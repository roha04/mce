#include "GlobalSystem.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace
{
std::atomic<float> g_calculationProgress{0.0f};

constexpr double kCoordTolerance = 1.0e-6;

std::vector<int> fixedDofs;

void applyEssentialBCsToDense(const int dofCount, std::vector<double>& matrix, std::vector<double>& rhs)
{
    const std::size_t n = static_cast<std::size_t>(dofCount);
    for (const int dof : fixedDofs)
    {
        if (dof < 0 || dof >= dofCount)
            continue;

        const std::size_t d = static_cast<std::size_t>(dof);
        for (std::size_t j = 0; j < n; ++j)
        {
            matrix[j * n + d] = 0.0;
            matrix[d * n + j] = 0.0;
        }

        matrix[d * n + d] = kPenaltyNumber;
        rhs[d] = 0.0;
    }
}

void gatherElementNodes(const Mesh& mesh, const Element& element, Node elementNodes[kHex20NodeCount])
{
    const std::vector<Node>& nodes = mesh.getNodes();
    for (int i = 0; i < kHex20NodeCount; ++i)
        elementNodes[i] = nodes[static_cast<std::size_t>(element.nodes[i])];
}
} // namespace

float getCalculationProgress()
{
    return g_calculationProgress.load(std::memory_order_relaxed);
}

void resetCalculationProgress()
{
    g_calculationProgress.store(0.0f, std::memory_order_relaxed);
}

void setCalculationProgress(const float value)
{
    g_calculationProgress.store(value, std::memory_order_relaxed);
}

void GlobalSystem::allocate(const int dofCount, const int halfBandwidth)
{
    dofCount_ = dofCount;
    halfBandwidth_ = halfBandwidth;
    kBand_.assign(static_cast<std::size_t>(dofCount_) * static_cast<std::size_t>(halfBandwidth_), 0.0);
    f_.assign(static_cast<std::size_t>(dofCount_), 0.0);
    u_.assign(static_cast<std::size_t>(dofCount_), 0.0);
}

int GlobalSystem::bandIndex(const int i, const int j) const
{
    if (i < 0 || j < i || j - i >= halfBandwidth_)
        return -1;
    return i * halfBandwidth_ + (j - i);
}

double GlobalSystem::getK(const int i, const int j) const
{
    if (i < 0 || j < 0 || i >= dofCount_ || j >= dofCount_)
        return 0.0;

    if (i <= j)
    {
        const int index = bandIndex(i, j);
        return index >= 0 ? kBand_[static_cast<std::size_t>(index)] : 0.0;
    }

    const int index = bandIndex(j, i);
    return index >= 0 ? kBand_[static_cast<std::size_t>(index)] : 0.0;
}

void GlobalSystem::setK(const int i, const int j, const double value)
{
    if (i > j)
        return setK(j, i, value);

    const int index = bandIndex(i, j);
    if (index < 0)
        return;
    kBand_[static_cast<std::size_t>(index)] = value;
}

void GlobalSystem::addK(const int i, const int j, const double value)
{
    if (i > j)
        return addK(j, i, value);

    const int index = bandIndex(i, j);
    if (index < 0)
        return;
    kBand_[static_cast<std::size_t>(index)] += value;
}

void GlobalSystem::addF(const int dof, const double value)
{
    if (dof < 0 || dof >= dofCount_)
        return;
    f_[static_cast<std::size_t>(dof)] += value;
}

int GlobalSystem::globalDof(const int nodeId, const int component)
{
    return nodeId * 3 + component;
}

void GlobalSystem::assembleElement(const Element& element,
                                   const ElementStiffnessMatrix& Ke,
                                   const ElementLoadVector& Fe)
{
    int globalDofs[kElementDofCount];
    for (int localNode = 0; localNode < kHex20NodeCount; ++localNode)
    {
        const int globalNode = element.nodes[localNode];
        globalDofs[3 * localNode + 0] = globalDof(globalNode, 0);
        globalDofs[3 * localNode + 1] = globalDof(globalNode, 1);
        globalDofs[3 * localNode + 2] = globalDof(globalNode, 2);
    }

    for (int i = 0; i < kElementDofCount; ++i)
    {
        addF(globalDofs[i], Fe[static_cast<std::size_t>(i)]);

        for (int j = i; j < kElementDofCount; ++j)
            addK(globalDofs[i], globalDofs[j], Ke[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]);
    }
}

void GlobalSystem::applyPenaltyToDof(const int dof)
{
    if (dof < 0 || dof >= dofCount_)
        return;

    addK(dof, dof, kPenaltyNumber);
    f_[static_cast<std::size_t>(dof)] = 0.0;
}

void GlobalSystem::applyFixedFaceY0(const Mesh& mesh)
{
    fixedDofs.clear();
    const std::vector<Node>& nodes = mesh.getNodes();

    for (std::size_t i = 0; i < nodes.size(); ++i)
    {
        if (std::abs(nodes[i].y) >= kCoordTolerance)
            continue;

        const int nodeId = static_cast<int>(i);
        for (int component = 0; component < 3; ++component)
        {
            const int dof = globalDof(nodeId, component);
            applyPenaltyToDof(dof);
            fixedDofs.push_back(dof);
        }
    }
}

void GlobalSystem::applyOrphanConstraints(const Mesh& mesh)
{
    std::vector<bool> used(mesh.nodeCount(), false);

    for (const Element& element : mesh.getElements())
    {
        for (int local = 0; local < kHex20NodeCount; ++local)
            used[static_cast<std::size_t>(element.nodes[local])] = true;
    }

    orphanNodeCount_ = 0;
    for (std::size_t nodeId = 0; nodeId < used.size(); ++nodeId)
    {
        if (used[nodeId])
            continue;

        ++orphanNodeCount_;
        applyPenaltyToDof(globalDof(static_cast<int>(nodeId), 0));
        applyPenaltyToDof(globalDof(static_cast<int>(nodeId), 1));
        applyPenaltyToDof(globalDof(static_cast<int>(nodeId), 2));
    }
}

void GlobalSystem::expandBandToDense(std::vector<double>& dense) const
{
    dense.assign(static_cast<std::size_t>(dofCount_) * static_cast<std::size_t>(dofCount_), 0.0);

    for (int i = 0; i < dofCount_; ++i)
    {
        const int colLimit = std::min(i + halfBandwidth_ - 1, dofCount_ - 1);
        for (int j = i; j <= colLimit; ++j)
        {
            const double value = getK(i, j);
            dense[static_cast<std::size_t>(i) * static_cast<std::size_t>(dofCount_) + static_cast<std::size_t>(j)] = value;
            if (i != j)
                dense[static_cast<std::size_t>(j) * static_cast<std::size_t>(dofCount_) + static_cast<std::size_t>(i)] = value;
        }
    }
}

void GlobalSystem::solveBandedGauss()
{
    std::vector<double> matrix;
    expandBandToDense(matrix);
    u_ = f_;
    applyEssentialBCsToDense(dofCount_, matrix, u_);

    const std::size_t n = static_cast<std::size_t>(dofCount_);
    auto at = [&](std::size_t row, std::size_t col) -> double& {
        return matrix[row * n + col];
    };

    for (std::size_t i = 0; i < n; ++i)
    {
        g_calculationProgress.store(
            0.60f + 0.30f * static_cast<float>(i) / static_cast<float>(n > 0 ? n : 1),
            std::memory_order_relaxed);

        const double pivot = at(i, i);
        if (std::abs(pivot) < 1.0e-14)
            throw std::runtime_error("Zero pivot in banded Gauss solver");

        const std::size_t rowLimit = std::min(i + static_cast<std::size_t>(halfBandwidth_ - 1), n - 1);
        for (std::size_t k = i + 1; k <= rowLimit; ++k)
        {
            const double coeff = at(k, i) / pivot;
            if (std::abs(coeff) < 1.0e-30)
                continue;

            const std::size_t colLimit = std::min(k + static_cast<std::size_t>(halfBandwidth_ - 1), n - 1);
            for (std::size_t j = i; j <= colLimit; ++j)
                at(k, j) -= coeff * at(i, j);

            u_[k] -= coeff * u_[i];
        }
    }

    for (int i = dofCount_ - 1; i >= 0; --i)
    {
        const std::size_t ui = static_cast<std::size_t>(i);
        g_calculationProgress.store(
            0.90f + 0.09f * static_cast<float>(dofCount_ - i) / static_cast<float>(dofCount_ > 0 ? dofCount_ : 1),
            std::memory_order_relaxed);

        const double pivot = at(ui, ui);
        const std::size_t colLimit = std::min(ui + static_cast<std::size_t>(halfBandwidth_ - 1), n - 1);

        double sum = u_[ui];
        for (std::size_t j = ui + 1; j <= colLimit; ++j)
            sum -= at(ui, j) * u_[j];

        u_[ui] = sum / pivot;
    }
}

void GlobalSystem::buildFromMesh(const Mesh& mesh,
                                 const MaterialProperties& material,
                                 const double tractionY)
{
    g_calculationProgress.store(0.05f, std::memory_order_relaxed);

    const int halfBandwidth = mesh.computeHalfBandwidth();
    allocate(static_cast<int>(mesh.nodeCount()) * 3, halfBandwidth);

    const GaussShapeDerivativeCache& gaussCache = GaussShapeDerivativeCache::instance();
    Node elementNodes[kHex20NodeCount];

    const std::size_t elementCount = mesh.getElements().size();
    std::size_t elementIndex = 0;

    for (const Element& element : mesh.getElements())
    {
        gatherElementNodes(mesh, element, elementNodes);

        const ElementStiffnessMatrix Ke = assembleElementStiffnessMatrix(elementNodes, material, gaussCache);

        ElementLoadVector Fe{};
        Fe.fill(0.0);
        if (isElementOnTopFace(elementNodes, mesh.getLy()))
            Fe = assembleElementLoadVectorTopFace(elementNodes, -tractionY);

        assembleElement(element, Ke, Fe);

        ++elementIndex;
        if (elementCount > 0)
        {
            const float assemblyProgress =
                0.05f + 0.50f * static_cast<float>(elementIndex) / static_cast<float>(elementCount);
            g_calculationProgress.store(assemblyProgress, std::memory_order_relaxed);
        }
    }

    g_calculationProgress.store(0.58f, std::memory_order_relaxed);
    applyFixedFaceY0(mesh);
    applyOrphanConstraints(mesh);
    g_calculationProgress.store(0.60f, std::memory_order_relaxed);
    solveBandedGauss();
    g_calculationProgress.store(0.99f, std::memory_order_relaxed);
}

double GlobalSystem::maxDisplacementY() const
{
    double maxValue = 0.0;
    for (int node = 0; node < dofCount_ / 3; ++node)
    {
        const double value = std::abs(u_[static_cast<std::size_t>(globalDof(node, 1))]);
        maxValue = std::max(maxValue, value);
    }
    return maxValue;
}

double GlobalSystem::totalReactionForceY() const
{
    double sum = 0.0;
    for (int i = 0; i < dofCount_; i += 3)
        sum += f_[static_cast<std::size_t>(i + 1)];
    return sum;
}

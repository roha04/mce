#include "GlobalSystem.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace
{
std::atomic<float> g_calculationProgress{0.0f};

constexpr double kCoordTolerance = 1.0e-6;

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

void GlobalSystem::applyBoundaryConditions(const Mesh& mesh)
{
    applyFixedFaceZU(mesh);
}

void GlobalSystem::applyFixedFaceZU(const Mesh& mesh)
{
    for (const int nodeId : mesh.getZU())
    {
        if (nodeId < 0 || static_cast<std::size_t>(nodeId) >= mesh.nodeCount())
            continue;

        for (int component = 0; component < 3; ++component)
            applyPenaltyToDof(globalDof(nodeId, component));
    }
}

void GlobalSystem::solveBandedGauss()
{
    u_ = f_;
    const int n = dofCount_;
    if (n <= 0)
        return;

    const int bandHalfWidth = halfBandwidth_ - 1;

    for (int i = 0; i < n; ++i)
    {
        g_calculationProgress.store(
            0.60f + 0.30f * static_cast<float>(i) / static_cast<float>(n),
            std::memory_order_relaxed);

        const double pivot = getK(i, i);
        if (std::abs(pivot) < 1.0e-14)
            throw std::runtime_error("Zero pivot in banded Gauss solver");

        const int rowLimit = std::min(i + bandHalfWidth, n - 1);
        for (int k = i + 1; k <= rowLimit; ++k)
        {
            const double coeff = getK(i, k) / pivot;
            if (std::abs(coeff) < 1.0e-30)
                continue;

            const int colLimit = std::min(k + bandHalfWidth, n - 1);
            for (int j = k; j <= colLimit; ++j)
                addK(k, j, -coeff * getK(i, j));

            u_[static_cast<std::size_t>(k)] -= coeff * u_[static_cast<std::size_t>(i)];
        }
    }

    for (int i = n - 1; i >= 0; --i)
    {
        g_calculationProgress.store(
            0.90f + 0.09f * static_cast<float>(n - i) / static_cast<float>(n),
            std::memory_order_relaxed);

        const double pivot = getK(i, i);
        const int colLimit = std::min(i + bandHalfWidth, n - 1);

        double sum = u_[static_cast<std::size_t>(i)];
        for (int j = i + 1; j <= colLimit; ++j)
            sum -= getK(i, j) * u_[static_cast<std::size_t>(j)];

        u_[static_cast<std::size_t>(i)] = sum / pivot;
    }
}

void GlobalSystem::assembleFromMesh(const Mesh& mesh,
                                    const MaterialProperties& material,
                                    const bool includeSurfaceLoads)
{
    const int halfBandwidth = mesh.computeHalfBandwidth();
    allocate(static_cast<int>(mesh.nodeCount()) * 3, halfBandwidth);

    const GaussShapeDerivativeCache& dfiabg = GaussShapeDerivativeCache::instance();
    Node elementNodes[kHex20NodeCount];

    const std::vector<Element>& elements = mesh.getElements();
    const std::size_t elementCount = elements.size();

    for (const Element& element : elements)
    {
        gatherElementNodes(mesh, element, elementNodes);

        const ElementStiffnessMatrix Ke = assembleElementStiffnessMatrix(elementNodes, material, dfiabg);

        ElementLoadVector Fe{};
        Fe.fill(0.0);
        assembleElement(element, Ke, Fe);
    }

    if (!includeSurfaceLoads)
        return;

    const std::vector<Mesh::ZPEntry>& zp = mesh.getZP();
    for (const Mesh::ZPEntry& entry : zp)
    {
        if (entry.elementIndex < 0 || static_cast<std::size_t>(entry.elementIndex) >= elementCount)
            continue;

        const Element& element = elements[static_cast<std::size_t>(entry.elementIndex)];
        gatherElementNodes(mesh, element, elementNodes);

        const ElementLoadVector Fe =
            assembleElementLoadVectorTopFace(elementNodes, -entry.pressure, entry.faceIndex);

        ElementStiffnessMatrix KeZero{};
        for (auto& row : KeZero)
            row.fill(0.0);
        assembleElement(element, KeZero, Fe);
    }
}

void GlobalSystem::buildFromMesh(const Mesh& mesh,
                                 const MaterialProperties& material,
                                 const double /*tractionY*/)
{
    g_calculationProgress.store(0.05f, std::memory_order_relaxed);

    assembleFromMesh(mesh, material, true);

    g_calculationProgress.store(0.58f, std::memory_order_relaxed);
    applyBoundaryConditions(mesh);
    g_calculationProgress.store(0.60f, std::memory_order_relaxed);
    solveBandedGauss();
    g_calculationProgress.store(0.99f, std::memory_order_relaxed);
}

void GlobalSystem::setForceAndSolve(const std::vector<double>& force)
{
    if (static_cast<int>(force.size()) != dofCount_)
        throw std::runtime_error("Force vector size mismatch in setForceAndSolve");

    f_ = force;
    solveBandedGauss();
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

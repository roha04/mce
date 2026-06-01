#include "StressAnalysis.h"

#include <algorithm>
#include <cmath>

namespace
{
void computeStrainFromDerivatives(const double dNGlobal[kHex20NodeCount][3],
                                  const double elementDisplacements[kElementDofCount],
                                  double strain[6])
{
    for (int c = 0; c < 6; ++c)
        strain[c] = 0.0;

    for (int k = 0; k < kHex20NodeCount; ++k)
    {
        const double u = elementDisplacements[3 * k + 0];
        const double v = elementDisplacements[3 * k + 1];
        const double w = elementDisplacements[3 * k + 2];

        strain[0] += dNGlobal[k][0] * u;
        strain[1] += dNGlobal[k][1] * v;
        strain[2] += dNGlobal[k][2] * w;
        strain[3] += dNGlobal[k][1] * u + dNGlobal[k][0] * v;
        strain[4] += dNGlobal[k][2] * v + dNGlobal[k][1] * w;
        strain[5] += dNGlobal[k][2] * u + dNGlobal[k][0] * w;
    }
}

StressTensor applyHookeLaw(const double strain[6], const double D[6][6])
{
    StressTensor stress{};
    const double stressVector[6] = {
        D[0][0] * strain[0] + D[0][1] * strain[1] + D[0][2] * strain[2],
        D[1][0] * strain[0] + D[1][1] * strain[1] + D[1][2] * strain[2],
        D[2][0] * strain[0] + D[2][1] * strain[1] + D[2][2] * strain[2],
        D[3][3] * strain[3],
        D[4][4] * strain[4],
        D[5][5] * strain[5]};

    stress.sxx = stressVector[0];
    stress.syy = stressVector[1];
    stress.szz = stressVector[2];
    stress.sxy = stressVector[3];
    stress.syz = stressVector[4];
    stress.sxz = stressVector[5];
    return stress;
}
} // namespace

double computeMaxPrincipalStress(const StressTensor& stress)
{
    const double I1 = stress.sxx + stress.syy + stress.szz;
    const double I2 = stress.sxx * stress.syy + stress.syy * stress.szz + stress.szz * stress.sxx
                    - stress.sxy * stress.sxy - stress.syz * stress.syz - stress.sxz * stress.sxz;
    const double I3 = stress.sxx * stress.syy * stress.szz
                    + 2.0 * stress.sxy * stress.syz * stress.sxz
                    - stress.sxx * stress.syz * stress.syz
                    - stress.syy * stress.sxz * stress.sxz
                    - stress.szz * stress.sxy * stress.sxy;

    const double a = -I1;
    const double b = I2;
    const double c = -I3;

    const double shift = a / 3.0;
    const double p = b - a * a / 3.0;
    const double q = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c;

    const double discriminant = q * q / 4.0 + p * p * p / 27.0;
    if (discriminant >= 0.0)
    {
        const double sqrtDisc = std::sqrt(discriminant);
        const double u = std::cbrt(-q / 2.0 + sqrtDisc);
        const double v = std::cbrt(-q / 2.0 - sqrtDisc);
        return u + v - shift;
    }

    const double r = std::sqrt(-p * p * p / 27.0);
    const double phi = std::acos(std::clamp(-q / (2.0 * r), -1.0, 1.0));
    const double m0 = 2.0 * std::cbrt(r);

    const double s1 = m0 * std::cos(phi / 3.0) - shift;
    const double s2 = m0 * std::cos((phi + 2.0 * 3.14159265358979323846) / 3.0) - shift;
    const double s3 = m0 * std::cos((phi + 4.0 * 3.14159265358979323846) / 3.0) - shift;
    return (std::max)({s1, s2, s3});
}

double computeMinPrincipalStress(const StressTensor& stress)
{
    const double I1 = stress.sxx + stress.syy + stress.szz;
    const double I2 = stress.sxx * stress.syy + stress.syy * stress.szz + stress.szz * stress.sxx
                    - stress.sxy * stress.sxy - stress.syz * stress.syz - stress.sxz * stress.sxz;
    const double I3 = stress.sxx * stress.syy * stress.szz
                    + 2.0 * stress.sxy * stress.syz * stress.sxz
                    - stress.sxx * stress.syz * stress.syz
                    - stress.syy * stress.sxz * stress.sxz
                    - stress.szz * stress.sxy * stress.sxy;

    const double a = -I1;
    const double b = I2;
    const double c = -I3;

    const double shift = a / 3.0;
    const double p = b - a * a / 3.0;
    const double q = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c;

    const double discriminant = q * q / 4.0 + p * p * p / 27.0;
    if (discriminant >= 0.0)
    {
        const double sqrtDisc = std::sqrt(discriminant);
        const double u = std::cbrt(-q / 2.0 + sqrtDisc);
        const double v = std::cbrt(-q / 2.0 - sqrtDisc);
        return u + v - shift;
    }

    const double r = std::sqrt(-p * p * p / 27.0);
    const double phi = std::acos(std::clamp(-q / (2.0 * r), -1.0, 1.0));
    const double m0 = 2.0 * std::cbrt(r);

    const double s1 = m0 * std::cos(phi / 3.0) - shift;
    const double s2 = m0 * std::cos((phi + 2.0 * 3.14159265358979323846) / 3.0) - shift;
    const double s3 = m0 * std::cos((phi + 4.0 * 3.14159265358979323846) / 3.0) - shift;
    return (std::min)({s1, s2, s3});
}

StressTensor StressAnalyzer::computeStressAtElementNode(
    const Node elementNodes[kHex20NodeCount],
    const double elementDisplacements[kElementDofCount],
    const int localNode,
    const MaterialProperties& material)
{
    double xi = 0.0;
    double eta = 0.0;
    double zeta = 0.0;
    getHex20LocalNodeCoordinates(localNode, xi, eta, zeta);

    double dNLocal[kHex20NodeCount][3];
    computeHex20ShapeDerivativesLocal(xi, eta, zeta, dNLocal);

    Jacobian3x3 jacobian = computeJacobian(dNLocal, elementNodes);
    double invJ[3][3];
    if (!jacobian.invert(invJ))
        return {};

    double dNGlobal[kHex20NodeCount][3];
    transformShapeDerivativesToGlobal(dNLocal, invJ, dNGlobal);

    double strain[6];
    computeStrainFromDerivatives(dNGlobal, elementDisplacements, strain);

    double D[6][6];
    buildIsotropicElasticityMatrix6x6(material, D);
    return applyHookeLaw(strain, D);
}

void StressAnalyzer::compute(const Mesh& mesh,
                             const GlobalSystem& system,
                             const MaterialProperties& material)
{
    nodeResults_.assign(mesh.nodeCount(), NodeStressResult{});
    contributionCount_.assign(mesh.nodeCount(), 0);
    maxPrincipal_ = 0.0;
    minPrincipal_ = 0.0;
    bool hasPrincipal = false;

    const std::vector<double>& u = system.displacements();
    Node elementNodes[kHex20NodeCount];
    double elementDisplacements[kElementDofCount];

    for (const Element& element : mesh.getElements())
    {
        for (int i = 0; i < kHex20NodeCount; ++i)
            elementNodes[i] = mesh.getNodes()[static_cast<std::size_t>(element.nodes[i])];

        for (int i = 0; i < kElementDofCount; ++i)
        {
            const int localNode = i / 3;
            const int component = i % 3;
            const int globalNode = element.nodes[localNode];
            elementDisplacements[i] = u[static_cast<std::size_t>(globalNode * 3 + component)];
        }

        for (int localNode = 0; localNode < kHex20NodeCount; ++localNode)
        {
            const int globalNode = element.nodes[localNode];
            const StressTensor localStress =
                computeStressAtElementNode(elementNodes, elementDisplacements, localNode, material);

            NodeStressResult& result = nodeResults_[static_cast<std::size_t>(globalNode)];
            result.stress.sxx += localStress.sxx;
            result.stress.syy += localStress.syy;
            result.stress.szz += localStress.szz;
            result.stress.sxy += localStress.sxy;
            result.stress.syz += localStress.syz;
            result.stress.sxz += localStress.sxz;
            ++contributionCount_[static_cast<std::size_t>(globalNode)];
        }
    }

    for (std::size_t node = 0; node < nodeResults_.size(); ++node)
    {
        if (contributionCount_[node] == 0)
            continue;

        const double inv = 1.0 / static_cast<double>(contributionCount_[node]);
        NodeStressResult& result = nodeResults_[node];
        result.stress.sxx *= inv;
        result.stress.syy *= inv;
        result.stress.szz *= inv;
        result.stress.sxy *= inv;
        result.stress.syz *= inv;
        result.stress.sxz *= inv;
        result.principalMax = computeMaxPrincipalStress(result.stress);
        result.principalMin = computeMinPrincipalStress(result.stress);
        result.valid = true;
        maxPrincipal_ = std::max(maxPrincipal_, result.principalMax);
        if (!hasPrincipal)
        {
            minPrincipal_ = result.principalMin;
            hasPrincipal = true;
        }
        else
            minPrincipal_ = std::min(minPrincipal_, result.principalMin);
    }
}

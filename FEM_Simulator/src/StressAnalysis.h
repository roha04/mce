#pragma once

// =============================================================================
// StressAnalysis — етап 4 (заняття 12): напруження та головні напруження
// =============================================================================
// За знайденими U обчислюється тензор напружень у вузлах (закон Гука),
// далі — головні напруження σ1, σ2, σ3 з кубічного рівняння інваріантів (49).
// =============================================================================

#include "GlobalSystem.h"
#include "MathCore.h"

#include <vector>

/// Компоненти симетричного тензора напружень у декартовій системі.
struct StressTensor
{
    double sxx = 0.0;
    double syy = 0.0;
    double szz = 0.0;
    double sxy = 0.0;
    double syz = 0.0;
    double sxz = 0.0;
};

/// Результат у вузлі після усереднення по суміжних елементах.
struct NodeStressResult
{
    StressTensor stress{};
    double principalMax = 0.0;  // σ1 — максимальне головне (розтяг)
    double principalMin = 0.0;  // σ3 — мінімальне головне (стиск)
    bool valid = false;
};

class StressAnalyzer
{
public:
    /// Обходить усі СЕ, накопичує напруження в вузлах, усереднює, рахує σ1 та σ3.
    void compute(const Mesh& mesh, const GlobalSystem& system, const MaterialProperties& material);

    const std::vector<NodeStressResult>& nodeResults() const { return nodeResults_; }
    double maxPrincipalStress() const { return maxPrincipal_; }
    double minPrincipalStress() const { return minPrincipal_; }

private:
    /// Напруження в локальному вузлі СЕ за U елемента (похідні ψ у вузлі, не в точках Гауса).
    static StressTensor computeStressAtElementNode(const Node elementNodes[kHex20NodeCount],
                                                 const double elementDisplacements[kElementDofCount],
                                                 int localNode,
                                                 const MaterialProperties& material);

    std::vector<NodeStressResult> nodeResults_;
    std::vector<int> contributionCount_;
    double maxPrincipal_ = 0.0;
    double minPrincipal_ = 0.0;
};

/// Корені кубічного рівняння інваріантів: σ³ - I1·σ² + I2·σ - I3 = 0.
double computeMaxPrincipalStress(const StressTensor& stress);
double computeMinPrincipalStress(const StressTensor& stress);

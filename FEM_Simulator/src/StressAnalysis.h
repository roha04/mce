#pragma once

#include "GlobalSystem.h"
#include "MathCore.h"

#include <vector>

struct StressTensor
{
    double sxx = 0.0;
    double syy = 0.0;
    double szz = 0.0;
    double sxy = 0.0;
    double syz = 0.0;
    double sxz = 0.0;
};

struct NodeStressResult
{
    StressTensor stress{};
    double principalMax = 0.0;
    bool valid = false;
};

// Обчислення напружень у вузлах (Етап 4, Заняття 12).
class StressAnalyzer
{
public:
    void compute(const Mesh& mesh, const GlobalSystem& system, const MaterialProperties& material);

    const std::vector<NodeStressResult>& nodeResults() const { return nodeResults_; }
    double maxPrincipalStress() const { return maxPrincipal_; }

private:
    static StressTensor computeStressAtElementNode(const Node elementNodes[kHex20NodeCount],
                                                 const double elementDisplacements[kElementDofCount],
                                                 int localNode,
                                                 const MaterialProperties& material);

    std::vector<NodeStressResult> nodeResults_;
    std::vector<int> contributionCount_;
    double maxPrincipal_ = 0.0;
};

double computeMaxPrincipalStress(const StressTensor& stress);

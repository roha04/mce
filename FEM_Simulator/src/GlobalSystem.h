#pragma once

#include "MathCore.h"

#include <vector>

constexpr double kPenaltyNumber = 1.0e30;

// Глобальна система K·U = F у стрічковому форматі (Етап 3).
class GlobalSystem
{
public:
    void buildFromMesh(const Mesh& mesh, const MaterialProperties& material, double tractionY);

    int dofCount() const { return dofCount_; }
    int halfBandwidth() const { return halfBandwidth_; }
    int orphanNodeCount() const { return orphanNodeCount_; }

    const std::vector<double>& displacements() const { return u_; }
    const std::vector<double>& forceVector() const { return f_; }
    double maxDisplacementY() const;
    double totalReactionForceY() const;

private:
    void allocate(int dofCount, int halfBandwidth);

    // Конвертація (i, j) -> 1D індекс стрічкового масиву (лише j >= i).
    int bandIndex(int i, int j) const;
    double getK(int i, int j) const;
    void setK(int i, int j, double value);
    void addK(int i, int j, double value);
    void addF(int dof, double value);

    static int globalDof(int nodeId, int component);

    void assembleElement(const Element& element,
                         const ElementStiffnessMatrix& Ke,
                         const ElementLoadVector& Fe);

    void applyFixedFaceY0(const Mesh& mesh);
    void applyOrphanConstraints(const Mesh& mesh);
    void applyPenaltyToDof(int dof);

    void expandBandToDense(std::vector<double>& dense) const;
    void solveBandedGauss();

    int dofCount_ = 0;
    int halfBandwidth_ = 0;
    int orphanNodeCount_ = 0;

    std::vector<double> kBand_;
    std::vector<double> f_;
    std::vector<double> u_;
};

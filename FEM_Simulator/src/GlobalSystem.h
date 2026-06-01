#pragma once

#include "MathCore.h"

#include <vector>

constexpr double kPenaltyNumber = 1.0e30;

float getCalculationProgress();
void setCalculationProgress(float value);
void resetCalculationProgress();

// Глобальна система K·U = F у стрічковому форматі MG (заняття 9–11).
class GlobalSystem
{
public:
    void buildFromMesh(const Mesh& mesh, const MaterialProperties& material, double tractionY);

    // Збірка MG та F без крайових умов і без розв'язання (тест зан. 13).
    void assembleFromMesh(const Mesh& mesh, const MaterialProperties& material, bool includeSurfaceLoads);

    void applyBoundaryConditions(const Mesh& mesh);

    int dofCount() const { return dofCount_; }
    int halfBandwidth() const { return halfBandwidth_; }

    const std::vector<double>& displacements() const { return u_; }
    const std::vector<double>& forceVector() const { return f_; }
    const std::vector<double>& stiffnessBand() const { return kBand_; }

    double maxDisplacementY() const;
    double totalReactionForceY() const;

    // Перезапуск розв'язувача з новим F (налагодження, заняття 13).
    void setForceAndSolve(const std::vector<double>& force);

private:
    void allocate(int dofCount, int halfBandwidth);

    int bandIndex(int i, int j) const;
    double getK(int i, int j) const;
    void setK(int i, int j, double value);
    void addK(int i, int j, double value);
    void addF(int dof, double value);

    static int globalDof(int nodeId, int component);

    void assembleElement(const Element& element,
                         const ElementStiffnessMatrix& Ke,
                         const ElementLoadVector& Fe);

    void applyFixedFaceZU(const Mesh& mesh);
    void applyPenaltyToDof(int dof);

    void solveBandedGauss();

    int dofCount_ = 0;
    int halfBandwidth_ = 0;

    std::vector<double> kBand_;
    std::vector<double> f_;
    std::vector<double> u_;
};

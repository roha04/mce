#pragma once

// =============================================================================
// GlobalSystem — етап 3 методички: глобальна збірка MG, F та розв'язання K·U = F
// =============================================================================
// MG зберігається у стрічковому (ленточному) вигляді: лише коефіцієнти з
// j - i < ng (півширина стрічки). Розв'язувач — власний метод Гауса під цей формат.
// Крайові умови: метод штрафу (заняття 10) для вузлів ZU.
// =============================================================================

#include "MathCore.h"

#include <vector>

/// Велике число на діагональ для закріплення DOF (метод 2, заняття 10).
constexpr double kPenaltyNumber = 1.0e30;

float getCalculationProgress();
void setCalculationProgress(float value);
void resetCalculationProgress();

class GlobalSystem
{
public:
    /// Повний цикл: збірка K і F → штраф на ZU → розв'язання → U.
    void buildFromMesh(const Mesh& mesh, const MaterialProperties& material, double tractionY);

    /// Лише збірка глобальної матриці MG та вектора F (без КУ і без розв'язання).
    /// includeSurfaceLoads: додати внешнє навантаження з масиву ZP.
    void assembleFromMesh(const Mesh& mesh, const MaterialProperties& material, bool includeSurfaceLoads);

    /// Накладає штрафні КУ для всіх DOF вузлів з ZU (Ux=Uy=Uz=0).
    void applyBoundaryConditions(const Mesh& mesh);

    int dofCount() const { return dofCount_; }
    int halfBandwidth() const { return halfBandwidth_; }

    const std::vector<double>& displacements() const { return u_; }
    const std::vector<double>& forceVector() const { return f_; }
    const std::vector<double>& stiffnessBand() const { return kBand_; }

    double maxDisplacementY() const;
    double totalReactionForceY() const;

    /// Підставляє новий F і повторно вирішує K·U = F (налагодження, заняття 13).
    void setForceAndSolve(const std::vector<double>& force);

private:
    void allocate(int dofCount, int halfBandwidth);

    /// Індекс у 1D-масиві kBand_ для пари (i,j), j >= i; -1 якщо поза стрічкою.
    int bandIndex(int i, int j) const;
    double getK(int i, int j) const;
    void setK(int i, int j, double value);
    void addK(int i, int j, double value);
    void addF(int dof, double value);

    /// Глобальний номер ступеня свободи: вузол nodeId, компонента 0..2 (x,y,z).
    static int globalDof(int nodeId, int component);

    /// Розкидає K^e та F^e елемента в глобальні MG і F за NT.
    void assembleElement(const Element& element,
                         const ElementStiffnessMatrix& Ke,
                         const ElementLoadVector& Fe);

    void applyFixedFaceZU(const Mesh& mesh);
    void applyPenaltyToDof(int dof);

    /// Прямий і зворотний хід Гауса без розгортання MG у повну N×N матрицю.
    void solveBandedGauss();

    int dofCount_ = 0;
    int halfBandwidth_ = 0;

    std::vector<double> kBand_;
    std::vector<double> f_;
    std::vector<double> u_;
};

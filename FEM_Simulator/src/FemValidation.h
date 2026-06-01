#pragma once

// =============================================================================
// FemValidation — автоматичні перевірки з заняття 13 методички
// =============================================================================

#include <string>

struct FemValidationReport
{
    bool jacobianUnitCubeOk = false;
    double jacobianUnitCubeDet = 0.0;

    bool solverUnitRhsOk = false;
    double solverUnitRhsMaxError = 0.0;

    bool linearLoadScalingNote = true;
};

/// Запускає тест det(J) та тест розв'язувача F_i = Σ_j K_ij → U_i ≈ 1.
FemValidationReport runFemValidationChecks();

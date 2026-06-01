#pragma once

#include <string>

struct FemValidationReport
{
    bool jacobianUnitCubeOk = false;
    double jacobianUnitCubeDet = 0.0;

    bool solverUnitRhsOk = false;
    double solverUnitRhsMaxError = 0.0;

    bool linearLoadScalingNote = true;
};

FemValidationReport runFemValidationChecks();

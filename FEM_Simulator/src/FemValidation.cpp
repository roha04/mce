#include "FemValidation.h"

#include "GlobalSystem.h"
#include "MathCore.h"

#include <algorithm>
#include <cmath>
#include <vector>

// Автотести заняття 13: коректність Якобіана та стрічкового розв'язувача.

FemValidationReport runFemValidationChecks()
{
    FemValidationReport report{};

    // Тест 1: одиничний куб 4×4×4 м, det(J) у центрі (ξ=η=ζ=0) має бути ≈ 8.
    report.jacobianUnitCubeDet = verifyUnitCubeJacobianDeterminant();
    report.jacobianUnitCubeOk = std::abs(report.jacobianUnitCubeDet - 8.0) < 0.05;

    try
    {
        // Тест 2: мінімальна сітка 1×1×1, зібрати K без поверхневого навантаження, накласти ZU.
        Mesh mesh;
        mesh.generateRectangularParallelepiped(1.0, 1.0, 1.0, 1, 1, 1);
        mesh.buildBoundaryData(0.0);

        MaterialProperties material;
        material.youngModulus = 2.0e11;
        material.poissonRatio = 0.3;

        GlobalSystem system;
        system.assembleFromMesh(mesh, material, false);
        system.applyBoundaryConditions(mesh);

        const int n = system.dofCount();
        const int hb = system.halfBandwidth();
        const std::vector<double>& band = system.stiffnessBand();

        auto bandGet = [&](int i, int j) -> double {
            if (i > j)
                std::swap(i, j);
            if (j - i >= hb)
                return 0.0;
            return band[static_cast<std::size_t>(i * hb + (j - i))];
        };

        // F_i = сума елементів рядка K (для симетричної K це дає K·1); очікуємо U≈1 на вільних DOF.
        std::vector<double> unitRhs(static_cast<std::size_t>(n), 0.0);
        for (int i = 0; i < n; ++i)
        {
            double rowSum = 0.0;
            const int jMax = std::min(i + hb - 1, n - 1);
            for (int j = 0; j <= jMax; ++j)
                rowSum += bandGet(i, j);
            unitRhs[static_cast<std::size_t>(i)] = rowSum;
        }

        system.setForceAndSolve(unitRhs);

        const std::vector<double>& u = system.displacements();
        double maxError = 0.0;
        for (int i = 0; i < n; ++i)
            maxError = std::max(maxError, std::abs(u[static_cast<std::size_t>(i)] - 1.0));

        report.solverUnitRhsMaxError = maxError;
        report.solverUnitRhsOk = maxError < 1.0e-9;
    }
    catch (...)
    {
        report.solverUnitRhsOk = false;
        report.solverUnitRhsMaxError = 1.0;
    }

    return report;
}

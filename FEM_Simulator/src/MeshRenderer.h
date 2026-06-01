#pragma once

#include "MathCore.h"
#include "StressAnalysis.h"

#include <glad/glad.h>

#include <vector>

class GlobalSystem;

class MeshRenderer
{
public:
    void build(const Mesh& mesh, const GlobalSystem& system, const StressAnalyzer& stressAnalyzer);

    void updateIndices(const Mesh& mesh, bool showSectionX, bool showSectionZ);

    void render(float scaleFactor,
                float aspectRatio,
                float cameraDist,
                float camAngleX,
                float camAngleY,
                float camPanX,
                float camPanY) const;
    void shutdown();

    float minStress() const { return minStress_; }
    float maxStress() const { return maxStress_; }

private:
    struct SurfaceVertex
    {
        float baseX, baseY, baseZ; // BasePosition (оригінальні x, y, z)
        float dispX, dispY, dispZ; // Displacement (Ux, Uy, Uz) — без масштабу
        float stress;              // головне напруження σ₁
    };

    void uploadIndices();

    GLuint vao_ = 0;
    GLuint vbo_ = 0;
    GLuint ebo_ = 0;
    GLuint shaderProgram_ = 0;

    GLint locProjection_ = -1;
    GLint locView_ = -1;
    GLint locModel_ = -1;
    GLint locScaleFactor_ = -1;
    GLint locMinStress_ = -1;
    GLint locMaxStress_ = -1;
    GLint locIsWireframe_ = -1;

    float minStress_ = 0.0f;
    float maxStress_ = 1.0f;

    float modelCenterX_ = 0.0f;
    float modelCenterY_ = 0.0f;
    float modelCenterZ_ = 0.0f;
    float modelExtent_ = 1.0f;

    std::vector<SurfaceVertex> vertices_;
    std::vector<unsigned int> indices_;
    int indexCount_ = 0;
};

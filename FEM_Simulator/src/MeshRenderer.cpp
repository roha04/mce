#include "MeshRenderer.h"

#include "GlobalSystem.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace
{
constexpr double kPlaneTolerance = 1.0e-6;

constexpr const char* kVertexShader = R"(#version 330 core
layout(location = 0) in vec3 aBasePosition;
layout(location = 1) in vec3 aDisplacement;
layout(location = 2) in float aStress;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_Model;
uniform float u_ScaleFactor;

out float vStress;

void main()
{
    vec3 final_pos = aBasePosition + aDisplacement * u_ScaleFactor;
    gl_Position = u_Projection * u_View * u_Model * vec4(final_pos, 1.0);
    vStress = aStress;
}
)";

constexpr const char* kFragmentShader = R"(#version 330 core
in float vStress;
out vec4 FragColor;

uniform float u_MinStress;
uniform float u_MaxStress;
uniform bool u_IsWireframe;

void main()
{
    if (u_IsWireframe)
    {
        FragColor = vec4(0.1, 0.1, 0.1, 1.0);
        return;
    }

    float normalized_stress = 0.0;
    if (u_MaxStress > u_MinStress)
        normalized_stress = clamp((vStress - u_MinStress) / (u_MaxStress - u_MinStress), 0.0, 1.0);

    vec3 color = mix(vec3(0.0, 0.0, 1.0), vec3(1.0, 0.0, 0.0), normalized_stress);
    FragColor = vec4(color, 1.0);
}
)";

// 6 граней HEX20 через 4 кутові вузли (локальні індекси 0..7).
constexpr int kHexCornerFaces[6][4] = {
    {0, 1, 2, 3}, // ζ = -1 (нижня)
    {4, 5, 6, 7}, // ζ = +1 (верхня)
    {0, 3, 7, 4}, // ξ = -1
    {1, 2, 6, 5}, // ξ = +1
    {0, 1, 5, 4}, // η = -1 (Y = 0)
    {3, 2, 6, 7}  // η = +1 (Y = Ly)
};

GLuint compileShader(const GLenum type, const char* source)
{
    const GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &source, nullptr);
    glCompileShader(shader);

    GLint ok = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &ok);
    if (!ok)
    {
        char log[1024];
        glGetShaderInfoLog(shader, sizeof(log), nullptr, log);
        throw std::runtime_error(std::string("Shader compile error: ") + log);
    }
    return shader;
}

GLuint createProgram()
{
    const GLuint vs = compileShader(GL_VERTEX_SHADER, kVertexShader);
    const GLuint fs = compileShader(GL_FRAGMENT_SHADER, kFragmentShader);
    const GLuint program = glCreateProgram();
    glAttachShader(program, vs);
    glAttachShader(program, fs);
    glLinkProgram(program);
    glDeleteShader(vs);
    glDeleteShader(fs);

    GLint ok = 0;
    glGetProgramiv(program, GL_LINK_STATUS, &ok);
    if (!ok)
        throw std::runtime_error("Shader link error");
    return program;
}

void makeIdentity(float out[16])
{
    for (int i = 0; i < 16; ++i)
        out[i] = 0.0f;
    out[0] = out[5] = out[10] = out[15] = 1.0f;
}

void makeTranslate(const float tx, const float ty, const float tz, float out[16])
{
    makeIdentity(out);
    out[12] = tx;
    out[13] = ty;
    out[14] = tz;
}

void makeRotateX(const float angle, float out[16])
{
    makeIdentity(out);
    const float c = std::cos(angle);
    const float s = std::sin(angle);
    out[5] = c;
    out[6] = s;
    out[9] = -s;
    out[10] = c;
}

void makeRotateY(const float angle, float out[16])
{
    makeIdentity(out);
    const float c = std::cos(angle);
    const float s = std::sin(angle);
    out[0] = c;
    out[2] = -s;
    out[8] = s;
    out[10] = c;
}

void multiplyMat4(const float a[16], const float b[16], float out[16])
{
    float temp[16]{};
    for (int col = 0; col < 4; ++col)
    {
        for (int row = 0; row < 4; ++row)
        {
            double sum = 0.0;
            for (int k = 0; k < 4; ++k)
                sum += static_cast<double>(a[k * 4 + row]) * static_cast<double>(b[col * 4 + k]);
            temp[col * 4 + row] = static_cast<float>(sum);
        }
    }
    for (int i = 0; i < 16; ++i)
        out[i] = temp[i];
}

void makePerspective(const float fovY, const float aspect, const float zNear, const float zFar, float out[16])
{
    for (int i = 0; i < 16; ++i)
        out[i] = 0.0f;

    const float f = 1.0f / std::tan(fovY * 0.5f);
    out[0] = f / aspect;
    out[5] = f;
    out[10] = (zFar + zNear) / (zNear - zFar);
    out[11] = -1.0f;
    out[14] = (2.0f * zFar * zNear) / (zNear - zFar);
}

void makeLookAt(const float eye[3], const float center[3], const float up[3], float out[16])
{
    float f[3] = {center[0] - eye[0], center[1] - eye[1], center[2] - eye[2]};
    const float fLen = std::sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
    f[0] /= fLen;
    f[1] /= fLen;
    f[2] /= fLen;

    float s[3] = {
        f[1] * up[2] - f[2] * up[1],
        f[2] * up[0] - f[0] * up[2],
        f[0] * up[1] - f[1] * up[0]};
    const float sLen = std::sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
    s[0] /= sLen;
    s[1] /= sLen;
    s[2] /= sLen;

    float u[3] = {
        s[1] * f[2] - s[2] * f[1],
        s[2] * f[0] - s[0] * f[2],
        s[0] * f[1] - s[1] * f[0]};

    out[0] = s[0];
    out[4] = s[1];
    out[8] = s[2];
    out[12] = -(s[0] * eye[0] + s[1] * eye[1] + s[2] * eye[2]);
    out[1] = u[0];
    out[5] = u[1];
    out[9] = u[2];
    out[13] = -(u[0] * eye[0] + u[1] * eye[1] + u[2] * eye[2]);
    out[2] = -f[0];
    out[6] = -f[1];
    out[10] = -f[2];
    out[14] = f[0] * eye[0] + f[1] * eye[1] + f[2] * eye[2];
    out[3] = 0.0f;
    out[7] = 0.0f;
    out[11] = 0.0f;
    out[15] = 1.0f;
}

void addQuadIndices(const Element& element, const int localCorners[4], std::vector<unsigned int>& indices)
{
    const unsigned int i0 = static_cast<unsigned int>(element.nodes[localCorners[0]]);
    const unsigned int i1 = static_cast<unsigned int>(element.nodes[localCorners[1]]);
    const unsigned int i2 = static_cast<unsigned int>(element.nodes[localCorners[2]]);
    const unsigned int i3 = static_cast<unsigned int>(element.nodes[localCorners[3]]);

    indices.push_back(i0);
    indices.push_back(i1);
    indices.push_back(i2);
    indices.push_back(i0);
    indices.push_back(i2);
    indices.push_back(i3);
}

void addElementFaceIndices(const Element& element, const int faceId, std::vector<unsigned int>& indices)
{
    addQuadIndices(element, kHexCornerFaces[faceId], indices);
}

void addElementAllFaces(const Element& element, std::vector<unsigned int>& indices)
{
    for (int face = 0; face < 6; ++face)
        addElementFaceIndices(element, face, indices);
}

void computeElementBounds(const Element& element, const Mesh& mesh,
                          double& minX, double& maxX,
                          double& minY, double& maxY,
                          double& minZ, double& maxZ)
{
    const std::vector<Node>& nodes = mesh.getNodes();
    minX = maxX = nodes[static_cast<std::size_t>(element.nodes[0])].x;
    minY = maxY = nodes[static_cast<std::size_t>(element.nodes[0])].y;
    minZ = maxZ = nodes[static_cast<std::size_t>(element.nodes[0])].z;

    for (int corner = 1; corner < 8; ++corner)
    {
        const Node& node = nodes[static_cast<std::size_t>(element.nodes[corner])];
        minX = (std::min)(minX, node.x);
        maxX = (std::max)(maxX, node.x);
        minY = (std::min)(minY, node.y);
        maxY = (std::max)(maxY, node.y);
        minZ = (std::min)(minZ, node.z);
        maxZ = (std::max)(maxZ, node.z);
    }
}

bool elementCrossesPlaneX(const Element& element, const Mesh& mesh, const double planeX)
{
    double minX, maxX, minY, maxY, minZ, maxZ;
    computeElementBounds(element, mesh, minX, maxX, minY, maxY, minZ, maxZ);
    (void)minY;
    (void)maxY;
    (void)minZ;
    (void)maxZ;
    return minX <= planeX + kPlaneTolerance && maxX >= planeX - kPlaneTolerance;
}

bool elementCrossesPlaneZ(const Element& element, const Mesh& mesh, const double planeZ)
{
    double minX, maxX, minY, maxY, minZ, maxZ;
    computeElementBounds(element, mesh, minX, maxX, minY, maxY, minZ, maxZ);
    (void)minX;
    (void)maxX;
    (void)minY;
    (void)maxY;
    return minZ <= planeZ + kPlaneTolerance && maxZ >= planeZ - kPlaneTolerance;
}
} // namespace

void MeshRenderer::build(const Mesh& mesh,
                         const GlobalSystem& system,
                         const StressAnalyzer& stressAnalyzer)
{
    shutdown();

    const std::vector<Node>& nodes = mesh.getNodes();
    const std::vector<double>& u = system.displacements();
    const auto& stressResults = stressAnalyzer.nodeResults();

    vertices_.resize(nodes.size());
    minStress_ = 0.0f;
    maxStress_ = 0.0f;
    bool hasValidStress = false;

    for (std::size_t nodeId = 0; nodeId < nodes.size(); ++nodeId)
    {
        SurfaceVertex& vertex = vertices_[nodeId];

        // BasePosition — лише оригінальні координати, без деформації.
        vertex.baseX = static_cast<float>(nodes[nodeId].x);
        vertex.baseY = static_cast<float>(nodes[nodeId].y);
        vertex.baseZ = static_cast<float>(nodes[nodeId].z);

        const NodeStressResult& stress = stressResults[nodeId];
        if (stress.valid)
        {
            // Displacement: U[i*3], U[i*3+1], U[i*3+2] — без множення на scale.
            vertex.dispX = static_cast<float>(u[nodeId * 3 + 0]);
            vertex.dispY = static_cast<float>(u[nodeId * 3 + 1]);
            vertex.dispZ = static_cast<float>(u[nodeId * 3 + 2]);
            vertex.stress = static_cast<float>(stress.principalMax);

            if (!hasValidStress)
            {
                minStress_ = maxStress_ = vertex.stress;
                hasValidStress = true;
            }
            else
            {
                minStress_ = (std::min)(minStress_, vertex.stress);
                maxStress_ = (std::max)(maxStress_, vertex.stress);
            }
        }
        else
        {
            // Orphan / невалідні вузли: нульове переміщення, щоб не ламати топологію.
            vertex.dispX = 0.0f;
            vertex.dispY = 0.0f;
            vertex.dispZ = 0.0f;
            vertex.stress = 0.0f;
        }
    }

    if (!hasValidStress)
        maxStress_ = minStress_ + 1.0f;

    modelCenterX_ = static_cast<float>(mesh.getLx() * 0.5);
    modelCenterY_ = static_cast<float>(mesh.getLy() * 0.5);
    modelCenterZ_ = static_cast<float>(mesh.getLz() * 0.5);
    modelExtent_ = static_cast<float>(
        std::sqrt(mesh.getLx() * mesh.getLx() + mesh.getLy() * mesh.getLy() + mesh.getLz() * mesh.getLz()));

    shaderProgram_ = createProgram();
    locProjection_ = glGetUniformLocation(shaderProgram_, "u_Projection");
    locView_ = glGetUniformLocation(shaderProgram_, "u_View");
    locModel_ = glGetUniformLocation(shaderProgram_, "u_Model");
    locScaleFactor_ = glGetUniformLocation(shaderProgram_, "u_ScaleFactor");
    locMinStress_ = glGetUniformLocation(shaderProgram_, "u_MinStress");
    locMaxStress_ = glGetUniformLocation(shaderProgram_, "u_MaxStress");
    locIsWireframe_ = glGetUniformLocation(shaderProgram_, "u_IsWireframe");
    glGenVertexArrays(1, &vao_);
    glGenBuffers(1, &vbo_);
    glGenBuffers(1, &ebo_);

    glBindVertexArray(vao_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(vertices_.size() * sizeof(SurfaceVertex)),
                 vertices_.data(),
                 GL_STATIC_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(SurfaceVertex), reinterpret_cast<void*>(0));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(SurfaceVertex), reinterpret_cast<void*>(3 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(SurfaceVertex), reinterpret_cast<void*>(6 * sizeof(float)));

    glBindVertexArray(0);

    updateIndices(mesh, false, false);
}

void MeshRenderer::updateIndices(const Mesh& mesh, const bool showSectionX, const bool showSectionZ)
{
    indices_.clear();

    const int Nx = mesh.getNx();
    const int Ny = mesh.getNy();
    const int Nz = mesh.getNz();
    const double sectionX = mesh.getLx() * 0.5;
    const double sectionZ = mesh.getLz() * 0.5;

    const bool anySection = showSectionX || showSectionZ;

    std::size_t elementIndex = 0;
    for (int ez = 0; ez < Nz; ++ez)
    {
        for (int ey = 0; ey < Ny; ++ey)
        {
            for (int ex = 0; ex < Nx; ++ex)
            {
                const Element& element = mesh.getElements()[elementIndex++];

                if (anySection)
                {
                    bool drawElement = false;
                    if (showSectionX && elementCrossesPlaneX(element, mesh, sectionX))
                        drawElement = true;
                    if (showSectionZ && elementCrossesPlaneZ(element, mesh, sectionZ))
                        drawElement = true;

                    if (drawElement)
                        addElementAllFaces(element, indices_);
                }
                else
                {
                    if (ez == 0)
                        addElementFaceIndices(element, 0, indices_);
                    if (ez == Nz - 1)
                        addElementFaceIndices(element, 1, indices_);
                    if (ex == 0)
                        addElementFaceIndices(element, 2, indices_);
                    if (ex == Nx - 1)
                        addElementFaceIndices(element, 3, indices_);
                    if (ey == 0)
                        addElementFaceIndices(element, 4, indices_);
                    if (ey == Ny - 1)
                        addElementFaceIndices(element, 5, indices_);
                }
            }
        }
    }

    uploadIndices();
}

void MeshRenderer::uploadIndices()
{
    indexCount_ = static_cast<int>(indices_.size());

    glBindVertexArray(vao_);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(indices_.size() * sizeof(unsigned int)),
                 indices_.empty() ? nullptr : indices_.data(),
                 GL_DYNAMIC_DRAW);
    glBindVertexArray(0);
}

void MeshRenderer::render(const float scaleFactor,
                          const float aspectRatio,
                          const float cameraDist,
                          const float camAngleX,
                          const float camAngleY,
                          const float camPanX,
                          const float camPanY) const
{
    if (vao_ == 0 || indexCount_ == 0)
        return;

    float projection[16]{};
    float view[16]{};
    float model[16]{};

    const float zNear = 0.1f;
    const float zFar = (std::max)(cameraDist + modelExtent_ * 4.0f, 100.0f);
    makePerspective(45.0f * 3.14159265f / 180.0f, aspectRatio, zNear, zFar, projection);

    // View = T(pan) * T(0,0,-dist) — віддалення та зсув у площині екрана.
    float viewDist[16]{};
    float viewPan[16]{};
    makeTranslate(0.0f, 0.0f, -cameraDist, viewDist);
    makeTranslate(camPanX, camPanY, 0.0f, viewPan);
    multiplyMat4(viewPan, viewDist, view);

    // Model = T(center) * R_y * R_x * T(-center)
    float toCenter[16]{};
    float fromCenter[16]{};
    float rotX[16]{};
    float rotY[16]{};
    float rotYrotX[16]{};
    float rotAtOrigin[16]{};

    makeTranslate(modelCenterX_, modelCenterY_, modelCenterZ_, toCenter);
    makeTranslate(-modelCenterX_, -modelCenterY_, -modelCenterZ_, fromCenter);
    makeRotateX(camAngleX, rotX);
    makeRotateY(camAngleY, rotY);

    multiplyMat4(rotY, rotX, rotYrotX);
    multiplyMat4(rotYrotX, fromCenter, rotAtOrigin);
    multiplyMat4(toCenter, rotAtOrigin, model);

    glUseProgram(shaderProgram_);
    glUniformMatrix4fv(locProjection_, 1, GL_FALSE, projection);
    glUniformMatrix4fv(locView_, 1, GL_FALSE, view);
    glUniformMatrix4fv(locModel_, 1, GL_FALSE, model);
    glUniform1f(locScaleFactor_, scaleFactor);
    glUniform1f(locMinStress_, minStress_);
    glUniform1f(locMaxStress_, maxStress_);

    glBindVertexArray(vao_);

    // 1) Кольорова заливка (поля напружень)
    glUniform1i(locIsWireframe_, 0);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0f, 1.0f);
    glDrawElements(GL_TRIANGLES, indexCount_, GL_UNSIGNED_INT, nullptr);

    // 2) Каркас сітки поверх заливки
    glDisable(GL_POLYGON_OFFSET_FILL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glUniform1i(locIsWireframe_, 1);
    glDrawElements(GL_TRIANGLES, indexCount_, GL_UNSIGNED_INT, nullptr);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glBindVertexArray(0);
}

void MeshRenderer::shutdown()
{
    if (ebo_ != 0)
    {
        glDeleteBuffers(1, &ebo_);
        ebo_ = 0;
    }
    if (vbo_ != 0)
    {
        glDeleteBuffers(1, &vbo_);
        vbo_ = 0;
    }
    if (vao_ != 0)
    {
        glDeleteVertexArrays(1, &vao_);
        vao_ = 0;
    }
    if (shaderProgram_ != 0)
    {
        glDeleteProgram(shaderProgram_);
        shaderProgram_ = 0;
    }
    vertices_.clear();
    indices_.clear();
    indexCount_ = 0;
}

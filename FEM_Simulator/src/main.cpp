#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "GlobalSystem.h"
#include "MathCore.h"
#include "MeshRenderer.h"
#include "StressAnalysis.h"

int main()
{
    if (!glfwInit())
        return -1;

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    GLFWwindow* window = glfwCreateWindow(1280, 720, u8"МСЕ Симулятор - Аналіз Міцності", nullptr, nullptr);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        glfwDestroyWindow(window);
        glfwTerminate();
        return -1;
    }

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    ImGui::StyleColorsDark();

#ifdef _WIN32
    ImFontConfig fontConfig;
    fontConfig.OversampleH = 2;
    fontConfig.OversampleV = 2;
    io.Fonts->AddFontFromFileTTF(
        "C:\\Windows\\Fonts\\segoeui.ttf",
        18.0f,
        &fontConfig,
        io.Fonts->GetGlyphRangesCyrillic());
#endif

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");

    Mesh mesh;
    mesh.generateRectangularParallelepiped(10.0, 5.0, 3.0, 4, 2, 3);

    MaterialProperties material;
    material.youngModulus = 2.0e11;
    material.poissonRatio = 0.3;

    GlobalSystem system;
    system.buildFromMesh(mesh, material, 1.0e7);

    StressAnalyzer stressAnalyzer;
    stressAnalyzer.compute(mesh, system, material);

    MeshRenderer renderer;
    renderer.build(mesh, system, stressAnalyzer);

    float scaleFactor = 100.0f;
    bool showSectionX = false;
    bool showSectionZ = false;
    bool prevSectionX = false;
    bool prevSectionZ = false;

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        if (showSectionX != prevSectionX || showSectionZ != prevSectionZ)
        {
            renderer.updateIndices(mesh, showSectionX, showSectionZ);
            prevSectionX = showSectionX;
            prevSectionZ = showSectionZ;
        }

        int width = 0;
        int height = 0;
        glfwGetFramebufferSize(window, &width, &height);
        const float aspect = height > 0 ? static_cast<float>(width) / static_cast<float>(height) : 1.0f;

        glViewport(0, 0, width, height);
        glEnable(GL_DEPTH_TEST);
        glClearColor(0.12f, 0.12f, 0.12f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        renderer.render(scaleFactor, aspect);

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::Begin(u8"Результати МСЕ");
        ImGui::Text(u8"Вузлів: %zu", mesh.nodeCount());
        ImGui::Text(u8"Елементів: %zu", mesh.elementCount());
        ImGui::Text(u8"DOF: %d", system.dofCount());
        ImGui::Text(u8"Півширина L: %d", system.halfBandwidth());
        ImGui::Text(u8"Orphan-вузлів: %d", system.orphanNodeCount());
        ImGui::Separator();
        ImGui::Text(u8"max |U_y| = %.6e м", system.maxDisplacementY());
        ImGui::Text(u8"max σ₁ = %.6e Па", stressAnalyzer.maxPrincipalStress());
        ImGui::SliderFloat(u8"Масштаб деформації", &scaleFactor, 1.0f, 5000.0f);
        ImGui::Separator();
        ImGui::Checkbox(u8"Показати січення X", &showSectionX);
        ImGui::Checkbox(u8"Показати січення Z", &showSectionZ);
        ImGui::End();

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    renderer.shutdown();
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}

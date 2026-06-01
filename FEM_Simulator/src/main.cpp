#include <glad/glad.h>

#include <GLFW/glfw3.h>



#include "imgui.h"

#include "imgui_impl_glfw.h"

#include "imgui_impl_opengl3.h"



#include "GlobalSystem.h"

#include "MathCore.h"

#include "MeshRenderer.h"

#include "StressAnalysis.h"



#include <algorithm>
#include <chrono>

#include <fstream>

#include <future>

#include <iomanip>

#include <memory>

#include <stdexcept>

#include <string>

#include <utility>

#include <vector>



struct SimulationParams

{

    float Lx = 10.0f;

    float Ly = 5.0f;

    float Lz = 3.0f;



    int Nx = 4;

    int Ny = 2;

    int Nz = 3;



    float youngModulus = 2.0e11f;

    float poissonRatio = 0.3f;



    float pressure = 1.0e7f; // P — тиск на верхній грані Y = Ly

};



struct SimulationContext

{

    Mesh mesh;

    GlobalSystem system;

    StressAnalyzer stressAnalyzer;

    std::vector<int> validNodeIds;

    bool hasResults = false;

    std::string lastError;

};



struct SimulationComputeResult

{

    bool success = false;

    std::string error;

    SimulationContext context;

};



namespace

{

void buildValidNodeIds(const StressAnalyzer& stressAnalyzer, std::vector<int>& validNodeIds)

{

    validNodeIds.clear();

    const auto& results = stressAnalyzer.nodeResults();

    validNodeIds.reserve(results.size());



    for (std::size_t nodeId = 0; nodeId < results.size(); ++nodeId)

    {

        if (results[nodeId].valid)

            validNodeIds.push_back(static_cast<int>(nodeId));

    }

}



bool exportResultsToCsv(const SimulationContext& ctx)

{

    if (!ctx.hasResults)

        return false;



    std::ofstream file("results.csv");

    if (!file)

        return false;



    file << std::setprecision(10);

    file << "ID,X,Y,Z,U_x,U_y,U_z,Sigma_1\n";



    const std::vector<Node>& nodes = ctx.mesh.getNodes();

    const std::vector<double>& u = ctx.system.displacements();

    const auto& stressResults = ctx.stressAnalyzer.nodeResults();



    for (const int nodeId : ctx.validNodeIds)

    {

        const std::size_t id = static_cast<std::size_t>(nodeId);

        const Node& node = nodes[id];



        file << nodeId << ','

             << node.x << ','

             << node.y << ','

             << node.z << ','

             << u[id * 3 + 0] << ','

             << u[id * 3 + 1] << ','

             << u[id * 3 + 2] << ','

             << stressResults[id].principalMax << '\n';

    }



    return true;

}



bool validateParams(const SimulationParams& params, std::string& error)

{

    if (params.Lx <= 0.0f || params.Ly <= 0.0f || params.Lz <= 0.0f)

    {

        error = u8"Lx, Ly, Lz must be positive";

        return false;

    }

    if (params.Nx < 1 || params.Ny < 1 || params.Nz < 1)

    {

        error = u8"Nx, Ny, Nz must be >= 1";

        return false;

    }

    if (params.youngModulus <= 0.0f)

    {

        error = u8"E must be positive";

        return false;

    }

    if (params.poissonRatio <= -1.0f || params.poissonRatio >= 0.5f)

    {

        error = u8"nu must be in (-1, 0.5)";

        return false;

    }

    return true;

}



SimulationComputeResult runSimulationCompute(const SimulationParams& params)

{

    SimulationComputeResult result{};

    resetCalculationProgress();



    std::string validationError;

    if (!validateParams(params, validationError))

    {

        result.success = false;

        result.error = validationError;

        return result;

    }



    try

    {

        result.context.mesh.clear();

        result.context.mesh.generateRectangularParallelepiped(

            static_cast<double>(params.Lx),

            static_cast<double>(params.Ly),

            static_cast<double>(params.Lz),

            params.Nx,

            params.Ny,

            params.Nz);



        MaterialProperties material;

        material.youngModulus = static_cast<double>(params.youngModulus);

        material.poissonRatio = static_cast<double>(params.poissonRatio);



        result.context.system.buildFromMesh(

            result.context.mesh, material, static_cast<double>(params.pressure));



        result.context.stressAnalyzer.compute(result.context.mesh, result.context.system, material);

        buildValidNodeIds(result.context.stressAnalyzer, result.context.validNodeIds);

        result.context.hasResults = true;
        result.context.lastError.clear();
        result.success = true;
        setCalculationProgress(1.0f);
        return result;

    }

    catch (const std::exception& ex)

    {

        result.success = false;

        result.error = ex.what();

        return result;

    }

    catch (...)

    {

        result.success = false;

        result.error = u8"Unknown simulation error";

        return result;

    }

}



void applyRendererFromContext(SimulationContext& ctx,

                              MeshRenderer& renderer,

                              bool showSectionX,

                              bool showSectionZ)

{

    renderer.build(ctx.mesh, ctx.system, ctx.stressAnalyzer);

    renderer.updateIndices(ctx.mesh, showSectionX, showSectionZ);

}

void applyAutoFitCamera(const SimulationParams& params,
                        float& cameraDist,
                        float& camPanX,
                        float& camPanY)
{
    const float maxDim = (std::max)(params.Lx, (std::max)(params.Ly, params.Lz));
    cameraDist = maxDim * 2.5f;
    camPanX = -params.Lx * 0.2f;
    camPanY = 0.0f;
}

struct AsyncSimulationState

{

    bool running = false;

    std::future<SimulationComputeResult> future;

};



void startSimulationAsync(AsyncSimulationState& asyncState, const SimulationParams& params)

{

    if (asyncState.running)

        return;



    resetCalculationProgress();

    asyncState.running = true;

    asyncState.future = std::async(std::launch::async, [params]() {

        return runSimulationCompute(params);

    });

}



bool pollSimulationAsync(AsyncSimulationState& asyncState,

                         SimulationContext& ctx,

                         MeshRenderer& renderer,

                         bool showSectionX,

                         bool showSectionZ)

{

    if (!asyncState.running || !asyncState.future.valid())

        return false;



    if (asyncState.future.wait_for(std::chrono::seconds(0)) != std::future_status::ready)

        return false;



    SimulationComputeResult result = asyncState.future.get();

    asyncState.running = false;

    asyncState.future = {};



    if (result.success)

    {

        ctx = std::move(result.context);

        applyRendererFromContext(ctx, renderer, showSectionX, showSectionZ);
        setCalculationProgress(1.0f);
        return true;

    }



    ctx.hasResults = false;

    ctx.lastError = result.error;

    return false;

}



void drawCalculationProgress(bool computeRunning)

{

    if (!computeRunning)

        return;



    const float progress = getCalculationProgress();

    ImGui::ProgressBar(progress, ImVec2(-1.0f, 0.0f));

    ImGui::Text(u8"Обчислення: %.0f%%", progress * 100.0f);

}



void drawSettingsWindow(SimulationParams& params, bool& runRequested, bool computeRunning)

{

    ImGui::Begin(u8"Налаштування");



    drawCalculationProgress(computeRunning);



    ImGui::Text(u8"Геометрія (м)");

    ImGui::InputFloat(u8"L_x", &params.Lx);

    ImGui::InputFloat(u8"L_y", &params.Ly);

    ImGui::InputFloat(u8"L_z", &params.Lz);



    ImGui::Separator();

    ImGui::Text(u8"Сітка (СЕ)");

    ImGui::InputInt(u8"N_x", &params.Nx);

    ImGui::InputInt(u8"N_y", &params.Ny);

    ImGui::InputInt(u8"N_z", &params.Nz);

    if (params.Nx < 1) params.Nx = 1;

    if (params.Ny < 1) params.Ny = 1;

    if (params.Nz < 1) params.Nz = 1;



    ImGui::Separator();

    ImGui::Text(u8"Матеріал");

    ImGui::InputFloat(u8"E (Па)", &params.youngModulus);

    ImGui::InputFloat(u8"nu", &params.poissonRatio);



    ImGui::Separator();

    ImGui::Text(u8"Закріплення");

    {

        static const char* boundaryItems[] = {u8"Нижня грань (Y = 0) — жорстко закріплена"};

        static int boundaryIndex = 0;

        ImGui::BeginDisabled();

        ImGui::Combo(u8"Умови закріплення", &boundaryIndex, boundaryItems, 1);

        ImGui::EndDisabled();

    }



    ImGui::Separator();

    ImGui::Text(u8"Навантаження");

    ImGui::InputFloat(u8"P (Па)", &params.pressure);



    ImGui::Separator();

    ImGui::BeginDisabled(computeRunning);

    if (ImGui::Button(u8"Обчислити (Перерахувати)", ImVec2(-1.0f, 0.0f)))

        runRequested = true;

    ImGui::EndDisabled();



    ImGui::End();

}



void drawResultsWindow(SimulationContext& ctx,

                       float& scaleFactor,

                       bool& showSectionX,

                       bool& showSectionZ,

                       float& cameraDist,

                       float& camAngleX,

                       float& camAngleY,
                       float& camPanX,
                       float& camPanY,
                       bool computeRunning)

{

    ImGui::Begin(u8"Результати МСЕ");



    drawCalculationProgress(computeRunning);



    if (!ctx.hasResults && !computeRunning)

    {

        ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.4f, 1.0f), u8"Розрахунок не виконано");

        if (!ctx.lastError.empty())

            ImGui::TextWrapped("%s", ctx.lastError.c_str());

        ImGui::End();

        return;

    }



    if (computeRunning)

    {

        ImGui::Text(u8"Очікуйте завершення розрахунку...");

        ImGui::End();

        return;

    }



    ImGui::Text(u8"Вузлів: %zu", ctx.mesh.nodeCount());

    ImGui::Text(u8"Елементів: %zu", ctx.mesh.elementCount());

    ImGui::Text(u8"DOF: %d", ctx.system.dofCount());

    ImGui::Text(u8"Півширина L: %d", ctx.system.halfBandwidth());

    ImGui::Text(u8"Orphan-вузлів: %d", ctx.system.orphanNodeCount());

    ImGui::Separator();

    ImGui::Text(u8"max |U_y| = %.6e м", ctx.system.maxDisplacementY());

    ImGui::Text(u8"max σ₁ = %.6e Па", ctx.stressAnalyzer.maxPrincipalStress());



    ImGui::SliderFloat(u8"Масштаб деформації", &scaleFactor, 1.0f, 5000.0f, "%.0f",

                       ImGuiSliderFlags_Logarithmic);



    ImGui::Separator();

    ImGui::Text(u8"Камера");

    ImGui::SliderFloat(u8"Віддалення (Zoom)", &cameraDist, 5.0f, 300.0f);

    ImGui::SliderFloat(u8"Обертання X", &camAngleX, -3.14f, 3.14f);

    ImGui::SliderFloat(u8"Обертання Y", &camAngleY, -3.14f, 3.14f);

    ImGui::SliderFloat(u8"Зсув X (Pan)", &camPanX, -50.0f, 50.0f);

    ImGui::SliderFloat(u8"Зсув Y (Pan)", &camPanY, -50.0f, 50.0f);



    ImGui::Separator();

    ImGui::Checkbox(u8"Показати січення X", &showSectionX);

    ImGui::Checkbox(u8"Показати січення Z", &showSectionZ);

    ImGui::End();

}



void drawTableWindow(SimulationContext& ctx, bool computeRunning)

{

    ImGui::Begin(u8"Таблиця результатів");



    if (computeRunning)

    {

        ImGui::Text(u8"Таблиця буде доступна після завершення розрахунку");

        ImGui::End();

        return;

    }



    if (!ctx.hasResults)

    {

        ImGui::Text(u8"Натисніть «Обчислити» у вікні налаштувань");

        ImGui::End();

        return;

    }



    if (ImGui::Button(u8"Експорт у CSV"))

    {

        if (exportResultsToCsv(ctx))

            ImGui::OpenPopup(u8"ExportOk");

        else

            ImGui::OpenPopup(u8"ExportFail");

    }



    if (ImGui::BeginPopup(u8"ExportOk"))

    {

        ImGui::Text(u8"Збережено: results.csv");

        ImGui::EndPopup();

    }

    if (ImGui::BeginPopup(u8"ExportFail"))

    {

        ImGui::Text(u8"Помилка запису results.csv");

        ImGui::EndPopup();

    }



    ImGui::Text(u8"Валідних вузлів: %d", static_cast<int>(ctx.validNodeIds.size()));



    const ImGuiTableFlags tableFlags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg

                                      | ImGuiTableFlags_Resizable | ImGuiTableFlags_ScrollY;

    if (ImGui::BeginTable(u8"ResultsTable", 8, tableFlags, ImVec2(0.0f, 400.0f)))

    {

        ImGui::TableSetupColumn(u8"ID");

        ImGui::TableSetupColumn(u8"X");

        ImGui::TableSetupColumn(u8"Y");

        ImGui::TableSetupColumn(u8"Z");

        ImGui::TableSetupColumn(u8"U_x");

        ImGui::TableSetupColumn(u8"U_y");

        ImGui::TableSetupColumn(u8"U_z");

        ImGui::TableSetupColumn(u8"Sigma_1");

        ImGui::TableSetupScrollFreeze(0, 1);

        ImGui::TableHeadersRow();



        const std::vector<Node>& nodes = ctx.mesh.getNodes();

        const std::vector<double>& u = ctx.system.displacements();

        const auto& stressResults = ctx.stressAnalyzer.nodeResults();



        ImGuiListClipper clipper;

        clipper.Begin(static_cast<int>(ctx.validNodeIds.size()));

        while (clipper.Step())

        {

            for (int row = clipper.DisplayStart; row < clipper.DisplayEnd; ++row)

            {

                const int nodeId = ctx.validNodeIds[static_cast<std::size_t>(row)];

                const std::size_t id = static_cast<std::size_t>(nodeId);

                const Node& node = nodes[id];



                ImGui::TableNextRow();

                ImGui::TableNextColumn();

                ImGui::Text("%d", nodeId);

                ImGui::TableNextColumn();

                ImGui::Text("%.6f", node.x);

                ImGui::TableNextColumn();

                ImGui::Text("%.6f", node.y);

                ImGui::TableNextColumn();

                ImGui::Text("%.6f", node.z);

                ImGui::TableNextColumn();

                ImGui::Text("%.6e", u[id * 3 + 0]);

                ImGui::TableNextColumn();

                ImGui::Text("%.6e", u[id * 3 + 1]);

                ImGui::TableNextColumn();

                ImGui::Text("%.6e", u[id * 3 + 2]);

                ImGui::TableNextColumn();

                ImGui::Text("%.6e", stressResults[id].principalMax);

            }

        }



        ImGui::EndTable();

    }



    ImGui::End();

}

} // namespace



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



    SimulationParams params;

    SimulationContext simContext;

    MeshRenderer renderer;

    AsyncSimulationState asyncState;



    float scaleFactor = 100.0f;

    float cameraDist = 40.0f;

    float camAngleX = 0.5f;

    float camAngleY = -0.5f;

    float camPanX = 0.0f;

    float camPanY = 0.0f;

    bool showSectionX = false;

    bool showSectionZ = false;

    bool prevSectionX = false;

    bool prevSectionZ = false;



    startSimulationAsync(asyncState, params);



    while (!glfwWindowShouldClose(window))

    {

        glfwPollEvents();



        if (pollSimulationAsync(asyncState, simContext, renderer, showSectionX, showSectionZ))
            applyAutoFitCamera(params, cameraDist, camPanX, camPanY);



        ImGui_ImplOpenGL3_NewFrame();

        ImGui_ImplGlfw_NewFrame();

        ImGui::NewFrame();



        bool runRequested = false;

        drawSettingsWindow(params, runRequested, asyncState.running);



        if (runRequested && !asyncState.running)

            startSimulationAsync(asyncState, params);



        if (showSectionX != prevSectionX || showSectionZ != prevSectionZ)

        {

            if (simContext.hasResults && !asyncState.running)

                renderer.updateIndices(simContext.mesh, showSectionX, showSectionZ);

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



        if (simContext.hasResults && !asyncState.running)

            renderer.render(scaleFactor, aspect, cameraDist, camAngleX, camAngleY, camPanX, camPanY);



        drawResultsWindow(simContext,

                          scaleFactor,

                          showSectionX,

                          showSectionZ,

                          cameraDist,

                          camAngleX,

                          camAngleY,
                          camPanX,
                          camPanY,
                          asyncState.running);

        drawTableWindow(simContext, asyncState.running);



        ImGui::Render();

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());



        glfwSwapBuffers(window);

    }



    if (asyncState.future.valid())

        asyncState.future.wait();



    renderer.shutdown();

    ImGui_ImplOpenGL3_Shutdown();

    ImGui_ImplGlfw_Shutdown();

    ImGui::DestroyContext();



    glfwDestroyWindow(window);

    glfwTerminate();



    return 0;

}


#include <glad/glad.h>

#include <GLFW/glfw3.h>



#include "imgui.h"

#include "imgui_impl_glfw.h"

#include "imgui_impl_opengl3.h"



#include "GlobalSystem.h"

#include "MathCore.h"

#include "MeshRenderer.h"

#include "StressAnalysis.h"

#include "FemValidation.h"



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

    file << "ID,X,Y,Z,U_x,U_y,U_z,Sigma_1,Sigma_3\n";



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

             << stressResults[id].principalMax << ','

             << stressResults[id].principalMin << '\n';

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

        result.context.mesh.buildBoundaryData(static_cast<double>(params.pressure));

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



void rebuildMeshPreview(const SimulationParams& params, Mesh& previewMesh, MeshRenderer& renderer)
{
    previewMesh.clear();
    previewMesh.generateRectangularParallelepiped(
        static_cast<double>(params.Lx),
        static_cast<double>(params.Ly),
        static_cast<double>(params.Lz),
        params.Nx,
        params.Ny,
        params.Nz);
    previewMesh.buildBoundaryData(static_cast<double>(params.pressure));
    renderer.buildMeshPreview(previewMesh);
}

void drawWorkingArraysWindow(const SimulationContext& ctx, const Mesh& previewMesh, bool hasPreview)
{
    ImGui::Begin(u8"Робочі масиви (п.20 viii)");

    const Mesh& mesh = ctx.hasResults ? ctx.mesh : previewMesh;
    const bool hasMesh = ctx.hasResults || hasPreview;

    if (!hasMesh)
    {
        ImGui::Text(u8"Згенеруйте сітку кнопкою «Оновити сітку» або виконайте розрахунок.");
        ImGui::End();
        return;
    }

    const std::size_t nqp = mesh.nodeCount();
    const std::size_t nel = mesh.elementCount();
    ImGui::Text(u8"nqp = %zu, nel = %zu", nqp, nel);
    ImGui::Text(u8"|ZU| = %zu, |ZP| = %zu", mesh.getZU().size(), mesh.getZP().size());

    if (ctx.hasResults)
    {
        ImGui::Text(u8"MG (стрічка): %zu елементів, L = %d",
                    ctx.system.stiffnessBand().size(),
                    ctx.system.halfBandwidth());
        ImGui::Text(u8"|F| = %zu, |U| = %zu",
                    ctx.system.forceVector().size(),
                    ctx.system.displacements().size());
    }

    if (ImGui::CollapsingHeader(u8"AKT (7) — перші 8 вузлів"))
    {
        const std::vector<double>& akt = mesh.getAKT();
        if (ImGui::BeginTable(u8"AktTable", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg))
        {
            ImGui::TableSetupColumn(u8"ID");
            ImGui::TableSetupColumn(u8"X");
            ImGui::TableSetupColumn(u8"Y");
            ImGui::TableSetupColumn(u8"Z");
            ImGui::TableHeadersRow();
            const int rows = static_cast<int>((std::min)(nqp, std::size_t{8}));
            for (int j = 0; j < rows; ++j)
            {
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::Text("%d", j);
                ImGui::TableNextColumn();
                ImGui::Text("%.4f", akt[static_cast<std::size_t>(j)]);
                ImGui::TableNextColumn();
                ImGui::Text("%.4f", akt[nqp + static_cast<std::size_t>(j)]);
                ImGui::TableNextColumn();
                ImGui::Text("%.4f", akt[2 * nqp + static_cast<std::size_t>(j)]);
            }
            ImGui::EndTable();
        }
    }

    if (ImGui::CollapsingHeader(u8"NT (21) — перший СЕ"))
    {
        if (!mesh.getNT().empty())
        {
            const auto& nt0 = mesh.getNT()[0];
            for (int i = 0; i < kHex20NodeCount; ++i)
                ImGui::Text(u8"NT[%d,0] = %d", i, nt0[static_cast<std::size_t>(i)]);
        }
    }

    if (ImGui::CollapsingHeader(u8"ZU (30)"))
    {
        for (std::size_t i = 0; i < mesh.getZU().size() && i < 16; ++i)
            ImGui::Text(u8"ZU[%zu] = %d", i, mesh.getZU()[i]);
    }

    if (ImGui::CollapsingHeader(u8"ZP (31)"))
    {
        for (std::size_t i = 0; i < mesh.getZP().size() && i < 8; ++i)
        {
            const Mesh::ZPEntry& e = mesh.getZP()[i];
            ImGui::Text(u8"ZP[%zu]: elem=%d face=%d P=%.3e", i, e.elementIndex, e.faceIndex, e.pressure);
        }
    }

    if (ctx.hasResults && ImGui::CollapsingHeader(u8"MG, F, U — зріз"))
    {
        const auto& band = ctx.system.stiffnessBand();
        const auto& f = ctx.system.forceVector();
        const auto& u = ctx.system.displacements();
        const int show = static_cast<int>((std::min)(f.size(), std::size_t{12}));
        for (int i = 0; i < show; ++i)
            ImGui::Text(u8"F[%d]=%.3e  U[%d]=%.3e", i, f[static_cast<std::size_t>(i)], i, u[static_cast<std::size_t>(i)]);
        if (!band.empty())
            ImGui::Text(u8"MG[0]=%.3e", band[0]);
    }

    ImGui::End();
}

void drawCalculationProgress(bool computeRunning)

{

    if (!computeRunning)

        return;



    const float progress = getCalculationProgress();

    ImGui::ProgressBar(progress, ImVec2(-1.0f, 0.0f));

    ImGui::Text(u8"Обчислення: %.0f%%", progress * 100.0f);

}



void drawSettingsWindow(SimulationParams& params,
                        bool& runRequested,
                        bool& meshPreviewRequested,
                        bool& validationRequested,
                        bool computeRunning)

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

    if (ImGui::Button(u8"Оновити сітку (тріангуляція)", ImVec2(-1.0f, 0.0f)))
        meshPreviewRequested = true;

    if (ImGui::Button(u8"Обчислити (Перерахувати)", ImVec2(-1.0f, 0.0f)))
        runRequested = true;

    if (ImGui::Button(u8"Перевірки методички (зан. 13)", ImVec2(-1.0f, 0.0f)))
        validationRequested = true;

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

    ImGui::Text(u8"Півширина L (ng): %d", ctx.system.halfBandwidth());

    ImGui::Text(u8"|ZU| = %zu, |ZP| = %zu", ctx.mesh.getZU().size(), ctx.mesh.getZP().size());

    ImGui::Separator();

    ImGui::Text(u8"max |U_y| = %.6e м", ctx.system.maxDisplacementY());

    ImGui::Text(u8"σ₁ (max) = %.6e Па", ctx.stressAnalyzer.maxPrincipalStress());

    ImGui::Text(u8"σ₃ (min) = %.6e Па", ctx.stressAnalyzer.minPrincipalStress());



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

    if (ImGui::BeginTable(u8"ResultsTable", 9, tableFlags, ImVec2(0.0f, 400.0f)))

    {

        ImGui::TableSetupColumn(u8"ID");

        ImGui::TableSetupColumn(u8"X");

        ImGui::TableSetupColumn(u8"Y");

        ImGui::TableSetupColumn(u8"Z");

        ImGui::TableSetupColumn(u8"U_x");

        ImGui::TableSetupColumn(u8"U_y");

        ImGui::TableSetupColumn(u8"U_z");

        ImGui::TableSetupColumn(u8"Sigma_1");

        ImGui::TableSetupColumn(u8"Sigma_3");

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

                ImGui::TableNextColumn();

                ImGui::Text("%.6e", stressResults[id].principalMin);

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

    Mesh previewMesh;
    bool meshPreviewReady = false;
    FemValidationReport validationReport{};
    bool showValidationPopup = false;

    rebuildMeshPreview(params, previewMesh, renderer);
    meshPreviewReady = true;



    while (!glfwWindowShouldClose(window))

    {

        glfwPollEvents();



        if (pollSimulationAsync(asyncState, simContext, renderer, showSectionX, showSectionZ))
            applyAutoFitCamera(params, cameraDist, camPanX, camPanY);



        ImGui_ImplOpenGL3_NewFrame();

        ImGui_ImplGlfw_NewFrame();

        ImGui::NewFrame();



        bool runRequested = false;
        bool meshPreviewRequested = false;
        bool validationRequested = false;

        drawSettingsWindow(params, runRequested, meshPreviewRequested, validationRequested, asyncState.running);

        if (meshPreviewRequested && !asyncState.running)
        {
            rebuildMeshPreview(params, previewMesh, renderer);
            meshPreviewReady = true;
            applyAutoFitCamera(params, cameraDist, camPanX, camPanY);
        }

        if (validationRequested && !asyncState.running)
        {
            validationReport = runFemValidationChecks();
            showValidationPopup = true;
        }

        if (runRequested && !asyncState.running)
            startSimulationAsync(asyncState, params);

        if (showValidationPopup)
            ImGui::OpenPopup(u8"ValidationReport");
        if (ImGui::BeginPopupModal(u8"ValidationReport", nullptr, ImGuiWindowFlags_AlwaysAutoResize))
        {
            ImGui::Text(u8"det(J) у центрі куба 4x4x4 м = %.6f (очікується 8)", validationReport.jacobianUnitCubeDet);
            ImGui::TextColored(
                validationReport.jacobianUnitCubeOk ? ImVec4(0.3f, 1.0f, 0.3f, 1.0f) : ImVec4(1.0f, 0.4f, 0.4f, 1.0f),
                validationReport.jacobianUnitCubeOk ? u8"OK" : u8"FAIL");
            ImGui::Text(u8"Тест СЛАР (K без КУ): F_i=sum_j K_ij, max|U-1| = %.3e", validationReport.solverUnitRhsMaxError);
            ImGui::TextColored(
                validationReport.solverUnitRhsOk ? ImVec4(0.3f, 1.0f, 0.3f, 1.0f) : ImVec4(1.0f, 0.4f, 0.4f, 1.0f),
                validationReport.solverUnitRhsOk ? u8"OK" : u8"FAIL");
            if (ImGui::Button(u8"Закрити"))
            {
                showValidationPopup = false;
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndPopup();
        }



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
        else if (meshPreviewReady && !asyncState.running)
            renderer.render(1.0f, aspect, cameraDist, camAngleX, camAngleY, camPanX, camPanY);



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

        drawWorkingArraysWindow(simContext, previewMesh, meshPreviewReady);

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


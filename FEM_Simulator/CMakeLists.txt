cmake_minimum_required(VERSION 3.15)
project(FEM_Simulator VERSION 1.0 LANGUAGES C CXX)

# Встановлюємо стандарт C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

# --- 1. АВТОМАТИЧНЕ ЗАВАНТАЖЕННЯ GLFW ---
FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG 3.3.8
)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(glfw)

# --- 2. АВТОМАТИЧНЕ ЗАВАНТАЖЕННЯ ImGui ---
FetchContent_Declare(
    imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG master
)
FetchContent_MakeAvailable(imgui)

# --- 3. GLAD (локальна копія для IDE та збірки) ---
set(GLAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/glad)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --- 4. ЗБІРКА ПРОЕКТУ ---
# Збираємо всі твої .cpp файли з папки src
file(GLOB MY_SOURCES "src/*.cpp")

# Знаходимо файли ImGui (основні + бекенд для GLFW/OpenGL)
set(IMGUI_SRC_DIR ${imgui_SOURCE_DIR})
file(GLOB IMGUI_SOURCES
    ${IMGUI_SRC_DIR}/*.cpp
    ${IMGUI_SRC_DIR}/backends/imgui_impl_glfw.cpp
    ${IMGUI_SRC_DIR}/backends/imgui_impl_opengl3.cpp
)

# Створюємо фінальний виконуваний файл
add_executable(FEM_App 
    ${MY_SOURCES}
    ${IMGUI_SOURCES}
    ${GLAD_DIR}/src/glad.c
)

# --- 5. ПІДКЛЮЧЕННЯ ШЛЯХІВ ТА ЛІНКОВКА ---
target_include_directories(FEM_App PRIVATE
    src
    ${GLAD_DIR}/include
    ${IMGUI_SRC_DIR}
    ${IMGUI_SRC_DIR}/backends
)

# Лінкуємо GLFW
target_link_libraries(FEM_App PRIVATE glfw)

# Лінкуємо системний OpenGL
if (APPLE)
    find_library(OPENGL_LIBRARY OpenGL)
    target_link_libraries(FEM_App PRIVATE ${OPENGL_LIBRARY})
elseif (WIN32)
    target_link_libraries(FEM_App PRIVATE opengl32)
    target_compile_options(FEM_App PRIVATE /utf-8)
else()
    find_package(OpenGL REQUIRED)
    target_link_libraries(FEM_App PRIVATE OpenGL::GL)
endif()
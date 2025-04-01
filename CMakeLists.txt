cmake_minimum_required(VERSION 3.21)
project(team_10_2025 CXX)

option(BUILD_GUI "Build Vulkan/ImGui executable" OFF)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose build type: Debug  Release" FORCE)
endif ()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Build app: ${BUILD_GUI}")

message("Build type: " ${CMAKE_BUILD_TYPE})

set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -O3")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_SOURCE_DIR "src")
include_directories(${CMAKE_PROJECT_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}") 

# Dependencies
include(FetchContent)

# fmt library
FetchContent_Declare(fmt GIT_REPOSITORY https://github.com/fmtlib/fmt.git GIT_TAG 10.2.1)
FetchContent_MakeAvailable(fmt)
list(APPEND LIBRARIES fmt)

# glm library
FetchContent_Declare(glm GIT_REPOSITORY https://github.com/g-truc/glm.git GIT_TAG 1.0.1)
FetchContent_MakeAvailable(glm)
list(APPEND LIBRARIES glm::glm)

# HighFive library 
FetchContent_Declare(highfive GIT_REPOSITORY https://github.com/highfive-devs/highfive.git GIT_TAG v3.0.0-beta2)
FetchContent_MakeAvailable(highfive)
list(APPEND LIBRARIES HighFive::HighFive)

# OpenMP library 
find_package(OpenMP REQUIRED)
list(APPEND LIBRARIES OpenMP::OpenMP_CXX)

# CLI executable
file(GLOB_RECURSE SOURCES_CORE "${CMAKE_SOURCE_DIR}/core/*.cpp")
file(GLOB_RECURSE SOURCES_CLI "${CMAKE_SOURCE_DIR}/cli/*.cpp")
add_executable(${CMAKE_PROJECT_NAME} ${SOURCES_CORE} ${SOURCES_CLI})
target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/core)
target_link_libraries(${CMAKE_PROJECT_NAME} ${LIBRARIES})

if (BUILD_GUI)
    # Vulkan
    find_package(Vulkan REQUIRED)
    list(APPEND LIBRARIES Vulkan::Vulkan)

    # GLFW
    option(GLFW_BUILD_DOCS OFF)
    option(GLFW_BUILD_EXAMPLES OFF)
    option(GLFW_BUILD_TESTS OFF)
    option(GLFW_INSTALL OFF)
    FetchContent_Declare(glfw GIT_REPOSITORY https://github.com/glfw/glfw.git GIT_TAG 3.3.9)
    FetchContent_MakeAvailable(glfw)
    set(GLFW_INCLUDE_DIR "${glfw_SOURCE_DIR}/include")
    list(APPEND LIBRARIES glfw)

    # ImGui
    FetchContent_Declare(imgui GIT_REPOSITORY https://github.com/ocornut/imgui.git GIT_TAG v1.91.9b)
    FetchContent_MakeAvailable(imgui)
    file(GLOB IMGUI_SOURCES ${imgui_SOURCE_DIR}/*.cpp)
    set(IMGUI_BACKENDS_DIR ${imgui_SOURCE_DIR}/backends)
    set(IMGUI_BACKEND_SOURCES ${IMGUI_BACKENDS_DIR}/imgui_impl_glfw.cpp ${IMGUI_BACKENDS_DIR}/imgui_impl_vulkan.cpp)
    set(IMGUI_BACKEND_HEADERS ${IMGUI_BACKENDS_DIR}/imgui_impl_glfw.h ${IMGUI_BACKENDS_DIR}/imgui_impl_vulkan.h)
    add_library(imgui STATIC ${IMGUI_SOURCES} ${IMGUI_BACKEND_HEADERS} ${IMGUI_BACKEND_SOURCES})
    target_include_directories(imgui PUBLIC ${IMGUI_BACKENDS_DIR} ${Vulkan_INCLUDE_DIR} ${imgui_SOURCE_DIR} ${GLFW_INCLUDE_DIR})
    set_target_properties(imgui PROPERTIES LINKER_LANGUAGE CXX)
    list(APPEND LIBRARIES imgui)

    # ImPlot
    FetchContent_Declare(implot GIT_REPOSITORY https://github.com/epezent/implot.git)
    FetchContent_MakeAvailable(implot)
    file(GLOB IMPLOT_SOURCES ${implot_SOURCE_DIR}/*.cpp)
    add_library(implot STATIC ${IMPLOT_SOURCES} ${IMPLOT_HEADERS})
    target_include_directories(implot PUBLIC ${implot_SOURCE_DIR} ${imgui_SOURCE_DIR})
    list(APPEND LIBRARIES implot)

    # GUI executable
    file(GLOB_RECURSE SOURCES_GUI  "${CMAKE_SOURCE_DIR}/gui/*.cpp")
    add_executable(${CMAKE_PROJECT_NAME}_app ${SOURCES_CORE} ${SOURCES_GUI})
    target_include_directories(${CMAKE_PROJECT_NAME}_app PUBLIC ${CMAKE_SOURCE_DIR}/core)
    target_link_libraries(${CMAKE_PROJECT_NAME}_app ${LIBRARIES})
endif()

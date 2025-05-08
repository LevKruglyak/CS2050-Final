cmake_minimum_required(VERSION 3.21)
project(lkxpm CXX)

option(BUILD_APP "Build Vulkan/ImGui executable" OFF)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose build type: Debug  Release" FORCE)
endif ()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Build app: ${BUILD_APP}")

message("Build type: " ${CMAKE_BUILD_TYPE})

set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_DEBUG "-ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_SOURCE_DIR "src")
include_directories(${CMAKE_PROJECT_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}") 

# Dependencies
include(FetchContent)

# glm library
FetchContent_Declare(glm GIT_REPOSITORY https://github.com/g-truc/glm.git GIT_TAG 1.0.1)
FetchContent_MakeAvailable(glm)
list(APPEND LIBRARIES glm::glm)
list(APPEND INCLUDES ${glm_SOURCE_DIR})

# json library
FetchContent_Declare(nlohmann_json GIT_REPOSITORY https://github.com/nlohmann/json.git GIT_TAG v3.12.0)
FetchContent_MakeAvailable(nlohmann_json)
list(APPEND LIBRARIES nlohmann_json::nlohmann_json)

# Perlin noise library
FetchContent_Declare(perlin GIT_REPOSITORY https://github.com/Reputeless/PerlinNoise.git GIT_TAG v3.0.0)
FetchContent_MakeAvailable(perlin)
list(APPEND INCLUDES ${perlin_SOURCE_DIR})

# stb library
FetchContent_Declare(stb GIT_REPOSITORY https://github.com/nothings/stb.git GIT_TAG master)
FetchContent_MakeAvailable(stb)
list(APPEND INCLUDES ${stb_SOURCE_DIR})

# OpenMP library 
find_package(OpenMP REQUIRED)
list(APPEND LIBRARIES OpenMP::OpenMP_CXX)

# MPI library
find_package(MPI REQUIRED)
list(APPEND LIBRARIES MPI::MPI_CXX)

# FFTW library
FetchContent_Declare(findfftw GIT_REPOSITORY https://github.com/egpbos/findFFTW.git GIT_TAG master)
FetchContent_MakeAvailable(findfftw)
list(APPEND CMAKE_MODULE_PATH "${findfftw_SOURCE_DIR}")
find_package(FFTW REQUIRED COMPONENTS DOUBLE_LIB DOUBLE_OPENMP_LIB DOUBLE_MPI_LIB)
list(APPEND LIBRARIES ${FFTW_DOUBLE_LIB} ${FFTW_DOUBLE_OPENMP_LIB} ${FFTW_DOUBLE_MPI_LIB})

# CLI executable
file(GLOB_RECURSE SOURCES_CORE "${CMAKE_SOURCE_DIR}/core/*.cpp")
set(SOURCE_MAIN_CLI "${CMAKE_SOURCE_DIR}/main_cli.cpp")

add_executable(${CMAKE_PROJECT_NAME} ${SOURCES_CORE} ${SOURCE_MAIN_CLI})
target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/core ${INCLUDES})
target_link_libraries(${CMAKE_PROJECT_NAME} ${LIBRARIES})

if (BUILD_APP)
    # ImGui Bundle
    FetchContent_Declare(imgui_bundle GIT_REPOSITORY https://github.com/pthom/imgui_bundle.git GIT_TAG main)
    FetchContent_MakeAvailable(imgui_bundle)

    # App executable
    set(SOURCE_MAIN_APP "${CMAKE_SOURCE_DIR}/main_app.cpp")
    imgui_bundle_add_app(${CMAKE_PROJECT_NAME}_app ${SOURCES_CORE} ${SOURCE_MAIN_APP})
    target_link_libraries(${CMAKE_PROJECT_NAME}_app PRIVATE ${LIBRARIES})
    target_include_directories(${CMAKE_PROJECT_NAME}_app PRIVATE ${CMAKE_SOURCE_DIR}/core ${INCLUDES})
endif()

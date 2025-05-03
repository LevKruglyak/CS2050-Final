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
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O3")
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

# small_vec library
FetchContent_Declare(small_vec GIT_REPOSITORY https://github.com/gharveymn/small_vector.git GIT_TAG v0.10.2)
FetchContent_MakeAvailable(small_vec)
list(APPEND LIBRARIES gch::small_vector)

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
    # ImGui Bundle
    FetchContent_Declare(imgui_bundle GIT_REPOSITORY https://github.com/pthom/imgui_bundle.git GIT_TAG main)
    FetchContent_MakeAvailable(imgui_bundle)

    # GUI executable
    file(GLOB_RECURSE SOURCES_GUI  "${CMAKE_SOURCE_DIR}/gui/*.cpp")
    imgui_bundle_add_app(${CMAKE_PROJECT_NAME}_app ${SOURCES_CORE} ${SOURCES_GUI})
    target_link_libraries(${CMAKE_PROJECT_NAME}_app PRIVATE ${LIBRARIES})
    target_include_directories(${CMAKE_PROJECT_NAME}_app PRIVATE ${CMAKE_SOURCE_DIR}/core ${glm_SOURCE_DIR} ${legit_profiler_SOURCE_DIR})
endif()

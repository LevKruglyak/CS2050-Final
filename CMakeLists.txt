cmake_minimum_required(VERSION 3.21)
project(team_10_2025 CXX)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose build type: Debug  Release" FORCE)
endif ()

message("Build type: " ${CMAKE_BUILD_TYPE})

set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -O3")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_SOURCE_DIR "src")

include(FetchContent)

# fmt library
FetchContent_Declare(fmt GIT_REPOSITORY https://github.com/fmtlib/fmt.git GIT_TAG 10.2.1)
FetchContent_MakeAvailable(fmt)

list(APPEND LIBRARIES fmt)

# Executable
file(GLOB_RECURSE SOURCES "${CMAKE_SOURCE_DIR}/*.cpp")
add_executable(${CMAKE_PROJECT_NAME} ${SOURCES})
target_link_libraries(${CMAKE_PROJECT_NAME} ${LIBRARIES})

cmake_minimum_required(VERSION 3.5.2)

project(catch_builder CXX)
include(ExternalProject)
find_package(Git REQUIRED)

PROJECT(${GLOBAL_PROJECT_NAME})
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fdiagnostics-color=always")
include_directories("include")
add_subdirectory(sources)
include_directories(${CATCH_INCLUDE_DIR})
enable_testing(true)
add_subdirectory(tests)

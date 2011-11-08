## Minimum required version of CMake to configure the project
cmake_minimum_required (VERSION 2.8)

## Name of the project
project (k3match)

## Adjust the module path 
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

## Project release version
set (K3MATCH_VERSION_MAJOR 0)
set (K3MATCH_VERSION_MINOR 0)
set (K3MATCH_VERSION_PATCH 1)
set (K3MATCH_VERSION "${K3MATCH_VERSION_MAJOR}.${K3MATCH_VERSION_MINOR}.${K3MATCH_VERSION_PATCH}")

## Add definitions
add_definitions (-O2)
add_definitions (-fPIC)
add_definitions (-Wall)
add_definitions (-Werror)
add_definitions (-g)

## Options
option (K3MATCH_SHARED_LIBRARY "Build k3match as shared library." ON)
option (K3MATCH_SINGLE_PRECISION "Use single precision for all floating point calculations." OFF)
option (K3MATCH_WITH_PYTHON_BINDINGS "Build Python bindings." ON)

## Set include directories
include_directories (BEFORE ${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR})

## Configuration
configure_file (
    ${PROJECT_SOURCE_DIR}/k3match_config.h.in
    ${PROJECT_BINARY_DIR}/k3match_config.h
)

## Install header files
install (FILES ${PROJECT_SOURCE_DIR}/k3match.h ${PROJECT_BINARY_DIR}/k3match_config.h
  DESTINATION include
  )

## Library
add_subdirectory (k3match)

## Bindings
add_subdirectory (python)

## Test
add_subdirectory (test)


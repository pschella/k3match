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
add_definitions (-O3)
add_definitions (-Wall)
add_definitions (-Werror)

## Options
option (K3MATCH_SHARED_LIBRARY "Build k3match as shared library." ON)

## Set include directories
include_directories (BEFORE ${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR})

## Install header files
install (FILES ${PROJECT_SOURCE_DIR}/k3match.h ${PROJECT_BINARY_DIR}/k3match_config.h
  DESTINATION include
  )

## Library
add_subdirectory (k3match)

## Test
add_subdirectory (test)


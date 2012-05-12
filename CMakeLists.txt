## Minimum required version of CMake to configure the project
cmake_minimum_required (VERSION 2.8)

## Name of the project
project (k3match)

## Adjust the module path 
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

## Project release version
set (K3MATCH_VERSION_MAJOR 1)
set (K3MATCH_VERSION_MINOR 1)
set (K3MATCH_VERSION "${K3MATCH_VERSION_MAJOR}.${K3MATCH_VERSION_MINOR}")

## Add definitions
add_definitions (-O3)
add_definitions (-Wall)
add_definitions (-Werror)

## Options
option (K3MATCH_SHARED_LIBRARY "Build k3match as shared library." ON)
option (K3MATCH_BUILD_TESTS "Build k3match test programs." OFF)

## Set include directories
include_directories (BEFORE ${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR})

## Install header files
install (FILES ${PROJECT_SOURCE_DIR}/k3match.h
  DESTINATION include
  )

## Library
add_subdirectory (k3match)

## Test
if (K3MATCH_BUILD_TESTS)
 add_subdirectory (test)
endif (K3MATCH_BUILD_TESTS)


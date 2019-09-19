# Project: Matrixgen
# Author: Stanislaw Hüll
#
# This file contains the basic CMake setup of the Matrixgen project. Targets
# are defined in 'src/CMakeLists.txt'.

## 1. Prologue
## ---------------------------------
cmake_minimum_required(VERSION 3.14)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

## 2. Project & Subdirectories
## --------------------------
project(Matrixgen
  VERSION 0.0.1
  LANGUAGES CXX)
set(CMAKE_CXX_EXTENTIONS OFF)

add_subdirectory(src)

option(BUILD_TESTS "Build tests" ON)
if(BUILD_TESTS)
  add_subdirectory(unittests)
endif()

option(BUILD_EXAMPLES "Build examples" ON)
if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

## 3. Install
## -----------------------------
find_package(CMakeshift 3.4 REQUIRED)
include(CMakeshift/InstallSimpleLibraryPackage)
cmakeshift_install_simple_library_package(
  EXPORT MatrixgenTargets
  VERSION_COMPATIBILITY SameMajorVersion
  CONFIG_TEMPLATE "cmake/MatrixgenConfig.cmake.in")

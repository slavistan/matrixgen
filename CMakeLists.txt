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
  VERSION 1.0.0
  LANGUAGES CXX)

add_subdirectory(src)

option(BUILD_TESTS "Build tests" OFF)
if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(unittests)
endif()

option(BUILD_EXAMPLES "Build examples" OFF)
if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

## 3. Install
## -----------------------------
find_package(CMakeshift 3.7 REQUIRED)
include(CMakeshift/InstallBasicPackageFiles)
cmakeshift_install_basic_package_files(Matrixgen
  EXPORT MatrixgenTargets
  COMPATIBILITY SameMajorVersion
  DEPENDENCIES
    "CMakeshift 3.7"
    "gsl-lite 0.35"
    "Eigen3 3.3"
    "TBB"
  )

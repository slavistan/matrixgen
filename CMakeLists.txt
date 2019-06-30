# Project: ASCMatrixgen
# Author: Stanislaw Hüll
#
# This file contains the basic CMake setup of the ASCMatrixgen project. Targets
# are defined in 'src/CMakeLists.txt'.

## 1. Prologue
## ---------------------------------
cmake_minimum_required(VERSION 3.12)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # generate a compilation database
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

## 2. Project & Subdirectories
## --------------------------
project(ASCMatrixgen
  VERSION 0.0.1
  LANGUAGES CXX)

add_subdirectory(src)

option(BUILD_TESTS "Build tests" ON)
if(BUILD_TESTS)
    add_subdirectory(test)
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
  EXPORT ASCMatrixgenTargets
  VERSION_COMPATIBILITY SameMajorVersion
  CONFIG_TEMPLATE "cmake/ASCMatrixgenConfig.cmake.in")

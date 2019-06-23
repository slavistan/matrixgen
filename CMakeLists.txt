cmake_minimum_required(VERSION 3.12)

project(matrixgen VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 20)            # set target's default C++ standard
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # generate a compilation database

#############
# Project Targets
############

find_package(Eigen3)

add_library(matrixgen INTERFACE)
target_include_directories(matrixgen INTERFACE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(matrixgen INTERFACE
  Eigen3::Eigen
)


# Miscellaneous targets related to this project but not required for its usage (such as unit-tests).
option(BUILD_TESTS "Build tests" ON)
if(BUILD_TESTS)
    add_subdirectory(test)
endif()

# include external tools
include(cmake/external-tools.cmake)

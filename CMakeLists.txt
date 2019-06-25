cmake_minimum_required(VERSION 3.12)

project(matrixgen VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 20)            # set target's default C++ standard
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # generate a compilation database

#############
# Project Targets
############

find_package(Eigen3)

add_library(matrixgen INTERFACE)
target_include_directories(matrixgen 
  INTERFACE
    "${PROJECT_SOURCE_DIR}/include"
)
target_compile_features(matrixgen INTERFACE cxx_std_20)
target_link_libraries(matrixgen INTERFACE
  Eigen3::Eigen
)

option(BUILD_TESTS "Build tests" ON)
if(BUILD_TESTS)
    add_subdirectory(test)
endif()

option(BUILD_EXAMPLES "Build examples" ON)
if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# define alias so we can use the qualified name everywhere
add_library(${PROJECT_NAME}::matrixgen ALIAS matrixgen)

install(
    DIRECTORY "${PROJECT_SOURCE_DIR}/include/"
    DESTINATION "include")

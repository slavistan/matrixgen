## Project: Matrixgen
## Target: matrixgen
## Author: Stanislaw Hüll

cmake_minimum_required(VERSION 3.14)

# Find dependencies
find_package(CMakeshift 3.4 REQUIRED)
find_package(Eigen3)

add_library(matrixgen INTERFACE)

# Target directories
target_include_directories(matrixgen 
  INTERFACE
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>")

# Target properties
target_compile_features(matrixgen
  INTERFACE
    cxx_std_20)

# Target dependencies
target_link_libraries(matrixgen
  INTERFACE
    Eigen3::Eigen
)


# Install header files and target
install(
    DIRECTORY "${PROJECT_SOURCE_DIR}/include/"
    DESTINATION "include")
install(
    TARGETS matrixgen
    EXPORT MatrixgenTargets
    ARCHIVE DESTINATION "lib"
    LIBRARY DESTINATION "lib"
    RUNTIME DESTINATION "bin"
    INCLUDES DESTINATION "include")

## Project: Matrixgen
## Target: matrixgen
## Author: Stanislaw Hüll

# TODO: 'Compile' headers into build directory 'build/generated/' and install
#       them from there. (I.e. introduce a compilation phase for the headers)
# TODO: Implement a version comment into all header files (and maybe generate
#       a version #define).

cmake_minimum_required(VERSION 3.14)

find_package(CMakeshift 3.4 REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(matrixgen INTERFACE)

target_include_directories(matrixgen 
  INTERFACE
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>")

target_compile_features(matrixgen
  INTERFACE
    cxx_std_20)

target_link_libraries(matrixgen
  INTERFACE
    Eigen3::Eigen
)

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

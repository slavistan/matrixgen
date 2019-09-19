## Project: Matrixgen
## Target: matrixgen
## Author: Stanislaw Hüll

cmake_minimum_required(VERSION 3.14)

find_package( CMakeshift 3.4 REQUIRED )
find_package( MS-GSL REQUIRED )
find_package( Eigen3 3.3 REQUIRED )
# TODO: Debug-mode does not link without TBB. Why?
#       The source code never even mentions TBB :(
find_package( TBB REQUIRED )

add_library(matrixgen INTERFACE)

# Configure MS-GSL's contracts. This macro is required to generate a proper
# runtime error message for a contract violation. Otherwise a simple
# std::terminate will be reported without context information.
# TODO: Make sure that absolutely no code is generated in release mode.
#       The MS-GSL's implementation suggests that something is always left
#       irrespective of any macro.
target_compile_options(matrixgen
  INTERFACE
  $<$<CONFIG:Debug>:-DGSL_THROW_ON_CONTRACT_VIOLATION>)

target_include_directories(matrixgen 
  INTERFACE
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>")

target_compile_features(matrixgen
  INTERFACE
    cxx_std_20)

target_link_libraries(matrixgen
  INTERFACE
    Eigen3::Eigen
    MS-GSL::MS-GSL
    TBB::tbb
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

## Project: Matrixgen
## Target: matrixgen
## Author: Stanislaw Hüll

cmake_minimum_required(VERSION 3.14)

find_package( CMakeshift 3.4 REQUIRED )
find_package( gsl-lite 0.35 REQUIRED )
find_package( Eigen3 3.3 REQUIRED )
find_package( TBB REQUIRED ) # TBB is required for gcc's C++20 features

add_library(matrixgen INTERFACE)

# Disable run-time contract checking for Release builds
target_compile_options(matrixgen
  INTERFACE
    $<$<CONFIG:Release>:-Dgsl_CONFIG_CONTRACT_CHECKING_OFF>)

target_include_directories(matrixgen 
  INTERFACE
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>")

target_compile_features(matrixgen
  INTERFACE
    cxx_std_20)

target_link_libraries(matrixgen
  INTERFACE
    Eigen3::Eigen
    gsl::gsl-lite
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

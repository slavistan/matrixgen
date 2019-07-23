# Project: Matrixgen
# Target: unittests
# Author: Stanislaw Hüll

find_package(CMakeshift 3.4 REQUIRED)
find_package(Catch2 REQUIRED)

add_executable(unittests unittests.cpp)
target_link_libraries(unittests
  PRIVATE
    matrixgen
    Catch2::Catch2)

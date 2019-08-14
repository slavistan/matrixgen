# Project: Matrixgen
# Target:  unittests
# Author:  Stanislaw Hüll

find_package(CMakeshift 3.4   REQUIRED)
find_package(doctest    2.3.4 REQUIRED)

add_executable(unittests unittests.cpp)
target_link_libraries(unittests
  PRIVATE
    matrixgen
    doctest::doctest)

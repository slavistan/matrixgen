# Project: Matrixgen
# Target: example[1-n]
# Author: Stanislaw Hüll

cmake_minimum_required(VERSION 3.10)

add_executable(example1 example1.cpp)
target_link_libraries(example1
  PRIVATE
    matrixgen)

add_executable(example2 example2.cpp)
target_link_libraries(example2
  PRIVATE
    matrixgen)

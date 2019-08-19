# Project: Matrixgen
# Target: example[1-n]
# Author: Stanislaw Hüll

cmake_minimum_required(VERSION 3.10)

add_executable(ex1-adjmat ex1-adjmat.cpp)
target_link_libraries(ex1-adjmat
  PRIVATE
    matrixgen)
#
# add_executable(example2 example2.cpp)
# target_link_libraries(example2
#   PRIVATE
#     matrixgen)
#
# add_executable(example3 example3.cpp)
# target_link_libraries(example3
#   PRIVATE
#     matrixgen)

add_executable(ex4-adjmat ex4-adjmat.cpp)
target_link_libraries(ex4-adjmat
  PRIVATE
    matrixgen)

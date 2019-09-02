# Project: Matrixgen
# Target: ex[1-n]-adjmat
# Author: Stanislaw Hüll

cmake_minimum_required(VERSION 3.10)

add_executable(ex1-adjmat ex1-adjmat.cpp)
target_link_libraries(ex1-adjmat
  PRIVATE
    matrixgen)

add_executable(ex2-adjmat ex2-adjmat.cpp)
target_link_libraries(ex2-adjmat
  PRIVATE
    matrixgen)

add_executable(ex3-adjmat ex3-adjmat.cpp)
target_link_libraries(ex3-adjmat
  PRIVATE
    matrixgen)

add_executable(ex4-adjmat ex4-adjmat.cpp)
target_link_libraries(ex4-adjmat
  PRIVATE
    matrixgen)

add_executable(ex1-assemble ex1-assemble.cpp)
target_link_libraries(ex1-assemble
  PRIVATE
    matrixgen)

add_executable(ex2-assemble ex2-assemble.cpp)
target_link_libraries(ex2-assemble
  PRIVATE
    matrixgen)

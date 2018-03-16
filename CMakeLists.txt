cmake_minimum_required(VERSION 2.8.3)
project(matrix_math_comparison)

set (CMAKE_CXX_STANDARD 11)

#find_package(GSL)


include_directories(include
  lib/Eigen
  lib/TNT
  lib/TNT/TNT
  lib/TNT/JAMA
  lib/GSL
)

#add_executable(matrix_test
  #src/main.cpp
#)

add_executable(eigen_decomp
  src/eigen_decomp.cpp
)

#target_link_libraries(matrix_test gsl gslcblas)

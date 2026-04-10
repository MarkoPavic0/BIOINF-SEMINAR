cmake_minimum_required(VERSION 3.10)
project(bsc_alignment VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(bsc_alignment src/alignment.cpp)

set(GENERATED_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated_include)
set(GENERATED_ALIGNMENT_HEADER ${GENERATED_INCLUDE_DIR}/alignment/alignment.hpp)
file(MAKE_DIRECTORY ${GENERATED_INCLUDE_DIR}/alignment)
file(WRITE ${GENERATED_ALIGNMENT_HEADER} "#pragma once\n#include \"bsc/alignment.hpp\"\nnamespace bsc = ::alignment;\n")

target_include_directories(bsc_alignment PUBLIC
  include
  ${GENERATED_INCLUDE_DIR}
)

add_subdirectory(test)
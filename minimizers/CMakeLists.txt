cmake_minimum_required(VERSION 3.10)
project(bsc_minimizers VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(bsc_minimizers src/minimizers.cpp)

target_include_directories(bsc_minimizers PUBLIC include)
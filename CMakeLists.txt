cmake_minimum_required(VERSION 3.10)
project(bsc_mapper VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add alignment library
add_subdirectory(alignment)

# Add minimizers library
add_subdirectory(minimizers)

# Main executable
add_executable(bsc_mapper src/main.cpp)

# Link libraries
target_link_libraries(bsc_mapper bsc_alignment bsc_minimizers)

# Include directories
target_include_directories(bsc_mapper PRIVATE alignment/include minimizers/include)
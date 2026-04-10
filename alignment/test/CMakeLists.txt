add_executable(alignment_test alignment_test.cpp)

target_link_libraries(alignment_test
  PRIVATE
    bsc_alignment
    bioparser)

target_include_directories(alignment_test
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_test(NAME alignment_test COMMAND alignment_test)

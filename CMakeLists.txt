cmake_minimum_required(VERSION 3.4 FATAL_ERROR)
project(mask-generator CXX)

add_subdirectory(extern/pybind11)
include_directories(extern)

include_directories(${PROJECT_SOURCE_DIR}/maskgen)
add_subdirectory(maskgen)
add_subdirectory(example)

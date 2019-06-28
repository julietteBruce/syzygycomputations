cmake_minimum_required(VERSION 3.10)
project (SyzygyComputation)
#this is a bit broken, but will work for now, in generaly we need to tell the python scripts where the binaries are
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)
add_subdirectory(src/CreateMatrices)

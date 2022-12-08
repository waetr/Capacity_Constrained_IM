cmake_minimum_required(VERSION 3.10)
project(exp)
set(CMAKE_CXX_STANDARD 14)
MATH(EXPR stack_size "1024*1024*1024")
set(CMAKE_EXE_LINKER_FLAGS "-O3,-Wl,--stack,${stack_size}")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(exp src/large.cpp)
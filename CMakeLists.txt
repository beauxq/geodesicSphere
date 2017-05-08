cmake_minimum_required(VERSION 3.7)
project(geodesicSphere)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp gs.h)
add_executable(geodesicSphere ${SOURCE_FILES})
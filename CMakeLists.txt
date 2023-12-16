cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 20)
project(FEM-Study)


set(SOURCES
    src/utils.h
    )

add_executable(FEM ${SOURCES} src/main.cpp)

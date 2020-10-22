cmake_minimum_required(VERSION 3.7.0)
project(Shinsung VERSION 0.1.1)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "-std=c++11")

include_directories(
    "header"
)

add_executable(Shinsung
    main.cpp
    source/function.cpp
    source/matrix.cpp
    source/Initialize.cpp
    source/ReadIMU.cpp
)


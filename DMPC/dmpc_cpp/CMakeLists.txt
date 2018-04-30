cmake_minimum_required(VERSION 3.9)
project(dmpc_cpp)

set(CMAKE_CXX_STANDARD 11)

add_executable(dmpc_cpp
        main.cpp)

# Declare a C++ library
add_library(dmpc
        dmpc.cpp
        )

target_link_libraries(dmpc_cpp
        dmpc
        )

cmake_minimum_required(VERSION 3.9)
project(dmpc_cpp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
# a CMake module named "FindGUROBI.cmake" is available in cmake/modules/
set(GUROBI_LIBRARY  "/home/carlos/opt/gurobi752/linux64/lib/libgurobi75.so")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
find_package(GUROBI REQUIRED)

if(GUROBI_FOUND)
    message(STATUS "GUROBI FOUND")
else(GUROBI_FOUND)
    message(FATAL_ERROR "!!! GUROBI NOT FOUND !!!")
endif(GUROBI_FOUND)

include_directories(${GUROBI_INCLUDE_DIRS})

set(OOQP_LIBS ooqpbase ooqpdense ooqpgensparse ooqpmehrotra ooqpbound
        ooqpgendense ooqpgondzio ooqpsparse blas ma27 gfortran)

add_executable(dmpc_cpp
        main.cpp)

# Declare a C++ library
add_library(dmpc
        dmpc.cpp)

target_link_libraries(dmpc_cpp
        dmpc
        eigen-quadprog
        eigen-gurobi
        pthread
        ooqpei
        ${GUROBI_LIBRARIES}
        ${OOQP_LIBS}
        )


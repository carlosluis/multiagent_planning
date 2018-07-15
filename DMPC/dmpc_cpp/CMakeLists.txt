cmake_minimum_required(VERSION 3.9)
project(dmpc_cpp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

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
        pthread
        ooqpei
        ${OOQP_LIBS}
        )


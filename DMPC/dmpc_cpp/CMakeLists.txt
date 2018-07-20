cmake_minimum_required(VERSION 3.9)
project(dmpc_cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake-modules")

find_package(cplex REQUIRED)

include_directories(${CPLEX_INCLUDE_DIR} ${CPLEX_CONCERT_INCLUDE_DIR})

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
        ${CPLEX_LIBRARIES}
        ${CMAKE_DL_LIBS}
        )

add_executable(test_cplex
        iloqpex1.cpp)

set_property(
        TARGET test_cplex
        PROPERTY COMPILE_DEFINITIONS IL_STD
)
target_link_libraries(test_cplex
        ${CPLEX_LIBRARIES}
        ${CMAKE_DL_LIBS})


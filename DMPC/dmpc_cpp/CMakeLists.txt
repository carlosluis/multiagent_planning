cmake_minimum_required(VERSION 3.9)
project(dmpc_cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake-modules")

find_package(cplex REQUIRED)

find_package(qpoases REQUIRED)

include_directories(${CPLEX_INCLUDE_DIR} ${CPLEX_CONCERT_INCLUDE_DIR} ${QPOASES_INCLUDE_DIR})

MESSAGE("-- Check for LAPACK")
IF (APPLE)
    # Mac OS X has LAPACK/BLAS bundled in a framework called "vecLib". Search for
    # that instead of for the normal "lapack" library.
    FIND_LIBRARY(LAPACK_LIB NAMES vecLib)
ELSE (APPLE)
    FIND_LIBRARY(LAPACK_LIB NAMES lapack)
ENDIF (APPLE)

IF (NOT EXISTS ${LAPACK_LIB})
    MESSAGE(FATAL_ERROR
            "Can't find LAPACK "
            "Please specify -DLAPACK_LIB=...")
ENDIF (NOT EXISTS ${LAPACK_LIB})
MESSAGE ("-- Found LAPACK library: ${LAPACK_LIB}")

# On Apple BLAS is linked with vecLib, so don't add it separately.
IF (NOT APPLE)
    MESSAGE("-- Check for BLAS")
    FIND_LIBRARY(BLAS_LIB NAMES blas)
    IF (NOT EXISTS ${BLAS_LIB})
        MESSAGE(FATAL_ERROR
                "Can't find BLAS (needed for LAPACK and SuiteSparse.)"
                "Please specify -DBLAS_LIB=...")
    ENDIF (NOT EXISTS ${BLAS_LIB})
    MESSAGE ("-- Found BLAS library: ${BLAS_LIB}")
ENDIF (NOT APPLE)


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
        ${LAPACK_LIB}
        ${BLAS_LIB}
        ${QPOASES_LIBRARY}
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

add_executable(test_qpoases
        example1.cpp)


target_link_libraries(test_qpoases
        ${LAPACK_LIB}
        ${BLAS_LIB}
        ${QPOASES_LIBRARY}
        )


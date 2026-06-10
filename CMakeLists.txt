cmake_minimum_required(VERSION 3.15)


project(MolAlignLib VERSION 1.0 LANGUAGES Fortran C)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Fortran compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "7.0")
        message(FATAL_ERROR "GNU Fortran compiler version 7.0 or higher is required.")
    endif()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Wunused")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -stand f08")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback -fpe0 -check bounds -warn all")
endif()

# C compiler flags
if(CMAKE_C_COMPILER_ID MATCHES "GNU")
    set(CMAKE_C_FLAGS_RELEASE "-O3")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g -Wall -Wextra")
endif()

add_subdirectory(fortran)
add_subdirectory(bindings)

if(SKBUILD)
    message(STATUS "Detected scikit-build environment; enabling Python extension build")
    add_subdirectory(python)
endif()

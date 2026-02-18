cmake_minimum_required(VERSION 3.15)

project(MolAlignLib VERSION 1.0 LANGUAGES Fortran C)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Fortran compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
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

# Set module output directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Create the static library
add_subdirectory(modules)

# Build the original Fortran programs
add_subdirectory(programs)

# Installation rules for the Fortran library
install(TARGETS molalignlib
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
    DESTINATION include/molalignlib
    FILES_MATCHING PATTERN "*.mod"
)

# Build Python extension only when called from scikit-build
if(SKBUILD)
    message(STATUS "Detected scikit-build environment; enabling Python extension build")
    add_subdirectory(bindings)
endif()

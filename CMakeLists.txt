cmake_minimum_required(VERSION 3.15)


# Only Fortran is needed for the library, the programs and the C binding
# (whose wrappers are Fortran bind(C) code; bindings/c only installs the
# header). C is enabled further down, only for the Python extension.
project(MolAlignLib VERSION 1.0 LANGUAGES Fortran)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Fortran compiler flags
#
# -fvisibility=hidden is added here globally so that every Fortran object
# (molalignlib's core modules as well as the per-binding static libraries
# built in fortran/CMakeLists.txt) defaults to hidden symbol visibility.
# This only affects the *dynamic* symbol table used when a .so is loaded —
# ordinary compile-time/link-time calls between object files bundled into
# the same final shared object are unaffected. The purpose is to stop the
# private copies of molalignlib that get statically linked into the
# atormsd/conformsd bindings from colliding in the global symbol table if
# both end up loaded into the same process (e.g. both Cython extensions
# imported in one Python interpreter).
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "7.0")
        message(FATAL_ERROR "GNU Fortran compiler version 7.0 or higher is required.")
    endif()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008 -fvisibility=hidden")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Wunused")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -stand f08 -fvisibility=hidden")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback -fpe0 -check bounds -warn all")
endif()

add_subdirectory(fortran)
add_subdirectory(bindings/c)

if(SKBUILD)
    message(STATUS "Detected scikit-build environment; enabling Python extension build")

    # The Cython-generated extension is the only C code in the build.
    # enable_language() must run here, in the top-level directory, because
    # it is the highest directory common to the extension and to the
    # molalignlib library it links against.
    enable_language(C)

    # C compiler flags
    if(CMAKE_C_COMPILER_ID MATCHES "GNU")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fvisibility=hidden")
        set(CMAKE_C_FLAGS_RELEASE "-O3")
        set(CMAKE_C_FLAGS_DEBUG "-O0 -g -Wall -Wextra")
    endif()

    add_subdirectory(bindings/python)
endif()

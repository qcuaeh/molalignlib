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
# of the molalign library defaults to hidden symbol visibility. This only
# affects the *dynamic* symbol table used when a shared object is loaded;
# ordinary link-time calls between object files bundled into the same
# shared object are unaffected. When the static library is linked into the
# Python extension, it keeps the library's internal symbols out of the
# extension's dynamic symbol table, so the only exported symbol is
# PyInit_molalignlib and nothing can clash with other libraries loaded into
# the same process.
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "4.8")
        message(FATAL_ERROR
            "GNU Fortran 4.8 or newer is required "
            "(found ${CMAKE_Fortran_COMPILER_VERSION}).")
    endif()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008 -fvisibility=hidden")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Wunused")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -stand f08 -fvisibility=hidden")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback -fpe0 -check bounds -warn all")
endif()

add_subdirectory(source)
add_subdirectory(bindings/c)

if(SKBUILD)
    message(STATUS "Detected scikit-build environment; enabling Python extension build")

    # The Cython-generated extension is the only C code in the build.
    # enable_language() must run here, in the top-level directory, because
    # it is the highest directory common to the extension and to the
    # molalign library it links against.
    enable_language(C)

    # C compiler flags and version check (the C compiler is only used for
    # the Cython-generated extension, so it is only checked here)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU")
        if(CMAKE_C_COMPILER_VERSION VERSION_LESS "5.3")
            message(FATAL_ERROR
                "GCC 5.3 or newer is required to build the Python extension "
                "(found ${CMAKE_C_COMPILER_VERSION}).")
        endif()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fvisibility=hidden")
        set(CMAKE_C_FLAGS_RELEASE "-O3")
        set(CMAKE_C_FLAGS_DEBUG "-O0 -g -Wall -Wextra")
    endif()

    add_subdirectory(bindings/python)
endif()

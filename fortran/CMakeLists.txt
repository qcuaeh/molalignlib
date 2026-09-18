set(FORTRAN_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# molalignlib is consumed by the C/Cython bindings, so it needs the
# ISO-C-compatible kinds module.
set(KINDS_MODULE kinds_iso_c.f90)

# ---------------------------------------------------------------------------
# Helper: append one module (by filename, relative to FORTRAN_SRC_DIR)
# to a source list variable. A macro runs in the caller's scope, so it
# can update out_var directly with no PARENT_SCOPE plumbing needed.
# ---------------------------------------------------------------------------
macro(add_module_source out_var filename)
    list(APPEND ${out_var} ${FORTRAN_SRC_DIR}/${filename})
endmacro()

# ---------------------------------------------------------------------------
# Collect sources: the kinds module first, then the rest of modules.txt,
# then clibs.txt (one bare filename per line, same as modules.txt — just
# the C-binding sources such as clib_atormsd.f90 and clib_conformsd.f90).
# Blank lines and comments are stripped from both files as they're read.
#
# There is only one library target below, so both C bindings' Fortran
# sources are compiled straight into it rather than into per-binding
# static libraries; clibs.txt just keeps them listed separately from the
# core modules for clarity, without needing a separate library per entry.
# ---------------------------------------------------------------------------
set(LIB_SOURCES "")
add_module_source(LIB_SOURCES ${KINDS_MODULE})

file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/modules.txt _module_lines)
foreach(line ${_module_lines})
    string(REGEX MATCH "^#" _comment "${line}")
    if(NOT _comment)
        string(STRIP "${line}" line)
        if(line)
            add_module_source(LIB_SOURCES ${line})
        endif()
    endif()
endforeach()

file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/clibs.txt _clib_lines)
foreach(line ${_clib_lines})
    string(REGEX MATCH "^#" _comment "${line}")
    if(NOT _comment)
        string(STRIP "${line}" line)
        if(line)
            add_module_source(LIB_SOURCES ${line})
        endif()
    endif()
endforeach()

# ---------------------------------------------------------------------------
# Target: molalignlib
#   Uses kinds_iso_c.f90 (for C bindings). Part of the default build.
#   All module sources live here, including clib_atormsd.f90 and
#   clib_conformsd.f90 — a single library, not one per binding. Programs
#   below and the Cython extensions (see python/CMakeLists.txt) all link
#   against this one target so nothing is compiled twice.
#
#   Symbol visibility defaults to hidden (-fvisibility=hidden, set
#   globally in the top-level CMakeLists.txt). Each Cython extension
#   still statically links this library independently, so each .so ends
#   up with its own private copy of molalignlib's symbols; hidden
#   visibility keeps those private copies out of the dynamic symbol
#   table so they don't collide if both extensions are imported into the
#   same Python process. This is unaffected by using one library instead
#   of two, since visibility only governs the *dynamic* symbol table
#   used when a .so is loaded, not which static library the objects came
#   from.
# ---------------------------------------------------------------------------
add_library(molalignlib STATIC ${LIB_SOURCES})

# OUTPUT_NAME is set explicitly because CMake would otherwise derive the
# archive name from the target name itself (libmolalignlib.a). The CMake
# target stays named "molalignlib" everywhere else in the build (programs,
# the Cython extensions in python/CMakeLists.txt, install() below) — only
# the file that lands on disk is renamed, to libmolalign.a.
set_target_properties(molalignlib PROPERTIES
    OUTPUT_NAME molalign
    POSITION_INDEPENDENT_CODE ON
    Fortran_MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/modules/molalignlib
)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/molalignlib)

target_include_directories(molalignlib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules/molalignlib>
    $<INSTALL_INTERFACE:include/molalignlib>
)

install(TARGETS molalignlib
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

# ---------------------------------------------------------------------------
# Targets: programs (skipped when called from scikit-build)
#   Each executable links against molalignlib directly — no intermediate
#   OBJECT library needed, so no sources are compiled a second time.
# ---------------------------------------------------------------------------
if(NOT SKBUILD)

    file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/programs.txt _prog_lines)

    foreach(line ${_prog_lines})
        string(REGEX MATCH "^#" _comment "${line}")
        if(NOT _comment)
            string(STRIP "${line}" line)
            if(line)
                string(REPLACE " " ";" pair "${line}")
                list(LENGTH pair pair_len)
                if(pair_len EQUAL 2)
                    list(GET pair 0 src_rel)
                    list(GET pair 1 exec_name)

                    add_executable(${exec_name}
                        ${FORTRAN_SRC_DIR}/${src_rel}
                    )

                    target_link_libraries(${exec_name} PRIVATE molalignlib)

                    set_target_properties(${exec_name} PROPERTIES
                        Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/molalignlib
                    )

                    install(TARGETS ${exec_name} RUNTIME DESTINATION bin)

                endif()
            endif()
        endif()
    endforeach()

endif()
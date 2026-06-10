set(FORTRAN_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# ---------------------------------------------------------------------------
# Helper: read a source-list file, strip blank lines and comments,
# and prepend the Fortran directory to each entry.
# ---------------------------------------------------------------------------
function(read_source_list filename out_var)
    file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/${filename} _lines)
    set(_sources "")
    foreach(line ${_lines})
        string(REGEX MATCH "^#" _comment "${line}")
        if(NOT _comment)
            string(STRIP "${line}" line)
            if(line)
                list(APPEND _sources ${FORTRAN_SRC_DIR}/${line})
            endif()
        endif()
    endforeach()
    set(${out_var} "${_sources}" PARENT_SCOPE)
endfunction()

# ---------------------------------------------------------------------------
# Collect sources
# ---------------------------------------------------------------------------
read_source_list(library_modules.txt LIB_SOURCES)

# ---------------------------------------------------------------------------
# Target: molalignlib
#   Uses kinds_iso_c.f90 (for C bindings). Part of the default build.
#   All module sources live here; programs link against this library so
#   nothing is compiled twice.
# ---------------------------------------------------------------------------
add_library(molalignlib STATIC ${LIB_SOURCES})

set_target_properties(molalignlib PROPERTIES
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

    file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/fortran_programs.txt _prog_lines)

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

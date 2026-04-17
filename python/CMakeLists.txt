# ---------------------------------------------------------------------------
# Python extension build (only reached when called from scikit-build)
# ---------------------------------------------------------------------------

find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)

# Query the correct extension suffix (.cpython-312-x86_64-linux-gnu.so etc.)
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c
        "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"
    OUTPUT_VARIABLE Python3_EXT_SUFFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "Python extension suffix: ${Python3_EXT_SUFFIX}")

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# ---------------------------------------------------------------------------
# Helper macro: Cythonize a .pyx file, build a MODULE library, install it.
# ---------------------------------------------------------------------------
macro(add_cython_extension EXT_NAME PYX_FILE)
    set(_pyx ${CMAKE_CURRENT_SOURCE_DIR}/molalignlib/${PYX_FILE})
    set(_c   ${CMAKE_CURRENT_BINARY_DIR}/${EXT_NAME}.c)

    add_custom_command(
        OUTPUT  ${_c}
        COMMAND ${Python3_EXECUTABLE} -m cython
                -3
                -I ${CMAKE_SOURCE_DIR}/bindings
                -o ${_c}
                ${_pyx}
        DEPENDS ${_pyx}
        COMMENT "Cythonizing molalignlib/${PYX_FILE}"
        VERBATIM
    )

    add_library(${EXT_NAME}_ext MODULE ${_c})
    target_link_libraries(${EXT_NAME}_ext PRIVATE molalignlib)

    # Set the output name to include the correct Python ABI suffix so that
    # Python's import system finds the freshly built .so and not any stale
    # installed copy with a plain name.
    set_target_properties(${EXT_NAME}_ext PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        PREFIX        ""
        OUTPUT_NAME   "${EXT_NAME}"
        SUFFIX        "${Python3_EXT_SUFFIX}"
    )

    target_include_directories(${EXT_NAME}_ext PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/molalignlib
        ${CMAKE_SOURCE_DIR}/bindings
        ${Python3_INCLUDE_DIRS}
        ${Python3_NumPy_INCLUDE_DIRS}
    )

    if(Python3_LIBRARIES)
        target_link_libraries(${EXT_NAME}_ext PRIVATE ${Python3_LIBRARIES})
    endif()

    install(TARGETS ${EXT_NAME}_ext
        LIBRARY DESTINATION molalignlib
    )
endmacro()

# ---------------------------------------------------------------------------
# Extensions
# ---------------------------------------------------------------------------
add_cython_extension(atormsd   atormsd.pyx)
add_cython_extension(conformsd conformsd.pyx)

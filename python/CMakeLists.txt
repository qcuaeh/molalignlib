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
#   The extension links directly against the single molalignlib static
#   library built in fortran/CMakeLists.txt (which contains
#   clib_atormsd.f90 and clib_conformsd.f90 alongside the core modules,
#   rather than each binding getting its own static library).
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

    # Belt-and-suspenders: also hide symbols by default in the extension's
    # own compiled C code (the Cython-generated wrapper). PyInit_<name>
    # stays exported regardless, since Python.h marks it with default
    # visibility explicitly via PyMODINIT_FUNC. This keeps the extension
    # from exposing anything beyond its module init function.
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${EXT_NAME}_ext PRIVATE -fvisibility=hidden)
    endif()

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
add_cython_extension(molalign  molalign.pyx)
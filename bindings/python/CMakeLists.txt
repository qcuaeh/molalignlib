# ---------------------------------------------------------------------------
# Python extension build (only reached when called from scikit-build)
#
# The whole Python interface is a single Cython source, molalignlib.pyx,
# compiled into a single top-level extension module (molalignlib.<abi>.so):
# no package directory and no separate .py file.
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
# Helper macro: Cythonize a .pyx file, build a MODULE library, install it
# at the top level of the wheel (scikit-build-core sets the install prefix
# to the wheel's platlib directory, so DESTINATION . is site-packages).
#   EXT_NAME must match the .pyx file name, since Cython derives the
#   module name (and PyInit_<name>) from it. The extension links against
#   the molalign static library built in source/CMakeLists.txt.
# ---------------------------------------------------------------------------
macro(add_cython_extension EXT_NAME PYX_FILE)
    set(_pyx ${CMAKE_CURRENT_SOURCE_DIR}/${PYX_FILE})
    set(_c   ${CMAKE_CURRENT_BINARY_DIR}/${EXT_NAME}.c)

    add_custom_command(
        OUTPUT  ${_c}
        COMMAND ${Python3_EXECUTABLE} -m cython
                -3
                -I ${CMAKE_SOURCE_DIR}/bindings/c
                -o ${_c}
                ${_pyx}
        DEPENDS ${_pyx}
        COMMENT "Cythonizing ${PYX_FILE}"
        VERBATIM
    )

    add_library(${EXT_NAME}_ext MODULE ${_c})
    target_link_libraries(${EXT_NAME}_ext PRIVATE MolAlignLib::molalign)

    # The top-level CMakeLists adds -fvisibility=hidden to CMAKE_C_FLAGS.
    # On Python < 3.9, PyMODINIT_FUNC has no visibility attribute, so that
    # flag hides PyInit_<name> and the module can't be imported. Target
    # options come after CMAKE_C_FLAGS on the command line, so this wins.
    # All other Cython-generated functions are static, so nothing else
    # gets exported.
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${EXT_NAME}_ext PRIVATE -fvisibility=default)
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
        ${CMAKE_SOURCE_DIR}/bindings/c
        ${Python3_INCLUDE_DIRS}
        ${Python3_NumPy_INCLUDE_DIRS}
    )

    if(Python3_LIBRARIES)
        target_link_libraries(${EXT_NAME}_ext PRIVATE ${Python3_LIBRARIES})
    endif()

    install(TARGETS ${EXT_NAME}_ext
        LIBRARY DESTINATION .
    )
endmacro()

# ---------------------------------------------------------------------------
# Extension
# ---------------------------------------------------------------------------
add_cython_extension(molalignlib molalignlib.pyx)

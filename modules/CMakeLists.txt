# Read module files from module_files.txt
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/module_files.txt MODULE_FILES)

# Filter out commented lines
set(FORTRAN_MODULE_SOURCES "")
foreach(line ${MODULE_FILES})
    string(REGEX MATCH "^#" is_comment ${line})
    if(NOT is_comment)
        string(STRIP ${line} line)
        if(line)
            list(APPEND FORTRAN_MODULE_SOURCES ${line})
        endif()
    endif()
endforeach()

# Create the static library
add_library(molalignlib STATIC ${FORTRAN_MODULE_SOURCES})

# Set library properties
set_target_properties(molalignlib PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Make sure module directory is created
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# Add module directory to include path
target_include_directories(molalignlib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<INSTALL_INTERFACE:include/molalignlib>
)

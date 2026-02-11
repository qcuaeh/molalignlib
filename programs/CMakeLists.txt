# Read program files from program_files.txt
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/program_files.txt PROGRAM_LINES)

# Parse each line and create executables
foreach(line ${PROGRAM_LINES})
    string(REGEX MATCH "^#" is_comment ${line})
    if(NOT is_comment)
        string(STRIP ${line} line)
        if(line)
            # Split line into source file and executable name
            string(REPLACE " " ";" line_list ${line})
            list(LENGTH line_list list_length)
            
            if(list_length EQUAL 2)
                list(GET line_list 0 source_file)
                list(GET line_list 1 exec_name)
                
                # Create executable
                add_executable(${exec_name} ${source_file})
                
                # Link with static library
                target_link_libraries(${exec_name} PRIVATE molalignlib)
                
                # Set module directory
                set_target_properties(${exec_name} PROPERTIES
                    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
                )
                
                # Installation
                install(TARGETS ${exec_name} RUNTIME DESTINATION bin)
            endif()
        endif()
    endif()
endforeach()

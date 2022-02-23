if ( NOT WIN32 )
    # Symlink the test data
    list( APPEND test_data
        gsw_data_v3_0.nc
    )

    foreach( FILENAME ${test_data} )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                         ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
                         ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
    endforeach()
endif()

# Get a list of test sources
list( APPEND test_srcs
    gsw_poly_check.f90
    gsw_check_functions.f90
)

# If both static and shared libraries are built,
# use the 1st one to link the tests
list( GET LIB_TARGETS 0 _LIB_TARGET )

# Create executable for each test program, and add test
foreach( src ${test_srcs} )

    if ( WIN32 )
       string( REPLACE ".f90" "" exec_name ${src} )
    else()
       string( REPLACE ".f90" ".x" exec_name ${src} )
    endif()
    add_executable( ${exec_name} ${src} )
    add_dependencies( ${exec_name} ${LIB_TARGETS} )
    target_link_libraries( ${exec_name} PUBLIC ${_LIB_TARGET} )
    set_target_properties( ${exec_name} PROPERTIES BUILD_WITH_INSTALL_RPATH FALSE )
    set_target_properties( ${exec_name} PROPERTIES SKIP_BUILD_RPATH         FALSE )
    target_include_directories( ${exec_name}
        PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${INCLUDE_DIR}> )
    set_property(TARGET ${exec_name} PROPERTY FOLDER gsw)


    string( REPLACE ".f90" ""   test_name ${src} )
    add_test( NAME ${test_name} COMMAND ${exec_name} )

endforeach()

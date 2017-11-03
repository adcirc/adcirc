
IF(BUILD_HPXADCIRC)

    add_executable(
        hpxadcirc
        src/hpxadcirc.cpp
    )
    target_link_libraries(hpxadcirc hpx hpx_init)
    hpx_setup_target(
        COMPONENT_DEPENDENCIES iostreams
    )
ENDIF(BUILD_HPXADCIRC)
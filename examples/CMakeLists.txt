set( example_programs
    allencahn
    cahnhilliard
    degcahnhilliard
    loretimarch
    )

include_directories(
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/vendor/sundials-2.4.0/include
    ${PROJECT_BINARY_DIR}/vendor/sundials-2.4.0/include
    )

foreach( program ${example_programs} )
    add_executable( ${program} ${program}.cpp )
    target_link_libraries( ${program}
        odeity
        sundials_cvode_static
        sundials_nvecserial_static
        ${extraLibs}
        )
endforeach()

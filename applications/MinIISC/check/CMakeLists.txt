include(CTest)

#
# instances in the Coloring data set and their optimal objective value
#
set(instances
    "prob.10.30.100.0\;2"
    "prob.15.40.100.1\;3"
    "prob.20.50.100.0\;2"
    "prob.5.030.100.0\;3"
)

#
# add a test that builds the executable
#
add_test(NAME applications-miniisc-build
        COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config $<CONFIG> --target miniisc
        )
#
# avoid that several build jobs try to concurrently build the SCIP library
# note that this ressource lock name is not the actual libscip target
#
set_tests_properties(applications-miniisc-build
                    PROPERTIES
                        RESOURCE_LOCK libscip
                    )

#
# add a test for every instance
#
foreach(instance ${instances})
    list(GET instance 0 basename)
    list(GET instance 1 optval)
    math(EXPR optvalminusone "${optval} - 1")

    #
    # call the Binpacking binary and validate the solve with the given objective value
    #
    add_test(NAME "applications-miniisc-${basename}"
            COMMAND $<TARGET_FILE:miniisc> ${CMAKE_CURRENT_SOURCE_DIR}/../data/${basename}.lp
            )
    set_tests_properties("applications-miniisc-${basename}"
                        PROPERTIES
                            PASS_REGULAR_EXPRESSION "Primal Bound       : \\+(${optval}\\.000000*|${optvalminusone}\\.999999*)e\\+00"
                            FAIL_REGULAR_EXPRESSION "ERROR"
                            DEPENDS applications-miniisc-build
                            RESOURCE_LOCK libscip
                        )
    if(WIN32)
        # on windows we need to execute the application and examples executables from the directory containing the libscip.dll,
        # on other systems this directory does not exist.
        set_tests_properties("applications-miniisc-${basename}"
                PROPERTIES
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin/$<CONFIG>
            )
    endif()
endforeach()

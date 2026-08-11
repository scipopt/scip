include(CTest)

#
# varying number of particles and optimal solutions
#
set(instances
    "2\;-1"
    "3\;-3"
    "4\;-6"
)

add_test(NAME examples-lj-build
         COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config $<CONFIG> --target lj
        )
#
# avoid that several build jobs try to concurrently build the SCIP library
# note that this resource lock name is not the actual libscip target
#
set_tests_properties(examples-lj-build
                     PROPERTIES
                        RESOURCE_LOCK libscip
                    )

#
# add a test for every instance
#
foreach(instance ${instances})
    list(GET instance 0 nparticles)
    list(GET instance 1 optval)
    #
    # call the Lennard-Jones binary and validate the solve with the given objective value
    #
    add_test(NAME "examples-lj-${nparticles}"
             COMMAND $<TARGET_FILE:lj> ${nparticles}
            )
    set_tests_properties("examples-lj-${nparticles}"
                         PROPERTIES
                            PASS_REGULAR_EXPRESSION "(optimal solution found|gap limit reached).*Primal Bound *: ${optval}\\.0000.*e\\+00"
                            FAIL_REGULAR_EXPRESSION "ERROR"
                            DEPENDS examples-lj-build
                            RESOURCE_LOCK libscip
                        )
    if(WIN32)
        # on windows we need to execute the application and examples executables from the directory containing the libscip.dll,
        # on other systems this directory does not exist.
        set_tests_properties("examples-lj-${nparticles}"
                PROPERTIES
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin/$<CONFIG>
            )
    endif()
endforeach()

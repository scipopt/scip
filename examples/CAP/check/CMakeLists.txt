include(CTest)

#
# instances of varying size in the CAP data set.
#
set(instances
    "cap81\;9"
    "cap82\;17"
    "cap83\;27" #TODO fix objective values, add other instances
)

add_test(NAME examples-cap-build
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target cap
        )
#
# avoid that several build jobs try to concurrently build the SCIP library
# note that this ressource lock name is not the actual libscip target
#
set_tests_properties(examples-cap-build
                    PROPERTIES
                        RESOURCE_LOCK libscip
                    )

#
# add a test for every instance
#
foreach(instance ${instances})
    list(GET instance 0 basename)
    list(GET instance 1 optval)
    #
    # call the CAP binary and validate the solve with the given objective value
    #
    add_test(NAME "examples-cap-${basename}"
            COMMAND $<TARGET_FILE:cap> -f ${CMAKE_CURRENT_SOURCE_DIR}/../data/${basename}.cap -o ${optval} ${optval}
            )
    set_tests_properties("examples-cap-${basename}"
                        PROPERTIES
                            PASS_REGULAR_EXPRESSION "Validation         : Success"
                            DEPENDS examples-cap-build
                        )
endforeach()
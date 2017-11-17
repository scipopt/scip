include(CTest)

#
# instances of varying size in the Binpacking data set.
#
set(instances
    "u20_00\;9"
    "u40_00\;17"
    "u60_00\;27"
#    "u80_00\;33"
#    "u100_00\;40"
#    "u120_00\;48"
# the remaining 3 from this set take too long.
)

add_test(NAME examples-binpacking-build
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target binpacking
        )
#
# avoid that several build jobs try to concurrently build the SCIP library
# note that this ressource lock name is not the actual libscip target
#
set_tests_properties(examples-binpacking-build
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
    # call the Binpacking binary and validate the solve with the given objective value
    #
    add_test(NAME "examples-binpacking-${basename}"
            COMMAND $<TARGET_FILE:binpacking> -f ${CMAKE_CURRENT_SOURCE_DIR}/../data/${basename}.bpa -o ${optval} ${optval}
            )
    set_tests_properties("examples-binpacking-${basename}"
                        PROPERTIES
                            PASS_REGULAR_EXPRESSION "Validation         : Success"
                            DEPENDS examples-binpacking-build
                        )
endforeach()
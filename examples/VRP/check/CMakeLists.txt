include(CTest)
#
# define the instance sets
#
# semicolon '\;' is used to split an instance and its optimal objective value
# For infeasible instances, '+infinity' is used (or '-infinity' in case of maximization)
#
set(instances
    "eil13\;247"
#    "eil22\;375"
#    this instance takes too long
    "eil7\;104"
)

#
# add test to build the VRP executable
#
add_test(NAME examples-vrp-build
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target vrp
        )
#
# avoid that several build jobs try to concurrently build the SCIP library
# note that this ressource lock name is not the actual libscip target
#
set_tests_properties(examples-vrp-build
                    PROPERTIES
                        RESOURCE_LOCK libscip
                    )
#
# loop over the instances
#
foreach(instance ${instances})
    #
    # treat the instance as a tuple (list) of two values
    #
    list(GET instance 0 basename)
    list(GET instance 1 optval)
    add_test(NAME examples-vrp-${basename}
            COMMAND $<TARGET_FILE:vrp> ${CMAKE_CURRENT_SOURCE_DIR}/../data/${basename}.vrp
            )
    set_tests_properties(examples-vrp-${basename}
                        PROPERTIES
                            DEPENDS examples-vrp-build
                            )
endforeach(instance)
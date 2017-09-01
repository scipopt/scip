include(CTest)

#
# instances in the Coloring data set and their optimal objective value
#
set(instances
"1-FullIns_3\;4"
"myciel3\;4"
"queen9_9\;10"
"will199GPIA\;7"
)

#
# add a test that builds the executable
#
add_test(NAME applications-coloring-build
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target coloring
        )
#
# avoid that several build jobs try to concurrently build the SCIP library
# note that this ressource lock name is not the actual libscip target
#
set_tests_properties(applications-coloring-build
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
    add_test(NAME "applications-coloring-${basename}"
            COMMAND $<TARGET_FILE:coloring> -f ${CMAKE_CURRENT_SOURCE_DIR}/../data/${basename}.col -o ${optval} ${optval}
            )
    set_tests_properties("applications-coloring-${basename}"
                        PROPERTIES
                            PASS_REGULAR_EXPRESSION "Validation         : Success"
                            DEPENDS applications-coloring-build
                        )
endforeach()
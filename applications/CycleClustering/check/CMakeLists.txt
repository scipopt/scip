include(CTest)

#
# instances in the Coloring data set and their optimal objective value
#
set(instances
"Paper/Pot3/Pot3_T_20.txt_3cluster\;0.005242"
"Paper/Pot4/Pot4_T_20.txt_4cluster\;0.010748"
"Paper/Pot6/Pot6_T_20.txt_6cluster\;0.004931"
)

#
# add a test that builds the executable
#
add_test(NAME applications-cycleclustering-build
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target cycleclustering
        )
#
# avoid that several build jobs try to concurrently build the SCIP library
# note that this ressource lock name is not the actual libscip target
#
set_tests_properties(applications-cycleclustering-build
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
    add_test(NAME "applications-cycleclustering-${basename}"
            COMMAND $<TARGET_FILE:cycleclustering> -f ${CMAKE_CURRENT_SOURCE_DIR}/../instances/${basename}.spa -o ${optval} ${optval}
            )
    set_tests_properties("applications-cycleclustering-${basename}"
                        PROPERTIES
                            PASS_REGULAR_EXPRESSION "Validation         : Success"
                            DEPENDS applications-cycleclustering-build
                        )
endforeach()
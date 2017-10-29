include(CTest)

#
# instances in the stp data set and their optimal objective value
#
set(instances
# "bipe2u\;54"
"cc3-4u\;23"
"d18\;223"
# "e18\;564"
"hc6u\;39"
# "i640-111\;6167"
# "i640-112\;6304"
)

#
# add a test that builds the application
#
add_test(NAME applications-scipstp-build
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target scipstp
        )
#
# avoid that several build jobs try to concurrently build the SCIP library
# note that this ressource lock name is not the actual libscip target
#
set_tests_properties(applications-scipstp-build
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
    # call the STP binary and validate the solve with the given objective value
    #
    add_test(NAME "applications-scipstp-${basename}"
            COMMAND $<TARGET_FILE:scipstp> -f ${CMAKE_CURRENT_SOURCE_DIR}/../data/short/${basename}.stp -o ${optval} ${optval}
            )
    set_tests_properties("applications-scipstp-${basename}"
                        PROPERTIES
                            PASS_REGULAR_EXPRESSION "Validation         : Success"
                            DEPENDS applications-scipstp-build
                        )
endforeach()
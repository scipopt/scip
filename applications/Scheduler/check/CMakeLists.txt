include(CTest)

#
# instances in the Scheduler data set and their optimal objective value
#
set(instances
"j301_1.sm\;43"
"j301_2.sm\;47"
"j301_3.sm\;47"
"j301_4.sm\;62"
"j301_5.sm\;39"
"j301_6.sm\;48"
"j301_7.sm\;60"
"j301_8.sm\;53"
"j301_9.sm\;49"
"j301_10.sm\;45"
"c10j4m3.cmin\;190"
"c10j4m1.cmin\;193"
"c10j2m4.cmin\;167"
"c10j2m5.cmin\;207"
"c10j2m1.cmin\;183"
"c10j2m2.cmin\;169"
"c10j2m3.cmin\;144"
"c10j4m4.cmin\;260"
"c10j4m5.cmin\;217"
"c10j4m2.cmin\;244"
"c10j3m1.cmin\;237"
"c10j3m3.cmin\;216"
"c10j3m2.cmin\;220"
"c10j3m5.cmin\;193"
"c10j3m4.cmin\;227"
)

#
# add a test to build the application
#
add_test(NAME applications-scheduler-build
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target scheduler
        )
#
# avoid that several build jobs try to concurrently build the SCIP library
# note that this ressource lock name is not the actual libscip target
#
set_tests_properties(applications-scheduler-build
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
    # call the Scheduler binary and validate the solve with the given objective value
    #
    add_test(NAME "applications-scheduler-${basename}"
            COMMAND $<TARGET_FILE:scheduler> -f ${CMAKE_CURRENT_SOURCE_DIR}/../data/${basename} -o ${optval} ${optval}
            )
    set_tests_properties("applications-scheduler-${basename}"
                        PROPERTIES
                            PASS_REGULAR_EXPRESSION "Validation         : Success"
                            DEPENDS applications-scheduler-build
                        )
endforeach()
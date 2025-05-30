include(CTest)
#
# define the instance sets
#
# semicolon '\;' is used to split an instance and its optimal objective value
# For infeasible instances, '+infinity' is used (or '-infinity' in case of maximization)
#
set(instances
    "att48\;10628"
    "berlin52\;7542"
    "burma14\;3323"
    "eil51\;426"
    "ulysses16\;6859"
    "ulysses22\;7013"
)

#
# add test that builds the executable
#
add_test(NAME examples-tsp-build
        COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config $<CONFIG> --target sciptsp
        )
#
# avoid that several build jobs try to concurrently build the SCIP library
# note that this ressource lock name is not the actual libscip target
#
set_tests_properties(examples-tsp-build
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
    add_test(NAME examples-tsp-${basename}
            COMMAND $<TARGET_FILE:sciptsp> -f ${CMAKE_CURRENT_SOURCE_DIR}/../tspdata/${basename}.tsp -o ${optval} ${optval}
            )
    set_tests_properties(examples-tsp-${basename}
                        PROPERTIES
                            PASS_REGULAR_EXPRESSION "Validation         : Success"
                            FAIL_REGULAR_EXPRESSION "ERROR"
                            DEPENDS examples-tsp-build
                            RESOURCE_LOCK libscip
                        )
    if(WIN32)
        # on windows we need to execute the application and examples executables from the directory containing the libscip.dll,
        # on other systems this directory does not exist.
        set_tests_properties(examples-tsp-${basename}
                PROPERTIES
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin/$<CONFIG>
            )
    endif()
endforeach(instance)

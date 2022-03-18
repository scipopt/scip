#
# define the instance sets
#
# semicolon '\;' is used to split an instance and its optimal objective value
# For infeasible instances, '+infinity' is used (or '-infinity' in case of maximization)
#

set(instances
    "ex1\;3100"
#     "ex2\;3224"
# this example takes 700 seconds
    "t65i11xx\;16389651"
    "t70x11xx\;343471236"
)

#
# add a test to build the lop example
#
add_test(NAME examples-lop-build
        COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config $<CONFIG> --target lop
        )
#
# avoid that several build jobs try to concurrently build the SCIP library
# note that this ressource lock name is not the actual libscip target
#
set_tests_properties(examples-lop-build
                    PROPERTIES
                        RESOURCE_LOCK libscip
                    )

#
#
# loop over the instances
#
foreach(instance ${instances})
    list(GET instance 0 basename)
    list(GET instance 1 optval)
        #
        # treat the instance as a tuple (list) of two values
        #
        add_test(NAME examples-lop-${basename}
                COMMAND $<TARGET_FILE:lop>  -f ${CMAKE_CURRENT_SOURCE_DIR}/../data/${basename}.lop -o ${optval} ${optval}
                )
        set_tests_properties(examples-lop-${basename}
                            PROPERTIES
                                PASS_REGULAR_EXPRESSION "Validation         : Success"
                                DEPENDS examples-lop-build
                                RESOURCE_LOCK libscip
                            )
        if(WIN32)
            # on windows we need to execute the application and examples executables from the directory containing the libscip.dll,
            # on other systems this directory does not exist.
            set_tests_properties(examples-lop-${basename}
                PROPERTIES
                        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin/$<CONFIG>
                )
        endif()
endforeach(instance)

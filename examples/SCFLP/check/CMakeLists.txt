include(CTest)

#
# instances of varying size in the CAP data set.
#
set(instances
    "cap41\;+1.22951678000000e+06"
    "cap42\;+1.29781129000000e+06"
    "cap43\;+1.36281129000000e+06"
    "cap44\;+1.46031129000000e+06"
)

set(quadinstances
    "cap41\;+1.78141e+05"
    "cap42\;+2.43141e+05"
    "cap43\;+3.08141e+05"
    "cap44\;+4.05641e+05"
)

add_test(NAME examples-scflp-build
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target scflp
        )
#
# avoid that several build jobs try to concurrently build the SCIP library
# note that this ressource lock name is not the actual libscip target
#
set_tests_properties(examples-scflp-build
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
    # call the SCFLP binary and validate the solve with the given objective value
    #
    add_test(NAME "examples-scflp-${basename}"
            COMMAND $<TARGET_FILE:scflp> -s ${CMAKE_CURRENT_SOURCE_DIR}/../settings/nobenders.set -f ${CMAKE_CURRENT_SOURCE_DIR}/../data/${basename}.cap -o ${optval} ${optval}
            )
    set_tests_properties("examples-scflp-${basename}"
                        PROPERTIES
                            PASS_REGULAR_EXPRESSION "Validation         : Success"
                            DEPENDS examples-scflp-build
                        )

    #
    # call the SCFLP binary and validate the solve with the given objective value
    #
    add_test(NAME "examples-scflpbenders-${basename}"
            COMMAND $<TARGET_FILE:scflp> -s ${CMAKE_CURRENT_SOURCE_DIR}/../settings/benders.set -f ${CMAKE_CURRENT_SOURCE_DIR}/../data/${basename}.cap -o ${optval} ${optval}
            )
    set_tests_properties("examples-scflpbenders-${basename}"
                        PROPERTIES
                            PASS_REGULAR_EXPRESSION "Validation         : Success"
                            DEPENDS examples-scflp-build
                        )
endforeach()

#
# add a test for every quadratic instance
#
foreach(instance ${quadinstances})
    list(GET instance 0 basename)
    list(GET instance 1 optval)
    #
    # call the SCFLP binary and validate the solve with the given objective value
    #
    add_test(NAME "examples-scflp-quad-${basename}"
            COMMAND $<TARGET_FILE:scflp> -s ${CMAKE_CURRENT_SOURCE_DIR}/../settings/nobenders-quad.set -f ${CMAKE_CURRENT_SOURCE_DIR}/../data/${basename}.cap -o ${optval} ${optval}
            )
    set_tests_properties("examples-scflp-${basename}"
                        PROPERTIES
                            PASS_REGULAR_EXPRESSION "Validation         : Success"
                            DEPENDS examples-scflp-build
                        )

    #
    # call the SCFLP binary and validate the solve with the given objective value
    #
    add_test(NAME "examples-scflpbenders-quad-${basename}"
            COMMAND $<TARGET_FILE:scflp> -s ${CMAKE_CURRENT_SOURCE_DIR}/../settings/benders-quad.set -f ${CMAKE_CURRENT_SOURCE_DIR}/../data/${basename}.cap -o ${optval} ${optval}
            )
    set_tests_properties("examples-scflpbenders-${basename}"
                        PROPERTIES
                            PASS_REGULAR_EXPRESSION "Validation         : Success"
                            DEPENDS examples-scflp-build
                        )
endforeach()

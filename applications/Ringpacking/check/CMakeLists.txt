include(CTest)

#
# simple test instances for Ringpacking
#
set(instances
#  "circle6\;1"
  "ring1\;1"
  "ring2\;1"
  "ring3\;1"
  )

add_test(NAME applications-ringpacking-build
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ringpacking
  )

#
# avoid that several build jobs try to concurrently build the SCIP library
# note that this ressource lock name is not the actual libscip target
#
set_tests_properties(applications-ringpacking-build
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
  # call the Ringpacking binary and validate the solve with the given objective value
  #
  add_test(NAME "applications-ringpacking-${basename}"
    COMMAND $<TARGET_FILE:ringpacking> -f ${CMAKE_CURRENT_SOURCE_DIR}/../data/${basename}.rpa -o ${optval} ${optval}
    )
  set_tests_properties("applications-ringpacking-${basename}"
    PROPERTIES
    PASS_REGULAR_EXPRESSION "Validation         : Success"
    DEPENDS applications-ringpacking-build
    )
endforeach()

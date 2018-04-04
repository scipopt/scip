include(CTest)

#
# simple test instances for Ringpacking
#
set(instances
  #"test_circle0\;2"
  "test_circle1\;2"
  "test_circle2\;1"
  "test_circle3\;1"
  "test_circle4\;2"
  "test_circle5\;1"
  "test_circle6\;1"
  #"test_circle16\;1"
  "test_ring1\;1"
  )

add_test(NAME examples-ringpacking-build
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ringpacking
  )

#
# avoid that several build jobs try to concurrently build the SCIP library
# note that this ressource lock name is not the actual libscip target
#
set_tests_properties(examples-ringpacking-build
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
  add_test(NAME "examples-ringpacking-${basename}"
    COMMAND $<TARGET_FILE:ringpacking> -f ${CMAKE_CURRENT_SOURCE_DIR}/../data/test/${basename}.rpa -o ${optval} ${optval}
    )
  set_tests_properties("examples-ringpacking-${basename}"
    PROPERTIES
    PASS_REGULAR_EXPRESSION "Validation         : Success"
    DEPENDS examples-ringpacking-build
    )
endforeach()

include(CTest)

#
# simple test instances for Ringpacking
#
set(instances
  "circle6\;1\;2" # the optimal value is 1, but without Ipopt it may not be found switch to a heuristic price-and-branch, so we test only for 2
  "ring1\;1\;1"
  "ring2\;1\;1"
  "ring3\;1\;1"
  )

add_test(NAME applications-ringpacking-build
  COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config $<CONFIG> --target ringpacking
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
  list(GET instance 1 dualval)
  list(GET instance 2 primval)

  #
  # call the Ringpacking binary and validate the solve with the given objective value
  #
  add_test(NAME "applications-ringpacking-${basename}"
    COMMAND $<TARGET_FILE:ringpacking> -f ${CMAKE_CURRENT_SOURCE_DIR}/../data/${basename}.rpa -o ${primval} ${dualval}
    )
  set_tests_properties("applications-ringpacking-${basename}"
    PROPERTIES
    PASS_REGULAR_EXPRESSION "Validation         : Success"
    DEPENDS applications-ringpacking-build
    RESOURCE_LOCK libscip
    )
  if(WIN32)
    # on windows we need to execute the application and examples executables from the directory containing the libscip.dll,
    # on other systems this directory does not exist.
    set_tests_properties("applications-ringpacking-${basename}"
                PROPERTIES
                      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin/$<CONFIG>
      )
  endif()
endforeach()

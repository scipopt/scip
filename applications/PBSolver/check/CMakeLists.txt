include(CTest)

#
# simple test instances for PBSolver
#
set(instances
  "normalized-bsg_10_4_5.opb\;-4"
  "normalized-mds_10_4_3.opb\;2"
  "normalized-mds_50_25_5.opb\;3"
  "normalized-90_rounds_0_errors.opb\;0"
  "normalized-t2001.13queen13.1110976730--soft-33-100-0.wbo\;4"
  )

add_test(NAME applications-pbsolver-build
  COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config $<CONFIG> --target pbsolver
  )

#
# avoid that several build jobs try to concurrently build the SCIP library
# note that this ressource lock name is not the actual libscip target
#
set_tests_properties(applications-pbsolver-build
  PROPERTIES
  RESOURCE_LOCK libscip
  )

#
# add a test for every instance
#
foreach(instance ${instances})
  list(GET instance 0 name)
  list(GET instance 1 optval)
  get_filename_component(basename ${name} NAME_WE)
  #
  # call the PBSolver binary and validate the solve with the given objective value
  #
  add_test(NAME "applications-pbsolver-${basename}"
    COMMAND $<TARGET_FILE:pbsolver> -f ${CMAKE_CURRENT_SOURCE_DIR}/../../../check/instances/PseudoBoolean/${name}
    )
  set_tests_properties("applications-pbsolver-${basename}"
    PROPERTIES
    PASS_REGULAR_EXPRESSION "o ${optval}"
    PASS_REGULAR_EXPRESSION "s OPTIMUM FOUND"
    DEPENDS applications-pbsolver-build
    )
  if(WIN32)
    # on windows we need to execute the application and examples executables from the directory containing the libscip.dll,
    # on other systems this directory does not exist.
    set_tests_properties("applications-pbsolver-${basename}"
                PROPERTIES
                      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin/$<CONFIG>
      )
  endif()
endforeach()

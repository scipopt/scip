cmake_minimum_required(VERSION 3.9)

project(Eventhdlr)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(eventhdlr
   src/cmain.c
   src/event_bestsol.c
   src/event_boundwriting.c)

target_link_libraries(eventhdlr ${SCIP_LIBRARIES})

if( TARGET examples )
    add_dependencies( examples eventhdlr )
endif()

#
# add some trivial tests
#
include(CTest)

#
# define the instance sets
#
# semicolon '\;' is used to split an instance and its optimal objective value
# For infeasible instances, '+infinity' is used (or '-infinity' in case of maximization)
#

#
# CP instances
#
set(instances
    "instances/CP/linking.cip\;2"
    "instances/CP/j301_2.cip\;47"
    "instances/MIP/flugpl.mps\;1201500"
    "instances/MIP/gt2.mps\;21166"
    "instances/MINLP/ex1266.mps\;16.3"
    "instances/MINLP/m3.osil\;37.8"
    "instances/Symmetry/partorb_1-FullIns_3.cip\;4"
    )

#
#
#
add_test(NAME examples-eventhdlr-build
        COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config $<CONFIG> --target eventhdlr
        )
#
# avoid that several build jobs try to concurrently build the SCIP library
# note that this ressource lock name is not the actual libscip target
#
set_tests_properties(examples-eventhdlr-build
                    PROPERTIES
                        RESOURCE_LOCK libscip
                    )
#
# loop over the instances
#
foreach(instance ${instances})
    list(GET instance 0 path)
    list(GET instance 1 optval)
    get_filename_component(basename ${path} NAME)
    #
    # add a test that depends on the build job and executes the example binary for this instance
    #
    add_test(NAME examples-eventhdlr-${basename}
            COMMAND $<TARGET_FILE:eventhdlr> -f ${CMAKE_CURRENT_SOURCE_DIR}/../../check/${path} -o ${optval} ${optval}
            )
    set_tests_properties(examples-eventhdlr-${basename}
                        PROPERTIES
                            PASS_REGULAR_EXPRESSION "Validation         : Success"
                            DEPENDS examples-eventhdlr-build
                            )
    if(WIN32)
        # on windows we need to execute the application and examples executables from the directory containing the libscip.dll,
        # on other systems this directory does not exist.
        set_tests_properties(examples-eventhdlr-${basename}
                PROPERTIES
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin/$<CONFIG>
            )
    endif()
endforeach(instance)

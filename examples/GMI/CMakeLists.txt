cmake_minimum_required(VERSION 3.0)

project(GMI)

find_package(SCIP REQUIRED)
include_directories(${SCIP_INCLUDE_DIRS})

add_executable(gmi
   src/cmain.c
   src/sepa_gmi.c)

# link to math library if it is available
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()

target_link_libraries(gmi ${SCIP_LIBRARIES} ${LIBM})

if( TARGET examples )
    add_dependencies( examples gmi )
endif()

include(CTest)

#
# define the instance sets
#
# semicolon '\;' is used to split an instance and its optimal objective value
# For infeasible instances, '+infinity' is used (or '-infinity' in case of maximization)
#

set(instances
    "instances/MIP/dcmulti.mps\;188182"
    "instances/MIP/misc03.mps\;3360"
    "instances/MIP/stein27.fzn\;18"
    "instances/MINLP/parincQuadratic.osil\;49920.5564"
    "instances/MINLP/tltr.mps\;48.0666666667"
    "instances/PseudoBoolean/factor-mod-size=9-P0=67-P1=349-P2=67-P3=499-P4=79-P5=347-P6=307-B.opb\;3"
    "instances/PseudoBoolean/normalized-bsg_10_4_5.opb\;-4"
    "instances/PseudoBoolean/normalized-mds_10_4_3.opb\;2"
    )
set(settings
    "gmionly"
    "gmirep"
    "gmi"
    "scipdefault")

#
# add the build of the binary as a test
#
add_test(NAME examples-gmi-build
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target gmi
        )
#
# avoid that several build jobs try to concurrently build the SCIP library
# note that this ressource lock name is not the actual libscip target
#
set_tests_properties(examples-gmi-build
                    PROPERTIES
                        RESOURCE_LOCK libscip
                    )
#
#
# loop over the instances
#
foreach(instance ${instances})
    list(GET instance 0 path)
    list(GET instance 1 optval)
    get_filename_component(basename ${path} NAME)
    foreach(setting ${settings})
        #
        # treat the instance as a tuple (list) of two values
        #
        add_test(NAME examples-gmi-${setting}-${basename}
                COMMAND $<TARGET_FILE:gmi> ${CMAKE_CURRENT_SOURCE_DIR}/../../check/${path} ${CMAKE_CURRENT_SOURCE_DIR}/settings/${setting}.set -o ${optval} ${optval}
                )
        set_tests_properties(examples-gmi-${setting}-${basename}
                            PROPERTIES
                                DEPENDS examples-gmi-build
                            )
    endforeach(setting)
endforeach(instance)
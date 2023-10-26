cmake_minimum_required(VERSION 3.9)

project(Queens)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(queens
   src/queens.cpp
   src/queens_main.cpp)

target_link_libraries(queens ${SCIP_LIBRARIES})

if( TARGET examples )
    add_dependencies( examples queens )
endif()

include(CTest)
#
# add a test to build the executable
#
add_test(NAME examples-queens-build
        COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config $<CONFIG> --target queens
        )
#
# avoid that several build jobs try to concurrently build the SCIP library
# note that this ressource lock name is not the actual libscip target
#
set_tests_properties(examples-queens-build
                    PROPERTIES
                        RESOURCE_LOCK libscip
                    )

#
# test the queens example on checker boards of different size
#

set(nSet
    1
    2
    4
    8
    16
    )
foreach(n ${nSet})
    add_test(NAME examples-queens-${n}
            COMMAND $<TARGET_FILE:queens> ${n}
            )
    set_tests_properties(examples-queens-${n}
                        PROPERTIES
                            DEPENDS examples-queens-build
                        )
    if(WIN32)
        # on windows we need to execute the application and examples executables from the directory containing the libscip.dll,
        # on other systems this directory does not exist.
        set_tests_properties(examples-queens-${n}
                PROPERTIES
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin/$<CONFIG>
            )
    endif()
endforeach(n)

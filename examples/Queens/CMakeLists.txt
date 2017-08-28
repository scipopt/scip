project(Queens)

find_package(SCIP REQUIRED)
include_directories(${SCIP_INCLUDE_DIRS})

add_executable(queens
   src/queens.cpp
   src/queens_main.cpp)

target_link_libraries(queens ${SCIP_LIBRARIES})
add_dependencies(examples queens)

include(CTest)
#
# add a test to build the executable
#
add_test(NAME examples-queens-build
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target queens
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
endforeach(n)

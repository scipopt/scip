cmake_minimum_required(VERSION 3.0)

project(CallableLibrary)
include(CTest)

find_package(SCIP REQUIRED)
include_directories(${SCIP_INCLUDE_DIRS})

add_executable(circle
   src/circle.c)

add_executable(gastrans
   src/gastrans.c)

add_executable(string
   src/string.c)

add_executable(brachistochrone
   src/brachistochrone.c)

add_executable(circlepacking
   src/circlepacking.c)

# link to math library if it is available
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()

target_link_libraries(circle ${SCIP_LIBRARIES})
target_link_libraries(gastrans ${SCIP_LIBRARIES} ${LIBM})
target_link_libraries(string ${SCIP_LIBRARIES})
target_link_libraries(brachistochrone ${SCIP_LIBRARIES} ${LIBM})
target_link_libraries(circlepacking ${SCIP_LIBRARIES} ${LIBM})

if( TARGET examples )
    add_dependencies( examples circle gastrans string brachistochrone circlepacking )
endif()

#
# add one test that builds the executable and one that runs it
#
macro(addbuildandtest target)
    add_test(NAME examples-callablelibrary_${target}-build
            COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${target}
            )
    #
    # avoid that several build jobs try to concurrently build the SCIP library
    # note that this ressource lock name is not the actual libscip target
    #
    set_tests_properties(examples-callablelibrary_${target}-build
                        PROPERTIES
                            RESOURCE_LOCK libscip
                        )

    add_test(NAME examples-callablelibrary_${target}
            COMMAND $<TARGET_FILE:${target}>
            )
    set_tests_properties(examples-callablelibrary_${target}
                        PROPERTIES
                            DEPENDS examples-callablelibrary_${target}-build
                        )
endmacro(addbuildandtest)

addbuildandtest(circle)
addbuildandtest(gastrans)
addbuildandtest(string)
addbuildandtest(brachistochrone)
addbuildandtest(circlepacking)

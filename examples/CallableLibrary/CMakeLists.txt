cmake_minimum_required(VERSION 3.11)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

project(CallableLibrary)
include(CTest)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(BEFORE PUBLIC ${SCIP_INCLUDE_DIRS})

add_executable(gastrans
   src/gastrans.c)

add_executable(spring
   src/spring.c)

add_executable(brachistochrone
   src/brachistochrone.c)

add_executable(circlepacking
   src/circlepacking.c)

# link to math library if it is available
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()

target_link_libraries(gastrans ${SCIP_LIBRARIES} ${LIBM})
target_compile_options(gastrans PRIVATE ${SCIP_COMPILE_FLAGS})
target_link_libraries(spring ${SCIP_LIBRARIES})
target_compile_options(spring PRIVATE ${SCIP_COMPILE_FLAGS})
target_link_libraries(brachistochrone ${SCIP_LIBRARIES} ${LIBM})
target_compile_options(brachistochrone PRIVATE ${SCIP_COMPILE_FLAGS})
target_link_libraries(circlepacking ${SCIP_LIBRARIES} ${LIBM})
target_compile_options(circlepacking PRIVATE ${SCIP_COMPILE_FLAGS})

if( TARGET examples )
    add_dependencies( examples gastrans spring brachistochrone circlepacking )
endif()

#
# add one test that builds the executable and one that runs it
#
macro(addbuildandtest target)
    add_test(NAME examples-callablelibrary_${target}-build
            COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config $<CONFIG> --target ${target}
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
			    FAIL_REGULAR_EXPRESSION "ERROR"
                            DEPENDS examples-callablelibrary_${target}-build
                        )
    if(WIN32)
        # on windows we need to execute the application and examples executables from the directory containing the libscip.dll,
        # on other systems this directory does not exist.
        set_tests_properties(examples-callablelibrary_${target}
                PROPERTIES
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin/$<CONFIG>
            )
    endif()
endmacro(addbuildandtest)

addbuildandtest(gastrans)
addbuildandtest(spring)
addbuildandtest(brachistochrone)
addbuildandtest(circlepacking)

cmake_minimum_required(VERSION 3.11)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

project(CallableLibrary)

include(CTest)

# Find SCIP unless we already have it because we build it.
if(NOT TARGET SCIP::SCIP)
  find_package(SCIP
     REQUIRED
     CONFIG
     HINTS "${SCIP_DIR}"
     )
endif()

# Find math library.
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()


add_executable(gastrans
   src/gastrans.c)
target_link_libraries(gastrans SCIP::SCIP ${LIBM})

add_executable(spring
   src/spring.c)
target_link_libraries(spring SCIP::SCIP ${LIBM})

add_executable(brachistochrone
   src/brachistochrone.c)
target_link_libraries(brachistochrone SCIP::SCIP ${LIBM})

add_executable(circlepacking
   src/circlepacking.c)
target_link_libraries(circlepacking SCIP::SCIP ${LIBM})

if(SCIP_HAS_LTO)
   set_target_properties(gastrans PROPERTIES
      INTERPROCEDURAL_OPTIMIZATION 1)
   set_target_properties(spring PROPERTIES
      INTERPROCEDURAL_OPTIMIZATION 1)
   set_target_properties(brachistochrone PROPERTIES
      INTERPROCEDURAL_OPTIMIZATION 1)
   set_target_properties(circlepacking PROPERTIES
      INTERPROCEDURAL_OPTIMIZATION 1)
endif()

if(TARGET examples)
   add_dependencies(examples gastrans spring brachistochrone circlepacking)
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
   set_tests_properties(examples-callablelibrary_${target}-build PROPERTIES
      RESOURCE_LOCK libscip)

   add_test(NAME examples-callablelibrary_${target}
      COMMAND $<TARGET_FILE:${target}>)
   set_tests_properties(examples-callablelibrary_${target} PROPERTIES
			FAIL_REGULAR_EXPRESSION "ERROR"
      DEPENDS examples-callablelibrary_${target}-build
      )
   if(WIN32)
      # on windows we need to execute the application and examples executables from the directory containing the libscip.dll,
      # on other systems this directory does not exist.
      set_tests_properties(examples-callablelibrary_${target}
         PROPERTIES
         WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/bin/$<CONFIG>"
         )
   endif()
endmacro(addbuildandtest)

addbuildandtest(gastrans)
addbuildandtest(spring)
addbuildandtest(brachistochrone)
addbuildandtest(circlepacking)


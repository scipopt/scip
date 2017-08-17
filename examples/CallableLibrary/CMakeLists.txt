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

# link to math library if it is available
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()

target_link_libraries(circle ${SCIP_LIBRARIES})
target_link_libraries(gastrans ${SCIP_LIBRARIES} ${LIBM})
target_link_libraries(string ${SCIP_LIBRARIES})

add_dependencies(examples circle gastrans string)

#
# add tests for the three executables
#
add_test(NAME examples-callablelibrary-circle COMMAND $<TARGET_FILE:circle>)
add_test(NAME examples-callablelibrary-gastrans COMMAND $<TARGET_FILE:gastrans>)
add_test(NAME examples-callablelibrary-string COMMAND $<TARGET_FILE:string>)

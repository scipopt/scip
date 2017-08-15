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

# check if we can use log10() with just including math.h
unset(HAVE_LOG10 CACHE)
check_symbol_exists(log10 "math.h" HAVE_LOG10)
if(NOT HAVE_LOG10)
    # if not we link libm and check again
    unset(HAVE_LOG10 CACHE)
    list(APPEND CMAKE_REQUIRED_LIBRARIES m)
    check_symbol_exists(log10 "math.h" HAVE_LOG10)
    if(HAVE_LOG10)
        set(libm m)
    else()
        message(FATAL_ERROR "No log10() found")
    endif()
endif()

target_link_libraries(circle ${SCIP_LIBRARIES})
target_link_libraries(gastrans ${SCIP_LIBRARIES} ${libm})
target_link_libraries(string ${SCIP_LIBRARIES})

add_dependencies(examples circle gastrans string)

#
# add tests for the three executables
#
add_test(NAME examples-callablelibrary-circle COMMAND $<TARGET_FILE:circle>)
add_test(NAME examples-callablelibrary-gastrans COMMAND $<TARGET_FILE:gastrans>)
add_test(NAME examples-callablelibrary-string COMMAND $<TARGET_FILE:string>)

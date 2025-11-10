cmake_minimum_required(VERSION 3.11)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

project(SCFLP)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(BEFORE PUBLIC ${SCIP_INCLUDE_DIRS})

add_executable(scflp
   src/cmain.c
   src/probdata_scflp.c
   src/reader_scflp.c )
set_property(TARGET scflp PROPERTY INTERPROCEDURAL_OPTIMIZATION ${LTO_AVAILABLE})

# link to math library if it is available
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()

target_link_libraries(scflp ${SCIP_LIBRARIES} ${LIBM})
target_compile_options(scflp PRIVATE ${SCIP_COMPILE_FLAGS})

if( TARGET examples )
    add_dependencies( examples scflp )
endif()

#
# add check subdirectory for tests
#
add_subdirectory(check)

enable_testing()

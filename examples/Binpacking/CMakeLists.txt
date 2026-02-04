cmake_minimum_required(VERSION 3.11)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

project(Binpacking)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG NO_DEFAULT_PATH PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(BEFORE PUBLIC ${SCIP_INCLUDE_DIRS})

add_executable(binpacking
   src/branch_ryanfoster.c
   src/cmain.c
   src/cons_samediff.c
   src/pricer_binpacking.c
   src/probdata_binpacking.c
   src/reader_bpa.c
   src/vardata_binpacking.c)
set_property(TARGET binpacking PROPERTY INTERPROCEDURAL_OPTIMIZATION ${LTO_AVAILABLE})

# link to math library if it is available
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()

target_link_libraries(binpacking ${SCIP_LIBRARIES} ${LIBM})
target_compile_options(binpacking PRIVATE ${SCIP_COMPILE_FLAGS})

if( TARGET examples )
    add_dependencies( examples binpacking )
endif()

#
# add check subdirectory for tests
#
add_subdirectory(check)

enable_testing()

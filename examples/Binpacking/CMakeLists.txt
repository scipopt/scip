cmake_minimum_required(VERSION 3.11)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

project(Binpacking)

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

add_executable(binpacking
   src/branch_ryanfoster.c
   src/cmain.c
   src/cons_samediff.c
   src/pricer_binpacking.c
   src/probdata_binpacking.c
   src/reader_bpa.c
   src/vardata_binpacking.c)
if(SCIP_HAS_LTO)
   set_target_properties(binpacking PROPERTIES
      INTERPROCEDURAL_OPTIMIZATION 1)
endif()
target_link_libraries(binpacking SCIP::SCIP ${LIBM})

if(TARGET examples)
   add_dependencies(examples binpacking)
endif()

# Add check subdirectory for tests.
add_subdirectory(check)

enable_testing()

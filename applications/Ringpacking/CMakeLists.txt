cmake_minimum_required(VERSION 3.11)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

project(Rinpacking)

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

add_executable(ringpacking
   src/cmain.c
   src/reader_rpa.c
   src/probdata_rpa.c
   src/pricer_rpa.c
   src/cons_rpa.c
   src/pattern.c)
target_link_libraries(ringpacking SCIP::SCIP ${LIBM})

if(SCIP_HAS_LTO)
   set_target_properties(ringpacking PROPERTIES
      INTERPROCEDURAL_OPTIMIZATION 1)
endif()

if(TARGET applications)
    add_dependencies(applications ringpacking)
endif()

#
# add check subdirectory for tests
#
add_subdirectory(check)

#
# add tests subdirectory for unit tests
#
add_subdirectory(tests)

enable_testing()

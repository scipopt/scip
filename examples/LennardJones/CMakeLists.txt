cmake_minimum_required(VERSION 3.11)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

project(LennardJones)

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

add_executable(lj
   src/main.c
   src/nlhdlr_lj.c
   )
if(SCIP_HAS_LTO)
   set_target_properties(lj PROPERTIES
      INTERPROCEDURAL_OPTIMIZATION 1)
endif()
target_link_libraries(lj SCIP::SCIP ${LIBM})

if(TARGET examples)
   add_dependencies(examples lj)
endif()

add_subdirectory(check)
enable_testing()

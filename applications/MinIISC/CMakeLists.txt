cmake_minimum_required(VERSION 3.11)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

project(MinIISC)

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


add_executable(miniisc
   src/benders.c
   src/miniisc.c
   src/readargs.c
   )
target_link_libraries(miniisc SCIP::SCIP ${LIBM})

if(SCIP_HAS_LTO)
   set_target_properties(miniisc PROPERTIES
      INTERPROCEDURAL_OPTIMIZATION 1)
endif()

if(TARGET applications)
   add_dependencies(applications miniisc)
endif()

add_subdirectory(check)

enable_testing()

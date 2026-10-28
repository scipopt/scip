cmake_minimum_required(VERSION 3.11)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

project(Coloring)

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

add_executable(coloring
   src/branch_coloring.c
   src/branch_strongcoloring.c
   src/coloringplugins.c
   src/cons_storeGraph.c
   src/heur_init.c
   src/main.c
   src/pricer_coloring.c
   src/probdata_coloring.c
   src/reader_col.c
   src/reader_csol.c
   )
target_link_libraries(coloring SCIP::SCIP ${LIBM})

if(SCIP_HAS_LTO)
   set_target_properties(coloring PROPERTIES
      INTERPROCEDURAL_OPTIMIZATION 1)
endif()

if(TARGET applications)
   add_dependencies(applications coloring)
endif()

add_subdirectory(check)
enable_testing()

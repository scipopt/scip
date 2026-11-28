cmake_minimum_required(VERSION 3.11)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

project(CycleClustering)

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

add_executable(cycleclustering
   src/main.c
   src/cycplugins.c
   src/sepa_partition.c
   src/sepa_edge.c
   src/sepa_subtour.c
   src/heur_cycgreedy.c
   src/heur_fuzzyround.c
   src/heur_cyckerlin.c
   src/heur_redsize.c
   src/branch_multinode.c
   src/event_newsol.c
   src/reader_cyc.c
   src/probdata_cyc.c
   )
target_link_libraries(cycleclustering SCIP::SCIP ${LIBM})

if(SCIP_HAS_LTO)
   set_target_properties(cycleclustering PROPERTIES
      INTERPROCEDURAL_OPTIMIZATION 1)
endif()

if(TARGET applications)
   add_dependencies(applications cycleclustering)
endif()

add_subdirectory(check)
enable_testing()

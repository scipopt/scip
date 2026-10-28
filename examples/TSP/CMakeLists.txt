cmake_minimum_required(VERSION 3.11)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

project(TSP)

# Find SCIP unless we already have it because we build it.
if(NOT TARGET SCIP::SCIP)
  find_package(SCIP
     REQUIRED
     CONFIG
     HINTS "${SCIP_DIR}"
     )
endif()

if(ZLIB)
   find_package(ZLIB)
endif()

add_executable(sciptsp
   src/ConshdlrSubtour.cpp
   src/cppmain.cpp
   src/EventhdlrNewSol.cpp
   src/GomoryHuTree.cpp
   src/Heur2opt.cpp
   src/HeurFarthestInsert.cpp
   src/HeurFrats.cpp
   src/ProbDataTSP.cpp
   src/ReaderTSP.cpp)
target_link_libraries(sciptsp
   SCIP::SCIP
   $<$<BOOL:${ZLIB_FOUND}>:ZLIB::ZLIB>
   )

if(SCIP_HAS_LTO)
   set_target_properties(sciptsp PROPERTIES
      INTERPROCEDURAL_OPTIMIZATION 1)
endif()

if(TARGET examples)
    add_dependencies(examples sciptsp)
endif()

add_subdirectory(check)
enable_testing()

cmake_minimum_required(VERSION 3.11)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

project(TSP)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG NO_DEFAULT_PATH PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(BEFORE PUBLIC ${SCIP_INCLUDE_DIRS})

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
set_property(TARGET sciptsp PROPERTY INTERPROCEDURAL_OPTIMIZATION ${LTO_AVAILABLE})

target_link_libraries(sciptsp ${SCIP_LIBRARIES})
target_compile_options(sciptsp PRIVATE ${SCIP_COMPILE_FLAGS})

if(ZLIB)
    find_package(ZLIB)
    if(ZLIB_FOUND)
        target_link_libraries(sciptsp ${ZLIB_LIBRARIES})
    endif()
endif()

if( TARGET examples )
    add_dependencies( examples sciptsp )
endif()

add_subdirectory(check)
enable_testing()

cmake_minimum_required(VERSION 3.9)

project(TSP)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(sciptsp
   src/ConshdlrSubtour.cpp
   src/cppmain.cpp
   src/EventhdlrNewSol.cpp
   src/GomoryHuTree.cpp
   src/Heur2opt.cpp
   src/HeurFarthestInsert.cpp
   src/HeurFrats.cpp
   src/ProbDataTSP.cpp
   src/ReaderTSP.cpp
   src/gzstream.cpp)

target_link_libraries(sciptsp ${SCIP_LIBRARIES})

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

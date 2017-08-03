project(TSP)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../src)

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

target_link_libraries(sciptsp libscip)
add_dependencies(examples sciptsp)

add_subdirectory(check)
enable_testing()

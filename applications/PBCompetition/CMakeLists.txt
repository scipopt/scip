cmake_minimum_required(VERSION 3.9)

project(PBCompetition)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(pbscip
    src/heur_smallcard.c
    src/heur_indtwoopt.c
    src/heur_indrounding.c
    src/heur_indcoefdiving.c
    src/heur_indoneopt.c
    src/main.cpp
    src/message_pbscip.c
    src/event_bestsol.c
)

target_link_libraries(pbscip ${SCIP_LIBRARIES})

if( TARGET applications )
    add_dependencies( applications pbscip )
endif()

enable_testing()

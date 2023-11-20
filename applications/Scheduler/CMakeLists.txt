cmake_minimum_required(VERSION 3.9)

project(Scheduler)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(scheduler
    src/main.cpp
    src/cons_optcumulative.c
    src/heur_listscheduling.c
    src/heur_optcumulative.c
    src/reader_cmin.c
    src/reader_rcp.c
    src/reader_sch.c
    src/reader_sm.c
)

target_link_libraries(scheduler ${SCIP_LIBRARIES})

if( TARGET applications )
    add_dependencies( applications scheduler )
endif()

add_subdirectory(check)

enable_testing()

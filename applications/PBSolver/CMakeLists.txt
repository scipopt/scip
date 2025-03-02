cmake_minimum_required(VERSION 3.11)

project(PBSolver)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(pbsolver
    src/main.c
    src/message_pb.c
    src/event_bestsol.c
)

target_link_libraries(pbsolver ${SCIP_LIBRARIES})

if( TARGET applications )
    add_dependencies( applications pbsolver )
endif()

enable_testing()

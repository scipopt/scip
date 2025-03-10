cmake_minimum_required(VERSION 3.11)

project(PBSolver)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(BEFORE PUBLIC ${SCIP_INCLUDE_DIRS})

add_executable(pbsolver
    src/main.c
    src/message_pb.c
    src/event_bestsol.c
)

# link to math library if it is available
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()

target_link_libraries(pbsolver ${SCIP_LIBRARIES} ${LIBM})
target_compile_options(pbsolver PRIVATE ${SCIP_COMPILE_FLAGS})

if( TARGET applications )
    add_dependencies( applications pbsolver )
endif()

#
# add check subdirectory for tests
#
add_subdirectory(check)

enable_testing()

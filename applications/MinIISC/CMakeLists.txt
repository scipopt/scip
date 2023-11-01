cmake_minimum_required(VERSION 3.9)

project(MinIISC)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(miniisc
    src/benders.c
    src/miniisc.c
    src/readargs.c
)

# link to math library if it is available
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()

target_link_libraries(miniisc ${SCIP_LIBRARIES} ${LIBM})

if( TARGET applications )
    add_dependencies( applications miniisc )
endif()

add_subdirectory(check)

enable_testing()

cmake_minimum_required(VERSION 3.11)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

project(MinIISC)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG NO_DEFAULT_PATH PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(BEFORE PUBLIC ${SCIP_INCLUDE_DIRS})

add_executable(miniisc
    src/benders.c
    src/miniisc.c
    src/readargs.c
)
set_property(TARGET miniisc PROPERTY INTERPROCEDURAL_OPTIMIZATION ${LTO_AVAILABLE})

# link to math library if it is available
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()

target_link_libraries(miniisc ${SCIP_LIBRARIES} ${LIBM})
target_compile_options(miniisc PRIVATE ${SCIP_COMPILE_FLAGS})

if( TARGET applications )
    add_dependencies( applications miniisc )
endif()

add_subdirectory(check)

enable_testing()

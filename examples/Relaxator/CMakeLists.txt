cmake_minimum_required(VERSION 3.11)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

project(Relaxator)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG NO_DEFAULT_PATH PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(BEFORE PUBLIC ${SCIP_INCLUDE_DIRS})

add_executable(relaxator
   src/cmain.c
   src/relax_lp.c
   src/relax_nlp.c)
set_property(TARGET relaxator PROPERTY INTERPROCEDURAL_OPTIMIZATION ${LTO_AVAILABLE})

target_link_libraries(relaxator ${SCIP_LIBRARIES})
target_compile_options(relaxator PRIVATE ${SCIP_COMPILE_FLAGS})

if( TARGET examples )
    add_dependencies( examples relaxator )
endif()

add_subdirectory(check)
enable_testing()

cmake_minimum_required(VERSION 3.11)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

project(SCFLP)

# Find SCIP unless we already have it because we build it.
if(NOT TARGET SCIP::SCIP)
  find_package(SCIP
     REQUIRED
     CONFIG
     HINTS "${SCIP_DIR}"
     )
endif()

# Find math library.
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()

add_executable(scflp
   src/cmain.c
   src/probdata_scflp.c
   src/reader_scflp.c )
target_link_libraries(scflp SCIP::SCIP ${LIBM})

if(SCIP_HAS_LTO)
   set_target_properties(scflp PROPERTIES
      INTERPROCEDURAL_OPTIMIZATION 1)
endif()

if(TARGET examples)
   add_dependencies(examples scflp)
endif()

# Add check subdirectory for tests.
add_subdirectory(check)

enable_testing()

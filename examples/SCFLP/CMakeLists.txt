cmake_minimum_required(VERSION 3.9)

project(SCFLP)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(scflp
   src/cmain.c
   src/probdata_scflp.c
   src/reader_scflp.c )

# link to math library if it is available
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()

target_link_libraries(scflp ${SCIP_LIBRARIES} ${LIBM})
if( TARGET examples )
    add_dependencies( examples scflp )
endif()

#
# add check subdirectory for tests
#
add_subdirectory(check)

enable_testing()

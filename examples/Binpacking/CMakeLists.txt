cmake_minimum_required(VERSION 3.3)

project(Binpacking)

message("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
message("SCIP_BINARY_DIR: ${SCIP_BINARY_DIR}")


if(TARGET SCIP::SCIP)
  message("Found SCIP")
  # find package by SCIP PATH
  if (EXISTS "${SCIP_BINARY_DIR}/scip-config.cmake")
    message("File found")
  else()
    message("Somethings wrong...")
  endif()
  
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)

else()
  message("Couldn't find scip, searching")
  find_package(SCIP REQUIRED)
endif()

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(binpacking
   src/branch_ryanfoster.c
   src/cmain.c
   src/cons_samediff.c
   src/pricer_binpacking.c
   src/probdata_binpacking.c
   src/reader_bpa.c
   src/vardata_binpacking.c)

# link to math library if it is available
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()

target_link_libraries(binpacking ${SCIP_LIBRARIES} ${LIBM})
if( TARGET examples )
    add_dependencies( examples binpacking )
endif()

#
# add check subdirectory for tests
#
add_subdirectory(check)

enable_testing()

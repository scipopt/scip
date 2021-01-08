cmake_minimum_required(VERSION 3.3)

project(Binpacking)

if(SCIP::SCIP)
  message("SCIP::SCIP is defined")
else()
  message("SCIP::SCIP is NOT defined")
endif()

if(TARGET scip)
  message("scip is defined")
else()
  message("scip is NOT defined")
endif()

if(TARGET SCIP::SCIP)
  message("TARGET SCIP::SCIP is defined")
else()
  message("TARGET SCIP::SCIP is NOT defined")
endif()


if(TARGET SCIP::SCIP)
  message("Found SCIP")
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

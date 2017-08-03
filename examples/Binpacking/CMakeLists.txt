project(Binpacking)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../src)

add_executable(binpacking
   src/branch_ryanfoster.c
   src/cmain.c
   src/cons_samediff.c
   src/pricer_binpacking.c
   src/probdata_binpacking.c
   src/reader_bpa.c
   src/vardata_binpacking.c)

target_link_libraries(binpacking libscip)
add_dependencies(examples binpacking)

#
# add check subdirectory for tests
#
add_subdirectory(check)

enable_testing()

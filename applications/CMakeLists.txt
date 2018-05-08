add_custom_target(applications)

add_subdirectory(Coloring EXCLUDE_FROM_ALL)
add_subdirectory(MinIISC EXCLUDE_FROM_ALL)
add_subdirectory(Ringpacking EXCLUDE_FROM_ALL)
add_subdirectory(Scheduler EXCLUDE_FROM_ALL)
add_subdirectory(STP EXCLUDE_FROM_ALL)
add_subdirectory(CycleClustering EXCLUDE_FROM_ALL)

# disable PolySCIP when compiling with MSVC (TODO enable it again after #1994 has been fixed)
if(NOT MSVC)
  add_subdirectory(PolySCIP EXCLUDE_FROM_ALL)
endif()

enable_testing()

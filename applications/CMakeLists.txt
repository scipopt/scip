add_custom_target(applications)

add_subdirectory(Coloring EXCLUDE_FROM_ALL)
add_subdirectory(MinIISC EXCLUDE_FROM_ALL)
add_subdirectory(PolySCIP EXCLUDE_FROM_ALL)
add_subdirectory(Scheduler EXCLUDE_FROM_ALL)
add_subdirectory(STP EXCLUDE_FROM_ALL)

enable_testing()

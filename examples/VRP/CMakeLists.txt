cmake_minimum_required(VERSION 3.11)

project(VRP)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(BEFORE PUBLIC ${SCIP_INCLUDE_DIRS})

add_executable(vrp
   src/main_vrp.cpp
   src/pricer_vrp.cpp)

target_link_libraries(vrp ${SCIP_LIBRARIES})
target_compile_options(vrp PRIVATE ${SCIP_COMPILE_FLAGS})

if( TARGET examples )
    add_dependencies( examples vrp )
endif()


add_subdirectory(check)
enable_testing()

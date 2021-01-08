cmake_minimum_required(VERSION 3.3)

project(SCIPampl)

if(TARGET SCIP::SCIP)
  message("Found SCIP")
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()

include_directories(${SCIP_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/solvers)

add_executable(scipampl
   src/cmain.c
   src/reader_nl.c)

target_link_libraries(scipampl ${SCIP_LIBRARIES} ${PROJECT_SOURCE_DIR}/solvers/amplsolver.a ${LIBM})

if(CMAKE_DL_LIBS)
    target_link_libraries(scipampl ${CMAKE_DL_LIBS})
endif()

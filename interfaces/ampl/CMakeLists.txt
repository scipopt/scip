cmake_minimum_required(VERSION 3.3)

project(SCIPampl)

find_package(SCIP REQUIRED)
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

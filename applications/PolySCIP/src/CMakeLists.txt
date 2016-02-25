set(HEADERS_FILES
  LiftedWeightSpaceSolver.h 
  ProbDataObjectives.h
  ReaderMOP.h
  Skeleton.h
  WeightedSolver.h
  WeightSpaceVertex.h
  )

set(SOURCE_FILES
  main.cpp 
  LiftedWeightSpaceSolver.cpp 
  ProbDataObjectives.cpp 
  ReaderMOP.cpp 
  Skeleton.cpp 
  WeightedSolver.cpp 
  WeightSpaceVertex.cpp 
)

add_executable(${EXECUTABLE_NAME} ${SOURCE_FILES})
target_link_libraries(${EXECUTABLE_NAME} ${SCIPOPT_LIB} ${SCIP_LIB} ${OBJSCIP_LIB} ${NLPI_LIB} ${LPS_LIB} ${SOPLEX_LIB})

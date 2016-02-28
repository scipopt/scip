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

# check for libsoplex in ${SCIP_ROOT_DIR}/lib directory
file(GLOB SOPLEX_SHARED_LIST ${SCIP_ROOT_DIR}/lib/libsoplex.*so)
list(LENGTH SOPLEX_SHARED_LIST SOPLEX_SHARED_COUNT)
file(GLOB SOPLEX_STATIC_LIST ${SCIP_ROOT_DIR}/lib/libsoplex.*a)
list(LENGTH SOPLEX_STATIC_LIST SOPLEX_STATIC_COUNT)
if(SOPLEX_SHARED_COUNT)
  list(GET SOPLEX_SHARED_LIST 0 HEAD_SHARED)
  get_filename_component(SOPLEX_LIB_SHARED ${HEAD_SHARED} NAME)
  message("USING soplex lib: ${SOPLEX_LIB_SHARED}")
  execute_process(COMMAND readlink ${HEAD_SHARED} OUTPUT_VARIABLE SOPLEX_FULLPATH_REL)
elseif(SOPLEX_STATIC_COUNT) # check for static libsoplex in ${SCIP_ROOT_DIR}/lib directory
  list(GET SOPLEX_STATIC_LIST 0 HEAD_STATIC)
  get_filename_component(SOPLEX_LIB_STATIC ${HEAD_STATIC} NAME)
  message("USING static soplex lib: ${SOPLEX_LIB_STATIC}")
  execute_process(COMMAND readlink ${HEAD_STATIC} OUTPUT_VARIABLE SOPLEX_FULLPATH_REL)
else(SOPLEX_SHARED_COUNT)
  message("No SoPlex library found.")
endif(SOPLEX_SHARED_COUNT)
get_filename_component(SOPLEX_FULLPATH_ABS ${SOPLEX_FULLPATH_REL} ABSOLUTE)
get_filename_component(SOPLEX_LIBDIR ${SOPLEX_FULLPATH_ABS} DIRECTORY)
message("SOPLEX_LIBDIR: ${SOPLEX_LIBDIR}")
link_directories(${SOPLEX_LIBDIR})

# check for libzimpl* in ${SCIP_ROOT_DIR}/lib dic
file(GLOB ZIMPL_LIST ${SCIP_ROOT_DIR}/lib/libzimpl.*)
list(LENGTH ZIMPL_LIST ZIMPL_COUNT)
if(ZIMPL_COUNT)
  list(GET ZIMPL_LIST 0 HEAD_ZIMPL)
  get_filename_component(ZIMPL_LIB ${HEAD_ZIMPL} NAME)
  message("USING zimpl lib: ${ZIMPL_LIB}")
  execute_process(COMMAND readlink ${HEAD_ZIMPL} OUTPUT_VARIABLE ZIMPL_FULLPATH_REL)
  get_filename_component(ZIMPL_FULLPATH_ABS ${ZIMPL_FULLPATH_REL} ABSOLUTE)
  get_filename_component(ZIMPL_LIBDIR ${ZIMPL_FULLPATH_ABS} DIRECTORY)
  message("ZIMPL_LIBDIR: ${ZIMPL_LIBDIR}")
  link_directories(${ZIMPL_LIBDIR})
endif(ZIMPL_COUNT)

add_executable(${EXECUTABLE_NAME} ${SOURCE_FILES})

# add target link libraries depending on whether libzimpl* and which libsoplex* was found
if(ZIMPL_COUNT) 
  if(SOPLEX_SHARED_COUNT)
    target_link_libraries(${EXECUTABLE_NAME} ${SCIP_LIB} ${OBJSCIP_LIB} ${NLPI_LIB} ${LPS_LIB} ${GMP_LIBRARIES} ${Readline_LIBRARY} ${ZLIB_LIBRARIES} ${ZIMPL_LIB} ${SOPLEX_LIB_SHARED})
  elseif(SOPLEX_STATIC_COUNT)
    target_link_libraries(${EXECUTABLE_NAME} ${SCIP_LIB} ${OBJSCIP_LIB} ${NLPI_LIB} ${LPS_LIB} ${GMP_LIBRARIES} ${Readline_LIBRARY} ${ZLIB_LIBRARIES} ${ZIMPL_LIB} ${SOPLEX_LIB_STATIC})
  endif(SOPLEX_SHARED_COUNT)
else(ZIMPL_COUNT) 
  if(SOPLEX_SHARED_COUNT)
    target_link_libraries(${EXECUTABLE_NAME} ${SCIP_LIB} ${OBJSCIP_LIB} ${NLPI_LIB} ${LPS_LIB} ${GMP_LIBRARIES} ${Readline_LIBRARY} ${ZLIB_LIBRARIES} ${SOPLEX_LIB_SHARED})
  elseif(SOPLEX_STATIC_COUNT)
    target_link_libraries(${EXECUTABLE_NAME} ${SCIP_LIB} ${OBJSCIP_LIB} ${NLPI_LIB} ${LPS_LIB} ${GMP_LIBRARIES} ${Readline_LIBRARY} ${ZLIB_LIBRARIES} ${SOPLEX_LIB_STATIC})
  endif(SOPLEX_SHARED_COUNT)
endif(ZIMPL_COUNT)
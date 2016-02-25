SET(LEMON_ROOT_DIR "" CACHE PATH "Lemon root directory")

find_path(LEMON_INC lemon/list_graph.h
  HINTS /usr/include/lemon/include /usr/local/include/lemon/include ${LEMON_ROOT_DIR}/include
)

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(lemon REQUIRED LEMON_INC)
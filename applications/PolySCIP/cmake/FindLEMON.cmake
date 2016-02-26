find_path(LEMON_INC lemon/list_graph.h
  PATHS /usr/local/include/lemon/include /usr/include/lemon/include
  )

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(lemon REQUIRED LEMON_INC)
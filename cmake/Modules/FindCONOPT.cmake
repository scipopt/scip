find_path(CONOPT_INCLUDE_DIRS
    NAMES conopt.h
    HINTS ${CONOPT_DIR} $ENV{CONOPT_DIR}
    PATH_SUFFIXES include)

find_library(CONOPT_LIBRARY
    NAMES conopt
    HINTS ${CONOPT_DIR} $ENV{CONOPT_DIR}
    PATH_SUFFIXES lib)

set(CONOPT_LIBRARIES ${CONOPT_LIBRARY})

if(DEFINED ENV{CONOPT_LICENSE_INT_1} AND DEFINED ENV{CONOPT_LICENSE_INT_2} AND DEFINED ENV{CONOPT_LICENSE_INT_3} AND DEFINED ENV{CONOPT_LICENSE_TEXT})
   message("CONOPT license found")
   add_compile_definitions(CONOPT_LICENSE_INT_1=$ENV{CONOPT_LICENSE_INT_1})
   add_compile_definitions(CONOPT_LICENSE_INT_2=$ENV{CONOPT_LICENSE_INT_2})
   add_compile_definitions(CONOPT_LICENSE_INT_3=$ENV{CONOPT_LICENSE_INT_3})
   add_compile_definitions(CONOPT_LICENSE_TEXT=\"$ENV{CONOPT_LICENSE_TEXT}\")
else()
   message("CONOPT license not found. The license needs to be specified in environment variables CONOPT_LICENSE_INT_1, CONOPT_LICENSE_INT_2, CONOPT_LICENSE_INT_3 and CONOPT_LICENSE_TEXT.")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CONOPT DEFAULT_MSG CONOPT_INCLUDE_DIRS CONOPT_LIBRARIES)

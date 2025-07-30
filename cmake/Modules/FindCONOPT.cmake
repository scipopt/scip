find_path(CONOPT_INCLUDE_DIRS
    NAMES coidef.f90
    HINTS ${CONOPT_DIR} $ENV{CONOPT_DIR}
    PATH_SUFFIXES include)

find_library(CONOPT_LIBRARY
    NAMES conopt4
    HINTS ${CONOPT_DIR} $ENV{CONOPT_DIR}
    PATH_SUFFIXES lib)

set(CONOPT_LIBRARIES ${CONOPT_LIBRARY})

if(DEFINED ENV{LICENSE_INT_1} AND DEFINED ENV{LICENSE_INT_2} AND DEFINED ENV{LICENSE_INT_3} AND DEFINED ENV{LICENSE_TEXT})
   message("CONOPT license found")
   add_compile_definitions(LICENSE_INT_1=$ENV{LICENSE_INT_1})
   add_compile_definitions(LICENSE_INT_2=$ENV{LICENSE_INT_2})
   add_compile_definitions(LICENSE_INT_3=$ENV{LICENSE_INT_3})
   add_compile_definitions(LICENSE_TEXT=\"$ENV{LICENSE_TEXT}\")
else()
   message("CONOPT license not found. The license needs to be specified in environment variables LICENSE_INT_1, LICENSE_INT_2, LICENSE_INT_3 and LICENSE_TEXT.")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CONOPT DEFAULT_MSG CONOPT_INCLUDE_DIRS CONOPT_LIBRARIES)

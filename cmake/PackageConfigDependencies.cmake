
function(package_config_dependency_dir result package_found package_dir)
   # Writes CMake code to ${result} for checking ${package_found} and ${package_dir} variables
   # for whether ${package_dir} at build time is different from a variable specified by the user.

   # Emit a warning if ${package}_DIR is different from value at build time.
   set(code "")
   if(DEFINED ENV{${package_dir}})
      set(current_value $ENV{${package_dir}})
      list(APPEND code "if(DEFINED ENV{${package_dir}} AND NOT (\"\$ENV{${package_dir}}\" STREQUAL \"${current_value}\"))")
      list(APPEND code "   message(WARNING \"Environment variable ${package_dir} is \\\"\$ENV{${package_dir}}\\\", but was \\\"${current_value}\\\" at build time.\")")
      list(APPEND code "endif()")
   elseif(${package_dir})
      set(current_value ${${package_dir}})
      list(APPEND code "if(${package_dir} AND NOT (${package_dir} STREQUAL \"${current_value}\"))")
      list(APPEND code "   message(WARNING \"CMake variable ${package_dir} is \\\"\$${package_dir}\\\", but was \\\"${current_value}\\\" at build time.\")")
      list(APPEND code "endif()")
   else()
      set(${result} "" PARENT_SCOPE)
      return()
   endif()
   # If ${package}_DIR is not given then we set it to the build time value.
   list(APPEND code "if(NOT ${package_found} AND NOT DEFINED ENV{${package_dir}} AND NOT ${package_dir})")
   list(APPEND code "   set(${package_dir} \"${current_value}\")")
   list(APPEND code "endif()")
   string(REPLACE ";" "\n      " code "${code}")
   set(${result} "${code}" PARENT_SCOPE)
endfunction(package_config_dependency_dir)


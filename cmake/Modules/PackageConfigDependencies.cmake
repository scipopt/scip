
function(package_config_dependency_dir result package_found package_dir)
   # Writes CMake code to ${result} for checking ${package_found} and ${package_dir} variables
   # for whether ${package_dir} at build time is different from a variable specified by the user.

   # Emit a warning if ${package}_DIR is different from value at build time.
   set(code "if(${package_dir} AND NOT ${package_dir} STREQUAL \"${${package_dir}}\")")
   set(code "${code}\n   message(WARNING \"${package_dir} is \\\"\$\{${package_dir}\}\\\", but was \\\"${${package_dir}}\\\" at build time.\")")
   set(code "${code}\nendif()")

   # If ${package}_DIR is not given then we set it to the build time value.
   set(code "${code}\nif(NOT ${package_found} AND NOT ${package_dir})")
   set(code "${code}\n   set(${package_dir} \"${${package_dir}}\")")
   set(code "${code}\nendif()")

   set(${result} "${code}" PARENT_SCOPE)
endfunction(package_config_dependency_dir)


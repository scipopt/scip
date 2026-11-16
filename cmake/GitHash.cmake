# The MIT License (MIT)
#
# Copyright (c)
#   2026 Matthias Walter
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

function(githash_query RETURN_VAR FILE MACRO APPEND_BRANCH)
   # Search for git.
   find_program(GIT git)
   if(GIT AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
      execute_process(
         COMMAND ${GIT} describe --always --dirty
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
         OUTPUT_VARIABLE _githash_query_result OUTPUT_STRIP_TRAILING_WHITESPACE
         )
      string(REGEX REPLACE "^.*-g" "" _githash_query_result ${_githash_query_result})
      if(APPEND_BRANCH)
         execute_process(
            COMMAND ${GIT} branch --show-current
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            OUTPUT_VARIABLE _githash_branch_result OUTPUT_STRIP_TRAILING_WHITESPACE
            )
         set(${RETURN_VAR} ${_githash_query_result}-${_githash_branch_result} PARENT_SCOPE)
      else()
         set(${RETURN_VAR} ${_githash_query_result} PARENT_SCOPE)
      endif()
   else()
      # Otherwise we check FILE for a C/C++ macro definition of MACRO.
      file(STRINGS "${FILE}" _githash_query_define REGEX "define ${MACRO} .*")
      if(_githash_query_define)
         string(REGEX MATCH "\\\"(.*)\\\"" _githash_query_matched ${_githash_query_define})
         if(_githash_query_matched)
            set(${RETURN_VAR} ${CMAKE_MATCH_1} PARENT_SCOPE)
         endif()
      endif()
   endif()
endfunction()

function(githash_add_target TARGET FILE MACRO APPEND_BRANCH)
   # Create a CMake script that runs `git describe --always --dirty` and
   # compares the resulting git hash with that defined by MACRO in FILE.
   # If they differ, FILE is updated.
   file(WRITE ${CMAKE_BINARY_DIR}/${TARGET}.cmake "
find_program(GIT git)
if(EXISTS \"${FILE}\")
   file(STRINGS \"${FILE}\" _githash_add_target_define REGEX \"define ${MACRO} .*\")
   string(REGEX MATCH \"\\\"(.*)\\\"\" _githash_add_target_matched \${_githash_add_target_define})
   if(_githash_add_target_matched)
      set(_githash_add_target_old \${CMAKE_MATCH_1})
   endif()
endif()
if(GIT AND EXISTS \"${CMAKE_CURRENT_SOURCE_DIR}/.git\")
   execute_process(
      COMMAND \${GIT} describe --always --dirty
      WORKING_DIRECTORY \"${CMAKE_CURRENT_SOURCE_DIR}\"
      OUTPUT_VARIABLE _githash_add_target_githash OUTPUT_STRIP_TRAILING_WHITESPACE
      )
   if(${APPEND_BRANCH})
      execute_process(
         COMMAND \${GIT} branch --show-current
         WORKING_DIRECTORY \"${CMAKE_CURRENT_SOURCE_DIR}\"
         OUTPUT_VARIABLE _githash_branch_result OUTPUT_STRIP_TRAILING_WHITESPACE
         )
      string(REGEX REPLACE \"^.*-g\" \"\" _githash_query_result \${_githash_add_target_githash})
      set(_githash_add_target_githash \${_githash_add_target_githash}-\${_githash_branch_result})
   else()
      string(REGEX REPLACE \"^.*-g\" \"\" _githash_add_target_githash \${_githash_add_target_githash})
   endif()
else()
   set(_githash_add_target_githash \${_githash_add_target_old})
endif()
if(NOT _githash_add_target_githash)
   message(WARNING \"No git information available!\")
elseif(\"\${_githash_add_target_githash}\" STREQUAL \"\${_githash_add_target_old}\")
   message(VERBOSE \"No update of githash in \\\"${FILE}\\\" needed.\" )
else()
   file(WRITE \"${FILE}\" \"#define ${MACRO} \\\"\${_githash_add_target_githash}\\\"\n\")
   message(STATUS \"Updating githash \\\"\${_githash_add_target_githash}\\\" in \\\"${FILE}\\\".\")
endif()
#")

   # Add a target that executes the created file. 
   add_custom_target(${TARGET}
      COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/${TARGET}.cmake)
endfunction()

# Cmake file for comparing two status files
#
# - We read the files FILE1 and FILE2 into cmake variables.
# - Then all floating point numbers are eliminated to get rid of variances in running time etc.
# - The results are writen back to files and compared.
#
# message(STATUS "Comparing ${FILE1} and ${FILE2}")
# read FILE1 into variable R1
file(STRINGS ${FILE1} R1)
# remove all occurence of floating point numbers
string(REGEX REPLACE "[0-9]+[.][0-9]+" "" R1OUT "${R1}")
# write result into file
file(WRITE ${FILE1}.s ${R1OUT})
#
# similarly for FILE2
file(STRINGS ${FILE2} R2)
string(REGEX REPLACE "[0-9]+[.][0-9]+" "" R2OUT "${R2}")
file(WRITE ${FILE2}.s ${R2OUT})
#
# finally compare the two resulting files
execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files ${FILE1}.s ${FILE2}.s)

#
# define the instance sets
#
# semicolon '\;' is used to split an instance and its optimal objective value
# For infeasible instances, '+infinity' is used (or '-infinity' in case of maximization)
#

set(instances
    "ex1\;3100"
#     "ex2\;3224"
# this example takes 700 seconds
    "t65i11xx\;16389651"
    "t70x11xx\;343471236"
)

#
# loop over the instances
#
foreach(instance ${instances})
    list(GET instance 0 basename)
    list(GET instance 1 optval)
        #
        # treat the instance as a tuple (list) of two values
        #
        add_test(NAME examples-lop-${basename} COMMAND $<TARGET_FILE:lop>  -f ${CMAKE_CURRENT_SOURCE_DIR}/../data/${basename}.lop -o ${optval} ${optval})
        set_tests_properties(examples-lop-${basename} PROPERTIES PASS_REGULAR_EXPRESSION "Validation         : Success")
endforeach(instance)
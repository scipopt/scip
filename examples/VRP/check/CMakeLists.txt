include(CTest)
#
# define the instance sets
#
# semicolon '\;' is used to split an instance and its optimal objective value
# For infeasible instances, '+infinity' is used (or '-infinity' in case of maximization)
#
set(instances
    "eil13\;247"
    "eil22\;375"
    "eil7\;104"
)

#
# loop over the instances
#
foreach(instance ${instances})
    #
    # treat the instance as a tuple (list) of two values
    #
    list(GET instance 0 basename)
    list(GET instance 1 optval)
    add_test(NAME examples-vrp-${basename} COMMAND $<TARGET_FILE:vrp> ${CMAKE_CURRENT_SOURCE_DIR}/../data/${basename}.vrp)
endforeach(instance)
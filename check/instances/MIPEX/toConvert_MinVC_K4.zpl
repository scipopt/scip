# Minimum vertex cover for a K4 graph
# We have a set of edges and a set of vertex
# Each variable is represented by a vertex
#
var x_v1 binary;
var x_v2 binary;
var x_v3 binary;
var x_v4 binary;

#Objective Function:
#
minimize cost: x_v1 + x_v2 + x_v3 + x_v4;

#Constraints:
#
subto E1: x_v1 + x_v2 >= 1;
subto E2: x_v1 + x_v3 >= 1;
subto E3: x_v2 + x_v3 >= 1;
subto E4: x_v1 + x_v4 >= 1;
subto E5: x_v2 + x_v4 >= 1;
subto E6: x_v3 + x_v4 >= 1; 


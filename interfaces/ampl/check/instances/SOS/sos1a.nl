g3 2 1 0	# problem sos1a
 3 1 1 0 0	# vars, constraints, objectives, ranges, eqns
 0 0	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 0 0 0	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 3 3	# nonzeros in Jacobian, gradients
 3 2	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
S0 3 sosno
0 1
1 1
2 1
S0 3 ref
0 1
1 2
2 3
C0
n0
O0 1
n0
r
1 1
b
0 0 0.8
0 0 0.6
0 0 0.6
k2
1
2
J0 3
0 1
1 1
2 1
G0 3
0 0.9
1 1
2 1.1

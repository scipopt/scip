g3 0 1 0	# problem suffix1
 3 2 1 0 0	# vars, constraints, objectives, ranges, eqns
 0 0	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 0 0 0	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 5 3	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
S0 1 initial
0 2
S0 1 removable
0 1
S0 2 sosno
0 1
1 1
S1 1 initial
1 2
S1 1 separate
1 2
S1 1 enforce
1 2
S1 1 check
1 2
S1 1 propagate
1 2
S1 1 dynamic
1 1
S1 1 removable
1 1
C0
n0
C1
n0
O0 1
n0
r
1 10
1 0
b
2 0
2 0
0 0 100
k2
2
4
J0 2
0 1
1 1
J1 3
0 1
1 1
2 -1
G0 3
0 1
1 1
2 1

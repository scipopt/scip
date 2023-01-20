g3 0 1 0	# problem sos2a
 7 5 1 0 3	# vars, constraints, objectives, ranges, eqns
 0 0	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 0 0 0	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 15 2	# nonzeros in Jacobian, gradients
 3 2	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
S0 3 sosno
0 -1
1 -1
2 -1
S0 3 ref
0 1
1 2
2 3
C0
n0
C1
n0
C2
n0
C3
n0
C4
n0
O0 0
n0
r
4 1
4 0
4 0
2 -1.3
2 1.3
b
2 0.8
2 0
2 0
2 0
2 0
3
3
k6
3
6
9
10
11
12
J0 3
0 1
1 1
2 1
J1 4
0 1
1 2
2 3
5 -1
J2 4
0 1
1 2
2 3
6 -1
J3 2
3 1
6 -1
J4 2
4 1
6 1
G0 2
3 1
4 1

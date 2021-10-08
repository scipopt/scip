*  This test problem has 3 constraints and 3 variables
* the constraints are:
* c0: x1 + 9 x1^2 + 0.5 x1 * x2 + 0.5 x2 * x1 = 0.1
* c1: x2 + 2 x3^2 >= 0.2
* c2: 0.5 x1 + x1^2 + x2^2 - x3^2 <= 0.3
*
* the objective function is quadratic as well 2x1*x2 + 0.1x2*x3 -0.1 x3*x1
* It gets transformed into a constraint 2x1*x2 + 0.1x2*x3 -0.1 x3*x1 <= qmatrixvar
NAME          test.mps
ROWS
 N  Obj
 E  c0
 L  c1
 G  c2
COLUMNS
    x1        c0                               1  c2                             0.5
    x2        c1                               1
RHS
    RHS       c0                             0.1  c1                             0.2
    RHS       c2                             0.3
BOUNDS
 FR Bound     x1
 FR Bound     x2
 FR Bound     x3
QMATRIX
    x1        x2                               2
    x2        x3                             0.1
    x3        x1                            -0.1
QCMATRIX c0
    x1        x1                               9
    x1        x2                             0.5
    x2        x1                             0.5
QCMATRIX c1
    x3        x3                               2
QCMATRIX c2
    x1        x1                               1
    x2        x2                               1
    x3        x3                              -1
ENDATA

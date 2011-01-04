* SCIP STATISTICS
*   Problem name     : sep1.zpl
*   Variables        : 30 (2 binary, 0 integer, 0 implicit integer, 28 continuous)
*   Constraints      : 32
*   Obj. scale       : 1
*   Obj. offset      : 0
NAME          sep1.zpl
ROWS
 N  Obj 
 E  c0 
 E  c1 
 E  c2 
 E  c3 
 E  c4 
 E  c5 
 E  c6 
 E  c7 
 E  c8 
 E  c9 
 E  c10 
 E  c11 
 E  c12 
 E  c13 
 E  c14 
 E  c15 
 E  c16 
 E  c17 
 E  c18 
 E  c19 
 E  c20 
 E  c21 
 G  c22 
 L  c23 
 G  c24 
 L  c25 
 G  c26 
 G  c27 
 L  c28 
 L  c29 
 G  c30 
 E  c31 
COLUMNS
    x27       c18                              1  c28                              1 
    x27       c26                              1  c31                             35 
    x28       c19                              1  c26                             -4 
    x28       c28                              1 
    x2        c20                              1  c27                             -3 
    x2        c29                              1 
    x3        c27                              1  c21                              1 
    x3        c31                             30  c29                              1 
    x4        c0                           -0.55  c1                           -0.45 
    x4        c31                            -10 
    x5        c1                            -0.5  c0                            -0.5 
    x5        c31                             -8 
    x6        c8                              -1  c0                               1 
    x7        c1                               1  c9                              -1 
    x8        c31                             -1  c2                               1 
    x8        c23                              1  c8                               1 
    x8        c12                          -0.15  c10                          -0.85 
    x8        c22                              1 
    x9        c9                               1  c3                               1 
    x9        c31                             -1  c13                           -0.8 
    x9        c11                           -0.2  c22                              1 
    x9        c23                              1 
    x10       c8                               1  c24                              1 
    x10       c25                              1  c16                         -0.025 
    x10       c14                         -0.975  c4                               1 
    x10       c31                             -4 
    x11       c15                          -0.05  c24                              1 
    x11       c5                               1  c17                          -0.95 
    x11       c9                               1  c31                             -4 
    x11       c25                              1 
    x12       c6                               1  c18                             -1 
    x12       c8                               1 
    x13       c19                             -1  c9                               1 
    x13       c7                               1 
    x14       c8                               1  c20                             -1 
    x15       c9                               1  c21                             -1 
    x16       c18                             -1  c10                              1 
    x17       c19                             -1  c11                              1 
    x18       c20                             -1  c12                              1 
    x19       c21                             -1  c13                              1 
    x20       c14                              1  c18                             -1 
    x21       c19                             -1  c15                              1 
    x22       c16                              1  c20                             -1 
    x23       c21                             -1  c17                              1 
    INTSTART  'MARKER'                            'INTORG'                           
    x0        c30                              1  c25                            -25 
    x0        c24                           -2.5  c31                            -50 
    x1        c30                              1  c23                            -25 
    x1        c22                           -2.5  c31                             -2 
    INTEND    'MARKER'                            'INTEND'                           
    x29       Obj                              1  c31                              1 
RHS
    RHS       c28                             15  c29                             18 
    RHS       c30                              1 
BOUNDS
 BV Bound     x0                                 
 BV Bound     x1                                 
 UP Bound     x2                              50 
 UP Bound     x3                              50 
 UP Bound     x4                              25 
 UP Bound     x5                              25 
 UP Bound     x6                              50 
 UP Bound     x7                              50 
 UP Bound     x8                              50 
 UP Bound     x9                              50 
 UP Bound     x10                             50 
 UP Bound     x11                             50 
 UP Bound     x12                             50 
 UP Bound     x13                             50 
 UP Bound     x14                             50 
 UP Bound     x15                             50 
 UP Bound     x16                             50 
 UP Bound     x17                             50 
 UP Bound     x18                             50 
 UP Bound     x19                             50 
 UP Bound     x20                             50 
 UP Bound     x21                             50 
 UP Bound     x22                             50 
 UP Bound     x23                             50 
 UP Bound     x24                              1 
 UP Bound     x25                              1 
 UP Bound     x26                              1 
 UP Bound     x27                             50 
 UP Bound     x28                             50 
 FR Bound     x29                                
QCMATRIX c2
    x24       x6                            -0.5 
    x6        x24                           -0.5 
QCMATRIX c3
    x24       x7                            -0.5 
    x7        x24                           -0.5 
QCMATRIX c4
    x25       x6                            -0.5 
    x6        x25                           -0.5 
QCMATRIX c5
    x25       x7                            -0.5 
    x7        x25                           -0.5 
QCMATRIX c6
    x26       x6                            -0.5 
    x6        x26                           -0.5 
QCMATRIX c7
    x26       x7                            -0.5 
    x7        x26                           -0.5 
ENDATA
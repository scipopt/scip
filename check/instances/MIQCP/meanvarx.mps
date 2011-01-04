* SCIP STATISTICS
*   Problem name     : meanvarx.zpl
*   Variables        : 36 (14 binary, 0 integer, 0 implicit integer, 22 continuous)
*   Constraints      : 45
*   Obj. scale       : 1
*   Obj. offset      : 0
NAME          meanvarx.zpl
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
 L  c8 
 L  c9 
 L  c10 
 L  c11 
 L  c12 
 L  c13 
 L  c14 
 L  c15 
 G  c16 
 G  c17 
 G  c18 
 G  c19 
 G  c20 
 G  c21 
 G  c22 
 L  c23 
 L  c24 
 L  c25 
 L  c26 
 L  c27 
 L  c28 
 L  c29 
 G  c30 
 G  c31 
 G  c32 
 G  c33 
 G  c34 
 G  c35 
 G  c36 
 L  c37 
 L  c38 
 L  c39 
 L  c40 
 L  c41 
 L  c42 
 L  c43 
 E  c44 
COLUMNS
    x22       Obj                              1  c44                              1 
    x23       c0                               1  c44                        0.06435 
    x23       c1                               1 
    x24       c2                               1  c0                               1 
    x24       c44                         0.0548 
    x25       c0                               1  c3                               1 
    x25       c44                        0.02505 
    x26       c4                               1  c0                               1 
    x26       c44                         0.0762 
    x27       c0                               1  c44                        0.03815 
    x27       c5                               1 
    x28       c0                               1  c44                         0.0927 
    x28       c6                               1 
    x29       c0                               1  c44                          0.031 
    x29       c7                               1 
    x30       c1                              -1  c16                              1 
    x30       c9                               1  c8                               1 
    x31       c17                              1  c10                              1 
    x31       c8                               1  c2                              -1 
    x32       c18                              1  c8                               1 
    x32       c3                              -1  c11                              1 
    x33       c8                               1  c19                              1 
    x33       c12                              1  c4                              -1 
    x34       c8                               1  c20                              1 
    x34       c5                              -1  c13                              1 
    x35       c6                              -1  c21                              1 
    x35       c8                               1  c14                              1 
    x14       c22                              1  c15                              1 
    x14       c7                              -1  c8                               1 
    x15       c1                               1  c30                              1 
    x15       c23                              1 
    x16       c31                              1  c2                               1 
    x16       c24                              1 
    x17       c3                               1  c32                              1 
    x17       c25                              1 
    x18       c4                               1  c33                              1 
    x18       c26                              1 
    x19       c27                              1  c34                              1 
    x19       c5                               1 
    x20       c6                               1  c35                              1 
    x20       c28                              1 
    x21       c7                               1  c36                              1 
    x21       c29                              1 
    INTSTART  'MARKER'                            'INTORG'                           
    x0        c9                           -0.11  c37                              1 
    x0        c16                          -0.03 
    x1        c10                           -0.1  c17                          -0.04 
    x1        c38                              1 
    x2        c39                              1  c11                          -0.07 
    x2        c18                          -0.04 
    x3        c12                          -0.11  c19                          -0.03 
    x3        c40                              1 
    x4        c41                              1  c20                          -0.03 
    x4        c13                           -0.2 
    x5        c14                           -0.1  c21                          -0.03 
    x5        c42                              1 
    x6        c15                           -0.1  c22                          -0.03 
    x6        c43                              1 
    x7        c23                           -0.2  c37                              1 
    x7        c30                          -0.02 
    x8        c38                              1  c24                          -0.15 
    x8        c31                          -0.02 
    x9        c32                          -0.04  c39                              1 
    x10       c33                          -0.04  c40                              1 
    x11       c41                              1  c34                          -0.04 
    x11       c27                           -0.1 
    x12       c28                          -0.15  c42                              1 
    x12       c35                          -0.04 
    x13       c43                              1  c36                          -0.04 
    x13       c29                           -0.2 
RHS
    RHS       c0                               1  c1                             0.2 
    RHS       c2                             0.2  c5                             0.2 
    RHS       c6                             0.2  c7                             0.2 
    RHS       c8                             0.3  c37                              1 
    RHS       c38                              1  c39                              1 
    RHS       c40                              1  c41                              1 
    RHS       c42                              1  c43                              1 
BOUNDS
 BV Bound     x0                                 
 BV Bound     x1                                 
 BV Bound     x2                                 
 BV Bound     x3                                 
 BV Bound     x4                                 
 BV Bound     x5                                 
 BV Bound     x6                                 
 BV Bound     x7                                 
 BV Bound     x8                                 
 BV Bound     x9                                 
 BV Bound     x10                                
 BV Bound     x11                                
 BV Bound     x12                                
 BV Bound     x13                                
 PL Bound     x14                                
 PL Bound     x15                                
 PL Bound     x16                                
 PL Bound     x17                                
 PL Bound     x18                                
 PL Bound     x19                                
 PL Bound     x20                                
 PL Bound     x21                                
 FR Bound     x22                                
 PL Bound     x23                                
 PL Bound     x24                                
 PL Bound     x25                                
 PL Bound     x26                                
 PL Bound     x27                                
 PL Bound     x28                                
 PL Bound     x29                                
 PL Bound     x30                                
 PL Bound     x31                                
 PL Bound     x32                                
 PL Bound     x33                                
 PL Bound     x34                                
 PL Bound     x35                                
QCMATRIX c44
    x23       x23                         -42.18 
    x24       x24                         -70.89 
    x25       x25                         -25.51 
    x26       x26                         -22.33 
    x27       x27                         -30.01 
    x28       x28                         -42.23 
    x29       x29                         -16.42 
    x24       x23                         -20.18 
    x23       x24                         -20.18 
    x25       x23                         -10.88 
    x23       x25                         -10.88 
    x26       x23                           -5.3 
    x23       x26                           -5.3 
    x27       x23                         -12.32 
    x23       x27                         -12.32 
    x28       x23                         -23.84 
    x23       x28                         -23.84 
    x29       x23                         -17.41 
    x23       x29                         -17.41 
    x25       x24                         -21.58 
    x24       x25                         -21.58 
    x26       x24                         -15.41 
    x24       x26                         -15.41 
    x27       x24                         -23.24 
    x24       x27                         -23.24 
    x28       x24                          -23.8 
    x24       x28                          -23.8 
    x29       x24                         -12.62 
    x24       x29                         -12.62 
    x26       x25                           -9.6 
    x25       x26                           -9.6 
    x27       x25                         -22.63 
    x25       x27                         -22.63 
    x28       x25                         -13.22 
    x25       x28                         -13.22 
    x29       x25                           -4.7 
    x25       x29                           -4.7 
    x27       x26                         -10.32 
    x26       x27                         -10.32 
    x28       x26                         -10.46 
    x26       x28                         -10.46 
    x29       x26                             -1 
    x26       x29                             -1 
    x28       x27                         -16.36 
    x27       x28                         -16.36 
    x29       x27                           -7.2 
    x27       x29                           -7.2 
    x29       x28                           -9.9 
    x28       x29                           -9.9 
ENDATA
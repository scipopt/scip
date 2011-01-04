* SCIP STATISTICS
*   Problem name     : tltr.zpl
*   Variables        : 49 (12 binary, 36 integer, 0 implicit integer, 1 continuous)
*   Constraints      : 55
*   Obj. scale       : 1
*   Obj. offset      : 0
NAME          tltr.zpl
ROWS
 N  Obj 
 E  c0 
 G  c1 
 G  c2 
 G  c3 
 L  c4 
 L  c5 
 L  c6 
 L  c7 
 L  c8 
 L  c9 
 L  c10 
 L  c11 
 L  c12 
 L  c13 
 L  c14 
 L  c15 
 L  c16 
 L  c17 
 L  c18 
 L  c19 
 L  c20 
 L  c21 
 L  c22 
 L  c23 
 L  c24 
 L  c25 
 L  c26 
 L  c27 
 L  c28 
 L  c29 
 L  c30 
 L  c31 
 L  c32 
 L  c33 
 L  c34 
 L  c35 
 L  c36 
 L  c37 
 L  c38 
 L  c39 
 L  c40 
 L  c41 
 L  c42 
 L  c43 
 L  c44 
 L  c45 
 L  c46 
 L  c47 
 L  c48 
 L  c49 
 L  c50 
 L  c51 
 L  c52 
 L  c53 
 L  c54 
COLUMNS
    INTSTART  'MARKER'                            'INTORG'                           
    x0        c40                           -500 
    x1        c0                             -35  c41                          -1270 
    x2        c42                          -1270  c0                             -35 
    x30       c13                             -1  c4                              12 
    x30       c31                              1 
    x31       c14                             -1  c32                              1 
    x31       c5                              12 
    x32       c6                              12  c33                              1 
    x32       c15                             -1 
    x33       c16                             -1  c7                              12 
    x33       c34                              1 
    x34       c35                              1  c17                             -1 
    x34       c8                              12 
    x35       c9                              12  c36                              1 
    x35       c18                             -1 
    x36       c19                             -1  c10                             12 
    x36       c37                              1 
    x37       c20                             -1  c11                             12 
    x37       c38                              1 
    x38       c12                             12  c21                             -1 
    x38       c39                              1 
    x12       c13                             -1  c31                              1 
    x12       c4                              24 
    x13       c5                              24  c32                              1 
    x13       c14                             -1 
    x14       c33                              1  c15                             -1 
    x14       c6                              24 
    x15       c7                              24  c34                              1 
    x15       c16                             -1 
    x16       c35                              1  c17                             -1 
    x16       c8                              24 
    x17       c9                              24  c18                             -1 
    x17       c36                              1 
    x18       c10                             24  c19                             -1 
    x18       c37                              1 
    x19       c20                             -1  c11                             24 
    x19       c38                              1 
    x20       c12                             24  c39                              1 
    x20       c21                             -1 
    x21       c13                             -1  c4                              36 
    x21       c31                              1 
    x22       c5                              36  c14                             -1 
    x22       c32                              1 
    x23       c15                             -1  c33                              1 
    x23       c6                              36 
    x24       c7                              36  c34                              1 
    x24       c16                             -1 
    x25       c8                              36  c17                             -1 
    x25       c35                              1 
    x26       c9                              36  c18                             -1 
    x26       c36                              1 
    x27       c19                             -1  c37                              1 
    x27       c10                             36 
    x28       c38                              1  c20                             -1 
    x28       c11                             36 
    x29       c21                             -1  c39                              1 
    x29       c12                             36 
    x3        c22                            -72  c43                             -1 
    x3        c13                              1  c31                             -5 
    x3        c0               -6.53333333333333  c4                             -48 
    x4        c23                           -182  c32                             -5 
    x4        c0               -6.53333333333333  c5                             -48 
    x4        c14                              1  c44                             -1 
    x5        c0                         -6.7375  c45                             -1 
    x5        c6                             -62  c15                              1 
    x5        c24                           -182  c33                             -5 
    x6        c43                              1  c0               -6.53333333333333 
    x6        c16                              1  c25                            -72 
    x6        c46                             -1  c7                             -48 
    x6        c34                             -5 
    x7        c47                             -1  c44                              1 
    x7        c17                              1  c26                           -182 
    x7        c8                             -48  c35                             -5 
    x7        c0               -6.53333333333333 
    x8        c36                             -5  c9                             -62 
    x8        c18                              1  c48                             -1 
    x8        c27                           -182  c45                              1 
    x8        c0                         -6.7375 
    x9        c46                              1  c0               -6.53333333333333 
    x9        c37                             -5  c28                            -72 
    x9        c10                            -48  c19                              1 
    x10       c29                           -182  c38                             -5 
    x10       c0               -6.53333333333333  c11                            -48 
    x10       c47                              1  c20                              1 
    x11       c39                             -5  c30                           -182 
    x11       c0                         -6.7375  c12                            -62 
    x11       c21                              1  c48                              1 
    x39       c49                             -1  c40                              7 
    x39       c22                              1 
    x40       c41                              7  c50                             -1 
    x40       c23                              1 
    x41       c42                              7  c24                              1 
    x41       c51                             -1 
    x42       c40                              7  c49                              1 
    x42       c52                             -1  c25                              1 
    x43       c41                              7  c26                              1 
    x43       c50                              1  c53                             -1 
    x44       c42                              7  c27                              1 
    x44       c51                              1  c54                             -1 
    x45       c40                              7  c52                              1 
    x45       c28                              1 
    x46       c41                              7  c29                              1 
    x46       c53                              1 
    x47       c54                              1  c30                              1 
    x47       c42                              7 
    INTEND    'MARKER'                            'INTEND'                           
    x48       c0                               1  Obj                              1 
RHS
    RHS       c1                               9  c2                              15 
    RHS       c3                              80 
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
 UP Bound     x12                              5 
 UP Bound     x13                              5 
 UP Bound     x14                              5 
 UP Bound     x15                              5 
 UP Bound     x16                              5 
 UP Bound     x17                              5 
 UP Bound     x18                              5 
 UP Bound     x19                              5 
 UP Bound     x20                              5 
 UP Bound     x21                              5 
 UP Bound     x22                              5 
 UP Bound     x23                              5 
 UP Bound     x24                              5 
 UP Bound     x25                              5 
 UP Bound     x26                              5 
 UP Bound     x27                              5 
 UP Bound     x28                              5 
 UP Bound     x29                              5 
 UP Bound     x30                              5 
 UP Bound     x31                              5 
 UP Bound     x32                              5 
 UP Bound     x33                              5 
 UP Bound     x34                              5 
 UP Bound     x35                              5 
 UP Bound     x36                              5 
 UP Bound     x37                              5 
 UP Bound     x38                              5 
 UP Bound     x39                            100 
 UP Bound     x40                            100 
 UP Bound     x41                            100 
 UP Bound     x42                            100 
 UP Bound     x43                            100 
 UP Bound     x44                            100 
 UP Bound     x45                            100 
 UP Bound     x46                            100 
 UP Bound     x47                            100 
 FR Bound     x48                                
QCMATRIX c1
    x39       x30                            0.5 
    x30       x39                            0.5 
    x42       x33                            0.5 
    x33       x42                            0.5 
    x45       x36                            0.5 
    x36       x45                            0.5 
    x40       x31                            0.5 
    x31       x40                            0.5 
    x43       x34                            0.5 
    x34       x43                            0.5 
    x46       x37                            0.5 
    x37       x46                            0.5 
    x41       x32                            0.5 
    x32       x41                            0.5 
    x44       x35                            0.5 
    x35       x44                            0.5 
    x47       x38                            0.5 
    x38       x47                            0.5 
QCMATRIX c2
    x39       x12                            0.5 
    x12       x39                            0.5 
    x42       x15                            0.5 
    x15       x42                            0.5 
    x45       x18                            0.5 
    x18       x45                            0.5 
    x40       x13                            0.5 
    x13       x40                            0.5 
    x43       x16                            0.5 
    x16       x43                            0.5 
    x46       x19                            0.5 
    x19       x46                            0.5 
    x41       x14                            0.5 
    x14       x41                            0.5 
    x44       x17                            0.5 
    x17       x44                            0.5 
    x47       x20                            0.5 
    x20       x47                            0.5 
QCMATRIX c3
    x39       x21                            0.5 
    x21       x39                            0.5 
    x42       x24                            0.5 
    x24       x42                            0.5 
    x45       x27                            0.5 
    x27       x45                            0.5 
    x40       x22                            0.5 
    x22       x40                            0.5 
    x43       x25                            0.5 
    x25       x43                            0.5 
    x46       x28                            0.5 
    x28       x46                            0.5 
    x41       x23                            0.5 
    x23       x41                            0.5 
    x44       x26                            0.5 
    x26       x44                            0.5 
    x47       x29                            0.5 
    x29       x47                            0.5 
ENDATA
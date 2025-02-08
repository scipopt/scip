Minimize
 obj: cost
Subject To
 c1:  - cost + 80 x1 + 60 x2 + 55 x3 + 30 x4 + 25 x5 + 80 x6 + 60 x7 + 35 x8
      + 80 x9 + 55 x10  = 0
 c2:  x1 + x2 + 0.8 x3 + 0.6 x4 + 0.4 x5 >= 2.1
 c3:  x6 + 0.9 x7 + 0.5 x8 >= 1.2
 c4:  x9 + 0.9 x10 >= 0.8
 c5:  0.2 x2 + x3 + 0.5 x4 + 0.5 x5 + 0.2 x7 + 0.5 x8 + x10 - service  = 0
 c6:  x1 + x6 + x9 >= 1
 c7:  x1 + x2 + x3 + x6 + x7 + x9 >= 2
 c8:  x2 + x3 + x4 + x5 <= 0
 c9:  x4 + x5 + x8 <= 1
 c10: x1 + x10 <= 1
Bounds
 service >= 3.2
Binaries
 x1  x2  x3  x4  x5  x6  x7  x8  x9  x10
End

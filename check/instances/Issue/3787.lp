\ SCIP STATISTICS
\   Problem name     : sparse2
\   Variables        : 6 (2 binary, 0 integer, 0 implicit integer, 4 continuous)
\   Constraints      : 5
Maximize
 Obj: +1 z#8 +1 z#10
Subject to
 rhs_2: +0.25 x#28 +0.25 x#18 >= +0.01821599166123
 rhs_8: -0.25 x#28 +0.25 x#18 -1 x#8 = +1.08811956920877
 lhs_10: -0.25 x#28 +0.25 x#18 -1 x#10 <= +0.23636262424061
Bounds
 0 <= z#8 <= 1
 0 <= z#10 <= 1
 x#28 free
 x#18 free
 x#8 free
 x#10 free
Binaries
 z#8 z#10
SOS
 SOS8: S1:: z#8:1 x#8:2
 SOS10: S1:: z#10:1 x#10:2
End

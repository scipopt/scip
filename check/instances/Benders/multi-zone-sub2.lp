\ Ring-Structured Capacity Expansion
\ Master Variables: y1, y2, y3 (Facility Investments)

Minimize
  obj: 4 x22 + 9 x32  \ Zone 2 served by Facility 2 and 3

Subject To
  \ Demand 2 can be met by Facility 2 OR Facility 3
  demand_z2: x22 + x32 >= 70

  \ Facility 2 serves Zone 1 and Zone 2
  cap_f2: x21 + x22 - 100 y2 <= 0

Bounds
  0 <= x22
  0 <= x32
  0 <= y2 <= 1

Generals
  y2

End

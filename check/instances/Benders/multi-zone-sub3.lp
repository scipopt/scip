\ Ring-Structured Capacity Expansion
\ Master Variables: y1, y2, y3 (Facility Investments)

Minimize
  obj: 6 x33 + 7 x13  \ Zone 3 served by Facility 3 and 1

Subject To
  \ Demand 3 can be met by Facility 3 OR Facility 1
  demand_z3: x33 + x13 >= 40

  \ Facility 3 serves Zone 2 and Zone 3
  cap_f3: x32 + x33 - 100 y3 <= 0

Bounds
  0 <= x33
  0 <= x13
  0 <= y3 <= 1

Generals
  y3

End

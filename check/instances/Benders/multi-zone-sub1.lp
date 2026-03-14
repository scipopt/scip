\ Ring-Structured Capacity Expansion
\ Master Variables: y1, y2, y3 (Facility Investments)

Minimize
  obj: 5 x11 + 8 x21  \ Zone 1 served by Facility 1 and 2

Subject To
  \ --- DEMAND CONSTRAINTS ---
  \ Demand 1 can be met by Facility 1 OR Facility 2
  demand_z1: x11 + x21 >= 50

  \ --- CAPACITY CONSTRAINTS (Linking Master y to Subproblem x) ---
  \ Facility 1 serves Zone 1 and Zone 3
  cap_f1: x11 + x13 - 100 y1 <= 0

Bounds
  0 <= x11
  0 <= x21
  0 <= y1 <= 1

Generals
  y1

End

\ Ring-Structured Capacity Expansion
\ Master Variables: y1, y2, y3 (Facility Investments)

Minimize
  obj: 100 y1 + 120 y2 + 110 y3

Subject To
  \ Master Budget Constraint (Still only allows 2 facilities)
  budget: y1 + y2 + y3 <= 2

Bounds
  0 <= y1 <= 1
  0 <= y2 <= 1
  0 <= y3 <= 1

Generals
  y1 y2 y3

End

\ Multi-Region Facility Location Problem
\ Master Variables: y_north, y_central, y_south

Minimize
  obj: 1000 y_north + 1200 y_central + 1100 y_south

Subject To
  \ Optional: Policy constraint (Must open at least 2 facilities)
  min_open: y_north + y_central + y_south >= 2

Bounds
  y_north <= 1
  y_central <= 1
  y_south <= 1

Binaries
  y_north y_central y_south

End

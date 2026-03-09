\ Multi-Region Facility Location Problem
\ Master Variables: y_north, y_central, y_south

Minimize
  obj: 4 x_n_c1 + 5 x_n_c2 + 3 x_c_c1 + 6 x_c_c2 + 5 x_s_c1 + 2 x_s_c2

Subject To
  \ Demand constraints for Customers
  demand_c1: x_n_c1 + x_c_c1 + x_s_c1 >= 80
  demand_c2: x_n_c2 + x_c_c2 + x_s_c2 >= 70

  \ Capacity constraints (Linking Master and Subproblem)
  \ Facility North capacity = 100
  cap_north:   x_n_c1 + x_n_c2 - 100 y_north   <= 0
  \ Facility Central capacity = 150
  cap_central: x_c_c1 + x_c_c2 - 150 y_central <= 0
  \ Facility South capacity = 90
  cap_south:   x_s_c1 + x_s_c2 - 90 y_south   <= 0

Bounds
  0 <= x_n_c1
  0 <= x_n_c2
  0 <= x_c_c1
  0 <= x_c_c2
  0 <= x_s_c1
  0 <= x_s_c2
  y_north <= 1
  y_central <= 1
  y_south <= 1

Binaries
  y_north y_central y_south

End

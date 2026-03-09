* SCIP STATISTICS
*   Problem name     : sub.lp
*   Variables        : 9 (3 binary, 0 integer, 0 implicit integer, 6 continuous)
*   Constraints      : 5
NAME          sub.lp
OBJSENSE
  MIN
ROWS
 N  Obj
 G  demand_c1
 G  demand_c2
 L  cap_north
 L  cap_central
 L  cap_south
COLUMNS
    INTSTART     'MARKER'                               'INTORG'
    y_north      Obj                                 0  cap_north                        -100
    y_central    Obj                                 0  cap_central                      -150
    y_south      cap_south                         -90  Obj                                 0
    INTEND       'MARKER'                               'INTEND'
    x_n_c1       demand_c1                           1  cap_north                           1
    x_n_c1       Obj                                 4
    x_n_c2       Obj                                 5  cap_north                           1
    x_n_c2       demand_c2                           1
    x_c_c1       Obj                                 3  cap_central                         1
    x_c_c1       demand_c1                           1
    x_c_c2       Obj                                 6  cap_central                         1
    x_c_c2       demand_c2                           1
    x_s_c1       demand_c1                           1  cap_south                           1
    x_s_c1       Obj                                 5
    x_s_c2       demand_c2                           1  Obj                                 2
    x_s_c2       cap_south                           1
RHS
    RHS          demand_c1                          80  demand_c2                          70
    RHS          cap_north                           0  cap_central                         0
    RHS          cap_south                           0
BOUNDS
 BV Bound        y_north
 BV Bound        y_central
 BV Bound        y_south
 PL Bound        x_c_c2
 PL Bound        x_s_c1
 PL Bound        x_s_c2
 PL Bound        x_n_c2
 PL Bound        x_n_c1
 PL Bound        x_c_c1
ENDATA

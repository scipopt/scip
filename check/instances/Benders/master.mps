* SCIP STATISTICS
*   Problem name     : master.lp
*   Variables        : 3 (3 binary, 0 integer, 0 implicit integer, 0 continuous)
*   Constraints      : 1
NAME          master.lp
OBJSENSE
  MIN
ROWS
 N  Obj
 G  min_open
COLUMNS
    INTSTART   'MARKER'                             'INTORG'
    y_north    Obj                            1000  min_open                          1
    y_central  Obj                            1200  min_open                          1
    y_south    Obj                            1100  min_open                          1
    INTEND     'MARKER'                             'INTEND'
RHS
    RHS        min_open                          2
BOUNDS
 BV Bound      y_north
 BV Bound      y_central
 BV Bound      y_south
ENDATA

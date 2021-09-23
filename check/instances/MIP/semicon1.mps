*  MIP written by GAMS Convert at 09/13/10 21:35:46
*  
*  Equation counts
*      Total        E        G        L        N        X        C
*          4        2        1        1        0        0        0
*  
*  Variable counts
*                   x        b        i      s1s      s2s       sc       si
*      Total     cont   binary  integer     sos1     sos2    scont     sint
*          5        4        0        0        0        0        1        0
*  FX      0        0        0        0        0        0        0        0
*  
*  Nonzero counts
*      Total    const       NL      DLL
*          9        9        0        0
*
NAME          Semicon1
*
* original model was minimizing
*
ROWS
 N  obj     
 E  c1      
 L  c2      
 G  c3      
 E  c4      
COLUMNS
    x1        obj                  1
    x1        c1                   1
    sc2       c4                   1
    x3        c1                  -1
    x3        c2                  -1
    x4        c1                  -1
    x4        c3                   1
    x5        c2                   1
    x5        c3                   1
    x5        c4                   1
RHS
    rhs       c2                 8.9
    rhs       c3                 8.9
    rhs       c4                  10
BOUNDS
 FR bnd       x1      
 LO bnd       sc2                2.8
 SC bnd       sc2                 10
ENDATA

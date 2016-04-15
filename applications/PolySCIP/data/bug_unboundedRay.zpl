set I := {1..2};
param c1[I] := <1> 1, <2> 0.8;
param c2[I] := <1> 0, <2> -1;
var x[I] integer >= 0;

minimize Obj1: sum <i> in I: c1[i]*x[i];
Obj2: sum <i> in I: c2[i]*x[2];
subto Cons1: x[1] + x[2] >= 2;
subto Cons2: -1*x[1] + x[2] <= 0;
subto Cons3: -0.5*x[1] + x[2] <= 1.5;
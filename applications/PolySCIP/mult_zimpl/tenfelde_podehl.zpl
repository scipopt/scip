set I := {1..4};

param obj1[I*I] := |1, 2, 3, 4|
                 |1|3, 6, 4, 5|
                 |2|2, 3, 5, 4|
                 |3|3, 5, 4, 2|
                 |4|4, 5, 3, 6|;

param obj2[I*I] := |1, 2, 3, 4|
                 |1|2, 3, 5, 4|
                 |2|5, 3, 4, 3|
                 |3|5, 2, 6, 4|
                 |4|4, 5, 2, 5|;

param obj3[I*I] := |1, 2, 3, 4|
                 |1|4, 2, 4, 2|
                 |2|4, 2, 4, 6|
                 |3|4, 2, 6, 3|
                 |4|2, 4, 5, 3|;

var x[I*I] binary;

minimize Obj1: sum <i,j> in I*I: obj1[i,j]*x[i,j];
Obj2: sum <i,j> in I*I: obj2[i,j]*x[i,j];
Obj3: sum <i,j> in I*I: obj3[i,j]*x[i,j];

subto rows: forall <i> in I do
      sum <j> in I: x[i,j] == 1;

subto col: forall <i> in I do
      sum <j> in I: x[j,i] == 1;
      
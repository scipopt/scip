set I := { 1 .. 10 };

var x[I];

maximize xyz: sum <i> in I : x[i];

subto c1: x[1] x[2] == x[3] x[4];
subto c2: x[5] x[6] == x[7] x[8];
subto c3: sum <i> in I : x[i] <= 20;

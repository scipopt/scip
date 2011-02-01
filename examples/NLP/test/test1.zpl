set I := { 1 .. 10 };

var x[I];

maximize xyz: sum <i> in I : i*x[i];

subto c1: x[1] x[2] == x[3] x[4];
subto c2: sum <i> in I : i*x[i] <= 20;

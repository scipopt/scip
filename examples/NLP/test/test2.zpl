set I := { 1 .. 4 };

var x[I];

maximize xyz: sum <i> in I : x[i];

subto c1: 5 * x[1] x[1] == x[3] x[4];

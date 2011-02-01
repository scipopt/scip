set I := { 1 .. 4 };

var x[I] <= 100;
var f <= 1000;

maximize xyz: f;

subto c1: f == x[1] x[2] + x[3] x[4];
subto c2: x[1] x[3] <= 17;
subto c3: x[2] x[4] <= 22;

var x >= 0;
var y >= 0;

var z binary;
var z2 binary;

minimize obj: -2*x - y;
subto c1: x <= 1;
subto c2: y <= 1;

#subto i: vif z == 1  then x == 0 else y == 0 end, indicator;
subto i1: vif z == 0 then x <= 0 end, indicator;
subto i2: vif z == 1 then y <= 0 end, indicator;

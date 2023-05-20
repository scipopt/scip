# ampl -oglogic3 -P logic3.mod

var x binary;
var y binary;
var z binary;

minimize obj: z - 0.5*x;

subject to
  e1: (x and y) or (x == y) or (not z);

# x = 1, z = 0 should be optimal

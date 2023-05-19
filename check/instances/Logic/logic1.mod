# ampl -oglogic1 -P logic1.mod

var x binary;
var y binary;

minimize obj: y ;

subject to
  e1: x or y;
  e2: not x;
  e3: not x == y;

# x = 0, y = 1 is the only feasible solution here

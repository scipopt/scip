# ampl -oglogic2 -P logic2.mod

var x binary;
var y binary;

minimize obj: x + y;

subject to
  e1: x and y;
  e2: (not x) or y;
  e3: x == y;

# x = 1, y = 1 is the only feasible solution here

# create nl file via  ampl -P -ognlobj nlobj.mod

var x >= 0 := 1;
var y >= 0 := 2;

minimize obj: x + 2*y*y + 5 ;

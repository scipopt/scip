# create nl file via  ampl -P -oglp1 lp1.mod

var x >= 0 ;
var y >= 0 ;
var z >= 0 ;

maximize obj: x + 2*y + z ;

subject to
  e1: x + y <= 10;
  e2: -y - z >= -50;

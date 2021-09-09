# create nl file via  ampl -P -ogcommonexpr1 commonexpr1.mod

var x >= 0 ;
var xsqr = x^2 ;

minimize obj: x ;

subject to
  e1: xsqr - exp(xsqr) + x >= 1 ;

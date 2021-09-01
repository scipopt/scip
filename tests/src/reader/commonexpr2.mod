# create nl file via  ampl -P -ogcommonexpr2 commonexpr2.mod

var x >= 0 ;
var xsqr = x^2 ;
var xcube = x^3 ;
var xsqrt = sqrt(x) ;

minimize obj: x ;

subject to
  e1: xsqr >= 1;
  e2: xsqr + xsqrt >= 1;
  e3: xsqr + xsqrt + xcube >= 1;

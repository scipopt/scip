# create nl file via  ampl -P -ogsuffix1 suffix1.mod

var x >= 0 ;
var y >= 0 ;
var z >= 0 , <= 100 ;

maximize obj: x + y + z ;

subject to
  e1: x + y <= 10;
  e2: x + y <= z;   # constraint that doesn't really matter in an optimal solution (z=100)

suffix initial   integer IN;
suffix separate  integer IN;
suffix enforce   integer IN;
suffix check     integer IN;
suffix propagate integer IN;
suffix dynamic   integer IN;
suffix removable integer IN;
suffix sosno     integer IN;

let x.initial    := 2;
let x.removable  := 1;
let e2.initial   := 2;
let e2.separate  := 2;
let e2.enforce   := 2;
let e2.check     := 2;
let e2.propagate := 2;
let e2.dynamic   := 1;
let e2.removable := 1;

# (x,y) are in a SOS1; we omit .ref here, as should be optional
let x.sosno := 1;
let y.sosno := 1;

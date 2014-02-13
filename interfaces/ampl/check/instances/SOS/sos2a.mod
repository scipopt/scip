suffix sosno integer IN;
suffix ref integer IN; 

var s1 >= 0.8;
var s2 >= 0;
var s3 >= 0;
var x4 >= 0;
var x5 >= 0;
var x7;
var x8;

minimize obj:    x4 + x5;

subject to

e1:  s1 + s2 + s3 = 1;

e2:  s1 + 2*s2 + 3*s3 = x7;

e3:  s1 + 2*s2 + 3*s3 = x8;

e6:  x4 - x8 >= -1.3;

e7:  x5 + x8 >= 1.3;

let s1.sosno := -1;
let s2.sosno := -1;
let s3.sosno := -1;

let s1.ref := 1;
let s2.ref := 2;
let s3.ref := 3;

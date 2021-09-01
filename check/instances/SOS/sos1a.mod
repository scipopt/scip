suffix sosno integer IN;
suffix ref integer IN; 

var s1 >= 0, <= 0.8;
var s2 >= 0, <= 0.6;
var s3 >= 0, <= 0.6;

maximize obj:    0.9*s1 + s2 + 1.1*s3;

subject to

e2:    s1 + s2 + s3 <= 1;

let s1.sosno := 1;
let s2.sosno := 1;
let s3.sosno := 1;

let s1.ref := 1;
let s2.ref := 2;
let s3.ref := 3;

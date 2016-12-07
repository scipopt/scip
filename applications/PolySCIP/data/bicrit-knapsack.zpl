# instances from http://l1.lamsade.dauphine.fr/~hugot/instances.html
param n := read "Type_A/100_items/2KP100-TA-0.dat" as "2n" skip 2 use 1;
set I := {1..n};
param weight[I] := read "Type_A/100_items/2KP100-TA-0.dat" as "2n" skip 4 use 100;
param cost1[I] := read "Type_A/100_items/2KP100-TA-0.dat" as "3n" skip 4 use 100;
param cost2[I] := read "Type_A/100_items/2KP100-TA-0.dat" as "4n" skip 4 use 100;
param max_weight := read "Type_A/100_items/2KP100-TA-0.dat" as "2n" skip 106 use 1;
var x[I] binary;

maximize Obj: sum <i> in I: cost1[i]*x[i];
Obj2: sum <i> in I: cost2[i]*x[i];

subto Weight: sum <i> in I: weight[i]*x[i] <= max_weight;
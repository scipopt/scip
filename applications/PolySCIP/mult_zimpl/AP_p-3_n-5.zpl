# generates tri-criteria assignment problem
param data_file := "../data/AP_p-3_n-5.dat";
param no_objs := read data_file as "1n" use 1;
param no_vars := read data_file as "1n" use 1 skip 1;

set I := {1..no_vars};

set T := {1..no_objs*no_vars*no_vars};
param coeffs[T] := read data_file as "n+" match "[0-9]+" skip 2;

param Obj1[<i,j> in I*I] := coeffs[(i-1)*no_vars + j];

param offset := no_vars*no_vars;
param Obj2[<i,j> in I*I] := coeffs[(i-1)*no_vars + j + offset];

param Obj3[<i,j> in I*I] := coeffs[(i-1)*no_vars + j + 2*offset];

var x[I*I] binary;

minimize Obj1: sum <i,j> in I*I: Obj1[i,j]*x[i,j];
Obj2: sum <i,j> in I*I: Obj2[i,j]*x[i,j];
Obj3: sum <i,j> in I*I: Obj3[i,j]*x[i,j];

subto row: forall <i> in I do
      sum <j> in I: x[i,j] == 1;	

subto col: forall <i> in I do
      sum <j> in I: x[j,i] == 1;
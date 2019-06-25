param objsen := 1;
set NZO := {
 <"OBJECTIV","x_v1", 1>
,
 <"E4_1","x_v1", 1>
,
 <"E2_1","x_v1", 1>
,
 <"E1_1","x_v1", 1>
,
 <"OBJECTIV","x_v2", 1>
,
 <"E5_1","x_v2", 1>
,
 <"E3_1","x_v2", 1>
,
 <"E1_1","x_v2", 1>
,
 <"OBJECTIV","x_v3", 1>
,
 <"E6_1","x_v3", 1>
,
 <"E3_1","x_v3", 1>
,
 <"E2_1","x_v3", 1>
,
 <"OBJECTIV","x_v4", 1>
,
 <"E6_1","x_v4", 1>
,
 <"E5_1","x_v4", 1>
,
 <"E4_1","x_v4", 1>
};
set IVAR := {
<"x_v1">
,<"x_v2">
,<"x_v3">
,<"x_v4">
};
set CVAR := {
};
set V := IVAR + CVAR;
param lower[V] := 
<"x_v1">0
,<"x_v2">0
,<"x_v3">0
,<"x_v4">0
;
param upper[V] := 
<"x_v1">1
,<"x_v2">1
,<"x_v3">1
,<"x_v4">1
;
set RANGEROW := {
<" dummy ">
};
set ROW := { 
<"E3_1","G">
,<"E6_1","G">
,<"E1_1","G">
,<"E4_1","G">
,<"OBJECTIV","N">
,<"E2_1","G">
,<"E5_1","G">
};
set R := proj(ROW, <1>);
param rhs[R] := 
<"E3_1">1
,<"E6_1">1
,<"E1_1">1
,<"E4_1">1
,<"OBJECTIV">0
,<"E2_1">1
,<"E5_1">1
;
param range[RANGEROW] := 
<" dummy ">0
;
var x[<v> in IVAR] integer >= if lower[v] > -1e20 then lower[v] else -infinity end <= if upper[v] < 1e20 then upper[v] else infinity end;
minimize obj: sum <r,"N"> in ROW : (
sum <r,c,a> in NZO with <c> in IVAR: objsen * a * x[c] 
);
subto ce: forall <r,"E"> in ROW do
if <r> in RANGEROW and sgn(range[r]) > 0 then 
rhs[r] <= 
sum <r,c,a> in NZO with <c> in IVAR : a * x[c]
<= rhs[r] + abs(range[r])
else if <r> in RANGEROW and sgn(range[r]) < 0 then 
rhs[r] - abs(range[r]) <= 
sum <r,c,a> in NZO with <c> in IVAR : a * x[c]
<= rhs[r]
 else 
sum <r,c,a> in NZO with <c> in IVAR : a * x[c]
== rhs[r]
 end end;
subto cl: forall <r,"L"> in ROW do 
if <r> in RANGEROW then 
rhs[r] - abs(range[r]) <= 
sum <r,c,a> in NZO with <c> in IVAR : a * x[c] 
<= rhs[r]
 else 
sum <r,c,a> in NZO with <c> in IVAR : a * x[c] 
<= rhs[r]
 end;
subto cg: forall <r,"G"> in ROW do 
if <r> in RANGEROW then 
rhs[r] + abs(range[r]) >= 
sum <r,c,a> in NZO with <c> in IVAR : a * x[c] 
>= rhs[r]
 else 
sum <r,c,a> in NZO with <c> in IVAR : a * x[c] 
>= rhs[r]
 end;

param objsen := 1;
set NZO := {
 <"KOSTEN","STM1", 2700>, <"ANZ1","STM1", 1>
,
 <"STD1","STM1", 150>, <"UEB1","STM1", -20>
,
 <"ANZ2","STM1", 0.9>
,
 <"KOSTEN","ANM1", 1500>, <"STD1","ANM1", -100>
,
 <"ANZ2","ANM1", 1>
,
 <"KOSTEN","UE1", 30>, <"STD1","UE1", 1>
,
 <"UEB1","UE1", 1>
,
 <"KOSTEN","STM2", 2700>, <"ANZ2","STM2", -1>
,
 <"STD2","STM2", 150>, <"UEB2","STM2", -20>
,
 <"ANZ3","STM2", 0.9>
,
 <"KOSTEN","ANM2", 1500>, <"STD2","ANM2", -100>
,
 <"ANZ3","ANM2", 1>
,
 <"KOSTEN","UE2", 30>, <"STD2","UE2", 1>
,
 <"UEB2","UE2", 1>
,
 <"KOSTEN","STM3", 2700>, <"ANZ3","STM3", -1>
,
 <"STD3","STM3", 150>, <"UEB3","STM3", -20>
,
 <"ANZ4","STM3", 0.9>
,
 <"KOSTEN","ANM3", 1500>, <"STD3","ANM3", -100>
,
 <"ANZ4","ANM3", 1>
,
 <"KOSTEN","UE3", 30>, <"STD3","UE3", 1>
,
 <"UEB3","UE3", 1>
,
 <"KOSTEN","STM4", 2700>, <"ANZ4","STM4", -1>
,
 <"STD4","STM4", 150>, <"UEB4","STM4", -20>
,
 <"ANZ5","STM4", 0.9>
,
 <"KOSTEN","ANM4", 1500>, <"STD4","ANM4", -100>
,
 <"ANZ5","ANM4", 1>
,
 <"KOSTEN","UE4", 30>, <"STD4","UE4", 1>
,
 <"UEB4","UE4", 1>
,
 <"KOSTEN","STM5", 2700>, <"ANZ5","STM5", -1>
,
 <"STD5","STM5", 150>, <"UEB5","STM5", -20>
,
 <"ANZ6","STM5", 0.9>
,
 <"KOSTEN","ANM5", 1500>, <"STD5","ANM5", -100>
,
 <"ANZ6","ANM5", 1>
,
 <"KOSTEN","UE5", 30>, <"STD5","UE5", 1>
,
 <"UEB5","UE5", 1>
,
 <"KOSTEN","STM6", 2700>, <"ANZ6","STM6", -1>
,
 <"STD6","STM6", 150>, <"UEB6","STM6", -20>
,
 <"KOSTEN","ANM6", 1500>, <"STD6","ANM6", -100>
,
 <"KOSTEN","UE6", 30>, <"STD6","UE6", 1>
,
 <"UEB6","UE6", 1>
};
set IVAR := {
<"ANM5">
,<"ANM6">
,<"STM2">
,<"STM3">
,<"ANM1">
,<"STM4">
,<"ANM2">
,<"STM5">
,<"ANM3">
,<"STM6">
,<"ANM4">
};
set CVAR := {
<"UE6">
,<"STM1">
,<"UE1">
,<"UE2">
,<"UE3">
,<"UE4">
,<"UE5">
};
set V := IVAR + CVAR;
param lower[V] := 
<"UE6">0
,<"ANM5">0
,<"ANM6">0
,<"STM1">0
,<"STM2">57
,<"UE1">0
,<"STM3">57
,<"UE2">0
,<"ANM1">0
,<"STM4">57
,<"UE3">0
,<"ANM2">0
,<"STM5">57
,<"UE4">0
,<"ANM3">0
,<"STM6">57
,<"UE5">0
,<"ANM4">0
;
param upper[V] := 
<"UE6">1e21
,<"ANM5">18
,<"ANM6">18
,<"STM1">1e21
,<"STM2">75
,<"UE1">1e21
,<"STM3">75
,<"UE2">1e21
,<"ANM1">18
,<"STM4">75
,<"UE3">1e21
,<"ANM2">18
,<"STM5">75
,<"UE4">1e21
,<"ANM3">18
,<"STM6">75
,<"UE5">1e21
,<"ANM4">18
;
set RANGEROW := {
<" dummy ">
};
set ROW := { 
<"UEB5","L">
,<"ANZ5","E">
,<"STD3","G">
,<"UEB6","L">
,<"ANZ6","E">
,<"STD4","G">
,<"STD5","G">
,<"STD6","G">
,<"KOSTEN","N">
,<"UEB1","L">
,<"ANZ1","E">
,<"UEB2","L">
,<"ANZ2","E">
,<"UEB3","L">
,<"ANZ3","E">
,<"STD1","G">
,<"UEB4","L">
,<"ANZ4","E">
,<"STD2","G">
};
set R := proj(ROW, <1>);
param rhs[R] := 
<"UEB5">0
,<"ANZ5">0
,<"STD3">8000
,<"UEB6">0
,<"ANZ6">0
,<"STD4">10000
,<"STD5">9000
,<"STD6">12000
,<"KOSTEN">0
,<"UEB1">0
,<"ANZ1">60
,<"UEB2">0
,<"ANZ2">0
,<"UEB3">0
,<"ANZ3">0
,<"STD1">8000
,<"UEB4">0
,<"ANZ4">0
,<"STD2">9000
;
param range[RANGEROW] := 
<" dummy ">0
;
var x[<v> in IVAR] integer >= if lower[v] > -1e20 then lower[v] else -infinity end <= if upper[v] < 1e20 then upper[v] else infinity end;
var y[<v> in CVAR] real >= if lower[v] > -1e20 then lower[v] else -infinity end <= if upper[v] < 1e20 then upper[v] else infinity end;
minimize obj: sum <r,"N"> in ROW : (
sum <r,c,a> in NZO with <c> in IVAR: objsen * a * x[c] 
+ sum <r,c,a> in NZO with <c> in CVAR: objsen * a * y[c]
);
subto ce: forall <r,"E"> in ROW do
if <r> in RANGEROW and sgn(range[r]) > 0 then 
rhs[r] <= 
sum <r,c,a> in NZO with <c> in IVAR : a * x[c]
+ sum <r,c,a> in NZO with <c> in CVAR : a * y[c]
<= rhs[r] + abs(range[r])
else if <r> in RANGEROW and sgn(range[r]) < 0 then 
rhs[r] - abs(range[r]) <= 
sum <r,c,a> in NZO with <c> in IVAR : a * x[c]
+ sum <r,c,a> in NZO with <c> in CVAR : a * y[c]
<= rhs[r]
 else 
sum <r,c,a> in NZO with <c> in IVAR : a * x[c]
+ sum <r,c,a> in NZO with <c> in CVAR : a * y[c]
== rhs[r]
 end end;
subto cl: forall <r,"L"> in ROW do 
if <r> in RANGEROW then 
rhs[r] - abs(range[r]) <= 
sum <r,c,a> in NZO with <c> in IVAR : a * x[c] 
+ sum <r,c,a> in NZO with <c> in CVAR : a * y[c] 
<= rhs[r]
 else 
sum <r,c,a> in NZO with <c> in IVAR : a * x[c] 
+ sum <r,c,a> in NZO with <c> in CVAR : a * y[c] 
<= rhs[r]
 end;
subto cg: forall <r,"G"> in ROW do 
if <r> in RANGEROW then 
rhs[r] + abs(range[r]) >= 
sum <r,c,a> in NZO with <c> in IVAR : a * x[c] 
+ sum <r,c,a> in NZO with <c> in CVAR : a * y[c] 
>= rhs[r]
 else 
sum <r,c,a> in NZO with <c> in IVAR : a * x[c] 
+ sum <r,c,a> in NZO with <c> in CVAR : a * y[c] 
>= rhs[r]
 end;

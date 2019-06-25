param objsen := 1;
set NZO := {
 <"SOS0","A0", 1>, <"SOSA","A0", 1>
,
 <"OBJECT","A1", 1>, <"SOS1","A1", 1>
,
 <"BILANCIO","A1", 202>, <"SOSA","A1", 1>
,
 <"OBJECT","A2", 2>, <"SOS2","A2", 1>
,
 <"BILANCIO","A2", 404>, <"SOSA","A2", 1>
,
 <"OBJECT","A3", 3>, <"SOS3","A3", 1>
,
 <"BILANCIO","A3", 606>, <"SOSA","A3", 1>
,
 <"OBJECT","A4", 4>, <"SOS4","A4", 1>
,
 <"BILANCIO","A4", 808>, <"SOSA","A4", 1>
,
 <"OBJECT","A5", 5>, <"SOS5","A5", 1>
,
 <"BILANCIO","A5", 1010>, <"SOSA","A5", 1>
,
 <"OBJECT","A6", 6>, <"SOS6","A6", 1>
,
 <"BILANCIO","A6", 1212>, <"SOSA","A6", 1>
,
 <"OBJECT","A7", 7>, <"SOS7","A7", 1>
,
 <"BILANCIO","A7", 1414>, <"SOSA","A7", 1>
,
 <"OBJECT","A8", 8>, <"SOS8","A8", 1>
,
 <"BILANCIO","A8", 1616>, <"SOSA","A8", 1>
,
 <"OBJECT","A9", 9>, <"SOS9","A9", 1>
,
 <"BILANCIO","A9", 1818>
,
 <"SOS0","B0", 1>, <"SOSB","B0", 1>
,
 <"SOS1","B1", 1>, <"BILANCIO","B1", -79>
,
 <"SOSB","B1", 1>
,
 <"SOS2","B2", 1>, <"BILANCIO","B2", -158>
,
 <"SOSB","B2", 1>
,
 <"SOS3","B3", 1>, <"BILANCIO","B3", -237>
,
 <"SOSB","B3", 1>
,
 <"SOS4","B4", 1>, <"BILANCIO","B4", -316>
,
 <"SOSB","B4", 1>
,
 <"SOS5","B5", 1>, <"BILANCIO","B5", -395>
,
 <"SOSB","B5", 1>
,
 <"SOS6","B6", 1>, <"BILANCIO","B6", -474>
,
 <"SOSB","B6", 1>
,
 <"SOS7","B7", 1>, <"BILANCIO","B7", -553>
,
 <"SOSB","B7", 1>
,
 <"SOS8","B8", 1>, <"BILANCIO","B8", -632>
,
 <"SOSB","B8", 1>
,
 <"SOS9","B9", 1>, <"BILANCIO","B9", -711>
,
 <"SOSB","B9", 1>
,
 <"SOS0","C0", 1>, <"SOSC","C0", 1>
,
 <"SOS1","C1", 1>, <"BILANCIO","C1", 100023>
,
 <"SOSC","C1", 1>
,
 <"SOS2","C2", 1>, <"BILANCIO","C2", 200046>
,
 <"SOSC","C2", 1>
,
 <"SOS3","C3", 1>, <"BILANCIO","C3", 300069>
,
 <"SOSC","C3", 1>
,
 <"SOS4","C4", 1>, <"BILANCIO","C4", 400092>
,
 <"SOSC","C4", 1>
,
 <"SOS5","C5", 1>, <"BILANCIO","C5", 500115>
,
 <"SOSC","C5", 1>
,
 <"SOS6","C6", 1>, <"BILANCIO","C6", 600138>
,
 <"SOSC","C6", 1>
,
 <"SOS7","C7", 1>, <"BILANCIO","C7", 700161>
,
 <"SOSC","C7", 1>
,
 <"SOS8","C8", 1>, <"BILANCIO","C8", 800184>
,
 <"SOSC","C8", 1>
,
 <"SOS9","C9", 1>, <"BILANCIO","C9", 900207>
,
 <"SOSC","C9", 1>
,
 <"SOS0","D0", 1>, <"SOSD","D0", 1>
,
 <"SOS1","D1", 1>, <"BILANCIO","D1", -89810>
,
 <"SOSD","D1", 1>
,
 <"SOS2","D2", 1>, <"BILANCIO","D2", -179620>
,
 <"SOSD","D2", 1>
,
 <"SOS3","D3", 1>, <"BILANCIO","D3", -269430>
,
 <"SOSD","D3", 1>
,
 <"SOS4","D4", 1>, <"BILANCIO","D4", -359240>
,
 <"SOSD","D4", 1>
,
 <"SOS5","D5", 1>, <"BILANCIO","D5", -449050>
,
 <"SOSD","D5", 1>
,
 <"SOS6","D6", 1>, <"BILANCIO","D6", -538860>
,
 <"SOSD","D6", 1>
,
 <"SOS7","D7", 1>, <"BILANCIO","D7", -628670>
,
 <"SOSD","D7", 1>
,
 <"SOS8","D8", 1>, <"BILANCIO","D8", -718480>
,
 <"SOSD","D8", 1>
,
 <"SOS9","D9", 1>, <"BILANCIO","D9", -808290>
,
 <"SOSD","D9", 1>
,
 <"SOS0","E_0", 1>, <"SOSE","E_0", 1>
,
 <"SOS1","E_1", 1>, <"BILANCIO","E_1", -9980>
,
 <"SOSE","E_1", 1>
,
 <"SOS2","E_2", 1>, <"BILANCIO","E_2", -19960>
,
 <"SOSE","E_2", 1>
,
 <"SOS3","E_3", 1>, <"BILANCIO","E_3", -29940>
,
 <"SOSE","E_3", 1>
,
 <"SOS4","E_4", 1>, <"BILANCIO","E_4", -39920>
,
 <"SOSE","E_4", 1>
,
 <"SOS5","E_5", 1>, <"BILANCIO","E_5", -49900>
,
 <"SOSE","E_5", 1>
,
 <"SOS6","E_6", 1>, <"BILANCIO","E_6", -59880>
,
 <"SOSE","E_6", 1>
,
 <"SOS7","E_7", 1>, <"BILANCIO","E_7", -69860>
,
 <"SOSE","E_7", 1>
,
 <"SOS8","E_8", 1>, <"BILANCIO","E_8", -79840>
,
 <"SOSE","E_8", 1>
,
 <"SOS9","E_9", 1>, <"BILANCIO","E_9", -89820>
,
 <"SOSE","E_9", 1>
,
 <"SOS0","F0", 1>, <"SOSF","F0", 1>
,
 <"SOS1","F1", 1>, <"BILANCIO","F1", 1000>
,
 <"SOSF","F1", 1>
,
 <"SOS2","F2", 1>, <"BILANCIO","F2", 2000>
,
 <"SOSF","F2", 1>
,
 <"SOS3","F3", 1>, <"BILANCIO","F3", 3000>
,
 <"SOSF","F3", 1>
,
 <"SOS4","F4", 1>, <"BILANCIO","F4", 4000>
,
 <"SOSF","F4", 1>
,
 <"SOS5","F5", 1>, <"BILANCIO","F5", 5000>
,
 <"SOSF","F5", 1>
,
 <"SOS6","F6", 1>, <"BILANCIO","F6", 6000>
,
 <"SOSF","F6", 1>
,
 <"SOS7","F7", 1>, <"BILANCIO","F7", 7000>
,
 <"SOSF","F7", 1>
,
 <"SOS8","F8", 1>, <"BILANCIO","F8", 8000>
,
 <"SOSF","F8", 1>
,
 <"SOS9","F9", 1>, <"BILANCIO","F9", 9000>
,
 <"SOSF","F9", 1>
,
 <"SOS0","G0", 1>, <"SOSG","G0", 1>
,
 <"SOS1","G1", 1>, <"BILANCIO","G1", 100>
,
 <"SOSG","G1", 1>
,
 <"SOS2","G2", 1>, <"BILANCIO","G2", 200>
,
 <"SOSG","G2", 1>
,
 <"SOS3","G3", 1>, <"BILANCIO","G3", 300>
,
 <"SOSG","G3", 1>
,
 <"SOS4","G4", 1>, <"BILANCIO","G4", 400>
,
 <"SOSG","G4", 1>
,
 <"SOS5","G5", 1>, <"BILANCIO","G5", 500>
,
 <"SOSG","G5", 1>
,
 <"SOS6","G6", 1>, <"BILANCIO","G6", 600>
,
 <"SOSG","G6", 1>
,
 <"SOS7","G7", 1>, <"BILANCIO","G7", 700>
,
 <"SOSG","G7", 1>
,
 <"SOS8","G8", 1>, <"BILANCIO","G8", 800>
,
 <"SOSG","G8", 1>
,
 <"SOS9","G9", 1>, <"BILANCIO","G9", 900>
,
 <"SOSG","G9", 1>
,
 <"SOS0","H0", 1>, <"SOSH","H0", 1>
,
 <"SOS1","H1", 1>, <"BILANCIO","H1", 10000>
,
 <"SOSH","H1", 1>
,
 <"SOS2","H2", 1>, <"BILANCIO","H2", 20000>
,
 <"SOSH","H2", 1>
,
 <"SOS3","H3", 1>, <"BILANCIO","H3", 30000>
,
 <"SOSH","H3", 1>
,
 <"SOS4","H4", 1>, <"BILANCIO","H4", 40000>
,
 <"SOSH","H4", 1>
,
 <"SOS5","H5", 1>, <"BILANCIO","H5", 50000>
,
 <"SOSH","H5", 1>
,
 <"SOS6","H6", 1>, <"BILANCIO","H6", 60000>
,
 <"SOSH","H6", 1>
,
 <"SOS7","H7", 1>, <"BILANCIO","H7", 70000>
,
 <"SOSH","H7", 1>
,
 <"SOS8","H8", 1>, <"BILANCIO","H8", 80000>
,
 <"SOSH","H8", 1>
,
 <"SOS9","H9", 1>, <"BILANCIO","H9", 90000>
,
 <"SOSH","H9", 1>
,
 <"SOS0","I0", 1>, <"SOSI","I0", 1>
,
 <"SOS1","I1", 1>, <"BILANCIO","I1", 100>
,
 <"SOSI","I1", 1>
,
 <"SOS2","I2", 1>, <"BILANCIO","I2", 200>
,
 <"SOSI","I2", 1>
,
 <"SOS3","I3", 1>, <"BILANCIO","I3", 300>
,
 <"SOSI","I3", 1>
,
 <"SOS4","I4", 1>, <"BILANCIO","I4", 400>
,
 <"SOSI","I4", 1>
,
 <"SOS5","I5", 1>, <"BILANCIO","I5", 500>
,
 <"SOSI","I5", 1>
,
 <"SOS6","I6", 1>, <"BILANCIO","I6", 600>
,
 <"SOSI","I6", 1>
,
 <"SOS7","I7", 1>, <"BILANCIO","I7", 700>
,
 <"SOSI","I7", 1>
,
 <"SOS8","I8", 1>, <"BILANCIO","I8", 800>
,
 <"SOSI","I8", 1>
,
 <"SOS9","I9", 1>, <"BILANCIO","I9", 900>
,
 <"SOSI","I9", 1>
,
 <"SOS0","L0", 1>, <"SOSL","L0", 1>
,
 <"SOS1","L1", 1>, <"BILANCIO","L1", -1>
,
 <"SOSL","L1", 1>
,
 <"SOS2","L2", 1>, <"BILANCIO","L2", -2>
,
 <"SOSL","L2", 1>
,
 <"SOS3","L3", 1>, <"BILANCIO","L3", -3>
,
 <"SOSL","L3", 1>
,
 <"SOS4","L4", 1>, <"BILANCIO","L4", -4>
,
 <"SOSL","L4", 1>
,
 <"SOS5","L5", 1>, <"BILANCIO","L5", -5>
,
 <"SOSL","L5", 1>
,
 <"SOS6","L6", 1>, <"BILANCIO","L6", -6>
,
 <"SOSL","L6", 1>
,
 <"SOS7","L7", 1>, <"BILANCIO","L7", -7>
,
 <"SOSL","L7", 1>
,
 <"SOS8","L8", 1>, <"BILANCIO","L8", -8>
,
 <"SOSL","L8", 1>
,
 <"SOS9","L9", 1>, <"BILANCIO","L9", -9>
,
 <"SOSL","L9", 1>
};
set IVAR := {
<"B0">
,<"B1">
,<"B2">
,<"B3">
,<"B4">
,<"B5">
,<"B6">
,<"D0">
,<"B7">
,<"B8">
,<"D1">
,<"B9">
,<"D2">
,<"D3">
,<"D4">
,<"D5">
,<"D6">
,<"F0">
,<"D7">
,<"F1">
,<"D8">
,<"F2">
,<"D9">
,<"F3">
,<"F4">
,<"F5">
,<"F6">
,<"H0">
,<"F7">
,<"H1">
,<"F8">
,<"H2">
,<"F9">
,<"H3">
,<"H4">
,<"H5">
,<"E_0">
,<"H6">
,<"E_1">
,<"H7">
,<"E_2">
,<"H8">
,<"E_3">
,<"H9">
,<"E_4">
,<"E_5">
,<"E_6">
,<"E_7">
,<"E_8">
,<"L0">
,<"E_9">
,<"L1">
,<"L2">
,<"L3">
,<"L4">
,<"L5">
,<"L6">
,<"L7">
,<"L8">
,<"L9">
,<"A0">
,<"A1">
,<"A2">
,<"A3">
,<"A4">
,<"A5">
,<"A6">
,<"C0">
,<"A7">
,<"C1">
,<"A8">
,<"A9">
,<"C2">
,<"C3">
,<"C4">
,<"C5">
,<"C6">
,<"C7">
,<"C8">
,<"C9">
,<"G0">
,<"G1">
,<"G2">
,<"G3">
,<"G4">
,<"G5">
,<"G6">
,<"I0">
,<"G7">
,<"I1">
,<"G8">
,<"I2">
,<"G9">
,<"I3">
,<"I4">
,<"I5">
,<"I6">
,<"I7">
,<"I8">
,<"I9">
};
set CVAR := {
};
set V := IVAR + CVAR;
param lower[V] := 
<"B0">0
,<"B1">0
,<"B2">0
,<"B3">0
,<"B4">0
,<"B5">0
,<"B6">0
,<"D0">0
,<"B7">0
,<"B8">0
,<"D1">0
,<"B9">0
,<"D2">0
,<"D3">0
,<"D4">0
,<"D5">0
,<"D6">0
,<"F0">0
,<"D7">0
,<"F1">0
,<"D8">0
,<"F2">0
,<"D9">0
,<"F3">0
,<"F4">0
,<"F5">0
,<"F6">0
,<"H0">0
,<"F7">0
,<"H1">0
,<"F8">0
,<"H2">0
,<"F9">0
,<"H3">0
,<"H4">0
,<"H5">0
,<"E_0">0
,<"H6">0
,<"E_1">0
,<"H7">0
,<"E_2">0
,<"H8">0
,<"E_3">0
,<"H9">0
,<"E_4">0
,<"E_5">0
,<"E_6">0
,<"E_7">0
,<"E_8">0
,<"L0">0
,<"E_9">0
,<"L1">0
,<"L2">0
,<"L3">0
,<"L4">0
,<"L5">0
,<"L6">0
,<"L7">0
,<"L8">0
,<"L9">0
,<"A0">0
,<"A1">0
,<"A2">0
,<"A3">0
,<"A4">0
,<"A5">0
,<"A6">0
,<"C0">0
,<"A7">0
,<"C1">0
,<"A8">0
,<"A9">0
,<"C2">0
,<"C3">0
,<"C4">0
,<"C5">0
,<"C6">0
,<"C7">0
,<"C8">0
,<"C9">0
,<"G0">0
,<"G1">0
,<"G2">0
,<"G3">0
,<"G4">0
,<"G5">0
,<"G6">0
,<"I0">0
,<"G7">0
,<"I1">0
,<"G8">0
,<"I2">0
,<"G9">0
,<"I3">0
,<"I4">0
,<"I5">0
,<"I6">0
,<"I7">0
,<"I8">0
,<"I9">0
;
param upper[V] := 
<"B0">1
,<"B1">1
,<"B2">1
,<"B3">1
,<"B4">1
,<"B5">1
,<"B6">1
,<"D0">1
,<"B7">1
,<"B8">1
,<"D1">1
,<"B9">1
,<"D2">1
,<"D3">1
,<"D4">1
,<"D5">1
,<"D6">1
,<"F0">1
,<"D7">1
,<"F1">1
,<"D8">1
,<"F2">1
,<"D9">1
,<"F3">1
,<"F4">1
,<"F5">1
,<"F6">1
,<"H0">1
,<"F7">1
,<"H1">1
,<"F8">1
,<"H2">1
,<"F9">1
,<"H3">1
,<"H4">1
,<"H5">1
,<"E_0">1
,<"H6">1
,<"E_1">1
,<"H7">1
,<"E_2">1
,<"H8">1
,<"E_3">1
,<"H9">1
,<"E_4">1
,<"E_5">1
,<"E_6">1
,<"E_7">1
,<"E_8">1
,<"L0">1
,<"E_9">1
,<"L1">1
,<"L2">1
,<"L3">1
,<"L4">1
,<"L5">1
,<"L6">1
,<"L7">1
,<"L8">1
,<"L9">1
,<"A0">1
,<"A1">1
,<"A2">1
,<"A3">1
,<"A4">1
,<"A5">1
,<"A6">1
,<"C0">1
,<"A7">1
,<"C1">1
,<"A8">1
,<"A9">1
,<"C2">1
,<"C3">1
,<"C4">1
,<"C5">1
,<"C6">1
,<"C7">1
,<"C8">1
,<"C9">1
,<"G0">1
,<"G1">1
,<"G2">1
,<"G3">1
,<"G4">1
,<"G5">1
,<"G6">1
,<"I0">1
,<"G7">1
,<"I1">1
,<"G8">1
,<"I2">1
,<"G9">1
,<"I3">1
,<"I4">1
,<"I5">1
,<"I6">1
,<"I7">1
,<"I8">1
,<"I9">1
;
set RANGEROW := {
<" dummy ">
};
set ROW := { 
<"SOSC","E">
,<"SOS6","E">
,<"SOSD","E">
,<"SOS7","E">
,<"OBJECT","N">
,<"SOSE","E">
,<"SOS8","E">
,<"SOSF","E">
,<"SOS9","E">
,<"SOSG","E">
,<"BILANCIO","E">
,<"SOSH","E">
,<"SOSI","E">
,<"SOS0","E">
,<"SOS1","E">
,<"SOSL","E">
,<"SOS2","E">
,<"SOS3","E">
,<"SOSA","E">
,<"SOS4","E">
,<"SOSB","E">
,<"SOS5","E">
};
set R := proj(ROW, <1>);
param rhs[R] := 
<"SOSC">1
,<"SOS6">1
,<"SOSD">1
,<"SOS7">1
,<"OBJECT">0
,<"SOSE">1
,<"SOS8">1
,<"SOSF">1
,<"SOS9">1
,<"SOSG">1
,<"BILANCIO">0
,<"SOSH">1
,<"SOSI">1
,<"SOS0">1
,<"SOS1">1
,<"SOSL">1
,<"SOS2">1
,<"SOS3">1
,<"SOSA">1
,<"SOS4">1
,<"SOSB">1
,<"SOS5">1
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

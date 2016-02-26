objective value:                       0.001032940184
x_1_1                                               1 	(obj:0)
x_2_2                                               1 	(obj:0)
x_3_1                                               1 	(obj:0)
x_4_2                                               1 	(obj:0)
x_5_3                                               1 	(obj:0)
x_6_1                                               1 	(obj:0)
x_7_2                                               1 	(obj:0)
x_8_2                                               1 	(obj:0)
x_9_1                                               1 	(obj:0)
x_10_1                                              1 	(obj:0)
x_11_3                                              1 	(obj:0)
x_12_3                                              1 	(obj:0)
x_13_1                                              1 	(obj:0)
x_14_2                                              1 	(obj:0)
x_15_2                                              1 	(obj:0)
x_16_1                                              1 	(obj:0)
x_17_3                                              1 	(obj:0)
x_18_2                                              1 	(obj:0)
x_19_3                                              1 	(obj:0)
x_20_1                                              1 	(obj:0)
ind_1                                               1 	(obj:0)
abs_1_3                                             1 	(obj:0)
ind_2                                               1 	(obj:0)
abs_2_1                                             1 	(obj:0)
ind_3                                               1 	(obj:0)
abs_3_2                                             1 	(obj:0)
epsI                                   0.001032940184 	(obj:1)
SCIP Status        : solving was interrupted [user interrupt]
Total Time         :      59.82
  solving          :      59.82
  presolving       :       0.05 (included in solving)
  reading          :       0.00
  copying          :       0.01 (1 #copies) (minimal 0.01, maximal 0.01, average 0.01)
Original Problem   :
  Problem name     : instances/Pot6Cycle_T_20B_1t.txt.spa
  Variables        : 101 (100 binary, 0 integer, 0 implicit integer, 1 continuous)
  Constraints      : 154 initial, 154 maximal
  Objective sense  : maximize
Presolved Problem  :
  Problem name     : t_instances/Pot6Cycle_T_20B_1t.txt.spa
  Variables        : 559 (87 binary, 0 integer, 16 implicit integer, 456 continuous)
  Constraints      : 1873 initial, 3749 maximal
Presolvers         :   ExecTime  SetupTime  Calls  FixedVars   AggrVars   ChgTypes  ChgBounds   AddHoles    DelCons    AddCons   ChgSides   ChgCoefs
  boundshift       :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  components       :       0.00       0.00      1          0          0          0          0          0          0          0          0          0
  convertinttobin  :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  domcol           :       0.01       0.00      1          0          0          0          0          0          0          0          0          0
  dualagg          :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  dualinfer        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  gateextraction   :       0.00       0.00      1          0          0          0          0          0        108         34          0          0
  implfree         :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  implics          :       0.00       0.00      6          0          0          0          0          0          0          0          0          0
  inttobinary      :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  redvub           :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  stuffing         :       0.00       0.00      1          0          0          0          0          0          0          0          0          0
  trivial          :       0.00       0.00     11          2          0          0          0          0          0          0          0          0
  tworowbnd        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  dualfix          :       0.00       0.00     11          4          0          0          0          0          0          0          0          0
  genvbounds       :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  probing          :       0.02       0.00      1          0          0          0          0          0          0          0          0          0
  pseudoobj        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  varbound         :       0.00       0.00     15          0          0          0         78          0        119         28          0          0
  setppc           :       0.00       0.00     19          0          0          0          0          0          0          0          0          0
  and              :       0.00       0.00      4          0         12          0          0          0         12          0          0          0
  linear           :       0.01       0.00     14         42          8          0          1          0        132          0         25         25
  orbitope         :       0.00       0.00      6          6          0          0          0          0          0          0          0          0
  logicor          :       0.00       0.00      2          0          0          0          0          0          0          0          0          0
  bounddisjunction :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  quadratic        :       0.01       0.00      2          0          0          0          0          0          0       2028          0          0
  root node        :          -          -      -          0          -          -         70          -          -          -          -          -
Constraints        :     Number  MaxNumber  #Separate #Propagate    #EnfoLP    #EnfoPS     #Check   #ResProp    Cutoffs    DomReds       Cuts    Applied      Conss   Children
  integral         :          0          0          0          0      13664          0      37673          0          0       4307          0          0          0      26248
  varbound         :        845        845      19785      94943        115          0      20081      22966        323     239603          0          0          0          0
  setppc           :         81+        86      19785      94619        115          0        355      37739       1013      18493         33         11          0          0
  and              :         19         19        515       6703          0          0        144       1429          1        401         26          6          0          0
  linear           :        927        927      19785      93606        115          0        140      23569        784     461740     250781     116428          0          0
  orbitope         :          1          1      19867          0          0          0          0          0        115          0       2200        488          0          0
  logicor          :          0+       988      18569      22409          0          0          0       1456         56        264         19         16          0          0
  bounddisjunction :          0+       907          0      28461          0          0          0        311          2       1877          0          0          0          0
  countsols        :          0          0          0          0          0          0      17446          0          0          0          0          0          0          0
Constraint Timings :  TotalTime  SetupTime   Separate  Propagate     EnfoLP     EnfoPS      Check    ResProp    SB-Prop
  integral         :       2.28       0.00       0.00       0.00       2.28       0.00       0.00       0.00       0.00
  varbound         :       0.70       0.00       0.19       0.36       0.00       0.00       0.13       0.00       0.02
  setppc           :       0.74       0.00       0.00       0.67       0.00       0.00       0.00       0.01       0.06
  and              :       0.01       0.00       0.01       0.00       0.00       0.00       0.00       0.00       0.00
  linear           :       2.53       0.00       0.74       1.70       0.01       0.00       0.00       0.02       0.06
  orbitope         :       0.05       0.00       0.05       0.00       0.00       0.00       0.00       0.00       0.00
  logicor          :       0.33       0.00       0.13       0.20       0.00       0.00       0.00       0.00       0.00
  bounddisjunction :       0.30       0.00       0.00       0.29       0.00       0.00       0.00       0.00       0.01
  countsols        :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
Propagators        : #Propagate   #ResProp    Cutoffs    DomReds
  dualfix          :          2          0          0          0
  genvbounds       :          0          0          0          0
  obbt             :          0          0          0          0
  probing          :          0          0          0          0
  pseudoobj        :      66817          0          0         70
  redcost          :      46330          0          0        216
  rootredcost      :          1          0          0          0
  vbounds          :     108186          2          0         15
Propagator Timings :  TotalTime  SetupTime   Presolve  Propagate    ResProp    SB-Prop
  dualfix          :       0.00       0.00       0.00       0.00       0.00       0.00
  genvbounds       :       0.01       0.00       0.00       0.01       0.00       0.00
  obbt             :       0.00       0.00       0.00       0.00       0.00       0.00
  probing          :       0.02       0.00       0.02       0.00       0.00       0.00
  pseudoobj        :       0.02       0.00       0.00       0.02       0.00       0.00
  redcost          :       0.05       0.00       0.00       0.05       0.00       0.00
  rootredcost      :       0.01       0.00       0.00       0.01       0.00       0.00
  vbounds          :       0.03       0.00       0.00       0.03       0.00       0.00
Conflict Analysis  :       Time      Calls    Success    DomReds  Conflicts   Literals    Reconvs ReconvLits   LP Iters
  propagation      :       0.12       1541       1393          -      21796       40.1        107       48.7          -
  infeasible LP    :       0.03        233        152          -       1793       35.8          0        0.0          0
  bound exceed. LP :       0.00          0          0          -          0        0.0          0        0.0          0
  strong branching :       0.00          0          0          -          0        0.0          0        0.0          0
  pseudo solution  :       0.00          0          0          -          0        0.0          0        0.0          -
  applied globally :       0.06          -          -          0       3093       21.0          -          -          -
  applied locally  :          -          -          -          0         28       33.8          -          -          -
Separators         :   ExecTime  SetupTime      Calls    Cutoffs    DomReds       Cuts    Applied      Conss
  cut pool         :       0.03                 19780          -          -      17495          -          -    (maximal pool size: 842)
  cgmip            :       0.00       0.00          0          0          0          0          0          0
  clique           :       0.02       0.00        173          0          0          6          6          0
  closecuts        :       0.05       0.00          5          0          0        260          0          0
  cmir             :       0.03       0.00         15          0          0        440      15944          0
  disjunctive      :       0.00       0.00          0          0          0          0          0          0
  eccuts           :       0.00       0.00          0          0          0          0          0          0
  flowcover        :       0.14       0.00         22          0          0        320         26          0
  gomory           :       0.07       0.00         15          0          0          1          1          0
  impliedbounds    :       0.01       0.00        173          0          0       1336        236          0
  intobj           :       0.00       0.00          0          0          0          0          0          0
  mcf              :       0.00       0.00          1          0          0          0          0          0
  oddcycle         :       0.00       0.00          8          0          0          2          2          0
  rapidlearning    :       0.00       0.00          0          0          0          0          0          0
  strongcg         :       0.05       0.00         30          0          0          3          1          0
  zerohalf         :       0.04       0.00         15          0          0          0          0          0
Pricers            :   ExecTime  SetupTime      Calls       Vars
  problem variables:       0.00          -          0          0
Branching Rules    :   ExecTime  SetupTime   BranchLP  BranchExt   BranchPS    Cutoffs    DomReds       Cuts      Conss   Children
  allfullstrong    :       0.00       0.00          0          0          0          0          0          0          0          0
  cloud            :       0.00       0.00          0          0          0          0          0          0          0          0
  distribution     :       0.00       0.00          0          0          0          0          0          0          0          0
  fullstrong       :       0.00       0.00          0          0          0          0          0          0          0          0
  inference        :       0.00       0.00          0          0          0          0          0          0          0          0
  leastinf         :       0.00       0.00          0          0          0          0          0          0          0          0
  mostinf          :       0.00       0.00          0          0          0          0          0          0          0          0
  multaggr         :       0.00       0.00          0          0          0          0          0          0          0          0
  nodereopt        :       0.00       0.00          0          0          0          0          0          0          0          0
  pscost           :       0.00       0.00          0          0          0          0          0          0          0          0
  random           :       0.00       0.00          0          0          0          0          0          0          0          0
  relpscost        :       2.27       0.00      13549          0          0          0       4307          0        -18      26248
Primal Heuristics  :   ExecTime  SetupTime      Calls      Found       Best
  LP solutions     :       0.00          -          -          0          -
  pseudo solutions :       0.00          -          -          0          -
  strong branching :       0.00          -          -          0          -
  SpaGreedy        :      20.60       0.00      17445          2          1
  actconsdiving    :       0.00       0.00          0          0          0
  bound            :       0.00       0.00          0          0          0
  clique           :       0.00       0.00          0          0          0
  coefdiving       :       1.72       0.00        743          0          0
  crossover        :       0.00       0.00          0          0          0
  dins             :       0.00       0.00          0          0          0
  distributiondivin:       0.61       0.00        146          0          0
  dualval          :       0.00       0.00          0          0          0
  feaspump         :       0.31       0.00         50          0          0
  fixandinfer      :       0.00       0.00          0          0          0
  fracdiving       :       0.51       0.00        146          0          0
  guideddiving     :       0.00       0.00          0          0          0
  indicator        :       0.00       0.00          0          0          0
  intdiving        :       0.00       0.00          0          0          0
  intshifting      :       0.13       0.00         59          0          0
  linesearchdiving :       1.01       0.00        266          0          0
  localbranching   :       0.00       0.00          0          0          0
  mutation         :       0.00       0.00          0          0          0
  nlpdiving        :       0.00       0.00          0          0          0
  objpscostdiving  :       0.32       0.00         21          0          0
  octane           :       0.00       0.00          0          0          0
  ofins            :       0.00       0.00          0          0          0
  oneopt           :       0.00       0.00          0          0          0
  proximity        :       0.00       0.00          0          0          0
  pscostdiving     :       1.41       0.00        494          0          0
  randrounding     :       0.67       0.00      19901          0          0
  rens             :       0.70       0.00          1          0          0
  reoptsols        :       0.00       0.00          0          0          0
  rins             :       0.00       0.00          0          0          0
  rootsoldiving    :       0.25       0.00         23          0          0
  rounding         :       0.08       0.00       2550          0          0
  shiftandpropagate:       0.00       0.00          0          0          0
  shifting         :       0.53       0.00       2030          0          0
  simplerounding   :       0.00       0.00          0          0          0
  subnlp           :       0.00       0.00          0          0          0
  trivial          :       0.00       0.00          2          0          0
  trivialnegation  :       0.00       0.00          0          0          0
  trysol           :       0.00       0.00          0          0          0
  twoopt           :       0.00       0.00          0          0          0
  undercover       :       0.00       0.00          1          0          0
  vbounds          :       0.00       0.00          0          0          0
  veclendiving     :       0.39       0.00        143          0          0
  zeroobj          :       0.00       0.00          0          0          0
  zirounding       :       0.01       0.00       1000          0          0
  other solutions  :          -          -          -          0          -
Diving Statistics  :      Calls      Nodes   LP Iters Backtracks   MinDepth   MaxDepth   AvgDepth  NLeafSols  MinSolDpt  MaxSolDpt  AvgSolDpt
  actconsdiving    :          0          -          -          -          -          -          -          -          -          -          -
  coefdiving       :        743       4909      41153        838         17         43       26.3          -          -          -          -
  distributiondivin:        146       1145      13989        166         20         37       25.8          -          -          -          -
  fracdiving       :        146       1265      13575        151         20         48       27.1          -          -          -          -
  guideddiving     :          0          -          -          -          -          -          -          -          -          -          -
  linesearchdiving :        266       1436      27459        252         19         34       24.4          -          -          -          -
  pscostdiving     :        494       2795      41465        479         20         53       25.6          -          -          -          -
  veclendiving     :        143        914      11085        163         19         31       24.8          -          -          -          -
LP                 :       Time      Calls Iterations  Iter/call   Iter/sec  Time-0-It Calls-0-It
  primal LP        :       0.02        218         19       1.36     950.00       0.02        204
  dual LP          :      17.03      39459     790949      21.22   46444.45       0.26       2188
  lex dual LP      :       0.00          0          0       0.00          -
  barrier LP       :       0.00          0          0       0.00          -       0.00          0
  diving/probing LP:       4.06       5081     185534      36.52   45698.03
  strong branching :       1.70        857      58550      68.32   34441.18
    (at root node) :          -         32       3381     105.66          -
  conflict analysis:       0.00          0          0       0.00          -
B&B Tree           :
  number of runs   :          1
  nodes            :      17445 (13124 internal, 4321 leaves)
  nodes (total)    :      17445 (13124 internal, 4321 leaves)
  nodes left       :       8792
  max depth        :         24
  max depth (total):         24
  backtracks       :       4811 (27.6%)
  delayed cutoffs  :         13
  repropagations   :        211 (2379 domain reductions, 9 cutoffs)
  avg switch length:       5.24
  switching time   :       3.03
Root Node          :
  First LP value   : +3.00000000000000e+00
  First LP Iters   :        354 (35400.00 Iter/sec)
  First LP Time    :       0.01
  Final Dual Bound : +1.18800689568374e-02
  Final Root Iters :       6136
Solution           :
  Solutions found  :          2 (1 improvements)
  First Solution   : +1.03294018400000e-03   (in run 1, after 1 nodes, 0.05 seconds, depth 0, found by <SpaGreedy>)
  Gap First Sol.   :   infinite
  Gap Last Sol.    :   infinite
  Primal Bound     : +1.03294018400000e-03   (in run 1, after 1 nodes, 0.05 seconds, depth 0, found by <SpaGreedy>)
  Dual Bound       : +3.23416882769396e-03
  Gap              :     213.10 %
  Avg. Gap         :      91.49 % (5472.94 primal-dual integral)

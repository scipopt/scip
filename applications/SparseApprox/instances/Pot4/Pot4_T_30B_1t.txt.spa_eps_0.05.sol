objective value:                 0.000637853732999988
x_1_1                                               1 	(obj:0)
x_2_2                                               1 	(obj:0)
x_3_4                                               1 	(obj:0)
x_4_2                                               1 	(obj:0)
x_5_1                                               1 	(obj:0)
x_6_4                                               1 	(obj:0)
x_7_2                                               1 	(obj:0)
x_8_1                                               1 	(obj:0)
x_9_1                                               1 	(obj:0)
x_10_4                                              1 	(obj:0)
x_11_2                                              1 	(obj:0)
x_12_4                                              1 	(obj:0)
x_13_1                                              1 	(obj:0)
x_14_2                                              1 	(obj:0)
x_15_1                                              1 	(obj:0)
x_16_4                                              1 	(obj:0)
x_17_2                                              1 	(obj:0)
x_18_2                                              1 	(obj:0)
x_19_4                                              1 	(obj:0)
x_20_1                                              1 	(obj:0)
x_21_2                                              1 	(obj:0)
x_22_1                                              1 	(obj:0)
x_23_4                                              1 	(obj:0)
x_24_1                                              1 	(obj:0)
x_25_1                                              1 	(obj:0)
x_26_2                                              1 	(obj:0)
x_27_4                                              1 	(obj:0)
x_28_2                                              1 	(obj:0)
x_29_2                                              1 	(obj:0)
x_30_4                                              1 	(obj:0)
ind_1                                               1 	(obj:0)
abs_1_4                                             1 	(obj:0)
ind_2                                               1 	(obj:0)
abs_2_1                                             1 	(obj:0)
ind_4                                               1 	(obj:0)
abs_4_2                                             1 	(obj:0)
epsI                             0.000637853732999988 	(obj:1)
SCIP Status        : solving was interrupted [node limit reached]
Total Time         :       5.33
  solving          :       5.32
  presolving       :       0.09 (included in solving)
  reading          :       0.01
  copying          :       0.01 (1 #copies) (minimal 0.01, maximal 0.01, average 0.01)
Original Problem   :
  Problem name     : ../instances/Pot4/Pot4_T_30B_1t.txt.spa
  Variables        : 141 (140 binary, 0 integer, 0 implicit integer, 1 continuous)
  Constraints      : 263 initial, 263 maximal
  Objective sense  : maximize
Presolved Problem  :
  Problem name     : t_../instances/Pot4/Pot4_T_30B_1t.txt.spa
  Variables        : 910 (132 binary, 0 integer, 16 implicit integer, 762 continuous)
  Constraints      : 3165 initial, 3165 maximal
Presolvers         :   ExecTime  SetupTime  Calls  FixedVars   AggrVars   ChgTypes  ChgBounds   AddHoles    DelCons    AddCons   ChgSides   ChgCoefs
  boundshift       :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  components       :       0.00       0.00      1          0          0          0          0          0          0          0          0          0
  convertinttobin  :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  domcol           :       0.00       0.00      1          0          0          0          0          0          0          0          0          0
  dualagg          :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  dualinfer        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  gateextraction   :       0.00       0.00      3          0          0          0          0          0        108         34          0          0
  implfree         :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  implics          :       0.00       0.00      4          0          0          0          0          0          0          0          0          0
  inttobinary      :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  redvub           :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  stuffing         :       0.00       0.00      1          0          0          0          0          0          0          0          0          0
  trivial          :       0.00       0.00      7          0          0          0          0          0          0          0          0          0
  tworowbnd        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  dualfix          :       0.00       0.00      7          4          0          0          0          0          0          0          0          0
  genvbounds       :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  probing          :       0.06       0.00      1          0          0          0          4          0          0          0          0          0
  pseudoobj        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  varbound         :       0.00       0.00      9          0          0          0          0          0         28         28          0          0
  setppc           :       0.01       0.00     12          0          0          0          0          0          0          0          0          0
  and              :       0.00       0.00      4          0         12          0          0          0         12          0          0          0
  linear           :       0.01       0.00      9          3          1          0          4          0          7          0         12         12
  logicor          :       0.00       0.00      7          0          0          0          0          0          0          0          0          0
  bounddisjunction :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  quadratic        :       0.01       0.00      2          0          0          0          0          0          0       2995          0          0
  root node        :          -          -      -          1          -          -         61          -          -          -          -          -
Constraints        :     Number  MaxNumber  #Separate #Propagate    #EnfoLP    #EnfoPS     #Check   #ResProp    Cutoffs    DomReds       Cuts    Applied      Conss   Children
  integral         :          0          0          0          0        359          0         76          0          0       2328          0          0          0        198
  varbound         :       1389       1389         75       6392          0          0         11        109          6       3053          0          0          0          0
  setppc           :        160        160         75       6386          0          0          6        189          4        306         10          0          0          0
  and              :         22         22        160       4886          0          0          4          9          0         92         28          9          0          0
  linear           :       1539       1539         75       6382          0          0          2         79          6       2550        681        263          0          0
  logicor          :         55+        59         65        929          0          0          0         27          5         34          0          0          0          0
  bounddisjunction :          0+         6          0        376          0          0          0          0          0          0          0          0          0          0
  countsols        :          0          0          0          0          0          0          9          0          0          0          0          0          0          0
Constraint Timings :  TotalTime  SetupTime   Separate  Propagate     EnfoLP     EnfoPS      Check    ResProp    SB-Prop
  integral         :       3.13       0.00       0.00       0.00       3.13       0.00       0.00       0.00       0.00
  varbound         :       0.05       0.00       0.01       0.03       0.00       0.00       0.00       0.00       0.01
  setppc           :       0.06       0.00       0.00       0.04       0.00       0.00       0.00       0.00       0.02
  and              :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
  linear           :       0.23       0.00       0.01       0.16       0.00       0.00       0.00       0.00       0.06
  logicor          :       0.01       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.01
  bounddisjunction :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
  countsols        :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
Propagators        : #Propagate   #ResProp    Cutoffs    DomReds
  dualfix          :          2          0          0          0
  genvbounds       :          0          0          0          0
  obbt             :          0          0          0          0
  probing          :          0          0          0          0
  pseudoobj        :       1051          0          0         60
  redcost          :        796          0          0          1
  rootredcost      :          1          0          0          0
  vbounds          :       6195          0          0          0
Propagator Timings :  TotalTime  SetupTime   Presolve  Propagate    ResProp    SB-Prop
  dualfix          :       0.00       0.00       0.00       0.00       0.00       0.00
  genvbounds       :       0.00       0.00       0.00       0.00       0.00       0.00
  obbt             :       0.00       0.00       0.00       0.00       0.00       0.00
  probing          :       0.06       0.00       0.06       0.00       0.00       0.00
  pseudoobj        :       0.00       0.00       0.00       0.00       0.00       0.00
  redcost          :       0.00       0.00       0.00       0.00       0.00       0.00
  rootredcost      :       0.00       0.00       0.00       0.00       0.00       0.00
  vbounds          :       0.00       0.00       0.00       0.00       0.00       0.00
Conflict Analysis  :       Time      Calls    Success    DomReds  Conflicts   Literals    Reconvs ReconvLits   LP Iters
  propagation      :       0.00         15         14          -        100       54.8          3       60.3          -
  infeasible LP    :       0.00         10          0          -          0        0.0          0        0.0          0
  bound exceed. LP :       0.00          0          0          -          0        0.0          0        0.0          0
  strong branching :       0.00          0          0          -          0        0.0          0        0.0          0
  pseudo solution  :       0.00          0          0          -          0        0.0          0        0.0          -
  applied globally :       0.00          -          -          0         14        9.5          -          -          -
  applied locally  :          -          -          -          0          2       40.5          -          -          -
Separators         :   ExecTime  SetupTime      Calls    Cutoffs    DomReds       Cuts    Applied      Conss
  cut pool         :       0.00                    73          -          -       1512          -          -    (maximal pool size: 911)
  cgmip            :       0.00       0.00          0          0          0          0          0          0
  clique           :       0.00       0.00         74          0          0         68         34          0
  closecuts        :       0.00       0.00          0          0          0          0          0          0
  cmir             :       0.02       0.00         10          0          0        481        192          0
  disjunctive      :       0.00       0.00          0          0          0          0          0          0
  eccuts           :       0.00       0.00          0          0          0          0          0          0
  flowcover        :       0.38       0.00         15          0          0        278         37          0
  gomory           :       0.05       0.00         10          0          0          4          4          0
  impliedbounds    :       0.00       0.00         74          0          0        798        155          0
  intobj           :       0.00       0.00          0          0          0          0          0          0
  mcf              :       0.00       0.00          1          0          0          0          0          0
  oddcycle         :       0.00       0.00          0          0          0          0          0          0
  rapidlearning    :       0.00       0.00          0          0          0          0          0          0
  strongcg         :       0.08       0.00         20          0          0         10          7          0
  zerohalf         :       0.00       0.00          0          0          0          0          0          0
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
  relpscost        :       3.13       0.00        359          0          0          0       2328          0          0        198
Primal Heuristics  :   ExecTime  SetupTime      Calls      Found       Best
  LP solutions     :       0.00          -          -          0          -
  pseudo solutions :       0.00          -          -          0          -
  strong branching :       0.00          -          -          0          -
  SpaGreedy        :       0.00       0.00          8          2          1
  actconsdiving    :       0.00       0.00          0          0          0
  bound            :       0.00       0.00          0          0          0
  clique           :       0.00       0.00          0          0          0
  coefdiving       :       0.01       0.00          1          0          0
  crossover        :       0.00       0.00          0          0          0
  dins             :       0.00       0.00          0          0          0
  distributiondivin:       0.01       0.00          1          0          0
  dualval          :       0.00       0.00          0          0          0
  feaspump         :       0.02       0.00          1          0          0
  fixandinfer      :       0.00       0.00          0          0          0
  fracdiving       :       0.02       0.00          1          0          0
  guideddiving     :       0.00       0.00          0          0          0
  indicator        :       0.00       0.00          0          0          0
  intdiving        :       0.00       0.00          0          0          0
  intshifting      :       0.01       0.00          2          0          0
  linesearchdiving :       0.00       0.00          0          0          0
  localbranching   :       0.00       0.00          0          0          0
  mutation         :       0.00       0.00          0          0          0
  nlpdiving        :       0.00       0.00          0          0          0
  objpscostdiving  :       0.04       0.00          1          0          0
  octane           :       0.00       0.00          0          0          0
  ofins            :       0.00       0.00          0          0          0
  oneopt           :       0.00       0.00          0          0          0
  proximity        :       0.00       0.00          0          0          0
  pscostdiving     :       0.02       0.00          1          0          0
  randrounding     :       0.37       0.00         74          0          0
  rens             :       0.05       0.00          1          0          0
  reoptsols        :       0.00       0.00          0          0          0
  rins             :       0.00       0.00          0          0          0
  rootsoldiving    :       0.00       0.00          1          0          0
  rounding         :       0.01       0.00        207          0          0
  shiftandpropagate:       0.00       0.00          0          0          0
  shifting         :       0.09       0.00        104          0          0
  simplerounding   :       0.00       0.00          0          0          0
  subnlp           :       0.00       0.00          0          0          0
  trivial          :       0.00       0.00          2          0          0
  trivialnegation  :       0.00       0.00          0          0          0
  trysol           :       0.00       0.00          0          0          0
  twoopt           :       0.00       0.00          0          0          0
  undercover       :       0.00       0.00          1          0          0
  vbounds          :       0.00       0.00          0          0          0
  veclendiving     :       0.02       0.00          1          0          0
  zeroobj          :       0.00       0.00          0          0          0
  zirounding       :       0.00       0.00         98          0          0
  other solutions  :          -          -          -          0          -
Diving Statistics  :      Calls      Nodes   LP Iters Backtracks   MinDepth   MaxDepth   AvgDepth  NLeafSols  MinSolDpt  MaxSolDpt  AvgSolDpt
  actconsdiving    :          0          -          -          -          -          -          -          -          -          -          -
  coefdiving       :          1         21        221          2         35         35       35.0          -          -          -          -
  distributiondivin:          1         15        320          1         34         34       34.0          -          -          -          -
  fracdiving       :          1         18        435          2         33         33       33.0          -          -          -          -
  guideddiving     :          0          -          -          -          -          -          -          -          -          -          -
  linesearchdiving :          0          -          -          -          -          -          -          -          -          -          -
  pscostdiving     :          1         14        218          0         30         30       30.0          -          -          -          -
  veclendiving     :          1         22        414          4         36         36       36.0          -          -          -          -
LP                 :       Time      Calls Iterations  Iter/call   Iter/sec  Time-0-It Calls-0-It
  primal LP        :       0.00         11          0       0.00          -       0.00         11
  dual LP          :       0.79        540      16244      43.09   20562.03       0.02        163
  lex dual LP      :       0.00          0          0       0.00          -
  barrier LP       :       0.00          0          0       0.00          -       0.00          0
  diving/probing LP:       0.21        101       3475      34.41   16547.62
  strong branching :       2.72        525      59872     114.04   22011.76
    (at root node) :          -         23       2036      88.52          -
  conflict analysis:       0.00          0          0       0.00          -
B&B Tree           :
  number of runs   :          1
  nodes            :        100 (99 internal, 1 leaves)
  nodes (total)    :        100 (99 internal, 1 leaves)
  nodes left       :         99
  max depth        :         19
  max depth (total):         19
  backtracks       :          7 (7.0%)
  delayed cutoffs  :          0
  repropagations   :          2 (0 domain reductions, 0 cutoffs)
  avg switch length:       1.95
  switching time   :       0.00
Root Node          :
  First LP value   : +3.00000000000000e+00
  First LP Iters   :        622 (62200.00 Iter/sec)
  First LP Time    :       0.01
  Final Dual Bound : +5.94623284467778e-03
  Final Root Iters :       3797
Solution           :
  Solutions found  :          2 (1 improvements)
  First Solution   : +6.37853732999988e-04   (in run 1, after 1 nodes, 0.09 seconds, depth 0, found by <SpaGreedy>)
  Gap First Sol.   :   infinite
  Gap Last Sol.    :   infinite
  Primal Bound     : +6.37853732999988e-04   (in run 1, after 1 nodes, 0.09 seconds, depth 0, found by <SpaGreedy>)
  Dual Bound       : +5.31056205004122e-03
  Gap              :     732.57 %
  Avg. Gap         :      91.67 % (487.71 primal-dual integral)

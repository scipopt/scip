objective value:                 0.000787434402999998
x_1_1                                               1 	(obj:0)
x_2_2                                               1 	(obj:0)
x_3_4                                               1 	(obj:0)
x_4_1                                               1 	(obj:0)
x_5_1                                               1 	(obj:0)
x_6_2                                               1 	(obj:0)
x_7_2                                               1 	(obj:0)
x_8_4                                               1 	(obj:0)
x_9_1                                               1 	(obj:0)
x_10_2                                              1 	(obj:0)
x_11_4                                              1 	(obj:0)
x_12_1                                              1 	(obj:0)
x_13_2                                              1 	(obj:0)
x_14_4                                              1 	(obj:0)
x_15_1                                              1 	(obj:0)
x_16_4                                              1 	(obj:0)
x_17_2                                              1 	(obj:0)
x_18_2                                              1 	(obj:0)
x_19_4                                              1 	(obj:0)
x_20_4                                              1 	(obj:0)
ind_1                                               1 	(obj:0)
abs_1_2                                             1 	(obj:0)
ind_2                                               1 	(obj:0)
abs_2_4                                             1 	(obj:0)
ind_4                                               1 	(obj:0)
abs_4_1                                             1 	(obj:0)
epsI                             0.000787434402999998 	(obj:1)
SCIP Status        : solving was interrupted [node limit reached]
Total Time         :       2.88
  solving          :       2.88
  presolving       :       0.05 (included in solving)
  reading          :       0.00
  copying          :       0.01 (1 #copies) (minimal 0.01, maximal 0.01, average 0.01)
Original Problem   :
  Problem name     : ../instances/Pot4/Pot4_T_20B_alpha=1.txt.spa
  Variables        : 101 (100 binary, 0 integer, 0 implicit integer, 1 continuous)
  Constraints      : 193 initial, 193 maximal
  Objective sense  : maximize
Presolved Problem  :
  Problem name     : t_../instances/Pot4/Pot4_T_20B_alpha=1.txt.spa
  Variables        : 590 (92 binary, 0 integer, 16 implicit integer, 482 continuous)
  Constraints      : 2015 initial, 2037 maximal
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
  probing          :       0.03       0.00      1          0          0          0          4          0          0          0          0          0
  pseudoobj        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  varbound         :       0.00       0.00      9          0          0          0          0          0         28         28          0          0
  setppc           :       0.01       0.00     12          0          0          0          0          0          0          0          0          0
  and              :       0.00       0.00      4          0         12          0          0          0         12          0          0          0
  linear           :       0.00       0.00      9          3          1          0          4          0          7          0         12         12
  logicor          :       0.00       0.00      7          0          0          0          0          0          0          0          0          0
  bounddisjunction :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  quadratic        :       0.01       0.00      2          0          0          0          0          0          0       1915          0          0
  root node        :          -          -      -          1          -          -         57          -          -          -          -          -
Constraints        :     Number  MaxNumber  #Separate #Propagate    #EnfoLP    #EnfoPS     #Check   #ResProp    Cutoffs    DomReds       Cuts    Applied      Conss   Children
  integral         :          0          0          0          0        239          0         52          0          0       1423          0          0          0        180
  varbound         :        869        869         84       5019          7          0         15        349         13       1563          0          0          0          0
  setppc           :        110        110         84       5006          7          0         11        311          9        143         10          0          0          0
  and              :         22         22        125       3902          0          0          9         29          3         56         24          7          0          0
  linear           :        979        979         84       4994          7          0          6        171         38       2306        671        466          0          0
  logicor          :         35+        47         78       1412          0          0          0         35          2         13          0          0          0          0
  bounddisjunction :          0+        31          0       1318          0          0          0          2          0          4          0          0          0          0
  countsols        :          0          0          0          0          0          0         18          0          0          0          0          0          0          0
Constraint Timings :  TotalTime  SetupTime   Separate  Propagate     EnfoLP     EnfoPS      Check    ResProp    SB-Prop
  integral         :       1.53       0.00       0.00       0.00       1.53       0.00       0.00       0.00       0.00
  varbound         :       0.03       0.00       0.00       0.03       0.00       0.00       0.00       0.00       0.00
  setppc           :       0.03       0.00       0.00       0.03       0.00       0.00       0.00       0.00       0.00
  and              :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
  linear           :       0.12       0.00       0.01       0.10       0.00       0.00       0.00       0.00       0.01
  logicor          :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
  bounddisjunction :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
  countsols        :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
Propagators        : #Propagate   #ResProp    Cutoffs    DomReds
  dualfix          :          2          0          0          0
  genvbounds       :          0          0          0          0
  obbt             :          0          0          0          0
  probing          :          0          0          0          0
  pseudoobj        :        844          0          0         56
  redcost          :        562          0          0          1
  rootredcost      :          1          0          0          0
  vbounds          :       4862          0          0          0
Propagator Timings :  TotalTime  SetupTime   Presolve  Propagate    ResProp    SB-Prop
  dualfix          :       0.00       0.00       0.00       0.00       0.00       0.00
  genvbounds       :       0.00       0.00       0.00       0.00       0.00       0.00
  obbt             :       0.00       0.00       0.00       0.00       0.00       0.00
  probing          :       0.03       0.00       0.03       0.00       0.00       0.00
  pseudoobj        :       0.00       0.00       0.00       0.00       0.00       0.00
  redcost          :       0.02       0.00       0.00       0.01       0.00       0.01
  rootredcost      :       0.00       0.00       0.00       0.00       0.00       0.00
  vbounds          :       0.00       0.00       0.00       0.00       0.00       0.00
Conflict Analysis  :       Time      Calls    Success    DomReds  Conflicts   Literals    Reconvs ReconvLits   LP Iters
  propagation      :       0.00         27         27          -        180       33.3          4       30.5          -
  infeasible LP    :       0.00          6          3          -         33       32.6          0        0.0          0
  bound exceed. LP :       0.00          0          0          -          0        0.0          0        0.0          0
  strong branching :       0.00          0          0          -          0        0.0          0        0.0          0
  pseudo solution  :       0.00          0          0          -          0        0.0          0        0.0          -
  applied globally :       0.00          -          -          0         45       14.7          -          -          -
  applied locally  :          -          -          -          0          0        0.0          -          -          -
Separators         :   ExecTime  SetupTime      Calls    Cutoffs    DomReds       Cuts    Applied      Conss
  cut pool         :       0.00                    83          -          -       1102          -          -    (maximal pool size: 781)
  cgmip            :       0.00       0.00          0          0          0          0          0          0
  clique           :       0.00       0.00         84          0          0         38         36          0
  closecuts        :       0.00       0.00          0          0          0          0          0          0
  cmir             :       0.02       0.00         10          0          0        189        133          0
  disjunctive      :       0.00       0.00          0          0          0          0          0          0
  eccuts           :       0.00       0.00          0          0          0          0          0          0
  flowcover        :       0.13       0.00         15          0          0        301         18          0
  gomory           :       0.05       0.00         10          0          0          4          4          0
  impliedbounds    :       0.02       0.00         84          0          0        701        143          0
  intobj           :       0.00       0.00          0          0          0          0          0          0
  mcf              :       0.00       0.00          1          0          0          0          0          0
  oddcycle         :       0.00       0.00          0          0          0          0          0          0
  rapidlearning    :       0.00       0.00          0          0          0          0          0          0
  strongcg         :       0.04       0.00         20          0          0         11          3          0
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
  relpscost        :       1.53       0.00        232          0          0          0       1423          0          0        180
Primal Heuristics  :   ExecTime  SetupTime      Calls      Found       Best
  LP solutions     :       0.00          -          -          0          -
  pseudo solutions :       0.00          -          -          0          -
  strong branching :       0.00          -          -          0          -
  SpaGreedy        :       0.00       0.00         17          2          1
  actconsdiving    :       0.00       0.00          0          0          0
  bound            :       0.00       0.00          0          0          0
  clique           :       0.00       0.00          0          0          0
  coefdiving       :       0.01       0.00          1          0          0
  crossover        :       0.00       0.00          0          0          0
  dins             :       0.00       0.00          0          0          0
  distributiondivin:       0.00       0.00          1          0          0
  dualval          :       0.00       0.00          0          0          0
  feaspump         :       0.01       0.00          1          0          0
  fixandinfer      :       0.00       0.00          0          0          0
  fracdiving       :       0.01       0.00          1          0          0
  guideddiving     :       0.00       0.00          0          0          0
  indicator        :       0.00       0.00          0          0          0
  intdiving        :       0.00       0.00          0          0          0
  intshifting      :       0.00       0.00          3          0          0
  linesearchdiving :       0.00       0.00          0          0          0
  localbranching   :       0.00       0.00          0          0          0
  mutation         :       0.00       0.00          0          0          0
  nlpdiving        :       0.00       0.00          0          0          0
  objpscostdiving  :       0.01       0.00          1          0          0
  octane           :       0.00       0.00          0          0          0
  ofins            :       0.00       0.00          0          0          0
  oneopt           :       0.00       0.00          0          0          0
  proximity        :       0.00       0.00          0          0          0
  pscostdiving     :       0.00       0.00          1          0          0
  randrounding     :       0.18       0.00         84          0          0
  rens             :       0.37       0.00          1          0          0
  reoptsols        :       0.00       0.00          0          0          0
  rins             :       0.00       0.00          0          0          0
  rootsoldiving    :       0.00       0.00          0          0          0
  rounding         :       0.00       0.00        171          0          0
  shiftandpropagate:       0.00       0.00          0          0          0
  shifting         :       0.07       0.00        100          0          0
  simplerounding   :       0.00       0.00          0          0          0
  subnlp           :       0.00       0.00          0          0          0
  trivial          :       0.00       0.00          2          0          0
  trivialnegation  :       0.00       0.00          0          0          0
  trysol           :       0.00       0.00          0          0          0
  twoopt           :       0.00       0.00          0          0          0
  undercover       :       0.00       0.00          1          0          0
  vbounds          :       0.00       0.00          0          0          0
  veclendiving     :       0.00       0.00          1          0          0
  zeroobj          :       0.00       0.00          0          0          0
  zirounding       :       0.00       0.00         90          0          0
  other solutions  :          -          -          -          0          -
Diving Statistics  :      Calls      Nodes   LP Iters Backtracks   MinDepth   MaxDepth   AvgDepth  NLeafSols  MinSolDpt  MaxSolDpt  AvgSolDpt
  actconsdiving    :          0          -          -          -          -          -          -          -          -          -          -
  coefdiving       :          1          6        214          1         22         22       22.0          -          -          -          -
  distributiondivin:          1          1         10          0         21         21       21.0          -          -          -          -
  fracdiving       :          1          5         64          1         23         23       23.0          -          -          -          -
  guideddiving     :          0          -          -          -          -          -          -          -          -          -          -
  linesearchdiving :          0          -          -          -          -          -          -          -          -          -          -
  pscostdiving     :          1          2         32          0         21         21       21.0          -          -          -          -
  veclendiving     :          1          2         28          0         22         22       22.0          -          -          -          -
LP                 :       Time      Calls Iterations  Iter/call   Iter/sec  Time-0-It Calls-0-It
  primal LP        :       0.00         13          0       0.00          -       0.00         13
  dual LP          :       0.26        397       7689      29.92   29573.08       0.03        140
  lex dual LP      :       0.00          0          0       0.00          -
  barrier LP       :       0.00          0          0       0.00          -       0.00          0
  diving/probing LP:       0.03         38       1199      31.55   39966.67
  strong branching :       1.33        346      42129     121.76   31675.94
    (at root node) :          -         25       1686      67.44          -
  conflict analysis:       0.00          0          0       0.00          -
B&B Tree           :
  number of runs   :          1
  nodes            :        100 (90 internal, 10 leaves)
  nodes (total)    :        100 (90 internal, 10 leaves)
  nodes left       :         81
  max depth        :         20
  max depth (total):         20
  backtracks       :         11 (11.0%)
  delayed cutoffs  :          0
  repropagations   :          2 (1 domain reductions, 0 cutoffs)
  avg switch length:       2.85
  switching time   :       0.00
Root Node          :
  First LP value   : +3.00000000000000e+00
  First LP Iters   :        377
  First LP Time    :       0.00
  Final Dual Bound : +7.03075231743873e-03
  Final Root Iters :       2883
Solution           :
  Solutions found  :          2 (1 improvements)
  First Solution   : +7.87434402999998e-04   (in run 1, after 1 nodes, 0.05 seconds, depth 0, found by <SpaGreedy>)
  Gap First Sol.   :   infinite
  Gap Last Sol.    :   infinite
  Primal Bound     : +7.87434402999998e-04   (in run 1, after 1 nodes, 0.05 seconds, depth 0, found by <SpaGreedy>)
  Dual Bound       : +6.49097573427210e-03
  Gap              :     724.32 %
  Avg. Gap         :      91.17 % (262.56 primal-dual integral)

objective value:                 0.000789625601999987
x_1_1                                               1 	(obj:0)
x_2_2                                               1 	(obj:0)
x_3_4                                               1 	(obj:0)
x_4_1                                               1 	(obj:0)
x_5_1                                               1 	(obj:0)
x_6_2                                               1 	(obj:0)
x_7_2                                               1 	(obj:0)
x_8_4                                               1 	(obj:0)
x_9_1                                               1 	(obj:0)
x_10_1                                              1 	(obj:0)
x_11_4                                              1 	(obj:0)
x_12_4                                              1 	(obj:0)
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
epsI                             0.000789625601999987 	(obj:1)
SCIP Status        : solving was interrupted [node limit reached]
Total Time         :       2.44
  solving          :       2.44
  presolving       :       0.05 (included in solving)
  reading          :       0.00
  copying          :       0.00 (1 #copies) (minimal 0.00, maximal 0.00, average 0.00)
Original Problem   :
  Problem name     : ../instances/Pot4/Pot4_T_20B_alpha=14.txt.spa
  Variables        : 101 (100 binary, 0 integer, 0 implicit integer, 1 continuous)
  Constraints      : 193 initial, 193 maximal
  Objective sense  : maximize
Presolved Problem  :
  Problem name     : t_../instances/Pot4/Pot4_T_20B_alpha=14.txt.spa
  Variables        : 590 (92 binary, 0 integer, 16 implicit integer, 482 continuous)
  Constraints      : 2015 initial, 2055 maximal
Presolvers         :   ExecTime  SetupTime  Calls  FixedVars   AggrVars   ChgTypes  ChgBounds   AddHoles    DelCons    AddCons   ChgSides   ChgCoefs
  boundshift       :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  components       :       0.00       0.00      1          0          0          0          0          0          0          0          0          0
  convertinttobin  :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  domcol           :       0.00       0.00      1          0          0          0          0          0          0          0          0          0
  dualagg          :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  dualinfer        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  gateextraction   :       0.01       0.00      3          0          0          0          0          0        108         34          0          0
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
  setppc           :       0.00       0.00     12          0          0          0          0          0          0          0          0          0
  and              :       0.00       0.00      4          0         12          0          0          0         12          0          0          0
  linear           :       0.00       0.00      9          3          1          0          4          0          7          0         12         12
  logicor          :       0.00       0.00      7          0          0          0          0          0          0          0          0          0
  bounddisjunction :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  quadratic        :       0.01       0.00      2          0          0          0          0          0          0       1915          0          0
  root node        :          -          -      -          1          -          -         57          -          -          -          -          -
Constraints        :     Number  MaxNumber  #Separate #Propagate    #EnfoLP    #EnfoPS     #Check   #ResProp    Cutoffs    DomReds       Cuts    Applied      Conss   Children
  integral         :          0          0          0          0        259          0         33          0          0       1525          0          0          0        196
  varbound         :        869        869         96       5441          1          0          9        433         16       2100          0          0          0          0
  setppc           :        110        110         96       5425          1          0          6        484         19        225         10          1          0          0
  and              :         22         22        158       4244          0          0          4         40          0         91         19          7          0          0
  linear           :        979        979         96       5406          1          0          2        281         49       2605        451        258          0          0
  logicor          :         35+        48         95       1483          0          0          0         28          7         27          1          1          0          0
  bounddisjunction :          0+        44          0       1722          0          0          0          0          0          3          0          0          0          0
  countsols        :          0          0          0          0          0          0          9          0          0          0          0          0          0          0
Constraint Timings :  TotalTime  SetupTime   Separate  Propagate     EnfoLP     EnfoPS      Check    ResProp    SB-Prop
  integral         :       1.20       0.00       0.00       0.00       1.20       0.00       0.00       0.00       0.00
  varbound         :       0.03       0.00       0.01       0.01       0.00       0.00       0.00       0.00       0.01
  setppc           :       0.05       0.00       0.00       0.03       0.00       0.00       0.00       0.00       0.02
  and              :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
  linear           :       0.12       0.00       0.00       0.10       0.00       0.00       0.00       0.00       0.02
  logicor          :       0.03       0.00       0.00       0.02       0.00       0.00       0.00       0.00       0.01
  bounddisjunction :       0.01       0.00       0.00       0.01       0.00       0.00       0.00       0.00       0.00
  countsols        :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
Propagators        : #Propagate   #ResProp    Cutoffs    DomReds
  dualfix          :          2          0          0          0
  genvbounds       :          0          0          0          0
  obbt             :          0          0          0          0
  probing          :          0          0          0          0
  pseudoobj        :        874          0          0         56
  redcost          :        615          0          0          1
  rootredcost      :          1          0          0          0
  vbounds          :       5304          0          0          0
Propagator Timings :  TotalTime  SetupTime   Presolve  Propagate    ResProp    SB-Prop
  dualfix          :       0.00       0.00       0.00       0.00       0.00       0.00
  genvbounds       :       0.00       0.00       0.00       0.00       0.00       0.00
  obbt             :       0.00       0.00       0.00       0.00       0.00       0.00
  probing          :       0.03       0.00       0.03       0.00       0.00       0.00
  pseudoobj        :       0.00       0.00       0.00       0.00       0.00       0.00
  redcost          :       0.00       0.00       0.00       0.00       0.00       0.00
  rootredcost      :       0.00       0.00       0.00       0.00       0.00       0.00
  vbounds          :       0.00       0.00       0.00       0.00       0.00       0.00
Conflict Analysis  :       Time      Calls    Success    DomReds  Conflicts   Literals    Reconvs ReconvLits   LP Iters
  propagation      :       0.00         42         42          -        318       36.2          7       28.6          -
  infeasible LP    :       0.00          7          0          -          0        0.0          0        0.0          0
  bound exceed. LP :       0.00          0          0          -          0        0.0          0        0.0          0
  strong branching :       0.00          0          0          -          0        0.0          0        0.0          0
  pseudo solution  :       0.00          0          0          -          0        0.0          0        0.0          -
  applied globally :       0.00          -          -          0         63       15.7          -          -          -
  applied locally  :          -          -          -          0          1       31.0          -          -          -
Separators         :   ExecTime  SetupTime      Calls    Cutoffs    DomReds       Cuts    Applied      Conss
  cut pool         :       0.00                    95          -          -       1187          -          -    (maximal pool size: 771)
  cgmip            :       0.00       0.00          0          0          0          0          0          0
  clique           :       0.02       0.00         96          0          0         44         35          0
  closecuts        :       0.00       0.00          0          0          0          0          0          0
  cmir             :       0.04       0.00         10          0          0        213        186          0
  disjunctive      :       0.00       0.00          0          0          0          0          0          0
  eccuts           :       0.00       0.00          0          0          0          0          0          0
  flowcover        :       0.10       0.00         15          0          0        217         17          0
  gomory           :       0.02       0.00         10          0          0         10         10          0
  impliedbounds    :       0.00       0.00         96          0          0        727        145          0
  intobj           :       0.00       0.00          0          0          0          0          0          0
  mcf              :       0.00       0.00          1          0          0          0          0          0
  oddcycle         :       0.00       0.00          0          0          0          0          0          0
  rapidlearning    :       0.00       0.00          0          0          0          0          0          0
  strongcg         :       0.05       0.00         20          0          0          4          2          0
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
  relpscost        :       1.20       0.00        258          0          0          0       1525          0          0        196
Primal Heuristics  :   ExecTime  SetupTime      Calls      Found       Best
  LP solutions     :       0.00          -          -          0          -
  pseudo solutions :       0.00          -          -          0          -
  strong branching :       0.00          -          -          0          -
  SpaGreedy        :       0.00       0.00          8          2          1
  actconsdiving    :       0.00       0.00          0          0          0
  bound            :       0.00       0.00          0          0          0
  clique           :       0.00       0.00          0          0          0
  coefdiving       :       0.00       0.00          1          0          0
  crossover        :       0.00       0.00          0          0          0
  dins             :       0.00       0.00          0          0          0
  distributiondivin:       0.01       0.00          2          0          0
  dualval          :       0.00       0.00          0          0          0
  feaspump         :       0.01       0.00          1          0          0
  fixandinfer      :       0.00       0.00          0          0          0
  fracdiving       :       0.01       0.00          2          0          0
  guideddiving     :       0.00       0.00          0          0          0
  indicator        :       0.00       0.00          0          0          0
  intdiving        :       0.00       0.00          0          0          0
  intshifting      :       0.04       0.00          7          0          0
  linesearchdiving :       0.00       0.00          1          0          0
  localbranching   :       0.00       0.00          0          0          0
  mutation         :       0.00       0.00          0          0          0
  nlpdiving        :       0.00       0.00          0          0          0
  objpscostdiving  :       0.00       0.00          1          0          0
  octane           :       0.00       0.00          0          0          0
  ofins            :       0.00       0.00          0          0          0
  oneopt           :       0.00       0.00          0          0          0
  proximity        :       0.00       0.00          0          0          0
  pscostdiving     :       0.01       0.00          2          0          0
  randrounding     :       0.30       0.00         96          0          0
  rens             :       0.26       0.00          1          0          0
  reoptsols        :       0.00       0.00          0          0          0
  rins             :       0.00       0.00          0          0          0
  rootsoldiving    :       0.00       0.00          1          0          0
  rounding         :       0.00       0.00        204          0          0
  shiftandpropagate:       0.00       0.00          0          0          0
  shifting         :       0.06       0.00        110          0          0
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
  zirounding       :       0.00       0.00         97          0          0
  other solutions  :          -          -          -          0          -
Diving Statistics  :      Calls      Nodes   LP Iters Backtracks   MinDepth   MaxDepth   AvgDepth  NLeafSols  MinSolDpt  MaxSolDpt  AvgSolDpt
  actconsdiving    :          0          -          -          -          -          -          -          -          -          -          -
  coefdiving       :          1          4          0          1         21         21       21.0          -          -          -          -
  distributiondivin:          2         11        104          2         21         25       23.0          -          -          -          -
  fracdiving       :          2         12         75          2         20         21       20.5          -          -          -          -
  guideddiving     :          0          -          -          -          -          -          -          -          -          -          -
  linesearchdiving :          1          7         31          1         25         25       25.0          -          -          -          -
  pscostdiving     :          2         17        156          2         20         28       24.0          -          -          -          -
  veclendiving     :          1          4          0          1         21         21       21.0          -          -          -          -
LP                 :       Time      Calls Iterations  Iter/call   Iter/sec  Time-0-It Calls-0-It
  primal LP        :       0.00         16          0       0.00          -       0.00         16
  dual LP          :       0.18        418       9388      29.80   52155.56       0.00        103
  lex dual LP      :       0.00          0          0       0.00          -
  barrier LP       :       0.00          0          0       0.00          -       0.00          0
  diving/probing LP:       0.07         45       2536      56.36   36228.57
  strong branching :       0.97        364      29689      81.56   30607.22
    (at root node) :          -         19       1125      59.21          -
  conflict analysis:       0.00          0          0       0.00          -
B&B Tree           :
  number of runs   :          1
  nodes            :        100 (98 internal, 2 leaves)
  nodes (total)    :        100 (98 internal, 2 leaves)
  nodes left       :         97
  max depth        :         19
  max depth (total):         19
  backtracks       :          7 (7.0%)
  delayed cutoffs  :          0
  repropagations   :          2 (0 domain reductions, 0 cutoffs)
  avg switch length:       2.06
  switching time   :       0.00
Root Node          :
  First LP value   : +3.00000000000000e+00
  First LP Iters   :        357
  First LP Time    :       0.00
  Final Dual Bound : +7.03745560075667e-03
  Final Root Iters :       3658
Solution           :
  Solutions found  :          2 (1 improvements)
  First Solution   : +7.89625601999987e-04   (in run 1, after 1 nodes, 0.05 seconds, depth 0, found by <SpaGreedy>)
  Gap First Sol.   :   infinite
  Gap Last Sol.    :   infinite
  Primal Bound     : +7.89625601999987e-04   (in run 1, after 1 nodes, 0.05 seconds, depth 0, found by <SpaGreedy>)
  Dual Bound       : +6.49952438901749e-03
  Gap              :     723.11 %
  Avg. Gap         :      91.64 % (223.60 primal-dual integral)

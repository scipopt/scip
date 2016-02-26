objective value:                 0.000275616363000039
x_1_2                                               1 	(obj:0)
x_2_1                                               1 	(obj:0)
x_3_4                                               1 	(obj:0)
x_4_1                                               1 	(obj:0)
x_5_4                                               1 	(obj:0)
x_6_2                                               1 	(obj:0)
x_7_1                                               1 	(obj:0)
x_8_4                                               1 	(obj:0)
x_9_2                                               1 	(obj:0)
x_10_1                                              1 	(obj:0)
x_11_2                                              1 	(obj:0)
x_12_4                                              1 	(obj:0)
x_13_4                                              1 	(obj:0)
x_14_4                                              1 	(obj:0)
x_15_1                                              1 	(obj:0)
x_16_2                                              1 	(obj:0)
x_17_1                                              1 	(obj:0)
x_18_4                                              1 	(obj:0)
x_19_2                                              1 	(obj:0)
x_20_1                                              1 	(obj:0)
x_21_4                                              1 	(obj:0)
x_22_2                                              1 	(obj:0)
x_23_2                                              1 	(obj:0)
x_24_2                                              1 	(obj:0)
x_25_1                                              1 	(obj:0)
x_26_4                                              1 	(obj:0)
x_27_1                                              1 	(obj:0)
x_28_4                                              1 	(obj:0)
x_29_1                                              1 	(obj:0)
x_30_4                                              1 	(obj:0)
ind_1                                               1 	(obj:0)
abs_1_4                                             1 	(obj:0)
ind_2                                               1 	(obj:0)
abs_2_1                                             1 	(obj:0)
ind_4                                               1 	(obj:0)
abs_4_2                                             1 	(obj:0)
epsI                             0.000275616363000039 	(obj:1)
SCIP Status        : solving was interrupted [time limit reached]
Total Time         :    3600.00
  solving          :    3600.00
  presolving       :       0.15 (included in solving)
  reading          :       0.00
  copying          :       0.01 (1 #copies) (minimal 0.01, maximal 0.01, average 0.01)
Original Problem   :
  Problem name     : ../instances/Pot3NonCycle/Pot3NonCycle_T_30B_1t.txt.spa
  Variables        : 141 (140 binary, 0 integer, 0 implicit integer, 1 continuous)
  Constraints      : 263 initial, 263 maximal
  Objective sense  : maximize
Presolved Problem  :
  Problem name     : t_../instances/Pot3NonCycle/Pot3NonCycle_T_30B_1t.txt.spa
  Variables        : 910 (132 binary, 0 integer, 16 implicit integer, 762 continuous)
  Constraints      : 3141 initial, 9545 maximal
Presolvers         :   ExecTime  SetupTime  Calls  FixedVars   AggrVars   ChgTypes  ChgBounds   AddHoles    DelCons    AddCons   ChgSides   ChgCoefs
  boundshift       :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  components       :       0.01       0.00      1          0          0          0          0          0          0          0          0          0
  convertinttobin  :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  domcol           :       0.00       0.00      1          0          0          0          0          0          0          0          0          0
  dualagg          :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  dualinfer        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  gateextraction   :       0.00       0.00      3          0          0          0          0          0        108         34          0          0
  implfree         :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  implics          :       0.00       0.00      4          0          0          0          0          0          0          0          0          0
  inttobinary      :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  redvub           :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  stuffing         :       0.01       0.00      1          0          0          0          0          0          0          0          0          0
  trivial          :       0.00       0.00      7          0          0          0          0          0          0          0          0          0
  tworowbnd        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  dualfix          :       0.00       0.00      7          4          0          0          0          0          0          0          0          0
  genvbounds       :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  probing          :       0.09       0.00      1          0          0          0          4          0          0          0          0          0
  pseudoobj        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  varbound         :       0.00       0.00      9          0          0          0          0          0         28         28          0          0
  setppc           :       0.00       0.00     12          0          0          0          0          0          0          0          0          0
  and              :       0.00       0.00      4          0         12          0          0          0         12          0          0          0
  linear           :       0.01       0.00      9          3          1          0          4          0          7          0         12         12
  logicor          :       0.01       0.00      7          0          0          0          0          0          0          0          0          0
  bounddisjunction :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  quadratic        :       0.01       0.00      2          0          0          0          0          0          0       2971          0          0
  root node        :          -          -      -          1          -          -         66          -          -          -          -          -
Constraints        :     Number  MaxNumber  #Separate #Propagate    #EnfoLP    #EnfoPS     #Check   #ResProp    Cutoffs    DomReds       Cuts    Applied      Conss   Children
  integral         :          0          0          0          0     592679          0     221606          0          0      73065          0          0        933    1173800
  varbound         :       1365       1365         79    3337831        176          2      78063    1716411         34   18989832          0          0          0          0
  setppc           :        160        160         79    3337793        176          2         15    1845759        317    1206152         10          0          0          0
  and              :         22         22      17941     108963          0          0         11     195796          6       8014         78         62          0          0
  linear           :       1539       1539         79    3337474        132          0          9     287705        110   32681778      13354      12840          0          0
  logicor          :         55+      1974         78     826353          0          0          0     262210         27      32950          2          2          0          0
  bounddisjunction :          0+      4516          0     986016          0          0          0      19776          6      98049          0          0          0          0
  countsols        :          0          0          0          0          0          1     143541          0          0          0          0          0          0          0
Constraint Timings :  TotalTime  SetupTime   Separate  Propagate     EnfoLP     EnfoPS      Check    ResProp    SB-Prop
  integral         :      42.38       0.00       0.00       0.00      42.34       0.00       0.04       0.00       0.00
  varbound         :      17.17       0.00       0.01      16.26       0.01       0.00       0.31       0.36       0.22
  setppc           :      18.23       0.00       0.00      17.67       0.00       0.00       0.00       0.27       0.29
  and              :       0.09       0.00       0.00       0.04       0.00       0.00       0.00       0.05       0.00
  linear           :      83.22       0.00       0.05      82.08       0.00       0.00       0.00       0.11       0.98
  logicor          :       4.87       0.00       0.00       4.70       0.00       0.00       0.00       0.17       0.00
  bounddisjunction :      17.14       0.00       0.00      17.09       0.00       0.00       0.00       0.00       0.05
  countsols        :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
Propagators        : #Propagate   #ResProp    Cutoffs    DomReds
  dualfix          :          2          0          0          0
  genvbounds       :          0          0          0          0
  obbt             :          0          0          0          0
  probing          :          0          0          0          0
  pseudoobj        :    3522720          0          0         65
  redcost          :    1263551          0          0     128884
  rootredcost      :          1          0          0          0
  vbounds          :    3929965       6175          0        228
Propagator Timings :  TotalTime  SetupTime   Presolve  Propagate    ResProp    SB-Prop
  dualfix          :       0.00       0.00       0.00       0.00       0.00       0.00
  genvbounds       :       0.44       0.00       0.00       0.44       0.00       0.00
  obbt             :       0.00       0.00       0.00       0.00       0.00       0.00
  probing          :       0.09       0.00       0.09       0.00       0.00       0.00
  pseudoobj        :       1.39       0.00       0.00       1.39       0.00       0.00
  redcost          :       6.33       0.00       0.00       6.33       0.00       0.00
  rootredcost      :       0.77       0.00       0.00       0.77       0.00       0.00
  vbounds          :       1.46       0.00       0.00       1.43       0.01       0.02
Conflict Analysis  :       Time      Calls    Success    DomReds  Conflicts   Literals    Reconvs ReconvLits   LP Iters
  propagation      :       0.03        390        387          -       1822       46.0         11       46.8          -
  infeasible LP    :      16.27      99545      76613          -    1032232       39.3        186        9.5          0
  bound exceed. LP :       0.00          0          0          -          0        0.0          0        0.0          0
  strong branching :       0.00          0          0          -          0        0.0          0        0.0          0
  pseudo solution  :       0.00          0          0          -          0        0.0          0        0.0          -
  applied globally :       1.84          -          -          0     142339       20.0          -          -          -
  applied locally  :          -          -          -          0       3569       47.6          -          -          -
Separators         :   ExecTime  SetupTime      Calls    Cutoffs    DomReds       Cuts    Applied      Conss
  cut pool         :       0.00                    77          -          -       1535          -          -    (maximal pool size: 930)
  cgmip            :       0.00       0.00          0          0          0          0          0          0
  clique           :       0.04       0.00         78          0          0         67         36          0
  closecuts        :       0.00       0.00          0          0          0          0          0          0
  cmir             :       0.04       0.00         10          0          0        341        225          0
  disjunctive      :       0.00       0.00          0          0          0          0          0          0
  eccuts           :       0.00       0.00          0          0          0          0          0          0
  flowcover        :       0.57       0.00         15          0          0        281         37          0
  gomory           :       0.13       0.00         10          0          0          8          8          0
  impliedbounds    :       0.00       0.00         78          0          0        736        137          0
  intobj           :       0.00       0.00          0          0          0          0          0          0
  mcf              :       0.00       0.00          1          0          0          0          0          0
  oddcycle         :       0.00       0.00          0          0          0          0          0          0
  rapidlearning    :       0.00       0.00          0          0          0          0          0          0
  strongcg         :       0.19       0.00         20          0          0         11          3          0
  zerohalf         :       0.00       0.00          0          0          0          0          0          0
Pricers            :   ExecTime  SetupTime      Calls       Vars
  problem variables:       0.00          -          0          0
Branching Rules    :   ExecTime  SetupTime   BranchLP  BranchExt   BranchPS    Cutoffs    DomReds       Cuts      Conss   Children
  allfullstrong    :       0.00       0.00          0          0          0          0          0          0          0          0
  cloud            :       0.00       0.00          0          0          0          0          0          0          0          0
  distribution     :       0.00       0.00          0          0          0          0          0          0          0          0
  fullstrong       :       0.00       0.00          0          0          0          0          0          0          0          0
  inference        :       0.00       0.00          0          0          1          0          0          0          0          2
  leastinf         :       0.00       0.00          0          0          0          0          0          0          0          0
  mostinf          :       0.00       0.00          0          0          0          0          0          0          0          0
  multaggr         :       0.00       0.00          0          0          0          0          0          0          0          0
  nodereopt        :       0.00       0.00          0          0          0          0          0          0          0          0
  pscost           :       0.00       0.00          0          0          0          0          0          0          0          0
  random           :       0.00       0.00          0          0          0          0          0          0          0          0
  relpscost        :      42.03       0.00     592503          0          0          0      73065          0        407    1173800
Primal Heuristics  :   ExecTime  SetupTime      Calls      Found       Best
  LP solutions     :       0.00          -          -          0          -
  pseudo solutions :       0.00          -          -          0          -
  strong branching :       0.00          -          -          0          -
  SpaGreedy        :     187.40       0.00     143540          2          1
  actconsdiving    :       0.00       0.00          0          0          0
  bound            :       0.00       0.00          0          0          0
  clique           :       0.00       0.00          0          0          0
  coefdiving       :      41.34       0.00       2764          0          0
  crossover        :       0.24       0.00          0          0          0
  dins             :       0.00       0.00          0          0          0
  distributiondivin:      69.19       0.00       2765          0          0
  dualval          :       0.00       0.00          0          0          0
  feaspump         :      26.74       0.00        430          0          0
  fixandinfer      :       0.00       0.00          0          0          0
  fracdiving       :      50.33       0.00       2765          0          0
  guideddiving     :       0.00       0.00          0          0          0
  indicator        :       0.07       0.00          0          0          0
  intdiving        :       0.00       0.00          0          0          0
  intshifting      :       0.54       0.00         53          0          0
  linesearchdiving :      56.13       0.00       2764          0          0
  localbranching   :       0.00       0.00          0          0          0
  mutation         :       0.00       0.00          0          0          0
  nlpdiving        :       0.01       0.00          0          0          0
  objpscostdiving  :      25.82       0.00        346          0          0
  octane           :       0.00       0.00          0          0          0
  ofins            :       0.00       0.00          0          0          0
  oneopt           :       0.24       0.00          0          0          0
  proximity        :       0.00       0.00          0          0          0
  pscostdiving     :      46.14       0.00       2764          0          0
  randrounding     :       1.14       0.00      78125          0          0
  rens             :       0.45       0.00          1          0          0
  reoptsols        :       0.00       0.00          0          0          0
  rins             :       0.33       0.00          0          0          0
  rootsoldiving    :      25.82       0.00        744          0          0
  rounding         :       0.46       0.00      10753          0          0
  shiftandpropagate:       0.00       0.00          0          0          0
  shifting         :       1.92       0.00       4070          0          0
  simplerounding   :       0.05       0.00          0          0          0
  subnlp           :       0.21       0.00          0          0          0
  trivial          :       0.00       0.00          2          0          0
  trivialnegation  :       0.00       0.00          0          0          0
  trysol           :       0.15       0.00          0          0          0
  twoopt           :       0.00       0.00          0          0          0
  undercover       :       0.00       0.00          1          0          0
  vbounds          :       0.00       0.00          0          0          0
  veclendiving     :      46.45       0.00       2764          0          0
  zeroobj          :       0.00       0.00          0          0          0
  zirounding       :       0.21       0.00       1000          0          0
  other solutions  :          -          -          -          0          -
Diving Statistics  :      Calls      Nodes   LP Iters Backtracks   MinDepth   MaxDepth   AvgDepth  NLeafSols  MinSolDpt  MaxSolDpt  AvgSolDpt
  actconsdiving    :          0          -          -          -          -          -          -          -          -          -          -
  coefdiving       :       2764      25314     502322       3914         17         55       29.9          -          -          -          -
  distributiondivin:       2765      19193     605770       3373         13         45       27.9          -          -          -          -
  fracdiving       :       2765      25185     641466       3958         15         49       29.9          -          -          -          -
  guideddiving     :          0          -          -          -          -          -          -          -          -          -          -
  linesearchdiving :       2764      21241     707116       4315         15         43       28.3          -          -          -          -
  pscostdiving     :       2764      30189     581934       3357         17         51       31.9          -          -          -          -
  veclendiving     :       2764      21734     565894       4841         17         41       28.3          -          -          -          -
LP                 :       Time      Calls Iterations  Iter/call   Iter/sec  Time-0-It Calls-0-It
  primal LP        :       0.16         18       4181     696.83   26131.25       0.01         12
  dual LP          :    2389.67    1354195   36379554      30.23   15223.67      38.68     150688
  lex dual LP      :       0.00          0          0       0.00          -
  barrier LP       :       0.00          0          0       0.00          -       0.00          0
  diving/probing LP:     325.40      75618    4786814      63.30   14710.55
  strong branching :      27.35       8752     418223      47.79   15291.52
    (at root node) :          -         23       1959      85.17          -
  conflict analysis:       0.00          0          0       0.00          -
B&B Tree           :
  number of runs   :          1
  nodes            :    1098637 (586901 internal, 511736 leaves)
  nodes (total)    :    1098637 (586901 internal, 511736 leaves)
  nodes left       :      75151
  max depth        :         42
  max depth (total):         42
  backtracks       :     272829 (24.8%)
  delayed cutoffs  :         15
  repropagations   :      41899 (457108 domain reductions, 8 cutoffs)
  avg switch length:       5.39
  switching time   :     219.03
Root Node          :
  First LP value   : +3.00000000000000e+00
  First LP Iters   :        581 (29050.00 Iter/sec)
  First LP Time    :       0.02
  Final Dual Bound : +3.59814281799942e-03
  Final Root Iters :       4973
Solution           :
  Solutions found  :          2 (1 improvements)
  First Solution   : +2.75616363000039e-04   (in run 1, after 1 nodes, 0.15 seconds, depth 0, found by <SpaGreedy>)
  Gap First Sol.   :   infinite
  Gap Last Sol.    :   infinite
  Primal Bound     : +2.75616363000039e-04   (in run 1, after 1 nodes, 0.15 seconds, depth 0, found by <SpaGreedy>)
  Dual Bound       : +3.51431817724514e-04
  Gap              :      27.51 %
  Avg. Gap         :      92.34 % (332439.81 primal-dual integral)

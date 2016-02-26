objective value:                  0.00103971965299997
x_1_1                                               1 	(obj:0)
x_2_2                                               1 	(obj:0)
x_3_4                                               1 	(obj:0)
x_4_2                                               1 	(obj:0)
x_5_4                                               1 	(obj:0)
x_6_4                                               1 	(obj:0)
x_7_1                                               1 	(obj:0)
x_8_2                                               1 	(obj:0)
x_9_4                                               1 	(obj:0)
x_10_1                                              1 	(obj:0)
x_11_2                                              1 	(obj:0)
x_12_4                                              1 	(obj:0)
x_13_1                                              1 	(obj:0)
x_14_1                                              1 	(obj:0)
x_15_2                                              1 	(obj:0)
x_16_4                                              1 	(obj:0)
x_17_2                                              1 	(obj:0)
x_18_2                                              1 	(obj:0)
x_19_2                                              1 	(obj:0)
x_20_1                                              1 	(obj:0)
x_21_4                                              1 	(obj:0)
x_22_1                                              1 	(obj:0)
x_23_1                                              1 	(obj:0)
x_24_4                                              1 	(obj:0)
x_25_4                                              1 	(obj:0)
x_26_4                                              1 	(obj:0)
x_27_4                                              1 	(obj:0)
x_28_4                                              1 	(obj:0)
x_29_2                                              1 	(obj:0)
x_30_1                                              1 	(obj:0)
ind_1                                               1 	(obj:0)
abs_1_4                                             1 	(obj:0)
ind_2                                               1 	(obj:0)
abs_2_1                                             1 	(obj:0)
ind_4                                               1 	(obj:0)
abs_4_2                                             1 	(obj:0)
epsI                              0.00103971965299997 	(obj:1)
SCIP Status        : solving was interrupted [time limit reached]
Total Time         :    3600.01
  solving          :    3600.00
  presolving       :       0.27 (included in solving)
  reading          :       0.01
  copying          :       0.00 (0 times copied the problem)
Original Problem   :
  Problem name     : ../instances/Pot6Cycle/Pot6Cycle_T_30B_1t.txt.spa
  Variables        : 141 (140 binary, 0 integer, 0 implicit integer, 1 continuous)
  Constraints      : 263 initial, 263 maximal
  Objective sense  : maximize
Presolved Problem  :
  Problem name     : t_../instances/Pot6Cycle/Pot6Cycle_T_30B_1t.txt.spa
  Variables        : 910 (132 binary, 0 integer, 16 implicit integer, 762 continuous)
  Constraints      : 3189 initial, 4264 maximal
Presolvers         :   ExecTime  SetupTime  Calls  FixedVars   AggrVars   ChgTypes  ChgBounds   AddHoles    DelCons    AddCons   ChgSides   ChgCoefs
  boundshift       :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  components       :       0.01       0.00      1          0          0          0          0          0          0          0          0          0
  convertinttobin  :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  domcol           :       0.01       0.00      1          0          0          0          0          0          0          0          0          0
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
  probing          :       0.16       0.00      1          0          0          0          0          0          0          0          0          0
  pseudoobj        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  varbound         :       0.01       0.00      9          0          0          0          0          0         28         28          0          0
  setppc           :       0.00       0.00     12          0          0          0          0          0          0          0          0          0
  and              :       0.00       0.00      4          0         12          0          0          0         12          0          0          0
  linear           :       0.04       0.00      9          3          1          0          4          0          7          0         12         12
  logicor          :       0.01       0.00      7          0          0          0          0          0          0          0          0          0
  bounddisjunction :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  quadratic        :       0.03       0.00      2          0          0          0          0          0          0       3019          0          0
  root node        :          -          -      -          1          -          -         68          -          -          -          -          -
Constraints        :     Number  MaxNumber  #Separate #Propagate    #EnfoLP    #EnfoPS     #Check   #ResProp    Cutoffs    DomReds       Cuts    Applied      Conss   Children
  integral         :          0          0          0          0    1761261          0     519773          0          0      46051          0          0          0    2831250
  varbound         :       1413       1413         86    8711528     343231          0      60724    1031560      29564   33958767          0          0          0          0
  setppc           :        160        160         86    8681864     343231          0      14529    2128344      11977    2065639         10          0          0          0
  and              :         22         22       1109      20743          0          0      14527      72729          2        920         24          9          0          0
  linear           :       1539       1539         86    8669985     343231          0      14472    1169402      37919   57814806   21057416   21057080          0          0
  logicor          :         55+       755         83    1415975          0          0          0     376502        224      17093          0          0          0          0
  bounddisjunction :          0+       437          0    2260497          0          0          0      59880       1014     819238          0          0          0          0
  countsols        :          0          0          0          0          0          0     459043          0          0          0          0          0          0          0
Constraint Timings :  TotalTime  SetupTime   Separate  Propagate     EnfoLP     EnfoPS      Check    ResProp    SB-Prop
  integral         :      21.02       0.00       0.00       0.00      20.98       0.00       0.04       0.00       0.00
  varbound         :      41.72       0.00       0.00      37.47       2.83       0.00       1.14       0.17       0.11
  setppc           :      25.66       0.00       0.00      24.89       0.08       0.00       0.00       0.49       0.20
  and              :       0.05       0.00       0.01       0.00       0.00       0.00       0.01       0.03       0.00
  linear           :     164.45       0.00       0.04     143.71      18.40       0.00       0.50       1.19       0.61
  logicor          :       3.41       0.00       0.00       3.14       0.00       0.00       0.00       0.26       0.01
  bounddisjunction :       4.62       0.00       0.00       4.61       0.00       0.00       0.00       0.00       0.01
  countsols        :       0.01       0.00       0.00       0.00       0.00       0.00       0.01       0.00       0.00
Propagators        : #Propagate   #ResProp    Cutoffs    DomReds
  dualfix          :          2          0          0          0
  genvbounds       :          0          0          0          0
  obbt             :          0          0          0          0
  probing          :          0          0          0          0
  pseudoobj        :    8931412          0         15         67
  redcost          :    3813595          0          0     402422
  rootredcost      :          1          0          0          0
  vbounds          :   10470483          0          0        107
Propagator Timings :  TotalTime  SetupTime   Presolve  Propagate    ResProp    SB-Prop
  dualfix          :       0.00       0.00       0.00       0.00       0.00       0.00
  genvbounds       :       1.73       0.00       0.00       1.72       0.00       0.01
  obbt             :       0.00       0.00       0.00       0.00       0.00       0.00
  probing          :       0.16       0.00       0.16       0.00       0.00       0.00
  pseudoobj        :       3.92       0.00       0.00       3.92       0.00       0.00
  redcost          :      19.70       0.00       0.00      19.70       0.00       0.00
  rootredcost      :       2.01       0.00       0.00       1.99       0.00       0.02
  vbounds          :       3.50       0.00       0.00       3.50       0.00       0.00
Conflict Analysis  :       Time      Calls    Success    DomReds  Conflicts   Literals    Reconvs ReconvLits   LP Iters
  propagation      :       7.61      42816      42526          -     979837       66.7       2979       80.6          -
  infeasible LP    :       2.32      19793       3928          -      71967       54.3          7       13.6          0
  bound exceed. LP :       0.00          0          0          -          0        0.0          0        0.0          0
  strong branching :       0.00          0          0          -          0        0.0          0        0.0          0
  pseudo solution  :       0.00          0          0          -          0        0.0          0        0.0          -
  applied globally :       0.44          -          -          0      27346       26.1          -          -          -
  applied locally  :          -          -          -          0      21102       32.1          -          -          -
Separators         :   ExecTime  SetupTime      Calls    Cutoffs    DomReds       Cuts    Applied      Conss
  cut pool         :       0.00                    84          -          -        973          -          -    (maximal pool size: 620)
  cgmip            :       0.00       0.00          0          0          0          0          0          0
  clique           :       0.02       0.00         85          0          0         50         46          0
  closecuts        :       0.00       0.00          0          0          0          0          0          0
  cmir             :       0.03       0.00         10          0          0         62         11          0
  disjunctive      :       0.00       0.00          0          0          0          0          0          0
  eccuts           :       0.00       0.00          0          0          0          0          0          0
  flowcover        :       0.53       0.00         15          0          0        250         48          0
  gomory           :       0.15       0.00         10          0          0         23         18          0
  impliedbounds    :       0.01       0.00         85          0          0        649        156          0
  intobj           :       0.00       0.00          0          0          0          0          0          0
  mcf              :       0.00       0.00          1          0          0          0          0          0
  oddcycle         :       0.00       0.00          0          0          0          0          0          0
  rapidlearning    :       0.00       0.00          0          0          0          0          0          0
  strongcg         :       0.24       0.00         20          0          0          3          2          0
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
  relpscost        :      20.21       0.00    1418030          0          0          0      46051          0        -84    2831250
Primal Heuristics  :   ExecTime  SetupTime      Calls      Found       Best
  LP solutions     :       0.00          -          -          0          -
  pseudo solutions :       0.00          -          -          0          -
  strong branching :       0.00          -          -          0          -
  SpaGreedy        :     591.87       0.00     459042          2          1
  actconsdiving    :       0.00       0.00          0          0          0
  bound            :       0.00       0.00          0          0          0
  clique           :       0.00       0.00          0          0          0
  coefdiving       :      45.30       0.00       8981          0          0
  crossover        :       0.29       0.00          0          0          0
  dins             :       0.00       0.00          0          0          0
  distributiondivin:      87.91       0.00       8981          0          0
  dualval          :       0.00       0.00          0          0          0
  feaspump         :      14.99       0.00       4490          0          0
  fixandinfer      :       0.00       0.00          0          0          0
  fracdiving       :      49.08       0.00       8981          0          0
  guideddiving     :       0.06       0.00          0          0          0
  indicator        :       0.27       0.00          0          0          0
  intdiving        :       0.00       0.00          0          0          0
  intshifting      :       0.39       0.00         94          0          0
  linesearchdiving :      48.37       0.00       8981          0          0
  localbranching   :       0.00       0.00          0          0          0
  mutation         :       0.00       0.00          0          0          0
  nlpdiving        :       0.00       0.00          0          0          0
  objpscostdiving  :      23.36       0.00       3289          0          0
  octane           :       0.00       0.00          0          0          0
  ofins            :       0.00       0.00          0          0          0
  oneopt           :       0.61       0.00          0          0          0
  proximity        :       0.00       0.00          0          0          0
  pscostdiving     :      52.85       0.00       8981          0          0
  randrounding     :       1.25       0.00      43424          0          0
  rens             :       0.00       0.00          0          0          0
  reoptsols        :       0.00       0.00          0          0          0
  rins             :       0.43       0.00          0          0          0
  rootsoldiving    :      22.16       0.00       4491          0          0
  rounding         :       0.87       0.00      16692          0          0
  shiftandpropagate:       0.00       0.00          0          0          0
  shifting         :       1.90       0.00       6065          0          0
  simplerounding   :       0.18       0.00          0          0          0
  subnlp           :       0.53       0.00          0          0          0
  trivial          :       0.00       0.00          2          0          0
  trivialnegation  :       0.00       0.00          0          0          0
  trysol           :       0.54       0.00          0          0          0
  twoopt           :       0.00       0.00          0          0          0
  undercover       :       0.00       0.00          1          0          0
  vbounds          :       0.00       0.00          0          0          0
  veclendiving     :      46.31       0.00       8981          0          0
  zeroobj          :       0.00       0.00          0          0          0
  zirounding       :       0.38       0.00       1000          0          0
  other solutions  :          -          -          -          0          -
Diving Statistics  :      Calls      Nodes   LP Iters Backtracks   MinDepth   MaxDepth   AvgDepth  NLeafSols  MinSolDpt  MaxSolDpt  AvgSolDpt
  actconsdiving    :          0          -          -          -          -          -          -          -          -          -          -
  coefdiving       :       8981      92101    1181326       7531         20         70       32.1          -          -          -          -
  distributiondivin:       8981      84025    1323824       8447         21         54       31.2          -          -          -          -
  fracdiving       :       8981      89789    1349994       7836         25         71       31.8          -          -          -          -
  guideddiving     :          0          -          -          -          -          -          -          -          -          -          -
  linesearchdiving :       8981      89586    1293740       7706         21         60       31.8          -          -          -          -
  pscostdiving     :       8981      88428    1288307       8252         21         76       31.7          -          -          -          -
  veclendiving     :       8981      91914    1237155       8460         26         43       32.0          -          -          -          -
LP                 :       Time      Calls Iterations  Iter/call   Iter/sec  Time-0-It Calls-0-It
  primal LP        :       0.00         13          0       0.00          -       0.00         13
  dual LP          :    1432.22    3191574   38619860      13.83   26965.03      75.27     399882
  lex dual LP      :       0.00          0          0       0.00          -
  barrier LP       :       0.00          0          0       0.00          -       0.00          0
  diving/probing LP:     238.26     247068    8734594      35.35   36659.93
  strong branching :      10.20       3644     191862      52.65   18810.00
    (at root node) :          -         21       1971      93.86          -
  conflict analysis:       0.00          0          0       0.00          -
B&B Tree           :
  number of runs   :          1
  nodes            :    2569646 (1415625 internal, 1154021 leaves)
  nodes (total)    :    2569646 (1415625 internal, 1154021 leaves)
  nodes left       :     260016
  max depth        :         35
  max depth (total):         35
  backtracks       :     671939 (26.1%)
  delayed cutoffs  :       1590
  repropagations   :      25799 (802391 domain reductions, 1573 cutoffs)
  avg switch length:       4.34
  switching time   :     358.17
Root Node          :
  First LP value   : +3.00000000000000e+00
  First LP Iters   :        617 (61700.00 Iter/sec)
  First LP Time    :       0.01
  Final Dual Bound : +1.09496686887221e-02
  Final Root Iters :       5268
Solution           :
  Solutions found  :          2 (1 improvements)
  First Solution   : +1.03971965299997e-03   (in run 1, after 1 nodes, 0.27 seconds, depth 0, found by <SpaGreedy>)
  Gap First Sol.   :   infinite
  Gap Last Sol.    :   infinite
  Primal Bound     : +1.03971965299997e-03   (in run 1, after 1 nodes, 0.27 seconds, depth 0, found by <SpaGreedy>)
  Dual Bound       : +2.92464550895063e-03
  Gap              :     181.29 %
  Avg. Gap         :      90.51 % (325837.89 primal-dual integral)

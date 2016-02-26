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
  presolving       :       0.29 (included in solving)
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
  Constraints      : 3189 initial, 4615 maximal
Presolvers         :   ExecTime  SetupTime  Calls  FixedVars   AggrVars   ChgTypes  ChgBounds   AddHoles    DelCons    AddCons   ChgSides   ChgCoefs
  boundshift       :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  components       :       0.00       0.00      1          0          0          0          0          0          0          0          0          0
  convertinttobin  :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  domcol           :       0.01       0.00      1          0          0          0          0          0          0          0          0          0
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
  probing          :       0.17       0.00      1          0          0          0          0          0          0          0          0          0
  pseudoobj        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  varbound         :       0.00       0.00      9          0          0          0          0          0         28         28          0          0
  setppc           :       0.02       0.00     12          0          0          0          0          0          0          0          0          0
  and              :       0.00       0.00      4          0         12          0          0          0         12          0          0          0
  linear           :       0.04       0.00      9          3          1          0          4          0          7          0         12         12
  logicor          :       0.00       0.00      7          0          0          0          0          0          0          0          0          0
  bounddisjunction :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  quadratic        :       0.04       0.00      2          0          0          0          0          0          0       3019          0          0
  root node        :          -          -      -          1          -          -         66          -          -          -          -          -
Constraints        :     Number  MaxNumber  #Separate #Propagate    #EnfoLP    #EnfoPS     #Check   #ResProp    Cutoffs    DomReds       Cuts    Applied      Conss   Children
  integral         :          0          0          0          0    1698959          0     250321          0          0      42852          0          0          0    3050176
  varbound         :       1413       1413         79   10382020     171986          0      58360    2364631      33376   53509916          0          0          0          0
  setppc           :        160+       265         79   10348632     171986          0      10265    3064093      23474    3320982         10          0          0          0
  and              :         22         22       1752      31204          0          0      10262     161770         13       1506         23          8          0          0
  linear           :       1539       1539         79   10325157     136508          0      10208    2037006      53791  105943440    8375381    8375044          0          0
  logicor          :         55+       878         56    2846674          0          0          0    1061699       9900     553724          0          0          0          0
  bounddisjunction :          0+       522          0    2865032          0          0          0      89430       5077    1037204          0          0          0          0
  countsols        :          0          0          0          0          0          0     191944          0          0          0          0          0          0          0
Constraint Timings :  TotalTime  SetupTime   Separate  Propagate     EnfoLP     EnfoPS      Check    ResProp    SB-Prop
  integral         :      27.28       0.00       0.00       0.00      27.27       0.00       0.01       0.00       0.00
  varbound         :      55.03       0.00       0.02      52.06       1.63       0.00       0.86       0.39       0.07
  setppc           :      50.68       0.00       0.00      49.55       0.19       0.00       0.03       0.63       0.28
  and              :       0.05       0.00       0.00       0.01       0.00       0.00       0.01       0.03       0.00
  linear           :     256.12       0.00       0.03     244.96       8.13       0.00       0.39       2.03       0.58
  logicor          :      13.65       0.00       0.00      13.01       0.00       0.00       0.00       0.63       0.01
  bounddisjunction :       7.05       0.00       0.00       7.03       0.00       0.00       0.00       0.02       0.00
  countsols        :       0.01       0.00       0.00       0.00       0.00       0.00       0.01       0.00       0.00
Propagators        : #Propagate   #ResProp    Cutoffs    DomReds
  dualfix          :          2          0          0          0
  genvbounds       :          0          0          0          0
  obbt             :          0          0          0          0
  probing          :          0          0          0          0
  pseudoobj        :   10367970          0         49         65
  redcost          :    3713907          0        253     511570
  rootredcost      :          1          0          0          0
  vbounds          :   12079774         53          0        559
Propagator Timings :  TotalTime  SetupTime   Presolve  Propagate    ResProp    SB-Prop
  dualfix          :       0.00       0.00       0.00       0.00       0.00       0.00
  genvbounds       :       2.06       0.00       0.00       2.04       0.00       0.02
  obbt             :       0.00       0.00       0.00       0.00       0.00       0.00
  probing          :       0.17       0.00       0.17       0.00       0.00       0.00
  pseudoobj        :       5.08       0.00       0.00       5.07       0.00       0.01
  redcost          :      25.58       0.00       0.00      25.58       0.00       0.00
  rootredcost      :       2.74       0.00       0.00       2.74       0.00       0.00
  vbounds          :       3.83       0.00       0.00       3.83       0.00       0.00
Conflict Analysis  :       Time      Calls    Success    DomReds  Conflicts   Literals    Reconvs ReconvLits   LP Iters
  propagation      :      12.48      72031      66333          -    1458070       64.8       3080       71.9          -
  infeasible LP    :       2.60      20896       7263          -     120977       49.9          0        0.0          0
  bound exceed. LP :       0.00          0          0          -          0        0.0          0        0.0          0
  strong branching :       0.00          0          0          -          0        0.0          0        0.0          0
  pseudo solution  :       0.00          0          0          -          0        0.0          0        0.0          -
  applied globally :       0.61          -          -          0      35245       25.3          -          -          -
  applied locally  :          -          -          -          0      36060       38.1          -          -          -
Separators         :   ExecTime  SetupTime      Calls    Cutoffs    DomReds       Cuts    Applied      Conss
  cut pool         :       0.00                    77          -          -        908          -          -    (maximal pool size: 612)
  cgmip            :       0.00       0.00          0          0          0          0          0          0
  clique           :       0.01       0.00         78          0          0         49         46          0
  closecuts        :       0.00       0.00          0          0          0          0          0          0
  cmir             :       0.02       0.00         10          0          0         62         12          0
  disjunctive      :       0.00       0.00          0          0          0          0          0          0
  eccuts           :       0.00       0.00          0          0          0          0          0          0
  flowcover        :       0.50       0.00         15          0          0        250         48          0
  gomory           :       0.14       0.00         10          0          0         23         18          0
  impliedbounds    :       0.01       0.00         78          0          0        617        148          0
  intobj           :       0.00       0.00          0          0          0          0          0          0
  mcf              :       0.00       0.00          1          0          0          0          0          0
  oddcycle         :       0.00       0.00          0          0          0          0          0          0
  rapidlearning    :       0.00       0.00          0          0          0          0          0          0
  strongcg         :       0.20       0.00         20          0          0          3          2          0
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
  relpscost        :      26.61       0.00    1526973          0          0          0      42852          0       -200    3050176
Primal Heuristics  :   ExecTime  SetupTime      Calls      Found       Best
  LP solutions     :       0.00          -          -          0          -
  pseudo solutions :       0.00          -          -          0          -
  strong branching :       0.00          -          -          0          -
  SpaGreedy        :     251.06       0.00     191943          2          1
  actconsdiving    :       0.00       0.00          0          0          0
  bound            :       0.00       0.00          0          0          0
  clique           :       0.00       0.00          0          0          0
  coefdiving       :      53.97       0.00       9805          0          0
  crossover        :       0.43       0.00          0          0          0
  dins             :       0.00       0.00          0          0          0
  distributiondivin:      97.62       0.00       9804          0          0
  dualval          :       0.00       0.00          0          0          0
  feaspump         :      19.26       0.00       4412          0          0
  fixandinfer      :       0.00       0.00          0          0          0
  fracdiving       :      56.72       0.00       9804          0          0
  guideddiving     :       0.09       0.00          0          0          0
  indicator        :       0.21       0.00          0          0          0
  intdiving        :       0.00       0.00          0          0          0
  intshifting      :       0.49       0.00        124          0          0
  linesearchdiving :      57.58       0.00       9805          0          0
  localbranching   :       0.00       0.00          0          0          0
  mutation         :       0.00       0.00          0          0          0
  nlpdiving        :       0.02       0.00          0          0          0
  objpscostdiving  :      22.85       0.00       2671          0          0
  octane           :       0.00       0.00          0          0          0
  ofins            :       0.00       0.00          0          0          0
  oneopt           :       0.65       0.00          0          0          0
  proximity        :       0.00       0.00          0          0          0
  pscostdiving     :      60.47       0.00       9805          0          0
  randrounding     :       1.18       0.00      47414          0          0
  rens             :       0.00       0.00          0          0          0
  reoptsols        :       0.00       0.00          0          0          0
  rins             :       0.64       0.00          0          0          0
  rootsoldiving    :      28.31       0.00       4903          0          0
  rounding         :       0.71       0.00      17310          0          0
  shiftandpropagate:       0.00       0.00          0          0          0
  shifting         :       1.98       0.00       5323          0          0
  simplerounding   :       0.27       0.00          0          0          0
  subnlp           :       0.50       0.00          0          0          0
  trivial          :       0.00       0.00          2          0          0
  trivialnegation  :       0.00       0.00          0          0          0
  trysol           :       0.54       0.00          0          0          0
  twoopt           :       0.00       0.00          0          0          0
  undercover       :       0.00       0.00          1          0          0
  vbounds          :       0.00       0.00          0          0          0
  veclendiving     :      55.67       0.00       9804          0          0
  zeroobj          :       0.00       0.00          0          0          0
  zirounding       :       0.60       0.00       1000          0          0
  other solutions  :          -          -          -          0          -
Diving Statistics  :      Calls      Nodes   LP Iters Backtracks   MinDepth   MaxDepth   AvgDepth  NLeafSols  MinSolDpt  MaxSolDpt  AvgSolDpt
  actconsdiving    :          0          -          -          -          -          -          -          -          -          -          -
  coefdiving       :       9805     101132    1273014       9066         25         75       31.3          -          -          -          -
  distributiondivin:       9804      92536    1425970       9952         21         55       30.5          -          -          -          -
  fracdiving       :       9804      96869    1448217       9664         23         78       30.9          -          -          -          -
  guideddiving     :          0          -          -          -          -          -          -          -          -          -          -
  linesearchdiving :       9805      96247    1411463       9155         26         58       30.8          -          -          -          -
  pscostdiving     :       9805      95676    1388175       9511         21         78       30.8          -          -          -          -
  veclendiving     :       9804      98997    1327373       9674         26         39       31.1          -          -          -          -
LP                 :       Time      Calls Iterations  Iter/call   Iter/sec  Time-0-It Calls-0-It
  primal LP        :       0.00         12          0       0.00          -       0.00         12
  dual LP          :    1461.99    3257159   38599126      14.45   26401.77     108.70     585692
  lex dual LP      :       0.00          0          0       0.00          -
  barrier LP       :       0.00          0          0       0.00          -       0.00          0
  diving/probing LP:     270.62     280759    9558829      34.05   35321.96
  strong branching :      14.96       4219     248079      58.80   16582.82
    (at root node) :          -         21       1736      82.67          -
  conflict analysis:       0.00          0          0       0.00          -
B&B Tree           :
  number of runs   :          1
  nodes            :    2767774 (1525088 internal, 1242686 leaves)
  nodes (total)    :    2767774 (1525088 internal, 1242686 leaves)
  nodes left       :     275649
  max depth        :         34
  max depth (total):         34
  backtracks       :     720467 (26.0%)
  delayed cutoffs  :       6754
  repropagations   :      35632 (2570679 domain reductions, 6694 cutoffs)
  avg switch length:       4.47
  switching time   :     485.47
Root Node          :
  First LP value   : +3.00000000000000e+00
  First LP Iters   :        617
  First LP Time    :       0.00
  Final Dual Bound : +1.09496686887221e-02
  Final Root Iters :       4919
Solution           :
  Solutions found  :          2 (1 improvements)
  First Solution   : +1.03971965299997e-03   (in run 1, after 1 nodes, 0.29 seconds, depth 0, found by <SpaGreedy>)
  Gap First Sol.   :   infinite
  Gap Last Sol.    :   infinite
  Primal Bound     : +1.03971965299997e-03   (in run 1, after 1 nodes, 0.29 seconds, depth 0, found by <SpaGreedy>)
  Dual Bound       : +2.85192014639320e-03
  Gap              :     174.30 %
  Avg. Gap         :      90.51 % (325837.02 primal-dual integral)

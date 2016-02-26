objective value:                  0.00063135529199998
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
x_26_1                                              1 	(obj:0)
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
epsI                              0.00063135529199998 	(obj:1)
SCIP Status        : solving was interrupted [time limit reached]
Total Time         :    3600.00
  solving          :    3600.00
  presolving       :       0.09 (included in solving)
  reading          :       0.00
  copying          :       0.00 (0 times copied the problem)
Original Problem   :
  Problem name     : ../instances/Pot4/Pot4_T_30B_1t.txt.spa
  Variables        : 141 (140 binary, 0 integer, 0 implicit integer, 1 continuous)
  Constraints      : 263 initial, 263 maximal
  Objective sense  : maximize
Presolved Problem  :
  Problem name     : t_../instances/Pot4/Pot4_T_30B_1t.txt.spa
  Variables        : 910 (132 binary, 0 integer, 16 implicit integer, 762 continuous)
  Constraints      : 3165 initial, 4724 maximal
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
  stuffing         :       0.00       0.00      1          0          0          0          0          0          0          0          0          0
  trivial          :       0.00       0.00      7          0          0          0          0          0          0          0          0          0
  tworowbnd        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  dualfix          :       0.00       0.00      7          4          0          0          0          0          0          0          0          0
  genvbounds       :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  probing          :       0.05       0.00      1          0          0          0          4          0          0          0          0          0
  pseudoobj        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  varbound         :       0.00       0.00      9          0          0          0          0          0         28         28          0          0
  setppc           :       0.00       0.00     12          0          0          0          0          0          0          0          0          0
  and              :       0.00       0.00      4          0         12          0          0          0         12          0          0          0
  linear           :       0.02       0.00      9          3          1          0          4          0          7          0         12         12
  logicor          :       0.00       0.00      7          0          0          0          0          0          0          0          0          0
  bounddisjunction :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  quadratic        :       0.01       0.00      2          0          0          0          0          0          0       2995          0          0
  root node        :          -          -      -          1          -          -         72          -          -          -          -          -
Constraints        :     Number  MaxNumber  #Separate #Propagate    #EnfoLP    #EnfoPS     #Check   #ResProp    Cutoffs    DomReds       Cuts    Applied      Conss   Children
  integral         :          0          0          0          0    1691349          0     194196          0          0      59335          0          0         26    2333480
  varbound         :       1389       1389         82    7489694     521091          0      84819    1569072       7283   35155214          0          0          0          0
  setppc           :        160        160         82    7481684     521091          0      15782    3262284       8543    2628688         12          0          0          0
  and              :         22         22       2378      25703        113          0      13691     118726          2       1903         25          9          0          0
  linear           :       1539       1539         82    7473866     372020          0      13645     686896      19032   97636978   32374133   32373723          0          0
  logicor          :         55+       749         74    2101036          0          0          0     363822       2285     270403          3          3          0          0
  bounddisjunction :          0+       883          0    2163994          0          0          0      90519        528     628221          0          0          0          0
  countsols        :          0          0          0          0          0          0     109366          0          0          0          0          0          0          0
Constraint Timings :  TotalTime  SetupTime   Separate  Propagate     EnfoLP     EnfoPS      Check    ResProp    SB-Prop
  integral         :      25.31       0.00       0.00       0.00      25.27       0.00       0.04       0.00       0.00
  varbound         :      36.87       0.00       0.01      29.30       6.06       0.00       1.06       0.27       0.17
  setppc           :      40.16       0.00       0.00      38.75       0.51       0.00       0.02       0.67       0.21
  and              :       0.06       0.00       0.00       0.00       0.00       0.00       0.03       0.03       0.00
  linear           :     201.16       0.00       0.01     168.11      31.81       0.00       0.20       0.43       0.60
  logicor          :      10.45       0.00       0.00      10.30       0.00       0.00       0.00       0.14       0.01
  bounddisjunction :      16.14       0.00       0.00      16.11       0.00       0.00       0.00       0.02       0.01
  countsols        :       0.01       0.00       0.00       0.00       0.00       0.00       0.01       0.00       0.00
Propagators        : #Propagate   #ResProp    Cutoffs    DomReds
  dualfix          :          2          0          0          0
  genvbounds       :          0          0          0          0
  obbt             :          0          0          0          0
  probing          :          0          0          0          0
  pseudoobj        :    7995048          0          1         71
  redcost          :    3535318          0         38     296595
  rootredcost      :          1          0          0          0
  vbounds          :    9148354         67          0        415
Propagator Timings :  TotalTime  SetupTime   Presolve  Propagate    ResProp    SB-Prop
  dualfix          :       0.00       0.00       0.00       0.00       0.00       0.00
  genvbounds       :       1.31       0.00       0.00       1.31       0.00       0.00
  obbt             :       0.00       0.00       0.00       0.00       0.00       0.00
  probing          :       0.05       0.00       0.05       0.00       0.00       0.00
  pseudoobj        :       2.61       0.00       0.00       2.61       0.00       0.00
  redcost          :      16.20       0.00       0.00      16.20       0.00       0.00
  rootredcost      :       1.46       0.00       0.00       1.46       0.00       0.00
  vbounds          :       3.15       0.00       0.00       3.15       0.00       0.00
Conflict Analysis  :       Time      Calls    Success    DomReds  Conflicts   Literals    Reconvs ReconvLits   LP Iters
  propagation      :       2.30      18743      18692          -     406097       59.3        833       71.2          -
  infeasible LP    :      20.49     118761      70536          -    1166807       45.5         12       10.5          0
  bound exceed. LP :       0.00          0          0          -          0        0.0          0        0.0          0
  strong branching :       0.00          0          0          -          0        0.0          0        0.0          0
  pseudo solution  :       0.00          0          0          -          0        0.0          0        0.0          -
  applied globally :       1.63          -          -          0      90073       22.6          -          -          -
  applied locally  :          -          -          -          0       8207       38.4          -          -          -
Separators         :   ExecTime  SetupTime      Calls    Cutoffs    DomReds       Cuts    Applied      Conss
  cut pool         :       0.00                    80          -          -       1007          -          -    (maximal pool size: 684)
  cgmip            :       0.00       0.00          0          0          0          0          0          0
  clique           :       0.00       0.00         81          0          0         45         38          0
  closecuts        :       0.00       0.00          0          0          0          0          0          0
  cmir             :       0.01       0.00         10          0          0         60         10          0
  disjunctive      :       0.00       0.00          0          0          0          0          0          0
  eccuts           :       0.00       0.00          0          0          0          0          0          0
  flowcover        :       0.44       0.00         15          0          0        259         60          0
  gomory           :       0.06       0.00         10          0          0         15         14          0
  impliedbounds    :       0.03       0.00         81          0          0        679        154          0
  intobj           :       0.00       0.00          0          0          0          0          0          0
  mcf              :       0.00       0.00          1          0          0          0          0          0
  oddcycle         :       0.00       0.00          0          0          0          0          0          0
  rapidlearning    :       0.00       0.00          0          0          0          0          0          0
  strongcg         :       0.14       0.00         20          0          0          8          4          0
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
  relpscost        :      24.81       0.00    1170258          0          0          0      59335          0       -395    2333480
Primal Heuristics  :   ExecTime  SetupTime      Calls      Found       Best
  LP solutions     :       0.00          -          -          0          -
  pseudo solutions :       0.00          -          -          0          -
  strong branching :       0.01          -          -          0          -
  SpaGreedy        :     442.56       0.00     109365          2          1
  actconsdiving    :       0.00       0.00          0          0          0
  bound            :       0.00       0.00          0          0          0
  clique           :       0.00       0.00          0          0          0
  coefdiving       :      40.46       0.00       7968          0          0
  crossover        :       0.32       0.00          0          0          0
  dins             :       0.00       0.00          0          0          0
  distributiondivin:      63.37       0.00       7967          0          0
  dualval          :       0.00       0.00          0          0          0
  feaspump         :      21.23       0.00       3984          0          0
  fixandinfer      :       0.00       0.00          0          0          0
  fracdiving       :      38.82       0.00       7967          0          0
  guideddiving     :       0.04       0.00          0          0          0
  indicator        :       0.19       0.00          0          0          0
  intdiving        :       0.00       0.00          0          0          0
  intshifting      :       0.57       0.00        154          0          0
  linesearchdiving :      39.93       0.00       7967          0          0
  localbranching   :       0.00       0.00          0          0          0
  mutation         :       0.00       0.00          0          0          0
  nlpdiving        :       0.03       0.00          0          0          0
  objpscostdiving  :      29.07       0.00       3821          0          0
  octane           :       0.00       0.00          0          0          0
  ofins            :       0.00       0.00          0          0          0
  oneopt           :       0.50       0.00          0          0          0
  proximity        :       0.00       0.00          0          0          0
  pscostdiving     :      37.88       0.00       7967          0          0
  randrounding     :       1.17       0.00      71069          0          0
  rens             :       0.00       0.00          0          0          0
  reoptsols        :       0.00       0.00          0          0          0
  rins             :       0.44       0.00          0          0          0
  rootsoldiving    :      23.86       0.00       3984          0          0
  rounding         :       0.49       0.00      15180          0          0
  shiftandpropagate:       0.00       0.00          0          0          0
  shifting         :       0.85       0.00       3461          0          0
  simplerounding   :       0.24       0.00          0          0          0
  subnlp           :       0.40       0.00          0          0          0
  trivial          :       0.00       0.00          2          0          0
  trivialnegation  :       0.00       0.00          0          0          0
  trysol           :       0.41       0.00          0          0          0
  twoopt           :       0.00       0.00          0          0          0
  undercover       :       0.00       0.00          1          0          0
  vbounds          :       0.00       0.00          0          0          0
  veclendiving     :      39.05       0.00       7967          0          0
  zeroobj          :       0.00       0.00          0          0          0
  zirounding       :       0.33       0.00       1000          0          0
  other solutions  :          -          -          -          0          -
Diving Statistics  :      Calls      Nodes   LP Iters Backtracks   MinDepth   MaxDepth   AvgDepth  NLeafSols  MinSolDpt  MaxSolDpt  AvgSolDpt
  actconsdiving    :          0          -          -          -          -          -          -          -          -          -          -
  coefdiving       :       7968      73754     891139       9574         19         60       30.4          -          -          -          -
  distributiondivin:       7967      58865     853251       7728         18         50       29.0          -          -          -          -
  fracdiving       :       7967      67258     976296       7717         18         63       29.9          -          -          -          -
  guideddiving     :          0          -          -          -          -          -          -          -          -          -          -
  linesearchdiving :       7967      65580     975746       7572         19         54       29.6          -          -          -          -
  pscostdiving     :       7967      67295     854930       7350         19         65       29.9          -          -          -          -
  veclendiving     :       7967      68051     889105       8171         19         47       29.9          -          -          -          -
LP                 :       Time      Calls Iterations  Iter/call   Iter/sec  Time-0-It Calls-0-It
  primal LP        :       0.00         12          0       0.00          -       0.00         12
  dual LP          :    1606.00    2865677   64309378      27.53   40043.20      77.22     529380
  lex dual LP      :       0.00          0          0       0.00          -
  barrier LP       :       0.00          0          0       0.00          -       0.00          0
  diving/probing LP:     202.79     251179    7183119      28.60   35421.47
  strong branching :      15.43       5023     364634      72.59   23631.50
    (at root node) :          -         19       2090     110.00          -
  conflict analysis:       0.00          0          0       0.00          -
B&B Tree           :
  number of runs   :          1
  nodes            :    2115151 (1166740 internal, 948411 leaves)
  nodes (total)    :    2115151 (1166740 internal, 948411 leaves)
  nodes left       :     216133
  max depth        :         35
  max depth (total):         35
  backtracks       :     550906 (26.0%)
  delayed cutoffs  :       2198
  repropagations   :      66726 (1966249 domain reductions, 1400 cutoffs)
  avg switch length:       5.37
  switching time   :     425.31
Root Node          :
  First LP value   : +3.00000000000000e+00
  First LP Iters   :        622 (62200.00 Iter/sec)
  First LP Time    :       0.01
  Final Dual Bound : +6.06855291545550e-03
  Final Root Iters :       3913
Solution           :
  Solutions found  :          2 (1 improvements)
  First Solution   : +6.31355291999980e-04   (in run 1, after 1 nodes, 0.11 seconds, depth 0, found by <SpaGreedy>)
  Gap First Sol.   :   infinite
  Gap Last Sol.    :   infinite
  Primal Bound     : +6.31355291999980e-04   (in run 1, after 1 nodes, 0.11 seconds, depth 0, found by <SpaGreedy>)
  Dual Bound       : +1.16435967982734e-03
  Gap              :      84.42 %
  Avg. Gap         :      89.60 % (322560.73 primal-dual integral)

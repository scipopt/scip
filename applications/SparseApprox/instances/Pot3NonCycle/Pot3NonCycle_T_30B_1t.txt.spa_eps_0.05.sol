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
Total Time         :    3600.01
  solving          :    3600.00
  presolving       :       0.29 (included in solving)
  reading          :       0.01
  copying          :       0.02 (1 #copies) (minimal 0.02, maximal 0.02, average 0.02)
Original Problem   :
  Problem name     : ../instances/Pot3NonCycle/Pot3NonCycle_T_30B_1t.txt.spa
  Variables        : 141 (140 binary, 0 integer, 0 implicit integer, 1 continuous)
  Constraints      : 263 initial, 263 maximal
  Objective sense  : maximize
Presolved Problem  :
  Problem name     : t_../instances/Pot3NonCycle/Pot3NonCycle_T_30B_1t.txt.spa
  Variables        : 910 (132 binary, 0 integer, 16 implicit integer, 762 continuous)
  Constraints      : 3141 initial, 4454 maximal
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
  probing          :       0.17       0.00      1          0          0          0          4          0          0          0          0          0
  pseudoobj        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  varbound         :       0.01       0.00      9          0          0          0          0          0         28         28          0          0
  setppc           :       0.00       0.00     12          0          0          0          0          0          0          0          0          0
  and              :       0.00       0.00      4          0         12          0          0          0         12          0          0          0
  linear           :       0.05       0.00      9          3          1          0          4          0          7          0         12         12
  logicor          :       0.00       0.00      7          0          0          0          0          0          0          0          0          0
  bounddisjunction :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  quadratic        :       0.04       0.00      2          0          0          0          0          0          0       2971          0          0
  root node        :          -          -      -          1          -          -         65          -          -          -          -          -
Constraints        :     Number  MaxNumber  #Separate #Propagate    #EnfoLP    #EnfoPS     #Check   #ResProp    Cutoffs    DomReds       Cuts    Applied      Conss   Children
  integral         :          0          0          0          0    1591157          0     227062          0          0      78153          0          0         61    3152558
  varbound         :       1365       1365         83    8887476      10111          0     109468    1687998       2804   39850773          0          0          0          0
  setppc           :        160        160         83    8884672      10111          0        402    2864492       5916    2979675          6          0          0          0
  and              :         22         22       1062      15815          0          0        378      94583          9        973         14          9          0          0
  linear           :       1539       1539         83    8878747       3593          0        375     689516       5675   86874195     308259     307823          0          0
  logicor          :         55+       737         82    2108077          0          0          0    1011620        587      97220          2          2          0          0
  bounddisjunction :          0+       745          0    1958828          0          0          0      64390        161     390542          0          0          0          0
  countsols        :          0          0          0          0          0          0     117593          0          0          0          0          0          0          0
Constraint Timings :  TotalTime  SetupTime   Separate  Propagate     EnfoLP     EnfoPS      Check    ResProp    SB-Prop
  integral         :      35.74       0.00       0.00       0.00      35.72       0.00       0.02       0.00       0.00
  varbound         :      40.30       0.00       0.01      39.20       0.18       0.00       0.57       0.18       0.16
  setppc           :      39.49       0.00       0.00      38.65       0.02       0.00       0.00       0.65       0.17
  and              :       0.02       0.00       0.00       0.01       0.00       0.00       0.00       0.01       0.00
  linear           :     186.74       0.00       0.00     184.66       0.55       0.00       0.02       0.60       0.91
  logicor          :       9.19       0.00       0.00       8.55       0.00       0.00       0.00       0.63       0.01
  bounddisjunction :      10.81       0.00       0.00      10.80       0.00       0.00       0.00       0.00       0.01
  countsols        :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
Propagators        : #Propagate   #ResProp    Cutoffs    DomReds
  dualfix          :          2          0          0          0
  genvbounds       :          0          0          0          0
  obbt             :          0          0          0          0
  probing          :          0          0          0          0
  pseudoobj        :    9575250          0          0         64
  redcost          :    3230014          0          0      67613
  rootredcost      :          1          0          0          0
  vbounds          :   10478090        253          0        319
Propagator Timings :  TotalTime  SetupTime   Presolve  Propagate    ResProp    SB-Prop
  dualfix          :       0.00       0.00       0.00       0.00       0.00       0.00
  genvbounds       :       1.66       0.00       0.00       1.65       0.00       0.01
  obbt             :       0.00       0.00       0.00       0.00       0.00       0.00
  probing          :       0.17       0.00       0.17       0.00       0.00       0.00
  pseudoobj        :       4.62       0.00       0.00       4.62       0.00       0.00
  redcost          :      12.17       0.00       0.00      12.17       0.00       0.00
  rootredcost      :       2.33       0.00       0.00       2.33       0.00       0.00
  vbounds          :       3.26       0.00       0.00       3.24       0.00       0.02
Conflict Analysis  :       Time      Calls    Success    DomReds  Conflicts   Literals    Reconvs ReconvLits   LP Iters
  propagation      :       1.69       9480       9440          -     191663       63.2        562       75.4          -
  infeasible LP    :      15.41      89882      53204          -     714277       52.4         18        8.9          0
  bound exceed. LP :       0.00          0          0          -          0        0.0          0        0.0          0
  strong branching :       0.00          0          0          -          0        0.0          0        0.0          0
  pseudo solution  :       0.00          0          0          -          0        0.0          0        0.0          -
  applied globally :       0.58          -          -          0      36152       22.9          -          -          -
  applied locally  :          -          -          -          0      20960       51.4          -          -          -
Separators         :   ExecTime  SetupTime      Calls    Cutoffs    DomReds       Cuts    Applied      Conss
  cut pool         :       0.00                    81          -          -        925          -          -    (maximal pool size: 695)
  cgmip            :       0.00       0.00          0          0          0          0          0          0
  clique           :       0.02       0.00         82          0          0         51         36          0
  closecuts        :       0.00       0.00          0          0          0          0          0          0
  cmir             :       0.03       0.00         10          0          0         84         40          0
  disjunctive      :       0.00       0.00          0          0          0          0          0          0
  eccuts           :       0.00       0.00          0          0          0          0          0          0
  flowcover        :       0.63       0.00         15          0          0        263         37          0
  gomory           :       0.15       0.00         10          0          0         19         16          0
  impliedbounds    :       0.00       0.00         82          0          0        694        131          0
  intobj           :       0.00       0.00          0          0          0          0          0          0
  mcf              :       0.00       0.00          1          0          0          0          0          0
  oddcycle         :       0.00       0.00          0          0          0          0          0          0
  rapidlearning    :       0.00       0.00          0          0          0          0          0          0
  strongcg         :       0.22       0.00         20          0          0          5          3          0
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
  relpscost        :      34.95       0.00    1581046          0          0          0      78153          0       -151    3152558
Primal Heuristics  :   ExecTime  SetupTime      Calls      Found       Best
  LP solutions     :       0.00          -          -          0          -
  pseudo solutions :       0.00          -          -          0          -
  strong branching :       0.02          -          -          0          -
  SpaGreedy        :     156.51       0.00     117592          2          1
  actconsdiving    :       0.00       0.00          0          0          0
  bound            :       0.00       0.00          0          0          0
  clique           :       0.00       0.00          0          0          0
  coefdiving       :      30.64       0.00       8635          0          0
  crossover        :       0.48       0.00          0          0          0
  dins             :       0.00       0.00          0          0          0
  distributiondivin:      71.69       0.00       8636          0          0
  dualval          :       0.00       0.00          0          0          0
  feaspump         :      23.24       0.00       1460          0          0
  fixandinfer      :       0.00       0.00          0          0          0
  fracdiving       :      31.93       0.00       8636          0          0
  guideddiving     :       0.02       0.00          0          0          0
  indicator        :       0.29       0.00          0          0          0
  intdiving        :       0.00       0.00          0          0          0
  intshifting      :       0.22       0.00         61          0          0
  linesearchdiving :      33.18       0.00       8636          0          0
  localbranching   :       0.00       0.00          0          0          0
  mutation         :       0.00       0.00          0          0          0
  nlpdiving        :       0.03       0.00          0          0          0
  objpscostdiving  :      24.38       0.00       1841          0          0
  octane           :       0.00       0.00          0          0          0
  ofins            :       0.00       0.00          0          0          0
  oneopt           :       0.71       0.00          0          0          0
  proximity        :       0.00       0.00          0          0          0
  pscostdiving     :      30.65       0.00       8636          0          0
  randrounding     :       1.56       0.00     109042          0          0
  rens             :       1.00       0.00          1          0          0
  reoptsols        :       0.00       0.00          0          0          0
  rins             :       0.48       0.00          0          0          0
  rootsoldiving    :      29.06       0.00       2071          0          0
  rounding         :       0.87       0.00      17636          0          0
  shiftandpropagate:       0.00       0.00          0          0          0
  shifting         :       2.08       0.00       4388          0          0
  simplerounding   :       0.18       0.00          0          0          0
  subnlp           :       0.45       0.00          0          0          0
  trivial          :       0.00       0.00          2          0          0
  trivialnegation  :       0.00       0.00          0          0          0
  trysol           :       0.61       0.00          0          0          0
  twoopt           :       0.00       0.00          0          0          0
  undercover       :       0.00       0.00          1          0          0
  vbounds          :       0.00       0.00          0          0          0
  veclendiving     :      30.63       0.00       8636          0          0
  zeroobj          :       0.00       0.00          0          0          0
  zirounding       :       0.46       0.00       1000          0          0
  other solutions  :          -          -          -          0          -
Diving Statistics  :      Calls      Nodes   LP Iters Backtracks   MinDepth   MaxDepth   AvgDepth  NLeafSols  MinSolDpt  MaxSolDpt  AvgSolDpt
  actconsdiving    :          0          -          -          -          -          -          -          -          -          -          -
  coefdiving       :       8635      57266     594141       9151         19         57       27.6          -          -          -          -
  distributiondivin:       8636      44848     579162       9343         15         40       26.2          -          -          -          -
  fracdiving       :       8636      54783     648408       9385         14         66       27.3          -          -          -          -
  guideddiving     :          0          -          -          -          -          -          -          -          -          -          -
  linesearchdiving :       8636      55894     631306       9336         17         38       27.4          -          -          -          -
  pscostdiving     :       8636      48035     587744       8813         17         49       26.6          -          -          -          -
  veclendiving     :       8636      55162     597554       9234         16         41       27.3          -          -          -          -
LP                 :       Time      Calls Iterations  Iter/call   Iter/sec  Time-0-It Calls-0-It
  primal LP        :       0.01         12          0       0.00       0.00       0.01         12
  dual LP          :    1712.13    3224716   46988517      15.74   27444.48      55.13     240187
  lex dual LP      :       0.00          0          0       0.00          -
  barrier LP       :       0.00          0          0       0.00          -       0.00          0
  diving/probing LP:     190.57     167434    5129578      30.64   26917.03
  strong branching :      18.81       7736     324295      41.92   17240.56
    (at root node) :          -         20       2036     101.80          -
  conflict analysis:       0.00          0          0       0.00          -
B&B Tree           :
  number of runs   :          1
  nodes            :    2965385 (1576279 internal, 1389106 leaves)
  nodes (total)    :    2965385 (1576279 internal, 1389106 leaves)
  nodes left       :     187151
  max depth        :         44
  max depth (total):         44
  backtracks       :     775907 (26.2%)
  delayed cutoffs  :         23
  repropagations   :      18286 (316020 domain reductions, 17 cutoffs)
  avg switch length:       5.60
  switching time   :     584.81
Root Node          :
  First LP value   : +3.00000000000000e+00
  First LP Iters   :        581 (29050.00 Iter/sec)
  First LP Time    :       0.02
  Final Dual Bound : +3.59814281799942e-03
  Final Root Iters :       3737
Solution           :
  Solutions found  :          2 (1 improvements)
  First Solution   : +2.75616363000039e-04   (in run 1, after 1 nodes, 0.29 seconds, depth 0, found by <SpaGreedy>)
  Gap First Sol.   :   infinite
  Gap Last Sol.    :   infinite
  Primal Bound     : +2.75616363000039e-04   (in run 1, after 1 nodes, 0.29 seconds, depth 0, found by <SpaGreedy>)
  Dual Bound       : +4.05315119886364e-04
  Gap              :      47.06 %
  Avg. Gap         :      92.34 % (332441.61 primal-dual integral)

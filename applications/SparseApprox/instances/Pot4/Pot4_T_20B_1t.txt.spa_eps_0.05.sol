objective value:                 0.000608190768000005
x_1_1                                               1 	(obj:0)
x_2_2                                               1 	(obj:0)
x_3_2                                               1 	(obj:0)
x_4_4                                               1 	(obj:0)
x_5_1                                               1 	(obj:0)
x_6_2                                               1 	(obj:0)
x_7_4                                               1 	(obj:0)
x_8_1                                               1 	(obj:0)
x_9_1                                               1 	(obj:0)
x_10_2                                              1 	(obj:0)
x_11_4                                              1 	(obj:0)
x_12_2                                              1 	(obj:0)
x_13_4                                              1 	(obj:0)
x_14_2                                              1 	(obj:0)
x_15_1                                              1 	(obj:0)
x_16_2                                              1 	(obj:0)
x_17_4                                              1 	(obj:0)
x_18_4                                              1 	(obj:0)
x_19_2                                              1 	(obj:0)
x_20_1                                              1 	(obj:0)
ind_1                                               1 	(obj:0)
abs_1_2                                             1 	(obj:0)
ind_2                                               1 	(obj:0)
abs_2_4                                             1 	(obj:0)
ind_4                                               1 	(obj:0)
abs_4_1                                             1 	(obj:0)
epsI                             0.000608190768000005 	(obj:1)
SCIP Status        : solving was interrupted [user interrupt]
Total Time         :      14.01
  solving          :      14.00
  presolving       :       0.09 (included in solving)
  reading          :       0.01
  copying          :       0.05 (7 #copies) (minimal 0.00, maximal 0.01, average 0.01)
Original Problem   :
  Problem name     : instances/Pot4/Pot4_T_20B_1t.txt.spa
  Variables        : 101 (100 binary, 0 integer, 0 implicit integer, 1 continuous)
  Constraints      : 193 initial, 193 maximal
  Objective sense  : maximize
Presolved Problem  :
  Problem name     : t_instances/Pot4/Pot4_T_20B_1t.txt.spa
  Variables        : 590 (92 binary, 0 integer, 16 implicit integer, 482 continuous)
  Constraints      : 2039 initial, 2333 maximal
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
  pseudoobj        :       0.00       0.00      2          0          0          0          0          0          0          0          0          0
  varbound         :       0.00       0.00      9          0          0          0          0          0         28         28          0          0
  setppc           :       0.00       0.00     12          0          0          0          0          0          0          0          0          0
  and              :       0.00       0.00      4          0         12          0          0          0         12          0          0          0
  linear           :       0.01       0.00      9          3          1          0          4          0          7          0         12         12
  logicor          :       0.00       0.00      7          0          0          0          0          0          0          0          0          0
  bounddisjunction :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  quadratic        :       0.00       0.00      2          0          0          0          0          0          0       1939          0          0
  root node        :          -          -      -          1          -          -         59          -          -          -          -          -
Constraints        :     Number  MaxNumber  #Separate #Propagate    #EnfoLP    #EnfoPS     #Check   #ResProp    Cutoffs    DomReds       Cuts    Applied      Conss   Children
  integral         :          0          0          0          0       2938          0       2226          0          0       1988          0          0          0       5460
  varbound         :        893        893       1352      29389          0          0       1423       3015         32      52527          0          0          0          0
  setppc           :        110        110       1352      29353          0          0         16       3530         50       4709          5          0          0          0
  and              :         22         22       1040       9434          0          0         14        400          5        707         29         17          0          0
  linear           :        979        979       1352      29301          0          0          7       1534         81     108363      20261       9672          0          0
  logicor          :         35+       156       1352       7676          0          0          0        343         10         65         48         17          0          0
  bounddisjunction :          0+       195          0       6841          0          0          0          8          0         30          0          0          0          0
  countsols        :          0          0          0          0          0          0        609          0          0          0          0          0          0          0
Constraint Timings :  TotalTime  SetupTime   Separate  Propagate     EnfoLP     EnfoPS      Check    ResProp    SB-Prop
  integral         :       1.50       0.00       0.00       0.00       1.50       0.00       0.00       0.00       0.00
  varbound         :       0.18       0.00       0.04       0.13       0.00       0.00       0.00       0.00       0.01
  setppc           :       0.20       0.00       0.00       0.18       0.00       0.00       0.00       0.00       0.02
  and              :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
  linear           :       0.64       0.00       0.19       0.42       0.00       0.00       0.00       0.00       0.03
  logicor          :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
  bounddisjunction :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
  countsols        :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
Propagators        : #Propagate   #ResProp    Cutoffs    DomReds
  dualfix          :          2          0          0          0
  genvbounds       :          0          0          0          0
  obbt             :          0          0          0          0
  probing          :          0          0          0          0
  pseudoobj        :      15614          0          0         58
  redcost          :       7258          0          0         37
  rootredcost      :          1          0          0          0
  vbounds          :      31934          2          0        246
Propagator Timings :  TotalTime  SetupTime   Presolve  Propagate    ResProp    SB-Prop
  dualfix          :       0.00       0.00       0.00       0.00       0.00       0.00
  genvbounds       :       0.01       0.00       0.00       0.01       0.00       0.00
  obbt             :       0.00       0.00       0.00       0.00       0.00       0.00
  probing          :       0.03       0.00       0.03       0.00       0.00       0.00
  pseudoobj        :       0.01       0.00       0.00       0.00       0.00       0.01
  redcost          :       0.01       0.00       0.00       0.01       0.00       0.00
  rootredcost      :       0.01       0.00       0.00       0.01       0.00       0.00
  vbounds          :       0.03       0.00       0.00       0.03       0.00       0.00
Conflict Analysis  :       Time      Calls    Success    DomReds  Conflicts   Literals    Reconvs ReconvLits   LP Iters
  propagation      :       0.02        100         96          -       1229       41.2          6       37.2          -
  infeasible LP    :       0.01        116         81          -       1047       35.0          0        0.0          0
  bound exceed. LP :       0.00          0          0          -          0        0.0          0        0.0          0
  strong branching :       0.00          0          0          -          0        0.0          0        0.0          0
  pseudo solution  :       0.00          0          0          -          0        0.0          0        0.0          -
  applied globally :       0.00          -          -          0        322       18.8          -          -          -
  applied locally  :          -          -          -          0          8       35.1          -          -          -
Separators         :   ExecTime  SetupTime      Calls    Cutoffs    DomReds       Cuts    Applied      Conss
  cut pool         :       0.00                  1346          -          -      19063          -          -    (maximal pool size: 771)
  cgmip            :       0.00       0.00          0          0          0          0          0          0
  clique           :       0.03       0.00         82          0          0         34         26          0
  closecuts        :       0.09       0.00          7          0          0        251          0          0
  cmir             :       0.07       0.00         15          0          0        351       9766          0
  disjunctive      :       0.00       0.00          0          0          0          0          0          0
  eccuts           :       0.00       0.00          0          0          0          0          0          0
  flowcover        :       0.18       0.00         22          0          0        275         30          0
  gomory           :       0.07       0.00         15          0          0        143        143          0
  impliedbounds    :       0.00       0.00         82          0          0        539        141          0
  intobj           :       0.00       0.00          0          0          0          0          0          0
  mcf              :       0.00       0.00          1          0          0          0          0          0
  oddcycle         :       0.01       0.00          8          0          0          2          1          0
  rapidlearning    :       0.00       0.00          0          0          0          0          0          0
  strongcg         :       0.09       0.00         30          0          0         31         16          0
  zerohalf         :       0.08       0.00         15          0          0          0          0          0
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
  relpscost        :       1.50       0.00       2938          0          0          0       1988          0          0       5460
Primal Heuristics  :   ExecTime  SetupTime      Calls      Found       Best
  LP solutions     :       0.00          -          -          0          -
  pseudo solutions :       0.00          -          -          0          -
  strong branching :       0.00          -          -          0          -
  SpaGreedy        :       0.73       0.00        607          2          1
  actconsdiving    :       0.18       0.00         15          0          0
  bound            :       0.00       0.00          1          0          0
  clique           :       0.01       0.00          1          0          0
  coefdiving       :       0.36       0.00         60          0          0
  crossover        :       0.04       0.00          1          0          0
  dins             :       0.01       0.00          1          0          0
  distributiondivin:       0.43       0.00         59          0          0
  dualval          :       0.00       0.00          0          0          0
  feaspump         :       0.25       0.00         16          0          0
  fixandinfer      :       0.00       0.00          0          0          0
  fracdiving       :       0.43       0.00         60          0          0
  guideddiving     :       0.00       0.00          0          0          0
  indicator        :       0.00       0.00          0          0          0
  intdiving        :       0.13       0.00         16          0          0
  intshifting      :       0.03       0.00         22          0          0
  linesearchdiving :       0.68       0.00         41          0          0
  localbranching   :       0.00       0.00          0          0          0
  mutation         :       0.00       0.00          0          0          0
  nlpdiving        :       0.00       0.00          0          0          0
  objpscostdiving  :       0.26       0.00          9          0          0
  octane           :       0.00       0.00          0          0          0
  ofins            :       0.00       0.00          0          0          0
  oneopt           :       0.00       0.00          0          0          0
  proximity        :       0.00       0.00          0          0          0
  pscostdiving     :       0.65       0.00         55          0          0
  randrounding     :       0.27       0.00       1607          0          0
  rens             :       0.07       0.00          1          0          0
  reoptsols        :       0.00       0.00          0          0          0
  rins             :       0.00       0.00          0          0          0
  rootsoldiving    :       0.22       0.00         15          0          0
  rounding         :       0.04       0.00        892          0          0
  shiftandpropagate:       0.00       0.00          0          0          0
  shifting         :       0.13       0.00        532          0          0
  simplerounding   :       0.00       0.00          0          0          0
  subnlp           :       0.01       0.00          0          0          0
  trivial          :       0.00       0.00          2          0          0
  trivialnegation  :       0.00       0.00          0          0          0
  trysol           :       0.01       0.00          0          0          0
  twoopt           :       0.00       0.00          0          0          0
  undercover       :       0.00       0.00          1          0          0
  vbounds          :       0.01       0.00          1          0          0
  veclendiving     :       0.39       0.00         60          0          0
  zeroobj          :       0.04       0.00          1          1          1
  zirounding       :       0.04       0.00       1000          0          0
  other solutions  :          -          -          -          0          -
Diving Statistics  :      Calls      Nodes   LP Iters Backtracks   MinDepth   MaxDepth   AvgDepth  NLeafSols  MinSolDpt  MaxSolDpt  AvgSolDpt
  actconsdiving    :         15        206       4916         28         19         26       21.3          -          -          -          -
  coefdiving       :         60        923      12802         71         19         48       25.4          -          -          -          -
  distributiondivin:         59        658      11161         77         15         44       21.5          -          -          -          -
  fracdiving       :         60        746      12769         93         18         53       22.7          -          -          -          -
  guideddiving     :          0          -          -          -          -          -          -          -          -          -          -
  linesearchdiving :         41        449      17807         69         15         27       20.4          -          -          -          -
  pscostdiving     :         55        867      17885         55         18         53       26.1          -          -          -          -
  veclendiving     :         60        658      10867         87         18         26       20.9          -          -          -          -
LP                 :       Time      Calls Iterations  Iter/call   Iter/sec  Time-0-It Calls-0-It
  primal LP        :       0.03        248         15       2.50     500.00       0.03        242
  dual LP          :       4.78       6380     172397      28.41   36066.32       0.04        312
  lex dual LP      :       0.00          0          0       0.00          -
  barrier LP       :       0.00          0          0       0.00          -       0.00          0
  diving/probing LP:       3.23       2737     119573      43.69   37019.50
  strong branching :       1.16        428      40527      94.69   34937.07
    (at root node) :          -         19       1659      87.32          -
  conflict analysis:       0.00          0          0       0.00          -
B&B Tree           :
  number of runs   :          1
  nodes            :       4445 (2730 internal, 1715 leaves)
  nodes (total)    :       4445 (2730 internal, 1715 leaves)
  nodes left       :       1016
  max depth        :         19
  max depth (total):         19
  backtracks       :       1272 (28.6%)
  delayed cutoffs  :          0
  repropagations   :         21 (298 domain reductions, 0 cutoffs)
  avg switch length:       4.05
  switching time   :       0.47
Root Node          :
  First LP value   : +3.00000000000000e+00
  First LP Iters   :        402
  First LP Time    :       0.00
  Final Dual Bound : +5.88514652887867e-03
  Final Root Iters :       3521
Solution           :
  Solutions found  :          3 (2 improvements)
  First Solution   : +0.00000000000000e+00   (in run 1, after 0 nodes, 0.04 seconds, depth 0, found by <zeroobj>)
  Gap First Sol.   :   infinite
  Gap Last Sol.    :   infinite
  Primal Bound     : +6.08190768000005e-04   (in run 1, after 1 nodes, 0.09 seconds, depth 0, found by <SpaGreedy>)
  Dual Bound       : +1.46701540314461e-03
  Gap              :     141.21 %
  Avg. Gap         :      90.47 % (1266.65 primal-dual integral)

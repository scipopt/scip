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
Total Time         :      10.68
  solving          :      10.68
  presolving       :       0.08 (included in solving)
  reading          :       0.00
  copying          :       0.42 (59 #copies) (minimal 0.00, maximal 0.01, average 0.01)
Original Problem   :
  Problem name     : instances/Pot4/Pot4_T_20B_1t.txt.spa
  Variables        : 101 (100 binary, 0 integer, 0 implicit integer, 1 continuous)
  Constraints      : 193 initial, 193 maximal
  Objective sense  : maximize
Presolved Problem  :
  Problem name     : t_instances/Pot4/Pot4_T_20B_1t.txt.spa
  Variables        : 590 (92 binary, 0 integer, 16 implicit integer, 482 continuous)
  Constraints      : 2039 initial, 2165 maximal
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
  pseudoobj        :       0.00       0.00      2          0          0          0          0          0          0          0          0          0
  varbound         :       0.00       0.00      9          0          0          0          0          0         28         28          0          0
  setppc           :       0.00       0.00     12          0          0          0          0          0          0          0          0          0
  and              :       0.00       0.00      4          0         12          0          0          0         12          0          0          0
  linear           :       0.01       0.00      9          3          1          0          4          0          7          0         12         12
  logicor          :       0.00       0.00      7          0          0          0          0          0          0          0          0          0
  bounddisjunction :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  quadratic        :       0.01       0.00      2          0          0          0          0          0          0       1939          0          0
  root node        :          -          -      -          1          -          -         23          -          -          -          -          -
Constraints        :     Number  MaxNumber  #Separate #Propagate    #EnfoLP    #EnfoPS     #Check   #ResProp    Cutoffs    DomReds       Cuts    Applied      Conss   Children
  integral         :          0          0          0          0       1484          0       1966          0          0       2152          0          0          0       2550
  varbound         :        893        893       1014      13521          0          0       1679       1335         18      27193          0          0          0          0
  setppc           :        110        110       1014      13503          0          0         83       1366         17       2265          4          0          0          0
  and              :         22         22        508       5286          0          0         16        161          0        395         25         16          0          0
  linear           :        979        979       1014      13486          0          0         10        631         62      53206      11170       5437          0          0
  logicor          :         35+       101        921       2745          0          0          0         30          3        114         17         11          0          0
  bounddisjunction :          0+        81          0       2276          0          0          0          1          0        119          0          0          0          0
  countsols        :          0          0          0          0          0          0        286          0          0          0          0          0          0          0
Constraint Timings :  TotalTime  SetupTime   Separate  Propagate     EnfoLP     EnfoPS      Check    ResProp    SB-Prop
  integral         :       1.60       0.00       0.00       0.00       1.60       0.00       0.00       0.00       0.00
  varbound         :       0.09       0.00       0.00       0.08       0.00       0.00       0.01       0.00       0.00
  setppc           :       0.08       0.00       0.00       0.06       0.00       0.00       0.00       0.00       0.02
  and              :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
  linear           :       0.24       0.00       0.09       0.13       0.00       0.00       0.00       0.00       0.02
  logicor          :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
  bounddisjunction :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
  countsols        :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
Propagators        : #Propagate   #ResProp    Cutoffs    DomReds
  dualfix          :          2          0          0          0
  genvbounds       :          0          0          0          0
  obbt             :          0          0          0          0
  probing          :          0          0          0          0
  pseudoobj        :       7041          0          0         22
  redcost          :       3989          0          0          6
  rootredcost      :          1          0          0          0
  vbounds          :      14614          7          0         44
Propagator Timings :  TotalTime  SetupTime   Presolve  Propagate    ResProp    SB-Prop
  dualfix          :       0.00       0.00       0.00       0.00       0.00       0.00
  genvbounds       :       0.00       0.00       0.00       0.00       0.00       0.00
  obbt             :       0.00       0.00       0.00       0.00       0.00       0.00
  probing          :       0.03       0.00       0.03       0.00       0.00       0.00
  pseudoobj        :       0.00       0.00       0.00       0.00       0.00       0.00
  redcost          :       0.01       0.00       0.00       0.01       0.00       0.00
  rootredcost      :       0.01       0.00       0.00       0.01       0.00       0.00
  vbounds          :       0.00       0.00       0.00       0.00       0.00       0.00
Conflict Analysis  :       Time      Calls    Success    DomReds  Conflicts   Literals    Reconvs ReconvLits   LP Iters
  propagation      :       0.00         38         38          -        289       40.3          3       37.7          -
  infeasible LP    :       0.00         70         43          -        564       34.9          0        0.0          0
  bound exceed. LP :       0.00          0          0          -          0        0.0          0        0.0          0
  strong branching :       0.00          0          0          -          0        0.0          0        0.0          0
  pseudo solution  :       0.00          0          0          -          0        0.0          0        0.0          -
  applied globally :       0.00          -          -          0        149       16.7          -          -          -
  applied locally  :          -          -          -          0          6       33.0          -          -          -
Separators         :   ExecTime  SetupTime      Calls    Cutoffs    DomReds       Cuts    Applied      Conss
  cut pool         :       0.00                  1010          -          -       8575          -          -    (maximal pool size: 576)
  cgmip            :       0.00       0.00          0          0          0          0          0          0
  clique           :       0.01       0.00         27          0          0         27         22          0
  closecuts        :       0.07       0.00          5          0          0        247          0          0
  cmir             :       0.04       0.00         15          0          0        378       5037          0
  disjunctive      :       0.00       0.00          0          0          0          0          0          0
  eccuts           :       0.00       0.00          0          0          0          0          0          0
  flowcover        :       0.24       0.00         22          0          0        327         34          0
  gomory           :       0.04       0.00         15          0          0         98         98          0
  impliedbounds    :       0.00       0.00         27          0          0        237         62          0
  intobj           :       0.00       0.00          0          0          0          0          0          0
  mcf              :       0.00       0.00          1          0          0          0          0          0
  oddcycle         :       0.00       0.00          7          0          0          2          2          0
  rapidlearning    :       0.00       0.00          0          0          0          0          0          0
  strongcg         :       0.13       0.00         22          0          0          5          3          0
  zerohalf         :       0.07       0.00         15          0          0          0          0          0
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
  relpscost        :       1.60       0.00       1484          0          0          0       2152          0          0       2550
Primal Heuristics  :   ExecTime  SetupTime      Calls      Found       Best
  LP solutions     :       0.00          -          -          0          -
  pseudo solutions :       0.00          -          -          0          -
  strong branching :       0.00          -          -          0          -
  SpaGreedy        :       0.33       0.00        284          2          1
  actconsdiving    :       0.05       0.00          6          0          0
  bound            :       0.00       0.00        152          0          0
  clique           :       0.01       0.00          2          0          0
  coefdiving       :       0.09       0.00         23          0          0
  crossover        :       0.03       0.00          1          0          0
  dins             :       0.45       0.00         51          0          0
  distributiondivin:       0.09       0.00         23          0          0
  dualval          :       0.00       0.00          0          0          0
  feaspump         :       0.17       0.00         11          0          0
  fixandinfer      :       0.00       0.00          0          0          0
  fracdiving       :       0.18       0.00         23          0          0
  guideddiving     :       0.00       0.00          0          0          0
  indicator        :       0.00       0.00          0          0          0
  intdiving        :       0.01       0.00          6          0          0
  intshifting      :       0.04       0.00         15          0          0
  linesearchdiving :       0.28       0.00         23          0          0
  localbranching   :       0.01       0.00          0          0          0
  mutation         :       0.00       0.00          0          0          0
  nlpdiving        :       0.00       0.00          0          0          0
  objpscostdiving  :       0.22       0.00          5          0          0
  octane           :       0.00       0.00          0          0          0
  ofins            :       0.00       0.00          0          0          0
  oneopt           :       0.00       0.00          0          0          0
  proximity        :       0.00       0.00          0          0          0
  pscostdiving     :       0.30       0.00         23          0          0
  randrounding     :       0.06       0.00       1218          0          0
  rens             :       1.81       0.00          4          0          0
  reoptsols        :       0.00       0.00          0          0          0
  rins             :       0.00       0.00          0          0          0
  rootsoldiving    :       0.08       0.00         12          0          0
  rounding         :       0.01       0.00        653          0          0
  shiftandpropagate:       0.00       0.00          0          0          0
  shifting         :       0.22       0.00        495          0          0
  simplerounding   :       0.00       0.00          0          0          0
  subnlp           :       0.00       0.00          0          0          0
  trivial          :       0.02       0.00        153          0          0
  trivialnegation  :       0.00       0.00          0          0          0
  trysol           :       0.00       0.00          0          0          0
  twoopt           :       0.00       0.00          0          0          0
  undercover       :       0.00       0.00          2          0          0
  vbounds          :       0.02       0.00        152          0          0
  veclendiving     :       0.18       0.00         23          0          0
  zeroobj          :       0.03       0.00          1          1          1
  zirounding       :       0.04       0.00       1000          0          0
  other solutions  :          -          -          -          0          -
Diving Statistics  :      Calls      Nodes   LP Iters Backtracks   MinDepth   MaxDepth   AvgDepth  NLeafSols  MinSolDpt  MaxSolDpt  AvgSolDpt
  actconsdiving    :          6         72       1238         15         17         28       22.8          -          -          -          -
  coefdiving       :         23        272       2950         37         17         43       24.6          -          -          -          -
  distributiondivin:         23        209       3290         24         16         36       22.4          -          -          -          -
  fracdiving       :         23        258       5210         40         17         35       24.8          -          -          -          -
  guideddiving     :          0          -          -          -          -          -          -          -          -          -          -
  linesearchdiving :         23        235       7862         37         17         27       23.0          -          -          -          -
  pscostdiving     :         23        466       9273         25         16         52       32.3          -          -          -          -
  veclendiving     :         23        266       4756         47         15         26       22.2          -          -          -          -
LP                 :       Time      Calls Iterations  Iter/call   Iter/sec  Time-0-It Calls-0-It
  primal LP        :       0.00        189          2       1.00          -       0.00        187
  dual LP          :       2.72       3296      99461      31.67   36566.54       0.02        155
  lex dual LP      :       0.00          0          0       0.00          -
  barrier LP       :       0.00          0          0       0.00          -       0.00          0
  diving/probing LP:       1.43       1286      51937      40.39   36319.58
  strong branching :       1.39        411      45874     111.62   33002.88
    (at root node) :          -         24       1788      74.50          -
  conflict analysis:       0.00          0          0       0.00          -
B&B Tree           :
  number of runs   :          1
  nodes            :       1931 (1275 internal, 656 leaves)
  nodes (total)    :       1931 (1275 internal, 656 leaves)
  nodes left       :        620
  max depth        :         24
  max depth (total):         24
  backtracks       :        571 (29.6%)
  delayed cutoffs  :          0
  repropagations   :          8 (40 domain reductions, 0 cutoffs)
  avg switch length:       5.55
  switching time   :       0.27
Root Node          :
  First LP value   : +3.00000000000000e+00
  First LP Iters   :        402
  First LP Time    :       0.00
  Final Dual Bound : +2.18820857024034e-01
  Final Root Iters :       2261
Solution           :
  Solutions found  :          3 (2 improvements)
  First Solution   : +0.00000000000000e+00   (in run 1, after 0 nodes, 0.03 seconds, depth 0, found by <zeroobj>)
  Gap First Sol.   :   infinite
  Gap Last Sol.    :   infinite
  Primal Bound     : +6.08190768000005e-04   (in run 1, after 1 nodes, 0.08 seconds, depth 0, found by <SpaGreedy>)
  Dual Bound       : +1.63721326371650e-03
  Gap              :     169.19 %
  Avg. Gap         :      99.73 % (1065.17 primal-dual integral)

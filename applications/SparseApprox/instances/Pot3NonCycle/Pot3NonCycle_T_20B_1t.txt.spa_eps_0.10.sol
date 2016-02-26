objective value:                 0.000288402607999991
x_1_1                                               1 	(obj:0)
x_2_3                                               1 	(obj:0)
x_3_2                                               1 	(obj:0)
x_4_3                                               1 	(obj:0)
x_5_2                                               1 	(obj:0)
x_6_1                                               1 	(obj:0)
x_7_3                                               1 	(obj:0)
x_8_2                                               1 	(obj:0)
x_9_1                                               1 	(obj:0)
x_10_3                                              1 	(obj:0)
x_11_1                                              1 	(obj:0)
x_12_2                                              1 	(obj:0)
x_13_2                                              1 	(obj:0)
x_14_2                                              1 	(obj:0)
x_15_3                                              1 	(obj:0)
x_16_1                                              1 	(obj:0)
x_17_3                                              1 	(obj:0)
x_18_2                                              1 	(obj:0)
x_19_1                                              1 	(obj:0)
x_20_3                                              1 	(obj:0)
ind_1                                               1 	(obj:0)
abs_1_3                                             1 	(obj:0)
ind_2                                               1 	(obj:0)
abs_2_1                                             1 	(obj:0)
ind_3                                               1 	(obj:0)
abs_3_2                                             1 	(obj:0)
epsI                             0.000288402607999991 	(obj:1)
SCIP Status        : problem is solved [optimal solution found]
Total Time         :      40.50
  solving          :      40.50
  presolving       :       0.04 (included in solving)
  reading          :       0.00
  copying          :       0.00 (0 times copied the problem)
Original Problem   :
  Problem name     : instances/Pot3NonCycle/Pot3NonCycle_T_20B_1t.txt.spa
  Variables        : 73 (72 binary, 0 integer, 0 implicit integer, 1 continuous)
  Constraints      : 128 initial, 128 maximal
  Objective sense  : maximize
Presolved Problem  :
  Problem name     : t_instances/Pot3NonCycle/Pot3NonCycle_T_20B_1t.txt.spa
  Variables        : 340 (67 binary, 0 integer, 9 implicit integer, 264 continuous)
  Constraints      : 1124 initial, 2910 maximal
Presolvers         :   ExecTime  SetupTime  Calls  FixedVars   AggrVars   ChgTypes  ChgBounds   AddHoles    DelCons    AddCons   ChgSides   ChgCoefs
  boundshift       :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  components       :       0.00       0.00      1          0          0          0          0          0          0          0          0          0
  convertinttobin  :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  domcol           :       0.01       0.00      1          0          0          0          0          0          0          0          0          0
  dualagg          :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  dualinfer        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  gateextraction   :       0.00       0.00      3          0          0          0          0          0         57         18          0          0
  implfree         :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  implics          :       0.00       0.00      4          0          0          0          0          0          0          0          0          0
  inttobinary      :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  redvub           :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  stuffing         :       0.00       0.00      1          0          0          0          0          0          0          0          0          0
  trivial          :       0.00       0.00      7          0          0          0          0          0          0          0          0          0
  tworowbnd        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  dualfix          :       0.00       0.00      7          3          0          0          0          0          0          0          0          0
  genvbounds       :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  probing          :       0.01       0.00      1          0          0          0          0          0          0          0          0          0
  pseudoobj        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  varbound         :       0.00       0.00      9          0          0          0          0          0         15         15          0          0
  setppc           :       0.00       0.00     12          0          0          0          0          0          0          0          0          0
  and              :       0.00       0.00      4          0          6          0          0          0          6          0          0          0
  linear           :       0.01       0.00      9          1          1          0          2          0          3          0          6          6
  logicor          :       0.00       0.00      7          0          0          0          0          0          0          0          0          0
  bounddisjunction :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  quadratic        :       0.01       0.00      2          0          0          0          0          0          0       1044          0          0
  root node        :          -          -      -          5          -          -        144          -          -          -          -          -
Constraints        :     Number  MaxNumber  #Separate #Propagate    #EnfoLP    #EnfoPS     #Check   #ResProp    Cutoffs    DomReds       Cuts    Applied      Conss   Children
  integral         :          0          0          0          0      28099          0       2226          0          0       3827          0          0         13      55272
  varbound         :        473        473         69     184791          3          0        242      58698        132     698223          0          0          0          0
  setppc           :         85         85         69     184657          3          0         11      43962        232      58263          0          0          0          0
  and              :         12         12        232       3875          0          0          8      13170          3        241          9          5          0          0
  linear           :        536        536         69     184426          3          0          6      15946        367    1290602        356        242          0          0
  logicor          :         18+       591         56      52007          0          0          0       1734        106       5289          0          0          0          0
  bounddisjunction :          0+      1490          0      74228          0          0          0       1501          4       5154          0          0          0          0
  countsols        :          0          0          0          0          0          0       1983          0          0          0          0          0          0          0
Constraint Timings :  TotalTime  SetupTime   Separate  Propagate     EnfoLP     EnfoPS      Check    ResProp    SB-Prop
  integral         :       0.95       0.00       0.00       0.00       0.95       0.00       0.00       0.00       0.00
  varbound         :       0.70       0.00       0.00       0.69       0.00       0.00       0.00       0.01       0.00
  setppc           :       0.44       0.00       0.00       0.43       0.00       0.00       0.00       0.01       0.00
  and              :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
  linear           :       1.99       0.00       0.00       1.95       0.00       0.00       0.00       0.02       0.02
  logicor          :       0.14       0.00       0.00       0.14       0.00       0.00       0.00       0.00       0.00
  bounddisjunction :       0.49       0.00       0.00       0.48       0.00       0.00       0.00       0.00       0.01
  countsols        :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
Propagators        : #Propagate   #ResProp    Cutoffs    DomReds
  dualfix          :          2          0          0          0
  genvbounds       :          0          0          0          0
  obbt             :          0          0          0          0
  probing          :          0          0          0          0
  pseudoobj        :     182060          0          0         63
  redcost          :      60935          0          0       8655
  rootredcost      :          1          0          0          0
  vbounds          :     212619         91          0         42
Propagator Timings :  TotalTime  SetupTime   Presolve  Propagate    ResProp    SB-Prop
  dualfix          :       0.00       0.00       0.00       0.00       0.00       0.00
  genvbounds       :       0.07       0.00       0.00       0.07       0.00       0.00
  obbt             :       0.00       0.00       0.00       0.00       0.00       0.00
  probing          :       0.01       0.00       0.01       0.00       0.00       0.00
  pseudoobj        :       0.02       0.00       0.00       0.02       0.00       0.00
  redcost          :       0.17       0.00       0.00       0.17       0.00       0.00
  rootredcost      :       0.03       0.00       0.00       0.03       0.00       0.00
  vbounds          :       0.06       0.00       0.00       0.06       0.00       0.00
Conflict Analysis  :       Time      Calls    Success    DomReds  Conflicts   Literals    Reconvs ReconvLits   LP Iters
  propagation      :       0.04        479        477          -       5772       34.4         15       41.3          -
  infeasible LP    :       0.28       4527       3961          -      33049       22.7          1        7.0          0
  bound exceed. LP :       0.00          0          0          -          0        0.0          0        0.0          0
  strong branching :       0.00          0          0          -          0        0.0          0        0.0          0
  pseudo solution  :       0.00          0          0          -          0        0.0          0        0.0          -
  applied globally :       0.06          -          -          0      10955       15.4          -          -          -
  applied locally  :          -          -          -          0         22       34.1          -          -          -
Separators         :   ExecTime  SetupTime      Calls    Cutoffs    DomReds       Cuts    Applied      Conss
  cut pool         :       0.00                    67          -          -        805          -          -    (maximal pool size: 681)
  cgmip            :       0.00       0.00          0          0          0          0          0          0
  clique           :       0.01       0.00         68          0          0          3          3          0
  closecuts        :       0.00       0.00          0          0          0          0          0          0
  cmir             :       0.04       0.00         10          0          0        396        364          0
  disjunctive      :       0.00       0.00          0          0          0          0          0          0
  eccuts           :       0.00       0.00          0          0          0          0          0          0
  flowcover        :       0.14       0.00         15          0          0        168          4          0
  gomory           :       0.01       0.00         10          0          0         14         10          0
  impliedbounds    :       0.00       0.00         68          0          0        318         73          0
  intobj           :       0.00       0.00          0          0          0          0          0          0
  mcf              :       0.00       0.00          1          0          0          0          0          0
  oddcycle         :       0.00       0.00          0          0          0          0          0          0
  rapidlearning    :       0.00       0.00          0          0          0          0          0          0
  strongcg         :       0.04       0.00         20          0          0          0          0          0
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
  relpscost        :       0.93       0.00      28096          0          0          0       3827          0        -15      55272
Primal Heuristics  :   ExecTime  SetupTime      Calls      Found       Best
  LP solutions     :       0.00          -          -          0          -
  pseudo solutions :       0.00          -          -          0          -
  strong branching :       0.00          -          -          0          -
  SpaGreedy        :       0.42       0.00       1982          2          1
  actconsdiving    :       0.00       0.00          0          0          0
  bound            :       0.00       0.00          0          0          0
  clique           :       0.00       0.00          0          0          0
  coefdiving       :       0.61       0.00        196          0          0
  crossover        :       0.03       0.00          0          0          0
  dins             :       0.00       0.00          0          0          0
  distributiondivin:       0.81       0.00        196          0          0
  dualval          :       0.00       0.00          0          0          0
  feaspump         :       0.31       0.00         41          0          0
  fixandinfer      :       0.00       0.00          0          0          0
  fracdiving       :       0.56       0.00        197          0          0
  guideddiving     :       0.00       0.00          0          0          0
  indicator        :       0.00       0.00          0          0          0
  intdiving        :       0.00       0.00          0          0          0
  intshifting      :       0.04       0.00         14          0          0
  linesearchdiving :       0.58       0.00        196          0          0
  localbranching   :       0.00       0.00          0          0          0
  mutation         :       0.00       0.00          0          0          0
  nlpdiving        :       0.00       0.00          0          0          0
  objpscostdiving  :       0.23       0.00         25          0          0
  octane           :       0.00       0.00          0          0          0
  ofins            :       0.00       0.00          0          0          0
  oneopt           :       0.01       0.00          0          0          0
  proximity        :       0.00       0.00          0          0          0
  pscostdiving     :       0.50       0.00        197          0          0
  randrounding     :       0.12       0.00        292          0          0
  rens             :       0.00       0.00          0          0          0
  reoptsols        :       0.00       0.00          0          0          0
  rins             :       0.01       0.00          0          0          0
  rootsoldiving    :       0.16       0.00         26          0          0
  rounding         :       0.04       0.00       2281          0          0
  shiftandpropagate:       0.00       0.00          0          0          0
  shifting         :       0.12       0.00        437          0          0
  simplerounding   :       0.01       0.00          0          0          0
  subnlp           :       0.01       0.00          0          0          0
  trivial          :       0.00       0.00          2          0          0
  trivialnegation  :       0.00       0.00          0          0          0
  trysol           :       0.00       0.00          0          0          0
  twoopt           :       0.00       0.00          0          0          0
  undercover       :       0.00       0.00          1          0          0
  vbounds          :       0.00       0.00          0          0          0
  veclendiving     :       0.71       0.00        197          0          0
  zeroobj          :       0.00       0.00          0          0          0
  zirounding       :       0.01       0.00       1000          0          0
  other solutions  :          -          -          -          0          -
Diving Statistics  :      Calls      Nodes   LP Iters Backtracks   MinDepth   MaxDepth   AvgDepth  NLeafSols  MinSolDpt  MaxSolDpt  AvgSolDpt
  actconsdiving    :          0          -          -          -          -          -          -          -          -          -          -
  coefdiving       :        196       1799      16473        252         12         31       21.0          -          -          -          -
  distributiondivin:        196       1349      17947        253         13         30       18.8          -          -          -          -
  fracdiving       :        197       1624      19423        259         14         33       20.3          -          -          -          -
  guideddiving     :          0          -          -          -          -          -          -          -          -          -          -
  linesearchdiving :        196       1388      21864        303         14         28       18.6          -          -          -          -
  pscostdiving     :        197       1889      16783        218         14         35       21.7          -          -          -          -
  veclendiving     :        197       1658      19987        353         14         27       19.6          -          -          -          -
LP                 :       Time      Calls Iterations  Iter/call   Iter/sec  Time-0-It Calls-0-It
  primal LP        :       0.00         10          0       0.00          -       0.00         10
  dual LP          :      20.59      63973     937316      16.88   45522.88       0.70       8448
  lex dual LP      :       0.00          0          0       0.00          -
  barrier LP       :       0.00          0          0       0.00          -       0.00          0
  diving/probing LP:       3.05       6678     143763      21.53   47135.41
  strong branching :       0.64        804      34074      42.38   53240.62
    (at root node) :          -         19       1910     100.53          -
  conflict analysis:       0.00          0          0       0.00          -
B&B Tree           :
  number of runs   :          1
  nodes            :      55244 (27636 internal, 27608 leaves)
  nodes (total)    :      55244 (27636 internal, 27608 leaves)
  nodes left       :          0
  max depth        :         29
  max depth (total):         29
  backtracks       :      16122 (29.2%)
  delayed cutoffs  :         29
  repropagations   :       2482 (23199 domain reductions, 22 cutoffs)
  avg switch length:       4.83
  switching time   :       4.29
Root Node          :
  First LP value   : +2.00000000000000e+00
  First LP Iters   :        198
  First LP Time    :       0.00
  Final Dual Bound : +3.68148498333573e-03
  Final Root Iters :       2391
Solution           :
  Solutions found  :          2 (1 improvements)
  First Solution   : +2.88402607999991e-04   (in run 1, after 1 nodes, 0.04 seconds, depth 0, found by <SpaGreedy>)
  Gap First Sol.   :   infinite
  Gap Last Sol.    :   infinite
  Primal Bound     : +2.88402607999991e-04   (in run 1, after 1 nodes, 0.04 seconds, depth 0, found by <SpaGreedy>)
  Dual Bound       : +2.88402607999991e-04
  Gap              :       0.00 %
  Avg. Gap         :      92.26 % (3736.42 primal-dual integral)

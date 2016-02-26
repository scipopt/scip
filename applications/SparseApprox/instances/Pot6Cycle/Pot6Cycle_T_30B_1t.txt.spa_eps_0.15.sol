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
Total Time         :    3600.00
  solving          :    3600.00
  presolving       :       0.09 (included in solving)
  reading          :       0.00
  copying          :       0.00 (0 times copied the problem)
Original Problem   :
  Problem name     : ../instances/Pot6Cycle/Pot6Cycle_T_30B_1t.txt.spa
  Variables        : 141 (140 binary, 0 integer, 0 implicit integer, 1 continuous)
  Constraints      : 263 initial, 263 maximal
  Objective sense  : maximize
Presolved Problem  :
  Problem name     : t_../instances/Pot6Cycle/Pot6Cycle_T_30B_1t.txt.spa
  Variables        : 910 (132 binary, 0 integer, 16 implicit integer, 762 continuous)
  Constraints      : 3189 initial, 4743 maximal
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
  probing          :       0.06       0.00      1          0          0          0          0          0          0          0          0          0
  pseudoobj        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  varbound         :       0.00       0.00      9          0          0          0          0          0         28         28          0          0
  setppc           :       0.00       0.00     12          0          0          0          0          0          0          0          0          0
  and              :       0.00       0.00      4          0         12          0          0          0         12          0          0          0
  linear           :       0.01       0.00      9          3          1          0          4          0          7          0         12         12
  logicor          :       0.00       0.00      7          0          0          0          0          0          0          0          0          0
  bounddisjunction :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  quadratic        :       0.02       0.00      2          0          0          0          0          0          0       3019          0          0
  root node        :          -          -      -          1          -          -         69          -          -          -          -          -
Constraints        :     Number  MaxNumber  #Separate #Propagate    #EnfoLP    #EnfoPS     #Check   #ResProp    Cutoffs    DomReds       Cuts    Applied      Conss   Children
  integral         :          0          0          0          0    1892826          0     232596          0          0      41359          0          0          7    2630484
  varbound         :       1413       1413         75    7801625     574743          0      80849    1824313      23559   30427699          0          0          0          0
  setppc           :        160+       213         75    7777765     574743          0      34526    4398638       7916    2096432          8          0          0          0
  and              :         22         22       1855      25629          0          0      34445      99429          2       1246         23          9          0          0
  linear           :       1539       1539         75    7770148     574743          0      34379     905919      36977   71624698   41405970   41405647          0          0
  logicor          :         55+       536         64    1331431          0          0          0     960894        283      11086          1          1          0          0
  bounddisjunction :          0+      1096          0    2488739          0          0          0     183943       1226     905271          0          0          0          0
  countsols        :          0          0          0          0          0          0     151744          0          0          0          0          0          0          0
Constraint Timings :  TotalTime  SetupTime   Separate  Propagate     EnfoLP     EnfoPS      Check    ResProp    SB-Prop
  integral         :      17.60       0.00       0.00       0.00      17.54       0.00       0.06       0.00       0.00
  varbound         :      34.23       0.00       0.00      28.30       3.77       0.00       1.79       0.26       0.11
  setppc           :      31.29       0.00       0.00      30.16       0.07       0.00       0.03       0.78       0.25
  and              :       0.09       0.00       0.01       0.00       0.00       0.00       0.05       0.03       0.00
  linear           :     164.56       0.00       0.03     127.79      34.42       0.00       1.09       0.93       0.30
  logicor          :       5.13       0.00       0.00       4.68       0.00       0.00       0.00       0.45       0.00
  bounddisjunction :      19.83       0.00       0.00      19.78       0.00       0.00       0.00       0.04       0.01
  countsols        :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
Propagators        : #Propagate   #ResProp    Cutoffs    DomReds
  dualfix          :          2          0          0          0
  genvbounds       :          0          0          0          0
  obbt             :          0          0          0          0
  probing          :          0          0          0          0
  pseudoobj        :    8516353          0         12         68
  redcost          :    3962913          0        582     322823
  rootredcost      :          1          0          0          0
  vbounds          :    9651867        697          0        299
Propagator Timings :  TotalTime  SetupTime   Presolve  Propagate    ResProp    SB-Prop
  dualfix          :       0.00       0.00       0.00       0.00       0.00       0.00
  genvbounds       :       1.24       0.00       0.00       1.24       0.00       0.00
  obbt             :       0.00       0.00       0.00       0.00       0.00       0.00
  probing          :       0.06       0.00       0.06       0.00       0.00       0.00
  pseudoobj        :       2.09       0.00       0.00       2.09       0.00       0.00
  redcost          :      22.08       0.00       0.00      22.08       0.00       0.00
  rootredcost      :       1.18       0.00       0.00       1.17       0.00       0.01
  vbounds          :       2.71       0.00       0.00       2.70       0.00       0.01
Conflict Analysis  :       Time      Calls    Success    DomReds  Conflicts   Literals    Reconvs ReconvLits   LP Iters
  propagation      :       5.19      33075      28735          -     527037       63.9       1912       74.1          -
  infeasible LP    :      21.57     141160      70911          -    1132382       51.9         55       15.4          0
  bound exceed. LP :       0.00          0          0          -          0        0.0          0        0.0          0
  strong branching :       0.00          0          0          -          0        0.0          0        0.0          0
  pseudo solution  :       0.00          0          0          -          0        0.0          0        0.0          -
  applied globally :       1.35          -          -          0      97563       23.6          -          -          -
  applied locally  :          -          -          -          0      12124       41.0          -          -          -
Separators         :   ExecTime  SetupTime      Calls    Cutoffs    DomReds       Cuts    Applied      Conss
  cut pool         :       0.00                    73          -          -        769          -          -    (maximal pool size: 613)
  cgmip            :       0.00       0.00          0          0          0          0          0          0
  clique           :       0.03       0.00         74          0          0         47         45          0
  closecuts        :       0.00       0.00          0          0          0          0          0          0
  cmir             :       0.02       0.00         10          0          0         62         13          0
  disjunctive      :       0.00       0.00          0          0          0          0          0          0
  eccuts           :       0.00       0.00          0          0          0          0          0          0
  flowcover        :       0.32       0.00         15          0          0        244         54          0
  gomory           :       0.06       0.00         10          0          0         34         21          0
  impliedbounds    :       0.01       0.00         74          0          0        572        146          0
  intobj           :       0.00       0.00          0          0          0          0          0          0
  mcf              :       0.00       0.00          1          0          0          0          0          0
  oddcycle         :       0.00       0.00          0          0          0          0          0          0
  rapidlearning    :       0.00       0.00          0          0          0          0          0          0
  strongcg         :       0.15       0.00         20          0          0          3          2          0
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
  relpscost        :      16.96       0.00    1318083          0          0          0      41359          0       -652    2630484
Primal Heuristics  :   ExecTime  SetupTime      Calls      Found       Best
  LP solutions     :       0.00          -          -          0          -
  pseudo solutions :       0.00          -          -          0          -
  strong branching :       0.02          -          -          0          -
  SpaGreedy        :     610.32       0.00     151743          2          1
  actconsdiving    :       0.00       0.00          0          0          0
  bound            :       0.00       0.00          0          0          0
  clique           :       0.00       0.00          0          0          0
  coefdiving       :      26.95       0.00       8113          0          0
  crossover        :       0.26       0.00          0          0          0
  dins             :       0.00       0.00          0          0          0
  distributiondivin:      57.68       0.00       8113          0          0
  dualval          :       0.00       0.00          0          0          0
  feaspump         :       5.94       0.00       4057          0          0
  fixandinfer      :       0.00       0.00          0          0          0
  fracdiving       :      29.10       0.00       8113          0          0
  guideddiving     :       0.05       0.00          0          0          0
  indicator        :       0.27       0.00          0          0          0
  intdiving        :       0.00       0.00          0          0          0
  intshifting      :       0.27       0.00         99          0          0
  linesearchdiving :      29.60       0.00       8113          0          0
  localbranching   :       0.00       0.00          0          0          0
  mutation         :       0.00       0.00          0          0          0
  nlpdiving        :       0.00       0.00          0          0          0
  objpscostdiving  :      13.26       0.00       4057          0          0
  octane           :       0.00       0.00          0          0          0
  ofins            :       0.00       0.00          0          0          0
  oneopt           :       0.42       0.00          0          0          0
  proximity        :       0.00       0.00          0          0          0
  pscostdiving     :      31.13       0.00       8113          0          0
  randrounding     :       1.09       0.00      45386          0          0
  rens             :       0.00       0.00          0          0          0
  reoptsols        :       0.00       0.00          0          0          0
  rins             :       0.39       0.00          0          0          0
  rootsoldiving    :      10.67       0.00       4057          0          0
  rounding         :       0.59       0.00      16104          0          0
  shiftandpropagate:       0.00       0.00          0          0          0
  shifting         :       0.85       0.00       2616          0          0
  simplerounding   :       0.26       0.00          0          0          0
  subnlp           :       0.36       0.00          0          0          0
  trivial          :       0.00       0.00          2          0          0
  trivialnegation  :       0.00       0.00          0          0          0
  trysol           :       0.50       0.00          0          0          0
  twoopt           :       0.00       0.00          0          0          0
  undercover       :       0.00       0.00          1          0          0
  vbounds          :       0.00       0.00          0          0          0
  veclendiving     :      29.76       0.00       8113          0          0
  zeroobj          :       0.00       0.00          0          0          0
  zirounding       :       0.36       0.00       1000          0          0
  other solutions  :          -          -          -          0          -
Diving Statistics  :      Calls      Nodes   LP Iters Backtracks   MinDepth   MaxDepth   AvgDepth  NLeafSols  MinSolDpt  MaxSolDpt  AvgSolDpt
  actconsdiving    :          0          -          -          -          -          -          -          -          -          -          -
  coefdiving       :       8113      66078     718391       4553         20         68       30.3          -          -          -          -
  distributiondivin:       8113      62242     767080       5761         23         52       29.8          -          -          -          -
  fracdiving       :       8113      65091     800714       4803         23         55       30.2          -          -          -          -
  guideddiving     :          0          -          -          -          -          -          -          -          -          -          -
  linesearchdiving :       8113      65473     789195       4544         24         51       30.1          -          -          -          -
  pscostdiving     :       8113      66141     750799       5405         21         69       30.2          -          -          -          -
  veclendiving     :       8113      68859     785994       6786         24         48       30.4          -          -          -          -
LP                 :       Time      Calls Iterations  Iter/call   Iter/sec  Time-0-It Calls-0-It
  primal LP        :       0.00         11          0       0.00          -       0.00         11
  dual LP          :    1575.36    3178074   67821842      23.87   43051.65      45.51     336571
  lex dual LP      :       0.00          0          0       0.00          -
  barrier LP       :       0.00          0          0       0.00          -       0.00          0
  diving/probing LP:     124.11     187565    5180352      27.62   41740.00
  strong branching :       9.69       4077     214873      52.70   22174.72
    (at root node) :          -         20       2273     113.65          -
  conflict analysis:       0.00          0          0       0.00          -
B&B Tree           :
  number of runs   :          1
  nodes            :    2351802 (1315242 internal, 1036560 leaves)
  nodes (total)    :    2351802 (1315242 internal, 1036560 leaves)
  nodes left       :     275391
  max depth        :         33
  max depth (total):         33
  backtracks       :     578431 (24.6%)
  delayed cutoffs  :       3292
  repropagations   :      63975 (922893 domain reductions, 2778 cutoffs)
  avg switch length:       5.00
  switching time   :     398.75
Root Node          :
  First LP value   : +3.00000000000000e+00
  First LP Iters   :        617
  First LP Time    :       0.00
  Final Dual Bound : +1.09496686887221e-02
  Final Root Iters :       4655
Solution           :
  Solutions found  :          2 (1 improvements)
  First Solution   : +1.03971965299997e-03   (in run 1, after 1 nodes, 0.09 seconds, depth 0, found by <SpaGreedy>)
  Gap First Sol.   :   infinite
  Gap Last Sol.    :   infinite
  Primal Bound     : +1.03971965299997e-03   (in run 1, after 1 nodes, 0.09 seconds, depth 0, found by <SpaGreedy>)
  Dual Bound       : +2.74356097724507e-03
  Gap              :     163.88 %
  Avg. Gap         :      90.51 % (325828.20 primal-dual integral)

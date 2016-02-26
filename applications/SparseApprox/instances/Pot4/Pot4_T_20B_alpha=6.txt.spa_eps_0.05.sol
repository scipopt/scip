objective value:                  0.00077680452900003
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
x_14_2                                              1 	(obj:0)
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
epsI                              0.00077680452900003 	(obj:1)
SCIP Status        : solving was interrupted [node limit reached]
Total Time         :       3.57
  solving          :       3.57
  presolving       :       0.10 (included in solving)
  reading          :       0.00
  copying          :       0.00 (0 times copied the problem)
Original Problem   :
  Problem name     : ../instances/Pot4/Pot4_T_20B_alpha=6.txt.spa
  Variables        : 101 (100 binary, 0 integer, 0 implicit integer, 1 continuous)
  Constraints      : 193 initial, 193 maximal
  Objective sense  : maximize
Presolved Problem  :
  Problem name     : t_../instances/Pot4/Pot4_T_20B_alpha=6.txt.spa
  Variables        : 590 (92 binary, 0 integer, 16 implicit integer, 482 continuous)
  Constraints      : 2015 initial, 2051 maximal
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
  logicor          :       0.00       0.01      7          0          0          0          0          0          0          0          0          0
  bounddisjunction :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  quadratic        :       0.01       0.00      2          0          0          0          0          0          0       1915          0          0
  root node        :          -          -      -          1          -          -         57          -          -          -          -          -
Constraints        :     Number  MaxNumber  #Separate #Propagate    #EnfoLP    #EnfoPS     #Check   #ResProp    Cutoffs    DomReds       Cuts    Applied      Conss   Children
  integral         :          0          0          0          0        237          0         23          0          0       1311          0          0          0        172
  varbound         :        869        869         86       5229          0          0         11        677         27       2053          0          0          0          0
  setppc           :        110        110         86       5202          0          0          5        611         15        184          8          1          0          0
  and              :         22         22        135       4079          0          0          3         54          3         46         28          8          0          0
  linear           :        979        979         86       5184          0          0          1        410         49       2290        473        225          0          0
  logicor          :         35+        43         85       1274          0          0          0         51          3         14          2          2          0          0
  bounddisjunction :          0+        48          0       1718          0          0          0          0          0          1          0          0          0          0
  countsols        :          0          0          0          0          0          0         11          0          0          0          0          0          0          0
Constraint Timings :  TotalTime  SetupTime   Separate  Propagate     EnfoLP     EnfoPS      Check    ResProp    SB-Prop
  integral         :       1.86       0.00       0.00       0.00       1.86       0.00       0.00       0.00       0.00
  varbound         :       0.14       0.00       0.01       0.10       0.00       0.00       0.00       0.00       0.03
  setppc           :       0.09       0.00       0.00       0.08       0.00       0.00       0.00       0.00       0.01
  and              :       0.01       0.00       0.00       0.01       0.00       0.00       0.00       0.00       0.00
  linear           :       0.28       0.00       0.01       0.23       0.00       0.00       0.00       0.00       0.04
  logicor          :       0.01       0.01       0.00       0.00       0.00       0.00       0.00       0.00       0.00
  bounddisjunction :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
  countsols        :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
Propagators        : #Propagate   #ResProp    Cutoffs    DomReds
  dualfix          :          2          0          0          0
  genvbounds       :          0          0          0          0
  obbt             :          0          0          0          0
  probing          :          0          0          0          0
  pseudoobj        :        853          0          0         56
  redcost          :        561          0          0          1
  rootredcost      :          1          0          0          0
  vbounds          :       5070          0          0          0
Propagator Timings :  TotalTime  SetupTime   Presolve  Propagate    ResProp    SB-Prop
  dualfix          :       0.00       0.00       0.00       0.00       0.00       0.00
  genvbounds       :       0.00       0.00       0.00       0.00       0.00       0.00
  obbt             :       0.00       0.00       0.00       0.00       0.00       0.00
  probing          :       0.05       0.00       0.05       0.00       0.00       0.00
  pseudoobj        :       0.00       0.00       0.00       0.00       0.00       0.00
  redcost          :       0.01       0.00       0.00       0.01       0.00       0.00
  rootredcost      :       0.01       0.00       0.00       0.01       0.00       0.00
  vbounds          :       0.01       0.00       0.00       0.01       0.00       0.00
Conflict Analysis  :       Time      Calls    Success    DomReds  Conflicts   Literals    Reconvs ReconvLits   LP Iters
  propagation      :       0.00         48         48          -        623       39.7         10       38.9          -
  infeasible LP    :       0.00          5          0          -          0        0.0          0        0.0          0
  bound exceed. LP :       0.00          0          0          -          0        0.0          0        0.0          0
  strong branching :       0.00          0          0          -          0        0.0          0        0.0          0
  pseudo solution  :       0.00          0          0          -          0        0.0          0        0.0          -
  applied globally :       0.00          -          -          0         59       17.6          -          -          -
  applied locally  :          -          -          -          0          4       32.2          -          -          -
Separators         :   ExecTime  SetupTime      Calls    Cutoffs    DomReds       Cuts    Applied      Conss
  cut pool         :       0.01                    85          -          -       1152          -          -    (maximal pool size: 782)
  cgmip            :       0.00       0.00          0          0          0          0          0          0
  clique           :       0.04       0.00         86          0          0         35         30          0
  closecuts        :       0.00       0.00          0          0          0          0          0          0
  cmir             :       0.03       0.00         10          0          0        197        124          0
  disjunctive      :       0.00       0.00          0          0          0          0          0          0
  eccuts           :       0.00       0.00          0          0          0          0          0          0
  flowcover        :       0.21       0.00         15          0          0        240         27          0
  gomory           :       0.04       0.00         10          0          0          6          6          0
  impliedbounds    :       0.01       0.00         86          0          0        707        146          0
  intobj           :       0.00       0.00          0          0          0          0          0          0
  mcf              :       0.00       0.00          1          0          0          0          0          0
  oddcycle         :       0.00       0.00          0          0          0          0          0          0
  rapidlearning    :       0.00       0.00          0          0          0          0          0          0
  strongcg         :       0.07       0.00         20          0          0         82         22          0
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
  relpscost        :       1.86       0.00        237          0          0          0       1311          0          0        172
Primal Heuristics  :   ExecTime  SetupTime      Calls      Found       Best
  LP solutions     :       0.00          -          -          0          -
  pseudo solutions :       0.00          -          -          0          -
  strong branching :       0.00          -          -          0          -
  SpaGreedy        :       0.01       0.00         10          2          1
  actconsdiving    :       0.00       0.00          0          0          0
  bound            :       0.00       0.00          0          0          0
  clique           :       0.00       0.00          0          0          0
  coefdiving       :       0.02       0.00          2          0          0
  crossover        :       0.00       0.00          0          0          0
  dins             :       0.00       0.00          0          0          0
  distributiondivin:       0.00       0.00          1          0          0
  dualval          :       0.00       0.00          0          0          0
  feaspump         :       0.01       0.00          1          0          0
  fixandinfer      :       0.00       0.00          0          0          0
  fracdiving       :       0.00       0.00          1          0          0
  guideddiving     :       0.00       0.00          0          0          0
  indicator        :       0.00       0.00          0          0          0
  intdiving        :       0.00       0.00          0          0          0
  intshifting      :       0.02       0.00          2          0          0
  linesearchdiving :       0.01       0.00          1          0          0
  localbranching   :       0.00       0.00          0          0          0
  mutation         :       0.00       0.00          0          0          0
  nlpdiving        :       0.00       0.00          0          0          0
  objpscostdiving  :       0.01       0.00          1          0          0
  octane           :       0.00       0.00          0          0          0
  ofins            :       0.00       0.00          0          0          0
  oneopt           :       0.00       0.00          0          0          0
  proximity        :       0.00       0.00          0          0          0
  pscostdiving     :       0.00       0.00          1          0          0
  randrounding     :       0.61       0.00         86          0          0
  rens             :       0.00       0.00          0          0          0
  reoptsols        :       0.00       0.00          0          0          0
  rins             :       0.00       0.00          0          0          0
  rootsoldiving    :       0.00       0.00          1          0          0
  rounding         :       0.03       0.00        187          0          0
  shiftandpropagate:       0.00       0.00          0          0          0
  shifting         :       0.04       0.00         94          0          0
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
  zirounding       :       0.00       0.00         85          0          0
  other solutions  :          -          -          -          0          -
Diving Statistics  :      Calls      Nodes   LP Iters Backtracks   MinDepth   MaxDepth   AvgDepth  NLeafSols  MinSolDpt  MaxSolDpt  AvgSolDpt
  actconsdiving    :          0          -          -          -          -          -          -          -          -          -          -
  coefdiving       :          2         31        415          4         28         32       30.0          -          -          -          -
  distributiondivin:          1          3          0          1         20         20       20.0          -          -          -          -
  fracdiving       :          1          3          0          1         20         20       20.0          -          -          -          -
  guideddiving     :          0          -          -          -          -          -          -          -          -          -          -
  linesearchdiving :          1          5         38          0         21         21       21.0          -          -          -          -
  pscostdiving     :          1          3          0          1         20         20       20.0          -          -          -          -
  veclendiving     :          1          4          0          1         21         21       21.0          -          -          -          -
LP                 :       Time      Calls Iterations  Iter/call   Iter/sec  Time-0-It Calls-0-It
  primal LP        :       0.01         14          0       0.00       0.00       0.01         14
  dual LP          :       0.32        392       7108      26.13   22212.50       0.05        120
  lex dual LP      :       0.00          0          0       0.00          -
  barrier LP       :       0.00          0          0       0.00          -       0.00          0
  diving/probing LP:       0.09         52       1185      22.79   13166.67
  strong branching :       1.54        332      31051      93.53   20162.99
    (at root node) :          -         24       1516      63.17          -
  conflict analysis:       0.00          0          0       0.00          -
B&B Tree           :
  number of runs   :          1
  nodes            :        100 (86 internal, 14 leaves)
  nodes (total)    :        100 (86 internal, 14 leaves)
  nodes left       :         73
  max depth        :         20
  max depth (total):         20
  backtracks       :         11 (11.0%)
  delayed cutoffs  :          0
  repropagations   :          2 (0 domain reductions, 0 cutoffs)
  avg switch length:       2.14
  switching time   :       0.01
Root Node          :
  First LP value   : +3.00000000000000e+00
  First LP Iters   :        375
  First LP Time    :       0.00
  Final Dual Bound : +7.05325623049629e-03
  Final Root Iters :       3452
Solution           :
  Solutions found  :          2 (1 improvements)
  First Solution   : +7.76804529000030e-04   (in run 1, after 1 nodes, 0.10 seconds, depth 0, found by <SpaGreedy>)
  Gap First Sol.   :   infinite
  Gap Last Sol.    :   infinite
  Primal Bound     : +7.76804529000030e-04   (in run 1, after 1 nodes, 0.10 seconds, depth 0, found by <SpaGreedy>)
  Dual Bound       : +6.51541054660569e-03
  Gap              :     738.75 %
  Avg. Gap         :      92.52 % (330.31 primal-dual integral)

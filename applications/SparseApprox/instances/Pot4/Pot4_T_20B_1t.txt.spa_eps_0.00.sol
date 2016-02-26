objective value:                  0.00061994677600001
x_1_2                                               1 	(obj:0)
x_2_1                                               1 	(obj:0)
x_3_3                                               1 	(obj:0)
x_4_1                                               1 	(obj:0)
x_5_2                                               1 	(obj:0)
x_6_3                                               1 	(obj:0)
x_7_1                                               1 	(obj:0)
x_8_2                                               1 	(obj:0)
x_9_2                                               1 	(obj:0)
x_10_3                                              1 	(obj:0)
x_11_1                                              1 	(obj:0)
x_12_3                                              1 	(obj:0)
x_13_1                                              1 	(obj:0)
x_14_3                                              1 	(obj:0)
x_15_2                                              1 	(obj:0)
x_16_2                                              1 	(obj:0)
x_17_1                                              1 	(obj:0)
x_18_1                                              1 	(obj:0)
x_19_3                                              1 	(obj:0)
x_20_2                                              1 	(obj:0)
ind_1                                               1 	(obj:0)
abs_1_2                                             1 	(obj:0)
ind_2                                               1 	(obj:0)
abs_2_3                                             1 	(obj:0)
ind_3                                               1 	(obj:0)
abs_3_1                                             1 	(obj:0)
epsI                              0.00061994677600001 	(obj:1)
SCIP Status        : solving was interrupted [user interrupt]
Total Time         :     331.61
  solving          :     331.60
  presolving       :       0.16 (included in solving)
  reading          :       0.01
  copying          :       0.15 (8 #copies) (minimal 0.01, maximal 0.03, average 0.02)
Original Problem   :
  Problem name     : instances/Pot4/Pot4_T_20B_1t.txt.spa
  Variables        : 73 (72 binary, 0 integer, 0 implicit integer, 1 continuous)
  Constraints      : 128 initial, 128 maximal
  Objective sense  : maximize
Presolved Problem  :
  Problem name     : t_instances/Pot4/Pot4_T_20B_1t.txt.spa
  Variables        : 340 (67 binary, 0 integer, 9 implicit integer, 264 continuous)
  Constraints      : 1124 initial, 1910 maximal
Presolvers         :   ExecTime  SetupTime  Calls  FixedVars   AggrVars   ChgTypes  ChgBounds   AddHoles    DelCons    AddCons   ChgSides   ChgCoefs
  boundshift       :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  components       :       0.00       0.00      1          0          0          0          0          0          0          0          0          0
  convertinttobin  :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  domcol           :       0.00       0.00      1          0          0          0          0          0          0          0          0          0
  dualagg          :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  dualinfer        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  gateextraction   :       0.00       0.00      3          0          0          0          0          0         57         18          0          0
  implfree         :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  implics          :       0.00       0.00      4          0          0          0          0          0          0          0          0          0
  inttobinary      :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  redvub           :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  stuffing         :       0.01       0.00      1          0          0          0          0          0          0          0          0          0
  trivial          :       0.00       0.00      7          0          0          0          0          0          0          0          0          0
  tworowbnd        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  dualfix          :       0.00       0.00      7          3          0          0          0          0          0          0          0          0
  genvbounds       :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  probing          :       0.08       0.00      1          0          0          0          0          0          0          0          0          0
  pseudoobj        :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  varbound         :       0.00       0.00      9          0          0          0          0          0         15         15          0          0
  setppc           :       0.01       0.00     12          0          0          0          0          0          0          0          0          0
  and              :       0.00       0.00      4          0          6          0          0          0          6          0          0          0
  linear           :       0.04       0.00      9          1          1          0          2          0          3          0          6          6
  logicor          :       0.00       0.00      7          0          0          0          0          0          0          0          0          0
  bounddisjunction :       0.00       0.00      0          0          0          0          0          0          0          0          0          0
  quadratic        :       0.02       0.00      2          0          0          0          0          0          0       1044          0          0
  root node        :          -          -      -          1          -          -         54          -          -          -          -          -
Constraints        :     Number  MaxNumber  #Separate #Propagate    #EnfoLP    #EnfoPS     #Check   #ResProp    Cutoffs    DomReds       Cuts    Applied      Conss   Children
  integral         :          0          0          0          0     124977          0       1068          0          0       5153          0          0          4     248944
  varbound         :        473        473         85     802003         20          0       1035      76197        796    1987598          0          0          0          0
  setppc           :         85+        86         85     800805         20          0         41      65404       1242     238041          0          0          0          0
  and              :         12         12        772       7438          0          0         30       5552          2        312         21          6          0          0
  linear           :        536        536         85     800001         20          0         17      45782       3853    4734152        590        435          0          0
  logicor          :         18+       295         85     191695          0          0          0        243         42       1476          2          2          0          0
  bounddisjunction :          0+       514          0     197915          0          0          0       2952          9      38070          0          0          0          0
  countsols        :          0          0          0          0          7          0         14          0          0          0          0          0          0          0
Constraint Timings :  TotalTime  SetupTime   Separate  Propagate     EnfoLP     EnfoPS      Check    ResProp    SB-Prop
  integral         :       4.79       0.00       0.00       0.00       4.79       0.00       0.00       0.00       0.00
  varbound         :      12.23       0.00       0.01      12.08       0.00       0.00       0.00       0.05       0.09
  setppc           :       7.75       0.00       0.00       7.63       0.00       0.00       0.00       0.09       0.03
  and              :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
  linear           :      63.96       0.00       0.03      63.48       0.01       0.00       0.00       0.18       0.26
  logicor          :       0.87       0.00       0.23       0.64       0.00       0.00       0.00       0.00       0.00
  bounddisjunction :       1.90       0.00       0.00       1.90       0.00       0.00       0.00       0.00       0.00
  countsols        :       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00       0.00
Propagators        : #Propagate   #ResProp    Cutoffs    DomReds
  dualfix          :          1          0          0          0
  genvbounds       :          0          0          0          0
  obbt             :          0          0          0          0
  probing          :          0          0          0          0
  pseudoobj        :     774095          0        110         50
  redcost          :     253145          0          0       3810
  rootredcost      :         17          0          0          4
  vbounds          :     926751        239          0        131
Propagator Timings :  TotalTime  SetupTime   Presolve  Propagate    ResProp    SB-Prop
  dualfix          :       0.00       0.00       0.00       0.00       0.00       0.00
  genvbounds       :       0.19       0.00       0.00       0.19       0.00       0.00
  obbt             :       0.00       0.00       0.00       0.00       0.00       0.00
  probing          :       0.08       0.00       0.08       0.00       0.00       0.00
  pseudoobj        :       0.61       0.00       0.00       0.61       0.00       0.00
  redcost          :       1.30       0.00       0.00       1.30       0.00       0.00
  rootredcost      :       0.44       0.00       0.00       0.44       0.00       0.00
  vbounds          :       0.40       0.00       0.00       0.40       0.00       0.00
Conflict Analysis  :       Time      Calls    Success    DomReds  Conflicts   Literals    Reconvs ReconvLits   LP Iters
  propagation      :       0.58       2227       1743          -      21806       44.8        138       47.3          -
  infeasible LP    :       0.65       1900       1428          -      18868       35.2          1       15.0          0
  bound exceed. LP :       0.00          0          0          -          0        0.0          0        0.0          0
  strong branching :       0.00          0          0          -          0        0.0          0        0.0          0
  pseudo solution  :       0.00          0          0          -          0        0.0          0        0.0          -
  applied globally :       0.17          -          -          0       4018       21.0          -          -          -
  applied locally  :          -          -          -          0        361       35.1          -          -          -
Separators         :   ExecTime  SetupTime      Calls    Cutoffs    DomReds       Cuts    Applied      Conss
  cut pool         :       0.01                    83          -          -       1577          -          -    (maximal pool size: 1049)
  cgmip            :       0.00       0.00          0          0          0          0          0          0
  clique           :       0.04       0.00         84          0          0        120         16          0
  closecuts        :       0.00       0.00          0          0          0          0          0          0
  cmir             :       0.13       0.00         10          0          0        594        413          0
  disjunctive      :       0.00       0.00          0          0          0          0          0          0
  eccuts           :       0.00       0.00          0          0          0          0          0          0
  flowcover        :       1.08       0.00         15          0          0        486         17          0
  gomory           :       0.13       0.00         10          0          0          7          7          0
  impliedbounds    :       0.01       0.00         84          0          0        781        150          0
  intobj           :       0.00       0.00          0          0          0          0          0          0
  mcf              :       0.00       0.00          1          0          0          0          0          0
  oddcycle         :       0.00       0.00          0          0          0          0          0          0
  rapidlearning    :       0.00       0.00          0          0          0          0          0          0
  strongcg         :       0.27       0.00         20          0          0          3          1          0
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
  relpscost        :       4.65       0.00     124950          0          0          0       5153          0         -2     248944
Primal Heuristics  :   ExecTime  SetupTime      Calls      Found       Best
  LP solutions     :       0.00          -          -          7          -
  pseudo solutions :       0.00          -          -          0          -
  strong branching :       0.00          -          -          0          -
  SpaGreedy        :       0.00       0.00          0          0          0
  actconsdiving    :       0.00       0.00          0          0          0
  bound            :       0.00       0.00          0          0          0
  clique           :       0.00       0.00          0          0          0
  coefdiving       :       9.71       0.00        993          1          1
  crossover        :       6.77       0.00          8          0          0
  dins             :       0.00       0.00          0          0          0
  distributiondivin:      15.81       0.00        991          0          0
  dualval          :       0.00       0.00          0          0          0
  feaspump         :       0.02       0.00          2          0          0
  fixandinfer      :       0.00       0.00          0          0          0
  fracdiving       :       7.01       0.00        907          0          0
  guideddiving     :       0.00       0.00          0          0          0
  indicator        :       0.01       0.00          0          0          0
  intdiving        :       0.00       0.00          0          0          0
  intshifting      :       0.02       0.00         14          0          0
  linesearchdiving :       7.29       0.00        932          0          0
  localbranching   :       0.00       0.00          0          0          0
  mutation         :       0.00       0.00          0          0          0
  nlpdiving        :       0.00       0.00          0          0          0
  objpscostdiving  :       0.91       0.00         59          0          0
  octane           :       0.00       0.00          0          0          0
  ofins            :       0.00       0.00          0          0          0
  oneopt           :       0.11       0.00          7          0          0
  proximity        :       0.00       0.00          0          0          0
  pscostdiving     :       8.14       0.00        993          2          2
  randrounding     :       0.96       0.00       1040          0          0
  rens             :       0.00       0.00          0          0          0
  reoptsols        :       0.00       0.00          0          0          0
  rins             :       0.12       0.00          0          0          0
  rootsoldiving    :       0.89       0.00         97          0          0
  rounding         :       0.59       0.00       4870          0          0
  shiftandpropagate:       0.01       0.00          1          1          1
  shifting         :       0.74       0.00        523          0          0
  simplerounding   :       0.07       0.00          0          0          0
  spaswitch        :       0.07       0.00         19          9          9
  subnlp           :       0.01       0.00          0          0          0
  trivial          :       0.00       0.00          2          0          0
  trivialnegation  :       0.00       0.00          0          0          0
  trysol           :       0.04       0.00          0          0          0
  twoopt           :       0.00       0.00          0          0          0
  undercover       :       0.00       0.00          1          0          0
  vbounds          :       0.00       0.00          0          0          0
  veclendiving     :       8.34       0.00        976          0          0
  zeroobj          :       0.00       0.00          0          0          0
  zirounding       :       0.20       0.00       1000          0          0
  other solutions  :          -          -          -          0          -
Diving Statistics  :      Calls      Nodes   LP Iters Backtracks   MinDepth   MaxDepth   AvgDepth  NLeafSols  MinSolDpt  MaxSolDpt  AvgSolDpt
  actconsdiving    :          0          -          -          -          -          -          -          -          -          -          -
  coefdiving       :        993      11444     121869       1200         16         37       25.5          -          -          -          -
  distributiondivin:        991       8297     110740       1217         15         36       22.9          -          -          -          -
  fracdiving       :        907       8447     117967       1256         16         37       23.5          -          -          -          -
  guideddiving     :          0          -          -          -          -          -          -          -          -          -          -
  linesearchdiving :        932       7569     118112       1347         16         36       22.3          -          -          -          -
  pscostdiving     :        993      11813     131188       1125         17         43       26.6          -          -          -          -
  veclendiving     :        976       8878     117931       1424         16         35       23.5          -          -          -          -
LP                 :       Time      Calls Iterations  Iter/call   Iter/sec  Time-0-It Calls-0-It
  primal LP        :       0.00         13          0       0.00          -       0.00         13
  dual LP          :      54.34     244196    2312586       9.79   42557.71       0.84       8027
  lex dual LP      :       0.00          0          0       0.00          -
  barrier LP       :       0.00          0          0       0.00          -       0.00          0
  diving/probing LP:      16.36      31314     773374      24.70   47272.25
  strong branching :       0.78        854      36757      43.04   47124.36
    (at root node) :          -         19       1189      62.58          -
  conflict analysis:       0.00          0          0       0.00          -
B&B Tree           :
  number of runs   :          1
  nodes            :     238539 (124472 internal, 114067 leaves)
  nodes (total)    :     238539 (124472 internal, 114067 leaves)
  nodes left       :      10353
  max depth        :         34
  max depth (total):         34
  backtracks       :      64639 (27.1%)
  delayed cutoffs  :         51
  repropagations   :       1781 (13707 domain reductions, 45 cutoffs)
  avg switch length:       4.41
  switching time   :      65.99
Root Node          :
  First LP value   : +2.00000000000000e+00
  First LP Iters   :        211
  First LP Time    :       0.00
  Final Dual Bound : +1.62415875355763e-01
  Final Root Iters :       4166
Solution           :
  Solutions found  :         20 (20 improvements)
  First Solution   : +0.00000000000000e+00   (in run 1, after 1 nodes, 0.17 seconds, depth 17, found by <shiftandpropagate>)
  Gap First Sol.   :   infinite
  Gap Last Sol.    :      76.77 %
  Primal Bound     : +6.19946776000010e-04   (in run 1, after 99784 nodes, 124.80 seconds, depth 16, found by <spaswitch>)
  Dual Bound       : +1.03648230094114e-03
  Gap              :      67.19 %
  Avg. Gap         :      45.57 % (15112.58 primal-dual integral)

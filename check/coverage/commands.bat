disp branch
disp conflict
disp cons
disp
disp
disp heur
disp nlpis
disp nodesel
disp param
disp presol
disp pricers
disp problem
disp dualsol
disp prop
disp reader
disp relax
disp sepa
disp sols
disp solution
disp stat
disp trans
transproblem
disp transsolution
disp value
disp varbranchstatistics
disp ..
set emph counter
set emph cpsolver
set emph easycip
set emph feasibility
set emph hardlp
set emph optimality
set heur emph aggr
set heur emph fast
set heur emph off
set sepa emph aggr
set sepa emph fast
set sepa emph off
set presol emph aggr
set presol emph fast
set presol emph off
set misc printreason 0
set misc printreason T
set diffsave tmp.set
set save tmp.set
set default
set load tmp.set
set default
read check/instances/MIP/bell5.mps
set lim no 10000
count
write allsolutions
read check/instances/MIP/bell5.mps
set cons countsol collect T
set lim no 10000
count
write allsolutions
disp lpsolqual
read check/instances/MIP/bell5.mps
opt
disp solution
disp sols 1
disp sols -1
disp val c1
help
free
read check/instances/MIP/bell5.mps
set lim no 1
set lim obj 1e+07
o
set lim obj 1e+08
write lp tmp.lp
write mip tmp.lp
write nlp tmp.lp
write problem tmp.lp
write genproblem tmp.lp
write transproblem tmp.lp
write gentransproblem tmp.lp
newstart
set lim obj 9e+06
set lp pricing d
set lp pricing
d
fix lp pricing T
set lp pricing s
fix lp pricing F
set branch prio
c1
1
set branch direction
c1
1
set branch direction c1 -1
o
read check/instances/MIP/bell5.mps
change bounds c1
1
1
change add [linear]\ <new>:\ <c1>\ >=\ 1;
change objsense
max
change objsense
min
presolve
disp mem
write cliquegraph tmp.gml
disp problem
disp stat
change freetrans
change minuc
q
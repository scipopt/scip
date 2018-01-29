#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#@file    check.awk
#@brief   SCIP Check Report Generator
#@author  Thorsten Koch
#@author  Tobias Achterberg
#@author  Alexander Martin
#@author  Timo Berthold
#@author  Robert Waniek
#@author  Gregor Hendel
#@author  Marc Pfetsch
#
function abs(x)
{
   return x < 0 ? -x : x;
}
function min(x,y)
{
   return (x) < (y) ? (x) : (y);
}
function max(x,y)
{
   return (x) > (y) ? (x) : (y);
}

# check if a limit was reached that counts as timeout
function isLimitReached()
{
   if( timeout || gapreached || sollimitreached || memlimitreached || nodelimitreached )
      return 1;
   else
      return 0;
}
# set a limit status and increase the timers and counters accordingly
function setStatusToLimit()
{
   if( timeout )
      status = "timeout";
   else if( gapreached )
      status = "gaplimit";
   else if( sollimitreached )
      status = "sollimit";
   else if( memlimitreached )
      status = "memlimit";
   else if( nodelimitreached )
      status = "nodelimit";
   else
      printf("Warning: Trying to set limit status for problem %s, but no limit was reached.\n", prob);

   timeouttime += tottime;
   timeouts++;
}

# set a file status and increase timers and counters accordingly
function setStatusToFail(statusstr)
{
   status = statusstr;
   failtime += tottime;
   fail++;
}

# is the dual bound better than the optimal or best known solution according to the given solu file?
# 'Better' means larger for minimization problems, else 'smaller'.
function isDualBoundBetter()
{
   if( !(prob in sol) )
      return 0;
   # objective sense of 1 means minimization
   if( (objsense == 1 && db-sol[prob] > reltol) || ( objsense == -1 && sol[prob]-db > reltol) )
      return 1;
   else
      return 0;
}

# is the primal bound better than the optimal or best known solution according to the given solu file?
# 'Better' means smaller for minimization problems, else 'larger'.
function isPrimalBoundBetter()
{
   if( !(prob in sol) )
      return 0;
   # objective sense of 1 means minimization
   if( (objsense == 1 && sol[prob] - pb > reltol) || (objsense == -1 && pb - sol[prob] > reltol) )
      return 1;
   else
      return 0;
}

# are primal and dual bound equal with respect to tolerances?
function isPrimalDualBoundEqual()
{
   if( abs(pb - db) <= max(abstol, reltol) )
      return 1;
   else
      return 0;
}

function isPrimalBoundBetterThanBestDual()
{
   if( !(prob in bestdual) )
      return 0;
   # objective sense of 1 means minimization
   if( (objsense == 1 && bestdual[prob] - pb > reltol) || (objsense == -1 && pb - bestdual[prob] > reltol) )
      return 1;
   return 0;
}

# write status and bound to solu file in case that the status requires a bound
function writeToSoluFile(status, thebound)
{
   # write status and problem name
   printf("=%s= %-18s",status, prob) > NEWSOLUFILE;

   #write bound for appropriate statusses
   if( status == "opt" || status == "best" || status == "bestdual" )
   {
      printf(" %16.9g", thebound) > NEWSOLUFILE;
   }

   printf("\n") > NEWSOLUFILE;
}

# does this status string represent a failed status?
function isStatusFail(thisstatus)
{
   # check if the first four characters are 'fail'
   if( substr(thisstatus, 1, 4) == "fail" )
      return 1;
   else
      return 0;
}

BEGIN {
   timegeomshift = 10.0;
   nodegeomshift = 100.0;
   sblpgeomshift = 0.0;
   pavshift = 0.0;
   onlyinsolufile = 0;          # should only instances be reported that are included in the .solu file?
   onlyintestfile = 0;          # should only instances be reported that are included in the .test file?  TEMPORARY HACK!
   onlypresolvereductions = 0;  # should only instances with presolve reductions be shown?
   if (useshortnames == "") {
       useshortnames = 1;       # should problem name be truncated to fit into column?
   }
   writesolufile = 0;           # should a solution file be created from the results? Use '1' for writing a new solution file, or '2' for writing an update
                                # respecting the previous solu file information and updating it by better solution values for previously unsolved instances
   printsoltimes = 0;           # should the times until first and best solution be shown
   checksol = 1;                # should the solution check of SCIP be parsed and counted as a fail if best solution is not feasible?
   analyseconf = 0;             # should conflict analysis be reported?
   NEWSOLUFILE = "new_solufile.solu";
   infty = +1e+20;
   headerprinted = 0;
   namelength = 35;             # maximal length of instance names (can be increased)
   usetimestamps = 0;

   nprobs = 0;
   sbab = 0;
   slp = 0;
   ssim = 0;
   ssblp = 0;
   stottime = 0.0;
   stimetofirst = 0.0;
   stimetobest = 0.0;
   nodegeom = 0.0;
   timegeom = 0.0;
   timetofirstgeom = 0.0;
   timetobestgeom = 0.0;
   sblpgeom = 0.0;
   conftimegeom = 0.0;
   confgeom = 0.0;
   basictimegeom = 0.0;
   overheadtimegeom = 0.0;
   shiftednodegeom = nodegeomshift;
   shiftedtimegeom = timegeomshift;
   shiftedsblpgeom = sblpgeomshift;
   shiftedconftimegeom = timegeomshift;
   shiftedconfgeom = timegeomshift;
   shiftedbasictimegeom = timegeomshift;
   shiftedoverheadtimegeom = timegeomshift;
   shiftedtimetofirstgeom = timegeomshift;
   shiftedtimetobestgeom = timegeomshift;
   timeouttime = 0.0;
   timeouts = 0;
   failtime = 0.0;
   fail = 0;
   pass = 0;
   settings = "default";
   lpsname = "?";
   lpsversion = "?";
   scipversion = "?";
   githash = "?";
   conftottime = 0.0;
   sumconfs = 0;
   overheadtottime = 0.0;

   #initialize paver input file
   if( PAVFILE != "" )
   {
      printf("* Trace Record Definition\n") > PAVFILE;
      printf("* InputFileName,ModelType,SolverName,Direction,ModelStatus,SolverStatus,ObjectiveValue,ObjectiveValueEstimate,SolverTime\n") > PAVFILE;
      printf("* NumberOfNodes,NumberOfIterations,NumberOfEquations,NumberOfVariables\n") > PAVFILE;
   }
}
/^IP\// || /^MINLP\// {  # HACK to parse .test files
   n = split ($1, a, "/");
   m = split(a[n], b, ".");
   prob = b[1];
   if( b[m] == "gz" || b[m] == "z" || b[m] == "GZ" || b[m] == "Z" )
      m--;
   for( i = 2; i < m; ++i )
      prob = prob "." b[i];
   intestfile[prob] = 1;
}
/=opt=/  {  # get optimum
   solstatus[$2] = "opt";
   sol[$2] = $3;
}
/=inf=/  {  # problem infeasible (no feasible solution exists)
   solstatus[$2] = "inf";
   sol[$2] = +infty;
}
/=best=/ {  # get best known primal bound
   solstatus[$2] = "best";
   sol[$2] = $3;
}
/=bestdual=/ {  # get best known dual bound
   bestdual[$2] = $3;
}
/=unkn=/ {  # no feasible solution known
   solstatus[$2] = "unkn";
}
#
# problem name
#
/^@01/ {
   filename = $2;
   grepresult = ""

   n = split ($2, a, "/");
   m = split(a[n], b, ".");
   prob = b[1];
   if( b[m] == "gz" || b[m] == "z" || b[m] == "GZ" || b[m] == "Z" )
      m--;
   for( i = 2; i < m; ++i )
      prob = prob "." b[i];

   if( useshortnames && length(prob) > namelength )
      shortprob = substr(prob, length(prob)-namelength-1, namelength);
   else
      shortprob = prob;

   # Escape _ for TeX
   n = split(prob, a, "_");
   pprob = a[1];
   for( i = 2; i <= n; i++ )
      pprob = pprob "\\_" a[i];
   vars = 0;
   binvars = 0;
   intvars = 0;
   implvars = 0;
   contvars = 0;
   cons = 0;
   lincons = 0;
   quadcons = 0;
   nonlincons = 0;
   origvars = 0;
   origcons = 0;
   objsense = 0;
   timeout = 0;
   feasible = 0;
   pb = +infty;
   objectivelimit = +infty;
   firstpb = +infty;
   db = -infty;
   dbset = 0;
   dbforobjsense = -infty;
   simpiters = 0;
   bbnodes = 0;
   primlps = 0;
   primiter = 0;
   duallps = 0;
   dualiter = 0;
   sblps = 0;
   sbiter = 0;
   tottime = 0.0;
   timetofirst = -1.0;
   timetobest = -1.0;
   inconflict = 0;
   inconstime = 0;
   confclauses = 0;
   confliterals = 0.0;
   conftime = 0.0;
   conf_prop  = 0;
   conf_infLP = 0;
   conf_bndEx = 0;
   conf_strbr = 0;
   conf_pseud = 0;
   overheadtime = 0.0;
   aborted = 1;
   readerror = 0;
   gapreached = 0;
   sollimitreached = 0;
   memlimitreached = 0;
   nodelimitreached = 0;
   objlimitreached = 0;
   starttime = 0.0;
   endtime = 0.0;
   timelimit = 0.0;
   inoriginalprob = 1;
   incons = 0;
   valgrinderror = 0;
   valgrindleaks = 0;
   bestsolfeas = 1;
   reoptimization = 0;
   niter = 0;
}

/@03/ {
   starttime = $2;
}
/@04/ {
   endtime = $2;
}

/^SCIP version/ {
   # get SCIP version
   scipversion = $3;

   # get name of LP solver
   if( $13 == "SoPlex" )
      lpsname = "spx";
   else if( $13 == "SoPlex2" )
      lpsname = "spx2";
   else if( $13 == "SoPlex1" )
      lpsname = "spx1";
   else if( $13 == "CPLEX" )
      lpsname = "cpx";
   else if( $13 == "NONE]" )
   {
      lpsname = "none";
      lpsversion = "";
   }
   else if( $13 == "Clp" )
      lpsname = "clp";
   else if( $13 == "MOSEK" )
      lpsname = "msk";
   else if( $13 == "Gurobi" )
      lpsname = "grb";
   else if( $13 == "QSopt" )
      lpsname = "qso";
   else if( $13 == "Xpress" )
      lpsname = "xprs";

    # get LP solver version
   if( NF >= 16 )
   {
      split($14, v, "]");
      lpsversion = v[1];
   }

   # get git hash
   if( $(NF-1) == "[GitHash:" ) {
      split($NF, v, "]");
      githash = v[1];
   }
}
/^SCIP> SCIP> / {
   $0 = substr($0, 13, length($0)-12);
}
/^SCIP> / {
   $0 = substr($0, 7, length($0)-6);
}
/^loaded parameter file/ {
   settings = $4;
   sub(/<.*settings\//, "", settings);
   sub(/\.set>/, "", settings);
}
/^reading user parameter file/ {
   settings = $5;
   sub(/<.*settings\//, "", settings);
   sub(/\.set>/, "", settings);
}
/^parameter <limits\/time> set to/ {
   timelimit = $5;
}
/^limits\/time =/ {
   timelimit = $3;
}
/set limits objective/ {
   objectivelimit = $4;
}

/^read problem/ { niter += 1; }

# check if reoptimization is enabled
/^Solving Nodes      :/ {
   if( $6 == "reactivated)" )
      reoptimization = 1;
}

# check for primalbound before statistic printing
/^Primal Bound       :/ {
   pb = $4
}

#
# get objective sense
#
/^  Objective sense  :/ {
   if ( $4 == "minimize" || $4 == "minimize\r")
      objsense = 1;
   if ( $4 == "maximize" || $4 == "maximize\r" )
      objsense = -1;

   # objsense is 0 otherwise
}

# SCIP API version >= 9
/^  Objective        :/ {
   if( objsense == 0 )
   {
      if ( $3 == "minimize," || $3 == "minimize,\r")
         objsense = 1;
      if ( $3 == "maximize," || $3 == "maximize,\r" )
         objsense = -1;

      # objsense is 0 otherwise
   }
}

#
# conflict analysis
#
/^Conflict Analysis  :/ {
   inconflict = 1;
}
/^  propagation      :/ {
   if( inconflict == 1 ) {
      conftime += $3; #confclauses += $5 + $7; confliterals += $5 * $6 + $7 * $8;
      conf_prop += $7;
   }
}
/^  infeasible LP    :/ {
   if( inconflict == 1 ) {
      conftime += $4; #confclauses += $6 + $8; confliterals += $6 * $7 + $8 * $9;
      conf_infLP += $8;
   }
}
/^  bound exceed. LP :/ {
   if( inconflict == 1 ) {
      conftime += $4; #confclauses += $6 + $8; confliterals += $6 * $7 + $8 * $9;
      conf_bndEx += $9;
   }
}
/^  strong branching :/ {
   if( inconflict == 1 ) {
      conftime += $4; #confclauses += $6 + $8; confliterals += $6 * $7 + $8 * $9;
      conf_strbr += $8;
   }
}
/^  pseudo solution  :/ {
   if( inconflict == 1 ) {
      conftime += $4; #confclauses += $6 + $8; confliterals += $6 * $7 + $8 * $9;
      conf_pseud += $8;
   }
}
/^  applied globally :/ {
   if( inconflict == 1 ) {
      confclauses += $8; confliterals += $8 * $9;
   }
}
/^  applied locally  :/ {
   if( inconflict == 1 ) {
      confclauses += $8; confliterals += $8 * $9;
   }
}
/^Separators         :/ {
   inconflict = 0;
}
/^Constraint Timings :/ {
   inconstime = 1;
}
#/^  logicor          :/ {
#   if( inconstime == 1 )
#      overheadtime += $3;
#}
/^Propagators        :/ {
   inconstime = 0;
}
/^  switching time   :/ {
   overheadtime += $4;
}
#
# problem size
#
/^Presolved Problem  :/ {
   inoriginalprob = 0;
}
/^  Variables        :/ {
   if( inoriginalprob )
      origvars = $3;
   else
   {
      vars = $3;
      intvars = $6;
      implvars = $8;
      contvars = $11;
      binvars = vars - intvars - implvars - contvars;
   }
}
/^  Constraints      :/ {
   if( inoriginalprob )
      origcons = $3;
   else
      cons = $3;
}

#
# count number of linear constraints
#
/^Constraints        :/ {
   incons = 1;
   lincons = 0;
   quadcons = 0;
   nonlincons = 0;
}
/^  knapsack         :/ || /^  setppc           :/ || /^  linear           :/ || /^  logicor          :/ || /^  varbound         :/ {
   if( incons == 1 )
   {
      n = split ($3, a, "+");
      lincons += a[1];
   }
}
/^  quadratic        :/ || /^  soc              :/ {
   if( incons == 1 )
   {
      n = split ($3, a, "+");
      quadcons += a[1];
   }
}
/^  bivariate        :/ || /^  nonlinear        :/ || /^  abspower         :/ {
   if( incons == 1 )
   {
      n = split ($3, a, "+");
      nonlincons += a[1];
   }
}
/^Constraint Timings :/ {
   incons = 0;
}

#
# solution
#
/^Original Problem   : no problem exists./ {
   readerror = 1;
}
/^SCIP Status        :/ {
   # replace / by \/ in filename
   fname = filename
   gsub(/\//, "\\/",fname);

   #grep between filename and next @01 for an error
   command = "sed -n '/"fname"/,/@01/p' "ERRFILE" | grep 'returned with error code'";
   command | getline grepresult;

   # set aborted flag correctly
   if( grepresult == "" )
      aborted = 0;

   close(command)
}

# grep for an objective limit induced infeasibility
/^SCIP Status        : problem is solved \[infeasible\] \(objective limit reached\)/ {
    objlimitreached = 1;
}
/solving was interrupted/ {
   timeout = 1;
}
/gap limit reached/ {
   gapreached = 1;
}
/solution limit reached/ {
   sollimitreached = 1;
}
/memory limit reached/ {
   memlimitreached = 1;
}
/node limit reached/ {
   nodelimitreached = 1;
}
/problem is solved/ {
   timeout = 0;
}
/best solution is not feasible in original problem/ {
   bestsolfeas = 0;
}

/Check SOL:/ {
   intcheck = $4;
   conscheck = $6;
   objcheck = $8;
   if( !intcheck || !conscheck || !objcheck )
      bestsolfeas = 0;
}

/^  First Solution   :/ {
   timetofirst = $11;
   firstpb = $4;
}
/^  Primal Bound     :/ {
   if( $4 == "infeasible" || $4 == "infeasible\r" )
   {
      pb = +infty;
      db = +infty;
      dbset = 1;
      feasible = 0;
   }
   else if( $4 == "-"  || $4 == "-\r")
   {
      pb = +infty;
      feasible = 0;
   }
   else if( $5 == "(user" && $6 == "objective" && $7 == "limit)" )
   {
      pb = $4;
      feasible = 0;
   }
   else
   {
      pb = $4;
      feasible = 1;
      timetobest = $11;
   }
}
/^  Dual Bound       :/ {
   if( $4 != "-" && $4 != "-\r" )
   {
      db = $4;
      dbset = 1;
   }
}
/^Dual Bound         :/ {
    dbforobjsense = $4; # in old SCIP log files, this value is used to detect the objective sense
}
#
# iterations
#
/^  primal LP        :/ {
   simpiters += $6;
}
/^  dual LP          :/ {
   simpiters += $6;
}
/^  barrier LP       :/ {
   simpiters += $6;
}
/^  nodes \(total\)    :/ {
   if( reoptimization == 1 )
      bbnodes += $4;
   else
      bbnodes = $4;
}
/^  primal LP        :/ {
   if( reoptimization == 1 )
   {
      primlps += $5;
      primiter += $6;
   }
   else
   {
      primlps = $5;
      primiter = $6;
   }
}
/^  dual LP          :/ {
   if( reoptimization == 1 )
   {
      duallps += $5;
      dualiter += $6;
   }
   else
   {
      duallps = $5;
      dualiter = $6;
   }
}
/^  strong branching :/ {
   if( reoptimization == 1 )
   {
      sblps += $5;
      sbiter += $6;
   }
   else
   {
      sblps = $5;
      sbiter = $6;
   }
}
#
# time
#
/^Solving Time       :/ {  # for older scip version ( < 2.0.1.3 )
   tottime = $4;
}
/^  solving          :/ {
   tottime = $3;
}
#
# valgrind check
#
/^==[0-9]*== ERROR SUMMARY:/       {
   valgrinderror = $4
}
/^==[0-9]*==    definitely lost:/  {
   valgrindleaks += $4
}
/^==[0-9]*==    indirectly lost:/  {
   valgrindleaks += $4
}
/^==[0-9]*==    possibly lost:/    {
   valgrindleaks += $4
}
#
# solver status overview (in order of priority):
# 1) solver broke before returning solution => abort
# 2) solver cut off the optimal solution (solu-file-value is not between primal and dual bound) => fail
#    (especially if problem is claimed to be solved but solution is not the optimal solution)
# 3) solver solved problem with the value in solu-file (if existing) => ok
# 4) solver solved problem which has no (optimal) value in solu-file => solved
#    (since we here don't detect the direction of optimization, it is possible
#     that a solver claims an optimal solution which contradicts a known feasible solution)
# 5) solver found solution better than known best solution (or no solution was known so far) => better
# 7) solver reached gaplimit or limit of number of solutions => gaplimit, sollimit
# 8) solver reached any other limit (like time or nodes) => timeout
# 9) otherwise => unknown
#
/^=ready=/ {

   #since the header depends on the parameter printsoltimes and settings it is no longer possible to print it in the BEGIN section
   if( !headerprinted )
   {
      ntexcolumns = 8 + (2 * printsoltimes);

      if( TEXFILE != "" )
      {
         #print header of tex file table
         printf("\\documentclass[leqno]{article}\n")                      >TEXFILE;
         printf("\\usepackage{a4wide}\n")                                 >TEXFILE;
         printf("\\usepackage{amsmath,amsfonts,amssymb,booktabs}\n")      >TEXFILE;
         printf("\\usepackage{supertabular}\n")                           >TEXFILE;
         printf("\\pagestyle{empty}\n\n")                                 >TEXFILE;
         printf("\\begin{document}\n\n")                                  >TEXFILE;
         printf("\\begin{center}\n")                                      >TEXFILE;
         printf("\\setlength{\\tabcolsep}{2pt}\n")                        >TEXFILE;
         printf("\\newcommand{\\g}{\\raisebox{0.25ex}{\\tiny $>$}}\n")    >TEXFILE;
         printf("\\tablehead{\n\\toprule\n")                              >TEXFILE;
         printf("Name                &  Conss &   Vars &     Dual Bound &   Primal Bound &  Gap\\%% &     Nodes &     Time ") >TEXFILE;
         if( printsoltimes )
            printf(" &     To First      &    To Last   ") > TEXFILE;
         printf("\\\\\n") > TEXFILE;
         printf("\\midrule\n}\n")                                         >TEXFILE;
         printf("\\tabletail{\n\\midrule\n")                              >TEXFILE;
         printf("\\multicolumn{%d}{r} \\; continue next page \\\\\n", ntexcolumns) >TEXFILE;
         printf("\\bottomrule\n}\n")                                      >TEXFILE;
         printf("\\tablelasttail{\\bottomrule}\n")                        >TEXFILE;
         printf("\\tablecaption{SCIP with %s settings}\n",settings)       >TEXFILE;
         printf("\\begin{supertabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lrrrrrrr") >TEXFILE;
         if( printsoltimes )
            printf("rr") > TEXFILE;
         printf("@{}}\n") > TEXFILE;
      }

      # print header of table when this regular expression is matched for the first time

      # prepare header
      hyphenstr = "";
      for (i = 0; i < namelength; ++i)
         hyphenstr = sprintf("%s-", hyphenstr);

      # first part: name of given length
      tablehead1 = hyphenstr;
      tablehead2 = sprintf("Name%*s", namelength-4, " ");
      tablehead3 = hyphenstr;

      # append rest of header
      if( reoptimization == 0 )
      {
         tablehead1 = tablehead1"+------+--- Original --+-- Presolved --+----------------+----------------+------+---------+--------+-------+";
         tablehead2 = tablehead2"| Type | Conss |  Vars | Conss |  Vars |   Dual Bound   |  Primal Bound  | Gap%% |  Iters  |  Nodes |  Time |";
         tablehead3 = tablehead3"+------+-------+-------+-------+-------+----------------+----------------+------+---------+--------+-------+";
      }
      else
      {
         tablehead1 = tablehead1"+------+--- Original --+-- Presolved --+-----------+----------------+----------------+---------+--------+-------+";
         tablehead2 = tablehead2"| Type | Conss |  Vars | Conss |  Vars | Reopt Its | last Dual Bnd  | last Primal Bnd|  Iters  |  Nodes |  Time |";
         tablehead3 = tablehead3"+------+-------+-------+-------+-------+-----------+----------------+----------------+---------+--------+-------+";
      }

      if( analyseconf == 1 )
      {
         tablehead1 = tablehead1"---------------Conflict Analysis---------------+";
         tablehead2 = tablehead2"  inf  |  bnd  |  sbr  |  pro  |  psol |  time |";
         tablehead3 = tablehead3"-------+-------+-------+-------+---------------+";
      }

      if( printsoltimes == 1 )
      {
         tablehead1 = tablehead1"----------+---------+";
         tablehead2 = tablehead2" To First | To Best |";
         tablehead3 = tablehead3"----------+---------+";
      }

      tablehead1 = tablehead1"--------\n";
      tablehead2 = tablehead2"        \n";
      tablehead3 = tablehead3"--------\n";

      printf(tablehead1);
      printf(tablehead2);
      printf(tablehead3);

      headerprinted = 1;
   }

   if( (!onlyinsolufile || prob in solstatus) &&
       (!onlyintestfile || intestfile[prob]) ) #&& conf_prop + conf_infLP + conf_bndEx + conf_pseud + conf_strbr > 0 )
   {
      # if sol file could not be read, fix status to be "unknown"
      if( !(prob in solstatus) )
         solstatus[prob] = "unkn";

      #avoid problems when comparing floats and integer (make everything float)
      temp = pb;
      pb = 1.0*temp;
      temp = db;
      db = 1.0*temp;
      temp = dbforobjsense;
      dbforobjsense = 1.0*temp;

      # if objsense could not be determined so far (output is maybe too old)
      if( objsense == 0 )
      {
         reltol = 1e-5 * max(abs(pb),1.0);
         abstol = 1e-4;

         # firstpb and dbforobjsense are used to detect the direction of optimization (min or max)
         if( timetofirst < 0.0 )
            temp = pb;
         else
            temp = firstpb;
         firstpb = 1.0*temp;

         if ( firstpb - dbforobjsense > max(abstol,reltol) )
            objsense = 1;   # minimize
         else
            objsense = -1;  # maximize
      }

      # treat primal and dual bound differently if objective limit was reached
      if( objlimitreached && objectivelimit < +infty )
      {
          pb = objectivelimit;
          db = objectivelimit;
      }

      # modify primal bound for maximization problems without primal solution
      if( objsense == -1 && pb >= +infty )
         pb = -1.0 * pb;

      # modify dual bound for infeasible maximization problems
      if( objsense == -1 && (db >= +infty || dbset == 0) )
         db = -1.0 * db;

      nprobs++;

      markersym = "\\g";
      if( abs(pb - db) < 1e-06 && pb < infty )
      {
         gap = 0.0;
         markersym = "  ";
      }
      else if( abs(db) < 1e-06 )
         gap = -1.0;
      else if( abs(pb) < 1e-06 )
         gap = -1.0;
      else if( pb*db < 0.0 )
         gap = -1.0;
      else if( abs(db) >= +infty )
         gap = -1.0;
      else if( abs(pb) >= +infty )
         gap = -1.0;
      else
         gap = 100.0*abs((pb-db)/min(abs(db),abs(pb)));
      if( gap < 0.0 )
         gapstr = "  --  ";
      else if( gap < 1e+04 )
         gapstr = sprintf("%6.1f", gap);
      else
         gapstr = " Large";

      if( vars == 0 )
         probtype = "   --";
      else if( lincons < cons )
      {
         if( cons == lincons+quadcons )
         {
            if( binvars == 0 && intvars == 0 )
               probtype = "  QCP";
            else
               probtype = "MIQCP";
         }
         else if( cons == lincons+quadcons+nonlincons )
         {
            if( binvars == 0 && intvars == 0 )
               probtype = "  NLP";
            else
               probtype = "MINLP";
         }
         else
            probtype = "  CIP";
      }
      else if( binvars == 0 && intvars == 0 )
         probtype = "   LP";
      else if( contvars == 0 ) {
         if( intvars == 0 && implvars == 0 )
            probtype = "   BP";
         else
            probtype = "   IP";
      }
      else
      {
         if( intvars == 0 )
            probtype = "  MBP";
         else
            probtype = "  MIP";
      }

      if( aborted && endtime - starttime > timelimit && timelimit > 0.0 )
      {
         timeout = 1;
         tottime = endtime - starttime;
      }
      else if( gapreached || sollimitreached || memlimitreached || nodelimitreached )
      {
         timeout = 0;
         if( memlimitreached )
            tottime = max(endtime - starttime, timelimit);
      }

      if( aborted && tottime == 0.0 )
         tottime = timelimit;
      if( timelimit > 0.0 )
         tottime = min(tottime, timelimit);

      if( usetimestamps != 0 )
         tottime = endtime - starttime;

      if( aborted || timetobest < 0.0 )
      {
         timetofirst = tottime;
         timetobest = tottime;
      }


      lps = primlps + duallps;
      simplex = primiter + dualiter;
      stottime += tottime;
      stimetofirst += timetofirst;
      stimetobest += timetobest;
      sbab += bbnodes;
      slp += lps;
      ssim += simplex;
      ssblp += sblps;
      conftottime += conftime;
      sumconfs += (conf_infLP + conf_bndEx + conf_prop + conf_pseud + conf_strbr);
      overheadtottime += overheadtime;
      basictime = tottime - conftime - overheadtime;

      nodegeom = nodegeom^((nprobs-1)/nprobs) * max(bbnodes, 1.0)^(1.0/nprobs);
      sblpgeom = sblpgeom^((nprobs-1)/nprobs) * max(sblps, 1.0)^(1.0/nprobs);
      timegeom = timegeom^((nprobs-1)/nprobs) * max(tottime, 1.0)^(1.0/nprobs);
      overheadtimegeom = overheadtimegeom^((nprobs-1)/nprobs) * max(overheadtime, 1.0)^(1.0/nprobs);
      basictimegeom = basictimegeom^((nprobs-1)/nprobs) * max(basictime, 1.0)^(1.0/nprobs);

      shiftednodegeom = shiftednodegeom^((nprobs-1)/nprobs) * max(bbnodes+nodegeomshift, 1.0)^(1.0/nprobs);
      shiftedsblpgeom = shiftedsblpgeom^((nprobs-1)/nprobs) * max(sblps+sblpgeomshift, 1.0)^(1.0/nprobs);
      shiftedtimegeom = shiftedtimegeom^((nprobs-1)/nprobs) * max(tottime+timegeomshift, 1.0)^(1.0/nprobs);
      shiftedoverheadtimegeom = shiftedoverheadtimegeom^((nprobs-1)/nprobs) * max(overheadtime+timegeomshift, 1.0)^(1.0/nprobs);
      shiftedbasictimegeom = shiftedbasictimegeom^((nprobs-1)/nprobs) * max(basictime+timegeomshift, 1.0)^(1.0/nprobs);

      shiftedconftimegeom = shiftedconftimegeom^((nprobs-1)/nprobs) * max(conftime+timegeomshift, 1.0)^(1.0/nprobs);
      shiftedconfgeom = shiftedconfgeom^((nprobs-1)/nprobs) * max((conf_infLP + conf_bndEx + conf_prop + conf_pseud + conf_strbr)+timegeomshift, 1.0)^(1.0/nprobs);
      conftimegeom = conftimegeom^((nprobs-1)/nprobs) * max(conftime, 1.0)^(1.0/nprobs);
      confgeom = confgeom^((nprobs-1)/nprobs) * max((conf_infLP + conf_bndEx + conf_prop + conf_pseud + conf_strbr), 1.0)^(1.0/nprobs);

      timetobestgeom = timetobestgeom^((nprobs-1)/nprobs) * max(timetobest,1.0)^(1.0/nprobs);
      timetofirstgeom = timetofirstgeom^((nprobs-1)/nprobs) * max(timetofirst,1.0)^(1.0/nprobs);
      shiftedtimetofirstgeom = shiftedtimetofirstgeom^((nprobs-1)/nprobs) * max(timetofirst + timegeomshift, 1.0)^(1.0/nprobs);
      shiftedtimetobestgeom = shiftedtimetobestgeom^((nprobs-1)/nprobs) * max(timetobest + timegeomshift, 1.0)^(1.0/nprobs);

      status = "";
      reltol = 1e-5 * max(abs(pb),1.0);
      abstol = 1e-4;

      if( readerror )
      {
         setStatusToFail("fail (readerror)");
      }
      else if( aborted )
      {
         setStatusToFail("fail (abort)");
      }
      else if( checksol && !bestsolfeas )
      {
         setStatusToFail("fail (solution infeasible)");
      }
      else if( solstatus[prob] == "opt" )
      {
         # in case a solution was found we compare primal and dual bound
         if( feasible && ( isPrimalBoundBetter() || isDualBoundBetter() ) )
         {
            setStatusToFail("fail (objective value)");
         }
         else if( !feasible && objlimitreached )
         {
            # if the objective limit was at least as tight as the optimal solution value, we accept the infeasibility
            if( (objsense == 1 && sol[prob]-objectivelimit >= -reltol) || (objsense == -1 && (objectivelimit - sol[prob] >= -reltol)) )
            {
               status = "ok";
               pass++;
            }
            else
            {
               setStatusToFail("fail (objective value)")
            }
         }
         else if( isLimitReached() )
         {
            setStatusToLimit();
         }
         else if( (db == -infty || isPrimalDualBoundEqual()) && !isPrimalBoundBetter() && !isDualBoundBetter() )
         {
            status = "ok";
            pass++;
         }
         else
         {
            setStatusToFail("fail");
         }
      }
      else if( solstatus[prob] == "best" || prob in bestdual )
      {
         # we failed if the dual bound was higher/lower than the best known primal bound
         if( isDualBoundBetter() )
         {
            setStatusToFail("fail (dual bound)");
         }
         else if( isPrimalBoundBetterThanBestDual() )
         {
            setStatusToFail("fail (primal bound)");
         }
         else if( isLimitReached() )
         {
            setStatusToLimit();

            if( isPrimalBoundBetter() )
               status = "better";
         }
         else if( isPrimalDualBoundEqual() )
         {
               status = "solved not verified";
               pass++;
         }
         else
         {
             setStatusToFail("fail (stopped unsolved before limit)");
         }
      }
      else if( solstatus[prob] == "unkn" )
      {
         if( isLimitReached() )
         {
            setStatusToLimit();

            if( abs(pb) < infty )
               status = "better";
         }
         else if( isPrimalDualBoundEqual() )
         {
            status = "solved not verified";
            pass++;
         }
         else
         {
            status = "unknown";
         }
      }
      else if( solstatus[prob] == "inf" )
      {
         if( !feasible )
         {
            if( timeout || memlimitreached || nodelimitreached )
            {
               setStatusToLimit();
            }
            else
            {
               status = "ok";
               pass++;
            }
         }
         else
         {
            setStatusToFail("fail (solution on infeasible instance)");
         }
      }
      else
      {
         if( isLimitReached() )
         {
            setStatusToLimit();
         }
         else if( isPrimalDualBoundEqual() )
         {
            status = "solved not verified";
            pass++;
         }
         else
         {
            status = "unknown";
         }
      }

      if( valgrinderror > 0 || valgrindleaks > 0 )
      {
         setStatusToFail("fail (valgrind)")
      }

      # write some solu file information
      if( writesolufile )
      {
         newsolpb = pb;

         if( pb == +infty && db == +infty )
             newsolstatus = "inf";
         else if( isPrimalDualBoundEqual() )
            newsolstatus = "opt";
         else if( abs(pb) < +infty )
            newsolstatus = "best";
         else
            newsolstatus = "unkn";

         # in update mode, use values from this run only if the primal bound is strictly better than the current best known
         # in case of a tie, we trust the solu-file more than the current results
         if( writesolufile == 2 )
         {
            # an opt or inf status from the solu-file always wins. In case of a best status in the solu-file, we use the value of this run
            # if the primal bound is strictly better or the instance was solved to proven optimality respecting the value from the solu file
            if( solstatus[prob] == "opt" ||
                solstatus[prob] == "inf" ||
                (solstatus[prob] == "best" && !isPrimalBoundBetter() && (!isPrimalDualBoundEqual() || isDualBoundBetter())) )
            {
               newsolstatus = solstatus[prob];
               newsolpb = sol[prob];
            }
         }

         writeToSoluFile(newsolstatus, newsolpb);

         # write dual bound, if we wrote a best known primal bound
         if( newsolstatus == "best" && abs(db) < +infty )
         {
            # in update mode, use value from solu-file if present and better
            if( writesolufile == 2 && (prob in bestdual) )
               newsoldb = (objsense == 1) ? max(bestdual[prob], db) : min(bestdual[prob], db);
            else
               newsoldb = db;
            writeToSoluFile("bestdual", newsoldb);
         }
      }

      #write output to both the tex file and the console depending on whether printsoltimes is activated or not
      if( !onlypresolvereductions || origcons > cons || origvars > vars )
      {
         if (TEXFILE != "")
         {
            if( reoptimization == 0 )
            {
               printf("%-*s & %6d & %6d & %16.9g & %16.9g & %6s &%s%8d &%s%7.1f",
                      namelength, pprob, cons, vars, db, pb, gapstr, markersym, bbnodes, markersym, tottime)  >TEXFILE;
            }
            else
            {
               printf("%-*s & %6d & %6d & %10d & %16.9g & %16.9g &%s%8d &%s%7.1f",
                      namelength, pprob, cons, vars, niter, db, pb, markersym, bbnodes, markersym, tottime)  >TEXFILE;
            }

	    if( analyseconf == 1 )
            {
               printf(" & %7d & %7d & %7d & %7d & %7d & %7.1f ", conf_infLP, conf_bndEx, conf_strbr, conf_prop, conf_pseud, conftime)  >TEXFILE;
            }

	    if( printsoltimes )
               printf(" & %7.1f & %7.1f", timetofirst, timetobest) > TEXFILE;
            printf("\\\\\n") > TEXFILE;
         }

         # note: probtype has length 5, but field width is 6
         if( reoptimization == 0 )
         {
            printf("%-*s  %-5s %7d %7d %7d %7d %16.9g %16.9g %6s %9d %8d %7.1f ",
                   namelength, shortprob, probtype, origcons, origvars, cons, vars, db, pb, gapstr, simpiters, bbnodes, tottime);
         }
         else
         {
            printf("%-*s  %-5s %7d %7d %7d %7d %11d %16.9g %16.9g %9d %8d %7.1f ",
                   namelength, shortprob, probtype, origcons, origvars, cons, vars, niter, db, pb, simpiters, bbnodes, tottime);
	 }

	 if( analyseconf == 1 )
         {
            printf("%7d %7d %7d %7d %7d %7.1f ", conf_infLP, conf_bndEx, conf_strbr, conf_prop, conf_pseud, conftime);
         }

	 if( printsoltimes )
            printf(" %9.1f %9.1f ", timetofirst, timetobest);

         printf("%s\n", status);
      }

      if( PAVFILE != "" )
      {
         #PAVER output: see http://www.gamsworld.org/performance/paver/pprocess_submit.htm
         if( status == "fail (abort)" )
         {
            modelstat = 13;
            solverstat = 13;
         }
         else if( isStatusFail(status) || status == "unknown" )
         {
            modelstat = 7;
            solverstat = 1;
         }
         else if( status == "timeout" )
         {
            modelstat = abs(pb) < infty ? 8 : 14;
            solverstat = 3;
         }
         else if( status == "nodelimit" || status == "memlimit" || status == "sollimit" )
         {
            modelstat = abs(pb) < infty ? 8 : 14;
            solverstat = 2;  # GAMS does not have a status for these limits, so let's report iteration limit
         }
         else if( status == "gaplimit" || status == "better" )
         {
            modelstat = 8;
            solverstat = 1;
         }
         else if( status == "ok" || status == "solved not verified" )
         {
            modelstat = 1;
            solverstat = 1;
         }
         else
         {
            modelstat = 13;
            solverstat = 13;
         }
         pavprob = prob;
         if( length(pavprob) > 25 )
            pavprob = substr(pavprob, length(pavprob)-24,25);
         if( vars == 0 )
            gamsprobtype = "LP";
         else if( lincons < cons && binvars == 0 && intvars == 0 )
            gamsprobtype = "NLP";
         else if( lincons < cons )
            gamsprobtype = "MINLP";
         else if( binvars == 0 && intvars == 0 )
            gamsprobtype = "LP";
         else
            gamsprobtype = "MIP";
         #InputFileName,ModelType,SolverName,Direction,ModelStatus,SolverStatus,ObjectiveValue,ObjectiveValueEstimate,SolverTime
         #NumberOfNodes,NumberOfIterations,NumberOfEquations,NumberOfVariables
         printf("%s,%s,SCIP_%s,%d,%d,%d,%.8g,%.8g,%g,", pavprob, gamsprobtype, settings, objsense == 1 ? 0 : 1, modelstat, solverstat, pb, db, tottime+pavshift) > PAVFILE;
         printf("%d,%d,%d,%d\n", bbnodes, simpiters, cons, vars) > PAVFILE;
      }
   }
}
END {
   shiftednodegeom -= nodegeomshift;
   shiftedsblpgeom -= sblpgeomshift;
   shiftedtimegeom -= timegeomshift;
   shiftedconftimegeom -= timegeomshift;
   shiftedoverheadtimegeom -= timegeomshift;
   shiftedbasictimegeom -= timegeomshift;
   shiftedtimetofirstgeom -= timegeomshift;
   shiftedtimetobestgeom -= timegeomshift;

   if( TEXFILE != "" )
   {
      printf("\\midrule\n")                                                 >TEXFILE;
      printf("%-14s (%2d) &        &        &                &                &        & %9d & %8.1f",
             "Total", nprobs, sbab, stottime) >TEXFILE;
      if( analyseconf )
         printf(" & %8.1f & %8.1f", allconfs, conftime) > TEXFILE;
      if( printsoltimes )
         printf(" & %8.1f & %8.1f", stimetofirst, stimetobest) > TEXFILE;
      printf("\\\\\n") > TEXFILE;
      printf("%-14s      &        &        &                &                &        & %9d & %8.1f",
             "Geom. Mean", nodegeom, timegeom) >TEXFILE;
      if( analyseconf )
         printf(" & %8.1f & %8.1f", confgeom, conftimegeom) > TEXFILE;
      if( printsoltimes )
         printf(" & %8.1f & %8.1f", timetofirstgeom, timetobestgeom) > TEXFILE;
      printf("\\\\\n") > TEXFILE;
      printf("%-14s      &        &        &                &                &        & %9d & %8.1f ",
             "Shifted Geom.", shiftednodegeom, shiftedtimegeom) >TEXFILE;
      if( analyseconf )
         printf(" & %8.1f & %8.1f", shiftedconfgeom, shiftedconftimegeom) > TEXFILE;
      if( printsoltimes )
         printf(" & %8.1f & %8.1f", shiftedtimetofirstgeom, shiftedtimetobestgeom) > TEXFILE;
      printf("\\\\\n") > TEXFILE;
   }
   printf(tablehead3);
   printf("\n");

   tablefooter1 = "------------------------------[Nodes]---------------[Time]------";
   tablefooter2 = "  Cnt  Pass  Time  Fail  total(k)     geom.     total     geom.";
   tablefooter3 = "----------------------------------------------------------------";

   if( analyseconf == 1 )
   {
      tablebottom1 = tablebottom1"--------[NConf]-----------[ConfTime]-----";
      tablebottom2 = tablebottom2"     total     geom.     total     geom.";
      tablebottom3 = tablebottom3"-----------------------------------------";
   }

   if( printsoltimes ) {
      tablebottom1 = tablebottom1"--------[ToFirst]-----------[ToLast]-----";
      tablebottom2 = tablebottom2"     total     geom.     total     geom.";
      tablebottom3 = tablebottom3"-----------------------------------------";
   }

   tablefooter1 = tablefooter1"\n";
   tablefooter2 = tablefooter2"\n";
   tablefooter3 = tablefooter3"\n";

   printf(tablefooter1);
   printf(tablefooter2);
   printf(tablefooter3);

   printf("%5d %5d %5d %5d %9d %9.1f %9.1f %9.1f ",
          nprobs, pass, timeouts, fail, sbab / 1000, nodegeom, stottime, timegeom);

   if( analyseconf == 1 )
     printf("%9d %9.1f %9.1f %9.1f", sumconfs, confgeom, conftottime, conftimegeom);

   if( printsoltimes )
      printf("%9.1f %9.1f %9.1f %9.1f", stimetofirst, timetofirstgeom, stimetobest, timetobestgeoconftimem);

   printf("\n");
   printf(" shifted geom. [%5d/%5.1f]      %9.1f           %9.1f ",
          nodegeomshift, timegeomshift, shiftednodegeom, shiftedtimegeom);
   if( analyseconf )
      printf("          %9.1f           %9.1f ", shiftedconfgeom, shiftedconftimegeom);
   if( printsoltimes )
      printf("          %9.1f           %9.1f ", shiftedtimetofirstgeom, shiftedtimetobestgeom);
   printf("\n");
   printf(tablefooter3);

   if( TEXFILE != "" )
   {
      printf("\\noalign{\\vspace{6pt}}\n")                                  >TEXFILE;
      printf("\\end{supertabular*}\n")                                      >TEXFILE;
      printf("\\end{center}\n")                                             >TEXFILE;
      printf("\\end{document}\n")                                           >TEXFILE;
   }

   printf("@02 timelimit: %g\n", timelimit);
   printf("@01 SCIP(%s)%s(%s):%s", scipversion, lpsname, lpsversion, settings);
   if( githash != "?" )
      printf(" [GitHash: %s]\n", githash);
   else
      printf("\n");

   if( TEXFILE != "" )
      close(TEXFILE);
   if( PAVFILE != "" )
      close(PAVFILE);
}

#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#@file    check_gams.awk
#@brief   GAMS Tracefile Check Report Generator
#@author  Robert Waniek
#@author  Stefan Vigerske
#@author  check.awk (copied large portions from there)
#
function max(x,y)
{
    return (x) > (y) ? (x) : (y);
}
function min(x,y)
{
   return (x) < (y) ? (x) : (y);
}
function abs(a)
{
  if (a>0) return a;
  else return -a;
}
function isEQ(a, b)
{
  scale = abs(a);
  if ( scale >= 1e+10 )
    scale = abs(b);
  if ( scale < 0.0001 )
    scale = 1;
  return abs(a-b) < 0.0002 * scale;
}

BEGIN  { 
   timegeomshift = 10.0;
   nodegeomshift = 100.0;
   onlyinsolufile = 0;       # should only instances be reported that are included in the .solu file?
   useshortnames = 1;        # should problem name be truncated to fit into column?
   namelength = 18;          # maximal length of instance names (can be increased)
   writesolufile = 0;        # should a solution file be created from the results
   NEWSOLUFILE = "new_solufile.solu";
   infty = 1e+20;

   if (solver == "") solver="SCIP";
   if (penaltytime == "") penaltytime=3600.0;

   nprobs = 0;
   sbab = 0;
   ssim = 0;
   stottime = 0.0;
   nodegeom = 1.0;
   timegeom = 1.0;
   shiftednodegeom = 1.0;
   shiftedtimegeom = 1.0;
   timeouttime = 0.0;
   timeouts = 0;
   failtime = 0.0;
   fail = 0;
   pass = 0;
   timelimit = 0;
   settings = "default";
}
/=opt=/  { solstatus[$2] = "opt"; sol[$2] = $3; }   # get optimum
/=inf=/  { solstatus[$2] = "inf"; }                 # problem infeasible (no feasible solution exists)
/=best=/ { solstatus[$2] = "best"; sol[$2] = $3; }  # get best known solution value
/=unkn=/ { solstatus[$2] = "unkn"; }                # no feasible solution known
/^\*/    { FS="," }  # a start at the beginnnig of a line we take as start of tracefile, so change FS to ','
/^\* SOLVER/     { solver=$2; }
/^\* TIMELIMIT/  { timelimit=$2; }
/^\* SETTINGS/   { settings=$2; }
/^\*/  { next; } # skip other comments and invalid problems
/^ *$/ { next; } # skip empty lines

#These need to coincide with those in in check_gams.sh
#TODO make this more flexible (see readtrace.awk from ptools)
#01 InputFileName
#02 ModelType
#03 SolverName
#04 OptionFile
#05 Direction
#06 NumberOfEquations
#07 NumberOfVariables
#08 NumberOfDiscreteVariables
#09 NumberOfNonZeros
#10 NumberOfNonlinearNonZeros
#11 ModelStatus
#12 SolverStatus
#13 ObjectiveValue
#14 ObjectiveValueEstimate
#15 SolverTime
#16 ETSolver
#17 NumberOfIterations
#18 NumberOfNodes

/.*/ {
  if( $3 == solver || $3 == "EXAMINER2" )
  {
    model[nprobs] = $1;
    type[nprobs] = $2;
    maxobj[nprobs] = ( $5 == 1 ? 1 : 0 );
    cons[nprobs] = $6;
    vars[nprobs] = $7;
    modstat[nprobs] = $11;
    solstat[nprobs] = $12;
    dualbnd[nprobs] = $14;
    primalbnd[nprobs] = $13;
    time[nprobs] = $15; # use max($15,$16) for max of solver reported time and wallclock time
    iters[nprobs] = $17;
    nodes[nprobs] = $18;
    nprobs++;
  }
}

END {
  # prepare header
  hyphenstr = "";
  for (i = 0; i < namelength; ++i)
     hyphenstr = sprintf("%s-", hyphenstr);

  # first part: name of given length
  tablehead1 = hyphenstr;
  tablehead2 = sprintf("Name%*s", namelength-4, " ");
  tablehead3 = hyphenstr;

  # append rest of header
  tablehead1 = tablehead1"+------+--- Original --+-- Presolved --+----------------+----------------+------+---------+--------+-------+-------\n";
  tablehead2 = tablehead2"| Type | Conss |  Vars | Conss |  Vars |   Dual Bound   |  Primal Bound  | Gap%% |  Iters  |  Nodes |  Time |       \n";
  tablehead3 = tablehead3"+------+-------+-------+-------+-------+----------------+----------------+------+---------+--------+-------+-------\n";

  # print header
  printf(tablehead1);
  printf(tablehead2);
  printf(tablehead3);

  # initialize paver input file
  printf("* Trace Record Definition\n") > PAVFILE;
  printf("* InputFileName,ModelType,SolverName,Direction,ModelStatus,SolverStatus,ObjectiveValue,ObjectiveValueEstimate,SolverTime\n") > PAVFILE;
  printf("* NumberOfNodes,NumberOfIterations,NumberOfEquations,NumberOfVariables\n") > PAVFILE;

  for (m = 0; m < nprobs; m++)
  {
     prob = model[m];
     if( useshortnames && length(prob) > namelength )
       shortprob = substr(prob, length(prob)-namelength-1, namelength);
     else
       shortprob = prob;

     if (primalbnd[m] == "NA")
       primalbnd[m] = ( maxobj[m] ? -infty : +infty );
       
     # do not trust primal bound in trace file if model status indicates that no feasible solution was found
     if( modstat[m] != 1 && modstat[m] != 2 && modstat[m] != 3 && modstat[m] != 7 && modstat[m] != 8 )
       primalbnd[m] = ( maxobj[m] ? -infty : +infty );

     # if dual bound is not given, but solver claimed model status "optimal", then we set dual bound to primal bound
     if (dualbnd[m] == "NA")
       dualbnd[m] = ( modstat[m] == 1 ? primalbnd[m] : ( maxobj[m] ? +infty : -infty ) );

     db = dualbnd[m];
     pb = primalbnd[m];
     
     # we consider every solver status between 4 and 7 and above 8 as unusal interrupt, i.e., abort
     # (8 is user interrupt, which is likely to be due to schulz stopping a solver on the hard timelimit)
     aborted = 0;
     if( solstat[m] >= 4 && solstat[m] != 8 )
       aborted = 1;
       
     if( nodes[m] == "NA" )
       nodes[m] = 0;

     # we consider iteration limit reached as node limit reached, because there is no extra status for node limits
     nodelimreached = (solstat[m] == 2);

     # TODO consider gaplimit
     gapreached = 0;
     
     timeout = 0;

     if( !onlyinsolufile || solstatus[prob] != "" )  {

       #avoid problems when comparing floats and integer (make everything float)
       temp = pb;
       pb = 1.0*temp;
       temp = db;
       db = 1.0*temp;
      
       optimal = 0;
       markersym = "\\g";
       if( abs(pb - db) < 1e-06 && pb < infty ) {
         gap = 0.0;
         optimal = 1;
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
      
       probtype = type[m];

       tottime = time[m];
       if( time[m] >= timelimit && timelimit > 0.0 )
         timeout = 1;
       else if( gapreached || nodelimreached )
         timeout = 0;
       if( tottime == 0.0 )
         tottime = timelimit;
       if( timelimit > 0.0 )
         tottime = min(tottime, timelimit);

       simplex = iters[m];
       stottime += tottime;
       sbab += nodes[m];
       ssim += simplex;

       nodegeom = nodegeom * max(nodes[m], 1.0)^(1.0/nprobs);
       timegeom = timegeom * max(tottime, 1.0)^(1.0/nprobs);

       shiftednodegeom = shiftednodegeom * max(nodes[m]+nodegeomshift, 1.0)^(1.0/nprobs);
       shiftedtimegeom = shiftedtimegeom * max(tottime+timegeomshift, 1.0)^(1.0/nprobs);

       status = "";
       if( aborted ) {
         status = "abort";
         failtime += tottime;
         fail++;
       }
       else if( solstatus[prob] == "opt" ) {
         reltol = 1e-4 * max(abs(pb),1.0);
         abstol = 1e-4;

         if( ( !maxobj[m] && (db-sol[prob] > reltol || sol[prob]-pb > reltol) ) || ( maxobj[m] && (sol[prob]-db > reltol || pb-sol[prob] > reltol) ) ) {
            status = "fail";
            failtime += tottime;
            fail++;
         }
         else {
           if( timeout || gapreached || nodelimreached ) {
             if( timeout )
                status = "timeout";
             else if( gapreached )
                status = "gaplimit";
             else if( nodelimreached )
                status = "nodelimit";
             timeouttime += tottime;
             timeouts++;
           }
           else {
             if( (abs(pb - db) <= max(abstol, reltol)) && abs(pb - sol[prob]) <= reltol ) {
               status = "ok";
               pass++;
             }
             else {
               status = "fail";
               failtime += tottime;
               fail++;
             }
           }
         }
       }
       else if( solstatus[prob] == "best" ) {
         reltol = 1e-4 * max(abs(pb),1.0);
         abstol = 1e-4;

         if( ( !maxobj[m] && db-sol[prob] > reltol) || ( maxobj[m] && sol[prob]-db > reltol) ) {
           status = "fail";
           failtime += tottime;
           fail++;
         }
         else {
           if( timeout || gapreached || nodelimreached ) {
             if( (!maxobj[m] && sol[prob]-pb > reltol) || (maxobj[m] && pb-sol[prob] > reltol) ) {
               status = "better";
               timeouttime += tottime;
               timeouts++;
             }
             else {
               if( timeout )
                 status = "timeout";
               else if( gapreached )
                 status = "gaplimit";
               else if( nodelimreached )
                 status = "nodelimit";
               timeouttime += tottime;
               timeouts++;
             }
           }
           else {
             if( abs(pb - db) <= max(abstol, reltol) ) {
               status = "solved not verified";
               pass++;
             }
             else {
               status = "fail";
               failtime += tottime;
               fail++;
             }
           }
         }
       }
       else if( solstatus[prob] == "unkn" ) {
	  reltol = 1e-4 * max(abs(pb),1.0);
	  abstol = 1e-4;

	  if( timeout || gapreached || nodelimreached ) {
	     if( abs(pb) < infty ) {
		status = "better";
		timeouttime += tottime;
		timeouts++;
	     }
	     else {
		if( timeout )
		   status = "timeout";
		else if( gapreached )
		   status = "gaplimit";
		else if( nodelimreached )
			status = "nodelimit";
		timeouttime += tottime;
		timeouts++;
	     }
	  }
	  else if( abs(pb - db) <= max(abstol, reltol) ) {
	     status = "solved not verified";
	     pass++;
	  }
	  else {
	     status = "unknown";
	  }
       }
       else if( solstatus[prob] == "inf" ) {
	  if( !feasible ) {
	     if( timeout ) {
		status = "timeout";
		timeouttime += tottime;
		timeouts++;
	     }
	     else {
		status = "ok";
		pass++;
	     }
	  }
	  else {
	     status = "fail";
	     failtime += tottime;
	     fail++;
	  }
       }
       else {
	  reltol = 1e-4 * max(abs(pb),1.0);
	  abstol = 1e-4;

	  if( timeout || gapreached || nodelimreached ) {
	     if( timeout )
		status = "timeout";
	     else if( gapreached )
		status = "gaplimit";
		  else if( nodelimreached )
		status = "nodelimit";
	     timeouttime += tottime;
	     timeouts++;
	  }
	  else if( abs(pb - db) < max(abstol,reltol) ) {
	     status = "solved not verified";
	     pass++;
	  }
	  else {
	     status = "unknown";
	  }
       }

       if( writesolufile ) {
         if( pb == +infty && db == +infty )
           printf("=inf= %-18s\n",prob)>NEWSOLUFILE;
         else if( pb == db )
           printf("=opt= %-18s %16.9g\n",prob,pb)>NEWSOLUFILE;
         else if( pb < +infty )
           printf("=best= %-18s %16.9g\n",prob,pb)>NEWSOLUFILE;
         else
           printf("=unkn= %-18s\n",prob)>NEWSOLUFILE;
       }

       #write output to both the tex file and the console depending on whether printsoltimes is activated or not
       printf("%-*s & %6d & %6d & %16.9g & %16.9g & %6s &%s%8d &%s%7.1f",
              namelength, pprob, cons[m], vars[m], db, pb, gapstr, markersym, nodes[m], markersym, time[m]) > TEXFILE;
       printf("\\\\\n") > TEXFILE;

       # note: probtype has length 5, but field width is 6
       printf("%-*s  %-5s %7d %7d      ??      ?? %16.9g %16.9g %6s %9d %8d %7.1f %s (%2d - %2d)\n",
              namelength, shortprob, probtype, cons[m], vars[m], db, pb, gapstr, iters[m], nodes[m], time[m], status, modstat[m], solstat[m]);

       #PAVER output: see http://www.gamsworld.org/performance/paver/pprocess_submit.htm
       if( status == "abort" ) {
         modelstat = 13;
         solverstat = 13;
       } else if( status == "fail" || status == "unknown" ) {
         modelstat = 7;
         solverstat = 1;
       } else if( status == "timeout" ) {
         modelstat = abs(pb) < infty ? 8 : 14;
         solverstat = 3;
       } else if( status == "nodelimit" ) {
         modelstat = abs(pb) < infty ? 8 : 14;
         solverstat = 2;  # GAMS does not have a status for these limits, so let's report iteration limit
       } else if( status == "gaplimit" || status == "better" ) {
         modelstat = 8;
         solverstat = 1;
       } else if( status == "ok" || status == "solved not verified" ) {
         modelstat = 1;
         solverstat = 1;
       } else {
         modelstat = 13;
         solverstat = 13;
       }
       pavprob = prob;
       if( length(pavprob) > 25 )
         pavprob = substr(pavprob, length(pavprob)-24,25);
       #InputFileName,ModelType,SolverName,Direction,ModelStatus,SolverStatus,ObjectiveValue,ObjectiveValueEstimate,SolverTime
       #NumberOfNodes,NumberOfIterations,NumberOfEquations,NumberOfVariables
       printf("%s,MINLP,%s_%s,%d,%d,%d,%g,%g,%g,", pavprob, solver, settings, maxobj[m] ? 1 : 0, modelstat, solverstat, pb, db, time[m]) > PAVFILE;
       printf("%d,%d,%d,%d\n", nodes[m], iters[m], cons[m], vars[m]) > PAVFILE;
     }
   }
   shiftednodegeom -= nodegeomshift;
   shiftedtimegeom -= timegeomshift;
   
   printf("------------------+------+-------+-------+-------+-------+----------------+----------------+------+--------+-------+-------+-------\n");
   printf("\n");
   printf("------------------------------[Nodes]---------------[Time]------\n");
   printf("  Cnt  Pass  Time  Fail  total(k)     geom.     total     geom. \n");
   printf("----------------------------------------------------------------\n");
   printf("%5d %5d %5d %5d %9d %9.1f %9.1f %9.1f\n",
     nprobs, pass, timeouts, fail, sbab / 1000, nodegeom, stottime, timegeom);
   printf(" shifted geom. [%5d/%5.1f]      %9.1f           %9.1f\n",
     nodegeomshift, timegeomshift, shiftednodegeom, shiftedtimegeom);
   printf("----------------------------------------------------------------\n");
   if( timelimit > 0 )
     printf("@02 timelimit: %g\n", timelimit);
   printf("@01 %s(%s)\n", solver, settings);
   printf("\n");
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   disp_default.c
 * @brief  default display columns
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "disp_default.h"


#define DISP_NAME_SOLFOUND      "solfound"
#define DISP_DESC_SOLFOUND      "letter that indicates the heuristic, that found the solution"
#define DISP_HEAD_SOLFOUND      " "
#define DISP_WIDT_SOLFOUND      1
#define DISP_PRIO_SOLFOUND      80000
#define DISP_POSI_SOLFOUND      0
#define DISP_STRI_SOLFOUND      FALSE

#define DISP_NAME_TIME          "time"
#define DISP_DESC_TIME          "total solution time"
#define DISP_HEAD_TIME          "time"
#define DISP_WIDT_TIME          5
#define DISP_PRIO_TIME          1000
#define DISP_POSI_TIME          50
#define DISP_STRI_TIME          TRUE

#define DISP_NAME_NODENUM       "nodenum"
#define DISP_DESC_NODENUM       "number of processed nodes"
#define DISP_HEAD_NODENUM       "node"
#define DISP_WIDT_NODENUM       7
#define DISP_PRIO_NODENUM       100000
#define DISP_POSI_NODENUM       100
#define DISP_STRI_NODENUM       TRUE

#define DISP_NAME_NODESLEFT     "nodesleft"
#define DISP_DESC_NODESLEFT     "number of unprocessed nodes"
#define DISP_HEAD_NODESLEFT     "left"
#define DISP_WIDT_NODESLEFT     7
#define DISP_PRIO_NODESLEFT     90000
#define DISP_POSI_NODESLEFT     200
#define DISP_STRI_NODESLEFT     TRUE

#define DISP_NAME_LPITERATIONS  "lpiterations"
#define DISP_DESC_LPITERATIONS  "number of simplex iterations"
#define DISP_HEAD_LPITERATIONS  "LP iter"
#define DISP_WIDT_LPITERATIONS  7
#define DISP_PRIO_LPITERATIONS  30000
#define DISP_POSI_LPITERATIONS  1000
#define DISP_STRI_LPITERATIONS  TRUE

#define DISP_NAME_MEMUSED       "memused"
#define DISP_DESC_MEMUSED       "total number of bytes used in block memory"
#define DISP_HEAD_MEMUSED       "mem"
#define DISP_WIDT_MEMUSED       5
#define DISP_PRIO_MEMUSED       20000
#define DISP_POSI_MEMUSED       1500
#define DISP_STRI_MEMUSED       TRUE

#define DISP_NAME_ACTDEPTH      "actdepth"
#define DISP_DESC_ACTDEPTH      "depth of actual node"
#define DISP_HEAD_ACTDEPTH      "depth"
#define DISP_WIDT_ACTDEPTH      5
#define DISP_PRIO_ACTDEPTH      500
#define DISP_POSI_ACTDEPTH      2000
#define DISP_STRI_ACTDEPTH      TRUE

#define DISP_NAME_MAXDEPTH      "maxdepth"
#define DISP_DESC_MAXDEPTH      "maximal depth of all processed nodes"
#define DISP_HEAD_MAXDEPTH      "mdpt"
#define DISP_WIDT_MAXDEPTH      5
#define DISP_PRIO_MAXDEPTH      2000
#define DISP_POSI_MAXDEPTH      2100
#define DISP_STRI_MAXDEPTH      TRUE

#define DISP_NAME_ACTVARS       "actvars"
#define DISP_DESC_ACTVARS       "number of variables in actual node"
#define DISP_HEAD_ACTVARS       "vars"
#define DISP_WIDT_ACTVARS       5
#define DISP_PRIO_ACTVARS       120
#define DISP_POSI_ACTVARS       3000
#define DISP_STRI_ACTVARS       TRUE

#define DISP_NAME_ACTCONSS      "actconss"
#define DISP_DESC_ACTCONSS      "number of enabled constraints in actual node"
#define DISP_HEAD_ACTCONSS      "cons"
#define DISP_WIDT_ACTCONSS      5
#define DISP_PRIO_ACTCONSS      130
#define DISP_POSI_ACTCONSS      3100
#define DISP_STRI_ACTCONSS      TRUE

#define DISP_NAME_ACTCOLS       "actcols"
#define DISP_DESC_ACTCOLS       "number of LP columns in actual node"
#define DISP_HEAD_ACTCOLS       "cols"
#define DISP_WIDT_ACTCOLS       5
#define DISP_PRIO_ACTCOLS       100
#define DISP_POSI_ACTCOLS       3200
#define DISP_STRI_ACTCOLS       TRUE

#define DISP_NAME_ACTROWS       "actrows"
#define DISP_DESC_ACTROWS       "number of LP rows in actual node"
#define DISP_HEAD_ACTROWS       "rows"
#define DISP_WIDT_ACTROWS       5
#define DISP_PRIO_ACTROWS       110
#define DISP_POSI_ACTROWS       3300
#define DISP_STRI_ACTROWS       TRUE

#define DISP_NAME_CUTS          "cuts"
#define DISP_DESC_CUTS          "total number of cuts applied to the LPs"
#define DISP_HEAD_CUTS          "cuts"
#define DISP_WIDT_CUTS          5
#define DISP_PRIO_CUTS          90
#define DISP_POSI_CUTS          3400
#define DISP_STRI_CUTS          TRUE

#define DISP_NAME_SEPAROUNDS    "separounds"
#define DISP_DESC_SEPAROUNDS    "number of separation rounds performed at the actual node"
#define DISP_HEAD_SEPAROUNDS    "sepa"
#define DISP_WIDT_SEPAROUNDS    4
#define DISP_PRIO_SEPAROUNDS    10
#define DISP_POSI_SEPAROUNDS    3500
#define DISP_STRI_SEPAROUNDS    TRUE

#define DISP_NAME_POOLSIZE      "poolsize"
#define DISP_DESC_POOLSIZE      "number of LP rows in the cut pool"
#define DISP_HEAD_POOLSIZE      "pool"
#define DISP_WIDT_POOLSIZE      5
#define DISP_PRIO_POOLSIZE      70
#define DISP_POSI_POOLSIZE      3600
#define DISP_STRI_POOLSIZE      TRUE

#define DISP_NAME_ACTDUALBOUND  "actdualbound"
#define DISP_DESC_ACTDUALBOUND  "dual bound of actual node"
#define DISP_HEAD_ACTDUALBOUND  "actdualbound"
#define DISP_WIDT_ACTDUALBOUND  14
#define DISP_PRIO_ACTDUALBOUND  50
#define DISP_POSI_ACTDUALBOUND  7000
#define DISP_STRI_ACTDUALBOUND  TRUE

#define DISP_NAME_AVGDUALBOUND  "avgdualbound"
#define DISP_DESC_AVGDUALBOUND  "average dual bound of all unprocessed nodes"
#define DISP_HEAD_AVGDUALBOUND  "avgdualbound"
#define DISP_WIDT_AVGDUALBOUND  14
#define DISP_PRIO_AVGDUALBOUND  50000
#define DISP_POSI_AVGDUALBOUND  8000
#define DISP_STRI_AVGDUALBOUND  TRUE

#define DISP_NAME_DUALBOUND     "dualbound"
#define DISP_DESC_DUALBOUND     "current dual bound"
#define DISP_HEAD_DUALBOUND     "dualbound"
#define DISP_WIDT_DUALBOUND     14
#define DISP_PRIO_DUALBOUND     70000
#define DISP_POSI_DUALBOUND     9000
#define DISP_STRI_DUALBOUND     TRUE

#define DISP_NAME_PRIMALBOUND   "primalbound"
#define DISP_DESC_PRIMALBOUND   "current primal bound"
#define DISP_HEAD_PRIMALBOUND   "primalbound"
#define DISP_WIDT_PRIMALBOUND   14
#define DISP_PRIO_PRIMALBOUND   80000
#define DISP_POSI_PRIMALBOUND   10000
#define DISP_STRI_PRIMALBOUND   TRUE

#define DISP_NAME_GAP           "gap"
#define DISP_DESC_GAP           "current relative primal-dual gap"
#define DISP_HEAD_GAP           "gap"
#define DISP_WIDT_GAP           8
#define DISP_PRIO_GAP           60000
#define DISP_POSI_GAP           20000
#define DISP_STRI_GAP           TRUE




/*
 * Callback methods
 */

static
DECL_DISPOUTPUT(SCIPdispOutputSolfound)
{
   SOL* sol;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_SOLFOUND) == 0);
   assert(scip != NULL);

   sol = SCIPgetBestSol(scip);
   if( sol != NULL
      && SCIPsolGetNodenum(sol) == SCIPgetNodenum(scip)
      && SCIPisEQ(scip, SCIPsolGetObj(sol), SCIPgetTransUpperBound(scip)) )
   {
      fprintf(file, "%c", SCIPheurGetDispchar(SCIPsolGetHeur(sol)));
   }
   else
      fprintf(file, " ");

   return SCIP_OKAY;
}

static
DECL_DISPOUTPUT(SCIPdispOutputTime)
{
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_TIME) == 0);
   assert(scip != NULL);

   SCIPdispTime(file, SCIPgetSolvingTime(scip), DISP_WIDT_TIME);

   return SCIP_OKAY;
}

static
DECL_DISPOUTPUT(SCIPdispOutputNodenum)
{
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_NODENUM) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, SCIPgetNodenum(scip), DISP_WIDT_NODENUM);

   return SCIP_OKAY;
}

static
DECL_DISPOUTPUT(SCIPdispOutputNodesleft)
{
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_NODESLEFT) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, SCIPgetNNodesLeft(scip), DISP_WIDT_NODESLEFT);

   return SCIP_OKAY;
}

static
DECL_DISPOUTPUT(SCIPdispOutputLpiterations)
{
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_LPITERATIONS) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, SCIPgetNLPIterations(scip), DISP_WIDT_LPITERATIONS);

   return SCIP_OKAY;
}

static
DECL_DISPOUTPUT(SCIPdispOutputActdepth)
{
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_ACTDEPTH) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, SCIPgetActDepth(scip), DISP_WIDT_ACTDEPTH);

   return SCIP_OKAY;
}

static
DECL_DISPOUTPUT(SCIPdispOutputMemused)
{
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_MEMUSED) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, SCIPgetMemUsed(scip), DISP_WIDT_MEMUSED);

   return SCIP_OKAY;
}

static
DECL_DISPOUTPUT(SCIPdispOutputMaxdepth)
{
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_MAXDEPTH) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, SCIPgetMaxDepth(scip), DISP_WIDT_MAXDEPTH);

   return SCIP_OKAY;
}

static
DECL_DISPOUTPUT(SCIPdispOutputActvars)
{
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_ACTVARS) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, SCIPgetNVars(scip), DISP_WIDT_ACTVARS);

   return SCIP_OKAY;
}

static
DECL_DISPOUTPUT(SCIPdispOutputActconss)
{
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_ACTCONSS) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, SCIPgetNEnabledConss(scip), DISP_WIDT_ACTCONSS);

   return SCIP_OKAY;
}

static
DECL_DISPOUTPUT(SCIPdispOutputActcols)
{
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_ACTCOLS) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, SCIPgetNLPCols(scip), DISP_WIDT_ACTCOLS);

   return SCIP_OKAY;
}

static
DECL_DISPOUTPUT(SCIPdispOutputActrows)
{
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_ACTROWS) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, SCIPgetNLPRows(scip), DISP_WIDT_ACTROWS);

   return SCIP_OKAY;
}

static
DECL_DISPOUTPUT(SCIPdispOutputCuts)
{
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_CUTS) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, SCIPgetNCutsApplied(scip), DISP_WIDT_CUTS);

   return SCIP_OKAY;
}

static
DECL_DISPOUTPUT(SCIPdispOutputSeparounds)
{
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_SEPAROUNDS) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, SCIPgetNSepaRounds(scip), DISP_WIDT_SEPAROUNDS);

   return SCIP_OKAY;
}

static
DECL_DISPOUTPUT(SCIPdispOutputPoolsize)
{
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_POOLSIZE) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, SCIPgetPoolsize(scip), DISP_WIDT_POOLSIZE);

   return SCIP_OKAY;
}

static
DECL_DISPOUTPUT(SCIPdispOutputActdualbound)
{
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_ACTDUALBOUND) == 0);
   assert(scip != NULL);

   fprintf(file, "%13.6e ", SCIPgetActDualBound(scip));

   return SCIP_OKAY;
}

static
DECL_DISPOUTPUT(SCIPdispOutputAvgdualbound)
{
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_AVGDUALBOUND) == 0);
   assert(scip != NULL);

   fprintf(file, "%13.6e ", SCIPgetAvgDualBound(scip));

   return SCIP_OKAY;
}

static
DECL_DISPOUTPUT(SCIPdispOutputDualbound)
{
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_DUALBOUND) == 0);
   assert(scip != NULL);

   fprintf(file, "%13.6e ", SCIPgetDualBound(scip));

   return SCIP_OKAY;
}

static
DECL_DISPOUTPUT(SCIPdispOutputPrimalbound)
{
   Real primalbound;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_PRIMALBOUND) == 0);
   assert(scip != NULL);

   primalbound = SCIPgetPrimalBound(scip);
   if( SCIPisInfinity(scip, ABS(primalbound)) )
      fprintf(file, "      --      ");
   else
      fprintf(file, "%13.6e ", primalbound);
   return SCIP_OKAY;
}

static
DECL_DISPOUTPUT(SCIPdispOutputGap)
{
   Real gap;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_GAP) == 0);
   assert(scip != NULL);

   gap = SCIPgetGap(scip);

   if( SCIPisInfinity(scip, gap) )
      fprintf(file, "    Inf ");
   else if( gap >= 10000.00 )
      fprintf(file, "  Large ");
   else
      fprintf(file, "%7.2f%%", 100.0*gap);

   return SCIP_OKAY;
}




/*
 * default display columns specific interface methods
 */

/** includes the default display columns in SCIP */
RETCODE SCIPincludeDispDefault(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_SOLFOUND, DISP_DESC_SOLFOUND, DISP_HEAD_SOLFOUND,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputSolfound, NULL, 
                  DISP_WIDT_SOLFOUND, DISP_PRIO_SOLFOUND, DISP_POSI_SOLFOUND, DISP_STRI_SOLFOUND) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_TIME, DISP_DESC_TIME, DISP_HEAD_TIME,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputTime, NULL, 
                  DISP_WIDT_TIME, DISP_PRIO_TIME, DISP_POSI_TIME, DISP_STRI_TIME) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_NODENUM, DISP_DESC_NODENUM, DISP_HEAD_NODENUM,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputNodenum, NULL, 
                  DISP_WIDT_NODENUM, DISP_PRIO_NODENUM, DISP_POSI_NODENUM, DISP_STRI_NODENUM) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_NODESLEFT, DISP_DESC_NODESLEFT, DISP_HEAD_NODESLEFT,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputNodesleft, NULL, 
                  DISP_WIDT_NODESLEFT, DISP_PRIO_NODESLEFT, DISP_POSI_NODESLEFT, DISP_STRI_NODESLEFT) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_LPITERATIONS, DISP_DESC_LPITERATIONS, DISP_HEAD_LPITERATIONS,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputLpiterations, NULL, 
                  DISP_WIDT_LPITERATIONS, DISP_PRIO_LPITERATIONS, DISP_POSI_LPITERATIONS, DISP_STRI_LPITERATIONS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_MEMUSED, DISP_DESC_MEMUSED, DISP_HEAD_MEMUSED,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputMemused, NULL, 
                  DISP_WIDT_MEMUSED, DISP_PRIO_MEMUSED, DISP_POSI_MEMUSED, DISP_STRI_MEMUSED) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_ACTDEPTH, DISP_DESC_ACTDEPTH, DISP_HEAD_ACTDEPTH,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputActdepth, NULL, 
                  DISP_WIDT_ACTDEPTH, DISP_PRIO_ACTDEPTH, DISP_POSI_ACTDEPTH, DISP_STRI_ACTDEPTH) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_MAXDEPTH, DISP_DESC_MAXDEPTH, DISP_HEAD_MAXDEPTH,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputMaxdepth, NULL, 
                  DISP_WIDT_MAXDEPTH, DISP_PRIO_MAXDEPTH, DISP_POSI_MAXDEPTH, DISP_STRI_MAXDEPTH) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_ACTVARS, DISP_DESC_ACTVARS, DISP_HEAD_ACTVARS,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputActvars, NULL, 
                  DISP_WIDT_ACTVARS, DISP_PRIO_ACTVARS, DISP_POSI_ACTVARS, DISP_STRI_ACTVARS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_ACTCONSS, DISP_DESC_ACTCONSS, DISP_HEAD_ACTCONSS,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputActconss, NULL, 
                  DISP_WIDT_ACTCONSS, DISP_PRIO_ACTCONSS, DISP_POSI_ACTCONSS, DISP_STRI_ACTCONSS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_ACTCOLS, DISP_DESC_ACTCOLS, DISP_HEAD_ACTCOLS,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputActcols, NULL, 
                  DISP_WIDT_ACTCOLS, DISP_PRIO_ACTCOLS, DISP_POSI_ACTCOLS, DISP_STRI_ACTCOLS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_ACTROWS, DISP_DESC_ACTROWS, DISP_HEAD_ACTROWS,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputActrows, NULL, 
                  DISP_WIDT_ACTROWS, DISP_PRIO_ACTROWS, DISP_POSI_ACTROWS, DISP_STRI_ACTROWS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_CUTS, DISP_DESC_CUTS, DISP_HEAD_CUTS,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputCuts, NULL, 
                  DISP_WIDT_CUTS, DISP_PRIO_CUTS, DISP_POSI_CUTS, DISP_STRI_CUTS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_SEPAROUNDS, DISP_DESC_SEPAROUNDS, DISP_HEAD_SEPAROUNDS,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputSeparounds, NULL, 
                  DISP_WIDT_SEPAROUNDS, DISP_PRIO_SEPAROUNDS, DISP_POSI_SEPAROUNDS, DISP_STRI_SEPAROUNDS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_POOLSIZE, DISP_DESC_POOLSIZE, DISP_HEAD_POOLSIZE,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputPoolsize, NULL, 
                  DISP_WIDT_POOLSIZE, DISP_PRIO_POOLSIZE, DISP_POSI_POOLSIZE, DISP_STRI_POOLSIZE) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_ACTDUALBOUND, DISP_DESC_ACTDUALBOUND, DISP_HEAD_ACTDUALBOUND,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputActdualbound, NULL, 
                  DISP_WIDT_ACTDUALBOUND, DISP_PRIO_ACTDUALBOUND, DISP_POSI_ACTDUALBOUND, DISP_STRI_ACTDUALBOUND) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_AVGDUALBOUND, DISP_DESC_AVGDUALBOUND, DISP_HEAD_AVGDUALBOUND,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputAvgdualbound, NULL, 
                  DISP_WIDT_AVGDUALBOUND, DISP_PRIO_AVGDUALBOUND, DISP_POSI_AVGDUALBOUND, DISP_STRI_AVGDUALBOUND) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_DUALBOUND, DISP_DESC_DUALBOUND, DISP_HEAD_DUALBOUND,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputDualbound, NULL, 
                  DISP_WIDT_DUALBOUND, DISP_PRIO_DUALBOUND, DISP_POSI_DUALBOUND, DISP_STRI_DUALBOUND) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_PRIMALBOUND, DISP_DESC_PRIMALBOUND, DISP_HEAD_PRIMALBOUND,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputPrimalbound, NULL, 
                  DISP_WIDT_PRIMALBOUND, DISP_PRIO_PRIMALBOUND, DISP_POSI_PRIMALBOUND, DISP_STRI_PRIMALBOUND) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_GAP, DISP_DESC_GAP, DISP_HEAD_GAP,
                  SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, SCIPdispOutputGap, NULL, 
                  DISP_WIDT_GAP, DISP_PRIO_GAP, DISP_POSI_GAP, DISP_STRI_GAP) );

   return SCIP_OKAY;
}


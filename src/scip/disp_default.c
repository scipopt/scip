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

#define DISP_NAME_ACTDEPTH      "actdepth"
#define DISP_DESC_ACTDEPTH      "depth of actual node"
#define DISP_HEAD_ACTDEPTH      "depth"
#define DISP_WIDT_ACTDEPTH      5
#define DISP_PRIO_ACTDEPTH      1000
#define DISP_POSI_ACTDEPTH      2000
#define DISP_STRI_ACTDEPTH      TRUE

#define DISP_NAME_MAXDEPTH      "maxdepth"
#define DISP_DESC_MAXDEPTH      "maximal depth of all processed nodes"
#define DISP_HEAD_MAXDEPTH      "mdpt"
#define DISP_WIDT_MAXDEPTH      5
#define DISP_PRIO_MAXDEPTH      2000
#define DISP_POSI_MAXDEPTH      2100
#define DISP_STRI_MAXDEPTH      TRUE

#define DISP_NAME_ACTCOLS       "actcols"
#define DISP_DESC_ACTCOLS       "number of LP columns in actual node"
#define DISP_HEAD_ACTCOLS       "cols"
#define DISP_WIDT_ACTCOLS       6
#define DISP_PRIO_ACTCOLS       100
#define DISP_POSI_ACTCOLS       3000
#define DISP_STRI_ACTCOLS       TRUE

#define DISP_NAME_ACTROWS       "actrows"
#define DISP_DESC_ACTROWS       "number of LP rows in actual node"
#define DISP_HEAD_ACTROWS       "rows"
#define DISP_WIDT_ACTROWS       6
#define DISP_PRIO_ACTROWS       110
#define DISP_POSI_ACTROWS       3100
#define DISP_STRI_ACTROWS       TRUE

#define DISP_NAME_POOLSIZE      "poolsize"
#define DISP_DESC_POOLSIZE      "number of LP rows in the cut pool"
#define DISP_HEAD_POOLSIZE      "pool"
#define DISP_WIDT_POOLSIZE      6
#define DISP_PRIO_POOLSIZE      80
#define DISP_POSI_POOLSIZE      3200
#define DISP_STRI_POOLSIZE      TRUE

#define DISP_NAME_ACTDUALBOUND  "actdualbound"
#define DISP_DESC_ACTDUALBOUND  "dual bound of actual node"
#define DISP_HEAD_ACTDUALBOUND  "actdualbound"
#define DISP_WIDT_ACTDUALBOUND  14
#define DISP_PRIO_ACTDUALBOUND  9000
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
#define DISP_PRIO_GAP           50000
#define DISP_POSI_GAP           20000
#define DISP_STRI_GAP           TRUE




/*
 * Callback methods
 */

static
DECL_DISPOUTP(SCIPdispOutpSolfound)
{
   SOL* sol;
   Longint nodenum;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_SOLFOUND) == 0);
   assert(scip != NULL);

   CHECK_OKAY( SCIPgetBestSol(scip, &sol) );
   CHECK_OKAY( SCIPgetNodenum(scip, &nodenum) );
   if( sol != NULL && SCIPsolGetNodenum(sol) == nodenum )
   {
      fprintf(file, "%c", SCIPheurGetDispchar(SCIPsolGetHeur(sol)));
   }
   else
      fprintf(file, " ");

   return SCIP_OKAY;
}

static
DECL_DISPOUTP(SCIPdispOutpNodenum)
{
   Longint nodenum;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_NODENUM) == 0);
   assert(scip != NULL);

   CHECK_OKAY( SCIPgetNodenum(scip, &nodenum) );
   SCIPdispDecimal(file, nodenum, DISP_WIDT_NODENUM);

   return SCIP_OKAY;
}

static
DECL_DISPOUTP(SCIPdispOutpNodesleft)
{
   int nnodes;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_NODESLEFT) == 0);
   assert(scip != NULL);

   CHECK_OKAY( SCIPgetNNodesLeft(scip, &nnodes) );
   SCIPdispDecimal(file, nnodes, DISP_WIDT_NODESLEFT);

   return SCIP_OKAY;
}

static
DECL_DISPOUTP(SCIPdispOutpLpiterations)
{
   int lpiterations;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_LPITERATIONS) == 0);
   assert(scip != NULL);

   CHECK_OKAY( SCIPgetNLPIterations(scip, &lpiterations) );
   SCIPdispDecimal(file, lpiterations, DISP_WIDT_LPITERATIONS);

   return SCIP_OKAY;
}

static
DECL_DISPOUTP(SCIPdispOutpActdepth)
{
   int actdepth;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_ACTDEPTH) == 0);
   assert(scip != NULL);

   CHECK_OKAY( SCIPgetActDepth(scip, &actdepth) );
   SCIPdispDecimal(file, actdepth, DISP_WIDT_ACTDEPTH);

   return SCIP_OKAY;
}

static
DECL_DISPOUTP(SCIPdispOutpMaxdepth)
{
   int maxdepth;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_MAXDEPTH) == 0);
   assert(scip != NULL);

   CHECK_OKAY( SCIPgetMaxDepth(scip, &maxdepth) );
   SCIPdispDecimal(file, maxdepth, DISP_WIDT_MAXDEPTH);

   return SCIP_OKAY;
}

static
DECL_DISPOUTP(SCIPdispOutpActcols)
{
   int actcols;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_ACTCOLS) == 0);
   assert(scip != NULL);

   CHECK_OKAY( SCIPgetLPCols(scip, NULL, &actcols) );
   SCIPdispDecimal(file, actcols, DISP_WIDT_ACTCOLS);

   return SCIP_OKAY;
}

static
DECL_DISPOUTP(SCIPdispOutpActrows)
{
   int actrows;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_ACTROWS) == 0);
   assert(scip != NULL);

   CHECK_OKAY( SCIPgetLPRows(scip, NULL, &actrows) );
   SCIPdispDecimal(file, actrows, DISP_WIDT_ACTROWS);

   return SCIP_OKAY;
}

static
DECL_DISPOUTP(SCIPdispOutpPoolsize)
{
   int poolsize;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_POOLSIZE) == 0);
   assert(scip != NULL);

   CHECK_OKAY( SCIPgetPoolsize(scip, &poolsize) );
   SCIPdispDecimal(file, poolsize, DISP_WIDT_POOLSIZE);

   return SCIP_OKAY;
}

static
DECL_DISPOUTP(SCIPdispOutpActdualbound)
{
   Real actdualbound;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_ACTDUALBOUND) == 0);
   assert(scip != NULL);

   CHECK_OKAY( SCIPgetActDualBound(scip, &actdualbound) );
   fprintf(file, "%13.6e ", actdualbound);

   return SCIP_OKAY;
}

static
DECL_DISPOUTP(SCIPdispOutpAvgdualbound)
{
   Real avgdualbound;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_AVGDUALBOUND) == 0);
   assert(scip != NULL);

   CHECK_OKAY( SCIPgetAvgDualBound(scip, &avgdualbound) );
   fprintf(file, "%13.6e ", avgdualbound);

   return SCIP_OKAY;
}

static
DECL_DISPOUTP(SCIPdispOutpDualbound)
{
   Real dualbound;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_DUALBOUND) == 0);
   assert(scip != NULL);

   CHECK_OKAY( SCIPgetDualBound(scip, &dualbound) );
   fprintf(file, "%13.6e ", dualbound);

   return SCIP_OKAY;
}

static
DECL_DISPOUTP(SCIPdispOutpPrimalbound)
{
   Real primalbound;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_PRIMALBOUND) == 0);
   assert(scip != NULL);

   CHECK_OKAY( SCIPgetPrimalBound(scip, &primalbound) );
   if( SCIPisInfinity(scip, ABS(primalbound)) )
      fprintf(file, "      --      ");
   else
      fprintf(file, "%13.6e ", primalbound);

   return SCIP_OKAY;
}

static
DECL_DISPOUTP(SCIPdispOutpGap)
{
   Real dualbound;
   Real primalbound;
   Real gap;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_GAP) == 0);
   assert(scip != NULL);

   CHECK_OKAY( SCIPgetDualBound(scip, &dualbound) );
   CHECK_OKAY( SCIPgetPrimalBound(scip, &primalbound) );

   if( SCIPisZero(scip, dualbound) || SCIPisInfinity(scip, ABS(primalbound)) )
      fprintf(file, "    Inf ");
   else
   {
      gap = 100.0 * ABS((dualbound - primalbound)/dualbound);
      
      if( gap >= 10000.00 )
         fprintf(file, "  Large ");
      else
         fprintf(file, "%7.2f%%", gap);
   }

   return SCIP_OKAY;
}




/*
 * default display columns specific interface methods
 */

RETCODE SCIPincludeDispDefault(         /**< includes the default display columns in SCIP */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_SOLFOUND, DISP_DESC_SOLFOUND, DISP_HEAD_SOLFOUND,
                  NULL, NULL, NULL, SCIPdispOutpSolfound, NULL, 
                  DISP_WIDT_SOLFOUND, DISP_PRIO_SOLFOUND, DISP_POSI_SOLFOUND, DISP_STRI_SOLFOUND) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_NODENUM, DISP_DESC_NODENUM, DISP_HEAD_NODENUM,
                  NULL, NULL, NULL, SCIPdispOutpNodenum, NULL, 
                  DISP_WIDT_NODENUM, DISP_PRIO_NODENUM, DISP_POSI_NODENUM, DISP_STRI_NODENUM) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_NODESLEFT, DISP_DESC_NODESLEFT, DISP_HEAD_NODESLEFT,
                  NULL, NULL, NULL, SCIPdispOutpNodesleft, NULL, 
                  DISP_WIDT_NODESLEFT, DISP_PRIO_NODESLEFT, DISP_POSI_NODESLEFT, DISP_STRI_NODESLEFT) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_LPITERATIONS, DISP_DESC_LPITERATIONS, DISP_HEAD_LPITERATIONS,
                  NULL, NULL, NULL, SCIPdispOutpLpiterations, NULL, 
                  DISP_WIDT_LPITERATIONS, DISP_PRIO_LPITERATIONS, DISP_POSI_LPITERATIONS, DISP_STRI_LPITERATIONS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_ACTDEPTH, DISP_DESC_ACTDEPTH, DISP_HEAD_ACTDEPTH,
                  NULL, NULL, NULL, SCIPdispOutpActdepth, NULL, 
                  DISP_WIDT_ACTDEPTH, DISP_PRIO_ACTDEPTH, DISP_POSI_ACTDEPTH, DISP_STRI_ACTDEPTH) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_MAXDEPTH, DISP_DESC_MAXDEPTH, DISP_HEAD_MAXDEPTH,
                  NULL, NULL, NULL, SCIPdispOutpMaxdepth, NULL, 
                  DISP_WIDT_MAXDEPTH, DISP_PRIO_MAXDEPTH, DISP_POSI_MAXDEPTH, DISP_STRI_MAXDEPTH) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_ACTCOLS, DISP_DESC_ACTCOLS, DISP_HEAD_ACTCOLS,
                  NULL, NULL, NULL, SCIPdispOutpActcols, NULL, 
                  DISP_WIDT_ACTCOLS, DISP_PRIO_ACTCOLS, DISP_POSI_ACTCOLS, DISP_STRI_ACTCOLS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_ACTROWS, DISP_DESC_ACTROWS, DISP_HEAD_ACTROWS,
                  NULL, NULL, NULL, SCIPdispOutpActrows, NULL, 
                  DISP_WIDT_ACTROWS, DISP_PRIO_ACTROWS, DISP_POSI_ACTROWS, DISP_STRI_ACTROWS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_POOLSIZE, DISP_DESC_POOLSIZE, DISP_HEAD_POOLSIZE,
                  NULL, NULL, NULL, SCIPdispOutpPoolsize, NULL, 
                  DISP_WIDT_POOLSIZE, DISP_PRIO_POOLSIZE, DISP_POSI_POOLSIZE, DISP_STRI_POOLSIZE) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_ACTDUALBOUND, DISP_DESC_ACTDUALBOUND, DISP_HEAD_ACTDUALBOUND,
                  NULL, NULL, NULL, SCIPdispOutpActdualbound, NULL, 
                  DISP_WIDT_ACTDUALBOUND, DISP_PRIO_ACTDUALBOUND, DISP_POSI_ACTDUALBOUND, DISP_STRI_ACTDUALBOUND) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_AVGDUALBOUND, DISP_DESC_AVGDUALBOUND, DISP_HEAD_AVGDUALBOUND,
                  NULL, NULL, NULL, SCIPdispOutpAvgdualbound, NULL, 
                  DISP_WIDT_AVGDUALBOUND, DISP_PRIO_AVGDUALBOUND, DISP_POSI_AVGDUALBOUND, DISP_STRI_AVGDUALBOUND) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_DUALBOUND, DISP_DESC_DUALBOUND, DISP_HEAD_DUALBOUND,
                  NULL, NULL, NULL, SCIPdispOutpDualbound, NULL, 
                  DISP_WIDT_DUALBOUND, DISP_PRIO_DUALBOUND, DISP_POSI_DUALBOUND, DISP_STRI_DUALBOUND) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_PRIMALBOUND, DISP_DESC_PRIMALBOUND, DISP_HEAD_PRIMALBOUND,
                  NULL, NULL, NULL, SCIPdispOutpPrimalbound, NULL, 
                  DISP_WIDT_PRIMALBOUND, DISP_PRIO_PRIMALBOUND, DISP_POSI_PRIMALBOUND, DISP_STRI_PRIMALBOUND) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_GAP, DISP_DESC_GAP, DISP_HEAD_GAP,
                  NULL, NULL, NULL, SCIPdispOutpGap, NULL, 
                  DISP_WIDT_GAP, DISP_PRIO_GAP, DISP_POSI_GAP, DISP_STRI_GAP) );

   return SCIP_OKAY;
}


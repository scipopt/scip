/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: disp_default.c,v 1.51 2005/03/14 16:10:19 bzfpfend Exp $"

/**@file   disp_default.c
 * @brief  default display columns
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/disp_default.h"


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
#define DISP_PRIO_TIME          4000
#define DISP_POSI_TIME          50
#define DISP_STRI_TIME          TRUE

#define DISP_NAME_NNODES        "nnodes"
#define DISP_DESC_NNODES        "number of processed nodes"
#define DISP_HEAD_NNODES        "node"
#define DISP_WIDT_NNODES        7
#define DISP_PRIO_NNODES        100000
#define DISP_POSI_NNODES        100
#define DISP_STRI_NNODES        TRUE

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

#define DISP_NAME_DEPTH         "depth"
#define DISP_DESC_DEPTH         "depth of current node"
#define DISP_HEAD_DEPTH         "depth"
#define DISP_WIDT_DEPTH         5
#define DISP_PRIO_DEPTH         500
#define DISP_POSI_DEPTH         2000
#define DISP_STRI_DEPTH         TRUE

#define DISP_NAME_MAXDEPTH      "maxdepth"
#define DISP_DESC_MAXDEPTH      "maximal depth of all processed nodes"
#define DISP_HEAD_MAXDEPTH      "mdpt"
#define DISP_WIDT_MAXDEPTH      5
#define DISP_PRIO_MAXDEPTH      5000
#define DISP_POSI_MAXDEPTH      2100
#define DISP_STRI_MAXDEPTH      TRUE

#define DISP_NAME_PLUNGEDEPTH   "plungedepth"
#define DISP_DESC_PLUNGEDEPTH   "current plunging depth"
#define DISP_HEAD_PLUNGEDEPTH   "pdpt"
#define DISP_WIDT_PLUNGEDEPTH   5
#define DISP_PRIO_PLUNGEDEPTH   10
#define DISP_POSI_PLUNGEDEPTH   2200
#define DISP_STRI_PLUNGEDEPTH   TRUE

#define DISP_NAME_NFRAC         "nfrac"
#define DISP_DESC_NFRAC         "number of fractional variables in the current solution"
#define DISP_HEAD_NFRAC         "frac"
#define DISP_WIDT_NFRAC         5
#define DISP_PRIO_NFRAC         40
#define DISP_POSI_NFRAC         2500
#define DISP_STRI_NFRAC         TRUE

#define DISP_NAME_VARS          "vars"
#define DISP_DESC_VARS          "number of variables in the problem"
#define DISP_HEAD_VARS          "vars"
#define DISP_WIDT_VARS          5
#define DISP_PRIO_VARS          3000
#define DISP_POSI_VARS          3000
#define DISP_STRI_VARS          TRUE

#define DISP_NAME_CONSS         "conss"
#define DISP_DESC_CONSS         "number of globally valid constraints in the problem"
#define DISP_HEAD_CONSS         "cons"
#define DISP_WIDT_CONSS         5
#define DISP_PRIO_CONSS         3100
#define DISP_POSI_CONSS         3100
#define DISP_STRI_CONSS         TRUE

#define DISP_NAME_CURCONSS      "curconss"
#define DISP_DESC_CURCONSS      "number of enabled constraints in current node"
#define DISP_HEAD_CURCONSS      "ccons"
#define DISP_WIDT_CURCONSS      5
#define DISP_PRIO_CURCONSS      600
#define DISP_POSI_CURCONSS      3200
#define DISP_STRI_CURCONSS      TRUE

#define DISP_NAME_CURCOLS       "curcols"
#define DISP_DESC_CURCOLS       "number of LP columns in current node"
#define DISP_HEAD_CURCOLS       "cols"
#define DISP_WIDT_CURCOLS       5
#define DISP_PRIO_CURCOLS       800
#define DISP_POSI_CURCOLS       3300
#define DISP_STRI_CURCOLS       TRUE

#define DISP_NAME_CURROWS       "currows"
#define DISP_DESC_CURROWS       "number of LP rows in current node"
#define DISP_HEAD_CURROWS       "rows"
#define DISP_WIDT_CURROWS       5
#define DISP_PRIO_CURROWS       900
#define DISP_POSI_CURROWS       3400
#define DISP_STRI_CURROWS       TRUE

#define DISP_NAME_CUTS          "cuts"
#define DISP_DESC_CUTS          "total number of cuts applied to the LPs"
#define DISP_HEAD_CUTS          "cuts"
#define DISP_WIDT_CUTS          5
#define DISP_PRIO_CUTS          2100
#define DISP_POSI_CUTS          3500
#define DISP_STRI_CUTS          TRUE

#define DISP_NAME_SEPAROUNDS    "separounds"
#define DISP_DESC_SEPAROUNDS    "number of separation rounds performed at the current node"
#define DISP_HEAD_SEPAROUNDS    "sepa"
#define DISP_WIDT_SEPAROUNDS    4
#define DISP_PRIO_SEPAROUNDS    100
#define DISP_POSI_SEPAROUNDS    3600
#define DISP_STRI_SEPAROUNDS    TRUE

#define DISP_NAME_POOLSIZE      "poolsize"
#define DISP_DESC_POOLSIZE      "number of LP rows in the cut pool"
#define DISP_HEAD_POOLSIZE      "pool"
#define DISP_WIDT_POOLSIZE      5
#define DISP_PRIO_POOLSIZE      50
#define DISP_POSI_POOLSIZE      3700
#define DISP_STRI_POOLSIZE      TRUE

#define DISP_NAME_CONFLICTS     "conflicts"
#define DISP_DESC_CONFLICTS     "total number of conflicts found in conflict analysis"
#define DISP_HEAD_CONFLICTS     "confs"
#define DISP_WIDT_CONFLICTS     5
#define DISP_PRIO_CONFLICTS     2000
#define DISP_POSI_CONFLICTS     4000
#define DISP_STRI_CONFLICTS     TRUE

#define DISP_NAME_STRONGBRANCHS "strongbranchs"
#define DISP_DESC_STRONGBRANCHS "total number of strong branching calls"
#define DISP_HEAD_STRONGBRANCHS "strbr"
#define DISP_WIDT_STRONGBRANCHS 5
#define DISP_PRIO_STRONGBRANCHS 1000
#define DISP_POSI_STRONGBRANCHS 5000
#define DISP_STRI_STRONGBRANCHS TRUE

#define DISP_NAME_CURDUALBOUND  "curdualbound"
#define DISP_DESC_CURDUALBOUND  "dual bound of current node"
#define DISP_HEAD_CURDUALBOUND  "curdualbound"
#define DISP_WIDT_CURDUALBOUND  14
#define DISP_PRIO_CURDUALBOUND  400
#define DISP_POSI_CURDUALBOUND  7000
#define DISP_STRI_CURDUALBOUND  TRUE

#define DISP_NAME_AVGDUALBOUND  "avgdualbound"
#define DISP_DESC_AVGDUALBOUND  "average dual bound of all unprocessed nodes"
#define DISP_HEAD_AVGDUALBOUND  "avgdualbound"
#define DISP_WIDT_AVGDUALBOUND  14
#define DISP_PRIO_AVGDUALBOUND  700
#define DISP_POSI_AVGDUALBOUND  8000
#define DISP_STRI_AVGDUALBOUND  TRUE

#define DISP_NAME_DUALBOUND     "dualbound"
#define DISP_DESC_DUALBOUND     "current global dual bound"
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

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputSolfound)
{  /*lint --e{715}*/
   SOL* sol;
   DISPDATA* dispdata;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_SOLFOUND) == 0);
   assert(scip != NULL);

   dispdata = SCIPdispGetData(disp);
   sol = SCIPgetBestSol(scip);
   if( sol != (SOL*)dispdata )
   {
      fprintf(file, "%c", SCIPheurGetDispchar(SCIPgetSolHeur(scip, sol)));
      SCIPdispSetData(disp, (DISPDATA*)sol);
   }
   else
      fprintf(file, " ");

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputTime)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_TIME) == 0);
   assert(scip != NULL);

   SCIPdispTime(file, SCIPgetSolvingTime(scip), DISP_WIDT_TIME);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputNNodes)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_NNODES) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, SCIPgetNNodes(scip), DISP_WIDT_NNODES);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputNodesleft)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_NODESLEFT) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, (Longint)SCIPgetNNodesLeft(scip), DISP_WIDT_NODESLEFT);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputLpiterations)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_LPITERATIONS) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, SCIPgetNLPIterations(scip), DISP_WIDT_LPITERATIONS);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputDepth)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_DEPTH) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, (Longint)SCIPgetDepth(scip), DISP_WIDT_DEPTH);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputMemused)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_MEMUSED) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, SCIPgetMemUsed(scip), DISP_WIDT_MEMUSED);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputMaxdepth)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_MAXDEPTH) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, (Longint)SCIPgetMaxDepth(scip), DISP_WIDT_MAXDEPTH);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputPlungedepth)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_PLUNGEDEPTH) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, (Longint)SCIPgetPlungeDepth(scip), DISP_WIDT_PLUNGEDEPTH);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputNfrac)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_NFRAC) == 0);
   assert(scip != NULL);

   if( SCIPhasCurrentNodeLP(scip) && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
      SCIPdispDecimal(file, (Longint)SCIPgetNLPBranchCands(scip), DISP_WIDT_NFRAC);
   else
      fprintf(file, "   - ");

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputVars)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_VARS) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, (Longint)SCIPgetNVars(scip), DISP_WIDT_VARS);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputConss)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_CONSS) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, (Longint)SCIPgetNGlobalConss(scip), DISP_WIDT_CONSS);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputCurconss)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_CURCONSS) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, (Longint)SCIPgetNEnabledConss(scip), DISP_WIDT_CURCONSS);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputCurcols)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_CURCOLS) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, (Longint)SCIPgetNLPCols(scip), DISP_WIDT_CURCOLS);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputCurrows)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_CURROWS) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, (Longint)SCIPgetNLPRows(scip), DISP_WIDT_CURROWS);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputCuts)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_CUTS) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, (Longint)SCIPgetNCutsApplied(scip), DISP_WIDT_CUTS);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputSeparounds)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_SEPAROUNDS) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, (Longint)SCIPgetNSepaRounds(scip), DISP_WIDT_SEPAROUNDS);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputPoolsize)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_POOLSIZE) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, (Longint)SCIPgetNPoolCuts(scip), DISP_WIDT_POOLSIZE);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputConflicts)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_CONFLICTS) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, SCIPgetNConflictClausesFound(scip), DISP_WIDT_CONFLICTS);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputStrongbranchs)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_STRONGBRANCHS) == 0);
   assert(scip != NULL);

   SCIPdispDecimal(file, SCIPgetNStrongbranchs(scip), DISP_WIDT_STRONGBRANCHS);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputCurdualbound)
{  /*lint --e{715}*/
   Real curdualbound;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_CURDUALBOUND) == 0);
   assert(scip != NULL);

   curdualbound = SCIPgetLocalDualbound(scip);
   if( SCIPisInfinity(scip, REALABS(curdualbound)) )
      fprintf(file, "      --      ");
   else
      fprintf(file, "%13.6e ", curdualbound);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputAvgdualbound)
{  /*lint --e{715}*/
   Real avgdualbound;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_AVGDUALBOUND) == 0);
   assert(scip != NULL);

   avgdualbound = SCIPgetAvgDualbound(scip);
   if( SCIPisInfinity(scip, REALABS(avgdualbound)) )
      fprintf(file, "      --      ");
   else
      fprintf(file, "%13.6e ", avgdualbound);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputDualbound)
{  /*lint --e{715}*/
   Real dualbound;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_DUALBOUND) == 0);
   assert(scip != NULL);

   dualbound = SCIPgetDualbound(scip);
   if( SCIPisInfinity(scip, REALABS(dualbound)) )
      fprintf(file, "      --      ");
   else
      fprintf(file, "%13.6e ", dualbound);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputPrimalbound)
{  /*lint --e{715}*/
   Real primalbound;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_PRIMALBOUND) == 0);
   assert(scip != NULL);

   primalbound = SCIPgetPrimalbound(scip);
   if( SCIPisInfinity(scip, REALABS(primalbound)) )
      fprintf(file, "      --      ");
   else
      fprintf(file, "%13.6e%c", primalbound, SCIPisPrimalboundSol(scip) ? ' ' : '*');

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
DECL_DISPOUTPUT(SCIPdispOutputGap)
{  /*lint --e{715}*/
   Real gap;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_GAP) == 0);
   assert(scip != NULL);

   gap = SCIPgetGap(scip);

   if( SCIPisInfinity(scip, gap) )
      fprintf(file, "    Inf ");
   else if( gap >= 100.00 )
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
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputSolfound, NULL, 
         DISP_WIDT_SOLFOUND, DISP_PRIO_SOLFOUND, DISP_POSI_SOLFOUND, DISP_STRI_SOLFOUND) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_TIME, DISP_DESC_TIME, DISP_HEAD_TIME,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputTime, NULL, 
         DISP_WIDT_TIME, DISP_PRIO_TIME, DISP_POSI_TIME, DISP_STRI_TIME) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_NNODES, DISP_DESC_NNODES, DISP_HEAD_NNODES,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputNNodes, NULL, 
         DISP_WIDT_NNODES, DISP_PRIO_NNODES, DISP_POSI_NNODES, DISP_STRI_NNODES) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_NODESLEFT, DISP_DESC_NODESLEFT, DISP_HEAD_NODESLEFT,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputNodesleft, NULL, 
         DISP_WIDT_NODESLEFT, DISP_PRIO_NODESLEFT, DISP_POSI_NODESLEFT, DISP_STRI_NODESLEFT) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_LPITERATIONS, DISP_DESC_LPITERATIONS, DISP_HEAD_LPITERATIONS,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputLpiterations, NULL, 
         DISP_WIDT_LPITERATIONS, DISP_PRIO_LPITERATIONS, DISP_POSI_LPITERATIONS, DISP_STRI_LPITERATIONS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_MEMUSED, DISP_DESC_MEMUSED, DISP_HEAD_MEMUSED,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputMemused, NULL, 
         DISP_WIDT_MEMUSED, DISP_PRIO_MEMUSED, DISP_POSI_MEMUSED, DISP_STRI_MEMUSED) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_DEPTH, DISP_DESC_DEPTH, DISP_HEAD_DEPTH,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputDepth, NULL, 
         DISP_WIDT_DEPTH, DISP_PRIO_DEPTH, DISP_POSI_DEPTH, DISP_STRI_DEPTH) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_MAXDEPTH, DISP_DESC_MAXDEPTH, DISP_HEAD_MAXDEPTH,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputMaxdepth, NULL, 
         DISP_WIDT_MAXDEPTH, DISP_PRIO_MAXDEPTH, DISP_POSI_MAXDEPTH, DISP_STRI_MAXDEPTH) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_PLUNGEDEPTH, DISP_DESC_PLUNGEDEPTH, DISP_HEAD_PLUNGEDEPTH,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputPlungedepth, NULL, 
         DISP_WIDT_PLUNGEDEPTH, DISP_PRIO_PLUNGEDEPTH, DISP_POSI_PLUNGEDEPTH, DISP_STRI_PLUNGEDEPTH) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_NFRAC, DISP_DESC_NFRAC, DISP_HEAD_NFRAC,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputNfrac, NULL, 
         DISP_WIDT_NFRAC, DISP_PRIO_NFRAC, DISP_POSI_NFRAC, DISP_STRI_NFRAC) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_VARS, DISP_DESC_VARS, DISP_HEAD_VARS,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputVars, NULL, 
         DISP_WIDT_VARS, DISP_PRIO_VARS, DISP_POSI_VARS, DISP_STRI_VARS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_CONSS, DISP_DESC_CONSS, DISP_HEAD_CONSS,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputConss, NULL, 
         DISP_WIDT_CONSS, DISP_PRIO_CONSS, DISP_POSI_CONSS, DISP_STRI_CONSS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_CURCONSS, DISP_DESC_CURCONSS, DISP_HEAD_CURCONSS,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputCurconss, NULL, 
         DISP_WIDT_CURCONSS, DISP_PRIO_CURCONSS, DISP_POSI_CURCONSS, DISP_STRI_CURCONSS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_CURCOLS, DISP_DESC_CURCOLS, DISP_HEAD_CURCOLS,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputCurcols, NULL, 
         DISP_WIDT_CURCOLS, DISP_PRIO_CURCOLS, DISP_POSI_CURCOLS, DISP_STRI_CURCOLS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_CURROWS, DISP_DESC_CURROWS, DISP_HEAD_CURROWS,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputCurrows, NULL, 
         DISP_WIDT_CURROWS, DISP_PRIO_CURROWS, DISP_POSI_CURROWS, DISP_STRI_CURROWS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_CUTS, DISP_DESC_CUTS, DISP_HEAD_CUTS,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputCuts, NULL, 
         DISP_WIDT_CUTS, DISP_PRIO_CUTS, DISP_POSI_CUTS, DISP_STRI_CUTS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_SEPAROUNDS, DISP_DESC_SEPAROUNDS, DISP_HEAD_SEPAROUNDS,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputSeparounds, NULL, 
         DISP_WIDT_SEPAROUNDS, DISP_PRIO_SEPAROUNDS, DISP_POSI_SEPAROUNDS, DISP_STRI_SEPAROUNDS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_POOLSIZE, DISP_DESC_POOLSIZE, DISP_HEAD_POOLSIZE,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputPoolsize, NULL, 
         DISP_WIDT_POOLSIZE, DISP_PRIO_POOLSIZE, DISP_POSI_POOLSIZE, DISP_STRI_POOLSIZE) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_CONFLICTS, DISP_DESC_CONFLICTS, DISP_HEAD_CONFLICTS,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputConflicts, NULL, 
         DISP_WIDT_CONFLICTS, DISP_PRIO_CONFLICTS, DISP_POSI_CONFLICTS, DISP_STRI_CONFLICTS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_STRONGBRANCHS, DISP_DESC_STRONGBRANCHS, DISP_HEAD_STRONGBRANCHS,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputStrongbranchs, NULL, 
         DISP_WIDT_STRONGBRANCHS, DISP_PRIO_STRONGBRANCHS, DISP_POSI_STRONGBRANCHS, DISP_STRI_STRONGBRANCHS) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_CURDUALBOUND, DISP_DESC_CURDUALBOUND, DISP_HEAD_CURDUALBOUND,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputCurdualbound, NULL, 
         DISP_WIDT_CURDUALBOUND, DISP_PRIO_CURDUALBOUND, DISP_POSI_CURDUALBOUND, DISP_STRI_CURDUALBOUND) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_AVGDUALBOUND, DISP_DESC_AVGDUALBOUND, DISP_HEAD_AVGDUALBOUND,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputAvgdualbound, NULL, 
         DISP_WIDT_AVGDUALBOUND, DISP_PRIO_AVGDUALBOUND, DISP_POSI_AVGDUALBOUND, DISP_STRI_AVGDUALBOUND) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_DUALBOUND, DISP_DESC_DUALBOUND, DISP_HEAD_DUALBOUND,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputDualbound, NULL, 
         DISP_WIDT_DUALBOUND, DISP_PRIO_DUALBOUND, DISP_POSI_DUALBOUND, DISP_STRI_DUALBOUND) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_PRIMALBOUND, DISP_DESC_PRIMALBOUND, DISP_HEAD_PRIMALBOUND,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputPrimalbound, NULL, 
         DISP_WIDT_PRIMALBOUND, DISP_PRIO_PRIMALBOUND, DISP_POSI_PRIMALBOUND, DISP_STRI_PRIMALBOUND) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_GAP, DISP_DESC_GAP, DISP_HEAD_GAP,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputGap, NULL, 
         DISP_WIDT_GAP, DISP_PRIO_GAP, DISP_POSI_GAP, DISP_STRI_GAP) );

   return SCIP_OKAY;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
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
#define DISP_WIDT_NODENUM       8
#define DISP_PRIO_NODENUM       100000
#define DISP_POSI_NODENUM       100
#define DISP_STRI_NODENUM       TRUE

#define DISP_NAME_NODESLEFT     "nodesleft"
#define DISP_DESC_NODESLEFT     "number of unprocessed nodes"
#define DISP_HEAD_NODESLEFT     "left"
#define DISP_WIDT_NODESLEFT     8
#define DISP_PRIO_NODESLEFT     90000
#define DISP_POSI_NODESLEFT     200
#define DISP_STRI_NODESLEFT     TRUE

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
   int nodenum;

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
   int nodenum;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_NODENUM) == 0);
   assert(scip != NULL);

   CHECK_OKAY( SCIPgetNodenum(scip, &nodenum) );
   fprintf(file, "%7d ", nodenum);

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
   fprintf(file, "%7d ", nnodes);

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
                  NULL, NULL, SCIPdispOutpSolfound, NULL, 
                  DISP_WIDT_SOLFOUND, DISP_PRIO_SOLFOUND, DISP_POSI_SOLFOUND, DISP_STRI_SOLFOUND) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_NODENUM, DISP_DESC_NODENUM, DISP_HEAD_NODENUM,
                  NULL, NULL, SCIPdispOutpNodenum, NULL, 
                  DISP_WIDT_NODENUM, DISP_PRIO_NODENUM, DISP_POSI_NODENUM, DISP_STRI_NODENUM) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_NODESLEFT, DISP_DESC_NODESLEFT, DISP_HEAD_NODESLEFT,
                  NULL, NULL, SCIPdispOutpNodesleft, NULL, 
                  DISP_WIDT_NODESLEFT, DISP_PRIO_NODESLEFT, DISP_POSI_NODESLEFT, DISP_STRI_NODESLEFT) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_DUALBOUND, DISP_DESC_DUALBOUND, DISP_HEAD_DUALBOUND,
                  NULL, NULL, SCIPdispOutpDualbound, NULL, 
                  DISP_WIDT_DUALBOUND, DISP_PRIO_DUALBOUND, DISP_POSI_DUALBOUND, DISP_STRI_DUALBOUND) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_PRIMALBOUND, DISP_DESC_PRIMALBOUND, DISP_HEAD_PRIMALBOUND,
                  NULL, NULL, SCIPdispOutpPrimalbound, NULL, 
                  DISP_WIDT_PRIMALBOUND, DISP_PRIO_PRIMALBOUND, DISP_POSI_PRIMALBOUND, DISP_STRI_PRIMALBOUND) );
   CHECK_OKAY( SCIPincludeDisp(scip, DISP_NAME_GAP, DISP_DESC_GAP, DISP_HEAD_GAP,
                  NULL, NULL, SCIPdispOutpGap, NULL, 
                  DISP_WIDT_GAP, DISP_PRIO_GAP, DISP_POSI_GAP, DISP_STRI_GAP) );

   return SCIP_OKAY;
}


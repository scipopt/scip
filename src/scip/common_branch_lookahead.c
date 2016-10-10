/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   common_branch_Lookahead.c
 * @brief  LookaheadAbbreviated branching rule
 * @author Christoph Schubert
 */

#include "scip/common_branch_lookahead.h"

/**
 * Allocates buffer memory for the given ValidDomRedData and the contained arrays.
 */
SCIP_RETCODE allocValidBoundData(
   SCIP*                 scip,               /**< SCIP data structure */
   VALIDDOMREDDATA**     validbounddata      /**< The struct to be allocated. */
)
{
   int ntotalvars;
   int i;

   ntotalvars = SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBuffer(scip, validbounddata) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*validbounddata)->upperbounds, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*validbounddata)->lowerbounds, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*validbounddata)->boundstatus, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*validbounddata)->boundedvars, ntotalvars) );

   for( i = 0; i < ntotalvars; i++ )
   {
      (*validbounddata)->boundstatus[i] = 0;
   }
   (*validbounddata)->nboundedvars = 0;

   return SCIP_OKAY;
}

/**
 * Frees the buffer memory of the given ValidDomRedData and the contained arrays.
 */
void freeValidBoundData(
   SCIP*                 scip,               /**< SCIP data structure */
   VALIDDOMREDDATA**     validbounddata      /**< The struct that should be freed. */
)
{
   SCIPfreeBufferArray(scip, &(*validbounddata)->boundedvars);
   SCIPfreeBufferArray(scip, &(*validbounddata)->boundstatus);
   SCIPfreeBufferArray(scip, &(*validbounddata)->lowerbounds);
   SCIPfreeBufferArray(scip, &(*validbounddata)->upperbounds);
   SCIPfreeBuffer(scip, validbounddata);
}

SCIP_RETCODE branchOnVar(
   SCIP*                 scip                /**< SCIP data structure */,
   SCIP_VAR*             var,                /**< the variable to branch on */
   SCIP_Real             val                 /**< the value to branch on */
)
{
   SCIP_NODE* downchild = NULL;
   SCIP_NODE* upchild = NULL;

   SCIPdebugMessage("Effective branching on var <%s> with value <%g>. Old domain: [%g..%g].\n",
      SCIPvarGetName(var), val, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));

   SCIP_CALL( SCIPbranchVarVal(scip, var, val, &downchild, NULL, &upchild) );

   assert(downchild != NULL);
   assert(upchild != NULL);

   return SCIP_OKAY;
}

/**
 * Handles the assignment of ne bounds (valid and supposed). Therefore the new bound is written directly over the old
 * bound. Analog the new bound status is written directly over the old one.
 *
 * @return TRUE, if no upper and lower bound for the given var were yet set; FALSE, else.
 */
SCIP_Bool addBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar,          /**< The variable for which a bound should be added. */
   SCIP_Real             newbound,           /**< The value of the new bound. */
   SCIP_Bool             keepminbound,       /**< In case there is already a bound for the variable, this Indicates
                                              *   whether the min or the max value of the new and the old bound should
                                              *   be kept. */
   BOUNDSTATUS           boundtype,          /**< The type of the new bound. Must be BOUNDSTATUS_UPPERBOUND or
                                              *   BOUNDSTATUS_LOWERBOUND. */
   SCIP_Real*            oldbound,           /**< Pointer to the old bound. Depending on the oldboundstatus this may contain
                                              *   no meaningful data. Also gets the new bound set */
   BOUNDSTATUS*          oldboundstatus      /**< Pointer to the old boundstatus. Also gets the new status set*/
)
{
   SCIP_Bool newboundadded = FALSE;

   assert(scip != NULL);
   assert(branchvar != NULL);
   assert(boundtype == BOUNDSTATUS_LOWERBOUND || boundtype == BOUNDSTATUS_UPPERBOUND);
   assert(oldbound != NULL);
   assert(oldboundstatus != NULL);

   if( *oldboundstatus == boundtype || *oldboundstatus == BOUNDSTATUS_BOTH )
   {
      /* we already have a valid bound with a fitting type set, so we can take min/max of this and the "newbound". */
      if (keepminbound)
      {
         *oldbound = MIN(newbound, *oldbound);
      }
      else
      {
         *oldbound = MAX(newbound, *oldbound);
      }
      SCIPdebugMessage("Updating an existent new bound. var <%s> type <%d> oldbound <%g> newbound <%g>.\n",
         SCIPvarGetName(branchvar), boundtype, *oldbound, newbound);
   }
   else
   {
      /* We either have no new bounds or only a bound with the other type for our var, so we can set the new bound directly. */
      SCIPdebugMessage("Adding new bound. var <%s> type <%d> newbound <%g>.\n", SCIPvarGetName(branchvar), boundtype,
         newbound);
      *oldbound = newbound;

      if( *oldboundstatus == BOUNDSTATUS_NONE )
      {
         *oldboundstatus = boundtype;
         newboundadded = TRUE;
      }
      else
      {
         *oldboundstatus = BOUNDSTATUS_BOTH;
      }
   }
   return newboundadded;
}

/**
 * Adds the given upper bound as a valid bound to the ValidDomRedData container.
 * A valid bound comes from a cutoff on the first level.
 */
void addValidUpperBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar,          /**< the var to assign the new bound to */
   SCIP_Real             newupperbound,      /**< the new upper bound */
   VALIDDOMREDDATA*      validbounds         /**< the container to a add the bound to */
)
{
   int varindex;
   SCIP_Real* oldupperbound;
   BOUNDSTATUS* oldboundstatus;
   SCIP_Bool newboundadded;

   assert(scip != NULL);
   assert(branchvar != NULL);
   assert(validbounds != NULL);

   /* the container is indexed with the probindex */
   varindex = SCIPvarGetProbindex( branchvar );
   oldupperbound = &validbounds->upperbounds[varindex];
   oldboundstatus = &validbounds->boundstatus[varindex];

   /* Add the given bound to the container and keep the min of the new and the old bound.
    * We want to keep the min bound, as the minimum of both valid upper bounds is tighter and still valid. */
   newboundadded = addBound(scip, branchvar, newupperbound, TRUE, BOUNDSTATUS_UPPERBOUND, oldupperbound, oldboundstatus);

   if( newboundadded )
   {
      /* in case of a new entry add the varindex to the container */
      int nboundedvars = validbounds->nboundedvars;
      validbounds->boundedvars[nboundedvars] = varindex;
      validbounds->nboundedvars = nboundedvars + 1;
   }
}

/**
 * Adds the given lower bound as a valid bound to the ValidDomRedData container.
 * A valid bound comes from a cutoff on the first level.
 */
void addValidLowerBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar,          /**< the var to assign the new bound to */
   SCIP_Real             newlowerbound,      /**< the new lower bound */
   VALIDDOMREDDATA*      validbounds         /**< the container to a add the bound to */
)
{
   int varindex;
   SCIP_Real* oldlowerbound;
   BOUNDSTATUS* oldboundstatus;
   SCIP_Bool newboundadded;

   assert(scip != NULL);
   assert(branchvar != NULL);
   assert(validbounds != NULL);

   /* the container is indexed with the probindex */
   varindex = SCIPvarGetProbindex( branchvar );
   oldlowerbound = &validbounds->lowerbounds[varindex];
   oldboundstatus = &validbounds->boundstatus[varindex];

   /* Add the given bound to the container and keep the max of the new and the old bound.
    * We want to keep the max bound, as the maximum of both valid lower bounds is tighter and still valid. */
   newboundadded = addBound(scip, branchvar, newlowerbound, FALSE, BOUNDSTATUS_LOWERBOUND, oldlowerbound, oldboundstatus);

   if( newboundadded )
   {
      /* in case of a new entry add the varindex to the container */
      int nboundedvars = validbounds->nboundedvars;
      validbounds->boundedvars[nboundedvars] = varindex;
      validbounds->nboundedvars = nboundedvars + 1;
   }
}

static const char* names[18] = { "", "SCIP_DIDNOTRUN", "SCIP_DELAYED", "SCIP_DIDNOTFIND", "SCIP_FEASIBLE", "SCIP_INFEASIBLE",
   "SCIP_UNBOUNDED", "SCIP_CUTOFF", "SCIP_SEPARATED", "SCIP_NEWROUND", "SCIP_REDUCEDDOM", "SCIP_CONSADDED",
   "SCIP_CONSCHANGED", "SCIP_BRANCHED", "SCIP_SOLVELP", "SCIP_FOUNDSOL", "SCIP_SUSPENDED", "SCIP_SUCCESS" };

const char* getStatusString(
   SCIP_RESULT           result
)
{
   return names[result];
}

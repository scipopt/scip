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
/*#define SCIP_DEBUG*/
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
   SCIP_CALL( SCIPallocBufferArray(scip, &(*validbounddata)->violatedbybaselp, ntotalvars) );

   for( i = 0; i < ntotalvars; i++ )
   {
      (*validbounddata)->boundstatus[i] = 0;
      (*validbounddata)->violatedbybaselp[i] = FALSE;
   }
   (*validbounddata)->nboundedvars = 0;
   (*validbounddata)->nviolatedbybaselp = 0;

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
   SCIPfreeBufferArray(scip, &(*validbounddata)->violatedbybaselp);
   SCIPfreeBufferArray(scip, &(*validbounddata)->boundedvars);
   SCIPfreeBufferArray(scip, &(*validbounddata)->boundstatus);
   SCIPfreeBufferArray(scip, &(*validbounddata)->lowerbounds);
   SCIPfreeBufferArray(scip, &(*validbounddata)->upperbounds);
   SCIPfreeBuffer(scip, validbounddata);
}

SCIP_RETCODE branchOnVar(
   SCIP*                 scip                /**< SCIP data structure */,
   SCIP_VAR*             var,                /**< the variable to branch on */
   SCIP_Real             val,                /**< the value to branch on */
   SCIP_Real             bestdown,
   SCIP_Bool             bestdownvalid,
   SCIP_Real             bestup,
   SCIP_Real             bestupvalid,
   SCIP_Real             provedbound
)
{
   SCIP_NODE* downchild = NULL;
   SCIP_NODE* upchild = NULL;

   SCIPdebugMessage("Effective branching on var <%s> with value <%g>. Old domain: [%g..%g].\n",
      SCIPvarGetName(var), val, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));

   SCIP_CALL( SCIPbranchVarVal(scip, var, val, &downchild, NULL, &upchild) );

   assert(downchild != NULL);
   assert(upchild != NULL);

   /* update the lower bounds in the children; we must not do this if columns are missing in the LP
    * (e.g., because we are doing branch-and-price) or the problem should be solved exactly
    */
   if( SCIPallColsInLP(scip) && !SCIPisExactSolve(scip) )
   {
      SCIP_CALL( SCIPupdateNodeLowerbound(scip, downchild, bestdownvalid ? MAX(bestdown, provedbound) : provedbound) );
      SCIP_CALL( SCIPupdateNodeLowerbound(scip, upchild, bestupvalid ? MAX(bestup, provedbound) : provedbound) );
   }
   SCIPdebugMessage(" -> down child's lowerbound: %g\n", SCIPnodeGetLowerbound(downchild));
   SCIPdebugMessage(" -> up child's lowerbound: %g\n", SCIPnodeGetLowerbound(upchild));

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
   SCIP_Real             baselpsolval,       /**< the lp solution of the base node */
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
      int nboundedvars;

      /* in case of a new entry add the varindex to the container */
      nboundedvars = validbounds->nboundedvars;
      validbounds->boundedvars[nboundedvars] = varindex;
      validbounds->nboundedvars = nboundedvars + 1;

      if( !validbounds->violatedbybaselp[varindex] && SCIPisFeasGT(scip, baselpsolval, newupperbound) )
      {
         validbounds->violatedbybaselp[varindex] = TRUE;
         validbounds->nviolatedbybaselp++;
      }
   }
}

/**
 * Adds the given lower bound as a valid bound to the ValidDomRedData container.
 * A valid bound comes from a cutoff on the first level.
 */
void addValidLowerBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             baselpsolval,       /**< the lp solution of the base node */
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
      int nboundedvars;

      /* in case of a new entry add the varindex to the container */
      nboundedvars = validbounds->nboundedvars;

      validbounds->boundedvars[nboundedvars] = varindex;
      validbounds->nboundedvars = nboundedvars + 1;

      if( !validbounds->violatedbybaselp[varindex] && SCIPisFeasLT(scip, baselpsolval, newlowerbound) )
      {
         validbounds->violatedbybaselp[varindex] = TRUE;
         validbounds->nviolatedbybaselp++;
      }
   }
}

/**
 * Adds the domain reductions found throughout the execution of the branching rule.
 * Domain reductions of a variable occur if:
 * - one branch on the first level is cutoff (directly or because both branches of a second level variable were cutoff)
 * - both second level branches in the same direction for the same first level variable are cutoff
 *
 * @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 * @see VALIDDOMREDDATA
 */
SCIP_RETCODE addDomainReductions(
   SCIP*                 scip,               /**< SCIP data structure */
   VALIDDOMREDDATA*      validbounds,        /**< The struct containing all bounds that should be added. */
   SCIP_Bool*            domredcutoff,
   SCIP_Bool*            domred
   )
{
   int i;
   int nboundedvars;
   SCIP_VAR** probvars;
   int nboundsadded = 0;

   assert(scip != NULL);
   assert(validbounds != NULL);

   nboundedvars = validbounds->nboundedvars;
   probvars = SCIPgetVars(scip);

   /* Instead of iterating over all problem vars, we iterate only over those that have a new bound set. */
   SCIPdebugMessage("Trying to add domain reductions for <%d> variables.\n", nboundedvars);
   for( i = 0; i < nboundedvars; i++ )
   {
      int probvarindex;
      BOUNDSTATUS boundstatus;
      SCIP_VAR* branchvar;

      /* The iteration index is only valid as an index for the boundedvars array. Every other array in validbounds is
       * indexed by the probindex contained in boundedvars. */
      probvarindex = validbounds->boundedvars[i];
      boundstatus = validbounds->boundstatus[probvarindex];
      branchvar = probvars[probvarindex];

      /* Handle the new lower bound (if present) */
      if( !*domredcutoff && (boundstatus == BOUNDSTATUS_LOWERBOUND || boundstatus == BOUNDSTATUS_BOTH) )
      {
         SCIP_Bool infeasible;
         SCIP_Bool tightened;
         SCIP_Real oldlowerbound;
         SCIP_Real proposedlowerbound;
         SCIP_Real newlowerbound;

         /* get the old and the new lower bound */
         oldlowerbound = SCIPvarGetLbLocal(branchvar);
         proposedlowerbound = validbounds->lowerbounds[probvarindex];

         /* add the new bound */
         SCIP_CALL( SCIPtightenVarLb(scip, branchvar, proposedlowerbound, FALSE, &infeasible, &tightened) );

         newlowerbound = SCIPvarGetLbLocal(branchvar);
         SCIPdebugMessage("Variable <%s>, old lower bound <%g>, proposed lower bound <%g>, new lower bound <%g>\n",
            SCIPvarGetName(branchvar), oldlowerbound, proposedlowerbound, newlowerbound);

         if( infeasible )
         {
            /* the domain reduction may result in an empty model (ub < lb) */
            *domredcutoff = TRUE;
            SCIPdebugMessage("The domain reduction of variable <%s> resulted in an empty model.\n",
               SCIPvarGetName(branchvar));
         }
         else if( tightened )
         {
            /* the lb is now strictly greater than before */
            *domred = TRUE;
            SCIPdebugMessage("The lower bound of variable <%s> was successfully tightened.\n", SCIPvarGetName(branchvar));
            nboundsadded++;
         }
      }

      /* Handle the new upper bound (if present) */
      if( !*domredcutoff && (boundstatus == BOUNDSTATUS_UPPERBOUND || boundstatus == BOUNDSTATUS_BOTH) )
      {
         SCIP_Bool infeasible;
         SCIP_Bool tightened;
         SCIP_Real oldupperbound;
         SCIP_Real proposedupperbound;
         SCIP_Real newupperbound;

         /* get the old and the new upper bound */
         oldupperbound = SCIPvarGetUbLocal(branchvar);
         proposedupperbound = validbounds->upperbounds[probvarindex];

         /* add the new bound */
         SCIP_CALL( SCIPtightenVarUb(scip, branchvar, proposedupperbound, FALSE, &infeasible, &tightened) );

         newupperbound = SCIPvarGetUbLocal(branchvar);
         SCIPdebugMessage("Variable <%s>, old upper bound <%g>, proposed upper bound <%g>, new upper bound <%g>\n",
            SCIPvarGetName(branchvar), oldupperbound, proposedupperbound, newupperbound);

         if( infeasible )
         {
            /* the domain reduction may result in an empty model (ub < lb) */
            *domredcutoff = TRUE;
            SCIPdebugMessage("The upper bound of variable <%s> could not be tightened.\n", SCIPvarGetName(branchvar));
         }
         else if( tightened )
         {
            /* the ub is now strictly smaller than before */
            *domred = TRUE;
            SCIPdebugMessage("The upper bound of variable <%s> was successfully tightened.\n", SCIPvarGetName(branchvar));
            nboundsadded++;
         }
      }

      if( boundstatus != BOUNDSTATUS_NONE )
      {
         /* Reset the array s.t. it only contains zero values. Necessary for the CleanBufferArray usage. */
         validbounds->boundstatus[probvarindex] = BOUNDSTATUS_NONE;
      }
   }

   SCIPdebugMessage("Added <%d> real domain reductions to the problem.\n", nboundsadded);

   return SCIP_OKAY;
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

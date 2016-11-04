/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_completesol.c
 * @brief  COMPLETESOL - primal heuristic trying to complete given partial solutions
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "scip/heur_completesol.h"
#include "scip/scipdefplugins.h"       /* needed for the secondary SCIP instance */
#include "scip/pub_misc.h"
#include "scip/def.h"

#define HEUR_NAME             "completesol"
#define HEUR_DESC             "primal heuristic trying to complete given partial solutions"
#define HEUR_DISPCHAR         'h'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         0
#define HEUR_TIMING           SCIP_HEURTIMING_BEFOREPRESOL
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

/* default values for heuristic plugins */
#define DEFAULT_MAXNODES      5000LL    /**< maximum number of nodes to regard in the subproblem */
#define DEFAULT_MAXUNKRATE    0.85      /**< maximum percentage of unknown solution values */
#define DEFAULT_ADDALLSOLS   FALSE      /**< should all subproblem solutions be added to the original SCIP? */
#define DEFAULT_MINNODES      50LL      /**< minimum number of nodes to regard in the subproblem */
#define DEFAULT_NODESOFS      500LL     /**< number of nodes added to the contingent of the total nodes */
#define DEFAULT_NODESQUOT     0.1       /**< subproblem nodes in relation to nodes of the original problem */
#define DEFAULT_LPLIMFAC      2.0       /**< factor by which the limit on the number of LP depends on the node limit */
#define DEFAULT_OBJWEIGHT     1.0       /**< weight of the original objective function (1: only original objective) */
#define DEFAULT_MINIMPROVE    0.01      /**< factor by which the incumbent should be improved at least */
#define DEFAULT_MINOBJWEIGHT 1e-3       /**< minimal weight for original objective function (zero could lead to infinite solutions) */
#define DEFAULT_IGNORECONT  FALSE       /**< should solution values for continuous variables be ignored? */
#define DEFAULT_BESTSOLS        5       /**< heuristic stops, if the given number of improving solutions were found (-1: no limit) */
#define DEFAULT_MAXPROPROUNDS  10       /**< maximal number of iterations in proagation (-1: no limit) */

/* event handler properties */
#define EVENTHDLR_NAME         "Completesol"
#define EVENTHDLR_DESC         "LP event handler for " HEUR_NAME " heuristic"


/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Real             maxunknownrate;     /**< maximal rate of changed coefficients in the objective function */
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem */
   SCIP_Bool             addallsols;         /**< should all subproblem solutions be added to the original SCIP? */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem */
   SCIP_Real             nodelimit;          /**< the nodelimit employed in the current sub-SCIP, for the event handler*/
   SCIP_Real             lplimfac;           /**< factor by which the limit on the number of LP depends on the node limit */
   SCIP_Real             objweight;          /**< weight of the original objective function (1: only original obj, 0: try to keep to given solution) */
   SCIP_Real             minimprove;         /**< factor by which the incumbent should be improved at least */
   SCIP_Bool             ignorecont;         /**< should solution values for continuous variables be ignored? */
   int                   bestsols;           /**< heuristic stops, if the given number of improving solutions were found (-1: no limit) */
   int                   maxproprounds;      /**< maximal number of iterations in proagation (-1: no limit) */
};

/* ---------------- Callback methods of event handler ---------------- */

/* exec the event handler
 *
 * we interrupt the solution process
 */
static
SCIP_DECL_EVENTEXEC(eventExecCompletesol)
{
   SCIP_HEURDATA* heurdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_LPSOLVED);

   heurdata = (SCIP_HEURDATA*)eventdata;
   assert(heurdata != NULL);

   /* interrupt solution process of sub-SCIP */
   if( SCIPgetNLPs(scip) > heurdata->lplimfac * heurdata->nodelimit )
   {
      SCIPdebugMsg(scip, "interrupt after %" SCIP_LONGINT_FORMAT " LPs\n",SCIPgetNLPs(scip));
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   return SCIP_OKAY;
}

/** creates a subproblem by fixing a number of variables */
static
SCIP_RETCODE createSubproblem(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem */
   SCIP_HEURDATA*        heurdata,           /**< heuristic's private data structure */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_SOL*             partialsol,         /**< partial solution */
   SCIP_Bool*            tightened,          /**< array to store for which variables we have found bound tightenings */
   SCIP_Bool*            success             /**< pointer to store whether the creation was successful */
   )
{
   SCIP_VAR** vars;
   SCIP_CONS* objcons;
   SCIP_Real epsobj;
   SCIP_Real cutoff;
   SCIP_Real upperbound;
   char consobjname[SCIP_MAXSTRLEN];
   int nvars;
   int i;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(heurdata != NULL);

   *success = TRUE;

   /* if there is already a solution, add an objective cutoff */
   if( SCIPgetNSols(scip) > 0 )
   {
      cutoff = SCIPinfinity(scip);
      assert(!SCIPisInfinity(scip, SCIPgetUpperbound(scip)));

      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);

      if( !SCIPisInfinity(scip, -1.0 * SCIPgetLowerbound(scip)) )
         cutoff = (1 - heurdata->minimprove) * SCIPgetUpperbound(scip) + heurdata->minimprove * SCIPgetLowerbound(scip);
      else
      {
         if( SCIPgetUpperbound(scip) >= 0 )
            cutoff = (1 - heurdata->minimprove) * SCIPgetUpperbound(scip);
         else
            cutoff = (1 + heurdata->minimprove) * SCIPgetUpperbound(scip);
      }
      cutoff = MIN(upperbound, cutoff);
      SCIPdebugMsg(scip, "set cutoff=%g for sub-SCIP\n", cutoff);
   }
   else
      cutoff = SCIPinfinity(scip);

   /* calculate objective coefficients for all potential epsilons */
   if( SCIPisEQ(scip, heurdata->objweight, 1.0) )
      return SCIP_OKAY;
   else if( !SCIPisInfinity(scip, cutoff) )
      epsobj = 1.0;
   else
   {
      /* divide by objweight to avoid changing objective coefficient of original problem variables */
      epsobj = (1.0 - heurdata->objweight)/heurdata->objweight;

      /* scale with -1 if we have a maximization problem */
      if( SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE )
         epsobj *= -1.0;
   }

   /* get active variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   objcons = NULL;

   /* add constraints to messure the distance to the given partial solution */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_Real solval;
      int idx;

      assert(SCIPvarIsActive(vars[i]));

      /* add objective function as a constraint, if a primal bound exists */
      if( SCIPisInfinity(scip, cutoff) )
      {
         /* create the constraints */
         if( objcons == NULL )
         {
            SCIP_Real lhs;
            SCIP_Real rhs;

            if( SCIPgetObjsense(subscip) == SCIP_OBJSENSE_MINIMIZE )
            {
               lhs = -SCIPinfinity(subscip);
               rhs = cutoff;
            }
            else
            {
               lhs = cutoff;
               rhs = SCIPinfinity(subscip);
            }

            (void)SCIPsnprintf(consobjname, SCIP_MAXSTRLEN, "obj");
            SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &objcons, consobjname, 0, NULL, NULL, lhs, rhs) );
         }

         /* add the variable to the constraints */
         SCIP_CALL( SCIPaddCoefLinear(subscip, objcons, subvars[i], SCIPvarGetObj(subvars[i])) );

         /* set objective coefficient to 0.0 */
         SCIP_CALL( SCIPchgVarObj(subscip, subvars[i], 0.0) );
      }

      solval = SCIPgetSolVal(scip, partialsol, vars[i]);

      /* skip variables with unknown solution value */
      if( solval == SCIP_UNKNOWN ) /*lint !e777*/
         continue;

      idx = SCIPvarGetProbindex(vars[i]);
      assert(idx >= 0);

      /* skip variables where we already found some bound tightenings */
      if( tightened[idx] == FALSE )
      {
         /* special case: vars[i] is binary; we do not add an extra variable, but we mimic the behaviour we would get with it.
          * E.g., if the solval is 0.3, setting the variable to 0 would give a cost of 0.3 * epsobj, setting it to 1 gives
          * 0.7 * epsobj. Thus, 0.3 * epsobj can be treated as a constant in the objective function and the variable gets
          * an objective coefficient of 0.4 * epsobj.
          */
         if( SCIPvarIsBinary(vars[i]) )
         {
            SCIP_Real frac = SCIPfeasFrac(scip, solval);
            SCIP_Real objcoef;

            frac = MIN(frac, 1-frac);
            objcoef = (1 - 2*frac) * epsobj * (int)SCIPgetObjsense(scip);

            if( solval > 0.5 )
            {
               SCIP_CALL( SCIPchgVarObj(scip, vars[i], -objcoef) );
            }
            else
            {
               SCIP_CALL( SCIPchgVarObj(scip, vars[i], objcoef) );
            }
         }
         else
         {
            SCIP_CONS* conspos;
            SCIP_CONS* consneg;
            SCIP_VAR* eps;
            char consnamepos[SCIP_MAXSTRLEN];
            char consnameneg[SCIP_MAXSTRLEN];
            char epsname[SCIP_MAXSTRLEN];

            /* create two new variables */
            (void)SCIPsnprintf(epsname, SCIP_MAXSTRLEN, "eps_%s", SCIPvarGetName(subvars[i]));

            SCIP_CALL( SCIPcreateVarBasic(subscip, &eps, epsname, 0.0, SCIPinfinity(scip), epsobj, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(subscip, eps) );

            /* create two constraints */
            (void)SCIPsnprintf(consnamepos, SCIP_MAXSTRLEN, "cons_%s_pos", SCIPvarGetName(subvars[i]));
            (void)SCIPsnprintf(consnameneg, SCIP_MAXSTRLEN, "cons_%s_neq", SCIPvarGetName(subvars[i]));

            /* x_{i} - s_{i} <= e_{i}   <==>   x_{i} - e_{i} <= s_{i} */
            SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &conspos, consnamepos, 0, NULL, NULL, -SCIPinfinity(scip), solval) );
            SCIP_CALL( SCIPaddCoefLinear(subscip, conspos, subvars[i], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(subscip, conspos, eps, -1.0) );
            SCIP_CALL( SCIPaddCons(subscip, conspos) );
            SCIP_CALL( SCIPreleaseCons(subscip, &conspos) );

            /* s_{i} - x_{i} <= e_{i}   <==>   e_{i} - x_{i} >= s_{i} */
            SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &consneg, consnameneg, 0, NULL, NULL, solval, SCIPinfinity(scip)) );
            SCIP_CALL( SCIPaddCoefLinear(subscip, consneg, subvars[i], -1.0) );
            SCIP_CALL( SCIPaddCoefLinear(subscip, consneg, eps, 1.0) );
            SCIP_CALL( SCIPaddCons(subscip, consneg) );
            SCIP_CALL( SCIPreleaseCons(subscip, &consneg) );

            /* release the variables */
            SCIP_CALL( SCIPreleaseVar(subscip, &eps) );
         }
      }
   }

   /* add and release the constraint representing the original objective function */
   if( objcons != NULL )
   {
      SCIP_CALL( SCIPaddCons(subscip, objcons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &objcons) );
   }

   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_HEUR*            heur,               /**< Completesol heuristic structure */
   SCIP_SOL*             subsol,             /**< solution of the subproblem or the partial */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables */
   int nvars;                                /* the original problem's number of variables */
   SCIP_SOL* newsol;                         /* solution to be created for the original problem */
   int v;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */
   assert(nvars <= SCIPgetNOrigVars(subscip));

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );

   for( v = 0; v < nvars; v++ )
   {
      SCIP_Real solval = SCIPgetSolVal(subscip, subsol, subvars[v]);

      assert(!SCIPisInfinity(subscip, solval) && !SCIPisInfinity(subscip, -solval));
      assert(solval != SCIP_UNKNOWN); /*lint !e777*/

      SCIP_CALL( SCIPsetSolVal(scip, newsol, vars[v], solval) );
   }

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, success) );

   return SCIP_OKAY;
}

/** perform a probing bound change or fixes the variable */
static
SCIP_RETCODE chgProbingBound(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             newval,             /**< new bound */
   SCIP_BRANCHDIR        branchdir           /**< bound change direction */
   )
{
   SCIP_Real ub;
   SCIP_Real lb;

   assert(scip != NULL);
   assert(var != NULL);

   ub = SCIPvarGetUbLocal(var);
   lb = SCIPvarGetLbLocal(var);

   switch (branchdir) {
   case SCIP_BRANCHDIR_DOWNWARDS:
      if( SCIPisLT(scip, newval, ub) && SCIPisGE(scip, newval, lb) )
      {
         SCIP_CALL( SCIPchgVarUbProbing(scip, var, newval) );
      }
      break;
   case SCIP_BRANCHDIR_UPWARDS:
      if( SCIPisLE(scip, newval, ub) && SCIPisGT(scip, newval, lb) )
      {
         SCIP_CALL( SCIPchgVarLbProbing(scip, var, newval) );
      }
      break;
   case SCIP_BRANCHDIR_FIXED:
      if( SCIPisLE(scip, newval, ub) && SCIPisGE(scip, newval, lb) )
      {
         SCIP_CALL( SCIPfixVarProbing(scip, var, newval) );
      }
      break;
   default:
      return SCIP_INVALIDDATA;
   }/*lint !e788*/

   return SCIP_OKAY;
}

/** tries variables bound changes guided by the given solution */
static
SCIP_RETCODE tightenVariables(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic's private data structure */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nvars,              /**< number of problem variables */
   SCIP_SOL*             sol,                /**< solution to guide the bound changes */
   SCIP_Bool*            tightened           /**< array to store if variable bound could be tightened */
   )
{
#ifndef NDEBUG
   SCIP_Bool incontsection;
#endif
   SCIP_Bool cutoff;
   SCIP_Longint ndomreds;
   SCIP_Longint ndomredssum;
   int nbndtightenings;
   int v;

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(vars != NULL);
   assert(nvars >= 0);
   assert(sol != NULL);
   assert(tightened != NULL);

   assert(SCIPsolGetOrigin(sol) == SCIP_SOLORIGIN_PARTIAL);

   SCIPdebugMsg(scip, "> start probing along the solution values\n");

   nbndtightenings = 0;
   ndomredssum = 0;
#ifndef NDEBUG
   incontsection = FALSE;
#endif

   /* there is at least one integral variable; open one probing node for all non-continuous variables */
   if( nvars - SCIPgetNContVars(scip) > 0 )
   {
      SCIP_CALL( SCIPnewProbingNode(scip) );
   }

   for( v = 0; v < nvars; v++ )
   {
      SCIP_Real solval;

      assert(SCIPvarIsActive(vars[v]));

#ifndef NDEBUG
      incontsection |= (!SCIPvarIsIntegral(vars[v])); /*lint !e514*/
      assert(!incontsection || !SCIPvarIsIntegral(vars[v]));
#endif

      /* return if continuous variables should ignored */
      if( heurdata->ignorecont && SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS )
         break;

      /* return if we have found enough domain reductions tightenings */
      if( ndomredssum > 0.3*nvars )
         break;

      solval = SCIPgetSolVal(scip, sol, vars[v]);

      /* skip unknows variables */
      if( solval == SCIP_UNKNOWN ) /*lint !e777*/
         continue;
      assert(!SCIPisInfinity(scip, solval) && !SCIPisInfinity(scip, -solval));

      cutoff = FALSE;
      ndomreds = 0;

      /* variable is binary or integer */
      if( SCIPvarIsIntegral(vars[v]) )
      {
         /* the solution value is integral, try to fix them */
         if( SCIPisIntegral(scip, solval) )
         {
            SCIP_CALL( chgProbingBound(scip, vars[v], solval, SCIP_BRANCHDIR_FIXED) );
            tightened[SCIPvarGetProbindex(vars[v])] = TRUE;
            ++nbndtightenings;

#ifdef SCIP_MORE_DEBUG
               SCIPdebugMsg(scip, "> fix variable <%s> = [%g,%g] to %g \n", SCIPvarGetName(vars[v]),
                     SCIPvarGetLbGlobal(vars[v]), SCIPvarGetUbGlobal(vars[v]), solval);
#endif
         }
         else
         {
            SCIP_Real ub = SCIPceil(scip, solval) + 1.0;
            SCIP_Real lb = SCIPfloor(scip, solval) - 1.0;

            /* try tightening of upper bound */
            if( SCIPisLT(scip, ub, SCIPvarGetUbLocal(vars[v])) )
            {
               SCIP_CALL( chgProbingBound(scip, vars[v], solval, SCIP_BRANCHDIR_DOWNWARDS) );
               tightened[SCIPvarGetProbindex(vars[v])] = TRUE;
               ++nbndtightenings;

#ifdef SCIP_MORE_DEBUG
               SCIPdebugMsg(scip, "> tighten upper bound of variable <%s>: %g to %g\n", SCIPvarGetName(vars[v]),
                     SCIPvarGetUbGlobal(vars[v]), ub);
#endif
            }

            /* try tightening of lower bound */
            if( SCIPisGT(scip, lb, SCIPvarGetLbLocal(vars[v])) )
            {
               SCIP_CALL( chgProbingBound(scip, vars[v], solval, SCIP_BRANCHDIR_UPWARDS) );
               tightened[SCIPvarGetProbindex(vars[v])] = TRUE;
               ++nbndtightenings;

#ifdef SCIP_MORE_DEBUG
               SCIPdebugMsg(scip, "> tighten lower bound of variable <%s>: %g to %g\n", SCIPvarGetName(vars[v]),
                     SCIPvarGetLbGlobal(vars[v]), ub);
#endif
            }
         }
      }
      /* variable is continuous */
      else
      {
         /* fix to lb or ub */
         if( SCIPisEQ(scip, solval, SCIPvarGetLbLocal(vars[v])) || SCIPisEQ(scip, solval, SCIPvarGetUbLocal(vars[v])) )
         {
            /* open a new probing node */
            if( SCIPgetProbingDepth(scip) < SCIP_MAXTREEDEPTH - 1 )
            {
               SCIP_CALL( SCIPnewProbingNode(scip) );
            }
            SCIP_CALL( chgProbingBound(scip, vars[v], solval, SCIP_BRANCHDIR_FIXED) );

            SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, &cutoff, &ndomreds) );

            if( cutoff || ndomreds == 0 )
            {
               SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip)-1) );
            }
            else
            {
               assert(SCIPvarGetProbindex(vars[v]) >= 0);
               tightened[SCIPvarGetProbindex(vars[v])] = TRUE;
               ++nbndtightenings;
#ifdef SCIP_MORE_DEBUG
               SCIPdebugMsg(scip, "> fix variable <%s> = [%g,%g] to %g (ndomreds=%lld)\n", SCIPvarGetName(vars[v]),
                     SCIPvarGetLbGlobal(vars[v]), SCIPvarGetUbGlobal(vars[v]), solval, ndomreds);
#endif
            }
         }
         else
         {
            SCIP_Real offset;
            SCIP_Real ub = SCIPvarGetUbGlobal(vars[v]);
            SCIP_Real lb = SCIPvarGetLbGlobal(vars[v]);

            /* both bound are finite; add 10% of domain range to the new bounds */
            if( !SCIPisInfinity(scip, -lb) && SCIPisInfinity(scip, ub) )
               offset = REALABS(0.1 * (ub-lb));
            else
            {
               /* if one bound is finite, add 10% of range between solval and finite bound the new bounds */
               if( !SCIPisInfinity(scip, -lb) )
                  offset = REALABS(0.1 * (solval-lb));
               else
               {
                  assert(!SCIPisInfinity(scip, ub));
                  offset = REALABS(0.1 * (ub-solval));
               }
            }

            /* update bounds */
            ub = SCIPceil(scip, solval) + offset;
            lb = SCIPfloor(scip, solval) - offset;

            /* try tightening of upper bound */
            if( SCIPisLT(scip, ub, SCIPvarGetUbLocal(vars[v])) )
            {
               /* open a new probing node */
               if( SCIPgetProbingDepth(scip) < SCIP_MAXTREEDEPTH-10 )
               {
                  SCIP_CALL( SCIPnewProbingNode(scip) );
               }
               SCIP_CALL( chgProbingBound(scip, vars[v], solval, SCIP_BRANCHDIR_DOWNWARDS) );

               SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, &cutoff, &ndomreds) );

               if( cutoff || ndomreds == 0 )
               {
                  SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip)-1) );
               }
               else
               {
                  assert(SCIPvarGetProbindex(vars[v]) >= 0);
                  tightened[SCIPvarGetProbindex(vars[v])] = TRUE;
                  ++nbndtightenings;
#ifdef SCIP_MORE_DEBUG
                  SCIPdebugMsg(scip, "> tighten upper bound of variable <%s>: %g to %g (ndomreds=%lld)\n",
                        SCIPvarGetName(vars[v]), SCIPvarGetUbGlobal(vars[v]), ub, ndomreds);
#endif
               }
            }

            /* try tightening of lower bound */
            if( SCIPisGT(scip, lb, SCIPvarGetLbLocal(vars[v])) )
            {
               /* open a new probing node */
               if( SCIPgetProbingDepth(scip) < SCIP_MAXTREEDEPTH-10 )
               {
                  SCIP_CALL( SCIPnewProbingNode(scip) );
               }
               SCIP_CALL( chgProbingBound(scip, vars[v], solval, SCIP_BRANCHDIR_UPWARDS) );

               SCIP_CALL( SCIPpropagateProbing(scip, -1, &cutoff, &ndomreds) );

               if( cutoff || ndomreds == 0 )
               {
                  SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip)-1) );
               }
               else
               {
                  assert(SCIPvarGetProbindex(vars[v]) >= 0);
                  tightened[SCIPvarGetProbindex(vars[v])] = TRUE;
                  ++nbndtightenings;
#ifdef SCIP_MORE_DEBUG
                  SCIPdebugMsg(scip, "> tighten lower bound of variable <%s>: %g to %g (ndomreds=%lld)\n",
                        SCIPvarGetName(vars[v]), SCIPvarGetLbGlobal(vars[v]), lb, ndomreds);
#endif
               }
            }
         }
      }

      ndomredssum += ndomreds;
   }

   SCIPdebugMsg(scip, "> found %d bound tightenings and %lld induced domain reductions.\n", nbndtightenings, ndomredssum);

   return SCIP_OKAY;
}

/** main procedure of the completesol heuristic, creates and solves a sub-SCIP */
static
SCIP_RETCODE applyCompletesol(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic's private data structure */
   SCIP_RESULT*          result,             /**< result data structure */
   SCIP_Longint          nstallnodes,        /**< number of stalling nodes for the subproblem */
   SCIP_SOL*             partialsol          /**< partial solutions */
   )
{
   SCIP* subscip;
   SCIP_HASHMAP* varmapf;
   SCIP_VAR** vars;
   SCIP_VAR** subvars;
   SCIP_Bool* tightened;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_Real timelimit;
   SCIP_Real memorylimit;
   int nvars;
   int i;

   SCIP_SOL** subsols;
   int nsubsols;

   SCIP_Bool valid;
   SCIP_Bool success;
   SCIP_RETCODE retcode;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);
   assert(partialsol != NULL);

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMsg(scip, "+---+ Start Completesol heuristic +---+\n");

   /* check whether there is enough time and memory left */
   timelimit = 0.0;
   memorylimit = 0.0;
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( ! SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );

   /* substract the memory already used by the main SCIP and the estimated memory usage of external software */
   if( ! SCIPisInfinity(scip, memorylimit) )
   {
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
      memorylimit -= SCIPgetMemExternEstim(scip)/1048576.0;
   }

   /* abort if no time is left or not enough memory to create a copy of SCIP, including external memory usage */
   if( timelimit <= 0.0 || memorylimit <= 2.0*SCIPgetMemExternEstim(scip)/1048576.0 )
   {
      SCIPdebugMsg(scip, "-> not enough memory left\n");
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   /* get variable data */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* get buffer memory and initialize it to FALSE */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &tightened, nvars) );

   SCIP_CALL( SCIPstartProbing(scip) );

   SCIP_CALL( tightenVariables(scip, heurdata, vars, nvars, partialsol, tightened) );

   /* initialize the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapf, SCIPblkmem(subscip), SCIPcalcHashtableSize(5 * nvars)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

   eventhdlr = NULL;
   valid = FALSE;

   /* copy complete SCIP instance */
   SCIP_CALL( SCIPcopyConsCompression(scip, subscip, varmapf, NULL, "completesol", NULL, NULL, 0, FALSE, FALSE, TRUE, &valid) );
   SCIPdebugMsg(scip, "Copying the SCIP instance was %s complete.\n", valid ? "" : "not ");

   /* create event handler for LP events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecCompletesol, NULL) );
   if( eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for " HEUR_NAME " heuristic not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* map all variables */
   for( i = 0; i < nvars; i++ )
   {
     subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapf, vars[i]);
     assert(subvars[i] != NULL);
   }
   /* free hash map */
   SCIPhashmapFree(&varmapf);

   /* create a new problem, which fixes variables with same value in bestsol and LP relaxation */
   SCIP_CALL( createSubproblem(scip, subscip, heurdata, subvars, partialsol, tightened, &success) );
   if( !success )
   {
      SCIPdebugMsg(scip, "Error while creating completesol subproblem wrt partial solurion <%p>.\n", (void*)partialsol);
      goto TERMINATE;
   }
   SCIPdebugMsg(scip, "Completesol subproblem: %d vars, %d cons\n", SCIPgetNVars(subscip), SCIPgetNConss(subscip));

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

   /* set limits for the subproblem */
   heurdata->nodelimit = heurdata->maxnodes;
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", nstallnodes) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", heurdata->maxnodes) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", heurdata->bestsols) );

   /* forbid recursive call of heuristics and separators solving sub-SCIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* use best estimate node selection */
   if( SCIPfindNodesel(subscip, "estimate") != NULL && !SCIPisParamFixed(subscip, "nodeselection/estimate/stdpriority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/estimate/stdpriority", INT_MAX/4) );
   }

   /* use inference branching */
   if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
   }

   /* disable conflict analysis */
   if( !SCIPisParamFixed(subscip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/enable", FALSE) );
   }

#ifdef SCIP_DEBUG
   /* for debugging OFINS, enable MIP output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000000) );
#endif

   SCIP_CALL( SCIPtransformProb(subscip) );
   SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, NULL) );

   /* solve the subproblem */
   SCIPdebugMsg(scip, "solving subproblem: nstallnodes=%" SCIP_LONGINT_FORMAT ", maxnodes=%" SCIP_LONGINT_FORMAT "\n", nstallnodes, heurdata->maxnodes);
   retcode = SCIPsolve(subscip);

   SCIP_CALL( SCIPdropEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, -1) );

   /* errors in solving the subproblem should not kill the overall solving process;
    * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL( retcode );
#endif
      SCIPwarningMessage(scip, "Error while solving subproblem in completesol heuristic; sub-SCIP terminated with code <%d>\n", retcode);
   }

   /* print solving statistics of subproblem if we are in SCIP's debug mode */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );

   /* check, whether a solution was found;
    * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
    */
   nsubsols = SCIPgetNSols(subscip);
   subsols = SCIPgetSols(subscip);
   success = FALSE;
   for( i = 0; i < nsubsols && (!success || heurdata->addallsols); i++ )
   {
      SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], &success) );
      if( success )
         *result = SCIP_FOUNDSOL;
   }

   SCIPstatisticPrintf("%s statistic: fixed %6.3f integer variables, needed %6.1f seconds, %" SCIP_LONGINT_FORMAT " nodes, solution %10.4f found at node %" SCIP_LONGINT_FORMAT "\n",
      HEUR_NAME, 0.0, SCIPgetSolvingTime(subscip), SCIPgetNNodes(subscip), success ? SCIPgetPrimalbound(scip) : SCIPinfinity(scip),
      nsubsols > 0 ? SCIPsolGetNodenum(SCIPgetBestSol(subscip)) : -1 );

  TERMINATE:
   /* free subproblem */
   SCIPfreeBufferArray(scip, &subvars);
   SCIPfreeBufferArray(scip, &tightened);
   SCIP_CALL( SCIPfree(&subscip) );

   SCIP_CALL( SCIPendProbing(scip) );

   return SCIP_OKAY;
}




/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyCompletesol)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurCompletesol(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeCompletesol)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecCompletesol)
{/*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_VAR** vars;
   SCIP_SOL** partialsols;
   SCIP_Longint nstallnodes;
   int npartialsols;
   int nunknown;
   int nfracints;
   int nvars;
   int s;
   int v;

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DELAYED;

   /* do not call heuristic of node was already detected to be infeasible */
   if( nodeinfeasible )
      return SCIP_OKAY;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   *result = SCIP_DIDNOTRUN;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* do not run after restart */
   if( SCIPgetNRuns(scip) > 1 )
      return SCIP_OKAY;

   /* get variable data and return of no variables are left in the problem */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   if( nvars == 0 )
      return SCIP_OKAY;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward Completesol if it succeeded often */
   nstallnodes = (SCIP_Longint)(nstallnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nstallnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-SCIP as 100 nodes */
   nstallnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nstallnodes < heurdata->minnodes )
   {
      SCIPdebugMsg(scip, "skipping Complete: nstallnodes=%" SCIP_LONGINT_FORMAT ", minnodes=%" SCIP_LONGINT_FORMAT "\n",
         nstallnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

   /* check the number of variables with unknown value and continuous variables with fractional value */
   nunknown = 0;
   nfracints = 0;

   /* get all partial sols */
   npartialsols = SCIPgetNPartialSols(scip);
   partialsols = SCIPgetPartialSols(scip);

   /* loop over all partial solutions */
   for( s = 0; s < npartialsols; s++ )
   {
      SCIP_SOL* sol;
      SCIP_Real solval;
      SCIP_Real unknownrate;

      sol = partialsols[s];
      assert(sol != NULL);
      assert(SCIPsolIsPartial(sol));

      nunknown = 0;
      /* loop over all variables */
      for( v = 0; v < nvars; v++ )
      {
         assert(SCIPvarIsActive(vars[v]));

         /* skip continuous variables if they should ignored */
         if( SCIPvarIsIntegral(vars[v]) && heurdata->ignorecont )
            continue;

         solval = SCIPgetSolVal(scip, sol, vars[v]);

         /* we only want to count variables that are unfixed after the presolving */
         if( solval == SCIP_UNKNOWN ) /*lint !e777*/
            ++nunknown;
         else if( SCIPvarIsIntegral(vars[v]) && !SCIPisIntegral(scip, solval) )
            ++nfracints;
      }

      if( heurdata->ignorecont )
         unknownrate = nunknown/((SCIP_Real)SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) + SCIPgetNImplVars(scip));
      else
         unknownrate = nunknown/((SCIP_Real)nvars);
      SCIPdebugMsg(scip, "%d (rate %.4f) unknown solution values\n", nunknown, unknownrate);

      /* run the heuristic, if not too many unknown variables exist */
      if( unknownrate > heurdata->maxunknownrate )
         continue;

      /* all variables have a finite/known solution value and the all integer variables have an integral solution value,
       * create a new solution without solving a sub-SCIP
       */
      if( nunknown == 0 && nfracints == 0 )
      {
         SCIP_VAR** origvars;
         SCIP_SOL* newsol;
         SCIP_Bool stored;
         int norigvars;

         origvars = SCIPgetOrigVars(scip);
         norigvars = SCIPgetNOrigVars(scip);

         SCIP_CALL( SCIPcreateOrigSol(scip, &newsol, heur) );

         for( v = 0; v < norigvars; v++ )
         {
            solval = SCIPgetSolVal(scip, sol, origvars[v]);
            assert(solval != SCIP_UNKNOWN); /*lint !e777*/

            SCIP_CALL( SCIPsetSolVal(scip, newsol, origvars[v], solval) );
         }

         SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &stored) );
         if( stored )
            *result = SCIP_FOUNDSOL;
      }
      else
      {
         /* run the heuristic */
         SCIP_CALL( applyCompletesol(scip, heur, heurdata, result, nstallnodes, sol) );
      }
   }

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the completesol primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurCompletesol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create completesol primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   assert(heurdata != NULL);

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecCompletesol, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyCompletesol) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeCompletesol) );

   /* add completesol primal heuristic parameters */

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/maxunkownrate",
         "maximal rate of changed coefficients",
         &heurdata->maxunknownrate, FALSE, DEFAULT_MAXUNKRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/addallsols",
         "should all subproblem solutions be added to the original SCIP?",
         &heurdata->addallsols, TRUE, DEFAULT_ADDALLSOLS, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/lplimfac",
         "factor by which the limit on the number of LP depends on the node limit",
         &heurdata->lplimfac, TRUE, DEFAULT_LPLIMFAC, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/objweight",
         "weight of the original objective function (1: only original objective)",
         &heurdata->objweight, TRUE, DEFAULT_OBJWEIGHT, DEFAULT_MINOBJWEIGHT, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which the incumbent should be improved at least",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/ignorecont",
         "should solution values for continuous variables be ignored?",
         &heurdata->ignorecont, FALSE, DEFAULT_IGNORECONT, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/solutions",
         "heuristic stops, if the given number of improving solutions were found (-1: no limit)",
         &heurdata->bestsols, FALSE, DEFAULT_BESTSOLS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxproprounds",
         "maximal number of iterations in proagation (-1: no limit)",
         &heurdata->maxproprounds, FALSE, DEFAULT_MAXPROPROUNDS, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

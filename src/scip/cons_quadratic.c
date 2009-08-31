/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_quadratic.c,v 1.32 2009/08/31 22:32:36 bzfviger Exp $"

/**@file   cons_quadratic.c
 * @ingroup CONSHDLRS
 * @brief  constraint handler for quadratic constraints
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/** TODO list:
 * - SCIP coding style guidelines
 * - SCIP might fix variables on +/- infty; remove them in presolve and take care later
 * - constraints that are always feasible w.r.t. local/global bounds should be enabled/disabled (see logicor, setppc)
 * - round constraint bounds to integers if all coefficients and variables are (impl.) integer
 * - constraints in one variable should be replaced by linear variable or similar
 * - recognize and reformulate complementarity constraints (x*y = 0)
 * - what is a good sepalp freq? is it 1?
 * - check if some quadratic terms appear in several constraints and try to simplify (e.g., nous1)
 * - skip separation in enfolp if for current LP (check LP id) was already separated
 * - don't iterate over hash map, use array additionally
 * - watch unbounded variables to enable/disable propagation
 * - sort order in bilinvar1/bilinvar2 such that the var which is involved in more terms is in bilinvar1, and use this info propagate and AddLinearReform
 */

#include <assert.h>
#include <string.h>

#include "scip/cons_quadratic.h"
#include "scip/cons_linear.h"
#include "scip/cons_and.h"
#include "scip/cons_varbound.h"
#include "scip/intervalarith.h"
#ifdef WITH_CONSBRANCHNL
#include "cons_branchnonlinear.h"
#endif
#ifdef WITH_SOC3
#include "cons_soc3.h"
#include "cons_soc.h"
#endif
#include "scip/heur_nlp.h"
#include "scip/nlpi.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "quadratic"
#define CONSHDLR_DESC          "quadratic constraints of the form lhs <= b^T x + x^T A x <= rhs"
#define CONSHDLR_SEPAPRIORITY        10 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY       -50 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -4000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             2 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            10 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */


/*
 * Data structures
 */

/** Eventdata for variable bound change events. */
struct SCIP_EventData
{
   SCIP_CONSDATA* consdata;   /**< the constraint data */
   int            varidx;     /**< the index of a linear variable which bound change is catched, or negative for quadratic variables */
};

/** Data of a quadratic constraint. */
struct SCIP_ConsData
{
   SCIP_Real       lhs;           /**<  left hand side of constraint */
   SCIP_Real       rhs;           /**< right hand side of constraint */

   int             n_linvar;      /**< number of linear variables */
   SCIP_VAR**      linvar;        /**< linear variables */
   SCIP_Real*      lincoeff;      /**< coefficients of linear variables */

   int             n_quadvar;     /**< number of variables in quadratic terms */
   SCIP_VAR**      quadvar;       /**< variables in quadratic terms */
   SCIP_Real*      quadlincoeff;  /**< linear coefficients of quadratic variables */
   SCIP_Real*      quadsqrcoeff;  /**< coefficients of square terms of quadratic variables */
   int*            n_adjbilin;    /**< number of bilinear terms where the variable is involved */
   int**           adjbilin;      /**< indices of bilinear terms in which variable is involved */

   int             n_bilin;       /**< number of bilinear terms */
   SCIP_VAR**      bilinvar1;     /**< first variable in bilinear term */
   SCIP_VAR**      bilinvar2;     /**< second variable in bilinear term */
   SCIP_Real*      bilincoeff;    /**< coefficient of bilinear term */

   unsigned char   is_convex:1;   /**< whether quadratic function is convex */
   unsigned char   is_concave:1;  /**< whether quadratic function is concave */
   unsigned char   is_removedfixings:1; /**< whether we have removed fixed/aggr/multiaggr variables */
   unsigned char   is_propagated:1; /**< whether the constraint was propagated with respect to the current bounds */
   unsigned char   is_presolved:1;/**< whether we have checked for possibilities of upgrading or implicit integer variables */
   unsigned char   soc_added:1;   /**< whether a SOC constraint has been added for this constraint */

   SCIP_Real       lhsviol;       /**< violation of lower bound by current solution (used temporarily inside constraint handler) */
   SCIP_Real       rhsviol;       /**< violation of lower bound by current solution (used temporarily inside constraint handler) */

   SCIP_EVENTDATA* linbndchgeventdata;   /**< eventdata for bound change of linear variable */
   SCIP_EVENTDATA* quadbndchgeventdata;  /**< eventdata for bound change on quadratic variable */
   SCIP_INTERVAL*  linrange;      /**< range of each linear term */
   SCIP_INTERVAL   quadrange;     /**< range of quadratic term as used in isIntervalFeasible */
   SCIP_INTERVAL*  quadrangevar;  /**< range of quadratic term except one variable as used in propagation */
   SCIP_INTERVAL*  bilinrange;    /**< range of bilinear term as used in propagation */
};

/** quadratic term as used during presolve */
typedef struct PresolveQuadTerm
{
   SCIP_Real     lincoeff;   /**< linear coefficient of a variable */
   SCIP_Real     sqrcoeff;   /**< square coefficient of a variable */
   SCIP_HASHMAP* bilin;      /**< bilinear terms involving a variable: mapping from SCIP_VAR* to PresolveBilinItem* */
   int           component;  /**< component(=block) number in block separable form */
} PresolveQuadTerm;

/** bilinear term as used inside PresolveQuadTerm */
typedef struct PresolveBilinItem
{
   SCIP_Real     coeff;      /**< coefficient of bilinear term */
   int           bilinidx;   /**< index in bilin array */
} PresolveBilinItem;

#ifndef WITH_CONSBRANCHNL
/** Stores information about the infeasibility assigned to a variable */
struct VarInfeasibility
{
   SCIP_Real                 min;
   SCIP_Real                 max;
   SCIP_Real                 sum;
   struct VarInfeasibility*  next;
};
/** Infeasibility of a variable. */
typedef struct VarInfeasibility VARINFEASIBILITY;
#endif

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Bool       replace_sqrbinary;   /**< whether squares of binary variables are replaced by the variable itself */
   int             replace_binaryprod_length; /**< length of linear term which when multiplied with a binary variable is replaced by an auxiliary variable and an equivalent linear formulation */
   SCIP_Bool       disaggregation;      /**< whether we should disaggregate block separable quadratic constraints */
#ifdef WITH_SOC3
   SCIP_Bool       upgrade_soc;         /**< whether we should upgrade to SOC constraints, if possible */
   int             soc3_nr_auxvars;     /**< how many auxiliary variables to use when upgrading to SOC3 constraints */
#endif
   SCIP_Real       mincutefficacy;      /**< minimal efficacy of a cut in order to add it to relaxation */
   SCIP_Bool       do_scaling;          /**< whether constraints should be scaled in the feasibility check */
   SCIP_Bool       fast_propagate;      /**< whether a faster but maybe less effective propagation should be used */
   SCIP_Real       defaultbound;        /**< a bound to set for variables that are unbounded and in a nonconvex term after presolve */
   SCIP_Real       cutmaxrange;         /**< maximal range (maximal coef / minimal coef) of a cut in order to be added to LP */

   SCIP_HEUR*      nlpheur;             /**< a pointer to the NLP heuristic */
   SCIP_EVENTHDLR* eventhdlr;           /**< our handler for variable bound change events */

#ifndef WITH_CONSBRANCHNL
   SCIP_HASHMAP*      branchcand;             /**< branching candidates */
   VARINFEASIBILITY*  varinfeas;              /**< list of variable infeasibilities */
   
   char               strategy;               /**< branching strategy */
   SCIP_Real          mindistbrpointtobound;  /**< minimal (fractional) distance of branching point to bound */
#else
   SCIP_CONSHDLR*  branchnl;            /**< a pointer to the constraint handler for branching on nonlinear variables */
#endif
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

#ifndef WITH_CONSBRANCHNL
/** clears list of branching candidates */
static
SCIP_RETCODE clearBranchingCandidates(
   SCIP*            scip,         /**< SCIP data structure */
   SCIP_CONSHDLR*   conshdlr      /**< constraint handler  */
   )
{
   SCIP_CONSHDLRDATA* data;
   VARINFEASIBILITY* v;
   VARINFEASIBILITY* w = NULL;
  
   assert(scip != NULL);
   assert(conshdlr != NULL);
  
   data = SCIPconshdlrGetData(conshdlr);
   assert(data != NULL);
   assert(data->branchcand != NULL);
   
   v = data->varinfeas;
   while (v)
   {
      w = v->next;
      SCIPfreeBlockMemory(scip, &v);
      v = w;
   }
   data->varinfeas = NULL;
   
   SCIP_CALL( SCIPhashmapRemoveAll(data->branchcand) );
   
   return SCIP_OKAY;
}

/** determines branching point for a variable */
static
SCIP_RETCODE selectBranchingPoint(
   SCIP*             scip,           /**< SCIP data structure */
   SCIP_CONSHDLR*    conshdlr,       /**< constraint handler for branching on nonlinear variables */
   SCIP_VAR*         var,            /**< branching variables */
   SCIP_Real*        leftub,         /**< buffer to store new upper bound of variable in left  branch */
   SCIP_Real*        rightlb         /**< buffer to store new lower bound of variable in right branch */
   )
{
   SCIP_CONSHDLRDATA*   data;
   SCIP_Real            branchpoint;
   SCIP_Real            lb, ub;

   assert(scip != NULL);
   assert(var  != NULL);
   assert(leftub  != NULL);
   assert(rightlb != NULL);
     
   assert(scip != NULL);
   assert(conshdlr != NULL);
  
   data = SCIPconshdlrGetData(conshdlr);
   assert(data != NULL);

   branchpoint = SCIPgetVarSol(scip, var);
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);

   if (SCIPisInfinity(scip, branchpoint))
   {
      if (SCIPisPositive(scip, lb))
         branchpoint = lb + 1000;
      else
         branchpoint = 0.0;
   }
   else if (SCIPisInfinity(scip, -branchpoint))
   {
      if (SCIPisNegative(scip, ub))
         branchpoint = ub - 1000;
      else
         branchpoint = 0.0;
   }

   if (SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS)
   {
      if (!SCIPisInfinity(scip, -lb) && !SCIPisInfinity(scip, ub))
      {
         /* branch on value of LP solution
          * if it is too close to the bounds, move more into the middle of the interval */
         if (branchpoint < (1-data->mindistbrpointtobound) * lb + data->mindistbrpointtobound * ub)
            branchpoint = (1-data->mindistbrpointtobound) * lb + data->mindistbrpointtobound * ub;
         else if (branchpoint > data->mindistbrpointtobound * lb + (1-data->mindistbrpointtobound) * ub)
            branchpoint = data->mindistbrpointtobound * lb + (1-data->mindistbrpointtobound) * ub;

         /* for very tiny intervals we set it into the middle */
         if (!SCIPisGT(scip, branchpoint, lb) || !SCIPisLT(scip, branchpoint, ub))
            branchpoint = (lb+ub) * .5;
      }
      else if (!SCIPisLT(scip, lb, branchpoint))
      {
         assert(SCIPisInfinity(scip, ub));
         branchpoint = lb + MAX(0.5*ABS(lb), 1000);
      }
      else if (!SCIPisGT(scip, ub, branchpoint))
      {
         assert(SCIPisInfinity(scip, -lb));
         branchpoint = ub - MAX(0.5*ABS(ub), 1000);
      }

      *leftub = *rightlb = branchpoint;
   }
   else
   {
      if (branchpoint > ub)
         branchpoint = ub;
      else if (branchpoint < lb)
         branchpoint = lb;
      if (SCIPisIntegral(scip, branchpoint))
      {
         if (branchpoint < .5*(lb+ub))
            branchpoint += .5;
         else
            branchpoint -= .5;
      }
      *rightlb = SCIPceil(scip, branchpoint);
      *leftub  = SCIPfloor(scip, branchpoint);
   }
   
   return SCIP_OKAY;
}

/** selects a branching variable and decides about branching point */
static
SCIP_RETCODE selectBranchingVariable(
   SCIP*            scip,           /**< SCIP data structure */
   SCIP_CONSHDLR*   conshdlr,       /**< constraint handler for branching on nonlinear variables */
   SCIP_VAR**       var,            /**< buffer to store branching variable */
   SCIP_Real*       leftub,         /**< buffer to store new upper bound of branching variable in  left branch */
   SCIP_Real*       rightlb         /**< buffer to store new lower bound of branching varialbe in right branch */
   )
{
   SCIP_CONSHDLRDATA*   data;
   int                  listidx;
   SCIP_HASHMAPLIST*    candlist;
   SCIP_VAR*            cand;
   VARINFEASIBILITY*    infeas;
   SCIP_Real            candleftub, candrightlb;
   SCIP_Real            deltaminus, deltaplus;
   SCIP_Real            pscostdown, pscostup;
   SCIP_Real            score, bestscore;
  
   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(var      != NULL);
   assert(leftub   != NULL);
   assert(rightlb  != NULL);
  
   data = SCIPconshdlrGetData(conshdlr);
   assert(data != NULL);
   assert(data->branchcand != NULL);
   
   *var = NULL;
   bestscore = -1.0;
   
   for (listidx = 0; listidx < SCIPhashmapGetNLists(data->branchcand); ++listidx)
   {
      for (candlist = SCIPhashmapGetList(data->branchcand, listidx); candlist; candlist = SCIPhashmapListGetNext(candlist))
      {
         cand   = (SCIP_VAR*) SCIPhashmapListGetOrigin(candlist);
         infeas = (VARINFEASIBILITY*) SCIPhashmapListGetImage(candlist);

         switch (data->strategy)
         {
            case 'b':
               SCIP_CALL( selectBranchingPoint(scip, conshdlr, cand, &candleftub, &candrightlb) );
               assert(candleftub == candrightlb || SCIPvarGetType(cand) <= SCIP_VARTYPE_IMPLINT);
              
               if (SCIPisInfinity(scip, -SCIPvarGetLbLocal(cand)))
                  deltaminus = SCIPisInfinity(scip, infeas->max) ? SCIPinfinity(scip) : 0.1 * infeas->sum + 0.8 * infeas->min + 1.3 * infeas->max;
               else
                  deltaminus = candleftub - SCIPvarGetLbLocal(cand);
              
               if (SCIPisInfinity(scip,  SCIPvarGetUbLocal(cand)))
                  deltaplus  = SCIPisInfinity(scip, infeas->max) ? SCIPinfinity(scip) : 0.1 * infeas->sum + 0.8 * infeas->min + 1.3 * infeas->max;
               else
                  deltaplus  = SCIPvarGetUbLocal(cand) - candrightlb;
              
               break;
              
            case 'r':
               SCIP_CALL( selectBranchingPoint(scip, conshdlr, cand, &candleftub, &candrightlb) );
               assert(candleftub == candrightlb || SCIPvarGetType(cand) <= SCIP_VARTYPE_IMPLINT);
            
               if (SCIPisInfinity(scip, -SCIPvarGetLbLocal(cand)))
                  deltaplus  = SCIPisInfinity(scip, infeas->max) ? SCIPinfinity(scip) : 0.1 * infeas->sum + 0.8 * infeas->min + 1.3 * infeas->max;
               else
                  deltaplus  = candleftub - SCIPvarGetLbLocal(cand);
            
               if (SCIPisInfinity(scip,  SCIPvarGetUbLocal(cand)))
                  deltaminus = SCIPisInfinity(scip, infeas->max) ? SCIPinfinity(scip) : 0.1 * infeas->sum + 0.8 * infeas->min + 1.3 * infeas->max;
               else
                  deltaminus = SCIPvarGetUbLocal(cand) - candrightlb;
            
               break;

            case 'i':
               deltaminus = deltaplus = 0.1 * infeas->sum + 0.8 * infeas->min + 1.3 * infeas->max;
               break;

            default :
               SCIPerrorMessage("branching strategy %c unknown\n", data->strategy);
               return SCIP_ERROR;
         }

         if (SCIPisInfinity(scip, deltaminus) || SCIPisInfinity(scip, deltaplus))
            score = SCIPinfinity(scip);
         else
         {
            pscostdown = SCIPgetVarPseudocost(scip, cand, -deltaminus);
            pscostup   = SCIPgetVarPseudocost(scip, cand,  deltaplus);
            score      = SCIPgetBranchScore(scip, cand, pscostdown, pscostup);
         }
         SCIPdebugMessage("branching score variable %s = %g; \tinfeas = %g; \ttype=%d  bestscore=%g\n", SCIPvarGetName(cand), score, 0.1 * infeas->sum + 0.8 * infeas->min + 1.3 * infeas->max, SCIPvarGetType(cand), bestscore);

         if ( SCIPisSumGT(scip, score, bestscore) )
         {
            bestscore = score;
            *var      = cand;
            *leftub   = candleftub;
            *rightlb  = candrightlb;
         }
         else if ( SCIPisSumEQ(scip, score, bestscore) && !(SCIPisInfinity(scip, -SCIPvarGetLbLocal(*var)) && SCIPisInfinity(scip, SCIPvarGetUbLocal(*var))))
         { /* if best candidate so far is bounded or unbounded at atmost one side, maybe take new candidate */
            if ( (SCIPisInfinity(scip, -SCIPvarGetLbLocal(cand)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(cand))) &&
                 (SCIPisInfinity(scip, -SCIPvarGetLbLocal(*var)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(*var))) )
            { /* if both variables are unbounded but one of them is bounded on one side, take the one with the larger bound on this side (hope that this avoids branching on always the same variable) */
               if (SCIPvarGetUbLocal(cand) > SCIPvarGetUbLocal(*var) ||
                  SCIPvarGetLbLocal(cand) < SCIPvarGetLbLocal(*var))
               {
                  *var = cand;
                  *leftub   = candleftub;
                  *rightlb  = candrightlb;
               }
            }
            else if (SCIPvarGetType(*var) == SCIPvarGetType(cand))
            { /* if both have the same type, take the one with larger diameter */
               if (SCIPisLT(scip, SCIPvarGetUbLocal(*var) - SCIPvarGetLbLocal(*var), SCIPvarGetUbLocal(cand) - SCIPvarGetLbLocal(cand)))
               {
                  *var = cand;
                  *leftub   = candleftub;
                  *rightlb  = candrightlb;
               }
            }
            else if (SCIPvarGetType(*var) > SCIPvarGetType(cand))
            { /* take the one with better type ("more discrete") */
               *var = cand;
               *leftub   = candleftub;
               *rightlb  = candrightlb;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_RETCODE enforceByBranching(
	SCIP*          scip,       /**< SCIP data structure */
	SCIP_CONSHDLR* conshdlr,   /**< constraint handler */
	SCIP_RESULT*   result      /**< buffer where to store result */
	)
{
   SCIP_CONSHDLRDATA*   data;
   SCIP_VAR*            brvar = NULL;
   SCIP_Real            leftub=0.0, rightlb=0.0;
   SCIP_Real            leftobjest, rightobjest; 
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(result != NULL);
   
   data = SCIPconshdlrGetData(conshdlr);
   assert(data != NULL);
   
   if (data->varinfeas == NULL)
   { /* have no candidates for branching */
      return SCIP_OKAY;
   }

   SCIP_CALL( selectBranchingVariable(scip, conshdlr, &brvar, &leftub, &rightlb) );
   SCIP_CALL( clearBranchingCandidates(scip, conshdlr) );
   
   if (!brvar)
   {
      SCIPwarningMessage("branching variable selection failed to select a variable\n");
      return SCIP_OKAY;
   }
   
   leftobjest = SCIPcalcChildEstimate(scip, brvar, leftub);
   rightobjest = (leftub != rightlb) ? SCIPcalcChildEstimate(scip, brvar, rightlb) : leftobjest;
   if (leftobjest > SCIPinfinity(scip))
      leftobjest = SCIPinfinity(scip)/5.;
   if (rightobjest > SCIPinfinity(scip))
      rightobjest = leftobjest;

   if (SCIPvarGetStatus(brvar) == SCIP_VARSTATUS_MULTAGGR)
   {
      SCIP_NODE* node;
      SCIP_CONS* cons;
      SCIP_Real  val = 1.0;
      SCIPdebugMessage("branching on multiaggregated variable %s: new intervals: [%g, %g] [%g, %g]\n", SCIPvarGetName(brvar), SCIPvarGetLbLocal(brvar), leftub, rightlb, SCIPvarGetUbLocal(brvar));

      SCIP_CALL( SCIPcreateChild(scip, &node, 0.0, leftobjest) );
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "branch", 1, &brvar, &val, SCIPvarGetLbLocal(brvar), leftub, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      SCIP_CALL( SCIPcreateChild(scip, &node, 0.0, rightobjest) );
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "branch", 1, &brvar, &val, rightlb, SCIPvarGetUbLocal(brvar), TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }
   else
   {
      if (SCIPvarGetType(brvar) != SCIP_VARTYPE_CONTINUOUS)
      {
         SCIPdebugMessage("branching on discrete variable %s\n", SCIPvarGetName(brvar));
         SCIP_CALL( SCIPbranchVar(scip, brvar, NULL, NULL, NULL) );
      }
      else
      {
         SCIP_NODE* node;
         SCIPdebugMessage("branching on continuous variable %s: new intervals: [%g, %g] [%g, %g]\n", SCIPvarGetName(brvar), SCIPvarGetLbLocal(brvar), leftub, rightlb, SCIPvarGetUbLocal(brvar));

         SCIP_CALL( SCIPcreateChild(scip, &node, 0.0, leftobjest) );
         SCIP_CALL( SCIPchgVarUbNode(scip, node, brvar, leftub) );

         SCIP_CALL( SCIPcreateChild(scip, &node, 0.0, rightobjest) );
         SCIP_CALL( SCIPchgVarLbNode(scip, node, brvar, rightlb) );
      }
   }

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/** Updates or initializes the infeasibility of a variable.
 * If called the first time for some variable, then this variable is added to the list of branching candidates.
 */
static
SCIP_RETCODE updateVarInfeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             varinfeasibility    /**< infeasibility of variable */
   )
{
   SCIP_CONSHDLRDATA*  data;
   VARINFEASIBILITY*   varinfeas;
  
   SCIPdebugMessage("register infeasibility %g for variable %s  [%g, %g]\n", varinfeasibility, SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(varinfeasibility >= 0.0);
   assert(!SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));
  
   data = SCIPconshdlrGetData(conshdlr);
   assert(data != NULL);
   assert(data->branchcand != NULL);
   
   varinfeas = (VARINFEASIBILITY*)SCIPhashmapGetImage(data->branchcand, (void*)var);
   
   if (varinfeas == NULL)
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &varinfeas) );
      varinfeas->min  = varinfeasibility;
      varinfeas->max  = varinfeasibility;
      varinfeas->sum  = varinfeasibility;
      varinfeas->next = data->varinfeas;
      data->varinfeas = varinfeas;
      
      SCIP_CALL( SCIPhashmapInsert(data->branchcand, (void*)var, varinfeas) );
   }
   else
   {
      varinfeas->sum += varinfeasibility;
      if (varinfeasibility < varinfeas->min)
         varinfeas->min = varinfeasibility;
      if (varinfeasibility > varinfeas->max)
         varinfeas->max = varinfeasibility;
   }
   
   return SCIP_OKAY;
}
#endif

/** translate from one value of infinity to another
 * if val is >= infty1, then give infty2, else give val */
#define infty2infty(infty1, infty2, val) (val >= infty1 ? infty2 : val)

#ifndef WITH_LAPACK

/** if WITH_LAPACK not set, but WITH_IPOPT, then we can use Lapack from Ipopt and also get the Fortran naming convention from it */
#ifdef WITH_IPOPT
#define WITH_LAPACK
#include "IpoptConfig.h"
#endif

#else

/** if WITH_LAPACK is set, then also F77_FUNC should be set, otherwise we try a default that works on common systems */
#ifndef F77_FUNC
#warning "do not know about fortran naming convention for using Lapack; please consider defining F77_FUNC"
/* this is compiler and machine dependent; the following just assumes a Linux/gcc system */
#define F77_FUNC(name,NAME) name ## _
/* #define F77_FUNC_(name,NAME) name ## _ */
#endif

#endif /* ifndef/else WITH_LAPACK */

#ifdef WITH_LAPACK

/** LAPACK Fortran subroutine DSYEV */
void F77_FUNC(dsyev,DSYEV)(
   char*   jobz,    /**< 'N' to compute eigenvalues only, 'V' to compute eigenvalues and eigenvectors */
   char*   uplo,    /**< 'U' if upper triangle of A is stored, 'L' if lower triangle of A is stored */
   int*    n,       /**< dimension */
   double* A,       /**< matrix A on entry; orthonormal eigenvectors on exit, if jobz == 'V' and info == 0; if jobz == 'N', then the matrix data is destroyed */
   int*    ldA,     /**< leading dimension, probably equal to n */ 
   double* W,       /**< buffer for the eigenvalues in ascending order */
   double* WORK,    /**< workspace array */
   int*    LWORK,   /**< length of WORK; if LWORK = -1, then the optimal workspace size is calculated and returned in WORK(1) */
   int*    info     /**< == 0: successful exit; < 0: illegal argument at given position; > 0: failed to converge */
);

static
SCIP_RETCODE LapackDsyev(
   SCIP*      scip,                    /**< SCIP data structure */
   SCIP_Bool  compute_eigenvectors,    /**< whether also eigenvectors should be computed */
   int        N,                       /**< dimension */
   SCIP_Real* a,                       /**< matrix data on input (size N*N); eigenvectors on output if compute_eigenvectors == TRUE */
   SCIP_Real* w                        /**< buffer to store eigenvalues (size N) */
   )
{
   int     INFO;
   char    JOBZ = compute_eigenvectors ? 'V' : 'N';
   char    UPLO = 'L';
   int     LDA  = N;
   double* WORK = NULL;
   int     LWORK;
   double  WORK_PROBE;
   int     i;

   /* First we find out how large LWORK should be */
   LWORK = -1;
   F77_FUNC(dsyev,DSYEV)(&JOBZ, &UPLO, &N, a, &LDA, w, &WORK_PROBE, &LWORK, &INFO);
   if (INFO)
   {
      SCIPerrorMessage("There was an error when calling DSYEV. INFO = %d\n", INFO);
      return SCIP_ERROR;
   }

   LWORK = (int) WORK_PROBE;
   assert(LWORK > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &WORK, LWORK) );
   for (i = 0; i < LWORK; ++i)
      WORK[i] = i;
   F77_FUNC(dsyev,DSYEV)(&JOBZ, &UPLO, &N, a, &LDA, w, WORK, &LWORK, &INFO);
   SCIPfreeBufferArray(scip, &WORK);
   if (INFO)
   {
       SCIPerrorMessage("There was an error when calling DSYEV. INFO = %d\n", INFO);
       return SCIP_ERROR;
   }

   return SCIP_OKAY;
}
#endif

static
SCIP_DECL_EVENTEXEC(processVarEvent)
{
   SCIP_CONSDATA* consdata;
   
   assert(scip != NULL);
   assert(event != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), CONSHDLR_NAME) == 0);

   consdata = eventdata->consdata;
   assert(consdata != NULL);

   if (SCIPeventGetType(event) & SCIP_EVENTTYPE_VARFIXED)
   {
      consdata->is_removedfixings = FALSE;
      /* return SCIP_OKAY; */
   }

   if (consdata->is_propagated)
   { /* if we have a bound tightening, we might not have propagated bounds anymore */
      if (SCIPeventGetType(event) == SCIP_EVENTTYPE_UBTIGHTENED)
      {
         consdata->is_propagated =
            eventdata->varidx >= 0 &&
            (SCIPisInfinity(scip, -consdata->lhs) || consdata->lincoeff[eventdata->varidx] < 0) &&
            (SCIPisInfinity(scip,  consdata->rhs) || consdata->lincoeff[eventdata->varidx] > 0);
      }
      else if (SCIPeventGetType(event) == SCIP_EVENTTYPE_LBTIGHTENED)
      {
         consdata->is_propagated =
            eventdata->varidx >= 0 &&
            (SCIPisInfinity(scip, -consdata->lhs) || consdata->lincoeff[eventdata->varidx] > 0) &&
            (SCIPisInfinity(scip,  consdata->rhs) || consdata->lincoeff[eventdata->varidx] < 0);
      }
   }

   if (eventdata->varidx >= 0)
   { /* make linrange[varidx] invalid */
      SCIPintervalSetEmpty(SCIPinfinity(scip), &consdata->linrange[eventdata->varidx]);
   }
   else
   { /* make quadrange, quadrangevar[-varidx-1], and some bilinrange invalid */
      assert(consdata->n_quadvar > 0);
      assert(-eventdata->varidx-1 < consdata->n_quadvar);
      SCIPintervalSetEmpty(SCIPinfinity(scip), &consdata->quadrange);
      if (!SCIPintervalIsEmpty(consdata->quadrangevar[-eventdata->varidx-1]))
      {
         int i;
         SCIPintervalSetEmpty(SCIPinfinity(scip), &consdata->quadrangevar[-eventdata->varidx-1]);
         for (i = 0; i < consdata->n_adjbilin[-eventdata->varidx-1]; ++i)
            SCIPintervalSetEmpty(SCIPinfinity(scip), &consdata->bilinrange[consdata->adjbilin[-eventdata->varidx-1][i]]);
      }
   }

   return SCIP_OKAY;
}

/* TODO store index of event ? */
/** catch variable events */
static
SCIP_RETCODE catchVarEvents(
   SCIP*           scip,        /**< SCIP data structure */
   SCIP_EVENTHDLR* eventhdlr,   /**< event handler */
   SCIP_CONS*      cons         /**< constraint for which to catch bound change events */      
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->linbndchgeventdata, consdata->n_linvar) );
   for( i = 0; i < consdata->n_linvar; ++i )
   {
      SCIPintervalSetEmpty(SCIPinfinity(scip), &consdata->linrange[i]);
      consdata->linbndchgeventdata[i].consdata = consdata;
      consdata->linbndchgeventdata[i].varidx   = i;
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->linvar[i],  SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, eventhdlr, &consdata->linbndchgeventdata[i], NULL) );
   }
   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->quadbndchgeventdata, consdata->n_quadvar) );
   SCIPintervalSetEmpty(SCIPinfinity(scip), &consdata->quadrange);
   for( i = 0; i < consdata->n_quadvar; ++i )
   {
      SCIPintervalSetEmpty(SCIPinfinity(scip), &consdata->quadrangevar[i]);
      consdata->quadbndchgeventdata[i].consdata = consdata;
      consdata->quadbndchgeventdata[i].varidx   = -i-1;
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->quadvar[i], SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, eventhdlr, &consdata->quadbndchgeventdata[i], NULL) );
   }
   for( i = 0; i < consdata->n_bilin; ++i )
      SCIPintervalSetEmpty(SCIPinfinity(scip), &consdata->bilinrange[i]);
   consdata->is_propagated = FALSE;

   return SCIP_OKAY;
}

/** drop variable events */
static
SCIP_RETCODE dropVarEvents(
   SCIP*           scip,        /**< SCIP data structure */
   SCIP_EVENTHDLR* eventhdlr,   /**< event handler */
   SCIP_CONS*      cons         /**< constraint for which to catch bound change events */      
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   for (i = 0; i < consdata->n_linvar; ++i)
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->linvar[i],  SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, eventhdlr, &consdata->linbndchgeventdata[i], -1) );
   
   for (i = 0; i < consdata->n_quadvar; ++i)
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->quadvar[i], SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, eventhdlr, &consdata->quadbndchgeventdata[i], -1) );

   return SCIP_OKAY;
}

/** sets or replaces function data of constraint
 * Takes care of and release/capture of variables, but not of unlock/lock.
 */
static
SCIP_RETCODE consdataSetFunctionData(
   SCIP*          scip,      /**< SCIP data structure */
   SCIP_CONSDATA* consdata,  /**< the constraint data where to store function */
   SCIP_HASHMAP*  terms,     /**< linear and quadratic terms of function */
   int            component  /**< which component of the terms, or -1 to use all */
   )
{
   int                i, j;
   SCIP_HASHMAPLIST*  list;
   SCIP_HASHMAPLIST*  list2;
   SCIP_VAR*          var;
   SCIP_VAR*          bvar;
   PresolveQuadTerm*  term;
   PresolveBilinItem* bitem;
   int                n_linvar  = 0;
   int                n_quadvar = 0;
   int                n_bilin   = 0;
   int                i_lin   = 0;
   int                i_quad  = 0;
   int                i_bilin = 0;
   int                i_adjbilin = 0;
   
   assert(scip     != NULL);
   assert(consdata != NULL);
   assert(terms    != NULL);

   for (i = 0; i < consdata->n_linvar; ++i)
      SCIP_CALL( SCIPreleaseVar(scip, &consdata->linvar[i]) );
   for (i = 0; i < consdata->n_quadvar; ++i)
      SCIP_CALL( SCIPreleaseVar(scip, &consdata->quadvar[i]) );

   /* collect statistics */
   for (i = 0; i < SCIPhashmapGetNLists(terms); ++i)
      for (list = SCIPhashmapGetList(terms, i); list; list = SCIPhashmapListGetNext(list))
      {
         term = (PresolveQuadTerm*) SCIPhashmapListGetImage(list);
         if (component >= 0 && term->component != component)
            continue;

         if (!term->sqrcoeff && !term->bilin)
         { /* linear variable */
            if (term->lincoeff)
               ++n_linvar;
            continue;
         }
         ++n_quadvar;

         if (term->bilin)
            for (j = 0; j < SCIPhashmapGetNLists(term->bilin); ++j)
               for (list2 = SCIPhashmapGetList(term->bilin, j); list2; list2 = SCIPhashmapListGetNext(list2))
                  ++n_bilin;
      }
   n_bilin /= 2; /* because we counted each bilinear term twice */

   /* realloc memory */
   consdata->n_linvar = n_linvar;
   if (n_linvar)
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->linvar,   n_linvar) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->lincoeff, n_linvar) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->linrange, n_linvar) );
   }
   else
   {
      SCIPfreeMemoryArrayNull(scip, &consdata->linvar);
      SCIPfreeMemoryArrayNull(scip, &consdata->lincoeff);
      SCIPfreeMemoryArrayNull(scip, &consdata->linrange);
   }

   if (consdata->adjbilin)
      for (i = 0; i < consdata->n_quadvar; ++i)
         SCIPfreeMemoryArrayNull(scip, &consdata->adjbilin[i]);

   consdata->n_quadvar = n_quadvar;
   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->quadvar,      n_quadvar) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->quadlincoeff, n_quadvar) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->quadsqrcoeff, n_quadvar) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->n_adjbilin,   n_quadvar) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->adjbilin,     n_quadvar) );

   consdata->n_bilin = n_bilin;
   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->bilinvar1,  n_bilin) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->bilinvar2,  n_bilin) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->bilincoeff, n_bilin) );

   /* set constraint data */
   for (i = 0; i < SCIPhashmapGetNLists(terms); ++i)
      for (list = SCIPhashmapGetList(terms, i); list; list = SCIPhashmapListGetNext(list))
      {
         term = (PresolveQuadTerm*) SCIPhashmapListGetImage(list);
         if (component >= 0 && term->component != component)
            continue;
         var  = (SCIP_VAR*) SCIPhashmapListGetOrigin(list);

         if (!term->sqrcoeff && !term->bilin)
         { /* linear variable */
            if (term->lincoeff)
            {
               assert(i_lin < n_linvar);
               assert(consdata->linvar   != NULL);
               assert(consdata->lincoeff != NULL);
               consdata->linvar  [i_lin] = var;
               consdata->lincoeff[i_lin] = term->lincoeff;
               ++i_lin;
            }
            continue;
         }

         /* quadratic variable */
         assert(i_quad < n_quadvar);
         consdata->quadvar[i_quad]      = var;
         consdata->quadlincoeff[i_quad] = term->lincoeff;
         consdata->quadsqrcoeff[i_quad] = term->sqrcoeff;
         consdata->n_adjbilin[i_quad]   = 0;
         if (term->bilin)
            for (j = 0; j < SCIPhashmapGetNLists(term->bilin); ++j)
               for (list2 = SCIPhashmapGetList(term->bilin, j); list2; list2 = SCIPhashmapListGetNext(list2))
                  ++consdata->n_adjbilin[i_quad];

         if (consdata->n_adjbilin[i_quad])
         { /* bilinear terms of quadratic variable */
            SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->adjbilin[i_quad], consdata->n_adjbilin[i_quad]) );
            i_adjbilin = 0;
            for (j = 0; j < SCIPhashmapGetNLists(term->bilin); ++j)
               for (list2 = SCIPhashmapGetList(term->bilin, j); list2; list2 = SCIPhashmapListGetNext(list2))
               {
                  bvar  = (SCIP_VAR*)          SCIPhashmapListGetOrigin(list2);
                  bitem = (PresolveBilinItem*) SCIPhashmapListGetImage(list2);

                  if (bitem->bilinidx >= 0)
                  { /* bilinear term has been created already, just store its index in adjbilin */
                     consdata->adjbilin[i_quad][i_adjbilin] = bitem->bilinidx;
                     ++i_adjbilin;
                  }
                  else
                  { /* bilinear term needs to be created here */
                     PresolveQuadTerm*  bterm;
                     PresolveBilinItem* bbitem;
                     bterm  = (PresolveQuadTerm*)  SCIPhashmapGetImage(terms, bvar);
                     assert(bterm != NULL);
                     assert(bterm->bilin != NULL);
                     bbitem = (PresolveBilinItem*) SCIPhashmapGetImage(bterm->bilin, var);
                     assert(bbitem != NULL);
                     assert(bbitem->bilinidx < 0);
                     assert(bbitem->coeff == bitem->coeff);

                     assert(i_bilin < n_bilin);
                     consdata->bilinvar1 [i_bilin] = var;
                     consdata->bilinvar2 [i_bilin] = bvar;
                     consdata->bilincoeff[i_bilin] = bitem->coeff;
                     consdata->adjbilin[i_quad][i_adjbilin] = i_bilin;

                     bitem ->bilinidx = i_bilin;
                     bbitem->bilinidx = i_bilin;

                     ++i_adjbilin;
                     ++i_bilin;
                  }
               }
         }
         else
         {
            consdata->adjbilin[i_quad]   = NULL;
         }

         ++i_quad;
      }

   assert(i_lin   == n_linvar);
   assert(i_quad  == n_quadvar);
   assert(i_bilin == n_bilin);

   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->quadrangevar, n_quadvar) );
   if (n_bilin)
      SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->bilinrange, n_bilin) );
   else
      SCIPfreeMemoryArrayNull(scip, &consdata->bilinrange);

   for (i = 0; i < consdata->n_linvar; ++i)
      SCIP_CALL( SCIPcaptureVar(scip, consdata->linvar[i]) );
   for (i = 0; i < consdata->n_quadvar; ++i)
      SCIP_CALL( SCIPcaptureVar(scip, consdata->quadvar[i]) );

   return SCIP_OKAY;
}

static
void consdataFree(SCIP* scip, SCIP_CONSDATA** consdata)
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   SCIPfreeMemoryArrayNull(scip, &(*consdata)->linvar);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->lincoeff);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->linrange);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->linbndchgeventdata);

   SCIPfreeMemoryArrayNull(scip, &(*consdata)->quadvar);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->quadlincoeff);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->quadsqrcoeff);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->quadrangevar);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->quadbndchgeventdata);
   if ((*consdata)->adjbilin)
   {
      int i;
      assert((*consdata)->n_adjbilin);
      for (i = 0; i < (*consdata)->n_quadvar; ++i)
         SCIPfreeMemoryArrayNull(scip, &(*consdata)->adjbilin[i]);
      SCIPfreeMemoryArray(scip, &(*consdata)->adjbilin);
      SCIPfreeMemoryArray(scip, &(*consdata)->n_adjbilin);
   }

   SCIPfreeMemoryArrayNull(scip, &(*consdata)->bilinvar1);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->bilinvar2);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->bilincoeff);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->bilinrange);

   SCIPfreeMemory(scip, consdata);
   *consdata = NULL;
}

/** Adds a variable with coefficients to the presolve data structure. */
static
SCIP_RETCODE presolveQuadTermAdd(
   SCIP*         scip,        /**< SCIP data structure */
   SCIP_HASHMAP* terms,       /**< the terms where to add something */
   SCIP_VAR*     var,         /**< the variable to add */
   SCIP_Real     lincoeff,    /**< the linear coefficient of the variable */
   SCIP_Real     sqrcoeff,    /**< the coefficient of the square term of the variable */
   SCIP_VAR*     bilinvar,    /**< the other variable in a bilinear term, or NULL */  
   SCIP_Real     bilincoeff   /**< the coefficient of the bilinear term */
   )
{
   PresolveQuadTerm* term;
   
   assert(scip  != NULL);
   assert(terms != NULL);
   assert(var   != NULL);
   
   term = (PresolveQuadTerm*) SCIPhashmapGetImage(terms, var);
   if( !term )
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &term) );
      term->lincoeff  = lincoeff;
      term->sqrcoeff  = sqrcoeff;
      term->bilin     = NULL;
      term->component = -1;
      SCIPhashmapInsert( terms, var, term );
   }
   else
   {
      term->lincoeff += lincoeff;
      term->sqrcoeff += sqrcoeff;
   }

   if( bilinvar )
   {
      PresolveBilinItem* bitem;
      if( term->bilin == NULL )
      {
         SCIP_CALL( SCIPhashmapCreate(&term->bilin, SCIPblkmem(scip), 10) ); /* TODO 10 a good value? */
         SCIP_CALL( SCIPallocBlockMemory(scip, &bitem) );
         SCIP_CALL( SCIPhashmapInsert(term->bilin, bilinvar, bitem) );
         bitem->coeff    = bilincoeff;
         bitem->bilinidx = -1;
      }
      else
      {
         bitem = (PresolveBilinItem*) SCIPhashmapGetImage(term->bilin, bilinvar);
         if( !bitem )
         {
            SCIP_CALL( SCIPallocBlockMemory(scip, &bitem) );
            SCIP_CALL( SCIPhashmapInsert(term->bilin, bilinvar, bitem) );
            bitem->coeff    = bilincoeff;
            bitem->bilinidx = -1;
         }
         else
         {
            bitem->coeff   += bilincoeff;
         }
      }
   }

   return SCIP_OKAY;
}

/** frees the presolve data */
static
void presolveQuadTermFree(
   SCIP*            scip,   /**< SCIP data structure */
   SCIP_HASHMAP**   terms   /**< quadratic terms to free */
   )
{
   int i, j;
   SCIP_HASHMAPLIST*  list1;
   SCIP_HASHMAPLIST*  list2;
   PresolveQuadTerm*  term;
   PresolveBilinItem* bitem;
   
   assert(scip  != NULL);
   assert(terms != NULL);

   if (!*terms)
      return;

   for (i = 0; i < SCIPhashmapGetNLists(*terms); ++i)
   {
      list1 = SCIPhashmapGetList(*terms, i);
      for (; list1; list1 = SCIPhashmapListGetNext(list1) )
      {
         term = (PresolveQuadTerm*)SCIPhashmapListGetImage(list1);
         if (term->bilin)
         {
            for (j = 0; j < SCIPhashmapGetNLists(term->bilin); ++j)
            {
               list2 = SCIPhashmapGetList(term->bilin, j);
               for (; list2; list2 = SCIPhashmapListGetNext(list2) )
               {
                  bitem = (PresolveBilinItem*)SCIPhashmapListGetImage(list2);
                  SCIPfreeBlockMemory(scip, &bitem);
               }
            }
            SCIPhashmapFree(&term->bilin);
         }
         SCIPfreeBlockMemory(scip, &term);
      }
   }
   SCIPhashmapFree(terms);
   assert(*terms == NULL);
}

/* Prints the quadratic function stored in the presolve data structure. */
#ifdef SCIP_DEBUG
static
void presolveQuadTermPrint(
   SCIP*         scip,   /**< SCIP data structure */
   SCIP_HASHMAP* terms,  /**< quadratic terms */
   FILE*         file    /**< file to print to, or NULL for stdout */
   )
{
   SCIP_HASHMAPLIST*  list1;
   SCIP_HASHMAPLIST*  list2;
   SCIP_VAR*          var;
   SCIP_VAR*          bvar;
   PresolveQuadTerm*  term;
   PresolveBilinItem* bitem;
   int                i, j;
   
   assert(scip  != NULL);
   assert(terms != NULL);

   for (i = 0; i < SCIPhashmapGetNLists(terms); ++i)
   {
      list1 = SCIPhashmapGetList(terms, i);
      for (; list1; list1 = SCIPhashmapListGetNext(list1) )
      {
         var  = (SCIP_VAR*)        SCIPhashmapListGetOrigin(list1);
         term = (PresolveQuadTerm*)SCIPhashmapListGetImage (list1);
         if (term->lincoeff)
            SCIPinfoMessage(scip, file, "%+g*%s ",      term->lincoeff, SCIPvarGetName(var));
         if (term->sqrcoeff)
            SCIPinfoMessage(scip, file, "%+g*sqr(%s) ", term->sqrcoeff, SCIPvarGetName(var));
         if (term->bilin)
         {
            SCIPinfoMessage(scip, file, "+%s*(", SCIPvarGetName(var));
            for (j = 0; j < SCIPhashmapGetNLists(term->bilin); ++j)
            {
               list2 = SCIPhashmapGetList(term->bilin, j);
               for (; list2; list2 = SCIPhashmapListGetNext(list2) )
               {
                  bvar  = (SCIP_VAR*)         SCIPhashmapListGetOrigin(list2);
                  bitem = (PresolveBilinItem*)SCIPhashmapListGetImage (list2);
                  if (var < bvar)
                     SCIPinfoMessage(scip, file, "%+g*%s ", bitem->coeff, SCIPvarGetName(bvar));
               }
            }
            SCIPinfoMessage(scip, file, ")");
         }
      }
   }
}
#endif

/** Adds a linear term into the presolve data structure. */
static
SCIP_RETCODE presolveAddLinearTerm(
   SCIP*          scip,     /**< SCIP data structure */
   SCIP_HASHMAP*  terms,    /**< presolve quadratic terms */
   SCIP_Real*     constant, /**< presolve constant term */
   SCIP_Real      coeff,    /**< coefficient of variable */
   SCIP_VAR*      var       /**< variable to add */
   )
{
   assert(scip  != NULL);
   assert(terms != NULL);
   assert(constant != NULL);
   assert(var   != NULL);
   assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL);

   if (!coeff)
      return SCIP_OKAY;

   switch (SCIPvarGetStatus(var))
   {
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
      {
         if (SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var))) /* var is fixed */
            *constant += coeff * SCIPvarGetLbGlobal(var);
         else
            SCIP_CALL( presolveQuadTermAdd(scip, terms, var, coeff, 0.0, NULL, 0.0) );
         break;
      }
      case SCIP_VARSTATUS_FIXED:
      {
         assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)));
         *constant += coeff * SCIPvarGetLbGlobal(var);
         break;
      }
      case SCIP_VARSTATUS_AGGREGATED:
      { /* var is replaced by scalar * aggrvar + constant */
         *constant += coeff * SCIPvarGetAggrConstant(var);
         SCIP_CALL( presolveAddLinearTerm(scip, terms, constant, coeff * SCIPvarGetAggrScalar(var), SCIPvarGetAggrVar(var)) );
         break;
      }
      case SCIP_VARSTATUS_MULTAGGR:
      { /* var is replaced by sum_i scalar_i * aggrvar_i + constant */
         int        navars   = SCIPvarGetMultaggrNVars(var);
         SCIP_VAR** avars    = SCIPvarGetMultaggrVars(var);
         SCIP_Real* ascalars = SCIPvarGetMultaggrScalars(var);
         int        j        = 0;
         SCIPdebugMessage("replace multiaggregated variable %s\n", SCIPvarGetName(var));
         *constant += coeff * SCIPvarGetMultaggrConstant(var);
         for (; j < navars; ++j)
            SCIP_CALL( presolveAddLinearTerm(scip, terms, constant, coeff * ascalars[j], avars[j]) );
         break;
      }
      case SCIP_VARSTATUS_NEGATED:
      { /* var is replaced by constant - negvar */
         *constant += coeff * SCIPvarGetNegationConstant(var);
         SCIP_CALL( presolveAddLinearTerm(scip, terms, constant, -coeff, SCIPvarGetNegationVar(var)) );
         break;
      }
      default:
         SCIPerrorMessage("unexpected variable status: %d\n", SCIPvarGetStatus(var));
         return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE presolveAddBilinearTerm(
   SCIP*          scip,     /**< SCIP data structure */
   SCIP_HASHMAP*  terms,    /**< presolve quadratic terms */
   SCIP_Real*     constant, /**< presolve constant term */
   SCIP_Real      coeff,    /**< coefficient of bilinear variable */
   SCIP_VAR*      var1,     /**< first  variable in bilinear term to add */
   SCIP_VAR*      var2      /**< second variable in bilinear term to add */
   )
{
   assert(scip  != NULL);
   assert(terms != NULL);
   assert(constant != NULL);
   assert(var1  != NULL);
   assert(var2  != NULL);
   assert(SCIPvarGetStatus(var1) != SCIP_VARSTATUS_ORIGINAL);
   assert(SCIPvarGetStatus(var2) != SCIP_VARSTATUS_ORIGINAL);

   if (!coeff)
      return SCIP_OKAY;

   switch (SCIPvarGetStatus(var2))
   {
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
      {
         if (SCIPisEQ(scip, SCIPvarGetLbGlobal(var2), SCIPvarGetUbGlobal(var2)))
         { /* var2 is fixed */
            SCIP_CALL( presolveAddLinearTerm(scip, terms, constant, coeff*SCIPvarGetLbGlobal(var2), var1) );
         }
         else if (var1 == var2)
         { /* var2 is not fixed but equal to var1 */
            SCIP_CALL( presolveQuadTermAdd(scip, terms, var1, 0.0, coeff, NULL, 0.0) );
         }
         else
         { /* var2 is not fixed */
            assert(var1 != var2);
            if (SCIPvarIsActive(var1))
            {
               SCIP_CALL( presolveQuadTermAdd(scip, terms, var1, 0.0, 0.0, var2, coeff) );
               SCIP_CALL( presolveQuadTermAdd(scip, terms, var2, 0.0, 0.0, var1, coeff) );
            }
            else
            {
               SCIP_CALL( presolveAddBilinearTerm(scip, terms, constant, coeff, var2, var1) );
            }
         }
         break;
      }
      case SCIP_VARSTATUS_FIXED:
      {
         assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(var2), SCIPvarGetUbGlobal(var2)));
         SCIP_CALL( presolveAddLinearTerm(scip, terms, constant, coeff*SCIPvarGetLbGlobal(var2), var1) );
         break;
      }
      case SCIP_VARSTATUS_AGGREGATED:
      { /* var2 is replaced by scalar * aggrvar + constant */
         SCIP_CALL( presolveAddLinearTerm  (scip, terms, constant, coeff*SCIPvarGetAggrConstant(var2), var1) );
         SCIP_CALL( presolveAddBilinearTerm(scip, terms, constant, coeff*SCIPvarGetAggrScalar(var2),   var1, SCIPvarGetAggrVar(var2)) );
         break;
      }
      case SCIP_VARSTATUS_MULTAGGR:
      { /* var2 is replaced by sum_i scalar_i * aggrvar_i + constant */
         int        navars   = SCIPvarGetMultaggrNVars(var2);
         SCIP_VAR** avars    = SCIPvarGetMultaggrVars(var2);
         SCIP_Real* ascalars = SCIPvarGetMultaggrScalars(var2);
         int        j        = 0;
         SCIPdebugMessage("replace multiaggregated variable %s in term %g*%s*%s by %g", SCIPvarGetName(var2), coeff, SCIPvarGetName(var1), SCIPvarGetName(var2), SCIPvarGetMultaggrConstant(var2));
         for (; j < navars; ++j)
         {
            SCIPdebugPrintf(" + %g%s", ascalars[j], SCIPvarGetName(avars[j]));
            SCIP_CALL( presolveAddBilinearTerm(scip, terms, constant, coeff * ascalars[j], var1, avars[j]) );
         }
         SCIPdebugPrintf("\n");
         SCIP_CALL( presolveAddLinearTerm(scip, terms, constant, coeff * SCIPvarGetMultaggrConstant(var2), var1) );
         break;
      }
      case SCIP_VARSTATUS_NEGATED:
      { /* var2 is replaced by constant - negvar */
         SCIP_CALL( presolveAddLinearTerm  (scip, terms, constant, coeff * SCIPvarGetNegationConstant(var2), var1) );
         SCIP_CALL( presolveAddBilinearTerm(scip, terms, constant, -coeff, var1, SCIPvarGetNegationVar(var2)) );
         break;
      }
      default:
         SCIPerrorMessage("unexpected variable status: %d\n", SCIPvarGetStatus(var2));
         return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE presolveAddSquareTerm(
   SCIP*          scip,     /**< SCIP data structure */
   SCIP_HASHMAP*  terms,    /**< presolve quadratic terms */
   SCIP_Real*     constant, /**< presolve constant term */
   SCIP_Real      coeff,    /**< coefficient of variable */
   SCIP_VAR*      var       /**< variable to add */
   )
{
   assert(scip  != NULL);
   assert(terms != NULL);
   assert(constant != NULL);
   assert(var   != NULL);
   assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL);

   if (!coeff)
      return SCIP_OKAY;

   switch (SCIPvarGetStatus(var))
   {
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
      {
         if (SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)))
         { /* variable is fixed now */
            SCIP_Real val = SCIPvarGetLbGlobal(var);
            *constant += coeff * val * val;
         }
         else
         { /* variable is not fixed */
            SCIP_CALL( presolveQuadTermAdd(scip, terms, var, 0.0, coeff, NULL, 0.0) );
         }
         break;
      }
      case SCIP_VARSTATUS_FIXED:
      {
         SCIP_Real val = SCIPvarGetLbGlobal(var);
         assert(SCIPisEQ(scip, val, SCIPvarGetUbGlobal(var)));
         *constant += coeff * val * val;
         break;
      }
      case SCIP_VARSTATUS_AGGREGATED:
      { /* var is replaced by scalar * aggrvar + constant */
         SCIP_Real aggconst  = SCIPvarGetAggrConstant(var);
         SCIP_Real aggscalar = SCIPvarGetAggrScalar(var);
         SCIP_CALL( presolveAddSquareTerm(scip, terms, constant,   coeff*aggscalar*aggscalar, SCIPvarGetAggrVar(var)) );
         SCIP_CALL( presolveAddLinearTerm(scip, terms, constant, 2*coeff*aggscalar*aggconst , SCIPvarGetAggrVar(var)) );
         *constant += coeff * aggconst * aggconst;
         break;
      }
      case SCIP_VARSTATUS_MULTAGGR:
      { /* var is replaced by sum_i scalar_i * aggrvar_i + constant */
         int        navars    = SCIPvarGetMultaggrNVars(var);
         SCIP_VAR** avars     = SCIPvarGetMultaggrVars(var);
         SCIP_Real* ascalars  = SCIPvarGetMultaggrScalars(var);
         SCIP_Real  aconstant = SCIPvarGetMultaggrConstant(var);
         int        jj        = 0;
         int        k;
         SCIPdebugMessage("replace multiaggregated variable %s by %g", SCIPvarGetName(var), aconstant);
         *constant += coeff * aconstant * aconstant;
         for (; jj < navars; ++jj)
         {
            SCIPdebugPrintf("+ %g%s", ascalars[jj], SCIPvarGetName(avars[jj]));
            SCIP_CALL( presolveAddLinearTerm(scip, terms, constant, 2*coeff*ascalars[jj]*aconstant, avars[jj]) );
            for (k = 0; k < jj; ++k)
               SCIP_CALL( presolveAddBilinearTerm(scip, terms, constant, 2*coeff*ascalars[jj]*ascalars[k], avars[k], avars[jj]) );
            SCIP_CALL( presolveAddSquareTerm(scip, terms, constant, coeff*ascalars[jj]*ascalars[jj], avars[jj]) );
         }
         SCIPdebugPrintf("\n");
         break;
      }
      case SCIP_VARSTATUS_NEGATED:
      { /* var is replaced by constant - negvar */
         SCIP_Real negconst = SCIPvarGetNegationConstant(var);
         *constant += coeff * negconst * negconst;
         SCIP_CALL( presolveAddLinearTerm(scip, terms, constant, -2*coeff*negconst, SCIPvarGetNegationVar(var)) );
         SCIP_CALL( presolveAddSquareTerm(scip, terms, constant, coeff, SCIPvarGetNegationVar(var)) );
         break;
      }
      default:
         SCIPerrorMessage("unexpected variable status: %d\n", SCIPvarGetStatus(var));
         return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** gets constraint function data as set of PresolveQuadTerm's for reformulation and preprocessing
 * replaces fixed and aggregated variables
 * replaces squares of binary variables
 */
static
SCIP_RETCODE presolveCreateQuadTerm(
   SCIP*           scip,              /**< SCIP data structure */
   SCIP_HASHMAP**  terms,             /**< storage for PresolveQuadTerm's */
   SCIP_Real*      constant,          /**< storage for constant */
   SCIP_Bool*      have_change,       /**< whether a change in the function has been identified (variables fixed, aggregated,...) */
   SCIP_CONS*      cons,              /**< constraint */
   SCIP_Bool       replace_sqrbinary  /**< whether to replace squares of binary variables */
   )
{
   SCIP_CONSDATA*     consdata;
   SCIP_VAR*          var;
   int                i, j;
   SCIP_Real          coeff;
   SCIP_HASHMAPLIST*  listitem;
   SCIP_HASHMAPLIST*  blistitem;
   PresolveQuadTerm*  term;
   PresolveBilinItem* bilinitem;
   SCIP_Bool          have_bilin;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(terms != NULL);
   assert(constant != NULL);
   assert(have_change != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( SCIPhashmapCreate(terms, SCIPblkmem(scip), consdata->n_quadvar + consdata->n_linvar) );
   *constant    = 0.;
   *have_change = FALSE;

   for (i = 0; i < consdata->n_linvar; ++i)
   {
      var = consdata->linvar[i];
      *have_change |= (SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN) && (SCIPvarGetStatus(var) != SCIP_VARSTATUS_LOOSE);
      *have_change |= SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
      SCIP_CALL( presolveAddLinearTerm(scip, *terms, constant, consdata->lincoeff[i], var) );
   }

   for (j = 0; j < consdata->n_quadvar; ++j)
   {
      if (!consdata->quadlincoeff[j] && !consdata->quadsqrcoeff[j])
         continue;

      var = consdata->quadvar[j];
      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL); /* should not happen after transformation */
      switch (SCIPvarGetStatus(var))
      {
         case SCIP_VARSTATUS_LOOSE:
         case SCIP_VARSTATUS_COLUMN:
         {
            if (SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)))
            { /* variable is fixed now */
               SCIP_Real val = SCIPvarGetLbGlobal(var);
               *constant += consdata->quadlincoeff[j] * val;
               *constant += consdata->quadsqrcoeff[j] * val*val;
               *have_change = TRUE;
            }
            else
            { /* variable is not fixed */
               SCIP_CALL( presolveQuadTermAdd(scip, *terms, var, consdata->quadlincoeff[j], consdata->quadsqrcoeff[j], NULL, 0.0) );
            }
            break;
         }
         case SCIP_VARSTATUS_FIXED:
         {
            SCIP_Real val = SCIPvarGetLbGlobal(var);
            assert(SCIPisEQ(scip, val, SCIPvarGetUbGlobal(var)));
            *constant += consdata->quadlincoeff[j] * val;
            *constant += consdata->quadsqrcoeff[j] * val*val;
            *have_change = TRUE;
            break;
         }
         case SCIP_VARSTATUS_AGGREGATED:
         case SCIP_VARSTATUS_MULTAGGR:
         case SCIP_VARSTATUS_NEGATED:
         {
            SCIP_CALL( presolveAddLinearTerm(scip, *terms, constant, consdata->quadlincoeff[j], var) );
            SCIP_CALL( presolveAddSquareTerm(scip, *terms, constant, consdata->quadsqrcoeff[j], var) );
            *have_change = TRUE;
            break;
         }
         default:
            SCIPerrorMessage("unexpected variable status: %d\n", SCIPvarGetStatus(var));
            return SCIP_ERROR;
      }
   }

   for (j = 0; j < consdata->n_bilin; ++j)
   { /* go through bilinear terms */
      var   = consdata->bilinvar1[j];
      coeff = consdata->bilincoeff[j];

      if (!coeff)
         continue;

      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL);
      switch (SCIPvarGetStatus(var))
      {
         case SCIP_VARSTATUS_LOOSE:
         case SCIP_VARSTATUS_COLUMN:
         {
            if (SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)))
            { /* variable is now fixed */
               SCIP_CALL( presolveAddLinearTerm(scip, *terms, constant, coeff*SCIPvarGetLbGlobal(var), consdata->bilinvar2[j]) );
               *have_change = TRUE;
            }
            else
            { /* variable is not fixed */
               *have_change |= (SCIPvarGetStatus(consdata->bilinvar2[j]) != SCIP_VARSTATUS_COLUMN) && (SCIPvarGetStatus(consdata->bilinvar2[j]) != SCIP_VARSTATUS_LOOSE);
               *have_change |= SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->bilinvar2[j]), SCIPvarGetUbGlobal(consdata->bilinvar2[j]));
               SCIP_CALL( presolveAddBilinearTerm(scip, *terms, constant, coeff, var, consdata->bilinvar2[j]) );
            }
            break;
         }
         case SCIP_VARSTATUS_FIXED:
         case SCIP_VARSTATUS_AGGREGATED:
         case SCIP_VARSTATUS_MULTAGGR:
         case SCIP_VARSTATUS_NEGATED:
         {
            SCIP_CALL( presolveAddBilinearTerm(scip, *terms, constant, consdata->bilincoeff[j], consdata->bilinvar2[j], var) );
            *have_change = TRUE;
            break;
         }
         default:
            SCIPerrorMessage("unexpected variable status: %d\n", SCIPvarGetStatus(var));
            return SCIP_ERROR;
      }
   }

   /* clean 0's */
   for (i = 0; i < SCIPhashmapGetNLists(*terms); ++i)
   {
      listitem = SCIPhashmapGetList(*terms, i);
      while ( listitem )
      {
         var = (SCIP_VAR*) SCIPhashmapListGetOrigin(listitem);
         term = (PresolveQuadTerm*) SCIPhashmapListGetImage(listitem);
         if (replace_sqrbinary && (SCIPvarGetType(var) == SCIP_VARTYPE_BINARY) && term->sqrcoeff)
         { /* replace square of binary variable by variable itself (x^2 = x) */
            term->lincoeff += term->sqrcoeff;
            term->sqrcoeff  = 0.0;
            *have_change = TRUE;
         }

         if (SCIPisZero(scip, term->lincoeff))
            term->lincoeff = 0.0;
         if (SCIPisZero(scip, term->sqrcoeff))
            term->sqrcoeff = 0.0;

         have_bilin = FALSE;
         if (term->bilin)
         {
            for (j = 0; j < SCIPhashmapGetNLists(term->bilin); ++j)
            {
               blistitem = SCIPhashmapGetList(term->bilin, j);
               while ( blistitem )
               {
                  bilinitem = (PresolveBilinItem*) SCIPhashmapListGetImage(blistitem);
                  if( SCIPisZero(scip, bilinitem->coeff) )
                  {
                     SCIPfreeBlockMemory(scip, &bilinitem);
                     SCIP_CALL( SCIPhashmapRemove(term->bilin, SCIPhashmapListGetOrigin(blistitem)) );
                     blistitem = SCIPhashmapGetList(term->bilin, j); /* restart from the beginning of this list again */
                  }
                  else
                  {
                     blistitem = SCIPhashmapListGetNext(blistitem);
                     have_bilin = TRUE;
                  }
               }
            }
            if (!have_bilin)
               SCIPhashmapFree(&term->bilin);
         }

         if (!term->lincoeff && !term->sqrcoeff && !have_bilin)
         {
            SCIPfreeBlockMemory(scip, &term);
            SCIP_CALL( SCIPhashmapRemove(*terms, SCIPhashmapListGetOrigin(listitem)) );
            listitem = SCIPhashmapGetList(*terms, i); /* restart from the beginning of this list again */
         }
         else
         {
            listitem = SCIPhashmapListGetNext(listitem);
         }           
      }
   }

   return SCIP_OKAY;
}

/** assigns component (block) number to all variables associated to a given one */
static
void presolveAssignComponent(
   SCIP_HASHMAP* terms,
   SCIP_VAR*     var,
   int           component
   )
{
   PresolveQuadTerm* term;
   SCIP_HASHMAPLIST* list;
   int               j;
   
   assert(terms != NULL);
   assert(var   != NULL);
   
   term = (PresolveQuadTerm*)SCIPhashmapGetImage(terms, var);
   assert(term != NULL);
   assert(term->component == component || term->component == -1);

   if (term->component >= 0)
      return;
   
   term->component = component;

   if (!term->bilin)
      return;
   for (j = 0; j < SCIPhashmapGetNLists(term->bilin); ++j)
      for (list = SCIPhashmapGetList(term->bilin, j); list; list = SCIPhashmapListGetNext(list))
         presolveAssignComponent(terms, (SCIP_VAR*)SCIPhashmapListGetOrigin(list), component);
}

/** find components (blocks) of connected quadratic variables in constraint function data
 */
static
SCIP_RETCODE presolveFindComponents(
   SCIP*           scip,        /**< SCIP data structure */
   SCIP_HASHMAP*   terms,       /**< constraint function in form of PresolveQuadTerm's */
   int*            n_components /**< buffer where to store number of components found */
   )
{
   SCIP_HASHMAPLIST* list;
   SCIP_HASHMAPLIST* list2;
   PresolveQuadTerm* term;
   int               i, j;
   
   assert(scip != NULL);
   assert(terms != NULL);
   assert(n_components != NULL);

   *n_components = 0;
   for (i = 0; i < SCIPhashmapGetNLists(terms); ++i)
      for (list = SCIPhashmapGetList(terms, i); list; list = SCIPhashmapListGetNext(list))
      {
         term = (PresolveQuadTerm*) SCIPhashmapListGetImage(list);
         if (term->sqrcoeff == 0.0 && term->bilin == NULL) /* linear variable */
            continue;
         if (term->component >= 0) /* component already assigned */
            continue;
         
         term->component = *n_components;
         ++*n_components;
         
         if (!term->bilin)
            continue;
         for (j = 0; j < SCIPhashmapGetNLists(term->bilin); ++j)
            for (list2 = SCIPhashmapGetList(term->bilin, j); list2; list2 = SCIPhashmapListGetNext(list2))
               presolveAssignComponent(terms, (SCIP_VAR*)SCIPhashmapListGetOrigin(list2), term->component);
      }
   
   return SCIP_OKAY;
}

/** disaggregates a block separable constraint into several quadratic constraints
 */
static
SCIP_RETCODE presolveDisaggregate(
   SCIP*           scip,        /**< SCIP data structure */
   SCIP_CONSHDLR*  conshdlr,    /**< the constraint handler itself */
   SCIP_CONS*      cons,        /**< constraint */
   SCIP_HASHMAP*   terms,       /**< constraint function in form of PresolveQuadTerm's */
   SCIP_Real       constant,    /**< constant part of constraint function */
   int             n_components /**< number of components (blocks) */
)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_HASHMAPLIST*  list;
   PresolveQuadTerm*  term;
   SCIP_VAR**         newlinvar;
   SCIP_Real*         newlincoeff;
   SCIP_CONS*         blockcons;
   SCIP_CONSDATA*     blockconsdata;
   SCIP_VAR*          auxvar;
   SCIP_VAR*          var;
   int                i, k;
   int                n_newlinvar;
   int                i_newlinvar;
   char               name[255];

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(terms != NULL);
   assert(n_components > 1);
   assert(SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("disaggregate constraint %s\n", SCIPconsGetName(cons));

   n_newlinvar = n_components - 1;
   for (i = 0; i < SCIPhashmapGetNLists(terms); ++i)
      for (list = SCIPhashmapGetList(terms, i); list; list = SCIPhashmapListGetNext(list))
      {
         term = (PresolveQuadTerm*) SCIPhashmapListGetImage(list);
         if (term->component == -1)
         {
            assert(SCIPisZero(scip, term->sqrcoeff));
            assert(term->bilin == NULL);
            if (SCIPisZero(scip, term->lincoeff))
               continue;
            ++n_newlinvar;
         }
      }

   SCIP_CALL( SCIPallocMemoryArray(scip, &newlinvar,   n_newlinvar) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &newlincoeff, n_newlinvar) );

   i_newlinvar = 0;
   for (k = 1; k < n_components; ++k)
   {
      SCIP_CALL( SCIPallocMemory( scip, &blockconsdata) );

      /* we need to enforce only one bound here
       * we cannot make auxvar implicit integer then
       * hope it does not have bad effect on bound tightening
       */
      blockconsdata->lhs          = SCIPisInfinity(scip, -consdata->lhs) ? -SCIPinfinity(scip) : 0.;
      blockconsdata->rhs          = SCIPisInfinity(scip,  consdata->rhs) ?  SCIPinfinity(scip) : 0.;

      blockconsdata->n_linvar     = 0;
      blockconsdata->linvar       = NULL;
      blockconsdata->lincoeff     = NULL;

      blockconsdata->n_quadvar    = 0;
      blockconsdata->quadvar      = NULL;
      blockconsdata->quadlincoeff = NULL;
      blockconsdata->quadsqrcoeff = NULL;
      blockconsdata->n_adjbilin   = NULL;
      blockconsdata->adjbilin     = NULL;

      blockconsdata->n_bilin      = 0;
      blockconsdata->bilincoeff   = NULL;
      blockconsdata->bilinvar1    = NULL;
      blockconsdata->bilinvar2    = NULL;

      blockconsdata->linrange     = NULL;
      blockconsdata->quadrangevar = NULL;
      blockconsdata->linbndchgeventdata = NULL;
      blockconsdata->quadbndchgeventdata = NULL;
      blockconsdata->bilinrange   = NULL;

      SCIP_CALL( consdataSetFunctionData(scip, blockconsdata, terms, k) );
      assert(!blockconsdata->n_linvar);
      assert(!blockconsdata->linvar);
      assert(!blockconsdata->lincoeff);
      assert(!blockconsdata->linrange);
      /* we scale the new constraints by n_components, so that the sum of the feasibility violations stays below the feasibility tolerance in a solution */
#if 0
      for (i = 0; i < blockconsdata->n_quadvar; ++i)
      {
         blockconsdata->quadlincoeff[i] *= n_components;
         blockconsdata->quadsqrcoeff[i] *= n_components;
      }
      for (i = 0; i < blockconsdata->n_bilin; ++i)
         blockconsdata->bilincoeff[i] *= n_components;
#endif

      SCIPsnprintf(name, 255, "%s#%u", SCIPconsGetName(cons), k);
      SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0., SCIP_VARTYPE_CONTINUOUS, SCIPconsIsInitial(cons), SCIPconsIsRemovable(cons), NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, auxvar) );

      SCIP_CALL( SCIPallocMemoryArray(scip, &blockconsdata->linvar,   1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &blockconsdata->lincoeff, 1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &blockconsdata->linrange, 1) );

      blockconsdata->n_linvar    = 1;
      blockconsdata->linvar[0]   = auxvar;
      blockconsdata->lincoeff[0] = -1;

      blockconsdata->is_convex     = FALSE;
      blockconsdata->is_concave    = FALSE;
      blockconsdata->is_removedfixings = TRUE;
      blockconsdata->is_propagated = FALSE;
      blockconsdata->is_presolved  = FALSE; /* so that in the next presolve round maybe auxvar is made implicit integer */
      blockconsdata->soc_added     = FALSE;

      SCIP_CALL( SCIPcreateCons(scip, &blockcons, name, conshdlr, blockconsdata, 
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
         SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
         SCIPconsIsStickingAtNode(cons)) );

      SCIP_CALL( SCIPaddCons(scip, blockcons) );
      SCIPdebugMessage("created new constraint %s: ", SCIPconsGetName(blockcons));
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintCons(scip, blockcons, NULL) );
#endif
      if (SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVING)
         SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, blockcons) );

      SCIP_CALL( SCIPreleaseCons(scip, &blockcons) );

      assert(i_newlinvar < n_newlinvar);
      newlinvar  [i_newlinvar] = auxvar;
#if 0
      newlincoeff[i_newlinvar] = 1./n_components;
#else
      newlincoeff[i_newlinvar] = 1.;
#endif
      ++i_newlinvar;
      SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   }

   SCIP_CALL( consdataSetFunctionData(scip, consdata, terms, 0) );
   assert(!consdata->n_linvar);
   assert(!consdata->linvar);
   assert(!consdata->lincoeff);
   assert(!consdata->linrange);

   for (i = 0; i < SCIPhashmapGetNLists(terms); ++i)
      for (list = SCIPhashmapGetList(terms, i); list; list = SCIPhashmapListGetNext(list))
      {
         term = (PresolveQuadTerm*) SCIPhashmapListGetImage(list);
         if (term->component == -1 && !SCIPisZero(scip, term->lincoeff))
         {
            var = (SCIP_VAR*) SCIPhashmapListGetOrigin(list);
            assert(i_newlinvar < n_newlinvar);
            newlinvar  [i_newlinvar] = var;
            newlincoeff[i_newlinvar] = term->lincoeff;
            ++i_newlinvar;
            SCIP_CALL( SCIPcaptureVar(scip, var) );
         }
      }
   assert(i_newlinvar == n_newlinvar);
   consdata->n_linvar = n_newlinvar;
   consdata->linvar   = newlinvar;
   consdata->lincoeff = newlincoeff;
   if (!SCIPisInfinity(scip, -consdata->lhs))
      consdata->lhs -= constant;
   if (!SCIPisInfinity(scip,  consdata->rhs))
      consdata->rhs -= constant;
   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->linrange, consdata->n_linvar) );

   SCIPdebugMessage("modified constraint %s to: ", SCIPconsGetName(cons));
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
#endif

   return SCIP_OKAY;
}
#if 0
/** Reformulates products of binary variables as AND constraint.
 * For a product y*x, with x and y binary variables, the product is replaced by a new auxiliary variable z and the constraint z = {x and y} is added.
 */
static
SCIP_RETCODE presolveTryAddAND(
   SCIP*           scip,        /**< SCIP data structure */
   SCIP_CONS*      cons,        /**< constraint */
   SCIP_HASHMAP*   terms,       /**< constraint function in form of PresolveQuadTerm's */
   int*            n_consadded  /**< buffer where to add the number of AND constraints added */
)
{
   SCIP_VAR*          x;
   SCIP_VAR*          y;
   PresolveQuadTerm*  xterm;
   PresolveQuadTerm*  yterm;
   PresolveBilinItem* bitem;
   char               name[255];
   int                i, j;
   SCIP_HASHMAPLIST*  list;
   SCIP_HASHMAPLIST*  blist;
   SCIP_VAR*          auxvar;
   SCIP_CONS*         andcons;
   SCIP_VAR*          vars[2];

   assert(scip != NULL);
   assert(cons != NULL);
   assert(terms != NULL);
   assert(n_consadded != NULL);
   
   *n_consadded = 0;
   
   for (i = 0; i < SCIPhashmapGetNLists(terms); ++i)
      for (list = SCIPhashmapGetList(terms, i); list; list = SCIPhashmapListGetNext(list))
      {
         x = (SCIP_VAR*) SCIPhashmapListGetOrigin(list);
         if (SCIPvarGetType(x) != SCIP_VARTYPE_BINARY)
            continue;

         xterm = (PresolveQuadTerm*) SCIPhashmapListGetImage(list);
         if (!xterm->bilin)
            continue;

         for (j = 0; j < SCIPhashmapGetNLists(xterm->bilin); ++j)
         {
            blist = SCIPhashmapGetList(xterm->bilin, j);
            while (blist)
            {
               y = (SCIP_VAR*) SCIPhashmapListGetOrigin(blist);
               if (SCIPvarGetType(y) != SCIP_VARTYPE_BINARY)
               {
                  blist = SCIPhashmapListGetNext(blist);
                  continue;
               }

               SCIPsnprintf(name, 255, "prod%s*%s", SCIPvarGetName(x), SCIPvarGetName(y));
               SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, 0., 1., 0., SCIP_VARTYPE_BINARY, TRUE, TRUE, NULL, NULL, NULL, NULL) );
               SCIP_CALL( SCIPaddVar(scip, auxvar) );

               vars[0] = x;
               vars[1] = y;
               SCIPsnprintf(name, 255, "%sAND%s", SCIPvarGetName(x), SCIPvarGetName(y));
               SCIP_CALL( SCIPcreateConsAnd(scip, &andcons, name, auxvar, 2, vars,
                  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                  SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                  SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),   SCIPconsIsStickingAtNode(cons)) );
               SCIP_CALL( SCIPaddCons(scip, andcons) );
               SCIPdebugMessage("added AND constraint: ");
#ifdef SCIP_DEBUG
               SCIPprintCons(scip, andcons, NULL);
#endif
               SCIP_CALL( SCIPreleaseCons(scip, &andcons) );

               bitem = (PresolveBilinItem*) SCIPhashmapListGetImage(blist);
               SCIP_CALL( presolveQuadTermAdd(scip, terms, auxvar, bitem->coeff, 0.0, NULL, 0.0) );
               SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

               /* delete BilinItem and entry from xterm->bilin */
               SCIPfreeBlockMemory(scip, &bitem);
               SCIP_CALL( SCIPhashmapRemove(xterm->bilin, y) );
               blist = SCIPhashmapGetList(xterm->bilin, j);

               /* find and delete BilinItem in yterm and entry from yterm->bilin */
               yterm = (PresolveQuadTerm*) SCIPhashmapGetImage(terms, y);
               assert(yterm != NULL);
               assert(yterm->bilin != NULL);
               bitem = (PresolveBilinItem*) SCIPhashmapGetImage(yterm->bilin, x);
               assert(bitem != NULL);
               SCIPfreeBlockMemory(scip, &bitem);
               SCIP_CALL( SCIPhashmapRemove(yterm->bilin, x) );

               ++*n_consadded;
            }
         }
      }

   if (*n_consadded)
   {
      /* clean empty bilin's and empty quad terms */
      for (i = 0; i < SCIPhashmapGetNLists(terms); ++i)
      {
         list = SCIPhashmapGetList(terms, i);
         while (list)
         {
            yterm = (PresolveQuadTerm*) SCIPhashmapListGetImage(list);
            if (yterm->bilin && SCIPhashmapIsEmpty(yterm->bilin))
               SCIPhashmapFree(&yterm->bilin);

            if (!yterm->bilin && !yterm->lincoeff && !yterm->sqrcoeff)
            {
               SCIPfreeBlockMemory(scip, &yterm);
               SCIP_CALL( SCIPhashmapRemove(terms, (SCIP_VAR*)SCIPhashmapListGetOrigin(list)) );
               list = SCIPhashmapGetList(terms, i); /* restart from beginning of list */
            }
            else
               list = SCIPhashmapListGetNext(list);
         }
      }
   }

   return SCIP_OKAY;
}
#endif

/** Reformulates products of binary times bounded continuous variables as system of linear inequalities (plus auxiliary variable).
 * A product x*y, with y a binary variable and x a continous variable with finite bounds,
 * an auxiliary variable z and the inequalities x^L * y <= z <= x^U * y and x - (1-y)*x^L <= z <= x - (1-y)*x^U are added.
 */
static
SCIP_RETCODE presolveTryAddLinearReform(
   SCIP*           scip,        /**< SCIP data structure */
   SCIP_CONS*      cons,        /**< constraint */
   SCIP_HASHMAP*   terms,       /**< constraint function in form of PresolveQuadTerm's */
   int*            n_varsadded, /**< buffer where to store the number of auxiliary variables added */
   int             maxnrvar     /**< maximal number of variables in linear term to consider when replacing by one auxiliary variable */
)
{
   SCIP_VAR**         xvars  = NULL;
   SCIP_Real*         xcoeff = NULL;
   SCIP_INTERVAL      xbnds;
   SCIP_INTERVAL      tmp;
   int                n_xvars;
   SCIP_VAR*          y;
   SCIP_VAR*          bvar;
   SCIP_HASHMAPLIST*  list;
   SCIP_HASHMAPLIST*  blist;
   PresolveQuadTerm*  yterm;
   PresolveBilinItem* bitem;
   PresolveQuadTerm*  bterm;
   char               name[255];
   int                i, j;
   int                n_bilin;
   SCIP_VAR*          auxvar;
   SCIP_CONS*         auxcons;
   SCIP_Bool          maxnrvar_full = FALSE; /* indicates whether we stopped collecting xvars because the maxnrvar limit was reached */

   assert(scip != NULL);
   assert(cons != NULL);
   assert(terms != NULL);
   assert(n_varsadded != NULL);
   
   *n_varsadded = 0;

   if (!maxnrvar)
      return SCIP_OKAY;
   
   for (i = 0; i < SCIPhashmapGetNLists(terms); i = maxnrvar_full ? i : i+1)
      for (list = SCIPhashmapGetList(terms, i); list; list = maxnrvar_full ? list : SCIPhashmapListGetNext(list))
      {
         maxnrvar_full = FALSE;
         
         y = (SCIP_VAR*) SCIPhashmapListGetOrigin(list);
         if (SCIPvarGetType(y) != SCIP_VARTYPE_BINARY)
            continue;
         
         yterm = (PresolveQuadTerm*) SCIPhashmapListGetImage(list);
         if (!yterm->bilin)
            continue;
         
         n_bilin = SCIPhashmapGetNEntries(yterm->bilin);
         if (!n_bilin)
            continue;
         
         n_xvars = 0;
         SCIPintervalSet(&xbnds, 0.);
         SCIP_CALL( SCIPreallocBufferArray(scip, &xvars,  MIN(maxnrvar, n_bilin)+2) ); /* add 2 for later use when creating linear constraints */
         SCIP_CALL( SCIPreallocBufferArray(scip, &xcoeff, MIN(maxnrvar, n_bilin)+2) );
         
         /* setup list of variables x_i with coefficients a_i that are multiplied with binary y: y*(sum_i a_i*x_i)
          * and compute range of sum_i a_i*x_i
          */
         for (j = 0; j < SCIPhashmapGetNLists(yterm->bilin) && !maxnrvar_full; ++j)
            for (blist = SCIPhashmapGetList(yterm->bilin, j); blist; blist = SCIPhashmapListGetNext(blist))
            {
               if (n_xvars >= maxnrvar)
               {
                  maxnrvar_full = TRUE;
                  break;
               }
               
               bvar  = (SCIP_VAR*)          SCIPhashmapListGetOrigin(blist);
               if (SCIPisInfinity(scip, -SCIPvarGetLbGlobal(bvar)) || SCIPisInfinity(scip, SCIPvarGetUbGlobal(bvar)))
                  continue;
               
               bitem = (PresolveBilinItem*) SCIPhashmapListGetImage(blist);
               xvars [n_xvars] = bvar;
               xcoeff[n_xvars] = bitem->coeff;
               
               SCIPintervalSetBounds(&tmp, MIN(SCIPvarGetLbGlobal(bvar), SCIPvarGetUbGlobal(bvar)), MAX(SCIPvarGetLbGlobal(bvar), SCIPvarGetUbGlobal(bvar)));
               SCIPintervalMulScalar(SCIPinfinity(scip), &tmp, tmp, bitem->coeff);
               SCIPintervalAdd(SCIPinfinity(scip), &xbnds, xbnds, tmp);

               ++n_xvars;
               
               /* free PresolveBilinItem structure in bilin map at y; entry from bilin map is removed later */
               SCIPfreeBlockMemory(scip, &bitem);

               /* free PresolveBilinItem structure in bilin map at bvar and delete entry from this bilin map */
               bterm = (PresolveQuadTerm*)SCIPhashmapGetImage(terms, bvar);
               assert(bterm != NULL);
               bitem = (PresolveBilinItem*) SCIPhashmapGetImage(bterm->bilin, y);
               assert(bitem != NULL);
               SCIPfreeBlockMemory(scip, &bitem);
               SCIP_CALL( SCIPhashmapRemove(bterm->bilin, y) );
            }
         
         if (!n_xvars) /* all x_i seem to be unbounded */
            continue;
         
         /* remove entries from bilin map at y */
         for (j = 0; j < n_xvars; ++j)
            SCIP_CALL( SCIPhashmapRemove(yterm->bilin, xvars[j]) );
         
         assert(!SCIPisInfinity(scip, -SCIPintervalGetInf(xbnds)));
         assert(!SCIPisInfinity(scip,  SCIPintervalGetSup(xbnds)));
         
         if (n_xvars == 1 && SCIPvarGetType(xvars[0]) == SCIP_VARTYPE_BINARY)
         { /* product of two binary variables, replace by auxvar and AND constraint */
            /* add auxiliary variable z */
            SCIPsnprintf(name, 255, "prod%s*%s", SCIPvarGetName(y), SCIPvarGetName(xvars[0]));
            SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, 0., 1., 0., SCIP_VARTYPE_BINARY, TRUE, TRUE, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, auxvar) );
            
            xvars[1] = y;
            SCIPsnprintf(name, 255, "%sAND%s", SCIPvarGetName(y), SCIPvarGetName(xvars[0]));
            SCIP_CALL( SCIPcreateConsAnd(scip, &auxcons, name, auxvar, 2, xvars,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
               SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),   SCIPconsIsStickingAtNode(cons)) );
            SCIP_CALL( SCIPaddCons(scip, auxcons) );
            SCIPdebugMessage("added AND constraint: ");
#ifdef SCIP_DEBUG
            SCIPprintCons(scip, auxcons, NULL);
#endif
            SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
            
            SCIP_CALL( presolveQuadTermAdd(scip, terms, auxvar, xcoeff[0], 0.0, NULL, 0.0) );

            SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );
         }
         else
         { /* product of binary avariable with more than one binary or with continuous variables, replace by auxvar and linear constraints */
            /* add auxiliary variable z */
            if (n_xvars == 1)
               SCIPsnprintf(name, 255, "prod%s*%s", SCIPvarGetName(y), SCIPvarGetName(xvars[0]));
            else
               SCIPsnprintf(name, 255, "prod%s*%s*more", SCIPvarGetName(y), SCIPvarGetName(xvars[0]));
            SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, MIN(0., SCIPintervalGetInf(xbnds)), MAX(0., SCIPintervalGetSup(xbnds)), 0., SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, auxvar) );

            if (!SCIPisZero(scip, SCIPintervalGetInf(xbnds)))
            { /* add 0 <= z - xbnds.inf * y constraint (as varbound constraint) */
               SCIPsnprintf(name, 255, "linreform%s_1", SCIPvarGetName(y));
               SCIP_CALL( SCIPcreateConsVarbound(scip, &auxcons, name, auxvar, y, -SCIPintervalGetInf(xbnds), 0, SCIPinfinity(scip),
                  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                  SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                  SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),   SCIPconsIsStickingAtNode(cons)) );
               SCIP_CALL( SCIPaddCons(scip, auxcons) );
               SCIPdebugMessage("added varbound constraint: ");
#ifdef SCIP_DEBUG
               SCIPprintCons(scip, auxcons, NULL);
#endif
               SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
            }
            if (!SCIPisZero(scip, SCIPintervalGetSup(xbnds)))
            { /* add z - xbnds.sup * y <= 0 constraint (as varbound constraint) */
               SCIPsnprintf(name, 255, "linreform%s_2", SCIPvarGetName(y));
               SCIP_CALL( SCIPcreateConsVarbound(scip, &auxcons, name, auxvar, y, -SCIPintervalGetSup(xbnds), -SCIPinfinity(scip), 0,
                  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                  SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                  SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),   SCIPconsIsStickingAtNode(cons)) );
               SCIP_CALL( SCIPaddCons(scip, auxcons) );
               SCIPdebugMessage("added varbound constraint: ");
#ifdef SCIP_DEBUG
               SCIPprintCons(scip, auxcons, NULL);
#endif
               SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
            }

            /* add xbnds.inf <= sum_i a_i*x_i + xbnds.inf * y - z constraint */
            xvars[n_xvars]    = y;
            xvars[n_xvars+1]  = auxvar;
            xcoeff[n_xvars]   = SCIPintervalGetInf(xbnds);
            xcoeff[n_xvars+1] = -1;

            SCIPsnprintf(name, 255, "linreform%s_3", SCIPvarGetName(y));
            SCIP_CALL( SCIPcreateConsLinear(scip, &auxcons, name, n_xvars+2, xvars, xcoeff, SCIPintervalGetInf(xbnds), SCIPinfinity(scip),
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
               SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),   SCIPconsIsStickingAtNode(cons)) );
            SCIP_CALL( SCIPaddCons(scip, auxcons) );
            SCIPdebugMessage("added linear constraint: ");
#ifdef SCIP_DEBUG
            SCIPprintCons(scip, auxcons, NULL);
#endif
            SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );

            /* add sum_i a_i*x_i + xbnds.sup * y - z <= xbnds.sup constraint */
            xcoeff[n_xvars]   = SCIPintervalGetSup(xbnds);

            SCIPsnprintf(name, 255, "linreform%s_4", SCIPvarGetName(y));
            SCIP_CALL( SCIPcreateConsLinear(scip, &auxcons, name, n_xvars+2, xvars, xcoeff, -SCIPinfinity(scip), SCIPintervalGetSup(xbnds),
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
               SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),   SCIPconsIsStickingAtNode(cons)) );
            SCIP_CALL( SCIPaddCons(scip, auxcons) );
            SCIPdebugMessage("added linear constraint: ");
#ifdef SCIP_DEBUG
            SCIPprintCons(scip, auxcons, NULL);
#endif
            SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );

            /* add z to this constraint */
            SCIP_CALL( presolveQuadTermAdd(scip, terms, auxvar, 1.0, 0.0, NULL, 0.0) );
            
            SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );
         }
         
         ++*n_varsadded;
      }
      
   if (*n_varsadded)
   {
      printf("added %d sets of constraints to reformulate product with binary variable in constraint %s\n", *n_varsadded, SCIPconsGetName(cons));

      /* clean empty bilin's and empty quad terms */
      for (i = 0; i < SCIPhashmapGetNLists(terms); ++i)
      {
         list = SCIPhashmapGetList(terms, i);
         while (list)
         {
            yterm = (PresolveQuadTerm*) SCIPhashmapListGetImage(list);
            if (yterm->bilin && SCIPhashmapIsEmpty(yterm->bilin))
               SCIPhashmapFree(&yterm->bilin);
            
            if (!yterm->bilin && !yterm->lincoeff && !yterm->sqrcoeff)
            {
               SCIPfreeBlockMemory(scip, &yterm);
               SCIP_CALL( SCIPhashmapRemove(terms, (SCIP_VAR*)SCIPhashmapListGetOrigin(list)) );
               list = SCIPhashmapGetList(terms, i); /* restart from beginning of list */
            }
            else
               list = SCIPhashmapListGetNext(list);
         }
      }
   }

   SCIPfreeBufferArrayNull(scip, &xvars);
   SCIPfreeBufferArrayNull(scip, &xcoeff);

   return SCIP_OKAY;
}

#ifdef WITH_SOC3
/** creates SOC3 constraints for input sum_{i=1}^n alpha_i x_i^2 <= beta y^2.
 * if n>2, calls same function recursively
 */
static
SCIP_RETCODE presolveCreateSOC3(
   SCIP*        scip,       /**< SCIP data structure */
   SCIP_VAR**   lhsvar,     /**< variables on left hand side (x_i) */
   SCIP_Real*   lhscoeff,   /**< coefficients of variables on left hand side (alpha_i) */
   int          nlhs,       /**< number of variables on left hand side (n) */
   SCIP_VAR*    rhsvar,     /**< variable on right hand side (y) */
   SCIP_Real    rhscoeff,   /**< coefficient of variable on right hand side (beta) */
   const char*  basename,   /**< prefix for variable and constraint name */
   SCIP_CONS*   origcons,   /**< original constraint for which this SOC3 set is added */
   int          soc3_nr_auxvars /**< number of auxiliary variables to use for a SOC3 constraint, or 0 if automatic */
   )
{
   char       name[255];
   SCIP_CONS* newcons;
   SCIP_VAR*  auxvar1;
   SCIP_VAR*  auxvar2;

   assert(scip     != NULL);
   assert(lhsvar   != NULL);
   assert(lhscoeff != NULL);
   assert(nlhs     >= 2);
   assert(rhsvar   != NULL);
   assert(basename != NULL);
   
   if (nlhs == 2)
   { /* end of recursion */
      assert(lhscoeff[0] > 0);
      assert(lhscoeff[1] > 0);
      assert(rhscoeff    > 0);
      assert(lhsvar[0] != NULL);
      assert(lhsvar[1] != NULL);
      SCIPsnprintf(name, 255, "%s_soc3", basename); /* TODO think of better name */
      if (soc3_nr_auxvars)
      {
         SCIP_CALL( SCIPcreateConsSOC3(scip, &newcons, name, soc3_nr_auxvars,
            sqrt(lhscoeff[0]), lhsvar[0], sqrt(lhscoeff[1]), lhsvar[1], sqrt(rhscoeff), rhsvar,
            SCIPconsIsInitial(origcons), SCIPconsIsSeparated(origcons), SCIPconsIsEnforced(origcons),
            SCIPconsIsChecked(origcons), SCIPconsIsPropagated(origcons),  SCIPconsIsLocal(origcons),
            SCIPconsIsDynamic(origcons), SCIPconsIsRemovable(origcons), SCIPconsIsStickingAtNode(origcons)) );
      }
      else
      {
         SCIP_CALL( SCIPcreateConsSOC3byTol(scip, &newcons, name, SCIPfeastol(scip),
            sqrt(lhscoeff[0]), lhsvar[0], sqrt(lhscoeff[1]), lhsvar[1], sqrt(rhscoeff), rhsvar,
            SCIPconsIsInitial(origcons), SCIPconsIsSeparated(origcons), SCIPconsIsEnforced(origcons),
            SCIPconsIsChecked(origcons), SCIPconsIsPropagated(origcons),  SCIPconsIsLocal(origcons),
            SCIPconsIsDynamic(origcons), SCIPconsIsRemovable(origcons), SCIPconsIsStickingAtNode(origcons)) );
      }
      SCIP_CALL( SCIPaddCons(scip, newcons) );
      SCIPdebugMessage("added SOC3 constraint %s: ", name);
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintCons(scip, newcons, NULL) );
#endif
      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
      
      return SCIP_OKAY;
   }

   if (nlhs == 3)
   { /* a bit special case too */
      /* for first two variables on lhs, create a new aux.var and a new SOC3 */
      SCIPsnprintf(name, 255, "%s#z1", basename);
      SCIP_CALL( SCIPcreateVar(scip, &auxvar1, name, 0., SCIPinfinity(scip), 0., SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, auxvar1) );

      /* constraint alpha_0 x_0^2 + alpha_1 x_1^2 <= auxvar^2 */
      SCIP_CALL( presolveCreateSOC3(scip, lhsvar, lhscoeff, 2, auxvar1, 1., name, origcons, soc3_nr_auxvars) );

      /* create new constraint alpha_2 x_2^2 + auxvar^2 <= rhscoeff * rhsvar^2 */
      SCIPsnprintf(name, 255, "%s_soc3", basename); /* TODO think of better name */
      if (soc3_nr_auxvars)
      {
         SCIP_CALL( SCIPcreateConsSOC3(scip, &newcons, name, soc3_nr_auxvars,
            sqrt(lhscoeff[2]), lhsvar[2], 1., auxvar1, sqrt(rhscoeff), rhsvar,
            SCIPconsIsInitial(origcons), SCIPconsIsSeparated(origcons), SCIPconsIsEnforced(origcons),
            SCIPconsIsChecked(origcons), SCIPconsIsPropagated(origcons),  SCIPconsIsLocal(origcons),
            SCIPconsIsDynamic(origcons), SCIPconsIsRemovable(origcons), SCIPconsIsStickingAtNode(origcons)) );
      }
      else
      {
         SCIP_CALL( SCIPcreateConsSOC3byTol(scip, &newcons, name, SCIPfeastol(scip),
            sqrt(lhscoeff[2]), lhsvar[2], 1., auxvar1, sqrt(rhscoeff), rhsvar,
            SCIPconsIsInitial(origcons), SCIPconsIsSeparated(origcons), SCIPconsIsEnforced(origcons),
            SCIPconsIsChecked(origcons), SCIPconsIsPropagated(origcons),  SCIPconsIsLocal(origcons),
            SCIPconsIsDynamic(origcons), SCIPconsIsRemovable(origcons), SCIPconsIsStickingAtNode(origcons)) );
      }
      SCIP_CALL( SCIPaddCons(scip, newcons) );
      SCIPdebugMessage("added SOC3 constraint %s: ", name);
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintCons(scip, newcons, NULL) );
#endif
      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
            
      SCIP_CALL( SCIPreleaseVar(scip, &auxvar1) );
      
      return SCIP_OKAY;
   }
   
   /* nlhs >= 4 */
   
   SCIPsnprintf(name, 255, "%s#z1", basename);
   SCIP_CALL( SCIPcreateVar(scip, &auxvar1, name, 0., SCIPinfinity(scip), 0., SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, auxvar1) );

   /* constraints for left half of lhs */
   SCIP_CALL( presolveCreateSOC3(scip, lhsvar, lhscoeff, nlhs/2, auxvar1, 1., name, origcons, soc3_nr_auxvars) );

   SCIPsnprintf(name, 255, "%s#z2", basename);
   SCIP_CALL( SCIPcreateVar(scip, &auxvar2, name, 0., SCIPinfinity(scip), 0., SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, auxvar2) );

   SCIP_CALL( presolveCreateSOC3(scip, &lhsvar[nlhs/2], &lhscoeff[nlhs/2], nlhs-nlhs/2, auxvar2, 1., name, origcons, soc3_nr_auxvars) );
   
   /* SOC constraint binding both auxvar's */
   SCIPsnprintf(name, 255, "%s_soc3", basename); /* TODO think of better name */
   if (soc3_nr_auxvars)
   {
      SCIP_CALL( SCIPcreateConsSOC3(scip, &newcons, name, soc3_nr_auxvars,
         1., auxvar1, 1., auxvar2, sqrt(rhscoeff), rhsvar,
         SCIPconsIsInitial(origcons), SCIPconsIsSeparated(origcons), SCIPconsIsEnforced(origcons),
         SCIPconsIsChecked(origcons), SCIPconsIsPropagated(origcons),  SCIPconsIsLocal(origcons),
         SCIPconsIsDynamic(origcons), SCIPconsIsRemovable(origcons), SCIPconsIsStickingAtNode(origcons)) );
   }
   else
   {
      SCIP_CALL( SCIPcreateConsSOC3byTol(scip, &newcons, name, SCIPfeastol(scip),
         1., auxvar1, 1., auxvar2, sqrt(rhscoeff), rhsvar,
         SCIPconsIsInitial(origcons), SCIPconsIsSeparated(origcons), SCIPconsIsEnforced(origcons),
         SCIPconsIsChecked(origcons), SCIPconsIsPropagated(origcons),  SCIPconsIsLocal(origcons),
         SCIPconsIsDynamic(origcons), SCIPconsIsRemovable(origcons), SCIPconsIsStickingAtNode(origcons)) );
   }
   SCIP_CALL( SCIPaddCons(scip, newcons) );
   SCIPdebugMessage("added SOC3 constraint %s: ", name);
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintCons(scip, newcons, NULL) );
#endif
   SCIP_CALL( SCIPreleaseCons(scip, &newcons) );

   SCIP_CALL( SCIPreleaseVar(scip, &auxvar1) );
   SCIP_CALL( SCIPreleaseVar(scip, &auxvar2) );
   
   return SCIP_OKAY;
}

/** adds for a constraint sum_i alpha_i x_i^2 <= beta * y a set of SOC3 constraints */
static
SCIP_RETCODE presolveTryAddSOC(
   SCIP*           scip,        /**< SCIP data structure */
   SCIP_CONS*      cons,        /**< constraint */
   SCIP_HASHMAP*   terms,       /**< constraint function in form of PresolveQuadTerm's */
   SCIP_Real       constant,    /**< constant part of constraint function */
   int             soc3_nr_auxvars /**< number of auxiliary variables to use for a SOC3 constraint, or 0 if automatic */
)
{
   SCIP_CONSDATA* consdata;
   int            nterms;
   SCIP_VAR**     lhsvar;
   SCIP_Real*     lhscoeff;
   int            lhscount = 0;
   SCIP_VAR*      rhsvar = NULL; 
   SCIP_Real      rhscoeff = 0.0;
   SCIP_HASHMAPLIST* list;
   SCIP_VAR*         var;
   PresolveQuadTerm* term;
   int            i;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(terms != NULL);
   assert(!SCIPisInfinity(scip, ABS(constant)));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if (!SCIPisZero(scip, consdata->rhs - constant) && !SCIPisZero(scip, consdata->lhs - constant))
      return SCIP_OKAY; /* do not even need to try */

   nterms = SCIPhashmapGetNEntries(terms);

   SCIP_CALL( SCIPallocBufferArray(scip, &lhsvar,   nterms-1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lhscoeff, nterms-1) );

   if (SCIPisZero(scip, consdata->rhs - constant))
   { /* try whether constraint on upper bound is SOC */
      for (i = 0; i < SCIPhashmapGetNLists(terms); ++i)
         for (list = SCIPhashmapGetList(terms, i); list; list = SCIPhashmapListGetNext(list))
         {
            var  = (SCIP_VAR*)         SCIPhashmapListGetOrigin(list);
            term = (PresolveQuadTerm*) SCIPhashmapListGetImage(list);
            if (term->lincoeff || term->bilin)
            { /* variable with linear part or variable in bilinear term -> both sides cannot be SOC */
               SCIPfreeBufferArray(scip, &lhsvar);
               SCIPfreeBufferArray(scip, &lhscoeff);
               return SCIP_OKAY;
            }
            if (!term->sqrcoeff) /* loose variable */
               continue;

            if (term->sqrcoeff > 0)
            {
               if (lhscount >= nterms-1)
               { /* too many variables on lhs, i.e., all variables seem to have positive coefficient */
                  rhsvar = NULL;
                  i = SCIPhashmapGetNLists(terms);
                  break;
               }
               lhsvar  [lhscount] = var;
               lhscoeff[lhscount] = term->sqrcoeff;
               ++lhscount;
            }
            else if (rhsvar || SCIPisNegative(scip, SCIPvarGetLbGlobal(var)))
            { /* second variable with negative coefficient or neg. coefficient but var not >= 0 -> cannot be SOC */
               rhsvar = NULL;
               i = SCIPhashmapGetNLists(terms);
               break;
            }
            else
            {
               rhsvar   = var;
               rhscoeff = -term->sqrcoeff;
            }
         }
   }
   
   if (rhsvar && lhscount >= 2)
   { /* found SOC constraint, so upgrade to SOC constraint(s) (below) and relax right hand side */
      consdata->rhs = SCIPinfinity(scip);
   }
   else if (SCIPisZero(scip, consdata->lhs - constant))
   { /* if the first failed, try if constraint on lower bound is SOC (using negated coefficients) */
      rhsvar = NULL;
      lhscount = 0;
      for (i = 0; i < SCIPhashmapGetNLists(terms); ++i)
         for (list = SCIPhashmapGetList(terms, i); list; list = SCIPhashmapListGetNext(list))
         {
            var  = (SCIP_VAR*)         SCIPhashmapListGetOrigin(list);
            term = (PresolveQuadTerm*) SCIPhashmapListGetImage(list);

            if (term->lincoeff || term->bilin)
            { /* variable with linear part or variable in bilinear term -> both sides cannot be SOC */
               SCIPfreeBufferArray(scip, &lhsvar);
               SCIPfreeBufferArray(scip, &lhscoeff);
               return SCIP_OKAY;
            }

            if (!term->sqrcoeff) /* loose variable */
               continue;

            if (term->sqrcoeff < 0)
            {
               if (lhscount >= nterms-1)
               { /* too many variables on lhs, i.e., all variables seem to have positive coefficient */
                  rhsvar = NULL;
                  i = SCIPhashmapGetNLists(terms);
                  break;
               }
               lhsvar  [lhscount] =  var;
               lhscoeff[lhscount] = -term->sqrcoeff;
               ++lhscount;
            }
            else if (rhsvar || SCIPisNegative(scip, SCIPvarGetLbGlobal(var)))
            { /* second variable with negative coefficient or neg. coefficient but var not >= 0 -> cannot be SOC */
               rhsvar = NULL;
               i = SCIPhashmapGetNLists(terms);
               break;
            }
            else
            {
               rhsvar   = var;
               rhscoeff = term->sqrcoeff;
            }
         }
      
      if (rhsvar && lhscount >= 2)
      { /* found SOC constraint, so upgrade to SOC constraint(s) (below) and relax left hand side */
         consdata->lhs = -SCIPinfinity(scip);         
      }
   }

   if (rhsvar && lhscount >= 2)
   { /* have SOC constraint */
      SCIPdebugMessage("found constraint %s to be SOC\n", SCIPconsGetName(cons));
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
#endif
      SCIP_CALL( presolveCreateSOC3(scip, lhsvar, lhscoeff, lhscount, rhsvar, rhscoeff, SCIPconsGetName(cons), cons, soc3_nr_auxvars) );
      consdata->soc_added = TRUE;
   }

   SCIPfreeBufferArray(scip, &lhsvar);
   SCIPfreeBufferArray(scip, &lhscoeff);

   return SCIP_OKAY; 
}

/** checks if constraint can be formulated as SOC constraint and creates corresponding SOC.
 * Depending on which side could be upgraded, consdata->lhs or rhs is set to -/+infinity.
 */
static
SCIP_RETCODE presolveTryUpgradeSOC(
   SCIP*           scip,        /**< SCIP data structure */
   SCIP_CONS*      cons,        /**< constraint */
   SCIP_HASHMAP*   terms,       /**< constraint function in form of PresolveQuadTerm's */
   SCIP_Real       constant     /**< constant part of constraint function */
   )
{
   SCIP_CONSDATA* consdata;
   int            nterms;
   SCIP_VAR**     lhsvar;
   SCIP_Real*     lhscoeff = NULL;
   SCIP_Real*     lhsoffset = NULL;
   SCIP_Real      lhsconstant = 0.0;
   int            lhscount = 0;
   SCIP_VAR*      rhsvar = NULL; 
   SCIP_Real      rhscoeff = 0.0;
   SCIP_Real      rhsoffset = 0.0;
   SCIP_HASHMAPLIST* list;
   SCIP_VAR*         var;
   PresolveQuadTerm* term;
   int            i, j;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(terms != NULL);
   assert(!SCIPisInfinity(scip, ABS(constant)));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   nterms = SCIPhashmapGetNEntries(terms);

   SCIP_CALL( SCIPallocBufferArray(scip, &lhsvar,   nterms-1) );

   if (!SCIPisInfinity(scip, consdata->rhs))
   { /* try whether constraint on right hand side is SOC */
      lhsconstant = constant - consdata->rhs;

      for (i = 0; i < SCIPhashmapGetNLists(terms); ++i)
         for (list = SCIPhashmapGetList(terms, i); list; list = SCIPhashmapListGetNext(list))
         {
            var  = (SCIP_VAR*)         SCIPhashmapListGetOrigin(list);
            term = (PresolveQuadTerm*) SCIPhashmapListGetImage(list);
            if (term->bilin)
            { /* variable with bilinear term -> both sides cannot be SOC */
               SCIPfreeBufferArray(scip, &lhsvar);
               SCIPfreeBufferArrayNull(scip, &lhscoeff);
               SCIPfreeBufferArrayNull(scip, &lhsoffset);
               return SCIP_OKAY;
            }
            
            if (!term->sqrcoeff)
            {
               if (term->lincoeff)
               {  /* linear variable, not supported yet @TODO */
                  rhsvar = NULL;
                  i = SCIPhashmapGetNLists(terms);
                  break;
               }
               else
               {  /* loose variable */
                  continue;
               }
            }
            
            if (term->sqrcoeff > 0)
            {
               if (lhscount >= nterms-1)
               { /* too many variables on lhs, i.e., all variables seem to have positive coefficient */
                  rhsvar = NULL;
                  i = SCIPhashmapGetNLists(terms);
                  break;
               }
               
               lhsvar  [lhscount] = var;

               if (term->sqrcoeff != 1.0)
               {
                  if (!lhscoeff)
                  {
                     SCIP_CALL( SCIPallocBufferArray(scip, &lhscoeff, nterms-1) );
                     for (j = 0; j < lhscount; ++j)
                        lhscoeff[j] = 1.0;
                  }
                  lhscoeff[lhscount] = sqrt(term->sqrcoeff);
               }
               else if (lhscoeff)
                  lhscoeff[lhscount] = 1.0;
               
               if (term->lincoeff != 0.0)
               {
                  if (!lhsoffset)
                  {
                     SCIP_CALL( SCIPallocBufferArray(scip, &lhsoffset, nterms-1) );
                     for (j = 0; j < lhscount; ++j)
                        lhsoffset[j] = 0.0;
                  }
                  lhsoffset[lhscount] = term->lincoeff / (2 * term->sqrcoeff);
                  lhsconstant -= term->lincoeff * term->lincoeff / (4 * term->sqrcoeff);
               }
               else if (lhsoffset)
                  lhsoffset[lhscount] = 0.0;
               
               ++lhscount;
            }
            else if (rhsvar || SCIPisLT(scip, SCIPvarGetLbGlobal(var), term->lincoeff / (2*term->sqrcoeff)))
            { /* second variable with negative coefficient -> cannot be SOC */
              /* or lower bound of variable does not ensure positivity of right hand side */
               rhsvar = NULL;
               i = SCIPhashmapGetNLists(terms);
               break;
            }
            else
            { 
               rhsvar    = var;
               rhscoeff  = sqrt(-term->sqrcoeff);
               rhsoffset = term->lincoeff / (2*term->sqrcoeff);
               lhsconstant -= term->lincoeff * term->lincoeff / (4 * term->sqrcoeff);
            }
         }
   }
   
   if (rhsvar && lhscount >= 2 && !SCIPisNegative(scip, lhsconstant))
   { /* found SOC constraint, so upgrade to SOC constraint(s) (below) and relax right hand side */
      consdata->rhs = SCIPinfinity(scip);
   }
   else if (!SCIPisInfinity(scip, -consdata->lhs))
   { /* if the first failed, try if constraint on left hand side is SOC (using negated coefficients) */
      SCIPfreeBufferArrayNull(scip, &lhscoeff);
      SCIPfreeBufferArrayNull(scip, &lhsoffset);
      lhscount = 0;
      rhsvar = NULL;
      
      lhsconstant = consdata->lhs - constant;

      for (i = 0; i < SCIPhashmapGetNLists(terms); ++i)
         for (list = SCIPhashmapGetList(terms, i); list; list = SCIPhashmapListGetNext(list))
         {
            var  = (SCIP_VAR*)         SCIPhashmapListGetOrigin(list);
            term = (PresolveQuadTerm*) SCIPhashmapListGetImage(list);
            if (term->bilin)
            { /* variable with bilinear term -> both sides cannot be SOC */
               SCIPfreeBufferArray(scip, &lhsvar);
               SCIPfreeBufferArrayNull(scip, &lhscoeff);
               SCIPfreeBufferArrayNull(scip, &lhsoffset);
               return SCIP_OKAY;
            }
            
            if (!term->sqrcoeff)
            {
               if (term->lincoeff)
               {  /* linear variable, not supported yet @TODO */
                  rhsvar = NULL;
                  i = SCIPhashmapGetNLists(terms);
                  break;
               }
               else
               {  /* loose variable */
                  continue;
               }
            }
            
            if (term->sqrcoeff < 0)
            {
               if (lhscount >= nterms-1)
               { /* too many variables on lhs, i.e., all variables seem to have negative coefficient */
                  rhsvar = NULL;
                  i = SCIPhashmapGetNLists(terms);
                  break;
               }
               
               lhsvar  [lhscount] = var;

               if (term->sqrcoeff != -1.0)
               {
                  if (!lhscoeff)
                  {
                     SCIP_CALL( SCIPallocBufferArray(scip, &lhscoeff, nterms-1) );
                     for (j = 0; j < lhscount; ++j)
                        lhscoeff[j] = 1.0;
                  }
                  lhscoeff[lhscount] = sqrt(-term->sqrcoeff);
               }
               else if (lhscoeff)
                  lhscoeff[lhscount] = 1.0;
               
               if (term->lincoeff != 0.0)
               {
                  if (!lhsoffset)
                  {
                     SCIP_CALL( SCIPallocBufferArray(scip, &lhsoffset, nterms-1) );
                     for (j = 0; j < lhscount; ++j)
                        lhsoffset[j] = 0.0;
                  }
                  lhsoffset[lhscount] = term->lincoeff / (2 * term->sqrcoeff);
                  lhsconstant += term->lincoeff * term->lincoeff / (4 * term->sqrcoeff);
               }
               else if (lhsoffset)
                  lhsoffset[lhscount] = 0.0;
               
               ++lhscount;
            }
            else if (rhsvar || SCIPisLT(scip, SCIPvarGetLbGlobal(var), term->lincoeff / (2*term->sqrcoeff)))
            { /* second variable with positive coefficient -> cannot be SOC */
              /* or lower bound of variable does not ensure positivity of right hand side */
               rhsvar = NULL;
               i = SCIPhashmapGetNLists(terms);
               break;
            }
            else
            { 
               rhsvar    = var;
               rhscoeff  = sqrt(term->sqrcoeff);
               rhsoffset = term->lincoeff / (2*term->sqrcoeff);
               lhsconstant += term->lincoeff * term->lincoeff / (4 * term->sqrcoeff);
            }
         }
      
      if (rhsvar && lhscount >= 2 && !SCIPisNegative(scip, lhsconstant))
      { /* found SOC constraint, so upgrade to SOC constraint(s) (below) and relax left hand side */
         consdata->lhs = -SCIPinfinity(scip);
      }
   }

   if (rhsvar && lhscount >= 2 && !SCIPisNegative(scip, lhsconstant))
   { /* have SOC constraint */
      SCIP_CONS* soccons;
      SCIPdebugMessage("found constraint %s to be SOC\n", SCIPconsGetName(cons));
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
#endif
  
      SCIP_CALL( SCIPcreateConsSOC(scip, &soccons, SCIPconsGetName(cons),
         lhscount, lhsvar, lhscoeff, lhsoffset, lhsconstant,
         rhsvar, rhscoeff, rhsoffset,
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
         SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintCons(scip, soccons, NULL) );
#endif
      SCIP_CALL( SCIPaddCons(scip, soccons) );
      SCIP_CALL( SCIPreleaseCons(scip, &soccons) );
   }

   SCIPfreeBufferArray(scip, &lhsvar);
   SCIPfreeBufferArrayNull(scip, &lhscoeff);
   SCIPfreeBufferArrayNull(scip, &lhsoffset);
   
   return SCIP_OKAY;
}
#endif

static
SCIP_RETCODE checkCurvature(
   SCIP*           scip,       /**< SCIP data structure */
   SCIP_CONSHDLR*  conshdlr,   /**< quadratic constraint handler */
   SCIP_CONS*      cons        /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;
   int            i, n, nn;
   double*        matrix;
   SCIP_HASHMAP*  var2index;
#ifdef WITH_LAPACK
   int            row, col;
   double*        alleigval;
#endif

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   n = consdata->n_quadvar;

   if (n == 1)
   {
      assert(consdata->n_bilin == 0);
      consdata->is_convex  = !SCIPisNegative(scip, consdata->quadsqrcoeff[0]);
      consdata->is_concave = !SCIPisPositive(scip, consdata->quadsqrcoeff[0]);
      return SCIP_OKAY;
   }

   if (n == 0)
   {
      consdata->is_convex  = TRUE;
      consdata->is_concave = TRUE;
      return SCIP_OKAY;
   }

   if (consdata->n_bilin == 0)
   {
      consdata->is_convex  = TRUE;
      consdata->is_concave = TRUE;
      for (i = 0; i < n; ++i)
      {
         consdata->is_convex  &= !SCIPisNegative(scip, consdata->quadsqrcoeff[i]);
         consdata->is_concave &= !SCIPisPositive(scip, consdata->quadsqrcoeff[i]);
      }
      return SCIP_OKAY;
   }

   if (n == 2)
   { /* compute eigenvalues by hand */
      SCIP_Real d;
      assert(consdata->n_bilin == 1);
      d = (consdata->quadsqrcoeff[0] - consdata->quadsqrcoeff[1]);
      d = sqrt(d*d + consdata->bilincoeff[0]*consdata->bilincoeff[0]);
      consdata->is_convex  = !SCIPisNegative(scip, -(consdata->quadsqrcoeff[0] + consdata->quadsqrcoeff[1]) - d);
      consdata->is_concave = !SCIPisPositive(scip, -(consdata->quadsqrcoeff[0] + consdata->quadsqrcoeff[1]) + d);
      return SCIP_OKAY;
   }

   /* lower triangular of quadratic term matrix, scaled by box diameter */
   nn = n * n;
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix, nn) );
   memset(matrix, 0, nn*sizeof(double));

   consdata->is_convex  = TRUE;
   consdata->is_concave = TRUE;

   SCIP_CALL( SCIPhashmapCreate(&var2index, SCIPblkmem(scip), n) );
   for (i = 0; i < n; ++i)
   {
      if (consdata->n_adjbilin[i])
      {
         SCIP_CALL( SCIPhashmapInsert(var2index, consdata->quadvar[i], (void*)(size_t)i) );
         matrix[i*n + i] = consdata->quadsqrcoeff[i];
      }
      else
      {
         consdata->is_convex  &= !SCIPisNegative(scip, consdata->quadsqrcoeff[i]);
         consdata->is_concave &= !SCIPisPositive(scip, consdata->quadsqrcoeff[i]);
      }
   }

   if (!consdata->is_convex && !consdata->is_concave)
   {
      SCIPfreeBufferArray(scip, &matrix);
      SCIPhashmapFree(&var2index);
      return SCIP_OKAY;
   }

#ifdef WITH_LAPACK
   for (i = 0; i < consdata->n_bilin; ++i)
   {
      assert(SCIPhashmapExists(var2index, consdata->bilinvar1[i]));
      assert(SCIPhashmapExists(var2index, consdata->bilinvar2[i]));
      row = (int)(size_t)SCIPhashmapGetImage(var2index, consdata->bilinvar1[i]);
      col = (int)(size_t)SCIPhashmapGetImage(var2index, consdata->bilinvar2[i]);
      if (row < col)
         matrix[row * n + col] = consdata->bilincoeff[i]/2;
      else
         matrix[col * n + row] = consdata->bilincoeff[i]/2;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &alleigval, n) );
   /* TODO can we compute only min and max eigval?
      TODO can we estimate the numerical error? */
   if (LapackDsyev(scip, FALSE, n, matrix, alleigval) != SCIP_OKAY)
   {
      SCIPwarningMessage("Failed to compute eigenvalues of quadratic coefficient matrix of constraint %s. Assuming matrix is indefinite.\n", SCIPconsGetName(cons));
      consdata->is_convex  = FALSE;
      consdata->is_concave = FALSE;
   }
   else
   {
      consdata->is_convex  &= !SCIPisNegative(scip, alleigval[0]);
      consdata->is_concave &= !SCIPisPositive(scip, alleigval[n-1]);
   }
   SCIPfreeBufferArray(scip, &alleigval);
   
#else
   consdata->is_convex  = FALSE;
   consdata->is_concave = FALSE;
#endif
   
   SCIPhashmapFree(&var2index);
   SCIPfreeBufferArray(scip, &matrix);

   return SCIP_OKAY;
}

/** locks a linear variable in a constraint */
static
SCIP_RETCODE lockLinearVariable(
   SCIP*       scip,   /**< SCIP data structure */
   SCIP_CONS*  cons,   /**< constraint where to lock a variable */
   SCIP_VAR*   var,    /**< variable to lock */
   SCIP_Real   coeff   /**< coefficient of variable in constraint */
   )
{
   SCIP_CONSDATA* consdata;
   
   assert(scip  != NULL);
   assert(cons  != NULL);
   assert(var   != NULL);
   assert(coeff != 0.0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if (coeff > 0)
      SCIP_CALL( SCIPlockVarCons(scip, var, cons, !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip,  consdata->rhs)) );
   else
      SCIP_CALL( SCIPlockVarCons(scip, var, cons, !SCIPisInfinity(scip,  consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );

   return SCIP_OKAY;
}

/** unlocks a linear variable in a constraint */
static
SCIP_RETCODE unlockLinearVariable(
   SCIP*       scip,   /**< SCIP data structure */
   SCIP_CONS*  cons,   /**< constraint where to unlock a variable */
   SCIP_VAR*   var,    /**< variable to unlock */
   SCIP_Real   coeff   /**< coefficient of variable in constraint */
   )
{
   SCIP_CONSDATA* consdata;
   
   assert(scip  != NULL);
   assert(cons  != NULL);
   assert(var   != NULL);
   assert(coeff != 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if (coeff > 0)
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons, !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip,  consdata->rhs)) );
   else
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons, !SCIPisInfinity(scip,  consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );

   return SCIP_OKAY;
}

/** Sets bounds for variables in not evidently convex terms to some predefined value.
 */
static
SCIP_RETCODE boundUnboundedVars(
   SCIP*           scip,        /**< SCIP data structure */
   SCIP_CONS*      cons,        /**< constraint */
   SCIP_Real       bound,       /**< value to use for bound */
   int*            n_chgbnds    /**< buffer where to add the number of bound changes, or NULL */
)
{
   SCIP_Bool      infeasible;
   SCIP_CONSDATA* consdata;
   int            i;
   
   assert(scip != NULL);
   assert(cons != NULL);

   if (SCIPisInfinity(scip, bound))
      return SCIP_OKAY;

   consdata =  SCIPconsGetData(cons);
   assert(consdata != NULL);

   for (i = 0; i < consdata->n_quadvar; ++i)
   {
      if (!consdata->n_adjbilin[i] &&
         (SCIPisInfinity(scip,  consdata->rhs) || consdata->quadsqrcoeff[i] > 0) &&
         (SCIPisInfinity(scip, -consdata->lhs) || consdata->quadsqrcoeff[i] < 0))
         continue; /* skip evidently convex terms */

      if (SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->quadvar[i])))
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "set lower bound of %s to %g\n", SCIPvarGetName(consdata->quadvar[i]), -bound);
         SCIP_CALL( SCIPtightenVarLb(scip, consdata->quadvar[i], -bound, FALSE, &infeasible, NULL) );
         assert(!infeasible);
         if (n_chgbnds)
            ++*n_chgbnds;
      }

      if (SCIPisInfinity(scip,  SCIPvarGetUbLocal(consdata->quadvar[i])))
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "set upper bound of %s to %g\n", SCIPvarGetName(consdata->quadvar[i]),  bound);
         SCIP_CALL( SCIPtightenVarUb(scip, consdata->quadvar[i],  bound, FALSE, &infeasible, NULL) );
         assert(!infeasible);
         if (n_chgbnds) 
            ++*n_chgbnds;
      }
   }

   return SCIP_OKAY;
}

static
SCIP_Real getGradientNorm(
   SCIP*       scip,     /**< SCIP data structure */
   SCIP_CONS*  cons,     /**< constraint */
   SCIP_SOL*   sol       /**< solution or NULL if LP solution should be used */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real      norm = 0.0;
   SCIP_Real      g;
   int            i, j, k;
   SCIP_VAR*      var;
   
   assert(scip != NULL);
   assert(cons != NULL);
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   /* TODO allow also other norms than euclidean, maybe read separating/efficacynorm */
   
   for (i = 0; i < consdata->n_linvar; ++i)
      norm += consdata->lincoeff[i] * consdata->lincoeff[i];
   
   for (i = 0; i < consdata->n_quadvar; ++i)
   {
      var = consdata->quadvar[i];
      assert(!SCIPisInfinity(scip, ABS(SCIPgetSolVal(scip, sol, var))));
      g  =     consdata->quadlincoeff[i];
      g += 2 * consdata->quadsqrcoeff[i] * SCIPgetSolVal(scip, sol, var);
      for (j = 0; j < consdata->n_adjbilin[i]; ++j)
      {
         k = consdata->adjbilin[i][j];
         if (consdata->bilinvar1[k] == var)
            g += consdata->bilincoeff[k] * SCIPgetSolVal(scip, sol, consdata->bilinvar2[k]);
         else
            g += consdata->bilincoeff[k] * SCIPgetSolVal(scip, sol, consdata->bilinvar1[k]);
      }
      norm += g*g;
   }
   
   return sqrt(norm);
}

static
SCIP_RETCODE computeViolation(
   SCIP*       scip,       /**< SCIP data structure */
   SCIP_CONS*  cons,       /**< constraint */
   SCIP_SOL*   sol,        /**< solution or NULL if LP solution should be used */
   SCIP_Bool   do_scaling  /**< whether we should scale the violation by the gradient of the quadratic function */ 
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real      val = 0.0;
   SCIP_Real      varval = 0.0;
   int            i, j;
   
   assert(scip != NULL);
   assert(cons != NULL);
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   /* TODO take better care of variables at +/- infinity: e.g., run instance waste in debug mode with a short timelimit (30s) */
   for (i = 0; i < consdata->n_linvar; ++i)
   {
      if (SCIPisInfinity(scip, ABS(SCIPgetSolVal(scip, sol, consdata->linvar[i]))))
      {
         consdata->lhsviol = consdata->rhsviol = SCIPinfinity(scip);
         return SCIP_OKAY;
      }
      val += consdata->lincoeff[i] * SCIPgetSolVal(scip, sol, consdata->linvar[i]);
   }

   for (j = 0; j < consdata->n_quadvar; ++j)
   {
      varval = SCIPgetSolVal(scip, sol, consdata->quadvar[j]);
      if (SCIPisInfinity(scip, ABS(varval)))
      {
         consdata->lhsviol = consdata->rhsviol = SCIPinfinity(scip);
         return SCIP_OKAY;
      }
      val   += (consdata->quadlincoeff[j] + consdata->quadsqrcoeff[j] * varval) * varval;
   }
   
   for (j = 0; j < consdata->n_bilin; ++j)
      val += consdata->bilincoeff[j] * SCIPgetSolVal(scip, sol, consdata->bilinvar1[j]) * SCIPgetSolVal(scip, sol, consdata->bilinvar2[j]);

   if ( val < consdata->lhs && !SCIPisInfinity(scip, -consdata->lhs) )
      consdata->lhsviol = consdata->lhs - val;
   else
      consdata->lhsviol = 0.0;
   
   if ( val > consdata->rhs && !SCIPisInfinity(scip,  consdata->rhs) )
      consdata->rhsviol = val - consdata->rhs;
   else
      consdata->rhsviol = 0.0;
   
   if (do_scaling && (consdata->lhsviol || consdata->rhsviol))
   {
      SCIP_Real norm = getGradientNorm(scip, cons, sol);
      if (norm > 1.)
      { /* TODO scale only if > 1., or should it be larger SCIPsumepsilon? */
         consdata->lhsviol /= norm;
         consdata->rhsviol /= norm;
      }
   }
   
   return SCIP_OKAY;
}

static
SCIP_RETCODE computeViolations(
   SCIP*        scip,        /**< SCIP data structure */
   SCIP_CONS**  conss,       /**< constraints */
   int          nconss,      /**< number of constraints */
   SCIP_SOL*    sol,         /**< solution or NULL if LP solution should be used */
   SCIP_Bool    do_scaling,  /**< whether to do scaling when computing violation */
   SCIP_CONS**  maxviolcon   /**< buffer to store constraint with largest violation, or NULL if solution is feasible */
)
{
   SCIP_CONSDATA* consdata;
   SCIP_Real      viol;
   SCIP_Real      maxviol = 0.0;
   int            c;

   assert(scip != NULL);
   assert(nconss == 0 || conss != NULL);
   assert(maxviolcon != NULL);
   
   *maxviolcon = NULL;
   
   for (c = 0; c < nconss; ++c)
   {
      assert(conss[c] != NULL);
      
      SCIP_CALL( computeViolation(scip, conss[c], sol, do_scaling) );

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      
      viol = MAX(consdata->lhsviol, consdata->rhsviol);
      if (viol > maxviol && SCIPisFeasPositive(scip, viol))
      {
         maxviol = viol;
         *maxviolcon = conss[c];
      }
   }
   
   return SCIP_OKAY;
}

/** Adds range of quadratic term w.r.t. local bounds to given interval.
 */
static
SCIP_RETCODE addQuadRange(
   SCIP*            scip,          /**< SCIP data structure */
   SCIP_CONS*       cons,          /**< constraint */
   SCIP_Real        intervalinfty, /**< value of infinity to use for interval operations */
   SCIP_INTERVAL*   resultant,     /**< interval where to add to */
   SCIP_VAR*        except         /**< a variable to skip in evaluation, NULL if nothing should be skipped */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_INTERVAL  lincoeff; /* linear interval coefficient in quadratic form */
   SCIP_INTERVAL  xrng;
   SCIP_INTERVAL  tmp;
   int            i, j, k;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(resultant != NULL);
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   for (i = 0; i < consdata->n_quadvar; ++i)
   {
      if (consdata->quadvar[i] == except)
         continue;
      
      SCIPintervalSetBounds(&xrng, 
         -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->quadvar[i]), SCIPvarGetUbLocal(consdata->quadvar[i]))),
          infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->quadvar[i]), SCIPvarGetUbLocal(consdata->quadvar[i])))
         );
      SCIPintervalSet(&lincoeff, consdata->quadlincoeff[i]);

      for (j = 0; j < consdata->n_adjbilin[i]; ++j)
      {
         k = consdata->adjbilin[i][j];
         if (consdata->bilinvar1[k] != consdata->quadvar[i])
            continue; /* handle this term later */
         if (consdata->bilinvar2[k] == except)
            continue; /* variable is skipped */
         
         SCIPintervalSetBounds(&tmp, 
            -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinvar2[k]), SCIPvarGetUbLocal(consdata->bilinvar2[k]))),
             infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinvar2[k]), SCIPvarGetUbLocal(consdata->bilinvar2[k])))
            );
         SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilincoeff[k]);
         SCIPintervalAdd(intervalinfty, &lincoeff, lincoeff, tmp);
      }
      
      SCIPintervalQuad(intervalinfty, &tmp, consdata->quadsqrcoeff[i], lincoeff, xrng);
      assert(SCIPintervalGetSup(tmp) > -intervalinfty);
      assert(SCIPintervalGetInf(tmp) <  intervalinfty);
      SCIPintervalAdd(intervalinfty, resultant, *resultant, tmp);
   }
   
   return SCIP_OKAY;
}

/** checks by interval analysis whether a violated constraint is infeasible
 * If lbviol and ubviol is below feasibility tolerance, the check is skipped.
 * If isfeasible is set to false, then constraint is infeasible w.r.t. current local bounds.
 * If isfeasible is set to true, then this gives no information.
 */
static
SCIP_RETCODE isIntervalFeasible(
   SCIP*       scip,          /**< SCIP data structure */
   SCIP_CONS*  cons,          /**< constraint to check */
   SCIP_Bool*  isfeasible     /**< buffer to store the result */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_INTERVAL  val;
   SCIP_Real      intervalinfty;
   int            i;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(isfeasible != NULL);
   assert(SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVING && SCIPgetStage(scip) < SCIP_STAGE_SOLVED);
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   *isfeasible = TRUE;
   intervalinfty = 1000 * SCIPinfinity(scip) * SCIPinfinity(scip);
   
   if (!SCIPisPositive(scip, consdata->lhsviol) && !SCIPisPositive(scip, consdata->rhsviol))
      return SCIP_OKAY; /* obviously we have a feasible point */

   if (SCIPintervalIsEmpty(consdata->quadrange))
   { /* need to update quadrange */
      SCIPintervalSet(&consdata->quadrange, 0.0);
      SCIP_CALL( addQuadRange(scip, cons, intervalinfty, &consdata->quadrange, NULL) );
   }

   val = consdata->quadrange;
   
   for (i = 0; i < consdata->n_linvar; ++i)
   {
      if (SCIPintervalIsEmpty(consdata->linrange[i]))
      { /* need to update linrange for var. i */
         SCIPintervalSetBounds(&consdata->linrange[i], 
            -infty2infty(SCIPinfinity(scip), intervalinfty, -SCIPvarGetLbLocal(consdata->linvar[i])), 
             infty2infty(SCIPinfinity(scip), intervalinfty,  SCIPvarGetUbLocal(consdata->linvar[i]))
            );
         SCIPintervalMulScalar(intervalinfty, &consdata->linrange[i], consdata->linrange[i], consdata->lincoeff[i]);
      }
      SCIPintervalAdd(intervalinfty, &val, val, consdata->linrange[i]);
   }
   
   if (SCIPisFeasGT(scip, consdata->lhs, SCIPintervalGetSup(val)) || SCIPisFeasLT(scip, consdata->rhs, SCIPintervalGetInf(val)))
   {
      SCIPdebugMessage("interval arithmetic found constraint %s infeasible: bounds = [%g, %g], interval = [%g, %g]\n", SCIPconsGetName(cons), consdata->lhs, consdata->rhs, SCIPintervalGetInf(val), SCIPintervalGetSup(val));
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      *isfeasible = FALSE;
   }
   
   return SCIP_OKAY;
}

/** checks by interval analysis whether a set of constraints is infeasible
 * If isfeasible is set to false, then one constraint is infeasible w.r.t. current local bounds.
 * If isfeasible is set to true, then this gives no information.
 */
static
SCIP_RETCODE areIntervalFeasible(
   SCIP*        scip,          /**< SCIP data structure */
   SCIP_CONS**  conss,         /**< constraints to check */
   int          nconss,        /**< number of constraints to check */
   SCIP_CONS*   firstcons,     /**< constraint to check first, can be NULL */
   SCIP_Bool*   isfeasible     /**< buffer to store the result */
   )
{
   int c;
   
   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   assert(isfeasible != NULL);

   if (firstcons)
   {
      SCIP_CALL( isIntervalFeasible(scip, firstcons, isfeasible) );
      if (!*isfeasible)
         return SCIP_OKAY;
   }

   for (c = 0; *isfeasible && c < nconss; ++c)
   {
      if (conss[c] == firstcons)
         continue;
      SCIP_CALL( isIntervalFeasible(scip, conss[c], isfeasible) );
   }

   return SCIP_OKAY;
}


/** generates a cut based on linearization (if convex) or McCormick (if nonconvex)
 */
static
SCIP_RETCODE generateCut(
   SCIP*           scip,       /**< SCIP data structure */
   SCIP_CONS*      cons,       /**< constraint */
   SCIP_SOL*       sol,        /**< solution to separate, or NULL if LP solution should be used */
   SCIP_BOUNDTYPE  violbound,  /**< for which bound a cut should be generated */
   SCIP_ROW**      row,        /**< storage for cut */
   SCIP_Real       maxrange    /**< maximal range allowed */
)
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool      is_convex;
   SCIP_Bool      is_global;
   SCIP_Real      coeff, rowcoeff, xcoeff, ycoeff;
   SCIP_Real      bnd, bnd_;
   SCIP_VAR*      x;
   SCIP_VAR*      y;
   SCIP_Real      xlb, xub, xval;
   SCIP_Real      ylb, yub, yval;
   int            j;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(row  != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   is_convex = (violbound == SCIP_BOUNDTYPE_LOWER) ? consdata->is_concave : consdata->is_convex;
   is_global = SCIPconsIsGlobal(cons) && is_convex;

   SCIP_CALL( SCIPcreateEmptyRow(scip, row, "cut", -SCIPinfinity(scip), SCIPinfinity(scip), !is_global /* locally */, TRUE /* modifiable */, TRUE /* removable */ ) );
   bnd = (violbound == SCIP_BOUNDTYPE_LOWER) ? consdata->lhs : consdata->rhs;
   assert(!SCIPisInfinity(scip, ABS(bnd)));

   /* add linear part */
   SCIP_CALL( SCIPaddVarsToRow(scip, *row, consdata->n_linvar, consdata->linvar, consdata->lincoeff) );
   /* TODO should we buffer the coefficients of the quadratic variables before adding them to the row? */

   if (is_convex)
   {  /* do first-order taylor for each term */
      for (j = 0; j < consdata->n_quadvar; ++j)
      { /* linear term + linearization of square term */
         rowcoeff = consdata->quadlincoeff[j];

         x    = consdata->quadvar[j];
         xval = SCIPgetSolVal(scip, sol, x);
         /* can happen when called from initlp */
         if (xval < SCIPvarGetLbLocal(x))
            xval = SCIPvarGetLbLocal(x);
         else if (xval > SCIPvarGetUbLocal(x))
            xval = SCIPvarGetUbLocal(x);
         if (SCIPisInfinity(scip, ABS(xval)))
         {
            SCIPdebugMessage("skip linearization of square term in constraint %s because var %s is at infinity\n", SCIPconsGetName(cons), SCIPvarGetName(x));
            SCIP_CALL( SCIPreleaseRow(scip, row) );
            return SCIP_OKAY;
         }

         if (consdata->n_bilin || SCIPvarGetType(x) == SCIP_VARTYPE_CONTINUOUS || SCIPisIntegral(scip, xval))
         {
            rowcoeff += 2*consdata->quadsqrcoeff[j]*xval;
            bnd      +=   consdata->quadsqrcoeff[j]*xval*xval;
         }
         else
         { /* if variable is discrete but fractional and there are no bilinear terms, try to be more clever */
            /* TODO: could we do something similar even if there are bilinear terms? */
            SCIP_Real f = SCIPfloor(scip, xval);
            rowcoeff += consdata->quadsqrcoeff[j] * (2*f+1);
            bnd      += consdata->quadsqrcoeff[j] * f * (f+1);
         }

         SCIP_CALL( SCIPaddVarToRow(scip, *row, x, rowcoeff) );
      }

      for (j = 0; j < consdata->n_bilin; ++j)
      { /* linearization of bilinear terms */
         x    = consdata->bilinvar1[j];
         xval = SCIPgetSolVal(scip, sol, x);
         if (xval < SCIPvarGetLbLocal(x))
            xval = SCIPvarGetLbLocal(x);
         else if (xval > SCIPvarGetUbLocal(x))
            xval = SCIPvarGetUbLocal(x);
         if (SCIPisInfinity(scip, ABS(xval)))
         {
            SCIPdebugMessage("skip linearization of bilinear term in constraint %s because var %s is at infinity\n", SCIPconsGetName(cons), SCIPvarGetName(x));
            SCIP_CALL( SCIPreleaseRow(scip, row) );
            return SCIP_OKAY;
         }

         y    = consdata->bilinvar2[j];
         yval = SCIPgetSolVal(scip, sol, y);
         if (yval < SCIPvarGetLbLocal(y))
            yval = SCIPvarGetLbLocal(y);
         else if (yval > SCIPvarGetUbLocal(y))
            yval = SCIPvarGetUbLocal(y);
         if (SCIPisInfinity(scip, ABS(yval)))
         {
            SCIPdebugMessage("skip linearization of bilinear term in constraint %s because var %s is at infinity\n", SCIPconsGetName(cons), SCIPvarGetName(y));
            SCIP_CALL( SCIPreleaseRow(scip, row) );
            return SCIP_OKAY;
         }

         coeff = consdata->bilincoeff[j];

         SCIP_CALL( SCIPaddVarToRow(scip, *row, x, coeff * yval) );
         SCIP_CALL( SCIPaddVarToRow(scip, *row, y, coeff * xval) );
         bnd += coeff * xval * yval;
      }
   }
   else
   { /* underestimate and linearize each term separately -> McCormick */
      for (j = 0; j < consdata->n_quadvar; ++j)
      {
         rowcoeff = consdata->quadlincoeff[j];

         x    = consdata->quadvar[j];
         xval = SCIPgetSolVal(scip, sol, x);
         xlb  = SCIPvarGetLbLocal(x);
         xub  = SCIPvarGetUbLocal(x);
         if (xval < xlb)
            xval = xlb;
         else if (xval > xub)
            xval = xub;
         if (SCIPisInfinity(scip, ABS(xval)))
         {
            SCIPdebugMessage("skip underestimator of square term in constraint %s because var %s is at infinity\n", SCIPconsGetName(cons), SCIPvarGetName(x));
            SCIP_CALL( SCIPreleaseRow(scip, row) );
            return SCIP_OKAY;
         }

         /* linearization of square term */
         coeff = consdata->quadsqrcoeff[j];

         if ((violbound == SCIP_BOUNDTYPE_LOWER && coeff <= 0) ||
             (violbound == SCIP_BOUNDTYPE_UPPER && coeff >  0))
         { /* convex -> linearize */
            if (SCIPvarGetType(x) == SCIP_VARTYPE_CONTINUOUS || SCIPisIntegral(scip, xval))
            {
               rowcoeff += 2*coeff*xval;
               bnd      +=   coeff*xval*xval;
            }
            else
            { /* if variable is discrete but fractional, try to be more clever */
               SCIP_Real f = SCIPfloor(scip, xval);
               rowcoeff += coeff*(2*f+1);
               bnd      += coeff*f*(f+1);
            }
         }
         else
         { /* not convex -> secand approximation */
            if (SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub))
            {
               SCIPdebugMessage("skip secand approx of square term in constraint %s because var %s is at infinity\n", SCIPconsGetName(cons), SCIPvarGetName(x));
               SCIP_CALL( SCIPreleaseRow(scip, row) );
               return SCIP_OKAY;
            }

            rowcoeff += coeff * (xlb+xub);
            bnd      += coeff * xlb * xub;
         }

         SCIP_CALL( SCIPaddVarToRow(scip, *row, x, rowcoeff) );
      }

      for (j = 0; j < consdata->n_bilin; ++j)
      {
         x    = consdata->bilinvar1[j];
         xval = SCIPgetSolVal(scip, sol, x);
         xlb  = SCIPvarGetLbLocal(x);
         xub  = SCIPvarGetUbLocal(x);
         if (xval < xlb)
            xval = xlb;
         else if (xval > xub)
            xval = xub;
         if (SCIPisInfinity(scip, ABS(xval)))
         {
            SCIPdebugMessage("skip underestimator of bilinear term in constraint %s because var %s is at infinity\n", SCIPconsGetName(cons), SCIPvarGetName(x));
            SCIP_CALL( SCIPreleaseRow(scip, row) );
            return SCIP_OKAY;
         }

         y    = consdata->bilinvar2[j];
         yval = SCIPgetSolVal(scip, sol, y);
         ylb  = SCIPvarGetLbLocal(y);
         yub  = SCIPvarGetUbLocal(y);
         if (yval < ylb)
            yval = ylb;
         else if (yval > yub)
            yval = yub;
         if (SCIPisInfinity(scip, ABS(xval)))
         {
            SCIPdebugMessage("skip underestimator of bilinear term in constraint %s because var %s is at infinity\n", SCIPconsGetName(cons), SCIPvarGetName(y));
            SCIP_CALL( SCIPreleaseRow(scip, row) );
            return SCIP_OKAY;
         }

         coeff = consdata->bilincoeff[j];
         if (violbound == SCIP_BOUNDTYPE_LOWER)
            coeff = -coeff;

         if (coeff > 0)
         {
            if (!SCIPisInfinity(scip, -xlb) && !SCIPisInfinity(scip, -ylb) &&
               (SCIPisInfinity(scip,  xub) ||  SCIPisInfinity(scip,  yub) ||
                  (xub-xlb)*yval + (yub-ylb)*xval <= xub*yub - xlb*ylb) )
            {
               xcoeff = coeff * ylb;
               ycoeff = coeff * xlb;
               bnd_   = coeff * xlb * ylb;
            }
            else if (!SCIPisInfinity(scip, xub) && !SCIPisInfinity(scip, yub))
            {
               xcoeff = coeff * yub;
               ycoeff = coeff * xub;
               bnd_   = coeff * xub * yub;
            }
            else
            {
               SCIPdebugMessage("skip underestimator of bilinear term in constraint %s because var %s or %s is unbounded\n", SCIPconsGetName(cons), SCIPvarGetName(x), SCIPvarGetName(y));
               SCIP_CALL( SCIPreleaseRow(scip, row) );
               return SCIP_OKAY;
            }
         }
         else
         { /* coeff < 0 */
            if (!SCIPisInfinity(scip,  xub) && !SCIPisInfinity(scip, -ylb) &&
               (SCIPisInfinity(scip, -xlb) ||  SCIPisInfinity(scip,  yub) ||
                  (xub-xlb)*yval - (yub-ylb)*xval <= xub*ylb - xlb*yub) )
            {
               xcoeff = coeff * ylb;
               ycoeff = coeff * xub;
               bnd_   = coeff * xub * ylb;
            }
            else if (!SCIPisInfinity(scip, -xlb) && !SCIPisInfinity(scip, yub))
            {
               xcoeff = coeff * yub;
               ycoeff = coeff * xlb;
               bnd_   = coeff * xlb * yub;
            }
            else
            {
               SCIPdebugMessage("skip underestimator of bilinear term in constraint %s because var %s or %s is unbounded\n", SCIPconsGetName(cons), SCIPvarGetName(x), SCIPvarGetName(y));
               SCIP_CALL( SCIPreleaseRow(scip, row) );
               return SCIP_OKAY;
            }
         }
         if (violbound == SCIP_BOUNDTYPE_LOWER)
         {
            xcoeff = -xcoeff;
            ycoeff = -ycoeff;
            bnd_   = -bnd_;
         }
         SCIP_CALL( SCIPaddVarToRow(scip, *row, x, xcoeff) );
         SCIP_CALL( SCIPaddVarToRow(scip, *row, y, ycoeff) );
         bnd += bnd_;
      }
   }

   SCIPdebugMessage(" -> found cut rhs=%f, min=%f, max=%f range=%g\n",
       ABS(bnd),
       SCIPgetRowMinCoef(scip, *row), SCIPgetRowMaxCoef(scip, *row),
       SCIPgetRowMaxCoef(scip, *row)/SCIPgetRowMinCoef(scip, *row));

   if (SCIPisInfinity(scip, ABS(bnd)))
   { /* seems to be a numerically bad cut */
      SCIPdebugMessage("skip cut for constraint %s because of very large left or right hand side: %g\n", SCIPconsGetName(cons), bnd);
      SCIP_CALL( SCIPreleaseRow(scip, row) );
      return SCIP_OKAY;
   }

   if (SCIPisGT(scip, SCIPgetRowMaxCoef(scip, *row)/SCIPgetRowMinCoef(scip, *row), maxrange))
   { /* seems to be a numerically bad cut */
      SCIPdebugMessage("skip cut for constraint %s because of very large range: %g\n", SCIPconsGetName(cons), SCIPgetRowMaxCoef(scip, *row)/SCIPgetRowMinCoef(scip, *row));
      SCIP_CALL( SCIPreleaseRow(scip, row) );
      return SCIP_OKAY;
   }

   if (violbound == SCIP_BOUNDTYPE_LOWER)
      SCIP_CALL( SCIPchgRowLhs(scip, *row, bnd) );
   else
      SCIP_CALL( SCIPchgRowRhs(scip, *row, bnd) );

   return SCIP_OKAY;
}

/** tries to separate solution or LP solution by a linear cut
 * assumes that constraint violations have been computed 
 */
static
SCIP_RETCODE separatePoint(
   SCIP*          scip,          /**< SCIP data structure */
   SCIP_CONSHDLR* conshdlr,      /**< quadratic constraints handler */
   SCIP_CONS**    conss,         /**< constraints */
   int            nconss,        /**< number of constraints */
   int            nusefulconss,  /**< number of constraints that seem to be useful */
   SCIP_SOL*      sol,           /**< solution to separate, or NULL if LP solution should be used */
   SCIP_RESULT*   result,        /**< result of separation */
   SCIP_Bool      addweakcuts    /**< whether also weak (only slightly violated) cuts should be added in a nonconvex constraint */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          efficacy;
   SCIP_BOUNDTYPE     violbound;
   int                c;
   SCIP_ROW*          row;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nusefulconss <= nconss);
   assert(result != NULL);
   
   *result = SCIP_FEASIBLE;
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for (c = 0; c < nconss; ++c)
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if (SCIPisFeasPositive(scip, consdata->lhsviol) || SCIPisFeasPositive(scip, consdata->rhsviol))
      {
         /* we are not feasible anymore */
         if (*result == SCIP_FEASIBLE)
            *result = SCIP_DIDNOTFIND;

         violbound = SCIPisFeasPositive(scip, consdata->lhsviol) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;

         /* generate cut */
         SCIP_CALL( generateCut(scip, conss[c], sol, violbound, &row, conshdlrdata->cutmaxrange) );
         if (!row) /* failed to generate cut */
            continue;

         efficacy = SCIPgetCutEfficacy(scip, sol, row);

         if (efficacy > conshdlrdata->mincutefficacy ||  /* ''strong'' cut */
            (SCIPisFeasPositive(scip, efficacy) &&       /* ''weak'' cut, use only if */ 
               (addweakcuts ||   /* flag is set, or */
                  (violbound == SCIP_BOUNDTYPE_UPPER && consdata->is_convex ) || /* convex  constraint, or */
                  (violbound == SCIP_BOUNDTYPE_LOWER && consdata->is_concave)    /* concave constraint */
               )))
         { /* cut cuts off solution */
            SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE /* forcecut */) );
            *result = SCIP_SEPARATED;
            SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );
         }

         SCIP_CALL( SCIPreleaseRow (scip, &row) );
      }

      /* enforce only useful constraints
       * others are only checked and enforced if we are still feasible or have not found a separating cut yet
       */ 
      if (c >= nusefulconss && *result == SCIP_SEPARATED)
         break;
   }

   return SCIP_OKAY;
}

/** computes the infeasibilities of variables from the convexification gaps in the constraints and notifies the branching rule about them
 */
static
SCIP_RETCODE registerVariableInfeasibilities(
   SCIP*          scip,          /**< SCIP data structure */
   SCIP_CONSHDLR* conshdlr,      /**< constraint handler */
   SCIP_CONS**    conss,         /**< constraints to check */
   int            nconss,        /**< number of constraints to check */
   int*           nnotify        /**< counter for number of notifications performed */
)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          viol = 0.0;
   int                c, j;
   SCIP_Real          xlb, xub, xval;
   SCIP_Real          ylb, yub, yval;
   SCIP_Real          gap;
   SCIP_Real          coeff_;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
#ifdef WITH_CONSBRANCHNL
   assert(conshdlrdata->branchnl != NULL);
#endif
   
   *nnotify = 0;

   for (c = 0; c < nconss; ++c)
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      SCIPdebugMessage("con %s violation: %g %g  convex: %d %d\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol, consdata->is_convex, consdata->is_concave);
      
      if (!consdata->n_quadvar)
         continue;
      
      if (SCIPisFeasPositive(scip, consdata->lhsviol) && !consdata->is_concave)
         viol = consdata->lhsviol;
      else if (SCIPisFeasPositive(scip, consdata->rhsviol) && !consdata->is_convex)
         viol = consdata->rhsviol;
      else
         continue;
      SCIPdebugMessage("con %s violation: %g %g  convex: %d %d\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol, consdata->is_convex, consdata->is_concave);
      
      for (j = 0; j < consdata->n_quadvar; ++j)
      { /* square terms */
         if ((SCIPisFeasPositive(scip, consdata->rhsviol) && consdata->quadsqrcoeff[j] < 0) ||
             (SCIPisFeasPositive(scip, consdata->lhsviol) && consdata->quadsqrcoeff[j] > 0))
         {
            xlb  = SCIPvarGetLbLocal(consdata->quadvar[j]);
            xub  = SCIPvarGetUbLocal(consdata->quadvar[j]);
            if (SCIPisEQ(scip, xlb, xub))
               continue;
            
            xval = SCIPgetSolVal(scip, NULL, consdata->quadvar[j]);
            
            if (SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub))
               gap = SCIPinfinity(scip);
            else if (xval < xlb || xval > xub)
               continue;
            else
               gap = (xval-xlb)*(xub-xval)/(1+2*ABS(xval));
            assert(!SCIPisNegative(scip, gap));
#ifdef WITH_CONSBRANCHNL
            SCIP_CALL( SCIPconshdlrBranchNonlinearUpdateVarInfeasibility(scip, conshdlrdata->branchnl, consdata->quadvar[j], MAX(gap, 0.)) );
#else
            SCIP_CALL( updateVarInfeasibility(scip, conshdlr, consdata->quadvar[j], MAX(gap, 0.)) );
#endif
            ++*nnotify;
         }
      }

      for (j = 0; j < consdata->n_bilin; ++j)
      { /* bilinear terms */
         xlb  = SCIPvarGetLbLocal(consdata->bilinvar1[j]);
         xub  = SCIPvarGetUbLocal(consdata->bilinvar1[j]);
         if (SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub))
         {
#ifdef WITH_CONSBRANCHNL
            SCIP_CALL(  SCIPconshdlrBranchNonlinearUpdateVarInfeasibility(scip, conshdlrdata->branchnl, consdata->bilinvar1[j], SCIPinfinity(scip)) );
#else
            SCIP_CALL( updateVarInfeasibility(scip, conshdlr, consdata->bilinvar1[j], SCIPinfinity(scip)) );
#endif
            ++*nnotify;
            continue;
         }

         ylb  = SCIPvarGetLbLocal(consdata->bilinvar2[j]);
         yub  = SCIPvarGetUbLocal(consdata->bilinvar2[j]);
         if (SCIPisInfinity(scip, -ylb) || SCIPisInfinity(scip, yub))
         {
#ifdef WITH_CONSBRANCHNL
            SCIP_CALL(  SCIPconshdlrBranchNonlinearUpdateVarInfeasibility(scip, conshdlrdata->branchnl, consdata->bilinvar2[j], SCIPinfinity(scip)) );
#else
            SCIP_CALL(  updateVarInfeasibility(scip, conshdlr, consdata->bilinvar2[j], SCIPinfinity(scip)) );
#endif
            ++*nnotify;
            continue;
         }

         xval = SCIPgetSolVal(scip, NULL, consdata->bilinvar1[j]);
         if (xval < xlb)
            xval = xlb;
         else if (xval > xub)
            xval = xub;
         
         yval = SCIPgetSolVal(scip, NULL, consdata->bilinvar2[j]);
         if (yval < ylb)
            yval = ylb;
         else if (yval > yub)
            yval = yub;
         
         coeff_ = SCIPisFeasPositive(scip, consdata->lhsviol) ? -consdata->bilincoeff[j] : consdata->bilincoeff[j];
         if (coeff_ > 0)
         {
            if ((xub-xlb)*yval + (yub-ylb)*xval <= xub*yub - xlb*ylb)
            {
               gap = (xval*yval - xlb*yval - ylb*xval + xlb*ylb) / (1+sqrt(xval*xval + yval*yval));
            }
            else
            {
               gap = (xval*yval - xval*yub - yval*xub + xub*yub) / (1+sqrt(xval*xval + yval*yval));
            }
         }
         else
         { /* coeff_ < 0 */
            if ((xub-xlb)*yval - (yub-ylb)*xval <= xub*ylb - xlb*yub)
            {
               gap = -(xval*yval - xval*ylb - yval*xub + xub*ylb) / (1+sqrt(xval*xval + yval*yval));
            }
            else
            {
               gap = -(xval*yval - xval*yub - yval*xlb + xlb*yub) / (1+sqrt(xval*xval + yval*yval));
            }
         }
         
         assert(!SCIPisNegative(scip, gap));
         if (gap < 0.) gap = 0.;
         
         if (!SCIPisEQ(scip, xlb, xub))
         {
#ifdef WITH_CONSBRANCHNL
            SCIP_CALL( SCIPconshdlrBranchNonlinearUpdateVarInfeasibility(scip, conshdlrdata->branchnl, consdata->bilinvar1[j], gap) );
#else
            SCIP_CALL( updateVarInfeasibility(scip, conshdlr, consdata->bilinvar1[j], gap) );
#endif            
            ++*nnotify;
         }
         if (!SCIPisEQ(scip, ylb, yub))
         {
#ifdef WITH_CONSBRANCHNL
            SCIP_CALL( SCIPconshdlrBranchNonlinearUpdateVarInfeasibility(scip, conshdlrdata->branchnl, consdata->bilinvar2[j], gap) );
#else
            SCIP_CALL( updateVarInfeasibility(scip, conshdlr, consdata->bilinvar2[j], gap) );
#endif
            ++*nnotify;
         }
      }
   }

   return SCIP_OKAY;
}

/** registers a variable from a violated constraint as branching candidate that has a large absolute value in the LP relaxation */
static
SCIP_RETCODE registerLargeLPValueVariableForBranching(
   SCIP*          scip,       /**< SCIP data structure */
   SCIP_CONSHDLR* conshdlr,   /**< constraint handler */
   SCIP_CONS**    conss,      /**< constraints */
   int            nconss,     /**< number of constraints */
   SCIP_VAR**     brvar       /**< buffer to store branching variable */
   )
{
   SCIP_CONSDATA*      consdata;
   SCIP_Real           val, brvarval = 0.0;
   int                 i, c;
   
   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   
   *brvar = NULL;
   
   for (c = 0; c < nconss; ++c)
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      
      if (!SCIPisFeasPositive(scip, consdata->lhsviol) && !SCIPisFeasPositive(scip, consdata->rhsviol))
         continue;
      
      for (i = 0; i < consdata->n_quadvar; ++i)
      {
         val = SCIPgetSolVal(scip, NULL, consdata->quadvar[i]);
         if (ABS(val) > brvarval)
         {
            brvarval = ABS(val);
            *brvar = consdata->quadvar[i];
         }
      }
   }
   
   if (*brvar)
   {
#ifdef WITH_CONSBRANCHNL
      assert(SCIPconshdlrGetData(conshdlr) != NULL);
      assert(SCIPconshdlrGetData(conshdlr)->branchnl != NULL);
      SCIP_CALL(  SCIPconshdlrBranchNonlinearUpdateVarInfeasibility(scip, SCIPconshdlrGetData(conshdlr)->branchnl, *brvar, brvarval) );
#else
      SCIP_CALL(  updateVarInfeasibility(scip, conshdlr, *brvar, brvarval) );
#endif
   }
   
   return SCIP_OKAY;
}

/** Solves a linear equation b*x in rhs and reduces bounds on x or deduces infeasibility if possible.
 */
static
SCIP_RETCODE propagateBoundsLinearVar(
   SCIP*          scip,          /**< SCIP data structure */ 
   SCIP_CONS*     cons,          /**< constraint where we currently propagate */
   SCIP_Real      intervalinfty, /**< infinity value used in interval operations */
   SCIP_VAR*      var,           /**< variable which domain we might reduce */
   SCIP_Real      b,             /**< linear coefficient of variable */
   SCIP_INTERVAL  rhs,           /**< right hand side */
   SCIP_RESULT*   result,        /**< result of propagation */
   int*           nchgbds        /**< buffer where to add number of tightened bounds */
   )
{
   SCIP_Bool infeas;
   SCIP_Bool tightened;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(var  != NULL);
   assert(result != NULL);
   assert(nchgbds != NULL);

   SCIPintervalDivScalar(intervalinfty, &rhs, rhs, b);
   
   if (SCIPisInfinity(scip, SCIPintervalGetInf(rhs)) || SCIPisInfinity(scip, -SCIPintervalGetSup(rhs)))
   { /* domain outside [-infty, +infty] -> declare node infeasible */
      *result = SCIP_CUTOFF;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      return SCIP_OKAY;
   }

   if (!SCIPisInfinity(scip, -SCIPintervalGetInf(rhs)))
   {
      SCIP_CALL( SCIPtightenVarLb(scip, var, SCIPintervalGetInf(rhs), FALSE, &infeas, &tightened) );
      if (infeas)
      {
         SCIPdebugMessage("found %s infeasible due to domain propagation for linear variable %s\n", SCIPconsGetName(cons), SCIPvarGetName(var));
         *result = SCIP_CUTOFF;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         return SCIP_OKAY;
      }
      if (tightened)
      {
         SCIPdebugMessage("tightened lower bound of linear variable %s in constraint %s to %g\n", SCIPvarGetName(var), SCIPconsGetName(cons), SCIPvarGetLbLocal(var));
         ++nchgbds;
         *result = SCIP_REDUCEDDOM;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
   }

   if (!SCIPisInfinity(scip, SCIPintervalGetSup(rhs)))
   {
      SCIP_CALL( SCIPtightenVarUb(scip, var, SCIPintervalGetSup(rhs), FALSE, &infeas, &tightened) );
      if (infeas)
      {
         SCIPdebugMessage("found %s infeasible due to domain propagation for linear variable %s\n", SCIPconsGetName(cons), SCIPvarGetName(var));
         *result = SCIP_CUTOFF;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         return SCIP_OKAY;
      }
      if (tightened)
      {
         SCIPdebugMessage("tightened upper bound of linear variable %s in constraint %s to %g\n", SCIPvarGetName(var), SCIPconsGetName(cons), SCIPvarGetUbLocal(var));
         ++nchgbds;
         *result = SCIP_REDUCEDDOM;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
   }

   return SCIP_OKAY;
}

/** Solves a quadratic equation a*x^2 + b*x in rhs (with b an interval) and reduces bounds on x or deduces infeasibility if possible.
 */
static
SCIP_RETCODE propagateBoundsQuadVar(
   SCIP*          scip,          /**< SCIP data structure */ 
   SCIP_CONS*     cons,          /**< constraint where we currently propagate */
   SCIP_Real      intervalinfty, /**< infinity value used in interval operations */
   SCIP_VAR*      var,           /**< variable which bounds with might tighten */
   SCIP_Real      a,             /**< coefficient in square term */
   SCIP_INTERVAL  b,             /**< coefficient in linear term */
   SCIP_INTERVAL  rhs,           /**< right hand side of quadratic equation */
   SCIP_RESULT*   result,        /**< result of propagation */
   int*           nchgbds        /**< buffer where to add number of tightened bounds */
   )
{
   SCIP_INTERVAL newrange;
   SCIP_Bool     infeas;
   SCIP_Bool     tightened;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(var  != NULL);
   assert(result != NULL);
   assert(nchgbds != NULL);

   if (SCIPvarGetLbLocal(var) >= 0.0)
   { /* need only positive solutions */
      SCIPintervalSolveUnivariateQuadExpressionPositive(intervalinfty, &newrange, a, b, rhs);
   }
   else if (SCIPvarGetUbLocal(var) <= 0.)
   { /* need only negative solutions */
      SCIP_INTERVAL tmp;
      SCIPintervalSetBounds(&tmp, -SCIPintervalGetSup(b), -SCIPintervalGetInf(b));
      SCIPintervalSolveUnivariateQuadExpressionPositive(intervalinfty, &tmp, a, tmp, rhs);
      if (SCIPintervalIsEmpty(tmp))
      {
         SCIPdebugMessage("found %s infeasible due to domain propagation for quadratic variable %s\n", SCIPconsGetName(cons), SCIPvarGetName(var));
         *result = SCIP_CUTOFF;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         return SCIP_OKAY;
      }
      SCIPintervalSetBounds(&newrange, -SCIPintervalGetSup(tmp), -SCIPintervalGetInf(tmp));
   }
   else
   {
      SCIPintervalSolveUnivariateQuadExpression(intervalinfty, &newrange, a, b, rhs);
   }

   if (SCIPisInfinity(scip, SCIPintervalGetInf(newrange)) || SCIPisInfinity(scip, -SCIPintervalGetSup(newrange)))
   { /* domain outside [-infty, +infty] -> declare node infeasible */
      SCIPdebugMessage("found %s infeasible because propagated domain of quadratic variable %s is outside of (-infty, +infty)\n", SCIPconsGetName(cons), SCIPvarGetName(var));
      *result = SCIP_CUTOFF;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      return SCIP_OKAY;
   }

   if (SCIPintervalIsEmpty(newrange))
   {
      SCIPdebugMessage("found %s infeasible due to domain propagation for quadratic variable %s\n", SCIPconsGetName(cons), SCIPvarGetName(var));
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   if (!SCIPisInfinity(scip, -SCIPintervalGetInf(newrange)))
   {
      SCIP_CALL( SCIPtightenVarLb(scip, var, SCIPintervalGetInf(newrange), FALSE, &infeas, &tightened) );
      if (infeas)
      {
         SCIPdebugMessage("found %s infeasible due to domain propagation for quadratic variable %s\n", SCIPconsGetName(cons), SCIPvarGetName(var));
         *result = SCIP_CUTOFF;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         return SCIP_OKAY;
      }
      if (tightened)
      {
         SCIPdebugMessage("tightened lower bound of quadratic variable %s in constraint %s to %g\n", SCIPvarGetName(var), SCIPconsGetName(cons), SCIPvarGetLbLocal(var));
         ++nchgbds;
         *result = SCIP_REDUCEDDOM;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
   }

   if (!SCIPisInfinity(scip, SCIPintervalGetSup(newrange)))
   {
      SCIP_CALL( SCIPtightenVarUb(scip, var, SCIPintervalGetSup(newrange), FALSE, &infeas, &tightened) );
      if (infeas)
      {
         SCIPdebugMessage("found %s infeasible due to domain propagation for quadratic variable %s\n", SCIPconsGetName(cons), SCIPvarGetName(var));
         *result = SCIP_CUTOFF;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         return SCIP_OKAY;
      }
      if (tightened)
      {
         SCIPdebugMessage("tightened upper bound of quadratic variable %s in constraint %s to %g -> %g\n", SCIPvarGetName(var), SCIPconsGetName(cons), SCIPintervalGetSup(newrange), SCIPvarGetUbLocal(var));
         ++nchgbds;
         *result = SCIP_REDUCEDDOM;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
   }

   return SCIP_OKAY;
}

/** updates the ranges for linear variables in a constraint data;
 * adds up all ranges;
 * returns whether there is a variable which range is entire
 */
static
void propagateBoundsUpdateLinRange(
   SCIP*             scip,              /**< SCIP data structure */
   SCIP_CONSDATA*    consdata,          /**< constraint data */
   SCIP_Real         intervalinfty,     /**< infinity value used in interval operations */
   SCIP_INTERVAL*    linrangesum,       /**< for summing up ranges of linear terms */
   int*              entire_var_index   /**< buffer to store index of single variable which domain is entire, or -1 if there is none, or -2 if there are at least two */ 
   )
{
   SCIP_VAR* var;
   int       i;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(linrangesum != NULL);
   assert(entire_var_index != NULL);

   SCIPintervalSet(linrangesum,  0.);
   *entire_var_index = -1;
   for (i = 0; i < consdata->n_linvar; ++i)
   {
      if (SCIPintervalIsEmpty(consdata->linrange[i]))
      {
         var = consdata->linvar[i];
         SCIPintervalSetBounds(&consdata->linrange[i],
            -infty2infty(SCIPinfinity(scip), intervalinfty, -SCIPvarGetLbLocal(var)),
             infty2infty(SCIPinfinity(scip), intervalinfty,  SCIPvarGetUbLocal(var))
            );
         SCIPintervalMulScalar(intervalinfty, &consdata->linrange[i], consdata->linrange[i], consdata->lincoeff[i]);
      }
#ifndef NDEBUG
      else {
         SCIP_INTERVAL tmp;
         var = consdata->linvar[i];
         SCIPintervalSetBounds(&tmp,
            -infty2infty(SCIPinfinity(scip), intervalinfty, -SCIPvarGetLbLocal(var)),
             infty2infty(SCIPinfinity(scip), intervalinfty,  SCIPvarGetUbLocal(var))
            );
         SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->lincoeff[i]);
         assert(SCIPintervalIsSubsetEQ(intervalinfty, tmp, consdata->linrange[i]));
      }
#endif

      if (SCIPintervalIsEntire(intervalinfty, consdata->linrange[i]))
      {
         if (*entire_var_index >= 0)
         {
            *entire_var_index = -2;
            return;
         }
         *entire_var_index = i;
      }
      
      SCIPintervalAdd(intervalinfty, linrangesum, *linrangesum, consdata->linrange[i]);
   }
}

/** updates the ranges for quadratic terms associated to each variable in a constraint data;
 * adds up all ranges;
 * returns whether there is a term which range is entire
 */
static
void propagateBoundsUpdateQuadRangeVar(
   SCIP*             scip,              /**< SCIP data structure */ 
   SCIP_CONSDATA*    consdata,          /**< constraint data */
   SCIP_Real         intervalinfty,     /**< infinity value used in interval operations */
   SCIP_INTERVAL*    quadrangesum,      /**< for summing up ranges of quadratic terms */
   int*              entire_var_index   /**< buffer to store index of single variable which domain is entire, or -1 if there is none, or -2 if there are at least two */ 
   )
{
   SCIP_VAR*     var;
   SCIP_INTERVAL tmp;
   int           i;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(quadrangesum != NULL);
   assert(entire_var_index != NULL);

   SCIPintervalSet(quadrangesum, 0.0);
   *entire_var_index = -1; /* index of quadratic variable which range is entire */
   for (i = 0; i < consdata->n_quadvar; ++i)
   {
      if (SCIPintervalIsEmpty(consdata->quadrangevar[i]))
      {
         var = consdata->quadvar[i];
         SCIPintervalSetBounds(&consdata->quadrangevar[i], 
            -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var))), 
             infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)))
             );

         if (consdata->quadlincoeff[i])
         {
            SCIPintervalSet(&tmp, consdata->quadlincoeff[i]);
            SCIPintervalQuad(intervalinfty, &consdata->quadrangevar[i], consdata->quadsqrcoeff[i], tmp, consdata->quadrangevar[i]);
         }
         else
         {
            SCIPintervalSquare(intervalinfty, &consdata->quadrangevar[i], consdata->quadrangevar[i]);
            assert(SCIPintervalGetInf(consdata->quadrangevar[i]) < intervalinfty);
            SCIPintervalMulScalar(intervalinfty, &consdata->quadrangevar[i], consdata->quadrangevar[i], consdata->quadsqrcoeff[i]);
         }
      }
#ifndef NDEBUG
      else
      {
         var = consdata->quadvar[i];
         SCIPintervalSetBounds(&tmp, 
            -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var))), 
             infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)))
             );

         if (consdata->quadlincoeff[i])
         {
            SCIP_INTERVAL tmp2;
            SCIPintervalSet(&tmp2, consdata->quadlincoeff[i]);
            SCIPintervalQuad(intervalinfty, &tmp, consdata->quadsqrcoeff[i], tmp2, tmp);
         }
         else
         {
            SCIPintervalSquare(intervalinfty, &tmp, tmp);
            SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->quadsqrcoeff[i]);
         }
         assert(SCIPintervalIsSubsetEQ(intervalinfty, tmp, consdata->quadrangevar[i]));
      }
#endif

      if (SCIPintervalIsEntire(intervalinfty, consdata->quadrangevar[i]))
      {
         if (*entire_var_index >= 0)
         {
            *entire_var_index = -2;
            return; /* cannot reduce bounds on any variable if there are more than one variable with an entire range */
         }
         *entire_var_index = i;
      }
      
      SCIPintervalAdd(intervalinfty, quadrangesum, *quadrangesum, consdata->quadrangevar[i]);
   }
}

/** updates the ranges for bilinear terms in a constraint data;
 * adds up all ranges
 */
static
void propagateBoundsUpdateBilinRange(
   SCIP*             scip,              /**< SCIP data structure */ 
   SCIP_CONSDATA*    consdata,          /**< constraint data */
   SCIP_Real         intervalinfty,     /**< infinity value used in interval operations */
   SCIP_INTERVAL*    bilinrangesum      /**< for summing up ranges of bilinear terms */
   )
{
   SCIP_VAR*      var;
   SCIP_INTERVAL  tmp;
   int            i;

   SCIPintervalSet(bilinrangesum, 0.);
   for (i = 0; i < consdata->n_bilin; ++i)
   { /* check if a bilinrange need an update */
      if (SCIPintervalIsEmpty(consdata->bilinrange[i]))
      {
         var = consdata->bilinvar1[i];
         SCIPintervalSetBounds(&consdata->bilinrange[i],
            -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var))), 
             infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)))
            );
         var = consdata->bilinvar2[i];
         SCIPintervalSetBounds(&tmp,
            -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var))),
             infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var)))
            );
         SCIPintervalMul(intervalinfty, &consdata->bilinrange[i], consdata->bilinrange[i], tmp);
         assert(SCIPintervalGetInf(consdata->bilinrange[i]) <  intervalinfty);
         assert(SCIPintervalGetSup(consdata->bilinrange[i]) > -intervalinfty);
         SCIPintervalMulScalar(intervalinfty, &consdata->bilinrange[i], consdata->bilinrange[i], consdata->bilincoeff[i]);
      }
#ifndef NDEBUG
      else 
      {
         SCIP_INTERVAL tmp2;
         var = consdata->bilinvar1[i];
         SCIPintervalSetBounds(&tmp2, 
            -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var))),
             infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var)))
            );
         var = consdata->bilinvar2[i];
         SCIPintervalSetBounds(&tmp,
            -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var))),
             infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var)))
            );
         SCIPintervalMul(intervalinfty, &tmp, tmp2, tmp);
         SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilincoeff[i]);
         assert(SCIPintervalIsSubsetEQ(intervalinfty, tmp, consdata->bilinrange[i]));
      }
#endif
      SCIPintervalAdd(intervalinfty, bilinrangesum, *bilinrangesum, consdata->bilinrange[i]);
   }
}

/** Propagates bounds on a constraint */
static
SCIP_RETCODE propagateBounds(
   SCIP*              scip,     /**< SCIP data structure */
   SCIP_CONSHDLR*     conshdlr, /**< constraint handler */
   SCIP_CONS*         cons,     /**< constraint to process */
   SCIP_RESULT*       result,   /**< pointer to store the result of the propagation call */
   int*               nchgbds   /**< buffer where to add the the number of changed bounds */
  )
{
   SCIP_CONSDATA*     consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_INTERVAL      consbounds;    /* lower and upper bounds of constraint */
   SCIP_Real          intervalinfty; /* infinity used for interval computation */  
   int                i, j, k, l;
   int                entire_var_index;  /* index of a variable which domain is entire */
   SCIP_INTERVAL      linrangesum;   /* range of linear part */
   SCIP_INTERVAL      quadrangesum;  /* sum of ranges of quadratic variable parts (ax^2+bx) */
   SCIP_INTERVAL      bilinrangesum; /* range of complete bilinear part */
   SCIP_VAR*          var;
   SCIP_Real          a;   /* quadratic coefficient of quadratic equation */
   SCIP_INTERVAL      b;   /* linear coefficient of quadratic equation */ 
   SCIP_INTERVAL      rhs; /* right hand side of quadratic equation */
   SCIP_INTERVAL      tmp;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(cons     != NULL);
   assert(result   != NULL);
   assert(nchgbds  != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->linrange     || consdata->n_linvar  == 0);
   assert(consdata->quadrangevar || consdata->n_quadvar == 0);
   assert(consdata->bilinrange   || consdata->n_bilin   == 0);

   *result = SCIP_DIDNOTRUN;

   if (consdata->is_propagated)
      return SCIP_OKAY;

#ifndef NDEBUG
   /* assert that there are no variables that are fixed to -/+ infinity */
   for (i = 0; i < consdata->n_linvar; ++i)
   {
      assert(!SCIPisInfinity(scip,  SCIPvarGetLbLocal(consdata->linvar[i])));
      assert(!SCIPisInfinity(scip, -SCIPvarGetUbLocal(consdata->linvar[i])));
   }

   for (i = 0; i < consdata->n_quadvar; ++i)
   {
      assert(!SCIPisInfinity(scip,  SCIPvarGetLbLocal(consdata->quadvar[i])));
      assert(!SCIPisInfinity(scip, -SCIPvarGetUbLocal(consdata->quadvar[i])));
   }
#endif

   consdata->is_propagated = TRUE;

   *result = SCIP_DIDNOTFIND;
   intervalinfty = 1000 * SCIPinfinity(scip) * SCIPinfinity(scip);

   SCIPintervalSetBounds(&consbounds,
      -infty2infty(SCIPinfinity(scip), intervalinfty, -consdata->lhs),
       infty2infty(SCIPinfinity(scip), intervalinfty,  consdata->rhs)
      );
   
   propagateBoundsUpdateLinRange(scip, consdata, intervalinfty, &linrangesum, &entire_var_index);
   if (entire_var_index == -2)
      return SCIP_OKAY; /* at least two variables that are completely unbounded; cannot propagate anything */

   if (SCIPintervalIsEmpty(consdata->quadrange))
   { /* quadrange needs update */
      SCIPintervalSet(&consdata->quadrange, 0.0);
      SCIP_CALL( addQuadRange(scip, cons, intervalinfty, &consdata->quadrange, NULL) );
   }

   if (entire_var_index >= 0)
   {
      assert(entire_var_index < consdata->n_linvar);
      if (SCIPintervalIsEntire(intervalinfty, consdata->quadrange))
         return SCIP_OKAY; /* nothing we can do */

      var = consdata->linvar[entire_var_index];
      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR);

      rhs = consbounds;
      for (i = 0; i < consdata->n_linvar; ++i)
         if (i != entire_var_index)
            SCIPintervalSub(intervalinfty, &rhs, rhs, consdata->linrange[i]);
      SCIPintervalSub(intervalinfty, &rhs, rhs, consdata->quadrange);

      SCIP_CALL( propagateBoundsLinearVar(scip, cons, intervalinfty, var, consdata->lincoeff[entire_var_index], rhs, result, nchgbds) );

      return SCIP_OKAY;
   }
   if (SCIPintervalIsEntire(intervalinfty, linrangesum)) /* can still happen if two half-unbounded were added -> bad luck */
      return SCIP_OKAY;

   /* intersects linrangesum+quadrange with consbounds */
   SCIPintervalAdd(intervalinfty, &tmp, linrangesum, consdata->quadrange);
#if 1 /* hopefully ok now where bug in interval arith. fixed */
   SCIPintervalIntersect(&consbounds, consbounds, tmp);
   if (SCIPintervalIsEmpty(consbounds))
   {
      SCIPdebugMessage("found %s infeasible due to forward propagation\n", SCIPconsGetName(cons));
      *result = SCIP_CUTOFF;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      return SCIP_OKAY;
   }
#else
   SCIPintervalIntersect(&tmp, consbounds, tmp);
   if (SCIPintervalIsEmpty(tmp))
   { /* check again with slightly larger bounds: workaround for instance product */
      SCIP_INTERVAL tmp2;
      SCIPintervalSetBounds(&tmp2, -SCIPfeastol(scip)/2, SCIPfeastol(scip)/2);
      SCIPintervalAdd(intervalinfty, &tmp, linrangesum, consdata->quadrange);
      SCIPintervalAdd(intervalinfty, &tmp, tmp, tmp2);
      SCIPintervalIntersect(&consbounds, consbounds, tmp);
      if (SCIPintervalIsEmpty(consbounds))
      {
         SCIPdebugMessage("found %s infeasible due to forward propagation\n", SCIPconsGetName(cons));
         *result = SCIP_CUTOFF;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         return SCIP_OKAY;
      }
      else
      {
         SCIPdebugMessage("node only feasible after widening bounds by %g\n", SCIPfeastol(scip));
      }
   }
   else
      consbounds = tmp;
#endif
   
   /* domain propagation for linear variables */
   for (i = 0; i < consdata->n_linvar; ++i)
   {
      var = consdata->linvar[i];
      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR);
      assert(!SCIPisZero(scip, consdata->lincoeff[i]));

      SCIPintervalSub(intervalinfty, &rhs, consbounds, linrangesum);
      SCIPintervalUndoSub(intervalinfty, &rhs, rhs, consdata->linrange[i]);
      SCIPintervalSub(intervalinfty, &rhs, rhs, consdata->quadrange);

      SCIP_CALL( propagateBoundsLinearVar(scip, cons, intervalinfty, var, consdata->lincoeff[i], rhs, result, nchgbds) );
      if (*result == SCIP_CUTOFF)
         return SCIP_OKAY;
   }

   if (!consdata->n_quadvar)
      return SCIP_OKAY;

   propagateBoundsUpdateQuadRangeVar(scip, consdata, intervalinfty, &quadrangesum, &entire_var_index);
   if (entire_var_index == -2)
      return SCIP_OKAY; /* cannot reduce bounds on any variable if there are more than one variable which domain is entire */

   if (entire_var_index >= 0)
   { /* there is exactly one quadratic variable which domain is entire; this one is the only chance where we can achieve a bound tightening then */ 
      var = consdata->quadvar[entire_var_index];
      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR);

      SCIPintervalSub(intervalinfty, &rhs, consbounds, linrangesum);
      
      SCIPintervalSet(&tmp, 0.0);
      SCIP_CALL( addQuadRange(scip, cons, intervalinfty, &tmp, var) );
      SCIPintervalSub(intervalinfty, &rhs, rhs, tmp);

      /* add up coefficient interval in linear term of expression */
      SCIPintervalSet(&b, consdata->quadlincoeff[entire_var_index]);
      for (k = 0; k < consdata->n_adjbilin[entire_var_index]; ++k)
      {
         l = consdata->adjbilin[entire_var_index][k];
         if (consdata->bilinvar1[l] == var)
         {
            assert(consdata->bilinvar2[l] != var);
            SCIPintervalSetBounds(&tmp, 
               -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinvar2[l]), SCIPvarGetUbLocal(consdata->bilinvar2[l]))),
                infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinvar2[l]), SCIPvarGetUbLocal(consdata->bilinvar2[l])))
               );
         }
         else
         {
            assert(consdata->bilinvar2[l] == var);
            SCIPintervalSetBounds(&tmp, 
               -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinvar1[l]), SCIPvarGetUbLocal(consdata->bilinvar1[l]))),
                infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinvar1[l]), SCIPvarGetUbLocal(consdata->bilinvar1[l])))
               );
         }
         SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilincoeff[l]);
         SCIPintervalAdd(intervalinfty, &b, b, tmp);
      }
      if (SCIPintervalIsEntire(SCIPinfinity(scip), b))
         return SCIP_OKAY; /* no hope to reduce a bound */

      a = consdata->quadsqrcoeff[entire_var_index];

      SCIP_CALL( propagateBoundsQuadVar(scip, cons, intervalinfty, var, a, b, rhs, result, nchgbds) );    
      return SCIP_OKAY;
   }

   propagateBoundsUpdateBilinRange(scip, consdata, intervalinfty, &bilinrangesum);
   if (SCIPintervalIsEntire(intervalinfty, bilinrangesum)) /* propagation on quad. vars makes no sense */
      return SCIP_OKAY;

   /* move everything into consbounds */
   SCIPintervalSub(intervalinfty, &consbounds, consbounds, linrangesum);
   if (conshdlrdata->fast_propagate && consdata->n_quadvar > 2)
   {
      SCIPintervalSub(intervalinfty, &consbounds, consbounds, bilinrangesum);
      SCIPintervalSub(intervalinfty, &consbounds, consbounds, quadrangesum);
   }

   if (SCIPintervalIsEntire(intervalinfty, consbounds))
      return SCIP_OKAY;

   for (j = 0; j < consdata->n_quadvar; ++j)
   {
      var = consdata->quadvar[j];
      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR);

      /* skip fixed variables */
      if (SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)))
         continue;

      /* setup equation a*x^2 + b*x \in rhs, x is var */
      a = consdata->quadsqrcoeff[j];

      SCIPintervalSet(&b, consdata->quadlincoeff[j]);
      if (conshdlrdata->fast_propagate)
      {
         if (consdata->n_quadvar > 2)
            SCIPintervalUndoSub(intervalinfty, &rhs, consbounds, consdata->quadrangevar[j]);
         else if (consdata->n_quadvar == 2)
         {
            if (SCIPintervalIsEmpty(consdata->quadrangevar[1-j]))
            { /* this can happen if j==1 and we just improved the bound for j==0 */
               SCIPintervalSetBounds(&tmp,
                  -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->quadvar[1-j]), SCIPvarGetUbLocal(consdata->quadvar[1-j]))),
                   infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->quadvar[1-j]), SCIPvarGetUbLocal(consdata->quadvar[1-j])))
                  );
               if (consdata->quadlincoeff[1-j])
               {
                  SCIP_INTERVAL tmp2;
                  SCIPintervalSet(&tmp2, consdata->quadlincoeff[1-j]);
                  SCIPintervalQuad(intervalinfty, &tmp, consdata->quadsqrcoeff[1-j], tmp2, tmp);
               }
               else
               {
                  SCIPintervalSquare(intervalinfty, &tmp, tmp);
                  assert(SCIPintervalGetInf(tmp) < intervalinfty);
                  SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->quadsqrcoeff[1-j]);
               }
               SCIPintervalSub(intervalinfty, &rhs, consbounds, tmp);
            }
            else
            {
               SCIPintervalSub(intervalinfty, &rhs, consbounds, consdata->quadrangevar[1-j]);
            }
         }
         else
            rhs = consbounds;

         /* we should not just put all bilinear terms into the right hand side, that would be fatal for equations like x*y \in [...]
          * since recomputing the best quad range is too expensive,
          * we undo all the substractions of bilinrange where this var is involved and setup an appropriate linear term (b)
          */
         if (consdata->n_adjbilin[j])
         {
            for (k = 0; k < consdata->n_adjbilin[j]; ++k)
            {
               l = consdata->adjbilin[j][k];
               if (consdata->n_quadvar > 2)
               {
                  if (SCIPintervalIsEmpty(consdata->bilinrange[l]))
                  { /* this might happen, if we just reduced a bound on a variable in this bilinear term
                     * however, we should not update bilinrange[l] here since it is not invalidated in future bound changes as long as the corresponding quadrangevar's are invalid
                     */
                     SCIP_INTERVAL tmp2;
                     SCIPintervalSetBounds(&tmp2,
                        -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinvar1[l]), SCIPvarGetUbLocal(consdata->bilinvar1[l]))),
                         infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinvar1[l]), SCIPvarGetUbLocal(consdata->bilinvar1[l])))
                        );
                     SCIPintervalSetBounds(&tmp,
                        -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinvar2[l]), SCIPvarGetUbLocal(consdata->bilinvar2[l]))),
                         infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinvar2[l]), SCIPvarGetUbLocal(consdata->bilinvar2[l])))
                        );
                     SCIPintervalMul(intervalinfty, &tmp, tmp, tmp2);
                     assert(SCIPintervalGetInf(tmp) <  intervalinfty);
                     assert(SCIPintervalGetSup(tmp) > -intervalinfty);
                     SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilincoeff[l]);
                     SCIPintervalUndoSub(intervalinfty, &rhs, rhs, tmp);
                  }
                  else
                     SCIPintervalUndoSub(intervalinfty, &rhs, rhs, consdata->bilinrange[l]);
               }

               if (consdata->bilinvar1[l] == var)
               {
                  assert(consdata->bilinvar2[l] != var);
                  SCIPintervalSetBounds(&tmp, 
                     -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinvar2[l]), SCIPvarGetUbLocal(consdata->bilinvar2[l]))),
                      infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinvar2[l]), SCIPvarGetUbLocal(consdata->bilinvar2[l])))
                     );
               }
               else
               {
                  assert(consdata->bilinvar2[l] == var);
                  SCIPintervalSetBounds(&tmp, 
                     -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinvar1[l]), SCIPvarGetUbLocal(consdata->bilinvar1[l]))),
                      infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinvar1[l]), SCIPvarGetUbLocal(consdata->bilinvar1[l])))
                     );
               }
               SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilincoeff[l]);
               SCIPintervalAdd(intervalinfty, &b, b, tmp);
            }
            if (SCIPintervalIsEntire(intervalinfty, b))
               continue; /* no hope to reduce a bound */
         }
      }
      else
      {
         SCIPintervalSet(&tmp, 0.);
         SCIP_CALL( addQuadRange(scip, cons, intervalinfty, &tmp, var) ); /* add up everything in quad.part not belonging to var */
         SCIPintervalSub(intervalinfty, &rhs, consbounds, tmp);   /* put bounds - linrangesum - quadratic_except_var into rhs */
         for (k = 0; k < consdata->n_adjbilin[j]; ++k)
         {
            l = consdata->adjbilin[j][k];
            if (consdata->bilinvar1[l] == var)
            {
               assert(consdata->bilinvar2[l] != var);
               SCIPintervalSetBounds(&tmp,
                  -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinvar2[l]), SCIPvarGetUbLocal(consdata->bilinvar2[l]))),
                   infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinvar2[l]), SCIPvarGetUbLocal(consdata->bilinvar2[l])))
                  );
            }
            else
            {
               assert(consdata->bilinvar2[l] == var);
               SCIPintervalSetBounds(&tmp, 
                  -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinvar1[l]), SCIPvarGetUbLocal(consdata->bilinvar1[l]))),
                   infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinvar1[l]), SCIPvarGetUbLocal(consdata->bilinvar1[l])))
                  );
            }
            SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilincoeff[l]);
            SCIPintervalAdd(intervalinfty, &b, b, tmp);
         }
         if (SCIPintervalIsEntire(intervalinfty, b))
            continue; /* no hope to reduce a bound */
      }

      if (SCIPintervalIsEntire(intervalinfty, rhs))
         continue;

      SCIP_CALL( propagateBoundsQuadVar(scip, cons, intervalinfty, var, a, b, rhs, result, nchgbds) );
      if (*result == SCIP_CUTOFF)
         break;
   }

   return SCIP_OKAY;
}

/** NLPI initialization method of constraint handler
 * 
 * The constraint handler should create an NLPI representation of the constraints in the provided NLPI.
 */
SCIP_RETCODE SCIPconsInitnlpiQuadratic(
   SCIP*           scip,     /**< SCIP data structure */
   SCIP_CONSHDLR*  conshdlr, /**< quadratic constraint handler - it's C wrapper*/
   SCIP_NLPI*      nlpi,     /**< NLPI where to add constraints */
   int             nconss,   /**< number of constraints */
   SCIP_CONS**     conss,    /**< quadratic constraints */
   SCIP_HASHMAP*   var_scip2nlp /**< mapping from SCIP variables to variable indices in NLPI */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real*     lhs;
   SCIP_Real*     rhs;
   int*           linoffset;   /* row offsets */
   int*           linindex;    /* column indices */
   SCIP_Real*     lincoeff;    /* coefficients of linear variables */
   int            linnnz;
   int*           nquadrows;
   int**          quadrowidx;
   int**          quadoffset;
   int**          quadindex;
   SCIP_Real**    quadcoeff;
   int            quadnnz;
   int            i, j, k, l;
   int            lincnt;      /* how many linear cofficients written so far */
   SCIP_HASHMAP*  var2rowidx;
   SCIP_VAR*      othervar;
      
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(nlpi != NULL);
   assert(conss != NULL || nconss == 0);

   if (!nconss)
      return SCIP_OKAY;

   assert(var_scip2nlp != NULL);

   linnnz = 0;
   for (i = 0; i < nconss; ++i)
   {
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      
      linnnz += consdata->n_linvar;
      for (j = 0; j < consdata->n_quadvar; ++j)
         if (consdata->quadlincoeff[j])
            ++linnnz; /* we assume here that the quadratic variables are disjunct from the linear variables; but it might not harm if this assumption does not hold */
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &lhs, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhs, nconss) );

   SCIP_CALL( SCIPallocBufferArray(scip, &linoffset, nconss+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linindex,  linnnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lincoeff,  linnnz) );

   SCIP_CALL( SCIPallocBufferArray(scip, &nquadrows,  nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadrowidx, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadoffset, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadindex,  nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadcoeff,  nconss) );

   lincnt = 0;
   for (i = 0; i < nconss; ++i)
   {
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      lhs[i] = consdata->lhs;
      rhs[i] = consdata->rhs;

      linoffset[i] = lincnt;
      memcpy(&lincoeff[lincnt], consdata->lincoeff, consdata->n_linvar*sizeof(SCIP_Real));
      for (j = 0; j < consdata->n_linvar; ++j, ++lincnt)
      {
         assert(SCIPhashmapExists(var_scip2nlp, consdata->linvar[j]));
         linindex[lincnt] = (int) (size_t) SCIPhashmapGetImage(var_scip2nlp, consdata->linvar[j]);
      }

      quadnnz = consdata->n_bilin;
      for (j = 0; j < consdata->n_quadvar; ++j)
         if (consdata->quadsqrcoeff[j])
            ++quadnnz;

      if (!quadnnz)
      {
         nquadrows[i]  = 0;
         quadrowidx[i] = NULL;
         quadoffset[i] = NULL;
         quadindex[i]  = NULL;
         quadcoeff[i]  = NULL;
         continue;
      }

      nquadrows[i] = consdata->n_quadvar;
      SCIPallocBufferArray(scip, &quadrowidx[i], consdata->n_quadvar);
      SCIPallocBufferArray(scip, &quadoffset[i], consdata->n_quadvar+1);
      SCIPallocBufferArray(scip, &quadindex[i],  quadnnz);
      SCIPallocBufferArray(scip, &quadcoeff[i],  quadnnz);

      SCIP_CALL( SCIPhashmapCreate(&var2rowidx, SCIPblkmem(scip), consdata->n_quadvar) );
      k = 0;
      for (j = 0; j < consdata->n_quadvar; ++j)
      {
         assert(SCIPhashmapExists(var_scip2nlp, consdata->quadvar[j]));
         if (consdata->quadlincoeff[j])
         {
            linindex[lincnt] = (int) (size_t) SCIPhashmapGetImage(var_scip2nlp, consdata->quadvar[j]);
            lincoeff[lincnt] = consdata->quadlincoeff[j];
            ++lincnt;
         }

         assert( !SCIPhashmapExists(var2rowidx, consdata->quadvar[j]) );
         quadrowidx[i][j] = (int) (size_t) SCIPhashmapGetImage(var_scip2nlp, consdata->quadvar[j]);
         SCIP_CALL( SCIPhashmapInsert(var2rowidx, consdata->quadvar[j], (void*) (size_t) j) );

         quadoffset[i][j] = k;
         if (consdata->quadsqrcoeff[j])
         {
            assert(k < quadnnz);
            quadindex[i][k] = j;
            quadcoeff[i][k] = consdata->quadsqrcoeff[j];
            ++k;
         }

         for (l = 0; l < consdata->n_adjbilin[j]; ++l)
         {
            othervar = consdata->bilinvar1[consdata->adjbilin[j][l]];
            if (othervar == consdata->quadvar[j])
               othervar = consdata->bilinvar2[consdata->adjbilin[j][l]];
            assert(othervar != consdata->quadvar[j]);

            if (SCIPhashmapExists(var2rowidx, othervar)) /* processed the other var already, so now its time to add the corresponding bilin term */
            {
               assert(k < quadnnz);
               quadindex[i][k] = (int) (size_t) SCIPhashmapGetImage(var2rowidx, othervar);
               quadcoeff[i][k] = consdata->bilincoeff[consdata->adjbilin[j][l]];
               ++k;
            }
            /* otherwise we leave this bilinear term for later */
         }
      }
      quadoffset[i][consdata->n_quadvar] = k;
      assert(k == quadnnz);
      SCIPhashmapFree(&var2rowidx);
   }
   linoffset[nconss] = lincnt;
   assert(lincnt == linnnz);

   SCIP_CALL( SCIPnlpiAddConstraints(scip, nlpi, nconss,
      lhs, rhs,
      linoffset, linindex, lincoeff,
      nquadrows, quadrowidx, quadoffset, quadindex, quadcoeff,
      NULL, NULL, NULL) );

   for (i = 0; i < nconss; ++i)
   {
      SCIPfreeBufferArrayNull(scip, &quadrowidx[i]);
      SCIPfreeBufferArrayNull(scip, &quadoffset[i]);
      SCIPfreeBufferArrayNull(scip, &quadindex[i]);
      SCIPfreeBufferArrayNull(scip, &quadcoeff[i]);
   }

   SCIPfreeBufferArray(scip, &lhs);
   SCIPfreeBufferArray(scip, &rhs);

   SCIPfreeBufferArray(scip, &linoffset);
   SCIPfreeBufferArray(scip, &linindex);
   SCIPfreeBufferArray(scip, &lincoeff);

   SCIPfreeBufferArray(scip, &nquadrows);
   SCIPfreeBufferArray(scip, &quadrowidx);
   SCIPfreeBufferArray(scip, &quadoffset);
   SCIPfreeBufferArray(scip, &quadindex);
   SCIPfreeBufferArray(scip, &quadcoeff);

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#if 1
static
SCIP_DECL_CONSFREE(consFreeQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   
   SCIPfreeMemory(scip, &conshdlrdata);
   
   return SCIP_OKAY;
}
#else
#define consFreeQuadratic NULL
#endif


/** initialization method of constraint handler (called after problem was transformed) */
#if 1
static
SCIP_DECL_CONSINIT(consInitQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

#ifdef WITH_CONSBRANCHNL
   conshdlrdata->branchnl = SCIPfindConshdlr(scip, "branchnonlinear");
   if (conshdlrdata->branchnl == NULL && nconss > 0)
   {
      SCIPerrorMessage("cannot find constraint handler for branching on nonlinear variables");
      return SCIP_PLUGINNOTFOUND;
   }
#else
   /* TODO: what is a good estimate for the hashmap size? should the constraint handler notify about the number of potential candidates? */
   SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->branchcand, SCIPblkmem(scip), MAX(1, SCIPgetNVars(scip))) );
#endif
   
#ifndef WITH_LAPACK
   if (nconss)
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Quadratic constraint handler does not have LAPACK for eigenvalue computation. Will assume that matrices (with size > 2x2) are indefinite.\n");
   }
#endif
   
   return SCIP_OKAY;
}
#else
#define consInitQuadratic NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 1
static
SCIP_DECL_CONSEXIT(consExitQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

#ifdef WITH_CONSBRANCHNL
   conshdlrdata->branchnl = NULL;
#else
   SCIP_CALL( clearBranchingCandidates(scip, conshdlr) );
   if (conshdlrdata->branchcand != NULL)
      SCIPhashmapFree(&conshdlrdata->branchcand);
#endif
   
   return SCIP_OKAY;
}
#else
#define consExitQuadratic NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 1
static
SCIP_DECL_CONSINITPRE(consInitpreQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   
   *result = SCIP_FEASIBLE;

   for (c = 0; c < nconss; ++c)
       SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );

   return SCIP_OKAY;
}
#else
#define consInitpreQuadratic NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 1
static
SCIP_DECL_CONSEXITPRE(consExitpreQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_HASHMAP*      terms = NULL;
   SCIP_Real          constant;
   SCIP_Bool          have_change;
   int                i, c;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_FEASIBLE;

   for (c = 0; c < nconss; ++c)
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if (!consdata->is_removedfixings)
      {
         SCIP_CALL( presolveCreateQuadTerm(scip, &terms, &constant, &have_change, conss[c], conshdlrdata->replace_sqrbinary) );

         if (have_change)
         {
            SCIPdebugMessage("exitpre found changed variable in constraint %s:\n", SCIPconsGetName(conss[c]));
            /* unlock all variables */
            for (i = 0; i < consdata->n_linvar; ++i)
               SCIP_CALL( unlockLinearVariable(scip, conss[c], consdata->linvar[i], consdata->lincoeff[i]) );
            for (i = 0; i < consdata->n_quadvar; ++i)
               SCIP_CALL( SCIPunlockVarCons(scip, consdata->quadvar[i], conss[c], TRUE, TRUE) );
            SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );

            /* reset function */
            SCIP_CALL( consdataSetFunctionData(scip, consdata, terms, -1) );
            if (!SCIPisInfinity(scip, -consdata->lhs))
               consdata->lhs -= constant;
            if (!SCIPisInfinity(scip,  consdata->rhs))
               consdata->rhs -= constant;

            /* lock all variables */
            for (i = 0; i < consdata->n_linvar; ++i)
               SCIP_CALL( lockLinearVariable(scip, conss[c], consdata->linvar[i], consdata->lincoeff[i]) );
            for (i = 0; i < consdata->n_quadvar; ++i)
               SCIP_CALL( SCIPlockVarCons(scip, consdata->quadvar[i], conss[c], TRUE, TRUE) );
            SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );
         }

         presolveQuadTermFree(scip, &terms);
      }

#ifndef NDEBUG
      for (i = 0; i < consdata->n_linvar; ++i)
         assert(SCIPvarIsActive(consdata->linvar[i]));

      for (i = 0; i < consdata->n_quadvar; ++i)
         assert(SCIPvarIsActive(consdata->quadvar[i]));
#endif

      SCIP_CALL( boundUnboundedVars(scip, conss[c], conshdlrdata->defaultbound, NULL) );
   }

   return SCIP_OKAY;
}
#else
#define consExitpreQuadratic NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 1
static
SCIP_DECL_CONSINITSOL(consInitsolQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Bool          convex = TRUE;
   int                c;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   
   for (c = 0; c < nconss; ++c)
   {
       consdata = SCIPconsGetData(conss[c]);
       assert(consdata != NULL);
       
       SCIP_CALL( checkCurvature(scip, conshdlr, conss[c]) );
       
       if (!SCIPisInfinity(scip,  consdata->rhs) && !consdata->is_convex)
       {
           SCIPdebugMessage("nonconvex because of upper bound in con %s\n", SCIPconsGetName(conss[c]));
           convex = FALSE;
       }
       if (!SCIPisInfinity(scip, -consdata->lhs) && !consdata->is_concave)
       {
           SCIPdebugMessage("nonconvex because of lower bound in con %s\n", SCIPconsGetName(conss[c]));
           convex = FALSE;
       }
   }
   
   if (nconss)
   {
       if (convex)
           SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "All quadratic constraints are convex\n");
       else
           SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "There are nonconvex quadratic constraints.\n");
   }
       
   conshdlrdata->nlpheur = SCIPfindHeur(scip, "nlp");
   
   return SCIP_OKAY;
}
#else
#define consInitsolQuadratic NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 1
static
SCIP_DECL_CONSEXITSOL(consExitsolQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int                c;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for (c = 0; c < nconss; ++c)
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );
   
   conshdlrdata->nlpheur = NULL;

   return SCIP_OKAY;
}
#else
#define consExitsolQuadratic NULL
#endif


/** frees specific constraint data */
#if 1
static
SCIP_DECL_CONSDELETE(consDeleteQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int                i;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if (SCIPgetStage(scip) == SCIP_STAGE_SOLVING)
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );

   for (i = 0; i < (*consdata)->n_linvar; ++i)
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->linvar[i]) );
   for (i = 0; i < (*consdata)->n_quadvar; ++i)
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->quadvar[i]) );

   consdataFree(scip, consdata);

   assert(*consdata == NULL);

   return SCIP_OKAY;
}
#else
#define consDeleteQuadratic NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */
#if 1
static
SCIP_DECL_CONSTRANS(consTransQuadratic)
{  
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;
   int i;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   SCIP_CALL( SCIPduplicateMemory(scip, &targetdata, sourcedata) );

   if (targetdata->n_linvar)
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &targetdata->lincoeff, sourcedata->lincoeff, targetdata->n_linvar) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &targetdata->linvar, targetdata->n_linvar) );
      SCIP_CALL( SCIPgetTransformedVars(scip, targetdata->n_linvar, sourcedata->linvar, targetdata->linvar) );
      for (i = 0; i < targetdata->n_linvar; ++i)
         SCIP_CALL( SCIPcaptureVar(scip, targetdata->linvar[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &targetdata->linrange, targetdata->n_linvar) );
   }
   else
   {
      targetdata->linvar   = NULL;
      targetdata->lincoeff = NULL;
      targetdata->linrange = NULL;
   }
   targetdata->linbndchgeventdata  = NULL;
   targetdata->quadbndchgeventdata = NULL;

   if (targetdata->n_quadvar)
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &targetdata->quadlincoeff, sourcedata->quadlincoeff, targetdata->n_quadvar) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &targetdata->quadsqrcoeff, sourcedata->quadsqrcoeff, targetdata->n_quadvar) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &targetdata->n_adjbilin,   sourcedata->n_adjbilin,   targetdata->n_quadvar) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &targetdata->quadvar,  targetdata->n_quadvar) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &targetdata->adjbilin, targetdata->n_quadvar) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &targetdata->quadrangevar, targetdata->n_quadvar) );
      SCIP_CALL( SCIPgetTransformedVars(scip, targetdata->n_quadvar, sourcedata->quadvar, targetdata->quadvar) );
      for (i = 0; i < targetdata->n_quadvar; ++i)
      {
         SCIP_CALL( SCIPcaptureVar(scip, targetdata->quadvar[i]) );
         if (targetdata->n_adjbilin[i])
            SCIP_CALL( SCIPduplicateMemoryArray(scip, &targetdata->adjbilin[i], sourcedata->adjbilin[i], sourcedata->n_adjbilin[i]) );
         else
            targetdata->adjbilin[i] = NULL;
      }

      if (targetdata->n_bilin)
      {
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &targetdata->bilincoeff, sourcedata->bilincoeff, targetdata->n_bilin) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &targetdata->bilinvar1, targetdata->n_bilin) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &targetdata->bilinvar2, targetdata->n_bilin) );
         SCIP_CALL( SCIPgetTransformedVars(scip, targetdata->n_bilin, sourcedata->bilinvar1, targetdata->bilinvar1) );
         SCIP_CALL( SCIPgetTransformedVars(scip, targetdata->n_bilin, sourcedata->bilinvar2, targetdata->bilinvar2) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &targetdata->bilinrange, targetdata->n_bilin) );
      }
      else
      {
         targetdata->bilincoeff = NULL;
         targetdata->bilinvar1  = NULL;
         targetdata->bilinvar2  = NULL;
         targetdata->bilinrange = NULL;
      }
   }
   else
   {
      targetdata->quadvar      = NULL;
      targetdata->quadlincoeff = NULL;
      targetdata->quadsqrcoeff = NULL;
      targetdata->n_adjbilin   = NULL;
      targetdata->adjbilin     = NULL;
      targetdata->bilincoeff   = NULL;
      targetdata->bilinvar1    = NULL;
      targetdata->bilinvar2    = NULL;
      targetdata->bilinrange   = NULL;
      targetdata->quadrangevar = NULL;
   }

   targetdata->is_removedfixings = FALSE;
   targetdata->is_propagated     = FALSE;
   targetdata->is_presolved      = FALSE;
   targetdata->soc_added         = FALSE;

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
      SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
      SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),  SCIPconsIsLocal(sourcecons),
      SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
      SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}
#else
#define consTransQuadratic NULL
#endif


/** LP initialization method of constraint handler */
#if 1
static
SCIP_DECL_CONSINITLP(consInitlpQuadratic)
{  
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_ROW*          row;
   int                i;

   assert(scip  != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for (i = 0; i < nconss; ++i)
   {
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      if (!SCIPisInfinity(scip, -consdata->lhs))
      {
         SCIP_CALL( generateCut(scip, conss[i], NULL, SCIP_BOUNDTYPE_LOWER, &row, conshdlrdata->cutmaxrange) );
         if (row)
         {
            SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE /* forcecut */) );
#ifdef SCIP_DEBUG
            SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif
            SCIP_CALL( SCIPreleaseRow (scip, &row) );
         }
      }
      if (!SCIPisInfinity(scip, consdata->rhs))
      {
         SCIP_CALL( generateCut(scip, conss[i], NULL, SCIP_BOUNDTYPE_UPPER, &row, conshdlrdata->cutmaxrange) );
         if (row)
         {
            SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE /* forcecut */) );
#ifdef SCIP_DEBUG
            SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif
            SCIP_CALL( SCIPreleaseRow (scip, &row) );
         }
      }
   }

   return SCIP_OKAY;
}
#else
#define consInitlpQuadratic NULL
#endif


/** separation method of constraint handler for LP solutions */
#if 1
static
SCIP_DECL_CONSSEPALP(consSepalpQuadratic)
{  
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;
   SCIP_Bool          intervalfeas;

   assert(scip   != NULL);
   assert(conshdlr != NULL);
   assert(conss  != NULL || nconss == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conss, nconss, NULL, conshdlrdata->do_scaling, &maxviolcon) );
   if (!maxviolcon)
      return SCIP_OKAY;

   SCIP_CALL( areIntervalFeasible(scip, conss, nconss, maxviolcon, &intervalfeas) );
   if (!intervalfeas)
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, result, TRUE) );
   if (*result == SCIP_SEPARATED)
      return SCIP_OKAY;

   return SCIP_OKAY;
}
#else
#define consSepalpQuadratic NULL
#endif


/** separation method of constraint handler for arbitrary primal solutions */
#if 1
static
SCIP_DECL_CONSSEPASOL(consSepasolQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;
   SCIP_Bool          intervalfeas;
   
   assert(scip   != NULL);
   assert(conshdlr != NULL);
   assert(conss  != NULL || nconss == 0);
   assert(sol    != NULL);
   assert(result != NULL);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( computeViolations(scip, conss, nconss, sol, conshdlrdata->do_scaling, &maxviolcon) );
   if (!maxviolcon)
      return SCIP_OKAY;

   SCIP_CALL( areIntervalFeasible(scip, conss, nconss, maxviolcon, &intervalfeas) );
   if (!intervalfeas)
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, sol, result, TRUE) );

   return SCIP_OKAY;
}
#else
#define consSepasolQuadratic NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;
   SCIP_Bool          intervalfeas;
   SCIP_RESULT        separateresult;
   int                nnotify;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(nconss == 0 || conss != NULL);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   
   SCIP_CALL( computeViolations(scip, conss, nconss, NULL, conshdlrdata->do_scaling, &maxviolcon) );
   if (!maxviolcon)
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }
   
   *result = SCIP_INFEASIBLE;

   SCIP_CALL( areIntervalFeasible(scip, conss, nconss, maxviolcon, &intervalfeas) );
   if (!intervalfeas)
   {
       *result = SCIP_CUTOFF;
       return SCIP_OKAY;
   }
   
   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, &separateresult, FALSE) );
   if (separateresult == SCIP_SEPARATED)
   {
      *result = SCIP_SEPARATED;
      return SCIP_OKAY;
   }

   /* we are not feasible, the whole node is not infeasible, and we cannot find a cut
    * -> collect variables for branching
    */
   
   SCIPdebugMessage("separation failed; max viol: %g+%g\n", SCIPconsGetData(maxviolcon)->lhsviol, SCIPconsGetData(maxviolcon)->rhsviol);

   SCIP_CALL( registerVariableInfeasibilities(scip, conshdlr, conss, nconss, &nnotify) );
   if (!nnotify)
   { /* fallback: separation probably failed because of numerical difficulties with a convex constraint; try to resolve by branching */ 
      SCIP_VAR* brvar = NULL;
      SCIP_CALL( registerLargeLPValueVariableForBranching(scip, conshdlr, conss, nconss, &brvar) );
      if (!brvar)
      {
         SCIPwarningMessage("Could not find any branching variable candidate. Cutting off node.\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      else
      {
         SCIPdebugMessage("Could not find any usual branching variable candidate. Proposed variable %s with LP value %g for branching.\n", SCIPvarGetName(brvar), SCIPgetSolVal(scip, NULL, brvar));
      }
   }
#ifndef WITH_CONSBRANCHNL
   SCIP_CALL( enforceByBranching(scip, conshdlr, result) );
   assert(*result == SCIP_BRANCHED);
#endif
   
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;
   SCIP_CONSDATA*     consdata;
   SCIP_Bool          intervalfeas;
   int                c, i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(nconss == 0 || conss != NULL);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   
   SCIP_CALL( computeViolations(scip, conss, nconss, NULL, conshdlrdata->do_scaling, &maxviolcon) );
   if (!maxviolcon)
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }
   
   *result = SCIP_INFEASIBLE;
      
   SCIP_CALL( areIntervalFeasible(scip, conss, nconss, maxviolcon, &intervalfeas) );
   if (!intervalfeas)
   {
       *result = SCIP_CUTOFF;
       return SCIP_OKAY;
   }
   
   /* we are not feasible and the whole node is not infeasible
    * -> collect variables for branching
    */

   for (c = 0; c < nconss; ++c)
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      SCIPdebugMessage("con %s violation: %g %g  convex: %d %d\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol, consdata->is_convex, consdata->is_concave);
      
      if (!consdata->n_quadvar)
         continue;
      
      if (!SCIPisFeasPositive(scip, consdata->lhsviol) && !SCIPisFeasPositive(scip, consdata->rhsviol))
         continue;
      
      SCIPdebugMessage("con %s violation: %g %g\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol);
      
      for (i = 0; i < consdata->n_quadvar; ++i)
         if (!SCIPisEQ(scip, SCIPvarGetLbLocal(consdata->quadvar[i]), SCIPvarGetUbLocal(consdata->quadvar[i])))
#ifdef WITH_CONSBRANCHNL
            SCIP_CALL( SCIPconshdlrBranchNonlinearUpdateVarInfeasibility(scip, conshdlrdata->branchnl, consdata->quadvar[i], consdata->lhsviol + consdata->rhsviol) );
#else
            SCIP_CALL( updateVarInfeasibility(scip, conshdlr, consdata->quadvar[i], consdata->lhsviol + consdata->rhsviol) );
#endif      
   }
#ifndef WITH_CONSBRANCHNL
   SCIP_CALL( enforceByBranching(scip, conshdlr, result) );
   assert(*result == SCIP_BRANCHED);
#endif
   
   return SCIP_OKAY;
}



/** domain propagation method of constraint handler */
#if 1
static
SCIP_DECL_CONSPROP(consPropQuadratic)
{
   SCIP_RESULT propresult;
   int         c;
   int         nchgbds = 0;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL);

   *result = SCIP_DIDNOTFIND;

   for (c = 0; c < nconss && *result != SCIP_CUTOFF; ++c)
   {
      SCIP_CALL( propagateBounds(scip, conshdlr, conss[c], &propresult, &nchgbds) );
      if (propresult != SCIP_DIDNOTFIND && propresult != SCIP_DIDNOTRUN)
         *result = propresult;
   }

   return SCIP_OKAY;
}
#else
#define consPropQuadratic NULL
#endif


/** presolving method of constraint handler */
#if 1
static
SCIP_DECL_CONSPRESOL(consPresolQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   int                c;
   SCIP_Bool          have_change;
   SCIP_HASHMAP*      terms;
   SCIP_Real          constant;
   
   assert(scip   != NULL);
   assert(conshdlr != NULL);
   assert(conss  != NULL || nconss == 0);
   assert(result != NULL);
   
   *result = SCIP_DIDNOTFIND;
   
   if (nrounds && !nnewfixedvars && !nnewupgdconss && !nnewchgbds && !nnewaggrvars && !nnewchgvartypes)
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for (c = 0; c < nconss; ++c)
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      have_change = FALSE;
      if (!consdata->is_removedfixings || !consdata->is_presolved || (!consdata->is_propagated && conshdlrdata->replace_binaryprod_length))
      {
         SCIPdebugMessage("process constraint %s\n", SCIPconsGetName(conss[c]));
         terms = NULL;
         SCIP_CALL( presolveCreateQuadTerm(scip, &terms, &constant, &have_change, conss[c], conshdlrdata->replace_sqrbinary) );
#ifdef SCIP_DEBUG
         SCIPdebugMessage("%d; %g ", have_change, constant);
         presolveQuadTermPrint(scip, terms, NULL);
         SCIPdebugPrintf("\n");
#endif

         if (conshdlrdata->replace_binaryprod_length)
         {
            int n_consadded = 0;
#if 0
            if (conshdlrdata->replace_binaryprod_forceAND)
            {
               SCIP_CALL( presolveTryAddAND(scip, conss[c], terms, &n_consadded) );
               if (n_consadded)
               { /* does this count as an upgrade? */
                  *result = SCIP_SUCCESS;
                  have_change = TRUE;
               }
            }
#endif
            SCIP_CALL( presolveTryAddLinearReform(scip, conss[c], terms, &n_consadded, conshdlrdata->replace_binaryprod_length) );
            if (n_consadded)
            { /* does this count as an upgrade? */
               *result = SCIP_SUCCESS;
               have_change = TRUE;
            }
         }

         if (have_change || !consdata->is_presolved)
         {
            int n_components; /* number of quadratic blocks */
            int i;

            if (SCIPhashmapIsEmpty(terms))
            { /* all variables fixed or removed, constraint now constant */
               presolveQuadTermFree(scip, &terms);
               if ( (!SCIPisInfinity(scip, -consdata->lhs) && SCIPisFeasLT(scip, consdata->lhs, constant    )) ||
                    (!SCIPisInfinity(scip,  consdata->rhs) && SCIPisFeasGT(scip, constant,     consdata->rhs)) )
               { /* constant is out of bounds */
                  SCIPdebugMessage("constraint %s is constant and infeasible\n", SCIPconsGetName(conss[c]));
                  *result = SCIP_CUTOFF;
                  break;
               }
               else
               { /* constant inside bounds, i.e., constraint always feasible */
                  SCIPdebugMessage("constraint %s is constant and feasible, removing\n", SCIPconsGetName(conss[c]));
                  SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) ); /* because only the presolved form is constant, in consdata we might still have variables */
                  SCIP_CALL( SCIPdelCons(scip, conss[c]) );
                  ++*ndelconss;
                  *result = SCIP_SUCCESS;
                  continue;
               }
            }
#ifdef WITH_SOC3
            if (conshdlrdata->upgrade_soc)
            {
               SCIP_CALL( presolveTryUpgradeSOC(scip, conss[c], terms, constant) );
               if (SCIPisInfinity(scip, -consdata->lhs) && SCIPisInfinity(scip, consdata->rhs))
               { /* can delete constraint here since it was upgraded to SOC3 */
                  SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );
                  SCIP_CALL( SCIPdelCons(scip, conss[c]) );
                  presolveQuadTermFree(scip, &terms);
                  continue;
               }
            }
            
            if (conshdlrdata->soc3_nr_auxvars >= 0 && !consdata->soc_added)
            { /* try to add outer-approximation by SOC constraint */
               SCIP_CALL( presolveTryAddSOC(scip, conss[c], terms, constant, conshdlrdata->soc3_nr_auxvars) );
               if (consdata->soc_added)
               {
                  ++*nupgdconss;
                  *result = SCIP_SUCCESS;
                  
                  if (SCIPisInfinity(scip, -consdata->lhs) && SCIPisInfinity(scip, consdata->rhs))
                  { /* can delete constraint here since it was upgraded to SOC3 */
                     SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );
                     SCIP_CALL( SCIPdelCons(scip, conss[c]) );
                     presolveQuadTermFree(scip, &terms);
                     continue;
                  }
               }
            }
#endif
            
            SCIP_CALL( presolveFindComponents(scip, terms, &n_components) );
            assert(n_components >= 0);

            if (n_components == 0)
            { /* upgrade to linear constraint */
               SCIP_HASHMAPLIST* list;
               SCIP_VAR*         var;
               PresolveQuadTerm* term;
               SCIP_CONS*        lincon;
               
               SCIPdebugMessage("upgrade to linear constraint\n");
               SCIP_CALL( SCIPcreateConsLinear(scip, &lincon, SCIPconsGetName(conss[c]), 0, NULL, NULL,
                  SCIPisInfinity(scip, -consdata->lhs) ? consdata->lhs : consdata->lhs-constant,
                  SCIPisInfinity(scip,  consdata->rhs) ? consdata->rhs : consdata->rhs-constant,
                  SCIPconsIsInitial(conss[c]), SCIPconsIsSeparated(conss[c]), SCIPconsIsEnforced(conss[c]),
                  SCIPconsIsChecked(conss[c]), SCIPconsIsPropagated(conss[c]),  SCIPconsIsLocal(conss[c]),
                  SCIPconsIsModifiable(conss[c]), SCIPconsIsDynamic(conss[c]), SCIPconsIsRemovable(conss[c]),
                  SCIPconsIsStickingAtNode(conss[c])) );
               for (i = 0; i < SCIPhashmapGetNLists(terms); ++i)
                  for (list = SCIPhashmapGetList(terms, i); list; list = SCIPhashmapListGetNext(list))
                  {
                     var  = (SCIP_VAR*)        SCIPhashmapListGetOrigin(list);
                     term = (PresolveQuadTerm*)SCIPhashmapListGetImage(list);
                     assert(SCIPisZero(scip, term->sqrcoeff));
                     assert(term->bilin == NULL);
                     if (SCIPisZero(scip, term->lincoeff))
                        continue;
                     SCIP_CALL( SCIPaddCoefLinear(scip, lincon, var, term->lincoeff) );
                  }

#ifdef SCIP_DEBUG
               SCIP_CALL( SCIPprintCons(scip, lincon, NULL) );
#endif
               SCIP_CALL( SCIPaddCons(scip, lincon) );
               SCIP_CALL( SCIPreleaseCons(scip, &lincon) );
               SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );
               SCIP_CALL( SCIPdelCons(scip, conss[c]) );
               ++*nupgdconss;
               *result = SCIP_SUCCESS;
               presolveQuadTermFree(scip, &terms);
               continue;
            }

            if (have_change || (n_components > 1 && conshdlrdata->disaggregation))
            {
               /* unlock all variables */
               for (i = 0; i < consdata->n_linvar; ++i)
                  SCIP_CALL( unlockLinearVariable(scip, conss[c], consdata->linvar[i], consdata->lincoeff[i]) );
               for (i = 0; i < consdata->n_quadvar; ++i)
                  SCIP_CALL( SCIPunlockVarCons(scip, consdata->quadvar[i], conss[c], TRUE, TRUE) );
               SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );

               if (n_components == 1 || !conshdlrdata->disaggregation)
               {
                  SCIP_CALL( consdataSetFunctionData(scip, consdata, terms, -1) );
                  if (!SCIPisInfinity(scip, -consdata->lhs))
                     consdata->lhs -= constant;
                  if (!SCIPisInfinity(scip,  consdata->rhs))
                     consdata->rhs -= constant;
               }
               else
               {
                  SCIP_CALL( presolveDisaggregate(scip, conshdlr, conss[c], terms, constant, n_components) );
                  have_change = TRUE;
               }

               /* lock all variables */
               for (i = 0; i < consdata->n_linvar; ++i)
                  SCIP_CALL( lockLinearVariable(scip, conss[c], consdata->linvar[i], consdata->lincoeff[i]) );
               for (i = 0; i < consdata->n_quadvar; ++i)
                  SCIP_CALL( SCIPlockVarCons(scip, consdata->quadvar[i], conss[c], TRUE, TRUE) );
               SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );
            }
         }

         presolveQuadTermFree(scip, &terms);

         consdata->is_removedfixings = TRUE;
      }

      if (!consdata->is_propagated)
      { /* try domain propagation if there were bound changes or constraint has changed (in which case, catchVarEvents has set is_propagated to false) */
         SCIP_RESULT propresult;
         SCIP_CALL( propagateBounds(scip, conshdlr, conss[c], &propresult, nchgbds) );
         switch (propresult)
         {
            case SCIP_REDUCEDDOM:
               *result = SCIP_SUCCESS;
               break;
            case SCIP_CUTOFF:
               SCIPdebugMessage("propagation on constraint %s says problem is infeasible in presolve\n", SCIPconsGetName(conss[c]));
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
               break;
            default:
               assert(propresult == SCIP_DIDNOTFIND || propresult == SCIP_DIDNOTRUN);
         }

         if (propresult != SCIP_REDUCEDDOM && !SCIPisInfinity(scip, conshdlrdata->defaultbound))
         {
            if (nrounds > 0)
            {
               int nboundchanges = 0;
               SCIP_CALL( boundUnboundedVars(scip, conss[c], conshdlrdata->defaultbound, &nboundchanges) );
               if (nboundchanges)
               {
                  *nchgbds += nboundchanges;
                  *result   = SCIP_SUCCESS;
               }
            }
            else
            { /* wait for next round (or do in exitpre if no next round) */
               consdata->is_propagated = FALSE;
            }
         }
      }

      if ((nnewchgvartypes || have_change || !consdata->is_presolved)
         && (SCIPisEQ(scip, consdata->lhs, consdata->rhs) && SCIPisIntegral(scip, consdata->lhs)))
      { /* check if we have a single linear continuous variable that we can make implicit integer 
        * TODO allow for coefficient != +/-1 in front of linear var
        */
         int       ncontvar  = 0;
         SCIP_VAR* candidate = NULL;
         SCIP_Bool fail      = FALSE;
         int       i;
         
         for (i = 0; !fail && i < consdata->n_linvar; ++i)
            if (!SCIPisIntegral(scip, consdata->lincoeff[i]))
               fail = TRUE;
            else if (SCIPvarGetType(consdata->linvar[i]) == SCIP_VARTYPE_CONTINUOUS)
            {
               if (ncontvar) /* now at 2nd continuous variable */
                  fail = TRUE;
               else if (SCIPisEQ(scip, ABS(consdata->lincoeff[i]), 1.0))
                  candidate = consdata->linvar[i];
               ++ncontvar;
            }
         for (i = 0; !fail && i < consdata->n_quadvar; ++i)
            fail = SCIPvarGetType(consdata->quadvar[i]) == SCIP_VARTYPE_CONTINUOUS ||
                  !SCIPisIntegral(scip, consdata->quadlincoeff[i]) ||
                  !SCIPisIntegral(scip, consdata->quadsqrcoeff[i]);
         for (i = 0; !fail && i < consdata->n_bilin; ++i)
            fail = !SCIPisIntegral(scip, consdata->bilincoeff[i]);

         if (!fail && candidate)
         {
            SCIPdebugMessage("make variable %s implicit integer due to constraint %s\n", SCIPvarGetName(candidate), SCIPconsGetName(conss[c]));
            SCIP_CALL( SCIPchgVarType(scip, candidate, SCIP_VARTYPE_IMPLINT) );
            ++*nchgvartypes;
            *result = SCIP_SUCCESS;
         }
      }

      consdata->is_presolved = TRUE;
   }

   return SCIP_OKAY;
}
#else
#define consPresolQuadratic NULL
#endif


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropQuadratic)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of quadratic constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropQuadratic NULL
#endif

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockQuadratic)
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool      haslb, hasub;
   int            i;
   
   assert(scip != NULL);
   assert(cons != NULL);
  
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   haslb = !SCIPisInfinity(scip, -consdata->lhs);
   hasub = !SCIPisInfinity(scip,  consdata->rhs);

   for (i = 0; i < consdata->n_linvar; ++i)
   {
      if (consdata->lincoeff[i] > 0)
      {
         if( haslb )
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvar[i], nlockspos, nlocksneg) );
         if( hasub )
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvar[i], nlocksneg, nlockspos) );
      }
      else
      {
         if( haslb )
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvar[i], nlocksneg, nlockspos) );
         if( hasub )
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvar[i], nlockspos, nlocksneg) );
      }
   }

   for (i = 0; i < consdata->n_quadvar; ++i)
   { /* TODO try to be more clever */
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->quadvar[i], nlockspos+nlocksneg, nlockspos+nlocksneg) );
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveQuadratic)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of quadratic constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveQuadratic NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveQuadratic)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of quadratic constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveQuadratic NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableQuadratic)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of quadratic constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableQuadratic NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableQuadratic)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of quadratic constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableQuadratic NULL
#endif


/** constraint display method of constraint handler */
#if 1
static
SCIP_DECL_CONSPRINT(consPrintQuadratic)
{  
   SCIP_CONSDATA* consdata;
   int            j;
   
   assert(scip != NULL);
   assert(cons != NULL);
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* print left hand side for ranged rows */
   if( !SCIPisInfinity(scip, -consdata->lhs)
      && !SCIPisInfinity(scip, consdata->rhs)
      && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, "%.15g <= ", consdata->lhs);

   /* print coefficients and variables */
   if( consdata->n_linvar == 0 && consdata->n_quadvar == 0)
      SCIPinfoMessage(scip, file, "0 ");
   else
   {
      for (j = 0; j < consdata->n_linvar; ++j)
         SCIPinfoMessage(scip, file, "%+.15g<%s>[%c] ", consdata->lincoeff[j], SCIPvarGetName(consdata->linvar[j]),
            SCIPvarGetType(consdata->linvar[j]) == SCIP_VARTYPE_BINARY ? 'B' :
            SCIPvarGetType(consdata->linvar[j]) == SCIP_VARTYPE_INTEGER ? 'I' :
            SCIPvarGetType(consdata->linvar[j]) == SCIP_VARTYPE_IMPLINT ? 'I' : 'C');

      for (j = 0; j < consdata->n_quadvar; ++j)
      {
         if (consdata->quadlincoeff[j])
            SCIPinfoMessage(scip, file, "%+.15g<%s>[%c]", consdata->quadlincoeff[j], SCIPvarGetName(consdata->quadvar[j]),
               SCIPvarGetType(consdata->quadvar[j]) == SCIP_VARTYPE_BINARY ? 'B' :
               SCIPvarGetType(consdata->quadvar[j]) == SCIP_VARTYPE_INTEGER ? 'I' :
               SCIPvarGetType(consdata->quadvar[j]) == SCIP_VARTYPE_IMPLINT ? 'I' : 'C');
         if (consdata->quadsqrcoeff[j])
            SCIPinfoMessage(scip, file, "%+.15gsqr(<%s>[%c])", consdata->quadsqrcoeff[j], SCIPvarGetName(consdata->quadvar[j]),
               SCIPvarGetType(consdata->quadvar[j]) == SCIP_VARTYPE_BINARY ? 'B' :
               SCIPvarGetType(consdata->quadvar[j]) == SCIP_VARTYPE_INTEGER ? 'I' :
               SCIPvarGetType(consdata->quadvar[j]) == SCIP_VARTYPE_IMPLINT ? 'I' : 'C');
      }

      for (j = 0; j < consdata->n_bilin; ++j)
      {
         SCIPinfoMessage(scip, file, "%+.15g<%s>[%c]*<%s>[%c]", consdata->bilincoeff[j],
            SCIPvarGetName(consdata->bilinvar1[j]),
            SCIPvarGetType(consdata->bilinvar1[j]) == SCIP_VARTYPE_BINARY ? 'B' :
            SCIPvarGetType(consdata->bilinvar1[j]) == SCIP_VARTYPE_INTEGER ? 'I' :
            SCIPvarGetType(consdata->bilinvar1[j]) == SCIP_VARTYPE_IMPLINT ? 'I' : 'C',
            SCIPvarGetName(consdata->bilinvar2[j]),
            SCIPvarGetType(consdata->bilinvar2[j]) == SCIP_VARTYPE_BINARY ? 'B' :
            SCIPvarGetType(consdata->bilinvar2[j]) == SCIP_VARTYPE_INTEGER ? 'I' :
            SCIPvarGetType(consdata->bilinvar2[j]) == SCIP_VARTYPE_IMPLINT ? 'I' : 'C');
      }
   }

   /* print right hand side */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, " == %.15g\n", consdata->rhs);
   else if( !SCIPisInfinity(scip, consdata->rhs) )
      SCIPinfoMessage(scip, file, " <= %.15g\n", consdata->rhs);
   else if( !SCIPisInfinity(scip, -consdata->lhs) )
      SCIPinfoMessage(scip, file, " >= %.15g\n", consdata->lhs);
   else
      SCIPinfoMessage(scip, file, " [free]\n");

   return SCIP_OKAY;
}
#else
#define consPrintQuadratic NULL
#endif

/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          maxviol = 0.0;
   int                c;

   assert(scip != NULL);
   assert(nconss == 0 || conss != NULL);
   assert(sol != NULL);
   assert(result != NULL);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_FEASIBLE;

   for (c = 0; c < nconss; ++c)
   {
      SCIP_CALL( computeViolation(scip, conss[c], sol, conshdlrdata->do_scaling) );
      
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      
      if (SCIPisFeasPositive(scip, consdata->lhsviol) || SCIPisFeasPositive(scip, consdata->rhsviol))
      {
         *result = SCIP_INFEASIBLE;
         if (printreason)
         {
            SCIPinfoMessage(scip, NULL, "quadratic constraint %s violated by %g+%g\n\t", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol);
            SCIP_CALL( consPrintQuadratic(scip, conshdlr, conss[c], NULL) );
         }
         if (!conshdlrdata->nlpheur)
            return SCIP_OKAY;
         if (consdata->lhsviol > maxviol || consdata->rhsviol > maxviol)
            maxviol = consdata->lhsviol + consdata->rhsviol;
      }
   }
   
   if (*result == SCIP_INFEASIBLE && conshdlrdata->nlpheur)
      SCIP_CALL( SCIPheurNlpUpdateStartpoint(scip, conshdlrdata->nlpheur, sol, maxviol) );

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
#if 1
static
SCIP_DECL_CONSCOPY(consCopyQuadratic)
{  
   SCIP_CONSDATA* consdata;
   SCIP_VAR**     linvars   = NULL;
   SCIP_VAR**     quadvars  = NULL;
   SCIP_VAR**     bilinvar1 = NULL;
   SCIP_VAR**     bilinvar2 = NULL;
   int            i, j;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(sourcescip != NULL);
   assert(sourcecons != NULL);
   assert(varmap != NULL);
   assert(succeed != NULL);
   
   consdata = SCIPconsGetData(sourcecons);
   assert(consdata != NULL);
   
   *succeed = TRUE; /* think positive */
   
   if (consdata->n_linvar)
   {
      SCIP_CALL( SCIPallocBufferArray(sourcescip, &linvars, consdata->n_linvar) );
      for (i = 0; i < consdata->n_linvar; ++i)
      {
         linvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, consdata->linvar[i]);
         if (linvars[i] == NULL)
         {
            *succeed = FALSE;
            break;
         }
      }
   }
   
   if (consdata->n_bilin && *succeed)
   {
      SCIP_CALL( SCIPallocBufferArray(sourcescip, &bilinvar1, consdata->n_bilin) );
      SCIP_CALL( SCIPallocBufferArray(sourcescip, &bilinvar2, consdata->n_bilin) );
   }

   if (consdata->n_quadvar && *succeed)
   {
      SCIP_CALL( SCIPallocBufferArray(sourcescip, &quadvars, consdata->n_quadvar) );
      for (i = 0; i < consdata->n_quadvar; ++i)
      {
         quadvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, consdata->quadvar[i]);
         if (quadvars[i] == NULL)
         {
            *succeed = FALSE;
            break;
         }
         
         assert(consdata->n_bilin || consdata->n_adjbilin[i] == 0);
         
         for (j = 0; j < consdata->n_adjbilin[i]; ++j)
         {
            if (consdata->bilinvar1[consdata->adjbilin[i][j]] == consdata->quadvar[i])
            {
               assert(consdata->bilinvar2[consdata->adjbilin[i][j]] != consdata->quadvar[i]);
               bilinvar1[consdata->adjbilin[i][j]] = quadvars[i];
            }
            else
            {
               assert(consdata->bilinvar2[consdata->adjbilin[i][j]] == consdata->quadvar[i]);
               bilinvar2[consdata->adjbilin[i][j]] = quadvars[i];
            }
         }
      }
   }

   if (*succeed)
   {
      assert(stickingatnode == FALSE);
      SCIP_CALL( SCIPcreateConsQuadratic2(scip, cons, name ? name : SCIPconsGetName(sourcecons),
         consdata->n_linvar, linvars, consdata->lincoeff,
         consdata->n_quadvar, quadvars, consdata->quadlincoeff, consdata->quadsqrcoeff,
         consdata->n_adjbilin, consdata->adjbilin,
         consdata->n_bilin, bilinvar1, bilinvar2, consdata->bilincoeff,
         consdata->lhs, consdata->rhs,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );
   }
   else
      *cons = NULL;
   
   SCIPfreeBufferArrayNull(sourcescip, &linvars);
   SCIPfreeBufferArrayNull(sourcescip, &quadvars);
   SCIPfreeBufferArrayNull(sourcescip, &bilinvar1);
   SCIPfreeBufferArrayNull(sourcescip, &bilinvar2);
   
   return SCIP_OKAY;
}
#else
#define consCopyQuadratic NULL
#endif


/** constraint parsing method of constraint handler */
#if 0
static
SCIP_DECL_CONSPARSE(consParseQuadratic)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of quadratic constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consParseQuadratic NULL
#endif


/*
 * constraint specific interface methods
 */

/** creates the handler for quadratic constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrQuadratic(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create quadratic constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   memset(conshdlrdata, 0, sizeof(SCIP_CONSHDLRDATA));

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         consFreeQuadratic, consInitQuadratic, consExitQuadratic,
         consInitpreQuadratic, consExitpreQuadratic, consInitsolQuadratic, consExitsolQuadratic,
         consDeleteQuadratic, consTransQuadratic, consInitlpQuadratic,
         consSepalpQuadratic, consSepasolQuadratic, consEnfolpQuadratic, consEnfopsQuadratic, consCheckQuadratic,
         consPropQuadratic, consPresolQuadratic, consRespropQuadratic, consLockQuadratic,
         consActiveQuadratic, consDeactiveQuadratic,
         consEnableQuadratic, consDisableQuadratic,
         consPrintQuadratic, consCopyQuadratic, consParseQuadratic,
         conshdlrdata) );

   /* add quadratic constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/replace_sqrbinary",           "whether a square of a binary variables should be replaced by the binary variable",                                                                           &conshdlrdata->replace_sqrbinary,           FALSE, TRUE,                                        NULL, NULL) );
/*   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/replace_binaryprod",          "whether products with a binary variable should be replaced by a new variable and an AND constraint (if possible) or a set of equivalent linear constraints", &conshdlrdata->replace_binaryprod,          FALSE, TRUE,                                        NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/replace_binaryprod_forceAND", "whether products of binary variables should be replaced always by an AND constraint",                                                                        &conshdlrdata->replace_binaryprod_forceAND, FALSE, FALSE,                                       NULL, NULL) );
*/
   SCIP_CALL( SCIPaddIntParam (scip, "constraints/"CONSHDLR_NAME"/replace_binaryprod",          "max. length of linear term which when multiplied with a binary variables is replaced by an auxiliary variable and a linear reformulation (0 to turn off)",   &conshdlrdata->replace_binaryprod_length,   FALSE, INT_MAX, 0, INT_MAX,                         NULL, NULL) );
#ifdef WITH_SOC3
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/upgrade_soc",                 "whether constraints should be upgraded to SOC constraints, if possible",                                                                                     &conshdlrdata->upgrade_soc,                 FALSE, FALSE,                                       NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam (scip, "constraints/"CONSHDLR_NAME"/soc3nrauxvar",                "number of auxiliary variables when adding a SOC3 constraint; 0 to determine automatically (could become very large!); -1 to turn off",                       &conshdlrdata->soc3_nr_auxvars,             FALSE, 8, -1, INT_MAX,                              NULL, NULL) );
#endif
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/disaggregate",                "whether quadratic constraints consisting of several quadratic blocks should be disaggregated in several constraints",                                        &conshdlrdata->disaggregation,              FALSE, TRUE,                                        NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/minefficacy",                 "minimal efficacy for a cut to be added to the LP; overwrites separating/efficacy",                                                                           &conshdlrdata->mincutefficacy,              FALSE, 0.0001, 0.0, SCIPinfinity(scip),             NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/scaling",                     "whether a quadratic constraint should be scaled w.r.t. the current gradient norm when checking for feasibility",                                             &conshdlrdata->do_scaling,                  FALSE, TRUE,                                        NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/fastpropagate",               "whether a propagation should be used that is faster in case of bilinear term, but also less efficient",                                                      &conshdlrdata->fast_propagate,              FALSE, TRUE,                                        NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/defaultbound",                "a default bound to impose on unbounded variables in quadratic terms (-defaultbound is used for missing lower bounds)",                                       &conshdlrdata->defaultbound,                TRUE,  SCIPinfinity(scip), 0.0, SCIPinfinity(scip), NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/cutmaxrange",                 "maximal range of a cut (maximal coefficient divided by minimal coefficient) in order to be added to LP relaxation",                                          &conshdlrdata->cutmaxrange,                 FALSE, 1e+10, 0.0, SCIPinfinity(scip),              NULL, NULL) );
#ifndef WITH_CONSBRANCHNL
   SCIP_CALL( SCIPaddCharParam(scip, "constraints/"CONSHDLR_NAME"/strategy",                    "strategy to use for selecting branching variable: b: rb-int-br, r: rb-int-br-rev, i: rb-inf",                                                                &conshdlrdata->strategy,                    FALSE, 'r', "bri",                                  NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/mindistbrpointtobound",       "minimal fractional distance of branching point to variable bounds; a value of 0.5 leads to branching always in the middle of a bounded domain",              &conshdlrdata->mindistbrpointtobound,       FALSE, 0.2, 0.0001, 0.5,                            NULL, NULL) );
#endif
   
   SCIP_CALL( SCIPincludeEventhdlr(scip, CONSHDLR_NAME, "signals a bound change to a quadratic constraint",
      NULL, NULL, NULL, NULL, NULL, NULL, processVarEvent, NULL) );
   conshdlrdata->eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME);

   return SCIP_OKAY;
}

/** creates and captures a quadratic constraint */
SCIP_RETCODE SCIPcreateConsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   n_linvars,          /**< number of linear terms */
   SCIP_VAR**            linvars,            /**< variables in linear part */
   SCIP_Real*            lincoeff,           /**< coefficients of variables in linear part */
   int                   n_quadterm,         /**< number of quadratic terms */
   SCIP_VAR**            quadvars1,          /**< index of first variable in quadratic terms */
   SCIP_VAR**            quadvars2,          /**< index of second variable in quadratic terms */
   SCIP_Real*            quadcoeff,          /**< coefficients of quadratic terms */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation */
   SCIP_Real             rhs,                /**< right hand side of quadratic equation */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are seperated as constraints. */
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_HASHMAP*  terms;
   int i;

   /* find the quadratic constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("quadratic constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocMemory( scip, &consdata) );
   memset(consdata, 0, sizeof(SCIP_CONSDATA));

   consdata->lhs         = lhs;
   consdata->rhs         = rhs;

   SCIP_CALL( SCIPhashmapCreate(&terms, SCIPblkmem(scip), n_linvars + n_quadterm) );
   for( i = 0; i < n_linvars; ++i )
   {
      if( SCIPisZero(scip, lincoeff[i]) )
         continue;
      SCIP_CALL( presolveQuadTermAdd(scip, terms, linvars[i], lincoeff[i], 0., NULL, 0.) );
   }
   for( i = 0; i < n_quadterm; ++i )
   {
      if( SCIPisZero(scip, quadcoeff[i]) )
         continue;
      assert(SCIPgetStage(scip) < SCIP_STAGE_PRESOLVED || SCIPvarGetStatus(quadvars1[i]) != SCIP_VARSTATUS_MULTAGGR);
      assert(SCIPgetStage(scip) < SCIP_STAGE_PRESOLVED || SCIPvarGetStatus(quadvars2[i]) != SCIP_VARSTATUS_MULTAGGR);
      if( quadvars1[i] == quadvars2[i] )
      {
         SCIP_CALL( presolveQuadTermAdd(scip, terms, quadvars1[i], 0., quadcoeff[i], NULL, 0.) );
      }
      else
      {
         SCIP_CALL( presolveQuadTermAdd(scip, terms, quadvars1[i], 0., 0., quadvars2[i], quadcoeff[i]) );
         SCIP_CALL( presolveQuadTermAdd(scip, terms, quadvars2[i], 0., 0., quadvars1[i], quadcoeff[i]) );
      }
   }
   SCIP_CALL( consdataSetFunctionData(scip, consdata, terms, -1) );
   presolveQuadTermFree(scip, &terms);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, FALSE) );

   if( SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVING )
   {
      SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);
      
      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, *cons) );

      if( SCIPgetStage(scip) > SCIP_STAGE_INITSOLVE )
      {
         SCIP_CALL( checkCurvature(scip, conshdlr, *cons) );
         SCIPdebugMessage("new quadratic constraint %s is %sconvex and %sconcave\n", name, consdata->is_convex ? "" : "not ", consdata->is_concave ? "" : "not ");
      }
   }

   return SCIP_OKAY;
}

/** creates and captures a quadratic constraint */
SCIP_RETCODE SCIPcreateConsQuadratic2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   n_linvar,           /**< number of linear terms */
   SCIP_VAR**            linvar,             /**< variables in linear part */ 
   SCIP_Real*            lincoeff,           /**< coefficients of variables in linear part */ 
   int                   n_quadvar,          /**< number of quadratic terms */
   SCIP_VAR**            quadvar,            /**< variables in quadratic terms */
   SCIP_Real*            quadlincoeff,       /**< linear coefficients of quadratic variables */
   SCIP_Real*            quadsqrcoeff,       /**< coefficients of square terms of quadratic variables */
   int*                  n_adjbilin,         /**< number of bilinear terms where the variable is involved */
   int**                 adjbilin,           /**< indices of bilinear terms in which variable is involved */
   int                   n_bilin,            /**< number of bilinear terms */
   SCIP_VAR**            bilinvar1,          /**< first variable in bilinear term */
   SCIP_VAR**            bilinvar2,          /**< second variable in bilinear term */
   SCIP_Real*            bilincoeff,         /**< coefficient of bilinear term */
   SCIP_Real             lhs,                /**< constraint  left hand side */
   SCIP_Real             rhs,                /**< constraint right hand side */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint dynamic? */
   SCIP_Bool             removable           /**< should the constraint be removed from the LP due to aging or cleanup? */
)
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   int i;
   
   assert( n_linvar  == 0 || (linvar    != NULL && lincoeff     != NULL) );
   assert( n_quadvar == 0 || (quadvar   != NULL && quadlincoeff != NULL && quadsqrcoeff != NULL) );
   assert( n_bilin   == 0 || (bilinvar1 != NULL && bilinvar2    != NULL && bilincoeff   != NULL && n_quadvar > 0) );
   
   /* find the quadratic constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("quadratic constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocMemory( scip, &consdata) );
   memset(consdata, 0, sizeof(SCIP_CONSDATA));

   consdata->lhs = lhs;
   consdata->rhs = rhs;

   if (n_linvar)
   {
      consdata->n_linvar = n_linvar;
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->linvar,   linvar,   n_linvar) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->lincoeff, lincoeff, n_linvar) );
      for( i = 0; i < n_linvar; ++i )
         SCIP_CALL( SCIPcaptureVar(scip, linvar[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->linrange, n_linvar) );
   }

   if (n_quadvar)
   {
      consdata->n_quadvar = n_quadvar;
      SCIP_CALL(    SCIPduplicateMemoryArray(scip, &consdata->quadvar,      quadvar,      n_quadvar) );
      if (quadlincoeff)
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->quadlincoeff, quadlincoeff, n_quadvar) );
      else
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->quadlincoeff, n_quadvar) );
         memset(consdata->quadlincoeff, n_quadvar, n_quadvar * sizeof(SCIP_Real));
      }
      if (quadsqrcoeff)
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->quadsqrcoeff, quadsqrcoeff, n_quadvar) );
      else
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->quadsqrcoeff, n_quadvar) );
         memset(consdata->quadsqrcoeff, n_quadvar, n_quadvar * sizeof(SCIP_Real));
      }
      if (n_adjbilin)
      {
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->n_adjbilin,   n_adjbilin,   n_quadvar) );
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->adjbilin,     adjbilin,     n_quadvar) );
      }
      else
      {
         assert(n_bilin == 0);
         SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->n_adjbilin, n_quadvar) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->adjbilin, n_quadvar) );
         memset(consdata->n_adjbilin, 0, n_quadvar * sizeof(int ));
         memset(consdata->adjbilin,   0, n_quadvar * sizeof(int*));
      }
      for (i = 0; i < n_quadvar; ++i)
      {
         SCIP_CALL( SCIPcaptureVar(scip, quadvar[i]) );
         if (consdata->n_adjbilin[i])
         {
            assert(adjbilin[i] != NULL);
            SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->adjbilin[i], adjbilin[i], n_adjbilin[i]) );
         }
         else
         {
            assert(consdata->adjbilin[i] == NULL);
         }
      }
      SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->quadrangevar, n_quadvar) );
   }

   if (n_bilin)
   {
      consdata->n_bilin = n_bilin;
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->bilincoeff, bilincoeff, n_bilin) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->bilinvar1,  bilinvar1,  n_bilin) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->bilinvar2,  bilinvar2,  n_bilin) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->bilinrange, n_bilin) );
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
      local, modifiable, dynamic, removable, FALSE) );

   if( SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVING )
   {
      SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);
      
      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, *cons) );

      if( SCIPgetStage(scip) > SCIP_STAGE_INITSOLVE )
      {
         SCIP_CALL( checkCurvature(scip, conshdlr, *cons) );
         SCIPdebugMessage("new quadratic constraint %s is %sconvex and %sconcave\n", name, consdata->is_convex ? "" : "not ", consdata->is_concave ? "" : "not ");
      }
   }

   return SCIP_OKAY;
}

/** Gets the number of variables in the linear term of a quadratic constraint.
 */
int SCIPgetNLinearVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->n_linvar;
}

/** Gets the variables in the linear part of a quadratic constraint.
 * Length is given by SCIPgetNLinearVarsQuadratic.
 */
SCIP_VAR** SCIPgetLinearVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->linvar;
}

/** Gets the number of variables in the quadratic term of a quadratic constraint.
 */
int SCIPgetNQuadVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->n_quadvar;
}

/** Gets the variables in the quadratic part of a quadratic constraint.
 * Length is given by SCIPgetNQuadVarsQuadratic.
 */
SCIP_VAR** SCIPgetQuadVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->quadvar;
}

/** Gets the coefficients in the linear part of a quadratic constraint.
 * Length is given by SCIPgetNLinearVarsQuadratic.
 */
SCIP_Real* SCIPgetCoeffLinearVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->lincoeff;
}

/** Gets the linear coefficients in the quadratic part of a quadratic constraint.
 * Length is given by SCIPgetNQuadVarsQuadratic.
 */
SCIP_Real* SCIPgetLinearCoeffQuadVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->quadlincoeff;
}

/** Gets the square coefficients in the quadratic part of a quadratic constraint.
 * Length is given by SCIPgetNQuadVarsQuadratic.
 */
SCIP_Real* SCIPgetSqrCoeffQuadVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->quadsqrcoeff;
}

/** Gets the number of bilinear terms in a quadratic constraint.
 */
int SCIPgetNBilinTermQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->n_bilin;
}

/** Gets the first variables in the bilinear terms in a quadratic constraint.
 * Length is given by SCIPgetNBilinTermQuadratic.
 */
SCIP_VAR** SCIPgetBilinVar1Quadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->bilinvar1;
}

/** Gets the second variables in the bilinear terms in a quadratic constraint.
 * Length is given by SCIPgetNBilinTermQuadratic.
 */
SCIP_VAR** SCIPgetBilinVar2Quadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->bilinvar2;
}

/** Gets the coefficients of the bilinear terms in a quadratic constraint.
 * Length is given by SCIPgetNBilinTermQuadratic.
 */
SCIP_Real* SCIPgetBilinCoeffQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->bilincoeff;
}

/** Gets the left hand side of a quadratic constraint.
 */
SCIP_Real SCIPgetLhsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->lhs;
}

/** Gets the right hand side of a quadratic constraint.
 */
SCIP_Real SCIPgetRhsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->rhs;
}

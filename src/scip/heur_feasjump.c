/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_feasjump.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief  Feasibility Jump heuristic
 * @author Gioni Mexi
 * @author Dominik Kamp
 *
 * Feasibility Jump repeatedly selects a variable and sets its value to minimize the weighted sum of constraint
 * violations. The weights are adjusted for constraints that remain violated in local minima.
 *
 * Based on:
 *  Bjørnar Luteberget and Giorgio Sartor.
 *  "Feasibility Jump: an LP-free Lagrangian MIP heuristic"
 *  Mathematical Programming Computation 15.2 (2023): 365-388.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <limits.h>
#include <string.h>
#include <math.h>

#include "blockmemshell/memory.h"
#include "scip/heur_feasjump.h"
#include "scip/pub_cons.h"
#include "scip/pub_heur.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_var.h"
#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_general.h"
#include "scip/scip_heur.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_randnumgen.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_timing.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"

#define HEUR_NAME             "feasjump"
#define HEUR_DESC             "feasibility jump"
#define HEUR_DISPCHAR         SCIP_HEURDISPCHAR_LNS
#define HEUR_PRIORITY         900
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFOREPRESOL | SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      FALSE

#define DEFAULT_RANDSEED       29        /**< default random seed */
#define DEFAULT_WEIGHTDECAY    1.0       /**< decay factor for the weight update */
#define DEFAULT_RANDPROBGOOD   0.001     /**< probability of random selection from good variables */
#define DEFAULT_RANDPROBVIOL   0.01      /**< probability of random selection from violated constraint */
#define DEFAULT_SAMPLESIZE     25        /**< number of candidates to sample in tournament selection */
#define DEFAULT_VERBOSITY      0         /**< verbosity level of the feasibility jump solver */
#define DEFAULT_MAXSOLS        1         /**< maximum number of solutions to find (-1: unlimited) */
#define DEFAULT_MINDECREASEREL 0.1       /**< minimum relative decrease to reset the effort budget */
#define DEFAULT_MAXEFFORTFAC   128       /**< maximum effort factor per nonzero without improvement (-1: unlimited) */
#define DEFAULT_RECOMPUTEFREQ  0         /**< activity recomputation frequency (-1: off, 0: near boundary only) */
#define DEFAULT_USEINITIALSOL  FALSE     /**< should best known solution be used? */
#define DEFAULT_ONLYMAINSCIP   FALSE     /**< should run in subscip be skipped? */
#define DEFAULT_ONLYBEFORENODE TRUE      /**< should run before presolving be skipped? */
#define DEFAULT_ONLYWITHOUTSOL FALSE     /**< should run with solution be skipped? */

/*
 * Data structures
 */

/* variable types */
enum FJ_Vartype
{
   FJ_CONTINUOUS = 0,
   FJ_INTEGER    = 1
};
typedef enum FJ_Vartype FJ_VARTYPE;

/* index-coefficient pair */
struct FJ_IdxCoeff
{
   int                   idx;                /**< variable or constraint index */
   SCIP_Real             coeff;              /**< coefficient */
};
typedef struct FJ_IdxCoeff FJ_IDXCOEFF;

/* variable structure */
struct FJ_Var
{
   FJ_VARTYPE            vartype;            /**< variable type */
   SCIP_Real             lb;                 /**< lower bound */
   SCIP_Real             ub;                 /**< upper bound */
   SCIP_Real             objcoeff;           /**< objective coefficient */
   FJ_IDXCOEFF*          coeffs;             /**< coefficients in constraints */
   int                   ncoeffs;            /**< number of coefficients */
};
typedef struct FJ_Var FJ_VAR;

/* constraint structure */
struct FJ_Constraint
{
   SCIP_Real             lhs;                /**< left hand side, -infinity if single */
   SCIP_Real             rhs;                /**< right hand side, infinity if inactive */
   FJ_IDXCOEFF*          coeffs;             /**< coefficients */
   int                   ncoeffs;            /**< number of coefficients */
};
typedef struct FJ_Constraint FJ_CONSTRAINT;

/* move structure */
struct FJ_Move
{
   SCIP_Real             value;              /**< variable value */
   SCIP_Real             score;              /**< move score */
};
typedef struct FJ_Move FJ_MOVE;

/* shift buffer entry */
struct FJ_ShiftBuffer
{
   SCIP_Real             value;              /**< shift value */
   SCIP_Real             weight;             /**< weight */
};
typedef struct FJ_ShiftBuffer FJ_SHIFTBUFFER;

/* feasibility jump problem data */
struct FJ_Problem
{
   FJ_VAR*               vars;               /**< variables */
   int                   nvars;              /**< number of variables */
   FJ_CONSTRAINT*        constraints;        /**< constraints */
   int                   nconstraints;       /**< number of constraints */
   int                   constraintssize;    /**< size of constraints array */
   int                   nnonzeros;          /**< number of nonzeros */
   SCIP_Bool             feasibility;        /**< problem without objective? */
   SCIP_Real             locallowerbound;    /**< local lower bound */
};
typedef struct FJ_Problem FJ_PROBLEM;

/* feasibility jump solver data
 *
 * @todo Currently only jump moves are implemented (minimizing weighted violations).
 *       The design allows for additional move types (e.g., incremental up/down moves), which would be stored in
 *       separate arrays and evaluated together for variable selection.
 */
struct FJ_Solver
{
   FJ_PROBLEM*           problem;            /**< problem data */
   FJ_MOVE*              jumpmoves;          /**< jump moves for each variable */
   FJ_SHIFTBUFFER**      shiftbuffer;        /**< buffer for computing jumps (array of pointers) */
   int                   shiftbuffersize;    /**< size of shift buffer */
   int                   nshiftbuffer;       /**< number of entries in shift buffer */
   int*                  goodvarsset;        /**< indices of variables with positive score */
   int                   ngoodvars;          /**< number of good variables */
   int                   goodvarssize;       /**< size of good vars array */
   int*                  goodvarssetidx;     /**< mapping from variable to index in goodvarsset */
   SCIP_Real*            incumbentassignment;/**< incumbent assignment */
   SCIP_Real*            act;                /**< per-constraint current activity */
   int*                  unsatidxs;          /**< violated constraint indices */
   int                   nunsat;             /**< number of violated constraints */
   int                   unsatidxssize;      /**< size of violated array */
   int*                  consnmoves;         /**< per-constraint moves since the last exact activity recomputation */
   SCIP_Real*            weight;             /**< per-constraint penalty weight */
   int*                  unsatidx;           /**< per-constraint position in unsatidxs, -1 if satisfied */
   SCIP_Real             objectiveweight;    /**< weight for objective function */
   SCIP_Real             weightupdateincrement; /**< weight update increment */
   int                   nsols;              /**< number of solutions found */
   int                   previolations;      /**< number of violations at last improvement */
   SCIP_Longint          effort;             /**< effort spent since last improvement */
};
typedef struct FJ_Solver FJ_SOLVER;

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   SCIP_Real             weightupdatedecay;  /**< decay factor for the weight update */
   SCIP_Real             randprobgood;       /**< probability of random selection from good variables */
   SCIP_Real             randprobviol;       /**< probability of random selection from violated constraint */
   int                   samplesize;         /**< number of candidates to sample in tournament selection */
   int                   verbosity;          /**< verbosity level of the feasibility jump solver */
   int                   maxsols;            /**< maximum number of solutions to find (-1: unlimited) */
   SCIP_Real             mindecreaserel;     /**< minimum relative decrease to reset the effort budget */
   int                   maxeffortfac;       /**< maximum effort factor per nonzero without improvement (-1: unlimited) */
   int                   recomputefreq;      /**< activity recomputation frequency (-1: off, 0: near boundary only) */
   SCIP_Bool             useinitialsol;      /**< should best known solution be used? */
   SCIP_Bool             onlymainscip;       /**< should run in subscip be skipped? */
   SCIP_Bool             onlybeforenode;     /**< should run before presolving be skipped? */
   SCIP_Bool             onlywithoutsol;     /**< should run with solution be skipped? */
   SCIP_Longint          maxeffort;          /**< absolute effort budget without improvement (maxeffortfac * nnonzeros) */
};


/*
 * Local methods for FJ_Problem
 */

/** computes the constraint's contribution to the feasibility score */
static
SCIP_Real fjConstraintScore(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_CONSTRAINT*        constraint,         /**< constraint */
   SCIP_Real             act                 /**< activity */
   )
{
   assert(scip != NULL);
   assert(constraint != NULL);

   if( !SCIPisInfinity(scip, constraint->rhs) && act >= constraint->rhs )
      return constraint->rhs - act;

   if( !SCIPisInfinity(scip, -constraint->lhs) && act <= constraint->lhs )
      return act - constraint->lhs;

   return 0.0;
}

/** creates a problem */
static
SCIP_RETCODE fjProblemCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_PROBLEM**          problem             /**< pointer to problem */
   )
{
   FJ_PROBLEM* prob;

   assert(scip != NULL);
   assert(problem != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &prob) );
   prob->vars = NULL;
   prob->nvars = 0;
   prob->constraints = NULL;
   prob->nconstraints = 0;
   prob->constraintssize = 0;
   prob->nnonzeros = 0;
   prob->feasibility = TRUE;
   prob->locallowerbound = 0.0;

   *problem = prob;

   return SCIP_OKAY;
}

/** frees a problem */
static
SCIP_RETCODE fjProblemFree(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_PROBLEM**          problem             /**< pointer to problem */
   )
{
   FJ_PROBLEM* prob;
   int i;

   assert(scip != NULL);
   assert(problem != NULL);

   prob = *problem;
   assert(prob != NULL);

   /* free constraint data */
   if( prob->constraints != NULL )
   {
      for( i = prob->nconstraints - 1; i >= 0; --i )
      {
         if( prob->constraints[i].coeffs != NULL )
         {
            SCIPfreeBlockMemoryArray(scip, &prob->constraints[i].coeffs, prob->constraints[i].ncoeffs);
         }
      }
      SCIPfreeBlockMemoryArray(scip, &prob->constraints, prob->constraintssize);
   }

   /* free variable data */
   if( prob->vars != NULL )
   {
      for( i = prob->nvars - 1; i >= 0; --i )
      {
         if( prob->vars[i].coeffs != NULL )
         {
            SCIPfreeBlockMemoryArray(scip, &prob->vars[i].coeffs, prob->vars[i].ncoeffs);
         }
      }
      SCIPfreeBlockMemoryArray(scip, &prob->vars, prob->nvars);
   }

   SCIPfreeBlockMemory(scip, problem);

   return SCIP_OKAY;
}

/** adds a variable to the problem */
static
SCIP_RETCODE fjProblemAddVar(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_PROBLEM*           problem,            /**< problem */
   FJ_VARTYPE            vartype,            /**< variable type */
   SCIP_Real             lb,                 /**< lower bound */
   SCIP_Real             ub,                 /**< upper bound */
   SCIP_Real             objcoeff,           /**< objective coefficient */
   int*                  idx                 /**< pointer to store variable index */
   )
{
   FJ_VAR* var;

   assert(scip != NULL);
   assert(problem != NULL);
   assert(problem->vars != NULL);
   assert(idx != NULL);

   *idx = problem->nvars;
   var = problem->vars + problem->nvars;
   var->vartype = vartype;
   var->lb = lb;
   var->ub = ub;
   var->objcoeff = objcoeff;
   var->coeffs = NULL;
   var->ncoeffs = 0;

   problem->nvars++;

   return SCIP_OKAY;
}

/** adds a constraint to the problem */
static
SCIP_RETCODE fjProblemAddConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_PROBLEM*           problem,            /**< problem */
   SCIP_Real*            consvals,           /**< coefficient values */
   int*                  consinds,           /**< variable indices */
   int                   consnvars,          /**< number of variables */
   SCIP_Real             constant,           /**< row activity constant */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool             normalize,          /**< should be normalized? */
   SCIP_Bool*            success             /**< is addition successful? */
   )
{
   FJ_CONSTRAINT* constraint;
   SCIP_Real abslhs;
   SCIP_Real absrhs;
   SCIP_Real scalar;
   int i;

   assert(scip != NULL);
   assert(problem != NULL);
   assert(problem->constraints != NULL);
   assert(consvals != NULL);
   assert(consinds != NULL);
   assert(consnvars >= 1);
   assert(success != NULL);

   /* normalize the sides so that right hand side is single or non-negative */
   if( normalize )
   {
      assert(!SCIPisInfinity(scip, -lhs) || !SCIPisInfinity(scip, rhs));
      abslhs = REALABS(lhs);
      absrhs = REALABS(rhs);
      scalar = MIN(abslhs, absrhs);
      if( scalar < 1.0 )
         scalar = 1.0;

      if( SCIPisInfinity(scip, rhs) || ( !SCIPisInfinity(scip, -lhs) && rhs < 0.0 ) )
      {
         SCIPswapReals(&lhs, &rhs);
         scalar *= -1.0;
      }
   }
   else
      scalar = 1.0;

   lhs = !SCIPisInfinity(scip, scalar > 0.0 ? -lhs : lhs) ? (lhs - constant) / scalar : -SCIPinfinity(scip);
   rhs = (rhs - constant) / scalar;

   /* terminate if a normalized side exceeds infinity threshold */
   if( SCIPisInfinity(scip, lhs) || SCIPisInfinity(scip, -rhs) )
   {
      *success = FALSE;

      return SCIP_OKAY;
   }
   assert(lhs <= rhs);

   *success = TRUE;

   constraint = problem->constraints + problem->nconstraints;
   constraint->lhs = lhs;
   constraint->rhs = rhs;
   constraint->ncoeffs = consnvars;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &constraint->coeffs, consnvars) );

   /* store the normalized coefficients */
   for( i = 0; i < consnvars; ++i )
   {
      constraint->coeffs[i].idx = consinds[i];
      constraint->coeffs[i].coeff = consvals[i] / scalar;
   }

   problem->nnonzeros += consnvars;
   problem->nconstraints++;

   return SCIP_OKAY;
}

/** allocates and fills variable cross-references */
static
SCIP_RETCODE fjProblemFillReferences(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_PROBLEM*           problem,            /**< problem */
   int*                  nvarcoeffs          /**< number of coefficients per variable */
   )
{
   FJ_VAR* var;
   FJ_CONSTRAINT* constraint;
   int varidx;
   int i;
   int j;

   assert(scip != NULL);
   assert(problem != NULL);
   assert(nvarcoeffs != NULL);

   /* allocate the variable coefficients */
   for( i = 0; i < problem->nvars; ++i )
   {
      assert(problem->vars[i].coeffs == NULL);
      assert(problem->vars[i].ncoeffs == 0);

      if( nvarcoeffs[i] > 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->vars[i].coeffs, nvarcoeffs[i]) );
      }
   }

   /* fill the variable cross-references from the constraints */
   for( i = 0; i < problem->nconstraints; ++i )
   {
      constraint = problem->constraints + i;

      for( j = 0; j < constraint->ncoeffs; ++j )
      {
         varidx = constraint->coeffs[j].idx;
         var = problem->vars + varidx;
         assert(var->ncoeffs < nvarcoeffs[varidx]);
         var->coeffs[var->ncoeffs].idx = i;
         var->coeffs[var->ncoeffs].coeff = constraint->coeffs[j].coeff;
         var->ncoeffs++;
      }
   }

#ifndef NDEBUG
   for( i = 0; i < problem->nvars; ++i )
      assert(problem->vars[i].ncoeffs == nvarcoeffs[i]);
#endif

   return SCIP_OKAY;
}


/*
 * Local methods for FJ_Solver
 */

/** comparison function for shift buffer sorting */
static
SCIP_DECL_SORTPTRCOMP(fjShiftBufferComp)
{
   FJ_SHIFTBUFFER* a = (FJ_SHIFTBUFFER*)elem1;
   FJ_SHIFTBUFFER* b = (FJ_SHIFTBUFFER*)elem2;

   /* compare by value first, then by weight as tiebreaker (matches C++ std::pair behavior) */
   if( a->value < b->value )
      return -1;
   else if( a->value > b->value )
      return 1;
   else if( a->weight < b->weight )
      return -1;
   else if( a->weight > b->weight )
      return 1;
   else
      return 0;
}

/** recomputes a constraint's incumbent activity exactly from its coefficients and resets its move counter */
static
void fjSolverRecomputeConstraint(
   FJ_SOLVER*            solver,             /**< solver */
   int                   constridx           /**< constraint index */
   )
{
   FJ_CONSTRAINT* constraint;
   int i;

   assert(solver != NULL);
   assert(constridx >= 0 && constridx < solver->problem->nconstraints);

   constraint = solver->problem->constraints + constridx;
   solver->act[constridx] = 0.0;

   for( i = 0; i < constraint->ncoeffs; ++i )
      solver->act[constridx] += constraint->coeffs[i].coeff * solver->incumbentassignment[constraint->coeffs[i].idx];

   solver->consnmoves[constridx] = 0;
}

/** creates a solver */
static
SCIP_RETCODE fjSolverCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_SOLVER**           solver,             /**< pointer to solver */
   FJ_PROBLEM*           problem             /**< problem */
   )
{
   FJ_SOLVER* solv;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(problem != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &solv) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->jumpmoves, problem->nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->goodvarsset, problem->nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->goodvarssetidx, problem->nvars) );

   solv->problem = problem;

   for( i = 0; i < problem->nvars; ++i )
   {
      solv->jumpmoves[i].value = 0.0;
      solv->jumpmoves[i].score = -INFINITY;
   }

   solv->shiftbuffer = NULL;
   solv->shiftbuffersize = 0;
   solv->nshiftbuffer = 0;
   solv->ngoodvars = 0;
   solv->goodvarssize = problem->nvars;

   for( i = 0; i < problem->nvars; ++i )
      solv->goodvarssetidx[i] = -1;

   solv->incumbentassignment = NULL;
   solv->act = NULL;
   solv->unsatidxs = NULL;
   solv->nunsat = 0;
   solv->unsatidxssize = 0;
   solv->consnmoves = NULL;
   solv->weight = NULL;
   solv->unsatidx = NULL;
   solv->objectiveweight = 0.0;
   solv->weightupdateincrement = 1.0;
   solv->nsols = 0;
   solv->previolations = 0;
   solv->effort = 0;
   *solver = solv;

   return SCIP_OKAY;
}

/** resets the incumbent solution */
static
SCIP_RETCODE fjSolverResetIncumbent(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_SOLVER*            solver,             /**< solver */
   SCIP_Real*            initialvalues       /**< initial values, or NULL */
   )
{
   FJ_PROBLEM* problem;
   FJ_CONSTRAINT* constraint;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);

   problem = solver->problem;

   /* allocate arrays if needed */
   if( solver->incumbentassignment == NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->incumbentassignment, problem->nvars) );
   }
   if( solver->act == NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->act, problem->nconstraints) );
   }
   if( solver->unsatidxs == NULL )
   {
      solver->unsatidxssize = SCIPcalcMemGrowSize(scip, problem->nconstraints);
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->unsatidxs, solver->unsatidxssize) );
   }
   if( solver->consnmoves == NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->consnmoves, problem->nconstraints) );
   }
   if( solver->weight == NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->weight, problem->nconstraints) );
      i = 0;
      if( !problem->feasibility )
         solver->weight[i++] = 0.5;
      while( i < problem->nconstraints )
         solver->weight[i++] = 1.0;
   }
   if( solver->unsatidx == NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->unsatidx, problem->nconstraints) );
   }

   /* set initial values */
   if( initialvalues != NULL )
   {
      for( i = 0; i < problem->nvars; ++i )
         solver->incumbentassignment[i] = initialvalues[i];
   }
   else
   {
      for( i = 0; i < problem->nvars; ++i )
         solver->incumbentassignment[i] = !SCIPisInfinity(scip, -problem->vars[i].lb) ? problem->vars[i].lb
               : !SCIPisInfinity(scip, problem->vars[i].ub) ? problem->vars[i].ub : 0.0;
   }

   /* reset constraint activities and violated constraints list; the checked violation is doubled to classify within
    * half the feasibility tolerance, a margin for the normalization deviation
    */
   solver->nunsat = 0;
   for( i = 0; i < problem->nconstraints; ++i )
   {
      constraint = problem->constraints + i;

      fjSolverRecomputeConstraint(solver, i);

      if( SCIPisFeasNegative(scip, 2.0 * fjConstraintScore(scip, constraint, solver->act[i])) )
      {
         solver->unsatidx[i] = solver->nunsat;
         solver->unsatidxs[solver->nunsat] = i;
         solver->nunsat++;
      }
      else
      {
         solver->unsatidx[i] = -1;
      }
   }

   return SCIP_OKAY;
}

/** frees a solver */
static
SCIP_RETCODE fjSolverFree(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_SOLVER**           solver              /**< pointer to solver */
   )
{
   FJ_SOLVER* solv;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);

   solv = *solver;
   assert(solv != NULL);

   if( solv->shiftbuffer != NULL )
   {
      for( i = solv->shiftbuffersize - 1; i >= 0; --i )
         SCIPfreeBlockMemory(scip, solv->shiftbuffer + i); /*lint !e866*/
      SCIPfreeBlockMemoryArray(scip, &solv->shiftbuffer, solv->shiftbuffersize);
   }
   if( solv->unsatidx != NULL )
      SCIPfreeBlockMemoryArray(scip, &solv->unsatidx, solv->problem->nconstraints);
   if( solv->weight != NULL )
      SCIPfreeBlockMemoryArray(scip, &solv->weight, solv->problem->nconstraints);
   if( solv->consnmoves != NULL )
      SCIPfreeBlockMemoryArray(scip, &solv->consnmoves, solv->problem->nconstraints);
   if( solv->unsatidxs != NULL )
      SCIPfreeBlockMemoryArray(scip, &solv->unsatidxs, solv->unsatidxssize);
   if( solv->act != NULL )
      SCIPfreeBlockMemoryArray(scip, &solv->act, solv->problem->nconstraints);
   if( solv->incumbentassignment != NULL )
      SCIPfreeBlockMemoryArray(scip, &solv->incumbentassignment, solv->problem->nvars);
   SCIPfreeBlockMemoryArray(scip, &solv->goodvarssetidx, solv->problem->nvars);
   SCIPfreeBlockMemoryArray(scip, &solv->goodvarsset, solv->goodvarssize);
   SCIPfreeBlockMemoryArray(scip, &solv->jumpmoves, solv->problem->nvars);
   SCIPfreeBlockMemory(scip, solver);

   return SCIP_OKAY;
}


/*
 * Core operations
 */

/** updates the jump move value for a variable */
static
SCIP_RETCODE fjSolverUpdateJumpValue(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_SOLVER*            solver,             /**< solver */
   int                   varidx              /**< variable index */
   )
{
   FJ_PROBLEM* problem;
   FJ_VAR* var;
   FJ_CONSTRAINT* constraint;
   SCIP_Real value;
   SCIP_Real bestscore;
   SCIP_Real bestvalue;
   SCIP_Real currentscore;
   SCIP_Real currentvalue;
   SCIP_Real currentslope;
   SCIP_Real currentvaluestore;
   SCIP_Real currentslopestore;
   SCIP_Real residual;
   SCIP_Real validrangelb;
   SCIP_Real validrangeub;
   SCIP_Real coeff;
   int currentindexstore;
   int cstridx;
   int oldsize;
   int newsize;
   int i;
   int k;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(varidx >= 0 && varidx < solver->problem->nvars);

   problem = solver->problem;
   var = problem->vars + varidx;
   value = solver->incumbentassignment[varidx];

   solver->nshiftbuffer = 0;
   currentvalue = var->lb;
   currentslope = 0.0;

   /* compute shift buffer */
   for( i = 0; i < var->ncoeffs; ++i )
   {
      coeff = var->coeffs[i].coeff;
      cstridx = var->coeffs[i].idx;
      constraint = problem->constraints + cstridx;

      /* skip scoring of inactive objective constraint */
      if( SCIPisInfinity(scip, constraint->rhs) )
         continue;

      assert(coeff != 0.0); /*lint !e777*/
      residual = solver->act[cstridx] - coeff * value;

      if( coeff >= 0.0 )
      {
         if( !SCIPisInfinity(scip, -constraint->lhs) )
         {
            validrangelb = (constraint->lhs - residual) / coeff;

            if( var->vartype == FJ_INTEGER )
               validrangelb = round(validrangelb);

            if( SCIPisFeasNegative(scip, 2.0 * (validrangelb * coeff + residual - constraint->lhs)) )
            {
               validrangelb = var->vartype == FJ_INTEGER ? validrangelb + 1.0
                     : nextafter(validrangelb, (SCIP_Real)INFINITY);
            }
         }
         else
            validrangelb = constraint->lhs;

         if( !SCIPisInfinity(scip, constraint->rhs) )
         {
            validrangeub = (constraint->rhs - residual) / coeff;

            if( var->vartype == FJ_INTEGER )
               validrangeub = round(validrangeub);

            if( SCIPisFeasPositive(scip, 2.0 * (validrangeub * coeff + residual - constraint->rhs)) )
            {
               validrangeub = var->vartype == FJ_INTEGER ? validrangeub - 1.0
                     : nextafter(validrangeub, -(SCIP_Real)INFINITY);
            }
         }
         else
            validrangeub = constraint->rhs;
      }
      else
      {
         if( !SCIPisInfinity(scip, constraint->rhs) )
         {
            validrangelb = (constraint->rhs - residual) / coeff;

            if( var->vartype == FJ_INTEGER )
               validrangelb = round(validrangelb);

            if( SCIPisFeasPositive(scip, 2.0 * (validrangelb * coeff + residual - constraint->rhs)) )
            {
               validrangelb = var->vartype == FJ_INTEGER ? validrangelb + 1.0
                     : nextafter(validrangelb, (SCIP_Real)INFINITY);
            }
         }
         else
            validrangelb = -constraint->rhs;

         if( !SCIPisInfinity(scip, -constraint->lhs) )
         {
            validrangeub = (constraint->lhs - residual) / coeff;

            if( var->vartype == FJ_INTEGER )
               validrangeub = round(validrangeub);

            if( SCIPisFeasNegative(scip, 2.0 * (validrangeub * coeff + residual - constraint->lhs)) )
            {
               validrangeub = var->vartype == FJ_INTEGER ? validrangeub - 1.0
                     : nextafter(validrangeub, -(SCIP_Real)INFINITY);
            }
         }
         else
            validrangeub = -constraint->lhs;
      }

      SCIPdebugMsg(scip, "cstridx=%d coeff=%g weight=%g residual=%g range=[%g,%g]\n", cstridx, coeff,
         solver->weight[cstridx], residual, validrangelb, validrangeub);

      if( validrangelb > currentvalue )
      {
         currentslope -= solver->weight[cstridx];
         if( validrangelb < var->ub )
         {
            /* ensure capacity */
            if( solver->nshiftbuffer >= solver->shiftbuffersize )
            {
               oldsize = solver->shiftbuffersize;
               newsize = SCIPcalcMemGrowSize(scip, solver->nshiftbuffer + 1);
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &solver->shiftbuffer, oldsize, newsize) );
               for( k = oldsize; k < newsize; ++k )
               {
                  SCIP_CALL( SCIPallocBlockMemory(scip, solver->shiftbuffer + k) ); /*lint !e866*/
               }
               solver->shiftbuffersize = newsize;
            }
            solver->shiftbuffer[solver->nshiftbuffer]->value = validrangelb;
            solver->shiftbuffer[solver->nshiftbuffer]->weight = solver->weight[cstridx];
            solver->nshiftbuffer++;
         }
      }

      if( validrangeub <= currentvalue )
         currentslope += solver->weight[cstridx];
      else if( validrangeub < var->ub )
      {
         /* ensure capacity */
         if( solver->nshiftbuffer >= solver->shiftbuffersize )
         {
            oldsize = solver->shiftbuffersize;
            newsize = SCIPcalcMemGrowSize(scip, solver->nshiftbuffer + 1);
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &solver->shiftbuffer, oldsize, newsize) );
            for( k = oldsize; k < newsize; ++k )
            {
               SCIP_CALL( SCIPallocBlockMemory(scip, solver->shiftbuffer + k) ); /*lint !e866*/
            }
            solver->shiftbuffersize = newsize;
         }
         solver->shiftbuffer[solver->nshiftbuffer]->value = validrangeub;
         solver->shiftbuffer[solver->nshiftbuffer]->weight = solver->weight[cstridx];
         solver->nshiftbuffer++;
      }
   }

   /* add bounds to shift buffer */
   if( solver->nshiftbuffer + 2 >= solver->shiftbuffersize )
   {
      oldsize = solver->shiftbuffersize;
      newsize = SCIPcalcMemGrowSize(scip, solver->nshiftbuffer + 2);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &solver->shiftbuffer, oldsize, newsize) );
      for( k = oldsize; k < newsize; ++k )
      {
         SCIP_CALL( SCIPallocBlockMemory(scip, solver->shiftbuffer + k) ); /*lint !e866*/
      }
      solver->shiftbuffersize = newsize;
   }
   solver->shiftbuffer[solver->nshiftbuffer]->value = var->lb;
   solver->shiftbuffer[solver->nshiftbuffer]->weight = 0.0;
   solver->nshiftbuffer++;
   solver->shiftbuffer[solver->nshiftbuffer]->value = var->ub;
   solver->shiftbuffer[solver->nshiftbuffer]->weight = 0.0;
   solver->nshiftbuffer++;

   /* sort shift buffer */
   SCIPsortPtr((void**)solver->shiftbuffer, fjShiftBufferComp, solver->nshiftbuffer);

   /* find best shift */
   i = 0;

   do
   {
      currentslope += solver->shiftbuffer[i]->weight;
      currentvalue = solver->shiftbuffer[i]->value;

      SCIPdebugMsg(scip, "breakpoint i=%d value=%g slope=%g\n", i, currentvalue, currentslope);

      if( i + 1 == solver->nshiftbuffer
         || SCIPisPositive(scip, currentslope) || ( SCIPisZero(scip, currentslope)
         && ( SCIPisPositive(scip, var->objcoeff) || ( SCIPisZero(scip, var->objcoeff)
         && -solver->shiftbuffer[i]->value < solver->shiftbuffer[i + 1]->value ) ) ) )
         break;

      ++i;
   }
   while( TRUE ); /*lint !e506*/

   bestvalue = currentvalue;

   /* force significant shift */
   if( SCIPisEQ(scip, currentvalue, value) )
   {
      SCIPdebugMsg(scip, "min at incumbent, forcing shift\n");

      bestscore = INFINITY;
      currentscore = 0.0;
      currentvaluestore = currentvalue;
      currentslopestore = currentslope;
      currentindexstore = i;

      if( SCIPisGT(scip, solver->shiftbuffer[solver->nshiftbuffer - 1]->value, value) )
      {
         do
         {
            ++i;
            currentscore += (solver->shiftbuffer[i]->value - currentvalue) * currentslope;
            currentslope += solver->shiftbuffer[i]->weight;
            currentvalue = solver->shiftbuffer[i]->value;

            SCIPdebugMsg(scip, "right i=%d value=%g slope=%g score=%g\n", i, currentvalue, currentslope, currentscore);
         }
         while( SCIPisLE(scip, currentvalue, value) );

         bestscore = currentscore;
         bestvalue = currentvalue;
      }

      currentscore = 0.0;
      currentvalue = currentvaluestore;
      currentslope = currentslopestore;
      i = currentindexstore;

      if( SCIPisLT(scip, solver->shiftbuffer[0]->value, value) )
      {
         do
         {
            currentslope -= solver->shiftbuffer[i]->weight;
            --i;
            currentscore += (solver->shiftbuffer[i]->value - currentvalue) * currentslope;
            currentvalue = solver->shiftbuffer[i]->value;

            SCIPdebugMsg(scip, "left i=%d value=%g slope=%g score=%g\n", i, currentvalue, currentslope, currentscore);
         }
         while( SCIPisGE(scip, currentvalue, value) );

         if( bestscore > currentscore )
            bestvalue = currentvalue;
      }
   }

   SCIPdebugMsg(scip, "varidx=%d best=%g\n", varidx, bestvalue);

   solver->jumpmoves[varidx].value = bestvalue;

   return SCIP_OKAY;
}

/** updates whether a variable has positive score moves */
static
SCIP_RETCODE fjSolverUpdateGoodMoves(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_SOLVER*            solver,             /**< solver */
   int                   varidx              /**< variable index */
   )
{
   SCIP_Bool anygoodmoves;
   int newsize;
   int lastsetidx;
   int lastvaridx;
   int thissetidx;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(varidx >= 0 && varidx < solver->problem->nvars);

   anygoodmoves = (solver->jumpmoves[varidx].score > 0.0);

   if( anygoodmoves && solver->goodvarssetidx[varidx] == -1 )
   {
      /* became good, add to good set */
      if( solver->ngoodvars >= solver->goodvarssize )
      {
         newsize = SCIPcalcMemGrowSize(scip, solver->ngoodvars + 1);
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &solver->goodvarsset, solver->goodvarssize, newsize) );
         solver->goodvarssize = newsize;
      }
      solver->goodvarssetidx[varidx] = solver->ngoodvars;
      solver->goodvarsset[solver->ngoodvars] = varidx;
      solver->ngoodvars++;
   }
   else if( !anygoodmoves && solver->goodvarssetidx[varidx] != -1 )
   {
      /* became bad, remove from good set */
      lastsetidx = solver->ngoodvars - 1;
      lastvaridx = solver->goodvarsset[lastsetidx];
      thissetidx = solver->goodvarssetidx[varidx];

      solver->goodvarsset[thissetidx] = lastvaridx;
      solver->goodvarssetidx[lastvaridx] = thissetidx;
      solver->goodvarssetidx[varidx] = -1;
      solver->ngoodvars--;
   }

   return SCIP_OKAY;
}

/** sets a variable value and updates violated constraints */
static
SCIP_RETCODE fjSolverSetValue(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int                   varidx,             /**< variable index */
   SCIP_Real             newvalue,           /**< new value */
   SCIP_Longint*         effort              /**< pointer to effort counter */
   )
{
   FJ_PROBLEM* problem;
   FJ_CONSTRAINT* constraint;
   FJ_MOVE* move;
   SCIP_Real oldvalue;
   SCIP_Real delta;
   SCIP_Real coeff;
   SCIP_Real oldact;
   SCIP_Real newact;
   SCIP_Real newscore;
   SCIP_Real otherincumbent;
   SCIP_Real othercoeff;
   SCIP_Real oldcandidateact;
   SCIP_Real oldscoreterm;
   SCIP_Real newcandidateact;
   SCIP_Real newscoreterm;
   SCIP_Bool recompute;
   int cstridx;
   int lastunsatidx;
   int lastconstraintidx;
   int thisunsatidx;
   int otheridx;
   int i;
   int j;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(effort != NULL);

   problem = solver->problem;
   assert(varidx >= 0 && varidx < problem->nvars);
   oldvalue = solver->incumbentassignment[varidx];
   delta = newvalue - oldvalue;
   solver->incumbentassignment[varidx] = newvalue;
   newscore = SCIP_INVALID;

   /* update activities of all involved constraints */
   for( i = 0; i < problem->vars[varidx].ncoeffs; ++i )
   {
      cstridx = problem->vars[varidx].coeffs[i].idx;
      coeff = problem->vars[varidx].coeffs[i].coeff;
      constraint = problem->constraints + cstridx;
      oldact = solver->act[cstridx];
      ++solver->consnmoves[cstridx];
      recompute = heurdata->recomputefreq > 0 && solver->consnmoves[cstridx] >= heurdata->recomputefreq;

      if( !recompute )
      {
         solver->act[cstridx] = oldact + coeff * delta;
         newscore = fjConstraintScore(scip, constraint, solver->act[cstridx]);
         ++(*effort);

         if( heurdata->recomputefreq >= 0 && newscore < 0.0 && !SCIPisFeasNegative(scip, newscore) )
            recompute = TRUE;
      }

      if( recompute )
      {
         fjSolverRecomputeConstraint(solver, cstridx);
         newscore = fjConstraintScore(scip, constraint, solver->act[cstridx]);
         *effort += constraint->ncoeffs;
      }

      assert(newscore <= 0.0);

      /* skip scoring of inactive objective constraint */
      if( SCIPisInfinity(scip, constraint->rhs) )
         continue;

      newact = solver->act[cstridx];

      /* update the violated constraint list */
      if( solver->unsatidx[cstridx] == -1 && SCIPisFeasNegative(scip, 2.0 * newscore) )
      {
         /* became violated */
         solver->unsatidx[cstridx] = solver->nunsat;
         solver->unsatidxs[solver->nunsat] = cstridx;
         solver->nunsat++;
      }
      if( solver->unsatidx[cstridx] != -1 && !SCIPisFeasNegative(scip, 2.0 * newscore) )
      {
         /* became satisfied */
         lastunsatidx = solver->nunsat - 1;
         lastconstraintidx = solver->unsatidxs[lastunsatidx];
         thisunsatidx = solver->unsatidx[cstridx];

         solver->unsatidxs[thisunsatidx] = lastconstraintidx;
         solver->unsatidx[lastconstraintidx] = thisunsatidx;
         solver->unsatidx[cstridx] = -1;
         solver->nunsat--;
      }

      /* update moves for other variables in this constraint */
      *effort += constraint->ncoeffs;
      for( j = 0; j < constraint->ncoeffs; ++j )
      {
         otheridx = constraint->coeffs[j].idx;
         if( otheridx != varidx )
         {
            move = solver->jumpmoves + otheridx;
            otherincumbent = solver->incumbentassignment[otheridx];
            othercoeff = constraint->coeffs[j].coeff;
            oldcandidateact = oldact + othercoeff * (move->value - otherincumbent);
            oldscoreterm = solver->weight[cstridx] * (fjConstraintScore(scip, constraint, oldcandidateact)
                  - fjConstraintScore(scip, constraint, oldact));
            newcandidateact = newact + othercoeff * (move->value - otherincumbent);
            newscoreterm = solver->weight[cstridx] * (fjConstraintScore(scip, constraint, newcandidateact)
                  - fjConstraintScore(scip, constraint, newact));

            move->score += newscoreterm - oldscoreterm;
            SCIP_CALL( fjSolverUpdateGoodMoves(scip, solver, otheridx) );
         }
      }
   }

   return SCIP_OKAY;
}

/** resets moves for a variable */
static
SCIP_RETCODE fjSolverResetMoves(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_SOLVER*            solver,             /**< solver */
   int                   varidx              /**< variable index */
   )
{
   FJ_PROBLEM* problem;
   FJ_VAR* var;
   FJ_MOVE* move;
   FJ_CONSTRAINT* constraint;
   SCIP_Real coeff;
   SCIP_Real candidateact;
   int cstridx;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(varidx >= 0 && varidx < solver->problem->nvars);

   problem = solver->problem;
   var = problem->vars + varidx;
   move = solver->jumpmoves + varidx;

   solver->effort += var->ncoeffs;
   SCIP_CALL( fjSolverUpdateJumpValue(scip, solver, varidx) );

   move->score = solver->objectiveweight * var->objcoeff * (solver->incumbentassignment[varidx] - move->value);

   for( i = 0; i < var->ncoeffs; ++i )
   {
      cstridx = var->coeffs[i].idx;
      coeff = var->coeffs[i].coeff;
      constraint = problem->constraints + cstridx;

      /* skip scoring of inactive objective constraint */
      if( SCIPisInfinity(scip, constraint->rhs) )
         continue;

      candidateact = solver->act[cstridx] + coeff * (move->value - solver->incumbentassignment[varidx]);

      move->score += solver->weight[cstridx] * (fjConstraintScore(scip, constraint, candidateact)
            - fjConstraintScore(scip, constraint, solver->act[cstridx]));
   }

   SCIP_CALL( fjSolverUpdateGoodMoves(scip, solver, varidx) );

   return SCIP_OKAY;
}

/** initializes the solver */
static
SCIP_RETCODE fjSolverInit(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_SOLVER*            solver,             /**< solver */
   SCIP_Real*            initialvalues       /**< initial values, or NULL */
   )
{
   int i;

   assert(scip != NULL);
   assert(solver != NULL);

   SCIP_CALL( fjSolverResetIncumbent(scip, solver, initialvalues) );

   /* reset variable scores */
   solver->ngoodvars = 0;
   for( i = 0; i < solver->problem->nvars; ++i )
   {
      SCIP_CALL( fjSolverResetMoves(scip, solver, i) );
   }

   /* force objective violation to submit better initial */
   solver->previolations = MAX(solver->nunsat, 1);
   solver->effort = 0;

   return SCIP_OKAY;
}


/*
 * Weight management
 */

/** updates weights at local minimum */
static
SCIP_RETCODE fjSolverUpdateWeights(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   FJ_PROBLEM* problem;
   FJ_CONSTRAINT* constraint;
   FJ_MOVE* move;
   SCIP_Real coeff;
   SCIP_Real candidateact;
   SCIP_Real diff;
   int cstridx;
   int varidx;
   int i;
   int j;

   assert(scip != NULL);
   assert(solver != NULL);

   problem = solver->problem;

   if( solver->nunsat == 0 )
   {
      solver->objectiveweight += solver->weightupdateincrement;

      for( i = 0; i < problem->nvars; ++i )
      {
         move = solver->jumpmoves + i;
         move->score += solver->weightupdateincrement * problem->vars[i].objcoeff
               * (solver->incumbentassignment[i] - move->value);
      }
   }
   else
   {
      for( i = 0; i < solver->nunsat; ++i )
      {
         cstridx = solver->unsatidxs[i];
         constraint = problem->constraints + cstridx;

         solver->weight[cstridx] += solver->weightupdateincrement;

         for( j = 0; j < constraint->ncoeffs; ++j )
         {
            varidx = constraint->coeffs[j].idx;
            coeff = constraint->coeffs[j].coeff;
            move = solver->jumpmoves + varidx;
            candidateact = solver->act[cstridx] + coeff * (move->value
                  - solver->incumbentassignment[varidx]);
            diff = solver->weightupdateincrement * (fjConstraintScore(scip, constraint, candidateact)
                  - fjConstraintScore(scip, constraint, solver->act[cstridx]));

            move->score += diff;

            SCIP_CALL( fjSolverUpdateGoodMoves(scip, solver, varidx) );
         }
      }
   }

   solver->weightupdateincrement /= heurdata->weightupdatedecay;

   return SCIP_OKAY;
}


/*
 * Move operators
 */

/** selects a variable to modify */
static
SCIP_RETCODE fjSolverSelectVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int*                  varidx              /**< pointer to store selected variable */
   )
{
   FJ_PROBLEM* problem;
   FJ_CONSTRAINT* constraint;
   SCIP_Real randval;
   SCIP_Real bestscore;
   SCIP_Real score;
   int cstridx;
   int randomidx;
   int samplesize;
   int bestvar;
   int setidx;
   int candidatevar;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(varidx != NULL);

   problem = solver->problem;

   if( solver->ngoodvars > 0 )
   {
      randval = SCIPrandomGetReal(heurdata->randnumgen, 0.0, 1.0);

      if( randval < heurdata->randprobgood )
      {
         randomidx = SCIPrandomGetInt(heurdata->randnumgen, 0, solver->ngoodvars - 1);
         *varidx = solver->goodvarsset[randomidx];
         return SCIP_OKAY;
      }
      else
      {
         samplesize = MIN(heurdata->samplesize, solver->ngoodvars);
         bestscore = -INFINITY;
         bestvar = -1;

         for( i = 0; i < samplesize; ++i )
         {
            setidx = SCIPrandomGetInt(heurdata->randnumgen, 0, solver->ngoodvars - 1);
            candidatevar = solver->goodvarsset[setidx];
            score = solver->jumpmoves[candidatevar].score;

            if( score > bestscore )
            {
               bestscore = score;
               bestvar = candidatevar;
            }
         }

         assert(bestvar >= 0);
         *varidx = bestvar;
         return SCIP_OKAY;
      }
   }

   /* local minimum, update weights */
   SCIP_CALL( fjSolverUpdateWeights(scip, solver, heurdata) );

   if( solver->nunsat > 0 )
   {
      cstridx = solver->unsatidxs[SCIPrandomGetInt(heurdata->randnumgen, 0, solver->nunsat - 1)];
      constraint = problem->constraints + cstridx;
      randval = SCIPrandomGetReal(heurdata->randnumgen, 0.0, 1.0);

      if( randval < heurdata->randprobviol )
      {
         randomidx = SCIPrandomGetInt(heurdata->randnumgen, 0, constraint->ncoeffs - 1);
         *varidx = constraint->coeffs[randomidx].idx;
         return SCIP_OKAY;
      }
      else
      {
         bestscore = -INFINITY;
         bestvar = -1;

         for( i = 0; i < constraint->ncoeffs; ++i )
         {
            candidatevar = constraint->coeffs[i].idx;
            score = solver->jumpmoves[candidatevar].score;

            if( score > bestscore )
            {
               bestscore = score;
               bestvar = candidatevar;
            }
         }

         assert(bestvar >= 0);
         *varidx = bestvar;
         return SCIP_OKAY;
      }
   }

   /* fallback to random choice */
   *varidx = SCIPrandomGetInt(heurdata->randnumgen, 0, problem->nvars - 1);

   return SCIP_OKAY;
}

/** performs a variable move */
static
SCIP_RETCODE fjSolverDoVariableMove(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int                   varidx              /**< variable index */
   )
{
   SCIP_Real newvalue;
   SCIP_Longint effort = 0;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(varidx >= 0 && varidx < solver->problem->nvars);

   newvalue = solver->jumpmoves[varidx].value;

   SCIP_CALL( fjSolverSetValue(scip, solver, heurdata, varidx, newvalue, &effort) );
   solver->effort += effort;

   SCIP_CALL( fjSolverResetMoves(scip, solver, varidx) );

   return SCIP_OKAY;
}


/*
 * SCIP integration
 */

/** adds the objective cutoff constraint as the first constraint
 *
 *  The objective is folded into a single unnormalized constraint whose right hand side is the objective
 *  target: the floored cutoff bound for an integral objective, the cutoff bound minus the minimal improvement otherwise,
 *  or infinity without a cutoff bound so the side stays infinite until the first acceptance retightens it. An empty
 *  objective adds no constraint and leaves the local lower bound at zero so the first feasible solution counts as
 *  optimal; otherwise the local lower bound is stored for the optimality check.
 */
static
SCIP_RETCODE addObjCutoff(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_PROBLEM*           problem,            /**< problem */
   int*                  nvarcoeffs,         /**< number of coefficients per variable */
   SCIP_Bool             global              /**< is problem global? */
   )
{
   SCIP_Real* vals;
   int* inds;
   SCIP_Real cutoffbound;
   SCIP_Real rhs;
   SCIP_Bool success;
   int objnnzs;
   int i;

   assert(scip != NULL);
   assert(problem != NULL);
   assert(problem->nconstraints == 0);

   cutoffbound = SCIPgetCutoffbound(scip);

   inds = NULL;
   vals = NULL;
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, problem->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, problem->nvars) );

   /* collect the objective variables from the already extracted variable coefficients */
   objnnzs = 0;
   for( i = 0; i < problem->nvars; ++i )
   {
      if( SCIPisZero(scip, problem->vars[i].objcoeff) )
         continue;

      inds[objnnzs] = i;
      vals[objnnzs] = problem->vars[i].objcoeff;
      objnnzs++;
      ++nvarcoeffs[i];
   }

   /* add the objective cutoff for a non-empty objective */
   if( objnnzs > 0 )
   {
      problem->feasibility = FALSE;
      problem->locallowerbound = global ? SCIPgetLowerbound(scip) : SCIPgetLocalLowerbound(scip);

      if( SCIPisInfinity(scip, cutoffbound) )
         rhs = SCIPinfinity(scip);
      else if( SCIPisObjIntegral(scip) )
         rhs = floor(cutoffbound); /*lint !e644*/
      else
         rhs = cutoffbound - SCIPcutoffbounddelta(scip);

      SCIP_CALL( fjProblemAddConstraint(scip, problem, vals, inds, objnnzs, 0.0, -SCIPinfinity(scip), rhs, FALSE, &success) );
      assert(success);
   }

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &inds);

   return SCIP_OKAY;
}

/** extracts problem data before presolve */
static
SCIP_RETCODE extractProblemDataBeforePresolve(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_PROBLEM*           problem,            /**< problem */
   SCIP_Bool*            success             /**< was the problem successfully extracted? */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR** consvars;
   SCIP_CONS** conss;
   SCIP_CONS* cons;
   SCIP_CONSHDLR* conshdlr;
   SCIP_Real* consvals;
   int* consinds;
   int* nvarcoeffs;
   FJ_VARTYPE vartype;
   SCIP_Real constant;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real obj;
   int nvars;
   int nconss;
   int nlinconss;
   int consnvars;
   int requiredsize;
   int idx;
   int i;
   int j;

   assert(scip != NULL);
   assert(problem != NULL);
   assert(success != NULL);

   *success = FALSE;
   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);
   nlinconss = 0;

   conshdlr = SCIPfindConshdlr(scip, "linear");

   if( conshdlr != NULL )
      nlinconss += SCIPconshdlrGetNConss(conshdlr);

   conshdlr = SCIPfindConshdlr(scip, "setppc");

   if( conshdlr != NULL )
      nlinconss += SCIPconshdlrGetNConss(conshdlr);

   conshdlr = SCIPfindConshdlr(scip, "logicor");

   if( conshdlr != NULL )
      nlinconss += SCIPconshdlrGetNConss(conshdlr);

   conshdlr = SCIPfindConshdlr(scip, "knapsack");

   if( conshdlr != NULL )
      nlinconss += SCIPconshdlrGetNConss(conshdlr);

   conshdlr = SCIPfindConshdlr(scip, "varbound");

   if( conshdlr != NULL )
      nlinconss += SCIPconshdlrGetNConss(conshdlr);

   /* only run on linear constraint types */
   if( nlinconss < nconss )
      return SCIP_OKAY;

   vars = SCIPgetVars(scip);
   conss = SCIPgetConss(scip);
   consinds = NULL;
   consvals = NULL;
   consvars = NULL;
   SCIP_CALL( SCIPallocBufferArray(scip, &consinds, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nvarcoeffs, nvars) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->vars, nvars) );
   problem->constraintssize = nconss + 1;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->constraints, problem->constraintssize) );

   /* create feasjump solver variables */
   for( j = 0; j < nvars; ++j )
   {
      lb = SCIPvarGetLbLocal(vars[j]);
      ub = SCIPvarGetUbLocal(vars[j]);
      obj = SCIPvarGetObj(vars[j]);
      vartype = SCIPvarIsIntegral(vars[j]) ? FJ_INTEGER : FJ_CONTINUOUS;
      SCIP_CALL( fjProblemAddVar(scip, problem, vartype, lb, ub, obj, &idx) );
   }

   /* add non-empty cutoff before the constraints */
   SCIP_CALL( addObjCutoff(scip, problem, nvarcoeffs, TRUE) );

   /* create feasjump solver constraints */
   for( i = 0; i < nconss; ++i )
   {
      cons = conss[i];

      SCIP_CALL( SCIPgetConsNVars(scip, cons, &consnvars, success) );
      if( !(*success) )
         goto TERMINATE;

      if( consnvars == 0 )
         continue;

      SCIP_CALL( SCIPgetConsVals(scip, cons, consvals, nvars, success) );
      if( !(*success) )
         goto TERMINATE;

      SCIP_CALL( SCIPgetConsVars(scip, cons, consvars, nvars, success) );
      if( !(*success) )
         goto TERMINATE;

      lhs = SCIPconsGetLhs(scip, cons, success);
      if( !(*success) )
         goto TERMINATE;

      rhs = SCIPconsGetRhs(scip, cons, success);
      if( !(*success) )
         goto TERMINATE;

      constant = 0.0;
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, consvars, consvals, &consnvars, nvars, &constant, &requiredsize) );
      assert(requiredsize <= nvars);

      if( consnvars == 0 )
         continue;

      for( j = 0; j < consnvars; ++j )
      {
         assert(SCIPvarGetProbindex(consvars[j]) >= 0);
         assert(SCIPvarGetProbindex(consvars[j]) < nvars);
         consinds[j] = SCIPvarGetProbindex(consvars[j]);
         ++nvarcoeffs[consinds[j]];
      }

      SCIP_CALL( fjProblemAddConstraint(scip, problem, consvals, consinds, consnvars, constant, lhs, rhs, TRUE, success) );
      if( !(*success) )
         goto TERMINATE;
   }

   SCIP_CALL( fjProblemFillReferences(scip, problem, nvarcoeffs) );

   /* require constraint apart from objective */
   *success = problem->feasibility || problem->nconstraints > 1;

TERMINATE:
   SCIPfreeBufferArray(scip, &nvarcoeffs);
   SCIPfreeBufferArray(scip, &consvars);
   SCIPfreeBufferArray(scip, &consvals);
   SCIPfreeBufferArray(scip, &consinds);

   return SCIP_OKAY;
}

/** extracts problem data from LP */
static
SCIP_RETCODE extractProblemData(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_PROBLEM*           problem,            /**< problem */
   SCIP_COL**            cols,               /**< LP columns */
   SCIP_ROW**            rows,               /**< LP rows */
   int                   ncols,              /**< number of LP columns */
   int                   nrows,              /**< number of LP rows */
   SCIP_Bool*            success             /**< was the problem successfully extracted? */
   )
{
   SCIP_Real* vals;
   int* inds;
   int* nvarcoeffs;
   SCIP_COL* col;
   SCIP_VAR* var;
   SCIP_ROW* row;
   SCIP_COL** rowcols;
   SCIP_Real* rowvals;
   FJ_VARTYPE vartype;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real obj;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real constant;
   int idx;
   int nnonz;
   int nvals;
   int i;
   int j;

   assert(scip != NULL);
   assert(problem != NULL);
   assert(success != NULL);

   *success = TRUE;
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, ncols) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nvarcoeffs, ncols) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->vars, ncols) );
   problem->constraintssize = nrows + 1;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->constraints, problem->constraintssize) );

   /* add the columns */
   for( i = 0; i < ncols; ++i )
   {
      col = cols[i];
      var = SCIPcolGetVar(col);
      lb = SCIPcolGetLb(col);
      ub = SCIPcolGetUb(col);
      obj = SCIPcolGetObj(col);

      vartype = SCIPvarIsIntegral(var) ? FJ_INTEGER : FJ_CONTINUOUS;
      SCIP_CALL( fjProblemAddVar(scip, problem, vartype, lb, ub, obj, &idx) );
   }

   /* add non-empty cutoff before the rows */
   SCIP_CALL( addObjCutoff(scip, problem, nvarcoeffs, FALSE) );

   /* add the rows */
   for( i = 0; i < nrows; ++i )
   {
      row = rows[i];
      rowcols = SCIProwGetCols(row);
      rowvals = SCIProwGetVals(row);
      constant = SCIProwGetConstant(row);
      lhs = SCIProwGetLhs(row);
      rhs = SCIProwGetRhs(row);
      nnonz = SCIProwGetNNonz(row);
      nvals = 0;

      for( j = 0; j < nnonz; ++j )
      {
         if( SCIPcolIsInLP(rowcols[j]) )
         {
            assert(nvals < ncols);
            inds[nvals] = SCIPcolGetLPPos(rowcols[j]);
            assert(inds[nvals] >= 0);
            assert(inds[nvals] < ncols);
            vals[nvals] = rowvals[j];
            assert(vals[nvals] != 0.0); /*lint !e777*/
            ++nvarcoeffs[inds[nvals]];
            ++nvals;
         }
      }

      if( nvals == 0 )
         continue;

      SCIP_CALL( fjProblemAddConstraint(scip, problem, vals, inds, nvals, constant, lhs, rhs, TRUE, success) );
      if( !(*success) )
         goto TERMINATE;
   }

   SCIP_CALL( fjProblemFillReferences(scip, problem, nvarcoeffs) );

   /* require constraint apart from objective */
   *success = problem->feasibility || problem->nconstraints > 1;

TERMINATE:
   SCIPfreeBufferArray(scip, &nvarcoeffs);
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &inds);

   return SCIP_OKAY;
}

/** tries to add incumbent solution to SCIP */
static
SCIP_RETCODE checkIncumbentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic pointer */
   SCIP_HEURTIMING       heurtiming,         /**< heuristic timing mask */
   SCIP_Real*            incumbent,          /**< incumbent solution */
   SCIP_Bool*            stored              /**< pointer to store whether the solution was accepted */
   )
{
   SCIP_SOL* sol;
   SCIP_VAR** vars;
   SCIP_COL** cols;
   SCIP_VAR* var;
   int nvars;
   int ncols;
   int i;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(incumbent != NULL);
   assert(stored != NULL);

   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

   if( heurtiming == SCIP_HEURTIMING_BEFOREPRESOL )
   {
      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);

      for( i = 0; i < nvars; ++i )
      {
         SCIP_CALL( SCIPsetSolVal(scip, sol, vars[i], incumbent[i]) );
      }
   }
   else
   {
      SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

      for( i = 0; i < ncols; ++i )
      {
         var = SCIPcolGetVar(cols[i]);
         SCIP_CALL( SCIPsetSolVal(scip, sol, var, incumbent[i]) );
      }
   }

   SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, FALSE, FALSE, TRUE, stored) );

   return SCIP_OKAY;
}

/** reports a solution and checks the termination conditions once per iteration */
static
SCIP_RETCODE fjCheckTermination(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_HEUR*            heur,               /**< heuristic */
   SCIP_HEURTIMING       heurtiming,         /**< heuristic timing */
   SCIP_Bool*            terminate           /**< pointer to store whether to terminate */
   )
{
   FJ_PROBLEM* problem;
   FJ_CONSTRAINT* objcutoff;
   SCIP_Real objective;
   SCIP_Bool stored;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(heurdata != NULL);
   assert(heur != NULL);
   assert(terminate != NULL);

   /* check the effort budget without a significant improvement */
   if( solver->nunsat >= solver->previolations * (1.0 - heurdata->mindecreaserel) )
      *terminate = heurdata->maxeffort >= 0 && solver->effort >= heurdata->maxeffort;
   /* report a feasible solution and reset the effort state */
   else
   {
      if( solver->nunsat == 0 )
      {
         SCIP_CALL( checkIncumbentSol(scip, heur, heurtiming, solver->incumbentassignment, &stored) );
         problem = solver->problem;
         objective = !problem->feasibility ? solver->act[0] : 0.0;

         /* terminate when rejected */
         if( !stored )
         {
            if( heurdata->verbosity >= 1 )
               SCIPinfoMessage(scip, NULL, "Feasibility jump: rejected objective %g\n", objective);

            *terminate = TRUE;

            return SCIP_OKAY;
         }

         if( heurdata->verbosity >= 1 )
            SCIPinfoMessage(scip, NULL, "Feasibility jump: accepted objective %g\n", objective);

         ++solver->nsols;

         if( (heurdata->maxsols >= 0 && solver->nsols >= heurdata->maxsols)
            || SCIPisInfinity(scip, -objective)
            || SCIPisLE(scip, objective, problem->locallowerbound) )
         {
            *terminate = TRUE;

            return SCIP_OKAY;
         }

         objcutoff = problem->constraints;
         assert(SCIPisInfinity(scip, -objcutoff->lhs));
         objcutoff->rhs = SCIPisObjIntegral(scip)
               ? floor(SCIPgetCutoffbound(scip)) : SCIPgetCutoffbound(scip) - SCIPcutoffbounddelta(scip);

         /* the activity drifted while the cutoff was inactive, so recompute it exactly before reclassifying */
         fjSolverRecomputeConstraint(solver, 0);

         /* the tightened cutoff makes the incumbent violate it, so improvement re-enters as a violation */
         if( solver->unsatidx[0] == -1 && SCIPisFeasNegative(scip, 2.0 * fjConstraintScore(scip, objcutoff,
                  solver->act[0])) )
         {
            solver->unsatidx[0] = solver->nunsat;
            solver->unsatidxs[solver->nunsat] = 0;
            solver->nunsat++;

            /* refresh the moves of the objective variables against the tightened cutoff */
            for( i = 0; i < objcutoff->ncoeffs; ++i )
            {
               SCIP_CALL( fjSolverResetMoves(scip, solver, objcutoff->coeffs[i].idx) );
            }
         }
      }

      *terminate = FALSE;

      /* reset the effort budget using the retightened violation count so the baseline follows the objective cutoff
       * that the next better solution has to re-satisfy
       */
      assert(solver->nunsat > 0);
      solver->previolations = solver->nunsat;
      solver->effort = 0;
   }

   return SCIP_OKAY;
}

/** runs the feasibility jump algorithm */
static
SCIP_RETCODE runFeasjump(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic pointer */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_RESULT*          result,             /**< pointer to store the result of the heuristic call */
   SCIP_HEURTIMING       heurtiming          /**< heuristic timing mask */
   )
{
   FJ_PROBLEM* problem;
   FJ_SOLVER* solver;
   SCIP_COL** cols;
   SCIP_ROW** rows;
   SCIP_Real* initial;
   SCIP_SOL* bestsol;
   SCIP_Bool success;
   SCIP_Bool terminate;
   int ncols;
   int nrows;
   int varidx;
   int i;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( fjProblemCreate(scip, &problem) );

   if( heurtiming == SCIP_HEURTIMING_BEFOREPRESOL )
   {
      SCIP_CALL( extractProblemDataBeforePresolve(scip, problem, &success) );
   }
   else
   {
      SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
      SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
      SCIP_CALL( extractProblemData(scip, problem, cols, rows, ncols, nrows, &success) );
   }

   if( !success || problem->nvars == 0 || problem->nconstraints == 0 )
   {
      SCIP_CALL( fjProblemFree(scip, &problem) );
      return SCIP_OKAY;
   }

   /* set effort budget */
   heurdata->maxeffort = heurdata->maxeffortfac >= 0
         ? heurdata->maxeffortfac * (SCIP_Longint)problem->nnonzeros : (SCIP_Longint)-1;

   SCIP_CALL( fjSolverCreate(scip, &solver, problem) );

   /* initialize with best solution if requested */
   if( SCIPgetBestSol(scip) != NULL && heurdata->useinitialsol && heurtiming != SCIP_HEURTIMING_BEFOREPRESOL )
   {
      initial = NULL;
      bestsol = SCIPgetBestSol(scip);

      SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &initial, ncols) );

      for( i = 0; i < ncols; ++i )
      {
         initial[i] = SCIPgetSolVal(scip, bestsol, SCIPcolGetVar(cols[i]));
      }

      SCIP_CALL( fjSolverInit(scip, solver, initial) );
      SCIPfreeBufferArray(scip, &initial);
   }
   else
   {
      SCIP_CALL( fjSolverInit(scip, solver, NULL) );
   }

   /* main loop */
   terminate = FALSE;

   while( !terminate && !SCIPisStopped(scip) )
   {
      /* select and perform move */
      SCIP_CALL( fjSolverSelectVariable(scip, solver, heurdata, &varidx) );
      SCIP_CALL( fjSolverDoVariableMove(scip, solver, heurdata, varidx) );

      SCIPdebugMsg(scip, "varidx=%d value=%g nunsat=%d good=%d effort=%" SCIP_LONGINT_FORMAT " obj=%g\n", varidx,
         solver->incumbentassignment[varidx], solver->nunsat, solver->ngoodvars, solver->effort,
         !problem->feasibility ? solver->act[0] : 0.0);

      /* report the solution and check the termination conditions */
      SCIP_CALL( fjCheckTermination(scip, solver, heurdata, heur, heurtiming, &terminate) );
   }

   /* report found solution */
   if( solver->nsols > 0 )
      *result = SCIP_FOUNDSOL;

   SCIP_CALL( fjSolverFree(scip, &solver) );
   SCIP_CALL( fjProblemFree(scip, &problem) );

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyFeasjump)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   SCIP_CALL( SCIPincludeHeurFeasjump(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeFeasjump)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitFeasjump)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen, DEFAULT_RANDSEED, TRUE) );

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitFeasjump)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeRandom(scip, &heurdata->randnumgen);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecFeasjump)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_Bool cutoff;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   if( ( heurdata->onlymainscip && SCIPgetSubscipDepth(scip) > 0 )
      || ( heurdata->onlybeforenode && heurtiming == SCIP_HEURTIMING_BEFOREPRESOL )
      || ( heurdata->onlywithoutsol && SCIPgetNSols(scip) > 0 ) )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Executing Feasjump\n");

   if( heurtiming == SCIP_HEURTIMING_BEFORENODE )
   {
      if( !SCIPhasCurrentNodeLP(scip) )
         return SCIP_OKAY;

      if( !SCIPisLPConstructed(scip) )
      {
         SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

         if( cutoff )
         {
            SCIP_CALL( SCIPcutoffNode(scip, SCIPgetCurrentNode(scip)) );
            return SCIP_OKAY;
         }

         SCIP_CALL( SCIPflushLP(scip) );
      }
   }

   SCIP_CALL( runFeasjump(scip, heur, heurdata, result, heurtiming) );

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the feasjump primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurFeasjump(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   heurdata = NULL;

   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
   assert(heurdata != NULL);

   heurdata->randnumgen = NULL;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecFeasjump, heurdata) );

   assert(heur != NULL);

   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyFeasjump) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeFeasjump) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitFeasjump) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitFeasjump) );

   /* add heuristic parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/weightupdatedecay",
         "decay factor for the weight update",
         &heurdata->weightupdatedecay, FALSE, DEFAULT_WEIGHTDECAY, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/randprobgood",
         "probability of random selection from good variables",
         &heurdata->randprobgood, FALSE, DEFAULT_RANDPROBGOOD, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/randprobviol",
         "probability of random selection from violated constraint",
         &heurdata->randprobviol, FALSE, DEFAULT_RANDPROBVIOL, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/samplesize",
         "number of candidates to sample in tournament selection",
         &heurdata->samplesize, FALSE, DEFAULT_SAMPLESIZE, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/verbosity",
         "verbosity level of the feasibility jump solver",
         &heurdata->verbosity, FALSE, DEFAULT_VERBOSITY, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/maxsols",
         "maximum number of solutions to find (-1: unlimited)",
         &heurdata->maxsols, FALSE, DEFAULT_MAXSOLS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/mindecreaserel",
         "minimum relative decrease to reset the effort budget",
         &heurdata->mindecreaserel, FALSE, DEFAULT_MINDECREASEREL, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/maxeffortfac",
         "maximum effort factor per nonzero without improvement (-1: unlimited)",
         &heurdata->maxeffortfac, FALSE, DEFAULT_MAXEFFORTFAC, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/recomputefreq",
         "activity recomputation frequency (-1: off, 0: near boundary only)",
         &heurdata->recomputefreq, FALSE, DEFAULT_RECOMPUTEFREQ, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/useinitialsol",
         "should best known solution be used?",
         &heurdata->useinitialsol, FALSE, DEFAULT_USEINITIALSOL, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/onlymainscip",
         "should run in subscip be skipped?",
         &heurdata->onlymainscip, FALSE, DEFAULT_ONLYMAINSCIP, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/onlybeforenode",
         "should run before presolving be skipped?",
         &heurdata->onlybeforenode, FALSE, DEFAULT_ONLYBEFORENODE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/onlywithoutsol",
         "should run with solution be skipped?",
         &heurdata->onlywithoutsol, FALSE, DEFAULT_ONLYWITHOUTSOL, NULL, NULL) );

   return SCIP_OKAY;
}

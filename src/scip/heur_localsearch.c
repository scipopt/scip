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

/**@file   heur_localsearch.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief  Local search heuristic
 * @author Dominik Kamp
 * @author Gennesaret Tjusila
 *
 * Local search heuristic based on weighted constraint satisfaction with tabu control.
 *
 * Based on:
 *  Peng Lin, Shaowei Cai, Mengchuan Zou, Jinkun Lin.
 *  "Local-MIP: Efficient Local Search for Mixed Integer Programming"
 *  Artificial Intelligence, Volume 348, 2025, 104405.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <limits.h>
#include <string.h>
#include <math.h>

#include "blockmemshell/memory.h"
#include "scip/heur_localsearch.h"
#include "scip/pub_cons.h"
#include "scip/pub_heur.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_var.h"
#include "scip/scip_cons.h"
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

#define HEUR_NAME                  "localsearch"
#define HEUR_DESC                  "local search heuristic"
#define HEUR_DISPCHAR              SCIP_HEURDISPCHAR_LNS
#define HEUR_PRIORITY              800
#define HEUR_FREQ                  0
#define HEUR_FREQOFS               0
#define HEUR_MAXDEPTH              -1
#define HEUR_TIMING                SCIP_HEURTIMING_BEFOREPRESOL | SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP           FALSE

#define DEFAULT_RANDSEED           2832      /**< random seed */
#define DEFAULT_SMOOTHPROB         3         /**< smooth probability (out of 10000) */
#define DEFAULT_TABUBASE           3         /**< minimum tabu tenure */
#define DEFAULT_TABUVARIATION      10        /**< random tabu variation */
#define DEFAULT_SAMPLEUNSAT        12        /**< unsat constraints to sample */
#define DEFAULT_BMSUNSATINFEAS     2000      /**< max moves scored (infeasible) */
#define DEFAULT_BMSUNSATFEAS       3000      /**< max moves scored (feasible) */
#define DEFAULT_SAMPLESAT          20        /**< sat constraints to sample */
#define DEFAULT_BMSSAT             190       /**< max moves scored (sat) */
#define DEFAULT_BMSFLIP            20        /**< binary flip candidates */
#define DEFAULT_BMSRANDOM          150       /**< random move candidates */
#define DEFAULT_MAXEFFORT          100000000 /**< max effort without improvement */
#define DEFAULT_VERBOSITY          0         /**< verbosity level */
#define DEFAULT_ONLYWITHOUTSOL     FALSE     /**< only run without existing solution? */
#define DEFAULT_USEINITIALSOL      TRUE      /**< use best known solution as start? */
#define DEFAULT_BEFOREPRESOL       FALSE     /**< run before presolving? */
#define DEFAULT_MAXSOLS            1         /**< max solutions to find */
#define DEFAULT_ITERATIONS         1000      /**< iterations to check progress */
#define DEFAULT_MINPERCENTDECREASE 10        /**< minimum percent decrease required */
#define DEFAULT_CALLBACKEFFORT     1000000   /**< effort between termination checks */


/*
 * Data structures
 */

/** index-coefficient pair for cross-references */
struct LS_IdxCoeff
{
   int                   idx;                /**< variable or constraint index */
   SCIP_Real             coeff;              /**< normalized coefficient */
   int                   pos;                /**< position in the other direction's array */
};
typedef struct LS_IdxCoeff LS_IDXCOEFF;

/** variable structure (static model data) */
struct LS_Var
{
   SCIP_Real             lb;                 /**< lower bound */
   SCIP_Real             ub;                 /**< upper bound */
   SCIP_Real             obj;                /**< objective coefficient */
   LS_IDXCOEFF*          coeffs;             /**< constraints this var appears in */
   int                   ncoeffs;            /**< number of constraints */
   int                   objidx;             /**< position in objvaridxs[], -1 if not in obj */
};
typedef struct LS_Var LS_VAR;

/** constraint structure (static model data, normalized) */
struct LS_Constraint
{
   SCIP_Real             lhs;                /**< normalized left hand side, -infinity if single */
   SCIP_Real             rhs;                /**< normalized right hand side */
   SCIP_Real             mhs;                /**< normalized middle hand side, -infinity if single */
   LS_IDXCOEFF*          coeffs;             /**< normalized variable coefficients */
   int                   ncoeffs;            /**< number of coefficients */
};
typedef struct LS_Constraint LS_CONSTRAINT;

/** problem data (static model) */
struct LS_Problem
{
   LS_VAR*               vars;               /**< variables */
   int*                  varidxs;            /**< variable index map */
   int                   nbinvars;           /**< number of binary variables */
   int                   nintvars;           /**< number of integral variables */
   int                   nvars;              /**< number of variables */
   LS_CONSTRAINT*        conss;              /**< regular constraints */
   int                   nconss;             /**< number of constraints */
   int                   consssize;          /**< constraint array size */
   int*                  objvaridxs;         /**< variable indices with nonzero obj coeff */
   int                   nobjvars;           /**< number of objective variables */
   int                   nnonzeros;          /**< total number of nonzero entries */
};
typedef struct LS_Problem LS_PROBLEM;

/** solver data (dynamic solving state) */
struct LS_Solver
{
   LS_PROBLEM*           problem;            /**< problem data */
   int                   curstep;            /**< current iteration step */
   SCIP_Bool             optimality;         /**< was feasibility reached once before? */
   SCIP_Bool             iskeepfeas;         /**< currently in lift-move chain? */
   SCIP_Real             objtarget;          /**< objective constraint side, infinite if none */
   int                   objweight;          /**< objective scoring weight */
   SCIP_Real*            incumbentassignment; /**< current variable assignment */
   SCIP_Real             incumbentobjective; /**< current objective value */
   int*                  unsatidxs;          /**< violated constraint indices */
   int                   nunsat;             /**< number of violated constraints */
   SCIP_Real*            act;                /**< per-constraint current activity */
   int*                  weight;             /**< per-constraint penalty weight */
   int*                  unsatidx;           /**< per-constraint position in unsatidxs */
   int*                  allowincstep;       /**< per-variable tabu: earliest allow increase */
   int*                  allowdecstep;       /**< per-variable tabu: earliest allow decrease */
   int*                  lastincstep;        /**< per-variable step of last increase */
   int*                  lastdecstep;        /**< per-variable step of last decrease */
   int*                  neighborvaridxs;    /**< neighbor candidate var indices */
   SCIP_Real*            neighborvalues;     /**< neighbor candidate move values */
   int                   neighborsize;       /**< current number of neighbors */
   int                   neighborcap;        /**< allocated neighbor capacity */
   SCIP_Bool*            scoretable;         /**< per-variable binary dedup flag */
   int*                  scoreidxs;          /**< stack of marked score indices */
   int                   nscoreidxs;         /**< number of marked score entries */
   int*                  tempunsatidxs;      /**< temp copy for partial shuffle */
   SCIP_Bool*            sampledconstrs;     /**< per-constraint dedup flag for sat sampling */
   int*                  sampledidxs;        /**< stack of sampled constraint indices */
   int                   nsampled;           /**< number of sampled constraints */
   SCIP_Real*            liftbound;          /**< per obj-var: feasible bound in improving direction */
   SCIP_Bool*            liftvalid;          /**< are liftbounds valid */
   SCIP_Longint          subscore;           /**< secondary score (set by tightScore) */
   int                   effortatlastcallback; /**< effort at last termination check */
   int                   effortatlastimprovement; /**< effort at last improvement */
   int                   totaleffort;        /**< total effort spent */
   int                   bestviolations;     /**< best violation count seen */
   SCIP_Real             bestobjective;      /**< best objective value seen */
   int                   prevviolations;     /**< violations at last progress check */
   SCIP_Real             prevobjective;      /**< objective at last progress check */
   SCIP_Bool             isimproved;         /**< improved at last progress check? */
};
typedef struct LS_Solver LS_SOLVER;

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   int                   maxeffort;          /**< max effort without improvement */
   int                   iterations;         /**< iterations to check progress */
   int                   mindecrease;        /**< minimum percent decrease required */
   int                   verbosity;          /**< verbosity level */
   int                   maxsols;            /**< max solutions to find */
   int                   callbackeffort;     /**< effort between termination checks */
   int                   tabubase;           /**< minimum tabu tenure */
   int                   tabuvariation;      /**< random tabu variation */
   int                   smoothprob;         /**< smooth probability (out of 10000) */
   int                   sampleunsat;        /**< unsat constraints to sample */
   int                   bmsunsatinfeas;     /**< max moves scored (infeasible) */
   int                   bmsunsatfeas;       /**< max moves scored (feasible) */
   int                   samplesat;          /**< sat constraints to sample */
   int                   bmssat;             /**< max moves scored (sat) */
   int                   bmsflip;            /**< binary flip candidates */
   int                   bmsrandom;          /**< random move candidates */
   SCIP_Bool             onlywithoutsol;     /**< only run without existing solution? */
   SCIP_Bool             useinitialsol;      /**< use best known solution as start? */
   SCIP_Bool             beforepresol;       /**< run before presolving? */
};


/*
 * Local methods for LS_Problem
 */

/** creates a problem */
static
SCIP_RETCODE lsProblemCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_PROBLEM**          problem             /**< pointer to problem */
   )
{
   LS_PROBLEM* prob;

   assert(scip != NULL);
   assert(problem != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &prob) );
   BMSclearMemory(prob);

   *problem = prob;

   return SCIP_OKAY;
}

/** frees a problem */
static
SCIP_RETCODE lsProblemFree(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_PROBLEM**          problem             /**< pointer to problem */
   )
{
   LS_PROBLEM* prob;
   int i;

   assert(scip != NULL);
   assert(problem != NULL);

   prob = *problem;
   assert(prob != NULL);

   if( prob->varidxs != NULL )
      SCIPfreeBlockMemoryArray(scip, &prob->varidxs, prob->nvars);

   if( prob->objvaridxs != NULL )
      SCIPfreeBlockMemoryArray(scip, &prob->objvaridxs, prob->nobjvars);

   if( prob->conss != NULL )
   {
      for( i = prob->nconss - 1; i >= 0; --i )
      {
         if( prob->conss[i].coeffs != NULL )
            SCIPfreeBlockMemoryArray(scip, &prob->conss[i].coeffs, prob->conss[i].ncoeffs);
      }
      SCIPfreeBlockMemoryArray(scip, &prob->conss, prob->consssize);
   }

   if( prob->vars != NULL )
   {
      for( i = prob->nvars - 1; i >= 0; --i )
      {
         if( prob->vars[i].coeffs != NULL )
            SCIPfreeBlockMemoryArray(scip, &prob->vars[i].coeffs, prob->vars[i].ncoeffs);
      }
      SCIPfreeBlockMemoryArray(scip, &prob->vars, prob->nvars);
   }

   SCIPfreeBlockMemory(scip, problem);

   return SCIP_OKAY;
}

/** stores a variable at given index */
static
SCIP_RETCODE lsProblemAddVar(
   LS_PROBLEM*           problem,            /**< problem */
   int                   idx,                /**< variable index */
   SCIP_Real             lb,                 /**< lower bound */
   SCIP_Real             ub,                 /**< upper bound */
   SCIP_Real             obj                 /**< objective coefficient */
   )
{
   LS_VAR* var;

   assert(problem != NULL);
   assert(problem->vars != NULL);
   assert(problem->varidxs != NULL);

   var = problem->vars + idx;
   var->lb = lb;
   var->ub = ub;
   var->obj = obj;
   var->objidx = -1;
   var->coeffs = NULL;
   var->ncoeffs = 0;
   problem->varidxs[problem->nvars++] = idx;

   return SCIP_OKAY;
}

/** tries to normalize and add constraint to the problem with scalar max(min(|lhs|,|rhs|),1) */
static
SCIP_RETCODE lsProblemAddConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_PROBLEM*           problem,            /**< problem */
   SCIP_Real*            consvals,           /**< coefficient values */
   int*                  consinds,           /**< variable indices */
   int                   consnvars,          /**< number of variables */
   SCIP_Real             constant,           /**< row activity constant */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool*            success             /**< is addition successful? */
   )
{
   LS_CONSTRAINT* cons;
   SCIP_Real abslhs;
   SCIP_Real absrhs;
   SCIP_Real scalar;
   int i;

   assert(scip != NULL);
   assert(problem != NULL);
   assert(problem->conss != NULL);
   assert(consvals != NULL);
   assert(consinds != NULL);
   assert(consnvars >= 1);
   assert(!SCIPisInfinity(scip, -lhs) || !SCIPisInfinity(scip, rhs));
   assert(success != NULL);

   /* normalize the sides so that right hand side is single or non-negative */
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

   /* normalize coefficients and store the constraint */
   cons = problem->conss + problem->nconss;
   cons->lhs = lhs;
   cons->rhs = rhs;
   cons->mhs = !SCIPisInfinity(scip, -lhs) ? 0.5 * (lhs + rhs) : lhs;
   cons->ncoeffs = consnvars;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cons->coeffs, consnvars) );

   for( i = 0; i < consnvars; ++i )
   {
      cons->coeffs[i].idx = consinds[i];
      cons->coeffs[i].coeff = consvals[i] / scalar;
   }

   problem->nnonzeros += consnvars;
   ++problem->nconss;

   return SCIP_OKAY;
}

/** allocates and fills variable cross-references and builds the binary index */
static
SCIP_RETCODE lsProblemFillReferences(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_PROBLEM*           problem,            /**< problem */
   int*                  nvarcoeffs          /**< number of coefficients per variable */
   )
{
   LS_VAR* var;
   LS_CONSTRAINT* cons;
   int varidx;
   int i;
   int j;

   assert(scip != NULL);
   assert(problem != NULL);
   assert(nvarcoeffs != NULL);

   /* allocate variable coefficients */
   for( i = 0; i < problem->nvars; ++i )
   {
      assert(problem->vars[i].coeffs == NULL);
      assert(problem->vars[i].ncoeffs == 0);

      if( nvarcoeffs[i] > 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->vars[i].coeffs, nvarcoeffs[i]) );
      }
   }

   /* fill variable cross-references */
   for( i = 0; i < problem->nconss; ++i )
   {
      cons = problem->conss + i;

      for( j = 0; j < cons->ncoeffs; ++j )
      {
         varidx = cons->coeffs[j].idx;
         var = problem->vars + varidx;
         assert(var->ncoeffs < nvarcoeffs[varidx]);
         cons->coeffs[j].pos = var->ncoeffs;
         var->coeffs[var->ncoeffs].idx = i;
         var->coeffs[var->ncoeffs].coeff = cons->coeffs[j].coeff;
         var->coeffs[var->ncoeffs].pos = j;
         ++var->ncoeffs;
      }
   }

#ifndef NDEBUG
   for( i = 0; i < problem->nvars; ++i )
      assert(problem->vars[i].ncoeffs == nvarcoeffs[i]);
#endif

   return SCIP_OKAY;
}

/** inserts a constraint into the unsatisfied set */
static
void lsSolverInsertUnsat(
   LS_SOLVER*            solver,             /**< solver */
   int                   constridx           /**< constraint index */
   )
{
   assert(solver != NULL);
   assert(constridx >= 0);
   assert(constridx < solver->problem->nconss);
   assert(solver->unsatidx[constridx] == -1);
   assert(solver->nunsat < solver->problem->nconss);

   solver->unsatidx[constridx] = solver->nunsat;
   solver->unsatidxs[solver->nunsat++] = constridx;
}

/** removes a constraint from the unsatisfied set */
static
void lsSolverRemoveUnsat(
   LS_SOLVER*            solver,             /**< solver */
   int                   constridx           /**< constraint index */
   )
{
   int lastconstridx;
   int pos;

   assert(solver != NULL);
   assert(constridx >= 0);
   assert(constridx < solver->problem->nconss);
   assert(solver->unsatidx[constridx] >= 0);
   assert(solver->nunsat > 0);

   pos = solver->unsatidx[constridx];
   lastconstridx = solver->unsatidxs[--solver->nunsat];
   solver->unsatidxs[pos] = lastconstridx;
   solver->unsatidx[lastconstridx] = pos;
   solver->unsatidx[constridx] = -1;
}


/*
 * Local methods for LS_Solver
 */

/** recomputes incumbent objective from current assignment */
static
void lsSolverRecomputeObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver              /**< solver */
   )
{
   LS_PROBLEM* problem;
   int i;

   problem = solver->problem;

   solver->incumbentobjective = 0.0;
   for( i = 0; i < problem->nobjvars; ++i )
   {
      solver->incumbentobjective += problem->vars[problem->objvaridxs[i]].obj
            * solver->incumbentassignment[problem->objvaridxs[i]];
   }

   solver->totaleffort += problem->nobjvars;
}

/** creates and initializes a solver */
static
SCIP_RETCODE lsSolverCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER**           solverptr,          /**< pointer to solver */
   LS_PROBLEM*           problem,            /**< problem */
   SCIP_SOL*             sol,                /**< SCIP solution, or NULL for bounds */
   SCIP_COL**            cols                /**< LP columns, or NULL */
   )
{
   LS_SOLVER* solver;
   LS_CONSTRAINT* cons;
   int nvars;
   int nobjvars;
   int nnonzeros;
   int nconss;
   int varidx;
   int i;
   int j;

   assert(scip != NULL);
   assert(solverptr != NULL);
   assert(problem != NULL);

   nvars = problem->nvars;
   nobjvars = problem->nobjvars;
   nnonzeros = problem->nnonzeros + nobjvars;
   nconss = problem->nconss;

   SCIP_CALL( SCIPallocBlockMemory(scip, &solver) );

   solver->problem = problem;
   solver->curstep = 0;
   solver->optimality = FALSE;
   solver->iskeepfeas = FALSE;
   solver->objtarget = SCIPgetCutoffbound(scip) - SCIPcutoffbounddelta(scip);
   solver->objweight = 1;
   solver->effortatlastcallback = 0;
   solver->effortatlastimprovement = 0;
   solver->totaleffort = 0;
   solver->bestviolations = INT_MAX;
   solver->bestobjective = SCIPinfinity(scip);
   solver->prevviolations = INT_MAX;
   solver->prevobjective = SCIPinfinity(scip);
   solver->isimproved = TRUE;
   solver->subscore = 0;
   solver->neighborsize = 0;
   solver->nscoreidxs = 0;
   solver->nsampled = 0;

   /* allocate per-variable parallel arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->incumbentassignment, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->allowincstep, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->allowdecstep, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->lastincstep, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->lastdecstep, nvars) );

   /* allocate per-constraint parallel arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->act, nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->weight, nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->unsatidx, nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->unsatidxs, nconss) );

   /* allocate work arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->neighborvaridxs, nnonzeros) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->neighborvalues, nnonzeros) );
   solver->neighborcap = nnonzeros;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->scoretable, nvars) );
   BMSclearMemoryArray(solver->scoretable, nvars);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->scoreidxs, problem->nbinvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->tempunsatidxs, nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->sampledconstrs, nconss) );
   BMSclearMemoryArray(solver->sampledconstrs, nconss);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->sampledidxs, nconss) );
   if( nobjvars > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->liftbound, nobjvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->liftvalid, nobjvars) );
   }
   else
   {
      solver->liftbound = NULL;
      solver->liftvalid = NULL;
   }

   /* initialize variable values */
   if( sol != NULL )
   {
      assert(cols != NULL);

      for( i = 0; i < nvars; ++i )
      {
         varidx = problem->varidxs[i];
         solver->incumbentassignment[varidx] = SCIPgetSolVal(scip, sol, SCIPcolGetVar(cols[i]));
         if( solver->incumbentassignment[varidx] < problem->vars[varidx].lb )
            solver->incumbentassignment[varidx] = problem->vars[varidx].lb;
         else if( solver->incumbentassignment[varidx] > problem->vars[varidx].ub )
            solver->incumbentassignment[varidx] = problem->vars[varidx].ub;
         else if( varidx < problem->nintvars )
            solver->incumbentassignment[varidx] = SCIPfeasRound(scip, solver->incumbentassignment[varidx]);
      }
   }
   else
   {
      for( i = 0; i < nvars; ++i )
      {
         solver->incumbentassignment[i] = !SCIPisInfinity(scip, -problem->vars[i].lb) ? problem->vars[i].lb
               : !SCIPisInfinity(scip, problem->vars[i].ub) ? problem->vars[i].ub : 0.0;
      }
   }

   /* init per-variable dynamic arrays */
   BMSclearMemoryArray(solver->allowincstep, nvars);
   BMSclearMemoryArray(solver->allowdecstep, nvars);
   BMSclearMemoryArray(solver->lastincstep, nvars);
   BMSclearMemoryArray(solver->lastdecstep, nvars);

   /* init per-constraint dynamic arrays */
   solver->nunsat = 0;
   for( i = 0; i < nconss; ++i )
   {
      solver->weight[i] = 1;
      solver->unsatidx[i] = -1;
   }

   /* compute objective value */
   lsSolverRecomputeObjective(scip, solver);

   /* compute constraint activity values */
   for( i = 0; i < nconss; ++i )
   {
      cons = problem->conss + i;

      solver->act[i] = 0.0;
      for( j = 0; j < cons->ncoeffs; ++j )
         solver->act[i] += cons->coeffs[j].coeff * solver->incumbentassignment[cons->coeffs[j].idx];

      if( SCIPisFeasPositive(scip, !SCIPisInfinity(scip, -cons->mhs) && solver->act[i] <= cons->mhs
         ? cons->lhs - solver->act[i] : solver->act[i] - cons->rhs) )
         lsSolverInsertUnsat(solver, i);
   }

   solver->totaleffort += problem->nnonzeros;

   /* enter optimality mode if initial assignment feasible */
   if( solver->nunsat == 0 )
      solver->optimality = TRUE;

   *solverptr = solver;

   return SCIP_OKAY;
}

/** frees a solver */
static
SCIP_RETCODE lsSolverFree(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER**           solverptr           /**< pointer to solver */
   )
{
   LS_PROBLEM* problem;
   LS_SOLVER* solver;
   int nvars;
   int nobjvars;
   int nnonzeros;
   int nconss;

   assert(scip != NULL);
   assert(solverptr != NULL);
   solver = *solverptr;
   assert(solver != NULL);
   problem = solver->problem;
   nvars = problem->nvars;
   nobjvars = problem->nobjvars;
   nnonzeros = problem->nnonzeros + nobjvars;
   nconss = problem->nconss;

   if( nobjvars > 0 )
   {
      SCIPfreeBlockMemoryArray(scip, &solver->liftvalid, nobjvars);
      SCIPfreeBlockMemoryArray(scip, &solver->liftbound, nobjvars);
   }
   SCIPfreeBlockMemoryArray(scip, &solver->sampledidxs, nconss);
   SCIPfreeBlockMemoryArray(scip, &solver->sampledconstrs, nconss);
   SCIPfreeBlockMemoryArray(scip, &solver->tempunsatidxs, nconss);
   SCIPfreeBlockMemoryArray(scip, &solver->scoreidxs, problem->nbinvars);
   SCIPfreeBlockMemoryArray(scip, &solver->scoretable, nvars);
   SCIPfreeBlockMemoryArray(scip, &solver->neighborvalues, nnonzeros);
   SCIPfreeBlockMemoryArray(scip, &solver->neighborvaridxs, nnonzeros);
   SCIPfreeBlockMemoryArray(scip, &solver->unsatidxs, nconss);
   SCIPfreeBlockMemoryArray(scip, &solver->unsatidx, nconss);
   SCIPfreeBlockMemoryArray(scip, &solver->weight, nconss);
   SCIPfreeBlockMemoryArray(scip, &solver->act, nconss);
   SCIPfreeBlockMemoryArray(scip, &solver->lastdecstep, nvars);
   SCIPfreeBlockMemoryArray(scip, &solver->lastincstep, nvars);
   SCIPfreeBlockMemoryArray(scip, &solver->allowdecstep, nvars);
   SCIPfreeBlockMemoryArray(scip, &solver->allowincstep, nvars);
   SCIPfreeBlockMemoryArray(scip, &solver->incumbentassignment, nvars);

   SCIPfreeBlockMemory(scip, solverptr);

   return SCIP_OKAY;
}

/*
 * Core operations
 */

/** applies a variable move */
static
void lsSolverApplyMove(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int                   varidx,             /**< variable index */
   SCIP_Real             newvalue            /**< new variable value */
   )
{
   LS_PROBLEM* problem;
   LS_VAR* var;
   SCIP_Real oldvalue;
   int constridx;
   LS_CONSTRAINT* cons;
   SCIP_Real oldact;
   SCIP_Real newact;
   SCIP_Real previol;
   SCIP_Real newviol;
   int i;
   int j;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(heurdata != NULL);

   problem = solver->problem;
   var = problem->vars + varidx;

   /* update variable value */
   oldvalue = solver->incumbentassignment[varidx];
   solver->incumbentassignment[varidx] = newvalue;

   assert(solver->incumbentassignment[varidx] >= var->lb);
   assert(solver->incumbentassignment[varidx] <= var->ub);
   assert(varidx >= problem->nintvars || SCIPrealIsExactlyIntegral(solver->incumbentassignment[varidx]));

   /* recompute incumbent objective */
   if( var->objidx >= 0 )
      lsSolverRecomputeObjective(scip, solver);

   /* recompute activity for each affected constraint */
   for( i = 0; i < var->ncoeffs; ++i )
   {
      constridx = var->coeffs[i].idx;
      cons = problem->conss + constridx;
      oldact = solver->act[constridx];
      newact = 0.0;

      /* full recompute */
      for( j = 0; j < cons->ncoeffs; ++j )
         newact += cons->coeffs[j].coeff * solver->incumbentassignment[cons->coeffs[j].idx];

      previol = !SCIPisInfinity(scip, -cons->mhs) && oldact <= cons->mhs ? cons->lhs - oldact : oldact - cons->rhs;
      newviol = !SCIPisInfinity(scip, -cons->mhs) && newact <= cons->mhs ? cons->lhs - newact : newact - cons->rhs;
      solver->act[constridx] = newact;

      /* update violated constraint */
      if( !SCIPisFeasPositive(scip, previol) && SCIPisFeasPositive(scip, newviol) )
         lsSolverInsertUnsat(solver, constridx);
      /* update satisfied constraint */
      else if( SCIPisFeasPositive(scip, previol) && !SCIPisFeasPositive(scip, newviol) )
         lsSolverRemoveUnsat(solver, constridx);

      solver->totaleffort += cons->ncoeffs;
   }

   /* enter optimality mode if current assignment feasible */
   if( solver->nunsat == 0 )
      solver->optimality = TRUE;

   /* set tabu */
   if( newvalue > oldvalue )
   {
      solver->lastincstep[varidx] = solver->curstep;
      solver->allowdecstep[varidx] = solver->curstep + heurdata->tabubase
            + SCIPrandomGetInt(heurdata->randnumgen, 0, heurdata->tabuvariation);
   }
   else
   {
      solver->lastdecstep[varidx] = solver->curstep;
      solver->allowincstep[varidx] = solver->curstep + heurdata->tabubase
            + SCIPrandomGetInt(heurdata->randnumgen, 0, heurdata->tabuvariation);
   }
}

/** computes tight move value for an objective term */
static
SCIP_Real getTargetMove(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   int                   termidx             /**< term index within the objective */
   )
{
   LS_PROBLEM* problem;
   LS_VAR* var;
   SCIP_Real value;
   SCIP_Real residual;
   SCIP_Real coeff;
   SCIP_Real movevalue;
   SCIP_Bool decrease;
   int varidx;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(!SCIPisInfinity(scip, solver->objtarget));

   problem = solver->problem;
   varidx = problem->objvaridxs[termidx];
   var = problem->vars + varidx;
   coeff = var->obj;
   assert(coeff != 0.0); /*lint !e777*/
   value = solver->incumbentassignment[varidx];
   residual = solver->incumbentobjective - value * coeff;
   movevalue = (solver->objtarget - residual) / coeff;
   decrease = coeff >= 0.0;

   /* bound and round */
   if( movevalue < var->lb )
      movevalue = var->lb;
   else if( movevalue > var->ub )
      movevalue = var->ub;
   else if( varidx < problem->nintvars )
      movevalue = round(movevalue);

   /* at best bound or objective target reached */
   if( (decrease ? movevalue == var->lb : movevalue == var->ub) /*lint !e777*/
      || ( movevalue != value && SCIPisLE(scip, movevalue * coeff + residual, solver->objtarget) ) ) /*lint !e777*/
      return movevalue;

   /* improve a step */
   return varidx < problem->nintvars
         ? movevalue + (decrease ? -1.0 : 1.0)
         : nextafter(movevalue, decrease ? -(SCIP_Real)INFINITY : (SCIP_Real)INFINITY);
}

/** computes tight move value for a constraint term to the selected side */
static
SCIP_Real getTightMove(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   int                   constridx,          /**< constraint index */
   int                   termidx,            /**< term index within constraint */
   SCIP_Bool             tolhs               /**< move to left hand side instead of right hand side? */
   )
{
   LS_PROBLEM* problem;
   LS_VAR* var;
   LS_CONSTRAINT* cons;
   SCIP_Real value;
   SCIP_Real residual;
   SCIP_Real coeff;
   SCIP_Real movevalue;
   SCIP_Real moveactivity;
   SCIP_Real stepvalue;
   SCIP_Real stepactivity;
   SCIP_Bool decrease;
   int varidx;

   assert(scip != NULL);
   assert(solver != NULL);

   problem = solver->problem;
   cons = problem->conss + constridx;
   assert(!tolhs || !SCIPisInfinity(scip, -cons->lhs));
   varidx = cons->coeffs[termidx].idx;
   var = problem->vars + varidx;
   coeff = cons->coeffs[termidx].coeff;
   assert(coeff != 0.0); /*lint !e777*/
   value = solver->incumbentassignment[varidx];
   residual = solver->act[constridx] - value * coeff;
   movevalue = ((tolhs ? cons->lhs : cons->rhs) - residual) / coeff;
   decrease = tolhs == (coeff < 0.0);

   /* bound and round */
   if( movevalue < var->lb )
      movevalue = var->lb;
   else if( movevalue > var->ub )
      movevalue = var->ub;
   else if( varidx < problem->nintvars )
      movevalue = round(movevalue);

   /* at nearest bound */
   if( decrease ? movevalue == var->lb : movevalue == var->ub ) /*lint !e777*/
      return movevalue;

   /* get current activity */
   if( movevalue == value ) /*lint !e777*/
   {
      if( solver->unsatidx[constridx] == -1 )
         return movevalue;

      moveactivity = solver->act[constridx];
   }
   /* get move activity */
   else
   {
      moveactivity = movevalue * coeff + residual;

      if( tolhs ? !SCIPisFeasNegative(scip, moveactivity - cons->lhs)
                : !SCIPisFeasPositive(scip, moveactivity - cons->rhs) )
         return movevalue;
   }

   /* satisfy a step */
   stepvalue = varidx < problem->nintvars
         ? movevalue + (decrease ? -1.0 : 1.0)
         : nextafter(movevalue, decrease ? -(SCIP_Real)INFINITY : (SCIP_Real)INFINITY);

   /* single violation decreased */
   if( SCIPisInfinity(scip, -cons->lhs) )
      return stepvalue;

   /* get step activity */
   stepactivity = stepvalue * coeff + residual;

   /* step smaller violation */
   if( tolhs
      ? SCIPisLT(scip, stepactivity - cons->rhs, cons->lhs - moveactivity)
      : SCIPisLT(scip, cons->lhs - stepactivity, moveactivity - cons->rhs) )
      return stepvalue;

   /* move nearest otherwise */
   return movevalue;
}

/** computes move score */
static
SCIP_Longint getTightScore(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   int                   varidx,             /**< variable index */
   SCIP_Real             newvalue            /**< proposed move value */
   )
{
   LS_PROBLEM* problem;
   LS_VAR* var;
   SCIP_Real value;
   SCIP_Real residual;
   SCIP_Real objresidual;
   SCIP_Longint score;
   int constridx;
   SCIP_Real coeff;
   SCIP_Real act;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real mhs;
   SCIP_Real newact;
   SCIP_Real previol;
   SCIP_Real newviol;
   SCIP_Bool presat;
   SCIP_Bool nowsat;
   SCIP_Bool prestable;
   SCIP_Bool nowstable;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);

   problem = solver->problem;
   var = problem->vars + varidx;
   value = solver->incumbentassignment[varidx];
   score = 0;
   solver->subscore = 0;

   /* objective scoring in optimality mode */
   if( solver->optimality && var->objidx >= 0 )
   {
      objresidual = solver->incumbentobjective - value * var->obj;
      newact = newvalue * var->obj + objresidual;

      /* score: objective change */
      if( ( !SCIPisInfinity(scip, newact) || !SCIPisInfinity(scip, solver->incumbentobjective) )
         && ( !SCIPisInfinity(scip, -newact) || !SCIPisInfinity(scip, -solver->incumbentobjective) )
         && !SCIPisEQ(scip, newact, solver->incumbentobjective) )
      {
         if( solver->incumbentobjective > newact )
            score += solver->objweight;
         else
            score -= solver->objweight;
      }

      /* subscore: stable objective transition */
      if( !SCIPisInfinity(scip, solver->objtarget) )
      {
         prestable = SCIPisLT(scip, solver->incumbentobjective, solver->objtarget);
         nowstable = SCIPisLT(scip, newact, solver->objtarget);
      }
      else
      {
         prestable = !SCIPisInfinity(scip, solver->incumbentobjective);
         nowstable = !SCIPisInfinity(scip, newact);
      }

      if( !prestable && nowstable )
         solver->subscore += solver->objweight;
      else if( prestable && !nowstable )
         solver->subscore -= solver->objweight;

      ++solver->totaleffort;
   }

   /* regular constraint scoring */
   for( i = 0; i < var->ncoeffs; ++i )
   {
      constridx = var->coeffs[i].idx;
      coeff = var->coeffs[i].coeff;
      act = solver->act[constridx];
      lhs = problem->conss[constridx].lhs;
      rhs = problem->conss[constridx].rhs;
      mhs = problem->conss[constridx].mhs;
      residual = act - value * coeff;
      newact = newvalue * coeff + residual;
      previol = !SCIPisInfinity(scip, -mhs) && act <= mhs ? lhs - act : act - rhs;
      newviol = !SCIPisInfinity(scip, -mhs) && newact <= mhs ? lhs - newact : newact - rhs;

      /* score: feasibility transition */
      presat = !SCIPisFeasPositive(scip, previol);
      nowsat = !SCIPisFeasPositive(scip, newviol);

      if( !presat && nowsat )
         score += solver->weight[constridx];
      else if( presat && !nowsat )
         score -= solver->weight[constridx];
      else if( !presat && !nowsat && ( !SCIPisInfinity(scip, previol) || !SCIPisInfinity(scip, newviol) )
         && !SCIPisEQ(scip, previol, newviol) )
      {
         if( previol > newviol )
            score += solver->weight[constridx] / 2;
         else
            score -= solver->weight[constridx] / 2;
      }

      /* subscore: stable activity transition */
      prestable = SCIPisFeasNegative(scip, previol);
      nowstable = SCIPisFeasNegative(scip, newviol);

      if( !prestable && nowstable )
         solver->subscore += solver->weight[constridx];
      else if( prestable && !nowstable )
         solver->subscore -= solver->weight[constridx];
   }

   solver->totaleffort += var->ncoeffs;

   return score;
}


/*
 * Weight management
 */

/** updates weights of violated constraints */
static
void lsSolverUpdateWeight(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver              /**< solver */
   )
{
   int constridx;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);

   /* increase constraint weights if violated */
   for( i = 0; i < solver->nunsat; ++i )
   {
      constridx = solver->unsatidxs[i];
      ++solver->weight[constridx];
   }

   /* increase objective weight if feasible */
   if( solver->nunsat == 0 )
      ++solver->objweight;
}

/** smooths weights of satisfied constraints */
static
void lsSolverSmoothWeight(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver              /**< solver */
   )
{
   LS_PROBLEM* problem;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);

   problem = solver->problem;

   /* decrease constraint weights if feasible */
   for( i = 0; i < problem->nconss; ++i )
   {
      if( solver->unsatidx[i] == -1 && solver->weight[i] > 0 )
         --solver->weight[i];
   }

   /* decrease objective weight if better */
   if( solver->optimality && solver->objweight > 0 && !SCIPisInfinity(scip, solver->incumbentobjective)
      && SCIPisLE(scip, solver->incumbentobjective, solver->objtarget) )
      --solver->objweight;
}


/*
 * Move operators
 */

/** collects neighbor candidates from a constraint */
static
void collectConstraintNeighbors(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   int                   constridx           /**< constraint index */
   )
{
   LS_CONSTRAINT* cons;
   SCIP_Real value;
   SCIP_Real movevalue;
   SCIP_Bool tolhs;
   int varidx;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);

   /* move to nearest side if infeasible and farthest side if satisfied */
   cons = solver->problem->conss + constridx;
   tolhs = !SCIPisInfinity(scip, -cons->mhs)
         && SCIPisLE(scip, solver->act[constridx], cons->mhs) == (solver->unsatidx[constridx] >= 0);

   for( i = 0; i < cons->ncoeffs; ++i )
   {
      varidx = cons->coeffs[i].idx;
      value = solver->incumbentassignment[varidx];
      movevalue = getTightMove(scip, solver, constridx, i, tolhs);

      /* tabu check */
      if( ( movevalue < value && solver->allowdecstep[varidx] > solver->curstep )
         || ( movevalue > value && solver->allowincstep[varidx] > solver->curstep ) )
         continue;

      /* skip zero moves */
      if( SCIPisEQ(scip, movevalue, value) )
         continue;

      assert(solver->neighborsize < solver->neighborcap);
      solver->neighborvaridxs[solver->neighborsize] = varidx;
      solver->neighborvalues[solver->neighborsize] = movevalue;
      ++solver->neighborsize;
   }

   solver->totaleffort += cons->ncoeffs;
}

/** collects neighbor candidates from objective terms */
static
void collectObjectiveNeighbors(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_Bool             softtabu            /**< TRUE: use soft aspiration criterion */
   )
{
   LS_PROBLEM* problem;
   SCIP_Real value;
   SCIP_Real movevalue;
   int varidx;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);

   problem = solver->problem;

   for( i = 0; i < problem->nobjvars; ++i )
   {
      varidx = problem->objvaridxs[i];
      value = solver->incumbentassignment[varidx];
      movevalue = getTargetMove(scip, solver, i);

      /* tabu check */
      if( softtabu )
      {
         if( ( movevalue < value && solver->lastincstep[varidx] == solver->curstep - 1 )
            || ( movevalue > value && solver->lastdecstep[varidx] == solver->curstep - 1 ) )
            continue;
      }
      else
      {
         if( ( movevalue < value && solver->allowdecstep[varidx] > solver->curstep )
            || ( movevalue > value && solver->allowincstep[varidx] > solver->curstep ) )
            continue;
      }

      /* skip zero moves */
      if( SCIPisEQ(scip, movevalue, value) )
         continue;

      assert(solver->neighborsize < solver->neighborcap);
      solver->neighborvaridxs[solver->neighborsize] = varidx;
      solver->neighborvalues[solver->neighborsize] = movevalue;
      ++solver->neighborsize;
   }

   solver->totaleffort += problem->nobjvars;
}

/** partial Fisher-Yates shuffle to subsample neighbors */
static
void subsampleNeighbors(
   LS_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int                   budget              /**< max neighbors to keep */
   )
{
   int randidx;
   int tmpvar;
   SCIP_Real tmpvalue;
   int i;

   assert(solver != NULL);
   assert(heurdata != NULL);

   if( solver->neighborsize <= budget )
      return;

   for( i = 0; i < budget; ++i )
   {
      randidx = SCIPrandomGetInt(heurdata->randnumgen, 0, solver->neighborsize - 1 - i) + i;

      tmpvar = solver->neighborvaridxs[randidx];
      tmpvalue = solver->neighborvalues[randidx];
      solver->neighborvaridxs[randidx] = solver->neighborvaridxs[i];
      solver->neighborvalues[randidx] = solver->neighborvalues[i];
      solver->neighborvaridxs[i] = tmpvar;
      solver->neighborvalues[i] = tmpvalue;
   }
   solver->neighborsize = budget;
}

/** scores neighbors and selects the best move */
static
SCIP_Bool selectBestNeighbor(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_Longint          minscore,           /**< minimum acceptable score */
   int*                  bestvaridx,         /**< pointer to store best variable */
   SCIP_Real*            bestvalue           /**< pointer to store best move value */
   )
{
   LS_PROBLEM* problem;
   SCIP_Longint bestscore;
   SCIP_Longint bestsubscore;
   int varidx;
   SCIP_Real movevalue;
   SCIP_Longint score;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(bestvaridx != NULL);
   assert(bestvalue != NULL);

   problem = solver->problem;
   bestscore = minscore;
   bestsubscore = -SCIP_LONGINT_MAX;
   *bestvaridx = -1;
   *bestvalue = 0.0;

   for( i = 0; i < solver->neighborsize; ++i )
   {
      varidx = solver->neighborvaridxs[i];
      movevalue = solver->neighborvalues[i];

      /* dedup binaries */
      if( varidx < problem->nbinvars )
      {
         if( solver->scoretable[varidx] )
            continue;
         solver->scoretable[varidx] = TRUE;
         solver->scoreidxs[solver->nscoreidxs++] = varidx;
      }

      score = getTightScore(scip, solver, varidx, movevalue);

      if( score > bestscore
         || (score == bestscore && solver->subscore > bestsubscore) )
      {
         bestscore = score;
         bestsubscore = solver->subscore;
         *bestvaridx = varidx;
         *bestvalue = movevalue;
      }
   }

   return (*bestvaridx >= 0);
}

/** cleans up the score table */
static
void cleanScoreTable(
   LS_SOLVER*            solver              /**< solver */
   )
{
   int i;

   assert(solver != NULL);

   for( i = 0; i < solver->nscoreidxs; ++i )
      solver->scoretable[solver->scoreidxs[i]] = FALSE;
   solver->nscoreidxs = 0;
}

/** performs the unsat tight move */
static
SCIP_RETCODE lsSolverUnsatTightMove(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_Bool*            result              /**< is move applied? */
   )
{
   int nsample;
   int budget;
   int bestvaridx;
   SCIP_Real bestvalue;
   int randidx;
   int tmp;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);

   *result = FALSE;
   solver->neighborsize = 0;

   /* sample unsat constraints */
   if( solver->nunsat > 0 )
   {
      nsample = MIN(heurdata->sampleunsat, solver->nunsat);

      if( nsample < solver->nunsat )
      {
         /* partial Fisher-Yates shuffle on temp copy */
         BMScopyMemoryArray(solver->tempunsatidxs, solver->unsatidxs, solver->nunsat);
         for( i = 0; i < nsample; ++i )
         {
            randidx = SCIPrandomGetInt(heurdata->randnumgen, 0, solver->nunsat - 1 - i) + i;
            tmp = solver->tempunsatidxs[i];

            solver->tempunsatidxs[i] = solver->tempunsatidxs[randidx];
            solver->tempunsatidxs[randidx] = tmp;
         }

         for( i = 0; i < nsample; ++i )
            collectConstraintNeighbors(scip, solver, solver->tempunsatidxs[i]);
      }
      else
      {
         for( i = 0; i < solver->nunsat; ++i )
            collectConstraintNeighbors(scip, solver, solver->unsatidxs[i]);
      }
   }

   /* collect objective terms when the objective constraint is violated */
   if( solver->optimality && !SCIPisInfinity(scip, solver->objtarget)
      && SCIPisGT(scip, solver->incumbentobjective, solver->objtarget) )
      collectObjectiveNeighbors(scip, solver, FALSE);

   /* subsample to budget */
   budget = solver->optimality ? heurdata->bmsunsatfeas : heurdata->bmsunsatinfeas;
   subsampleNeighbors(solver, heurdata, budget);

   /* score and select */
   if( selectBestNeighbor(scip, solver, (SCIP_Longint)0, &bestvaridx, &bestvalue) )
   {
      lsSolverApplyMove(scip, solver, heurdata, bestvaridx, bestvalue);
      *result = TRUE;
   }

   return SCIP_OKAY;
}

/** performs the sat tight move */
static
SCIP_RETCODE lsSolverSatTightMove(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_Bool*            result              /**< is move applied? */
   )
{
   LS_PROBLEM* problem;
   int bestvaridx;
   SCIP_Real bestvalue;
   int constridx;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);

   *result = FALSE;
   problem = solver->problem;
   assert(problem->nconss > 0);
   solver->neighborsize = 0;

   if( !solver->optimality )
      return SCIP_OKAY;

   /* sample SAT constraints (skip duplicates, unsat, and cutoff; matching Local's sampleSet) */
   solver->nsampled = 0;
   for( i = 0; i < heurdata->samplesat; ++i )
   {
      constridx = SCIPrandomGetInt(heurdata->randnumgen, 0, problem->nconss - 1);

      /* skip already sampled */
      if( solver->sampledconstrs[constridx] )
         continue;

      /* skip UNSAT */
      if( solver->unsatidx[constridx] >= 0 )
         continue;

      /* mark as sampled */
      solver->sampledconstrs[constridx] = TRUE;
      solver->sampledidxs[solver->nsampled++] = constridx;

      collectConstraintNeighbors(scip, solver, constridx);
   }

   /* clean up sampled flags */
   for( i = 0; i < solver->nsampled; ++i )
      solver->sampledconstrs[solver->sampledidxs[i]] = FALSE;
   solver->nsampled = 0;

   /* subsample to budget */
   subsampleNeighbors(solver, heurdata, heurdata->bmssat);

   /* score and select */
   if( selectBestNeighbor(scip, solver, (SCIP_Longint)0, &bestvaridx, &bestvalue) )
   {
      lsSolverApplyMove(scip, solver, heurdata, bestvaridx, bestvalue);
      *result = TRUE;
   }

   return SCIP_OKAY;
}

/** performs the flip move */
static
SCIP_RETCODE lsSolverFlipMove(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_Bool*            result              /**< is move applied? */
   )
{
   LS_PROBLEM* problem;
   SCIP_Longint bestscore;
   SCIP_Longint bestsubscore;
   int bestvaridx;
   SCIP_Real bestvalue;
   int varidx;
   LS_VAR* var;
   SCIP_Real movevalue;
   SCIP_Longint score;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);

   *result = FALSE;
   problem = solver->problem;

   if( problem->nbinvars == 0 )
      return SCIP_OKAY;

   bestscore = 0;
   bestsubscore = -SCIP_LONGINT_MAX;
   bestvaridx = -1;
   bestvalue = 0.0;

   for( i = 0; i < heurdata->bmsflip; ++i )
   {
      varidx = SCIPrandomGetInt(heurdata->randnumgen, 0, problem->nbinvars - 1);
      var = problem->vars + varidx;

      /* dedup */
      if( solver->scoretable[varidx] )
         continue;
      solver->scoretable[varidx] = TRUE;
      solver->scoreidxs[solver->nscoreidxs++] = varidx;

      /* flip to opposite */
      if( solver->incumbentassignment[varidx] != var->lb ) /*lint !e777*/
      {
         /* check lower tabu */
         if( solver->allowdecstep[varidx] > solver->curstep )
            continue;

         movevalue = var->lb;
      }
      else
      {
         /* check upper tabu */
         if( solver->allowincstep[varidx] > solver->curstep )
            continue;

         movevalue = var->ub;
      }

      score = getTightScore(scip, solver, varidx, movevalue);

      if( score > bestscore
         || (score == bestscore && solver->subscore > bestsubscore) )
      {
         bestscore = score;
         bestsubscore = solver->subscore;
         bestvaridx = varidx;
         bestvalue = movevalue;
      }
   }

   if( bestvaridx >= 0 )
   {
      lsSolverApplyMove(scip, solver, heurdata, bestvaridx, bestvalue);
      *result = TRUE;
   }

   return SCIP_OKAY;
}

/** performs the random tight move */
static
SCIP_RETCODE lsSolverRandomTightMove(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   LS_CONSTRAINT* cons;
   SCIP_Real bestvalue;
   SCIP_Real value;
   SCIP_Real movevalue;
   SCIP_Bool tolhs;
   int bestvaridx;
   int constridx;
   int varidx;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(heurdata != NULL);

   solver->neighborsize = 0;

   /* pick 1 random unsat constraint */
   if( solver->nunsat > 0 )
   {
      constridx = solver->unsatidxs[SCIPrandomGetInt(heurdata->randnumgen, 0, solver->nunsat - 1)];
      cons = solver->problem->conss + constridx;
      tolhs = !SCIPisInfinity(scip, -cons->mhs) && solver->act[constridx] <= cons->mhs;

      /* collect terms with soft tabu */
      for( i = 0; i < cons->ncoeffs; ++i )
      {
         varidx = cons->coeffs[i].idx;
         value = solver->incumbentassignment[varidx];
         movevalue = getTightMove(scip, solver, constridx, i, tolhs);

         /* soft aspiration: only block immediate reversal */
         if( ( movevalue < value && solver->lastincstep[varidx] == solver->curstep - 1 )
            || ( movevalue > value && solver->lastdecstep[varidx] == solver->curstep - 1 ) )
            continue;

         if( SCIPisEQ(scip, movevalue, value) )
            continue;

         assert(solver->neighborsize < solver->neighborcap);
         solver->neighborvaridxs[solver->neighborsize] = varidx;
         solver->neighborvalues[solver->neighborsize] = movevalue;
         ++solver->neighborsize;
      }

      solver->totaleffort += cons->ncoeffs;
   }

   /* collect objective terms when the objective constraint is violated */
   if( solver->optimality && !SCIPisInfinity(scip, solver->objtarget)
      && SCIPisGT(scip, solver->incumbentobjective, solver->objtarget) )
      collectObjectiveNeighbors(scip, solver, TRUE);

   /* subsample */
   subsampleNeighbors(solver, heurdata, heurdata->bmsrandom);

   /* score ALL (accept negative scores) */
   if( selectBestNeighbor(scip, solver, -SCIP_LONGINT_MAX, &bestvaridx, &bestvalue) )
      lsSolverApplyMove(scip, solver, heurdata, bestvaridx, bestvalue);

   return SCIP_OKAY;
}

/** performs the lift move (feasible objective improvement) */
static
SCIP_RETCODE lsSolverLiftMove(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_Bool*            result              /**< is move applied? */
   )
{
   LS_PROBLEM* problem;
   LS_VAR* bestvar;
   LS_VAR* var;
   LS_CONSTRAINT* cons;
   SCIP_Real objresidual;
   SCIP_Real objcoeff;
   SCIP_Real objvalue;
   SCIP_Real bestobjvalue;
   SCIP_Real bestvalue;
   SCIP_Real value;
   SCIP_Real movevalue;
   SCIP_Bool tolhs;
   int objpos;
   int bestlaststep;
   int bestvaridx;
   int laststep;
   int varidx;
   int i;
   int j;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);

   *result = FALSE;
   problem = solver->problem;

   /* select best objective-improving move */
   bestobjvalue = solver->incumbentobjective;
   bestvaridx = -1;
   bestvalue = 0.0;
   bestlaststep = INT_MAX;

   for( i = 0; i < problem->nobjvars; ++i )
   {
      varidx = problem->objvaridxs[i];
      var = problem->vars + varidx;
      objcoeff = var->obj;
      value = solver->incumbentassignment[varidx];

      /* validate lift bound */
      if( !solver->iskeepfeas || !solver->liftvalid[i] )
      {
         if( objcoeff >= 0.0 )
         {
            /* at best bound */
            if( SCIPisGE(scip, var->lb, value) )
            {
               solver->liftbound[i] = value;
               solver->liftvalid[i] = TRUE;
               continue;
            }

            solver->liftbound[i] = var->lb;

            for( j = 0; j < var->ncoeffs; ++j )
            {
               cons = problem->conss + var->coeffs[j].idx;
               tolhs = var->coeffs[j].coeff > 0.0;

               /* no improvement restriction */
               if( tolhs && SCIPisInfinity(scip, -cons->lhs) )
                  continue;

               movevalue = getTightMove(scip, solver, var->coeffs[j].idx, var->coeffs[j].pos, tolhs);

               if( solver->liftbound[i] < movevalue )
               {
                  /* lift is blocked */
                  if( SCIPisGE(scip, movevalue, value) )
                  {
                     solver->liftbound[i] = value;
                     break;
                  }

                  solver->liftbound[i] = movevalue;
               }
            }
         }
         else
         {
            /* at best bound */
            if( SCIPisLE(scip, var->ub, value) )
            {
               solver->liftbound[i] = value;
               solver->liftvalid[i] = TRUE;
               continue;
            }

            solver->liftbound[i] = var->ub;

            for( j = 0; j < var->ncoeffs; ++j )
            {
               cons = problem->conss + var->coeffs[j].idx;
               tolhs = var->coeffs[j].coeff < 0.0;

               /* no improvement restriction */
               if( tolhs && SCIPisInfinity(scip, -cons->lhs) )
                  continue;

               movevalue = getTightMove(scip, solver, var->coeffs[j].idx, var->coeffs[j].pos, tolhs);

               if( solver->liftbound[i] > movevalue )
               {
                  /* lift is blocked */
                  if( SCIPisLE(scip, movevalue, value) )
                  {
                     solver->liftbound[i] = value;
                     break;
                  }

                  solver->liftbound[i] = movevalue;
               }
            }
         }

         solver->liftvalid[i] = TRUE;
         solver->totaleffort += j;
      }

      movevalue = solver->liftbound[i];

      /* skip blocked move */
      if( objcoeff >= 0.0 ? SCIPisGE(scip, movevalue, value) : SCIPisLE(scip, movevalue, value) )
         continue;

      objresidual = solver->incumbentobjective - value * objcoeff;
      objvalue = movevalue * objcoeff + objresidual;
      laststep = (movevalue < value) ? solver->lastdecstep[varidx] : solver->lastincstep[varidx];

      if( ( SCIPisInfinity(scip, -objvalue) || SCIPisInfinity(scip, bestobjvalue)
         || SCIPisLE(scip, objvalue, bestobjvalue) )
         && ( laststep < bestlaststep || ( !SCIPisInfinity(scip, objvalue) && !SCIPisInfinity(scip, -bestobjvalue)
         && SCIPisLT(scip, objvalue, bestobjvalue) ) ) )
      {
         bestobjvalue = objvalue;
         bestvaridx = varidx;
         bestvalue = movevalue;
         bestlaststep = laststep;
      }
   }

   solver->totaleffort += problem->nobjvars;

   if( bestvaridx >= 0 )
   {
      lsSolverApplyMove(scip, solver, heurdata, bestvaridx, bestvalue);
      bestvar = problem->vars + bestvaridx;

      for( i = 0; i < bestvar->ncoeffs; ++i )
      {
         cons = problem->conss + bestvar->coeffs[i].idx;

         /* invalidate other bounds */
         for( j = 0; j < cons->ncoeffs; ++j )
         {
            varidx = cons->coeffs[j].idx;

            if( varidx != bestvaridx )
            {
               objpos = problem->vars[varidx].objidx;

               if( objpos >= 0 )
                  solver->liftvalid[objpos] = FALSE;
            }
         }
      }

      solver->iskeepfeas = TRUE;
      *result = TRUE;
   }
   else
      solver->iskeepfeas = FALSE;

   return SCIP_OKAY;
}


/*
 * SCIP integration
 */

/** sets up objective data on variables */
static
SCIP_RETCODE setupObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_PROBLEM*           problem             /**< problem */
   )
{
   int nobjvars;
   int pos;
   int i;

   assert(scip != NULL);
   assert(problem != NULL);

   /* count objective variables */
   nobjvars = 0;
   for( i = 0; i < problem->nvars; ++i )
   {
      if( !SCIPisZero(scip, problem->vars[i].obj) )
         ++nobjvars;
   }

   problem->nobjvars = nobjvars;

   if( nobjvars > 0 )
   {
      pos = 0;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->objvaridxs, nobjvars) );
      for( i = 0; i < problem->nvars; ++i )
      {
         if( !SCIPisZero(scip, problem->vars[i].obj) )
         {
            problem->vars[i].objidx = pos;
            problem->objvaridxs[pos++] = i;
         }
      }
   }

   return SCIP_OKAY;
}

/** extracts problem data before presolve */
static
SCIP_RETCODE extractProblemDataBeforePresolve(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_PROBLEM*           problem,            /**< problem */
   SCIP_Bool*            success             /**< is the problem successfully extracted? */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR** consvars;
   SCIP_CONS** conss;
   SCIP_CONSHDLR* conshdlr;
   SCIP_Real* consvals;
   int* nvarcoeffs;
   int* consinds;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real constant;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int binpos;
   int intpos;
   int contpos;
   int varidx;
   int nvars;
   int nconss;
   int nlinconss;
   int consnvars;
   int requiredsize;
   int i;
   int j;

   assert(scip != NULL);
   assert(problem != NULL);
   assert(problem->vars == NULL);
   assert(problem->nvars == 0);
   assert(problem->conss == NULL);
   assert(problem->nconss == 0);
   assert(success != NULL);

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

   if( nlinconss < nconss )
   {
      *success = FALSE;

      return SCIP_OKAY;
   }

   *success = TRUE;
   vars = SCIPgetVars(scip);
   conss = SCIPgetConss(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &consinds, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nvarcoeffs, nvars) );

   /* allocate variable, index, and constraint arrays */
   if( nvars > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->vars, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->varidxs, nvars) );
   }
   if( nconss > 0 )
   {
      problem->consssize = nconss;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->conss, problem->consssize) );
   }

   /* count variables for each type */
   for( i = 0; i < nvars; ++i )
   {
      if( SCIPvarIsIntegral(vars[i]) )
      {
         lb = SCIPfeasCeil(scip, SCIPvarGetLbLocal(vars[i]));
         ub = SCIPfeasFloor(scip, SCIPvarGetUbLocal(vars[i]));

         if( ub - lb == 1.0 ) /*lint !e777*/
            ++problem->nbinvars;
         else if( lb != ub ) /*lint !e777*/
            ++problem->nintvars;
      }
   }

   /* add variables to the problem */
   problem->nintvars += problem->nbinvars;
   binpos = 0;
   intpos = problem->nbinvars;
   contpos = problem->nintvars;

   for( i = 0; i < nvars; ++i )
   {
      if( SCIPvarIsIntegral(vars[i]) )
      {
         lb = SCIPfeasCeil(scip, SCIPvarGetLbLocal(vars[i]));
         ub = SCIPfeasFloor(scip, SCIPvarGetUbLocal(vars[i]));

         if( ub - lb == 1.0 ) /*lint !e777*/
            varidx = binpos++;
         else if( lb != ub ) /*lint !e777*/
            varidx = intpos++;
         else
            varidx = contpos++;
      }
      else
      {
         lb = SCIPvarGetLbLocal(vars[i]);
         ub = SCIPvarGetUbLocal(vars[i]);

         varidx = contpos++;
      }

      SCIP_CALL( lsProblemAddVar(problem, varidx, lb, ub, SCIPvarGetObj(vars[i])) );
   }

   /* setup objective */
   SCIP_CALL( setupObjective(scip, problem) );

   /* resolve problem constraints */
   for( i = 0; i < nconss; ++i )
   {
      SCIP_CALL( SCIPgetConsNVars(scip, conss[i], &consnvars, success) );
      if( !(*success) )
         goto TERMINATE;

      if( consnvars == 0 )
         continue;

      SCIP_CALL( SCIPgetConsVals(scip, conss[i], consvals, nvars, success) );
      if( !(*success) )
         goto TERMINATE;

      SCIP_CALL( SCIPgetConsVars(scip, conss[i], consvars, nvars, success) );
      if( !(*success) )
         goto TERMINATE;

      lhs = SCIPconsGetLhs(scip, conss[i], success);
      if( !(*success) )
         goto TERMINATE;

      rhs = SCIPconsGetRhs(scip, conss[i], success);
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
         consinds[j] = problem->varidxs[SCIPvarGetProbindex(consvars[j])];
         ++nvarcoeffs[consinds[j]];
      }

      SCIP_CALL( lsProblemAddConstraint(scip, problem, consvals, consinds, consnvars, constant, lhs, rhs, success) );
      if( !(*success) )
         goto TERMINATE;
   }

   SCIP_CALL( lsProblemFillReferences(scip, problem, nvarcoeffs) );

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
   LS_PROBLEM*           problem,            /**< problem */
   SCIP_COL**            cols,               /**< LP columns */
   SCIP_ROW**            rows,               /**< LP rows */
   int                   ncols,              /**< number of LP columns */
   int                   nrows,              /**< number of LP rows */
   SCIP_Bool*            success             /**< is the problem successfully extracted? */
   )
{
   SCIP_COL** rowcols;
   SCIP_COL* col;
   SCIP_ROW* row;
   SCIP_Real* vals;
   SCIP_Real* rowvals;
   int* nvarcoeffs;
   int* inds;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real constant;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int binpos;
   int intpos;
   int contpos;
   int varidx;
   int nnonz;
   int nvals;
   int i;
   int j;

   assert(scip != NULL);
   assert(problem != NULL);
   assert(problem->vars == NULL);
   assert(problem->nvars == 0);
   assert(problem->conss == NULL);
   assert(problem->nconss == 0);
   assert(success != NULL);

   *success = TRUE;
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, ncols) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nvarcoeffs, ncols) );

   /* allocate variable, index, and constraint arrays */
   if( ncols > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->vars, ncols) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->varidxs, ncols) );
   }
   if( nrows > 0 )
   {
      problem->consssize = nrows;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->conss, problem->consssize) );
   }

   /* count variables for each type */
   for( i = 0; i < ncols; ++i )
   {
      col = cols[i];

      if( SCIPvarIsIntegral(SCIPcolGetVar(col)) )
      {
         lb = SCIPfeasCeil(scip, SCIPcolGetLb(col));
         ub = SCIPfeasFloor(scip, SCIPcolGetUb(col));

         if( ub - lb == 1.0 ) /*lint !e777*/
            ++problem->nbinvars;
         else if( lb != ub ) /*lint !e777*/
            ++problem->nintvars;
      }
   }

   /* add variables to the problem */
   problem->nintvars += problem->nbinvars;
   binpos = 0;
   intpos = problem->nbinvars;
   contpos = problem->nintvars;

   for( i = 0; i < ncols; ++i )
   {
      col = cols[i];

      if( SCIPvarIsIntegral(SCIPcolGetVar(col)) )
      {
         lb = SCIPfeasCeil(scip, SCIPcolGetLb(col));
         ub = SCIPfeasFloor(scip, SCIPcolGetUb(col));

         if( ub - lb == 1.0 ) /*lint !e777*/
            varidx = binpos++;
         else if( lb != ub ) /*lint !e777*/
            varidx = intpos++;
         else
            varidx = contpos++;
      }
      else
      {
         lb = SCIPcolGetLb(col);
         ub = SCIPcolGetUb(col);

         varidx = contpos++;
      }

      SCIP_CALL( lsProblemAddVar(problem, varidx, lb, ub, SCIPcolGetObj(col)) );
   }

   /* setup objective */
   SCIP_CALL( setupObjective(scip, problem) );

   /* resolve problem constraints */
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
            assert(SCIPcolGetLPPos(rowcols[j]) >= 0);
            assert(SCIPcolGetLPPos(rowcols[j]) < ncols);
            inds[nvals] = problem->varidxs[SCIPcolGetLPPos(rowcols[j])];
            vals[nvals] = rowvals[j];
            ++nvarcoeffs[inds[nvals]];
            ++nvals;
         }
      }

      if( nvals == 0 )
         continue;

      SCIP_CALL( lsProblemAddConstraint(scip, problem, vals, inds, nvals, constant, lhs, rhs, success) );
      if( !(*success) )
         goto TERMINATE;
   }

   SCIP_CALL( lsProblemFillReferences(scip, problem, nvarcoeffs) );

TERMINATE:
   SCIPfreeBufferArray(scip, &nvarcoeffs);
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &inds);

   return SCIP_OKAY;
}

/** checks incumbent solution against SCIP */
static
SCIP_RETCODE checkIncumbentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic pointer */
   SCIP_HEURTIMING       heurtiming,         /**< heuristic timing mask */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   LS_PROBLEM* problem;
   SCIP_VAR** vars;
   SCIP_COL** cols;
   SCIP_SOL* sol;
   SCIP_Bool stored;
   int nvars;
   int ncols;
   int i;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(solver != NULL);
   assert(result != NULL);

   problem = solver->problem;
   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

   if( heurtiming == SCIP_HEURTIMING_BEFOREPRESOL )
   {
      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);

      for( i = 0; i < nvars; ++i )
      {
         SCIP_CALL( SCIPsetSolVal(scip, sol, vars[i], solver->incumbentassignment[problem->varidxs[i]]) );
      }
   }
   else
   {
      SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

      for( i = 0; i < ncols; ++i )
      {
         SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPcolGetVar(cols[i]), solver->incumbentassignment[problem->varidxs[i]]) );
      }
   }

   SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, FALSE, TRUE, &stored) );

   if( stored )
   {
      SCIPdebugMsg(scip, "Local search found feasible solution\n");
      *result = SCIP_FOUNDSOL;
   }

   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   return SCIP_OKAY;
}

/** checks for solution reporting and termination conditions */
static
SCIP_RETCODE lsCheckTermination(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_HEUR*            heur,               /**< heuristic */
   SCIP_HEURTIMING       heurtiming,         /**< heuristic timing */
   SCIP_Real*            solution,           /**< solution to report, or NULL */
   int*                  nsols,              /**< pointer to number of solutions found */
   SCIP_RESULT*          result,             /**< pointer to store result */
   SCIP_Bool*            terminate           /**< pointer to store termination flag */
   )
{
   assert(scip != NULL);
   assert(solver != NULL);
   assert(heurdata != NULL);
   assert(heur != NULL);
   assert(nsols != NULL);
   assert(result != NULL);
   assert(terminate != NULL);

   *terminate = FALSE;

   if( solution != NULL || solver->totaleffort - solver->effortatlastcallback > heurdata->callbackeffort )
   {
      SCIP_Bool quitnumsol;
      SCIP_Bool quiteffort;
      SCIP_Bool quitunbounded;
      int nfoundsols;

      nfoundsols = *nsols;
      solver->effortatlastcallback = solver->totaleffort;

      if( solution != NULL )
      {
         SCIP_CALL( checkIncumbentSol(scip, heur, heurtiming, solver, result) );

         if( *result == SCIP_FOUNDSOL )
         {
            solver->objtarget = SCIPgetCutoffbound(scip) - SCIPcutoffbounddelta(scip);
            ++nfoundsols;

            if( heurdata->verbosity >= 1 )
            {
               SCIPinfoMessage(scip, NULL, "Local search: found solution #%d obj=%g\n",
                  nfoundsols, solver->incumbentobjective);
            }
         }
         else
            *terminate = TRUE;
      }

      if( !(*terminate) )
      {
         quitnumsol = (nfoundsols >= heurdata->maxsols);
         if( quitnumsol )
         {
            SCIPdebugMsg(scip, "Local: quitting, solutions %d >= %d\n", nfoundsols, heurdata->maxsols);
         }

         quiteffort = (!solver->isimproved
               || solver->totaleffort - solver->effortatlastimprovement > heurdata->maxeffort);
         if( quiteffort )
         {
            SCIPdebugMsg(scip, "Local: quitting, no improvement\n");
         }

         quitunbounded = (solver->nunsat == 0 && SCIPisInfinity(scip, -solver->incumbentobjective));
         if( quitunbounded )
         {
            SCIPdebugMsg(scip, "Local: quitting, unbounded\n");
         }

         *terminate = quitnumsol || quiteffort || quitunbounded;
      }

      *nsols = nfoundsols;
   }

   return SCIP_OKAY;
}

/** runs the local search algorithm */
static
SCIP_RETCODE runLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic pointer */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_RESULT*          result,             /**< pointer to store the result */
   SCIP_HEURTIMING       heurtiming          /**< heuristic timing mask */
   )
{
   LS_PROBLEM* problem;
   LS_SOLVER* solver;
   SCIP_COL** cols = NULL;
   SCIP_ROW** rows = NULL;
   int ncols = 0;
   int nrows = 0;
   int nsols;
   SCIP_Bool terminate;
   SCIP_Bool success;
   SCIP_Real* solution;
   SCIP_SOL* bestsol;
   SCIP_Bool liftresult;
   SCIP_Bool moveresult;
   int violationsdecrease;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( lsProblemCreate(scip, &problem) );

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

   if( !success || problem->nvars == 0 || problem->nconss == 0 )
   {
      SCIP_CALL( lsProblemFree(scip, &problem) );
      return SCIP_OKAY;
   }

   /* initialize solution */
   bestsol = (heurdata->useinitialsol && heurtiming != SCIP_HEURTIMING_BEFOREPRESOL) ? SCIPgetBestSol(scip) : NULL;

   /* create solver */
   SCIP_CALL( lsSolverCreate(scip, &solver, problem, bestsol, cols) );

   if( heurdata->verbosity >= 1 )
   {
      SCIPinfoMessage(scip, NULL, "Local search: starting. nvars=%d nconss=%d nnonzeros=%d\n",
         problem->nvars, problem->nconss, problem->nnonzeros);
   }

   /* main loop */
   nsols = 0;
   terminate = FALSE;

   while( !terminate && !SCIPisStopped(scip) )
   {
      solution = NULL;

      /* feasible phase */
      if( solver->nunsat == 0 )
      {
         if( !SCIPisInfinity(scip, solver->incumbentobjective)
            && SCIPisLE(scip, solver->incumbentobjective, solver->objtarget) )
         {
            solver->effortatlastimprovement = solver->totaleffort;
            solver->isimproved = TRUE;
            solution = solver->incumbentassignment;
         }

         SCIP_CALL( lsCheckTermination(scip, solver, heurdata, heur, heurtiming, solution, &nsols, result,
               &terminate) );
         if( terminate )
            break;

         liftresult = FALSE;

         SCIP_CALL( lsSolverLiftMove(scip, solver, heurdata, &liftresult) );
         ++solver->curstep;

         SCIPdebugMsg(scip, "Local: step=%d nunsat=%d effort=%d obj=%g phase=lift\n",
               solver->curstep, solver->nunsat, solver->totaleffort, solver->incumbentobjective);

         if( liftresult )
            continue;
      }

      /* infeasible phase */
      moveresult = FALSE;

      SCIP_CALL( lsSolverUnsatTightMove(scip, solver, heurdata, &moveresult) );

      if( !moveresult )
      {
         SCIP_CALL( lsSolverSatTightMove(scip, solver, heurdata, &moveresult) );
      }

      if( !moveresult )
      {
         SCIP_CALL( lsSolverFlipMove(scip, solver, heurdata, &moveresult) );
      }

      if( !moveresult )
      {
         cleanScoreTable(solver);

         if( SCIPrandomGetInt(heurdata->randnumgen, 0, 9999) > heurdata->smoothprob )
            lsSolverUpdateWeight(scip, solver);
         else
            lsSolverSmoothWeight(scip, solver);

         SCIP_CALL( lsSolverRandomTightMove(scip, solver, heurdata) );
      }

      cleanScoreTable(solver);

      /* track improvement */
      if( solver->nunsat < solver->bestviolations
         || ( solver->nunsat == solver->bestviolations
         && !SCIPisInfinity(scip, solver->incumbentobjective) && !SCIPisInfinity(scip, -solver->bestobjective)
         && SCIPisLT(scip, solver->incumbentobjective, solver->bestobjective) ) )
      {
         solver->effortatlastimprovement = solver->totaleffort;
         solver->bestviolations = solver->nunsat;
         solver->bestobjective = solver->incumbentobjective;
      }

      /* track progress */
      if( solver->curstep == 0 )
      {
         solver->prevviolations = solver->nunsat;
         solver->prevobjective = solver->incumbentobjective;
      }
      else if( solver->curstep % heurdata->iterations == 0 )
      {
         violationsdecrease = solver->prevviolations - solver->nunsat;
         if( solver->prevviolations == 0 || heurdata->mindecrease * solver->prevviolations <= 100 * violationsdecrease )
            solver->isimproved = TRUE;
         else if( violationsdecrease < 0 )
            solver->isimproved = FALSE;
         else
            solver->isimproved = (!SCIPisInfinity(scip, solver->incumbentobjective)
                  && !SCIPisInfinity(scip, -solver->prevobjective)
                  && SCIPisLT(scip, solver->incumbentobjective, solver->prevobjective));
         solver->prevviolations = solver->nunsat;
         solver->prevobjective = solver->incumbentobjective;
      }

      ++solver->curstep;

      SCIPdebugMsg(scip, "Local: step=%d nunsat=%d effort=%d obj=%g phase=search\n",
            solver->curstep, solver->nunsat, solver->totaleffort, solver->incumbentobjective);

      SCIP_CALL( lsCheckTermination(scip, solver, heurdata, heur, heurtiming, NULL, &nsols, result,
            &terminate) );
   }

   SCIP_CALL( lsSolverFree(scip, &solver) );
   SCIP_CALL( lsProblemFree(scip, &problem) );

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyLocalsearch)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   SCIP_CALL( SCIPincludeHeurLocalsearch(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeLocalsearch)
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
SCIP_DECL_HEURINIT(heurInitLocalsearch)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen, DEFAULT_RANDSEED, TRUE) );

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitLocalsearch)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeRandom(scip, &heurdata->randnumgen);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecLocalsearch)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_Bool cutoff;

   heurdata = SCIPheurGetData(heur);

   assert(result != NULL);
   assert(heurdata != NULL);

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMsg(scip, "Executing Local search\n");

   if( !heurdata->beforepresol && heurtiming == SCIP_HEURTIMING_BEFOREPRESOL )
      return SCIP_OKAY;

   if( heurtiming == SCIP_HEURTIMING_BEFORENODE )
   {
      if( SCIPgetBestSol(scip) != NULL && heurdata->onlywithoutsol )
         return SCIP_OKAY;

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

   SCIP_CALL( runLocal(scip, heur, heurdata, result, heurtiming) );

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the localsearch primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurLocalsearch(
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
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecLocalsearch, heurdata) );

   assert(heur != NULL);

   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyLocalsearch) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeLocalsearch) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitLocalsearch) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitLocalsearch) );

   /* add parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/maxeffort",
         "maximum effort spent in local search",
         &heurdata->maxeffort, FALSE, DEFAULT_MAXEFFORT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/iterations",
         "number of iterations to check progress",
         &heurdata->iterations, FALSE, DEFAULT_ITERATIONS, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/mindecrease",
         "minimum percentage decrease for violated constraints",
         &heurdata->mindecrease, FALSE, DEFAULT_MINPERCENTDECREASE, 0, 100, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/verbosity",
         "verbosity level of the local search solver",
         &heurdata->verbosity, FALSE, DEFAULT_VERBOSITY, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/maxsols",
         "maximum number of solutions to find",
         &heurdata->maxsols, FALSE, DEFAULT_MAXSOLS, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/callbackeffort",
         "effort between callbacks",
         &heurdata->callbackeffort, FALSE, DEFAULT_CALLBACKEFFORT, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/tabubase",
         "minimum tabu tenure",
         &heurdata->tabubase, FALSE, DEFAULT_TABUBASE, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/tabuvariation",
         "random tabu variation",
         &heurdata->tabuvariation, FALSE, DEFAULT_TABUVARIATION, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/smoothprob",
         "smooth probability (out of 10000)",
         &heurdata->smoothprob, FALSE, DEFAULT_SMOOTHPROB, 0, 10000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/sampleunsat",
         "unsat constraints to sample",
         &heurdata->sampleunsat, FALSE, DEFAULT_SAMPLEUNSAT, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/bmsunsatinfeas",
         "max moves scored during infeasible phase",
         &heurdata->bmsunsatinfeas, FALSE, DEFAULT_BMSUNSATINFEAS, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/bmsunsatfeas",
         "max moves scored during feasible phase",
         &heurdata->bmsunsatfeas, FALSE, DEFAULT_BMSUNSATFEAS, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/samplesat",
         "sat constraints to sample",
         &heurdata->samplesat, FALSE, DEFAULT_SAMPLESAT, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/bmssat",
         "max moves scored from sat constraints",
         &heurdata->bmssat, FALSE, DEFAULT_BMSSAT, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/bmsflip",
         "binary flip candidates",
         &heurdata->bmsflip, FALSE, DEFAULT_BMSFLIP, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/bmsrandom",
         "random move candidates",
         &heurdata->bmsrandom, FALSE, DEFAULT_BMSRANDOM, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/onlywithoutsol",
         "should local search be called only if no solution exists?",
         &heurdata->onlywithoutsol, FALSE, DEFAULT_ONLYWITHOUTSOL, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/beforepresol",
         "should the heuristic be called before presolving?",
         &heurdata->beforepresol, FALSE, DEFAULT_BEFOREPRESOL, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/useinitialsol",
         "should the heuristic use the best known solution as initialization?",
         &heurdata->useinitialsol, FALSE, DEFAULT_USEINITIALSOL, NULL, NULL) );

   return SCIP_OKAY;
}

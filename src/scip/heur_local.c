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

/**@file   heur_local.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief  Local search heuristic
 * @author Dominik Kamp
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
#include "scip/heur_local.h"
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

#define HEUR_NAME                  "local"
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

/** variable types */
enum LS_Vartype
{
   LS_BINARY     = 0,
   LS_INTEGER    = 1,
   LS_CONTINUOUS = 2
};
typedef enum LS_Vartype LS_VARTYPE;

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
   LS_VARTYPE            vartype;            /**< variable type */
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
   SCIP_Real             rhs;                /**< normalized right hand side */
   LS_IDXCOEFF*          coeffs;             /**< normalized variable coefficients */
   int                   ncoeffs;            /**< number of coefficients */
};
typedef struct LS_Constraint LS_CONSTRAINT;

/** problem data (static model) */
struct LS_Problem
{
   LS_VAR*               vars;               /**< variables */
   int                   nvars;              /**< number of variables */
   LS_CONSTRAINT*        conss;              /**< regular constraints */
   int                   nconss;             /**< number of constraints */
   int*                  objvaridxs;         /**< variable indices with nonzero obj coeff */
   int                   nobjvars;           /**< number of objective variables */
   int*                  binaryidxs;         /**< binary variable indices */
   int                   nbinary;            /**< number of binary variables */
   int                   nnonzeros;          /**< total number of nonzero entries */
};
typedef struct LS_Problem LS_PROBLEM;

/** solver data (dynamic solving state) */
struct LS_Solver
{
   LS_PROBLEM*           problem;            /**< problem data */
   int                   curstep;            /**< current iteration step */
   SCIP_Bool             iskeepfeas;         /**< currently in lift-move chain? */
   SCIP_Real             objcutoff;          /**< objective cutoff bound */
   int                   objweight;          /**< objective scoring weight */
   SCIP_Real*            incumbentassignment; /**< current variable assignment */
   SCIP_Real             incumbentobjective; /**< current objective value */
   int*                  unsatidxs;          /**< violated constraint indices */
   int                   nunsat;             /**< number of violated constraints */
   SCIP_Real*            lhs;                /**< per-constraint current LHS */
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
   SCIP_Real*            liftlowervalue;     /**< per obj-var: feasible lower value */
   SCIP_Real*            liftuppervalue;     /**< per obj-var: feasible upper value */
   SCIP_Bool*            affectedset;        /**< per-variable affected flag */
   int*                  affectedvars;       /**< stack of affected var indices */
   int                   naffected;          /**< number of affected variables */
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

   if( prob->binaryidxs != NULL )
      SCIPfreeBlockMemoryArray(scip, &prob->binaryidxs, prob->nbinary);

   if( prob->objvaridxs != NULL )
      SCIPfreeBlockMemoryArray(scip, &prob->objvaridxs, prob->nobjvars);

   if( prob->conss != NULL )
   {
      for( i = prob->nconss - 1; i >= 0; --i )
      {
         if( prob->conss[i].coeffs != NULL )
            SCIPfreeBlockMemoryArray(scip, &prob->conss[i].coeffs, prob->conss[i].ncoeffs);
      }
      SCIPfreeBlockMemoryArray(scip, &prob->conss, prob->nconss);
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

/** adds a variable to the problem */
static
SCIP_RETCODE lsProblemAddVar(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_PROBLEM*           problem,            /**< problem */
   LS_VARTYPE            vartype,            /**< variable type */
   SCIP_Real             lb,                 /**< lower bound */
   SCIP_Real             ub,                 /**< upper bound */
   SCIP_Real             obj                 /**< objective coefficient */
   )
{
   LS_VAR* var;

   assert(scip != NULL);
   assert(problem != NULL);
   assert(problem->vars != NULL);

   var = problem->vars + problem->nvars;
   var->vartype = vartype;
   var->lb = lb;
   var->ub = ub;
   var->obj = obj;
   var->coeffs = NULL;
   var->ncoeffs = 0;
   var->objidx = -1;

   ++problem->nvars;

   return SCIP_OKAY;
}

/** adds a constraint to the problem (normalizes by max(|rhs|,1)) */
static
SCIP_RETCODE lsProblemAddConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_PROBLEM*           problem,            /**< problem */
   SCIP_Real             constant,           /**< row activity constant */
   SCIP_Real             rhs,                /**< right hand side */
   int                   ncoeffs,            /**< number of coefficients */
   int*                  varidxs,            /**< variable indices */
   SCIP_Real*            coeffs           /**< coefficient values */
   )
{
   LS_VAR* var;
   LS_CONSTRAINT* cons;
   SCIP_Real scalar;
   SCIP_Real normcoeff;
   int i;

   assert(scip != NULL);
   assert(problem != NULL);
   assert(problem->conss != NULL);
   assert(varidxs != NULL || ncoeffs == 0);
   assert(coeffs != NULL || ncoeffs == 0);

   /* skip empty constraints */
   if( ncoeffs == 0 )
      return SCIP_OKAY;

   scalar = MAX3(rhs, -rhs, 1.0);
   cons = problem->conss + problem->nconss;
   cons->rhs = (rhs - constant) / scalar;
   cons->ncoeffs = 0;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cons->coeffs, ncoeffs) );

   for( i = 0; i < ncoeffs; ++i )
   {
      var = problem->vars + varidxs[i];
      assert(var->coeffs != NULL);
      normcoeff = coeffs[i] / scalar;

      /* set cross-positions: constraint entry points to var entry and vice versa */
      cons->coeffs[cons->ncoeffs].idx = varidxs[i];
      cons->coeffs[cons->ncoeffs].coeff = normcoeff;
      cons->coeffs[cons->ncoeffs].pos = var->ncoeffs;

      var->coeffs[var->ncoeffs].idx = problem->nconss;
      var->coeffs[var->ncoeffs].coeff = normcoeff;
      var->coeffs[var->ncoeffs].pos = cons->ncoeffs;

      ++cons->ncoeffs;
      ++var->ncoeffs;
   }

   problem->nnonzeros += cons->ncoeffs;
   ++problem->nconss;

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

/** creates a solver */
static
SCIP_RETCODE lsSolverCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER**           solver,             /**< pointer to solver */
   LS_PROBLEM*           problem             /**< problem */
   )
{
   LS_SOLVER* solv;
   int nvars;
   int nobjvars;
   int nnonzeros;
   int nconss;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(problem != NULL);

   nvars = problem->nvars;
   nobjvars = problem->nobjvars;
   nnonzeros = problem->nnonzeros + nobjvars;
   nconss = problem->nconss;

   SCIP_CALL( SCIPallocBlockMemory(scip, &solv) );

   solv->problem = problem;
   solv->curstep = 0;
   solv->iskeepfeas = FALSE;
   solv->objcutoff = SCIPgetCutoffbound(scip);
   solv->objweight = 1;
   solv->incumbentassignment = NULL;
   solv->incumbentobjective = 0.0;
   solv->nunsat = 0;
   solv->effortatlastcallback = 0;
   solv->effortatlastimprovement = 0;
   solv->totaleffort = 0;
   solv->bestviolations = INT_MAX;
   solv->bestobjective = SCIPinfinity(scip);
   solv->prevviolations = INT_MAX;
   solv->prevobjective = SCIPinfinity(scip);
   solv->isimproved = TRUE;
   solv->subscore = 0;
   solv->naffected = 0;
   solv->neighborsize = 0;
   solv->nscoreidxs = 0;
   solv->nsampled = 0;

   /* allocate per-variable parallel arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->allowincstep, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->allowdecstep, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->lastincstep, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->lastdecstep, nvars) );

   /* allocate per-constraint parallel arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->lhs, nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->weight, nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->unsatidx, nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->unsatidxs, nconss) );

   /* allocate work arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->neighborvaridxs, nnonzeros) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->neighborvalues, nnonzeros) );
   solv->neighborcap = nnonzeros;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->scoretable, nvars) );
   BMSclearMemoryArray(solv->scoretable, nvars);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->scoreidxs, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->tempunsatidxs, nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->sampledconstrs, nconss) );
   BMSclearMemoryArray(solv->sampledconstrs, nconss);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->sampledidxs, nconss) );
   if( nobjvars > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->liftlowervalue, nobjvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->liftuppervalue, nobjvars) );
   }
   else
   {
      solv->liftlowervalue = NULL;
      solv->liftuppervalue = NULL;
   }
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->affectedset, nvars) );
   BMSclearMemoryArray(solv->affectedset, nvars);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->affectedvars, nvars) );

   *solver = solv;

   return SCIP_OKAY;
}

/** frees a solver */
static
SCIP_RETCODE lsSolverFree(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER**           solver              /**< pointer to solver */
   )
{
   LS_SOLVER* solv;
   int nvars;
   int nobjvars;
   int nnonzeros;
   int nconss;

   assert(scip != NULL);
   assert(solver != NULL);

   solv = *solver;
   assert(solv != NULL);

   nvars = solv->problem->nvars;
   nobjvars = solv->problem->nobjvars;
   nnonzeros = solv->problem->nnonzeros + nobjvars;
   nconss = solv->problem->nconss;

   SCIPfreeBlockMemoryArray(scip, &solv->affectedvars, nvars);
   SCIPfreeBlockMemoryArray(scip, &solv->affectedset, nvars);
   if( nobjvars > 0 )
   {
      SCIPfreeBlockMemoryArray(scip, &solv->liftuppervalue, nobjvars);
      SCIPfreeBlockMemoryArray(scip, &solv->liftlowervalue, nobjvars);
   }
   SCIPfreeBlockMemoryArray(scip, &solv->sampledidxs, nconss);
   SCIPfreeBlockMemoryArray(scip, &solv->sampledconstrs, nconss);
   SCIPfreeBlockMemoryArray(scip, &solv->tempunsatidxs, nconss);
   SCIPfreeBlockMemoryArray(scip, &solv->scoreidxs, nvars);
   SCIPfreeBlockMemoryArray(scip, &solv->scoretable, nvars);
   SCIPfreeBlockMemoryArray(scip, &solv->neighborvalues, nnonzeros);
   SCIPfreeBlockMemoryArray(scip, &solv->neighborvaridxs, nnonzeros);
   SCIPfreeBlockMemoryArray(scip, &solv->unsatidxs, nconss);
   SCIPfreeBlockMemoryArray(scip, &solv->unsatidx, nconss);
   SCIPfreeBlockMemoryArray(scip, &solv->weight, nconss);
   SCIPfreeBlockMemoryArray(scip, &solv->lhs, nconss);
   SCIPfreeBlockMemoryArray(scip, &solv->lastdecstep, nvars);
   SCIPfreeBlockMemoryArray(scip, &solv->lastincstep, nvars);
   SCIPfreeBlockMemoryArray(scip, &solv->allowdecstep, nvars);
   SCIPfreeBlockMemoryArray(scip, &solv->allowincstep, nvars);
   if( solv->incumbentassignment != NULL )
      SCIPfreeBlockMemoryArray(scip, &solv->incumbentassignment, nvars);

   SCIPfreeBlockMemory(scip, solver);

   return SCIP_OKAY;
}

/** initializes variable values */
static
SCIP_RETCODE lsSolverInitSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_SOL*             sol,                /**< SCIP solution, or NULL for bounds */
   SCIP_COL**            cols,               /**< LP columns, or NULL */
   int                   ncols               /**< number of LP columns, or 0 */
   )
{
   LS_PROBLEM* problem;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);

   problem = solver->problem;

   if( solver->incumbentassignment == NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->incumbentassignment, problem->nvars) );
   }

   if( sol != NULL )
   {
      assert(cols != NULL);
      assert(ncols == problem->nvars);

      for( i = 0; i < ncols; ++i )
      {
         solver->incumbentassignment[i] = SCIPgetSolVal(scip, sol, SCIPcolGetVar(cols[i]));
         if( problem->vars[i].vartype != LS_CONTINUOUS )
            solver->incumbentassignment[i] = round(solver->incumbentassignment[i]);
         if( solver->incumbentassignment[i] < problem->vars[i].lb )
            solver->incumbentassignment[i] = problem->vars[i].lb;
         else if( solver->incumbentassignment[i] > problem->vars[i].ub )
            solver->incumbentassignment[i] = problem->vars[i].ub;
      }
   }
   else
   {
      for( i = 0; i < problem->nvars; ++i )
      {
         solver->incumbentassignment[i] = !SCIPisInfinity(scip, -problem->vars[i].lb) ? problem->vars[i].lb
               : !SCIPisInfinity(scip, problem->vars[i].ub) ? problem->vars[i].ub : 0.0;
      }
   }

   return SCIP_OKAY;
}

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
      solver->incumbentobjective += problem->vars[problem->objvaridxs[i]].obj
            * solver->incumbentassignment[problem->objvaridxs[i]];

   solver->totaleffort += problem->nobjvars;
}

/** initializes constraint LHS values and violation tracking */
static
SCIP_RETCODE lsSolverInitState(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver              /**< solver */
   )
{
   LS_PROBLEM* problem;
   LS_CONSTRAINT* cons;
   SCIP_Real violation;
   int i;
   int j;

   assert(scip != NULL);
   assert(solver != NULL);

   problem = solver->problem;

   /* init per-variable dynamic arrays */
   BMSclearMemoryArray(solver->allowincstep, problem->nvars);
   BMSclearMemoryArray(solver->allowdecstep, problem->nvars);
   BMSclearMemoryArray(solver->lastincstep, problem->nvars);
   BMSclearMemoryArray(solver->lastdecstep, problem->nvars);

   /* init per-constraint dynamic arrays */
   solver->nunsat = 0;
   for( i = 0; i < problem->nconss; ++i )
   {
      solver->weight[i] = 1;
      solver->unsatidx[i] = -1;
   }

   /* compute constraint LHS values */
   for( i = 0; i < problem->nconss; ++i )
   {
      cons = problem->conss + i;

      solver->lhs[i] = 0.0;
      for( j = 0; j < cons->ncoeffs; ++j )
      {
         solver->lhs[i] += cons->coeffs[j].coeff * solver->incumbentassignment[cons->coeffs[j].idx];
      }

      violation = solver->lhs[i] - cons->rhs;
      if( SCIPisFeasPositive(scip, violation) )
         lsSolverInsertUnsat(solver, i);
   }

   solver->totaleffort += problem->nnonzeros;

   /* compute objective value */
   lsSolverRecomputeObjective(scip, solver);

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
   SCIP_Real oldlhs;
   SCIP_Real newlhs;
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
   assert(var->vartype == LS_CONTINUOUS || SCIPrealIsExactlyIntegral(solver->incumbentassignment[varidx]));

   /* recompute objective */
   if( !SCIPisInfinity(scip, solver->objcutoff) && var->objidx >= 0 )
      lsSolverRecomputeObjective(scip, solver);

   /* recompute LHS for each affected constraint */
   for( i = 0; i < var->ncoeffs; ++i )
   {
      constridx = var->coeffs[i].idx;
      cons = problem->conss + constridx;
      oldlhs = solver->lhs[constridx];
      newlhs = 0.0;

      /* full recompute */
      for( j = 0; j < cons->ncoeffs; ++j )
      {
         newlhs += cons->coeffs[j].coeff * solver->incumbentassignment[cons->coeffs[j].idx];
      }

      previol = oldlhs - cons->rhs;
      newviol = newlhs - cons->rhs;
      solver->lhs[constridx] = newlhs;

      /* update violation tracking */
      if( !SCIPisFeasPositive(scip, previol) && SCIPisFeasPositive(scip, newviol) )
      {
         /* became violated */
         lsSolverInsertUnsat(solver, constridx);
      }
      else if( SCIPisFeasPositive(scip, previol) && !SCIPisFeasPositive(scip, newviol) )
      {
         /* became satisfied */
         lsSolverRemoveUnsat(solver, constridx);
      }

      solver->totaleffort += cons->ncoeffs;
   }

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

/** computes tight move value for a constraint term */
static
SCIP_Bool lsSolverTightMove(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   int                   constridx,          /**< constraint index */
   int                   termidx,            /**< term index within constraint */
   SCIP_Real*            result              /**< pointer to store move value */
   )
{
   LS_PROBLEM* problem;
   LS_VAR* var;
   LS_CONSTRAINT* cons;
   SCIP_Real value;
   SCIP_Real residual;
   SCIP_Real coeff;
   SCIP_Real movevalue;
   int varidx;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(result != NULL);

   problem = solver->problem;
   cons = problem->conss + constridx;
   coeff = cons->coeffs[termidx].coeff;
   varidx = cons->coeffs[termidx].idx;
   var = problem->vars + varidx;

   assert(coeff != 0.0); /*lint !e777*/

   value = solver->incumbentassignment[varidx];
   residual = solver->lhs[constridx] - value * coeff;
   movevalue = (cons->rhs - residual) / coeff;

   /* round, verify, and adjust */
   if( var->vartype == LS_CONTINUOUS )
   {
      if( SCIPisFeasPositive(scip, movevalue * coeff + residual - cons->rhs) )
      {
         if( coeff < 0.0 )
            movevalue = nextafter(movevalue, SCIP_DEFAULT_INFINITY);
         else
            movevalue = nextafter(movevalue, -SCIP_DEFAULT_INFINITY);
      }
   }
   else
   {
      movevalue = round(movevalue);

      if( SCIPisFeasPositive(scip, movevalue * coeff + residual - cons->rhs) )
      {
         if( coeff < 0.0 )
            movevalue += 1.0;
         else
            movevalue -= 1.0;
      }
   }

   /* clamp move value to variable bounds */
   if( movevalue < var->lb )
      movevalue = var->lb;
   else if( movevalue > var->ub )
      movevalue = var->ub;

   *result = movevalue;

   /* check whether clamped move value satisfies the constraint */
   return !SCIPisFeasPositive(scip, movevalue * coeff + residual - cons->rhs);
}

/** computes move score */
static
SCIP_Longint lsSolverTightScore(
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
   SCIP_Real newobj;
   SCIP_Bool prebetter;
   SCIP_Bool nowbetter;
   int constridx;
   SCIP_Real coeff;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real newlhs;
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

   /* objective scoring */
   if( !SCIPisInfinity(scip, solver->objcutoff) && var->objidx >= 0 )
   {
      objresidual = solver->incumbentobjective - value * var->obj;
      newobj = newvalue * var->obj + objresidual;

      if( ( !SCIPisInfinity(scip, newobj) || !SCIPisInfinity(scip, solver->incumbentobjective) )
         && ( !SCIPisInfinity(scip, -newobj) || !SCIPisInfinity(scip, -solver->incumbentobjective) )
         && !SCIPisEQ(scip, newobj, solver->incumbentobjective) )
      {
         if( newobj < solver->incumbentobjective )
            score += solver->objweight;
         else
            score -= solver->objweight;
      }

      /* subscore: crossing objcutoff boundary */
      prebetter = SCIPisLT(scip, solver->incumbentobjective, solver->objcutoff);
      nowbetter = SCIPisLT(scip, newobj, solver->objcutoff);

      if( !prebetter && nowbetter )
         solver->subscore += solver->objweight;
      else if( prebetter && !nowbetter )
         solver->subscore -= solver->objweight;
   }

   /* regular constraint scoring */
   for( i = 0; i < var->ncoeffs; ++i )
   {
      constridx = var->coeffs[i].idx;
      coeff = var->coeffs[i].coeff;
      lhs = solver->lhs[constridx];
      rhs = problem->conss[constridx].rhs;
      residual = lhs - value * coeff;
      newlhs = newvalue * coeff + residual;
      previol = lhs - rhs;
      newviol = newlhs - rhs;
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

      /* subscore: stable transitions */
      prestable = SCIPisFeasNegative(scip, previol);
      nowstable = SCIPisFeasNegative(scip, newviol);

      if( !prestable && nowstable )
         solver->subscore += solver->weight[constridx];
      else if( prestable && !nowstable )
         solver->subscore -= solver->weight[constridx];
   }

   solver->totaleffort += var->ncoeffs + 1;

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

   for( i = 0; i < solver->nunsat; ++i )
   {
      constridx = solver->unsatidxs[i];
      ++solver->weight[constridx];
   }

   if( solver->nunsat == 0 )
      ++solver->objweight;

   solver->totaleffort += solver->nunsat + 1;
}

/** smooths weights of satisfied constraints */
static
void lsSolverSmoothWeight(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver              /**< solver */
   )
{
   LS_PROBLEM* problem;
   SCIP_Real violation;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);

   problem = solver->problem;

   for( i = 0; i < problem->nconss; ++i )
   {
      violation = solver->lhs[i] - problem->conss[i].rhs;

      if( !SCIPisFeasPositive(scip, violation) && solver->weight[i] > 0 )
         --solver->weight[i];
   }

   if( solver->objweight > 0 && !SCIPisInfinity(scip, solver->incumbentobjective)
      && SCIPisLT(scip, solver->incumbentobjective, solver->objcutoff) )
      --solver->objweight;

   solver->totaleffort += problem->nconss + 1;
}


/*
 * Move operators
 */

/** performs the unsat tight move */
static
SCIP_RETCODE lsSolverUnsatTightMove(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_Bool*            result              /**< pointer to store whether move was applied */
   );

/** performs the sat tight move */
static
SCIP_RETCODE lsSolverSatTightMove(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_Bool*            result              /**< pointer to store whether move was applied */
   );

/** performs the flip move */
static
SCIP_RETCODE lsSolverFlipMove(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_Bool*            result              /**< pointer to store whether move was applied */
   );

/** performs the random tight move */
static
SCIP_RETCODE lsSolverRandomTightMove(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   );

/** performs the lift move */
static
SCIP_RETCODE lsSolverLiftMove(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_Bool*            result              /**< pointer to store whether move was applied */
   );

/** collects neighbor candidates from a constraint */
static
void collectConstraintNeighbors(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   int                   constridx,          /**< constraint index */
   SCIP_Bool             towardsatisfy       /**< TRUE: fallback toward satisfaction,
                                              *   FALSE: fallback toward boundary */
   )
{
   LS_PROBLEM* problem;
   LS_VAR* var;
   LS_CONSTRAINT* cons;
   SCIP_Real value;
   SCIP_Real movevalue;
   int varidx;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);

   problem = solver->problem;
   cons = problem->conss + constridx;

   for( i = 0; i < cons->ncoeffs; ++i )
   {
      varidx = cons->coeffs[i].idx;
      var = problem->vars + varidx;
      value = solver->incumbentassignment[varidx];

      if( !lsSolverTightMove(scip, solver, constridx, i, &movevalue) && !towardsatisfy )
      {
         /* move to opposite bound that tightens constraint */
         if( cons->coeffs[i].coeff < 0.0 )
            movevalue = var->lb;
         else
            movevalue = var->ub;
      }

      /* tabu check */
      if( (movevalue < value && solver->curstep < solver->allowdecstep[varidx])
         || (movevalue > value && solver->curstep < solver->allowincstep[varidx]) )
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
   LS_VAR* var;
   SCIP_Real value;
   SCIP_Real residual;
   SCIP_Real movevalue;
   SCIP_Real objtarget;
   int varidx;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);

   problem = solver->problem;
   objtarget = solver->objcutoff - SCIPcutoffbounddelta(scip);

   for( i = 0; i < problem->nobjvars; ++i )
   {
      varidx = problem->objvaridxs[i];
      var = problem->vars + varidx;

      assert(var->obj != 0.0); /*lint !e777*/

      value = solver->incumbentassignment[varidx];
      residual = solver->incumbentobjective - value * var->obj;
      movevalue = (objtarget - residual) / var->obj;

      /* round, verify, and adjust */
      if( var->vartype == LS_CONTINUOUS )
      {
         if( SCIPisFeasPositive(scip, movevalue * var->obj + residual - objtarget) )
         {
            if( var->obj < 0.0 )
               movevalue = nextafter(movevalue, SCIP_DEFAULT_INFINITY);
            else
               movevalue = nextafter(movevalue, -SCIP_DEFAULT_INFINITY);
         }
      }
      else
      {
         movevalue = round(movevalue);

         if( SCIPisFeasPositive(scip, movevalue * var->obj + residual - objtarget) )
         {
            if( var->obj < 0.0 )
               movevalue += 1.0;
            else
               movevalue -= 1.0;
         }
      }

      /* clamp move value to variable bounds */
      if( movevalue < var->lb )
         movevalue = var->lb;
      else if( movevalue > var->ub )
         movevalue = var->ub;

      /* tabu check */
      if( softtabu )
      {
         if( (movevalue < value && solver->curstep == solver->lastincstep[varidx] + 1)
            || (movevalue > value && solver->curstep == solver->lastdecstep[varidx] + 1) )
            continue;
      }
      else
      {
         if( (movevalue < value && solver->curstep < solver->allowdecstep[varidx])
            || (movevalue > value && solver->curstep < solver->allowincstep[varidx]) )
            continue;
      }

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
   LS_VAR* var;
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
      var = problem->vars + varidx;

      /* dedup binaries */
      if( var->vartype == LS_BINARY )
      {
         if( solver->scoretable[varidx] )
            continue;
         solver->scoretable[varidx] = TRUE;
         solver->scoreidxs[solver->nscoreidxs++] = varidx;
      }

      score = lsSolverTightScore(scip, solver, varidx, movevalue);

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
   SCIP_Bool*            result              /**< pointer to store whether move was applied */
   )
{
   int nsample;
   int budget;
   int bestvaridx;
   SCIP_Real bestvalue;
   int randidx;
   int tmp;
   SCIP_Bool satresult;
   SCIP_Bool flipresult;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);

   *result = FALSE;
   solver->neighborsize = 0;
   solver->nscoreidxs = 0;

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
         {
            collectConstraintNeighbors(scip, solver, solver->tempunsatidxs[i], TRUE);
         }
      }
      else
      {
         for( i = 0; i < solver->nunsat; ++i )
         {
            collectConstraintNeighbors(scip, solver, solver->unsatidxs[i], TRUE);
         }
      }
   }

   /* collect objective terms if needed */
   if( SCIPisInfinity(scip, solver->incumbentobjective)
      || SCIPisGE(scip, solver->incumbentobjective, solver->objcutoff) )
      collectObjectiveNeighbors(scip, solver, FALSE);

   /* subsample to budget */
   budget = !SCIPisInfinity(scip, solver->objcutoff) ? heurdata->bmsunsatfeas : heurdata->bmsunsatinfeas;
   subsampleNeighbors(solver, heurdata, budget);

   /* score and select */
   if( selectBestNeighbor(scip, solver, (SCIP_Longint)0, &bestvaridx, &bestvalue) )
   {
      lsSolverApplyMove(scip, solver, heurdata, bestvaridx, bestvalue);
      cleanScoreTable(solver);
      *result = TRUE;
      return SCIP_OKAY;
   }

   /* fallback: try SatTightMove (if feasible), then FlipMove */
   if( !SCIPisInfinity(scip, solver->objcutoff) )
   {
      satresult = FALSE;

      SCIP_CALL( lsSolverSatTightMove(scip, solver, heurdata, &satresult) );
      if( satresult )
      {
         cleanScoreTable(solver);
         *result = TRUE;
         return SCIP_OKAY;
      }
   }

   flipresult = FALSE;

   SCIP_CALL( lsSolverFlipMove(scip, solver, heurdata, &flipresult) );
   cleanScoreTable(solver);
   *result = flipresult;

   return SCIP_OKAY;
}

/** performs the sat tight move */
static
SCIP_RETCODE lsSolverSatTightMove(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_Bool*            result              /**< pointer to store whether move was applied */
   )
{
   LS_PROBLEM* problem;
   int nsaved;
   int bestvaridx;
   SCIP_Real bestvalue;
   int constridx;
   SCIP_Real violation;
   int newcount;
   int randidx;
   int tmpvar;
   SCIP_Real tmpvalue;
   SCIP_Longint bestscore;
   SCIP_Longint bestsubscore;
   int varidx;
   SCIP_Real movevalue;
   LS_VAR* var;
   SCIP_Longint score;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);

   *result = FALSE;
   problem = solver->problem;
   nsaved = solver->neighborsize;

   if( problem->nconss == 0 )
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
      violation = solver->lhs[constridx] - problem->conss[constridx].rhs;
      if( SCIPisFeasPositive(scip, violation) )
         continue;

      /* mark as sampled */
      solver->sampledconstrs[constridx] = TRUE;
      solver->sampledidxs[solver->nsampled++] = constridx;

      collectConstraintNeighbors(scip, solver, constridx, FALSE);
   }

   /* clean up sampled flags */
   for( i = 0; i < solver->nsampled; ++i )
      solver->sampledconstrs[solver->sampledidxs[i]] = FALSE;
   solver->nsampled = 0;

   /* subsample to budget (only new entries) */
   if( solver->neighborsize - nsaved > heurdata->bmssat )
   {
      /* shuffle only the new entries */
      newcount = solver->neighborsize - nsaved;

      for( i = 0; i < heurdata->bmssat; ++i )
      {
         randidx = SCIPrandomGetInt(heurdata->randnumgen, 0, newcount - 1 - i) + i;

         tmpvar = solver->neighborvaridxs[nsaved + randidx];
         tmpvalue = solver->neighborvalues[nsaved + randidx];
         solver->neighborvaridxs[nsaved + randidx] = solver->neighborvaridxs[nsaved + i];
         solver->neighborvalues[nsaved + randidx] = solver->neighborvalues[nsaved + i];
         solver->neighborvaridxs[nsaved + i] = tmpvar;
         solver->neighborvalues[nsaved + i] = tmpvalue;
      }
      solver->neighborsize = nsaved + heurdata->bmssat;
   }

   /* score from nsaved onward using shared scoretable */
   bestscore = 0;
   bestsubscore = -SCIP_LONGINT_MAX;

   bestvaridx = -1;
   bestvalue = 0.0;

   for( i = nsaved; i < solver->neighborsize; ++i )
   {
      varidx = solver->neighborvaridxs[i];
      movevalue = solver->neighborvalues[i];
      var = problem->vars + varidx;

      if( var->vartype == LS_BINARY )
      {
         if( solver->scoretable[varidx] )
            continue;
         solver->scoretable[varidx] = TRUE;
         solver->scoreidxs[solver->nscoreidxs++] = varidx;
      }

      score = lsSolverTightScore(scip, solver, varidx, movevalue);

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

/** performs the flip move */
static
SCIP_RETCODE lsSolverFlipMove(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_Bool*            result              /**< pointer to store whether move was applied */
   )
{
   LS_PROBLEM* problem;
   SCIP_Longint bestscore;
   SCIP_Longint bestsubscore;
   int bestvaridx;
   SCIP_Real bestvalue;
   int binidx;
   int varidx;
   LS_VAR* var;
   SCIP_Real value;
   SCIP_Real movevalue;
   SCIP_Longint score;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);

   *result = FALSE;
   problem = solver->problem;

   if( problem->nbinary == 0 )
      return SCIP_OKAY;

   bestscore = 0;
   bestsubscore = -SCIP_LONGINT_MAX;
   bestvaridx = -1;
   bestvalue = 0.0;

   for( i = 0; i < heurdata->bmsflip; ++i )
   {
      binidx = SCIPrandomGetInt(heurdata->randnumgen, 0, problem->nbinary - 1);
      varidx = problem->binaryidxs[binidx];
      var = problem->vars + varidx;

      /* dedup */
      if( solver->scoretable[varidx] )
         continue;
      solver->scoretable[varidx] = TRUE;
      solver->scoreidxs[solver->nscoreidxs++] = varidx;

      value = solver->incumbentassignment[varidx];
      movevalue = (value > 0.5) ? 0.0 : 1.0;

      /* bounds and tabu check */
      if( ( movevalue < value
         && ( SCIPisFeasPositive(scip, var->lb) || solver->curstep < solver->allowdecstep[varidx] ) )
         || ( movevalue > value
         && ( SCIPisFeasNegative(scip, var->ub - 1.0) || solver->curstep < solver->allowincstep[varidx] ) ) )
         continue;

      score = lsSolverTightScore(scip, solver, varidx, movevalue);

      if( score > bestscore
         || (score == bestscore && solver->subscore > bestsubscore) )
      {
         bestscore = score;
         bestsubscore = solver->subscore;
         bestvaridx = varidx;
         bestvalue = movevalue;
      }
   }

   solver->totaleffort += heurdata->bmsflip;

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
   LS_PROBLEM* problem;
   LS_CONSTRAINT* cons;
   int bestvaridx;
   SCIP_Real bestvalue;
   int constridx;
   int varidx;
   SCIP_Real value;
   SCIP_Real movevalue;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(heurdata != NULL);

   problem = solver->problem;
   solver->neighborsize = 0;

   /* pick 1 random unsat constraint */
   if( solver->nunsat > 0 )
   {
      constridx = solver->unsatidxs[SCIPrandomGetInt(heurdata->randnumgen, 0, solver->nunsat - 1)];
      cons = problem->conss + constridx;

      /* collect terms with soft tabu */
      for( i = 0; i < cons->ncoeffs; ++i )
      {
         varidx = cons->coeffs[i].idx;
         value = solver->incumbentassignment[varidx];

         (void)lsSolverTightMove(scip, solver, constridx, i, &movevalue);

         /* soft aspiration: only block immediate reversal */
         if( (movevalue < value && solver->curstep == solver->lastincstep[varidx] + 1)
            || (movevalue > value && solver->curstep == solver->lastdecstep[varidx] + 1) )
            continue;

         if( SCIPisEQ(scip, movevalue, value) )
            continue;

         assert(solver->neighborsize < solver->neighborcap);
         solver->neighborvaridxs[solver->neighborsize] = varidx;
         solver->neighborvalues[solver->neighborsize] = movevalue;
         ++solver->neighborsize;
      }
   }

   /* collect objective terms if needed */
   if( SCIPisInfinity(scip, solver->incumbentobjective)
      || SCIPisGE(scip, solver->incumbentobjective, solver->objcutoff) )
      collectObjectiveNeighbors(scip, solver, TRUE);

   /* subsample */
   subsampleNeighbors(solver, heurdata, heurdata->bmsrandom);

   /* score ALL (accept negative scores) */
   if( selectBestNeighbor(scip, solver, -SCIP_LONGINT_MAX, &bestvaridx, &bestvalue) )
      lsSolverApplyMove(scip, solver, heurdata, bestvaridx, bestvalue);
   cleanScoreTable(solver);

   return SCIP_OKAY;
}

/** performs the lift move (feasible objective improvement) */
static
SCIP_RETCODE lsSolverLiftMove(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_Bool*            result              /**< pointer to store whether move was applied */
   )
{
   LS_PROBLEM* problem;
   SCIP_Real bestobjvalue;
   int bestvaridx;
   SCIP_Real bestvalue;
   int bestlaststep;
   LS_VAR* var;
   LS_VAR* bestvar;
   LS_VAR* affvar;
   LS_CONSTRAINT* cons;
   SCIP_Real value;
   SCIP_Real objresidual;
   SCIP_Real objcoeff;
   SCIP_Real objvalue;
   SCIP_Real varvalue;
   SCIP_Real lv;
   SCIP_Real uv;
   SCIP_Real gap;
   SCIP_Real gap2;
   SCIP_Real coeff;
   SCIP_Real coeff2;
   SCIP_Real movevalue;
   SCIP_Real movevalue2;
   int varidx;
   int constridx;
   int cidx;
   int affidx;
   int objpos;
   int laststep;
   int i;
   int j;
   int k;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);

   *result = FALSE;
   problem = solver->problem;

   if( SCIPisInfinity(scip, solver->objcutoff) || problem->nobjvars == 0 )
      return SCIP_OKAY;

   /* compute feasible value ranges */
   if( !solver->iskeepfeas )
   {
      for( i = 0; i < problem->nobjvars; ++i )
      {
         varidx = problem->objvaridxs[i];
         var = problem->vars + varidx;

         solver->liftlowervalue[i] = var->lb;
         solver->liftuppervalue[i] = var->ub;

         for( j = 0; j < var->ncoeffs; ++j )
         {
            constridx = var->coeffs[j].idx;
            gap = solver->lhs[constridx] - problem->conss[constridx].rhs;
            coeff = var->coeffs[j].coeff;

            if( SCIPisFeasZero(scip, gap) )
            {
               /* constraint is tight */
               if( coeff > 0.0 )
                  solver->liftuppervalue[i] = solver->incumbentassignment[varidx];
               else
                  solver->liftlowervalue[i] = solver->incumbentassignment[varidx];
            }
            else if( lsSolverTightMove(scip, solver, constridx, var->coeffs[j].pos, &movevalue) )
            {
               if( coeff > 0.0 )
               {
                  if( movevalue < solver->liftuppervalue[i] )
                     solver->liftuppervalue[i] = movevalue;
               }
               else
               {
                  if( movevalue > solver->liftlowervalue[i] )
                     solver->liftlowervalue[i] = movevalue;
               }
            }

            if( solver->liftlowervalue[i] >= solver->liftuppervalue[i] )
               break;
         }

         solver->totaleffort += var->ncoeffs;
      }
   }

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
      lv = solver->liftlowervalue[i];
      uv = solver->liftuppervalue[i];

      if( lv == uv ) /*lint !e777*/
         continue;

      /* pick direction improving objective */
      objresidual = solver->incumbentobjective - value * objcoeff;
      if( objcoeff > 0.0 )
      {
         objvalue = lv * objcoeff + objresidual;
         varvalue = lv;
      }
      else
      {
         objvalue = uv * objcoeff + objresidual;
         varvalue = uv;
      }

      /* skip zero move */
      if( SCIPisEQ(scip, varvalue, value) )
         continue;

      laststep = (varvalue < value) ? solver->lastdecstep[varidx] : solver->lastincstep[varidx];

      if( ( SCIPisInfinity(scip, -objvalue) || SCIPisInfinity(scip, bestobjvalue)
         || SCIPisLE(scip, objvalue, bestobjvalue) )
         && ( laststep < bestlaststep || ( !SCIPisInfinity(scip, objvalue) && !SCIPisInfinity(scip, -bestobjvalue)
         && SCIPisLT(scip, objvalue, bestobjvalue) ) ) )
      {
         bestobjvalue = objvalue;
         bestvaridx = varidx;
         bestvalue = varvalue;
         bestlaststep = laststep;
      }
   }

   if( bestvaridx >= 0 )
   {
      lsSolverApplyMove(scip, solver, heurdata, bestvaridx, bestvalue);
      solver->iskeepfeas = TRUE;
      *result = TRUE;

      /* recompute value ranges for affected variables */
      bestvar = problem->vars + bestvaridx;
      solver->naffected = 0;

      for( j = 0; j < bestvar->ncoeffs; ++j )
      {
         constridx = bestvar->coeffs[j].idx;
         cons = problem->conss + constridx;

         for( k = 0; k < cons->ncoeffs; ++k )
         {
            affidx = cons->coeffs[k].idx;

            if( !solver->affectedset[affidx] )
            {
               solver->affectedset[affidx] = TRUE;
               solver->affectedvars[solver->naffected++] = affidx;
            }
         }
      }

      /* recompute value ranges for affected objective variables */
      for( j = 0; j < solver->naffected; ++j )
      {
         affidx = solver->affectedvars[j];
         affvar = problem->vars + affidx;
         objpos = affvar->objidx;

         solver->affectedset[affidx] = FALSE;

         if( objpos == -1 )
            continue;

         solver->liftlowervalue[objpos] = affvar->lb;
         solver->liftuppervalue[objpos] = affvar->ub;

         for( k = 0; k < affvar->ncoeffs; ++k )
         {
            cidx = affvar->coeffs[k].idx;
            gap2 = solver->lhs[cidx] - problem->conss[cidx].rhs;
            coeff2 = affvar->coeffs[k].coeff;

            if( SCIPisFeasZero(scip, gap2) )
            {
               if( coeff2 > 0.0 )
                  solver->liftuppervalue[objpos] = solver->incumbentassignment[affidx];
               else
                  solver->liftlowervalue[objpos] = solver->incumbentassignment[affidx];
            }
            else if( lsSolverTightMove(scip, solver, cidx, affvar->coeffs[k].pos, &movevalue2) )
            {
               if( coeff2 > 0.0 )
               {
                  if( movevalue2 < solver->liftuppervalue[objpos] )
                     solver->liftuppervalue[objpos] = movevalue2;
               }
               else
               {
                  if( movevalue2 > solver->liftlowervalue[objpos] )
                     solver->liftlowervalue[objpos] = movevalue2;
               }
            }

            if( solver->liftlowervalue[objpos] >= solver->liftuppervalue[objpos] )
               break;
         }

         solver->totaleffort += affvar->ncoeffs;
      }
   }
   else
   {
      solver->iskeepfeas = FALSE;
   }

   return SCIP_OKAY;
}


/*
 * SCIP integration
 */

/** adds a row to the local search problem */
static
SCIP_RETCODE addRowInLocalSolver(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_PROBLEM*           problem,            /**< problem */
   SCIP_Real*            consvals,           /**< coefficient values */
   int*                  consinds,           /**< variable indices */
   int                   consnvars,          /**< number of variables */
   SCIP_Real             constant,           /**< row activity constant */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs                 /**< right hand side */
   )
{
   SCIP_Real* negvals;
   int i;

   assert(scip != NULL);
   assert(problem != NULL);

   /* ranged row */
   if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) && !SCIPisEQ(scip, lhs, rhs) )
   {
      SCIP_CALL( lsProblemAddConstraint(scip, problem, constant, rhs, consnvars, consinds, consvals) );

      /* negated constraint for GTE side */
      SCIP_CALL( SCIPallocBufferArray(scip, &negvals, consnvars) );
      for( i = 0; i < consnvars; ++i )
         negvals[i] = -consvals[i];
      SCIP_CALL( lsProblemAddConstraint(scip, problem, -constant, -lhs, consnvars, consinds, negvals) );
      SCIPfreeBufferArray(scip, &negvals);
   }
   else if( SCIPisEQ(scip, rhs, lhs) )
   {
      /* equality: two constraints */
      SCIP_CALL( lsProblemAddConstraint(scip, problem, constant, rhs, consnvars, consinds, consvals) );

      SCIP_CALL( SCIPallocBufferArray(scip, &negvals, consnvars) );
      for( i = 0; i < consnvars; ++i )
         negvals[i] = -consvals[i];
      SCIP_CALL( lsProblemAddConstraint(scip, problem, -constant, -rhs, consnvars, consinds, negvals) );
      SCIPfreeBufferArray(scip, &negvals);
   }
   else if( SCIPisInfinity(scip, rhs) )
   {
      /* GTE */
      SCIP_CALL( SCIPallocBufferArray(scip, &negvals, consnvars) );
      for( i = 0; i < consnvars; ++i )
         negvals[i] = -consvals[i];
      SCIP_CALL( lsProblemAddConstraint(scip, problem, -constant, -lhs, consnvars, consinds, negvals) );
      SCIPfreeBufferArray(scip, &negvals);
   }
   else
   {
      /* LTE */
      assert(SCIPisInfinity(scip, -lhs));
      SCIP_CALL( lsProblemAddConstraint(scip, problem, constant, rhs, consnvars, consinds, consvals) );
   }

   return SCIP_OKAY;
}

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
   SCIP_Bool*            success             /**< was the problem successfully extracted? */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR** consvars;
   SCIP_CONS** conss;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR* var;
   SCIP_CONS* cons;
   SCIP_Real* consvals;
   int* nvarcoeffs;
   int* consinds;
   SCIP_Real constant;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int nvars;
   int nconss;
   int nlinconss;
   int consnvars;
   int requiredsize;
   int nbinary;
   int pos;
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

   if( nlinconss < nconss )
      return SCIP_OKAY;

   vars = SCIPgetVars(scip);
   conss = SCIPgetConss(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &consinds, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nvarcoeffs, nvars) );

   /* count resolved constraints and coefficients */
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

      /* resolve to active variables */
      constant = 0.0;
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, consvars, consvals, &consnvars, nvars, &constant, &requiredsize) );
      assert(requiredsize <= nvars);

      if( consnvars == 0 )
         continue;

      if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) )
      {
         problem->nconss += 2;
         for( j = 0; j < consnvars; ++j )
            nvarcoeffs[SCIPvarGetProbindex(consvars[j])] += 2;
      }
      else
      {
         ++problem->nconss;
         for( j = 0; j < consnvars; ++j )
            ++nvarcoeffs[SCIPvarGetProbindex(consvars[j])];
      }
   }

   /* allocate variable and constraint arrays */
   if( nvars > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->vars, nvars) );
   }
   if( problem->nconss > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->conss, problem->nconss) );

      /* zero-initialize for safe cleanup on early termination */
      BMSclearMemoryArray(problem->conss, problem->nconss);
   }

   /* add variables to problem */
   for( i = 0; i < nvars; ++i )
   {
      var = vars[i];
      SCIP_CALL( lsProblemAddVar(scip, problem,
            SCIPvarIsIntegral(var) ? (SCIPvarIsBinary(var) ? LS_BINARY : LS_INTEGER) : LS_CONTINUOUS,
            SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPvarGetObj(var)) );
   }

   /* allocate coefficient arrays */
   for( i = 0; i < nvars; ++i )
   {
      if( nvarcoeffs[i] > 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->vars[i].coeffs, nvarcoeffs[i]) );
      }
   }

   /* setup objective */
   SCIP_CALL( setupObjective(scip, problem) );

   /* reset nconss before adding constraints */
   problem->nconss = 0;

   /* add constraints to problem */
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

      /* resolve to active variables */
      constant = 0.0;
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, consvars, consvals, &consnvars, nvars, &constant, &requiredsize) );
      assert(requiredsize <= nvars);

      if( consnvars == 0 )
         continue;

      for( j = 0; j < consnvars; ++j )
      {
         consinds[j] = SCIPvarGetProbindex(consvars[j]);
         assert(consinds[j] >= 0);
      }

      SCIP_CALL( addRowInLocalSolver(scip, problem, consvals, consinds, consnvars, constant, lhs, rhs) );
   }

#ifndef NDEBUG
   for( i = 0; i < nvars; ++i )
      assert(problem->vars[i].ncoeffs == nvarcoeffs[i]);
#endif

   /* build binary index */
   nbinary = 0;

   for( i = 0; i < problem->nvars; ++i )
   {
      if( problem->vars[i].vartype == LS_BINARY )
         ++nbinary;
   }
   problem->nbinary = nbinary;

   if( nbinary > 0 )
   {
      pos = 0;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->binaryidxs, nbinary) );
      for( i = 0; i < problem->nvars; ++i )
      {
         if( problem->vars[i].vartype == LS_BINARY )
            problem->binaryidxs[pos++] = i;
      }
   }

   *success = TRUE;

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
   int                   nrows               /**< number of LP rows */
   )
{
   SCIP_COL** rowcols;
   SCIP_VAR* var;
   SCIP_COL* col;
   SCIP_ROW* row;
   SCIP_Real* vals;
   SCIP_Real* rowvals;
   int* nvarcoeffs;
   int* inds;
   SCIP_Real constant;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int nnonz;
   int nvals;
   int nbinary;
   int pos;
   int i;
   int j;

   assert(scip != NULL);
   assert(problem != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &inds, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, ncols) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nvarcoeffs, ncols) );

   /* count constraints and coefficients */
   for( i = 0; i < nrows; ++i )
   {
      row = rows[i];
      rowcols = SCIProwGetCols(row);
      nnonz = SCIProwGetNNonz(row);

      /* count columns for this row */
      nvals = 0;
      for( j = 0; j < nnonz; ++j )
      {
         if( SCIPcolIsInLP(rowcols[j]) )
            ++nvals;
      }

      /* skip rows without columns */
      if( nvals == 0 )
         continue;

      lhs = SCIProwGetLhs(row);
      rhs = SCIProwGetRhs(row);

      if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) )
      {
         problem->nconss += 2;
         for( j = 0; j < nnonz; ++j )
         {
            if( SCIPcolIsInLP(rowcols[j]) )
               nvarcoeffs[SCIPcolGetLPPos(rowcols[j])] += 2;
         }
      }
      else
      {
         ++problem->nconss;
         for( j = 0; j < nnonz; ++j )
         {
            if( SCIPcolIsInLP(rowcols[j]) )
               ++nvarcoeffs[SCIPcolGetLPPos(rowcols[j])];
         }
      }
   }

   /* allocate variable and constraint arrays */
   if( ncols > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->vars, ncols) );
   }
   if( problem->nconss > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->conss, problem->nconss) );
   }

   /* add variables to problem */
   for( i = 0; i < ncols; ++i )
   {
      col = cols[i];
      var = SCIPcolGetVar(col);
      SCIP_CALL( lsProblemAddVar(scip, problem,
            SCIPvarIsIntegral(var) ? (SCIPvarIsBinary(var) ? LS_BINARY : LS_INTEGER) : LS_CONTINUOUS,
            SCIPcolGetLb(col), SCIPcolGetUb(col), SCIPcolGetObj(col)) );
   }

   /* allocate coefficient arrays */
   for( i = 0; i < ncols; ++i )
   {
      if( nvarcoeffs[i] > 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->vars[i].coeffs, nvarcoeffs[i]) );
      }
   }

   /* setup objective */
   SCIP_CALL( setupObjective(scip, problem) );

   /* reset nconss before adding constraints */
   problem->nconss = 0;

   /* add constraints to problem */
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
            ++nvals;
         }
      }

      if( nvals == 0 )
         continue;

      SCIP_CALL( addRowInLocalSolver(scip, problem, vals, inds, nvals, constant, lhs, rhs) );
   }

#ifndef NDEBUG
   for( i = 0; i < ncols; ++i )
      assert(problem->vars[i].ncoeffs == nvarcoeffs[i]);
#endif

   SCIPfreeBufferArray(scip, &nvarcoeffs);
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &inds);

   /* build binary index */
   nbinary = 0;

   for( i = 0; i < problem->nvars; ++i )
   {
      if( problem->vars[i].vartype == LS_BINARY )
         ++nbinary;
   }
   problem->nbinary = nbinary;

   if( nbinary > 0 )
   {
      pos = 0;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->binaryidxs, nbinary) );
      for( i = 0; i < problem->nvars; ++i )
      {
         if( problem->vars[i].vartype == LS_BINARY )
            problem->binaryidxs[pos++] = i;
      }
   }

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
   SCIP_SOL* sol;
   SCIP_Bool stored;
   SCIP_VAR** vars;
   int nvars;
   SCIP_COL** cols;
   int ncols;
   int i;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(solver != NULL);
   assert(result != NULL);

   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

   if( heurtiming == SCIP_HEURTIMING_BEFOREPRESOL )
   {
      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);

      for( i = 0; i < nvars; ++i )
      {
         SCIP_CALL( SCIPsetSolVal(scip, sol, vars[i], solver->incumbentassignment[i]) );
      }
   }
   else
   {
      SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

      for( i = 0; i < ncols; ++i )
      {
         SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPcolGetVar(cols[i]), solver->incumbentassignment[i]) );
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
            solver->objcutoff = SCIPgetCutoffbound(scip);
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
      if( !success )
      {
         SCIP_CALL( lsProblemFree(scip, &problem) );
         return SCIP_OKAY;
      }
   }
   else
   {
      SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
      SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

      SCIP_CALL( extractProblemData(scip, problem, cols, rows, ncols, nrows) );
   }

   if( problem->nvars == 0 || problem->nconss == 0 )
   {
      SCIP_CALL( lsProblemFree(scip, &problem) );
      return SCIP_OKAY;
   }

   /* create solver */
   SCIP_CALL( lsSolverCreate(scip, &solver, problem) );

   /* initialize solution */
   bestsol = (heurdata->useinitialsol && heurtiming != SCIP_HEURTIMING_BEFOREPRESOL) ? SCIPgetBestSol(scip) : NULL;
   SCIP_CALL( lsSolverInitSolution(scip, solver, bestsol, cols, ncols) );
   SCIP_CALL( lsSolverInitState(scip, solver) );

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
            && SCIPisLT(scip, solver->incumbentobjective, solver->objcutoff) )
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
         if( SCIPrandomGetInt(heurdata->randnumgen, 0, 9999) > heurdata->smoothprob )
            lsSolverUpdateWeight(scip, solver);
         else
            lsSolverSmoothWeight(scip, solver);

         SCIP_CALL( lsSolverRandomTightMove(scip, solver, heurdata) );
      }

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
SCIP_DECL_HEURCOPY(heurCopyLocal)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   SCIP_CALL( SCIPincludeHeurLocal(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeLocal)
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
SCIP_DECL_HEURINIT(heurInitLocal)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen, DEFAULT_RANDSEED, TRUE) );

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitLocal)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeRandom(scip, &heurdata->randnumgen);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecLocal)
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

/** creates the local search primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurLocal(
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
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecLocal, heurdata) );

   assert(heur != NULL);

   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyLocal) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeLocal) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitLocal) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitLocal) );

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

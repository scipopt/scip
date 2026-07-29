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
#define DEFAULT_BMSUNSATINFEAS     2000      /**< max moves scored during infeasible phase */
#define DEFAULT_BMSUNSATFEAS       3000      /**< max moves scored during feasible phase */
#define DEFAULT_SAMPLESAT          20        /**< sat constraints to sample */
#define DEFAULT_BMSSAT             190       /**< max moves scored from sat constraints */
#define DEFAULT_BMSFLIP            20        /**< binary flip candidates */
#define DEFAULT_BMSRANDOM          150       /**< random move candidates */
#define DEFAULT_VERBOSITY          0         /**< verbosity level of the local search solver */
#define DEFAULT_MAXSOLS            1         /**< maximum number of solutions to find (-1: unlimited) */
#define DEFAULT_MINDECREASEREL     0.1       /**< minimum relative decrease to reset the effort budget */
#define DEFAULT_MAXEFFORTFAC       256       /**< maximum effort factor per nonzero without improvement (-1: unlimited) */
#define DEFAULT_RECOMPUTEFREQ      0         /**< activity recomputation frequency (-1: off, 0: near boundary only) */
#define DEFAULT_USEINITIALSOL      TRUE      /**< should best known solution be used? */
#define DEFAULT_USETERMINALSOL     TRUE      /**< should previous terminal solution be used? */
#define DEFAULT_ONLYMAINSCIP       FALSE     /**< should run in subscip be skipped? */
#define DEFAULT_ONLYBEFORENODE     TRUE      /**< should run before presolving be skipped? */
#define DEFAULT_ONLYWITHOUTSOL     FALSE     /**< should run with solution be skipped? */


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   SCIP_SOL*             terminalsol;        /**< terminal solution to warm-start the next call, or NULL */
   int                   smoothprob;         /**< smooth probability (out of 10000) */
   int                   tabubase;           /**< minimum tabu tenure */
   int                   tabuvariation;      /**< random tabu variation */
   int                   sampleunsat;        /**< unsat constraints to sample */
   int                   bmsunsatinfeas;     /**< max moves scored during infeasible phase */
   int                   bmsunsatfeas;       /**< max moves scored during feasible phase */
   int                   samplesat;          /**< sat constraints to sample */
   int                   bmssat;             /**< max moves scored from sat constraints */
   int                   bmsflip;            /**< binary flip candidates */
   int                   bmsrandom;          /**< random move candidates */
   int                   verbosity;          /**< verbosity level of the local search solver */
   int                   maxsols;            /**< maximum number of solutions to find (-1: unlimited) */
   SCIP_Real             mindecreaserel;     /**< minimum relative decrease to reset the effort budget */
   int                   maxeffortfac;       /**< maximum effort factor per nonzero without improvement (-1: unlimited) */
   int                   recomputefreq;      /**< activity recomputation frequency (-1: off, 0: near boundary only) */
   SCIP_Bool             useinitialsol;      /**< should best known solution be used? */
   SCIP_Bool             useterminalsol;     /**< should previous terminal solution be used? */
   SCIP_Bool             onlymainscip;       /**< should run in subscip be skipped? */
   SCIP_Bool             onlybeforenode;     /**< should run before presolving be skipped? */
   SCIP_Bool             onlywithoutsol;     /**< should run with solution be skipped? */
   SCIP_Longint          maxeffort;          /**< absolute effort budget */
};

/**@todo split into simple data arrays */
/** index-coefficient pair for cross-references */
struct LS_IdxCoeff
{
   int                   idx;                /**< variable or constraint index */
   int                   pos;                /**< position in the other direction's array */
   SCIP_Real             coeff;              /**< normalized coefficient */
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

/** constraint structure (static model data) */
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
   SCIP_Real             locallowerbound;    /**< local lower bound */
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
   SCIP_Bool             objunsat;           /**< is the objective target unreached? */
   int                   objnmoves;          /**< objective moves since its last computation */
   SCIP_Real*            act;                /**< per-constraint current activity */
   int*                  unsatidxs;          /**< violated constraint indices */
   int                   nunsat;             /**< number of violated constraints */
   int*                  satidxs;            /**< satisfied constraint indices */
   int                   nsat;               /**< number of satisfied constraints */
   int*                  idxs;               /**< per-constraint index: >=0 unsatidxs position, <0 satidxs position */
   int*                  consnmoves;         /**< per-constraint moves since their last computation */
   int*                  weight;             /**< per-constraint penalty weight */
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
   SCIP_Bool*            sampledconstrs;     /**< per-constraint dedup flag for sat sampling */
   int*                  sampledidxs;        /**< stack of sampled constraint indices */
   int                   nsampled;           /**< number of sampled constraints */
   SCIP_Real*            liftbound;          /**< per obj-var: feasible bound in improving direction */
   int*                  liftidxs;           /**< per obj-var: lift coefficient position, -1 if free */
   SCIP_Longint          subscore;           /**< secondary score (set by tightScore) */
   int                   nsols;              /**< number of solutions */
   int                   previolations;      /**< violation count at last improvement */
   SCIP_Real             preobjective;       /**< objective value at last improvement */
   SCIP_Longint          effort;             /**< effort spent since last improvement */
};
typedef struct LS_Solver LS_SOLVER;


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

/** allocates and fills variable cross-references */
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


/*
 * Local methods for LS_Solver
 */

/** declares constraint violated: moves from satisfied to violated set */
static
void lsSolverMakeUnsat(
   LS_SOLVER*            solver,             /**< solver */
   int                   constridx           /**< constraint index */
   )
{
   int lastconstridx;
   int pos;

   assert(solver != NULL);
   assert(constridx >= 0);
   assert(constridx < solver->problem->nconss);
   assert(solver->idxs[constridx] < 0);
   assert(solver->nsat > 0);
   assert(solver->nunsat < solver->problem->nconss);

   lastconstridx = solver->satidxs[--solver->nsat];
   pos = -solver->idxs[constridx] - 1;
   solver->satidxs[pos] = lastconstridx;
   solver->idxs[lastconstridx] = -(pos + 1);
   solver->unsatidxs[solver->nunsat] = constridx;
   solver->idxs[constridx] = solver->nunsat++;
}

/** declares constraint satisfied: moves from violated to satisfied set */
static
void lsSolverMakeSat(
   LS_SOLVER*            solver,             /**< solver */
   int                   constridx           /**< constraint index */
   )
{
   int lastconstridx;
   int pos;

   assert(solver != NULL);
   assert(constridx >= 0);
   assert(constridx < solver->problem->nconss);
   assert(solver->idxs[constridx] >= 0);
   assert(solver->nunsat > 0);
   assert(solver->nsat < solver->problem->nconss);

   lastconstridx = solver->unsatidxs[--solver->nunsat];
   pos = solver->idxs[constridx];
   solver->unsatidxs[pos] = lastconstridx;
   solver->idxs[lastconstridx] = pos;
   solver->satidxs[solver->nsat] = constridx;
   solver->idxs[constridx] = -(++solver->nsat);
}

/** recomputes objective value from current incumbent assignment */
static
void lsSolverRecomputeObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver              /**< solver */
   )
{
   LS_PROBLEM* problem = solver->problem;
   int i;

   solver->incumbentobjective = 0.0;

   for( i = 0; i < problem->nobjvars; ++i )
   {
      /**@todo use quad precision if necessary */
      solver->incumbentobjective += problem->vars[problem->objvaridxs[i]].obj
            * solver->incumbentassignment[problem->objvaridxs[i]];
   }

   solver->objnmoves = 0;
}

/** updates target flag from current objective value */
static
void lsSolverUpdateObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver              /**< solver */
   )
{
   solver->objunsat = !SCIPisInfinity(scip, solver->objtarget)
         && SCIPisGT(scip, solver->incumbentobjective, solver->objtarget);
}

/** recomputes constraint activity from current incumbent assignment */
static
void lsSolverRecomputeConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   int                   constridx           /**< constraint index */
   )
{
   LS_CONSTRAINT* cons = solver->problem->conss + constridx;
   int i;

   solver->act[constridx] = 0.0;

   for( i = 0; i < cons->ncoeffs; ++i )
   {
      /**@todo use quad precision if necessary */
      solver->act[constridx] += cons->coeffs[i].coeff * solver->incumbentassignment[cons->coeffs[i].idx];
   }

   solver->consnmoves[constridx] = 0;
}

/** updates violation indicator from current constraint activity */
static
void lsSolverUpdateConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   int                   constridx           /**< constraint index */
   )
{
   LS_CONSTRAINT* cons = solver->problem->conss + constridx;

   if( SCIPisFeasPositive(scip, 2.0 * (!SCIPisInfinity(scip, -cons->mhs) && solver->act[constridx] <= cons->mhs
      ? cons->lhs - solver->act[constridx] : solver->act[constridx] - cons->rhs)) )
   {
      if( solver->idxs[constridx] < 0 )
         lsSolverMakeUnsat(solver, constridx);
   }
   else
   {
      if( solver->idxs[constridx] >= 0 )
         lsSolverMakeSat(solver, constridx);
   }
}

/** creates and initializes a solver */
static
SCIP_RETCODE lsSolverCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURTIMING       heurtiming,         /**< heuristic timing mask */
   SCIP_SOL*             sol,                /**< SCIP solution, or NULL for bounds */
   LS_PROBLEM*           problem,            /**< problem */
   LS_SOLVER**           solverptr           /**< pointer to solver */
   )
{
   LS_SOLVER* solver;
   SCIP_VAR** vars;
   SCIP_COL** cols;
   SCIP_Real upperbound;
   int nvars;
   int nobjvars;
   int nnonzeros;
   int nconss;
   int i;

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
   solver->iskeepfeas = FALSE;
   solver->objtarget = SCIPisObjIntegral(scip)
         ? floor(SCIPgetCutoffbound(scip)) : SCIPgetCutoffbound(scip) - SCIPcutoffbounddelta(scip);
   solver->objweight = 1;
   solver->subscore = 0;
   solver->neighborsize = 0;
   solver->nscoreidxs = 0;
   solver->nsampled = 0;
   solver->nsols = 0;
   solver->nunsat = 0;
   solver->nsat = 0;

   /* allocate incumbent assignment */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->incumbentassignment, nvars) );

   /* allocate per-constraint parallel arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->act, nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->unsatidxs, nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->satidxs, nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->idxs, nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->consnmoves, nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->weight, nconss) );

   /* allocate per-variable tabu arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->allowincstep, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->allowdecstep, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->lastincstep, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->lastdecstep, nvars) );

   /* allocate work arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->neighborvaridxs, nnonzeros) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->neighborvalues, nnonzeros) );
   solver->neighborcap = nnonzeros;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->scoretable, nvars) );
   BMSclearMemoryArray(solver->scoretable, nvars);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->scoreidxs, problem->nbinvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->sampledconstrs, nconss) );
   BMSclearMemoryArray(solver->sampledconstrs, nconss);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->sampledidxs, nconss) );
   if( nobjvars > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->liftbound, nobjvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solver->liftidxs, nobjvars) );
   }
   else
   {
      solver->liftbound = NULL;
      solver->liftidxs = NULL;
   }

   /* initialize bound solution */
   if( sol == NULL )
   {
      for( i = 0; i < nvars; ++i )
      {
         solver->incumbentassignment[i] = !SCIPisInfinity(scip, -problem->vars[i].lb) ? problem->vars[i].lb
               : !SCIPisInfinity(scip, problem->vars[i].ub) ? problem->vars[i].ub : 0.0;
      }
   }
   /* initialize best solution */
   else
   {
      if( heurtiming == SCIP_HEURTIMING_BEFOREPRESOL )
      {
         vars = SCIPgetVars(scip);

         for( i = 0; i < nvars; ++i )
            solver->incumbentassignment[problem->varidxs[i]] = SCIPgetSolVal(scip, sol, vars[i]);
      }
      else
      {
         SCIP_CALL( SCIPgetLPColsData(scip, &cols, NULL) );

         for( i = 0; i < nvars; ++i )
            solver->incumbentassignment[problem->varidxs[i]] = SCIPgetSolVal(scip, sol, SCIPcolGetVar(cols[i]));
      }

      for( i = 0; i < nvars; ++i )
      {
         if( solver->incumbentassignment[i] < problem->vars[i].lb )
            solver->incumbentassignment[i] = problem->vars[i].lb;
         else if( solver->incumbentassignment[i] > problem->vars[i].ub )
            solver->incumbentassignment[i] = problem->vars[i].ub;
         else if( i < problem->nintvars )
            solver->incumbentassignment[i] = SCIPfeasRound(scip, solver->incumbentassignment[i]);
      }
   }

   /* init per-variable dynamic arrays */
   BMSclearMemoryArray(solver->allowincstep, nvars);
   BMSclearMemoryArray(solver->allowdecstep, nvars);
   BMSclearMemoryArray(solver->lastincstep, nvars);
   BMSclearMemoryArray(solver->lastdecstep, nvars);

   /* compute objective value with target flag */
   lsSolverRecomputeObjective(scip, solver);
   lsSolverUpdateObjective(scip, solver);

   /* compute constraint activities with violation indicators */
   for( i = 0; i < nconss; ++i )
   {
      solver->weight[i] = 1;
      solver->satidxs[solver->nsat] = i;
      solver->idxs[i] = -(++solver->nsat);
      lsSolverRecomputeConstraint(scip, solver, i);
      lsSolverUpdateConstraint(scip, solver, i);
   }

   /* enter optimality mode if initial assignment feasible */
   solver->optimality = solver->nunsat == 0;

   /* initialize the effort state */
   solver->previolations = solver->nunsat;
   solver->preobjective = solver->incumbentobjective;
   solver->effort = 0;

   /* use upper bound to allow initial improvement */
   if( solver->optimality )
   {
      upperbound = SCIPgetUpperbound(scip);

      if( solver->preobjective < upperbound )
         solver->preobjective = upperbound;
   }

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
      SCIPfreeBlockMemoryArray(scip, &solver->liftidxs, nobjvars);
      SCIPfreeBlockMemoryArray(scip, &solver->liftbound, nobjvars);
   }
   SCIPfreeBlockMemoryArray(scip, &solver->sampledidxs, nconss);
   SCIPfreeBlockMemoryArray(scip, &solver->sampledconstrs, nconss);
   SCIPfreeBlockMemoryArray(scip, &solver->scoreidxs, problem->nbinvars);
   SCIPfreeBlockMemoryArray(scip, &solver->scoretable, nvars);
   SCIPfreeBlockMemoryArray(scip, &solver->neighborvalues, nnonzeros);
   SCIPfreeBlockMemoryArray(scip, &solver->neighborvaridxs, nnonzeros);
   SCIPfreeBlockMemoryArray(scip, &solver->lastdecstep, nvars);
   SCIPfreeBlockMemoryArray(scip, &solver->lastincstep, nvars);
   SCIPfreeBlockMemoryArray(scip, &solver->allowdecstep, nvars);
   SCIPfreeBlockMemoryArray(scip, &solver->allowincstep, nvars);
   SCIPfreeBlockMemoryArray(scip, &solver->weight, nconss);
   SCIPfreeBlockMemoryArray(scip, &solver->consnmoves, nconss);
   SCIPfreeBlockMemoryArray(scip, &solver->idxs, nconss);
   SCIPfreeBlockMemoryArray(scip, &solver->satidxs, nconss);
   SCIPfreeBlockMemoryArray(scip, &solver->unsatidxs, nconss);
   SCIPfreeBlockMemoryArray(scip, &solver->act, nconss);
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
   SCIP_Bool recompute;
   int i;

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

   /* update objective value with target flag */
   if( var->objidx >= 0 )
   {
      ++solver->objnmoves;
      recompute = heurdata->recomputefreq > 0 && solver->objnmoves >= heurdata->recomputefreq;

      if( !recompute )
      {
         solver->incumbentobjective += var->obj * (newvalue - oldvalue);
         ++solver->effort;

         if( heurdata->recomputefreq >= 0
            && solver->incumbentobjective > solver->objtarget
            && !SCIPisPositive(scip, 0.5 * (solver->incumbentobjective - solver->objtarget)) )
            recompute = TRUE;
      }

      if( recompute )
      {
         lsSolverRecomputeObjective(scip, solver);
         solver->effort += problem->nobjvars;
      }

      lsSolverUpdateObjective(scip, solver);
   }

   /* update constraint activities with violation indicators */
   for( i = 0; i < var->ncoeffs; ++i )
   {
      constridx = var->coeffs[i].idx;
      cons = problem->conss + constridx;
      ++solver->consnmoves[constridx];
      recompute = heurdata->recomputefreq > 0 && solver->consnmoves[constridx] >= heurdata->recomputefreq;

      if( !recompute )
      {
         solver->act[constridx] += var->coeffs[i].coeff * (newvalue - oldvalue);
         ++solver->effort;

         if( heurdata->recomputefreq >= 0
            && ( !SCIPisInfinity(scip, -cons->mhs) && solver->act[constridx] <= cons->mhs
            ? solver->act[constridx] < cons->lhs
            && !SCIPisFeasNegative(scip, solver->act[constridx] - cons->lhs)
            : solver->act[constridx] > cons->rhs
            && !SCIPisFeasPositive(scip, solver->act[constridx] - cons->rhs) ) )
            recompute = TRUE;
      }

      if( recompute )
      {
         lsSolverRecomputeConstraint(scip, solver, constridx);
         solver->effort += cons->ncoeffs;
      }

      lsSolverUpdateConstraint(scip, solver, constridx);
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

/** computes break move value for an objective term */
static
SCIP_Real getBreakMove(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   int                   termidx             /**< term index within the objective */
   )
{
   LS_PROBLEM* problem;
   LS_VAR* var;
   LS_CONSTRAINT* liftcons;
   SCIP_Real value;
   SCIP_Real coeff;
   SCIP_Real residual;
   SCIP_Real movevalue;
   SCIP_Real liftcoeff;
   SCIP_Real liftresidual;
   SCIP_Real liftvalue;
   SCIP_Bool decrease;
   SCIP_Bool lifttolhs;
   int varidx;
   int liftidx;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(!SCIPisInfinity(scip, solver->objtarget));

   problem = solver->problem;
   assert(termidx >= 0);
   assert(termidx < problem->nobjvars);
   varidx = problem->objvaridxs[termidx];
   assert(varidx >= 0);
   assert(varidx < problem->nvars);
   var = problem->vars + varidx;
   coeff = var->obj;
   assert(coeff != 0.0); /*lint !e777*/
   value = solver->incumbentassignment[varidx];
   residual = solver->incumbentobjective - value * coeff;
   movevalue = (solver->objtarget - residual) / coeff;
   decrease = coeff >= 0.0;
   liftidx = solver->liftidxs[termidx];
   assert(liftidx >= 0);
   assert(liftidx < var->ncoeffs);
   liftcons = problem->conss + var->coeffs[liftidx].idx;
   liftcoeff = var->coeffs[liftidx].coeff;
   assert(liftcoeff != 0.0); /*lint !e777*/
   liftresidual = solver->act[var->coeffs[liftidx].idx] - value * liftcoeff;
   lifttolhs = decrease == (liftcoeff >= 0.0);
   assert(!lifttolhs || !SCIPisInfinity(scip, -liftcons->lhs));
   liftvalue = ((lifttolhs ? liftcons->lhs - SCIPfeastol(scip) : liftcons->rhs + SCIPfeastol(scip)) - liftresidual) / liftcoeff;

   /* better and breaking */
   if( decrease ? movevalue > liftvalue : movevalue < liftvalue )
      movevalue = liftvalue;

   /* bound and round */
   if( movevalue < var->lb )
      movevalue = var->lb;
   else if( movevalue > var->ub )
      movevalue = var->ub;
   else if( varidx < problem->nintvars )
      movevalue = round(movevalue);

   /* improve a step if target unreached */
   if( decrease ? movevalue != var->lb : movevalue != var->ub ) /*lint !e777*/
   {
      if( movevalue == value || SCIPisGT(scip, movevalue * coeff + residual, solver->objtarget) /*lint !e777*/
         || ( lifttolhs ? !SCIPisFeasNegative(scip, 2.0 * (movevalue * liftcoeff + liftresidual - liftcons->lhs))
                        : !SCIPisFeasPositive(scip, 2.0 * (movevalue * liftcoeff + liftresidual - liftcons->rhs)) ) )
         return varidx < problem->nintvars
               ? movevalue + (decrease ? -1.0 : 1.0)
               : nextafter(movevalue, decrease ? -(SCIP_Real)INFINITY : (SCIP_Real)INFINITY);
   }
   /* only move if lift is broken */
   else
   {
      if( lifttolhs ? !SCIPisFeasNegative(scip, 2.0 * (movevalue * liftcoeff + liftresidual - liftcons->lhs))
                    : !SCIPisFeasPositive(scip, 2.0 * (movevalue * liftcoeff + liftresidual - liftcons->rhs)) )
         return value;
   }

   return movevalue;
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
   assert(termidx >= 0);
   assert(termidx < cons->ncoeffs);
   assert(!tolhs || !SCIPisInfinity(scip, -cons->lhs));
   varidx = cons->coeffs[termidx].idx;
   assert(varidx >= 0);
   assert(varidx < problem->nvars);
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
      if( solver->idxs[constridx] < 0 )
         return movevalue;

      moveactivity = solver->act[constridx];
   }
   /* get move activity */
   else
   {
      moveactivity = movevalue * coeff + residual;

      if( tolhs ? !SCIPisFeasNegative(scip, 2.0 * (moveactivity - cons->lhs))
                : !SCIPisFeasPositive(scip, 2.0 * (moveactivity - cons->rhs)) )
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

      ++solver->effort;
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
      previol = 2.0 * (!SCIPisInfinity(scip, -mhs) && act <= mhs ? lhs - act : act - rhs);
      newviol = 2.0 * (!SCIPisInfinity(scip, -mhs) && newact <= mhs ? lhs - newact : newact - rhs);

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

   solver->effort += var->ncoeffs;

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

   /* increase objective weight if feasible */
   if( solver->nunsat == 0 )
      ++solver->objweight;

   /* increase constraint weights if violated */
   for( i = 0; i < solver->nunsat; ++i )
   {
      constridx = solver->unsatidxs[i];
      ++solver->weight[constridx];
   }
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

   /* decrease objective weight if better */
   if( solver->optimality && solver->objweight > 0 && !solver->objunsat
      && !SCIPisInfinity(scip, solver->incumbentobjective) )
      --solver->objweight;

   /* decrease constraint weights if feasible */
   for( i = 0; i < problem->nconss; ++i )
   {
      if( solver->idxs[i] < 0 && solver->weight[i] > 0 )
         --solver->weight[i];
   }
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
         && SCIPisLE(scip, solver->act[constridx], cons->mhs) == (solver->idxs[constridx] >= 0);

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

   solver->effort += cons->ncoeffs;
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
      /* skip free variable */
      if( solver->liftidxs[i] == -1 )
         continue;

      varidx = problem->objvaridxs[i];
      value = solver->incumbentassignment[varidx];
      movevalue = getBreakMove(scip, solver, i);

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

   solver->effort += problem->nobjvars;
}

/** scores neighbors and selects the best move */
static
SCIP_Bool selectBestNeighbor(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_Longint          minscore,           /**< minimum acceptable score */
   int                   budget,             /**< max neighbors to score */
   int*                  bestvaridx,         /**< pointer to store best variable */
   SCIP_Real*            bestvalue           /**< pointer to store best move value */
   )
{
   LS_PROBLEM* problem;
   SCIP_Real movevalue;
   SCIP_Longint bestscore;
   SCIP_Longint bestsubscore;
   SCIP_Longint score;
   int nsample;
   int idx;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(heurdata != NULL);
   assert(bestvaridx != NULL);
   assert(bestvalue != NULL);

   problem = solver->problem;
   bestscore = minscore;
   bestsubscore = 0;
   *bestvaridx = -1;
   *bestvalue = 0.0;

   /* select smaller side */
   nsample = solver->neighborsize - budget;
   if( nsample > budget )
      nsample = budget;
   else if( nsample < 0 )
      nsample = 0;

   /* partial Fisher-Yates sampling */
   for( i = 0; i < nsample; ++i )
   {
      idx = SCIPrandomGetInt(heurdata->randnumgen, i, solver->neighborsize - 1);
      SCIPswapInts(&solver->neighborvaridxs[i], &solver->neighborvaridxs[idx]);
      SCIPswapReals(&solver->neighborvalues[i], &solver->neighborvalues[idx]);
   }

   /* score shuffled range */
   if( nsample == budget )
      i = 0;
   else
      nsample = solver->neighborsize;
   while( i < nsample )
   {
      idx = solver->neighborvaridxs[i];
      movevalue = solver->neighborvalues[i];
      ++i;

      /* dedup binaries */
      if( idx < problem->nbinvars )
      {
         if( solver->scoretable[idx] )
            continue;
         solver->scoretable[idx] = TRUE;
         solver->scoreidxs[solver->nscoreidxs++] = idx;
      }

      score = getTightScore(scip, solver, idx, movevalue);

      if( score > bestscore
         || ( score == bestscore && solver->subscore > bestsubscore ) )
      {
         bestscore = score;
         bestsubscore = solver->subscore;
         *bestvaridx = idx;
         *bestvalue = movevalue;
      }
   }

   return *bestvaridx >= 0;
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
   SCIP_Real bestvalue;
   int bestvaridx;
   int nsample;
   int idx;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);

   *result = FALSE;
   solver->neighborsize = 0;

   /* select smaller side */
   nsample = solver->nunsat - heurdata->sampleunsat;
   if( nsample > heurdata->sampleunsat )
      nsample = heurdata->sampleunsat;
   else if( nsample < 0 )
      nsample = 0;

   /* partial Fisher-Yates sampling */
   for( i = 0; i < nsample; ++i )
   {
      idx = SCIPrandomGetInt(heurdata->randnumgen, i, solver->nunsat - 1);
      SCIPswapInts(&solver->unsatidxs[i], &solver->unsatidxs[idx]);
      solver->idxs[solver->unsatidxs[i]] = i;
      solver->idxs[solver->unsatidxs[idx]] = idx;
   }

   /* collect shuffled range */
   if( nsample == heurdata->sampleunsat )
      i = 0;
   else
      nsample = solver->nunsat;
   while( i < nsample )
      collectConstraintNeighbors(scip, solver, solver->unsatidxs[i++]);

   /* collect objective terms when the objective constraint is violated */
   if( solver->optimality && solver->objunsat )
      collectObjectiveNeighbors(scip, solver, FALSE);

   /* score and select within budget */
   if( selectBestNeighbor(scip, solver, heurdata, (SCIP_Longint)0,
      solver->optimality ? heurdata->bmsunsatfeas : heurdata->bmsunsatinfeas, &bestvaridx, &bestvalue) )
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
   int bestvaridx;
   SCIP_Real bestvalue;
   int constridx;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);

   *result = FALSE;
   solver->neighborsize = 0;

   if( !solver->optimality || solver->nsat == 0 )
      return SCIP_OKAY;

   /* sample satisfied constraints from satisfied set and skip duplicates */
   solver->nsampled = 0;
   for( i = 0; i < heurdata->samplesat; ++i )
   {
      /**@todo sample without replacement like unsat */
      constridx = solver->satidxs[SCIPrandomGetInt(heurdata->randnumgen, 0, solver->nsat - 1)];
      assert(constridx >= 0);
      assert(constridx < solver->problem->nconss);

      /* skip already sampled */
      if( solver->sampledconstrs[constridx] )
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

   /* score and select within budget */
   if( selectBestNeighbor(scip, solver, heurdata, (SCIP_Longint)0, heurdata->bmssat, &bestvaridx, &bestvalue) )
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
   bestsubscore = 0;
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
         || ( score == bestscore && solver->subscore > bestsubscore ) )
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

      solver->effort += cons->ncoeffs;
   }

   /* collect objective terms when the objective constraint is violated */
   if( solver->optimality && solver->objunsat )
      collectObjectiveNeighbors(scip, solver, TRUE);

   /* score and select within budget accepting negative scores */
   if( selectBestNeighbor(scip, solver, heurdata, -SCIP_LONGINT_MAX, heurdata->bmsrandom, &bestvaridx, &bestvalue) )
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

   if( solver->nunsat > 0 )
   {
      *result = FALSE;

      return SCIP_OKAY;
   }

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
      if( !solver->iskeepfeas || solver->liftidxs[i] == -1 )
      {
         solver->liftidxs[i] = -1;

         if( objcoeff >= 0.0 )
         {
            /* at best bound */
            if( SCIPisGE(scip, var->lb, value) )
            {
               solver->liftbound[i] = value;
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
                  if( SCIPisLT(scip, solver->liftbound[i], movevalue) )
                     solver->liftidxs[i] = j;

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
                  if( SCIPisGT(scip, solver->liftbound[i], movevalue) )
                     solver->liftidxs[i] = j;

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

         solver->effort += j;
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

   solver->effort += problem->nobjvars;

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
                  solver->liftidxs[objpos] = -1;
            }
         }
      }

      solver->iskeepfeas = TRUE;
      *result = TRUE;
   }
   else
   {
      solver->iskeepfeas = FALSE;
      *result = FALSE;
   }

   return SCIP_OKAY;
}


/*
 * SCIP integration
 */

/** sets up objective data on variables */
static
SCIP_RETCODE setupObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   LS_PROBLEM*           problem,            /**< problem */
   SCIP_Bool             global              /**< is problem global? */
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
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->objvaridxs, nobjvars) );
      problem->locallowerbound = global ? SCIPgetLowerbound(scip) : SCIPgetLocalLowerbound(scip);
      pos = 0;

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
   SCIP_CALL( setupObjective(scip, problem, TRUE) );

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
   SCIP_CALL( setupObjective(scip, problem, FALSE) );

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

/** stores incumbent in given solution and tries in SCIP */
static
SCIP_RETCODE checkIncumbentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURTIMING       heurtiming,         /**< heuristic timing mask */
   LS_SOLVER*            solver,             /**< solver */
   SCIP_SOL*             sol,                /**< solution to store the incumbent assignment */
   SCIP_Bool*            stored              /**< pointer to store whether the solution was accepted */
   )
{
   LS_PROBLEM* problem;
   SCIP_VAR** vars;
   SCIP_COL** cols;
   int nvars;
   int ncols;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(sol != NULL);
   assert(stored != NULL);

   problem = solver->problem;

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

   SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, FALSE, TRUE, stored) );

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
   SCIP_Bool*            terminate           /**< pointer to store termination flag */
   )
{
   LS_PROBLEM* problem;
   SCIP_SOL* sol;
   SCIP_Bool stored;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(heurdata != NULL);
   assert(heur != NULL);
   assert(terminate != NULL);

   /* check effort budget if without significant improvement */
   if( solver->nunsat >= (1.0 - heurdata->mindecreaserel) * solver->previolations
      && ( solver->nunsat > solver->previolations
      || SCIPisInfinity(scip, solver->incumbentobjective) || SCIPisInfinity(scip, -solver->preobjective)
      || SCIPisGE(scip, solver->incumbentobjective, solver->preobjective) ) )
   {
      if( heurdata->maxeffort >= 0 && solver->effort >= heurdata->maxeffort )
      {
         SCIPdebugMsg(scip, "Local: quitting, no improvement\n");
         *terminate = TRUE;
      }
      else
         *terminate = FALSE;
   }
   /* submit feasible solution and reset effort state */
   else
   {
      if( solver->nunsat == 0 && !solver->objunsat && !SCIPisInfinity(scip, solver->incumbentobjective) )
      {
         SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
         SCIP_CALL( checkIncumbentSol(scip, heurtiming, solver, sol, &stored) );
         problem = solver->problem;

         /* recover rejection by recomputation or restart */
         if( !stored )
         {
            if( heurdata->verbosity >= 1 )
               SCIPinfoMessage(scip, NULL, "Local search: rejected objective %g\n", solver->incumbentobjective);

            if( solver->objnmoves != 0 )
            {
               lsSolverRecomputeObjective(scip, solver);
               lsSolverUpdateObjective(scip, solver);
               solver->effort += problem->nobjvars;
            }

            for( i = 0; i < problem->nconss; ++i )
            {
               if( solver->consnmoves[i] != 0 )
               {
                  lsSolverRecomputeConstraint(scip, solver, i);
                  lsSolverUpdateConstraint(scip, solver, i);
                  solver->effort += problem->conss[i].ncoeffs;
               }
            }

            *terminate = solver->nunsat == 0 && !solver->objunsat && !SCIPisInfinity(scip, solver->incumbentobjective);

            if( heurdata->useterminalsol && *terminate )
            {
               assert(heurdata->terminalsol == NULL);
               heurdata->terminalsol = sol;
            }
            else
            {
               SCIP_CALL( SCIPfreeSol(scip, &sol) );
            }

            return SCIP_OKAY;
         }

         if( heurdata->verbosity >= 1 )
            SCIPinfoMessage(scip, NULL, "Local search: accepted objective %g\n", solver->incumbentobjective);

         SCIP_CALL( SCIPfreeSol(scip, &sol) );
         ++solver->nsols;

         if( (heurdata->maxsols >= 0 && solver->nsols >= heurdata->maxsols)
            || SCIPisInfinity(scip, -solver->incumbentobjective)
            || SCIPisLE(scip, solver->incumbentobjective, problem->locallowerbound) )
         {
            *terminate = TRUE;

            return SCIP_OKAY;
         }

         solver->objtarget = SCIPisObjIntegral(scip)
               ? floor(SCIPgetCutoffbound(scip)) : SCIPgetCutoffbound(scip) - SCIPcutoffbounddelta(scip);
         lsSolverUpdateObjective(scip, solver);
      }

      *terminate = FALSE;

      /* reset effort state */
      solver->previolations = solver->nunsat;
      solver->preobjective = solver->incumbentobjective;
      solver->effort = 0;
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
   SCIP_Bool terminate;
   SCIP_Bool success;
   SCIP_SOL* bestsol;
   SCIP_Bool moveresult;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* create problem structure */
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

   /* set effort budget */
   heurdata->maxeffort = heurdata->maxeffortfac >= 0
         ? heurdata->maxeffortfac * ((SCIP_Longint)problem->nnonzeros + (SCIP_Longint)problem->nobjvars)
         : (SCIP_Longint)-1;

   /* get initial solution */
   bestsol = heurdata->terminalsol != NULL ? heurdata->terminalsol
         : heurdata->useinitialsol ? SCIPgetBestSol(scip) : NULL;

   /* create solver structure */
   SCIP_CALL( lsSolverCreate(scip, heurtiming, bestsol, problem, &solver) );

   /* free terminal solution */
   if( heurdata->terminalsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &heurdata->terminalsol) );
   }

   if( heurdata->verbosity >= 1 )
   {
      SCIPinfoMessage(scip, NULL, "Local search: starting. nvars=%d nconss=%d nnonzeros=%d\n",
         problem->nvars, problem->nconss, problem->nnonzeros);
   }

   /* main search loop */
   while( !SCIPisStopped(scip) )
   {
      /* lift feasible solution to local optimality */
      SCIP_CALL( lsSolverLiftMove(scip, solver, heurdata, &moveresult) );
      if( moveresult )
         goto STEP;

      /* reset effort budget, report feasible solution, and determine termination */
      SCIP_CALL( lsCheckTermination(scip, solver, heurdata, heur, heurtiming, &terminate) );
      if( terminate )
         break;

      /* search for feasible or better solution */
      SCIP_CALL( lsSolverUnsatTightMove(scip, solver, heurdata, &moveresult) );
      if( moveresult )
         goto STEP;

      SCIP_CALL( lsSolverSatTightMove(scip, solver, heurdata, &moveresult) );
      if( moveresult )
         goto STEP;

      SCIP_CALL( lsSolverFlipMove(scip, solver, heurdata, &moveresult) );
      if( moveresult )
         goto STEP;

      cleanScoreTable(solver);

      if( SCIPrandomGetInt(heurdata->randnumgen, 0, 9999) > heurdata->smoothprob )
         lsSolverUpdateWeight(scip, solver);
      else
         lsSolverSmoothWeight(scip, solver);

      SCIP_CALL( lsSolverRandomTightMove(scip, solver, heurdata) );

   STEP:
      cleanScoreTable(solver);

      SCIPdebugMsg(scip, "Local: step=%d nunsat=%d effort=%" SCIP_LONGINT_FORMAT " obj=%g phase=search\n",
            solver->curstep, solver->nunsat, solver->effort, solver->incumbentobjective);

      ++solver->curstep;
   }

   /* report found solution */
   if( solver->nsols > 0 )
      *result = SCIP_FOUNDSOL;

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

   if( heurdata->terminalsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &heurdata->terminalsol) );
   }

   SCIPfreeRandom(scip, &heurdata->randnumgen);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecLocalsearch)
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

   SCIPdebugMsg(scip, "Executing Local Search\n");

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
   heurdata->terminalsol = NULL;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecLocalsearch, heurdata) );

   assert(heur != NULL);

   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyLocalsearch) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeLocalsearch) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitLocalsearch) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitLocalsearch) );

   /* add heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/smoothprob",
         "smooth probability (out of 10000)",
         &heurdata->smoothprob, FALSE, DEFAULT_SMOOTHPROB, 0, 10000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/tabubase",
         "minimum tabu tenure",
         &heurdata->tabubase, FALSE, DEFAULT_TABUBASE, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/tabuvariation",
         "random tabu variation",
         &heurdata->tabuvariation, FALSE, DEFAULT_TABUVARIATION, 0, INT_MAX, NULL, NULL) );
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
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/verbosity",
         "verbosity level of the local search solver",
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
         "heuristics/" HEUR_NAME "/useterminalsol",
         "should previous terminal solution be used?",
         &heurdata->useterminalsol, FALSE, DEFAULT_USETERMINALSOL, NULL, NULL) );
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

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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
 * For more details see:
 *  Luteberget, Bjørnar, and Giorgio Sartor. "Feasibility Jump: an LP-free Lagrangian MIP heuristic."
 *  Mathematical Programming Computation 15.2 (2023): 365-388.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <math.h>

#include "blockmemshell/memory.h"
#include "scip/heur_feasjump.h"
#include "scip/pub_event.h"
#include "scip/pub_heur.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_sol.h"
#include "scip/pub_var.h"
#include "scip/pub_cons.h"
#include "scip/scip_event.h"
#include "scip/scip_general.h"
#include "scip/scip_heur.h"
#include "scip/scip_cons.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_randnumgen.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"
#include "scip/scip_timing.h"

#define HEUR_NAME             "feasjump"
#define HEUR_DESC             "feasibility jump"
#define HEUR_DISPCHAR         SCIP_HEURDISPCHAR_LNS
#define HEUR_PRIORITY         900
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFOREPRESOL | SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_RANDSEED       29        /**< default random seed */
#define DEFAULT_WEIGHTDECAY    1.0       /**< decay factor for the weight update */
#define DEFAULT_MAXEFFORT      50000000  /**< maximum effort spend in feasibility jump */
#define DEFAULT_VERBOSITY      0         /**< verbosity level of the feasibility jump solver */
#define DEFAULT_ONLYWITHOUTSOL FALSE     /**< should the feasibility jump be called only if no solution exists? */
#define DEFAULT_USEINITIALSOL  FALSE     /**< should the feasibility jump use the best known solution as
                                          *   initialization? */
#define DEFAULT_ITERATIONS     1000      /**< number of iterations to check if the heuristic reduces the total
                                          *   number of violated constraints */
#define DEFAULT_MINPERCENTDECREASE   10  /**< minimum percentage decrease for the number of violated constraints */
#define DEFAULT_MAXSOLS        1         /**< maximum number of solutions to find */
#define DEFAULT_BEFORE_PRESOL  FALSE     /**< should the heuristic be called before presolving? */
#define DEFAULT_SAMPLESIZE     25        /**< number of candidates to sample in tournament selection */
#define DEFAULT_RANDPROBGOOD   0.001     /**< probability of random selection from good variables */
#define DEFAULT_RANDPROBVIOL   0.01      /**< probability of random selection from violated constraint */
#define DEFAULT_CALLBACKEFFORT 500000    /**< effort between callbacks */
#define DEFAULT_WEIGHTSCALETHRESHOLD 1.0e20 /**< threshold for rescaling weights to avoid overflow */

/* row types */
enum FJ_Rowtype
{
   FJ_EQUAL = 0,
   FJ_LTE   = 1,
   FJ_GTE   = 2
};
typedef enum FJ_Rowtype FJ_ROWTYPE;

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
   int                   coeffssize;         /**< size of coeffs array */
};
typedef struct FJ_Var FJ_VAR;

/* constraint structure */
struct FJ_Constraint
{
   FJ_ROWTYPE            sense;              /**< constraint sense */
   SCIP_Real             rhs;                /**< right hand side */
   FJ_IDXCOEFF*          coeffs;             /**< coefficients */
   int                   ncoeffs;            /**< number of coefficients */
   int                   coeffssize;         /**< size of coeffs array */
   SCIP_Real             weight;             /**< constraint weight */
   SCIP_Real             incumbentlhs;       /**< LHS value at incumbent */
   int                   violatedidx;        /**< index in violated constraints list, -1 if satisfied */
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
   SCIP*                 scip;               /**< SCIP data structure */
   FJ_VAR*               vars;               /**< variables */
   int                   nvars;              /**< number of variables */
   int                   varssize;           /**< size of vars array */
   FJ_CONSTRAINT*        constraints;        /**< constraints */
   int                   nconstraints;       /**< number of constraints */
   int                   constraintssize;    /**< size of constraints array */
   SCIP_Real*            incumbentassignment;/**< incumbent assignment */
   int*                  violatedconstraints;/**< violated constraints indices */
   int                   nviolated;          /**< number of violated constraints */
   int                   violatedsize;       /**< size of violated array */
   SCIP_Bool             usedrelaxcontinuous;/**< was continuous variable relaxation used? */
   int                   nnonzeros;          /**< number of nonzeros */
   SCIP_Real             incumbentobjective; /**< incumbent objective value */
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
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   SCIP_Real             bestobjective;      /**< best objective value found */
   SCIP_Real             objectiveweight;    /**< weight for objective function */
   int                   bestviolationscore; /**< best violation score */
   int                   effortatlastcallback;/**< effort at last callback */
   int                   effortatlastimprovement;/**< effort at last improvement */
   int                   totaleffort;        /**< total effort */
   int                   violationsdecrease; /**< decrease in violations */
   int                   prevviolations;     /**< previous number of violations */
   int                   percentdecrease;    /**< percent decrease in violations */
   int                   iterations;         /**< iteration count */
   SCIP_Real             weightupdatedecay;  /**< weight update decay */
   SCIP_Real             weightupdateincrement;/**< weight update increment */
   int                   nbumps;             /**< number of weight bumps */
   int                   verbosity;          /**< verbosity level */
   int                   checkiterations;    /**< number of iterations to check progress */
   int                   samplesize;         /**< number of candidates to sample in tournament selection */
   SCIP_Real             randprobgood;       /**< probability of random selection from good variables */
   SCIP_Real             randprobviol;       /**< probability of random selection from violated constraint */
   int                   callbackeffort;     /**< effort between callbacks */
   SCIP_Real             weightscalethreshold; /**< threshold for rescaling weights to avoid overflow */
};
typedef struct FJ_Solver FJ_SOLVER;

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   SCIP_Real             weightupdatedecay;  /**< decay factor for the weight update */
   int                   maxeffort;          /**< max effort in feasibility jump */
   int                   verbosity;          /**< verbosity level of the feasibility jump solver */
   SCIP_Bool             onlywithoutsol;     /**< should the heuristic only be called, if no primal solution exists? */
   SCIP_Bool             useinitialsol;      /**< should the heuristic use the best known solution as initialization? */
   SCIP_Bool             beforepresol;       /**< should the heuristic be called before presolving? */
   int                   iterations;         /**< number of iterations to check if the heuristic reduces the
                                              *   total number of violated constraints */
   int                   mindecrease;        /**< minimum percentage decrease for the number of violated
                                              *   constraints */
   int                   maxsols;            /**< maximum number of solutions to find */
   int                   samplesize;         /**< number of candidates to sample in tournament selection */
   SCIP_Real             randprobgood;       /**< probability of random selection from good variables */
   SCIP_Real             randprobviol;       /**< probability of random selection from violated constraint */
   int                   callbackeffort;     /**< effort between callbacks */
   SCIP_Real             weightscalethreshold; /**< threshold for rescaling weights to avoid overflow */
};


/*
 * Local methods for FJ_Problem
 */

/** computes the constraint's contribution to the feasibility score */
static
SCIP_Real fjConstraintScore(
   FJ_CONSTRAINT*        constraint,         /**< constraint */
   SCIP_Real             lhs                 /**< LHS value */
   )
{
   assert(constraint != NULL);

   lhs -= constraint->rhs;

   if( constraint->sense == FJ_EQUAL )
      return -REALABS(lhs);
   else if( constraint->sense == FJ_LTE )
      return -MAX(lhs, 0.0);
   else
      return MIN(lhs, 0.0);
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
   prob->scip = scip;
   prob->vars = NULL;
   prob->nvars = 0;
   prob->varssize = 0;
   prob->constraints = NULL;
   prob->nconstraints = 0;
   prob->constraintssize = 0;
   prob->incumbentassignment = NULL;
   prob->violatedconstraints = NULL;
   prob->nviolated = 0;
   prob->violatedsize = 0;
   prob->usedrelaxcontinuous = FALSE;
   prob->nnonzeros = 0;
   prob->incumbentobjective = 0.0;

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

   /* free variable data */
   if( prob->vars != NULL )
   {
      for( i = 0; i < prob->nvars; ++i )
      {
         if( prob->vars[i].coeffs != NULL )
         {
            SCIPfreeBlockMemoryArray(scip, &prob->vars[i].coeffs, prob->vars[i].coeffssize);
         }
      }
      SCIPfreeBlockMemoryArray(scip, &prob->vars, prob->varssize);
   }

   /* free constraint data */
   if( prob->constraints != NULL )
   {
      for( i = 0; i < prob->nconstraints; ++i )
      {
         if( prob->constraints[i].coeffs != NULL )
         {
            SCIPfreeBlockMemoryArray(scip, &prob->constraints[i].coeffs, prob->constraints[i].coeffssize);
         }
      }
      SCIPfreeBlockMemoryArray(scip, &prob->constraints, prob->constraintssize);
   }

   if( prob->incumbentassignment != NULL )
      SCIPfreeBlockMemoryArray(scip, &prob->incumbentassignment, prob->nvars);

   if( prob->violatedconstraints != NULL )
      SCIPfreeBlockMemoryArray(scip, &prob->violatedconstraints, prob->violatedsize);

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
   assert(idx != NULL);

   /* ensure capacity */
   if( problem->nvars >= problem->varssize )
   {
      int newsize = SCIPcalcMemGrowSize(scip, problem->nvars + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &problem->vars, problem->varssize, newsize) );
      problem->varssize = newsize;
   }

   *idx = problem->nvars;
   var = &problem->vars[problem->nvars];
   var->vartype = vartype;
   var->lb = lb;
   var->ub = ub;
   var->objcoeff = objcoeff;
   var->coeffs = NULL;
   var->ncoeffs = 0;
   var->coeffssize = 0;

   problem->nvars++;

   return SCIP_OKAY;
}

/** adds a constraint to the problem */
static
SCIP_RETCODE fjProblemAddConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_PROBLEM*           problem,            /**< problem */
   FJ_ROWTYPE            sense,              /**< constraint sense */
   SCIP_Real             rhs,                /**< right hand side */
   int                   ncoeffs,            /**< number of coefficients */
   int*                  rowinds,            /**< variable indices */
   SCIP_Real*            rowcoeffs,          /**< coefficients */
   SCIP_Bool             relaxcontinuous,    /**< should continuous variables be relaxed? */
   int*                  idx                 /**< pointer to store constraint index */
   )
{
   FJ_CONSTRAINT* constraint;
   SCIP_Real newrhs;
   SCIP_Real scalar;
   int i;
   int nnewcoeffs;

   assert(scip != NULL);
   assert(problem != NULL);
   assert(rowinds != NULL || ncoeffs == 0);
   assert(rowcoeffs != NULL || ncoeffs == 0);

   if( relaxcontinuous )
      problem->usedrelaxcontinuous = TRUE;

   /* handle continuous variable relaxation for equality constraints */
   if( relaxcontinuous && sense == FJ_EQUAL )
   {
      SCIP_Bool hascontinuous = FALSE;
      for( i = 0; i < ncoeffs; ++i )
      {
         if( problem->vars[rowinds[i]].vartype == FJ_CONTINUOUS )
         {
            hascontinuous = TRUE;
            break;
         }
      }

      if( hascontinuous )
      {
         int idx1, idx2;
         SCIP_CALL( fjProblemAddConstraint(scip, problem, FJ_GTE, rhs, ncoeffs, rowinds, rowcoeffs,
               relaxcontinuous, &idx1) );
         SCIP_CALL( fjProblemAddConstraint(scip, problem, FJ_LTE, rhs, ncoeffs, rowinds, rowcoeffs,
               relaxcontinuous, &idx2) );
         if( idx != NULL )
            *idx = INT_MAX;
         return SCIP_OKAY;
      }
   }

   newrhs = rhs;
   scalar = MAX3(rhs, -rhs, 1.0);
   nnewcoeffs = 0;

   /* adjust rhs if relaxing continuous variables */
   if( relaxcontinuous )
   {
      for( i = 0; i < ncoeffs; ++i )
      {
         if( problem->vars[rowinds[i]].vartype == FJ_CONTINUOUS )
         {
            if( sense == FJ_LTE )
            {
               if( rowcoeffs[i] >= 0.0 )
               {
                  if( !SCIPisInfinity(scip, -problem->vars[rowinds[i]].lb) )
                     newrhs -= rowcoeffs[i] * problem->vars[rowinds[i]].lb;
                  else
                  {
                     if( idx != NULL )
                        *idx = INT_MAX;
                     return SCIP_OKAY;
                  }
               }
               else
               {
                  if( !SCIPisInfinity(scip, problem->vars[rowinds[i]].ub) )
                     newrhs -= rowcoeffs[i] * problem->vars[rowinds[i]].ub;
                  else
                  {
                     if( idx != NULL )
                        *idx = INT_MAX;
                     return SCIP_OKAY;
                  }
               }
            }
            else if( sense == FJ_GTE )
            {
               if( rowcoeffs[i] >= 0.0 )
               {
                  if( !SCIPisInfinity(scip, problem->vars[rowinds[i]].ub) )
                     newrhs -= rowcoeffs[i] * problem->vars[rowinds[i]].ub;
                  else
                  {
                     if( idx != NULL )
                        *idx = INT_MAX;
                     return SCIP_OKAY;
                  }
               }
               else
               {
                  if( !SCIPisInfinity(scip, -problem->vars[rowinds[i]].lb) )
                     newrhs -= rowcoeffs[i] * problem->vars[rowinds[i]].lb;
                  else
                  {
                     if( idx != NULL )
                        *idx = INT_MAX;
                     return SCIP_OKAY;
                  }
               }
            }
            else
            {
               if( idx != NULL )
                  *idx = INT_MIN;
               return SCIP_OKAY;
            }
         }
         else
            nnewcoeffs++;
      }
   }
   else
      nnewcoeffs = ncoeffs;

   /* check if constraint is trivially satisfied or infeasible */
   if( nnewcoeffs == 0 )
   {
      SCIP_Bool ok;
      if( sense == FJ_LTE )
         ok = !SCIPisFeasNegative(scip, newrhs / scalar);
      else if( sense == FJ_GTE )
         ok = !SCIPisFeasPositive(scip, newrhs / scalar);
      else
         ok = SCIPisFeasZero(scip, newrhs / scalar);

      if( idx != NULL )
         *idx = ok ? INT_MAX : INT_MIN;
      return SCIP_OKAY;
   }

   /* ensure capacity */
   if( problem->nconstraints >= problem->constraintssize )
   {
      int newsize = SCIPcalcMemGrowSize(scip, problem->nconstraints + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &problem->constraints, problem->constraintssize, newsize) );
      problem->constraintssize = newsize;
   }

   if( idx != NULL )
      *idx = problem->nconstraints;

   constraint = &problem->constraints[problem->nconstraints];
   constraint->sense = sense;
   constraint->rhs = newrhs / scalar;
   constraint->ncoeffs = 0;
   constraint->coeffssize = nnewcoeffs;
   constraint->weight = 1.0;
   constraint->incumbentlhs = 0.0;
   constraint->violatedidx = -1;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &constraint->coeffs, constraint->coeffssize) );

   /* add coefficients and update variable structures */
   for( i = 0; i < ncoeffs; ++i )
   {
      FJ_VAR* var;

      if( relaxcontinuous && problem->vars[rowinds[i]].vartype == FJ_CONTINUOUS )
         continue;

      constraint->coeffs[constraint->ncoeffs].idx = rowinds[i];
      constraint->coeffs[constraint->ncoeffs].coeff = rowcoeffs[i] / scalar;
      constraint->ncoeffs++;

      var = &problem->vars[rowinds[i]];
      if( var->ncoeffs >= var->coeffssize )
      {
         int newsize = SCIPcalcMemGrowSize(scip, var->ncoeffs + 1);
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &var->coeffs, var->coeffssize, newsize) );
         var->coeffssize = newsize;
      }
      var->coeffs[var->ncoeffs].idx = problem->nconstraints;
      var->coeffs[var->ncoeffs].coeff = rowcoeffs[i] / scalar;
      var->ncoeffs++;
   }

   problem->nnonzeros += constraint->ncoeffs;
   problem->nconstraints++;

   return SCIP_OKAY;
}

/** resets the incumbent solution */
static
SCIP_RETCODE fjProblemResetIncumbent(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_PROBLEM*           problem,            /**< problem */
   SCIP_Real*            initialvalues       /**< initial values, or NULL */
   )
{
   int i;
   int j;

   assert(scip != NULL);
   assert(problem != NULL);

   /* allocate arrays if needed */
   if( problem->incumbentassignment == NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->incumbentassignment, problem->nvars) );
   }
   if( problem->violatedconstraints == NULL )
   {
      problem->violatedsize = SCIPcalcMemGrowSize(scip, problem->nconstraints);
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->violatedconstraints, problem->violatedsize) );
   }

   /* set initial values */
   if( initialvalues != NULL )
   {
      for( i = 0; i < problem->nvars; ++i )
         problem->incumbentassignment[i] = initialvalues[i];
   }
   else
   {
      for( i = 0; i < problem->nvars; ++i )
         problem->incumbentassignment[i] = !SCIPisInfinity(scip, -problem->vars[i].lb) ? problem->vars[i].lb
               : !SCIPisInfinity(scip, problem->vars[i].ub) ? problem->vars[i].ub : 0.0;
   }

   /* reset incumbent objective */
   problem->incumbentobjective = 0.0;
   for( i = 0; i < problem->nvars; ++i )
      problem->incumbentobjective += problem->vars[i].objcoeff * problem->incumbentassignment[i];

   /* reset constraint LHSs and violated constraints list */
   problem->nviolated = 0;
   for( i = 0; i < problem->nconstraints; ++i )
   {
      FJ_CONSTRAINT* constraint = &problem->constraints[i];
      constraint->incumbentlhs = 0.0;

      for( j = 0; j < constraint->ncoeffs; ++j )
      {
         int varidx = constraint->coeffs[j].idx;
         constraint->incumbentlhs += constraint->coeffs[j].coeff * problem->incumbentassignment[varidx];
      }

      if( SCIPisFeasNegative(scip, fjConstraintScore(constraint, constraint->incumbentlhs)) )
      {
         constraint->violatedidx = problem->nviolated;
         problem->violatedconstraints[problem->nviolated] = i;
         problem->nviolated++;
      }
      else
      {
         constraint->violatedidx = -1;
      }
   }

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
   SCIP_Real varincumbent;
   SCIP_Real bestscore;
   SCIP_Real bestvalue;
   SCIP_Real currentscore;
   SCIP_Real currentvalue;
   SCIP_Real currentslope;
   SCIP_Real currentvaluestore;
   SCIP_Real currentslopestore;
   int currentindexstore;
   int i;
   int j;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(varidx >= 0 && varidx < solver->problem->nvars);

   problem = solver->problem;
   var = &problem->vars[varidx];
   varincumbent = problem->incumbentassignment[varidx];

   solver->nshiftbuffer = 0;
   currentvalue = var->lb;
   currentslope = 0.0;

   /* compute shift buffer */
   for( i = 0; i < var->ncoeffs; ++i )
   {
      int cstridx = var->coeffs[i].idx;
      SCIP_Real cellcoeff = var->coeffs[i].coeff;
      FJ_CONSTRAINT* constraint = &problem->constraints[cstridx];
      int nbounds = 0;
      SCIP_Real bounds[6];

      assert(cellcoeff != 0.0); /*lint !e777*/

      /* determine bounds based on constraint sense */
      if( constraint->sense == FJ_LTE )
      {
         bounds[nbounds++] = -SCIPinfinity(scip);
         bounds[nbounds++] = constraint->rhs;
      }
      else if( constraint->sense == FJ_GTE )
      {
         bounds[nbounds++] = constraint->rhs;
         bounds[nbounds++] = SCIPinfinity(scip);
      }
      else
      {
         bounds[nbounds++] = -SCIPinfinity(scip);
         bounds[nbounds++] = constraint->rhs;
         bounds[nbounds++] = constraint->rhs;
         bounds[nbounds++] = constraint->rhs;
         bounds[nbounds++] = constraint->rhs;
         bounds[nbounds++] = SCIPinfinity(scip);
      }

      for( j = 0; j < nbounds; j += 2 )
      {
         SCIP_Real residualincumbent = constraint->incumbentlhs - cellcoeff * varincumbent;
         SCIP_Real validrangelb;
         SCIP_Real validrangeub;

         if( cellcoeff >= 0.0 )
         {
            if( !SCIPisInfinity(scip, -bounds[j]) )
            {
               validrangelb = (bounds[j] - residualincumbent) / cellcoeff;
               if( var->vartype == FJ_INTEGER )
               {
                  validrangelb = round(validrangelb);
                  if( SCIPisFeasNegative(scip, validrangelb * cellcoeff + residualincumbent - bounds[j]) )
                     validrangelb += 1.0;
               }
               else
               {
                  if( SCIPisFeasNegative(scip, validrangelb * cellcoeff + residualincumbent - bounds[j]) )
                     validrangelb = nextafter(validrangelb, (SCIP_Real)INFINITY);
               }
            }
            else
               validrangelb = bounds[j];

            if( !SCIPisInfinity(scip, bounds[j + 1]) )
            {
               validrangeub = (bounds[j + 1] - residualincumbent) / cellcoeff;
               if( var->vartype == FJ_INTEGER )
               {
                  validrangeub = round(validrangeub);
                  if( SCIPisFeasPositive(scip, validrangeub * cellcoeff + residualincumbent - bounds[j + 1]) )
                     validrangeub -= 1.0;
               }
               else
               {
                  if( SCIPisFeasPositive(scip, validrangeub * cellcoeff + residualincumbent - bounds[j + 1]) )
                     validrangeub = nextafter(validrangeub, -(SCIP_Real)INFINITY);
               }
            }
            else
               validrangeub = bounds[j + 1];
         }
         else
         {
            if( !SCIPisInfinity(scip, bounds[j + 1]) )
            {
               validrangelb = (bounds[j + 1] - residualincumbent) / cellcoeff;
               if( var->vartype == FJ_INTEGER )
               {
                  validrangelb = round(validrangelb);
                  if( SCIPisFeasPositive(scip, validrangelb * cellcoeff + residualincumbent - bounds[j + 1]) )
                     validrangelb += 1.0;
               }
               else
               {
                  if( SCIPisFeasPositive(scip, validrangelb * cellcoeff + residualincumbent - bounds[j + 1]) )
                     validrangelb = nextafter(validrangelb, (SCIP_Real)INFINITY);
               }
            }
            else
               validrangelb = -bounds[j + 1];

            if( !SCIPisInfinity(scip, -bounds[j]) )
            {
               validrangeub = (bounds[j] - residualincumbent) / cellcoeff;
               if( var->vartype == FJ_INTEGER )
               {
                  validrangeub = round(validrangeub);
                  if( SCIPisFeasNegative(scip, validrangeub * cellcoeff + residualincumbent - bounds[j]) )
                     validrangeub -= 1.0;
               }
               else
               {
                  if( SCIPisFeasNegative(scip, validrangeub * cellcoeff + residualincumbent - bounds[j]) )
                     validrangeub = nextafter(validrangeub, -(SCIP_Real)INFINITY);
               }
            }
            else
               validrangeub = -bounds[j];
         }

         if( validrangelb > validrangeub )
            continue;

         if( validrangelb > currentvalue )
         {
            currentslope -= constraint->weight;
            if( validrangelb < var->ub )
            {
               /* ensure capacity */
               if( solver->nshiftbuffer >= solver->shiftbuffersize )
               {
                  int oldsize = solver->shiftbuffersize;
                  int newsize = SCIPcalcMemGrowSize(scip, solver->nshiftbuffer + 1);
                  int k;
                  SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &solver->shiftbuffer, oldsize, newsize) );
                  for( k = oldsize; k < newsize; ++k )
                  {
                     SCIP_CALL( SCIPallocBlockMemory(scip, &solver->shiftbuffer[k]) ); /*lint !e866*/
                  }
                  solver->shiftbuffersize = newsize;
               }
               solver->shiftbuffer[solver->nshiftbuffer]->value = validrangelb;
               solver->shiftbuffer[solver->nshiftbuffer]->weight = constraint->weight;
               solver->nshiftbuffer++;
            }
         }

         if( validrangeub <= currentvalue )
            currentslope += constraint->weight;
         else if( validrangeub < var->ub )
         {
            /* ensure capacity */
            if( solver->nshiftbuffer >= solver->shiftbuffersize )
            {
               int oldsize = solver->shiftbuffersize;
               int newsize = SCIPcalcMemGrowSize(scip, solver->nshiftbuffer + 1);
               int k;
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &solver->shiftbuffer, oldsize, newsize) );
               for( k = oldsize; k < newsize; ++k )
               {
                  SCIP_CALL( SCIPallocBlockMemory(scip, &solver->shiftbuffer[k]) ); /*lint !e866*/
               }
               solver->shiftbuffersize = newsize;
            }
            solver->shiftbuffer[solver->nshiftbuffer]->value = validrangeub;
            solver->shiftbuffer[solver->nshiftbuffer]->weight = constraint->weight;
            solver->nshiftbuffer++;
         }
      }
   }

   /* add bounds to shift buffer */
   if( solver->nshiftbuffer + 2 >= solver->shiftbuffersize )
   {
      int oldsize = solver->shiftbuffersize;
      int newsize = SCIPcalcMemGrowSize(scip, solver->nshiftbuffer + 2);
      int k;
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &solver->shiftbuffer, oldsize, newsize) );
      for( k = oldsize; k < newsize; ++k )
      {
         SCIP_CALL( SCIPallocBlockMemory(scip, &solver->shiftbuffer[k]) ); /*lint !e866*/
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
   if( SCIPisEQ(scip, currentvalue, varincumbent) )
   {
      bestscore = INFINITY;
      currentscore = 0.0;
      currentvaluestore = currentvalue;
      currentslopestore = currentslope;
      currentindexstore = i;

      if( SCIPisGT(scip, solver->shiftbuffer[solver->nshiftbuffer - 1]->value, varincumbent) )
      {
         do
         {
            ++i;
            currentscore += (solver->shiftbuffer[i]->value - currentvalue) * currentslope;
            currentslope += solver->shiftbuffer[i]->weight;
            currentvalue = solver->shiftbuffer[i]->value;
         }
         while( SCIPisLE(scip, currentvalue, varincumbent) );

         bestscore = currentscore;
         bestvalue = currentvalue;
      }

      currentscore = 0.0;
      currentvalue = currentvaluestore;
      currentslope = currentslopestore;
      i = currentindexstore;

      if( SCIPisLT(scip, solver->shiftbuffer[0]->value, varincumbent) )
      {
         do
         {
            currentslope -= solver->shiftbuffer[i]->weight;
            --i;
            currentscore += (solver->shiftbuffer[i]->value - currentvalue) * currentslope;
            currentvalue = solver->shiftbuffer[i]->value;
         }
         while( SCIPisGE(scip, currentvalue, varincumbent) );

         if( bestscore > currentscore )
            bestvalue = currentvalue;
      }
   }

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

   assert(scip != NULL);
   assert(solver != NULL);
   assert(varidx >= 0 && varidx < solver->problem->nvars);

   anygoodmoves = (solver->jumpmoves[varidx].score > 0.0);

   if( anygoodmoves && solver->goodvarssetidx[varidx] == -1 )
   {
      /* became good, add to good set */
      if( solver->ngoodvars >= solver->goodvarssize )
      {
         int newsize = SCIPcalcMemGrowSize(scip, solver->ngoodvars + 1);
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
      int lastsetidx = solver->ngoodvars - 1;
      int lastvaridx = solver->goodvarsset[lastsetidx];
      int thissetidx = solver->goodvarssetidx[varidx];

      solver->goodvarsset[thissetidx] = lastvaridx;
      solver->goodvarssetidx[lastvaridx] = thissetidx;
      solver->goodvarssetidx[varidx] = -1;
      solver->ngoodvars--;
   }

   return SCIP_OKAY;
}

/** sets a variable value and updates violated constraints */
static
SCIP_RETCODE fjProblemSetValue(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_PROBLEM*           problem,            /**< problem */
   FJ_SOLVER*            solver,             /**< solver (for move updates) */
   int                   varidx,             /**< variable index */
   SCIP_Real             newvalue,           /**< new value */
   int*                  effort              /**< pointer to effort counter */
   )
{
   SCIP_Real oldvalue;
   SCIP_Real delta;
   int i;
   int j;

   assert(scip != NULL);
   assert(problem != NULL);
   assert(varidx >= 0 && varidx < problem->nvars);
   assert(effort != NULL);

   oldvalue = problem->incumbentassignment[varidx];
   delta = newvalue - oldvalue;
   problem->incumbentassignment[varidx] = newvalue;
   problem->incumbentobjective += problem->vars[varidx].objcoeff * delta;

   /* update LHSs of all involved constraints */
   for( i = 0; i < problem->vars[varidx].ncoeffs; ++i )
   {
      int cstridx = problem->vars[varidx].coeffs[i].idx;
      SCIP_Real coeff = problem->vars[varidx].coeffs[i].coeff;
      FJ_CONSTRAINT* constraint = &problem->constraints[cstridx];
      SCIP_Real oldlhs = constraint->incumbentlhs;
      SCIP_Real newlhs = oldlhs + coeff * delta;
      SCIP_Real newcost;

      constraint->incumbentlhs = newlhs;
      newcost = fjConstraintScore(constraint, newlhs);

      /* add/remove from violated constraints list */
      if( constraint->violatedidx == -1 && SCIPisFeasNegative(scip, newcost) )
      {
         /* became violated */
         constraint->violatedidx = problem->nviolated;
         problem->violatedconstraints[problem->nviolated] = cstridx;
         problem->nviolated++;
      }
      if( constraint->violatedidx != -1 && !SCIPisFeasNegative(scip, newcost) )
      {
         /* became satisfied */
         int lastviolatedidx = problem->nviolated - 1;
         int lastconstraintidx = problem->violatedconstraints[lastviolatedidx];
         int thisviolatedidx = constraint->violatedidx;

         problem->violatedconstraints[thisviolatedidx] = lastconstraintidx;
         problem->constraints[lastconstraintidx].violatedidx = thisviolatedidx;
         constraint->violatedidx = -1;
         problem->nviolated--;
      }

      /* update moves for other variables in this constraint */
      if( solver != NULL )
      {
         *effort += constraint->ncoeffs;
         for( j = 0; j < constraint->ncoeffs; ++j )
         {
            int otheridx = constraint->coeffs[j].idx;
            if( otheridx != varidx )
            {
               FJ_MOVE* move = &solver->jumpmoves[otheridx];
               SCIP_Real otherincumbent = problem->incumbentassignment[otheridx];
               SCIP_Real othercoeff = constraint->coeffs[j].coeff;
               SCIP_Real oldmodifiedlhs = oldlhs + othercoeff * (move->value - otherincumbent);
               SCIP_Real oldscoreterm = constraint->weight * (fjConstraintScore(constraint, oldmodifiedlhs)
                     - fjConstraintScore(constraint, oldlhs));
               SCIP_Real newmodifiedlhs = newlhs + othercoeff * (move->value - otherincumbent);
               SCIP_Real newscoreterm = constraint->weight * (fjConstraintScore(constraint, newmodifiedlhs)
                     - fjConstraintScore(constraint, newlhs));

               move->score += newscoreterm - oldscoreterm;
               SCIP_CALL( fjSolverUpdateGoodMoves(scip, solver, otheridx) );
            }
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
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(varidx >= 0 && varidx < solver->problem->nvars);

   problem = solver->problem;
   var = &problem->vars[varidx];
   move = &solver->jumpmoves[varidx];

   solver->totaleffort += var->ncoeffs;
   SCIP_CALL( fjSolverUpdateJumpValue(scip, solver, varidx) );

   move->score = solver->objectiveweight * var->objcoeff * (move->value - problem->incumbentassignment[varidx]);

   for( i = 0; i < var->ncoeffs; ++i )
   {
      int cstridx = var->coeffs[i].idx;
      SCIP_Real cellcoeff = var->coeffs[i].coeff;
      FJ_CONSTRAINT* constraint = &problem->constraints[cstridx];
      SCIP_Real candidatelhs = constraint->incumbentlhs + cellcoeff * (move->value
            - problem->incumbentassignment[varidx]);

      move->score += constraint->weight * (fjConstraintScore(constraint, candidatelhs)
            - fjConstraintScore(constraint, constraint->incumbentlhs));
   }

   SCIP_CALL( fjSolverUpdateGoodMoves(scip, solver, varidx) );

   return SCIP_OKAY;
}

/** updates weights at local minimum */
static
SCIP_RETCODE fjSolverUpdateWeights(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_SOLVER*            solver              /**< solver */
   )
{
   FJ_PROBLEM* problem;
   SCIP_Bool rescaleallweights = FALSE;
   int i;
   int j;

   assert(scip != NULL);
   assert(solver != NULL);

   problem = solver->problem;

   if( solver->verbosity >= 2 )
      SCIPinfoMessage(scip, NULL, "Feasibility Jump: Reached a local minimum.\n");

   solver->nbumps++;

   if( problem->nviolated == 0 )
   {
      solver->objectiveweight += solver->weightupdateincrement;
      if( solver->objectiveweight > solver->weightscalethreshold )
         rescaleallweights = TRUE;

      solver->totaleffort += problem->nvars;
      for( i = 0; i < problem->nvars; ++i )
      {
         FJ_MOVE* move = &solver->jumpmoves[i];
         move->score += solver->weightupdateincrement * problem->vars[i].objcoeff
               * (move->value - problem->incumbentassignment[i]);
      }
   }
   else
   {
      for( i = 0; i < problem->nviolated; ++i )
      {
         int cstridx = problem->violatedconstraints[i];
         FJ_CONSTRAINT* constraint = &problem->constraints[cstridx];

         constraint->weight += solver->weightupdateincrement;
         if( constraint->weight > solver->weightscalethreshold )
            rescaleallweights = TRUE;

         solver->totaleffort += constraint->ncoeffs;
         for( j = 0; j < constraint->ncoeffs; ++j )
         {
            int cellidx = constraint->coeffs[j].idx;
            SCIP_Real cellcoeff = constraint->coeffs[j].coeff;
            FJ_MOVE* move = &solver->jumpmoves[cellidx];
            SCIP_Real candidatelhs = constraint->incumbentlhs + cellcoeff * (move->value
                  - problem->incumbentassignment[cellidx]);
            SCIP_Real diff = solver->weightupdateincrement * (fjConstraintScore(constraint, candidatelhs)
                  - fjConstraintScore(constraint, constraint->incumbentlhs));

            move->score += diff;

            SCIP_CALL( fjSolverUpdateGoodMoves(scip, solver, cellidx) );
         }
      }
   }

   solver->weightupdateincrement /= solver->weightupdatedecay;

   if( rescaleallweights )
   {
      /* TODO check if rescaling the increment makes sense */ 
      solver->weightupdateincrement /= solver->weightscalethreshold;
      solver->objectiveweight /= solver->weightscalethreshold;

      for( i = 0; i < problem->nconstraints; ++i )
         problem->constraints[i].weight /= solver->weightscalethreshold;

      solver->totaleffort += problem->nconstraints;

      for( i = 0; i < problem->nvars; ++i )
      {
         SCIP_CALL( fjSolverResetMoves(scip, solver, i) );
      }
   }

   return SCIP_OKAY;
}

/** selects a variable to modify */
static
SCIP_RETCODE fjSolverSelectVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_SOLVER*            solver,             /**< solver */
   int*                  varidx              /**< pointer to store selected variable */
   )
{
   FJ_PROBLEM* problem;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(varidx != NULL);

   problem = solver->problem;

   if( solver->ngoodvars > 0 )
   {
      SCIP_Real randval = SCIPrandomGetReal(solver->randnumgen, 0.0, 1.0);

      if( randval < solver->randprobgood )
      {
         int randomidx = SCIPrandomGetInt(solver->randnumgen, 0, solver->ngoodvars - 1);
         *varidx = solver->goodvarsset[randomidx];
         return SCIP_OKAY;
      }
      else
      {
         int samplesize = MIN(solver->samplesize, solver->ngoodvars);
         SCIP_Real bestscore = -INFINITY;
         int bestvar = -1;

         solver->totaleffort += samplesize;

         for( i = 0; i < samplesize; ++i )
         {
            int setidx = SCIPrandomGetInt(solver->randnumgen, 0, solver->ngoodvars - 1);
            int candidatevar = solver->goodvarsset[setidx];
            SCIP_Real score = solver->jumpmoves[candidatevar].score;

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
   SCIP_CALL( fjSolverUpdateWeights(scip, solver) );

   if( problem->nviolated > 0 )
   {
      int cstridx = problem->violatedconstraints[SCIPrandomGetInt(solver->randnumgen, 0, problem->nviolated - 1)];
      FJ_CONSTRAINT* constraint = &problem->constraints[cstridx];
      SCIP_Real randval = SCIPrandomGetReal(solver->randnumgen, 0.0, 1.0);

      if( randval < solver->randprobviol )
      {
         int randomidx = SCIPrandomGetInt(solver->randnumgen, 0, constraint->ncoeffs - 1);
         *varidx = constraint->coeffs[randomidx].idx;
         return SCIP_OKAY;
      }
      else
      {
         SCIP_Real bestscore = -INFINITY;
         int bestvar = -1;

         for( i = 0; i < constraint->ncoeffs; ++i )
         {
            int cellidx = constraint->coeffs[i].idx;
            SCIP_Real score = solver->jumpmoves[cellidx].score;

            if( score > bestscore )
            {
               bestscore = score;
               bestvar = cellidx;
            }
         }

         assert(bestvar >= 0);
         *varidx = bestvar;
         return SCIP_OKAY;
      }
   }

   /* fallback to random choice */
   *varidx = SCIPrandomGetInt(solver->randnumgen, 0, problem->nvars - 1);

   return SCIP_OKAY;
}

/** performs a variable move */
static
SCIP_RETCODE fjSolverDoVariableMove(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_SOLVER*            solver,             /**< solver */
   int                   varidx              /**< variable index */
   )
{
   SCIP_Real newvalue;
   int effort = 0;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(varidx >= 0 && varidx < solver->problem->nvars);

   newvalue = solver->jumpmoves[varidx].value;

   SCIP_CALL( fjProblemSetValue(scip, solver->problem, solver, varidx, newvalue, &effort) );
   solver->totaleffort += effort;

   SCIP_CALL( fjSolverResetMoves(scip, solver, varidx) );

   return SCIP_OKAY;
}

/** creates a solver */
static
SCIP_RETCODE fjSolverCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_SOLVER**           solver,             /**< pointer to solver */
   FJ_PROBLEM*           problem,            /**< problem */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   int                   verbosity,          /**< verbosity level */
   SCIP_Real             weightupdatedecay,  /**< weight update decay */
   SCIP_Real             bestobjective,      /**< best objective value */
   int                   checkiterations,    /**< iterations to check progress */
   int                   samplesize,         /**< number of candidates to sample in tournament selection */
   SCIP_Real             randprobgood,       /**< probability of random selection from good variables */
   SCIP_Real             randprobviol,       /**< probability of random selection from violated constraint */
   int                   callbackeffort,     /**< effort between callbacks */
   SCIP_Real             weightscalethreshold /**< threshold for rescaling weights to avoid overflow */
   )
{
   FJ_SOLVER* solv;
   int i;

   assert(scip != NULL);
   assert(solver != NULL);
   assert(problem != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &solv) );
   solv->problem = problem;
   solv->verbosity = verbosity;
   solv->weightupdatedecay = weightupdatedecay;
   solv->bestobjective = bestobjective;
   solv->checkiterations = checkiterations;
   solv->samplesize = samplesize;
   solv->randprobgood = randprobgood;
   solv->randprobviol = randprobviol;
   solv->callbackeffort = callbackeffort;
   solv->weightscalethreshold = weightscalethreshold;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->jumpmoves, problem->nvars) );
   for( i = 0; i < problem->nvars; ++i )
   {
      solv->jumpmoves[i].value = 0.0;
      solv->jumpmoves[i].score = -INFINITY;
   }

   solv->shiftbuffer = NULL;
   solv->shiftbuffersize = 0;
   solv->nshiftbuffer = 0;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->goodvarsset, problem->nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solv->goodvarssetidx, problem->nvars) );
   solv->ngoodvars = 0;
   solv->goodvarssize = problem->nvars;
   for( i = 0; i < problem->nvars; ++i )
      solv->goodvarssetidx[i] = -1;

   solv->randnumgen = randnumgen;

   solv->objectiveweight = 0.0;
   solv->bestviolationscore = INT_MAX;
   solv->effortatlastcallback = 0;
   solv->effortatlastimprovement = 0;
   solv->totaleffort = 0;
   solv->violationsdecrease = INT_MAX;
   solv->prevviolations = INT_MAX;
   solv->percentdecrease = 100;
   solv->iterations = 0;
   solv->weightupdateincrement = 1.0;
   solv->nbumps = 0;

   *solver = solv;

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

   assert(scip != NULL);
   assert(solver != NULL);

   solv = *solver;
   assert(solv != NULL);

   SCIPfreeBlockMemoryArray(scip, &solv->goodvarssetidx, solv->problem->nvars);
   SCIPfreeBlockMemoryArray(scip, &solv->goodvarsset, solv->goodvarssize);
   if( solv->shiftbuffer != NULL )
   {
      int i;
      for( i = 0; i < solv->shiftbuffersize; ++i )
      {
         SCIPfreeBlockMemory(scip, &solv->shiftbuffer[i]); /*lint !e866*/
      }
      SCIPfreeBlockMemoryArray(scip, &solv->shiftbuffer, solv->shiftbuffersize);
   }
   SCIPfreeBlockMemoryArray(scip, &solv->jumpmoves, solv->problem->nvars);
   SCIPfreeBlockMemory(scip, solver);

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

   SCIP_CALL( fjProblemResetIncumbent(scip, solver->problem, initialvalues) );

   solver->totaleffort += solver->problem->nnonzeros;

   /* reset variable scores */
   solver->ngoodvars = 0;
   for( i = 0; i < solver->problem->nvars; ++i )
   {
      SCIP_CALL( fjSolverResetMoves(scip, solver, i) );
   }

   return SCIP_OKAY;
}


/*
 * Local methods for heuristic
 */

/** adds a row to the feasibility jump solver */
static
SCIP_RETCODE addRowInFeasjumpSolver(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_PROBLEM*           problem,            /**< problem */
   SCIP_Real*            consvals,           /**< array with coefficients of the row */
   int*                  consinds,           /**< array with column indices of the row */
   int                   consnvars,          /**< number of variables in the row */
   SCIP_Real             lhs,                /**< left hand side of the row */
   SCIP_Real             rhs                 /**< right hand side of the row */
   )
{
   FJ_ROWTYPE rowtype;
   SCIP_Real solverrhs;
   int idx;

   assert(scip != NULL);
   assert(problem != NULL);

   /* ranged row */
   if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) && !SCIPisEQ(scip, lhs, rhs) )
   {
      solverrhs = lhs;
      SCIP_CALL( fjProblemAddConstraint(scip, problem, FJ_GTE, solverrhs, consnvars, consinds, consvals, FALSE, &idx) );
      solverrhs = rhs;
      SCIP_CALL( fjProblemAddConstraint(scip, problem, FJ_LTE, solverrhs, consnvars, consinds, consvals, FALSE, &idx) );
   }
   else
   {
      if( SCIPisEQ(scip, rhs, lhs) )
      {
         rowtype = FJ_EQUAL;
         solverrhs = rhs;
      }
      else if( SCIPisInfinity(scip, rhs) )
      {
         rowtype = FJ_GTE;
         solverrhs = lhs;
      }
      else
      {
         assert(SCIPisInfinity(scip, -lhs));
         rowtype = FJ_LTE;
         solverrhs = rhs;
      }
      SCIP_CALL( fjProblemAddConstraint(scip, problem, rowtype, solverrhs, consnvars, consinds, consvals,
            FALSE, &idx) );
   }

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
   SCIP_VAR** consvars;
   SCIP_CONSHDLR* conshdlrsos1;
   SCIP_CONSHDLR* conshdlrsos2;
   SCIP_CONSHDLR* conshdlrnonlinear;
   SCIP_Real* consvals;
   int* consinds;
   int nconss;
   int nvars;
   int i;
   int j;

   assert(scip != NULL);
   assert(problem != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* we do not use the heuristic before presolving if sos or non-linear constraints are present */
   conshdlrsos1 = SCIPfindConshdlr(scip, "SOS1");
   conshdlrsos2 = SCIPfindConshdlr(scip, "SOS2");
   conshdlrnonlinear = SCIPfindConshdlr(scip, "nonlinear");

   if( (conshdlrsos1 != NULL && SCIPconshdlrGetNConss(conshdlrsos1) > 0) ||
       (conshdlrsos2 != NULL && SCIPconshdlrGetNConss(conshdlrsos2) > 0) ||
       (conshdlrnonlinear != NULL && SCIPconshdlrGetNConss(conshdlrnonlinear) > 0) )
      return SCIP_OKAY;

   nconss = SCIPgetNConss(scip);
   nvars = SCIPgetNVars(scip);

   SCIP_CONS** conss = SCIPgetConss(scip);
   SCIP_VAR** probvars = SCIPgetVars(scip);

   consinds = NULL;
   consvals = NULL;
   consvars = NULL;
   SCIP_CALL( SCIPallocBufferArray(scip, &consinds, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nvars) );

   /* create feasjump solver variables */
   for( j = 0; j < nvars; ++j )
   {
      FJ_VARTYPE vartype;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real obj;
      int idx;

      lb = SCIPvarGetLbLocal(probvars[j]);
      ub = SCIPvarGetUbLocal(probvars[j]);
      obj = SCIPvarGetObj(probvars[j]);
      vartype = SCIPvarIsIntegral(probvars[j]) ? FJ_INTEGER : FJ_CONTINUOUS;
      SCIP_CALL( fjProblemAddVar(scip, problem, vartype, lb, ub, obj, &idx) );
   }

   /* create feasjump solver constraints */
   for( i = 0; i < nconss; ++i )
   {
      SCIP_CONS* cons;
      SCIP_Real lhs;
      SCIP_Real rhs;
      int consnvars;

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

      SCIP_CALL( SCIPgetConsNVars(scip, cons, &consnvars, success) );
      if( !(*success) )
         goto TERMINATE;

      /* check for non-active vars */
      for( j = 0; j < consnvars; ++j )
      {
         if( SCIPvarGetProbindex(consvars[j]) == -1 )
         {
            *success = FALSE;
            goto TERMINATE;
         }
         consinds[j] = SCIPvarGetProbindex(consvars[j]);
      }
      SCIP_CALL( addRowInFeasjumpSolver(scip, problem, consvals, consinds, consnvars, lhs, rhs) );
   }

   *success = TRUE;

TERMINATE:
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
   int                   nrows               /**< number of LP rows */
   )
{
   SCIP_Real* vals;
   int* inds;
   int i;
   int j;

   assert(scip != NULL);
   assert(problem != NULL);

   inds = NULL;
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, ncols) );

   /* add the columns */
   for( i = 0; i < ncols; ++i )
   {
      SCIP_COL* col;
      SCIP_VAR* var;
      FJ_VARTYPE vartype;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real obj;
      int idx;

      col = cols[i];
      var = SCIPcolGetVar(col);
      lb = SCIPcolGetLb(col);
      ub = SCIPcolGetUb(col);
      obj = SCIPcolGetObj(col);

      vartype = SCIPvarIsIntegral(var) ? FJ_INTEGER : FJ_CONTINUOUS;
      SCIP_CALL( fjProblemAddVar(scip, problem, vartype, lb, ub, obj, &idx) );
   }

   /* add the rows */
   for( i = 0; i < nrows; ++i )
   {
      SCIP_ROW* row;
      SCIP_COL** rowcols;
      SCIP_Real* rowvals;
      SCIP_Real lhs;
      SCIP_Real rhs;
      int nnonz;
      int nvals;

      row = rows[i];
      lhs = SCIProwGetLhs(row) - SCIProwGetConstant(row);
      rhs = SCIProwGetRhs(row) - SCIProwGetConstant(row);
      rowcols = SCIProwGetCols(row);
      rowvals = SCIProwGetVals(row);
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
            assert(vals[nvals] != 0.0);
            ++nvals;
         }
      }
      SCIP_CALL( addRowInFeasjumpSolver(scip, problem, vals, inds, nvals, lhs, rhs) );
   }

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &inds);

   return SCIP_OKAY;
}

/** adds objective cutoff constraint */
static
SCIP_RETCODE addObjCutoff(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_PROBLEM*           problem,            /**< problem */
   SCIP_COL**            cols,               /**< LP columns */
   int                   ncols               /**< number of LP columns */
   )
{
   SCIP_Real* vals;
   int* inds;
   SCIP_Real rhs;
   int nnonz;
   int objnnzs;
   int i;
   int idx;

   assert(scip != NULL);
   assert(problem != NULL);

   nnonz = ncols;

   inds = NULL;
   vals = NULL;
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, nnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nnonz) );

   objnnzs = 0;
   for( i = 0; i < ncols; ++i )
   {
      SCIP_COL* col;
      SCIP_Real obj;

      col = cols[i];
      obj = SCIPcolGetObj(col);

      if( SCIPisZero(scip, obj) )
         continue;

      inds[objnnzs] = SCIPcolGetLPPos(col);
      assert(inds[objnnzs] < ncols);

      vals[objnnzs] = obj;
      objnnzs++;
   }

   assert(!SCIPisInfinity(scip, SCIPgetCutoffbound(scip)));
   rhs = SCIPgetCutoffbound(scip) - SCIPcutoffbounddelta(scip);
   SCIP_CALL( fjProblemAddConstraint(scip, problem, FJ_LTE, rhs, objnnzs, inds, vals, TRUE, &idx) );

   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &vals);

   return SCIP_OKAY;
}

/** tries to add incumbent solution to SCIP */
static
SCIP_RETCODE checkIncumbentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic pointer */
   SCIP_HEURTIMING       heurtiming,         /**< heuristic timing mask */
   SCIP_Real*            incumbent,          /**< incumbent solution */
   SCIP_RESULT*          result              /**< pointer to store the result of the heuristic call */
   )
{
   SCIP_SOL* sol;
   SCIP_Bool stored;
   int i;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(incumbent != NULL);
   assert(result != NULL);

   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

   if( heurtiming == SCIP_HEURTIMING_BEFOREPRESOL )
   {
      SCIP_VAR** vars;
      int nvars;

      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);

      for( i = 0; i < nvars; ++i )
      {
         SCIP_CALL( SCIPsetSolVal(scip, sol, vars[i], incumbent[i]) );
      }
   }
   else
   {
      SCIP_COL** cols;
      int ncols;

      SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

      for( i = 0; i < ncols; ++i )
      {
         SCIP_VAR* var;

         var = SCIPcolGetVar(cols[i]);
         SCIP_CALL( SCIPsetSolVal(scip, sol, var, incumbent[i]) );
      }
   }

   SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, FALSE, TRUE, &stored) );

   if( stored )
   {
      SCIPdebugMsg(scip, "Feasjump found feasible solution\n");
      *result = SCIP_FOUNDSOL;
   }

   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   return SCIP_OKAY;
}

/** checks termination conditions and performs callback checks */
static
SCIP_RETCODE fjCheckTermination(
   SCIP*                 scip,               /**< SCIP data structure */
   FJ_SOLVER*            solver,             /**< solver */
   FJ_PROBLEM*           problem,            /**< problem */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_HEUR*            heur,               /**< heuristic */
   SCIP_HEURTIMING       heurtiming,         /**< heuristic timing */
   SCIP_Real*            solution,           /**< solution to report, or NULL */
   int*                  nsols,              /**< pointer to number of solutions found so far */
   SCIP_RESULT*          result,             /**< pointer to store result */
   SCIP_Real             starttime,          /**< start time */
   SCIP_Bool*            terminate           /**< pointer to store whether to terminate */
   )
{
   assert(scip != NULL);
   assert(solver != NULL);
   assert(problem != NULL);
   assert(heurdata != NULL);
   assert(heur != NULL);
   assert(nsols != NULL);
   assert(result != NULL);
   assert(terminate != NULL);

   *terminate = FALSE;

   /* check if callback should be performed */
   if( solution != NULL || solver->totaleffort - solver->effortatlastcallback > solver->callbackeffort )
   {
      SCIP_Bool quitnumsol;
      SCIP_Bool quiteffort;
      SCIP_Bool quitnoimprove;
      int nfoundsols = *nsols;

      if( heurdata->verbosity >= 2 )
         SCIPinfoMessage(scip, NULL, "Feasibility Jump: calling user termination.\n");

      solver->effortatlastcallback = solver->totaleffort;

      /* if we received a solution, check and report it */
      if( solution != NULL )
      {
         SCIP_CALL( checkIncumbentSol(scip, heur, heurtiming, solution, result) );
         if( *result == SCIP_FOUNDSOL )
            nfoundsols++;
         /* if the heuristic returns a solution that is not feasible we have to stop */
         else
            *terminate = TRUE;
      }

      /* check termination conditions */
      if( !(*terminate) )
      {
         quitnumsol = (nfoundsols >= heurdata->maxsols);
         if( quitnumsol )
         {
            SCIPdebugMsg(scip, "Feasibility Jump: quitting because number of solutions %d >= %d.\n",
               nfoundsols, heurdata->maxsols);
         }

         quiteffort = (solver->totaleffort - solver->effortatlastimprovement > heurdata->maxeffort);
         if( quiteffort )
         {
            SCIPdebugMsg(scip, "Feasibility Jump: quitting because effort %d > %d.\n",
               solver->totaleffort - solver->effortatlastimprovement, heurdata->maxeffort);
         }

         quitnoimprove = (solver->percentdecrease < heurdata->mindecrease);
         if( quitnoimprove )
         {
            SCIPdebugMsg(scip, "Feasibility Jump: Percentage decrease of violated constraints (%d) "
               "is smaller than the minimum of %d percent decrease.\n", solver->percentdecrease,
               heurdata->mindecrease);
         }

         *terminate = quitnumsol || quiteffort || quitnoimprove;
      }

      if( *terminate )
      {
         SCIP_Real time = SCIPgetSolvingTime(scip) - starttime;
         SCIPdebugMsg(scip, "Feasibility Jump: effort rate: %g Mops/sec\n",
            solver->totaleffort / time / 1.0e6);
         if( heurdata->verbosity >= 2 )
            SCIPinfoMessage(scip, NULL, "Feasibility Jump: quitting.\n");
      }
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
   int ncols;
   int nrows;
   int nsols;
   int step;
   SCIP_Real starttime;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( fjProblemCreate(scip, &problem) );

   if( heurtiming == SCIP_HEURTIMING_BEFOREPRESOL )
   {
      SCIP_Bool success;

      SCIP_CALL( extractProblemDataBeforePresolve(scip, problem, &success) );
      if( !success )
      {
         SCIP_CALL( fjProblemFree(scip, &problem) );
         return SCIP_OKAY;
      }
   }
   else
   {
      SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
      SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

      SCIP_CALL( extractProblemData(scip, problem, cols, rows, ncols, nrows) );

      /* add objective cutoff */
      if( SCIPgetBestSol(scip) != NULL )
      {
         SCIP_CALL( addObjCutoff(scip, problem, cols, ncols) );
      }
   }

   SCIP_CALL( fjSolverCreate(scip, &solver, problem, heurdata->randnumgen, heurdata->verbosity,
         heurdata->weightupdatedecay, SCIPinfinity(scip), heurdata->iterations, heurdata->samplesize,
         heurdata->randprobgood, heurdata->randprobviol, heurdata->callbackeffort,
         heurdata->weightscalethreshold) );

   starttime = SCIPgetSolvingTime(scip);
   nsols = 0;

   /* initialize with best solution if requested */
   if( SCIPgetBestSol(scip) != NULL && heurdata->useinitialsol && heurtiming != SCIP_HEURTIMING_BEFOREPRESOL )
   {
      SCIP_Real* initial = NULL;
      SCIP_SOL* bestsol = SCIPgetBestSol(scip);
      int i;

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

   if( heurdata->verbosity >= 1 )
      SCIPinfoMessage(scip, NULL, "Feasibility Jump: starting solve. weightUpdateDecay=%g, relaxContinuous=%d  \n",
         heurdata->weightupdatedecay, problem->usedrelaxcontinuous);

   /* main loop */
   for( step = 0; step < INT_MAX; ++step )
   {
      int varidx;
      SCIP_Bool terminate = FALSE;

      /* check termination at start of iteration */
      SCIP_CALL( fjCheckTermination(scip, solver, problem, heurdata, heur, heurtiming, NULL, &nsols, result,
            starttime, &terminate) );
      if( terminate )
         break;

      /* periodic output */
      if( step % 100000 == 0 && heurdata->verbosity >= 1 )
      {
         SCIPinfoMessage(scip, NULL, "Feasibility Jump: step %d viol %d good %d bumps %d\n",
            step, problem->nviolated, solver->ngoodvars, solver->nbumps);
      }

      /* check for improvement */
      if( problem->nviolated < solver->bestviolationscore )
      {
         solver->effortatlastimprovement = solver->totaleffort;
         solver->bestviolationscore = problem->nviolated;
      }

      if( problem->nviolated == 0 && problem->incumbentobjective < solver->bestobjective )
      {
         solver->effortatlastimprovement = solver->totaleffort;
         solver->bestobjective = problem->incumbentobjective;
         solver->percentdecrease = 100;

         /* check termination after finding solution */
         SCIP_CALL( fjCheckTermination(scip, solver, problem, heurdata, heur, heurtiming,
               problem->incumbentassignment, &nsols, result, starttime, &terminate) );
         if( terminate || *result != SCIP_FOUNDSOL )
            break;
      }

      if( SCIPisStopped(scip) || problem->nvars == 0 )
         break;

      /* select and perform move */
      SCIP_CALL( fjSolverSelectVariable(scip, solver, &varidx) );
      SCIP_CALL( fjSolverDoVariableMove(scip, solver, varidx) );

      /* check progress every N iterations */
      if( step % heurdata->iterations == 0 && step <= heurdata->iterations )
      {
         solver->prevviolations = problem->nviolated;
      }
      else if( step % heurdata->iterations == 0 )
      {
         if( heurdata->verbosity >= 1 )
            SCIPinfoMessage(scip, NULL, "Feasibility Jump: viol %d prev viol %d \n",
               problem->nviolated, solver->prevviolations);

         solver->violationsdecrease = solver->prevviolations - problem->nviolated;
         if( solver->prevviolations == 0 )
            solver->percentdecrease = 100;
         else
            solver->percentdecrease = (100 * solver->violationsdecrease) / solver->prevviolations;

         solver->prevviolations = problem->nviolated;

         /* check termination after checking progress */
         SCIP_CALL( fjCheckTermination(scip, solver, problem, heurdata, heur, heurtiming, NULL, &nsols, result,
               starttime, &terminate) );
         if( terminate )
            break;
      }
   }

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

   /* call inclusion method of primal heuristic */
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

   /* free heuristic data */
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

   /* create random number generator */
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

   /* free random number generator */
   SCIPfreeRandom(scip, &heurdata->randnumgen);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecFeasjump)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_Bool cutoff;

   heurdata = SCIPheurGetData(heur);

   assert(result != NULL);
   assert(heurdata != NULL);

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMsg(scip, "Executing Feasjump\n");

   if( !heurdata->beforepresol && heurtiming == SCIP_HEURTIMING_BEFOREPRESOL )
      return SCIP_OKAY;

   if( heurtiming == SCIP_HEURTIMING_BEFORENODE )
   {
      /* stop execution method if we call the heuristic only when no primal feasible solution is found yet */
      if( SCIPgetBestSol(scip) != NULL && heurdata->onlywithoutsol )
         return SCIP_OKAY;

      /* stop if there is no LP available */
      if( !SCIPhasCurrentNodeLP(scip) )
         return SCIP_OKAY;

      if( !SCIPisLPConstructed(scip) )
      {
         SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

         /* manually cut off the node if the LP construction detected infeasibility */
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

   /* create Feasjump primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
   assert(heurdata != NULL);

   heurdata->onlywithoutsol = DEFAULT_ONLYWITHOUTSOL;
   heurdata->weightupdatedecay = DEFAULT_WEIGHTDECAY;
   heurdata->maxeffort = DEFAULT_MAXEFFORT;
   heurdata->verbosity = DEFAULT_VERBOSITY;
   heurdata->useinitialsol = DEFAULT_USEINITIALSOL;
   heurdata->beforepresol = DEFAULT_BEFORE_PRESOL;
   heurdata->iterations = DEFAULT_ITERATIONS;
   heurdata->mindecrease = DEFAULT_MINPERCENTDECREASE;
   heurdata->maxsols = DEFAULT_MAXSOLS;
   heurdata->samplesize = DEFAULT_SAMPLESIZE;
   heurdata->randprobgood = DEFAULT_RANDPROBGOOD;
   heurdata->randprobviol = DEFAULT_RANDPROBVIOL;
   heurdata->callbackeffort = DEFAULT_CALLBACKEFFORT;
   heurdata->weightscalethreshold = DEFAULT_WEIGHTSCALETHRESHOLD;
   heurdata->randnumgen = NULL;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecFeasjump, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyFeasjump) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeFeasjump) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitFeasjump) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitFeasjump) );

   /* add feasjump primal heuristic parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/weightupdatedecay",
         "decay factor for the weight update",
         &heurdata->weightupdatedecay, FALSE, DEFAULT_WEIGHTDECAY, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/maxeffort",
         "maximum effort spend in feasibility jump",
         &heurdata->maxeffort, FALSE, DEFAULT_MAXEFFORT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/verbosity",
         "verbosity level of the feasibility jump solver",
         &heurdata->verbosity, FALSE, DEFAULT_VERBOSITY, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/onlywithoutsol",
         "should the feasibility jump be called only if no solution exists?",
         &heurdata->onlywithoutsol, FALSE, DEFAULT_ONLYWITHOUTSOL, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/beforepresol",
         "should the heuristic be called before presolving?",
         &heurdata->beforepresol, FALSE, DEFAULT_BEFORE_PRESOL, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/useinitialsol",
         "should the feasibility jump use the best known solution as initialization?",
         &heurdata->useinitialsol, FALSE, DEFAULT_USEINITIALSOL, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/iterations",
         "number of iterations to check if the heuristic reduces the total number of violated constraints",
         &heurdata->iterations, FALSE, DEFAULT_ITERATIONS, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/mindecrease",
         "minimum percentage decrease for the number of violated constraints",
         &heurdata->mindecrease, FALSE, DEFAULT_MINPERCENTDECREASE, 0, 100, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/maxsols",
         "maximum number of solutions to find",
         &heurdata->maxsols, FALSE, DEFAULT_MAXSOLS, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/samplesize",
         "number of candidates to sample in tournament selection",
         &heurdata->samplesize, FALSE, DEFAULT_SAMPLESIZE, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/randprobgood",
         "probability of random selection from good variables",
         &heurdata->randprobgood, FALSE, DEFAULT_RANDPROBGOOD, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/randprobviol",
         "probability of random selection from violated constraint",
         &heurdata->randprobviol, FALSE, DEFAULT_RANDPROBVIOL, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/callbackeffort",
         "effort between callbacks",
         &heurdata->callbackeffort, FALSE, DEFAULT_CALLBACKEFFORT, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/weightscalethreshold",
         "threshold for rescaling weights to avoid overflow",
         &heurdata->weightscalethreshold, FALSE, DEFAULT_WEIGHTSCALETHRESHOLD, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

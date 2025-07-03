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

/* NLP interface for the CONOPT solver */

/**@file    nlpi_conopt.c
 * @ingroup DEFPLUGINS_NLPI
 * @brief   CONOPT NLP interface
 * @author  Ksenia Bestuzheva
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/nlpi_conopt.h"
#include "scip/nlpioracle.h"
#include "scip/scip_mem.h"
#include "scip/scip_numerics.h"
#include "scip/scip_nlp.h"
#include "scip/scip_nlpi.h"
#include "scip/scip_randnumgen.h"
#include "scip/pub_misc.h"
#include "scip/pub_message.h"
#include "scip/type_clock.h"

#include "coiheader.h"
#include "scip_message.h"
#include "scip_timing.h"

#define NLPI_NAME              "conopt"                    /**< short concise name of solver */
#define NLPI_DESC              "solver interface template" /**< description of solver */
#define NLPI_PRIORITY          0                           /**< priority of NLP solver */

#define DEFAULT_RANDSEED       107                         /**< initial random seed */
#define MAXPERTURB             0.01                        /**< maximal perturbation of bounds in starting point heuristic */

/*
 * Data structures
 */

struct SCIP_NlpiData
{
   SCIP_CLOCK*           solvetime;          /**< clock for measuring solving time */
};

struct SCIP_NlpiProblem
{
   SCIP*                 scip;               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle;             /**< Oracle-helper to store and evaluate NLP */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */

   SCIP_Bool             firstrun;           /**< whether the next NLP solve will be the first one (with the current problem structure) */
   SCIP_Real*            initguess;          /**< initial values for primal variables, or NULL if not known */

   SCIP_NLPSOLSTAT       solstat;            /**< solution status from last NLP solve */
   SCIP_NLPTERMSTAT      termstat;           /**< termination status from last NLP solve */
   SCIP_Real             solvetime;          /**< time spend for last NLP solve */
   int                   niterations;        /**< number of iterations for last NLP solve */
   SCIP_Real             objval;             /**< objective value from last run */

   SCIP_Real*            lastprimal;         /**< primal solution from last run, if available */
   SCIP_Real*            lastdualcons;       /**< dual solution from last run, if available */
   SCIP_Real*            lastduallb;         /**< dual solution for lower bounds from last run, if available */
   SCIP_Real*            lastdualub;         /**< dual solution for upper bounds from last run, if available */

   coiHandle_t           CntVect;            /**< pointer to CONOPT Control Vector */

   /* options that will be passed to CONOPT via callbacks */
   int                   verblevel;          /**< verbosity level, same as in SCIP_NlpParam: 0 off, 1 normal, 2 debug, > 2 more debug */
   SCIP_Real             opttol;
};


/*
 * Local methods
 */

/** Implementations of CONOPT callbacks */

/** CONOPT callback to pass solution back to SCIP
 *
 * The dual values are multiplied by -1, since CONOPT's definition of marginal values is 'by
 * how much the objective increases if we increase the constraint side by +eps', which, for
 * minimisation problems (which is always the case here), is the negative of what SCIP expects.
 */
static int COI_CALLCONV Solution(
   const double          XVAL[],             /**< solution values of the variables (provided by CONOPT) */
   const double          XMAR[],             /**< marginal values (provided by CONOPT) */
   const int             XBAS[],             /**< basis indicators for the variables (provided by CONOPT) */
   const int             XSTA[],             /**< status values for the variables (provided by CONOPT) */
   const double          YVAL[],             /**< values of the left hand sides of all rows in the optimal solution (provided by CONOPT) */
   const double          YMAR[],             /**< marginal values corresponding to rows (provided by CONOPT) */
   const int             YBAS[],             /**< basis indicators for the rows or constraints (provided by CONOPT) */
   const int             YSTA[],             /**< status values for the rows or constraints (provided by CONOPT) */
   int                   NUMVAR,             /**< number of variables (provided by CONOPT) */
   int                   NUMCON,             /**< number of constraints (provided by CONOPT) */
   void*                 USRMEM              /**< user memory pointer (i.e. pointer to SCIP_NLPIPROBLEM) */
   )
{
   SCIP_NLPIPROBLEM* problem = (SCIP_NLPIPROBLEM*)USRMEM;
   int noraclevars;
   int noracleconss;

   assert(problem != NULL);
   assert(problem->oracle != NULL);

   noraclevars = SCIPnlpiOracleGetNVars(problem->oracle);
   noracleconss = SCIPnlpiOracleGetNConstraints(problem->oracle);

   /* number of SCIP variables always less or equal, since CONOPT variables can also contain slack variables */
   assert(NUMVAR >= noraclevars);

   /* CONOPT has one more constraint than oracle (the objective) */
   assert(NUMCON == noracleconss + 1);

   /* copy values from CONOPT into SCIP arrays. Note that in CONOPT, there is one
    * extra constraint (the objective) and slack variables that are not explicitly
    * present in SCIP, which the code below does not copy */
   if( problem->lastprimal == NULL )
   {
      if( NUMVAR > 0 )
      {
         if( SCIPduplicateBlockMemoryArray(problem->scip, &problem->lastprimal, XVAL, noraclevars) != SCIP_OKAY ||
               SCIPallocClearBlockMemoryArray(problem->scip, &problem->lastduallb, noraclevars) != SCIP_OKAY ||
               SCIPallocClearBlockMemoryArray(problem->scip, &problem->lastdualub, noraclevars) != SCIP_OKAY )
         {
            SCIPerrorMessage("Failed to allocate memory for a solution from CONOPT\n");
            return SCIP_NOMEMORY;
         }
      }
      if( NUMCON > 0 )
      {
         if( SCIPduplicateBlockMemoryArray(problem->scip, &problem->lastdualcons, YMAR, noracleconss) != SCIP_OKAY )
         {
            SCIPerrorMessage("Failed to allocate memory for a solution from CONOPT\n");
            return SCIP_NOMEMORY;
         }
      }
   }
   else
   {
      BMScopyMemoryArray(problem->lastprimal, XVAL, noraclevars);
      BMSclearMemoryArray(problem->lastduallb, noraclevars);
      BMSclearMemoryArray(problem->lastdualub, noraclevars);
      BMScopyMemoryArray(problem->lastdualcons, YMAR, noracleconss);
   }

   /* replace initial guess with the obtained solution */
   if( problem->lastprimal != NULL )
   {
      if( problem->initguess == NULL )
      {
         if( SCIPduplicateMemoryArray(problem->scip, &problem->initguess, problem->lastprimal, noraclevars) != SCIP_OKAY )
         {
            SCIPerrorMessage("Failed to allocate memory for an initial guess from\n");
            return SCIP_NOMEMORY;
         }
      }
      else
         BMScopyMemoryArray(problem->initguess, problem->lastprimal, noraclevars);
   }

   /* get dual multipliers for variable bounds */
   for( int i = 0; i < noraclevars; i++ )
   {
      if( XBAS[i] == 0 ) /* Xi is at lower bound */
         problem->lastduallb[i] = -XMAR[i];
      else if( XBAS[i] == 1 ) /* Xi is at upper bound */
         problem->lastdualub[i] = -XMAR[i];
   }

   for( int i = 0; i < NUMCON-1; i++ )
      (problem->lastdualcons[i]) *= -1;

   return 0;
}

/** CONOPT callback to pass variable bounds, constraint types and sides
 * and Jacobian structure and constant values to CONOPT
 */
static int COI_CALLCONV ReadMatrix(
   double                LOWER[],            /**< lower bounds on the variables (set to CONOPT's minus infinity before callback is issued */
   double                CURR[],             /**< initial values of the variables (set to zero before callback is issued */
   double                UPPER[],            /**< upper bounds on the variables (set to CONOPT's plus infinity before callback is issued */
   int                   VSTA[],             /**< initial status values for the variable (used if coidef_inistat() was called with IniStat = 1 or 2) */
   int                   TYPE[],             /**< constraint types (equation, inequality, free) */
   double                RHS[],              /**< right hand sides values of constraints (default is zero) */
   int                   ESTA[],             /**<initial status values for constraint slacks (used if coidef_inistat() was called with IniStat = 1 or 2) */
   int                   COLSTA[],           /**< starting indices of Jacobian columns in ROWNO */
   int                   ROWNO[],            /**< row numbers of Jacobian nonzeros */
   double                VALUE[],            /**< values of the Jacobian elements (defined for all constant Jacobian elements) */
   int                   NLFLAG[],           /**< flags marking nonlinear Jacobian elements (following the same sparse representation as ROWNO) */
   int                   NUMVAR,             /**< number of variables as defined in COIDEF_NumVar() */
   int                   NUMCON,             /**< number of constraints as defined in COIDEF_NumCon() */
   int                   NUMNZ,              /**< number of Jacobian elements as defined in COIDEF_NumNz() */
   void*                 USRMEM              /**< user memory pointer (i.e. pointer to SCIP_NLPIPROBLEM) */
   )
{
   SCIP_NLPIPROBLEM* problem = (SCIP_NLPIPROBLEM*)USRMEM;
   SCIP* scip;
   SCIP_NLPIORACLE* oracle;
   const SCIP_Real* lbs;
   const SCIP_Real* ubs;
   int norigvars;
   int nslackvars = 0;
   int njacnlnnz;
   int* rangeconsidxs = NULL;
   const SCIP_Bool* jacrownlflags;
   const int* jaccoloffsets;
   const int* jacrows;
   int nobjnz;
   const int* objnz;
   const SCIP_Bool* objnlflags;
   int nobjnlnz;
   int objnzcnt = 0;
   int* nrownz;
   SCIP_RETCODE retcode;

   assert(problem != NULL);

   oracle = problem->oracle;
   assert(oracle != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   SCIPdebugMsg(scip, "NLP sizes passed to CONOPT: NUMVAR = %d, NUMCON = %d, NUMNZ = %d\n", NUMVAR, NUMCON, NUMNZ);

   norigvars = SCIPnlpiOracleGetNVars(oracle);
   lbs = SCIPnlpiOracleGetVarLbs(oracle);
   ubs = SCIPnlpiOracleGetVarUbs(oracle);

   /* save indices of range constraints if there are any */
   if( NUMVAR - norigvars > 0 )
   {
      retcode = SCIPallocBufferArray(scip, &rangeconsidxs, NUMVAR - norigvars);
      if( retcode != SCIP_OKAY )
      {
         SCIPerrorMessage("No memory in a callback of CONOPT\n");
         return retcode;
      }
   }

   /* add all 'normal' (i.e. non-slack) variables here */
   for( int i = 0; i < norigvars; i++ )
   {
      if( !SCIPisInfinity(scip, -lbs[i]) )
         LOWER[i] = lbs[i];
      if( !SCIPisInfinity(scip, ubs[i]) )
         UPPER[i] = ubs[i];
   }

   /* specify initial values of original variables,
    * (if there is a previous solution, it is stored in initguess) */
   if( problem->initguess != NULL )
   {
      BMScopyMemoryArray(CURR, problem->initguess, norigvars);

      /* some values may have been set outside the bounds - project them */
      for( int i = 0; i < norigvars; ++i )
      {
         if( lbs[i] > CURR[i] )
            CURR[i] = SCIPrandomGetReal(problem->randnumgen, lbs[i], lbs[i] + MAXPERTURB*MIN(1.0, ubs[i]-lbs[i]));
         else if( ubs[i] < CURR[i] )
            CURR[i] = SCIPrandomGetReal(problem->randnumgen, ubs[i] - MAXPERTURB*MIN(1.0, ubs[i]-lbs[i]), ubs[i]);
      }
   }
   else
   {
      /* if no initial guess given, project 0 onto variable bounds */
      assert(problem->randnumgen != NULL);

      SCIPdebugMsg(scip, "CONOPT started without initial primal values; make up starting guess by projecting 0 onto variable bounds\n");

      for( int i = 0; i < norigvars; ++i )
      {
         if( lbs[i] > 0.0 )
            CURR[i] = SCIPrandomGetReal(problem->randnumgen, lbs[i], lbs[i] + MAXPERTURB*MIN(1.0, ubs[i]-lbs[i]));
         else if( ubs[i] < 0.0 )
            CURR[i] = SCIPrandomGetReal(problem->randnumgen, ubs[i] - MAXPERTURB*MIN(1.0, ubs[i]-lbs[i]), ubs[i]);
         else
            CURR[i] = SCIPrandomGetReal(problem->randnumgen,
               MAX(lbs[i], -MAXPERTURB*MIN(1.0, ubs[i]-lbs[i])), MIN(ubs[i], MAXPERTURB*MIN(1.0, ubs[i]-lbs[i])));
      }
   }

   for( int i = 0; i < NUMCON-1; i++ )
   {
      SCIP_Real lhs = SCIPnlpiOracleGetConstraintLhs(oracle, i);
      SCIP_Real rhs = SCIPnlpiOracleGetConstraintRhs(oracle, i);

      assert(lhs <= rhs);

      if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) )
      {
         TYPE[i] = 0; /* an equality or a ranged row modelled as equality */

         if( !SCIPisEQ(scip, lhs, rhs) )
         {
            assert(rangeconsidxs != NULL);

            /* range constraint lhs <= g(x) <= rhs: reformulate as g(x) - s = 0 and lhs <= s <= rhs */
            RHS[i] = 0.0;
            LOWER[norigvars + nslackvars] = lhs;
            UPPER[norigvars + nslackvars] = rhs;

            /* set initial value of slack variable */
            if( lhs > 0.0 )
               CURR[norigvars + nslackvars] = SCIPrandomGetReal(problem->randnumgen, lhs, lhs + MAXPERTURB*MIN(1.0, rhs-lhs));
            else if( rhs < 0.0 )
               CURR[norigvars + nslackvars] = SCIPrandomGetReal(problem->randnumgen, rhs - MAXPERTURB*MIN(1.0, rhs-lhs), rhs);
            else
               CURR[norigvars + nslackvars] = SCIPrandomGetReal(problem->randnumgen,
                  MAX(lhs, -MAXPERTURB*MIN(1.0, rhs-lhs)), MIN(rhs, MAXPERTURB*MIN(1.0, rhs-lhs)));

            rangeconsidxs[nslackvars] = i;
            nslackvars++;
         }
         else
            RHS[i] = lhs;
      }
      else if( !SCIPisInfinity(scip, -lhs) )
      {
         TYPE[i] = 1;
         RHS[i] = lhs;
      }
      else if( !SCIPisInfinity(scip, rhs) )
      {
         TYPE[i] = 2;
         RHS[i] = rhs;
      }
      else
         TYPE[i] = 3;
   }
   assert(norigvars + nslackvars == NUMVAR);

   /* the last constraint is the objective */
   TYPE[NUMCON-1] = 3; /* objective must be a free row */
   RHS[NUMCON-1] = 0.0;

   /* Jacobian information */

   /* TODO make it so that the column representation is only composed when asked for */
   SCIP_CALL( SCIPnlpiOracleGetJacobianColSparsity(scip, oracle, &jaccoloffsets, &jacrows, &jacrownlflags, &njacnlnnz) );
   assert(jaccoloffsets == NULL || jaccoloffsets[norigvars] <= NUMNZ);

   /* move structure info into COLSTA and ROWNO; while doing so, also add nonzeroes for the objective
    * (which CONOPT sees as the last constraint, i.e. constraint with index NUMCON-1) */
   retcode = SCIPallocCleanBufferArray(scip, &nrownz, NUMCON);
   if( retcode != SCIP_OKAY )
   {
      SCIPerrorMessage("No memory in a callback of CONOPT\n");
      return retcode;
   }

   retcode = SCIPnlpiOracleGetObjGradientNnz(scip, oracle, &objnz, &objnlflags, &nobjnz, &nobjnlnz);
   if( retcode != SCIP_OKAY )
   {
      SCIPerrorMessage("Error in the ReadMatrix callback of CONOPT\n");
      return retcode;
   }

   for( int i = 0; i < norigvars; i++ )
   {
      /* starts of columns get shifted by how many objective nonzeros were added */
      COLSTA[i] = jaccoloffsets != NULL ? jaccoloffsets[i] + objnzcnt : objnzcnt;

      if( jaccoloffsets != NULL )
      {
         /* nonzeroes of constraints */
         for( int j = jaccoloffsets[i]; j < jaccoloffsets[i+1]; j++ )
         {
            ROWNO[j+objnzcnt] = jacrows[j];
            NLFLAG[j+objnzcnt] = jacrownlflags[j] ? 1 : 0;
            if( NLFLAG[j+objnzcnt] == 0 )
            {
               VALUE[j+objnzcnt] = SCIPnlpiOracleGetConstraintCoef(oracle, jacrows[j], nrownz[jacrows[j]]);
               ++(nrownz[jacrows[j]]);
            }
         }
      }

      /* nonzeroes of objective */
      if( i == objnz[objnzcnt] )
      {
         int idx = jaccoloffsets != NULL ? jaccoloffsets[i+1] + objnzcnt : objnzcnt;

         ROWNO[idx] = NUMCON - 1;
         NLFLAG[idx] = objnlflags[objnzcnt] ? 1 : 0;
         if( NLFLAG[idx] == 0 )
         {
            /* in the oracle, index -1 is used for the objective */
            VALUE[idx] = SCIPnlpiOracleGetConstraintCoef(oracle, -1, nrownz[NUMCON-1]);
            ++(nrownz[NUMCON-1]);
         }
         ++objnzcnt;
      }
   }
   assert(COLSTA[0] == 0);
   COLSTA[norigvars] = jaccoloffsets != NULL ? jaccoloffsets[norigvars] + objnzcnt : objnzcnt;
   BMSclearMemoryArray(nrownz, NUMCON);
   SCIPfreeCleanBufferArray(scip, &nrownz);

   if( nslackvars > 0 )
   {
      /* add a nonzero for each slack variable */
      for( int i = 0; i < nslackvars; i++ )
      {
         COLSTA[norigvars+i] = COLSTA[norigvars] + i; /* for each slack var, only one nonzero is added */
         ROWNO[COLSTA[norigvars+i]] = rangeconsidxs[i];
         NLFLAG[COLSTA[norigvars+i]] = 0;
         VALUE[COLSTA[norigvars+i]] = -1.0;
      }
      SCIPfreeBufferArray(scip, &rangeconsidxs);
      COLSTA[NUMVAR] = NUMNZ;
   }

#ifndef NDEBUG
   for( int i = 0; i < NUMNZ; ++i )
      assert(ROWNO[i] >= 0 && ROWNO[i] < NUMCON);
#endif

#ifdef STRUCTURE_DEBUG
   SCIPdebugMsg(scip, "Jacobian structure information:\n");
   SCIPdebugMsg(scip, "COLSTA =\n");
   for( int i = 0; i <= NUMVAR; i++ )
      printf("%d, ", COLSTA[i]);
   printf("\n");

   SCIPdebugMsg(scip, "ROWNO =\n");
   for( int i = 0; i < NUMNZ; i++ )
      printf("%d, ", ROWNO[i]);
   printf("\n");

   SCIPdebugMsg(scip, "NLFLAG =\n");
   for( int i = 0; i < NUMNZ; i++ )
      printf("%d, ", NLFLAG[i]);
   printf("\n");

   SCIPdebugMsg(scip, "VALUE =\n");
   for( int i = 0; i < NUMNZ; i++ )
      if( SCIPisInfinity(scip, VALUE[i]) )
         printf("inf, ");
      else
         printf("%g, ", VALUE[i]);
   printf("\n");
   fflush(stdout);
#endif

   return 0;
}

/** callback for CONOPT's standard output */
static int COI_CALLCONV Message(
   int                   SMSG,               /**< number of lines in the message that should go to the Screen file, between 0 and 30 */
   int                   DMSG,               /**< number of lines in the message that should go to the Status file, between 0 and 30 */
   int                   NMSG,               /**< number of lines in the message that should go to the Documentation file, between 0 and 30 */
   char*                 MSGV[],             /**< array with the lengths of the individual message lines, the lengths are between 1 and 132 */
   void*                 USRMEM              /**< user memory pointer (i.e. pointer to SCIP_NLPIPROBLEM) */
   )
{
   /* TODO how to handle status and documentation messages? */
   SCIP_NLPIPROBLEM* problem = (SCIP_NLPIPROBLEM*)USRMEM;

   assert(problem != NULL);
   assert(problem->scip != NULL);

   if( problem->verblevel > 0 )
   {
      for( int i = 0; i < SMSG; i++ )
         SCIPinfoMessage(problem->scip, NULL, "%s\n", MSGV[i]);
   }

   return 0;
}

/** callback for CONOPT's standard error output */
static int COI_CALLCONV ErrMsg(
   int                   ROWNO,              /**< number of row */
   int                   COLNO,              /**< number of column */
   int                   POSNO,              /**< number of a Jacobian element */
   const char*           MSG,                /**< the message */
   void*                 USRMEM              /**< user memory pointer (i.e. pointer to SCIP_NLPIPROBLEM) */
   )
{
   SCIP_NLPIPROBLEM* problem = (SCIP_NLPIPROBLEM*)USRMEM;

   if( problem->verblevel > 0 )
   {
      if( ROWNO == -1 && COLNO == -1 )
         SCIPinfoMessage(problem->scip, NULL, "\nCONOPT error/warning: Jacobian element %d: ", POSNO);
      else if( ROWNO == -1 )
         SCIPinfoMessage(problem->scip, NULL, "\nCONOPT error/warning: variable %d: ", COLNO);
      else if( COLNO == -1 )
         SCIPinfoMessage(problem->scip, NULL, "\nCONOPT error/warning: constraint %d: ", ROWNO);
      else
         SCIPinfoMessage(problem->scip, NULL, "\nCONOPT error/warning: variable %d appearing in constraint %d: ", COLNO, ROWNO);
      SCIPinfoMessage(problem->scip, NULL, "%s\n", MSG);
   }

   return 0;
}

/** callback for CONOPT to report the solving statuses */
static int COI_CALLCONV Status(
   int                   MODSTA,             /**< model status (corresponds to SCIP's solstat) */
   int                   SOLSTA,             /**< solver status (corresponds to SCIP's termstat) */
   int                   ITER,               /**< number of iterations */
   double                OBJVAL,             /**< objective value */
   void*                 USRMEM              /**< user memory pointer (i.e. pointer to SCIP_NLPIPROBLEM) */
   )
{
   SCIP* scip;
   SCIP_NLPIPROBLEM* problem = (SCIP_NLPIPROBLEM*)USRMEM;

   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   SCIPdebugMsg(scip, "CONOPT has finished optimizing\n");
   SCIPdebugMsg(scip, "Iteration count = %8d\n", ITER);
   SCIPdebugMsg(scip, "Objective value = %10f\n", OBJVAL);

   problem->niterations = ITER;
   problem->objval = OBJVAL;

   switch( MODSTA )
   {
      case 1:
         SCIPdebugMsg(scip, "NLP problem solved to global optimality\n");
         problem->solstat = SCIP_NLPSOLSTAT_GLOBOPT;
         break;
      case 2:
         SCIPdebugMsg(scip, "NLP problem solved to local optimality\n");
         problem->solstat = SCIP_NLPSOLSTAT_LOCOPT;
         break;
      case 3:
         SCIPdebugMsg(scip, "NLP problem unbounded\n");
         problem->solstat = SCIP_NLPSOLSTAT_UNBOUNDED;
         break;
      case 4:
         SCIPdebugMsg(scip, "NLP problem infeasible\n");
         problem->solstat = SCIP_NLPSOLSTAT_GLOBINFEASIBLE;
         break;
      case 5:
         SCIPdebugMsg(scip, "NLP problem locally infeasible\n");
         problem->solstat = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
         break;
      case 6: /* intermediate infeasible */
      case 7: /* intermediate non-optimal */
      case 12: /* unknown error */
      case 13:
         SCIPdebugMsg(scip, "NLP problem status unknown (CONOPT status %d)\n", MODSTA);
         problem->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
         break;
      default:
         SCIPerrorMessage("CONOPT returned an unexpected solution status %d\n", MODSTA);
         problem->solstat  = SCIP_NLPSOLSTAT_UNKNOWN;
   }

   switch( SOLSTA )
   {
      case 1:
         problem->termstat = SCIP_NLPTERMSTAT_OKAY;
         SCIPdebugMsg(scip, "CONOPT terminated with normal status.\n");
         break;
      case 2:
         problem->termstat = SCIP_NLPTERMSTAT_ITERLIMIT;
         SCIPdebugMsg(scip, "CONOPT terminated due to an iteration limit.\n");
         break;
      case 3:
         problem->termstat = SCIP_NLPTERMSTAT_TIMELIMIT;
         SCIPdebugMsg(scip, "CONOPT terminated due to a time limit.\n");
         break;
      case 5:
         problem->termstat = SCIP_NLPTERMSTAT_EVALERROR;
         SCIPdebugMsg(scip, "CONOPT terminated due to evaluation errors.\n");
         break;
      case 8:
         problem->termstat = SCIP_NLPTERMSTAT_INTERRUPT;
         SCIPdebugMsg(scip, "CONOPT interrupted by user.\n");
         break;
      case 4: /* terminated by solver */
      case 6: /* unknown */
      case 9: /* error: setup failure */
      case 10: /* error: solver failure */
      case 11: /* error: internal solver error */
         SCIPdebugMsg(scip, "CONOPT terminated with status %d\n", SOLSTA);
         problem->termstat = SCIP_NLPTERMSTAT_OTHER;
         break;
      default:
         SCIPerrorMessage("CONOPT returned an unexpected termination status %d\n", SOLSTA);
         problem->termstat = SCIP_NLPTERMSTAT_OTHER;
   }

   return 0;
}

/** CONOPT callback for function and Jacobian evaluation
 *
 *  The callback has three modes, indicated by MODE:
 *
 *  1: Only evaluate the sum of the nonlinear and linear terms in row ROWNO and return the value in G.
 *  2: Only evaluate the nonlinear Jacobian elements in row ROWNO and return them in JAC.
 *  3: Perform both option 1 and 2.
 */
static int COI_CALLCONV FDEval(
   const double          X[],                /**< point of evaluation (provided by CONOPT) */
   double*               G,                  /**< value of the function */
   double                JAC[],              /**< vector of Jacobian values */
   int                   ROWNO,              /**< number of the row for which nonlinearities are to be evaluated (provided by CONOPT) */
   const int             JACNUM[],           /**< list of column numbers for the nonlinear nonzero Jacobian elements in the current row (provided by CONOPT when MODE = 2 or 3) */
   int                   MODE,               /**< indicator for mode of evaluation (provided by CONOPT) */
   int                   IGNERR,             /**< indicator whether CONOPT assumes the point to be safe (0) or potentially unsafe (1) */
   int*                  ERRCNT,             /**< scalar function evaluation error indicator (set to 1 if a function value cannot be computed */
   int                   NUMVAR,             /**< number of variables (provided by CONOPT) */
   int                   NUMJAC,             /**< number of nonlinear nonzero Jacobian elements in the current row */
   int                   THREAD,             /**< only relevant in multi-threading environments */
   void*                 USRMEM              /**< user memory pointer (i.e. pointer to SCIP_NLPIPROBLEM) */
   )
{
   SCIP_RETCODE retcode;
   SCIP_NLPIPROBLEM* problem = (SCIP_NLPIPROBLEM*)USRMEM;

   assert(problem != NULL);
   assert(ROWNO <= SCIPnlpiOracleGetNConstraints(problem->oracle));

   if( MODE == 1 || MODE == 3 )
   {
      /* the last 'constraint' is the objective */
      retcode = ROWNO < SCIPnlpiOracleGetNConstraints(problem->oracle) ?
            SCIPnlpiOracleEvalConstraintValue(problem->scip, problem->oracle, ROWNO, X, G) :
            SCIPnlpiOracleEvalObjectiveValue(problem->scip, problem->oracle, X, G);
      if( retcode != SCIP_OKAY || *G == SCIP_INVALID )
         *ERRCNT = 1;
   }

   if( MODE == 2 || MODE == 3 )
   {
      SCIP_Real conval;

      /* the last 'constraint' is the objective */
      retcode = ROWNO < SCIPnlpiOracleGetNConstraints(problem->oracle) ?
            SCIPnlpiOracleEvalConstraintGradient(problem->scip, problem->oracle, ROWNO, X, TRUE, &conval, JAC) :
            SCIPnlpiOracleEvalObjectiveGradient(problem->scip, problem->oracle, X, TRUE, &conval, JAC);
      if( retcode != SCIP_OKAY )
         *ERRCNT = 1;
   }

   return 0;
}

/** CONOPT callback to pass some of the options (not supported via COIDEF_* functions) to CONOPT */
static int COI_CALLCONV Option(
   int                   NCALL,              /** number of callback call (provided by CONOPT) */
   double*               RVAL,               /** pointer to set the value of a real option */
   int*                  IVAL,               /** pointer to set the value of an integer option */
   int*                  LVAL,               /** pointer to set the value of a binary option */
   char*                 NAME,               /** pointer to set the name of the option */
   void*                 USRMEM              /**< user memory pointer (i.e. pointer to SCIP_NLPIPROBLEM) */
   )
{
   SCIP_NLPIPROBLEM* problem = (SCIP_NLPIPROBLEM*)USRMEM;

   assert(problem != NULL);

   switch( NCALL )
   {
      case 0: /* information about stopping criteria */
         strcpy(NAME, "LOCNPT");
         if( problem->verblevel >= 2 )
            *IVAL = 1;
         else
            *IVAL = 0;
         break;
      case 1: /* optimality tolerance */
         strcpy(NAME, "RTREDG");
         *RVAL = problem->opttol;
         break;
      default:
         strcpy(NAME,"");
   }

   return 0;
}

/** CONOPT callback to define structure of the Hessian of the Lagrangian */
static int COI_CALLCONV LagrStr(
   int                   HSRW[],             /**< row numbers of the lower triangular part of the Hessian */
   int                   HSCL[],             /**< column numbers of the lower triangular part of the Hessian; elements must be
                                                  sorted column-wise, and within each column, row-wise */
   int*                  NODRV,              /**< can be set to 1 if the derivatives could not be computed */
   int                   NUMVAR,             /**< number of variables as defined in coidef_numvar() (provided by CONOPT) */
   int                   NUMCON,             /**< number of constraints as defined in coidef_numcon() (provided by CONOPT) */
   int                   NHESS,              /**< number of nonzero elements in the Hessian (provided by CONOPT) */
   void*                 USRMEM              /**< user memory pointer (i.e. pointer to SCIP_NLPIPROBLEM) */
   )
{
   SCIP_NLPIPROBLEM* problem = (SCIP_NLPIPROBLEM*)USRMEM;
   const int* hessoffsets;
   const int* rows;
   int col = 0;

   assert(problem != NULL);

   if( SCIPnlpiOracleGetHessianLagSparsity(problem->scip, problem->oracle, &hessoffsets, &rows, TRUE) != SCIP_OKAY )
   {
      *NODRV = 1;
      return 0;
   }

   for( int i = 0; i < NHESS; i++ )
   {
      /* check if it's time to switch to the next row (this will also skip empty rows) */
      while( i == hessoffsets[col+1] )
         col++;

      HSRW[i] = rows[i];
      HSCL[i] = col;
   }

   return 0;
}

/** CONOPT callback to compute the Hessian of the Lagrangian
 *
 *  the Lagrangian is written as L = sum_{r in rows}U[r] * function(r)
 */
static int COI_CALLCONV LagrVal(
   const double          X[],                /**< point in which the Hessian should be computed (provided by CONOPT) */
   const double          U[],                /**< vector of weights on the individual constraints (provided by CONOPT) */
   const int             HSRW[],             /**< row numbers of the lower triangular part of the Hessian (provided by CONOPT) */
   const int             HSCL[],             /**< column numbers of the lower triangular part of the Hessian; elements must be
                                                  sorted column-wise, and within each column, row-wise (provided by CONOPT) */
   double                HSVL[],             /**< values of Hessian entries */
   int*                  NODRV,              /**< can be set to 1 if the derivatives could not be computed */
   int                   NUMVAR,             /**< number of variables as defined in coidef_numvar() (provided by CONOPT) */
   int                   NUMCON,             /**< number of constraints as defined in coidef_numcon() (provided by CONOPT) */
   int                   NHESS,              /**< number of nonzero elements in the Hessian (provided by CONOPT) */
   void*                 USRMEM              /**< user memory pointer (i.e. pointer to SCIP_NLPIPROBLEM) */
   )
{
   SCIP_NLPIPROBLEM* problem = (SCIP_NLPIPROBLEM*)USRMEM;
   SCIP_Real* hessvals;

   assert(problem != NULL);

   if( SCIPallocBlockMemoryArray(problem->scip, &hessvals, NHESS) != SCIP_OKAY )
   {
      *NODRV = 1;
      return 0;
   }

   /* TODO better handling for isnew */
   if( SCIPnlpiOracleEvalHessianLag(problem->scip, problem->oracle, X, TRUE, TRUE, U[NUMCON-1], U, hessvals, TRUE)
         != SCIP_OKAY )
   {
      SCIPfreeBlockMemoryArray(problem->scip, &hessvals, NHESS);
      *NODRV = 1;
      return 0;
   }

   for( int i = 0; i < NHESS; i++ )
      HSVL[i] = hessvals[i];

   SCIPfreeBlockMemoryArray(problem->scip, &hessvals, NHESS);

   return 0;
}

#ifdef DISABLED_CODE
/** CONOPT callback to provide interval information
 *
 *  The callback has three modes (indicated by the argument MODE):
 *    1: Only evaluate the bounds on the sum of the nonlinear terms in row ROWNO and return the interval in [GMIN, GMAX].
 *    2: Only evaluate the bounds on the nonlinear Jacobian elements in row ROWNO and return the intervals in [JMIN, JMAX].
 *    3: Perform both option 1 and 2. (Currently, MODE = 3 is not used).
 *
 *  TODO this is work in progress
 */
static int COI_CALLCONV FDInterval(
   const double*         XMIN,               /**< vector of lower bounds on variables (provided by CONOPT) */
   const double*         XMAX,               /**< vector of upper bounds on variables (provided by CONOPT) */
   double*               GMIN,               /**< pointer to store a lower bound on the function value in constraint ROWNO if MODE = 1 or 3 */
   double*               GMAX,               /**< pointer to store an upper bound on the function value in constraint ROWNO if MODE = 1 or 3 */
   double*               JMIN,               /**< if MODE = 2 or 3, lower bound on the derivative of row ROWNO with respect to variable I may be returned in JMIN(I) */
   double*               JAMX,               /**< if MODE = 2 or 3, upper bound on the derivative of row ROWNO with respect to variable I may be returned in JMIN(I) */
   int                   ROWNO,              /**< number of the current (nonlinear) constraint (provided by CONOPT) */
   const int*            JCNM,               /**< list of column numbers of nonlinear nonzeroes of the Jacobian in the current row, only defined when MODE = 2 or 3 */
   int                   MODE,               /**< indicator for mode of evaluation (provided by CONOPT) */
   double                PINF,               /**< CONOPT's plus infinity (provided by CONOPT) */
   int                   N,                  /**< number of variables in the model (provided by CONOPT) */
   int                   NJ,                 /**< number of nonlinear nonzeroes in the Jacobian in current row, only defined when MODE = 2 or 3 */
   void*                 USRMEM              /**< user memory pointer (i.e. pointer to SCIP_NLPIPROBLEM) */
   )
{
   SCIP_NLPIPROBLEM* problem = (SCIP_NLPIPROBLEM*)USRMEM;
   SCIP_NLPIORACLE* oracle;

   assert(problem != NULL);

   oracle = problem->oracle;
   assert(oracle != NULL);



   if( MODE == 1 )
   { /* evaluate the bounds on the sum of the nonlinear terms in row ROWNO and return the interval in [GMIN, GMAX] */

   }

   return 0;
}
#endif

/* NLPI local methods */

/** frees solution arrays and sets the solstat and termstat to unknown and other, resp. */
static
void invalidateSolution(
   SCIP_NLPIPROBLEM*     problem             /**< data structure of problem */
   )
{
   int nvars;
   int nconss;

   assert(problem != NULL);

   nvars = SCIPnlpiOracleGetNVars(problem->oracle);
   nconss = SCIPnlpiOracleGetNConstraints(problem->oracle);

   SCIPfreeBlockMemoryArrayNull(problem->scip, &(problem->lastprimal), nvars);
   SCIPfreeBlockMemoryArrayNull(problem->scip, &(problem->lastdualcons), nconss);
   SCIPfreeBlockMemoryArrayNull(problem->scip, &(problem->lastduallb), nvars);
   SCIPfreeBlockMemoryArrayNull(problem->scip, &(problem->lastdualub), nvars);

   problem->solstat  = SCIP_NLPSOLSTAT_UNKNOWN;
   problem->termstat = SCIP_NLPTERMSTAT_OTHER;
}

static SCIP_RETCODE initConopt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIDATA*        data,               /**< pointer to NLPIDATA structure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to NLPI problem structure */
)
{
   int COI_Error = 0; /* CONOPT error counter */
   int nrangeconss = 0;
   int nconss;
   int nvars;
   const int* jacoffsets;
   const int* hessoffsets;
   int nnlnz;
   const int* objgradnz;
   const SCIP_Bool* objnl;
   int nobjgradnz;
   int nobjgradnls;

   assert(data != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   // printf("\nInitialising conopt for problem:\n");
   // SCIPnlpiOraclePrintProblem(scip, problem->oracle, NULL);
   // printf("\n");

   if( data->solvetime == NULL )
      SCIP_CALL( SCIPcreateClock(scip, &(data->solvetime)) );

   nconss = SCIPnlpiOracleGetNConstraints(problem->oracle);
   nvars = SCIPnlpiOracleGetNVars(problem->oracle);

   /* count range constraints: because CONOPT doesn't support them directly, will need to add a slack variable for each ranged constraint */
   for( int i = 0; i < nconss; i++ )
   {
      SCIP_Real lhs = SCIPnlpiOracleGetConstraintLhs(problem->oracle, i);
      SCIP_Real rhs = SCIPnlpiOracleGetConstraintRhs(problem->oracle, i);

      if( !SCIPisInfinity(problem->scip, -lhs) && !SCIPisInfinity(problem->scip, rhs) && !SCIPisEQ(scip, lhs, rhs) )
         nrangeconss++;
   }

   COI_Error += COIDEF_EmptyCol(problem->CntVect, 1);

   /* inform CONOPT about problem sizes */
   COI_Error += COIDEF_NumVar(problem->CntVect, nvars + nrangeconss);
   COI_Error += COIDEF_NumCon(problem->CntVect, nconss + 1); /* objective counts as another constraint here */

   /* jacobian information */
   SCIP_CALL( SCIPnlpiOracleGetJacobianColSparsity(scip, problem->oracle, &jacoffsets, NULL, NULL, &nnlnz) );
   SCIP_CALL( SCIPnlpiOracleGetObjGradientNnz(scip, problem->oracle, &objgradnz, &objnl, &nobjgradnz, &nobjgradnls) );

   /* each slack var adds a Jacobian nnz; objective also counts as constraint */
   COI_Error += COIDEF_NumNz(problem->CntVect, jacoffsets != NULL ? jacoffsets[nvars] + nrangeconss + nobjgradnz :
         nrangeconss + nobjgradnz);

   /* Jacobian nonzeroes include those of constraints and objective */
   COI_Error += COIDEF_NumNlNz(problem->CntVect, nnlnz + nobjgradnls);

   /* hessian sparsity information */
   SCIP_CALL( SCIPnlpiOracleGetHessianLagSparsity(scip, problem->oracle, &hessoffsets, NULL, TRUE) );
   COI_Error += COIDEF_NumHess(problem->CntVect, hessoffsets[nvars]);

   /* tell CONOPT to minimise the objective (the oracle always gives a minimisation problem) */
   COI_Error += COIDEF_OptDir(problem->CntVect, -1);

   /* oracle gives objective as a constraint, hence use ObjCon (not ObjVar) here;
    * we treat objective as the last constraint, hence index nconss */
   COI_Error += COIDEF_ObjCon(problem->CntVect, nconss);

   /* register callback routines */
   COI_Error += COIDEF_Message(problem->CntVect, &Message);
   COI_Error += COIDEF_ErrMsg(problem->CntVect, &ErrMsg);
   COI_Error += COIDEF_Status(problem->CntVect, &Status);
   COI_Error += COIDEF_Solution(problem->CntVect, &Solution);
   COI_Error += COIDEF_ReadMatrix(problem->CntVect, &ReadMatrix);
   COI_Error += COIDEF_FDEval(problem->CntVect, &FDEval);
   COI_Error += COIDEF_Option(problem->CntVect, &Option);
   COI_Error += COIDEF_2DLagrStr(problem->CntVect, &LagrStr);
   COI_Error += COIDEF_2DLagrVal(problem->CntVect, &LagrVal);
   /* TODO complete and add the interval callback */

   COI_Error += COIDEF_FVincLin(problem->CntVect, 1);

   /* pass the problem pointer to CONOPT, so that it may be used in CONOPT callbacks */
   COI_Error += COIDEF_UsrMem(problem->CntVect, (void*)problem);

   /* register license */
   /* TODO make sure the license error gets printed by CONOPT if no valid license (even with reduced output) */
#if defined(LICENSE_INT_1) && defined(LICENSE_INT_2) && defined(LICENSE_INT_3) && defined(LICENSE_TEXT)
   COI_Error += COIDEF_License(problem->CntVect, LICENSE_INT_1, LICENSE_INT_2, LICENSE_INT_3, LICENSE_TEXT);
#endif

   if( COI_Error )
      SCIPinfoMessage(scip, NULL, "Error %d encountered during initialising CONOPT", COI_Error);

   return SCIP_OKAY;
}

static void handleConoptParam(
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   const SCIP_NLPPARAM   param               /**< NLP solve parameters */
)
{
   int COI_Error = 0; /* CONOPT error counter */

   assert(problem != NULL);

   if( param.warmstart )
   {
      SCIPdebugMsg(problem->scip, "warmstart parameter not supported by CONOPT interface yet. Ignored.\n");
   }
   if( param.lobjlimit > -SCIP_REAL_MAX )
   {
      SCIPwarningMessage(problem->scip, "lobjlimit parameter not supported by CONOPT interface yet. Ignored.\n");
   }
   if( param.fastfail )
   {
      SCIPdebugMsg(problem->scip, "fastfail parameter not supported by CONOPT interface yet. Ignored.\n");
   }

   /* options that we can set directly */
   COI_Error += COIDEF_ItLim(problem->CntVect, param.iterlimit);
   COI_Error += COIDEF_ResLim(problem->CntVect, param.timelimit);

   /* options that need to be handled in callbacks */
   problem->verblevel = param.verblevel;
   problem->opttol = param.opttol;

   /* TODO enable memory limits in NLPIs? */

   if( COI_Error )
      SCIPinfoMessage(problem->scip, NULL, "Errors encountered during setting parameters, %d", COI_Error);
}


/*
 * Callback methods of NLP solver interface
 */

#ifdef SCIP_DISABLED_CODE
/** copy method of NLP interface (called when SCIP copies plugins) */
static
SCIP_DECL_NLPICOPY(nlpiCopyXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/
#else
#define nlpiCopyConopt NULL
#endif

/** destructor of NLP interface to free nlpi data */
static
SCIP_DECL_NLPIFREE(nlpiFreeConopt)
{
   assert(nlpi != NULL);
   assert(nlpidata != NULL);
   assert(*nlpidata != NULL);

   if( (*nlpidata)->solvetime != NULL )
      SCIP_CALL( SCIPfreeClock(scip, &((*nlpidata)->solvetime)) );

   SCIPfreeBlockMemory(scip, nlpidata);
   assert(*nlpidata == NULL);

   return SCIP_OKAY;
}  /*lint !e715*/

#ifdef SCIP_DISABLED_CODE
/** gets pointer for NLP solver */
static
SCIP_DECL_NLPIGETSOLVERPOINTER(nlpiGetSolverPointerXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return NULL;  /*lint !e527*/
}  /*lint !e715*/
#else
#define nlpiGetSolverPointerConopt NULL
#endif

/** create a problem instance */
static
SCIP_DECL_NLPICREATEPROBLEM(nlpiCreateProblemConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   SCIP_CALL( SCIPallocClearBlockMemory(scip, problem) );
   assert(*problem != NULL);

   (*problem)->firstrun = TRUE;
   (*problem)->scip = scip;

   /* initialize oracle */
   SCIP_CALL( SCIPnlpiOracleCreate(scip, &(*problem)->oracle) );
   SCIP_CALL( SCIPnlpiOracleSetProblemName(scip, (*problem)->oracle, name) );

   coiCreate(&((*problem)->CntVect));

   if( (*problem)->CntVect == NULL )
   {
      SCIPerrorMessage("Could not create CONOPT control vector\n");
      return SCIP_ERROR;
   }

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &(*problem)->randnumgen, DEFAULT_RANDSEED, TRUE) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** free a problem instance */
static
SCIP_DECL_NLPIFREEPROBLEM(nlpiFreeProblemConopt)
{
   assert(nlpi     != NULL);
   assert(problem  != NULL);
   assert(*problem != NULL);

   invalidateSolution(*problem);

   if( (*problem)->oracle != NULL )
   {
      SCIP_CALL( SCIPnlpiOracleFree(scip, &(*problem)->oracle) );
   }

   SCIPfreeRandom(scip, &(*problem)->randnumgen);
   SCIPfreeMemoryArrayNull(scip, &(*problem)->initguess);

   coiFree(&((*problem)->CntVect));

   SCIPfreeBlockMemory(scip, problem);
   *problem = NULL;

   return SCIP_OKAY;
}  /*lint !e715*/

#ifdef SCIP_DISABLED_CODE
/** gets pointer to solver-internal problem instance */
static
SCIP_DECL_NLPIGETPROBLEMPOINTER(nlpiGetProblemPointerXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return NULL;  /*lint !e527*/
}  /*lint !e715*/
#else
#define nlpiGetProblemPointerConopt NULL
#endif

/** add variables */
static
SCIP_DECL_NLPIADDVARS(nlpiAddVarsConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleAddVars(scip, problem->oracle, nvars, lbs, ubs, varnames) );

   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/


/** add constraints */
static
SCIP_DECL_NLPIADDCONSTRAINTS(nlpiAddConstraintsConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleAddConstraints(scip, problem->oracle, nconss, lhss, rhss,
         nlininds, lininds, linvals, exprs, names) );

   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets or overwrites objective, a minimization problem is expected */
static
SCIP_DECL_NLPISETOBJECTIVE(nlpiSetObjectiveConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   if( expr != NULL || SCIPnlpiOracleIsConstraintNonlinear(problem->oracle, -1) )
      problem->firstrun = TRUE;

   SCIP_CALL( SCIPnlpiOracleSetObjective(scip, problem->oracle, constant, nlins, lininds, linvals, expr) );

   invalidateSolution(problem);

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change variable bounds */
static
SCIP_DECL_NLPICHGVARBOUNDS(nlpiChgVarBoundsConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgVarBounds(scip, problem->oracle, nvars, indices, lbs, ubs) );

   invalidateSolution(problem);

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change constraint bounds */
static
SCIP_DECL_NLPICHGCONSSIDES(nlpiChgConsSidesConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgConsSides(scip, problem->oracle, nconss, indices, lhss, rhss) );

   invalidateSolution(problem);

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** delete a set of variables */
static
SCIP_DECL_NLPIDELVARSET(nlpiDelVarSetConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleDelVarSet(scip, problem->oracle, dstats) );

   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** delete a set of constraints */
static
SCIP_DECL_NLPIDELCONSSET(nlpiDelConstraintSetConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleDelConsSet(scip, problem->oracle, dstats) );

   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** changes (or adds) linear coefficients in a constraint or objective */
static
SCIP_DECL_NLPICHGLINEARCOEFS(nlpiChgLinearCoefsConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgLinearCoefs(scip, problem->oracle, idx, nvals, varidxs, vals) );

   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** replaces the expression tree of a constraint or objective */
static
SCIP_DECL_NLPICHGEXPR(nlpiChgExprConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgExpr(scip, problem->oracle, idxcons, expr) );

   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change the constant offset in the objective */
static
SCIP_DECL_NLPICHGOBJCONSTANT(nlpiChgObjConstantConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgObjConstant(scip, problem->oracle, objconstant) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets initial guess */
static
SCIP_DECL_NLPISETINITIALGUESS(nlpiSetInitialGuessConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   if( primalvalues != NULL )
   {
      if( !problem->initguess )
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &problem->initguess, primalvalues, SCIPnlpiOracleGetNVars(problem->oracle)) );
      else
         BMScopyMemoryArray(problem->initguess, primalvalues, SCIPnlpiOracleGetNVars(problem->oracle));
   }
   else
      SCIPfreeMemoryArrayNull(scip, &problem->initguess);

   return SCIP_OKAY;
}  /*lint !e715*/

/** try to solve NLP
 *
 * Note that SCIP will already have reset a timelimit of SCIP_REAL_MAX to the time remaining for the SCIP solve in SCIPnlpiSolve().
 */
static
SCIP_DECL_NLPISOLVE(nlpiSolveConopt)
{
   SCIP_NLPIDATA* data;
#ifdef SCIP_DEBUG
   int COI_Error = 0; /* CONOPT error counter */
#endif

   assert(nlpi != NULL);
   assert(problem != NULL);

   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);

   SCIPdebugMsg(scip, "solve with parameters " SCIP_NLPPARAM_PRINT(param));

   SCIP_CALL( SCIPnlpiOracleResetEvalTime(scip, problem->oracle) );

   if( param.timelimit == 0.0 )
   {
      /* there is nothing we can do if we are not given any time */
      problem->niterations = 0;
      problem->solvetime = 0.0;
      problem->termstat = SCIP_NLPTERMSTAT_TIMELIMIT;
      problem->solstat = SCIP_NLPSOLSTAT_UNKNOWN;

      return SCIP_OKAY;
   }

   problem->niterations = -1;
   problem->solvetime  = -1.0;

   /* set CONOPT parameters */
   handleConoptParam(problem, param);

   /* initialize Conopt data if necessary */
   if( problem->firstrun )
   {
      initConopt(scip, data, problem);
      problem->firstrun = FALSE;
   }

   /* measure time */
   SCIPresetClock(scip, data->solvetime);
   SCIPstartClock(scip, data->solvetime);

   /* optimize */
#ifdef SCIP_DEBUG
   COI_Error = COI_Solve(problem->CntVect);

   /* CONOPT may return either a positive error code, which is one of its own error codes,
    * or a negative error code that is a SCIP_RETCODE returned from one of the callbacks */
   if( COI_Error )
      SCIPinfoMessage(scip, NULL, "Errors encountered in CONOPT during solution, %d", COI_Error);
#else
   COI_Solve(problem->CntVect);
#endif

   /* store statistics (some statistics are passed back to SCIP by the Status callback) */
   problem->solvetime = SCIPgetClockTime(scip, data->solvetime);

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gives solution status */
static
SCIP_DECL_NLPIGETSOLSTAT(nlpiGetSolstatConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   return problem->solstat;
}  /*lint !e715*/

/** gives termination reason */
static
SCIP_DECL_NLPIGETTERMSTAT(nlpiGetTermstatConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   return problem->termstat;
}  /*lint !e715*/

/** gives primal and dual solution values */
static
SCIP_DECL_NLPIGETSOLUTION(nlpiGetSolutionConopt)
{
   assert(problem != NULL);

   if( primalvalues != NULL )
      *primalvalues = problem->lastprimal;

   if( consdualvalues != NULL )
      *consdualvalues = problem->lastdualcons;

   if( varlbdualvalues != NULL )
      *varlbdualvalues = problem->lastduallb;

   if( varubdualvalues != NULL )
      *varubdualvalues = problem->lastdualub;

   if( objval != NULL )
      *objval = problem->objval;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gives solve statistics */
static
SCIP_DECL_NLPIGETSTATISTICS(nlpiGetStatisticsConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(statistics != NULL);

   statistics->niterations = problem->niterations;
   statistics->totaltime = problem->solvetime;
   statistics->evaltime = SCIPnlpiOracleGetEvalTime(scip, problem->oracle);
   statistics->consviol = SCIP_INVALID; /* TODO currently unavailable in CONOPT; might change */
   statistics->boundviol = 0.0;

   return SCIP_OKAY;
}  /*lint !e715*/

/*
 * NLP solver interface specific interface methods
 */

/** create solver interface for Conopt solver and includes it into SCIP */
SCIP_RETCODE SCIPincludeNlpSolverConopt(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NLPIDATA* nlpidata;

   /* create Conopt solver interface data */
   SCIP_CALL( SCIPallocClearBlockMemory(scip, &nlpidata) );

   /* create and include solver interface */
   SCIP_CALL( SCIPincludeNlpi(scip,
         NLPI_NAME, NLPI_DESC, NLPI_PRIORITY,
         nlpiCopyConopt, nlpiFreeConopt, nlpiGetSolverPointerConopt,
         nlpiCreateProblemConopt, nlpiFreeProblemConopt, nlpiGetProblemPointerConopt,
         nlpiAddVarsConopt, nlpiAddConstraintsConopt, nlpiSetObjectiveConopt,
         nlpiChgVarBoundsConopt, nlpiChgConsSidesConopt, nlpiDelVarSetConopt, nlpiDelConstraintSetConopt,
         nlpiChgLinearCoefsConopt, nlpiChgExprConopt, nlpiChgObjConstantConopt,
         nlpiSetInitialGuessConopt, nlpiSolveConopt,
         nlpiGetSolstatConopt, nlpiGetTermstatConopt, nlpiGetSolutionConopt, nlpiGetStatisticsConopt,
         nlpidata) );

   return SCIP_OKAY;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   miniisc.c
 * @brief  find a minimum IIS cover
 * @author Marc Pfetsch
 */

#include <string.h>
#include <scip/scipdefplugins.h>
#include <lpi/lpi.h>

#include "benders.h"
#include "readargs.h"

/* default parameters */
#define DEFAULT_SOLVEMASTERAPPROX      FALSE      /**< Solve master problem approximately? */
#define DEFAULT_MASTERGAPLIMIT         0.1        /**< gap bound for approximately solving the master problem */
#define DEFAULT_REOPTIMIZATION         FALSE      /**< Use reoptimization to solve master problem? */
#define DEFAULT_MASTERSTALLNODES       5000L      /**< stall nodes for the master problem */

/** data needed for cut generation */
struct BENDERS_Data
{
   SCIP_LPI*             lp;                 /**< alternative polyhedron */
   int                   m;                  /**< number of constraints considered */
};


/* Macro for setting parameters in LPI */
#define SCIP_CALL_PARAM(x) /*lint -e527 */ do                                                   \
{                                                                                               \
   SCIP_RETCODE _restat_;                                                                       \
   if ( (_restat_ = (x)) != SCIP_OKAY && (_restat_ != SCIP_PARAMETERUNKNOWN) )                  \
   {                                                                                            \
      SCIPerrorMessage("[%s:%d] Error <%d> in function call\n", __FILE__, __LINE__, _restat_);  \
      SCIPABORT();                                                                              \
      return _restat_;                                                                          \
   }                                                                                            \
}                                                                                               \
while ( FALSE )


/** Fix variable @a ind to 0 */
static
SCIP_RETCODE fixAltLPVariable(
   SCIP_LPI*             lp,                 /**< alternative LP */
   int                   ind                 /**< variable that should be fixed to 0 */
   )
{
   SCIP_Real lb = 0.0;
   SCIP_Real ub = 0.0;

   /* change bounds */
   SCIP_CALL( SCIPlpiChgBounds(lp, 1, &ind, &lb, &ub) );

   return SCIP_OKAY;
}


/** fix variables in @a S to 0 */
static
SCIP_RETCODE fixAltLPVariables(
   SCIP*                 masterscip,         /**< SCIP pointer */
   int                   nmastervars,        /**< number of variables in master */
   SCIP_Bool*            S,                  /**< indices to fix */
   SCIP_LPI*             lp                  /**< alternative LP */
   )
{
   SCIP_Real* lb = NULL;
   SCIP_Real* ub = NULL;
   int* indices = NULL;
   int cnt = 0;
   int j;

   assert( masterscip != NULL );
   assert( S != NULL );
   assert( lp != NULL );

   SCIP_CALL( SCIPallocBufferArray(masterscip, &lb, nmastervars) );
   SCIP_CALL( SCIPallocBufferArray(masterscip, &ub, nmastervars) );
   SCIP_CALL( SCIPallocBufferArray(masterscip, &indices, nmastervars) );

   /* collect bounds to be changed */
   for (j = 0; j < nmastervars; ++j)
   {
      if ( S[j] )
      {
         indices[cnt] = j;
         lb[cnt] = 0.0;
         ub[cnt] = 0.0;
         ++cnt;
      }
   }

   /* change bounds */
   if ( cnt > 0 )
   {
      SCIP_CALL( SCIPlpiChgBounds(lp, cnt, indices, lb, ub) );
   }

   SCIPfreeBufferArray(masterscip, &indices);
   SCIPfreeBufferArray(masterscip, &ub);
   SCIPfreeBufferArray(masterscip, &lb);

   return SCIP_OKAY;
}

/** unfix variables in @a S */
static
SCIP_RETCODE unfixAltLPVariables(
   SCIP*                 masterscip,         /**< SCIP pointer */
   int                   nmastervars,        /**< number of variables in master */
   SCIP_Bool*            S,                  /**< indices to fix */
   SCIP_LPI*             lp                  /**< alternative LP */
   )
{
   SCIP_Real* lb = NULL;
   SCIP_Real* ub = NULL;
   int* indices = NULL;
   int cnt = 0;
   int j;

   assert( masterscip != NULL );
   assert( S != NULL );
   assert( lp != NULL );

   SCIP_CALL( SCIPallocBufferArray(masterscip, &lb, nmastervars) );
   SCIP_CALL( SCIPallocBufferArray(masterscip, &ub, nmastervars) );
   SCIP_CALL( SCIPallocBufferArray(masterscip, &indices, nmastervars) );

   /* collect bounds to be changed */
   for (j = 0; j < nmastervars; ++j)
   {
      if ( S[j] )
      {
         indices[cnt] = j;
         lb[cnt] = 0.0;
         ub[cnt] = SCIPlpiInfinity(lp);
         ++cnt;
      }
   }

   /* change bounds */
   if ( cnt > 0 )
   {
      SCIP_CALL( SCIPlpiChgBounds(lp, cnt, indices, lb, ub) );
   }

   SCIPfreeBufferArray(masterscip, &indices);
   SCIPfreeBufferArray(masterscip, &ub);
   SCIPfreeBufferArray(masterscip, &lb);

   return SCIP_OKAY;
}

/** Check whether the given LP is infeasible
 *
 *  If @a primal is false we assume that the problem is <em>dual feasible</em>, e.g., the problem
 *  was only changed by fixing bounds!
 *
 *  This is the workhorse for all methods that have to solve the alternative LP. We try in several
 *  ways to recover from possible stability problems.
 *
 *  @pre It is assumed that all parameters for the alternative LP are set.
 */
static
SCIP_RETCODE checkAltLPInfeasible(
   SCIP*                 masterscip,         /**< SCIP pointer */
   SCIP_LPI*             lp,                 /**< LP */
   SCIP_Bool             primal,             /**< whether we are using the primal or dual simplex */
   SCIP_Bool*            infeasible,         /**< output: whether the LP is infeasible */
   SCIP_Bool*            error               /**< output: whether an error occured */
   )
{
   SCIP_RETCODE retcode;

   assert( masterscip != NULL );
   assert( lp != NULL );
   assert( infeasible != NULL );
   assert( error != NULL );

   *error = FALSE;

   /* solve LP */
   if ( primal )
      retcode = SCIPlpiSolvePrimal(lp);  /* use primal simplex */
   else
      retcode = SCIPlpiSolveDual(lp);    /* use dual simplex */

   if ( retcode == SCIP_LPERROR )
   {
      *error = TRUE;
      return SCIP_OKAY;
   }
   SCIP_CALL( retcode );

   /* resolve if LP is not stable */
   if ( ! SCIPlpiIsStable(lp) )
   {
      SCIP_CALL_PARAM( SCIPlpiSetIntpar(lp, SCIP_LPPAR_FROMSCRATCH, TRUE) );
      SCIP_CALL_PARAM( SCIPlpiSetIntpar(lp, SCIP_LPPAR_PRESOLVING, FALSE) );
      SCIPwarningMessage(masterscip, "Numerical problems, retrying ...\n");

      /* re-solve LP */
      if ( primal )
         retcode = SCIPlpiSolvePrimal(lp);  /* use primal simplex */
      else
         retcode = SCIPlpiSolveDual(lp);    /* use dual simplex */

      /* reset parameters */
      SCIP_CALL_PARAM( SCIPlpiSetIntpar(lp, SCIP_LPPAR_FROMSCRATCH, FALSE) );
      SCIP_CALL_PARAM( SCIPlpiSetIntpar(lp, SCIP_LPPAR_PRESOLVING, TRUE) );

      if ( retcode == SCIP_LPERROR )
      {
         *error = TRUE;
         return SCIP_OKAY;
      }
      SCIP_CALL( retcode );
   }

   /* check whether we are in the paradoxical situation that
    * - the primal is not infeasible
    * - the primal is not unbounded
    * - the LP is not optimal
    * - we have a primal ray
    *
    * If we ran the dual simplex algorithm, then we run again with the primal simplex
    */
   if ( ! SCIPlpiIsPrimalInfeasible(lp) && ! SCIPlpiIsPrimalUnbounded(lp) && ! SCIPlpiIsOptimal(lp) && SCIPlpiExistsPrimalRay(lp) && ! primal )
   {
      SCIPwarningMessage(masterscip, "The dual simplex produced a primal ray. Retrying with primal ...\n");

      /* the following settings might be changed: */
      SCIP_CALL_PARAM( SCIPlpiSetIntpar(lp, SCIP_LPPAR_FROMSCRATCH, TRUE) );
      SCIP_CALL_PARAM( SCIPlpiSetIntpar(lp, SCIP_LPPAR_PRESOLVING, TRUE) );
      SCIP_CALL_PARAM( SCIPlpiSetIntpar(lp, SCIP_LPPAR_SCALING, TRUE) );

      SCIP_CALL( SCIPlpiSolvePrimal(lp) );   /* use primal simplex */

      /* reset parameters */
      SCIP_CALL_PARAM( SCIPlpiSetIntpar(lp, SCIP_LPPAR_FROMSCRATCH, FALSE) );
      SCIP_CALL_PARAM( SCIPlpiSetIntpar(lp, SCIP_LPPAR_PRESOLVING, TRUE) );
      SCIP_CALL_PARAM( SCIPlpiSetIntpar(lp, SCIP_LPPAR_SCALING, TRUE) );
   }

   /* examine LP solution status */
   if ( SCIPlpiIsPrimalInfeasible(lp) )     /* the LP is provably infeasible */
   {
      assert( ! SCIPlpiIsPrimalUnbounded(lp) );   /* can't be unbounded or optimal */
      assert( ! SCIPlpiIsOptimal(lp) );           /* if it is infeasible! */
      *infeasible = TRUE;                         /* LP is infeasible */
      return SCIP_OKAY;
   }
   else
   {
      /* By assumption the dual is feasible if the dual simplex is run, therefore
       * the status has to be primal unbounded or optimal. */
      if ( ! SCIPlpiIsPrimalUnbounded(lp) && ! SCIPlpiIsOptimal(lp) )
      {
         /* We have a status different from unbounded or optimal. This should not be the case ... */
         if (primal)
            SCIPwarningMessage(masterscip, "Primal simplex returned with unknown status: %d\n", SCIPlpiGetInternalStatus(lp));
         else
            SCIPwarningMessage(masterscip, "Dual simplex returned with unknown status: %d\n", SCIPlpiGetInternalStatus(lp));

         /* SCIP_CALL( SCIPlpiWriteLP(lp, "debug.lp") ); */
         *error = TRUE;
         return SCIP_OKAY;
      }
   }

   /* at this point we have a feasible solution */
   *infeasible = FALSE;
   return SCIP_OKAY;
}


/** produce Benders cuts from the alternative polyhedron
 *
 *  input:
 *   - masterscip:       SCIP pointer of Benders master problem
 *   - nmastervars:      number of variables in master problem
 *   - mastervars:       variables in master problem
 *   - mastersolution:   solution of Benders master problem
 *   - data:             user data for oracle
 *   - timelimit:        time limit for subproblem
 *   - ntotalcuts:       total number of cuts
 *  output:
 *   - ncuts:            number of cuts added
 *   - status:           status
 *
 *  @todo apply time limit
 */
static
BENDERS_CUTORACLE(cutoracle)
{  /*lint --e{715}*/
#ifdef SCIP_DEBUG
   char name[SCIP_MAXSTRLEN];
#endif
   SCIP_LPI* lp;
   SCIP_Real* primsol;
   SCIP_Real value = 0.0;
   SCIP_Bool* S;
   int size = 0;
   int step = 0;
   int ncols;
   int j;

   assert( masterscip != NULL );
   assert( data != NULL );
   assert( mastersolution != NULL );
   assert( ncuts != NULL );
   assert( status != NULL );
   assert( data->lp != NULL );
   assert( data->m == nmastervars );

   lp = data->lp;

   *ncuts = 0;
   *status = BENDERS_STATUS_UNKNOWN;

   SCIP_CALL( SCIPlpiGetNCols(lp, &ncols) );
   SCIP_CALL( SCIPallocBufferArray(masterscip, &primsol, ncols) );
   assert( nmastervars <= ncols );

   /* init set S */
   SCIP_CALL( SCIPallocClearBufferArray(masterscip, &S, nmastervars) );
   for (j = 0; j < nmastervars; ++j)
   {
      assert( SCIPisFeasIntegral(masterscip, mastersolution[j]) );
      if ( mastersolution[j] > 0.5 )
      {
         S[j] = TRUE;
         ++size;
         value += SCIPvarGetObj(mastervars[j]);
      }
   }
   SCIP_CALL( fixAltLPVariables(masterscip, nmastervars, S, lp) );

   do
   {
      SCIP_CONS* cons;
      SCIP_VAR** vars;
      SCIP_Bool infeasible;
      SCIP_Real candobj = -1.0;
      SCIP_Bool error;
      int sizeIIS = 0;
      int candidate = -1;
      int cnt = 0;

      if ( step == 0 )
      {
         /* the first LP is solved without warm start, after that we use a warmstart. */
         SCIP_CALL_PARAM( SCIPlpiSetIntpar(lp, SCIP_LPPAR_FROMSCRATCH, TRUE) );
         SCIP_CALL( checkAltLPInfeasible(masterscip, lp, TRUE, &infeasible, &error) );
         SCIP_CALL_PARAM( SCIPlpiSetIntpar(lp, SCIP_LPPAR_FROMSCRATCH, FALSE) );
      }
      else
         SCIP_CALL( checkAltLPInfeasible(masterscip, lp, FALSE, &infeasible, &error) );

      if ( error )
      {
         *status = BENDERS_STATUS_ERROR;
         break;
      }

      /* if the alternative polyhedron is infeasible, we found a cover */
      if ( infeasible )
      {
         /* if the problem is infeasible in the first step, we are successful */
         if ( step == 0 )
            *status = BENDERS_STATUS_SUCESS;

         SCIPdebugMessage("   size: %4d  produced possible cover with objective value %f.\n", size, value);
         break;
      }

      /* get solution of alternative LP */
      SCIP_CALL( SCIPlpiGetSol(lp, NULL, primsol, NULL, NULL, NULL) );

      /* find candidate for variable to add */
      for (j = 0; j < nmastervars; ++j)
      {
         /* check support of the solution, i.e., the corresponding IIS */
         if ( ! SCIPisFeasZero(masterscip, primsol[j]) )
         {
            assert( ! S[j] );
            ++sizeIIS;

            /* take first element */
            if ( candidate < 0 )
            {
               candidate = j;
               candobj = SCIPvarGetObj(mastervars[j]);
            }
         }
      }

      /* check for error */
      if ( candidate < 0 )
      {
         /* Because of numerical problem it might happen that the solution primsol above is zero
          * within the tolerances. In this case we quit. */
         break;
      }
      assert( candidate >= 0 );
      assert( ! S[candidate] );
      assert( sizeIIS > 0 );

      SCIPdebugMessage("   size: %4d  add %4d with objective value %6g and alt-LP solution value %-8.4g  (IIS size: %4d).\n",
         size, candidate, candobj, primsol[candidate], sizeIIS);

      /* update new set S */
      S[candidate] = TRUE;
      ++size;
      value += candobj;

      SCIP_CALL( SCIPallocBufferArray(masterscip, &vars, nmastervars) );

      /* collect variables corresponding to support to cut */
      for (j = 0; j < nmastervars; ++j)
      {
         /* check support of the solution, i.e., the corresponding IIS */
         if ( ! SCIPisFeasZero(masterscip, primsol[j]) )
            vars[cnt++] = mastervars[j];
      }
      assert( cnt == sizeIIS );

#ifdef SCIP_DEBUG
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "iis%d", (int) ntotalcuts + *ncuts);
      SCIP_CALL( SCIPcreateConsLogicor(masterscip, &cons, name, cnt, vars, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
#else
      SCIP_CALL( SCIPcreateConsLogicor(masterscip, &cons, "", cnt, vars, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
#endif

#ifdef SCIP_OUTPUT
      SCIP_CALL( SCIPprintCons(masterscip, cons, NULL) );
      SCIPinfoMessage(masterscip, NULL, ";\n");
#endif

      SCIP_CALL( SCIPaddCons(masterscip, cons) );
      SCIP_CALL( SCIPreleaseCons(masterscip, &cons) );

      SCIPfreeBufferArray(masterscip, &vars);

      ++(*ncuts);
      *status = BENDERS_STATUS_ADDEDCUT;

      /* fix chosen variable to 0 */
      SCIP_CALL( fixAltLPVariable(lp, candidate) );

      ++step;
   }
   while (step < nmastervars);

   SCIP_CALL( unfixAltLPVariables(masterscip, nmastervars, S, lp) );

   SCIPfreeBufferArray(masterscip, &S);
   SCIPfreeBufferArray(masterscip, &primsol);

   return SCIP_OKAY;
}


/** creates column in alternative polyhedron */
static
SCIP_RETCODE createAltLPColumn(
   SCIP*                 origscip,           /**< SCIP pointer */
   SCIP_LPI*             lp,                 /**< alternative LP */
   int                   nvars,              /**< number of variables in column */
   SCIP_VAR**            vars,               /**< variables for column */
   SCIP_Real*            vals,               /**< values for column */
   SCIP_Real             rhscoef,            /**< coefficient for first row */
   SCIP_Real             sign                /**< sign (+1,-1) for column */
   )
{
   SCIP_Real obj = 1.0;
   SCIP_Real lb = 0.0;
   SCIP_Real ub;
   SCIP_Real* matval;
   int* matind;
   int matbeg = 0;
   int cnt = 0;
   int v;

   assert( origscip != NULL );
   assert( vars != NULL );
   assert( vals != NULL );
   assert( SCIPisEQ(origscip, sign, 1.0) || SCIPisEQ(origscip, sign, -1.0) );

   if ( SCIPisInfinity(origscip, rhscoef) || SCIPisInfinity(origscip, -rhscoef) )
      return SCIP_OKAY;

   /* set up data for construction */
   SCIP_CALL( SCIPallocBufferArray(origscip, &matind, nvars + 1) );
   SCIP_CALL( SCIPallocBufferArray(origscip, &matval, nvars + 1) );

   /* handle first row */
   if ( ! SCIPisFeasZero(origscip, rhscoef) )
   {
      matind[cnt] = 0;
      matval[cnt++] = sign * rhscoef;
   }

   /* set up column */
   for (v = 0; v < nvars; ++v)
   {
      assert( vars[v] != NULL );
      if ( vals != NULL )
         matval[cnt] = vals[v] * sign;
      else
         matval[cnt] = sign;
      matind[cnt++] = SCIPvarGetIndex(vars[v]) + 1;
   }

   /* now add column */
   ub = SCIPlpiInfinity(lp);

   SCIP_CALL( SCIPlpiAddCols(lp, 1, &obj, &lb, &ub, NULL, cnt, &matbeg, matind, matval) );

   SCIPfreeBufferArray(origscip, &matval);
   SCIPfreeBufferArray(origscip, &matind);

   return SCIP_OKAY;
}


/** create alternative polyhedron */
static
SCIP_RETCODE createAltLP(
   SCIP*                 origscip,           /**< original SCIP instance */
   SCIP_LPI*             lp                  /**< alternative polyhedron */
   )
{
   SCIP_CONS** origconss;
   int norigconss;
   int c;
   int v;

   assert( origscip != NULL );
   assert( lp != NULL );

   origconss = SCIPgetConss(origscip);
   norigconss = SCIPgetNConss(origscip);

   for (c = 0; c < norigconss; ++c)
   {
      const char* origconshdlrname;
      SCIP_CONSHDLR* origconshdlr;
      SCIP_VAR** origconsvars;
      SCIP_CONS* origcons;
      int norigconsvars;

      origcons = origconss[c];
      assert( origcons != NULL );

      origconshdlr = SCIPconsGetHdlr(origcons);
      assert( origconshdlr != NULL );
      origconshdlrname = SCIPconshdlrGetName(origconshdlr);

      if ( strcmp(origconshdlrname, "linear") == 0 )
      {
         origconsvars = SCIPgetVarsLinear(origscip, origcons);
         norigconsvars = SCIPgetNVarsLinear(origscip, origcons);

         SCIP_CALL( createAltLPColumn(origscip, lp, norigconsvars, origconsvars, SCIPgetValsLinear(origscip, origcons), SCIPgetRhsLinear(origscip, origcons), 1.0) );
         SCIP_CALL( createAltLPColumn(origscip, lp, norigconsvars, origconsvars, SCIPgetValsLinear(origscip, origcons), SCIPgetLhsLinear(origscip, origcons), -1.0) );
      }
      else if ( strcmp(origconshdlrname, "setppc") == 0 )
      {
         origconsvars = SCIPgetVarsSetppc(origscip, origcons);
         norigconsvars = SCIPgetNVarsSetppc(origscip, origcons);

         switch ( SCIPgetTypeSetppc(origscip, origcons) )
         {
         case SCIP_SETPPCTYPE_PARTITIONING :
            SCIP_CALL( createAltLPColumn(origscip, lp, norigconsvars, origconsvars, NULL, 1.0, 1.0) );
            SCIP_CALL( createAltLPColumn(origscip, lp, norigconsvars, origconsvars, NULL, 1.0, -1.0) );
            break;
         case SCIP_SETPPCTYPE_PACKING :
            SCIP_CALL( createAltLPColumn(origscip, lp, norigconsvars, origconsvars, NULL, 1.0, 1.0) );
            break;
         case SCIP_SETPPCTYPE_COVERING :
            SCIP_CALL( createAltLPColumn(origscip, lp, norigconsvars, origconsvars, NULL, 1.0, -1.0) );
            break;
         }
      }
      else if ( strcmp(origconshdlrname, "logicor") == 0 )
      {
         origconsvars = SCIPgetVarsLogicor(origscip, origcons);
         norigconsvars = SCIPgetNVarsLogicor(origscip, origcons);

         SCIP_CALL( createAltLPColumn(origscip, lp, norigconsvars, origconsvars, NULL, 1.0, -1.0) );
      }
      else if ( strcmp(origconshdlrname, "knapsack") == 0 )
      {
	 SCIP_Longint* origweights;
	 SCIP_Real* consvals;

         origconsvars = SCIPgetVarsKnapsack(origscip, origcons);
         norigconsvars = SCIPgetNVarsKnapsack(origscip, origcons);

         /* copy Longint array to SCIP_Real array */
         origweights = SCIPgetWeightsKnapsack(origscip, origcons);
         SCIP_CALL( SCIPallocBufferArray(origscip, &consvals, norigconsvars) );

         for ( v = 0; v < norigconsvars; ++v )
            consvals[v] = (SCIP_Real) origweights[v];

         SCIP_CALL( createAltLPColumn(origscip, lp, norigconsvars, origconsvars, consvals, (SCIP_Real) SCIPgetCapacityKnapsack(origscip, origcons), 1.0) );

         SCIPfreeBufferArray(origscip, &consvals);
      }
      else if ( strcmp(origconshdlrname, "varbound") == 0 )
      {
	 SCIP_VAR* consvars[2];
	 SCIP_Real consvals[2];

         consvars[0] = SCIPgetVarVarbound(origscip, origcons);
         consvars[1] = SCIPgetVbdvarVarbound(origscip, origcons);

         consvals[0] = 1.0;
         consvals[1] = SCIPgetVbdcoefVarbound(origscip, origcons);

         SCIP_CALL( createAltLPColumn(origscip, lp, 2, consvars, consvals, SCIPgetRhsVarbound(origscip, origcons), 1.0) );
         SCIP_CALL( createAltLPColumn(origscip, lp, 2, consvars, consvals, SCIPgetLhsVarbound(origscip, origcons), -1.0) );
      }
      else
      {
         SCIPwarningMessage(origscip, "Cannot handle constraints of type <%s>.\n", origconshdlrname);
      }
   }
   return SCIP_OKAY;
}


/** solve minimum IIS cover problem */
static
SCIP_RETCODE solveMinIISC(
   const char*           filename,           /**< problem name */
   const char*           settingsname,       /**< name of parameter file (or NULL) */
   SCIP_Real             timelimit,          /**< time limit read from arguments */
   SCIP_Real             memlimit,           /**< memory limit read from arguments */
   int                   dispfreq            /**< display frequency */
   )
{
   char name[SCIP_MAXSTRLEN];
   BENDERS_DATA data;
   SCIP* masterscip;
   SCIP* origscip;
   SCIP_STATUS status;
   SCIP_LPI* lp;
   SCIP_Real lhs = -1.0;
   SCIP_Real rhs = -1.0;
   SCIP_VAR** origvars;
   SCIP_Real obj = 0.0;
   SCIP_Real lb = 0.0;
   SCIP_Real ub;
   int norigvars;
   int nrows = 0;
   int m = 0;
   int v;

   /* parameters */
   SCIP_Bool solvemasterapprox;
   SCIP_Longint masterstallnodes;
   SCIP_Real mastergaplimit;
   SCIP_Bool reoptimization;

   /* create master SCIP */
   SCIP_CALL( SCIPcreate(&masterscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(masterscip) );
   if ( getProblemName(filename, name, SCIP_MAXSTRLEN) == 0 )
   {
      SCIPerrorMessage("Cannot extract problem name for filename <%s>.\n", filename);
      return SCIP_ERROR;
   }
   SCIP_CALL( SCIPcreateProb(masterscip, name, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPsetObjsense(masterscip, SCIP_OBJSENSE_MINIMIZE) );

   SCIPinfoMessage(masterscip, NULL, "Finding a minimum IIS cover using a set covering approach.\n");
   SCIPinfoMessage(masterscip, NULL, "Implemented by Marc Pfetsch, 2015\n\n");

   SCIPprintVersion(masterscip, NULL);
   SCIPinfoMessage(masterscip, NULL, "\n");

   /* add parameters */
   SCIP_CALL( SCIPaddBoolParam(masterscip,
         "miniisc/solvemasterapprox",
         "Solve master problem approximately?",
         &solvemasterapprox, TRUE, DEFAULT_SOLVEMASTERAPPROX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(masterscip,
         "miniisc/mastergaplimit",
         "gap bound for approximately solving the master problem",
         &mastergaplimit, TRUE, DEFAULT_MASTERGAPLIMIT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(masterscip,
         "miniisc/masterstallnodes",
         "stall nodes for the master problem",
         &masterstallnodes, TRUE, DEFAULT_MASTERSTALLNODES, 0L, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(masterscip,
         "miniisc/reoptimization",
         "Use reoptimization to solve master problem?",
         &reoptimization, TRUE, DEFAULT_REOPTIMIZATION, NULL, NULL) );

   /* read parameters if required */
   if ( settingsname != NULL )
   {
      if ( SCIPfileExists(settingsname) )
      {
         SCIPinfoMessage(masterscip, NULL, "\nreading user parameter file <%s> ...\n\n", settingsname);
         SCIP_CALL( SCIPreadParams(masterscip, settingsname) );
         SCIP_CALL( SCIPwriteParams(masterscip, NULL, FALSE, TRUE) );
      }
      else
      {
         SCIPwarningMessage(masterscip, NULL, "\nparameter file <%s> not found - using default parameters.\n", settingsname);
      }
   }

   if ( ! SCIPisInfinity(masterscip, timelimit) )
      SCIPinfoMessage(masterscip, NULL, "limits/time = %f\n\n", timelimit);

   SCIPinfoMessage(masterscip, NULL, "Input file:\t%s\n", filename);
   SCIPinfoMessage(masterscip, NULL, "Problem name:\t%s\n\n", name);

   /* ----------------------------------------------------------------------------------------*/

   /* read instance to create alternative polyhedron */
   SCIP_CALL( SCIPcreate(&origscip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(origscip) );

   /* read problem */
   SCIP_CALL( SCIPreadProb(origscip, filename, NULL) );

   /* check that we have an LP */
   if ( SCIPgetNOrigBinVars(origscip) + SCIPgetNOrigIntVars(origscip) > 0 )
   {
      SCIPinfoMessage(masterscip, NULL, "ERROR: input file contains integer variables. The code only works for LPs.\n");
      return SCIP_ERROR;
   }

   /* ----------------------------------------------------------------------------------------*/

   /* init alternative polyhedron */
   SCIP_CALL( SCIPlpiCreate(&lp, SCIPgetMessagehdlr(masterscip), "altlp", SCIP_OBJSEN_MINIMIZE) );

   /* init parameters */
   SCIP_CALL_PARAM( SCIPlpiSetIntpar(lp, SCIP_LPPAR_FROMSCRATCH, FALSE) );
   SCIP_CALL_PARAM( SCIPlpiSetIntpar(lp, SCIP_LPPAR_PRESOLVING, TRUE) );
   SCIP_CALL_PARAM( SCIPlpiSetIntpar(lp, SCIP_LPPAR_SCALING, TRUE) );
   SCIP_CALL_PARAM( SCIPlpiSetIntpar(lp, SCIP_LPPAR_FASTMIP, FALSE) );

   /* add first row */
   SCIP_CALL( SCIPlpiAddRows(lp, 1, &lhs, &rhs, NULL, 0, NULL, NULL, NULL) );

   norigvars = SCIPgetNOrigVars(origscip);
   origvars = SCIPgetOrigVars(origscip);

   /* add rows for each variable */
   lhs = 0.0;
   rhs = 0.0;
   for (v = 0; v < norigvars; ++v)
   {
      SCIP_CALL( SCIPlpiAddRows(lp, 1, &lhs, &rhs, NULL, 0, NULL, NULL, NULL) );
   }
   SCIP_CALL( SCIPlpiGetNRows(lp, &nrows) );

   /* create alternative polyhedron */
   SCIP_CALL( createAltLP(origscip, lp) );

   /* get number of constraints */
   SCIP_CALL( SCIPlpiGetNCols(lp, &m) );

   /* add columns for bounds */
   ub = SCIPlpiInfinity(lp);
   for (v = 0; v < norigvars; ++v)
   {
      SCIP_Real val;
      SCIP_VAR* var;
      SCIP_Real matval[2];
      int matind[2];
      int matbeg = 0;
      int cnt = 0;

      var = origvars[v];
      assert( var != NULL );
      assert( 0 <= SCIPvarGetIndex(var) && SCIPvarGetIndex(var) < nrows );

      /* if the lower bound is finite */
      val = SCIPvarGetLbGlobal(var);
      if ( ! SCIPisInfinity(origscip, -val) )
      {
         if ( ! SCIPisZero(origscip, val) )
         {
            matind[cnt] = 0;
            matval[cnt++] = -val;
         }
         matind[cnt] = SCIPvarGetIndex(var) + 1;
         matval[cnt++] = -1.0;
         SCIP_CALL( SCIPlpiAddCols(lp, 1, &obj, &lb, &ub, NULL, cnt, &matbeg, matind, matval) );
      }

      /* if the upper bound is finite */
      cnt = 0;
      val = SCIPvarGetUbGlobal(var);
      if ( ! SCIPisInfinity(origscip, val) )
      {
         if ( ! SCIPisZero(origscip, val) )
         {
            matind[cnt] = 0;
            matval[cnt++] = val;
         }
         matind[cnt] = SCIPvarGetIndex(var) + 1;
         matval[cnt++] = 1.0;
         SCIP_CALL( SCIPlpiAddCols(lp, 1, &obj, &lb, &ub, NULL, cnt, &matbeg, matind, matval) );
      }
   }

   /* free SCIP instance */
   SCIP_CALL( SCIPfree(&origscip) );

#ifdef SCIP_OUTPUT
   SCIP_CALL( SCIPlpiWriteLP(lp, "alt.lp") );
#endif

   /* ----------------------------------------------------------------------------------------*/
   /* initialize master problem */
   for (v = 0; v < m; ++v)
   {
      SCIP_VAR* var;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "y%d", v);
      SCIP_CALL( SCIPcreateVar(masterscip, &var, name, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(masterscip, var) );
      SCIP_CALL( SCIPreleaseVar(masterscip, &var) );
   }

   /* run Benders algorithm */
   data.lp = lp;
   data.m = m;
   SCIP_CALL( runBenders(masterscip, cutoracle, &data, timelimit, memlimit, dispfreq, reoptimization, solvemasterapprox,
         masterstallnodes, mastergaplimit, SCIP_VERBLEVEL_NORMAL, &status) );

   SCIP_CALL( SCIPlpiFree(&lp) );

   SCIP_CALL( SCIPfree(&masterscip) );

   return SCIP_OKAY;
}




/** main function */
int
main(
   int                   argc,               /**< number of shell parameters */
   char**                argv                /**< array with shell parameters */
   )
{
   SCIP_RETCODE retcode;
   const char* filename;
   const char* settingsname;
   SCIP_Real timelimit;
   SCIP_Real memlimit;
   SCIP_Longint nodelimit;
   int dispfreq;

   retcode = readArguments(argc, argv, &filename, &settingsname, &timelimit, &memlimit, &nodelimit, &dispfreq);
   if ( retcode != SCIP_OKAY )
      return -1;
   assert( filename != NULL );

   /* read file */
   if ( ! SCIPfileExists(filename) )
   {
      SCIPerrorMessage("file <%s> does not exist.\n", filename);
      return -1;
   }

   retcode = solveMinIISC(filename, settingsname, timelimit, memlimit, dispfreq);
   if ( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   BMScheckEmptyMemory();

   return 0;
}

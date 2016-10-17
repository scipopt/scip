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

/**@file   prop_nlobbt.c
 * @brief  nlobbt propagator
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/prop_nlobbt.h"
#include "nlpi/nlpi.h"

#define PROP_NAME              "nlobbt"
#define PROP_DESC              "propagator template"
#define PROP_PRIORITY          -1100000
#define PROP_FREQ                     0
#define PROP_DELAY                 TRUE
#define PROP_TIMING            SCIP_PROPTIMING_AFTERLPLOOP

#define DEFAULT_MINNONCONVEXFRAC   0.20      /**< default minimum (nconvex nlrows)/(nonconvex nlrows) threshold to apply NLOBBT */
#define DEFAULT_MINLINEARFRAC      0.02      /**< default minimum (nconvex nlrows)/(linear nlrows) threshold to apply NLOBBT */
#define DEFAULT_FEASTOLFAC         0.01      /**< default factor for NLP feasibility tolerance */
#define DEFAULT_RELOBJTOLFAC       0.01      /**< default factor for NLP relative objective tolerance */
#define DEFAULT_ADDLPROWS          TRUE      /**< should (non-initial) LP rows be used? */

/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   SCIP_Bool             skipprop;           /**< should the propagator be skipped? */
   SCIP_Real             feastolfac;         /**< factor for NLP feasibility tolerance */
   SCIP_Real             relobjtolfac;       /**< factor for NLP relative objective tolerance */
   SCIP_Real             minnonconvexfrac;   /**< minimum (nconvex nlrows)/(nonconvex nlrows) threshold to apply NLOBBT */
   SCIP_Real             minlinearfrac;      /**< minimum (nconvex nlrows)/(linear nlrows) threshold to apply NLOBBT */
   SCIP_Bool             addlprows;          /**< should (non-initial) LP rows be used? */
};

/*
 * Local methods
 */

/** clears the propagator data */
static
SCIP_RETCODE propdataClear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert(propdata != NULL);

   propdata->skipprop = FALSE;

   return SCIP_OKAY;
}

/** clears and frees propagator data */
static
SCIP_RETCODE propdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA**       propdata            /**< propagator data */
   )
{
   assert(propdata != NULL);
   assert(*propdata != NULL);

   SCIP_CALL( propdataClear(scip, *propdata) );
   SCIPfreeBlockMemory(scip, propdata);
   *propdata = NULL;

   return SCIP_OKAY;
}

/** adds LP rows to a given NLPI problem */
static
SCIP_RETCODE addLpRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< nlpi problem */
   SCIP_HASHMAP*         var2idx,            /**< empty hash map to store mapping between variables and indices in nlpi
                                              *   problem */
   SCIP_ROW**            rows,               /**< rows to add */
   int                   nrows               /**< total number of rows to add */
   )
{
   int* lininds;
   int i;

   assert(nlpi != NULL);
   assert(nlpiprob != NULL);
   assert(var2idx != NULL);
   assert(rows != NULL || nrows == 0);

   if( nrows <= 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &lininds, SCIPgetNVars(scip)) );

   for( i = 0; i < nrows; ++i )
   {
      SCIP_Real* linvals;
      SCIP_Real lhs;
      SCIP_Real rhs;
      const char* name;
      int nlininds;
      int k;

      assert(rows[i] != NULL);
      assert(SCIProwGetNNonz(rows[i]) <= SCIPgetNVars(scip));

      for( k = 0; k < SCIProwGetNNonz(rows[i]); ++k )
      {
         SCIP_VAR* var;

         var = SCIPcolGetVar(SCIProwGetCols(rows[i])[k]);
         assert(var != NULL);
         assert(SCIPhashmapExists(var2idx, (void*)var));

         lininds[k] = (int)(size_t)SCIPhashmapGetImage(var2idx, (void*)var);
         assert(lininds[k] >= 0 && lininds[k] < SCIPgetNVars(scip));
      }

      linvals = SCIProwGetVals(rows[i]);
      nlininds = SCIProwGetNNonz(rows[i]);
      lhs = SCIProwGetLhs(rows[i]) - SCIProwGetConstant(rows[i]);
      rhs = SCIProwGetRhs(rows[i]) - SCIProwGetConstant(rows[i]);
      name = SCIProwGetName(rows[i]);

      SCIP_CALL( SCIPnlpiAddConstraints(nlpi, nlpiprob, 1, &lhs, &rhs, &nlininds, &lininds, &linvals, NULL,
            NULL, NULL, NULL, &name) );
   }

   SCIPfreeBufferArray(scip, &lininds);

   return SCIP_OKAY;
}

/** creates a convex NLP relaxation and stores it in a given NLPI problem;
 *
 *  @note function does not copy the objective
 **/
static
SCIP_RETCODE createNlpRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLROW**          nlrows,             /**< nonlinear rows */
   int                   nnlrows,            /**< total number of nonlinear rows */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< empty nlpi problem */
   SCIP_HASHMAP*         var2idx,            /**< empty hash map to store mapping between variables and indices in nlpi
                                              *   problem */
   int*                  nlcount             /**< array to store the number of occurrences of variables in convex and
                                              *   nonlinear rows */
   )
{
   SCIP_EXPRTREE** exprtrees;
   int** exprvaridxs;
   SCIP_QUADELEM** quadelems;
   int* nquadelems;
   SCIP_Real** linvals;
   int** lininds;
   int* nlininds;
   SCIP_Real* lhss;
   SCIP_Real* rhss;
   const char** names;
   SCIP_VAR** vars;
   int nvars;
   SCIP_Real* lbs;
   SCIP_Real* ubs;
   const char** varnames;
   int nconss;
   int i;

   assert(nlpiprob != NULL);
   assert(var2idx != NULL);
   assert(nlrows != NULL);
   assert(nnlrows > 0);
   assert(nlpi != NULL);
   assert(nlcount != NULL);

   BMSclearMemoryArray(nlcount, SCIPgetNVars(scip));
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   nconss = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &exprtrees, nnlrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &exprvaridxs, nnlrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, nnlrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nquadelems, nnlrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals, nnlrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lininds, nnlrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlininds, nnlrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &names, nnlrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lhss, nnlrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhss, nnlrows) );

   SCIP_CALL( SCIPallocBufferArray(scip, &lbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varnames, nvars) );

   /* create an unique mapping between variables and {0,..,nvars-1} */
   for( i = 0; i < nvars; ++i )
   {
      assert(vars[i] != NULL);
      SCIP_CALL( SCIPhashmapInsert(var2idx, (void*)vars[i], (void*)(size_t)i) );

      lbs[i] = SCIPvarGetLbLocal(vars[i]);
      ubs[i] = SCIPvarGetUbLocal(vars[i]);
      varnames[i] = SCIPvarGetName(vars[i]);
  }

   /** add variables */
   SCIP_CALL( SCIPnlpiAddVars(nlpi, nlpiprob, nvars, lbs, ubs, varnames) );
   SCIPfreeBufferArray(scip, &varnames);
   SCIPfreeBufferArray(scip, &ubs);
   SCIPfreeBufferArray(scip, &lbs);

   /* add convex nonlinear rows to NLPI problem */
   for( i = 0; i < nnlrows; ++i )
   {
      SCIP_Bool userhs;
      SCIP_Bool uselhs;
      int k;

      assert(nlrows[i] != NULL);

      uselhs = FALSE;
      userhs = FALSE;

      /* check curvature together with constraint sides of a nonlinear row */
      if( SCIPnlrowGetNQuadElems(nlrows[i]) == 0 && SCIPnlrowGetExprtree(nlrows[i]) == NULL )
      {
         uselhs = TRUE;
         userhs = TRUE;
      }
      else if( SCIPnlrowGetCurvature(nlrows[i]) == SCIP_EXPRCURV_CONVEX && !SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrows[i])) )
         userhs = TRUE;
      else if( SCIPnlrowGetCurvature(nlrows[i]) == SCIP_EXPRCURV_CONCAVE && !SCIPisInfinity(scip, SCIPnlrowGetLhs(nlrows[i])) )
         uselhs = TRUE;
      else
         continue;
      assert(uselhs || userhs);

      lhss[nconss] = uselhs ? SCIPnlrowGetLhs(nlrows[i]) - SCIPnlrowGetConstant(nlrows[i]) : -SCIPinfinity(scip);
      rhss[nconss] = userhs ? SCIPnlrowGetRhs(nlrows[i]) - SCIPnlrowGetConstant(nlrows[i]) :  SCIPinfinity(scip);
      names[nconss] = SCIPnlrowGetName(nlrows[i]);
      nlininds[nconss] = 0;
      lininds[nconss] = NULL;
      linvals[nconss] = NULL;
      nquadelems[nconss] = 0;
      quadelems[nconss] = NULL;
      exprtrees[nconss] = NULL;
      exprvaridxs[nconss] = NULL;

      /* copy linear part */
      if( SCIPnlrowGetNLinearVars(nlrows[i]) > 0 )
      {
         SCIP_VAR* var;

         nlininds[nconss] = SCIPnlrowGetNLinearVars(nlrows[i]);

         SCIP_CALL( SCIPallocBufferArray(scip, &lininds[nconss], nlininds[nconss]) );
         SCIP_CALL( SCIPallocBufferArray(scip, &linvals[nconss], nlininds[nconss]) );

         for( k = 0; k < nlininds[nconss]; ++k )
         {
            var = SCIPnlrowGetLinearVars(nlrows[i])[k];
            assert(var != NULL);
            assert(SCIPhashmapExists(var2idx, (void*)var));

            lininds[nconss][k] = (int)(size_t)SCIPhashmapGetImage(var2idx, (void*)var);
            assert(var == vars[lininds[nconss][k]]);
            linvals[nconss][k] = SCIPnlrowGetLinearCoefs(nlrows[i])[k];
         }
      }

      /* copy quadratic part */
      if( SCIPnlrowGetNQuadElems(nlrows[i]) > 0 )
      {
         SCIP_QUADELEM quadelem;
         SCIP_VAR* var1;
         SCIP_VAR* var2;

         nquadelems[nconss] = SCIPnlrowGetNQuadElems(nlrows[i]);
         SCIP_CALL( SCIPallocBufferArray(scip, &quadelems[nconss], nquadelems[nconss]) );

         for( k = 0; k < nquadelems[nconss]; ++k )
         {
            quadelem = SCIPnlrowGetQuadElems(nlrows[i])[k];

            var1 = SCIPnlrowGetQuadVars(nlrows[i])[quadelem.idx1];
            assert(var1 != NULL);
            assert(SCIPhashmapExists(var2idx, (void*)var1));

            var2 = SCIPnlrowGetQuadVars(nlrows[i])[quadelem.idx2];
            assert(var2 != NULL);
            assert(SCIPhashmapExists(var2idx, (void*)var2));

            quadelems[nconss][k].coef = quadelem.coef;
            quadelems[nconss][k].idx1 = (int)(size_t)SCIPhashmapGetImage(var2idx, (void*)var1);
            quadelems[nconss][k].idx2 = (int)(size_t)SCIPhashmapGetImage(var2idx, (void*)var2);

            /* update nlcount */
            ++nlcount[quadelems[nconss][k].idx1];
            if( quadelems[nconss][k].idx1 != quadelems[nconss][k].idx2 )
               ++nlcount[quadelems[nconss][k].idx2];
         }
      }

      /* copy expression tree */
      if( SCIPnlrowGetExprtree(nlrows[i]) != NULL )
      {
         SCIP_VAR* var;

         /* note that we don't need to copy the expression trees since they are copied internally when adding a
          * nonlinear constraint to the NLP solver through the NLPI
          */
         exprtrees[nconss] = SCIPnlrowGetExprtree(nlrows[i]);

         SCIP_CALL( SCIPallocBufferArray(scip, &exprvaridxs[nconss], SCIPexprtreeGetNVars(exprtrees[nconss])) );

         for( k = 0; k < SCIPexprtreeGetNVars(exprtrees[nconss]); ++k )
         {
            var = SCIPexprtreeGetVars(exprtrees[nconss])[k];
            assert(var != NULL);
            assert(SCIPhashmapExists(var2idx, (void*)var));

            exprvaridxs[nconss][k] = (int)(size_t)SCIPhashmapGetImage(var2idx, (void*)var);

            /* update nlcount */
            ++nlcount[exprvaridxs[nconss][k]];
         }
      }

      ++nconss;
   }
   assert(nconss > 0);

   /* pass all constraint information to nlpi */
   SCIP_CALL( SCIPnlpiAddConstraints(nlpi, nlpiprob, nconss, lhss, rhss, nlininds, lininds, linvals, nquadelems,
         quadelems, exprvaridxs, exprtrees, names) );

   /* free memory */
   for( i = nconss - 1; i >= 0; --i )
   {
      if( exprtrees[i] != NULL )
      {
         assert(exprvaridxs[i] != NULL);
         SCIPfreeBufferArray(scip, &exprvaridxs[i]);
      }

      if( nquadelems[i] > 0 )
      {
         assert(quadelems[i] != NULL);
         SCIPfreeBufferArray(scip, &quadelems[i]);
      }

      if( nlininds[i] > 0 )
      {
         assert(linvals[i] != NULL);
         assert(lininds[i] != NULL);
         SCIPfreeBufferArray(scip, &linvals[i]);
         SCIPfreeBufferArray(scip, &lininds[i]);
      }
   }

   SCIPfreeBufferArray(scip, &rhss);
   SCIPfreeBufferArray(scip, &lhss);
   SCIPfreeBufferArray(scip, &names);
   SCIPfreeBufferArray(scip, &nlininds);
   SCIPfreeBufferArray(scip, &lininds);
   SCIPfreeBufferArray(scip, &linvals);
   SCIPfreeBufferArray(scip, &nquadelems);
   SCIPfreeBufferArray(scip, &quadelems);
   SCIPfreeBufferArray(scip, &exprvaridxs);
   SCIPfreeBufferArray(scip, &exprtrees);

   return SCIP_OKAY;
}

/** checks whether it is worth to call nonlinear OBBT procedure */
static
SCIP_Bool isNlobbtApplicable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagation data */
   )
{
   SCIP_NLROW** nlrows;
   int nnonconvexnlrows;
   int nconvexnlrows;
   int nlinearnlrows;
   int nnlrows;
   int i;

   nlrows = SCIPgetNLPNlRows(scip);
   nnlrows = SCIPgetNNLPNlRows(scip);
   nnonconvexnlrows = 0;
   nconvexnlrows = 0;
   nlinearnlrows = 0;

   for( i = 0; i < nnlrows; ++i )
   {
      if( SCIPnlrowGetNQuadElems(nlrows[i]) == 0 && SCIPnlrowGetExprtree(nlrows[i]) == NULL )
         ++nlinearnlrows;
      else if( SCIPnlrowGetCurvature(nlrows[i]) == SCIP_EXPRCURV_CONVEX )
      {
         if( !SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrows[i])) )
            ++nconvexnlrows;
         if( !SCIPisInfinity(scip, -SCIPnlrowGetLhs(nlrows[i])) )
            ++nnonconvexnlrows;
      }
      else if( SCIPnlrowGetCurvature(nlrows[i]) == SCIP_EXPRCURV_CONCAVE )
      {
         if( !SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrows[i])) )
            ++nnonconvexnlrows;
         if( !SCIPisInfinity(scip, -SCIPnlrowGetLhs(nlrows[i])) )
            ++nconvexnlrows;
      }
      else
      {
         if( !SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrows[i])) )
            ++nnonconvexnlrows;
         if( !SCIPisInfinity(scip, -SCIPnlrowGetLhs(nlrows[i])) )
            ++nnonconvexnlrows;
      }
   }

   SCIPdebugMsg(scip, "nconvex=%d nnonconvex=%d nlinear=%d\n", nconvexnlrows, nnonconvexnlrows, nlinearnlrows);

   return nconvexnlrows > 0
      && (SCIPisGE(scip, nconvexnlrows, nnonconvexnlrows * propdata->minnonconvexfrac))
      && (SCIPisGE(scip, nconvexnlrows, nlinearnlrows * propdata->minlinearfrac));
}

/** filters variables which achieve their lower or dual bound in the current NLP solution */
static
SCIP_RETCODE filterCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< nlpi problem */
   SCIP_HASHMAP*         var2idx,            /**< hash map to store mapping between variables and indices in nlpi
                                              *   problem */
   int*                  nlcount,            /**< array to store the number of occurrences of variables in convex and
                                              *   nonlinear rows */
   int*                  status,             /**< array to store the status of each candidate (0: unfiltered, 1:
                                              *   filtered lb, 2 filtered ub, 3: filtered lb and ub)
                                              */
   SCIP_VAR**            cands,              /**< candidate array */
   int                   pos,                /**< current position in the cands array */
   SCIP_NLPSOLSTAT       nlpsolstat          /**< NLP solution status */
   )
{
   SCIP_Real* primal;
   int i;

   assert(nlpi != NULL);
   assert(nlpiprob != NULL);
   assert(cands != NULL);
   assert(var2idx != NULL);
   assert(nlcount != NULL);
   assert(status != NULL);
   assert(pos >= 0 && pos < SCIPgetNVars(scip));

   /* at least a feasible solution is needed */
   if( nlpsolstat > SCIP_NLPSOLSTAT_FEASIBLE )
      return SCIP_OKAY;

   SCIP_CALL( SCIPnlpiGetSolution(nlpi, nlpiprob, &primal, NULL, NULL, NULL) );

   for( i = pos + 1; i < SCIPgetNVars(scip) && nlcount[i] > 0; ++i )
   {
      SCIP_VAR* var;
      SCIP_Real val;
      int varidx;

      /* candidate has been filtered already */
      if( status[i] == 3 )
         continue;

      var = cands[i];
      assert(var != NULL && SCIPhashmapExists(var2idx, (void*)var));

      varidx = (int)(size_t)SCIPhashmapGetImage(var2idx, (void*)var);
      assert(SCIPgetVars(scip)[varidx] == var);
      val = primal[varidx];

      if( status[i] != 1 && !SCIPisInfinity(scip, -val) && SCIPisRelLE(scip, val, SCIPvarGetLbLocal(var)) )
      {
         SCIPdebugMsg(scip, "filter LB of %s in [%g,%g] with %g \n", SCIPvarGetName(var), SCIPvarGetLbLocal(var),
            SCIPvarGetUbLocal(var), val);
         status[i] = status[i] == 0 ? 1 : 3;
      }

      if( status[i] != 2 && !SCIPisInfinity(scip, val) && SCIPisRelGE(scip, val, SCIPvarGetUbLocal(var)) )
      {
         SCIPdebugMsg(scip, "filter UB of %s in [%g,%g] with %g \n", SCIPvarGetName(var), SCIPvarGetLbLocal(var),
            SCIPvarGetUbLocal(var), val);
         status[i] = status[i] == 0 ? 2 : 3;
      }
   }

   return SCIP_OKAY;
}

/** sets the objective function, solves the NLP, and tightens the given variable; after calling this function, the
 *  objective function is set to zero
 *
 *  @note function assumes that objective function is zero
 */
static
SCIP_RETCODE solveNlp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< nlpi problem */
   SCIP_VAR*             var,                /**< variable to propagate */
   int                   varidx,             /**< variable index in the nlpi */
   SCIP_BOUNDTYPE        boundtype,          /**< direction in which to optimize */
   SCIP_HASHMAP*         var2idx,            /**< hash map to store mapping between variables and indices in nlpi
                                              *   problem */
   int*                  nlcount,            /**< array to store the number of occurrences of variables in convex and
                                              *   nonlinear rows */
   int*                  status,             /**< array to store the status of each candidate (0: unfiltered, 1:
                                              *   filtered lb, 2 filtered ub, 3: filtered lb and ub)
                                              */
   SCIP_VAR**            cands,              /**< candidate array */
   int                   pos,                /**< current position in cands array */
   SCIP_RESULT*          result              /**< pointer to store result */
   )
{
   SCIP_Real* primal;
   SCIP_Real obj;

#ifdef SCIP_DEBUG
   SCIP_Real oldlb;
   SCIP_Real oldub;

   oldlb = SCIPvarGetLbLocal(var);
   oldub = SCIPvarGetUbLocal(var);
#endif

   assert(nlpi != NULL);
   assert(nlpiprob != NULL);
   assert(var != NULL);
   assert(varidx >= 0 && varidx < SCIPgetNVars(scip));
   assert(result != NULL && *result != SCIP_CUTOFF);

   /* set corresponding objective coefficient and solve NLP */
   obj = boundtype == SCIP_BOUNDTYPE_LOWER ? 1.0 : -1.0;
   SCIP_CALL( SCIPnlpiSetObjective(nlpi, nlpiprob, 1, &varidx, &obj, 0, NULL, NULL, NULL, 0.0) );
   SCIP_CALL( SCIPnlpiSolve(nlpi, nlpiprob) );

   /* filter bound candidates first, otherwise we do not have access to the primal solution values */
   SCIP_CALL( filterCands(scip, nlpi, nlpiprob, var2idx, nlcount, status, cands, pos, SCIPnlpiGetSolstat(nlpi, nlpiprob)) );

   /* try to tighten variable bound */
   if( SCIPnlpiGetSolstat(nlpi, nlpiprob) <= SCIP_NLPSOLSTAT_LOCOPT )
   {
      SCIP_Bool tightened;
      SCIP_Bool infeasible;

      SCIP_CALL( SCIPnlpiGetSolution(nlpi, nlpiprob, &primal, NULL, NULL, NULL) );

      if( boundtype == SCIP_BOUNDTYPE_LOWER )
      {
         SCIP_CALL( SCIPtightenVarLb(scip, var, primal[varidx], FALSE, &infeasible, &tightened) );
      }
      else
      {
         SCIP_CALL( SCIPtightenVarUb(scip, var, primal[varidx], FALSE, &infeasible, &tightened) );
      }

      if( infeasible )
      {
         SCIPdebugMsg(scip, "detect infeasibility after propagating %s\n", SCIPvarGetName(var));
         *result = SCIP_CUTOFF;
      }
      else if( tightened )
      {
         SCIP_Real lb;
         SCIP_Real ub;

         *result = SCIP_REDUCEDDOM;

         /* update bounds in NLP */
         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);
         SCIP_CALL( SCIPnlpiChgVarBounds(nlpi, nlpiprob, 1, &varidx, &lb, &ub) );

#ifdef SCIP_DEBUG
         SCIPdebugMsg(scip, "tightened bounds of %s from [%g,%g] to [%g,%g]\n", SCIPvarGetName(var), oldlb, oldub, lb, ub);
#endif
      }
   }

   /* reset objective function */
   obj = 0.0;
   SCIP_CALL( SCIPnlpiSetObjective(nlpi, nlpiprob, 1, &varidx, &obj, 0, NULL, NULL, NULL, 0.0) );

   return SCIP_OKAY;
}

/** main method of the propagator
 *
 *  creates a convex NLP relaxation and solves the OBBT-NLPs for each possible candidate;
 *  binary and variables with a small domain will be ignored to reduce the computational cost of the propagator; after
 *  solving each NLP we filter out all variable candidates which are on their lower or upper bound; candidates with a
 *  larger number of occurrences are preferred
 *
 *  @note propagor will be disabled if it did not find any reductions
 */
static
SCIP_RETCODE applyNlobbt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagation data */
   SCIP_RESULT*          result              /**< pointer to store result */
   )
{
   SCIP_VAR** cands;
   SCIP_NLPIPROBLEM* nlpiprob;
   SCIP_HASHMAP* var2idx;
   SCIP_NLPI* nlpi;
   int* status;
   int* nlcount;
   int nlprows;
   int i;

   assert(result != NULL);
   assert(!propdata->skipprop);
   assert(SCIPgetNNlpis(scip) > 0);

   *result = SCIP_DIDNOTRUN;

   if( !isNlobbtApplicable(scip, propdata) )
   {
      /* do not call the propagator anymore (excepted after a restart) */
      SCIPdebugMsg(scip, "nlobbt propagator is not applicable\n");
      propdata->skipprop = TRUE;
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   nlpi = SCIPgetNlpis(scip)[0];
   SCIP_CALL( SCIPnlpiCreateProblem(nlpi, &nlpiprob, "nlobbt-nlp") );
   SCIP_CALL( SCIPhashmapCreate(&var2idx, SCIPblkmem(scip), SCIPcalcHashtableSize(SCIPgetNVars(scip))) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlcount, SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &status, SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &cands, SCIPgetVars(scip), SCIPgetNVars(scip)) );

   BMSclearMemoryArray(status, SCIPgetNVars(scip));

   /* compute convex NLP relaxation */
   SCIP_CALL( createNlpRelax(scip, nlpi, SCIPgetNLPNlRows(scip), SCIPgetNNLPNlRows(scip), nlpiprob, var2idx, nlcount) );
   SCIPnlpiSetRealPar(nlpi, nlpiprob, SCIP_NLPPAR_FEASTOL, SCIPfeastol(scip) * propdata->feastolfac);
   SCIPnlpiSetRealPar(nlpi, nlpiprob, SCIP_NLPPAR_RELOBJTOL, SCIPfeastol(scip) * propdata->relobjtolfac);

   /* add rows which have been found after creating the initial root LP */
   nlprows = SCIPgetNLPRows(scip) - SCIPgetNLPInitialRows(scip);
   if( nlprows > 0 )
   {
      SCIPdebugMsg(scip, "add %d (%d) LP rows\n", nlprows, SCIPgetNLPRows(scip));
      SCIP_CALL( addLpRows(scip, nlpi, nlpiprob, var2idx, &SCIPgetLPRows(scip)[SCIPgetNLPInitialRows(scip)], nlprows) );
   }

   /* use nlcount to decide on candidates */
   SCIPsortDownIntPtr(nlcount, (void*)cands, SCIPgetNVars(scip));

   /* main propagation loop */
   for( i = 0; i < SCIPgetNVars(scip) && nlcount[i] > 0 && *result != SCIP_CUTOFF && !SCIPisStopped(scip); ++i )
   {
      int varidx;

      /* skip binary or almost fixed variables */
      if( SCIPvarGetType(cands[i]) == SCIP_VARTYPE_BINARY
         || SCIPisRelEQ(scip, SCIPvarGetLbLocal(cands[i]), SCIPvarGetUbLocal(cands[i])) )
         continue;

      /* get index of var in the nlpi */
      assert(SCIPhashmapExists(var2idx, (void*)cands[i]) );
      varidx = (int)(size_t)SCIPhashmapGetImage(var2idx, (void*)cands[i]);
      assert(cands[i] == SCIPgetVars(scip)[varidx]);

      /* case: minimize var */
      if( status[i] != 1 && status[i] != 3 )
      {
         SCIP_CALL( solveNlp(scip, nlpi, nlpiprob, cands[i], varidx, SCIP_BOUNDTYPE_LOWER, var2idx, nlcount, status, cands, i, result) );
      }

      /* case: maximize var */
      if( *result != SCIP_CUTOFF && status[i] != 2 && status[i] != 3 )
      {
         SCIP_CALL( solveNlp(scip, nlpi, nlpiprob, cands[i], varidx, SCIP_BOUNDTYPE_UPPER, var2idx, nlcount, status, cands, i, result) );
      }
   }

   /* stop calling propagator if it did not find a reduction */
   if( *result == SCIP_DIDNOTFIND )
   {
      SCIPdebugMsg(scip, "no reductions found -> stop calling nlobbt propagator\n");
      propdata->skipprop = TRUE;
   }

   /* free the memory */
   SCIPfreeBufferArray(scip, &cands);
   SCIPfreeBufferArray(scip, &status);
   SCIPfreeBufferArray(scip, &nlcount);
   SCIPhashmapFree(&var2idx);
   SCIP_CALL( SCIPnlpiFreeProblem(nlpi, &nlpiprob) );

   return SCIP_OKAY;
}

/*
 * Callback methods of propagator
 */

/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeNlobbt)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIP_CALL( propdataFree(scip, &propdata) );
   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}

/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
static
SCIP_DECL_PROPEXITSOL(propExitsolNlobbt)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIP_CALL( propdataClear(scip, propdata) );

   return SCIP_OKAY;
}

/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecNlobbt)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   if( propdata->skipprop || SCIPgetStage(scip) != SCIP_STAGE_SOLVING || SCIPinRepropagation(scip)
      || SCIPinProbing(scip) || SCIPinDive(scip) || !SCIPallowObjProp(scip) || SCIPgetNNlpis(scip) == 0 )
   {
      SCIPdebugMsg(scip, "skip nlobbt propagator\n");
      return SCIP_OKAY;
   }

   /* do not run if SCIP does not have constructed an NLP */
   if( !SCIPisNLPConstructed(scip) )
   {
      SCIPdebugMsg(scip, "NLP not constructed, skipping nlobbt\n");
      return SCIP_OKAY;
   }

   /* call main procedure of nonlinear OBBT propagator */
   SCIP_CALL( applyNlobbt(scip, propdata, result) );

   return SCIP_OKAY;
}

/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropNlobbt)
{  /*lint --e{715}*/
   *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/*
 * propagator specific interface methods
 */

/** creates the nlobbt propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropNlobbt(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   propdata = NULL;
   prop = NULL;

   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );
   assert(propdata != NULL);
   BMSclearMemory(propdata);

   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecNlobbt, propdata) );
   assert(prop != NULL);

   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeNlobbt) );
   SCIP_CALL( SCIPsetPropExitsol(scip, prop, propExitsolNlobbt) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropNlobbt) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/"PROP_NAME"/feastolfac",
         "factor for NLP feasibility tolerance",
         &propdata->feastolfac, TRUE, DEFAULT_FEASTOLFAC, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/"PROP_NAME"/relobjtolfac",
         "factor for NLP relative objective tolerance",
         &propdata->relobjtolfac, TRUE, DEFAULT_RELOBJTOLFAC, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/"PROP_NAME"/minnonconvexfrac",
         "(nconvex nlrows)/(nonconvex nlrows) threshold to apply propagator",
         &propdata->minnonconvexfrac, TRUE, DEFAULT_MINNONCONVEXFRAC, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/"PROP_NAME"/minlinearfrac",
         "minimum (nconvex nlrows)/(linear nlrows) threshold to apply propagator",
         &propdata->minlinearfrac, TRUE, DEFAULT_MINLINEARFRAC, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/"PROP_NAME"/addlprows",
         "should non-initial LP rows be used?",
         &propdata->addlprows, FALSE, DEFAULT_ADDLPROWS, NULL, NULL) );

   return SCIP_OKAY;
}

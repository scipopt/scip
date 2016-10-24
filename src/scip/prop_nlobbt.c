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
#include "scip/prop_genvbounds.h"
#include "nlpi/nlpi.h"

#define PROP_NAME              "nlobbt"
#define PROP_DESC              "propagator template"
#define PROP_PRIORITY          -1100000
#define PROP_FREQ                     0
#define PROP_DELAY                 TRUE
#define PROP_TIMING            SCIP_PROPTIMING_AFTERLPLOOP

#define DEFAULT_MINNONCONVEXFRAC   0.20      /**< default minimum (#convex nlrows)/(#nonconvex nlrows) threshold to apply propagator */
#define DEFAULT_MINLINEARFRAC      0.02      /**< default minimum (#convex nlrows)/(#linear nlrows) threshold to apply propagator */
#define DEFAULT_FEASTOLFAC         0.01      /**< default factor for NLP feasibility tolerance */
#define DEFAULT_RELOBJTOLFAC       0.01      /**< default factor for NLP relative objective tolerance */
#define DEFAULT_ADDLPROWS          TRUE      /**< should (non-initial) LP rows be used? */
#define DEFAULT_ITLIMITFACTOR       2.0      /**< multiple of root node LP iterations used as total LP iteration
                                              *   limit for nlobbt (<= 0: no limit ) */
#define DEFAULT_NLPITERLIMIT          0      /**< default iteration limit of NLP solver; 0 for no limit */
#define DEFAULT_NLPTIMELIMIT        0.0      /**< default time limit of NLP solver; 0.0 for no limit */
#define DEFAULT_NLPVERLEVEL           0      /**< verbosity level of NLP solver */
#define DEFAULT_RANDSEED             79      /**< initial random seed */

/*
 * Data structures
 */

/** status of bound candidates */
enum BoundStatus
{
   UNSOLVED = 1,                             /**< did not solve LB or UB problem */
   SOLVEDLB = 2,                             /**< solved LB problem */
   SOLVEDUB = 4,                             /**< solved UB problem */
   SOLVED   = SOLVEDLB | SOLVEDUB            /**< solved LB and UB problem */
};
typedef enum BoundStatus BOUNDSTATUS;

/** propagator data */
struct SCIP_PropData
{
   SCIP_NLPI*            nlpi;               /**< nlpi used to create the nlpi problem */
   SCIP_NLPIPROBLEM*     nlpiprob;           /**< nlpi problem representing the convex NLP relaxation */
   SCIP_HASHMAP*         var2nlpiidx;        /**< mapping between variables and nlpi indices */
   SCIP_VAR**            nlpivars;           /**< array containing all variables of the nlpi */
   int                   nlpinvars;          /**< total number of nlpi variables */
   SCIP_Real*            nlscore;            /**< score for each nonlinear variable */
   int*                  status;             /**< array containing a bound status for each candidate (type int* is
                                              *   necessary to use sort functions) */
   SCIP_PROP*            genvboundprop;      /**< genvbound propagator */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   SCIP_Bool             skipprop;           /**< should the propagator be skipped? */
   SCIP_Longint          lastnode;           /**< number of last node where obbt was performed */
   int                   lastidx;            /**< index to store the last undone and unfiltered bound */

   int                   nlpiterlimit;       /**< iteration limit of NLP solver; 0 for no limit */
   SCIP_Real             nlptimelimit;       /**< time limit of NLP solver; 0.0 for no limit */
   int                   nlpverblevel;       /**< verbosity level of NLP solver */
   SCIP_NLPSTATISTICS*   nlpstatistics;      /**< statistics from NLP solver */

   SCIP_Real             feastolfac;         /**< factor for NLP feasibility tolerance */
   SCIP_Real             relobjtolfac;       /**< factor for NLP relative objective tolerance */
   SCIP_Real             minnonconvexfrac;   /**< minimum (#convex nlrows)/(#nonconvex nlrows) threshold to apply propagator */
   SCIP_Real             minlinearfrac;      /**< minimum (#convex nlrows)/(#linear nlrows) threshold to apply propagator */
   SCIP_Bool             addlprows;          /**< should (non-initial) LP rows be used? */
   SCIP_Real             itlimitfactor;      /**< LP iteration limit for nlobbt will be this factor times total LP
                                              *   iterations in root node */
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

   if( propdata->nlpiprob != NULL )
   {
      assert(propdata->nlpi != NULL);

      SCIPfreeMemoryArray(scip, &propdata->status);
      SCIPfreeMemoryArray(scip, &propdata->nlscore);
      SCIPfreeMemoryArray(scip, &propdata->nlpivars);
      SCIPhashmapFree(&propdata->var2nlpiidx);
      SCIPnlpiFreeProblem(propdata->nlpi, &propdata->nlpiprob);

      propdata->nlpinvars = 0;
   }
   assert(propdata->nlpinvars == 0);

   propdata->skipprop = FALSE;
   propdata->lastidx = 0;
   propdata->lastnode = -1;

   return SCIP_OKAY;
}

/** adds linear rows to the convex NLP relaxation */
static
SCIP_RETCODE nlpRelaxAddRows(
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

   SCIPdebugMsg(scip, "call nlpRelaxAddRows() with %d rows\n", nrows);

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

/** creates a convex NLP relaxation and stores it in a given NLPI problem; the function computes for each variable which
 *  the number of non-linearly occurrences and stores it in the nlscore array; the values are used to sort all variables
 *  in applyNlobbt() and process the 'best' candidates first; to reduce performance variability we store the scores as
 *  floating points to perturb them later
 *
 *  @note function does not copy the objective
 *
 *  @note the first row corresponds always to the cutoff row (even if cutoffbound is SCIPinfinity(scip))
 **/
static
SCIP_RETCODE nlpRelaxCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLROW**          nlrows,             /**< nonlinear rows */
   int                   nnlrows,            /**< total number of nonlinear rows */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< empty nlpi problem */
   SCIP_HASHMAP*         var2idx,            /**< empty hash map to store mapping between variables and indices in nlpi
                                              *   problem */
   SCIP_Real*            nlscore,            /**< array to store the score of each nonlinear variable (NULL if not
                                              *   needed) */
   SCIP_Real             cutoffbound         /**< cutoff bound */
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

   SCIPdebugMsg(scip, "call nlpRelaxCreate() with cutoffbound %g\n", cutoffbound);

   if( nlscore != NULL )
   {
      BMSclearMemoryArray(nlscore, SCIPgetNVars(scip));
   }
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   nconss = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &exprtrees, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &exprvaridxs, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nquadelems, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lininds, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlininds, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &names, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lhss, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhss, nnlrows + 1) );

   SCIP_CALL( SCIPallocBufferArray(scip, &lbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varnames, nvars) );

   /* create a unique mapping between variables and {0,..,nvars-1} */
   for( i = 0; i < nvars; ++i )
   {
      assert(vars[i] != NULL);
      SCIP_CALL( SCIPhashmapInsert(var2idx, (void*)vars[i], (void*)(size_t)i) );

      lbs[i] = SCIPvarGetLbLocal(vars[i]);
      ubs[i] = SCIPvarGetUbLocal(vars[i]);
      varnames[i] = SCIPvarGetName(vars[i]);
  }

   /* add variables */
   SCIP_CALL( SCIPnlpiAddVars(nlpi, nlpiprob, nvars, lbs, ubs, varnames) );
   SCIPfreeBufferArray(scip, &varnames);
   SCIPfreeBufferArray(scip, &ubs);
   SCIPfreeBufferArray(scip, &lbs);

   /* add row for cutoff bound even if cutoffbound == SCIPinfinity() */
   lhss[0] = -SCIPinfinity(scip);
   rhss[0] = cutoffbound;
   names[0] = "objcutoff";
   lininds[0] = NULL;
   linvals[0] = NULL;
   nlininds[0] = 0;
   nquadelems[0] = 0;
   quadelems[0] = NULL;
   exprtrees[0] = NULL;
   exprvaridxs[0] = NULL;

   SCIP_CALL( SCIPallocBufferArray(scip, &lininds[0], nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals[0], nvars) );

   for( i = 0; i < nvars; ++i )
   {
      if( !SCIPisZero(scip, SCIPvarGetObj(vars[i])) )
      {
         linvals[nconss][nlininds[0]] = SCIPvarGetObj(vars[i]);
         lininds[nconss][nlininds[0]] = i;
         ++nlininds[0];
      }
   }
   ++nconss;

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

            /* update nlscore */
            if( nlscore != NULL )
            {
               ++nlscore[quadelems[nconss][k].idx1];
               if( quadelems[nconss][k].idx1 != quadelems[nconss][k].idx2 )
                  ++nlscore[quadelems[nconss][k].idx2];
            }
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

            /* update nlscore */
            if( nlscore != NULL )
               ++nlscore[exprvaridxs[nconss][k]];
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

      if( linvals[i] != NULL )
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

/** updates bounds of each variable and the cutoff row in the nlpiproblem */
static
SCIP_RETCODE nlpRelaxUpdate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_Real* lbs;
   SCIP_Real* ubs;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int* inds;
   int i;

   SCIPdebugMsg(scip, "call nlpRelaxUpdate()\n");

   /* update variable bounds */
   SCIP_CALL( SCIPallocBufferArray(scip, &lbs, propdata->nlpinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubs, propdata->nlpinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, propdata->nlpinvars) );

   for( i = 0; i < propdata->nlpinvars; ++i )
   {
      assert(propdata->nlpivars[i] != NULL);
      assert(SCIPhashmapExists(propdata->var2nlpiidx, (void*)propdata->nlpivars[i]));

      lbs[i] = SCIPvarGetLbLocal(propdata->nlpivars[i]);
      ubs[i] = SCIPvarGetUbLocal(propdata->nlpivars[i]);
      inds[i] = (int)(size_t)SCIPhashmapGetImage(propdata->var2nlpiidx, (void*)propdata->nlpivars[i]);
      assert(inds[i] >= 0 && inds[i] < propdata->nlpinvars);
   }

   SCIP_CALL( SCIPnlpiChgVarBounds(propdata->nlpi, propdata->nlpiprob, propdata->nlpinvars, inds, lbs, ubs) );

   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &ubs);
   SCIPfreeBufferArray(scip, &lbs);

   /* update cutoff row */
   lhs = -SCIPinfinity(scip);
   rhs = SCIPgetCutoffbound(scip);
   i = 0;

   SCIP_CALL( SCIPnlpiChgConsSides(propdata->nlpi, propdata->nlpiprob, 1, &i, &lhs, &rhs) );

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
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   int                   pos,                /**< current position in the propdata->nlpivars array */
   SCIP_NLPSOLSTAT       nlpsolstat          /**< NLP solution status */
   )
{
   SCIP_Real* primal;

   assert(pos >= 0 && pos < propdata->nlpinvars);

   /* at least a feasible solution is needed */
   if( nlpsolstat > SCIP_NLPSOLSTAT_FEASIBLE )
      return SCIP_OKAY;

   SCIP_CALL( SCIPnlpiGetSolution(propdata->nlpi, propdata->nlpiprob, &primal, NULL, NULL, NULL) );
   assert(primal != NULL);

   for( ; pos < propdata->nlpinvars; ++pos )
   {
      SCIP_VAR* var;
      SCIP_Real val;
      int varidx;

      /* only uninteresting variables left -> stop filtering */
      if( SCIPisLE(scip, propdata->nlscore[pos], 0.0) )
         break;

      var = propdata->nlpivars[pos];
      assert(var != NULL && SCIPhashmapExists(propdata->var2nlpiidx, (void*)var));

      varidx = (int)(size_t)SCIPhashmapGetImage(propdata->var2nlpiidx, (void*)var);
      assert(SCIPgetVars(scip)[varidx] == var);
      val = primal[varidx];

      if( (propdata->status[pos] & SOLVEDLB) == 0 && !SCIPisInfinity(scip, -val) && SCIPisFeasLE(scip, val, SCIPvarGetLbLocal(var)) )
      {
         SCIPdebugMsg(scip, "filter LB of %s in [%g,%g] with %g\n", SCIPvarGetName(var), SCIPvarGetLbLocal(var),
            SCIPvarGetUbLocal(var), val);
         propdata->status[pos] |= SOLVEDLB;
         assert((propdata->status[pos] & SOLVEDLB) != 0);
      }

      if( (propdata->status[pos] & SOLVEDUB) == 0 && !SCIPisInfinity(scip, val) && SCIPisFeasGE(scip, val, SCIPvarGetUbLocal(var)) )
      {
         SCIPdebugMsg(scip, "filter UB of %s in [%g,%g] with %g\n", SCIPvarGetName(var), SCIPvarGetLbLocal(var),
            SCIPvarGetUbLocal(var), val);
         propdata->status[pos] |= SOLVEDUB;
         assert((propdata->status[pos] & SOLVEDUB) != 0);
      }
   }

   return SCIP_OKAY;
}

/** tries to add a generalized variable bound by exploiting the dual solution of the last NLP solve */
static
SCIP_RETCODE addGenVBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var,                /**< variable used in last NLP solve */
   int                   varidx,             /**< variable index in the propdata->nlpivars array */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound provided by the genvbound */
   SCIP_Real             cutoffbound         /**< cutoff bound */
   )
{
   SCIP_VAR** lvbvars;
   SCIP_Real* lvbcoefs;
   SCIP_Real* primal;
   SCIP_Real* dual;
   SCIP_Real* alpha;
   SCIP_Real* beta;
   SCIP_Real constant;
   SCIP_Real gamma;
   int nlvbvars;
   int i;

   assert(propdata->genvboundprop != NULL);
   assert(var != NULL);
   assert(varidx >= 0 && varidx < propdata->nlpinvars);

   if( SCIPnlpiGetSolstat(propdata->nlpi, propdata->nlpiprob) > SCIP_NLPSOLSTAT_LOCOPT )
      return SCIP_OKAY;

   SCIP_CALL( SCIPnlpiGetSolution(propdata->nlpi, propdata->nlpiprob, &primal, &dual, &alpha, &beta) );

   /* not possible to generate genvbound if the duals for the propagated variable do not disappear */
   if( !SCIPisFeasZero(scip, alpha[varidx] + beta[varidx]) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &lvbcoefs, propdata->nlpinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lvbvars, propdata->nlpinvars) );
   constant = boundtype == SCIP_BOUNDTYPE_LOWER ? primal[varidx] : -primal[varidx];
   gamma = 0.0;
   nlvbvars = 0;

   /* collect coefficients of genvbound */
   for( i = 0; i < propdata->nlpinvars; ++i )
   {
      if( !SCIPisZero(scip, beta[i] - alpha[i]) )
      {
         lvbvars[nlvbvars] = propdata->nlpivars[i];
         lvbcoefs[nlvbvars] = beta[i] - alpha[i];
         ++nlvbvars;

         constant += (alpha[i] - beta[i]) * primal[i];
      }
   }

   /* first dual multiplier corresponds to the cutoff row if cutoffbound < SCIPinfinity() */
   if( !SCIPisInfinity(scip, cutoffbound) && SCIPisGT(scip, dual[0], 0.0) )
   {
      gamma = dual[0];
      constant += gamma * cutoffbound;
   }

   /* add genvbound to genvbounds propagator */
   if( !SCIPisInfinity(scip, REALABS(constant)) && (nlvbvars > 0 || SCIPisFeasGT(scip, gamma, 0.0)) )
   {
      SCIP_CALL( SCIPgenVBoundAdd(scip, propdata->genvboundprop, lvbvars, var, lvbcoefs, nlvbvars, -gamma, constant,
            boundtype) );
      SCIPdebugMsg(scip, "add genvbound for %s\n", SCIPvarGetName(var));
   }

   SCIPfreeBufferArray(scip, &lvbvars);
   SCIPfreeBufferArray(scip, &lvbcoefs);

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
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var,                /**< variable to propagate */
   int                   varidx,             /**< variable index in the propdata->nlpivars array */
   SCIP_BOUNDTYPE        boundtype,          /**< minimize or maximize var? */
   int                   pos,                /**< current position in propdata->nlpivars array */
   int*                  nlpiter,            /**< buffer to store the total number of nlp iterations */
   SCIP_RESULT*          result              /**< pointer to store result */
   )
{
   SCIP_Real timelimit;
   SCIP_Real* primal;
   SCIP_Real obj;
   int iterlimit;

#ifdef SCIP_DEBUG
   SCIP_Real oldlb;
   SCIP_Real oldub;

   oldlb = SCIPvarGetLbLocal(var);
   oldub = SCIPvarGetUbLocal(var);
#endif

   assert(var != NULL);
   assert(varidx >= 0 && varidx < propdata->nlpinvars);
   assert(result != NULL && *result != SCIP_CUTOFF);

   *nlpiter = 0;

   /* set time and iteration limit */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
   {
      timelimit -= SCIPgetSolvingTime(scip);
      if( timelimit <= 0.0 )
      {
         SCIPdebugMsg(scip, "skip NLP solve; no time left\n");
         return SCIP_OKAY;
      }
   }
   if( propdata->nlptimelimit > 0.0 )
      timelimit = MIN(propdata->nlptimelimit, timelimit);
   iterlimit = propdata->nlpiterlimit > 0 ? propdata->nlpiterlimit : INT_MAX;
   SCIP_CALL( SCIPnlpiSetRealPar(propdata->nlpi, propdata->nlpiprob, SCIP_NLPPAR_TILIM, timelimit) );
   SCIP_CALL( SCIPnlpiSetIntPar(propdata->nlpi, propdata->nlpiprob, SCIP_NLPPAR_ITLIM, iterlimit) );

   /* set corresponding objective coefficient and solve NLP */
   obj = boundtype == SCIP_BOUNDTYPE_LOWER ? 1.0 : -1.0;
   SCIP_CALL( SCIPnlpiSetObjective(propdata->nlpi, propdata->nlpiprob, 1, &varidx, &obj, 0, NULL, NULL, NULL, 0.0) );

   SCIPdebugMsg(scip, "solve var=%s boundtype=%d nlscore=%g\n", SCIPvarGetName(var), boundtype,
      propdata->nlscore[pos]);
   SCIP_CALL( SCIPnlpiSolve(propdata->nlpi, propdata->nlpiprob) );
   SCIPdebugMsg(scip, "NLP solstat = %d\n", SCIPnlpiGetSolstat(propdata->nlpi, propdata->nlpiprob));

   /* collect NLP statistics */
   assert(propdata->nlpstatistics != NULL);
   SCIP_CALL( SCIPnlpiGetStatistics(propdata->nlpi, propdata->nlpiprob, propdata->nlpstatistics) );
   *nlpiter = SCIPnlpStatisticsGetNIterations(propdata->nlpstatistics);
   SCIPdebugMsg(scip, "iterations %d time %g\n", *nlpiter, SCIPnlpStatisticsGetTotalTime(propdata->nlpstatistics));

   /* try to add a genvbound in the root node */
   if( propdata->genvboundprop != NULL && SCIPgetDepth(scip) == 0 )
   {
      SCIP_CALL( addGenVBound(scip, propdata, var, varidx, boundtype, SCIPgetCutoffbound(scip)) );
   }

   /* filter bound candidates first, otherwise we do not have access to the primal solution values */
   SCIP_CALL( filterCands(scip, propdata, pos+1, SCIPnlpiGetSolstat(propdata->nlpi, propdata->nlpiprob)) );

   /* try to tighten variable bound */
   if( SCIPnlpiGetSolstat(propdata->nlpi, propdata->nlpiprob) <= SCIP_NLPSOLSTAT_LOCOPT )
   {
      SCIP_Bool tightened;
      SCIP_Bool infeasible;

      SCIP_CALL( SCIPnlpiGetSolution(propdata->nlpi, propdata->nlpiprob, &primal, NULL, NULL, NULL) );

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
         SCIP_CALL( SCIPnlpiChgVarBounds(propdata->nlpi, propdata->nlpiprob, 1, &varidx, &lb, &ub) );

#ifdef SCIP_DEBUG
         SCIPdebugMsg(scip, "tightened bounds of %s from [%g,%g] to [%g,%g]\n", SCIPvarGetName(var), oldlb, oldub, lb, ub);
#endif
      }
   }

   /* reset objective function */
   obj = 0.0;
   SCIP_CALL( SCIPnlpiSetObjective(propdata->nlpi, propdata->nlpiprob, 1, &varidx, &obj, 0, NULL, NULL, NULL, 0.0) );

   return SCIP_OKAY;
}

/** main method of the propagator
 *
 *  creates a convex NLP relaxation and solves the OBBT-NLPs for each possible candidate;
 *  binary and variables with a small domain will be ignored to reduce the computational cost of the propagator; after
 *  solving each NLP we filter out all variable candidates which are on their lower or upper bound; candidates with a
 *  larger number of occurrences are preferred
 */
static
SCIP_RETCODE applyNlobbt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagation data */
   SCIP_RESULT*          result              /**< pointer to store result */
   )
{
   int nlpiterleft;
   int i;

   assert(result != NULL);
   assert(!propdata->skipprop);
   assert(SCIPgetNNlpis(scip) > 0);

   *result = SCIP_DIDNOTRUN;

   if( propdata->nlpiprob == NULL && !isNlobbtApplicable(scip, propdata) )
   {
      /* do not call the propagator anymore (except after a restart) */
      SCIPdebugMsg(scip, "nlobbt propagator is not applicable\n");
      propdata->skipprop = TRUE;
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   /* compute NLP iteration limit */
   if( propdata->itlimitfactor > 0.0 )
      nlpiterleft = propdata->itlimitfactor * SCIPgetNRootLPIterations(scip);
   else
      nlpiterleft = INT_MAX;

   /* recompute NLP relaxation if the variable set changed */
   if( propdata->nlpiprob != NULL && SCIPgetNVars(scip) != propdata->nlpinvars )
   {
      SCIP_CALL( propdataClear(scip, propdata) );
      assert(propdata->nlpiprob == NULL);
   }

   /* create or update NLP relaxation */
   if( propdata->nlpiprob == NULL )
   {
      propdata->nlpinvars = SCIPgetNVars(scip);
      propdata->nlpi = SCIPgetNlpis(scip)[0];
      assert(propdata->nlpi != NULL);

      SCIP_CALL( SCIPnlpiCreateProblem(propdata->nlpi, &propdata->nlpiprob, "nlobbt-nlp") );
      SCIP_CALL( SCIPhashmapCreate(&propdata->var2nlpiidx, SCIPblkmem(scip),
            SCIPcalcHashtableSize(propdata->nlpinvars)) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &propdata->nlpivars, SCIPgetVars(scip), propdata->nlpinvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &propdata->nlscore, propdata->nlpinvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &propdata->status, propdata->nlpinvars) );

      SCIP_CALL( nlpRelaxCreate(scip, propdata->nlpi, SCIPgetNLPNlRows(scip), SCIPgetNNLPNlRows(scip), propdata->nlpiprob,
            propdata->var2nlpiidx, propdata->nlscore, SCIPgetCutoffbound(scip)) );

      /* initialize bound status; perturb nlscores by a factor which ensures sure that zero scores remain zero */
      assert(propdata->randnumgen != NULL);
      for( i = 0; i < propdata->nlpinvars; ++i )
      {
         propdata->status[i] = UNSOLVED;
         propdata->nlscore[i] *= 1.0 + SCIPrandomGetReal(propdata->randnumgen, SCIPfeastol(scip), 2.0 * SCIPfeastol(scip));
      }

      /* add rows of the LP */
      if( SCIPgetDepth(scip) == 0 )
      {
         SCIP_CALL( nlpRelaxAddRows(scip, propdata->nlpi, propdata->nlpiprob, propdata->var2nlpiidx, SCIPgetLPRows(scip),
               SCIPgetNLPRows(scip)) );
      }
   }
   else
   {
      SCIP_CALL( nlpRelaxUpdate(scip, propdata) );
   }

   assert(propdata->nlpiprob != NULL);
   assert(propdata->var2nlpiidx != NULL);
   assert(propdata->nlpivars != NULL);
   assert(propdata->nlscore != NULL);

   /* sort variables w.r.t. their nlscores if we did not solve any NLP for this node */
   if( propdata->lastidx == 0 )
   {
      SCIPsortDownRealIntPtr(propdata->nlscore, propdata->status, (void*)propdata->nlpivars, propdata->nlpinvars);
   }

   /* set parameters of NLP solver */
   SCIP_CALL( SCIPnlpiSetRealPar(propdata->nlpi, propdata->nlpiprob, SCIP_NLPPAR_FEASTOL,
         SCIPfeastol(scip) * propdata->feastolfac) );
   SCIP_CALL( SCIPnlpiSetRealPar(propdata->nlpi, propdata->nlpiprob, SCIP_NLPPAR_FEASTOL,
         SCIPfeastol(scip) * propdata->feastolfac) );
   SCIP_CALL( SCIPnlpiSetRealPar(propdata->nlpi, propdata->nlpiprob, SCIP_NLPPAR_RELOBJTOL,
         SCIPfeastol(scip) * propdata->relobjtolfac) );
   SCIP_CALL( SCIPnlpiSetIntPar(propdata->nlpi, propdata->nlpiprob, SCIP_NLPPAR_VERBLEVEL, propdata->nlpverblevel) );

   /* main propagation loop */
   for( i = propdata->lastidx; i < propdata->nlpinvars
           && nlpiterleft > 0
           && SCIPisGT(scip, propdata->nlscore[i], 0.0)
           && *result != SCIP_CUTOFF
           && !SCIPisStopped(scip);
        ++i )
   {
      SCIP_VAR* var;
      int varidx;
      int iters;

      var = propdata->nlpivars[i];
      assert(var != NULL);

      /* skip binary or almost fixed variables */
      if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY
         || SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         continue;

      SCIPdebugMsg(scip, "iterations left %d\n", nlpiterleft);

      /* get index of var in the nlpi */
      assert(SCIPhashmapExists(propdata->var2nlpiidx, (void*)var) );
      varidx = (int)(size_t)SCIPhashmapGetImage(propdata->var2nlpiidx, (void*)var);
      assert(var == SCIPgetVars(scip)[varidx]);

      /* case: minimize var */
      if( (propdata->status[i] & SOLVEDLB) == 0 )
      {
         SCIP_CALL( solveNlp(scip, propdata, var, varidx, SCIP_BOUNDTYPE_LOWER, i, &iters, result) );
         nlpiterleft -= iters;
      }

      /* case: maximize var */
      if( *result != SCIP_CUTOFF && (propdata->status[i] & SOLVEDUB) == 0 )
      {
         SCIP_CALL( solveNlp(scip, propdata, var, varidx, SCIP_BOUNDTYPE_UPPER, i, &iters, result) );
         nlpiterleft -= iters;
      }
   }

   /* remember last position */
   propdata->lastidx = i;

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

   SCIP_CALL( propdataClear(scip, propdata) );
   SCIPfreeBlockMemory(scip, &propdata);
   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}

/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
static
SCIP_DECL_PROPINITSOL(propInitsolNlobbt)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert(scip != NULL);
   assert(prop != NULL);

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* if genvbounds propagator is not available, we cannot create genvbounds */
   propdata->genvboundprop = SCIPfindProp(scip, "genvbounds");

   SCIP_CALL( SCIPrandomCreate(&propdata->randnumgen, SCIPblkmem(scip),
         SCIPinitializeRandomSeed(scip, DEFAULT_RANDSEED)) );
   SCIP_CALL( SCIPnlpStatisticsCreate(&propdata->nlpstatistics) );
   propdata->lastnode = -1;

   return SCIP_OKAY;
}

/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
static
SCIP_DECL_PROPEXITSOL(propExitsolNlobbt)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIPnlpStatisticsFree(&propdata->nlpstatistics);
   SCIPrandomFree(&propdata->randnumgen);

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

   /* only run if LP all columns are in the LP, i.e., the LP is a relaxation; e.g., do not run if pricers are active
    * since pricing is not performed in probing mode
    */
   if( !SCIPallColsInLP(scip) )
   {
      SCIPdebugMsg(scip, "not all columns in LP, skipping obbt\n");
      return SCIP_OKAY;
   }

   /* do not run if SCIP does not have constructed an NLP */
   if( !SCIPisNLPConstructed(scip) )
   {
      SCIPdebugMsg(scip, "NLP not constructed, skipping nlobbt\n");
      return SCIP_OKAY;
   }

   if( SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) != propdata->lastnode )
   {
      propdata->lastnode = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
      propdata->lastidx = 0;
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
   SCIP_CALL( SCIPsetPropInitsol(scip, prop, propInitsolNlobbt) );
   SCIP_CALL( SCIPsetPropExitsol(scip, prop, propExitsolNlobbt) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropNlobbt) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/"PROP_NAME"/feastolfac",
         "factor for NLP feasibility tolerance",
         &propdata->feastolfac, TRUE, DEFAULT_FEASTOLFAC, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/"PROP_NAME"/relobjtolfac",
         "factor for NLP relative objective tolerance",
         &propdata->relobjtolfac, TRUE, DEFAULT_RELOBJTOLFAC, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/"PROP_NAME"/minnonconvexfrac",
         "(#convex nlrows)/(#nonconvex nlrows) threshold to apply propagator",
         &propdata->minnonconvexfrac, TRUE, DEFAULT_MINNONCONVEXFRAC, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/"PROP_NAME"/minlinearfrac",
         "minimum (#convex nlrows)/(#linear nlrows) threshold to apply propagator",
         &propdata->minlinearfrac, TRUE, DEFAULT_MINLINEARFRAC, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/"PROP_NAME"/addlprows",
         "should non-initial LP rows be used?",
         &propdata->addlprows, FALSE, DEFAULT_ADDLPROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "propagating/"PROP_NAME"/nlpiterlimit",
         "iteration limit of NLP solver; 0 for no limit",
         &propdata->nlpiterlimit, TRUE, DEFAULT_NLPTIMELIMIT, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/"PROP_NAME"/nlptimelimit",
         "time limit of NLP solver; 0.0 for no limit",
         &propdata->nlptimelimit, TRUE, DEFAULT_NLPTIMELIMIT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "propagating/"PROP_NAME"/nlpverblevel",
         "verbosity level of NLP solver",
         &propdata->nlpverblevel, TRUE, DEFAULT_NLPVERLEVEL, 0, 5, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/"PROP_NAME"/itlimitfactor",
         "LP iteration limit for nlobbt will be this factor times total LP iterations in root node",
         &propdata->itlimitfactor, TRUE, DEFAULT_ITLIMITFACTOR, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

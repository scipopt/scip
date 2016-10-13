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

/*
 * Data structures
 */


/** bound data */
struct Bound
{
   SCIP_VAR*             var;                /**< variable */
   SCIP_Real             newval;             /**< stores a probably tighter value for this bound */
   SCIP_BOUNDTYPE        boundtype;          /**< type of bound */
   unsigned int          found:1;            /**< stores whether a probably tighter value for this bound was found */
   unsigned int          done:1;             /**< has this bound been processed already? */
};
typedef struct Bound BOUND;

/** propagator data */
struct SCIP_PropData
{
   BOUND**               bounds;             /**< array of interesting bounds */
   int                   nbounds;            /**< number of interesting bounds */
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
   assert(scip != NULL);
   assert(propdata != NULL);

   /* free the bounds candidates */
   if( propdata->bounds != NULL )
   {
      int i;

      for( i = propdata->nbounds - 1; i >= 0; --i )
      {
         BOUND* bound;
         bound = propdata->bounds[i];
         SCIPfreeMemory(scip, &bound);
      }

      SCIPfreeMemoryArray(scip, &propdata->bounds);
   }

   propdata->nbounds = 0;

   return SCIP_OKAY;
}

/** clears and frees propagator data */
static
SCIP_RETCODE propdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA**       propdata            /**< propagator data */
   )
{
   assert(scip != NULL);
   assert(propdata != NULL);
   assert(*propdata != NULL);

   SCIP_CALL( propdataClear(scip, *propdata) );
   SCIPfreeBlockMemory(scip, propdata);
   *propdata = NULL;

   return SCIP_OKAY;
}

/** creates a convex NLP relaxation and stores it in a given NLPI problem */
static
SCIP_RETCODE createNlpRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLROW**          nlrows,             /**< nonlinear rows */
   int                   nnlrows,            /**< total number of nonlinear rows */
   SCIP_NLPIPROBLEM**    nlpiprob,           /**< pointer to store the nlpi problem */
   SCIP_HASHMAP*         var2idx             /**< empty hash map to store mapping between variables and indices in nlpi
                                              *   problem */
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
   SCIP_Real* objvals;
   SCIP_Real* lbs;
   SCIP_Real* ubs;
   const char** varnames;
   int* objinds;
   int nconss;
   int i;

   assert(scip != NULL);
   assert(nlpiprob != NULL);
   assert(var2idx != NULL);
   assert(nlrows != NULL);
   assert(nnlrows > 0);
   assert(nlpi != NULL);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   nconss = 0;

   SCIP_CALL( SCIPnlpiCreateProblem(nlpi, nlpiprob, "convex_NLP") );

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

   SCIP_CALL( SCIPallocBufferArray(scip, &objvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &objinds, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varnames, nvars) );

   /* create an unique mapping from variables to {0,..,nvars-1} */
   for( i = 0; i < nvars; ++i )
   {
      assert(vars[i] != NULL);
      SCIP_CALL( SCIPhashmapInsert(var2idx, (void*)vars[i], (void*)(size_t)i) );

      objinds[i] = i;
      objvals[i] = SCIPvarGetObj(vars[i]);

      lbs[i] = SCIPvarGetLbLocal(vars[i]);
      ubs[i] = SCIPvarGetUbLocal(vars[i]);
      varnames[i] = SCIPvarGetName(vars[i]);
  }

   /** add variables */
   SCIP_CALL( SCIPnlpiAddVars(nlpi, *nlpiprob, nvars, lbs, ubs, varnames) );

   /* set objective */
   SCIP_CALL( SCIPnlpiSetObjective(nlpi, *nlpiprob, nvars, objinds, objvals, 0, NULL, NULL, NULL, 0.0) );
   SCIPfreeBufferArray(scip, &objinds);
   SCIPfreeBufferArray(scip, &objvals);

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
         }
      }

      ++nconss;
   }
   assert(nconss > 0);

   SCIP_CALL( SCIPnlpiAddConstraints(nlpi, *nlpiprob, nconss, lhss, rhss, nlininds, lininds, linvals, nquadelems,
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

   return SCIP_OKAY;
}

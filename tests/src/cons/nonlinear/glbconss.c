/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
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

/**@file   glbconss.c
 * @brief  test for generating globally valid nonlinear constraints
 * @author Benjamin Mueller
 *
 * This test adds a constraint handler to separate globally valid nonlinear constraints during the tree search. For a
 * given set of variables \f$\{(x_i,y_j)\}\f$, the constraint handler enforces the following quadratic constraints
 *
 * \f[
 *     (x_i - x_j)^2 + (y_i - y_j)^2 \ge (r_i + r_j)^2,
 * \f]
 *
 * for every pair \f$i,j\f$ with \f$i < j\f$. The quadratic constraints are motivated by non-overlapping conditions in
 * the circle packing problem. During the enforcement, the constraint handler checks all of the quadratic constraints
 * and adds the first one that is violated and has not been separated so far.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "cpp"
#define CONSHDLR_DESC          "constraint handler for circle packing problem structure"
#define CONSHDLR_ENFOPRIORITY       -50 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */
#define CONSHDLR_MAXPREROUNDS         0 /**< maximal rounds in presolve; set to 1 to enable upgrade to nonlinear cons in presolve */

/*
 * forward declarations
 */

/*
 * constraint specific interface methods
 */
static
SCIP_RETCODE SCIPincludeConshdlrCpp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a circle packing constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
static
SCIP_RETCODE SCIPcreateConsCpp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR**            xs,                 /**< variables representing the x-coordinate of each circle */
   SCIP_VAR**            ys,                 /**< variables representing the y-coordinate of each circle */
   SCIP_Real*            rs,                 /**< radius of each circle */
   int                   ncircles,           /**< total number of circles */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

/** creates and captures a circle packing constraint in its most basic form */
static
SCIP_RETCODE SCIPcreateConsBasicCpp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR**            xs,                 /**< variables representing the x-coordinate of each circle */
   SCIP_VAR**            ys,                 /**< variables representing the y-coordinate of each circle */
   SCIP_Real*            rs,                 /**< radius of each circle */
   int                   ncircles            /**< total number of circles */
   );

/*
 * Data structures
 */

/** constraint data for cpp constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            xs;                 /**< variables representing the x-coordinate of each circle */
   SCIP_VAR**            ys;                 /**< variables representing the y-coordinate of each circle */
   SCIP_Real*            rs;                 /**< radius of each circle */
   int                   ncircles;           /**< total number of circles */
   SCIP_HASHMAP*         map;                /**< hashmap to store which nonlinear constraints have been generated */
};

/*
 * Local methods
 */

/** creates constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to constraint data */
   SCIP_VAR**            xs,                 /**< variables representing the x-coordinate of each circle */
   SCIP_VAR**            ys,                 /**< variables representing the y-coordinate of each circle */
   SCIP_Real*            rs,                 /**< radius of each circle */
   int                   ncircles            /**< total number of circles */
   )
{
   int i;

   assert(consdata != NULL);
   assert(xs != NULL);
   assert(ys != NULL);
   assert(rs != NULL);
   assert(ncircles > 0);

   SCIP_CALL( SCIPallocClearBlockMemory(scip, consdata) );

   /* duplicate data */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->xs, xs, ncircles) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->ys, ys, ncircles) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->rs, rs, ncircles) );
   (*consdata)->ncircles = ncircles;

   /* transform variables, if necessary */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->ncircles, (*consdata)->xs, (*consdata)->xs) );
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->ncircles, (*consdata)->ys, (*consdata)->ys) );
   }

   /* capture variables */
   for( i = 0; i < ncircles; ++i )
   {
      SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->xs[i]) );
      SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->ys[i]) );
   }

   return SCIP_OKAY;
}

/** frees constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to constraint data */
   )
{
   int i;

   /* release variables */
   for( i = 0; i < (*consdata)->ncircles; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->xs[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->ys[i]) );
   }

   /* free memory */
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->xs, (*consdata)->ncircles);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->ys, (*consdata)->ncircles);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->rs, (*consdata)->ncircles);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** helper method to check whether a solution satisfies a circle packing constraint */
static
SCIP_RESULT checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to check (NULL for the LP solution) */
   SCIP_CONS*            cons,               /**< circle packing constraint */
   int*                  violi,              /**< index of the first circle that violates the constraint (might be NULL) */
   int*                  violj               /**< index of the second circle that violates the constraint (might be NULL) */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   for( i = 0; i < consdata->ncircles - 1; ++i )
   {
      SCIP_Real xi;
      SCIP_Real yi;
      SCIP_Real ri;
      int j;

      xi = SCIPgetSolVal(scip, sol, consdata->xs[i]);
      yi = SCIPgetSolVal(scip, sol, consdata->ys[i]);
      ri = consdata->rs[i];

      for( j = i+1; j < consdata->ncircles; ++j )
      {
         SCIP_Real xj;
         SCIP_Real yj;
         SCIP_Real rj;

         xj = SCIPgetSolVal(scip, sol, consdata->xs[j]);
         yj = SCIPgetSolVal(scip, sol, consdata->ys[j]);
         rj = consdata->rs[j];

         if( SCIPisFeasLT(scip, xi*xi + xj*xj + yi*yi + yj*yj -2.0*xi*xj -2.0*yi*yj, SQR(ri+rj)) )
         {
            /* store the index of the circles that violate the constraint */
            if( violi != NULL )
               *violi = i;
            if( violj != NULL )
               *violj = j;

            SCIPdebugMsg(scip, "solution is infeasible for constraint %s (%d,%d)\n", SCIPconsGetName(cons), i, j);
            return SCIP_INFEASIBLE;
         }
      }
   }

   return SCIP_FEASIBLE;
}

/** helper method to enforce a given set of constraints
 *
 * Note that the circle packing constraint handler has a smaller enforce priority than the nonlinear constraint handler.
 */
static
SCIP_RETCODE enfoConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_CONS**           conss,              /**< circle packing constraints */
   int                   nconss,             /**< total number of circle packing constraints */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   int c;

   assert(conss != NULL);
   assert(nconss > 0);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;
      int violi;
      int violj;

      assert(conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      assert(consdata->map != NULL);

      /* only consider violated constraints */
      if( checkCons(scip, sol, conss[c], &violi, &violj) == SCIP_INFEASIBLE )
      {
         SCIP_VAR* xi;
         SCIP_VAR* yi;
         SCIP_VAR* xj;
         SCIP_VAR* yj;
         SCIP_Real ri;
         SCIP_Real rj;
         int hash;

         *result = SCIP_INFEASIBLE;

         assert(violi >= 0 && violi < consdata->ncircles);
         assert(violj >= 0 && violj < consdata->ncircles);

         xi = consdata->xs[violi];
         yi = consdata->ys[violi];
         xj = consdata->xs[violj];
         yj = consdata->ys[violj];
         ri = consdata->rs[violi];
         rj = consdata->rs[violj];
         hash = violi * consdata->ncircles + violj;

         /* make sure that every nonlinear constraint is only added once */
         if( !SCIPhashmapExists(consdata->map, (void*)(size_t)hash) )
         {
            SCIP_VAR* quadvars1[6] = {xi, xj, xi, yi, yj, yi};
            SCIP_VAR* quadvars2[6] = {xi, xj, xj, yi, yj, yj};
            SCIP_Real quadcoefs[6] = {1.0, 1.0, -2.0, 1.0, 1.0, -2.0};
            SCIP_Real lhs = SQR(ri + rj);
            SCIP_CONS* cons;
            char name[SCIP_MAXSTRLEN];

            SCIP_CALL( SCIPhashmapInsert(consdata->map, (void*)(size_t)hash, NULL) );

            /* create a nonlinear constraint */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d_%d", SCIPconsGetName(conss[c]), violi, violj);
            SCIP_CALL( SCIPcreateConsQuadraticNonlinear(scip, &cons, name, 0, NULL, NULL, 6, quadvars1, quadvars2,
               quadcoefs, lhs, SCIPinfinity(scip),
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );

            SCIPdebugMsg(scip, "add constraint %s in node %lld\n", name, SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

            /* add and release constraint */
            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );

            *result = SCIP_CONSADDED;

            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteCpp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* free constraint data */
   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransCpp)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* create constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->xs, sourcedata->ys, sourcedata->rs, sourcedata->ncircles) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpCpp)
{  /*lint --e{715}*/

   SCIP_CALL( enfoConss(scip, NULL, conss, nconss, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsCpp)
{  /*lint --e{715}*/

   SCIP_CALL( enfoConss(scip, NULL, conss, nconss, result) );

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckCpp)
{  /*lint --e{715}*/
   int c;

   *result = SCIP_FEASIBLE;

   for( c = 0; c < nconss && (*result == SCIP_FEASIBLE); ++c )
   {
      *result = checkCons(scip, sol, conss[c], NULL, NULL);
   }

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockCpp)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   for( i = 0; i < consdata->ncircles; ++i )
   {
      SCIP_CALL( SCIPaddVarLocksType(scip, consdata->xs[i], locktype, nlocksneg + nlockspos, nlocksneg + nlockspos) );
      SCIP_CALL( SCIPaddVarLocksType(scip, consdata->ys[i], locktype, nlocksneg + nlockspos, nlocksneg + nlockspos) );
   }

   return SCIP_OKAY;
}

/** presolve that replaces cpp constraint by all nonlinear constraints
 * (for debug purposes)
 */
static
SCIP_DECL_CONSPRESOL(consPresolCpp)
{
   int c;

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;
      int i, j;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      assert(consdata->map == NULL);

      for( i = 0; i < consdata->ncircles; ++i )
      {
         for( j = i+1; j < consdata->ncircles; ++j )
         {
            SCIP_VAR* xi;
            SCIP_VAR* yi;
            SCIP_VAR* xj;
            SCIP_VAR* yj;
            SCIP_Real ri;
            SCIP_Real rj;

            xi = consdata->xs[i];
            yi = consdata->ys[i];
            xj = consdata->xs[j];
            yj = consdata->ys[j];
            ri = consdata->rs[i];
            rj = consdata->rs[j];

            SCIP_VAR* quadvars1[6] = {xi, xj, xi, yi, yj, yi};
            SCIP_VAR* quadvars2[6] = {xi, xj, xj, yi, yj, yj};
            SCIP_Real quadcoefs[6] = {1.0, 1.0, -2.0, 1.0, 1.0, -2.0};
            SCIP_Real lhs = SQR(ri + rj);
            SCIP_CONS* cons;
            char name[SCIP_MAXSTRLEN];

            /* create a nonlinear constraint */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d_%d", SCIPconsGetName(conss[c]), i, j);
            SCIP_CALL( SCIPcreateConsQuadraticNonlinear(scip, &cons, name, 0, NULL, NULL, 6, quadvars1, quadvars2,
               quadcoefs, lhs, SCIPinfinity(scip),
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );

            SCIPdebugMsg(scip, "add constraint %s in node %lld\n", name, SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

            /* add and release constraint */
            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }
      }

      SCIP_CALL( SCIPdelCons(scip, conss[c]) );
      ++*nupgdconss;
   }

   return SCIP_OKAY;
}

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolCpp)
{  /*lint --e{715}*/
   int c;

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      assert(consdata->map == NULL);

      SCIP_CALL( SCIPhashmapCreate(&consdata->map, SCIPblkmem(scip), SQR(consdata->ncircles)) );
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolCpp)
{  /*lint --e{715}*/
   int c;

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      assert(consdata->map != NULL);

      SCIPhashmapFree(&consdata->map);
   }

   return SCIP_OKAY;

}

/*
 * constraint specific interface methods
 */
static
SCIP_RETCODE SCIPincludeConshdlrCpp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{

   SCIP_CONSHDLR* conshdlr;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpCpp, consEnfopsCpp, consCheckCpp, consLockCpp,
         NULL) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransCpp) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteCpp) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolCpp, CONSHDLR_MAXPREROUNDS, SCIP_PRESOLTIMING_ALWAYS) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolCpp) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolCpp) );

   return SCIP_OKAY;
}

/** creates and captures a circle packing constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
static
SCIP_RETCODE SCIPcreateConsCpp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR**            xs,                 /**< variables representing the x-coordinate of each circle */
   SCIP_VAR**            ys,                 /**< variables representing the y-coordinate of each circle */
   SCIP_Real*            rs,                 /**< radius of each circle */
   int                   ncircles,           /**< total number of circles */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the cpp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("cpp constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, xs, ys, rs, ncircles) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a circle packing constraint in its most basic form */
static
SCIP_RETCODE SCIPcreateConsBasicCpp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR**            xs,                 /**< variables representing the x-coordinate of each circle */
   SCIP_VAR**            ys,                 /**< variables representing the y-coordinate of each circle */
   SCIP_Real*            rs,                 /**< radius of each circle */
   int                   ncircles            /**< total number of circles */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcreateConsCpp(scip, cons, name, xs, ys, rs, ncircles,
      TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/*
 * TESTS
 */

static SCIP* scip;

/** setup of test run */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include circle packing constraint handler */
   SCIP_CALL( SCIPincludeConshdlrCpp(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );
}

/** deinitialization method */
static
void teardown(void)
{
   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

/** helper method to create x- and y-variables */
static
SCIP_RETCODE createVars(
   SCIP_Real             lbx,                /**< lower bound of each x-variable */
   SCIP_Real             ubx,                /**< upper bound of each x-variable */
   SCIP_Real             objx,               /**< objective coefficient of each x-variable */
   SCIP_Real             lby,                /**< lower bound of each y-variable */
   SCIP_Real             uby,                /**< upper bound of each y-variable */
   SCIP_Real             objy,               /**< objective coefficient of each y-variable */
   int                   nvars,              /**< total number of x- and y-variables */
   SCIP_VAR**            xs,                 /**< array to store x-variables */
   SCIP_VAR**            ys                  /**< array to store y-variables */
   )
{
   char name[SCIP_MAXSTRLEN];
   int i;

   for( i = 0; i < nvars; ++i )
   {
      (void) SCIPsnprintf(name , SCIP_MAXSTRLEN, "x_%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &xs[i], name, lbx, ubx, objx, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, xs[i]) );

      (void) SCIPsnprintf(name , SCIP_MAXSTRLEN, "y_%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &ys[i], name, lby, uby, objy, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, ys[i]) );
   }

   return SCIP_OKAY;
}

/* define test suite */
TestSuite(glbconss, .init = setup, .fini = teardown);

/*
 * define tests
 */

Test(glbconss, cpp1)
{
   SCIP_CONS* cons;
   SCIP_VAR* xs[4];
   SCIP_VAR* ys[4];
   SCIP_Real rs[4] = {1.0, 1.0, 1.0, 1.0};
   SCIP_SOL* sol;
   int ncircles = 4;
   int i;

   /* create variables */
   SCIP_CALL( createVars(0.0, 10.0, -1.0, 0.0, 2.0, 0.0, ncircles, xs, ys) );

   /* create circle packing constraint */
   SCIP_CALL( SCIPcreateConsBasicCpp(scip, &cons, "cpp", xs, ys, rs, ncircles) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* set a node limit so it has a chance to fail instead of running forever
    * currently (2020-06-19) solves in 1620 nodes
    */
   SCIP_CALL( SCIPsetLongintParam(scip, "limits/totalnodes", 10000L) );

   /* solve problem */
   SCIP_CALL( SCIPsolve(scip) );
   cr_expect(SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL);
   cr_expect_float_eq(SCIPgetPrimalbound(scip), -3.60000003185305e+01, 1e-4);
   sol = SCIPgetBestSol(scip);
   cr_assert(sol != NULL);

   /* check whether all circles satisfy the distance requirements */
   for( i = 0; i < ncircles - 1; ++i )
   {
      SCIP_Real xi = SCIPgetSolVal(scip, sol, xs[i]);
      SCIP_Real yi = SCIPgetSolVal(scip, sol, ys[i]);
      int j;

      for( j = i + 1; j < ncircles; ++j )
      {
         SCIP_Real xj = SCIPgetSolVal(scip, sol, xs[j]);
         SCIP_Real yj = SCIPgetSolVal(scip, sol, ys[j]);
         SCIP_Real dist;

         dist = SQR(xi - xj) + SQR(yi - yj);
         cr_expect(SCIPisFeasGE(scip, dist, SQR(rs[i] + rs[j])));
      }
   }

   /* release variables */
   for( i = 0; i < ncircles; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &xs[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &ys[i]) );
   }
}

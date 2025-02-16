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

/**@file   cons_orbitope.c
 * @ingroup DEFPLUGINS_CONS
 * @brief  interface for constraint handlers of type partitioning, packing, and full to ensure backwards compatibility
 * @author Christopher Hojny
 *
 * This interface ensures backwards compatibility to be able to add packing, partitioning, and full
 * orbitopes via the same function call.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/cons_orbitope.h"
#include "scip/cons_orbitope_full.h"
#include "scip/cons_orbitope_pp.h"
#include "scip/symmetry.h"
#include <symmetry/type_symmetry.h>

/*
 * Local methods
 */

/** strengthen full orbitopes to packing/partitioning orbitopes if possible */
static
SCIP_RETCODE strengthenOrbitopeConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< variable matrix of orbitope constraint */
   int*                  nrows,              /**< pointer to number of rows of variable matrix */
   int                   ncols,              /**< number of columns of variable matrix */
   SCIP_ORBITOPETYPE*    type                /**< pointer to store type of orbitope constraint after strengthening */
   )
{
   SCIP_Bool* pprows = NULL;
   int npprows;
   int nrowsorig;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( vars != NULL );
   assert( *nrows > 0 );
   assert( ncols > 0 );
   assert( type != NULL );

   nrowsorig = *nrows;
   SCIP_CALL( SCIPisPackingPartitioningOrbitope(scip, vars, *nrows, ncols, &pprows, &npprows, type) );

   /* If only some rows are contained in set packing/partitioning constraints, it may still be worth it
    * to exploit the packing/partitioning structure on these rows, because packing/partitioning orbitopes
    * are more restrictive than full orbitopes. If at least three rows have this property, we discard
    * all rows not contained in set packing/partitioning constraints and add the smaller packing sub-orbitope.
    */
   if ( npprows >= 3 )
   {
      int r = *nrows - 1;
      int i;

      assert( pprows != NULL );

      while ( r >= 0 )
      {
         if ( ! pprows[r] )
         {
            for (i = r; i < *nrows - 1; ++i)
            {
               SCIP_VAR** varrow;
               varrow = vars[i];
               vars[i] = vars[i+1];
               vars[i+1] = varrow;
            }
            --(*nrows);
         }
         --r;
      }
      *type = SCIP_ORBITOPETYPE_PACKING;
   }

   /* pprows might not have been initialized if there are no setppc conss */
   if ( pprows != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &pprows, nrowsorig);
   }

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** includes the orbitope constraint handlers used by this interface */
SCIP_RETCODE SCIPincludeConshdlrOrbitope(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPincludeConshdlrOrbitopeFull(scip) );
   SCIP_CALL( SCIPincludeConshdlrOrbitopePP(scip) );

   return SCIP_OKAY;
}


/** creates and captures an orbitope constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR***           vars,               /**< matrix of variables on which the symmetry acts */
   SCIP_ORBITOPETYPE     orbitopetype,       /**< type of orbitope constraint */
   int                   nrows,              /**< number of rows of variable matrix */
   int                   ncols,              /**< number of columns of variable matrix */
   SCIP_Bool             resolveprop,        /**< should propagation be resolved? */
   SCIP_Bool             ismodelcons,        /**< whether the orbitope is a model constraint */
   SCIP_Bool             checkpporbitope,    /**< Check if full orbitope constraints can be upgraded to pp-orbitope? */
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
   assert( nrows > 0 );
   assert( ncols > 0 );

   if ( checkpporbitope && orbitopetype != SCIP_ORBITOPETYPE_PARTITIONING
      && orbitopetype != SCIP_ORBITOPETYPE_PACKING )
   {
      SCIP_CALL( strengthenOrbitopeConstraint(scip, vars, &nrows, ncols, &orbitopetype) );
   }

   if ( orbitopetype == SCIP_ORBITOPETYPE_FULL )
   {
      SCIP_CALL( SCIPcreateConsOrbitopeFull(scip, cons, name, vars, nrows, ncols, resolveprop, ismodelcons,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }
   else
   {
      SCIP_CALL( SCIPcreateConsOrbitopePP(scip, cons, name, vars, orbitopetype, nrows, ncols,
            resolveprop, ismodelcons, initial, separate, enforce, check, propagate,
            local, modifiable, dynamic, removable, stickingatnode) );
   }

   return SCIP_OKAY;
}

/** creates and captures an orbitope constraint in its most basic variant, i. e., with all constraint flags set to their
 *  default values, which can be set afterwards using SCIPsetConsFLAGNAME()
 *
 *  @see SCIPcreateConsOrbitope() for the default constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR***           vars,               /**< matrix of variables on which the symmetry acts */
   SCIP_ORBITOPETYPE     orbitopetype,       /**< type of orbitope constraint */
   int                   nrows,              /**< number of rows of variable matrix */
   int                   ncols,              /**< number of columns of variable matrix */
   SCIP_Bool             resolveprop,        /**< should propagation be resolved? */
   SCIP_Bool             ismodelcons,        /**< whether the orbitope is a model constraint */
   SCIP_Bool             checkpporbitope     /**< Check if full orbitope constraints can be upgraded to pp-orbitope? */
   )
{
   SCIP_CALL( SCIPcreateConsOrbitope(scip, cons, name, vars, orbitopetype, nrows, ncols,
         resolveprop, ismodelcons, checkpporbitope,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

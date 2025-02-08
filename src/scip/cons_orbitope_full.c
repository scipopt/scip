/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

/**@file   cons_orbitope_full.c
 * @ingroup DEFPLUGINS_CONS
 * @brief  constraint handler for full orbitope constraints w.r.t. the full symmetric group
 * @author Timo Berthold
 * @author Marc Pfetsch
 * @author Christopher Hojny
 *
 * The type of constraints of this constraint handler is described in cons_orbitope_full.h.
 * When creating the constraint, users can decide whether it is a constraint defining the model
 * or "just" use to handle symmetries. In the latter case, symmetry reductions are only performed
 * by the constraint handler if strong dual reductions are permitted.
 *
 * The details of the method implemented here are described in the following papers.
 *
 * Orbitopal fixing for the full (sub-)orbitope and application to the Unit Commitment Problem@n
 * Pascale Bendotti, Pierre Fouilhoux, and Cecile Rottner,@n
 * Optimization Online: http://www.optimization-online.org/DB_HTML/2017/10/6301.html
 *
 * Two linear time propagation algorithms for full orbitopes are described in this paper, a static
 * version and a dynamic one. While the static version uses a fixed variable order, the dynamic
 * version determines the variable order during the solving process via branching descisions.
 * We only implemented the static version, because constraints should define the model and should
 * not be changed during the solving process. Instead, a dynamic version of orbitopal fixing has
 * been implemented as a routine in prop_symmetry.c.
 *
 * Polytopes associated with symmetry handling@n
 * Christopher Hojny and Marc E. Pfetsch,@n
 * Math. Program. (2018)
 *
 * In this paper, a linear time separation algorithm for orbisacks (full orbitopes with two columnes)
 * is described. We use this algorithm for every pair of adjacent columns within the orbitope.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/cons_orbisack.h"
#include "scip/cons_orbitope_full.h"
#include "scip/cons_setppc.h"
#include "scip/pub_cons.h"
#include "scip/pub_message.h"
#include "scip/pub_var.h"
#include "scip/scip.h"
#include "scip/scip_branch.h"
#include "scip/scip_conflict.h"
#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_cut.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_sol.h"
#include "scip/scip_var.h"
#include "scip/symmetry.h"
#include <symmetry/type_symmetry.h>

/* constraint handler properties */
#define CONSHDLR_NAME          "orbitope_full"
#define CONSHDLR_DESC          "symmetry breaking constraint handler relying on full orbitopes"
#define CONSHDLR_SEPAPRIORITY    +40100 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -1005200 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -1005200 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler */
#define CONSHDLR_PRESOLTIMING            SCIP_PRESOLTIMING_MEDIUM /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */

#define DEFAULT_FORCECONSCOPY     FALSE /**< whether orbitope constraints should be forced to be copied to sub SCIPs */

/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Bool             forceconscopy;      /**< whether orbitope constraints should be forced to be copied to sub SCIPs */
};

/** constraint data for orbitope constraints */
struct SCIP_ConsData
{
   SCIP_VAR***           vars;               /**< matrix of variables on which the symmetry acts */
   int                   nrows;              /**< number of rows */
   int                   ncols;              /**< number of columns*/
   SCIP_Bool             resolveprop;        /**< should propagation be resolved? */
   SCIP_Bool             ismodelcons;        /**< whether the orbitope is a model constraint */
};


/*
 * Local methods
 */

/** frees an orbitope constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to orbitope constraint data */
   )
{
   int i;
   int j;
   int p;
   int q;

   assert( consdata != NULL );
   assert( *consdata != NULL );

   p = (*consdata)->nrows;
   q = (*consdata)->ncols;
   for (i = 0; i < p; ++i)
   {
      /* release variables in vars array */
      for (j = 0; j < q; ++j)
      {
         assert( (*consdata)->vars[i] != NULL );
         SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->vars[i][j]) );
      }

      SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vars[i]), q);     /*lint !e866*/
   }

   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vars), p);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** creates orbitope constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store constraint data */
   SCIP_VAR***           vars,               /**< variables array, must have size nrows x ncols*/
   int                   nrows,              /**< number of rows */
   int                   ncols,              /**< number of columns */
   SCIP_Bool             resolveprop,        /**< should propagation be resolved? */
   SCIP_Bool             ismodelcons         /**< whether the orbitope is a model constraint */
   )
{
   int i;
   int j;

   assert(consdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->vars, nrows) );

   for (i = 0; i < nrows; ++i)
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars[i], vars[i], ncols) );    /*lint !e866*/
   }

   (*consdata)->nrows = nrows;
   (*consdata)->ncols = ncols;
   (*consdata)->resolveprop = resolveprop;
   (*consdata)->ismodelcons = ismodelcons;

   /* get transformed variables, if we are in the transformed problem */
   if ( SCIPisTransformed(scip) )
   {
      for (i = 0; i < nrows; ++i)
      {
         /* make sure that no variable gets multiaggregated (cannot be handled by cons_orbitope, since one cannot easily
          * eliminate single variables from an orbitope constraint).
          */
         for (j = 0; j < ncols; ++j)
         {
            SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->vars[i][j], &(*consdata)->vars[i][j]) );
            SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, (*consdata)->vars[i][j]) );
         }
      }
   }

   /* capture vars contained in vars array */
   for (i = 0; i < nrows; ++i)
   {
      for (j = 0; j < ncols; ++j)
      {
         assert( (*consdata)->vars[i][j] != NULL );
         SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->vars[i][j]) );
      }
   }

   return SCIP_OKAY;
}


/** Compute lexicographically minimal face of the hypercube w.r.t. some coordinate fixing */
static
SCIP_RETCODE findLexMinFace(
   SCIP_VAR***           vars,               /**< variable matrix */
   int**                 lexminfixes,        /**< fixings characterzing lex-min face */
   int*                  minfixedrowlexmin,  /**< index of minimum fixed row for each column or NULL (if in prop) */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility has been detected or NULL (if in resprop) */
   int                   m,                  /**< number of rows in vars */
   int                   n,                  /**< number of columns in vars */
   SCIP_Bool             resprop             /**< whether we are in resprop (TRUE) or prop (FALSE) */
   )
{
   int i;
   int j;

   *infeasible = FALSE;

   assert( vars != NULL );
   assert( lexminfixes != NULL );
   assert( !resprop || minfixedrowlexmin != NULL );
   assert( m > 0 );
   assert( n > 0 );
   assert( infeasible != NULL );

   /* iterate over columns in reverse order and find the lexicographically minimal face
    * of the hypercube containing lexminfixes
    */
   for (j = n - 2; j >= 0; --j)
   {
      int maxdiscriminating = m;
      int minfixed = -1;

      /* fix free entries in column j to the corresponding value in column j + 1 and collect some information */
      for (i = 0; i < m; ++i)
      {
         /* is row i j-discriminating? */
         if ( minfixed == -1 && lexminfixes[i][j] != 0 && lexminfixes[i][j + 1] != 1 )
         {
            assert( lexminfixes[i][j + 1] == 0 );

            maxdiscriminating = i;
         }

         /* is row i j-fixed? */
         if ( minfixed == -1 && lexminfixes[i][j] != lexminfixes[i][j + 1] && lexminfixes[i][j] != 2 )
         {
            assert( lexminfixes[i][j + 1] != 2 );

            minfixed = i;

            /* detect infeasibility */
            if ( maxdiscriminating > minfixed )
            {
               *infeasible = TRUE;

               return SCIP_OKAY;
            }
         }
      }

      /* ensure that column j is lexicographically not smaller than column j + 1 */
      for (i = 0; i < m; ++i)
      {
         if ( lexminfixes[i][j] == 2 )
         {
            if ( i < maxdiscriminating || minfixed == -1 )
               lexminfixes[i][j] = lexminfixes[i][j + 1];
            else if ( i == maxdiscriminating )
               lexminfixes[i][j] = 1;
            else
               lexminfixes[i][j] = 0;
         }
      }

      if ( resprop )
      {
         assert( minfixedrowlexmin != NULL );

         /* store minimum fixed row */
         if ( minfixed == -1 )
            minfixedrowlexmin[j] = m - 1;
         else
            minfixedrowlexmin[j] = minfixed;

         /* columns 1, ..., n-2 are contained in two columns (take the minimum) and
          * the minimum fixed row of column n-1 is determined by column n-2 */
         if ( minfixedrowlexmin[j + 1] < minfixedrowlexmin[j] )
            minfixedrowlexmin[j + 1] = minfixedrowlexmin[j];
      }
   }

   return SCIP_OKAY;
}


/** Compute lexicographically maximal face of the hypercube w.r.t. some coordinate fixing */
static
SCIP_RETCODE findLexMaxFace(
   SCIP_VAR***           vars,               /**< variable matrix */
   int**                 lexmaxfixes,        /**< fixings characterzing lex-max face */
   int*                  minfixedrowlexmax,  /**< index of minimum fixed row for each column or NULL (if in prop) */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility has been detected or NULL (if in resprop) */
   int                   m,                  /**< number of rows in vars */
   int                   n,                  /**< number of columns in vars */
   SCIP_Bool             resprop             /**< whether we are in resprop (TRUE) or prop (FALSE) */
   )
{
   int i;
   int j;

   *infeasible = FALSE;

   assert( vars != NULL );
   assert( lexmaxfixes != NULL );
   assert( !resprop || minfixedrowlexmax != NULL );
   assert( m > 0 );
   assert( n > 0 );
   assert( infeasible != NULL );

   for (j = 1; j < n; ++j)
   {
      int maxdiscriminating = m;
      int minfixed = -1;

      /* fix free entries in column j to the corresponding value in column j - 1 and collect some information */
      for (i = 0; i < m; ++i)
      {
         /* is row i j-discriminating? */
         if ( minfixed == -1 && lexmaxfixes[i][j - 1] != 0 && lexmaxfixes[i][j] != 1 )
         {
            assert( lexmaxfixes[i][j - 1] == 1 );

            maxdiscriminating = i;
         }

         /* is row i j-fixed? */
         if ( minfixed == -1 && lexmaxfixes[i][j - 1] != lexmaxfixes[i][j] && lexmaxfixes[i][j] != 2 )
         {
            assert( lexmaxfixes[i][j - 1] != 2 );

            minfixed = i;

            /* detect infeasibility */
            if ( maxdiscriminating > minfixed )
            {
               *infeasible = TRUE;

               return SCIP_OKAY;
            }
         }
      }

      /* ensure that column j is lexicographically not greater than column j - 1 */
      for (i = 0; i < m; ++i)
      {
         if ( lexmaxfixes[i][j] == 2 )
         {
            if ( i < maxdiscriminating || minfixed == -1 )
               lexmaxfixes[i][j] = lexmaxfixes[i][j - 1];
            else if ( i == maxdiscriminating )
               lexmaxfixes[i][j] = 0;
            else
               lexmaxfixes[i][j] = 1;
         }
      }

      if ( resprop )
      {
         assert( minfixedrowlexmax != NULL );

         /* store minimum fixed row */
         if ( minfixed == -1 )
            minfixedrowlexmax[j] = m - 1;
         else
            minfixedrowlexmax[j] = minfixed;

         /* columns 1, ..., n-2 are contained in two columns (take the minimum) and
          * the minimum fixed row of column 0 is determined by column 1 */
         if ( minfixedrowlexmax[j - 1] < minfixedrowlexmax[j] )
            minfixedrowlexmax[j - 1] = minfixedrowlexmax[j];
      }
   }

   return SCIP_OKAY;
}


/** propagation method for a single full orbitope constraint */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be processed */
   SCIP_Bool*            infeasible,         /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars          /**< pointer to add up the number of found domain reductions */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int** lexminfixes;
   int** lexmaxfixes;
   int i;
   int j;
   int m;
   int n;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );
   assert( nfixedvars != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   *nfixedvars = 0;
   *infeasible = FALSE;

   /* if the constraint is not a model constraint, check whether symmetry reductions are permitted */
   if ( !consdata->ismodelcons && !SCIPallowStrongDualReds(scip) )
      return SCIP_OKAY;

   assert( consdata->nrows > 0 );
   assert( consdata->ncols > 0 );
   assert( consdata->vars != NULL );

   m = consdata->nrows;
   n = consdata->ncols;
   vars = consdata->vars;

   /* Initialize lexicographically minimal matrix by fixed entries at the current node.
    * Free entries in the last column are set to 0.
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &lexminfixes, m) );
   for (i = 0; i < m; ++i)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &lexminfixes[i], n) ); /*lint !e866*/
   }

   for (i = 0; i < m; ++i)
   {
      for (j = 0; j < n; ++j)
      {
         if ( SCIPvarGetLbLocal(vars[i][j]) > 0.5 )
            lexminfixes[i][j] = 1;
         else if ( SCIPvarGetUbLocal(vars[i][j]) < 0.5 || j == n - 1 )
            lexminfixes[i][j] = 0;
         else
            lexminfixes[i][j] = 2;
      }
   }

   /* find lexicographically minimal face of hypercube containing lexmin fixes */
   SCIP_CALL( findLexMinFace(vars, lexminfixes, NULL, infeasible, m, n, FALSE) );

   if ( *infeasible == TRUE )
      goto FREELEXMIN;

   /* Initialize lexicographically maximal matrix by fixed entries at the current node.
    * Free entries in the first column are set to 1.
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &lexmaxfixes, m) );
   for (i = 0; i < m; ++i)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &lexmaxfixes[i], n) ); /*lint !e866*/
   }

   for (i = 0; i < m; ++i)
   {
      for (j = 0; j < n; ++j)
      {
         if ( SCIPvarGetUbLocal(vars[i][j]) < 0.5 )
            lexmaxfixes[i][j] = 0;
         else if ( SCIPvarGetLbLocal(vars[i][j]) > 0.5 || j == 0 )
            lexmaxfixes[i][j] = 1;
         else
            lexmaxfixes[i][j] = 2;
      }
   }

   /* find lexicographically maximal face of hypercube containing lexmax fixes */
   SCIP_CALL( findLexMaxFace(vars, lexmaxfixes, NULL, infeasible, m, n, FALSE) );

   if ( *infeasible )
      goto FREELEXMAX;

   /* Find for each column j the minimal row in which lexminfixes and lexmaxfixes differ. Fix all entries above this
    * row to the corresponding value in lexminfixes (or lexmaxfixes).
    */
   for (j = 0; j < n; ++j)
   {
      for (i = 0; i < m; ++i)
      {
         if ( lexminfixes[i][j] != lexmaxfixes[i][j] )
            break;

         if ( SCIPvarGetLbLocal(vars[i][j]) < 0.5 && SCIPvarGetUbLocal(vars[i][j]) > 0.5 )
         {
            SCIP_Bool success;

            SCIP_CALL( SCIPinferBinvarCons(scip, vars[i][j], (SCIP_Bool) lexminfixes[i][j],
                  cons, 0, infeasible, &success) );

            if ( success )
               *nfixedvars += 1;
         }
      }
   }

 FREELEXMAX:
   for (i = 0; i < m; ++i)
      SCIPfreeBufferArray(scip, &lexmaxfixes[i]);
   SCIPfreeBufferArray(scip, &lexmaxfixes);

 FREELEXMIN:
   for (i = 0; i < m; ++i)
      SCIPfreeBufferArray(scip, &lexminfixes[i]);
   SCIPfreeBufferArray(scip, &lexminfixes);

   return SCIP_OKAY;
}


/** Propagation conflict resolving method of propagator
 *
 *  In this function we use that all variable reductions that can be found by the propagation algorithm
 *  are only due to the fixed variables that are in or above the minimum fixed row of each pair of adjacent
 *  columns of the lexmin and lexmax matrices.
 *
 *  Since the storage of an integer is not enough to store the complete information about the fixing,
 *  we have to use the linear time algorithm for finding the lexmin and lexmax
 *  matrices and determine from this the minimum fixed rows.
 */
static
SCIP_RETCODE resolvePropagationFullOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler of the corresponding constraint */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   int                   inferinfo,          /**< inference information */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation conflict resolving call */
   )
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int** lexminfixes;
   int** lexmaxfixes;
   int* minfixedrowlexmin;
   int* minfixedrowlexmax;
   int i;
   int j;
   int m;
   int n;
   SCIP_Bool terminate;

   *result = SCIP_DIDNOTFIND;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( result != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->nrows > 0 );
   assert( consdata->ncols > 0 );
   assert( consdata->vars != NULL );

   m = consdata->nrows;
   n = consdata->ncols;
   vars = consdata->vars;

   assert( inferinfo <= consdata->nrows );

   /* Initialize lexicographically minimal matrix by fixed entries at the current node.
    * Free entries in the last column are set to 0.
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &lexminfixes, m) );
   for (i = 0; i < m; ++i)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &lexminfixes[i], n) ); /*lint !e866*/
   }

   /* store minimum fixed row for each column */
   SCIP_CALL( SCIPallocBufferArray(scip, &minfixedrowlexmin, n) );
   minfixedrowlexmin[n - 1] = -1;

   for (i = 0; i < m; ++i)
   {
      for (j = 0; j < n; ++j)
      {
         if ( SCIPvarGetLbAtIndex(vars[i][j], bdchgidx, FALSE) > 0.5 )
            lexminfixes[i][j] = 1;
         else if ( SCIPvarGetUbAtIndex(vars[i][j], bdchgidx, FALSE) < 0.5 || j == n - 1 )
            lexminfixes[i][j] = 0;
         else
            lexminfixes[i][j] = 2;
      }
   }

   /* find lexicographically minimal face of hypercube containing lexmin fixes */
   SCIP_CALL( findLexMinFace(vars, lexminfixes, minfixedrowlexmin, &terminate, m, n, TRUE) );

   if ( terminate )
      goto FREELEXMIN;

   /* Initialize lexicographically maximal matrix by fixed entries at the current node.
    * Free entries in the first column are set to 1.
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &lexmaxfixes, m) );
   for (i = 0; i < m; ++i)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &lexmaxfixes[i], n) ); /*lint !e866*/
   }

   /* store minimum fixed row for each column */
   SCIP_CALL( SCIPallocBufferArray(scip, &minfixedrowlexmax, n) );
   minfixedrowlexmax[0] = -1;

   for (i = 0; i < m; ++i)
   {
      for (j = 0; j < n; ++j)
      {
         if ( SCIPvarGetUbAtIndex(vars[i][j], bdchgidx, FALSE) < 0.5 )
            lexmaxfixes[i][j] = 0;
         else if ( SCIPvarGetLbAtIndex(vars[i][j], bdchgidx, FALSE) > 0.5 || j == 0 )
            lexmaxfixes[i][j] = 1;
         else
            lexmaxfixes[i][j] = 2;
      }
   }

   /* find lexicographically maximal face of hypercube containing lexmax fixes */
   SCIP_CALL( findLexMaxFace(vars, lexmaxfixes, minfixedrowlexmax, &terminate, m, n, TRUE) );

   if ( terminate )
      goto FREELEXMAX;

   /* Find for each column j the minimal row in which lexminfixes and lexmaxfixes differ. Fix all entries above this
    * row to the corresponding value in lexminfixes (or lexmaxfixes).
    */
   for (j = 0; j < n; ++j)
   {
      int ub = MAX(minfixedrowlexmin[j], minfixedrowlexmax[j]);

      for (i = 0; i <= ub; ++i)
      {
         if ( SCIPvarGetLbAtIndex(vars[i][j], bdchgidx, FALSE) > 0.5 ||
            SCIPvarGetUbAtIndex(vars[i][j], bdchgidx, FALSE) < 0.5 )
         {
            SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i][j]) );
            *result = SCIP_SUCCESS;
         }
      }
   }

 FREELEXMAX:
   SCIPfreeBufferArray(scip, &minfixedrowlexmax);
   for (i = 0; i < m; ++i)
      SCIPfreeBufferArray(scip, &lexmaxfixes[i]);
   SCIPfreeBufferArray(scip, &lexmaxfixes);

 FREELEXMIN:
   SCIPfreeBufferArray(scip, &minfixedrowlexmin);
   for (i = 0; i < m; ++i)
      SCIPfreeBufferArray(scip, &lexminfixes[i]);
   SCIPfreeBufferArray(scip, &lexminfixes);

   return SCIP_OKAY;
}


/** check full orbitope solution for feasibility */
static
SCIP_RETCODE checkFullOrbitopeSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to process */
   SCIP_SOL*             sol,                /**< solution to be checked */
   SCIP_Bool             printreason,        /**< whether reason for infeasibility should be printed */
   SCIP_Bool*            feasible            /**< memory address to store whether solution is feasible */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   int nrows;
   int ncols;
   int j;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( feasible != NULL );

   consdata = SCIPconsGetData(cons);

   assert( consdata != NULL );
   assert( consdata->vars != NULL );
   assert( consdata->nrows > 0 );
   assert( consdata->ncols > 0 );
   assert( ! consdata->ismodelcons ); /* non-model constraints are never checked */

   vars = consdata->vars;
   nrows = consdata->nrows;
   ncols = consdata->ncols;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars1, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars2, nrows) );

   /* iterate over adjacent columns of orbitope and check whether the first column in this
    * column pair is lexicographically not smaller than the second column in the pair */
   *feasible = TRUE;
   for (j = 1; j < ncols && *feasible; ++j)
   {
      for (i = 0; i < nrows; ++i)
      {
         vars1[i] = vars[i][j - 1];
         vars2[i] = vars[i][j];
      }

      SCIP_CALL( SCIPcheckSolutionOrbisack(scip, sol, vars1, vars2, nrows, printreason, feasible) );
   }

   SCIPfreeBufferArray(scip, &vars2);
   SCIPfreeBufferArray(scip, &vars1);

   return SCIP_OKAY;
}


/** separate orbisack cover inequalities */
static
SCIP_RETCODE separateCoversOrbisack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to process */
   SCIP_SOL*             sol,                /**< solution to separate (NULL for the LP solution) */
   int*                  ngen,               /**< pointer to store number of generated cuts */
   SCIP_Bool*            infeasible          /**< pointer to store whether infeasibility has been detected */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int nrows;
   int ncols;
   int i;
   int j;
   SCIP_Real rhs;
   SCIP_Real lhs;
   SCIP_Real* coeffs1;
   SCIP_Real* coeffs2;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( ngen != NULL );
   assert( infeasible != NULL );

   *ngen = 0;
   *infeasible = FALSE;

   /* get basic data */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   vars = consdata->vars;
   nrows = consdata->nrows;
   ncols = consdata->ncols;

   /* allocate memory for cover inequalities */
   SCIP_CALL( SCIPallocBufferArray(scip, &coeffs1, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coeffs2, nrows) );

   lhs = 0.0;
   rhs = 0.0;

   /* separate orbisack cover inequalities for adjacent columns */
   for (j = 0; j < ncols - 1 && ! *infeasible; ++j)
   {
      SCIP_Real rowval;

      for (i = 0; i < nrows; ++i)
      {
         rowval = SCIPgetSolVal(scip, sol, vars[i][j + 1]) - SCIPgetSolVal(scip, sol, vars[i][j]);

         /* check whether cover inequality is violated */
         if ( SCIPisEfficacious(scip, rowval + lhs - rhs) )
         {
            SCIP_ROW* row;
            int k;

            /* set coefficients for current inequality */
            coeffs1[i] = -1.0;
            coeffs2[i] = 1.0;

            /* add violated orbisack cover inequality */
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, "orbisackcover", -SCIPinfinity(scip), rhs,
                  FALSE, FALSE, TRUE) );
            SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

            for (k = 0; k <= i; ++k)
            {
               SCIP_CALL( SCIPaddVarToRow(scip, row, vars[k][j], coeffs1[k]) );
               SCIP_CALL( SCIPaddVarToRow(scip, row, vars[k][j + 1], coeffs2[k]) );
            }
            SCIP_CALL( SCIPflushRowExtensions(scip, row) );

            SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
#ifdef SCIP_DEBUG
            SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif
            SCIP_CALL( SCIPreleaseRow(scip, &row) );

            *ngen += 1;
            if ( *infeasible )
               break;

            /* reset coefficients for next inequality */
            coeffs1[i] = 0.0;
            coeffs2[i] = 0.0;
         }

         /* add argmax( 1 - vals[i][0], vals[i][1] ) as coefficient and ensure that both vars1[0] and vars2[0] are
          * contained in the LIFTED cover inequality */
         rowval = SCIPgetSolVal(scip, sol, vars[i][j]) + SCIPgetSolVal(scip, sol, vars[i][j + 1]);
         if ( SCIPisEfficacious(scip, 1.0 - rowval) )
         {
            coeffs1[i] = -1.0;
            coeffs2[i] = 0.0;
            lhs -= SCIPgetSolVal(scip, sol, vars[i][j]);

            /* apply lifting? */
            if ( i == 0 )
            {
               coeffs2[i] = 1.0;
               lhs += SCIPgetSolVal(scip, sol, vars[i][j + 1]);
            }
         }
         else
         {
            coeffs1[i] = 0.0;
            coeffs2[i] = 1.0;
            lhs += SCIPgetSolVal(scip, sol, vars[i][j]);
            rhs += 1.0;

            /* apply lifting? */
            if ( i == 0 )
            {
               coeffs1[i] = -1.0;
               lhs -= SCIPgetSolVal(scip, sol, vars[i][j]);
               rhs -= 1.0;
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &coeffs1);
   SCIPfreeBufferArray(scip, &coeffs2);

   return SCIP_OKAY;
}


/** separate or enforce constraints */
static
SCIP_RETCODE separateConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   int                   nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
   SCIP_SOL*             sol,                /**< solution to separate (NULL for the LP solution) */
   SCIP_RESULT*          result,             /**< pointer to store the result (should be initialized) */
   SCIP_Bool             enforce             /**< whether we enforce orbitope constraints */
   )
{
   SCIP_Bool infeasible = FALSE;
   int nfixedvars = 0;
   int ncuts = 0;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   /* loop through constraints */
   for (c = 0; c < nconss && ! infeasible; c++)
   {
      SCIP_CONSDATA* consdata;
      int nconscuts = 0;

      assert( conss[c] != NULL );

      /* get data of constraint */
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* skip non-model constraints if strong dual reductions are not permitted */
      if ( !consdata->ismodelcons && !SCIPallowStrongDualReds(scip) )
         continue;

      /* do not enforce non-model constraints */
      if ( enforce && !consdata->ismodelcons )
         continue;

      SCIP_CALL( separateCoversOrbisack(scip, conss[c], sol, &nconscuts, &infeasible) );
      ncuts += nconscuts;

      /* stop after the useful constraints if we found cuts of fixed variables */
      if ( c >= nusefulconss && (ncuts > 0 || nfixedvars > 0) )
         break;
   }

   if ( infeasible )
   {
      SCIPdebugMsg(scip, "Infeasible node.\n");
      *result = SCIP_CUTOFF;
   }
   else if ( nfixedvars > 0 )
   {
      SCIPdebugMsg(scip, "Fixed %d variables.\n", nfixedvars);
      *result = SCIP_REDUCEDDOM;
   }
   else if ( ncuts > 0 )
   {
      SCIPdebugMsg(scip, "Separated %dinequalities.\n", ncuts);
      *result = SCIP_SEPARATED;
   }
   else
   {
      SCIPdebugMsg(scip, "No violated inequality found during separation.\n");
   }

   return SCIP_OKAY;
}


/** check whether all variables in an orbitope constraint are fixed */
static
SCIP_RETCODE checkRedundantCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be processed */
   SCIP_Bool*            redundant           /**< pointer to store whether constraint is redundant (contains no active vars) */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int i;
   int j;
   int nrows;
   int ncols;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( redundant != NULL );

   *redundant = FALSE;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->vars != NULL );
   assert( consdata->nrows > 0 );
   assert( consdata->ncols > 0 );

   vars = consdata->vars;
   nrows = consdata->nrows;
   ncols = consdata->ncols;

   /* check whether there exists an active variable in the orbitope */
   for (i = 0; i < nrows; ++i)
   {
      for (j = 0; j < ncols; ++j)
      {
         if ( SCIPvarIsActive(vars[i][j]) )
            return SCIP_OKAY;
      }
   }
   *redundant = TRUE;

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyOrbitopeFull)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrOrbitopeFull(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** frees constraint handler */
static
SCIP_DECL_CONSFREE(consFreeOrbitopeFull)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != 0 );
   assert( conshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteOrbitopeFull)
{  /*lint --e{715}*/
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransOrbitopeFull)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* create linear constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->vars, sourcedata->nrows, sourcedata->ncols,
         sourcedata->resolveprop, sourcedata->ismodelcons) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpOrbitopeFull)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Separation of full orbitope constraint handler <%s> for LP solution.\n",
      SCIPconshdlrGetName(conshdlr));

   *result = SCIP_DIDNOTRUN;

   /* if solution is integer, skip separation */
   if ( SCIPgetNLPBranchCands(scip) <= 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* separate constraints */
   SCIP_CALL( separateConstraints(scip, conshdlr, conss, nconss, nusefulconss, NULL, result, FALSE) );

   return SCIP_OKAY;
}

/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolOrbitopeFull)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Separation of full orbitope constraint handler <%s> for primal solution.\n",
      SCIPconshdlrGetName(conshdlr));

   *result = SCIP_DIDNOTFIND;

   /* separate constraints */
   SCIP_CALL( separateConstraints(scip, conshdlr, conss, nconss, nusefulconss, sol, result, FALSE) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpOrbitopeFull)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( result != NULL );

   /* we have a negative priority, so we should come after the integrality conshdlr */
   assert( SCIPgetNLPBranchCands(scip) == 0 );

   SCIPdebugMsg(scip, "Enforcement for full orbitope constraint handler <%s> for LP solution.\n",
      SCIPconshdlrGetName(conshdlr));

   *result = SCIP_FEASIBLE;

   /* separate constraints */
   SCIP_CALL( separateConstraints(scip, conshdlr, conss, nconss, nusefulconss, NULL, result, TRUE) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxOrbitopeFull)
{  /*lint --e{715}*/
   assert( result != NULL );
   assert( scip != NULL );

   SCIPdebugMsg(scip, "Enforcement for full orbitope constraint handler <%s> for relaxation solution.\n",
      SCIPconshdlrGetName(conshdlr));

   *result = SCIP_FEASIBLE;

   /* separate constraints */
   SCIP_CALL( separateConstraints(scip, conshdlr, conss, nconss, nusefulconss, sol, result, TRUE) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsOrbitopeFull)
{  /*lint --e{715}*/
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_FEASIBLE;
   if ( objinfeasible || solinfeasible )
      return SCIP_OKAY;

   /* loop through constraints */
   for (c = 0; c < nconss && *result == SCIP_FEASIBLE; ++c)
   {
      SCIP_CONS* cons;
      SCIP_CONSDATA* consdata;
      SCIP_Bool feasible;

      /* get data of constraint */
      cons = conss[c];
      assert( cons != 0 );
      consdata = SCIPconsGetData(cons);

      assert( consdata != NULL );

      /* do not enforce non-model constraints */
      if ( !consdata->ismodelcons )
         continue;

      SCIP_CALL( checkFullOrbitopeSolution(scip, cons, NULL, FALSE, &feasible) );

      *result = SCIP_INFEASIBLE;
   }

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckOrbitopeFull)
{  /*lint --e{715}*/
   int c;
   SCIP_CONSDATA* consdata;
   SCIP_Bool feasible;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_FEASIBLE;

   /* loop through constraints */
   for( c = 0; c < nconss && (*result == SCIP_FEASIBLE || completely); ++c )
   {
      assert( conss[c] != 0 );
      consdata = SCIPconsGetData(conss[c]);

      assert( consdata != NULL );

      /* do not check non-model constraints */
      if ( !consdata->ismodelcons )
         continue;

      SCIP_CALL( checkFullOrbitopeSolution(scip, conss[c], sol, printreason, &feasible) );

      if ( ! feasible )
         *result = SCIP_INFEASIBLE;
   }

   if( feasible )
   {
      SCIPdebugMsg(scip, "Solution is feasible.\n");
   }
   else
   {
      SCIPdebugMsg(scip, "Solution is infeasible.\n");
   }

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropOrbitopeFull)
{  /*lint --e{715}*/
   SCIP_Bool infeasible = FALSE;
   int nfixedvars = 0;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* propagate all useful constraints */
   for (c = 0; c < nusefulconss && !infeasible; ++c)
   {
      assert( conss[c] != 0 );

      SCIPdebugMsg(scip, "Propagation of full orbitope constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      SCIP_CALL( propagateCons(scip, conss[c], &infeasible, &nfixedvars) );
   }

   /* return the correct result */
   if ( infeasible )
   {
      *result = SCIP_CUTOFF;
      SCIPdebugMsg(scip, "Propagation via orbitopal fixing proved node to be infeasible.\n");
   }
   else if ( nfixedvars > 0 )
   {
      *result = SCIP_REDUCEDDOM;
      SCIPdebugMsg(scip, "Propagated %d variables via orbitopal fixing.\n", nfixedvars);
   }
   else if ( nusefulconss > 0 )
   {
      *result = SCIP_DIDNOTFIND;
      SCIPdebugMsg(scip, "Propagation via orbitopal fixing did not find anything.\n");
   }

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolOrbitopeFull)
{  /*lint --e{715}*/
   SCIP_Bool infeasible = FALSE;
   int noldfixedvars;
   int c;
   SCIP_Bool redundant;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;
   noldfixedvars = *nfixedvars;

   /* propagate all useful constraints
    *
    * @todo use an event handler to only propagate if a variable in the orbitope has been fixed
    */
   for (c = 0; c < nconss && !infeasible; ++c)
   {
      int nfixed = 0;

      assert( conss[c] != 0 );

      SCIPdebugMsg(scip, "Presolving of full orbitope constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      /* first propagate */
      SCIP_CALL( propagateCons(scip, conss[c], &infeasible, &nfixed) );
      *nfixedvars += nfixed;

      if ( ! infeasible )
      {
         SCIP_CALL( checkRedundantCons(scip, conss[c], &redundant) );

         if ( redundant )
         {
            SCIPdebugMsg(scip, "Full orbitope constraint <%s> is redundant: it does not contain active variables\n",
               SCIPconsGetName(conss[c]));
            SCIP_CALL( SCIPdelCons(scip, conss[c]) );
            assert( ! SCIPconsIsActive(conss[c]) );
            (*ndelconss)++;
            continue;
         }
      }
   }

   if ( infeasible )
   {
      *result = SCIP_CUTOFF;
      SCIPdebugMsg(scip, "Presolving detected infeasibility.\n");
   }
   else if ( *nfixedvars > noldfixedvars )
   {
      *result = SCIP_SUCCESS;
   }
   else if ( nconss > 0 )
   {
      *result = SCIP_DIDNOTFIND;
      SCIPdebugMsg(scip, "Presolving via orbitopal fixing did not find anything.\n");
   }

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropOrbitopeFull)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( cons != NULL );
   assert( infervar != NULL );
   assert( bdchgidx != NULL );
   assert( result != NULL );

   SCIP_CALL( resolvePropagationFullOrbitope(scip, conshdlr, cons, inferinfo, bdchgidx, result) );

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockOrbitopeFull)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int i;
   int j;
   int nrows;
   int ncols;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( locktype == SCIP_LOCKTYPE_MODEL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->nrows > 0 );
   assert( consdata->ncols > 0 );
   assert( consdata->vars != NULL );

   SCIPdebugMsg(scip, "Locking method for full orbitope constraint handler\n");

   nrows = consdata->nrows;
   ncols = consdata->ncols;
   vars = consdata->vars;

   /* add up locks and down locks on each variable */
   for (i = 0; i < nrows; ++i)
   {
      for (j = 0; j < ncols; ++j)
         SCIP_CALL( SCIPaddVarLocksType(scip, vars[i][j], locktype, nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   return SCIP_OKAY;
}


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintOrbitopeFull)
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int i;
   int j;
   int nrows;
   int ncols;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->nrows > 0 );
   assert( consdata->ncols > 0 );
   assert( consdata->vars != NULL );

   nrows = consdata->nrows;
   ncols = consdata->ncols;
   vars = consdata->vars;

   SCIPdebugMsg(scip, "Printing method for full orbitope constraint handler\n");

   SCIPinfoMessage(scip, file, "fullOrbitope(");

   for (i = 0; i < nrows; ++i)
   {
      for (j = 0; j < ncols; ++j)
      {
         if ( j > 0 )
            SCIPinfoMessage(scip, file, ",");
         SCIP_CALL( SCIPwriteVarName(scip, file, vars[i][j], TRUE) );
      }
      if ( i < nrows-1 )
         SCIPinfoMessage(scip, file, ".");
   }
   SCIPinfoMessage(scip, file, ")");

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyOrbitopeFull)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   SCIP_VAR*** sourcevars;
   SCIP_VAR*** vars;
   int nrows;
   int ncols;
   int i;
   int k;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( sourcescip != NULL );
   assert( sourceconshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(sourceconshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );
   assert( varmap != NULL );
   assert( valid != NULL );

   *valid = TRUE;

   SCIPdebugMsg(scip, "Copying method for full orbitope constraint handler.\n");

   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );
   assert( sourcedata->nrows > 0 );
   assert( sourcedata->ncols > 0 );
   assert( sourcedata->vars != NULL );

   conshdlrdata = SCIPconshdlrGetData(sourceconshdlr);
   assert( conshdlrdata != NULL );

   /* do not copy non-model constraints */
   if ( !sourcedata->ismodelcons && !conshdlrdata->forceconscopy )
   {
      *valid = FALSE;

      return SCIP_OKAY;
   }

   nrows = sourcedata->nrows;
   ncols = sourcedata->ncols;
   sourcevars = sourcedata->vars;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nrows) );
   for (i = 0; i < nrows && *valid; ++i)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(vars[i]), ncols) );  /*lint !e866*/

      for (j = 0; j < ncols && *valid; ++j)
      {
         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[i][j], &(vars[i][j]), varmap, consmap, global, valid) );
         assert( !(*valid) || vars[i][j] != NULL );
      }
   }

   /* only create the target constraint, if all variables could be copied */
   if ( *valid )
   {
      /* create copied constraint */
      if ( name == NULL )
         name = SCIPconsGetName(sourcecons);

      SCIP_CALL( SCIPcreateConsOrbitopeFull(scip, cons, name, vars, nrows, ncols,
            sourcedata->resolveprop, sourcedata->ismodelcons,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   /* free space; only up to row i if copying failed */
   assert( 0 <= i && i <= nrows );
   for (k = i - 1; k >= 0; --k)
      SCIPfreeBufferArray(scip, &vars[k]);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseOrbitopeFull)
{  /*lint --e{715}*/
   const char* s;
   char* endptr;
   SCIP_VAR*** vars;
   SCIP_VAR* var;
   int nrows;
   int maxnrows;
   int ncols;
   int maxncols;
   int k;
   int j;

   assert( success != NULL );

   *success = TRUE;
   s = str;

   /* skip white space */
   SCIP_CALL( SCIPskipSpace((char**)&s) );

   if( strncmp(s, "fullOrbitope(", 13) != 0 )
   {
      SCIPerrorMessage("Syntax error - expected \"fullOrbitope(\": %s\n", s);
      *success = FALSE;
      return SCIP_OKAY;
   }
   s += 13;

   /* loop through string */
   nrows = 0;
   ncols = 0;
   maxnrows = 10;
   maxncols = 10;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, maxnrows) );

   j = 0;
   do
   {
      /* parse variable name */
      SCIP_CALL( SCIPparseVarName(scip, s, &var, &endptr) );

      if( var == NULL )
      {
         endptr = strchr(endptr, ')');

         if( endptr == NULL || j > 0 )
         {
            SCIPerrorMessage("not enough variables.\n");
            *success = FALSE;
         }

         break;
      }

      s = endptr;
      assert( s != NULL );

      /* skip white space */
      SCIP_CALL( SCIPskipSpace((char**)&s) );

      /* begin new row if required */
      if( j == 0 )
      {
         ++nrows;

         if( nrows > maxnrows )
         {
            maxnrows = SCIPcalcMemGrowSize(scip, nrows);
            SCIP_CALL( SCIPreallocBufferArray(scip, &vars, maxnrows) );
            assert( nrows <= maxnrows );
         }

         SCIP_CALL( SCIPallocBufferArray(scip, &(vars[nrows-1]), nrows == 1 ? maxncols : ncols) ); /*lint !e866*/
      }

      /* determine number of columns */
      if( nrows == 1 )
      {
         ncols = j+1;

         if( *s == '.' || *s == ')' )
            SCIP_CALL( SCIPreallocBufferArray(scip, &(vars[nrows-1]), ncols) ); /*lint !e866*/
         else if( ncols > maxncols )
         {
            maxncols = SCIPcalcMemGrowSize(scip, ncols);
            SCIP_CALL( SCIPreallocBufferArray(scip, &(vars[nrows-1]), maxncols) ); /*lint !e866*/
            assert( ncols <= maxncols );
         }
      }
      else if( ( j < ncols-1 ) == ( *s == '.' || *s == ')' ) )
      {
         SCIPerrorMessage("variables per row do not match.\n");
         *success = FALSE;
         break;
      }

      vars[nrows-1][j] = var;

      if( *s == '.' )
         j = 0;
      else
         ++j;

      /* skip ',' or '.' */
      if( *s == ',' || *s == '.' )
         ++s;
   }
   while( *s != ')' );

   if( *success )
      SCIP_CALL( SCIPcreateConsOrbitopeFull(scip, cons, name, vars, nrows, ncols, TRUE, TRUE,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

   for( k = nrows - 1; k >= 0; --k )
      SCIPfreeBufferArray(scip, &vars[k]);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsOrbitopeFull)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );
   assert( success != NULL );
   assert( vars != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   if ( varssize < consdata->ncols * consdata->nrows )
      (*success) = FALSE;
   else
   {
      int cnt = 0;
      int i;
      int j;

      for (i = 0; i < consdata->nrows; ++i)
      {
         for (j = 0; j < consdata->ncols; ++j)
            vars[cnt++] = consdata->vars[i][j];
      }
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsOrbitopeFull)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   (*nvars) = consdata->ncols * consdata->nrows;
   (*success) = TRUE;

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for orbitope constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrOrbitopeFull(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create orbitope constraint handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpOrbitopeFull, consEnfopsOrbitopeFull, consCheckOrbitopeFull, consLockOrbitopeFull,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyOrbitopeFull, consCopyOrbitopeFull) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeOrbitopeFull) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteOrbitopeFull) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsOrbitopeFull) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsOrbitopeFull) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseOrbitopeFull) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolOrbitopeFull,
         CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintOrbitopeFull) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropOrbitopeFull, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropOrbitopeFull) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpOrbitopeFull, consSepasolOrbitopeFull, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransOrbitopeFull) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxOrbitopeFull) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/forceconscopy",
         "Whether orbitope constraints should be forced to be copied to sub SCIPs.",
         &conshdlrdata->forceconscopy, TRUE, DEFAULT_FORCECONSCOPY, NULL, NULL) );

   return SCIP_OKAY;
}


/** creates and captures a full orbitope constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsOrbitopeFull(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR***           vars,               /**< matrix of variables on which the symmetry acts */
   int                   nrows,            /**< number of set partitioning/packing constraints  <=> p */
   int                   ncols,            /**< number of symmetric variable blocks             <=> q */
   SCIP_Bool             resolveprop,        /**< should propagation be resolved? */
   SCIP_Bool             ismodelcons,        /**< whether the orbitope is a model constraint */
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

   /* find the orbitope constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("full orbitope constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   assert( nrows > 0 );
   assert( ncols > 0 );

   /* run some checks */
#ifndef NDEBUG
   {
      SCIP_Real obj;
      int i;
      int j;
      for (i = 0; i < nrows; ++i)
      {
         /* initialize obj to infinity */
         obj = SCIPinfinity(scip);
         for (j = 0; j < ncols; ++j)
         {
            SCIP_Bool fixedZero;
            SCIP_VAR* var;

            var = vars[i][j];
            assert(var != NULL);

            if ( SCIPvarIsNegated(var) )
               var = SCIPvarGetNegatedVar(var);

            /* all variables need to be binary */
            assert( SCIPvarIsBinary(var) );

            /* fixed variables have obj = 0; for variables fixed to 0, we assume that there is no
               problem (but we cannot always check it, e.g., when in the original problem
               variables were fixed and this problem was copied.) */
            fixedZero = ( SCIPisZero(scip, SCIPvarGetLbGlobal(var)) && SCIPisZero(scip, SCIPvarGetUbGlobal(var)) );

            /* @todo adapt correctness of the following check for sub-scips */
            if ( SCIPgetSubscipDepth(scip) == 0 )
            {
               /* check whether all variables in a row have the same objective */
               if ( ! fixedZero && SCIPisInfinity(scip, obj) )
                  obj = SCIPvarGetObj(var);
               else
               {
                  assert( fixedZero || ! SCIPvarIsActive(var) || SCIPisEQ(scip, obj, SCIPvarGetObj(var)) );
               }
            }
         }
      }
   }
#endif

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, vars, nrows, ncols, resolveprop, ismodelcons) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a full orbitope constraint in its most basic variant, i. e., with all constraint flags set to
 *  their default values, which can be set afterwards using SCIPsetConsFLAGNAME()
 *
 *  @see SCIPcreateConsOrbitopeFull() for the default constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicOrbitopeFull(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR***           vars,               /**< matrix of variables on which the symmetry acts */
   int                   nrows,              /**< number of set partitioning/packing constraints  <=> p */
   int                   ncols,              /**< number of symmetric variable blocks             <=> q */
   SCIP_Bool             resolveprop,        /**< should propagation be resolved? */
   SCIP_Bool             ismodelcons         /**< whether the orbitope is a model constraint */
   )
{
   SCIP_CALL( SCIPcreateConsOrbitopeFull(scip, cons, name, vars, nrows, ncols,
         resolveprop, ismodelcons, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

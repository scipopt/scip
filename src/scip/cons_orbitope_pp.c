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

/**@file   cons_orbitope_pp.c
 * @ingroup DEFPLUGINS_CONS
 * @brief  constraint handler for partitioning/packing orbitope constraints w.r.t. the full symmetric group
 * @author Timo Berthold
 * @author Marc Pfetsch
 * @author Christopher Hojny
 *
 * The type of constraints of this constraint handler is described in cons_orbitope_pp.h.
 * When creating the constraint, users can decide whether it is a constraint defining the model
 * or "just" use to handle symmetries. In the latter case, symmetry reductions are only performed
 * by the constraint handler if strong dual reductions are permitted.
 *
 * The details of the method implemented here are described in the following papers.
 *
 * Packing and Partitioning Orbitopes@n
 * Volker Kaibel and Marc E. Pfetsch,@n
 * Math. Program. 114, No. 1, 1-36 (2008)
 *
 * Among other things, this paper describes so-called shifted column inequalities of the following
 * form \f$x(S) \leq x(B)\f$, where \f$S\f$ is a so-called shifted column and \f$B\f$ is a so-called
 * bar. These inequalities can be used to handle symmetry and they are separated in this constraint
 * handler. We use the linear time separation algorithm of the paper.@par
 *
 * Orbitopal Fixing@n
 * Volker Kaibel, Matthias Peinhardt, and Marc E. Pfetsch,@n
 * Discrete Optimization 8, No. 4, 595-610 (2011)
 * (A preliminary version appears in Proc. IPCO 2007.)
 *
 * In this paper a linear time propagation algorithm is described, a variant of which is implemented
 * here. The implemented variant does not run in linear time, but is very fast in practice.
 *
 * <table>
 *   <caption>translation table</caption>
 *   <tr><td>here</td><td>paper</td></tr>
 *   <tr><td></td><td></td></tr>
 *   <tr><td>nrows        </td><td>p       </td></tr>
 *   <tr><td>ncols        </td><td>q       </td></tr>
 *   <tr><td>vars         </td><td>x       </td></tr>
 *   <tr><td>vals         </td><td>A^\\star</td></tr>
 *   <tr><td>weights      </td><td>\\omega </td></tr>
 *   <tr><td>cases        </td><td>\\tau   </td></tr>
 *   <tr><td>fixtriangle  </td><td>--      </td></tr>
 *   <tr><td>resolveprop  </td><td>--      </td></tr>
 *   <tr><td>firstnonzeros</td><td>\\mu    </td></tr>
 *   <tr><td>lastones     </td><td>\\alpha </td></tr>
 *   <tr><td>frontiersteps</td><td>\\Gamma </td></tr>
 * </table>
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/cons_orbitope_pp.h"
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
#define CONSHDLR_NAME          "orbitope_pp"
#define CONSHDLR_DESC          "symmetry breaking constraint handler relying on partitioning/packing orbitopes"
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
   SCIP_VAR***           vars;               /**< matrix of variables on which the symmetry acts                     */
   SCIP_VAR**            tmpvars;            /**< temporary storage for variables                                    */
   SCIP_Real**           vals;               /**< LP-solution for those variables                                    */
   SCIP_Real*            tmpvals;            /**< temporary storage for values                                       */
   SCIP_Real**           weights;            /**< SC weight table                                                    */
   int**                 cases;              /**< indicator of the SC cases                                          */
   int                   nrows;              /**< number of rows in orbitope matrix               <=> p              */
   int                   ncols;              /**< number of columns in orbitope matrix            <=> q              */
   SCIP_ORBITOPETYPE     orbitopetype;       /**< type of orbitope constraint                                        */
   SCIP_Bool             resolveprop;        /**< should propagation be resolved?                                    */
   SCIP_Bool             istrianglefixed;    /**< has the upper right triangle already globally been fixed to zero?  */
   SCIP_Bool             ismodelcons;        /**< whether the orbitope is a model constraint                         */
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
   int nrows;
   int ncols;

   assert( consdata != NULL );
   assert( *consdata != NULL );

   nrows = (*consdata)->nrows;
   ncols = (*consdata)->ncols;
   for (i = 0; i < nrows; ++i)
   {
      /* release variables in vars array */
      for (j = 0; j < ncols; ++j)
      {
         assert( (*consdata)->vars[i] != NULL );
         SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->vars[i][j]) );
      }

      SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->cases[i]), ncols);    /*lint !e866*/
      SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vars[i]), ncols);     /*lint !e866*/
      SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->weights[i]), ncols);  /*lint !e866*/
      SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vals[i]), ncols);     /*lint !e866*/
   }

   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->cases), nrows);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vars), nrows);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->weights), nrows);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vals), nrows);

   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->tmpvals), nrows + ncols);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->tmpvars), nrows + ncols);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** creates orbitope constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure                                     */
   SCIP_CONSDATA**       consdata,           /**< pointer to store constraint data                        */
   SCIP_VAR***           vars,               /**< variables array, must have size nspcons x nblocks       */
   int                   nrows,              /**< number of rows in orbitope matrix               <=> p   */
   int                   ncols,              /**< number of columns in orbitope matrix            <=> q   */
   SCIP_ORBITOPETYPE     orbitopetype,       /**< type of orbitope constraint                             */
   SCIP_Bool             resolveprop,        /**< should propagation be resolved?                         */
   SCIP_Bool             ismodelcons         /**< whether the orbitope is a model constraint              */
   )
{
   int i;
   int j;

   assert(consdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->vals, nrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->weights, nrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->vars, nrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->cases, nrows) );

   for (i = 0; i < nrows; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->vals[i], ncols) );                 /*lint !e866*/
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->weights[i], ncols) );              /*lint !e866*/
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars[i], vars[i], ncols) );    /*lint !e866*/
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->cases[i], ncols) );                /*lint !e866*/
   }
   (*consdata)->tmpvals = NULL;
   (*consdata)->tmpvars = NULL;
   (*consdata)->nrows = nrows;
   (*consdata)->ncols = ncols;
   (*consdata)->orbitopetype = orbitopetype;
   (*consdata)->resolveprop = resolveprop;
   (*consdata)->istrianglefixed = FALSE;
   (*consdata)->ismodelcons = ismodelcons;

   /* get transformed variables, if we are in the transformed problem */
   if ( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->tmpvals, nrows + ncols) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->tmpvars, nrows + ncols) );

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


#ifdef PRINT_MATRIX
/** debug method, prints variable matrix */
static
void printMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< the constraint data */
   )
{
   int i;
   int j;

   assert( consdata != NULL );
   assert( consdata->nrows > 0 );
   assert( consdata->ncols > 0 );
   assert( consdata->vars != NULL );

   for (j = 0; j < consdata->ncols; ++j)
      SCIPinfoMessage(scip, NULL, "-");

   SCIPinfoMessage(scip, NULL, "\n");
   for (i = 0; i < consdata->nrows; ++i)
   {
      for (j = 0; j < consdata->ncols; ++j)
      {
         if ( SCIPvarGetUbLocal(consdata->vars[i][j]) - SCIPvarGetLbLocal(consdata->vars[i][j]) < 0.5 )
            SCIPinfoMessage(scip, NULL, "%1.0f", REALABS(SCIPvarGetUbLocal(consdata->vars[i][j])));
         else
            SCIPinfoMessage(scip, NULL, " ");
      }
      SCIPinfoMessage(scip, NULL, "|\n");
   }
   for (j = 0; j < consdata->ncols; ++j)
      SCIPinfoMessage(scip, NULL, "-");
   SCIPinfoMessage(scip, NULL, "\n");
}
#endif


#ifdef SHOW_SCI
/** Print SCI in nice form for debugging */
static
SCIP_RETCODE printSCI(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   nrows,              /**< number of rows */
   int                   ncols,              /**< number of columns */
   int**                 cases,              /**< SCI dynamic programming table */
   int                   i,                  /**< row position of bar */
   int                   j                   /**< column position of bar */
   )
{
   int k;
   int l;
   int** M;
   int p1;
   int p2;

   SCIP_CALL( SCIPallocBufferArray(scip, &M, nrows) );
   for (k = 0; k < nrows; ++k)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &M[k], ncols) ); /*lint !e866*/
      for (l = 0; l < ncols; ++l)
         M[k][l] = 0;
   }

   /* first add bar */
   for (l = j; l < ncols; ++l)
   {
      assert( M[i][l] == 0 );
      M[i][l] = 1;
   }

   /* then add shifted column */
   p1 = i-1;
   p2 = j-1;
   do
   {
      assert( cases[p1][p2] != -1 );
      assert( p1 >= 0 && p1 < i );
      assert( p2 >= 0 && p2 < j );

      /* if case 1 */
      if ( cases[p1][p2] == 1 )
         --p2;   /* decrease column */
      else
      {
         /* case 2 or 3: */
         assert( cases[p1][p2] == 2 || cases[p1][p2] == 3 );
         assert( M[p1][p2] == 0 );
         M[p1][p2] = -1;
         if ( cases[p1][p2] == 3 )
            break;
      }
      --p1;  /* decrease row */
   }
   while ( p1 >= 0 );   /* should always be true, i.e. the break should end the loop */
   assert( cases[p1][p2] == 3 );

   /* now output matrix M */
   for (l = 0; l < ncols; ++l)
      SCIPinfoMessage(scip, NULL, "-");
   SCIPinfoMessage(scip, NULL, "\n");

   for (k = 0; k < nrows; ++k)
   {
      for (l = 0; l < ncols; ++l)
      {
         if ( l > k )
            SCIPinfoMessage(scip, NULL, "*");
         else
         {
            switch (M[k][l])
            {
            case 1:
               SCIPinfoMessage(scip, NULL, "+");
               break;
            case -1:
               SCIPinfoMessage(scip, NULL, "-");
               break;
            case 0:
               SCIPinfoMessage(scip, NULL, "#");
               break;
            default:
               SCIPerrorMessage("unexpected matrix entry <%d>: should be -1, 0 or +1\n", M[k][l]);
               SCIPABORT();
            }
         }
      }
      SCIPinfoMessage(scip, NULL, "\n");
   }

   for (l = 0; l < ncols; ++l)
      SCIPinfoMessage(scip, NULL, "-");
   SCIPinfoMessage(scip, NULL, "\n");

   for (k = 0; k < nrows; ++k)
      SCIPfreeBufferArray(scip, &M[k]);
   SCIPfreeBufferArray(scip, &M);

   return SCIP_OKAY;
}
#endif


/** copies the variables values from the solution to the constraint data structure */
static
void copyValues(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< the constraint data     */
   SCIP_SOL*             sol                 /**< a primal solution or NULL for the current LP optimum */
   )
{
   int i;
   int j;

   assert( scip != NULL );
   assert( consdata != NULL );
   assert( consdata->nrows > 0 );
   assert( consdata->ncols > 0 );
   assert( consdata->vars != NULL );
   assert( consdata->vals != NULL );

   for (i = 0; i < consdata->nrows; ++i)
   {
      for (j = 0; j < consdata->ncols; ++j)
         consdata->vals[i][j] = SCIPgetSolVal(scip, sol, consdata->vars[i][j]);
   }
}


/** compute the dynamic programming table for SC
 *
 *  Build up dynamic programming table in order to find SCs with minimum weight.
 *
 *  The values of the minimal SCIs are stored in @a weights.
 *  The array @a cases[i][j] stores which of the cases were applied to get @a weights[i][j].
 *  Here, 3 means that we have reached the upper limit.
 *
 *  We assume that the upper right triangle is fixed to 0. Hence we can perform the computation a
 *  bit more efficient.
 */
static
void computeSCTable(
   SCIP*                 scip,               /**< SCIP pointer                                              */
   int                   nrows,              /**< number of rows in orbitope matrix                 <=> p   */
   int                   ncols,              /**< number of columns in orbitope matrix              <=> q   */
   SCIP_Real**           weights,            /**< SC weight table                                           */
   int**                 cases,              /**< indicator of the SC cases                                 */
   SCIP_Real**           vals                /**< current solution                                          */
   )
{
   SCIP_Real minvalue;
   int diagsize;
   int i;
   int j;

   assert( weights != NULL );
   assert( cases != NULL );
   assert( vals != NULL );

#ifndef NDEBUG
   /* for debugging */
   for (i = 0; i < nrows; ++i)
   {
      for (j = 0; j < ncols; ++j)
      {
         if ( i >= j )
         {
            weights[i][j] = -1.0;
            cases[i][j] = -1;
         }
      }
   }
#endif

   /* initialize diagonal */
   minvalue = vals[0][0];
   weights[0][0] = minvalue;
   cases[0][0] = 3;

   /* get last row of triangle */
   diagsize = ncols;
   if ( nrows < ncols )
      diagsize = nrows;

   for (j = 1; j < diagsize; ++j)
   {
      /* use LT to move entry as far to the left as possible */
      if ( SCIPisLT(scip, vals[j][j], minvalue) )
      {
         minvalue = vals[j][j];
         cases[j][j] = 3;
      }
      else
         cases[j][j] = 1;
      weights[j][j] = minvalue;
   }

   /* initialize first column */
   for (i = 1; i < nrows; ++i)
   {
      weights[i][0] = weights[i-1][0] + vals[i][0];
      cases[i][0] = 2;  /* second case */
   }

   /* build the table */
   for (i = 2; i < nrows; ++i)
   {
      for (j = 1; j < ncols && j < i; ++j)
      {
         SCIP_Real weightleft;
         SCIP_Real weightright;

         assert( cases[i-1][j] != -1 );
         assert( cases[i-1][j-1] != -1 );

         weightleft = weights[i-1][j-1];
         weightright = vals[i][j] + weights[i-1][j];

         /* For first column: cannot take left possibility */
         if ( SCIPisLT(scip, weightleft, weightright) )
         {
            weights[i][j] = weightleft;
            cases[i][j] = 1;
         }
         else
         {
            weights[i][j] = weightright;
            cases[i][j] = 2;
         }
      }
   }
}


/** fix upper right triangle if necessary */
static
SCIP_RETCODE fixTriangle(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be processed */
   SCIP_Bool*            infeasible,         /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars          /**< pointer to add up the number of found domain reductions */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   SCIP_Bool fixedglobal;
   SCIP_Bool fixed;
   int diagsize;
   int nrows;
   int ncols;
   int i;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );
   assert( nfixedvars != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->nrows > 0 );
   assert( consdata->ncols > 0 );
   assert( consdata->vars != NULL );

   *infeasible = FALSE;
   *nfixedvars = 0;

   if ( consdata->istrianglefixed )
      return SCIP_OKAY;

   nrows = consdata->nrows;
   ncols = consdata->ncols;
   vars = consdata->vars;
   fixedglobal = TRUE;

   /* get last row of triangle */
   diagsize = ncols;
   if ( nrows < ncols )
      diagsize = nrows;

   /* fix variables to 0 */
   for (i = 0; i < diagsize; ++i)
   {
      for (j = i+1; j < ncols; ++j)
      {
         /* fix variable, if not in the root the fixation is local */
         SCIP_CALL( SCIPfixVar(scip, vars[i][j], 0.0, infeasible, &fixed) );

         if ( *infeasible )
         {
            SCIPdebugMsg(scip, "The problem is infeasible: some variable in the upper right triangle is fixed to 1.\n");
            return SCIP_OKAY;
         }

         if ( fixed )
            ++(*nfixedvars);

         if ( SCIPvarGetUbGlobal(vars[i][j]) > 0.5 )
            fixedglobal = FALSE;
      }
   }
   if ( *nfixedvars > 0 )
   {
      SCIPdebugMsg(scip, "<%s>: %s fixed upper right triangle to 0 (fixed vars: %d).\n",
         SCIPconsGetName(cons), fixedglobal ? "globally" : "locally", *nfixedvars);
   }
   else
   {
      SCIPdebugMsg(scip, "<%s>: Upper right triangle already fixed to 0.\n", SCIPconsGetName(cons));
   }

   if ( fixedglobal )
      consdata->istrianglefixed = TRUE;

   return SCIP_OKAY;
}


/** separates shifted column inequalities according to the solution stored in consdata->vals */
static
SCIP_RETCODE separateSCIs(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< the constraint data */
   SCIP_Bool*            infeasible,         /**< whether we detected infeasibility */
   int*                  nfixedvars,         /**< pointer to store the number of variables fixed */
   int*                  ncuts               /**< pointer to store number of separated SCIs */
   )
{
   SCIP_Real** vals;
   SCIP_Real** weights;
   SCIP_Real* tmpvals;
   SCIP_VAR*** vars;
   SCIP_VAR** tmpvars;
   int** cases;
   int nrows;
   int ncols;
   int i;
   int j;
   int l;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL);
   assert( nfixedvars != NULL );
   assert( ncuts != NULL );

   assert( consdata != NULL );
   assert( consdata->nrows > 0 );
   assert( consdata->ncols > 0 );
   assert( consdata->vars != NULL );
   assert( consdata->vals != NULL );
   assert( consdata->tmpvars != NULL );
   assert( consdata->tmpvals != NULL );
   assert( consdata->weights != NULL );
   assert( consdata->cases != NULL );

   *infeasible = FALSE;
   *nfixedvars = 0;
   *ncuts = 0;

   nrows = consdata->nrows;
   ncols = consdata->ncols;
   vars = consdata->vars;
   vals = consdata->vals;
   tmpvars = consdata->tmpvars;
   tmpvals = consdata->tmpvals;
   weights = consdata->weights;
   cases = consdata->cases;

   /* check for upper right triangle */
   if ( ! consdata->istrianglefixed )
   {
      SCIP_CALL( fixTriangle(scip, cons, infeasible, nfixedvars) );
      if ( *infeasible )
         return SCIP_OKAY;
      if ( *nfixedvars > 0 )
         return SCIP_OKAY;
   }

   /* compute table if necessary (i.e., not computed before) */
   computeSCTable(scip, nrows, ncols, weights, cases, vals);

   /* loop through rows */
   for (i = 1; i < nrows && ! (*infeasible); ++i)
   {
      SCIP_Real bar;       /* value of bar: */
      int lastcolumn;      /* last column considered as part of the bar */

      bar = 0.0;
      lastcolumn = ncols - 1;
      if ( lastcolumn > i )
         lastcolumn = i;

      /* traverse row from right to left: */
      /* j >= 1, since for j = 0, i.e., the bar is a complete row, there does not exist an SCI */
      for (j = lastcolumn; j > 0; --j)
      {
         bar += vals[i][j];

         /* check whether weights[i-1][j-1] < bar  (<=> bar - weights[i-1][j-1] > 0), i.e. cut is violated) */
         if ( SCIPisEfficacious(scip, bar - weights[i-1][j-1]) )
         {
#ifndef NDEBUG
            SCIP_Real weight = 0.0;
#endif
            SCIP_ROW* row;
#ifdef SCIP_DEBUG
            char name[SCIP_MAXSTRLEN];
#endif
            int nvars;
            int p1;
            int p2;

            nvars = 0;
            p1 = i-1;
            p2 = j-1;

            /* first add bar */
            for (l = j; l <= lastcolumn; ++l)
            {
               tmpvars[nvars] = vars[i][l];
               tmpvals[nvars] = 1.0;
               nvars++;
            }

            /* then add shifted column */
            do
            {
               assert( cases[p1][p2] != -1 );
               assert( p1 >= 0 && p1 < i );
               assert( p2 >= 0 && p2 < j );

               /* if case 1 */
               if (cases[p1][p2] == 1)
                  p2--;   /* decrease column */
               else
               {
                  /* case 2 or 3: */
                  assert( cases[p1][p2] == 2 || cases[p1][p2] == 3 );
                  tmpvars[nvars] = vars[p1][p2];
                  tmpvals[nvars] = -1.0;
                  nvars++;
#ifndef NDEBUG
                  weight += vals[p1][p2];
#endif
                  if ( cases[p1][p2] == 3 )
                     break;
               }
               p1--;  /* decrease row */
            }
            while ( p1 >= 0 );   /* should always be true, i.e. the break should end the loop */
            assert( cases[p1][p2] == 3 );

            /* generate cut */
#ifdef SCIP_DEBUG
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sci_%d_%d", i, j);
            SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, name, -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
#else
            SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, "", -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
#endif
            SCIP_CALL( SCIPaddVarsToRow(scip, row, nvars, tmpvars, tmpvals) );
            /*SCIP_CALL( SCIPprintRow(scip, row, NULL) ); */
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
            SCIP_CALL( SCIPreleaseRow(scip, &row) );
            ++(*ncuts);

#ifdef SHOW_SCI
            SCIP_CALL( printSCI(scip, nrows, ncols, cases, i, j) );
#endif

            assert( SCIPisSumEQ(scip, weights[i-1][j-1], weight) );
         }
      }
   }
   return SCIP_OKAY;
}


/** propagation method for a single orbitope constraint */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be processed */
   SCIP_Bool*            infeasible,         /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars          /**< pointer to add up the number of found domain reductions */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_ORBITOPETYPE orbitopetype;
   SCIP_VAR*** vars;
   int* firstnonzeros;
   int* lastones;
   int* frontiersteps;
   int lastoneprevrow;
   int nrows;
   int ncols;
   int nsteps;
   int i;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );
   assert( nfixedvars != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   *nfixedvars = 0;

   /* if the constraint is not a model constraint, check whether symmetry reductions are permitted */
   if( !consdata->ismodelcons && !SCIPallowStrongDualReds(scip) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->nrows > 0 );
   assert( consdata->ncols > 0 );
   assert( consdata->vars != NULL );

   nrows = consdata->nrows;
   ncols = consdata->ncols;
   vars = consdata->vars;
   orbitopetype = consdata->orbitopetype;

   assert( orbitopetype == SCIP_ORBITOPETYPE_PACKING || orbitopetype == SCIP_ORBITOPETYPE_PARTITIONING );

   /* fix upper right triangle if still necessary */
   if ( ! consdata->istrianglefixed )
   {
      int nfixed = 0;
      SCIP_CALL( fixTriangle(scip, cons, infeasible, &nfixed) );
      *nfixedvars += nfixed;
   }

   /* prepare further propagation */
   SCIP_CALL( SCIPallocBufferArray(scip, &firstnonzeros, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lastones, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &frontiersteps, ncols) );

#ifdef PRINT_MATRIX
   SCIPdebugMsg(scip, "Matrix:\n");
   printMatrix(scip, consdata);
#endif

   /* propagate */
   lastoneprevrow = 0;
   lastones[0] = 0;

   if ( orbitopetype == SCIP_ORBITOPETYPE_PACKING )
   {
      /* packing case: if entry (0,0) is fixed to 0 */
      if ( SCIPvarGetUbLocal(vars[0][0]) < 0.5 )
      {
         lastoneprevrow = -1;
         lastones[0] = -1;
      }
   }
   nsteps = 0;

   for (i = 1; i < nrows; ++i)
   {
      int lastcolumn;
      int firstnonzeroinrow;
      int lastoneinrow;
      SCIP_Bool infrontier;

      /* last column considered as part of the bar: */
      lastcolumn = ncols - 1;
      if ( lastcolumn > i )
         lastcolumn = i;

      /* find first position not fixed to 0 (partitioning) or fixed to 1 (packing) */
      firstnonzeroinrow = -1;
      for (j = 0; j <= lastcolumn; ++j)
      {
         if ( orbitopetype == SCIP_ORBITOPETYPE_PARTITIONING )
         {
            /* partitioning case: if variable is not fixed to 0 */
            if ( SCIPvarGetUbLocal(vars[i][j]) > 0.5 )
            {
               firstnonzeroinrow = j;
               break;
            }
         }
         else
         {
            /* packing case: if variable is fixed to 1 */
            if ( SCIPvarGetLbLocal(vars[i][j]) > 0.5 )
            {
               firstnonzeroinrow = j;
               break;
            }
         }
      }
      /* if all variables are fixed to 0 in the partitioning case - should not happen */
      if ( firstnonzeroinrow == -1 && orbitopetype == SCIP_ORBITOPETYPE_PARTITIONING )
      {
         SCIPdebugMsg(scip, " -> Infeasible node: all variables in row %d are fixed to 0.\n", i);
         *infeasible = TRUE;
         /* conflict should be analyzed by setppc constraint handler */
         goto TERMINATE;
      }
      firstnonzeros[i] = firstnonzeroinrow;
      assert( orbitopetype == SCIP_ORBITOPETYPE_PACKING || firstnonzeroinrow >= 0 );
      assert( -1 <= firstnonzeroinrow && firstnonzeroinrow <= lastcolumn );

      /* compute rightmost possible position for a 1 */
      assert( orbitopetype == SCIP_ORBITOPETYPE_PACKING || 0 <= lastoneprevrow );
      assert( lastoneprevrow <= lastcolumn );

      /* if we are at right border or if entry in column lastoneprevrow+1 is fixed to 0 */
      infrontier = FALSE;
      assert( lastoneprevrow + 1 >= 0 );
      if ( lastoneprevrow == ncols-1 || SCIPvarGetUbLocal(vars[i][lastoneprevrow+1]) < 0.5 ) /*lint !e679*/
         lastoneinrow = lastoneprevrow;
      else
      {
         lastoneinrow = lastoneprevrow + 1;
         frontiersteps[nsteps++] = i;
         infrontier = TRUE;
      }

      /* store lastoneinrow */
      assert( orbitopetype == SCIP_ORBITOPETYPE_PACKING || 0 <= lastoneinrow );
      assert( lastoneinrow <= lastcolumn );
      lastones[i] = lastoneinrow;

      /* check whether we are infeasible */
      if ( firstnonzeroinrow > lastoneinrow )
      {
         int k;

#ifdef SCIP_DEBUG
         if ( orbitopetype == SCIP_ORBITOPETYPE_PARTITIONING )
         {
            SCIPdebugMsg(scip, " -> Infeasible node: row %d, leftmost nonzero at %d, rightmost 1 at %d\n",
               i, firstnonzeroinrow, lastoneinrow);
         }
         else
         {
            SCIPdebugMsg(scip, " -> Infeasible node: row %d, 1 at %d, rightmost position for 1 at %d\n",
               i, firstnonzeroinrow, lastoneinrow);
         }
#endif

         /* check if conflict analysis is applicable */
         if ( SCIPisConflictAnalysisApplicable(scip) )
         {
            /* conflict analysis only applicable in SOLVING stage */
            assert( SCIPgetStage(scip) == SCIP_STAGE_SOLVING || SCIPinProbing(scip) );

            /* perform conflict analysis */
            SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

            if ( orbitopetype == SCIP_ORBITOPETYPE_PARTITIONING )
            {
               /* add bounds (variables fixed to 0) that result in the first nonzero entry */
               for (j = 0; j <= lastcolumn; ++j)
               {
                  /* add varaibles in row up to the first variable fixed to 0 */
                  if ( SCIPvarGetUbLocal(vars[i][j]) > 0.5 )
                     break;

                  assert( SCIPvarGetUbLocal(vars[i][j]) < 0.5 );
                  SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i][j]) );
               }
            }
            else
            {
               /* add bounds that result in the last one - check top left entry for packing case */
               if ( lastones[0] == -1 )
               {
                  assert( SCIPvarGetUbLocal(vars[0][0]) < 0.5 );
                  SCIP_CALL( SCIPaddConflictBinvar(scip, vars[0][0]) );
               }

               /* mark variable fixed to 1 */
               assert( SCIPvarGetLbLocal(vars[i][firstnonzeroinrow]) > 0.5 );
               SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i][firstnonzeroinrow]) );
            }

            /* add bounds that result in the last one - pass through rows */
            for (k = 1; k < i; ++k)
            {
               int l;
               l = lastones[k] + 1;

               /* if the frontier has not moved and we are not beyond the matrix boundaries */
               if ( l <= ncols-1 && l <= k && lastones[k-1] == lastones[k] )
               {
                  assert( SCIPvarGetUbLocal(vars[k][l]) < 0.5 );
                  SCIP_CALL( SCIPaddConflictBinvar(scip, vars[k][l]) );
               }
            }
            SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
         }

         *infeasible = TRUE;
         goto TERMINATE;
      }

      /* fix entries beyond the last possible position for a 1 in the row to 0 (see Lemma 1 in the paper) */
      for (j = lastoneinrow+1; j <= lastcolumn; ++j)
      {
         /* if the entry is not yet fixed to 0 */
         if ( SCIPvarGetUbLocal(vars[i][j]) > 0.5 )
         {
            SCIP_Bool tightened;
            int inferInfo;

            SCIPdebugMsg(scip, " -> Fixing entry (%d,%d) to 0.\n", i, j);

            tightened = FALSE;

            /* fix variable to 0 and store position of (i,lastoneinrow+1) for conflict resolution */
            inferInfo = i * ncols + lastoneinrow + 1;
            /* correction according to Lemma 1 in the paper (second part): store (i,lastoneinrow+2) */
            if ( !infrontier )
               ++inferInfo;
            SCIP_CALL( SCIPinferBinvarCons(scip, vars[i][j], FALSE, cons, inferInfo, infeasible, &tightened) );

            /* if entry is fixed to one -> infeasible node */
            if ( *infeasible )
            {
               SCIPdebugMsg(scip, " -> Infeasible node: row %d, 1 in column %d beyond rightmost position %d\n", i, j, lastoneinrow);
               /* check if conflict analysis is applicable */
               if( SCIPisConflictAnalysisApplicable(scip) )
               {
                  int k;

                  /* conflict analysis only applicable in SOLVING stage */
                  assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING || SCIPinProbing(scip));

                  /* perform conflict analysis */
                  SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

                  /* add current bound */
                  SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i][j]) );

                  /* add bounds that result in the last one - check top left entry for packing case */
                  if ( orbitopetype == SCIP_ORBITOPETYPE_PACKING && lastones[0] == -1 )
                  {
                     assert( SCIPvarGetUbLocal(vars[0][0]) < 0.5 );
                     SCIP_CALL( SCIPaddConflictBinvar(scip, vars[0][0]) );
                  }

                  /* add bounds that result in the last one - pass through rows */
                  for (k = 1; k < i; ++k)
                  {
                     int l;
                     l = lastones[k] + 1;

                     /* if the frontier has not moved and we are not beyond the matrix boundaries */
                     if ( l <= ncols-1 && l <= k && lastones[k-1] == lastones[k] )
                     {
                        assert( SCIPvarGetUbLocal(vars[k][l]) < 0.5 );
                        SCIP_CALL( SCIPaddConflictBinvar(scip, vars[k][l]) );
                     }
                  }
                  SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
               }

               goto TERMINATE;
            }
            if ( tightened )
               ++(*nfixedvars);
         }
      }

      lastoneprevrow = lastoneinrow;
   }

   /* check whether fixing any entry to 0 results in a contradiction -> loop through rows in frontiersteps (a.k.a. gamma) */
   for (j = 0; j < nsteps; ++j)
   {
      int s;
      int lastoneinrow;

      s = frontiersteps[j];
      lastoneinrow = lastones[s];
      /* note for packing case: if we are in a frontier step then lastoneinrow >= 0 */
      assert( 0 <= lastoneinrow && lastoneinrow < ncols );

      /* if entry is not fixed */
      if ( SCIPvarGetLbLocal(vars[s][lastoneinrow]) < 0.5 && SCIPvarGetUbLocal(vars[s][lastoneinrow]) > 0.5 )
      {
         int betaprev;
         betaprev = lastoneinrow - 1;

         /* loop through rows below s */
         for (i = s+1; i < nrows; ++i)
         {
            int beta;

            assert( betaprev + 1 >= 0 );
            if ( betaprev == ncols-1 || SCIPvarGetUbLocal(vars[i][betaprev+1]) < 0.5 ) /*lint !e679*/
               beta = betaprev;
            else
               beta = betaprev + 1;
            assert( -1 <= beta && beta < ncols );

            if ( firstnonzeros[i] > beta )
            {
               SCIP_Bool tightened = FALSE;
               int inferInfo;

               /* can fix (s,lastoneinrow) (a.k.a (s,alpha)) to 1 (do not need to fix other entries to 0, since they
                * will be automatically fixed by SCIPtightenVarLb.)
                */
               assert( SCIPvarGetLbLocal(vars[s][lastoneinrow]) < 0.5 );
               SCIPdebugMsg(scip, " -> Fixing entry (%d,%d) to 1.\n", s, lastoneinrow);

               /* store position (i,firstnonzeros[i]) */
               inferInfo = ncols * nrows + i * ncols + firstnonzeros[i];
               SCIP_CALL( SCIPinferBinvarCons(scip, vars[s][lastoneinrow], TRUE, cons, inferInfo, infeasible, &tightened) );

               assert( !(*infeasible) );
               if ( tightened )
                  ++(*nfixedvars);
               break;
            }
            betaprev = beta;
         }
      }
   }

 TERMINATE:
   SCIPfreeBufferArray(scip, &frontiersteps);
   SCIPfreeBufferArray(scip, &lastones);
   SCIPfreeBufferArray(scip, &firstnonzeros);

   return SCIP_OKAY;
}


/** Propagation conflict resolving method of propagator
 *
 *  In this function we use that the propagation method above implicitly propagates SCIs, i.e., every
 *  fixing can also be gotten via an SCI-fixing.
 *
 *  Since the storage of an integer is not enough to store the complete information about the fixing
 *  nor a complete shifted column, we have to use the linear time algorithm for SCIs.
 *
 *  The inferinfo integer is set as follows:
 *
 *  - If a shifted column is fixed to 0 and the corresponding bar does not necessarily has value 1
 *    then we fix these entries to 0 and inferinfo is i * ncols + j, where (i,j) is the leader of the
 *    bar. The SCI depends on whether i is in Gamma or not (see Lemma 1 in the paper and the comments
 *    above).
 *
 *  - If a bar has value 1 and the shifted column has one entry that is not fixed, it can be fixed to
 *    1 and inferinfo is (nrows*ncols) + i * ncols + j, where (i,j) is the leader of the bar; see
 *    Proposition 1 (2c).
 */
static
SCIP_RETCODE resolvePropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   int                   inferinfo,          /**< inference information */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation conflict resolving call */
   )
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Real** vals;
   SCIP_Real** weights;
   SCIP_VAR*** vars;
   SCIP_ORBITOPETYPE orbitopetype;
   int** cases;

   int i;
   int j;
   int nrows;
   int ncols;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( result != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->nrows > 0 );
   assert( consdata->ncols > 0 );
   assert( consdata->vars != NULL );
   assert( consdata->vals != NULL );
   assert( consdata->weights != NULL );
   assert( consdata->cases != NULL );
   assert( consdata->istrianglefixed );

   *result = SCIP_DIDNOTFIND;
   if ( ! consdata->resolveprop )
      return SCIP_OKAY;

   nrows = consdata->nrows;
   ncols = consdata->ncols;
   vars = consdata->vars;
   vals = consdata->vals;
   weights = consdata->weights;
   orbitopetype = consdata->orbitopetype;
   cases = consdata->cases;

   SCIPdebugMsg(scip, "Propagation resolution method of orbitope constraint using orbitopal fixing\n");

   /* fill table */
   for (i = 0; i < nrows; ++i)
   {
      int lastcolumn;

      /* last column considered as part of the bar: */
      lastcolumn = ncols - 1;
      if ( lastcolumn > i )
         lastcolumn = i;
      for (j = 0; j <= lastcolumn; ++j)
      {
         /* if the variable was fixed to zero at conflict time */
         if ( SCIPgetVarUbAtIndex(scip, vars[i][j], bdchgidx, FALSE) < 0.5 )
            vals[i][j] = 0.0;
         else
         {
            /* if the variable was fixed to one at conflict time */
            if ( SCIPgetVarLbAtIndex(scip, vars[i][j], bdchgidx, FALSE) > 0.5 )
               vals[i][j] = 2.0;
            else
               vals[i][j] = 1.0;
         }
      }
   }

#ifdef PRINT_MATRIX
   SCIPdebugMsg(scip, "Matrix:\n");
   printMatrix(scip, consdata);
#endif

   /* computation of table: this now minimizes the value of the shifted column */
   assert( consdata->istrianglefixed );
   computeSCTable(scip, nrows, ncols, weights, cases, vals);

   /* if we fixed variables in the bar to zero */
   assert( inferinfo >= 0 && inferinfo < 2 * nrows * ncols );
   if ( inferinfo < nrows * ncols )
   {
      int p1;
      int p2;
#ifdef SCIP_DEBUG
      char str[SCIP_MAXSTRLEN];
      char tmpstr[SCIP_MAXSTRLEN];
#endif

      i = (int) (inferinfo / ncols);
      j = inferinfo % ncols;
      assert( 0 <= i && i < nrows );
      assert( 0 <= j && j < ncols );

      /* find SCI with value 0 */
      assert( weights[i-1][j-1] < 0.5 );

      SCIPdebugMsg(scip, " -> reason for x[%d][%d] = ... = x[%d][%d] = 0 was the following SC:\n", i, j, i, MIN(i,ncols));
#ifdef SCIP_DEBUG
      str[0] = '\0';
#endif

      p1 = i-1;
      p2 = j-1;
      do
      {
         assert( cases[p1][p2] != -1 );
         assert( p1 >= 0 && p1 < i );
         assert( p2 >= 0 && p2 < j );

         /* if case 1 */
         if ( cases[p1][p2] == 1 )
            --p2;   /* decrease column */
         else
         {
            /* case 2 or 3: */
            assert( cases[p1][p2] == 2 || cases[p1][p2] == 3 );
            assert( SCIPgetVarUbAtIndex(scip, vars[p1][p2], bdchgidx, FALSE) < 0.5 );
            SCIP_CALL( SCIPaddConflictUb(scip, vars[p1][p2], bdchgidx) );
            *result = SCIP_SUCCESS;

#ifdef SCIP_DEBUG
            (void) SCIPsnprintf(tmpstr, SCIP_MAXSTRLEN, " (%d,%d)", p1, p2);
            (void) strncat(str, tmpstr, SCIP_MAXSTRLEN-1);
#endif

            if ( cases[p1][p2] == 3 )
               break;
         }
         --p1;  /* decrease row */
      }
      while ( p1 >= 0 );   /* should always be true, i.e. the break should end the loop */
      assert( cases[p1][p2] == 3 );

#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, "%s\n", str);
#endif
   }
   else
   {
      int k;
      int p1;
      int p2;
#ifndef NDEBUG
      int pos1;
      int pos2;
#endif
#ifdef SCIP_DEBUG
      char str[SCIP_MAXSTRLEN];
      char tmpstr[SCIP_MAXSTRLEN];
#endif

      /* if we fixed a variable in the SC to 1 */
      inferinfo -= nrows * ncols;
      i = (int) inferinfo / ncols;
      j = inferinfo % ncols;
      assert( 0 <= i && i < nrows );
      assert( 0 <= j && j < ncols );

      /* In rare cases it might happen that we fixed a variable to 1, but the node later becomes infeasible by globally
       * fixing variables to 0. In this case, it might happen that we find a SC with value 0 instead of 1. We then
       * cannot use this SC to repropagate (and do not know how to reconstruct the original reasoning). */
      if ( weights[i-1][j-1] > 0.5 && weights[i-1][j-1] < 1.5 )
      {
         SCIPdebugMsg(scip, " -> reason for x[%d][%d] = 1 was the following SC:\n", i, j);
#ifdef SCIP_DEBUG
         (void) SCIPsnprintf(str, SCIP_MAXSTRLEN, "SC:");
#endif

         p1 = i-1;
         p2 = j-1;
#ifndef NDEBUG
         pos1 = -1;
         pos2 = -1;
#endif
         do
         {
            assert( cases[p1][p2] != -1 );
            assert( p1 >= 0 && p1 < i );
            assert( p2 >= 0 && p2 < j );

            /* if case 1 */
            if ( cases[p1][p2] == 1 )
               --p2;   /* decrease column */
            else
            {
               /* case 2 or 3: reason are formed by variables in SC fixed to 0 */
               assert( cases[p1][p2] == 2 || cases[p1][p2] == 3 );
               if ( SCIPgetVarUbAtIndex(scip, vars[p1][p2], bdchgidx, FALSE) < 0.5 )
               {
                  SCIP_CALL( SCIPaddConflictUb(scip, vars[p1][p2], bdchgidx) );
                  *result = SCIP_SUCCESS;

#ifdef SCIP_DEBUG
                  (void) SCIPsnprintf(tmpstr, SCIP_MAXSTRLEN, " (%d,%d)", p1, p2);
                  (void) strncat(str, tmpstr, SCIP_MAXSTRLEN-1);
#endif
               }
#ifndef NDEBUG
               else
               {
                  assert( SCIPgetVarLbAtIndex(scip, vars[p1][p2], bdchgidx, FALSE) < 0.5 );
                  assert( pos1 == -1 && pos2 == -1 );
                  pos1 = p1;
                  pos2 = p2;
               }
#endif
               if ( cases[p1][p2] == 3 )
                  break;
            }
            --p1;  /* decrease row */
         }
         while ( p1 >= 0 );   /* should always be true, i.e., the break should end the loop */
         assert( cases[p1][p2] == 3 );
         assert( pos1 >= 0 && pos2 >= 0 );

         /* distinguish partitioning/packing */
         if ( orbitopetype == SCIP_ORBITOPETYPE_PARTITIONING )
         {
            /* partitioning case */
#ifdef SCIP_DEBUG
            (void) SCIPsnprintf(tmpstr, SCIP_MAXSTRLEN, "  before bar: ");
            (void) strncat(str, tmpstr, SCIP_MAXSTRLEN-1);
#endif
            /* add variables before the bar in the partitioning case */
            for (k = 0; k < j; ++k)
            {
               assert( SCIPgetVarUbAtIndex(scip, vars[i][k], bdchgidx, FALSE) < 0.5 );
               SCIP_CALL( SCIPaddConflictUb(scip, vars[i][k], bdchgidx) );
               *result = SCIP_SUCCESS;
#ifdef SCIP_DEBUG
               (void) SCIPsnprintf(tmpstr, SCIP_MAXSTRLEN, " (%d,%d)", i, k);
               (void) strncat(str, tmpstr, SCIP_MAXSTRLEN-1);
#endif
            }

#ifdef SCIP_DEBUG
            SCIPdebugMsg(scip, "%s\n", str);
#endif
         }
         else
         {
            /* packing case */
            int lastcolumn;

            /* last column considered as part of the bar: */
            lastcolumn = ncols - 1;
            if ( lastcolumn > i )
               lastcolumn = i;

            /* search for variable in the bar that is fixed to 1 in the packing case */
            for (k = j; k <= lastcolumn; ++k)
            {
               if ( SCIPgetVarLbAtIndex(scip, vars[i][k], bdchgidx, FALSE) > 0.5 )
               {
                  SCIP_CALL( SCIPaddConflictLb(scip, vars[i][k], bdchgidx) );
                  *result = SCIP_SUCCESS;
                  SCIPdebugMsg(scip, "   and variable x[%d][%d] fixed to 1.\n", i, k);
                  break;
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** check packing/partitioning orbitope solution for feasibility */
static
SCIP_RETCODE enfopsPackingPartitioningOrbitopeSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< pointer to orbitope constraint */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real** weights;
   SCIP_Real** vals;
   int** cases;
   int nrows;
   int ncols;
   int i;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* do not enforce non-model constraints if strong dual reductions are not permitted */
   if ( !consdata->ismodelcons && !SCIPallowStrongDualReds(scip) )
      return SCIP_OKAY;

   assert( consdata->nrows > 0 );
   assert( consdata->ncols > 0 );
   assert( consdata->vals != NULL );
   assert( consdata->weights != NULL );
   assert( consdata->cases != NULL );

   /* check for upper right triangle */
   if ( ! consdata->istrianglefixed )
   {
      SCIP_Bool infeasible = FALSE;
      int nfixedvars = 0;

      SCIP_CALL( fixTriangle(scip, cons, &infeasible, &nfixedvars) );
      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if ( nfixedvars > 0 )
      {
         *result = SCIP_REDUCEDDOM;
         return SCIP_OKAY;
      }
   }

   nrows = consdata->nrows;
   ncols = consdata->ncols;
   vals = consdata->vals;
   weights = consdata->weights;
   cases = consdata->cases;

   /* get solution */
   copyValues(scip, consdata, NULL);
   SCIPdebugMsg(scip, "Enforcing (pseudo solutions) for orbitope constraint <%s>\n", SCIPconsGetName(cons));

   /* compute table */
   assert( consdata->istrianglefixed );
   computeSCTable(scip, nrows, ncols, weights, cases, vals);

   /* loop through rows */
   for (i = 1; i < nrows; ++i)
   {
      SCIP_Real bar = 0.0;
      int lastcolumn;

      lastcolumn = ncols - 1;

      /* last column considered as part of the bar: */
      if ( lastcolumn > i )
         lastcolumn = i;

      /* traverse row from right to left */
      for (j = lastcolumn; j > 0; --j)
      {
         bar += vals[i][j];
         assert( SCIPisIntegral(scip, vals[i][j]) );

         /* check whether weights[i-1][j-1] < bar  (<=> bar - weights[i-1][j-1] > 0), i.e. cut is violated) */
         if ( SCIPisGT(scip, bar - weights[i-1][j-1], 0.0) )
         {
            SCIPdebugMsg(scip, "Solution is infeasible.\n");
            *result = SCIP_INFEASIBLE;
            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}


/** check packing/partitioning orbitope solution for feasibility */
static
SCIP_RETCODE checkPackingPartitioningOrbitopeSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< pointer to orbitope constraint */
   SCIP_SOL*             sol,                /**< solution to be checked */
   SCIP_RESULT*          result,             /**< pointer to store the result of the enforcing call */
   SCIP_Bool             printreason         /**< whether reason for infeasibility should be printed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   SCIP_Real** vals;
   SCIP_Real** weights;
   int** cases;
   int nrows;
   int ncols;
   int i;
   int j;

   /* get data of constraint */
   assert( cons != 0 );
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->nrows > 0 );
   assert( consdata->ncols > 0 );
   assert( consdata->vars != NULL );
   assert( consdata->vals != NULL );
   assert( consdata->weights != NULL );
   assert( consdata->cases != NULL );

   nrows = consdata->nrows;
   ncols = consdata->ncols;
   vars = consdata->vars;
   vals = consdata->vals;
   weights  = consdata->weights;
   cases = consdata->cases;

   /* get solution */
   copyValues(scip, consdata, sol);
   SCIPdebugMsg(scip, "Checking orbitope constraint <%s> ...\n", SCIPconsGetName(cons));

   /* check upper right triangle (if not yet fixed to zero or in debug mode */
#ifdef NDEBUG
   if ( ! consdata->istrianglefixed )
#endif
   {
      int diagsize;

      /* get last row of triangle */
      diagsize = ncols;
      if ( nrows < ncols )
         diagsize = nrows;

      /* check variables */
      for (i = 0; i < diagsize; ++i)
      {
         for (j = i+1; j < ncols; ++j)
         {
            if ( ! SCIPisFeasZero(scip, vals[i][j]) )
            {
               if ( printreason )
                  SCIPinfoMessage(scip, NULL, "variable x[%d][%d] = %f on upper right nonzero.\n", i, j, vals[i][j]);
               *result = SCIP_INFEASIBLE;

               return SCIP_OKAY;
            }
         }
      }
   }

   /* compute table */
   computeSCTable(scip, nrows, ncols, weights, cases, vals);

   /* loop through rows */
   for (i = 1; i < nrows; ++i)
   {
      SCIP_Real bar;
      int lastcolumn;

      lastcolumn = ncols - 1;
      bar = 0.0;
      /* last column considered as part of the bar: */
      if ( lastcolumn > i )
         lastcolumn = i;

      /* traverse row from right to left */
      for (j = lastcolumn; j > 0; --j)
      {
         bar += vals[i][j];
         assert( SCIPisFeasIntegral(scip, vals[i][j]) );

         /* check whether weights[i-1][j-1] < bar  (<=> bar - weights[i-1][j-1] > 0), i.e. cut is violated) */
         if ( SCIPisGT(scip, bar - weights[i-1][j-1], 0.0) )
         {
            SCIPdebugMsg(scip, "Solution is infeasible.\n");
            *result = SCIP_INFEASIBLE;

            if ( printreason )
            {
               int l;
               int p1;
               int p2;

               SCIPinfoMessage(scip, NULL, "violated SCI: bar(");

               /* first output bar */
               for (l = j; l < ncols; ++l)
                  SCIPinfoMessage(scip, NULL, "<%s> (%f)", SCIPvarGetName(vars[i][l]), consdata->vals[i][l]);

               SCIPinfoMessage(scip, NULL, ")  SC(");

               /* output shifted column */
               p1 = i-1;
               p2 = j-1;
               do
               {
                  assert( cases[p1][p2] != -1 );
                  assert( p1 >= 0 && p1 < i );
                  assert( p2 >= 0 && p2 < j );

                  /* if case 1 */
                  if (cases[p1][p2] == 1)
                     --p2;   /* decrease column */
                  else
                  {
                     /* case 2 or 3: */
                     assert( cases[p1][p2] == 2 || cases[p1][p2] == 3 );
                     SCIPinfoMessage(scip, NULL, "<%s> (%f)", SCIPvarGetName(vars[p1][p2]), consdata->vals[p1][p2]);
                     if ( cases[p1][p2] == 3 )
                        break;
                  }
                  --p1;  /* decrease row */
               }
               while ( p1 >= 0 );   /* should always be true, i.e. the break should end the loop */
               assert( cases[p1][p2] == 3 );

               SCIPinfoMessage(scip, NULL, ")");
            }

            return SCIP_OKAY;
         }
      }
   }

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
      int nconsfixedvars = 0;
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

      /* get solution */
      copyValues(scip, consdata, sol);

      /* separate */
      SCIP_CALL( separateSCIs(scip, conshdlr, conss[c], consdata, &infeasible, &nconsfixedvars, &nconscuts) );
      nfixedvars += nconsfixedvars;
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

/** replace aggregated variables by active variables */
static
SCIP_RETCODE replaceAggregatedVarsOrbitopePP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to be processed */
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

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->vars != NULL );
   assert( consdata->nrows > 0 );
   assert( consdata->ncols > 0 );

   vars = consdata->vars;
   nrows = consdata->nrows;
   ncols = consdata->ncols;

   /* check whether there exists an aggregated variable in the orbitope */
   for (i = 0; i < nrows; ++i)
   {
      for (j = 0; j < ncols; ++j)
      {
         SCIP_VAR* var;
         SCIP_Bool negated;

         assert( SCIPvarGetStatus(vars[i][j]) != SCIP_VARSTATUS_MULTAGGR ); /* variables are marked as not to be multi-aggregated */

         SCIP_CALL( SCIPgetBinvarRepresentative(scip, vars[i][j], &var, &negated) );
         SCIP_UNUSED( negated );
         assert( SCIPvarIsActive(var) || SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED || SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED );
         if ( var != vars[i][j] )
         {
            SCIP_CALL( SCIPreleaseVar(scip, &vars[i][j]) );
            vars[i][j] = var;
            SCIP_CALL( SCIPcaptureVar(scip, var) );
         }
      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyOrbitopePP)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(valid != NULL);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrOrbitopePP(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** frees constraint handler */
static
SCIP_DECL_CONSFREE(consFreeOrbitopePP)
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
SCIP_DECL_CONSDELETE(consDeleteOrbitopePP)
{  /*lint --e{715}*/
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransOrbitopePP)
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
         sourcedata->orbitopetype, sourcedata->resolveprop, sourcedata->ismodelcons) );

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
SCIP_DECL_CONSSEPALP(consSepalpOrbitopePP)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Separation of packing/partitioning orbitope constraint handler <%s> for LP solution.\n",
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
SCIP_DECL_CONSSEPASOL(consSepasolOrbitopePP)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Separation of packing/partitioning orbitope constraint handler <%s> for primal solution.\n",
      SCIPconshdlrGetName(conshdlr));

   *result = SCIP_DIDNOTFIND;

   /* separate constraints */
   SCIP_CALL( separateConstraints(scip, conshdlr, conss, nconss, nusefulconss, sol, result, FALSE) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpOrbitopePP)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( result != NULL );

   /* we have a negative priority, so we should come after the integrality conshdlr */
   assert( SCIPgetNLPBranchCands(scip) == 0 );

   SCIPdebugMsg(scip, "Enforcement for packing/partitioning orbitope constraint handler <%s> for LP solution.\n",
      SCIPconshdlrGetName(conshdlr));

   *result = SCIP_FEASIBLE;

   /* separate constraints */
   SCIP_CALL( separateConstraints(scip, conshdlr, conss, nconss, nusefulconss, NULL, result, TRUE) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxOrbitopePP)
{  /*lint --e{715}*/
   assert( result != NULL );
   assert( scip != NULL );

   SCIPdebugMsg(scip,
      "Enforcement for packing/partitioning orbitope constraint handler <%s> for relaxation solution.\n",
      SCIPconshdlrGetName(conshdlr));

   *result = SCIP_FEASIBLE;

   /* separate constraints */
   SCIP_CALL( separateConstraints(scip, conshdlr, conss, nconss, nusefulconss, sol, result, TRUE) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsOrbitopePP)
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
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONS* cons;
      SCIP_CONSDATA* consdata;

      /* get data of constraint */
      cons = conss[c];
      assert( cons != 0 );
      consdata = SCIPconsGetData(cons);

      assert( consdata != NULL );

      /* do not enforce non-model constraints */
      if ( ! consdata->ismodelcons )
         continue;

      SCIP_CALL( enfopsPackingPartitioningOrbitopeSolution(scip, cons, result) );

      if ( *result == SCIP_INFEASIBLE )
         break;
   }

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckOrbitopePP)
{  /*lint --e{715}*/
   int c;
   SCIP_CONSDATA* consdata;

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

      SCIP_CALL( checkPackingPartitioningOrbitopeSolution(scip, conss[c], sol, result, printreason) );
   }

   if( *result == SCIP_FEASIBLE )
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
SCIP_DECL_CONSPROP(consPropOrbitopePP)
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
      int nfixed;

      assert( conss[c] != 0 );

      SCIPdebugMsg(scip, "Propagation of packing/partitioning orbitope constraint <%s> ...\n",
         SCIPconsGetName(conss[c]));

      SCIP_CALL( propagateCons(scip, conss[c], &infeasible, &nfixed) );
      nfixedvars += nfixed;
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
SCIP_DECL_CONSPRESOL(consPresolOrbitopePP)
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

      SCIPdebugMsg(scip, "Presolving of packing/partitioning orbitope constraint <%s> ...\n",
         SCIPconsGetName(conss[c]));

      /* first propagate */
      SCIP_CALL( propagateCons(scip, conss[c], &infeasible, &nfixed) );
      *nfixedvars += nfixed;

      if ( ! infeasible )
      {
         SCIP_CALL( checkRedundantCons(scip, conss[c], &redundant) );

         if ( redundant )
         {
            SCIPdebugMsg(scip,
               "Packing/Partitioning orbitope constraint <%s> is redundant: it does not contain active variables\n",
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
SCIP_DECL_CONSRESPROP(consRespropOrbitopePP)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( cons != NULL );
   assert( infervar != NULL );
   assert( bdchgidx != NULL );
   assert( result != NULL );

   SCIP_CALL( resolvePropagation(scip, cons, inferinfo, bdchgidx, result) );

   return SCIP_OKAY;
}


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreOrbitopePP)
{
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   for (c = 0; c < nconss; ++c)
   {
      /* replace aggregated variables by active variables */
      SCIP_CALL( replaceAggregatedVarsOrbitopePP(scip, conss[c]) );
   }
   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockOrbitopePP)
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

   SCIPdebugMsg(scip, "Locking method for packing/partitioning orbitope constraint handler\n");

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
SCIP_DECL_CONSPRINT(consPrintOrbitopePP)
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int i;
   int j;
   int nrows;
   int ncols;
   SCIP_ORBITOPETYPE orbitopetype;

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
   orbitopetype = consdata->orbitopetype;

   SCIPdebugMsg(scip, "Printing method for packing/partitioning orbitope constraint handler\n");

   switch ( orbitopetype )
   {
   case SCIP_ORBITOPETYPE_PARTITIONING:
      SCIPinfoMessage(scip, file, "partOrbitope(");
      break;
   case SCIP_ORBITOPETYPE_PACKING:
      SCIPinfoMessage(scip, file, "packOrbitope(");
      break;
   default:
      SCIPABORT();
   } /*lint !e788*/

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
SCIP_DECL_CONSCOPY(consCopyOrbitopePP)
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

   SCIPdebugMsg(scip, "Copying method for packing/partOrbitope orbitope constraint handler.\n");

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

      SCIP_CALL( SCIPcreateConsOrbitopePP(scip, cons, name, vars,
            sourcedata->orbitopetype, nrows, ncols, sourcedata->resolveprop, sourcedata->ismodelcons,
            initial, separate, enforce, check, propagate,
            local, modifiable, dynamic, removable, stickingatnode) );
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
SCIP_DECL_CONSPARSE(consParseOrbitopePP)
{  /*lint --e{715}*/
   const char* s;
   char* endptr;
   SCIP_ORBITOPETYPE orbitopetype;
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

   if( strncmp(s, "partOrbitope(", 13) == 0 )
      orbitopetype = SCIP_ORBITOPETYPE_PARTITIONING;
   else
   {
      if( strncmp(s, "packOrbitope(", 13) != 0 )
      {
         SCIPerrorMessage("Syntax error - expected \"partOrbitope\" or \"packOrbitope\": %s\n", s);
         *success = FALSE;
         return SCIP_OKAY;
      }
      orbitopetype = SCIP_ORBITOPETYPE_PACKING;
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
         ncols = j + 1;

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
   {
      SCIP_CALL( SCIPcreateConsOrbitopePP(scip, cons, name, vars, orbitopetype, nrows, ncols, TRUE, TRUE,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   for( k = nrows - 1; k >= 0; --k )
      SCIPfreeBufferArray(scip, &vars[k]);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsOrbitopePP)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );
   assert( success != NULL );
   assert( vars != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   if ( varssize < consdata->ncols * consdata->nrows )
      *success = FALSE;
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
      *success = TRUE;
   }

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsOrbitopePP)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   *nvars = consdata->ncols * consdata->nrows;
   *success = TRUE;

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for packing/partitioning orbitope constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrOrbitopePP(
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
         consEnfolpOrbitopePP, consEnfopsOrbitopePP, consCheckOrbitopePP, consLockOrbitopePP,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyOrbitopePP, consCopyOrbitopePP) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeOrbitopePP) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteOrbitopePP) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsOrbitopePP) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsOrbitopePP) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseOrbitopePP) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolOrbitopePP, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreOrbitopePP) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintOrbitopePP) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropOrbitopePP, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropOrbitopePP) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpOrbitopePP, consSepasolOrbitopePP, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransOrbitopePP) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxOrbitopePP) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/forceconscopy",
         "Whether orbitope constraints should be forced to be copied to sub SCIPs.",
         &conshdlrdata->forceconscopy, TRUE, DEFAULT_FORCECONSCOPY, NULL, NULL) );

   return SCIP_OKAY;
}


/** creates and captures a packing/partitioning orbitope constraint
 *
 *  @pre This constraint handler assumes that constraints which enforce the packing/partitioning constraints are
 *  contained in the problem. It does not implement, e.g., separation and propagation of set packing/partitioning
 *  constraints, since this would just copy large parts of the code of the setppc constraint handler.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsOrbitopePP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR***           vars,               /**< matrix of variables on which the symmetry acts */
   SCIP_ORBITOPETYPE     orbitopetype,       /**< type of orbitope constraint */
   int                   nrows,              /**< number of rows in orbitope matrix               <=> p   */
   int                   ncols,              /**< number of columns in orbitope matrix            <=> q   */
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
      SCIPerrorMessage("packing/partitioning orbitope constraint handler not found\n");
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
   SCIP_CALL( consdataCreate(scip, &consdata, vars, nrows, ncols, orbitopetype, resolveprop, ismodelcons) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a packing/partitioning orbitope constraint in its most basic variant, i. e., with all
 *  constraint flags set to their default values, which can be set afterwards using SCIPsetConsFLAGNAME()
 *
 *  @see SCIPcreateConsOrbitopePP() for the default constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicOrbitopePP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR***           vars,               /**< matrix of variables on which the symmetry acts */
   SCIP_ORBITOPETYPE     orbitopetype,       /**< type of orbitope constraint */
   int                   nrows,              /**< number of rows in orbitope matrix               <=> p   */
   int                   ncols,              /**< number of columns in orbitope matrix            <=> q   */
   SCIP_Bool             resolveprop,        /**< should propagation be resolved? */
   SCIP_Bool             ismodelcons         /**< whether the orbitope is a model constraint */
   )
{
   SCIP_CALL( SCIPcreateConsOrbitopePP(scip, cons, name, vars, orbitopetype, nrows, ncols,
         resolveprop, ismodelcons, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

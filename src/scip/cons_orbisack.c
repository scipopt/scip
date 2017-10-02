/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_orbisack.c
 * @brief  constraint handler for orbisack constraints
 * @author Christopher Hojny
 *
 * The constraint works on two vectors of variables, which are interpreted as columns of a matrix such that the first
 * column is lexicographically not smaller than the second.\
 *
 * @TODO: Currently, the copy methods of the constraint handler are deactivated, as otherwise, we would reduce the
 * effect of heuristics. The trade-off for this is that we cannot transfer the symmetry information to the sub-scips
 * of the components presolver.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "scip/cons_orbisack.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "orbisack"
#define CONSHDLR_DESC          "symmetry breaking constraint handler relying on orbisacks"
#define CONSHDLR_SEPAPRIORITY    +40100 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -1005200 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -1005200 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             5 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             5 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP
#define CONSHDLR_PRESOLTIMING      SCIP_PRESOLTIMING_EXHAUSTIVE

/* default parameters: */
#define DEFAULT_PPORBISACK                   FALSE     /** whether orbisack is of packing/partitioning type */

/* default parameters for separation routines: */
#define DEFAULT_ORBISEPARATION               FALSE     /**< whether orbisack inequalities should be separated */
#define DEFAULT_COVERSEPARATION              FALSE     /**< whether cover inequalities should be separated */
#define DEFAULT_ENFORCING                    FALSE     /**< whether we use enforcing methods of the constraint handler */
#define DEFAULT_CHECK                        FALSE     /**< whether we use check methods of the constraint handler */

/* default parameters for constraints */
#define DEFAULT_COEFFBOUND               1000000.0     /**< maximum size of coefficients in orbisack inequalities */


/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Bool             coverSeparation_;   /**< whether only cover inequalities should be separated */
   SCIP_Bool             orbiSeparation_;    /**< whether orbisack as well as cover inequalities should be separated */
   SCIP_Bool             orbiEnforcing_;     /**< whether we use enforcing methods of constraint hanlder */
   SCIP_Bool             orbiCheck_;         /**< whether we use check methods of constraint hanlder */
   SCIP_Real             coeffBound_;        /**< maximum size of coefficients in orbisack inequalities */
};

/** constraint data for orbisack constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars1_;             /**< first column of variable matrix */
   SCIP_VAR**            vars2_;             /**< second column of variable matrix */
   unsigned int          nRows_;             /**< number of rows of variable matrix */
   SCIP_Real*            vals1_;             /**< LP-solution for those variables in first column */
   SCIP_Real*            vals2_;             /**< LP-solution for those variables in second column */
   SCIP_Real             isPporbisack_;      /**< whether the orbisack is a packing/partitioning orbisack */
};


/*
 * Local methods
 */

/** frees an orbisack constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to orbisack constraint data */
   )
{
   unsigned int nrows;

   assert( consdata != NULL );
   assert( *consdata != NULL );

   nrows = (*consdata)->nRows_;
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vars2_), nrows);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vars1_), nrows);

   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vals2_), nrows);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vals1_), nrows);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** creates orbisack constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store constraint data */
   SCIP_VAR*const*       vars1,              /**< first column of variable matrix */
   SCIP_VAR*const*       vars2,              /**< second column of variable matrix */
   unsigned int          nrows,              /**< number of rows in variable matrix */
   SCIP_Bool             isPporbisack        /**< whether orbisack is of packing/partitioning type */
   )
{
   assert(consdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->vals1_, nrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->vals2_, nrows) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars1_, vars1, nrows) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars2_, vars2, nrows) );

#ifndef NDEBUG
   {
      unsigned int i;
      for (i = 0; i < nrows; ++i)
      {
         assert( SCIPvarIsBinary(vars1[i]) );
         assert( SCIPvarIsBinary(vars2[i]) );
      }
   }
#endif

   (*consdata)->nRows_ = nrows;
   (*consdata)->isPporbisack_ = isPporbisack;

   return SCIP_OKAY;
}


/** get the current lp solution from SCIP solution @p sol
 *
 *  We store the solution in the data of the constraint.
 */
static
SCIP_RETCODE getValues(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_SOL*             sol,                /**< solution (may be NULL) */
   const SCIP_CONSDATA*  consdata            /**< constraint data */
   )
{
   unsigned int nRows;
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   SCIP_Real* vals1;
   SCIP_Real* vals2;

   assert( scip != NULL );
   assert( consdata != NULL );
   assert( consdata->nRows_ > 0 );
   assert( consdata->vars1_ != NULL );
   assert( consdata->vars2_ != NULL );
   assert( consdata->vals1_ != NULL );
   assert( consdata->vals2_ != NULL );

   nRows = consdata->nRows_;
   vars1 = consdata->vars1_;
   vars2 = consdata->vars2_;
   vals1 = consdata->vals1_;
   vals2 = consdata->vals2_;

   SCIP_CALL( SCIPgetSolVals(scip, sol, nRows, vars1, vals1) );
   SCIP_CALL( SCIPgetSolVals(scip, sol, nRows, vars2, vals2) );

   return SCIP_OKAY;
}


/** generate initial LP cut
 *
 *  We generate the inequality of the orbisack on the elements of the first row, i.e.,
 *  the inequality \f$-x_{1,1} + x_{1,2} \leq 0\f$.
 */
static
SCIP_RETCODE initLP(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool*            infeasible          /**< whether we detected infeasibility */
   )
{
   SCIP_CONSDATA* consdata;
   unsigned int nRows;
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   SCIP_VAR* tmpVars[2];
   SCIP_ROW* row;
   SCIP_Bool isPporbisack;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );

   *infeasible = FALSE;

   consdata = SCIPconsGetData(cons);
   assert( consdata != 0 );
   assert( consdata->vars1_ != NULL );
   assert( consdata->vars2_ != NULL );
   assert( consdata->vals1_ != NULL );
   assert( consdata->vals2_ != NULL );

   nRows = consdata->nRows_;
   vars1 = consdata->vars1_;
   vars2 = consdata->vars2_;

   /* avoid stupid problems */
   if ( nRows == 0 )
      return SCIP_OKAY;

   tmpVars[0] = vars1[0];
   tmpVars[1] = vars2[0];

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), "orbisack0#0", -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, tmpVars[0], -1.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, tmpVars[1], 1.0) );

   SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, infeasible) );

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif
   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   isPporbisack = consdata->isPporbisack_;

   if ( isPporbisack )
   {
      SCIP_VAR** tempVars;
      SCIP_Real* tempVals;
      unsigned int i;
      char name[SCIP_MAXSTRLEN];

      SCIP_CALL( SCIPallocBufferArray(scip, &tempVars, nRows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &tempVals, nRows) );

      tempVals[0] = 1.0;

      for (i = 0; i < nRows - 1; ++i)
      {
         tempVars[0] = vars2[i];

         if ( i > 0 )
         {
            tempVars[i] = vars1[i - 1];
            tempVals[i] = -1.0;
         }

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pporbisack#%d", i);
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), name, -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
         SCIP_CALL( SCIPaddVarsToRow(scip, row, i + 1, tempVars, tempVals) );

         SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, infeasible) );

#ifdef SCIP_DEBUG
         SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif
         SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }

      SCIPfreeBufferArray(scip, &tempVals);
      SCIPfreeBufferArray(scip, &tempVars);
   }

   return SCIP_OKAY;
}


/** add orbisack cover inequality */
static
SCIP_RETCODE addOrbisackCover(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   const SCIP_CONSDATA*  consdata,           /**< constraint data */
   SCIP_Real*            coeffs1,            /**< coefficients of the variables of the first column of the inequality to be added */
   SCIP_Real*            coeffs2,            /**< coefficients of the variables of the second column of the inequality to be added */
   SCIP_Real             rhs,                /**< right-hand side of inequality to be added */
   SCIP_Bool*            infeasible          /**< whether we detected infeasibility */
   )
{
   SCIP_ROW* row;
   unsigned int i;

   assert( scip != NULL );
   assert( consdata != NULL );
   assert( consdata->nRows_ > 0 );
   assert( consdata->vars1_ != NULL );
   assert( consdata->vars2_ != NULL );
   assert( consdata->vals1_ != NULL );
   assert( consdata->vals2_ != NULL );
   assert( coeffs1 != NULL );
   assert( coeffs2 != NULL );
   assert( infeasible != NULL );

   *infeasible = FALSE;

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), "orbisackcover", -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
   for (i = 0; i < consdata->nRows_; ++i)
   {
      SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->vars1_[i], coeffs1[i]) );
      SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->vars2_[i], coeffs2[i]) );
   }
   SCIP_CALL( SCIPflushRowExtensions(scip, row) );

   SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, infeasible) );

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif
   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   return SCIP_OKAY;
}


/** Separate orbisack cover inequalities
 *
 *  We currently do NOT enter cuts into the pool.
 *
 *  We iterate over the nRow-many cover inequalities which are potentially
 *  maximal w.r.t. their violation.
 */
static
SCIP_RETCODE separateOrbisackCovers(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   const SCIP_CONSDATA*  consdata,           /**< constraint data */
   unsigned int*         nGen,               /**< number of separated covers */
   SCIP_Bool*            infeasible          /**< whether we detected infeasibility */
   )
{
   SCIP_Real rhs = 0.0;
   SCIP_Real lhs = 0.0;
   unsigned int nRows;
   SCIP_Real* vals1;
   SCIP_Real* vals2;
   unsigned int i;
   SCIP_Real* coeff1;
   SCIP_Real* coeff2;

   assert( scip != NULL );
   assert( consdata != NULL );
   assert( consdata->nRows_ > 0 );
   assert( consdata->vars1_ != NULL );
   assert( consdata->vars2_ != NULL );
   assert( consdata->vals1_ != NULL );
   assert( consdata->vals2_ != NULL );
   assert( infeasible != NULL );
   assert( nGen != NULL );

   *infeasible = FALSE;
   *nGen = 0;

   nRows = consdata->nRows_;
   vals1 = consdata->vals1_;
   vals2 = consdata->vals2_;

   /* allocate memory for inequality coefficients */
   SCIP_CALL( SCIPallocBufferArray(scip, &coeff1, nRows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coeff2, nRows) );

   /* initialize coefficient matrix */
   for (i = 0; i < nRows; ++i)
   {
      coeff1[i] = 0.0;
      coeff2[i] = 0.0;
   }

   /* detect violated covers */
   for (i = 0; i < nRows; ++i)
   {
      /* cover inequality is violated */
      if ( SCIPisEfficacious(scip, -vals1[i] + vals2[i] + lhs - rhs) )
      {
         /* set coefficients for inequality */
         coeff1[i] = -1.0;
         coeff2[i] = 1.0;
         SCIP_CALL( addOrbisackCover(scip, cons, consdata, coeff1, coeff2, rhs, infeasible) );
         ++(*nGen);
         if ( *infeasible )
            break;

         /* reset coefficients for next inequality */
         coeff1[i] = 0.0;
         coeff2[i] = 0.0;
      }

      /* add argmax( 1 - vals[i][0], vals[i][1] ) as coefficient */
      if ( SCIPisEfficacious(scip, 1.0 - vals1[i] - vals2[i]) )
      {
         coeff1[i] = -1.0;
         lhs = lhs - vals1[i];
      }
      else
      {
         coeff2[i] = 1.0;
         rhs += 1.0;
         lhs = lhs + vals2[i];
      }
   }

   /* free coefficient matrix */
   SCIPfreeBufferArray(scip, &coeff2);
   SCIPfreeBufferArray(scip, &coeff1);

   return SCIP_OKAY;
}


/** add orbisack inequality */
static
SCIP_RETCODE addOrbisackInequality(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   const SCIP_CONSDATA*  consdata,           /**< constraint data */
   SCIP_Real*            coeffs1,            /**< first column of coefficient matrix of inequality to be added */
   SCIP_Real*            coeffs2,            /**< second column of coefficient matrix of inequality to be added */
   SCIP_Real             rhs,                /**< right-hand side of inequality to be added */
   SCIP_Bool*            infeasible          /**< whether we detected infeasibility */
   )
{
   SCIP_ROW* row;
   unsigned int i;

   assert( scip != NULL );
   assert( consdata != NULL );
   assert( consdata->nRows_ > 0 );
   assert( consdata->vars1_ != NULL );
   assert( consdata->vars2_ != NULL );
   assert( consdata->vals1_ != NULL );
   assert( consdata->vals2_ != NULL );
   assert( coeffs1 != NULL );
   assert( coeffs2 != NULL );
   assert( infeasible != NULL );

   *infeasible = FALSE;

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), "orbisack", -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

   for (i = 0; i < consdata->nRows_; ++i)
   {
      SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->vars1_[i], coeffs1[i]) );
      SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->vars2_[i], coeffs2[i]) );
   }
   SCIP_CALL( SCIPflushRowExtensions(scip, row) );

   SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, infeasible) );

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif
   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   return SCIP_OKAY;
}


/** Separate orbisack inequalities
 *
 *  We currently do NOT enter cuts into the pool.
 *
 *  We stop if we checked for each possible basement row, whether a cut
 *  could be added. If the coefficients grow too large, we start
 *  separating cover inequalities.
 *
 *  We implement the separation algorithm for orbisacks described in@n
 *  A. Loos. Describing Orbitopes by Linear Inequalities and Projection Based Tools.
 *  PhD thesis, Otto-von-Guericke-Universit\"at Magdeburg, 2010.
 */
static
SCIP_RETCODE separateOrbisack(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            conss,              /**< constraint */
   const SCIP_CONSDATA*  consdata,           /**< constraint data */
   unsigned int*         nGen,               /**< number of generated cuts */
   SCIP_Bool*            infeasible,         /**< whether we detected infeasibility */
   SCIP_Bool             coverSeparation,    /**< whether we separate cover inequalities */
   SCIP_Real             coeffBound          /**< maximum size of coefficients in orbisack inequalities */
   )
{
   unsigned int nRows;
   SCIP_Real* vals1;
   SCIP_Real* vals2;
   SCIP_Real* coeff1;
   SCIP_Real* coeff2;
   unsigned int i;
   unsigned int basement;
   SCIP_Real rhs;
   SCIP_Real lhs;
   SCIP_Real valueA;
   SCIP_Real valueB;
   SCIP_Real valueC;

   assert( scip != NULL );
   assert( consdata != NULL );
   assert( consdata->nRows_ > 0 );
   assert( consdata->vars1_ != NULL );
   assert( consdata->vars2_ != NULL );
   assert( consdata->vals1_ != NULL );
   assert( consdata->vals2_ != NULL );
   assert( nGen != NULL );
   assert( infeasible != NULL );
   assert( coeffBound >= 0.0 );

   *infeasible = FALSE;
   *nGen = 0;

   nRows = consdata->nRows_;
   vals1 = consdata->vals1_;
   vals2 = consdata->vals2_;

   /* avoid trivial problems */
   if ( nRows < 2 )
      return SCIP_OKAY;

   /* allocate memory for inequality coefficients */
   SCIP_CALL( SCIPallocBufferArray(scip, &coeff1, nRows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coeff2, nRows) );

   /* initialize coefficient matrix row 0 */
   coeff1[0] = -1.0;
   coeff2[0] = 1.0;
   for (i = 2; i < nRows; ++i)
   {
      coeff1[i] = 0.0;
      coeff2[i] = 0.0;
   }

   /* initialize right-hand side and left-hand side (lhs for row 0) */
   rhs = 0.0;
   lhs = - vals1[0] + vals2[0];

   /* basement row of orbisack */
   basement = 1;

   /* update value of left-hand side and coefficients for basement row = 1 */
   lhs += - vals1[1] + vals2[1];
   coeff1[1] = -1.0;
   coeff2[1] = 1.0;

   while ( basement < nRows )
   {
      /* to avoid numerical troubles, we bound the size of coefficients and rhs */
      if ( rhs > coeffBound || -coeff1[0] > coeffBound || coeff2[0] > coeffBound )
      {
         /* avoid separating cover inequalities twice */
         if ( ! coverSeparation )
         {
            unsigned int ncuts;
            SCIP_CALL( separateOrbisackCovers(scip, conss, consdata, &ncuts, infeasible) );
            *nGen += ncuts;
         }
         break;
      }

      /* if current inequality is violated */
      if ( SCIPisEfficacious(scip, lhs - rhs) )
      {
         SCIP_CALL( addOrbisackInequality(scip, conss, consdata, coeff1, coeff2, rhs, infeasible) );
         ++(*nGen);
         if ( *infeasible )
            break;
      }

      /* if we cannot extend the inequality by a further row */
      if ( basement >= nRows - 1 )
         break;

      valueA = lhs + vals1[basement] - vals1[basement + 1] + vals2[basement + 1] - rhs - 1.0; /*lint !e679, !e834*/
      valueB = lhs - vals2[basement] - vals1[basement + 1] + vals2[basement + 1] - rhs; /*lint !e679, !e834*/
      valueC = 2.0 * lhs + vals1[basement] - vals2[basement] - vals1[basement + 1] + vals2[basement + 1] - 2.0 * rhs; /*lint !e679, !e834*/

      /* update inequality */
      if ( valueA >= valueB && valueA >= valueC )
      {
         ++rhs;
         coeff1[basement] = 0.0;
         lhs += vals1[basement++];
         coeff1[basement] = -1.0;
         coeff2[basement] = 1.0;
         lhs += - vals1[basement] + vals2[basement];
      }
      else if ( valueB >= valueA && valueB >= valueC )
      {
         coeff2[basement] = 0.0;
         lhs -= vals2[basement++];
         coeff1[basement] = -1.0;
         coeff2[basement] = 1.0;
         lhs += - vals1[basement] + vals2[basement];
      }
      else
      {
         rhs *= 2.0;
         lhs = 0.0;
         for (i = 0; i < basement; ++i)
         {
            coeff1[i] = 2.0 * coeff1[i];
            coeff2[i] = 2.0 * coeff2[i];
            lhs += coeff1[i] * vals1[i] + coeff2[i] * vals2[i];
         }
         coeff1[basement] = -1.0;
         coeff2[basement] = 1.0;
         lhs -= vals1[basement];
         lhs += vals2[basement++];
         coeff1[basement] = -1.0;
         coeff2[basement] = 1.0;
         lhs -= vals1[basement];
         lhs += vals2[basement];
      }
   }

   /* free allocated memory */
   SCIPfreeBufferArray(scip, &coeff2);
   SCIPfreeBufferArray(scip, &coeff1);

   return SCIP_OKAY;
}


/** propagation */
static
SCIP_RETCODE propVariables(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint to be propagated */
   SCIP_Bool*            infeasible,         /**< whether it was detected that the node is infeasible */
   SCIP_Bool*            found,              /**< whether a new propagation could be found */
   unsigned int*         nGen                /**< number of generated bound strengthenings */
   )
{
   SCIP_CONSDATA* consdata;
   unsigned int nRows;
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   SCIP_Bool tightened;
   unsigned int* solu1;
   unsigned int* solu2;
   unsigned int i, r;
   SCIP_VAR* var;
   SCIP_Real ub;
   SCIP_Real lb;
   SCIP_VAR* var1;
   SCIP_VAR* var2;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );
   assert( nGen != NULL );
   assert( found != NULL );

   SCIPdebugMessage("Propagating variables of constraint <%s>.\n", SCIPconsGetName(cons));

   *nGen = 0;
   *infeasible = FALSE;
   *found = FALSE;

   /* get data of constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->vars1_ != NULL );
   assert( consdata->vars2_ != NULL );
   assert( consdata->nRows_ > 0 );

   nRows = consdata->nRows_;
   vars1 = consdata->vars1_;
   vars2 = consdata->vars2_;
   tightened = FALSE;

   if ( nRows < 2 )
      return SCIP_OKAY;

   /* determine current solution */
   SCIP_CALL( SCIPallocBufferArray(scip, &solu1, (int) nRows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solu2, (int) nRows) );

   for (i = 0; i < nRows; ++i)
   {
      /* determine value in first column */
      var = vars1[i];

      /* possibly get non-negated variables (variables might have been aggregated in between) */
      if ( SCIPvarIsNegated(var) )
         var = SCIPvarGetNegatedVar(var);
      assert( ! SCIPvarIsNegated(var) );
      assert( SCIPvarIsBinary(var) );

      /* get local upper and lower bound on variable */
      ub = SCIPvarGetUbLocal(var);
      lb = SCIPvarGetLbLocal(var);

      /* if variable is fixed to 1 -> solu[i][j] = 1,
       * if it is fixed to 0 -> solu[i][j] = 0,
       * else, -> solu[i][j] = 2 */
      if ( lb > 0.5 )
         solu1[i] = 1;
      else if (ub < 0.5)
         solu1[i] = 0;
      else
         solu1[i] = 2;

      /* determine value in second column */
      var = vars2[i];

      /* possibly get non-negated variables (variables might have been aggregated in between) */
      if ( SCIPvarIsNegated(var) )
         var = SCIPvarGetNegatedVar(var);
      assert( ! SCIPvarIsNegated(var) );
      assert( SCIPvarIsBinary(var) );

      /* get local upper and lower bound on variable */
      ub = SCIPvarGetUbLocal(var);
      lb = SCIPvarGetLbLocal(var);

      /* if variable is fixed to 1 -> solu[i][j] = 1,
       * if it is fixed to 0 -> solu[i][j] = 0,
       * else, -> solu[i][j] = 2 */
      if ( lb > 0.5 )
         solu2[i] = 1;
      else if (ub < 0.5)
         solu2[i] = 0;
      else
         solu2[i] = 2;
   }

   /* loop through all variables */
   for (i = 0; i < nRows; ++i)
   {
      /* get variables of first and second column */
      var1 = vars1[i];
      var2 = vars2[i];
      assert( var1 != NULL );
      assert( var2 != NULL );

      /* possibly get non-negated variables (variables might have been aggregated in between) */
      if ( SCIPvarIsNegated(var1) )
         var1 = SCIPvarGetNegatedVar(var1);
      assert( ! SCIPvarIsNegated(var1) );

      if ( SCIPvarIsNegated(var2) )
         var2 = SCIPvarGetNegatedVar(var2);
      assert( ! SCIPvarIsNegated(var2) );

      /* if variable in first column is fixed to 0 and variable in second column is fixed to 1 */
      if ( SCIPvarGetUbLocal(var1) < 0.5 && SCIPvarGetLbLocal(var2) > 0.5 )
      {
         SCIP_Bool noCritical = TRUE;

         SCIPdebugMessage("Check variable pair (%u,0) and (%u,1).\n", i, i);

         /* check whether there is a critical row above row i, otherwise the solution is infeasible
          * if there is a row without fixing (2) above the current row, we cannot obtain a result */
         for (r = 0; r < i; ++r)
         {
            if ( (solu1[r] == 1 && solu2[r] == 0) || solu1[r] == 2 || solu2[r] == 2 )
            {
               noCritical = FALSE;
               break;
            }
         }

         if ( noCritical )
         {
            SCIPdebugMessage(" -> node infeasible (row was fixed to 0,1 but there was no critical row above).\n");

            /* perform conflict analysis */
            if ( SCIPisConflictAnalysisApplicable(scip) )
            {
#if ( SCIP_VERSION > 321 || (SCIP_VERSION == 321 && SCIP_SUBVERSION >= 2 ))
               SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );
#else
               SCIP_CALL( SCIPinitConflictAnalysis(scip) );
#endif

               for (r = 0; r <= i; ++r)
               {
                  SCIP_CALL( SCIPaddConflictBinvar(scip, vars1[r]) );
                  SCIP_CALL( SCIPaddConflictBinvar(scip, vars2[r]) );
               }

               SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

               *infeasible = TRUE;

               /* free current solution */
               SCIPfreeBufferArray(scip, &solu2);
               SCIPfreeBufferArray(scip, &solu1);

               return SCIP_OKAY;
            }
         }
      }
      /* if variable in the first column is fixed to 0 and the variable in the second column is free */
      else if ( SCIPvarGetUbLocal(var1) < 0.5 && SCIPvarGetUbLocal(var2) > 0.5 )
      {
         SCIP_Bool allConstant = TRUE;

         SCIPdebugMessage("Check variable pair (%u,0) and (%u,1).\n", i, i);

         /* Check whether all rows above row i are constant. In this case, the variable in the second */
         /* column can be fixed to 0. If an entry above row i is unfixed or a row is not constant, we cannot */
         /* fix the second entry in row i. */
         for (r = 0; r < i; ++r)
         {
            if ( solu1[r] == 2 || solu2[r] == 2 || solu1[r] != solu2[r] )
            {
               allConstant = FALSE;
               break;
            }
         }

         /* fix variable in the second column to 0 */
         if ( allConstant )
         {
            assert( SCIPvarGetLbLocal(var2) < 0.5 );
            SCIP_CALL( SCIPinferVarUbCons(scip, var2, 0.0, cons, i, FALSE, infeasible, &tightened) ); /*lint !e713*/
            assert( ! *infeasible );

            *found = *found || tightened;
            if ( tightened )
               ++(*nGen);
         }
      }
   }

   /* free current solution */
   SCIPfreeBufferArray(scip, &solu2);
   SCIPfreeBufferArray(scip, &solu1);

   return SCIP_OKAY;
}


/** separate orbisack and cover inequalities */
static
SCIP_RETCODE separateInequalities(
   SCIP*                 scip,               /**< pointer to scip */
   SCIP_RESULT*          result,             /**< pointer to store the result of separation */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool infeasible = FALSE;
   unsigned int nGen1 = 0;
   unsigned int nGen2 = 0;

   assert( scip != NULL );
   assert( result != NULL );

   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert( conshdlrdata != NULL );

   if ( conshdlrdata->orbiSeparation_ )
   {
      SCIP_CALL( separateOrbisack(scip, cons, consdata, &nGen1, &infeasible, FALSE, conshdlrdata->coeffBound_) );
   }

   if ( ! infeasible && conshdlrdata->coverSeparation_ )
   {
      SCIP_CALL( separateOrbisackCovers(scip, cons, consdata, &nGen2, &infeasible) );
   }

   if ( infeasible )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   if ( nGen1 + nGen2 > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


/*--------------------------------------------------------------------------------------------
 *--------------------------------- SCIP functions -------------------------------------------
 *--------------------------------------------------------------------------------------------*/

#if 0
/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyOrbisack)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrOrbisack(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}
#endif


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteOrbisack)
{
   /*lint --e{715}*/
   assert( scip != 0 );
   assert( conshdlr != 0 );
   assert( consdata != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}


/** frees constraint handler */
static
SCIP_DECL_CONSFREE(consFreeOrbisack)
{   /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;

   assert( scip != 0 );
   assert( conshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransOrbisack)
{
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* consdata = NULL;
   unsigned int nRows;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );
   assert( targetcons != NULL );

   SCIPdebugMessage("Transforming constraint.\n");

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );
   assert( sourcedata->nRows_ > 0 );
   assert( sourcedata->vars1_ != NULL );
   assert( sourcedata->vars2_ != NULL );
   assert( sourcedata->vals1_ != NULL );
   assert( sourcedata->vals2_ != NULL );

   /* create transformed constraint data (copy data where necessary) */
   nRows = sourcedata->nRows_;

   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->nRows_ = nRows;
   consdata->isPporbisack_ = sourcedata->isPporbisack_;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars1_, nRows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars2_, nRows) );

   SCIP_CALL( SCIPgetTransformedVars(scip, nRows, sourcedata->vars1_, consdata->vars1_) );
   SCIP_CALL( SCIPgetTransformedVars(scip, nRows, sourcedata->vars2_, consdata->vars2_) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vals1_, nRows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vals2_, nRows) );

   /* create transformed constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, consdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
         SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpOrbisack)
{
   int c;
#if ( SCIP_VERSION < 321 || ( SCIP_VERSION == 321 && SCIP_SUBVERSION == 0 ) )
   SCIP_Bool infeasible = FALSE;
#else
   assert( infeasible != NULL );
   *infeasible = FALSE;
#endif

   assert( scip != 0 );
   assert( conshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      /* get data of constraint */
      assert( conss[c] != 0 );

      SCIPdebugMessage("Generating initial orbisack cut for constraint <%s> ...\n", SCIPconsGetName(conss[c]));

#if ( SCIP_VERSION < 321 || ( SCIP_VERSION == 321 && SCIP_SUBVERSION == 0 ) )
      SCIP_CALL( initLP(scip, conss[c], &infeasible) );
      if ( infeasible )
         break;
#else
      SCIP_CALL( initLP(scip, conss[c], infeasible) );
      if ( *infeasible )
         break;
#endif
      SCIPdebugMessage("Generated initial orbisack cut.\n");
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solution */
static
SCIP_DECL_CONSSEPALP(consSepalpOrbisack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMessage("Separation method for orbisack constraints\n");

   *result = SCIP_DIDNOTRUN;

   /* if solution is not integer */
   if ( SCIPgetNLPBranchCands(scip) > 0 )
   {
      *result = SCIP_DIDNOTFIND;

      /* loop through constraints */
      for (c = 0; c < nconss; ++c)
      {
         /* get data of constraint */
         assert( conss[c] != NULL );
         consdata = SCIPconsGetData(conss[c]);

         /* get solution */
         SCIP_CALL( getValues(scip, NULL, consdata) );

         SCIPdebugMessage("Separating orbisack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

         SCIP_CALL( separateInequalities(scip, result, conss[c], consdata) );
         if ( *result == SCIP_CUTOFF )
            break;
      }
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solution */
static
SCIP_DECL_CONSSEPASOL(consSepasolOrbisack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMessage("Separation method for orbisack constraints\n");

   *result = SCIP_DIDNOTFIND;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      /* get data of constraint */
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);

      /* get solution */
      SCIP_CALL( getValues(scip, sol, consdata) );

      SCIPdebugMessage("Separating orbisack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      SCIP_CALL( separateInequalities(scip, result, conss[c], consdata) );
      if ( *result == SCIP_CUTOFF )
         break;
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions
 *
 *  @pre It is assumed that the solution is integral (this can be ensured by appropriate priorities).
 */
static
SCIP_DECL_CONSENFOLP(consEnfolpOrbisack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   unsigned int nGen = 0;
   SCIP_Bool infeasible = FALSE;
   int c;

   assert( scip != 0 );
   assert( conshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != 0 );

   SCIPdebugMessage("Enfolp method for orbisack constraints\n");

   /* we have a negative priority, so we should come after the integrality conshdlr. */
   assert( SCIPgetNLPBranchCands(scip) == 0 );

   *result = SCIP_FEASIBLE;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      /* get data of constraint */
      assert( conss[c] != 0 );
      consdata = SCIPconsGetData(conss[c]);

      /* get solution */
      SCIP_CALL( getValues(scip, NULL, consdata) );

      SCIPdebugMessage("Enforcing orbisack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      /* Separate only cover inequalities to ensure that enforcing works correctly. */
      /* Otherwise, it may happen that infeasible solutions cannot be detected since */
      /* we bound the size of the coefficients for the orbisack inequalities. */
      SCIP_CALL( separateOrbisackCovers(scip, conss[c], consdata, &nGen, &infeasible) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      SCIPdebugMessage("Generated orbisack inequalities for <%s>: %u\n", SCIPconsGetName(conss[c]), nGen);

      if ( nGen > 0 )
         *result = SCIP_SEPARATED;
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsOrbisack)
{
   SCIP_Bool terminated = FALSE;
   SCIP_CONSDATA* consdata;
   unsigned int nRows;
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   SCIP_Real* solu1;
   SCIP_Real* solu2;
   unsigned int i;
   SCIP_Real val;
   int c;

   /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMessage("Enforcing method for orbisack constraints (pseudo solutions) ...\n");

   *result = SCIP_FEASIBLE;

   if ( objinfeasible || solinfeasible )
      return SCIP_OKAY;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      /* get data of constraint */
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL);
      assert( consdata->nRows_ > 0 );
      assert( consdata->vars1_ != NULL );
      assert( consdata->vars2_ != NULL );
      nRows = consdata->nRows_;
      vars1 = consdata->vars1_;
      vars2 = consdata->vars2_;

      /* determine solution */
      SCIP_CALL( SCIPallocBufferArray(scip, &solu1, (int) nRows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &solu2, (int) nRows) );

      for (i = 0; i < nRows; ++i)
      {
         /* get value of first column */
         val = SCIPgetSolVal(scip, NULL, vars1[i]);
         assert( SCIPisFeasIntegral(scip, val) );

         /* if variable is fixes to 1 -> solu[i][j] = 1, else = 0 */
         if ( val > 0.5 )
            solu1[i] = 1;
         else
            solu1[i] = 0;

         /* get value of first column */
         val = SCIPgetSolVal(scip, NULL, vars2[i]);
         assert( SCIPisFeasIntegral(scip, val) );

         /* if variable is fixes to 1 -> solu[i][j] = 1, else = 0 */
         if ( val > 0.5 )
            solu2[i] = 1;
         else
            solu2[i] = 0;
      }

      /* check whether columns of solution are not lexicographically maximal */
      for (i = 0; i < nRows; ++i)
      {
         /* if row i is constant */
         if ( SCIPisEQ(scip, solu1[i],solu2[i]) )
            continue;

         /* if first non-constant row is [1,0]: feasible */
         if ( SCIPisEQ(scip, solu1[i], 1.0) )
            break;
         else /* infeasible */
         {
            SCIPdebugMessage("Solution is infeasible.\n");
            *result = SCIP_INFEASIBLE;
            terminated = TRUE;
            break;
         }
      }

      /* free buffers */
      SCIPfreeBufferArrayNull(scip, &solu2);
      SCIPfreeBufferArrayNull(scip, &solu1);

      if ( terminated )
         break;
   }

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckOrbisack)
{
   SCIP_Bool terminated = FALSE;
   SCIP_CONSDATA* consdata;
   unsigned int nRows;
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   SCIP_Real* solu1;
   SCIP_Real* solu2;
   unsigned int i;
   SCIP_Real val;
   int c;

   /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_FEASIBLE;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      /* get data of constraint */
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL);
      assert( consdata->nRows_ > 0 );
      assert( consdata->vars1_ != NULL );
      assert( consdata->vars2_ != NULL );

      SCIPdebugMessage("Check method for orbisack constraint <%s> (%u rows) ...\n", SCIPconsGetName(conss[c]), consdata->nRows_);

      nRows = consdata->nRows_;
      vars1 = consdata->vars1_;
      vars2 = consdata->vars2_;

      /* determine solution */
      SCIP_CALL( SCIPallocBufferArray(scip, &solu1, (int) nRows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &solu2, (int) nRows) );

      for (i = 0; i < nRows; ++i)
      {
         /* get value of first column */
         val = SCIPgetSolVal(scip, sol, vars1[i]);
         assert( SCIPisFeasIntegral(scip, val) );

         /* if variable is fixes to 1 -> solu[i][j] = 1, else = 0 */
         if ( val > 0.5 )
            solu1[i] = 1;
         else
            solu1[i] = 0;

         /* get value of second column */
         val = SCIPgetSolVal(scip, sol, vars2[i]);
         assert( SCIPisFeasIntegral(scip, val) );

         /* if variable is fixes to 1 -> solu[i][j] = 1, else = 0 */
         if ( val > 0.5 )
            solu2[i] = 1;
         else
            solu2[i] = 0;
      }

      /* check whether columns of solution are not lexicographically maximal */
      for (i = 0; i < nRows; ++i)
      {
         /* if row i is constant */
         if ( SCIPisEQ(scip, solu1[i], solu2[i]) )
            continue;

         /* if first non-constant row is [1,0]: feasible */
         if ( SCIPisEQ(scip, solu1[i], 1.0) )
         {
            break;
         }
         else /* infeasible */
         {
            SCIPdebugMessage("Solution is infeasible.\n");
            *result = SCIP_INFEASIBLE;
            terminated = TRUE;

            if ( printreason )
               SCIPinfoMessage(scip, NULL, "First non-constant row (%u) of orbisack is (0,1).\n", i);

            break;
         }
      }

      /* free buffers */
      SCIPfreeBufferArrayNull(scip, &solu2);
      SCIPfreeBufferArrayNull(scip, &solu1);

      if ( terminated )
         break;
   }

   if ( ! terminated )
      SCIPdebugMessage("Solution is feasible.\n");

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropOrbisack)
{
   int c;

   /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMessage("Propagation method of orbisack constraint handler.\n");

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      SCIP_Bool found = FALSE;
      unsigned int nGen = 0;

      assert( conss[c] != NULL );

      SCIP_CALL( propVariables(scip, conss[c], &infeasible, &found, &nGen) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      if ( found )
         *result = SCIP_REDUCEDDOM;
   }

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolOrbisack)
{  /*lint --e{715}*/
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMessage("Presolving method of orbisack constraint handler. Propagating orbisack inequalities.\n");
   *result = SCIP_DIDNOTFIND;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      SCIP_Bool found = FALSE;
      unsigned int nGen = 0;

      assert( conss[c] != NULL );
      SCIP_CALL( propVariables(scip, conss[c], &infeasible, &found, &nGen) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         break;
      }

      if ( nGen > 0 )
      {
         *nfixedvars += (int) nGen;
         *result = SCIP_SUCCESS;
      }
   }

   return SCIP_OKAY;
}


/** Propagation resolution for conflict analysis */
static
SCIP_DECL_CONSRESPROP(consRespropOrbisack)
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
#ifndef NDEBUG
   unsigned int nRows;
#endif
   unsigned int i;

   /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( infervar != NULL );
   assert( bdchgidx != NULL );
   assert( result != NULL );

   SCIPdebugMessage("Propagation resolution method of orbisack constraint handler.\n");

   *result = SCIP_DIDNOTFIND;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL);
   assert( consdata->nRows_ > 0 );
   assert( consdata->vars1_ != NULL );
   assert( consdata->vars2_ != NULL );

   vars1 = consdata->vars1_;
   vars2 = consdata->vars2_;

#ifndef NDEBUG
   nRows = consdata->nRows_;
   assert( 0 <= inferinfo && inferinfo < (int) nRows );
#endif

   assert( vars2[inferinfo] == infervar );
   assert( SCIPvarGetUbAtIndex(vars2[inferinfo], bdchgidx, FALSE) > 0.5 && SCIPvarGetUbAtIndex(vars2[inferinfo], bdchgidx, TRUE) < 0.5 );

   if ( SCIPvarGetUbAtIndex(vars2[inferinfo], bdchgidx, FALSE) > 0.5 && SCIPvarGetUbAtIndex(vars2[inferinfo], bdchgidx, TRUE) < 0.5 )
   {
      SCIPdebugMessage(" -> reason for setting x[%i][1] = 0 was fixing x[%i][0] to 0 and each row above is fixed to the same value.\n",
         inferinfo, inferinfo);

      for (i = 0; (int) i < inferinfo; ++i)
      {
         SCIP_CALL( SCIPaddConflictUb(scip, vars1[i], bdchgidx) );
         SCIP_CALL( SCIPaddConflictLb(scip, vars1[i], bdchgidx) );
         SCIP_CALL( SCIPaddConflictUb(scip, vars2[i], bdchgidx) );
         SCIP_CALL( SCIPaddConflictLb(scip, vars2[i], bdchgidx) );
      }
      SCIP_CALL( SCIPaddConflictUb(scip, vars1[inferinfo], bdchgidx) );

      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}


/** Lock variables
 *
 *  We assume we have only one global (void) constraint and lock all variables.
 *
 * - Orbisack constraints may get violated if the variables of the first column
 *   are rounded down, we therefor call SCIPaddVarLocks(..., nlockspos, nlocksneg).
 * - Orbisack constraints may get violated if the variables of the second column
 *   are rounded up , we therefor call SCIPaddVarLocks(..., nlocksneg, nlockspo ).
 */
static
SCIP_DECL_CONSLOCK(consLockOrbisack)
{
   SCIP_CONSDATA* consdata;
   unsigned int nRows;
   unsigned int i;
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;

   /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   SCIPdebugMessage("Locking method for orbisack constraint handler.\n");

   /* get data of original constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL);
   assert( consdata->nRows_ > 0 );
   assert( consdata->vars1_ != NULL );
   assert( consdata->vars2_ != NULL );

   nRows = consdata->nRows_;
   vars1 = consdata->vars1_;
   vars2 = consdata->vars2_;

   for (i = 0; i < nRows; ++i)
   {
      SCIP_CALL( SCIPaddVarLocks(scip, vars1[i], nlockspos, nlocksneg) );
      SCIP_CALL( SCIPaddVarLocks(scip, vars2[i], nlocksneg, nlockspos) );
   }

   return SCIP_OKAY;
}


/** constraint display method of constraint handler
 *
 *  The constraint handler should output a representation of the constraint into the given text file.
 */
static
SCIP_DECL_CONSPRINT(consPrintOrbisack)
{
   /*lint --e{715}*/
   return SCIP_OKAY;
}


#if 0
/** constraint copying method of constraint handler
 *
 *  The constraint handler can provide a copy method which copy a constraint from one SCIP data structure into another
 *  SCIP data structure.
 */
static
SCIP_DECL_CONSCOPY(consCopyOrbisack)
{
   SCIP_CONSDATA* sourcedata;
   unsigned int nRows;
   SCIP_VAR** sourceVars1;
   SCIP_VAR** sourceVars2;
   SCIP_VAR** vars1 = NULL;
   SCIP_VAR** vars2 = NULL;
   unsigned int i;

   /*lint --e{715}*/
   assert( scip != NULL );
   assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(sourcecons)), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( sourcescip != NULL );
   assert( sourcecons != NULL );
   assert( varmap != NULL );
   assert( valid != NULL );

   *valid = TRUE;

   SCIPdebugMessage("Copying method for orbisack constraint handler.\n");

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL);
   assert( sourcedata->vars1_ != NULL );
   assert( sourcedata->vars2_ != NULL );
   assert( sourcedata->vals1_ != NULL );
   assert( sourcedata->vals2_ != NULL );
   assert( sourcedata->nRows_  > 0 );

   /* copy constraint data */
   nRows = sourcedata->nRows_;
   sourceVars1 = sourcedata->vars1_;
   sourceVars2 = sourcedata->vars2_;

   /* separately allocate space to account for unsuccessful copying */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars1, (int) nRows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars2, (int) nRows) );

   /* get copies */
   for (i = 0; i < nRows; ++i)
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourceVars1[i], &(vars1[i]), varmap, consmap, global, valid) );
      assert( *valid );
      assert( vars1[i] != NULL );

      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourceVars2[i], &(vars2[i]), varmap, consmap, global, valid) );
      assert( *valid );
      assert( vars2[i] != NULL );
   }

   /* create copied constraint */
   if ( name == NULL )
      name = SCIPconsGetName(sourcecons);

   SCIP_CALL( SCIPcreateConsOrbisack(scip, cons, name, vars1, vars2, nRows, sourcedata->isPporbisack_,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic,
         removable, stickingatnode) );

   SCIPfreeBufferArrayNull(scip, &vars2);
   SCIPfreeBufferArrayNull(scip, &vars1);

   return SCIP_OKAY;
}
#endif


/** creates the handler for orbisack constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrOrbisack(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;
   SCIP_CONSHDLR* conshdlr;

   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpOrbisack, consEnfopsOrbisack, consCheckOrbisack, consLockOrbisack,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
#if 0
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyOrbisack, consCopyOrbisack) );
#endif
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeOrbisack) );

   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteOrbisack) );
#if ( SCIP_VERSION >= 320 )
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolOrbisack, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
#else
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolOrbisack, CONSHDLR_MAXPREROUNDS, CONSHDLR_DELAYPRESOL) );
#endif
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintOrbisack) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropOrbisack, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropOrbisack) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpOrbisack, consSepasolOrbisack, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransOrbisack) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpOrbisack) );

   /* separation methods */
   SCIP_CALL( SCIPaddBoolParam(scip, "cons/orbisack/coverSeparation", "Separate cover inequalities for orbisacks?", &conshdlrdata->coverSeparation_, TRUE,
         DEFAULT_COVERSEPARATION, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "cons/orbisack/orbiSeparation", "Separate orbisack inequalities?", &conshdlrdata->orbiSeparation_, TRUE,
         DEFAULT_ORBISEPARATION, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "cons/orbisack/enforcing", "Enforce orbisack constraints?", &conshdlrdata->orbiEnforcing_, TRUE,
         DEFAULT_ENFORCING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "cons/orbisack/check", "Check orbisack constraints?", &conshdlrdata->orbiCheck_, TRUE,
         DEFAULT_CHECK, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "cons/orbisack/coeffBound", "Maximum size of coefficients for orbisack inequalities", &conshdlrdata->coeffBound_, TRUE,
         DEFAULT_COEFFBOUND, 0.0, DBL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates and captures a orbisack constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsOrbisack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*const*       vars1,              /**< first column of matrix of variables on which the symmetry acts */
   SCIP_VAR*const*       vars2,              /**< second column of matrix of variables on which the symmetry acts */
   unsigned int          nrows,              /**< number of rows in variable matrix */
   SCIP_Bool             isPporbisack,       /**< whether the orbisack is a packing/partitioning orbisack */
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

   /* find the orbisack constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("orbisack constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   assert( nrows > 0 );

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, vars1, vars2, nrows, isPporbisack) );

   /* do we use enforcing? */
   SCIP_CALL( SCIPgetBoolParam(scip, "cons/orbisack/enforcing", &enforce) );

   /* do we use check routines? */
   SCIP_CALL( SCIPgetBoolParam(scip, "cons/orbisack/check", &check) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}


/** creates and captures an orbisack constraint in its most basic variant
 *
 *  All constraint flags set to their default values, which can be set afterwards using SCIPsetConsFLAGNAME() in scip.h.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicOrbisack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR**            vars1,              /**< first column of matrix of variables on which the symmetry acts */
   SCIP_VAR**            vars2,              /**< second column of matrix of variables on which the symmetry acts */
   unsigned int          nRows               /**< number of rows in constraint matrix */
   )
{
   SCIP_CALL( SCIPcreateConsOrbisack(scip, cons, name, vars1, vars2, nRows, DEFAULT_PPORBISACK,
         TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

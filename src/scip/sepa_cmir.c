/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: sepa_cmir.c,v 1.21 2004/10/05 11:01:38 bzfpfend Exp $"

/**@file   sepa_cmir.c
 * @brief  complemented mixed integer rounding cuts separator (Marchand's version)
 * @author Kati Wolter
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "sepa_cmir.h"


#define SEPA_NAME              "cmir"
#define SEPA_DESC              "complemented mixed integer rounding cuts separator (Marchand's version)"
#define SEPA_PRIORITY             -1000
#define SEPA_FREQ                    20

#define DEFAULT_MAXROUNDS             3 /**< maximal number of cmir separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        -1 /**< maximal number of cmir separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXTRIES            100 /**< maximal number of rows to start aggregation with per separation round
                                         *   (-1: unlimited) */
#define DEFAULT_MAXTRIESROOT         -1 /**< maximal number of rows to start aggregation with per round in the root node
                                         *   (-1: unlimited) */
#define DEFAULT_MAXAGGRS              4 /**< maximal number of aggregations for each row per separation round */
#define DEFAULT_MAXAGGRSROOT          8 /**< maximal number of aggreagtions for each row per round in the root node */
#define DEFAULT_MAXSEPACUTS          50 /**< maximal number of cmir cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT     500 /**< maximal number of cmir cuts separated per separation round in root node */
#define DEFAULT_MAXSLACK            0.0 /**< maximal slack of rows to be used in aggregation */
#define DEFAULT_MAXSLACKROOT        0.1 /**< maximal slack of rows to be used in aggregation in the root node */
#define DEFAULT_SLACKSCORE         1e-3 /**< weight of slack in the aggregation scoring of the rows */
#define DEFAULT_MAXROWFAC          1e+4 /**< maximal row aggregation factor */
#define DEFAULT_MINROWFAC         -1e+4 /**< minimal row aggregation factor */
#define DEFAULT_MAXTESTDELTA         20	/**< maximal number of different deltas to try */
#define DEFAULT_MAXTESTDELTAROOT    100 /**< maximal number of different deltas to try in the root node */
#define DEFAULT_MAXCONTS             10 /**< maximal number of active continuous variables in aggregated row */
#define DEFAULT_MAXCONTSROOT         10 /**< maximal number of active continuous variables in aggregated row in the root */
#define DEFAULT_DYNAMICCUTS        TRUE /**< should generated cuts be removed from the LP if they are no longer tight? */

#define BOUNDSWITCH                 0.5
#define USEVBDS                    TRUE
#define ALLOWLOCAL                 TRUE


/*
 * Data structures
 */

/** separator data */
struct SepaData
{
   Real             maxslack;           /**< maximal slack of rows to be used in aggregation */
   Real             maxslackroot;       /**< maximal slack of rows to be used in aggregation in the root node */
   Real             slackscore;         /**< weight of slack in the aggregation scoring of the rows */
   int              maxrounds;          /**< maximal number of cmir separation rounds per node (-1: unlimited) */
   int              maxroundsroot;      /**< maximal number of cmir separation rounds in the root node (-1: unlimited) */
   int              maxtries;           /**< maximal number of rows to start aggregation with per separation round
                                         *   (-1: unlimited) */
   int              maxtriesroot;       /**< maximal number of rows to start aggregation with per round in the root node
                                         *   (-1: unlimited) */
   int              maxaggrs;           /**< maximal number of aggregations for each row per separation round */
   int              maxaggrsroot;       /**< maximal number of aggreagtions for each row per round in the root node */
   int              maxsepacuts;        /**< maximal number of cmir cuts separated per separation round */
   int              maxsepacutsroot;    /**< maximal number of cmir cuts separated per separation round in root node */
   int              maxrowfac;          /**< maximal row aggregation factor */
   int              minrowfac;          /**< minimal row aggregation factor */
   int              maxtestdelta;	/**< maximal number of different deltas to try */
   int              maxtestdeltaroot;   /**< maximal number of different deltas to try in the root node */
   int              maxconts;	        /**< maximal number of active continuous variables in aggregated row */
   int              maxcontsroot;       /**< maximal number of active continuous variables in aggregated row in the root */
   Bool             dynamiccuts;        /**< should generated cuts be removed from the LP if they are no longer tight? */
};




/*
 * Local methods
 */

/** adds given cut to LP if violated */
static
RETCODE addCut(
   SCIP*            scip,               /**< SCIP data structure */
   SEPADATA*        sepadata,           /**< separator data */
   Real*            varsolvals,         /**< solution values of active variables */
   Real*            cutcoefs,           /**< coefficients of active variables in cut */
   Real             cutrhs,             /**< right hand side of cut */
   Bool             cutislocal,         /**< is the cut only locally valid? */
   int*             ncuts               /**< pointer to count the number of added cuts */
   )
{
   VAR** vars;
   COL** cutcols;
   Real* cutvals;
   Real cutact;
   Real cutsqrnorm;
   Real cutnorm;
   Real val;
   int nvars;
   int cutlen;
   int v;
   Bool success;
   
   assert(scip != NULL);
   assert(sepadata != NULL);      
   assert(varsolvals != NULL);
   assert(cutcoefs != NULL);
   assert(ncuts != NULL);

   /* get active problem variables */
   CHECK_OKAY( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   assert(nvars == 0 || vars != NULL);

   /* get temporary memory for storing the cut as sparse row */
   CHECK_OKAY( SCIPallocBufferArray(scip, &cutcols, nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &cutvals, nvars) );
   
   /* store the cut as sparse row, calculate activity of cut */
   cutlen = 0;
   cutact = 0.0;
   cutsqrnorm = 0.0;
   for( v = 0; v < nvars; ++v )
   {
      val = cutcoefs[v];
      if( !SCIPisZero(scip, val) )
      {
         assert(SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_COLUMN);
         cutact += val * varsolvals[v];
         cutsqrnorm += SQR(val);
         cutcols[cutlen] = SCIPvarGetCol(vars[v]);
         cutvals[cutlen] = val;
         cutlen++;
      }
   }
   cutnorm = SQRT(cutsqrnorm);
   
   if( SCIPisPositive(scip, cutnorm) && SCIPisEfficacious(scip, (cutact - cutrhs)/cutnorm) )
   {
      ROW* cut;
      char cutname[MAXSTRLEN];
      
      /* create the cut */
      sprintf(cutname, "cmir%d_%d", SCIPgetNLPs(scip), *ncuts);
      CHECK_OKAY( SCIPcreateRow(scip, &cut, cutname, cutlen, cutcols, cutvals, -SCIPinfinity(scip), cutrhs, 
            cutislocal, FALSE, sepadata->dynamiccuts) );

      debugMessage(" -> found potential c-mir cut <%s>: activity=%f, rhs=%f, norm=%f, eff=%f\n",
         cutname, cutact, cutrhs, cutnorm, SCIPgetCutEfficacy(scip, cut));
      debug(SCIPprintRow(scip, cut, NULL));
      
#if 0
      /* try to scale the cut to integral values */
      CHECK_OKAY( SCIPmakeRowIntegral(scip, cut, 1000, 10000.0, &success) );
      if( success && !SCIPisCutEfficacious(scip, cut) )
      {
         debugMessage(" -> c-mir cut <%s> no longer efficacious: act=%f, rhs=%f, norm=%f, eff=%f\n",
            cutname, cutact, cutrhs, cutnorm, SCIPgetCutEfficacy(scip, cut));
         debug(SCIPprintRow(scip, cut, NULL));
         success = FALSE;
      }
#else
      success = TRUE;
#endif

      /* if scaling was successful, add the cut */
      if( success )
      {
         debugMessage(" -> found c-mir cut <%s>: act=%f, rhs=%f, norm=%f, eff=%f\n",
            cutname, cutact, cutrhs, cutnorm, SCIPgetCutEfficacy(scip, cut));
         debug(SCIPprintRow(scip, cut, NULL));
         CHECK_OKAY( SCIPaddCut(scip, cut, FALSE) );
         if( !cutislocal )
         {
            CHECK_OKAY( SCIPaddPoolCut(scip, cut) );
         }
         (*ncuts)++;
      }
      
      /* release the row */
      CHECK_OKAY( SCIPreleaseRow(scip, &cut) );
   }
   
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &cutvals);
   SCIPfreeBufferArray(scip, &cutcols);

   return SCIP_OKAY;   
}

/** adds delta to continuous variables counters */
static
void updateNConts(
   SCIP*            scip,               /**< SCIP data structure */ 
   COL*             col,                /**< column of continuous variable */
   int              delta,              /**< delta value of counters */
   int*             nconts,             /**< pointer to count number of continuous variabls */
   int*             nactiveconts        /**< pointer to count number of active continuous variabls */
   )
{
   Real primsol;
   Real lb;
   Real ub;

   assert(nconts != NULL);
   assert(nactiveconts != NULL);

   (*nconts) += delta;

   primsol = SCIPcolGetPrimsol(col);
   lb = SCIPcolGetLb(col);
   ub = SCIPcolGetUb(col);
   assert(SCIPisEQ(scip, primsol, SCIPvarGetLPSol(SCIPcolGetVar(col))));
   assert(SCIPisEQ(scip, lb, SCIPvarGetLbLocal(SCIPcolGetVar(col))));
   assert(SCIPisEQ(scip, ub, SCIPvarGetUbLocal(SCIPcolGetVar(col))));

   if( SCIPisLT(scip, lb, primsol) && SCIPisLT(scip, primsol, ub) )
      (*nactiveconts) += delta;
}

/** decreases the score of a row in order to not aggregate it again too soon */
static
void decreaseRowScore(
   SCIP*            scip,               /**< SCIP data structure */ 
   Real*            rowlhsscores,       /**< aggregation scores for left hand sides of row */
   Real*            rowrhsscores,       /**< aggregation scores for right hand sides of row */
   int              rowidx              /**< index of row to decrease score for */
   )
{
   assert(rowlhsscores != NULL);
   assert(rowrhsscores != NULL);
   assert(rowlhsscores[rowidx] < 0.0);
   assert(rowrhsscores[rowidx] < 0.0);

   if( !SCIPisInfinity(scip, -rowlhsscores[rowidx]) )
      rowlhsscores[rowidx] *= 1.1;
   if( !SCIPisInfinity(scip, -rowrhsscores[rowidx]) )
      rowrhsscores[rowidx] *= 1.1;
}

/** aggregates different single mixed integer constraints by taking linear combinations of the rows of the LP  */
static
RETCODE aggregation(
   SCIP*            scip,               /**< SCIP data structure */ 
   SEPADATA*        sepadata,           /**< separator data */
   Real*            varsolvals,         /**< LP solution value of all variables in LP */
   Real*            rowlhsscores,       /**< aggregation scores for left hand sides of row */
   Real*            rowrhsscores,       /**< aggregation scores for right hand sides of row */
   int              startrow,           /**< index of row to start aggregation */ 
   int              maxaggrs,           /**< maximal number of aggregations */
   Real             maxslack,           /**< maximal slack of rows to be used in aggregation */
   int              maxtestdelta,	/**< maximal number of different deltas to try */
   int              maxconts,	        /**< maximal number of active continuous variables in aggregated row */
   int*             ncuts               /**< pointer to count the number of generated cuts */
   )
{
   Real* aggrcoefs;       /* coefficients of all variables in aggregated row */
   Real* rowweights;      /* weight of rows in all aggregations */ 
   Real* testeddeltas;

   int nstartnonzcols;    /* number of nonzero columns of startrow */
   COL** startnonzcols;   /* columns with nonzero coefficients of startrow */
   Real* startnonzcoefs;  /* nonzero coefficients of startrow */    
   Real startrowact;      /* activity of startrow */

   Real* cutcoefs;         /* coefficients of variables in cut */
   Real cutrhs;            /* right hand side of the cut */
   Bool success;
   Bool cutislocal;

   int naggrs;
   int nconts;
   int nactiveconts;
   COL* bestcol;          

   VAR** vars;
   int nvars;
   COL** cols;
   int ncols;
   ROW** rows;
   int nrows;
   int c;
   int r;

   assert(scip != NULL);
   assert(sepadata != NULL);      
   assert(varsolvals != NULL);
   assert(rowlhsscores != NULL);
   assert(rowrhsscores != NULL);
   assert(ncuts != NULL);

   /* get active problem variables and LP data */
   CHECK_OKAY( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   assert(nvars == 0 || vars != NULL);
   CHECK_OKAY( SCIPgetLPColsData(scip, &cols, &ncols) );
   assert(ncols == 0 || cols != NULL);
   CHECK_OKAY( SCIPgetLPRowsData(scip, &rows, &nrows) );
   assert(nrows == 0 || rows != NULL);
   assert(0 <= startrow && startrow < nrows);

   debugMessage("start c-MIR aggregation with row <%s> (%d/%d)\n", SCIProwGetName(rows[startrow]), startrow, nrows);

   /* get temporary memory */
   CHECK_OKAY( SCIPallocBufferArray(scip, &aggrcoefs, ncols) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &rowweights, nrows) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &cutcoefs, nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &testeddeltas, ncols) );

   /* initialize weights of rows in aggregation */
   for( r = 0; r < nrows; r++ )
      rowweights[r] = 0.0;
   startrowact = SCIPgetRowActivity(scip, rows[startrow]);
   if( startrowact <= 0.5 * SCIProwGetLhs(rows[startrow]) + 0.5 * SCIProwGetRhs(rows[startrow]) )
      rowweights[startrow] = -1.0;
   else 
      rowweights[startrow] = 1.0;

   /* decrease score of startrow in order to not aggregate it again too soon */
   decreaseRowScore(scip, rowlhsscores, rowrhsscores, startrow);

   /* get nonzero columns and coefficients of startrow */
   startnonzcols =  SCIProwGetCols(rows[startrow]);
   nstartnonzcols = SCIProwGetNLPNonz(rows[startrow]);
   startnonzcoefs = SCIProwGetVals(rows[startrow]);

   /* for all columns of startrow store coefficient as coefficient in aggregated row */ 
   clearMemoryArray(aggrcoefs, ncols);
   nconts = 0;
   nactiveconts = 0;
   for( c = 0; c < nstartnonzcols; c++ )
   {
      VAR* var;
      int pos;

      var = SCIPcolGetVar(startnonzcols[c]);
      pos = SCIPcolGetLPPos(startnonzcols[c]);
      assert(pos >= 0); 
      aggrcoefs[pos] = rowweights[startrow] * startnonzcoefs[c];
      if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
         updateNConts(scip, startnonzcols[c], +1, &nconts, &nactiveconts);
   }

   /* try to generate cut from the current aggregated row 
    * add cut if found, otherwise add another row to aggregated row 
    * in order to get rid of a continuous variable
    */
   naggrs = 0;
   while( nactiveconts <= maxconts && naggrs <= maxaggrs )
   {
      Real bestdelta;
      Real bestviolation; 
      int ntesteddeltas;

      ROW* bestrow;          
      COL** bestrownonzcols;     /* columns with nonzero coefficients in best row to add */
      Real* bestrownonzcoefs;    /* nonzero coefficients of columns in best row to add */
      int nbestrownonzcols;      /* number of columns with nonzero coefficients in best row to add */
      Real bestrowact;           /* activity of best row to add */
      Real bestbounddist;
      Real bestscore;
      Real aggrfact;         

#ifdef DEBUG
      debugMessage("aggregation of startrow %d and %d additional rows with %d continuous variables (%d active):\n",
         startrow, naggrs, nconts, nactiveconts);
      for( c = 0; c < ncols; ++c )
      {
         if( aggrcoefs[c] != 0.0 )
            printf(" %+g<%s>(%g)", aggrcoefs[c], SCIPvarGetName(SCIPcolGetVar(cols[c])),
               SCIPvarGetLPSol(SCIPcolGetVar(cols[c])));
      }
      printf("\n");
#endif

      /* Step 1: try to generate a MIR cut out of the current aggregation */

      /* search delta for generating a cut with maximum violation: 
       * delta = coefficient of integer variable, which lies between its bounds
       */ 
      ntesteddeltas = 0;
      bestdelta = 0.0;
      bestviolation = 0.0;
      for( c = 0; c < ncols && ntesteddeltas < maxtestdelta; c++ )
      {
         VAR* var;
         Real primsol;
         Real lb;
         Real ub;
         Real delta;
         Real cutact;
         Real violation;
         Bool tested;
         int i;

         /* ignore variables not existing in current aggregation (or with close to zero coefficients) */
         if( SCIPisZero(scip, aggrcoefs[c]) )
            continue;

         /* ignore continuous variables */
         var = SCIPcolGetVar(cols[c]);
         if( !SCIPvarIsIntegral(var) )
            continue;

         assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS);
         primsol = SCIPcolGetPrimsol(cols[c]);
         lb = SCIPcolGetLb(cols[c]);
         ub = SCIPcolGetUb(cols[c]);
         assert(SCIPisEQ(scip, primsol, varsolvals[SCIPvarGetProbindex(var)]));
         assert(SCIPisEQ(scip, lb, SCIPvarGetLbLocal(var)));
         assert(SCIPisEQ(scip, ub, SCIPvarGetUbLocal(var)));

         /* ignore variables with current solution value on its bounds */
         if( SCIPisEQ(scip, primsol, lb) || SCIPisEQ(scip, primsol, ub) )
            continue;

         /* try to divide the aggregation by this coefficient */
         delta = 1 / ABS(aggrcoefs[c]);

         /* check, if delta was already tested */
         tested = FALSE;
         for( i = 0; i < ntesteddeltas && !tested; i++ )
            tested = SCIPisEQ(scip, testeddeltas[i], delta);
         if( tested )
            continue;
         
         testeddeltas[ntesteddeltas] = delta;
         ntesteddeltas++;
         
         /* create a MIR cut out of the weighted LP rows */
         CHECK_OKAY( SCIPcalcMIR(scip, BOUNDSWITCH, USEVBDS, ALLOWLOCAL, 0.05, rowweights, delta, 
               cutcoefs, &cutrhs, &cutact, &success, &cutislocal) );
         assert(ALLOWLOCAL || !cutislocal);
         debugMessage("delta = %g -> success: %d\n", delta, success);

         /* delta generates cut which is more violated */
         if( success )
         {
            violation = cutact - cutrhs;
            debugMessage("act = %g  rhs = %g  viol = %g, old bestviol = %g\n", 
               cutact, cutrhs, violation, bestviolation);
            if( violation > bestviolation )
            {
               bestdelta = delta;
               bestviolation = violation;
            }
         }
      }
  
      /* delta found */
      if( SCIPisFeasPositive(scip, bestviolation) )
      {
         Real cutact;
         Real violation;
         Real delta;
         Bool tested;
         int i;
         int j;
         int oldncuts;

         assert(bestdelta != 0.0);

         /* Try to improve violation by multiplying delta with 2, 4 and 8 */
         for( i = 0, delta = bestdelta; i < 3; i++, delta *= 2.0 )
         {
            /* check, if delta was already tested */
            tested = FALSE;
            for( j = 0; j < ntesteddeltas && !tested; j++ )
               tested = SCIPisEQ(scip, testeddeltas[j], delta);
            if( tested )
               continue;

            /* create a MIR cut out of the weighted LP rows */
            CHECK_OKAY( SCIPcalcMIR(scip, BOUNDSWITCH, USEVBDS, ALLOWLOCAL, 0.05, rowweights, delta, 
                  cutcoefs, &cutrhs, &cutact, &success, &cutislocal) );
            assert(ALLOWLOCAL || !cutislocal);
            debugMessage("delta = %g -> success: %d\n", delta, success);
            if( success )
            {
               violation = cutact - cutrhs;
               if( violation > bestviolation )
               {
                  bestdelta = delta;
                  bestviolation = violation;
               }
            }
         }
         
         /* generate cut with bestdelta */
         oldncuts = *ncuts;
         CHECK_OKAY( SCIPcalcMIR(scip, BOUNDSWITCH, USEVBDS, ALLOWLOCAL, 0.05, rowweights, bestdelta, 
               cutcoefs, &cutrhs, &cutact, &success, &cutislocal) );
         assert(ALLOWLOCAL || !cutislocal);
         CHECK_OKAY( addCut(scip, sepadata, varsolvals, cutcoefs, cutrhs, cutislocal, ncuts) );

         /* if the cut was successfully added, abort the aggregation of further rows */
         if( *ncuts > oldncuts )
            break;
      }

      /* abort, if no more active continuous variable is left or if we reached the maximal number of aggregations */
      if( nactiveconts == 0 || naggrs == maxaggrs )
         break;


      /* Step 2: aggregate an additional row in order to remove a continuous variable */
      debugMessage(" -> search column to eliminate\n");

      /* search for "best" continuous variable in aggregated row:
       * - solution value is strictly between lower and upper bound
       * - it exists a not yet aggregated row with nonzero coefficient in this column
       * out of these variables:
       * - prefer variables with larger distance of current solution value to its bounds
       * - of those with large bound distance, prefer variables that can be eliminated with a row of high score
       */
      bestcol = NULL;
      bestbounddist = 0.0;
      bestscore = -SCIPinfinity(scip);
      bestrow = NULL;
      aggrfact = 0.0;
      for( c = 0; c < ncols; c++ )
      {
         VAR* var;
         COL* col;
         Real primsol;
         Real lb;
         Real ub;
         Real distlower;
         Real distupper;
         Real bounddist;
         
         /* ignore columns not in current aggregation */
         if( aggrcoefs[c] == 0.0 )
            continue;

         col = cols[c];
         var = SCIPcolGetVar(col);
         
         /* ignore integral variables */
         if( SCIPvarIsIntegral(var) )
            continue;

         /* get minimum distance of LP solution value of variable to its bounds */
         primsol = SCIPcolGetPrimsol(col);
         lb = SCIPcolGetLb(col);
         ub = SCIPcolGetUb(col);
         distlower = primsol - lb;
         distupper = ub - primsol;
         bounddist = MIN(distlower, distupper);

         /* check, if variable is candidate to be the new best variable */
         if( SCIPisPositive(scip, bounddist) && bounddist >= bestbounddist - 0.1 )
         {
            ROW** nonzrows;
            Real* nonzcoefs;
            int nnonzrows;

            debugMessage("     -> col <%s>[%g,%g]: sol=%g, dist=%g\n", SCIPvarGetName(var), lb, ub, primsol, bounddist);

            /* look for "best" row to add (minimal slack), but don't add rows again,
             * that are already involved in aggregation
             */
            nnonzrows = SCIPcolGetNLPNonz(col);
            nonzrows = SCIPcolGetRows(col);
            nonzcoefs = SCIPcolGetVals(col);
                  
            for( r = 0; r < nnonzrows; r++ )
            {
               Real score;
               Real fact;
               Real activity;
               int lppos;

               lppos = SCIProwGetLPPos(nonzrows[r]);
               assert(0 <= lppos && lppos < nrows);

               debugMessage("        -> row <%s>: weight=%g, pos=%d, fact=%g, %g <= %g <= %g\n",
                  SCIProwGetName(nonzrows[r]), rowweights[lppos], lppos, - aggrcoefs[c] / nonzcoefs[r],
                  SCIProwGetLhs(nonzrows[r]), SCIPgetRowLPActivity(scip, nonzrows[r]), SCIProwGetRhs(nonzrows[r]));

               /* take only unmodifiable LP rows, that are not yet aggregated */
               if( rowweights[lppos] != 0.0 || SCIProwIsModifiable(nonzrows[r]) )
                  continue;

               /* don't aggregate rows that would lead to a too extreme aggregation factor */
               fact = - aggrcoefs[c] / nonzcoefs[r]; 
               if( fact < sepadata->minrowfac || fact > sepadata->maxrowfac || SCIPisZero(scip, fact) )
                  continue;

               /* check, if the row's slack multiplied with the aggregation factor is too large */
               activity = SCIPgetRowLPActivity(scip, nonzrows[r]);
               if( (fact < 0.0 && SCIPisGT(scip, fact * (SCIProwGetLhs(nonzrows[r]) - activity), maxslack))
                  || (fact > 0.0 && SCIPisGT(scip, fact * (SCIProwGetRhs(nonzrows[r]) - activity), maxslack)) )
                  continue;

               /* choose row with best aggregation score */
               assert(!SCIPisInfinity(scip, -SCIProwGetLhs(nonzrows[r])) || SCIPisInfinity(scip, -rowlhsscores[lppos]));
               assert(!SCIPisInfinity(scip, SCIProwGetRhs(nonzrows[r])) || SCIPisInfinity(scip, -rowrhsscores[lppos]));
               score = (fact < 0.0 ? rowlhsscores[lppos] : rowrhsscores[lppos]);
               if( !SCIPisInfinity(scip, -score)
                  && (bounddist > bestbounddist + 0.1 || score > bestscore) )
               {
                  bestbounddist = bounddist;
                  bestscore = score; 
                  bestcol = col;
                  bestrow = nonzrows[r];
                  aggrfact = fact;
                  debugMessage("     -> col <%s>: %g * row <%s>, bounddist=%g, score=%g\n",
                     SCIPvarGetName(SCIPcolGetVar(bestcol)), aggrfact, SCIProwGetName(bestrow), bestbounddist, score);
               }
            }
         }
      }
      assert((bestcol == NULL) == (bestrow == NULL));
         
      /* abort, if no row can be added to remove an additional active continuous variable */
      if( bestcol == NULL )
         break;

            
      /* Step 3: add row to aggregation */
      debugMessage(" -> adding %+g<%s> to eliminate variable <%s> (aggregation %d)\n", 
         aggrfact, SCIProwGetName(bestrow), SCIPvarGetName(SCIPcolGetVar(bestcol)), naggrs+1);
      assert(rowweights[SCIProwGetLPPos(bestrow)] == 0.0);
      assert(!SCIPisZero(scip, aggrfact));

      /* change row's aggregation weight */
      rowweights[SCIProwGetLPPos(bestrow)] = aggrfact;
               
      /* decrease score of aggregation row in order to not aggregate it again too soon */
      decreaseRowScore(scip, rowlhsscores, rowrhsscores, SCIProwGetLPPos(bestrow));

      /* change coefficients of aggregation and update the number of continuous variables */
      bestrownonzcols = SCIProwGetCols(bestrow);
      bestrownonzcoefs = SCIProwGetVals(bestrow);
      nbestrownonzcols = SCIProwGetNLPNonz(bestrow);
      for( c = 0; c < nbestrownonzcols; c++ )
      {
         VAR* var;
         int pos;

         var = SCIPcolGetVar(bestrownonzcols[c]);
         pos = SCIPcolGetLPPos(bestrownonzcols[c]);
         assert(pos >= 0);

         if( aggrcoefs[pos] != 0.0 && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
            updateNConts(scip, bestrownonzcols[c], -1, &nconts, &nactiveconts);
         aggrcoefs[pos] += bestrownonzcoefs[c] * aggrfact;
         if( SCIPisZero(scip, aggrcoefs[pos]) )
            aggrcoefs[pos] = 0.0;
         else if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
            updateNConts(scip, bestrownonzcols[c], +1, &nconts, &nactiveconts);
      }
      naggrs++;

      debugMessage(" -> %d continuous variables left (%d/%d active), %d/%d aggregations\n", 
         nconts, nactiveconts, maxconts, naggrs, maxaggrs);
   }

   /* free datastructures */
   SCIPfreeBufferArray(scip, &testeddeltas);
   SCIPfreeBufferArray(scip, &cutcoefs);
   SCIPfreeBufferArray(scip, &rowweights);
   SCIPfreeBufferArray(scip, &aggrcoefs);

   return SCIP_OKAY; 
}



/*
 * Callback methods of separator
 */

/** destructor of separator to free user data (called when SCIP is exiting) */
static
DECL_SEPAFREE(sepaFreeCmir)
{  /*lint --e{715}*/
   SEPADATA* sepadata;

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;

}


/** initialization method of separator (called when problem solving starts) */
#define sepaInitCmir NULL


/** deinitialization method of separator (called when problem solving exits) */
#define sepaExitCmir NULL


/** execution method of separator */
static
DECL_SEPAEXEC(sepaExecCmir)
{  /*lint --e{715}*/
   SEPADATA* sepadata;
   VAR** vars;
   Real* varsolvals;
   ROW** rows;     
   Real* rowlhsscores;
   Real* rowrhsscores;
   Real* rowscores;
   int* roworder;
   Real maxslack;
   int nvars;
   int nrows;
   int ntries;
   int depth;
   int ncalls;
   int maxtries;
   int maxaggrs;
   int maxsepacuts;
   int maxtestdelta;
   int maxconts;
   int ncuts;
   int r;

   assert(sepa != NULL);
   assert(scip != NULL);
 
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   depth = SCIPgetDepth(scip);
   ncalls = SCIPsepaGetNCallsAtNode(sepa);

   /* only call the cmir cut separator a given number of times at each node */
   if( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;

   /* get all rows and number of columns */
   CHECK_OKAY( SCIPgetLPRowsData(scip, &rows, &nrows) ); 
   assert(nrows == 0 || rows != NULL);

   /* nothing to do, if LP is empty */
   if( nrows == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* get active problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   assert(nvars == 0 || vars != NULL);

   /* get data structure */
   CHECK_OKAY( SCIPallocBufferArray(scip, &rowlhsscores, nrows) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &rowrhsscores, nrows) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &rowscores, nrows) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &roworder, nrows) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &varsolvals, nvars) );
  
   /* get the LP solution for all active variables */
   CHECK_OKAY( SCIPgetVarSols(scip, nvars, vars, varsolvals) );

   /* get the maximal number of cuts allowed in a separation round */
   if( depth == 0 )
   {
      maxtries = sepadata->maxtriesroot;
      maxaggrs = sepadata->maxaggrsroot;
      maxsepacuts = sepadata->maxsepacutsroot;
      maxslack = sepadata->maxslackroot;
      maxtestdelta = sepadata->maxtestdeltaroot;
      maxconts = sepadata->maxcontsroot;
   }   
   else
   {
      maxtries = sepadata->maxtries;
      maxaggrs = sepadata->maxaggrs;
      maxsepacuts = sepadata->maxsepacuts;
      maxslack = sepadata->maxslack;
      maxtestdelta = sepadata->maxtestdelta;
      maxconts = sepadata->maxconts;
   }

   /* calculate aggregation scores for both sides of all rows, and sort rows by nonincreasing maximal score */
   for( r = 0; r < nrows; r++ )
   {
      Real activity;
      Real lhs;
      Real rhs;
      Real lencoef;
      Real rownorm;
      Real slack;
      int i;

      assert(SCIProwGetLPPos(rows[r]) == r);

      activity = SCIPgetRowLPActivity(scip, rows[r]);
      lhs = SCIProwGetLhs(rows[r]);
      rhs = SCIProwGetRhs(rows[r]);
      rownorm = SCIProwGetNorm(rows[r]);
      rownorm = MAX(rownorm, 0.1);
      lencoef = (Real)SCIProwGetNNonz(rows[r])/(Real)nvars;
      assert(SCIPisPositive(scip, rownorm));

      slack = (activity - lhs)/rownorm;
      if( !SCIPisInfinity(scip, -lhs) && SCIPisLE(scip, slack, maxslack)
         && (ALLOWLOCAL || !SCIProwIsLocal(rows[r])) )
         rowlhsscores[r] = -lencoef - sepadata->slackscore * slack;
      else
         rowlhsscores[r] = -SCIPinfinity(scip);

      slack = (rhs - activity)/rownorm;
      if( !SCIPisInfinity(scip, rhs) && SCIPisLE(scip, slack, maxslack)
         && (ALLOWLOCAL || !SCIProwIsLocal(rows[r])) )
         rowrhsscores[r] = -lencoef - sepadata->slackscore * slack;
      else
         rowrhsscores[r] = -SCIPinfinity(scip);

      rowscores[r] = MAX(rowlhsscores[r], rowrhsscores[r]);
      for( i = r; i > 0 && rowscores[r] > rowscores[roworder[i-1]]; --i )
         roworder[i] = roworder[i-1];
      assert(0 <= i && i <= r);
      roworder[i] = r;
   }
 
   /* start aggregation heuristic for each row in the LP */
   ncuts = 0;
   ntries = (maxtries >= 0 ? MIN(maxtries, nrows) : nrows);
   for( r = 0; r < ntries && ncuts < maxsepacuts && !SCIPisInfinity(scip, -rowscores[roworder[r]]); r++ )
   {
      CHECK_OKAY( aggregation(scip, sepadata, varsolvals, rowlhsscores, rowrhsscores, roworder[r], 
            maxaggrs, maxslack, maxtestdelta, maxconts, &ncuts) );
   }

   /* free data structure */
   SCIPfreeBufferArray(scip, &varsolvals);
   SCIPfreeBufferArray(scip, &roworder);
   SCIPfreeBufferArray(scip, &rowscores);
   SCIPfreeBufferArray(scip, &rowrhsscores);
   SCIPfreeBufferArray(scip, &rowlhsscores);

   if( ncuts > 0 )
      *result = SCIP_SEPARATED;
    
   return SCIP_OKAY;
}




/*
 * separator specific interface methods
 */

/** creates the cmir separator and includes it in SCIP */
RETCODE SCIPincludeSepaCmir(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   SEPADATA* sepadata;

   /* create cmir separator data */
   CHECK_OKAY( SCIPallocMemory(scip, &sepadata) );

   /* include separator */
   CHECK_OKAY( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ,
         sepaFreeCmir, sepaInitCmir, sepaExitCmir, sepaExecCmir,
         sepadata) );

   /* add cmir separator parameters */
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/cmir/maxrounds",
         "maximal number of cmir separation rounds per node (-1: unlimited)",
         &sepadata->maxrounds, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/cmir/maxroundsroot",
         "maximal number of cmir separation rounds in the root node (-1: unlimited)",
         &sepadata->maxroundsroot, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/cmir/maxtries",
         "maximal number of rows to start aggregation with per separation round (-1: unlimited)",
         &sepadata->maxtries, DEFAULT_MAXTRIES, -1, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/cmir/maxtriesroot",
         "maximal number of rows to start aggregation with per separation round in the root node (-1: unlimited)",
         &sepadata->maxtriesroot, DEFAULT_MAXTRIESROOT, -1, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/cmir/maxaggrs",
         "maximal number of aggregations for each row per separation round",
         &sepadata->maxaggrs, DEFAULT_MAXAGGRS, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/cmir/maxaggrsroot",
         "maximal number of aggregations for each row per separation round in the root node",
         &sepadata->maxaggrsroot, DEFAULT_MAXAGGRSROOT, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/cmir/maxsepacuts",
         "maximal number of cmir cuts separated per separation round",
         &sepadata->maxsepacuts, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/cmir/maxsepacutsroot",
         "maximal number of cmir cuts separated per separation round in the root node",
         &sepadata->maxsepacutsroot, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
         "separating/cmir/maxslack",
         "maximal slack of rows to be used in aggregation",
         &sepadata->maxslack, DEFAULT_MAXSLACK, 0.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
         "separating/cmir/maxslackroot",
         "maximal slack of rows to be used in aggregation in the root node",
         &sepadata->maxslackroot, DEFAULT_MAXSLACKROOT, 0.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
         "separating/cmir/slackscore",
         "weight of slack in the aggregation scoring of the rows",
         &sepadata->slackscore, DEFAULT_SLACKSCORE, 0.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/cmir/maxrowfac",
         "maximal row aggregation factor",
         &sepadata->maxrowfac, DEFAULT_MAXROWFAC, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/cmir/minrowfac",
         "minimal row aggregation factor",
         &sepadata->minrowfac, DEFAULT_MINROWFAC, INT_MIN, 0, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/cmir/maxtestdelta",
         "maximal number of different deltas to try",
         &sepadata->maxtestdelta, DEFAULT_MAXTESTDELTA, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/cmir/maxtestdeltaroot",
         "maximal number of different deltas to try in the root node",
         &sepadata->maxtestdeltaroot, DEFAULT_MAXTESTDELTAROOT, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/cmir/maxconts",
         "maximal number of active continuous variables in aggregated row", 
         &sepadata->maxconts, DEFAULT_MAXCONTS, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/cmir/maxcontsroot",
         "maximal number of active continuous variables in aggregated row in the root node", 
         &sepadata->maxcontsroot, DEFAULT_MAXCONTSROOT, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddBoolParam(scip,
         "separating/cmir/dynamiccuts",
         "should generated cuts be removed from the LP if they are no longer tight?",
         &sepadata->dynamiccuts, DEFAULT_DYNAMICCUTS, NULL, NULL) );

   return SCIP_OKAY;
}

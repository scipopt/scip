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
#pragma ident "@(#) $Id: heur_rounding.c,v 1.24 2004/04/27 15:50:00 bzfpfend Exp $"

/**@file   heur_rounding.c
 * @brief  LP rounding heuristic that tries to recover from intermediate infeasibilities
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_rounding.h"


#define HEUR_NAME         "rounding"
#define HEUR_DESC         "LP rounding heuristic with infeasibility recovering"
#define HEUR_DISPCHAR     'R'
#define HEUR_PRIORITY     -1000
#define HEUR_FREQ         5
#define HEUR_FREQOFS      0
#define HEUR_MAXDEPTH     -1
#define HEUR_PSEUDONODES  FALSE         /** call heuristic at nodes where only a pseudo solution exist? */


/* locally defined heuristic data */
struct HeurData
{
   SOL*             sol;                /**< working solution */
};




/*
 * local methods
 */

/** update row violation arrays after a row's activity value changed */
static
void updateViolations(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   ROW**            violrows,           /**< array with currently violated rows */
   int*             violrowpos,         /**< position of LP rows in violrows array */
   int*             nviolrows,          /**< pointer to the number of currently violated rows */
   Real             oldactivity,        /**< old activity value of LP row */
   Real             newactivity         /**< new activity value of LP row */
   )
{
   Real lhs;
   Real rhs;
   Bool oldviol;
   Bool newviol;

   assert(violrows != NULL);
   assert(violrowpos != NULL);
   assert(nviolrows != NULL);

   lhs = SCIProwGetLhs(row);
   rhs = SCIProwGetRhs(row);
   oldviol = (SCIPisFeasLT(scip, oldactivity, lhs) || SCIPisFeasGT(scip, oldactivity, rhs));
   newviol = (SCIPisFeasLT(scip, newactivity, lhs) || SCIPisFeasGT(scip, newactivity, rhs));
   if( oldviol != newviol )
   {
      int rowpos;
      int violpos;
      
      rowpos = SCIProwGetLPPos(row);
      assert(rowpos >= 0);

      if( oldviol )
      {
         /* the row violation was repaired: remove row from violrows array, decrease violation count */
         violpos = violrowpos[rowpos];
         assert(0 <= violpos && violpos < *nviolrows);
         assert(violrows[violpos] == row);
         violrowpos[rowpos] = -1;
         if( violpos != *nviolrows-1 )
         {
            violrows[violpos] = violrows[*nviolrows-1];
            violrowpos[SCIProwGetLPPos(violrows[violpos])] = violpos;
         }
         (*nviolrows)--;
      }
      else
      {
         /* the row is now violated: add row to violrows array, increase violation count */
         assert(violrowpos[rowpos] == -1);
         violrows[*nviolrows] = row;
         violrowpos[rowpos] = *nviolrows;
         (*nviolrows)++;
      }
   }
}

/** update row activities after a variable's solution value changed */
static
RETCODE updateActivities(
   SCIP*            scip,               /**< SCIP data structure */
   Real*            activities,         /**< LP row activities */
   ROW**            violrows,           /**< array with currently violated rows */
   int*             violrowpos,         /**< position of LP rows in violrows array */
   int*             nviolrows,          /**< pointer to the number of currently violated rows */
   int              nlprows,            /**< number of rows in current LP */
   VAR*             var,                /**< variable that has been changed */
   Real             oldsolval,          /**< old solution value of variable */
   Real             newsolval           /**< new solution value of variable */
   )
{
   COL* col;
   ROW** colrows;
   Real* colvals;
   Real delta;
   int ncolrows;
   int r;

   assert(activities != NULL);
   assert(nviolrows != NULL);
   assert(0 <= *nviolrows && *nviolrows <= nlprows);

   delta = newsolval - oldsolval;
   col = SCIPvarGetCol(var);
   colrows = SCIPcolGetRows(col);
   colvals = SCIPcolGetVals(col);
   ncolrows = SCIPcolGetNNonz(col);
   assert(ncolrows == 0 || (colrows != NULL && colvals != NULL));
   /**@todo This is awfully slow, because the loop runs over ALL non-zero coefficients of the column, not only the
    *       ones that are currently in the LP. Solution: implement a SCIPcolGetLPRowInds() call, that returns an
    *       int array with the positions of the LP rows in the column's rows- and vals- arrays. This array
    *       should be cached and declared invalid if the LP is modified. (analoguosly for SCIProwGetLPColInds())
    */
   for( r = 0; r < ncolrows; ++r )
   {
      ROW* row;
      int rowpos;
      
      row = colrows[r];
      rowpos = SCIProwGetLPPos(row);
      assert(-1 <= rowpos && rowpos < nlprows);
      
      if( rowpos >= 0 && !SCIProwIsLocal(row) )
      {
         Real oldactivity;
         Real newactivity;
         
         assert(SCIProwIsInLP(row));
         
         /* update row activity */
         oldactivity = activities[rowpos];
         newactivity = oldactivity + delta * colvals[r];
         activities[rowpos] = newactivity;

         /* update row violation arrays */
         updateViolations(scip, row, violrows, violrowpos, nviolrows, oldactivity, newactivity);
      }            
   }

   return SCIP_OKAY;
}

/** returns a variable, that pushes activity of the row in the given direction with minimal negative impact on other rows;
 *  if variables have equal impact, chooses the one with best objective value improvement in corresponding direction;
 *  rounding in a direction is forbidden, if this forces the objective value over the upper bound
 */
static
RETCODE selectRounding(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   Real             minobj,             /**< minimal objective value possible after rounding remaining fractional vars */
   ROW*             row,                /**< LP row */
   int              direction,          /**< should the activity be increased (+1) or decreased (-1)? */
   VAR**            roundvar,           /**< pointer to store the rounding variable, returns NULL if impossible */
   Real*            oldsolval,          /**< old (fractional) solution value of rounding variable */
   Real*            newsolval           /**< new (rounded) solution value of rounding variable */
   )
{
   COL* col;
   VAR* var;
   Real val;
   COL** rowcols;
   Real* rowvals;
   Real solval;
   Real roundval;
   Real obj;
   Real deltaobj;
   Real bestdeltaobj;
   VARTYPE vartype;
   int nrowcols;
   int nlocks;
   int minnlocks;
   int c;

   assert(direction == +1 || direction == -1);
   assert(roundvar != NULL);
   assert(oldsolval != NULL);
   assert(newsolval != NULL);

   /* get row entries */
   rowcols = SCIProwGetCols(row);
   rowvals = SCIProwGetVals(row);
   nrowcols = SCIProwGetNNonz(row);

   /* select rounding variable */
   minnlocks = INT_MAX;
   bestdeltaobj = SCIPinfinity(scip);
   *roundvar = NULL;
   for( c = 0; c < nrowcols; ++c )
   {
      col = rowcols[c];
      var = SCIPcolGetVar(col);
      
      vartype = SCIPvarGetType(var);
      if( vartype == SCIP_VARTYPE_BINARY || vartype == SCIP_VARTYPE_INTEGER )
      {
         solval = SCIPgetSolVal(scip, sol, var);
      
         if( !SCIPisIntegral(scip, solval) )
         {
            val = rowvals[c];
            obj = SCIPvarGetObj(var);

            if( direction * val < 0.0 )
            {
               /* rounding down */
               nlocks = SCIPvarGetNLocksDown(var);
               if( nlocks <= minnlocks )
               {
                  roundval = SCIPfloor(scip, solval);
                  deltaobj = obj * (roundval - solval);
                  if( (nlocks < minnlocks || deltaobj < bestdeltaobj) && minobj - obj < SCIPgetUpperbound(scip) )
                  {
                     minnlocks = nlocks;
                     bestdeltaobj = deltaobj;
                     *roundvar = var;
                     *oldsolval = solval;
                     *newsolval = roundval;
                  }
               }
            }
            else
            {
               /* rounding up */
               assert(direction * val > 0.0);
               nlocks = SCIPvarGetNLocksUp(var);
               if( nlocks <= minnlocks )
               {
                  roundval = SCIPceil(scip, solval);
                  deltaobj = obj * (roundval - solval);
                  if( (nlocks < minnlocks || deltaobj < bestdeltaobj) && minobj + obj < SCIPgetUpperbound(scip) )
                  {
                     minnlocks = nlocks;
                     bestdeltaobj = deltaobj;
                     *roundvar = var;
                     *oldsolval = solval;
                     *newsolval = roundval;
                  }
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** returns a variable, that increases the activity of the row */
static
RETCODE selectIncreaseRounding(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   Real             minobj,             /**< minimal objective value possible after rounding remaining fractional vars */
   ROW*             row,                /**< LP row */
   VAR**            roundvar,           /**< pointer to store the rounding variable, returns NULL if impossible */
   Real*            oldsolval,          /**< old (fractional) solution value of rounding variable */
   Real*            newsolval           /**< new (rounded) solution value of rounding variable */
   )
{
   return selectRounding(scip, sol, minobj, row, +1, roundvar, oldsolval, newsolval);
}

/** returns a variable, that decreases the activity of the row */
static
RETCODE selectDecreaseRounding(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   Real             minobj,             /**< minimal objective value possible after rounding remaining fractional vars */
   ROW*             row,                /**< LP row */
   VAR**            roundvar,           /**< pointer to store the rounding variable, returns NULL if impossible */
   Real*            oldsolval,          /**< old (fractional) solution value of rounding variable */
   Real*            newsolval           /**< new (rounded) solution value of rounding variable */
   )
{
   return selectRounding(scip, sol, minobj, row, -1, roundvar, oldsolval, newsolval);
}

/** returns a fractional variable, that has most impact on rows in opposite direction, i.e. that is most crucial to
 *  fix in the other direction;
 *  if variables have equal impact, chooses the one with best objective value improvement in corresponding direction;
 *  rounding in a direction is forbidden, if this forces the objective value over the upper bound
 */
static
RETCODE selectEssentialRounding(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   Real             minobj,             /**< minimal objective value possible after rounding remaining fractional vars */
   VAR**            lpcands,            /**< fractional variables in LP */
   int              nlpcands,           /**< number of fractional variables in LP */
   VAR**            roundvar,           /**< pointer to store the rounding variable, returns NULL if impossible */
   Real*            oldsolval,          /**< old (fractional) solution value of rounding variable */
   Real*            newsolval           /**< new (rounded) solution value of rounding variable */
   )
{
   VAR* var;
   Real solval;
   Real roundval;
   Real obj;
   Real deltaobj;
   Real bestdeltaobj;
   int maxnlocks;
   int nlocks;
   int v;

   assert(roundvar != NULL);
   assert(oldsolval != NULL);
   assert(newsolval != NULL);

   /* select rounding variable */
   maxnlocks = -1;
   bestdeltaobj = SCIPinfinity(scip);
   *roundvar = NULL;
   for( v = 0; v < nlpcands; ++v )
   {
      var = lpcands[v];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER);
      
      solval = SCIPgetSolVal(scip, sol, var);
      if( !SCIPisIntegral(scip, solval) )
      {
         obj = SCIPvarGetObj(var);

         /* rounding down */
         nlocks = SCIPvarGetNLocksUp(var);
         if( nlocks >= maxnlocks )
         {
            roundval = SCIPfloor(scip, solval);
            deltaobj = obj * (roundval - solval);
            if( (nlocks > maxnlocks || deltaobj < bestdeltaobj) && minobj - obj < SCIPgetUpperbound(scip) )
            {
               maxnlocks = nlocks;
               bestdeltaobj = deltaobj;
               *roundvar = var;
               *oldsolval = solval;
               *newsolval = roundval;
            }
         }

         /* rounding up */
         nlocks = SCIPvarGetNLocksDown(var);
         if( nlocks >= maxnlocks )
         {
            roundval = SCIPceil(scip, solval);
            deltaobj = obj * (roundval - solval);
            if( (nlocks > maxnlocks || deltaobj < bestdeltaobj) && minobj + obj < SCIPgetUpperbound(scip) )
            {
               maxnlocks = nlocks;
               bestdeltaobj = deltaobj;
               *roundvar = var;
               *oldsolval = solval;
               *newsolval = roundval;
            }
         }
      }
   }

   return SCIP_OKAY;
}




/*
 * Callback methods
 */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
#define heurFreeRounding NULL


/** initialization method of primal heuristic (called after problem was transformed) */
static
DECL_HEURINIT(heurInitRounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(SCIPheurGetData(heur) == NULL);
   assert(scip != NULL);

   /* create heuristic data */
   CHECK_OKAY( SCIPallocMemory(scip, &heurdata) );
   CHECK_OKAY( SCIPcreateSol(scip, &heurdata->sol, heur) );
   SCIPheurSetData(heur, heurdata);

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
DECL_HEUREXIT(heurExitRounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   CHECK_OKAY( SCIPfreeSol(scip, &heurdata->sol) );
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
DECL_HEUREXEC(heurExecRounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   HEURDATA* heurdata;
   SOL* sol;
   VAR** lpcands;
   Real* lpcandssol;
   ROW** lprows;
   Real* activities;
   ROW** violrows;
   int* violrowpos;
   Real obj;
   Real bestroundval;
   Real minobj;
   Bool aborted;
   int nlpcands;
   int nlprows;
   int nfrac;
   int nviolrows;
   int c;
   int r;

   /**@todo try to shift continuous variables to stay feasible */
   /**@todo improve rounding heuristic */

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPhasActNodeLP(scip));

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* get fractional variables, that should be integral */
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands, NULL) );
   nfrac = nlpcands;

   /* only call heuristic, if LP solution is fractional */
   if( nfrac == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* get LP rows */
   CHECK_OKAY( SCIPgetLPRowsData(scip, &lprows, &nlprows) );

   debugMessage("executing rounding heuristic: %d LP rows, %d fractionals\n", nlprows, nfrac);

   /* get memory for activities, violated rows, and row violation positions */
   CHECK_OKAY( SCIPallocBufferArray(scip, &activities, nlprows) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &violrows, nlprows) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &violrowpos, nlprows) );

   /* get the activities for all globally valid rows */
   for( r = 0; r < nlprows; ++r )
   {
      ROW* row;

      row = lprows[r];
      assert(SCIProwGetLPPos(row) == r);

      if( !SCIProwIsLocal(row) )
      {
         activities[r] = SCIPgetRowActivity(scip, row);
         assert(SCIPisFeasGE(scip, activities[r], SCIProwGetLhs(row)));
         assert(SCIPisFeasLE(scip, activities[r], SCIProwGetRhs(row)));
         violrowpos[r] = -1;
      }
   }
   nviolrows = 0;

   /* get the working solution from heuristic's local data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   sol = heurdata->sol;
   assert(sol != NULL);

   /* copy the current LP solution to the working solution */
   CHECK_OKAY( SCIPlinkLPSol(scip, sol) );

   /* calculate the minimal objective value possible after rounding fractional variables */
   minobj = SCIPsolGetObj(sol);
   assert(minobj < SCIPgetUpperbound(scip));
   for( c = 0; c < nlpcands; ++c )
   {
      obj = SCIPvarGetObj(lpcands[c]);
      bestroundval = obj > 0.0 ? SCIPfloor(scip, lpcandssol[c]) : SCIPceil(scip, lpcandssol[c]);
      minobj += obj * (bestroundval - lpcandssol[c]);
   }

   /* try to round remaining variables in order to become/stay feasible */
   aborted = FALSE;
   while( nfrac > 0 && !aborted )
   {
      VAR* roundvar;
      Real oldsolval;
      Real newsolval;
         
      debugMessage("rounding heuristic: nfrac=%d, nviolrows=%d, obj=%g (best possible obj: %g)\n",
         nfrac, nviolrows, SCIPretransformObj(scip, SCIPsolGetObj(sol)), SCIPretransformObj(scip, minobj));

      assert(minobj < SCIPgetUpperbound(scip)); /* otherwise, the rounding variable selection should have returned NULL */

      /* choose next variable to process:
       *  - if a violated row exists, round a variable decreasing the violation, that has least impact on other rows
       *  - otherwise, round a variable, that has strongest devastating impact on rows in opposite direction
       */
      if( nviolrows > 0 )
      {
         ROW* row;
         int rowpos;

         row = violrows[nviolrows-1];
         rowpos = SCIProwGetLPPos(row);
         assert(0 <= rowpos && rowpos < nlprows);
         assert(violrowpos[rowpos] == nviolrows-1);

         debugMessage("rounding heuristic: try to fix violated row <%s>: %g <= %g <= %g\n",
            SCIProwGetName(row), SCIProwGetLhs(row), activities[rowpos], SCIProwGetRhs(row));
         if( SCIPisFeasLT(scip, activities[rowpos], SCIProwGetLhs(row)) )
         {
            /* lhs is violated: select a variable rounding, that increases the activity */
            CHECK_OKAY( selectIncreaseRounding(scip, sol, minobj, row, &roundvar, &oldsolval, &newsolval) );
         }
         else
         {
            assert(SCIPisFeasGT(scip, activities[rowpos], SCIProwGetRhs(row)));
            /* rhs is violated: select a variable rounding, that decreases the activity */
            CHECK_OKAY( selectDecreaseRounding(scip, sol, minobj, row, &roundvar, &oldsolval, &newsolval) );
         }
      }
      else
      {
         debugMessage("rounding heuristic: search rounding variable and try to stay feasible\n");
         CHECK_OKAY( selectEssentialRounding(scip, sol, minobj, lpcands, nlpcands, &roundvar, &oldsolval, &newsolval) );
      }

      /* check, whether rounding was possible */
      if( roundvar == NULL )
      {
         debugMessage("rounding heuristic:  -> didn't find a rounding variable\n");
         aborted = TRUE;
      }
      else
      {
         debugMessage("rounding heuristic:  -> round var <%s>, oldval=%g, newval=%g, obj=%g\n",
            SCIPvarGetName(roundvar), oldsolval, newsolval, SCIPvarGetObj(roundvar));
         
         /* update row activities of globally valid rows */
         CHECK_OKAY( updateActivities(scip, activities, violrows, violrowpos, &nviolrows, nlprows, 
                        roundvar, oldsolval, newsolval) );
         
         /* store new solution value and decrease fractionality counter */
         CHECK_OKAY( SCIPsetSolVal(scip, sol, roundvar, newsolval) );
         nfrac--;

         /* update minimal objective value possible after rounding remaining variables */
         obj = SCIPvarGetObj(roundvar);
         if( obj > 0.0 && newsolval > oldsolval )
            minobj += obj;
         else if( obj < 0.0 && newsolval < oldsolval )
            minobj -= obj;

         debugMessage("rounding heuristic:  -> nfrac=%d, nviolrows=%d, obj=%g (best possible obj: %g)\n",
            nfrac, nviolrows, SCIPretransformObj(scip, SCIPsolGetObj(sol)), SCIPretransformObj(scip, minobj));
      }
   }

   /* check, if the new solution is feasible */
   if( nfrac == 0 && nviolrows == 0 )
   {
      Bool stored;

      assert(!aborted);

      /* check solution for feasibility, and add it to solution store if possible
       * neither integrality nor feasibility of LP rows has to be checked, because this is already
       * done in the rounding heuristic itself
       */
      CHECK_OKAY( SCIPtrySol(scip, sol, FALSE, FALSE, &stored) );

      if( stored )
      {
#ifdef DEBUG
         printf("found feasible rounded solution:\n");
         SCIPprintSol(scip, sol, NULL);
#endif
         *result = SCIP_FOUNDSOL;
      }
   }

   /* free memory buffers */
   CHECK_OKAY( SCIPfreeBufferArray(scip, &violrowpos) );
   CHECK_OKAY( SCIPfreeBufferArray(scip, &violrows) );
   CHECK_OKAY( SCIPfreeBufferArray(scip, &activities) );
   
   return SCIP_OKAY;
}




/*
 * heuristic specific interface methods
 */

/** creates the rounding heuristic with infeasibility recovering and includes it in SCIP */
RETCODE SCIPincludeHeurRounding(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   /* include heuristic */
   CHECK_OKAY( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
                  HEUR_MAXDEPTH, HEUR_PSEUDONODES,
                  heurFreeRounding, heurInitRounding, heurExitRounding, heurExecRounding,
                  NULL) );

   return SCIP_OKAY;
}


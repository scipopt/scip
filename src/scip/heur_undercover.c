/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_undercover.c
 * @ingroup PRIMALHEURISTICS
 * @brief  undercover primal heuristic for MIQCPs
 *
 * @author Timo Berthold
 * @author Ambros Gleixner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/heur_undercover.h"

#define HEUR_NAME               "undercover"
#define HEUR_DESC               "solves a sub-CIP determined by a set covering approach"
#define HEUR_DISPCHAR           'U'
#define HEUR_PRIORITY           -1110000
#define HEUR_FREQ               0
#define HEUR_FREQOFS            0
#define HEUR_MAXDEPTH           -1
#define HEUR_TIMING             SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP        TRUE         /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_FIXINGALTS      "li"         /**< sequence of fixing values used: 'l'p relaxation, 'n'lp relaxation, 'i'ncumbent solution */
#define DEFAULT_MAXNODES        (SCIP_Longint)500/**< maximum number of nodes to regard in the subproblem */
#define DEFAULT_MINNODES        (SCIP_Longint)500/**< minimum number of nodes to regard in the subproblem */
#define DEFAULT_NODESOFS        (SCIP_Longint)500/**< number of nodes added to the contingent of the total nodes */
#define DEFAULT_CONFLICTWEIGHT  1000.0       /**< weight for conflict score in fixing order */
#define DEFAULT_CUTOFFWEIGHT    1.0          /**< weight for cutoff score in fixing order */
#define DEFAULT_INFERENCEWEIGHT 1.0          /**< weight for inference score in foxing order */
#define DEFAULT_MAXCOVERSIZE    1.0          /**< maximum coversize (as fraction of total number of variables) */
#define DEFAULT_MINIMPROVE      0.0          /**< factor by which heuristic should at least improve the incumbent */
#define DEFAULT_NODESQUOT       0.1          /**< subproblem nodes in relation to nodes of the original problem */
#define DEFAULT_RECOVERDIV      0.9          /**< fraction of covering variables in the last cover which need to change their value when re-covering */
#define DEFAULT_BEFORECUTS      TRUE         /**< should undercover called at root node before cut separation? */
#define DEFAULT_FIXINTFIRST     FALSE        /**< should integer variables in the cover be fixed first? */
#define DEFAULT_LOCKSROUNDING   TRUE         /**< shall LP values for integer vars be rounded according to locks? */
#define DEFAULT_ONLYCONVEXIFY   FALSE        /**< should we only fix/dom.red. variables creating nonconvexity? */
#define DEFAULT_POSTNLP         TRUE         /**< should the nlp heuristic be called to polish a feasible solution? */
#define DEFAULT_MAXBACKTRACKS   6            /**< maximum number of backtracks */
#define DEFAULT_MAXRECOVERS     1            /**< maximum number of re-coverings */
#define DEFAULT_MAXREORDERS     1            /**< maximum number of reorderings of the fixing order */
#define DEFAULT_COVERINGOBJ     'u'          /**< objective function of the covering problem */
#define DEFAULT_COPYCUTS        TRUE         /**< should all active cuts from the cutpool of the original scip be copied
                                              *   to constraints of the subscip
					      */

#define COVERINGOBJS            "cdlmtu"     /**< list of objective functions of the covering problem */
#define MAXNLPFAILS             1            /**< maximum number of fails after which we give up solving the nlp relaxation */
#define MAXPOSTNLPFAILS         1            /**< maximum number of fails after which we give up calling nlp local search */
#define MINTIMELEFT             1.0          /**< don't start expensive parts of the heuristics if less than this amount of time left */
#define SUBMIPSETUPCOSTS        200          /**< number of nodes equivalent for the costs for setting up the sub-CIP */


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_CONSHDLR**       nlconshdlrs;        /**< array of nonlinear constraint handlers */
   SCIP_HEUR*            nlpheur;            /**< pointer to nlp local search heuristics */
   char*                 fixingalts;         /**< sequence of fixing values used: 'l'p relaxation, 'n'lp relaxation, 'i'ncumbent solution */
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes */
   SCIP_Longint          nusednodes;         /**< nodes already used by heuristic in earlier calls */
   SCIP_Real             conflictweight;     /**< weight for conflict score in fixing order */
   SCIP_Real             cutoffweight;       /**< weight for cutoff score in fixing order */
   SCIP_Real             inferenceweight;    /**< weight for inference score in foxing order */
   SCIP_Real             maxcoversize;       /**< maximum coversize (as fraction of total number of variables) */
   SCIP_Real             minimprove;         /**< factor by which heuristic should at least improve the incumbent */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem */
   SCIP_Real             recoverdiv;         /**< fraction of covering variables in the last cover which need to change their value when re-covering */
   SCIP_Bool             beforecuts;         /**< should undercover be called at root node before cut separation? */
   SCIP_Bool             fixintfirst;        /**< should integer variables in the cover be fixed first? */
   SCIP_Bool             globalbounds;       /**< should global bounds on variables be used instead of local bounds at focus node? */
   SCIP_Bool             locksrounding;      /**< shall LP values for integer vars be rounded according to locks? */
   SCIP_Bool             nlpsolved;          /**< has current nlp relaxation already been solved successfully? */
   SCIP_Bool             nlpfailed;          /**< has solving the nlp relaxation failed? */
   SCIP_Bool             onlyconvexify;      /**< should we only fix/dom.red. variables creating nonconvexity? */
   SCIP_Bool             postnlp;            /**< should the nlp heuristic be called to polish a feasible solution? */
   int                   maxbacktracks;      /**< maximum number of backtracks */
   int                   maxrecovers;        /**< maximum number of re-coverings */
   int                   maxreorders;        /**< maximum number of reorderings of the fixing order */
   int                   nnlpfails;          /**< number of fails when solving the nlp relaxation after last success */
   int                   npostnlpfails;      /**< number of fails of the nlp local search after last success */
   int                   nnlconshdlrs;       /**< number of nonlinear constraint handlers */
   char                  coveringobj;        /**< objective function of the covering problem */
   SCIP_Bool             copycuts;           /**< should all active cuts from cutpool be copied to constraints in
                                              *   subproblem?
                                              */
};


/*
 * Local methods
 */


/** determines, whether a variable is fixed to the given value */
static
SCIP_Bool varIsFixed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to check */
   SCIP_Real             val,                /**< value to check */
   SCIP_Bool             global              /**< should global bounds be used? */
   )
{
   SCIP_Bool isfixed;

   if( global )
      isfixed = SCIPisFeasEQ(scip, val, SCIPvarGetLbGlobal(var)) && SCIPisFeasEQ(scip, val, SCIPvarGetUbGlobal(var));
   else
      isfixed = SCIPisFeasEQ(scip, val, SCIPvarGetLbLocal(var)) && SCIPisFeasEQ(scip, val, SCIPvarGetUbLocal(var));

   return isfixed;
}


/** determines, whether a term is already constant, because the variable is fixed or the coefficient is zero */
static
SCIP_Bool termIsConstant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to check */
   SCIP_Real             coeff,              /**< coefficient to check */
   SCIP_Bool             global              /**< should global bounds be used? */
   )
{
   /* if the variable has zero coefficient in the original problem, the term is linear */
   if( SCIPisZero(scip, coeff) )
      return TRUE;

   /* if the variable is fixed in the original problem, the term is linear */
   if( global )
      return SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
   else
      return SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));

}


/** determines, whether a term is convex */
static
SCIP_Bool termIsConvex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             lhs,                /**< left hand side of the constraint */
   SCIP_Real             rhs,                /**< right hand side of the constraint */
   SCIP_Bool             sign                /**< signature of the term */
   )
{
   return sign ? SCIPisInfinity(scip, -lhs) : SCIPisInfinity(scip, rhs);
}


/** increases counters */
static
void  incCounters(
   int*                  termcounter,        /**< array to count in how many nonlinear terms a variable appears */
   int*                  conscounter,        /**< array to count in how many constraints a variable appears */
   SCIP_Bool*            consmarker,         /**< was this variable already counted for this constraint? */
   int                   idx                 /**< problem index of the variable */
   )
{
   termcounter[idx]++;
   if( !consmarker[idx] )
   {
      conscounter[idx]++;
      consmarker[idx] = TRUE;
   }
   return;
}


/** update time limit */
static
SCIP_RETCODE updateTimelimit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck,               /**< clock timer */
   SCIP_Real*            timelimit           /**< time limit */
   )
{
   *timelimit -= SCIPgetClockTime(scip, clck);
   SCIP_CALL( SCIPresetClock(scip, clck) );
   SCIP_CALL( SCIPstartClock(scip, clck) );

   return SCIP_OKAY;
}


/** analyzes a nonlinear row and adds constraints and fixings to the covering problem */
static
SCIP_RETCODE processNlRow(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_NLROW*           nlrow,              /**< nonlinear row representation of a nonlinear constraint */
   SCIP*                 coveringscip,       /**< SCIP data structure for the covering problem */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            coveringvars,       /**< array to store the covering problem's variables */
   int*                  termcounter,        /**< counter array for number of nonlinear nonzeros per variable */
   int*                  conscounter,        /**< counter array for number of nonlinear constraints per variable */
   SCIP_Bool*            consmarker,         /**< marker array if constraint has been counted in conscounter */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully */
   )
{
   SCIP_EXPRTREE* exprtree;
   SCIP_Bool infeas;
   SCIP_Bool fixed;
   int t;
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(nlrow != NULL);
   assert(coveringscip != NULL);
   assert(nvars >= 1);
   assert(coveringvars != NULL);
   assert(termcounter != NULL);
   assert(conscounter != NULL);
   assert(consmarker != NULL);
   assert(success != NULL);

   *success = FALSE;
   BMSclearMemoryArray(consmarker, nvars);

   /* go through expression tree */
   exprtree = SCIPnlrowGetExprtree(nlrow);
   if( exprtree != NULL )
   {
      SCIP_VAR** exprtreevars;
      int nexprtreevars;
      int probidx;

      /* get variables in expression tree */
      nexprtreevars = SCIPexprtreeGetNVars(exprtree);
      exprtreevars = SCIPexprtreeGetVars(exprtree);

      /* currently, we fix all variables contained in the expression tree; this is possibly unnecessarily restrictive,
       * but ensures a linear resp. convex subproblem */
      if( exprtreevars != NULL )
      {
         int i;

         for( i = nexprtreevars-1; i >= 0; i-- )
         {
            assert(exprtreevars[i] != NULL);

            /* if constraints with inactive variables are present, we will have difficulties creating the sub-CIP later */
            probidx = SCIPvarGetProbindex(exprtreevars[i]);
            if( probidx == -1 )
            {
               SCIPdebugMessage("inactive variables detected in nonlinear row <%s>\n", SCIPnlrowGetName(nlrow));
               return SCIP_OKAY;
            }

            /* term is constant, nothing to do */
            if( termIsConstant(scip, exprtreevars[i], 1.0, heurdata->globalbounds) )
               continue;

            /* otherwise fix variable */
            SCIP_CALL( SCIPfixVar(coveringscip, coveringvars[probidx], 1.0, &infeas, &fixed) );
            assert(!infeas);
            assert(fixed);

            /* update counters */
            incCounters(termcounter, conscounter, consmarker, probidx);

            SCIPdebugMessage("fixing var <%s> in covering problem to 1\n", SCIPvarGetName(coveringvars[probidx]));
         }
      }
   }

   /* go through all quadratic terms */
   for( t = SCIPnlrowGetNQuadElems(nlrow)-1; t >= 0; t-- )
   {
      SCIP_QUADELEM* quadelem;
      SCIP_VAR* bilinvar1;
      SCIP_VAR* bilinvar2;
      int probidx1;
      int probidx2;

      /* get quadratic term */
      quadelem = &SCIPnlrowGetQuadElems(nlrow)[t];

      /* get involved variables */
      bilinvar1 = SCIPnlrowGetQuadVars(nlrow)[quadelem->idx1];
      bilinvar2 = SCIPnlrowGetQuadVars(nlrow)[quadelem->idx2];
      assert(bilinvar1 != NULL);
      assert(bilinvar2 != NULL);

      /* if constraints with inactive variables are present, we will have difficulties creating the sub-CIP later */
      probidx1 = SCIPvarGetProbindex(bilinvar1);
      probidx2 = SCIPvarGetProbindex(bilinvar2);
      if( probidx1 == -1 || probidx2 == -1 )
      {
         SCIPdebugMessage("inactive variables detected in nonlinear row <%s>\n", SCIPnlrowGetName(nlrow));
         return SCIP_OKAY;
      }

      /* we have a square term */
      if( bilinvar1 == bilinvar2 )
      {
         /* term is constant, nothing to do */
         if( termIsConstant(scip, bilinvar1, quadelem->coef, heurdata->globalbounds) )
            continue;

         /* if we only convexify and term is convex considering the bounds of the nlrow, nothing to do */
         if( heurdata->onlyconvexify && termIsConvex(scip, SCIPnlrowGetLhs(nlrow), SCIPnlrowGetRhs(nlrow), quadelem->coef >= 0) )
            continue;

         /* otherwise variable has to be in the cover */
         SCIP_CALL( SCIPfixVar(coveringscip, coveringvars[probidx1], 1.0, &infeas, &fixed) );
         assert(!infeas);
         assert(fixed);

         /* update counters */
         incCounters(termcounter, conscounter, consmarker, probidx1);

         SCIPdebugMessage("fixing var <%s> in covering problem to 1\n", SCIPvarGetName(coveringvars[probidx1]));
      }
      /* we have a bilinear term */
      else
      {
         SCIP_CONS* coveringcons;
         SCIP_VAR* coveringconsvars[2];

         /* if the term is linear because one of the variables is fixed or the coefficient is zero, nothing to do */
         if( termIsConstant(scip, bilinvar1, quadelem->coef, heurdata->globalbounds)
            || termIsConstant(scip, bilinvar2, quadelem->coef, heurdata->globalbounds) )
            continue;

         /* create covering constraint */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_covering%d", SCIPnlrowGetName(nlrow), t);
         coveringconsvars[0] = coveringvars[probidx1];
         coveringconsvars[1] = coveringvars[probidx2];
         SCIP_CALL( SCIPcreateConsSetcover(coveringscip, &coveringcons, name, 2, coveringconsvars,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

         if( coveringcons == NULL )
         {
            SCIPdebugMessage("failed to create set covering constraint <%s>\n", name);
            return SCIP_OKAY;
         }

         /* add and release covering constraint */
         SCIP_CALL( SCIPaddCons(coveringscip, coveringcons) );
         SCIP_CALL( SCIPreleaseCons(coveringscip, &coveringcons) );

         /* update counters for both variables */
         incCounters(termcounter, conscounter, consmarker, probidx1);
         incCounters(termcounter, conscounter, consmarker, probidx2);
      }
   }

   *success = TRUE;

   return SCIP_OKAY;
}


/** creates the covering problem to determine a number of variables to be fixed */
static
SCIP_RETCODE createCoveringProblem(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 coveringscip,       /**< SCIP data structure for the covering problem */
   SCIP_VAR**            coveringvars,       /**< array to store the covering problem's variables */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully */
   )
{
   SCIP_VAR** vars;
   SCIP_CONSHDLR* conshdlr;
   SCIP_HASHMAP* nlrowmap;
   SCIP_Bool* consmarker;
   int* conscounter;
   int* termcounter;

   int nlocksup;
   int nlocksdown;
   int nvars;
   int i;
   int probindex;
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(coveringscip != NULL);
   assert(coveringvars != NULL);
   assert(heurdata != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* create problem data structure */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_covering", SCIPgetProbName(scip));
   SCIP_CALL( SCIPcreateProb(coveringscip, name, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* allocate and initialize to zero counter arrays for weighted objectives */
   SCIP_CALL( SCIPallocBufferArray(scip, &consmarker, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conscounter, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termcounter, nvars) );
   BMSclearMemoryArray(conscounter, nvars);
   BMSclearMemoryArray(termcounter, nvars);

   /* create covering variable for each variable in the original problem (fix it or not?) in the same order as in the
      original problem */
   for( i = 0; i < nvars; i++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_covering", SCIPvarGetName(vars[i]));
      SCIP_CALL( SCIPcreateVar(coveringscip, &coveringvars[i], name, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY,
            TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      assert(coveringvars[i] != NULL);      
      SCIP_CALL( SCIPaddVar(coveringscip, coveringvars[i]) );
   }

   /* first, go through some special constraint handlers which we do not want to treat by looking at their nlrow
    * representation; we store these in a hash map and afterwards process all nlrows which are not found in the hash map */
   nlrowmap = NULL;
   if( SCIPisNLPConstructed(scip) )
   {
      int nnlprows;

      assert(SCIPgetNLP(scip) != NULL);
 
      nnlprows = SCIPnlpGetNNlRows(SCIPgetNLP(scip));
      if( nnlprows > 0 )
      {
         int mapsize;

         /* calculate size of hash map */
         conshdlr = SCIPfindConshdlr(scip, "and");
         mapsize = SCIPconshdlrGetNConss(conshdlr);
         conshdlr = SCIPfindConshdlr(scip, "quadratic");
         mapsize += SCIPconshdlrGetNConss(conshdlr);
         conshdlr = SCIPfindConshdlr(scip, "soc");
         mapsize += SCIPconshdlrGetNConss(conshdlr);
         mapsize = MAX(mapsize, nnlprows);
         mapsize = SCIPcalcHashtableSize(2*mapsize);
         assert(mapsize > 0);

         /* create hash map */
         SCIP_CALL( SCIPhashmapCreate(&nlrowmap, SCIPblkmem(scip), mapsize) );
         assert(nlrowmap != NULL);
      }
   }

   /* go through all "and" constraints in the original problem */
   conshdlr = SCIPfindConshdlr(scip, "and");
   if( conshdlr != NULL )
   {
      int c;

      for( c = SCIPconshdlrGetNConss(conshdlr)-1; c >= 0; c-- )
      {
         SCIP_CONS* andcons;
         SCIP_CONS* coveringcons;
         SCIP_VAR** andvars;
         SCIP_VAR** coveringconsvars;
         SCIP_Real* coveringconsvals;

         int ntofix;
         int v;

         /* get original constraint and variables */
         andcons = SCIPconshdlrGetConss(conshdlr)[c];
         assert(andcons != NULL);
         andvars = SCIPgetVarsAnd(scip, andcons);
         assert(andvars != NULL);

         /* "and" constraints are not passed to the nlp, hence nothing to store in the hash map */

         /* allocate memory for covering constraint */
         SCIP_CALL( SCIPallocBufferArray(coveringscip, &coveringconsvars, SCIPgetNVarsAnd(scip, andcons)+1) );
         SCIP_CALL( SCIPallocBufferArray(coveringscip, &coveringconsvals, SCIPgetNVarsAnd(scip, andcons)+1) );

         /* collect unfixed variables */
         BMSclearMemoryArray(consmarker, nvars);
         ntofix = 0;
         for( v = SCIPgetNVarsAnd(scip, andcons)-1; v >= 0; v-- )
         {
            assert(andvars[v] != NULL);

            /* if constraints with inactive variables are present, we will have difficulties creating the sub-CIP later */
            probindex = SCIPvarGetProbindex(andvars[v]);
            if( probindex == -1 )
            {
               SCIPdebugMessage("inactive variables detected in constraint <%s>\n", SCIPconsGetName(andcons));
               SCIPfreeBufferArray(coveringscip, &coveringconsvals);
               SCIPfreeBufferArray(coveringscip, &coveringconsvars);
               goto TERMINATE;
            }

            /* if variable is fixed to 0, entire constraint can be linearized */
            if( varIsFixed(scip, andvars[v], 0.0, heurdata->globalbounds) )
            {
               ntofix = 0;
               break;
            }

            /* if variable is fixed, nothing to do */
            if( termIsConstant(scip, andvars[v], 1.0, heurdata->globalbounds) )
            {
               continue;
            }

            /* add covering variable for unfixed original variable */
            coveringconsvars[ntofix] = coveringvars[probindex];
            coveringconsvals[ntofix] = 1.0;
            ntofix++;
         }

         /* if constraints with inactive variables are present, we will have difficulties creating the sub-CIP later */
         probindex = SCIPvarGetProbindex(SCIPgetResultantAnd(scip, andcons));
         if( probindex == -1 )
         {
            SCIPdebugMessage("inactive variables detected in constraint <%s>\n", SCIPconsGetName(andcons));
            SCIPfreeBufferArray(coveringscip, &coveringconsvals);
            SCIPfreeBufferArray(coveringscip, &coveringconsvars);
            goto TERMINATE;
         }

         /* if less than 2 variables are unfixed or the resultant variable is fixed, the entire constraint can be linearized anyway */
         if( ntofix >= 2 && !termIsConstant(scip, vars[probindex], 1.0, heurdata->globalbounds) )
         {
            assert(ntofix <= SCIPgetNVarsAnd(scip, andcons));

            /* add covering variable for unfixed resultant */
            coveringconsvars[ntofix] = coveringvars[probindex];
            coveringconsvals[ntofix] = (SCIP_Real)(ntofix - 1);
            ntofix++;

            /* create covering constraint */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_covering", SCIPconsGetName(andcons));
            SCIP_CALL( SCIPcreateConsLinear(coveringscip, &coveringcons, name, ntofix, coveringconsvars, coveringconsvals,
                  (SCIP_Real)(ntofix - 2), SCIPinfinity(coveringscip),
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

            if( coveringcons == NULL )
            {
               SCIPdebugMessage("failed to create linear constraint <%s>\n", name);
               SCIPfreeBufferArray(coveringscip, &coveringconsvals);
               SCIPfreeBufferArray(coveringscip, &coveringconsvars);
               goto TERMINATE;
            }

            /* add and release covering constraint */
            SCIP_CALL( SCIPaddCons(coveringscip, coveringcons) );
            SCIP_CALL( SCIPreleaseCons(coveringscip, &coveringcons) );

            /* update counters */
            for( v = ntofix-1; v >= 0; v-- )
               incCounters(termcounter, conscounter, consmarker, SCIPvarGetProbindex(coveringconsvars[v]));
         }

         /* free memory for covering constraint */
         SCIPfreeBufferArray(coveringscip, &coveringconsvals);
         SCIPfreeBufferArray(coveringscip, &coveringconsvars);
      }
   }

   /* go through all "quadratic" constraints in the original problem */
   conshdlr = SCIPfindConshdlr(scip, "quadratic");
   if( conshdlr != NULL )
   {
      int c;

      for( c = SCIPconshdlrGetNConss(conshdlr)-1; c >= 0; c-- )
      {
         SCIP_CONS* quadcons;
         SCIP_NLROW* nlrow;

         /* get original constraint */
         quadcons = SCIPconshdlrGetConss(conshdlr)[c];
         assert(quadcons != NULL);

         /* get nlrow representation and store it in hash map */
         SCIP_CALL( SCIPgetNlRowQuadratic(scip, quadcons, &nlrow) );
         assert(nlrow != NULL);
         if( nlrowmap != NULL )
         {
            assert(!SCIPhashmapExists(nlrowmap, nlrow));
            SCIP_CALL( SCIPhashmapInsert(nlrowmap, nlrow, quadcons) );
         }

         /* if we only want to convexify and curvature and bounds prove already convexity, nothing to do */
         if( heurdata->onlyconvexify
            && ((SCIPisInfinity(scip, -SCIPgetLhsQuadratic(scip, quadcons)) && SCIPisConvexQuadratic(scip, quadcons))
               || (SCIPisInfinity(scip, SCIPgetRhsQuadratic(scip, quadcons)) && SCIPisConcaveQuadratic(scip, quadcons))) )
            continue;

         /* process nlrow */
         *success = FALSE;
         SCIP_CALL( processNlRow(scip, heurdata, nlrow, coveringscip, nvars, coveringvars, termcounter, conscounter, consmarker, success) );

         if( *success == FALSE )
            goto TERMINATE;
      }

      *success = FALSE;
   }

   /* go through all "soc" constraints in the original problem */
   conshdlr = SCIPfindConshdlr(scip, "soc");
   if( conshdlr != NULL && !heurdata->onlyconvexify )
   {
      int c;

      for( c = SCIPconshdlrGetNConss(conshdlr)-1; c >= 0; c-- )
      {
         SCIP_CONS* soccons;
         SCIP_CONS* coveringcons;
         SCIP_VAR** soclhsvars;
         SCIP_VAR* socrhsvar;
         SCIP_VAR** coveringconsvars;
         SCIP_NLROW* nlrow;

         int ntofix;
         int v;

         /* get original constraints and variables */
         soccons = SCIPconshdlrGetConss(conshdlr)[c];
         assert(soccons != NULL);
         socrhsvar = SCIPgetRhsVarSOC(scip, soccons);
         assert(socrhsvar != NULL);
         soclhsvars = SCIPgetLhsVarsSOC(scip, soccons);
         assert(soclhsvars != NULL);

         /* get nlrow representation and store it in hash map */
         SCIP_CALL( SCIPgetNlRowSOC(scip, soccons, &nlrow) );
         assert(nlrow != NULL);
         if( nlrowmap != NULL )
         {
            assert(!SCIPhashmapExists(nlrowmap, nlrow));
            SCIP_CALL( SCIPhashmapInsert(nlrowmap, nlrow, soccons) );
         }

         /* allocate memory for covering constraint */
         SCIP_CALL( SCIPallocBufferArray(coveringscip, &coveringconsvars, SCIPgetNLhsVarsSOC(scip, soccons)+1) );

         /* collect unfixed variables */
         BMSclearMemoryArray(consmarker, nvars);
         ntofix = 0;

         /* if constraints with inactive variables are present, we will have difficulties creating the sub-CIP later */
         probindex = SCIPvarGetProbindex(socrhsvar);
         if( probindex == -1 )
         {
            SCIPdebugMessage("inactive variables detected in constraint <%s>\n", SCIPconsGetName(soccons));
            SCIPfreeBufferArray(coveringscip, &coveringconsvars);
            goto TERMINATE;
         }

         /* add covering variable for unfixed rhs variable */
         if( !termIsConstant(scip, socrhsvar, SCIPgetRhsCoefSOC(scip, soccons), heurdata->globalbounds) )
         {
            SCIP_CALL( SCIPgetNegatedVar(coveringscip, coveringvars[probindex], &coveringconsvars[ntofix]) );
            ntofix++;
         }

         /* go through lhs variables */
         for( v = SCIPgetNLhsVarsSOC(scip, soccons)-1; v >= 0; v-- )
         {
            assert(soclhsvars[v] != NULL);

            /* if constraints with inactive variables are present, we will have difficulties creating the sub-CIP later */
            probindex = SCIPvarGetProbindex(soclhsvars[v]);
            if( probindex == -1 )
            {
               SCIPdebugMessage("inactive variables detected in constraint <%s>\n", SCIPconsGetName(soccons));
               SCIPfreeBufferArray(coveringscip, &coveringconsvars);
               goto TERMINATE;
            }

            /* add covering variable for unfixed lhs variable */
            if( !termIsConstant(scip, soclhsvars[v], SCIPgetLhsCoefsSOC(scip, soccons)[v], heurdata->globalbounds) )
            {
               SCIP_CALL( SCIPgetNegatedVar(coveringscip, coveringvars[probindex], &coveringconsvars[ntofix]) );
               ntofix++;
            }
         }

         if( ntofix >= 2 )
         {
            /* create covering constraint */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_covering", SCIPconsGetName(soccons));
            SCIP_CALL( SCIPcreateConsSetpack(coveringscip, &coveringcons, name, ntofix, coveringconsvars,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

            if( coveringcons == NULL )
            {
               SCIPdebugMessage("failed to create set packing constraint <%s>\n", name);
               SCIPfreeBufferArray(coveringscip, &coveringconsvars);
               goto TERMINATE;
            }

            /* add and release covering constraint */
            SCIP_CALL( SCIPaddCons(coveringscip, coveringcons) );
            SCIP_CALL( SCIPreleaseCons(coveringscip, &coveringcons) );

            /* update counters */
            for( v = ntofix-1; v >= 0; v-- )
               incCounters(termcounter, conscounter, consmarker, SCIPvarGetProbindex(SCIPvarGetNegatedVar(coveringconsvars[v])));
         }

         /* free memory for covering constraint */
         SCIPfreeBufferArray(coveringscip, &coveringconsvars);
      }
   }

   /* go through all yet unprocessed nlrows */
   if( nlrowmap != NULL )
   {
      SCIP_NLP* nlp;
      SCIP_NLROW** nlrows;
      int nnlrows;

      assert(SCIPisNLPConstructed(scip));

      /* get nlp */
      nlp = SCIPgetNLP(scip);
      assert(nlp != NULL);

      /* get nlrows */
      nnlrows = SCIPnlpGetNNlRows(nlp);
      nlrows = SCIPnlpGetNlRows(nlp);

      for( i = nnlrows-1; i >= 0; i-- )
      {
         assert(nlrows[i] != NULL);

         /* nlrow or corresponding constraint already processed */
         if( SCIPhashmapExists(nlrowmap, nlrows[i]) )
            continue;
         
         /* process nlrow */
         *success = FALSE;
         SCIP_CALL( processNlRow(scip, heurdata, nlrows[i], coveringscip, nvars, coveringvars, termcounter, conscounter, consmarker, success) );

         if( *success == FALSE )
            goto TERMINATE;
      }
   }

   /* set objective function of covering problem */
   switch( heurdata->coveringobj )
   {
   case 'c': /* number of influenced nonlinear constraints */
      for( i = nvars-1; i >= 0; i-- )
      {
         SCIP_CALL( SCIPchgVarObj(coveringscip, coveringvars[i], (SCIP_Real) conscounter[i]) );
      }
      break;
   case 'd': /* domain size */
      for( i = nvars-1; i >= 0; i-- )
      {
         SCIP_CALL( SCIPchgVarObj(coveringscip, coveringvars[i],
               (heurdata->globalbounds ? SCIPvarGetUbGlobal(vars[i]) - SCIPvarGetLbGlobal(vars[i]) : SCIPvarGetUbLocal(vars[i]) - SCIPvarGetLbLocal(vars[i]))) );
      }
      break;
   case 'l': /* number of locks */
      for( i = nvars-1; i >= 0; i-- )
      {
         nlocksup = SCIPvarGetNLocksUp(vars[i]);
         nlocksdown = SCIPvarGetNLocksDown(vars[i]);
         SCIP_CALL( SCIPchgVarObj(coveringscip, coveringvars[i], (SCIP_Real) (nlocksup+nlocksdown+1)) );
      }
      break;
   case 'm': /* min(up locks, down locks)+1 */
      for( i = nvars-1; i >= 0; i-- )
      {
         nlocksup = SCIPvarGetNLocksUp(vars[i]);
         nlocksdown = SCIPvarGetNLocksDown(vars[i]);
         SCIP_CALL( SCIPchgVarObj(coveringscip, coveringvars[i], (SCIP_Real) (MIN(nlocksup, nlocksdown)+1)) );
      }
      break;
   case 't': /* number of influenced nonlinear terms */
      for( i = nvars-1; i >= 0; i-- )
      {
         SCIP_CALL( SCIPchgVarObj(coveringscip, coveringvars[i], (SCIP_Real) termcounter[i]) );
      }
      break;
   case 'u': /* unit penalties */
      for( i = nvars-1; i >= 0; i-- )
      {
         SCIP_CALL( SCIPchgVarObj(coveringscip, coveringvars[i], 1.0) );
      }
      break;
   default:
      SCIPerrorMessage("invalid choice <%c> for covering objective\n", heurdata->coveringobj);
      goto TERMINATE;
   }

   /* covering problem successfully set up */
   *success = TRUE;

 TERMINATE:
   /* free nlrow hash map */
   if( nlrowmap != NULL )
   {
      SCIPhashmapFree(&nlrowmap);
   }

   /* free counter arrays for weighted objectives */
   SCIPfreeBufferArray(scip, &termcounter); 
   SCIPfreeBufferArray(scip, &conscounter);
   SCIPfreeBufferArray(scip, &consmarker);
 
   return SCIP_OKAY;
}


/** adds a constraint to the covering problem to forbid the given cover */
static
SCIP_RETCODE forbidCover(
   SCIP*                 scip,               /**< SCIP data structure of the covering problem */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            vars,               /**< variable array */
   int                   coversize,          /**< size of the cover */
   int*                  cover,              /**< problem indices of the variables in the cover */
   int                   diversification,    /**< how many unfixed variables have to change their value? */
   SCIP_Bool*            success,            /**< pointer to store whether the cutoff constraint was created successfully */
   SCIP_Bool*            infeas              /**< pointer to store whether the cutoff proves (local or global) infeasibility */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR** consvars;

   char name[SCIP_MAXSTRLEN];
   int nconsvars;
   int i;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nvars >= 1);
   assert(cover != NULL);
   assert(coversize >= 1);
   assert(coversize <= nvars);
   assert(diversification >= 1);
   assert(success != NULL);
   assert(infeas != NULL);

   *success = FALSE;
   *infeas = FALSE;

   /* allocate memory for constraint variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, coversize) );
   nconsvars = 0;
   cons = NULL;

   /* create constraint name */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "forbid_cover_assignment");

   /* if all variables in the cover are binary and we require only one variable to change its value, then we create a
    * set covering constraint */
   if( diversification == 1 )
   {
      /* build up constraint */
      for( i = coversize-1; i >= 0; i-- )
      {
         if( !SCIPisFeasGE(scip, SCIPvarGetLbLocal(vars[cover[i]]), 1.0) )
         {
            SCIP_CALL( SCIPgetNegatedVar(scip, vars[cover[i]], &consvars[nconsvars]) );
            nconsvars++;
         }
      }

      /* if all covering variables are fixed to one, the constraint cannot be satisfied */
      if( nconsvars == 0 )
      {
         *infeas = TRUE;
      }
      else
      {
         /* create constraint */
         SCIP_CALL( SCIPcreateConsSetcover(scip, &cons, name, nconsvars, consvars,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
      }
   }
   /* if all variables in the cover are binary and we require more variables to change their value, then we create a
    * linear constraint */
   else
   {
      SCIP_Real* consvals;
      SCIP_Real rhs;

      /* build up constraint */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvals, coversize) );
      for( i = coversize-1; i >= 0; i-- )
      {
         if( !SCIPisFeasGE(scip, SCIPvarGetLbLocal(vars[cover[i]]), 1.0) )
         {
            consvars[nconsvars] = vars[cover[i]];
            consvals[nconsvars] = 1.0;
            nconsvars++;
         }
      }
      rhs = (SCIP_Real) (nconsvars-diversification);

      /* if too many covering variables are fixed to 1, the constraint cannot be sataisfied */
      if( rhs < 0 )
      {
         *infeas = TRUE;
      }
      else
      {
         /* create constraint */
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name,
               nconsvars, consvars, consvals, -SCIPinfinity(scip), rhs,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
      }

      /* free memory */
      SCIPfreeBufferArray(scip, &consvals);
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &consvars);

   /* if proven infeasible, we do not even add the constraint; otherwise we add and release the constraint if created
    * successfully */
   if( !(*infeas) && cons != NULL )
   {
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      *success = TRUE;
   }

   return SCIP_OKAY;
}


/** adds a set covering or bound disjunction constraint to the original problem */
static
SCIP_RETCODE createNogood(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   bdlen,              /**< length of bound disjunction */
   SCIP_VAR**            bdvars,             /**< array of variables in bound disjunction */
   SCIP_BOUNDTYPE*       bdtypes,            /**< array of bound types in bound disjunction */
   SCIP_Real*            bdbounds,           /**< array of bounds in bound disjunction */
   SCIP_Bool             local,              /**< is constraint valid only locally? */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   SCIP_Bool*            success             /**< pointer to store whether the cutoff constraint was created successfully */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR** consvars;
   SCIP_Bool isbinary;
   char name[SCIP_MAXSTRLEN];
   int i;

   assert(scip != NULL);
   assert(bdlen >= 1);
   assert(bdvars != NULL);
   assert(bdtypes != NULL);
   assert(bdbounds != NULL);
   assert(success != NULL);

   /* initialize */
   *success = FALSE;
   cons = NULL;
   consvars = NULL;

   /* create constraint name */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "undercover_cutoff");

   /* check if all variables in the cover are binary */
   isbinary = TRUE;
   for( i = bdlen-1; i >= 0 && isbinary; i-- )
   {
      isbinary = isbinary && SCIPvarIsBinary(bdvars[i]);
   }

   /* if all variables in the cover are binary, then we create a logicor constraint */
   if( isbinary )
   {
      /* allocate memory for constraint variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, bdlen) );

      /* build up constraint */
      for( i = bdlen-1; i >= 0; i-- )
      {
         assert(bdtypes[i] == SCIP_BOUNDTYPE_LOWER || SCIPisFeasZero(scip, bdbounds[i]));
         assert(bdtypes[i] == SCIP_BOUNDTYPE_UPPER || SCIPisFeasEQ(scip, bdbounds[i], 1.0));

         if( bdtypes[i] == SCIP_BOUNDTYPE_LOWER )
         {
            consvars[i] = bdvars[i];
         }
         else
         {
            assert(bdtypes[i] == SCIP_BOUNDTYPE_UPPER);
            SCIP_CALL( SCIPgetNegatedVar(scip, bdvars[i], &consvars[i]) );
         }
      }

      /* create nogood constraint */
      SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, name, bdlen, consvars,
            FALSE, TRUE, FALSE, FALSE, TRUE, local, FALSE, dynamic, removable, FALSE) );
   }
   /* otherwise we create a bound disjunction constraint as given */
   else
   {
      /* create nogood constraint */
      SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &cons, name, bdlen, bdvars, bdtypes, bdbounds,
            FALSE, TRUE, FALSE, FALSE, TRUE, local, FALSE, dynamic, removable, FALSE) );
   }

   /* add and release constraint if created successfully */
   if( cons != NULL )
   {
      if( local )
      {
         SCIP_CALL( SCIPaddConsLocal(scip, cons, NULL) );
      }
      else
      {
         SCIP_CALL( SCIPaddCons(scip, cons) );
      }

      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      *success = TRUE;
   }

   /* free memory */
   SCIPfreeBufferArrayNull(scip, &consvars);

   return SCIP_OKAY;
}


/** solve covering problem */
static
SCIP_RETCODE solveCoveringProblem(
   SCIP*                 coveringscip,       /**< SCIP data structure for the covering problem */
   int                   ncoveringvars,      /**< number of the covering problem's variables */
   SCIP_VAR**            coveringvars,       /**< array of the covering problem's variables */
   int*                  coversize,          /**< size of the computed cover */
   int*                  cover,              /**< array to store indices of the variables in the computed cover
                                              *   (should be ready to hold ncoveringvars entries) */
   SCIP_Real             timelimit,          /**< time limit */
   SCIP_Real             memorylimit,        /**< memory limit */
   SCIP_Bool*            success             /**< feasible cover found? */
   )
{
   SCIP_Real* solvals;
   SCIP_Real totalpenalty;
   SCIP_RETCODE retcode;
   int i;

   assert(coveringscip != NULL);
   assert(coveringvars != NULL);
   assert(cover != NULL);
   assert(coversize != NULL);
   assert(timelimit > 0.0);
   assert(memorylimit > 0.0);
   assert(success != NULL);

   *success = FALSE;

   /* forbid call of heuristics and separators solving sub-CIPs */
   SCIP_CALL( SCIPsetSubscipsOff(coveringscip, TRUE) );

   /* set time and memory limit */
   SCIP_CALL( SCIPsetRealParam(coveringscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(coveringscip, "limits/memory", memorylimit) );

   /* do not abort on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(coveringscip, "misc/catchctrlc", FALSE) );
 
   /* disable output to console in optimized mode, enable in SCIP's debug mode */
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPsetIntParam(coveringscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(coveringscip, "display/freq", 100000) );
#else
   SCIP_CALL( SCIPsetIntParam(coveringscip, "display/verblevel", 0) );
#endif
 
   /* solve covering problem */
   retcode = SCIPsolve(coveringscip);
   
   /* Errors in solving the covering problem should not kill the overall solving process 
    * Hence, the return code is catched and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   { 
#ifndef NDEBUG
      SCIP_CALL( retcode );     
#endif
      SCIPwarningMessage("Error while solving covering problem in Undercover heuristic; subSCIP terminated with code <%d>\n",retcode);
      return SCIP_OKAY;
   }

   /* check, whether a feasible cover was found */
   if( SCIPgetNSols(coveringscip) == 0 )
      return SCIP_OKAY;

   /* store solution */
   SCIP_CALL( SCIPallocBufferArray(coveringscip, &solvals, ncoveringvars) );
   SCIP_CALL( SCIPgetSolVals(coveringscip, SCIPgetBestSol(coveringscip), ncoveringvars, coveringvars, solvals) );

   *coversize = 0;
   totalpenalty = 0.0;
   for( i = 0; i < ncoveringvars; i++ )
   {
      if( solvals[i] > 0.5 )
      {
         cover[*coversize] = i;
         (*coversize)++;
      }
      totalpenalty += SCIPvarGetObj(coveringvars[i]);
   }

   /* print solution if we are in SCIP's debug mode */
   assert(SCIPgetBestSol(coveringscip) != NULL);
   SCIPdebugMessage("found a feasible cover: %d/%d variables fixed, normalized penalty=%g\n\n",
      *coversize, SCIPgetNOrigVars(coveringscip), SCIPgetSolOrigObj(coveringscip, SCIPgetBestSol(coveringscip))/(totalpenalty+SCIPsumepsilon(coveringscip)));
   SCIPdebug( SCIP_CALL( SCIPprintSol(coveringscip, SCIPgetBestSol(coveringscip), NULL, FALSE) ) );
   SCIPdebugMessage("\r                                                  \n");

   *success = TRUE;

   /* free array of solution values */
   SCIPfreeBufferArray(coveringscip, &solvals);

   return SCIP_OKAY;
}


/** computes fixing order and returns whether order has really been changed */
static
SCIP_RETCODE computeFixingOrder(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int                   nvars,              /**< number of variables in the original problem */
   SCIP_VAR**            vars,               /**< variables in the original problem */
   int                   coversize,          /**< size of the cover */
   int*                  cover,              /**< problem indices of the variables in the cover */
   int                   lastfailed,         /**< position in cover array of the variable the fixing of which yielded
                                              *   infeasibility in last dive (or >= coversize, in which case *success
                                              *   is always TRUE) */
   SCIP_Bool*            success             /**< has order really been changed? */
   )
{
   SCIP_Real* scores;
   SCIP_Real bestscore;
   int i;

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(nvars >= 1);
   assert(vars != NULL);
   assert(coversize >= 1);
   assert(cover != NULL);
   assert(lastfailed >= 0);

   *success = FALSE;

   /* if fixing the first variable had failed, do not try with another order */
   if( lastfailed == 0 )
      return SCIP_OKAY;

   /* allocate buffer array for score values */
   SCIP_CALL( SCIPallocBufferArray(scip, &scores, coversize) );

   /* compute score values */
   bestscore = SCIPinfinity(scip);
   for( i = coversize-1; i >= 0; i-- )
   {
      SCIP_VAR* var;

      /* get variable in the cover */
      assert(cover[i] >= 0);
      assert(cover[i] < nvars);
      var = vars[cover[i]];

      /* compute score; switch sign because we will sort in non-decreasing order */
      scores[i] = -(heurdata->conflictweight * SCIPgetVarConflictScore(scip, var)
         + heurdata->inferenceweight * SCIPgetVarAvgInferenceCutoffScore(scip, var, heurdata->cutoffweight));
      assert(scores[i] <= 0.0);

      /* update maximum score */
      bestscore = MIN(bestscore, scores[i]);
   }
   bestscore -= 1.0;

   /* put integers to the front */
   if( heurdata->fixintfirst )
   {
      for( i = coversize-1; i >= 0; i-- )
      {
         if( SCIPvarIsIntegral(vars[cover[i]]) )
            scores[i] += bestscore;
      }
   }

   /* put last failed variable to the front */
   if( lastfailed < coversize )
   {
      scores[lastfailed] += 2*bestscore;
      i = cover[lastfailed];
   }

   /* sort by non-decreasing (negative) score */
   SCIPsortRealInt(scores, cover, coversize);
   assert(lastfailed >= coversize || cover[0] == i);

   /* free buffer memory */
   SCIPfreeBufferArray(scip, &scores);

   *success = TRUE;

   return SCIP_OKAY;
}


/** gets fixing value */
static
SCIP_RETCODE getFixingValue(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_VAR*             var,                /**< variable in the original problem */
   SCIP_Real*            val,                /**< buffer for returning fixing value */
   char                  fixalt,             /**< source of the fixing value: 'l'p relaxation, 'n'lp relaxation, 'i'ncumbent solution */
   SCIP_Bool*            success,            /**< could value be retrieved successfully? */
   int                   bdlen,              /**< current length of probing path */
   SCIP_VAR**            bdvars,             /**< array of variables with bound changes along probing path */
   SCIP_BOUNDTYPE*       bdtypes,            /**< array of bound types in bound disjunction */
   SCIP_Real*            oldbounds           /**< array of bounds before fixing */
   )
{
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(var != NULL);
   assert(val != NULL);
   assert(success != NULL);

   *success = FALSE;

   switch( fixalt )
   {
   case 'l':
      /* get the last LP solution value */
      *val = SCIPvarGetLPSol(var);
      *success = TRUE;
      break;
   case 'n':
      /* only call this function if nlp relaxation is available */
      assert(SCIPisNLPConstructed(scip));

      /* the solution values are already available */
      if( heurdata->nlpsolved )
      {
         assert(!heurdata->nlpfailed);

         /* retrieve nlp solution value */
         assert(SCIPgetNLP(scip) != NULL);
         SCIP_CALL( SCIPnlpGetVarSolVal(SCIPgetNLP(scip), var, val) );
         *success = TRUE;
      }
      /* solve nlp relaxation unless it has not failed too often before */
      else if( !heurdata->nlpfailed )
      {
         SCIP_NLPSOLSTAT stat;
         int i;

         /* restore bounds at start of probing, since otherwise, if in backtrack mode, nlp solver becomes most likely
          * locally infeasible */
         SCIP_CALL( SCIPstartDiveNLP(scip) );

         for( i = bdlen-1; i >= 0; i-- )
         {  /*lint --e{850}*/
            SCIP_VAR* relaxvar;
            SCIP_Real lb;
            SCIP_Real ub;

            relaxvar = bdvars[i];

            /* both bounds were tightened */
            if( i > 0 && bdvars[i-1] == relaxvar )
            {
               assert(bdtypes[i] != bdtypes[i-1]);

               lb = bdtypes[i] == SCIP_BOUNDTYPE_UPPER ? oldbounds[i] : oldbounds[i-1];
               ub = bdtypes[i] == SCIP_BOUNDTYPE_UPPER ? oldbounds[i-1] : oldbounds[i];
               i--;
            }
            /* lower bound was tightened */            
            else if( bdtypes[i] == SCIP_BOUNDTYPE_UPPER )
            {
               lb = oldbounds[i];
               ub = SCIPvarGetUbLocal(relaxvar);
            }
            /* upper bound was tightened */            
            else
            {
               lb = SCIPvarGetLbLocal(relaxvar);
               ub = oldbounds[i];
            }

            assert(SCIPisLE(scip, lb, SCIPvarGetLbLocal(relaxvar)));
            assert(SCIPisGE(scip, ub, SCIPvarGetUbLocal(relaxvar)));

            /* relax bounds */
            SCIP_CALL( SCIPchgVarBoundsDiveNLP(scip, relaxvar, lb, ub) );
         }

         /* activate nlp solver output if we are in SCIP's debug mode */
         SCIPdebug( SCIP_CALL( SCIPnlpSetIntPar(SCIPgetNLP(scip), SCIP_NLPPAR_VERBLEVEL, 1) ) );

         /* set starting point to lp solution */
         SCIP_CALL( SCIPsetNLPInitialGuessSol(scip, NULL) );

         /* solve nlp relaxation */
         SCIP_CALL( SCIPsolveDiveNLP(scip) );
         stat = SCIPgetNLPSolstat(scip);
         *success = stat == SCIP_NLPSOLSTAT_GLOBOPT || stat == SCIP_NLPSOLSTAT_LOCOPT || stat == SCIP_NLPSOLSTAT_FEASIBLE;

         SCIPdebugMessage("solving nlp relaxation to obtain fixing values %s (stat=%d)\n", *success ? "successful" : "failed", stat);

         if( *success )
         {
            /* solving succeeded */
            heurdata->nnlpfails = 0;
            heurdata->nlpsolved = TRUE;

            /* retrieve nlp solution value */
            assert(SCIPgetNLP(scip) != NULL);
            SCIP_CALL( SCIPnlpGetVarSolVal(SCIPgetNLP(scip), var, val) );
         }
         else
         {
            /* solving failed */
            heurdata->nnlpfails++;
            heurdata->nlpfailed = TRUE;
            heurdata->nlpsolved = FALSE;

            SCIPdebugMessage("solving nlp relaxation failed (%d time%s%s)\n",
               heurdata->nnlpfails, heurdata->nnlpfails > 1 ? "s" : "", heurdata->nnlpfails >= MAXNLPFAILS ? ", will not be called again" : "");
         }
      }
      break;
   case 'i':
      /* only call this function if a feasible solution is available */
      assert(SCIPgetBestSol(scip) != NULL);

      /* get value in the incumbent solution */
      *val = SCIPgetSolVal(scip, SCIPgetBestSol(scip), var);
      *success = TRUE;
      break;
   default:
      break;
   }

   return SCIP_OKAY;
}


/** calculates up to four alternative values for backtracking, if fixing the variable failed. 
 * The alternatives are the two bounds of the variable, and the averages of the bounds and the fixing value.
 * For infinite bounds, fixval +/- abs(fixval) will be used instead. 
 */
static
void calculateAlternatives(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_VAR*             var,                /**< variable to calculate alternatives for */
   SCIP_Real             fixval,             /**< reference fixing value */
   int*                  nalternatives,      /**< number of fixing values computed */
   SCIP_Real*            alternatives        /**< array to store the alternative fixing values */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   /* for binary variables, there is only two possible fixing values */
   if( SCIPvarIsBinary(var) )
   {
      if( SCIPisFeasEQ(scip, fixval, 0.0) || SCIPisFeasEQ(scip, fixval, 1.0) )
      {
         alternatives[0] = 1.0 - fixval;
         *nalternatives = 1;
      }
      else
      {
         alternatives[0] = 0.0;
         alternatives[1] = 1.0;
         *nalternatives = 2;
      }
      return;
   }

   /* get bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);

   /* if lower bound is infinite, use x'-|x'|; if x' is zero, use -1.0 instead */
   if( SCIPisInfinity(scip, -lb) )
      lb = SCIPisFeasZero(scip, fixval) ? -1.0 : fixval - ABS(fixval);

   /* if upper bound is infinite, use x'+|x'|; if x' is zero, use 1.0 instead */
   if( SCIPisInfinity(scip, ub) )
      ub = SCIPisFeasZero(scip, fixval) ? 1.0 : fixval + ABS(fixval);

   assert(!SCIPisEQ(scip, lb, ub));

   /* collect alternatives */
   *nalternatives = 0;

   /* use lower bound if not equal to x' */
   if( !SCIPisFeasEQ(scip, lb, fixval) )
   {
      alternatives[*nalternatives] = lb;
      (*nalternatives)++;
   }
   
   /* use upper bound if not equal to x' */
   if( !SCIPisFeasEQ(scip, ub, fixval) )
   {
      alternatives[*nalternatives] = ub;
      (*nalternatives)++;
   }
   
   /* use the average of x' and lower bound as alternative value, if this is not equal to any of the other values */
   if( !SCIPisFeasEQ(scip, lb, fixval) && (!SCIPvarIsIntegral(var) || !SCIPisFeasEQ(scip, lb, fixval-1)) )
   {
      alternatives[*nalternatives] = (lb+fixval)/2.0;
      (*nalternatives)++; 
   }

   /* use the average of x' and upper bound as alternative value, if this is not equal to any of the other values */
   if( !SCIPisFeasEQ(scip, ub, fixval) && (!SCIPvarIsIntegral(var) || !SCIPisFeasEQ(scip, ub, fixval+1)) )
   {
      alternatives[*nalternatives] = (ub+fixval)/2.0;
      (*nalternatives)++; 
   }

   assert(*nalternatives <= 4);

   return;
}


/** rounds the given fixing value */
static
SCIP_RETCODE roundFixingValue(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_Real*            val,                /**< fixing value to be rounded */
   SCIP_VAR*             var,                /**< corresponding variable */
   SCIP_Bool             locksrounding       /**< shall we round according to locks? (otherwise to nearest integer) */
   )
{
   SCIP_Real x;

   x = *val;

   /* if integral within feasibility tolerance, only shift to nearest integer */
   if( SCIPisFeasIntegral(scip, x) )
      x = SCIPfeasFrac(scip, x) < 0.5 ? SCIPfeasFloor(scip, x) : SCIPfeasCeil(scip, x);

   /* round in the direction of least locks with fractionality as tie breaker */
   else if( locksrounding )
   {
      if( SCIPvarGetNLocksDown(var) < SCIPvarGetNLocksUp(var) ) 
         x = SCIPfeasFloor(scip, x);
      else if( SCIPvarGetNLocksDown(var) > SCIPvarGetNLocksUp(var) )
         x = SCIPfeasCeil(scip, x);
      else 
         x = SCIPfeasFrac(scip, x) < 0.5 ? SCIPfeasFloor(scip, x) : SCIPfeasCeil(scip, x);
   }
   /* round in the direction of least fractionality with locks as tie breaker */
   else
   {
      if( SCIPfeasFrac(scip, x) < 0.5)
         x = SCIPfeasFloor(scip, x);
      else if( SCIPfeasFrac(scip, x) > 0.5 )
         x = SCIPfeasCeil(scip, x);
      else
         x = SCIPvarGetNLocksDown(var) < SCIPvarGetNLocksUp(var) ? SCIPfeasFloor(scip, x) : SCIPfeasCeil(scip, x);
   }

   /* return rounded fixing value */
   *val = x;

   return SCIP_OKAY;
}


/** copy the solution of the subproblem to newsol */
static
SCIP_RETCODE copySol(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_SOL*             subsol,             /**< solution of the subproblem */
   SCIP_SOL**            newsol              /**< solution to the original problem */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables */
   SCIP_Real* subsolvals;                    /* solution values of the subproblem */
   int        nvars;
        
   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);
   assert(newsol != NULL);
   assert(*newsol != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   /* subSCIP may have more variable than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */ 
   assert(nvars <= SCIPgetNOrigVars(subscip)); 
 
   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );
   SCIP_CALL( SCIPsetSolVals(scip, *newsol, nvars, vars, subsolvals) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}


/** solve subproblem and pass best feasible solution to original SCIP instance */
static
SCIP_RETCODE solveSubproblem(
   SCIP*                 scip,               /**< SCIP data structure of the original problem */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   int                   coversize,          /**< size of the cover */
   int*                  cover,              /**< problem indices of the variables in the cover */
   SCIP_Real*            fixingvals,         /**< fixing values for the variables in the cover */
   SCIP_Real             timelimit,          /**< time limit */
   SCIP_Real             memorylimit,        /**< memory limit */
   SCIP_Longint          nodelimit,          /**< node limit */
   SCIP_Longint          nstallnodes,        /**< number of stalling nodes for the subproblem */
   SCIP_Bool*            validsolved,        /**< was problem constructed from a valid copy and solved (proven optimal or infeasible)? */
   SCIP_SOL**            sol,                /**< best solution found in subproblem (if feasible); *sol must be NULL, solution will be created */
   SCIP_Longint*         nusednodes          /**< number of nodes used for solving the subproblem */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP* subscip;
   SCIP_VAR** subvars;
   SCIP_VAR** vars;
   SCIP_HASHMAP* varmap;

   SCIP_RETCODE retcode;

   int nvars;
   int i;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(cover != NULL);
   assert(fixingvals != NULL);
   assert(coversize >= 1);
   assert(timelimit > 0.0);
   assert(memorylimit > 0.0);
   assert(nodelimit >= 1);
   assert(nstallnodes >= 1);
   assert(validsolved != NULL);
   assert(sol != NULL);
   assert(*sol == NULL);
   assert(nusednodes != NULL);

   *validsolved = FALSE;
   *nusednodes = 0;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* create subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(subscip), SCIPcalcHashtableSize(5 * nvars)) );

   /* copy original problem to subproblem; do not copy pricers */
   SCIP_CALL( SCIPcopy(scip, subscip, varmap, NULL, "undercoversub", heurdata->globalbounds, FALSE, validsolved) );

   if( heurdata->copycuts )
   {
      /** copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
      SCIP_CALL( SCIPcopyCuts(scip, subscip, varmap, NULL, heurdata->globalbounds) );
   }

   SCIPdebugMessage("problem copied, copy %svalid\n", *validsolved ? "" : "in");

   /* store subproblem variables */
   for( i = nvars-1; i >= 0; i-- )
   {
      subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, vars[i]);
      assert(subvars[i] != NULL);
   }

   /* fix subproblem variables in the cover */
   SCIPdebugMessage("fixing variables\n");
   for( i = coversize-1; i >= 0; i-- )
   {
      assert(cover[i] >= 0);
      assert(cover[i] < nvars);

      SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[cover[i]], fixingvals[i]) );
      SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[cover[i]], fixingvals[i]) );
   }

   /* set the parameters such that good solutions are found fast */
   SCIPdebugMessage("setting subproblem parameters\n");
   SCIP_CALL( SCIPsetEmphasis(subscip, SCIP_PARAMSETTING_FEASIBILITY, TRUE) );
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );
   SCIP_CALL( SCIPsetHeuristics(subscip, SCIP_PARAMSETTING_AGGRESSIVE, TRUE) );

   /* deactivate expensive pre-root heuristics, since it may happen that the lp relaxation of the subproblem is already
      infeasible; in this case, we do not want to waste time on heuristics before solving the root lp */
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/shiftandpropagate/freq", -1) );

   /* forbid recursive call of undercover heuristic */
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/"HEUR_NAME"/freq", -1) );
   
   /* set time, memory and node limits */
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nodelimit) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", nstallnodes) );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );
 
   /* disable output to console in optimized mode, enable in SCIP's debug mode */
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000) );
#else
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
#endif

   /* if there is already a solution, add an objective cutoff; note: this does not affect the validity of the subproblem
    * if we find solutions later, thus we do not set *validsolved to FALSE */
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_Real cutoff;
      SCIP_Real upperbound;

      assert(!SCIPisInfinity(scip, SCIPgetUpperbound(scip)));
      upperbound = SCIPgetUpperbound(scip);

      if( SCIPisInfinity(scip, -SCIPgetLowerbound(scip)) )
         cutoff = (upperbound >= 0 ? 1.0 - heurdata->minimprove : 1.0 + heurdata->minimprove) * upperbound;
      else
         cutoff = (1.0 - heurdata->minimprove) * upperbound + heurdata->minimprove * SCIPgetLowerbound(scip);

      cutoff = MIN(upperbound, cutoff);
      SCIP_CALL( SCIPsetObjlimit(subscip, cutoff) );

      SCIPdebugMessage("adding objective cutoff=%g (minimprove=%g)\n", cutoff, heurdata->minimprove);
   }
      
   /* solve subproblem */
   SCIPdebugMessage("solving subproblem started\n");
   retcode = SCIPsolve(subscip);
   
   /* Errors in solving the subproblem should not kill the overall solving process 
    * Hence, the return code is catched and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   { 
#ifndef NDEBUG
      SCIP_CALL( retcode );     
#endif
      SCIPwarningMessage("Error while solving subproblem in Undercover heuristic; subSCIP terminated with code <%d>\n",retcode);
      return SCIP_OKAY;
   }

   /* print solving statistics of subproblem if we are in SCIP's debug mode */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );

   /* store solving status; note: if we proved infeasibility in presence of an objective cutoff beyond the primal bound,
    * the subproblem was not a valid copy */
   *validsolved = *validsolved && (SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL
      || (SCIPgetStatus(subscip) == SCIP_STATUS_INFEASIBLE && (SCIPgetNSols(scip) == 0 || heurdata->minimprove <= 0.0)));
   *nusednodes = SCIPgetNNodes(subscip);

   /* if a solution was found for the subproblem, create corresponding solution in the original problem */
   if( SCIPgetNSols(subscip) > 0 && (SCIPgetStatus(subscip) != SCIP_STATUS_INFEASIBLE || heurdata->minimprove > 0.0) )
   {
      SCIP_SOL** subsols;
      SCIP_Bool success;
      int nsubsols;

      /* create solution */
      SCIP_CALL( SCIPcreateSol(scip, sol, heur) );

      /* check, whether a solution was found;
       * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted */
      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      assert(subsols != NULL);

      success = FALSE;
      for( i = 0; i < nsubsols && !success; i++ )
      {
         /* transform solution to original problem */
         SCIP_CALL( copySol(scip, subscip, subvars, subsols[i], sol) );

         /* try to add new solution to scip */
         SCIP_CALL( SCIPtrySol(scip, *sol, FALSE, TRUE, TRUE, TRUE, &success) );
      }

      if( success )
      {   
         assert(i >= 1);
         SCIPdebugMessage("heuristic found %d solutions in subproblem; solution %d feasible in original problem\n", nsubsols, i);
      }
      else
      {
         /* free solution structure, since we found no feasible solution */
         SCIP_CALL( SCIPfreeSol(scip, sol) );
         *sol = NULL;
      }

      /* if the best subproblem solution was not accepted in the original problem, we do not trust the solving status */
      *validsolved = *validsolved && i == 1;
   }

   /* free variable mapping hash map, array of subproblem variables, and subproblem */
   SCIPhashmapFree(&varmap);
   SCIPfreeBufferArray(scip, &subvars);
   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}


/** perform fixing of a variable and record bound disjunction information */
static
SCIP_RETCODE performFixing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to fix */
   SCIP_Real             val,                /**< fixing value */
   SCIP_Bool*            infeas,             /**< pointer to store whether the fixing lead to infeasibility */
   int*                  bdlen,              /**< current length of bound disjunction */
   SCIP_VAR**            bdvars,             /**< array of variables in bound disjunction */
   SCIP_BOUNDTYPE*       bdtypes,            /**< array of bound types in bound disjunction */
   SCIP_Real*            bdbounds,           /**< array of bounds in bound disjunction */
   SCIP_Real*            oldbounds           /**< array of bounds before fixing */
   )
{
   SCIP_Longint ndomredsfound;
   SCIP_Real oldlb;
   SCIP_Real oldub;
   int oldbdlen;

   assert(scip != NULL);
   assert(var != NULL);
   assert(val >= SCIPvarGetLbLocal(var));
   assert(val <= SCIPvarGetUbLocal(var));
   assert(infeas != NULL);
   assert(bdlen != NULL);
   assert(*bdlen >= 0);
   assert(*bdlen < 2*SCIPgetNVars(scip)-1);
   assert(bdvars != NULL);
   assert(bdtypes != NULL);
   assert(bdbounds != NULL);

   assert(!SCIPvarIsIntegral(var) || SCIPisFeasIntegral(scip, val));

   /* remember length of probing path */
   oldbdlen = *bdlen;

   /* get bounds of the variable to fix */
   oldlb = SCIPvarGetLbLocal(var);
   oldub = SCIPvarGetUbLocal(var);

   /* decrease upper bound to fixing value */
   *infeas = FALSE;
   if( SCIPisUbBetter(scip, val, oldlb, oldub) )
   {
      /* create next probing node */
      SCIP_CALL( SCIPnewProbingNode(scip) );
      SCIP_CALL( SCIPchgVarUbProbing(scip, var, val) );

      SCIPdebugMessage("tentatively decreasing upper bound of variable <%s> to %g for probing\n",
         SCIPvarGetName(var), val);

      /* store bound disjunction information */
      bdvars[*bdlen] = var;
      bdtypes[*bdlen] = SCIP_BOUNDTYPE_LOWER;
      bdbounds[*bdlen] = SCIPvarIsIntegral(var) ? SCIPfeasCeil(scip, val)+1.0 : val;
      oldbounds[*bdlen] = oldub;
      (*bdlen)++;

      /* propagate the bound change; conflict analysis is performed automatically */
      SCIP_CALL( SCIPpropagateProbing(scip, 0, infeas, &ndomredsfound) );
      SCIPdebugMessage("  --> propagation reduced %lld further domains\n", ndomredsfound);

      /* if propagation led to a cutoff, we backtrack immediately */
      if( *infeas )
      {
         *bdlen = oldbdlen;
         return SCIP_OKAY;
      }

      /* store bound before propagation */
      oldbounds[*bdlen] = oldlb;

      /* move fixing value into the new domain, since it may be outside due to numerical issues or previous propagation */
      oldlb = SCIPvarGetLbLocal(var);
      oldub = SCIPvarGetUbLocal(var);
      val = MIN(val, oldub);
      val = MAX(val, oldlb);

      assert(!SCIPvarIsIntegral(var) || SCIPisFeasIntegral(scip, val));
   }

   /* update lower bound to fixing value */
   *infeas = FALSE;
   if( SCIPisLbBetter(scip, val, oldlb, oldub) )
   {
      /* create next probing node */
      SCIP_CALL( SCIPnewProbingNode(scip) );
      SCIP_CALL( SCIPchgVarLbProbing(scip, var, val) );

      SCIPdebugMessage("tentatively increasing lower bound of variable <%s> to %g for probing\n",
         SCIPvarGetName(var), val);

      /* store bound disjunction information */
      bdvars[*bdlen] = var;
      bdtypes[*bdlen] = SCIP_BOUNDTYPE_UPPER;
      bdbounds[*bdlen] = SCIPvarIsIntegral(var) ? SCIPfeasCeil(scip, val)-1.0 : val;
      (*bdlen)++;

      /* propagate the bound change */
      SCIP_CALL( SCIPpropagateProbing(scip, 0, infeas, &ndomredsfound) );
      SCIPdebugMessage("  --> propagation reduced %lld further domains\n", ndomredsfound);

      /* if propagation led to a cutoff, we backtrack immediately */
      if( *infeas )
      {
         *bdlen = oldbdlen;
         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
}


/** main procedure of the undercover heuristic */
static
SCIP_RETCODE SCIPapplyUndercover(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_RESULT*          result,             /**< result data structure */
   SCIP_Real             timelimit,          /**< time limit */
   SCIP_Real             memorylimit,        /**< memory limit */
   SCIP_Longint          nstallnodes         /**< number of stalling nodes for the subproblem */
   )
{
   SCIP_HEURDATA* heurdata;                  /* heuristic data */
   SCIP_VAR** vars;                          /* original problem's variables */
   SCIP_CLOCK* clock;                        /* clock for updating time limit */

   SCIP* coveringscip;                       /* SCIP data structure for covering problem */
   SCIP_VAR** coveringvars;                  /* covering variables */
   SCIP_Real* fixingvals;                    /* array for storing fixing values used */
   int* cover;                               /* array to store problem indices of variables in the computed cover */

   SCIP_VAR** bdvars;                        /* array of variables in bound disjunction along the probing path */
   SCIP_BOUNDTYPE* bdtypes;                  /* array of bound types in bound disjunction along the probing path */
   SCIP_Real* bdbounds;                      /* array of bounds in bound disjunction along the probing path */
   SCIP_Real* oldbounds;                     /* array of bounds before fixing along the probing path */

   SCIP_Bool success;
   int bdlen;                                /* current length of bound disjunction along the probing path */
   int coversize;
   int nvars;
   int ncovers;
   int nunfixeds;
   int i;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(result != NULL);
   assert(*result == SCIP_DIDNOTFIND);

   /* create and start timing */
   SCIP_CALL( SCIPcreateClock(scip, &clock) );
   SCIP_CALL( SCIPstartClock(scip, clock) );

   /* initialize */
   fixingvals = NULL;
   cover = NULL;
   bdvars = NULL;
   bdtypes = NULL;
   bdbounds = NULL;
   oldbounds = NULL;
   bdlen = 0;
   coversize = 0;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* nlp relaxation has not been solved yet (only solve once, not again for each cover or dive, because it is expensive) */
   heurdata->nlpsolved = FALSE;

   /* if solving the nlp relaxation has failed too often in previous runs, or nlp and nlp solver is not available, we do
    * not even try */
   heurdata->nlpfailed = heurdata->nnlpfails >= MAXNLPFAILS || !SCIPisNLPConstructed(scip) || SCIPgetNNlpis(scip) == 0;

   /* get variable data of original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* create covering problem */
   success = FALSE;
   SCIP_CALL( SCIPcreate(&coveringscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(coveringscip) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coveringvars, nvars) );
   SCIP_CALL( createCoveringProblem(scip, coveringscip, coveringvars, heurdata, &success) );

   if( !success )
   {
      SCIPdebugMessage("creating covering problem failed, terminating\n");
      goto TERMINATE;
   }
   else
      SCIPdebugMessage("covering problem created successfully\n");

   /* count number of unfixed covering variables */
   nunfixeds = 0;
   for( i = nvars-1; i >= 0; i-- )
   {
      if( SCIPisFeasEQ(coveringscip, SCIPvarGetLbGlobal(coveringvars[i]), 1.0) )
         nunfixeds++;
   }

   /* update time limit */
   SCIP_CALL( updateTimelimit(scip, clock, &timelimit) );

   if( timelimit <= MINTIMELEFT )
   {
      SCIPdebugMessage("time limit hit, terminating\n");
      goto TERMINATE;
   }

   /* update memory left */
   memorylimit -= SCIPgetMemUsed(coveringscip)/1048576.0;

   /* allocate memory for storing bound disjunction information along probing path */
   SCIP_CALL( SCIPallocBufferArray(scip, &bdvars, 2*nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bdtypes, 2*nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bdbounds, 2*nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oldbounds, 2*nvars) );

   /* re-covering loop */
   SCIP_CALL( SCIPallocBufferArray(scip, &cover, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixingvals, nvars) );
   ncovers = 0;
   success = FALSE;
   while( ncovers <= heurdata->maxrecovers && !success )
   {
      int lastfailed;
      int ndives;
      int nfixingalts;
      int nfixedints;
      int nfixedconts;

      /* solve covering problem; free transformed covering problem immediately */
      SCIPdebugMessage("solving covering problem\n\n");
      success = FALSE;
      SCIP_CALL( solveCoveringProblem(coveringscip, nvars, coveringvars, &coversize, cover,
            timelimit, memorylimit + SCIPgetMemUsed(coveringscip)/1048576.0, &success) );
      assert(coversize >= 0);
      assert(coversize <= nvars);

      SCIP_CALL( SCIPfreeTransform(coveringscip) );
      ncovers++;

      /* terminate if no feasible cover was found */
      if( !success )
      {
         SCIPdebugMessage("no feasible cover found in covering problem %d, terminating\n", ncovers);
         goto TERMINATE;
      }

      /* terminate, if cover is empty or too large */
      if( coversize == 0 || coversize > nvars*(heurdata->maxcoversize) )
      {
         SCIPdebugMessage("terminating due to coversize=%d\n", coversize);
         goto TERMINATE;
      }

      /* round-fix-propagate-analyze-backtrack-reorder */
      nfixingalts = (int) strlen(heurdata->fixingalts);
      assert(nfixingalts >= 1);

      /* reordering loop */
      ndives = 0;
      nfixedints = 0;
      nfixedconts = 0;
      success = FALSE;
      lastfailed = coversize;
      while( ndives <= heurdata->maxreorders && !success )
      {
         SCIP_Bool infeas;
         SCIP_Bool lpsolved;
         SCIP_Bool reordered;

         /* compute fixing order */
         SCIP_CALL( computeFixingOrder(scip, heurdata, nvars, vars, coversize, cover, lastfailed, &reordered) );
         reordered = reordered || ndives == 0;
         SCIPdebugMessage("%sordering variables in cover %s\n", ndives == 0 ? "" : "re", reordered ? "" : "failed");

         /* stop if order has not changed */
         if( !reordered )
            break;

         /* start probing in original problem */
         lpsolved = SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL;
         SCIP_CALL( SCIPstartProbing(scip) );
         ndives++;

         /* round-fix-propagate-analyze-backtrack for each variable in the cover */
         nfixedints = 0;
         nfixedconts = 0;
         bdlen = 0;
         infeas = FALSE;
         for( i = 0; i < coversize && !infeas; i++ )
         {
            SCIP_Real* boundalts;
            SCIP_Real* usedvals;
            SCIP_Real val;
            int nbacktracks;
            int nboundalts;
            int nfailedvals;
            int nusedvals;
            int probingdepth;
            int idx;

            /* get probindex of next variable in the cover */
            idx = cover[i];

            /* nothing to do if the variable was already fixed, e.g., by propagation */
            if( SCIPisEQ(scip, SCIPvarGetLbLocal(vars[idx]), SCIPvarGetUbLocal(vars[idx])) )
            {
               fixingvals[i] = SCIPvarGetLbLocal(vars[idx]);
               continue;
            }

            /* we will store the fixing values already used to avoid try the same value twice */
            SCIP_CALL( SCIPallocBufferArray(scip, &usedvals, heurdata->maxbacktracks+1) );
            nusedvals = 0;

            /* backtracking loop */
            infeas = TRUE;
            nfailedvals = 0;
            nboundalts = 0;
            boundalts = NULL;
            val = 0.0;
            for( nbacktracks = 0; nbacktracks <= heurdata->maxbacktracks+nfailedvals && infeas; nbacktracks++ )
            {
               SCIP_Real oldlb;
               SCIP_Real oldub;
               SCIP_Bool usedbefore;
               int j;

               probingdepth = SCIPgetProbingDepth(scip);

               /* get fixing value */
               if( nbacktracks < nfixingalts )
               {
                  /* if the lp relaxation is not solved, we do not even try to retrieve the lp solution value;
                   * if the nlp relaxation is not constructed, we do not even try to retrieve the nlp solution value; 
                   * if there is no feasible solution yet, we do not even try to obtain the value in the incumbent */
                  success = FALSE;
                  if( (heurdata->fixingalts[nbacktracks] != 'l' || lpsolved)
                     && (heurdata->fixingalts[nbacktracks] != 'n' || !heurdata->nlpfailed)
                     && (heurdata->fixingalts[nbacktracks] != 'i' || SCIPgetBestSol(scip) != NULL) )
                  {
                     SCIP_CALL( getFixingValue(scip, heurdata, vars[idx], &val, heurdata->fixingalts[nbacktracks], &success, bdlen, bdvars, bdtypes, oldbounds) );
                  }

                  if( !success )
                  {
                     SCIPdebugMessage("retrieving fixing value '%c' for variable <%s> failed, trying next in the list\n",
                        heurdata->fixingalts[nbacktracks], SCIPvarGetName(vars[idx]));
                     nfailedvals++;
                     continue;
                  }

                  /* for the first (successfully retrieved) fixing value, compute (at most 4) bound dependent
                   * alternative fixing values */
                  if( boundalts == NULL )
                  {
                     SCIP_CALL( SCIPallocBufferArray(scip, &boundalts, 4) );
                     nboundalts = 0;
                     calculateAlternatives(scip, vars[idx], val, &nboundalts, boundalts);
                     assert(nboundalts >= 0);
                     assert(nboundalts <= 4);
                  }
               }
               /* get alternative fixing value */
               else if( boundalts != NULL && nbacktracks < nfixingalts+nboundalts )
               {
                  assert(nbacktracks-nfixingalts >= 0);
                  val = boundalts[nbacktracks-nfixingalts];
               }
               else
                  break;

               /* round fixing value */
               if( SCIPvarIsIntegral(vars[idx]) && !SCIPisIntegral(scip, val) )
               {
                  SCIP_CALL( roundFixingValue(scip, &val, vars[idx], heurdata->locksrounding) );
                  assert(SCIPisIntegral(scip, val));
               }

               /* move value into the domain, since it may be outside due to numerical issues or previous propagation */
               oldlb = SCIPvarGetLbLocal(vars[idx]);
               oldub = SCIPvarGetUbLocal(vars[idx]);
               val = MIN(val, oldub);
               val = MAX(val, oldlb);

               assert(!SCIPvarIsIntegral(vars[idx]) || SCIPisFeasIntegral(scip, val));

               /* check if this fixing value was already used */
               usedbefore = FALSE;
               for( j = nusedvals-1; j >= 0 && !usedbefore; j-- )
                  usedbefore = SCIPisFeasEQ(scip, val, usedvals[j]);

               if( usedbefore )
               {
                  nfailedvals++;
                  continue;
               }

               /* store fixing value */
               assert(nusedvals < heurdata->maxbacktracks);
               usedvals[nusedvals] = val;
               nusedvals++;

               /* fix-propagate-analyze */
               SCIP_CALL( performFixing(scip, vars[idx], val, &infeas, &bdlen, bdvars, bdtypes, bdbounds, oldbounds) );

               /* if infeasible, backtrack and try alternative fixing value */
               if( infeas )
               {
                  SCIPdebugMessage("  --> cutoff detected - backtracking\n");
                  SCIP_CALL( SCIPbacktrackProbing(scip, probingdepth) );
               }
            }

            /* free array of alternative backtracking values */
            if( boundalts != NULL)
               SCIPfreeBufferArray(scip, &boundalts);
            SCIPfreeBufferArray(scip, &usedvals);

            /* backtracking loop unsuccessful */
            if( infeas )
            {
               SCIPdebugMessage("no feasible fixing value found for variable <%s> in fixing order %d\n",
                  SCIPvarGetName(vars[idx]), ndives);
               break;
            }
            /* fixing successful */
            else
            {
               /* store successful fixing value */
               fixingvals[i] = val;

               /* statistics */
               if( SCIPvarGetType(vars[idx]) == SCIP_VARTYPE_CONTINUOUS )
                  nfixedconts++;
               else
                  nfixedints++;
            }
         }

         /* end of dive */
         SCIP_CALL( SCIPendProbing(scip) );
         success = !infeas;
         lastfailed = i;
         assert(infeas || i == coversize);
         assert(!infeas || i < coversize);
      }

      /* update time limit */
      SCIPdebugMessage("%d dive%s of fix-and-propagate for cover %d took %.1f seconds\n", ndives, ndives > 1 ? "s" : "", ncovers, SCIPgetClockTime(scip, clock));
      SCIP_CALL( updateTimelimit(scip, clock, &timelimit) );

      if( timelimit <= MINTIMELEFT )
      {
         SCIPdebugMessage("time limit hit, terminating\n");
         goto TERMINATE;
      }

      /* no feasible fixing could be found for the current cover */
      if( !success )
      {
         SCIPdebugMessage("no feasible fixing values found for cover %d\n", ncovers);
      }
      else
      {
         SCIP_SOL* sol;
         SCIP_Longint nsubnodes;
         SCIP_Bool validsolved;

         SCIPdebugMessage("heuristic successfully fixed %d variables (%d integral, %d continuous) during probing\n",
            nfixedints+nfixedconts, nfixedints, nfixedconts); /*lint !e771*/

         /* solve sub-CIP and pass feasible solutions to original problem */
         success = FALSE;
         validsolved = FALSE;
         sol = NULL;
         nsubnodes = 0;

         SCIP_CALL( solveSubproblem(scip, heur, coversize, cover, fixingvals,
               timelimit, memorylimit, heurdata->maxnodes, nstallnodes, &validsolved, &sol, &nsubnodes) );

         /* update number of sub-CIP nodes used by heuristic so far */
         heurdata->nusednodes += nsubnodes;

         /* if the subproblem was constructed from a valid copy and solved, try to forbid the assignment of fixing
          * values to variables in the cover */
         if( validsolved )
         {
            SCIP_Real maxvarsfac;
            SCIP_Bool useconf;

            if( SCIPgetRealParam(scip, "conflict/maxvarsfac", &maxvarsfac) != SCIP_OKAY )
               maxvarsfac = 0.1;

            useconf = bdlen > 0 && bdlen < maxvarsfac*nvars;
            if( useconf )
            {
               /* even if we had reset the global bounds at start of probing, the constraint might be only locally valid due to local constraints/cuts */ 
               SCIP_CALL( createNogood(scip, bdlen, bdvars, bdtypes, bdbounds, SCIPgetDepth(scip) > 0, TRUE, TRUE, &success) );
            }

            SCIPdebugMessage("subproblem solved (%s), forbidding assignment in original problem %s, %sconflict length=%d\n",
               sol == NULL ? "infeasible" : "optimal",
               success ? "successful" : "failed", useconf ? "" : "skipped due to ", bdlen);
         }

         /* heuristic succeeded */
         success = sol != NULL;
         if( success )
         {
            *result = SCIP_FOUNDSOL;
            success = TRUE;

            /* update time limit */
            SCIP_CALL( updateTimelimit(scip, clock, &timelimit) );

            /* call nlp local search heuristic unless it has failed too often */
            if( heurdata->postnlp && heurdata->npostnlpfails < MAXPOSTNLPFAILS )
            {
               if( nfixedconts == 0 && validsolved )
               {
                  SCIPdebugMessage("subproblem solved to optimality while all covering variables are integral, hence skipping nlp local search\n");
               }
               else if( timelimit <= MINTIMELEFT )
               {
                  SCIPdebugMessage("time limit hit, skipping nlp local search\n");
               }
               else if( heurdata->nlpheur == NULL )
               {
                  SCIPdebugMessage("nlp heuristic not found, skipping nlp local search\n");
               }
               else
               {
                  SCIP_RESULT nlpresult;

                  SCIP_CALL( SCIPapplyHeurSubNlp(scip, heurdata->nlpheur, &nlpresult, sol, -1LL, timelimit, heurdata->minimprove, NULL) );
                  SCIPdebugMessage("nlp local search %s\n", nlpresult == SCIP_FOUNDSOL ? "successful" : "failed");

                  if( nlpresult == SCIP_FOUNDSOL )
                     heurdata->npostnlpfails = 0;
                  else
                     heurdata->npostnlpfails++;
               }
            }

            /* free solution */
            SCIP_CALL( SCIPfreeSol(scip, &sol) );
         }
      }

      /* heuristic failed but we have another re-covering try, hence we forbid the current cover in the covering problem */
      if( !success && ncovers <= heurdata->maxrecovers )
      {
         SCIP_Bool infeas;
         int diversification;

         /* compute minimal number of unfixed covering variables (in the cover) which have to change their value */
         diversification = (int) SCIPfeasCeil(scip, (heurdata->recoverdiv) * (SCIP_Real) nunfixeds);
         diversification = MAX(diversification, 1);

         /* forbid unsuccessful cover globally in covering problem */
         SCIP_CALL( forbidCover(coveringscip, nvars, coveringvars, coversize, cover, diversification, &success, &infeas) );

         if( infeas )
         {
            SCIPdebugMessage("re-covering problem infeasible (diversification=%d), terminating\n", diversification);
            goto TERMINATE;
         }
         else if( !success )
         {
            SCIPdebugMessage("failed to forbid current cover in the covering problem, terminating\n");
            goto TERMINATE;
         }
         else
         {
            SCIPdebugMessage("added constraint to the covering problem in order to forbid current cover\n");
            success = FALSE;
         }
      }
   }

 TERMINATE:
   if( *result != SCIP_FOUNDSOL && *result != SCIP_DELAYED )
   {
      SCIPdebugMessage("heuristic terminating unsuccessfully\n");
   }

   /* we must remain in nlp diving mode until here to be able to retrieve nlp solution values easily */
   assert((SCIPisNLPConstructed(scip) == FALSE && heurdata->nlpsolved == FALSE) ||
      (SCIPisNLPConstructed(scip) == TRUE && heurdata->nlpsolved == SCIPnlpIsDiving(SCIPgetNLP(scip))));
   if( heurdata->nlpsolved )
   {
      SCIP_CALL( SCIPendDiveNLP(scip) );
   }

   /* free arrays for storing the cover */
   SCIPfreeBufferArrayNull(scip, &fixingvals);
   SCIPfreeBufferArrayNull(scip, &cover);

   /* free arrays for storing bound disjunction information along probing path */
   SCIPfreeBufferArrayNull(scip, &oldbounds);
   SCIPfreeBufferArrayNull(scip, &bdbounds);
   SCIPfreeBufferArrayNull(scip, &bdtypes);
   SCIPfreeBufferArrayNull(scip, &bdvars);

   /* free covering problem */
   SCIPfreeBufferArray(scip, &coveringvars);
   SCIP_CALL( SCIPfree(&coveringscip) );

   /* free clock */
   SCIP_CALL( SCIPfreeClock(scip, &clock) );

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyUndercover)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurUndercover(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeUndercover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);
   assert(heur != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitUndercover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* we use the nlp for building up the covering problem, i.e., even if we do not solve the nlp relaxation or perform
    * nlp local search;
    * however, if we want to use nlp fixing values exclusively and no nlp solver is available,
    * heuristic will not run anyway */
   if( strcmp(heurdata->fixingalts, "n") != 0 || SCIPgetNNlpis(scip) > 0 )
      SCIPmarkRequireNLP(scip);

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#define heurExitUndercover NULL


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolUndercover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   int h;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* initialize counters to zero */
   heurdata->nusednodes = 0;
   heurdata->npostnlpfails = 0;
   heurdata->nnlpfails = 0;

   /* if the heuristic is called at the root node, we may want to be called directly after the initial root LP solve */
   if( heurdata->beforecuts && SCIPheurGetFreqofs(heur) == 0 )
      SCIPheurSetTimingmask(heur, SCIP_HEURTIMING_DURINGLPLOOP);

   /* find nonlinear constraint handlers */
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->nlconshdlrs, 3) );
   h = 0;
   heurdata->nlconshdlrs[h] = SCIPfindConshdlr(scip, "and");
   if( heurdata->nlconshdlrs[h] != NULL )
      h++;
   heurdata->nlconshdlrs[h] = SCIPfindConshdlr(scip, "quadratic");
   if( heurdata->nlconshdlrs[h] != NULL )
      h++;
   heurdata->nlconshdlrs[h] = SCIPfindConshdlr(scip, "soc");
   if( heurdata->nlconshdlrs[h] != NULL )
      h++;
   heurdata->nnlconshdlrs = h;

   /* find nlp local search heuristic */
   heurdata->nlpheur = SCIPfindHeur(scip, "subnlp");

   /* add global linear constraints to nlp relaxation */
   if( SCIPisNLPConstructed(scip) && heurdata->nlpheur != NULL )
   {
      SCIP_CALL( SCIPaddLinearConsToNlpHeurSubNlp(scip, heurdata->nlpheur, TRUE, TRUE) );
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolUndercover)
{
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free array of nonlinear constraint handlers */
   SCIPfreeMemoryArray(scip, &heurdata->nlconshdlrs);

   /* reset timing, if it was changed temporary (at the root node) */
   SCIPheurSetTimingmask(heur, HEUR_TIMING);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecUndercover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;                  /* heuristic data */
   SCIP_Real timelimit;                      /* time limit for the subproblem */
   SCIP_Real memorylimit;                    /* memory limit for the subproblem */
   SCIP_Longint nstallnodes;                 /* number of stalling nodes for the subproblem */
   SCIP_Bool run;

   int h;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* only call heuristic once at the root */
   if( SCIPgetDepth(scip) == 0 && SCIPheurGetNCalls(heur) > 0 )
      return SCIP_OKAY;
   
   /* if we want to use nlp fixing values exclusively and no nlp solver is available, we cannot run */
   if( strcmp(heurdata->fixingalts, "n") == 0 && SCIPgetNNlpis(scip) == 0 )
   {
      SCIPdebugMessage("skipping undercover heuristic: want to use nlp fixing values exclusively, but no nlp solver available\n");
      return SCIP_OKAY;
   }

   /* calculate stallnode limit */
   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));
   
   /* reward heuristic if it succeeded often */
   nstallnodes = (SCIP_Longint)(nstallnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur) + 1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nstallnodes -= SUBMIPSETUPCOSTS * SCIPheurGetNCalls(heur);  /* account for the setup costs of the sub-CIP */
   nstallnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nstallnodes -= heurdata->nusednodes;
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);
   nstallnodes = MAX(nstallnodes, 1);

   /* only call heuristics if we have enough nodes left to call sub-CIP solving */
   if( nstallnodes < heurdata->minnodes )
   {
      SCIPdebugMessage("skipping undercover heuristic: nstallnodes=%"SCIP_LONGINT_FORMAT", minnodes=%"SCIP_LONGINT_FORMAT"\n", nstallnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

   /* only call heuristics if we have enough time left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   if( timelimit <= 2*MINTIMELEFT )
   {
      SCIPdebugMessage("skipping undercover heuristic: time left=%g\n", timelimit);
      return SCIP_OKAY;
   }

   /* only call heuristics if we have enough memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(scip, memorylimit) )
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
   if( memorylimit <= 0.0 )
   {
      SCIPdebugMessage("skipping undercover heuristic: too little memory\n");
      return SCIP_OKAY;
   }

   /* only call heuristic if nonlinear constraints are present */
   run = FALSE;
   for( h = heurdata->nnlconshdlrs-1; h >= 0 && !run; h-- )
   {
      run = SCIPconshdlrGetNConss(heurdata->nlconshdlrs[h]) > 0;
   }

   /* go through all nlrows and check for general nonlinearities */
   if( SCIPisNLPConstructed(scip) )
   {
      SCIP_NLP* nlp;
      SCIP_NLROW** nlrows;
      int nnlrows;
      int i;

      /* get nlp */
      nlp = SCIPgetNLP(scip);
      assert(nlp != NULL);

      /* get nlrows */
      nnlrows = SCIPnlpGetNNlRows(nlp);
      nlrows = SCIPnlpGetNlRows(nlp);

      /* check for an nlrow with nontrivial expression tree or quadratic terms; start from 0 since we expect the linear
       * nlrows at the end */
      for( i = nnlrows-1; i >= 0 && !run; i-- )
      {
         assert(nlrows[i] != NULL);
         run = SCIPnlrowGetExprtree(nlrows[i]) != NULL && SCIPexprtreeGetNVars(SCIPnlrowGetExprtree(nlrows[i])) > 0;
         run = run || SCIPnlrowGetNQuadVars(nlrows[i]) > 0;
      }
   }

   if( !run )
   {
      SCIPdebugMessage("skipping undercover heuristic: no nonlinear constraints found\n");
      return SCIP_OKAY;
   }

   /* only call heuristics if solving has not stopped yet */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* reset timing, if it was changed temporary (at the root node) */
   if( heurtiming != HEUR_TIMING )
      SCIPheurSetTimingmask(heur, HEUR_TIMING);

   /* call heuristic */
   *result = SCIP_DIDNOTFIND;
   SCIPdebugMessage("calling undercover heuristic for <%s> at depth %d\n", SCIPgetProbName(scip), SCIPgetDepth(scip));

   SCIP_CALL( SCIPapplyUndercover(scip, heur, result, timelimit, memorylimit, nstallnodes) );

   return SCIP_OKAY;
}





/*
 * primal heuristic specific interface methods
 */

/** creates the undercover primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurUndercover(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create undercover primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* always use local bounds */
   heurdata->globalbounds = FALSE;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyUndercover,
         heurFreeUndercover, heurInitUndercover, heurExitUndercover, 
         heurInitsolUndercover, heurExitsolUndercover, heurExecUndercover,
         heurdata) );

   /* add string parameters */
   heurdata->fixingalts = NULL;
   SCIP_CALL( SCIPaddStringParam(scip, "heuristics/"HEUR_NAME"/fixingalts",
         "prioritized sequence of fixing values used ('l'p relaxation, 'n'lp relaxation, 'i'ncumbent solution)",
         &heurdata->fixingalts, FALSE, DEFAULT_FIXINGALTS, NULL, NULL) );

   /* add longint parameters */
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
 
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
 
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   /* add real parameters */
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/conflictweight", 
         "weight for conflict score in fixing order",
         &heurdata->conflictweight, TRUE, DEFAULT_CONFLICTWEIGHT, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/cutoffweight", 
         "weight for cutoff score in fixing order",
         &heurdata->cutoffweight, TRUE, DEFAULT_CUTOFFWEIGHT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/inferenceweight", 
         "weight for inference score in fixing order",
         &heurdata->inferenceweight, TRUE, DEFAULT_INFERENCEWEIGHT, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/maxcoversize",
         "maximum coversize (as fraction of total number of variables)",
         &heurdata->maxcoversize, TRUE, DEFAULT_MAXCOVERSIZE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minimprove",
         "factor by which the heuristic should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, -1.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/recoverdiv",
         "fraction of covering variables in the last cover which need to change their value when re-covering",
         &heurdata->recoverdiv, TRUE, DEFAULT_RECOVERDIV, 0.0, 1.0, NULL, NULL) );

   /* add bool parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/beforecuts",
         "should the heuristic be called at root node before cut separation?",
         &heurdata->beforecuts, TRUE, DEFAULT_BEFORECUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/fixintfirst",
         "should integer variables in the cover be fixed first?",
         &heurdata->fixintfirst, TRUE, DEFAULT_FIXINTFIRST, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/locksrounding",
         "shall LP values for integer vars be rounded according to locks?",
         &heurdata->locksrounding, TRUE, DEFAULT_LOCKSROUNDING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/onlyconvexify",
         "should we only fix variables in order to obtain a convex problem?",
         &heurdata->onlyconvexify, FALSE, DEFAULT_ONLYCONVEXIFY, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/postnlp",
         "should the nlp heuristic be called to polish a feasible solution?",
         &heurdata->postnlp, FALSE, DEFAULT_POSTNLP, NULL, NULL) );

   /* add int parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxbacktracks",
         "maximum number of backtracks in fix-and-propagate",
         &heurdata->maxbacktracks, TRUE, DEFAULT_MAXBACKTRACKS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxrecovers",
         "maximum number of re-coverings",
         &heurdata->maxrecovers, TRUE, DEFAULT_MAXRECOVERS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxreorders",
         "maximum number of reorderings of the fixing order",
         &heurdata->maxreorders, TRUE, DEFAULT_MAXREORDERS, 0, INT_MAX, NULL, NULL) );

   /* add char parameters */
   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/"HEUR_NAME"/coveringobj",
         "objective function of the covering problem ('b'ranching status, influenced nonlinear 'c'onstraints/'t'erms, 'd'omain size, 'l'ocks, 'm'in of up/down locks, 'u'nit penalties, constraint 'v'iolation)",
         &heurdata->coveringobj, TRUE, DEFAULT_COVERINGOBJ, COVERINGOBJS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/copycuts",
         "should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   return SCIP_OKAY;
}

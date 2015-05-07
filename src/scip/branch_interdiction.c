/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_interdiction.c
 * @brief  interdiction branching rule by Lodi, Ralphs, Rossi, and Smriglio
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/branch_interdiction.h"


#define BRANCHRULE_NAME            "interdiction"
#define BRANCHRULE_DESC            "branching rule for (mixed) binary programs by Lodi, Ralphs, Rossi, and Smriglio."
#define BRANCHRULE_PRIORITY     -20000
#define BRANCHRULE_MAXDEPTH          3
#define BRANCHRULE_MAXBOUNDDIST      1.0

#define DEFAULT_COVERSIZE            8
#define DEFAULT_MINLPITER          100LL
#define DEFAULT_NCOVERCANDS         10

#define DEFAULT_CONFLICTWEIGHT       0.01   /**< weight in score calculations for conflict score */
#define DEFAULT_CONFLENGTHWEIGHT     0.0    /**< weight in score calculations for conflict length score*/
#define DEFAULT_INFERENCEWEIGHT      0.1    /**< weight in score calculations for inference score */
#define DEFAULT_CUTOFFWEIGHT         0.0001 /**< weight in score calculations for cutoff score */
#define DEFAULT_PSCOSTWEIGHT         1.0    /**< weight in score calculations for pseudo cost score */

#define DEFAULT_PROPAGATION         TRUE
#define DEFAULT_NPROPROUNDS           -1

/*
 * Data structures
 */

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Longint          maxlpiter;               /**< maximal number of LP iterations per node */
   SCIP_Longint          minlpiter;               /**< minimal number of LP iterations per node */
   SCIP_SOL*             sol;
   SCIP_Bool             solfound;                /**< found an improving solution */
   SCIP_Bool             propagate;               /**< propagate in probing */
   int                   ncovercands;             /**< number of candidates to check in each iteration */
   int                   nproprounds;             /**< maximal number of propagation in probing */
   int                   fails;                   /**< number of fails */
   int                   maxcoversize;            /**< maximal size of an improving solution cover */
   char*                 calcsolcover;            /**< method for calculating the solution cover */
   SCIP_Real             conflictweight;          /**< weight in score calculations for conflict score */
   SCIP_Real             conflengthweight;        /**< weight in score calculations for conflict length score*/
   SCIP_Real             cutoffweight;            /**< weight in score calculations for cutoff score */
   SCIP_Real             inferenceweight;         /**< weight in score calculations for inference score */
   SCIP_Real             pscostweight;            /**< weight in score calculations for pseudo cost score */
};

/* calculate the scores */
static
void calcScore(
    SCIP*                scip,                    /**< SCIP data structure */
    SCIP_BRANCHRULEDATA* branchruledata,          /**< data of branching rule */
    SCIP_Real*           score,                   /**< array to store scores */
    SCIP_VAR**           branchcands,             /**< branching candidates */
    SCIP_Real*           branchcandssol,          /**< solution values */
    SCIP_Real*           branchcandssfrac,        /**< fractional value of variable in current solution */
    int                  nbranchcands             /**< number of branching candidates */
   )
{
   SCIP_Real conflictscore;
   SCIP_Real avgconflictscore;
   SCIP_Real conflengthscore;
   SCIP_Real avgconflengthscore;
   SCIP_Real inferencescore;
   SCIP_Real avginferencescore;
   SCIP_Real cutoffscore;
   SCIP_Real avgcutoffscore;
   SCIP_Real pscostscore;
   SCIP_Real avgpscostscore;
   int c;

   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(score != NULL);
   assert(branchcands != NULL);
   assert(branchcandssol != NULL);

   /* get average conflict, inference, and pseudocost scores */
   avgconflictscore = SCIPgetAvgConflictScore(scip);
   avgconflictscore = MAX(avgconflictscore, 0.1);
   avgconflengthscore = SCIPgetAvgConflictlengthScore(scip);
   avgconflengthscore = MAX(avgconflengthscore, 0.1);
   avginferencescore = SCIPgetAvgInferenceScore(scip);
   avginferencescore = MAX(avginferencescore, 0.1);
   avgcutoffscore = SCIPgetAvgCutoffScore(scip);
   avgcutoffscore = MAX(avgcutoffscore, 0.1);
   avgpscostscore = SCIPgetAvgPseudocostScore(scip);
   avgpscostscore = MAX(avgpscostscore, 0.1);

   for( c = 0; c < nbranchcands; c++ )
   {
      /* get conflict, inference, cutoff, and pseudo cost scores for candidate */
      conflictscore = SCIPgetVarConflictScore(scip, branchcands[c]);
      conflengthscore = SCIPgetVarConflictlengthScore(scip, branchcands[c]);
      inferencescore = SCIPgetVarAvgInferenceScore(scip, branchcands[c]);
      cutoffscore = SCIPgetVarAvgCutoffScore(scip, branchcands[c]);
      pscostscore = SCIPgetVarPseudocostScore(scip, branchcands[c], branchcandssol[c]);

      score[c] = branchruledata->conflictweight * (1.0 - 1.0/(1.0+conflictscore/avgconflictscore))
            + branchruledata->conflengthweight * (1.0 - 1.0/(1.0+conflengthscore/avgconflengthscore))
            + branchruledata->inferenceweight * (1.0 - 1.0/(1.0+inferencescore/avginferencescore))
            + branchruledata->cutoffweight * (1.0 - 1.0/(1.0+cutoffscore/avgcutoffscore))
            + branchruledata->pscostweight * (1.0 - 1.0/(1.0+pscostscore/avgpscostscore));

      /* avoid close to integral variables */
      if( MIN(branchcandssfrac[c], 1.0 - branchcandssfrac[c]) < 10.0*SCIPfeastol(scip) )
         score[c] *= 1e-6;
   }

   return;
}

/* sort candidates in non-increasing score order */
static
void sortCands(
   SCIP*                 scip,
   SCIP_VAR**            branchcands,
   SCIP_Real*            score,
   int                   nbranchcands
   )
{
   int c;

   assert(scip != NULL);
   assert(branchcands != NULL);
   assert(score != NULL);
   assert(nbranchcands >= 0);

   SCIPsortDownRealPtr(score, (void*) branchcands, nbranchcands);

   return;
}

/* calculate an exact solution cover */
static
SCIP_RETCODE calcExactCover(
   SCIP*                 subscip,                 /**< sub-SCIP data structure */
   SCIP_VAR**            vars,                    /**< SCIP variables */
   int                   nvars,                   /**< number of variables */
   SCIP_RESULT*          result                   /**< pointer to store the result */
   )
{
   SCIPdebugMessage("method 'calcExactCover' not implemented yet.\n");

   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}

/* calculate an relaxed solution cover */
static
SCIP_RETCODE calcRelaxCover(
   SCIP*                 subscip,                 /**< sub-SCIP data structure */
   SCIP_VAR**            vars,                    /**< SCIP variables */
   int                   nvars,                   /**< number of variables */
   SCIP_RESULT*          result                   /**< pointer to store the result */
   )
{
   SCIPdebugMessage("method 'calcRelaxCover' not implemented yet.\n");

   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}

/* calculate an heuristically solution cover */
static
SCIP_RETCODE calcHeurCover(
   SCIP*                 scip,                    /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata,          /**< branchrule data */
   SCIP_Real             lpobj,                   /**< objectiv value of the LP relaxation */
   SCIP_VAR**            cover,                   /**< found solution cover */
   int*                  coversize,               /**< size of the cover */
   SCIP_Bool*            cutoff,                  /**< pointer to store whether the all-zero node can be pruned */
   SCIP_Bool*            downinf,                 /**< last down branch is infeasible */
   SCIP_Bool*            upinf,                   /**< last up branch is infeasibke */
   SCIP_RESULT*          result                   /**< pointer to store the result */
   )
{
   SCIP_VAR** branchcands;
   SCIP_Real* branchcandssol;
   SCIP_Real* branchcandsfrac;
   SCIP_Real bestlpobj;
   SCIP_Bool lperror;
   int nbranchcands;
   int niter;

   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(cover != NULL);
   assert(result != NULL);

   niter = 0;
   *coversize = 0;
   *cutoff = FALSE;
   *downinf = FALSE;
   *upinf = FALSE;

   SCIPdebugMessage("start constructing cover, lp threshold = %g\n", lpobj);

   /* start the probing mode */
   SCIP_CALL( SCIPstartProbing(scip) );

   while( *coversize < branchruledata->maxcoversize && !branchruledata->solfound )
   {
      SCIP_CALL( SCIPgetLPBranchCands(scip, &branchcands, &branchcandssol, &branchcandsfrac, &nbranchcands, NULL, NULL) );

      if( nbranchcands > 0 )
      {
         SCIP_Real* score;
         SCIP_Real bestbranchscore;
         SCIP_Real bestlpobj;
         SCIP_Longint lpiters;
         SCIP_Longint lpitersleft;
         SCIP_Longint ndomredsfound;
         SCIP_DOMCHG* domchgs;
         SCIP_NODE* probingnode;
         int bestcand;
         int ncands;
         int c;
         int d;

         /* allocate memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &score, nbranchcands) );

         /* calculate the scores */
         calcScore(scip, branchruledata, score, branchcands, branchcandssol, branchcandsfrac, nbranchcands);

         /* sort candidates in non-increasing score order */
         sortCands(scip, branchcands, score, nbranchcands);

         lpiters = SCIPgetNLPIterations(scip);
         lpitersleft = branchruledata->maxlpiter;

         /* open a new probing node */
         SCIP_CALL( SCIPnewProbingNode(scip) );

         probingnode = SCIPgetCurrentNode(scip);
         assert(SCIPnodeGetType(probingnode) == SCIP_NODETYPE_PROBINGNODE);

         ncands = MIN(nbranchcands, branchruledata->ncovercands);
         bestbranchscore = SCIP_REAL_MIN;
         bestcand = -1;

         for( c = 0; c < ncands && lpitersleft > 0; c++ )
         {
            SCIP_Real dwbranchobj;
            SCIP_Real upbranchobj;
            SCIP_Real branchscore;
            int ndomchgs;
            int b;

            for( b = 0; b <= 1; b++ )
            {
               if( b == 0 )
               {
                  /** test the 0-branch */
                  SCIP_CALL( SCIPfixVarProbing(scip, branchcands[c], 0.0) );
               }
               else
               {
                  SCIP_CALL( SCIPfixVarProbing(scip, branchcands[c], 1.0) );
               }

               /* run propagation */
               ndomredsfound = 0;
               if( branchruledata->propagate )
               {
                  SCIP_CALL( SCIPpropagateProbing(scip, branchruledata->nproprounds, cutoff, &ndomredsfound) );
               }

               if( *cutoff )
               {
                  assert(branchruledata->propagate);

                  SCIPdebugMessage("   cutoff after propagation on branch\n");

                  if( b == 0 )
                     *downinf = TRUE;
                  else
                     *upinf = TRUE;

                  /* we add the variable but do not create the all-zero node */
                  cover[*coversize] = branchcands[c];
                  (*coversize)++;

                  SCIPfreeBufferArray(scip, &score);
                  goto TERMINATE;
               }

               /* solve the probing LP */
               SCIP_CALL( SCIPsolveProbingLP(scip, lpitersleft, &lperror, cutoff) );

               if( b == 0 )
                  dwbranchobj = SCIPgetLPObjval(scip);
               else
                  upbranchobj = SCIPgetLPObjval(scip);

               /* we stop at lperror and cutoff */
               if( lperror )
               {
                  SCIPfreeBufferArray(scip, &score);
                  goto TERMINATE;
               }

               if( *cutoff )
               {
                  if( b == 0 )
                     *downinf = TRUE;
                  else
                     *upinf = TRUE;

                  /* we add the variable but do not create the all-zero node */
                  cover[*coversize] = branchcands[c];
                  (*coversize)++;

                  SCIPfreeBufferArray(scip, &score);
                  goto TERMINATE;
               }

               /* undo the bound change */
               SCIP_CALL( SCIPbacktrackProbing(scip, *coversize) );

               lpitersleft = lpitersleft - (SCIPgetNLPIterations(scip) - lpiters);
               lpiters = SCIPgetNLPIterations(scip);

            }

            /* calculate the score and check whether the candidate is better */
            branchscore = SCIPgetBranchScore(scip, NULL, dwbranchobj, upbranchobj);

            if( lpobj < dwbranchobj && branchscore > bestbranchscore )
            {
               bestbranchscore = branchscore;
               bestcand = c;
            }

         }

         if( bestcand != -1 )
         {
            SCIPdebugMessage("   found cover variable <%s> with branchscore %g\n", SCIPvarGetName(branchcands[bestcand]), bestbranchscore);

            /** open 0-branch */
            SCIP_CALL( SCIPfixVarProbing(scip, branchcands[c], 0.0) );

            cover[*coversize] = branchcands[c];
            (*coversize)++;
            break;
         }
         else
         {
            /* we can terminate since no candidate is good enough */
            goto TERMINATE;
         }

         SCIPfreeBufferArray(scip, &score);

         if( c == nbranchcands || lpitersleft <= 0 )
         {
            *result = SCIP_DIDNOTFIND;
            branchruledata->fails++;
            goto TERMINATE;
         }
      }
      else
      {
         SCIP_Bool stored;

         /* get the LP solution */
         SCIP_CALL( SCIPlinkLPSol(scip, branchruledata->sol) );

         /* try the solution */
         SCIP_CALL( SCIPtrySol(scip, branchruledata->sol, FALSE, FALSE, FALSE, FALSE, &stored) );

         if( stored )
         {
            SCIPdebugMessage("   found new incumbent\n");
            branchruledata->solfound = TRUE;
            *cutoff = TRUE;
            break;
         }
      }

      niter++;
   }

   *result = SCIP_SUCCESS;

  TERMINATE:

   /* end probing mode */
   SCIP_CALL( SCIPendProbing(scip) );

   return SCIP_OKAY;
}

/* perform the branching */
static
SCIP_RETCODE branchInterdiction(
   SCIP*                 scip,                    /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata,          /**< data of branching rule */
   SCIP_VAR**            cover,                   /**< solution cover */
   int                   coversize,               /**< size of the cover */
   SCIP_Bool             cutoff,                  /**< true, if the all-zero node can be pruned */
   SCIP_RESULT*          result                   /**< pointer to store the result */
   )
{
   int c;

   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(cover != NULL);
   assert(coversize > 0);
   assert(result != NULL);

   SCIPdebugMessage("perform interdiction branching: create %d nodes\n", coversize + 1 - cutoff);

   if( coversize == 1 && cutoff )
   {
      SCIP_VAR* var;

      assert(cover[0] != NULL);
      var = SCIPfindVar(scip, SCIPvarGetName(cover[0]));
      assert(var != NULL);

      /* fix variable to 1.0 */
      SCIP_CALL( SCIPchgVarLbNode(scip, SCIPgetCurrentNode(scip), var, 1.0) );

      *result = SCIP_REDUCEDDOM;
   }
   else
   {
      c = 0;

      while( c <= coversize )
      {
         SCIP_NODE* child;
         SCIP_VAR* var;
         int v;

         SCIP_CALL( SCIPcreateChild(scip, &child, 0.0, SCIPgetLocalTransEstimate(scip)) );
         assert(child != NULL);

         /* change fist c-1 bounds to 0 */
         for( v = 0; v < c && (c < coversize || !cutoff); v++ )
         {
            int d;

            assert(cover[v] != NULL);

            var = SCIPfindVar(scip, SCIPvarGetName(cover[v]));
            assert(var != NULL);

            SCIP_CALL( SCIPchgVarUbNode(scip, child, var, 0.0) );
         }

         if( c < coversize )
         {
            /* change the c-th bound to 1 */
            assert(cover[c] != NULL);

            var = SCIPfindVar(scip, SCIPvarGetName(cover[c]));
            assert(var != NULL);

            SCIP_CALL( SCIPchgVarLbNode(scip, child, var, 1.0) );
         }
         c++;

         if( c == coversize && cutoff )
            break;
      }

      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}

/** perform interdiction branching for 'normal' nodes */
static
SCIP_RETCODE execInterdiction(
   SCIP*                 scip,                    /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,              /**< branching rule */
   SCIP_NODE*            curnode,                 /**< current focus node*/
   SCIP_RESULT*          result                   /**< pointer to store the result */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** cover;
   SCIP_VAR** vars;
   SCIP_Real lpobj;
   SCIP_Bool cutoff;
   SCIP_Bool downinf;
   SCIP_Bool upinf;
   int nfracbins;
   int coversize;
   int nvars;
   int d;
   int v;

   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(curnode != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   nfracbins = 0;
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   SCIP_CALL( SCIPgetLPBranchCands(scip, NULL, NULL, NULL, &nfracbins, NULL, NULL) );

   if( nfracbins <= 1 )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* get the number of LP iterations */
   branchruledata->maxlpiter = (SCIP_Longint)(SCIPgetNRootLPIterations(scip)/10);
   branchruledata->maxlpiter /= MAX(1, branchruledata->fails);

   if( branchruledata->solfound )
      branchruledata->maxlpiter *= 2;

   if( branchruledata->maxlpiter < branchruledata->minlpiter )
      return SCIP_OKAY;

   branchruledata->solfound = FALSE;

   SCIPdebugMessage("%d fractional binaries, LP iteration limit = %"SCIP_LONGINT_FORMAT", min. LP iterations = %"SCIP_LONGINT_FORMAT", max cover size = %d\n",
         nfracbins, branchruledata->maxlpiter, DEFAULT_MINLPITER, branchruledata->maxcoversize);

   if( branchruledata->maxlpiter < DEFAULT_MINLPITER )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* allocate memory for the solution cover and bound changes found during propagation */
   SCIP_CALL( SCIPallocBufferArray(scip, &cover, branchruledata->maxcoversize) );

   lpobj = MAX(SCIPgetAvgLowerbound(scip), SCIPnodeGetLowerbound(curnode));
   cutoff = FALSE;

   if( strcmp(branchruledata->calcsolcover, "e") == 0 )
   {
      SCIP_CALL( calcExactCover(scip, vars, nvars, result) );
   }
   else if( strcmp(branchruledata->calcsolcover, "r") == 0 )
   {
      SCIP_CALL( calcRelaxCover(scip, vars, nvars, result) );
   }
   else
   {
      assert(strcmp(branchruledata->calcsolcover, "h") == 0);
      SCIP_CALL( calcHeurCover(scip, branchruledata, lpobj, cover, &coversize, &cutoff, &downinf, &upinf, result) );
   }
   assert(*result == SCIP_DIDNOTFIND || *result == SCIP_SUCCESS);

   if( coversize > 0 )
   {
      SCIP_CALL( branchInterdiction(scip, branchruledata, cover, coversize, cutoff, result) );
      assert(*result == SCIP_BRANCHED || *result == SCIP_REDUCEDDOM);
   }
   else
      *result = SCIP_DIDNOTFIND;

   /* free subproblem */
   SCIPfreeBufferArray(scip, &cover);

   return SCIP_OKAY;
}

/** perform interdiction branching for reoptimization nodes */
static
SCIP_RETCODE execReoptInterdiction(
   SCIP*                 scip,                    /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,              /**< branching rule */
   SCIP_NODE*            curnode,                 /**< current focus node*/
   SCIP_RESULT*          result                   /**< pointer to store the result */
   )
{
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(curnode != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* do not run is the constraint need not be split */
   if( !SCIPnodeSplit(scip, curnode) )
      return SCIP_OKAY;

   // HIER WEITER ! ! !

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitInterdiction)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   if( strcmp(branchruledata->calcsolcover, "e") != 0
    && strcmp(branchruledata->calcsolcover, "r") != 0
    && strcmp(branchruledata->calcsolcover, "h") != 0 )
   {
      SCIPerrorMessage("unknown parameter: %s\n", branchruledata->calcsolcover);
      SCIPABORT();
   }

   branchruledata->fails = 0;
   branchruledata->maxlpiter = 0;
   branchruledata->solfound = FALSE;
   SCIP_CALL( SCIPcreateSol(scip, &branchruledata->sol, NULL) );

   return SCIP_OKAY;
}

/** deinitialization method of branching rule (called before transformed problem is freed) */
static
SCIP_DECL_BRANCHEXIT(branchExitInterdiction)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);

   SCIP_CALL( SCIPfreeSol(scip, &branchruledata->sol) );

   return SCIP_OKAY;
}

/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#define branchInitsolInterdiction NULL

/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#define branchExitsolInterdiction NULL

/** branching execution method for external candidates */
#define branchExecextInterdiction NULL

/** branching execution method for not completely fixed pseudo solutions */
#define branchExecpsInterdiction NULL

/** copy method for branchrule plugins (called when SCIP copies plugins) */
#define branchCopyInterdiction NULL

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeInterdiction)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);

   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpInterdiction)
{  /*lint --e{715}*/
   SCIP_NODE* curnode;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   curnode = SCIPgetCurrentNode(scip);
   assert(curnode != NULL);

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* we do not run if any variable is binary */
   if( SCIPgetNBinVars(scip) == 0 )
      return SCIP_OKAY;

   /* we do not run if integer variables exist */
   if( SCIPgetNIntVars(scip) != 0 )
      return SCIP_OKAY;

   SCIPdebugMessage("Execlp method of %s branching in node %llu\n", BRANCHRULE_NAME, SCIPnodeGetNumber(curnode));

   if( SCIPnodeGetReoptID(curnode) >= 0 )
   {
      SCIP_CALL( execReoptInterdiction(scip, branchrule, curnode, result) );
   }

   if( *result == SCIP_DIDNOTRUN )
   {
      SCIP_CALL( execInterdiction(scip, branchrule, curnode, result) );
   }

   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the interdiction branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleInterdiction(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create interdiction branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitInterdiction) );
   SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitInterdiction) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeInterdiction) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpInterdiction) );

   /* parameters */
   branchruledata->calcsolcover = NULL;
   SCIP_CALL( SCIPaddStringParam(scip, "branching/"BRANCHRULE_NAME"/calcsolcover", "calculate an 'e'xect, 'r'elaxed, or 'h'euristically solution cover",
         &branchruledata->calcsolcover, FALSE, "h", NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "branching/"BRANCHRULE_NAME"/ncovercands", "number of candidates to check in each iteration",
         &branchruledata->ncovercands, FALSE, DEFAULT_NCOVERCANDS, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "branching/"BRANCHRULE_NAME"/maxcoversize", "maximal size of a solution cover",
         &branchruledata->maxcoversize, FALSE, DEFAULT_COVERSIZE, 2, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "branching/"BRANCHRULE_NAME"/minlpiter", "minimal number of LP iterations for executing the branching rule",
         &branchruledata->minlpiter, TRUE, DEFAULT_MINLPITER, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "branching/"BRANCHRULE_NAME"/conflictweight", "weight for conflict score",
         &branchruledata->conflictweight, TRUE, DEFAULT_CONFLICTWEIGHT, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "branching/"BRANCHRULE_NAME"/conflengthweight", "weight for conflict score",
         &branchruledata->conflengthweight, TRUE, DEFAULT_CONFLENGTHWEIGHT, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "branching/"BRANCHRULE_NAME"/cutoffweight", "weight for cutoff length score",
         &branchruledata->cutoffweight, TRUE, DEFAULT_CUTOFFWEIGHT, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "branching/"BRANCHRULE_NAME"/inferenceweight", "weight for inference score",
         &branchruledata->inferenceweight, TRUE, DEFAULT_INFERENCEWEIGHT, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "branching/"BRANCHRULE_NAME"/pscostweight", "weight for pseudo-cost score",
         &branchruledata->pscostweight, TRUE, DEFAULT_PSCOSTWEIGHT, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/"BRANCHRULE_NAME"/propagate", "use propagation while constructing the cover",
         &branchruledata->propagate, TRUE, DEFAULT_PROPAGATION, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "branching/"BRANCHRULE_NAME"/nproprounds", "maximal number of propagation rounds",
         &branchruledata->nproprounds, TRUE, DEFAULT_NPROPROUNDS, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

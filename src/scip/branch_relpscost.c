/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: branch_relpscost.c,v 1.36 2005/07/15 17:20:04 bzfpfend Exp $"

/**@file   branch_relpscost.c
 * @brief  reliable pseudo costs branching rule
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_relpscost.h"


#define BRANCHRULE_NAME          "relpscost"
#define BRANCHRULE_DESC          "reliability branching on pseudo cost values"
#define BRANCHRULE_PRIORITY      10000
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

#define DEFAULT_MINRELIABLE      1.0    /**< minimal value for minimum pseudo cost size to regard pseudo cost value as reliable */
#define DEFAULT_MAXRELIABLE      8.0    /**< maximal value for minimum pseudo cost size to regard pseudo cost value as reliable */
#define DEFAULT_SBITERQUOT       0.5    /**< maximal fraction of strong branching LP iterations compared to normal iters */
#define DEFAULT_SBITEROFS   100000      /**< additional number of allowed strong branching LP iterations */
#define DEFAULT_MAXLOOKAHEAD     8      /**< maximal number of further variables evaluated without better score */
#define DEFAULT_INITCAND       100      /**< maximal number of candidates initialized with strong branching per node */
#define DEFAULT_INITITER         0      /**< iteration limit for strong branching init of pseudo cost entries (0: auto) */
#define DEFAULT_MAXBDCHGS        5      /**< maximal number of bound tightenings before the node is reevaluated (-1: unlimited) */


/** branching rule data */
struct BranchruleData
{
   Real             minreliable;        /**< minimal value for minimum pseudo cost size to regard pseudo cost value as reliable */
   Real             maxreliable;        /**< maximal value for minimum pseudo cost size to regard pseudo cost value as reliable */
   Real             sbiterquot;         /**< maximal fraction of strong branching LP iterations compared to normal iters */
   int              sbiterofs;          /**< additional number of allowed strong branching LP iterations */
   int              maxlookahead;       /**< maximal number of further variables evaluated without better score */
   int              initcand;           /**< maximal number of candidates initialized with strong branching per node */
   int              inititer;           /**< iteration limit for strong branching init of pseudo cost entries (0: auto) */
   int              maxbdchgs;          /**< maximal number of bound tightenings before the node is reevaluated (-1: unlimited) */
};



/*
 * local methods
 */

/** adds given index and direction to bound change arrays */
static
RETCODE addBdchg(
   SCIP*            scip,               /**< SCIP data structure */
   int**            bdchginds,          /**< pointer to bound change index array */
   Bool**           bdchgdowninfs,      /**< pointer to bound change direction array */
   int*             nbdchgs,            /**< pointer to number of bound changes */
   int              ind,                /**< index to store in bound change array */
   Bool             downinf             /**< is the down branch infeasible? */
   )
{
   assert(bdchginds != NULL);
   assert(bdchgdowninfs != NULL);
   assert(nbdchgs != NULL);

   CHECK_OKAY( SCIPreallocBufferArray(scip, bdchginds, (*nbdchgs) + 1) );
   CHECK_OKAY( SCIPreallocBufferArray(scip, bdchgdowninfs, (*nbdchgs) + 1) );
   assert(*bdchginds != NULL);
   assert(*bdchgdowninfs != NULL);

   (*bdchginds)[*nbdchgs] = ind;
   (*bdchgdowninfs)[*nbdchgs] = downinf;
   (*nbdchgs)++;

   return SCIP_OKAY;
}

/** frees bound change arrays */
static
void freeBdchgs(
   SCIP*            scip,               /**< SCIP data structure */
   int**            bdchginds,          /**< pointer to bound change index array */
   Bool**           bdchgdowninfs,      /**< pointer to bound change direction array */
   int*             nbdchgs             /**< pointer to number of bound changes */
   )
{
   assert(bdchginds != NULL);
   assert(bdchgdowninfs != NULL);
   assert(nbdchgs != NULL);

   SCIPfreeBufferArrayNull(scip, bdchgdowninfs);
   SCIPfreeBufferArrayNull(scip, bdchginds);
   *nbdchgs = 0;
}

/** applies bound changes stored in bound change arrays */
static
RETCODE applyBdchgs(
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            lpcands,            /**< fractional branching candidates */
   Real*            lpcandssol,         /**< LP solution array of branching candidates */
   int*             bdchginds,          /**< bound change index array */
   Bool*            bdchgdowninfs,      /**< bound change direction array */
   int              nbdchgs             /**< number of bound changes */
   )
{
   int i;

   debugMessage("applying %d bound changes\n", nbdchgs);

   for( i = 0; i < nbdchgs; ++i )
   {
      int c;

      c = bdchginds[i];
      if( bdchgdowninfs[i] )
      {
         /* downwards rounding is infeasible -> change lower bound of variable to upward rounding */
         CHECK_OKAY( SCIPchgVarLb(scip, lpcands[c], SCIPfeasCeil(scip, lpcandssol[c])) );
      }
      else
      {
         /* upwards rounding is infeasible -> change upper bound of variable to downward rounding */
         CHECK_OKAY( SCIPchgVarUb(scip, lpcands[c], SCIPfeasFloor(scip, lpcandssol[c])) );
      }
      debugMessage(" -> <%s> (sol:%g) -> [%g,%g]\n",
         SCIPvarGetName(lpcands[c]), lpcandssol[c], SCIPvarGetLbLocal(lpcands[c]), SCIPvarGetUbLocal(lpcands[c]));
   }

   return SCIP_OKAY;
}




/*
 * Callback methods
 */

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
DECL_BRANCHFREE(branchFreeRelpscost)
{  /*lint --e{715}*/
   BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
#define branchInitRelpscost NULL


/** deinitialization method of branching rule (called before transformed problem is freed) */
#define branchExitRelpscost NULL


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#define branchInitsolRelpscost NULL


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#define branchExitsolRelpscost NULL


/** branching execution method for fractional LP solutions */
static
DECL_BRANCHEXECLP(branchExeclpRelpscost)
{  /*lint --e{715}*/
   BRANCHRULEDATA* branchruledata;
   VAR** lpcands;
   Real* lpcandssol;
   Real* lpcandsfrac;
   Real lpobjval;
   Real cutoffbound;
   Real bestsbdown;
   Real bestsbup;
   Real provedbound;
   Bool bestsbdownvalid;
   Bool bestsbupvalid;
   Bool bestisstrongbranch;
   Bool allcolsinlp;
   Bool exactsolve;
   int nlpcands;
   int ninitcands;
   int bestcand;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   debugMessage("Execlp method of relpscost branching\n");

   *result = SCIP_DIDNOTRUN;

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get current LP objective bound of the local sub problem and global cutoff bound */
   lpobjval = SCIPgetLPObjval(scip);
   cutoffbound = SCIPgetCutoffbound(scip);

   /* check, if we want to solve the problem exactly, meaning that strong branching information is not useful
    * for cutting off sub problems and improving lower bounds of children
    */
   exactsolve = SCIPisExactSolve(scip);

   /* check, if all existing columns are in LP, and thus the strong branching results give lower bounds */
   allcolsinlp = SCIPallColsInLP(scip);

   /* get branching candidates */
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, NULL, &nlpcands) );
   assert(nlpcands > 0);

   bestcand = -1;
   bestisstrongbranch = FALSE;
   bestsbdown = SCIP_INVALID;
   bestsbup = SCIP_INVALID;
   bestsbdownvalid = FALSE;
   bestsbupvalid = FALSE;
   provedbound = lpobjval;

   if( nlpcands == 1 )
   {
      /* only one candidate: nothing has to be done */
      bestcand = 0;
      ninitcands = 0;
   }
   else
   {
      int* initcands;
      Real* initcandscores;
      int* bdchginds;
      Bool* bdchgdowninfs;
      int maxninitcands;
      int nuninitcands;
      int nbdchgs;
      int nbdconflicts;
      Real bestpsscore;
      Real bestpsfracscore;
      Real bestpsdomainscore;
      Real bestsbscore;
      Real bestuninitsbscore;
      Real bestsbfracscore;
      Real bestsbdomainscore;
      Real prio;
      Real reliable;
      Real maxlookahead;
      Real lookahead;
      Longint nodenum;
      Longint nsblpiterations;
      Longint maxnsblpiterations;
      int maxbdchgs;
      int bestpscand;
      int bestsbcand;
      int inititer;
      int i;
      int c;

      /* get maximal number of candidates to initialize with strong branching; if the current solutions is not basic,
       * we cannot apply the simplex algorithm and therefore don't initialize any candidates
       */
      maxninitcands = MIN(nlpcands, branchruledata->initcand);
      if( !SCIPisLPSolBasic(scip) )
         maxninitcands = 0;

      /* calculate maximal number of strong branching LP iterations; if we used too many, don't apply strong branching
       * any more
       */
      maxnsblpiterations = (Longint)(branchruledata->sbiterquot * SCIPgetNNodeLPIterations(scip))
         + branchruledata->sbiterofs + SCIPgetNRootStrongbranchLPIterations(scip);
      if( SCIPgetNStrongbranchLPIterations(scip) > maxnsblpiterations )
         maxninitcands = 0;

      /* get buffer for storing the unreliable candidates */
      CHECK_OKAY( SCIPallocBufferArray(scip, &initcands, maxninitcands+1) ); /* allocate one additional slot for convenience */
      CHECK_OKAY( SCIPallocBufferArray(scip, &initcandscores, maxninitcands+1) );
      ninitcands = 0;

      /* get current node number */
      nodenum = SCIPgetNNodes(scip);

      /* initialize bound change arrays */
      bdchginds = NULL;
      bdchgdowninfs = NULL;
      nbdchgs = 0;
      nbdconflicts = 0;
      maxbdchgs = branchruledata->maxbdchgs;

      /* calculate value used as reliability */
      nsblpiterations = SCIPgetNStrongbranchLPIterations(scip);
      prio = (maxnsblpiterations - nsblpiterations)/(nsblpiterations + 1.0);
      reliable = (1.0-prio) * branchruledata->minreliable + prio * branchruledata->maxreliable;

      /* search for the best pseudo cost candidate, while remembering unreliable candidates in a sorted buffer */
      nuninitcands = 0;
      bestpscand = -1;
      bestpsscore = -SCIPinfinity(scip);
      bestpsfracscore = -SCIPinfinity(scip);
      bestpsdomainscore = -SCIPinfinity(scip);
      for( c = 0; c < nlpcands; ++c )
      {
         Real score;
         Bool usesb;

         assert(lpcands[c] != NULL);
         assert(!SCIPisFeasIntegral(scip, lpcandssol[c]));

         /* get pseudo cost score of candidate */
         score = SCIPgetVarPseudocostScoreCurrentRun(scip, lpcands[c], lpcandssol[c]);
         usesb = FALSE;

         /* don't use strong branching on variables that have already been initialized at the current node */
         if( SCIPgetVarStrongbranchNode(scip, lpcands[c]) == nodenum )
         {
            Real down;
            Real up;
            Real lastlpobjval;
            Real downgain;
            Real upgain;

            /* use the score of the strong branching call at the current node */
            CHECK_OKAY( SCIPgetVarStrongbranchLast(scip, lpcands[c], &down, &up, NULL, NULL, NULL, &lastlpobjval) );
            downgain = MAX(down - lastlpobjval, 0.0);
            upgain = MAX(up - lastlpobjval, 0.0);
            score = SCIPgetBranchScore(scip, lpcands[c], downgain, upgain);

            debugMessage(" -> strong branching on variable <%s> already performed (down=%g (%+g), up=%g (%+g), score=%g)\n",
               SCIPvarGetName(lpcands[c]), down, downgain, up, upgain, score);
         }
         else if( maxninitcands > 0 )
         {
            Real downsize;
            Real upsize;
            Real size;

            /* check, if the pseudo cost score of the variable is reliable */
            downsize = SCIPgetVarPseudocostCountCurrentRun(scip, lpcands[c], SCIP_BRANCHDIR_DOWNWARDS);
            upsize = SCIPgetVarPseudocostCountCurrentRun(scip, lpcands[c], SCIP_BRANCHDIR_UPWARDS);
            size = MIN(downsize, upsize);

            /* use strong branching on variables with unreliable pseudo cost scores */
            usesb = (size < reliable);

            /* count the number of variables that are completely uninitialized */
            if( size < 0.1 )
               nuninitcands++;
         }

         if( usesb )
         {
            int j;

            /* pseudo cost of variable is not reliable: insert candidate in initcands buffer */
            for( j = ninitcands; j > 0 && score > initcandscores[j-1]; --j )
            {
               initcands[j] = initcands[j-1];
               initcandscores[j] = initcandscores[j-1];
            }
            initcands[j] = c;
            initcandscores[j] = score;
            ninitcands++;
            ninitcands = MIN(ninitcands, maxninitcands);
         }
         else
         {
            /* variable will keep it's pseudo cost value: check for better pseudo cost score of candidate */
            if( SCIPisSumGE(scip, score, bestpsscore) )
            {
               Real fracscore;
               Real domainscore;

               fracscore = MIN(lpcandsfrac[c], 1.0 - lpcandsfrac[c]);
               domainscore = -(SCIPvarGetUbLocal(lpcands[c]) - SCIPvarGetLbLocal(lpcands[c]));
               if( SCIPisSumGT(scip, score, bestpsscore)
                  || SCIPisSumGT(scip, fracscore, bestpsfracscore)
                  || (SCIPisSumGE(scip, fracscore, bestpsfracscore) && domainscore > bestpsdomainscore) )
               {
                  bestpscand = c;
                  bestpsscore = score;
                  bestpsfracscore = fracscore;
                  bestpsdomainscore = domainscore;
               }
            }
         }
      }

      /* initialize unreliable candidates with strong branching until maxlookahead is reached,
       * search best strong branching candidate
       */
      maxlookahead = (Real)branchruledata->maxlookahead * (1.0 + (Real)nuninitcands/(Real)nlpcands);
      inititer = branchruledata->inititer;
      if( inititer == 0 )
      {
         Longint nlpiterations;
         int nlps;

         /* iteration limit is set to twice the average number of iterations spent to resolve a dual feasible LP;
          * at the first few nodes, this average is not very exact, so we better increase the iteration limit on
          * these very important nodes
          */
         nlpiterations = SCIPgetNDualResolveLPIterations(scip);
         nlps = SCIPgetNDualResolveLPs(scip);
         if( nlps == 0 )
         {
            nlpiterations = SCIPgetNNodeInitLPIterations(scip);
            nlps = SCIPgetNNodeInitLPs(scip);
            if( nlps == 0 )
            {
               nlpiterations = 1000;
               nlps = 1;
            }
         }
         assert(nlps >= 1);
         inititer = (int)(2*nlpiterations / nlps);
         inititer = (int)((Real)inititer * (1.0 + 20.0/nodenum));
         inititer = MAX(inititer, 10);
         inititer = MIN(inititer, 10000);
      }
      
      debugMessage("strong branching (reliable=%g, %d/%d cands, %d uninit, maxcands=%d, maxlookahead=%g, inititer=%d, iters:%lld/%lld, basic:%d)\n",
         reliable, ninitcands, nlpcands, nuninitcands, maxninitcands, maxlookahead, inititer, 
         SCIPgetNStrongbranchLPIterations(scip), maxnsblpiterations, SCIPisLPSolBasic(scip));

      bestsbcand = -1;
      bestsbscore = -SCIPinfinity(scip);
      bestsbfracscore = -SCIPinfinity(scip);
      bestsbdomainscore = -SCIPinfinity(scip);
      lookahead = 0.0;
      for( i = 0; i < ninitcands && lookahead < maxlookahead
              && (i < maxlookahead || SCIPgetNStrongbranchLPIterations(scip) < maxnsblpiterations); ++i )
      {
         Real down;
         Real up;
         Real downgain;
         Real upgain;
         Real score;
         Bool downvalid;
         Bool upvalid;
         Bool lperror;
         Bool downinf;
         Bool upinf;
         Bool downconflict;
         Bool upconflict;

         /* get candidate number to initialize */
         c = initcands[i];
         assert(!SCIPisFeasIntegral(scip, lpcandssol[c]));

         debugMessage("init pseudo cost (%g/%g) of <%s> at %g (score:%g) with strong branching (%d iters) -- %lld/%lld iterations\n",
            SCIPgetVarPseudocostCountCurrentRun(scip, lpcands[c], SCIP_BRANCHDIR_DOWNWARDS), 
            SCIPgetVarPseudocostCountCurrentRun(scip, lpcands[c], SCIP_BRANCHDIR_UPWARDS), 
            SCIPvarGetName(lpcands[c]), lpcandssol[c], initcandscores[i],
            inititer, SCIPgetNStrongbranchLPIterations(scip), maxnsblpiterations);

         /* use strong branching on candidate */
         CHECK_OKAY( SCIPgetVarStrongbranch(scip, lpcands[c], inititer, 
               &down, &up, &downvalid, &upvalid, &downinf, &upinf, &downconflict, &upconflict, &lperror) );

         /* check for an error in strong branching */
         if( lperror )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "(node %lld) error in strong branching call for variable <%s> with solution %g\n", 
               SCIPgetNNodes(scip), SCIPvarGetName(lpcands[c]), lpcandssol[c]);
            break;
         }

         /* evaluate strong branching */
         down = MAX(down, lpobjval);
         up = MAX(up, lpobjval);
         downgain = down - lpobjval;
         upgain = up - lpobjval;
         assert(!allcolsinlp || exactsolve || !downvalid || downinf == SCIPisGE(scip, down, cutoffbound));
         assert(!allcolsinlp || exactsolve || !upvalid || upinf == SCIPisGE(scip, up, cutoffbound));
         assert(downinf || !downconflict);
         assert(upinf || !upconflict);

         /* the minimal lower bound of both children is a proved lower bound of the current subtree */
         if( allcolsinlp && !exactsolve && downvalid && upvalid )
         {
            Real minbound;
            
            minbound = MIN(down, up);
            provedbound = MAX(provedbound, minbound);
         }

         /* check if there are infeasible roundings */
         if( downinf || upinf )
         {
            assert(allcolsinlp);
            assert(!exactsolve);
            
            /* if for both infeasibilities, a conflict clause was created, we don't need to fix the variable by hand,
             * but better wait for the next propagation round to fix them as an inference, and potentially produce a
             * cutoff that can be analyzed
             */
            if( allowaddcons && downinf == downconflict && upinf == upconflict )
            {
               debugMessage(" -> variable <%s> is infeasible in %s: conflict constraint added\n",
                  SCIPvarGetName(lpcands[c]), 
                  downinf && upinf ? "both directions" : (downinf ? "downward branch" : "upwardbranch"));
               *result = SCIP_CONSADDED;
               nbdconflicts++;
               if( (downinf && upinf)
                  || (maxbdchgs >= 0 && nbdchgs + nbdconflicts >= maxbdchgs) )
                  break; /* terminate initialization loop, because enough roundings are performed or a cutoff was found */
            }
            else if( downinf && upinf )
            {
               /* both roundings are infeasible -> node is infeasible */
               debugMessage(" -> variable <%s> is infeasible in both directions (conflict: %d/%d)\n",
                  SCIPvarGetName(lpcands[c]), downconflict, upconflict);
               *result = SCIP_CUTOFF;
               break; /* terminate initialization loop, because node is infeasible */
            }
            else
            {
               /* rounding is infeasible in one direction -> round variable in other direction */
               debugMessage(" -> variable <%s> is infeasible in %s branch (conflict: %d/%d)\n",
                  SCIPvarGetName(lpcands[c]), downinf ? "downward" : "upward", downconflict, upconflict);
               CHECK_OKAY( addBdchg(scip, &bdchginds, &bdchgdowninfs, &nbdchgs, c, downinf) );
               if( maxbdchgs >= 0 && nbdchgs + nbdconflicts >= maxbdchgs )
                  break; /* terminate initialization loop, because enough roundings are performed */
            }
         }
         else
         {
            /* check for a better score */
            score = SCIPgetBranchScore(scip, lpcands[c], downgain, upgain);
            if( SCIPisSumGE(scip, score, bestsbscore) )
            {
               Real fracscore;
               Real domainscore;

               fracscore = MIN(lpcandsfrac[c], 1.0 - lpcandsfrac[c]);
               domainscore = -(SCIPvarGetUbLocal(lpcands[c]) - SCIPvarGetLbLocal(lpcands[c]));
               if( SCIPisSumGT(scip, score, bestsbscore)
                  || SCIPisSumGT(scip, fracscore, bestsbfracscore)
                  || (SCIPisSumGE(scip, fracscore, bestsbfracscore) && domainscore > bestsbdomainscore) )
               {
                  bestsbcand = c;
                  bestsbscore = score;
                  bestsbdown = down;
                  bestsbup = up;
                  bestsbdownvalid = downvalid;
                  bestsbupvalid = upvalid;
                  bestsbfracscore = fracscore;
                  bestsbdomainscore = domainscore;
                  lookahead = 0.0;
               }
               else
                  lookahead += 0.5;
            }
            else
               lookahead += 1.0;
         
            /* update pseudo cost values */
            CHECK_OKAY( SCIPupdateVarPseudocost(scip, lpcands[c], 0.0-lpcandsfrac[c], downgain, 1.0) );
            CHECK_OKAY( SCIPupdateVarPseudocost(scip, lpcands[c], 1.0-lpcandsfrac[c], upgain, 1.0) );

            debugMessage(" -> variable <%s> (solval=%g, down=%g (%+g), up=%g (%+g), score=%g) -- best: <%s> (%g / %g / %g), lookahead=%g/%g\n",
               SCIPvarGetName(lpcands[c]), lpcandssol[c], down, downgain, up, upgain, score, 
               SCIPvarGetName(lpcands[bestsbcand]), bestsbscore, bestsbfracscore, bestsbdomainscore, 
               lookahead, maxlookahead);
         }
      }

      /* get the score of the best uninitialized strong branching candidate */
      if( i < ninitcands )
         bestuninitsbscore = initcandscores[i];
      else
         bestuninitsbscore = -SCIPinfinity(scip);

      /* if the best pseudo cost candidate is better than the best uninitialized strong branching candidate,
       * compare it to the best initialized strong branching candidate
       */
      if( bestpsscore > bestuninitsbscore && SCIPisSumGT(scip, bestpsscore, bestsbscore) )
      {
         bestcand = bestpscand;
         bestisstrongbranch = FALSE;
      }
      else if( bestsbcand >= 0 )
      {
         bestcand = bestsbcand;
         bestisstrongbranch = TRUE;
      }
      else
      {
         /* no candidate was initialized, and the best score is the one of the first candidate in the initialization
          * queue
          */
         assert(ninitcands >= 1);
         bestcand = initcands[0];
         bestisstrongbranch = FALSE;
      }

      /* apply domain reductions */
      if( nbdchgs > 0 )
      {
         if( *result != SCIP_CUTOFF )
         {
            CHECK_OKAY( applyBdchgs(scip, lpcands, lpcandssol, bdchginds, bdchgdowninfs, nbdchgs) );
            *result = SCIP_REDUCEDDOM;
         }
         freeBdchgs(scip, &bdchginds, &bdchgdowninfs, &nbdchgs);
      }

      /* free buffer for the unreliable candidates */
      SCIPfreeBufferArray(scip, &initcandscores);
      SCIPfreeBufferArray(scip, &initcands);
   }

   /* if no domain could be reduced, create the branching */
   if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM && *result != SCIP_CONSADDED )
   {
      NODE* node;
      Real proveddown;
      Real provedup;
      Real downprio;

      assert(*result == SCIP_DIDNOTRUN);
      assert(0 <= bestcand && bestcand < nlpcands);
      assert(!SCIPisFeasIntegral(scip, lpcandssol[bestcand]));
      assert(SCIPisLT(scip, provedbound, cutoffbound));

      /* perform the branching */
      debugMessage(" -> %d (%d) cands, sel cand %d: var <%s> (sol=%g, down=%g (%+g), up=%g (%+g), sb=%d, psc=%g/%g [%g])\n",
         nlpcands, ninitcands, bestcand, SCIPvarGetName(lpcands[bestcand]), lpcandssol[bestcand],
         bestsbdown, bestsbdown - lpobjval, bestsbup, bestsbup - lpobjval, bestisstrongbranch,
         SCIPgetVarPseudocostCurrentRun(scip, lpcands[bestcand], 
            SCIPfeasFloor(scip, lpcandssol[bestcand]) - lpcandssol[bestcand]),
         SCIPgetVarPseudocostCurrentRun(scip, lpcands[bestcand], 
            SCIPfeasCeil(scip, lpcandssol[bestcand]) - lpcandssol[bestcand]),
         SCIPgetVarPseudocostScoreCurrentRun(scip, lpcands[bestcand], lpcandssol[bestcand]));

      /* choose preferred branching direction */
      switch( SCIPvarGetBranchDirection(lpcands[bestcand]) )
      {
      case SCIP_BRANCHDIR_DOWNWARDS:
         downprio = 1.0;
         break;
      case SCIP_BRANCHDIR_UPWARDS:
         downprio = -1.0;
         break;
      case SCIP_BRANCHDIR_AUTO:
         downprio = SCIPvarGetRootSol(lpcands[bestcand]) - lpcandssol[bestcand];
         break;
      default:
         errorMessage("invalid preferred branching direction <%d> of variable <%s>\n", 
            SCIPvarGetBranchDirection(lpcands[bestcand]), SCIPvarGetName(lpcands[bestcand]));
         return SCIP_INVALIDDATA;
      }

      /* calculate proved lower bounds for children */
      proveddown = provedbound;
      provedup = provedbound;
      if( bestisstrongbranch )
      {
         if( bestsbdownvalid )
            proveddown = MAX(provedbound, bestsbdown);
         if( bestsbupvalid )
            provedup = MAX(provedbound, bestsbup);
      }
      
      /* create child node with x <= floor(x') */
      debugMessage(" -> creating child: <%s> <= %g\n",
         SCIPvarGetName(lpcands[bestcand]), SCIPfeasFloor(scip, lpcandssol[bestcand]));
      CHECK_OKAY( SCIPcreateChild(scip, &node, downprio) );
      CHECK_OKAY( SCIPchgVarUbNode(scip, node, lpcands[bestcand], SCIPfeasFloor(scip, lpcandssol[bestcand])) );
      if( allcolsinlp && !exactsolve )
      {
         assert(SCIPisLT(scip, proveddown, cutoffbound));
         CHECK_OKAY( SCIPupdateNodeLowerbound(scip, node, proveddown) );
      }
      debugMessage(" -> child's lowerbound: %g\n", SCIPnodeGetLowerbound(node));
      
      /* create child node with x >= ceil(x') */
      debugMessage(" -> creating child: <%s> >= %g\n", 
         SCIPvarGetName(lpcands[bestcand]), SCIPfeasCeil(scip, lpcandssol[bestcand]));
      CHECK_OKAY( SCIPcreateChild(scip, &node, -downprio) );
      CHECK_OKAY( SCIPchgVarLbNode(scip, node, lpcands[bestcand], SCIPfeasCeil(scip, lpcandssol[bestcand])) );
      if( allcolsinlp && !exactsolve )
      {
         assert(SCIPisLT(scip, provedup, cutoffbound));
         CHECK_OKAY( SCIPupdateNodeLowerbound(scip, node, provedup) );
      }
      debugMessage(" -> child's lowerbound: %g\n", SCIPnodeGetLowerbound(node));

      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}


/** branching execution method for not completely fixed pseudo solutions */
#define branchExecpsRelpscost NULL




/*
 * branching specific interface methods
 */

/** creates the reliable pseudo cost braching rule and includes it in SCIP */
RETCODE SCIPincludeBranchruleRelpscost(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   BRANCHRULEDATA* branchruledata;

   /* create relpscost branching rule data */
   CHECK_OKAY( SCIPallocMemory(scip, &branchruledata) );
   
   /* include branching rule */
   CHECK_OKAY( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, 
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST,
         branchFreeRelpscost, branchInitRelpscost, branchExitRelpscost, branchInitsolRelpscost, branchExitsolRelpscost, 
         branchExeclpRelpscost, branchExecpsRelpscost,
         branchruledata) );

   /* relpscost branching rule parameters */
   CHECK_OKAY( SCIPaddRealParam(scip,
         "branching/relpscost/minreliable", 
         "minimal value for minimum pseudo cost size to regard pseudo cost value as reliable",
         &branchruledata->minreliable, DEFAULT_MINRELIABLE, 0.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
         "branching/relpscost/maxreliable", 
         "maximal value for minimum pseudo cost size to regard pseudo cost value as reliable",
         &branchruledata->maxreliable, DEFAULT_MAXRELIABLE, 0.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
         "branching/relpscost/sbiterquot", 
         "maximal fraction of strong branching LP iterations compared to node relaxation LP iterations",
         &branchruledata->sbiterquot, DEFAULT_SBITERQUOT, 0.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "branching/relpscost/sbiterofs", 
         "additional number of allowed strong branching LP iterations",
         &branchruledata->sbiterofs, DEFAULT_SBITEROFS, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "branching/relpscost/maxlookahead", 
         "maximal number of further variables evaluated without better score",
         &branchruledata->maxlookahead, DEFAULT_MAXLOOKAHEAD, 1, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "branching/relpscost/initcand", 
         "maximal number of candidates initialized with strong branching per node",
         &branchruledata->initcand, DEFAULT_INITCAND, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "branching/relpscost/inititer", 
         "iteration limit for strong branching initializations of pseudo cost entries (0: auto)",
         &branchruledata->inititer, DEFAULT_INITITER, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "branching/relpscost/maxbdchgs", 
         "maximal number of bound tightenings before the node is reevaluated (-1: unlimited)",
         &branchruledata->maxbdchgs, DEFAULT_MAXBDCHGS, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

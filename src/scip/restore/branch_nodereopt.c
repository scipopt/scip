/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_nodereopt.c
 * @brief  nodereopt branching rule
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>
#include <string.h>
#include <math.h>

#include "scip/branch_nodereopt.h"
#include "scip/branch_relpscost.h"
#include "scip/cons_logicor.h"
#include "scip/scip.h"
#include "scip/tree.h"
#include "scip/reopt.h"

#define BRANCHRULE_NAME            "nodereopt"
#define BRANCHRULE_DESC            "branching rule for node reoptimization"
#define BRANCHRULE_PRIORITY        536870911
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

/*
 * Data structures
 */

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Bool             initialized;             /**< is the data structure initialized? */
   SCIP_Bool             reopt;                   /**< is reoptimization enabled?  */
   SCIP_Bool             strongbranchinginit;     /**< run a strong branching initialization? */
   SCIP_Real             simsolverootlp;          /**< threshold to skip solving the root LP */

   /** Statistic stuff */
   int                   nsplits;                 /**< number of nodes split by the branching rule */
   int                   nrevivednodes;           /**< number of nodes reoptimized by the branching rule */
   SCIP_Real             time;
};

/*
 *  static methods
 */

#ifdef SCIP_DISABLED_CODE
/**
 * add a infeasible constraint, because the problem is proven to be infeasible
 */
static
SCIP_RETCODE addInfeasibleConstraint(
   SCIP*                 scip
)
{
   SCIP_CONS* cons;

   assert(scip != NULL);

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "infeasibleproblem", 0, NULL, NULL, 1, 0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}
#endif

/*
 * heuristic that calculates a largest induced subtree (LR)
 */
#ifdef SCIP_DISABLED_CODE
static
SCIP_RETCODE findLR(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_VAR**            varsLR1,
   SCIP_VAR**            varsLR2,
   SCIP_Real*            valsLR1,
   SCIP_Real*            valsLR2,
   SCIP_BOUNDTYPE*       boundsLR1,
   SCIP_BOUNDTYPE*       boundsLR2,
   int                   nallocs,
   int*                  nvarsLR1,
   int*                  nvarsLR2,
   SCIP_Real*            lossLR1,
   SCIP_Real*            lossLR2,
   int*                  SetToCompress,
   int                   nnodesToCompress,
   int*                  LR1,
   int*                  LR2,
   int*                  nnodesLR1,
   int*                  nnodesLR2
)
{
   SCIP_Bool* commonvars1;
   SCIP_Bool* commonvars2;
   SCIP_Bool* commonvars1_best;
   SCIP_Bool* commonvars2_best;
   SCIP_Bool* nodevars;
   const char** varnames;
   int* nvars;
   int node;
   int parentID;
   int var;

   SCIP_Real tmp_lossLR1;
   SCIP_Real tmp_lossLR2;
   SCIP_Real tmp_minloss;
   SCIP_Real minloss;
   int startnode;

   int tmp_nnodesLR1;
   int tmp_nnodesLR2;
   int tmp_nvarsLR1;
   int tmp_nvarsLR2;

   assert(scip != NULL);
   assert(branchruledata != NULL);

   /* start time */
   SCIP_CALL( SCIPstartClock(scip, branchruledata->lrtime) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &commonvars1, 2*SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &commonvars2, 2*SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &commonvars1_best, 2*SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &commonvars2_best, 2*SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nodevars, 2*SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &varnames, SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nvars, nnodesToCompress) );

   *nnodesLR1 = 0;
   *nnodesLR2 = 0;
   *nvarsLR1 = 0;
   *nvarsLR2 = 0;
   *lossLR1 = SCIPinfinity(scip);
   *lossLR2 = SCIPinfinity(scip);
   minloss = SCIPinfinity(scip);

   for(startnode = 0; startnode < nnodesToCompress; startnode++)
   {
      tmp_nnodesLR1 = 0;
      tmp_nnodesLR2 = 0;
      tmp_nvarsLR1 = 0;
      tmp_nvarsLR2 = 0;

      /*
       * clear bool arrays
       * the array has structure: [x1=0, x1=1, x2=0, x2=1, ... , xn=0, xn=1]
       * */
      for(var = 0; var < 2*SCIPgetNVars(scip); var++)
      {
         commonvars1[var] = FALSE;
         commonvars2[var] = FALSE;
      }

      for(node = 0; node < nnodesToCompress; node++)
         nvars[node] = 0;

      /* initialize array with the first node */
      parentID = SetToCompress[startnode];
      assert(branchruledata->nodedata[parentID] != NULL);
      assert(branchruledata->nodedata[parentID]->nvars >= 1);

      while( parentID != 0 )
      {
         assert( branchruledata->nodedata[parentID] != NULL);
         for(var = 0; var < branchruledata->nodedata[parentID]->nvars; var++)
         {
            if( SCIPisFeasEQ(scip, branchruledata->nodedata[parentID]->varbounds[var], 0))
            {
               commonvars1[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])] = TRUE;
               (tmp_nvarsLR1)++;
            }
            else
            {
               commonvars1[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])+1] = TRUE;
               (tmp_nvarsLR1)++;
            }
            varnames[SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])] = SCIPvarGetName(branchruledata->nodedata[parentID]->vars[var]);
         }

         /* get strong branched vars */
         if( branchruledata->nodedata[parentID]->dualreds )
         {
            LOGICORDATA* consdata;
            int nvarscons;

            assert(branchruledata->nodedata[parentID]->dualconscur != NULL);
            assert(branchruledata->nodedata[parentID]->dualconscur->nvars > 0);

            for(var = 0; var < branchruledata->nodedata[parentID]->dualconscur->nvars; var++)
            {
               if( SCIPisFeasEQ(scip, branchruledata->nodedata[parentID]->dualconscur->vals[var], 0) )
               {
                  commonvars1[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->dualconscur->vars[var])] = TRUE;
                  (tmp_nvarsLR1)++;
               }
               else
               {
                  commonvars1[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->dualconscur->vars[var])+1] = TRUE;
                  (tmp_nvarsLR1)++;
               }
               varnames[SCIPvarGetIndex(branchruledata->nodedata[parentID]->dualconscur->vars[var])] = SCIPvarGetName(branchruledata->nodedata[parentID]->dualconscur->vars[var]);
            }
         }

         /* go to the next saved node above */
         parentID = branchruledata->nodedata[parentID]->parentID;
      }

      /* add the first node to LR1 */
      LR1[tmp_nnodesLR1] = SetToCompress[startnode];
      (tmp_nnodesLR1)++;

      /* set length of root path */
      nvars[startnode] = (tmp_nvarsLR1);

      tmp_lossLR1 = nvars[startnode];
      tmp_lossLR2 = nvars[startnode];

      /* iterate over all other nodes */
      for(node = 0; node < nnodesToCompress; node++)
      {
         if( node != startnode )
         {
            int ncommonvars;

            parentID = SetToCompress[node];
            ncommonvars = 0;

            /* set length of root path */
            nvars[node] = 0;

            /* clear bool array */
            for(var = 0; var < 2*SCIPgetNVars(scip); var++)
               nodevars[var] = FALSE;

            /* check if the intersection is empty and delete variables which are set to a different bound in the current node */
            while( parentID != 0 )
            {
               assert( branchruledata->nodedata[parentID] != NULL);
               for(var = 0; var < branchruledata->nodedata[parentID]->nvars; var++)
               {
                  assert(branchruledata->nodedata[parentID]->vars != NULL);
                  assert(branchruledata->nodedata[parentID]->varbounds != NULL);

                  if( !nodevars[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])]
                   && !nodevars[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])+1] )
                  {
                     if( SCIPisFeasEQ(scip, branchruledata->nodedata[parentID]->varbounds[var], 0) )
                     {
                        /* add the variable to the var array */
                        nodevars[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])] = TRUE;

                        if( commonvars1[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])] )
                           ncommonvars++;
                     }
                     else
                     {
                        /* add the variable to the var array */
                        nodevars[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])+1] = TRUE;

                        if( commonvars1[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])+1] )
                           ncommonvars++;
                     }

                     /* add variables name that are not part of the branching path of the first seen node */
                     if( varnames[SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])] == NULL )
                        varnames[SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])] = SCIPvarGetName(branchruledata->nodedata[parentID]->vars[var]);

                     /* set length of root path */
                     nvars[node]++;
                  }
               }

               /* get strong branched vars */
               if( branchruledata->nodedata[parentID]->dualreds )
               {
                  LOGICORDATA* consdata;
                  int nvarscons;

                  assert(branchruledata->nodedata[parentID]->dualconscur != NULL);
                  assert(branchruledata->nodedata[parentID]->dualconscur->nvars > 0);

                  for(var = 0; var < branchruledata->nodedata[parentID]->dualconscur->nvars; var++)
                  {
                     if( SCIPisFeasEQ(scip, branchruledata->nodedata[parentID]->dualconscur->vals[var], 0) )
                     {
                        /* add the variable to the var array */
                        nodevars[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->dualconscur->vars[var])] = TRUE;

                        if( commonvars1[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->dualconscur->vars[var])] )
                           ncommonvars++;
                     }
                     else
                     {
                        /* add the variable to the var array */
                        nodevars[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->dualconscur->vars[var])+1] = TRUE;

                        if( commonvars1[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->dualconscur->vars[var])+1] )
                           ncommonvars++;
                     }
                     if( varnames[SCIPvarGetIndex(branchruledata->nodedata[parentID]->dualconscur->vars[var])] == NULL )
                        varnames[SCIPvarGetIndex(branchruledata->nodedata[parentID]->dualconscur->vars[var])] = SCIPvarGetName(branchruledata->nodedata[parentID]->dualconscur->vars[var]);

                     /* set length of root path */
                     nvars[node]++;
                  }
               }

               /* go to the next saved node above */
               parentID = branchruledata->nodedata[parentID]->parentID;
            }

            /* delete all variable in commonvars1 and commonvars2 that are not part of vars, respectively. */
            for(var = 0; var < 2*SCIPgetNVars(scip); var++)
            {
               if( ncommonvars > 0 && commonvars1[var] && !nodevars[var] )
               {
                  commonvars1[var] = FALSE;
                  tmp_nvarsLR1--;
                  assert(ncommonvars <= tmp_nvarsLR1);
               }
               /* this is the first node in LR2 */
               else if( ncommonvars == 0 && tmp_nnodesLR2 == 0 && nodevars[var])
               {
                  commonvars2[var] = TRUE;
                  (tmp_nvarsLR2)++;
               }
               else if( ncommonvars == 0 && commonvars2[var] && !nodevars[var] )
               {
                  commonvars2[var] = FALSE;
                  tmp_nvarsLR2--;
                  assert(tmp_nnodesLR2 >= 0);
               }
            }

            assert(ncommonvars == 0 || ncommonvars == tmp_nvarsLR1);

            /* the intersection is empty */
            if( ncommonvars == 0 )
            {
               LR2[tmp_nnodesLR2] = SetToCompress[node];
               (tmp_nnodesLR2)++;
            }
            else
            {
               LR1[tmp_nnodesLR1] = SetToCompress[node];
               (tmp_nnodesLR1)++;
            }

            tmp_lossLR1 += nvars[node];
            tmp_lossLR2 += nvars[node];
         }
      }

      assert(tmp_nnodesLR1 + tmp_nnodesLR2 == nnodesToCompress);

      /* calculate loss of LR1 and LR2 */
      for(node = 0; node < tmp_nnodesLR1; node++)
         (tmp_lossLR1) -= (tmp_nvarsLR1);

      for(node = 0; node < *nnodesLR2; node++)
         (tmp_lossLR2) -= (tmp_nvarsLR2);

      if( tmp_nnodesLR1 > 0 && tmp_nnodesLR2 > 0 )
         tmp_minloss = MIN(tmp_lossLR1/tmp_nnodesLR1, tmp_lossLR2/tmp_nnodesLR2);
      else if( tmp_nnodesLR1 > 0 )
         tmp_minloss = tmp_lossLR1/tmp_nnodesLR1;
      else
         tmp_minloss = tmp_lossLR2/tmp_nnodesLR2;

      if( tmp_minloss < minloss )
      {
         *lossLR1 = tmp_lossLR1;
         *lossLR2 = tmp_lossLR2;

         *nnodesLR1 = tmp_nnodesLR1;
         *nnodesLR2 = tmp_nnodesLR2;

         /* update minloss */
         minloss = tmp_minloss;

         /* collect data for LR1 */
         *nvarsLR1 = 0;
         *nvarsLR2 = 0;
         for(var = 0; var < 2*SCIPgetNVars(scip); var++)
         {
            if( commonvars1[var] )
            {
               varsLR1[*nvarsLR1] = SCIPfindVar(scip, varnames[(int)((SCIP_Real)var/2)]);
               valsLR1[*nvarsLR1] =  var%2 == 0 ? 0 : 1;
               boundsLR1[*nvarsLR1] = (SCIP_BOUNDTYPE) (SCIPisFeasEQ(scip, valsLR1[*nvarsLR1], 0) ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER);

               assert(varsLR1[*nvarsLR1] != NULL);
               (*nvarsLR1)++;
            }

            if( varsLR2 != NULL && commonvars2[var] )
            {
               varsLR2[*nvarsLR2] = SCIPfindVar(scip, varnames[(int)((SCIP_Real)var/2)]);
               valsLR2[*nvarsLR2] =  var%2 == 0 ? 0 : 1;
               boundsLR2[*nvarsLR2] = (SCIP_BOUNDTYPE) (SCIPisFeasEQ(scip, valsLR2[*nvarsLR2], 0) ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER);

               assert(varsLR2[*nvarsLR2] != NULL);
               (*nvarsLR2)++;
            }
         }
      }
   }

   /* free memory */
   SCIPfreeBlockMemoryArray(scip, &commonvars1, 2*SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &commonvars2, 2*SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &nodevars, 2*SCIPgetNVars(scip));
   SCIPfreeMemoryArray(scip, &varnames);
   SCIPfreeBlockMemoryArray(scip, &nvars, nnodesToCompress);

   branchruledata->lrcalls++;

   /* stop time */
   SCIP_CALL( SCIPstopClock(scip, branchruledata->lrtime) );

   return SCIP_OKAY;
}

/*
 * generate weak compression
 */
static
SCIP_RETCODE genWC(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_VAR**            vars,
   SCIP_Real*            vals,
   SCIP_BOUNDTYPE*       bounds,
   LOGICORDATA**         conss,
   int                   nvars,
   int                   nconss
)
{
   LOGICORDATA* consdata;
   int nodeIDtobranch;
   int nodeIDleaf;
   int var;
   int c;
   int parentID;

   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(bounds != NULL);
   assert(conss != NULL);
   assert(nvars > 0);
   assert(nconss >= 0);

   /*
    * fix variables and create node with a corresponding logic-or constraint
    */
   assert(!SCIPqueueIsEmpty(branchruledata->openIDs));

   nodeIDtobranch = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
   assert(branchruledata->nodedata[nodeIDtobranch] == NULL || branchruledata->nodedata[nodeIDtobranch]->nvars == 0);

   /* initialize node */
   SCIP_CALL( initNode(scip, branchruledata, nodeIDtobranch) );
   branchruledata->nodedata[nodeIDtobranch]->reopttype = SCIP_REOPTTYPE_TRANSIT;
   branchruledata->nodedata[nodeIDtobranch]->parentID = 0;

   /* allocate memory for branching decisions */
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &branchruledata->nodedata[nodeIDtobranch]->vars, vars, nvars) );
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &branchruledata->nodedata[nodeIDtobranch]->varbounds, vals, nvars) );
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &branchruledata->nodedata[nodeIDtobranch]->varboundtypes, bounds, nvars) );
   branchruledata->nodedata[nodeIDtobranch]->allocvarmem = nvars;
   branchruledata->nodedata[nodeIDtobranch]->nvars = nvars;

   /* allocate queue for child nodes */
   SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[nodeIDtobranch]->nodechilds, 2, 2) );

   nodeIDleaf = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
   assert(branchruledata->nodedata[nodeIDleaf] == NULL || branchruledata->nodedata[nodeIDleaf]->nvars == 0);

   /* initialize node */
   SCIP_CALL( initNode(scip, branchruledata, nodeIDleaf) );
   branchruledata->nodedata[nodeIDleaf]->reopttype = nvars > 1 ? SCIP_REOPTTYPE_LOGICORNODE : SCIP_REOPTTYPE_LEAF;
   branchruledata->nodedata[nodeIDleaf]->parentID = 0;

   if( nvars > 1 )
   {
      /* create queue for local constaints */
      SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[nodeIDleaf]->conss, 1, 1) );

      /* create and add the constraint */
      SCIP_CALL( SCIPallocMemory(scip, &consdata) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->vars, vars, nvars) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->vals, vals, nvars) );
      consdata->nvars = nvars;
      consdata->constype = REOPT_CONSTYPE_STRBRANCHED;
      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[nodeIDleaf]->conss, (void*) consdata) );
   }
   else
   {
      assert(nvars == 1);

      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeIDleaf]->vars, 1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeIDleaf]->varbounds, 1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeIDleaf]->varboundtypes, 1) );
      branchruledata->nodedata[nodeIDleaf]->allocvarmem = 1;
      branchruledata->nodedata[nodeIDleaf]->nvars = 1;

      branchruledata->nodedata[nodeIDleaf]->vars[0] = vars[0];
      branchruledata->nodedata[nodeIDleaf]->varbounds[0] = 1 - vals[0];
      branchruledata->nodedata[nodeIDleaf]->varboundtypes[0] = (SCIP_BOUNDTYPE) ( SCIP_BOUNDTYPE_UPPER - bounds[0]);
   }

   /* add both nodes as child nodes to the root node */
   if( branchruledata->nodedata[0]->nodechilds == NULL )
   {
      SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[0]->nodechilds, 2, 2) );
   }
   SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[0]->nodechilds, (void*) (size_t) nodeIDtobranch) );
   SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[0]->nodechilds, (void*) (size_t) nodeIDleaf) );

   parentID = nodeIDtobranch;

   /*
    * create child nodes with constraints and corresponding fixings
    */
   for(c = 0; c < nconss; c++)
   {
      LOGICORDATA* consdatabranch;

      /*
       * create node fixings (will be leaf node)
       */
      nodeIDleaf = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
      assert(branchruledata->nodedata[nodeIDleaf] == NULL || branchruledata->nodedata[nodeIDleaf]->nvars == 0);

      /* initialize the node */
      SCIP_CALL( initNode(scip, branchruledata, nodeIDleaf) );
      branchruledata->nodedata[nodeIDleaf]->reopttype = SCIP_REOPTTYPE_LEAF;
      branchruledata->nodedata[nodeIDleaf]->parentID = parentID;

      /* allocate memory for the fixings */
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeIDleaf]->vars, conss[c]->nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeIDleaf]->varbounds, conss[c]->nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeIDleaf]->varboundtypes, conss[c]->nvars) );
      branchruledata->nodedata[nodeIDleaf]->allocvarmem = conss[c]->nvars;
      branchruledata->nodedata[nodeIDleaf]->nvars = conss[c]->nvars;

      /* copy the variables, bound, and boundtypes */
      for(var = 0; var < conss[c]->nvars; var++)
      {
         branchruledata->nodedata[nodeIDleaf]->vars[var] = conss[c]->vars[var];
         branchruledata->nodedata[nodeIDleaf]->varbounds[var] = conss[c]->vals[var];
         branchruledata->nodedata[nodeIDleaf]->varboundtypes[var] = ( SCIPisFeasEQ(scip, conss[c]->vals[var], 1) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER );
      }

      /* add the node as a child of the last branched node */
      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[parentID]->nodechilds, (void*) (size_t) nodeIDleaf) );

      /*
       * create node with constraint (need to be branched)
       */
      nodeIDtobranch = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
      assert(branchruledata->nodedata[nodeIDtobranch] == NULL || branchruledata->nodedata[nodeIDtobranch]->nvars == 0);

      /* initialize node */
      SCIP_CALL( initNode(scip, branchruledata, nodeIDtobranch) );
      if( branchruledata->nodedata[nodeIDtobranch]->nodechilds == NULL )
      {
         SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[nodeIDtobranch]->nodechilds, 2, 2) );
      }
      branchruledata->nodedata[nodeIDtobranch]->reopttype = SCIP_REOPTTYPE_LOGICORNODE;
      branchruledata->nodedata[nodeIDtobranch]->parentID = parentID;

      /* create queue for local constraints */
      SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[nodeIDtobranch]->conss, 1, 1) );

      /* add the constraint */
      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[nodeIDtobranch]->conss, (void*) conss[c]) );

      /* add the node as a child of the last branched node */
      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[parentID]->nodechilds, (void*) (size_t) nodeIDtobranch) );

      /* set the current node as the next node to branch */
      parentID = nodeIDtobranch;
   }

   branchruledata->wclastnodeID = nodeIDtobranch;

   return SCIP_OKAY;
}



/*
 * reduce the tree to feasibility
 */
static
SCIP_RETCODE genLC(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_VAR***           vars,
   SCIP_Real**           bounds,
   SCIP_BOUNDTYPE**      boundtypes,
   int                   nnodes,
   int*                  nvars,
   SCIP_QUEUE**          conss,
   int                   newparentID,
   int*                  nodeID_cons,
   SCIP_Real*            avgdepth
)
{
   int consID;
   int nodeID;
   int var;
   int newID;

   assert(scip != NULL);
   assert(branchruledata != NULL);

   /* start time */
   SCIP_CALL( SCIPstartClock(scip, branchruledata->wctime) );

   branchruledata->wccalls++;
   (*avgdepth) = 0;

   /* ensure that parent node is allocated */
   assert(branchruledata->nodedata[newparentID] != NULL);
   if( branchruledata->nodedata[newparentID]->nodechilds == NULL )
   {
      SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[newparentID]->nodechilds, nnodes+1, 2) );
   }

   /* save node with added constraints */
   if (SCIPqueueIsEmpty(branchruledata->openIDs))
   {
      SCIP_CALL(reallocNodedata(scip, branchruledata));
   }

   newID = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
   assert(newID >= 1);
   assert(newID < branchruledata->allocmemsizenodedata);
   assert(branchruledata->nodedata[newID] == NULL
      || (branchruledata->nodedata[newID]->nvars == 0
          && (branchruledata->nodedata[newID]->conss == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[newID]->conss))) );

   if( branchruledata->nodedata[newID] == NULL )
   {
      SCIP_CALL( initNode(scip, branchruledata, newID) );
   }

   /* set the root as parent node*/
   branchruledata->nodedata[newID]->parentID = newparentID;

   /* set the reopttype */
   branchruledata->nodedata[newID]->reopttype = SCIP_REOPTTYPE_LOGICORNODE;

   *nodeID_cons = newID;
   consID = newID;

   /* increase number of saved nodes */
   branchruledata->nsavednodes++;

   if( branchruledata->nodedata[newID]->conss == NULL )
   {
      SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[newID]->conss, nnodes, 2) );
   }

   /* add the node as a child of the root node */
   assert(branchruledata->nodedata[newparentID]->nodechilds != NULL);
   SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[newparentID]->nodechilds, (void*) (size_t) newID) );

   /** save the new generated nodes */

   /* store nodes */
   for(nodeID = 0; nodeID < nnodes; nodeID++)
   {
      if (SCIPqueueIsEmpty(branchruledata->openIDs))
      {
         SCIP_CALL(reallocNodedata(scip, branchruledata));
      }

      newID = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
      assert(newID >= 1);
      assert(newID < branchruledata->allocmemsizenodedata);
      assert(nvars[nodeID] > 0);
      assert(vars[nodeID] != NULL);
      assert(bounds[nodeID] != NULL);
      assert(boundtypes[nodeID] != NULL);
      assert(branchruledata->nodedata[newID] == NULL
         || (branchruledata->nodedata[newID]->nvars == 0
             && (branchruledata->nodedata[newID]->conss == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[newID]->conss))) );

      if( branchruledata->nodedata[newID] == NULL )
      {
         SCIP_CALL( initNode(scip, branchruledata, newID) );
      }

      /* set the parent node*/
      branchruledata->nodedata[newID]->parentID = newparentID;

      /* set the reopttype */
      branchruledata->nodedata[newID]->reopttype = SCIP_REOPTTYPE_LEAF;

      /* increase number of saved nodes */
      branchruledata->nsavednodes++;

      /* allocate memory */
      if( branchruledata->nodedata[newID]->allocvarmem == 0 )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[newID]->vars, nvars[nodeID]) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[newID]->varbounds, nvars[nodeID]) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[newID]->varboundtypes, nvars[nodeID]) );
         branchruledata->nodedata[newID]->allocvarmem = nvars[nodeID];
      }

      assert(branchruledata->nodedata[newID]->vars != NULL);
      assert(branchruledata->nodedata[newID]->varbounds != NULL);
      assert(branchruledata->nodedata[newID]->varboundtypes != NULL);

      /* copy variable data */
      for(var = 0; var < nvars[nodeID]; var++)
      {
         branchruledata->nodedata[newID]->vars[var] = vars[nodeID][var];
         branchruledata->nodedata[newID]->varbounds[var] = bounds[nodeID][var];
         branchruledata->nodedata[newID]->varboundtypes[var] = boundtypes[nodeID][var];
         branchruledata->nodedata[newID]->nvars++;
      }

      (*avgdepth) += nvars[nodeID];

      /* copy constraint data */
      if( branchruledata->nodedata[newID]->conss == NULL && !SCIPqueueIsEmpty(conss[nodeID]) )
      {
         SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[newID]->conss, SCIPqueueNElems(conss[nodeID])+1, 2) );
      }
      assert(branchruledata->nodedata[newID]->conss == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[newID]->conss));

      while( !SCIPqueueIsEmpty(conss[nodeID]) )
      {
         SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[newID]->conss, SCIPqueueRemove(conss[nodeID])) );
      }

      /* create a constraint and add them to the node at consID */
      if( nvars[nodeID] == 1 )
      {
         /* allocate memory */
         if( branchruledata->nodedata[consID]->allocvarmem == 0 )
         {
            SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[consID]->vars, nvars[nodeID]) );
            SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[consID]->varbounds, nvars[nodeID]) );
            SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[consID]->varboundtypes, nvars[nodeID]) );
            branchruledata->nodedata[consID]->allocvarmem = SCIPgetNOrigBinVars(scip);
         }

         /* fix the variable to the negated values */
         branchruledata->nodedata[consID]->vars[branchruledata->nodedata[consID]->nvars] = vars[nodeID][0];
         branchruledata->nodedata[consID]->varbounds[branchruledata->nodedata[consID]->nvars] = 1 - bounds[nodeID][0];
         branchruledata->nodedata[consID]->varboundtypes[branchruledata->nodedata[consID]->nvars] = ( 1 - boundtypes[nodeID][0] == 0 ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER );
         branchruledata->nodedata[consID]->nvars++;
      }
      else
      {
         LOGICORDATA* consdata;

         SCIP_CALL( SCIPallocMemory(scip, &consdata) );
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->vars, vars[nodeID], nvars[nodeID]) );
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->vals, bounds[nodeID], nvars[nodeID]) );
         consdata->nvars = nvars[nodeID];
         consdata->constype = REOPT_CONSTYPE_INFSUBTREE;

         /* add the constraint to the node */
         SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[consID]->conss, (void*) consdata) );
      }

      /* add the node as a child node */
      assert(branchruledata->nodedata[newparentID]->nodechilds != NULL);
      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[newparentID]->nodechilds, (void*) (size_t) newID) );
   }

   (*avgdepth) /= nnodes;

   /* for each constraint that consists of exactly one variable, we fixed the variable to the negated value,
    * thus, the number of constraints + number of variables has to be nfeasnodeIDs */
   assert(nnodes == SCIPqueueNElems(branchruledata->nodedata[consID]->conss) + branchruledata->nodedata[consID]->nvars);

   /* stop time */
   SCIP_CALL( SCIPstopClock(scip, branchruledata->wctime) );

   return SCIP_OKAY;
}

/*
 * find weak compression
 */
static
SCIP_RETCODE findWC(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nodeID,
   SCIP_VAR**            vars,
   SCIP_Real*            vals,
   SCIP_BOUNDTYPE*       bounds,
   int*                  nvars,
   LOGICORDATA**         conss,
   int*                  nconss
)
{
   int parentID;
   int var;

   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(bounds != NULL);
   assert(nvars != NULL);


   parentID = nodeID;
   *nvars = 0;
   *nconss = 0;

   /* go up to the root */
   while( parentID != 0 )
   {
      assert(branchruledata->nodedata[parentID] != NULL);

      /* copy branching information */
      for(var = 0; var < branchruledata->nodedata[parentID]->nvars; var++)
      {
         vars[*nvars] = branchruledata->nodedata[parentID]->vars[var];
         vals[*nvars] = branchruledata->nodedata[parentID]->varbounds[var];
         bounds[*nvars] = branchruledata->nodedata[parentID]->varboundtypes[var];
         (*nvars) += 1;
      }

      /* collect all added constraints along the root path */
      if( branchruledata->nodedata[parentID]->conss != NULL && !SCIPqueueIsEmpty(branchruledata->nodedata[parentID]->conss) )
      {
         int ncons;
         int cons;

         ncons = SCIPqueueNElems(branchruledata->nodedata[parentID]->conss);
         cons = 0;
         while( cons < ncons )
         {
            LOGICORDATA* consdata;
            LOGICORDATA* consdataCopy;

            consdata = (LOGICORDATA*) SCIPqueueRemove(branchruledata->nodedata[parentID]->conss);

            if( consdata->constype == REOPT_CONSTYPE_STRBRANCHED )
            {
               SCIP_CALL( SCIPallocMemory(scip, &consdataCopy) );
               SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdataCopy->vars, consdata->vars, consdata->nvars) );
               SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdataCopy->vals, consdata->vals, consdata->nvars) );
               consdataCopy->nvars = consdata->nvars;
               consdataCopy->constype = consdata->constype;

               conss[*nconss] = consdataCopy;
               (*nconss) += 1;
            }
            cons++;
            SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[parentID]->conss, (void*) consdata) );
         }
      }

      parentID = branchruledata->nodedata[parentID]->parentID;
   }

   return SCIP_OKAY;
}

/*
 * find a lazy compression
 */
static
SCIP_RETCODE findLC(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int*                  SetToCompress,
   int                   nnodesToCompress,
   SCIP_VAR***           vars,
   SCIP_Real**           vals,
   SCIP_BOUNDTYPE**      bounds,
   int*                  nvars,
   SCIP_QUEUE**          conss
)
{
   int nodeID;
   int parentID;
   int var;

   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(nnodesToCompress >= 0);
   assert(SetToCompress != NULL || nnodesToCompress == 0);
   assert(vars != NULL || nnodesToCompress == 0);
   assert(vals != NULL || nnodesToCompress == 0);
   assert(bounds != NULL || nnodesToCompress == 0);
   assert(nvars != NULL || nnodesToCompress == 0);

   /* collect the root path */
   for(nodeID = 0; nodeID < nnodesToCompress; nodeID++)
   {
      assert(vars[nodeID] != NULL);
      assert(vals[nodeID] != NULL);
      assert(bounds[nodeID] != NULL);

      parentID = SetToCompress[nodeID];
      nvars[nodeID] = 0;

      /* go up to the root */
      while( parentID != 0 )
      {
         assert(branchruledata->nodedata[parentID] != NULL);

         /* copy branching information */
         for(var = 0; var < branchruledata->nodedata[parentID]->nvars; var++)
         {
            vars[nodeID][nvars[nodeID]] = branchruledata->nodedata[parentID]->vars[var];
            vals[nodeID][nvars[nodeID]] = branchruledata->nodedata[parentID]->varbounds[var];
            bounds[nodeID][nvars[nodeID]] = branchruledata->nodedata[parentID]->varboundtypes[var];
            nvars[nodeID]++;
         }

         /* collect all added constaints along the root path */
         if( branchruledata->nodedata[parentID]->conss != NULL && !SCIPqueueIsEmpty(branchruledata->nodedata[parentID]->conss) )
         {
            int ncons;
            int cons;

            ncons = SCIPqueueNElems(branchruledata->nodedata[parentID]->conss);
            cons = 0;
            while( cons < ncons )
            {
               LOGICORDATA* consdata;
               LOGICORDATA* consdataCopy;

               consdata = (LOGICORDATA*) SCIPqueueRemove(branchruledata->nodedata[parentID]->conss);

               if( consdata->constype == REOPT_CONSTYPE_STRBRANCHED )
               {
                  SCIP_CALL( SCIPallocMemory(scip, &consdataCopy) );
                  SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdataCopy->vars, consdata->vars, consdata->nvars) );
                  SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdataCopy->vals, consdata->vals, consdata->nvars) );
                  consdataCopy->nvars = consdata->nvars;
                  consdataCopy->constype = consdata->constype;
               }
               cons++;
               SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[parentID]->conss, (void*) consdata) );
            }
         }

         /* collect all dual branching information */
         if( branchruledata->nodedata[parentID]->dualreds )
         {
            LOGICORDATA* consdata;
            int nvarscons;

            nvarscons = SCIPbranchrulePseudoGetNPseudoVars(scip, parentID);
            SCIP_CALL( SCIPallocMemory(scip, &consdata) );
            SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->vars, nvarscons) );
            SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->vals, nvarscons) );
            consdata->nvars = -1;
            SCIP_CALL( SCIPbranchrulePseudoGenerateCons(scip, consdata, &(consdata->nvars), nvarscons, parentID, FALSE, FALSE) );

            assert(consdata->nvars <= nvarscons);

            for(var = 0; var < consdata->nvars; var++)
            {
               vars[nodeID][nvars[nodeID]] = consdata->vars[var];
               vals[nodeID][nvars[nodeID]] = consdata->vals[var];
               bounds[nodeID][nvars[nodeID]] = ( SCIPisFeasEQ(scip, consdata->vals[var], 0) ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER );
               nvars[nodeID]++;
            }

            SCIPfreeMemoryArray(scip, &consdata->vals);
            SCIPfreeMemoryArray(scip, &consdata->vars);
            SCIPfreeMemory(scip, &consdata);
         }

         parentID = branchruledata->nodedata[parentID]->parentID;
      }

      assert(nvars[nodeID] > 0);
   }

   return SCIP_OKAY;
}

/*
 * generate the LR
 */
static
SCIP_RETCODE genLR(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_VAR**            vars,
   SCIP_Real*            vals,
   SCIP_BOUNDTYPE*       bounds,
   int                   nvars,
   int                   parentID,
   int*                  consnodeID
)
{
   LOGICORDATA* consdata;
   int LRnodeID;
   int notLRnodeID;
   int var;

   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);
   assert(nvars == 0 || bounds != NULL);
   assert(!SCIPqueueIsEmpty(branchruledata->openIDs));

   /* start time */
   SCIP_CALL( SCIPstartClock(scip, branchruledata->lrtime) );

   /** ensure that the parent node is allocated */
   assert(0 <= parentID && parentID < branchruledata->allocmemsizenodedata);
   assert(branchruledata->nodedata[parentID] != NULL);

   /** generate LR child */
   if (SCIPqueueIsEmpty(branchruledata->openIDs))
   {
      SCIP_CALL(reallocNodedata(scip, branchruledata));
   }

   LRnodeID = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
   assert(branchruledata->nodedata[LRnodeID] == NULL
      || (branchruledata->nodedata[LRnodeID]->nvars == 0
          && (branchruledata->nodedata[LRnodeID]->conss == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[LRnodeID]->conss))) );
   assert(LRnodeID > 0);
   assert(branchruledata->nodedata[LRnodeID] == NULL || branchruledata->nodedata[LRnodeID]->nvars == 0);

   if( branchruledata->nodedata[LRnodeID] == NULL )
   {
      SCIP_CALL( initNode(scip, branchruledata, LRnodeID) );
   }

   if(branchruledata->nodedata[LRnodeID]->allocvarmem == 0)
   {
      assert(branchruledata->nodedata[LRnodeID]->vars == NULL );
      assert(branchruledata->nodedata[LRnodeID]->varbounds == NULL );
      assert(branchruledata->nodedata[LRnodeID]->varboundtypes == NULL );

      /** Allocate memory for node information */
      branchruledata->nodedata[LRnodeID]->allocvarmem = SCIPgetNOrigVars(scip);
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[LRnodeID]->vars), branchruledata->nodedata[LRnodeID]->allocvarmem) );
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[LRnodeID]->varbounds), branchruledata->nodedata[LRnodeID]->allocvarmem) );
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[LRnodeID]->varboundtypes), branchruledata->nodedata[LRnodeID]->allocvarmem) );
   }

   /* check memory */
   SCIP_CALL( checkMemory(scip, branchruledata, LRnodeID, nvars) );
   assert(branchruledata->nodedata[LRnodeID]->allocvarmem >= nvars);

   /* copy bounds */
   for(var = 0; var < nvars; var++)
   {
      assert(vars[var] != NULL);
      branchruledata->nodedata[LRnodeID]->vars[var] = vars[var];
      branchruledata->nodedata[LRnodeID]->varbounds[var] = vals[var];
      branchruledata->nodedata[LRnodeID]->varboundtypes[var] = bounds[var];
      branchruledata->nodedata[LRnodeID]->nvars++;
   }

   /* set reopttype */
   branchruledata->nodedata[LRnodeID]->reopttype = SCIP_REOPTTYPE_LEAF;

   /* set parentID */
   branchruledata->nodedata[LRnodeID]->parentID = parentID;

   /* increase number of saved nodes */
   branchruledata->nsavednodes++;

   /** generate child with constraint */
   if (SCIPqueueIsEmpty(branchruledata->openIDs))
   {
      SCIP_CALL(reallocNodedata(scip, branchruledata));
   }

   notLRnodeID = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
   assert(branchruledata->nodedata[notLRnodeID] == NULL
      || (branchruledata->nodedata[notLRnodeID]->nvars == 0
          && (branchruledata->nodedata[notLRnodeID]->conss == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[notLRnodeID]->conss))) );
   assert(branchruledata->nodedata[notLRnodeID] == NULL || branchruledata->nodedata[notLRnodeID]->nvars == 0);
   assert(notLRnodeID > 0);

   if( branchruledata->nodedata[notLRnodeID] == NULL )
   {
      SCIP_CALL( initNode(scip, branchruledata, notLRnodeID) );
   }

   if( nvars > 1 )
   {
      if (branchruledata->nodedata[notLRnodeID]->conss == NULL)
      {
         SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[notLRnodeID]->conss, 2, 2) );
      }

      SCIP_CALL( SCIPallocMemory(scip, &consdata) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->vars, vars, nvars) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->vals, vals, nvars) );
      consdata->nvars = nvars;
      consdata->constype = REOPT_CONSTYPE_INFSUBTREE;

      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[notLRnodeID]->conss, (void*) consdata) );
   }
   else
   {
      if(branchruledata->nodedata[notLRnodeID]->allocvarmem == 0)
      {
         assert(branchruledata->nodedata[notLRnodeID]->vars == NULL );
         assert(branchruledata->nodedata[notLRnodeID]->varbounds == NULL );
         assert(branchruledata->nodedata[notLRnodeID]->varboundtypes == NULL );

         /** Allocate memory for node information */
         branchruledata->nodedata[notLRnodeID]->allocvarmem = SCIPgetNOrigVars(scip);
         SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[notLRnodeID]->vars), branchruledata->nodedata[notLRnodeID]->allocvarmem) );
         SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[notLRnodeID]->varbounds), branchruledata->nodedata[notLRnodeID]->allocvarmem) );
         SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[notLRnodeID]->varboundtypes), branchruledata->nodedata[notLRnodeID]->allocvarmem) );
      }

      /* check memory */
      SCIP_CALL( checkMemory(scip, branchruledata, notLRnodeID, 1) );

      assert(vars[0] != NULL);
      branchruledata->nodedata[notLRnodeID]->vars[0] = vars[0];
      branchruledata->nodedata[notLRnodeID]->varbounds[0] = 1 - vals[0];
      branchruledata->nodedata[notLRnodeID]->varboundtypes[0] = (SCIP_BOUNDTYPE) (1 - bounds[0]);
      branchruledata->nodedata[notLRnodeID]->nvars++;
   }

   /* set reopttype */
   branchruledata->nodedata[notLRnodeID]->reopttype = SCIP_REOPTTYPE_LOGICORNODE;

   /* set parentID */
   branchruledata->nodedata[notLRnodeID]->parentID = parentID;

   /* set the ID of the node with the added constraint */
   (*consnodeID) = notLRnodeID;

   /* add this two nodes as child nodes below the root */
   if( branchruledata->nodedata[parentID]->nodechilds == NULL )
   {
      SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[parentID]->nodechilds, 2, 2) );
   }

   SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[parentID]->nodechilds, (void*) (size_t) LRnodeID) );
   SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[parentID]->nodechilds, (void*) (size_t) notLRnodeID) );

   /* stop time */
   SCIP_CALL( SCIPstopClock(scip, branchruledata->lrtime) );

   return SCIP_OKAY;
}

/*
 * run heuristics to compress the search tree
 */
static
SCIP_RETCODE runHeuristics(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_Bool*            success
)
{
   SCIP_VAR** varsLR1;
   SCIP_VAR** varsLR2;
   SCIP_Real* valsLR1;
   SCIP_Real* valsLR2;
   SCIP_BOUNDTYPE* boundsLR1;
   SCIP_BOUNDTYPE* boundsLR2;
   SCIP_Real lossLR1;
   SCIP_Real lossLR2;
   SCIP_Real minLoss;
   SCIP_Bool successLR1;
   SCIP_Bool successLR2;
   int depth;
   int nnodesLR1;
   int nnodesLR2;
   int nodeID_cons;
   int nodeID;
   int nnodesToCompress;
   int nallocvars;
   int nvarsLR1;
   int nvarsLR2;
   int nrepresentatives;
   int allocsize;
   int* SetToCompress;
   int* LR1;
   int* LR2;
   int nodeWC;
   int fixedvarsWC;

   SCIP_VAR** varsWC;
   SCIP_Real* valsWC;
   SCIP_BOUNDTYPE* boundsWC;
   int nvarsWC;
   LOGICORDATA** conssWC;
   int nconssWC;

   assert(scip != NULL);
   assert(branchruledata != NULL);

   depth = 0;
   nnodesLR1 = 0;
   nnodesLR2 = 0;
   nodeID_cons = 0;
   nodeID = 0;
   nnodesToCompress = 0;
   nallocvars = 0;
   nvarsLR1 = 0;
   nvarsLR2 = 0;
   nrepresentatives = 0;
   allocsize = 0;

   lossLR1 = SCIPinfinity(scip);
   lossLR2 = SCIPinfinity(scip);

   nodeWC = -1;
   nvarsWC = 0;
   nconssWC = 0;
   fixedvarsWC = 1;

   /* allocate general memory */
   allocsize = branchruledata->allocmemsizenodedata - SCIPqueueNElems(branchruledata->openIDs);
   nallocvars = SCIPgetNBinVars(scip);
   SCIP_CALL( SCIPallocMemoryArray(scip, &SetToCompress, allocsize) );

   /* allocate memory for LR*/
   SCIP_CALL( SCIPallocMemoryArray(scip, &varsLR1, nallocvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &valsLR1, nallocvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &boundsLR1, nallocvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &boundsLR2, nallocvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &varsLR2, nallocvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &valsLR2, nallocvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &LR1, allocsize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &LR2, allocsize) );

   /* allocate memory for WC */
   SCIP_CALL( SCIPallocMemoryArray(scip, &varsWC, nallocvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &valsWC, nallocvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &boundsWC, nallocvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &conssWC, 20) );

   /* collect nodeIDs */
   nnodesToCompress= 0;
   for(nodeID = 1; nodeID < branchruledata->allocmemsizenodedata; nodeID++)
   {
      if( branchruledata->nodedata[nodeID] != NULL
       && lengthBranchPathByID(branchruledata, nodeID) >= 1
       && (branchruledata->nodedata[nodeID]->nodechilds == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds)) )
      {
         if( branchruledata->cpressnodes == 0 )
         {
            SetToCompress[nnodesToCompress] = nodeID;
            nnodesToCompress++;
         }
         else if( branchruledata->cpressnodes == 1 && branchruledata->nodedata[nodeID]->reopttype != SCIP_REOPTTYPE_FEASIBLE )
         {
            SetToCompress[nnodesToCompress] = nodeID;
            nnodesToCompress++;
         }
         else if( branchruledata->cpressnodes == 2 && branchruledata->nodedata[nodeID]->reopttype == SCIP_REOPTTYPE_FEASIBLE )
         {
            SetToCompress[nnodesToCompress] = nodeID;
            nnodesToCompress++;
         }
      }
   }

   nrepresentatives = 0;
   nodeID_cons = 0;
   *success = FALSE;
   successLR1 = FALSE;
   successLR2 = FALSE;
   depth = 0;
   minLoss = SCIPinfinity(scip);

   /* find LR nodes */
   if( nnodesToCompress > 0 && branchruledata->lrenable )
   {
      SCIP_CALL( findLR(scip, branchruledata, varsLR1, varsLR2, valsLR1, valsLR2, boundsLR1, boundsLR2,
            nallocvars, &nvarsLR1, &nvarsLR2, &lossLR1, &lossLR2, SetToCompress, nnodesToCompress, LR1, LR2,
            &nnodesLR1, &nnodesLR2) );

      minLoss = MIN(lossLR1/nnodesLR1, lossLR2/nnodesLR2);
      depth = lossLR1/nnodesLR1 <= lossLR2/nnodesLR2 ? nvarsLR1 : nvarsLR2;

      if( branchruledata->lrenable && (minLoss < branchruledata->lrloss || depth > branchruledata->lrdepth) )
      {
         if( nnodesLR1 >= 1 && nvarsLR1 > 1 ) /* if nnnodesLR1feas = 1 then LR is equivalent to LC */
            successLR1 = TRUE;
         if( nnodesLR2 >= 1 && nvarsLR2 > 1 ) /* if nnnodesLR1feas = 1 then LR is equivalent to LC */
            successLR2 = TRUE;
      }
   }

   if( !successLR1 && !successLR2 && !branchruledata->wcenable )
      goto SKIP;

   assert(!branchruledata->lrenable || nnodesToCompress == nnodesLR1 + nnodesLR2);

   if( !successLR1 && !successLR2 && branchruledata->lastprunedID > 0 && branchruledata->wcenable)
   {
      int nfixedvars;
      nfixedvars = lengthBranchPathByID(branchruledata, branchruledata->lastprunedID);

      branchruledata->wccalls++;

      if( nfixedvars > branchruledata->wcdepth && nfixedvars < SCIPgetNBinVars(scip) )
      {
         fixedvarsWC = nfixedvars;
         nodeWC = branchruledata->lastprunedID;
      }
   }

   if( nodeWC == branchruledata->wclastnodeID )
      nodeWC = -1;


   /* find weak compression of nodes */
   if( !successLR1 && !successLR2 && branchruledata->wcenable && nodeWC > 0 )
   {
      SCIP_CALL( findWC(scip, branchruledata, nodeWC, varsWC, valsWC, boundsWC, &nvarsWC, conssWC, &nconssWC) );
   }

   assert(successLR1 + (nodeWC > 0) <= 1);
   assert(successLR2 + (nodeWC > 0) <= 1);

   if( !successLR1 && !successLR2 && nodeWC <= 0 )
      goto SKIP;

   /* if at least one heuristic was successful we clear the node data */
   if( successLR1 || successLR2 || nodeWC > 0 )
   {
      /* reset the saved data */
      SCIP_CALL( clearNodes(scip, branchruledata, FALSE) );
      SCIP_CALL( SCIPbranchrulePseudoReset(scip, TRUE, FALSE) );

      /* initialize the root data */
      SCIP_CALL( initNode(scip, branchruledata, 0) );
   }

   /* compress the nodes */
   if( successLR1 || successLR2 || nodeWC > 0 )
   {
      if( successLR1 || successLR2 )
      {
         assert(branchruledata->lrenable);

         /* do the compression */
         if( successLR1 && lossLR1 <= lossLR2 )
         {
            /* generate LR1 */
            assert(branchruledata->nodedata[0] != NULL);
            assert(successLR1);

            SCIP_CALL( genLR(scip, branchruledata, varsLR1, valsLR1, boundsLR1, nvarsLR1, 0, &nodeID_cons) );

            assert(0 < nodeID_cons && nodeID_cons < branchruledata->allocmemsizenodedata);
            nrepresentatives = 2;

            /* generate LR2 */
            if( successLR2 )
            {
               assert(branchruledata->nodedata[nodeID_cons] != NULL);
               assert(successLR2);

               SCIP_CALL( genLR(scip, branchruledata, varsLR2, valsLR2, boundsLR2, nvarsLR2, nodeID_cons, &nodeID_cons) );

               assert(0 < nodeID_cons && nodeID_cons < branchruledata->allocmemsizenodedata);
               nrepresentatives++;
            }

            branchruledata->lrsuccess++;
            branchruledata->lrloss = minLoss;
            branchruledata->lrdepth = MAX(branchruledata->lrdepth, depth);

            *success = TRUE;
         }
         else if( successLR2 && lossLR1 > lossLR2 )
         {
            /* generate LR2 */
            assert(branchruledata->nodedata[0] != NULL);
            SCIP_CALL( genLR(scip, branchruledata, varsLR2, valsLR2, boundsLR2, nvarsLR2, nodeID_cons, &nodeID_cons) );
            assert(0 < nodeID_cons && nodeID_cons < branchruledata->allocmemsizenodedata);
            nrepresentatives = 2;

            /* generate LR1 if nnodeLR1 > 0 */
            if( successLR2 )
            {
               assert(branchruledata->nodedata[nodeID_cons] != NULL);
               SCIP_CALL( genLR(scip, branchruledata, varsLR1, valsLR1, boundsLR1, nvarsLR1, nodeID_cons, &nodeID_cons) );
               assert(0 < nodeID_cons && nodeID_cons < branchruledata->allocmemsizenodedata);
               nrepresentatives++;
            }

            branchruledata->lrsuccess++;
            branchruledata->lrloss = minLoss;
            branchruledata->lrdepth = MAX(branchruledata->lrdepth, depth);

            *success = TRUE;
         }

         branchruledata->lrk += nrepresentatives;
      }
      else
      {
         if( nodeWC > 0)
         {
            SCIP_CALL( genWC(scip, branchruledata, varsWC, valsWC, boundsWC, conssWC, nvarsWC, nconssWC) );

            branchruledata->wcsuccess++;
            nrepresentatives = nconssWC+2;
            branchruledata->wck += nrepresentatives;
            branchruledata->wcdepth = nvarsWC + 1 >= 0.75*SCIPgetNBinVars(scip) ? 1 : nvarsWC;
            *success = TRUE;
         }
      }
   }

   if( *success )
   {
      printf("** reoptimization ** heuristic compression of the search frontier:\n");
      printf("*               nnodes       loss loss/nodes      depth\n");

      printf("* nodes    %10d %10s %10s %10s\n", nnodesToCompress, "", "", "");

      if( successLR1 )
         printf("*    LR1 : %10d %10.2f %10.2f %10d\n", nnodesLR1, lossLR1, lossLR1/nnodesLR1, nvarsLR1);
      else
         printf("*    LR1 : %10d %10s %10s %10s\n", 0, "-", "-", "-");

      if( successLR2 )
         printf("*    LR2 : %10d %10.2f %10.2f %10d\n", nnodesLR2, lossLR2, lossLR2/nnodesLR2, nvarsLR2);
      else
         printf("*    LR2 : %10d %10s %10s %10s\n", 0, "-", "-", "-");

      if( nodeWC > 0 )
         printf("*     WC : %10d %10s %10s %10d (nvars)\n", 1, "-", "-", nvarsWC);
      else
         printf("*     WC : %10d %10s %10s %10s\n", 0, "-", "-", "-");

      printf("* \n");
      printf("* compression of size: %5d\n", nrepresentatives);
      printf("************************************************************\n");
   }

   SKIP:

   /* allocate memory for WC */
   SCIPfreeMemoryArray(scip, &varsWC);
   SCIPfreeMemoryArray(scip, &valsWC);
   SCIPfreeMemoryArray(scip, &boundsWC);
   SCIPfreeMemoryArray(scip, &conssWC);

   /* free memory for LR */
   SCIPfreeMemoryArray(scip, &LR2);
   SCIPfreeMemoryArray(scip, &LR1);
   SCIPfreeMemoryArray(scip, &valsLR2);
   SCIPfreeMemoryArray(scip, &varsLR2);
   SCIPfreeMemoryArray(scip, &boundsLR2);
   SCIPfreeMemoryArray(scip, &boundsLR1);
   SCIPfreeMemoryArray(scip, &valsLR1);
   SCIPfreeMemoryArray(scip, &varsLR1);

   /* free general memory */
   SCIPfreeMemoryArray(scip, &SetToCompress);

   return SCIP_OKAY;
}
#endif

/*
 * Execute the branching of nodes with additional constraints.
 */
static
SCIP_RETCODE Exec(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_RESULT*          result
)
{
   SCIP_Bool localrestart;
   SCIP_NODE* node;
   SCIP_CONS* cons;

   int* childids;
   int nchilds;
   int nodeID;
   int childID;
   int c;
   int ncreatedchilds;

   assert(scip != NULL );
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);

   node = SCIPgetCurrentNode(scip);

   nodeID = SCIPgetRootNode(scip) == node ? 0 : SCIPnodeGetReoptID(node);

   /* calculate local similarity and delete the induced subtree if
    * the similarity is to low */
   localrestart = FALSE;
   SCIP_CALL( SCIPcheckLocalRestart(scip, node, &localrestart) );

   if( localrestart )
   {
      *result = SCIP_DIDNOTRUN;
      goto TERMINATE;
   }

   SCIPdebugMessage("current node is %lld, ID %d:\n", SCIPnodeGetNumber(node), nodeID);

   /**
    * current node is equal to the root and the root was pseudo-branched
    * we have to create two child nodes; one with the pseudo-constraint and
    * one with the negated fixings.
    */
   if(SCIPgetRootNode(scip) == node && SCIPnodeSplit(scip, node) )
   {
      LOGICORDATA* consdata;
      int v;

      consdata = NULL;

      /* allocate buffer memory */
      SCIP_CALL( SCIPallocBuffer(scip, &consdata) );
      SCIP_CALL( SCIPallocBufferArray(scip, &consdata->vars, SCIPgetNBinVars(scip)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &consdata->vals, SCIPgetNBinVars(scip)) );
      consdata->allocmem = SCIPgetNBinVars(scip);
      consdata->nvars = 0;

      /* get the constraint which splits the node */
      SCIP_CALL( SCIPnodeGetSplitCons(scip, nodeID, consdata) );

      /**
       * if the bound changes based on dual information induces an
       * infeasible subtree, we add this constraint as a global one.
       * otherwise, we split the root into two dummy nodes and
       * proceed as like the root is a transit node.
       */
      switch (consdata->constype) {
         case REOPT_CONSTYPE_STRBRANCHED:

            /* split the root into two dummy nodes */
            SCIP_CALL( SCIPsplitReoptRoot(scip) );

         break;

         case REOPT_CONSTYPE_INFSUBTREE:
            assert( SCIPgetNReoptNodeIDs(scip, node) == 0 );

            /* the logic-or constraint induces an infeasible subtree
             * and at the beginning of the next iteration the constraint can
             * be added globally.
             */
            SCIP_CALL( SCIPaddReoptGlbCons(scip, consdata) );

            /* create a logic-or constraint and add them to the current root */
            for(v = 0; v < consdata->nvars; v++)
            {
               if( SCIPisFeasEQ(scip, consdata->vals[v], 1) )
               {
                  assert(SCIPvarIsOriginal(consdata->vars[v]));
                  SCIP_CALL( SCIPgetNegatedVar(scip, consdata->vars[v], &consdata->vars[v]) );
                  assert(SCIPvarIsNegated(consdata->vars[v]));
               }
            }

            SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, "glbinfsub", consdata->nvars, consdata->vars,
                  TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE) );

            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );

            /* remove all dual information leading to split the node */
            SCIP_CALL( SCIPnodeReoptResetDualcons(scip, node) );

            /* we can return because nothing is to do */
            *result = SCIP_DIDNOTRUN;

            /* free buffer memory */
            SCIPfreeBufferArray(scip, &consdata->vals);
            SCIPfreeBufferArray(scip, &consdata->vars);
            SCIPfreeBuffer(scip, &consdata);

            goto TERMINATE;

            break;

         default:
            break;
      }

      /* free buffer memory */
      SCIPfreeBufferArray(scip, &consdata->vals);
      SCIPfreeBufferArray(scip, &consdata->vars);
      SCIPfreeBuffer(scip, &consdata);

      nodeID = 0;
      goto REVIVE;
   }
   else if(SCIPgetRootNode(scip) == SCIPgetCurrentNode(scip) && !SCIPnodeSplit(scip, node))
   {
      nodeID = 0;
      goto REVIVE;
   }

   assert(nodeID >= 1);

   REVIVE:

   /* get the IDs of all child nodes */
   SCIP_CALL( SCIPallocBufferArray(scip, &childids, SCIPgetNReoptNodeIDs(scip, node)) );
   SCIP_CALL( SCIPgetReoptNodeIDs(scip, node, childids, SCIPgetNReoptNodeIDs(scip, node), &nchilds) );
   assert(SCIPgetNReoptNodeIDs(scip, node) == nchilds);

   ncreatedchilds = 0;

   for(c = 0; c < nchilds; c++)
   {
      SCIP_NODE* child1;
      SCIP_NODE* child2;
      SCIP_REOPTTYPE reopttype;

      child1 = NULL;
      child2 = NULL;

      childID = childids[c];
      assert(childID >= 1);

      reopttype = SCIPreoptGetNodeType(scip, childID);

      SCIPdebugMessage("process child at ID %d\n", childID);

      /** check weather the constraints contains variable and if so, check the type of the constraint */
      if( reopttype == SCIP_REOPTTYPE_STRBRANCHED || reopttype == SCIP_REOPTTYPE_INFSUBTREE )
      {
        SCIPdebugMessage(" -> node %s\n", reopttype == SCIP_REOPTTYPE_STRBRANCHED ? "need to split" : "includes infeasible subtree");

         /* split the node into two disjoint node or cut off an infeasible subset */
         switch( reopttype ) {
            case SCIP_REOPTTYPE_STRBRANCHED:

               /** the constraint split the node into two new nodes */
               SCIP_CALL( SCIPcreateChild(scip, &child1, 0.0, SCIPgetLocalTransEstimate(scip)) );
               SCIP_CALL( SCIPcreateChild(scip, &child2, 0.0, SCIPgetLocalTransEstimate(scip)) );
               assert(child1 != NULL );
               assert(child2 != NULL );

               SCIPdebugMessage(" -> create 2 nodes: #%lld and #%lld\n", SCIPnodeGetNumber(child1), SCIPnodeGetNumber(child2));

               /* reoptimize the node:
                *  - apply the bound changes along the stored branching path in both nodes
                *  - apply the all bound changes based on primal information caught between the
                *    first and second bound change based on dual information to child1
                *  - fix all bound changes based on dual information in child1 and add the
                *    corresponding constraint to child 2
                *  - add all local constraints to both nodes */
               SCIP_CALL( SCIPapplyReopt(scip, child1, child2, childID) );

               ncreatedchilds += 2;

               break;

            case SCIP_REOPTTYPE_INFSUBTREE:

               /* create only one child node */
               SCIP_CALL( SCIPcreateChild(scip, &child1, 0.0, SCIPgetLocalTransEstimate(scip)) );
               assert(child1 != NULL);
               assert(child2 == NULL);

               SCIPdebugMessage(" -> create 1 node: #%lld\n", SCIPnodeGetNumber(child1));

               /* reoptimize the node:
                *  - apply the bound changes along the stored branching path to child1
                *  - add constraint to child 1 to cut off the infeasible subtree
                *  - add all local constraints to child1 */
               SCIP_CALL( SCIPapplyReopt(scip, NULL, child1, childID) );

               ncreatedchilds++;

               break;

            default:
               assert(reopttype == SCIP_REOPTTYPE_STRBRANCHED || reopttype == SCIP_REOPTTYPE_STRBRANCHED);
               break;
         }
      }
      else
      {
         /**
          * node at position childID includes no bound changes based on dual information
          */
         SCIP_CALL( SCIPcreateChild(scip, &child1, 0.0, SCIPgetLocalTransEstimate(scip)) );
         assert(child1 != NULL );
         assert(child2 == NULL );

         /* reoptimize the node:
          *  - apply the bound changes along the stored branching path to child1
          *  - add all local constraints to child1 */
         SCIP_CALL( SCIPapplyReopt(scip, child1, NULL, childID) );

         ncreatedchilds++;
      }

#ifdef SCIP_DISABLED_CODE
      /** set LPI state */
      if(branchruledata->savelpbasis
      && branchruledata->nodedata[childID]->lpistate != NULL
      && branchruledata->nodedata[childID]->conss == NULL
      && SCIPreoptGetSimToPrevious(scip->reopt) > 0.9)
      {
         SCIP_CALL( SCIPchildSetLpistate(child1, branchruledata->nodedata[childID]->lpistate) );
         SCIPdebugMessage("use LPI from previous round in node %lld (sim = %.4f)\n",SCIPnodeGetNumber(child1), SCIPreoptGetSimToPrevious(scip->reopt));
         printf("use LPI from previous round in node %lld (sim = %.4f)\n",SCIPnodeGetNumber(child1), SCIPreoptGetSimToPrevious(scip->reopt));
      }

      /** add local constraint from an iteration before (if some exists) to child1 */
      if(branchruledata->nodedata[childID]->conss != NULL && SCIPqueueNElems(branchruledata->nodedata[childID]->conss) - savedconsdata > 0 )
      {
         /** generate all local constraints and add them to child1 and child2 (if exists) */
         SCIP_CALL(genLocalCons(scip, branchruledata, childID, child1, child2, savedconsdata));
      }

      /** remove flag 'pseudobranched' */
      branchruledata->nodedata[childID]->dualreds = (branchruledata->nodedata[childID]->dualconscur != NULL);
#endif

      /* update the reoptid */
      SCIPnodeSetReoptID(child1, childID);

      /* set the REOPTTYPE */
      assert(reopttype >= SCIP_REOPTTYPE_TRANSIT);
      if( reopttype == SCIP_REOPTTYPE_STRBRANCHED )
         reopttype = SCIP_REOPTTYPE_TRANSIT;
      SCIPnodeSetReopttype(child1, reopttype);
      SCIPdebugMessage(" -> set reopttype: %d\n", reopttype);

      /** check if child2 includes some added constraints and save the node */
      if(child2 != NULL && SCIPnodeGetNAddedcons(child2) > 0)
      {
         SCIP_CALL( SCIPaddReoptnode(scip, child2, SCIP_REOPTTYPE_LEAF, TRUE) );
      }
   }

   if( ncreatedchilds == 0 )
      *result = SCIP_DIDNOTRUN;
   else
   {
      /* increase the counter for reoptimized nodes */
      while( ncreatedchilds > 0 )
      {
         SCIPenforceNReoptnodes(scip);
         --ncreatedchilds;
      }
      *result = SCIP_BRANCHED;
   }

   /* free the buffer memory */
   SCIPfreeBufferArray(scip, &childids);

   TERMINATE:

   SCIPdebugMessage("**** finish reoptimizing %d child nodes of node %lld ****\n", ncreatedchilds, SCIPnodeGetNumber(node));

   return SCIP_OKAY;
}

/***********************
 * non-static methods
 ***********************/

/*
 * checks if reoptimization detects infeasibility
 */
#ifdef SCIP_DISABLED_CODE
SCIP_RETCODE SCIPbranchnodereoptCheckFeasibility(
   SCIP*                 scip
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL );
   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert(branchruledata != NULL );
   assert(branchruledata->reopt);

   if( !branchruledata->initialized || !branchruledata->reopt)
      return SCIP_OKAY;

   if( (branchruledata->nodedata[0]->nodechilds == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[0]->nodechilds))
     && branchruledata->nodedata[0]->reopttype == SCIP_REOPTTYPE_TRANSIT)
   {
      branchruledata->infeasibleproblem = TRUE;
      SCIP_CALL( addInfeasibleConstraint(scip) );
      branchruledata->infeasibleconsadded = TRUE;
   }

   return SCIP_OKAY;
}
#endif

#ifdef SCIP_DISABLED_CODE
SCIP_RETCODE SCIPbranchruleNodereoptSetRootLPI(
   SCIP*                 scip
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL );
   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert(branchruledata != NULL );
   assert(branchruledata->reopt);
   assert(branchruledata->savelpbasis);
   assert(SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip));

   /** set the LPI state for the root node (if some exists) */
   if( branchruledata->nodedata[0] != NULL
    && branchruledata->nodedata[0]->lpistate != NULL
    && SCIPreoptGetSimToPrevious(scip->reopt) > 0.5)
   {
#ifdef DEBUG_MODE
      printf("** reoptimization ** use LPI from previous round in root node (sim = %.4f)\n", SCIPreoptGetSimToPrevious(scip->reopt));
#endif
      scip->tree->focuslpistate = branchruledata->nodedata[0]->lpistate;
   }

   return SCIP_OKAY;
}
#endif

/*
 * set the candidate for the WC heuristic
 */
#ifdef SCIP_DISABLED_CODE
void SCIPbranchruleNodereoptSetWCCand(
   SCIP*                 scip,
   SCIP_NODE*            node
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   int nodeID;

   assert(scip != NULL);
   assert(node != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   if( SCIPnodeGetLowerbound(node) > branchruledata->lastpruneddualbound
    && SCIPisFeasLT(scip, fabs(SCIPnodeGetLowerbound(node)), SCIPinfinity(scip)) )
   {
      branchruledata->lastpruneddualbound = SCIPnodeGetLowerbound(node);
      branchruledata->lastprunednr = SCIPnodeGetNumber(node);
      branchruledata->lastprunedID = SCIPnodeGetReoptID(node);
   }
}
#endif

/*
 * save unexplored nodes
 */
#ifdef SCIP_DISABLED_CODE
SCIP_RETCODE SCIPbranchruleNodereoptSaveOpenNodes(
   SCIP*                 scip
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_NODE** opennodes;
   int nodeID;
   int nopennodes;
   int opennode;

   assert(scip != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   SCIPdebugMessage("-> save %d open nodes.\n", SCIPgetNNodesLeft(scip));

   /* save open leaves */
   SCIP_CALL( SCIPgetLeaves(scip, &opennodes, &nopennodes) );
   for(opennode = 0; opennode < nopennodes; opennode++)
   {
      SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, opennodes[opennode], SCIP_REOPTTYPE_PRUNED, TRUE) );
   }

   /* save open children */
   SCIP_CALL( SCIPgetChildren(scip, &opennodes, &nopennodes) );
   for(opennode = 0; opennode < nopennodes; opennode++)
   {
      SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, opennodes[opennode], SCIP_REOPTTYPE_PRUNED, TRUE) );
   }

   /* save open siblings */
   SCIP_CALL( SCIPgetSiblings(scip, &opennodes, &nopennodes) );
   for(opennode = 0; opennode < nopennodes; opennode++)
   {
      SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, opennodes[opennode], SCIP_REOPTTYPE_PRUNED, TRUE) );
   }

   return SCIP_OKAY;
}
#endif


#ifdef SCIP_DISABLED_CODE
SCIP_RETCODE SCIPbranchruleNodereoptSolveLP(
   SCIP*                 scip,
   SCIP_NODE*            node,
   SCIP_Bool*            solvelp
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   int nodeID;

   assert(scip != NULL );

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL );

   /* get the ID */
   nodeID = (SCIPgetRootNode(scip) == node ? 0 : SCIPnodeGetReoptID(node));

   (*solvelp) = TRUE;

   return SCIP_OKAY;

   if( nodeID == 0 )
   {
      if( branchruledata->nodedata[0]->nodechilds != NULL && !SCIPqueueIsEmpty(branchruledata->nodedata[0]->nodechilds) )
      {
         if( SCIPreoptGetSimToPrevious(scip->reopt) >= branchruledata->objsimrootLP )
            (*solvelp) = FALSE;
      }
   }
   else
      switch (branchruledata->solvelp) {
         /* solve all LPs */
         case 0:
            if( SCIPnodeGetReopttype(node) < SCIP_REOPTTYPE_LEAF )
            {
               if( branchruledata->nodedata[nodeID]->nvars < branchruledata->solvelpdiff)
                  (*solvelp) = FALSE;
            }
            break;

         default:
            if( branchruledata->nodedata[nodeID]->nodechilds != NULL && !SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds) )
            {
               if( branchruledata->nodedata[nodeID]->nvars < branchruledata->solvelpdiff && (int) SCIPnodeGetReopttype(node) < branchruledata->solvelp )
                  (*solvelp) = FALSE;
            }
            break;
      }

   assert(*solvelp || (branchruledata->nodedata[nodeID]->nodechilds != NULL && !SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds)));

   return SCIP_OKAY;
}
#endif

#ifdef SCIP_DISABLED_CODE
SCIP_RETCODE SCIPbranchruledataNodereoptLPTimeStart(
   SCIP*                 scip
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->reopt);

   SCIP_CALL( SCIPstartClock(scip, branchruledata->lpclock) );
   SCIP_CALL( SCIPstartClock(scip, branchruledata->threadclock[branchruledata->nlp%8]) );

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPbranchruledataNodereoptLPTimeStop(
   SCIP*                 scip
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->reopt);

   SCIP_CALL( SCIPstopClock(scip, branchruledata->lpclock) );
   SCIP_CALL( SCIPstopClock(scip, branchruledata->threadclock[branchruledata->nlp%8]) );
   branchruledata->nlp++;

   return SCIP_OKAY;
}
#endif

/*
 * Callback methods of branching rule
 */

#define branchCopynodereopt NULL;
#define branchExitnodereopt NULL;

static
SCIP_DECL_BRANCHINIT(branchInitnodereopt)
{
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL );
   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL );

   /** check if reoptimization is enabled */
   if (!branchruledata->initialized)
   {
      SCIP_CALL(SCIPgetBoolParam(scip, "reoptimization/enable", &branchruledata->reopt));
   }

   /** initialize the data and change parameters */
   if (!branchruledata->initialized && branchruledata->reopt)
   {
      /** statistic */
      branchruledata->nrevivednodes = 0;
      branchruledata->nsplits = 0;

      /** mark data structure initialized */
      branchruledata->initialized = TRUE;

      /** get parameters  */
      SCIP_CALL( SCIPgetBoolParam(scip, "reoptimization/strongbranchinginit", &branchruledata->strongbranchinginit) );
   }

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreenodereopt)
{
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL );
   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL );

   /** free data structure only if reoptimization is enabled */
   if (branchruledata->initialized )
   {
      assert(branchruledata->reopt);

      branchruledata->initialized = FALSE;
      branchruledata->reopt = FALSE;
   }
   assert(!branchruledata->initialized);

   SCIPfreeMemory(scip, &branchruledata);

   return SCIP_OKAY;
}

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpnodereopt)
{
   SCIP_BRANCHRULEDATA* branchruledata;
#ifdef SCIP_DISABLED_CODE
   SCIP_VAR** branchcands;
   SCIP_Real* branchcandssol;
   SCIP_Real* branchcandsfrac;
   int nbranchcands;
#endif

   assert(branchrule != NULL );
   assert(*result != SCIP_BRANCHED);

   branchruledata = SCIPbranchruleGetData(branchrule);

   *result = SCIP_DIDNOTRUN;

   if ( branchruledata->reopt && SCIPreoptimizeNode(scip, SCIPgetCurrentNode(scip)) )
   {

#ifdef SCIP_DISABLED_CODE
      if( branchruledata->strongbranchinginit
       && SCIPreoptGetSimToPrevious(scip->reopt) < branchruledata->simsolverootlp
       && SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip))
      {
         /* get branching candidates */
         SCIP_CALL( SCIPgetLPBranchCands(scip, &branchcands, &branchcandssol, &branchcandsfrac, NULL, &nbranchcands, NULL) );

         /* run strong branching initialization */
         if( nbranchcands > 0 )
         {
            /* select only some 'good' candidates */
            SCIP_CALL( checkLPBranchCands(scip, branchruledata, 0, branchcands, branchcandssol, branchcandsfrac, &nbranchcands) );

            if( nbranchcands > 0 )
            {
               SCIP_CALL( SCIPexecRelpscostBranching(scip, TRUE, branchcands, branchcandssol, branchcandsfrac, nbranchcands, FALSE, result) );
            }
            assert(*result == SCIP_DIDNOTRUN || *result == SCIP_CUTOFF || *result == SCIP_REDUCEDDOM);
         }
      }
#endif

      if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM)
      {
         if( SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip) )
         {
            SCIP_CALL( Exec(scip, branchruledata, result) );
         }
         else
         {
            assert(1 <= SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)));

            SCIP_CALL( Exec(scip, branchruledata, result) );
         }
      }
   }

   return SCIP_OKAY;
}

/** branching execution method for external candidates */
static SCIP_DECL_BRANCHEXECEXT(branchExecextnodereopt)
{
   SCIP_BRANCHRULEDATA* branchruledata;
#ifdef SCIP_DISABALED_CODE
   SCIP_VAR** branchcands;
   SCIP_Real* branchcandssol;
   SCIP_Real* branchcandsfrac;
   int nbranchcands;
#endif

   assert(branchrule != NULL );
   assert(*result != SCIP_BRANCHED);

   branchruledata = SCIPbranchruleGetData(branchrule);

   *result = SCIP_DIDNOTRUN;

   if ( branchruledata->reopt && SCIPreoptimizeNode(scip, SCIPgetCurrentNode(scip)) )
   {
#ifdef SCIP_DISABLED_CODE
      if( branchruledata->strongbranchinginit
       && SCIPreoptGetSimToPrevious(scip->reopt) < branchruledata->simsolverootlp
       && SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip))
      {
         /* get branching candidates */
         SCIP_CALL( SCIPgetLPBranchCands(scip, &branchcands, &branchcandssol, &branchcandsfrac, NULL, &nbranchcands, NULL) );

         /* run strong branching initialization */
         if( nbranchcands > 0 )
         {
            /* select only some 'good' candidates */
            SCIP_CALL( checkLPBranchCands(scip, branchruledata, 0, branchcands, branchcandssol, branchcandsfrac, &nbranchcands) );

            if( nbranchcands > 0 )
            {
               SCIP_CALL( SCIPexecRelpscostBranching(scip, TRUE, branchcands, branchcandssol, branchcandsfrac, nbranchcands, FALSE, result) );
            }
            assert(*result == SCIP_DIDNOTRUN || *result == SCIP_CUTOFF || *result == SCIP_REDUCEDDOM);
         }
      }
#endif

      if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM)
      {
         if( SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip) )
         {
            SCIP_CALL( Exec(scip, branchruledata, result) );
         }
         else
         {
            assert(1 <= SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)));

            SCIP_CALL( Exec(scip, branchruledata, result) );
         }
      }
   }

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static SCIP_DECL_BRANCHEXECPS(branchExecpsnodereopt)
{
   SCIP_BRANCHRULEDATA* branchruledata;
#ifdef SCIP_DISABLED_CODE
   SCIP_VAR** branchcands;
   SCIP_Real* branchcandssol;
   SCIP_Real* branchcandsfrac;
   int nbranchcands;
#endif

   assert(branchrule != NULL );
   assert(*result != SCIP_BRANCHED);

   branchruledata = SCIPbranchruleGetData(branchrule);

   *result = SCIP_DIDNOTRUN;

   if ( branchruledata->reopt && SCIPreoptimizeNode(scip, SCIPgetCurrentNode(scip)) )
   {
#ifdef SCIP_DISABLED_CODE
      if( branchruledata->strongbranchinginit
       && SCIPreoptGetSimToPrevious(scip->reopt) < branchruledata->simsolverootlp
       && SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip))
      {
         /* get branching candidates */
         SCIP_CALL( SCIPgetLPBranchCands(scip, &branchcands, &branchcandssol, &branchcandsfrac, NULL, &nbranchcands, NULL) );

         /* run strong branching initialization */
         if( nbranchcands > 0 )
         {
            /* select only some 'good' candidates */
            SCIP_CALL( checkLPBranchCands(scip, branchruledata, 0, branchcands, branchcandssol, branchcandsfrac, &nbranchcands) );

            if( nbranchcands > 0 )
            {
               SCIP_CALL( SCIPexecRelpscostBranching(scip, TRUE, branchcands, branchcandssol, branchcandsfrac, nbranchcands, FALSE, result) );
            }
            assert(*result == SCIP_DIDNOTRUN || *result == SCIP_CUTOFF || *result == SCIP_REDUCEDDOM);
         }
      }
#endif

      if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM )
      {
         if( SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip) )
         {
            SCIP_CALL( Exec(scip, branchruledata, result) );
         }
         else
         {
            assert(1 <= SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)));

            SCIP_CALL( Exec(scip, branchruledata, result) );
         }
      }
   }

   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the nodereopt branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleNodereopt(
   SCIP*                 scip                     /**< SCIP data structure */
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL );

   /* create nodereopt branching rule data */
   SCIP_CALL(SCIPallocMemory(scip, &branchruledata));
   branchruledata->initialized = FALSE;

   /* include nodereopt branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC,
         BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata));

   assert(branchrule != NULL );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL(SCIPsetBranchruleFree(scip, branchrule, branchFreenodereopt));
   SCIP_CALL(SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpnodereopt));
   SCIP_CALL(SCIPsetBranchruleExecExt(scip, branchrule, branchExecextnodereopt));
   SCIP_CALL(SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsnodereopt));
   SCIP_CALL(SCIPsetBranchruleInit(scip, branchrule, branchInitnodereopt));

   return SCIP_OKAY;
}

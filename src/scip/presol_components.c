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
#define SCIP_DEBUG

/**@file   presol_components.c
 * @brief  components presolver
 * @author Dieter Weninger
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/presol_components.h"


#define PRESOL_NAME            "components"
#define PRESOL_DESC            "components presolver"
#define PRESOL_PRIORITY        -9200000      /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS              0      /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY               TRUE      /**< should presolver be delayed, if other presolvers found reductions? */

#define DEFAULT_SEARCH             TRUE      /**< should be searched for components? */
#define DEFAULT_WRITEPROBLEMS     FALSE      /**< should the single components be written as an .lp-file? */
#define DEFAULT_MAXINTVARS           20      /**< maximum number of integer (or binary) variables to solve a subproblem directly (-1: no solving) */
#define DEFAULT_NODELIMIT         10000      /**< maximum number of nodes to be solved in subproblems */
#define DEFAULT_INTFACTOR             1      /**< the weight of an integer variable compared to binary variables */

#define START_LIST_SIZE              10      /**< first size of constraint list per variable */
#define LIST_MEMORY_GAIN             10      /**< memory extension factor */

/*
 * Data structures
 */

/** control parameters */
struct SCIP_PresolData
{
   SCIP_Bool             dosearch;           /** should be searched for components? */
   SCIP_Bool             writeproblems;      /** should the single components be written as an .lp-file? */
   int                   maxintvars;         /** maximum number of integer (or binary) variables to solve a subproblem directly (-1: no solving) */
   SCIP_Longint          nodelimit;          /** maximum number of nodes to be solved in subproblems */
   SCIP_Real             intfactor;          /** the weight of an integer variable compared to binary variables */

   SCIP**                components;         /** sub-SCIPs storing the components */
   SCIP_HASHMAP**        varmaps;            /** hashmaps mapping from original variables to variables in the sub-SCIPs */
   SCIP_HASHMAP*         consmap;            /** hashmaps mapping from original constraints to constraints in the sub-SCIPs
                                              *  (needed only for performance reasons)
                                              */
   int                   componentssize;     /** size of arrays components and varmaps */
   int                   ncomponents;        /** number of components */
};

/** struct controlling memory usage of list */
struct MemoryInfo
{
   int                   navail;             /** number of available elements */
   int                   nused;              /** number of consumed elements */
};
typedef struct MemoryInfo MEMORYINFO;

/** struct saving which variable is present in which constraint */
struct List
{
   int                   nvars;              /** number of variables present in list */
   int**                 conss;              /** list of constraints for one variable */
   MEMORYINFO*           meminfo;            /** memory information of constraint list per variable */
};
typedef struct List LIST;


/*
 * Local methods
 */

/** initializes presolver data */
static
void initPresoldata(
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   assert(presoldata != NULL);

   presoldata->dosearch = 0;
   presoldata->components = NULL;
   presoldata->varmaps = NULL;
   presoldata->consmap = NULL;
   presoldata->componentssize = 0;
   presoldata->ncomponents = 0;
}

/** initialize list */
static
SCIP_RETCODE initList(
   SCIP*                 scip,
   LIST*                 list
   )
{
   int v;

   assert(scip != NULL);
   assert(list != NULL);

   list->nvars = SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &list->conss, list->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &list->meminfo, list->nvars) );
   for( v = 0; v < list->nvars; v++ )
   {
      list->conss[v] = NULL;
      list->meminfo[v].navail = 0;
      list->meminfo[v].nused = 0;
   }

   return SCIP_OKAY;
}

/** begin constraint filling for one variable */
static
SCIP_RETCODE startList(
   SCIP*                 scip,
   LIST*                 list,
   int                   varidx,
   int                   considx
   )
{
   SCIP_CALL( SCIPallocMemoryArray(scip, &list->conss[varidx], START_LIST_SIZE) );
   list->meminfo[varidx].navail = START_LIST_SIZE;
   assert(list->meminfo[varidx].navail > 1);
   list->conss[varidx][0] = considx;
   list->conss[varidx][1] = -1;
   list->meminfo[varidx].nused = 1;

   return SCIP_OKAY;
}

/** continue constraint filling for one variable */
static
SCIP_RETCODE appendList(
   SCIP*                 scip,
   LIST*                 list,
   int                   varidx,
   int                   considx
   )
{
   if( list->meminfo[varidx].nused >= (list->meminfo[varidx].navail - 1) )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &list->conss[varidx], list->meminfo[varidx].navail * LIST_MEMORY_GAIN) );
      list->meminfo[varidx].navail = list->meminfo[varidx].navail * LIST_MEMORY_GAIN;
   }
   list->conss[varidx][list->meminfo[varidx].nused] = considx;
   list->meminfo[varidx].nused++;
   list->conss[varidx][list->meminfo[varidx].nused] = -1;

   return SCIP_OKAY;
}

/** fill constraint into list for this variable */
static
SCIP_RETCODE fillConsIndex(
   SCIP*                 scip,
   LIST*                 list,
   int                   varidx,
   int                   considx
   )
{
   assert(scip != NULL);
   assert(list != NULL);
   assert(varidx < list->nvars);

   if( list->conss[varidx] == NULL )
   {
      SCIP_CALL( startList(scip, list, varidx, considx) );
   }
   else
   {
      SCIP_CALL( appendList(scip, list, varidx, considx) );
   }

   return SCIP_OKAY;
}

/** delete list */
static
void freeList(
   SCIP*                 scip,
   LIST**                list
   )
{
   int v;

   assert(scip != NULL);
   assert(list != NULL);
   assert(*list != NULL);

   if( (*list)->nvars > 0 )
   {
      assert((*list)->conss != NULL);
      assert((*list)->meminfo != NULL);

      for( v = 0; v < (*list)->nvars; v++ )
      {
         if( ((*list)->conss[v]) != NULL )
         {
            SCIPfreeMemoryArray(scip, &((*list)->conss[v]));
         }
      }

      SCIPfreeBufferArray(scip, &((*list)->conss));
      SCIPfreeBufferArray(scip, &((*list)->meminfo));

      (*list)->nvars = 0;
   }

   SCIPfreeBuffer(scip, list);
}

/** initializes the data for storing connected components */
static
SCIP_RETCODE initComponentData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   int                   ncomponents         /**< number of independent components */
   )
{
   assert(scip != NULL);
   assert(presoldata != NULL);
   assert(ncomponents > 0);
   assert(presoldata->ncomponents == 0);
   assert(presoldata->componentssize == 0);
   assert(presoldata->components == NULL);
   assert(presoldata->varmaps == NULL);
   assert(presoldata->consmap == NULL);

   /* allocate memory for sub-SCIPs and variable maps */
   SCIP_CALL( SCIPallocMemoryArray(scip, &presoldata->components, ncomponents) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &presoldata->varmaps, ncomponents) );
   SCIP_CALL( SCIPhashmapCreate(&presoldata->consmap, SCIPblkmem(scip), 10 * SCIPgetNConss(scip)) );
   presoldata->componentssize = ncomponents;

   return SCIP_OKAY;
}

/** free the data for storing connected components */
static
SCIP_RETCODE freeComponentData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   int c;

   assert(scip != NULL);
   assert(presoldata != NULL);

   /* free sub-SCIPs and variable hash maps */
   for( c = 0; c < presoldata->ncomponents; ++c )
   {
      if( presoldata->components[c] != NULL )
      {
         SCIP_CALL( SCIPfree(&presoldata->components[c]) );
      }
      if( presoldata->varmaps[c] != NULL )
      {
         SCIPhashmapFree(&presoldata->varmaps[c]);
      }
   }

   SCIPhashmapFree(&presoldata->consmap);

   SCIPfreeMemoryArray(scip, &presoldata->components);
   SCIPfreeMemoryArray(scip, &presoldata->varmaps);
   presoldata->ncomponents = 0;
   presoldata->componentssize = 0;
   presoldata->components = NULL;
   presoldata->varmaps = NULL;

   return SCIP_OKAY;
}

/** copies a connected component given by a set of constraints into a sub-SCIP */
static
SCIP_RETCODE buildComponentSubscip(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   SCIP_CONS**           conss,              /**< constraints contained in this component */
   int                   nconss,             /**< number of constraints contained in this component */
   SCIP_VAR**            vars,               /**< variables contained in this component */
   int                   nvars,              /**< number of variables contained in this component */
   int*                  nsolvedprobs,       /**< pointer to increase, if the subproblem was solved */
   int*                  ndeletedvars,       /**< pointer to store the number of deleted variables */
   int*                  ndeletedcons,       /**< pointer to store the number of deleted constraints */
   SCIP_Real*            subsolvetime        /**< pointer to store time needed to solve the subproblem */
   )
{
   char probname[SCIP_MAXSTRLEN];
   char consname[SCIP_MAXSTRLEN];
   SCIP* subscip;
   SCIP_CONS* newcons;
   SCIP_Real timelimit;
   SCIP_Real memorylimit;
   SCIP_Bool success;
   int c;
   int i;

   assert(scip != NULL);
   assert(presoldata != NULL);
   assert(conss != NULL);
   assert(nconss > 0);
   assert(nvars > 0);

   c = presoldata->ncomponents;

   //printf("build sub-SCIP for component %d (%d vars, %d conss)\n", c, nvars, nconss);

   /* create sub-SCIP */
   SCIP_CALL( SCIPcreate(&(presoldata->components[c])) );
   subscip = presoldata->components[c];

   /* create variable hashmap */
   SCIP_CALL( SCIPhashmapCreate(&presoldata->varmaps[c], SCIPblkmem(scip), 10 * nvars) );

   /* copy plugins */
   success = TRUE;
   SCIP_CALL( SCIPcopyPlugins(scip, subscip,
         TRUE, /* readers */
         TRUE, /* pricers */
         TRUE, /* conshdlrs */
         TRUE, /* conflicthdlrs */
         TRUE, /* presolvers */
         TRUE, /* relaxators */
         TRUE, /* separators */
         TRUE, /* propagators */
         TRUE, /* heuristics */
         TRUE, /* eventhandler */
         TRUE, /* nodeselectors (SCIP gives an error if there is none) */
         TRUE, /* branchrules */
         TRUE, /* displays */
         FALSE, /* dialogs */
         TRUE, /* nlpis */
         &success) );

   if( success )
   {
      /* copy parameter settings */
      SCIP_CALL( SCIPcopyParamSettings(scip, subscip) );

#if 1
      /* reduce the effort spent for hash tables */
      SCIP_CALL( SCIPsetBoolParam(subscip, "misc/usevartable", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(subscip, "misc/useconstable", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(subscip, "misc/usesmalltables", TRUE) );
#endif

      /* set gap limit to 0 */
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/gap", 0.0) );

      /* do not catch control-C */
      SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );


      /* check whether there is enough time and memory left */
      timelimit = 0.0;
      memorylimit = 0.0;
      SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
      if( !SCIPisInfinity(scip, timelimit) )
         timelimit -= SCIPgetSolvingTime(scip);
      SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
      if( !SCIPisInfinity(scip, memorylimit) )
         memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
      if( timelimit <= 0.0 || memorylimit <= 0.0 )
         goto TERMINATE;

      /* set limits for the subproblem */
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );

      /* set node limit */
      if( presoldata->nodelimit != -1 )
      {
         SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", presoldata->nodelimit) );
      }
      //SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", 20000) );

      /* disable output */
      SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", SCIP_VERBLEVEL_NONE) );

      /* create problem in sub-SCIP */
      /* get name of the original problem and add "comp_nr" */
      (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_comp_%d", SCIPgetProbName(scip), c);
      SCIP_CALL( SCIPcreateProb(subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

      for( i = 0; i < nconss; ++i )
      {
         /* copy the constraint */
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s", SCIPconsGetName(conss[i]));
         SCIP_CALL( SCIPgetConsCopy(scip, subscip, conss[i], &newcons, SCIPconsGetHdlr(conss[i]),
               presoldata->varmaps[c], presoldata->consmap, consname,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, &success) );

         /* break if constraint was not successfully copied */
         if( !success )
            break;

         SCIP_CALL( SCIPaddCons(subscip, newcons) );
         SCIP_CALL( SCIPreleaseCons(subscip, &newcons) );
      }
   }

   /* ignore this component, if a problem relevant plugin or a constraint could not be copied */
   if( success )
   {
      presoldata->ncomponents++;

      //printf("++++++++++++++ sub-SCIP for component %d: %d vars (%d bin, %d int, %d impl, %d cont), %d conss\n",
      //   c, nvars, SCIPgetNBinVars(subscip), SCIPgetNIntVars(subscip), SCIPgetNImplVars(subscip), SCIPgetNContVars(subscip), nconss);

      if( presoldata->writeproblems )
      {
         (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_comp_%d.lp", SCIPgetProbName(scip), c);
         printf("write problem to file %s\n", probname);
         SCIP_CALL( SCIPwriteOrigProblem(subscip, probname, NULL, FALSE) );
      }

      if( SCIPgetNBinVars(subscip) + SCIPgetNIntVars(subscip) <= presoldata->maxintvars )
      {
         //SCIP_CALL( SCIPpresolve(subscip) ); // TODO: simulation of presolving without solve

         SCIP_CALL( SCIPsolve(subscip) );

         //printf("solved subproblem %d: status = %d, time = %.2f\n", c, SCIPgetStatus(subscip), SCIPgetSolvingTime(subscip));
         *subsolvetime += SCIPgetSolvingTime(subscip);

         if( SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL )
         {
            SCIP_SOL* sol;
            SCIP_Bool infeasible;
            SCIP_Bool fixed;

            ++(*nsolvedprobs);

            sol = SCIPgetBestSol(subscip);

            /* fix variables contained in the sub-scip */
            for( i = 0; i < nvars; ++i )
            {
               assert( SCIPhashmapExists(presoldata->varmaps[c], vars[i]) );
               SCIP_CALL( SCIPfixVar(scip, vars[i], SCIPgetSolVal(subscip, sol, SCIPhashmapGetImage(presoldata->varmaps[c], vars[i])), &infeasible, &fixed) );
               assert(!infeasible);
               assert(fixed);
            }
            (*ndeletedvars) += nvars;

            /* delete constraints contained in the sub-scip */
            for( i = 0; i < nconss; ++i )
            {
               SCIP_CALL( SCIPdelCons(scip, conss[i]) );
            }
            (*ndeletedcons) += nconss;
         }
         else
         {
            printf("++++++++++++++ sub-SCIP for component %d not solved (status=%d, time=%.2f): %d vars (%d bin, %d int, %d impl, %d cont), %d conss\n",
               c, SCIPgetStatus(subscip), SCIPgetSolvingTime(subscip), nvars, SCIPgetNBinVars(subscip), SCIPgetNIntVars(subscip), SCIPgetNImplVars(subscip), SCIPgetNContVars(subscip), nconss);
         }
      }
      else
      {
         printf("++++++++++++++ sub-SCIP for component %d not solved: %d vars (%d bin, %d int, %d impl, %d cont), %d conss\n",
            c, nvars, SCIPgetNBinVars(subscip), SCIPgetNIntVars(subscip), SCIPgetNImplVars(subscip), SCIPgetNContVars(subscip), nconss);
      }
   }

 TERMINATE:
   SCIP_CALL( SCIPfree(&presoldata->components[c]) );
   SCIPhashmapFree(&presoldata->varmaps[c]);
   presoldata->components[c] = NULL;
   presoldata->varmaps[c] = NULL;

   return SCIP_OKAY;
}


/** loop over constraints, get active variables and fill adjacency list */
static
SCIP_RETCODE fillAdjList(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ADJLIST*         adjlist,            /**< adjacency list */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   LIST*                 conslist,           /**< list for filling which variable is in which constraints */
   SCIP_Bool*            success             /**< flag indicating successful adjacency list filling */
   )
{
   int c;
   int v;
   SCIP_VAR** consvars;
   int nconsvars;
   int varssize;
   int requiredsize;
   int probindex1;
   int probindex2;

   assert(scip != NULL);
   assert(adjlist != NULL);
   assert(conss != NULL);
   assert(conslist != NULL);
   assert(success != NULL);

   *success = TRUE;

   nconsvars = 0;
   requiredsize = 0;
   varssize = SCIPgetNVars(scip);

   /* use big buffer for storing variables
      TODO: maybe we can do better with less memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, varssize) );

   for( c = 0; c < nconss; c++ )
   {
      SCIP_CALL( SCIPgetConsNVars(scip, conss[c], &nconsvars, success) );

      if( !(*success) )
      {
         break;
      }

      SCIP_CALL( SCIPgetConsVars(scip, conss[c], consvars, varssize, success) );

      if( !(*success) )
      {
         break;
      }

      SCIP_CALL( SCIPgetActiveVars(scip, consvars, &nconsvars, varssize, &requiredsize) );

      if( requiredsize > varssize )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &consvars, requiredsize) );
         varssize = requiredsize;

         /* call get active variables a second time */
         SCIP_CALL( SCIPgetActiveVars(scip, consvars, &nconsvars, varssize, &requiredsize) );
      }

      probindex1 = SCIPvarGetProbindex(consvars[0]);
      assert(probindex1 >= 0);
      SCIP_CALL( fillConsIndex(scip, conslist, probindex1, c) );

      if( nconsvars > 1 )
      {
         /* create sparse adjacency list
            sparse means, to add only those edges necessary for component calculation */
         for( v = 1; v < nconsvars; v++ )
         {
            probindex2 = SCIPvarGetProbindex(consvars[v]);
            assert(probindex2 >= 0);
            SCIP_CALL( fillConsIndex(scip, conslist, probindex2, c) );

            /* we add a directed edge in both directions */
            SCIP_CALL( SCIPadjlistAddEdgeSave(adjlist, probindex1, probindex2) );
            SCIP_CALL( SCIPadjlistAddEdgeSave(adjlist, probindex2, probindex1) );
         }
      }
   }

   SCIPfreeBufferArray(scip, &consvars);

   return SCIP_OKAY;
}

/** use components to detect which vars and cons belong to one subscip
 *  and try to solve all subscips having not too much integer vars
 */
static
SCIP_RETCODE createSubScipsAndSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int*                  components,         /**< array with component number for every variable */
   int*                  ncomponents,        /**< size of array with component numbers */
   LIST*                 conslist,           /**< list containing information which variable occur in which constraint */
   int*                  ndeletedcons,       /**< number of deleted constraints by component presolver */
   int*                  ndeletedvars,       /**< number of fixed variables by component presolver */
   int*                  nsolvedprobs        /**< number of solved subproblems */
   )
{
   int comp;
   int v;
   int c;
   SCIP_VAR** vars;
   int nvars;
   SCIP_CONS** tmpconss;
   int ntmpconss;
   SCIP_VAR** tmpvars;
   int ntmpvars;
   int nbinvars;
   int nintvars;
   SCIP_Real subsolvetime;
   SCIP_Bool* consincomponent;

   assert(scip != NULL);
   assert(presoldata != NULL);
   assert(conss != NULL);
   assert(components != NULL);
   assert(ncomponents != NULL);
   assert(conslist != NULL);
   assert(ndeletedcons != NULL);
   assert(ndeletedvars != NULL);
   assert(nsolvedprobs != NULL);

   *ndeletedcons = 0;
   *ndeletedvars = 0;
   *nsolvedprobs = 0;
   subsolvetime = 0.0;

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   initComponentData(scip, presoldata, *ncomponents);

   SCIP_CALL( SCIPallocBufferArray(scip, &tmpconss, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consincomponent, nconss) );

   /* loop over all components
      start loop from 1 because components are numbered form 1..n */
   for( comp = 1; comp <= *ncomponents; comp++ )
   {
      /* initially, no constraint is in this component */
      BMSclearMemoryArray(consincomponent, nconss);

      ntmpconss = 0;
      ntmpvars = 0;
      nbinvars = 0;
      nintvars = 0;

      for( v = 0; v < nvars; v++ )
      {
         if( components[v] == comp )
         {
            /* variable is present in this component */
            tmpvars[ntmpvars] = vars[v];

            /* check whether variable is of binary or integer type */
            if( SCIPvarGetType(tmpvars[ntmpvars]) == SCIP_VARTYPE_BINARY )
               nbinvars++;
            else if( SCIPvarGetType(tmpvars[ntmpvars]) == SCIP_VARTYPE_INTEGER )
               nintvars++;

            ++ntmpvars;

            if( conslist->conss[v] != NULL )
            {
               c = 0;
               while( conslist->conss[v][c] != -1 )
               {
                  consincomponent[conslist->conss[v][c]] = TRUE;
                  c++;
               }
            }
         }
      }

      /* get constraints for this component */
      ntmpconss = 0;
      for( c = 0; c < nconss; c++ )
      {
         if( consincomponent[c] == TRUE )
         {
            tmpconss[ntmpconss] = conss[c];
            ntmpconss++;
         }
      }

      if( (nbinvars + presoldata->intfactor * nintvars <= presoldata->maxintvars) || presoldata->writeproblems )
      {
         /* build subscip for one component */
         SCIP_CALL( buildComponentSubscip(scip, presoldata, tmpconss, ntmpconss,
               tmpvars, ntmpvars, nsolvedprobs, ndeletedvars, ndeletedcons, &subsolvetime) );
      }
      else
      {
         printf("++++++++++++++ sub-SCIP for component %d not created: %d vars (%d bin, %d int, %d cont), %d conss\n",
            presoldata->ncomponents, ntmpvars, nbinvars, nintvars, ntmpvars - nintvars - nbinvars, ntmpconss);
      }
   }

   SCIPfreeBufferArray(scip, &consincomponent);
   SCIPfreeBufferArray(scip, &tmpvars);
   SCIPfreeBufferArray(scip, &tmpconss);

   freeComponentData(scip, presoldata);

   return SCIP_OKAY;
}

/** performs presolving by searching for components */
static
SCIP_RETCODE presolComponents(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_PRESOL*          presol,             /**< the presolver itself */
   int*                  nfixedvars,         /**< pointer to total number of variables fixed of all presolvers */
   int*                  ndelconss,          /**< pointer to total number of deleted constraints of all presolvers */
   SCIP_RESULT*          result              /**< pointer to store the result of the presolving call */
   )
{
   SCIP_CONSHDLR** conshdlrs;
   int nconshdlrs;
   SCIP_CONS** conss;
   int nconss;
   int nvars;
   SCIP_PRESOLDATA* presoldata;
   SCIP_ADJLIST* adjlist;
   int* components;
   int ncomponents;
   int ndeletedcons;
   int ndeletedvars;
   int nsolvedprobs;
   int i;
   int c;
   LIST* conslist;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(presol != NULL);
   //assert(nfixedvars != NULL);
   //assert(ndelconss != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING || SCIPinProbing(scip) )
      return SCIP_OKAY;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   if( !presoldata->dosearch )
   {
      /* do not search for components */
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   ncomponents = 0;
   ndeletedvars = 0;
   ndeletedcons = 0;
   nsolvedprobs = 0;

   /* collect number of enforced constraints */
   nconshdlrs = SCIPgetNConshdlrs(scip);
   conshdlrs = SCIPgetConshdlrs(scip);
   nconss = 0;
   for( i = 0; i < nconshdlrs; ++i )
   {
      nconss += SCIPconshdlrGetNEnfoConss(conshdlrs[i]);
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &conss, nconss) );

   /* copy the contraints */
   nconss = 0;
   for( i = 0; i < nconshdlrs; ++i )
   {
      SCIP_CONS** conshdlrconss;
      int nconshdlrconss;

      conshdlrconss = SCIPconshdlrGetEnfoConss(conshdlrs[i]);
      nconshdlrconss = SCIPconshdlrGetNEnfoConss(conshdlrs[i]);

      for( c = 0; c < nconshdlrconss; ++c )
      {
         conss[nconss] = conshdlrconss[c];
         nconss++;
      }
   }

   nvars = SCIPgetNVars(scip);

   if( nvars > 1 && nconss > 1 )
   {
      SCIP_CALL( SCIPallocBuffer(scip, &conslist) );
      initList(scip, conslist);

      SCIP_CALL( SCIPadjlistCreate(&adjlist, nvars) );

      SCIP_CALL( fillAdjList(scip, adjlist, conss, nconss, conslist, &success) );

      if( success )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &components, nvars) );

         SCIP_CALL( SCIPadjlistComputeComponents(adjlist, components, &ncomponents) );

         SCIP_CALL( createSubScipsAndSolve(scip, presoldata, conss, nconss,
               components, &ncomponents, conslist, &ndeletedcons, &ndeletedvars, &nsolvedprobs ) );

         SCIPfreeBufferArray(scip, &components);
      }

      SCIPadjlistFree(&adjlist);
      freeList(scip, &conslist);
   }

   SCIPfreeBufferArray(scip, &conss);

   //(*nfixedvars) = (*nfixedvars) + ndeletedvars;
   //(*ndelconss) = (*ndelconss) + ndeletedcons;

   SCIPdebugMessage("### %d components, %d solved, %d delcons, %d delvars\n",
      ncomponents, nsolvedprobs, ndeletedcons, ndeletedvars);

   return SCIP_OKAY;
}


/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_PRESOLCOPY(presolCopyComponents)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of components presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolCopyComponents NULL
#endif


/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeComponents)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   SCIPfreeMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}


/** initialization method of presolver (called after problem was transformed) */
#if 0
static
SCIP_DECL_PRESOLINIT(presolInitComponents)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of components presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolInitComponents NULL
#endif


/** deinitialization method of presolver (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_PRESOLEXIT(presolExitComponents)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of components presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolExitComponents NULL
#endif


/** presolving initialization method of presolver (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_PRESOLINITPRE(presolInitpreComponents)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of components presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolInitpreComponents NULL
#endif


/** presolving deinitialization method of presolver (called after presolving has been finished) */
static
SCIP_DECL_PRESOLEXITPRE(presolExitpreComponents)
{  /*lint --e{715}*/

   SCIP_CALL( presolComponents(scip, presol, NULL, NULL, result) );

   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecComponents)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the components presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolComponents(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;

   /* create components presolver data */
   SCIP_CALL( SCIPallocMemory(scip, &presoldata) );
   initPresoldata(presoldata);

   /* include presolver */
   SCIP_CALL( SCIPincludePresol(scip,
         PRESOL_NAME,
         PRESOL_DESC,
         PRESOL_PRIORITY,
         PRESOL_MAXROUNDS,
         PRESOL_DELAY,
         presolCopyComponents,
         presolFreeComponents,
         presolInitComponents,
         presolExitComponents,
         presolInitpreComponents,
         presolExitpreComponents,
         presolExecComponents,
         presoldata) );

   /* add presolver parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/components/dosearch",
         "search for components (0: no search, 1: do search)",
         &presoldata->dosearch, FALSE, DEFAULT_SEARCH, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/components/writeproblems",
         "should the single components be written as an .lp-file?",
         &presoldata->writeproblems, FALSE, DEFAULT_WRITEPROBLEMS, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/components/maxintvars",
         "maximum number of integer (or binary) variables to solve a subproblem directly (-1: no solving)",
         &presoldata->maxintvars, FALSE, DEFAULT_MAXINTVARS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip,
         "presolving/components/nodelimit",
         "maximum number of nodes to be solved in subproblems",
         &presoldata->nodelimit, FALSE, DEFAULT_NODELIMIT, -1, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/components/intfactor",
         "the weight of an integer variable compared to binary variables",
         &presoldata->intfactor, FALSE, DEFAULT_INTFACTOR, 0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

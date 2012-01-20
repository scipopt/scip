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
/*#define SCIP_DEBUG*/

/**@file    presol_concomp.c
 * @ingroup PRESOLVERS
 * @brief   small connected components presolver
 * @author  Dieter Weninger
 * @author  Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "matrix/matrix.h"
#include "scip/presol_concomp.h"

#define PRESOL_NAME             "concomp"
#define PRESOL_DESC             "connected components presolver"
#define PRESOL_PRIORITY         -9200000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS              -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY                TRUE /**< should presolver be delayed, if other presolvers found reductions? */

#define DEFAULT_SEARCH              TRUE /**< should be searched for connected components? */
#define DEFAULT_SHOW_SIZES         FALSE /**< should the sizes of the single components be shown? */
#define DEFAULT_WRITEPROBLEMS      FALSE /**< should the single components be written as an .lp-file? */
#define DEFAULT_MAXINTVARS            20 /**< maximum number of integer (or binary) variables to solve a subproblem directly (-1: no solving) */
#define DEFAULT_NODELIMIT          10000 /**< maximum number of nodes to be solved in subproblems */
#define DEFAULT_INTFACTOR              1 /**< the weight of an integer variable compared to binary variables */

#define START_ADJSIZE                 10 /**< first size of adjacency list */
#define ADJLIST_MEMORY_GAIN           10 /**< memory extension factor */
#define USE_NO_REDUNDANT_EDGES      TRUE /**< TRUE := do not add edges redundant, FALSE := redundant edge adding is allowed */

/*
 * Data structures
 */

/** control parameter for concomp */
struct SCIP_PresolData
{
   SCIP_Bool             dosearch;           /** should be searched for connected components? */
   SCIP_Bool             showsizes;          /** should the sizes of the single components be shown? */
   SCIP_Bool             writeproblems;      /** should the single components be written as an .lp-file? */
   int                   maxintvars;         /** maximum number of integer (or binary) variables to solve a subproblem directly (-1: no solving) */
   SCIP_Longint          nodelimit;          /** maximum number of nodes to be solved in subproblems */
   SCIP_Real             intfactor;          /** the weight of an integer variable compared to binary variables */

   SCIP**                components;         /** sub-SCIPs storing the connected components */
   SCIP_HASHMAP**        varmaps;            /** hashmaps mapping from original variables to variables in the sub-SCIPs */
   SCIP_HASHMAP*         consmap;            /** hashmaps mapping from original constraints to constraints in the sub-SCIPs
                                              *  (needed only for performance reasons)
                                              */
   int                   componentssize;     /** size of arrays components and varmaps */
   int                   ncomponents;        /** number of connected components */
};


/** struct controlling memory usage of adjacency list */
struct ListInfo
{
   int                   nmemory;            /** amount of memory allocated */
   int                   nmemused;           /** amount of memory currently in use */
};
typedef struct ListInfo LISTINFO;

/** graph data structure */
struct Graph
{
   int                   nvertices;          /** number of vertices */
   int                   nedges;             /** number of edges doubled */
   int                   ncomponents;        /** number of connected components */

   int*                  visited;            /** status of vertices: -1 := not initialized, 0 := initialized, 1 := visited with DFS */
   int**                 adjlist;            /** adjacency list to every vertice */
   LISTINFO*             adjlistinfo;        /** additional info to every adjlist */
   int*                  marks;              /** counter for connected components */
};
typedef struct Graph GRAPH;


/*
 * Local methods
 */

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
   SCIP*                 scip,
   SCIP_PRESOLDATA*      presoldata
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
   int*                  ndeletedconss,      /**< pointer to store the number of deleted constraints */
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
            (*ndeletedconss) += nconss;
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

/** initializes presolver data */
static
void initPresoldata(
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   assert(presoldata != NULL);

   presoldata->dosearch = 0;
   presoldata->showsizes = 0;

   presoldata->components = NULL;
   presoldata->varmaps = NULL;
   presoldata->consmap = NULL;
   presoldata->componentssize = 0;
   presoldata->ncomponents = 0;
}

/** use graph data structure to create sub-SCIPs */
static
SCIP_RETCODE createSubscips(
   SCIP*                 scip,               /**< SCIP data structure */
   CONSTRAINTMATRIX*     matrix,
   GRAPH*                graph,
   SCIP_PRESOLDATA*      presoldata
   )
{
   SCIP_CONS** conss;
   int nconss;
   int nvars;
   SCIP_CONS** tmpconss;
   SCIP_VAR** tmpvars;
   int ntmpconss;
   int ntmpvars;
   SCIP_Bool* consincomp;
   int comp;
   int v;
   int c;
   int* colpnt;
   int* colend;
   int nbinvars;
   int nintvars;
   int nsolvedprobs;
   int ndeletedvars;
   int ndeletedconss;
   SCIP_Real subsolvetime;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(graph != NULL);
   assert(graph->ncomponents >= 1);

   conss = matrix->conss;
   nconss = matrix->nrows;
   nvars = matrix->ncols;
   nsolvedprobs = 0;
   ndeletedvars = 0;
   ndeletedconss = 0;
   subsolvetime = 0.0;

   assert(matrix->nrows == nconss);

   initComponentData(scip, presoldata, graph->ncomponents);

   SCIP_CALL( SCIPallocBufferArray(scip, &tmpconss, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consincomp, nconss) );

   /* loop over all connected components */
   for( comp = 0; comp < graph->ncomponents; comp++ )
   {
      /* initially, no constraint is in this component */
      BMSclearMemoryArray(consincomp, nconss);

      ntmpvars = 0;
      nbinvars = 0;
      nintvars = 0;

      /* loop over all vertices (variables) */
      for( v = 0; v < graph->nvertices; v++ )
      {
         if( graph->marks[v] == comp )
         {
            /* variable is present in this component */
            tmpvars[ntmpvars] = matrix->vars[v];

            /* check whether variable is of binary or integer type */
            if( SCIPvarGetType(tmpvars[ntmpvars]) == SCIP_VARTYPE_BINARY )
               nbinvars++;
            else if( SCIPvarGetType(tmpvars[ntmpvars]) == SCIP_VARTYPE_INTEGER )
               nintvars++;

            ++ntmpvars;

            colpnt = matrix->colmatind + matrix->colmatbeg[v];
            colend = colpnt + matrix->colmatcnt[v];
            for( ; colpnt < colend; colpnt++ )
            {
               assert(colpnt != NULL);
               assert(*colpnt < matrix->nrows);
               consincomp[*colpnt] = TRUE;
            }
         }
      }

      /* get constraints for this component */
      ntmpconss = 0;
      for( c = 0; c < nconss; c++ )
      {
         if( consincomp[c] == TRUE )
         {
            tmpconss[ntmpconss] = conss[c];
            ntmpconss++;
         }
      }

      //printf("++++++++++++++ sub-SCIP for component %d: %d vars (%d bin, %d int, %d cont), %d conss\n",
      //   presoldata->ncomponents, ntmpvars, nbinvars, nintvars, ntmpvars - nintvars - nbinvars, ntmpconss);

      if( (nbinvars + presoldata->intfactor * nintvars <= presoldata->maxintvars) || presoldata->writeproblems )
      {
         /* build subscip for one component */
         SCIP_CALL( buildComponentSubscip(scip, presoldata, tmpconss, ntmpconss, tmpvars, ntmpvars, &nsolvedprobs, &ndeletedvars, &ndeletedconss, &subsolvetime) );
      }
      else
      {
         printf("++++++++++++++ sub-SCIP for component %d not created: %d vars (%d bin, %d int, %d cont), %d conss\n",
            presoldata->ncomponents, ntmpvars, nbinvars, nintvars, ntmpvars - nintvars - nbinvars, ntmpconss);
      }
   }

   printf("solved %d/%d subproblems: deleted %d vars, %d conss, subproblem solving time %.2f\n",
      nsolvedprobs, graph->ncomponents, ndeletedvars, ndeletedconss, subsolvetime);

   SCIPfreeBufferArray(scip, &consincomp);
   SCIPfreeBufferArray(scip, &tmpvars);
   SCIPfreeBufferArray(scip, &tmpconss);

   freeComponentData(scip, presoldata);

   return SCIP_OKAY;
}


/*
 * functions for graph datastruct
 */
/** start adjacency list for head vertex and fill one neighboring vertex */
static
SCIP_RETCODE startLine(
   SCIP*                 scip,
   GRAPH*                graph,
   int                   head,
   int                   neighbor
   )
{
   assert(scip != NULL);
   assert(graph != NULL);

   SCIP_CALL( SCIPallocMemoryArray(scip, &graph->adjlist[head], START_ADJSIZE) );
   graph->adjlistinfo[head].nmemory = START_ADJSIZE;
   assert(graph->adjlistinfo[head].nmemory > 1);
   graph->adjlist[head][0] = neighbor;
   graph->adjlist[head][1] = -1;
   graph->adjlistinfo[head].nmemused = 1;
   graph->nedges++;

   return SCIP_OKAY;
}

/** continue adjacency list for head vertex and fill another neighboring vertex */
static
SCIP_RETCODE appendLine(
   SCIP*                 scip,
   GRAPH*                graph,
   int                   head,
   int                   neighbor
   )
{
   int* node;
   int present;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(0 <= head);
   assert(head < graph->nvertices);
   assert(0 <= neighbor);
   assert(neighbor < graph->nvertices);

   node = graph->adjlist[head];
   present = 0;

   while( USE_NO_REDUNDANT_EDGES && *node != -1 )
   {
      if( *node == neighbor )
      {
         present = 1;
         break;
      }
      node++;
   }

   if( !present && (*node != neighbor) )
   {
      if( graph->adjlistinfo[head].nmemused >= (graph->adjlistinfo[head].nmemory - 1) )
      {
         SCIP_CALL( SCIPreallocMemoryArray(scip, &graph->adjlist[head], graph->adjlistinfo[head].nmemory * ADJLIST_MEMORY_GAIN) );
         graph->adjlistinfo[head].nmemory = graph->adjlistinfo[head].nmemory * ADJLIST_MEMORY_GAIN;
      }
      graph->adjlist[head][graph->adjlistinfo[head].nmemused] = neighbor;
      graph->adjlistinfo[head].nmemused++;
      graph->adjlist[head][graph->adjlistinfo[head].nmemused] = -1;
      graph->nedges++;
   }

   return SCIP_OKAY;
}

/** add two edges to the graph: vertex1 -> vertex2 and vertex2 -> vertex1 */
static
SCIP_RETCODE addEdgeGraph(
   SCIP*                 scip,
   GRAPH*                graph,
   int                   vertex1,
   int                   vertex2
   )
{
   assert(scip != NULL);
   assert(graph != NULL);
   assert(vertex1 < graph->nvertices);
   assert(vertex2 < graph->nvertices);
   assert(vertex1 != vertex2);

   /* edge v1 - v2 */
   if( graph->adjlist[vertex1] == NULL )
   {
      startLine(scip, graph, vertex1, vertex2);
   }
   else
   {
      appendLine(scip, graph, vertex1, vertex2);
   }

   /* edge v2 - v1 */
   if( graph->adjlist[vertex2] == NULL )
   {
      startLine(scip, graph, vertex2, vertex1);
   }
   else
   {
      appendLine(scip, graph, vertex2, vertex1);
   }

   return SCIP_OKAY;
}

/** initialize graph out of the matrix */
static
SCIP_RETCODE initGraph(
   SCIP*                 scip,
   CONSTRAINTMATRIX*     matrix,
   GRAPH*                graph,
   SCIP_Bool*            initialized
   )
{
   int v;
   int* rowpnt1;
   int* rowpnt2;
   int* rowend;
   assert(graph != NULL);

   SCIPdebugMessage("Entering init graph: rows=%d, cols=%d, nonzs=%d.\n",
      matrix->nrows, matrix->ncols, matrix->nnonzs);

   *initialized = FALSE;
   rowpnt1 = NULL;
   rowpnt2 = NULL;
   rowend = NULL;

   graph->nvertices = matrix->ncols;
   graph->nedges = 0;
   graph->ncomponents = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &graph->visited, graph->nvertices) );
   SCIP_CALL( SCIPallocBufferArray(scip, &graph->adjlist, graph->nvertices) );
   SCIP_CALL( SCIPallocBufferArray(scip, &graph->adjlistinfo, graph->nvertices) );
   SCIP_CALL( SCIPallocBufferArray(scip, &graph->marks, graph->nvertices) );
   for( v = 0; v < graph->nvertices; v++ )
   {
      graph->visited[v] = -1;
      graph->adjlist[v] = NULL;
      graph->adjlistinfo[v].nmemory = 0;
      graph->adjlistinfo[v].nmemused = 0;
      graph->marks[v] = -1;
   }

   /* add edges  */
   for( v = 0; v < matrix->nrows; v++ )
   {
      if( matrix->rowmatcnt[v] > 1 )
      {
         /* components with minimum two variables */
         rowpnt1 = matrix->rowmatind + matrix->rowmatbeg[v];
         rowend = rowpnt1 + matrix->rowmatcnt[v];
         for(; rowpnt1 < rowend - 1; rowpnt1++ )
         {
            /* create sparse-edge graph: i.e. we only add edges
               in a manner for finding the connected components */
            rowpnt2=rowpnt1+1;
            graph->visited[*rowpnt1] = 0;
            graph->visited[*rowpnt2] = 0;
            addEdgeGraph(scip, graph, *rowpnt1, *rowpnt2);
         }
      }
      else if( matrix->rowmatcnt[v] == 1 )
      {
         /* components containing only one variable */
         rowpnt1 = matrix->rowmatind + matrix->rowmatbeg[v];
         graph->visited[*rowpnt1] = 0;
      }
   }
   *initialized = TRUE;

   SCIPdebugMessage("Leaving init graph: vertices=%d, edges=%d.\n",
      graph->nvertices, graph->nedges);

   return SCIP_OKAY;
}

/** recursive depth-first-search algorithm */
static
void dfsGraph(
   GRAPH*                graph,
   int                   v
   )
{
   int* node;

   assert(graph != NULL);

   graph->marks[v] = graph->ncomponents;
   graph->visited[v] = 1;

   if( graph->adjlist[v] != NULL)
   {
      node = graph->adjlist[v];

      while( *node != -1 )
      {
         if( graph->visited[*node] == 0 )
         {
            dfsGraph(graph,*node);
         }
         node++;
      }
   }
}

/** calculate the connected components with DFS */
static
void calcGraphConComp(
   GRAPH*                graph
   )
{
   int v;

   assert(graph != NULL);
   SCIPdebugMessage("Entering graph calculated connected components.\n");

   graph->ncomponents = 0;
   for( v = 0; v < graph->nvertices; v++ )
   {
      if( graph->visited[v] == 0 )
      {
         dfsGraph(graph,v);
         graph->ncomponents++;
      }
   }

   SCIPdebugMessage("Leaving graph calculated connected components.\n");
}

#if 0
/* debug function */
static
void writeGraph(
   SCIP*                 scip,
   CONSTRAINTMATRIX*     matrix,
   GRAPH*                graph
   )
{
   int v;
   int* node;
   FILE* fid;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(graph != NULL);
   fid = NULL;
   fid = fopen("graph.txt","wt");

   for( v = 0; v < graph->nvertices; v++ )
   {
      fprintf(fid,"[v%d/%s/%d] => ",v,SCIPvarGetName(matrix->vars[v]),
         /* representation within objective function */
         (SCIPisEQ(scip,SCIPvarGetObj(matrix->vars[v]),0.0))>0?0:1);
      if( graph->adjlist[v] != NULL)
      {
         node = graph->adjlist[v];
         while( *node != -1 )
         {
            fprintf(fid, "[%d]/%s ",*node,SCIPvarGetName(matrix->vars[*node]));
            node++;
         }
      }
      fprintf(fid,"\n");
   }
   fclose(fid);
}
#endif

#if 0
/* debug function */
static
void writeConComp(
   SCIP*                 scip,
   CONSTRAINTMATRIX*     matrix,
   GRAPH*                graph
   )
{
   int c;
   int v;
   FILE* fid;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(graph != NULL);
   fid = NULL;
   fid = fopen("coco.txt","wt");

   for( c = 0; c < graph->ncomponents; c++ )
   {
      fprintf(fid,"[c%d] => ",c);
      for( v = 0; v < graph->nvertices; v++ )
      {
         if( graph->marks[v] == c )
         {
            fprintf(fid,"[%d]/%s ",v,SCIPvarGetName(matrix->vars[v]));
         }
      }
      fprintf(fid,"\n");
   }
   fclose(fid);
}
#endif

/** calculate the number of vertices for component idx */
static
int sizeGraphComponent(
   GRAPH*                graph,
   int                   idx
   )
{
   int v;
   int size;

   assert(graph != NULL);
   assert(idx < graph->ncomponents);

   size = 0;
   for( v = 0; v < graph->nvertices; v++ )
   {
      if( graph->marks[v] == idx )
      {
         size++;
      }
   }
   return size;
}

/** calculate the number of vertices for component idx present within obj function */
static
int sizeObjGraphComponent(
   SCIP*                 scip,
   CONSTRAINTMATRIX*     matrix,
   GRAPH*                graph,
   int                   idx
   )
{
   int v;
   int size;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(idx < graph->ncomponents);

   size = 0;
   for( v = 0; v < graph->nvertices; v++ )
   {
      if( graph->marks[v] == idx &&
         !SCIPisEQ(scip,SCIPvarGetObj(matrix->vars[v]),0.0) )
      {
         size++;
      }
   }
   return size;
}

/** calculate the number of contraints present within component idx */
static
int sizeGraphConstraints(
   SCIP*                 scip,
   CONSTRAINTMATRIX*     matrix,
   GRAPH*                graph,
   int                   idx
   )
{
   int i;
   int size;
   int* rowflags;
   int* colpnt;
   int* colend;

   assert(matrix != NULL);
   assert(graph != NULL);
   assert(idx < graph->ncomponents);

   SCIP_CALL( SCIPallocBufferArray(scip, &rowflags, matrix->nrows) );
   for( i = 0; i < matrix->nrows; i++ )
   {
      rowflags[i] = 0;
   }

   for( i = 0; i < graph->nvertices; i++ )
   {
      if( graph->marks[i] == idx )
      {
         colpnt = matrix->colmatind + matrix->colmatbeg[i];
         colend = colpnt + matrix->colmatcnt[i];
         for(; colpnt < colend; colpnt++ )
         {
            rowflags[*colpnt] = 1;
         }
      }
   }

   size = 0;
   for( i = 0; i < matrix->nrows; i++ )
   {
      if( rowflags[i] == 1 )
      {
         size++;
      }
   }

   SCIPfreeBufferArray(scip, &rowflags);

   return size;
}

/** delete graph data structure */
static
void freeGraph(
   SCIP*                 scip,
   GRAPH**               graph
   )
{
   int v;

   assert(scip != NULL);
   assert(graph != NULL);

   if( (*graph)->nvertices > 0 )
   {
      assert(*graph != NULL);
      assert((*graph)->visited != NULL);
      assert((*graph)->adjlist != NULL);
      assert((*graph)->adjlistinfo != NULL);
      assert((*graph)->marks != NULL);

      for( v = 0; v < (*graph)->nvertices; v++ )
      {
         if( ((*graph)->adjlist[v]) != NULL )
         {
            SCIPfreeMemoryArray(scip, &((*graph)->adjlist[v]));
         }
      }

      SCIPfreeBufferArray(scip, &((*graph)->marks));
      SCIPfreeBufferArray(scip, &((*graph)->adjlistinfo));
      SCIPfreeBufferArray(scip, &((*graph)->adjlist));
      SCIPfreeBufferArray(scip, &((*graph)->visited));

      (*graph)->nvertices = 0;
      (*graph)->nedges = 0;
      (*graph)->ncomponents = 0;
   }

   SCIPfreeBuffer(scip,graph);
}

/*
 * Callback methods of presolver
 */

/* TODO: Implement all necessary presolver methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_PRESOLCOPY(presolCopyConcomp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of small connected components presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolCopyConcomp NULL
#endif


/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeConcomp)
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
SCIP_DECL_PRESOLINIT(presolInitConcomp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of small connected components presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolInitConcomp NULL
#endif


/** deinitialization method of presolver (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_PRESOLEXIT(presolExitConcomp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of small connected components presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolExitConcomp NULL
#endif

/** performs presolving by searching for connected components */
static
SCIP_RETCODE presolComponents(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_PRESOL*          presol,             /**< the presolver itself */
   int*                  nfixedvars,         /**< pointer to total number of variables fixed of all presolvers */
   int*                  ndelconss,          /**< pointer to total number of deleted constraints of all presolvers */
   SCIP_RESULT*          result              /**< pointer to store the result of the presolving call */
   )
{
   CONSTRAINTMATRIX* matrix;
   SCIP_PRESOLDATA* presoldata;
   SCIP_Bool matrixInitialized;
   GRAPH* graph;
   SCIP_Bool graphInitialized;
   int i;

   matrix = NULL;
   matrixInitialized = FALSE;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING || SCIPinProbing(scip) )
      return SCIP_OKAY;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   if( !presoldata->dosearch )
   {
      /* do not search for connected components  */
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPallocBuffer(scip, &matrix) );
   initMatrix(scip, matrix, &matrixInitialized);

   if( matrixInitialized )
   {
      graphInitialized = FALSE;
      SCIP_CALL( SCIPallocBuffer(scip, &graph) );
      initGraph(scip, matrix, graph, &graphInitialized);
      /* writeGraph(scip,matrix,graph); */

      if( graphInitialized )
      {
         calcGraphConComp(graph);
         /* writeConComp(scip,matrix,graph); */

         printf("### GRAPH (array): %d connected components (%d cons, %d vars):\n",
            graph->ncomponents,matrix->nrows,matrix->ncols);

         if( presoldata->showsizes )
         {
            for( i = 0; i < graph->ncomponents; i++ )
            {
               printf("### [CoCo%d]-> %d cons, %d vars (%d vars obj).\n",
                  i,
                  sizeGraphConstraints(scip,matrix,graph,i),
                  sizeGraphComponent(graph, i),
                  sizeObjGraphComponent(scip,matrix,graph,i));
            }
         }
         createSubscips(scip, matrix, graph, presoldata);

      }
      freeGraph(scip, &graph);
   }

   freeMatrix(scip, &matrix);

   return SCIP_OKAY;
}


/** presolving initialization method of presolver (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_PRESOLINITPRE(presolInitpreConcomp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of small connected components presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolInitpreConcomp NULL
#endif


/** presolving deinitialization method of presolver (called after presolving has been finished) */
static
SCIP_DECL_PRESOLEXITPRE(presolExitpreConcomp)
{  /*lint --e{715}*/

   SCIP_CALL( presolComponents(scip, presol, NULL, NULL, result) );

   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecConcomp)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the small connected components presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolConcomp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;

   /* create concomp presolver data */
   SCIP_CALL( SCIPallocMemory(scip, &presoldata) );
   initPresoldata(presoldata);

   /* include presolver */
   SCIP_CALL( SCIPincludePresol(scip,
         PRESOL_NAME,
         PRESOL_DESC,
         PRESOL_PRIORITY,
         PRESOL_MAXROUNDS,
         PRESOL_DELAY,
         presolCopyConcomp,
         presolFreeConcomp,
         presolInitConcomp,
         presolExitConcomp,
         presolInitpreConcomp,
         presolExitpreConcomp,
         presolExecConcomp,
         presoldata) );

   /* add presolver parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/concomp/dosearch",
         "search for connected components (0: no search, 1: do search)",
         &presoldata->dosearch, FALSE, DEFAULT_SEARCH, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/concomp/showsizes",
         "show single sizes (0: no, 1: show)",
         &presoldata->showsizes, FALSE, DEFAULT_SHOW_SIZES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/concomp/writeproblems",
         "should the single components be written as an .lp-file?",
         &presoldata->writeproblems, FALSE, DEFAULT_WRITEPROBLEMS, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/concomp/maxintvars",
         "maximum number of integer (or binary) variables to solve a subproblem directly (-1: no solving)",
         &presoldata->maxintvars, FALSE, DEFAULT_MAXINTVARS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip,
         "presolving/concomp/nodelimit",
         "maximum number of nodes to be solved in subproblems",
         &presoldata->nodelimit, FALSE, DEFAULT_NODELIMIT, -1, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/concomp/intfactor",
         "the weight of an integer variable compared to binary variables",
         &presoldata->intfactor, FALSE, DEFAULT_INTFACTOR, 0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

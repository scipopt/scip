/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   sepa_multilinear.c
 * @ingroup DEFPLUGINS_SEPA
 * @brief  multilinear separator
 * @author Matthias Walter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/sepa_multilinear.h"
#include "scip/cons_and.h"
#include "scip/cons_nonlinear.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/set.h"
#include "scip/hypergraph.h"


#define SEPA_NAME                      "multilinear"
#define SEPA_DESC                      "multilinear cut separator"
#define SEPA_PRIORITY                  100000
#define SEPA_FREQ                      1
#define SEPA_MAXBOUNDDIST              1.0
#define SEPA_USESSUBSCIP               FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                     FALSE /**< should separation method be delayed if other separators found cuts? */

#define DEFAULT_MIN_OVERLAPS           1
#define DEFAULT_SCAN_AND               TRUE
#define DEFAULT_SCAN_PRODUCT           FALSE
#define DEFAULT_MAX_STANDARD           0
#define DEFAULT_MAX_ONEFLOWER          10000000
#define DEFAULT_MAX_TWOFLOWER          10000000
#define DEFAULT_DELAY_STANDARD         FALSE
#define DEFAULT_MAX_USELESS_ONEFLOWER  2 /**< Number of useless separation rounds after which we stop separating. */
#define DEFAULT_MAX_USELESS_TWOFLOWER  2 /**< Number of useless separation rounds after which we stop separating. */

/* Define this for quickly testing whether instances are affected at all. */
// #define PRINT_HYPERGRAPH_AND_EXIT

/* TODO: These old codes shall be removed once the computational comparison with the new ones is published. */
// #define USE_OLD_ONEFLOWER_SEPARATION
// #define USE_OLD_TWOFLOWER_SEPARATION

#ifdef USE_OLD_ONEFLOWER_SEPARATION
#define MAXNSEPA_ONEFLOWER_PER_BASE 2 /* Maximum number of inequalities to be separated that have the same base. */
#endif /* USE_OLD_ONEFLOWER_SEPARATION */

/** data associated with each hypergraph node */
struct SCIP_Hypergraph_NodeData
{
   SCIP_VAR*             var;                /**< Variable associated with this node. */
   SCIP_Real             coefscale;          /**< Factor to scale a cut's coefficient with; equals 1/ub. */
   SCIP_Real             solval;             /**< Solution value in [0,1]-transformed space. */
};

/** data associated with each hypergraph edge */
struct SCIP_Hypergraph_EdgeData
{
   SCIP_VAR*             var;                /**< Variable of the product. */
   SCIP_Real             coefficient;        /**< Coefficient of the product. */
   SCIP_Real             coefscale;          /**< Factor to scale a cut's coefficient with;
                                              **  equals 1/(product of ubs * coefficient). */
   SCIP_Real             solval;             /**< Value of the variable in [0,1]-transformed space. */
   SCIP_Real             slackval;           /**< The slack value of that edge. */
};

/** data associated with each overlap, i.e., intersection of two edges. */
struct SCIP_Hypergraph_OverlapData
{
   SCIP_Real             sumnodecomplements; /**< Sum of 1-z_v for all v in this overlap. */
   SCIP_Real             minedgecomplement;  /**< 1 minus the maximum solution value of any edge incident to this overlap. */
   SCIP_HYPERGRAPH_EDGE  minedge;            /**< An edge for which minedgecomplement is attained. */
};


/** Separator data. */
struct SCIP_SepaData
{
   int                   lastrun;            /**< Last run for which we constructed a hypergraph. */
   SCIP_Bool             scanand;            /**< Whether to scan AND constraints when constructing a hypergraph. */
   SCIP_Bool             scanproduct;        /**< Whether to scan product expressions when constructing a hypergraph. */
   SCIP_HYPERGRAPH*      hypergraph;         /**< The hypergraph. */
   int                   nsepacuts;          /**< Total number of generated cuts. */

   SCIP_Real             timehypercreation;  /**< Total time spent on constructing hypergraphs. */
   SCIP_Real             timehyperoverlaps;  /**< Total time spent on computing hypergraphs' overlaps. */
   SCIP_Real             timepreparation;    /**< Time spent on joint preparation for all separation algorithms. */
   int                   nsepastandard;      /**< Total number of generated standard relaxation inequalities. */
   SCIP_Real             timesepaoneflower;  /**< Total time spent on separation problem for 1-flower inequalities. */
   int                   nsepaoneflower;     /**< Total number of generated 1-flower inequalities. */
   SCIP_Real             timesepatwoflower;  /**< Total time spent on separation problem for 1-flower inequalities. */
   int                   nsepatwoflower;     /**< Total number of generated 2-flower inequalities. */

   int                   minnoverlaps;       /**< Minimum number of overlaps needed to actually try separation. */
   int                   maxstandard;        /**< Maximum number of standard relaxation inequalities per round. */
   int                   maxoneflower;       /**< Maximum number of 1-flower inequalities per round. */
   int                   maxtwoflower;       /**< Maximum number of 2-flower inequalities per round. */
   SCIP_Bool             delaystandard;      /**< Whether to only generate standard inequalities if also flower. */
   int                   maxuselessoneflower;/**< Number of useless separation rounds after which we stop separating. */
   int                   maxuselesstwoflower;/**< Number of useless separation rounds after which we stop separating. */
   int                   nuselessoneflower;  /**< Number of recent useless separation rounds for 1-flowers. */
   int                   nuselesstwoflower;  /**< Number of recent useless separation rounds for 2-flowers. */
};

/*
 * Local methods
 */

/** @brief constructs the hypergraph from transformed problem */
static
SCIP_RETCODE constructHypergraph(
   SCIP*                 scip,               /**< SCIP data structure. */
   SCIP_SEPADATA*        sepadata            /**< Sepadata. */
   )
{
   SCIP_CLOCK* clock;
   int nvars;
   SCIP_HASHMAP* varsnodes = NULL;
   SCIP_HYPERGRAPH_VERTEX* vertices = NULL;
   int nvertices;
   int memvertices = 32;
   SCIP_CONSHDLR* conshdlr;

   assert(scip);
   assert(sepadata);
   assert(sepadata->hypergraph == NULL);

   SCIPdebugMsg(scip, "scanning constraints to construct hypergraph.\n");

   SCIP_CALL( SCIPcreateClock(scip, &clock) );
   SCIP_CALL( SCIPstartClock(scip, clock) );

   nvars = SCIPgetNVars(scip);
   if( nvars == 0 )
      ++nvars;
   SCIP_CALL( SCIPhypergraphCreate(&sepadata->hypergraph, SCIPblkmem(scip), nvars, nvars, nvars, 4,
      sizeof(SCIP_HYPERGRAPH_VERTEXDATA), sizeof(SCIP_HYPERGRAPH_EDGEDATA), sizeof(SCIP_HYPERGRAPH_OVERLAPDATA) ) );
   SCIP_HYPERGRAPH* hypergraph = sepadata->hypergraph;

   /* Create map from variables to nodes and vertex array. */
   SCIP_CALL( SCIPhashmapCreate(&varsnodes, SCIPblkmem(scip), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vertices, memvertices) );

   if( sepadata->scanand )
   {
      /* Scan AND constraints. */
      conshdlr = SCIPfindConshdlr(scip, "and");
      if( conshdlr != NULL )
      {
         SCIP_CONS** conss;
         int nconss;

         conss = SCIPconshdlrGetConss(conshdlr);
         nconss = SCIPconshdlrGetNConss(conshdlr);
         SCIPdebugMsg(scip, "  processing %d AND constraints.\n", nconss);
         for( int i = 0; i < nconss; ++i )
         {
            SCIP_CONS* cons;
            SCIP_VAR** vars;
            SCIP_HYPERGRAPH_EDGE edge;
            SCIP_HYPERGRAPH_EDGEDATA* edgedata = NULL;

            /* Get number of variables. */
            cons = conss[i];
            nvertices = SCIPgetNVarsAnd(scip, cons);
            if( nvertices > memvertices )
            {
               int newcapacity;

               newcapacity = SCIPcalcMemGrowSize(scip, nvertices);
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &vertices, memvertices, newcapacity) );
               memvertices = newcapacity;
            }

            vars = SCIPgetVarsAnd(scip, cons);
            for( int j = 0; j < nvertices; ++j )
            {
               SCIP_HYPERGRAPH_VERTEX v = SCIPhashmapGetImageInt(varsnodes, vars[j]);
               if( v == INT_MAX )
               {
                  SCIP_HYPERGRAPH_VERTEXDATA* vertexdata = NULL;
                  SCIP_CALL( SCIPhypergraphAddVertex(sepadata->hypergraph, &v, &vertexdata) );
                  vertexdata->var = vars[j];
                  SCIP_CALL( SCIPhashmapInsertInt(varsnodes, vars[j], v) );
               }
               vertices[j] = v;
            }
            SCIP_CALL( SCIPhypergraphAddEdge(sepadata->hypergraph, nvertices, vertices, &edge, &edgedata) );
            edgedata->var = SCIPgetResultantAnd(scip, cons);
            edgedata->coefficient = 1.0;
         }
      }
   }

   if( sepadata->scanproduct )
   {
      /* Scan nonlinear constraints for product expressions. */
      conshdlr = SCIPfindConshdlr(scip, "nonlinear");
      if( conshdlr != NULL )
      {
         int nconss;
         SCIP_CONS** conss;
         SCIP_EXPRITER* it = NULL;

         nconss = SCIPconshdlrGetNConss(conshdlr);
         conss = SCIPconshdlrGetConss(conshdlr);
         assert( conss != NULL || (nconss == 0) );
         SCIPdebugMsg(scip, "  processing %d nonlinear constraints.\n", nconss);

         /* Prepare iteration such that we visit every expression only once. */
         SCIP_CALL( SCIPcreateExpriter(scip, &it) );
         SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );
         SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ENTEREXPR);

         for( int c = 0; c < nconss; ++c )
         {
            SCIP_EXPR* expr;

            /* Iterate through all expressions of the nonlinear constraint that we haven't seen so far. */
            expr = SCIPgetExprNonlinear(conss[c]);
            for( expr = SCIPexpriterRestartDFS(it, expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) ) /*lint !e441 *//*lint !e440 */
            {
               if( SCIPisExprProduct(scip, expr) )
               {
                  SCIP_VAR* exprvar;
                  SCIP_EXPR** children;
                  int j;

                  exprvar = SCIPgetExprAuxVarNonlinear(expr);
                  if( exprvar == NULL )
                     continue;

                  children = SCIPexprGetChildren(expr);
                  nvertices = SCIPexprGetNChildren(expr);

                  if( nvertices > memvertices )
                  {
                     int newcapacity;

                     newcapacity = SCIPcalcMemGrowSize(scip, nvertices);
                     SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &vertices, memvertices, newcapacity) );
                     memvertices = newcapacity;
                  }

                  for( j = 0; j < nvertices; ++j )
                  {
                     SCIP_VAR* auxvar;
                     SCIP_HYPERGRAPH_VERTEX v;

                     auxvar = SCIPgetExprAuxVarNonlinear(children[j]);
                     if( auxvar == NULL || SCIPisLT(scip, SCIPvarGetLbGlobal(auxvar), 0.0) )
                        break;
                     v = SCIPhashmapGetImageInt(varsnodes, auxvar);
                     if( v == INT_MAX )
                     {
                        SCIP_HYPERGRAPH_VERTEXDATA* vertexdata = NULL;
                        SCIP_CALL( SCIPhypergraphAddVertex(sepadata->hypergraph, &v, &vertexdata) );
                        vertexdata->var = auxvar;
                        SCIP_CALL( SCIPhashmapInsertInt(varsnodes, auxvar, v) );
                     }
                     vertices[j] = v;
                  }
                  if( j == nvertices )
                  {
                     /* If the lower bound is nonnegative, then we can consider its relaxed domain [0, ub]. */
                     SCIP_HYPERGRAPH_EDGEDATA* edgedata = NULL;
                     SCIP_HYPERGRAPH_EDGE edge;
                     SCIP_CALL( SCIPhypergraphAddEdge(sepadata->hypergraph, nvertices, vertices, &edge, &edgedata) );
                     edgedata->var = exprvar;
                     edgedata->coefficient = SCIPgetCoefExprProduct(expr);
                  }
               }
            }
         }

         SCIPfreeExpriter(&it);
      }
   }

   SCIPfreeBlockMemoryArray(scip, &vertices, memvertices);
   SCIPhashmapFree(&varsnodes);

   /* Compute each node's incident edges. */
   SCIP_CALL( SCIPhypergraphComputeVerticesEdges(hypergraph) );

   assert( SCIPhypergraphIsValid(hypergraph, stdout) );

   SCIP_Real time = SCIPgetClockTime(scip, clock);
   sepadata->timehypercreation += time;

   /* Find all pair-wise intersections of edges (of size at least 2). */
   SCIP_CALL( SCIPhypergraphComputeOverlaps(hypergraph, NULL, NULL) );

   sepadata->timehyperoverlaps += SCIPgetClockTime(scip, clock) - time;
   SCIP_CALL( SCIPfreeClock(scip, &clock) );

   assert( SCIPhypergraphIsValid(hypergraph, stdout) );

   SCIPdebugMsg(scip, "the hypergraph has %d proper overlaps.\n", SCIPhypergraphGetNOverlaps(hypergraph));

   SCIP_CALL( SCIPhypergraphComputeOverlapsEdges(hypergraph) );

#ifdef PRINT_HYPERGRAPH_AND_EXIT
   SCIPmessagePrintInfo(SCIPgetMessagehdlr(scip), "The hypergraph has %d vertices, %d edges and %d overlaps and was "
      "computed in %f+%f=%f seconds.\n", SCIPhypergraphGetNVertices(hypergraph), SCIPhypergraphGetNEdges(hypergraph),
      SCIPhypergraphGetNOverlaps(hypergraph), sepadata->timehypercreation, sepadata->timehyperoverlaps,
      sepadata->timehypercreation + sepadata->timehyperoverlaps);
   return SCIP_ERROR;
#endif /* PRINT_HYPERGRAPH_AND_EXIT */

   return SCIP_OKAY;
}

/* @brief prepare the separation for all cutting plane types by storing relevant data with the hypergraph */
static
SCIP_RETCODE prepareSeparation(
   SCIP*                 scip,               /**< SCIP data structure. */
   SCIP_SEPA*            sepa,               /**< Separator. */
   SCIP_SOL*             sol                 /**< Solution to be separated. */
   )
{
   SCIP_SEPADATA* sepadata;
   int nedges;
   int noverlaps;
   SCIP_Real feastol;
   SCIP_CLOCK* clock = NULL;

   assert(scip);
   assert(sepa);

   feastol = SCIPfeastol(scip);
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata);
   if( sepadata->hypergraph == NULL)
      return SCIP_OKAY;

   SCIP_CALL( SCIPcreateClock(scip, &clock) );
   SCIP_CALL( SCIPstartClock(scip, clock) );

   SCIPdebugMsg(scip, "separating a solution of value %g...\n", sol == NULL ? SCIPgetLPObjval(scip) : SCIPsolGetOrigObj(sol));

   /* We extract the solution values for the hypergraph's vertices. */
   for( SCIP_HYPERGRAPH_VERTEX vertex = 0; vertex < SCIPhypergraphGetNVertices(sepadata->hypergraph); ++vertex )
   {
      SCIP_HYPERGRAPH_VERTEXDATA* vertexdata;
      SCIP_Real value;
      SCIP_Real ub;

      vertexdata = SCIPhypergraphVertexData(sepadata->hypergraph, vertex);
      assert( SCIPisGE(scip, SCIPvarGetLbLocal(vertexdata->var), 0.0) );
      value = SCIPgetSolVal(scip, sol, vertexdata->var);
      ub = SCIPvarGetUbLocal(vertexdata->var);

      /* Compute the scaling factor for scaling the domain [lb,ub] to [0,1]. */
      vertexdata->coefscale = 1.0 / MAX(ub, feastol);
      vertexdata->solval = value * vertexdata->coefscale;
   }

   /* We extract the solution values for the hypergraph's edges. */
   nedges = SCIPhypergraphGetNEdges(sepadata->hypergraph);
   for( SCIP_HYPERGRAPH_EDGE edge = 0; edge < nedges; ++edge )
   {
      SCIP_HYPERGRAPH_EDGEDATA* edgedata;
      SCIP_Real ubprod = 1.0;
      int size;
      SCIP_HYPERGRAPH_VERTEX* vertices;

      edgedata = SCIPhypergraphEdgeData(sepadata->hypergraph, edge);
      size = SCIPhypergraphEdgeSize(sepadata->hypergraph, edge);
      vertices = SCIPhypergraphEdgeVertices(sepadata->hypergraph, edge);
      for( int i = 0; i < size; ++i )
         ubprod /= SCIPhypergraphVertexData(sepadata->hypergraph, vertices[i])->coefscale;
      SCIP_Real value = SCIPgetSolVal(scip, sol, edgedata->var);

      /* Compute the scaling factor for scaling the domain [0, product of ubs * edge-coefficient ] to [0,1]. */
      edgedata->coefscale = 1.0 / (edgedata->coefficient * ubprod);
      edgedata->coefscale = MIN(edgedata->coefscale, 1.0 / feastol );
      edgedata->solval = value * edgedata->coefscale;
      edgedata->slackval = edgedata->solval - 1.0;
      for( int i = 0; i < size; ++i )
         edgedata->slackval += 1.0 - SCIPhypergraphVertexData(sepadata->hypergraph, vertices[i])->solval;
   }

   /* Initialize overlap's sumnodecomplements to \sum_{v in o} (1-z_v).
    * Together with minedgecomplement, which is the minimum value of 1-z_e over all incident edges e,
    * the difference sumnodecomplements - minedgecomplement represents the violation increase if the sum of (1-z_v)
    * is replaced by (1-z_e) in k-flower inequalities.
    */
   noverlaps = SCIPhypergraphGetNOverlaps(sepadata->hypergraph);
   for( SCIP_HYPERGRAPH_OVERLAP overlap = 0; overlap < noverlaps; ++overlap )
   {
      SCIP_HYPERGRAPH_OVERLAPDATA* overlapdata;
      int nvertices;
      SCIP_HYPERGRAPH_VERTEX* vertices;

      overlapdata = SCIPhypergraphOverlapData(sepadata->hypergraph, overlap);
      overlapdata->minedgecomplement = 2.0;
      overlapdata->sumnodecomplements = 0.0;
      nvertices = SCIPhypergraphOverlapSize(sepadata->hypergraph, overlap);
      vertices = SCIPhypergraphOverlapVertices(sepadata->hypergraph, overlap);
      for( int i = 0; i < nvertices; ++i )
      {
         SCIP_HYPERGRAPH_VERTEXDATA* vertexdata;

         vertexdata = SCIPhypergraphVertexData(sepadata->hypergraph, vertices[i]);
         overlapdata->sumnodecomplements += 1.0 - vertexdata->solval;
      }
   }

   /* Iterate over all edge-overlap incidences to compute each overlap's maximum-value edge.
    * Note that the initial value of minedgecomplement is guaranteed to be larger than 1 - z_e for any incident edge e,
    * so it will be updated in the loop below. */
   for( SCIP_HYPERGRAPH_EDGE edge = 0; edge < SCIPhypergraphGetNEdges(sepadata->hypergraph); ++edge )
   {
      SCIP_HYPERGRAPH_EDGEDATA* edgedata;
      int beyond;

      edgedata = SCIPhypergraphEdgeData(sepadata->hypergraph, edge);
      beyond = SCIPhypergraphEdgesOverlapsBeyond(sepadata->hypergraph, edge);
      for( int i = SCIPhypergraphEdgesOverlapsFirst(sepadata->hypergraph, edge); i < beyond; ++i )
      {
         SCIP_HYPERGRAPH_OVERLAP overlap;
         SCIP_HYPERGRAPH_OVERLAPDATA* overlapdata;
         SCIP_Real edgecomplement;

         overlap = SCIPhypergraphEdgesOverlapsGetAtIndex(sepadata->hypergraph, i);
         overlapdata = SCIPhypergraphOverlapData(sepadata->hypergraph, overlap);
         edgecomplement = 1.0 - edgedata->solval;
         if( edgecomplement < overlapdata->minedgecomplement )
         {
            overlapdata->minedgecomplement = edgecomplement;
            overlapdata->minedge = edge;
         }
      }
   }

   sepadata->timepreparation += SCIPgetClockTime(scip, clock);
   SCIP_CALL( SCIPfreeClock(scip, &clock) );

   return SCIP_OKAY;
}

/**
 * @brief add a generated cut row to the cut pool (for the root node) or as a row (otherwise)
 */
static
SCIP_RETCODE addCut(
   SCIP*                 scip,               /**< SCIP datastructure. */
   SCIP_SOL*             sol,                /**< Solution. */
   SCIP_ROW*             row,                /**< Cutting plane. */
   int*                  pnumseparated,      /**< Pointer to store number of separated cuts. */
   SCIP_RESULT*          presult,            /**< Pointer to store result. */
   SCIP_Bool*            padded              /**< Pointer for storing whether it was added. */
   )
{
   assert(scip);
   assert(row);
   assert(pnumseparated);
   assert(presult);

   if( SCIPisCutEfficacious(scip, sol, row) )
   {
      if( SCIPgetDepth(scip) == 0 )
      {
         SCIP_CALL( SCIPaddPoolCut(scip, row) );
      }
      else
      {
         SCIP_Bool infeasible;
         SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
         if( infeasible )
            *presult = SCIP_CUTOFF;
      }
      if( *presult != SCIP_CUTOFF )
         *presult = SCIP_SEPARATED;
      (*pnumseparated)++;
      if( padded )
         *padded = TRUE;
   }
   else if( padded )
      *padded = FALSE;

   return SCIP_OKAY;
}

/**
 * @brief separate missing inequalities from the standard relaxation
 */
static
SCIP_RETCODE separateStandard(
   SCIP*                 scip,               /**< SCIP datastructure. */
   SCIP_SEPA*            sepa,               /**< Separator. */
   SCIP_SOL*             sol,                /**< Solution to be separated. */
   int                   maxnsepa,           /**< Maximum number of separated inequalities. */
   SCIP_RESULT*          presult             /**< Pointer to store result. */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_CLOCK* clock = NULL;
   int nseparated = 0;

   assert(scip);
   assert(sepa);
   assert(presult);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata);

   if( sepadata->hypergraph == NULL)
      return SCIP_OKAY;

#ifdef SCIP_DEBUG
   SCIPdebugMessage("  ...standard relaxation inequalities: ");
   fflush(stdout);
#endif

   SCIP_CALL( SCIPcreateClock(scip, &clock) );
   SCIP_CALL( SCIPstartClock(scip, clock) );

   for( SCIP_HYPERGRAPH_EDGE edge = 0; edge < SCIPhypergraphGetNEdges(sepadata->hypergraph) && *presult != SCIP_CUTOFF
      && nseparated < maxnsepa; ++edge )
   {
      SCIP_HYPERGRAPH_EDGEDATA* edgedata;
      SCIP_HYPERGRAPH_VERTEX* vertices;
      int size;

      edgedata = SCIPhypergraphEdgeData(sepadata->hypergraph, edge);
      vertices = SCIPhypergraphEdgeVertices(sepadata->hypergraph, edge);
      size = SCIPhypergraphEdgeSize(sepadata->hypergraph, edge);
      for( int i = 0; i < size && *presult != SCIP_CUTOFF && nseparated < maxnsepa; ++i )
      {
         SCIP_HYPERGRAPH_VERTEXDATA* vertexdata;

         vertexdata = SCIPhypergraphVertexData(sepadata->hypergraph, vertices[i]);
         if( SCIPisEfficacious(scip, edgedata->solval - vertexdata->solval) )
         {
            SCIP_VAR* vars[2];
            SCIP_Real coefs[2];
            char name[SCIP_MAXSTRLEN];
            SCIP_ROW* row = NULL;

            /* 0 <= z_v - z_e */

            vars[0] = vertexdata->var;
            vars[1] = edgedata->var;
            coefs[0] = 1.0 * vertexdata->coefscale;
            coefs[1] = -1.0 * edgedata->coefscale;
            SCIPsnprintf(name, SCIP_MAXSTRLEN, "multilinear_%05d_standard", ++sepadata->nsepacuts); /*lint !e534 */
            SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &row, sepa, name, 0.0, SCIPinfinity(scip), SCIPgetDepth(scip) > 0,
               FALSE, TRUE) );
            SCIP_CALL( SCIPaddVarsToRow( scip, row, 2, vars, coefs) );

            SCIP_CALL( addCut(scip, sol, row, &nseparated, presult, NULL) );
            SCIP_CALL( SCIPreleaseRow(scip, &row) );
         }
      }
   }

   sepadata->nsepastandard += nseparated;

#ifdef SCIP_DEBUG
   printf("found %d in %.3fs.\n", nseparated, SCIPgetClockTime(scip, clock));
   fflush(stdout);
#endif

   SCIP_CALL( SCIPfreeClock(scip, &clock) );

   return SCIP_OKAY;
}

#ifdef USE_OLD_ONEFLOWER_SEPARATION

/* Separate 1-flower inequalities. */
static
SCIP_RETCODE separateOneFlowerOld(
   SCIP*                 scip,               /**< SCIP data structure. */
   SCIP_SEPA*            sepa,               /**< Separator. */
   SCIP_SOL*             sol,                /**< Solution to be separated. */
   int                   maxnsepa,           /**< Maximum number of separated inequalities. */
   SCIP_RESULT*          presult             /**< Pointer to store result. */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_HYPERGRAPH* hypergraph;
   SCIP_CLOCK* clock = NULL;
   int nseparated = 0;
   SCIP_HYPERGRAPH_ITER iter;

   assert(scip);
   assert(sepa);
   assert(presult);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata);

   hypergraph = sepadata->hypergraph;
   if( hypergraph == NULL )
      return SCIP_OKAY;

#ifdef SCIP_DEBUG
   SCIPdebugMessage("  ...flower inequalities with 1 neighbor: ");
   fflush(stdout);
#endif

   SCIP_CALL( SCIPcreateClock(scip, &clock) );
   SCIP_CALL( SCIPstartClock(scip, clock) );

   SCIP_CALL( SCIPhypergraphIterInit(hypergraph, &iter) );
   for( SCIP_HYPERGRAPH_EDGE base = 0; base < SCIPhypergraphGetNEdges(hypergraph) && *presult != SCIP_CUTOFF
      && nseparated < maxnsepa; ++base )
   {
      SCIP_HYPERGRAPH_EDGEDATA* basedata;
      SCIP_HYPERGRAPH_VERTEX* basevertices;
      int nbasevertices;
      int nbaseseparated = 0;  /**< Number of cuts separated for this base edge. */

      basedata = SCIPhypergraphEdgeData(hypergraph, base);

      /* If the slack value is at least 1 then the inequality can never be violated. */
      if( SCIPisFeasGE(scip, basedata->slackval, 1.0) )
         continue;

      basevertices = SCIPhypergraphEdgeVertices(hypergraph, base);
      nbasevertices = SCIPhypergraphEdgeSize(hypergraph, base);
      for( SCIPhypergraphIterStart(hypergraph, &iter, base, 2, TRUE, FALSE);
         SCIPhypergraphIterValid(&iter); SCIPhypergraphIterNext(hypergraph, &iter) )
      {
         SCIP_HYPERGRAPH_EDGEDATA* adjacentdata;
         SCIP_Real gain;

         /* The standard inequality has slack stdslack:         x_e - 1 + \sum_{v \in e} (1-x_v).
          * The 1-flower with base e and neighbor f has slack:  x_e - 1 + \sum_{v \in e \setminus f} (1-x_v) + 1-x_f.
          * The gain (difference standard - 1-flower) is thus: \sum_{v \in e \cap f} (1-x_v) - 1 + x_f
          * If gain > stdslack then the 1-flower inequality has negative slack.
          */

         adjacentdata = SCIPhypergraphEdgeData(hypergraph, iter.adjacent);
         gain = adjacentdata->solval - 1.0;
         for( int i = 0; i < iter.ncommonvertices; ++i )
         {
            gain += 1.0 - SCIPhypergraphVertexData(hypergraph, iter.commonvertices[i])->solval;
         }

         if( SCIPisEfficacious(scip, gain - basedata->slackval) )
         {
            SCIP_ROW* row = NULL;
            char name[SCIP_MAXSTRLEN];

            SCIPsnprintf(name, SCIP_MAXSTRLEN, "multilinear_%05d_1flower", ++sepadata->nsepacuts);
            SCIP_CALL( SCIPcreateRowSepa(scip, &row, sepa, name, 0, NULL, NULL, iter.ncommonvertices - nbasevertices,
               SCIPinfinity(scip), SCIPgetDepth(scip) > 0, FALSE, TRUE) );
            SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
            SCIP_CALL( SCIPaddVarToRow(scip, row, basedata->var, 1.0 * basedata->coefscale ) );
            SCIP_CALL( SCIPaddVarToRow(scip, row, adjacentdata->var, -1.0  * adjacentdata->coefscale) );
            for( int i = 0; i < nbasevertices; ++i )
            {
               SCIP_HYPERGRAPH_VERTEXDATA* vertexdata;

               vertexdata = SCIPhypergraphVertexData(hypergraph, basevertices[i]);
               SCIP_CALL( SCIPaddVarToRow(scip, row, vertexdata->var , -1.0 * vertexdata->coefscale) );
            }
            for( int i = 0; i < iter.ncommonvertices; ++i )
            {
               SCIP_HYPERGRAPH_VERTEXDATA* vertexdata;

               vertexdata = SCIPhypergraphVertexData(hypergraph, iter.commonvertices[i]);
               SCIP_CALL( SCIPaddVarToRow(scip, row, vertexdata->var , +1.0 * vertexdata->coefscale ) );
            }
            SCIP_CALL( SCIPflushRowExtensions(scip, row) );
            SCIP_CALL( addCut(scip, sol, row, &nseparated, presult, NULL) );
            SCIP_CALL( SCIPreleaseRow(scip, &row) );

            if( nbaseseparated > MAXNSEPA_ONEFLOWER_PER_BASE )
               break;
         }
      }
   }

   SCIPhypergraphIterClear(hypergraph, &iter);
   sepadata->nsepaoneflower += nseparated;

#ifdef SCIP_DEBUG
   printf("found %d in %.3fs.\n", nseparated, SCIPgetClockTime(scip, clock));
   fflush(stdout);
#endif

   sepadata->timesepaoneflower += SCIPgetClockTime(scip, clock);
   SCIP_CALL( SCIPfreeClock(scip, &clock) );

   return SCIP_OKAY;
}

#else /* !USE_OLD_ONEFLOWER_SEPARATION */

/* Separate 1-flower inequalities. */
static
SCIP_RETCODE separateOneFlower(
   SCIP*                 scip,               /**< SCIP data structure. */
   SCIP_SEPA*            sepa,               /**< Separator. */
   SCIP_SOL*             sol,                /**< Solution to be separated. */
   int                   maxnsepa,           /**< Maximum number of separated inequalities. */
   SCIP_RESULT*          presult             /**< Pointer to store result. */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_HYPERGRAPH* hypergraph;
   SCIP_CLOCK* clock = NULL;
   int nseparated = 0;

   assert(scip);
   assert(sepa);
   assert(presult);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata);

   hypergraph = sepadata->hypergraph;
   if( hypergraph == NULL)
      return SCIP_OKAY;

#ifdef SCIP_DEBUG
   SCIPdebugMessage("  ...flower inequalities with 1 neighbor: ");
   fflush(stdout);
#endif

   SCIP_CALL( SCIPcreateClock(scip, &clock) );
   SCIP_CALL( SCIPstartClock(scip, clock) );

   for( SCIP_HYPERGRAPH_EDGE base = 0; (base < SCIPhypergraphGetNEdges(hypergraph)) && (nseparated < maxnsepa); ++base )
   {
      SCIP_HYPERGRAPH_EDGEDATA* basedata;
      int i;
      int beyond;
      SCIP_Real bestgain = -1.0;
      SCIP_HYPERGRAPH_OVERLAP bestoverlap = -1;

      basedata = SCIPhypergraphEdgeData(hypergraph, base);
      i = SCIPhypergraphEdgesOverlapsFirst(hypergraph, base);
      beyond = SCIPhypergraphEdgesOverlapsBeyond(hypergraph, base);

      /*
       * We need to test if z_base + (1-z_overlap) + sum_{v in base \ overlap} (1-z_v) < 1.
       * This is done by minimizing the change (1-z_overlap) - sum_{v in overlap) (1-z)
       * w.r.t. the standard inequality over all overlaps that are incident to base.
       */

      for( ; i < beyond; ++i )
      {
         SCIP_HYPERGRAPH_OVERLAP overlap;
         SCIP_HYPERGRAPH_OVERLAPDATA* overlapdata;
         SCIP_Real gain;

         overlap = SCIPhypergraphEdgesOverlapsGetAtIndex(hypergraph, i);
         overlapdata = SCIPhypergraphOverlapData(hypergraph, overlap);
         gain = overlapdata->sumnodecomplements - overlapdata->minedgecomplement;
         if( gain > bestgain )
         {
            bestgain = gain;
            bestoverlap = overlap;
         }
      }

      if( bestgain > 0 && SCIPisEfficacious(scip, bestgain - basedata->slackval) )
      {
         SCIP_ROW* row = NULL;
         char name[SCIP_MAXSTRLEN];
         SCIP_HYPERGRAPH_VERTEX* basevertices;
         int nbasevertices;
         SCIP_HYPERGRAPH_VERTEX* overlapvertices;
         int noverlapvertices;
         SCIP_HYPERGRAPH_EDGE adjacent;
         SCIP_HYPERGRAPH_EDGEDATA* adjacentdata;

         basevertices = SCIPhypergraphEdgeVertices(hypergraph, base);
         nbasevertices = SCIPhypergraphEdgeSize(hypergraph, base);
         overlapvertices = SCIPhypergraphOverlapVertices(hypergraph, bestoverlap);
         noverlapvertices = SCIPhypergraphOverlapSize(hypergraph, bestoverlap);
         adjacent = SCIPhypergraphOverlapData(hypergraph, bestoverlap)->minedge;
         adjacentdata = SCIPhypergraphEdgeData(hypergraph, adjacent);
         SCIPsnprintf(name, SCIP_MAXSTRLEN, "multilinear_%05d_1flower", ++sepadata->nsepacuts); /*lint !e534*/
         SCIP_CALL( SCIPcreateRowSepa(scip, &row, sepa, name, 0, NULL, NULL, 1.0 * (noverlapvertices - nbasevertices),
            SCIPinfinity(scip), SCIPgetDepth(scip) > 0, FALSE, TRUE) );
         SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
         SCIP_CALL( SCIPaddVarToRow(scip, row, basedata->var, 1.0 * basedata->coefscale ) );
         SCIP_CALL( SCIPaddVarToRow(scip, row, adjacentdata->var, -1.0  * adjacentdata->coefscale) );
         for( i = 0; i < nbasevertices; ++i )
         {
            SCIP_HYPERGRAPH_VERTEXDATA* vertexdata;

            vertexdata = SCIPhypergraphVertexData(hypergraph, basevertices[i]);
            SCIP_CALL( SCIPaddVarToRow(scip, row, vertexdata->var , -1.0 * vertexdata->coefscale) );
         }
         for( i = 0; i < noverlapvertices; ++i )
         {
            SCIP_HYPERGRAPH_VERTEXDATA* vertexdata;

            vertexdata = SCIPhypergraphVertexData(hypergraph, overlapvertices[i]);
            SCIP_CALL( SCIPaddVarToRow(scip, row, vertexdata->var , +1.0 * vertexdata->coefscale ) );
         }
         SCIP_CALL( SCIPflushRowExtensions(scip, row) );

         SCIP_CALL( addCut(scip, sol, row, &nseparated, presult, NULL) );
         SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }
   }

   sepadata->nsepaoneflower += nseparated;

#ifdef SCIP_DEBUG
   printf("found %d in %.3fs.\n", nseparated, SCIPgetClockTime(scip, clock));
   fflush(stdout);
#endif

   sepadata->timesepaoneflower += SCIPgetClockTime(scip, clock);
   SCIP_CALL( SCIPfreeClock(scip, &clock) );

   return SCIP_OKAY;
}

#endif /* USE_OLD_ONEFLOWER_SEPARATION */


#ifdef USE_OLD_TWOFLOWER_SEPARATION

/* Separate 2-flower inequalities. */
static
SCIP_RETCODE separateTwoFlowerOld(
   SCIP*                 scip,               /**< SCIP data structure. */
   SCIP_SEPA*            sepa,               /**< Separator. */
   SCIP_SOL*             sol,                /**< Solution to be separated. */
   int                   maxnsepa,           /**< Maximum number of separated inequalities. */
   SCIP_RESULT*          presult             /**< Pointer to store result. */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_HYPERGRAPH* hypergraph;
   SCIP_CLOCK* clock = NULL;
   int nedges;
   int nseparated = 0;
   int* markedvertices = NULL;
   SCIP_HYPERGRAPH_ITER iter1;
   SCIP_HYPERGRAPH_ITER iter2;

   assert(scip);
   assert(sepa);
   assert(presult);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata);

   hypergraph = sepadata->hypergraph;
   if( hypergraph == NULL)
      return SCIP_OKAY;

#ifdef SCIP_DEBUG
   SCIPdebugMessage("  ...flower inequalities with 2 neighbors: ");
   fflush(stdout);
#endif

   SCIP_CALL( SCIPcreateClock(scip, &clock) );
   SCIP_CALL( SCIPstartClock(scip, clock) );

   nedges = SCIPhypergraphGetNEdges(hypergraph);
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &markedvertices, SCIPhypergraphGetNVertices(hypergraph)) );
   SCIP_CALL( SCIPhypergraphIterInit(hypergraph, &iter1) );
   SCIP_CALL( SCIPhypergraphIterInit(hypergraph, &iter2) );
   for( SCIP_HYPERGRAPH_EDGE base = 0; (base < nedges) && (*presult != SCIP_CUTOFF) && (nseparated < maxnsepa); ++base )
   {
      SCIP_HYPERGRAPH_EDGEDATA* basedata;
      SCIP_HYPERGRAPH_VERTEX* basevertices;
      int nbasevertices;
      int nbaseseparated = 0;

      basedata = SCIPhypergraphEdgeData(hypergraph, base);

      /* If the slack value is at least 1 then the inequality can never be violated. */
      if( SCIPisFeasGE(scip, basedata->slackval, 1.0) )
         continue;

      basevertices = SCIPhypergraphEdgeVertices(hypergraph, base);
      nbasevertices = SCIPhypergraphEdgeSize(hypergraph, base);
      for( SCIPhypergraphIterStart(hypergraph, &iter1, base, 2, FALSE, FALSE); SCIPhypergraphIterValid(&iter1);
         SCIPhypergraphIterNext(hypergraph, &iter1) )
      {
         SCIP_HYPERGRAPH_EDGEDATA* adjacent1data;
         SCIP_Real gain1;

         adjacent1data = SCIPhypergraphEdgeData(hypergraph, iter1.adjacent);
         gain1 = adjacent1data->solval - 1.0;
         for( int i = 0; i < iter1.ncommonvertices; ++i )
         {
            SCIP_HYPERGRAPH_VERTEX v;

            v = iter1.commonvertices[i];
            gain1 += 1.0 - SCIPhypergraphVertexData(hypergraph, v)->solval;
            markedvertices[v] = 1;
         }

         for( SCIPhypergraphIterStart(hypergraph, &iter2, base, 2, FALSE, FALSE); SCIPhypergraphIterValid(&iter2);
            SCIPhypergraphIterNext(hypergraph, &iter2) )
         {
            SCIP_Bool disjoint;
            SCIP_HYPERGRAPH_EDGEDATA* adjacent2data;
            SCIP_Real gain2;

            /* Ordering of neighbor edges. */
            if( iter2.adjacent <= iter1.adjacent )
               continue;
            
            disjoint = TRUE;
            for( int i = 0; i < iter2.ncommonvertices; ++i)
            {
               if( markedvertices[iter2.commonvertices[i]] )
               {
                  disjoint = FALSE;
                  break;
               }
            }
            if( !disjoint )
               continue;
            
            /* The standard inequality has slack: x_e - 1 + \sum_{v \in e} (1-x_v).
             * The 2-flower with base e and neighbors f,g has slack:
             *    x_e - 1 + \sum_{v \in e \setminus (f \cup g)} (1-x_v) + 1-x_f + 1 - x_g.
             * The gain (difference standard - 2-flower) is thus: \sum_{v \in e \cap (f \cup g)} (1-x_v) + x_f + x_g - 2
             * If gain > stdslack then the 2-flower inequality has negative slack.
             */

            adjacent2data = SCIPhypergraphEdgeData(hypergraph, iter2.adjacent);
            gain2 = gain1 + adjacent2data->solval - 1.0;
            for( int j = 0; j < iter2.ncommonvertices; ++j )
            {
               SCIP_HYPERGRAPH_VERTEX w;

               w = iter2.commonvertices[j];
               gain2 += 1.0 - SCIPhypergraphVertexData(hypergraph, w)->solval;
            }

            if( SCIPisEfficacious(scip, gain2 - basedata->slackval) )
            {
               SCIP_ROW* row = NULL;
               char name[SCIP_MAXSTRLEN];
               SCIP_Bool separated;

               SCIPsnprintf(name, SCIP_MAXSTRLEN, "multilinear_%05d_2flower", ++sepadata->nsepacuts);
               SCIP_CALL( SCIPcreateRowSepa(scip, &row, sepa, name, 0, NULL, NULL,
                     iter1.ncommonvertices + iter2.ncommonvertices - nbasevertices - 1.0, SCIPinfinity(scip),
                     SCIPgetDepth(scip) > 0, FALSE, TRUE) );
               SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
               SCIP_CALL( SCIPaddVarToRow(scip, row, basedata->var, 1.0 * basedata->coefscale) );
               SCIP_CALL( SCIPaddVarToRow(scip, row, adjacent1data->var, -1.0 * adjacent1data->coefscale) );
               SCIP_CALL( SCIPaddVarToRow(scip, row, adjacent2data->var, -1.0 * adjacent2data->coefscale) );
               for( int i = 0; i < nbasevertices; ++i )
               {
                  SCIP_HYPERGRAPH_VERTEXDATA* vertexdata;

                  vertexdata = SCIPhypergraphVertexData(hypergraph, basevertices[i]);
                  SCIP_CALL( SCIPaddVarToRow(scip, row, vertexdata->var, -1.0 * vertexdata->coefscale) );
               }
               for( int i = 0; i < iter1.ncommonvertices; ++i )
               {
                  SCIP_HYPERGRAPH_VERTEXDATA* vertexdata;

                  vertexdata = SCIPhypergraphVertexData(hypergraph, iter1.commonvertices[i]);
                  SCIP_CALL( SCIPaddVarToRow(scip, row, vertexdata->var, +1.0 * vertexdata->coefscale) );
               }
               for( int j = 0; j < iter2.ncommonvertices; ++j )
               {
                  SCIP_HYPERGRAPH_VERTEXDATA* vertexdata;

                  vertexdata = SCIPhypergraphVertexData(hypergraph, iter2.commonvertices[j]);
                  SCIP_CALL( SCIPaddVarToRow(scip, row, vertexdata->var, +1.0 * vertexdata->coefscale) );
               }
               SCIP_CALL( SCIPflushRowExtensions(scip, row) );
               SCIP_CALL( addCut(scip, sol, row, &nseparated, presult, &separated) );
               if( separated )
                  ++nbaseseparated;
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
               if( nbaseseparated > 2 )
                  break;
            }
         }

         for( int i = 0; i < iter1.ncommonvertices; ++i )
            markedvertices[iter1.commonvertices[i]] = 0;
      }
   }

   SCIPhypergraphIterClear(hypergraph, &iter2);
   SCIPhypergraphIterClear(hypergraph, &iter1);
   SCIPfreeCleanBuffer(scip, &markedvertices );
   sepadata->nsepatwoflower += nseparated;

#ifdef SCIP_DEBUG
   printf("found %d in %.3fs.\n", nseparated, SCIPgetClockTime(scip, clock));
   fflush(stdout);
#endif

   sepadata->timesepatwoflower += SCIPgetClockTime(scip, clock);
   SCIP_CALL( SCIPfreeClock(scip, &clock) );

   return SCIP_OKAY;
}

#else /* !USE_OLD_TWOFLOWER_SEPARATION */

/* Separate 2-flower inequalities. */
static
SCIP_RETCODE separateTwoFlower(
   SCIP*                 scip,               /**< SCIP datastructure. */
   SCIP_SEPA*            sepa,               /**< Separator. */
   SCIP_SOL*             sol,                /**< Solution to be separated. */
   int                   maxnsepa,           /**< Maximum number of separated inequalities. */
   SCIP_RESULT*          presult             /**< Pointer to store result. */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_HYPERGRAPH* hypergraph;
   SCIP_CLOCK* clock = NULL;
   SCIP_HYPERGRAPH_EDGE base;
   int nseparated = 0;

   assert(scip);
   assert(sepa);
   assert(presult);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata);

   hypergraph = sepadata->hypergraph;
   if( hypergraph == NULL)
      return SCIP_OKAY;

#ifdef SCIP_DEBUG
   SCIPdebugMessage("  ...flower inequalities with 2 neighbors: ");
   fflush(stdout);
#endif

   SCIP_CALL( SCIPcreateClock(scip, &clock) );
   SCIP_CALL( SCIPstartClock(scip, clock) );

   for( base = 0; (base < SCIPhypergraphGetNEdges(hypergraph)) && (nseparated < maxnsepa); ++base )
   {
      SCIP_HYPERGRAPH_EDGEDATA* basedata;
      int basesize;
      int first;
      int beyond;
      int i1;
      int i2;
      SCIP_HYPERGRAPH_OVERLAP bestoverlap1 = -1;
      SCIP_HYPERGRAPH_OVERLAP bestoverlap2 = -1;
      SCIP_Real bestgain = -1.0;

      basedata = SCIPhypergraphEdgeData(hypergraph, base);
      basesize = SCIPhypergraphEdgeSize(hypergraph, base);
      first = SCIPhypergraphEdgesOverlapsFirst(hypergraph, base);
      beyond = SCIPhypergraphEdgesOverlapsBeyond(hypergraph, base);

      /*
       * We need to test if
       *   z_base + (1-z_overlap1) + (1-z_overlap2) + sum_{v in base \ (overlap1 \cup overlap2)} (1-z_v) < 1.
       * This is done by minimizing the change (1-z_overlap1) - sum_{v in overlap1) (1-z) and
       * of (1-z_overlap2) - sum_{v in overlap2) (1-z), both w.r.t. the standard inequality over all overlaps that are
       * incident to base. However, this only correctly reflects the changes violation if the two overlaps are disjoint.
       */

      for( i1 = first; i1 < beyond; ++i1 )
      {
         SCIP_HYPERGRAPH_OVERLAP overlap1;
         SCIP_HYPERGRAPH_OVERLAPDATA* overlap1data;
         SCIP_Real gain1;
         int overlap1size;

         overlap1 = SCIPhypergraphEdgesOverlapsGetAtIndex(hypergraph, i1);
         overlap1size = SCIPhypergraphOverlapSize(hypergraph, overlap1);

         /* If |o1| > |e|/2 then for sure |o2| > holds for sure, so we can stop. */
         if( overlap1size > basesize / 2)
            break;

         overlap1data = SCIPhypergraphOverlapData(hypergraph, overlap1);
         gain1 = overlap1data->sumnodecomplements - overlap1data->minedgecomplement;

         for( i2 = i1 + 1; i2 < beyond; ++i2 )
         {
            SCIP_HYPERGRAPH_OVERLAP overlap2;
            int overlap2size;

            overlap2 = SCIPhypergraphEdgesOverlapsGetAtIndex(hypergraph, i2);
            overlap2size = SCIPhypergraphOverlapSize(hypergraph, overlap2);
            if( (overlap1size + overlap2size > basesize)
               || ((overlap1size == overlap2size) && (overlap1 >= overlap2))
               || !SCIPhypergraphOverlapsDisjoint(hypergraph, overlap1, overlap2) )
            {
               continue;
            }

            SCIP_HYPERGRAPH_OVERLAPDATA* overlap2data;
            SCIP_Real gain2;
            SCIP_Real totalgain;

            overlap2data = SCIPhypergraphOverlapData(hypergraph, overlap2);
            gain2 = overlap2data->sumnodecomplements - overlap2data->minedgecomplement;
            totalgain = gain1 + gain2;

            if( totalgain > bestgain )
            {
               bestgain = totalgain;
               bestoverlap1 = overlap1;
               bestoverlap2 = overlap2;
            }
         }
      }

      if( bestgain > 0 && SCIPisEfficacious(scip, bestgain - basedata->slackval) )
      {
         SCIP_HYPERGRAPH_VERTEX* basevertices;
         int nbasevertices;
         SCIP_HYPERGRAPH_VERTEX* overlap1vertices;
         int noverlap1vertices;
         SCIP_HYPERGRAPH_EDGE adjacent1;
         SCIP_HYPERGRAPH_EDGEDATA* adjacent1data;
         SCIP_HYPERGRAPH_VERTEX* overlap2vertices;
         int noverlap2vertices;
         SCIP_HYPERGRAPH_EDGE adjacent2;
         SCIP_HYPERGRAPH_EDGEDATA* adjacent2data;
         char name[SCIP_MAXSTRLEN];
         SCIP_ROW* row = NULL;

         basevertices = SCIPhypergraphEdgeVertices(sepadata->hypergraph, base);
         nbasevertices = SCIPhypergraphEdgeSize(sepadata->hypergraph, base);
         overlap1vertices = SCIPhypergraphOverlapVertices(sepadata->hypergraph, bestoverlap1);
         noverlap1vertices = SCIPhypergraphOverlapSize(sepadata->hypergraph, bestoverlap1);
         adjacent1 = SCIPhypergraphOverlapData(sepadata->hypergraph, bestoverlap1)->minedge;
         adjacent1data = SCIPhypergraphEdgeData(sepadata->hypergraph, adjacent1);
         overlap2vertices = SCIPhypergraphOverlapVertices(sepadata->hypergraph, bestoverlap2);
         noverlap2vertices = SCIPhypergraphOverlapSize(sepadata->hypergraph, bestoverlap2);
         adjacent2 = SCIPhypergraphOverlapData(sepadata->hypergraph, bestoverlap2)->minedge;
         adjacent2data = SCIPhypergraphEdgeData(sepadata->hypergraph, adjacent2);

         SCIPsnprintf(name, SCIP_MAXSTRLEN, "multilinear_%05d_2flower", ++sepadata->nsepacuts); /*lint !e534*/
         SCIP_CALL( SCIPcreateRowSepa(scip, &row, sepa, name, 0, NULL, NULL,
               noverlap1vertices + noverlap2vertices - nbasevertices - 1.0, SCIPinfinity(scip), SCIPgetDepth(scip) > 0,
               FALSE, TRUE) );
         SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
         SCIP_CALL( SCIPaddVarToRow(scip, row, basedata->var, 1.0 * basedata->coefscale ) );
         SCIP_CALL( SCIPaddVarToRow(scip, row, adjacent1data->var, -1.0  * adjacent1data->coefscale) );
         SCIP_CALL( SCIPaddVarToRow(scip, row, adjacent2data->var, -1.0  * adjacent2data->coefscale) );
         for( int i = 0; i < nbasevertices; ++i )
         {
            SCIP_HYPERGRAPH_VERTEXDATA* data;

            data = SCIPhypergraphVertexData(sepadata->hypergraph, basevertices[i]);
            SCIP_CALL( SCIPaddVarToRow(scip, row, data->var , -1.0 * data->coefscale) );
         }
         for( int i = 0; i < noverlap1vertices; ++i )
         {
            SCIP_HYPERGRAPH_VERTEXDATA* data;

            data = SCIPhypergraphVertexData(sepadata->hypergraph, overlap1vertices[i]);
            SCIP_CALL( SCIPaddVarToRow(scip, row, data->var , +1.0 * data->coefscale ) );
         }
         for( int i = 0; i < noverlap2vertices; ++i )
         {
            SCIP_HYPERGRAPH_VERTEXDATA* data;

            data = SCIPhypergraphVertexData(sepadata->hypergraph, overlap2vertices[i]);
            SCIP_CALL( SCIPaddVarToRow(scip, row, data->var , +1.0 * data->coefscale ) );
         }
         SCIP_CALL( SCIPflushRowExtensions(scip, row) );
         SCIP_CALL( addCut(scip, sol, row, &nseparated, presult, NULL) );
         SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }
   }

   sepadata->nsepatwoflower += nseparated;

#ifdef SCIP_DEBUG
   printf("found %d in %.3fs.\n", nseparated, SCIPgetClockTime(scip, clock));
   fflush(stdout);
#endif

   sepadata->timesepatwoflower += SCIPgetClockTime(scip, clock);
   SCIP_CALL( SCIPfreeClock(scip, &clock) );

   return SCIP_OKAY;
}

#endif /* USE_OLD_TWOFLOWER_SEPARATION */

/**
 * @brief Main separation function.
 */
static
SCIP_RETCODE separate(
   SCIP*                 scip,               /**< SCIP data structure. */
   SCIP_SEPA*            sepa,               /**< Separator. */
   SCIP_SOL*             sol,                /**< Solution to be separated (or \c NULL for the LP solution). */
   SCIP_RESULT*          result              /**< Pointer for storing the result. */
   )
{
   SCIP_SEPADATA* sepadata;

   assert(scip);

   *result = SCIP_DIDNOTFIND;

   /* We do not apply the cuts for subscips. */
   if( SCIPgetSubscipDepth(scip) > 0 )
      return SCIP_OKAY;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata);

   if( sepadata->lastrun != SCIPgetNRuns(scip) )
   {
      /* If we still have a hypergraph, we delete it. */
      if( sepadata->hypergraph )
      {
         SCIPdebugMsg(scip, "recreating hypergraph.\n");
         SCIP_CALL( SCIPhypergraphFree(&sepadata->hypergraph) );
      }
      else
         SCIPdebugMsg(scip, "creating hypergraph.\n");

      SCIP_CALL( constructHypergraph(scip, sepadata) );
      assert(sepadata->hypergraph != NULL);

      /* If the hypergraph is not interesting then we just delete it again. */
      if( SCIPhypergraphGetNOverlaps(sepadata->hypergraph) < sepadata->minnoverlaps )
      {
         SCIPdebugMsg(scip, "deleting hypergraph due to only %d overlaps.\n",
            SCIPhypergraphGetNOverlaps(sepadata->hypergraph));
         SCIP_CALL( SCIPhypergraphFree(&sepadata->hypergraph) );
         assert(sepadata->hypergraph == NULL);
      }

      sepadata->lastrun = SCIPgetNRuns(scip);
   }

   if( sepadata->hypergraph
      && ((sepadata->nuselessoneflower <= sepadata->maxuselessoneflower)
      || (sepadata->nuselesstwoflower <= sepadata->maxuselesstwoflower)) )
   {
      *result = SCIP_DIDNOTFIND;

      SCIP_CALL( prepareSeparation(scip, sepa, sol) );

      if( sepadata->maxstandard && !sepadata->delaystandard )
      {
         SCIP_CALL( separateStandard(scip, sepa, sol, sepadata->maxstandard, result) );
         if( *result == SCIP_CUTOFF )
            return SCIP_OKAY;
      }

      if( sepadata->maxoneflower && (sepadata->nuselessoneflower <= sepadata->maxuselessoneflower) )
      {
         int oldnsepaoneflower = sepadata->nsepaoneflower;
#ifdef USE_OLD_ONEFLOWER_SEPARATION
         SCIP_CALL( separateOneFlowerOld(scip, sepa, sol, sepadata->maxoneflower, result) );
#else /* !USE_OLD_ONEFLOWER_SEPARATION */
         SCIP_CALL( separateOneFlower(scip, sepa, sol, sepadata->maxoneflower, result) );
#endif /* USE_OLD_ONEFLOWER_SEPARATION */
         if( sepadata->nsepaoneflower > oldnsepaoneflower )
            sepadata->nuselessoneflower = 0;
         else
            sepadata->nuselessoneflower++;
         if( *result == SCIP_CUTOFF )
            return SCIP_OKAY;
      }

      if( sepadata->maxtwoflower && (sepadata->nuselesstwoflower <= sepadata->maxuselesstwoflower) )
      {
         int oldnsepatwoflower = sepadata->nsepatwoflower;
#ifdef USE_OLD_TWOFLOWER_SEPARATION
         SCIP_CALL( separateTwoFlowerOld(scip, sepa, sol, sepadata->maxtwoflower, result) );
#else /* !USE_OLD_ONEFLOWER_SEPARATION */
         SCIP_CALL( separateTwoFlower(scip, sepa, sol, sepadata->maxtwoflower, result) );
#endif /* USE_OLD_ONEFLOWER_SEPARATION */
         if( sepadata->nsepatwoflower > oldnsepatwoflower )
            sepadata->nuselesstwoflower = 0;
         else
            sepadata->nuselesstwoflower++;
         if( *result == SCIP_CUTOFF )
            return SCIP_OKAY;
      }

      if( *result == SCIP_SEPARATED && sepadata->maxstandard && sepadata->delaystandard )
      {
         SCIP_CALL( separateStandard(scip, sepa, sol, sepadata->maxstandard, result) );
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of separator.
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyMultilinear)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of separator */
   SCIP_CALL( SCIPincludeSepaMultilinear(scip) );

   return SCIP_OKAY;
}

/** Destructor of separator to free user data (called when SCIP is exiting). */
static
SCIP_DECL_SEPAFREE(sepaFreeMultilinear)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   /* Free the separator data. */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata);
   assert(sepadata->hypergraph == NULL);

   SCIPfreeBlockMemory(scip, &sepadata);
   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}

/** Initialization method of separator (called after problem was transformed). */
static
SCIP_DECL_SEPAINIT(sepaInitMultilinear)
{  /*lint --e{715}*/

   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata);

   sepadata->timehypercreation = 0.0;
   sepadata->timehyperoverlaps = 0.0;

   return SCIP_OKAY;
}

/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolMultilinear)
{
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata);

   if( SCIPgetSubscipDepth(scip) > 0 )
      return SCIP_OKAY;

   if( sepadata->hypergraph != NULL )
   {
      SCIP_CALL( SCIPhypergraphFree(&sepadata->hypergraph) );
   }

   return SCIP_OKAY;
}

/** LP solution separation method of separator. */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpMultilinear)
{  /*lint --e{715}*/

   SCIP_CALL( separate(scip, sepa, NULL, result) );

   return SCIP_OKAY;
}

/** arbitrary primal solution separation method of separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolMultilinear)
{  /*lint --e{715}*/

   SCIP_CALL( separate(scip, sepa, sol, result) );

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the multilinear separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaMultilinear(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata = NULL;
   SCIP_SEPA* sepa = NULL;

   /* create multilinear separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   sepadata->lastrun = -1;
   sepadata->hypergraph = NULL;
   sepadata->nsepacuts = 0;
   sepadata->timesepaoneflower = 0.0;
   sepadata->timesepatwoflower = 0.0;
   sepadata->timepreparation = 0.0;
   sepadata->nsepastandard = 0;
   sepadata->nsepaoneflower = 0;
   sepadata->nsepatwoflower = 0;
   sepadata->nuselessoneflower = 0;
   sepadata->nuselesstwoflower = 0;

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY, sepaExeclpMultilinear, sepaExecsolMultilinear, sepadata) );

   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyMultilinear) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeMultilinear) );
   SCIP_CALL( SCIPsetSepaInit(scip, sepa, sepaInitMultilinear) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitsolMultilinear) );

   /* add multilinear separator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/multilinear/scanand",
         "Whether to scan AND constraints when constructing hypergraph", &sepadata->scanand, FALSE, DEFAULT_SCAN_AND, 0,
         0) );
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/multilinear/scanproduct",
         "Whether to scan product expressions when constructing hypergraph", &sepadata->scanproduct, FALSE,
         DEFAULT_SCAN_PRODUCT, 0, 0) );
   SCIP_CALL( SCIPaddIntParam(scip, "separating/multilinear/maxstandard",
         "Maximum number of standard relaxation inequalities per cut round", &sepadata->maxstandard, FALSE,
         DEFAULT_MAX_STANDARD, 0, INT_MAX, 0, 0) );
   SCIP_CALL( SCIPaddIntParam(scip, "separating/multilinear/maxoneflower",
         "Maximum number of 1-flower inequalities per cut round", &sepadata->maxoneflower, FALSE, DEFAULT_MAX_ONEFLOWER,
         0, INT_MAX, 0, 0) );
   SCIP_CALL( SCIPaddIntParam(scip, "separating/multilinear/maxtwoflower",
         "Maximum number of 2-flower inequalities per cut round", &sepadata->maxtwoflower, FALSE, DEFAULT_MAX_TWOFLOWER,
         0, INT_MAX, 0, 0) );
   SCIP_CALL( SCIPaddIntParam(scip, "separating/multilinear/minnoverlaps",
         "Minimum number of overlaps necessary to try separation", &sepadata->minnoverlaps, FALSE, DEFAULT_MIN_OVERLAPS,
         0, INT_MAX, 0, 0) );
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/multilinear/delaystandard",
         "Whether to only generate standard inequalities if also flowers were generated", &sepadata->delaystandard,
         FALSE, DEFAULT_DELAY_STANDARD, 0, 0) );
   SCIP_CALL( SCIPaddIntParam(scip, "separating/multilinear/maxuselessoneflower",
         "Number of useless separation rounds after which we stop separating 1-flowers", &sepadata->maxuselessoneflower,
         FALSE, DEFAULT_MAX_USELESS_ONEFLOWER, 0, INT_MAX, 0, 0) );
   SCIP_CALL( SCIPaddIntParam(scip, "separating/multilinear/maxuselesstwoflower",
         "Number of useless separation rounds after which we stop separating 2-flowers", &sepadata->maxuselesstwoflower,
         FALSE, DEFAULT_MAX_USELESS_TWOFLOWER, 0, INT_MAX, 0, 0) );

   return SCIP_OKAY;
}

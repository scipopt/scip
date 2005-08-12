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
#pragma ident "@(#) $Id: sepa_clique.c,v 1.13 2005/08/12 12:36:22 bzfpfend Exp $"

/**@file   sepa_clique.c
 * @brief  clique separator
 * @author Kati Wolter
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/sepa_clique.h"
#include "scip/tclique.h"


#define SEPA_NAME              "clique"
#define SEPA_DESC              "clique separator of stable set relaxation"
#define SEPA_PRIORITY             -5000
#define SEPA_FREQ                    10
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_SCALEVAL         1000.0 /**< factor for scaling weights */
#define DEFAULT_MAXTREENODES         -1 /**< maximal number of nodes in branch and bound tree (-1: no limit) */




/*
 * Data structures
 */

/** separator data */
struct SepaData
{
   TCLIQUEGRAPH*    tcliquegraph;       /**< tclique graph data structure */
   SCIP*            scip;               /**< SCIP data structure */
   Real*            varsolvals;         /**< LP solution of binary variables (contained in a 3-clique in implgraph) */
   Real             scaleval;           /**< factor for scaling weights */
   int              maxtreenodes;       /**< maximal number of nodes in branch and bound tree (-1: no limit) */
   int              ncuts;              /**< number of cuts found */
   Bool             tcliquegraphloaded; /**< TRUE if tcliquegraph is allready loaded (tcliquegraph can be NULL),
                                         *   FALSE otherwise */ 
};

/** tclique graph data */
struct TcliqueGraph
{
   VAR**            vars;               /**< active problem variables (or negated variables) the nodes belong to */
   WEIGHT*          weights;	        /**< weight of nodes */
   int*             adjnodesidxs;       /**< indices in adjnodes array of first adjacent nodes for each node */
   int*             cliqueidsidxs;      /**< indices in cliqueids array of first clique the node is contained in */
   int*             adjnodes;	        /**< adjacent nodes of edges */
   int*             cliqueids;          /**< unique ids of cliques */
   int              adjnodessize;       /**< size of adjnodes array */
   int              cliqueidssize;      /**< size of cliqueids array */
   int              nnodes;		/**< number of nodes in graph */
};




/*
 * Methods for tclique graph
 */

/** creates an empty tclique graph data structure */
static
RETCODE tcliquegraphCreate(
   SCIP*            scip,               /**< SCIP data structure */
   TCLIQUEGRAPH**   tcliquegraph        /**< pointer to tclique graph data */
   )
{
   int maxnnodes;

   assert(tcliquegraph != NULL);

   CHECK_OKAY( SCIPallocMemory(scip, tcliquegraph) );
   
   /* there are at most 2*nbinvars nodes in the graph */
   maxnnodes = 2*SCIPgetNBinVars(scip);
   assert(maxnnodes > 0);

   /* allocate memory for tclique graph arrays */
   CHECK_OKAY( SCIPallocMemoryArray(scip, &(*tcliquegraph)->vars, maxnnodes) );
   CHECK_OKAY( SCIPallocMemoryArray(scip, &(*tcliquegraph)->weights, maxnnodes) );
   CHECK_OKAY( SCIPallocMemoryArray(scip, &(*tcliquegraph)->adjnodesidxs, maxnnodes+1) );
   CHECK_OKAY( SCIPallocMemoryArray(scip, &(*tcliquegraph)->cliqueidsidxs, maxnnodes+1) );
   (*tcliquegraph)->adjnodesidxs[0] = 0;  /* the last slot defines the end of the last node */
   (*tcliquegraph)->cliqueidsidxs[0] = 0; /* the last slot defines the end of the last node */
   (*tcliquegraph)->adjnodes = NULL;
   (*tcliquegraph)->cliqueids = NULL;
   (*tcliquegraph)->adjnodessize = 0;
   (*tcliquegraph)->cliqueidssize = 0;
   (*tcliquegraph)->nnodes = 0;

   return SCIP_OKAY;
}

/** frees the tclique graph data structure and releases all captured variables */
static
RETCODE tcliquegraphFree(
   SCIP*            scip,               /**< SCIP data structure */
   TCLIQUEGRAPH**   tcliquegraph        /**< pointer to tclique graph data */
   )
{
   int v;

   assert(tcliquegraph != NULL);
   assert(*tcliquegraph != NULL);
   
   /* release variables */
   for( v = 0; v < (*tcliquegraph)->nnodes; ++v )
   {
      CHECK_OKAY( SCIPreleaseVar(scip, &(*tcliquegraph)->vars[v]) );
   }

   /* free memory */
   SCIPfreeMemoryArray(scip, &(*tcliquegraph)->vars);
   SCIPfreeMemoryArray(scip, &(*tcliquegraph)->weights);
   SCIPfreeMemoryArray(scip, &(*tcliquegraph)->adjnodesidxs);
   SCIPfreeMemoryArray(scip, &(*tcliquegraph)->cliqueidsidxs);
   SCIPfreeMemoryArrayNull(scip, &(*tcliquegraph)->adjnodes);
   SCIPfreeMemoryArrayNull(scip, &(*tcliquegraph)->cliqueids);
   SCIPfreeMemory(scip, tcliquegraph);

   return SCIP_OKAY;
}

/** ensures that the adjnodes array can store at least num entries */
static
RETCODE tcliquegraphEnsureAdjnodesSize(
   SCIP*            scip,               /**< SCIP data structure */
   TCLIQUEGRAPH*    tcliquegraph,       /**< tclique graph data */
   int              num                 /**< minimal number of adjacent nodes to be able to store in the array */
   )
{
   assert(tcliquegraph != NULL);

   if( num > tcliquegraph->adjnodessize )
   {
      tcliquegraph->adjnodessize = SCIPcalcMemGrowSize(scip, num);
      CHECK_OKAY( SCIPreallocMemoryArray(scip, &tcliquegraph->adjnodes, tcliquegraph->adjnodessize) );
   }
   assert(num <= tcliquegraph->adjnodessize);

   return SCIP_OKAY;
}

/** ensures that the cliqueids array can store at least num entries */
static
RETCODE tcliquegraphEnsureCliqueidsSize(
   SCIP*            scip,               /**< SCIP data structure */
   TCLIQUEGRAPH*    tcliquegraph,       /**< tclique graph data */
   int              num                 /**< minimal number of adjacent nodes to be able to store in the array */
   )
{
   assert(tcliquegraph != NULL);

   if( num > tcliquegraph->cliqueidssize )
   {
      tcliquegraph->cliqueidssize = SCIPcalcMemGrowSize(scip, num);
      CHECK_OKAY( SCIPreallocMemoryArray(scip, &tcliquegraph->cliqueids, tcliquegraph->cliqueidssize) );
   }
   assert(num <= tcliquegraph->cliqueidssize);

   return SCIP_OKAY;
}

/** adds a node to the tclique graph defined as a variable-value pair; adds all cliques to the cliqueids array the
 *  variable is contained in with the given value
 */
static
RETCODE tcliquegraphAddNode(
   SCIP*            scip,               /**< SCIP data structure */
   TCLIQUEGRAPH**   tcliquegraph,       /**< pointer to tclique graph data */
   VAR*             var,                /**< active binary problem variable */
   Bool             value,              /**< value of the variable in the node */
   int*             nodeidx             /**< pointer to store the index of the new node */
   )
{
   VAR* nodevar;
   int* cliqueids;
   CLIQUE** cliques;
   int ncliques;
   int nadjnodes;
   int ncliqueids;
   int i;

   assert(tcliquegraph != NULL);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
   assert(SCIPvarIsActive(var));
   assert(nodeidx != NULL);

   /* create tclique graph data if not yet existing */
   if( *tcliquegraph == NULL )
   {
      CHECK_OKAY( tcliquegraphCreate(scip, tcliquegraph) );
   }
   assert(*tcliquegraph != NULL);
   assert((*tcliquegraph)->nnodes < 2*SCIPgetNBinVars(scip) - 1);

   /* if the value is FALSE, use the negated variable for the node */
   if( !value )
   {
      CHECK_OKAY( SCIPgetNegatedVar(scip, var, &nodevar) );
   }
   else
      nodevar = var;

   /* get the current number of used entries in adjnodes and cliqueids arrays */
   nadjnodes = (*tcliquegraph)->adjnodesidxs[(*tcliquegraph)->nnodes];
   ncliqueids = (*tcliquegraph)->cliqueidsidxs[(*tcliquegraph)->nnodes];

   /* insert the variable into the tclique graph */
   *nodeidx = (*tcliquegraph)->nnodes;
   CHECK_OKAY( SCIPcaptureVar(scip, nodevar) );
   (*tcliquegraph)->vars[*nodeidx] = nodevar;
   (*tcliquegraph)->weights[*nodeidx] = 0;
   (*tcliquegraph)->nnodes++;

   /* store the ids of the variable's cliques in the cliqueids array */
   ncliques = SCIPvarGetNCliques(var, value);
   cliques = SCIPvarGetCliques(var, value);
   CHECK_OKAY( tcliquegraphEnsureCliqueidsSize(scip, *tcliquegraph, ncliqueids + ncliques) );
   cliqueids = (*tcliquegraph)->cliqueids;
   for( i = 0; i < ncliques; ++i )
   {
      assert(ncliqueids < (*tcliquegraph)->cliqueidssize);
      cliqueids[ncliqueids] = SCIPcliqueGetId(cliques[i]);
      assert(i == 0 || cliqueids[ncliqueids-1] <= cliqueids[ncliqueids]);
      ncliqueids++;
   }

   /* store the new number of used entries in adjnodes and cliqueids arrays */
   (*tcliquegraph)->adjnodesidxs[(*tcliquegraph)->nnodes] = nadjnodes;
   (*tcliquegraph)->cliqueidsidxs[(*tcliquegraph)->nnodes] = ncliqueids;

   return SCIP_OKAY;
}

/** adds all variable/value pairs to the tclique graph that are contained in an existing 3-clique */
static
RETCODE tcliquegraphAddCliqueVars(
   SCIP*            scip,               /**< SCIP data structure */
   TCLIQUEGRAPH**   tcliquegraph,       /**< pointer to tclique graph data */
   int**            cliquegraphidx      /**< array to store tclique graph node index of variable/value pairs */
   )
{
   VAR** vars;
   int nvars;
   int i;

   assert(tcliquegraph != NULL);
   assert(cliquegraphidx != NULL);
   assert(cliquegraphidx[0] != NULL);
   assert(cliquegraphidx[1] != NULL);

   /* get binary problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNBinVars(scip);

   for( i = 0; i < nvars; ++i )
   {
      VAR* var;
      int value;

      var = vars[i];
      for( value = 0; value < 2; ++value )
      {
         assert(cliquegraphidx[value][i] == -1);

         if( SCIPvarGetNCliques(var, (Bool)value) >= 1 )
         {
            /* all cliques stored in the clique table are at least 3-cliques */
            assert(SCIPcliqueGetNVars(SCIPvarGetCliques(var, (Bool)value)[0]) >= 3);
            CHECK_OKAY( tcliquegraphAddNode(scip, tcliquegraph, var, (Bool)value, &cliquegraphidx[value][i]) );
         }
      }
   }

   return SCIP_OKAY;
}

/* adds all variable/value pairs to the tclique graph that are contained in a 3-clique in the implication graph */
static
RETCODE tcliquegraphAddImplicsVars(
   SCIP*            scip,               /**< SCIP data structure */
   TCLIQUEGRAPH**   tcliquegraph,       /**< pointer to tclique graph data */
   int**            cliquegraphidx      /**< array to store tclique graph node index of variable/value pairs */
   )
{
   VAR** vars;
   int nvars;
   int xi;

   assert(tcliquegraph != NULL);
   assert(cliquegraphidx != NULL);
   assert(cliquegraphidx[0] != NULL);
   assert(cliquegraphidx[1] != NULL);

   /* get binary variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNBinVars(scip);
   assert(nvars > 0);

   /* detect 3-cliques in the clique graph: triples (x,y,z) in the implication graph with x -> y', x -> z', and y -> z';
    * in order to avoid double checks, we only check triplets with variable indices xindex < yindex < zindex;
    * we don't have to check triples with x == y, x == z, or y == z, because cliques with the same variable occuring
    * twice lead to fixings of this or all other variables which should already have been detected in presolving
    */
   for( xi = 0; xi < nvars-2; ++xi ) /* at least two variables must be left over for y and z */
   {
      VAR* x;
      int xindex;
      int xvalue;

      x = vars[xi];
      xindex = SCIPvarGetIndex(x);

      for( xvalue = 0; xvalue < 2; xvalue++ )
      {
         VAR** ximplvars;
         BOUNDTYPE* ximpltypes; 
         int xnbinimpls;
         int i;

         /* scan the implications of x == xvalue for potential y candidates */
         xnbinimpls = SCIPvarGetNBinImpls(x, (Bool)xvalue);
         ximplvars = SCIPvarGetImplVars(x, (Bool)xvalue);
         ximpltypes = SCIPvarGetImplTypes(x, (Bool)xvalue);
      
         /* ignore implicants with yindex <= xindex */
         for( i = 0; i < xnbinimpls-1 && SCIPvarGetIndex(ximplvars[i]) <= xindex; ++i )
         {}

         /* loop over all y > x */
         for( ; i < xnbinimpls-1; ++i ) /* at least one variable must be left over for z */
         {
            VAR* y;
            int yindex;
            int yvalue;
            int yi;
            VAR** yimplvars;
            BOUNDTYPE* yimpltypes; 
            int ynbinimpls;
            int xk;
            int yk;

            y = ximplvars[i];
            yi = SCIPvarGetProbindex(y);
            yindex = SCIPvarGetIndex(y);
            assert(0 <= yi && yi < nvars);
            assert(xindex < yindex); /* the implied variables are sorted by increasing variable index */
            assert(yindex < SCIPvarGetIndex(ximplvars[i+1]));

            /* check, whether the implicant is y == 0 or y == 1 (yi conflicts with x == xvalue) */
            yvalue = (int)(ximpltypes[i] == SCIP_BOUNDTYPE_UPPER);
            assert(0 <= yvalue && yvalue <= 1);

            /* scan the remaining implications of x == xvalue and the implications of y == yvalue for equal entries */
            ynbinimpls = SCIPvarGetNBinImpls(y, (Bool)yvalue);
            yimplvars = SCIPvarGetImplVars(y, (Bool)yvalue);
            yimpltypes = SCIPvarGetImplTypes(y, (Bool)yvalue);

            /* we simultaneously scan both implication arrays for candidates z with x < y < z */
            xk = i+1;
            for( yk = 0; yk < ynbinimpls && SCIPvarGetIndex(yimplvars[yk]) <= yindex; ++yk )
            {}
            while( xk < xnbinimpls && yk < ynbinimpls )
            {
               int zindex;
               
               assert(yindex < SCIPvarGetIndex(ximplvars[xk]));
               assert(yindex < SCIPvarGetIndex(yimplvars[yk]));

               /* scan the implications of x */
               zindex = SCIPvarGetIndex(yimplvars[yk]);
               while( xk < xnbinimpls && SCIPvarGetIndex(ximplvars[xk]) < zindex )
                  xk++;
               if( xk >= xnbinimpls )
                  break;

               /* scan the implications of y */
               zindex = SCIPvarGetIndex(ximplvars[xk]);
               while( yk < ynbinimpls && SCIPvarGetIndex(yimplvars[yk]) < zindex )
                  yk++;
               if( yk >= ynbinimpls )
                  break;

               /* check, whether we reached a common implied variable */
               if( ximplvars[xk] != yimplvars[yk] )
                  continue;

               assert(SCIPvarGetIndex(ximplvars[xk]) == zindex);
               assert(SCIPvarGetIndex(yimplvars[yk]) == zindex);

               /* check, whether both implications are of the same type */
               if( ximpltypes[xk] == yimpltypes[yk] )
               {
                  VAR* z;
                  int zi;
                  int zvalue;

                  /* we found z with xindex < yindex < zindex and x + y + z <= 1 */
                  z = ximplvars[xk];
                  zi = SCIPvarGetProbindex(z);
                  assert(0 <= zi && zi < nvars);
                  assert(SCIPvarGetIndex(z) == zindex);
                  assert(xindex < yindex && yindex < zindex);

                  /* check, whether the implicant is z == 0 or z == 1 (z conflicts with x == xvalue) */
                  zvalue = (int)(ximpltypes[xk] == SCIP_BOUNDTYPE_UPPER);
                  assert(0 <= zvalue && zvalue <= 1);

                  /* add nodes x == xvalue, y == yvalue, and z == zvalue to clique graph (if not yet existing) */
                  if( cliquegraphidx[xvalue][xi] == -1 )
                  {
                     CHECK_OKAY( tcliquegraphAddNode(scip, tcliquegraph, x, (Bool)xvalue, &cliquegraphidx[xvalue][xi]) );
                  }
                  if( cliquegraphidx[yvalue][yi] == -1 )
                  {
                     CHECK_OKAY( tcliquegraphAddNode(scip, tcliquegraph, y, (Bool)yvalue, &cliquegraphidx[yvalue][yi]) );
                  }
                  if( cliquegraphidx[zvalue][zi] == -1 )
                  {
                     CHECK_OKAY( tcliquegraphAddNode(scip, tcliquegraph, z, (Bool)zvalue, &cliquegraphidx[zvalue][zi]) );
                  }
               }
               
               /* proceed with the next pair of implications */
               xk++;
               yk++;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/* adds all variable/value pairs to the tclique graph that have implications to two variables of the same existing
 * clique
 */
static
RETCODE tcliquegraphAddImplicsCliqueVars(
   SCIP*            scip,               /**< SCIP data structure */
   TCLIQUEGRAPH**   tcliquegraph,       /**< pointer to tclique graph data */
   int**            cliquegraphidx      /**< array to store tclique graph node index of variable/value pairs */
   )
{
   VAR** vars;
   int nvars;
   int xi;

   assert(tcliquegraph != NULL);
   assert(cliquegraphidx != NULL);
   assert(cliquegraphidx[0] != NULL);
   assert(cliquegraphidx[1] != NULL);

   /* get binary variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNBinVars(scip);
   assert(nvars > 0);

   /* detect triples (x,y,z) with implications x -> y' and x -> z' and a clique that contains y and z;
    * because all cliques stored in the clique table are at least 3-cliques, y and z are already member of the
    * tclique graph due to tcliquegraphAddCliqueVars(); therefore, we only have to process pairs x == xvalue
    * that are currently not member of the tclique graph
    */
   for( xi = 0; xi < nvars; ++xi )
   {
      VAR* x;
      int xvalue;

      x = vars[xi];
      for( xvalue = 0; xvalue < 2; xvalue++ )
      {
         VAR** ximplvars;
         BOUNDTYPE* ximpltypes; 
         int xnbinimpls;
         int yk;

         /* we only need to process pairs x == xvalue that are not yet member of the tclique graph */
         if( cliquegraphidx[xvalue][xi] >= 0 )
            continue;

         /* scan the implications of x == xvalue for potential y and z candidates */
         xnbinimpls = SCIPvarGetNBinImpls(x, (Bool)xvalue);
         ximplvars = SCIPvarGetImplVars(x, (Bool)xvalue);
         ximpltypes = SCIPvarGetImplTypes(x, (Bool)xvalue);
      
         /* loop over all y */
         for( yk = 0; yk < xnbinimpls-1; ++yk ) /* at least one variable must be left over for z */
         {
            VAR* y;
            int yvalue;
            int yi;
            int zk;
            
            y = ximplvars[yk];
            yi = SCIPvarGetProbindex(y);

            /* check, whether the implicant is y == 0 or y == 1 (y conflicts with x == xvalue) */
            yvalue = (int)(ximpltypes[yk] == SCIP_BOUNDTYPE_UPPER);
            if( SCIPvarGetNCliques(y, yvalue) == 0 )
               continue;

            /* loop over all z > y */
            for( zk = yk+1; zk < xnbinimpls; ++zk )
            {
               VAR* z;
               int zvalue;
               int zi;
               
               z = ximplvars[zk];
               zi = SCIPvarGetProbindex(z);
               
               /* check, whether the implicant is z == 0 or z == 1 (z conflicts with x == xvalue) */
               zvalue = (int)(ximpltypes[zk] == SCIP_BOUNDTYPE_UPPER);

               /* check whether y and z have a common clique */
               if( SCIPvarsHaveCommonClique(y, yvalue, z, zvalue, FALSE) )
               {
                  /* add nodes x == xvalue to clique graph */
                  assert(cliquegraphidx[xvalue][xi] == -1);
                  assert(cliquegraphidx[yvalue][yi] >= 0);
                  assert(cliquegraphidx[zvalue][zi] >= 0);
                  CHECK_OKAY( tcliquegraphAddNode(scip, tcliquegraph, x, (Bool)xvalue, &cliquegraphidx[xvalue][xi]) );
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** adds all implications between nodes in the tclique graph to the tclique graph */
static
RETCODE tcliquegraphAddImplics(
   SCIP*            scip,               /**< SCIP data structure */
   TCLIQUEGRAPH*    tcliquegraph,       /**< tclique graph data */
   int**            cliquegraphidx      /**< array to store tclique graph node index of variable/value pairs */
   )
{
   int nadjnodes;
   int i;

   assert(cliquegraphidx != NULL);
   assert(cliquegraphidx[0] != NULL);
   assert(cliquegraphidx[1] != NULL);

   /* there is nothing to do, if the graph is empty */
   if( tcliquegraph == NULL )
      return SCIP_OKAY;
   
   /* the tclique graph should currently contain no implications */
   assert(tcliquegraph->adjnodes == NULL);
   assert(tcliquegraph->adjnodessize == 0);

   nadjnodes = 0;
   for( i = 0; i < tcliquegraph->nnodes; ++i )
   {
      VAR* var;
      Bool value;
      VAR** implvars;
      BOUNDTYPE* impltypes;
      int nbinimpls;
      int adjnodesidx;
      int j;
      
      /* store start adjnodes index of the node */
      adjnodesidx = nadjnodes;
      tcliquegraph->adjnodesidxs[i] = adjnodesidx;

      /* get active problem variable */
      var = tcliquegraph->vars[i];
      value = TRUE;
      CHECK_OKAY( SCIPvarGetProbvarBinary(&var, &value) );
      assert(0 <= SCIPvarGetProbindex(var) && SCIPvarGetProbindex(var) < SCIPgetNBinVars(scip));
      assert(cliquegraphidx[value][SCIPvarGetProbindex(var)] == i);

      /* get implications on binary variables */
      nbinimpls = SCIPvarGetNBinImpls(var, value);
      implvars = SCIPvarGetImplVars(var, value);
      impltypes = SCIPvarGetImplTypes(var, value);
      for( j = 0; j < nbinimpls; ++j )
      {
         int probidx;
         int graphidx;
         Bool implvalue;

         probidx = SCIPvarGetProbindex(implvars[j]);
         implvalue = (impltypes[j] == SCIP_BOUNDTYPE_UPPER);
         assert(0 <= probidx && probidx < SCIPgetNBinVars(scip));

         graphidx = cliquegraphidx[implvalue][probidx];
         if( graphidx >= 0 )
         {
            int pos;

            assert(graphidx < tcliquegraph->nnodes);
            assert((implvalue == TRUE && tcliquegraph->vars[graphidx] == implvars[j])
               || (implvalue == FALSE && SCIPvarGetNegationVar(tcliquegraph->vars[graphidx]) == implvars[j]));

            /* allocate memory for addional arc */
            CHECK_OKAY( tcliquegraphEnsureAdjnodesSize(scip, tcliquegraph, nadjnodes+1) );

            /* store the adjacent node in the tclique graph data structure, sorted by index */
            for( pos = nadjnodes; pos > adjnodesidx && tcliquegraph->adjnodes[pos-1] > graphidx; --pos )
               tcliquegraph->adjnodes[pos] = tcliquegraph->adjnodes[pos-1];
            tcliquegraph->adjnodes[pos] = graphidx;
            nadjnodes++;
         }
      }
   }
   assert(nadjnodes/2 == (nadjnodes+1)/2); /* we must have an even number of arcs */

   /* store final adjnodes index */
   tcliquegraph->adjnodesidxs[tcliquegraph->nnodes] = nadjnodes;

   return SCIP_OKAY;
}

/* creates tclique data structure using the implication graph;
 * only variables that are contained in a 3-clique are added as nodes to the clique graph
 */
static 
RETCODE loadTcliquegraph(
   SCIP*            scip,               /**< SCIP data structure */
   SEPADATA*        sepadata            /**< separator data */
   )
{
   int* cliquegraphidx[2];
   int nvars;
   int i;

   assert(sepadata != NULL);
   assert(sepadata->tcliquegraph == NULL);

   /* there is nothing to do, if no binary variables are present in the problem */
   nvars = SCIPgetNBinVars(scip);
   if( nvars == 0 )
      return SCIP_OKAY;

   /* get temporary memory for mapping variable/value pairs to clique graph nodes */
   CHECK_OKAY( SCIPallocBufferArray(scip, &cliquegraphidx[0], nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &cliquegraphidx[1], nvars) );
   for( i = 0; i < nvars; ++i )
   {
      cliquegraphidx[0][i] = -1;
      cliquegraphidx[1][i] = -1;
   }

   /* insert all variable/value pairs that are contained in an existing 3-clique */
   CHECK_OKAY( tcliquegraphAddCliqueVars(scip, &sepadata->tcliquegraph, cliquegraphidx) );

   /* insert all variable/value pairs that are contained in a 3-clique in the implication graph */
   CHECK_OKAY( tcliquegraphAddImplicsVars(scip, &sepadata->tcliquegraph, cliquegraphidx) );

   /* insert all variable/value pairs that have implications to two variables of the same existing clique */
   CHECK_OKAY( tcliquegraphAddImplicsCliqueVars(scip, &sepadata->tcliquegraph, cliquegraphidx) );

   /* add all implications between used variables to the tclique graph */
   CHECK_OKAY( tcliquegraphAddImplics(scip, sepadata->tcliquegraph, cliquegraphidx) );

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &cliquegraphidx[1]);
   SCIPfreeBufferArray(scip, &cliquegraphidx[0]);

   return SCIP_OKAY;
}

/* updates the weights in the tclique graph data structure */
static 
void updateTcliquegraph(
   SCIP*            scip,               /**< SCIP data structure */
   SEPADATA*        sepadata            /**< separator data */
   )
{
   TCLIQUEGRAPH* tcliquegraph;
   int i;

   assert(sepadata != NULL);
   assert(sepadata->varsolvals != NULL);

   tcliquegraph = sepadata->tcliquegraph;
   assert(tcliquegraph != NULL);

   /* updates weight of all nodes in tclique data structure */
   for( i = 0; i < tcliquegraph->nnodes; i++ )
   {
      int weight;

      weight = (WEIGHT)SCIPfeasFloor(scip, sepadata->varsolvals[i] * sepadata->scaleval);
      tcliquegraph->weights[i] = MAX(weight, 0);
   }
}




/*
 * TClique Graph Callbacks
 */

/** gets number of nodes in the graph */
static
TCLIQUE_GETNNODES(tcliqueGetnnodesClique)
{
   assert(tcliquegraph != NULL);

   return tcliquegraph->nnodes;
}

/** gets weight of nodes in the graph */
static
TCLIQUE_GETWEIGHTS(tcliqueGetweightsClique)
{
   assert(tcliquegraph != NULL);

   return tcliquegraph->weights;
}

/** returns whether the nodes are member of a common clique */
static
Bool nodesHaveCommonClique(
   TCLIQUEGRAPH*    tcliquegraph,       /**< tclique graph data */
   int              node1,              /**< first node */
   int              node2               /**< second node */
   )
{
   int* cliqueids;
   int i1;
   int i2;
   int endi1;
   int endi2;

   assert(tcliquegraph != NULL);

   cliqueids = tcliquegraph->cliqueids;
   i1 = tcliquegraph->cliqueidsidxs[node1];
   endi1 = tcliquegraph->cliqueidsidxs[node1+1];
   i2 = tcliquegraph->cliqueidsidxs[node2];
   endi2 = tcliquegraph->cliqueidsidxs[node2+1];
   while( i1 < endi1 && i2 < endi2 )
   {
      while( i1 < endi1 && cliqueids[i1] < cliqueids[i2] )
         i1++;
      if( i1 == endi1 )
         break;

      while( i2 < endi2 && cliqueids[i2] < cliqueids[i1] )
         i2++;
      if( i2 == endi2 )
         break;

      if( cliqueids[i1] == cliqueids[i2] )
         return TRUE;
   }

   return FALSE;
}

/** returns, whether the edge (node1, node2) is in the graph */
static
TCLIQUE_ISEDGE(tcliqueIsedgeClique)
{
   int left;
   int right;

   assert(tcliquegraph != NULL);
   assert(0 <= node1 && node1 < tcliquegraph->nnodes);
   assert(0 <= node2 && node2 < tcliquegraph->nnodes);

   /* check if node2 is contained in adjacency list of node1 (list is ordered by adjacent nodes) */
   left = tcliquegraph->adjnodesidxs[node1];
   right = tcliquegraph->adjnodesidxs[node1+1]-1;
   while( left <= right )
   {
      int middle;
      int node;

      middle = (left+right)/2;
      node = tcliquegraph->adjnodes[middle];
      if( node < node2 )
         left = middle+1;
      else if( node > node2 )
         right = middle-1;
      else
         return TRUE;
   }

   /* check if the nodes are member of a common clique */
   return nodesHaveCommonClique(tcliquegraph, node1, node2);
}

/** selects all nodes from a given set of nodes which are adjacent to a given node
 *  and returns the number of selected nodes
 */
static
TCLIQUE_SELECTADJNODES(tcliqueSelectadjnodesClique)
{
   int* graphadjnodes;
   int nadjnodes;
   int nodeadjindex;
   int nodeadjend;
   int i;

   assert(tcliquegraph != NULL);
   assert(0 <= node && node < tcliquegraph->nnodes);
   assert(nnodes == 0 || nodes != NULL);
   assert(adjnodes != NULL);

   nadjnodes = 0;

   /* check for each node in given nodes set, if it is adjacent to the given node or shares a common clique */
   graphadjnodes = tcliquegraph->adjnodes;
   nodeadjindex = tcliquegraph->adjnodesidxs[node];
   nodeadjend = tcliquegraph->adjnodesidxs[node+1];
   for( i = 0; i < nnodes; i++ )
   {
      /* check if the node is adjacent to the given node (nodes and adjacent nodes are ordered by node index) */
      assert(0 <= nodes[i] && nodes[i] < tcliquegraph->nnodes);
      assert(i == 0 || nodes[i-1] < nodes[i]);
      while( nodeadjindex < nodeadjend && graphadjnodes[nodeadjindex] < nodes[i] )
         nodeadjindex++;
      if( nodeadjindex < nodeadjend && graphadjnodes[nodeadjindex] == nodes[i] )
      {
         /* current node is adjacent to given node */
         adjnodes[nadjnodes] = nodes[i]; 
         nadjnodes++;
      }
      else
      {
         /* current node is not adjacent to given node: check if they share a common clique */
         if( nodesHaveCommonClique(tcliquegraph, node, nodes[i]) )
         {
            adjnodes[nadjnodes] = nodes[i]; 
            nadjnodes++;
         }         
      }
   }
   
   return nadjnodes;
}

/** generates cuts using a clique found by algorithm for maximum weight clique 
 *  and decides whether to stop generating cliques with the algorithm for maximum weight clique
 */
static
TCLIQUE_NEWSOL(tcliqueNewsolClique)
{
   SEPADATA* sepadata;
   WEIGHT minweightinc;

   assert(acceptsol != NULL);
   assert(stopsolving != NULL);

   sepadata = (SEPADATA*)tcliquedata;
   assert(sepadata != NULL);
   assert(sepadata->scip != NULL);
   assert(sepadata->tcliquegraph != NULL);
   assert(sepadata->ncuts >= 0);
   
   /* we don't accept the solution as new incumbent, because we want to find many violated clique inequalities */
   *acceptsol = FALSE;
   *stopsolving = FALSE;

   /* slightly increase the minimal weight for additional cliques */
   minweightinc = (cliqueweight - *minweight)/MIN(1, (int)(sepadata->scaleval/10.0));
   minweightinc = MIN(minweightinc, 1);
   *minweight += minweightinc;

   /* adds cut if weight of the clique is greater than 1 */
   if( cliqueweight > sepadata->scaleval )
   {
      SCIP* scip;
      Real* varsolvals;
      Real unscaledweight;
      int i;

      scip = sepadata->scip;
      varsolvals = sepadata->varsolvals; 
      assert(varsolvals != NULL);

      /* calculate the weight of the clique in unscaled fractional variable space */
      unscaledweight = 0.0;
      for( i = 0; i < ncliquenodes; i++ )
         unscaledweight += varsolvals[cliquenodes[i]];
      
      if( SCIPisEfficacious(scip, unscaledweight - 1.0) )
      {
         VAR** vars;
         int nvars;
         ROW* cut;
         char cutname[MAXSTRLEN];
      
         nvars = sepadata->tcliquegraph->nnodes; 
         vars = sepadata->tcliquegraph->vars; 
         assert(nvars > 0);
         assert(vars != NULL);

         /* create the cut */
         sprintf(cutname, "clique%d_%d", SCIPgetNLPs(scip), sepadata->ncuts);
         CHECK_ABORT( SCIPcreateEmptyRow(scip, &cut, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );

         CHECK_ABORT( SCIPcacheRowExtensions(scip, cut) );
         assert(ncliquenodes <= nvars);
         debugMessage("clique in graph:");
         for( i = 0; i < ncliquenodes; ++i )
         {
            assert(cliquenodes[i] < nvars);
            CHECK_ABORT( SCIPaddVarToRow(scip, cut, vars[cliquenodes[i]], 1.0) );
            debugPrintf(" [%d]<%s>", cliquenodes[i], SCIPvarGetName(vars[cliquenodes[i]]));
         }
         debugPrintf("\n");
         CHECK_ABORT( SCIPflushRowExtensions(scip, cut) );

         debugMessage("found clique cut (act=%g): ", unscaledweight);
         debug(SCIPprintRow(scip, cut, NULL));

         CHECK_ABORT( SCIPaddCut(scip, cut, FALSE) );
         CHECK_ABORT( SCIPaddPoolCut(scip, cut) );
         sepadata->ncuts++;
         
         /* release the row */
         CHECK_ABORT( SCIPreleaseRow(scip, &cut) );
      }
   }
}





/*
 * Callback methods of separator
 */

/** destructor of separator to free user data (called when SCIP is exiting) */
static
DECL_SEPAFREE(sepaFreeClique)
{  /*lint --e{715}*/
   SEPADATA* sepadata;
   
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->tcliquegraph == NULL);

   SCIPfreeMemory(scip, &sepadata);
   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** initialization method of separator (called after problem was transformed) */
#define sepaInitClique NULL


/** deinitialization method of separator (called before transformed problem is freed) */
#define sepaExitClique NULL


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#define sepaInitsolClique NULL


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
DECL_SEPAEXITSOL(sepaExitsolClique)
{  /*lint --e{715}*/
   SEPADATA* sepadata;
   
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* free tclique data */
   if( sepadata->tcliquegraph != NULL )
   {
      CHECK_OKAY( tcliquegraphFree(scip, &sepadata->tcliquegraph) );
   }
   assert(sepadata->tcliquegraph == NULL);
   sepadata->tcliquegraphloaded = FALSE;

   return SCIP_OKAY;
}

/** execution method of separator */
static
DECL_SEPAEXEC(sepaExecClique)
{  /*lint --e{715}*/
   SEPADATA* sepadata;
   TCLIQUEGRAPH* tcliquegraph;
   int* cliquenodes; 	        
   WEIGHT cliqueweight;    
   int ncliquenodes;	        
   int maxtreenodes;

   assert(scip != NULL);
   
   *result = SCIP_DIDNOTRUN;

   /* get separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   
   sepadata->ncuts = 0;

   /* if we already detected that no implications between binary variables exist, nothing has to be done */
   if( sepadata->tcliquegraph == NULL && sepadata->tcliquegraphloaded )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;
      
   /* load tclique data structure */
   if( !sepadata->tcliquegraphloaded )
   {
      assert(sepadata->tcliquegraph == NULL);

      CHECK_OKAY( loadTcliquegraph(scip, sepadata) );
      sepadata->tcliquegraphloaded = TRUE;
   
      if( sepadata->tcliquegraph == NULL )
      {
         /* we did not find any variables that are contained in a clique with at least 3 variables in the
          * implication graph or in the clique table -> nothing has to be done
          */
         return SCIP_OKAY;
      }
   }
   tcliquegraph = sepadata->tcliquegraph;
   assert(tcliquegraph != NULL);
   
   /* store LP-solution in sepadata and update weights in tclique data structure */
   CHECK_OKAY( SCIPallocBufferArray(scip, &sepadata->varsolvals, tcliquegraph->nnodes) );
   CHECK_OKAY( SCIPgetVarSols(scip, tcliquegraph->nnodes, tcliquegraph->vars, sepadata->varsolvals) );
   updateTcliquegraph(scip, sepadata);

   /* get maximal number of tree nodes */
   maxtreenodes = (sepadata->maxtreenodes == -1 ? INT_MAX : sepadata->maxtreenodes);

   /* finds maximum weight clique in tclique */
   CHECK_OKAY( SCIPallocBufferArray(scip, &cliquenodes, tcliquegraph->nnodes) );
   tcliqueMaxClique(tcliqueGetnnodesClique, tcliqueGetweightsClique, tcliqueIsedgeClique, tcliqueSelectadjnodesClique, 
      tcliquegraph, tcliqueNewsolClique, (TCLIQUEDATA*)sepadata,
      cliquenodes, &ncliquenodes, &cliqueweight, (int)sepadata->scaleval-1, (int)sepadata->scaleval+1, maxtreenodes);

   /* frees data structures */
   SCIPfreeBufferArray(scip, &cliquenodes);
   SCIPfreeBufferArray(scip, &sepadata->varsolvals);

   /* adjust result code */
   if( sepadata->ncuts > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}




/*
 * separator specific interface methods
 */

/** creates the clique separator and includes it in SCIP */
RETCODE SCIPincludeSepaClique(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   SEPADATA* sepadata;

   /* create clique separator data */
   CHECK_OKAY( SCIPallocMemory(scip, &sepadata) );
   sepadata->tcliquegraph = NULL;
   sepadata->scip = scip;
   sepadata->varsolvals = NULL; 
   sepadata->ncuts = 0;
   sepadata->tcliquegraphloaded = FALSE;

   /* include separator */
   CHECK_OKAY( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_DELAY,
         sepaFreeClique, sepaInitClique, sepaExitClique, sepaInitsolClique, sepaExitsolClique, sepaExecClique,
         sepadata) );

   /* add clique separator parameters */
   CHECK_OKAY( SCIPaddRealParam(scip,
         "separating/clique/scaleval",
         "factor for scaling weights",
         &sepadata->scaleval, DEFAULT_SCALEVAL, 1.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/clique/maxtreenodes",
         "maximal number of nodes in branch and bound tree (-1: no limit)",
         &sepadata->maxtreenodes, DEFAULT_MAXTREENODES, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

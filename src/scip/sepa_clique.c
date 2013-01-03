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

/**@file   sepa_clique.c
 * @brief  clique separator
 * @author Kati Wolter
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_clique.h"
#include "tclique/tclique.h"
#include "scip/pub_misc.h"


#define SEPA_NAME              "clique"
#define SEPA_DESC              "clique separator of stable set relaxation"
#define SEPA_PRIORITY             -5000
#define SEPA_FREQ                     0
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_SCALEVAL         1000.0 /**< factor for scaling weights */
#define DEFAULT_MAXTREENODES      10000 /**< maximal number of nodes in branch and bound tree (-1: no limit) */
#define DEFAULT_BACKTRACKFREQ      1000 /**< frequency for premature backtracking up to tree level 1 (0: no backtracking) */
#define DEFAULT_MAXSEPACUTS          10 /**< maximal number of clique cuts separated per separation round (-1: no limit) */
#define DEFAULT_MAXZEROEXTENSIONS  1000 /**< maximal number of zero-valued variables extending the clique (-1: no limit) */
#define DEFAULT_CLIQUETABLEMEM  20000.0 /**< maximal memory size of dense clique table (in kb) */
#define DEFAULT_CLIQUEDENSITY      0.05 /**< minimal density of cliques to use a dense clique table */


/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   TCLIQUE_GRAPH*        tcliquegraph;       /**< tclique graph data structure */
   SCIP*                 scip;               /**< SCIP data structure */
   SCIP_SEPA*            sepa;               /**< separator */
   SCIP_SOL*             sol;                /**< primal solution that is currently separated */
   SCIP_Real*            varsolvals;         /**< LP solution of binary variables (contained in a 3-clique in implgraph) */
   SCIP_Real             scaleval;           /**< factor for scaling weights */
   SCIP_Longint          ncalls;             /**< number of calls to the clique separator */
   int                   maxtreenodes;       /**< maximal number of nodes in branch and bound tree (-1: no limit) */
   int                   backtrackfreq;      /**< frequency for premature backtracking up to tree level 1 (0: no backtracking) */
   int                   maxsepacuts;        /**< maximal number of clique cuts separated per separation round (-1: no limit) */
   int                   maxzeroextensions;  /**< maximal number of zero-valued variables extending the clique (-1: no limit) */
   SCIP_Real             cliquetablemem;     /**< maximal memory size of dense clique table (in kb) */
   SCIP_Real             cliquedensity;      /**< minimal density of cliques to use a dense clique table */
   int                   ncuts;              /**< number of cuts found */
   SCIP_Bool             tcliquegraphloaded; /**< TRUE if tcliquegraph is already loaded (tcliquegraph can be NULL),
                                              *   FALSE otherwise */
   SCIP_RETCODE          retcode;            /**< error code which might occur during the maximal clique algorithm */
};

/** tclique graph data */
struct TCLIQUE_Graph
{
   SCIP_VAR**            vars;               /**< active problem variables (or negated variables) the nodes belong to */
   TCLIQUE_WEIGHT*       weights;            /**< weight of nodes */
   int*                  adjnodesidxs;       /**< indices in adjnodes array of first adjacent nodes for each node */
   int*                  cliqueidsidxs;      /**< indices in cliqueids array of first clique the node is contained in */
   int*                  adjnodes;           /**< adjacent nodes of edges */
   int*                  cliqueids;          /**< unique ids of cliques */
   unsigned int*         cliquetable;        /**< dense bitvector clique table (array stored as a vector) */
   int                   adjnodessize;       /**< size of adjnodes array */
   int                   cliqueidssize;      /**< size of cliqueids array */
   int                   nnodes;             /**< number of nodes in graph */
   int                   tablewidth;         /**< number of unsigned ints per row in the table */
};


/*
 * Methods for tclique graph
 */

/** creates an empty tclique graph data structure */
static
SCIP_RETCODE tcliquegraphCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   TCLIQUE_GRAPH**       tcliquegraph        /**< pointer to tclique graph data */
   )
{
   int maxnnodes;

   assert(tcliquegraph != NULL);

   SCIP_CALL( SCIPallocMemory(scip, tcliquegraph) );

   /* there are at most 2*nbinvars nodes in the graph */
   maxnnodes = 2*SCIPgetNBinVars(scip);
   assert(maxnnodes > 0);

   /* allocate memory for tclique graph arrays */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*tcliquegraph)->vars, maxnnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*tcliquegraph)->weights, maxnnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*tcliquegraph)->adjnodesidxs, maxnnodes+1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*tcliquegraph)->cliqueidsidxs, maxnnodes+1) );
   (*tcliquegraph)->adjnodesidxs[0] = 0;  /* the last slot defines the end of the last node */
   (*tcliquegraph)->cliqueidsidxs[0] = 0; /* the last slot defines the end of the last node */
   (*tcliquegraph)->adjnodes = NULL;
   (*tcliquegraph)->cliqueids = NULL;
   (*tcliquegraph)->cliquetable = NULL;
   (*tcliquegraph)->adjnodessize = 0;
   (*tcliquegraph)->cliqueidssize = 0;
   (*tcliquegraph)->nnodes = 0;
   (*tcliquegraph)->tablewidth = 0;

   return SCIP_OKAY;
}

/** frees the tclique graph data structure and releases all captured variables */
static
SCIP_RETCODE tcliquegraphFree(
   SCIP*                 scip,               /**< SCIP data structure */
   TCLIQUE_GRAPH**       tcliquegraph        /**< pointer to tclique graph data */
   )
{
   int v;

   assert(tcliquegraph != NULL);
   assert(*tcliquegraph != NULL);

   /* release variables */
   for( v = 0; v < (*tcliquegraph)->nnodes; ++v )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*tcliquegraph)->vars[v]) );
   }

   /* free memory */
   SCIPfreeMemoryArray(scip, &(*tcliquegraph)->vars);
   SCIPfreeMemoryArray(scip, &(*tcliquegraph)->weights);
   SCIPfreeMemoryArray(scip, &(*tcliquegraph)->adjnodesidxs);
   SCIPfreeMemoryArray(scip, &(*tcliquegraph)->cliqueidsidxs);
   SCIPfreeMemoryArrayNull(scip, &(*tcliquegraph)->adjnodes);
   SCIPfreeMemoryArrayNull(scip, &(*tcliquegraph)->cliqueids);
   SCIPfreeMemoryArrayNull(scip, &(*tcliquegraph)->cliquetable);
   SCIPfreeMemory(scip, tcliquegraph);

   return SCIP_OKAY;
}

/** ensures that the adjnodes array can store at least num entries */
static
SCIP_RETCODE tcliquegraphEnsureAdjnodesSize(
   SCIP*                 scip,               /**< SCIP data structure */
   TCLIQUE_GRAPH*        tcliquegraph,       /**< tclique graph data */
   int                   num                 /**< minimal number of adjacent nodes to be able to store in the array */
   )
{
   assert(tcliquegraph != NULL);

   if( num > tcliquegraph->adjnodessize )
   {
      tcliquegraph->adjnodessize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &tcliquegraph->adjnodes, tcliquegraph->adjnodessize) );
   }
   assert(num <= tcliquegraph->adjnodessize);

   return SCIP_OKAY;
}

/** ensures that the cliqueids array can store at least num entries */
static
SCIP_RETCODE tcliquegraphEnsureCliqueidsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   TCLIQUE_GRAPH*        tcliquegraph,       /**< tclique graph data */
   int                   num                 /**< minimal number of adjacent nodes to be able to store in the array */
   )
{
   assert(tcliquegraph != NULL);

   if( num > tcliquegraph->cliqueidssize )
   {
      tcliquegraph->cliqueidssize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &tcliquegraph->cliqueids, tcliquegraph->cliqueidssize) );
   }
   assert(num <= tcliquegraph->cliqueidssize);

   return SCIP_OKAY;
}

/** adds a node to the tclique graph defined as a variable-value pair; adds all cliques to the cliqueids array the
 *  variable is contained in with the given value
 */
static
SCIP_RETCODE tcliquegraphAddNode(
   SCIP*                 scip,               /**< SCIP data structure */
   TCLIQUE_GRAPH**       tcliquegraph,       /**< pointer to tclique graph data */
   SCIP_VAR*             var,                /**< active binary problem variable */
   SCIP_Bool             value,              /**< value of the variable in the node */
   int*                  nodeidx             /**< pointer to store the index of the new node */
   )
{
   SCIP_VAR* nodevar;
   int* cliqueids;
   SCIP_CLIQUE** cliques;
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
      SCIP_CALL( tcliquegraphCreate(scip, tcliquegraph) );
   }
   assert(*tcliquegraph != NULL);
   assert((*tcliquegraph)->nnodes < 2*SCIPgetNBinVars(scip));

   /* if the value is FALSE, use the negated variable for the node */
   if( !value )
   {
      SCIP_CALL( SCIPgetNegatedVar(scip, var, &nodevar) );
   }
   else
      nodevar = var;

   /* get the current number of used entries in adjnodes and cliqueids arrays */
   nadjnodes = (*tcliquegraph)->adjnodesidxs[(*tcliquegraph)->nnodes];
   ncliqueids = (*tcliquegraph)->cliqueidsidxs[(*tcliquegraph)->nnodes];

   /* insert the variable into the tclique graph */
   *nodeidx = (*tcliquegraph)->nnodes;
   SCIP_CALL( SCIPcaptureVar(scip, nodevar) );
   (*tcliquegraph)->vars[*nodeidx] = nodevar;
   (*tcliquegraph)->weights[*nodeidx] = 0;
   (*tcliquegraph)->nnodes++;

   /* store the ids of the variable's cliques in the cliqueids array */
   ncliques = SCIPvarGetNCliques(var, value);
   cliques = SCIPvarGetCliques(var, value);
   SCIP_CALL( tcliquegraphEnsureCliqueidsSize(scip, *tcliquegraph, ncliqueids + ncliques) );
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
SCIP_RETCODE tcliquegraphAddCliqueVars(
   SCIP*                 scip,               /**< SCIP data structure */
   TCLIQUE_GRAPH**       tcliquegraph,       /**< pointer to tclique graph data */
   int**                 cliquegraphidx      /**< array to store tclique graph node index of variable/value pairs */
   )
{
   SCIP_VAR** vars;
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
      SCIP_VAR* var;
      int value;

      var = vars[i];

      for( value = 0; value < 2; ++value )
      {
         assert(cliquegraphidx[value][i] == -1);
         
         if( SCIPvarGetNCliques(var, (SCIP_Bool)value) >= 1 )
         {
            /* all cliques stored in the clique table are at least 3-cliques */
            SCIP_CALL( tcliquegraphAddNode(scip, tcliquegraph, var, (SCIP_Bool)value, &cliquegraphidx[value][i]) );
         }
      }
   }

   return SCIP_OKAY;
}

/* adds all variable/value pairs to the tclique graph that are contained in a 3-clique in the implication graph */
static
SCIP_RETCODE tcliquegraphAddImplicsVars(
   SCIP*                 scip,               /**< SCIP data structure */
   TCLIQUE_GRAPH**       tcliquegraph,       /**< pointer to tclique graph data */
   int**                 cliquegraphidx      /**< array to store tclique graph node index of variable/value pairs */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int xi;

   assert(tcliquegraph != NULL);
   assert(cliquegraphidx != NULL);
   assert(cliquegraphidx[0] != NULL);
   assert(cliquegraphidx[1] != NULL);

   /**@todo search for 3-clique in the implication graph w.r.t. to implicit binary variable (SCIPvarIsBinary()) */
   /* get binary variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNBinVars(scip);

   if( nvars == 0 )
      return SCIP_OKAY;

   assert(nvars > 0);

   /* detect 3-cliques in the clique graph: triples (x,y,z) in the implication graph with x -> y', x -> z', and y -> z';
    * in order to avoid double checks, we only check triplets with variable indices xindex < yindex < zindex;
    * we don't have to check triples with x == y, x == z, or y == z, because cliques with the same variable occuring
    * twice lead to fixings of this or all other variables which should already have been detected in presolving
    */
   for( xi = 0; xi < nvars-2 && !SCIPisStopped(scip); ++xi ) /* at least two variables must be left over for y and z */
   {
      SCIP_VAR* x;
      int xvalue;
#ifndef NDEBUG
      int xindex;
#endif

      x = vars[xi];
#ifndef NDEBUG
      xindex = SCIPvarGetIndex(x);
#endif
      assert(SCIPvarGetType(x) == SCIP_VARTYPE_BINARY);

      for( xvalue = 0; xvalue < 2; xvalue++ )
      {
         SCIP_VAR** ximplvars;
         SCIP_BOUNDTYPE* ximpltypes;
         int xnbinimpls;
         int i;

         /**@todo search for 3-clique in the implication graph w.r.t. to implicit binary variable (SCIPvarIsBinary()) */
         /* scan the implications of x == xvalue for potential y candidates */
         xnbinimpls = SCIPvarGetNBinImpls(x, (SCIP_Bool)xvalue);
         ximplvars = SCIPvarGetImplVars(x, (SCIP_Bool)xvalue);
         ximpltypes = SCIPvarGetImplTypes(x, (SCIP_Bool)xvalue);

	 assert(xnbinimpls >= 0);
	 if( xnbinimpls > 0 )
	 {
#ifndef NDEBUG
	    SCIP_Bool found;

	    /* search for the position in the sorted array (via binary search) */
	    found = SCIPsortedvecFindPtr((void**)ximplvars, SCIPvarComp, (void*)x, xnbinimpls, &i);
	    assert(!found); /* no implication to itself */
	    assert(i >= 0 && i <= xnbinimpls);
#else
	    /* search for the position in the sorted array (via binary search) */
	    (void) SCIPsortedvecFindPtr((void**)ximplvars, SCIPvarComp, (void*)x, xnbinimpls, &i);
#endif
	 }
	 else
	    i = 0;

         /* loop over all y > x */
         for( ; i < xnbinimpls-1 && !SCIPisStopped(scip); ++i ) /* at least one variable must be left over for z */
         {
            SCIP_VAR* y;
            int yvalue;
            int yi;
            SCIP_VAR** yimplvars;
            SCIP_BOUNDTYPE* yimpltypes;
            int ynbinimpls;
            int xk;
            int yk;
#ifndef NDEBUG
            int yindex;
#endif

            y = ximplvars[i];
            yi = SCIPvarGetProbindex(y);
#ifndef NDEBUG
            yindex = SCIPvarGetIndex(y);
#endif
            assert(yi < nvars);

            /* consider only implications with active implvar y */
            if( yi < 0 )
               continue;

            assert(xindex < yindex); /* the implied variables are sorted by increasing variable index */
            assert(yindex < SCIPvarGetIndex(ximplvars[i+1]));

            /* check, whether the implicant is y == 0 or y == 1 (yi conflicts with x == xvalue) */
            yvalue = (int)(ximpltypes[i] == SCIP_BOUNDTYPE_UPPER);
            assert(0 <= yvalue && yvalue <= 1);

            /**@todo search for 3-clique in the implication graph w.r.t. to implicit binary variable (SCIPvarIsBinary()) */
            /* scan the remaining implications of x == xvalue and the implications of y == yvalue for equal entries */
            ynbinimpls = SCIPvarGetNBinImpls(y, (SCIP_Bool)yvalue);
            yimplvars = SCIPvarGetImplVars(y, (SCIP_Bool)yvalue);
            yimpltypes = SCIPvarGetImplTypes(y, (SCIP_Bool)yvalue);

            /* we simultaneously scan both implication arrays for candidates z with x < y < z */
            xk = i+1;

	    assert(ynbinimpls >= 0);
	    if( ynbinimpls > 0 )
	    {
#ifndef NDEBUG
	       SCIP_Bool found;

	       /* search for the position in the sorted array (via binary search) */
	       found = SCIPsortedvecFindPtr((void**)yimplvars, SCIPvarComp, (void*)y, ynbinimpls, &yk);
	       assert(!found); /* no implication to itself */
	       assert(yk >= 0 && yk <= ynbinimpls);
#else
	       (void) SCIPsortedvecFindPtr((void**)yimplvars, SCIPvarComp, (void*)y, ynbinimpls, &yk);
#endif
	    }
	    else
	       yk = 0;

            while( xk < xnbinimpls && yk < ynbinimpls )
            {
               int zindex;

               assert(yindex < SCIPvarGetIndex(ximplvars[xk]));
               assert(yindex < SCIPvarGetIndex(yimplvars[yk]));

               /* scan the implications of x */
               zindex = SCIPvarGetIndex(yimplvars[yk]);

               while ( xk < xnbinimpls && SCIPvarGetIndex(ximplvars[xk]) < zindex )
                  ++xk;

               if( xk >= xnbinimpls )
                  break;

               /* scan the implications of y */
               zindex = SCIPvarGetIndex(ximplvars[xk]);

               while ( yk < ynbinimpls && SCIPvarGetIndex(yimplvars[yk]) < zindex )
                  ++yk;

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
                  SCIP_VAR* z;
                  int zi;
                  int zvalue;
                  
                  /* we found z with xindex < yindex < zindex and x + y + z <= 1 */
                  z = ximplvars[xk];
                  zi = SCIPvarGetProbindex(z);
                  assert(zi < nvars);
                  
                  /* consider only implications with active implvar z */
                  if( zi < 0 )
                     continue;
                  
                  assert(SCIPvarGetIndex(z) == zindex);
                  assert(xindex < yindex && yindex < zindex);

                  /* check, whether the implicant is z == 0 or z == 1 (z conflicts with x == xvalue) */
                  zvalue = (int)(ximpltypes[xk] == SCIP_BOUNDTYPE_UPPER);
                  assert(0 <= zvalue && zvalue <= 1);

                  /* add nodes x == xvalue, y == yvalue, and z == zvalue to clique graph (if not yet existing) */
                  if( cliquegraphidx[xvalue][xi] == -1 )
                  {
                     SCIP_CALL( tcliquegraphAddNode(scip, tcliquegraph, x, (SCIP_Bool)xvalue, &cliquegraphidx[xvalue][xi]) );
                  }
                  if( cliquegraphidx[yvalue][yi] == -1 )
                  {
                     SCIP_CALL( tcliquegraphAddNode(scip, tcliquegraph, y, (SCIP_Bool)yvalue, &cliquegraphidx[yvalue][yi]) );
                  }
                  if( cliquegraphidx[zvalue][zi] == -1 )
                  {
                     SCIP_CALL( tcliquegraphAddNode(scip, tcliquegraph, z, (SCIP_Bool)zvalue, &cliquegraphidx[zvalue][zi]) );
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
SCIP_RETCODE tcliquegraphAddImplicsCliqueVars(
   SCIP*                 scip,               /**< SCIP data structure */
   TCLIQUE_GRAPH**       tcliquegraph,       /**< pointer to tclique graph data */
   int**                 cliquegraphidx      /**< array to store tclique graph node index of variable/value pairs */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int xi;

   assert(tcliquegraph != NULL);
   assert(cliquegraphidx != NULL);
   assert(cliquegraphidx[0] != NULL);
   assert(cliquegraphidx[1] != NULL);

   /**@todo search for 3-clique in the implication graph w.r.t. to implicit binary variable (SCIPvarIsBinary()) */
   /* get binary variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNBinVars(scip);
   assert(nvars > 0);

   /* detect triples (x,y,z) with implications x -> y' and x -> z' and a clique that contains y and z;
    * because all cliques stored in the clique table are at least 3-cliques, y and z are already member of the
    * tclique graph due to tcliquegraphAddCliqueVars(); therefore, we only have to process pairs x == xvalue
    * that are currently not member of the tclique graph
    */
   for( xi = 0; xi < nvars && !SCIPisStopped(scip); ++xi )
   {
      SCIP_VAR* x;
      SCIP_Bool xvalue;

      x = vars[xi];

      for( xvalue = 0; xvalue < 2; xvalue++ )
      {
         SCIP_VAR** ximplvars;
         SCIP_BOUNDTYPE* ximpltypes;
         int xnbinimpls;
         int yk;

         /* we only need to process pairs x == xvalue that are not yet member of the tclique graph */
         if( cliquegraphidx[xvalue][xi] >= 0 )
            continue;

         /**@todo search for 3-clique in the implication graph w.r.t. to implicit binary variable (SCIPvarIsBinary()) */
         /* scan the implications of x == xvalue for potential y and z candidates */
         xnbinimpls = SCIPvarGetNBinImpls(x, (SCIP_Bool)xvalue);
         ximplvars = SCIPvarGetImplVars(x, (SCIP_Bool)xvalue);
         ximpltypes = SCIPvarGetImplTypes(x, (SCIP_Bool)xvalue);

         /* loop over all y, at least one variable must be left over for z */
         for( yk = 0; yk < xnbinimpls-1 && !SCIPisStopped(scip); ++yk ) 
         {
            SCIP_VAR* y;
            SCIP_Bool yvalue;
            int yi;
            int zk;

            y = ximplvars[yk];
            yi = SCIPvarGetProbindex(y);
            assert(yi < nvars);

            /* consider only implications with active implvar y */
            if( yi < 0 )
               continue;

            assert(SCIPvarGetType(y) == SCIP_VARTYPE_BINARY);
            
            /* check, whether the implicant is y == 0 or y == 1 (y conflicts with x == xvalue) */
            yvalue = (ximpltypes[yk] == SCIP_BOUNDTYPE_UPPER);
            if( SCIPvarGetNCliques(y, yvalue) == 0 )
               continue;

            /* loop over all z > y */
            for( zk = yk+1; zk < xnbinimpls; ++zk )
            {
               SCIP_VAR* z;
               SCIP_Bool zvalue;
               int zi;

               z = ximplvars[zk];
               zi = SCIPvarGetProbindex(z);
               assert(zi < nvars);

               /* consider only implications with active implvar z */
               if( zi < 0 )
                  continue;

               assert(SCIPvarGetType(z) == SCIP_VARTYPE_BINARY);
               
               /* check, whether the implicant is z == 0 or z == 1 (z conflicts with x == xvalue) */
               zvalue = (ximpltypes[zk] == SCIP_BOUNDTYPE_UPPER);

               /* check whether y and z have a common clique stored in the clique table
                * (the ones in the implication graph have already been processed by tcliquegraphAddImplicsVars())
                */
               if( SCIPvarsHaveCommonClique(y, yvalue, z, zvalue, FALSE) )
               {
                  /* add nodes x == xvalue to clique graph */
                  assert(cliquegraphidx[xvalue][xi] == -1);
                  assert(cliquegraphidx[yvalue][yi] >= 0);
                  assert(cliquegraphidx[zvalue][zi] >= 0);
                  SCIP_CALL( tcliquegraphAddNode(scip, tcliquegraph, x, (SCIP_Bool)xvalue, &cliquegraphidx[xvalue][xi]) );
                  assert(cliquegraphidx[xvalue][xi] >= 0);
                  break;
               }
            }
            if( cliquegraphidx[xvalue][xi] >= 0 )
               break;
         }
      }
   }

   return SCIP_OKAY;
}

/** adds all implications between nodes in the tclique graph to the tclique graph */
static
SCIP_RETCODE tcliquegraphAddImplics(
   SCIP*                 scip,               /**< SCIP data structure */
   TCLIQUE_GRAPH*        tcliquegraph,       /**< tclique graph data */
   int**                 cliquegraphidx      /**< array to store tclique graph node index of variable/value pairs */
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
      SCIP_VAR* var;
      SCIP_Bool value;
      SCIP_VAR** implvars;
      SCIP_BOUNDTYPE* impltypes;
      int nbinimpls;
      int adjnodesidx;
      int j;

      /* store start adjnodes index of the node */
      adjnodesidx = nadjnodes;
      tcliquegraph->adjnodesidxs[i] = adjnodesidx;

      /* get active problem variable */
      var = tcliquegraph->vars[i];
      value = TRUE;
      SCIP_CALL( SCIPvarGetProbvarBinary(&var, &value) );
      assert(0 <= SCIPvarGetProbindex(var) && SCIPvarGetProbindex(var) < SCIPgetNBinVars(scip));
      assert(cliquegraphidx[value][SCIPvarGetProbindex(var)] == i);

      /**@todo search for 3-clique in the implication graph w.r.t. to implicit binary variable (SCIPvarIsBinary()) */
      /* get implications on binary variables */
      nbinimpls = SCIPvarGetNBinImpls(var, value);
      implvars = SCIPvarGetImplVars(var, value);
      impltypes = SCIPvarGetImplTypes(var, value);
      for( j = 0; j < nbinimpls; ++j )
      {
         int probidx;
         int graphidx;
         SCIP_Bool implvalue;

         probidx = SCIPvarGetProbindex(implvars[j]);
         implvalue = (impltypes[j] == SCIP_BOUNDTYPE_UPPER);
         assert(probidx < SCIPgetNVars(scip));
                  
         /* consider only implications with active implvar */
         if( probidx < 0 )
            continue;

         graphidx = cliquegraphidx[implvalue][probidx];
         if( graphidx >= 0 )
         {
            int pos;

            assert(graphidx < tcliquegraph->nnodes);
            assert((implvalue == TRUE && tcliquegraph->vars[graphidx] == implvars[j])
               || (implvalue == FALSE && SCIPvarGetNegationVar(tcliquegraph->vars[graphidx]) == implvars[j]));

            /* allocate memory for additional arc */
            SCIP_CALL( tcliquegraphEnsureAdjnodesSize(scip, tcliquegraph, nadjnodes+1) );

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

/** constructs dense clique incidence matrix */
static
SCIP_RETCODE tcliquegraphConstructCliqueTable(
   SCIP*                 scip,               /**< SCIP data structure */
   TCLIQUE_GRAPH*        tcliquegraph,       /**< tclique graph data */
   SCIP_Real             cliquetablemem,     /**< maximal memory size of dense clique table (in kb) */
   SCIP_Real             cliquedensity       /**< minimal density of cliques to store as dense table */
   )
{
   SCIP_CLIQUE** cliques;
   int* varids;
   unsigned int* cliquetable;
   SCIP_Real density;
   int nbits;
   int tablesize;
   int tablewidth;
   int ncliques;
   int nelems;
   int i;

   /* get clique table */
   cliques = SCIPgetCliques(scip);
   ncliques = SCIPgetNCliques(scip);
   if( ncliques == 0 )
      return SCIP_OKAY;

   assert(tcliquegraph != NULL);

   /* calculate size of dense clique table */
   nbits = 8*sizeof(unsigned int);
   tcliquegraph->tablewidth = (tcliquegraph->nnodes + nbits-1) / nbits; /* number of ints needed */
   
   /* check if dense clique table is too large (calculate as Reals to avoid overflow) */
   if( (SCIP_Real)tcliquegraph->nnodes * (SCIP_Real)tcliquegraph->tablewidth/1024.0 > cliquetablemem )
      return SCIP_OKAY;

   /* calculate clique entry density */
   nelems = 0;
   for( i = 0; i < ncliques; ++i )
      nelems += SCIPcliqueGetNVars(cliques[i]);
   density = (SCIP_Real)nelems / ((SCIP_Real)ncliques * (SCIP_Real)tcliquegraph->nnodes);
   if( density < cliquedensity )
      return SCIP_OKAY;

   /* allocate memory */
   tablesize = tcliquegraph->nnodes * tcliquegraph->tablewidth;
   SCIPdebugMessage("clique separator: constructing dense clique table (%d kb, %d cliques, %d nodes, density: %.2f)\n",
      tablesize/1024, SCIPgetNCliques(scip), tcliquegraph->nnodes, density);

   SCIP_CALL( SCIPallocMemoryArray(scip, &tcliquegraph->cliquetable, tablesize) );
   BMSclearMemoryArray(tcliquegraph->cliquetable, tablesize);

   /* insert the cliques as complete graphs to the incidence matrix */
   SCIP_CALL( SCIPallocBufferArray(scip, &varids, tcliquegraph->nnodes) );
   cliquetable = tcliquegraph->cliquetable;
   tablewidth = tcliquegraph->tablewidth;
   for( i = 0; i < ncliques && !SCIPisStopped(scip); ++i )
   {
      SCIP_VAR** vars;
      SCIP_Bool* vals;
      int nvars;
      int u;
      int v;

      vars = SCIPcliqueGetVars(cliques[i]);
      vals = SCIPcliqueGetValues(cliques[i]);
      nvars = SCIPcliqueGetNVars(cliques[i]);
#if 0  /**@todo this assert is currently not valid since implicit binary variables in cliques are ignored, 
        * i.e., corresponding nodes and edges are not added to the tclique graph. Enable assert again if 
        * this feature it incorporated. 
        */
      assert(nvars <= tcliquegraph->nnodes);
#endif

      /* get the node numbers of the variables */
      for( u = 0; u < nvars && !SCIPisStopped(scip); ++u )
      {
         SCIP_VAR* var;

         if( SCIPvarGetType(vars[u]) != SCIP_VARTYPE_BINARY )
            continue;

         var = (vals[u] ? vars[u] : SCIPvarGetNegatedVar(vars[u]));
         assert(var != NULL); /* var must exist even if negated, since it is stored in the tcliquegraph */
         for( v = 0; v < tcliquegraph->nnodes && var != tcliquegraph->vars[v]; ++v )
         {}
         assert(v < tcliquegraph->nnodes);
         varids[u] = v;
      }

      /* flag the edges in the incidence matrix (excluding diagonal entries) */
      for( u = 0; u < nvars-1 && !SCIPisStopped(scip); ++u )
      {
         int nu;
         int rowstart;
         int colofs;
         unsigned int colmask;

         if( SCIPvarGetType(vars[u]) != SCIP_VARTYPE_BINARY )
            continue;

         nu = varids[u];
         rowstart = nu*tablewidth;
         colofs = nu/nbits;
         colmask = 1 << (nu % nbits); /*lint !e701*/
         for( v = u+1; v < nvars; ++v )
         {
            int nv;
            unsigned int mask;

            if( SCIPvarGetType(vars[v]) != SCIP_VARTYPE_BINARY )
               continue;

            nv = varids[v];
            mask = 1 << (nv % nbits); /*lint !e701*/
            cliquetable[rowstart+nv/nbits] |= mask;
            cliquetable[nv*tablewidth+colofs] |= colmask;
         }
      }
   }
   SCIPfreeBufferArray(scip, &varids);

   SCIPdebugMessage("clique separator: finished constructing dense clique table\n");

   return SCIP_OKAY;
}

/* creates tclique data structure using the implication graph;
 * only variables that are contained in a 3-clique are added as nodes to the clique graph
 */
static
SCIP_RETCODE loadTcliquegraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< separator data */
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
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquegraphidx[0], nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquegraphidx[1], nvars) );
   for( i = 0; i < nvars; ++i )
   {
      cliquegraphidx[0][i] = -1;
      cliquegraphidx[1][i] = -1;
   }

   /* insert all variable/value pairs that are contained in an existing 3-clique */
   SCIP_CALL( tcliquegraphAddCliqueVars(scip, &sepadata->tcliquegraph, cliquegraphidx) );

   /* insert all variable/value pairs that are contained in a 3-clique in the implication graph */
   SCIP_CALL( tcliquegraphAddImplicsVars(scip, &sepadata->tcliquegraph, cliquegraphidx) );

   /* insert all variable/value pairs that have implications to two variables of the same existing clique */
   SCIP_CALL( tcliquegraphAddImplicsCliqueVars(scip, &sepadata->tcliquegraph, cliquegraphidx) );

   /* add all implications between used variables to the tclique graph */
   SCIP_CALL( tcliquegraphAddImplics(scip, sepadata->tcliquegraph, cliquegraphidx) );

   /* it occurs that it might be that some cliques were not yet removed from the global clique array, so SCIPgetNClique
    * can be greater than 0, even if there is no clique with some variables left
    *
    * @todo clean up empty cliques
    */
   if( sepadata->tcliquegraph != NULL )
   {
      /* construct the dense clique table */
      SCIP_CALL( tcliquegraphConstructCliqueTable(scip, sepadata->tcliquegraph, sepadata->cliquetablemem, sepadata->cliquedensity) );
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &cliquegraphidx[1]);
   SCIPfreeBufferArray(scip, &cliquegraphidx[0]);
   if( SCIPisStopped(scip) && sepadata->tcliquegraph != NULL )
      SCIP_CALL( tcliquegraphFree(scip,&sepadata->tcliquegraph) );
   return SCIP_OKAY;
}

/* updates the weights in the tclique graph data structure */
static
void updateTcliquegraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< separator data */
   )
{
   TCLIQUE_GRAPH* tcliquegraph;
   int i;

   assert(sepadata != NULL);
   assert(sepadata->varsolvals != NULL);

   tcliquegraph = sepadata->tcliquegraph;
   assert(tcliquegraph != NULL);

   /* updates weight of all nodes in tclique data structure */
   for( i = 0; i < tcliquegraph->nnodes; i++ )
   {
      int weight;

      weight = (TCLIQUE_WEIGHT)SCIPfeasFloor(scip, sepadata->varsolvals[i] * sepadata->scaleval);
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
SCIP_Bool nodesHaveCommonClique(
   TCLIQUE_GRAPH*        tcliquegraph,       /**< tclique graph data */
   int                   node1,              /**< first node */
   int                   node2               /**< second node */
   )
{
   assert(tcliquegraph != NULL);

   /* return TRUE for equal nodes */
   if( node1 == node2 )
      return TRUE;

   /* check whether the dense clique table was constructed */
   if( tcliquegraph->cliquetable != NULL )
   {
      int nbits;
      unsigned int mask;
      int colofs;

      /* check entry in the table */
      nbits = 8*sizeof(unsigned int);
      mask = (1 << (node2 % nbits)); /*lint !e701*/
      colofs = node2 / nbits;
      assert(((tcliquegraph->cliquetable[node1*tcliquegraph->tablewidth + colofs] & mask) != 0)
         == ((tcliquegraph->cliquetable[node2*tcliquegraph->tablewidth + node1/nbits] & (1 << (node1 % nbits))) != 0)); /*lint !e701*/
      return ((tcliquegraph->cliquetable[node1*tcliquegraph->tablewidth + colofs] & mask) != 0);
   }
   else
   {
      int* cliqueids;
      int i1;
      int i2;
      int endi1;
      int endi2;
      
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

/** basic code for new cliques (needed because of error handling) */
static
SCIP_RETCODE newsolCliqueAddRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< the cut separator itself */
   SCIP_SEPADATA*        sepadata,           /**< data of separator */
   int                   ncliquenodes,       /**< number of nodes in clique */
   int*                  cliquenodes         /**< nodes in clique */
   )
{
   SCIP_VAR** vars;
   SCIP_ROW* cut;
   char cutname[SCIP_MAXSTRLEN];
   int i;

   vars = sepadata->tcliquegraph->vars;
   assert(sepadata->tcliquegraph->nnodes > 0);
   assert(vars != NULL);

   /* create the cut (handle retcode since we do not have a backtrace) */
   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "clique%"SCIP_LONGINT_FORMAT"_%d", sepadata->ncalls, sepadata->ncuts);
   SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );

   SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

   assert(ncliquenodes <= sepadata->tcliquegraph->nnodes);
   /*SCIPdebugMessage(" -> clique in graph:");*/
   for( i = 0; i < ncliquenodes; ++i )
   {
      assert(cliquenodes[i] < sepadata->tcliquegraph->nnodes);
      SCIP_CALL( SCIPaddVarToRow(scip, cut, vars[cliquenodes[i]], 1.0) );
      /*SCIPdebugPrintf(" [%d]<%s>", cliquenodes[i], SCIPvarGetName(vars[cliquenodes[i]]));*/
   }
   /*SCIPdebugPrintf("\n");*/
   SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

   /*SCIPdebug( SCIP_CALL(SCIPprintRow(scip, cut, NULL)) );*/

   SCIP_CALL( SCIPaddCut(scip, sepadata->sol, cut, FALSE) );
   SCIP_CALL( SCIPaddPoolCut(scip, cut) );

   /* release the row */
   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   return SCIP_OKAY;
}

/** generates cuts using a clique found by algorithm for maximum weight clique
 *  and decides whether to stop generating cliques with the algorithm for maximum weight clique
 */
static
TCLIQUE_NEWSOL(tcliqueNewsolClique)
{
   SCIP_SEPADATA* sepadata;
   TCLIQUE_WEIGHT minweightinc;

   assert(acceptsol != NULL);
   assert(stopsolving != NULL);

   sepadata = (SCIP_SEPADATA*)tcliquedata;
   assert(sepadata != NULL);
   assert(sepadata->scip != NULL);
   assert(sepadata->sepa != NULL);
   assert(sepadata->tcliquegraph != NULL);
   assert(sepadata->ncuts >= 0);

   /* we don't accept the solution as new incumbent, because we want to find many violated clique inequalities */
   *acceptsol = FALSE;
   *stopsolving = FALSE;

   /* slightly increase the minimal weight for additional cliques */
   minweightinc = (cliqueweight - *minweight)/10;
   minweightinc = MAX(minweightinc, 1);
   *minweight += minweightinc;

   /* adds cut if weight of the clique is greater than 1 */
   if( cliqueweight > sepadata->scaleval )
   {
      SCIP* scip;
      SCIP_SEPA* sepa;
      SCIP_Real* varsolvals;
      SCIP_Real unscaledweight;
      int i;

      scip = sepadata->scip;
      sepa = sepadata->sepa;
      varsolvals = sepadata->varsolvals;
      assert(varsolvals != NULL);

      /* calculate the weight of the clique in unscaled fractional variable space */
      unscaledweight = 0.0;
      for( i = 0; i < ncliquenodes; i++ )
         unscaledweight += varsolvals[cliquenodes[i]];

      if( SCIPisEfficacious(scip, unscaledweight - 1.0) )
      {
         SCIP_RETCODE retcode;

         /* explicitly handle return code */
         retcode = newsolCliqueAddRow(scip, sepa, sepadata, ncliquenodes, cliquenodes);
         if ( retcode == SCIP_OKAY )
         {
            SCIPdebugMessage(" -> found clique cut (act=%g)\n", unscaledweight);
            sepadata->ncuts++;

            /* if we found more than half the cuts we are allowed to generate, we accept the clique as new incumbent,
             * such that only more violated cuts are generated afterwards
             */
            if( sepadata->maxsepacuts >= 0 )
            {
               if( sepadata->ncuts > sepadata->maxsepacuts/2 )
                  *acceptsol = TRUE;
               if( sepadata->ncuts >= sepadata->maxsepacuts )
                  *stopsolving = TRUE;
            }
         }
         else
         {
            /* in case an internal SCIP error occurred we stop the algorithm and store the error code for later
             * evaluation
             */
            sepadata->retcode = retcode;
            *stopsolving = TRUE;
         }
      }
   }
}


/*
 * main separation method
 */

/** searches and adds clique cuts that separate the given primal solution */
static
SCIP_RETCODE separateCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< the cut separator itself */
   SCIP_SOL*             sol,                /**< primal solution that should be separated, or NULL for LP solution */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   TCLIQUE_GRAPH* tcliquegraph;
   int* cliquenodes;
   TCLIQUE_WEIGHT cliqueweight;
   TCLIQUE_STATUS tcliquestatus;
   int ncliquenodes;
   int maxtreenodes;
   int maxzeroextensions;

   assert(scip != NULL);
   assert(*result == SCIP_DIDNOTRUN);

   /* get separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   sepadata->sol = sol;
   sepadata->ncalls = SCIPsepaGetNCalls(sepa);
   sepadata->ncuts = 0;

   /* if we already detected that no implications between binary variables exist, nothing has to be done */
   if( sepadata->tcliquegraph == NULL && sepadata->tcliquegraphloaded )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* load tclique data structure */
   if( !sepadata->tcliquegraphloaded )
   {
      assert(sepadata->tcliquegraph == NULL);

      SCIPdebugMessage("loading implication and clique graph\n");
      SCIP_CALL( loadTcliquegraph(scip, sepadata) );
      sepadata->tcliquegraphloaded = TRUE;

      if( sepadata->tcliquegraph == NULL )
      {
         if( SCIPisStopped(scip) )
            sepadata->tcliquegraphloaded = FALSE;
         /* we did not find any variables that are contained in a clique with at least 3 variables in the
          * implication graph or in the clique table -> nothing has to be done
          */
         else
	 {
            SCIPdebugMessage("no 3-cliques found in implication graph\n");
         }

         return SCIP_OKAY;
      }
   }
   tcliquegraph = sepadata->tcliquegraph;
   assert(tcliquegraph != NULL);

   /* store LP-solution in sepadata and update weights in tclique data structure */
   SCIP_CALL( SCIPallocBufferArray(scip, &sepadata->varsolvals, tcliquegraph->nnodes) );
   SCIP_CALL( SCIPgetSolVals(scip, sol, tcliquegraph->nnodes, tcliquegraph->vars, sepadata->varsolvals) );
   updateTcliquegraph(scip, sepadata);

   /* get maximal number of tree nodes and maximal zero-extensions */
   maxtreenodes = (sepadata->maxtreenodes == -1 ? INT_MAX : sepadata->maxtreenodes);
   maxzeroextensions = (sepadata->maxzeroextensions == -1 ? INT_MAX : sepadata->maxzeroextensions);

   SCIPdebugMessage("searching for violated clique cuts\n");

   sepadata->retcode = SCIP_OKAY;

   /* finds maximum weight clique in tclique */
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquenodes, tcliquegraph->nnodes) );
   tcliqueMaxClique(tcliqueGetnnodesClique, tcliqueGetweightsClique, tcliqueIsedgeClique, tcliqueSelectadjnodesClique,
      tcliquegraph, tcliqueNewsolClique, (TCLIQUE_DATA*)sepadata,
      cliquenodes, &ncliquenodes, &cliqueweight, (int)sepadata->scaleval-1, (int)sepadata->scaleval+1,
      maxtreenodes, sepadata->backtrackfreq, maxzeroextensions, -1, NULL, &tcliquestatus);

   /* in case an internal error occurred during the maximal clique computation, evaluate that one */
   SCIP_CALL( sepadata->retcode );

   SCIPdebugMessage("finished searching clique cuts: found %d cuts\n", sepadata->ncuts);

   /* frees data structures */
   SCIPfreeBufferArray(scip, &cliquenodes);
   SCIPfreeBufferArray(scip, &sepadata->varsolvals);

   /* adjust result code */
   if( sepadata->ncuts > 0 )
      *result = SCIP_SEPARATED;

   /* better reset the sol pointer in sepadata to avoid having an invalid pointer */
   sepadata->sol = NULL;

   return SCIP_OKAY;
}


/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyClique)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaClique(scip) );
 
   return SCIP_OKAY;
}


/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeClique)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->tcliquegraph == NULL);

   SCIPfreeMemory(scip, &sepadata);
   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolClique)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* free tclique data */
   if( sepadata->tcliquegraph != NULL )
   {
      SCIP_CALL( tcliquegraphFree(scip, &sepadata->tcliquegraph) );
   }
   assert(sepadata->tcliquegraph == NULL);
   sepadata->tcliquegraphloaded = FALSE;

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpClique)
{
   /*lint --e{715}*/

   *result = SCIP_DIDNOTRUN;

   /* only call separator, if we are not close to terminating */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* only call separator, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call separator, if there are fractional variables */
   if( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   /* separate cuts on the LP solution */
   SCIP_CALL( separateCuts(scip, sepa, NULL, result) );

   return SCIP_OKAY;
}


/** arbitrary primal solution separation method of separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolClique)
{
   /*lint --e{715}*/

   *result = SCIP_DIDNOTRUN;

   /* separate cuts on the given primal solution */
   SCIP_CALL( separateCuts(scip, sepa, sol, result) );

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the clique separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaClique(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create clique separator data */
   SCIP_CALL( SCIPallocMemory(scip, &sepadata) );
   sepadata->tcliquegraph = NULL;
   sepadata->scip = scip;
   sepadata->sol = NULL;
   sepadata->varsolvals = NULL;
   sepadata->ncalls = 0;
   sepadata->ncuts = 0;
   sepadata->tcliquegraphloaded = FALSE;

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpClique, sepaExecsolClique,
         sepadata) );

   assert(sepa != NULL);
   sepadata->sepa = sepa;

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyClique) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeClique) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitsolClique) );

   /* add clique separator parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/clique/scaleval",
         "factor for scaling weights",
         &sepadata->scaleval, TRUE, DEFAULT_SCALEVAL, 1.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/clique/maxtreenodes",
         "maximal number of nodes in branch and bound tree (-1: no limit)",
         &sepadata->maxtreenodes, TRUE, DEFAULT_MAXTREENODES, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/clique/backtrackfreq",
         "frequency for premature backtracking up to tree level 1 (0: no backtracking)",
         &sepadata->backtrackfreq, TRUE, DEFAULT_BACKTRACKFREQ, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/clique/maxsepacuts",
         "maximal number of clique cuts separated per separation round (-1: no limit)",
         &sepadata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/clique/maxzeroextensions",
         "maximal number of zero-valued variables extending the clique (-1: no limit)",
         &sepadata->maxzeroextensions, TRUE, DEFAULT_MAXZEROEXTENSIONS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/clique/cliquetablemem",
         "maximal memory size of dense clique table (in kb)",
         &sepadata->cliquetablemem, TRUE, DEFAULT_CLIQUETABLEMEM, 0.0, (SCIP_Real)INT_MAX/1024.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/clique/cliquedensity",
         "minimal density of cliques to use a dense clique table",
         &sepadata->cliquedensity, TRUE, DEFAULT_CLIQUEDENSITY, 0.0, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}

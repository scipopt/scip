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
#pragma ident "@(#) $Id: sepa_clique.c,v 1.10 2005/08/08 11:16:35 bzfpfend Exp $"

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
   int              ncuts;              /**< number of cuts found */
   VAR**            vars;               /**< binary variables (contained in a 3-clique in implgraph) */    
   int              nvars;              /**< number of binary variables (contained in a 3-clique in implgraph)*/  
   Real*            varsolvals;         /**< LP solution of binary variables (contained in a 3-clique in implgraph) */
   Real             scaleval;           /**< factor for scaling weights */
   int              maxtreenodes;       /**< maximal number of nodes in branch and bound tree (-1: no limit) */
   Bool             tcliqueloaded;      /**< TRUE if tcliquegraph is allready loaded (tcliquegraph can be NULL),
                                         *   FALSE otherwise */ 
};





/*
 * Local methods
 */

/** gets index in tclique data structure of binary (active or negated) variable 
 *  with probindex of corresponding active binary variable */
static
int getVarTcliqueIndex(
   int              probindex,          /**< given probindex of cor. active binary variable */
   Bool             neg                 /**< TRUE if given variable is negated, FALSE if it is active */
   )
{
   /* index i of variable x in tclique data structure    
    *   i =  2 * probindex                 x is active  variable (neg = FALSE)
    *   i = (2 * probindex) + 1            x is negated varaible (neg = TRUE)
    */
   assert(FALSE == 0);
   assert(TRUE == 1);
   assert(neg == FALSE || neg == TRUE);

   return 2*probindex + (int)neg;
} 

/** gets probindex of active binary variable corresponding to given binary (active or negated) variable 
 *  with index i in tclique data structure and flag that states whether given variable is active or negated
 */ 
static
void getVarIndex(
   int              i,                  /**< given index of variable in tclique data structure */
   int*             probindex,          /**< pointer to store probindex of active variable cor. to given variable */
   Bool*            neg                 /**< pointer to store wh. given var is negated (TRUE if negated, FALSE if active) */
   )
{
   assert(neg != NULL);
   assert(probindex != NULL);

   /* probindex of corresponding active variable  
    *   probindex =  i   / 2               i is even (given variable is active)
    *   probindex = (i-1)/ 2               i is odd  (given variable is negated)
    */
   *neg = (i % 2 == 1);
   *probindex = i / 2;
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
   int* cliquegraphidx;
   Bool* implused;
   VAR** binvars;
   int nimplications;
   int nbinvars;
   int xi;

   assert(sepadata != NULL);
   assert(sepadata->tcliquegraph == NULL);
   assert(sepadata->vars == NULL);
   assert(sepadata->nvars == 0);

   /* get binary variables */
   binvars = SCIPgetVars(scip);
   nbinvars = SCIPgetNBinVars(scip);
   if( nbinvars == 0 )
      return SCIP_OKAY;

   /* get temporary memory for mapping variable indices to clique graph nodes and storing which implications are
    * already in the graph
    */
   nimplications = SCIPgetNImplications(scip);
   CHECK_OKAY( SCIPallocBufferArray(scip, &cliquegraphidx, 2*nbinvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &implused, nimplications) );
   for( xi = 0; xi < 2*nbinvars; ++xi )
      cliquegraphidx[xi] = -1;
   clearMemoryArray(implused, nimplications);
   sepadata->nvars = 0;

   debugMessage("initializing clique graph (%d implications)\n", nimplications);

   /* detect 3-cliques in the clique graph: triples (x,y,z) in the implication graph with x -> y', x -> z', and y -> z';
    * in order to avoid double checks, we only check triplets with variable indices xindex < yindex < zindex
    */
   for( xi = 0; xi < 2*(nbinvars-2); ++xi ) /* at least two variables must be left over for y and z */
   {
      VAR* x;
      int xindex;
      int xprobindex;
      Bool xneg;
      VAR** ximplvars;
      BOUNDTYPE* ximpltypes; 
      int* ximplids;
      int xnbinimpls;
      int i;

      getVarIndex(xi, &xprobindex, &xneg);
      assert(0 <= xprobindex && xprobindex < nbinvars);
      x = binvars[xprobindex];
      xindex = SCIPvarGetIndex(x);

      /* scan the implications of xi == 1 for potential y candidates */
      xnbinimpls = SCIPvarGetNBinImpls(x, !xneg);
      ximplvars = SCIPvarGetImplVars(x, !xneg);
      ximpltypes = SCIPvarGetImplTypes(x, !xneg);
      ximplids = SCIPvarGetImplIds(x, !xneg);
      
      /* ignore implicants with yindex <= xindex */
      for( i = 0; i < xnbinimpls && SCIPvarGetIndex(ximplvars[i]) <= xindex; ++i )
      {}

      /* loop over all y > x */
      for( ; i < xnbinimpls-1; ++i ) /* at least one variable must be left over for z */
      {
         VAR* y;
         int yindex;
         int yprobindex;
         Bool yneg;
         int yi;
         VAR** yimplvars;
         BOUNDTYPE* yimpltypes; 
         int* yimplids;
         int ynbinimpls;
         int xk;
         int yk;
         int xyid;

         y = ximplvars[i];
         yindex = SCIPvarGetIndex(y);
         assert(xindex < yindex); /* the implied variables are sorted by increasing variable index */
         xyid = ximplids[i];

         /* check, whether the implicant is y == 0 or y == 1 (yi is the variable that conflicts with xi == 1) */
         yneg = (ximpltypes[i] == SCIP_BOUNDTYPE_LOWER);
         yprobindex = SCIPvarGetProbindex(y);
         yi = getVarTcliqueIndex(yprobindex, yneg);

         /* scan the remaining implications of xi == 1 and the implications of yi == 1 for equal entries */
         ynbinimpls = SCIPvarGetNBinImpls(y, !yneg);
         yimplvars = SCIPvarGetImplVars(y, !yneg);
         yimpltypes = SCIPvarGetImplTypes(y, !yneg);
         yimplids = SCIPvarGetImplIds(y, !yneg);

         /* we simultaneously scan both implication arrays for candidates z with x < y < z */
         xk = i+1;
         assert(yindex < SCIPvarGetIndex(ximplvars[xk]));
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
               int zprobindex;
               Bool zneg;
               int zi;
               int xzid;
               int yzid;

               /* we found z with x < y < z and xi + yi + zi <= 1 */
               z = ximplvars[xk];
               xzid = ximplids[xk];
               yzid = yimplids[yk];
               assert(SCIPvarGetIndex(z) == zindex);

               /* check, whether the implicant is z == 0 or z == 1 (zi is the variable that conflicts with xi == 1) */
               zneg = (ximpltypes[xk] == SCIP_BOUNDTYPE_LOWER);
               zprobindex = SCIPvarGetProbindex(z);
               zi = getVarTcliqueIndex(zprobindex, zneg);

               /* create the tclique data structure, if not yet existing */
               if( sepadata->tcliquegraph == NULL )
               {
                  assert(sepadata->vars == NULL);
                  assert(sepadata->nvars == 0);

                  if( !tcliqueCreate(&sepadata->tcliquegraph) )
                  {
                     debugMessage("failed to create tclique data structure\n");
                     return SCIP_NOMEMORY;
                  }
                  CHECK_OKAY( SCIPallocMemoryArray(scip, &sepadata->vars, 2*nbinvars) );
               }
               assert(sepadata->tcliquegraph != NULL);
               assert(sepadata->vars != NULL);

               /* add nodes xi, yi, and zi to clique graph (if not yet existing) */
               if( cliquegraphidx[xi] == -1 )
               {
                  assert(sepadata->nvars < 2*nbinvars);
                  if( !tcliqueAddNode(sepadata->tcliquegraph, sepadata->nvars, 0) )
                  {
                     debugMessage("failed to add node %d to clique graph\n", sepadata->nvars);
                     return SCIP_NOMEMORY;
                  }
                  cliquegraphidx[xi] = sepadata->nvars;
                  if( xneg )
                  {
                     CHECK_OKAY( SCIPgetNegatedVar(scip, x, &sepadata->vars[sepadata->nvars]) ); 
                  }
                  else
                     sepadata->vars[sepadata->nvars] = x;
                  sepadata->nvars++;
               }
               if( cliquegraphidx[yi] == -1 )
               {
                  assert(sepadata->nvars < 2*nbinvars);
                  if( !tcliqueAddNode(sepadata->tcliquegraph, sepadata->nvars, 0) )
                  {
                     debugMessage("failed to add node %d to clique graph\n", sepadata->nvars);
                     return SCIP_NOMEMORY;
                  }
                  cliquegraphidx[yi] = sepadata->nvars;
                  if( yneg )
                  {
                     CHECK_OKAY( SCIPgetNegatedVar(scip, y, &sepadata->vars[sepadata->nvars]) ); 
                  }
                  else
                     sepadata->vars[sepadata->nvars] = y;
                  sepadata->nvars++;
               }
               if( cliquegraphidx[zi] == -1 )
               {
                  assert(sepadata->nvars < 2*nbinvars);
                  if( !tcliqueAddNode(sepadata->tcliquegraph, sepadata->nvars, 0) )
                  {
                     debugMessage("failed to add node %d to clique graph\n", sepadata->nvars);
                     return SCIP_NOMEMORY;
                  }
                  cliquegraphidx[zi] = sepadata->nvars;
                  if( zneg )
                  {
                     CHECK_OKAY( SCIPgetNegatedVar(scip, z, &sepadata->vars[sepadata->nvars]) ); 
                  }
                  else
                     sepadata->vars[sepadata->nvars] = z;
                  sepadata->nvars++;
               }

               /* add edges (xi,yi), (xi,zi), and (yi,zi) to clique graph (multiple edges are deleted afterwards) */
               if( !implused[xyid] )
               {
                  if( !tcliqueAddEdge(sepadata->tcliquegraph, cliquegraphidx[xi], cliquegraphidx[yi]) )
                  {
                     debugMessage("failed to add edge (%d,%d) to clique graph\n",
                        cliquegraphidx[xi], cliquegraphidx[yi]);
                     return SCIP_NOMEMORY;
                  }
                  implused[xyid] = TRUE;
               }
               if( !implused[xzid] )
               {
                  if( !tcliqueAddEdge(sepadata->tcliquegraph, cliquegraphidx[xi], cliquegraphidx[zi]) )
                  {
                     debugMessage("failed to add edge (%d,%d) to clique graph\n",
                        cliquegraphidx[xi], cliquegraphidx[zi]);
                     return SCIP_NOMEMORY;
                  }
                  implused[xzid] = TRUE;
               }
               if( !implused[yzid] )
               {
                  if( !tcliqueAddEdge(sepadata->tcliquegraph, cliquegraphidx[yi], cliquegraphidx[zi]) )
                  {
                     debugMessage("failed to add edge (%d,%d) to clique graph\n",
                        cliquegraphidx[yi], cliquegraphidx[zi]);
                     return SCIP_NOMEMORY;
                  }
                  implused[yzid] = TRUE;
               }
            }

            /* proceed with the next pair of implications */
            xk++;
            yk++;
         }
      }
   }
   assert((sepadata->nvars == 0) == (sepadata->vars == NULL));
   assert((sepadata->nvars == 0) == (sepadata->tcliquegraph == NULL));

   /* reduce the size of the vars array */
   if( 0 < sepadata->nvars && sepadata->nvars < 2*nbinvars )
   {
      CHECK_OKAY( SCIPreallocMemoryArray(scip, &sepadata->vars, sepadata->nvars) );
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &implused);
   SCIPfreeBufferArray(scip, &cliquegraphidx);

   /* flush the changes to the clique graph */
   if( sepadata->tcliquegraph != NULL )
   {
      debugMessage(" -> flushing clique graph\n");
      if( !tcliqueFlush(sepadata->tcliquegraph) )
         return SCIP_NOMEMORY;
      
      debugMessage(" -> clique graph constructed (%d nodes, %d edges)\n", 
         tcliqueGetNNodes(sepadata->tcliquegraph), tcliqueGetNEdges(sepadata->tcliquegraph));

#if 0
      debugPrintf("---v---- tcliquegraph ----v---\n");
      debug(tcliquePrintData(sepadata->tcliquegraph));
      debugPrintf("---^---- tcliquegraph ----^---\n");
#endif
   }

   return SCIP_OKAY;
}

/* updates tclique data structure (changes weight of all nodes) */
static 
void updateTcliquegraph(
   SCIP*            scip,               /**< SCIP data structure */
   SEPADATA*        sepadata            /**< separator data */
   )
{
   int i;

   assert(sepadata != NULL);
   assert(sepadata->tcliquegraph != NULL);
   assert(tcliqueGetNNodes(sepadata->tcliquegraph) == sepadata->nvars);

   /* updates weight of all nodes in tclique data structure */
   for( i = 0; i < sepadata->nvars; i++ )
   {
      int weight;

      weight = (int)SCIPfeasFloor(scip, sepadata->varsolvals[i] * sepadata->scaleval);
      weight = MAX(weight, 0);
      tcliqueChangeWeight(sepadata->tcliquegraph, i, weight);
   }
}

/** generates cuts using a clique found by algorithm for maximum weight clique 
 *  and decides whether to stop generating cliques with the algorithm for maximum weight clique
 */
static
TCLIQUE_NEWSOL(tcliqueNewClique)
{
   SEPADATA* sepadata;
   VAR** vars;
   Real* varsolvals;
   Real unscaledweight;
   WEIGHT minweightinc;
   int i;
   int nvars;

   assert(acceptsol != NULL);
   assert(stopsolving != NULL);

   sepadata = (SEPADATA*)tcliquedata;

   assert(sepadata != NULL);
   assert(sepadata->scip != NULL);
   assert(sepadata->ncuts >= 0);
   assert(sepadata->nvars > 0);
   assert(sepadata->vars != NULL);
   assert(sepadata->varsolvals != NULL);
   
   /* we don't accept the solution as new incumbent, because we want to find many violated clique inequalities */
   *acceptsol = FALSE;
   *stopsolving = FALSE;

   /* slightly increase the minimal weight for additional cliques */
   minweightinc = (cliqueweight - *minweight)/MIN(1, (int)(sepadata->scaleval/10.0));
   minweightinc = MIN(minweightinc, 1);
   *minweight += minweightinc;

   vars = sepadata->vars; 
   nvars = sepadata->nvars; 
   varsolvals = sepadata->varsolvals; 

   unscaledweight = 0.0;

   /* adds cut if weight of the clique is greater than 1 */
   if( cliqueweight > sepadata->scaleval )
   {
      for( i = 0; i < ncliquenodes; i++ )
         unscaledweight += varsolvals[cliquenodes[i]];
      
      /*debugMessage("->%g, ncuts=%d\n", unscaledweight, sepadata->ncuts);*/

      if( SCIPisEfficacious(sepadata->scip, unscaledweight - 1.0) )
      {
         ROW* cut;
         char cutname[MAXSTRLEN];
      
         /* create the cut */
         sprintf(cutname, "clique%d_%d", SCIPgetNLPs(sepadata->scip), sepadata->ncuts);
         CHECK_ABORT( SCIPcreateEmptyRow(sepadata->scip, &cut, cutname, -SCIPinfinity(sepadata->scip), 1.0, 
               FALSE, FALSE, TRUE) );

         CHECK_ABORT( SCIPcacheRowExtensions(sepadata->scip, cut) );
         assert(ncliquenodes <= nvars);
         debugMessage("clique in graph:");
         for( i = 0; i < ncliquenodes; ++i )
         {
            assert(cliquenodes[i] < nvars);
            CHECK_ABORT( SCIPaddVarToRow(sepadata->scip, cut, vars[cliquenodes[i]], 1.0) );
            debugPrintf(" [%d]<%s>", cliquenodes[i], SCIPvarGetName(vars[cliquenodes[i]]));
         }
         debugPrintf("\n");
         CHECK_ABORT( SCIPflushRowExtensions(sepadata->scip, cut) );

         debugMessage("found clique cut (act=%g): ", unscaledweight);
         debug(SCIPprintRow(sepadata->scip, cut, NULL));

         CHECK_ABORT( SCIPaddCut(sepadata->scip, cut, FALSE) );
         (sepadata->ncuts)++;
         
         /* release the row */
         CHECK_ABORT( SCIPreleaseRow(sepadata->scip, &cut) );
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
   assert((sepadata->tcliquegraph == NULL && sepadata->nvars == 0 && sepadata->vars == NULL) || 
      (sepadata->tcliquegraph != NULL && sepadata->nvars > 0 && sepadata->vars != NULL) );

   /* free tclique data */
   if( sepadata->tcliquegraph != NULL )
   {
      SCIPfreeMemoryArray(scip, &(sepadata->vars));
      tcliqueFree(&(sepadata->tcliquegraph));
      sepadata->nvars = 0;
   }
   sepadata->tcliqueloaded = FALSE;

   return SCIP_OKAY;
}

/** execution method of separator */
static
DECL_SEPAEXEC(sepaExecClique)
{  /*lint --e{715}*/
   SEPADATA* sepadata;
   int* cliquenodes; 	        
   WEIGHT cliqueweight;    
   int ncliquenodes;	        
   int maxtreenodes;

   assert(scip != NULL);
   
   /* get separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   
   assert((sepadata->tcliquegraph == NULL && sepadata->nvars == 0 && sepadata->vars == NULL)
      || (sepadata->tcliquegraph != NULL && sepadata->nvars > 0 && sepadata->vars != NULL));

   *result = SCIP_DIDNOTRUN;
   sepadata->ncuts = 0;

   /* if we already detected that no implications between binary variables exist, nothing has to be done */
   if( sepadata->tcliquegraph == NULL && sepadata->tcliqueloaded )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;
      
   /* load tclique data structure */
   if( !sepadata->tcliqueloaded )
   {
      assert(sepadata->tcliquegraph == NULL);

      CHECK_OKAY( loadTcliquegraph(scip, sepadata) );
      assert((sepadata->tcliquegraph == NULL) == (sepadata->vars == NULL));
      sepadata->tcliqueloaded = TRUE;
   
      if( sepadata->vars == NULL )
      {
         /* we did not find any variables that are contained in a clique with at least 3 variables in the
          * implication graph -> we disable the clique cut separator
          */
         assert(sepadata->nvars == 0);
         return SCIP_OKAY;
      }
   }
   
   /* updates LP-solution in sepadata and weights in tclique data structure */
   assert(sepadata->tcliquegraph != NULL);
   CHECK_OKAY( SCIPallocBufferArray(scip, &sepadata->varsolvals, sepadata->nvars) );
   CHECK_OKAY( SCIPgetVarSols(scip, sepadata->nvars, sepadata->vars, sepadata->varsolvals) );
   updateTcliquegraph(scip, sepadata);

   /* get maximal number of tree nodes */
   maxtreenodes = (sepadata->maxtreenodes == -1 ? INT_MAX : sepadata->maxtreenodes);

   /* finds maximum weight clique in tclique */
   CHECK_OKAY( SCIPallocBufferArray(scip, &cliquenodes, tcliqueGetNNodes(sepadata->tcliquegraph)) );
   tcliqueMaxClique(NULL, NULL, NULL, NULL, sepadata->tcliquegraph, tcliqueNewClique, (TCLIQUEDATA*)sepadata,
      cliquenodes, &ncliquenodes, &cliqueweight, (int)sepadata->scaleval-1, (int)sepadata->scaleval+1, maxtreenodes);

   /* frees data structures */
   SCIPfreeBufferArray(scip, &cliquenodes);
   SCIPfreeBufferArray(scip, &sepadata->varsolvals);

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
   sepadata->ncuts = 0;
   sepadata->vars = NULL;  
   sepadata->nvars=0;   
   sepadata->varsolvals = NULL; 
   sepadata->tcliqueloaded = FALSE;

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

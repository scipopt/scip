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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: sepa_clique.c,v 1.1 2005/03/10 17:11:15 bzfpfend Exp $"

/**@file   sepa_clique.c
 * @brief  clique separator
 * @author Tobias Achterberg
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "sepa_clique.h"
#include "tclique/graph.h"
#include "tclique/branch.h"


/* separator properties */
#define SEPA_NAME              "clique"
#define SEPA_DESC              "clique separator of stable set relaxation"
#define SEPA_PRIORITY             -5000
#define SEPA_FREQ                    10
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_SCALEVAL         1000.0 /**< factor for scaling weights */
#define DEFAULT_MAXTREENODES        100 /**< maximal number of nodes in branch and bound tree (-1: no limit) */




/*
 * Data structures
 */

/** separator data */
struct SepaData
{
   TCLIQUEDATA*     tcliquedata;        /**< tclique data structure */
   SCIP*            scip;               /**< SCIP data structure */
   int              ncuts;              /**< number of cuts found */
   VAR**            vars;               /**< binary variables (contained in a 3-clique in implgraph) */    
   int              nvars;              /**< number of binary variables (contained in a 3-clique in implgraph)*/  
   Real*            varsolvals;         /**< LP solution of binary variables (contained in a 3-clique in implgraph) */
   Real             scaleval;           /**< factor for scaling weights */
   int              maxtreenodes;       /**< maximal number of nodes in branch and bound tree (-1: no limit) */
   Bool             tcliqueloaded;      /**< TRUE if tcliquedata is allready loaded (tcliquedata can be NULL), FALSE otherwise */ 
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
   int tcliqueindex; 

   /* index i of variable x in tclique data structure    
    *   i =  2 * probindex                 x is active  variable (neg = FALSE)
    *   i = (2 * probindex) + 1            x is negated varaible (neg = TRUE) */
   if( !neg )
      tcliqueindex = 2 * probindex;
   else
      tcliqueindex = (2 * probindex) + 1;

   return tcliqueindex;
} 

/** gets probindex of active binary variable corresponding to given binary (active or negated) variable 
 *  with index i in tclique data structure and flag that states whether given variable is active or negated */ 
static
void getVarIndex(
   int              i,                  /**< given index of variable in tclique data structure */
   Bool*            neg,                /**< pointer to store wh. given var is negated (TRUE if negated, FALSE if active) */
   int*             probindex           /**< pointer to store probindex of active variable cor. to given variable */
   )
{
   assert(neg != NULL);
   assert(probindex != NULL);

   /* probindex of corresponding active variable  
    *   probindex =  i   / 2               i is eaven (given variable is active)
    *   probindex = (i-1)/ 2               i is odd   (given variable is negated) */
   if( i % 2 == 0 )
   {
      *neg = FALSE;
      *probindex = i / 2;
   }
   else
   {
      *neg = TRUE;
      *probindex = i / 2;
   }
} 

/* adds given variable contained in a clique with 3 variables in implication graph 
 * to the set of variables used for tclique data structure. */
static 
RETCODE addVar(
   SCIP*            scip,               /**< SCIP data structure */
   SEPADATA*        sepadata,           /**< separator data */
   VAR**            binvars,            /**< active binary problem variables */ 
   int*             inthreeclique,      /**< index of all binary variables (active or negated) in set of vars for tclique */
   int              i,                  /**< index of variable to add in inthreeclique array */
   int*             nvars,
   VAR**            vars
   )
{
   VAR* var;
   Bool neg; 
   int probindex;
   
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(i >= 0);
   assert(nvars != NULL);
   assert(*nvars >= 0);
   assert(vars != NULL);
   assert(inthreeclique[i] >= -1);
   assert(inthreeclique[i] < *nvars);

   /* variable is not added yet */
   if( inthreeclique[i] == -1 )
   {
      /* sets index of variable in set of vars for tclique */
      inthreeclique[i] = *nvars;

      /* adds variable to set of vars for tclique */
      getVarIndex(i, &neg, &probindex);
      if( neg )
      {
         CHECK_OKAY( SCIPgetNegatedVar(scip, binvars[probindex], &var) ); 
         vars[inthreeclique[i]] = var;
      }
      else
      {
         assert(!neg);
         vars[inthreeclique[i]] = binvars[probindex];
      }
      (*nvars)++;
   } 
   /* only asserts */
   else
   {
      getVarIndex(i, &neg, &probindex);
      if( neg )
      {
         CHECK_OKAY( SCIPgetNegatedVar(scip, binvars[probindex], &var) ); 
         assert(vars[inthreeclique[i]] == var);
      }
      else
      {
         assert(vars[inthreeclique[i]] == binvars[probindex]);
      }
   }

   return SCIP_OKAY;
}

/* loads tclique data structure using implication graph; all nodes have dummy weight 1s 
 * and are contained in a clique with three nodes in implication graph */
static 
RETCODE loadTcliquedata(
   SCIP*            scip,               /**< SCIP data structure */
   SEPADATA*        sepadata,           /**< separator data */
   VAR**            binvars,            /**< active binary problem variables */ 
   int              nbinvars,           /**< number of binary active problem variables */
   VAR**            vars,               /**< pointer to store all binary variables contained in a 3-clique in implgraph */
   int*             nvars               /**< pointer to store number of binary vars contained in a 3-clique in implgraph */
   )
{
   int i;
   Bool* marked;
   int* inthreeclique;
   
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(sepadata->tcliquedata == NULL );
   assert(sepadata->vars == NULL);
   assert(sepadata->nvars == 0);

   /* gets data structures */
   CHECK_OKAY( SCIPallocBufferArray(scip, &marked, 2 * nbinvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &inthreeclique, 2 * nbinvars) );

   *nvars = 0;

   /* initialises inthreeclique */
   for( i = 0; i < 2 * nbinvars; i++ )
      inthreeclique[i] = -1;

   /* gets all binary variables (active and negated) that are contained in a clique 
    * with 3 variables in implication graph */
   for( i = 0; i < 2 * nbinvars; i++ )
   {
      int j;
      Bool neg;
      int probindex;
      BOUNDTYPE* impltypes; 
      Real* implbounds; 
      VAR** implvars;
      int nbinimpl;
      VAR* maxmarked;

      assert(inthreeclique[i] >= -1);
      assert(inthreeclique[i] < *nvars);

      if( inthreeclique[i] > -1) 
         continue;

      /* unmarks all variables */
      clearMemoryArray(marked, 2 *nbinvars);

      /* gets corresponding active variable x for current variable x_i */
      getVarIndex(i, &neg, &probindex);
      
      /* gets all implications x >= 1 ==> y <= 0 or y >= 1, if current variable x_i is active */
      if( !neg )
      {
         nbinimpl = SCIPvarGetNBinImpls(binvars[probindex], 1);
         implvars = SCIPvarGetImplVars(binvars[probindex], 1);
         impltypes = SCIPvarGetImplTypes(binvars[probindex], 1);
         implbounds = SCIPvarGetImplBounds(binvars[probindex], 1);
      }   
      /* gets all implications x <= 0 ==> y <= 0 or y >= 1, if current variable x_i is negated */
      else
      {
         assert(neg);

         nbinimpl = SCIPvarGetNBinImpls(binvars[probindex], 0);
         implvars = SCIPvarGetImplVars(binvars[probindex], 0);
         impltypes = SCIPvarGetImplTypes(binvars[probindex], 0);
         implbounds = SCIPvarGetImplBounds(binvars[probindex], 0);
      }  
      
      if( nbinimpl <= 1 )
      {
         assert(inthreeclique[i] == -1); 
         continue;
      }

      /* marks variables that are adjacent to current variable x_i in implication graph */
      for( j = 0; j < nbinimpl; j++ )
      {
         assert((impltypes[j] == SCIP_BOUNDTYPE_UPPER && SCIPisEQ(scip, implbounds[j], 0.0)) || 
            (impltypes[j] == SCIP_BOUNDTYPE_LOWER && SCIPisEQ(scip, implbounds[j], 1.0)));

         /* marks y for implication x <= 0/>=1 ==> y <= 0 */
         if( impltypes[j] == SCIP_BOUNDTYPE_UPPER )
            marked[getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), FALSE)] = TRUE;
         /* marks y' for implication x <= 0/>=1 ==> y >= 1 */
         else
         {
            assert(impltypes[j] == SCIP_BOUNDTYPE_LOWER);
            marked[getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), TRUE)] = TRUE;
         }
      }
      maxmarked = implvars[nbinimpl-1];
      
      /* searches for two adjacent marked variables (y,z) --> clique = {x_i, y, z} */ 
      for( j = 0; j < nbinimpl && inthreeclique[i] == -1; j++ )
      {
         int k;
         BOUNDTYPE* impltypesmarked; 
         Real* implboundsmarked; 
         VAR** implvarsmarked;
         int nbinimplmarked;
         
         /* gets all implications x_j >= 1 ==> y <= 0 or y >= 1, if x_j is active (x_i <=0/>=1 ==> x_j <= 0) */
         if( impltypes[j] == SCIP_BOUNDTYPE_UPPER )
         {
            if( !marked[getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), FALSE)] )
               continue;
            
            nbinimplmarked = SCIPvarGetNBinImpls(implvars[j], 1);
            implvarsmarked = SCIPvarGetImplVars(implvars[j], 1);
            impltypesmarked = SCIPvarGetImplTypes(implvars[j], 1);
            implboundsmarked = SCIPvarGetImplBounds(implvars[j], 1);
         }   
         /* gets all implications x_j <= 0 ==> y <= 0 or y >= 1, if x_j is negated (x_i <=0/>=1 ==> x_j >= 1) */
         else
         {
            assert(impltypes[j] == SCIP_BOUNDTYPE_LOWER);

            if( !marked[getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), TRUE)] )
               continue;

            nbinimplmarked = SCIPvarGetNBinImpls(implvars[j], 0);
            implvarsmarked = SCIPvarGetImplVars(implvars[j], 0);
            impltypesmarked = SCIPvarGetImplTypes(implvars[j], 0);
            implboundsmarked = SCIPvarGetImplBounds(implvars[j], 0);
         }  
         
         /* searches for marked variables in implications for x_j */
         for( k = 0; k < nbinimplmarked && inthreeclique[i] == -1 && implvarsmarked[k] <= maxmarked; k++ )
         {
            /* implication x_j <= 0 or x_j >= 1 ==> x_k <= 0 && x_k is marked */
            if( impltypesmarked[k] == SCIP_BOUNDTYPE_UPPER 
               && marked[getVarTcliqueIndex(SCIPvarGetProbindex(implvarsmarked[k]), FALSE)])
            {
               addVar(scip, sepadata, binvars, inthreeclique, i, nvars, vars);

               if( impltypes[j] == SCIP_BOUNDTYPE_UPPER )
                  addVar(scip, sepadata, binvars, inthreeclique, 
                     getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), FALSE), nvars, vars);
               else
               {
                  assert(impltypes[j] == SCIP_BOUNDTYPE_LOWER);
                  addVar(scip, sepadata, binvars, inthreeclique, 
                     getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), TRUE), nvars, vars);
               }
               addVar(scip, sepadata, binvars, inthreeclique, 
                  getVarTcliqueIndex(SCIPvarGetProbindex(implvarsmarked[k]), FALSE), nvars, vars);
               break;
            }
            /* implication x_j <= 0 or x_j >= 1 ==> x_k >= 1 && x_k' is marked */
            if( impltypesmarked[k] == SCIP_BOUNDTYPE_LOWER 
               && marked[getVarTcliqueIndex(SCIPvarGetProbindex(implvarsmarked[k]), TRUE)])
            {
               addVar(scip, sepadata, binvars, inthreeclique, i, nvars, vars);

               if( impltypes[j] == SCIP_BOUNDTYPE_UPPER )
                  addVar(scip, sepadata, binvars, inthreeclique, 
                     getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), FALSE), nvars, vars);
               else
               {
                  assert(impltypes[j] == SCIP_BOUNDTYPE_LOWER);
                  addVar(scip, sepadata, binvars, inthreeclique, 
                     getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), TRUE), nvars, vars);
               }
               addVar(scip, sepadata, binvars, inthreeclique, 
                  getVarTcliqueIndex(SCIPvarGetProbindex(implvarsmarked[k]), TRUE), nvars, vars);
               break;
            }
         }
      }
   }

   if( *nvars > 0 )
   {
      /* adds a node with weight = 1 for all binary variables (active or negated) that are contained 
       * in a clique with 3 variables in implication graph  */
      for( i = 0; i < *nvars; i++ )
         tcliqueAddNode(&sepadata->tcliquedata, i, 0);

      /* adds the following edges for all binary variables x (active or negated) that are contained 
       * in a clique with 3 variables in implication graph:  
       *    (x , y )      for all implications x >= 1 ==> y <= 0
       *    (x , y')      for all implications x >= 1 ==> y >= 1,             y' = 1 - x
       *    (x', y )      for all implications x <= 0 ==> y <= 0, x' = 1 - x 
       *    (x', y')      for all implications x <= 0 ==> y >= 1, x' = 1 - x, y' = 1 - y */
      for( i = 0; i < nbinvars; i++ )
      {
         BOUNDTYPE* impltypes; 
         Real* implbounds; 
         VAR** implvars;
         int j;
         int nbinimpl;
      
         assert(inthreeclique[getVarTcliqueIndex(i, FALSE)] >= -1);

         /* active binary variable */
         if( inthreeclique[getVarTcliqueIndex(i, FALSE)] > -1)
         {
            assert(inthreeclique[getVarTcliqueIndex(i, FALSE)] < *nvars);

            /* gets all implications x >= 1 ==> y <= 0 or y >= 1 */
            nbinimpl = SCIPvarGetNBinImpls(binvars[i], 1);
            implvars = SCIPvarGetImplVars(binvars[i], 1);
            impltypes = SCIPvarGetImplTypes(binvars[i], 1);
            implbounds = SCIPvarGetImplBounds(binvars[i], 1);
         
            /* adds edges for all implication x >= 1 ==> y <= 0 / y >= 1 to tclique data structure */  
            for( j = 0; j < nbinimpl; j++ )
            {
               assert((impltypes[j] == SCIP_BOUNDTYPE_UPPER && SCIPisEQ(scip, implbounds[j], 0.0)) || 
                  (impltypes[j] == SCIP_BOUNDTYPE_LOWER && SCIPisEQ(scip, implbounds[j], 1.0)));
            
               if( impltypes[j] == SCIP_BOUNDTYPE_UPPER 
                  && inthreeclique[getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), FALSE)] > -1 )
               {
                  assert(inthreeclique[getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), FALSE)] < *nvars);
         
                  tcliqueAddEdge(&sepadata->tcliquedata, inthreeclique[getVarTcliqueIndex(i, FALSE)], 
                     inthreeclique[getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), FALSE)]);
               }
               if( impltypes[j] == SCIP_BOUNDTYPE_LOWER 
                  && inthreeclique[getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), TRUE)] > -1 )
               {
                  assert(inthreeclique[getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), TRUE)] < *nvars);

                  tcliqueAddEdge(&sepadata->tcliquedata, inthreeclique[getVarTcliqueIndex(i, FALSE)], 
                     inthreeclique[getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), TRUE)]);
               }
            }  
         }
      
         assert(inthreeclique[getVarTcliqueIndex(i, TRUE)] >= -1);

         /* negated binary variable */
         if( inthreeclique[getVarTcliqueIndex(i, TRUE)] > -1)
         {
            assert(inthreeclique[getVarTcliqueIndex(i, TRUE)] < *nvars);

            /* gets all implications x <= 0 ==> y <= 0 or y >= 1 */
            nbinimpl = SCIPvarGetNBinImpls(binvars[i], 0);
            implvars = SCIPvarGetImplVars(binvars[i], 0);
            impltypes = SCIPvarGetImplTypes(binvars[i], 0);
            implbounds = SCIPvarGetImplBounds(binvars[i], 0);
         
            /* adds edges for all implication x <= 0 ==> y <= 0 / y >= 1 to tclique data structure */  
            for( j = 0; j < nbinimpl; j++ )
            {
               assert((impltypes[j] == SCIP_BOUNDTYPE_UPPER && SCIPisEQ(scip, implbounds[j], 0.0)) || 
                  (impltypes[j] == SCIP_BOUNDTYPE_LOWER && SCIPisEQ(scip, implbounds[j], 1.0)));
            
               if( impltypes[j] == SCIP_BOUNDTYPE_UPPER 
                  && inthreeclique[getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), FALSE)] > -1)
               {
                  assert(inthreeclique[getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), FALSE)] < *nvars);
         
                  tcliqueAddEdge(&sepadata->tcliquedata, inthreeclique[getVarTcliqueIndex(i, TRUE)], 
                     inthreeclique[getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), FALSE)]);
               }
               if( impltypes[j] == SCIP_BOUNDTYPE_LOWER 
                  && inthreeclique[getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), TRUE)] > -1)
               {
                  assert(inthreeclique[getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), TRUE)] < *nvars);

                  tcliqueAddEdge(&sepadata->tcliquedata, inthreeclique[getVarTcliqueIndex(i, TRUE)], 
                     inthreeclique[getVarTcliqueIndex(SCIPvarGetProbindex(implvars[j]), TRUE)]);
               }
            }  
         }
      }
   }

#ifdef DEBUG
   if( sepadata->tcliquedata != NULL )
      printTcliquedata(sepadata->tcliquedata);
   else
      debugMessage("no tclique data structure created, because no variable is cont. in clique with 3 nodes!\n");
#endif

   /* frees data structures */
   SCIPfreeBufferArray(scip, &inthreeclique);
   SCIPfreeBufferArray(scip, &marked);

   return SCIP_OKAY;
}

/* updates tclique data structure (changes weight of all nodes) */
static 
void updateTcliquedata(
   SCIP*            scip,               /**< SCIP data structure */
   SEPADATA*        sepadata            /**< separator data */
   )
{
   int i;

   assert(sepadata != NULL);
   assert(sepadata->tcliquedata != NULL);
   assert(getNnodes(sepadata->tcliquedata) == sepadata->nvars);

   /* updates weight of all nodes in tclique data structure */
   for( i = 0; i < sepadata->nvars; i++ )
   {
      tcliqueNodeChangeWeight(&sepadata->tcliquedata, i, 
         (int)SCIPfeasFloor(scip, sepadata->varsolvals[i] * sepadata->scaleval));
   }
}

/** generates cuts using a clique found by algorithm for maximum weight clique 
 *  and decides whether to stop generating cliques with the algorithm for maximum weight clique */
static
TCLIQUE_USRCALLBACK(tcliqueNewClique)
{
   SEPADATA* sepadata;
   VAR** vars;
   Real* varsolvals;
   Real weightmwcunscaled;
   int i;
   int nvars;

   sepadata = (SEPADATA*)usrdata;

   assert(sepadata != NULL);
   assert(sepadata->scip != NULL);
   assert(sepadata->ncuts >= 0);
   assert(sepadata->nvars > 0);
   assert(sepadata->vars != NULL);
   assert(sepadata->varsolvals != NULL);
   
   vars = sepadata->vars; 
   nvars = sepadata->nvars; 
   varsolvals = sepadata->varsolvals; 

#ifdef DEBUG
   printf("found %der mwc=[ ", nmwc);
   for( i = 0; i < nmwc; i++ )
   {
      printf("%s(%g), ", SCIPvarGetName(vars[mwc[i]]), varsolvals[mwc[i]]);
   }
   printf("] wmwc=%d", weightmwc);
#endif

   weightmwcunscaled = 0.0;

   /* adds cut if weight of the clique is greater than 1 */
   if( weightmwc > sepadata->scaleval )
   {
      for( i = 0; i < nmwc; i++ )
         weightmwcunscaled += varsolvals[mwc[i]];
      
#ifdef DEBUG
      printf("->%g, ncuts=%d\n", weightmwcunscaled, sepadata->ncuts);
#endif
      if( SCIPisEfficacious(sepadata->scip, weightmwcunscaled - 1.0) )
      {
         ROW* cut;
         char cutname[MAXSTRLEN];
      
         /* creates the cut */
         sprintf(cutname, "clique%d_%d", SCIPgetNLPs(sepadata->scip), sepadata->ncuts);
         CHECK_OKAY( SCIPcreateEmptyRow(sepadata->scip, &cut, cutname, -SCIPinfinity(sepadata->scip), 1.0, 
               FALSE, FALSE, TRUE) );

         CHECK_OKAY( SCIPcacheRowExtensions(sepadata->scip, cut) );
         assert(nmwc <= nvars);
         for( i = 0; i < nmwc; ++i )
         {
            assert(mwc[i] < nvars);
            CHECK_OKAY( SCIPaddVarToRow(sepadata->scip, cut, vars[mwc[i]], 1.0) );
         }
         CHECK_OKAY( SCIPflushRowExtensions(sepadata->scip, cut) );

#ifdef DEBUG
         printf("-----> CUT FOUND: ");
         SCIPprintRow(sepadata->scip, cut, NULL);
#endif
         CHECK_OKAY( SCIPaddCut(sepadata->scip, cut, FALSE) );
         (sepadata->ncuts)++;
         
         /* release the row */
         CHECK_OKAY( SCIPreleaseRow(sepadata->scip, &cut) );
      }
   }

   /* decides whether to stop the search for a maximum weight clique */
   return FALSE;
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
#if 0
static
DECL_SEPAINIT(sepaInitClique)
{  /*lint --e{715}*/
   errorMessage("method of clique separator not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitClique NULL
#endif


/** deinitialization method of separator (called before transformed problem is freed) */
#if 0
static
DECL_SEPAEXIT(sepaExitClique)
{  /*lint --e{715}*/
   errorMessage("method of clique separator not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitClique NULL
#endif

/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#if 0
static
DECL_SEPAINITSOL(sepaInitsolClique)
{  /*lint --e{715}*/
   errorMessage("method of clique separator not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitsolClique NULL
#endif

/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
DECL_SEPAEXITSOL(sepaExitsolClique)
{ 
   SEPADATA* sepadata;
   
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert((sepadata->tcliquedata == NULL && sepadata->nvars == 0 && sepadata->vars == NULL) || 
      (sepadata->tcliquedata != NULL && sepadata->nvars > 0 && sepadata->vars != NULL) );

   /* free tclique data */
   if( sepadata->tcliquedata != NULL )
   {
      SCIPfreeMemoryArray(scip, &(sepadata->vars));
      tcliqueFree(&(sepadata->tcliquedata));
      sepadata->nvars = 0;
   }
   sepadata->tcliqueloaded = FALSE;

   return SCIP_OKAY;
}

/** execution method of separator */
static
DECL_SEPAEXEC(sepaExecClique)
{
   VAR** binvars;
   SEPADATA* sepadata;
   Real* varsolvals;
   int* mwc; 	        
   WEIGHT weightmwc;    
   int nmwc;	        
   int nbinvars;
   int maxtreenodes;

   assert(scip != NULL);
   
   /* gets separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   
   assert((sepadata->tcliquedata == NULL && sepadata->nvars == 0 && sepadata->vars == NULL)
      || (sepadata->tcliquedata != NULL && sepadata->nvars > 0 && sepadata->vars != NULL));

   *result = SCIP_DIDNOTRUN;
   sepadata->ncuts = 0;
   
   /* gets active problem variables */
   binvars = SCIPgetVars(scip);
   nbinvars = SCIPgetNBinVars(scip);
   assert(binvars != NULL);

   if( nbinvars == 0 )
      return SCIP_OKAY;

   if( sepadata->tcliquedata == NULL && sepadata->tcliqueloaded )
      return SCIP_OKAY;

   assert((!sepadata->tcliqueloaded && sepadata->tcliquedata == NULL)
      || (sepadata->tcliqueloaded && sepadata->tcliquedata != NULL));

   *result = SCIP_DIDNOTFIND;
      
   /* loads tclique data structure */
   if( !sepadata->tcliqueloaded )
   {
      VAR** vars;
      int nvars;
      
      assert(sepadata->tcliquedata == NULL);

      CHECK_OKAY( SCIPallocBufferArray(scip, &vars, 2 * nbinvars) );
      nvars = 0;
      
      CHECK_OKAY( loadTcliquedata(scip, sepadata, binvars, nbinvars, vars, &nvars ) );
      sepadata->tcliqueloaded = TRUE;

      assert((sepadata->tcliquedata == NULL && nvars == 0) || 
         (sepadata->tcliquedata != NULL && nvars > 0) );
   
      /* found variables that are contained in a clique with 3 variables in implication graph 
       * -> creats tclique data structure */
      if( sepadata->tcliquedata != NULL )
      {
         int i;
         
         /* sets sepadata */
         CHECK_OKAY( SCIPallocMemoryArray(scip, &(sepadata->vars), nvars) );

         for( i = 0; i < nvars; i++ )
            sepadata->vars[i] = vars[i];

         sepadata->nvars = nvars;

         SCIPfreeBufferArray(scip, &vars);
      }
      /* did not find any variables that are contained in a clique with 3 variables in implication graph 
       * -> does not create tclique data structure */
      else
      {
         SCIPfreeBufferArray(scip, &vars);

         /* sets sepadata */
         assert(sepadata->vars == NULL);
         assert(sepadata->nvars == 0);

         return SCIP_OKAY;
      }
   }
   
   /* updates LP-solution in sepadata and weights in tclique data structure */
   assert(sepadata->tcliquedata != NULL);
   CHECK_OKAY( SCIPallocBufferArray(scip, &varsolvals, sepadata->nvars) );
   CHECK_OKAY( SCIPgetVarSols(scip, sepadata->nvars, sepadata->vars, varsolvals) );
   sepadata->varsolvals = varsolvals;
   updateTcliquedata(scip, sepadata);
   
   /* get maximal number of tree nodes */
   maxtreenodes = (sepadata->maxtreenodes == -1 ? INT_MAX : sepadata->maxtreenodes);

   /* finds maximum weight clique in tclique */
   CHECK_OKAY( SCIPallocBufferArray(scip, &mwc, getNnodes(sepadata->tcliquedata)) );
   maxClique(sepadata->tcliquedata, tcliqueNewClique, (void*)sepadata, mwc, &nmwc, &weightmwc, 
      (int)sepadata->scaleval-1, maxtreenodes);

   /* frees data structures */
   SCIPfreeBufferArray(scip, &mwc);
   SCIPfreeBufferArray(scip, &varsolvals);

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
   sepadata->tcliquedata = NULL;
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
         &sepadata->scaleval, DEFAULT_SCALEVAL, 0.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/clique/maxtreenodes",
         "maximal number of nodes in branch and bound tree (-1: no limit)",
         &sepadata->maxtreenodes, DEFAULT_MAXTREENODES, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

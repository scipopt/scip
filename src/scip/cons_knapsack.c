/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_knapsack.c,v 1.26 2004/03/15 15:54:36 bzfwolte Exp $"

/**@file   cons_knapsack.c
 * @brief  constraint handler for knapsack constraints
 * @author Tobias Achterberg
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_knapsack.h"
#include "cons_linear.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "knapsack"
#define CONSHDLR_DESC          "knapsack constraint of the form  a^T x <= b, x binary"
#define CONSHDLR_SEPAPRIORITY   +600000
#define CONSHDLR_ENFOPRIORITY   +600000
#define CONSHDLR_CHECKPRIORITY  -850000
#define CONSHDLR_SEPAFREQ             5
#define CONSHDLR_PROPFREQ             1
#define CONSHDLR_NEEDSCONS         TRUE

#define EVENTHDLR_NAME         "knapsack"
#define EVENTHDLR_DESC         "bound change event handler for knapsack constraints"

#define LINCONSUPGD_PRIORITY    +100000




/*
 * Local methods
 */

/** constraint data for knapsack constraints */
struct ConsData
{
   Real             capacity;           /**< capacity of knapsack */
   Real             weightsum;          /**< sum of all weights */
   Real             onesweightsum;      /**< sum of weights of variables fixed to one */
   VAR**            vars;               /**< variables in knapsack */
   Real*            weights;            /**< weights of knapsack items */
   EVENTDATA**      eventdatas;         /**< event datas for bound change events of the variables */
   ROW*             row;                /**< corresponding LP row */
   int              nvars;              /**< number of variables in knapsack */
   int              varssize;           /**< size of vars, weights, and eventdatas arrays */
   unsigned int     sorted:1;           /**< are the knapsack items sorted by weight? */
   unsigned int     propagated:1;       /**< is the knapsack constraint already propagated? */
};

/** event data for bound changes events */
struct EventData
{
   CONSDATA*        consdata;           /**< knapsack constraint data to process the bound change for */
   Real             weight;             /**< weight of variable */
};

/** creates event data */
static
RETCODE eventdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   EVENTDATA**      eventdata,          /**< pointer to store event data */
   CONSDATA*        consdata,           /**< constraint data */
   Real             weight              /**< weight of variable */
   )
{
   assert(eventdata != NULL);

   CHECK_OKAY( SCIPallocBlockMemory(scip, eventdata) );
   (*eventdata)->consdata = consdata;
   (*eventdata)->weight = weight;

   return SCIP_OKAY;
}  

/** frees event data */
static
RETCODE eventdataFree(
   SCIP*            scip,               /**< SCIP data structure */
   EVENTDATA**      eventdata           /**< pointer to event data */
   )
{
   assert(eventdata != NULL);

   SCIPfreeBlockMemory(scip, eventdata);

   return SCIP_OKAY;
}

/** sorts items in knapsack with nondecreasing weights */
static 
void sortItems(
  CONSDATA*        consdata            /**< constraint data */
  )
{
   int i;
   int j;
   VAR* var;
   Real weight;
   EVENTDATA* eventdata;

   assert(consdata != NULL);
   
   if( !consdata->sorted )
   {
     for( i = 0; i < consdata->nvars; i++)
     {
        var = consdata->vars[i];
        weight = consdata->weights[i];
        eventdata = consdata->eventdatas[i];
        
        for( j = i; j > 0 && weight < consdata->weights[j-1]; j--)
        {
           consdata->weights[j] = consdata->weights[j-1];
           consdata->vars[j] = consdata->vars[j-1];
           consdata->eventdatas[j] = consdata->eventdatas[j-1];
        }
        consdata->weights[j] = weight;
        consdata->vars[j] = var;
        consdata->eventdatas[j] = eventdata;
     }
     consdata->sorted = TRUE;
   } 
}

/** catches bound change events for variables in knapsack */
static
RETCODE catchEvents(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata            /**< constraint data */
   )
{
   int i;
   EVENTHDLR* eventhdlr;
   
   assert(consdata != NULL);

   eventhdlr = SCIPfindEventhdlr (scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   for( i = 0; i < consdata->nvars; i++)
   {
      CHECK_OKAY( eventdataCreate(scip, &consdata->eventdatas[i], consdata, consdata->weights[i]) );
      CHECK_OKAY( SCIPcatchVarEvent(scip, consdata->vars[i], SCIP_EVENTTYPE_LBCHANGED, eventhdlr, 
                     consdata->eventdatas[i]) );
   }

   return SCIP_OKAY;
}

/** drops bound change events for variables in knapsack */
static
RETCODE dropEvents(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata            /**< constraint data */
   )
{
   int i;
   EVENTHDLR* eventhdlr;
   
   assert(consdata != NULL);

   eventhdlr = SCIPfindEventhdlr (scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   for( i = 0; i < consdata->nvars; i++)
   {
      CHECK_OKAY( SCIPdropVarEvent(scip, consdata->vars[i], SCIP_EVENTTYPE_LBCHANGED, eventhdlr, consdata->eventdatas[i]) );
      CHECK_OKAY( eventdataFree(scip, &consdata->eventdatas[i]) );
   }

   return SCIP_OKAY;
}

/** creates knapsack constraint data */
static
RETCODE consdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to store constraint data */
   int              nvars,              /**< number of variables in knapsack */
   VAR**            vars,               /**< variables of knapsack */
   Real*            weights,            /**< weights of knapsack items */
   Real             capacity            /**< capacity of knapsack */
   )
{
   int i;

   assert(consdata != NULL);

   CHECK_OKAY( SCIPallocBlockMemory(scip, consdata) );
   CHECK_OKAY( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
   CHECK_OKAY( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->weights, weights, nvars) );
   (*consdata)->nvars = nvars;
   (*consdata)->varssize = nvars;
   (*consdata)->capacity = capacity;
   (*consdata)->row = NULL;
   (*consdata)->eventdatas = NULL;
   (*consdata)->weightsum = 0.0;
   (*consdata)->onesweightsum = 0.0;
   (*consdata)->sorted = FALSE;
   (*consdata)->propagated = FALSE;

   /* get transformed variables, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      CHECK_OKAY( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );
      CHECK_OKAY( SCIPallocBlockMemoryArray(scip, &(*consdata)->eventdatas, nvars) );
      /* catch events for variables */
      CHECK_OKAY( catchEvents(scip, *consdata) );
   } 

   /* calculate sum of weights */
   for( i = 0; i < (*consdata)->nvars; ++i )
      (*consdata)->weightsum += (*consdata)->weights[i];

   return SCIP_OKAY;
}

/** frees knapsack constraint data */
static
RETCODE consdataFree(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata            /**< pointer to the constraint data */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);

   if( (*consdata)->row != NULL )
   {
      CHECK_OKAY( SCIPreleaseRow(scip, &(*consdata)->row) );
   }
   if( (*consdata)->eventdatas != NULL )
   {
      CHECK_OKAY( dropEvents(scip, *consdata) );
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->eventdatas, (*consdata)->varssize);
   }
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->weights, (*consdata)->varssize);

   SCIPfreeBlockMemory(scip, consdata);
 
   return SCIP_OKAY;
}

/** creates LP row corresponding to knapsack constraint */
static 
RETCODE createRelaxation(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to check */
   )
{
   CONSDATA* consdata;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);

   CHECK_OKAY( SCIPcreateEmptyRow(scip, &consdata->row, SCIPconsGetName(cons), -SCIPinfinity(scip), consdata->capacity,
                  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
   CHECK_OKAY( SCIPaddVarsToRow(scip, consdata->row, consdata->nvars, consdata->vars, consdata->weights) );

   return SCIP_OKAY;
}  

/** adds linear relaxation of knapsack constraint to the LP */
static 
RETCODE addRelaxation(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to check */
   )
{
   CONSDATA* consdata;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row == NULL )
   {
      CHECK_OKAY( createRelaxation(scip, cons) );
   }
   assert(consdata->row != NULL);

   debugMessage("adding relaxation of knapsack constraint <%s> (capacity %g): ", 
      SCIPconsGetName(cons), consdata->capacity);
   debug( SCIProwPrint(consdata->row, NULL) );
   CHECK_OKAY( SCIPaddCut(scip, consdata->row, 1.0/(SCIProwGetNNonz(consdata->row)+1)) );

   return SCIP_OKAY;
}

/** checks knapsack constraint for feasibility of given solution: returns TRUE iff constraint is feasible */
static
Bool checkCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint to check */
   SOL*             sol,                /**< solution to check, NULL for current solution */
   Bool             checklprows         /**< should LP rows be checked? */
   )
{
   CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( checklprows || consdata->row == NULL || !SCIProwIsInLP(consdata->row) )
   {
      Real sum;
      int i;

      sum = 0.0;
      for( i = 0; i < consdata->nvars && sum <= consdata->capacity+0.1; i++ )
      {
         sum += consdata->weights[i] * SCIPgetSolVal(scip, sol, consdata->vars[i]);
      }
      return SCIPisFeasLE(scip, sum, consdata->capacity);
   }
   else
      return TRUE;
}

#define IDX(j,d) ((j)*(capacity+1)+(d))

/** solves knapsack problem with dynamic programming;
 *  if needed, one can provide arrays to store all selected items and all not selected items
 */
static
RETCODE solveKnapsack(
   SCIP*            scip,               /**< SCIP data structure */
   int              nitems,             /**< number of available items */
   int*             weights,            /**< item weights */
   Real*            profits,            /**< item profits */
   int              capacity,           /**< capacity of knapsack */
   int*             items,              /**< item numbers, or NULL */
   int*             solitems,           /**< array to store items in solution, or NULL */
   int*             nonsolitems,        /**< array to store items not in solution, or NULL */
   int*             nsolitems,          /**< pointer to store number of items in solution, or NULL */
   int*             nnonsolitems,       /**< pointer to store number of items not in solution, or NULL */
   Real*            solval              /**< pointer to store optimal solution value */
) 
{
   Real* optvalues;
   int d;
   int j;

   assert(weights != NULL);
   assert(profits != NULL);
   assert(solval != NULL);
   assert(capacity >= 0);
   assert(nitems >= 0);

   CHECK_OKAY( SCIPallocBufferArray(scip, &optvalues, (nitems+1)*(capacity+1)) );
   
   /* fill dynamic programming table with optimal values */
   for( d = 0; d <= capacity; d++ )
      optvalues[IDX(0,d)] = 0.0;
   for( j = 1; j <= nitems; j++ )
   {
      for( d = 0; d < weights[j-1] && d <= capacity; d++)
         optvalues[IDX(j,d)] = optvalues[IDX(j-1,d)];
      for( d = weights[j-1]; d <= capacity; d++ )
      {
         if( optvalues[IDX(j-1,d-weights[j-1])]+profits[j-1] > optvalues[IDX(j-1,d)] )
            optvalues[IDX(j,d)] = optvalues[IDX(j-1,d-weights[j-1])] + profits[j-1];
         else
            optvalues[IDX(j,d)] = optvalues[IDX(j-1,d)];
      } 
   }

   /* produce optimal solution by following the table */
   if( solitems != NULL)
   {
      assert(items != NULL);
      assert(nsolitems != NULL);
      assert(nonsolitems != NULL);
      assert(nnonsolitems != NULL);

      *nnonsolitems = 0;
      *nsolitems = 0;
      d = capacity;
      
      for( j = nitems; j > 0; j-- )
      {
         if( optvalues[IDX(j,d)] > optvalues[IDX(j-1,d)] )
         {
            solitems[*nsolitems] = items[j-1];
            (*nsolitems)++;
            d -= weights[j-1];
         } 
         else
         { 
            nonsolitems[*nnonsolitems] = items[j-1];
            (*nnonsolitems)++;
         }
         assert(d >= 0);
      }
      assert(*nsolitems + *nnonsolitems == nitems);
   }

   *solval = optvalues[IDX(nitems,capacity)];

   CHECK_OKAY( SCIPfreeBufferArray(scip, &optvalues) );

   return SCIP_OKAY;
}

/** lifts given knapsack cover */
static
RETCODE liftCover(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< cover inequality */
   CONS*            cons,               /**< knapsack constraint */
   int*             covervars,          /**< cover elements */
   int*             noncovervars,       /**< non-cover elements */
   int              ncovervars,         /**< number of cover elements */
   int              nnoncovervars       /**< number of non-cover elements */
   )
{
   CONSDATA* consdata;
   int* minweight;
   int* weights;
   int weight;
   int rescapacity;
   int beta;
   int i;
   int j;
   int z;
   int left;
   int right;
   int middle;

   assert(covervars != NULL);
   assert(noncovervars != NULL);
   assert(ncovervars > 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   CHECK_OKAY( SCIPallocBufferArray(scip, &minweight, ncovervars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &weights, ncovervars) );
  
   /* sort weights of covervars */
   for( i = 0; i < ncovervars; i++ )
   {
      weight = consdata->weights[covervars[i]];
      for( j = i; j > 0 && weight < weights[j-1]; j--)
         weights[j] = weights[j-1];
      weights[j] = weight;   
   }

   /* calculate minweight vector: 
    * minweight[z] := minimal sum of weights s.t. activity of cut inequality equals z
    */
   minweight[0] = 0;
   for( z = 1; z < ncovervars; z++ )
      minweight[z] = minweight[z-1] + weights[z-1];

   /* calculate lifting coefficients beta for noncovervars:
    * for each noncovervar i: 
    * 1. calculate maximal activity z_max of current cut inequality s.t. the knapsack is still feasible with x_i = 1 
    * 2. add x_i with lifting coefficient beta_i = ncovervars - 1 - z_max to cut inequality
    * 3. update minweight table: calculate minimal sum of weights s.t. activity of current cut inequality equals z
    */
   for( i = 0; i < nnoncovervars; i++ )
   {
      /* binary search in sorted minweight array for the largest entry that is not greater then capacity - weight_i */
      weight = consdata->weights[noncovervars[i]];
      rescapacity = consdata->capacity - weight;
      left = 0;
      right = ncovervars;
      while( left < right - 1)
      {
         middle = (left + right) / 2;
         if( minweight[middle] <= rescapacity ) 
            left = middle;
         else
            right = middle;
      }
      assert(left == right - 1);
      
      /* now zmax_i = left: calculate beta_i */
      beta = ncovervars - 1 - left;

      if( beta == 0 )
         continue;

      /* insert variable with coefficient beta into cut */
      CHECK_OKAY( SCIPaddVarToRow(scip, row, consdata->vars[noncovervars[i]], (Real) beta) );
      
      /* update minweight table: all entries below beta keep the same */
      for( z = ncovervars - 1; z >= beta; z-- )
         minweight[z] = MIN(minweight[z], minweight[z - beta] + weight);
   }

   CHECK_OKAY( SCIPfreeBufferArray(scip, &weights) );
   CHECK_OKAY( SCIPfreeBufferArray(scip, &minweight) );

   return SCIP_OKAY;
}
       
/** separates cover inequalities for given knapsack constraint */
static
RETCODE separateCovers(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< knapsack constraint */
   Bool*            separated           /**< pointer to store whether a cut was found */
   )
{
   CONSDATA* consdata;
   int* items;
   int* weights;
   Real* profits;
   int i;
   int nitems;
   int* fixedones;
   int nfixedones;
   int* fixedzeros;
   int nfixedzeros;
   int capacity;
   int* covervars;
   int ncovervars;
   int* noncovervars;
   int nnoncovervars;
   Real infeasibility;
   Real solval;
   Real solvalsum;
   Real transsol;
   int coverweight;
   int fixedonesweightsum;

   assert(separated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *separated = FALSE;
   
   /* solve the following knapsack problem with dynamic programming:
    * max SUM (1 - x_i^*) z_i
    *     SUM w_i z_i >= SUM w_i - capacity - 1
    */
   
   /* set up datastructures */
   CHECK_OKAY( SCIPallocBufferArray(scip, &items, consdata->nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &weights, consdata->nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &profits, consdata->nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &fixedones, consdata->nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &fixedzeros, consdata->nvars) );
   
   /* use all variables with fractional LP value */
   nitems = 0;
   nfixedones = 0;
   nfixedzeros = 0;
   solvalsum = 0.0;
   capacity = -(int)(consdata->capacity+0.5) - 1;
   fixedonesweightsum = 0;
   for( i = 0; i < consdata->nvars; i++ )
   {
      solval = SCIPgetVarSol(scip, consdata->vars[i]);
      solvalsum += solval;
      if( !SCIPisIntegral(scip, solval) )
      {
         items[nitems] = i;
         weights[nitems] = (int)(consdata->weights[i]+0.5);
         profits[nitems] = 1.0 - solval;
         nitems++;
         capacity += (int)(consdata->weights[i]+0.5);
      }
      else if( solval > 0.5 )
      {
         fixedones[nfixedones] = i;
         nfixedones++;
         capacity += (int)(consdata->weights[i]+0.5);
         fixedonesweightsum += (int)(consdata->weights[i]+0.5);
      }
      else
      {
         fixedzeros[nfixedzeros] = i;
         nfixedzeros++;
      }

   }

   if( capacity >= 0)
   {
      int* covervars;
      int ncovervars;

      /* solve separation knapsack with dynamic programming */
      CHECK_OKAY( SCIPallocBufferArray(scip, &covervars, consdata->nvars) );
      CHECK_OKAY( SCIPallocBufferArray(scip, &noncovervars, consdata->nvars) );
      CHECK_OKAY( solveKnapsack(scip, nitems, weights, profits, capacity, 
                     items, noncovervars, covervars, &nnoncovervars, &ncovervars, &transsol) ); 
      
      /* generate cutting plane */
      infeasibility = transsol - nitems - nfixedones + 1.0 + solvalsum;
      if( SCIPisFeasPositive(scip, infeasibility) ) 
      {
         ROW* row;
         char name[MAXSTRLEN];

         /* remove unnecessary items from cover to get a minimal cover */
         coverweight = 0;
         for( i = 0; i < ncovervars; i++)
            coverweight += consdata->weights[covervars[i]];
                  
         for( i = 0; i < nfixedones; i++)
         {
            fixedonesweightsum -= consdata->weights[fixedones[i]];
            if( coverweight + fixedonesweightsum <= consdata->capacity )
            {
               covervars[ncovervars] = fixedones[i];
               ncovervars++;
               coverweight += consdata->weights[fixedones[i]];
            }  
         }
         assert(coverweight > consdata->capacity);

         /* add variables with LP value zero to non-cover vars */
         for( i = 0; i < nfixedzeros; i++)
         {
            noncovervars[nnoncovervars] = fixedzeros[i];
            nnoncovervars++;
         }  

         /* create LP row */
         sprintf(name, "%s_%lld", SCIPconsGetName(cons), SCIPconshdlrGetNCutsFound(SCIPconsGetHdlr(cons)));
         CHECK_OKAY( SCIPcreateEmptyRow (scip, &row, name, -SCIPinfinity(scip), ncovervars-1.0, 
                        SCIPconsIsLocal(cons), FALSE, SCIPconsIsRemoveable(cons)) );
         for( i = 0; i < ncovervars; i++)
         {
            CHECK_OKAY( SCIPaddVarToRow(scip, row, consdata->vars[covervars[i]], 1.0) );
         }
      
         debugMessage("found cover cut for knapsack constraint <%s>: ", SCIPconsGetName(cons));
         debug(SCIProwPrint(row, NULL));
        
         /* lift variables not in cover to inequality */
         CHECK_OKAY( liftCover(scip, row, cons, covervars, noncovervars, ncovervars, nnoncovervars) ); 
           
         debugMessage("lifted cover cut for knapsack constraint <%s>: ", SCIPconsGetName(cons));
         debug(SCIProwPrint(row, NULL));
         CHECK_OKAY( SCIPaddCut(scip, row, infeasibility/(SCIProwGetNNonz(row)+1)) );
         CHECK_OKAY( SCIPreleaseRow(scip, &row) );
         *separated = TRUE;
      }

      CHECK_OKAY( SCIPfreeBufferArray(scip, &noncovervars) );
      CHECK_OKAY( SCIPfreeBufferArray(scip, &covervars) );
   }

   /* free temporary data */
   CHECK_OKAY( SCIPfreeBufferArray(scip, &fixedzeros) );
   CHECK_OKAY( SCIPfreeBufferArray(scip, &fixedones) );
   CHECK_OKAY( SCIPfreeBufferArray(scip, &profits) );
   CHECK_OKAY( SCIPfreeBufferArray(scip, &weights) );
   CHECK_OKAY( SCIPfreeBufferArray(scip, &items) );

   return SCIP_OKAY;
}

/** propagation methode for knapsack constraint */
static
RETCODE propagateCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< knapsack constraint */
   Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   Bool*            redundant,          /**< pointer to store whether constraint is redundant */
   int*             nfixedvars          /**< pointer to count number of fixings */
   )
{
   CONSDATA* consdata;
   Bool infeasible;
   Bool tightened;
   Real zerosweightsum;
   int i;

   assert(cutoff != NULL);
   assert(redundant != NULL);
   assert(nfixedvars != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *cutoff = FALSE;
   *redundant = FALSE;

   /* check, if constraint is already propagated */
   if( consdata->propagated )
      return SCIP_OKAY;

   /* check, if weights of fixed variables already exceeds knapsack capacity */
   if( consdata->capacity < consdata->onesweightsum )
   {
      /* start conflict analysis with the fixed-to-one variables */
      CHECK_OKAY( SCIPinitConflictAnalysis(scip) );
      for( i = 0; i < consdata->nvars; i++ )
      {
         if( SCIPvarGetLbLocal(consdata->vars[i]) > 0.5)
         {
            CHECK_OKAY( SCIPaddConflictVar(scip, consdata->vars[i]) );
         }
      }

      CHECK_OKAY( SCIPanalyzeConflict(scip, NULL) );
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* make sure, the items are sorted by non-decreasing weight */
   sortItems(consdata);

   /* fix all variables to zero, that don't fit into the knapsack anymore */
   zerosweightsum = 0.0;
   for( i = consdata->nvars - 1; i >= 0; i-- )
   {
      if( consdata->weights[i] > consdata->capacity - consdata->onesweightsum )
      {
         if( SCIPvarGetLbLocal(consdata->vars[i]) < 0.5 && SCIPvarGetUbLocal(consdata->vars[i]) > 0.5 )
         {
            CHECK_OKAY( SCIPinferBinVar(scip, consdata->vars[i], FALSE, cons, 0, &infeasible, &tightened) );
            assert(!infeasible);
            assert(tightened);
            (*nfixedvars)++;
         }
         zerosweightsum += consdata->weights[i];
      }
      else
         break;
   }

   /* if the remaining (potentially unfixed) variables would fit all into the knapsack, the knapsack is now redundant */
   if( consdata->weightsum - zerosweightsum <= consdata->capacity )
   {
      CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
      *redundant = TRUE;
   }

   /* mark the constraint propagated */
   consdata->propagated = TRUE;

   return SCIP_OKAY;
}

/** deletes coefficient at given position from constraint data */
static
RETCODE delCoefPos(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< knapsack constraint */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< position of coefficient to delete */
   )
{
   CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(consdata->vars[pos]));

   /* if necessary, update the rounding locks of variable */
   if( SCIPconsIsLocked(cons) )
   {
      assert(SCIPconsIsTransformed(cons));
      SCIPvarUnlock(consdata->vars[pos], (int)SCIPconsIsLockedNeg(cons), (int)SCIPconsIsLockedPos(cons));
   }

   if( SCIPconsIsTransformed(cons) )
   {
      /* drop bound tighten events */
      CHECK_OKAY( SCIPdropVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_LBCHANGED, eventhdlr, 
                     consdata->eventdatas[pos]) );
      CHECK_OKAY( eventdataFree(scip, &consdata->eventdatas[pos]) );
   }

   /* update weight sums */
   consdata->weightsum -= consdata->weights[pos];
   if( SCIPvarGetLbLocal(consdata->vars[pos]) > 0.5 )
      consdata->onesweightsum -= consdata->weights[pos];

   /* move the last variable to the free slot */
   consdata->vars[pos] = consdata->vars[consdata->nvars-1];
   consdata->weights[pos] = consdata->weights[consdata->nvars-1];
   consdata->eventdatas[pos] = consdata->eventdatas[consdata->nvars-1];
   consdata->nvars--;

   consdata->propagated = FALSE;

   return SCIP_OKAY;
}

/** deletes all fixed variables from knapsack constraint */
static
RETCODE applyFixings(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< knapsack constraint */
   EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   CONSDATA* consdata;
   VAR* var;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);

   v = 0;
   while( v < consdata->nvars )
   {
      var = consdata->vars[v];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);

      if( SCIPvarGetLbGlobal(var) > 0.5 )
      {
         assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(var), 1.0));
         consdata->capacity -= consdata->weights[v];
         CHECK_OKAY( delCoefPos(scip, cons, eventhdlr, v) );
      }
      else if( SCIPvarGetUbGlobal(var) < 0.5 )
      {
         assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 0.0));
         CHECK_OKAY( delCoefPos(scip, cons, eventhdlr, v) );
      }
      else
         ++v;
   }

   return SCIP_OKAY;
}





/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with an #if 0 ... #else #define ... are optional */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#if 0
static
DECL_CONSFREE(consFreeKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consFreeKnapsack NULL
#endif


/** initialization method of constraint handler (called when problem solving starts) */
#if 0
static
DECL_CONSINIT(consInitKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitKnapsack NULL
#endif


/** deinitialization method of constraint handler (called when problem solving exits) */
#if 0
static
DECL_CONSEXIT(consExitKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitKnapsack NULL
#endif


/** solving start notification method of constraint handler (called when presolving was finished) */
#if 0
static
DECL_CONSSOLSTART(consSolstartKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSolstartKnapsack NULL
#endif


/** frees specific constraint data */
static
DECL_CONSDELETE(consDeleteKnapsack)
{  /*lint --e{715}*/
   CHECK_OKAY( consdataFree(scip, consdata) );
   
   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
DECL_CONSTRANS(consTransKnapsack)
{  /*lint --e{715}*/
   CONSDATA* sourcedata;
   CONSDATA* targetdata;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* create target constraint data */
   CHECK_OKAY( consdataCreate(scip, &targetdata, sourcedata->nvars, sourcedata->vars, sourcedata->weights, 
                  sourcedata->capacity) ); 

   /* create target constraint */
   CHECK_OKAY( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
                  SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
                  SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
                  SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), SCIPconsIsRemoveable(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
DECL_CONSINITLP(consInitlpKnapsack)
{  /*lint --e{715}*/
   int i;

   for( i = 0; i < nconss; i++ )
   {
      if( SCIPconsIsInitial(conss[i]) )
      {
         CHECK_OKAY( addRelaxation(scip, conss[i]) );
      }
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler */
static
DECL_CONSSEPA(consSepaKnapsack)
{  /*lint --e{715}*/
   int i;
   Bool separated;
   *result=SCIP_DIDNOTFIND;
   
   for(i=0; i<nusefulconss; i++)
   {
      CHECK_OKAY( separateCovers(scip, conss[i], &separated) );
      if( separated ) 
         *result=SCIP_SEPARATED;
   }
   
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
DECL_CONSENFOLP(consEnfolpKnapsack)
{  /*lint --e{715}*/
   int i;

   for( i = 0; i < nconss; i++ )
   {
      if( !checkCons(scip, conss[i], NULL, FALSE) )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   } 
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
DECL_CONSENFOPS(consEnfopsKnapsack)
{  /*lint --e{715}*/
   int i;

   for( i = 0; i < nconss; i++ )
   {
      if( !checkCons(scip, conss[i], NULL, TRUE) )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   } 
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;  
}


/** feasibility check method of constraint handler for integral solutions */
static
DECL_CONSCHECK(consCheckKnapsack)
{  /*lint --e{715}*/
   int i;

   for( i = 0; i < nconss; i++ )
   {
      if( !checkCons(scip, conss[i], sol, checklprows) )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   } 
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
DECL_CONSPROP(consPropKnapsack)
{  /*lint --e{715}*/
   Bool cutoff;
   Bool redundant;
   int nfixedvars;
   int i;

   cutoff = FALSE;
   nfixedvars = 0;

   for( i = 0; i < nusefulconss && !cutoff; i++ )
   {
      CHECK_OKAY( propagateCons(scip, conss[i], &cutoff, &redundant, &nfixedvars) );
   } 

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nfixedvars > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** presolving method of constraint handler */
static
DECL_CONSPRESOL(consPresolKnapsack)
{  /*lint --e{715}*/
   EVENTHDLR* eventhdlr;
   Bool cutoff;
   Bool redundant;
   int oldnfixedvars;
   int oldndelconss;
   int i;

   cutoff = FALSE;
   oldnfixedvars = *nfixedvars;
   oldndelconss = *ndelconss;

   /* get event handler for bound change events */
   eventhdlr = SCIPfindEventhdlr (scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);
   
   for( i = 0; i < nconss && !cutoff; i++ )
   {
      /* remove all variables that are fixed to zero, check redundancy due to fixed-to-one variable */
      CHECK_OKAY( applyFixings(scip, conss[i], eventhdlr) );

      /* propagate constraint */
      CHECK_OKAY( propagateCons(scip, conss[i], &cutoff, &redundant, nfixedvars) );
      if( redundant )
         (*ndelconss)++;
   } 

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( *nfixedvars > oldnfixedvars || *ndelconss > oldndelconss )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** conflict variable resolving method of constraint handler */
static
DECL_CONSRESCVAR(consRescvarKnapsack)
{  /*lint --e{715}*/
   CONSDATA* consdata;
   int i;
   consdata = SCIPconsGetData(cons);
 
   assert(SCIPvarGetUbLocal(infervar) < 0.5);
   assert(consdata != NULL);

   for( i = 0; i < consdata->nvars; i++ )
   {
      if( SCIPvarWasFixedEarlier(consdata->vars[i], infervar) && SCIPvarGetLbLocal(consdata->vars[i]) > 0.5)
      {
         CHECK_OKAY( SCIPaddConflictVar(scip, consdata->vars[i]) );
      }
   }

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
DECL_CONSLOCK(consLockKnapsack)
{  /*lint --e{715}*/
   CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   for( i = 0; i < consdata->nvars; i++)
   {
      SCIPvarLock(consdata->vars[i], nlocksneg, nlockspos);
   }

   return SCIP_OKAY;
}


/** variable rounding unlock method of constraint handler */
static
DECL_CONSUNLOCK(consUnlockKnapsack)
{  /*lint --e{715}*/
   CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   for( i = 0; i < consdata->nvars; i++)
   {
      SCIPvarUnlock(consdata->vars[i], nunlocksneg, nunlockspos);
   }   

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
DECL_CONSACTIVE(consActiveKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveKnapsack NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
DECL_CONSDEACTIVE(consDeactiveKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveKnapsack NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
DECL_CONSENABLE(consEnableKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableKnapsack NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
DECL_CONSDISABLE(consDisableKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableKnapsack NULL
#endif




/*
 * Linear constraint upgrading
 */


/** creates and captures a knapsack constraint out of a linear inequality */
static
RETCODE createNormalizedKnapsack(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Real*            vals,               /**< array with inequality coefficients */
   Real             lhs,                /**< left hand side of inequality */
   Real             rhs,                /**< right hand side of inequality */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   VAR** transvars;
   Real* transvals;
   Real capacity;
   Real val;
   int mult;
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);
   assert(SCIPisInfinity(scip, -lhs) != SCIPisInfinity(scip, rhs));

   /* get temporary memory */
   CHECK_OKAY( SCIPallocBufferArray(scip, &transvars, nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &transvals, nvars) );

   /* if the right hand side is non-infinite, we have to negate all variables with negative coefficient;
    * otherwise, we have to negate all variables with positive coefficient and multiply the row with -1
    */
   if( SCIPisInfinity(scip, rhs) )
   {
      mult = -1;
      capacity = -lhs;
   }
   else
   {
      mult = +1;
      capacity = rhs;
   }

   /* negate positive or negative variables */
   for( v = 0; v < nvars; ++v )
   {
      val = mult * vals[v];
      if( val > 0.0 )
      {
         transvars[v] = vars[v];
         transvals[v] = val;
      }
      else
      {
         CHECK_OKAY( SCIPgetNegatedVar(scip, vars[v], &transvars[v]) );
         transvals[v] = -val;
         capacity -= val;
      }
      assert(transvars[v] != NULL);
   }

   /* create the constraint */
   CHECK_OKAY( SCIPcreateConsKnapsack(scip, cons, name, nvars, transvars, transvals, capacity,
                  initial, separate, enforce, check, propagate, local, modifiable, removeable) );

   /* free temporary memory */
   CHECK_OKAY( SCIPfreeBufferArray(scip, &transvals) );
   CHECK_OKAY( SCIPfreeBufferArray(scip, &transvars) );

   return SCIP_OKAY;
}

#ifdef LINCONSUPGD_PRIORITY
static
DECL_LINCONSUPGD(linconsUpgdKnapsack)
{  /*lint --e{715}*/
   Bool upgrade;

   assert(upgdcons != NULL);
   
   /* check, if linear constraint can be upgraded to a knapsack constraint
    * - all variables must be binary
    * - all coefficients must be integral
    * - exactly one of the sides must be infinite
    */
   upgrade = (nposbin + nnegbin == nvars)
      && (ncoeffspone + ncoeffsnone + ncoeffspint + ncoeffsnint == nvars)
      && (SCIPisInfinity(scip, -lhs) != SCIPisInfinity(scip, rhs));

   if( upgrade )
   {
      debugMessage("upgrading constraint <%s> to knapsack constraint\n", SCIPconsGetName(cons));
      
      /* create the bin Knapsack constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      CHECK_OKAY( createNormalizedKnapsack(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, lhs, rhs,
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
                     SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
                     SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
   }

   return SCIP_OKAY;
}
#endif

/** execution methode of bound change event handler */
static
DECL_EVENTEXEC(eventExecKnapsack)
{
   assert(eventdata != NULL);
   assert(eventdata->consdata != NULL);
   
   switch( SCIPeventGetType(event) )
   {
   case SCIP_EVENTTYPE_LBTIGHTENED:
      eventdata->consdata->onesweightsum += eventdata->weight;
      eventdata->consdata->propagated = FALSE;
      break;
   case SCIP_EVENTTYPE_LBRELAXED:
      eventdata->consdata->onesweightsum -= eventdata->weight;
      break;
   default:
      errorMessage("Invalid event type %x\n", SCIPeventGetType(event));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/*
 * constraint specific interface methods
 */

/** creates the handler for knapsack constraints and includes it in SCIP */
RETCODE SCIPincludeConshdlrKnapsack(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;
   EVENTHDLRDATA* eventhdlrdata;

   /* create knapsack constraint handler data */
   conshdlrdata = NULL;

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
                  CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_NEEDSCONS,
                  consFreeKnapsack, consInitKnapsack, consExitKnapsack, consSolstartKnapsack,
                  consDeleteKnapsack, consTransKnapsack, consInitlpKnapsack,
                  consSepaKnapsack, consEnfolpKnapsack, consEnfopsKnapsack, consCheckKnapsack, 
                  consPropKnapsack, consPresolKnapsack, consRescvarKnapsack,
                  consLockKnapsack, consUnlockKnapsack,
                  consActiveKnapsack, consDeactiveKnapsack, 
                  consEnableKnapsack, consDisableKnapsack,
                  conshdlrdata) );
  
   /* include event handler for bound change events */
   eventhdlrdata = NULL;
   CHECK_OKAY( SCIPincludeEventhdlr (scip, EVENTHDLR_NAME, EVENTHDLR_DESC, 
                  NULL, NULL, NULL, NULL, eventExecKnapsack,
                  eventhdlrdata) );

#ifdef LINCONSUPGD_PRIORITY
   /* include the linear constraint upgrade in the linear constraint handler */
   CHECK_OKAY( SCIPincludeLinconsUpgrade(scip, linconsUpgdKnapsack, LINCONSUPGD_PRIORITY) );
#endif

   return SCIP_OKAY;
}

/** creates and captures a knapsack constraint */
RETCODE SCIPcreateConsKnapsack(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              len,                /**< number of nonzeros in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Real*            vals,               /**< array with coefficients of constraint entries */
   Real             rhs,                /**< right hand side of constraint */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   Bool             removeable          /**< should the constraint be removed from the LP due to aging or cleanup? */
   )
{
   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   /* find the knapsack constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("knapsack constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   CHECK_OKAY( consdataCreate(scip, &consdata, len, vars, vals, rhs) );
        
   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
                  local, modifiable, removeable) );

   return SCIP_OKAY;
}

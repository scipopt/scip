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
#pragma ident "@(#) $Id: cons_knapsack.c,v 1.36 2004/04/16 10:43:33 bzfwolte Exp $"

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
#define CONSHDLR_SEPAFREQ            10
#define CONSHDLR_PROPFREQ             1
#define CONSHDLR_NEEDSCONS         TRUE

#define EVENTHDLR_NAME         "knapsack"
#define EVENTHDLR_DESC         "bound change event handler for knapsack constraints"

#define LINCONSUPGD_PRIORITY    +100000

#define MAX_SEPA_CAPACITY         10000 /**< maximal capacity of knapsack to apply cover separation */

#define DEFAULT_MAXROUNDS             5 /**< maximal number of separation rounds per node */
#define DEFAULT_MAXROUNDSROOT        10 /**< maximal number of separation rounds in the root node */
#define DEFAULT_MAXSEPACUTS          50 /**< maximal number of cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT     200 /**< maximal number of cuts separated per separation round in root node */




/*
 * Local methods
 */

/** constraint handler data */
struct ConshdlrData
{
   int              maxrounds;          /**< maximal number of separation rounds per node */
   int              maxroundsroot;      /**< maximal number of separation rounds in the root node */
   int              maxsepacuts;        /**< maximal number of cuts separated per separation round */
   int              maxsepacutsroot;    /**< maximal number of cuts separated per separation round in root node */
};

/** constraint data for knapsack constraints */
struct ConsData
{
   VAR**            vars;               /**< variables in knapsack */
   Longint*         weights;            /**< weights of knapsack items */
   EVENTDATA**      eventdatas;         /**< event datas for bound change events of the variables */
   ROW*             row;                /**< corresponding LP row */
   int              nvars;              /**< number of variables in knapsack */
   int              varssize;           /**< size of vars, weights, and eventdatas arrays */
   Longint          capacity;           /**< capacity of knapsack */
   Longint          weightsum;          /**< sum of all weights */
   Longint          onesweightsum;      /**< sum of weights of variables fixed to one */
   unsigned int     sorted:1;           /**< are the knapsack items sorted by weight? */
   unsigned int     propagated:1;       /**< is the knapsack constraint already propagated? */
};

/** event data for bound changes events */
struct EventData
{
   CONSDATA*        consdata;           /**< knapsack constraint data to process the bound change for */
   Longint          weight;             /**< weight of variable */
};

/** creates event data */
static
RETCODE eventdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   EVENTDATA**      eventdata,          /**< pointer to store event data */
   CONSDATA*        consdata,           /**< constraint data */
   Longint          weight              /**< weight of variable */
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
   Longint weight;
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
   Longint*         weights,            /**< weights of knapsack items */
   Longint          capacity            /**< capacity of knapsack */
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
   (*consdata)->weightsum = 0;
   (*consdata)->onesweightsum = 0;
   (*consdata)->sorted = FALSE;
   (*consdata)->propagated = FALSE;

   /* get transformed variables, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      CHECK_OKAY( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );

      /* catch events for variables */
      CHECK_OKAY( SCIPallocBlockMemoryArray(scip, &(*consdata)->eventdatas, nvars) );
      CHECK_OKAY( catchEvents(scip, *consdata) );
   } 

   /* calculate sum of weights */
   for( i = 0; i < (*consdata)->nvars; ++i )
   {
      (*consdata)->weightsum += (*consdata)->weights[i];
      if( SCIPvarGetLbLocal((*consdata)->vars[i]) > 0.5 )
         (*consdata)->onesweightsum += (*consdata)->weights[i];
   }

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

/** changes a single weight in knapsack constraint data */
static
void consdataChgWeight(
   CONSDATA*        consdata,           /**< knapsack constraint data */
   int              item,               /**< item number */
   Longint          newweight           /**< new weight of item */
   )
{
   Longint oldweight;

   assert(consdata != NULL);
   assert(0 <= item && item < consdata->nvars);
   assert(newweight > 0);

   oldweight = consdata->weights[item];
   consdata->weights[item] = newweight;
   consdata->weightsum += (newweight - oldweight);

   if( SCIPvarGetLbLocal(consdata->vars[item]) > 0.5 )
      consdata->onesweightsum += (newweight - oldweight);
      
   if( consdata->eventdatas != NULL )
   {
      assert(consdata->eventdatas[item] != NULL);
      assert(consdata->eventdatas[item]->weight == oldweight);
      consdata->eventdatas[item]->weight = newweight;
   }

   consdata->propagated = FALSE;
   consdata->sorted = FALSE;
}

/** creates LP row corresponding to knapsack constraint */
static 
RETCODE createRelaxation(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to check */
   )
{
   CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);

   CHECK_OKAY( SCIPcreateEmptyRow(scip, &consdata->row, SCIPconsGetName(cons),
                  -SCIPinfinity(scip), (Real)consdata->capacity,
                  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
   for( i = 0; i < consdata->nvars; ++i )
   {
      CHECK_OKAY( SCIPaddVarToRow(scip, consdata->row, consdata->vars[i], (Real)consdata->weights[i]) );
   }

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

   debugMessage("adding relaxation of knapsack constraint <%s> (capacity %lld): ", 
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

   debugMessage("checking knapsack constraint <%s> for feasibility of solution %p (lprows=%d)\n",
      SCIPconsGetName(cons), sol, checklprows);

   if( checklprows || consdata->row == NULL || !SCIProwIsInLP(consdata->row) )
   {
      Real sum;
      int i;

      sum = 0.0;
      for( i = 0; i < consdata->nvars && sum <= consdata->capacity + 0.1; i++ )
      {
         sum += consdata->weights[i] * SCIPgetSolVal(scip, sol, consdata->vars[i]);
      }
      return SCIPisFeasLE(scip, sum, (Real)consdata->capacity);
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
   Longint*         weights,            /**< item weights */
   Real*            profits,            /**< item profits */
   Longint          capacity,           /**< capacity of knapsack */
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

/** lifts given cardinality inequality sum(x_i) <= c */
static
RETCODE liftCardinality(
   SCIP*            scip,               /**< SCIP data structure */
   int*             liftcoefs,          /**< to store lifting coefficient of non-set elements */
   CONS*            cons,               /**< knapsack constraint */
   int*             setvars,            /**< set elements */
   int*             nonsetvars,         /**< non-set elements */
   int              nsetvars,           /**< number of set elements */
   int              nnonsetvars,        /**< number of non-set elements */
   int              maxcardinality,     /**< maximal cardinality of selected subset in given set */ 
   Real*            liftlpval           /**< LP solution value of lifted elements */  
  )
{
   CONSDATA* consdata;
   Longint* minweight;
   Longint* weights;
   Longint weight;
   Longint rescapacity;
   int i;
   int j;
   int z;
   int left;
   int right;
   int middle;
   Real solval;

   assert(setvars != NULL);
   assert(nonsetvars != NULL);
   assert(nsetvars > 0);
   assert(liftcoefs != NULL);
   assert(liftlpval != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   CHECK_OKAY( SCIPallocBufferArray(scip, &minweight, maxcardinality + 1) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &weights, nsetvars) );
  
   /* sort weights of setvars */
   for( i = 0; i < nsetvars; i++ )
   {
      weight = consdata->weights[setvars[i]];
      for( j = i; j > 0 && weight < weights[j-1]; j--)
         weights[j] = weights[j-1];
      weights[j] = weight;   
   }

   /* calculate minweight vector: 
    * minweight[z] := minimal sum of weights s.t. activity of cut inequality equals z
    */
   minweight[0] = 0;
   for( z = 1; z <= maxcardinality; z++ )
      minweight[z] = minweight[z-1] + weights[z-1];

   /* calculate lifting coefficients beta for nonsetvars:
    * for each nonsetvar i: 
    * 1. calculate maximal activity z_max of current cut inequality s.t. the knapsack is still feasible with x_i = 1 
    * 2. add x_i with lifting coefficient beta_i = maxcardinality - z_max to cut inequality
    * 3. update minweight table: calculate minimal sum of weights s.t. activity of current cut inequality equals z
    */
   *liftlpval = 0.0;
   for( i = 0; i < nnonsetvars; i++ )
   {
      /* binary search in sorted minweight array for the largest entry that is not greater then capacity - weight_i */
      weight = consdata->weights[nonsetvars[i]];
      rescapacity = consdata->capacity - weight;
      left = 0;
      right = maxcardinality + 1;
      while( left < right - 1)
      {
         middle = (left + right) / 2;
         if( minweight[middle] <= rescapacity ) 
            left = middle;
         else
            right = middle;
      }
      assert(left == right - 1);
      
      /* now zmax_i = left: calculate lifting coefficient liftcoefs[i] = beta_i */
      liftcoefs[i] = maxcardinality - left;

      if( liftcoefs[i] == 0 )
         continue;

      /* update slack of cut inequality */
      solval = SCIPgetVarSol(scip, consdata->vars[nonsetvars[i]]);
      *liftlpval += liftcoefs[i] * solval;
      
      /* update minweight table: all entries below beta keep the same */
      for( z = maxcardinality; z >= liftcoefs[i]; z-- )
         minweight[z] = MIN(minweight[z], minweight[z - liftcoefs[i]] + weight);
   }

   CHECK_OKAY( SCIPfreeBufferArray(scip, &weights) );
   CHECK_OKAY( SCIPfreeBufferArray(scip, &minweight) );

   return SCIP_OKAY;
}

#if 0
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
   CHECK_OKAY( liftCardinality(scip, row, cons, covervars, noncovervars, ncovervars, nnoncovervars, ncovervars - 1) );

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
   Longint* weights;
   Real* profits;
   Longint capacity;
   int i;
   int nitems;
   int* fixedones;
   int nfixedones;
   int* fixedzeros;
   int nfixedzeros;
   int* covervars;
   int ncovervars;
   int* noncovervars;
   int nnoncovervars;
   Real infeasibility;
   Real solval;
   Real solvalsum;
   Real transsol;
   Longint coverweight;
   Longint fixedonesweightsum;

   assert(separated != NULL);

   debugMessage("separating knapsack constraint <%s>\n", SCIPconsGetName(cons));

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
         weights[nitems] = consdata->weights[i];
         profits[nitems] = 1.0 - solval;
         nitems++;
         capacity += consdata->weights[i];
      }
      else if( solval > 0.5 )
      {
         fixedones[nfixedones] = i;
         nfixedones++;
         capacity += consdata->weights[i];
         fixedonesweightsum += consdata->weights[i];
      }
      else
      {
         fixedzeros[nfixedzeros] = i;
         nfixedzeros++;
      }
   }

   if( capacity >= 0 && capacity <= MAX_SEPA_CAPACITY )
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
         Real cutnorm;
         Real cutfeas;

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
            else
            {
               noncovervars[nnoncovervars] = fixedones[i];
               nnoncovervars++;
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
         CHECK_OKAY( SCIPcreateEmptyRow (scip, &row, name, -SCIPinfinity(scip), (Real)ncovervars-1.0, 
                        SCIPconsIsLocal(cons), FALSE, SCIPconsIsRemoveable(cons)) );
         for( i = 0; i < ncovervars; i++)
         {
            CHECK_OKAY( SCIPaddVarToRow(scip, row, consdata->vars[covervars[i]], 1.0) );
         }
      
         debugMessage("found cover cut for knapsack constraint <%s>: ", SCIPconsGetName(cons));
         debug(SCIProwPrint(row, NULL));
        
         /* lift variables not in cover into inequality */
         CHECK_OKAY( liftCover(scip, row, cons, covervars, noncovervars, ncovervars, nnoncovervars) ); 

         /* check, if cut is violated enough */
         cutnorm = SCIProwGetNorm(row);
         cutfeas = SCIPgetRowLPFeasibility(scip, row);
         if( SCIPisFeasNegative(scip, cutfeas/cutnorm) )
         {         
            debugMessage("lifted cover       cut for knapsack constraint <%s>: ", SCIPconsGetName(cons));
            debug(SCIProwPrint(row, NULL));
            CHECK_OKAY( SCIPaddCut(scip, row, -cutfeas/cutnorm/(SCIProwGetNNonz(row)+1)) );
            *separated = TRUE;
         }
         CHECK_OKAY( SCIPreleaseRow(scip, &row) );
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
#endif

/** separates lifted cardinality inequalities for given knapsack constraint */
static
RETCODE separateCardinality(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< knapsack constraint */
   Bool*            separated           /**< pointer to store whether a cut was found */
   )
{
   CONSDATA* consdata;
   int* items;
   Longint* weights;
   Real* profits;
   Longint capacity;
   int i;
   int nitems;
   int* fixedones;
   int nfixedones;
   int* fixedzeros;
   int nfixedzeros;
   int* covervars;
   int ncovervars;
   int* noncovervars;
   int nnoncovervars;
   Real infeasibility;
   Real solval;
   Real solvalsum;
   Real transsol;
   Longint coverweight;
   Longint maxweight;

   assert(separated != NULL);

   debugMessage("separating knapsack constraint <%s>\n", SCIPconsGetName(cons));

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
   for( i = 0; i < consdata->nvars; i++ )
   {
      solval = SCIPgetVarSol(scip, consdata->vars[i]);
      solvalsum += solval;
      if( !SCIPisIntegral(scip, solval) )
      {
         items[nitems] = i;
         weights[nitems] = consdata->weights[i];
         profits[nitems] = 1.0 - solval;
         nitems++;
         capacity += consdata->weights[i];
      }
      else if( solval > 0.5 )
      {
         fixedones[nfixedones] = i;
         nfixedones++;
         capacity += consdata->weights[i];
      }
      else
      {
         fixedzeros[nfixedzeros] = i;
         nfixedzeros++;
      }
   }

   if( capacity >= 0 && capacity <= MAX_SEPA_CAPACITY )
   {
      //ROW* row;
      //char name[MAXSTRLEN];
      //Real cutnorm;
      //Real cutfeas;
      Real* solvals;
      int j; 
      int idx;
      int* covervars;
      int ncovervars;
      Real slack;
      Real liftlpval;
      int* liftcoefs;

      /* solve separation knapsack with dynamic programming */
      CHECK_OKAY( SCIPallocBufferArray(scip, &covervars, consdata->nvars) );
      CHECK_OKAY( SCIPallocBufferArray(scip, &noncovervars, consdata->nvars) );
      CHECK_OKAY( solveKnapsack(scip, nitems, weights, profits, capacity, 
                     items, noncovervars, covervars, &nnoncovervars, &ncovervars, &transsol) ); 
      
      CHECK_OKAY( SCIPallocBufferArray(scip, &solvals, consdata->nvars) );
      CHECK_OKAY( SCIPgetVarSols(scip, consdata->nvars, consdata->vars, solvals) );
      CHECK_OKAY( SCIPallocBufferArray(scip, &liftcoefs, consdata->nvars) );
      
      coverweight = 0;
      slack = 0.0;
      for( i = 0; i < ncovervars; i++)
      {
         coverweight += consdata->weights[covervars[i]];
         slack += 1 - solvals[covervars[i]];
      }

      /* add all variables with LP value one to cover vars */
      for( i = 0; i < nfixedones; i++)
      {
         covervars[ncovervars] = fixedones[i];
         coverweight += consdata->weights[fixedones[i]];
         ncovervars++;
      }
      assert(coverweight > consdata->capacity);

      /* add variables with LP value zero to non-cover vars */
      for( i = 0; i < nfixedzeros; i++ )
      {
         noncovervars[nnoncovervars] = fixedzeros[i];
         nnoncovervars++;
      }  

      /* get maximum weight of vars not in cover */
      maxweight = 0;
      for( i = 0; i < nnoncovervars; i++ )
      {
         if( consdata->weights[noncovervars[i]] > maxweight )
            maxweight = consdata->weights[noncovervars[i]];
      }

      /* sort cover vars and noncover vars by nonincreasing LP value */
      for( i = 0; i < ncovervars; i++ ) 
      {
         idx = covervars[i];
         solval = solvals[idx];
         for( j = i - 1; j >= 0 && solvals[covervars[j]] < solval; j-- )
            covervars[j + 1] = covervars[j];
         covervars[j + 1] = idx;
      }
      for( i = 0; i < nnoncovervars; i++ ) 
      {
         idx = noncovervars[i];
         solval = solvals[idx];
         for( j = i - 1; j >= 0 && solvals[noncovervars[j]] < solval; j-- )
            noncovervars[j + 1] = noncovervars[j];
         noncovervars[j + 1] = idx;
      }

      /* generate lifted cardinality inequalities from cover */
      for( i = ncovervars - 1; i >= 0; i-- )
      {
         /* update coverweight and maxweight */
         coverweight -= consdata->weights[covervars[i]];
         if( consdata->weights[covervars[i]] > maxweight )
            maxweight = consdata->weights[covervars[i]];

         /* delete current var from set */
         solval = solvals[covervars[i]];
         for( j = nnoncovervars - 1; j >= 0 && solvals[noncovervars[j]] < solval; j-- )
            noncovervars[j + 1] = noncovervars[j];
         noncovervars[j + 1] = covervars[i];
         slack -= 1 - solval;
         nnoncovervars++;
         ncovervars--;

         if( coverweight <= consdata->capacity - maxweight )
            break;
        
         /* lift variables not in cover into cardinality inequality */
         CHECK_OKAY( liftCardinality(scip, liftcoefs, cons, covervars, noncovervars, ncovervars, 
                        nnoncovervars, ncovervars, &liftlpval) );
         
         /* create LP row */
         if( SCIPisFeasNegative(scip, (slack - liftlpval)/sqrt(ncovervars + 1)) )
         {
            int v;
            ROW* row;
            char name[MAXSTRLEN];
            Real cutnorm;
            Real cutfeas;
            
            sprintf(name, "%s_%lld", SCIPconsGetName(cons), SCIPconshdlrGetNCutsFound(SCIPconsGetHdlr(cons)));
            CHECK_OKAY( SCIPcreateEmptyRow (scip, &row, name, -SCIPinfinity(scip), (Real)ncovervars, 
                           SCIPconsIsLocal(cons), FALSE, SCIPconsIsRemoveable(cons)) );
            
            /* add cover elements to cut */
            for( v = 0; v < ncovervars; v++ )
            {
               CHECK_OKAY( SCIPaddVarToRow(scip, row, consdata->vars[covervars[v]], 1.0) );
            }
            
            /* add lifted non cover elements to cut */
            for( v = 0; v < nnoncovervars; v++ )
            {
               if( liftcoefs[v] > 0 )
               {
                  CHECK_OKAY( SCIPaddVarToRow(scip, row, consdata->vars[noncovervars[v]], (Real)liftcoefs[v]) );
               }
            }
            
            /* check, if cut is violated enough */
            cutnorm = SCIProwGetNorm(row);
            cutfeas = SCIPgetRowLPFeasibility(scip, row);
            if( SCIPisFeasNegative(scip, cutfeas/cutnorm) )
            {         
               debugMessage("lifted cardinality1 cut for knapsack constraint <%s>: ", SCIPconsGetName(cons));
               debug(SCIProwPrint(row, NULL));
               CHECK_OKAY( SCIPaddCut(scip, row, -cutfeas/cutnorm/(SCIProwGetNNonz(row)+1)) );
               *separated = TRUE;
            }
            CHECK_OKAY( SCIPreleaseRow(scip, &row) );
         }
         else
            break;
         
#if 0
         /* ALT: */
         sprintf(name, "%s_%lld", SCIPconsGetName(cons), SCIPconshdlrGetNCutsFound(SCIPconsGetHdlr(cons)));
         CHECK_OKAY( SCIPcreateEmptyRow (scip, &row, name, -SCIPinfinity(scip), (Real)ncovervars, 
                        SCIPconsIsLocal(cons), FALSE, SCIPconsIsRemoveable(cons)) );
         for( i = 0; i < ncovervars; i++ )
         {
            CHECK_OKAY( SCIPaddVarToRow(scip, row, consdata->vars[covervars[i]], 1.0) );
         }
            
         /* lift variables not in set into cardinality inequality */
         CHECK_OKAY( liftCardinality(scip, row, cons, covervars, noncovervars, ncovervars, nnoncovervars, ncovervars) ); 

         /* check, if cut is violated enough */
         cutnorm = SCIProwGetNorm(row);
         cutfeas = SCIPgetRowLPFeasibility(scip, row);
         if( SCIPisFeasNegative(scip, cutfeas/cutnorm) )
         {         
            debugMessage("lifted cardinality cut for knapsack constraint <%s>: ", SCIPconsGetName(cons));
            debug(SCIProwPrint(row, NULL));
            CHECK_OKAY( SCIPaddCut(scip, row, -cutfeas/cutnorm/(SCIProwGetNNonz(row)+1)) );
            *separated = TRUE;
         }
         CHECK_OKAY( SCIPreleaseRow(scip, &row) );
#endif
      }
      CHECK_OKAY( SCIPfreeBufferArray(scip, &liftcoefs) );
      CHECK_OKAY( SCIPfreeBufferArray(scip, &solvals) );
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

/** separates lifted cardinality inequalities for given knapsack constraint */
static
RETCODE separateCardinality2(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< knapsack constraint */
   Bool*            separated           /**< pointer to store whether a cut was found */
   )
{
   CONSDATA* consdata;
   Real* solvals;
   int* perm;
   Longint setweight;
   int* setvars;
   int* nonsetvars;   
   int nsetvars;
   int nnonsetvars;
   int i;
   int j;
   int item;
   Real slack;
   Real liftlpval;
   int* liftcoefs;

   assert(separated != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* set up datastructures */
   CHECK_OKAY( SCIPallocBufferArray(scip, &solvals, consdata->nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &perm, consdata->nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &setvars, consdata->nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &nonsetvars, consdata->nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &liftcoefs, consdata->nvars) );

   sortItems(consdata);
  
   /* sort items by nonincreasing LP solution value */
   CHECK_OKAY( SCIPgetVarSols(scip, consdata->nvars, consdata->vars, solvals) );
   for( i = 0; i < consdata->nvars; ++i )
   {
      for( j = i - 1; j >= 0 && solvals[perm[j]] <= solvals[i]; j-- )
         perm[j + 1] = perm[j];
      perm[j + 1] = i;
   }

   nsetvars = 0;
   nnonsetvars = 0;
   setweight = 0;
   slack = 0.0;

   while( nsetvars < consdata->nvars && setweight <= consdata->capacity && slack < 1.0 )
   {
      item = perm[nsetvars];
      setweight += consdata->weights[item];
      slack += 1.0 - solvals[item];
      setvars[nsetvars] = item;
      nsetvars++;
   }

   for( j = nsetvars; j < consdata->nvars; j++ )
   {
      nonsetvars[nnonsetvars] = perm[j];
      nnonsetvars++;
   }

   while( nsetvars >= 2 )
   {
      int n;

      item = setvars[nsetvars-1];

      /* remove next item from set */
      setweight -= consdata->weights[item];
      slack -= 1.0 - solvals[item];
      nsetvars--;

      /* add item to nonsetvars (as the first element) */
      for( j = nnonsetvars; j > 0; j-- )
         nonsetvars[j] = nonsetvars[j - 1];
      nonsetvars[j] = item;
      nnonsetvars++;

      /* remove all items from nonset that have too small capacity */
      n = 0;
      for( j = 0; j < nnonsetvars; ++j )
      {
         if( setweight + consdata->weights[nonsetvars[j]] > consdata->capacity )
         {
            nonsetvars[n] = nonsetvars[j];
            n++;
         }
      }
      nnonsetvars = n;

      if( nnonsetvars == 0 )
         break;
 
      /* lift variables not in set into inequality */
      CHECK_OKAY( liftCardinality(scip, liftcoefs, cons, setvars, nonsetvars, nsetvars, 
                     nnonsetvars, nsetvars, &liftlpval) );
   
      if( SCIPisFeasNegative(scip, (slack - liftlpval)/sqrt(nsetvars)) )
      {
         int v;
         ROW* row;
         char name[MAXSTRLEN];
         Real cutnorm;
         Real cutfeas;

         /* create LP row */
         sprintf(name, "%s_%lld", SCIPconsGetName(cons), SCIPconshdlrGetNCutsFound(SCIPconsGetHdlr(cons)));
         CHECK_OKAY( SCIPcreateEmptyRow (scip, &row, name, -SCIPinfinity(scip), (Real)nsetvars, 
                        SCIPconsIsLocal(cons), FALSE, SCIPconsIsRemoveable(cons)) );

         /* add set elements to cut */
         for( v = 0; v < nsetvars; v++ )
         {
            CHECK_OKAY( SCIPaddVarToRow(scip, row, consdata->vars[setvars[v]], 1.0) );
         }
         
         /* add lifted non set elements to cut */
         for( v = 0; v < nnonsetvars; v++ )
         {
            if( liftcoefs[v] > 0 )
            {
               CHECK_OKAY( SCIPaddVarToRow(scip, row, consdata->vars[nonsetvars[v]], (Real)liftcoefs[v]) );
            }
         }
                    
         /* check, if cut is violated enough */
         cutnorm = SCIProwGetNorm(row);
         cutfeas = SCIPgetRowLPFeasibility(scip, row);
         if( SCIPisFeasNegative(scip, cutfeas/cutnorm) )
         {         
            debugMessage("lifted cardinality2 cut for knapsack constraint <%s>: ", SCIPconsGetName(cons));
            debug(SCIProwPrint(row, NULL));
            CHECK_OKAY( SCIPaddCut(scip, row, -cutfeas/cutnorm/(SCIProwGetNNonz(row)+1)) );
            *separated = TRUE;
         }
         CHECK_OKAY( SCIPreleaseRow(scip, &row) );
      }
      else
         break;
   }
      
   /* free temporary data */   
   CHECK_OKAY( SCIPfreeBufferArray(scip, &liftcoefs) );
   CHECK_OKAY( SCIPfreeBufferArray(scip, &nonsetvars) );
   CHECK_OKAY( SCIPfreeBufferArray(scip, &setvars) );
   CHECK_OKAY( SCIPfreeBufferArray(scip, &perm) );
   CHECK_OKAY( SCIPfreeBufferArray(scip, &solvals) );
  
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
   Longint zerosweightsum;
   int i;

   assert(cutoff != NULL);
   assert(redundant != NULL);
   assert(nfixedvars != NULL);

   debugMessage("propagating knapsack constraint <%s>\n", SCIPconsGetName(cons));

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
   zerosweightsum = 0;
   for( i = consdata->nvars - 1; i >= 0; i-- )
   {
      if( consdata->weights[i] > consdata->capacity - consdata->onesweightsum )
      {
         if( SCIPvarGetLbLocal(consdata->vars[i]) < 0.5 )
         {
            if( SCIPvarGetUbLocal(consdata->vars[i]) > 0.5 )
            {
               CHECK_OKAY( SCIPinferBinVar(scip, consdata->vars[i], FALSE, cons, 0, &infeasible, &tightened) );
               assert(!infeasible);
               assert(tightened);
               (*nfixedvars)++;
            }
            zerosweightsum += consdata->weights[i];
         }
      }
      else
         break;
   }

   /* if the remaining (potentially unfixed) variables would fit all into the knapsack, the knapsack is now redundant */
   if( consdata->weightsum - zerosweightsum <= consdata->capacity )
   {
      debugMessage("knapsack constraint <%s> is redundant: weightsum=%lld, zerosweightsum=%lld, capacity=%lld\n",
         SCIPconsGetName(cons), consdata->weightsum, zerosweightsum, consdata->capacity);
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
   assert(consdata->weightsum >= 0);
   assert(consdata->onesweightsum >= 0);

   /* move the last variable to the free slot */
   consdata->vars[pos] = consdata->vars[consdata->nvars-1];
   consdata->weights[pos] = consdata->weights[consdata->nvars-1];
   consdata->eventdatas[pos] = consdata->eventdatas[consdata->nvars-1];
   consdata->nvars--;

   consdata->propagated = FALSE;
   consdata->sorted = FALSE;

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
   assert(consdata->onesweightsum == 0);

   return SCIP_OKAY;
}

/** divides weights by their greatest common divisor and divides capacity by the same value, rounding down the result */
static
void normalizeWeights(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   int*             nchgcoefs,          /**< pointer to count total number of changed coefficients */
   int*             nchgsides           /**< pointer to count number of side changes */
   )
{
   CONSDATA* consdata;
   Longint gcd;
   int i;

   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(!SCIPconsIsModifiable(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL); /* we are in presolve, so no LP row exists */
   assert(consdata->onesweightsum == 0); /* all fixed variables should have been removed */
   assert(consdata->weightsum > consdata->capacity); /* otherwise, the constraint is redundant */
   assert(consdata->nvars >= 1);

   /* sort items, because we can stop earlier if the smaller weights are evaluated first */
   sortItems(consdata);

   gcd = (Longint)consdata->weights[0];
   for( i = 1; i < consdata->nvars && gcd >= 2; ++i )
   {
      assert(SCIPvarGetLbLocal(consdata->vars[i]) < 0.5);
      assert(SCIPvarGetUbLocal(consdata->vars[i]) > 0.5); /* all fixed variables should have been removed */

      gcd = SCIPcalcGreComDiv(gcd, (Longint)consdata->weights[i]);
   }

   if( gcd >= 2 )
   {
      debugMessage("knapsack constraint <%s>: dividing weights by %lld\n", SCIPconsGetName(cons), gcd);

      for( i = 0; i < consdata->nvars; ++i )
         consdataChgWeight(consdata, i, consdata->weights[i]/gcd);
      consdata->capacity /= gcd;
      (*nchgcoefs) += consdata->nvars;
      (*nchgsides)++;
   }
}

/** tightens item weights and capacity in presolving:
 *  - given a knapsack  w*x + wi*xi <= capacity
 *  - let weightsum := sum{w} + wi
 *  (1) if  weightsum - wi < capacity:  (this can only apply for the heaviest item)
 *      - not using item i would make the knapsack constraint redundant
 *      - wi and capacity can be changed to have the same redundancy effect and the same results for
 *        fixing xi to zero or one, but with a reduced wi and tightened capacity to tighten the LP relaxation
 *      - change coefficients:
 *          wi'       := weightsum - capacity
 *          capacity' := capacity - (wi - wi')
 *  (2) if  min{w} + wi > capacity:
 *      - using item i would force fix other items to zero
 *      - wi can increased to the capacity
 */
static
void tightenWeights(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   int*             nchgcoefs,          /**< pointer to count total number of changed coefficients */
   int*             nchgsides           /**< pointer to count number of side changes */
   )
{
   CONSDATA* consdata;
   Longint weight;
   Longint newweight;
   Longint minweight;
   int i;

   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(!SCIPconsIsModifiable(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL); /* we are in presolve, so no LP row exists */
   assert(consdata->onesweightsum == 0); /* all fixed variables should have been removed */
   assert(consdata->weightsum > consdata->capacity); /* otherwise, the constraint is redundant */
   assert(consdata->nvars > 0);

   /* sort items, s.t. the heaviest one is in the last position */
   sortItems(consdata);

   /* apply rule (1) */
   weight = consdata->weights[consdata->nvars-1];
   if( consdata->weightsum - weight < consdata->capacity )
   {
      newweight = consdata->weightsum - consdata->capacity;
      consdataChgWeight(consdata, consdata->nvars-1, newweight);
      consdata->capacity -= (weight - newweight);
      (*nchgcoefs)++;
      (*nchgsides)++;
      debugMessage("knapsack constraint <%s>: changed weight of <%s> from %lld to %lld, capacity from %lld to %lld\n",
         SCIPconsGetName(cons), SCIPvarGetName(consdata->vars[consdata->nvars-1]),
         weight, newweight, consdata->capacity + (weight-newweight), consdata->capacity);
   }

   /* apply rule (2) */
   minweight = consdata->weights[0];
   for( i = consdata->nvars-1; i >= 0; --i )
   {
      weight = consdata->weights[i];
      if( minweight + weight > consdata->capacity && weight < consdata->capacity )
      {
         debugMessage("knapsack constraint <%s>: changing weight of <%s> from %lld to %lld\n",
            SCIPconsGetName(cons), SCIPvarGetName(consdata->vars[i]), weight, consdata->capacity);
         consdataChgWeight(consdata, i, consdata->capacity);
         (*nchgcoefs)++;
      }
      else
         break;
   }
}




/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
DECL_CONSFREE(consFreeKnapsack)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeMemory(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


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
   CONSHDLRDATA* conshdlrdata;
   Bool separated;
   int depth;
   int nrounds;
   int maxsepacuts;
   int ncuts;
   int i;

   *result = SCIP_DIDNOTRUN;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   depth = SCIPgetDepth(scip);
   nrounds = SCIPgetNSepaRounds(scip);

   debugMessage("knapsack separation of %d/%d constraints, round %d (max %d/%d)\n",
      nusefulconss, nconss, nrounds, conshdlrdata->maxroundsroot, conshdlrdata->maxrounds);

   /* only call the separator a given number of times at each node */
   if( (depth == 0 && nrounds >= conshdlrdata->maxroundsroot)
      || (depth > 0 && nrounds >= conshdlrdata->maxrounds) )
      return SCIP_OKAY;

   /* get the maximal number of cuts allowed in a separation round */
   maxsepacuts = depth == 0 ? conshdlrdata->maxsepacutsroot : conshdlrdata->maxsepacuts;

   *result = SCIP_DIDNOTFIND;

   ncuts = 0;
   for( i = 0; i < nusefulconss && ncuts < maxsepacuts; i++ )
   {
#if 0
      CHECK_OKAY( separateCovers(scip, conss[i], &separated) );
      if( separated )
      {
         ncuts++;
         *result = SCIP_SEPARATED;
      }
#endif

      CHECK_OKAY( separateCardinality(scip, conss[i], &separated) );
      if( separated )
      {
         ncuts++;
         *result = SCIP_SEPARATED;
      }

#if 0
      CHECK_OKAY( separateCardinality2(scip, conss[i], &separated) );
      if( separated )
      {
         ncuts++;
         *result = SCIP_SEPARATED;
      }
#endif
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
   CONSDATA* consdata;
   Bool cutoff;
   Bool redundant;
   int oldnfixedvars;
   int oldndelconss;
   int oldnchgcoefs;
   int oldnchgsides;
   int i;

   cutoff = FALSE;
   oldnfixedvars = *nfixedvars;
   oldndelconss = *ndelconss;
   oldnchgcoefs = *nchgcoefs;
   oldnchgsides = *nchgsides;

   /* get event handler for bound change events */
   eventhdlr = SCIPfindEventhdlr (scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);
   
   for( i = 0; i < nconss && !cutoff; i++ )
   {
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      /* remove all variables that are fixed to zero, check redundancy due to fixed-to-one variable */
      CHECK_OKAY( applyFixings(scip, conss[i], eventhdlr) );

      if( consdata->propagated )
         continue;

      /* propagate constraint */
      CHECK_OKAY( propagateCons(scip, conss[i], &cutoff, &redundant, nfixedvars) );
      if( redundant )
      {
         (*ndelconss)++;
         continue;
      }

      if( !SCIPconsIsModifiable(conss[i]) )
      {
         /* divide weights by their greatest common divisor */
         normalizeWeights(scip, conss[i], nchgcoefs, nchgsides);

         /* tighten capacity and weights */
         tightenWeights(scip, conss[i], nchgcoefs, nchgsides);
      }
   } 

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( *nfixedvars > oldnfixedvars || *ndelconss > oldndelconss
      || *nchgcoefs > oldnchgcoefs || *nchgsides > oldnchgsides )
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

   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
 
   assert(SCIPvarGetUbLocal(infervar) < 0.5);

   for( i = 0; i < consdata->nvars; i++ )
   {
      if( SCIPvarGetLbLocal(consdata->vars[i]) > 0.5 && SCIPvarWasFixedEarlier(consdata->vars[i], infervar) )
      {
         CHECK_OKAY( SCIPaddConflictVar(scip, consdata->vars[i]) );
      }
   }
   *result = SCIP_SUCCESS;

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
   Longint* weights;
   Longint capacity;
   Longint weight;
   int mult;
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);
   assert(SCIPisInfinity(scip, -lhs) != SCIPisInfinity(scip, rhs));

   /* get temporary memory */
   CHECK_OKAY( SCIPallocBufferArray(scip, &transvars, nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &weights, nvars) );

   /* if the right hand side is non-infinite, we have to negate all variables with negative coefficient;
    * otherwise, we have to negate all variables with positive coefficient and multiply the row with -1
    */
   if( SCIPisInfinity(scip, rhs) )
   {
      mult = -1;
      capacity = (Longint)SCIPfloor(scip, -lhs);
   }
   else
   {
      mult = +1;
      capacity = (Longint)SCIPfloor(scip, rhs);
   }

   /* negate positive or negative variables */
   for( v = 0; v < nvars; ++v )
   {
      weight = mult * vals[v];
      if( weight > 0 )
      {
         transvars[v] = vars[v];
         weights[v] = weight;
      }
      else
      {
         CHECK_OKAY( SCIPgetNegatedVar(scip, vars[v], &transvars[v]) );
         weights[v] = -weight;
         capacity -= weight;
      }
      assert(transvars[v] != NULL);
   }

   /* create the constraint */
   CHECK_OKAY( SCIPcreateConsKnapsack(scip, cons, name, nvars, transvars, weights, capacity,
                  initial, separate, enforce, check, propagate, local, modifiable, removeable) );

   /* free temporary memory */
   CHECK_OKAY( SCIPfreeBufferArray(scip, &weights) );
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
   CHECK_OKAY( SCIPallocMemory(scip, &conshdlrdata) );

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
   CHECK_OKAY( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC, 
                  NULL, NULL, NULL, NULL, eventExecKnapsack,
                  eventhdlrdata) );

#ifdef LINCONSUPGD_PRIORITY
   /* include the linear constraint upgrade in the linear constraint handler */
   CHECK_OKAY( SCIPincludeLinconsUpgrade(scip, linconsUpgdKnapsack, LINCONSUPGD_PRIORITY) );
#endif

   /* add knapsack constraint handler parameters */
   CHECK_OKAY( SCIPaddIntParam(scip,
                  "constraints/knapsack/maxrounds",
                  "maximal number of separation rounds per node",
                  &conshdlrdata->maxrounds, DEFAULT_MAXROUNDS, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
                  "constraints/knapsack/maxroundsroot",
                  "maximal number of separation rounds per node in the root node",
                  &conshdlrdata->maxroundsroot, DEFAULT_MAXROUNDSROOT, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
                  "constraints/knapsack/maxsepacuts",
                  "maximal number of cuts separated per separation round",
                  &conshdlrdata->maxsepacuts, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
                  "constraints/knapsack/maxsepacutsroot",
                  "maximal number of cuts separated per separation round in the root node",
                  &conshdlrdata->maxsepacutsroot, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a knapsack constraint */
RETCODE SCIPcreateConsKnapsack(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nvars,              /**< number of items in the knapsack */
   VAR**            vars,               /**< array with item variables */
   Longint*         weights,            /**< array with item weights */
   Longint          capacity,           /**< capacity of knapsack */
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
   CHECK_OKAY( consdataCreate(scip, &consdata, nvars, vars, weights, capacity) );
        
   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
                  local, modifiable, removeable) );

   return SCIP_OKAY;
}

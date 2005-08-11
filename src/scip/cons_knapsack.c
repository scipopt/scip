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
#pragma ident "@(#) $Id: cons_knapsack.c,v 1.103 2005/08/11 07:53:59 bzfpfend Exp $"

/**@file   cons_knapsack.c
 * @brief  constraint handler for knapsack constraints
 * @author Tobias Achterberg
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "knapsack"
#define CONSHDLR_DESC          "knapsack constraint of the form  a^T x <= b, x binary"
#define CONSHDLR_SEPAPRIORITY   +600000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   +600000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -600000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define EVENTHDLR_NAME         "knapsack"
#define EVENTHDLR_DESC         "bound change event handler for knapsack constraints"

#define LINCONSUPGD_PRIORITY    +100000

#define MAX_DYNPROG_CAPACITY      10000 /**< maximal capacity of knapsack to apply dynamic programming */

#define DEFAULT_SEPACARDFREQ         10 /**< multiplier on separation frequency, how often cardinality cuts are separated */
#define DEFAULT_MAXROUNDS             5 /**< maximal number of separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        -1 /**< maximal number of separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS          50 /**< maximal number of cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT     200 /**< maximal number of cuts separated per separation round in the root node */
#define DEFAULT_MAXNUMCARDLIFT       -1 /**< maximal number of cardinality inequalities lifted per sepa round (-1: unlimited) */




/*
 * Data structures
 */

/** constraint handler data */
struct ConshdlrData
{
   EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
   int              sepacardfreq;       /**< multiplier on separation frequency, how often cardinality cuts are separated */
   int              maxrounds;          /**< maximal number of separation rounds per node (-1: unlimited) */
   int              maxroundsroot;      /**< maximal number of separation rounds in the root node (-1: unlimited) */
   int              maxsepacuts;        /**< maximal number of cuts separated per separation round */
   int              maxsepacutsroot;    /**< maximal number of cuts separated per separation round in the root node */
   int              maxnumcardlift;     /**< maximal number of cardinality inequ. lifted per sepa round (-1: unlimited) */
};


/** constraint data for knapsack constraints */
struct ConsData
{
   VAR**            vars;               /**< variables in knapsack constraint */
   Longint*         weights;            /**< weights of variables in knapsack constraint */
   EVENTDATA**      eventdatas;         /**< event datas for bound change events of the variables */
   ROW*             row;                /**< corresponding LP row */
   int              nvars;              /**< number of variables in knapsack constraint */
   int              varssize;           /**< size of vars, weights, and eventdatas arrays */
   Longint          capacity;           /**< capacity of knapsack */
   Longint          weightsum;          /**< sum of all weights */
   Longint          onesweightsum;      /**< sum of weights of variables fixed to one */
   unsigned int     propagated:1;       /**< is the knapsack constraint already propagated? */
   unsigned int     presolved:1;        /**< is the knapsack constraint already presolved? */
   unsigned int     sorted:1;           /**< are the knapsack items sorted by weight? */
   unsigned int     merged:1;           /**< are the constraint's equal variables already merged? */
};


/** event data for bound changes events */
struct EventData
{
   CONSDATA*        consdata;           /**< knapsack constraint data to process the bound change for */
   Longint          weight;             /**< weight of variable */
};




/*
 * Local methods
 */

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

/** installs rounding locks for the given variable in the given knapsack constraint */
static
RETCODE lockRounding(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< knapsack constraint */
   VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding up may violate the constraint */
   CHECK_OKAY( SCIPlockVarCons(scip, var, cons, FALSE, TRUE) );

   return SCIP_OKAY;
}

/** removes rounding locks for the given variable in the given knapsack constraint */
static
RETCODE unlockRounding(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< knapsack constraint */
   VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding up may violate the constraint */
   CHECK_OKAY( SCIPunlockVarCons(scip, var, cons, FALSE, TRUE) );

   return SCIP_OKAY;
}

/** catches bound change events for variables in knapsack */
static
RETCODE catchEvents(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< constraint data */
   EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;
   
   assert(consdata != NULL);

   for( i = 0; i < consdata->nvars; i++)
   {
      CHECK_OKAY( eventdataCreate(scip, &consdata->eventdatas[i], consdata, consdata->weights[i]) );
      CHECK_OKAY( SCIPcatchVarEvent(scip, consdata->vars[i], 
            SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_UBRELAXED | SCIP_EVENTTYPE_VARFIXED,
            eventhdlr, consdata->eventdatas[i], NULL) );
   }

   return SCIP_OKAY;
}

/** drops bound change events for variables in knapsack */
static
RETCODE dropEvents(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< constraint data */
   EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;
   
   assert(consdata != NULL);

   for( i = 0; i < consdata->nvars; i++)
   {
      CHECK_OKAY( SCIPdropVarEvent(scip, consdata->vars[i],
            SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_UBRELAXED | SCIP_EVENTTYPE_VARFIXED,
            eventhdlr, consdata->eventdatas[i], -1) );
      CHECK_OKAY( eventdataFree(scip, &consdata->eventdatas[i]) );
   }

   return SCIP_OKAY;
}

/** ensures, that vars and vals arrays can store at least num entries */
static
RETCODE consdataEnsureVarsSize(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint data */
   int              num,                /**< minimum number of entries to store */
   Bool             transformed         /**< is constraint from transformed problem? */
   )
{
   assert(consdata != NULL);
   assert(consdata->nvars <= consdata->varssize);
   
   if( num > consdata->varssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      CHECK_OKAY( SCIPreallocBlockMemoryArray(scip, &consdata->vars, consdata->varssize, newsize) );
      CHECK_OKAY( SCIPreallocBlockMemoryArray(scip, &consdata->weights, consdata->varssize, newsize) );
      if( transformed )
      {
         CHECK_OKAY( SCIPreallocBlockMemoryArray(scip, &consdata->eventdatas, consdata->varssize, newsize) );
      }
      else
         assert(consdata->eventdatas == NULL);
      consdata->varssize = newsize;
   }
   assert(num <= consdata->varssize);

   return SCIP_OKAY;
}

/** creates knapsack constraint data */
static
RETCODE consdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to store constraint data */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              nvars,              /**< number of variables in knapsack */
   VAR**            vars,               /**< variables of knapsack */
   Longint*         weights,            /**< weights of knapsack items */
   Longint          capacity            /**< capacity of knapsack */
   )
{
   int i;

   assert(consdata != NULL);

   CHECK_OKAY( SCIPallocBlockMemory(scip, consdata) );
   if( nvars > 0 )
   {
      CHECK_OKAY( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
      CHECK_OKAY( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->weights, weights, nvars) );
   }
   else
   {
      (*consdata)->vars = NULL;
      (*consdata)->weights = NULL;
   }
   (*consdata)->nvars = nvars;
   (*consdata)->varssize = nvars;
   (*consdata)->capacity = capacity;
   (*consdata)->row = NULL;
   (*consdata)->eventdatas = NULL;
   (*consdata)->weightsum = 0;
   (*consdata)->onesweightsum = 0;
   (*consdata)->propagated = FALSE;
   (*consdata)->presolved = FALSE;
   (*consdata)->sorted = FALSE;
   (*consdata)->merged = FALSE;

   /* get transformed variables, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      CHECK_OKAY( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );

      /* catch events for variables */
      CHECK_OKAY( SCIPallocBlockMemoryArray(scip, &(*consdata)->eventdatas, nvars) );
      CHECK_OKAY( catchEvents(scip, *consdata, eventhdlr) );
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
   CONSDATA**       consdata,           /**< pointer to the constraint data */
   EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
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
      CHECK_OKAY( dropEvents(scip, *consdata, eventhdlr) );
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
   consdata->presolved = FALSE;
   consdata->sorted = FALSE;
}

/** creates LP row corresponding to knapsack constraint */
static 
RETCODE createRelaxation(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< knapsack constraint */
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

   CHECK_OKAY( SCIPcacheRowExtensions(scip, consdata->row) );
   for( i = 0; i < consdata->nvars; ++i )
   {
      CHECK_OKAY( SCIPaddVarToRow(scip, consdata->row, consdata->vars[i], (Real)consdata->weights[i]) );
   }
   CHECK_OKAY( SCIPflushRowExtensions(scip, consdata->row) );

   return SCIP_OKAY;
}  

/** adds linear relaxation of knapsack constraint to the LP */
static 
RETCODE addRelaxation(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< knapsack constraint */
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
   assert(!SCIProwIsInLP(consdata->row));

   debugMessage("adding relaxation of knapsack constraint <%s> (capacity %lld): ", 
      SCIPconsGetName(cons), consdata->capacity);
   debug( SCIProwPrint(consdata->row, NULL) );
   CHECK_OKAY( SCIPaddCut(scip, consdata->row, FALSE) );

   return SCIP_OKAY;
}

/** checks knapsack constraint for feasibility of given solution: returns TRUE iff constraint is feasible */
static
RETCODE checkCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint to check */
   SOL*             sol,                /**< solution to check, NULL for current solution */
   Bool             checklprows,        /**< should LP rows be checked? */
   Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   CONSDATA* consdata;

   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   debugMessage("checking knapsack constraint <%s> for feasibility of solution %p (lprows=%d)\n",
      SCIPconsGetName(cons), sol, checklprows);

   *violated = FALSE;

   if( checklprows || consdata->row == NULL || !SCIProwIsInLP(consdata->row) )
   {
      Real sum;
      int i;

      /* increase age of constraint; age is reset to zero, if a violation was found */
      CHECK_OKAY( SCIPincConsAge(scip, cons) );

      sum = 0.0;
      for( i = 0; i < consdata->nvars && sum <= consdata->capacity + 0.1; i++ )
      {
         sum += consdata->weights[i] * SCIPgetSolVal(scip, sol, consdata->vars[i]);
      }

      if( SCIPisFeasGT(scip, sum, (Real)consdata->capacity) )
      {
         CHECK_OKAY( SCIPresetConsAge(scip, cons) );
         *violated = TRUE;
      }
   }

   return SCIP_OKAY;
}

#define IDX(j,d) ((j)*(intcap+1)+(d))

/** solves knapsack problem with dynamic programming;
 *  if needed, one can provide arrays to store all selected items and all not selected items
 */
RETCODE SCIPsolveKnapsack(
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
   Real*            solval              /**< pointer to store optimal solution value, or NULL */
   ) 
{
   Real* optvalues;
   int intcap;
   int d;
   int j;

   assert(weights != NULL);
   assert(profits != NULL);
   assert(capacity >= 0);
   assert(capacity < INT_MAX);
   assert(nitems >= 0);

   intcap = (int)capacity;
   assert(intcap >= 0);
   CHECK_OKAY( SCIPallocBufferArray(scip, &optvalues, (nitems+1)*(intcap+1)) );
   
   /* fill dynamic programming table with optimal values */
   for( d = 0; d <= intcap; d++ )
      optvalues[IDX(0,d)] = 0.0;
   for( j = 1; j <= nitems; j++ )
   {
      int intweight;

      assert(0 <= weights[j-1] && weights[j-1] < INT_MAX);
      intweight = (int)weights[j-1];
      assert(intweight >= 0);

      for( d = 0; d < intweight && d <= intcap; d++ )
         optvalues[IDX(j,d)] = optvalues[IDX(j-1,d)];
      for( d = intweight; d <= intcap; d++ )
      {
         Real sumprofit;

         sumprofit = optvalues[IDX(j-1,d-intweight)] + profits[j-1];
         optvalues[IDX(j,d)] = MAX(sumprofit, optvalues[IDX(j-1,d)]);
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
      d = intcap;
      
      for( j = nitems; j > 0; j-- )
      {
         if( optvalues[IDX(j,d)] > optvalues[IDX(j-1,d)] )
         {
            assert(0 <= weights[j-1] && weights[j-1] < INT_MAX);
            solitems[*nsolitems] = items[j-1];
            (*nsolitems)++;
            d -= (int)weights[j-1];
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

   if( solval != NULL )
      *solval = optvalues[IDX(nitems,intcap)];

   SCIPfreeBufferArray(scip, &optvalues);

   return SCIP_OKAY;
}

/** gets a most violated minimal cover C = C2 & C1 for a given knapsack constraint, considering that a all variables 
 *  in a given set C2 of variables in knapsack constraint are fixed to one (variables for downlifting)
 */
static
RETCODE getCover(
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            vars,               /**< variables in knapsack constraint */
   int              nvars,              /**< number of variables in knapsack constraint */
   Longint*         weights,            /**< weights of variables in knapsack constraint */
   Longint          capacity,           /**< capacity of knapsack */
   Real*            solvals,            /**< LP values of all problem variables */
   int*             varsc2,             /**< set C2 of variables in knapsack constraint (variables for downlifting) */
   int*             varscn2,            /**< set of variables in knapsack constraint that are not in C2 (cand. for C1) */
   int              nvarsc2,            /**< number of variables in C2 */
   int              nvarscn2,           /**< number of variables not in C2 */
   Longint          varsc2weight,       /**< sum of weights of variables in C2 */
   int*             covervars,          /**< pointer to store cover variables C = C2 & C1 */
   int*             noncovervars,       /**< pointer to store noncover variables */
   int*             ncovervars,         /**< pointer to store number of cover variables */
   int*             ncovervarsc1,       /**< pointer to store number of cover variables in C1 (at the end of covervars) */
   int*             ncovervarsc2,       /**< pointer to store number of cover variables in C2 (at the beg of covervars) */
   int*             nnoncovervars,      /**< pointer to store number of noncover variables */
   Longint*         coverweight,        /**< pointer to store weight of cover */
   Real*            covervarsc1activity,/**< pointer to store activity of cover variables in C1 */  
   Bool*            found               /**< pointer to store whether a cover was found */
   )
{
   Real* dpprofits;
   Real relslack;
   Longint* dpweights;
   Longint dpcapacity;
   int* items;
   int* fixedzeros;
   int* setvars;
   int* nonsetvars;
   int nitems;
   int nfixedzeros;
   int nsetvars;
   int nnonsetvars;
   int i;
   int j; 

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(weights != NULL);
   assert(capacity >= 0);
   assert(solvals != NULL);
   assert(varsc2 != NULL);
   assert(varscn2 != NULL);
   assert(nvarsc2 >= 0);
   assert(nvarscn2 >= 0);
   assert(nvarsc2 + nvarscn2 == nvars);
   assert(varsc2weight >= 0);
   assert(covervars != NULL);
   assert(noncovervars != NULL);
   assert(ncovervars != NULL);
   assert(ncovervarsc1 != NULL);
   assert(ncovervarsc2 != NULL);
   assert(nnoncovervars != NULL);
   assert(coverweight != NULL);
   assert(covervarsc1activity != NULL);
   assert(found != NULL);
   
   /* allocates temporary memory */
   CHECK_OKAY( SCIPallocBufferArray(scip, &items, nvarscn2) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &dpweights, nvarscn2) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &dpprofits, nvarscn2) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &fixedzeros, nvarscn2) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &setvars, nvarscn2) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &nonsetvars, nvarscn2) );

   *found = FALSE;
   *ncovervars = 0;
   *ncovervarsc1 = 0;
   *ncovervarsc2 = 0;
   *nnoncovervars = 0;
   *coverweight = 0;
   *covervarsc1activity = 0.0;
      
   /* solves the following knapsack problem with dynamic programming or greedy heuristic:
    * max SUM (1 - x_i^*) z_i
    *     SUM w_i z_i <= SUM w_i - capacity - 1
    *             z_i == 0                       forall variables x_i in C2
    * with the meaning: x_i is member of cover <==> z_i == 0
    */
   
   /* uses all variables with fractional LP value from the set of variables not in C2 */
   nitems = 0;
   nfixedzeros = 0;
   dpcapacity = varsc2weight - capacity - 1;
   for( i = 0; i < nvarscn2; i++ )
   {
      assert(SCIPvarGetType(vars[varscn2[i]]) == SCIP_VARTYPE_BINARY);
      
      /* variable has fractional LP value */
      if( SCIPisFeasGT(scip, solvals[varscn2[i]], 0.0) )
      {
         /* sorts items by non-decreasing relative slack value (1-x_i^*)/w_i */
         relslack = (1.0 - solvals[varscn2[i]])/weights[varscn2[i]];
         for( j = nitems; j > 0 && relslack < (1.0 - solvals[varscn2[j-1]])/weights[varscn2[j-1]]; --j )
            items[j] = items[j-1];
         items[j] = varscn2[i];
         nitems++;
         dpcapacity += weights[varscn2[i]];
      }
      /* variable has LP value 0 */
      else
      {
         fixedzeros[nfixedzeros] = varscn2[i];
         nfixedzeros++;
      }
   }
   assert(nitems + nfixedzeros == nvarscn2);

   /* gets weights and slacks of items */
   for( i = 0; i < nitems; ++i )
   {
      dpweights[i] = weights[items[i]];
      dpprofits[i] = 1.0 - solvals[items[i]];
   }
   
   /* if the right hand side of the separation knapsack is negative, there cannot be a cover 
    * (when variables from fixedzeros are not used)
    */
   if( dpcapacity >= 0 )
   {
      /* if the right hand side is too large, solving the dynamic program would be too expensive in space and time */
      if( dpcapacity <= MAX_DYNPROG_CAPACITY )
      {
         /* solves separation knapsack with dynamic programming */
         CHECK_OKAY( SCIPsolveKnapsack(scip, nitems, dpweights, dpprofits, dpcapacity, 
               items, nonsetvars, setvars, &nnonsetvars, &nsetvars, NULL) ); 
         assert(nsetvars + nnonsetvars == nitems);
      }
      else
      {
         Longint setweight;
         Longint minimumweight;
         int minimumweightidx;
         int removevar;
         Bool minimal;

         /* uses greedy heuristic on separation knapsack */
         nnonsetvars = 0;
         nsetvars = 0;
         for( i = nitems-1; i >= 0 && dpcapacity - dpweights[i] >= 0; --i )
         {
            dpcapacity -= dpweights[i];
            nonsetvars[nnonsetvars] = items[i];
            (nnonsetvars)++;
         }
         for( ; i >= 0; --i )
         {
            setvars[nsetvars] = items[i];
            (nsetvars)++;
         }
         
         /* calculates weight of set */
         setweight = 0;
         for( i = 0; i < nsetvars; i++)
            setweight += weights[setvars[i]];

         /* gets minimum weigth of set variables */
         minimal = FALSE;
         minimumweight = LONGINT_MAX;
         minimumweightidx = -1;
         for( i = 0; i < nsetvars; i++ )
         {
            if( weights[setvars[i]] <= minimumweight )
            {
               minimumweightidx = i; 
               minimumweight = weights[setvars[i]];
            }
         }
         assert(minimumweightidx >= 0 && minimumweightidx < nsetvars);
         assert(minimumweight == weights[setvars[minimumweightidx]]);
         assert(setweight > capacity - varsc2weight);

         /* tests if set is allready minimal */
         if( setweight - minimumweight <= capacity - varsc2weight )
            minimal = TRUE;

         /* makes set minimal by removing variables (in decreasing order of slack) */
         for( i = 0; i < nsetvars && !minimal; i++ )
         {
            /* set variable x_i can not be removed with respect to the set property */
            if( setweight - weights[setvars[i]] <= capacity - varsc2weight )
               continue;

            /* removes x_i from setvars */
            removevar = setvars[i];
            for( j = i; j < nsetvars - 1; j++ )
               setvars[j] = setvars[j + 1];
            nsetvars--;
            setweight -= weights[removevar];
            
            /* adds x_i to nonsetvars */
            nonsetvars[nnonsetvars] = removevar;
            nnonsetvars++;
            
            /* updates minimumweight of setvars */     
            if( minimumweightidx == i )
            {
               /* gets minimum weigth of setvars */
               minimumweight = LONGINT_MAX;
               minimumweightidx = -1;
               for( j = 0; j < nsetvars; j++ )
               {
                  if( weights[setvars[j]] <= minimumweight )
                  {
                     minimumweightidx = j; 
                     minimumweight = weights[setvars[j]];
                  }
               }
            }
            else
               minimumweightidx--;
            assert(minimumweightidx >= 0 && minimumweightidx < nsetvars);
            assert(minimumweight == weights[setvars[minimumweightidx]]);
         
            /* updates index of current setvar */
            i--;

            /* tests if set is now minimal */
            assert(setweight > capacity - varsc2weight);
            assert(!minimal);
            if( (setweight) - minimumweight <= capacity - varsc2weight )
               minimal = TRUE;
         }
         assert(minimal);
      }
      assert(*ncovervars == 0);
      assert(*nnoncovervars == 0);
      assert(nvarsc2 + nsetvars + nnonsetvars + nfixedzeros == nvars);

      /* gets covervars (variables from C2 are in the beginning of covervars) and noncovervars: */ 

      /* adds all variables from C2 to covervars */
      for( i = 0; i < nvarsc2; i++)
      {
         covervars[*ncovervars] = varsc2[i];
         *coverweight += weights[varsc2[i]];
         (*ncovervars)++;
         (*ncovervarsc2)++;
      }

      /* adds all variables from setvars (variables from set of variables not in C2 which have been chosen to be in the 
       * cover) to covervars 
       */
      for( i = 0; i < nsetvars; i++)
      {
         covervars[*ncovervars] = setvars[i];
         *coverweight += weights[setvars[i]];
         *covervarsc1activity += solvals[setvars[i]];
         (*ncovervars)++;
         (*ncovervarsc1)++;
      }
      assert(*coverweight > capacity);
      assert((*ncovervarsc1) + (*ncovervarsc2) == *ncovervars);

      /* adds all variables from nonsetvars (variables from set of variables not in C2 which have not been chosen to be 
       * in the cover) to noncovervars 
       */ 
      for( i = 0; i < nnonsetvars; i++ )
      {
         noncovervars[*nnoncovervars] = nonsetvars[i];
         (*nnoncovervars)++;
      }  

      /* adds all variables with LP value zero to noncovervars */
      for( i = 0; i < nfixedzeros; i++ )
      {
         noncovervars[*nnoncovervars] = fixedzeros[i];
         (*nnoncovervars)++;
      }  
      assert((*ncovervarsc2) + (*ncovervarsc1) == (*ncovervars));
      assert((*ncovervars) + (*nnoncovervars) == nvars);
      *found = TRUE;
   }
   else
   {
      /* puts all variables into noncovervars because no cover was found (variables from fixedzeros have not been used) */
      
      /* adds all variables from C2 to noncovervars */
      for( i = 0; i < nvarsc2; i++)
      {
         noncovervars[*nnoncovervars] = varsc2[i];
         (*nnoncovervars)++;
      }
      
      /* adds all variables from set of variables not in C2 to noncovervars */
      for( i = 0; i < nvarscn2; i++)
      {
         noncovervars[*nnoncovervars] = varscn2[i];
         (*nnoncovervars)++;
      }
      assert(!(*found));
      assert(*coverweight == 0 && *covervarsc1activity == 0.0); 
   }

   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &nonsetvars);
   SCIPfreeBufferArray(scip, &setvars);
   SCIPfreeBufferArray(scip, &fixedzeros);
   SCIPfreeBufferArray(scip, &dpprofits);
   SCIPfreeBufferArray(scip, &dpweights);
   SCIPfreeBufferArray(scip, &items);

   return SCIP_OKAY;
}

/** gets a partition of all variables in a given knapsack constraint into a set of variables for downlifting (C2) and 
 *  a set of all remaining variables    
 */
static
void getPartition(
   SCIP*            scip,               /**< SCIP data structure */
   int              nvars,              /**< number of variables in knapsack constraint */
   Longint*         weights,            /**< weights of variables in knapsack constraint */
   Real*            solvals,            /**< LP values of all problem variables */
   int*             varsc2,             /**< pointer to store variables of knapsack constraint in C2 */
   int*             varscn2,            /**< pointer to store variables of knapsack constraint not in C1 */
   int*             nvarsc2,            /**< pointer to store number of variables of knapsack constraint in C2 */
   int*             nvarscn2,           /**< pointer to store number of variables of knapsack constraint not in C2 */
   Longint*         varsc2weight        /**< pointer to store sum of weights of variables of knapsack constraint in C2 */
   )
{
   int i;

   assert(scip != NULL);
   assert(nvars > 0);
   assert(weights != NULL);
   assert(solvals != NULL);
   assert(varsc2 != NULL);
   assert(varscn2 != NULL);
   assert(nvarsc2 !=NULL);
   assert(nvarscn2 != NULL);
   assert(varsc2weight !=NULL);

   *nvarscn2 = 0;
   *nvarsc2 = 0;
   *varsc2weight = 0;

   /* choses all variables with LP value equal to one for the set of variables for downlifting (C2) */
   for( i = 0; i < nvars; i++ )
   {
      if( SCIPisFeasEQ(scip, solvals[i], 1.0) )
      {
         varsc2[*nvarsc2] = i;
         (*nvarsc2)++;
         *varsc2weight += weights[i];
      }
      else
      {
         assert(SCIPisFeasLT(scip, solvals[i], 1.0));
         varscn2[*nvarscn2] = i;
         (*nvarscn2)++;
      }
   }
}

/** sorts given set of variables in knapsack constraint by non-increasing LP value (variables with equal LP value 
 *  are ordered by non-increasing weight)
 */
static
void sortSetvars(
   SCIP*            scip,               /**< SCIP data structure */
   int*             setvars,            /**< set of variables in knapsack constraint to be sorted */
   int              nsetvars,           /**< number of variables in set of variables in knapsack cons. to be sorted */
   Real*            solvals             /**< LP values of all problem variables */
   )
{
   Real solval;
   int idx;
   int i;
   int j;

   assert(scip != NULL);
   assert(setvars != NULL);
   assert(nsetvars >= 0);
   assert(solvals != NULL);

   /* sorts setvars by non-increasing LP value */
   for( i = 0; i < nsetvars; i++ ) 
   {
      idx = setvars[i];
      solval = solvals[idx];
      for( j = i - 1; j >= 0 && solvals[setvars[j]] < solval; --j )
         setvars[j + 1] = setvars[j];
      setvars[j + 1] = idx;
   }
}

/** enlarges minweight table to at least the given length */
static
RETCODE enlargeMinweighttableSize(
   SCIP*            scip,               /**< SCIP data structure */
   Longint**        fullminweightptr,   /**< pointer to fullminweight table */
   int*             fulltablelen,       /**< pointer to store number of entries in fullminweight table (incl. z=0) */
   int*             fulltablesize,      /**< pointer to current size of fullminweight table */
   int              newlen              /**< new length of fullminweight table */
   )
{
   int j;

   assert(fullminweightptr != NULL);
   assert(*fullminweightptr != NULL);

   if( newlen > *fulltablesize )
   {
      int newsize;

      /* reallocate table memory */
      newsize = MAX(newlen, 2*(*fulltablesize));
      CHECK_OKAY( SCIPreallocBufferArray(scip, fullminweightptr, newsize) );
      *fulltablesize = newsize;
   }

   /* initialize new elements */
   for( j = *fulltablelen; j < newlen; ++j )
      (*fullminweightptr)[j] = LONGINT_MAX;
   *fulltablelen = newlen;

   return SCIP_OKAY;
}

/** lifts up given inequality sum(j in C1) x_j <= liftrhs which is valid for the knapsack polytop { x binary | 
 *  x is solution of knapsack constraint, x_j = 1 for all j in C2, x_j = 0 for all j in noncovervars } to a valid 
 *  inequality for the knapsack polytop { x binary | x is solution of knapsack constraint, x_j = 1 for all j in C2,
 *  x_j = 0 for all j in noncovervars with LP value = 0 }
 */
static
RETCODE liftupKnapsackCover(
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            vars,               /**< variables in knapsack constraint */
   Longint*         weights,            /**< weights of variables in knapsack constraint */
   Longint          capacity,           /**< capacity of knapsack */
   Real*            solvals,            /**< LP values of all problem variables */
   int*             covervars,          /**< cover variables C = C2 & C1 (C2, C1 sorted by non-incr LP val then weight) */
   int*             noncovervars,       /**< noncover variables (sorted by non-incr LP val then weight) */
   int              ncovervars,         /**< number of cover variables */
   int              ncovervarsc1,       /**< number of cover variables in C1 (at the end of covervars) */
   int              ncovervarsc2,       /**< number of cover variables in C2 (at the beginning of covervars) */
   int              nnoncovervars,      /**< number of noncover variables */
   Longint          fixedoneweight,     /**< sum of weights of variables from C2 which are still fixed to 1 */ 
   int*             liftcoefs,          /**< lifting coefficients of var in knapsack constraint in lifted inequality */
   int*             liftrhs,            /**< right hand side of lifted inequality */
   Real*            liftlpval,          /**< LP solution value of lifted variables (without C1) */  
   int*             lastlifted,         /**< pointer to store index of noncover var lifted last (-1 if no var lifted) */
   Longint*         initminweight,      /**< initial minweight table (aggregation of sorted weights) */
   Longint**        fullminweightptr,   /**< pointer to fullminweight table */
   int*             fulltablelen,       /**< pointer to store number of entries in fullminweight table (incl. z=0) */
   int*             fulltablesize       /**< pointer to current size of fullminweight table */
   )
{
   Longint* fullminweight;
   Longint weight;
   Longint rescapacity;
   int liftvar;
   int i;
   int z;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(weights != NULL);
   assert(capacity >= 0);
   assert(solvals != NULL);
   assert(covervars != NULL);
   assert(noncovervars != NULL);
   assert(ncovervars >= 0);
   assert(ncovervarsc1 >= 0);
   assert(ncovervarsc2 >= 0);
   assert(ncovervarsc1 + ncovervarsc2 == ncovervars);
   assert(nnoncovervars >= 0);
   assert(fixedoneweight >= 0);
   assert(liftcoefs != NULL);
   assert(liftrhs != NULL);
   assert(liftlpval != NULL);
   assert(lastlifted != NULL);
   assert(initminweight != NULL);
   assert(fullminweightptr != NULL);
   assert(fulltablelen != NULL);
   assert(fulltablesize != NULL);
   assert(ncovervarsc1 < *fulltablesize);
   assert(initminweight[0] == 0);
   assert(*liftrhs <= ncovervarsc1);

   fullminweight = *fullminweightptr;
   assert(fullminweight != NULL);

   /* initialize the minweight table
    *  fullminweight[z] := minimal sum of weights s.t. activity of cut inequality equals z
    */
   copyMemoryArray(fullminweight, initminweight, ncovervarsc1+1);
   *fulltablelen = ncovervarsc1 + 1;
   assert(*liftrhs < *fulltablelen);
   assert(0 <= *fulltablelen && *fulltablelen <= *fulltablesize);

   /* calculate lifting coefficients for noncover variables with LP value > 0:
    *  for each noncovervar x_i with LP-value > 0 (in non-increasing order of LP-value and weight): 
    *   1. calculates max activity z_max of current lifted inequality s.t. the knapsack is still feasible with 
    *      x_i = 1, x_j = 0 for all j in noncovervars with j > i, x_j = 1 for all j in C2
    *   2. adds x_i with lifting coefficient beta_i = liftrhs - z_max to current lifted inequality
    *   3. updates minweight table: calculates minimal sum of weights s.t. activity of current lifted inequality equals z 
    */
   for( i = 0; i < nnoncovervars && SCIPisFeasPositive(scip, solvals[noncovervars[i]]); i++ )
   {
#ifdef LIFTUPOUT
      printf("----------------------------------- up lifting --------------------------------------------------\n");
      printf("ncovervars=%d (nc2=%d, nc1=%d) && covervars=\n", ncovervars, ncovervarsc2, ncovervarsc1);
      for( z = 0; z < ncovervars; z++ )
         printf("%d: x_%d [w=%lld, lp=%g, lb_global=%g, lb_local=%g, %s]\n", z, covervars[z], weights[covervars[z]], 
            solvals[covervars[z]], SCIPvarGetLbGlobal(vars[covervars[z]]), SCIPvarGetLbLocal(vars[covervars[z]]),
            SCIPvarGetName(vars[covervars[z]]));
      printf("capacity=%lld\n\n", capacity);
      printf("nnoncovervars=%d && noncovervars=\n", nnoncovervars);
      for( z = 0; z < nnoncovervars; z++ )
         printf("%d: x_%d [w=%lld, lp=%g, lb_global=%g, lb_local=%g, %s]\n", z, noncovervars[z], weights[noncovervars[z]], 
            solvals[noncovervars[z]], SCIPvarGetLbGlobal(vars[noncovervars[z]]), SCIPvarGetLbLocal(vars[noncovervars[z]]),
            SCIPvarGetName(vars[noncovervars[z]]));
      printf("\n\n");
      printf("i=%d: noncovervar x_%d:\n\n", i, noncovervars[i]);
      printf("fullminweighttable(liftrhs=%d, fulltablelen=%d)=[ ", *liftrhs, *fulltablelen);
      for( z = 0; z < *fulltablelen; z++ )
         printf("%lld ", fullminweight[z]);
      printf("]\n");
#endif

      liftvar = noncovervars[i];
      assert(liftvar >= 0);
      weight = weights[liftvar];

      /* binary search in sorted fullminweight array for the largest entry that is not greater than 
       * capacity - fixedoneweight - weight_i 
       */
      assert(*liftrhs < *fulltablelen);
      rescapacity = capacity - fixedoneweight - weight;
      if( rescapacity < 0 )
      {
         /* lifted inequality cut is valid for any beta_i, but it may no longer be facet inducing: set beta = liftrhs+1 */
         liftcoefs[liftvar] = (*liftrhs) + 1;
      }
      else if( fullminweight[*liftrhs] <= rescapacity )
      {
         /* weight of variable is too small to get a positive lifting coefficient */
         liftcoefs[liftvar] = 0;
         continue;
      }
      else
      {
         int left;
         int right;
         int middle;
         
         left = 0;
         right = (*liftrhs) + 1;
         while( left < right - 1)
         {
            middle = (left + right) / 2;
            assert(0 <= middle && middle < *fulltablelen);
            if( fullminweight[middle] <= rescapacity ) 
               left = middle;
            else
               right = middle;
         }
         assert(left == right - 1);
         assert(0 <= left && left < *fulltablelen);
         assert(fullminweight[left] <= rescapacity);
         assert(left == (*fulltablelen)-1 || fullminweight[left+1] > rescapacity);
            
         /* now z_max = left: calculate lifting coefficient beta_i */
         liftcoefs[liftvar] = (*liftrhs) - left;
      }
      
      /* update LP solution value of lifted variables (without C1) */
      *liftlpval += liftcoefs[liftvar] * solvals[liftvar];
         
#ifdef LIFTUPOUT
      printf("liftvar=xn%d: lifting coefficient=%d, liftlpval=%g\n", liftvar, liftcoefs[liftvar], *liftlpval);
#endif

      /* if we have downlifting candidates, we have to increase the size of the minweight table;
       * otherwise, we only have to keep track of the elements up to z = liftrhs
       */
      if( ncovervarsc2 > 0 )
      {
         CHECK_OKAY( enlargeMinweighttableSize(scip, fullminweightptr, fulltablelen, fulltablesize,
               *fulltablelen + liftcoefs[liftvar]) );
         fullminweight = *fullminweightptr;
      }

      /* update fullminweight table: calculates minimal sum of weights s.t. activity of curr lifted inequality equals z */
      for( z = (*fulltablelen) - 1; z >= liftcoefs[liftvar]; z-- )
      {
         Longint min;
         min = MIN(fullminweight[z], fullminweight[z - liftcoefs[liftvar]] + weights[liftvar]);
         assert(min > 0);
         fullminweight[z] = min;
      }
   }

   /* update noncover variable lifted last */
   *lastlifted = i-1;

   return SCIP_OKAY;
}

/** lifts down given inequality sum(j in C1 & noncovervars with LP value > 0) beta_j * x_j <= liftrhs which is valid for 
 *  the knapsack polytop { x binary | x is solution of knapsack constraint, x_j = 1 for all j in C2, x_j = 0 for all j in 
 *  noncovervars with LP value 0} to a valid inequality for the knapsack polytop { x binary | x is solution of knapsack 
 *  constraint, x_j = 0 for all j in noncovervars with LP value 0}
 */
static
RETCODE liftdownKnapsackCover(
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            vars,               /**< variables in knapsack constraint */
   Longint*         weights,            /**< weights of variables in knapsack constraint */
   Longint          capacity,           /**< capacity of knapsack */
   Real*            solvals,            /**< LP values of all problem variables */
   int*             covervars,          /**< cover variables C = C2 & C1 (C2, C1 sorted by non-incr LP val then weight) */
   int*             noncovervars,       /**< noncover variables (sorted by non-incr LP val then weight) */
   int              ncovervars,         /**< number of cover variables */
   int              ncovervarsc1,       /**< number of cover variables in C1 (at the end of covervars) */
   int              ncovervarsc2,       /**< number of cover variables in C2 (at the beginning of covervars) */
   int              nnoncovervars,      /**< number of noncover variables */
   Longint          fixedoneweight,     /**< sum of weights of variables from C2 which are still fixed to 1 */ 
   int*             liftcoefs,          /**< lifting coefficients of var in knapsack constraint in lifted inequality */
   int*             liftrhs,            /**< right hand side of lifted inequality */
   Real*            liftlpval,          /**< LP solution value of lifted variables (without C1) */  
   Longint**        fullminweightptr,   /**< pointer to fullminweight table */
   int*             fulltablelen,       /**< pointer to store number of entries in fullminweight table (incl. z=0) */
   int*             fulltablesize       /**< pointer to current size of fullminweight table */
   )
{
   Longint* fullminweight;
   Longint rescapacity;
   int liftvar;
   int left;
   int right;
   int middle;
   int i;
   int z;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(weights != NULL);
   assert(capacity >= 0);
   assert(solvals != NULL);
   assert(covervars != NULL);
   assert(noncovervars != NULL);
   assert(ncovervars >= 0);
   assert(ncovervarsc1 >= 0);
   assert(ncovervarsc2 >= 0);
   assert(ncovervarsc1 + ncovervarsc2 == ncovervars);
   assert(nnoncovervars >= 0);
   assert(fixedoneweight >= 0);
   assert(liftcoefs != NULL);
   assert(liftrhs != NULL);
   assert(liftlpval != NULL);
   assert(fullminweightptr != NULL);
   assert(fulltablelen != NULL);
   assert(0 <= *fulltablelen && *fulltablelen <= *fulltablesize);

   fullminweight = *fullminweightptr;
   assert(fullminweight != NULL);

   /* calculates lifting coefficients for cover variables in C2:
    *  for each covervar x_i in C2 (in non-increasing order of LP-value and weight): 
    *   1. calculates max activity z_max of current lifted inequality s.t. the knapsack is still feasible with 
    *      x_i = 0, x_j = 1 for all covervars in C2 with j > i, x_j = 0 for all j in noncovervars with LP value 0
    *   2. adds x_i with lifting coefficient beta_i = z_max - liftrhs to current lifted inequality
    *   3. changes liftrhs: liftrhs = liftrhs + beta_i
    *   4. updates minweight table: calculates minimal sum of weights s.t. activity of current lifted inequality equals z 
    */
   for( i = 0; i < ncovervarsc2; i++ )
   {
#ifdef LIFTDOWNOUT
      printf("------------------------------ down lifting -------------------------------------------------------\n");
      printf("ncovervars=%d (nc2=%d, nc1=%d) && covervars=\n", ncovervars, ncovervarsc2, ncovervarsc1);
      for( z = 0; z < ncovervars; z++ )
         printf("%d: x_%d [w=%lld, lp=%g, lb_global=%g, lb_local=%g, %s]\n", z, covervars[z], weights[covervars[z]], 
            solvals[covervars[z]], SCIPvarGetLbGlobal(vars[covervars[z]]), SCIPvarGetLbLocal(vars[covervars[z]]),
            SCIPvarGetName(vars[covervars[z]]));
      printf("capacity=%lld\n\n", capacity);
      printf("nnonsetvars=%d && nonsetvars=\n", nnoncovervars);
      for( z = 0; z < nnoncovervars; z++ )
         printf("%d: x_%d [w=%lld, lp=%g, lb_global=%g, lb_local=%g, %s]\n", z, noncovervars[z], weights[noncovervars[z]], 
            solvals[noncovervars[z]], SCIPvarGetLbGlobal(vars[noncovervars[z]]), SCIPvarGetLbLocal(vars[noncovervars[z]]),
            SCIPvarGetName(vars[noncovervars[z]]));
      printf("\n\n");
      printf("i=%d: covervar x_%d, weight_%d=%lld:\n같같같같같같같같같같같같같같같같같같같같�\n", i, covervars[i], 
         covervars[i], weights[covervars[i]]);
      printf("fullminweighttable (fulltablelen=%d)=[ ", *fulltablelen);
      for( z = 0; z < *fulltablelen; z++ )
         printf("%lld ", fullminweight[z]);
      printf("]\n\n");
#endif

      assert(SCIPisFeasEQ(scip, solvals[covervars[i]], 1.0)); 
      liftvar = covervars[i];
      assert(liftvar >= 0);
      
      /* updates sum of weights of variables in knapsack still fixed to 1 */
      fixedoneweight -= weights[liftvar];
    
      /* binary search in sorted fullminweight array for the largest entry that is not greater then 
       * capacity - fixedoneweight
       */
      rescapacity = capacity - fixedoneweight;
      assert(rescapacity >= 0);
      left = 0;
      right = *fulltablelen;
      while( left < right - 1)
      {
         middle = (left + right) / 2;
         assert(0 <= middle && middle < *fulltablelen);
         if( fullminweight[middle] <= rescapacity ) 
            left = middle;
         else
            right = middle;
      }
      assert(left == right - 1);
      assert(0 <= left && left < *fulltablelen);
      assert(fullminweight[left] <= rescapacity);
      assert(left == *fulltablelen - 1 || fullminweight[left+1] > rescapacity);

      /* now z_max = left: calculates lifting coefficient beta_i */
      assert(left >= (*liftrhs));
      liftcoefs[liftvar] = left - (*liftrhs);
      
      /* updates LP solution value of lifted variables (without C1) */
      *liftlpval += liftcoefs[liftvar] * solvals[liftvar];
      
#ifdef LIFTDOWNOUT
      printf("fixedoneweight=%lld, rescapacity=%lld, liftrhs=%d  ", fixedoneweight, rescapacity, *liftrhs);
      printf("==> zmax=%d, lifting coefficient=%d, liftlpval=%g\n------------------------------------------------------------------------------------------------------------------------------\n\n", left, liftcoefs[liftvar], *liftlpval);
#endif

      /* changes liftrhs */
      (*liftrhs) += liftcoefs[liftvar];
      
      /* updates fullminweight table: calculates minimal sum of weights s.t. activity of curr lifted inequality equals z */
      if( liftcoefs[liftvar] == 0 )
         continue;

      CHECK_OKAY( enlargeMinweighttableSize(scip, fullminweightptr, fulltablelen, fulltablesize,
            *fulltablelen + liftcoefs[liftvar]) );
      fullminweight = *fullminweightptr;
      for( z = *fulltablelen - 1; z >= liftcoefs[liftvar]; z-- )
      {
         Longint min;
         min = MIN(fullminweight[z], fullminweight[z - liftcoefs[liftvar]] + weights[liftvar]);
         assert(min > 0);
         fullminweight[z] = min;
      }
   }

   return SCIP_OKAY;
}

/** lifts up given inequality sum(j in C1 & noncovervars with LP value > 0 & C2) beta_j * x_j <= liftrhs which is valid 
 *  for the knapsack polytop { x binary | x is solution of knapsack constraint, x_j = 0 for all j in noncovervars with 
 *  LP value 0} to a valid inequality for the knapsack polytop { x binary | x is solution of knapsack constraint }
 */
static
void liftupZerosKnapsackCover(
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            vars,               /**< variables in knapsack constraint */
   Longint*         weights,            /**< weights of variables in knapsack constraint */
   Longint          capacity,           /**< capacity of knapsack */
   Real*            solvals,            /**< LP values of all problem variables */
   int*             covervars,          /**< cover variables C = C2 & C1 (C2, C1 sorted by non-incr LP val then weight) */
   int*             noncovervars,       /**< noncover variables (sorted by non-incr LP value then weight) */
   int              ncovervars,         /**< number of cover variables */
   int              ncovervarsc1,       /**< number of cover variables in C1 (at the end of covervars) */
   int              ncovervarsc2,       /**< number of cover variables in C2 (at the beginning of covervars) */
   int              nnoncovervars,      /**< number of noncover variables */
   int*             liftcoefs,          /**< lifting coefficients of variables in knapsack cons in lifted inequality */
   int*             liftrhs,            /**< right hand side of lifted inequality */
   Real*            liftlpval,          /**< LP solution value of lifted variables (without C1) */  
   int              lastlifted,         /**< index of last noncover var with LP value > 0 lifted (-1 if no var lifted) */
   Longint*         fullminweight       /**< fullminweight table */
   )
{
   Longint weight;
   Longint rescapacity;
   int liftvar;
   int i;
   int z;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(weights != NULL);
   assert(capacity >= 0);
   assert(solvals != NULL);
   assert(covervars != NULL);
   assert(noncovervars != NULL);
   assert(ncovervars >= 0);
   assert(ncovervarsc1 >= 0);
   assert(ncovervarsc2 >= 0);
   assert(ncovervarsc1 + ncovervarsc2 == ncovervars);
   assert(nnoncovervars >= 0);
   assert(liftcoefs != NULL);
   assert(liftrhs != NULL);
   assert(liftlpval != NULL);
   assert(fullminweight != NULL);

   /* calculates lifting coefficients for noncover variables with LP value = 0:
    *  for each noncovervar x_i with LP-value = 0 (i > lastlifted) (in non-increasing order of LP-value and weight): 
    *   1. calculates max activity z_max of current lifted inequality s.t. the knapsack is still feasible with 
    *      x_i = 1, x_j = 0 for all j in noncovervars with j > i
    *   2. adds x_i with lifting coefficient beta_i = liftrhs - z_max to current lifted inequality
    *   3. updates minweight table: calculates minimal sum of weights s.t. activity of current lifted inequality equals z 
    */
   for( i = lastlifted + 1; i < nnoncovervars; i++ )
   {
#ifdef LIFTUPOUT
      printf("----------------------------------- up lifting zeros --------------------------------------------------\n");
      printf("ncovervars=%d (nc2=%d, nc1=%d) && covervars=\n", ncovervars, ncovervarsc2, ncovervarsc1);
      for( z = 0; z < ncovervars; z++ )
         printf("%d: x_%d [w=%lld, lp=%g, lb_global=%g, lb_local=%g, %s]\n", z, covervars[z], weights[covervars[z]], 
            solvals[covervars[z]], SCIPvarGetLbGlobal(vars[covervars[z]]), SCIPvarGetLbLocal(vars[covervars[z]]),
            SCIPvarGetName(vars[covervars[z]]));
      printf("capacity=%lld\n\n", capacity);
      printf("nnonsetvars=%d && nonsetvars=\n", nnoncovervars);
      for( z = 0; z < nnoncovervars; z++ )
         printf("%d: x_%d [w=%lld, lp=%g, lb_global=%g, lb_local=%g, %s]\n", z, noncovervars[z], weights[noncovervars[z]], 
            solvals[noncovervars[z]], SCIPvarGetLbGlobal(vars[noncovervars[z]]), SCIPvarGetLbLocal(vars[noncovervars[z]]),
            SCIPvarGetName(vars[noncovervars[z]]));
      printf("\n\n");
      printf("i=%d: noncovervar x_%d:\n\n", i, noncovervars[i]);
      printf("fullminweighttable(liftrhs=%d)=[ ", *liftrhs);
      for( z = 0; z < *liftrhs; z++ )
         printf("%lld ", fullminweight[z]);
      printf("]\n");
#endif

      assert(SCIPisFeasLE(scip, solvals[noncovervars[i]], 0.0));  

      liftvar = noncovervars[i];
      assert(liftvar >= 0);
      weight = weights[liftvar];
 
      /* binary search in sorted fullminweight array for the largest entry that is not greater then capacity - weight_i */
      rescapacity = capacity - weight;
      if( rescapacity < 0 )
      {
         /* lifted inequality cut is valid for any beta_i, but it may no longer be facet inducing, sets beta_i = liftrhs+1 */
         liftcoefs[liftvar] = (*liftrhs) + 1;
      }
      else if( fullminweight[*liftrhs] <= rescapacity )
      {
         /* weight of variable is too small to get a positive lifting coefficient */
         liftcoefs[liftvar] = 0;
         continue;
      }
      else
      {
         int left;
         int right;
         int middle;
         
         left = 0;
         right = (*liftrhs) + 1;
         while( left < right - 1)
         {
            middle = (left + right) / 2;
            if( fullminweight[middle] <= rescapacity ) 
               left = middle;
            else
               right = middle;
         }
         assert(left == right - 1);
         assert(fullminweight[left] <= rescapacity);
         assert(left == *liftrhs || fullminweight[left+1] > rescapacity);
            
         /* now z_max = left: calculates lifting coefficient beta_i */
         liftcoefs[liftvar] = (*liftrhs) - left;
      }

      /* updates LP solution value of lifted variables (without C1) */
      *liftlpval += liftcoefs[liftvar] * solvals[liftvar];
         
#ifdef LIFTUPOUT
      printf("liftvar=xn%d: lifting coefficient=%d, liftlpval=%g\n", liftvar, liftcoefs[liftvar], *liftlpval);
#endif

      /* update fullminweight table: calculates minimal sum of weights s.t. activity of curr lifted inequality equals z */
      for( z = *liftrhs; z >= liftcoefs[liftvar]; z-- )
      {
         Longint min;
         min = MIN(fullminweight[z], fullminweight[z - liftcoefs[liftvar]] + weights[liftvar]);
         assert(min > 0);
         fullminweight[z] = min;
      }
   }
}

/** lifts given cardinality inequality sum(j in C1) x_j <= |C1| to a valid inequality of the full dimensional knapsack 
 *  polytop by using uplifting for all variables not in the cover and downlifting for all variables in the cover that 
 *  are fixed to one (C2)
 */
static
RETCODE liftKnapsackCover(
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            vars,               /**< variables in knapsack constraint */
   int              nvars,              /**< number of variables in knapsack constraint */
   Longint*         weights,            /**< weights of variables in knapsack constraint */
   Longint          capacity,           /**< capacity of knapsack */
   Real*            solvals,            /**< LP values of all problem variables */
   int*             covervars,          /**< cover variables C = C2 & C1 (C2, C1 sorted by non-incr LP val then weight) */
   int*             noncovervars,       /**< noncover variables (sorted by non-incr LP val then weight) */
   int              ncovervars,         /**< number of cover variables */
   int              ncovervarsc1,       /**< number of cover variables in C1 (at the end of covervars) */
   int              ncovervarsc2,       /**< number of cover variables in C2 (at the beginning of covervars) */
   int              nnoncovervars,      /**< number of noncover variables */
   Longint*         initminweight,      /**< initial minweight table (aggregation of sorted weights) */
   int*             liftcoefs,          /**< pointer to store lifting coefficient of variables in knapsack constraint */
   int*             liftrhs,            /**< pointer to store right hand side of the lifted cover inequality */
   Real*            liftlpval           /**< pointer to store LP solution value of lifted variables */  
   )
{
   Longint fixedoneweight;
   Longint* fullminweight;
   int fulltablesize;
   int fulltablelen;
   int lastlifted;
   int i;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(weights != NULL);
   assert(capacity >= 0);
   assert(solvals != NULL);
   assert(covervars != NULL);
   assert(noncovervars != NULL);
   assert(ncovervars >= 0);
   assert(ncovervarsc1 >= 0);
   assert(ncovervarsc2 >= 0);
   assert(ncovervarsc1 + ncovervarsc2 == ncovervars);
   assert(nnoncovervars >= 0);
   assert(ncovervars + nnoncovervars <= nvars); /* some noncovervars might have been removed already */
   assert(liftcoefs != NULL);
   assert(liftrhs != NULL);
   assert(liftlpval != NULL);

   /* sets rights hand side of cut for cardinality inequality sum(j in C1) x_j <= |C1| */
   *liftrhs = ncovervarsc1; 

   /* allocates temporary memory */
   fulltablesize = ncovervarsc1+1;
   CHECK_OKAY( SCIPallocBufferArray(scip, &fullminweight, fulltablesize) );

   /* initializes data structures */
   *liftlpval = 0.0;
   clearMemoryArray(liftcoefs, nvars);

   /* sets lifting coefficient of cover variables of C1 */
   for( i = ncovervarsc2; i < ncovervars; i++ )
   {
      assert(liftcoefs[covervars[i]] == 0);
      liftcoefs[covervars[i]] = 1;
   }
   
   /* lifts up all noncover variables with LP value > 0 */
   fixedoneweight = 0;
   for( i = 0; i < ncovervarsc2; i++ )
      fixedoneweight += weights[covervars[i]];
   CHECK_OKAY( liftupKnapsackCover(scip, vars, weights, capacity, solvals, covervars, noncovervars, ncovervars, 
         ncovervarsc1, ncovervarsc2, nnoncovervars, fixedoneweight, liftcoefs, liftrhs, liftlpval, &lastlifted, 
         initminweight, &fullminweight, &fulltablelen, &fulltablesize) );
   assert(fulltablelen <= fulltablesize);
   assert(lastlifted < nnoncovervars);

   /* lifts down all cover variables from C2 */
   CHECK_OKAY( liftdownKnapsackCover(scip, vars, weights, capacity, solvals, covervars, noncovervars, ncovervars,
         ncovervarsc1, ncovervarsc2, nnoncovervars, fixedoneweight, liftcoefs, liftrhs, liftlpval,
         &fullminweight, &fulltablelen, &fulltablesize) );

   /* lifts up all remaining noncover variables (noncover variables with LP value 0) */
   liftupZerosKnapsackCover(scip, vars, weights, capacity, solvals, covervars, noncovervars, ncovervars, ncovervarsc1, 
      ncovervarsc2, nnoncovervars, liftcoefs, liftrhs, liftlpval, lastlifted, fullminweight);

   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &fullminweight);

   return SCIP_OKAY;
}

/** lifts given cardinality inequality sum(j in C1) x_j <= |C1| to a valid inequality of the full dimensional knapsack 
 *  polytop by using uplifting for all variables not in the cover and downlifting for all variables in the cover that 
 *  are fixed to one (C2)
 */
RETCODE SCIPliftKnapsackCover(
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            vars,               /**< variables in knapsack constraint */
   int              nvars,              /**< number of variables in knapsack constraint */
   Longint*         weights,            /**< weights of variables in knapsack constraint */
   Longint          capacity,           /**< capacity of knapsack */
   Real*            solvals,            /**< LP values of all problem variables */
   int*             covervars,          /**< cover variables C = C2 & C1 (C2, C1 sorted by non-incr LP val then weight) */
   int*             noncovervars,       /**< noncover variables (sorted by non-incr LP val then weight) */
   int              ncovervars,         /**< number of cover variables */
   int              ncovervarsc1,       /**< number of cover variables in C1 (at the end of covervars) */
   int              ncovervarsc2,       /**< number of cover variables in C2 (at the beginning of covervars) */
   int              nnoncovervars,      /**< number of noncover variables */
   int*             liftcoefs,          /**< pointer to store lifting coefficient of variables in knapsack constraint */
   int*             liftrhs,            /**< pointer to store right hand side of the lifted cover inequality */
   Real*            liftlpval           /**< pointer to store LP solution value of lifted variables */  
   )
{
   Longint* initminweight;
   int i;

   assert(weights != NULL);
   assert(covervars != NULL);

   /* allocate temporary memory for initial minweight table */
   CHECK_OKAY( SCIPallocBufferArray(scip, &initminweight, ncovervarsc1+1) );
   
   /* sort weights of variables in C1 (leave open position 0 to store a zero afterwards) */
   for( i = 0; i < ncovervarsc1; i++ )
   {
      Longint weight;
      int j;

      assert(ncovervarsc2 + i < ncovervars);
      weight = weights[covervars[ncovervarsc2 + i]];
      for( j = i+1; j > 1 && weight < initminweight[j-1]; --j )
         initminweight[j] = initminweight[j-1];
      initminweight[j] = weight;
   }
   
   /* the initial minweight table is the aggregation of the sorted weights */
   initminweight[0] = 0;
   for( i = 1; i <= ncovervarsc1; ++i )
      initminweight[i] += initminweight[i-1];

   /* calculate the lifted knapsack cover cut */
   CHECK_OKAY( liftKnapsackCover(scip, vars, nvars, weights, capacity, solvals, covervars, noncovervars,
         ncovervars, ncovervarsc1, ncovervarsc2, nnoncovervars, initminweight, liftcoefs, liftrhs, liftlpval) );

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &initminweight);

   return SCIP_OKAY;
}

/** separates lifted cover inequalities for given knapsack problem */
RETCODE SCIPseparateKnapsackCover(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint that originates the knapsack problem */
   VAR**            vars,               /**< variables in knapsack constraint */
   int              nvars,              /**< number of variables in knapsack constraint */
   Longint*         weights,            /**< weights of variables in knapsack constraint */
   Longint          capacity,           /**< capacity of knapsack */
   int              maxnumcardlift,     /**< maximal number of cardinality inequ. lifted per sepa round (-1: unlimited) */
   int*             ncuts               /**< pointer to add up the number of found cuts */
   )
{
   Real* solvals;
   Real solval;
   Real covervarsc1activity;
   Real liftlpval;
   Bool coverfound;
   Longint varsc2weight;
   Longint coverweight;
   Longint weight;
   int* varscn2;
   int* varsc2;
   int* covervars;
   int* noncovervars;
   int* liftcoefs;
   int nvarscn2;
   int nvarsc2;
   int ncovervars;
   int ncovervarsc1;
   int ncovervarsc2;
   int nnoncovervars;
   int liftrhs;
   int loopend;
   int i;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(weights != NULL);
   assert(capacity >= 0);
   assert(ncuts != NULL);

   /* increases age of constraint (age is reset to zero, if a cut was found) */
   CHECK_OKAY( SCIPincConsAge(scip, cons) );
   
   /* allocates temporary memory */
   CHECK_OKAY( SCIPallocBufferArray(scip, &solvals, nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &varscn2, nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &varsc2, nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &covervars, nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &noncovervars, nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &liftcoefs, nvars) );

   /* gets LP solution values of all problem variables */
   CHECK_OKAY( SCIPgetVarSols(scip, nvars, vars, solvals) );

#ifdef SEPARATEOUT
   printf("=================================== separating knapsack constraint <%s>: =================================\n", 
      SCIPconsGetName(cons));
   printf("vars (variables in knapsack constraint)(nvars=%d):\n", nvars);
   for( i = 0; i < nvars; i++ )
      printf("%d: x_%d [w=%lld, lp=%g] <%s>\n", i, i, weights[i], solvals[i], SCIPvarGetName(vars[i]));
   printf("capacity of knapsack = %lld\n", capacity);
#endif   

   /* gets a partition of the set of variables in knapsack constraint in a set of variables for downlifting (varsc2) and 
    * a set of the remaining variables (varscn2) 
    */
   getPartition(scip, nvars, weights, solvals, varsc2, varscn2, &nvarsc2, &nvarscn2, &varsc2weight);
   
   /* gets a most violated minimal cover C = C2 & C1, considering that the variables for downlifting (varsc2) are fixed 
    * to one 
    */
   CHECK_OKAY( getCover(scip, vars, nvars, weights, capacity, solvals, varsc2, varscn2, nvarsc2, nvarscn2, varsc2weight,
         covervars, noncovervars, &ncovervars, &ncovervarsc1, &ncovervarsc2, &nnoncovervars, 
         &coverweight, &covervarsc1activity, &coverfound) );
   assert(ncovervars + nnoncovervars == nvars);

   if( coverfound ) 
   {
      Longint* initminweight;

#if 0
      /* increase C1 if |C1| = 1 */
      if( ncovervarsc2 > 0 && ncovervarsc1 == 1 )
      {
         covervarsc1activity += solvals[covervars[ncovervarsc2-1]];
         ncovervarsc1++;
         ncovervarsc2--;
      }
#endif 
      
      /* allocate temporary memory for initial minweight table */
      CHECK_OKAY( SCIPallocBufferArray(scip, &initminweight, ncovervarsc1+1) );

      /* sort weights of variables in C1 (leave open position 0 to store a zero afterwards) */
      for( i = 0; i < ncovervarsc1; i++ )
      {
         assert(ncovervarsc2 + i < ncovervars);
         weight = weights[covervars[ncovervarsc2 + i]];
         for( j = i+1; j > 1 && weight < initminweight[j-1]; --j )
            initminweight[j] = initminweight[j-1];
         initminweight[j] = weight;
      }

      /* the initial minweight table is the aggregation of the sorted weights */
      initminweight[0] = 0;
      for( i = 1; i <= ncovervarsc1; ++i )
         initminweight[i] += initminweight[i-1];

      /* sort covervars and noncovervars by non-increasing LP value */
      sortSetvars(scip, covervars, ncovervars, solvals);
      sortSetvars(scip, noncovervars, nnoncovervars, solvals);
         
      /* generate lifted cardinality inequalities, consecutively removing variables from the cardinality inequality 
       * (only from C1), starting with the full cover 
       */
      if( maxnumcardlift == -1 )
         loopend = ncovervarsc2;
      else
         loopend = MAX(ncovervarsc2, ncovervars - maxnumcardlift - 1);
      assert(loopend >= ncovervarsc2 && loopend <= ncovervars);
      for( j = ncovervars - 1; j >= loopend; j-- )
      {
         Longint oldweight;
         Longint newweight;
         int n;

#ifdef SEPARATEOUT
         printf(".................................. cardinality round j=%d .......................................\n", j);
         printf("\nbefor card-remove:\n");
         printf("covervars (ncovervars=%d, nc2=%d, nc1=%d):\n", ncovervars, ncovervarsc2, ncovervarsc1);
         for( i = 0; i < ncovervars; i++ )
            printf("%d: x_%d [w=%lld, lp=%g]\n", i, covervars[i], weights[covervars[i]], 
               solvals[covervars[i]]);
         printf("coverweight=%lld, covervarsc1activity=%g\n", coverweight, covervarsc1activity);
         printf("noncovervars (nnoncovervars=%d):\n", nnoncovervars);
         for( i = 0; i < nnoncovervars; i++ )
            printf("%d: x_%d [w=%lld, lp=%g]\n", i, noncovervars[i], weights[noncovervars[i]], 
               solvals[noncovervars[i]]);
#endif
         /* delete current variable from C1 */
         weight = weights[covervars[j]];
         solval = solvals[covervars[j]];
         coverweight -= weight;
         covervarsc1activity -= solval;
         ncovervars--;
         ncovervarsc1--;
         assert(ncovervars == j);
         assert(ncovervarsc1 + ncovervarsc2 == ncovervars);
         assert(SCIPisFeasGE(scip, covervarsc1activity, 0.0));

         /* remove the entry from the initial minweight table */
         oldweight = initminweight[ncovervarsc1+1];
         for( i = ncovervarsc1; i > 0 && oldweight - initminweight[i] != weight; --i )
         {
            newweight = oldweight - weight;
            oldweight = initminweight[i];
            assert(newweight > oldweight);
            initminweight[i] = newweight;
         }

         /* add current variable to noncovervars */
         for( i = nnoncovervars - 1; i >= 0 && solvals[noncovervars[i]] < solval; i-- )
            noncovervars[i + 1] = noncovervars[i];
         noncovervars[i + 1] = covervars[j];
         nnoncovervars++;
         assert(nnoncovervars + ncovervars <= nvars); /* some variables might have been removed from noncovervars */

         /* remove variables from noncovervars, that would fit into the knapsack together with all covervars;
          * the corresponding lifting coefficients would be zero
          */
         for( i = 0, n = 0; i < nnoncovervars; i++ )
         {
            if( coverweight + weights[noncovervars[i]] > capacity )
            {
               noncovervars[n] = noncovervars[i];
               n++;
            }
         }
         nnoncovervars = n;
            
#ifdef SEPARATEOUT
         printf("\nafter card-remove:\n");
         printf("covervars (ncovervars=%d, nc2=%d, nc1=%d):\n", ncovervars, ncovervarsc2, ncovervarsc1);
         for( i = 0; i < ncovervars; i++ )
            printf("%d: x_%d [w=%lld, lp=%g, beta=%d]\n", i, covervars[i], weights[covervars[i]], 
               solvals[covervars[i]], liftcoefs[covervars[i]]);
         printf("coverweight=%lld, covervarsc1activity=%g\n", coverweight, covervarsc1activity);
         printf("noncovervars (nnoncovervars=%d):\n", nnoncovervars);
         for( i = 0; i < nnoncovervars; i++ )
            printf("%d: x_%d [w=%lld, lp=%g, beta=%d]\n", i, noncovervars[i], weights[noncovervars[i]], 
               solvals[noncovervars[i]], liftcoefs[noncovervars[i]]);
#endif

         /* lifts cardinality inequality sum(j in C1) x_j <= |C1| to a valid inequality of the full dimensional knapsack 
          * polytop: 
          *  1. uplifting of noncovervars with LP value > 0 
          *  2. downlifting of covervarsc2 (C2)
          *  3. uplifting of noncovervars with LP value = 0 
          */
         CHECK_OKAY( liftKnapsackCover(scip, vars, nvars, weights, capacity, solvals, covervars, noncovervars, 
               ncovervars, ncovervarsc1, ncovervarsc2, nnoncovervars, initminweight,
               liftcoefs, &liftrhs, &liftlpval) );
         
         /* checks, if lifting yielded a violated cut */
         if( SCIPisEfficacious(scip, (covervarsc1activity + liftlpval - liftrhs)/sqrt((Real)liftrhs)) )
         {
            ROW* row;
            char name[MAXSTRLEN];
            int v;
            
            /* creates LP row */
            sprintf(name, "%s_card%lld_%d", SCIPconsGetName(cons), SCIPconshdlrGetNCutsFound(SCIPconsGetHdlr(cons)), j);
            CHECK_OKAY( SCIPcreateEmptyRow (scip, &row, name, -SCIPinfinity(scip), (Real)liftrhs, 
                  SCIPconsIsLocal(cons), FALSE, SCIPconsIsRemoveable(cons)) );
            
            /* adds all variables in the knapsack constraint with calculated lifting coefficient to the cut */
            CHECK_OKAY( SCIPcacheRowExtensions(scip, row) );
            for( v = 0; v < ncovervars; v++ )
            {
               CHECK_OKAY( SCIPaddVarToRow(scip, row, vars[covervars[v]], (Real)liftcoefs[covervars[v]]) ); 
            }
            for( v = 0; v < nnoncovervars; v++ )
            {
               if( liftcoefs[noncovervars[v]] > 0 )
               {
                  CHECK_OKAY( SCIPaddVarToRow(scip, row, vars[noncovervars[v]], (Real)liftcoefs[noncovervars[v]]) );
               }
            }
            CHECK_OKAY( SCIPflushRowExtensions(scip, row) );
            
            /* checks, if cut is violated enough */
            if( SCIPisCutEfficacious(scip, row) )
            {         
#ifdef CUTOUT
               printf("lifted cover cut for knapsack constraint <%s> round j=%d: ", SCIPconsGetName(cons), j);
               SCIProwPrint(row, NULL);
               printf("violation = %g", covervarsc1activity + liftlpval - liftrhs);
               printf("\n");

#endif               
               CHECK_OKAY( SCIPresetConsAge(scip, cons) );
               CHECK_OKAY( SCIPaddCut(scip, row, FALSE) );
               (*ncuts)++;
            }
            CHECK_OKAY( SCIPreleaseRow(scip, &row) );
         }

         /* if no variable in noncovervar is left that does not fit in the knapsack, we can break the separation loop */
         if( nnoncovervars == 0 )
            break;
      }

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &initminweight);
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &liftcoefs);
   SCIPfreeBufferArray(scip, &noncovervars);
   SCIPfreeBufferArray(scip, &covervars);
   SCIPfreeBufferArray(scip, &varsc2);
   SCIPfreeBufferArray(scip, &varscn2);
   SCIPfreeBufferArray(scip, &solvals);

   return SCIP_OKAY;
}

/** separates given knapsack constraint */
static
RETCODE separateCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< knapsack constraint */
   Bool             sepacardinality,    /**< should knapsack cardinality cuts be separated? */
   int              maxnumcardlift,     /**< maximal number of cardinality inequ. lifted per sepa round (-1: unlimited) */
   int*             ncuts               /**< pointer to add up the number of found cuts */
   )
{
   CONSDATA* consdata;
   Bool violated;

   assert(ncuts != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   debugMessage("separating knapsack constraint <%s>\n", SCIPconsGetName(cons));
   
   /* check knapsack constraint itself for feasibility */
   CHECK_OKAY( checkCons(scip, cons, NULL, FALSE, &violated) );
   
   if( violated )
   {
      /* add knapsack constraint as LP row to the LP */
      CHECK_OKAY( addRelaxation(scip, cons) );
      (*ncuts)++;
   }
   else if( sepacardinality )
   {
      CHECK_OKAY( SCIPseparateKnapsackCover(scip, cons, consdata->vars, consdata->nvars, consdata->weights, 
            consdata->capacity, maxnumcardlift, ncuts) );
   }
   
   return SCIP_OKAY;
}

/** propagation method for knapsack constraints */
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
   Longint onesweightsum;
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

   debugMessage("propagating knapsack constraint <%s>\n", SCIPconsGetName(cons));

   /* increase age of constraint; age is reset to zero, if a conflict or a propagation was found */
   CHECK_OKAY( SCIPincConsAge(scip, cons) );

   do
   {
      /* store the sum of weights of fixed-to-one variables locally, because the constraint data's sum can
       * change due to event processing
       */
      onesweightsum = consdata->onesweightsum;

      /* check, if weights of fixed variables already exceeds knapsack capacity */
      if( consdata->capacity < onesweightsum )
      {
         debugMessage(" -> cutoff - fixed weight: %lld, capacity: %lld\n", onesweightsum, consdata->capacity);
            
         CHECK_OKAY( SCIPresetConsAge(scip, cons) );
         *cutoff = TRUE;

         if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
         {
            /* start conflict analysis with the fixed-to-one variables */
            CHECK_OKAY( SCIPinitConflictAnalysis(scip) );
            for( i = 0; i < consdata->nvars; i++ )
            {
               if( SCIPvarGetLbLocal(consdata->vars[i]) > 0.5)
               {
                  CHECK_OKAY( SCIPaddConflictBinvar(scip, consdata->vars[i]) );
               }
            }
         
            CHECK_OKAY( SCIPanalyzeConflictCons(scip, cons, NULL) );
         }

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
                  debugMessage(" -> fixing variable <%s> to 0\n", SCIPvarGetName(consdata->vars[i]));
                  CHECK_OKAY( SCIPresetConsAge(scip, cons) );
                  CHECK_OKAY( SCIPinferBinvarCons(scip, consdata->vars[i], FALSE, cons, i, &infeasible, &tightened) );
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
      assert(consdata->onesweightsum >= onesweightsum);
   }
   while( consdata->onesweightsum > onesweightsum );

   /* if the remaining (potentially unfixed) variables would fit all into the knapsack, the knapsack is now redundant */
   if( consdata->weightsum - zerosweightsum <= consdata->capacity )
   {
      debugMessage(" -> knapsack constraint <%s> is redundant: weightsum=%lld, zerosweightsum=%lld, capacity=%lld\n",
         SCIPconsGetName(cons), consdata->weightsum, zerosweightsum, consdata->capacity);
      CHECK_OKAY( SCIPdelConsLocal(scip, cons) );
      *redundant = TRUE;
   }

   /* mark the constraint propagated */
   consdata->propagated = TRUE;

   return SCIP_OKAY;
}

/** adds coefficient to constraint data */
static
RETCODE addCoef(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< knapsack constraint */
   VAR*             var,                /**< variable to add to knapsack */
   Longint          weight              /**< weight of variable in knapsack */
   )
{
   CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
   assert(weight > 0);

   /* add the new coefficient to the LP row */
   if( consdata->row != NULL )
   {
      CHECK_OKAY( SCIPaddVarToRow(scip, consdata->row, var, (Real)weight) );
   }

   /* check for fixed variable */
   if( SCIPvarGetLbGlobal(var) > 0.5 )
   {
      /* variable is fixed to one: reduce capacity */
      consdata->capacity -= weight;
   }
   else if( SCIPvarGetUbGlobal(var) > 0.5 )
   {
      Bool negated;

      /* get binary representative of variable */
      CHECK_OKAY( SCIPgetBinvarRepresentative(scip, var, &var, &negated) );

      /* insert coefficient */
      CHECK_OKAY( consdataEnsureVarsSize(scip, consdata, consdata->nvars+1, SCIPconsIsTransformed(cons)) );
      consdata->vars[consdata->nvars] = var;
      consdata->weights[consdata->nvars] = weight;
      consdata->nvars++;

      /* install the rounding locks of variable */
      CHECK_OKAY( lockRounding(scip, cons, var) );

      /* catch events */
      if( SCIPconsIsTransformed(cons) )
      {
         CONSHDLRDATA* conshdlrdata;

         conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
         assert(conshdlrdata != NULL);
         CHECK_OKAY( eventdataCreate(scip, &consdata->eventdatas[consdata->nvars-1], consdata, weight) );
         CHECK_OKAY( SCIPcatchVarEvent(scip, var,
               SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_UBRELAXED | SCIP_EVENTTYPE_VARFIXED,
               conshdlrdata->eventhdlr, consdata->eventdatas[consdata->nvars-1], NULL) );
      }

      /* update weight sums */
      consdata->weightsum += weight;
      if( SCIPvarGetLbLocal(var) > 0.5 )
         consdata->onesweightsum += weight;

      consdata->sorted = FALSE;
      consdata->merged = FALSE;
   }
   consdata->propagated = FALSE;
   consdata->presolved = FALSE;

   return SCIP_OKAY;
}

/** deletes coefficient at given position from constraint data */
static
RETCODE delCoefPos(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< knapsack constraint */
   int              pos                 /**< position of coefficient to delete */
   )
{
   CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(consdata->vars[pos]));

   /* delete the coefficient from the LP row */
   if( consdata->row != NULL )
   {
      CHECK_OKAY( SCIPaddVarToRow(scip, consdata->row, consdata->vars[pos], -(Real)consdata->weights[pos]) );
   }

   /* remove the rounding locks of variable */
   CHECK_OKAY( unlockRounding(scip, cons, consdata->vars[pos]) );

   /* drop events */
   if( SCIPconsIsTransformed(cons) )
   {
      CONSHDLRDATA* conshdlrdata;
      
      conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
      assert(conshdlrdata != NULL);
      CHECK_OKAY( SCIPdropVarEvent(scip, consdata->vars[pos],
            SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_UBRELAXED | SCIP_EVENTTYPE_VARFIXED,
            conshdlrdata->eventhdlr, consdata->eventdatas[pos], -1) );
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
   consdata->presolved = FALSE;
   consdata->sorted = FALSE;

   return SCIP_OKAY;
}

/** replaces multiple occurrences of a variable or its negation by a single coefficient */
static
RETCODE mergeMultiples(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< linear constraint */
   )
{
   CONSDATA* consdata;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->merged )
      return SCIP_OKAY;

   /* loop backwards through the items: deletion only affects rear items */
   for( v = consdata->nvars-1; v >= 0; --v )
   {
      VAR* var;
      VAR* negvar;
      int w;

      var = consdata->vars[v];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
      assert(SCIPvarIsActive(var) || SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED);
      negvar = SCIPvarGetNegatedVar(var);
      assert(negvar == NULL || SCIPvarIsActive(var) || SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED);

      for( w = v-1; w >= 0; --w )
      {
         if( consdata->vars[w] == var )
         {
            /* variables v and w are equal: add weight of v to w, and delete v */
            consdataChgWeight(consdata, w, consdata->weights[v] + consdata->weights[w]);
            CHECK_OKAY( delCoefPos(scip, cons, v) );
            break;
         }
         else if( consdata->vars[w] == negvar )
         {
            /* variables v and w are opposite: subtract smaller weight from larger weight, reduce capacity,
             * and delete item of smaller weight
             */
            if( consdata->weights[v] == consdata->weights[w] )
            {
               /* both variables eliminate themselves: w*x + w*(1-x) == w */
               consdata->capacity -= consdata->weights[v];
               CHECK_OKAY( delCoefPos(scip, cons, v) ); /* this does not affect w, because w < v */
               assert(consdata->vars[w] == negvar);
               CHECK_OKAY( delCoefPos(scip, cons, w) );
               v = MIN(v, consdata->nvars); /* we could have removed the last two coefficients */
            }
            else if( consdata->weights[v] < consdata->weights[w] )
            {
               consdata->capacity -= consdata->weights[v];
               consdataChgWeight(consdata, w, consdata->weights[w] - consdata->weights[v]);
               CHECK_OKAY( delCoefPos(scip, cons, v) ); /* this does not affect w, because w < v */
               assert(consdata->vars[w] == negvar);
               assert(consdata->weights[w] > 0);
            }
            else
            {
               consdata->capacity -= consdata->weights[w];
               consdataChgWeight(consdata, v, consdata->weights[v] - consdata->weights[w]);
               CHECK_OKAY( delCoefPos(scip, cons, w) ); /* this may move v, but does no harm */
            }
            break;
         }
      }
   }

   consdata->merged = TRUE;

   return SCIP_OKAY;
}

/** deletes all fixed variables from knapsack constraint, and replaces variables with binary representatives */
static
RETCODE applyFixings(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< knapsack constraint */
   )
{
   CONSDATA* consdata;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);

   v = 0;
   while( v < consdata->nvars )
   {
      VAR* var;

      var = consdata->vars[v];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);

      if( SCIPvarGetLbGlobal(var) > 0.5 )
      {
         assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(var), 1.0));
         consdata->capacity -= consdata->weights[v];
         CHECK_OKAY( delCoefPos(scip, cons, v) );
      }
      else if( SCIPvarGetUbGlobal(var) < 0.5 )
      {
         assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 0.0));
         CHECK_OKAY( delCoefPos(scip, cons, v) );
      }
      else
      {
         VAR* repvar;
         Longint weight;
         Bool negated;
         
         /* get binary representative of variable */
         CHECK_OKAY( SCIPgetBinvarRepresentative(scip, var, &repvar, &negated) );

         /* check, if the variable should be replaced with the representative */
         if( repvar != var )
         {
            weight = consdata->weights[v];

            /* delete old (aggregated) variable */
            CHECK_OKAY( delCoefPos(scip, cons, v) );

            /* add representative instead */
            CHECK_OKAY( addCoef(scip, cons, repvar, weight) );
         }
         else
            ++v;
      }
   }
   assert(consdata->onesweightsum == 0);

   debugMessage("after fixings:\n");
   debug(CHECK_OKAY( SCIPprintCons(scip, cons, NULL) ));

   /* if aggregated variables have been replaced, multiple entries of the same variable are possible and we have
    * to clean up the constraint
    */
   CHECK_OKAY( mergeMultiples(scip, cons) );
   
   debugMessage("after merging:\n");
   debug(CHECK_OKAY( SCIPprintCons(scip, cons, NULL) ));

   return SCIP_OKAY;
}

/** divides weights by their greatest common divisor and divides capacity by the same value, rounding down the result */
static
void normalizeWeights(
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


/** initialization method of constraint handler (called after problem was transformed) */
#define consInitKnapsack NULL


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitKnapsack NULL


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#define consInitpreKnapsack NULL


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreKnapsack NULL


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolKnapsack NULL


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
DECL_CONSEXITSOL(consExitsolKnapsack)
{  /*lint --e{715}*/
   CONSDATA* consdata;
   int c;

   /* release the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->row != NULL )
      {
         CHECK_OKAY( SCIPreleaseRow(scip, &consdata->row) );
      }
   }

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
DECL_CONSDELETE(consDeleteKnapsack)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);
   
   /* free linear constraint */
   CHECK_OKAY( consdataFree(scip, consdata, conshdlrdata->eventhdlr) );
   
   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */ 
static
DECL_CONSTRANS(consTransKnapsack)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* sourcedata;
   CONSDATA* targetdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->row == NULL);  /* in original problem, there cannot be LP rows */

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* create target constraint data */
   CHECK_OKAY( consdataCreate(scip, &targetdata, conshdlrdata->eventhdlr,
         sourcedata->nvars, sourcedata->vars, sourcedata->weights, sourcedata->capacity) ); 

   /* create target constraint */
   CHECK_OKAY( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), 
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemoveable(sourcecons)) );

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
   Bool sepacardinality;
   int depth;
   int nrounds;
   int sepafreq;
   int sepacardfreq;
   int ncuts;
   int maxsepacuts;
   int i;

   *result = SCIP_DIDNOTRUN;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   depth = SCIPgetDepth(scip);
   nrounds = SCIPgetNSepaRounds(scip);
   
   debugMessage("knapsack separation of %d/%d constraints, round %d (max %d/%d)\n",
      nusefulconss, nconss, nrounds, conshdlrdata->maxroundsroot, conshdlrdata->maxrounds);

   /* only call the separator a given number of times at each node */
   if( (depth == 0 && conshdlrdata->maxroundsroot >= 0 && nrounds >= conshdlrdata->maxroundsroot)
      || (depth > 0 && conshdlrdata->maxrounds >= 0 && nrounds >= conshdlrdata->maxrounds) )
      return SCIP_OKAY;

   /* check, if we should additionally separate cardinality cuts */
   sepafreq = SCIPconshdlrGetSepaFreq(conshdlr);
   sepacardfreq = sepafreq * conshdlrdata->sepacardfreq;
   sepacardinality = (conshdlrdata->sepacardfreq >= 0)
      && ((sepacardfreq == 0 && depth == 0) || (sepacardfreq >= 1 && (depth % sepacardfreq == 0)));

   /* get the maximal number of cuts allowed in a separation round */
   maxsepacuts = (depth == 0 ? conshdlrdata->maxsepacutsroot : conshdlrdata->maxsepacuts);

   *result = SCIP_DIDNOTFIND;
   ncuts = 0;

   /* separate useful constraints */
   for( i = 0; i < nusefulconss && ncuts < maxsepacuts; i++ )
   {
      CHECK_OKAY( separateCons(scip, conss[i], sepacardinality, conshdlrdata->maxnumcardlift, &ncuts) );
   }
   
   /* adjust return value */
   if( ncuts > 0 )
      *result = SCIP_SEPARATED;
   
   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
DECL_CONSENFOLP(consEnfolpKnapsack)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   Bool violated;
   int maxncuts;
   int ncuts;
   int i;

   *result = SCIP_FEASIBLE;

   debugMessage("knapsack enforcement of %d/%d constraints\n", nusefulconss, nconss);

   /* get maximal number of cuts per round */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   maxncuts = (SCIPgetDepth(scip) == 0 ? conshdlrdata->maxsepacutsroot : conshdlrdata->maxsepacuts);

   ncuts = 0;

   /* search for violated useful knapsack constraints */
   for( i = 0; i < nusefulconss && ncuts < maxncuts; i++ )
   {
      CHECK_OKAY( checkCons(scip, conss[i], NULL, FALSE, &violated) );
      if( violated )
      {
         /* add knapsack constraint as LP row to the LP */
         CHECK_OKAY( addRelaxation(scip, conss[i]) );
         ncuts++;
      }
   } 

   /* as long as no violations were found, search for violated obsolete knapsack constraints */
   for( i = nusefulconss; i < nconss && ncuts == 0; i++ )
   {
      CHECK_OKAY( checkCons(scip, conss[i], NULL, FALSE, &violated) );
      if( violated )
      {
         /* add knapsack constraint as LP row to the LP */
         CHECK_OKAY( addRelaxation(scip, conss[i]) );
         ncuts++;
      }
   } 

   /* adjust the result code */
   if( ncuts > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
DECL_CONSENFOPS(consEnfopsKnapsack)
{  /*lint --e{715}*/
   Bool violated;
   int i;

   for( i = 0; i < nconss; i++ )
   {
      CHECK_OKAY( checkCons(scip, conss[i], NULL, TRUE, &violated) );
      if( violated )
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
   Bool violated;
   int i;

   for( i = 0; i < nconss; i++ )
   {
      CHECK_OKAY( checkCons(scip, conss[i], sol, checklprows, &violated) );
      if( violated )
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

   /* process useful constraints */
   for( i = 0; i < nusefulconss && !cutoff; i++ )
   {
      CHECK_OKAY( propagateCons(scip, conss[i], &cutoff, &redundant, &nfixedvars) );
   } 

   /* adjust result code */
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
   CONSDATA* consdata;
   Bool cutoff;
   Bool redundant;
   int oldnfixedvars;
   int oldndelconss;
   int oldnchgcoefs;
   int oldnchgsides;
   int i;

   /* remember old preprocessing counters */
   cutoff = FALSE;
   oldnfixedvars = *nfixedvars;
   oldndelconss = *ndelconss;
   oldnchgcoefs = *nchgcoefs;
   oldnchgsides = *nchgsides;

   for( i = 0; i < nconss; i++ )
   {
      CONS* cons;
      int thisnfixedvars;

      cons = conss[i];
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* force presolving the constraint in the initial round */
      if( nrounds == 0 )
         consdata->presolved = FALSE;

      if( consdata->presolved )
         continue;

      debugMessage("presolving knapsack constraint <%s>\n", SCIPconsGetName(cons));
      debug(CHECK_OKAY( SCIPprintCons(scip, cons, NULL) ));

      consdata->presolved = TRUE;

      /* remove all fixed variables */
      if( nrounds == 0 || nnewfixedvars > 0 || nnewaggrvars > 0 || *nfixedvars > oldnfixedvars )
      {
         CHECK_OKAY( applyFixings(scip, cons) );
      }

      /* propagate constraint */
      thisnfixedvars = *nfixedvars;
      CHECK_OKAY( propagateCons(scip, cons, &cutoff, &redundant, nfixedvars) );
      if( cutoff )
         break;
      if( redundant )
      {
         (*ndelconss)++;
         continue;
      }

      /* remove again all fixed variables, if further fixings were found */
      if( *nfixedvars > thisnfixedvars )
      {
         CHECK_OKAY( applyFixings(scip, cons) );
      }

      if( !SCIPconsIsModifiable(cons) )
      {
         /* divide weights by their greatest common divisor */
         normalizeWeights(cons, nchgcoefs, nchgsides);

         /* tighten capacity and weights */
         tightenWeights(cons, nchgcoefs, nchgsides);
      }
   } 

   /**@todo preprocess pairs of knapsack constraints */

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( *nfixedvars > oldnfixedvars || *ndelconss > oldndelconss
      || *nchgcoefs > oldnchgcoefs || *nchgsides > oldnchgsides )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** propagation conflict resolving method of constraint handler */
static
DECL_CONSRESPROP(consRespropKnapsack)
{  /*lint --e{715}*/
   CONSDATA* consdata;
   Longint capsum;
   int i;

   assert(result != NULL);
   assert(inferinfo >= 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
 
   assert(SCIPvarGetUbLocal(infervar) < 0.5);

   /* locate the inference variable and calculate the capacity that has to be used up to conclude infervar == 0;
    * inferinfo stores the position of the inference variable (but maybe the variables were resorted)
    */
   if( inferinfo < consdata->nvars && consdata->vars[inferinfo] == infervar )
      capsum = consdata->weights[inferinfo];
   else
   {
      for( i = 0; i < consdata->nvars && consdata->vars[i] != infervar; ++i )
      {}
      assert(i < consdata->nvars);
      capsum = consdata->weights[i];
   }

   /* add fixed-to-one variables up to the point, that their weight plus the weight of the conflict variable exceeds
    * the capacity
    */
   if( capsum <= consdata->capacity )
   {
      for( i = 0; i < consdata->nvars; i++ )
      {
         if( SCIPvarGetLbAtIndex(consdata->vars[i], bdchgidx, FALSE) > 0.5 )
         {
            CHECK_OKAY( SCIPaddConflictBinvar(scip, consdata->vars[i]) );
            capsum += consdata->weights[i];
            if( capsum > consdata->capacity )
               break;
         }
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
      CHECK_OKAY( SCIPaddVarLocks(scip, consdata->vars[i], nlocksneg, nlockspos) );
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#define consActiveKnapsack NULL


/** constraint deactivation notification method of constraint handler */
#define consDeactiveKnapsack NULL


/** constraint enabling notification method of constraint handler */
#define consEnableKnapsack NULL


/** constraint disabling notification method of constraint handler */
#define consDisableKnapsack NULL

/** constraint display method of constraint handler */
static
DECL_CONSPRINT(consPrintKnapsack)
{  /*lint --e{715}*/
   CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   for( i = 0; i < consdata->nvars; ++i )
   {
      if( i > 0 )
         SCIPinfoMessage(scip, file, " ");
      SCIPinfoMessage(scip, file, "%+lld<%s>", consdata->weights[i], SCIPvarGetName(consdata->vars[i]));
   }
   SCIPinfoMessage(scip, file, " <= %lld\n", consdata->capacity);

   return SCIP_OKAY;
}




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
   Bool             dynamic,            /**< is constraint subject to aging? */
   Bool             removeable          /**< should the relaxation be removed from the LP due to aging or cleanup? */
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
      capacity = (Longint)SCIPfeasFloor(scip, -lhs);
   }
   else
   {
      mult = +1;
      capacity = (Longint)SCIPfeasFloor(scip, rhs);
   }

   /* negate positive or negative variables */
   for( v = 0; v < nvars; ++v )
   {
      assert(SCIPisFeasIntegral(scip, vals[v]));
      weight = mult * (Longint)SCIPfeasFloor(scip, vals[v]);
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
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removeable) );

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &weights);
   SCIPfreeBufferArray(scip, &transvars);

   return SCIP_OKAY;
}

/** tries to upgrade a linear constraint into a knapsack constraint */
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
      
      /* create the knapsack constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      CHECK_OKAY( createNormalizedKnapsack(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, lhs, rhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), 
            SCIPconsIsDynamic(cons), SCIPconsIsRemoveable(cons)) );
   }

   return SCIP_OKAY;
}




/*
 * Event handler
 */

/** execution methode of bound change event handler */
static
DECL_EVENTEXEC(eventExecKnapsack)
{  /*lint --e{715}*/
   assert(eventdata != NULL);
   assert(eventdata->consdata != NULL);
   
   switch( SCIPeventGetType(event) )
   {
   case SCIP_EVENTTYPE_LBTIGHTENED:
      eventdata->consdata->onesweightsum += eventdata->weight;
      eventdata->consdata->propagated = FALSE;
      eventdata->consdata->presolved = FALSE;
      break;
   case SCIP_EVENTTYPE_LBRELAXED:
      eventdata->consdata->onesweightsum -= eventdata->weight;
      break;
   case SCIP_EVENTTYPE_UBRELAXED:
      /* if a variable fixed to 0 is unfixed, it is possible, that it can be fixed to 0 again */
      eventdata->consdata->propagated = FALSE;
      break;
   case SCIP_EVENTTYPE_VARFIXED:
      /* the variable should be removed from the constraint in presolving */
      eventdata->consdata->presolved = FALSE;
      break;
   default:
      errorMessage("invalid event type %x\n", SCIPeventGetType(event));
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
   EVENTHDLRDATA* eventhdlrdata;
   CONSHDLRDATA* conshdlrdata;

   /* include event handler for bound change events */
   eventhdlrdata = NULL;
   CHECK_OKAY( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC, 
         NULL, NULL, NULL, NULL, NULL, NULL, eventExecKnapsack,
         eventhdlrdata) );

   /* create knapsack constraint handler data */
   CHECK_OKAY( SCIPallocMemory(scip, &conshdlrdata) );

   /* get event handler for bound change events */
   conshdlrdata->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( conshdlrdata->eventhdlr == NULL )
   {
      errorMessage("event handler for knapsack constraints not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         consFreeKnapsack, consInitKnapsack, consExitKnapsack, 
         consInitpreKnapsack, consExitpreKnapsack, consInitsolKnapsack, consExitsolKnapsack,
         consDeleteKnapsack, consTransKnapsack, consInitlpKnapsack,
         consSepaKnapsack, consEnfolpKnapsack, consEnfopsKnapsack, consCheckKnapsack, 
         consPropKnapsack, consPresolKnapsack, consRespropKnapsack, consLockKnapsack,
         consActiveKnapsack, consDeactiveKnapsack, 
         consEnableKnapsack, consDisableKnapsack,
         consPrintKnapsack,
         conshdlrdata) );
  
   /* include the linear constraint upgrade in the linear constraint handler */
   CHECK_OKAY( SCIPincludeLinconsUpgrade(scip, linconsUpgdKnapsack, LINCONSUPGD_PRIORITY) );

   /* add knapsack constraint handler parameters */
   CHECK_OKAY( SCIPaddIntParam(scip,
         "constraints/knapsack/sepacardfreq",
         "multiplier on separation frequency, how often cardinality cuts are separated (-1: never, 0: only at root)",
         &conshdlrdata->sepacardfreq, DEFAULT_SEPACARDFREQ, -1, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "constraints/knapsack/maxrounds",
         "maximal number of separation rounds per node (-1: unlimited)",
         &conshdlrdata->maxrounds, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "constraints/knapsack/maxroundsroot",
         "maximal number of separation rounds per node in the root node (-1: unlimited)",
         &conshdlrdata->maxroundsroot, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "constraints/knapsack/maxsepacuts",
         "maximal number of cuts separated per separation round",
         &conshdlrdata->maxsepacuts, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "constraints/knapsack/maxsepacutsroot",
         "maximal number of cuts separated per separation round in the root node",
         &conshdlrdata->maxsepacutsroot, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "constraints/knapsack/maxnumcardlift",
         "maximal number of cardinality inequalities lifted per separation round (-1: unlimited)",
         &conshdlrdata->maxnumcardlift, DEFAULT_MAXNUMCARDLIFT, -1, INT_MAX, NULL, NULL) );

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
   Bool             dynamic,            /**< is constraint subject to aging? */
   Bool             removeable          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   )
{
   CONSHDLRDATA* conshdlrdata;
   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   /* find the knapsack constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("knapsack constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* create constraint data */
   CHECK_OKAY( consdataCreate(scip, &consdata, conshdlrdata->eventhdlr, nvars, vars, weights, capacity) );
        
   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removeable) );

   return SCIP_OKAY;
}

/** adds new item to knapsack constraint */
RETCODE SCIPaddCoefKnapsack(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   VAR*             var,                /**< item variable */
   Longint          weight              /**< item weight */
   )
{
   assert(var != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not a knapsack constraint\n");
      return SCIP_INVALIDDATA;
   }
   
   CHECK_OKAY( addCoef(scip, cons, var, weight) );

   return SCIP_OKAY;
}

/** gets the dual solution of the knapsack constraint in the current LP */
Real SCIPgetDualsolKnapsack(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   )
{
   CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not a knapsack constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualsol(consdata->row);
   else
      return 0.0;
}

/** gets the dual farkas value of the knapsack constraint in the current infeasible LP */
Real SCIPgetDualfarkasKnapsack(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   )
{
   CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not a knapsack constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualfarkas(consdata->row);
   else
      return 0.0;
}


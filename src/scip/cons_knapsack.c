/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_knapsack.c,v 1.181 2009/09/10 13:47:15 bzfwolte Exp $"

/**@file   cons_knapsack.c
 * @ingroup CONSHDLRS 
 * @brief  constraint handler for knapsack constraints
 * @author Tobias Achterberg
 * @author Kati Wolter
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "scip/pub_misc.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "knapsack"
#define CONSHDLR_DESC          "knapsack constraint of the form  a^T x <= b, x binary and a >= 0"
#define CONSHDLR_SEPAPRIORITY   +600000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   -600000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -600000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             0 /**< frequency for separating cuts; zero means to separate only in the root node */
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

#define LINCONSUPGD_PRIORITY    +100000 /**< priority of the constraint handler for upgrading of linear constraints */

#define MAX_DYNPROG_CAPACITY      10000 /**< maximal capacity of knapsack to apply dynamic programming */
#define MAX_USECLIQUES_SIZE        1000 /**< maximal number of items in knapsack where clique information is used */
#define MAX_ZEROITEMS_SIZE        10000 /**< maximal number of items to store in the zero list in preprocessing */

#define KNAPSACKRELAX_MAXDELTA        0.1 /**< maximal allowed rounding distance for scaling in knapsack relaxation */
#define KNAPSACKRELAX_MAXDNOM      1000LL /**< maximal allowed denominator in knapsack rational relaxation */
#define KNAPSACKRELAX_MAXSCALE     1000.0 /**< maximal allowed scaling factor in knapsack rational relaxation */

#define DEFAULT_SEPACARDFREQ          1 /**< multiplier on separation frequency, how often cardinality cuts are separated */
#define DEFAULT_MAXROUNDS             5 /**< maximal number of separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        -1 /**< maximal number of separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS          50 /**< maximal number of cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT     200 /**< maximal number of cuts separated per separation round in the root node */
#define DEFAULT_MAXNUMCARDLIFT       -1 /**< maximal number of cardinality inequalities lifted per sepa round (-1: unlimited) */
#define DEFAULT_MAXCARDBOUNDDIST    0.0 /**< maximal relative distance from current node's dual bound to primal bound compared
                                         *   to best node's dual bound for separating knapsack cardinality cuts */
#define DEFAULT_DISAGGREGATION     TRUE /**< should disaggregation of knapsack constraints be allowed in preprocessing? */
#define DEFAULT_SIMPLIFYINEQUALITIES FALSE/**< should presolving try to simplify knapsacks */

#define MAXABSVBCOEF               1e+5 /**< maximal absolute coefficient in variable bounds used for knapsack relaxation */



/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   int*                  ints1;              /**< cleared memory array, all entries are set to zero in initpre, if you use this
                                              *   you have to clear it at the end, exists only in presolving stage */
   int*                  ints2;              /**< cleared memory array, all entries are set to zero in initpre, if you use this
                                              *   you have to clear it at the end, exists only in presolving stage */
   SCIP_Longint*         longints1;          /**< cleared memory array, all entries are set to zero in initpre, if you use this
                                              *   you have to clear it at the end, exists only in presolving stage */
   SCIP_Longint*         longints2;          /**< cleared memory array, all entries are set to zero in initpre, if you use this
                                              *   you have to clear it at the end, exists only in presolving stage */
   SCIP_Bool*            bools1;             /**< cleared memory array, all entries are set to zero in initpre, if you use this
                                              *   you have to clear it at the end, exists only in presolving stage */
   SCIP_Bool*            bools2;             /**< cleared memory array, all entries are set to zero in initpre, if you use this
                                              *   you have to clear it at the end, exists only in presolving stage */
   SCIP_Bool*            bools3;             /**< cleared memory array, all entries are set to zero in initpre, if you use this
                                              *   you have to clear it at the end, exists only in presolving stage */
   SCIP_Bool*            bools4;             /**< cleared memory array, all entries are set to zero in initpre, if you use this
                                              *   you have to clear it at the end, exists only in presolving stage */
   SCIP_Real*            reals1;             /**< cleared memory array, all entries are set to zero in consinit, if you use this
                                              *   you have to clear it at the end */
   int                   ints1size;          /**< size of ints1 array */
   int                   ints2size;          /**< size of ints2 array */
   int                   longints1size;      /**< size of longints1 array */
   int                   longints2size;      /**< size of longints2 array */
   int                   bools1size;         /**< size of bools1 array */
   int                   bools2size;         /**< size of bools2 array */
   int                   bools3size;         /**< size of bools3 array */
   int                   bools4size;         /**< size of bools4 array */
   int                   reals1size;         /**< size of reals1 array */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
   SCIP_Real             maxcardbounddist;   /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for separating knapsack cardinality cuts */
   int                   sepacardfreq;       /**< multiplier on separation frequency, how often cardinality cuts are separated */
   int                   maxrounds;          /**< maximal number of separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of separation rounds in the root node (-1: unlimited) */
   int                   maxsepacuts;        /**< maximal number of cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of cuts separated per separation round in the root node */
   int                   maxnumcardlift;     /**< maximal number of cardinality inequ. lifted per sepa round (-1: unlimited) */
   SCIP_Bool             disaggregation;     /**< should disaggregation of knapsack constraints be allowed in preprocessing? */
   SCIP_Bool             simplifyinequalities;/**< should presolving try to cancel down or delete coefficients in inequalities */
};


/** constraint data for knapsack constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars;               /**< variables in knapsack constraint */
   SCIP_Longint*         weights;            /**< weights of variables in knapsack constraint */
   SCIP_EVENTDATA**      eventdatas;         /**< event datas for bound change events of the variables */
   int*                  cliquepartition;    /**< clique indices of the clique partition */
   SCIP_ROW*             row;                /**< corresponding LP row */
   int                   nvars;              /**< number of variables in knapsack constraint */
   int                   varssize;           /**< size of vars, weights, and eventdatas arrays */
   SCIP_Longint          capacity;           /**< capacity of knapsack */
   SCIP_Longint          weightsum;          /**< sum of all weights */
   SCIP_Longint          onesweightsum;      /**< sum of weights of variables fixed to one */
   unsigned int          propagated:1;       /**< is the knapsack constraint already propagated? */
   unsigned int          presolved:1;        /**< is the knapsack constraint already presolved? */
   unsigned int          sorted:1;           /**< are the knapsack items sorted by weight? */
   unsigned int          cliquepartitioned:1;/**< is the clique partition valid? */
   unsigned int          merged:1;           /**< are the constraint's equal variables already merged? */
   unsigned int          cliquesadded:1;     /**< were the cliques of the knapsack already added to clique table? */
};


/** event data for bound changes events */
struct SCIP_EventData
{
   SCIP_CONSDATA*        consdata;           /**< knapsack constraint data to process the bound change for */
   SCIP_Longint          weight;             /**< weight of variable */
   int                   filterpos;          /**< position of event in variable's event filter */
};




/*
 * Local methods
 */


/** creates event data */
static
SCIP_RETCODE eventdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTDATA**      eventdata,          /**< pointer to store event data */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Longint          weight              /**< weight of variable */
   )
{
   assert(eventdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, eventdata) );
   (*eventdata)->consdata = consdata;
   (*eventdata)->weight = weight;

   return SCIP_OKAY;
}  

/** frees event data */
static
SCIP_RETCODE eventdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTDATA**      eventdata           /**< pointer to event data */
   )
{
   assert(eventdata != NULL);

   SCIPfreeBlockMemory(scip, eventdata);

   return SCIP_OKAY;
}

/** sorts items in knapsack with nonincreasing weights */
static 
void sortItems(
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);
   assert(consdata->nvars == 0 || consdata->weights != NULL);
   assert(consdata->nvars == 0 || consdata->eventdatas != NULL);
   assert(consdata->nvars == 0 || consdata->cliquepartition != NULL);

   if( !consdata->sorted )
   {
      /* sort of four joint arrays of Long/pointer/pointer/ints, 
       * sorted by first array in non-increasing order via sort template */
      SCIPsortDownLongPtrPtrInt(
         consdata->weights, 
         (void**)consdata->vars, 
         (void**)consdata->eventdatas, 
         consdata->cliquepartition, 
         consdata->nvars);
      
      consdata->sorted = TRUE;
   }
#ifndef NDEBUG
   {
      /* check if the weight array is sorted in a non-increasing way */ 
      int i;
      for( i = 0; i < consdata->nvars-1; ++i )
         assert(consdata->weights[i] >= consdata->weights[i+1]);
   }
#endif
}

/** calculates a partition of the variables into cliques */
static
SCIP_RETCODE calcCliquepartition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->cliquepartition != NULL);

   if( !consdata->cliquepartitioned )
   {
      SCIP_CALL( SCIPcalcCliquePartition(scip, consdata->vars, consdata->nvars, consdata->cliquepartition) );
      consdata->cliquepartitioned = TRUE;
   }

   return SCIP_OKAY;
}

/** installs rounding locks for the given variable in the given knapsack constraint */
static
SCIP_RETCODE lockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding up may violate the constraint */
   SCIP_CALL( SCIPlockVarCons(scip, var, cons, FALSE, TRUE) );

   return SCIP_OKAY;
}

/** removes rounding locks for the given variable in the given knapsack constraint */
static
SCIP_RETCODE unlockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding up may violate the constraint */
   SCIP_CALL( SCIPunlockVarCons(scip, var, cons, FALSE, TRUE) );

   return SCIP_OKAY;
}

/** catches bound change events for variables in knapsack */
static
SCIP_RETCODE catchEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;
   
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);
   assert(consdata->nvars == 0 || consdata->weights != NULL);
   assert(consdata->nvars == 0 || consdata->eventdatas != NULL);

   for( i = 0; i < consdata->nvars; i++)
   {
      SCIP_CALL( eventdataCreate(scip, &consdata->eventdatas[i], consdata, consdata->weights[i]) );
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[i], 
            SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_UBRELAXED | SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_IMPLADDED,
            eventhdlr, consdata->eventdatas[i], &consdata->eventdatas[i]->filterpos) );
   }

   return SCIP_OKAY;
}

/** drops bound change events for variables in knapsack */
static
SCIP_RETCODE dropEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;
   
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);
   assert(consdata->nvars == 0 || consdata->weights != NULL);
   assert(consdata->nvars == 0 || consdata->eventdatas != NULL);

   for( i = 0; i < consdata->nvars; i++)
   {
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[i],
            SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_UBRELAXED | SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_IMPLADDED,
            eventhdlr, consdata->eventdatas[i], consdata->eventdatas[i]->filterpos) );
      SCIP_CALL( eventdataFree(scip, &consdata->eventdatas[i]) );
   }

   return SCIP_OKAY;
}

/** ensures, that vars and vals arrays can store at least num entries */
static
SCIP_RETCODE consdataEnsureVarsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< knapsack constraint data */
   int                   num,                /**< minimum number of entries to store */
   SCIP_Bool             transformed         /**< is constraint from transformed problem? */
   )
{
   assert(consdata != NULL);
   assert(consdata->nvars <= consdata->varssize);
   
   if( num > consdata->varssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, consdata->varssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->weights, consdata->varssize, newsize) );
      if( transformed )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->eventdatas, consdata->varssize, newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->cliquepartition, consdata->varssize, newsize) );
      }
      else
      {
         assert(consdata->eventdatas == NULL);
         assert(consdata->cliquepartition == NULL);
      }
      consdata->varssize = newsize;
   }
   assert(num <= consdata->varssize);

   return SCIP_OKAY;
}

/** creates knapsack constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   nvars,              /**< number of variables in knapsack */
   SCIP_VAR**            vars,               /**< variables of knapsack */
   SCIP_Longint*         weights,            /**< weights of knapsack items */
   SCIP_Longint          capacity            /**< capacity of knapsack */
   )
{
   int i;

   assert(consdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   if( nvars > 0 )
   {
      int k;
      int v;

      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->weights, weights, nvars) );
      
      k = 0;
      for( v = 0; v < nvars; ++v )
      {
         assert(vars[v] != NULL);
         
         /* all weight have to be not negative */
         assert( weights[v] >= 0 );

         if( weights[v] != 0 )
         {
            vars[k] = vars[v];
            weights[k] = weights[v];
            k++;
         }
      }
      (*consdata)->nvars = k;
   }
   else
   {
      (*consdata)->vars = NULL;
      (*consdata)->weights = NULL;
   }

   /* capacity has to be greater or equal to zero */
   assert( capacity >= 0 );

   (*consdata)->nvars = nvars;
   (*consdata)->varssize = nvars;
   (*consdata)->capacity = capacity;
   (*consdata)->eventdatas = NULL;
   (*consdata)->cliquepartition = NULL;
   (*consdata)->row = NULL;
   (*consdata)->weightsum = 0;
   (*consdata)->onesweightsum = 0;
   (*consdata)->propagated = FALSE;
   (*consdata)->presolved = FALSE;
   (*consdata)->sorted = FALSE;
   (*consdata)->cliquepartitioned = FALSE;
   (*consdata)->merged = FALSE;
   (*consdata)->cliquesadded = FALSE;

   /* get transformed variables, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );

      /* allocate memory for additional data structures */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->eventdatas, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->cliquepartition, nvars) );

      /* catch events for variables */
      SCIP_CALL( catchEvents(scip, *consdata, eventhdlr) );
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
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to the constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);

   if( (*consdata)->row != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->row) );
   }
   if( (*consdata)->eventdatas != NULL )
   {
      SCIP_CALL( dropEvents(scip, *consdata, eventhdlr) );
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->eventdatas, (*consdata)->varssize);
   }
   if( (*consdata)->cliquepartition != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->cliquepartition, (*consdata)->varssize);
   }
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->weights, (*consdata)->varssize);

   SCIPfreeBlockMemory(scip, consdata);
 
   return SCIP_OKAY;
}

/** changes a single weight in knapsack constraint data */
static
void consdataChgWeight(
   SCIP_CONSDATA*        consdata,           /**< knapsack constraint data */
   int                   item,               /**< item number */
   SCIP_Longint          newweight           /**< new weight of item */
   )
{
   SCIP_Longint oldweight;

   assert(consdata != NULL);
   assert(0 <= item && item < consdata->nvars);

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
   consdata->cliquepartitioned = FALSE; /* recalculate the clique partition after a weight was changed */
}

/** creates LP row corresponding to knapsack constraint */
static 
SCIP_RETCODE createRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< knapsack constraint */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);

   SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->row, SCIPconsGetName(cons),
         -SCIPinfinity(scip), (SCIP_Real)consdata->capacity,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );

   SCIP_CALL( SCIPcacheRowExtensions(scip, consdata->row) );
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, consdata->vars[i], (SCIP_Real)consdata->weights[i]) );
   }
   SCIP_CALL( SCIPflushRowExtensions(scip, consdata->row) );

   return SCIP_OKAY;
}  

/** adds linear relaxation of knapsack constraint to the LP */
static 
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_SOL*             sol                 /**< primal CIP solution, NULL for current LP solution */
   )
{
   SCIP_CONSDATA* consdata;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert(consdata->row != NULL);

   /* insert LP row as cut */
   if( !SCIProwIsInLP(consdata->row) )
   {
      SCIPdebugMessage("adding relaxation of knapsack constraint <%s> (capacity %"SCIP_LONGINT_FORMAT"): ", 
         SCIPconsGetName(cons), consdata->capacity);
      SCIPdebug( SCIProwPrint(consdata->row, NULL) );
      SCIP_CALL( SCIPaddCut(scip, sol, consdata->row, FALSE) );
   }

   return SCIP_OKAY;
}

/** checks knapsack constraint for feasibility of given solution: returns TRUE iff constraint is feasible */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SOL*             sol,                /**< solution to check, NULL for current solution */
   SCIP_Bool             checklprows,        /**< should LP rows be checked? */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   SCIP_CONSDATA* consdata;

   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("checking knapsack constraint <%s> for feasibility of solution %p (lprows=%u)\n",
      SCIPconsGetName(cons), (void*)sol, checklprows);

   *violated = FALSE;

   if( checklprows || consdata->row == NULL || !SCIProwIsInLP(consdata->row) )
   {
      SCIP_Real sum;
      int i;

      /* increase age of constraint; age is reset to zero, if a violation was found */
      SCIP_CALL( SCIPincConsAge(scip, cons) );

      sum = 0.0;
      for( i = 0; i < consdata->nvars && sum <= consdata->capacity + 0.1; ++i )
      {
         sum += consdata->weights[i] * SCIPgetSolVal(scip, sol, consdata->vars[i]);
      }

      if( SCIPisFeasGT(scip, sum, (SCIP_Real)consdata->capacity) )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         *violated = TRUE;

         if( printreason )
         {
            SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
            SCIPinfoMessage(scip, NULL, "violation: ");

            /* complete the activity computation to discover the violation */
            for( ; i < consdata->nvars; ++i )
            {
               sum += consdata->weights[i] * SCIPgetSolVal(scip, sol, consdata->vars[i]);
            }
            SCIPinfoMessage(scip, NULL, "violation: the capacity is violated by %.15g\n", sum - consdata->capacity);
         }
      }
   }
   
   return SCIP_OKAY;
}

/* IDX computes the integer index for the optimal solution array */
#define IDX(j,d) ((j)*(intcap+1)+(d))

/** solves knapsack problem in maximization form exactly using dynamic programming;
 *  if needed, one can provide arrays to store all selected items and all not selected items
 */
SCIP_RETCODE SCIPsolveKnapsackExactly(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nitems,             /**< number of available items */
   SCIP_Longint*         weights,            /**< item weights */
   SCIP_Real*            profits,            /**< item profits */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   int*                  items,              /**< item numbers */
   int*                  solitems,           /**< array to store items in solution, or NULL */
   int*                  nonsolitems,        /**< array to store items not in solution, or NULL */
   int*                  nsolitems,          /**< pointer to store number of items in solution, or NULL */
   int*                  nnonsolitems,       /**< pointer to store number of items not in solution, or NULL */
   SCIP_Real*            solval              /**< pointer to store optimal solution value, or NULL */
   ) 
{
   SCIP_Real* tempsort;
   SCIP_Real* optvalues;
   int intcap;
   int d;
   int j;
   int i;

   assert(weights != NULL);
   assert(profits != NULL);
   assert(capacity >= 0);
   assert(capacity < INT_MAX);
   assert(items != NULL);
   assert(nitems >= 0);

   intcap = (int)capacity;
   assert(intcap >= 0);
   SCIP_CALL( SCIPallocBufferArray(scip, &optvalues, (nitems+1)*(intcap+1)) );
   
   /* allocate memory for temporary array used for sorting; array should contain profits devided by corresponding weights (p_1 / w_1 ... p_n / w_n )*/
   SCIP_CALL( SCIPallocBufferArray(scip, &tempsort, nitems) );
   /* initialize temporary array */ 
   for (i = nitems - 1; i >= 0; --i)
   {
      tempsort[i] = profits[i] / weights [i];
   }

   /* sort tempsort, items, weights and profits such that p_1 / w_1 >= p_2 / w_2 >= ... >= p_n / w_n */
   SCIPsortDownRealLongRealInt( tempsort, weights, profits, items, nitems);

   /* free temporary array */
   SCIPfreeBufferArray(scip, &tempsort);

   /* fills dynamic programming table with optimal values */
   for( d = 0; d <= intcap; d++ )
      optvalues[IDX(0,d)] = 0.0;
   for( j = 1; j <= nitems; j++ )
   {
      int intweight;

      assert(0 <= weights[j-1] && weights[j-1] < SCIP_LONGINT_MAX);
      if( weights[j-1] >= INT_MAX )
      {
         assert(weights[j-1] > capacity);
         intweight = INT_MAX;
      }
      else
         intweight = (int)weights[j-1];
      assert(intweight >= 0);

      for( d = 0; d < intweight && d <= intcap; d++ )
         optvalues[IDX(j,d)] = optvalues[IDX(j-1,d)];
      for( d = intweight; d <= intcap; d++ )
      {
         SCIP_Real sumprofit;

         sumprofit = optvalues[IDX(j-1,d-intweight)] + profits[j-1];
         optvalues[IDX(j,d)] = MAX(sumprofit, optvalues[IDX(j-1,d)]);
      } 
   }

   /* produces optimal solution by following the table */
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

/**todo: use SCIPsolveKnapsackExactly() instead of SCIPsolveKnapsack() in all methods of cons_knapsack */

/** solves knapsack problem with dynamic programming;
 *  if needed, one can provide arrays to store all selected items and all not selected items
 */
SCIP_RETCODE SCIPsolveKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nitems,             /**< number of available items */
   SCIP_Longint*         weights,            /**< item weights */
   SCIP_Real*            profits,            /**< item profits */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   int*                  items,              /**< item numbers, or NULL */
   int*                  solitems,           /**< array to store items in solution, or NULL */
   int*                  nonsolitems,        /**< array to store items not in solution, or NULL */
   int*                  nsolitems,          /**< pointer to store number of items in solution, or NULL */
   int*                  nnonsolitems,       /**< pointer to store number of items not in solution, or NULL */
   SCIP_Real*            solval              /**< pointer to store optimal solution value, or NULL */
   ) 
{
   SCIP_Real* optvalues;
   int intcap;
   int d;
   int j;

   assert(weights != NULL);
   assert(profits != NULL);
   assert(capacity >= 0);
   assert(capacity < INT_MAX);
   assert(nitems >= 0);

   if( capacity == 0 )
   {
#ifndef NDEBUG

      for( j = nitems - 1; j >= 0; --j )
         assert(weights[j] > 0);
#endif
      if( solitems != NULL)
      {
         assert(nonsolitems != NULL);
         assert(nnonsolitems != NULL);
         assert(items != NULL);
         *nsolitems = 0;
         
         BMScopyMemoryArray(nonsolitems, items, nitems);
         *nnonsolitems = nitems;
      }
      if( solval != NULL )
         *solval = 0.0;
      
      return SCIP_OKAY;
   }

   intcap = (int)capacity;
   assert(intcap >= 0);
   SCIP_CALL( SCIPallocBufferArray(scip, &optvalues, (nitems+1)*(intcap+1)) );
   
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
         SCIP_Real sumprofit;

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

/** gets a most violated minimal cover C = C2 & C1 for a given knapsack constraint, considering that all variables 
 *  in a given set C2 of variables in knapsack constraint are fixed to one (variables for downlifting)
 */
static
SCIP_RETCODE getCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  varsc2,             /**< set C2 of variables in knapsack constraint (variables for downlifting) */
   int*                  varscn2,            /**< set of variables in knapsack constraint that are not in C2 (cand. for C1) */
   int                   nvarsc2,            /**< number of variables in C2 */
   int                   nvarscn2,           /**< number of variables not in C2 */
   SCIP_Longint          varsc2weight,       /**< sum of weights of variables in C2 */
   int*                  covervars,          /**< pointer to store cover variables C = C2 & C1 */
   int*                  noncovervars,       /**< pointer to store noncover variables */
   int*                  ncovervars,         /**< pointer to store number of cover variables */
   int*                  ncovervarsc1,       /**< pointer to store number of cover variables in C1 (at the end of covervars) */
   int*                  ncovervarsc2,       /**< pointer to store number of cover variables in C2 (at the beg of covervars) */
   int*                  nnoncovervars,      /**< pointer to store number of noncover variables */
   SCIP_Longint*         coverweight,        /**< pointer to store weight of cover */
   SCIP_Real*            covervarsc1activity,/**< pointer to store activity of cover variables in C1 */  
   SCIP_Bool*            found               /**< pointer to store whether a cover was found */
   )
{
   SCIP_Real* dpprofits;
   SCIP_Real relslack;
   SCIP_Longint* dpweights;
   SCIP_Longint dpcapacity;
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
   
   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &items, nvarscn2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &dpweights, nvarscn2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &dpprofits, nvarscn2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedzeros, nvarscn2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &setvars, nvarscn2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nonsetvars, nvarscn2) );

   *found = FALSE;
   *ncovervars = 0;
   *ncovervarsc1 = 0;
   *ncovervarsc2 = 0;
   *nnoncovervars = 0;
   *coverweight = 0;
   *covervarsc1activity = 0.0;
   nsetvars = 0;
   nnonsetvars = 0;
      
   /* solve the following knapsack problem with dynamic programming or greedy heuristic:
    * max SUM (1 - x_i^*) z_i
    *     SUM w_i z_i <= SUM w_i - capacity - 1
    *             z_i == 0                       forall variables x_i in C2
    * with the meaning: x_i is member of cover <==> z_i == 0
    */
   
   /* use all variables with fractional solution value from the set of variables not in C2 */
   nitems = 0;
   nfixedzeros = 0;
   dpcapacity = varsc2weight - capacity - 1;
   for( i = 0; i < nvarscn2; i++ )
   {
      assert(SCIPvarGetType(vars[varscn2[i]]) == SCIP_VARTYPE_BINARY);
      
      /* variable has fractional solution value */
      if( SCIPisFeasGT(scip, solvals[varscn2[i]], 0.0) )
      {
         /* sort items by non-decreasing relative slack value (1-x_i^*)/w_i */
         relslack = (1.0 - solvals[varscn2[i]])/weights[varscn2[i]];
         for( j = nitems; j > 0 && relslack < (1.0 - solvals[items[j-1]])/weights[items[j-1]]; --j )
            items[j] = items[j-1];
         items[j] = varscn2[i];
         nitems++;
         dpcapacity += weights[varscn2[i]];
      }
      /* variable has solution value 0 */
      else
      {
         fixedzeros[nfixedzeros] = varscn2[i];
         nfixedzeros++;
      }
   }
   assert(nitems + nfixedzeros == nvarscn2);

   /* get weights and slacks of items */
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
         /* solve separation knapsack with dynamic programming */
         SCIP_CALL( SCIPsolveKnapsack(scip, nitems, dpweights, dpprofits, dpcapacity, 
               items, nonsetvars, setvars, &nnonsetvars, &nsetvars, NULL) ); 
         assert(nsetvars + nnonsetvars == nitems);
      }
      else
      {
         SCIP_Longint setweight;

         /* use greedy heuristic on separation knapsack */
         nnonsetvars = 0;
         nsetvars = 0;
         setweight = 0;
         for( i = nitems-1; i >= 0; --i )
         {
            assert(weights[items[i]] == dpweights[i]);
            if( dpcapacity - dpweights[i] >= 0 )
            {
               dpcapacity -= dpweights[i];
               nonsetvars[nnonsetvars] = items[i];
               nnonsetvars++;
            }
            else
            {
               setweight += dpweights[i];
               setvars[nsetvars] = items[i];
               nsetvars++;
            }
         }

         /* the cover might be empty, which means the variables in C2 alone already suffice to exceed the capacity;
          * this can happen even if C2 consists of the variables at 1.0 in the LP relaxation, because the linear
          * constraint handler replaces continuous and integer variables by variable lower/upper bounds, and the LP
          * relaxation may violate these variable bounds
          */
         if( nsetvars > 0 )
         {
            SCIP_Longint minimumweight;
            int minimumweightidx;
            int removevar;
            SCIP_Bool minimal;
            /* get minimum weigth of set variables */
            minimal = FALSE;
            minimumweight = SCIP_LONGINT_MAX;
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

            /* test if set is allready minimal */
            if( setweight - minimumweight <= capacity - varsc2weight )
               minimal = TRUE;

            /* make set minimal by removing variables (in decreasing order of slack) */
            for( i = 0; i < nsetvars && !minimal; i++ )
            {
               /* set variable x_i can not be removed with respect to the set property */
               if( setweight - weights[setvars[i]] <= capacity - varsc2weight )
                  continue;

               /* remove x_i from setvars */
               removevar = setvars[i];
               for( j = i; j < nsetvars - 1; j++ )
                  setvars[j] = setvars[j + 1];
               nsetvars--;
               setweight -= weights[removevar];
            
               /* add x_i to nonsetvars */
               nonsetvars[nnonsetvars] = removevar;
               nnonsetvars++;
            
               /* update minimumweight of setvars */     
               if( minimumweightidx == i )
               {
                  /* get minimum weigth of setvars */
                  minimumweight = SCIP_LONGINT_MAX;
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
         
               /* update index of current setvar */
               i--;

               /* test if set is now minimal */
               assert(setweight > capacity - varsc2weight);
               assert(!minimal);
               if( setweight - minimumweight <= capacity - varsc2weight )
                  minimal = TRUE;
            }
            assert(minimal);
         }
      }
      assert(*ncovervars == 0);
      assert(*nnoncovervars == 0);
      assert(nvarsc2 + nsetvars + nnonsetvars + nfixedzeros == nvars);

      /* get covervars (variables from C2 are in the beginning of covervars) and noncovervars: */ 

      /* add all variables from C2 to covervars */
      for( i = 0; i < nvarsc2; i++)
      {
         covervars[*ncovervars] = varsc2[i];
         *coverweight += weights[varsc2[i]];
         (*ncovervars)++;
         (*ncovervarsc2)++;
      }

      /* add all variables from setvars (variables from set of variables not in C2 which have been chosen to be in the 
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

      /* add all variables from nonsetvars (variables from set of variables not in C2 which have not been chosen to be 
       * in the cover) to noncovervars 
       */ 
      for( i = 0; i < nnonsetvars; i++ )
      {
         noncovervars[*nnoncovervars] = nonsetvars[i];
         (*nnoncovervars)++;
      }  

      /* add all variables with solution value zero to noncovervars */
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
      /* put all variables into noncovervars because no cover was found (variables from fixedzeros have not been used) */
      
      /* add all variables from C2 to noncovervars */
      for( i = 0; i < nvarsc2; i++)
      {
         noncovervars[*nnoncovervars] = varsc2[i];
         (*nnoncovervars)++;
      }
      
      /* add all variables from set of variables not in C2 to noncovervars */
      for( i = 0; i < nvarscn2; i++)
      {
         noncovervars[*nnoncovervars] = varscn2[i];
         (*nnoncovervars)++;
      }
      assert(!(*found));
      assert(*coverweight == 0 && *covervarsc1activity == 0.0); 
   }

   /* free temporary memory */
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
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables in knapsack constraint */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  varsc2,             /**< pointer to store variables of knapsack constraint in C2 */
   int*                  varscn2,            /**< pointer to store variables of knapsack constraint not in C2 */
   int*                  nvarsc2,            /**< pointer to store number of variables of knapsack constraint in C2 */
   int*                  nvarscn2,           /**< pointer to store number of variables of knapsack constraint not in C2 */
   SCIP_Longint*         varsc2weight        /**< pointer to store sum of weights of variables of knapsack constraint in C2 */
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

   /* choses all variables with solution value equal to one for the set of variables for downlifting (C2) */
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

/** sorts given set of variables in knapsack constraint by non-increasing solution value (variables with equal solution value 
 *  are ordered by non-increasing weight)
 */
static
void sortSetvars(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  setvars,            /**< set of variables in knapsack constraint to be sorted */
   int                   nsetvars,           /**< number of variables in set of variables in knapsack cons. to be sorted */
   SCIP_Real*            solvals             /**< solution values of all problem variables */
   )
{
   SCIP_Real solval;
   int idx;
   int i;
   int j;

   assert(scip != NULL);
   assert(setvars != NULL);
   assert(nsetvars >= 0);
   assert(solvals != NULL);

   /* sorts setvars by non-increasing solution value */
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
SCIP_RETCODE enlargeMinweights(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Longint**        minweightsptr,      /**< pointer to minweights table */
   int*                  minweightslen,      /**< pointer to store number of entries in minweights table (incl. z=0) */
   int*                  minweightssize,     /**< pointer to current size of minweights table */
   int                   newlen              /**< new length of minweights table */
   )
{
   int j;

   assert(minweightsptr != NULL);
   assert(*minweightsptr != NULL);

   if( newlen > *minweightssize )
   {
      int newsize;

      /* reallocate table memory */
      newsize = SCIPcalcMemGrowSize(scip, newlen);
      SCIP_CALL( SCIPreallocBufferArray(scip, minweightsptr, newsize) );
      *minweightssize = newsize;
   }
   assert(newlen <= *minweightssize);

   /* initialize new elements */
   for( j = *minweightslen; j < newlen; ++j )
      (*minweightsptr)[j] = SCIP_LONGINT_MAX;
   *minweightslen = newlen;

   return SCIP_OKAY;
}

/** lifts up given inequality sum(j in C1) x_j <= liftrhs which is valid for the knapsack polytop { x binary | 
 *  x is solution of knapsack constraint, x_j = 1 for all j in C2, x_j = 0 for all j in noncovervars } to a valid 
 *  inequality for the knapsack polytop { x binary | x is solution of knapsack constraint, x_j = 1 for all j in C2,
 *  x_j = 0 for all j in noncovervars with solution value = 0 }
 */
static
SCIP_RETCODE liftupKnapsackCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  covervars,          /**< cover variables C = C2 & C1 (C2, C1 sorted by non-incr LP val then weight) */
   int*                  noncovervars,       /**< noncover variables (sorted by non-incr LP val then weight) */
   int                   ncovervars,         /**< number of cover variables */
   int                   ncovervarsc1,       /**< number of cover variables in C1 (at the end of covervars) */
   int                   ncovervarsc2,       /**< number of cover variables in C2 (at the beginning of covervars) */
   int                   nnoncovervars,      /**< number of noncover variables */
   SCIP_Longint          fixedoneweight,     /**< sum of weights of variables from C2 which are still fixed to 1 */ 
   int*                  liftcoefs,          /**< lifting coefficients of var in knapsack constraint in lifted inequality */
   int*                  liftrhs,            /**< right hand side of lifted inequality */
   SCIP_Real*            liftlpval,          /**< LP solution value of lifted variables (without C1) */  
   int*                  lastlifted,         /**< pointer to store index of noncover var lifted last (-1 if no var lifted) */
   SCIP_Longint*         initminweight,      /**< initial minweight table (aggregation of sorted weights) */
   SCIP_Longint**        minweightsptr,      /**< pointer to minweights table */
   int*                  minweightslen,      /**< pointer to store number of entries in minweights table (incl. z=0) */
   int*                  minweightssize      /**< pointer to current size of minweights table */
   )
{
   SCIP_Longint* minweights;
   SCIP_Longint rescapacity;
   int i;

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
   assert(minweightsptr != NULL);
   assert(minweightslen != NULL);
   assert(minweightssize != NULL);
   assert(ncovervarsc1 < *minweightssize);
   assert(initminweight[0] == 0);
   assert(*liftrhs <= ncovervarsc1);

   minweights = *minweightsptr;
   assert(minweights != NULL);

   /* initialize the minweight table
    *  minweights[z] := minimal sum of weights s.t. activity of cut inequality equals z
    */
   BMScopyMemoryArray(minweights, initminweight, ncovervarsc1+1);
   *minweightslen = ncovervarsc1 + 1;
   assert(*liftrhs < *minweightslen);
   assert(0 <= *minweightslen && *minweightslen <= *minweightssize);

   /* calculate lifting coefficients for noncover variables with solution value > 0:
    *  for each noncovervar x_i with LP-value > 0 (in non-increasing order of LP-value and weight): 
    *   1. calculates max activity z_max of current lifted inequality s.t. the knapsack is still feasible with 
    *      x_i = 1, x_j = 0 for all j in noncovervars with j > i, x_j = 1 for all j in C2
    *   2. adds x_i with lifting coefficient beta_i = liftrhs - z_max to current lifted inequality
    *   3. updates minweight table: calculates minimal sum of weights s.t. activity of current lifted inequality equals z 
    */
   for( i = 0; i < nnoncovervars && SCIPisFeasPositive(scip, solvals[noncovervars[i]]); i++ )
   {
      SCIP_Longint weight;
      int liftvar;
      int liftcoef;
      int z;

      liftvar = noncovervars[i];
      assert(liftvar >= 0);
      weight = weights[liftvar];

      /* binary search in sorted minweights array for the largest entry that is not greater than 
       * capacity - fixedoneweight - weight_i 
       */
      assert(*liftrhs < *minweightslen);
      rescapacity = capacity - fixedoneweight - weight;
      if( rescapacity < 0 )
      {
         /* lifted inequality cut is valid for any beta_i, but it may no longer be facet inducing: set beta = liftrhs+1 */
         liftcoef = (*liftrhs) + 1;
      }
      else if( minweights[*liftrhs] <= rescapacity )
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
            assert(0 <= middle && middle < *minweightslen);
            if( minweights[middle] <= rescapacity ) 
               left = middle;
            else
               right = middle;
         }
         assert(left == right - 1);
         assert(0 <= left && left < *minweightslen);
         assert(minweights[left] <= rescapacity);
         assert(left == (*minweightslen)-1 || minweights[left+1] > rescapacity);
            
         /* now z_max = left: calculate lifting coefficient beta_i */
         liftcoef = (*liftrhs) - left;
      }
      liftcoefs[liftvar] = liftcoef;

      /* update LP solution value of lifted variables (without C1) */
      *liftlpval += liftcoef * solvals[liftvar];
         
      /* if we have downlifting candidates, we have to increase the size of the minweight table;
       * otherwise, we only have to keep track of the elements up to z = liftrhs
       */
      if( ncovervarsc2 > 0 )
      {
         SCIP_CALL( enlargeMinweights(scip, minweightsptr, minweightslen, minweightssize, *minweightslen + liftcoef) );
         minweights = *minweightsptr;
      }

      /* update minweights table: calculates minimal sum of weights s.t. activity of curr lifted inequality equals z */
      for( z = (*minweightslen) - 1; z >= liftcoef; z-- )
      {
         SCIP_Longint min;
         min = MIN(minweights[z], minweights[z - liftcoef] + weights[liftvar]);
         assert(min > 0);
         minweights[z] = min;
      }
   }

   /* update noncover variable lifted last */
   *lastlifted = i-1;

   return SCIP_OKAY;
}

/** lifts down given inequality sum(j in C1 & noncovervars with solution value > 0) beta_j * x_j <= liftrhs which is valid for 
 *  the knapsack polytop { x binary | x is solution of knapsack constraint, x_j = 1 for all j in C2, x_j = 0 for all j in 
 *  noncovervars with solution value 0} to a valid inequality for the knapsack polytop { x binary | x is solution of knapsack 
 *  constraint, x_j = 0 for all j in noncovervars with solution value 0}
 */
static
SCIP_RETCODE liftdownKnapsackCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  covervars,          /**< cover variables C = C2 & C1 (C2, C1 sorted by non-incr LP val then weight) */
   int*                  noncovervars,       /**< noncover variables (sorted by non-incr LP val then weight) */
   int                   ncovervars,         /**< number of cover variables */
   int                   ncovervarsc1,       /**< number of cover variables in C1 (at the end of covervars) */
   int                   ncovervarsc2,       /**< number of cover variables in C2 (at the beginning of covervars) */
   int                   nnoncovervars,      /**< number of noncover variables */
   SCIP_Longint          fixedoneweight,     /**< sum of weights of variables from C2 which are still fixed to 1 */ 
   int*                  liftcoefs,          /**< lifting coefficients of var in knapsack constraint in lifted inequality */
   int*                  liftrhs,            /**< right hand side of lifted inequality */
   SCIP_Real*            liftlpval,          /**< LP solution value of lifted variables (without C1) */  
   SCIP_Longint**        minweightsptr,      /**< pointer to minweights table */
   int*                  minweightslen,      /**< pointer to store number of entries in minweights table (incl. z=0) */
   int*                  minweightssize       /**< pointer to current size of minweights table */
   )
{
   SCIP_Longint* minweights;
   int i;

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
   assert(minweightsptr != NULL);
   assert(minweightslen != NULL);
   assert(0 <= *minweightslen && *minweightslen <= *minweightssize);

   minweights = *minweightsptr;
   assert(minweights != NULL);

   /* calculate lifting coefficients for cover variables in C2:
    *  for each covervar x_i in C2 (in non-increasing order of LP-value and weight): 
    *   1. calculate max activity z_max of current lifted inequality s.t. the knapsack is still feasible with 
    *      x_i = 0, x_j = 1 for all covervars in C2 with j > i, x_j = 0 for all j in noncovervars with solution value 0
    *   2. add x_i with lifting coefficient beta_i = z_max - liftrhs to current lifted inequality
    *   3. change liftrhs: liftrhs = liftrhs + beta_i
    *   4. update minweight table: calculate minimal sum of weights s.t. activity of current lifted inequality equals z 
    */
   for( i = 0; i < ncovervarsc2; i++ )
   {
      SCIP_Longint rescapacity;
      int liftvar;
      int liftcoef;
      int left;
      int right;
      int middle;
      int z;

      assert(SCIPisFeasEQ(scip, solvals[covervars[i]], 1.0)); 
      liftvar = covervars[i];
      assert(liftvar >= 0);
      
      /* update sum of weights of variables in knapsack still fixed to 1 */
      fixedoneweight -= weights[liftvar];
    
      /* binary search in sorted minweights array for the largest entry that is not greater then 
       * capacity - fixedoneweight
       */
      rescapacity = capacity - fixedoneweight;
      assert(rescapacity >= 0);
      left = 0;
      right = *minweightslen;
      while( left < right - 1)
      {
         middle = (left + right) / 2;
         assert(0 <= middle && middle < *minweightslen);
         if( minweights[middle] <= rescapacity ) 
            left = middle;
         else
            right = middle;
      }
      assert(left == right - 1);
      assert(0 <= left && left < *minweightslen);
      assert(minweights[left] <= rescapacity);
      assert(left == *minweightslen - 1 || minweights[left+1] > rescapacity);

      /* now z_max = left: calculate lifting coefficient beta_i */
      assert(left >= (*liftrhs));
      liftcoef = left - (*liftrhs);
      liftcoefs[liftvar] = liftcoef;
      if( liftcoef == 0 )
         continue;

      /* update LP solution value of lifted variables (without C1) */
      *liftlpval += liftcoef * solvals[liftvar];
      
      /* change liftrhs */
      (*liftrhs) += liftcoef;
      
      /* update minweights table: calculate minimal sum of weights s.t. activity of curr lifted inequality equals z */
      SCIP_CALL( enlargeMinweights(scip, minweightsptr, minweightslen, minweightssize, *minweightslen + liftcoef) );
      minweights = *minweightsptr;
      for( z = *minweightslen - 1; z >= liftcoef; z-- )
      {
         SCIP_Longint min;
         min = MIN(minweights[z], minweights[z - liftcoef] + weights[liftvar]);
         assert(min > 0);
         minweights[z] = min;
      }
   }

   return SCIP_OKAY;
}

/** lifts up given inequality sum(j in C1 & noncovervars with solution value > 0 & C2) beta_j * x_j <= liftrhs which is valid 
 *  for the knapsack polytop { x binary | x is solution of knapsack constraint, x_j = 0 for all j in noncovervars with 
 *  solution value 0} to a valid inequality for the knapsack polytop { x binary | x is solution of knapsack constraint }
 */
static
void liftupZerosKnapsackCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  covervars,          /**< cover variables C = C2 & C1 (C2, C1 sorted by non-incr LP val then weight) */
   int*                  noncovervars,       /**< noncover variables (sorted by non-incr solution value then weight) */
   int                   ncovervars,         /**< number of cover variables */
   int                   ncovervarsc1,       /**< number of cover variables in C1 (at the end of covervars) */
   int                   ncovervarsc2,       /**< number of cover variables in C2 (at the beginning of covervars) */
   int                   nnoncovervars,      /**< number of noncover variables */
   int*                  liftcoefs,          /**< lifting coefficients of variables in knapsack cons in lifted inequality */
   int*                  liftrhs,            /**< right hand side of lifted inequality */
   SCIP_Real*            liftlpval,          /**< LP solution value of lifted variables (without C1) */  
   int                   lastlifted,         /**< index of last noncover var with solution value > 0 lifted (-1 if no var lifted) */
   SCIP_Longint*         minweights          /**< minweights table */
   )
{
   int i;

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
   assert(minweights != NULL);

   /* calculates lifting coefficients for noncover variables with solution value = 0:
    *  for each noncovervar x_i with LP-value = 0 (i > lastlifted) (in non-increasing order of LP-value and weight): 
    *   1. calculates max activity z_max of current lifted inequality s.t. the knapsack is still feasible with 
    *      x_i = 1, x_j = 0 for all j in noncovervars with j > i
    *   2. adds x_i with lifting coefficient beta_i = liftrhs - z_max to current lifted inequality
    *   3. updates minweight table: calculates minimal sum of weights s.t. activity of current lifted inequality equals z 
    */
   for( i = lastlifted + 1; i < nnoncovervars; i++ )
   {
      SCIP_Longint weight;
      SCIP_Longint rescapacity;
      int liftvar;
      int liftcoef;
      int z;

      assert(SCIPisFeasLE(scip, solvals[noncovervars[i]], 0.0));  

      liftvar = noncovervars[i];
      assert(liftvar >= 0);
      weight = weights[liftvar];
 
      /* binary search in sorted minweights array for the largest entry that is not greater then capacity - weight_i */
      rescapacity = capacity - weight;
      if( rescapacity < 0 )
      {
         /* lifted inequality cut is valid for any beta_i, but it may no longer be facet inducing, sets beta_i = liftrhs+1 */
         liftcoef = (*liftrhs) + 1;
      }
      else if( minweights[*liftrhs] <= rescapacity )
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
            if( minweights[middle] <= rescapacity ) 
               left = middle;
            else
               right = middle;
         }
         assert(left == right - 1);
         assert(minweights[left] <= rescapacity);
         assert(left == *liftrhs || minweights[left+1] > rescapacity);
            
         /* now z_max = left: calculates lifting coefficient beta_i */
         liftcoef = (*liftrhs) - left;
      }
      liftcoefs[liftvar] = liftcoef;

      /* updates LP solution value of lifted variables (without C1) */
      *liftlpval += liftcoef * solvals[liftvar];
         
      /* update minweights table: calculates minimal sum of weights s.t. activity of curr lifted inequality equals z */
      for( z = *liftrhs; z >= liftcoef; z-- )
      {
         SCIP_Longint min;
         min = MIN(minweights[z], minweights[z - liftcoef] + weights[liftvar]);
         assert(min > 0);
         minweights[z] = min;
      }
   }
}

/** lifts given cardinality inequality sum(j in C1) x_j <= |C1| to a valid inequality of the full dimensional knapsack 
 *  polytop by using uplifting for all variables not in the cover and downlifting for all variables in the cover that 
 *  are fixed to one (C2)
 */
static
SCIP_RETCODE liftKnapsackCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  covervars,          /**< cover variables C = C2 & C1 (C2, C1 sorted by non-incr LP val then weight) */
   int*                  noncovervars,       /**< noncover variables (sorted by non-incr LP val then weight) */
   int                   ncovervars,         /**< number of cover variables */
   int                   ncovervarsc1,       /**< number of cover variables in C1 (at the end of covervars) */
   int                   ncovervarsc2,       /**< number of cover variables in C2 (at the beginning of covervars) */
   int                   nnoncovervars,      /**< number of noncover variables */
   SCIP_Longint*         initminweight,      /**< initial minweight table (aggregation of sorted weights) */
   int*                  liftcoefs,          /**< pointer to store lifting coefficient of variables in knapsack constraint */
   int*                  liftrhs,            /**< pointer to store right hand side of the lifted cover inequality */
   SCIP_Real*            liftlpval           /**< pointer to store LP solution value of lifted variables */  
   )
{
   SCIP_Longint fixedoneweight;
   SCIP_Longint* minweights;
   int minweightssize;
   int minweightslen;
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
   minweightssize = ncovervarsc1+1;
   SCIP_CALL( SCIPallocBufferArray(scip, &minweights, minweightssize) );

   /* initializes data structures */
   *liftlpval = 0.0;
   BMSclearMemoryArray(liftcoefs, nvars);

   /* sets lifting coefficient of cover variables of C1 */
   for( i = ncovervarsc2; i < ncovervars; i++ )
   {
      assert(liftcoefs[covervars[i]] == 0);
      liftcoefs[covervars[i]] = 1;
   }
   
   /* lifts up all noncover variables with solution value > 0 */
   fixedoneweight = 0;
   for( i = 0; i < ncovervarsc2; i++ )
      fixedoneweight += weights[covervars[i]];
   SCIP_CALL( liftupKnapsackCover(scip, vars, weights, capacity, solvals, covervars, noncovervars, ncovervars, 
         ncovervarsc1, ncovervarsc2, nnoncovervars, fixedoneweight, liftcoefs, liftrhs, liftlpval, &lastlifted, 
         initminweight, &minweights, &minweightslen, &minweightssize) );
   assert(minweightslen <= minweightssize);
   assert(lastlifted < nnoncovervars);

   /* lifts down all cover variables from C2 */
   SCIP_CALL( liftdownKnapsackCover(scip, vars, weights, capacity, solvals, covervars, noncovervars, ncovervars,
         ncovervarsc1, ncovervarsc2, nnoncovervars, fixedoneweight, liftcoefs, liftrhs, liftlpval,
         &minweights, &minweightslen, &minweightssize) );

   /* lifts up all remaining noncover variables (noncover variables with solution value 0) */
   liftupZerosKnapsackCover(scip, vars, weights, capacity, solvals, covervars, noncovervars, ncovervars, ncovervarsc1, 
      ncovervarsc2, nnoncovervars, liftcoefs, liftrhs, liftlpval, lastlifted, minweights);

   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &minweights);

   return SCIP_OKAY;
}

/** lifts given cardinality inequality sum(j in C1) x_j <= |C1| to a valid inequality of the full dimensional knapsack 
 *  polytop by using uplifting for all variables not in the cover and downlifting for all variables in the cover that 
 *  are fixed to one (C2)
 */
SCIP_RETCODE SCIPliftKnapsackCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  covervars,          /**< cover variables C = C2 & C1 (C2, C1 sorted by non-incr LP val then weight) */
   int*                  noncovervars,       /**< noncover variables (sorted by non-incr LP val then weight) */
   int                   ncovervars,         /**< number of cover variables */
   int                   ncovervarsc1,       /**< number of cover variables in C1 (at the end of covervars) */
   int                   ncovervarsc2,       /**< number of cover variables in C2 (at the beginning of covervars) */
   int                   nnoncovervars,      /**< number of noncover variables */
   int*                  liftcoefs,          /**< pointer to store lifting coefficient of variables in knapsack constraint */
   int*                  liftrhs,            /**< pointer to store right hand side of the lifted cover inequality */
   SCIP_Real*            liftlpval           /**< pointer to store LP solution value of lifted variables */  
   )
{
   SCIP_Longint* initminweight;
   int i;

   assert(weights != NULL);
   assert(covervars != NULL);

   /* allocate temporary memory for initial minweight table */
   SCIP_CALL( SCIPallocBufferArray(scip, &initminweight, ncovervarsc1+1) );
   
   /* sort weights of variables in C1 (leave open position 0 to store a zero afterwards) */
   for( i = 0; i < ncovervarsc1; i++ )
   {
      SCIP_Longint weight;
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
   SCIP_CALL( liftKnapsackCover(scip, vars, nvars, weights, capacity, solvals, covervars, noncovervars,
         ncovervars, ncovervarsc1, ncovervarsc2, nnoncovervars, initminweight, liftcoefs, liftrhs, liftlpval) );

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &initminweight);

   return SCIP_OKAY;
}

/** separates lifted cover inequalities for given knapsack problem */
SCIP_RETCODE SCIPseparateKnapsackCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that originates the knapsack problem, or NULL */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_SOL*             sol,                /**< primal CIP solution to separate, NULL for current LP solution */
   int                   maxnumcardlift,     /**< maximal number of cardinality inequ. lifted per sepa round (-1: unlimited) */
   int*                  ncuts               /**< pointer to add up the number of found cuts */
   )
{
   SCIP_Real* solvals;
   SCIP_Real solval;
   SCIP_Real covervarsc1activity;
   SCIP_Real liftlpval;
   SCIP_Bool coverfound;
   SCIP_Longint varsc2weight;
   SCIP_Longint coverweight;
   SCIP_Longint weight;
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
   assert(vars != NULL);
   assert(nvars > 0);
   assert(weights != NULL);
   assert(capacity >= 0);
   assert(ncuts != NULL);

   /* increase age of constraint (age is reset to zero, if a cut was found) */
   if( cons != NULL )
   {
      SCIP_CALL( SCIPincConsAge(scip, cons) );
   }

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varscn2, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsc2, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &covervars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &noncovervars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &liftcoefs, nvars) );

   /* get solution values of all problem variables */
   SCIP_CALL( SCIPgetSolVals(scip, sol, nvars, vars, solvals) );

   /* get a partition of the set of variables in knapsack constraint in a set of variables for downlifting (varsc2) and 
    * a set of the remaining variables (varscn2) 
    */
   getPartition(scip, nvars, weights, solvals, varsc2, varscn2, &nvarsc2, &nvarscn2, &varsc2weight);
   
   /* get a most violated minimal cover C = C2 & C1, considering that the variables for downlifting (varsc2) are fixed 
    * to one 
    */
   SCIP_CALL( getCover(scip, vars, nvars, weights, capacity, solvals, varsc2, varscn2, nvarsc2, nvarscn2, varsc2weight,
         covervars, noncovervars, &ncovervars, &ncovervarsc1, &ncovervarsc2, &nnoncovervars, 
         &coverweight, &covervarsc1activity, &coverfound) );
   assert(ncovervars + nnoncovervars == nvars);

   if( coverfound ) 
   {
      SCIP_Longint* initminweight;

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
      SCIP_CALL( SCIPallocBufferArray(scip, &initminweight, ncovervarsc1+1) );

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

      /* sort covervars and noncovervars by non-increasing solution value */
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
         SCIP_Longint oldweight;
         SCIP_Longint newweight;
         int n;

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
            
         /* lift cardinality inequality sum(j in C1) x_j <= |C1| to a valid inequality of the full dimensional knapsack 
          * polytop: 
          *  1. uplifting of noncovervars with solution value > 0 
          *  2. downlifting of covervarsc2 (C2)
          *  3. uplifting of noncovervars with solution value = 0 
          */
         SCIP_CALL( liftKnapsackCover(scip, vars, nvars, weights, capacity, solvals, covervars, noncovervars, 
               ncovervars, ncovervarsc1, ncovervarsc2, nnoncovervars, initminweight,
               liftcoefs, &liftrhs, &liftlpval) );

         /* check, if lifting yielded a violated cut */
         if( SCIPisEfficacious(scip, (covervarsc1activity + liftlpval - liftrhs)/sqrt((SCIP_Real)MAX(liftrhs, 1))) )
         {
            SCIP_ROW* row;
            char name[SCIP_MAXSTRLEN];
            int v;
            
            /* create LP row */
            if( cons != NULL )
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_card%"SCIP_LONGINT_FORMAT"_%d", SCIPconsGetName(cons),
                  SCIPconshdlrGetNCutsFound(SCIPconsGetHdlr(cons)), j);
            else
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "card_%d", j);

            SCIP_CALL( SCIPcreateEmptyRow(scip, &row, name, -SCIPinfinity(scip), (SCIP_Real)liftrhs, 
                  cons != NULL ? SCIPconsIsLocal(cons) : FALSE, FALSE,
                  cons != NULL ? SCIPconsIsRemovable(cons) : TRUE) );
            
            /* add all variables in the knapsack constraint with calculated lifting coefficient to the cut */
            SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
            for( v = 0; v < ncovervars; v++ )
            {
               SCIP_CALL( SCIPaddVarToRow(scip, row, vars[covervars[v]], (SCIP_Real)liftcoefs[covervars[v]]) ); 
            }
            for( v = 0; v < nnoncovervars; v++ )
            {
               if( liftcoefs[noncovervars[v]] > 0 )
               {
                  SCIP_CALL( SCIPaddVarToRow(scip, row, vars[noncovervars[v]], (SCIP_Real)liftcoefs[noncovervars[v]]) );
               }
            }
            SCIP_CALL( SCIPflushRowExtensions(scip, row) );
            
            /* check, if cut is violated enough */
            if( SCIPisCutEfficacious(scip, sol, row) )
            {         
               if( cons != NULL )
               {
                  SCIP_CALL( SCIPresetConsAge(scip, cons) );
               }
               SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE) );
               if( !SCIProwIsLocal(row) )
               {
                  SCIP_CALL( SCIPaddPoolCut(scip, row) );
               }
               (*ncuts)++;
            }
            SCIP_CALL( SCIPreleaseRow(scip, &row) );
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

/* relaxes given general linear constraint into a knapsack constraint and separates lifted knapsack cover inequalities */
SCIP_RETCODE SCIPseparateRelaxedKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that originates the linear constraint, or NULL */
   int                   nknapvars,          /**< number of variables in the continuous knapsack constraint */
   SCIP_VAR**            knapvars,           /**< variables in the continuous knapsack constraint */
   SCIP_Real*            knapvals,           /**< coefficientce of the variables in the continuous knapsack constraint */
   SCIP_Real             valscale,           /**< -1.0 if lhs of row is used as rhs of c. k. constraint, +1.0 otherwise */
   SCIP_Real             rhs,                /**< right hand side of the continuous knapsack constraint */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   int*                  ncuts               /**< pointer to add up the number of found cuts */
   )
{
   SCIP_VAR** binvars;
   SCIP_VAR** consvars;
   SCIP_Real* binvals;
   SCIP_Longint* consvals;
   SCIP_Longint maxact;
   SCIP_Real intscalar;
   SCIP_Bool success;
   int nbinvars;
   int nconsvars;
   int i;

   int* tmpindices;
   int tmp;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool noknapsackconshdlr;

   assert(nknapvars > 0);
   assert(knapvars != NULL);

   tmpindices = NULL;

   SCIPdebugMessage("separate linear constraint <%s> relaxed to knapsack\n", cons != NULL ? SCIPconsGetName(cons) : "-");
   if( cons != NULL )
   {
      SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
   }

   SCIP_CALL( SCIPgetVarsData(scip, &binvars, NULL, &nbinvars, NULL, NULL, NULL) );

   if( nbinvars == 0 )
      return SCIP_OKAY;

   /* set up data structures */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nbinvars) );

   /* get conshdlrdata to use cleared memory */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      noknapsackconshdlr = TRUE;
      SCIP_CALL( SCIPallocBufferArray(scip, &binvals, nbinvars) );
      BMSclearMemoryArray(binvals, nbinvars);
   }
   else
   {
      noknapsackconshdlr = FALSE;
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      SCIP_CALL( SCIPallocBufferArray(scip, &tmpindices, nknapvars) );

      assert(conshdlrdata->reals1size > 0);

      /* next if condition should normally not be true, because it means that presolving has created more binary 
       * variables than binary + integer variables existed at the constraint initialization method, but for example if you would
       * transform all integers into their binary representation then it maybe happens
       */ 
      if( conshdlrdata->reals1size < nbinvars )
      {
	 int oldsize;
	 oldsize = conshdlrdata->reals1size;
	 
	 while( conshdlrdata->reals1size < nbinvars )
            conshdlrdata->reals1size *= 2;
         SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->reals1, conshdlrdata->reals1size) );
         BMSclearMemoryArray(&conshdlrdata->reals1[oldsize], conshdlrdata->reals1size - oldsize);
	}
      binvals = conshdlrdata->reals1;

      /* check for cleared array, all entries have to be zero */
#ifndef NDEBUG
      for( tmp = nbinvars - 1; tmp >= 0; --tmp )
      {
         assert(binvals[tmp] == 0);
      }
#endif
   }

   tmp = 0;

   /* relax continuous knapsack constraint:
    * 1. make all variables binary:
    *    if x_j is continuous or integer variable substitute:
    *      - a_j < 0: x_j = lb  or  x_j = b*z + d with variable lower bound b*z + d with binary variable z
    *      - a_j > 0: x_j = ub  or  x_j = b*z + d with variable upper bound b*z + d with binary variable z
    * 2. convert coefficients of all variables to positive integers:
    *      - scale all coefficients a_j to a~_j integral
    *      - substitute  x~_j = 1 - x_j if a~_j < 0
    */

   /* replace integer and continuous variables with binary variables */
   for( i = 0; i < nknapvars; i++ )
   {
      SCIP_VAR* var;

      var = knapvars[i];

      if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
      {
         assert(0 <= SCIPvarGetProbindex(var) && SCIPvarGetProbindex(var) < nbinvars);
         binvals[SCIPvarGetProbindex(var)] += valscale * knapvals[i];
         if( !noknapsackconshdlr )
         {
            assert(tmpindices != NULL);

            tmpindices[tmp] = SCIPvarGetProbindex(var);
            ++tmp;
         }
         SCIPdebugMessage(" -> binary variable %+.15g<%s>(%.15g)\n", 
            valscale * knapvals[i], SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var));
      }
      else if( valscale * knapvals[i] > 0.0 )
      {
         SCIP_VAR** zvlb;
         SCIP_Real* bvlb;
         SCIP_Real* dvlb;
         SCIP_Real bestlbsol;
         int bestlbtype;
         int nvlb;
         int j;

         /* a_j > 0: substitution with lb or vlb */
         nvlb = SCIPvarGetNVlbs(var);
         zvlb = SCIPvarGetVlbVars(var);
         bvlb = SCIPvarGetVlbCoefs(var);
         dvlb = SCIPvarGetVlbConstants(var);

         /* search for lb or vlb with maximal bound value */
         bestlbsol = SCIPvarGetLbGlobal(var);
         bestlbtype = -1;
         for( j = 0; j < nvlb; j++ )
         {
            /* use only numerical stable vlb with binary variable z */
            if( SCIPvarGetType(zvlb[j]) == SCIP_VARTYPE_BINARY && SCIPvarIsActive(zvlb[j]) 
               && REALABS(bvlb[j]) <= MAXABSVBCOEF )
            {
               SCIP_Real vlbsol;

               assert(0 <= SCIPvarGetProbindex(zvlb[j]) && SCIPvarGetProbindex(zvlb[j]) < nbinvars);
               vlbsol = bvlb[j] * SCIPgetSolVal(scip, sol, zvlb[j]) + dvlb[j];
               if( SCIPisGE(scip, vlbsol, bestlbsol) )
               {
                  bestlbsol = vlbsol;
                  bestlbtype = j;
               }
            }
         }

         /* if no lb or vlb with binary variable was found, we have to abort */
         if( SCIPisInfinity(scip, -bestlbsol) )
            goto TERMINATE;

         if( bestlbtype == -1 )
         {
            rhs -= valscale * knapvals[i] * bestlbsol;
            SCIPdebugMessage(" -> non-binary variable %+.15g<%s>(%.15g) replaced with lower bound %.15g (rhs=%.15g)\n",
               valscale * knapvals[i], SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var), SCIPvarGetLbGlobal(var), rhs);
         }
         else
         {
            assert(0 <= SCIPvarGetProbindex(zvlb[bestlbtype]) && SCIPvarGetProbindex(zvlb[bestlbtype]) < nbinvars);
            rhs -= valscale * knapvals[i] * dvlb[bestlbtype];
            binvals[SCIPvarGetProbindex(zvlb[bestlbtype])] += valscale * knapvals[i] * bvlb[bestlbtype];

            if( SCIPisInfinity(scip, REALABS(binvals[SCIPvarGetProbindex(zvlb[bestlbtype])])) )
               goto TERMINATE;

            if( !noknapsackconshdlr )
            {
               assert(tmpindices != NULL);

               tmpindices[tmp] = SCIPvarGetProbindex(zvlb[bestlbtype]);
               ++tmp;
            }
            SCIPdebugMessage(" -> non-binary variable %+.15g<%s>(%.15g) replaced with variable lower bound %+.15g<%s>(%.15g) %+.15g (rhs=%.15g)\n",
               valscale * knapvals[i], SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var),
               bvlb[bestlbtype], SCIPvarGetName(zvlb[bestlbtype]),
               SCIPgetSolVal(scip, sol, zvlb[bestlbtype]), dvlb[bestlbtype], rhs);
         }
      }
      else
      {
         SCIP_VAR** zvub;
         SCIP_Real* bvub;
         SCIP_Real* dvub;
         SCIP_Real bestubsol;
         int bestubtype;
         int nvub;
         int j;

         assert(valscale * knapvals[i] < 0.0);

         /* a_j < 0: substitution with ub or vub */
         nvub = SCIPvarGetNVubs(var);
         zvub = SCIPvarGetVubVars(var);
         bvub = SCIPvarGetVubCoefs(var);
         dvub = SCIPvarGetVubConstants(var);

         /* search for ub or vub with minimal bound value */
         bestubsol = SCIPvarGetUbGlobal(var);
         bestubtype = -1;
         for( j = 0; j < nvub; j++ )
         {
            /* use only numerical stable vub with active binary variable z */
            if( SCIPvarGetType(zvub[j]) == SCIP_VARTYPE_BINARY && SCIPvarIsActive(zvub[j]) 
               && REALABS(bvub[j]) <= MAXABSVBCOEF )
            {
               SCIP_Real vubsol;

               assert(0 <= SCIPvarGetProbindex(zvub[j]) && SCIPvarGetProbindex(zvub[j]) < nbinvars);
               vubsol = bvub[j] * SCIPgetSolVal(scip, sol, zvub[j]) + dvub[j];
               if( SCIPisLE(scip, vubsol, bestubsol) )
               {
                  bestubsol = vubsol;
                  bestubtype = j;
               }
            }
         }

         /* if no ub or vub with binary variable was found, we have to abort */
         if( SCIPisInfinity(scip, bestubsol) )
            goto TERMINATE;

         if( bestubtype == -1 )
         {
            rhs -= valscale * knapvals[i] * bestubsol;
            SCIPdebugMessage(" -> non-binary variable %+.15g<%s>(%.15g) replaced with upper bound %.15g (rhs=%.15g)\n",
               valscale * knapvals[i], SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var), SCIPvarGetUbGlobal(var), rhs);
         }
         else
         {
            assert(0 <= SCIPvarGetProbindex(zvub[bestubtype]) && SCIPvarGetProbindex(zvub[bestubtype]) < nbinvars);
            rhs -= valscale * knapvals[i] * dvub[bestubtype];
            binvals[SCIPvarGetProbindex(zvub[bestubtype])] += valscale * knapvals[i] * bvub[bestubtype];

            if( SCIPisInfinity(scip, REALABS(binvals[SCIPvarGetProbindex(zvub[bestubtype])])) )
               goto TERMINATE;

            if( !noknapsackconshdlr )
            {
               assert(tmpindices != NULL);

               tmpindices[tmp] = SCIPvarGetProbindex(zvub[bestubtype]);
               ++tmp;
            }
            SCIPdebugMessage(" -> non-binary variable %+.15g<%s>(%.15g) replaced with variable upper bound %+.15g<%s>(%.15g) %+.15g (rhs=%.15g)\n",
               valscale * knapvals[i], SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var),
               bvub[bestubtype], SCIPvarGetName(zvub[bestubtype]),
               SCIPgetSolVal(scip, sol, zvub[bestubtype]), dvub[bestubtype], rhs);
         }
      }
   }

   /* convert coefficents of all (now binary) variables to positive integers:
    *   - make all coefficients integral
    *   - make all coefficients positive (substitute negated variable)
    */
   nconsvars = 0;

   /* calculate scalar which makes all coefficients integral */
   SCIP_CALL( SCIPcalcIntegralScalar(binvals, nbinvars, -SCIPepsilon(scip), KNAPSACKRELAX_MAXDELTA,
         KNAPSACKRELAX_MAXDNOM, KNAPSACKRELAX_MAXSCALE, &intscalar, &success) );
   SCIPdebugMessage(" -> intscalar = %.15g\n", intscalar);

   /* if coefficients can not be made integral, we have to use a scalar of 1.0 and only round fractional coefficients down */
   if( !success )
      intscalar = 1.0;

   /* make all coefficients integral and positive:
    *  - scale a~_j = a_j * intscalar
    *  - substitute x~_j = 1 - x_j if a~_j < 0
    */
   rhs = rhs*intscalar;

   SCIPdebugMessage(" -> rhs = %.15g\n", rhs);
   maxact = 0;
   for( i = 0; i < nbinvars; i++ )
   {
      SCIP_VAR* var;
      SCIP_Longint val;

      val = (SCIP_Longint)SCIPfloor(scip, binvals[i]*intscalar);
      if( val == 0 )
         continue;

      if( val > 0 )
      {
         var = binvars[i];
         SCIPdebugMessage(" -> positive scaled binary variable %+"SCIP_LONGINT_FORMAT"<%s> (unscaled %.15g): not changed (rhs=%.15g)\n",
            val, SCIPvarGetName(var), binvals[i], rhs);
      }
      else
      {
         assert(val < 0);

         SCIP_CALL( SCIPgetNegatedVar(scip, binvars[i], &var) );
         val = -val;
         rhs += val;
         SCIPdebugMessage(" -> negative scaled binary variable %+"SCIP_LONGINT_FORMAT"<%s> (unscaled %.15g): substituted by (1 - <%s>) (rhs=%.15g)\n",
            -val, SCIPvarGetName(binvars[i]), binvals[i], SCIPvarGetName(var), rhs);
      }

      maxact += val;
      consvals[nconsvars] = val;
      consvars[nconsvars] = var;
      nconsvars++;
   }

   if( nconsvars > 0 )
   {
      SCIP_Longint capacity;

      assert(consvars != NULL);
      assert(consvals != NULL);
      capacity = (SCIP_Longint)SCIPfeasFloor(scip, rhs);
      if( maxact > capacity )
      {
#ifdef SCIP_DEBUG
         SCIP_Real act;

         SCIPdebugMessage(" -> linear constraint <%s> relaxed to knapsack:", cons != NULL ? SCIPconsGetName(cons) : "-");
         act = 0.0;
         for( i = 0; i < nconsvars; ++i )
         {
            SCIPdebugPrintf(" %+"SCIP_LONGINT_FORMAT"<%s>(%.15g)", consvals[i], SCIPvarGetName(consvars[i]),
               SCIPgetSolVal(scip, sol, consvars[i]));
            act += consvals[i] * SCIPgetSolVal(scip, sol, consvars[i]);
         }
         SCIPdebugPrintf(" <= %"SCIP_LONGINT_FORMAT" (%.15g) [act: %.15g, max: %"SCIP_LONGINT_FORMAT"]\n",
            capacity, rhs, act, maxact);
#endif

         /* separate lifted cut from relaxed knapsack constraint */
         SCIP_CALL( SCIPseparateKnapsackCover(scip, cons, consvars, nconsvars, consvals, capacity, sol, -1, ncuts) );
      }
   }

 TERMINATE:
   /* free data structures */
   if( noknapsackconshdlr)
   {
      SCIPfreeBufferArray(scip, &binvals);
   }
   else
   {
      /* clear binvals */
      for( --tmp; tmp >= 0; --tmp)
      {
         assert(tmpindices != NULL);
         binvals[tmpindices[tmp]] = 0;
      }
      SCIPfreeBufferArray(scip, &tmpindices);
   }
   SCIPfreeBufferArray(scip, &consvals);
   SCIPfreeBufferArray(scip, &consvars);

   return SCIP_OKAY;
}

/** separates given knapsack constraint */
static
SCIP_RETCODE separateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool             sepacardinality,    /**< should knapsack cardinality cuts be separated? */
   int                   maxnumcardlift,     /**< maximal number of cardinality inequ. lifted per sepa round (-1: unlimited) */
   int*                  ncuts               /**< pointer to add up the number of found cuts */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool violated;

   assert(ncuts != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("separating knapsack constraint <%s>\n", SCIPconsGetName(cons));
   
   /* check knapsack constraint itself for feasibility */
   SCIP_CALL( checkCons(scip, cons, sol, (sol != NULL), FALSE, &violated) );
   
   if( violated )
   {
      /* add knapsack constraint as LP row to the LP */
      SCIP_CALL( addRelaxation(scip, cons, sol) );
      (*ncuts)++;
   }
   else if( sepacardinality )
   {
      SCIP_CALL( SCIPseparateKnapsackCover(scip, cons, consdata->vars, consdata->nvars, consdata->weights, 
            consdata->capacity, sol, maxnumcardlift, ncuts) );
   }
   
   return SCIP_OKAY;
}

/** propagation method for knapsack constraints */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   SCIP_Bool*            redundant,          /**< pointer to store whether constraint is redundant */
   int*                  nfixedvars          /**< pointer to count number of fixings */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;
   SCIP_Longint zerosweightsum;
   SCIP_Longint onesweightsum;
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

   SCIPdebugMessage("propagating knapsack constraint <%s>\n", SCIPconsGetName(cons));

   /* increase age of constraint; age is reset to zero, if a conflict or a propagation was found */
   if( !SCIPinRepropagation(scip) )
   {
      SCIP_CALL( SCIPincConsAge(scip, cons) );
   }

   do
   {
      /* store the sum of weights of fixed-to-one variables locally, because the constraint data's sum can
       * change due to event processing
       */
      onesweightsum = consdata->onesweightsum;

      /* check, if weights of fixed variables already exceeds knapsack capacity */
      if( consdata->capacity < onesweightsum )
      {
         SCIPdebugMessage(" -> cutoff - fixed weight: %"SCIP_LONGINT_FORMAT", capacity: %"SCIP_LONGINT_FORMAT"\n", 
            onesweightsum, consdata->capacity);
            
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         *cutoff = TRUE;

         if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
         {
            /* start conflict analysis with the fixed-to-one variables */
            SCIP_CALL( SCIPinitConflictAnalysis(scip) );
            for( i = 0; i < consdata->nvars; i++ )
            {
               if( SCIPvarGetLbLocal(consdata->vars[i]) > 0.5)
               {
                  SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[i]) );
               }
            }
         
            SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
         }

         return SCIP_OKAY;
      }

      /* make sure, the items are sorted by non-increasing weight */
      sortItems(consdata);

      /* fix all variables to zero, that don't fit into the knapsack anymore */
      zerosweightsum = 0;
      for( i = 0; i < consdata->nvars; ++i )
      {
         if( consdata->weights[i] > consdata->capacity - consdata->onesweightsum )
         {
            if( SCIPvarGetLbLocal(consdata->vars[i]) < 0.5 )
            {
               if( SCIPvarGetUbLocal(consdata->vars[i]) > 0.5 )
               {
                  SCIPdebugMessage(" -> fixing variable <%s> to 0\n", SCIPvarGetName(consdata->vars[i]));
                  SCIP_CALL( SCIPresetConsAge(scip, cons) );
                  SCIP_CALL( SCIPinferBinvarCons(scip, consdata->vars[i], FALSE, cons, i, &infeasible, &tightened) );
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
   if( !SCIPconsIsModifiable(cons) && consdata->weightsum - zerosweightsum <= consdata->capacity )
   {
      SCIPdebugMessage(" -> knapsack constraint <%s> is redundant: weightsum=%"SCIP_LONGINT_FORMAT", zerosweightsum=%"SCIP_LONGINT_FORMAT", capacity=%"SCIP_LONGINT_FORMAT"\n",
         SCIPconsGetName(cons), consdata->weightsum, zerosweightsum, consdata->capacity);
      SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      *redundant = TRUE;
   }

   /* mark the constraint propagated */
   consdata->propagated = TRUE;

   return SCIP_OKAY;
}

/** adds coefficient to constraint data */
static
SCIP_RETCODE addCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_VAR*             var,                /**< variable to add to knapsack */
   SCIP_Longint          weight              /**< weight of variable in knapsack */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
   assert(weight > 0);

   /* add the new coefficient to the LP row */
   if( consdata->row != NULL )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, var, (SCIP_Real)weight) );
   }

   /* check for fixed variable */
   if( SCIPvarGetLbGlobal(var) > 0.5 )
   {
      /* variable is fixed to one: reduce capacity */
      consdata->capacity -= weight;
   }
   else if( SCIPvarGetUbGlobal(var) > 0.5 )
   {
      SCIP_Bool negated;

      /* get binary representative of variable */
      SCIP_CALL( SCIPgetBinvarRepresentative(scip, var, &var, &negated) );

      /* insert coefficient */
      SCIP_CALL( consdataEnsureVarsSize(scip, consdata, consdata->nvars+1, SCIPconsIsTransformed(cons)) );
      consdata->vars[consdata->nvars] = var;
      consdata->weights[consdata->nvars] = weight;
      consdata->nvars++;

      /* install the rounding locks of variable */
      SCIP_CALL( lockRounding(scip, cons, var) );

      /* catch events */
      if( SCIPconsIsTransformed(cons) )
      {
         SCIP_CONSHDLRDATA* conshdlrdata;

         conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
         assert(conshdlrdata != NULL);
         SCIP_CALL( eventdataCreate(scip, &consdata->eventdatas[consdata->nvars-1], consdata, weight) );
         SCIP_CALL( SCIPcatchVarEvent(scip, var,
               SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_UBRELAXED | SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_IMPLADDED,
               conshdlrdata->eventhdlr, consdata->eventdatas[consdata->nvars-1],
               &consdata->eventdatas[consdata->nvars-1]->filterpos) );
      }

      /* update weight sums */
      consdata->weightsum += weight;
      if( SCIPvarGetLbLocal(var) > 0.5 )
         consdata->onesweightsum += weight;

      consdata->sorted = FALSE;
      consdata->cliquepartitioned = FALSE;
      consdata->merged = FALSE;
   }
   consdata->propagated = FALSE;
   consdata->presolved = FALSE;
   consdata->cliquesadded = FALSE; /* new coefficient might lead to larger cliques */

   return SCIP_OKAY;
}

/** deletes coefficient at given position from constraint data */
static
SCIP_RETCODE delCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int                   pos                 /**< position of coefficient to delete */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(consdata->vars[pos]));

   /* delete the coefficient from the LP row */
   if( consdata->row != NULL )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, consdata->vars[pos], -(SCIP_Real)consdata->weights[pos]) );
   }

   /* remove the rounding locks of variable */
   SCIP_CALL( unlockRounding(scip, cons, consdata->vars[pos]) );

   /* drop events */
   if( SCIPconsIsTransformed(cons) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      
      conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
      assert(conshdlrdata != NULL);
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[pos],
            SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_UBRELAXED | SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_IMPLADDED,
            conshdlrdata->eventhdlr, consdata->eventdatas[pos], consdata->eventdatas[pos]->filterpos) );
      SCIP_CALL( eventdataFree(scip, &consdata->eventdatas[pos]) );
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
   if( consdata->eventdatas != NULL )
      consdata->eventdatas[pos] = consdata->eventdatas[consdata->nvars-1];
#if 0 /* not neccessary, because clique partition will be made invalid */
   if( consdata->cliquepartition != NULL )
      consdata->cliquepartition[pos] = consdata->cliquepartition[consdata->nvars-1];
#endif
   consdata->nvars--;

   consdata->propagated = FALSE;
   consdata->presolved = FALSE;
   consdata->sorted = FALSE;
   consdata->cliquepartitioned = FALSE; /* recalculate the clique partition after a coefficient was removed */

   return SCIP_OKAY;
}

/** removes all items with weight zero from knapsack constraint */
static
SCIP_RETCODE removeZeroWeights(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< knapsack constraint */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   for( v = consdata->nvars-1; v >= 0; --v )
   {
      if( consdata->weights[v] == 0 )
      {
         SCIP_CALL( delCoefPos(scip, cons, v) );
      }
   }

   return SCIP_OKAY;
}

/** replaces multiple occurrences of a variable or its negation by a single coefficient */
static
SCIP_RETCODE mergeMultiples(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< knapsack constraint */
   )
{
   SCIP_CONSDATA* consdata;
   int v;
   int prev;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->merged ) 
      return SCIP_OKAY;

   if( consdata->nvars <= 1 )
   {
      consdata->merged = TRUE;
      return SCIP_OKAY;
   }

   assert(consdata->vars != NULL || consdata->nvars == 0);

   /* sorting array after indices of variables, that's only for faster merging */ 
   SCIPsortPtrPtrLongInt((void**)consdata->vars, (void**)consdata->eventdatas, consdata->weights, consdata->cliquepartition,
			 SCIPvarCompActiveAndNegated, consdata->nvars);
   /* knapsack-sorting (descreasing weights) now lost */ 
   consdata->sorted = FALSE;

   v = consdata->nvars - 1;
   /* loop backwards through the items: deletion only affects rear items */
   for( prev = v - 1; prev >= 0; --prev )
   {
      SCIP_VAR* var1;
      SCIP_VAR* var2;
      SCIP_Bool negated1;
      SCIP_Bool negated2;
      
      negated1 = FALSE;
      negated2 = FALSE;
      
      var1 = consdata->vars[v];
      assert(SCIPvarGetType(var1) == SCIP_VARTYPE_BINARY);
      assert(SCIPvarIsActive(var1) || SCIPvarGetStatus(var1) == SCIP_VARSTATUS_NEGATED);
      if( SCIPvarGetStatus(var1) == SCIP_VARSTATUS_NEGATED )
      {
         var1 = SCIPvarGetNegatedVar(var1);
         negated1 = TRUE;
      }
      assert(var1 != NULL);
      
      var2 = consdata->vars[prev];
      assert(SCIPvarGetType(var2) == SCIP_VARTYPE_BINARY);
      assert(SCIPvarIsActive(var2) || SCIPvarGetStatus(var2) == SCIP_VARSTATUS_NEGATED);
      if( SCIPvarGetStatus(var2) == SCIP_VARSTATUS_NEGATED )
      {
         var2 = SCIPvarGetNegatedVar(var2);
         negated2 = TRUE;
      }
      assert(var2 != NULL);
      
      if( var1 == var2 )
      {
         /* both variables are either active or negated */
         if( negated1 == negated2 )
         {
            /* variables var1 and var2 are equal: add weight of var1 to var2, and delete var1 */
            consdataChgWeight(consdata, prev, consdata->weights[v] + consdata->weights[prev]);
            SCIP_CALL( delCoefPos(scip, cons, v) );
         }
         /* variables var1 and var2 are opposite: subtract smaller weight from larger weight, reduce capacity,
          * and delete item of smaller weight
          */
         else if( consdata->weights[v] == consdata->weights[prev] )
         {
            /* both variables eliminate themselves: w*x + w*(1-x) == w */
            consdata->capacity -= consdata->weights[v];
            SCIP_CALL( delCoefPos(scip, cons, v) ); /* this does not affect var2, because var2 stands before var1 */
            SCIP_CALL( delCoefPos(scip, cons, prev) );

            --prev;
         }
         else if( consdata->weights[v] < consdata->weights[prev] )
         {
            consdata->capacity -= consdata->weights[v];
            consdataChgWeight(consdata, prev, consdata->weights[prev] - consdata->weights[v]);
            assert(consdata->weights[prev] > 0);
            SCIP_CALL( delCoefPos(scip, cons, v) ); /* this does not affect var2, because var2 stands before var1 */
         }
         else
         {
            consdata->capacity -= consdata->weights[prev];
            consdataChgWeight(consdata, v, consdata->weights[v] - consdata->weights[prev]);
            assert(consdata->weights[v] > 0);
            SCIP_CALL( delCoefPos(scip, cons, prev) ); /* attention: normally we lose our order */
            /* restore order iff necessary */
            if( consdata->nvars != v ) /* otherwise the order still stands */
            {
               /* either that was the last pair or both, the negated and "normal" variable in front doesn't match var1, so the order is irrelevant */
               if( prev == 0 || (var1 != consdata->vars[prev - 1] && var1 != SCIPvarGetNegatedVar(consdata->vars[prev - 1])) )
                  --prev;
               else /* we need to let v at the same position*/
               {
                  consdata->cliquesadded = FALSE; /* new coefficient might lead to larger cliques */
                  /* don't decrease v, the same variable exist in front */
                  continue;
               }
            }
         }
         consdata->cliquesadded = FALSE; /* new coefficient might lead to larger cliques */
      }
      v = prev;
     }

   consdata->merged = TRUE;

   return SCIP_OKAY;
}

/*  tries to simplify coefficients and delete variables in knapsack a^Tx <= rhs
 *  in case there is only one binary variable with an odd weight and the capacity
 *  is odd too, then:
 *    - if the odd coefficient is equal to 1, delete the variable and decrease rhs by 1
 *    - otherwise, decrease coefficient and rhs by 1
 *  Afterwards we use the normalize method to further simplify the inequality. 
 */
static
SCIP_RETCODE simplifyInequalities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int*                  ndelconss,          /**< pointer to store the amount of deleted constraints */
   int*                  nchgcoefs,          /**< pointer to store the amount of changed coefficients */
   int*                  nchgsides,          /**< pointer to store the amount of changed sides */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool success;
   int v;
   int nvars;
   SCIP_VAR** vars;
   SCIP_Longint* weights;
   SCIP_Longint oddbinval;
   SCIP_Longint capacity;
   SCIP_Longint gcd;
   
   SCIP_VAR* oddbinvar;
   int noddvals = 0;
   int pos = 0;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( ndelconss != NULL );
   assert( nchgcoefs != NULL );
   assert( nchgsides != NULL );
   assert( cutoff != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   *cutoff = FALSE;

   if( !consdata->merged )
      return SCIP_OKAY;

   /* check if capacity is odd */
   if( !(consdata->capacity % 2) )
      return SCIP_OKAY;


   /* try to delete variables and simplify constraint */
   do
   {
      success = FALSE;
      noddvals = 0;

      capacity = consdata->capacity;
      vars = consdata->vars;
      weights = consdata->weights;
      nvars = consdata->nvars;
      
      assert(nvars > 0);

      oddbinvar = NULL;
      oddbinval = 0;
      
      /* search for binary variables with an odd coefficient */
      for( v = 0; v < nvars; ++v )
      {
         /* check if the coefficient is odd */
	if( weights[v] % 2 )
         {
            oddbinvar = vars[v];
            oddbinval = weights[v];
            pos = v;
            noddvals++;
            if (noddvals >= 2)
               return SCIP_OKAY;
         }
      }
      
      /* now we found exactly one binary variables with an odd coefficient and all other variables have even
       * coefficients and are not of continuous type; furthermore, only one side is odd */
      if( noddvals )
      {
         assert(noddvals == 1);
         assert(oddbinvar != NULL);
         
         oddbinval--;
         SCIPdebugMessage("knapsack constraint <%s>: decreasing coefficient for variable <%s> to <%lld> and rhs to <%lld>\n", 
            SCIPconsGetName(cons), SCIPvarGetName(oddbinvar), oddbinval , capacity - 1);

	 if( consdata->capacity == 0 )
	 {
	    *cutoff = TRUE;
	    return SCIP_OKAY;
	 }

         --(consdata->capacity);

         if( oddbinval == 0 )
         {
            SCIP_CALL( delCoefPos( scip, cons, pos ) );
	    /* if the last variable was erased delete the constraint too */
	    if( nvars == 1 )
	    {
	       assert(consdata->capacity >= 0);
	       SCIP_CALL( SCIPdelConsLocal(scip, cons) );
	       ++(*ndelconss);
	       return SCIP_OKAY;
	    }
         }
         else
         {
            consdataChgWeight(consdata, pos, oddbinval);
         }

         ++(*nchgcoefs);
         ++(*nchgsides);

         capacity = consdata->capacity;
         vars = consdata->vars;
         weights = consdata->weights;
         nvars = consdata->nvars;

         /*
          * division by greatest common divisor
          */
         gcd = weights[nvars - 1];
         for( v = nvars - 2; v >= 0 && gcd > 1; --v )
         {
            gcd = SCIPcalcGreComDiv(gcd, ABS(weights[v]));
         }
         
         assert(gcd >= 2);
         
         for( v = nvars - 1; v >= 0 && gcd > 1; --v )
         {
            consdataChgWeight(consdata, v, weights[v]/gcd);
         }
         (*nchgcoefs) += nvars;
         
         consdata->capacity /= gcd;
         (*nchgsides)++;
         
	 /* capacity still odd, try to simplify again */
	 if( consdata->capacity % 2 )
	   success = TRUE;
      }
   }
   while( success );
   
   return SCIP_OKAY;
}

/** deletes all fixed variables from knapsack constraint, and replaces variables with binary representatives */
static
SCIP_RETCODE applyFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);
   
   *cutoff = FALSE;

   SCIPdebugMessage("apply fixings:\n");
   SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons, NULL) ));
   
   v = 0;
   while( v < consdata->nvars )
   {
      SCIP_VAR* var;

      var = consdata->vars[v];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);

      if( SCIPvarGetLbGlobal(var) > 0.5 )
      {
         assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(var), 1.0));
         consdata->capacity -= consdata->weights[v];
         SCIP_CALL( delCoefPos(scip, cons, v) );
         consdata->cliquesadded = FALSE; /* reduced capacity might lead to larger cliques */
      }
      else if( SCIPvarGetUbGlobal(var) < 0.5 )
      {
         assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 0.0));
         SCIP_CALL( delCoefPos(scip, cons, v) );
      }
      else
      {
         SCIP_VAR* repvar;
	 SCIP_VAR* negvar;
         SCIP_VAR* workvar;
         SCIP_Longint weight;
         SCIP_Bool negated;
	 
	 weight = consdata->weights[v];
	 
         /* get binary representative of variable */
         SCIP_CALL( SCIPgetBinvarRepresentative(scip, var, &repvar, &negated) );
	 assert(repvar != NULL);

	 /* check for multiaggregation */
	 if( SCIPvarIsNegated(repvar) )
	 {
	    workvar = SCIPvarGetNegatedVar(repvar);
	    assert(workvar != NULL);
	    negated = TRUE;
	 }
	 else
	 {
	    workvar = repvar;
	    negated = FALSE;
	 }

	 /* @todo maybe resolve the problem that the eliminating of the multi-aggreation leads to a non-knapsack
          * constraint (converting into a linear constraint), for example the multiaggregation consist of a non-binary
          * variable or due to resolving now their are non-integral coefficients or a non-integral capacity 
          *
          * If repvar is not negated so workwar = repvar, otherwise workvar = 1 - repvar. This means,
          * weight * workvar = weight * (a_1*y_1 + ... + a_n*y_n + c) 
          *
          * The explaination for  the following block:  
	  * 1a) If repvar is a multiaggregated variable weight * repvar should be replaced by 
	  *     weight * (a_1*y_1 + ... + a_n*y_n + c).
	  * 1b) If repvar is a negated variable of a multiaggregated variable weight * repvar should be replaced by 
	  *     weight - weight * (a_1*y_1 + ... + a_n*y_n + c), for better further use here we switch the sign of weight
	  *     so now we have the replacement -weight + weight * (a_1*y_1 + ... + a_n*y_n + c).
	  * 2)  For all replacement variable we check:
	  * 2a) weight * a_i < 0 than we add -weight * a_i * y_i_neg to the constraint and adjust the capacity through 
	  *     capacity -= weight * a_i caused by the negation of y_i.
	  * 2b) weight * a_i >= 0 than we add weight * a_i * y_i to the constraint.
	  * 3a) If repvar was not negated we need to substract weight * c from capacity.
	  * 3b) If repvar was negated we need to substract weight * (c - 1) from capacity(note we switched the sign of 
	  *     weight in this case.
	  */
	 if( SCIPvarGetStatus(workvar) == SCIP_VARSTATUS_MULTAGGR )
	 {
	    SCIP_VAR** aggrvars;
	    SCIP_Real* aggrscalars;
	    SCIP_Real aggrconst;
	    int naggrvars;
	    int i;

	    SCIP_CALL( SCIPflattenVarAggregationGraph(scip, workvar) );
	    naggrvars = SCIPvarGetMultaggrNVars(workvar);
	    aggrvars = SCIPvarGetMultaggrVars(workvar);
	    aggrscalars = SCIPvarGetMultaggrScalars(workvar);
	    aggrconst = SCIPvarGetMultaggrConstant(workvar);
	    assert((aggrvars != NULL && aggrscalars != NULL) || naggrvars == 0);

	    if( !SCIPisIntegral(scip, weight * aggrconst) )
            {
               SCIPerrorMessage("try to resolve a multiaggregation with a non-integral value for weight*aggrconst = %g\n", weight*aggrconst);
               return SCIP_ERROR;
	    }
            
	    /* if workvar was negated, we have to flip the weight */
	    if( negated )
	      weight *= -1;
            
	    for( i = naggrvars - 1; i >= 0; --i )
	    {
               assert(aggrvars != NULL);
               assert(aggrscalars != NULL);

	       if( SCIPvarGetType(aggrvars[i]) != SCIP_VARTYPE_BINARY )
               {
                  SCIPerrorMessage("try to resolve a multiaggregation with a non-binary variable <%s>\n", aggrvars[i]);
                  return SCIP_ERROR;
               }
	       if( !SCIPisIntegral(scip, weight * aggrscalars[i]) )
               {
                  SCIPerrorMessage("try to resolve a multiaggregation with a non-integral value for weight*aggrscalars = %g\n", weight*aggrscalars[i]);
                  return SCIP_ERROR;
               }
	       /* if the new coefficent is smaller than zero, we need to add the negated variable instead and adjust the capacity */
	       if( SCIPisNegative(scip, weight * aggrscalars[i]) )
	       {
		  SCIP_CALL( SCIPgetNegatedVar(scip, aggrvars[i], &negvar));
		  assert(negvar != NULL);
		  SCIP_CALL( addCoef(scip, cons, negvar, (SCIP_Longint)(SCIPfloor(scip, -weight * aggrscalars[i] + 0.5))) );
		  consdata->capacity -= (SCIP_Longint)(SCIPfloor(scip, weight * aggrscalars[i] + 0.5));
	       }
	       else 
	       {
		  SCIP_CALL( addCoef(scip, cons, aggrvars[i], (SCIP_Longint)(SCIPfloor(scip, weight * aggrscalars[i] + 0.5))) );
	       }
	    }
            /* delete old coefficient */
	    SCIP_CALL( delCoefPos(scip, cons, v) );

	    /* adjust the capacity with the aggregation constant and if necessary the extra weight through tne
             * negation */ 
	    if(negated)
	       consdata->capacity -= (SCIP_Longint)SCIPfloor(scip, weight * (aggrconst - 1) + 0.5);
	    else
	       consdata->capacity -= (SCIP_Longint)SCIPfloor(scip, weight * aggrconst + 0.5);

            if( consdata->capacity < 0 )
            {
               *cutoff = TRUE;
               break;
	    }
	 }
	 /* check, if the variable should be replaced with the representative */
         else if( repvar != var )
         {
	    /* delete old (aggregated) variable */
	    SCIP_CALL( delCoefPos(scip, cons, v) );
	    
	    /* add representative instead */
	    SCIP_CALL( addCoef(scip, cons, repvar, weight) );
	 }
         else
            ++v;
      }
   }
   assert(consdata->onesweightsum == 0);

   SCIPdebugMessage("after applyFixings, before merging:\n");
   SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons, NULL) ));
   
   /* if aggregated variables have been replaced, multiple entries of the same variable are possible and we have to
    * clean up the constraint
    */
   if( !(*cutoff) )
   {
      SCIP_CALL( mergeMultiples(scip, cons) );
      SCIPdebugMessage("after applyFixings and merging:\n");
      SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons, NULL) ));
   }
   
   return SCIP_OKAY;
}

/** divides weights by their greatest common divisor and divides capacity by the same value, rounding down the result */
static
void normalizeWeights(
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int*                  nchgcoefs,          /**< pointer to count total number of changed coefficients */
   int*                  nchgsides           /**< pointer to count number of side changes */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Longint gcd;
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

   gcd = (SCIP_Longint)consdata->weights[consdata->nvars-1];
   for( i = consdata->nvars-2; i >= 0 && gcd >= 2; --i )
   {
      assert(SCIPvarGetLbLocal(consdata->vars[i]) < 0.5);
      assert(SCIPvarGetUbLocal(consdata->vars[i]) > 0.5); /* all fixed variables should have been removed */

      gcd = SCIPcalcGreComDiv(gcd, (SCIP_Longint)consdata->weights[i]);
   }

   if( gcd >= 2 )
   {
      SCIPdebugMessage("knapsack constraint <%s>: dividing weights by %"SCIP_LONGINT_FORMAT"\n", SCIPconsGetName(cons), gcd);

      for( i = 0; i < consdata->nvars; ++i )
         consdataChgWeight(consdata, i, consdata->weights[i]/gcd);
      consdata->capacity /= gcd;
      (*nchgcoefs) += consdata->nvars;
      (*nchgsides)++;
   }
}

/** inserts an element into the list of binary zero implications */
static
SCIP_RETCODE insertZerolist(
   SCIP*                 scip,               /**< SCIP data structure */
   int**                 liftcands,          /**< array of the lifting candidates */
   int*                  nliftcands,         /**< number of lifting candidates */
   int**                 firstidxs,          /**< array of first zeroitems indices */
   SCIP_Longint**        zeroweightsums,     /**< array of sums of weights of the implied-to-zero items */
   int**                 zeroitems,          /**< pointer to zero items array */
   int**                 nextidxs,           /**< pointer to array of next zeroitems indeces */
   int*                  zeroitemssize,      /**< pointer to size of zero items array */
   int*                  nzeroitems,         /**< pointer to length of zero items array */
   int                   probindex,          /**< problem index of variable y in implication y == v -> x == 0 */
   SCIP_Bool             value,              /**< value v of variable y in implication */
   int                   knapsackidx,        /**< index of variable x in knapsack */
   SCIP_Longint          knapsackweight,     /**< weight of variable x in knapsack */
   SCIP_Bool*            memlimitreached     /**< pointer to store whether the memory limit was reached */
   )
{
   assert(liftcands != NULL);
   assert(liftcands[value] != NULL);
   assert(nliftcands != NULL);
   assert(firstidxs != NULL);
   assert(firstidxs[value] != NULL);
   assert(zeroweightsums != NULL);
   assert(zeroweightsums[value] != NULL);
   assert(zeroitems != NULL);
   assert(nextidxs != NULL);
   assert(zeroitemssize != NULL);
   assert(nzeroitems != NULL);
   assert(*nzeroitems <= *zeroitemssize);
   assert(0 <= probindex && probindex < SCIPgetNBinVars(scip));
   assert(memlimitreached != NULL);

   /* allocate enough memory */
   if( *nzeroitems == *zeroitemssize )
   {
      /* we explicitly construct the complete implication graph where the knapsack variables are involved;
       * this can be too huge - abort on memory limit
       */
      if( *zeroitemssize >= MAX_ZEROITEMS_SIZE )
      {
         SCIPdebugMessage("memory limit of %d bytes reached in knapsack preprocessing - abort collecting zero items\n",
            *zeroitemssize);
         *memlimitreached = TRUE;
         return SCIP_OKAY;
      }
      *zeroitemssize *= 2;
      *zeroitemssize = MIN(*zeroitemssize, MAX_ZEROITEMS_SIZE);
      SCIP_CALL( SCIPreallocBufferArray(scip, zeroitems, *zeroitemssize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, nextidxs, *zeroitemssize) );
   }
   assert(*nzeroitems < *zeroitemssize);
   *memlimitreached = FALSE;

   /* insert element */
   (*zeroitems)[*nzeroitems] = knapsackidx;
   (*nextidxs)[*nzeroitems] = firstidxs[value][probindex];
   if( firstidxs[value][probindex] == 0 )
   {
      liftcands[value][nliftcands[value]] = probindex;
      nliftcands[value]++;
   }
   firstidxs[value][probindex] = *nzeroitems;
   (*nzeroitems)++;
   zeroweightsums[value][probindex] += knapsackweight;

   return SCIP_OKAY;
}

/** applies rule (3) of the weight tightening procedure, which can lift other variables into the knapsack:
 *  (3) for a clique C let C(xi == v) := C \ {j: xi == v -> xj == 0}),
 *      let cliqueweightsum(xi == v) := sum(W(C(xi == v)))
 *      if cliqueweightsum(xi == v) < capacity:
 *      - fixing variable xi to v would make the knapsack constraint redundant
 *      - the weight of the variable or its negation (depending on v) can be increased as long as it has the same
 *        redundancy effect:
 *          wi'       := capacity - cliqueweightsum(xi == v)
 *      this rule can also be applied to binary variables not in the knapsack!
 */
static
SCIP_RETCODE tightenWeightsLift(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int*                  nchgcoefs           /**< pointer to count total number of changed coefficients */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** binvars;
   int nbinvars;
   int* liftcands[2];          /* binary variables that have at least one entry in zeroitems */
   int* firstidxs[2];          /* first index in zeroitems for each binary variable/value pair, or zero for empty list */
   SCIP_Longint* zeroweightsums[2]; /* sums of weights of the implied-to-zero items */
   int* zeroitems;             /* item number in knapsack that is implied to zero */
   int* nextidxs;              /* next index in zeroitems for the same binary variable, or zero for end of list */
   int zeroitemssize;
   int nzeroitems;
   SCIP_Bool* zeroiteminserted[2];
   SCIP_Bool memlimitreached;
   int nliftcands[2];
   SCIP_Bool* cliqueused;
   SCIP_Bool* itemremoved;
   SCIP_Longint maxcliqueweightsum;
   SCIP_VAR** addvars;
   SCIP_Longint* addweights;
   SCIP_Longint addweightsum;
   int naddvars;
   int val;
   int i;

   int* tmpindices;
   SCIP_Bool* tmpboolindices;
   int* tmpindices2;
   SCIP_Bool* tmpboolindices2;
   int* tmpindices3;
   SCIP_Bool* tmpboolindices3;
   int tmp;
   int tmp2;
   int tmp3;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(nchgcoefs != NULL);
   assert(!SCIPconsIsModifiable(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL); /* we are in presolve, so no LP row exists */
   assert(consdata->weightsum > consdata->capacity); /* otherwise, the constraint is redundant */
   assert(consdata->nvars > 0);
   assert(consdata->merged);

   /* check if the knapsack has too many items for applying this costly method */
   if( consdata->nvars > MAX_USECLIQUES_SIZE )
      return SCIP_OKAY;

   /* sort items, s.t. the heaviest one is in the first position */
   sortItems(consdata);

   binvars = SCIPgetVars(scip);
   nbinvars = SCIPgetNBinVars(scip);
   assert(nbinvars > 0);

   /* get conshdlrdata to use cleared memory */
   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* allocate temporary memory for the list of implied to zero variables */
   zeroitemssize = MIN(nbinvars, MAX_ZEROITEMS_SIZE); /* initial size of zeroitems buffer */
   SCIP_CALL( SCIPallocBufferArray(scip, &liftcands[0], nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &liftcands[1], nbinvars) );

   assert(conshdlrdata->ints1size > 0);
   assert(conshdlrdata->ints2size > 0);
   assert(conshdlrdata->longints1size > 0);
   assert(conshdlrdata->longints2size > 0);

   /* next if conditions should normally not be true, because it means that presolving has created more binary variables
    * than binary + integer variables existed at the presolving initialization method, but for example if you would transform all 
    * integers into their binary representation then it maybe happens
    */ 
   if( conshdlrdata->ints1size < nbinvars )
   {
      int oldsize;
      oldsize = conshdlrdata->ints1size;

      while( conshdlrdata->ints1size < nbinvars )
         conshdlrdata->ints1size *= 2;
      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->ints1, conshdlrdata->ints1size) );
      BMSclearMemoryArray(&conshdlrdata->ints1[oldsize], conshdlrdata->ints1size - oldsize);
   }
   if( conshdlrdata->ints2size < nbinvars )
   {
      int oldsize;
      oldsize = conshdlrdata->ints2size;

      while( conshdlrdata->ints2size < nbinvars )
         conshdlrdata->ints2size *= 2;
      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->ints2, conshdlrdata->ints2size) );
      BMSclearMemoryArray(&conshdlrdata->ints2[oldsize], conshdlrdata->ints2size - oldsize);
   }
   if( conshdlrdata->longints1size < nbinvars )
   {
      int oldsize;
      oldsize = conshdlrdata->longints1size;

      while( conshdlrdata->longints1size < nbinvars )
         conshdlrdata->longints1size *= 2;
      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->longints1, conshdlrdata->longints1size) );
      BMSclearMemoryArray(&conshdlrdata->longints1[oldsize], conshdlrdata->longints1size - oldsize);
   }
   if( conshdlrdata->longints2size < nbinvars )
   {
      int oldsize;
      oldsize = conshdlrdata->longints2size;

      while( conshdlrdata->longints2size < nbinvars )
         conshdlrdata->longints2size *= 2;
      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->longints2, conshdlrdata->longints2size) );
      BMSclearMemoryArray(&conshdlrdata->longints2[oldsize], conshdlrdata->longints2size - oldsize);
   }

   firstidxs[0] = conshdlrdata->ints1;
   firstidxs[1] = conshdlrdata->ints2;
   zeroweightsums[0] = conshdlrdata->longints1;
   zeroweightsums[1] = conshdlrdata->longints2;

   /* check for cleared arrays, all entries are zero */
#ifndef NDEBUG
   for( tmp = nbinvars - 1; tmp >= 0; --tmp )
   {
      assert(firstidxs[0][tmp] == 0);
      assert(firstidxs[1][tmp] == 0);
      assert(zeroweightsums[0][tmp] == 0);
      assert(zeroweightsums[1][tmp] == 0);
   }
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &zeroitems, zeroitemssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nextidxs, zeroitemssize) );

   zeroitems[0] = -1; /* dummy element */
   nextidxs[0] = -1;
   nzeroitems = 1;
   nliftcands[0] = 0;
   nliftcands[1] = 0;

   assert(conshdlrdata->bools1size > 0);
   assert(conshdlrdata->bools2size > 0);

   /* next if conditions should normally not be true, because it means that presolving has created more binary variables
    * than binary + integer variables existed at the presolving initialization method, but for example if you would transform all 
    * integers into their binary representation then it maybe happens
    */ 
   if( conshdlrdata->bools1size < nbinvars )
   {
      int oldsize;
      oldsize = conshdlrdata->bools1size;

      while( conshdlrdata->bools1size < nbinvars )
         conshdlrdata->bools1size *= 2;
      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->bools1, conshdlrdata->bools1size) );
      BMSclearMemoryArray(&conshdlrdata->bools1[oldsize], conshdlrdata->bools1size - oldsize);
   }
   if( conshdlrdata->bools2size < nbinvars )
   {
      int oldsize;
      oldsize = conshdlrdata->bools2size;

      while( conshdlrdata->bools2size < nbinvars )
         conshdlrdata->bools2size *= 2;
      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->bools2, conshdlrdata->bools2size) );
      BMSclearMemoryArray(&conshdlrdata->bools2[oldsize], conshdlrdata->bools2size - oldsize);
   }

   zeroiteminserted[0] = conshdlrdata->bools1;
   zeroiteminserted[1] = conshdlrdata->bools2;

   /* check for cleared arrays, all entries are zero */
#ifndef NDEBUG
   for( tmp = nbinvars - 1; tmp >= 0; --tmp )
   {
      assert(zeroiteminserted[0][tmp] == 0);
      assert(zeroiteminserted[1][tmp] == 0);
   }
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &tmpindices, 2 * nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpboolindices, 2 * nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpindices2, 2 * nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpboolindices2, 2 * nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpindices3, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpboolindices3, consdata->nvars) );
   tmp2 = 0;
   tmp3 = 0;

   memlimitreached = FALSE;
   for( i = 0; i < consdata->nvars && !memlimitreached; ++i )
   {
      SCIP_VAR* var;
      SCIP_Longint weight;
      SCIP_Bool value;
      int varprobindex;
      SCIP_VAR** implvars;
      SCIP_BOUNDTYPE* impltypes;
      int nimpls;
      SCIP_CLIQUE** cliques;
      int ncliques;
      int j;

      tmp = 0;

      /* get corresponding active problem variable */
      var = consdata->vars[i];
      weight = consdata->weights[i];
      value = TRUE;
      SCIP_CALL( SCIPvarGetProbvarBinary(&var, &value) );
      varprobindex = SCIPvarGetProbindex(var);
      assert(0 <= varprobindex && varprobindex < nbinvars);

      /* update the zeroweightsum */
      zeroweightsums[!value][varprobindex] += weight; /*lint !e514*/
      tmpboolindices3[tmp3] = !value;
      tmpindices3[tmp3] = varprobindex;
      ++tmp3;

      /* initialize the arrays of inserted zero items */
      /* get implications of the knapsack item fixed to one: x == 1 -> y == (1-v);
       * the negation of these implications (y == v -> x == 0) are the ones that we are interested in
       */
      nimpls = SCIPvarGetNBinImpls(var, value);
      implvars = SCIPvarGetImplVars(var, value);
      impltypes = SCIPvarGetImplTypes(var, value);
      for( j = 0; j < nimpls && !memlimitreached; ++j )
      {
         int probindex;
         SCIP_Bool implvalue;

         assert(SCIPvarGetType(implvars[j]) == SCIP_VARTYPE_BINARY);
         probindex = SCIPvarGetProbindex(implvars[j]);
         assert(0 <= probindex && probindex < nbinvars);
         implvalue = (impltypes[j] == SCIP_BOUNDTYPE_UPPER); /* the negation of the implication */
         
         /* insert the item into the list of the implied variable/value */
         if( !zeroiteminserted[implvalue][probindex] )
         {
            if( firstidxs[implvalue][probindex] == 0 )
            {
               tmpboolindices2[tmp2] = implvalue;
               tmpindices2[tmp2] = probindex;
               ++tmp2;
            }
            SCIP_CALL( insertZerolist(scip, liftcands, nliftcands, firstidxs, zeroweightsums,
                  &zeroitems, &nextidxs, &zeroitemssize, &nzeroitems, probindex, implvalue, i, weight,
                  &memlimitreached) );
            zeroiteminserted[implvalue][probindex] = TRUE;
	    tmpboolindices[tmp] = implvalue;
	    tmpindices[tmp] = probindex;
            ++tmp;
         }
      }

      /* get the cliques where the knapsack item is member of with value 1 */
      ncliques = SCIPvarGetNCliques(var, value);
      cliques = SCIPvarGetCliques(var, value);
      for( j = 0; j < ncliques && !memlimitreached; ++j )
      {
         SCIP_VAR** cliquevars;
         SCIP_Bool* cliquevalues;
         int ncliquevars;
         int k;
         
         ncliquevars = SCIPcliqueGetNVars(cliques[j]);
         cliquevars = SCIPcliqueGetVars(cliques[j]);
         cliquevalues = SCIPcliqueGetValues(cliques[j]);
         for( k = 0; k < ncliquevars && !memlimitreached; ++k )
         {
            if( cliquevars[k] != var )
            {
               int probindex;
               SCIP_Bool implvalue;
               
               probindex = SCIPvarGetProbindex(cliquevars[k]);
               assert(0 <= probindex && probindex < nbinvars);
               implvalue = cliquevalues[k];
               
               /* insert the item into the list of the clique variable/value */
               if( !zeroiteminserted[implvalue][probindex] )
               {
                  if( firstidxs[implvalue][probindex] == 0 )
                  {
                     tmpboolindices2[tmp2] = implvalue;
                     tmpindices2[tmp2] = probindex;
                     ++tmp2;
                  }

                  SCIP_CALL( insertZerolist(scip, liftcands, nliftcands, firstidxs, zeroweightsums,
                        &zeroitems, &nextidxs, &zeroitemssize, &nzeroitems, probindex, implvalue, i, weight,
                        &memlimitreached) );
                  zeroiteminserted[implvalue][probindex] = TRUE;
		  tmpboolindices[tmp] = implvalue;
		  tmpindices[tmp] = probindex;
                  ++tmp;
               }
            }
         }
      }
      /* clear zeroiteminserted */
      for( --tmp; tmp >= 0; --tmp)
	zeroiteminserted[tmpboolindices[tmp]][tmpindices[tmp]] = FALSE;
   }
   SCIPfreeBufferArray(scip, &tmpboolindices);

   /* calculate the clique partition and the maximal sum of weights using the clique information */
   assert(consdata->sorted);
   SCIP_CALL( calcCliquepartition(scip, consdata) );

   assert(conshdlrdata->bools3size > 0);

   /* next if condition should normally not be true, because it means that presolving has created more binary variables
    * in one constraint than binary + integer variables existed in the whole problem at the presolving initialization method, but
    * for example if you would transform all integers into their binary representation then it maybe happens
    */ 
   if( conshdlrdata->bools3size < consdata->nvars )
   {
      int oldsize;
      oldsize = conshdlrdata->bools3size;

      while( conshdlrdata->bools3size < consdata->nvars )
         conshdlrdata->bools3size *= 2;
      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->bools3, conshdlrdata->bools3size) );
      BMSclearMemoryArray(&conshdlrdata->bools3[oldsize], conshdlrdata->bools3size - oldsize);
   }

   cliqueused = conshdlrdata->bools3;

   /* check for cleared array, all entries are zero */
#ifndef NDEBUG
   for( tmp = consdata->nvars - 1; tmp >= 0; --tmp )
      assert(cliqueused[tmp] == 0);
#endif

   maxcliqueweightsum = 0;

   tmp = 0;

   for( i = 0; i < consdata->nvars; ++i )
   {
      assert(0 <= consdata->cliquepartition[i] && consdata->cliquepartition[i] < consdata->nvars);
      if( !cliqueused[consdata->cliquepartition[i]] )
      {
         maxcliqueweightsum += consdata->weights[i];
         cliqueused[consdata->cliquepartition[i]] = TRUE;
         tmpindices[tmp] = consdata->cliquepartition[i];
         ++tmp;
      }
   }
   /* clear cliqueused */
   for( --tmp; tmp >= 0; --tmp)
      cliqueused[tmp] = FALSE;

   assert(conshdlrdata->bools4size > 0);

   /* next if condition should normally not be true, because it means that presolving has created more binary variables
    * in one constraint than binary + integer variables existed in the whole problem at the presolving initialization method, but
    * for example if you would transform all integers into their binary representation then it maybe happens
    */ 
   if( conshdlrdata->bools4size < consdata->nvars )
   {
      int oldsize;
      oldsize = conshdlrdata->bools4size;

      while( conshdlrdata->bools4size < consdata->nvars )
         conshdlrdata->bools4size *= 2;
      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->bools4, conshdlrdata->bools4size) );
      BMSclearMemoryArray(&conshdlrdata->bools4[oldsize], conshdlrdata->bools4size - oldsize);
   }

   itemremoved = conshdlrdata->bools4;

   /* check for cleared array, all entries are zero */
#ifndef NDEBUG
   for( tmp = consdata->nvars - 1; tmp >= 0; --tmp )
      assert(itemremoved[tmp] == 0);
#endif

   /* for each binary variable xi and each fixing v, calculate the cliqueweightsum and update the weight of the
    * variable in the knapsack (this is sequence-dependent because the new or modified weights have to be
    * included in subsequent cliqueweightsum calculations)
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &addvars, 2*nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &addweights, 2*nbinvars) );
   naddvars = 0;
   addweightsum = 0;
   for( val = 0; val < 2 && addweightsum < consdata->capacity; ++val )
   {
      for( i = 0; i < nliftcands[val] && addweightsum < consdata->capacity; ++i )
      {
         SCIP_Longint cliqueweightsum;
         int probindex;
         int idx;
         int j;

	 tmp = 0;

         probindex = liftcands[val][i];
         assert(0 <= probindex && probindex < nbinvars);

         /* ignore empty zero lists and variables that cannot be lifted anyways */
         if( firstidxs[val][probindex] == 0
            || maxcliqueweightsum - zeroweightsums[val][probindex] + addweightsum >= consdata->capacity )
            continue;

         /* mark the items that are implied to zero by setting the current variable to the current value */
         for( idx = firstidxs[val][probindex]; idx != 0; idx = nextidxs[idx] )
         {
            assert(0 < idx && idx < nzeroitems);
            assert(0 <= zeroitems[idx] && zeroitems[idx] < consdata->nvars);
            itemremoved[zeroitems[idx]] = TRUE;
         }

         /* calculate the residual cliqueweight sum */
         cliqueweightsum = addweightsum; /* the previously added items are single-element cliques */
         for( j = 0; j < consdata->nvars; ++j )
         {
            assert(0 <= consdata->cliquepartition[j] && consdata->cliquepartition[j] < consdata->nvars);
            if( !itemremoved[j] && !cliqueused[consdata->cliquepartition[j]] )
            {
               cliqueweightsum += consdata->weights[j];
               cliqueused[consdata->cliquepartition[j]] = TRUE;
	       tmpindices[tmp] = consdata->cliquepartition[j];
               ++tmp;
               if( cliqueweightsum >= consdata->capacity )
                  break;
            }
         }

         /* check if the weight of the variable/value can be increased */
         if( cliqueweightsum < consdata->capacity )
         {
            SCIP_VAR* var;
            SCIP_Longint weight;

            /* insert the variable (with value TRUE) in the list of additional items */
            assert(naddvars < 2*nbinvars);
            var = binvars[probindex];
            if( val == FALSE )
            {
               SCIP_CALL( SCIPgetNegatedVar(scip, var, &var) );
            }
            weight = consdata->capacity - cliqueweightsum;
            addvars[naddvars] = var;
            addweights[naddvars] = weight;
            addweightsum += weight;
            naddvars++;
            SCIPdebugMessage("knapsack constraint <%s>: adding lifted item %"SCIP_LONGINT_FORMAT"<%s>\n",
               SCIPconsGetName(cons), weight, SCIPvarGetName(var));
         }

	 /* clear itemremoved */
         for( idx = firstidxs[val][probindex]; idx != 0; idx = nextidxs[idx] )
         {
            assert(0 < idx && idx < nzeroitems);
            assert(0 <= zeroitems[idx] && zeroitems[idx] < consdata->nvars);
            itemremoved[zeroitems[idx]] = FALSE;
	 }
	 /* clear cliqueused */
	 for( --tmp; tmp >= 0; --tmp)
	   cliqueused[tmpindices[tmp]] = FALSE;
      }
   }
   SCIPfreeBufferArray(scip, &tmpindices);

   /* clear part of zeroweightsums */
   for( --tmp3; tmp3 >= 0; --tmp3)
      zeroweightsums[tmpboolindices3[tmp3]][tmpindices3[tmp3]] = 0;

   /* clear rest of zeroweightsums and firstidxs */
   for( --tmp2; tmp2 >= 0; --tmp2)
   {
      zeroweightsums[tmpboolindices2[tmp2]][tmpindices2[tmp2]] = 0;
      firstidxs[tmpboolindices2[tmp2]][tmpindices2[tmp2]] = 0;
   }  
   
   SCIPfreeBufferArray(scip, &tmpindices2);
   SCIPfreeBufferArray(scip, &tmpindices3);
   SCIPfreeBufferArray(scip, &tmpboolindices2);
   SCIPfreeBufferArray(scip, &tmpboolindices3);

   /* add all additional item weights */
   for( i = 0; i < naddvars; ++i )
   {
      SCIP_CALL( addCoef(scip, cons, addvars[i], addweights[i]) );
   }
   *nchgcoefs += naddvars;

   /* if new items were added, multiple entries of the same variable are possible and we have to clean up the constraint */
   SCIP_CALL( mergeMultiples(scip, cons) );

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &addweights);
   SCIPfreeBufferArray(scip, &addvars);
   SCIPfreeBufferArray(scip, &nextidxs);
   SCIPfreeBufferArray(scip, &zeroitems);
   SCIPfreeBufferArray(scip, &liftcands[1]);
   SCIPfreeBufferArray(scip, &liftcands[0]);

   return SCIP_OKAY;
}

/** tightens item weights and capacity in presolving:
 *  given a knapsack sum(wi*xi) <= capacity
 *  (1) let weightsum := sum(wi)
 *      if weightsum - wi < capacity:
 *      - not using item i would make the knapsack constraint redundant
 *      - wi and capacity can be changed to have the same redundancy effect and the same results for
 *        fixing xi to zero or one, but with a reduced wi and tightened capacity to tighten the LP relaxation
 *      - change coefficients:
 *          wi'       := weightsum - capacity
 *          capacity' := capacity - (wi - wi')
 *  (2) let W(C) be the maximal weight of clique C,
 *          cliqueweightsum := sum(W(C))
 *      if cliqueweightsum - W(C) < capacity:
 *      - not using any item of C would make the knapsack constraint redundant
 *      - weights wi, i in C, and capacity can be changed to have the same redundancy effect and the same results for
 *        fixing xi, i in C, to zero or one, but with a reduced wi and tightened capacity to tighten the LP relaxation
 *      - change coefficients:
 *          delta     := capacity - (cliqueweightsum - W(C))
 *          wi'       := max(wi - delta, 0)
 *          capacity' := capacity - delta
 *      This rule has to add the used cliques in order to ensure they are enforced - otherwise, the reduction might
 *      introduce infeasible solutions.
 *  (3) for a clique C let C(xi == v) := C \ {j: xi == v -> xj == 0}),
 *      let cliqueweightsum(xi == v) := sum(W(C(xi == v)))
 *      if cliqueweightsum(xi == v) < capacity:
 *      - fixing variable xi to v would make the knapsack constraint redundant
 *      - the weight of the variable or its negation (depending on v) can be increased as long as it has the same
 *        redundancy effect:
 *          wi'       := capacity - cliqueweightsum(xi == v)
 *      This rule can also be applied to binary variables not in the knapsack!
 *  (4) if min{w} + wi > capacity:
 *      - using item i would force to fix other items to zero
 *      - wi can be increased to the capacity
 */
static
SCIP_RETCODE tightenWeights(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int*                  nchgcoefs,          /**< pointer to count total number of changed coefficients */
   int*                  nchgsides           /**< pointer to count number of side changes */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Longint minweight;
   int i;

   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(!SCIPconsIsModifiable(cons));

   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL); /* we are in presolve, so no LP row exists */
   assert(consdata->onesweightsum == 0); /* all fixed variables should have been removed */
   assert(consdata->weightsum > consdata->capacity); /* otherwise, the constraint is redundant */
   assert(consdata->nvars > 0);


   SCIP_CALL( mergeMultiples(scip, cons) );

   /* apply rule (1) */
   do
   {
      assert(consdata->merged);

      /* sort items, s.t. the heaviest one is in the first position */
      sortItems(consdata);

      for( i = 0; i < consdata->nvars; ++i )
      {
         SCIP_Longint weight;
         
         weight = consdata->weights[i];
         if( consdata->weightsum - weight < consdata->capacity )
         {
            SCIP_Longint newweight;
            
            newweight = consdata->weightsum - consdata->capacity;
            consdataChgWeight(consdata, i, newweight);
            consdata->capacity -= (weight - newweight);
            (*nchgcoefs)++;
            (*nchgsides)++;
            assert(!consdata->sorted);
            SCIPdebugMessage("knapsack constraint <%s>: changed weight of <%s> from %"SCIP_LONGINT_FORMAT" to %"SCIP_LONGINT_FORMAT", capacity from %"SCIP_LONGINT_FORMAT" to %"SCIP_LONGINT_FORMAT"\n",
               SCIPconsGetName(cons), SCIPvarGetName(consdata->vars[i]), weight, newweight,
               consdata->capacity + (weight-newweight), consdata->capacity);
         }
         else
            break;
      }
   }
   while( !consdata->sorted && consdata->weightsum > consdata->capacity );

   /* check for redundancy */
   if( consdata->weightsum <= consdata->capacity )
      return SCIP_OKAY;

   /* apply rule (2) (don't apply, if the knapsack has too many items for applying this costly method) */
   if( consdata->nvars <= MAX_USECLIQUES_SIZE )
   {
      SCIP_Longint* maxcliqueweights;
      SCIP_Longint* newweightvals;
      int* newweightidxs;
      SCIP_Longint cliqueweightsum;

      SCIP_CALL( SCIPallocBufferArray(scip, &maxcliqueweights, consdata->nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &newweightvals, consdata->nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &newweightidxs, consdata->nvars) );

      /* repeat as long as changes have been applied */
      do
      {
         int ncliques;
         SCIP_Bool zeroweights;

         assert(consdata->merged);

         /* sort items, s.t. the heaviest one is in the first position */
         sortItems(consdata);

         /* calculate a clique partition */
         SCIP_CALL( calcCliquepartition(scip, consdata) );

         /* if there are only single element cliques, rule (2) is equivalent to rule (1) */
         if( consdata->cliquepartition[consdata->nvars-1] == consdata->nvars-1 )
            break;

         /* calculate the maximal weight of the cliques */
         cliqueweightsum = 0;
         ncliques = 0;
         for( i = 0; i < consdata->nvars; ++i )
         {
            SCIP_Longint weight;
            int cliquenum;

            cliquenum = consdata->cliquepartition[i];
            weight = consdata->weights[i];
            assert(cliquenum <= ncliques);
            assert(weight > 0);
            if( cliquenum == ncliques )
            {
               maxcliqueweights[ncliques] = weight;
               cliqueweightsum += weight;
               ncliques++;
            }
            assert(maxcliqueweights[cliquenum] >= weight);
         }

         /* apply rule on every clique */
         zeroweights = FALSE;
         for( i = 0; i < ncliques; ++i )
         {
            SCIP_Longint delta;

            delta = consdata->capacity - (cliqueweightsum - maxcliqueweights[i]);
            if( delta > 0 )
            {
               SCIP_Longint newcapacity;
               SCIP_Longint newmincliqueweight;
               SCIP_Bool forceclique;
               int nnewweights;
               int j;

               SCIPdebugMessage("knapsack constraint <%s>: weights of clique %d (maxweight: %"SCIP_LONGINT_FORMAT") can be tightened: cliqueweightsum=%"SCIP_LONGINT_FORMAT", capacity=%"SCIP_LONGINT_FORMAT" -> delta: %"SCIP_LONGINT_FORMAT"\n",
                  SCIPconsGetName(cons), i, maxcliqueweights[i], cliqueweightsum, consdata->capacity, delta);
               newcapacity = consdata->capacity - delta;
               newmincliqueweight = newcapacity + 1;
               forceclique = FALSE;
               nnewweights = 0;
#ifndef NDEBUG
               for( j = 0; j < i; ++j )
                  assert(consdata->cliquepartition[j] < i); /* no element j < i can be in clique i */
#endif
               for( j = i; j < consdata->nvars; ++j )
               {
                  if( consdata->cliquepartition[j] == i )
                  {
                     SCIP_Longint newweight;
                  
                     newweight = consdata->weights[j] - delta;
                     newweight = MAX(newweight, 0);

                     /* check if this new weight and the minimal new weight of already processed items both fit into the
                      * knapsack; if this is true, we have to add a clique constraint in order to enforce the clique
                      * (otherwise, the knapsack might have been one of the reasons for the clique, and the weight
                      * reduction might be infeasible, i.e., allows additional solutions)
                      */
                     if( newmincliqueweight + newweight <= consdata->capacity )
                     {
                        forceclique = TRUE;
                        if( !conshdlrdata->disaggregation )
                           break;
                     }

                     /* cache the new weight */
                     assert(nnewweights < consdata->nvars);
                     newweightvals[nnewweights] = newweight;
                     newweightidxs[nnewweights] = j;
                     nnewweights++;

                     assert(newweight <= newmincliqueweight); /* items are sorted by non-increasing weight! */
                     newmincliqueweight = newweight;
                  }
               }

               /* check if we really want to apply the change */
               if( conshdlrdata->disaggregation || !forceclique )
               {
                  int k;

                  SCIPdebugMessage(" -> change capacity from %"SCIP_LONGINT_FORMAT" to %"SCIP_LONGINT_FORMAT" (forceclique:%u)\n",
                     consdata->capacity, newcapacity, forceclique);
                  consdata->capacity = newcapacity;
                  (*nchgsides)++;

                  for( k = 0; k < nnewweights; ++k )
                  {
                     j = newweightidxs[k];
                     assert(0 <= j && j < consdata->nvars);
                     assert(consdata->cliquepartition[j] == i);

                     /* apply the weight change */
                     SCIPdebugMessage(" -> change weight of <%s> from %"SCIP_LONGINT_FORMAT" to %"SCIP_LONGINT_FORMAT"\n",
                        SCIPvarGetName(consdata->vars[j]), consdata->weights[j], newweightvals[k]);
                     consdataChgWeight(consdata, j, newweightvals[k]);
                     (*nchgcoefs)++;
                     assert(!consdata->sorted);
                     zeroweights = zeroweights || (newweightvals[k] == 0);
                  }

                  /* if before the weight update at least one pair of weights did not fit into the knapsack and now fits,
                   * we have to make sure, the clique is enforced - the clique might have been constructed partially from
                   * this constraint, and by reducing the weights, this clique information is not contained anymore in the
                   * knapsack constraint
                   */
                  if( forceclique )
                  {
                     SCIP_CONS* cliquecons;
                     char name[SCIP_MAXSTRLEN];
                     SCIP_VAR** cliquevars;

                     SCIP_CALL( SCIPallocBufferArray(scip, &cliquevars, nnewweights) );
                     for( k = 0; k < nnewweights; ++k )
                        cliquevars[k] = consdata->vars[newweightidxs[k]];

                     (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_clq_%"SCIP_LONGINT_FORMAT"_%d", SCIPconsGetName(cons), consdata->capacity, i);
                     SCIP_CALL( SCIPcreateConsSetpack(scip, &cliquecons, name, nnewweights, cliquevars,
                           SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                           SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
                           SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
                           SCIPconsIsStickingAtNode(cons)) );
                     SCIPdebugMessage(" -> adding clique constraint: ");
                     SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cliquecons, NULL) ) );
                     SCIP_CALL( SCIPaddCons(scip, cliquecons) );
                     SCIP_CALL( SCIPreleaseCons(scip, &cliquecons) );
                     SCIPfreeBufferArray(scip, &cliquevars);
                  }
               }
            }
         }
         if( zeroweights )
         {
            SCIP_CALL( removeZeroWeights(scip, cons) );
         }
      }
      while( !consdata->sorted && consdata->weightsum > consdata->capacity );

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &newweightidxs);
      SCIPfreeBufferArray(scip, &newweightvals);
      SCIPfreeBufferArray(scip, &maxcliqueweights);

      /* check for redundancy */
      if( consdata->weightsum <= consdata->capacity )
         return SCIP_OKAY;
   }

   /* apply rule (3) */
   SCIP_CALL( tightenWeightsLift(scip, cons, nchgcoefs) );

   /* check for redundancy */
   if( consdata->weightsum <= consdata->capacity )
      return SCIP_OKAY;

   /* apply rule (4) (all but smallest weight) */
   assert(consdata->merged);
   sortItems(consdata);
   minweight = consdata->weights[consdata->nvars-1];
   for( i = 0; i < consdata->nvars-1; ++i )
   {
      SCIP_Longint weight;

      weight = consdata->weights[i];
      assert(weight >= minweight);
      if( minweight + weight > consdata->capacity )
      {
         if( weight < consdata->capacity )
         {
            SCIPdebugMessage("knapsack constraint <%s>: changed weight of <%s> from %"SCIP_LONGINT_FORMAT" to %"SCIP_LONGINT_FORMAT"\n",
               SCIPconsGetName(cons), SCIPvarGetName(consdata->vars[i]), weight, consdata->capacity);
            assert(consdata->sorted);
            consdataChgWeight(consdata, i, consdata->capacity); /* this does not destroy the weight order! */
            assert(i == 0 || consdata->weights[i-1] >= consdata->weights[i]);
            consdata->sorted = TRUE;
            (*nchgcoefs)++;
         }
      }
      else
         break;
   }

   /* apply rule (4) (smallest weight) */
   if( consdata->nvars >= 2 )
   {
      SCIP_Longint weight;

      minweight = consdata->weights[consdata->nvars-2];
      weight = consdata->weights[consdata->nvars-1];
      assert(minweight >= weight);
      if( minweight + weight > consdata->capacity && weight < consdata->capacity )
      {
         SCIPdebugMessage("knapsack constraint <%s>: changed weight of <%s> from %"SCIP_LONGINT_FORMAT" to %"SCIP_LONGINT_FORMAT"\n",
            SCIPconsGetName(cons), SCIPvarGetName(consdata->vars[consdata->nvars-1]), weight, consdata->capacity);
         assert(consdata->sorted);
         consdataChgWeight(consdata, consdata->nvars-1, consdata->capacity); /* this does not destroy the weight order! */
         assert(minweight >= consdata->weights[consdata->nvars-1]);
         consdata->sorted = TRUE;
         (*nchgcoefs)++;
      }
   }

   return SCIP_OKAY;
}

/** adds cliques of the knapsack constraint to the global clique table */
static
SCIP_RETCODE addCliques(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nbdchgs             /**< pointer to count the number of performed bound changes */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** cliquevars;
   int ncliquevars;
   int i;

   assert(cutoff != NULL);
   assert(nbdchgs != NULL);

   *cutoff = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check whether the cliques have already been added */
   if( consdata->cliquesadded || consdata->nvars == 0 )
      return SCIP_OKAY;

   /* make sure, the items are merged and sorted by non-increasing weight */
   SCIP_CALL( mergeMultiples(scip, cons) );
   sortItems(consdata);

   /* get temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquevars, consdata->nvars) );

   /* build a largest clique by using the items with the maximal weights */
   cliquevars[0] = consdata->vars[0];
   for( i = 1; i < consdata->nvars && consdata->weights[i-1] + consdata->weights[i] > consdata->capacity; ++i )
      cliquevars[i] = consdata->vars[i];
   ncliquevars = i;

   if( ncliquevars >= 2 )
   {
      SCIP_Longint cliqueminweight;
      int thisnbdchgs;

      /* add the clique to the clique table */
      SCIP_CALL( SCIPaddClique(scip, cliquevars, NULL, ncliquevars, cutoff, &thisnbdchgs) );
      *nbdchgs += thisnbdchgs;

      /* try to replace the last item in the clique by a different item to obtain a slightly different clique */
      cliqueminweight = consdata->weights[ncliquevars-2];
      for( i = ncliquevars; i < consdata->nvars && cliqueminweight + consdata->weights[i] > consdata->capacity; ++i )
      {
         cliquevars[ncliquevars-1] = consdata->vars[i];
         SCIP_CALL( SCIPaddClique(scip, cliquevars, NULL, ncliquevars, cutoff, &thisnbdchgs) );
         *nbdchgs += thisnbdchgs;
      }
   }

   /* free temporary memory and mark the constraint */
   SCIPfreeBufferArray(scip, &cliquevars);
   consdata->cliquesadded = TRUE;

   return SCIP_OKAY;
}




/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeMemory(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int nvars;

   assert( scip != NULL );
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   nvars = SCIPgetNBinVars(scip);
   nvars += SCIPgetNIntVars(scip); /* maybe some integer variables get upgraded to binary variables due to presolving, so we need to count them too */

   SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->reals1, nvars) );
   BMSclearMemoryArray(conshdlrdata->reals1, nvars);
   conshdlrdata->reals1size = nvars;

   return SCIP_OKAY;
}

/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeMemoryArrayNull(scip, &conshdlrdata->reals1);
   conshdlrdata->reals1size = 0;

   return SCIP_OKAY;
}


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int nvars;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( nconss == 0 || conss != NULL );
   assert( result != NULL );

   *result = SCIP_FEASIBLE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   nvars = SCIPgetNBinVars(scip);
   nvars += SCIPgetNIntVars(scip); /* maybe some integer variables get upgraded to binary variables due to presolving, so we need to count them too */

   SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->ints1, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->ints2, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->longints1, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->longints2, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->bools1, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->bools2, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->bools3, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->bools4, nvars) );

   BMSclearMemoryArray(conshdlrdata->ints1, nvars);
   BMSclearMemoryArray(conshdlrdata->ints2, nvars);
   BMSclearMemoryArray(conshdlrdata->longints1, nvars);
   BMSclearMemoryArray(conshdlrdata->longints2, nvars);
   BMSclearMemoryArray(conshdlrdata->bools1, nvars);
   BMSclearMemoryArray(conshdlrdata->bools2, nvars);
   BMSclearMemoryArray(conshdlrdata->bools3, nvars);
   BMSclearMemoryArray(conshdlrdata->bools4, nvars);

   conshdlrdata->ints1size = nvars;
   conshdlrdata->ints2size = nvars;
   conshdlrdata->longints1size = nvars;
   conshdlrdata->longints2size = nvars;
   conshdlrdata->bools1size = nvars;
   conshdlrdata->bools2size = nvars;
   conshdlrdata->bools3size = nvars;
   conshdlrdata->bools4size = nvars;

   return SCIP_OKAY;
}


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   
   assert( scip != NULL );
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeMemoryArrayNull(scip, &conshdlrdata->ints1);
   SCIPfreeMemoryArrayNull(scip, &conshdlrdata->ints2);
   SCIPfreeMemoryArrayNull(scip, &conshdlrdata->longints1);
   SCIPfreeMemoryArrayNull(scip, &conshdlrdata->longints2);
   SCIPfreeMemoryArrayNull(scip, &conshdlrdata->bools1);
   SCIPfreeMemoryArrayNull(scip, &conshdlrdata->bools2);
   SCIPfreeMemoryArrayNull(scip, &conshdlrdata->bools3);
   SCIPfreeMemoryArrayNull(scip, &conshdlrdata->bools4);

   conshdlrdata->ints1size = 0;
   conshdlrdata->ints2size = 0;
   conshdlrdata->longints1size = 0;
   conshdlrdata->longints2size = 0;
   conshdlrdata->bools1size = 0;
   conshdlrdata->bools2size = 0;
   conshdlrdata->bools3size = 0;
   conshdlrdata->bools4size = 0;

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolKnapsack NULL


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert( scip != NULL );

   /* release the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->row != NULL )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->row) );
      }
   }

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);
   
   /* free knapsack constraint */
   SCIP_CALL( consdataFree(scip, consdata, conshdlrdata->eventhdlr) );
   
   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */ 
static
SCIP_DECL_CONSTRANS(consTransKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

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
   SCIP_CALL( consdataCreate(scip, &targetdata, conshdlrdata->eventhdlr,
         sourcedata->nvars, sourcedata->vars, sourcedata->weights, sourcedata->capacity) ); 

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), 
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpKnapsack)
{  /*lint --e{715}*/
   int i;

   for( i = 0; i < nconss; i++ )
   {
      assert(SCIPconsIsInitial(conss[i]));
      SCIP_CALL( addRelaxation(scip, conss[i], NULL) );
   }

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool sepacardinality;
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
   
   SCIPdebugMessage("knapsack separation of %d/%d constraints, round %d (max %d/%d)\n",
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

   /* check dual bound to see if we want to produce knapsack cardinality cuts at this node */
   if( sepacardinality )
   {
      SCIP_Real loclowerbound;
      SCIP_Real glblowerbound;
      SCIP_Real cutoffbound;
      SCIP_Real maxbound;

      loclowerbound = SCIPgetLocalLowerbound(scip);
      glblowerbound = SCIPgetLowerbound(scip);
      cutoffbound = SCIPgetCutoffbound(scip);
      maxbound = glblowerbound + conshdlrdata->maxcardbounddist * (cutoffbound - glblowerbound);
      sepacardinality = SCIPisLE(scip, loclowerbound, maxbound);
   }

   /* get the maximal number of cuts allowed in a separation round */
   maxsepacuts = (depth == 0 ? conshdlrdata->maxsepacutsroot : conshdlrdata->maxsepacuts);

   *result = SCIP_DIDNOTFIND;
   ncuts = 0;

   /* separate useful constraints */
   for( i = 0; i < nusefulconss && ncuts < maxsepacuts && !SCIPisStopped(scip); i++ )
   {
      SCIP_CALL( separateCons(scip, conss[i], NULL, sepacardinality, conshdlrdata->maxnumcardlift, &ncuts) );
   }
   
   /* adjust return value */
   if( ncuts > 0 )
      *result = SCIP_SEPARATED;
   
   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool sepacardinality;
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
   
   SCIPdebugMessage("knapsack separation of %d/%d constraints, round %d (max %d/%d)\n",
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
   for( i = 0; i < nusefulconss && ncuts < maxsepacuts && !SCIPisStopped(scip); i++ )
   {
      SCIP_CALL( separateCons(scip, conss[i], sol, sepacardinality, conshdlrdata->maxnumcardlift, &ncuts) );
   }
   
   /* adjust return value */
   if( ncuts > 0 )
      *result = SCIP_SEPARATED;
   
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool violated;
   int maxncuts;
   int ncuts;
   int i;

   *result = SCIP_FEASIBLE;

   SCIPdebugMessage("knapsack enforcement of %d/%d constraints\n", nusefulconss, nconss);

   /* get maximal number of cuts per round */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   maxncuts = (SCIPgetDepth(scip) == 0 ? conshdlrdata->maxsepacutsroot : conshdlrdata->maxsepacuts);

   ncuts = 0;

   /* search for violated useful knapsack constraints */
   for( i = 0; i < nusefulconss && ncuts < maxncuts; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], NULL, FALSE, FALSE, &violated) );
      if( violated )
      {
         /* add knapsack constraint as LP row to the LP */
         SCIP_CALL( addRelaxation(scip, conss[i], NULL) );
         ncuts++;
      }
   } 

   /* as long as no violations were found, search for violated obsolete knapsack constraints */
   for( i = nusefulconss; i < nconss && ncuts == 0; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], NULL, FALSE, FALSE, &violated) );
      if( violated )
      {
         /* add knapsack constraint as LP row to the LP */
         SCIP_CALL( addRelaxation(scip, conss[i], NULL) );
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
SCIP_DECL_CONSENFOPS(consEnfopsKnapsack)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int i;

   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], NULL, TRUE, FALSE, &violated) );
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
SCIP_DECL_CONSCHECK(consCheckKnapsack)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int i;

   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], sol, checklprows, printreason, &violated) );
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
SCIP_DECL_CONSPROP(consPropKnapsack)
{  /*lint --e{715}*/
   SCIP_Bool cutoff;
   SCIP_Bool redundant;
   int nfixedvars;
   int i;

   cutoff = FALSE;
   nfixedvars = 0;

   /* process useful constraints */
   for( i = 0; i < nusefulconss && !cutoff; i++ )
   {
      SCIP_CALL( propagateCons(scip, conss[i], &cutoff, &redundant, &nfixedvars) );
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
SCIP_DECL_CONSPRESOL(consPresolKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Bool cutoff;
   SCIP_Bool redundant;
   int oldnfixedvars;
   int oldnchgbds;
   int oldndelconss;
   int oldnchgcoefs;
   int oldnchgsides;
   int i;

   /* remember old preprocessing counters */
   cutoff = FALSE;
   oldnfixedvars = *nfixedvars;
   oldnchgbds = *nchgbds;
   oldndelconss = *ndelconss;
   oldnchgcoefs = *nchgcoefs;
   oldnchgsides = *nchgsides;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   
   for( i = 0; i < nconss && !SCIPisStopped(scip); i++ )
   {
      SCIP_CONS* cons;
      int thisnfixedvars;
      int thisnchgbds;

      cons = conss[i];
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* force presolving the constraint in the initial round */
      if( nrounds == 0 )
         consdata->presolved = FALSE;

      if( consdata->presolved )
         continue;

      SCIPdebugMessage("presolving knapsack constraint <%s>\n", SCIPconsGetName(cons));
      SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );

      consdata->presolved = TRUE;

      /* remove all fixed variables */
      if( nrounds == 0 || nnewfixedvars > 0 || nnewaggrvars > 0 || nnewchgbds > 0
         || *nfixedvars > oldnfixedvars || *nchgbds > oldnchgbds )
      {
         SCIP_CALL( applyFixings(scip, cons, &cutoff) );
         if( cutoff )
            break;
      }
      thisnfixedvars = *nfixedvars;
      thisnchgbds = *nchgbds;

      /* add cliques in the knapsack to the clique table */
      SCIP_CALL( addCliques(scip, cons, &cutoff, nchgbds) );
      if( cutoff )
         break;

      /* propagate constraint */
      SCIP_CALL( propagateCons(scip, cons, &cutoff, &redundant, nfixedvars) );
      if( cutoff )
         break;
      if( redundant )
      {
         (*ndelconss)++;
         continue;
      }

      /* remove again all fixed variables, if further fixings were found */
      if( *nfixedvars > thisnfixedvars || *nchgbds > thisnchgbds )
      {
         SCIP_CALL( applyFixings(scip, cons, &cutoff) );
         if( cutoff )
            break;
      }

      if( !SCIPconsIsModifiable(cons) )
      {
         /* check again for redundancy (applyFixings() might have decreased weightsum due to fixed-to-zero vars) */
         if( consdata->weightsum <= consdata->capacity )
         {
            SCIPdebugMessage(" -> knapsack constraint <%s> is redundant: weightsum=%"SCIP_LONGINT_FORMAT", capacity=%"SCIP_LONGINT_FORMAT"\n",
               SCIPconsGetName(cons), consdata->weightsum, consdata->capacity);
            SCIP_CALL( SCIPdelConsLocal(scip, cons) );
            continue;
         }

         /* divide weights by their greatest common divisor */
         normalizeWeights(cons, nchgcoefs, nchgsides);

         /* tighten capacity and weights */
         SCIP_CALL( tightenWeights(scip, cons, nchgcoefs, nchgsides) );

         /* try to simplify inequalities */
         if( conshdlrdata->simplifyinequalities )
         {
	   SCIP_CALL( simplifyInequalities(scip, cons, ndelconss, nchgcoefs, nchgsides, &cutoff) );
	   if( cutoff )
	      break;
         }
      }
   } 

   /**@todo preprocess pairs of knapsack constraints */

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( *nfixedvars > oldnfixedvars || *nchgbds > oldnchgbds || *ndelconss > oldndelconss
      || *nchgcoefs > oldnchgcoefs || *nchgsides > oldnchgsides )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Longint capsum;
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
            SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[i]) );
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
SCIP_DECL_CONSLOCK(consLockKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   for( i = 0; i < consdata->nvars; i++)
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlocksneg, nlockspos) );
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
SCIP_DECL_CONSPRINT(consPrintKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   for( i = 0; i < consdata->nvars; ++i )
   {
      if( i > 0 )
         SCIPinfoMessage(scip, file, " ");
      SCIPinfoMessage(scip, file, "%+"SCIP_LONGINT_FORMAT"<%s>", consdata->weights[i], SCIPvarGetName(consdata->vars[i]));
   }
   SCIPinfoMessage(scip, file, " <= %"SCIP_LONGINT_FORMAT"", consdata->capacity);
   
   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyKnapsack)
{  /*lint --e{715}*/
   SCIP_VAR** sourcevars;
   SCIP_Longint* weights;
   SCIP_Real* coefs;
   const char* consname;
   int nvars;
   int v;
   
   /* get variables and coefficients of the source constraint */
   sourcevars = SCIPgetVarsKnapsack(sourcescip, sourcecons);
   nvars = SCIPgetNVarsKnapsack(sourcescip, sourcecons);
   weights = SCIPgetWeightsKnapsack(sourcescip, sourcecons);

   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
   for( v = 0; v < nvars; ++v )
      coefs[v] = (SCIP_Real) weights[v];
   
   if( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   /* copy the logic using the linear constraint copy method */
   SCIP_CALL( SCIPcopyConsLinear(scip, cons, sourcescip, consname, nvars, sourcevars, coefs,
         -SCIPinfinity(scip), (SCIP_Real) SCIPgetCapacityKnapsack(sourcescip, sourcecons), varmap,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, success) );

   SCIPfreeBufferArray(scip, &coefs);

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
#define consParseKnapsack NULL



/*
 * Linear constraint upgrading
 */

/** creates and captures a knapsack constraint out of a linear inequality */
static
SCIP_RETCODE createNormalizedKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            vals,               /**< array with inequality coefficients */
   SCIP_Real             lhs,                /**< left hand side of inequality */
   SCIP_Real             rhs,                /**< right hand side of inequality */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? 
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? 
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? 
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally? 
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? 
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which 
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? 
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_VAR** transvars;
   SCIP_Longint* weights;
   SCIP_Longint capacity;
   SCIP_Longint weight;
   int mult;
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);
   assert(SCIPisInfinity(scip, -lhs) != SCIPisInfinity(scip, rhs));

   /* get temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &transvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &weights, nvars) );

   /* if the right hand side is non-infinite, we have to negate all variables with negative coefficient;
    * otherwise, we have to negate all variables with positive coefficient and multiply the row with -1
    */
   if( SCIPisInfinity(scip, rhs) )
   {
      mult = -1;
      capacity = (SCIP_Longint)SCIPfeasFloor(scip, -lhs);
   }
   else
   {
      mult = +1;
      capacity = (SCIP_Longint)SCIPfeasFloor(scip, rhs);
   }

   /* negate positive or negative variables */
   for( v = 0; v < nvars; ++v )
   {
      assert(SCIPisFeasIntegral(scip, vals[v]));
      weight = mult * (SCIP_Longint)SCIPfeasFloor(scip, vals[v]);
      if( weight > 0 )
      {
         transvars[v] = vars[v];
         weights[v] = weight;
      }
      else
      {
         SCIP_CALL( SCIPgetNegatedVar(scip, vars[v], &transvars[v]) );
         weights[v] = -weight;
         capacity -= weight;
      }
      assert(transvars[v] != NULL);
   }

   /* create the constraint */
   SCIP_CALL( SCIPcreateConsKnapsack(scip, cons, name, nvars, transvars, weights, capacity,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &weights);
   SCIPfreeBufferArray(scip, &transvars);

   return SCIP_OKAY;
}

/** tries to upgrade a linear constraint into a knapsack constraint */
static
SCIP_DECL_LINCONSUPGD(linconsUpgdKnapsack)
{  /*lint --e{715}*/
   SCIP_Bool upgrade;

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
      SCIPdebugMessage("upgrading constraint <%s> to knapsack constraint\n", SCIPconsGetName(cons));
      
      /* create the knapsack constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      SCIP_CALL( createNormalizedKnapsack(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, lhs, rhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), 
            SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
   }

   return SCIP_OKAY;
}




/*
 * Event handler
 */

/** execution methode of bound change event handler */
static
SCIP_DECL_EVENTEXEC(eventExecKnapsack)
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
   case SCIP_EVENTTYPE_VARFIXED:  /* the variable should be removed from the constraint in presolving */
   case SCIP_EVENTTYPE_IMPLADDED: /* further preprocessing might be possible due to additional implications */
      eventdata->consdata->presolved = FALSE;
      break;
   default:
      SCIPerrorMessage("invalid event type %x\n", SCIPeventGetType(event));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for knapsack constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrKnapsack(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* include event handler for bound change events */
   eventhdlrdata = NULL;
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC, 
         NULL, NULL, NULL, NULL, NULL, NULL, eventExecKnapsack,
         eventhdlrdata) );

   /* create knapsack constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );

   /* get event handler for bound change events */
   conshdlrdata->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( conshdlrdata->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for knapsack constraints not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         consFreeKnapsack, consInitKnapsack, consExitKnapsack, 
         consInitpreKnapsack, consExitpreKnapsack, consInitsolKnapsack, consExitsolKnapsack,
         consDeleteKnapsack, consTransKnapsack, consInitlpKnapsack,
         consSepalpKnapsack, consSepasolKnapsack, consEnfolpKnapsack, consEnfopsKnapsack, consCheckKnapsack, 
         consPropKnapsack, consPresolKnapsack, consRespropKnapsack, consLockKnapsack,
         consActiveKnapsack, consDeactiveKnapsack, 
         consEnableKnapsack, consDisableKnapsack,
         consPrintKnapsack, consCopyKnapsack, consParseKnapsack,
         conshdlrdata) );
  
   /* include the linear constraint upgrade in the linear constraint handler */
   SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdKnapsack, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );

   /* add knapsack constraint handler parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/knapsack/sepacardfreq",
         "multiplier on separation frequency, how often cardinality cuts are separated (-1: never, 0: only at root)",
         &conshdlrdata->sepacardfreq, TRUE, DEFAULT_SEPACARDFREQ, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/knapsack/maxcardbounddist",
         "maximal relative distance from current node's dual bound to primal bound compared to best node's dual bound for separating knapsack cardinality cuts",
         &conshdlrdata->maxcardbounddist, TRUE, DEFAULT_MAXCARDBOUNDDIST, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/knapsack/maxrounds",
         "maximal number of separation rounds per node (-1: unlimited)",
         &conshdlrdata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/knapsack/maxroundsroot",
         "maximal number of separation rounds per node in the root node (-1: unlimited)",
         &conshdlrdata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/knapsack/maxsepacuts",
         "maximal number of cuts separated per separation round",
         &conshdlrdata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/knapsack/maxsepacutsroot",
         "maximal number of cuts separated per separation round in the root node",
         &conshdlrdata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/knapsack/maxnumcardlift",
         "maximal number of cardinality inequalities lifted per separation round (-1: unlimited)",
         &conshdlrdata->maxnumcardlift, TRUE, DEFAULT_MAXNUMCARDLIFT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/knapsack/disaggregation",
         "should disaggregation of knapsack constraints be allowed in preprocessing?",
         &conshdlrdata->disaggregation, TRUE, DEFAULT_DISAGGREGATION, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/knapsack/simplifyinequalities",
         "should presolving try to simplify knapsacks",
         &conshdlrdata->simplifyinequalities, TRUE, DEFAULT_SIMPLIFYINEQUALITIES, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a knapsack constraint */
SCIP_RETCODE SCIPcreateConsKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of items in the knapsack */
   SCIP_VAR**            vars,               /**< array with item variables */
   SCIP_Longint*         weights,            /**< array with item weights */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which 
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the knapsack constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("knapsack constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, conshdlrdata->eventhdlr, nvars, vars, weights, capacity) );
        
   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** adds new item to knapsack constraint */
SCIP_RETCODE SCIPaddCoefKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< item variable */
   SCIP_Longint          weight              /**< item weight */
   )
{
   assert(var != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a knapsack constraint\n");
      return SCIP_INVALIDDATA;
   }
   
   SCIP_CALL( addCoef(scip, cons, var, weight) );

   return SCIP_OKAY;
}

/** gets the capacity of the knapsack constraint */
SCIP_Longint SCIPgetCapacityKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a knapsack constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->capacity;
}

/** gets the number of items in the knapsack constraint */
int SCIPgetNVarsKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a knapsack constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nvars;
}

/** gets the array of variables in the knapsack constraint; the user must not modify this array! */
SCIP_VAR** SCIPgetVarsKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a knapsack constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vars;
}

/** gets the array of weights in the knapsack constraint; the user must not modify this array! */
SCIP_Longint* SCIPgetWeightsKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a knapsack constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->weights;
}

/** gets the dual solution of the knapsack constraint in the current LP */
SCIP_Real SCIPgetDualsolKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a knapsack constraint\n");
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
SCIP_Real SCIPgetDualfarkasKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a knapsack constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualfarkas(consdata->row);
   else
      return 0.0;
}

/** returns the linear relaxation of the given knapsack constraint; may return NULL if no LP row was yet created;
 *  the user must not modify the row!
 */
SCIP_ROW* SCIPgetRowKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a knapsack\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->row;
}


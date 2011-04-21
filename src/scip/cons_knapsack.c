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

#define DEFAULT_SEPACARDFREQ          1 /**< multiplier on separation frequency, how often knapsack cuts are separated */
#define DEFAULT_MAXROUNDS             5 /**< maximal number of separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        -1 /**< maximal number of separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS          50 /**< maximal number of cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT     200 /**< maximal number of cuts separated per separation round in the root node */
#define DEFAULT_MAXCARDBOUNDDIST    0.0 /**< maximal relative distance from current node's dual bound to primal bound compared
                                         *   to best node's dual bound for separating knapsack cuts */
#define DEFAULT_DISAGGREGATION     TRUE /**< should disaggregation of knapsack constraints be allowed in preprocessing? */
#define DEFAULT_SIMPLIFYINEQUALITIES FALSE/**< should presolving try to simplify knapsacks */
#define DEFAULT_NEGATEDCLIQUE      TRUE /**< should negated clique information be used in solving process */

#define MAXABSVBCOEF               1e+5 /**< maximal absolute coefficient in variable bounds used for knapsack relaxation */
#define USESUPADDLIFT             FALSE /**< should lifted minimal cover ineq using superadditive up-lifting be separated in addition */

#define DEFAULT_PRESOLUSEHASHING   TRUE /**< should hash table be used for detecting redundant constraints in advance */
#define HASHSIZE_KNAPSACKCONS    131101 /**< minimal size of hash table in linear constraint tables */

#define DEFAULT_PRESOLPAIRWISE     TRUE /**< should pairwise constraint comparison be performed in presolving? */
#define NMINCOMPARISONS          200000 /**< number for minimal pairwise presol comparisons */
#define MINGAINPERNMINCOMPARISONS 1e-06 /**< minimal gain per minimal pairwise presol comparisons to repeat pairwise 
                                         *   comparison round */
#define DEFAULT_DUALPRESOLVING     TRUE /**< should dual presolving steps be preformed? */


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
                                              *   to best node's dual bound for separating knapsack cuts */
   int                   sepacardfreq;       /**< multiplier on separation frequency, how often knapsack cuts are separated */
   int                   maxrounds;          /**< maximal number of separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of separation rounds in the root node (-1: unlimited) */
   int                   maxsepacuts;        /**< maximal number of cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of cuts separated per separation round in the root node */
   SCIP_Bool             disaggregation;     /**< should disaggregation of knapsack constraints be allowed in preprocessing? */
   SCIP_Bool             simplifyinequalities;/**< should presolving try to cancel down or delete coefficients in inequalities */
   SCIP_Bool             negatedclique;      /**< should negated clique information be used in solving process */
   SCIP_Bool             presolpairwise;     /**< should pairwise constraint comparison be performed in presolving? */
   SCIP_Bool             presolusehashing;   /**< should hash table be used for detecting redundant constraints in advance */
   SCIP_Bool             dualpresolving;     /**< should dual presolving steps be preformed? */
};


/** constraint data for knapsack constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars;               /**< variables in knapsack constraint */
   SCIP_Longint*         weights;            /**< weights of variables in knapsack constraint */
   SCIP_EVENTDATA**      eventdatas;         /**< event datas for bound change events of the variables */
   int*                  cliquepartition;    /**< clique indices of the clique partition */
   int*                  negcliquepartition; /**< clique indices of the negated clique partition */
   SCIP_ROW*             row;                /**< corresponding LP row */
   int                   nvars;              /**< number of variables in knapsack constraint */
   int                   varssize;           /**< size of vars, weights, and eventdatas arrays */
   int                   ncliques;           /**< number of cliques in the clique partition */
   int                   nnegcliques;        /**< number of cliques in the negated clique partition */
   SCIP_Longint          capacity;           /**< capacity of knapsack */
   SCIP_Longint          weightsum;          /**< sum of all weights */
   SCIP_Longint          onesweightsum;      /**< sum of weights of variables fixed to one */
   unsigned int          propagated:1;       /**< is the knapsack constraint already propagated? */
   unsigned int          presolved:1;        /**< is the knapsack constraint already presolved? */
   unsigned int          sorted:1;           /**< are the knapsack items sorted by weight? */
   unsigned int          cliquepartitioned:1;/**< is the clique partition valid? */
   unsigned int          negcliquepartitioned:1;/**< is the negated clique partition valid? */
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


/** data structure to combine two sorting key values */
struct sortkeypair
{
   SCIP_Real             key1;               /**< first sort key value */
   SCIP_Real             key2;               /**< second sort key value */
};
typedef struct sortkeypair SORTKEYPAIR;


/*
 * Local methods
 */

/** comparison method for two sorting key pairs */
static
SCIP_DECL_SORTPTRCOMP(compSortkeypairs)
{
   SORTKEYPAIR* sortkeypair1 = (SORTKEYPAIR*)elem1;
   SORTKEYPAIR* sortkeypair2 = (SORTKEYPAIR*)elem2;

   if( sortkeypair1->key1 < sortkeypair2->key1 )
      return -1;
   else if( sortkeypair1->key1 > sortkeypair2->key1 )
      return +1;
   else if( sortkeypair1->key2 < sortkeypair2->key2 )
      return -1;
   else if( sortkeypair1->key2 > sortkeypair2->key2 )
      return +1;
   else 
      return 0;
}

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
   assert(consdata->nvars == 0 || (consdata->cliquepartition != NULL && consdata->negcliquepartition != NULL));

   if( !consdata->sorted )
   {
      int pos;
      int lastcliquenum;
      int v;

      /* sort of five joint arrays of Long/pointer/pointer/ints/ints, 
       * sorted by first array in non-increasing order via sort template */
      SCIPsortDownLongPtrPtrIntInt(
         consdata->weights, 
         (void**)consdata->vars, 
         (void**)consdata->eventdatas, 
         consdata->cliquepartition,
         consdata->negcliquepartition,
         consdata->nvars);
      
      /* sort all items with same weight according to their variable index, used for hash value for fast pairwise comparison of all constraints */
      for(v = consdata->nvars - 2; v >= 0; --v )
      {
         SCIP_VAR* tmpvar;
         SCIP_EVENTDATA* tmpeventdata;
         int tmpclique;
         int tmpnegclique;

         if( consdata->weights[v] == consdata->weights[v+1]
            && SCIPvarGetIndex(consdata->vars[v]) > SCIPvarGetIndex(consdata->vars[v + 1]) )
         {
            tmpvar = consdata->vars[v];
            tmpeventdata = consdata->eventdatas[v];
            tmpclique = consdata->cliquepartition[v];
            tmpnegclique = consdata->negcliquepartition[v];
            
            consdata->vars[v] = consdata->vars[v + 1];
            consdata->eventdatas[v] = consdata->eventdatas[v + 1];
            consdata->cliquepartition[v] = consdata->cliquepartition[v + 1];
            consdata->negcliquepartition[v] = consdata->negcliquepartition[v + 1];
            
            consdata->vars[v + 1] = tmpvar;
            consdata->eventdatas[v + 1] = tmpeventdata;
            consdata->cliquepartition[v + 1] = tmpclique;
            consdata->negcliquepartition[v + 1] = tmpnegclique;
         }
      }

      /* we need to make sure that our clique numbers of our normal clique will be in increasing order without gaps */
      if( consdata->cliquepartitioned )
      {
	 lastcliquenum = 0;

         for( pos = 0; pos < consdata->nvars; ++pos )
         {
            /* if the clique number in the normal clique at position pos is greater than the last found cliquenumber the
             * partition is invalid */
	    if( consdata->cliquepartition[pos] > lastcliquenum )
            {
               consdata->cliquepartitioned = FALSE;
               break;
            }
	    else if( consdata->cliquepartition[pos] == lastcliquenum )
	       ++lastcliquenum; 
         }
      }
      /* we need to make sure that our clique numbers of our negated clique will be in increasing order without gaps */
      if( consdata->negcliquepartitioned )
      {
	 lastcliquenum = 0;

         for( pos = 0; pos < consdata->nvars; ++pos )
         {
            /* if the clique number in the negated clique at position pos is greater than the last found cliquenumber the
             * partition is invalid */
	    if( consdata->negcliquepartition[pos] > lastcliquenum )
            {
               consdata->negcliquepartitioned = FALSE;
               break;
            }
	    else if( consdata->negcliquepartition[pos] == lastcliquenum )
	       ++lastcliquenum; 
         }
      }

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
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Bool             normalclique,       /** should normal cliquepartition be created */
   SCIP_Bool             negatedclique       /** should negated cliquepartition be created */
   )
{
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || (consdata->cliquepartition != NULL && consdata->negcliquepartition != NULL));

   if( normalclique && !consdata->cliquepartitioned )
   {
      SCIP_CALL( SCIPcalcCliquePartition(scip, consdata->vars, consdata->nvars, consdata->cliquepartition, &consdata->ncliques) );
      consdata->cliquepartitioned = TRUE;
   }

   if( negatedclique && !consdata->negcliquepartitioned )
   {
      SCIP_CALL( SCIPcalcNegatedCliquePartition(scip, consdata->vars, consdata->nvars, consdata->negcliquepartition, &consdata->nnegcliques) );
      consdata->negcliquepartitioned = TRUE;
   }

   assert(!consdata->cliquepartitioned || consdata->ncliques <= consdata->nvars);
   assert(!consdata->negcliquepartitioned || consdata->nnegcliques <= consdata->nvars);

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
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->negcliquepartition, consdata->varssize, newsize) );
      }
      else
      {
         assert(consdata->eventdatas == NULL);
         assert(consdata->cliquepartition == NULL);
         assert(consdata->negcliquepartition == NULL);
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

         if( weights[v] > 0 )
         {
            (*consdata)->vars[k] = vars[v];
            (*consdata)->weights[k] = weights[v];
            ++k;
         }
      }
      assert(k >= 0);

      (*consdata)->nvars = k;
      if( k < nvars )
      {
         if( k > 0 )
         {
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(*consdata)->vars, nvars, k) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(*consdata)->weights, nvars, k) );
         }
         else
         {
            SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, nvars);
            SCIPfreeBlockMemoryArray(scip, &(*consdata)->weights, nvars);

            (*consdata)->vars = NULL;
            (*consdata)->weights = NULL;
            assert( (*consdata)->nvars == 0 );
         }
      }
   }
   else
   {
      (*consdata)->vars = NULL;
      (*consdata)->weights = NULL;
      (*consdata)->nvars = 0;
   }

   /* capacity has to be greater or equal to zero */
   assert( capacity >= 0 );

   (*consdata)->varssize = (*consdata)->nvars;
   (*consdata)->capacity = capacity;
   (*consdata)->eventdatas = NULL;
   (*consdata)->cliquepartition = NULL;
   (*consdata)->negcliquepartition = NULL;
   (*consdata)->row = NULL;
   (*consdata)->weightsum = 0;
   (*consdata)->onesweightsum = 0;
   (*consdata)->ncliques = 0;
   (*consdata)->nnegcliques = 0;
   (*consdata)->propagated = FALSE;
   (*consdata)->presolved = FALSE;
   (*consdata)->sorted = FALSE;
   (*consdata)->cliquepartitioned = FALSE;
   (*consdata)->negcliquepartitioned = FALSE;
   (*consdata)->merged = FALSE;
   (*consdata)->cliquesadded = FALSE;

   /* get transformed variables, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );

      /* allocate memory for additional data structures */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->eventdatas, (*consdata)->nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->cliquepartition, (*consdata)->nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->negcliquepartition, (*consdata)->nvars) );

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
   if( (*consdata)->negcliquepartition != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->negcliquepartition, (*consdata)->varssize);
   }
   if( (*consdata)->cliquepartition != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->cliquepartition, (*consdata)->varssize);
   }
   if( (*consdata)->vars != NULL )
   {
      assert( (*consdata)->weights != NULL );
      assert( (*consdata)->varssize > 0 );
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->varssize);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->weights, (*consdata)->varssize);
   }

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

   /* recalculate cliques extraction after a weight was increased */
   if( oldweight < newweight )
   {
      consdata->cliquesadded = FALSE;
   }
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
#define IDX(j,d) ((j)*(intcap)+(d))

/** solves knapsack problem in maximization form exactly using dynamic programming;
 *  if needed, one can provide arrays to store all selected items and all not selected items
 *
 * @note in case you provide the solitems or nonsolitems array you also have to provide the counter part as well
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
   SCIP_Real*            solval,             /**< pointer to store optimal solution value, or NULL */
   SCIP_Bool*            success             /**< pointer to store if an error occured during solving(normally a memory
                                              *   problem) */
   ) 
{
   SCIP_RETCODE retcode;
   SCIP_Real* tempsort;
   SCIP_Real* optvalues;
   int intcap;
   SCIP_Longint d;
   int j;
   SCIP_Longint weightsum;
   int* myitems;
   SCIP_Longint* myweights;
   SCIP_Real* myprofits;
   int nmyitems;
   SCIP_Longint gcd;
   SCIP_Longint minweight;
   SCIP_Longint currminweight;
   SCIP_Longint greedycap;
   SCIP_Longint greedysolweight;
   SCIP_Real greedysolvalue;
   SCIP_Bool eqweights;
   SCIP_Bool isoptimal;

   assert(weights != NULL);
   assert(profits != NULL);
   assert(capacity >= 0);
   assert(items != NULL);
   assert(nitems >= 0);
   assert(success != NULL);

   *success = TRUE;
   
#ifndef NDEBUG
   for( j = nitems - 1; j >= 0; --j )
      assert(weights[j] >= 0);
#endif

   SCIPdebugMessage("Solving knapsack exactly.\n");

   /* initializing solution value */
   if( solval != NULL )
      *solval = 0.0;

   /* produces optimal solution by following the table */
   if( solitems != NULL)
   {
      assert(items != NULL);
      assert(nsolitems != NULL);
      assert(nonsolitems != NULL);
      assert(nnonsolitems != NULL);

      *nnonsolitems = 0;
      *nsolitems = 0;
   }

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &myweights, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &myprofits, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &myitems, nitems) );
   nmyitems = 0;
   weightsum = 0;
   minweight = SCIP_LONGINT_MAX;

   /* remove unnecessary items */
   for( j = 0; j < nitems; ++j )
   {
      assert(0 <= weights[j] && weights[j] < SCIP_LONGINT_MAX);
      /* items does not fit */
      if( weights[j] > capacity )
      {
         if( solitems != NULL)
         {
            nonsolitems[*nnonsolitems] = items[j];
            ++(*nnonsolitems);
         }
      }
      /* items we does not want */
      else if( profits[j] <= 0.0 )
      {
         if( solitems != NULL)
         {
            nonsolitems[*nnonsolitems] = items[j];
            ++(*nnonsolitems);
         }
      }
      /* items which always fit */
      else if( weights[j] == 0 )
      {
         if( solitems != NULL)
         {
            solitems[*nsolitems] = items[j];
            ++(*nsolitems);
         }
         if( solval != NULL )
            *solval += profits[j];
      }
      /* all important items */
      else
      {
         myweights[nmyitems] = weights[j];
         myprofits[nmyitems] = profits[j];
         myitems[nmyitems] = items[j];

         /* remember smallest item */
         if( myweights[nmyitems] < minweight )
            minweight = myweights[nmyitems];

         weightsum += myweights[nmyitems];
         ++nmyitems;
      }
   }

   /* no item is left then goto end */
   if( nmyitems == 0 )
   {
      SCIPdebugMessage("After preprocessing no items are left.\n");

      goto TERMINATE;
   }
   /* if all items fit, we also do not need to do the expensive stuff later on */
   else if( weightsum > 0 && weightsum <= capacity )
   {
      SCIPdebugMessage("After preprocessing all items fit into knapsack.\n");
      
      for( j = nmyitems - 1; j >= 0; --j )
      {
         if( solitems != NULL )
         {
            solitems[*nsolitems] = myitems[j];
            ++(*nsolitems);
         }
         if( solval != NULL )
            *solval += myprofits[j];
      }
      
      goto TERMINATE;
   }

   /* determine greatest common divisor */
   gcd = myweights[nmyitems - 1];
   for( j = nmyitems - 2; j >= 0 && gcd >= 2; --j )
      gcd = SCIPcalcGreComDiv(gcd, myweights[j]);

   SCIPdebugMessage("Gcd is %"SCIP_LONGINT_FORMAT".\n", gcd);

   /* divide by greatest common divisor */
   if( gcd > 1 )
   {
      eqweights = TRUE;
      for( j = nmyitems - 1; j >= 0; --j )
      {
         myweights[j] /= gcd;
         eqweights &= (myweights[j] == 1);
      }
      capacity /= gcd;
      minweight /= gcd;
   }
   else
      eqweights = FALSE;

   assert(minweight <= capacity);

   /* only one item fits, than take the best */
   if( minweight > capacity / 2 )
   {
      int p;

      SCIPdebugMessage("Only one item fits into knapsack, so take the best.\n");

      p = nmyitems - 1;

      /* find best item */
      for( j = nmyitems - 2; j >= 0; --j )
         if( myprofits[j] > myprofits[p] )
            p = j;

      /* update solution information */
      if( solitems != NULL)
      {
         solitems[*nsolitems] = myitems[p];
         ++(*nsolitems);
         for( j = nmyitems - 1; j >= 0; --j )
            if( j != p )
            {
               nonsolitems[*nnonsolitems] = myitems[j];
               ++(*nnonsolitems);
            }
      }
      /* update solution value */
      if( solval != NULL )
         *solval += myprofits[p];
      
      goto TERMINATE;
   }

   /* all items have the same weight, than take the best */
   if( eqweights )
   {
      SCIP_Real addval;

      SCIPdebugMessage("All weights are equal, so take the best.\n");

      SCIPsortDownRealIntLong(myprofits, myitems, myweights, nmyitems);

      addval = 0.0;
      /* update solution information */
      if( solitems != NULL)
      {
         /* if all items would fit we had handled this case before */
         assert(nmyitems > capacity);

         /* take the first best items into the solution */
         for( d = capacity - 1; d >= 0; --d )
         {
            solitems[*nsolitems] = myitems[d];
            ++(*nsolitems);
            addval += myprofits[d];
         }

         /* the rest are not in the solution */
         for( j = nmyitems - 1; j >= capacity; --j )
         {
            nonsolitems[*nnonsolitems] = myitems[j];
            ++(*nnonsolitems);
         }
      }
      /* update solution value */
      if( solval != NULL )
      {
         assert(addval > 0.0);
         *solval += addval;
      }

      goto TERMINATE;
   }

   /* in the follwoing table we do not need the first minweight columns */
   capacity -= (minweight - 1);

   /* we can only handle integers */
   if( capacity >= INT_MAX )
   {
      SCIPdebugMessage("Capacity is to big, so we cannot handle it here.\n");

      *success = FALSE;
      goto TERMINATE;
   }
   assert(capacity < INT_MAX);

   intcap = (int)capacity;
   assert(intcap >= 0);

   /* this condition is only to check if the size of memory which will be allocated is still positiv( so no error occurs
    * ), which will be necessary after this if condition
    */ 
   if( (intcap) < 0 || (nmyitems) * (intcap) < 0 || (nmyitems) * (intcap) * ((int) sizeof(*optvalues)) < 0 )
   {
      SCIPdebugMessage("Too much memory will be consumed.\n");

      *success = FALSE;
      goto TERMINATE;
   }

   /* allocate temporary memory and check for memory exceeding */ 
   retcode = SCIPallocBufferArray(scip, &optvalues, (nmyitems)*(intcap));
   if( retcode == SCIP_NOMEMORY )
   {
      SCIPdebugMessage("Did not get enough memory.\n");

      *success = FALSE;
      goto TERMINATE;
   }
   else
   {
      SCIP_CALL( retcode );
   }

   /* sort myitems (plus corresponding arrays myweights and myprofits) such that 
    * p_1/w_1 >= p_2/w_2 >= ... >= p_n/w_n, this is only use for greedy solution
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &tempsort, nmyitems) );
   for( j = nmyitems - 1; j >= 0; --j )
      tempsort[j] = myprofits[j]/myweights[j];

   SCIPsortDownRealLongRealInt(tempsort, myweights, myprofits, myitems, nmyitems);

   /* initialize values for greedy solution information */
   greedysolweight = 0;
   greedysolvalue = 0.0;
   isoptimal = TRUE;
   greedycap = capacity + (minweight - 1);

   SCIPdebugMessage("Determine greedy solution.\n");

   /* determine greedy solution */
   for( j = 0; j < nmyitems; ++j )
   {
      assert(myweights[j] <= greedycap);

      /* take all fitting items */
      if( myweights[j] + greedysolweight <= greedycap )
      {
         /* update greedy solution weight and value */
         greedysolweight += myweights[j];
         greedysolvalue += myprofits[j];
         continue;
      }
      else if( greedysolweight < greedycap )
         isoptimal = FALSE;
      break;
   }
   assert(greedysolweight > 0);
   assert(greedysolvalue > 0.0);

   /* greedy solution is optimal */
   if( isoptimal )
   {
      assert(greedysolweight == greedycap);

      SCIPdebugMessage("Greedy solution is optimal.\n");

      greedysolweight = 0;

      /* update solution information */
      if( solitems != NULL)
      {
         /* take the first best items into the solution */
         for( j = 0; j < nmyitems; ++j )
         {
            /* take all fitting items */
            if( myweights[j] + greedysolweight <= greedycap )
            {
               solitems[*nsolitems] = myitems[j];
               ++(*nsolitems);
               greedysolweight += myweights[j];
            }
            else
            {
               nonsolitems[*nnonsolitems] = myitems[j];
               ++(*nnonsolitems);
            }
         }
      }
      /* update solution value */
      if( solval != NULL )
      {
         assert(greedysolvalue > 0.0);
         *solval += greedysolvalue;
      }

      SCIPfreeBufferArray(scip, &tempsort);
      SCIPfreeBufferArray(scip, &optvalues);

      goto TERMINATE;
   }

   SCIPdebugMessage("Start real exact algorithm.\n");

   currminweight = myweights[0] - minweight;
   /* fills first row of dynamic programming table with optimal values */
   for( d = currminweight; d < intcap; ++d )
      optvalues[d] = myprofits[0];
   /* fills dynamic programming table with optimal values */
   for( j = 1; j < nmyitems; ++j )
   {
      int intweight;

      /* compute important part of weight, which will be represented in the table */
      intweight = (int)(myweights[j] - minweight);
      assert(0 <= intweight && intweight < intcap);
      
      /* copy all nonzeros from row above */
      for( d = currminweight; d < intweight && d < intcap; ++d )
         optvalues[IDX(j,d)] = optvalues[IDX(j-1,d)];

      /* initialize last line if we found a new smallest weight which becomes necessary for later use */
      if( intweight < currminweight )
      {
         for( d = intweight; d < currminweight; ++d )
            optvalues[IDX(j-1,d)] = 0.0;
      }

      /* update corresponding row */
      for( d = intweight; d < intcap; ++d )
      {
         SCIP_Real sumprofit;

         if( d - myweights[j] < currminweight )
            sumprofit = myprofits[j];
         else
            sumprofit = optvalues[IDX(j-1,d-myweights[j])] + myprofits[j];
         optvalues[IDX(j,d)] = MAX(sumprofit, optvalues[IDX(j-1,d)]);
      }
      /* update currminweight */
      if( intweight < currminweight )
         currminweight = intweight;
   }

   /* update optimal solution by following the table */
   if( solitems != NULL)
   {
      d = intcap - 1;

      SCIPdebugMessage("Fill the solution vector after solving exactly.\n");

      /* insert all items in (non-) solution vector */
      for( j = nmyitems - 1; j > 0; --j )
      {
         /* if we cannot find any item which is in our solution stop */
         if( d < 0 || optvalues[IDX(j,d)] == 0.0 )
         {
            assert(d >= -minweight);
            break;
         }

         /* collect solution items */
         if( optvalues[IDX(j,d)] > optvalues[IDX(j-1,d)] )
         {
            solitems[*nsolitems] = myitems[j];
            ++(*nsolitems);
            d -= myweights[j];
         } 
         /* collect non-solution items */
         else
         { 
            nonsolitems[*nnonsolitems] = myitems[j];
            ++(*nnonsolitems);
         }
      }

      /* insert remaining items */
      if( d >= 0 && optvalues[IDX(j,d)] > 0.0 )
      {
         assert(j == 0);
         solitems[*nsolitems] = myitems[j];
         ++(*nsolitems);
      }
      else
      {
         assert(j >= 0);
         assert(d < 0 || optvalues[IDX(j,d)] == 0.0);
         for( ; j >= 0; --j )
         {
            nonsolitems[*nnonsolitems] = myitems[j];
            ++(*nnonsolitems);
         }
      }
         
      assert(*nsolitems + *nnonsolitems == nitems);
   }

   /* update solution value */
   if( solval != NULL )
      *solval += optvalues[IDX(nmyitems-1,intcap-1)];

   /* free all temporary memory */
   SCIPfreeBufferArray(scip, &tempsort);
   SCIPfreeBufferArray(scip, &optvalues);

 TERMINATE:
   SCIPfreeBufferArray(scip, &myitems);
   SCIPfreeBufferArray(scip, &myprofits);
   SCIPfreeBufferArray(scip, &myweights);

   return SCIP_OKAY;
}


/** solves knapsack problem in maximization form approximately by solving the LP-relaxation of the problem using Dantzig's
 *  method and rounding down the solution; if needed, one can provide arrays to store all selected items and all not 
 *  selected items
 */
SCIP_RETCODE SCIPsolveKnapsackApproximately(
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
   SCIP_Longint solitemsweight;
   int j;
   
   assert(weights != NULL);
   assert(profits != NULL);
   assert(capacity >= 0);
   assert(items != NULL);
   assert(nitems >= 0);

   if( solitems != NULL )
   {
      *nsolitems = 0;
      *nnonsolitems = 0;
   }
   if( solval != NULL )
      *solval = 0.0;

   /* sort items (plus corresponding arrays weights and profits) such that 
    * p_1/w_1 >= p_2/w_2 >= ... >= p_n/w_n 
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &tempsort, nitems) );
   for( j = nitems - 1; j >= 0; --j )
   {
      tempsort[j] = profits[j]/weights[j];
   }
   SCIPsortDownRealLongRealInt(tempsort, weights, profits, items, nitems);

   /* selects items as long as they fit into the knapsack */
   solitemsweight = 0;
   for( j = 0; j < nitems && solitemsweight + weights[j] <= capacity; j++ )
   {
      if( solitems != NULL )
      {
         solitems[*nsolitems] = items[j];
         (*nsolitems)++;
      }
      if( solval != NULL )
         (*solval) += profits[j];
      solitemsweight += weights[j];
   }
   for( ; j < nitems && solitems != NULL; j++ )
   {
      nonsolitems[*nnonsolitems] = items[j];
      (*nnonsolitems)++;
   }
   
   SCIPfreeBufferArray(scip, &tempsort);

   return SCIP_OKAY;
}

#ifndef NDEBUG
/** returns, whether the the arrays transweights, transprofits and items are sorted such that 
 *  p_1 / w_1 >= p_2 / w_2 >= ... >= p_n / w_n and these arrays are not changed
 */
static
SCIP_Bool checkSolveKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nitems,             /**< number of available items */
   SCIP_Longint*         transweights,       /**< item weights */
   SCIP_Real*            transprofits,       /**< item profits */
   int*                  items,              /**< item numbers */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   SCIP_Bool             modtransused        /**< TRUE for mod trans sepa prob was used to find cover */
   )
{
   int j;
   
   assert(scip != NULL);
   assert(nitems >= 0);
   assert(weights != NULL);
   assert(solvals != NULL);
   assert(transweights != NULL);
   assert(transprofits != NULL);
   
   for( j = 1; j < nitems; j++ )
   {
      assert(SCIPisFeasGE(scip, transprofits[j-1]/transweights[j-1], transprofits[j]/transweights[j]));
      if( SCIPisFeasLT(scip, transprofits[j-1]/transweights[j-1], transprofits[j]/transweights[j]) )
         return FALSE;
      assert(weights[items[j]] == transweights[j]);
      if( weights[items[j]] != transweights[j] )
         return FALSE;
      if( modtransused )
      {
         assert(SCIPisFeasEQ(scip, (1.0 - solvals[items[j]]) * weights[items[j]], transprofits[j]));
         if( !SCIPisFeasEQ(scip, (1.0 - solvals[items[j]]) * weights[items[j]], transprofits[j]) )
            return FALSE;
      }
      else
      {
         assert(SCIPisFeasEQ(scip, (1.0 - solvals[items[j]]), transprofits[j]));
         if( !SCIPisFeasEQ(scip, (1.0 - solvals[items[j]]), transprofits[j]) )
            return FALSE;

      }
   }
   return TRUE;
}
#endif

/** gets a most violated cover C (\f$\sum_{j \in C} a_j > a_0\f$) for a given knapsack constraint \f$\sum_{j \in N} a_j x_j \leq a_0\f$
 *  taking into consideration the following fixing: \f$j \in C\f$, if \f$j \in N_1 = \{j \in N : x^*_j = 1\}\f$ and 
 *  \f$j \in N \setminus C\f$, if \f$j \in N_0 = \{j \in N : x^*_j = 0\}\f$, if one exists.
 */
static
SCIP_RETCODE getCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  covervars,          /**< pointer to store cover variables */
   int*                  noncovervars,       /**< pointer to store noncover variables */
   int*                  ncovervars,         /**< pointer to store number of cover variables */
   int*                  nnoncovervars,      /**< pointer to store number of noncover variables */
   SCIP_Longint*         coverweight,        /**< pointer to store weight of cover */
   SCIP_Bool*            found,              /**< pointer to store whether a cover was found */
   SCIP_Bool             modtransused,       /**< should modified transformed separation problem be used to find cover */
   int*                  ntightened,         /**< pointer to store number of variables with tightened upper bound */
   SCIP_Bool*            fractional          /**< pointer to store whether the LP sol for knapsack vars is fractional */
   )
{
   SCIP_Longint* transweights;
   SCIP_Longint* modtransweights;
   SCIP_Real* transprofits;
   SCIP_Real* modtransprofits;
   SCIP_Longint transcapacity;
   SCIP_Longint modtranscapacity;
   SCIP_Longint fixedonesweight;
   SCIP_Longint itemsweight;
   SCIP_Bool infeasible;
   int* fixedones;
   int* fixedzeros;
   int* items;
   int nfixedones;
   int nfixedzeros;
   int nitems;
   int j; 

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(weights != NULL);
   assert(capacity >= 0);
   assert(solvals != NULL);
   assert(covervars != NULL);
   assert(noncovervars != NULL);
   assert(ncovervars != NULL);
   assert(nnoncovervars != NULL);
   assert(coverweight != NULL);
   assert(found != NULL);
   assert(ntightened != NULL);
   assert(fractional != NULL);

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &transweights, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transprofits, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedones, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedzeros, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &items, nvars) );
   if( modtransused )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &modtransweights, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &modtransprofits, nvars) );
   }
   else
   {
      modtransweights = NULL;
      modtransprofits = NULL;
   }

   *found = FALSE;
   *ncovervars = 0;
   *nnoncovervars = 0;
   *coverweight = 0;
   *fractional = TRUE;
      
   /* gets the following sets 
    *  N_1 = {j in N : x*_j = 1} (fixedones), 
    *  N_0 = {j in N : x*_j = 0} (fixedzeros) and 
    *  N\(N_0 & N_1) (items), 
    * where x*_j is the solution value of variable x_j
    */
   nfixedones = 0;
   nfixedzeros = 0;
   nitems = 0;
   fixedonesweight = 0;
   itemsweight = 0;
   *ntightened = 0;
   for( j = 0; j < nvars; j++ )
   {
      assert(SCIPvarIsBinary(vars[j]));
      
      /* tightens upper bound of x_j if weight of x_j is greater than capacity of knapsack */
      if( weights[j] > capacity ) 
      {
         SCIP_CALL( SCIPtightenVarUb(scip, vars[j], 0.0, FALSE, &infeasible, NULL) );
         assert(!infeasible);
         (*ntightened)++;
         continue;
      }

      /* variable x_j has solution value one */
      if( SCIPisFeasEQ(scip, solvals[j], 1.0) ) 
      {
         fixedones[nfixedones] = j;
         nfixedones++;
         fixedonesweight += weights[j];
      }
      /* variable x_j has solution value zero */
      else if( SCIPisFeasEQ(scip, solvals[j], 0.0) )
      {
         fixedzeros[nfixedzeros] = j;
         nfixedzeros++;
      }
      /* variable x_j has fractional solution value */
      else
      {
         assert( SCIPisFeasGT(scip, solvals[j], 0.0) && SCIPisFeasLT(scip, solvals[j], 1.0) );
         items[nitems] = j;
         nitems++;
         itemsweight += weights[j];
      }
   }
   assert(nfixedones + nfixedzeros + nitems == nvars - (*ntightened));
   
   /* sets whether the LP solution x* for the knapsack variables is fractional; if it is not fractional we stop 
    * the separation routine 
    */
   assert(nitems >= 0);
   if( nitems == 0 )
   {
      *fractional = FALSE;
      goto TERMINATE;
   }
   assert(*fractional);

   /* transforms the traditional separation problem (under consideration of the following fixing: 
    * z_j = 1 for all j in N_1, z_j = 0 for all j in N_0)
    *
    *   min sum_{j in N\(N_0 & N_1)} (1 - x*_j) z_j
    *       sum_{j in N\(N_0 & N_1)} a_j z_j >= (a_0 + 1) - sum_{j in N_1} a_j
    *                                    z_j in {0,1}, j in N\(N_0 & N_1)
    *
    * to a knapsack problem in maximization form by complementing the variables
    *  
    * sum_{j in N\(N_0 & N_1)} (1 - x*_j) - 
    *   max sum_{j in N\(N_0 & N_1)} (1 - x*_j) z°_j
    *       sum_{j in N\(N_0 & N_1)} a_j z°_j <= sum_{j in N\N_0} a_j - (a_0 + 1) 
    *                                    z°_j in {0,1}, j in N\(N_0 & N_1)
    */

   /* gets weight and profit of variables in transformed knapsack problem */
   for( j = 0; j < nitems; j++ )
   {
      transweights[j] = weights[items[j]];
      transprofits[j] = 1.0 - solvals[items[j]];
   }
   /* gets capacity of transformed knapsack problem */
   transcapacity = fixedonesweight + itemsweight - capacity - 1;

   /* if capacity of transformed knapsack problem is less than zero, there is no cover 
    * (when variables fixed to zero are not used)
    */
   if( transcapacity < 0 )
   {
      assert(!(*found));
      goto TERMINATE;
   }
   
   if( modtransused )
   {
      assert(modtransweights != NULL);
      assert(modtransprofits != NULL);
      /* transforms the modified separation problem (under consideration of the following fixing: 
       * z_j = 1 for all j in N_1, z_j = 0 for all j in N_0)
       *
       *   min sum_{j in N\(N_0 & N_1)} (1 - x*_j) a_j z_j
       *       sum_{j in N\(N_0 & N_1)} a_j z_j >= (a_0 + 1) - sum_{j in N_1} a_j
       *                                    z_j in {0,1}, j in N\(N_0 & N_1)
       *
       * to a knapsack problem in maximization form by complementing the variables
       *  
       * sum_{j in N\(N_0 & N_1)} (1 - x*_j) a_j - 
       *   max sum_{j in N\(N_0 & N_1)} (1 - x*_j) a_j z°_j
       *       sum_{j in N\(N_0 & N_1)} a_j z°_j <= sum_{j in N\N_0} a_j - (a_0 + 1) 
       *                                    z°_j in {0,1}, j in N\(N_0 & N_1)
       */

      /* gets weight and profit of variables in modified transformed knapsack problem */
      for( j = 0; j < nitems; j++ )
      {
         modtransweights[j] = weights[items[j]];
         modtransprofits[j] = (1.0 - solvals[items[j]]) * weights[items[j]];
         assert(SCIPisFeasGT(scip, modtransprofits[j], 0.0));
      }
      /* gets capacity of modified transformed knapsack problem */
      modtranscapacity = fixedonesweight + itemsweight - capacity - 1;

      /* solves modified transformed knapsack problem approximately by solving the LP-relaxation of the 
       * modified transformed knapsack problem using Dantzig's method and rounding down the solution.
       * let z°* be the solution, then 
       *   j in C,          if z°*_j = 0 and
       *   i in N\C,        if z°*_j = 1.
       */ 
      SCIP_CALL(SCIPsolveKnapsackApproximately(scip, nitems, modtransweights, modtransprofits, modtranscapacity, items, 
            noncovervars, covervars, nnoncovervars, ncovervars, NULL)); 
      assert(checkSolveKnapsack(scip, nitems, modtransweights, modtransprofits, items, weights, solvals, modtransused));
   }
   else
   {
      /* solves transformed knapsack problem approximately by solving the LP-relaxation of the 
       * transformed knapsack problem using Dantzig's method and rounding down the solution.
       * let z°* be the solution, then 
       *   j in C,          if z°*_j = 0 and
       *   i in N\C,        if z°*_j = 1.
       */ 
      SCIP_CALL(SCIPsolveKnapsackApproximately(scip, nitems, transweights, transprofits, transcapacity, items, 
            noncovervars, covervars, nnoncovervars, ncovervars, NULL));
      assert(checkSolveKnapsack(scip, nitems, transweights, transprofits, items, weights, solvals, modtransused));
   }

   /* constructs cover C (sum_{j in C} a_j > a_0) */
   for( j = 0; j < *ncovervars; j++ )
   {
      (*coverweight) += weights[covervars[j]];
   }
   
   /* adds all variables from N_1 to C */
   for( j = 0; j < nfixedones; j++ )
   {
      covervars[*ncovervars] = fixedones[j];
      (*ncovervars)++;
      (*coverweight) += weights[fixedones[j]];
   }

   /* adds all variables from N_0 to N\C */
   for( j = 0; j < nfixedzeros; j++ )
   {
      noncovervars[*nnoncovervars] = fixedzeros[j];
      (*nnoncovervars)++;
   }
   assert((*ncovervars) + (*nnoncovervars) == nvars - (*ntightened));
   assert((*coverweight) > capacity);
   *found = TRUE;
   
 TERMINATE:
   /* frees temporary memory */
   if( modtransused )
   {
      SCIPfreeBufferArray(scip, &modtransprofits);
      SCIPfreeBufferArray(scip, &modtransweights);
   }
   SCIPfreeBufferArray(scip, &items);
   SCIPfreeBufferArray(scip, &fixedzeros);
   SCIPfreeBufferArray(scip, &fixedones);
   SCIPfreeBufferArray(scip, &transprofits);
   SCIPfreeBufferArray(scip, &transweights);

   return SCIP_OKAY;
}

#ifndef NDEBUG
/** checks if minweightidx is set correctly
 */
static
SCIP_Bool checkMinweightidx(
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   int*                  covervars,          /**< pointer to store cover variables */
   int                   ncovervars,         /**< pointer to store number of cover variables */
   SCIP_Longint          coverweight,        /**< pointer to store weight of cover */
   int                   minweightidx,       /**< index of variable in cover variables with minimum weight */
   int                   j                   /**< current index in cover variables */      
  )
{
   SCIP_Longint minweight;
   int i;

   assert(weights != NULL);
   assert(covervars != NULL);
   assert(ncovervars > 0);

   minweight = weights[covervars[minweightidx]];

   /* checks if all cover variables before index j have weight greater than minweight */
   for( i = 0; i < j; i++ )
   {
      assert(weights[covervars[i]] > minweight);
      if( weights[covervars[i]] <= minweight )
         return FALSE;
   }

   /* checks if all variables before index j can not be removed, i.e. i cannot be the next minweightidx */
   for( i = 0; i < j; i++ )
   {
      assert(coverweight - weights[covervars[i]] <= capacity);
      if( coverweight - weights[covervars[i]] > capacity )
         return FALSE;
   }
   return TRUE;
}
#endif


/** gets partition \f$(C_1,C_2)\f$ of minimal cover \f$C\f$, i.e. \f$C_1 \cup C_2 = C\f$ and \f$C_1 \cap C_2 = \emptyset\f$, 
 *  with \f$C_1\f$ not empty; choses partition as follows \f$C_2 = \{ j \in C : x^*_j = 1 \}\f$ and \f$C_1 = C \setminus C_2\f$
 */
static
void getPartitionCovervars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  covervars,          /**< cover variables */
   int                   ncovervars,         /**< number of cover variables */
   int*                  varsC1,             /**< pointer to store variables in C1 */
   int*                  varsC2,             /**< pointer to store variables in C2 */
   int*                  nvarsC1,            /**< pointer to store number of variables in C1 */
   int*                  nvarsC2             /**< pointer to store number of variables in C2 */
   )
{
   int j;

   assert(scip != NULL);
   assert(ncovervars >= 0);
   assert(solvals != NULL);
   assert(covervars != NULL);
   assert(varsC1 != NULL);
   assert(varsC2 != NULL);
   assert(nvarsC1 != NULL);
   assert(nvarsC2 != NULL);
   
   *nvarsC1 = 0;
   *nvarsC2 = 0;
   for( j = 0; j < ncovervars; j++ )
   {
      assert(SCIPisFeasGT(scip, solvals[covervars[j]], 0.0));

      /* variable has solution value one */
      if( SCIPisGE(scip, solvals[covervars[j]], 1.0) )
      {
         varsC2[*nvarsC2] = covervars[j];
         (*nvarsC2)++;
      }
      /* variable has solution value less than one */
      else
      {
         assert(SCIPisLT(scip, solvals[covervars[j]], 1.0));
         varsC1[*nvarsC1] = covervars[j];
         (*nvarsC1)++;
      }
   }
   assert((*nvarsC1) + (*nvarsC2) == ncovervars);
}

/** changes given partition (C_1,C_2) of minimal cover C, if |C1| = 1, by moving one and two (if possoble) variables from 
 *  C2 to C1 if |C1| = 1 and |C1| = 0, respectively.
 */
static
SCIP_RETCODE changePartitionCovervars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   int*                  varsC1,             /**< pointer to store variables in C1 */
   int*                  varsC2,             /**< pointer to store variables in C2 */
   int*                  nvarsC1,            /**< pointer to store number of variables in C1 */
   int*                  nvarsC2             /**< pointer to store number of variables in C2 */
   )
{
   SCIP_Real* sortkeysC2;
   int j;

   assert(*nvarsC1 >= 0 && *nvarsC1 <= 1);
   assert(*nvarsC2 > 0);

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeysC2, *nvarsC2) );

   /* sorts variables in C2 such that a_1 >= .... >= a_|C2| */
   for( j = 0; j < *nvarsC2; j++ )
      sortkeysC2[j] = (SCIP_Real) weights[varsC2[j]];
   SCIPsortDownRealInt(sortkeysC2, varsC2, *nvarsC2);

   /* adds one or two variable from C2 with smallest weight to C1 and removes them from C2 */
   assert(*nvarsC2 == 1 || weights[varsC2[(*nvarsC2)-1]] <= weights[varsC2[(*nvarsC2)-2]]);
   while( *nvarsC1 < 2 && *nvarsC2 > 0 )
   {
      varsC1[*nvarsC1] = varsC2[(*nvarsC2)-1];
      (*nvarsC1)++;
      (*nvarsC2)--;
   }
   
   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &sortkeysC2);

   return SCIP_OKAY;
}

/** changes given partition (C_1,C_2) of feasible set C, if |C1| = 1, by moving one variable from C2 to C1 */
static
SCIP_RETCODE changePartitionFeasiblesetvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   int*                  varsC1,             /**< pointer to store variables in C1 */
   int*                  varsC2,             /**< pointer to store variables in C2 */
   int*                  nvarsC1,            /**< pointer to store number of variables in C1 */
   int*                  nvarsC2             /**< pointer to store number of variables in C2 */
   )
{
   SCIP_Real* sortkeysC2;
   int j;

   assert(*nvarsC1 >= 0 && *nvarsC1 <= 1);
   assert(*nvarsC2 > 0);

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeysC2, *nvarsC2) );

   /* sorts variables in C2 such that a_1 >= .... >= a_|C2| */
   for( j = 0; j < *nvarsC2; j++ )
      sortkeysC2[j] = (SCIP_Real) weights[varsC2[j]];
   SCIPsortDownRealInt(sortkeysC2, varsC2, *nvarsC2);

   /* adds variable from C2 with smallest weight to C1 and removes it from C2 */
   assert(*nvarsC2 == 1 || weights[varsC2[(*nvarsC2)-1]] <= weights[varsC2[(*nvarsC2)-2]]);
   varsC1[*nvarsC1] = varsC2[(*nvarsC2)-1];
   (*nvarsC1)++;
   (*nvarsC2)--;
   
   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &sortkeysC2);

   return SCIP_OKAY;
}


/** gets partition \f$(F,R)\f$ of \f$N \setminus C\f$ where \f$C\f$ is a minimal cover, i.e. \f$F \cup R = N \setminus C\f$
 *  and \f$F \cap R = \emptyset\f$; choses partition as follows \f$R = \{ j \in N \setminus C : x^*_j = 0 \}\f$ and 
 *  \f$F = (N \setminus C) \setminus F\f$
 */
static
void getPartitionNoncovervars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  noncovervars,       /**< noncover variables */
   int                   nnoncovervars,      /**< number of noncover variables */
   int*                  varsF,              /**< pointer to store variables in F */
   int*                  varsR,              /**< pointer to store variables in R */
   int*                  nvarsF,             /**< pointer to store number of variables in F */
   int*                  nvarsR              /**< pointer to store number of variables in R */
   )
{
   int j;

   assert(scip != NULL);
   assert(nnoncovervars >= 0);
   assert(solvals != NULL);
   assert(noncovervars != NULL);
   assert(varsF != NULL);
   assert(varsR != NULL);
   assert(nvarsF != NULL);
   assert(nvarsR != NULL);
   
   *nvarsF = 0;
   *nvarsR = 0;
   for( j = 0; j < nnoncovervars; j++ )
   {
      /* variable has solution value zero */
      if( SCIPisFeasEQ(scip, solvals[noncovervars[j]], 0.0) )
      {
         varsR[*nvarsR] = noncovervars[j];
         (*nvarsR)++;
      }
      /* variable has solution value greater than zero */
      else
      {
         assert(SCIPisFeasGT(scip, solvals[noncovervars[j]], 0.0));
         varsF[*nvarsF] = noncovervars[j];
         (*nvarsF)++;
      }
   }
   assert((*nvarsF) + (*nvarsR) == nnoncovervars);
}

/** sorts variables in F, C_2, and R according to the second level lifting sequence that will be used in the sequential 
 *  lifting procedure 
 */
static 
SCIP_RETCODE getLiftingSequence(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   int*                  varsF,              /**< pointer to store variables in F */
   int*                  varsC2,             /**< pointer to store variables in C2 */
   int*                  varsR,              /**< pointer to store variables in R */
   int                   nvarsF,             /**< number of variables in F */
   int                   nvarsC2,            /**< number of variables in C2 */
   int                   nvarsR              /**< number of variables in R */
   )
{
   SORTKEYPAIR** sortkeypairsF;
   SCIP_Real* sortkeysC2;
   SCIP_Real* sortkeysR;
   int j;

   assert(scip != NULL);
   assert(solvals != NULL);
   assert(weights != NULL);
   assert(varsF != NULL);
   assert(varsC2 != NULL);
   assert(varsR != NULL);
   assert(nvarsF >= 0);
   assert(nvarsC2 >= 0);
   assert(nvarsR >= 0);
   
   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sortkeypairsF, nvarsF) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeysC2, nvarsC2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeysR, nvarsR) );
   
   /* gets sorting key for variables in F corresponding to the following lifting sequence 
    *  sequence 1: non-increasing absolute difference between x*_j and the value the variable is fixed to, i.e. 
    *              x*_1 >= x*_2 >= ... >= x*_|F| 
    * in case of equality uses 
    *  sequence 4: non-increasing a_j, i.e. a_1 >= a_2 >= ... >= a_|C_2| 
    */
   for( j = 0; j < nvarsF; j++ )
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &sortkeypairsF[j]) );

      sortkeypairsF[j]->key1 = solvals[varsF[j]]; 
      sortkeypairsF[j]->key2 = (SCIP_Real) weights[varsF[j]]; 
   }
   
   /* gets sorting key for variables in C_2 corresponding to the following lifting sequence 
    *  sequence 4: non-increasing a_j, i.e. a_1 >= a_2 >= ... >= a_|C_2| 
    */
   for( j = 0; j < nvarsC2; j++ )
      sortkeysC2[j] = (SCIP_Real) weights[varsC2[j]]; 
   
   /* gets sorting key for variables in R corresponding to the following lifting sequence 
    *  sequence 4: non-increasing a_j, i.e. a_1 >= a_2 >= ... >= a_|R| 
    */
   for( j = 0; j < nvarsR; j++ )
      sortkeysR[j] = (SCIP_Real) weights[varsR[j]]; 
   
   /* sorts F, C2 and R */
   if( nvarsF > 0 )
   {
      SCIPsortDownPtrInt((void**)sortkeypairsF, varsF, compSortkeypairs, nvarsF);
   }   
   if( nvarsC2 > 0 )
   {
      SCIPsortDownRealInt(sortkeysC2, varsC2, nvarsC2);
   }   
   if( nvarsR > 0)
   {
      SCIPsortDownRealInt(sortkeysR, varsR, nvarsR);
   }
   
   /* frees temporary memory */
   for( j = nvarsF-1; j >= 0; j-- )
      SCIPfreeBlockMemory(scip, &sortkeypairsF[j]);
   SCIPfreeBufferArray(scip, &sortkeysR);
   SCIPfreeBufferArray(scip, &sortkeysC2);
   SCIPfreeBlockMemoryArray(scip, &sortkeypairsF, nvarsF);
   
   return SCIP_OKAY;
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
   assert(minweightslen != NULL);
   assert(*minweightslen >= 0);
   assert(minweightssize != NULL);
   assert(*minweightssize >= 0);

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

/** lifts given inequality 
 *    sum_{j in M_1} x_j <= alpha_0 
 *  valid for 
 *    S^0 = { x in {0,1}^|M_1| : sum_{j in M_1} a_j x_j <= a_0 - sum_{j in M_2} a_j } 
 *  to a valid inequality 
 *    sum_{j in M_1} x_j + sum_{j in F} alpha_j x_j + sum_{j in M_2} alpha_j x_j + sum_{j in R} alpha_j x_j 
 *    <= alpha_0 + sum_{j in M_2} alpha_j 
 *  for 
 *    S = { x in {0,1}^|N| : sum_{j in N} a_j x_j <= a_0 }; 
 *  uses sequential up-lifting for the variables in F, sequential down-lifting for the variable in M_2, and 
 *  sequential up-lifting for the variabels in R; procedure can be used to strengthen minimal cover inequalities and 
 *  extended weight inequalities.   
 */ 
static
SCIP_RETCODE sequentialUpAndDownLifting(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   int                   ntightened,         /**< number of variables with tightened upper bound */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  varsM1,             /**< variables in M_1 */
   int*                  varsM2,             /**< variables in M_2 */
   int*                  varsF,              /**< variables in F */
   int*                  varsR,              /**< variables in R */
   int                   nvarsM1,            /**< number of variables in M_1 */
   int                   nvarsM2,            /**< number of variables in M_2 */
   int                   nvarsF,             /**< number of variables in F */
   int                   nvarsR,             /**< number of variables in R */
   int                   alpha0,             /**< rights hand side of given valid inequality */
   int*                  liftcoefs,          /**< pointer to store lifting coefficient of vars in knapsack constraint */
   SCIP_Real*            cutact,             /**< pointer to store activity of lifted valid inequality */
   int*                  liftrhs             /**< pointer to store right hand side of the lifted valid inequality */
)
{
   SCIP_Longint* minweights;
   SCIP_Real* sortkeys;
   SCIP_Longint fixedonesweight;
   int minweightssize;
   int minweightslen;
   int j;
   int w;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nvars >= 0);
   assert(weights != NULL);
   assert(capacity >= 0);
   assert(solvals != NULL);
   assert(varsM1 != NULL);
   assert(varsM2 != NULL);
   assert(varsF != NULL);
   assert(varsR != NULL);
   assert(nvarsM1 >= 0 && nvarsM1 <= nvars - ntightened);
   assert(nvarsM2 >= 0 && nvarsM2 <= nvars - ntightened);
   assert(nvarsF >= 0 && nvarsF <= nvars - ntightened);
   assert(nvarsR >= 0 && nvarsR <= nvars - ntightened);
   assert(nvarsM1 + nvarsM2 + nvarsF + nvarsR == nvars  - ntightened);
   assert(alpha0 >= 0);
   assert(liftcoefs != NULL);
   assert(cutact != NULL);
   assert(liftrhs != NULL);

   /* allocates temporary memory */
   minweightssize = nvarsM1 + 1;
   SCIP_CALL( SCIPallocBufferArray(scip, &minweights, minweightssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeys, nvarsM1) );
   
   /* initializes data structures */
   BMSclearMemoryArray(liftcoefs, nvars);
   *cutact = 0.0;

   /* sets lifting coefficient of variables in M1, sorts variables in M1 such that a_1 <= a_2 <= ... <= a_|M1| 
    * and calculates activity of the current valid inequality 
    */
   for( j = 0; j < nvarsM1; j++ )
   {
      assert(liftcoefs[varsM1[j]] == 0);
      liftcoefs[varsM1[j]] = 1;
      sortkeys[j] = (SCIP_Real) (weights[varsM1[j]]);
      (*cutact) += solvals[varsM1[j]];
   }

   SCIPsortRealInt(sortkeys, varsM1, nvarsM1);
  
   /* initializes (i = 1) the minweight table, defined as: minweights_i[w] = 
    *   min   sum_{j in M_1} a_j x_j + sum_{k=1}^{i-1} a_{j_k}     x_{j_k}   
    *   s.t.  sum_{j in M_1}     x_j + sum_{k=1}^{i-1} alpha_{j_k} x_{j_k} >= w
    *                                    x_j in {0,1} for j in M_1 & {j_i,...,j_i-1},
    * for i = 1,...,t with t = |N\M1| and w = 0,...,|M1| + sum_{k=1}^{i-1} alpha_{j_k}; 
    */
   minweights[0] = 0;
   for( w = 1; w <= nvarsM1; w++ )
      minweights[w] = minweights[w-1] + weights[varsM1[w-1]];
   minweightslen = nvarsM1 + 1;

   /* gets sum of weights of variables fixed to one, i.e. sum of weights of variables in M_2 */
   fixedonesweight = 0;   
   for( j = 0; j < nvarsM2; j++ )
      fixedonesweight += weights[varsM2[j]];
   assert(fixedonesweight >= 0);

   /* initializes right hand side of lifted valid inequality */
   *liftrhs = alpha0;
   
   /* sequentially up-lifts all variables in F: */
   for( j = 0; j < nvarsF; j++ )
   {
      SCIP_Longint weight;
      int liftvar;
      int liftcoef;
      int z; 

      liftvar = varsF[j];
      weight = weights[liftvar];
      assert(liftvar >= 0 && liftvar < nvars);
      assert(SCIPisFeasGT(scip, solvals[liftvar], 0.0));
      assert(weight > 0);

      /* knapsack problem is infeasible: 
       *   sets z = 0 
       */
      if( capacity - fixedonesweight - weight < 0 )
      {
         z = 0; 
      }
      /* knapsack problem is feasible: 
       *   sets z = max { w : 0 <= w <= liftrhs, minweights_i[w] <= a_0 - fixedonesweight - a_{j_i}  } = liftrhs,
       *   if minweights_i[liftrhs] <= a_0 - fixedonesweight - a_{j_i}
       */
      else if( minweights[*liftrhs] <= capacity - fixedonesweight - weight )
      {
         z = *liftrhs;
      } 
      /* knapsack problem is feasible: 
       *   uses binary search to find z = max { w : 0 <= w <= liftrhs, minweights_i[w] <= a_0 - fixedonesweight - a_{j_i} }
       */
      else
      {
         int left;
         int right;
         int middle;
         
         assert((*liftrhs) + 1 >= minweightslen || minweights[(*liftrhs) + 1] > capacity - fixedonesweight - weight);
         left = 0;
         right = (*liftrhs) + 1;
         while( left < right - 1 )
         {
            middle = (left + right) / 2;
            assert(0 <= middle && middle < minweightslen);
            if( minweights[middle] <= capacity - fixedonesweight - weight ) 
               left = middle;
            else
               right = middle;
         }
         assert(left == right - 1);
         assert(0 <= left && left < minweightslen);
         assert(minweights[left] <= capacity - fixedonesweight - weight );
         assert(left == minweightslen - 1 || minweights[left+1] > capacity - fixedonesweight - weight);

         /* now z = left */
         z = left;
         assert(z <= *liftrhs);
      }
         
      /* calculates lifting coefficients alpha_{j_i} = liftrhs - z */
      liftcoef = (*liftrhs) - z;
      liftcoefs[liftvar] = liftcoef;
      assert(liftcoef >= 0 && liftcoef <= (*liftrhs) + 1);

      /* minweight table and activity of current valid inequality will not change, if alpha_{j_i} = 0 */
      if( liftcoef == 0 )
         continue;

      /* updates activity of current valid inequality */
      (*cutact) += liftcoef * solvals[liftvar];

      /* enlarges current minweight table:
       *  from minweightlen = |M1| + sum_{k=1}^{i-1} alpha_{j_k} + 1 entries 
       *  to                  |M1| + sum_{k=1}^{i  } alpha_{j_k} + 1 entries
       * and sets minweights_i[w] = infinity for 
       *  w = |M1| + sum_{k=1}^{i-1} alpha_{j_k} + 1 , ... , |M1| + sum_{k=1}^{i} alpha_{j_k} 
       */
      SCIP_CALL( enlargeMinweights(scip, &minweights, &minweightslen, &minweightssize, minweightslen + liftcoef) );
     
      /* updates minweight table: minweight_i+1[w] = 
       *   min{ minweights_i[w], a_{j_i}},                                 if w <  alpha_j_i
       *   min{ minweights_i[w], minweights_i[w - alpha_j_i] + a_j_i},     if w >= alpha_j_i
       */
      for( w = minweightslen - 1; w >= 0; w-- )
      {
         SCIP_Longint min;
         if( w < liftcoef )
         {
            min = MIN(minweights[w], weight);
            minweights[w] = min;
         }
         else
         {
            assert(w >= liftcoef);
            min = MIN(minweights[w], minweights[w - liftcoef] + weight);
            minweights[w] = min;
         }
      }
   }
   assert(minweights[0] == 0);

   /* sequentially down-lifts all variables in M_2: */
   for( j = 0; j < nvarsM2; j++ )
   {
      SCIP_Longint weight;
      int liftvar;
      int liftcoef;
      int left;
      int right;
      int middle;
      int z; 

      liftvar = varsM2[j];
      weight = weights[liftvar];
      assert(SCIPisFeasEQ(scip, solvals[liftvar], 1.0));
      assert(liftvar >= 0 && liftvar < nvars);
      assert(weight > 0);
      
      /* uses binary search to find 
       *   z = max { w : 0 <= w <= |M_1| + sum_{k=1}^{i-1} alpha_{j_k}, minweights_[w] <= a_0 - fixedonesweight + a_{j_i}}
       */
      left = 0;
      right = minweightslen;
      while( left < right - 1 )
      {
         middle = (left + right) / 2;
         assert(0 <= middle && middle < minweightslen);
         if( minweights[middle] <= capacity - fixedonesweight + weight ) 
            left = middle;
         else
            right = middle;
      }
      assert(left == right - 1);
      assert(0 <= left && left < minweightslen);
      assert(minweights[left] <= capacity - fixedonesweight + weight );
      assert(left == minweightslen - 1 || minweights[left+1] > capacity - fixedonesweight + weight);

      /* now z = left */
      z = left;
      assert(z >= *liftrhs);
         
      /* calculates lifting coefficients alpha_{j_i} = z - liftrhs */
      liftcoef = z - (*liftrhs);
      liftcoefs[liftvar] = liftcoef;
      assert(liftcoef >= 0);

      /* updates sum of weights of variables fixed to one */
      fixedonesweight -= weight;

      /* updates right-hand side of current valid inequality */
      (*liftrhs) += liftcoef;
      assert(*liftrhs >= alpha0);

      /* minweight table and activity of current valid inequality will not change, if alpha_{j_i} = 0 */
      if( liftcoef == 0 )
         continue;

      /* updates activity of current valid inequality */
      (*cutact) += liftcoef * solvals[liftvar];

      /* enlarges current minweight table:
       *  from minweightlen = |M1| + sum_{k=1}^{i-1} alpha_{j_k} + 1 entries 
       *  to                  |M1| + sum_{k=1}^{i  } alpha_{j_k} + 1 entries
       * and sets minweights_i[w] = infinity for 
       *  w = |M1| + sum_{k=1}^{i-1} alpha_{j_k} + 1 , ... , |M1| + sum_{k=1}^{i} alpha_{j_k} 
       */
      SCIP_CALL( enlargeMinweights(scip, &minweights, &minweightslen, &minweightssize, minweightslen + liftcoef) );
     
      /* updates minweight table: minweight_i+1[w] = 
       *   min{ minweights_i[w], a_{j_i}},                                 if w <  alpha_j_i
       *   min{ minweights_i[w], minweights_i[w - alpha_j_i] + a_j_i},     if w >= alpha_j_i
       */
      for( w = minweightslen - 1; w >= 0; w-- )
      {
         SCIP_Longint min;
         if( w < liftcoef )
         {
            min = MIN(minweights[w], weight);
            minweights[w] = min;
         }
         else
         {
            assert(w >= liftcoef);
            min = MIN(minweights[w], minweights[w - liftcoef] + weight);
            minweights[w] = min;
         }
      }
   }
   assert(fixedonesweight == 0);
   assert(*liftrhs >= alpha0);
   
   /* sequentially up-lifts all variables in R: */
   for( j = 0; j < nvarsR; j++ )
   {
      SCIP_Longint weight;
      int liftvar;
      int liftcoef;
      int z; 

      liftvar = varsR[j];
      weight = weights[liftvar];
      assert(liftvar >= 0 && liftvar < nvars);
      assert(SCIPisFeasEQ(scip, solvals[liftvar], 0.0));
      assert(weight > 0);
      assert(capacity - weight >= 0);
      assert((*liftrhs) + 1 >= minweightslen || minweights[(*liftrhs) + 1] > capacity - weight);
      
      /* sets z = max { w : 0 <= w <= liftrhs, minweights_i[w] <= a_0 - a_{j_i} } = liftrhs, 
       * if minweights_i[liftrhs] <= a_0 - a_{j_i}
       */
      if( minweights[*liftrhs] <= capacity - weight )
      {
         z = *liftrhs;
      } 
      /* uses binary search to find z = max { w : 0 <= w <= liftrhs, minweights_i[w] <= a_0 - a_{j_i} } 
       */
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
            assert(0 <= middle && middle < minweightslen);
            if( minweights[middle] <= capacity - weight ) 
               left = middle;
            else
               right = middle;
         }
         assert(left == right - 1);
         assert(0 <= left && left < minweightslen);
         assert(minweights[left] <= capacity - weight );
         assert(left == minweightslen - 1 || minweights[left+1] > capacity - weight);

         /* now z = left */
         z = left;
         assert(z <= *liftrhs);
      }
         
      /* calculates lifting coefficients alpha_{j_i} = liftrhs - z */
      liftcoef = (*liftrhs) - z;
      liftcoefs[liftvar] = liftcoef;
      assert(liftcoef >= 0 && liftcoef <= *liftrhs);

      /* minweight table and activity of current valid inequality will not change, if alpha_{j_i} = 0 */
      if( liftcoef == 0 )
         continue;

      /* updates activity of current valid inequality */
      (*cutact) += liftcoef * solvals[liftvar];

      /* updates minweight table: minweight_i+1[w] = 
       *   min{ minweight_i[w], a_{j_i}},                                if w <  alpha_j_i
       *   min{ minweight_i[w], minweight_i[w - alpha_j_i] + a_j_i},     if w >= alpha_j_i
       */
      for( w = *liftrhs; w >= 0; w-- )
      {
         SCIP_Longint min;
         if( w < liftcoef )
         {
            min = MIN(minweights[w], weight);
            minweights[w] = min;
         }
         else
         {
            assert(w >= liftcoef);
            min = MIN(minweights[w], minweights[w - liftcoef] + weight);
            minweights[w] = min;
         }
      }
   }
   
   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &sortkeys);
   SCIPfreeBufferArray(scip, &minweights);

   return SCIP_OKAY;
}

/** lifts given minimal cover inequality 
 *  \f[
 *    \sum_{j \in C} x_j \leq |C| - 1 
 *  \f]
 *  valid for 
 *  \f[
 *    S^0 = \{ x \in {0,1}^{|C|} : \sum_{j \in C} a_j x_j \leq a_0 \}
 *  \f]
 *  to a valid inequality 
 *  \f[
 *    \sum_{j \in C} x_j + \sum_{j \in N \setminus C} \alpha_j x_j \leq |C| - 1
 *  \f]
 *  for
 *  \f[ 
 *    S = \{ x \in {0,1}^{|N|} : \sum_{j \in N} a_j x_j \leq a_0 \}; 
 *  \f]
 *  uses superadditive up-lifting for the variables in \f$N \setminus C\f$.
 */ 
static
SCIP_RETCODE superadditiveUpLifting(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   int                   ntightened,         /**< number of variables with tightened upper bound */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  covervars,          /**< cover variables */
   int*                  noncovervars,       /**< noncover variables */
   int                   ncovervars,         /**< number of cover variables */
   int                   nnoncovervars,      /**< number of noncover variables */
   SCIP_Longint          coverweight,        /**< weight of cover */
   SCIP_Real*            liftcoefs,          /**< pointer to store lifting coefficient of vars in knapsack constraint */
   SCIP_Real*            cutact              /**< pointer to store activity of lifted valid inequality */
)
{
   SCIP_Longint* maxweightsums;
   SCIP_Longint* intervalends;
   SCIP_Longint* rhos;
   SCIP_Real* sortkeys;
   SCIP_Longint lambda;
   int j;
   int h;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nvars >= 0);
   assert(weights != NULL);
   assert(capacity >= 0);
   assert(solvals != NULL);
   assert(covervars != NULL);
   assert(noncovervars != NULL);
   assert(ncovervars > 0 && ncovervars <= nvars);
   assert(nnoncovervars >= 0 && nnoncovervars <= nvars - ntightened);
   assert(ncovervars + nnoncovervars == nvars - ntightened);
   assert(liftcoefs != NULL);
   assert(cutact != NULL);

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeys, ncovervars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxweightsums, ncovervars + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &intervalends, ncovervars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhos, ncovervars) );
   
   /* initializes data structures */
   BMSclearMemoryArray(liftcoefs, nvars);
   *cutact = 0.0;

   /* sets lifting coefficient of variables in C, sorts variables in C such that a_1 >= a_2 >= ... >= a_|C| 
    * and calculates activity of current valid inequality 
    */
   for( j = 0; j < ncovervars; j++ )
   {
      assert(liftcoefs[covervars[j]] == 0.0);
      liftcoefs[covervars[j]] = 1.0;
      sortkeys[j] = (SCIP_Real) weights[covervars[j]];
      (*cutact) += solvals[covervars[j]];
   }
   SCIPsortDownRealInt(sortkeys, covervars, ncovervars);
   
   /* calculates weight excess of cover C */
   lambda = coverweight - capacity;
   assert(lambda > 0);

   /* calculates A_h for h = 0,...,|C|, I_h for h = 1,...,|C| and rho_h for h = 1,...,|C| */
   maxweightsums[0] = 0;
   for( h = 1; h <= ncovervars; h++ )
   {
      maxweightsums[h] = maxweightsums[h-1] + weights[covervars[h-1]];
      intervalends[h-1] = maxweightsums[h] - lambda;
      rhos[h-1] = MAX(0, weights[covervars[h-1]] - weights[covervars[0]] + lambda);
   }
   
   /* sorts variables in N\C such that a_{j_1} <= a_{j_2} <= ... <= a_{j_t} */
   for( j = 0; j < nnoncovervars; j++ )
      sortkeys[j] = (SCIP_Real) (weights[noncovervars[j]]);
   SCIPsortRealInt(sortkeys, noncovervars, nnoncovervars);

   /* calculates lifting coefficient for all variables in N\C */
   h = 0;
   for( j = 0; j < nnoncovervars; j++ )
   {
      int liftvar;
      SCIP_Longint weight;
      SCIP_Real liftcoef;

      liftvar = noncovervars[j];
      weight = weights[liftvar];
      
      while( intervalends[h] < weight )
         h++;

      if( h == 0 )
         liftcoef = h;
      else 
      {
         if( weight <= intervalends[h-1] + rhos[h] )
         {
            SCIP_Real tmp1;
            SCIP_Real tmp2;
            tmp1 =  (SCIP_Real) (intervalends[h-1] + rhos[h] - weight);
            tmp2 =  (SCIP_Real) rhos[1];
            liftcoef = h - ( tmp1 / tmp2 );
         }
         else
            liftcoef = h;
      }      

      /* sets lifting coefficient */
      assert(liftcoefs[liftvar] == 0.0);
      liftcoefs[liftvar] = liftcoef;

      /* updates activity of current valid inequality */
      (*cutact) += liftcoef * solvals[liftvar];
   } 

   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &rhos);
   SCIPfreeBufferArray(scip, &intervalends);
   SCIPfreeBufferArray(scip, &maxweightsums);
   SCIPfreeBufferArray(scip, &sortkeys);

   return SCIP_OKAY;
}


/** separates lifted minimal cover inequalities using sequential up- and down-lifting for given knapsack problem */
static
SCIP_RETCODE separateSequLiftedMinimalCoverInequality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that originates the knapsack problem, or NULL */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   int                   ntightened,         /**< number of variables with tightened upper bound */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  mincovervars,       /**< mincover variables */
   int*                  nonmincovervars,    /**< nonmincover variables */
   int                   nmincovervars,      /**< number of mincover variables */
   int                   nnonmincovervars,   /**< number of nonmincover variables */
   SCIP_SOL*             sol,                /**< primal SCIP solution to separate, NULL for current LP solution */
   int*                  ncuts               /**< pointer to add up the number of found cuts */
   )
{
   int* varsC1;
   int* varsC2;
   int* varsF;
   int* varsR;
   int nvarsC1;
   int nvarsC2;
   int nvarsF;
   int nvarsR;
   SCIP_Real cutact;
   int* liftcoefs;
   int liftrhs;

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &varsC1, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsC2, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsF, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsR, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &liftcoefs, nvars) );

   /* gets partition (C_1,C_2) of C, i.e. C_1 & C_2 = C and C_1 cap C_2 = emptyset, with C_1 not empty; choses partition 
    * as follows
    *   C_2 = { j in C : x*_j = 1 } and
    *   C_1 = C\C_2
    */
   getPartitionCovervars(scip, solvals, mincovervars, nmincovervars, varsC1, varsC2, &nvarsC1, &nvarsC2);
   assert(nvarsC1 + nvarsC2 == nmincovervars);
   assert(nmincovervars > 0); 
   assert(nvarsC1 >= 0); /* nvarsC1 > 0 does not always hold, because relaxed knapsack conss may already be violated */

   /* changes partition (C_1,C_2) of minimal cover C, if |C1| = 1, by moving one variable from C2 to C1 */
   if( nvarsC1 < 2 && nvarsC2 > 0)
   {
      SCIP_CALL( changePartitionCovervars(scip, weights, varsC1, varsC2, &nvarsC1, &nvarsC2) );
      assert(nvarsC1 >= 1);
   }
   assert(nvarsC2 == 0 || nvarsC1 >= 1);
   
   /* gets partition (F,R) of N\C, i.e. F & R = N\C and F cap R = emptyset; choses partition as follows 
    *   R = { j in N\C : x*_j = 0 } and 
    *   F = (N\C)\F
    */
   getPartitionNoncovervars(scip, solvals, nonmincovervars, nnonmincovervars, varsF, varsR, &nvarsF, &nvarsR);
   assert(nvarsF + nvarsR == nnonmincovervars);
   assert(nvarsC1 + nvarsC2 + nvarsF + nvarsR == nvars - ntightened);

   /* sorts variables in F, C_2, and R according to the second level lifting sequence that will be used in the sequential 
    * lifting procedure 
    */
   SCIP_CALL( getLiftingSequence(scip, solvals, weights, varsF, varsC2, varsR, nvarsF, nvarsC2, nvarsR) );

   /* lifts minimal cover inequality sum_{j in C_1} x_j <= |C_1| - 1 valid for 
    *
    *    S^0 = { x in {0,1}^|C_1| : sum_{j in C_1} a_j x_j <= a_0 - sum_{j in C_2} a_j } 
    *
    * to a valid inequality sum_{j in C_1} x_j + sum_{j in N\C_1} alpha_j x_j <= |C_1| - 1 + sum_{j in C_2} alpha_j for 
    *
    *      S = { x in {0,1}^|N|   : sum_{j in N}   a_j x_j <= a_0 }, 
    *
    * uses sequential up-lifting for the variables in F, sequential down-lifting for the variable in C_2 and sequential 
    * up-lifting for the variabels in R according to the second level lifting sequence    
    */ 
   SCIP_CALL( sequentialUpAndDownLifting(scip, vars, nvars, ntightened, weights, capacity, solvals, varsC1, varsC2, varsF,
         varsR, nvarsC1, nvarsC2, nvarsF, nvarsR, nvarsC1 - 1, liftcoefs, &cutact, &liftrhs) );

   /* checks, if lifting yielded a violated cut */
   if( SCIPisEfficacious(scip, (cutact - liftrhs)/sqrt((SCIP_Real)MAX(liftrhs, 1))) )
   {
      SCIP_ROW* row;
      char name[SCIP_MAXSTRLEN];
      int j;
            
      /* creates LP row */
      if( cons != NULL )
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_mcseq%"SCIP_LONGINT_FORMAT"", SCIPconsGetName(cons),
            SCIPconshdlrGetNCutsFound(SCIPconsGetHdlr(cons)));
      else
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nn_mcseq_%"SCIP_LONGINT_FORMAT"", *ncuts); 
      
      SCIP_CALL( SCIPcreateEmptyRow(scip, &row, name, -SCIPinfinity(scip), (SCIP_Real)liftrhs, 
            cons != NULL ? SCIPconsIsLocal(cons) : FALSE, FALSE,
            cons != NULL ? SCIPconsIsRemovable(cons) : TRUE) );
      
      /* adds all variables in the knapsack constraint with calculated lifting coefficient to the cut */
      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
      assert(nvarsC1 + nvarsC2 + nvarsF + nvarsR == nvars - ntightened);
      for( j = 0; j < nvarsC1; j++ )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, vars[varsC1[j]], 1.0) ); 
      }
      for( j = 0; j < nvarsC2; j++ )
      {
         if( liftcoefs[varsC2[j]] > 0 )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, row, vars[varsC2[j]], (SCIP_Real)liftcoefs[varsC2[j]]) );
         }
      }
      for( j = 0; j < nvarsF; j++ )
      {
         if( liftcoefs[varsF[j]] > 0 )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, row, vars[varsF[j]], (SCIP_Real)liftcoefs[varsF[j]]) );
         }
      }
      for( j = 0; j < nvarsR; j++ )
      {
         if( liftcoefs[varsR[j]] > 0 )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, row, vars[varsR[j]], (SCIP_Real)liftcoefs[varsR[j]]) );
         }
      }
      SCIP_CALL( SCIPflushRowExtensions(scip, row) );
            
      /* checks, if cut is violated enough */
      if( SCIPisCutEfficacious(scip, sol, row) )
      {
         if( cons != NULL )
         {
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
         }
         SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE) );
         (*ncuts)++;
      }
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }
   
   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &liftcoefs);
   SCIPfreeBufferArray(scip, &varsR);
   SCIPfreeBufferArray(scip, &varsF);
   SCIPfreeBufferArray(scip, &varsC2);
   SCIPfreeBufferArray(scip, &varsC1);

   return SCIP_OKAY;
}

/** separates lifted extended weight inequalities using sequential up- and down-lifting for given knapsack problem */
static
SCIP_RETCODE separateSequLiftedExtendedWeightInequality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that originates the knapsack problem, or NULL */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   int                   ntightened,         /**< number of variables with tightened upper bound */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  feassetvars,        /**< variables in feasible set */
   int*                  nonfeassetvars,     /**< variables not in feasible set */
   int                   nfeassetvars,       /**< number of variables in feasible set */
   int                   nnonfeassetvars,    /**< number of variables not in feasible set */
   SCIP_SOL*             sol,                /**< primal SCIP solution to separate, NULL for current LP solution */
   int*                  ncuts               /**< pointer to add up the number of found cuts */
   )
{
   int* varsT1;
   int* varsT2;
   int* varsF;
   int* varsR;
   int* liftcoefs;
   SCIP_Real cutact;
   int nvarsT1;
   int nvarsT2;
   int nvarsF;
   int nvarsR;
   int liftrhs;
   int j;

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &varsT1, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsT2, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsF, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsR, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &liftcoefs, nvars) );

   /* gets partition (T_1,T_2) of T, i.e. T_1 & T_2 = T and T_1 cap T_2 = emptyset, with T_1 not empty; choses partition 
    * as follows
    *   T_2 = { j in T : x*_j = 1 } and
    *   T_1 = T\T_2
    */
   getPartitionCovervars(scip, solvals, feassetvars, nfeassetvars, varsT1, varsT2, &nvarsT1, &nvarsT2);
   assert(nvarsT1 + nvarsT2 == nfeassetvars);
   
   /* changes partition (T_1,T_2) of feasible set T, if |T1| = 0, by moving one variable from T2 to T1 */
   if( nvarsT1 == 0 && nvarsT2 > 0)
   {
      SCIP_CALL( changePartitionFeasiblesetvars(scip, weights, varsT1, varsT2, &nvarsT1, &nvarsT2) );
      assert(nvarsT1 == 1);
   }
   assert(nvarsT2 == 0 || nvarsT1 > 0);
   
   /* gets partition (F,R) of N\T, i.e. F & R = N\T and F cap R = emptyset; choses partition as follows 
    *   R = { j in N\T : x*_j = 0 } and 
    *   F = (N\T)\F
    */
   getPartitionNoncovervars(scip, solvals, nonfeassetvars, nnonfeassetvars, varsF, varsR, &nvarsF, &nvarsR);
   assert(nvarsF + nvarsR == nnonfeassetvars);
   assert(nvarsT1 + nvarsT2 + nvarsF + nvarsR == nvars - ntightened);

   /* sorts variables in F, T_2, and R according to the second level lifting sequence that will be used in the sequential 
    * lifting procedure (the variable removed last from the initial cover does not have to be lifted first, therefore it 
    * is included in the sorting routine)
    */
   SCIP_CALL( getLiftingSequence(scip, solvals, weights, varsF, varsT2, varsR, nvarsF, nvarsT2, nvarsR) );

   /* lifts extended weight inequality sum_{j in T_1} x_j <= |T_1| valid for 
    *
    *    S^0 = { x in {0,1}^|T_1| : sum_{j in T_1} a_j x_j <= a_0 - sum_{j in T_2} a_j } 
    *
    * to a valid inequality sum_{j in T_1} x_j + sum_{j in N\T_1} alpha_j x_j <= |T_1| + sum_{j in T_2} alpha_j for 
    *
    *      S = { x in {0,1}^|N|   : sum_{j in N}   a_j x_j <= a_0 }, 
    *
    * uses sequential up-lifting for the variables in F, sequential down-lifting for the variable in T_2 and sequential 
    * up-lifting for the variabels in R according to the second level lifting sequence    
    */ 
   SCIP_CALL( sequentialUpAndDownLifting(scip, vars, nvars, ntightened, weights, capacity, solvals, varsT1, varsT2, varsF, varsR, 
         nvarsT1, nvarsT2, nvarsF, nvarsR, nvarsT1, liftcoefs, &cutact, &liftrhs) );

   /* checks, if lifting yielded a violated cut */
   if( SCIPisEfficacious(scip, (cutact - liftrhs)/sqrt((SCIP_Real)MAX(liftrhs, 1))) )
   {
      SCIP_ROW* row;
      char name[SCIP_MAXSTRLEN];
            
      /* creates LP row */
      if( cons != NULL )
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_ewseq%"SCIP_LONGINT_FORMAT"", SCIPconsGetName(cons),
            SCIPconshdlrGetNCutsFound(SCIPconsGetHdlr(cons)));
      else
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nn_ewseq_%"SCIP_LONGINT_FORMAT"", *ncuts); 
      
      SCIP_CALL( SCIPcreateEmptyRow(scip, &row, name, -SCIPinfinity(scip), (SCIP_Real)liftrhs, 
            cons != NULL ? SCIPconsIsLocal(cons) : FALSE, FALSE,
            cons != NULL ? SCIPconsIsRemovable(cons) : TRUE) );
            
      /* adds all variables in the knapsack constraint with calculated lifting coefficient to the cut */
      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
      assert(nvarsT1 + nvarsT2 + nvarsF + nvarsR == nvars - ntightened);
      for( j = 0; j < nvarsT1; j++ )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, vars[varsT1[j]], 1.0) ); 
      }
      for( j = 0; j < nvarsT2; j++ )
      {
         if( liftcoefs[varsT2[j]] > 0 )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, row, vars[varsT2[j]], (SCIP_Real)liftcoefs[varsT2[j]]) );
         }
      }
      for( j = 0; j < nvarsF; j++ )
      {
         if( liftcoefs[varsF[j]] > 0 )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, row, vars[varsF[j]], (SCIP_Real)liftcoefs[varsF[j]]) );
         }
      }
      for( j = 0; j < nvarsR; j++ )
      {
         if( liftcoefs[varsR[j]] > 0 )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, row, vars[varsR[j]], (SCIP_Real)liftcoefs[varsR[j]]) );
         }
      }
      SCIP_CALL( SCIPflushRowExtensions(scip, row) );
            
      /* checks, if cut is violated enough */
      if( SCIPisCutEfficacious(scip, sol, row) )
      {
         if( cons != NULL )
         {
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
         }
         SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE) );
         (*ncuts)++;
      }
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &liftcoefs);
   SCIPfreeBufferArray(scip, &varsR);
   SCIPfreeBufferArray(scip, &varsF);
   SCIPfreeBufferArray(scip, &varsT2);
   SCIPfreeBufferArray(scip, &varsT1);
   
   return SCIP_OKAY;
}

/** separates lifted minimal cover inequalities using superadditive up-lifting for given knapsack problem */
static
SCIP_RETCODE separateSupLiftedMinimalCoverInequality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that originates the knapsack problem, or NULL */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   int                   ntightened,         /**< number of variables with tightened upper bound */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  mincovervars,       /**< mincover variables */
   int*                  nonmincovervars,    /**< nonmincover variables */
   int                   nmincovervars,      /**< number of mincover variables */
   int                   nnonmincovervars,   /**< number of nonmincover variables */
   SCIP_Longint          mincoverweight,     /**< weight of minimal cover */
   SCIP_SOL*             sol,                /**< primal SCIP solution to separate, NULL for current LP solution */
   int*                  ncuts               /**< pointer to add up the number of found cuts */
   )
{
   SCIP_Real* realliftcoefs;
   SCIP_Real cutact;
   int liftrhs;

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &realliftcoefs, nvars) );

   /* lifts minimal cover inequality sum_{j in C} x_j <= |C| - 1 valid for 
    * 
    *    S^0 = { x in {0,1}^|C| : sum_{j in C} a_j x_j <= a_0 } 
    *
    * to a valid inequality sum_{j in C} x_j + sum_{j in N\C} alpha_j x_j <= |C| - 1 for 
    *
    *      S = { x in {0,1}^|N|   : sum_{j in N}   a_j x_j <= a_0 }, 
    *
    * uses superadditive up-lifting for the variables in N\C.    
    */ 
   SCIP_CALL( superadditiveUpLifting(scip, vars, nvars, ntightened, weights, capacity, solvals, mincovervars, 
         nonmincovervars, nmincovervars, nnonmincovervars, mincoverweight, realliftcoefs, &cutact) );
   liftrhs = nmincovervars - 1;

   /* checks, if lifting yielded a violated cut */
   if( SCIPisEfficacious(scip, (cutact - liftrhs)/sqrt((SCIP_Real)MAX(liftrhs, 1))) )
   {
      SCIP_ROW* row;
      char name[SCIP_MAXSTRLEN];
      int j;
            
      /* creates LP row */
      if( cons != NULL )
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_mcsup%"SCIP_LONGINT_FORMAT"", SCIPconsGetName(cons),
            SCIPconshdlrGetNCutsFound(SCIPconsGetHdlr(cons)));
      else
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nn_mcsup_%"SCIP_LONGINT_FORMAT"", *ncuts); 

      SCIP_CALL( SCIPcreateEmptyRow(scip, &row, name, -SCIPinfinity(scip), (SCIP_Real)liftrhs, 
            cons != NULL ? SCIPconsIsLocal(cons) : FALSE, FALSE,
            cons != NULL ? SCIPconsIsRemovable(cons) : TRUE) );
            
      /* adds all variables in the knapsack constraint with calculated lifting coefficient to the cut */
      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
      assert(nmincovervars + nnonmincovervars == nvars - ntightened);
      for( j = 0; j < nmincovervars; j++ )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, vars[mincovervars[j]], 1.0) ); 
      }
      for( j = 0; j < nnonmincovervars; j++ )
      {
         assert(SCIPisFeasGE(scip, realliftcoefs[nonmincovervars[j]], 0.0));
         if( SCIPisFeasGT(scip, realliftcoefs[nonmincovervars[j]], 0.0) ) 
         {
            SCIP_CALL( SCIPaddVarToRow(scip, row, vars[nonmincovervars[j]], realliftcoefs[nonmincovervars[j]]) );
         }
      }
      SCIP_CALL( SCIPflushRowExtensions(scip, row) );
            
      /* checks, if cut is violated enough */
      if( SCIPisCutEfficacious(scip, sol, row) )
      {
         if( cons != NULL )
         {
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
         }
         SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE) );
         (*ncuts)++;
      }
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &realliftcoefs);

   return SCIP_OKAY;
}

/** coverts given cover C to a minimal cover by removing variables in the reverse order in which the variables were chosen
 *  to be in C, i.e. in the order of non-increasing (1 - x*_j)/a_j, if the transformed separation problem was used to find
 *  C and in the order of non-increasing (1 - x*_j), if the modified transformed separation problem was used to find C; 
 *  note that all variables with x*_j = 1 will be removed last
 */
static
SCIP_RETCODE makeCoverMinimal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  covervars,          /**< pointer to store cover variables */
   int*                  noncovervars,       /**< pointer to store noncover variables */
   int*                  ncovervars,         /**< pointer to store number of cover variables */
   int*                  nnoncovervars,      /**< pointer to store number of noncover variables */
   SCIP_Longint*         coverweight,        /**< pointer to store weight of cover */
   SCIP_Bool             modtransused        /**< TRUE if mod trans sepa prob was used to find cover */    
)
{
   SORTKEYPAIR** sortkeypairs;
   SCIP_Longint minweight;
   int nsortkeypairs;
   int minweightidx;
   int j;
   int k;

   assert(scip != NULL);
   assert(covervars != NULL);
   assert(noncovervars != NULL);
   assert(ncovervars != NULL);
   assert(*ncovervars > 0);
   assert(nnoncovervars != NULL);
   assert(*nnoncovervars >= 0);
   assert(coverweight != NULL);
   assert(*coverweight > 0);
   assert(*coverweight > capacity);

   /* allocates temporary memory */
   nsortkeypairs = *ncovervars;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sortkeypairs, nsortkeypairs) );

   /* sorts C in the reverse order in which the variables were chosen to be in the cover, i.e. 
    *   such that (1 - x*_1)/a_1 >= ... >= (1 - x*_|C|)/a_|C|,  if          trans separation problem was used to find C 
    *   such that (1 - x*_1)     >= ... >= (1 - x*_|C|),        if modified trans separation problem was used to find C 
    * note that all variables with x*_j = 1 are in the end of the sorted C, so they will be removed last from C
    */
   assert(*ncovervars == nsortkeypairs);
   if( modtransused )
   {
      for( j = 0; j < *ncovervars; j++ )
      {
         SCIP_CALL( SCIPallocBlockMemory(scip, &sortkeypairs[j]) );
         
         sortkeypairs[j]->key1 = solvals[covervars[j]]; 
         sortkeypairs[j]->key2 = (SCIP_Real) weights[covervars[j]]; 
      }
   }
   else
   {
      for( j = 0; j < *ncovervars; j++ )
      {
         SCIP_CALL( SCIPallocBlockMemory(scip, &sortkeypairs[j]) );
         
         sortkeypairs[j]->key1 = (solvals[covervars[j]] - 1.0) / weights[covervars[j]]; 
         sortkeypairs[j]->key2 = (SCIP_Real) (-weights[covervars[j]]); 
      }
   }
   SCIPsortPtrInt((void**)sortkeypairs, covervars, compSortkeypairs, *ncovervars);

   /* gets j' with a_j' = min{ a_j : j in C } */
   minweightidx = 0;
   minweight = weights[covervars[minweightidx]];
   for( j = 1; j < *ncovervars; j++ )
   {
      if( weights[covervars[j]] <= minweight )
      {
         minweightidx = j;
         minweight = weights[covervars[minweightidx]];
      }
   }
   assert(minweightidx >= 0 && minweightidx < *ncovervars);
   assert(minweight > 0 && minweight <= *coverweight);

   /* removes variables from C until the remaining variables form a minimal cover */
   for( j = 0; j < *ncovervars && ((*coverweight) - minweight > capacity); j++ )
   {
      assert(minweightidx >= j);
      assert(checkMinweightidx(weights, capacity, covervars, *ncovervars, *coverweight, minweightidx, j));

      /* if sum_{i in C} a_i - a_j <= a_0, j cannot be removed from C */
      if( (*coverweight) - weights[covervars[j]] <= capacity )
         continue;
      
      /* adds j to N\C */
      noncovervars[*nnoncovervars] = covervars[j];
      (*nnoncovervars)++;

      /* removes j from C */
      (*coverweight) -= weights[covervars[j]];
      for( k = j; k < (*ncovervars) - 1; k++ )
         covervars[k] = covervars[k+1];
      (*ncovervars)--;

      /* updates j' with a_j' = min{ a_j : j in C } */
      if( j == minweightidx )
      {
         minweightidx = 0;
         minweight = weights[covervars[minweightidx]];
         for( k = 1; k < *ncovervars; k++ )
         {
            if( weights[covervars[k]] <= minweight )
            {
               minweightidx = k;
               minweight = weights[covervars[minweightidx]];
            }
         }
         assert(minweight > 0 && minweight <= *coverweight);
         assert(minweightidx >= 0 && minweightidx < *ncovervars);
      }
      else
      {
         assert(minweightidx > j);
         minweightidx--;
      }
      /* updates j */
      j--;
   }
   assert((*coverweight) > capacity);
   assert((*coverweight) - minweight <= capacity);
   
   /* frees temporary memory */
   for( j = nsortkeypairs-1; j >= 0; j-- )
      SCIPfreeBlockMemory(scip, &sortkeypairs[j]);
   SCIPfreeBlockMemoryArray(scip, &sortkeypairs, nsortkeypairs);
   
   return SCIP_OKAY;
}

/** coverts given initial cover C_init to a feasible set by removing variables in the reverse order in which 
 *  they were chosen to be in C_init:
 *   non-increasing (1 - x*_j)/a_j,   if          transformed separation problem was used to find C_init 
 *   non-increasing (1 - x*_j),       if modified transformed separation problem was used to find C_init.
 *  separates lifted extended weight inequalities using sequential up- and down-lifting for this feasible set 
 *  and all subsequent feasible sets.
 */
static
SCIP_RETCODE getFeasibleSet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that originates the knapsack problem */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   int                   ntightened,         /**< number of variables with tightened upper bound */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  covervars,          /**< pointer to store cover variables */
   int*                  noncovervars,       /**< pointer to store noncover variables */
   int*                  ncovervars,         /**< pointer to store number of cover variables */
   int*                  nnoncovervars,      /**< pointer to store number of noncover variables */
   SCIP_Longint*         coverweight,        /**< pointer to store weight of cover */
   SCIP_Bool             modtransused,       /**< TRUE if mod trans sepa prob was used to find cover */    
   SCIP_SOL*             sol,                /**< primal SCIP solution to separate, NULL for current LP solution */
   int*                  ncuts               /**< pointer to add up the number of found cuts */
)
{
   SCIP_Real* sortkeys;
   int j;
   int k;

   assert(scip != NULL);
   assert(covervars != NULL);
   assert(noncovervars != NULL);
   assert(ncovervars != NULL);
   assert(*ncovervars > 0);
   assert(nnoncovervars != NULL);
   assert(*nnoncovervars >= 0);
   assert(coverweight != NULL);
   assert(*coverweight > 0);
   assert(*coverweight > capacity);
   assert(*ncovervars + *nnoncovervars == nvars - ntightened);

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeys, *ncovervars) );

   /* sorts C in the reverse order in which the variables were chosen to be in the cover, i.e. 
    *   such that (1 - x*_1)/a_1 >= ... >= (1 - x*_|C|)/a_|C|,  if          trans separation problem was used to find C 
    *   such that (1 - x*_1)     >= ... >= (1 - x*_|C|),        if modified trans separation problem was used to find C 
    * note that all variables with x*_j = 1 are in the end of the sorted C, so they will be removed last from C
    */
   if( modtransused )
   {
      for( j = 0; j < *ncovervars; j++ )
      {
         sortkeys[j] = solvals[covervars[j]];
         assert(SCIPisFeasGE(scip, sortkeys[j], 0.0));
      }
   }
   else
   {
      for( j = 0; j < *ncovervars; j++ )
      {
         sortkeys[j] = (solvals[covervars[j]] - 1.0) / weights[covervars[j]];
         assert(SCIPisFeasLE(scip, sortkeys[j], 0.0));
      }
   }
   SCIPsortRealInt(sortkeys, covervars, *ncovervars);

   /* removes variables from C_init and separates lifted extended weight inequalities using sequential up- and down-lifting;
    * in addition to an extended weight inequality this gives cardinality inequalities  
    */
   while( *ncovervars >= 2 ) 
   {
      /* adds first element of C_init to N\C_init */
      noncovervars[*nnoncovervars] = covervars[0];
      (*nnoncovervars)++;

      /* removes first element from C_init */
      (*coverweight) -= weights[covervars[0]];
      for( k = 0; k < (*ncovervars) - 1; k++ )
         covervars[k] = covervars[k+1];
      (*ncovervars)--;

      assert(*ncovervars + *nnoncovervars == nvars - ntightened);
      if( (*coverweight) <= capacity )
      {
         SCIP_CALL( separateSequLiftedExtendedWeightInequality(scip, cons, vars, nvars, ntightened, weights, capacity, solvals, 
               covervars, noncovervars, *ncovervars, *nnoncovervars, sol, ncuts) );
      }
   }

   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &sortkeys);
   
   return SCIP_OKAY;
}

/** separates different classes of valid inequalities for the 0-1 knapsack problem */ 
SCIP_RETCODE SCIPseparateKnapsackCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that originates the knapsack problem, or NULL */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_SOL*             sol,                /**< primal SCIP solution to separate, NULL for current LP solution */
   int*                  ncuts               /**< pointer to add up the number of found cuts */
   )
{
   SCIP_Real* solvals;
   int* covervars;
   int* noncovervars;
   SCIP_Bool coverfound;
   SCIP_Bool fractional;
   SCIP_Bool modtransused;
   SCIP_Longint coverweight;
   int ncovervars;
   int nnoncovervars;
   int ntightened;

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
   
   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &covervars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &noncovervars, nvars) );

   /* gets solution values of all problem variables */
   SCIP_CALL( SCIPgetSolVals(scip, sol, nvars, vars, solvals) );

#ifdef SCIP_DEBUG
   {
      int i;

      SCIPdebugMessage("separate cuts for knapsack constraint originated by cons <%s>:\n", 
         cons == NULL ? "-" : SCIPconsGetName(cons));
      for( i = 0; i < nvars; ++i )
      {
         SCIPdebugPrintf("%+"SCIP_LONGINT_FORMAT"<%s>(%g)", weights[i], SCIPvarGetName(vars[i]), solvals[i]);
      }
      SCIPdebugPrintf(" <= %"SCIP_LONGINT_FORMAT"\n", capacity);
   }
#endif

   /* LMCI1 (lifted minimal cover inequalities using sequential up- and down-lifting) 
    * (and LMCI2 (lifted minimal cover inequalities using superadditive up-lifting)) 
    */

   /* gets a most violated initial cover C_init ( sum_{j in C_init} a_j > a_0 ) by using the 
    * MODIFIED transformed separation problem and taking into account the following fixing: 
    *   j in C_init,           if j in N_1 = {j in N : x*_j = 1} and 
    *   j in N\C_init,         if j in N_0 = {j in N : x*_j = 0}, 
    * if one exists 
    */
   SCIPdebugMessage("separate LMCI1 cuts:\n");
   modtransused = TRUE;
   SCIP_CALL( getCover(scip, vars, nvars, weights, capacity, solvals, covervars, noncovervars, &ncovervars, 
         &nnoncovervars, &coverweight, &coverfound, modtransused, &ntightened, &fractional) );
   assert(!coverfound || !fractional || ncovervars + nnoncovervars == nvars - ntightened);

   /* if x* is not fractional we stop the separation routine */
   if( !fractional )
      goto TERMINATE;
   
   /* if no cover was found we stop the separation routine for lifted minimal cover inequality */
   if( coverfound )
   {
      /* coverts initial cover C_init to a minimal cover C by removing variables in the reverse order in which the 
       * variables were chosen to be in C_init; note that variables with x*_j = 1 will be removed last
       */
      SCIP_CALL( makeCoverMinimal(scip, weights, capacity, solvals, covervars, noncovervars, &ncovervars, 
            &nnoncovervars, &coverweight, modtransused) );

      /* separates lifted minimal cover inequalities using sequential up- and down-lifting */
      SCIP_CALL( separateSequLiftedMinimalCoverInequality(scip, cons, vars, nvars, ntightened, weights, capacity, solvals,
            covervars, noncovervars, ncovervars, nnoncovervars, sol, ncuts) );

      if( USESUPADDLIFT ) /*lint !e506 !e774*/
      { 
         SCIPdebugMessage("separate LMCI2 cuts:\n"); 
         /* separates lifted minimal cover inequalities using superadditive up-lifting */
         SCIP_CALL( separateSupLiftedMinimalCoverInequality(scip, cons, vars, nvars, ntightened, weights, capacity, 
               solvals, covervars, noncovervars, ncovervars, nnoncovervars, coverweight, sol, ncuts) );
      }
   }

   /* LEWI (lifted extended weight inequalities using sequential up- and down-lifting) */

   /* gets a most violated initial cover C_init ( sum_{j in C_init} a_j > a_0 ) by using the 
    * transformed separation problem and taking into account the following fixing: 
    *   j in C_init,           if j in N_1 = {j in N : x*_j = 1} and 
    *   j in N\C_init,         if j in N_0 = {j in N : x*_j = 0}, 
    * if one exists 
    */
   SCIPdebugMessage("separate LEWI cuts:\n"); 
   modtransused = FALSE; 
   SCIP_CALL( getCover(scip, vars, nvars, weights, capacity, solvals, covervars, noncovervars, &ncovervars, 
            &nnoncovervars, &coverweight, &coverfound, modtransused, &ntightened, &fractional) );
   assert(fractional);
   assert(!coverfound || ncovervars + nnoncovervars == nvars - ntightened);

   /* if no cover was found we stop the separation routine */
   if( !coverfound )
      goto TERMINATE;
      
   /* coverts initial cover C_init to a feasible set by removing variables in the reverse order in which 
    * they were chosen to be in C_init and separates lifted extended weight inequalities using sequential 
    * up- and down-lifting for this feasible set and all subsequent feasible sets.
    */
   SCIP_CALL( getFeasibleSet(scip, cons, vars, nvars, ntightened, weights, capacity, solvals, covervars, noncovervars, 
         &ncovervars, &nnoncovervars, &coverweight, modtransused, sol, ncuts) );

 TERMINATE:
   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &noncovervars);
   SCIPfreeBufferArray(scip, &covervars);
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

   binvars = SCIPgetVars(scip);

   /* all variables which are of integral type can be potentially of binary type; this can be checked via the method SCIPvarIsBinary(var) */
   nbinvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);

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

      if( SCIPvarIsBinary(var) && SCIPvarIsActive(var) )
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
            if( SCIPvarIsBinary(zvlb[j]) && SCIPvarIsActive(zvlb[j]) && REALABS(bvlb[j]) <= MAXABSVBCOEF )
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
            if( SCIPvarIsBinary(zvub[j]) && SCIPvarIsActive(zvub[j]) && REALABS(bvub[j]) <= MAXABSVBCOEF )
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
         SCIP_CALL( SCIPseparateKnapsackCuts(scip, cons, consvars, nconsvars, consvals, capacity, sol, ncuts) );
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
   SCIP_SOL*             sol,                /**< primal SCIP solution, NULL for current LP solution */
   SCIP_Bool             sepacuts,           /**< should knapsack cuts be separated? */
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
   else if( sepacuts )
   {
      SCIP_CALL( SCIPseparateKnapsackCuts(scip, cons, consdata->vars, consdata->nvars, consdata->weights, 
            consdata->capacity, sol, ncuts) );
   }
   
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
   assert(SCIPvarIsBinary(var));
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
      consdata->negcliquepartitioned = FALSE;
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

   consdata->propagated = FALSE;
   consdata->presolved = FALSE;
   consdata->sorted = FALSE;

   /* try to use old clique partitions */
   if( consdata->cliquepartitioned )
   {
      assert(consdata->cliquepartition != NULL);
      /* if the clique number is equal to the number of variables we have only cliques with one element, so we don't 
       * change the clique number */
      if( consdata->cliquepartition[consdata->nvars - 1] != consdata->nvars - 1 )
      {
         int oldcliqenum;
         
         oldcliqenum = consdata->cliquepartition[pos];
         consdata->cliquepartition[pos] = consdata->cliquepartition[consdata->nvars-1];
	 
	 /* the following if and else cases assure that we have increasing cliquenumbers */
	 if( consdata->cliquepartition[pos] > pos )
	    consdata->cliquepartitioned = FALSE; /* recalculate the clique partition after a coefficient was removed */
	 else
	 {
	    int i;
	    int cliquenumbefore;

            /* if the old clique number was greater than the new one we have to check that before a bigger cliquenumber 
             * occurs the same as the old one is still in the cliquepartition */
            if( oldcliqenum > consdata->cliquepartition[pos] )
            {
               for( i = 0; i < consdata->nvars; ++i )
                  if( oldcliqenum == consdata->cliquepartition[i] )
                     break;
                  else if( oldcliqenum < consdata->cliquepartition[i] )
                  {
                     consdata->cliquepartitioned = FALSE; /* recalculate the clique partition after a coefficient was removed */
                     break;
                  }
               /* if we reached the end in the for loop, it means we have deleted the last element of the clique with
                * the biggest index, so decrease the number of cliques 
                */
               if( i == consdata->nvars )
                  --(consdata->ncliques);
            }
            /* if the old clique number was smaller than the new one we have to check the front for an element with 
             * cliquenumber minus 1 */
            else if( oldcliqenum < consdata->cliquepartition[pos] )
            {
               cliquenumbefore = consdata->cliquepartition[pos] - 1;
               for( i = pos - 1; i >= 0 && i >= cliquenumbefore && consdata->cliquepartition[i] < cliquenumbefore; --i );

               if( i < cliquenumbefore )
                  consdata->cliquepartitioned = FALSE; /* recalculate the clique partition after a coefficient was removed */
            }
            /* if we deleted the last element of the clique with biggest index, we have to decrease the clique number */
            else if( pos == consdata->nvars - 1)
            {
               cliquenumbefore = consdata->cliquepartition[pos];
               for( i = pos - 1; i >= 0 && i >= cliquenumbefore && consdata->cliquepartition[i] < cliquenumbefore; --i );

               if( i < cliquenumbefore )
                  --(consdata->ncliques);
            }
            /* if the old clique number is equal to the new one the cliquepartition should be ok */
	 }
      }
      else
	 --(consdata->ncliques);
   }

   if( consdata->negcliquepartitioned )
   {
      assert(consdata->negcliquepartition != NULL);
      /* if the clique number is equal to the number of variables we have only cliques with one element, so we don't 
       * change the clique number */
      if( consdata->negcliquepartition[consdata->nvars-1] != consdata->nvars - 1 )
      {
         int oldcliqenum;
         
         oldcliqenum = consdata->negcliquepartition[pos];
         consdata->negcliquepartition[pos] = consdata->negcliquepartition[consdata->nvars-1];
	 
	 /* the following if and else cases assure that we have increasing cliquenumbers */
	 if( consdata->negcliquepartition[pos] > pos )
	    consdata->negcliquepartitioned = FALSE; /* recalculate the clique partition after a coefficient was removed */
	 else
	 {
	    int i;
	    int cliquenumbefore;

            /* if the old clique number was greater than the new one we have to check that, before a bigger cliquenumber
             * occurs, the same as the old one occurs */
            if( oldcliqenum > consdata->negcliquepartition[pos] )
            {
               for( i = 0; i < consdata->nvars; ++i )
                  if( oldcliqenum == consdata->negcliquepartition[i] )
                     break;
                  else if( oldcliqenum < consdata->negcliquepartition[i] )
                  {
                     consdata->negcliquepartitioned = FALSE; /* recalculate the negated clique partition after a coefficient was removed */
                     break;
                  }
               /* if we reached the end in the for loop, it means we have deleted the last element of the clique with
                * the biggest index, so decrease the number of negeated cliques 
                */
               if( i == consdata->nvars )
                  --(consdata->nnegcliques);
            }
            /* if the old clique number was smaller than the new one we have to check the front for an element with 
             * cliquenumber minus 1 */
            else if( oldcliqenum < consdata->negcliquepartition[pos] )
            {
               cliquenumbefore = consdata->negcliquepartition[pos] - 1;
               for( i = pos - 1; i >= 0 && i >= cliquenumbefore && consdata->negcliquepartition[i] < cliquenumbefore; --i );

               if( i < cliquenumbefore )
                  consdata->negcliquepartitioned = FALSE; /* recalculate the negated clique partition after a coefficient was removed */
            }
            /* if we deleted the last element of the clique with biggest index, we have to decrease the clique number */
            else if( pos == consdata->nvars - 1)
            {
               cliquenumbefore = consdata->negcliquepartition[pos];
               for( i = pos - 1; i >= 0 && i >= cliquenumbefore && consdata->negcliquepartition[i] < cliquenumbefore; --i );

               if( i < cliquenumbefore )
                  --(consdata->nnegcliques);
            }
            /* otherwise if the old clique number is equal to the new one the cliquepartition should be ok */
	 }
      }
      else
	 --(consdata->nnegcliques);
   }

   --(consdata->nvars);

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
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_CONSDATA* consdata;
   int v;
   int prev;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *cutoff = FALSE;

   if( consdata->merged ) 
      return SCIP_OKAY;

   if( consdata->nvars <= 1 )
   {
      consdata->merged = TRUE;
      return SCIP_OKAY;
   }

   assert(consdata->vars != NULL || consdata->nvars == 0);

   /* sorting array after indices of variables, that's only for faster merging */ 
   SCIPsortPtrPtrLongIntInt((void**)consdata->vars, (void**)consdata->eventdatas, consdata->weights, 
      consdata->cliquepartition, consdata->negcliquepartition, SCIPvarCompActiveAndNegated, consdata->nvars);

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
      assert(SCIPvarIsBinary(var1));
      assert(SCIPvarIsActive(var1) || SCIPvarGetStatus(var1) == SCIP_VARSTATUS_NEGATED);
      if( SCIPvarGetStatus(var1) == SCIP_VARSTATUS_NEGATED )
      {
         var1 = SCIPvarGetNegatedVar(var1);
         negated1 = TRUE;
      }
      assert(var1 != NULL);
      
      var2 = consdata->vars[prev];
      assert(SCIPvarIsBinary(var2));
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
	       assert(prev == 0 || ((prev > 0) && (SCIPvarIsActive(consdata->vars[prev - 1]) || SCIPvarGetStatus(consdata->vars[prev - 1]) == SCIP_VARSTATUS_NEGATED)) );
               /* either that was the last pair or both, the negated and "normal" variable in front doesn't match var1, so the order is irrelevant */
               if( prev == 0 || (var1 != consdata->vars[prev - 1] && var1 != SCIPvarGetNegatedVar(consdata->vars[prev - 1])) )
                  --prev;
               else /* we need to let v at the same position*/
               {
                  consdata->cliquesadded = FALSE; /* reduced capacity might lead to larger cliques */
                  /* don't decrease v, the same variable may exist up front */
                  continue;
               }
            }
         }
         consdata->cliquesadded = FALSE; /* reduced capacity might lead to larger cliques */
      }
      v = prev;
   }

   consdata->merged = TRUE;

   /* check infeasibility */
   if( consdata->onesweightsum > consdata->capacity )
   {
      SCIPdebugMessage("merge multiples detected cutoff.\n");
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** in case the knapsack constraint is independent of every else, solve the knapsack problem (exactly) and apply the
 *  fixings (dual reductions) 
 */
static
SCIP_RETCODE dualPresolving(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int*                  nfixedvars,         /**< pointer to count number of fixings */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints  */
   SCIP_Bool*            deleted             /**< pointer to store if the constraint is deleted */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* profits;
   int* solitems;
   int* nonsolitems;
   int* items;
   SCIP_Real solval;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;
   SCIP_Bool applicable;
   int nsolitems;      
   int nnonsolitems;      
   int nvars;
   int v;

   assert(!SCIPconsIsModifiable(cons));

   /* constraints for which the check flag is set to FALSE, did not contribute to the lock numbers; therefore, we cannot
    * use the locks to decide for a dual reduction using this constraint; for example after a restart the cuts which are
    * added to the problems have the check flag set to FALSE 
    */
   if( !SCIPconsIsChecked(cons) )
      return SCIP_OKAY;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;
   vars = consdata->vars;
      
   SCIP_CALL( SCIPallocBufferArray(scip, &profits, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &items, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solitems, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nonsolitems, nvars) );

   applicable = TRUE;

   /* check if we can apply the dual reduction; this can be done if the knapsack has the only locks on this constraint;
    * collect object values which are the profits of the knapsack problem 
    */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Bool negated;
      
      var = vars[v];
      assert(var != NULL);

      /* the variable should not be (globally) fixed */
      assert(SCIPvarGetLbGlobal(var) < 0.5 && SCIPvarGetUbGlobal(var) > 0.5);
      
      if( SCIPvarGetNLocksDown(var) > 0 || SCIPvarGetNLocksUp(var) > 1 ) 
      {
         applicable = FALSE;
         break;
      }

      negated = FALSE;
      
      /* get the active variable */
      SCIP_CALL( SCIPvarGetProbvarBinary(&var, &negated) );
      assert(SCIPvarIsActive(var));

      if( negated )
         profits[v] = SCIPvarGetObj(var);
      else
         profits[v] = -SCIPvarGetObj(var);

      SCIPdebugMessage("variable <%s> -> item size %"SCIP_LONGINT_FORMAT", profit <%g>\n", 
         SCIPvarGetName(vars[v]), consdata->weights[v], profits[v]);
      items[v] = v;
   }
   
   if( applicable )
   {
      SCIP_Bool success;

      SCIPdebugMessage("the knapsack constraint <%s> is independent to rest of the problem\n", SCIPconsGetName(cons));
      SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );

      /* solve knapsack problem exactly */
      SCIP_CALL( SCIPsolveKnapsackExactly(scip, consdata->nvars, consdata->weights, profits, consdata->capacity, 
            items, solitems, nonsolitems, &nsolitems, &nnonsolitems, &solval, &success) );

      if( !success )
         goto TERMINATE;
      
      /* apply solution of the knapsack as dual reductions */
      for( v = 0; v < nsolitems; ++v )
      {
         SCIP_VAR* var;
         
         var = vars[solitems[v]];
         assert(var != NULL);
         
         SCIPdebugMessage("variable <%s> only locked up in knapsack constraints: dual presolve <%s>[%.15g,%.15g] >= 1.0\n",
            SCIPvarGetName(var), SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
         SCIP_CALL( SCIPtightenVarLb(scip, var, 1.0, TRUE, &infeasible, &tightened) );
         assert(!infeasible);
         assert(tightened);
         (*nfixedvars)++;
      }
      
      for( v = 0; v < nnonsolitems; ++v )
      {
         SCIP_VAR* var;
         
         var = vars[nonsolitems[v]];
         assert(var != NULL);

         SCIPdebugMessage("variable <%s> has no down locks: dual presolve <%s>[%.15g,%.15g] <= 0.0\n",
            SCIPvarGetName(var), SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
         SCIP_CALL( SCIPtightenVarUb(scip, var, 0.0, TRUE, &infeasible, &tightened) );
         assert(!infeasible);
         assert(tightened);
         (*nfixedvars)++;
      }

      SCIP_CALL( SCIPdelCons(scip, cons) );
      (*ndelconss)++;
      (*deleted) = TRUE;
   }
   
 TERMINATE:
   SCIPfreeBufferArray(scip, &nonsolitems);
   SCIPfreeBufferArray(scip, &solitems);
   SCIPfreeBufferArray(scip, &items);
   SCIPfreeBufferArray(scip, &profits);

   return SCIP_OKAY;
}

/** sort the variables and weights w.r.t. the clique partition; thereby ensure the current order of the variables when a
 *  weight of one variable is greater or equal another weight and both variables are in the same cliques */
static
SCIP_RETCODE stableSort(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< knapsack constraint data */
   SCIP_VAR**            vars,               /**< array for sorted variables */
   SCIP_Longint*         weights,            /**< array for sorted weights */
   int*                  cliquestartposs,    /**< starting position array for each clique */
   SCIP_Bool             usenegatedclique    /**< should negated or normal clique partition be used */
   )
{ 
   SCIP_VAR** origvars;
   int norigvars;
   SCIP_Longint* origweights;
   int* cliquepartition;
   int ncliques;

   SCIP_VAR*** varpointers;
   SCIP_Longint** weightpointers;
   int* cliquecount;

   int nextpos;
   int c;
   int v;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(vars != NULL);
   assert(weights != NULL);
   assert(cliquestartposs != NULL);

   origweights = consdata->weights;
   origvars = consdata->vars;
   norigvars = consdata->nvars;

   assert(origvars != NULL || norigvars == 0);
   assert(origweights != NULL || norigvars == 0);

   if( norigvars == 0 )
      return SCIP_OKAY;

   if( usenegatedclique )
   {
      assert(consdata->negcliquepartitioned);
         
      cliquepartition = consdata->negcliquepartition;
      ncliques = consdata->nnegcliques;
   }
   else
   {
      assert(consdata->cliquepartitioned);

      cliquepartition = consdata->cliquepartition;
      ncliques = consdata->ncliques;
   }

   assert(cliquepartition != NULL);
   assert(ncliques > 0);

   /* we first count all clique items and alloc temporary memory for a bucket sort */
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquecount, ncliques) );
   BMSclearMemoryArray(cliquecount, ncliques);
      
   /* first we count for each clique the number of elements */
   for( v = norigvars - 1; v >= 0; --v )
   {
      assert(0 <= cliquepartition[v] && cliquepartition[v] < ncliques);
      ++(cliquecount[cliquepartition[v]]);
   }

   /*@todo: maybe it is better to put largest cliques up front */

#ifndef NDEBUG
   BMSclearMemoryArray(vars, norigvars);
   BMSclearMemoryArray(weights, norigvars);
#endif
   SCIP_CALL( SCIPallocBufferArray(scip, &varpointers, ncliques) );
   SCIP_CALL( SCIPallocBufferArray(scip, &weightpointers, ncliques) );
   
   nextpos = 0;
   /* now we initialize all start pointers for each clique, so they will be ordered */
   for( c = 0; c < ncliques; ++c )
   {
      /* to reach the goal that all variables of each clique will be standing next to each other we will initialize the
       * starting pointers for each clique by adding the number of each clique to the last clique starting pointer
       * e.g. clique1 has 4 elements and clique2 has 3 elements the the starting pointer for clique1 will be the pointer
       *      to vars[0], the starting pointer to clique2 will be the pointer to vars[4] and to clique3 it will be 
       *      vars[7]
       *
       */
      varpointers[c] = (SCIP_VAR**) (vars + nextpos);
      cliquestartposs[c] = nextpos;
      weightpointers[c] = (SCIP_Longint*) (weights + nextpos);
      assert(cliquecount[c] > 0);
      nextpos += cliquecount[c];
      assert(nextpos > 0);
   }
   assert(nextpos == norigvars);
   cliquestartposs[c] = nextpos;

   /* now we copy all variable and weights to the right order */
   for( v = 0; v < norigvars; ++v )
   {
      *(varpointers[cliquepartition[v]]) = origvars[v];
      ++(varpointers[cliquepartition[v]]);
      *(weightpointers[cliquepartition[v]]) = origweights[v];
      ++(weightpointers[cliquepartition[v]]);
   }
#ifndef NDEBUG
   for( v = 0; v < norigvars; ++v )
   {
      assert(vars[v] != NULL);
      assert(weights[v] > 0);
   }
#endif

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &weightpointers);        
   SCIPfreeBufferArray(scip, &varpointers);
   SCIPfreeBufferArray(scip, &cliquecount);

   return SCIP_OKAY;
}

/** propagation method for knapsack constraints */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   SCIP_Bool*            redundant,          /**< pointer to store whether constraint is redundant */
   int*                  nfixedvars,         /**< pointer to count number of fixings */
   SCIP_Bool             usenegatedclique    /**< should negated clique information be used */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;
   SCIP_Longint zerosweightsum;
   SCIP_Longint* secondmaxweights;
   SCIP_Longint minweightsum;

   int nvars;
   int i;
   int nnegcliques;

   SCIP_VAR** myvars;
   SCIP_Longint* myweights;
   int* cliquestartposs;
   SCIP_Longint localminweightsum;
   SCIP_Bool foundmax;
   int c;

   assert(scip != NULL);
   assert(cons != NULL);
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

   /* we need a merged constraint cause without it the negated clique information could be invalid */
   usenegatedclique = usenegatedclique && consdata->merged;
   secondmaxweights = NULL;
   nnegcliques = -1;

   do
   {
      zerosweightsum = 0;
      nvars = consdata->nvars;
      
      minweightsum = 0;
      localminweightsum = 0;

      /* make sure, the items are sorted by non-increasing weight */
      sortItems(consdata);
      
      /* (1) compute the minimum weight of the knapsack constraint using negated clique information;
       *     a negated clique means, that at most one of the clique variables can be zero
       *     - minweightsum = sum_{negated cliques C} ( sum(wi : i \in C) - W(C) ), where W(C) is the maximal weight of C
       *     
       *     if for i \in C (a negated clique) minweightsum - wi + W(C) > capacity => xi = 1 
       *     since replacing i with the element of maximal weight leads to infeasibility
       */
      if( usenegatedclique && nvars > 0 )
      {
         /* compute clique partitions */
         SCIP_CALL( calcCliquepartition(scip, consdata, FALSE, TRUE) );
         nnegcliques = consdata->nnegcliques;

         /* if we have no real negated cliques we can stop here */
         if( nnegcliques == nvars )
            break;

         /* allocate temporary memory and initialize it */
         SCIP_CALL( SCIPduplicateBufferArray(scip, &myvars, consdata->vars, nvars) );
         SCIP_CALL( SCIPduplicateBufferArray(scip, &myweights, consdata->weights, nvars) ) ;  
         SCIP_CALL( SCIPallocBufferArray(scip, &cliquestartposs, nnegcliques + 1) );
         SCIP_CALL( SCIPallocBufferArray(scip, &secondmaxweights, nnegcliques) );
         BMSclearMemoryArray(secondmaxweights, nnegcliques);

         /* resort variables to avoid quadratic algorithm later on */
         SCIP_CALL( stableSort(scip, consdata, myvars, myweights, cliquestartposs, TRUE) );

         c = 0;
         foundmax = FALSE;
         
         for( i = 0; i < nvars; ++i )
         {
            /* ignore variables of the negated clique which are fixed to one since these are counted in
             * consdata->onesweightsum 
             */

            /* if there are only one variable negated cliques left we can stop */
            if( nnegcliques - c == nvars - i )
            {
               minweightsum += localminweightsum;
               localminweightsum = 0;
               break;
            }

            /* for summing up the minimum active weights due to cliques we have to omit the biggest weights of each
             * clique, we can only skip this clique if this variables is not fixed to zero, otherwise we have to fix all
             * other clique variables to one 
             */
            if( cliquestartposs[c] == i )
            {
	       assert(myweights[i] > 0);
               ++c;
               minweightsum += localminweightsum;
               localminweightsum = 0;
               foundmax = TRUE;

               if( SCIPvarGetLbLocal(myvars[i]) > 0.5 )
                  foundmax = FALSE;
               
               if( SCIPvarGetUbLocal(myvars[i]) > 0.5 )
                  continue;
            }

            if( SCIPvarGetLbLocal(myvars[i]) < 0.5 )
            {
	       assert(myweights[i] > 0);

               if( SCIPvarGetUbLocal(myvars[i]) > 0.5 )
               {
                  assert(myweights[i] <= myweights[cliquestartposs[c - 1]]); 

                  if( !foundmax )
                  {
                     foundmax = TRUE;

                     /* overwrite cliquestartpos to the position of the first unfixed variable in this clique */
                     cliquestartposs[c - 1] = i;

                     continue;
                  }
                  /* memorize second max weight for each clique */
                  if( secondmaxweights[c - 1] == 0 )
                     secondmaxweights[c - 1] = myweights[i];

                  localminweightsum += myweights[i];
               }
               /* we found a fixed variable to zero so all other variables in this negated clique have to be fixed to one */
               else
               {
		  int v;
		  /* fix all other variables of the negated clique to 1 */
		  for( v = cliquestartposs[c - 1]; v < cliquestartposs[c]; ++v )
		  {
                     if( v != i && SCIPvarGetLbLocal(myvars[v]) < 0.5 )
                     {
                        SCIPdebugMessage(" -> fixing variable <%s> to 1, due to negated clique information\n", SCIPvarGetName(myvars[v]));
                        SCIP_CALL( SCIPinferBinvarCons(scip, myvars[v], TRUE, cons, SCIPvarGetIndex(myvars[i]), &infeasible, &tightened) );
                        if( infeasible )
                        {
                           assert( SCIPvarGetUbLocal(myvars[v]) < 0.5 );
                          
                           /* initialize the conflict analysis */
                           SCIP_CALL( SCIPinitConflictAnalysis(scip) );

                           /* add the two variables which are fixed to zero within a negated clique */
                           SCIP_CALL( SCIPaddConflictBinvar(scip, myvars[i]) );
                           SCIP_CALL( SCIPaddConflictBinvar(scip, myvars[v]) );

                           /* start the conflict analysis */
                           SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

                           *cutoff = TRUE;
                           break;
                        }	 
                        assert(tightened);
                        ++(*nfixedvars);
                        SCIP_CALL( SCIPresetConsAge(scip, cons) );
                     }
		  }
		  if( *cutoff )
		     break;

                  /* reset local minweightsum for clique because all fixed to one variables are now counted in consdata->onesweightsum */
                  localminweightsum = 0;
                  /* we can jump to the end of this clique */
                  i = cliquestartposs[c] - 1;
               }
            }
         }
         /* add last clique minweightsum */
         minweightsum += localminweightsum;

         SCIPdebugMessage("knapsack constraint <%s> has minimum weigth sum of <%"SCIP_LONGINT_FORMAT">\n", 
            SCIPconsGetName(cons), minweightsum + consdata->onesweightsum );

         /* check, if weights of fixed variables don't exceeds knapsack capacity */
         if( !(*cutoff) && consdata->capacity >= minweightsum + consdata->onesweightsum )
         {
            SCIP_Longint maxweight;

            c = 0;
            maxweight = -1;

            /* check for each negated clique if we can conclude a fixing of a variable to one */
            for( i = 0; i < nvars; ++i )
            {
               /* if there are only one variable negated cliques left we can stop */
               if( nnegcliques - c == nvars - i )
                  break;

               /* we cannot fix the biggest weight */
               if( cliquestartposs[c] == i )
               {
                  maxweight = myweights[i];
                  ++c;
                  continue;
               }

               /* only check variables of negated cliques for which no variable is locally fixed */
               if( SCIPvarGetLbLocal(myvars[i]) < 0.5 && SCIPvarGetUbLocal(myvars[i]) > 0.5  )
               {		 
                  assert(maxweight >= myweights[i]);

                  /* if for i \in C (a negated clique) minweightsum - wi + W(C) > capacity => xi = 1 
                   * since replacing i with the element of maximal weight leads to infeasibility */
                  if( consdata->onesweightsum + minweightsum - myweights[i] + maxweight > consdata->capacity  )
                  {
                     SCIPdebugMessage(" -> fixing variable <%s> to 1, due to negated clique information\n", SCIPvarGetName(myvars[i]));
                     SCIP_CALL( SCIPinferBinvarCons(scip, myvars[i], TRUE, cons, -i, &infeasible, &tightened) );
                     assert(!infeasible);
                     assert(tightened);
                     ++(*nfixedvars);
                     SCIP_CALL( SCIPresetConsAge(scip, cons) );

                     /* update minweightsum because now the variable is fixed to one and its weight is counted by
                      * consdata->onesweightsum 
                      */
                     minweightsum -= myweights[i];
                     assert(minweightsum >= 0);
                  }
               }
            }
         }

	 if( *cutoff )
         {
            SCIPfreeBufferArray(scip, &secondmaxweights);
            SCIPfreeBufferArray(scip, &cliquestartposs);
            SCIPfreeBufferArray(scip, &myweights);
            SCIPfreeBufferArray(scip, &myvars);

            return SCIP_OKAY;
         }
      }

      /* check, if weights of fixed variables already exceed knapsack capacity */
      if( consdata->capacity < minweightsum + consdata->onesweightsum )
      {
         SCIPdebugMessage(" -> cutoff - fixed weight: %"SCIP_LONGINT_FORMAT", capacity: %"SCIP_LONGINT_FORMAT", minimum weight sum: %"SCIP_LONGINT_FORMAT" \n", 
            consdata->onesweightsum, consdata->capacity, minweightsum);
            
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         *cutoff = TRUE;

         if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
         {
            /* start conflict analysis with the fixed-to-one variables, add only as many as need to exceed the capacity */
            SCIP_Longint weight;

            weight = minweightsum;

            SCIP_CALL( SCIPinitConflictAnalysis(scip) );
            for( i = 0; i < nvars && weight <= consdata->capacity; i++ )
            {
               if( SCIPvarGetLbLocal(consdata->vars[i]) > 0.5)
               {
                  SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[i]) );
                  weight += consdata->weights[i];
               }
            }
         
            SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
         }

         if( usenegatedclique && nvars > 0 )
         {
            SCIPfreeBufferArray(scip, &secondmaxweights);
            SCIPfreeBufferArray(scip, &cliquestartposs);
            SCIPfreeBufferArray(scip, &myweights);
            SCIPfreeBufferArray(scip, &myvars);
         }
         return SCIP_OKAY;
      }

      assert(consdata->negcliquepartitioned || minweightsum == 0);
      
      /* if the sum of all weights of fixed variables to one plus the minimalweightsum (minimal weight which is already
       * used in this knapsack due to negated cliques) plus any weight minus the second largest weight in this cliques
       * exceeds the capacity the variables have to be fixed to zero (these variables should only be variables in the
       * cliques which have maxweights)
       */
      if( usenegatedclique && nvars > 0 )
      {
         SCIP_VAR* var;

         assert(nnegcliques == consdata->nnegcliques);
         assert(myvars != NULL);
         assert(myweights != NULL);
         assert(secondmaxweights != NULL);
         assert(cliquestartposs != NULL);

         for( c = 0; c < nnegcliques; ++c )
         {
            if( consdata->onesweightsum + minweightsum + myweights[cliquestartposs[c]] - secondmaxweights[c] > consdata->capacity )
            {
               assert(myweights[cliquestartposs[c]] >= secondmaxweights[c]);
                              
               var = myvars[cliquestartposs[c]];
               if( SCIPvarGetLbLocal(var) < 0.5 && SCIPvarGetUbLocal(var) > 0.5 )
               {
                  SCIPdebugMessage(" -> fixing variable <%s> to 0\n", SCIPvarGetName(var));
                  SCIP_CALL( SCIPresetConsAge(scip, cons) );
                  SCIP_CALL( SCIPinferBinvarCons(scip, var, FALSE, cons, cliquestartposs[c], &infeasible, &tightened) );
                  assert(!infeasible);
                  assert(tightened);
                  (*nfixedvars)++;
               }
            }
         }

         SCIPfreeBufferArray(scip, &secondmaxweights);
         SCIPfreeBufferArray(scip, &cliquestartposs);
         SCIPfreeBufferArray(scip, &myweights);
         SCIPfreeBufferArray(scip, &myvars);
      }
   }
   while( FALSE );

   /* check, if weights of fixed variables already exceed knapsack capacity, this can only happen if 'usenegatedclique'
    * is FALSE, or 'nnegcliques == nvars', otherwise the stronger condition above should have led to a cutoff 
    */
   if( consdata->capacity < consdata->onesweightsum )
   {
      SCIPdebugMessage(" -> cutoff - fixed weight: %"SCIP_LONGINT_FORMAT", capacity: %"SCIP_LONGINT_FORMAT" \n", 
         consdata->onesweightsum, consdata->capacity);
            
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      *cutoff = TRUE;

      if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      {
         /* start conflict analysis with the fixed-to-one variables, add only as many as need to exceed the capacity */
         SCIP_Longint weight;

         weight = 0;

         SCIP_CALL( SCIPinitConflictAnalysis(scip) );
         for( i = 0; i < nvars && weight <= consdata->capacity; i++ )
         {
            if( SCIPvarGetLbLocal(consdata->vars[i]) > 0.5)
            {
               SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[i]) );
               weight += consdata->weights[i];
            }
         }
         
         SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
      }

      return SCIP_OKAY;
   }
      
   /* fix all variables to zero, that don't fit into the knapsack anymore */
   for( i = 0; i < nvars; ++i )
   {
      /* if all weights of fixed variables to one plus any weight exceeds the capacity the variables have to be fixed
       * to zero
       */
      if( SCIPvarGetLbLocal(consdata->vars[i]) < 0.5 )
      {
         if( SCIPvarGetUbLocal(consdata->vars[i]) > 0.5 )
         {
            if( consdata->onesweightsum + consdata->weights[i] > consdata->capacity )
            {
               SCIPdebugMessage(" -> fixing variable <%s> to 0\n", SCIPvarGetName(consdata->vars[i]));
               SCIP_CALL( SCIPresetConsAge(scip, cons) );
               SCIP_CALL( SCIPinferBinvarCons(scip, consdata->vars[i], FALSE, cons, i, &infeasible, &tightened) );
               assert(!infeasible);
               assert(tightened);
               (*nfixedvars)++;
               zerosweightsum += consdata->weights[i];
            }
         }
         else
            zerosweightsum += consdata->weights[i];
      }
   }

   assert(consdata->onesweightsum + zerosweightsum <= consdata->weightsum);

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

/*  tries to simplify weights and delete variables in knapsack a^Tx <= rhs in case there is only one binary variable
 *  with an odd weight and the capacity is odd too, then:
 *    - if the odd weight is equal to 1, delete the variable and decrease capacity by 1
 *    - otherwise, decrease the weight and the capacity by 1
 *  Afterwards we use the normalize method to further simplify the inequality. 
 *  @todo if all variables with odd weights are in one (not negated!) clique we can decrease all weights and the capacity
 *        by 1 but we need to enforce this clique information (because we maybe loose this information by reducing 
 *        the weights) by creating a setppc constraint, but this seems to be expensive?!
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

   SCIP_CALL( mergeMultiples(scip, cons, cutoff) );
   assert(consdata->merged);
   if( *cutoff )
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
            if( noddvals >= 2 )
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

   /* check infeasibility */
   if ( consdata->onesweightsum > consdata->capacity )
   {
      SCIPdebugMessage("apply fixings detected cutoff.\n");
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   v = 0;
   while( v < consdata->nvars )
   {
      SCIP_VAR* var;

      var = consdata->vars[v];
      assert(SCIPvarIsBinary(var));

      if( SCIPvarGetLbGlobal(var) > 0.5 )
      {
         assert(SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(var), 1.0));
         consdata->capacity -= consdata->weights[v];
         SCIP_CALL( delCoefPos(scip, cons, v) );
         consdata->cliquesadded = FALSE; /* reduced capacity might lead to larger cliques */
      }
      else if( SCIPvarGetUbGlobal(var) < 0.5 )
      {
         assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(var), 0.0));
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

	       if( !SCIPvarIsBinary(aggrvars[i]) )
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
      SCIP_CALL( mergeMultiples(scip, cons, cutoff) );
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
   int nzeros;

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
   assert(0 <= probindex && probindex < SCIPgetNVars(scip) - SCIPgetNContVars(scip));
   assert(memlimitreached != NULL);

   nzeros = *nzeroitems;

   /* allocate enough memory */
   if( nzeros == *zeroitemssize )
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
   assert(nzeros < *zeroitemssize);

   if( *memlimitreached )
      *memlimitreached = FALSE;

   /* insert element */
   (*zeroitems)[nzeros] = knapsackidx;
   (*nextidxs)[nzeros] = firstidxs[value][probindex];
   if( firstidxs[value][probindex] == 0 )
   {
      liftcands[value][nliftcands[value]] = probindex;
      ++nliftcands[value];
   }
   firstidxs[value][probindex] = nzeros;
   ++(*nzeroitems);
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
   int*                  nchgcoefs,          /**< pointer to count total number of changed coefficients */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
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
   int nvars;
   int cliquenum;
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

   nvars = consdata->nvars;
   
   /* check if the knapsack has too many items/cliques for applying this costly method */
   if( (!consdata->cliquepartitioned && nvars > MAX_USECLIQUES_SIZE) || consdata->ncliques > MAX_USECLIQUES_SIZE ) 
      return SCIP_OKAY;
   
   /* sort items, s.t. the heaviest one is in the first position */
   sortItems(consdata);

   if( !consdata->cliquepartitioned && nvars > MAX_USECLIQUES_SIZE )
      return SCIP_OKAY;
   
   /* we have to consider all integral variables since even integer and implicit integer variables can have binary bounds */
   nbinvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);
   assert(nbinvars > 0);
   binvars = SCIPgetVars(scip);

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
      /* first add the implications (~x == 1 -> x == 0) */
      {
         SCIP_Bool implvalue;
         int probindex;

         probindex = SCIPvarGetProbindex(var);
         assert(0 <= probindex && probindex < nbinvars);
         
	 implvalue = !value;
            
         /* insert the item into the list of the implied variable/value */
         assert( !zeroiteminserted[implvalue][probindex] );
            
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

         assert(SCIPvarIsBinary(implvars[j]));
         probindex = SCIPvarGetProbindex(implvars[j]);
         assert(probindex < nbinvars);

         /* this assert should hold, but if not there is a old continue later on */
         assert(probindex >= 0);

         /* consider only implications with active implvar */
         if( probindex < 0 )
            continue;

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
   SCIP_CALL( calcCliquepartition(scip, consdata, TRUE, FALSE) );

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

   /* calculates maximal weight of cliques */
   for( i = 0; i < consdata->nvars; ++i )
   {
      cliquenum = consdata->cliquepartition[i];
      assert(0 <= cliquenum && cliquenum < consdata->nvars);

      if( !cliqueused[cliquenum] )
      {
         maxcliqueweightsum += consdata->weights[i];
         cliqueused[cliquenum] = TRUE;
         tmpindices[tmp] = cliquenum;
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
            cliquenum = consdata->cliquepartition[j];
            assert(0 <= cliquenum && cliquenum < consdata->nvars);
	    if( !itemremoved[j] )
            {
               if( !cliqueused[cliquenum] )
               {
                  cliqueweightsum += consdata->weights[j];
                  cliqueused[cliquenum] = TRUE;
                  tmpindices[tmp] = cliquenum;
                  ++tmp;
               }
		
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

   if( naddvars > 0 )
   {
      /* if new items were added, multiple entries of the same variable are possible and we have to clean up the constraint */
      SCIP_CALL( mergeMultiples(scip, cons, cutoff) );
   }

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
   int*                  nchgsides,          /**< pointer to count number of side changes */
   int*                  naddconss,          /**< pointer to count number of added constraints */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
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

   SCIP_CALL( mergeMultiples(scip, cons, cutoff) );
   if( *cutoff )
      return SCIP_OKAY;

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
   if( consdata->nvars <= MAX_USECLIQUES_SIZE || (consdata->cliquepartitioned && consdata->ncliques <= MAX_USECLIQUES_SIZE) ) 
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
         int cliquenum;
         SCIP_Bool zeroweights;

         assert(consdata->merged);

         /* sort items, s.t. the heaviest one is in the first position */
         sortItems(consdata);

         /* calculate a clique partition */
         SCIP_CALL( calcCliquepartition(scip, consdata, TRUE, FALSE) );

         /* if there are only single element cliques, rule (2) is equivalent to rule (1) */
         if( consdata->cliquepartition[consdata->nvars - 1] == consdata->nvars - 1 )
            break;

         /* calculate the maximal weight of the cliques and store the clique type */
         cliqueweightsum = 0;
         ncliques = 0;

         for( i = 0; i < consdata->nvars; ++i )
         {
            SCIP_Longint weight;

            cliquenum = consdata->cliquepartition[i];
            assert(0 <= cliquenum && cliquenum <= ncliques);

            weight = consdata->weights[i];
            assert(weight > 0);

            if( cliquenum == ncliques )
            {
               maxcliqueweights[ncliques] = weight;
               cliqueweightsum += weight;
               ++ncliques;
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
               SCIP_Longint newminweightsuminclique;
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

                     /* cache the new weight */
                     assert(nnewweights < consdata->nvars);
                     newweightvals[nnewweights] = newweight;
                     newweightidxs[nnewweights] = j;
                     nnewweights++;

                     assert(newweight <= newmincliqueweight); /* items are sorted by non-increasing weight! */
                     newmincliqueweight = newweight;
                  }
               }

               /* check if our clique information results out of this knapsack constraint and if so check if we would loose the clique information */
               if( nnewweights > 1 )
               {
#ifndef NDEBUG
                  j = newweightidxs[nnewweights - 2];
                  assert(0 <= j && j < consdata->nvars);
                  assert(consdata->cliquepartition[j] == i);
                  j = newweightidxs[nnewweights - 1];
                  assert(0 <= j && j < consdata->nvars);
                  assert(consdata->cliquepartition[j] == i);
#endif     

                  newminweightsuminclique = newweightvals[nnewweights - 2];
                  newminweightsuminclique += newweightvals[nnewweights - 1];

                  /* check if these new two minimal weights both fit into the knapsack; 
                   * if this is true, we have to add a clique constraint in order to enforce the clique
                   * (otherwise, the knapsack might have been one of the reasons for the clique, and the weight
                   * reduction might be infeasible, i.e., allows additional solutions)
                   */
                  if( newminweightsuminclique <= newcapacity )
                     forceclique = TRUE;
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
		     (*naddconss)++;
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
   SCIP_CALL( tightenWeightsLift(scip, cons, nchgcoefs, cutoff) );

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

/** adds negated cliques of the knapsack constraint to the global clique table */
static
SCIP_RETCODE addNegatedCliques(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< knapsack constraint */
   SCIP_Bool*const       cutoff,             /**< pointer to store whether the node can be cut off */
   int*const             nbdchgs             /**< pointer to count the number of performed bound changes */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** poscliquevars;
   SCIP_VAR** cliquevars;
   SCIP_Longint* maxweights;
   SCIP_Longint* gainweights;
   int* gaincliquepartition;
   SCIP_Bool* cliqueused;
   SCIP_Longint minactduetonegcliques;
   SCIP_Longint freecapacity;
   SCIP_Longint lastweight;
   SCIP_Longint beforelastweight;
   int nposcliquevars;
   int ncliquevars;
   int nvars;
   int nnegcliques;
   int lastcliqueused;
   int thisnbdchgs;
   int v;
   int w;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nbdchgs != NULL);

   *cutoff = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   nvars = consdata->nvars;

   /* check whether the cliques have already been added */
   if( consdata->cliquesadded || nvars == 0 )
      return SCIP_OKAY;

   /* make sure, the items are merged */
   SCIP_CALL( mergeMultiples(scip, cons, cutoff) );
   if( *cutoff )
      return SCIP_OKAY;

   /* make sure, items are sorted by non-increasing weight */
   sortItems(consdata);

   assert(consdata->merged);

   /* calculate a clique partition */
   SCIP_CALL( calcCliquepartition(scip, consdata, FALSE, TRUE) );
   nnegcliques = consdata->nnegcliques;

   /* if we have no negated cliques, stop */
   if( nnegcliques == nvars )
      return SCIP_OKAY;
   
   /* get temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &poscliquevars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquevars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &gainweights, nvars) );
   BMSclearMemoryArray(gainweights, nvars);
   SCIP_CALL( SCIPallocBufferArray(scip, &gaincliquepartition, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxweights, nnegcliques) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cliqueused, nnegcliques) );
   BMSclearMemoryArray(cliqueused, nnegcliques);

   nnegcliques = 0;
   minactduetonegcliques = 0;

   /* determine maximal weights for all negated cliques and calculate minimal weightsum due to negated cliques */
   for( v = 0; v < nvars; ++v )
   {
      assert(0 <= consdata->negcliquepartition[v] && consdata->negcliquepartition[v] <= nnegcliques);
      assert(consdata->weights[v] > 0);

      if( consdata->negcliquepartition[v] == nnegcliques )
      {
         nnegcliques++;
         maxweights[consdata->negcliquepartition[v]] = consdata->weights[v];
      }
      else
         minactduetonegcliques += consdata->weights[v];
   }

   nposcliquevars = 0;
   
   /* add cliques, using negated cliques information */
   if( minactduetonegcliques > 0 )
   {
      /* free capacity is the rest of not used capacity if the smallest amount of weights due to negated cliques are used */
      freecapacity = consdata->capacity - minactduetonegcliques;

      SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
      SCIPdebugMessage("Try to add negated cliques in knapsack constraint handler for constraint %s; capacity = %"SCIP_LONGINT_FORMAT", minactivity(due to neg. cliques) = %"SCIP_LONGINT_FORMAT", freecapacity = %"SCIP_LONGINT_FORMAT".\n", 
         SCIPconsGetName(cons), consdata->capacity, minactduetonegcliques, freecapacity);

      /* calculate possible gain by switching choosen items in negated cliques */
      for( v = 0; v < nvars; ++v )
      {
         if( !cliqueused[consdata->negcliquepartition[v]] )
         {
            cliqueused[consdata->negcliquepartition[v]] = TRUE;
            for( w = v + 1; w < nvars; ++w )
            {
               /* if we would take the biggest weight instead of another what would we gain, take weight[i] instead of
                * weight[j] (which are both in a negated clique) */
               if( consdata->negcliquepartition[v] == consdata->negcliquepartition[w] 
                  && consdata->weights[v] > consdata->weights[w] )
               {
                  poscliquevars[nposcliquevars] = consdata->vars[w];
                  gainweights[nposcliquevars] = maxweights[consdata->negcliquepartition[v]] - consdata->weights[w];
                  gaincliquepartition[nposcliquevars] = consdata->negcliquepartition[v];
                  ++nposcliquevars;
               }
            }
         }
      }
      
      /* try to create negated cliques */
      if( nposcliquevars > 0 )
      {
         /* sort possible gain per substitution of the clique members */
         SCIPsortDownLongPtrInt(gainweights,(void**) poscliquevars, gaincliquepartition, nposcliquevars);
      
         for( v = 0; v < nposcliquevars; ++v )
         {
            SCIP_CALL( SCIPgetNegatedVar(scip, poscliquevars[v], &cliquevars[0]) );
            ncliquevars = 1;
            lastweight = gainweights[v];
            beforelastweight = -1;
            lastcliqueused = gaincliquepartition[v];
            /* clear cliqueused to get an unused array */
            BMSclearMemoryArray(cliqueused, nnegcliques);
            cliqueused[gaincliquepartition[v]] = TRUE;
            
            /* taking bigger weights make the knapsack redundant so we will create cliques, only take items which are not
             * in the same negated clique and by taking two of them would exceed the free capacity */
            for( w = v + 1; w < nposcliquevars && !cliqueused[gaincliquepartition[w]] && gainweights[w] + lastweight > freecapacity; ++w )
            {
               beforelastweight = lastweight;
               lastweight = gainweights[w];
               lastcliqueused = gaincliquepartition[w];
               cliqueused[gaincliquepartition[w]] = TRUE;
               SCIP_CALL( SCIPgetNegatedVar(scip, poscliquevars[w], &cliquevars[ncliquevars]) );
               ++ncliquevars;
            }
            
            if( ncliquevars > 1 )
            {
#ifdef SCIP_DEBUG
               int b;
               SCIPdebugMessage("adding new Clique: ");
               for( b = 0; b < ncliquevars; ++b )
                  SCIPdebugPrintf("%s ", SCIPvarGetName(cliquevars[b]));
               SCIPdebugPrintf("\n");
#endif
               
               assert(beforelastweight > 0);
               /* add the clique to the clique table */
               SCIP_CALL( SCIPaddClique(scip, cliquevars, NULL, ncliquevars, cutoff, &thisnbdchgs) );
               if( *cutoff )
                  goto TERMINATE;
               *nbdchgs += thisnbdchgs;
               
               /* reset last used clique to get slightly different cliques */
               cliqueused[lastcliqueused] = FALSE;
               
               /* try to replace the last item in the clique by a different item to obtain a slightly different clique */
               for( ++w; w < nposcliquevars && !cliqueused[gaincliquepartition[w]] && beforelastweight + gainweights[w] > freecapacity; ++w )
               {
                  SCIP_CALL( SCIPgetNegatedVar(scip, poscliquevars[w], &cliquevars[ncliquevars - 1]) );
#ifdef SCIP_DEBUG
                  {
                     SCIPdebugMessage("adding new Clique: ");
                     for( b = 0; b < ncliquevars; ++b )
                        SCIPdebugPrintf("%s ", SCIPvarGetName(cliquevars[b]));
                     SCIPdebugPrintf("\n");
                  }
#endif
                  SCIP_CALL( SCIPaddClique(scip, cliquevars, NULL, ncliquevars, cutoff, &thisnbdchgs) );
                  if( *cutoff )
                     goto TERMINATE;
                  *nbdchgs += thisnbdchgs;
               }
            }
         }
      }   
   }

 TERMINATE:
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &cliqueused);
   SCIPfreeBufferArray(scip, &gaincliquepartition);
   SCIPfreeBufferArray(scip, &maxweights);
   SCIPfreeBufferArray(scip, &gainweights);
   SCIPfreeBufferArray(scip, &cliquevars);
   SCIPfreeBufferArray(scip, &poscliquevars);

   return SCIP_OKAY;
}

/** adds cliques of the knapsack constraint to the global clique table */
static
SCIP_RETCODE addCliques(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< knapsack constraint */
   SCIP_Bool*const       cutoff,             /**< pointer to store whether the node can be cut off */
   int*const             nbdchgs             /**< pointer to count the number of performed bound changes */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** cliquevars;
   int ncliquevars;
   int i;
   SCIP_Longint minactduetonegcliques;
   SCIP_Longint freecapacity;
   SCIP_Longint lastweight;
   SCIP_Longint beforelastweight;
   int nnegcliques;
   int cliquenum;
   int thisnbdchgs;
   SCIP_VAR** poscliquevars;
   SCIP_Longint* gainweights;
   int nposcliquevars;
   SCIP_Longint* secondmaxweights;
   int nvars;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nbdchgs != NULL);

   *cutoff = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   nvars = consdata->nvars;

   /* check whether the cliques have already been added */
   if( consdata->cliquesadded || nvars == 0 )
      return SCIP_OKAY;

   /* make sure, the items are merged */
   SCIP_CALL( mergeMultiples(scip, cons, cutoff) );
   if( *cutoff )
      return SCIP_OKAY;

   /* make sure, the items are sorted by non-increasing weight */
   sortItems(consdata);

   assert(consdata->merged);

   /* calculate a clique partition */
   SCIP_CALL( calcCliquepartition(scip, consdata, FALSE, TRUE) );
   nnegcliques = consdata->nnegcliques;
   assert(nnegcliques <= nvars);

   /* get temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &poscliquevars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquevars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &gainweights, nvars) );
   BMSclearMemoryArray(gainweights, nvars);
   SCIP_CALL( SCIPallocBufferArray(scip, &secondmaxweights, nnegcliques) );
   BMSclearMemoryArray(secondmaxweights, nnegcliques);

   minactduetonegcliques = 0;
   
   /* calculate minimal activity due to negated cliques, and determine second maximal weight in each clique */
   if( nnegcliques < nvars )
   {
      nnegcliques = 0;
   
      for( i = 0; i < nvars; ++i )
      {
         SCIP_Longint weight;
         
         cliquenum = consdata->negcliquepartition[i];
         assert(0 <= cliquenum && cliquenum <= nnegcliques);
         
         weight = consdata->weights[i];
         assert(weight > 0);
         
         if( cliquenum == nnegcliques )
            nnegcliques++;
         else
         {
            minactduetonegcliques += weight;
            if( secondmaxweights[cliquenum] == 0 )
               secondmaxweights[cliquenum] = weight;
         }
      }
   }

   /* add cliques, using negated cliques information */
   if( minactduetonegcliques > 0 )
   {
      /* free capacity is the rest of not used capacity if the smallest amount of weights due to negated cliques are used */
      freecapacity = consdata->capacity - minactduetonegcliques;

      SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
      SCIPdebugMessage("Try to add cliques in knapsack constraint handler for constraint %s; capacity = %"SCIP_LONGINT_FORMAT", minactivity(due to neg. cliques) = %"SCIP_LONGINT_FORMAT", freecapacity = %"SCIP_LONGINT_FORMAT".\n", 
         SCIPconsGetName(cons), consdata->capacity, minactduetonegcliques, freecapacity);

      /* create negated cliques out of negated cliques, if we do not take the smallest weight of a cliques ... */
      SCIP_CALL( addNegatedCliques(scip, cons, cutoff, nbdchgs ) );

      nposcliquevars = 0;

      for( i = nvars - 1; i >= 0; --i )
      {
         /* if we would take the biggest weight instead of the second biggest */
         cliquenum = consdata->negcliquepartition[i];
         if( consdata->weights[i] > secondmaxweights[cliquenum] )
         {
            poscliquevars[nposcliquevars] = consdata->vars[i];
            gainweights[nposcliquevars] = consdata->weights[i] - secondmaxweights[cliquenum];
            ++nposcliquevars;
         }
      }
      
      if( nposcliquevars > 0 )
      {
         /* sort possible gain per substitution of the clique members */
         SCIPsortDownLongPtr(gainweights,(void**) poscliquevars, nposcliquevars);

         cliquevars[0] = poscliquevars[0];
         ncliquevars = 1;
         lastweight = gainweights[0];
         beforelastweight = 0;
         
         /* taking bigger weights make the knapsack redundant so we will create cliques */
         for( i = 1; i < nposcliquevars && gainweights[i] + lastweight > freecapacity; ++i )
         {
            beforelastweight = lastweight;
            lastweight = gainweights[i];
            cliquevars[ncliquevars] = poscliquevars[i];
            ++ncliquevars;
         }
         
         if( ncliquevars > 1 )
         {
#ifdef SCIP_DEBUG
            int b;
            SCIPdebugMessage("adding new Clique: ");
            for( b = 0; b < ncliquevars; ++b )
               SCIPdebugPrintf("%s ", SCIPvarGetName(cliquevars[b]));
            SCIPdebugPrintf("\n");
#endif

            assert(beforelastweight > 0);
            /* add the clique to the clique table */
            SCIP_CALL( SCIPaddClique(scip, cliquevars, NULL, ncliquevars, cutoff, &thisnbdchgs) );
            if( *cutoff )
               goto TERMINATE;
            *nbdchgs += thisnbdchgs;
            
            /* try to replace the last item in the clique by a different item to obtain a slightly different clique */
            for( ++i; i < nposcliquevars && beforelastweight + gainweights[i] > freecapacity; ++i )
            {
               cliquevars[ncliquevars - 1] = poscliquevars[i];
#ifdef SCIP_DEBUG
               {
                  SCIPdebugMessage("adding new Clique: ");
                  for( b = 0; b < ncliquevars; ++b )
                     SCIPdebugPrintf("%s ", SCIPvarGetName(cliquevars[b]));
                  SCIPdebugPrintf("\n");
               }
#endif
               SCIP_CALL( SCIPaddClique(scip, cliquevars, NULL, ncliquevars, cutoff, &thisnbdchgs) );
               if( *cutoff )
                  goto TERMINATE;
               *nbdchgs += thisnbdchgs;
            }
         }
      }
   }

   /* build a largest clique by using the items with the maximal weights */
   cliquevars[0] = consdata->vars[0];
   for( i = 1; i < nvars && consdata->weights[i-1] + consdata->weights[i] > consdata->capacity; ++i )
      cliquevars[i] = consdata->vars[i];
   
   ncliquevars = i;

   if( ncliquevars >= 2 )
   {
      SCIP_Longint cliqueminweight;
#ifdef SCIP_DEBUG
      int b;
      SCIPdebugMessage("adding new Clique: ");
      for( b = 0; b < ncliquevars; ++b )
         SCIPdebugPrintf("%s ", SCIPvarGetName(cliquevars[b]));
      SCIPdebugPrintf("\n");
#endif

      /* add the clique to the clique table */
      SCIP_CALL( SCIPaddClique(scip, cliquevars, NULL, ncliquevars, cutoff, &thisnbdchgs) );
      if( *cutoff )
         goto TERMINATE;
      *nbdchgs += thisnbdchgs;

      /* try to replace the last item in the clique by a different item to obtain a slightly different clique */
      cliqueminweight = consdata->weights[ncliquevars-2];
      for( i = ncliquevars; i < nvars && cliqueminweight + consdata->weights[i] > consdata->capacity; ++i )
      {
         cliquevars[ncliquevars - 1] = consdata->vars[i];
#ifdef SCIP_DEBUG
         {
            SCIPdebugMessage("adding new Clique: ");
            for( b = 0; b < ncliquevars; ++b )
               SCIPdebugPrintf("%s ", SCIPvarGetName(cliquevars[b]));
            SCIPdebugPrintf("\n");
         }
#endif
         SCIP_CALL( SCIPaddClique(scip, cliquevars, NULL, ncliquevars, cutoff, &thisnbdchgs) );
         if( *cutoff )
            goto TERMINATE;
         *nbdchgs += thisnbdchgs;
      }
   }

 TERMINATE:
   /* free temporary memory and mark the constraint */
   SCIPfreeBufferArray(scip, &secondmaxweights);
   SCIPfreeBufferArray(scip, &gainweights);
   SCIPfreeBufferArray(scip, &cliquevars);
   SCIPfreeBufferArray(scip, &poscliquevars);
   consdata->cliquesadded = TRUE;

   return SCIP_OKAY;
}


/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashGetKeyKnapsackcons)
{  /*lint --e{715}*/
   /* the key is the element itself */ 
   return elem;
}

/** returns TRUE iff both keys are equal; two constraints are equal if they have the same variables and the 
 * same coefficients 
 */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqKnapsackcons)
{
   SCIP* scip;
   SCIP_CONSDATA* consdata1;
   SCIP_CONSDATA* consdata2;
   int i;

   consdata1 = SCIPconsGetData((SCIP_CONS*)key1);
   consdata2 = SCIPconsGetData((SCIP_CONS*)key2);
   assert(consdata1->sorted);
   assert(consdata2->sorted);
   scip = (SCIP*)userptr; 
   assert(scip != NULL);
   
   /* checks trivial case */
   if( consdata1->nvars != consdata2->nvars )
      return FALSE;

   for( i = consdata1->nvars - 1; i >= 0; --i )
   {
      /* tests if variables are equal */
      if( consdata1->vars[i] != consdata2->vars[i] )
      {
         assert(SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == 1 || 
            SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == -1);
         return FALSE;
      }
      assert(SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == 0); 
      
      /* tests if weights are equal too */  
      if( consdata1->weights[i] != consdata2->weights[i] )
         return FALSE;
   } 
   
   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValKnapsackcons)
{
   SCIP* scip;
   SCIP_CONSDATA* consdata;
   unsigned int hashval;
   int minidx;
   int mididx;
   int maxidx;
   int maxabsval;

   consdata = SCIPconsGetData((SCIP_CONS*)key);
   assert(consdata != NULL);
   assert(consdata->nvars > 0);

   scip = (SCIP*)userptr; 
   assert(scip != NULL);

   /* sorts the constraints */
   sortItems(consdata);

   minidx = SCIPvarGetIndex(consdata->vars[0]);
   mididx = SCIPvarGetIndex(consdata->vars[consdata->nvars / 2]);
   maxidx = SCIPvarGetIndex(consdata->vars[consdata->nvars - 1]);
   assert(minidx >= 0 && mididx >= 0 && maxidx >= 0);

   if( consdata->weights[0] > INT_MAX )
      maxabsval = 0;
   else
      maxabsval = (int) consdata->weights[0];

   /* hash value depends on vectors of variable indices */
   hashval = (consdata->nvars << 29) + (minidx << 22) + (mididx << 11) + maxidx + maxabsval; /*lint !e701*/

   return hashval;
}

/** updates the flags of the first constraint according to the ones of the second constraint */
static
SCIP_RETCODE updateFlags(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons0,              /**< constraint that should stay */
   SCIP_CONS*            cons1               /**< constraint that should be deleted */
   )
{
   if( SCIPconsIsInitial(cons1) )
   {
      SCIP_CALL( SCIPsetConsInitial(scip, cons0, TRUE) );
   }
   if( SCIPconsIsSeparated(cons1) )
   {
      SCIP_CALL( SCIPsetConsSeparated(scip, cons0, TRUE) );
   }
   if( SCIPconsIsEnforced(cons1) )
   {
      SCIP_CALL( SCIPsetConsEnforced(scip, cons0, TRUE) );
   }
   if( SCIPconsIsChecked(cons1) )
   {
      SCIP_CALL( SCIPsetConsChecked(scip, cons0, TRUE) );
   }
   if( SCIPconsIsPropagated(cons1) )
   {
      SCIP_CALL( SCIPsetConsPropagated(scip, cons0, TRUE) );
   }
   if( !SCIPconsIsDynamic(cons1) )
   {
      SCIP_CALL( SCIPsetConsDynamic(scip, cons0, FALSE) );
   }
   if( !SCIPconsIsRemovable(cons1) )
   {
      SCIP_CALL( SCIPsetConsRemovable(scip, cons0, FALSE) );
   }
   if( SCIPconsIsStickingAtNode(cons1) )
   {
      SCIP_CALL( SCIPsetConsStickingAtNode(scip, cons0, TRUE) );
   }

   return SCIP_OKAY;
}

/** compares each constraint with all other constraints for possible redundancy and removes or changes constraint 
 *  accordingly; in contrast to preprocessConstraintPairs(), it uses a hash table 
 */
static
SCIP_RETCODE detectRedundantConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   nconss,             /**< number of constraints in constraint set */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the problem is infeasible */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
)
{
   SCIP_HASHTABLE* hashtable;
   int hashtablesize;
   int c;

   assert(scip != NULL);
   assert(blkmem != NULL);
   assert(conss != NULL);
   assert(ndelconss != NULL);

   /* create a hash table for the constraint set */
   hashtablesize = SCIPcalcHashtableSize(10*nconss);
   hashtablesize = MAX(hashtablesize, HASHSIZE_KNAPSACKCONS);
   SCIP_CALL( SCIPhashtableCreate(&hashtable, blkmem, hashtablesize,
         hashGetKeyKnapsackcons, hashKeyEqKnapsackcons, hashKeyValKnapsackcons, (void*) scip) );

   /* check all constraints in the given set for redundancy */
   for( c = nconss - 1; c >= 0; --c )
   {
      SCIP_CONS* cons0;
      SCIP_CONS* cons1;
      SCIP_CONSDATA* consdata0;

      cons0 = conss[c];

      if( !SCIPconsIsActive(cons0) || SCIPconsIsModifiable(cons0) )
         continue;
      
      consdata0 = SCIPconsGetData(cons0);
      assert(consdata0 != NULL);
      if( consdata0->nvars == 0 )
      {
         if( consdata0->capacity < 0 )
         {
            *cutoff = TRUE;
            return SCIP_OKAY;
         }
         else
         {
            SCIP_CALL( SCIPdelCons(scip, cons0) );
            ++(*ndelconss);
            continue;
         }
      }

      /* get constraint from current hash table with same variables and same weights as cons0 */
      cons1 = (SCIP_CONS*)(SCIPhashtableRetrieve(hashtable, (void*)cons0));
 
      if( cons1 != NULL )
      {
         SCIP_CONS* consstay;
         SCIP_CONS* consdel;
         SCIP_CONSDATA* consdata1;

         assert(SCIPconsIsActive(cons1));
         assert(!SCIPconsIsModifiable(cons1));
      
         /* constraint found: create a new constraint with same coeffients and best left and right hand side; 
          * delete old constraints afterwards
          */
         consdata1 = SCIPconsGetData(cons1);
         
         assert(consdata1 != NULL);
         assert(consdata0->nvars > 0 && consdata0->nvars == consdata1->nvars);
         
         assert(consdata0->sorted && consdata1->sorted);
         assert(consdata0->vars[0] == consdata1->vars[0]);
         assert(consdata0->weights[0] == consdata1->weights[0]);

         SCIPdebugMessage("knapsack constraints <%s> and <%s> with equal coefficients\n",
            SCIPconsGetName(cons0), SCIPconsGetName(cons1));

         /* check which constraint has to stay; */
         if( consdata0->capacity < consdata1->capacity )
         {
            consstay = cons0;
            consdel = cons1;
            
            /* exchange consdel with consstay in hashtable */
            SCIP_CALL( SCIPhashtableRemove(hashtable, (void*) consdel) );
            SCIP_CALL( SCIPhashtableInsert(hashtable, (void*) consstay) );
         }
         else
         {
            consstay = cons1; 
            consdel = cons0; 
         }
        
         /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
         SCIP_CALL( updateFlags(scip, consstay, consdel) ); 

         /* delete consdel */
         SCIP_CALL( SCIPdelCons(scip, consdel) );
         ++(*ndelconss);

         assert(SCIPconsIsActive(consstay));
      }
      else
      {
         /* no such constraint in current hash table: insert cons0 into hash table */  
         SCIP_CALL( SCIPhashtableInsert(hashtable, (void*) cons0) );
      }
   }

   /* free hash table */
   SCIPhashtableFree(&hashtable);

   return SCIP_OKAY;
}


/** compares constraint with all prior constraints for possible redundancy or aggregation,
 *  and removes or changes constraint accordingly
 */
static
SCIP_RETCODE preprocessConstraintPairs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   firstchange,        /**< first constraint that changed since last pair preprocessing round */
   int                   chkind,             /**< index of constraint to check against all prior indices upto startind */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONS* cons0;
   SCIP_CONSDATA* consdata0;
   int c;

   assert(scip != NULL);
   assert(conss != NULL);
   assert(firstchange <= chkind);
   assert(ndelconss != NULL);

   /* get the constraint to be checked against all prior constraints */
   cons0 = conss[chkind];
   assert(cons0 != NULL);
   assert(SCIPconsIsActive(cons0));
   assert(!SCIPconsIsModifiable(cons0));

   consdata0 = SCIPconsGetData(cons0);
   assert(consdata0 != NULL);
   assert(consdata0->nvars >= 1);
   assert(consdata0->merged);

   /* sort the constraint */
   sortItems(consdata0);

   /* check constraint against all prior constraints */
   for( c = (consdata0->presolved ? firstchange : 0); c < chkind; ++c )
   {
      SCIP_CONS* cons1;
      SCIP_CONSDATA* consdata1;
      SCIP_Bool iscons0incons1contained;
      SCIP_Bool iscons1incons0contained;
      SCIP_Real quotient;
      int v;
      int v0;
      int v1;

      cons1 = conss[c];
      assert(cons1 != NULL);
      if( !SCIPconsIsActive(cons1) || SCIPconsIsModifiable(cons1) )
         continue;

      consdata1 = SCIPconsGetData(cons1);
      assert(consdata1 != NULL);

      /* if both constraints didn't change since last pair processing, we can ignore the pair */
      if( consdata0->presolved && consdata0->presolved )
         continue;
      
      assert(consdata1->nvars >= 1);
      assert(consdata1->merged);

      /* sort the constraint */
      sortItems(consdata1);

      quotient = ((SCIP_Real) consdata0->capacity) / ((SCIP_Real) consdata1->capacity);
      
      if( consdata0->nvars > consdata1->nvars )
      {
         iscons0incons1contained = FALSE;
         iscons1incons0contained = TRUE;
         v = consdata0->nvars - 1;
      }
      else if( consdata0->nvars < consdata1->nvars )
      {
         iscons0incons1contained = TRUE;
         iscons1incons0contained = FALSE;
         v = consdata1->nvars - 1;
      }
      else
      {
         iscons0incons1contained = TRUE;
         iscons1incons0contained = TRUE;
         v = consdata0->nvars - 1;
      }

      SCIPdebugMessage("preprocess knapsack constraint pair <%s> and <%s>\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));

      /* check consdata0 against consdata1:
       * 1. if all variables var_i of cons1 are in cons0 and for each of these variables (consdata0->weights[i] / quotient) >= consdata1->weights[i] 
       *    cons1 is redundant
       * 2. if all variables var_i of cons0 are in cons1 and for each of these variables (consdata0->weights[i] / quotient) <= consdata1->weights[i] 
       *    cons0 is redundant
       */
      v0 = consdata0->nvars - 1;
      v1 = consdata1->nvars - 1;

      for( ; v >= 0; --v )
      {
         assert(iscons0incons1contained || iscons1incons0contained);
         
         /* now there are more variables in cons1 left */
         if( v1 > v0 )
         {
            iscons1incons0contained = FALSE;
            if( !iscons0incons1contained )
               break;
         }
         /* now there are more variables in cons0 left */
         else if( v1 < v0 )
         {
            iscons0incons1contained = FALSE;
            if( !iscons1incons0contained )
               break;
         }

         assert(v >= v0 && v >= v1);
         
         /* both variables are the same */
         if( consdata0->vars[v0] == consdata1->vars[v1] )
         {
            /* if cons1 is possible contained in cons0 (consdata0->weights[v0] / quotient) must be greater equals consdata1->weights[v1] */
            if( iscons1incons0contained && SCIPisLT(scip, ((SCIP_Real) consdata0->weights[v0]) / quotient, (SCIP_Real) consdata1->weights[v1]) )
            {
               iscons1incons0contained = FALSE;
               if( !iscons0incons1contained )
                  break;
            }
            /* if cons0 is possible contained in cons1 (consdata0->weight[v0] / quotient) must be less equals consdata1->weight[v1] */
            else if( iscons0incons1contained && SCIPisGT(scip, ((SCIP_Real) consdata0->weights[v0]) / quotient, (SCIP_Real) consdata1->weights[v1]) )
            {
               iscons0incons1contained = FALSE;
               if( !iscons1incons0contained )
                  break;
            }
            --v0;
            --v1;
         }
         else
         {
            /* both constraints have a variables which is not part of the other constraint, so stop */
            if( iscons0incons1contained && iscons1incons0contained )
            {
               iscons0incons1contained = FALSE;
               iscons1incons0contained = FALSE;
               break;
            }
            assert(iscons0incons1contained ? (v1 >= v0) : iscons1incons0contained);
            assert(iscons1incons0contained ? (v1 <= v0) : iscons0incons1contained);
            /* continue to the next variable */
            if( iscons0incons1contained )
               --v1;
            else
               --v0;
            ++v;
         }
      }
      if( iscons1incons0contained )
      {
         SCIPdebugMessage("knapsack constraint <%s> is redundant\n", SCIPconsGetName(cons1));
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ) );

         /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
         SCIP_CALL( updateFlags(scip, cons0, cons1) ); 

         SCIP_CALL( SCIPdelCons(scip, cons1) );
         ++(*ndelconss);
      }
      else if( iscons0incons1contained )
      {
         SCIPdebugMessage("knapsack constraint <%s> is redundant\n", SCIPconsGetName(cons0));
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ) );

         /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
         SCIP_CALL( updateFlags(scip, cons1, cons0) ); 

         SCIP_CALL( SCIPdelCons(scip, cons0) );
         ++(*ndelconss);
         break;
      }
   }
   
   return SCIP_OKAY;
}

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
                                              *   are separated as constraints. */
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
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyKnapsack)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrKnapsack(scip) );
 
   *valid = TRUE;

   return SCIP_OKAY;
}

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

   /* all variables which are of integral type can be binary; this can be checked via the method SCIPvarIsBinary(var) */
   nvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);
 
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

   /* all variables which are of integral type can be binary; this can be checked via the method SCIPvarIsBinary(var) */
   nvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);

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

   SCIP_Real loclowerbound;
   SCIP_Real glblowerbound;
   SCIP_Real cutoffbound;
   SCIP_Real maxbound;
  
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

   /* check, if we should additionally separate knapsack cuts */
   sepafreq = SCIPconshdlrGetSepaFreq(conshdlr);
   sepacardfreq = sepafreq * conshdlrdata->sepacardfreq;
   sepacardinality = (conshdlrdata->sepacardfreq >= 0)
      && ((sepacardfreq == 0 && depth == 0) || (sepacardfreq >= 1 && (depth % sepacardfreq == 0)));

   /* check dual bound to see if we want to produce knapsack cuts at this node */
   loclowerbound = SCIPgetLocalLowerbound(scip);
   glblowerbound = SCIPgetLowerbound(scip);
   cutoffbound = SCIPgetCutoffbound(scip);
   maxbound = glblowerbound + conshdlrdata->maxcardbounddist * (cutoffbound - glblowerbound);
   sepacardinality = sepacardinality && SCIPisLE(scip, loclowerbound, maxbound);
   sepacardinality = sepacardinality && (SCIPgetNLPBranchCands(scip) > 0);

   /* get the maximal number of cuts allowed in a separation round */
   maxsepacuts = (depth == 0 ? conshdlrdata->maxsepacutsroot : conshdlrdata->maxsepacuts);

   *result = SCIP_DIDNOTFIND;
   ncuts = 0;

   /* separate useful constraints */
   for( i = 0; i < nusefulconss && ncuts < maxsepacuts && !SCIPisStopped(scip); i++ )
   {
      SCIP_CALL( separateCons(scip, conss[i], NULL, sepacardinality, &ncuts) );
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

   /* check, if we should additionally separate knapsack cuts */
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
      SCIP_CALL( separateCons(scip, conss[i], sol, sepacardinality, &ncuts) );
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
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   SCIP_Bool redundant;
   int nfixedvars;
   int i;

   cutoff = FALSE;
   nfixedvars = 0;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* process useful constraints */
   for( i = 0; i < nusefulconss && !cutoff; i++ )
   {
      SCIP_CALL( propagateCons(scip, conss[i], &cutoff, &redundant, &nfixedvars, conshdlrdata->negatedclique) );
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
   SCIP_CONS* cons;
   SCIP_Bool cutoff;
   SCIP_Bool redundant;
   SCIP_Bool success;
   int oldnfixedvars;
   int oldnchgbds;
   int oldndelconss;
   int oldnaddconss;
   int oldnchgcoefs;
   int oldnchgsides;
   int firstchange;
   int c;

   /* remember old preprocessing counters */
   cutoff = FALSE;
   oldnfixedvars = *nfixedvars;
   oldnchgbds = *nchgbds;
   oldndelconss = *ndelconss;
   oldnaddconss = *naddconss;
   oldnchgcoefs = *nchgcoefs;
   oldnchgsides = *nchgsides;
   firstchange = INT_MAX;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   
   for( c = 0; c < nconss && !SCIPisStopped(scip); c++ )
   {
      int thisnfixedvars;
      int thisnchgbds;

      cons = conss[c];
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* force presolving the constraint in the initial round */
      if( nrounds == 0 )
         consdata->presolved = FALSE;
      else if( consdata->presolved )
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

      /* merge constraint, so propagation works better */
      SCIP_CALL( mergeMultiples(scip, cons, &cutoff) );
      if( cutoff )
         return SCIP_OKAY;
      
      /* add cliques in the knapsack to the clique table */
      SCIP_CALL( addCliques(scip, cons, &cutoff, nchgbds) );
      if( cutoff )
         break;
      
      /* propagate constraint */
      SCIP_CALL( propagateCons(scip, cons, &cutoff, &redundant, nfixedvars, TRUE) );
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
         SCIP_CALL( tightenWeights(scip, cons, nchgcoefs, nchgsides, naddconss, &cutoff) );
         if( cutoff )
            break;

         /* try to simplify inequalities */
         if( conshdlrdata->simplifyinequalities )
         {
	   SCIP_CALL( simplifyInequalities(scip, cons, ndelconss, nchgcoefs, nchgsides, &cutoff) );
	   if( cutoff )
	      break;
         }

         if( SCIPconsIsActive(cons) && conshdlrdata->dualpresolving )
         {
            /* in case the knapsack constraints is independent of everything else, solve the knapsack and apply the dual reduction */
            SCIP_CALL( dualPresolving(scip, cons, nchgbds, ndelconss, &redundant) );
            if( redundant )
               continue;
         }
      }
      /* remember the first changed constraint to begin the next aggregation round with */
      if( firstchange == INT_MAX && !consdata->presolved )
         firstchange = c;
   } 

   /* preprocess pairs of knapsack constraints */
   if( !cutoff && conshdlrdata->presolusehashing ) 
   {
      /* detect redundant constraints; fast version with hash table instead of pairwise comparison */
      SCIP_CALL( detectRedundantConstraints(scip, SCIPblkmem(scip), conss, nconss, &cutoff, ndelconss) );
   }

   if( (*ndelconss != oldndelconss) || (*nchgsides != oldnchgsides) || (*nchgcoefs != oldnchgcoefs) || (*naddconss != oldnaddconss) )
      success = TRUE;
   else
      success = FALSE;

   if( !cutoff && firstchange < nconss && conshdlrdata->presolpairwise )
   {
      SCIP_Longint npaircomparisons;
      
      npaircomparisons = 0;
      oldndelconss = *ndelconss;
      oldnchgsides = *nchgsides;
      oldnchgcoefs = *nchgcoefs;

      for( c = firstchange; c < nconss && !cutoff && !SCIPisStopped(scip); ++c )
      {
         cons = conss[c];
         if( !SCIPconsIsActive(cons) || SCIPconsIsModifiable(cons) )
            continue;
         
         npaircomparisons += (!SCIPconsGetData(cons)->presolved) ? c : (c - firstchange);
         
         SCIP_CALL( preprocessConstraintPairs(scip, conss, firstchange, c, ndelconss) );
         
         if( npaircomparisons > NMINCOMPARISONS )
         {

            if( (*ndelconss != oldndelconss) || (*nchgsides != oldnchgsides) || (*nchgcoefs != oldnchgcoefs) )
               success = TRUE;
            if( ((*ndelconss - oldndelconss) + (*nchgsides - oldnchgsides)/2 + (*nchgcoefs - oldnchgcoefs)/10) / (npaircomparisons + 0.0) < MINGAINPERNMINCOMPARISONS )
               break;
            oldndelconss = *ndelconss;
            oldnchgsides = *nchgsides;
            oldnchgcoefs = *nchgcoefs;
            npaircomparisons = 0;
         }
      }
   }

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( success || *nfixedvars > oldnfixedvars || *nchgbds > oldnchgbds )
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

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
 
   /* check if we fixed a binary variable to one (due to negated clique) */
   if( inferinfo >= 0 && SCIPvarGetLbLocal(infervar) > 0.5 )
   {
      for( i = 0; i < consdata->nvars; ++i )
      {
         if( SCIPvarGetIndex(consdata->vars[i]) == inferinfo ) 
         {
            assert( SCIPvarGetUbAtIndex(consdata->vars[i], bdchgidx, FALSE) < 0.5 );
            SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[i]) );
            break;
         }
      }
      assert(i < consdata->nvars);
   }
   else
   {
      /* according to negated cliques the minweightsum and all variables which are fixed to one which led to a fixing of
       * another negated clique variable to one, the inferinfo was choosen to be the negative of the position in the 
       * knapsack constraint, see one above call of SCIPinferBinvarCons
       */
      if( inferinfo < 0 )
         capsum = 0;
      else
      {
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
   }
   
   /* NOTE: It might be the case that capsum < consdata->capacity. This is due the fact that the fixing of the variable
    *       to zero can included negated clique information. A negated clique means, that at most one of the clique
    *       variables can be zero. These information can be used to compute a minimum activity of the constraint and
    *       used to fix variables to zero.
    *
    *       Even if capsum < consdata->capacity we still reported a complete reason since the minimum activity is based
    *       on global variable bounds. It might even be the case that we reported to many variables which are fixed to
    *       one.
    */
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
         -SCIPinfinity(scip), (SCIP_Real) SCIPgetCapacityKnapsack(sourcescip, sourcecons), varmap, consmap, 
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, global, valid) );
   assert(cons != NULL);

   SCIPfreeBufferArray(scip, &coefs);

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
#define consParseKnapsack NULL



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
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

   case SCIP_EVENTTYPE_LBTIGHTENED:
      eventdata->consdata->onesweightsum += eventdata->weight;
      eventdata->consdata->propagated = FALSE;
      eventdata->consdata->presolved = FALSE;
      break;
   case SCIP_EVENTTYPE_LBRELAXED:
      eventdata->consdata->onesweightsum -= eventdata->weight;
      
      conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
      assert(conshdlr != NULL);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
   
      if( conshdlrdata->negatedclique )
      {
         /* if a variable fixed to 1 is unfixed, it is possible, that it can be fixed to 1 again */
         eventdata->consdata->propagated = FALSE;
      }
  
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
         NULL,
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
         conshdlrCopyKnapsack,
         consFreeKnapsack, consInitKnapsack, consExitKnapsack, 
         consInitpreKnapsack, consExitpreKnapsack, consInitsolKnapsack, consExitsolKnapsack,
         consDeleteKnapsack, consTransKnapsack, consInitlpKnapsack,
         consSepalpKnapsack, consSepasolKnapsack, consEnfolpKnapsack, consEnfopsKnapsack, consCheckKnapsack, 
         consPropKnapsack, consPresolKnapsack, consRespropKnapsack, consLockKnapsack,
         consActiveKnapsack, consDeactiveKnapsack, 
         consEnableKnapsack, consDisableKnapsack,
         consPrintKnapsack, consCopyKnapsack, consParseKnapsack,
         conshdlrdata) );

   if( SCIPfindConshdlr(scip,"linear") != NULL )
   {
      /* include the linear constraint to knapsack constraint upgrade in the linear constraint handler */
      SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdKnapsack, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );
   }
  
   /* add knapsack constraint handler parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/knapsack/sepacardfreq",
         "multiplier on separation frequency, how often knapsack cuts are separated (-1: never, 0: only at root)",
         &conshdlrdata->sepacardfreq, TRUE, DEFAULT_SEPACARDFREQ, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/knapsack/maxcardbounddist",
         "maximal relative distance from current node's dual bound to primal bound compared to best node's dual bound for separating knapsack cuts",
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
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/knapsack/disaggregation",
         "should disaggregation of knapsack constraints be allowed in preprocessing?",
         &conshdlrdata->disaggregation, TRUE, DEFAULT_DISAGGREGATION, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/knapsack/simplifyinequalities",
         "should presolving try to simplify knapsacks",
         &conshdlrdata->simplifyinequalities, TRUE, DEFAULT_SIMPLIFYINEQUALITIES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/knapsack/negatedclique",
         "should negated clique information be used in solving process",
         &conshdlrdata->negatedclique, TRUE, DEFAULT_NEGATEDCLIQUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/knapsack/presolpairwise",
         "should pairwise constraint comparison be performed in presolving?",
         &conshdlrdata->presolpairwise, TRUE, DEFAULT_PRESOLPAIRWISE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/knapsack/presolusehashing",
         "should hash table be used for detecting redundant constraints in advance", 
         &conshdlrdata->presolusehashing, TRUE, DEFAULT_PRESOLUSEHASHING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/knapsack/dualpresolving",
         "should dual presolving steps be preformed?",
         &conshdlrdata->dualpresolving, TRUE, DEFAULT_DUALPRESOLVING, NULL, NULL) );

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
                                              *   are separated as constraints. */
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


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

/**@file   cons_pseudoboolean.c
 * @brief  constraint handler for pseudo Boolean constraints
 * @author Stefan Heinz
 * @author Michael Winkler
 *
 *
 * The constraint handler deals with pseudo Boolean constraints. These are constraints of the form 
 * \f[
 * \mbox{lhs} \leq \sum_{k=0}^m c_k \cdot x_k  +  \sum_{i=0}^n c_i \cdot \prod_{j \in I_i} x_j \leq \mbox{rhs}
 * \f]
 * where all x are binary and all c are integer
 *
 * @todo Add eventhandling.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/cons_and.h"
#include "scip/cons_indicator.h"
#include "scip/cons_linear.h"
#include "scip/cons_pseudoboolean.h"
#include "scip/pub_var.h"
#include "scip/debug.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "pseudoboolean"
#define CONSHDLR_DESC          "constraint handler template"
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define DEFAULT_DECOMPOSE          TRUE /**< decompose the pseudo boolean constraint into a "linear" constraint "and" constraints */
#define DEFAULT_SEPARATENONLINEAR  TRUE /**< if decomposed, should the nonlinear constraints be separated during LP processing */
#define DEFAULT_PROPAGATENONLINEAR TRUE /**< if decomposed, should the nonlinear constraints be propagated during node processing */
#define DEFAULT_REMOVABLENONLINEAR TRUE /**< if decomposed, should the nonlinear constraints be removable */
#define USEINDICATOR               TRUE

/*
 * Data structures
 */
#define HASHSIZE_PSEUDOBOOLEANNONLINEARTERMS 131101 /**< minimal size of hash table in and constraint tables */


/* struct used to find fast whether an and-constraint is new or already existing, and to memorize the non-linear parts
 * inside the opb constraints */
struct Term
{
   SCIP_VAR**            vars;               /**< non-linear variables */
   SCIP_VAR*             resultant;          /**< resultant variable in and constraint if transform, NULL otherwise */
   int                   nvars;              /**< number of variable inside this nonlinear term */
   int                   svars;              /**< size of variable array inside this nonlinear term */
   unsigned int          sorted:1;           /**< are all variables in this term sorted by index */
   unsigned int          merged:1;           /**< are all variables in this term merged */
   unsigned int          removedfixings:1;   /**< are all fixed variables removed */
};
typedef struct Term TERM;

/** constraint data for pseudoboolean constraints */
struct SCIP_ConsData
{
   SCIP_Real             lhs;                /**< left hand side of constraint */
   SCIP_Real             rhs;                /**< right hand side of constraint */
   
   SCIP_VAR**            linvars;            /**< linear variables */
   SCIP_Real*            lincoefs;           /**< coefficients of linear variables */
   int                   nlinvars;           /**< number of linear variables */
   int                   slinvars;           /**< size of the linvars- and lincoefs-arrays */

   TERM**                nonlinterms;        /**< array of nonlinear terms */
   SCIP_Real*            nonlincoefs;        /**< coefficients of nonlinear terms */
   int                   nnonlinterms;       /**< number of nonlinear terms */
   int                   snonlinterms;       /**< size of the nonlinterms- and nonlincoefs-arrays */
   
   SCIP_VAR*             intvar;             /**< a artificial variable which was added only for the objective function,
                                              *   if this variable is not NULL this constraint (without this integer
                                              *   variable) describes the objective function */

   SCIP_VAR*             indvar;             /**< indicator variable if it's a soft constraint, or NULL */
   SCIP_Real             weight;             /**< weight of the soft constraint, if it is one */

   unsigned int          sorted:1;           /**< are all variables of the linear part sorted by index */
   unsigned int          issoftcons:1;       /**< is this a soft constraint */
   unsigned int          propagated:1;       /**< is constraint already propagated? */
   unsigned int          boundstightened:1;  /**< is constraint already propagated with bound tightening? */
   unsigned int          presolved:1;        /**< is constraint already presolved? */
   unsigned int          removedfixings:1;   /**< are all fixed variables removed from the constraint? */
   unsigned int          changed:1;          /**< was constraint changed since last aggregation round in preprocessing? */
   unsigned int          normalized:1;       /**< is the constraint in normalized form? */
   unsigned int          upgradetried:1;     /**< was the constraint already tried to be upgraded? */
#if 0
   unsigned int          upgraded:1;         /**< is the constraint upgraded and will it be removed after preprocessing? */
#endif
   unsigned int          merged:1;           /**< are the constraint's equal linear variables already merged? */
   unsigned int          cliquesadded:1;     /**< were the cliques of the constraint already extracted? */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   TERM**                allnonlinterms;     /**< array of nonlinterms terms inside the whole problem */ 
   int                   nallnonlinterms;    /**< number of nonlinear terms inside the whole problem */
   int                   sallnonlinterms;    /**< size for all nonlinear terms inside the whole problem */
   SCIP_HASHTABLE*       hashtable;          /**< hash table of the nonlinear terms */
   int                   hashtablesize;      /**< size for hash table of the nonlinear terms */
   SCIP_Bool             decompose;          /**< decompose the pseudo boolean constraint into a "linear" constraint "and" constraints */
};

/*
 * Local methods
 */

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashGetKeyNonLinearTerms)
{  /*lint --e{715}*/
   /* the key is the element itself */ 
   return elem;
}

/** returns TRUE iff both keys are equal; two non-linear terms are equal if they have the same variables */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqNonLinearTerms)
{
#ifndef NDEBUG
   SCIP* scip;
#endif
   TERM* term1;
   TERM* term2;
   int v;

   term1 = (TERM*)key1;
   term2 = (TERM*)key2;
#ifndef NDEBUG
   scip = (SCIP*)userptr; 
#endif
   assert(scip != NULL);
   assert(term1 != NULL);
   assert(term2 != NULL);
   assert(term1->vars != NULL);
   assert(term1->nvars > 1);
#ifndef NDEBUG
   {
      /* check that these variables are sorted */
      for( v = term1->nvars - 1; v > 0; --v )
         assert(SCIPvarGetIndex(term1->vars[v]) >= SCIPvarGetIndex(term1->vars[v - 1]));
   }
#endif
   assert(term2->vars != NULL);
   assert(term2->nvars > 1);
#ifndef NDEBUG
   {
      /* check that these variables are sorted */
      for( v = term2->nvars - 1; v > 0; --v )
         assert(SCIPvarGetIndex(term2->vars[v]) >= SCIPvarGetIndex(term2->vars[v - 1]));
   }
#endif

   /* checks trivial case */
   if( term1->nvars != term2->nvars )
      return FALSE;

   /* check each variable in both terms for equality */
   for( v = term1->nvars - 1; v >= 0; --v )
   {
      assert(term1->vars[v] != NULL);
      assert(term2->vars[v] != NULL);

      /* tests if variables are equal */
      if( term1->vars[v] != term2->vars[v] )
      {
         assert(SCIPvarCompare(term1->vars[v], term2->vars[v]) == 1 || 
            SCIPvarCompare(term1->vars[v], term2->vars[v]) == -1);
         return FALSE;
      }
      assert(SCIPvarCompare(term1->vars[v], term2->vars[v]) == 0); 
   } 
   
   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValNonLinearTerms)
{  /*lint --e{715}*/
   TERM* term;
   int minidx;
   int mididx;
   int maxidx;
   
   term = (TERM*)key;
   
   assert(term != NULL);
   assert(term->vars != NULL);
   assert(term->nvars > 1);
#ifndef NDEBUG
   {
      /* check that these variables are sorted */
      int v;
      for( v = term->nvars - 1; v > 0; --v )
         assert(SCIPvarGetIndex(term->vars[v]) >= SCIPvarGetIndex(term->vars[v - 1]));
   }
#endif

   minidx = SCIPvarGetIndex(term->vars[0]);
   mididx = SCIPvarGetIndex(term->vars[term->nvars / 2]);
   maxidx = SCIPvarGetIndex(term->vars[term->nvars - 1]);
   assert(minidx >= 0 && minidx <= maxidx);

   return (term->nvars << 29) + (minidx << 22) + (mididx << 11) + maxidx; /*lint !e701*/
}

/** sort a non-linear term by index */
static
void sortTerm(
   TERM*                 term                /**< non-linear term to sort */
   )
{
   assert(term != NULL);
   assert(term->vars == NULL || term->nvars > 0);

   if( term->sorted )
      return;

   if( term->nvars < 2 )
   {
      term->sorted = TRUE;
      return;
   }

   /* sort vars to be able to check easily whether this non-linear term exists multiple times */
   SCIPsortPtr((void**)(term->vars), SCIPvarComp, term->nvars);
   term->sorted = TRUE;
}

/** creates constraint handler data for pseudo boolean constraint handler */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocMemory(scip, conshdlrdata) );

   return SCIP_OKAY;
}


/** frees constraint handler data for pseudo boolean constraint handler */
static
void conshdlrdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);

   SCIPfreeMemory(scip, conshdlrdata);
}

/** initializes constraint handler data for pseudo boolean constraint handler */
static
SCIP_RETCODE conshdlrdataInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);

   conshdlrdata->allnonlinterms = NULL;
   conshdlrdata->nallnonlinterms = 0;
   conshdlrdata->sallnonlinterms = 10;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(conshdlrdata->allnonlinterms), conshdlrdata->sallnonlinterms ) );

   /* create a hash table for nonlinear terms */
   conshdlrdata->hashtablesize = SCIPcalcHashtableSize(HASHSIZE_PSEUDOBOOLEANNONLINEARTERMS);
   SCIP_CALL( SCIPhashtableCreate(&(conshdlrdata->hashtable), SCIPblkmem(scip), conshdlrdata->hashtablesize,
         hashGetKeyNonLinearTerms, hashKeyEqNonLinearTerms, hashKeyValNonLinearTerms, (void*) scip) );

   return SCIP_OKAY;
}

/** deletes constraint handler data for pseudo boolean constraint handler if necessary */
static
void conshdlrdataExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->hashtable != NULL )
   {
      /* free hash table */
      SCIPhashtableFree(&(conshdlrdata->hashtable));
      conshdlrdata->hashtablesize = 0;
   }

   if( conshdlrdata->allnonlinterms != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->allnonlinterms), conshdlrdata->sallnonlinterms );

      conshdlrdata->allnonlinterms = NULL;
      conshdlrdata->nallnonlinterms = 0;
      conshdlrdata->sallnonlinterms = 0;
   }
}

/** creates a non-linear term data object */
static
SCIP_RETCODE createTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           terms,              /**< nonlinear terms of variables, or NULL */
   SCIP_Real*            coefs,              /**< array with coefficients for the nonlinear terms */
   int                   nterms,             /**< number of terms of variables of nonlinear term */
   int*                  ntermvars,          /**< number of variables in nonlinear terms, or NULL */
   TERM***               pbterms,            /**< pointer to store all created terms */
   SCIP_Real**           pbcoefs,            /**< pointer to store all important coefficients for the nonlinear terms */
   int*                  npbterms            /**< pointer to store the number of created terms */
   )
{
   TERM* term;
   int t;

   assert(scip != NULL);
   assert(nterms == 0 || (terms != NULL && ntermvars != NULL && coefs != NULL));
   assert(pbterms != NULL);
   assert(pbcoefs != NULL);
   assert(npbterms != NULL);

   *npbterms = 0;
   
   if( nterms == 0 )
   {
      *pbterms = NULL;
      *pbcoefs = NULL;
      return SCIP_OKAY;
   }

   /* allocate memory for the non-linear terms array */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, pbterms, nterms) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, pbcoefs, nterms) );

   /* create all non-linear terms */
   for( t = 0; t < nterms; ++t )
   {
      assert(terms[t] != NULL);
      assert(ntermvars[t] > 1);

      /* ignore nonlinear terms with a zero coefficient or zero variables */
      if( !SCIPisZero(scip, coefs[t]) && ntermvars[t] > 0 )
      {
         SCIP_Bool varsareactive;
         SCIP_Bool sorted;
         SCIP_Bool merged;
         int v;

         SCIP_CALL( SCIPallocBlockMemory(scip, &term) );
         
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(term->vars), terms[t], ntermvars[t]) );
         term->nvars = ntermvars[t];
         term->svars = ntermvars[t];
         term->resultant = NULL;

         /* compute some termdata booleans */
         varsareactive = SCIPvarIsActive(terms[t][ntermvars[t] - 1]);
         sorted = TRUE;
         merged = TRUE;
         
         for( v = ntermvars[t] - 2; v >= 0; --v )
         {
            varsareactive = varsareactive && SCIPvarIsActive(terms[t][v]);
            sorted = sorted && (SCIPvarCompare(terms[t][v + 1], terms[t][v]) >= 0);
            merged = merged && (SCIPvarCompare(terms[t][v + 1], terms[t][v]) > 0);
         }
         term->sorted = sorted;
         term->removedfixings = varsareactive;
         term->merged = merged && varsareactive;
        
         /* save term */
         (*pbterms)[*npbterms] = term;
         (*pbcoefs)[*npbterms] = coefs[t];
         ++(*npbterms);
      }
   }
   
   return SCIP_OKAY;
}

/** creates a pseudo boolean constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to linear constraint data */
   SCIP_VAR**            linvars,            /**< array with variables for the linear terms */
   SCIP_Real*            lincoefs,           /**< array with coefficients for the linear terms */
   int                   nlinvars,           /**< number of linear terms in the constraint */
   TERM**                nonlinterms,        /**< array with nonlinvars terms */
   SCIP_Real*            nonlincoefs,        /**< array with coefficients for the nonlinear terms */
   int                   nnonlinterms,       /**< number of nonlinear terms in the constraint */
   SCIP_VAR*             indvar,             /**< indicator variable if it's a soft constraint, or NULL */
   SCIP_Real             weight,             /**< weight of the soft constraint, if it is one */
   SCIP_Bool             issoftcons,         /**< is this a soft constraint */
   SCIP_VAR*             intvar,             /**< a artificial variable which was added only for the objective function,
                                              *   if this variable is not NULL this constraint (without this integer
                                              *   variable) describes the objective function */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs                 /**< right hand side of row */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(nlinvars == 0 || linvars != NULL);
   assert(nlinvars == 0 || lincoefs != NULL);
   assert(nnonlinterms == 0 || nonlinterms != NULL);
   assert(nnonlinterms == 0 || nonlincoefs != NULL);

   /* adjust right hand side */
   if( SCIPisInfinity(scip, rhs) )
      rhs = SCIPinfinity(scip);
   else if( SCIPisInfinity(scip, -rhs) )
      rhs = -SCIPinfinity(scip);
   
   /* adjust left hand side */
   if( SCIPisInfinity(scip, -lhs) )
      lhs = -SCIPinfinity(scip);
   else if( SCIPisInfinity(scip, lhs) )
      lhs = SCIPinfinity(scip);

   /* check left and right side */
   if( SCIPisGT(scip, lhs, rhs) )
   {
      SCIPerrorMessage("left hand side of pseudo boolean constraint greater than right hand side\n");
      SCIPerrorMessage(" -> lhs=%g, rhs=%g\n", lhs, rhs);
      return SCIP_INVALIDDATA;
   }

   /* allocate memory for the constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   (*consdata)->sorted = FALSE;

   /* initialize the weights for soft constraints */
   (*consdata)->issoftcons = issoftcons;
   if( issoftcons )
   {
      (*consdata)->weight = weight;
      if( SCIPisTransformed(scip) )
      {
         SCIP_CALL( SCIPgetTransformedVar(scip, indvar, &((*consdata)->indvar)) );
      }
      else
         (*consdata)->indvar = indvar;
   }
   else
      (*consdata)->indvar = NULL;
  
   /* copy artificial integer variable if it exist */
   if( intvar != NULL )
   { 
      if( SCIPisTransformed(scip) )
      {
         SCIP_CALL( SCIPgetTransformedVar(scip, intvar, &((*consdata)->intvar)) );
      }
      else
         (*consdata)->intvar = intvar;
   }
   else
      (*consdata)->intvar = NULL;

   /* copy linear terms */
   (*consdata)->slinvars = nlinvars;
   (*consdata)->nlinvars = nlinvars;
   if( nlinvars > 0 )
   {
      int k;
      int v;

      /* copy linear term variables and coefficients */
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->linvars, linvars, nlinvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->lincoefs, lincoefs, nlinvars) );

      /* remove linear term variables with a zero coefficient */
      k = 0;
      for( v = 0; v < nlinvars; ++v )
      {
         assert((*consdata)->linvars[v] != NULL);
         if( !SCIPisZero(scip, (*consdata)->lincoefs[v]) )
         {
            (*consdata)->linvars[k] = (*consdata)->linvars[v];
            (*consdata)->lincoefs[k] = (*consdata)->lincoefs[v];
            k++;
         }
      }
      (*consdata)->nlinvars = k;
   }
   else
   {
      (*consdata)->linvars = NULL;
      (*consdata)->lincoefs = NULL;
   }

   /* copy nonlinear terms */
   (*consdata)->snonlinterms = nnonlinterms;
   (*consdata)->nnonlinterms = nnonlinterms;
   if( nnonlinterms > 0 )
   {
      int t;
      int k;
      
      /* copy nonlinear coefficients */
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->nonlincoefs, nonlincoefs, nnonlinterms) );

      /* allocate memory for the nonlinear term array */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->nonlinterms, nnonlinterms) );

      assert((*consdata)->snonlinterms == nnonlinterms);

      /* copy each nonlinear terms */
      k = 0;
      for( t = 0; t < nnonlinterms; ++t )
      {
         TERM* term;
         
         /* ignore nonlinear terms with a zero coefficient or zero variables */
         if( !SCIPisZero(scip, (*consdata)->nonlincoefs[t]) && nonlinterms[t]->nvars > 0 )
         {
            SCIP_CALL( SCIPallocBlockMemory(scip, &term) );
         
            SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(term->vars), nonlinterms[t]->vars, nonlinterms[t]->nvars) );
            term->nvars = nonlinterms[t]->nvars; 
            term->svars = nonlinterms[t]->nvars; 
            term->resultant = nonlinterms[t]->resultant;
            term->sorted = FALSE;

            (*consdata)->nonlinterms[k] = term;
            (*consdata)->nonlincoefs[k] = (*consdata)->nonlincoefs[t];
            k++;
         }
      }
      (*consdata)->nnonlinterms = k;
   }
   else
   {
      (*consdata)->nonlinterms = NULL;
      (*consdata)->nonlincoefs = NULL;
   }
   
   /* copy left and right hand side */
   (*consdata)->lhs = lhs;
   (*consdata)->rhs = rhs;
   
   if( SCIPisTransformed(scip) )
   {
      int t;
      
      /* get transformed variables for the linear terms */
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nlinvars, (*consdata)->linvars, (*consdata)->linvars) );

      /* get transformed variables for the nonlinear terms */
      for( t = 0; t < nnonlinterms; ++t )
      {
         TERM* term;

         term = (*consdata)->nonlinterms[t];
         assert(term != NULL);
         
         SCIP_CALL( SCIPgetTransformedVars(scip, term->nvars, term->vars, term->vars) );
         if( term->resultant != NULL )
         {
            SCIP_CALL( SCIPgetTransformedVar(scip, term->resultant, &(term->resultant)) );
         }
      }
   }

   return SCIP_OKAY;
}

/** free a pseudo boolean constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to linear constraint data */
   )
{
   int t;

   assert(scip != NULL);
   assert(consdata != NULL);

   /* free each nonlinear term */
   for( t = 0; t < (*consdata)->nnonlinterms; ++t )
   {
      assert((*consdata)->nonlinterms[t] != NULL);
      assert(((*consdata)->nonlinterms[t]->vars) != NULL);

      SCIPfreeBlockMemoryArray(scip, &((*consdata)->nonlinterms[t]->vars), (*consdata)->nonlinterms[t]->nvars);
      SCIPfreeBlockMemory(scip, &((*consdata)->nonlinterms[t])); /*lint !e866 */
   }
   
   /* free linear term variable and coefficient arrays */
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->nonlinterms), (*consdata)->snonlinterms);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->nonlincoefs), (*consdata)->snonlinterms);

   /* free linear term variable and coefficient arrays */
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->lincoefs), (*consdata)->slinvars);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->linvars), (*consdata)->slinvars);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** installs rounding locks for the given variable associated to the given coefficient in the linear part */
static
SCIP_RETCODE lockRoundingVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisZero(scip, val));

   /* lock variable */
   if( SCIPisPositive(scip, val) )
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons,
            !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip, consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons,
            !SCIPisInfinity(scip, consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** removes rounding locks for the given variable associated to the given coefficient in the linear part */
static
SCIP_RETCODE unlockRoundingVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisZero(scip, val));

   /* unlock variable */
   if( SCIPisPositive(scip, val) )
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons,
            !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip, consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons,
            !SCIPisInfinity(scip, consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** installs rounding locks for the given term associated to the given coefficient in the nonlinear part */
static
SCIP_RETCODE lockRoundingTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   TERM*                 term,               /**< nonlinear term */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   SCIP_VAR** vars;
   SCIP_CONSDATA* consdata;
   SCIP_Bool islhsnotinf;
   SCIP_Bool isrhsnotinf;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(term != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisZero(scip, val));
   assert(term->nvars == 0 || term->vars != NULL);

   /* check which sites are infinity */
   islhsnotinf = !SCIPisInfinity(scip, -consdata->lhs);
   isrhsnotinf = !SCIPisInfinity(scip, consdata->rhs);

   vars = term->vars;

   /* locking variables */
   if( SCIPisPositive(scip, val) )
   {
      for( v = term->nvars - 1; v >= 0; --v )
      {
         SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, islhsnotinf, isrhsnotinf) );
      }
   }
   else
   {
      for( v = term->nvars - 1; v >= 0; --v )
      {
         SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, isrhsnotinf, islhsnotinf) );
      }
   }

   return SCIP_OKAY;
}

/** removes rounding locks for the given term associated to the given coefficient in the nonlinear part */
static
SCIP_RETCODE unlockRoundingTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   TERM*                 term,               /**< nonlinear term */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   SCIP_VAR** vars;
   SCIP_CONSDATA* consdata;
   SCIP_Bool islhsnotinf;
   SCIP_Bool isrhsnotinf;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(term != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisZero(scip, val));
   assert(term->nvars == 0 || term->vars != NULL);

   /* check which sites are infinity */
   islhsnotinf = !SCIPisInfinity(scip, -consdata->lhs);
   isrhsnotinf = !SCIPisInfinity(scip, consdata->rhs);

   vars = term->vars;

   /* unlock variables */
   if( SCIPisPositive(scip, val) )
   {
      for( v = term->nvars - 1; v >= 0; --v )
      {
         SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, islhsnotinf, isrhsnotinf) );
      }
   }
   else
   {
      for( v = term->nvars - 1; v >= 0; --v )
      {
         SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, isrhsnotinf, islhsnotinf) );
      }
   }

   return SCIP_OKAY;
}

/** adds a variable to the pseudo boolean constraint (if it is not zero) */
static
SCIP_RETCODE addCoefLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   if( SCIPisZero(scip, val) )
      return SCIP_OKAY;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* resize linear arrays if necessary */
   if( consdata->slinvars == consdata->nlinvars )
   {
      int oldsize;
      
      oldsize = consdata->slinvars;
      consdata->slinvars = SCIPcalcMemGrowSize(scip, consdata->slinvars + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->linvars), oldsize, consdata->slinvars) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->lincoefs), oldsize, consdata->slinvars) );
   }

   /* save all data */
   consdata->linvars[consdata->nlinvars] = var; 
   consdata->lincoefs[consdata->nlinvars] = val;
   ++(consdata->nlinvars);

   /* install rounding locks for new variable */
   SCIP_CALL( lockRoundingVar(scip, cons, var, val) );

   consdata->propagated = FALSE;
   consdata->boundstightened = FALSE;
   consdata->presolved = FALSE;
   consdata->removedfixings = consdata->removedfixings && SCIPvarIsActive(var);
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->cliquesadded = FALSE;
   if( consdata->nlinvars == 1 )
   {
      consdata->sorted = TRUE;
      consdata->merged = TRUE;
   }
   else
   {
      consdata->sorted = consdata->sorted
         && (SCIPvarCompare(consdata->linvars[consdata->nlinvars - 2], consdata->linvars[consdata->nlinvars - 1]) <= 0);
      consdata->merged = consdata->merged 
         && (SCIPvarCompare(consdata->linvars[consdata->nlinvars - 2], consdata->linvars[consdata->nlinvars - 1]) < 0);
   }

   return SCIP_OKAY;
}

/** create a non-linear term and adds it to the constraint */
static
SCIP_RETCODE addCoefTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_VAR**            vars,               /**< variables of the nonlinear term */
   int                   nvars,              /**< number of variables of the nonlinear term */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;
   TERM* term;
   SCIP_Bool varsareactive;
   SCIP_Bool sorted;
   int v;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(nvars == 0 || vars != NULL);

   if( nvars == 0 || SCIPisZero(scip, val) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* resize term arrays if necessary */
   if( consdata->snonlinterms == consdata->nnonlinterms )
   {
      int oldsize;
      
      oldsize = consdata->snonlinterms;
      consdata->snonlinterms = SCIPcalcMemGrowSize(scip, consdata->snonlinterms + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->nonlinterms), oldsize, consdata->snonlinterms) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->nonlincoefs), oldsize, consdata->snonlinterms) );
   }

   /* allocate memory for new term */
   SCIP_CALL( SCIPallocBlockMemory(scip, &term) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(term->vars), vars, nvars) );

   /* save all data */
   term->nvars = nvars;
   term->svars = nvars;
   term->resultant = NULL;

   /* compute some termdata booleans */
   varsareactive = SCIPvarIsActive(vars[nvars - 1]);
   sorted = TRUE;

   for( v = nvars - 2; v >= 0; --v )
   {
      varsareactive = varsareactive && SCIPvarIsActive(vars[v]);
      sorted = sorted && (SCIPvarCompare(vars[v + 1], vars[v]) >= 0);
   }
   term->sorted = sorted;
   term->removedfixings = varsareactive;
   term->merged = varsareactive;
        
   /* save term */
   consdata->nonlinterms[consdata->nnonlinterms] = term;
   consdata->nonlincoefs[consdata->nnonlinterms] = val;
   ++(consdata->nnonlinterms);

   /* install rounding locks for new variable */
   SCIP_CALL( lockRoundingTerm(scip, cons, term, val) );

   /* change flags */
   consdata->propagated = FALSE;
   consdata->boundstightened = FALSE;
   consdata->presolved = FALSE;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->cliquesadded = FALSE;

   return SCIP_OKAY;
}

/** adds a variable to a non-linear term */
static
SCIP_RETCODE addVarToTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   TERM*                 term,               /**< non-linear term */
   SCIP_VAR*             var                 /**< variable to add to non-linear part */
   )
{
   assert(scip != NULL);
   assert(term != NULL);
   assert(term->vars != NULL);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);

   if( term->nvars == term->svars )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, term->nvars + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(term->vars), term->svars, newsize) );
      term->svars = newsize;
   }
   assert(term->nvars < term->svars);

   term->vars[term->nvars] = var;
   ++(term->nvars);

   /* if new variables has a bigger index than the old last one, we don't need to change the sorted or merged status */
   if( term->nvars > 1 && SCIPvarGetIndex(var) <= SCIPvarGetIndex(term->vars[term->nvars - 2]) )
   {
      term->merged = FALSE;
      term->sorted = FALSE;
   }
           
   return SCIP_OKAY;            
}

/** deletes coefficient at given position from linear part of constraint */
static
SCIP_RETCODE delCoefLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< pseudoboolean constraint */
   int                   pos                 /**< position of coefficient of linear part to delete */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real val;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nlinvars);

   var = consdata->linvars[pos];
   val = consdata->lincoefs[pos];
   assert(var != NULL);

   /* remove rounding locks for deleted variable */
   SCIP_CALL( unlockRoundingVar(scip, cons, var, val) );

   /* move the last variable to the free slot */
   consdata->linvars[pos] = consdata->linvars[consdata->nlinvars-1];
   consdata->lincoefs[pos] = consdata->lincoefs[consdata->nlinvars-1];
   
   --(consdata->nlinvars);

   consdata->propagated = FALSE;
   consdata->boundstightened = FALSE;
   consdata->presolved = FALSE;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->cliquesadded = FALSE;

   consdata->sorted = consdata->sorted && (pos + 1 >= consdata->nlinvars || (SCIPvarCompare(consdata->linvars[pos], consdata->linvars[pos + 1]) <= 0));

   return SCIP_OKAY;
}

/** deletes a nonlinear term and at the given position from the nonlinear part of the constraint */
static
SCIP_RETCODE delCoefTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int                   pos                 /**< position of coefficient of nonlinear part to delete */
   )
{
   SCIP_CONSDATA* consdata;
   TERM* term;
   SCIP_Real val;
   
   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nnonlinterms);

   term = consdata->nonlinterms[pos];
   val = consdata->nonlincoefs[pos];
   assert(term != NULL);

   /* remove rounding locks for deleted variable */
   SCIP_CALL( unlockRoundingTerm(scip, cons, term, val) );
 
   /* remove locks on resultant */
   if( term->resultant != NULL )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, term->resultant, -1, -1) );
   }

   /* delete memory */
   SCIPfreeBlockMemoryArray(scip, &(term->vars), term->nvars);
   SCIPfreeBlockMemory(scip, &term);

   --(consdata->nnonlinterms);
   consdata->nonlinterms[pos] = consdata->nonlinterms[consdata->nnonlinterms];
   consdata->nonlincoefs[pos] = consdata->nonlincoefs[consdata->nnonlinterms];

   consdata->propagated = FALSE;
   consdata->boundstightened = FALSE;
   consdata->presolved = FALSE;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->cliquesadded = FALSE;

   return SCIP_OKAY;
}

/** deletes a variable in a non-linear term at the given position */
static
SCIP_RETCODE delVarInTermPos(
   SCIP*                 scip,               /**< SCIP data structure */
   TERM*                 term,               /**< non-linear term */
   int                   pos                 /**< position of variable in non-linear part to delete */
   )
{
   assert(scip != NULL);
   assert(term != NULL);
   assert(term->vars != NULL);
   assert(0 <= pos && pos < term->nvars);

   if( pos < term->nvars - 1 )
      term->vars[pos] = term->vars[term->nvars - 1];

   --(term->nvars);

   /* if deleted element was not one of the two last variables, we lost our sorting */
   if( pos < term->nvars )
      term->sorted = FALSE;
               
   return SCIP_OKAY;            
}

/** sets left hand side of linear constraint */
static
SCIP_RETCODE chgLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisInfinity(scip, lhs));

   /* adjust value to not be smaller than -inf */
   if( SCIPisInfinity(scip, -lhs) )
      lhs = -SCIPinfinity(scip);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisInfinity(scip, consdata->lhs));

   /* check whether the side is not changed */
   if( SCIPisEQ(scip, consdata->lhs, lhs) )
      return SCIP_OKAY;

   /* ensure that rhs >= lhs is satisfied without numerical tolerance */   
   if( SCIPisEQ(scip, lhs, consdata->rhs) )
      consdata->rhs = lhs;

   /* if necessary, update the rounding locks of variables */
   if( SCIPconsIsLocked(cons) )
   {
      if( SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, -lhs) )
      {
         SCIP_VAR** vars;
         SCIP_Real* lincoefs;
         int nlinvars;
         TERM** terms;
         SCIP_Real* termcoefs;
         int nterms;
         int v;
         int t;

         /* the right hand side switched from infinity to a non-infinite value -> install rounding locks */
         vars = consdata->linvars;
         lincoefs = consdata->lincoefs;
         nlinvars = consdata->nlinvars;
         terms = consdata->nonlinterms;
         termcoefs = consdata->nonlincoefs;
         nterms = consdata->nnonlinterms;

         /* linear part */
         for( v = nlinvars; v >= 0; --v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, lincoefs[v]));

            if( SCIPisPositive(scip, lincoefs[v]) )
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
            }
            else
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
            }
         }

         /* non-linear part */
         for( t = nterms - 1; t >= 0; --t )
         {
            assert(terms[t] != NULL);
            assert(terms[t]->nvars == 0 || terms[t]->vars != NULL);
            assert(!SCIPisZero(scip, termcoefs[t]));

            vars = terms[t]->vars;
            
            for( v = terms[t]->nvars - 1; v >= 0;--v )
            {
               assert(vars[v] != NULL);

               if( SCIPisPositive(scip, termcoefs[t]) )
               {
                  SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
               }
               else
               {
                  SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
               }
            }
         }
      }
      else if( !SCIPisInfinity(scip, -consdata->lhs) && SCIPisInfinity(scip, -lhs) )
      {
         SCIP_VAR** vars;
         SCIP_Real* lincoefs;
         int nlinvars;
         TERM** terms;
         SCIP_Real* termcoefs;
         int nterms;
         int v;
         int t;

         /* the right hand side switched from infinity to a non-infinite value -> install rounding locks */
         vars = consdata->linvars;
         lincoefs = consdata->lincoefs;
         nlinvars = consdata->nlinvars;
         terms = consdata->nonlinterms;
         termcoefs = consdata->nonlincoefs;
         nterms = consdata->nnonlinterms;

         /* linear part */
         for( v = nlinvars; v >= 0; --v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, lincoefs[v]));

            if( SCIPisPositive(scip, lincoefs[v]) )
            {
               SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
            }
            else
            {
               SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
            }
         }

         /* non-linear part */
         for( t = nterms - 1; t >= 0; --t )
         {
            assert(terms[t] != NULL);
            assert(terms[t]->nvars == 0 || terms[t]->vars != NULL);
            assert(!SCIPisZero(scip, termcoefs[t]));

            vars = terms[t]->vars;
            
            for( v = terms[t]->nvars - 1; v >= 0;--v )
            {
               assert(vars[v] != NULL);

               if( SCIPisPositive(scip, termcoefs[t]) )
               {
                  SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
               }
               else
               {
                  SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
               }
            }
         }
      }
   }

   /* check whether the left hand side is increased, if and only if that's the case we maybe can propagate, tighten and add more cliques */
   if( SCIPisLT(scip, consdata->lhs, lhs) )
   {
      consdata->propagated = FALSE;
      consdata->boundstightened = FALSE;
      consdata->cliquesadded = FALSE;
   }

   /* set new left hand side and update constraint data */
   consdata->lhs = lhs;
   consdata->presolved = FALSE;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;

   return SCIP_OKAY;
}

/** sets right hand side of linear constraint */
static
SCIP_RETCODE chgRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisInfinity(scip, -rhs));

   /* adjust value to not be larger than inf */
   if( SCIPisInfinity(scip, rhs) )
      rhs = SCIPinfinity(scip);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisInfinity(scip, -consdata->rhs));

   /* check whether the side is not changed */
   if( SCIPisEQ(scip, consdata->rhs, rhs) )
      return SCIP_OKAY;

   /* ensure that rhs >= lhs is satisfied without numerical tolerance */   
   if( SCIPisEQ(scip, rhs, consdata->lhs) )
      consdata->lhs = rhs;

   /* if necessary, update the rounding locks of variables */
   if( SCIPconsIsLocked(cons) )
   {
      assert(SCIPconsIsTransformed(cons));

      if( SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, rhs) )
      {
         SCIP_VAR** vars;
         SCIP_Real* lincoefs;
         int nlinvars;
         TERM** terms;
         SCIP_Real* termcoefs;
         int nterms;
         int v;
         int t;

         /* the right hand side switched from infinity to a non-infinite value -> install rounding locks */
         vars = consdata->linvars;
         lincoefs = consdata->lincoefs;
         nlinvars = consdata->nlinvars;
         terms = consdata->nonlinterms;
         termcoefs = consdata->nonlincoefs;
         nterms = consdata->nnonlinterms;

         /* linear part */
         for( v = nlinvars - 1; v >= 0;--v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, lincoefs[v]));

            if( SCIPisPositive(scip, lincoefs[v]) )
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
            }
            else
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
            }
         }

         /* non-linear part */
         for( t = nterms - 1; t >= 0; --t )
         {
            assert(terms[t] != NULL);
            assert(terms[t]->nvars == 0 || terms[t]->vars != NULL);
            assert(!SCIPisZero(scip, termcoefs[t]));

            vars = terms[t]->vars;
            
            for( v = terms[t]->nvars - 1; v >= 0;--v )
            {
               assert(vars[v] != NULL);

               if( SCIPisPositive(scip, termcoefs[t]) )
               {
                  SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
               }
               else
               {
                  SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
               }
            }
         }
      }
      else if( !SCIPisInfinity(scip, consdata->rhs) && SCIPisInfinity(scip, rhs) )
      {
         SCIP_VAR** vars;
         SCIP_Real* lincoefs;
         int nlinvars;
         TERM** terms;
         SCIP_Real* termcoefs;
         int nterms;
         int v;
         int t;

         /* the right hand side switched from infinity to a non-infinite value -> install rounding locks */
         vars = consdata->linvars;
         lincoefs = consdata->lincoefs;
         nlinvars = consdata->nlinvars;
         terms = consdata->nonlinterms;
         termcoefs = consdata->nonlincoefs;
         nterms = consdata->nnonlinterms;

         /* linear part */
         for( v = nlinvars - 1; v >= 0;--v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, lincoefs[v]));

            if( SCIPisPositive(scip, lincoefs[v]) )
            {
               SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
            }
            else
            {
               SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
            }
         }

         /* non-linear part */
         for( t = nterms - 1; t >= 0; --t )
         {
            assert(terms[t] != NULL);
            assert(terms[t]->nvars == 0 || terms[t]->vars != NULL);
            assert(!SCIPisZero(scip, termcoefs[t]));

            vars = terms[t]->vars;
            
            for( v = terms[t]->nvars - 1; v >= 0;--v )
            {
               assert(vars[v] != NULL);

               if( SCIPisPositive(scip, termcoefs[t]) )
               {
                  SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
               }
               else
               {
                  SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
               }
            }
         }
      }
   }

   /* check whether the right hand side is decreased, if and only if that's the case we maybe can propagate, tighten and add more cliques */
   if( SCIPisGT(scip, consdata->rhs, rhs) )
   {
      consdata->propagated = FALSE;
      consdata->boundstightened = FALSE;
      consdata->cliquesadded = FALSE;
   }

   /* set new right hand side and update constraint data */
   consdata->rhs = rhs;
   consdata->presolved = FALSE;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;

   return SCIP_OKAY;
}

/* create and-constraints and get all and-resultants */
static
SCIP_RETCODE createAnds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< pseudoboolean constraint */
   SCIP_CONSHDLR*        conshdlr,           /**< pseudoboolean constraint handler */
   SCIP_VAR***           andvars,            /**< pointer to store created and-resultant vars, should be as big as the
                                              *   size of consdata->nnonlinterms of cons */
   int*                  nandvars            /**< pointer to store number of created and-resultant vars */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int nmaxandres;
   int t;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(conshdlr != NULL);
   assert(andvars != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nmaxandres = consdata->nnonlinterms;
   *nandvars = 0;

   if( nmaxandres == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, andvars, nmaxandres) );
   BMSclearMemoryArray(*andvars, nmaxandres);

   for( t = 0; t < nmaxandres; ++t )
   {
      TERM* term;
      TERM* tmpterm;
      SCIP_CONS* andcons;
      SCIP_VAR* var;
      char name[SCIP_MAXSTRLEN];
      SCIP_Bool separate;
      SCIP_Bool propagate;
      SCIP_Bool removable;
      SCIP_Bool created;

      /* here we only allow terms with non zero coefficients */
      assert(!SCIPisZero(scip, consdata->nonlincoefs[t]));

      term = consdata->nonlinterms[t];
      assert(term != NULL);
      /* here we only allow non empty terms */
      assert(term->nvars > 0);

      /* sort variables in non-linear term */
      sortTerm(term);
      assert(term->sorted);

      assert(conshdlrdata->hashtable != NULL);

      /* get constraint from current hash table with same variables as cons0 */
      tmpterm = (TERM*)(SCIPhashtableRetrieve(conshdlrdata->hashtable, (void*)term));

      /* if there is already the same and constraint created use this resultant */
      if( tmpterm != NULL )
      {
         var = tmpterm->resultant;
         assert(var != NULL);
      }
      else 
         var = NULL;

      created = FALSE;

      /* create new and-resultant variable */
      if( var == NULL )
      {
         /* resize data for all and constraints if necessary */
         if(conshdlrdata->nallnonlinterms == conshdlrdata->sallnonlinterms)
         {
            int oldsize;

            oldsize = conshdlrdata->sallnonlinterms;
            conshdlrdata->sallnonlinterms = SCIPcalcMemGrowSize(scip, conshdlrdata->sallnonlinterms + 1);
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(conshdlrdata->allnonlinterms), oldsize, conshdlrdata->sallnonlinterms) );
         }

         conshdlrdata->allnonlinterms[conshdlrdata->nallnonlinterms] = term;
         ++(conshdlrdata->nallnonlinterms);

         /* no such non-linear term in current hash table: insert the term into hash table */
         SCIP_CALL( SCIPhashtableInsert(conshdlrdata->hashtable, (void*) term) );

         /* create auxiliary variable */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_andresultant_%d", SCIPconsGetName(cons), conshdlrdata->nallnonlinterms - 1);
         SCIP_CALL( SCIPcreateVar(scip, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY, 
               TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );

         /* change branching priority of artificial variable to -1 */
         SCIP_CALL( SCIPchgVarBranchPriority(scip, var, -1) );

         /* add auxiliary variable to the problem */
         SCIP_CALL( SCIPaddVar(scip, var) );
         created = TRUE;

#ifdef SCIP_DEBUG_SOLUTION
         {
            SCIP_Real val;
            int v;

            for( v = term->nvars - 1; v >= 0; --v )
            {
               SCIP_CALL( SCIPdebugGetSolVal(scip, term->vars[v], &val) );
               assert(SCIPisFeasZero(scip, val) || SCIPisFeasEQ(scip, val, 1.0));

               if( val < 0.5 )
                  break;
            }
            SCIP_CALL( SCIPdebugAddSolVal(scip, var, (val < 0.5) ? 0.0 : 1.0) );
         }
#endif

         SCIP_CALL( SCIPgetBoolParam(scip, "constraints/"CONSHDLR_NAME"/nlcseparate", &separate) );
         SCIP_CALL( SCIPgetBoolParam(scip, "constraints/"CONSHDLR_NAME"/nlcpropagate", &propagate) );
         SCIP_CALL( SCIPgetBoolParam(scip, "constraints/"CONSHDLR_NAME"/nlcremovable", &removable) );

         /* create, add, and release "and" constraint for the multiplication of the binary variables */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_andcons_%d", SCIPconsGetName(cons), conshdlrdata->nallnonlinterms - 1);
         SCIP_CALL( SCIPcreateConsAnd(scip, &andcons, name, var, term->nvars, term->vars,
               SCIPconsIsInitial(cons), separate, SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons), propagate, 
               SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), removable, SCIPconsIsStickingAtNode(cons)) );
         SCIP_CALL( SCIPaddCons(scip, andcons) );
         SCIPdebugPrintCons(scip, andcons, NULL);
         SCIP_CALL( SCIPreleaseCons(scip, &andcons) );
      }
      assert(var != NULL);
      term->resultant = var;

      (*andvars)[*nandvars] = var;
      ++(*nandvars);

      if( created )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &var) );
      }
   }

   return SCIP_OKAY;
}


/** checks pseudo boolean constraint for feasibility of given solution or current solution */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< pseudo boolean constraint */
   SCIP_SOL*             sol,                /**< solution to be checked, or NULL for current solution */
   SCIP_Bool             printreason,        /**< print reason of violation (if any) */
   SCIP_Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real activity;
   int t;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(violated != NULL);

   SCIPdebugMessage("checking pseudo boolean constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebugPrintCons(scip, cons, NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *violated = FALSE;
   activity = 0.0;

   /* compute linear activity */
   for( t = 0; t < consdata->nlinvars; ++t )
      activity += consdata->lincoefs[t] * SCIPgetSolVal(scip, sol, consdata->linvars[t]);

   /* add nonlinear activity */
   for( t = 0; t < consdata->nnonlinterms; ++t )
   {
      TERM* term;
      SCIP_Real factor;
      int v;

      term = consdata->nonlinterms[t];
      factor = consdata->nonlincoefs[t];

      /* maybe take value of resultant */
      if( term->resultant != NULL )
      {
#ifndef NDEBUG
         SCIP_Real andvalue;

         andvalue = 1;
         /* check if the term is active */
         for( v = 0; v < term->nvars; ++v )
         {
            andvalue *= SCIPgetSolVal(scip, sol, term->vars[v]);
            if( SCIPisFeasZero(scip, andvalue) )
               break;
         }
         assert(SCIPisFeasEQ(scip, andvalue, SCIPgetSolVal(scip, sol, term->resultant)));
#endif
         factor *= SCIPgetSolVal(scip, sol, term->resultant);
      }
      else
      {
         /* check if the term is active */
         for( v = 0; v < term->nvars && !SCIPisZero(scip, factor); ++v )
         {
            factor *= SCIPgetSolVal(scip, sol, term->vars[v]);
         }
      }

      activity += factor;
   }

   SCIPdebugMessage("  consdata activity=%.15g (lhs=%.15g, rhs=%.15g, sol=%p)\n",
      activity, consdata->lhs, consdata->rhs, (void*)sol);

   if( SCIPisFeasLT(scip, activity, consdata->lhs) || SCIPisFeasGT(scip, activity, consdata->rhs) )
   {
      *violated = TRUE;

      /* only reset constraint age if we are in enforcement */
      if( sol == NULL )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }

      if ( printreason )
      {
         SCIPinfoMessage(scip, NULL, "violation <%s>: %g <= ", SCIPconsGetName(cons), consdata->lhs);
         /* print linear part */
         for( t = 0; t < consdata->nlinvars; ++t )
         {
            if ( t > 0 )
            {
               if ( consdata->lincoefs[t] > 0 )
                  SCIPinfoMessage(scip, NULL, " +");
            }
            SCIPinfoMessage(scip, NULL, " %g <%s> [%g]", consdata->lincoefs[t], SCIPvarGetName(consdata->linvars[t]),
               SCIPgetSolVal(scip, sol, consdata->linvars[t]));
         }
         for( t = 0; t < consdata->nnonlinterms; ++t )
         {
            TERM* term;
            SCIP_Real factor;
            int v;

            term = consdata->nonlinterms[t];
            factor = consdata->nonlincoefs[t];

            if ( t > 0 )
               SCIPinfoMessage(scip, NULL, "+ ");
            SCIPinfoMessage(scip, NULL, "%g ", consdata->nonlincoefs[t]);

            /* maybe take value of resultant */
            if( term->resultant != NULL )
               SCIPinfoMessage(scip, NULL, " * <%s> [%g]", SCIPvarGetName(term->resultant), SCIPgetSolVal(scip, sol, term->resultant));
            else
            {
               for( v = 0; v < term->nvars && !SCIPisZero(scip, factor); ++v )
                  SCIPinfoMessage(scip, NULL, " * <%s> [%g]", SCIPvarGetName(term->vars[v]), SCIPgetSolVal(scip, sol, term->vars[v]));
            }
         }
         SCIPinfoMessage(scip, NULL, " <= %g\n", consdata->rhs);
      }
   }
   else
   {
      *violated = FALSE;
      SCIP_CALL( SCIPincConsAge(scip, cons) );
   }

   return SCIP_OKAY;
}

/** creates by copying and captures a linear constraint */
static
SCIP_RETCODE copyConsPseudoboolean(
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_CONS**           targetcons,         /**< pointer to store the created target constraint */
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP_CONS*            sourcecons,         /**< source constraint which will be copied */
   const char*           name,               /**< name of constraint */
   SCIP_HASHMAP*         varmap,             /**< a SCIP_HASHMAP mapping variables of the source SCIP to corresponding
                                              *   variables of the target SCIP */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   SCIP_Bool             stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool*            valid               /**< pointer to store if the copying was valid */
   )
{
   SCIP_VAR** linvars;
   SCIP_Real* lincoefs;
   int nlinvars;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_VAR* intvar;
   SCIP_VAR* indvar;

   SCIP_VAR*** termvars;
   SCIP_Real* termcoefs;
   int* ntermvars;
   int nterms;

   SCIP_CONSDATA* sourceconsdata;

   int v;

   assert(targetscip != NULL);
   assert(targetcons != NULL);
   assert(sourcescip != NULL);
   assert(sourcecons != NULL);

   sourceconsdata = SCIPconsGetData(sourcecons);
   assert(sourceconsdata != NULL);

   (*valid) = TRUE;

   linvars = NULL;
   lincoefs = NULL;
   nlinvars = 0;
   termvars = NULL;
   ntermvars = NULL;
   termcoefs = NULL;
   nterms = 0;

   /* first start with the linear part */
   if( sourceconsdata->nlinvars > 0 )
   {
      assert(sourceconsdata->linvars != NULL);
      assert(sourceconsdata->lincoefs != NULL);
      
      nlinvars = sourceconsdata->nlinvars;

      /* allocate temporary variable array */
      SCIP_CALL( SCIPduplicateBufferArray(targetscip, &linvars, sourceconsdata->linvars, nlinvars) );
      /* duplicate coefficient array */
      SCIP_CALL( SCIPduplicateBufferArray(targetscip, &lincoefs, sourceconsdata->lincoefs, nlinvars) );
   
      /* map variables of the source constraint to variables of the target SCIP */
      for( v = 0; v < nlinvars; ++v )
      {
         assert(linvars[v] != NULL);
         SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, linvars[v], &linvars[v], varmap, consmap, global, valid) );
         assert(!(*valid) || linvars[v] != NULL);

         /* we do not copy, if a variable is missing */
         if( !(*valid) )
            goto TERMINATE;
      }
   }

   /* second the nonlinear part */
   if( sourceconsdata->nnonlinterms > 0 )
   {
      TERM* term;
      int t;

      assert(sourceconsdata->nonlinterms != NULL);
      assert(sourceconsdata->nonlincoefs != NULL);
      assert(sourceconsdata->nonlincoefs != NULL);
      
      nterms = sourceconsdata->nnonlinterms;

      /* allocate nonlinear array of variable array */
      SCIP_CALL( SCIPallocBufferArray(targetscip, &termvars, nterms) );
      BMSclearMemoryArray(termvars, nterms);

      /* allocate array for length for all nonlinear variable array */
      SCIP_CALL( SCIPallocBufferArray(targetscip, &ntermvars, nterms) );
      /* duplicate coefficient array */
      SCIP_CALL( SCIPduplicateBufferArray(targetscip, &termcoefs, sourceconsdata->nonlincoefs, nterms) );
   
      /* copy all term variables */
      for( t = 0; t < nterms; ++t )
      { 
         term = sourceconsdata->nonlinterms[t];
         assert(term != NULL);

         assert(term->nvars > 1);
         ntermvars[t] = term->nvars;

         /* allocate nonlinear array of variable array */
         SCIP_CALL( SCIPallocBufferArray(targetscip, &(termvars[t]), ntermvars[t]) ); /*lint !e866 */

         /* map variables of the source constraint to variables of the target SCIP */
         for( v = 0; v < ntermvars[t]; ++v )
         {
            SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, term->vars[v], &(termvars[t][v]), varmap, consmap, global, valid) );
            assert(!(*valid) || termvars[t][v] != NULL);
            
            /* we do not copy, if a variable is missing */
            if( !(*valid) )
               goto TERMINATE;
         }
      }
   }

   /* third the indicator and artificial integer variable part */
   assert(sourceconsdata->issoftcons == (sourceconsdata->indvar != NULL));
   indvar = sourceconsdata->indvar;
   intvar = sourceconsdata->intvar;

   if( indvar != NULL )
   {
      assert(*valid);
      SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, indvar, &indvar, varmap, consmap, global, valid) );
      assert(!(*valid) || indvar != NULL);
   }
   if( intvar != NULL && *valid )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, intvar, &intvar, varmap, consmap, global, valid) );
      assert(!(*valid) || intvar != NULL);
   }

   /* the sides */
   lhs = sourceconsdata->lhs;   
   rhs = sourceconsdata->rhs;   

   if( *valid )
   {
      SCIP_CALL( SCIPcreateConsPseudoboolean(targetscip, targetcons, name, linvars, nlinvars, lincoefs, 
            termvars, nterms, ntermvars, termcoefs,
            indvar, sourceconsdata->weight, sourceconsdata->issoftcons, intvar, lhs, rhs, 
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }  

 TERMINATE:
   /* free buffer array */
   if( sourceconsdata->nnonlinterms > 0 )
   {
      int t;

      if( termvars != NULL )
      {
         for( t = sourceconsdata->nnonlinterms - 1; t >= 0; --t )
            SCIPfreeBufferArrayNull(targetscip, &(termvars[t]));
      }

      SCIPfreeBufferArrayNull(targetscip, &termcoefs);
      SCIPfreeBufferArrayNull(targetscip, &ntermvars);
      SCIPfreeBufferArrayNull(targetscip, &termvars);
   }
   if( sourceconsdata->nlinvars > 0 )
   {
      SCIPfreeBufferArrayNull(targetscip, &lincoefs);
      SCIPfreeBufferArrayNull(targetscip, &linvars);
   }
   
   return SCIP_OKAY;
}

/** replaces all fixed and aggregated variables in linear part of constraint */
static
SCIP_RETCODE applyFixingsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< pseudoboolean constraint */
   )
{
   SCIP_CONSDATA* consdata;
   int v;
      
   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* start removing fixings in linear part */
   if( !consdata->removedfixings && consdata->nlinvars > 0 )
   {
      SCIP_VAR* var;
      SCIP_Real val;
      SCIP_Real fixedval;
      SCIP_Real aggrconst;
      SCIP_Real lhssubtrahend;
      SCIP_Real rhssubtrahend;
#if 0 /* only for multi-aggregated variables, should not be needed here */
      SCIP_VAR** aggrvars;
      SCIP_Real* aggrscalars;
      int naggrvars;
      int i;
#endif
      lhssubtrahend = 0.0;
      rhssubtrahend = 0.0;

      SCIPdebugMessage("applying fixings linear:\n");
      SCIPdebugPrintCons(scip, cons, NULL);

      assert(consdata->linvars != NULL);
      assert(consdata->lincoefs != NULL);

      v = 0;
      /* iterate over all variables and replace all fixed variables */
      while( v < consdata->nlinvars )
      {
         var = consdata->linvars[v];
         val = consdata->lincoefs[v];

         assert(var != NULL);
         assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
         assert(SCIPvarIsTransformed(var));

         switch( SCIPvarGetStatus(var) )
         {
         case SCIP_VARSTATUS_ORIGINAL:
            SCIPerrorMessage("original variable in transformed linear constraint\n");
            return SCIP_INVALIDDATA;

         case SCIP_VARSTATUS_LOOSE:
         case SCIP_VARSTATUS_COLUMN:
            ++v;
            break;

         case SCIP_VARSTATUS_FIXED:
            assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)));
            fixedval = SCIPvarGetLbGlobal(var);
            
            if( !SCIPisInfinity(scip, -consdata->lhs) )
               lhssubtrahend += val * fixedval;
            
            if( !SCIPisInfinity(scip, consdata->rhs) )
                  rhssubtrahend += val * fixedval;

            SCIP_CALL( delCoefLinear(scip, cons, v) );
            break;

         case SCIP_VARSTATUS_AGGREGATED:
            SCIP_CALL( addCoefLinear(scip, cons, SCIPvarGetAggrVar(var), val * SCIPvarGetAggrScalar(var)) );
            aggrconst = SCIPvarGetAggrConstant(var);

            if( !SCIPisInfinity(scip, -consdata->lhs) )
               lhssubtrahend += val * aggrconst;
            if( !SCIPisInfinity(scip, consdata->rhs) )
               rhssubtrahend += val * aggrconst;

            SCIP_CALL( delCoefLinear(scip, cons, v) );
            break;

#if 0 /* binary variables should not be multi-aggregated */
         case SCIP_VARSTATUS_MULTAGGR:
            SCIP_CALL( SCIPflattenVarAggregationGraph(scip,var) );
            naggrvars = SCIPvarGetMultaggrNVars(var);
            aggrvars = SCIPvarGetMultaggrVars(var);
            aggrscalars = SCIPvarGetMultaggrScalars(var);
            for( i = 0; i < naggrvars; ++i )
            {
               SCIP_CALL( addCoef(scip, cons, aggrvars[i], val * aggrscalars[i]) );
            }
            aggrconst = SCIPvarGetMultaggrConstant(var);

            if( !SCIPisInfinity(scip, -consdata->lhs) )
               lhssubtrahend += val * aggrconst;
            if( !SCIPisInfinity(scip, consdata->rhs) )
               rhssubtrahend += val * aggrconst;

            SCIP_CALL( delCoefLinear(scip, cons, v) );
            break;
#else
         case SCIP_VARSTATUS_MULTAGGR:
            SCIPerrorMessage("binary variable should not be multi-aggregated\n");
            SCIPABORT();
            break;
#endif
         case SCIP_VARSTATUS_NEGATED:
            SCIP_CALL( addCoefLinear(scip, cons, SCIPvarGetNegationVar(var), -val) );
            aggrconst = SCIPvarGetNegationConstant(var);

            if( !SCIPisInfinity(scip, -consdata->lhs) )
               lhssubtrahend += val * aggrconst;
            if( !SCIPisInfinity(scip, consdata->rhs) )
               rhssubtrahend += val * aggrconst;

            SCIP_CALL( delCoefLinear(scip, cons, v) );
            break;

         default:
            SCIPerrorMessage("unknown variable status\n");
            SCIPABORT();
         }
      }
      
      /* update left hand side */
      if( !SCIPisInfinity(scip, -consdata->lhs) )
      {
         if( SCIPisFeasEQ(scip, lhssubtrahend, consdata->lhs ) )
         {
            SCIP_CALL( chgLhs(scip, cons, 0.0) );
         }
         else
         {
            assert(!SCIPisInfinity(scip, consdata->lhs - lhssubtrahend));
            SCIP_CALL( chgLhs(scip, cons, consdata->lhs - lhssubtrahend) );
         }
      }
      /* update right hand side */
      if( !SCIPisInfinity(scip, consdata->rhs) )
      {
         if( SCIPisFeasEQ(scip, rhssubtrahend, consdata->rhs ) )
         {
            SCIP_CALL( chgRhs(scip, cons, 0.0) );
         }
         else
         {
            assert(!SCIPisInfinity(scip, consdata->rhs - rhssubtrahend));
            SCIP_CALL( chgRhs(scip, cons, consdata->rhs - rhssubtrahend) );
         }
      }

      consdata->removedfixings = TRUE;

      SCIPdebugMessage("after fixings:\n");
      SCIPdebugPrintCons(scip, cons, NULL);
   }
   assert(consdata->removedfixings || consdata->nlinvars == 0);

#ifndef NDEBUG
   /* check, if all linear fixings are applied */
   for( v = 0; v < consdata->nlinvars; ++v )
      assert(SCIPvarIsActive(consdata->linvars[v]));
#endif

   return SCIP_OKAY;
}

/** replaces all fixed and aggregated variables in non-linear part */
static
SCIP_RETCODE applyFixingsNonLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< pseudoboolean constraint */
   int*                  nfixedvars,         /**< pointer to number of fixed variables */
   int*                  naggrvars,          /**< pointer to number of aggregated variables */
   SCIP_Bool*            infeasible          /**< pointer to store if infeasibility is detected */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real lhssubtrahend;
   SCIP_Real rhssubtrahend;
   int v;
   int t;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(infeasible != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   lhssubtrahend = 0.0;
   rhssubtrahend = 0.0;

   /* iterate over all non-linear parts */
   t = 0;
   while( t < consdata->nnonlinterms )
   {
      SCIP_VAR* resultant;
      SCIP_VAR** vars;
      int nvars;
      SCIP_VAR* var;
      SCIP_Real val;
      SCIP_Bool fixed;

      resultant = consdata->nonlinterms[t]->resultant;
      vars = consdata->nonlinterms[t]->vars;
      nvars = consdata->nonlinterms[t]->nvars;

      assert(vars != NULL);
      assert(nvars > 2);

      val = consdata->nonlincoefs[t];

      /* if the resultant exist and is fixed */
      if( resultant != NULL && SCIPvarGetStatus(resultant) == SCIP_VARSTATUS_FIXED )
      {
         SCIP_Real fixedval;
         
         assert(SCIPvarGetType(resultant) == SCIP_VARTYPE_BINARY);
         assert(SCIPvarIsTransformed(resultant));
         assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(resultant), SCIPvarGetUbGlobal(resultant)));
         
         fixedval = SCIPvarGetLbGlobal(resultant);
            
         /* fix non-linear part if not already fixed, if resultant is one than all term variables should be one too */
         if( fixedval > 0.5 )
         {
            if( !SCIPisInfinity(scip, -consdata->lhs) )
               lhssubtrahend += val;
            
            if( !SCIPisInfinity(scip, consdata->rhs) )
               rhssubtrahend += val;

            for( v = nvars - 1; v >= 0; --v )
            {
               var = vars[v];
               assert(var != NULL);
               assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);

               /* var could not be fixed to zero, otherwise the and-resultant couldn't be fixed to one */
               if( SCIPvarGetUbLocal(var) < 0.5 )
               {
                  SCIPdebugMessage("the and-resultant <%s> which is fixed to 1 is in conflict with the and-variable <%s> which is fixed to zero\n", SCIPvarGetName(resultant), SCIPvarGetName(var));
                  *infeasible = TRUE;
                  return SCIP_OKAY;
               }

               if( SCIPvarGetLbLocal(var) < 0.5 )
               {
                  SCIPdebugMessage("fixing variable <%s> to 1, due to fixed and-resultant\n", SCIPvarGetName(var));
                  SCIP_CALL( SCIPfixVar(scip, var, 1.0, infeasible, &fixed) );
                  assert(fixed);
                  assert(!(*infeasible));
                  ++(*nfixedvars);
               }
            }
         }
            
         SCIPdebugMessage("removing term at pos %d(in the non-linear part), due to fixed and-resultant\n", t);
         SCIP_CALL( delCoefTerm(scip, cons, t) );
      }
      else if( !consdata->nonlinterms[t]->removedfixings )
      {
         SCIP_Bool delterm;
         SCIP_Bool allone;

         delterm = FALSE;
         allone = TRUE;

         v = 0;
         /* delete all variables which are already fixed in that term and if an and-resultant exist fix it if possible */
         /* note that number of term variables can decrease in this loop due to deleting */
         while( v < consdata->nonlinterms[t]->nvars )
         {
            var = vars[v];
            assert(var != NULL);
            assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
            assert(SCIPvarIsTransformed(var));

            /* if one variable is fixed to zero the term is redundant and the and-resultant have to be fixed to zero too */
            if( SCIPvarGetUbLocal(var) < 0.5 )
            {
               if( resultant != NULL )
               {
                  SCIPdebugMessage("fixing and-resultant <%s> to 0\n", SCIPvarGetName(resultant));
                  SCIP_CALL( SCIPfixVar(scip, resultant, 0.0, infeasible, &fixed) );
                  assert(fixed);
                  assert(!(*infeasible));
                  ++(*nfixedvars);
               }
               allone = FALSE;
               delterm = TRUE;
               break;
            }
               
            if( SCIPvarGetLbGlobal(var) > 0.5 )
            {
               /* delete variable in non-linear part */
               SCIP_CALL( delVarInTermPos(scip, consdata->nonlinterms[t], v) );
            }
            else
            {
               SCIP_VAR* repvar;
               SCIP_Bool negated;

               /* get binary representative of variable */
               SCIP_CALL( SCIPgetBinvarRepresentative(scip, var, &repvar, &negated) );

               /* check, if the variable should be replaced with the representative */
               if( repvar != var )
               {
                  /* delete old aggregated variable in non-linear part */
                  SCIP_CALL( delVarInTermPos(scip, consdata->nonlinterms[t], v) );

                  /* add representative variable instead to non-linear part */
                  SCIP_CALL( addVarToTerm(scip, consdata->nonlinterms[t], repvar) );
               }
               else
               {
                  allone = FALSE;
                  ++v;
               }
            }
         }
         
         /* all variables inside the non-linear part are fixed to one then we can delete this part (and possible fix the and-resultant) */
         if( allone )
         {
            if( resultant != NULL )
            {
               SCIPdebugMessage("fixing and-resultant <%s> to 1, all and-variables are fixed to 1\n", SCIPvarGetName(resultant));
               SCIP_CALL( SCIPfixVar(scip, resultant, 1.0, infeasible, &fixed) );
               assert(fixed);
               assert(!(*infeasible));
               ++(*nfixedvars);
            }
            delterm = TRUE;

            if( !SCIPisInfinity(scip, -consdata->lhs) )
               lhssubtrahend += val;
         
            if( !SCIPisInfinity(scip, consdata->rhs) )
               rhssubtrahend += val;
         }
         
         /* if only one variable is left, convert the non-linear term into a linear part */
         if( !delterm && consdata->nonlinterms[t]->nvars == 1 )
         {
            SCIP_Bool redundant;
            SCIP_Bool aggregated;

            assert(consdata->nonlinterms[t]->vars == vars);
            assert(vars[0] != NULL);
            SCIP_CALL( addCoefLinear(scip, cons, vars[0], val) );

            if( resultant != NULL )
            { 
               /* @todo: check for the correct aggregation order to throw the resultant away */
               SCIPdebugMessage("aggregating last unfixed variable <%s> in a non-linear part with the and-resultant <%s>\n", SCIPvarGetName(vars[0]), SCIPvarGetName(resultant));
               SCIP_CALL( SCIPaggregateVars(scip, vars[0], resultant, 1.0, -1.0, 0.0, infeasible, &redundant, &aggregated) );
               assert(aggregated || SCIPdoNotAggr(scip));
               assert(!redundant);
               assert(!(*infeasible));
               if( aggregated )
                  ++(*naggrvars);
            }
            delterm = TRUE;
         }

         /* remove redundant term */
         if( delterm )
         {
            SCIPdebugMessage("removing term at pos %d(in the non-linear part), due to fixed and-resultant\n", v);
            SCIP_CALL( delCoefTerm(scip, cons, t) );
         }
         else
            ++t;
      }
      else
         ++t;
   }

   /* update left hand side */
   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      if( SCIPisFeasEQ(scip, lhssubtrahend, consdata->lhs ) )
      {
         SCIP_CALL( chgLhs(scip, cons, 0.0) );
      }
      else
      {
         assert(!SCIPisInfinity(scip, consdata->lhs - lhssubtrahend));
         SCIP_CALL( chgLhs(scip, cons, consdata->lhs - lhssubtrahend) );
      }
   }
   /* update right hand side */
   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      if( SCIPisFeasEQ(scip, rhssubtrahend, consdata->rhs ) )
      {
         SCIP_CALL( chgRhs(scip, cons, 0.0) );
      }
      else
      {
         assert(!SCIPisInfinity(scip, consdata->rhs - rhssubtrahend));
         SCIP_CALL( chgRhs(scip, cons, consdata->rhs - rhssubtrahend) );
      }
   }

#ifndef NDEBUG
   /* check, if all non-linear fixings are applied */
   for( t = 0; t < consdata->nnonlinterms; ++t )
   {
      SCIP_VAR** vars;
      int nvars;

      vars = consdata->nonlinterms[t]->vars;
      nvars = consdata->nonlinterms[t]->nvars;

      for( v = 0; v < nvars; ++v )
         assert(SCIPvarIsActive(vars[v]) || (SCIPvarIsNegated(vars[v]) && SCIPvarIsActive(SCIPvarGetNegationVar(vars[v]))) );
   }
#endif

   return SCIP_OKAY;
}

/** replaces all fixed and aggregated variables by their non-fixed counterparts */
static
SCIP_RETCODE applyFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< pseudoboolean constraint */
   int*                  nfixedvars,         /**< number of fixed variables */
   int*                  naggrvars,          /**< number of aggregated variables */
   SCIP_Bool*            infeasible          /**< pointer to store if infeasibility is detected */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(infeasible != NULL);

   *infeasible = FALSE;

   /* first the non-linear part, because maybe a non-linear part will be 'upgraded' into a linear */
   SCIP_CALL( applyFixingsNonLinear(scip, cons, nfixedvars, naggrvars, infeasible) );

   /* second the linear part */
   SCIP_CALL( applyFixingsLinear(scip, cons) );

   return SCIP_OKAY;
}


/** replaces multiple occurrences of a variable or its negation by a single coefficient */
static
SCIP_RETCODE mergeMultiples(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int*                  ndelconss,          /**< pointer to store number of deleted constraints */
   int*                  nchgcoefs,          /**< pointer to store number of changed coefficients */
   SCIP_Bool*            cutoff              /**< pointer to store whether a fixing leads to a cutoff */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(ndelconss != NULL);
   assert(nchgcoefs != NULL);
   assert(cutoff != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->merged ) 
      return SCIP_OKAY;

   if( consdata->nlinvars <= 1 )
   {
      consdata->merged = TRUE;
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyPseudoboolean)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   
   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrPseudoboolean(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreePseudoboolean)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdataFree(scip, &conshdlrdata);
   
   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}

/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitPseudoboolean)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->allnonlinterms = NULL;
   conshdlrdata->nallnonlinterms = 0;
   conshdlrdata->sallnonlinterms = 0;

   conshdlrdata->hashtable = NULL;
   conshdlrdata->hashtablesize = 0;

   return SCIP_OKAY;
}

/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitPseudoboolean)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* free constraint handler data */
   conshdlrdataExit(scip, conshdlrdata);

   return SCIP_OKAY;
}

/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitprePseudoboolean)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* decompose all pseudo boolean constraints into a "linear" constraint and "and" constraints */
   if( conshdlrdata->decompose )
   {
      int c;

      /* initializes constraint handler data for pseudo boolean constraint handler, if we have constraints */
      if( nconss > 0 )
      {
	 SCIP_CALL( conshdlrdataInit(scip, conshdlrdata) );
      }

      for( c = 0; c < nconss; ++c )
      {
         SCIP_CONS* cons;
         SCIP_CONS* lincons;
         SCIP_CONSDATA* consdata;
         SCIP_VAR** andvars;
         int nandvars;
         int t;

         cons = conss[c];
         assert(cons != NULL);

	 /* only added constraints can be upgraded */
	 if( !SCIPconsIsAdded(cons) )
	    continue;

         consdata = SCIPconsGetData(cons);
         assert(consdata != NULL);

         andvars = NULL;
	 nandvars = 0;

         if( consdata->issoftcons )
         {
            SCIP_VAR* negindvar;
            SCIP_Real* coefs;
            char name[SCIP_MAXSTRLEN];
            SCIP_Real lhs;
            SCIP_Real rhs;
            SCIP_Bool initial;
            SCIP_Bool gotandress;
            int v;
#if USEINDICATOR == FALSE
            SCIP_Real maxact;
            SCIP_Real minact;
            SCIP_Real lb;
            SCIP_Real ub;
#else
            SCIP_CONS* indcons;
#endif         

            assert(consdata->weight != 0);
            assert(consdata->indvar != NULL);

            /* if it is a soft constraint, there should be no integer variable */
            assert(consdata->intvar == NULL);

            /* get negation of indicator variable */
            SCIP_CALL( SCIPgetNegatedVar(scip, consdata->indvar, &negindvar) );
            assert(negindvar != NULL);

            lhs = consdata->lhs;
            rhs = consdata->rhs;

            gotandress = FALSE;

#if USEINDICATOR == FALSE
            maxact = 0.0;
            minact = 0.0;

            coefs = consdata->lincoefs;

            /* adding all linear coefficients up */
            for( v = consdata->nlinvars - 1; v >= 0; --v )
               if( coefs[v] > 0 )
                  maxact += coefs[v];
               else
                  minact += coefs[v];

            coefs = consdata->nonlincoefs;

            /* adding all nonlinear coefficients up */
            for( v = consdata->nnonlinterms - 1; v >= 0; --v )
               if( coefs[v] > 0 )
                  maxact += coefs[v];
               else
                  minact += coefs[v];

            if( SCIPisInfinity(scip, maxact) )
            {
               SCIPwarningMessage(scip, "maxactivity = %g exceed infinity value.\n", maxact);
            }
            if( SCIPisInfinity(scip, -minact) )
            {
               SCIPwarningMessage(scip, "minactivity = %g exceed -infinity value.\n", minact);
            }

            /* @todo check whether it's better to set the initial flag to false */         
            initial = SCIPconsIsInitial(cons); /* FALSE; */

            /* first soft constraints for lhs */
            if( !SCIPisInfinity(scip, -lhs) )
            {
               /* first we are modelling the feasibility of the soft constraint by adding a slack variable */
               /* we ensure that if indvar == 1 => (a^T*x + ub*indvar >= lhs) */
               ub = lhs - minact;

               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_lhs_part1", SCIPconsGetName(cons));

               SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, name, consdata->nlinvars, consdata->linvars, 
                     consdata->lincoefs, lhs, SCIPinfinity(scip),
                     initial, SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

               /* create all and-constraints and get all and-resultants */
               SCIP_CALL( createAnds(scip, cons, conshdlr, &andvars, &nandvars) );
               assert(nandvars == 0 || andvars != NULL);
               assert(nandvars == consdata->nnonlinterms);
               gotandress = TRUE;

               /* add all and-resultants */
               for( t = 0; t < nandvars; ++t )
               {
                  /* add auxiliary variables to linear constraint */
                  SCIP_CALL( SCIPaddCoefLinear(scip, lincons, andvars[t], consdata->nonlincoefs[t]) );
               }
               /* add artificial indicator variable */
               SCIP_CALL( SCIPaddCoefLinear(scip, lincons, consdata->indvar, ub) );

               SCIP_CALL( SCIPaddCons(scip, lincons) );
               SCIPdebugPrintCons(scip, lincons, NULL);
               SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

               /* second we are modelling the implication that if the slack variable is on( negation is off), the constraint
                * is disabled, so only the cost arise if the slack variable is necessary */
               /* indvar == 1 => (a^T*x (+ ub * negindvar) <= lhs - 1) */
               ub = lhs - maxact - 1;

               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_lhs_part2", SCIPconsGetName(cons));

               SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, name, consdata->nlinvars, consdata->linvars, 
                     consdata->lincoefs, -SCIPinfinity(scip), lhs - 1,
                     initial, SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

               /* add all and-resultants */
               for( t = 0; t < nandvars; ++t )
               {
                  /* add auxiliary variables to linear constraint */
                  SCIP_CALL( SCIPaddCoefLinear(scip, lincons, andvars[t], consdata->nonlincoefs[t]) );
               }
               /* add artificial indicator variable */
               SCIP_CALL( SCIPaddCoefLinear(scip, lincons, negindvar , ub) );

               SCIP_CALL( SCIPaddCons(scip, lincons) );
               SCIPdebugPrintCons(scip, lincons, NULL);
               SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
            }

            /* second soft constraints for rhs */
            if( !SCIPisInfinity(scip, rhs) )
            {
               /* first we are modelling the feasibility of the soft-constraint by adding a slack variable */
               /* indvar == 1 => (a^T*x + lb * indvar <= rhs) */
               lb = rhs - maxact;

               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_rhs_part1", SCIPconsGetName(cons));

               SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, name, consdata->nlinvars, consdata->linvars,
                     consdata->lincoefs, -SCIPinfinity(scip), rhs,
                     initial, SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

               if( !gotandress )
               {
                  /* create all and-constraints and get all and-resultants */
                  SCIP_CALL( createAnds(scip, cons, conshdlr, &andvars, &nandvars) );
               }
               assert(nandvars == 0 || andvars != NULL);
               assert(nandvars == consdata->nnonlinterms);

               /* add all and-resultants */
               for( t = 0; t < nandvars; ++t )
               {
                  /* add auxiliary variables to linear constraint */
                  SCIP_CALL( SCIPaddCoefLinear(scip, lincons, andvars[t], consdata->nonlincoefs[t]) );
               }
               /* add artificial indicator variable */
               SCIP_CALL( SCIPaddCoefLinear(scip, lincons, consdata->indvar, lb) );

               SCIP_CALL( SCIPaddCons(scip, lincons) );
               SCIPdebugPrintCons(scip, lincons, NULL);
               SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

               /* second we are modelling the implication that if the slack variable is on( negation is off), the constraint
                * is disabled, so only the cost arise if the slack variable is necessary */
               /* indvar == 1 => (a^T*x (+ lb * negindvar) >= rhs + 1) */
               lb = rhs - minact + 1;

               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_rhs_part2", SCIPconsGetName(cons));

               SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, name, consdata->nlinvars, consdata->linvars,
                     consdata->lincoefs, rhs + 1, SCIPinfinity(scip),
                     initial, SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

               /* add all and-resultants */
               for( t = 0; t < nandvars; ++t )
               {
                  /* add auxiliary variables to linear constraint */
                  SCIP_CALL( SCIPaddCoefLinear(scip, lincons, andvars[t], consdata->nonlincoefs[t]) );
               }
               /* add artificial indicator variable */
               SCIP_CALL( SCIPaddCoefLinear(scip, lincons, negindvar, lb) );

               SCIP_CALL( SCIPaddCons(scip, lincons) );
               SCIPdebugPrintCons(scip, lincons, NULL);
               SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
            }
#else /* with indicator */
            /* @todo check whether it's better to set the initial flag to false */         
            initial = SCIPconsIsInitial(cons); /* FALSE; */

            if( !SCIPisInfinity(scip, rhs) )
            {
               /* first we are modelling the implication that if the negation of the indicator variable is on, the constraint
                * is enabled */
               /* indvar == 0 => a^T*x <= rhs */

               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_rhs_ind", SCIPconsGetName(cons));

               SCIP_CALL( SCIPcreateConsIndicator(scip, &indcons, name, negindvar, consdata->nlinvars, 
                     consdata->linvars, consdata->lincoefs, rhs,
                     initial, SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), 
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

               /* create all and-constraints and get all and-resultants */
               SCIP_CALL( createAnds(scip, cons, conshdlr, &andvars, &nandvars) );
               assert(nandvars == 0 || andvars != NULL);
               assert(nandvars == consdata->nnonlinterms);
               gotandress = TRUE;

               /* get linear constraint corresponding to the indicator constraint */
               lincons = SCIPgetLinearConsIndicator(indcons);
               assert(lincons != NULL);

               /* add all and-resultants */
               for( t = 0; t < nandvars; ++t )
               {
                  /* add auxiliary variables to linear constraint */
                  SCIP_CALL( SCIPaddCoefLinear(scip, lincons, andvars[t], consdata->nonlincoefs[t]) );
               }

               SCIP_CALL( SCIPaddCons(scip, indcons) );
               SCIPdebugPrintCons(scip, indcons, NULL);
               SCIP_CALL( SCIPreleaseCons(scip, &indcons) );
            }

            if( !SCIPisInfinity(scip, -lhs) )
            {
               /* second we are modelling the implication that if the negation of the indicator variable is on, the constraint
                * is enabled */
               /* change the a^T*x >= lhs to -a^Tx<= -lhs, for indicator constraint */

               /* get temporary memory */
               SCIP_CALL( SCIPduplicateBufferArray(scip, &coefs, consdata->lincoefs, consdata->nlinvars) );

               for( v = consdata->nlinvars - 1; v >= 0; --v )
                  coefs[v] *= -1;

               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_lhs_ind", SCIPconsGetName(cons));

               SCIP_CALL( SCIPcreateConsIndicator(scip, &indcons, name, negindvar, consdata->nlinvars,
                     consdata->linvars, coefs, -lhs,
                     initial, SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

               if( !gotandress )
               {
                  /* create all and-constraints and get all and-resultants */
                  SCIP_CALL( createAnds(scip, cons, conshdlr, &andvars, &nandvars) );
               }
               assert(nandvars == 0 || andvars != NULL);
               assert(nandvars == consdata->nnonlinterms);

               /* get linear constraint corresponding to the indicator constraint */
               lincons = SCIPgetLinearConsIndicator(indcons);
               assert(lincons != NULL);

               /* add all and-resultants */
               for( t = 0; t < nandvars; ++t )
               {
                  /* add auxiliary variables to linear constraint */
                  SCIP_CALL( SCIPaddCoefLinear(scip, lincons, andvars[t], consdata->nonlincoefs[t]) ); /*lint !e613 */
               }

               SCIP_CALL( SCIPaddCons(scip, indcons) );
               SCIPdebugPrintCons(scip, indcons, NULL);
               SCIP_CALL( SCIPreleaseCons(scip, &indcons) );

               /* free temporary memory */
               SCIPfreeBufferArray(scip, &coefs);
            }
#endif
            /* free temporary memory */
            SCIPfreeBufferArrayNull(scip, &andvars);
         }
         /* no soft constraint */
         else
         {
            /* create linear constraint */
            SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, SCIPconsGetName(cons), 
                  consdata->nlinvars, consdata->linvars, consdata->lincoefs, consdata->lhs, consdata->rhs, 
                  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), 
                  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

            /* create all and-constraints and get all and-resultants */
            SCIP_CALL( createAnds(scip, cons, conshdlr, &andvars, &nandvars) );
            assert(nandvars == 0 || andvars != NULL);
            assert(nandvars == consdata->nnonlinterms);

            for( t = 0; t < nandvars; ++t )
            {
               /* add auxiliary variables to linear constraint */
               SCIP_CALL( SCIPaddCoefLinear(scip, lincons, andvars[t], consdata->nonlincoefs[t]) );
            }

            if( consdata->intvar != NULL )
            {
               /* add auxiliary integer variables to linear constraint */
               SCIP_CALL( SCIPaddCoefLinear(scip, lincons, consdata->intvar, -1.0) );
            }

            SCIP_CALL( SCIPaddCons(scip, lincons) );
            SCIPdebugPrintCons(scip, lincons, NULL);
            SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

            /* free temporary memory */
            SCIPfreeBufferArrayNull(scip, &andvars);
         }
         /* remove pseudo boolean constraint */
         SCIP_CALL( SCIPdelCons(scip, cons) );
      }
   }

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeletePseudoboolean)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* free pseudo boolean constraint */
   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */ 
static
SCIP_DECL_CONSTRANS(consTransPseudoboolean)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* create linear constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, &targetdata, 
         sourcedata->linvars, sourcedata->lincoefs, sourcedata->nlinvars,
         sourcedata->nonlinterms, sourcedata->nonlincoefs, sourcedata->nnonlinterms,
         sourcedata->indvar, sourcedata->weight, sourcedata->issoftcons, sourcedata->intvar, sourcedata->lhs, sourcedata->rhs) );
   
   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpPseudoboolean)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int c;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   violated = FALSE;
   
   /* check all pseudo boolean constraints for feasibility */
   for( c = 0; c < nconss && !violated; ++c )
   {
      SCIP_CALL( checkCons(scip, conss[c], NULL, FALSE, &violated) );
   }

   if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;
   
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsPseudoboolean)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int c;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   violated = FALSE;
   
   /* check all pseudo boolean constraints for feasibility */
   for( c = 0; c < nconss && !violated; ++c )
   {
      SCIP_CALL( checkCons(scip, conss[c], NULL, FALSE, &violated) );
   }

   if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;
   
   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckPseudoboolean)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int c;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(sol != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   violated = FALSE;

   /* check all pseudo boolean constraints for feasibility */
   for( c = 0; c < nconss && !violated; ++c )
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* does not need to check indictor constraints which are switched off */
      if( consdata->issoftcons )
      {
         assert(consdata->indvar != NULL);
         if( SCIPisEQ(scip, SCIPgetSolVal(scip, sol, consdata->indvar), 1.0) )
            continue;
      }
      SCIP_CALL( checkCons(scip, conss[c], sol, printreason, &violated) );
   }

   if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;
   
   return SCIP_OKAY;
}

/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolPseudoboolean)
{  /*lint --e{715}*/
   SCIP_Bool cutoff;
   int oldnfixedvars;
   int oldnaggrvars;
   int oldnchgbds;
   int oldndelconss;
   int oldnupgdconss;
   int oldnchgcoefs;
   int oldnchgsides;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* remember old preprocessing counters */
   cutoff = FALSE;
   oldnfixedvars = *nfixedvars;
   oldnaggrvars = *naggrvars;
   oldnchgbds = *nchgbds;
   oldndelconss = *ndelconss;
   oldnupgdconss = *nupgdconss;
   oldnchgcoefs = *nchgcoefs;
   oldnchgsides = *nchgsides;

#if 0
   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
#endif

   for( c = 0; c < nconss && !cutoff && !SCIPisStopped(scip); ++c )
   {
      SCIP_CONS* cons;
#ifndef NDEBUG
      SCIP_CONSDATA* consdata;
#endif
      SCIP_Bool infeasible;

      infeasible = FALSE;

      cons = conss[c];
      assert(SCIPconsIsActive(cons));
#ifndef NDEBUG
      consdata = SCIPconsGetData(cons);
#endif
      assert(consdata != NULL);

      /* incorporate fixings and aggregations in constraint */
      SCIP_CALL( applyFixings(scip, cons, nfixedvars, naggrvars, &infeasible) );

      if( infeasible )
      {
         SCIPdebugMessage(" -> cutoff due to infeasible fixing(s)\n");
         cutoff = TRUE;
         break;
      }

      assert(consdata->removedfixings);

      /* if aggregated variables have been replaced, multiple entries of the same variable are possible and we have
       * to clean up the constraint
       */
      SCIP_CALL( mergeMultiples(scip, cons, ndelconss, nchgcoefs, &cutoff) );

      SCIPdebugMessage("after merging:\n");
      SCIPdebugPrintCons(scip, cons, NULL);

      /* we can only presolve pseudoboolean constraints, that are not modifiable */
      if( SCIPconsIsModifiable(cons) )
         continue;
   }

   /* return the correct result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( *nfixedvars > oldnfixedvars || *naggrvars > oldnaggrvars || *nchgbds > oldnchgbds || *ndelconss > oldndelconss
      || *nupgdconss > oldnupgdconss || *nchgcoefs > oldnchgcoefs || *nchgsides > oldnchgsides )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockPseudoboolean)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool haslhs;
   SCIP_Bool hasrhs;
   int v;
   int t;

   assert(scip != NULL);
   assert(cons != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   haslhs = !SCIPisInfinity(scip, -consdata->lhs);
   hasrhs = !SCIPisInfinity(scip, consdata->rhs);

   /* update rounding locks of every single variable corresponding to the linear terms */
   for( v = 0; v < consdata->nlinvars; ++v )
   {
      if( SCIPisPositive(scip, consdata->lincoefs[v]) )
      {
         if( haslhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[v], nlockspos, nlocksneg) );
         }
         if( hasrhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[v], nlocksneg, nlockspos) );
         }
      }
      else
      {
         if( haslhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[v], nlocksneg, nlockspos) );
         }
         if( hasrhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[v], nlockspos, nlocksneg) );
         }
      }
   }

   /* update rounding locks of every single variable corresponding to the nonlinear terms */
   for( t = 0; t < consdata->nnonlinterms; ++t )
   {
      TERM* term;

      term = consdata->nonlinterms[t];

      if( SCIPisPositive(scip, consdata->nonlincoefs[t]) )
      {
         if( haslhs )
         {
            for( v = 0; v < term->nvars; ++v )
            {
               SCIP_CALL( SCIPaddVarLocks(scip, term->vars[v], nlockspos, nlocksneg) );
            }
         }
         if( hasrhs )
         {
            for( v = 0; v < term->nvars; ++v )
            {
               SCIP_CALL( SCIPaddVarLocks(scip, term->vars[v], nlocksneg, nlockspos) );
            }
         }
      }
      else
      {
         if( haslhs )
         {
            for( v = 0; v < term->nvars; ++v )
            {
               SCIP_CALL( SCIPaddVarLocks(scip, term->vars[v], nlocksneg, nlockspos) );
            }
         }
         if( hasrhs )
         {
            for( v = 0; v < term->nvars; ++v )
            {
               SCIP_CALL( SCIPaddVarLocks(scip, term->vars[v], nlockspos, nlocksneg) );
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyPseudoboolean)
{  /*lint --e{715}*/
   const char* consname;

   assert(scip != NULL);
   assert(sourcescip != NULL);
   assert(sourcecons != NULL);

   if( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   SCIP_CALL( copyConsPseudoboolean(scip, cons, sourcescip, sourcecons, consname, varmap, consmap, 
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, global, valid) );
   assert(cons != NULL || *valid == FALSE);

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for pseudoboolean constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrPseudoboolean(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create pseudoboolean constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpPseudoboolean, consEnfopsPseudoboolean, consCheckPseudoboolean, consLockPseudoboolean,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyPseudoboolean, consCopyPseudoboolean) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeletePseudoboolean) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitPseudoboolean) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreePseudoboolean) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitPseudoboolean) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitprePseudoboolean) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolPseudoboolean, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYPRESOL) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransPseudoboolean) );

   /* add pseudoboolean constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/decompose",
         "decompose the pseudo boolean constraint into a \"linear\" constraint \"and\" constraints",
         &conshdlrdata->decompose, TRUE, DEFAULT_DECOMPOSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/nlcseparate", "if decomposed, should the nonlinear constraints be separated during LP processing?",
         NULL, TRUE, DEFAULT_SEPARATENONLINEAR, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/nlcpropagate", "if decomposed, should the nonlinear constraints be propagated during node processing?",
         NULL, TRUE, DEFAULT_PROPAGATENONLINEAR, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/nlcremovable", "if decomposed, should the nonlinear constraints be removable?",
         NULL, TRUE, DEFAULT_REMOVABLENONLINEAR, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a pseudoboolean constraint 
 *
 *  @note linear and nonlinear terms can be added using SCIPaddCoefPseudoboolean() and SCIPaddTermPseudoboolean(),
 *        respectively
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR**            linvars,            /**< variables of the linear part, or NULL */
   int                   nlinvars,           /**< number of variables of the linear part */
   SCIP_Real*            linvals,            /**< coefficients of linear part, or NULL */
   SCIP_VAR***           terms,              /**< nonlinear terms of variables, or NULL */
   int                   nterms,             /**< number of terms of variables of nonlinear term */
   int*                  ntermvars,          /**< number of variables in nonlinear terms, or NULL */
   SCIP_Real*            termvals,           /**< coefficients of nonlinear parts, or NULL */
   SCIP_VAR*             indvar,             /**< indicator variable if it's a soft constraint, or NULL */
   SCIP_Real             weight,             /**< weight of the soft constraint, if it is one */
   SCIP_Bool             issoftcons,         /**< is this a soft constraint */
   SCIP_VAR*             intvar,             /**< a artificial variable which was added only for the objective function,
                                              *   if this variable is not NULL this constraint (without this integer
                                              *   variable) describes the objective function */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
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
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   TERM** pbterms;
   SCIP_Real* pbcoefs;
   int npbterms;
   int t;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nlinvars == 0 || (linvars != NULL && linvals != NULL));
   assert(nterms == 0 || (terms != NULL && termvals != NULL && ntermvars != NULL));
   assert(issoftcons == (indvar != NULL));

   /* find the pseudoboolean constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("pseudo boolean constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create non-linear term data */
   SCIP_CALL( createTerms(scip, terms, termvals, nterms, ntermvars, &pbterms, &pbcoefs, &npbterms) );

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, linvars, linvals, nlinvars, pbterms, pbcoefs, npbterms, indvar, weight, issoftcons, intvar, lhs, rhs) );
   assert(consdata != NULL);
   
   /* delete temporary terms */
   for( t = npbterms - 1; t >= 0; --t )
   {
      assert(pbterms[t] != NULL);
      assert(pbterms[t]->vars != NULL);

      SCIPfreeBlockMemoryArray(scip, &(pbterms[t]->vars), pbterms[t]->nvars);
      SCIPfreeBlockMemory(scip, &(pbterms[t]) ); /*lint !e866 */
   }
   /* clear temporary created array */
   SCIPfreeBlockMemoryArrayNull(scip, &pbcoefs, npbterms);
   SCIPfreeBlockMemoryArrayNull(scip, &pbterms, npbterms);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );
   
   return SCIP_OKAY;
}

/** creates and captures a pseudoboolean constraint
 *  in its most basic variant, i. e., with all constraint flags set to their default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR**            linvars,            /**< variables of the linear part, or NULL */
   int                   nlinvars,           /**< number of variables of the linear part */
   SCIP_Real*            linvals,            /**< coefficients of linear part, or NULL */
   SCIP_VAR***           terms,              /**< nonlinear terms of variables, or NULL */
   int                   nterms,             /**< number of terms of variables of nonlinear term */
   int*                  ntermvars,          /**< number of variables in nonlinear terms, or NULL */
   SCIP_Real*            termvals,           /**< coefficients of nonlinear parts, or NULL */
   SCIP_VAR*             indvar,             /**< indicator variable if it's a soft constraint, or NULL */
   SCIP_Real             weight,             /**< weight of the soft constraint, if it is one */
   SCIP_Bool             issoftcons,         /**< is this a soft constraint */
   SCIP_VAR*             intvar,             /**< a artificial variable which was added only for the objective function,
                                              *   if this variable is not NULL this constraint (without this integer
                                              *   variable) describes the objective function */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   )
{
   SCIP_CALL( SCIPcreateConsPseudoboolean(scip, cons, name, linvars, nlinvars, linvals,
         terms, nterms, ntermvars, termvals, indvar, weight, issoftcons, intvar, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** adds a variable to the pseudo boolean constraint (if it is not zero) */
SCIP_RETCODE SCIPaddCoefPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( addCoefLinear(scip, cons, var, val) );

   return SCIP_OKAY;
}


/** adds nonlinear term to pseudo boolean constraint (if it is not zero) */
SCIP_RETCODE SCIPaddTermPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR**            vars,               /**< variables of the nonlinear term */
   int                   nvars,              /**< number of variables of the nonlinear term */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(nvars == 0 || vars != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( addCoefTerm(scip, cons, vars, nvars, val) );

   return SCIP_OKAY;
}

/** gets left hand side of pseudoboolean constraint */
SCIP_Real SCIPgetLhsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->lhs;
}

/** gets right hand side of pseudoboolean constraint */
SCIP_Real SCIPgetRhsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->rhs;
}

/** gets indicator variable of pseudoboolean constraint, or NULL if there is no */
SCIP_VAR* SCIPgetIndVarPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->indvar;
}

/** gets linear variables of pseudoboolean constraint */
SCIP_VAR** SCIPgetLinearVarsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->linvars;
}

/** gets linear coefficients of pseudoboolean constraint */
SCIP_Real* SCIPgetLinearValsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->lincoefs;
}


/** gets number of linear variables of pseudoboolean constraint */
int SCIPgetNLinearVarsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nlinvars;
}


/** gets array with non-linear term variables and corresponding number of variable in each non-linear term iff
 *  'termvarssize' (size of array for all non-linear term variables) is big enough, or will return the needed size for
 *  the array for all non-linear variables in termvarssize otherwise, of a given pseudoboolean constraint 
 */
SCIP_RETCODE SCIPgetTermVarsDataPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR***           termvars,           /**< array for all arrays of term variables */
   int*                  ntermvars,          /**< array for number of variable in each term */
   int*                  termvarssize        /**< pointer to store size of array which is needed for all term variables */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( *termvarssize > 0 && consdata->nnonlinterms <= *termvarssize )
   {
      int t;

      assert(termvars != NULL);
      assert(ntermvars != NULL);

      for( t = 0; t < consdata->nnonlinterms; ++t )
      {
         termvars[t] = consdata->nonlinterms[t]->vars;
         ntermvars[t] = consdata->nonlinterms[t]->nvars;
      }
   }

   *termvarssize = consdata->nnonlinterms;

   return SCIP_OKAY;
}

/** gets non-linear coefficients of pseudoboolean constraint */
SCIP_Real* SCIPgetTermValsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nonlincoefs;
}

/** gets number of non-linear coefficients of pseudoboolean constraint */
int SCIPgetNTermValsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nnonlinterms;
}

/** changes left hand side of pseudoboolean constraint */
SCIP_RETCODE SCIPchgLhsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( chgLhs(scip, cons, lhs) );

   return SCIP_OKAY;
}

/** changes right hand side of pseudoboolean constraint */
SCIP_RETCODE SCIPchgRhsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not pseudo boolean\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( chgRhs(scip, cons, rhs) );

   return SCIP_OKAY;
}

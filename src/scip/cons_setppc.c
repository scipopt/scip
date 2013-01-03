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

/**@file   cons_setppc.c
 * @brief  Constraint handler for the set partitioning / packing / covering constraints \f$1^T x\ \{=, \le, \ge\}\ 1\f$.
 * @author Tobias Achterberg
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <ctype.h>

#include "scip/cons_setppc.h"
#include "scip/cons_linear.h"
#include "scip/pub_misc.h"


#define CONSHDLR_NAME          "setppc"
#define CONSHDLR_DESC          "set partitioning / packing / covering constraints"
#define CONSHDLR_SEPAPRIORITY   +700000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   -700000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -700000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             0 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP

#define LINCONSUPGD_PRIORITY    +700000 /**< priority of the constraint handler for upgrading of linear constraints */

#define EVENTHDLR_NAME         "setppc"
#define EVENTHDLR_DESC         "bound change event handler for set partitioning / packing / covering constraints"

#define CONFLICTHDLR_NAME      "setppc"
#define CONFLICTHDLR_DESC      "conflict handler creating set covering constraints"
#define CONFLICTHDLR_PRIORITY  LINCONSUPGD_PRIORITY

#define DEFAULT_PRESOLPAIRWISE     TRUE /**< should pairwise constraint comparison be performed in presolving? */

#define HASHSIZE_SETPPCCONS      131101 /**< minimal size of hash table in setppc constraint tables */
#define DEFAULT_PRESOLUSEHASHING   TRUE /**< should hash table be used for detecting redundant constraints in advance */
#define NMINCOMPARISONS          200000 /**< number for minimal pairwise presolving comparisons */
#define MINGAINPERNMINCOMPARISONS 1e-06 /**< minimal gain per minimal pairwise presolving comparisons to repeat pairwise comparison round */

/*#define VARUSES*/  /* activate variable usage counting, that is necessary for LP and pseudo branching */
/*#define BRANCHLP*/ /* BRANCHLP is only useful if the ENFOPRIORITY is set to a positive value */
#ifdef BRANCHLP
#define MINBRANCHWEIGHT             0.3 /**< minimum weight of both sets in binary set branching */
#define MAXBRANCHWEIGHT             0.7 /**< maximum weight of both sets in binary set branching */
#endif
#define DEFAULT_NPSEUDOBRANCHES       2 /**< number of children created in pseudo branching (0: disable branching) */
#define DEFAULT_DUALPRESOLVING     TRUE /**< should dual presolving steps be performed? */

#define DEFAULT_CLIQUELIFTING     FALSE /**< should we try to lift variables into other clique constraints, fix
					 *   variables, aggregate them, and also shrink the amount of variables in
					 *   clique constraints
					 */
#define DEFAULT_ADDVARIABLESASCLIQUES FALSE/**< should we try to generate extra clique constraint out of all binary
                                            *   variables to hopefully fasten the detection of redundant clique
                                            *   constraints */
#define DEFAULT_CLIQUESHRINKING    TRUE /**< should we try to shrink the number of variables in a clique constraints, by
					 *   replacing more than one variable by only one
					 */

/* @todo maybe use event SCIP_EVENTTYPE_VARUNLOCKED to decide for another dual-presolving run on a constraint */

/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
#ifdef VARUSES
   SCIP_INTARRAY*        varuses;            /**< number of times a var is used in the active setppc constraints */
#endif
   int                   npseudobranches;    /**< number of children created in pseudo branching (0 to disable branching) */
   int                   noldfixedvars;      /**< number of fixed variables after last clique lifting run */
   int                   noldimpls;          /**< number of implication before last clique lifting run */
   int                   noldcliques;        /**< number of cliques before last clique lifting run */
   int                   noldupgrs;          /**< number of setppc constraints since the last clique lifting run */
   int                   nclqpresolve;       /**< number of setppc clique lifting runs */
   SCIP_Bool             updatedsetppctype;  /**< remember whether we upgraded a constraint type */
   SCIP_Bool             cliquelifting;      /**< should we perform the clique lifting procedure */
   SCIP_Bool             enablecliquelifting;/**< check whether we have enough changes to run the lifting procedure again */
   SCIP_Bool             cliqueshrinking;    /**< should we try to shrink the number of variables in a clique
					      *   constraints, by replacing more than one variable by only one
					      */
   SCIP_Bool             addvariablesascliques;/**< should we try to generate extra clique constraint out of all binary
                                                *   variables to hopefully fasten the detection of redundant clique
                                                *   constraints */
   unsigned int          randseed;           /**< start value for randomly permutate binary variable array to create
                                              *   extra cliques, iff addvariablesascliques is set to TRUE */
   SCIP_Bool             presolpairwise;     /**< should pairwise constraint comparison be performed in presolving? */
   SCIP_Bool             presolusehashing;   /**< should hash table be used for detecting redundant constraints in advance */
   SCIP_Bool             dualpresolving;     /**< should dual presolving steps be performed? */
};

/** constraint data for set partitioning / packing / covering constraints */
struct SCIP_ConsData
{
   SCIP_Longint          signature;          /**< bit signature of vars array */
   SCIP_ROW*             row;                /**< LP row, if constraint is already stored in LP row format */
   SCIP_VAR**            vars;               /**< variables of the constraint */
   int                   varssize;           /**< size of vars array */
   int                   nvars;              /**< number of variables in the constraint */
   int                   nfixedzeros;        /**< current number of variables fixed to zero in the constraint */
   int                   nfixedones;         /**< current number of variables fixed to one in the constraint */
   unsigned int          setppctype:2;       /**< type of constraint: set partitioning, packing or covering */
   unsigned int          sorted:1;           /**< are the constraint's variables sorted? */
   unsigned int          cliqueadded:1;      /**< was the set partitioning / packing constraint already added as clique? */
   unsigned int          validsignature:1;   /**< is the bit signature valid? */
   unsigned int          changed:1;          /**< was constraint changed since last redundancy round in preprocessing? */
   unsigned int          varsdeleted:1;      /**< were variables deleted after last cleanup? */
   unsigned int          merged:1;           /**< are the constraint's equal/negated variables already merged? */
   unsigned int          propagated:1;       /**< does this constraint need to be propagated? */
};




/*
 * Local methods
 */

/** compares two active constraints of type set partitioning or set packing such that a "-1" is return if
 *    1. the first constraint is a set partitioning constraint and the second is a set packing or
 *    2. both constraints are set partitioning constraints and the second has more! variables than the first or
 *    3. both constraints are set packing constraints and the second has less! variables than the first
 *  a "0" is return if
 *    1. both constraint are of the same type and have the amount of variables or
 *  and a "1" is returned otherwise
 */
static
int setppcCompare(
   SCIP_CONS*const       cons1,              /**< first problem variable */
   SCIP_CONS*const       cons2               /**< second problem variable */
   )
{
   SCIP_CONSDATA* consdata1;
   SCIP_CONSDATA* consdata2;

   assert(cons1 != NULL);
   assert(cons2 != NULL);
   assert(SCIPconsIsActive(cons1));
   assert(SCIPconsIsActive(cons2));

   /* the partitioning type should be the smallest value and the packing the second smallest */
   assert(SCIP_SETPPCTYPE_PARTITIONING < SCIP_SETPPCTYPE_PACKING);

   consdata1 = SCIPconsGetData(cons1);
   assert(consdata1 != NULL);
   assert(consdata1->setppctype != SCIP_SETPPCTYPE_COVERING); /*lint !e641*/
   consdata2 = SCIPconsGetData(cons2);
   assert(consdata2 != NULL);
   assert(consdata2->setppctype != SCIP_SETPPCTYPE_COVERING); /*lint !e641*/

   if( consdata1->setppctype < consdata2->setppctype ||
      (consdata1->setppctype == SCIP_SETPPCTYPE_PARTITIONING && consdata1->nvars < consdata2->nvars) || /*lint !e641*/
      (consdata2->setppctype == SCIP_SETPPCTYPE_PACKING && consdata1->nvars > consdata2->nvars) ) /*lint !e641*/
      return -1;
   else if( (consdata1->setppctype == consdata2->setppctype && consdata1->nvars == consdata2->nvars) ) /*lint !e641*/
      return 0;
   else
   {
      assert(consdata1->setppctype > consdata2->setppctype || (consdata1->setppctype == SCIP_SETPPCTYPE_PARTITIONING && consdata1->setppctype == consdata2->setppctype && consdata1->nvars > consdata2->nvars) || (consdata1->setppctype == SCIP_SETPPCTYPE_PACKING && consdata1->setppctype == consdata2->setppctype && consdata1->nvars < consdata2->nvars)); /*lint !e641*/
      return +1;
   }
}

/** sort constraints first after type (partitioning before packing before covering) and second after number of
 *  variables such that the partitioning constraints have increasing number of variables and the packing constraints
 *  have decreasing number of variables */
static
SCIP_DECL_SORTPTRCOMP(setppcConssSort)
{
   return setppcCompare((SCIP_CONS*)elem1, (SCIP_CONS*)elem2);
}


/** installs rounding locks for the given variable in the given setppc constraint */
static
SCIP_RETCODE lockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< setppc constraint */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   switch( consdata->setppctype )
   {
   case SCIP_SETPPCTYPE_PARTITIONING:
      /* rounding in both directions may violate the constraint */
      SCIP_CALL( SCIPlockVarCons(scip, var, cons, TRUE, TRUE) );
      break;
   case SCIP_SETPPCTYPE_PACKING:
      /* rounding up may violate the constraint */
      SCIP_CALL( SCIPlockVarCons(scip, var, cons, FALSE, TRUE) );
      break;
   case SCIP_SETPPCTYPE_COVERING:
      /* rounding down may violate the constraint */
      SCIP_CALL( SCIPlockVarCons(scip, var, cons, TRUE, FALSE) );
      break;
   default:
      SCIPerrorMessage("unknown setppc type\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** removes rounding locks for the given variable in the given setppc constraint */
static
SCIP_RETCODE unlockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< setppc constraint */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   switch( consdata->setppctype )
   {
   case SCIP_SETPPCTYPE_PARTITIONING:
      /* rounding in both directions may violate the constraint */
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons, TRUE, TRUE) );
      break;
   case SCIP_SETPPCTYPE_PACKING:
      /* rounding up may violate the constraint */
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons, FALSE, TRUE) );
      break;
   case SCIP_SETPPCTYPE_COVERING:
      /* rounding down may violate the constraint */
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons, TRUE, FALSE) );
      break;
   default:
      SCIPerrorMessage("unknown setppc type\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** creates constraint handler data for set partitioning / packing / covering constraint handler */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata,       /**< pointer to store the constraint handler data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPallocMemory(scip, conshdlrdata) );
#ifdef VARUSES
   SCIP_CALL( SCIPcreateIntarray(scip, &(*conshdlrdata)->varuses) );
#endif
   (*conshdlrdata)->npseudobranches = DEFAULT_NPSEUDOBRANCHES;

   /* set event handler for bound change events */
   (*conshdlrdata)->eventhdlr = eventhdlr;

   return SCIP_OKAY;
}

/** frees constraint handler data for set partitioning / packing / covering constraint handler */
static
SCIP_RETCODE conshdlrdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

#ifdef VARUSES
   SCIP_CALL( SCIPfreeIntarray(scip, &(*conshdlrdata)->varuses) );
#endif
   SCIPfreeMemory(scip, conshdlrdata);

   return SCIP_OKAY;
}

#ifdef VARUSES
/** adds the given value to the usage counter of the given variable */
static
SCIP_RETCODE conshdlrdataAddVaruses(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_VAR*             var,                /**< variable to increase usage counter for */
   int                   addval              /**< value to add to the usage counter */
   )
{
   SCIP_INTARRAY* varuses;

   assert(conshdlrdata != NULL);
   assert(var != NULL);

   varuses = conshdlrdata->varuses;
   assert(varuses != NULL);

   /* if the variable is the negation of a problem variable, count the varuses in the problem variable */
   if( SCIPvarIsNegated(var) )
   {
      SCIP_VAR* negvar;
      int varindex;

      /* move the varuses value of the negated variable to the active problem variable */
      varindex = SCIPvarGetIndex(var);
      addval += SCIPgetIntarrayVal(scip, varuses, varindex);
      SCIP_CALL( SCIPsetIntarrayVal(scip, varuses, varindex, 0) );
      SCIP_CALL( SCIPgetNegatedVar(scip, var, &negvar) );
      var = negvar;
   }

   /* increase varuses counter */
   SCIP_CALL( SCIPincIntarrayVal(scip, varuses, SCIPvarGetIndex(var), addval) );

   SCIPdebugMessage("added %d to varuses of <%s>: %d\n",
      addval, SCIPvarGetName(var), SCIPgetIntarrayVal(scip, varuses, SCIPvarGetIndex(var)));

   return SCIP_OKAY;
}

/** increases the usage counter of the given variable */
static
SCIP_RETCODE conshdlrdataIncVaruses(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_VAR*             var                 /**< variable to increase usage counter for */
   )
{
   assert(conshdlrdata != NULL);

   SCIPdebugMessage("increasing varuses of <%s>: %d\n",
      SCIPvarGetName(var), SCIPgetIntarrayVal(scip, conshdlrdata->varuses, SCIPvarGetIndex(var)));

   SCIP_CALL( conshdlrdataAddVaruses(scip, conshdlrdata, var, +1) );

   return SCIP_OKAY;
}

/** decreases the usage counter of the given variable */
static
SCIP_RETCODE conshdlrdataDecVaruses(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_VAR*             var                 /**< variable to increase usage counter for */
   )
{
   assert(conshdlrdata != NULL);

   SCIPdebugMessage("decreasing varuses of <%s>: %d\n",
      SCIPvarGetName(var), SCIPgetIntarrayVal(scip, conshdlrdata->varuses, SCIPvarGetIndex(var)));

   SCIP_CALL( conshdlrdataAddVaruses(scip, conshdlrdata, var, -1) );

   return SCIP_OKAY;
}

/** increases the usage counter of all variable in the constraint */
static
SCIP_RETCODE consdataIncVaruses(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONSDATA*        consdata            /**< setppc constraint data */
   )
{
   int v;

   assert(consdata != NULL);

   for( v = 0; v < consdata->nvars; ++v )
   {
      SCIP_CALL( conshdlrdataIncVaruses(scip, conshdlrdata, consdata->vars[v]) );
   }

   return SCIP_OKAY;
}

/** decreases the usage counter of all variable in the constraint */
static
SCIP_RETCODE consdataDecVaruses(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONSDATA*        consdata            /**< setppc constraint data */
   )
{
   int v;

   assert(consdata != NULL);

   for( v = 0; v < consdata->nvars; ++v )
   {
      SCIP_CALL( conshdlrdataDecVaruses(scip, conshdlrdata, consdata->vars[v]) );
   }

   return SCIP_OKAY;
}
#endif

/** ensures, that the vars array can store at least num entries */
static
SCIP_RETCODE consdataEnsureVarsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< setppc constraint data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(consdata != NULL);
   assert(consdata->nvars <= consdata->varssize);

   if( num > consdata->varssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, consdata->varssize, newsize) );
      consdata->varssize = newsize;
   }
   assert(num <= consdata->varssize);

   return SCIP_OKAY;
}

/** creates a set partitioning / packing / covering constraint data object */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store the set partitioning / packing / covering constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< variables of the constraint */
   SCIP_SETPPCTYPE       setppctype          /**< type of constraint: set partitioning, packing, or covering constraint */
   )
{
   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   (*consdata)->signature = 0;
   (*consdata)->row = NULL;
   if( nvars > 0 )
   {
      int v;

      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
      (*consdata)->varssize = nvars;
      (*consdata)->nvars = nvars;

      if( SCIPisTransformed(scip) )
      {
         /* get transformed variables */
         SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );
      }


      /* capture variables */
      for( v = 0; v < (*consdata)->nvars; v++ )
      {
         assert((*consdata)->vars[v] != NULL);
         SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->vars[v]) );
      }

   }
   else
   {
      (*consdata)->vars = NULL;
      (*consdata)->varssize = 0;
      (*consdata)->nvars = 0;
   }
   (*consdata)->nfixedzeros = 0;
   (*consdata)->nfixedones = 0;
   (*consdata)->setppctype = setppctype; /*lint !e641*/
   (*consdata)->sorted = (nvars <= 1);
   (*consdata)->cliqueadded = FALSE;
   (*consdata)->validsignature = FALSE;
   (*consdata)->changed = TRUE;
   (*consdata)->varsdeleted = FALSE;
   (*consdata)->merged = FALSE;
   (*consdata)->propagated = FALSE;

   return SCIP_OKAY;
}

/** creates a transformed set partitioning / packing / covering constraint data object */
static
SCIP_RETCODE consdataCreateTransformed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store the set partitioning / packing / covering constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< variables of the constraint */
   SCIP_SETPPCTYPE       setppctype          /**< type of constraint: set partitioning, packing, or covering constraint */
   )
{
   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, consdata, nvars, vars, setppctype) );

   /* transform the variables */
   SCIP_CALL( SCIPgetTransformedVars(scip, nvars, (*consdata)->vars, (*consdata)->vars) );

   return SCIP_OKAY;
}

/** frees a set partitioning / packing / covering constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to store the set partitioning / packing / covering constraint */
   )
{
   int v;

   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* release the row */
   if( (*consdata)->row != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->row) );
   }

   /* release variables */
   for( v = 0; v < (*consdata)->nvars; v++ )
   {
      assert((*consdata)->vars[v] != NULL);
      SCIP_CALL( SCIPreleaseVar(scip, &((*consdata)->vars[v])) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** prints set partitioning / packing / covering constraint to file stream */
static
SCIP_RETCODE consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< set partitioning / packing / covering constraint data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert(consdata != NULL);

   /* print coefficients */
   if( consdata->nvars == 0 )
      SCIPinfoMessage(scip, file, "0 ");

   /* write linear sum */
   SCIP_CALL( SCIPwriteVarsLinearsum(scip, file, consdata->vars, NULL, consdata->nvars, TRUE) );
   
   /* print right hand side */
   switch( consdata->setppctype )
   {
   case SCIP_SETPPCTYPE_PARTITIONING:
      SCIPinfoMessage(scip, file, " == 1");
      break;
   case SCIP_SETPPCTYPE_PACKING:
      SCIPinfoMessage(scip, file, " <= 1");
      break;
   case SCIP_SETPPCTYPE_COVERING:
      SCIPinfoMessage(scip, file, " >= 1");
      break;
   default:
      SCIPerrorMessage("unknown setppc type\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** returns the signature bitmask for the given variable */
static
SCIP_Longint getVarSignature(
   SCIP_VAR*             var                 /**< variable */
   )
{
   int sigidx;

   sigidx = SCIPvarGetIndex(var) % (int)(8*sizeof(SCIP_Longint));
   return ((SCIP_Longint)1) << sigidx; /*lint !e703*/
}

/** returns the bit signature of the given constraint data */
static
SCIP_Longint consdataGetSignature(
   SCIP_CONSDATA*        consdata            /**< set partitioning / packing / covering constraint data */
   )
{
   assert(consdata != NULL);

   if( !consdata->validsignature )
   {
      int i;

      consdata->signature = 0;
      for( i = 0; i < consdata->nvars; ++i )
         consdata->signature |= getVarSignature(consdata->vars[i]);
      consdata->validsignature = TRUE;
   }

   return consdata->signature;
}

/** sorts setppc constraint's variables by non-decreasing variable index */
static
void consdataSort(
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   assert(consdata != NULL);

   if( !consdata->sorted )
   {
      if( consdata->nvars <= 1 )
	 consdata->sorted = TRUE;
      else
      {
	 SCIPsortPtr((void**)consdata->vars, SCIPvarComp, consdata->nvars);
	 consdata->sorted = TRUE;
      }
   }
   assert(consdata->sorted);
#ifdef SCIP_DEBUG
   /* check sorting */
   {
      int v;

      for( v = 0; v < consdata->nvars; ++v )
      {
         assert(v == consdata->nvars-1 || SCIPvarCompare(consdata->vars[v], consdata->vars[v+1]) <= 0);
      }
   }
#endif
}

/** changes the type of a setppc constraint */
static
SCIP_RETCODE setSetppcType(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< setppc constraint */
   SCIP_SETPPCTYPE       setppctype          /**< new type of constraint */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( (SCIP_SETPPCTYPE)consdata->setppctype == setppctype )
      return SCIP_OKAY;

   SCIPdebugMessage(" -> converting <%s> into setppc type %d\n", SCIPconsGetName(cons), setppctype);

   /* remove rounding locks */
   if( SCIPconsIsLocked(cons) )
   {
      int v;

      for( v = 0; v < consdata->nvars; ++v )
      {
         SCIP_CALL( unlockRounding(scip, cons, consdata->vars[v]) );
      }
   }

   /* change the constraint type */
   consdata->setppctype = setppctype; /*lint !e641*/

   /* reinstall rounding locks again */
   if( SCIPconsIsLocked(cons) )
   {
      int v;

      for( v = 0; v < consdata->nvars; ++v )
      {
         SCIP_CALL( lockRounding(scip, cons, consdata->vars[v]) );
      }
   }

   /* remember that we changed a constraint type for clique lifting procedure */
   if( setppctype != SCIP_SETPPCTYPE_COVERING )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      conshdlr = SCIPconsGetHdlr(cons);
      assert(conshdlr != NULL);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      conshdlrdata->updatedsetppctype = TRUE;
   }

   return SCIP_OKAY;
}

/** catches events for variable at given position */
static
SCIP_RETCODE catchEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< set partitioning / packing / covering constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< array position of variable to catch bound change events for */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;
   SCIP_VAR* var;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(consdata->vars != NULL);

   var = consdata->vars[pos];
   assert(var != NULL);

   /* we are catching the following events:
    *
    * - SCIP_EVENTTYPE_BOUNDCHANGED: Is used to count the number of variable fixed locally to zero and one. That helps
    *                                to speed up the propagation
    *
    * - SCIP_EVENTTYPE_VARDELETED: Is caught to remove a deleted variable from the constraint
    *
    * - SCIP_EVENTTYPE_VARFIXED: Is used to get informed if a variable of the constraint was aggregated which means was
    *                            detected to be equal or a negated variable of on other variable. in case of a negation
    *                            this could lead to a redundant constraint if the (other) active variable is also part
    *                            of the constraint.
    */
   eventtype =  SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARDELETED | SCIP_EVENTTYPE_VARFIXED;

   /* catch bound change events on variable */
   SCIP_CALL( SCIPcatchVarEvent(scip, var, eventtype, eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );

   /* update the fixed variables counters for this variable */
   if( SCIPisEQ(scip, SCIPvarGetUbLocal(var), 0.0) )
      consdata->nfixedzeros++;
   else if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), 1.0) )
      consdata->nfixedones++;

   return SCIP_OKAY;
}

/** drops events for variable at given position */
static
SCIP_RETCODE dropEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< set partitioning / packing / covering constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< array position of variable to catch bound change events for */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;
   SCIP_VAR* var;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(consdata->vars != NULL);

   var = consdata->vars[pos];
   assert(var != NULL);

   eventtype =  SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARDELETED | SCIP_EVENTTYPE_VARFIXED;

   /* drop events on variable */
   SCIP_CALL( SCIPdropVarEvent(scip, var, eventtype, eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );

   /* update the fixed variables counters for this variable */
   if( SCIPisEQ(scip, SCIPvarGetUbLocal(var), 0.0) )
      consdata->nfixedzeros--;
   else if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), 1.0) )
      consdata->nfixedones--;

   return SCIP_OKAY;
}

/** catches bound change events for all variables in transformed setppc constraint */
static
SCIP_RETCODE catchAllEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< set partitioning / packing / covering constraint */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* catch event for every single variable */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( catchEvent(scip, cons, eventhdlr, i) );
   }

   return SCIP_OKAY;
}

/** drops bound change events for all variables in transformed setppc constraint */
static
SCIP_RETCODE dropAllEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< set partitioning / packing / covering constraint */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* drop event of every single variable */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( dropEvent(scip, cons, eventhdlr, i) );
   }

   return SCIP_OKAY;
}

/** adds coefficient in setppc constraint */
static
SCIP_RETCODE addCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< setppc constraint */
   SCIP_VAR*             var                 /**< variable to add to the constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;

   assert(var != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if( transformed )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var) );
   }
   assert(var != NULL);
   assert(transformed == SCIPvarIsTransformed(var));

   SCIP_CALL( consdataEnsureVarsSize(scip, consdata, consdata->nvars+1) );
   consdata->vars[consdata->nvars] = var;
   consdata->nvars++;
   if( consdata->validsignature )
      consdata->signature |= getVarSignature(var);
   consdata->sorted = (consdata->nvars == 1);
   consdata->changed = TRUE;

   /* capture the variable */
   SCIP_CALL( SCIPcaptureVar(scip, var) );

   /* if we are in transformed problem, catch the variable's events */
   if( transformed )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      assert(conshdlr != NULL);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* catch bound change events of variable */
      SCIP_CALL( catchEvent(scip, cons, conshdlrdata->eventhdlr, consdata->nvars-1) );

#ifdef VARUSES
      /* if the constraint is currently active, increase the variable usage counter */
      if( SCIPconsIsActive(cons) )
      {
         SCIP_CALL( conshdlrdataIncVaruses(scip, conshdlrdata, var) );
      }
#endif
   }

   /* install the rounding locks for the new variable */
   SCIP_CALL( lockRounding(scip, cons, var) );

   /* add the new coefficient to the LP row */
   if( consdata->row != NULL )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, var, 1.0) );
   }

   consdata->merged = FALSE;

   return SCIP_OKAY;
}

/** deletes coefficient at given position from setppc constraint data */
static
SCIP_RETCODE delCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< set partitioning / packing / covering constraint */
   int                   pos                 /**< position of coefficient to delete */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   var = consdata->vars[pos];
   assert(var != NULL);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(var));

   /* remove the rounding locks for the deleted variable */
   SCIP_CALL( unlockRounding(scip, cons, var) );

   /* if we are in transformed problem, delete the event data of the variable */
   if( SCIPconsIsTransformed(cons) )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* drop bound change events of variable */
      SCIP_CALL( dropEvent(scip, cons, conshdlrdata->eventhdlr, pos) );
   }

   /* delete coefficient from the LP row */
   if( consdata->row != NULL )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, var, -1.0) );
   }

   /* move the last variable to the free slot */
   if( pos != consdata->nvars - 1 )
   {
      consdata->vars[pos] = consdata->vars[consdata->nvars-1];
      consdata->sorted = FALSE;
   }
   consdata->nvars--;
   consdata->validsignature = FALSE;
   consdata->changed = TRUE;

   /* release variable */
   SCIP_CALL( SCIPreleaseVar(scip, &var) );

   return SCIP_OKAY;
}

/** in case a part (more than one variable) in the setppc constraint is independent of every else, we can perform dual
 *  reductions;
 *
 *  (1) set covering
 *      - fix the variable with the smallest object coefficient to one if the constraint is not modifiable and all
 *        variable are independant
 *      - fix all independant variables with negative object coefficient to one
 *      - fix all remaining independant variables to zero 
 *  (2) set partitioning
 *      - fix the variable with the smallest object coefficient to one if the constraint is not modifiable and all
 *        variables are independant
 *      - fix all remaining independant variables to zero
 *  (3) set packing
 *      - fix the variable with the smallest object coefficient to one if the object coefficient is negative or zero
 *        otherwise fix it to one, but only if the constraint is not modifiable and all variables are independant
 *      - fix all remaining independant variables to zero
 *
 * Note: the following dual reduction for set covering and set packing constraints is already performed by the presolver
 *       "dualfix"
 *       (1) in case of a set covering constraint the following dual reduction can be performed:
 *           - if a variable in a set covering constraint is only locked by that constraint and has negative or zero
 *             objective coefficient than it can be fixed to one
 *       (2) in case of a set packing constraint the following dual reduction can be performed:
 *           - if a variable in a set covering constraint is only locked by that constraint and has position or zero 
 *             objective coefficient than it can be fixed to zero
 */
static
SCIP_RETCODE dualPresolving(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< setppc constraint */
   int*                  nfixedvars,         /**< pointer to count number of fixings */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints  */
   SCIP_RESULT*          result              /**< pointer to store the result SCIP_SUCCESS, if presolving was performed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_SETPPCTYPE setppctype;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Real bestobjval;
   SCIP_Real objval;
   SCIP_Real fixval;
   SCIP_Bool infeasible;
   SCIP_Bool fixed;
   SCIP_Bool negated;
   int nfixables;
   int nlockdowns;
   int nlockups;
   int nvars;
   int idx;
   int v;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(nfixedvars != NULL);
   assert(ndelconss != NULL);
   assert(result != NULL);

   /* constraints for which the check flag is set to FALSE, did not contribute to the lock numbers; therefore, we cannot
    * use the locks to decide for a dual reduction using this constraint; for example after a restart the cuts which are
    * added to the problems have the check flag set to FALSE 
    */
   if( !SCIPconsIsChecked(cons) )
      return SCIP_OKAY;
   
   assert(SCIPconsIsActive(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* all fixed variables should be removed at that point */
   assert(consdata->nfixedones == 0);
   assert(consdata->nfixedzeros == 0);
   
   nvars = consdata->nvars;
   
   /* we don't want to consider small constraints (note that the constraints can be modifiable, so we can't delete this
    * constraint) 
    */
   if( nvars < 2 )
      return SCIP_OKAY;

   setppctype = (SCIP_SETPPCTYPE)consdata->setppctype;
   vars = consdata->vars;
   idx = -1;
   bestobjval = SCIP_INVALID;
   
   /* collect the rounding locks depending on the setppc type */
   switch( setppctype )
   {
   case SCIP_SETPPCTYPE_PARTITIONING:
      nlockdowns = 1;
      nlockups = 1;
      break;
   case SCIP_SETPPCTYPE_PACKING:
      nlockdowns = 0;
      nlockups = 1;
      break;
   case SCIP_SETPPCTYPE_COVERING:
      nlockdowns = 1;
      nlockups = 0;
      break;
   default:
      SCIPerrorMessage("unknown setppc type\n");
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }

   nfixables = 0;

   /* check if we can apply the dual reduction; therefore count the number of variables where the setppc has the only
    * locks on this constraint 
    */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert(var != NULL);

      /* the variable should not be (globally) fixed */
      assert(SCIPvarGetLbGlobal(var) < 0.5 && SCIPvarGetUbGlobal(var) > 0.5);
      
      /* in case an other constraints has also locks on that variable we cannot perform a dual reduction on these
       * variables
       */
      if( SCIPvarGetNLocksDown(var) > nlockdowns || SCIPvarGetNLocksUp(var) > nlockups ) 
         continue;

      ++nfixables;
      negated = FALSE;
      
      /* get the active variable */
      SCIP_CALL( SCIPvarGetProbvarBinary(&var, &negated) );
      assert(SCIPvarIsActive(var));

      if( negated )
         objval = -SCIPvarGetObj(var);
      else
         objval = SCIPvarGetObj(var);
      
      /* check if the current variable has a smaller objective coefficient */
      if( SCIPisLT(scip, objval, bestobjval) )
      {
         idx = v;
         bestobjval = objval;
      }
   }

   if( nfixables < 2 )
      return SCIP_OKAY;

   assert(idx >= 0 && idx < nvars);
   assert(bestobjval < SCIPinfinity(scip));
   
   *result = SCIP_SUCCESS;

   /* in case of set packing and set partitioning we fix the remaining variables to zero; Note that this would be also
    * done by the "dualfix" presolver or in the next presolving round of the setppc constraint handler; for performance
    * reasons we do it directly
    */
   if( setppctype != SCIP_SETPPCTYPE_COVERING ) 
   { 
      /* first part of all variables */
      for( v = 0; v < idx; ++v ) 
      {
         var = vars[v];
         assert(var != NULL);

         /* in case an other constraints has also locks on that variable we cannot perform a dual reduction on these
          * variables
          */
         if( SCIPvarGetNLocksDown(var) > nlockdowns || SCIPvarGetNLocksUp(var) > nlockups ) 
            continue;

         SCIP_CALL( SCIPfixVar(scip, var, 0.0, &infeasible, &fixed) ); 
         assert(!infeasible); 
         assert(fixed);
   
         SCIPdebugMessage(" -> fixed <%s> == 0.0\n", SCIPvarGetName(var));
         ++(*nfixedvars);
      }

      /* second part of all variables */
      for( v = idx + 1; v < nvars; ++v )
      {
         var = vars[v];
         assert(var != NULL);

         /* in case an other constraints has also locks on that variable we cannot perform a dual reduction on these
          * variables
          */
         if( SCIPvarGetNLocksDown(var) > nlockdowns || SCIPvarGetNLocksUp(var) > nlockups ) 
            continue;

         SCIP_CALL( SCIPfixVar(scip, var, 0.0, &infeasible, &fixed) );
         assert(!infeasible);
         assert(fixed);
   
         SCIPdebugMessage(" -> fixed <%s> == 0.0\n", SCIPvarGetName(var));
         ++(*nfixedvars);
      }
   }
   /* if we got a set covering constraint and not all variables are locked from this constraint it might not get
    * redundant (which is case if it is not possible to fix at least one variable to one), we fix all redundant
    * variables to their best bound
    */
   else if( nfixables != nvars )
   {
      SCIP_VAR* activevar;

      /* first part of all variables */
      for( v = 0; v < idx; ++v ) 
      {
         var = vars[v];
         assert(var != NULL);

         /* in case an other constraints has also locks on that variable we cannot perform a dual reduction on these
          * variables
          */
         if( SCIPvarGetNLocksDown(var) > nlockdowns || SCIPvarGetNLocksUp(var) > nlockups ) 
            continue;

         activevar = var;
         negated = FALSE;

         /* get the active variable */
         SCIP_CALL( SCIPvarGetProbvarBinary(&activevar, &negated) );
         assert(SCIPvarIsActive(activevar));
         
         if( negated )
            objval = -SCIPvarGetObj(activevar);
         else
            objval = SCIPvarGetObj(activevar);

         if( objval > 0.0 )
            fixval = 0.0;
         else
            fixval = 1.0;
                 
         SCIP_CALL( SCIPfixVar(scip, var, fixval, &infeasible, &fixed) ); 
         assert(!infeasible); 
         assert(fixed);
   
         SCIPdebugMessage(" -> fixed <%s> == %g\n", SCIPvarGetName(var), fixval);
         ++(*nfixedvars);
      }

      /* second part of all variables */
      for( v = idx + 1; v < nvars; ++v )
      {
         var = vars[v];
         assert(var != NULL);

         /* in case an other constraints has also locks on that variable we cannot perform a dual reduction on these
          * variables
          */
         if( SCIPvarGetNLocksDown(var) > nlockdowns || SCIPvarGetNLocksUp(var) > nlockups ) 
            continue;

         activevar = var;
         negated = FALSE;

         /* get the active variable */
         SCIP_CALL( SCIPvarGetProbvarBinary(&activevar, &negated) );
         assert(SCIPvarIsActive(activevar));
         
         if( negated )
            objval = -SCIPvarGetObj(activevar);
         else
            objval = SCIPvarGetObj(activevar);

         if( objval > 0.0 )
            fixval = 0.0;
         else
            fixval = 1.0;
         
         SCIP_CALL( SCIPfixVar(scip, var, fixval, &infeasible, &fixed) );
         assert(!infeasible);
         assert(fixed);
   
         SCIPdebugMessage(" -> fixed <%s> == %g\n", SCIPvarGetName(var), fixval);
         ++(*nfixedvars);
      }
   }

   /* if all variable have our appreciated number of locks and the constraint is not modifiable, or if we have a set
    * covering constraint and the bestobjval is less than or equal to zero, we can fix the variable with the smallest
    * objective coefficient and the constraint gets redundant
    */
   if( (nfixables == nvars && !SCIPconsIsModifiable(cons)) || (setppctype == SCIP_SETPPCTYPE_COVERING && bestobjval <= 0.0) )
   {
      /* in case of a set packing constraint with position objective values, all variables can be fixed to zero; in all
       * other cases the variable with the smallest objective values is fixed to one 
       */
      if( setppctype == SCIP_SETPPCTYPE_PACKING && bestobjval > 0.0 )
         fixval = 0.0;
      else
         fixval = 1.0;
      
      SCIP_CALL( SCIPfixVar(scip, vars[idx], fixval, &infeasible, &fixed) );
      assert(!infeasible);
      assert(fixed);
      
      SCIPdebugMessage(" -> fixed <%s> == %g\n", SCIPvarGetName(vars[idx]), fixval);
      ++(*nfixedvars);

      /* remove constraint since i*/
      SCIP_CALL( SCIPdelCons(scip, cons) );
      ++(*ndelconss);
   }
   
   return SCIP_OKAY;
}

/** find pairs of negated variables in constraint:
 *  partitioning/packing: all other variables must be zero, constraint is redundant
 *  covering: constraint is redundant
 *
 *  find sets of equal variables in constraint:
 *  partitioning/packing: variable must be zero
 *  covering: multiple entries of variable can be replaced by single entry
 */
static
SCIP_RETCODE mergeMultiples(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int*                  nfixedvars,         /**< pointer to store number of fixed variables */
   int*                  ndelconss,          /**< pointer to store number of deleted constraints */
   int*                  nchgcoefs,          /**< pointer to store number of changed coefficients */
   SCIP_Bool*            cutoff              /**< pointer to store whether a fixing leads to a cutoff */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nfixedvars != NULL);
   assert(ndelconss != NULL);
   assert(nchgcoefs != NULL);
   assert(cutoff != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->merged || SCIPconsIsDeleted(cons) )
      return SCIP_OKAY;

   if( consdata->nvars <= 1 )
   {
      consdata->merged = TRUE;
      return SCIP_OKAY;
   }

   assert(consdata->vars != NULL || consdata->nvars == 0);

   /* sorting array after indices of variables, that's only for faster merging */ 
   SCIPsortPtr((void**)consdata->vars, SCIPvarCompActiveAndNegated, consdata->nvars);
   /* setppc sorting now lost */ 
   consdata->sorted = FALSE;

   /* loop backwards through the items: deletion only affects rear items */
   for( v = consdata->nvars - 1; v > 0; --v )
   {
      SCIP_VAR* var1;
      SCIP_VAR* var2;
      SCIP_Bool negated1;
      SCIP_Bool negated2;
      
      negated1 = FALSE;
      negated2 = FALSE;
      
      var1 = consdata->vars[v];
      assert(SCIPvarIsBinary(var1));
      assert(SCIPvarIsActive(var1) || SCIPvarGetStatus(var1) == SCIP_VARSTATUS_NEGATED || SCIPvarGetStatus(var1) == SCIP_VARSTATUS_FIXED);
      if( SCIPvarGetStatus(var1) == SCIP_VARSTATUS_NEGATED )
      {
         var1 = SCIPvarGetNegatedVar(var1);
         negated1 = TRUE;
      }
      assert(var1 != NULL);
      
      var2 = consdata->vars[v-1];
      assert(SCIPvarIsBinary(var2));
      assert(SCIPvarIsActive(var2) || SCIPvarGetStatus(var2) == SCIP_VARSTATUS_NEGATED || SCIPvarGetStatus(var2) == SCIP_VARSTATUS_FIXED);
      if( SCIPvarGetStatus(var2) == SCIP_VARSTATUS_NEGATED )
      {
         var2 = SCIPvarGetNegatedVar(var2);
         negated2 = TRUE;
      }
      assert(var2 != NULL);
      
      if( var1 == var2 )
      {
         SCIP_Bool infeasible;
         SCIP_Bool fixed;

         /* one variables is active and the other is the same negated variable */
         if( negated1 != negated2  )
         {
            /* all other variable have to be zero if it's a partitioning or packing constraint */
            if( consdata->setppctype != SCIP_SETPPCTYPE_COVERING ) /*lint !e641*/
            {
               int i;

               assert(consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING
                  || consdata->setppctype == SCIP_SETPPCTYPE_PACKING); /*lint !e641*/

               for( i = consdata->nvars - 1; i >= 0; --i )
                  if( i != v && i != (v-1) )
                  {
                     SCIP_CALL( SCIPfixVar(scip, consdata->vars[i], 0.0, &infeasible, &fixed) );
                     if( infeasible )
                     {
                        SCIPdebugMessage("setppc constraint <%s>: infeasible fixing <%s> == 0\n",
                           SCIPconsGetName(cons), SCIPvarGetName(consdata->vars[i]));
                        *cutoff = TRUE;
                        return SCIP_OKAY;
                     }

                     if( fixed )
			++(*nfixedvars);
                  }
            }
            /* all setppc-type constraints are redundant */
            SCIP_CALL( SCIPdelCons(scip, cons) );
            ++(*ndelconss);
            return SCIP_OKAY;
         }
         /* both variables are either active or negated */
         else
         {
            /* this variable can be fixed to zero if it's a partitioning or packing constraint */
            if( consdata->setppctype != SCIP_SETPPCTYPE_COVERING ) /*lint !e641*/
            {
               assert(consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING
                  || consdata->setppctype == SCIP_SETPPCTYPE_PACKING); /*lint !e641*/

               SCIP_CALL( SCIPfixVar(scip, var1, negated1 ? 1.0 : 0.0, &infeasible, &fixed) );
               if( infeasible )
               {
                  SCIPdebugMessage("setppc constraint <%s>: infeasible fixing <%s> == %g\n",
                     SCIPconsGetName(cons), SCIPvarGetName(var1), negated1 ? 1.0 : 0.0);
                  *cutoff = TRUE;
                  return SCIP_OKAY;
               }

               if( fixed )
		  ++(*nfixedvars);
            }
            /* multiple entries of variable can be replaced by single entry */
            else
            {
               SCIP_CALL( delCoefPos(scip, cons, v) ); /* only some changed behind position v-1, so it's okay */
               ++(*nchgcoefs);
            }
         }
         consdata->changed = TRUE;
      }
   }
   consdata->merged = TRUE;

   return SCIP_OKAY;
}

/** deletes all zero-fixed variables and replace aggregated variables */
static
SCIP_RETCODE applyFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< set partitioning / packing / covering constraint */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   v = 0;
   while( v < consdata->nvars )
   {
      SCIP_VAR* var;

      var = consdata->vars[v];
      assert(SCIPvarIsBinary(var));

      if( SCIPvarGetUbGlobal(var) < 0.5 )
      {
         assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(var), 0.0));
         SCIP_CALL( delCoefPos(scip, cons, v) );
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
            /* delete old (aggregated) variable */
            SCIP_CALL( delCoefPos(scip, cons, v) );

            /* add representative instead */
            SCIP_CALL( addCoef(scip, cons, repvar) );
         }
         else
            ++v;
      }
   }

   return SCIP_OKAY;
}

/** analyzes conflicting assignment on given constraint where all of the variables where assigned to zero,
 *  and adds conflict constraint to problem
 */
static
SCIP_RETCODE analyzeConflictZero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< set partitioning / packing / covering constraint that detected the conflict */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   /* conflict analysis can only be applied in solving stage and if it is applicable */
   if( (SCIPgetStage(scip) != SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) || !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING
      || consdata->setppctype == SCIP_SETPPCTYPE_COVERING); /*lint !e641*/

   /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
   SCIP_CALL( SCIPinitConflictAnalysis(scip) );
   for( v = 0; v < consdata->nvars; ++v )
   {
      SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[v]) );
   }

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** analyzes conflicting assignment on given constraint where two of the variables where assigned to one,
 *  and adds conflict constraint to problem
 */
static
SCIP_RETCODE analyzeConflictOne(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< set partitioning / packing / covering constraint that detected the conflict */
   )
{
   SCIP_CONSDATA* consdata;
   int v;
   int n;

   /* conflict analysis can only be applied in solving stage and if it is applicable */
   if( (SCIPgetStage(scip) != SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) || !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING
      || consdata->setppctype == SCIP_SETPPCTYPE_PACKING); /*lint !e641*/

   /* initialize conflict analysis, and add the two variables assigned to one to conflict candidate queue */
   SCIP_CALL( SCIPinitConflictAnalysis(scip) );
   n = 0;
   for( v = 0; v < consdata->nvars && n < 2; ++v )
   {
      if( SCIPvarGetLbLocal(consdata->vars[v]) > 0.5 )
      {
         SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[v]) );
         n++;
      }
   }
   assert(n == 2);

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** checks constraint for violation only looking at the fixed variables, applies further fixings if possible */
static
SCIP_RETCODE processFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< set partitioning / packing / covering constraint to be processed */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   SCIP_Bool*            reduceddom,         /**< pointer to store TRUE, if a domain reduction was found */
   SCIP_Bool*            addcut,             /**< pointer to store whether this constraint must be added as a cut */
   SCIP_Bool*            mustcheck           /**< pointer to store whether this constraint must be checked for feasibility */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(reduceddom != NULL);
   assert(addcut != NULL);
   assert(mustcheck != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);
   assert(0 <= consdata->nfixedzeros && consdata->nfixedzeros <= consdata->nvars);
   assert(0 <= consdata->nfixedones && consdata->nfixedones <= consdata->nvars);

   *addcut = FALSE;
   *mustcheck = TRUE;

   /*SCIPdebugMessage("processing constraint <%s> with respect to fixed variables (%d fixed to 0.0, %d fixed to 1.0)\n",
     SCIPconsGetName(cons), consdata->nfixedzeros, consdata->nfixedones);*/

   if( consdata->nfixedones == 1 )
   {
      /* exactly one variable is fixed to 1:
       * - a set covering constraint is feasible anyway and can be disabled
       * - all other variables in a set partitioning or packing constraint must be zero
       */
      if( consdata->setppctype == SCIP_SETPPCTYPE_COVERING ) /*lint !e641*/
      {
         SCIPdebugMessage(" -> disabling set covering constraint <%s>\n", SCIPconsGetName(cons));
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      }
      else
      {
         if( consdata->nfixedzeros < consdata->nvars - 1 )
         {
            SCIP_VAR** vars;
            SCIP_VAR* var;
#ifndef NDEBUG
            SCIP_Bool fixedonefound;
#endif
            SCIP_Bool fixed;
            SCIP_Bool infeasible;
            SCIP_Bool tightened;
            int nvars;
            int v;

            SCIPdebugMessage(" -> fixing all other variables to zero in set packing/partitioning constraint <%s>\n",
               SCIPconsGetName(cons));

            /* unfixed variables exist: fix them to zero;
             * this could result in additional variables fixed to one due to aggregations; in this case, the
             * constraint is infeasible in local bounds
             */
            vars = consdata->vars;
            nvars = consdata->nvars;
#ifndef NDEBUG
            fixedonefound = FALSE;
#endif
            fixed = FALSE;
            for( v = 0; v < nvars && consdata->nfixedones == 1; ++v )
            {
               var = vars[v];
               assert(SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) || SCIPisFeasEQ(scip, SCIPvarGetUbLocal(var), 1.0));
               if( SCIPvarGetLbLocal(var) < 0.5 )
               {
                  SCIP_CALL( SCIPinferBinvarCons(scip, var, FALSE, cons, 0, &infeasible, &tightened) );
                  assert(!infeasible);
                  fixed = fixed || tightened;
                  SCIPdebugMessage("   -> fixed <%s> to zero (tightened=%u)\n", SCIPvarGetName(var), tightened);
               }
#ifndef NDEBUG
               else
                  fixedonefound = TRUE;
#endif
            }
            /* the fixed to one variable must have been found, and at least one variable must have been fixed */
            assert(consdata->nfixedones >= 2 || (fixedonefound && fixed));

            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            *reduceddom = TRUE;
         }

         /* now all other variables are fixed to zero:
          * the constraint is feasible, and if it's not modifiable, it is redundant
          */
         if( !SCIPconsIsModifiable(cons) && consdata->nfixedones == 1 )
         {
            SCIPdebugMessage(" -> disabling set packing/partitioning constraint <%s>\n", SCIPconsGetName(cons));
            SCIP_CALL( SCIPdelConsLocal(scip, cons) );
         }
      }
      *mustcheck = FALSE;
   }

   if( consdata->nfixedones >= 2 )
   {
      /* at least two variables are fixed to 1:
       * - a set covering constraint is feasible anyway and can be disabled
       * - a set partitioning or packing constraint is infeasible
       */
      if( consdata->setppctype == SCIP_SETPPCTYPE_COVERING ) /*lint !e641*/
      {
         SCIPdebugMessage(" -> disabling set covering constraint <%s>\n", SCIPconsGetName(cons));
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      }
      else
      {
         SCIPdebugMessage(" -> conflict on set packing/partitioning constraint <%s>\n", SCIPconsGetName(cons));

         SCIP_CALL( SCIPresetConsAge(scip, cons) );

         /* use conflict analysis to get a conflict constraint out of the conflicting assignment */
         SCIP_CALL( analyzeConflictOne(scip, cons) );

         *cutoff = TRUE;
      }
      *mustcheck = FALSE;
   }
   else if( consdata->nfixedzeros == consdata->nvars )
   {
      /* all variables are fixed to zero:
       * - a set packing constraint is feasible anyway, and if it's unmodifiable, it can be disabled
       * - a set partitioning or covering constraint is infeasible, and if it's unmodifiable, the node
       *   can be cut off -- otherwise, the constraint must be added as a cut and further pricing must
       *   be performed
       */
      assert(consdata->nfixedones == 0);

      if( consdata->setppctype == SCIP_SETPPCTYPE_PACKING ) /*lint !e641*/
      {
         if( !SCIPconsIsModifiable(cons) )
         {
            SCIPdebugMessage(" -> disabling set packing constraint <%s>\n", SCIPconsGetName(cons));
            SCIP_CALL( SCIPdelConsLocal(scip, cons) );
         }
      }
      else
      {
         SCIPdebugMessage(" -> set covering/partitioning constraint <%s> is infeasible\n", SCIPconsGetName(cons));

         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         if( SCIPconsIsModifiable(cons) )
            *addcut = TRUE;
         else
         {
            /* use conflict analysis to get a conflict constraint out of the conflicting assignment */
            SCIP_CALL( analyzeConflictZero(scip, cons) );

            *cutoff = TRUE;
         }
      }
      *mustcheck = FALSE;
   }
   else if( consdata->nfixedzeros == consdata->nvars - 1 && consdata->nfixedones == 0 )
   {
      /* all variables except one are fixed to zero:
       * - a set packing constraint is feasible anyway, and if it's unmodifiable, it can be disabled
       * - an unmodifiable set partitioning or covering constraint is feasible and can be disabled after the
       *   remaining variable is fixed to one
       * - a modifiable set partitioning or covering constraint must be checked manually
       */
      if( consdata->setppctype == SCIP_SETPPCTYPE_PACKING ) /*lint !e641*/
      {
         if( !SCIPconsIsModifiable(cons) )
         {
            SCIPdebugMessage(" -> disabling set packing constraint <%s>\n", SCIPconsGetName(cons));
            SCIP_CALL( SCIPdelConsLocal(scip, cons) );
         }
         *mustcheck = FALSE;
      }
      else if( !SCIPconsIsModifiable(cons) )
      {
         SCIP_VAR** vars;
         SCIP_VAR* var;
         SCIP_Bool infeasible;
         SCIP_Bool tightened;
         int nvars;
         int v;

         /* search the single variable that can be fixed */
         vars = consdata->vars;
         nvars = consdata->nvars;
         for( v = 0; v < nvars; ++v )
         {
            var = vars[v];
            assert(SCIPisFeasZero(scip, SCIPvarGetLbLocal(var)));
            assert(SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) || SCIPisFeasEQ(scip, SCIPvarGetUbLocal(var), 1.0));
            if( SCIPvarGetUbLocal(var) > 0.5 )
            {
               SCIPdebugMessage(" -> fixing remaining variable <%s> to one in set covering/partitioning constraint <%s>\n",
                  SCIPvarGetName(var), SCIPconsGetName(cons));
               SCIP_CALL( SCIPinferBinvarCons(scip, var, TRUE, cons, 0, &infeasible, &tightened) );
               assert(!infeasible);
               assert(tightened);
               break;
            }
         }
         assert(v < nvars);
         assert(consdata->nfixedzeros == consdata->nvars - 1);
         assert(consdata->nfixedones == 1);

         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
         *reduceddom = TRUE;
         *mustcheck = FALSE;
      }
   }
   assert(consdata->nfixedzeros + consdata->nfixedones <= consdata->nvars);

   return SCIP_OKAY;
}

/** checks constraint for violation, returns TRUE iff constraint is feasible */
static
SCIP_Bool checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< set partitioning / packing / covering constraint to be checked */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   SCIP_VAR** vars;
   SCIP_Real solval;
   SCIP_Real sum;
   SCIP_Real sumbound;
   int nvars;
   int v;

   /* calculate the constraint's activity */
   vars = consdata->vars;
   nvars = consdata->nvars;
   sum = 0.0;
   sumbound = ((SCIP_SETPPCTYPE)consdata->setppctype == SCIP_SETPPCTYPE_COVERING ? 1.0 : 1.0 + 2*SCIPfeastol(scip));
   for( v = 0; v < nvars && sum < sumbound; ++v )  /* if sum >= sumbound, the feasibility is clearly decided */
   {
      assert(SCIPvarIsBinary(vars[v]));
      solval = SCIPgetSolVal(scip, sol, vars[v]);
      assert(SCIPisFeasGE(scip, solval, 0.0) && SCIPisFeasLE(scip, solval, 1.0));
      sum += solval;
   }

   switch( consdata->setppctype )
   {
   case SCIP_SETPPCTYPE_PARTITIONING:
      return SCIPisFeasEQ(scip, sum, 1.0);
   case SCIP_SETPPCTYPE_PACKING:
      return SCIPisFeasLE(scip, sum, 1.0);
   case SCIP_SETPPCTYPE_COVERING:
      return SCIPisFeasGE(scip, sum, 1.0);
   default:
      SCIPerrorMessage("unknown setppc type\n");
      SCIPABORT();
      return FALSE; /*lint !e527*/
   }
}

/** creates an LP row in a set partitioning / packing / covering constraint data object */
static
SCIP_RETCODE createRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< set partitioning / packing / covering constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real lhs;
   SCIP_Real rhs;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);

   switch( consdata->setppctype )
   {
   case SCIP_SETPPCTYPE_PARTITIONING:
      lhs = 1.0;
      rhs = 1.0;
      break;
   case SCIP_SETPPCTYPE_PACKING:
      lhs = -SCIPinfinity(scip);
      rhs = 1.0;
      break;
   case SCIP_SETPPCTYPE_COVERING:
      lhs = 1.0;
      rhs = SCIPinfinity(scip);
      break;
   default:
      SCIPerrorMessage("unknown setppc type\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->row, SCIPconsGetHdlr(cons), SCIPconsGetName(cons), lhs, rhs,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );

   SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, consdata->row, consdata->nvars, consdata->vars, 1.0) );

   return SCIP_OKAY;
}

/** adds setppc constraint as cut to the LP */
static
SCIP_RETCODE addCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< setppc constraint */
   SCIP_SOL*             sol                 /**< primal CIP solution, NULL for current LP solution */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row == NULL )
   {
      /* convert set partitioning constraint data into LP row */
      SCIP_CALL( createRow(scip, cons) );
   }
   assert(consdata->row != NULL);

   /* insert LP row as cut */
   if( !SCIProwIsInLP(consdata->row) )
   {
      SCIPdebugMessage("adding constraint <%s> as cut to the LP\n", SCIPconsGetName(cons));
      SCIP_CALL( SCIPaddCut(scip, sol, consdata->row, FALSE) );
   }

   return SCIP_OKAY;
}

/** checks constraint for violation, and adds it as a cut if possible */
static
SCIP_RETCODE separateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< set partitioning / packing / covering constraint to be separated */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   SCIP_Bool*            separated,          /**< pointer to store TRUE, if a cut was found */
   SCIP_Bool*            reduceddom          /**< pointer to store TRUE, if a domain reduction was found */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool addcut;
   SCIP_Bool mustcheck;

   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(separated != NULL);
   assert(reduceddom != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);
   assert(0 <= consdata->nfixedzeros && consdata->nfixedzeros <= consdata->nvars);
   assert(0 <= consdata->nfixedones && consdata->nfixedones <= consdata->nvars);

   /* skip constraints already in the LP */
   if( sol == NULL && consdata->row != NULL && SCIProwIsInLP(consdata->row) )
      return SCIP_OKAY;

   SCIPdebugMessage("separating constraint <%s>\n", SCIPconsGetName(cons));

   /* check constraint for violation only looking at the fixed variables, apply further fixings if possible */
   if( sol == NULL )
   {
      SCIP_CALL( processFixings(scip, cons, cutoff, reduceddom, &addcut, &mustcheck) );
   }
   else
   {
      mustcheck = TRUE;
      addcut = FALSE;
   }

   if( mustcheck )
   {
      assert(!addcut);

      /* variable's fixings didn't give us any information -> we have to check the constraint */
      if( sol == NULL && consdata->row != NULL )
      {
         SCIP_Real feasibility;

         assert(!SCIProwIsInLP(consdata->row));
         feasibility = SCIPgetRowLPFeasibility(scip, consdata->row);
         addcut = SCIPisFeasNegative(scip, feasibility);
      }
      else
         addcut = !checkCons(scip, consdata, sol);

      if( !addcut )
      {
         /* constraint was feasible -> increase age */
         SCIP_CALL( SCIPincConsAge(scip, cons) );
      }
   }

   if( addcut )
   {
      /* insert LP row as cut */
      SCIP_CALL( addCut(scip, cons, sol) );
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      *separated = TRUE;
   }

   return SCIP_OKAY;
}

/** enforces the pseudo solution on the given constraint */
static
SCIP_RETCODE enforcePseudo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< set partitioning / packing / covering constraint to be separated */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   SCIP_Bool*            infeasible,         /**< pointer to store TRUE, if the constraint was infeasible */
   SCIP_Bool*            reduceddom,         /**< pointer to store TRUE, if a domain reduction was found */
   SCIP_Bool*            solvelp             /**< pointer to store TRUE, if the LP has to be solved */
   )
{
   SCIP_Bool addcut;
   SCIP_Bool mustcheck;

   assert(!SCIPhasCurrentNodeLP(scip));
   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(infeasible != NULL);
   assert(reduceddom != NULL);
   assert(solvelp != NULL);

   /* check constraint for violation only looking at the fixed variables, apply further fixings if possible */
   SCIP_CALL( processFixings(scip, cons, cutoff, reduceddom, &addcut, &mustcheck) );

   if( mustcheck )
   {
      SCIP_CONSDATA* consdata;

      assert(!addcut);

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      if( checkCons(scip, consdata, NULL) )
      {
         /* constraint was feasible -> increase age */
         SCIP_CALL( SCIPincConsAge(scip, cons) );
      }
      else
      {
         /* constraint was infeasible -> reset age */
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         *infeasible = TRUE;
      }
   }

   if( addcut )
   {
      /* a cut must be added to the LP -> we have to solve the LP immediately */
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      *solvelp = TRUE;
   }

   return SCIP_OKAY;
}

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashGetKeySetppccons)
{  /*lint --e{715}*/
   /* the key is the element itself */
   return elem;
}

/** returns TRUE iff both keys are equal; two constraints are equal if they have the same variables */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqSetppccons)
{
#ifndef NDEBUG
   SCIP* scip;
#endif
   SCIP_CONSDATA* consdata1;
   SCIP_CONSDATA* consdata2;
   SCIP_Bool coefsequal;
   int i;

   consdata1 = SCIPconsGetData((SCIP_CONS*)key1);
   consdata2 = SCIPconsGetData((SCIP_CONS*)key2);
   assert(consdata1->sorted);
   assert(consdata2->sorted);
#ifndef NDEBUG
   scip = (SCIP*)userptr;
   assert(scip != NULL);
#endif

   /* checks trivial case */
   if( consdata1->nvars != consdata2->nvars )
      return FALSE;

   coefsequal = TRUE;

   for( i = 0; i < consdata1->nvars; ++i )
   {
      /* tests if variables are equal */
      if( consdata1->vars[i] != consdata2->vars[i] )
      {
         assert(SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == 1 ||
            SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == -1);
         coefsequal = FALSE;
         break;
      }
      assert(SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == 0);
   }

   return coefsequal;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValSetppccons)
{
   SCIP_CONSDATA* consdata;
   unsigned int hashval;
   int minidx;
   int mididx;
   int maxidx;
#ifndef NDEBUG
   SCIP* scip;

   scip = (SCIP*)userptr;
   assert(scip != NULL);
#endif

   consdata = SCIPconsGetData((SCIP_CONS*)key);
   assert(consdata != NULL);
   assert(consdata->nvars > 0);

   /* sorts the constraints */
   consdataSort(consdata);

   minidx = SCIPvarGetIndex(consdata->vars[0]);
   mididx = SCIPvarGetIndex(consdata->vars[consdata->nvars / 2]);
   maxidx = SCIPvarGetIndex(consdata->vars[consdata->nvars - 1]);
   assert(minidx >= 0 && minidx <= maxidx);

   hashval = (consdata->nvars << 29) + (minidx << 22) + (mididx << 11) + maxidx; /*lint !e701*/

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


/** add extra clique-constraints resulting from a given cliquepartition to SCIP */
static
SCIP_RETCODE addExtraCliques(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_VAR**const       binvars,            /**< binary variables to create clique constraints */
   int const             nbinvars,           /**< number of binary variables to create clique constraints */
   int*const             cliquepartition,    /**< clique partition of binary variables */
   int const             ncliques,           /**< number of cliques in cliquepartition */
   SCIP_CONS**const      usefulconss,        /**< storage for created constraints */
   int*const             nusefulconss,       /**< pointer to store number of useful created constraints */
   SCIP_VAR**const       usefulvars,         /**< storage for all found variables */
   int*const             nusefulvars,        /**< pointer to store number of added variables */
   SCIP_HASHMAP*const    vartoindex,         /**< hashmap mapping variables to indices */
   int*const             varnconss,          /**< storage for remembering the number of constraints a variable occurs */
   int*const             maxnvarconsidx,     /**< storage for the maximal number of occurances of a variable */
   int**const            varconsidxs,        /**< storage for constraint indices in which the corresponding variable exists */
   int*const             maxnvars,           /**< pointer to store maximal number of variables of a constraint */
   int const             nrounds,            /**< actual presolving round */
   int*const             nfixedvars,         /**< pointer to count number of deleted variables */
   int*const             naddconss,          /**< pointer to count number of added constraints */
   int*const             ndelconss,          /**< pointer to count number of deleted constraints */
   int*const             nchgcoefs,          /**< pointer to count number of deleted coefficients */
   SCIP_Bool*const       cutoff              /**< pointer to store if the problem is infeasible due to a fixing */
   )
{
   SCIP_CONS* cliquecons;
   char name[SCIP_MAXSTRLEN];
   int varindex;
   int lastclqidx;
   int nadded;
   int c;
   int v;

   assert(scip != NULL);
   assert(binvars != NULL || nbinvars == 0);
   assert(cliquepartition != NULL || nbinvars == 0);
   assert(ncliques >= 0 && ncliques <= nbinvars);
   assert(usefulconss != NULL);
   assert(nusefulconss != NULL);
   assert(usefulvars != NULL);
   assert(nusefulvars != NULL);
   assert(vartoindex != NULL);
   assert(varnconss != NULL);
   assert(maxnvarconsidx != NULL);
   assert(varconsidxs != NULL);
   assert(maxnvars != NULL);
   assert(nfixedvars != NULL);
   assert(naddconss != NULL);
   assert(ndelconss != NULL);
   assert(nchgcoefs != NULL);
   assert(cutoff != NULL);

   /* no given binary variables */
   if( nbinvars == 0 || ncliques == 0 )
      return SCIP_OKAY;

   assert(binvars != NULL);
   assert(cliquepartition != NULL);

   /* no useful clique information */
   if( ncliques == nbinvars )
      return SCIP_OKAY;

   lastclqidx = 0;

   /* @todo: maybe sort cliques and accordingly the variables so it will be faster to add the constraints */
   for( c = 0; c < ncliques - 1; ++c )
   {
      if( lastclqidx >= cliquepartition[c] )
	 continue;

      nadded = 0;

      /* name the clique constraint */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "extra_clq_%d_round_%d", cliquepartition[c], nrounds);
      SCIP_CALL( SCIPcreateConsSetpack(scip, &cliquecons, name, 0, NULL,
	    TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      /* add variables to clique constraint */
      for( v = c; v < ncliques - 1; ++v )
      {
	 if( cliquepartition[c] == cliquepartition[v] )
	 {
	    SCIP_CALL( addCoef(scip, cliquecons, binvars[v]) );
	    ++nadded;
	 }
      }

      /* @todo: try to find a good value for what are enough variables to create this constraint, maybe at least
       *        (nmaxvars(over all conss)-nminvars(over all conss))/2 */
      if( nadded >= 2 )
      {
	 SCIP_CONSDATA* cliqueconsdata;

	 SCIPdebugMessage(" -> adding clique constraint: ");
	 SCIPdebugPrintCons(scip, cliquecons, NULL);
	 SCIP_CALL( SCIPaddCons(scip, cliquecons) );
	 ++(*naddconss);

	 /* we only want to consider merged constraints */
	 SCIP_CALL( mergeMultiples(scip, cliquecons, nfixedvars, ndelconss, nchgcoefs, cutoff) );
	 if( *cutoff )
	 {
	    SCIP_CALL( SCIPreleaseCons(scip, &cliquecons) );

	    return SCIP_OKAY;
	 }

	 cliqueconsdata = SCIPconsGetData(cliquecons);
	 assert(cliqueconsdata != NULL);

	 /* the artificial constraints could be deleled while merging */
	 if( !SCIPconsIsDeleted(cliquecons) && nadded - cliqueconsdata->nfixedzeros >= 2 )
	 {
	    assert(cliqueconsdata->nfixedones == 0);

	    /* save maximal number of vars */
	    if( cliqueconsdata->nvars > *maxnvars )
	       *maxnvars = cliqueconsdata->nvars;

	    /* adding variables and information about occurances to local data structure */
	    for( v = cliqueconsdata->nvars - 1; v >= 0; --v )
	    {
	       SCIP_VAR* var;

	       var = cliqueconsdata->vars[v];

	       /* there should be no fixed variable to one in our new clique constraint */
	       assert(SCIPvarGetLbLocal(var) < 0.5);

	       if( SCIPvarGetUbLocal(var) < 0.5 )
		  continue;

	       if( !SCIPhashmapExists(vartoindex, (void*) var) )
	       {
		  SCIP_VAR* tmpvar;

		  usefulvars[*nusefulvars] = var;
		  ++(*nusefulvars);
		  varindex = *nusefulvars;
		  SCIP_CALL( SCIPhashmapInsert(vartoindex, (void*) var, (void*) (size_t) varindex) );

		  /* get the maximal number of occurances of this variable, if this variables  */
		  tmpvar = SCIPvarIsNegated(var) ? SCIPvarGetNegatedVar(var) : var;
		  maxnvarconsidx[varindex] = SCIPvarGetNLocksDown(tmpvar) + SCIPvarGetNLocksUp(tmpvar);
		  SCIP_CALL( SCIPallocBufferArray(scip, &(varconsidxs[varindex]), maxnvarconsidx[varindex]) ); /*lint !e866*/
	       }
	       else
	       {
		  assert(SCIPhashmapGetImage(vartoindex, (void*) var) != NULL);
		  varindex = (int) (size_t) SCIPhashmapGetImage(vartoindex, (void*) var);

		  /* the number of occurances of a variable is not limited by the locks (so maybe we have to increase
		   * memory), because for examples converted cuts are not check and therefore they have no locks on
		   * their variables */
		  if( varnconss[varindex] == maxnvarconsidx[varindex] )
		  {
		     maxnvarconsidx[varindex] = SCIPcalcMemGrowSize(scip, maxnvarconsidx[varindex] + 1);
		     SCIP_CALL( SCIPreallocBufferArray(scip, &(varconsidxs[varindex]), maxnvarconsidx[varindex]) ); /*lint !e866*/
		  }

		  assert(varnconss[varindex] < maxnvarconsidx[varindex]);
		  /* add the constraint number to the variable list */
		  varconsidxs[varindex][varnconss[varindex]] = *nusefulconss;
		  /* increase number of occurances for variables */
		  ++(varnconss[varindex]);
	       }
	    }

	    /* save the type and constraint */
	    usefulconss[*nusefulconss] = cliquecons;
	    ++(*nusefulconss);
	 }
	 SCIP_CALL( SCIPreleaseCons(scip, &cliquecons) );
      }
      else
      {
	 SCIP_CALL( SCIPreleaseCons(scip, &cliquecons) );
      }
      lastclqidx = cliquepartition[c];
   }

   return SCIP_OKAY;
}


/** start to collect setpartitioning and setpacking constraints, and try to remove fixed variables and merged these
 *  constraints
 */
static
SCIP_RETCODE collectCliqueConss(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS**const      conss,              /**< constraint set */
   int const             nconss,             /**< number of constraints in constraint set */
   SCIP_CONS**const      usefulconss,        /**< storage for created constraints */
   int*const             nusefulconss,       /**< pointer to store number of useful created constraints */
   int*const             nfixedvars,         /**< pointer to count number of deleted variables */
   int*const             ndelconss,          /**< pointer to count number of deleted constraints */
   int*const             nchgcoefs,          /**< pointer to count number of deleted coefficients */
   SCIP_Bool*const       cutoff              /**< pointer to store if the problem is infeasible due to a fixing */
   )
{
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_Bool reduceddom;
   SCIP_Bool addcut;
   SCIP_Bool mustcheck;
   int c;

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   assert(usefulconss != NULL);
   assert(nusefulconss != NULL);
   assert(nfixedvars != NULL);
   assert(ndelconss != NULL);
   assert(nchgcoefs != NULL);
   assert(cutoff != NULL);

   if( nconss == 0 )
      return SCIP_OKAY;

   assert(conss != NULL);

   for( c = nconss - 1; c >= 0; --c )
   {
      cons = conss[c];

      /* we only want to consider constraints with either active or negated of active variables, applyfixings removes
       * aggregated and fixed variables to zero, processFixings removes fixings to one but no aggregation
       *
       * @todo: maybe write a new method for deleting aggregations and all fixings
       */
      SCIP_CALL( applyFixings(scip, cons) );
      SCIP_CALL( processFixings(scip, cons, cutoff, &reduceddom, &addcut, &mustcheck) );
      if( *cutoff )
         return SCIP_OKAY;

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* we only want to consider merged constraints */
      SCIP_CALL( mergeMultiples(scip, cons, nfixedvars, ndelconss, nchgcoefs, cutoff) );
      if( *cutoff )
         return SCIP_OKAY;

      if( SCIPconsIsModifiable(cons) || !SCIPconsIsActive(cons) )
         continue;

      assert(consdata->nfixedones == 0);

      if( consdata->nvars == 0 )
         continue;

      /* @todo: check for covering constraints with only two variables which are equal to a packing constraint with
       * negated variables */
      if( consdata->setppctype != SCIP_SETPPCTYPE_COVERING ) /*lint !e641*/
      {
         assert(consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING || consdata->setppctype == SCIP_SETPPCTYPE_PACKING); /*lint !e641*/

         usefulconss[*nusefulconss] = cons;
         ++(*nusefulconss);
      }
   }

   return SCIP_OKAY;
}

/** creating all necessary data in array structure, collect all clique constraint variables and occurances,
 *  @note works only with merged and active not set-covering constraints
 */
static
SCIP_RETCODE collectCliqueData(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS**const      usefulconss,        /**< clique constraints */
   int const             nusefulconss,       /**< number of clique constraints */
   SCIP_VAR**const       usefulvars,         /**< storage for all found variables */
   int*const             nusefulvars,        /**< pointer to store number of added variables */
   SCIP_HASHMAP*const    vartoindex,         /**< hashmap mapping variables to indices */
   int*const             varnconss,          /**< storage for remembering the number of constraints a variable occurs */
   int*const             maxnvarconsidx,     /**< storage for the maximal number of occurances of a variable */
   int**const            varconsidxs,        /**< storage for constraint indices in which the corresponding variable exists */
   int*const             maxnvars            /**< pointer to store maximal number of variables of a constraint */
   )
{
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   int varindex;
   int c;
   int v;

   assert(scip != NULL);
   assert(usefulconss != NULL || nusefulconss == 0);
   assert(usefulvars != NULL);
   assert(nusefulvars != NULL);
   assert(vartoindex != NULL);
   assert(varnconss != NULL);
   assert(maxnvarconsidx != NULL);
   assert(varconsidxs != NULL);
   assert(maxnvars != NULL);

   if( nusefulconss == 0 )
      return SCIP_OKAY;

   assert(usefulconss != NULL);

   for( c = nusefulconss - 1; c >= 0; --c )
   {
      cons = usefulconss[c];

      assert(SCIPconsIsActive(cons));

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* here we should have no covering constraints anymore and the constraint data should be merged */
      assert(consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING || consdata->setppctype == SCIP_SETPPCTYPE_PACKING); /*lint !e641*/
      assert(consdata->merged);

      /* save maximal number of vars */
      if( consdata->nvars > *maxnvars )
         *maxnvars = consdata->nvars;

      /* adding variables and information about occurances to local data structure */
      for( v = consdata->nvars - 1; v >= 0; --v )
      {
         SCIP_VAR* var;

         var = consdata->vars[v];

         /* don't remember fixed vars */
         if( SCIPvarGetLbLocal(var) > 0.5 && SCIPvarGetUbLocal(var) < 0.5 )
            continue;

	 /* only collect active or negated active varibels */
	 assert(SCIPvarIsActive(var) || (SCIPvarIsNegated(var) && SCIPvarIsActive(SCIPvarGetNegationVar(var))));

         if( !SCIPhashmapExists(vartoindex, (void*) var) )
         {
            SCIP_VAR* tmpvar;

            usefulvars[*nusefulvars] = var;
            ++(*nusefulvars);
            varindex = *nusefulvars;
            SCIP_CALL( SCIPhashmapInsert(vartoindex, (void*) var, (void*) (size_t) varindex) );

            /* get the maximal number of occurances of this variable, if this variables  */
            tmpvar = SCIPvarIsNegated(var) ? SCIPvarGetNegatedVar(var) : var;
            maxnvarconsidx[varindex] = SCIPvarGetNLocksDown(tmpvar) + SCIPvarGetNLocksUp(tmpvar);
            SCIP_CALL( SCIPallocBufferArray(scip, &(varconsidxs[varindex]), maxnvarconsidx[varindex]) );  /*lint !e866*/
         }
         else
         {
            assert(SCIPhashmapGetImage(vartoindex, (void*) var) != NULL);
            varindex = (int) (size_t) SCIPhashmapGetImage(vartoindex, (void*) var);
         }

         /* the number of occurances of a variable is not limited by the locks (so maybe we have to increase memory),
          * because for examples converted cuts are not check and therefore they have no locks on their variables */
         if( varnconss[varindex] == maxnvarconsidx[varindex] )
         {
            maxnvarconsidx[varindex] = SCIPcalcMemGrowSize(scip, maxnvarconsidx[varindex] + 1);
            SCIP_CALL( SCIPreallocBufferArray(scip, &(varconsidxs[varindex]), maxnvarconsidx[varindex]) ); /*lint !e866*/
         }

         assert(varnconss[varindex] < maxnvarconsidx[varindex]);
         /* add the constraint number to the variable list */
         varconsidxs[varindex][varnconss[varindex]] = c;
         /* increase number of occurances for variables */
         ++(varnconss[varindex]);
      }
   } /* data structure created */

   return SCIP_OKAY;
}

/** correct clique data due to an aggregation */
static
void deleteCliqueDataEntry(
   SCIP_VAR*const        var,                /**< variable which appears less */
   int const             considx,            /**< constraint index which to remove */
   SCIP_HASHMAP*const    vartoindex,         /**< hashmap mapping variables to indices */
   int*const             varnconss,          /**< storage for remembering the number of constraints a variable occurs */
   int**const            varconsidxs         /**< storage for constraint indices in which the corresponding variable exists */
   )
{
   int varindex;
   int i;
#ifndef NDEBUG
   SCIP_Bool found = FALSE;
#endif

   assert(var != NULL);
   assert(SCIPvarGetLbLocal(var) < 0.5 && SCIPvarGetUbLocal(var) > 0.5);
   assert(considx >= 0);
   assert(vartoindex != NULL);
   assert(varnconss != NULL);
   assert(varconsidxs != NULL);

   assert(SCIPhashmapGetImage(vartoindex, (void*) var) != NULL);
   varindex = (int) (size_t) SCIPhashmapGetImage(vartoindex, (void*) var);

   /* remove entry of variable at the given position */
   for( i = 0; i < varnconss[varindex]; ++i )
   {
      if( varconsidxs[varindex][i] == considx )
      {
	 varconsidxs[varindex][i] = varconsidxs[varindex][varnconss[varindex] - 1];
#ifndef NDEBUG
	 found = TRUE;
#endif
	 --(varnconss[varindex]);
	 break;
      }
   }
   assert(found);
}

/* correct local data structure, add constraint entry to variable data  */
static
SCIP_RETCODE addCliqueDataEntry(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_VAR*const        addvar,             /**< variable which was added */
   int const             considx,            /**< constraint index which to add */
   SCIP_Bool const       maybenew,           /**< could be a new variables, a negated of an already existing */
   SCIP_VAR**const       usefulvars,         /**< storage for all found variables */
   int*const             nusefulvars,        /**< pointer to store number of added variables */
   SCIP_HASHMAP*const    vartoindex,         /**< hashmap mapping variables to indices */
   int*const             varnconss,          /**< storage for remembering the number of constraints a variable occurs */
   int*const             maxnvarconsidx,     /**< storage for the maximal number of occurances of a variable */
   int**const            varconsidxs         /**< storage for constraint indices in which the corresponding variable exists */
   )
{
   int varindex;

   assert(scip != NULL);
   assert(addvar != NULL);
   assert(SCIPvarGetLbLocal(addvar) < 0.5 && SCIPvarGetUbLocal(addvar) > 0.5);
   assert(usefulvars != NULL);
   assert(nusefulvars != NULL);
   assert(vartoindex != NULL);
   assert(varnconss != NULL);
   assert(maxnvarconsidx != NULL);
   assert(varconsidxs != NULL);

   /* we add the variable to the hashmap if its new */
   if( maybenew && !SCIPhashmapExists(vartoindex, (void*) addvar) )
   {
      assert(SCIPvarIsActive(addvar) || SCIPvarIsNegated(addvar));
      assert(SCIPvarGetNegatedVar(addvar) != NULL && SCIPhashmapExists(vartoindex, (void*) SCIPvarGetNegatedVar(addvar)));

      /* @note because we can only have created a negated variable, and we already alloacted enough memory for
       * all (even not existing) negated variables the usefulvars array should be big enough
       */
      SCIPsortedvecInsertDownPtr((void**)usefulvars, SCIPvarCompActiveAndNegated, addvar, nusefulvars, NULL);
      varindex = *nusefulvars;
      SCIP_CALL( SCIPhashmapInsert(vartoindex, (void*) addvar, (void*) (size_t) varindex) );

      assert(varconsidxs[varindex] == NULL);

      maxnvarconsidx[varindex] = 1;
      SCIP_CALL( SCIPallocBufferArray(scip, &(varconsidxs[varindex]), maxnvarconsidx[varindex]) ); /*lint !e866*/
      varnconss[varindex] = 0;
   }
   else
   {
      varindex = (int) (size_t) SCIPhashmapGetImage(vartoindex, (void*) addvar);

      /* grow the needed memory if we added a variable */
      if( varnconss[varindex] == maxnvarconsidx[varindex] )
      {
	 maxnvarconsidx[varindex] = SCIPcalcMemGrowSize(scip, maxnvarconsidx[varindex] + 1);
	 SCIP_CALL( SCIPreallocBufferArray(scip, &(varconsidxs[varindex]), maxnvarconsidx[varindex]) ); /*lint !e866*/
      }
   }
   assert(varnconss[varindex] < maxnvarconsidx[varindex]);
   varconsidxs[varindex][varnconss[varindex]] = considx;

   /* increase number of occurances for variables */
   ++(varnconss[varindex]);

   return SCIP_OKAY;
}


/** check if constraint is already redundant or infeasible due to fixings, fix or aggregate left over variables if
 *  possible
 */
static
SCIP_RETCODE presolvePropagateCons(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< constraint */
   SCIP_Bool const       aggregate,          /**< try to aggregate if possible */
   SCIP_VAR**            undoneaggrvars,     /**< array to store aggregation variables, if aggregation is not performed
					      *   yet; both variables are standing next to each other; or NULL if
					      *   aggregate == TRUE
					      */
   SCIP_Bool*            undoneaggrtypes,    /**< array to store aggregation type, if aggregation is not performed yet;
					      *   type FALSE means the aggregation is of the form x + y = 1; type TRUE means
					      *   the aggregation is of the form x = y; or NULL if aggregate == TRUE
					      */
   int*const             naggregations,      /**< pointer to store number of aggregations which are not yet performed;
					      *   or NULL if aggregate == TRUE
					      */
   int*const             saggregations,      /**< pointer to store size of the array for aggregation type and two times
					      *   the value is the size of the array for the aggregation variables which
					      *   are not yet performed; or NULL if aggregate == TRUE
					      */
   int*const             nfixedvars,         /**< pointer to count number of deleted variables */
   int*const             naggrvars,          /**< pointer to count number of aggregated variables */
   int*const             ndelconss,          /**< pointer to count number of deleted constraints */
   SCIP_Bool*const       cutoff              /**< pointer to store if the problem is infeasible due to a fixing */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int nvars;
   int v;
   SCIP_Bool fixed;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(ndelconss != NULL);
   assert(cutoff != NULL);

   if( !SCIPconsIsActive(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->propagated )
      return SCIP_OKAY;

   consdata->propagated = TRUE;

   vars = consdata->vars;
   nvars = consdata->nvars;

   /* no variables left */
   if( nvars == 0 && !SCIPconsIsModifiable(cons) )
   {
      if( consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING || consdata->setppctype == SCIP_SETPPCTYPE_COVERING ) /*lint !e641*/
      {
	 SCIPdebugMessage("empty set-partition/-covering constraint found -> cutoff\n");
	 *cutoff = TRUE;

	 return SCIP_OKAY;
      }
      else
      {
	 assert(consdata->setppctype == SCIP_SETPPCTYPE_PACKING); /*lint !e641*/

	 /* delete constraint */
	 SCIPdebugMessage(" -> deleting constraint <%s>, no variables left\n", SCIPconsGetName(cons));
	 SCIP_CALL( SCIPdelCons(scip, cons) );
	 ++(*ndelconss);

	 return SCIP_OKAY;
      }
   }

   /* more then two variables are fixed */
   if( consdata->nfixedones > 1 )
   {
      /* at least two variables are fixed to 1:
       * - a set covering constraint is feasible anyway and can be deleted
       * - a set partitioning or packing constraint is infeasible
       */
      if( consdata->setppctype == SCIP_SETPPCTYPE_COVERING ) /*lint !e641*/
      {
	 /* delete constraint */
	 SCIPdebugMessage(" -> deleting set-covering constraint <%s>, at least two variables are fixed to 1\n", SCIPconsGetName(cons));
	 SCIP_CALL( SCIPdelCons(scip, cons) );
	 ++(*ndelconss);

	 return SCIP_OKAY;
      }

      SCIPdebugMessage("set partitioning / packing constraint <%s> is infeasible\n", SCIPconsGetName(cons));
      *cutoff = TRUE;

      return SCIP_OKAY;
   }

   if( consdata->nfixedones == 1 )
   {
      /* exactly one variable is fixed to 1:
       * - a set covering constraint is feasible anyway and can be disabled
       * - all other variables in a set partitioning or packing constraint must be zero
       */
      if( consdata->setppctype != SCIP_SETPPCTYPE_COVERING && consdata->nfixedzeros < nvars - 1 ) /*lint !e641*/
      {
         assert(vars != NULL);

	 for( v = nvars - 1; v >= 0; --v )
	 {
	    if( SCIPvarGetLbLocal(vars[v]) + 0.5 < SCIPvarGetUbLocal(vars[v]) )
	    {
	       SCIPdebugMessage("trying to fix <%s> to 0 due to at least one variable is already fixed to 1\n", SCIPvarGetName(vars[v]));

	       /* fix all remaining variables to zero, constraint is already feasible or infeasible */
	       SCIP_CALL( SCIPfixVar(scip, vars[v], 0.0, cutoff, &fixed) );
	       if( *cutoff )
	       {
		  SCIPdebugMessage("setppc constraint <%s>: infeasible fixing <%s> == 0\n",
		     SCIPconsGetName(cons), SCIPvarGetName(vars[v]));

		  return SCIP_OKAY;
	       }

	       assert(fixed);
	       ++(*nfixedvars);
	    }
	 }
      }

      if( !SCIPconsIsModifiable(cons) || consdata->setppctype == SCIP_SETPPCTYPE_COVERING ) /*lint !e641*/
      {
	 /* delete constraint */
	 SCIPdebugMessage(" -> deleting constraint <%s>, all variable are fixed\n", SCIPconsGetName(cons));
	 assert(SCIPconsIsActive(cons));
	 SCIP_CALL( SCIPdelCons(scip, cons) );
	 ++(*ndelconss);
      }

      return SCIP_OKAY;
   }

   /* other propagations can only be done on not modifiable constraints */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   assert(vars != NULL);

   /* all variables were fixed to zero then either delete the constraint or stop with infeasibility */
   if( consdata->nfixedzeros == nvars )
   {
      assert(consdata->nfixedones == 0);

      /* all variables are fixed to zero:
       * - a set packing constraint is feasible anyway and can be deleted
       * - a set partitioning or covering constraint is infeasible, and so is the whole problem
       */
      if( consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING || consdata->setppctype == SCIP_SETPPCTYPE_COVERING ) /*lint !e641*/
      {
	 SCIPdebugMessage("set partitioning / covering constraint <%s> is infeasible\n", SCIPconsGetName(cons));
	 *cutoff = TRUE;

	 return SCIP_OKAY;
      }

      /* delete constraint */
      SCIPdebugMessage(" -> deleting set-packing constraint <%s>, all variable are fixed to zero\n", SCIPconsGetName(cons));
      assert(SCIPconsIsActive(cons));
      SCIP_CALL( SCIPdelCons(scip, cons) );
      ++(*ndelconss);

      return SCIP_OKAY;
   }

   /* all but one variable were fixed to zero then delete the constraint and for setpartition fix the remaining variable to 1 */
   if( consdata->nfixedzeros + 1 == nvars )
   {
      assert(consdata->nfixedones == 0);

      /* all variables except one are fixed to zero:
       * - a set packing constraint is feasible anyway, and can be deleted
       * - a set partitioning or covering constraint is feasible and can be deleted after the
       *   remaining variable is fixed to one
       */
      if( consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING || consdata->setppctype == SCIP_SETPPCTYPE_COVERING ) /*lint !e641*/
      {
	 fixed = FALSE;
	 for( v = nvars - 1; v >= 0; --v )
	 {
	    assert(SCIPvarGetLbLocal(vars[v]) < 0.5);
	    if( SCIPvarGetUbLocal(vars[v]) > 0.5 )
	    {
	       SCIPdebugMessage("trying to fix <%s> to 1 due to it's the last unfixed variable is the set-partitioning/covering constraint\n", SCIPvarGetName(vars[v]));

	       /* fix the remaining set partition variable */
	       SCIP_CALL( SCIPfixVar(scip, vars[v], 1.0, cutoff, &fixed) );
	       if( *cutoff )
	       {
                  SCIPdebugMessage("setppc constraint <%s>: infeasible fixing <%s> == 1\n",
                     SCIPconsGetName(cons), SCIPvarGetName(vars[v]));

		  return SCIP_OKAY;
	       }

	       assert(fixed);
	       ++(*nfixedvars);
	       break;
	    }
	 }
	 assert(fixed);
      }

      /* delete constraint */
      SCIPdebugMessage(" -> deleting constraint <%s>, all variable are fixed\n", SCIPconsGetName(cons));
      assert(SCIPconsIsActive(cons));
      SCIP_CALL( SCIPdelCons(scip, cons) );
      ++(*ndelconss);

      return SCIP_OKAY;
   }

   /* all but two variable were fixed to zero in a setpartitioning constraint then delete the constraint and
    * aggregate the remaining two variables
    */
   if( consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING && consdata->nfixedzeros + 2 == nvars ) /*lint !e641*/
   {
      SCIP_VAR* var;

      var = NULL;
      for( v = nvars - 1; v >= 0; --v )
      {
	 assert(SCIPvarGetLbLocal(vars[v]) < 0.5);

	 if( SCIPvarGetUbLocal(vars[v]) > 0.5 )
	 {
	    if( var == NULL )
	       var = vars[v];
	    else
	    {
	       SCIP_Bool redundant;
	       SCIP_Bool aggregated;
#ifdef VARUSES
	       SCIP_CONSHDLR* conshdlr;
	       SCIP_CONSHDLRDATA* conshdlrdata;

	       /* get event handler and event handler data */
	       conshdlr = SCIPconsGetHdlr(cons);
	       assert(conshdlr != NULL);
	       conshdlrdata = SCIPconshdlrGetData(conshdlr);
	       assert(conshdlrdata != NULL);
#endif
	       if( aggregate )
	       {
		  SCIPdebugMessage("trying to aggregate <%s> and <%s> due to they are the last two unfixed variables in the set partitionning constraint <%s>\n", SCIPvarGetName(var), SCIPvarGetName(vars[v]), SCIPconsGetName(cons));

#ifdef VARUSES
		  /* in order to not mess up the variable usage counting, we have to decrease usage counting, aggregate,
		   * and increase usage counting again
		   */
		  SCIP_CALL( conshdlrdataDecVaruses(scip, conshdlrdata, var) );
		  SCIP_CALL( conshdlrdataDecVaruses(scip, conshdlrdata, vars[v]) );
#endif

		  /* aggregate last remaining variables in the set partitioning constraint */
		  SCIP_CALL( SCIPaggregateVars(scip, var, vars[v], 1.0, 1.0, 1.0, cutoff, &redundant, &aggregated) );
		  if( *cutoff )
		  {
		     SCIPdebugMessage("set partitioning constraint <%s>: aggregate <%s> + <%s> == 1\n",
			SCIPconsGetName(cons), SCIPvarGetName(var), SCIPvarGetName(vars[v]));

		     return SCIP_OKAY;
		  }

#ifdef VARUSES
		  /* increase variable usage counting again */
		  SCIP_CALL( conshdlrdataIncVaruses(scip, conshdlrdata, var) );
		  SCIP_CALL( conshdlrdataIncVaruses(scip, conshdlrdata, vars[v]) );
#endif

		  if( aggregated )
		     ++(*naggrvars);

		  if( redundant )
		  {
		     /* delete constraint */
		     SCIPdebugMessage(" -> deleting constraint <%s>, all variable are fixed\n", SCIPconsGetName(cons));
		     assert(SCIPconsIsActive(cons));
		     SCIP_CALL( SCIPdelCons(scip, cons) );
		     ++(*ndelconss);
		  }
	       }
	       else
	       {
		  assert(undoneaggrvars != NULL);
		  assert(undoneaggrtypes != NULL);
		  assert(naggregations != NULL);
		  assert(saggregations != NULL);

		  SCIPdebugMessage("memorize the aggregation of <%s> + <%s> = 1, because they are the last two unfixed variable in the set partitioning constraints <%s>\n", SCIPvarGetName(var), SCIPvarGetName(vars[v]), SCIPconsGetName(cons));

		  /* resize the aggregation arrays if necessary */
		  if( *saggregations == *naggregations )
		  {
		     *saggregations = SCIPcalcMemGrowSize(scip, *naggregations + 1);
		     assert(*saggregations > *naggregations);
		     SCIP_CALL( SCIPreallocBufferArray(scip, &undoneaggrtypes, *saggregations) );
		     SCIP_CALL( SCIPreallocBufferArray(scip, &undoneaggrvars, 2 * (*saggregations)) );

		     /* clear the aggregation type array to set the default to the aggregation of the form x + y = 1 */
		     BMSclearMemoryArray(&(undoneaggrtypes[*naggregations]), *saggregations - *naggregations); /*lint !e866*/
		  }

		  /* memorize aagregation variables*/
		  assert(undoneaggrtypes[*naggregations] == FALSE);
		  undoneaggrvars[2 * (*naggregations)] = var;
		  undoneaggrvars[2 * (*naggregations) + 1] = vars[v];
		  ++(*naggregations);

		  if( !SCIPdoNotAggr(scip) )
		  {
		     /* delete constraint */
		     SCIPdebugMessage(" -> deleting constraint <%s>, all variable are fixed\n", SCIPconsGetName(cons));
		     assert(SCIPconsIsActive(cons));
		     SCIP_CALL( SCIPdelCons(scip, cons) );
		     ++(*ndelconss);
		  }
	       }

	       return SCIP_OKAY;
	    }
	 }
      }
      /* we should never be here, because the last to unfixed variables should have been either aggregated or a cutoff
       * should be applied
       */
      assert(FALSE);
   }

   return SCIP_OKAY;
}


/** check for overlapping constraint */
static
SCIP_RETCODE checkForOverlapping(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< constraint which may overlap */
   int const             considx,            /**< constriant index to avoid checking against itself */
   int const             endidx,             /**< end index to check against given constraint */
   SCIP_CONS**const      usefulconss,        /**< clique constraints */
   int const             nusefulconss,       /**< number of clique constraints */
   SCIP_VAR**const       usefulvars,         /**< storage for all found variables */
   int*const             nusefulvars,        /**< pointer to store number of added variables */
   SCIP_HASHMAP*const    vartoindex,         /**< hashmap mapping variables to indices */
   int*const             varnconss,          /**< storage for remembering the number of constraints a variable occurs */
   int*const             maxnvarconsidx,     /**< storage for the maximal number of occurances of a variable */
   int**const            varconsidxs,        /**< storage for constraint indices in which the corresponding variable exists */
   int*const             countofoverlapping, /**< the amount of variables of cons which overlap in all other constraint */
   SCIP_Bool const       shrinking,          /**< try to replace some variables with one variable */
   SCIP_Bool*const       chgcons,            /**< pointer to store if the given constraint was changed, due to
					      *   added/deleted variables
					      */
   SCIP_VAR**            undoneaggrvars,     /**< array to store aggregation variables, if aggregation is not performed
					      *   yet; both variables are standing next to each other;
					      */
   SCIP_Bool*            undoneaggrtypes,    /**< array to store aggregation type, if aggregation is not performed yet;
					      *   type FALSE means the aggregation is of the form x + y = 1; type TRUE means
					      *   the aggregation is of the form x = y;
					      */
   int*const             naggregations,      /**< pointer to store number of aggregations which are not yet performed; */
   int*const             saggregations,      /**< pointer to store size of the array for aggregation type and two times
					      *   the value is the size of the array for the aggregation variables which
					      *   are not yet performed;
					      */
   int*const             nfixedvars,         /**< pointer to count number of deleted variables */
   int*const             naggrvars,          /**< pointer to count number of aggregated variables */
   int*const             nchgcoefs,          /**< pointer to count number of changed coefficients */
   int*const             ndelconss,          /**< pointer to count number of deleted constraints */
   SCIP_Bool*const       cutoff              /**< pointer to store if the problem is infeasible due to a fixing */
   )
{
   SCIP_CONS* cons1;
   SCIP_CONSDATA* consdata1;
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_VAR** vars1;
   SCIP_VAR* var;
   SCIP_VAR* var1;
   SCIP_Bool fixed;
   SCIP_Bool overlapdestroyed;
   int nvars;
   int nvars1;
   int oldnfixedzeros;
   int c;
   int v;
   int v1;
#ifndef NDEBUG
   int oldnaggrvars;
#endif

   assert(scip != NULL);
   assert(cons != NULL);
   assert(usefulconss != NULL && nusefulconss > 0);
   assert(0 <= considx && considx < nusefulconss);
   assert(usefulconss[considx] == cons);
   assert(0 <= endidx && endidx <= nusefulconss);
   assert(countofoverlapping != NULL);
   assert(chgcons != NULL);
   assert(undoneaggrvars != NULL);
   assert(undoneaggrtypes != NULL);
   assert(naggregations != NULL);
   assert(saggregations != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(nchgcoefs != NULL);
   assert(ndelconss != NULL);
   assert(cutoff != NULL);

   if( !SCIPconsIsActive(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;

   if( nvars == 0 )
      return SCIP_OKAY;

   vars = consdata->vars;
   assert(vars != NULL);

   oldnfixedzeros = consdata->nfixedzeros;
   overlapdestroyed = FALSE;

   /* first check for redundancy for all unprocessed constraints with cons */
   for( c = endidx - 1; c >= 0; --c )
   {
      cons1 = usefulconss[c];

      if( !SCIPconsIsActive(cons1) )
	 continue;

      /* avoid checking constraint against itself */
      if( considx == c )
	 continue;

      assert(usefulconss[c] != cons);

#ifndef NDEBUG
      oldnaggrvars = *naggrvars;
#endif

      /* check if constraint is already redundant or infeasible due to fixings, fix or aggregate left over variables if
       * possible
       */
      SCIP_CALL( presolvePropagateCons(scip, cons1, FALSE, undoneaggrvars, undoneaggrtypes, naggregations, saggregations, nfixedvars, naggrvars, ndelconss, cutoff) );

      if( *cutoff )
	 return SCIP_OKAY;

      /* we can't handle aggregated variables later on so we should have saved them for later */
      assert(*naggrvars == oldnaggrvars);

      if( !SCIPconsIsActive(cons1) )
	 continue;

      consdata1 = SCIPconsGetData(cons1);
      assert(consdata1 != NULL);

      nvars1 = consdata1->nvars;

      if( nvars1 == 0 )
	 continue;

      /* no more variables from cons as nvars1 can overlap */
      assert(countofoverlapping[c] <= nvars1);

      /* constraint should not be redundant or infeasible */
      assert(consdata1->nfixedones == 0);

      SCIPdebugMessage("constraint <%s> overlaps with constraint <%s> by %d variables\n", SCIPconsGetName(cons), SCIPconsGetName(cons1), countofoverlapping[c]);

      /* cons1 includes cons */
      if( !overlapdestroyed && countofoverlapping[c] == nvars - consdata->nfixedzeros )
      {
	 if( consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING ) /*lint !e641*/
	 {
	    if( nvars - consdata->nfixedzeros < nvars1 )
	    {
#ifndef NDEBUG
	       SCIP_Bool negated0;
	       SCIP_Bool negated1;
#endif

	       /* both constraints should stay merged */
	       assert(consdata->merged);
	       assert(consdata1->merged);

	       vars1 = consdata1->vars;
	       assert(vars1 != NULL);

	       /* sorting array after indices of variables, negated and active counterparts would stand side by side */
	       SCIPsortDownPtr((void**)vars1, SCIPvarCompActiveAndNegated, nvars1);
	       /* standard setppc-sorting now lost */
	       consdata1->sorted = FALSE;

	       /* iterate over the both cliques variables the "same" time */
	       for( v = nvars - 1, v1 = nvars1 - 1; v >= 0 && v1 >= 0; )
	       {
		  if( SCIPvarGetLbLocal(vars1[v1]) > 0.5 || SCIPvarGetUbLocal(vars1[v1]) < 0.5 )
		  {
		     --v1;
		     continue;
		  }
		  if( SCIPvarGetLbLocal(vars[v]) > 0.5 || SCIPvarGetUbLocal(vars[v]) < 0.5 )
		  {
		     --v;
		     continue;
		  }

		  /* all variables inside the second clique constraint should be either active or negated of an active one */
		  assert(SCIPvarIsActive(vars1[v1]) || (SCIPvarGetStatus(vars1[v1]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars1[v1]))));

		  /* get not negated variable and clique value in cons */
		  if( SCIPvarGetStatus(vars[v]) != SCIP_VARSTATUS_NEGATED )
		  {
		     var = vars[v];
#ifndef NDEBUG
		     negated0 = FALSE;
#endif
		  }
		  else
		  {
		     var = SCIPvarGetNegationVar(vars[v]);
#ifndef NDEBUG
		     negated0 = TRUE;
#endif
		  }

		  /* get active variable and clique value of next variable */
		  if( SCIPvarIsActive(vars1[v1]) )
		  {
		     var1 = vars1[v1];
#ifndef NDEBUG
		     negated1 = FALSE;
#endif
		  }
		  else
		  {
		     assert(SCIPvarGetStatus(vars1[v1]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars1[v1])));
		     var1 = SCIPvarGetNegationVar(vars1[v1]);
#ifndef NDEBUG
		     negated1 = TRUE;
#endif
		  }

		  /* variable index in the constraint smaller than the other one, so go to the next variable in cons */
		  if( SCIPvarGetIndex(var) < SCIPvarGetIndex(var1)  )
		     --v;
		  /* variable index in the constraint is greater than the other one, so fix this variable */
		  else if( SCIPvarGetIndex(var) > SCIPvarGetIndex(var1)  )
		  {
		     SCIPdebugMessage("trying to fix <%s> to 0 because it is in the same clique with a complete set partitioning constraint\n", SCIPvarGetName(vars1[v1]));

		     /* fix all variables except the one which has the negated var in the clique to zero */
		     SCIP_CALL( SCIPfixVar(scip, vars1[v1], 0.0, cutoff, &fixed) );
		     if( *cutoff )
		     {
			SCIPdebugMessage("fixing led to cutoff\n");

			return SCIP_OKAY;
		     }

		     assert(fixed);
		     ++(*nfixedvars);
		     --v1;
		  }
		  else
		  {
		     /* because the constraint's are merged it is not possible that one constraint contains a negated
		      * variable of another and because all variables in cons are in cons1 this should be really the
		      * same variable here; so we can decrease v and v1
		      */
		     assert(negated0 == negated1);

		     --v;
		     --v1;
		  }
	       }
	       /* maybe we ended because of cons(v reached -1) so try to add rest of cons1 to cons */
	       for( ; v1 >= 0; --v1)
	       {
		  if( SCIPvarGetLbLocal(vars1[v1]) > 0.5 || SCIPvarGetUbLocal(vars1[v1]) < 0.5 )
		     continue;

		  SCIPdebugMessage("trying to fix <%s> to 0 because it is in the same clique with a complete set partitioning constraint\n", SCIPvarGetName(vars1[v1]));

		  /* fix all variables except the one which has the negated var in the clique to zero */
		  SCIP_CALL( SCIPfixVar(scip, vars1[v1], 0.0, cutoff, &fixed) );
		  if( *cutoff )
		  {
		     SCIPdebugMessage("fixing led to cutoff\n");

		     return SCIP_OKAY;
		  }

		  assert(fixed);
		  ++(*nfixedvars);
	       }
	    }

	    /* if caused by all fixings now this set partitioning constraint doesn't have any variable which was
	     * fixed to one, it's infeasible */
	    if( consdata1->setppctype == SCIP_SETPPCTYPE_PARTITIONING && consdata1->nfixedzeros == nvars1 && consdata1->nfixedones != 1 ) /*lint !e641*/
	    {
	       SCIPdebugMessage("all variable in the set-partitioning constraint <%s> are fixed to zero, this leads to a cutoff\n", SCIPconsGetName(cons1));
	       *cutoff = TRUE;

	       return SCIP_OKAY;
	    }

	    assert(SCIPconsIsActive(cons1));
	    /* delete second constraint */
	    SCIPdebugMessage(" -> deleting constraint <%s> number <%d> because it includes the setpartitioning constraint <%s> number <%d>\n", SCIPconsGetName(cons1), c, SCIPconsGetName(cons), considx);

	    SCIP_CALL( SCIPdelCons(scip, cons1) );
	    ++(*ndelconss);
	 }
	 /* could already be deleted because the constraint was included in another set partition constraint */
	 else if( SCIPconsIsActive(cons) )
	 {
	    /* delete cons due to redundancy to cons1 */
	    SCIPdebugMessage(" -> deleting constraint <%s> number <%d> due to inclusion in constraint <%s> number <%d>\n", SCIPconsGetName(cons), considx, SCIPconsGetName(cons1), c);
	    SCIP_CALL( SCIPdelCons(scip, cons) );
	    ++(*ndelconss);
	 }
      }
      /* cons includes cons1
       *
       * @note that zero fixations from above can only appear through a set-partitioning constraint, this means if
       *       cons was the set-partitioning constraint only variables which are not in this constraint could be fixed
       *       to zero, and this also means that the overlapping variables in this particular case are still active or
       *       fixed to 1
       *       later on it could be possible that even variables in cons are fixed to zero, which can lead to wrong
       *       results when checking if countofoverlapping[c] + consdata1->nfixedzeros == nvars1, because a fixed
       *       variable could be counted twice
       */
      else if( (!overlapdestroyed && countofoverlapping[c] + consdata1->nfixedzeros == nvars1) || countofoverlapping[c] == nvars1 )
      {
	 /* even in deleted constraints we may fix unfixed variables */
	 if( consdata1->setppctype == SCIP_SETPPCTYPE_PARTITIONING ) /*lint !e641*/
	 {
	    const int oldnfixedvars = *nfixedvars;
#ifndef NDEBUG
	    SCIP_Bool negated0;
	    SCIP_Bool negated1;
#endif
	    /* both constraints should stay merged */
	    assert(consdata->merged);
	    assert(consdata1->merged);

	    vars1 = consdata1->vars;

	    /* sorting array after indices of variables, negated and active counterparts would stand side by side */
	    SCIPsortDownPtr((void**)vars1, SCIPvarCompActiveAndNegated, nvars1);
	    /* standard setppc-sorting now lost */
	    consdata1->sorted = FALSE;

	    /* iterate over the both cliques variables the "same" time */
	    for( v = nvars - 1, v1 = nvars1 - 1; v >= 0 && v1 >= 0; )
	    {
	       if( SCIPvarGetLbLocal(vars1[v1]) > 0.5 || SCIPvarGetUbLocal(vars1[v1]) < 0.5 )
	       {
		  --v1;
		  continue;
	       }
	       if( SCIPvarGetLbLocal(vars[v]) > 0.5 || SCIPvarGetUbLocal(vars[v]) < 0.5 )
	       {
		  --v;
		  continue;
	       }

	       /* all variables inside the second clique constraint should be either active or negated of an active one */
	       assert(SCIPvarIsActive(vars1[v1]) || (SCIPvarGetStatus(vars1[v1]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars1[v1]))));
	       /* all variables inside the first clique constraint should be either active or negated of an active one */
	       assert(SCIPvarIsActive(vars[v]) || (SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars[v]))));

	       /* get not negated variable and clique value in cons */
	       if( SCIPvarIsActive(vars[v]) )
	       {
		  var = vars[v];
#ifndef NDEBUG
		  negated0 = FALSE;
#endif
	       }
	       else
	       {
		  assert(SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars[v])));
		  var = SCIPvarGetNegationVar(vars[v]);
#ifndef NDEBUG
		  negated0 = TRUE;
#endif
	       }

	       /* get active variable and clique value of next variable */
	       if( SCIPvarIsActive(vars1[v1]) )
	       {
		  var1 = vars1[v1];
#ifndef NDEBUG
		  negated1 = FALSE;
#endif
	       }
	       else
	       {
		  assert(SCIPvarGetStatus(vars1[v1]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars1[v1])));
		  var1 = SCIPvarGetNegationVar(vars1[v1]);
#ifndef NDEBUG
		  negated1 = TRUE;
#endif
	       }

	       /* variable index in the constraint smaller than the other one, so go to the next variable in cons */
	       if( SCIPvarGetIndex(var) < SCIPvarGetIndex(var1)  )
	       {
		  SCIPdebugMessage("trying to fix <%s> to 0 because it is in the same clique with a complete set partitioning constraint\n", SCIPvarGetName(var));

		  /* fix all variables except the one which has the negated var in the clique to zero */
		  SCIP_CALL( SCIPfixVar(scip, vars[v], 0.0, cutoff, &fixed) );
		  if( *cutoff )
		  {
		     SCIPdebugMessage("fixing led to cutoff\n");

		     return SCIP_OKAY;
		  }

		  assert(fixed);
		  ++(*nfixedvars);

		  --v;
	       }
	       /* variable index in the constraint is greater than the other one, so fix this variable */
	       else if( SCIPvarGetIndex(var) > SCIPvarGetIndex(var1)  )
		  --v1;
	       else
	       {
		  /* because the constraint's are merged it is not possible that one constraint contains a negated
		   * variable of another and because all variables in cons1 are in cons this should be really the same
		   * variable here; so we can decrease v and v1
		   */
		  assert(negated0 == negated1);

		  --v;
		  --v1;
	       }
	    }

	    /* maybe we ended because of cons1(v1 reached -1) so try to add rest of cons to cons1 */
	    for( ; v >= 0; --v)
	    {
	       if( SCIPvarGetLbLocal(vars[v]) > 0.5 || SCIPvarGetUbLocal(vars[v]) < 0.5 )
		  continue;

	       SCIPdebugMessage("trying to fix <%s> to 0 because it is in the same clique with a complete set partitioning constraint\n", SCIPvarGetName(vars[v]));

	       /* fix all variables except the one which has the negated var in the clique to zero */
	       SCIP_CALL( SCIPfixVar(scip, vars[v], 0.0, cutoff, &fixed) );
	       if( *cutoff )
	       {
		  SCIPdebugMessage("fixing led to cutoff\n");

		  return SCIP_OKAY;
	       }

	       assert(fixed);
	       ++(*nfixedvars);
	    }

	    /* if caused by all fixings now this set partitioning constraint doesn't have any variable which was
	     * fixed to one, it's infeasible */
	    if( consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING && consdata->nfixedzeros == nvars && consdata->nfixedones != 1 ) /*lint !e641*/
	    {
	       SCIPdebugMessage("all variable in the set-partitioning constraint <%s> are fixed to zero, this leads to a cutoff\n", SCIPconsGetName(cons1));
	       *cutoff = TRUE;

	       return SCIP_OKAY;
	    }

	    /* could already be deleted because the constraint was included in another set partition constraint */
	    if( SCIPconsIsActive(cons) )
	    {
	       /* delete cons because it include another set partitioning constraint */
	       SCIPdebugMessage(" -> deleting constraint <%s> number <%d> because it includes the setpartitioning constraint <%s> number <%d>\n", SCIPconsGetName(cons), considx, SCIPconsGetName(cons1), c);
	       assert(SCIPconsIsActive(cons));
	       SCIP_CALL( SCIPdelCons(scip, cons) );
	       ++(*ndelconss);
	    }

	    /* due to fixings in cons0 mark overlapping invalid for checking with fixedzero variables together */
	    if( oldnfixedvars < *nfixedvars )
	       overlapdestroyed = TRUE;
	 }
	 else
	 {
	    assert(consdata1->setppctype == SCIP_SETPPCTYPE_PACKING); /*lint !e641*/

	    /* delete cons1 due to redundancy to cons */
	    SCIPdebugMessage(" -> deleting constraint <%s> number <%d> due to inclusion in constraint <%s> number <%d>\n", SCIPconsGetName(cons1), c, SCIPconsGetName(cons), considx);
	    assert(SCIPconsIsActive(cons1));
	    SCIP_CALL( SCIPdelCons(scip, cons1) );
	    ++(*ndelconss);
	 }
      }
      /* if cons has only one unfixed variable which is not in cons1 and cons1 has one variable which does not appaer in
       * cons and both constraints are setpartitioning constraints we might aggregate both not overlapping variables and
       * delete one constraint
       */
      else if( consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING && consdata1->setppctype == SCIP_SETPPCTYPE_PARTITIONING && countofoverlapping[c] == nvars - oldnfixedzeros - 1 && countofoverlapping[c] == nvars1 - 1 ) /*lint !e641*/
      {
	 SCIP_VAR* aggvar1;
	 SCIP_VAR* aggvar2;
	 SCIP_Bool negated0;
	 SCIP_Bool negated1;

	 aggvar1 = NULL;
	 aggvar2 = NULL;

	 /* both constraints should stay merged */
	 assert(consdata->merged);
	 assert(consdata1->merged);

	 vars1 = consdata1->vars;

	 /* sorting array after indices of variables, negated and active counterparts would stand side by side */
	 SCIPsortDownPtr((void**)vars1, SCIPvarCompActiveAndNegated, nvars1);
	 /* standard setppc-sorting now lost */
	 consdata1->sorted = FALSE;

	 /* iterate over the both cliques variables the "same" time */
	 for( v = nvars - 1, v1 = nvars1 - 1; v >= 0 && v1 >= 0; )
	 {
	    if( SCIPvarGetLbLocal(vars1[v1]) > 0.5 || SCIPvarGetUbLocal(vars1[v1]) < 0.5 )
	    {
	       --v1;
	       continue;
	    }
	    if( SCIPvarGetLbLocal(vars[v]) > 0.5 || SCIPvarGetUbLocal(vars[v]) < 0.5 )
	    {
	       --v;
	       continue;
	    }

	    /* all variables inside the second clique constraint should be either active or negated of an active one */
	    assert(SCIPvarIsActive(vars1[v1]) || (SCIPvarGetStatus(vars1[v1]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars1[v1]))));
	    /* all variables inside the first clique constraint should be either active or negated of an active one */
	    assert(SCIPvarIsActive(vars[v]) || (SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars[v]))));

	    /* get not negated variable and clique value in cons */
	    if( SCIPvarIsActive(vars[v]) )
	    {
	       var = vars[v];
	       negated0 = FALSE;
	    }
	    else
	    {
	       assert(SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars[v])));
	       var = SCIPvarGetNegationVar(vars[v]);
	       negated0 = TRUE;
	    }

	    /* get active variable and clique value of next variable */
	    if( SCIPvarIsActive(vars1[v1]) )
	    {
	       var1 = vars1[v1];
	       negated1 = FALSE;
	    }
	    else
	    {
	       assert(SCIPvarGetStatus(vars1[v1]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars1[v1])));
	       var1 = SCIPvarGetNegationVar(vars1[v1]);
	       negated1 = TRUE;
	    }

	    /* variable index in the constraint smaller than the other one, so go to the next variable in cons */
	    if( SCIPvarGetIndex(var) < SCIPvarGetIndex(var1)  )
	    {
	       assert(aggvar1 == NULL);
	       aggvar1 = vars[v];

	       if( aggvar2 != NULL )
		  break;

	       --v;
	    }
	    /* variable index in the constraint is greater than the other one, so fix this variable */
	    else if( SCIPvarGetIndex(var) > SCIPvarGetIndex(var1)  )
	    {
	       assert(aggvar2 == NULL);
	       aggvar2 = vars1[v1];

	       if( aggvar1 != NULL )
		  break;

	       --v1;
	    }
	    else
	    {
	       /* because the constraint's are merged it is not possible that one constraint contains a negated variable
		* of another, but both variables in both constraints still can be negated to each other
		*/
	       if( negated0 != negated1 )
	       {
		  /* cons is except for one variable equal to cons1 and the unequal variable in cons is negated
		   * to the one in cons1, so the problem is infeasible
		   */
		  SCIPdebugMessage("two set-partitioning constraint <%s> and <%s> have only one variable not in common, but this variable <%s> appears in one constraint as the negated version as in the other constraint\n", SCIPconsGetName(cons), SCIPconsGetName(cons1), SCIPvarGetName(vars[v]));
		  *cutoff = TRUE;

		  return SCIP_OKAY;
	       }
	       --v;
	       --v1;
	    }
	 }

	 /* due to fixings, it is possible that there are no active variables left, we we did not recognize which variables we could aggregate */
	 if( aggvar1 == NULL && aggvar2 == NULL )
	    continue;

	 /* determine second aggregation var, if not yet done */
	 if( aggvar2 == NULL )
	 {
	    for( ; v1 >= 0; --v1)
	    {
	       if( SCIPvarGetLbLocal(vars1[v1]) > 0.5 || SCIPvarGetUbLocal(vars1[v1]) < 0.5 )
		  continue;

	       aggvar2 = vars1[v1];
	       break;
	    }
	 }
	 /* determine first aggregation var, if not yet done */
	 else if( aggvar1 == NULL )
	 {
	    /* maybe we ended because of cons1(v1 reached -1) so find the aggvar1 in cons */
	    for( ; v >= 0; --v)
	    {
	       if( SCIPvarGetLbLocal(vars[v]) > 0.5 || SCIPvarGetUbLocal(vars[v]) < 0.5 )
		  continue;

	       aggvar1 = vars[v];
	       break;
	    }
	 }

	 /* due to fixings, it is possible that there are no active variables left, we we did not recognize which variables we could aggregate */
	 if( aggvar1 == NULL || aggvar2 == NULL )
	    continue;

	 SCIPdebugMessage("memorize the aggregation of <%s> == <%s>, because they are the last two variable which are different in these two set partitioning constraints <%s> <%s>\n", SCIPvarGetName(aggvar1), SCIPvarGetName(aggvar2), SCIPconsGetName(cons), SCIPconsGetName(cons1));

	 /* resize the aggregation arrays if necessary */
	 if( *saggregations == *naggregations )
	 {
	    *saggregations = SCIPcalcMemGrowSize(scip, *naggregations + 1);
	    assert(*saggregations > *naggregations);
	    SCIP_CALL( SCIPreallocBufferArray(scip, &undoneaggrtypes, *saggregations) );
	    SCIP_CALL( SCIPreallocBufferArray(scip, &undoneaggrvars, 2 * (*saggregations)) );

	    /* clear the aggregation type array to set the default to the aggregation of the form x + y = 1 */
	    BMSclearMemoryArray(&(undoneaggrtypes[*naggregations]), *saggregations - *naggregations); /*lint !e866*/
	 }

	 /* memorize aagregation variables*/
	 undoneaggrtypes[*naggregations] = TRUE;
	 undoneaggrvars[2 * (*naggregations)] = aggvar1;
	 undoneaggrvars[2 * (*naggregations) + 1] = aggvar2;
	 ++(*naggregations);

	 if( !SCIPdoNotAggr(scip) )
	 {
	    /* delete constraint */
	    SCIPdebugMessage(" -> deleting constraint <%s> number <%d> because it is dominated by constraint <%s>\n", SCIPconsGetName(cons1), c, SCIPconsGetName(cons));
	    assert(SCIPconsIsActive(cons1));
	    SCIP_CALL( SCIPdelCons(scip, cons1) );
	    ++(*ndelconss);
	 }
      }
      /* w.l.o.g. cons is a setpartitioning constraint and countofoverlapping == nvars - oldnfixedzeros - 1 we can
       * delete all overlapping variables in cons1 and add the negated variable of the not overlapped variable to cons
       * 1; the result should be a shorter constraint with the same impact
       */
      else if( shrinking && !overlapdestroyed && countofoverlapping[c] > 1 && ((consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING && countofoverlapping[c] == nvars - oldnfixedzeros - 1) || (consdata1->setppctype == SCIP_SETPPCTYPE_PARTITIONING && countofoverlapping[c] == nvars1 - 1)) ) /*lint !e641*/
      {
	 SCIP_CONSDATA* consdatachange;
	 SCIP_VAR** varstostay;
	 SCIP_VAR** varstochange;
	 SCIP_CONS* constochange;
	 SCIP_CONS* constostay;
	 SCIP_VAR* addvar;
	 SCIP_Bool negated0;
	 SCIP_Bool negated1;
	 int nvarstostay;
	 int nvarstochange;
	 int constochangeidx;
#ifndef NDEBUG
	 const int oldnchgcoefs = *nchgcoefs;
#endif

	 addvar = NULL;

	 assert((consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING) != (consdata1->setppctype == SCIP_SETPPCTYPE_PARTITIONING) || countofoverlapping[c] != nvars - 1 || countofoverlapping[c] != nvars1 - 1); /*lint !e641*/

	 /* both constraints should stay merged */
	 assert(consdata->merged);
	 assert(consdata1->merged);

	 /* sorting array after indices of variables, negated and active counterparts would stand side by side */
	 SCIPsortDownPtr((void**)(consdata1->vars), SCIPvarCompActiveAndNegated, nvars1);
	 /* standard setppc-sorting now lost */
	 consdata1->sorted = FALSE;

	 /* initialize variables */
	 if( consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING && countofoverlapping[c] == nvars - oldnfixedzeros - 1) /*lint !e641*/
	 {
	    varstostay = vars;
	    varstochange = consdata1->vars;
	    nvarstostay = nvars;
	    nvarstochange = nvars1;
	    constostay = cons;
	    constochange = cons1;
	    consdatachange = consdata1;
	    constochangeidx = c;
	 }
	 else
	 {
	    varstostay = consdata1->vars;
	    varstochange = vars;
	    nvarstostay = nvars1;
	    nvarstochange = nvars;
	    constostay = cons1;
	    constochange = cons;
	    consdatachange = consdata;
	    constochangeidx = considx;

	    *chgcons = TRUE;
	 }

	 /* iterate over the both cliques variables the "same" time, here we need the backward loop, because we
	  * delete some variables and we don not want to loose order
	  */
	 for( v = nvarstostay - 1, v1 = nvarstochange - 1; v >= 0 && v1 >= 0; )
	 {
	    if( SCIPvarGetLbLocal(varstochange[v1]) > 0.5 || SCIPvarGetUbLocal(varstochange[v1]) < 0.5 )
	    {
	       --v1;
	       continue;
	    }
	    if( SCIPvarGetLbLocal(varstostay[v]) > 0.5 || SCIPvarGetUbLocal(varstostay[v]) < 0.5 )
	    {
	       --v;
	       continue;
	    }

	    /* all variables inside the second clique constraint should be either active or negated of an active one */
	    assert(SCIPvarIsActive(varstochange[v1]) || (SCIPvarGetStatus(varstochange[v1]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(varstochange[v1]))));
	    /* all variables inside the first clique constraint should be either active or negated of an active one */
	    assert(SCIPvarIsActive(varstostay[v]) || (SCIPvarGetStatus(varstostay[v]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(varstostay[v]))));

	    /* get not negated variable and clique value in constostay */
	    if( SCIPvarIsActive(varstostay[v]) )
	    {
	       var = varstostay[v];
	       negated0 = FALSE;
	    }
	    else
	    {
	       assert(SCIPvarGetStatus(varstostay[v]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(varstostay[v])));
	       var = SCIPvarGetNegationVar(varstostay[v]);
	       negated0 = TRUE;
	    }

	    /* get active variable and clique value of in constochange*/
	    if( SCIPvarIsActive(varstochange[v1]) )
	    {
	       var1 = varstochange[v1];
	       negated1 = FALSE;
	    }
	    else
	    {
	       assert(SCIPvarGetStatus(varstochange[v1]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(varstochange[v1])));
	       var1 = SCIPvarGetNegationVar(varstochange[v1]);
	       negated1 = TRUE;
	    }

	    /* variable index in the constraint smaller than the other one, so go to the next variable in cons */
	    if( SCIPvarGetIndex(var) < SCIPvarGetIndex(var1)  )
	    {
	       assert(addvar == NULL);
	       addvar = varstostay[v];
	       --v;
	    }
	    /* variable index in the constraint is greater than the other one, so fix this variable */
	    else if( SCIPvarGetIndex(var) > SCIPvarGetIndex(var1)  )
	    {
	       --v1;
	    }
	    else
	    {
	       /* because the constraint's are merged it is not possible that one constraint contains a negated variable
		* of another, but both constraint might have a variable in neagted form of the other
		*/
	       if( negated0 != negated1 )
	       {
		  assert(addvar == NULL);

		  SCIPdebugMessage("-> trying to fix <%s> to 0 because it would exist twice in a constraint\n", SCIPvarGetName(varstochange[v1]));

		  /* fix variable to zero */
		  SCIP_CALL( SCIPfixVar(scip, varstochange[v1], 0.0, cutoff, &fixed) );
		  if( *cutoff )
		  {
		     SCIPdebugMessage("fixing led to cutoff\n");

		     return SCIP_OKAY;
		  }

		  assert(fixed);
		  ++(*nfixedvars);

		  /* the above fixing is equal to the fixation of varstostay[v] to 1, so we can call presolvePropagateCons() for consstay */
		  SCIP_CALL( presolvePropagateCons(scip, constostay, FALSE, NULL, NULL, NULL, NULL, nfixedvars, naggrvars, ndelconss, cutoff) );

		  return SCIP_OKAY;
	       }
	       else
	       {
		  /* correct local data structure, remove variable from constraint entry where it will be removed */
		  deleteCliqueDataEntry(varstochange[v1], constochangeidx, vartoindex, varnconss, varconsidxs);

		  SCIPdebugMessage(" -> deleting variable <%s> in constraint <%s> number %d, because it will be replaced\n", SCIPvarGetName(varstochange[v1]), SCIPconsGetName(constochange), constochangeidx);
		  /* delete overlapping variabes in constochange */
		  SCIP_CALL( delCoefPos(scip, constochange, v1) );
		  ++(*nchgcoefs);
	       }

	       --v;
	       --v1;
	    }
	 }
	 assert(addvar != NULL || v >= 0);
	 /* we should have removed exactly countofoverlapping[c] variables from the constochange */
	 assert(*nchgcoefs - oldnchgcoefs == countofoverlapping[c]);

	 /* determine addvar if not yet found */
	 if( addvar == NULL )
	 {
	    for( ; v >= 0; --v)
	    {
	       if( SCIPvarGetLbLocal(varstostay[v]) > 0.5 || SCIPvarGetUbLocal(varstostay[v]) < 0.5 )
		  continue;

	       /* all variables inside the first clique constraint should be either active or negated of an active one */
	       assert(SCIPvarIsActive(varstostay[v]) || (SCIPvarGetStatus(varstostay[v]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(varstostay[v]))));

	       addvar = varstostay[v];
	       break;
	    }
	 }
	 assert(addvar != NULL);

	 /* get representative variable for all deleted variables */
	 SCIP_CALL( SCIPgetNegatedVar(scip, addvar, &addvar) );
	 assert(addvar != NULL);

	 SCIPdebugMessage(" -> adding variable <%s> to constraint <%s> number %d\n", SCIPvarGetName(addvar), SCIPconsGetName(constochange), constochangeidx);
	 /* add representative for overlapping instead */
	 SCIP_CALL( addCoef(scip, constochange, addvar) );
	 ++(*nchgcoefs);

	 /* constraint should be still merged because this added variable is new in this constraint */
	 consdatachange->merged = TRUE;
	 assert(constochangeidx == (cons == constochange ? considx : c));

	 /* correct local data structure, add constraint entry to variable data  */
	 SCIP_CALL( addCliqueDataEntry(scip, addvar, constochangeidx, TRUE, usefulvars, nusefulvars, vartoindex, varnconss, maxnvarconsidx, varconsidxs) );

	 /* cons changed so much, that it cannot be used for more overlapping checks */
	 if( *chgcons )
	    return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
}

/** try to lift variables to given constraint */
/** @todo try another variant by determine lifting variables as the intersection of all cliques variables of the
 *        constraint variables, note that the insection changes after one variable was added
 */
static
SCIP_RETCODE liftCliqueVariables(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< constraint which may overlap */
   int const             arraypos,           /**< position of constraint in global array */
   SCIP_VAR**const       usefulvars,         /**< possible variables to lift */
   int*const             nusefulvars,        /**< pointer to store number of added variables */
   int const             endidx,             /**< end index for possible lifting variables */
   SCIP_Bool**           cliquevalues,       /**< pointer to clique values of constraint-variables, either one if the
					      *   varibale is active or zero if the variable is negated
					      *   @note this array can be resized in this method
					      */
   SCIP_HASHMAP*const    vartoindex,         /**< hashmap mapping variables to indices */
   int*const             varnconss,          /**< array with number of constraints a variable occurs */
   int*const             maxnvarconsidx,     /**< array with the maximal number of occurances of a variable */
   int**const            varconsidxs,        /**< array with constraint indices in which the corresponding variable
					      *   exists
					      */
   int*const             maxnvars,           /**< pointer to store maximal number of variables of a constraint */
   int*const             nadded,             /**< pointer to store number of possible added variables */
   SCIP_Bool*const       chgcons,            /**< pointer to store if the constraint was changed, due to added
					      *   variables
					      */
   int*const             nfixedvars,         /**< pointer to count number of deleted variables */
   int*const             ndelconss,          /**< pointer to count number of deleted constraints */
   SCIP_Bool*const       cutoff              /**< pointer to store if the problem is infeasible due to a fixing */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_VAR* var1;
   SCIP_Bool fixed;
   SCIP_Bool value;
   int nvars;
   int nottocheck; /* will be the position for a variable in cons0 which is in negated form in the same clique */
   int v;
   int v1;
   int k;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(usefulvars != NULL);
   assert(cliquevalues != NULL);
   assert(*cliquevalues != NULL);
   assert(vartoindex != NULL);
   assert(varnconss != NULL);
   assert(maxnvarconsidx != NULL);
   assert(varconsidxs != NULL);
   assert(maxnvars != NULL);
   assert(nadded != NULL);
   assert(chgcons != NULL);
   assert(nfixedvars != NULL);
   assert(ndelconss != NULL);
   assert(cutoff != NULL);

   if( !SCIPconsIsActive(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;

   if( nvars == 0 )
      return SCIP_OKAY;

   assert(nvars <= *maxnvars);

   vars = consdata->vars;
   assert(vars != NULL);

   v1 = endidx;

   /* now we try to add variables with index prior to endidx to cons */
   for( v = nvars - 1; v >= 0 && v1 >= 0; )
   {
      if( SCIPvarGetLbLocal(usefulvars[v1]) > 0.5 || SCIPvarGetUbLocal(usefulvars[v1]) < 0.5 )
      {
	 --v1;
	 continue;
      }
      if( SCIPvarGetUbLocal(vars[v]) < 0.5 )
      {
	 --v;
	 continue;
      }

      /* check that constraint variables are still correctly sorted, indices of active variables should be decreasing */
      assert(v == 0 || SCIPvarCompareActiveAndNegated(vars[v], vars[v - 1]) <= 0);

      /* there should no variables fixed to one occur in our constraint */
      assert(SCIPvarGetLbLocal(vars[v]) < 0.5 && SCIPvarGetUbLocal(vars[v]) > 0.5);
      assert(SCIPvarGetLbLocal(usefulvars[v1]) < 0.5 && SCIPvarGetUbLocal(usefulvars[v1]) > 0.5);

      /* all variables which we have inside the clique constraint and which can possibly be added should be either active or negated */
      assert(SCIPvarIsActive(vars[v]) || (SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars[v]))));
      assert(SCIPvarIsActive(usefulvars[v1]) || (SCIPvarGetStatus(usefulvars[v1]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(usefulvars[v1]))));

      /* constraint should during adding of variables stay merged, because for each variable which is added holds that
       * the index of this corresponding active variable is pairwise different to all indices of all active
       * corresponding variables inside the constraint
       * @note it should not happen that we add one variable and the corresponding counterpart to the same constraint */
      assert(consdata->merged);

      /* get active variable and clique value in cons */
      if( (*cliquevalues)[v] )
	 var = vars[v];
      else
      {
	 assert(SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars[v])));
	 var = SCIPvarGetNegationVar(vars[v]);
      }

      /* get active variable and clique value of next variable */
      if( SCIPvarIsActive(usefulvars[v1]) )
      {
	 var1 = usefulvars[v1];
	 value = TRUE;
      }
      else
      {
	 assert(SCIPvarGetStatus(usefulvars[v1]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(usefulvars[v1])));
	 var1 = SCIPvarGetNegationVar(usefulvars[v1]);
	 value = FALSE;
      }

      nottocheck = -1;
      k = 0;

      /* variable index in the constraint smaller than the other one, so go to the next variable in cons */
      if( SCIPvarGetIndex(var) < SCIPvarGetIndex(var1)  )
      {
	 --v;
	 continue;
      }
      /* variable index in the constraint is greater than the other one, so check for possible inclusion of the variable */
      else if( SCIPvarGetIndex(var) > SCIPvarGetIndex(var1)  )
      {
	 assert(consdata == SCIPconsGetData(cons));

	 /* check if every variable in the actual clique is in clique with the new variable */
	 for( k = nvars - 1; k >= 0; --k )
	 {
	    if( SCIPvarGetUbLocal(vars[k]) > 0.5 )
	    {
	       /* there should no variables fixed to one occur in our constraint */
	       assert(SCIPvarGetLbLocal(vars[k]) < 0.5);
	       assert(SCIPvarIsActive(vars[k]) || (SCIPvarGetStatus(vars[k]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars[k]))));

	       if( (*cliquevalues)[k] )
	       {
		  assert(SCIPvarIsActive(vars[k]));
		  var = vars[k];
	       }
	       else
	       {
		  assert(SCIPvarGetStatus(vars[k]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars[k])));
		  var = SCIPvarGetNegationVar(vars[k]);
	       }
	       if( !SCIPhaveVarsCommonClique(scip, var1, value, var, (*cliquevalues)[k], TRUE) )
		  break;
	    }
	 }
	 --v1;
      }
      /* variable index in the constraint is equal to the index of the other variable, check if these variables are
       * negated of each other so memorize the position and check for possible inclusion of the new variable and if
       * possible decrease indices
       */
      else
      {
	 /* one clique contains the negated and the other clique the corresponding active var */
	 if( value != (*cliquevalues)[v] )
	 {
	    nottocheck = v;

	    assert(consdata == SCIPconsGetData(cons));
	    assert(nvars <= consdata->nvars);

	    /* check if every variable in the actual clique is in clique with the new variable */
	    for( k = nvars - 1; k >= 0; --k )
	    {
	       if( SCIPvarGetUbLocal(vars[k]) > 0.5 )
	       {
		  /* there should no variables fixed to one occur in our constraint */
		  assert(SCIPvarGetLbLocal(vars[k]) < 0.5);

		  assert(SCIPvarIsActive(vars[k]) || (SCIPvarGetStatus(vars[k]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars[k]))));

		  if( k == nottocheck )
		     continue;

		  if( (*cliquevalues)[k] )
		  {
		     assert(SCIPvarIsActive(vars[k]));
		     var = vars[k];
		  }
		  else
		  {
		     assert(SCIPvarGetStatus(vars[k]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars[k])));
		     var = SCIPvarGetNegationVar(vars[k]);
		  }

		  if( !SCIPhaveVarsCommonClique(scip, var1, value, var, (*cliquevalues)[k], TRUE) )
		     break;
	       }
	    }
	 }
	 /* don't decrease v because it might happen that the corresponding negated variable of var is next in
	  * usefulvars
	  */
	 --v1;
      }

      /* if k is smaller than 0 than the possible new variables is in the same clique with all variables of cons,
       * so we add the new variable to clique constraint or fix some variables */
      if( k < 0 )
      {
	 ++(*nadded);

	 /* we found a variable which is the negated variable of another one in this clique so we can fix all
	  * other variable to zero and if it's a partitioning constraint we can also fix the variable of the
	  * negated to one and we can delete the constraint too */
	 if( nottocheck >= 0 )
	 {
	    assert(consdata == SCIPconsGetData(cons));
	    assert(nvars <= consdata->nvars);
	    assert(consdata->merged);

	    /* process all vars for possible fixing */
	    for( k = consdata->nvars - 1; k >= 0; --k )
	    {
	       if( SCIPvarGetUbLocal(vars[k]) > 0.5 )
	       {
		  /* there should no variables fixed to one occur in our constraint */
		  assert(SCIPvarGetLbLocal(vars[v]) < 0.5);

		  assert(SCIPvarIsActive(vars[k]) || (SCIPvarGetStatus(vars[k]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars[k]))));

		  if( k != nottocheck )
		  {
		     SCIPdebugMessage("trying to fix <%s> to 0 because we could lift a negated variable of another constraint variable\n", SCIPvarGetName(vars[k]));
		     /* fix variable to zero */
		     SCIP_CALL( SCIPfixVar(scip, vars[k], 0.0, cutoff, &fixed) );

		     if( *cutoff )
		     {
			SCIPdebugMessage("fixing led to cutoff\n");

			return SCIP_OKAY;
		     }

		     assert(fixed);

		     ++(*nfixedvars);
		  }
	       }
	    }
	    if( consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING ) /*lint !e641*/
	    {
	       assert(SCIPvarIsActive(vars[nottocheck]) || (SCIPvarGetStatus(vars[nottocheck]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars[nottocheck]))));

	       SCIPdebugMessage("trying to fix <%s> to 1 due to this setpartitioning variable is with its negated in the same clique\n", SCIPvarGetName(vars[nottocheck]));
	       /* fix the remaining variable to one, due to it's the only one left to satisfy the constraint */
	       SCIP_CALL( SCIPfixVar(scip, vars[nottocheck], 1.0, cutoff, &fixed) );
	       if( *cutoff )
	       {
		  SCIPdebugMessage("fixing led to cutoff\n");

		  return SCIP_OKAY;
	       }

	       assert(fixed);
	       ++(*nfixedvars);
	    }

	    /* delete constraint */
	    SCIPdebugMessage(" -> deleting constraint <%s> number <%d> due to active and negated variable in the same clique constraint\n", SCIPconsGetName(cons), arraypos);
	    assert(SCIPconsIsActive(cons));
	    SCIP_CALL( SCIPdelCons(scip, cons) );
	    ++(*ndelconss);

	    break;
	 }
	 /* we found a variable which could be added to a partitioning constraint so we can fix it to zero */
	 else if( consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING ) /*lint !e641*/
	 {
	    SCIPdebugMessage("trying to fix <%s> to 0 because this variable is in the same clique with a set partition\n", SCIPvarGetName(usefulvars[v1 + 1]));
	    /* fix variable to zero */
	    SCIP_CALL( SCIPfixVar(scip, usefulvars[v1 + 1], 0.0, cutoff, &fixed) );

	    if( *cutoff )
	    {
	       SCIPdebugMessage("fixing led to cutoff\n");

	       return SCIP_OKAY;
	    }

	    assert(fixed);

	    ++(*nfixedvars);
	 }
	 /* we have found a new variable for a set packing constraint cons, so add the found variable to the first constraint */
	 else
	 {
	    SCIP_VAR* addvar;

	    assert(SCIPconsIsActive(cons));

	    addvar = usefulvars[v1 + 1];

	    assert(SCIPvarGetLbLocal(addvar) < 0.5 && SCIPvarGetUbLocal(addvar) > 0.5);

	    /* add representative instead */
	    SCIPdebugMessage(" -> adding variable <%s> to constraint <%s> number %d\n", SCIPvarGetName(usefulvars[v1 + 1]), SCIPconsGetName(cons), arraypos);
	    SCIP_CALL( addCoef(scip, cons, addvar) );
	    assert(consdata == SCIPconsGetData(cons));
	    /* we know that this constraint stays merged but later on we have to resort */
	    consdata->merged = TRUE;

	    /* second we add the constraint index to the list of indices where this variable occurs */
	    assert(SCIPhashmapExists(vartoindex, (void*) addvar));

	    /* correct local data structure, add constraint entry to variable data  */
	    SCIP_CALL( addCliqueDataEntry(scip, addvar, arraypos, FALSE, usefulvars, nusefulvars, vartoindex, varnconss, maxnvarconsidx, varconsidxs) );

	    /* we need the new pointer to the variables, because due to adding variables it is possible that we
	     * did reallocate the variables array inside the constraint, the index v should stay the same because the
	     * added variable was inserted at the end and we are decreasing v in our for loop
	     */
	    vars = consdata->vars;
	    nvars = consdata->nvars;

	    /* we need to update our data structure */

	    /* resize clique array if necessary, due to adding variables */
	    if( (*maxnvars) < nvars )
	    {
	       while( (*maxnvars) < nvars )
		  (*maxnvars) *= 2 ;
	       SCIP_CALL( SCIPreallocBufferArray(scip, cliquevalues, (*maxnvars)) );
	    }
	    (*cliquevalues)[nvars - 1] = SCIPvarIsActive(addvar) ? TRUE : FALSE;

	    (*chgcons) = TRUE;
	 }
      }
   }

   if( !SCIPconsIsActive(cons) )
      return SCIP_OKAY;

   /* maybe we stopped because of cons(v reached -1) so try to add rest in usefulvars */
   for( ; v1 >= 0; --v1)
   {
      if( SCIPvarGetLbLocal(usefulvars[v1]) > 0.5 || SCIPvarGetUbLocal(usefulvars[v1]) < 0.5 )
	 continue;

      /* get active variable and clique value */
      if( SCIPvarIsActive(usefulvars[v1]) )
      {
	 var1 = usefulvars[v1];
	 value = TRUE;
      }
      else
      {
	 assert(SCIPvarGetStatus(usefulvars[v1]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(usefulvars[v1])));
	 var1 = SCIPvarGetNegationVar(usefulvars[v1]);
	 value = FALSE;
      }

      assert(consdata == SCIPconsGetData(cons));
      assert(nvars <= consdata->nvars);

      /* check if every variable in the actual clique is in clique with the new variable */
      for( k = nvars - 1; k >= 0; --k )
      {
	 if( SCIPvarGetUbLocal(vars[k]) > 0.5 )
	 {
	    /* there should no variables fixed to one occur in our constraint */
	    assert(SCIPvarGetLbLocal(vars[k]) < 0.5);

	    assert(SCIPvarIsActive(vars[k]) || (SCIPvarGetStatus(vars[k]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars[k]))));

	    if( (*cliquevalues)[k] )
	    {
	       assert(SCIPvarIsActive(vars[k]));
	       var = vars[k];
	    }
	    else
	    {
	       assert(SCIPvarGetStatus(vars[k]) == SCIP_VARSTATUS_NEGATED && SCIPvarIsActive(SCIPvarGetNegationVar(vars[k])));
	       var = SCIPvarGetNegationVar(vars[k]);
	    }

	    if( !SCIPvarsHaveCommonClique(var1, value, var, (*cliquevalues)[k], TRUE) )
	       break;
	 }
      }

      /* add new variable to clique constraint or fix some variables */
      if( k < 0 )
      {
	 /* we found a variable which could be added to a partitioning constraint so we can fix it to zero */
	 if( consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING ) /*lint !e641*/
	 {
	    SCIPdebugMessage("trying to fix <%s> to 0 because this variable is in the same clique with a set partition\n", SCIPvarGetName(usefulvars[v1]));

	    /* fix variable to zero */
	    SCIP_CALL( SCIPfixVar(scip, usefulvars[v1], 0.0, cutoff, &fixed) );
	    if( *cutoff )
	    {
	       SCIPdebugMessage("fixing led to cutoff\n");

	       return SCIP_OKAY;
	    }
	    assert(fixed);

	    ++(*nfixedvars);
	    ++(*nadded);
	 }
	 /* add the found variable to the first constraint */
	 else
	 {
	    SCIP_VAR* addvar;

	    assert(SCIPconsIsActive(cons));

	    addvar = usefulvars[v1];

	    assert(SCIPvarGetLbLocal(addvar) < 0.5 && SCIPvarGetUbLocal(addvar) > 0.5);

	    /* add representative instead */
	    SCIPdebugMessage(" -> adding variable <%s> to constraint <%s> number %d\n", SCIPvarGetName(addvar), SCIPconsGetName(cons), arraypos);
	    SCIP_CALL( addCoef(scip, cons, addvar) );
	    assert(consdata == SCIPconsGetData(cons));
	    /* we know that this constraint stays merged but later on we have to resort */
	    consdata->merged = TRUE;

	    /* second we add the constraint index to the list of indices where this variable occurs */
	    assert(SCIPhashmapExists(vartoindex, (void*) addvar));

	    /* correct local data structure, add constraint entry to variable data  */
	    SCIP_CALL( addCliqueDataEntry(scip, addvar, arraypos, FALSE, usefulvars, nusefulvars, vartoindex, varnconss, maxnvarconsidx, varconsidxs) );

	    /* we need the new pointer to the variables, because due to adding variables it is possible that we
	     * did reallocate the variables array inside the constraint, the index v should stay the same because the
	     * added variable was inserted at the end and we are decreasing v in our for loop
	     */
	    vars = consdata->vars;
	    nvars = consdata->nvars;

	    /* we need to update our data structure */

	    /* resize clique array if necessary, due to adding variables */
	    if( (*maxnvars) < nvars )
	    {
	       while( (*maxnvars) < nvars )
		  (*maxnvars) *= 2 ;
	       SCIP_CALL( SCIPreallocBufferArray(scip, cliquevalues, (*maxnvars)) );
	    }
	    (*cliquevalues)[nvars - 1] = SCIPvarIsActive(addvar) ? TRUE : FALSE;

	    ++(*nadded);
	    (*chgcons) = TRUE;
	 }
      }
   }

   return SCIP_OKAY;
}

/** perform all collected aggregations */
static
SCIP_RETCODE performAggregations(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_VAR**const       undoneaggrvars,     /**< aggregation variables storage */
   SCIP_Bool*const       undoneaggrtypes,    /**< aggregation type storage, type FALSE means the aggregation is of the
					      *   form x + y = 1; type TRUE means the aggregation is of the form x = y;
					      */
   int const             naggregations,      /**< number of aggregations to performed */
   int*const             naggrvars,          /**< pointer to count number of aggregated variables */
   SCIP_Bool*const       cutoff              /**< pointer to store if the problem is infeasible due to a fixing */
   )
{  /*lint --e{715}*/
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_Bool aggregated;
   SCIP_Bool redundant;
   int a;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(undoneaggrvars != NULL);
   assert(undoneaggrtypes != NULL);
   assert(naggregations > 0);
   assert(naggrvars != NULL);
   assert(cutoff != NULL);

   /* loop over all open aggreagtions and try to aggregate them */
   for( a = 0; a < naggregations; ++a  )
   {
      var1 = undoneaggrvars[2 * a];
      var2 = undoneaggrvars[2 * a + 1];
      assert(var1 != NULL);
      assert(var2 != NULL);

      SCIPdebugMessage("trying to aggregate <%s> %s <%s>%s\n", SCIPvarGetName(var1), undoneaggrtypes[a] ? "=" : "+", SCIPvarGetName(var2), undoneaggrtypes[a] ? "" : " = 1");

#ifdef VARUSES
      /* in order to not mess up the variable usage counting, we have to decrease usage counting, aggregate,
       * and increase usage counting again
       */
      SCIP_CALL( conshdlrdataDecVaruses(scip, conshdlrdata, var1) );
      SCIP_CALL( conshdlrdataDecVaruses(scip, conshdlrdata, var2) );
#endif

      /* aggregate last remaining variables in the set partitioning constraint */
      if( undoneaggrtypes[a] )
      {
	 SCIP_CALL( SCIPaggregateVars(scip, var1, var2, 1.0, -1.0, 0.0, cutoff, &redundant, &aggregated) );
      }
      else
      {
	 SCIP_CALL( SCIPaggregateVars(scip, var1, var2, 1.0, 1.0, 1.0, cutoff, &redundant, &aggregated) );
      }

      if( *cutoff )
      {
	 SCIPdebugMessage("aggregation was infeasible\n");

	 return SCIP_OKAY;
      }
      /* binary variables should always be aggregated, or due to fixation the aggregation is redundant */
      assert(redundant);

      if( aggregated )
	 ++(*naggrvars);

#ifdef VARUSES
      /* increase variable usage counting again */
      SCIP_CALL( conshdlrdataIncVaruses(scip, conshdlrdata, var1) );
      SCIP_CALL( conshdlrdataIncVaruses(scip, conshdlrdata, var2) );
#endif
   }

   return SCIP_OKAY;
}

/** check whether we can combine or grow cliques so some constraints become redundant or we can fix variables */
/** @todo try another variant, by building up the clique graph and delete unnecessary (transitive closure) edges and do
 *        a bfs search to search for common ancestors to get all possible lifting variables
 */
static
SCIP_RETCODE preprocessCliques(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONS**const      conss,              /**< constraint set */
   int const             nconss,             /**< number of constraints in constraint set */
   int const             nrounds,            /**< actual presolving round */
   int*const             firstchange,        /**< pointer to store first changed constraint */
   int*const             nfixedvars,         /**< pointer to count number of deleted variables */
   int*const             naggrvars,          /**< pointer to count number of aggregated variables */
   int*const             ndelconss,          /**< pointer to count number of deleted constraints */
   int*const             nchgcoefs,          /**< pointer to count number of deleted coefficients */
   SCIP_Bool*const       cutoff              /**< pointer to store if the problem is infeasible due to a fixing */
)
{
   /* extend cliques/constraints by checking whether some variables are in the same clique, no pairwise clique lifting
    * which would be slower
    */
   SCIP_CONS** usefulconss;                  /* array with pointers of constraint of setpartitioning and setpacking type */
   SCIP_VAR** usefulvars;                    /* array with pointers of variables in setpartitioning and setpacking constraints */
   int** varconsidxs;                        /* array consisting of constraint indices in which the corresponding variable exists */
   int* varnconss;                           /* array consisting of number of constraints the variable occurs */
   int* maxnvarconsidx;                      /* maximal number of occurances of a variable */
   int* countofoverlapping;                  /* the amount of variables which are in another constraint */
   SCIP_Bool* cliquevalues;                  /* values of clique-variables, either one if the varibale is active or zero if the variable is negated */

   SCIP_HASHMAP* vartoindex;                 /* mapping of SCIP variables to indices */
   SCIP_CONSDATA* consdata;

   SCIP_Bool useadditionalmem;
   SCIP_Bool chgcons0;
   int nvars;
   int c;
   int v;
   int nusefulconss;
   int nusefulvars;
   int susefulvars;
   int maxnvars;
   int varindex;

   SCIP_VAR** undoneaggrvars;                /* storage for not yet performed aggregations */
   SCIP_Bool* undoneaggrtypes;               /* storage for not yet performed aggregation type (x = y or x + y = 1) */
   int saggregations;
   int naggregations;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(conss != NULL || nconss == 0);
   assert(firstchange != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(ndelconss != NULL);
   assert(nchgcoefs != NULL);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   if( nconss == 0 )
      return SCIP_OKAY;

   nvars = SCIPgetNVars(scip);

   if( nvars == 0 )
      return SCIP_OKAY;

   susefulvars = 2 * nvars; /* two times because of negated vars, maybe due to deleted variables we need to increase this */

   /* a hashmap from varindex to postion in varconsidxs array, because above is still too small */
   SCIP_CALL( SCIPhashmapCreate(&vartoindex, SCIPblkmem(scip), SCIPcalcHashtableSize(5 * nvars)) );

   /* get temporary memory for the aggregation storage, to memorize aggregations which will be performed later, otherwise we would destroy our local data structures */
   saggregations = nvars;
   SCIP_CALL( SCIPallocBufferArray(scip, &undoneaggrvars, 2 * saggregations) );
   SCIP_CALL( SCIPallocBufferArray(scip, &undoneaggrtypes, saggregations) );
   BMSclearMemoryArray(undoneaggrtypes, saggregations);
   naggregations = 0;

   /* get temporary memory for all clique constraints, all appearing variables and the mapping from variables to constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &usefulconss, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &usefulvars, susefulvars) );
   BMSclearMemoryArray(usefulvars, susefulvars);
   SCIP_CALL( SCIPallocBufferArray(scip, &varnconss, susefulvars + 1) );
   BMSclearMemoryArray(varnconss, susefulvars + 1);
   SCIP_CALL( SCIPallocBufferArray(scip, &maxnvarconsidx, susefulvars + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varconsidxs, susefulvars + 1) );
   BMSclearMemoryArray(varconsidxs, susefulvars + 1);
   nusefulvars = 0;
   nusefulconss = 0;
   maxnvars = 0;

   useadditionalmem = FALSE;

   /* @todo: check for round limit for adding extra clique constraints */
   /* adding clique constraints which arises from global clique information */
   if( conshdlrdata->nclqpresolve == 0 && conshdlrdata->addvariablesascliques )
   {
      SCIP_VAR** binvars;
      int* cliquepartition;
      int ncliques;
      int nbinvars;
      int naddconss;

      nbinvars = SCIPgetNBinVars(scip);
      SCIP_CALL( SCIPduplicateBufferArray(scip, &binvars, SCIPgetVars(scip), nbinvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &cliquepartition, nbinvars) );

      /* @todo: check for better permutations/don't permutate the first round
       * @todo: take binary variables which are not of vartype SCIP_VARTYPE_BINARY into account
       */
      SCIPpermuteArray((void**)binvars, 0, nbinvars, &(conshdlrdata->randseed));

      /* try to create a clique-partition over all binary variables and create these cliques as new setppc constraints
       * and add them to the usefulconss array and adjust all necessary data this will hopefully lead to faster
       * detection of redundant constraints
       */
      SCIP_CALL( SCIPcalcCliquePartition(scip, binvars, nbinvars, cliquepartition, &ncliques) );

      /* resize usefulconss array if necessary */
      SCIP_CALL( SCIPreallocBufferArray(scip, &usefulconss, nconss + ncliques) );

      naddconss = 0;

      /* add extra clique constraints resulting from the cliquepartition calculation to SCIP and to the local data structure */
      SCIP_CALL( addExtraCliques(scip, binvars, nbinvars, cliquepartition, ncliques, usefulconss, &nusefulconss,
	    usefulvars, &nusefulvars, vartoindex, varnconss, maxnvarconsidx, varconsidxs, &maxnvars,
	    nrounds, nfixedvars, &naddconss, ndelconss, nchgcoefs, cutoff) );

      /* bad hack, we don't want to count these artificial created constraints if they got deleted, so ndelconss
       * can become negative which will be change to zero at the end of this method if it's still negative
       */
      *ndelconss -= naddconss;

      SCIPfreeBufferArray(scip, &cliquepartition);
      SCIPfreeBufferArray(scip, &binvars);

      if( *cutoff )
	 goto TERMINATE;
   }

   /* start to collect setpartitioning and setpacking constraints, and try to remove fixed variables and merged these
    * constraints
    */
   SCIP_CALL( collectCliqueConss(scip, conss, nconss, usefulconss, &nusefulconss, nfixedvars, ndelconss, nchgcoefs, cutoff) );
   /* @Note: Even after the call above some constraints can have fixed variables, because it might happen that caused by
    * mergeMultiplies some variables were fixed which occured already in previous constraints
    */
   if( *cutoff )
      goto TERMINATE;

   /* no usefulconss found */
   if( nusefulconss <= 1 )
      goto TERMINATE;

   /* @todo: maybe sort them after biggest indices too, or another variant would be to restore the order as they were
    *        read in
    */
   /* sort constraints first after type (partitioning before packing) and second after number of variables such that the
    * partitioning constraints have increasing number of variables and the packing constraints have decreasing number of
    * variables, because we loop from back to front we sort them downwards, so they are the other way around
    */
   SCIPsortDownPtr((void**)usefulconss, setppcConssSort, nusefulconss);

   /* creating all necessary data in array structure, collect all clique constraint variables and occurances */
   SCIP_CALL( collectCliqueData(scip, usefulconss, nusefulconss, usefulvars, &nusefulvars, vartoindex, varnconss, maxnvarconsidx, varconsidxs, &maxnvars) );
   assert(maxnvars > 0);

   /* allocate temporary memory for actual clique */
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquevalues, maxnvars) );
   /* allocate temporary memory for counting an overlap of variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &countofoverlapping, nusefulconss) );
   useadditionalmem = TRUE;

   /* sort usefulvars after indices of variables, negated and active counterparts will stand side by side */
   SCIPsortDownPtr((void**)usefulvars, SCIPvarCompActiveAndNegated, nusefulvars);

   /* extend cliques/constraints by checking whether some variables of a second constraint are in the same clique */
   for( c = nusefulconss - 1; c >= 0 && !SCIPisStopped(scip); --c )
   {
      SCIP_VAR** cons0vars;                  /* these are the clique variables */
      SCIP_CONS* cons0;
      int ncons0vars;
      SCIP_VAR* var0;
      int v1;
      int nadded;     /* number of possible added variables to constraint */
      int cons0fixedzeros;
      int oldnchgcoefs;
#ifndef NDEBUG
      const int oldnaggrvars = *naggrvars;
#endif
      cons0 = usefulconss[c];

      if( !SCIPconsIsActive(cons0) )
         continue;

      /* check if constraint is already redundant or infeasible due to fixings, fix or aggregate left over variables if
       * possible
       */
      SCIP_CALL( presolvePropagateCons(scip, cons0, FALSE, undoneaggrvars, undoneaggrtypes, &naggregations, &saggregations, nfixedvars, naggrvars, ndelconss, cutoff) );

      if( *cutoff )
	 break;

      /* we can't handle aggregated variables later on so we should have saved them for later */
      assert(*naggrvars == oldnaggrvars);

      if( !SCIPconsIsActive(cons0) )
         continue;

      /* we need to determine the cliquedata in each iteration because we eventual will change it later */
      consdata = SCIPconsGetData(cons0);
      assert(consdata != NULL);

      cons0vars = consdata->vars;
      ncons0vars = consdata->nvars;

      /* sorting array after indices of variables, negated and active counterparts will stand side by side */
      SCIPsortDownPtr((void**)cons0vars, SCIPvarCompActiveAndNegated, ncons0vars);
      /* standard setppc-sorting now lost */
      consdata->sorted = FALSE;

      /* clique array should be long enough */
      assert(maxnvars >= ncons0vars);

      /* clear old entries in overlapping constraint */
      BMSclearMemoryArray(countofoverlapping, nusefulconss);

      /* calculate overlapping */
      for( v = ncons0vars - 1; v >= 0 ; --v )
      {
	 var0 = cons0vars[v];

	 /* fixed variables later to the count */
	 if( SCIPvarGetLbLocal(var0) > 0.5 || SCIPvarGetUbLocal(var0) < 0.5 )
	    continue;

	 assert(SCIPhashmapExists(vartoindex, (void*) var0));

	 varindex = (int) (size_t) SCIPhashmapGetImage(vartoindex, (void*) var0);
	 for( v1 = varnconss[varindex] - 1; v1 >= 0 ; --v1 )
	    ++(countofoverlapping[varconsidxs[varindex][v1]]);
      }

      oldnchgcoefs = *nchgcoefs;
      cons0fixedzeros = consdata->nfixedzeros;

      chgcons0 = FALSE;

      /* check for overlapping constraint before starting lifting */
      SCIP_CALL( checkForOverlapping(scip, cons0, c, c, usefulconss, nusefulconss, usefulvars, &nusefulvars, vartoindex,
	    varnconss, maxnvarconsidx, varconsidxs, countofoverlapping, conshdlrdata->cliqueshrinking, &chgcons0,
	    undoneaggrvars, undoneaggrtypes, &naggregations, &saggregations,
	    nfixedvars, naggrvars, nchgcoefs, ndelconss, cutoff) );

      if( *cutoff )
	 break;

      /* we can't handle aggregated variables later on so we should have saved them for later */
      assert(*naggrvars == oldnaggrvars);

      /* if cons0 changed, we need to reorder the variables  */
      if( chgcons0 && *nchgcoefs > oldnchgcoefs )
      {
	 consdata = SCIPconsGetData(cons0);
	 assert(consdata != NULL);

	 cons0vars = consdata->vars;
	 ncons0vars = consdata->nvars;

	 /* sorting array after indices of variables, negated and active counterparts will stand side by side */
	 SCIPsortDownPtr((void**)cons0vars, SCIPvarCompActiveAndNegated, ncons0vars);
	 /* standard setppc-sorting now lost */
	 consdata->sorted = FALSE;
      }

      /* check cons0 again for redundancy/fixings, because due to fixings in all other constraints it might happen that cons0 is redundant now */
      if( consdata->nfixedones > 0 || consdata->nfixedzeros > cons0fixedzeros )
      {
	 /* check if constraint is already redundant or infeasible due to fixings, fix or aggregate left over variables if
	  * possible
	  */
	 SCIP_CALL( presolvePropagateCons(scip, cons0, FALSE, undoneaggrvars, undoneaggrtypes, &naggregations, &saggregations, nfixedvars, naggrvars, ndelconss, cutoff) );

	 if( *cutoff )
	    break;

	 /* we can't handle aggregated variables later on so we should have saved them for later */
	 assert(*naggrvars == oldnaggrvars);

	 if( !SCIPconsIsActive(cons0) )
	    continue;
      }

      nadded = 0;

      /* iterate over the cliques variables and all possible new clique variables at the "same" time, determine starting
       * index
       *
       * @note: it might be better to start the first round with our computed v1, but maybe it's better to switch to
       *        trying to add all variables the second time for set packing constraints
       */

      /* we try to add all variables to the partitioning constraints, to try to fix as much as possible */
      if( consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING ) /*lint !e641*/
         v1 = nusefulvars - 1;
      else
      {
	 /* if we already ran a presolving round we want to try to add new variables */
	 if( conshdlrdata->nclqpresolve > 0 )
	    v1 = nusefulvars - 1;
	 else
	 {
	    /* find start position of variable which we will try to add to our constraint, so we will get better clique constraints */
	    (void) SCIPsortedvecFindDownPtr((void**)usefulvars, SCIPvarCompActiveAndNegated, (void*)cons0vars[ncons0vars - 1], nusefulvars, &v1);
	    assert(v1 >= 0 && v1 < nusefulvars);
	    /* if constraint is not merged and we found a variable which is negated the same as it's neighbour we have to
	     * increase v1 to make sure that we don't loose this important variable */
	    if( v1 + 1 < nusefulvars && ((SCIPvarIsNegated(usefulvars[v1 + 1]) && SCIPvarGetNegatedVar(usefulvars[v1 + 1]) == usefulvars[v1]) || (SCIPvarIsNegated(usefulvars[v1]) && SCIPvarGetNegatedVar(usefulvars[v1]) == usefulvars[v1 + 1])) )
	       ++v1;
	 }
      }

      assert(maxnvars >= ncons0vars);
      /* initialize the cliquevalues array */
      for( v = ncons0vars - 1; v >= 0; --v )
      {
         if( SCIPvarGetLbLocal(cons0vars[v]) < 0.5 && SCIPvarGetUbLocal(cons0vars[v]) > 0.5 )
         {
            /* variable has to be either active or a negated variable of an active one */
            assert(SCIPvarIsActive(cons0vars[v]) || (SCIPvarGetStatus(cons0vars[v]) == SCIP_VARSTATUS_NEGATED &&
                  SCIPvarIsActive(SCIPvarGetNegationVar(cons0vars[v]))));
            cliquevalues[v] = SCIPvarIsActive(cons0vars[v]) ? TRUE : FALSE;
         }
      }

      chgcons0 = FALSE;

      /* try to lift variables to cons0 */
      SCIP_CALL( liftCliqueVariables(scip, cons0, c, usefulvars, &nusefulvars, v1, &cliquevalues, vartoindex, varnconss, maxnvarconsidx,
	    varconsidxs, &maxnvars, &nadded, &chgcons0, nfixedvars, ndelconss, cutoff) );

      if( *cutoff )
	 break;

      if( !SCIPconsIsActive(cons0) )
	 continue;

      /* check for redundant constraints due to changing cons0 */
      if( chgcons0 )
      {
         int i;

	 /* variables array has changed due to lifting variables, so get new values */
	 assert(consdata == SCIPconsGetData(cons0));
	 cons0vars = consdata->vars;
	 ncons0vars = consdata->nvars;

	 /* resorting array, because we added new variables, in order of indices of variables, negated
	  * and active counterparts would stand side by side
	  */
	 SCIPsortDownPtr((void**)cons0vars, SCIPvarCompActiveAndNegated, ncons0vars);
	 /* standard setppc-sorting now lost */
	 consdata->sorted = FALSE;

	 /* clear old entries in overlapping constraint */
         BMSclearMemoryArray(countofoverlapping, nusefulconss);

         for( v = ncons0vars - 1; v >= 0 ; --v )
         {
            var0 = cons0vars[v];

            /* fixed variables later to the count */
            if( SCIPvarGetLbLocal(var0) > 0.5 || SCIPvarGetUbLocal(var0) < 0.5 )
               continue;

            assert(SCIPhashmapExists(vartoindex, (void*) var0));

            varindex = (int) (size_t) SCIPhashmapGetImage(vartoindex, (void*) var0);
            for( i = varnconss[varindex] - 1; i >= 0 ; --i )
               ++(countofoverlapping[varconsidxs[varindex][i]]);
         }

	 chgcons0 = FALSE;

	 /* check for overlapping constraint after lifting, in the first round we will only check up front */
	 SCIP_CALL( checkForOverlapping(scip, cons0, c, (conshdlrdata->nclqpresolve > 0) ? nusefulconss : c,
	       usefulconss, nusefulconss, usefulvars, &nusefulvars, vartoindex, varnconss, maxnvarconsidx, varconsidxs,
	       countofoverlapping, conshdlrdata->cliqueshrinking, &chgcons0,
	       undoneaggrvars, undoneaggrtypes, &naggregations, &saggregations,
	       nfixedvars, naggrvars, nchgcoefs, ndelconss, cutoff) );

	 if( *cutoff )
	    break;

	 /* we can't handle aggregated variables later on so we should have saved them for later */
	 assert(*naggrvars == oldnaggrvars);
      }
   }

 TERMINATE:
   if( useadditionalmem )
   {
      SCIPfreeBufferArray(scip, &countofoverlapping);
      SCIPfreeBufferArray(scip, &cliquevalues);
   }

   /* free temporary memory for constraints, variables and the mapping between them in reverse order as they were
    * allocated
    */
   for( c = nusefulvars; c > 0; --c )
   {
      if( varconsidxs[c] != NULL )
      {
         SCIPfreeBufferArrayNull(scip, &(varconsidxs[c]));
      }
   }

   SCIPfreeBufferArray(scip, &varconsidxs);
   SCIPfreeBufferArray(scip, &maxnvarconsidx);
   SCIPfreeBufferArray(scip, &varnconss);
   SCIPfreeBufferArray(scip, &usefulvars);
   SCIPfreeBufferArray(scip, &usefulconss);


   /* perform all collected aggregations */
   if( !*cutoff && naggregations > 0 && !SCIPdoNotAggr(scip) )
   {
      SCIP_CALL( performAggregations(scip, conshdlrdata, undoneaggrvars, undoneaggrtypes, naggregations, naggrvars, cutoff) );
   }

   /* free temporary memory for the aggregation storage */
   SCIPfreeBufferArray(scip, &undoneaggrtypes);
   SCIPfreeBufferArray(scip, &undoneaggrvars);

   /* free hashmap */
   SCIPhashmapFree(&vartoindex);

   if( *ndelconss < 0 )
      *ndelconss = 0;

   return SCIP_OKAY;
}

/** compares each constraint with all other constraints for possible redundancy and removes or changes constraint
 *  accordingly; in contrast to removeRedundantConstraints(), it uses a hash table
 */
static
SCIP_RETCODE detectRedundantConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   nconss,             /**< number of constraints in constraint set */
   int*                  firstchange,        /**< pointer to store first changed constraint */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints */
   int*                  nchgsides           /**< pointer to count number of changed left/right hand sides */
)
{
   SCIP_HASHTABLE* hashtable;
   int hashtablesize;
   int c;

   assert(scip != NULL);
   assert(blkmem != NULL);
   assert(conss != NULL || nconss == 0);
   assert(firstchange != NULL);
   assert(ndelconss != NULL);
   assert(nchgsides != NULL);

   if( nconss == 0 )
      return SCIP_OKAY;

   assert(conss != NULL);

   /* create a hash table for the constraint set */
   hashtablesize = SCIPcalcHashtableSize(10*nconss);
   hashtablesize = MAX(hashtablesize, HASHSIZE_SETPPCCONS);
   SCIP_CALL( SCIPhashtableCreate(&hashtable, blkmem, hashtablesize,
         hashGetKeySetppccons, hashKeyEqSetppccons, hashKeyValSetppccons, (void*) scip) );

   /* check all constraints in the given set for redundancy */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONS* cons0;
      SCIP_CONS* cons1;

      cons0 = conss[c];

      if( !SCIPconsIsActive(cons0) || SCIPconsIsModifiable(cons0) )
         continue;

      /* get constraint from current hash table with same variables as cons0 and with coefficients either equal or negated
       * to the ones of cons0 */
      cons1 = (SCIP_CONS*)(SCIPhashtableRetrieve(hashtable, (void*)cons0));

      if( cons1 != NULL )
      {
         SCIP_CONSDATA* consdata0;
         SCIP_CONSDATA* consdata1;

         assert(SCIPconsIsActive(cons1));
         assert(!SCIPconsIsModifiable(cons1));

         /* constraint found: create a new constraint with same coefficients and best left and right hand side;
          * delete old constraints afterwards
          */
         consdata0 = SCIPconsGetData(cons0);
         consdata1 = SCIPconsGetData(cons1);

         assert(consdata0 != NULL && consdata1 != NULL);
         assert(consdata0->nvars >= 1 && consdata0->nvars == consdata1->nvars);

         assert(consdata0->sorted && consdata1->sorted);
         assert(consdata0->vars[0] == consdata1->vars[0]);

         SCIPdebugMessage("setppc constraints <%s> and <%s> have identical variable sets\n",
            SCIPconsGetName(cons0), SCIPconsGetName(cons1));
         SCIPdebugPrintCons(scip, cons0, NULL);
         SCIPdebugPrintCons(scip, cons1, NULL);

         /* if necessary change type of setppc constraint */
         if( consdata1->setppctype != SCIP_SETPPCTYPE_PARTITIONING && consdata0->setppctype != consdata1->setppctype ) /*lint !e641*/
         {
            /* change the type of cons0 */
            SCIP_CALL( setSetppcType(scip, cons1, SCIP_SETPPCTYPE_PARTITIONING) );
            (*nchgsides)++;
         }

         /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
         SCIP_CALL( updateFlags(scip, cons1, cons0) );

         /* delete cons0 */
         SCIP_CALL( SCIPdelCons(scip, cons0) );
         (*ndelconss)++;

         /* update the first changed constraint to begin the next aggregation round with */
         if( consdata0->changed && SCIPconsGetPos(cons1) < *firstchange )
            *firstchange = SCIPconsGetPos(cons1);

         assert(SCIPconsIsActive(cons1));
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

/** removes the redundant second constraint and updates the flags of the first one */
static
SCIP_RETCODE removeRedundantCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons0,              /**< constraint that should stay */
   SCIP_CONS*            cons1,              /**< constraint that should be deleted */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   assert(ndelconss != NULL);

   SCIPdebugMessage(" -> removing setppc constraint <%s> which is redundant to <%s>\n",
      SCIPconsGetName(cons1), SCIPconsGetName(cons0));
   SCIPdebugPrintCons(scip, cons0, NULL);
   SCIPdebugPrintCons(scip, cons1, NULL);

   /* update flags of cons0 */
   SCIP_CALL( updateFlags(scip, cons0, cons1) );

   /* delete cons1 */
   SCIP_CALL( SCIPdelCons(scip, cons1) );
   (*ndelconss)++;

   return SCIP_OKAY;
}

/** for cons0 contained in cons1, fixes variables of cons1 that are not in cons0 to zero */
static
SCIP_RETCODE fixAdditionalVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons0,              /**< constraint that is contained in the other */
   SCIP_CONS*            cons1,              /**< constraint that is a superset of the other */
   SCIP_Bool*            cutoff,             /**< pointer to store whether a cutoff was found */
   int*                  nfixedvars          /**< pointer to count number of fixed variables */
   )
{
   SCIP_CONSDATA* consdata0;
   SCIP_CONSDATA* consdata1;
   int v0;
   int v1;

   assert(cutoff != NULL);
   assert(nfixedvars != NULL);

   *cutoff = FALSE;

   /* get constraint datas */
   consdata0 = SCIPconsGetData(cons0);
   consdata1 = SCIPconsGetData(cons1);
   assert(consdata0 != NULL);
   assert(consdata1 != NULL);
   assert(consdata0->nvars < consdata1->nvars);
   assert(consdata0->sorted);
   assert(consdata1->sorted);

   /* fix variables in the range of cons0 */
   for( v0 = 0, v1 = 0; v0 < consdata0->nvars && !(*cutoff); ++v0, ++v1 )
   {
      int index0;

      assert(v1 < consdata1->nvars);
      index0 = SCIPvarGetIndex(consdata0->vars[v0]);
      for( ; SCIPvarGetIndex(consdata1->vars[v1]) < index0 && !(*cutoff); ++v1 )
      {
         SCIP_Bool fixed;

         /* fix variable to zero */
         SCIP_CALL( SCIPfixVar(scip, consdata1->vars[v1], 0.0, cutoff, &fixed) );
         if( fixed )
         {
            SCIPdebugMessage(" -> fixed <%s> == 0\n", SCIPvarGetName(consdata1->vars[v1]));
            (*nfixedvars)++;
         }
         assert(v1 < consdata1->nvars-1);
      }
      assert(SCIPvarGetIndex(consdata1->vars[v1]) == index0 || *cutoff);
   }

   /* fix remaining variables of cons1 */
   for( ; v1 < consdata1->nvars && !(*cutoff); ++v1 )
   {
      SCIP_Bool fixed;

      assert(consdata0->nvars == 0
         || SCIPvarGetIndex(consdata1->vars[v1]) > SCIPvarGetIndex(consdata0->vars[consdata0->nvars-1]));

      /* fix variable to zero */
      SCIP_CALL( SCIPfixVar(scip, consdata1->vars[v1], 0.0, cutoff, &fixed) );
      if( fixed )
      {
         SCIPdebugMessage(" -> fixed <%s> == 0\n", SCIPvarGetName(consdata1->vars[v1]));
         (*nfixedvars)++;
      }
   }

   return SCIP_OKAY;
}

/** applies reductions for cons0 being contained in cons1 */
static
SCIP_RETCODE processContainedCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons0,              /**< constraint that is contained in the other */
   SCIP_CONS*            cons1,              /**< constraint that is a superset of the other */
   SCIP_Bool*            cutoff,             /**< pointer to store whether a cutoff was found */
   int*                  nfixedvars,         /**< pointer to count number of fixed variables */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints */
   int*                  nchgsides           /**< pointer to count number of changed left/right hand sides */
   )
{
   SCIP_CONSDATA* consdata0;
   SCIP_CONSDATA* consdata1;

   assert(cutoff != NULL);
   assert(nfixedvars != NULL);
   assert(ndelconss != NULL);
   assert(nchgsides != NULL);

   *cutoff = FALSE;

   /* get constraint datas */
   consdata0 = SCIPconsGetData(cons0);
   consdata1 = SCIPconsGetData(cons1);
   assert(consdata0 != NULL);
   assert(consdata1 != NULL);
   assert(consdata0->nvars < consdata1->nvars);
   assert(consdata0->sorted);
   assert(consdata1->sorted);

   switch( consdata0->setppctype )
   {
   case SCIP_SETPPCTYPE_PARTITIONING:
      switch( consdata1->setppctype )
      {
      case SCIP_SETPPCTYPE_PARTITIONING:
      case SCIP_SETPPCTYPE_PACKING:
         /* cons0: partitioning, cons1: partitioning or packing
          * -> fix additional variables in cons1 to zero, remove cons1
          */
         SCIP_CALL( fixAdditionalVars(scip, cons0, cons1, cutoff, nfixedvars) );
         SCIP_CALL( removeRedundantCons(scip, cons0, cons1, ndelconss) );
         break;

      case SCIP_SETPPCTYPE_COVERING:
         /* cons0: partitioning, cons1: covering
          * -> remove cons1
          */
         SCIP_CALL( removeRedundantCons(scip, cons0, cons1, ndelconss) );
         break;

      default:
         SCIPerrorMessage("invalid setppc type <%d> of constraint <%s>\n", consdata1->setppctype, SCIPconsGetName(cons1));
         return SCIP_INVALIDDATA;
      }
      break;

   case SCIP_SETPPCTYPE_PACKING:
      switch( consdata1->setppctype )
      {
      case SCIP_SETPPCTYPE_PARTITIONING:
      case SCIP_SETPPCTYPE_PACKING:
         /* cons0: packing, cons1: partitioning or packing
          * -> remove cons0
          */
         SCIP_CALL( removeRedundantCons(scip, cons1, cons0, ndelconss) );
         break;

      case SCIP_SETPPCTYPE_COVERING:
         /* cons0: packing, cons1: covering
          * -> nothing can be deduced
          */
         break;

      default:
         SCIPerrorMessage("invalid setppc type <%d> of constraint <%s>\n", consdata1->setppctype, SCIPconsGetName(cons1));
         return SCIP_INVALIDDATA;
      }
      break;

   case SCIP_SETPPCTYPE_COVERING:
      switch( consdata1->setppctype )
      {
      case SCIP_SETPPCTYPE_PARTITIONING:
      case SCIP_SETPPCTYPE_PACKING:
         /* cons0: covering, cons1: partitioning or packing
          * -> fix additional variables in cons1 to zero, remove cons1, convert cons0 into partitioning
          */
         SCIP_CALL( fixAdditionalVars(scip, cons0, cons1, cutoff, nfixedvars) );
         SCIP_CALL( setSetppcType(scip, cons0, SCIP_SETPPCTYPE_PARTITIONING) );
         SCIP_CALL( removeRedundantCons(scip, cons0, cons1, ndelconss) );
         (*nchgsides)++;
         break;

      case SCIP_SETPPCTYPE_COVERING:
         /* cons0: covering, cons1: covering
          * -> remove cons1
          */
         SCIP_CALL( removeRedundantCons(scip, cons0, cons1, ndelconss) );
         break;

      default:
         SCIPerrorMessage("invalid setppc type <%d> of constraint <%s>\n", consdata1->setppctype, SCIPconsGetName(cons1));
         return SCIP_INVALIDDATA;
      }
      break;

   default:
      SCIPerrorMessage("invalid setppc type <%d> of constraint <%s>\n", consdata0->setppctype, SCIPconsGetName(cons0));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** deletes redundant constraints */
static
SCIP_RETCODE removeRedundantConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   firstchange,        /**< first constraint that changed since last pair preprocessing round */
   int                   chkind,             /**< index of constraint to check against all prior indices upto startind */
   SCIP_Bool*            cutoff,             /**< pointer to store whether a cutoff was found */
   int*                  nfixedvars,         /**< pointer to count number of fixed variables */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints */
   int*                  nchgsides           /**< pointer to count number of changed left/right hand sides */
   )
{
   SCIP_CONS* cons0;
   SCIP_CONSDATA* consdata0;
   SCIP_Longint signature0;
   SCIP_Bool cons0changed;
   int c;

   assert(scip != NULL);
   assert(conss != NULL);
   assert(cutoff != NULL);
   assert(nfixedvars != NULL);
   assert(ndelconss != NULL);
   assert(nchgsides != NULL);

   *cutoff = FALSE;

   /* get the constraint to be checked against all prior constraints */
   cons0 = conss[chkind];
   assert(SCIPconsIsActive(cons0));
   assert(!SCIPconsIsModifiable(cons0));

   consdata0 = SCIPconsGetData(cons0);
   assert(consdata0 != NULL);
   assert(consdata0->nvars >= 1);

   /* sort the constraint cons0 */
   consdataSort(consdata0);

   /* get the bit signature of the constraint */
   signature0 = consdataGetSignature(consdata0);

   /* check constraint against all prior constraints */
   cons0changed = consdata0->changed;
   consdata0->changed = FALSE;
   for( c = (cons0changed ? 0 : firstchange); c < chkind && !(*cutoff) && SCIPconsIsActive(cons0); ++c )
   {
      SCIP_CONS* cons1;
      SCIP_CONSDATA* consdata1;
      SCIP_Longint signature1;
      SCIP_Longint jointsignature;
      SCIP_Bool cons0iscontained;
      SCIP_Bool cons1iscontained;
      int v0;
      int v1;

      cons1 = conss[c];

      /* ignore inactive and modifiable constraints */
      if( !SCIPconsIsActive(cons1) || SCIPconsIsModifiable(cons1) )
         continue;

      consdata1 = SCIPconsGetData(cons1);
      assert(consdata1 != NULL);

      /* sort the constraint cons1 */
      consdataSort(consdata1);

      /* get the bit signature of cons1 */
      signature1 = consdataGetSignature(consdata1);

      /* check (based on signature) if the two constraints are not included in each other */
      jointsignature = (signature0 | signature1);
      if( jointsignature != signature0 && jointsignature != signature1 )
         continue;

      /* check whether one constraint is really a subset of the other */
      cons0iscontained = (consdata0->nvars <= consdata1->nvars);
      cons1iscontained = (consdata1->nvars <= consdata0->nvars);
      v0 = 0;
      v1 = 0;
      while( v0 < consdata0->nvars && v1 < consdata1->nvars )
      {
         int index0;
         int index1;

         index0 = SCIPvarGetIndex(consdata0->vars[v0]);
         index1 = SCIPvarGetIndex(consdata1->vars[v1]);
         if( index0 < index1 )
         {
            cons0iscontained = FALSE;
            if( !cons1iscontained )
               break;
            for( v0++; v0 < consdata0->nvars && SCIPvarGetIndex(consdata0->vars[v0]) < index1; v0++ )
            {}
         }
         else if( index1 < index0 )
         {
            cons1iscontained = FALSE;
            if( !cons0iscontained )
               break;
            for( v1++; v1 < consdata1->nvars && SCIPvarGetIndex(consdata1->vars[v1]) < index0; v1++ )
            {}
         }
         else
         {
            v0++;
            v1++;
         }
      }
      cons0iscontained = cons0iscontained && (v0 == consdata0->nvars);
      cons1iscontained = cons1iscontained && (v1 == consdata1->nvars);

      if( cons0iscontained && cons1iscontained )
      {
         SCIPdebugMessage("setppc constraints <%s> and <%s> have identical variable sets\n",
            SCIPconsGetName(cons0), SCIPconsGetName(cons1));
         SCIPdebugPrintCons(scip, cons0, NULL);
         SCIPdebugPrintCons(scip, cons1, NULL);

         /* both constraints consists of the same variables */
         if( consdata0->setppctype == consdata1->setppctype )
         {
            /* both constraints are equal: update flags in cons0 and delete cons1 */
            SCIP_CALL( removeRedundantCons(scip, cons0, cons1, ndelconss) );
         }
         else if( consdata0->setppctype == SCIP_SETPPCTYPE_PARTITIONING ) /*lint !e641*/
         {
            /* the set partitioning constraint is stronger: remove the other one */
            SCIP_CALL( removeRedundantCons(scip, cons0, cons1, ndelconss) );
         }
         else if( consdata1->setppctype == SCIP_SETPPCTYPE_PARTITIONING ) /*lint !e641*/
         {
            /* the set partitioning constraint is stronger: remove the other one */
            SCIP_CALL( removeRedundantCons(scip, cons1, cons0, ndelconss) );
         }
         else
         {
            /* one is a covering, the other one a packing constraint: replace them by a single partitioning constraint */
            assert((consdata0->setppctype == SCIP_SETPPCTYPE_COVERING && consdata1->setppctype == SCIP_SETPPCTYPE_PACKING)
               || (consdata1->setppctype == SCIP_SETPPCTYPE_COVERING && consdata0->setppctype == SCIP_SETPPCTYPE_PACKING)); /*lint !e641*/

            /* change the type of cons0 */
            SCIP_CALL( setSetppcType(scip, cons0, SCIP_SETPPCTYPE_PARTITIONING) );
            (*nchgsides)++;

            /* delete cons1 */
            SCIP_CALL( removeRedundantCons(scip, cons0, cons1, ndelconss) );
         }
      }
      else if( cons0iscontained )
      {
         /* cons0 is contained in cons1 */
         SCIPdebugMessage("setppc constraint <%s> is contained in <%s>\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));
         SCIPdebugPrintCons(scip, cons0, NULL);
         SCIPdebugPrintCons(scip, cons1, NULL);
         SCIP_CALL( processContainedCons(scip, cons0, cons1, cutoff, nfixedvars, ndelconss, nchgsides) );
      }
      else if( cons1iscontained )
      {
         /* cons1 is contained in cons1 */
         SCIPdebugMessage("setppc constraint <%s> is contained in <%s>\n", SCIPconsGetName(cons1), SCIPconsGetName(cons0));
         SCIPdebugPrintCons(scip, cons0, NULL);
         SCIPdebugPrintCons(scip, cons1, NULL);
         SCIP_CALL( processContainedCons(scip, cons1, cons0, cutoff, nfixedvars, ndelconss, nchgsides) );
      }
   }

   return SCIP_OKAY;
}

/* perform deletion of variables in all constraints of the constraint handler */
static
SCIP_RETCODE performVarDeletions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss              /**< number of constraints */
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   int v;


   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(nconss >= 0);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* iterate over all constraints */
   for( i = 0; i < nconss; i++ )
   {
      consdata = SCIPconsGetData(conss[i]);

      /* constraint is marked, that some of its variables were deleted */
      if( consdata->varsdeleted )
      {
         /* iterate over all variables of the constraint and delete marked variables */
         for( v = consdata->nvars - 1; v >= 0; v-- )
         {
            if( SCIPvarIsDeleted(consdata->vars[v]) )
            {
               SCIP_CALL( delCoefPos(scip, conss[i], v) );
            }
         }
         consdata->varsdeleted = FALSE;
      }
   }

   return SCIP_OKAY;
}


/*
 * upgrading of linear constraints
 */


/** creates and captures a set partitioning / packing / covering constraint */
static
SCIP_RETCODE createConsSetppc(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_SETPPCTYPE       setppctype,         /**< type of constraint: set partitioning, packing, or covering constraint */
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

   assert(scip != NULL);

   /* find the set partitioning constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("set partitioning / packing / covering constraint handler not found\n");
      return SCIP_INVALIDCALL;
   }

   /* create the constraint specific data */
   if( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM )
   {
      /* create constraint in original problem */
      SCIP_CALL( consdataCreate(scip, &consdata, nvars, vars, setppctype) );
   }
   else
   {
      /* create constraint in transformed problem */
      SCIP_CALL( consdataCreateTransformed(scip, &consdata, nvars, vars, setppctype) );
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   if( SCIPgetStage(scip) != SCIP_STAGE_PROBLEM )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* catch bound change events of variables */
      SCIP_CALL( catchAllEvents(scip, *cons, conshdlrdata->eventhdlr) );
   }

   return SCIP_OKAY;
}

/** creates and captures a normalized (with all coefficients +1) setppc constraint */
static
SCIP_RETCODE createNormalizedSetppc(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            vals,               /**< array with coefficients (+1.0 or -1.0) */
   int                   mult,               /**< multiplier on the coefficients(+1 or -1) */
   SCIP_SETPPCTYPE       setppctype,         /**< type of constraint: set partitioning, packing, or covering constraint */
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
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);
   assert(mult == +1 || mult == -1);

   /* get temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &transvars, nvars) );

   /* negate positive or negative variables */
   for( v = 0; v < nvars; ++v )
   {
      if( mult * vals[v] > 0.0 )
         transvars[v] = vars[v];
      else
      {
         SCIP_CALL( SCIPgetNegatedVar(scip, vars[v], &transvars[v]) );
      }
      assert(transvars[v] != NULL);
   }

   /* create the constraint */
   SCIP_CALL( createConsSetppc(scip, cons, name, nvars, transvars, setppctype,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

   /* release temporary memory */
   SCIPfreeBufferArray(scip, &transvars);

   return SCIP_OKAY;
}

static
SCIP_DECL_LINCONSUPGD(linconsUpgdSetppc)
{  /*lint --e{715}*/
   assert(upgdcons != NULL);

   /* check, if linear constraint can be upgraded to set partitioning, packing, or covering constraint
    * - all set partitioning / packing / covering constraints consist only of binary variables with a
    *   coefficient of +1.0 or -1.0 (variables with -1.0 coefficients can be negated):
    *        lhs     <= x1 + ... + xp - y1 - ... - yn <= rhs
    * - negating all variables y = (1-Y) with negative coefficients gives:
    *        lhs + n <= x1 + ... + xp + Y1 + ... + Yn <= rhs + n
    * - negating all variables x = (1-X) with positive coefficients and multiplying with -1 gives:
    *        p - rhs <= X1 + ... + Xp + y1 + ... + yn <= p - lhs
    * - a set partitioning constraint has left hand side of +1.0, and right hand side of +1.0 : x(S) == 1.0
    *    -> without negations:  lhs == rhs == 1 - n  or  lhs == rhs == p - 1
    * - a set packing constraint has left hand side of -infinity, and right hand side of +1.0 : x(S) <= 1.0
    *    -> without negations:  (lhs == -inf  and  rhs == 1 - n)  or  (lhs == p - 1  and  rhs = +inf)
    * - a set covering constraint has left hand side of +1.0, and right hand side of +infinity: x(S) >= 1.0
    *    -> without negations:  (lhs == 1 - n  and  rhs == +inf)  or  (lhs == -inf  and  rhs = p - 1)
    */
   if( nposbin + nnegbin == nvars && ncoeffspone + ncoeffsnone == nvars )
   {
      int mult;

      if( SCIPisEQ(scip, lhs, rhs) && (SCIPisEQ(scip, lhs, 1.0 - ncoeffsnone) || SCIPisEQ(scip, lhs, ncoeffspone - 1.0)) )
      {
         SCIPdebugMessage("upgrading constraint <%s> to set partitioning constraint\n", SCIPconsGetName(cons));

         /* check, if we have to multiply with -1 (negate the positive vars) or with +1 (negate the negative vars) */
         mult = SCIPisEQ(scip, lhs, 1.0 - ncoeffsnone) ? +1 : -1;

         /* create the set partitioning constraint (an automatically upgraded constraint is always unmodifiable) */
         assert(!SCIPconsIsModifiable(cons));
         SCIP_CALL( createNormalizedSetppc(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, mult,
               SCIP_SETPPCTYPE_PARTITIONING,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
               SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
               SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
      }
      else if( (SCIPisInfinity(scip, -lhs) && SCIPisEQ(scip, rhs, 1.0 - ncoeffsnone))
         || (SCIPisEQ(scip, lhs, ncoeffspone - 1.0) && SCIPisInfinity(scip, rhs)) )
      {
         SCIPdebugMessage("upgrading constraint <%s> to set packing constraint\n", SCIPconsGetName(cons));

         /* check, if we have to multiply with -1 (negate the positive vars) or with +1 (negate the negative vars) */
         mult = SCIPisInfinity(scip, -lhs) ? +1 : -1;

         /* create the set packing constraint (an automatically upgraded constraint is always unmodifiable) */
         assert(!SCIPconsIsModifiable(cons));
         SCIP_CALL( createNormalizedSetppc(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, mult,
               SCIP_SETPPCTYPE_PACKING,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
               SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
               SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
      }
      else if( (SCIPisEQ(scip, lhs, 1.0 - ncoeffsnone) && SCIPisInfinity(scip, rhs))
         || (SCIPisInfinity(scip, -lhs) && SCIPisEQ(scip, rhs, ncoeffspone - 1.0)) )
      {
         SCIPdebugMessage("upgrading constraint <%s> to set covering constraint\n", SCIPconsGetName(cons));

         /* check, if we have to multiply with -1 (negate the positive vars) or with +1 (negate the negative vars) */
         mult = SCIPisInfinity(scip, rhs) ? +1 : -1;

         /* create the set covering constraint (an automatically upgraded constraint is always unmodifiable) */
         assert(!SCIPconsIsModifiable(cons));
         SCIP_CALL( createNormalizedSetppc(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, mult,
               SCIP_SETPPCTYPE_COVERING,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
               SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
               SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopySetppc)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrSetppc(scip) );
 
   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeSetppc)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( conshdlrdataFree(scip, &conshdlrdata) );

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitSetppc)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->noldfixedvars = 0;
   conshdlrdata->noldimpls = 0;
   conshdlrdata->noldcliques = 0;
   conshdlrdata->noldupgrs = 0;
   conshdlrdata->nclqpresolve = 0;
   conshdlrdata->updatedsetppctype = FALSE;
   conshdlrdata->enablecliquelifting = TRUE;
   conshdlrdata->randseed = 1010010007; /*@todo try to take a better/another value */

   return SCIP_OKAY;
}

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolSetppc)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

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
SCIP_DECL_CONSDELETE(consDeleteSetppc)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* if constraint belongs to transformed problem space, drop bound change events on variables */
   if( (*consdata)->nvars > 0 && SCIPvarIsTransformed((*consdata)->vars[0]) )
   {
      SCIP_CALL( dropAllEvents(scip, cons, conshdlrdata->eventhdlr) );
   }

   /* free setppc constraint data */
   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransSetppc)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   /*debugMessage("Trans method of setppc constraints\n");*/

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->row == NULL);  /* in original problem, there cannot be LP rows */

   /* create constraint data for target constraint */
   SCIP_CALL( consdataCreateTransformed(scip, &targetdata, sourcedata->nvars, sourcedata->vars,
         (SCIP_SETPPCTYPE)sourcedata->setppctype) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   /* catch bound change events of variables */
   SCIP_CALL( catchAllEvents(scip, *targetcons, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpSetppc)
{  /*lint --e{715}*/
   int c;

   for( c = 0; c < nconss; ++c )
   {
      assert(SCIPconsIsInitial(conss[c]));
      SCIP_CALL( addCut(scip, conss[c], NULL) );
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpSetppc)
{  /*lint --e{715}*/
   SCIP_Bool cutoff;
   SCIP_Bool separated;
   SCIP_Bool reduceddom;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   SCIPdebugMessage("separating %d/%d set partitioning / packing / covering constraints\n", nusefulconss, nconss);

   *result = SCIP_DIDNOTFIND;

   cutoff = FALSE;
   separated = FALSE;
   reduceddom = FALSE;

   /* check all useful set partitioning / packing / covering constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], NULL, &cutoff, &separated, &reduceddom) );
   }

   /* combine set partitioning / packing / covering constraints to get more cuts */
   /**@todo further cuts of set partitioning / packing / covering constraints */

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;
   else if( separated )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolSetppc)
{  /*lint --e{715}*/
   SCIP_Bool cutoff;
   SCIP_Bool separated;
   SCIP_Bool reduceddom;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   SCIPdebugMessage("separating %d/%d set partitioning / packing / covering constraints\n", nusefulconss, nconss);

   *result = SCIP_DIDNOTFIND;

   cutoff = FALSE;
   separated = FALSE;
   reduceddom = FALSE;

   /* check all useful set partitioning / packing / covering constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], sol, &cutoff, &separated, &reduceddom) );
   }

   /* combine set partitioning / packing / covering constraints to get more cuts */
   /**@todo further cuts of set partitioning / packing / covering constraints */

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;
   else if( separated )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


#ifdef VARUSES
#ifdef BRANCHLP
/** if fractional variables exist, chooses a set S of them and branches on (i) x(S) == 0, and (ii) x(S) >= 1 */
static
SCIP_RETCODE branchLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< set partitioning / packing / covering constraint handler */
   SCIP_RESULT*          result              /**< pointer to store the result SCIP_BRANCHED, if branching was applied */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_INTARRAY* varuses;
   SCIP_VAR** lpcands;
   SCIP_VAR** sortcands;
   SCIP_VAR* var;
   SCIP_Real branchweight;
   SCIP_Real solval;
   int* uses;
   int nlpcands;
   int nsortcands;
   int nselcands;
   int numuses;
   int i;
   int j;

   /**@todo use a better set partitioning / packing / covering branching on LP solution (use SOS branching) */

   assert(conshdlr != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   varuses = conshdlrdata->varuses;
   assert(varuses != NULL);

   /* get fractional variables */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, NULL, NULL, &nlpcands, NULL) );
   if( nlpcands == 0 )
      return SCIP_OKAY;

   assert(MINBRANCHWEIGHT <= MAXBRANCHWEIGHT);

   /* get temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &sortcands, nlpcands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &uses, nlpcands) );

   /* sort fractional variables by number of uses in enabled set partitioning / packing / covering constraints */
   nsortcands = 0;
   for( i = 0; i < nlpcands; ++i )
   {
      var = lpcands[i];
      numuses = SCIPgetIntarrayVal(scip, varuses, SCIPvarGetIndex(var));
      if( numuses > 0 )
      {
         for( j = nsortcands; j > 0 && numuses > uses[j-1]; --j )
         {
            sortcands[j] = sortcands[j-1];
            uses[j] = uses[j-1];
         }
         assert(0 <= j && j <= nsortcands);
         sortcands[j] = var;
         uses[j] = numuses;
         nsortcands++;
      }
   }
   assert(nsortcands <= nlpcands);

   /* if none of the fractional variables is member of a set partitioning / packing / covering constraint,
    * we are not responsible for doing the branching
    */
   if( nsortcands > 0 )
   {
      SCIP_Real cumprio = 0.0;
      SCIP_Real minprio = SCIP_INVALID;
      SCIP_Real minestzero = SCIP_INVALID;
      SCIP_Real minestone = SCIP_INVALID;
      SCIP_Real tmp;

      /* select the first variables from the sorted candidate list, until MAXBRANCHWEIGHT is reached;
       * then choose one less
       */
      branchweight = 0.0;
      solval = 0.0;
      for( nselcands = 0; nselcands < nsortcands; ++nselcands )
      {
         solval = SCIPgetVarSol(scip, sortcands[nselcands]);
         assert(SCIPisFeasGE(scip, solval, 0.0) && SCIPisFeasLE(scip, solval, 1.0));
         branchweight += solval;

	 /* did we exceed the maximal weight */
	 if( branchweight > MAXBRANCHWEIGHT )
	    break;

	 /* @todo instead of taking the minimum into account try other variants like the maximum and the average */
	 /* calculate priorities and estimates by adding up/taking the minimum of all single priorities/estimates */
	 cumprio += SCIPcalcNodeselPriority(scip, sortcands[nselcands], SCIP_BRANCHDIR_DOWNWARDS, 0.0);
	 tmp = SCIPcalcNodeselPriority(scip, sortcands[nselcands], SCIP_BRANCHDIR_UPWARDS, 1.0);
	 minprio = MIN(minprio, tmp);
	 tmp = SCIPcalcChildEstimate(scip, sortcands[nselcands], 0.0);;
	 minestzero = MIN(minestzero, tmp);
	 tmp = SCIPcalcChildEstimate(scip, sortcands[nselcands], 1.0);;
	 minestone = MIN(minestone, tmp);
      }
      assert(minestzero != SCIP_INVALID);
      assert(minestone != SCIP_INVALID);
      assert(minprio != SCIP_INVALID);
      assert(nselcands > 0);
      branchweight -= solval;

      /* check, if we accumulated at least MIN and at most MAXBRANCHWEIGHT weight */
      if( MINBRANCHWEIGHT <= branchweight && branchweight <= MAXBRANCHWEIGHT )
      {
         SCIP_NODE* node;

         /* perform the binary set branching on the selected variables */
         assert(1 <= nselcands && nselcands <= nlpcands);

         /* create left child, fix x_i = 0 for all i \in S */
         SCIP_CALL( SCIPcreateChild(scip, &node, cumprio, minestzero) );
         for( i = 0; i < nselcands; ++i )
         {
            SCIP_CALL( SCIPchgVarUbNode(scip, node, sortcands[i], 0.0) );
         }

         /* create right child: add constraint x(S) >= 1 */
         SCIP_CALL( SCIPcreateChild(scip, &node, minprio, minestone) );
         if( nselcands == 1 )
         {
            /* only one candidate selected: fix it to 1.0 */
            SCIPdebugMessage("fixing variable <%s> to 1.0 in right child node\n", SCIPvarGetName(sortcands[0]));
            SCIP_CALL( SCIPchgVarLbNode(scip, node, sortcands[0], 1.0) );
         }
         else
         {
            SCIP_CONS* newcons;
            char name[SCIP_MAXSTRLEN];

            /* add set covering constraint x(S) >= 1 */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "BSB%"SCIP_LONGINT_FORMAT, SCIPgetNTotalNodes(scip));

            SCIP_CALL( SCIPcreateConsSetcover(scip, &newcons, name, nselcands, sortcands,
                  FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE) );
            SCIP_CALL( SCIPaddConsNode(scip, node, newcons, NULL) );
            SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
         }

         *result = SCIP_BRANCHED;

#ifdef SCIP_DEBUG
         SCIPdebugMessage("binary set branching: nselcands=%d/%d, weight(S)=%g, A={", nselcands, nlpcands, branchweight);
         for( i = 0; i < nselcands; ++i )
            SCIPdebugPrintf(" %s[%g]", SCIPvarGetName(sortcands[i]), SCIPgetSolVal(scip, NULL, sortcands[i]));
         SCIPdebugPrintf(" }\n");
#endif
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &uses);
   SCIPfreeBufferArray(scip, &sortcands);

   return SCIP_OKAY;
}
#endif

/** if unfixed variables exist, chooses a set S of them and creates |S|+1 child nodes:
 *   - for each variable i from S, create child node with x_0 = ... = x_i-1 = 0, x_i = 1
 *   - create an additional child node x_0 = ... = x_n-1 = 0
 */
static
SCIP_RETCODE branchPseudo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< set partitioning / packing / covering constraint handler */
   SCIP_RESULT*          result              /**< pointer to store the result SCIP_BRANCHED, if branching was applied */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_INTARRAY* varuses;
   SCIP_VAR** pseudocands;
   SCIP_VAR** branchcands;
   SCIP_VAR* var;
   SCIP_NODE* node;
   int* canduses;
   int npseudocands;
   int maxnbranchcands;
   int nbranchcands;
   int uses;
   int i;
   int j;

   /**@todo use a better set partitioning / packing / covering branching on pseudo solution (use SOS branching) */

   assert(conshdlr != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* check, if pseudo branching is disabled */
   if( conshdlrdata->npseudobranches <= 1 )
      return SCIP_OKAY;

   /* get fractional variables */
   SCIP_CALL( SCIPgetPseudoBranchCands(scip, &pseudocands, NULL, &npseudocands) );
   if( npseudocands == 0 )
      return SCIP_OKAY;

   varuses = conshdlrdata->varuses;
   assert(varuses != NULL);

   /* choose the maximal number of branching variables */
   maxnbranchcands = conshdlrdata->npseudobranches-1;
   assert(maxnbranchcands >= 1);

   /* get temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &branchcands, maxnbranchcands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &canduses, maxnbranchcands) );

   /* sort unfixed variables by number of uses in enabled set partitioning / packing / covering constraints */
   nbranchcands = 0;
   for( i = 0; i < npseudocands; ++i )
   {
      var = pseudocands[i];
      uses = SCIPgetIntarrayVal(scip, varuses, SCIPvarGetIndex(var));
      if( uses > 0 )
      {
         if( nbranchcands < maxnbranchcands || uses > canduses[nbranchcands-1] )
         {
            for( j = MIN(nbranchcands, maxnbranchcands-1); j > 0 && uses > canduses[j-1]; --j )
            {
               branchcands[j] = branchcands[j-1];
               canduses[j] = canduses[j-1];
            }
            assert(0 <= j && j <= nbranchcands && j < maxnbranchcands);
            branchcands[j] = var;
            canduses[j] = uses;
            if( nbranchcands < maxnbranchcands )
               nbranchcands++;
         }
      }
   }
   assert(nbranchcands <= maxnbranchcands);

   /* if none of the unfixed variables is member of a set partitioning / packing / covering constraint,
    * we are not responsible for doing the branching
    */
   if( nbranchcands > 0 )
   {
      SCIP_Real* estone;
      SCIP_Real minestzero = SCIP_INVALID;
      SCIP_Real tmp;

      SCIP_CALL( SCIPallocBufferArray(scip, &estone, nbranchcands) );

      /* @todo instead of taking the minimum into account try other variants like the maximum and the average */
      /* @todo calculate priorities instead of setting it to the number of branching candidates */
      /* calculate estimates by taking the minimum over all single estimates */
      for( i = 0; i < nbranchcands; ++i )
      {
	 tmp = SCIPcalcChildEstimate(scip, branchcands[i], 0.0);;
	 minestzero = MIN(minestzero, tmp);
	 estone[i] = SCIPcalcChildEstimate(scip, branchcands[i], 1.0);
      }
      assert(minestzero != SCIP_INVALID);

      /* branch on the first part of the sorted candidates:
       * - for each of these variables i, create a child node x_0 = ... = x_i-1 = 0, x_i = 1
       * - create an additional child node x_0 = ... = x_n-1 = 0
       */
      for( i = 0; i < nbranchcands; ++i )
      {
         /* create child with x_0 = ... = x_i-1 = 0, x_i = 1 */
         SCIP_CALL( SCIPcreateChild(scip, &node, (SCIP_Real)nbranchcands, MIN(minestzero, estone[i])) );
         for( j = 0; j < i; ++j )
         {
            SCIP_CALL( SCIPchgVarUbNode(scip, node, branchcands[j], 0.0) );
         }
         SCIP_CALL( SCIPchgVarLbNode(scip, node, branchcands[i], 1.0) );
      }
      /* create child with x_0 = ... = x_n = 0 */
      SCIP_CALL( SCIPcreateChild(scip, &node, (SCIP_Real)nbranchcands, minestzero) );
      for( i = 0; i < nbranchcands; ++i )
      {
         SCIP_CALL( SCIPchgVarUbNode(scip, node, branchcands[i], 0.0) );
      }

      *result = SCIP_BRANCHED;

      SCIPfreeBufferArray(scip, &estone);

#ifdef SCIP_DEBUG
      {
         int nchildren;
         SCIP_CALL( SCIPgetChildren(scip, NULL, &nchildren) );
         SCIPdebugMessage("branched on pseudo solution: %d children\n", nchildren);
      }
#endif
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &canduses);
   SCIPfreeBufferArray(scip, &branchcands);

   return SCIP_OKAY;
}
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpSetppc)
{  /*lint --e{715}*/
   SCIP_Bool cutoff;
   SCIP_Bool separated;
   SCIP_Bool reduceddom;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   SCIPdebugMessage("LP enforcing %d set partitioning / packing / covering constraints\n", nconss);

   *result = SCIP_FEASIBLE;

   cutoff = FALSE;
   separated = FALSE;
   reduceddom = FALSE;

   /* check all useful set partitioning / packing / covering constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff && !reduceddom; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], NULL, &cutoff, &separated, &reduceddom) );
   }

   /* check all obsolete set partitioning / packing / covering constraints for feasibility */
   for( c = nusefulconss; c < nconss && !cutoff && !separated && !reduceddom; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], NULL, &cutoff, &separated, &reduceddom) );
   }

#ifdef VARUSES
#ifdef BRANCHLP
   if( !cutoff && !separated && !reduceddom )
   {
      /* if solution is not integral, choose a variable set to branch on */
      SCIP_CALL( branchLP(scip, conshdlr, result) );
      if( *result != SCIP_FEASIBLE )
         return SCIP_OKAY;
   }
#endif
#endif

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( separated )
      *result = SCIP_SEPARATED;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsSetppc)
{  /*lint --e{715}*/
   SCIP_Bool cutoff;
   SCIP_Bool infeasible;
   SCIP_Bool reduceddom;
   SCIP_Bool solvelp;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   /* if the solution is infeasible anyway due to objective value, skip the constraint processing and branch directly */
#ifdef VARUSES
   if( objinfeasible )
   {
      *result = SCIP_DIDNOTRUN;
      SCIP_CALL( branchPseudo(scip, conshdlr, result) );
      return SCIP_OKAY;
   }
#endif

   SCIPdebugMessage("pseudo enforcing %d set partitioning / packing / covering constraints\n", nconss);

   *result = SCIP_FEASIBLE;

   cutoff = FALSE;
   infeasible = FALSE;
   reduceddom = FALSE;
   solvelp = FALSE;

   /* check all set partitioning / packing / covering constraints for feasibility */
   for( c = 0; c < nconss && !cutoff && !reduceddom && !solvelp; ++c )
   {
      SCIP_CALL( enforcePseudo(scip, conss[c], &cutoff, &infeasible, &reduceddom, &solvelp) );
   }

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;
   else if( solvelp )
      *result = SCIP_SOLVELP;
   else if( infeasible )
   {
      *result = SCIP_INFEASIBLE;

#ifdef VARUSES
      /* at least one constraint is violated by pseudo solution and we didn't find a better way to resolve this:
       * -> branch on pseudo solution
       */
      SCIP_CALL( branchPseudo(scip, conshdlr, result) );
#endif
   }

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckSetppc)
{  /*lint --e{715}*/
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   /* check all set partitioning / packing / covering constraints for feasibility */
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);
      if( checklprows || consdata->row == NULL || !SCIProwIsInLP(consdata->row) )
      {
         if( !checkCons(scip, consdata, sol) )
         {
            /* constraint is violated */
            *result = SCIP_INFEASIBLE;

            if( printreason )
            {
               int v;
               SCIP_Real sum = 0.0;

               SCIP_CALL( SCIPprintCons(scip, cons, NULL) );

               for( v = 0; v < consdata->nvars; ++v )
               {
                  assert(SCIPvarIsBinary(consdata->vars[v]));
                  sum += SCIPgetSolVal(scip, sol, consdata->vars[v]);
               }
               SCIPinfoMessage(scip, NULL, ";\nviolation: the right hand side is violated by by %.15g\n", ABS(sum - 1));
            }
            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropSetppc)
{  /*lint --e{715}*/
   SCIP_Bool cutoff;
   SCIP_Bool reduceddom;
   SCIP_Bool addcut;
   SCIP_Bool mustcheck;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   SCIPdebugMessage("propagating %d/%d set partitioning / packing / covering constraints\n", nusefulconss, nconss);

   cutoff = FALSE;
   reduceddom = FALSE;

   /* propagate all useful set partitioning / packing / covering constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( processFixings(scip, conss[c], &cutoff, &reduceddom, &addcut, &mustcheck) );
   }

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolSetppc)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int oldnfixedvars;
   int oldnaggrvars;
   int oldndelconss;
   int oldnchgcoefs;
   int firstchange;
   int firstclique;
   int lastclique;
   int startdelconss;
   int c;
   SCIP_Bool cutoff;
   SCIP_Bool delay;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   delay = FALSE;
   *result = SCIP_DIDNOTFIND;
   oldnfixedvars = *nfixedvars;
   oldndelconss = *ndelconss;
   oldnaggrvars = *naggrvars;
   oldnchgcoefs = *nchgcoefs;
   cutoff = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* determine whether we want to run the clique lifting procedure */
   conshdlrdata->enablecliquelifting = conshdlrdata->enablecliquelifting || conshdlrdata->updatedsetppctype
      || conshdlrdata->noldfixedvars != SCIPgetNFixedVars(scip) || conshdlrdata->noldimpls != SCIPgetNImplications(scip)
      || conshdlrdata->noldcliques != SCIPgetNCliques(scip) || conshdlrdata->noldupgrs != nconss;

   /* remember old values */
   startdelconss = *ndelconss;
   conshdlrdata->noldimpls = SCIPgetNImplications(scip);
   conshdlrdata->noldcliques = SCIPgetNCliques(scip);
   conshdlrdata->updatedsetppctype = FALSE;

   /* process constraints */
   firstchange = INT_MAX;
   firstclique = INT_MAX;
   lastclique = -1;
   for( c = 0; c < nconss && !SCIPisStopped(scip); ++c )
   {
      SCIP_CONS* cons;
      SCIP_CONSDATA* consdata;

      assert(*result != SCIP_CUTOFF);

      cons = conss[c];
      assert(cons != NULL);
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /*SCIPdebugMessage("presolving set partitioning / packing / covering constraint <%s>\n", SCIPconsGetName(cons));*/

      /* remove all variables that are fixed to zero and replace all aggregated variables */
      if( consdata->nfixedzeros > 0 || nnewaggrvars > 0 || nnewaddconss > 0 || *naggrvars > oldnaggrvars || (nrounds == 0 && SCIPgetNRuns(scip) > 1) )
      {
         SCIP_CALL( applyFixings(scip, cons) );
      }

      /** find pairs of negated variables in constraint:
       *  partitioning/packing: all other variables must be zero, constraint is redundant
       *  covering: constraint is redundant
       *
       *  find sets of equal variables in constraint:
       *  partitioning/packing: variable must be zero
       *  covering: multiple entries of variable can be replaced by single entry
       */
      SCIP_CALL( mergeMultiples(scip, cons, nfixedvars, ndelconss, nchgcoefs, &cutoff) );

      if( cutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      /* if constraint was deleted while merging, go to the next constraint */
      if( !SCIPconsIsActive(cons) )
         continue;

      /* check if constraint is already redundant or infeasible due to fixings, fix or aggregate left over variables if
       * possible
       */
      SCIP_CALL( presolvePropagateCons(scip, cons, TRUE, NULL, NULL, NULL, NULL, nfixedvars, naggrvars, ndelconss, &cutoff) );

      if( cutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      /* if constraint was deleted while propagation, go to the next constraint */
      if( !SCIPconsIsActive(cons) )
         continue;

      /* perform dual reductions */
      if( conshdlrdata->dualpresolving )
      {
         SCIP_CALL( dualPresolving(scip, cons, nfixedvars, ndelconss, result) );

         /* if dual reduction deleted the constraint we take the next */
         if( !SCIPconsIsActive(cons) )
            continue;
      }

      /* remember the first changed constraint to begin the next redundancy round with */
      if( firstchange == INT_MAX && consdata->changed )
         firstchange = c;

      /* remember the first and last constraints for which we have to add the clique information */
      if( !consdata->cliqueadded && consdata->nvars >= 2 )
      {
         if( firstclique == INT_MAX )
            firstclique = c;
         lastclique = c;
      }
   }

   /* update result pointer */
   if( oldnfixedvars < *nfixedvars || oldnaggrvars < *naggrvars || oldndelconss < *ndelconss || oldnchgcoefs < *nchgcoefs )
      *result = SCIP_SUCCESS;

   if( firstchange < nconss && conshdlrdata->presolusehashing )
   {
      /* detect redundant constraints; fast version with hash table instead of pairwise comparison */
      SCIP_CALL( detectRedundantConstraints(scip, SCIPblkmem(scip), conss, nconss, &firstchange, ndelconss, nchgsides) );
      if( oldndelconss < *ndelconss )
         *result = SCIP_SUCCESS;
   }

   /* clique lifting */
   if( conshdlrdata->cliquelifting && conshdlrdata->enablecliquelifting )
   {
      if( *nfixedvars == oldnfixedvars && *naggrvars == oldnaggrvars && *ndelconss == oldndelconss && *nchgcoefs == oldnchgcoefs && SCIPisPresolveFinished(scip) )
      {
	 SCIP_CALL( preprocessCliques(scip, conshdlrdata, conss, nconss, nrounds, &firstchange, nfixedvars, naggrvars, ndelconss, nchgcoefs, &cutoff) );
	 ++(conshdlrdata->nclqpresolve);

	 if( cutoff )
	 {
	    *result = SCIP_CUTOFF;
	    return SCIP_OKAY;
	 }
	 else if( oldnfixedvars < *nfixedvars || oldnaggrvars < *naggrvars || oldndelconss < *ndelconss || oldnchgcoefs < *nchgcoefs )
	    *result = SCIP_SUCCESS;

	 /* remember the number of fixings */
	 conshdlrdata->noldfixedvars = *nfixedvars + *naggrvars;
      }
      else
         delay = TRUE;
   }

   if( oldndelconss == *ndelconss )
   {
      /* check constraints for redundancy */
      if( conshdlrdata->presolpairwise )
      {
         SCIP_Longint npaircomparisons;
         npaircomparisons = 0;
         oldndelconss = *ndelconss;
         oldnfixedvars = *nfixedvars;

         for( c = firstchange; c < nconss && !SCIPisStopped(scip); ++c )
         {
            assert(*result != SCIP_CUTOFF);

            if( SCIPconsIsActive(conss[c]) && !SCIPconsIsModifiable(conss[c]) )
            {
               npaircomparisons += (SCIPconsGetData(conss[c])->changed) ? (SCIP_Longint) c : ((SCIP_Longint) c - (SCIP_Longint) firstchange);

               SCIP_CALL( removeRedundantConstraints(scip, conss, firstchange, c, &cutoff, nfixedvars, ndelconss, nchgsides) );
               if( cutoff )
               {
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }

               if( npaircomparisons > NMINCOMPARISONS )
               {
                  if( (*ndelconss - oldndelconss + *nfixedvars - oldnfixedvars) / ((SCIP_Real)npaircomparisons) < MINGAINPERNMINCOMPARISONS )
                     break;
                  oldndelconss = *ndelconss;
                  oldnfixedvars = *nfixedvars;
                  npaircomparisons = 0;
                  *result = SCIP_SUCCESS;
               }
            }
         }
      }
   }

   /* add clique and implication information */
   for( c = firstclique; c < lastclique && !SCIPisStopped(scip); ++c )
   {
      SCIP_CONS* cons;
      SCIP_CONSDATA* consdata;

      assert(*result != SCIP_CUTOFF);

      cons = conss[c];
      assert(cons != NULL);

      /* ignore deleted constraints */
      if( !SCIPconsIsActive(cons) )
         continue;

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      if( !consdata->cliqueadded && consdata->nvars >= 2 )
      {
         /* add a set partitioning / packing constraint as clique */
         if( consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING || consdata->setppctype == SCIP_SETPPCTYPE_PACKING ) /*lint !e641*/
         {
            SCIP_Bool infeasible;
            int ncliquebdchgs;

            SCIP_CALL( SCIPaddClique(scip, consdata->vars, NULL, consdata->nvars, &infeasible, &ncliquebdchgs) );
            *nchgbds += ncliquebdchgs;
            if( infeasible )
            {
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
         }
         else if( consdata->nvars == 2 && !SCIPconsIsModifiable(cons) )
         {
            SCIP_Bool infeasible;
            int nimplbdchgs;

            /* a two-variable set covering constraint x + y >= 1 yields the implication x == 0 -> y == 1 */
            SCIP_CALL( SCIPaddVarImplication(scip, consdata->vars[0], FALSE, consdata->vars[1],
                  SCIP_BOUNDTYPE_LOWER, 1.0, &infeasible, &nimplbdchgs) );
            *nchgbds += nimplbdchgs;
            if( infeasible )
            {
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
         }
         consdata->cliqueadded = TRUE;
      }
   }

   conshdlrdata->enablecliquelifting = FALSE;
   conshdlrdata->noldupgrs = nconss - (*ndelconss - startdelconss);

   /* return the correct result code */
   if( *result != SCIP_CUTOFF && delay )
      *result = SCIP_DELAYED;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropSetppc)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int v;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);
   assert(infervar != NULL);
   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("conflict resolving method of set partitioning / packing / covering constraint handler\n");

   if( (SCIP_SETPPCTYPE)consdata->setppctype == SCIP_SETPPCTYPE_COVERING
      || ((SCIP_SETPPCTYPE)consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING
         && SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE) > 0.5) )
   {
#ifndef NDEBUG
      SCIP_Bool confvarfound;
#endif

      /* the inference constraint is a set partitioning or covering constraint with the inference variable infered to 1.0:
       * the reason for the deduction is the assignment of 0.0 to all other variables
       */
#ifndef NDEBUG
      confvarfound = FALSE;
#endif
      for( v = 0; v < consdata->nvars; ++v )
      {
         if( consdata->vars[v] != infervar )
         {
            /* the reason variable must be assigned to zero */
            assert(SCIPvarGetUbAtIndex(consdata->vars[v], bdchgidx, FALSE) < 0.5);
            SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[v]) );
         }
#ifndef NDEBUG
         else
         {
            assert(!confvarfound);
            confvarfound = TRUE;
         }
#endif
      }
      assert(confvarfound);
   }
   else
   {
      /* the inference constraint is a set partitioning or packing constraint with the inference variable infered to 0.0:
       * the reason for the deduction is the assignment of 1.0 to a single variable
       */
      assert(SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE) < 0.5);
      for( v = 0; v < consdata->nvars; ++v )
      {
         if( SCIPvarGetLbAtIndex(consdata->vars[v], bdchgidx, FALSE) > 0.5 )
         {
            SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[v]) );
            break;
         }
      }
      assert(v < consdata->nvars);
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockSetppc)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int nlocksdown;
   int nlocksup;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   switch( consdata->setppctype )
   {
   case SCIP_SETPPCTYPE_PARTITIONING:
      nlocksdown = nlockspos + nlocksneg;
      nlocksup = nlockspos + nlocksneg;
      break;
   case SCIP_SETPPCTYPE_PACKING:
      nlocksdown = nlocksneg;
      nlocksup = nlockspos;
      break;
   case SCIP_SETPPCTYPE_COVERING:
      nlocksdown = nlockspos;
      nlocksup = nlocksneg;
      break;
   default:
      SCIPerrorMessage("unknown setppc type\n");
      return SCIP_INVALIDDATA;
   }

   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlocksdown, nlocksup) );
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#ifdef VARUSES
static
SCIP_DECL_CONSACTIVE(consActiveSetppc)
{  /*lint --e{715}*/
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPconsIsTransformed(cons));

   SCIPdebugMessage("activation information for set partitioning / packing / covering constraint <%s>\n",
      SCIPconsGetName(cons));

   /* increase the number of uses for each variable in the constraint */
   SCIP_CALL( consdataIncVaruses(scip, SCIPconshdlrGetData(conshdlr), SCIPconsGetData(cons)) );

   return SCIP_OKAY;
}
#endif


/** constraint deactivation notification method of constraint handler */
#ifdef VARUSES
static
SCIP_DECL_CONSDEACTIVE(consDeactiveSetppc)
{  /*lint --e{715}*/
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPconsIsTransformed(cons));

   SCIPdebugMessage("deactivation information for set partitioning / packing / covering constraint <%s>\n",
      SCIPconsGetName(cons));

   /* decrease the number of uses for each variable in the constraint */
   SCIP_CALL( consdataDecVaruses(scip, SCIPconshdlrGetData(conshdlr), SCIPconsGetData(cons)) );

   return SCIP_OKAY;
}
#endif

/** variable deletion method of constraint handler */
static
SCIP_DECL_CONSDELVARS(consDelvarsSetppc)
{
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL || nconss == 0 );

   if( nconss > 0 )
   {
      SCIP_CALL( performVarDeletions(scip, conshdlr, conss, nconss) );
   }

   return SCIP_OKAY;
}



/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintSetppc)
{  /*lint --e{715}*/

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );

   SCIP_CALL( consdataPrint(scip, SCIPconsGetData(cons), file) );

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopySetppc)
{  /*lint --e{715}*/
   SCIP_VAR** sourcevars;
   const char* consname;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int nvars;
   SCIP_SETPPCTYPE type;

   /* get variables and coefficients of the source constraint */
   sourcevars = SCIPgetVarsSetppc(sourcescip, sourcecons);
   nvars = SCIPgetNVarsSetppc(sourcescip, sourcecons);

   /* get setppc type */
   type = SCIPgetTypeSetppc(sourcescip, sourcecons);
   lhs = -SCIPinfinity(scip);
   rhs = SCIPinfinity(scip);

   switch( type )
   {
   case SCIP_SETPPCTYPE_PARTITIONING:
      lhs = 1.0;
      rhs = 1.0;
      break;
   case SCIP_SETPPCTYPE_PACKING:
      rhs = 1.0;
      break;
   case SCIP_SETPPCTYPE_COVERING:
      lhs = 1.0;
      break;
   default:
      SCIPerrorMessage("unknown setppc type\n");
      return SCIP_INVALIDDATA;
   }

   if( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   /* copy the logic using the linear constraint copy method */
   SCIP_CALL( SCIPcopyConsLinear(scip, cons, sourcescip, consname, nvars, sourcevars, NULL,
         lhs, rhs, varmap, consmap,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, global, valid) );
   assert(cons != NULL);

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseSetppc)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   int nvars;

   assert(scip != NULL);
   assert(success != NULL);
   assert(str != NULL);
   assert(name != NULL);
   assert(cons != NULL);

   *success = TRUE;

   nvars = 0;
   vars = NULL;

   /* check if lhs is just 0 */
   if( str[0] == '0' )
   {
      assert(str[1] == ' ');
      str += 2;
   }
   else
   {
      SCIP_Real* coefs;
      char* endptr;
      int coefssize;
      int requsize;

      /* initialize buffers for storing the coefficients */
      coefssize = 100;
      SCIP_CALL( SCIPallocBufferArray(scip, &vars,  coefssize) );
      SCIP_CALL( SCIPallocBufferArray(scip, &coefs, coefssize) );

      /* parse linear sum to get variables and coefficients */
      SCIP_CALL( SCIPparseVarsLinearsum(scip, str, vars, coefs, &nvars, coefssize, &requsize, &endptr, success) );

      if( *success && requsize > coefssize )
      {
         /* realloc buffers and try again */
         coefssize = requsize;
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars,  coefssize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &coefs, coefssize) );

         SCIP_CALL( SCIPparseVarsLinearsum(scip, str, vars, coefs, &nvars, coefssize, &requsize, &endptr, success) );
         assert(!*success || requsize <= coefssize); /* if successful, then should have had enough space now */
      }

      if( !*success )
      {
         SCIPerrorMessage("no luck in parsing linear sum '%s'\n", str);
      }
      else
         str = endptr;

      /* free coefficient array */
      SCIPfreeBufferArray(scip, &coefs);
   }

   /* remove white spaces */
   while( isspace((unsigned char)*str) )
      str++;

   if( *success )
   {
      switch( *str )
      {
         case '=' :
            SCIP_CALL( SCIPcreateConsSetpart(scip, cons, name, nvars, vars,
                  initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
            break;
         case '<' :
            SCIP_CALL( SCIPcreateConsSetpack(scip, cons, name, nvars, vars,
                  initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
            break;
         case '>' :
            SCIP_CALL( SCIPcreateConsSetcover(scip, cons, name, nvars, vars,
                  initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
            break;
         default:
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "error parsing setppc type\n");
            *success = FALSE;
            break;
      }
   }

   /* free variable array */
   SCIPfreeBufferArrayNull(scip, &vars);

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsSetppc)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( varssize < consdata->nvars )
      (*success) = FALSE;
   else
   {
      assert(vars != NULL);

      BMScopyMemoryArray(vars, consdata->vars, consdata->nvars);
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsSetppc)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   (*nvars) = consdata->nvars;
   (*success) = TRUE;

   return SCIP_OKAY;
}

/*
 * Callback methods of event handler
 */

static
SCIP_DECL_EVENTEXEC(eventExecSetppc)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   /*debugMessage("Exec method of bound change event handler for set partitioning / packing / covering constraints\n");*/

   consdata = (SCIP_CONSDATA*)eventdata;
   assert(consdata != NULL);

   eventtype = SCIPeventGetType(event);

   switch( eventtype )
   {
   case SCIP_EVENTTYPE_LBTIGHTENED:
      consdata->nfixedones++;
      consdata->propagated = FALSE;
      break;
   case SCIP_EVENTTYPE_LBRELAXED:
      consdata->nfixedones--;
      break;
   case SCIP_EVENTTYPE_UBTIGHTENED:
      consdata->nfixedzeros++;
      consdata->propagated = FALSE;
      break;
   case SCIP_EVENTTYPE_UBRELAXED:
      consdata->nfixedzeros--;
      break;
   case SCIP_EVENTTYPE_VARDELETED:
      consdata->varsdeleted = TRUE;
      break;
   case SCIP_EVENTTYPE_VARFIXED:
      if( consdata->merged )
      {
	 SCIP_VAR* var = SCIPeventGetVar(event);

	 /* this event should only arise during the presolving stage */
	 assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING);
	 assert(var != NULL);

	 /* one variable was changed to a negated or aggregated variable, so maybe we can merge again */
	 if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_FIXED && SCIPvarGetLbGlobal(var) < 0.5 && SCIPvarGetUbGlobal(var) > 0.5 )
	    consdata->merged = FALSE;
      }
      break;
   default:
      SCIPerrorMessage("invalid event type\n");
      return SCIP_INVALIDDATA;
   }
   assert(0 <= consdata->nfixedzeros && consdata->nfixedzeros <= consdata->nvars);
   assert(0 <= consdata->nfixedones && consdata->nfixedones <= consdata->nvars);

   /*debugMessage(" -> constraint has %d zero-fixed and %d one-fixed of %d variables\n",
     consdata->nfixedzeros, consdata->nfixedones, consdata->nvars);*/

   return SCIP_OKAY;
}




/*
 * Callback methods of conflict handler
 */

static
SCIP_DECL_CONFLICTEXEC(conflictExecSetppc)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   int i;

   assert(conflicthdlr != NULL);
   assert(strcmp(SCIPconflicthdlrGetName(conflicthdlr), CONFLICTHDLR_NAME) == 0);
   assert(bdchginfos != NULL || nbdchginfos == 0);
   assert(result != NULL);

   /* don't process already resolved conflicts */
   if( resolved )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   /* create array of variables in conflict constraint */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nbdchginfos) );
   for( i = 0; i < nbdchginfos; ++i )
   {
      assert(bdchginfos != NULL);

      vars[i] = SCIPbdchginfoGetVar(bdchginfos[i]);

      /* we can only treat binary variables */
      if( !SCIPvarIsBinary(vars[i]) )
         break;

      /* if the variable is fixed to one in the conflict set, we have to use its negation */
      if( SCIPbdchginfoGetNewbound(bdchginfos[i]) > 0.5 )
      {
         SCIP_CALL( SCIPgetNegatedVar(scip, vars[i], &vars[i]) );
      }
   }

   if( i == nbdchginfos )
   {
      SCIP_CONS* cons;
      char consname[SCIP_MAXSTRLEN];

      /* create a constraint out of the conflict set */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "cf%d_%"SCIP_LONGINT_FORMAT, SCIPgetNRuns(scip), SCIPgetNConflictConssApplied(scip));
      SCIP_CALL( SCIPcreateConsSetcover(scip, &cons, consname, nbdchginfos, vars,
            FALSE, separate, FALSE, FALSE, TRUE, local, FALSE, dynamic, removable, FALSE) );
      SCIP_CALL( SCIPaddConsNode(scip, node, cons, validnode) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      *result = SCIP_CONSADDED;
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for set partitioning / packing / covering constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSetppc(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EVENTHDLR* eventhdlr;

   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecSetppc, NULL) );

   /* create conflict handler for setppc constraints */
   SCIP_CALL( SCIPincludeConflicthdlrBasic(scip, NULL, CONFLICTHDLR_NAME, CONFLICTHDLR_DESC, CONFLICTHDLR_PRIORITY,
         conflictExecSetppc, NULL) );

   /* create constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata, eventhdlr) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSetppc, consEnfopsSetppc, consCheckSetppc, consLockSetppc,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
#ifdef VARUSES
   SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveSetppc) );
   SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveSetppc) );
#endif
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopySetppc, consCopySetppc) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSetppc) );
   SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsSetppc) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolSetppc) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeSetppc) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsSetppc) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsSetppc) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitSetppc) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpSetppc) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseSetppc) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSetppc, CONSHDLR_MAXPREROUNDS, CONSHDLR_DELAYPRESOL) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintSetppc) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropSetppc, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropSetppc) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpSetppc, consSepasolSetppc, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSetppc) );

   if( SCIPfindConshdlr(scip,"linear") != NULL )
   {
      /* include the linear constraint to setppc constraint upgrade in the linear constraint handler */
      SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdSetppc, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );
   }

   /* set partitioning constraint handler parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/"CONSHDLR_NAME"/npseudobranches",
         "number of children created in pseudo branching (0: disable pseudo branching)",
         &conshdlrdata->npseudobranches, TRUE, DEFAULT_NPSEUDOBRANCHES, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/presolpairwise",
         "should pairwise constraint comparison be performed in presolving?",
         &conshdlrdata->presolpairwise, TRUE, DEFAULT_PRESOLPAIRWISE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/presolusehashing",
         "should hash table be used for detecting redundant constraints in advance",
         &conshdlrdata->presolusehashing, TRUE, DEFAULT_PRESOLUSEHASHING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/dualpresolving",
         "should dual presolving steps be performed?",
         &conshdlrdata->dualpresolving, TRUE, DEFAULT_DUALPRESOLVING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/cliquelifting",
         " should we try to lift variables into other clique constraints, fix variables, aggregate them, and also shrink the amount of variables in clique constraints",
         &conshdlrdata->cliquelifting, TRUE, DEFAULT_CLIQUELIFTING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/addvariablesascliques",
         "should we try to generate extra cliques out of all binary variables to maybe fasten redundant constraint detection",
         &conshdlrdata->addvariablesascliques, TRUE, DEFAULT_ADDVARIABLESASCLIQUES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/cliqueshrinking",
         "should we try to shrink the number of variables in a clique constraints, by replacing more than one variable by only one",
         &conshdlrdata->cliqueshrinking, TRUE, DEFAULT_CLIQUESHRINKING, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a set partitioning constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsSetpart(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
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
   return createConsSetppc(scip, cons, name, nvars, vars, SCIP_SETPPCTYPE_PARTITIONING,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode);
}

/** creates and captures a set partitioning constraint with all constraint flags set
 *  to their default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicSetpart(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars                /**< array with variables of constraint entries */
   )
{
   SCIP_CALL( SCIPcreateConsSetpart(scip, cons, name, nvars, vars,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** creates and captures a set packing constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsSetpack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
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
   return createConsSetppc(scip, cons, name, nvars, vars, SCIP_SETPPCTYPE_PACKING,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode);
}

/** creates and captures a set packing constraint with all constraint flags set
 *  to their default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicSetpack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars                /**< array with variables of constraint entries */
   )
{
   SCIP_CALL( SCIPcreateConsSetpack(scip, cons, name, nvars, vars,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;

}

/** creates and captures a set covering constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsSetcover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
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
   return createConsSetppc(scip, cons, name, nvars, vars, SCIP_SETPPCTYPE_COVERING,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode);
}

/** creates and captures a set covering constraint with all constraint flags set
 *  to their default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicSetcover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars                /**< array with variables of constraint entries */
   )
{
   SCIP_CALL( SCIPcreateConsSetcover(scip, cons, name, nvars, vars,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;

}

/** adds coefficient in set partitioning / packing / covering constraint */
SCIP_RETCODE SCIPaddCoefSetppc(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var                 /**< variable to add to the constraint */
   )
{
   assert(var != NULL);

   /*debugMessage("adding variable <%s> to setppc constraint <%s>\n",
     SCIPvarGetName(var), SCIPconsGetName(cons));*/

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a set partitioning / packing / covering constraint\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( addCoef(scip, cons, var) );

   return SCIP_OKAY;
}

/** gets number of variables in set partitioning / packing / covering constraint */
int SCIPgetNVarsSetppc(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a set partitioning / packing / covering constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nvars;
}

/** gets array of variables in set partitioning / packing / covering constraint */
SCIP_VAR** SCIPgetVarsSetppc(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a set partitioning / packing / covering constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vars;
}

/** gets type of set partitioning / packing / covering constraint */
SCIP_SETPPCTYPE SCIPgetTypeSetppc(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a set partitioning / packing / covering constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return (SCIP_SETPPCTYPE)(consdata->setppctype);
}

/** gets the dual solution of the set partitioning / packing / covering constraint in the current LP */
SCIP_Real SCIPgetDualsolSetppc(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a set partitioning / packing / covering constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualsol(consdata->row);
   else
      return 0.0;
}

/** gets the dual Farkas value of the set partitioning / packing / covering constraint in the current infeasible LP */
SCIP_Real SCIPgetDualfarkasSetppc(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a set partitioning / packing / covering constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualfarkas(consdata->row);
   else
      return 0.0;
}

/** returns the linear relaxation of the given set partitioning / packing / covering constraint; may return NULL if no
 *  LP row was yet created; the user must not modify the row!
 */
SCIP_ROW* SCIPgetRowSetppc(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a set partitioning / packing / covering constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->row;
}

/** returns current number of variables fixed to one in the constraint  */
int SCIPgetNFixedonesSetppc(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a set partitioning / packing / covering constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nfixedones;
}
   

/** returns current number of variables fixed to zero in the constraint  */
int SCIPgetNFixedzerosSetppc(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a set partitioning / packing / covering constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nfixedzeros;
}


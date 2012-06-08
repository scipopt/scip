/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
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
};




/*
 * Local methods
 */

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
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocMemory(scip, conshdlrdata) );
#ifdef VARUSES
   SCIP_CALL( SCIPcreateIntarray(scip, &(*conshdlrdata)->varuses) );
#endif
   (*conshdlrdata)->npseudobranches = DEFAULT_NPSEUDOBRANCHES;

   /* get event handler for bound change events */
   (*conshdlrdata)->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( (*conshdlrdata)->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for set partitioning / packing / covering constraints not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

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
SCIP_RETCODE consdataSort(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   assert(consdata != NULL);

   if( consdata->nvars == 0 )
      consdata->sorted = TRUE;
   else if( !consdata->sorted )
   {
      SCIPsortPtr((void**)consdata->vars, SCIPvarComp, consdata->nvars);
      consdata->sorted = TRUE;
   }
   assert(consdata->sorted);

   return SCIP_OKAY;
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

   if( consdata->merged ) 
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
                     assert(fixed);
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
               assert(fixed);
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

   SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->row, SCIPconsGetName(cons), lhs, rhs,
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
   SCIP* scip;
   SCIP_CONSDATA* consdata;
   unsigned int hashval;
   int minidx;
   int mididx;
   int maxidx;
   
   consdata = SCIPconsGetData((SCIP_CONS*)key);
   assert(consdata != NULL);
   assert(consdata->nvars > 0);

   scip = (SCIP*)userptr; 
   assert(scip != NULL);

   /* sorts the constraints */
   SCIP_CALL_ABORT( consdataSort(scip, consdata) );

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

   assert(conss != NULL);
   assert(ndelconss != NULL);
   assert(nchgsides != NULL);

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
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ) );
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ) );
         
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
   SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ) );
   SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ) );

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

   /* sort the constraint */
   SCIP_CALL( consdataSort(scip, consdata0) );

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
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ) );
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ) );

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
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ) );
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ) );
         SCIP_CALL( processContainedCons(scip, cons0, cons1, cutoff, nfixedvars, ndelconss, nchgsides) );
      }
      else if( cons1iscontained )
      {
         /* cons1 is contained in cons1 */
         SCIPdebugMessage("setppc constraint <%s> is contained in <%s>\n", SCIPconsGetName(cons1), SCIPconsGetName(cons0));
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ) );
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ) );
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
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, NULL, NULL, &nlpcands) );
   if( nlpcands == 0 )
      return SCIP_OKAY;

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
      /* select the first variables from the sorted candidate list, until MAXBRANCHWEIGHT is reached;
       * then choose one less
       */
      branchweight = 0.0;
      solval = 0.0;
      for( nselcands = 0; nselcands < nsortcands && branchweight <= MAXBRANCHWEIGHT; ++nselcands )
      {
         solval = SCIPgetVarSol(scip, sortcands[nselcands]);
         assert(SCIPisFeasGE(scip, solval, 0.0) && SCIPisFeasLE(scip, solval, 1.0));
         branchweight += solval;
      }
      assert(nselcands > 0);
      nselcands--;
      branchweight -= solval;

      /* check, if we accumulated at least MIN and at most MAXBRANCHWEIGHT weight */
      if( MINBRANCHWEIGHT <= branchweight && branchweight <= MAXBRANCHWEIGHT )
      {
         SCIP_NODE* node;
         SCIP_Real downprio;

         /* perform the binary set branching on the selected variables */
         assert(1 <= nselcands && nselcands <= nlpcands);

         /* choose preferred branching direction */
         switch( SCIPvarGetBranchDirection(branchcands[0]) )
         {
         case SCIP_BRANCHDIR_DOWNWARDS:
            downprio = 1.0;
            break;
         case SCIP_BRANCHDIR_UPWARDS:
            downprio = -1.0;
            break;
         case SCIP_BRANCHDIR_AUTO:
            downprio = SCIPvarGetRootSol(branchcands[0]) - SCIPgetVarSol(scip, branchcands[0]);
            break;
         default:
            SCIPerrorMessage("invalid preferred branching direction <%d> of variable <%s>\n",
               SCIPvarGetBranchDirection(branchcands[0]), SCIPvarGetName(branchcands[0]));
            return SCIP_INVALIDDATA;
         }

         /* create left child, fix x_i = 0 for all i \in S */
         SCIP_CALL( SCIPcreateChild(scip, &node, downprio) );
         for( i = 0; i < nselcands; ++i )
         {
            SCIP_CALL( SCIPchgVarUbNode(scip, node, sortcands[i], 0.0) );
         }

         /* create right child: add constraint x(S) >= 1 */
         SCIP_CALL( SCIPcreateChild(scip, &node, -downprio) );
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
                  FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE) );
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
      /* branch on the first part of the sorted candidates:
       * - for each of these variables i, create a child node x_0 = ... = x_i-1 = 0, x_i = 1
       * - create an additional child node x_0 = ... = x_n-1 = 0
       */
      for( i = 0; i < nbranchcands; ++i )
      {
         /* create child with x_0 = ... = x_i-1 = 0, x_i = 1 */
         SCIP_CALL( SCIPcreateChild(scip, &node, (SCIP_Real)nbranchcands) );
         for( j = 0; j < i; ++j )
         {
            SCIP_CALL( SCIPchgVarUbNode(scip, node, branchcands[j], 0.0) );
         }
         SCIP_CALL( SCIPchgVarLbNode(scip, node, branchcands[i], 1.0) );
      }
      /* create child with x_0 = ... = x_n = 0 */
      SCIP_CALL( SCIPcreateChild(scip, &node, (SCIP_Real)i) );
      for( i = 0; i < nbranchcands; ++i )
      {
         SCIP_CALL( SCIPchgVarUbNode(scip, node, branchcands[i], 0.0) );
      }

      *result = SCIP_BRANCHED;

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
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
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
               SCIPinfoMessage(scip, NULL, "violation: the right hand side is violated by by %.15g\n", ABS(sum - 1));
            }
            return SCIP_OKAY;
         }
         else
         {
            SCIP_CALL( SCIPincConsAge(scip, cons) );
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
   int firstchange;
   int firstclique;
   int lastclique;
   int c;
   SCIP_Bool cutoff;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;
   oldnfixedvars = *nfixedvars;
   oldndelconss = *ndelconss;
   oldnaggrvars = *naggrvars;
   cutoff = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

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

      if( consdata->nfixedones >= 2 )
      {
         /* at least two variables are fixed to 1:
          * - a set covering constraint is feasible anyway and can be deleted
          * - a set partitioning or packing constraint is infeasible
          */
         if( consdata->setppctype == SCIP_SETPPCTYPE_COVERING ) /*lint !e641*/
         {
            SCIPdebugMessage("set covering constraint <%s> is redundant\n", SCIPconsGetName(cons));
            SCIP_CALL( SCIPdelCons(scip, cons) );
            (*ndelconss)++;
            *result = SCIP_SUCCESS;
            continue;
         }
         else
         {
            SCIPdebugMessage("set partitioning / packing constraint <%s> is infeasible\n", SCIPconsGetName(cons));
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
      }
      else if( consdata->nfixedones == 1 )
      {
         /* exactly one variable is fixed to 1:
          * - a set covering constraint is feasible anyway and can be disabled
          * - all other variables in a set partitioning or packing constraint must be zero
          */
         if( consdata->setppctype == SCIP_SETPPCTYPE_COVERING ) /*lint !e641*/
         {
            SCIPdebugMessage("set covering constraint <%s> is redundant\n", SCIPconsGetName(cons));
            SCIP_CALL( SCIPdelCons(scip, cons) );
            (*ndelconss)++;
            *result = SCIP_SUCCESS;
            continue;
         }
         else
         {
            SCIP_VAR* var;
            int v;

            SCIPdebugMessage("set partitioning / packing constraint <%s> has a variable fixed to 1.0\n", SCIPconsGetName(cons));
            for( v = 0; v < consdata->nvars; ++v )
            {
               var = consdata->vars[v];
               if( SCIPvarGetLbGlobal(var) + 0.5 < SCIPvarGetUbGlobal(var) )
               {
                  SCIP_Bool infeasible;
                  SCIP_Bool fixed;

                  SCIP_CALL( SCIPfixVar(scip, var, 0.0, &infeasible, &fixed) );
                  if( infeasible )
                  {
                     SCIPdebugMessage("setppc constraint <%s>: infeasible fixing <%s> == 0\n",
                        SCIPconsGetName(cons), SCIPvarGetName(var));
                     *result = SCIP_CUTOFF;
                     return SCIP_OKAY;
                  }
                  assert(fixed);
                  (*nfixedvars)++;
                  *result = SCIP_SUCCESS;
               }
            }

            /* now all other variables are fixed to zero:
             * the constraint is feasible, and if it's not modifiable, it is redundant
             */
            if( !SCIPconsIsModifiable(cons) )
            {
               SCIPdebugMessage("set partitioning / packing constraint <%s> is redundant\n", SCIPconsGetName(cons));
               SCIP_CALL( SCIPdelCons(scip, cons) );
               (*ndelconss)++;
               *result = SCIP_SUCCESS;
               continue;
            }
         }
      }
      else if( !SCIPconsIsModifiable(cons) )
      {
         /* all other preprocessing steps can only be done on non-modifiable constraints */
         if( consdata->nfixedzeros == consdata->nvars )
         {
            /* all variables are fixed to zero:
             * - a set packing constraint is feasible anyway and can be deleted
             * - a set partitioning or covering constraint is infeasible, and so is the whole problem
             */
            assert(consdata->nfixedones == 0);

            if( consdata->setppctype == SCIP_SETPPCTYPE_PACKING ) /*lint !e641*/
            {
               SCIPdebugMessage("set packing constraint <%s> is redundant: all variables fixed to zero\n", SCIPconsGetName(cons));
               SCIP_CALL( SCIPdelCons(scip, cons) );
               (*ndelconss)++;
               *result = SCIP_SUCCESS;
               continue;
            }
            else
            {
               SCIPdebugMessage("set partitioning / covering constraint <%s> is infeasible\n", SCIPconsGetName(cons));
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
         }
         else if( consdata->nfixedzeros == consdata->nvars - 1 )
         {
            /* all variables except one are fixed to zero:
             * - a set packing constraint is feasible anyway, and can be deleted
             * - a set partitioning or covering constraint is feasible and can be deleted after the
             *   remaining variable is fixed to one
             */
            assert(consdata->nfixedones == 0);

            if( consdata->setppctype == SCIP_SETPPCTYPE_PACKING ) /*lint !e641*/
            {
               SCIPdebugMessage("set packing constraint <%s> is redundant: all but one variable fixed to zero\n", SCIPconsGetName(cons));
               SCIP_CALL( SCIPdelCons(scip, cons) );
               (*ndelconss)++;
               *result = SCIP_SUCCESS;
               continue;
            }
            else
            {
               SCIP_VAR* var;
               SCIP_Bool infeasible;
               SCIP_Bool fixed;
               SCIP_Bool found;
               int v;

               SCIPdebugMessage("set partitioning / covering constraint <%s> has only one variable not fixed to 0.0\n",
                  SCIPconsGetName(cons));

               /* search unfixed variable */
               found = FALSE;
               var = NULL;
               for( v = 0; v < consdata->nvars && !found; ++v )
               {
                  var = consdata->vars[v];
                  found = !SCIPisZero(scip, SCIPvarGetUbGlobal(var));
               }
               assert(found);

               SCIP_CALL( SCIPfixVar(scip, var, 1.0, &infeasible, &fixed) );
               if( infeasible )
               {
                  SCIPdebugMessage("setppc constraint <%s>: infeasible fixing <%s> == 1\n",
                     SCIPconsGetName(cons), SCIPvarGetName(var));
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }
               assert(fixed);
               (*nfixedvars)++;

               SCIP_CALL( SCIPdelCons(scip, cons) );
               (*ndelconss)++;
               *result = SCIP_SUCCESS;
               continue;
            }
         }
         else if( consdata->nfixedzeros == consdata->nvars - 2
            && consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING ) /*lint !e641*/
         {
            SCIP_VAR* var;
            SCIP_VAR* var1;
            SCIP_VAR* var2;
            SCIP_Bool infeasible;
            SCIP_Bool redundant;
            SCIP_Bool aggregated;
            int v;

            /* aggregate variable and delete constraint, if set partitioning constraint consists only of two
             * non-fixed variables
             */

            /* search unfixed variable */
            var1 = NULL;
            var2 = NULL;
            for( v = 0; v < consdata->nvars && var2 == NULL; ++v )
            {
               var = consdata->vars[v];
               if( !SCIPisZero(scip, SCIPvarGetUbGlobal(var)) )
               {
                  if( var1 == NULL )
                     var1 = var;
                  else
                     var2 = var;
               }
            }
            assert(var1 != NULL && var2 != NULL);

#ifdef VARUSES
            /* in order to not mess up the variable usage counting, we have to decrease usage counting, aggregate,
             * and increase usage counting again
             */
            SCIP_CALL( conshdlrdataDecVaruses(scip, conshdlrdata, var1) );
            SCIP_CALL( conshdlrdataDecVaruses(scip, conshdlrdata, var2) );
#endif

            /* aggregate binary equality var1 + var2 == 1 */
            SCIPdebugMessage("set partitioning constraint <%s>: aggregate <%s> + <%s> == 1\n",
               SCIPconsGetName(cons), SCIPvarGetName(var1), SCIPvarGetName(var2));
            SCIP_CALL( SCIPaggregateVars(scip, var1, var2, 1.0, 1.0, 1.0, &infeasible, &redundant, &aggregated) );

#ifdef VARUSES
            /* increase variable usage counting again */
            SCIP_CALL( conshdlrdataIncVaruses(scip, conshdlrdata, var1) );
            SCIP_CALL( conshdlrdataIncVaruses(scip, conshdlrdata, var2) );
#endif

            /* evaluate aggregation result */
            if( infeasible )
            {
               SCIPdebugMessage("set partitioning constraint <%s>: infeasible aggregation <%s> + <%s> == 1\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var1), SCIPvarGetName(var2));
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }

            if( aggregated )
               (*naggrvars)++;

            if( redundant )
            {
               SCIP_CALL( SCIPdelCons(scip, cons) );
               (*ndelconss)++;
            }
            *result = SCIP_SUCCESS;
            continue;
         }
      }

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

   if( oldndelconss == *ndelconss )
   {
      if( firstchange < nconss && conshdlrdata->presolusehashing ) 
      {
         /* detect redundant constraints; fast version with hash table instead of pairwise comparison */
         SCIP_CALL( detectRedundantConstraints(scip, SCIPblkmem(scip), conss, nconss, &firstchange, ndelconss, nchgsides) );
         if( oldndelconss < *ndelconss )
            *result = SCIP_SUCCESS;
      }

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
#else
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
#else
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
      break;
   case SCIP_EVENTTYPE_LBRELAXED:
      consdata->nfixedones--;
      break;
   case SCIP_EVENTTYPE_UBTIGHTENED:
      consdata->nfixedzeros++;
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
	 /* this event should only arise during the presolving stage */
	 assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING);
	 assert(SCIPeventGetVar(event) != NULL);

	 /* one variable was changed to a negated or aggregated variable, so maybe we can merge again */
	 if( SCIPvarGetStatus(SCIPeventGetVar(event)) != SCIP_VARSTATUS_FIXED )
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

   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, NULL, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecSetppc, NULL) );

   /* create conflict handler for setppc constraints */
   SCIP_CALL( SCIPincludeConflicthdlrBasic(scip, NULL, CONFLICTHDLR_NAME, CONFLICTHDLR_DESC, CONFLICTHDLR_PRIORITY,
         conflictExecSetppc, NULL) );

   /* create constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSetppc, consEnfopsSetppc, consCheckSetppc, consLockSetppc,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopySetppc, consCopySetppc) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSetppc) );
   SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsSetppc) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolSetppc) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeSetppc) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsSetppc) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsSetppc) );
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
         "constraints/setppc/npseudobranches",
         "number of children created in pseudo branching (0: disable pseudo branching)",
         &conshdlrdata->npseudobranches, TRUE, DEFAULT_NPSEUDOBRANCHES, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/setppc/presolpairwise",
         "should pairwise constraint comparison be performed in presolving?",
         &conshdlrdata->presolpairwise, TRUE, DEFAULT_PRESOLPAIRWISE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/setppc/presolusehashing",
         "should hash table be used for detecting redundant constraints in advance", 
         &conshdlrdata->presolusehashing, TRUE, DEFAULT_PRESOLUSEHASHING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/setppc/dualpresolving",
         "should dual presolving steps be performed?",
         &conshdlrdata->dualpresolving, TRUE, DEFAULT_DUALPRESOLVING, NULL, NULL) );

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


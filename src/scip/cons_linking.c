/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_linking.c,v 1.1 2010/06/24 10:53:50 bzfheinz Exp $"

/**@file   cons_linking.c
 * @brief  constraint handler for linking constraints
 * @author Stefan Heinz
 * @author Jens Schulz
 *
 * The constraints handler stores linking constraints between an integer variable and an array of binary variables. Such
 * a linking constraint has the form:
 *
 * intvar = sum_{i=1}^n {(offset+i) * binvars[i]}
 *
 * with the additional side condition that exactly one binary variable has to be one (set partitioning condition).
 *
 * This constraint can be created only with the integer variables. In this case the binary variables are only created on
 * demand. That is, at that point someone ask for the binary variables. Therefore, such constraints can be used to get a
 * "binary representation" of the domain of the integer variable which will be dynamically created.
 *
 *
 * @todo add pairwise comparison of constraints in presolving (fast hash table version and complete pairwise comparison)
 * @todo in case the integer variable is set to lower or upper bound it follows that only the corresponding binary
 * variable has a positive value which is one, this can be used to fasten the checking routine
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/cons_linear.h"
#include "scip/cons_linking.h"
#include "scip/cons_setppc.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "linking"
#define CONSHDLR_DESC          "linking constraint x = offset + sum_{i=1}^{n} i*y_i, y1+...+yn = 1, x integer, y's binary"

#define EVENTHDLR_NAME         "linking"
#define EVENTHDLR_DESC         "event handler for linking constraints"

#define CONSHDLR_SEPAPRIORITY    750000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -2050000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -750000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */

#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */

#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */


//#define LINCONSUPGD_PRIORITY   +750000 /**< priority of the constraint handler for upgrading of linear constraints */

#define HASHSIZE_BINVARSCONS     131101 /**< minimal size of hash table in linking constraint handler */

#define DEFAULT_LINEARIZE         FALSE /**< should the linking constraint be linearize after the binary variable are created */

/*
 * Data structures
 */

/** constraint data for linking constraints */
struct SCIP_ConsData
{
   SCIP_VAR*             intvar;             /**< integer variable which is linked */
   SCIP_VAR**            binvars;            /**< binary variables */
   SCIP_ROW*             row1;               /**< LP row for the linking itself */
   SCIP_ROW*             row2;               /**< LP row ensuring the set partitioning condition of the binary variables */
   int                   nbinvars;           /**< number of binary variables */
   int                   offset;             /**< offset of the binary representation */
   int                   nfixedzeros;        /**< current number of variables fixed to zero in the constraint */
   int                   nfixedones;         /**< current number of variables fixed to one in the constraint */
   unsigned int          cliqueadded:1;      /**< was the set partitioning condition already added as clique? */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events on binary variables */
   SCIP_HASHMAP*         varmap;             /**< hash map mapping an integer variable to its linking constraint */
   SCIP_Bool             linearize;          /**< should the linking constraint be linearize after the binary variable are created */
};

/*
 * Local methods
 */

/** returns for a given integer variable the corresponding hash map key */
static
void* getHashmapKey(
   SCIP_VAR*             var                 /**< variable to get the hash map key for */
   )
{
   /* return the unique variable index + 1 */
   return (void*)(size_t)(SCIPvarGetIndex(var) + 1);
}

/** installs rounding locks for the binary variables in the given linking constraint */
static
SCIP_RETCODE lockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint */
   SCIP_VAR**            binvars,            /**< binary variables  */
   int                   nbinvars            /**< number of binary variables */
   )
{
   int b;
   
   for( b = 0; b < nbinvars; ++b )
   {
      SCIP_CALL( SCIPlockVarCons(scip, binvars[b], cons, TRUE, TRUE) );
   }
   
   return SCIP_OKAY;
}

/** removes rounding locks from the given binary variables */
static
SCIP_RETCODE unlockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint */
   SCIP_VAR**            binvars,            /**< binary variables  */
   int                   nbinvars            /**< number of binary variables */
   )
{
   int b;

   for( b = 0; b < nbinvars; ++b )
   {
      SCIP_CALL( SCIPunlockVarCons(scip, binvars[b], cons, TRUE, TRUE) );
   }
   return SCIP_OKAY;
}

/** creates constaint handler data for the linking constraint handler */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocMemory(scip, conshdlrdata) );

   /* create hash map */
   SCIP_CALL( SCIPhashmapCreate(&(*conshdlrdata)->varmap, SCIPblkmem(scip), HASHSIZE_BINVARSCONS) );
   
   /* get event handler for bound change events on binary variables */
   (*conshdlrdata)->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( (*conshdlrdata)->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for "CONSHDLR_NAME" constraints not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   return SCIP_OKAY;
}

/** frees constraint handler data for linking constraint handler */
static
SCIP_RETCODE conshdlrdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);
   assert((*conshdlrdata)->varmap != NULL);

   /* free hash map */
   SCIPhashmapFree(&(*conshdlrdata)->varmap);
   
   /* free memory of constraint handler data */
   SCIPfreeMemory(scip, conshdlrdata);

   return SCIP_OKAY;
}

/** prints linking constraint to file stream */
static
void consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linking constraint data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_VAR** binvars;
   SCIP_VAR* intvar;
   int nbinvars;
   int offset;
   int b;

   assert(scip != NULL);
   assert(consdata != NULL);

   intvar = consdata->intvar;
   binvars = consdata->binvars;
   nbinvars = consdata->nbinvars;
   offset = consdata->offset;
   
   assert(intvar != NULL);
   assert(binvars != NULL || nbinvars == 0);

   /* print coefficients */
   SCIPinfoMessage(scip, file, "<%s> = ", SCIPvarGetName(intvar));

   if( nbinvars == 0 )
   {
      SCIPinfoMessage(scip, file, "no binary variables yet");
   }
   
   for( b = 0; b < nbinvars; ++b )
   {
      SCIPinfoMessage(scip, file, "%+d<%s> ", offset + b, SCIPvarGetName(binvars[b]));
   }
}

/** catches events for variable at given position */
static
SCIP_RETCODE catchEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linking constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< array position of variable to catch bound change events for */
   )
{
   SCIP_VAR* var;

   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nbinvars);
   assert(consdata->binvars != NULL);

   var = consdata->binvars[pos];
   assert(var != NULL);

   /* catch bound change events on variable */
   /**@todo do we have to add the event SCIP_EVENTTYPE_VARFIXED? */
   SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );

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
   SCIP_CONSDATA*        consdata,           /**< linking constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< array position of variable to catch bound change events for */
   )
{
   SCIP_VAR* var;
   
   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nbinvars);
   assert(consdata->binvars != NULL);

   var = consdata->binvars[pos];
   assert(var != NULL);

   /* drop events on variable */
   SCIP_CALL( SCIPdropVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );

   /* update the fixed variables counters for this variable */
   if( SCIPisEQ(scip, SCIPvarGetUbLocal(var), 0.0) )
      consdata->nfixedzeros--;
   else if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), 1.0) )
      consdata->nfixedones--;

   return SCIP_OKAY;
}

/** catches bound change events for all variables in transformed linking constraint */
static
SCIP_RETCODE catchAllEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linking constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(consdata != NULL);
   
   if( consdata->nbinvars <= 1 )
      return SCIP_OKAY;

   /* catch event for every single variable */
   for( i = 0; i < consdata->nbinvars; ++i )
   {
      SCIP_CALL( catchEvent(scip, consdata, eventhdlr, i) );
   }

   return SCIP_OKAY;
}

/** drops bound change events for all variables in transformed linking constraint */
static
SCIP_RETCODE dropAllEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linking constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(consdata != NULL);

   if( consdata->nbinvars <= 1 )
      return SCIP_OKAY;

   /* drop event of every single variable */
   for( i = 0; i < consdata->nbinvars; ++i )
   {
      SCIP_CALL( dropEvent(scip, consdata, eventhdlr, i) );
   }
   
   return SCIP_OKAY;
}

/** linearize the given linking constraint into a set partitioning constraint for the binary variables and a linear
 *  constraint for the linking between the integer variable and the binary variables */
static
SCIP_RETCODE consdataLinearize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint */
   SCIP_CONSDATA*        consdata            /**< linking constraint data */
   )
{
   SCIP_CONS* lincons;
   SCIP_Real offset;
   int b;
   
   /* create set partitioning constraint for the binary variables */
   SCIP_CALL( SCIPcreateConsSetpart(scip, &lincons, SCIPconsGetName(cons), consdata->nbinvars, consdata->binvars, 
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons), 
         SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
         SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
   SCIP_CALL( SCIPaddCons(scip, lincons) );
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

   /* create linear constraint for the linking between the binary variables and the integer variable */
   SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, SCIPconsGetName(cons), 0, NULL, NULL, 0.0, 0.0, 
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons), 
         SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
         SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

   offset = (SCIP_Real)consdata->offset;

   for( b = 0; b < consdata->nbinvars; ++b )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, lincons, consdata->binvars[b], offset+b) );
   }
   SCIP_CALL( SCIPaddCoefLinear(scip, lincons, consdata->intvar, -1.0) ); 

   SCIP_CALL( SCIPaddCons(scip, lincons) );
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

   return SCIP_OKAY;
}
 
/** creates the binary  variables */
static
SCIP_RETCODE consdataCreateBinvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint */
   SCIP_CONSDATA*        consdata,           /**< linking constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler for bound change events on binary variables */
   SCIP_Bool             linearize           /**< should the linking constraint be linearized */
   )
{
   SCIP_VAR* intvar;
   SCIP_VAR* binvar;
   char name[SCIP_MAXSTRLEN];
   int nbinvars;
   int lb;
   int ub;
   int b;
      
   assert(scip != NULL);
   assert(consdata != NULL);

   SCIPdebugMessage("create binary variables for integer variable <%s>\n", SCIPvarGetName(consdata->intvar));

   intvar = consdata->intvar;
   lb = (int)(SCIPvarGetLbGlobal(intvar) + 0.5);
   ub = (int)(SCIPvarGetUbGlobal(intvar) + 0.5);
   nbinvars = ub-lb+1;
   
   /* allocate block memory for the binary variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->binvars, nbinvars) );
   
   /* check if the integer variable is fixed */
   if( nbinvars == 1 )
   {
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s[%d]", SCIPvarGetName(intvar), lb);

      /* creates and captures a fixed binary variables */
      SCIP_CALL( SCIPcreateVar(scip, &binvar, name, 1.0, 1.0, 0.0, SCIP_VARTYPE_BINARY, 
            FALSE, FALSE, NULL, NULL, NULL, NULL) ); 
      SCIP_CALL( SCIPaddVar(scip, binvar) );
      
      /* change branching priority, such that the integer variable has a higher priority as the binary variable */
      SCIP_CALL( SCIPchgVarBranchPriority(scip, binvar, -1.0) );
      
      consdata->binvars[0] = binvar;
      SCIP_CALL( SCIPreleaseVar(scip, &binvar) );
   }
   else
   {
      for( b = lb; b <= ub; ++b)
      {
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s[%d]", SCIPvarGetName(intvar), b);
         
         /* creates and captures variables */
         SCIP_CALL( SCIPcreateVar(scip, &binvar, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
               TRUE, FALSE, NULL, NULL, NULL, NULL) );
         
         /* add variable to the problem */
         SCIP_CALL( SCIPaddVar(scip, binvar) );
         consdata->binvars[b-lb] = binvar;
         SCIP_CALL( SCIPreleaseVar(scip, &binvar) );
      }
   }
   
   consdata->nbinvars = nbinvars;
   consdata->offset = lb;

   assert(consdata->nfixedzeros == 0);
   assert(consdata->nfixedones == 0);

   if( SCIPisTransformed(scip) )
   {
      /* (rounding) lock binary variable */
      SCIP_CALL( lockRounding(scip, cons, consdata->binvars, consdata->nbinvars) );
      
      /* catch bound change events of variables */
      SCIP_CALL( catchAllEvents(scip, consdata, eventhdlr) );

      if( nbinvars > 1 )
      {
         if( linearize )
         {
            SCIP_CALL( consdataLinearize(scip, cons, consdata) );
         }
         else
         {
            /* enable constraint */
            SCIP_CALL( SCIPenableCons(scip, cons) );
         }
      }
   }

   return SCIP_OKAY;
}   


/** creates consdata */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure  */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_CONSDATA**       consdata,           /**< pointer to constraint data */
   SCIP_VAR*             intvar,             /**< integer variable which is linked */
   SCIP_VAR**            binvars,            /**< binary variables */
   int                   nbinvars,           /**< number of binary starting variables */
   int                   offset              /**< offset ot the binary variable representation */
   )
{
   assert(scip!= NULL);
   assert(consdata != NULL);
   assert(intvar != NULL);
   assert(binvars != NULL || nbinvars == 0);
   assert(SCIPvarGetType(intvar) != SCIP_VARTYPE_CONTINUOUS);
   
   /* allocate memory for consdata */
   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   (*consdata)->intvar = intvar;
   (*consdata)->row1 = NULL;
   (*consdata)->row2 = NULL;
   (*consdata)->cliqueadded = FALSE;
   (*consdata)->nbinvars = 0;
   (*consdata)->offset = 0;
   (*consdata)->nfixedzeros = 0;
   (*consdata)->nfixedones = 0;

   if( binvars == NULL )
   {
      (*consdata)->binvars = NULL;
   }
   else
   {
      /* copy binary variable array */
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->binvars, binvars, nbinvars) );
   }

   (*consdata)->nfixedones = 0;
   (*consdata)->nfixedzeros = 0;
   (*consdata)->nbinvars = nbinvars;
   (*consdata)->offset = offset;

   /* get transformed variable, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      if( binvars != NULL )
      {
         SCIP_CALL( SCIPgetTransformedVars(scip, nbinvars, (*consdata)->binvars, (*consdata)->binvars) );
      
         /* catch bound change events of variables */
         SCIP_CALL( catchAllEvents(scip, *consdata, eventhdlr) );
      }
   
      SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->intvar, &(*consdata)->intvar) );
   }

   return SCIP_OKAY;
}


/** free consdata */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure  */
   SCIP_CONSDATA**       consdata            /**< pointer to consdata */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);
   assert((*consdata)->nbinvars == 0 || (*consdata)->binvars != NULL);

   /* release the rows */
   if( (*consdata)->row1 != NULL )
   {
      assert((*consdata)->row2 != NULL);

      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->row1) );
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->row2) );
   }

   /* free binary variable array */
   if( (*consdata)->nbinvars > 0 )
   {
      /* if constraint belongs to transformed problem space, drop bound change events on variables */
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->binvars, (*consdata)->nbinvars);
   }
   
   /* check if the fixed counters a reseted */
   assert((*consdata)->nfixedzeros == 0);
   assert((*consdata)->nfixedones == 0);
   
   /* free constraint data */
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** analyzes conflicting assignment on given constraint where reason comes from the integer variable lower or upper
 *  bound  */
static
SCIP_RETCODE analyzeConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint to be processed */
   SCIP_VAR*             intvar,             /**< integer variable  */   
   SCIP_VAR*             binvar,             /**< binary variable is the reason */   
   SCIP_Bool             lb,                 /**< lower bound of integer variable is the reason */
   SCIP_Bool             ub                  /**< upper bound of integer variable is the reason */
   )
{
   assert(scip != NULL);

   /* conflict analysis can only be applied in solving stage */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;
   
   /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
   SCIP_CALL( SCIPinitConflictAnalysis(scip) );

   if( lb )
   {
      assert(intvar != NULL);
      SCIP_CALL( SCIPaddConflictLb(scip, intvar, NULL) );
   }

   if( ub )
   {
      assert(intvar != NULL);
      SCIP_CALL( SCIPaddConflictUb(scip, intvar, NULL) );
   }
   
   if( binvar != NULL )
   {
      SCIP_CALL( SCIPaddConflictBinvar(scip, binvar) );
   }

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}  

/** fix integer variable to offset + pos */
static
SCIP_RETCODE consFixInteger(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint to be processed */
   int                   pos,                /**< position of binary variable */
   SCIP_Bool*            cutoff              /**< pointer to store TRUE, if the node can be cut off */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* intvar;
   int offset;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;
 
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   intvar = consdata->intvar;
   offset = consdata->offset;
         
   /* change lower bound of the integer variable */
   SCIP_CALL( SCIPinferVarLbCons(scip, intvar, (SCIP_Real)(pos+offset), cons, pos, TRUE, &infeasible, &tightened) );
   
   if( infeasible )
   {
      assert(pos+offset > (int)(SCIPvarGetUbLocal(intvar)+0.5) );
      assert(pos+offset >= (int)(SCIPvarGetLbLocal(intvar)+0.5) );
      /* conflict analysis can only be applied in solving stage */
      SCIP_CALL( analyzeConflict(scip, cons, intvar, consdata->binvars[pos], FALSE, TRUE) );
                  
      *cutoff = TRUE;
      return SCIP_OKAY;          
   }
   assert(pos+offset <= (int)(SCIPvarGetUbLocal(intvar)+0.5));

   /* change upper bound of the integer variable */
   SCIP_CALL( SCIPinferVarUbCons(scip, intvar, (SCIP_Real)(pos+offset), cons, pos, TRUE, &infeasible, &tightened) );
   
   if( infeasible )
   {
      assert(pos+offset < (int)(SCIPvarGetLbLocal(intvar)+0.5) );
      assert(pos+offset <= (int)(SCIPvarGetUbLocal(intvar)+0.5) );
      /* conflict analysis can only be applied in solving stage */
      SCIP_CALL( analyzeConflict(scip, cons, intvar, consdata->binvars[pos], TRUE, FALSE) );
                  
      *cutoff = TRUE;
      return SCIP_OKAY;          
   }
   
   assert((int)(SCIPvarGetUbLocal(intvar)+0.5) == (int)(SCIPvarGetLbLocal(intvar)+0.5) );
               
   return SCIP_OKAY;
}

/** checks constraint for violation from the local bound of the integer variable, applies fixings to the binary
 *  variables if possible */
static
SCIP_RETCODE processIntegerBoundChg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint to be processed */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*                  nchgbds,            /**< pointer to store the number of changes (foxed) variable bounds */
   SCIP_Bool*            mustcheck           /**< pointer to store whether this constraint must be checked for feasibility */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_VAR* intvar;
   int offset;
   int nbinvars;
   int lblocal;
   int ublocal;
   int b;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(nchgbds != NULL);
   assert(mustcheck != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   nbinvars = consdata->nbinvars;

   /* in case there is only at most one binary variables, the constraints should already be disabled */
   assert(nbinvars > 1);
   
   /* if more than one binary variable is fixed to one return */
   if( consdata->nfixedones > 0  || consdata->nfixedzeros == nbinvars-1 )
      return  SCIP_OKAY;

   intvar = consdata->intvar;
   assert(intvar != NULL);

   vars = consdata->binvars;
   offset = consdata->offset;

   lblocal = (int)(SCIPvarGetLbLocal(intvar) + 0.5);
   ublocal = (int)(SCIPvarGetUbLocal(intvar) + 0.5);
   assert(lblocal <= ublocal);

   //   lbglobal = (int)(SCIPvarGetLbGlobal(intvar) + 0.5);
   // ubglobal = (int)(SCIPvarGetUbGlobal(intvar) + 0.5);
   // assert(lbglobal <= ubglobal);

   /* fix binary variables to zero if not yet fixed, until local lower bound */
   for( b = offset; b < lblocal; ++b )
      //   for( b = lblocal-1 ; b >= lbglobal; --b )
   {
      assert(b - offset >= 0); /*@repaired*/
      assert(b - offset < nbinvars);/*@repaired*/
      assert(vars[b-offset] != NULL );
      
      SCIPdebugMessage("fix variable <%s> to zero due to the lower bound of the integer variable <%s> [%g,%g]\n", 
         SCIPvarGetName(vars[b-offset]), SCIPvarGetName(intvar), SCIPvarGetLbLocal(intvar), SCIPvarGetUbLocal(intvar));

      SCIP_CALL( SCIPinferBinvarCons(scip, vars[b-offset], FALSE, cons, -2, &infeasible, &tightened) );
      
      if( infeasible )
      {
         SCIP_CALL( analyzeConflict(scip, cons, intvar, vars[b-offset], TRUE, FALSE) );
         *cutoff = TRUE;
         return SCIP_OKAY;      
      } 

      if( tightened )
         (*nchgbds)++;
   }

   /* fix binary variables to zero if not yet fixed, from local upper bound + 1*/
   for( b = ublocal+1; b < nbinvars; ++b )
      //   for( b = ublocal+1; b <= ubglobal; ++b )
   {
      assert(b - offset >= 0);
      assert(b - offset < nbinvars);
      assert(vars[b-offset] != NULL );

      SCIPdebugMessage("fix variable <%s> to zero due to the upper bound of the integer variable <%s> [%g,%g]\n", 
         SCIPvarGetName(vars[b-offset]), SCIPvarGetName(intvar), SCIPvarGetLbLocal(intvar), SCIPvarGetUbLocal(intvar));
 
      SCIP_CALL( SCIPinferBinvarCons(scip, vars[b-offset], FALSE, cons, -3, &infeasible, &tightened) );

      if( infeasible )
      {
         SCIP_CALL( analyzeConflict(scip, cons, intvar, vars[b-offset], FALSE, TRUE) );
         *cutoff = TRUE;
         return SCIP_OKAY;      
      } 

      if( tightened )
         (*nchgbds)++;
   }

   *mustcheck = (*nchgbds) == 0;
   
   /* if integer variables is fixed, fix the corresponding binary variable to one */
   if( lblocal == ublocal )
   {
      SCIPdebugMessage("fix variable <%s> to one due to the fixed  integer variable <%s> [%g,%g]\n", 
         SCIPvarGetName(vars[b-offset]), SCIPvarGetName(intvar), SCIPvarGetLbLocal(intvar), SCIPvarGetUbLocal(intvar));

      SCIP_CALL( SCIPinferBinvarCons(scip, vars[lblocal-offset], TRUE, cons, -6, &infeasible, &tightened) );

      if( infeasible )
      {
         SCIP_CALL( analyzeConflict(scip, cons, intvar, vars[lblocal-offset], TRUE, TRUE) );
         *cutoff = TRUE;
         return SCIP_OKAY;      
      } 
      
      if( tightened )
         (*nchgbds)++;
      
      SCIPdebugMessage(" -> disabling linking constraint <%s>\n", SCIPconsGetName(cons));
      SCIP_CALL( SCIPdelConsLocal(scip, cons) );

      *mustcheck = FALSE;
   }
   
   return SCIP_OKAY;
}

/** tightened the integer variable due to binary variables which are fixed to zero */
static
SCIP_RETCODE tightenedIntvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint to be processed */
   SCIP_CONSDATA*        consdata,           /**< linking constraint to be processed */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*                  nchgbds,            /**< pointer to store the number of changed variable bounds */
   SCIP_Bool             removefixings       /**< should the zero fixed binary variable be removed */
   )
{
   SCIP_VAR* intvar;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;
   int offset;
   int lb;
   int ub;
   int newlb;
   int newub;
   int b;
   
   /* check if the lower and upper bound of integer variable can be adjusted */
   if( consdata->nfixedones == 1 || consdata->nfixedzeros >= consdata->nbinvars-1 )
      return SCIP_OKAY;
   
   if( *cutoff )
      return SCIP_OKAY;

   intvar = consdata->intvar;
   offset = consdata->offset;
   lb = (int)(SCIPvarGetLbLocal(intvar) + 0.5);
   ub = (int)(SCIPvarGetUbLocal(intvar) + 0.5);
   assert(lb <= ub);
      
   /* check if we can tighten the upper bound of the integer variable */
   for( b = ub-offset; b >= MAX(0, lb-offset); --b ) 
   {
      assert(b  >= 0);
      assert(b < consdata->nbinvars);
      
      if( SCIPvarGetUbLocal(consdata->binvars[b]) > 0.5 )
         break;
   }
   
   newub = b + offset;
   
   SCIP_CALL( SCIPinferVarUbCons(scip, intvar, (SCIP_Real)newub, cons, -5, TRUE, &infeasible, &tightened) );
   
   if( infeasible )
   {
      /* conflict analysis can only be applied in solving stage */
      if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      {  
         int k;
         
         SCIPdebugMessage("conflict at <%s> due to bounds and fixed binvars: [lb,ub] = [%d,%d]; b = %d; b+offset = %d,\n", 
            SCIPvarGetName(intvar), lb, ub, b, b+offset);

         SCIP_CALL( SCIPinitConflictAnalysis(scip) );
         
         /* add conflicting variables */
         SCIP_CALL( SCIPaddConflictLb(scip, intvar, NULL) );         
         SCIP_CALL( SCIPaddConflictUb(scip, intvar, NULL) );
         for( k = b+1; k <= ub-offset; ++k )
         {
            SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->binvars[k]) );
         }
         
         /* analyze the conflict */
         SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
      }
      *cutoff = TRUE;         
      return SCIP_OKAY;
   }
    
   if( tightened )
   {
      (*nchgbds)++;
      
      if( removefixings )
      {
         int nvars;

         nvars = ub - newub;
         
         /* unlock the fixe binary which we remove */
         SCIP_CALL( unlockRounding(scip, cons, &consdata->binvars[newub - offset + 1], nvars) );
         consdata->nbinvars -= nvars;
         consdata->nfixedzeros -= nvars;
      }
   }

   /* check if we can tighten the lower bound of the integer variable */
   for( b = lb-offset; b < MIN(consdata->nbinvars, ub-offset+1); ++b )
   {
      assert(b >= 0);
      assert(b < consdata->nbinvars);
      
      if( SCIPvarGetUbLocal(consdata->binvars[b]) > 0.5 )
         break;
   }            
      
   newlb = b + offset;

   SCIP_CALL( SCIPinferVarLbCons(scip, intvar, (SCIP_Real)newlb, cons, -4, TRUE, &infeasible, &tightened) );
   
   /* start conflict analysis if infeasible */
   if( infeasible ) 
   {
      if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      {
         int k;
         
         SCIPdebugMessage("conflict at <%s> due to bounds and fixed binvars: [lb,ub] = [%d,%d]; b= %d; b+offset = %d \n", 
            SCIPvarGetName(intvar), lb, ub, b, b+offset);
         
         SCIP_CALL( SCIPinitConflictAnalysis(scip) );
         
         /* add conflicting variables */
         SCIP_CALL( SCIPaddConflictLb(scip, intvar, NULL) );         
         SCIP_CALL( SCIPaddConflictUb(scip, intvar, NULL) );
         for( k = b-1; k >= lb-offset; --k )
         {
            SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->binvars[k]) );
         }
         
         /* analyze the conflict */
         SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
      }
      *cutoff = TRUE;
      return SCIP_OKAY;
   }
   
   if( tightened )
   {
      (*nchgbds)++;
      
      if( removefixings )
      {
         
      }
   }
   
   return SCIP_OKAY;
}
   
/** checks constraint for violation only looking at the fixed binary variables, applies further fixings if possible */
static
SCIP_RETCODE processBinvarFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint to be processed */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*                  nchgbds,            /**< pointer to store the number of changed variable bounds */
   SCIP_Bool*            addcut,             /**< pointer to store whether this constraint must be added as a cut */
   SCIP_Bool*            mustcheck           /**< pointer to store whether this constraint must be checked for feasibility */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(nchgbds != NULL);
   assert(addcut != NULL);
   assert(mustcheck != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nbinvars == 0 || consdata->binvars != NULL);
   assert(0 <= consdata->nfixedzeros && consdata->nfixedzeros <= consdata->nbinvars);
   assert(0 <= consdata->nfixedones && consdata->nfixedones <= consdata->nbinvars);

   /*SCIPdebugMessage("processing constraint <%s> with respect to fixed variables (%d fixed to 0.0, %d fixed to 1.0)\n",
     SCIPconsGetName(cons), consdata->nfixedzeros, consdata->nfixedones);*/

   /* in case there is only at most one binary variables, the constraints should already be disabled */
   assert(consdata->nbinvars > 1);

   if( *cutoff )
      return SCIP_OKAY;

   if( consdata->nfixedones == 1 )
   {
      /* exactly one variable is fixed to 1:
       * - all other binary variables in a set partitioning must be zero
       * - integer variable are fixed to that binary variable
       */
      if( consdata->nfixedzeros < consdata->nbinvars - 1 || 
         SCIPisLT(scip, SCIPvarGetLbLocal(consdata->intvar), SCIPvarGetUbLocal(consdata->intvar)) )
      {
         SCIP_VAR** vars;
         SCIP_VAR* var;
         SCIP_Bool fixedonefound;
         int nvars;
         int v;

         SCIPdebugMessage(" -> fixing all other variables to zero due to the set partitioning condition <%s>\n",
            SCIPconsGetName(cons));

         /* unfixed variables exist: fix them to zero;
          * this could result in additional variables fixed to one due to aggregations; in this case, the
          * constraint is infeasible in local bounds
          */
         vars = consdata->binvars;
         nvars = consdata->nbinvars;
         fixedonefound = FALSE;
         for( v = 0; v < nvars && consdata->nfixedones == 1 && !(*cutoff); ++v )
         {
            var = vars[v];
            assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
            if( SCIPvarGetLbLocal(var) < 0.5 )
            {
               SCIP_CALL( SCIPinferBinvarCons(scip, var, FALSE, cons, -1, &infeasible, &tightened) );
               assert(!infeasible);
               SCIPdebugMessage("   -> fixed <%s> to zero (tightened=%u)\n", SCIPvarGetName(var), tightened);
            }
            else
            {
               fixedonefound = TRUE;
               /* fix integer variable */
               SCIP_CALL( consFixInteger(scip, cons, v, cutoff) );
            }
            
         }
         if( !(*cutoff) ) 
         {
            /* the fixed to one variable must have been found, and at least one variable must have been fixed */
            assert(consdata->nfixedones >= 1 || fixedonefound); 
            
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            (*nchgbds)++;
         }
      }

      /* now all other variables are fixed to zero:
       * the constraint is feasible, and if it's not modifiable, it is redundant
       */
      if( !SCIPconsIsModifiable(cons) && consdata->nfixedones == 1 )
      {
         SCIPdebugMessage(" -> disabling set linking constraint <%s>\n", SCIPconsGetName(cons));
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      }
   }
   else if( consdata->nfixedones >= 2 )
   {
      /* at least two variables are fixed to 1:
       * - the set partitioning condition is violated
       */
      SCIPdebugMessage(" -> conflict on "CONSHDLR_NAME" constraint <%s> due to the set partitioning condition\n", SCIPconsGetName(cons));

      SCIP_CALL( SCIPresetConsAge(scip, cons) );

      /* conflict analysis can only be applied in solving stage */
      if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      {
         SCIP_VAR** vars;
         int nvars;
         int n;
         int v;

         vars = consdata->binvars;
         nvars = consdata->nbinvars;

         /* initialize conflict analysis, and add the two variables assigned to one to conflict candidate queue */
         SCIP_CALL( SCIPinitConflictAnalysis(scip) );
         n = 0;

         for( v = 0; v < nvars && n < 2; ++v )
         {
            if( SCIPvarGetLbLocal(vars[v]) > 0.5 )
            {
               SCIP_CALL( SCIPaddConflictBinvar(scip, vars[v]) );
               n++;
            }
         }
         assert(n == 2);

         /* analyze the conflict */
         SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
      }

      *cutoff = TRUE;
   }
   else if( consdata->nfixedzeros == consdata->nbinvars )
   {
      /* all variables are fixed to zero:
       * - the set partitioning condition is violated, and if it's unmodifiable, the node
       *   can be cut off -- otherwise, the constraint must be added as a cut and further pricing must
       *   be performed
       */
      assert(consdata->nfixedones == 0);

      SCIPdebugMessage(" -> "CONSHDLR_NAME" constraint <%s> is infeasible due to the set partitioning condition\n", 
         SCIPconsGetName(cons));
      
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      if( SCIPconsIsModifiable(cons) )
         *addcut = TRUE;
      else 
      {
         /* conflict analysis can only be applied in solving stage */
         if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
         {
            SCIP_VAR** vars;
            int nvars;
            int v;

            vars = consdata->binvars;
            nvars = consdata->nbinvars;

            /* initialize conflict analysis, add all variables of infeasible constraint to conflict candidate queue */
            SCIP_CALL( SCIPinitConflictAnalysis(scip) );
            for( v = 0; v < nvars; ++v )
            {
               assert(SCIPvarGetUbLocal(vars[v]) < 0.5);
               SCIP_CALL( SCIPaddConflictBinvar(scip, vars[v]) );
            }
            
            /* analyze the conflict */
            SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
         }
         *cutoff = TRUE;
      }
   }
   else if( consdata->nfixedzeros == consdata->nbinvars - 1 && consdata->nfixedones == 0 )
   {
      /* all variables except one are fixed to zero:
       * - an unmodifiable set partitioning constraint is feasible and can be disabled after the
       *   remaining variable is fixed to one
       * - a modifiable set partitioning constraint must be checked manually
       */
      if( !SCIPconsIsModifiable(cons) )
      {
         SCIP_VAR** vars;
         SCIP_VAR* var;
         int nvars;
         int v;

         /* search the single variable that can be fixed */
         vars = consdata->binvars;
         nvars = consdata->nbinvars;
         for( v = 0; v < nvars && !(*cutoff); ++v )
         {
            var = vars[v];
            assert(SCIPisZero(scip, SCIPvarGetLbLocal(var)));
            assert(SCIPisZero(scip, SCIPvarGetUbLocal(var)) || SCIPisEQ(scip, SCIPvarGetUbLocal(var), 1.0));
            if( SCIPvarGetUbLocal(var) > 0.5 )
            {
               assert(SCIPvarGetLbLocal(var) < 0.5);
               SCIPdebugMessage(" -> fixing remaining binary variable <%s> to one in "CONSHDLR_NAME" constraint <%s>\n",
                  SCIPvarGetName(var), SCIPconsGetName(cons));
               SCIP_CALL( SCIPinferBinvarCons(scip, var, TRUE, cons, -1, &infeasible, &tightened) );
               assert(!infeasible);
               assert(tightened);
               
               /* fix integer variable */
               SCIP_CALL( consFixInteger(scip, cons, v, cutoff) );
               break;
            }
         }
         assert(v < nvars);
         assert(consdata->nfixedzeros == consdata->nbinvars - 1);
         assert(consdata->nfixedones == 1);

         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
         (*nchgbds)++;
      }
   }
   else
   {
      SCIP_CALL( tightenedIntvar(scip, cons, consdata, cutoff, nchgbds, FALSE) );
   }
   
   *mustcheck = (*nchgbds) == 0;

   assert(consdata->nfixedzeros + consdata->nfixedones <= consdata->nbinvars);
   
   return SCIP_OKAY;
}

/** returns whether the given solution is feasible for the given linking constraint */
static
SCIP_Bool checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint to be checked */
   SCIP_SOL*             sol                 /**< primal solution, or NULL for current LP/pseudo solution */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** binvars;
   SCIP_Real solval;
   SCIP_Real linksum;
   SCIP_Real linksumbound;
   SCIP_Real setpartsum;
   SCIP_Real setpartsumbound;
   int nbinvars;
   int offset;
   int b;
   
   assert(scip != NULL);
   assert(cons != NULL);

   SCIPdebugMessage("checking linking constraint <%s> for feasibility of solution %p\n", SCIPconsGetName(cons), (void*)sol);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->binvars != NULL || consdata->nbinvars == 0);
   
   /* in case there is only at most one binary variables, the constraints should already be disabled */
   assert(consdata->nbinvars > 1);
      
   /* calculate the constraint's activity for the linking part and the set partitioning part */
   binvars = consdata->binvars;
   nbinvars = consdata->nbinvars;
   offset = consdata->offset;
   linksum = 0.0;
   setpartsum = 0.0;
   linksumbound = SCIPgetSolVal(scip, sol, consdata->intvar)  + 2*SCIPfeastol(scip);
   setpartsumbound = 1.0 + 2*SCIPfeastol(scip);
  
   for( b = 0; b < nbinvars && setpartsum < setpartsumbound && linksum < linksumbound; ++b )  /* if sum >= sumbound, the feasibility is clearly decided */
   {
      assert(SCIPvarGetType(binvars[b]) == SCIP_VARTYPE_BINARY);
      solval = SCIPgetSolVal(scip, sol, binvars[b]);
      assert(SCIPisFeasGE(scip, solval, 0.0) && SCIPisFeasLE(scip, solval, 1.0));
      linksum += (offset + b) * solval;
      setpartsum += solval;
   } 

   /* check if the fixed binary variable match with the integer variable */
   return SCIPisFeasEQ(scip, linksum, SCIPgetSolVal(scip, sol, consdata->intvar)) && SCIPisFeasEQ(scip, setpartsum, 1.0);
}

/** transfer aggregated integer variables to the corresponding binary variables */
static
SCIP_RETCODE aggregateVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< hash map mapping a integer variables to its linking constraint */
   SCIP_CONS**           conss,              /**< array of linking constraint */
   int                   nconss,             /**< number of linking constraints */
   int*                  naggrvars,          /**< pointer to store the number of aggregate variables */
   SCIP_Bool*            cutoff              /**< pointer to store if a cutoff was detected */
   )
{
   SCIP_CONS* aggrcons;
   SCIP_CONSDATA* aggrconsdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR** binvars;
   SCIP_VAR** aggrbinvars;
   SCIP_VAR* intvar;
   SCIP_VAR* aggrvar;
   SCIP_Real aggrconst;
   SCIP_Real aggrscalar;
   SCIP_Bool infeasible;
   SCIP_Bool redundant;
   SCIP_Bool aggregated;
   int offset;
   int aggroffset;
   int nbinvars;
   int shift;
   int b;
   int c;

   assert(varmap != NULL);
   
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      
      intvar = consdata->intvar;
      assert(intvar != NULL);
      
      if( SCIPvarGetStatus(intvar) == SCIP_VARSTATUS_AGGREGATED )
      {
         aggrvar =  SCIPvarGetAggrVar(intvar);
         aggrcons = (SCIP_CONS*) SCIPhashmapGetImage(varmap, getHashmapKey(aggrvar));

         /* check if the aggregate variable belongs to a linking constraint */
         if( aggrcons != NULL )
         {
            aggrconsdata = SCIPconsGetData(aggrcons);
            assert(aggrcons != NULL);

            aggrconst = SCIPvarGetAggrConstant(intvar);
            aggrscalar = SCIPvarGetAggrScalar(intvar);
            
            /**@todo extend the aggregation for those cases were the aggrscalar is not equal to 1.0 */
            if( SCIPisEQ(scip, aggrscalar, 1.0 ) )
            {
               /* since both variables are integer variable and the aggrscalar is 1.0 the aggrconst should 
                * integral */
               assert(SCIPisIntegral(scip, aggrconst));
               shift = (int)(aggrconst + 0.5);

               offset = consdata->offset;
               binvars = consdata->binvars;
               aggroffset = aggrconsdata->offset;
               aggrbinvars = aggrconsdata->binvars;
               
               nbinvars = MIN(consdata->nbinvars + offset, aggrconsdata->nbinvars + shift + aggroffset);
               
               for( b = MAX(offset, aggroffset-shift); b < nbinvars; ++b )
               {
                  assert(b - offset >= 0);
                  assert(b + shift - aggroffset >= 0);
                  assert(b < consdata->nbinvars);
                  assert(b < aggrconsdata->nbinvars - shift);

                  /* add aggregation x - y  = 0.0 */
                  SCIP_CALL( SCIPaggregateVars(scip, binvars[b-offset], aggrbinvars[b+shift-aggroffset], 1.0, -1.0, 0.0, 
                        &infeasible, &redundant, &aggregated) );

                  if( infeasible )
                  {
                     (*cutoff) = TRUE;
                     return SCIP_OKAY;
                  }
                  
                  if( aggregated )
                     (*naggrvars)++;
               }
            }
         }
      }
   }
   
   return SCIP_OKAY;
}  

/** create two rows for the linking constraint 
 *
 *  - row1: {sum_{b=1}^n-1 b * binvars[b]} - intvar = -offset 
 *  - row2: {sum_{b=0}^n-1 binvars[b]} = 1.0 
 */
static
SCIP_RETCODE createRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linking constraint */
   )
{
   SCIP_CONSDATA* consdata;
   char rowname[SCIP_MAXSTRLEN];
   int b;
   
   assert( cons != NULL);

   /* get constraint data */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row1 == NULL);
   assert(consdata->row2 == NULL);
   assert(consdata->nbinvars > 1);

   /* create the LP row which capturs the linking between the integer and binary variables */
   (void)SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s[link]", SCIPconsGetName(cons));

   SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->row1, rowname, -(SCIP_Real)consdata->offset, -(SCIP_Real)consdata->offset,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
   
   /* add integer variable to the row */
   assert(consdata->intvar != NULL);
   SCIP_CALL( SCIPaddVarToRow(scip, consdata->row1, consdata->intvar, -1.0) );

   /* adding all except the first binary variable to the row */
   assert(consdata->binvars != NULL);
   for( b = 1; b < consdata->nbinvars; ++b )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->row1, consdata->binvars[b], (SCIP_Real)b) );
   }
   
   /* create the LP row which captures the set partitioning condition of the binary variables */
   (void)SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s[setppc]", SCIPconsGetName(cons));
   assert( consdata->nbinvars > 0 );   

   SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->row2, rowname, 1.0, 1.0,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
   
   SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, consdata->row2, consdata->nbinvars, consdata->binvars, 1.0) );

   return SCIP_OKAY;
}


/** adds linking constraint as cut to the LP */
static
SCIP_RETCODE addCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint */
   SCIP_SOL*             sol                 /**< primal CIP solution, NULL for current LP solution */
   )
{
   SCIP_CONSDATA* consdata;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   /* in case there is only at most one binary variables, the constraints should already be disabled */
   assert(consdata->nbinvars > 1);

   if( consdata->row1 == NULL )
   {
      assert(consdata->row2 == NULL);

      /* convert linking data into LP rows */
      SCIP_CALL( createRows(scip, cons) );
   }
   assert(consdata->row1 != NULL);
   assert(consdata->row2 != NULL);
   
   /* insert LP linking row as cut */
   if( !SCIProwIsInLP(consdata->row1) )
   {
      SCIPdebugMessage("adding linking row of constraint <%s> as cut to the LP\n", SCIPconsGetName(cons));
      SCIP_CALL( SCIPaddCut(scip, sol, consdata->row1, TRUE/*FALSE*/) );
   }

   /* insert LP set partitioning row as cut */
   if( !SCIProwIsInLP(consdata->row2) )
   {
      SCIPdebugMessage("adding set partitioning row of constraint <%s> as cut to the LP\n", SCIPconsGetName(cons));
      SCIP_CALL( SCIPaddCut(scip, sol, consdata->row2, TRUE/*FALSE*/) );
   }
   
   return SCIP_OKAY;
}


/** checks constraint for violation, and adds it as a cuts if possible */
static
SCIP_RETCODE separateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint to be separated */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   SCIP_Bool*            separated,          /**< pointer to store TRUE, if a cut was found */
   int*                  nchgbds             /**< pointer to store the number of changed variables bounds */
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
   assert(nchgbds != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* in case there is only at most one binary variables, the constraints should already be disabled */
   assert(consdata->nbinvars > 1);

   SCIPdebugMessage("separating constraint <%s>\n", SCIPconsGetName(cons));
   
   addcut = FALSE;
   mustcheck = TRUE;

   /* check constraint for violation only looking at the fixed variables, apply further fixings if possible */
   if( sol == NULL )
   {
      SCIP_CALL( processIntegerBoundChg(scip, cons, cutoff, nchgbds, &mustcheck) );
   }

   if( mustcheck && !(*cutoff) )
   {
      /* variable's fixings didn't give us any information -> we have to check the constraint */
      if( sol == NULL && consdata->row1 != NULL )
      {
         assert(consdata->row2 != NULL);
         /* skip constraints already in the LP */
         if( SCIProwIsInLP(consdata->row1) && SCIProwIsInLP(consdata->row2))
            return SCIP_OKAY;
         else
         {
            SCIP_Real feasibility;
           
            feasibility = 1.;
            
            assert(!SCIProwIsInLP(consdata->row1) || !SCIProwIsInLP(consdata->row2));

            /* check first row (linking) for feasibility */
            if( !SCIProwIsInLP(consdata->row1) )
            {
               SCIP_Real intsol;
               SCIP_Real binsol;
               int idx;

               intsol = SCIPgetVarSol(scip, consdata->intvar);
               idx = SCIPfeasFloor(scip, intsol) - consdata->offset;
               assert(idx < consdata->nbinvars);
               binsol = SCIPgetVarSol(scip, consdata->binvars[idx]);


               feasibility = MIN(feasibility, SCIPisFeasEQ(scip, intsol, (consdata->offset + idx) * binsol ) ? -1.0 : 1.0);
            }

            /* check second row (setppc) for feasibility */
            if( !SCIProwIsInLP(consdata->row2) )
               feasibility = MIN( feasibility, SCIPgetRowLPFeasibility(scip, consdata->row2) );
            
            addcut = SCIPisFeasNegative(scip, feasibility);
         }
      }
      else
         addcut = !checkCons(scip, cons, sol);
      
      
      if( !addcut )
      {
         /* constraint was feasible -> increase age */
         SCIP_CALL( SCIPincConsAge(scip, cons) );
      }
   }     
   
   if( addcut )
   {
      /* insert LP row as cut */
      assert(!(*cutoff));
      SCIP_CALL( addCuts(scip, cons, sol) );
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      *separated = TRUE;
   }
      
   return SCIP_OKAY;
}

/** enforces the pseudo solution on the given constraint */
static
SCIP_RETCODE enforcePseudo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint to be separated */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   SCIP_Bool*            infeasible,         /**< pointer to store TRUE, if the constraint was infeasible */
   int*                  nchgbds,            /**< pointer to store the number of changed variable bounds */
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
   assert(nchgbds != NULL);
   assert(solvelp != NULL);

   addcut = FALSE;
   mustcheck = TRUE;

   /* check constraint for violation only looking at the fixed variables, apply further fixings if possible */
   SCIP_CALL( processIntegerBoundChg(scip, cons, cutoff, nchgbds, &mustcheck) );
   SCIP_CALL( processBinvarFixings(scip, cons, cutoff, nchgbds, &addcut, &mustcheck) );

   if( mustcheck )
   {
      SCIP_CONSDATA* consdata;

      assert(!addcut);

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      if( checkCons(scip, cons, NULL) )
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
      assert(!(*cutoff));
      /* a cut must be added to the LP -> we have to solve the LP immediately */
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      *solvelp = TRUE;
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyLinking)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrLinking(scip) );
 
   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeLinking)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( conshdlrdataFree(scip, &conshdlrdata) );

   return SCIP_OKAY;
}

/** initialization method of constraint handler (called after problem was transformed) */
#define consInitLinking NULL

/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitLinking NULL

/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreLinking)
{
   SCIP_CONSDATA* consdata;
   int c;
   
   /* disable all linking constraints which at most one variable */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      
      if( consdata->nbinvars <= 1 )
      {
         SCIP_CALL( SCIPdisableCons(scip, conss[c]) );
      }
   }
   
   return SCIP_OKAY;
}

/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreLinking NULL

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolLinking NULL

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolLinking)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      
      /* release the rows of all constraints */
      if( consdata->row1 != NULL )
      {
         assert(consdata->row2 != NULL);

         SCIP_CALL( SCIPreleaseRow(scip, &consdata->row1) );
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->row2) );
      }

#if 0
      /* release binvars of all constraints */
      if( consdata->nbinvars  > 0 )
      {
         int b;
         
         /* drop bound change events of variables */
         assert(SCIPconsIsTransformed(conss[c]));
         SCIP_CALL( dropAllEvents(scip, consdata, conshdlrdata->eventhdlr) );
         
         SCIP_CALL( unlockRounding(scip, conss[c], consdata->binvars, consdata->nbinvars) );
         
         SCIPfreeBlockMemoryArrayNull(scip, &consdata->binvars, consdata->nbinvars);
      }
      consdata->nbinvars = 0;

      assert(consdata->nfixedzeros == 0);
      assert(consdata->nfixedones == 0);
#endif
   }

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteLinking)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* remove linking constraint form variable hash map */
   assert(SCIPhashmapExists(conshdlrdata->varmap, getHashmapKey((*consdata)->intvar)));
   SCIP_CALL( SCIPhashmapRemove(conshdlrdata->varmap, getHashmapKey((*consdata)->intvar)) );
      
   if( (*consdata)->nbinvars > 0 && SCIPisTransformed(scip) )
   {
      SCIP_CALL( dropAllEvents(scip, *consdata, conshdlrdata->eventhdlr) );
   }
   
   /* free consdata  */
   SCIP_CALL( consdataFree(scip, consdata) );
      
   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransLinking)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->row1 == NULL);  /* in original problem, there cannot be LP rows */
   assert(sourcedata->row2 == NULL);  /* in original problem, there cannot be LP rows */
   
   SCIPdebugMessage("transform linking constraint for variable <%s>\n", SCIPvarGetName(sourcedata->intvar));
   
   /* create constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, conshdlrdata->eventhdlr, &targetdata, 
         sourcedata->intvar, sourcedata->binvars, sourcedata->nbinvars, sourcedata->offset) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   /* insert (transformed) linking constraint into the hash map */
   SCIP_CALL( SCIPhashmapInsert(conshdlrdata->varmap, getHashmapKey(targetdata->intvar), *targetcons) );
   
   return SCIP_OKAY;
}

/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpLinking)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;
   
   for( c = 0; c < nconss; ++c )
   {
      assert(SCIPconsIsInitial(conss[c]));

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      
      if( consdata->nbinvars <= 1 )
         continue;
      
      SCIP_CALL( addCuts(scip, conss[c], NULL) );
   }
   
   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpLinking)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   SCIP_Bool separated;
   int nchgbds;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);
   
   SCIPdebugMessage("separating %d/%d linking constraints\n", nusefulconss, nconss);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   separated = FALSE;
   nchgbds = 0;
   
   /* check all useful linking constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], NULL, &cutoff, &separated, &nchgbds) );
   }
   
   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;
   else if( separated )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolLinking)
{  /*lint --e{715}*/
   SCIP_Bool cutoff;
   SCIP_Bool separated;
   int nchgbds;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   SCIPdebugMessage("separating %d/%d "CONSHDLR_NAME" constraints\n", nusefulconss, nconss);

   cutoff = FALSE;
   separated = FALSE;
   nchgbds = 0;

   /* check all useful set partitioning / packing / covering constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], sol, &cutoff, &separated, &nchgbds) );
   }

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;
   else if( separated )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpLinking)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   SCIP_Bool separated;
   int nchgbds;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   SCIPdebugMessage("LP enforcing %d linking constraints\n", nconss);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   separated = FALSE;
   nchgbds = 0;

   /* check all useful linking constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff && nchgbds == 0; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], NULL, &cutoff, &separated, &nchgbds) );
   }

   /* check all obsolete linking constraints for feasibility */
   for( c = nusefulconss; c < nconss && !cutoff && !separated && nchgbds == 0; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], NULL, &cutoff, &separated, &nchgbds) );
   }

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;
   else if( separated )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_FEASIBLE;
   
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsLinking)
{  /*lint --e{715}*/
   SCIP_Bool cutoff;
   SCIP_Bool infeasible;
   int nchgbds;
   SCIP_Bool solvelp;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   SCIPdebugMessage("pseudo enforcing %d "CONSHDLR_NAME" constraints\n", nconss);

   if( objinfeasible )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
   
   cutoff = FALSE;
   infeasible = FALSE;
   nchgbds = 0;
   solvelp = FALSE;

   /* check all linking constraint for domain reductions and feasibility */
   for( c = 0; c < nconss && !cutoff && !solvelp; ++c )
   {
      SCIP_CALL( enforcePseudo(scip, conss[c], &cutoff, &infeasible, &nchgbds, &solvelp) );
   }

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;
   else if( solvelp )
      *result = SCIP_SOLVELP;
   else if( infeasible )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckLinking)
{  /*lint --e{715}*/
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   int c;
   
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   /* check all linking constraints for feasibility */
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      if( consdata->nbinvars > 1 && (checklprows || consdata->row1 == NULL || !SCIProwIsInLP(consdata->row1)) )
      {
         if( !checkCons(scip, cons, sol) )
         {
            /* constraint is violated */
            *result = SCIP_INFEASIBLE;
            
            if( printreason )
            {
               int pos;
               int b;

               pos = -1;
               
#ifndef NDEBUG
               for( b = 0; b < consdata->nbinvars; ++b )
               {
                  assert( consdata->binvars[b] != NULL);
                  assert( SCIPvarGetType(consdata->binvars[b]) == SCIP_VARTYPE_BINARY );
               }
#endif

               SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
               
               /* check that at most one binary variable is fixed */
               for( b = 0; b < consdata->nbinvars; ++b )
               {
                  assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, consdata->binvars[b])) );
                  
                  /* check if binary variable is fixed */
                  if( SCIPgetSolVal(scip, sol, consdata->binvars[b]) > 0.5 )
                  {
                     if( pos != -1 )
                     {
                        SCIPinfoMessage(scip, NULL, "violation: more than one binary variable is set to one");
                        break;
                     }
                     pos = b ;
                  }
               }
               
               /* check that at least one binary variable is fixed */
               if( pos == -1 )
               {
                  SCIPinfoMessage(scip, NULL, "violation: none of the binary variables is set to one");
               }
               else if( !SCIPisFeasEQ(scip, pos + consdata->offset, SCIPgetSolVal(scip, sol, consdata->intvar)) )
               {
                  /* check if the fixed binary variable match with the integer variable */
                  SCIPinfoMessage(scip, NULL, "violation: <%s> = <%g> and <%s> is one\n",
                     SCIPvarGetName(consdata->intvar), SCIPgetSolVal(scip, sol, consdata->intvar),
                     SCIPvarGetName(consdata->binvars[pos]) );
               }
            }
            
            return SCIP_OKAY;
         }
      }
   }
   
   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropLinking)
{  
   SCIP_Bool cutoff;
   SCIP_Bool addcut;
   SCIP_Bool mustcheck;
   int nchgbds;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   SCIPdebugMessage("propagating %d/%d "CONSHDLR_NAME" constraints\n", nusefulconss, nconss);

   cutoff = FALSE;
   nchgbds = 0;
   addcut = FALSE;
   mustcheck = TRUE;

   /* propagate all useful set partitioning / packing / covering constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( processIntegerBoundChg(scip, conss[c], &cutoff, &nchgbds, &mustcheck) );
      SCIP_CALL( processBinvarFixings(scip, conss[c], &cutoff, &nchgbds, &addcut, &mustcheck) );
   }

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;
   
   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolLinking)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int oldnfixedvars;
   int oldnchgbds;
   int oldnaggrvars;
   int oldndelconss;
   int firstchange;
   int firstclique;
   int lastclique;
   int c;
   SCIP_Bool fixed;
   SCIP_Bool cutoff;
   SCIP_Bool infeasible;
   SCIP_Bool mustcheck;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("presolve %d linking constraints", nconss);

   (*result) = SCIP_DIDNOTFIND;

   oldnchgbds = *nchgbds;
   oldnaggrvars = *naggrvars;
   oldnfixedvars = *nfixedvars;
   oldndelconss = *ndelconss;
   cutoff = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* process constraints */
   firstchange = INT_MAX;
   firstclique = INT_MAX;
   lastclique = -1;

   /* check for each linking constraint the set partitioning condition */
   for( c = 0; c < nconss && !SCIPisStopped(scip); ++c )
   {
      SCIP_CONS* cons;
      SCIP_CONSDATA* consdata;

      assert(*result != SCIP_CUTOFF);

      cons = conss[c];
      assert(cons != NULL);
      assert(!SCIPconsIsModifiable(cons));

      SCIPdebugMessage("presolve linking constraints <%s>\n", SCIPconsGetName(cons));
 
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* in case there is only at most one binary variables, the constraints should already be disabled */
      if( consdata->nbinvars <= 1)
         continue;

      /*SCIPdebugMessage("presolving set partitioning / packing / covering constraint <%s>\n", SCIPconsGetName(cons));*/
      if( consdata->nfixedones >= 2 )
      {
         /* at least two variables are fixed to 1:
          * - a linkink constraint is infeasible due to the set partitioning condition
          */
         SCIPdebugMessage(""CONSHDLR_NAME" constraint <%s> is infeasible\n", SCIPconsGetName(cons));
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      
      if( consdata->nfixedones == 1 )
      {
         /* exactly one variable is fixed to 1:
          * - all other binary variables must be zero due to the set partitioning condition
          * - integer variable has to be fixed to corresponding binary variable which is fixed to one
          * - if constraint is not modifiable it can be removed 
          */
         SCIP_VAR* var;
         int v;

         SCIPdebugMessage(""CONSHDLR_NAME" constraint <%s> has a binary variable fixed to 1.0\n", SCIPconsGetName(cons));

         for( v = 0; v < consdata->nbinvars; ++v )
         {
            var = consdata->binvars[v];
            if( SCIPvarGetLbGlobal(var) < 0.5 && SCIPvarGetUbGlobal(var) > 0.5 )
            {
               SCIP_CALL( SCIPfixVar(scip, var, 0.0, &infeasible, &fixed) );
               if( infeasible )
               {
                  SCIPdebugMessage(""CONSHDLR_NAME" constraint <%s>: infeasible fixing <%s> == 0\n",
                     SCIPconsGetName(cons), SCIPvarGetName(var));
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }
               assert(fixed);
               (*nfixedvars)++;
            }
            else if( SCIPvarGetLbGlobal(var) > 0.5 )
            {
               /* fix integer variable */
               SCIP_CALL( SCIPfixVar(scip, consdata->intvar, (SCIP_Real)(v+consdata->offset), &infeasible, &fixed) );
               if( infeasible )
               {
                  SCIPdebugMessage(""CONSHDLR_NAME" constraint <%s>: infeasible fixing <%s> == %d\n",
                     SCIPconsGetName(cons), SCIPvarGetName(consdata->intvar), v+consdata->offset);
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }
               
               if( fixed )
                  (*nfixedvars)++;
            }
         }

	 /* now all other variables are fixed to zero:
          * the constraint is feasible, and if it's not modifiable, it is redundant
          */
         SCIPdebugMessage(""CONSHDLR_NAME" constraint <%s> is redundant\n", SCIPconsGetName(cons));
         SCIP_CALL( SCIPdelCons(scip, cons) );
         (*ndelconss)++;
         continue;
      }
      
      if( consdata->nfixedzeros == consdata->nbinvars )
      {
         /* all variables are fixed to zero:
          * - a linking constraint is infeasible due the set partitioning condition
          */
         assert(consdata->nfixedones == 0);

         SCIPdebugMessage("linking constraint <%s> is infeasible due to set partitioning condition\n", SCIPconsGetName(cons));
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      
      if( consdata->nfixedzeros == consdata->nbinvars - 1 )
      {
         /* all variables except one are fixed to zero:
          * - a linking constraint is feasible due the set partitioning condition
          * - the remaining binary variable can be fixed to one
          * - integer variable has to be fixed to corresponding binary variable which  is fixed  to one
          * - constraint can be deleted since it is not modifiable
          */
         SCIP_VAR* var;
         SCIP_Bool found;
         int v;

         assert(consdata->nfixedones == 0);

         SCIPdebugMessage(""CONSHDLR_NAME" constraint <%s> has only one binary variable not fixed to zero\n",
            SCIPconsGetName(cons));
            
         /* search unfixed variable */
         found = FALSE;
         var = NULL;
         for( v = 0; v < consdata->nbinvars && !found; ++v )
         {
            var = consdata->binvars[v];
            found = SCIPvarGetUbGlobal(var) > 0.5;
         }
         assert(found);

         /*fix remaining binary variable */
         SCIP_CALL( SCIPfixVar(scip, var, 1.0, &infeasible, &fixed) );
         if( infeasible )
         {
            SCIPdebugMessage(""CONSHDLR_NAME" constraint <%s>: infeasible fixing <%s> == 1\n",
               SCIPconsGetName(cons), SCIPvarGetName(var));
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         assert(fixed);
         (*nfixedvars)++;

         /* fix integer variable */
         SCIP_CALL( SCIPfixVar(scip, consdata->intvar, (SCIP_Real)(v+consdata->offset), &infeasible, &fixed) );
         if( infeasible )
         {
            SCIPdebugMessage(""CONSHDLR_NAME" constraint <%s>: infeasible fixing <%s> == %d\n",
               SCIPconsGetName(cons), SCIPvarGetName(consdata->intvar), v+consdata->offset);
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         assert(fixed);
         (*nfixedvars)++;
            
         /* delete constraint from  problem */
         SCIP_CALL( SCIPdelCons(scip, cons) );
         (*ndelconss)++;
         continue;
      }
      
      if( consdata->nfixedzeros == consdata->nbinvars - 2 ) /*lint !e641*/
      {
         SCIP_VAR* var;
         SCIP_VAR* var1;
         SCIP_VAR* var2;
         SCIP_Bool redundant;
         SCIP_Bool aggregated;
         int v;

         /* aggregate variable, if set partitioning condition consists only of two
          * non-fixed variables
          */

         /* search unfixed variable */
         var1 = NULL;
         var2 = NULL;
         for( v = 0; v < consdata->nbinvars && var2 == NULL; ++v )
         {
            var = consdata->binvars[v];
            if( SCIPvarGetUbGlobal(var) > 0.5 ) 
            {
               if( var1 == NULL )
                  var1 = var;
               else
                  var2 = var;
            }
         }
         assert(var1 != NULL && var2 != NULL);

         /* aggregate binary equality var1 + var2 == 1 */
         SCIPdebugMessage(""CONSHDLR_NAME" constraint <%s>: aggregate <%s> + <%s> == 1\n",
            SCIPconsGetName(cons), SCIPvarGetName(var1), SCIPvarGetName(var2));
         SCIP_CALL( SCIPaggregateVars(scip, var1, var2, 1.0, 1.0, 1.0, &infeasible, &redundant, &aggregated) );

         /* evaluate aggregation result */
         if( infeasible )
         {
            SCIPdebugMessage("linking constraint <%s>: infeasible aggregation <%s> + <%s> == 1\n",
               SCIPconsGetName(cons), SCIPvarGetName(var1), SCIPvarGetName(var2));
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         if( aggregated )
            (*naggrvars)++;
      }

      /* apply integer bound to binary variables */
      SCIP_CALL( processIntegerBoundChg(scip, cons, &cutoff, nchgbds, &mustcheck) );
      
      /* tightened integer variable */
      SCIP_CALL( tightenedIntvar(scip, cons, consdata, &cutoff, nchgbds, TRUE) );
            
      if( cutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      /* remember the first changed constraint to begin the next redundancy round with */
      if( firstchange == INT_MAX )
         firstchange = c;

      /* remember the first and last constraints for which we have to add the clique information */
      if( !consdata->cliqueadded && consdata->nbinvars >= 2 )
      {
         if( firstclique == INT_MAX )
            firstclique = c;
         lastclique = c;
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

      if( !consdata->cliqueadded && consdata->nbinvars >= 3 )
      {
         /* add set partitioning condition as clique */
         int ncliquebdchgs;

         SCIP_CALL( SCIPaddClique(scip, consdata->binvars, NULL, consdata->nbinvars, &infeasible, &ncliquebdchgs) );
         *nchgbds += ncliquebdchgs;

         if( infeasible )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         consdata->cliqueadded = TRUE;
      }
   }

   /* transfer aggregated integer variables to the corresponding binary variables */
   SCIP_CALL( aggregateVariables(scip, conshdlrdata->varmap, conss, nconss, naggrvars, &cutoff) );

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( oldndelconss < *ndelconss || oldnfixedvars < *nfixedvars || oldnchgbds < *nchgbds || oldnaggrvars < *naggrvars)
      *result = SCIP_SUCCESS;
   
   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropLinking)
{  /*lint --e{715}*/
   
   SCIP_CONSDATA* consdata;
   SCIP_VAR* intvar;
   int v;

   SCIPdebugMessage("conflict resolving method of "CONSHDLR_NAME" constraint handler\n");

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   intvar = consdata->intvar;
   assert( intvar != NULL);
   
   *result = SCIP_DIDNOTFIND;

   if( inferinfo == -1 )
   {
      /* we have to resolve a fixing of a binary variable which was done due to fixed binary variables */
      assert(SCIPvarGetType(infervar) == SCIP_VARTYPE_BINARY); 
      assert(SCIPisFeasEQ(scip, SCIPvarGetUbAtIndex(intvar, bdchgidx, TRUE), SCIPvarGetUbAtIndex(intvar, bdchgidx, FALSE)));
      assert(SCIPisFeasEQ(scip, SCIPvarGetLbAtIndex(intvar, bdchgidx, TRUE), SCIPvarGetLbAtIndex(intvar, bdchgidx, FALSE)));

      if( boundtype == SCIP_BOUNDTYPE_UPPER )
      {
         /* we fixed the binary variable to zero since one of the others was fixed to one */
         assert(SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE) < 0.5); 
         
         for( v = 0; v < consdata->nbinvars; ++v )
         {
            if( SCIPvarGetLbAtIndex(consdata->binvars[v], bdchgidx, FALSE) > 0.5 )
            {
               SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->binvars[v]) );
               break;
            }
         }
         assert(v < consdata->nbinvars); 
      }
      else
      {
         /* we fixed the binary variable to one since all other binary variable were fixed to zero */
         assert(boundtype == SCIP_BOUNDTYPE_LOWER);
         assert(SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE) > 0.5); 
         
         for( v = 0; v < consdata->nbinvars; ++v )
         {
            if( consdata->binvars[v] != infervar )
            {
               /* the reason variable must be assigned to zero */
               assert(SCIPvarGetUbAtIndex(consdata->binvars[v], bdchgidx, FALSE) < 0.5);
               SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->binvars[v]) );
            }
         }
      }
   }
   else if( inferinfo == -2 )
   {
      SCIP_Real lb;
      
      /* we have to resolve a fixing of a binary variable which was done due to the integer variable lower bound */
      assert(SCIPvarGetType(infervar) == SCIP_VARTYPE_BINARY); 
      assert(SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE) < 0.5); 
      assert(SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE) < 0.5); /*@repair: neu*/
      assert(SCIPvarGetUbAtIndex(infervar, bdchgidx, FALSE) > 0.5); /*@repair: neu*/
      assert( SCIPisFeasEQ(scip, SCIPvarGetUbAtIndex(intvar, bdchgidx, TRUE), SCIPvarGetUbAtIndex(intvar, bdchgidx, FALSE)) );
      assert( SCIPisFeasEQ(scip, SCIPvarGetLbAtIndex(intvar, bdchgidx, TRUE), SCIPvarGetLbAtIndex(intvar, bdchgidx, FALSE)) );    

      lb = SCIPvarGetLbAtIndex(intvar, bdchgidx, TRUE);
      assert(infervar != consdata->binvars[(int)(lb + 0.5) - consdata->offset] ); /*@repaired: vorher + offset && == */

      SCIP_CALL( SCIPaddConflictLb( scip, intvar, bdchgidx) );
   }
   else if( inferinfo == -3 )
   {
      SCIP_Real ub;
      
      /* we have to resolve a fixing of a binary variable which was done due to the integer variable upper bound */
      assert(SCIPvarGetType(infervar) == SCIP_VARTYPE_BINARY); 
      assert(SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE) < 0.5); 
      assert(SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE) < 0.5); 
      assert(SCIPvarGetUbAtIndex(infervar, bdchgidx, FALSE) > 0.5); 
      assert( SCIPisFeasEQ(scip, SCIPvarGetUbAtIndex(intvar, bdchgidx, TRUE), SCIPvarGetUbAtIndex(intvar, bdchgidx, FALSE)) );
      assert( SCIPisFeasEQ(scip, SCIPvarGetLbAtIndex(intvar, bdchgidx, TRUE), SCIPvarGetLbAtIndex(intvar, bdchgidx, FALSE)) );

      ub = SCIPvarGetUbAtIndex(intvar, bdchgidx, TRUE);
      assert(infervar != consdata->binvars[(int)(ub + 0.5) - consdata->offset] );/*@repaired: vorher +offset */
      
      SCIP_CALL( SCIPaddConflictUb( scip, intvar, bdchgidx) );
   }
   else if( inferinfo == -4 )
   {
      int oldlb;
      int newlb;
      int offset;
      int b;
      
      /* we tightened the lower bound of the integer variable due the fixing of the corresponding binary variable to zero */
      
      assert(infervar == intvar);
      assert(boundtype == SCIP_BOUNDTYPE_LOWER);
      
      /* get old and new lower bound */
      oldlb = (int)(SCIPvarGetLbAtIndex(intvar, bdchgidx, FALSE) + 0.5); /*@repair: were both lower bounds */
      newlb = (int)(SCIPvarGetLbAtIndex(intvar, bdchgidx, TRUE) + 0.5);
      assert(oldlb < newlb);

      /* add old lower bound of integer variable to conflict */
      SCIP_CALL( SCIPaddConflictLb( scip, intvar, bdchgidx) ); /* @repair: new, belongs to the conflict  */

      offset = consdata->offset;

      /* add binvars to conflict */
      for( b = oldlb-offset; b < newlb-offset; ++b ) 
      {
         assert(b  >= 0);
         assert(b < consdata->nbinvars);
         assert(SCIPvarGetUbLocal(consdata->binvars[b]) < 0.5);

         SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->binvars[b]) );
      }
   }
   else if( inferinfo == -5 )
   {
      int oldub;
      int newub;
      int offset;
      int b;
      
      /* we tightened the upper bound of the integer variable due the fixing of the corresponding binary variable two zero */
      
      assert(infervar == intvar);
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);

      /* get old and new upper bound */
      oldub = (int)(SCIPvarGetUbAtIndex(intvar, bdchgidx, FALSE) + 0.5);
      newub = (int)(SCIPvarGetUbAtIndex(intvar, bdchgidx, TRUE) + 0.5);
      assert(oldub > newub);

      /* add old upper bound of integer variable to conflict */
      SCIP_CALL( SCIPaddConflictUb( scip, intvar, bdchgidx) ); 
      
      offset = consdata->offset;

      /* resolve tightening of upper bound of the integer variable by binary variables */
      for( b = oldub-offset; b > newub-offset; --b ) 
      {
         assert(b  >= 0);
         assert(b < consdata->nbinvars);
         assert(SCIPvarGetUbLocal(consdata->binvars[b]) < 0.5 || b+offset == newub); /*@repair: zweite bedingung neu*/

         SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->binvars[b]) );
      }
   }
   else if( inferinfo == -6 )
   {
      /* we fixed a binary variable to one since the integer variable was fixed */
      assert(SCIPvarGetType(infervar) == SCIP_VARTYPE_BINARY); 
      assert(boundtype == SCIP_BOUNDTYPE_LOWER);
      assert( SCIPisFeasEQ(scip, SCIPvarGetUbAtIndex(intvar, bdchgidx, TRUE), SCIPvarGetUbAtIndex(intvar, bdchgidx, FALSE)) );
      assert( SCIPisFeasEQ(scip, SCIPvarGetLbAtIndex(intvar, bdchgidx, TRUE), SCIPvarGetUbAtIndex(intvar, bdchgidx, FALSE)) );
      assert( SCIPisFeasEQ(scip, SCIPvarGetUbAtIndex(intvar, bdchgidx, TRUE), SCIPvarGetLbAtIndex(intvar, bdchgidx, FALSE)) );
      assert( SCIPisFeasEQ(scip, SCIPvarGetLbAtIndex(intvar, bdchgidx, TRUE), SCIPvarGetLbAtIndex(intvar, bdchgidx, FALSE)) );
      
      assert( !SCIPisFeasEQ(scip, SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE), SCIPvarGetLbAtIndex(infervar, bdchgidx, FALSE))  );
      
      SCIP_CALL( SCIPaddConflictLb( scip, intvar, bdchgidx) );
      SCIP_CALL( SCIPaddConflictUb( scip, intvar, bdchgidx) );
   }
   else
   {
      /* we fixed the integer variable to (inferinfo + offset) since the corresponding binary variable was fixed to one */
      assert(infervar == intvar);
      assert(inferinfo >= 0);
      assert(inferinfo < consdata->nbinvars);
      assert(inferinfo + consdata->offset == (int)(SCIPvarGetUbAtIndex(consdata->intvar, bdchgidx, TRUE) + 0.5)
         || inferinfo + consdata->offset == (int)(SCIPvarGetLbAtIndex(consdata->intvar, bdchgidx, TRUE) + 0.5));
      /* @repair: possibly only one bound has changed at this point in time, not both */
      
      assert(SCIPvarGetLbAtIndex(consdata->binvars[inferinfo], bdchgidx, FALSE) > 0.5);
      SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->binvars[inferinfo]) );
   }

   *result = SCIP_SUCCESS;
   
   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockLinking)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int b;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   /* look integer variable in both directions */
   SCIP_CALL( SCIPaddVarLocks(scip, consdata->intvar, nlockspos + nlocksneg, nlockspos + nlocksneg) );

   /* look binary variables in both directions */
   for( b = 0; b < consdata->nbinvars; ++b )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->binvars[b], nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }
   
   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#define consActiveLinking NULL


/** constraint deactivation notification method of constraint handler */
#define consDeactiveLinking NULL


/** constraint enabling notification method of constraint handler */
#define consEnableLinking NULL


/** constraint disabling notification method of constraint handler */
#define consDisableLinking NULL


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintLinking)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );

   consdataPrint(scip, SCIPconsGetData(cons), file);

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyLinking)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourceconsdata;
   SCIP_VAR** binvars;
   SCIP_VAR* intvar;
   int offset;
   int nbinvars;
   int size;
   
   SCIP_Real coef;
   SCIP_Real constant;
   int requiredsize;
   int v;
   
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(sourcecons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a linking constraint\n");
      SCIPABORT();
   }

   (*success) = TRUE;
   
   sourceconsdata = SCIPconsGetData(sourcecons);
   assert(sourceconsdata != NULL);

   /* get number of binary variables, intger variables, and offset of the source constraint */
   nbinvars = sourceconsdata->nbinvars;
   intvar = sourceconsdata->intvar;
   offset = sourceconsdata->offset;
   
   /* duplicate variable array */
   if( nbinvars > 0 )
   {
      SCIP_CALL( SCIPduplicateBufferArray(scip, &binvars, sourceconsdata->binvars, nbinvars) );
   }   
   else
      binvars = NULL;
   
   coef = 1.0;
   constant = 0.0;
   size = 1;
   
   if( SCIPisTransformed(sourcescip) )
   {
      SCIP_Bool negated;

      /* convert the binary variable into active variables and map this active variable to counter part in the target
       * SCIP */
      for( v = 0; v < nbinvars && *success; ++v )
      {
         SCIP_CALL( SCIPgetBinvarRepresentative(sourcescip, binvars[v], &binvars[v], &negated)  );
         assert(!negated);
         
         binvars[v] = (SCIP_VAR*) (size_t) SCIPhashmapGetImage(varmap, binvars[v]);
         
         if( binvars[v] == NULL )
            (*success) = FALSE;
      }
      
      if( *success )
      {
         /* convert the integer variable into an active variable and map this variable to the counter part in the target SCIP */
         SCIP_CALL( SCIPgetProbvarLinearSum(sourcescip, &intvar, &coef, &size, size, &constant, &requiredsize, TRUE) );
      
         if(requiredsize == 1 && SCIPisEQ(scip, coef, 1.0) && SCIPisZero(scip, constant) )
         {
            assert(SCIPisZero(scip, constant));
            intvar = (SCIP_VAR*) (size_t) SCIPhashmapGetImage(varmap, intvar);

            if( intvar == NULL )
               (*success) = FALSE;
         }
         else
            (*success) = FALSE;
      }
   }
   else
   {
      /* convert the binary variable into active variables and map this active variable to counter part in the target SCIP */
      for( v = 0; v < nbinvars; ++v )
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&binvars[v], &coef, &constant) );

         if( !SCIPisEQ(scip, coef, 1.0) || !SCIPisZero(scip, constant) ) 
         {
            (*success) = FALSE;
            break;
         }

         binvars[v] = (SCIP_VAR*) (size_t) SCIPhashmapGetImage(varmap, binvars[v]);
         
         if( binvars[v] == NULL )
         {
            (*success) = FALSE;
            break;
         }
      }
      
      if( *success )
      {
         /* convert the integer variable into an active variable and map this variable to the counter part in the target SCIP */
         SCIP_CALL( SCIPvarGetOrigvarSum(&intvar, &coef, &constant) );

         if( !SCIPisEQ(scip, coef, 1.0) || !SCIPisZero(scip, constant) ) 
            (*success) = FALSE;
         else
         {
            binvars[v] = (SCIP_VAR*) (size_t) SCIPhashmapGetImage(varmap, binvars[v]);
            
            if( binvars[v] == NULL )
               (*success) = FALSE;
         }
      }
   }
   
   if( *success )
   {      
      const char* consname;

      if( name != NULL )
         consname = name;
      else
         consname = SCIPconsGetName(sourcecons);

      SCIP_CALL( SCIPcreateConsLinking(scip, cons, name, intvar, binvars, nbinvars, offset,  
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }
   
   /* free buffer array */
   if( nbinvars > 0 )
   {
      SCIPfreeBufferArrayNull(scip, &binvars);
   }
   
   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
#define consParseLinking NULL






// static
// SCIP_DECL_LINCONSUPGD(linconsUpgdLinking)
// {  /*lint --e{715}*/

   
//    SCIP_VAR** binvars; 
//    SCIP_VAR* intvar; 
//    int* coeffs;

//    SCIP_Bool islinking;
//    int nbinvars; 
//    int offset; 
//    int n; 
//    int i;

//    //printf("try upgrade to linking constraint\n");

//    n = 0;
//    /*@todo: there must also be a set partitioning constraint! */

//    assert(upgdcons != NULL);

//    //   *result = SCIP_DIDNOTFIND;

//    /* - lhs and rhs are equal 
//     * - no fractional coefficients
//     * - 
//     * 
//     * 
//     */
//    if( SCIPisEQ(scip, lhs, rhs) && SCIPisFeasIntegral(scip, lhs) && ncoeffspfrac + ncoeffsnfrac == 0 )
//    {
//       // printf("equality constraint\n");

//       if( !(nposint == 1 && nnegbin + nposint == nvars) && !(nnegint == 1 && nnegint + nposbin == nvars) )
//       {
//          return SCIP_OKAY;
//       }
   
//       for( i = 0; i < nvars; ++i )
//       {
//          if( !SCIPisFeasIntegral(scip, vals[i]) )
//             return SCIP_OKAY;
//       }

//       /* integer variable coefficient must be 1 */
//       if( !SCIPisFeasEQ(scip, vals [nvars-1], 1.) )
//          return SCIP_OKAY;


//       printf("one integer variable and binary variables with integral coefficients\n");   
   
//       /* set offset to rhs, has to be updated eventually */
//       offset = (int)(lhs + 0.5);
      
//       nbinvars = nposbin + nnegbin;
//       SCIP_CALL( SCIPallocBufferArray(scip, &binvars, nbinvars) );
//       SCIP_CALL( SCIPallocBufferArray(scip, &coeffs, nbinvars) ); 
      
//       /* add binvars with positiv coefficients to array */
//       for( i = 0; i < nposbin; ++i )
//       {
//          coeffs[n] = vals[i];
//          binvars[n] = vars[i];
//          n++;
//       }
//       /* add binvars with negativ coefficients to array */ 
//       for( i = 0; i < nposbin; ++i ) 
//       { 
//          coeffs[n] = - vals[i];
//          binvars[n] = vars[i]; 
//          n++; 
//       } 

//       assert(n == nposbin + nnegbin);
//       assert(n == nvars - 1);
      

//       /* sort coefficients */
//       SCIPsortIntPtr(coeffs, (void**)binvars, nbinvars);

//       /* update offset with smallest coefficient */
//       offset += coeffs[0];

//       islinking = TRUE;

//       /* check coefficients to be increasing with step length 1 */      
//       for( i = 1; i < nbinvars; ++i )
//       {
//          if( coeffs[i] != coeffs[i-1] )
//          {
//             islinking = FALSE;
//             break;
//          }
//       }
      

//       intvar = vars[nvars-1];
//       assert(SCIPvarGetType(intvar) == SCIP_VARTYPE_INTEGER);

//       if( islinking )
//       {

//          SCIPdebugMessage("upgrading constraint <%s> to linking constraint\n", SCIPconsGetName(cons)); 

//          assert(!SCIPconsIsModifiable(cons)); 
//          SCIP_CALL( SCIPcreateConsLinking(scip, upgdcons, SCIPconsGetName(cons), 
//                intvar, binvars, nbinvars, offset,   
//                SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
//                SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), 
//                SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), 
//                SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
//          //*result = SCIP_SUCCESS; 


//          SCIPinfoMessage(scip, NULL, "upgrading constraint <%s> to linking constraint\n", SCIPconsGetName(cons)); 
//       }

//    }

//    return SCIP_OKAY;

// }







/*
 * Callback methods of event handler
 */

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecBinvar)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

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
   default:
      SCIPerrorMessage("invalid event type\n");
      return SCIP_INVALIDDATA;
   }
   assert(0 <= consdata->nfixedzeros && consdata->nfixedzeros <= consdata->nbinvars);
   assert(0 <= consdata->nfixedones && consdata->nfixedones <= consdata->nbinvars);

   /*debugMessage(" -> constraint has %d zero-fixed and %d one-fixed of %d variables\n",
     consdata->nfixedzeros, consdata->nfixedones, consdata->nvars);*/

   return SCIP_OKAY;
}

/*
 * constraint specific interface methods
 */

/** creates the handler for linking constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrLinking(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, eventExecBinvar, NULL) );
   
   /* create linking constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         conshdlrCopyLinking,
         consFreeLinking, consInitLinking, consExitLinking,
         consInitpreLinking, consExitpreLinking, consInitsolLinking, consExitsolLinking,
         consDeleteLinking, consTransLinking, consInitlpLinking,
         consSepalpLinking, consSepasolLinking, consEnfolpLinking, consEnfopsLinking, consCheckLinking,
         consPropLinking, consPresolLinking, consRespropLinking, consLockLinking,
         consActiveLinking, consDeactiveLinking,
         consEnableLinking, consDisableLinking,
         consPrintLinking, consCopyLinking, consParseLinking,
         conshdlrdata) );


   /* include the linear constraint to linking constraint upgrade in the linear constraint handler */
   //  SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdLinking, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );

   
   /* add linking constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, 
         "constraints/"CONSHDLR_NAME"/linearize", "this constraint will not propagate or separate, linear and setppc are used?", 
         &conshdlrdata->linearize, FALSE, DEFAULT_LINEARIZE, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures an linking constraint */
SCIP_RETCODE SCIPcreateConsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             intvar,             /**< integer variable which should be linked */
   SCIP_VAR**            binvars,            /**< binary variables, or NULL */
   int                   nbinvars,           /**< number of binary variables */
   int                   offset,             /**< offset ot the binary variable representation */
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
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(!SCIPisInfinity(scip, -SCIPvarGetLbGlobal(intvar)));
   assert(!SCIPisInfinity(scip, SCIPvarGetUbGlobal(intvar)));

   /* find the linking constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("linking constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   SCIPdebugMessage("create linking constraint for variable <%s> (SCIP stage %d)\n", 
      SCIPvarGetName(intvar), SCIPgetStage(scip));

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* check if the linking for the requestes integer variable already exists */
   assert(!SCIPhashmapExists(conshdlrdata->varmap, getHashmapKey(intvar)));

   /* create the constraint specific data */
   SCIP_CALL( consdataCreate(scip, conshdlrdata->eventhdlr, &consdata, intvar, binvars, nbinvars, offset) );
   
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, 
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

   if( SCIPisTransformed(scip) && consdata->nbinvars <= 1 )
   {
      SCIP_CALL( SCIPdisableCons(scip, *cons) );
   }
   
   /* insert linking constraint into the hash map */
   SCIP_CALL( SCIPhashmapInsert(conshdlrdata->varmap, getHashmapKey(intvar), *cons) );
   assert(SCIPhashmapExists(conshdlrdata->varmap, getHashmapKey(intvar)));

   return SCIP_OKAY;
}


/** checks if for the given integer variable a linking constraint exists */
SCIP_Bool SCIPexistsConsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             intvar              /**< integer variable which should be linked */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   
   return SCIPhashmapExists(conshdlrdata->varmap, getHashmapKey(intvar));
}
 
/** returns the linking constraint belonging the given integer variable or NULL if it does not exist yet */
SCIP_CONS* SCIPgetConsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             intvar              /**< integer variable which should be linked */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   
   return (SCIP_CONS*) SCIPhashmapGetImage(conshdlrdata->varmap, getHashmapKey(intvar));
}
 
/** returns the integer variable of the linking constraint */
SCIP_VAR* SCIPgetIntvarLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linking constraint */
   )
{
   SCIP_CONSDATA* consdata;
   
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a "CONSHDLR_NAME" constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->intvar;

}

/** returns the integer variable if this binary variable is linked to an INTEGER variable, or NULL, inefficient method */
SCIP_VAR* SCIPgetIntvarFromBinvarLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             binvar              /**< binary variable */
   )
{
   SCIP_CONS** conss;

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   
   int nconss;
   int c;
   int b;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   /* get constraints */
   nconss = SCIPconshdlrGetNConss(conshdlr);
   conss = SCIPconshdlrGetConss(conshdlr);
   assert(nconss == 0 || conss != NULL);

   /* check each linking constraint */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* compae each binvar to the given one */
      for( b = 0; b < consdata->nbinvars; ++b )
      {
         if( consdata->binvars[b] == binvar )
         {
            /* found integer variable */
            return consdata->intvar;
         }
      }
   }
   
   /* this binvar is not in any linking constraint */
   return NULL;
}

/** returns the binary variables of the linking constraint */
SCIP_RETCODE SCIPgetBinvarsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint */
   SCIP_VAR***           binvars,            /**< pointer to store the binary variables array pointer */
   int*                  nbinvars            /**< pointer to store the number of returned binary variables */
   )
{
   SCIP_CONSDATA* consdata;
   
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a "CONSHDLR_NAME" constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->binvars == NULL )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;
      
      conshdlr = SCIPconsGetHdlr(cons);
      assert(conshdlr != NULL);

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      SCIP_CALL( consdataCreateBinvars(scip, cons, consdata, conshdlrdata->eventhdlr, conshdlrdata->linearize) );
   }

   assert(consdata->binvars != NULL);
   
   if( binvars != NULL )
      (*binvars) = consdata->binvars;
   if( nbinvars != NULL )
      (*nbinvars) = consdata->nbinvars;

   return SCIP_OKAY;
}

/** returns the number of binary variables of the linking constraint */
int SCIPgetNBinvarsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linking constraint */
   )
{
   SCIP_CONSDATA* consdata;
   
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a "CONSHDLR_NAME" constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nbinvars;
}

/** returns the offset of the linking constraint */
int SCIPgetOffsetLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linking constraint */
   )
{
   SCIP_CONSDATA* consdata;
   
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a "CONSHDLR_NAME" constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->offset;
}


/** returns whether a given solution satisfies intvar = \sum_t  (t * x_{jt}) && \exist t: x_{jt} = 1 */
SCIP_Bool SCIPcheckConsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint to be checked */
   SCIP_SOL*             sol                 /**< primal solution, or NULL for current LP/pseudo solution */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** binvars;
   SCIP_Real solval;
   
   int nbinvars;
   int pos;
   int offset;
   
   
   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a "CONSHDLR_NAME" constraint\n");
      SCIPABORT();
   }

   SCIPdebugMessage("checking precedence constraint <%s> for feasibility of solution %p\n", SCIPconsGetName(cons), (void*)sol);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->binvars != NULL || consdata->nbinvars == 0);
   
   solval = SCIPgetSolVal(scip, sol, consdata->intvar);

   /* return false for fractional start variables */
   if( !SCIPisFeasIntegral(scip, solval) )
   {
      // printf("solval is fractional\n");
      return FALSE;
   }
    
   binvars = consdata->binvars;
   nbinvars = consdata->nbinvars;
   offset = consdata->offset;

   if( consdata->nbinvars == 0 || consdata->nbinvars == 1 )
      return TRUE;

   pos = (int)(solval+0.5)-offset;
   assert(pos >= 0 && pos < nbinvars);

   /* the binvar for time point 'solval' has to be set to one */
   if( SCIPisFeasEQ(scip, SCIPgetSolVal(scip, sol,binvars[pos]), 1.) )
   {
      return TRUE;
   }
  
   return FALSE;
}


/** returns whether this constraint is able to branch or not */
SCIP_Bool SCIPisLPCandLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons               /**< linking constraint to be checked */
   )
{

   assert(scip != NULL);
   assert(cons != NULL);
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a "CONSHDLR_NAME" constraint\n");
      SCIPABORT();
   }

   return !SCIPcheckConsLinking(scip, cons, NULL);
}



/** returns value for fractionality of this constraint in an arbitrary solution (current LP or pseudo) */
SCIP_Real SCIPgetFractionalityLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons               /**< linking constraint to be checked */
   )
{
   SCIP_VAR* intvar;
   SCIP_CONSDATA* consdata;

   SCIP_Real solvalue;
   SCIP_Real frac;

   assert(scip != NULL);
   assert(cons != NULL);
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a "CONSHDLR_NAME" constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   frac = 0.;
   intvar = consdata->intvar;

   solvalue = SCIPgetVarSol(scip, intvar);

   frac = (SCIPvarGetUbLocal(intvar) - solvalue) / (SCIPvarGetUbLocal(intvar) - SCIPvarGetLbLocal(intvar));

   return frac;
}

/** returns array of constraints that are able to perform branching and probing */
SCIP_RETCODE SCIPgetLPCandsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< linking constraint handler */
   SCIP_CONS**           conscands,          /**< array to store the constraints that are able to branch */
   SCIP_Real*            conscandsfracs,     /**< array with fractionalities for each constraint */
   int*                  nconscands          /**< pointer to store number of constraint candidates */
   )
{

   SCIP_CONS** conss;
   int nconss;
   int noldconscands;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conscands != NULL);
   assert(nconscands != NULL);

   noldconscands = *nconscands;

   conss = SCIPconshdlrGetConss(conshdlr);
   nconss = SCIPconshdlrGetNConss(conshdlr);
   assert(conss != NULL);

   /* check each constraint whether it is a branching candidate */
   for( c = 0; c < nconss; ++c )
   {
      if( SCIPisLPCandLinking(scip, conss[c]) )
      {
         /* add constraint to candidate store */
         conscands[*nconscands] = conss[c];
         conscandsfracs[*nconscands] = SCIPgetFractionalityLinking(scip, conss[c]);
         ++(*nconscands);
      }
   }

   SCIPdebugMessage("added %d branching linking constraint candidates\n",*nconscands-noldconscands);
   
   return SCIP_OKAY;
}




// /** adds given variable and bound change to hash map and bound change arrays */
// static
// SCIP_RETCODE addBdchg(
//    SCIP*                 scip,               /**< SCIP data structure */
//    BDCHGDATA*            bdchgdata,          /**< structure to keep bound chage data */
//    SCIP_VAR*             var,                /**< variable to store bound change */
//    int                   newbound,           /**< new bound for given variable */
//    SCIP_BOUNDTYPE        boundtype,          /**< lower or upper bound change */
//    int*                  nbdchgs,            /**< total number of bound changes occured so far */
//    SCIP_Bool*            infeasible          /**< pointer to store whether bound change makes the node infeasible */
//    )
// {
//    int nvars; 
//    int pos;

//    assert(scip != NULL);
//    assert(bdchgdata != NULL);
//    assert(bdchgdata->varhashmap != NULL);
//    assert(bdchgdata->lbchgs != NULL);
//    assert(bdchgdata->ubchgs != NULL);
//    assert(var != NULL);
   
//    nvars = bdchgdata->nvars;

//    if( infeasible != NULL )
//       *infeasible = FALSE;

//    if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
//       return SCIP_OKAY;

//    /* if variable is not in hash map insert it and increase array sizes */
//    if( !SCIPhashmapExists(bdchgdata->varhashmap, var) )
//    {
//       printf("capture additional variable <%s>\n", SCIPvarGetName(var) );
//       //SCIP_CALL( SCIPcaptureVar(scip, var) );
//       SCIPhashmapInsert(bdchgdata->varhashmap, var, (void*) (size_t)nvars);
//       SCIP_CALL( SCIPreallocBufferArray(scip, &bdchgdata->lbchgs, nvars + 1) );
//       SCIP_CALL( SCIPreallocBufferArray(scip, &bdchgdata->ubchgs, nvars + 1) );

//       bdchgdata->lbchgs[nvars] = SCIPfeasCeil(scip, SCIPvarGetLbLocal(var));
//       bdchgdata->ubchgs[nvars] = SCIPfeasFloor(scip, SCIPvarGetUbLocal(var));
//       (bdchgdata->nvars)++;

//       assert(SCIPhashmapExists(bdchgdata->varhashmap, var) 
//          && (int)(size_t) SCIPhashmapGetImage(bdchgdata->varhashmap, var) == nvars);

//    }
   
//    /* get position of this variable */
//    pos = (int)(size_t) SCIPhashmapGetImage(bdchgdata->varhashmap, var);
   
//    /* update bounds if necessary */
//    if( boundtype == SCIP_BOUNDTYPE_LOWER )
//    {
//       if( bdchgdata->lbchgs[pos] < newbound )
//       {
//          bdchgdata->lbchgs[pos] = newbound;
//          (*nbdchgs)++;
//       }
      
//       if( newbound > bdchgdata->ubchgs[pos] )
//       {
//          *infeasible = TRUE;
//       }
      
//    } 
//    else 
//    {
//       if( newbound < bdchgdata->ubchgs[pos] )
//       {
//          bdchgdata->ubchgs[pos] = newbound;
//          (*nbdchgs)++;
//       }
//       if( newbound < bdchgdata->lbchgs[pos] )
//       {
//          *infeasible = TRUE;
//       }
//    }  

//    return SCIP_OKAY;
// }



// /* calculates variable bounds for an up-branch and a down-branch, supposig a LP or pseudo solution is given */
// static
// SCIP_RETCODE calculateBounds(
//    SCIP*                 scip,               /**< SCIP data structure */
//    SCIP_VAR*             branchvar,          /**< branching variable */
//    int*                  downlb,             /**< lower bound of variable in down branch */
//    int*                  downub,             /**< upper bound of variable in down branch */
//    int*                  uplb,               /**< lower bound of variable in up branch */
//    int*                  upub                /**< upper bound of variable in up branch */
//    )
// {
//    SCIP_Real varsol;
//    int lbdown;
//    int ubdown;
//    int lbup;
//    int ubup;

//    int lblocal;
//    int ublocal;

//    assert(scip != NULL);
//    assert(branchvar != NULL);

//    varsol = SCIPgetVarSol(scip, branchvar);

//    lblocal = SCIPfeasCeil(scip, SCIPvarGetLbLocal(branchvar));
//    ublocal = SCIPfeasFloor(scip, SCIPvarGetUbLocal(branchvar));

//    /* calculate bounds in down branch */
//    lbdown = lblocal;
   
//    /* in down branch: new upper bound is at most local upper bound - 1 */
//    ubdown = SCIPfeasFloor(scip, varsol) ;
//    if( ubdown == ublocal )
//       ubdown--;      

//    assert(lbdown <= ubdown);

//    /* calculate bounds in up branch */
//    ubup = ublocal;
      
//    /* in right branch: new lower bound is at least local lower bound + 1 */
//    lbup = SCIPfeasCeil(scip, varsol);
//    if( lbup == lblocal )
//       lbup++;

//    assert(lbup <= ubup);

//    /* ensure that both branches partition the domain */
//    if( lbup == ubdown )
//    {
//       int middle = (lblocal + ublocal) / 2; /* implicit rounding */
      
//       if( lbup <= middle )
//          ubdown--;
//       else 
//          lbup++;
//    }

//    /* ensure a real partition of the domain */
//    assert(ubdown < lbup);
//    assert(lbdown <= ubdown);
//    assert(lbup <= ubup);

//    /* set return values */
//    if( downlb != NULL )
//       *downlb = lbdown;
//    if( downub != NULL )
//       *downub = ubdown;
//    if( uplb != NULL )
//       *uplb = lbup;
//    if( upub != NULL )
//       *upub = ubup;
 
//    return SCIP_OKAY;
// }


// /** applies probing of a single variable in the given direction, and stores evaluation in given arrays */
// static
// SCIP_RETCODE applyProbing(
//    SCIP*                 scip,               /**< SCIP data structure */
//    SCIP_VAR**            vars,               /**< problem variables */
//    int                   nvars,              /**< number of problem variables */
//    SCIP_VAR*             probingvar,         /**< variable to perform probing on */
//    SCIP_Bool             probingdir,         /**< value to fix probing variable to */
//    SCIP_Bool             solvelp,            /**< value to decide whether pricing loop shall be performed */
//    SCIP_Real*            proplbs,            /**< array to store lower bounds after full propagation */
//    SCIP_Real*            propubs,            /**< array to store upper bounds after full propagation */
//    SCIP_Real*            lpobjvalue,         /**< pointer to store the lp obj value if lp was solved */
//    SCIP_Bool*            lpsolved,           /**< pointer to store whether the lp was solved */
//    SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occured or the
//                                              *   solving process should be stopped (e.g., due to a time limit) */
//    SCIP_Bool*            cutoff              /**< pointer to store whether the probing direction is infeasible */
//    )
// {


//    SCIP_LPSOLSTAT lpsolstat;
//    SCIP_Real varsol;
//    int leftlbprobing;
//    int leftubprobing;
//    int rightlbprobing;
//    int rightubprobing;

//    leftubprobing = -1;
//    leftlbprobing = -1;
//    rightlbprobing = -1;
//    rightubprobing = -1;

//    assert(proplbs != NULL);
//    assert(propubs != NULL);
//    assert(cutoff != NULL);
//    assert(SCIPvarGetLbLocal(probingvar) - 0.5 < SCIPvarGetUbLocal(probingvar));
//    assert(SCIPisFeasIntegral(scip, SCIPvarGetLbLocal(probingvar)));
//    assert(SCIPisFeasIntegral(scip, SCIPvarGetUbLocal(probingvar)));
         
//    assert(!solvelp || (lpsolved!=NULL && lpobjvalue!=NULL && lperror!=NULL));
   
//    varsol = SCIPgetVarSol(scip, probingvar);
   
//    if( probingdir == FALSE )
//    {
//       SCIP_CALL( calculateBounds(scip, probingvar, 
//             &leftlbprobing, &leftubprobing, NULL, NULL) );
//    }
//    else
//    {
//       SCIP_CALL( calculateBounds(scip, probingvar, 
//             NULL, NULL, &rightlbprobing, &rightubprobing) );
//    }

//    SCIPdebugMessage("applying probing on variable <%s> == %u [%d,%d] (nlocks=%d/%d, impls=%d/%d, clqs=%d/%d)\n",
//       SCIPvarGetName(probingvar), probingdir,
//       probingdir ? rightlbprobing : leftlbprobing, probingdir ? rightubprobing : leftubprobing,
//       SCIPvarGetNLocksDown(probingvar), SCIPvarGetNLocksUp(probingvar),
//       SCIPvarGetNImpls(probingvar, FALSE), SCIPvarGetNImpls(probingvar, TRUE),
//       SCIPvarGetNCliques(probingvar, FALSE), SCIPvarGetNCliques(probingvar, TRUE));

//    /* start probing mode */
//    SCIP_CALL( SCIPstartProbing(scip) );

//    *lpsolved = FALSE;
//    *lperror = FALSE;

//    /* change variable bounds for the probing directions*/
//    if( probingdir == FALSE )
//    {
//       SCIP_CALL( SCIPchgVarUbProbing(scip, probingvar, leftubprobing) );
//    }
//    else
//    {
//       SCIP_CALL( SCIPchgVarLbProbing(scip, probingvar, rightlbprobing) );
//    }

//    /* apply propagation of implication graph and clique table */
// //    SCIP_CALL( SCIPpropagateProbingImplications(scip, cutoff) );
// //    if( !(*cutoff) )
// //    {
// //       int i;

// //       for( i = 0; i < nvars; ++i )
// //       {
// //          impllbs[i] = SCIPvarGetLbLocal(vars[i]);
// //          implubs[i] = SCIPvarGetUbLocal(vars[i]);
// //       }
// //    }

//    /* apply propagation */
//    if( !(*cutoff) )
//    {
//       SCIP_CALL( SCIPpropagateProbing(scip, -1 /* maxproprounds */, cutoff, NULL) );
//    }

//    /* evaluate propagation */
//    if( !(*cutoff) )
//    {
//       int i;

//       for( i = 0; i < nvars; ++i )
//       {
//          proplbs[i] = SCIPvarGetLbLocal(vars[i]);
//          propubs[i] = SCIPvarGetUbLocal(vars[i]);
//       }
//    }


//    /* if parameter is set, we want to use the outcome of the LP relaxation */
//    if( !(*cutoff) && solvelp )
//    {
//       /* @todo: parameter with or without pricing */
//       //SCIP_CALL( SCIPsolveProbingLP( scip, -1 /*iterlimit*/, lperror ) );

//       SCIP_CALL( SCIPsolveProbingLPWithPricing( scip, FALSE/* pretendroot */, FALSE /*displayinfo*/,
//            -1 /*maxpricerounds*/, lperror ) );
//       lpsolstat = SCIPgetLPSolstat(scip);

//       if( !(*lperror) )
//       {
//          /* get LP solution status, objective value */
//          *cutoff = *cutoff || (lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT || lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE);
//          if( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
//          {
//             *lpobjvalue = SCIPgetLPObjval(scip);
//             *lpsolved = TRUE;
//          }
//       }
//       else 
//       {
//          SCIPinfoMessage(scip, NULL, "something went wrong, an lp error occured\n");         
//       }
//    }


//    /* exit probing mode */
//    SCIP_CALL( SCIPendProbing(scip) );

//    SCIPdebugMessage("probing results in cutoff/lpsolved/lpobj: %s / %s / %g\n", 
//       *cutoff?"cutoff":"no cutoff", *lpsolved?"lpsolved":"lp not solved", *lpobjvalue);

//    return SCIP_OKAY;
// }





// /** performs probing according to a fixed branching scheme */
// SCIP_RETCODE SCIPperformProbingLinking(
//    SCIP*                 scip,               /**< SCIP data structure */
//    SCIP_CONS*            cons,               /**< linking constraint to be checked */
//    BDCHGDATA*            bdchgdata,          /**< structure containing bound changes for almost all variables */
//    int                   itlim,              /**< iteration limit for strong branchings */
//    SCIP_Bool             solvelp,            /**< value to decide whether pricing loop shall be performed */
//    SCIP_Real*            down,               /**< stores dual bound after branching column down */
//    SCIP_Real*            up,                 /**< stores dual bound after branching column up */
//    SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound, or NULL;
//                                               *   otherwise, it can only be used as an estimate value */
//    SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound, or NULL;
//                                               *   otherwise, it can only be used as an estimate value */
//    SCIP_Bool*            downinf,            /**< pointer to store whether the downwards branch is infeasible, or NULL */
//    SCIP_Bool*            upinf,              /**< pointer to store whether the upwards branch is infeasible, or NULL */
//    SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occured or the
//                                               *   solving process should be stopped (e.g., due to a time limit) */
//    int*                  nbdchgs             /**< pointer to store number of total bound changes */
// )
// {

//    /* data for variables and bdchg arrays */
//    SCIP_VAR** probvars;
//    SCIP_VAR** vars;
//    SCIP_CONSDATA* consdata;
//    SCIP_VAR* probingvar;
//    int nvars;
//    int nintvars;
//    int nbinvars;

//    SCIP_Real* leftproplbs;
//    SCIP_Real* leftpropubs;
//    SCIP_Real* rightproplbs;
//    SCIP_Real* rightpropubs;
//    SCIP_VARTYPE vartype;

//    SCIP_Real leftlpbound;
//    SCIP_Real rightlpbound;
//    SCIP_Bool leftlpsolved;
//    SCIP_Bool rightlpsolved;
//    SCIP_Bool leftlperror;
//    SCIP_Bool rightlperror;
//    SCIP_Bool leftcutoff;
//    SCIP_Bool rightcutoff;
   
//    SCIP_Bool delay;
//    SCIP_Bool cutoff;

//    int i;

//    SCIP_Real varsol;
//    int j;

//    assert(lperror != NULL);

//    if( downvalid != NULL )
//       *downvalid = FALSE;
//    if( upvalid != NULL )
//       *upvalid = FALSE;
//    if( downinf != NULL )
//       *downinf = FALSE;
//    if( upinf != NULL )
//       *upinf = FALSE;
   
//    /* get constraint data and integer start time variable */
//    consdata = SCIPconsGetData(cons);
//    assert(consdata != NULL);
//    probingvar = consdata->intvar;
//    assert(consdata->intvar != NULL);


//    vartype = SCIPvarGetType(probingvar);

//     if( SCIPisStopped(scip) )
//    {
//       SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
//          "   (%.1fs) probing aborted: solving stopped\n", SCIPgetSolvingTime(scip));
//       return SCIP_OKAY;
//    }   
   
//    /* get lp solution value of last run */
//    varsol = SCIPgetVarSol(scip, probingvar);

//    /* get all variables to store branching deductions of variable bounds */
//    /* get all variables and store them in array 'vars' */
//    SCIP_CALL( SCIPgetVarsData(scip, &probvars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   
//    SCIP_CALL( SCIPduplicateMemoryArray(scip, &vars, probvars, nvars) );
   
//    /* capture variables to make sure, the variables are not deleted */
//    for( i = 0; i < nvars; ++i )
//    {
//       SCIP_CALL( SCIPcaptureVar(scip, vars[i]) );
//    }
   
//    /* get temporary memory for storing probing results */
//    SCIP_CALL( SCIPallocBufferArray(scip, &leftproplbs, nvars) );
//    SCIP_CALL( SCIPallocBufferArray(scip, &leftpropubs, nvars) );
//    SCIP_CALL( SCIPallocBufferArray(scip, &rightproplbs, nvars) );
//    SCIP_CALL( SCIPallocBufferArray(scip, &rightpropubs, nvars) );

//    /* for each binary variable, probe fixing the variable to left and right */
//    delay = FALSE;
//    cutoff = FALSE;
//    leftcutoff = FALSE;
//    rightcutoff = FALSE;
  
//    /* left branch: apply probing for setting ub to LP solution value  */
//    SCIP_CALL( applyProbing(scip, vars, nvars, probingvar, FALSE, solvelp, 
//          leftproplbs, leftpropubs,
//          &leftlpbound, &leftlpsolved, &leftlperror, &leftcutoff) );
   
//    if( leftcutoff )
//    {
//       int newbound;

//       SCIP_CALL( calculateBounds(scip, probingvar, 
//             NULL, NULL, &newbound, NULL) );

//       // newbound = SCIPfeasCeil(scip, varsol);
// //       if( SCIPisFeasEQ(scip, newbound, SCIPvarGetLbLocal(probingvar)) )
// //          newbound++;

//       /* lower bound can be updated */
//       SCIPdebugMessage("change lower bound of probing variable <%s> from %g to %d, nlocks=(%d/%d)\n", 
//          SCIPvarGetName(probingvar), SCIPvarGetLbLocal(probingvar), newbound,
//          SCIPvarGetNLocksDown(probingvar), SCIPvarGetNLocksUp(probingvar));

//       SCIP_CALL( addBdchg(scip, bdchgdata, probingvar, newbound, SCIP_BOUNDTYPE_LOWER, nbdchgs, &cutoff) );
//    }
      
//    if( !cutoff )
//    {
//       /* right branch: apply probing for setting lb to LP solution value  */
//       SCIP_CALL( applyProbing(scip, vars, nvars, probingvar, TRUE, solvelp, 
//             rightproplbs, rightpropubs,
//             &rightlpbound, &rightlpsolved, &rightlperror, &rightcutoff) );
      
//       if( rightcutoff )
//       {
//          int newbound;

//          SCIP_CALL( calculateBounds(scip, probingvar, 
//                NULL, &newbound, NULL, NULL) );

//          // newbound = SCIPfeasFloor(scip, varsol);
// //          if( SCIPisFeasEQ(scip, newbound, SCIPvarGetUbLocal(probingvar)) )
// //          newbound--;

//          /* upper bound can be updated */
//          SCIPdebugMessage("change probing variable <%s> upper bound from %g to %d, nlocks=(%d/%d)\n",
//             SCIPvarGetName(probingvar), SCIPvarGetUbLocal(probingvar), newbound,
//             SCIPvarGetNLocksDown(probingvar), SCIPvarGetNLocksUp(probingvar));
         
//          SCIP_CALL( addBdchg(scip, bdchgdata, probingvar, newbound, SCIP_BOUNDTYPE_UPPER, nbdchgs, &cutoff) );
         
//       }
//    }

//    /* set return value of lperror */
//    cutoff = cutoff || (leftcutoff && rightcutoff);
//    *lperror = leftlperror || rightlperror;
   

//    /* analyze probing deductions */
   
//    /* 1. dualbounds */
//    if( leftlpsolved  )
//       *down = leftlpbound;
//    if( rightlpsolved )
//       *up = rightlpbound;
   
//    /* 2. update bounds */
//    for( j = 0; j < nvars && !cutoff; ++j )
//    {
//       SCIP_Real newlb;
//       SCIP_Real newub;
      
//       if( vars[j] == probingvar )
//          continue;
      
//       /* new bounds of the variable is the union of the propagated bounds of the left and right case */
//       newlb = MIN(leftproplbs[j], rightproplbs[j]);
//       newub = MAX(leftpropubs[j], rightpropubs[j]);
      
//       /* check for fixed variables */
//       if( SCIPisFeasEQ(scip, newlb, newub) )
//       {
//          /* in both probings, variable j is deduced to a fixed value */
//          SCIP_CALL( addBdchg(scip, bdchgdata, vars[j], newlb, SCIP_BOUNDTYPE_LOWER, nbdchgs, &cutoff) );
//          SCIP_CALL( addBdchg(scip, bdchgdata, vars[j], newub, SCIP_BOUNDTYPE_UPPER, nbdchgs, &cutoff) );
//          continue;
//       } 
//       else 
//       {
      
//          SCIP_Real oldlb;
//          SCIP_Real oldub;

//          assert(SCIPvarGetType(vars[j]) == SCIP_VARTYPE_BINARY || SCIPvarGetType(vars[j]) == SCIP_VARTYPE_INTEGER);

//          /* check for bound tightenings */
//          oldlb = SCIPvarGetLbLocal(vars[j]);
//          oldub = SCIPvarGetUbLocal(vars[j]);
//          if( SCIPisLbBetter(scip, newlb, oldlb, oldub) )
//          {
//             /* in both probings, variable j is deduced to be at least newlb: tighten lower bound */
//             SCIP_CALL( addBdchg(scip, bdchgdata, vars[j], newlb, SCIP_BOUNDTYPE_LOWER, nbdchgs, &cutoff) );
//          }
//          if( SCIPisUbBetter(scip, newub, oldlb, oldub) && !cutoff )
//          {
//             /* in both probings, variable j is deduced to be at most newub: tighten upper bound */
//             SCIP_CALL( addBdchg(scip, bdchgdata, vars[j], newub, SCIP_BOUNDTYPE_UPPER, nbdchgs, &cutoff) );
//          }
         
//       } 
            
//    } /* end check for deductions */

//    /* set correct return values */
//    if( down != NULL && leftlpsolved )
//       *down = leftlpbound;
//    if( up != NULL && rightlpsolved )
//       *up = rightlpbound;
//    if( downvalid != NULL && leftlpsolved )
//       *downvalid = TRUE;
//    if( downvalid != NULL && !leftlpsolved )
//       *downvalid = FALSE;
//    if( upvalid != NULL && rightlpsolved )
//       *upvalid = TRUE;
//    if( upvalid != NULL && !rightlpsolved )
//       *upvalid = FALSE;
//    if( downinf != NULL )
//       *downinf = leftcutoff;
//    if( upinf != NULL )
//       *upinf = rightcutoff;

//    if( cutoff )
//    {
//       *downinf = cutoff;
//       *upinf = cutoff;
//    }

//    /* free temporary memory */
//    SCIPfreeBufferArray(scip, &rightpropubs);
//    SCIPfreeBufferArray(scip, &rightproplbs);
//    SCIPfreeBufferArray(scip, &leftpropubs);
//    SCIPfreeBufferArray(scip, &leftproplbs);

//    /* capture variables to make sure, the variables are not deleted */
//    for( i = 0; i < nvars; ++i )
//    {
//       SCIP_CALL( SCIPreleaseVar(scip, &vars[i]) );
//    }
   
//    SCIPfreeMemoryArray(scip, &vars);

//    return SCIP_OKAY;

// }



// /** performs branching on the given constraint */
// SCIP_RETCODE SCIPperformBranchingLinking(
//    SCIP*                 scip,               /**< SCIP data structure */
//    SCIP_CONS*            cons,               /**< linking constraint to be checked */
//    SCIP_Bool             bestisstrongbranch,
//    SCIP_Real             bestsbdown,
//    SCIP_Real             bestsbup,
//    SCIP_Bool             allcolsinlp, 
//    SCIP_Bool             exactsolve,   
//    SCIP_Real             cutoffbound,
//    SCIP_RESULT*          result              /**< pointer to store result */
//    )
// {
//    SCIP_CONSDATA* consdata;
//    SCIP_NODE* downchild;
//    SCIP_NODE* upchild;
//    SCIP_VAR* var;
//    SCIP_Real proveddown;
//    SCIP_Real provedup;
//    int downbound;
//    int upbound;
     
//    assert(scip != NULL);
//    assert(cons != NULL);
//    assert(*result == SCIP_DIDNOTRUN);
   

//    assert(SCIPisLT(scip, provedbound, cutoffbound));

//    /* get constraint data */
//    consdata = SCIPconsGetData(cons);
//    var = consdata->intvar;
   
   
//    /* recalculate bounds of branching variable in up and down branch */
//    SCIP_CALL( calculateBounds(scip, var, 
//          NULL, &downbound, &upbound, NULL) );
   

//    SCIPdebugMessage("perform branching on var <%s> with leftbranch %d and rightbranch %d\n", 
//       SCIPvarGetName(var), downbound, upbound);
      
//    /* create children and change bounds of branching variable */
//    SCIP_CALL( SCIPcreateChild(scip, &downchild,  SCIPcalcNodeselPriority(scip, var, downbound), 
//          SCIPcalcChildEstimate(scip, var, downbound ) ) );
//    SCIP_CALL( SCIPcreateChild(scip, &upchild, SCIPcalcNodeselPriority(scip, var, upbound),
//          SCIPcalcChildEstimate(scip, var, upbound ) ) );
   
//    SCIP_CALL( SCIPchgVarUbNode(scip, downchild, var, downbound) );
//    SCIP_CALL( SCIPchgVarLbNode(scip, upchild, var, upbound) );
   
//    SCIP_CALL( applyBdchgs(scip, bdchgdata, downchild) );
//    SCIP_CALL( applyBdchgs(scip, bdchgdata, upchild) );
   
//    assert(downchild != NULL);
//    assert(upchild != NULL);
   
//    /* calculate proved lower bounds for children */
//    proveddown = provedbound;
//    provedup = provedbound;
//    if( bestisstrongbranch )
//    {
//       if( bestsbdownvalid )
//          proveddown = MAX(provedbound, bestsbdown);
//       if( bestsbupvalid )
//          provedup = MAX(provedbound, bestsbup);
//    }
   
//    /* update the lower bounds in the children */
//    if( allcolsinlp && !exactsolve )
//    {
//       assert(SCIPisLT(scip, proveddown, cutoffbound));
//       assert(SCIPisLT(scip, provedup, cutoffbound));
//       SCIP_CALL( SCIPupdateNodeLowerbound(scip, downchild, proveddown) );
//       SCIP_CALL( SCIPupdateNodeLowerbound(scip, upchild, provedup) );
//    }
//    SCIPdebugMessage(" -> down child's lowerbound: %g\n", SCIPnodeGetLowerbound(downchild));
//    SCIPdebugMessage(" -> up child's lowerbound  : %g\n", SCIPnodeGetLowerbound(upchild));
   
//    *result = SCIP_BRANCHED;
   
//    return SCIP_OKAY;
// }

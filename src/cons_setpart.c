/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_setpart.c
 * @brief  constraint handler for set partitioning constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_setpart.h"
#include "cons_setcover.h"
#include "cons_linear.h"


#define CONSHDLR_NAME          "setpart"
#define CONSHDLR_DESC          "set partitioning constraint"
#define CONSHDLR_SEPAPRIORITY   +700200
#define CONSHDLR_ENFOPRIORITY   +700200
#define CONSHDLR_CHECKPRIORITY  -700200
#define CONSHDLR_SEPAFREQ             4
#define CONSHDLR_PROPFREQ            -1
#define CONSHDLR_NEEDSCONS         TRUE

#define EVENTHDLR_NAME         "setpart"
#define EVENTHDLR_DESC         "bound change event handler for set partitioning constraints"

#define LINCONSUPGD_PRIORITY    +700200

#define DEFAULT_NPSEUDOBRANCHES       2  /**< number of children created in pseudo branching */
#define MINBRANCHWEIGHT             0.3  /**< minimum weight of both sets in set partitioning branching */
#define MAXBRANCHWEIGHT             0.7  /**< maximum weight of both sets in set partitioning branching */


/** constraint handler data */
struct ConsHdlrData
{
   INTARRAY*        varuses;            /**< number of times a var is used in the active set partitioning constraints */
   int              npseudobranches;    /**< number of children created in pseudo branching */
};

/** set partitioning constraint data */
struct SetpartCons
{
   VAR**            vars;               /**< variables of the constraint */
   int              varssize;           /**< size of vars array */
   int              nvars;              /**< number of variables in the constraint */
   int              nfixedzeros;        /**< actual number of variables fixed to zero in the constraint */
   int              nfixedones;         /**< actual number of variables fixed to one in the constraint */
   unsigned int     local:1;            /**< is constraint only valid locally? */
   unsigned int     modifiable:1;       /**< is data modifiable during node processing (subject to column generation)? */
   unsigned int     removeable:1;       /**< should the row be removed from the LP due to aging or cleanup? */
   unsigned int     transformed:1;      /**< does the constraint data belongs to the transformed problem? */
   unsigned int     changed:1;          /**< was constraint changed since last preprocess/propagate call? */
};
typedef struct SetpartCons SETPARTCONS; /**< set partitioning constraint data */

/** constraint data for set partitioning constraints */
struct ConsData
{
   SETPARTCONS*     setpartcons;        /**< set partitioning constraint data */
   ROW*             row;                /**< LP row, if constraint is already stored in LP row format */
};




/*
 * Local methods
 */

/** creates constaint handler data for set partitioning constraint handler */
static
RETCODE conshdlrdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);

   CHECK_OKAY( SCIPallocMemory(scip, conshdlrdata) );
   CHECK_OKAY( SCIPcreateIntarray(scip, &(*conshdlrdata)->varuses) );
   (*conshdlrdata)->npseudobranches = DEFAULT_NPSEUDOBRANCHES;

   return SCIP_OKAY;
}

/** frees constraint handler data for set partitioning constraint handler */
static
RETCODE conshdlrdataFree(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   int i;

   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   CHECK_OKAY( SCIPfreeIntarray(scip, &(*conshdlrdata)->varuses) );
   SCIPfreeMemory(scip, conshdlrdata);

   return SCIP_OKAY;
}

/** increases the usage counter of the given variable */
static
RETCODE conshdlrdataIncVaruses(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   VAR*             var                 /**< variable to increase usage counter for */
   )
{
   INTARRAY* varuses;

   assert(conshdlrdata != NULL);
   assert(var != NULL);

   varuses = conshdlrdata->varuses;
   assert(varuses != NULL);

   CHECK_OKAY( SCIPincIntarrayVal(scip, varuses, SCIPvarGetIndex(var), +1) );

   /*debugMessage("varuses of <%s>: %d\n", SCIPvarGetName(var), SCIPgetIntarrayVal(scip, varuses, SCIPvarGetIndex(var)));*/

   return SCIP_OKAY;
}

/** decreases the usage counter of the given variable */
static
RETCODE conshdlrdataDecVaruses(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   VAR*             var                 /**< variable to increase usage counter for */
   )
{
   INTARRAY* varuses;

   assert(conshdlrdata != NULL);
   assert(var != NULL);

   varuses = conshdlrdata->varuses;
   assert(varuses != NULL);

   CHECK_OKAY( SCIPincIntarrayVal(scip, varuses, SCIPvarGetIndex(var), -1) );
   assert(SCIPgetIntarrayVal(scip, varuses, SCIPvarGetIndex(var)) >= 0);

   /*debugMessage("varuses of <%s>: %d\n", SCIPvarGetName(var), SCIPgetIntarrayVal(scip, varuses, SCIPvarGetIndex(var)));*/

   return SCIP_OKAY;
}

/** creates event data for variable at given position, and catches events */
static
RETCODE setpartconsCatchEvent(
   SCIP*            scip,               /**< SCIP data structure */
   SETPARTCONS*     setpartcons,        /**< set partitioning constraint object */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< array position of variable to catch bound change events for */
   )
{
   VAR* var;

   assert(setpartcons != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < setpartcons->nvars);
   assert(setpartcons->vars != NULL);

   var = setpartcons->vars[pos];
   assert(var != NULL);

   /* catch bound change events on variables */
   CHECK_OKAY( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, (EVENTDATA*)setpartcons) );
   
   /* update the fixed variables counters for this variable */
   if( SCIPisEQ(scip, SCIPvarGetUbLocal(var), 0.0) )
      setpartcons->nfixedzeros++;
   else if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), 1.0) )
      setpartcons->nfixedones++;

   return SCIP_OKAY;
}

/** deletes event data for variable at given position, and drops events */
static
RETCODE setpartconsDropEvent(
   SCIP*            scip,               /**< SCIP data structure */
   SETPARTCONS*     setpartcons,        /**< set partitioning constraint object */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< array position of variable to catch bound change events for */
   )
{
   VAR* var;

   assert(setpartcons != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < setpartcons->nvars);
   assert(setpartcons->vars != NULL);

   var = setpartcons->vars[pos];
   assert(var != NULL);
   
   CHECK_OKAY( SCIPdropVarEvent(scip, var, eventhdlr, (EVENTDATA*)setpartcons) );

   /* update the fixed variables counters for this variable */
   if( SCIPisEQ(scip, SCIPvarGetUbLocal(var), 0.0) )
      setpartcons->nfixedzeros--;
   else if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), 1.0) )
      setpartcons->nfixedones--;

   return SCIP_OKAY;
}

/** catches bound change events and locks rounding for variable at given position in transformed set partitioning constraint */
static
RETCODE setpartconsLockCoef(
   SCIP*            scip,               /**< SCIP data structure */
   SETPARTCONS*     setpartcons,        /**< set partitioning constraint object */
   EVENTHDLR*       eventhdlr,          /**< event handler for bound change events, or NULL */
   int              pos                 /**< position of variable in set partitioning constraint */
   )
{
   VAR* var;
      
   assert(scip != NULL);
   assert(setpartcons != NULL);
   assert(setpartcons->transformed);
   assert(0 <= pos && pos < setpartcons->nvars);

   var = setpartcons->vars[pos];
   
   /*debugMessage("locking coefficient <%s> in set partitioning constraint\n", SCIPvarGetName(var));*/

   if( eventhdlr == NULL )
   {
      /* get event handler for updating set partitioning constraint activity bounds */
      eventhdlr = SCIPfindEventHdlr(scip, EVENTHDLR_NAME);
      if( eventhdlr == NULL )
      {
         errorMessage("event handler for set partitioning constraints not found");
         return SCIP_PLUGINNOTFOUND;
      }
   }

   /* catch bound change events on variable */
   assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL);
   CHECK_OKAY( setpartconsCatchEvent(scip, setpartcons, eventhdlr, pos) );

   /* forbid rounding of variable */
   if( !setpartcons->local )
      SCIPvarForbidRound(var);

   return SCIP_OKAY;
}

/** drops bound change events and unlocks rounding for variable at given position in transformed set partitioning constraint */
static
RETCODE setpartconsUnlockCoef(
   SCIP*            scip,               /**< SCIP data structure */
   SETPARTCONS*     setpartcons,        /**< set partitioning constraint object */
   EVENTHDLR*       eventhdlr,          /**< event handler for bound change events, or NULL */
   int              pos                 /**< position of variable in set partitioning constraint */
   )
{
   VAR* var;

   assert(scip != NULL);
   assert(setpartcons != NULL);
   assert(setpartcons->transformed);
   assert(0 <= pos && pos < setpartcons->nvars);

   var = setpartcons->vars[pos];

   /*debugMessage("unlocking coefficient <%s> in set partitioning constraint\n", SCIPvarGetName(var));*/

   if( eventhdlr == NULL )
   {
      /* get event handler for updating set partitioning constraint activity bounds */
      eventhdlr = SCIPfindEventHdlr(scip, EVENTHDLR_NAME);
      if( eventhdlr == NULL )
      {
         errorMessage("event handler for set partitioning constraints not found");
         return SCIP_PLUGINNOTFOUND;
      }
   }
   
   /* drop bound change events on variable */
   assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL);
   CHECK_OKAY( setpartconsDropEvent(scip, setpartcons, eventhdlr, pos) );

   /* allow rounding of variable */
   if( !setpartcons->local )
      SCIPvarAllowRound(var);

   return SCIP_OKAY;
}

/** catches bound change events and locks rounding for all variables in transformed set partitioning constraint */
static
RETCODE setpartconsLockAllCoefs(
   SCIP*            scip,               /**< SCIP data structure */
   SETPARTCONS*     setpartcons        /**< set partitioning constraint object */
   )
{
   EVENTHDLR* eventhdlr;
   int i;

   assert(scip != NULL);
   assert(setpartcons != NULL);
   assert(setpartcons->transformed);

   /* get event handler for updating set partitioning constraint activity bounds */
   eventhdlr = SCIPfindEventHdlr(scip, EVENTHDLR_NAME);
   if( eventhdlr == NULL )
   {
      errorMessage("event handler for set partitioning constraints not found");
      return SCIP_PLUGINNOTFOUND;
   }

   /* lock every single coefficient */
   for( i = 0; i < setpartcons->nvars; ++i )
   {
      CHECK_OKAY( setpartconsLockCoef(scip, setpartcons, eventhdlr, i) );
   }

   return SCIP_OKAY;
}

/** drops bound change events and unlocks rounding for all variables in transformed set partitioning constraint */
static
RETCODE setpartconsUnlockAllCoefs(
   SCIP*            scip,               /**< SCIP data structure */
   SETPARTCONS*     setpartcons        /**< set partitioning constraint object */
   )
{
   EVENTHDLR* eventhdlr;
   int i;

   assert(scip != NULL);
   assert(setpartcons != NULL);
   assert(setpartcons->transformed);

   /* get event handler for updating set partitioning constraint activity bounds */
   eventhdlr = SCIPfindEventHdlr(scip, EVENTHDLR_NAME);
   if( eventhdlr == NULL )
   {
      errorMessage("event handler for set partitioning constraints not found");
      return SCIP_PLUGINNOTFOUND;
   }

   /* unlock every single coefficient */
   for( i = 0; i < setpartcons->nvars; ++i )
   {
      CHECK_OKAY( setpartconsUnlockCoef(scip, setpartcons, eventhdlr, i) );
   }

   return SCIP_OKAY;
}

/** deletes coefficient at given position from set partitioning constraint object */
static
RETCODE setpartconsDelCoefPos(
   SCIP*            scip,               /**< SCIP data structure */
   SETPARTCONS*     setpartcons,        /**< set partitioning constraint object */
   int              pos                 /**< position of coefficient to delete */
   )
{
   VAR* var;

   assert(setpartcons != NULL);
   assert(0 <= pos && pos < setpartcons->nvars);

   var = setpartcons->vars[pos];
   assert(var != NULL);
   assert(setpartcons->transformed ^ (SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL));

   if( setpartcons->transformed )
   {
      /* drop bound change events and unlock the rounding of variable */
      CHECK_OKAY( setpartconsUnlockCoef(scip, setpartcons, NULL, pos) );
   }

   /* move the last variable to the free slot */
   setpartcons->vars[pos] = setpartcons->vars[setpartcons->nvars-1];
   setpartcons->nvars--;

   setpartcons->changed = TRUE;

   return SCIP_OKAY;
}

/** creates a set partitioning constraint data object */
static
RETCODE setpartconsCreate(
   SCIP*            scip,               /**< SCIP data structure */
   SETPARTCONS**    setpartcons,        /**< pointer to store the set partitioning constraint */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< variables of the constraint */
   Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   assert(setpartcons != NULL);
   assert(nvars == 0 || vars != NULL);

   CHECK_OKAY( SCIPallocBlockMemory(scip, setpartcons) );
   if( nvars > 0 )
   {
      VAR* var;
      int v;

      CHECK_OKAY( SCIPduplicateBlockMemoryArray(scip, &(*setpartcons)->vars, vars, nvars) );
      (*setpartcons)->varssize = nvars;
      (*setpartcons)->nvars = nvars;
   }
   else
   {
      (*setpartcons)->vars = NULL;
      (*setpartcons)->varssize = 0;
      (*setpartcons)->nvars = 0;
   }
   (*setpartcons)->nfixedzeros = 0;
   (*setpartcons)->nfixedones = 0;
   (*setpartcons)->local = FALSE;
   (*setpartcons)->modifiable = modifiable;
   (*setpartcons)->removeable = removeable;
   (*setpartcons)->transformed = FALSE;
   (*setpartcons)->changed = TRUE;

   return SCIP_OKAY;
}   

/** creates a transformed set partitioning constraint data object */
static
RETCODE setpartconsCreateTransformed(
   SCIP*            scip,               /**< SCIP data structure */
   SETPARTCONS**    setpartcons,        /**< pointer to store the set partitioning constraint */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< variables of the constraint */
   Bool             local,              /**< is constraint only locally valid? */
   Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   EVENTHDLR* eventhdlr;
   VAR* var;
   int i;

   assert(setpartcons != NULL);
   assert(nvars == 0 || vars != NULL);

   CHECK_OKAY( setpartconsCreate(scip, setpartcons, nvars, vars, modifiable, removeable) );
   (*setpartcons)->local = local;
   (*setpartcons)->transformed = TRUE;

   /* transform the variables */
   for( i = 0; i < (*setpartcons)->nvars; ++i )
   {
      var = (*setpartcons)->vars[i];
      assert(var != NULL);
      assert(SCIPisLE(scip, 0.0, SCIPvarGetLbLocal(var)));
      assert(SCIPisLE(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));
      assert(SCIPisLE(scip, SCIPvarGetUbLocal(var), 1.0));
      assert(SCIPisIntegral(scip, SCIPvarGetLbLocal(var)));
      assert(SCIPisIntegral(scip, SCIPvarGetUbLocal(var)));

      /* use transformed variables in constraint instead original ones */
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
      {
	 var = SCIPvarGetTransformed(var);
         (*setpartcons)->vars[i] = var;
         assert(var != NULL);
      }
      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL);
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
   }

   /* catch bound change events and lock the rounding of variables */
   CHECK_OKAY( setpartconsLockAllCoefs(scip, *setpartcons) );

   return SCIP_OKAY;
}

/** frees a set partitioning constraint data */
static
RETCODE setpartconsFree(
   SCIP*            scip,               /**< SCIP data structure */
   SETPARTCONS**    setpartcons        /**< pointer to store the set partitioning constraint */
   )
{
   assert(setpartcons != NULL);
   assert(*setpartcons != NULL);

   if( (*setpartcons)->transformed )
   {
      /* drop bound change events and unlock the rounding of variables */
      CHECK_OKAY( setpartconsUnlockAllCoefs(scip, *setpartcons) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*setpartcons)->vars, (*setpartcons)->varssize);
   SCIPfreeBlockMemory(scip, setpartcons);

   return SCIP_OKAY;
}

/** creates an LP row from a set partitioning constraint data object */
static
RETCODE setpartconsToRow(
   SCIP*            scip,               /**< SCIP data structure */
   SETPARTCONS*     setpartcons,        /**< set partitioning constraint data */
   const char*      name,               /**< name of the constraint */
   ROW**            row                 /**< pointer to an LP row data object */
   )
{
   int v;

   assert(setpartcons != NULL);
   assert(row != NULL);

   CHECK_OKAY( SCIPcreateRow(scip, row, name, 0, NULL, NULL, 1.0, 1.0,
                  setpartcons->local, setpartcons->modifiable, setpartcons->removeable) );
   
   for( v = 0; v < setpartcons->nvars; ++v )
   {
      CHECK_OKAY( SCIPaddVarToRow(scip, *row, setpartcons->vars[v], 1.0) );
   }

   return SCIP_OKAY;
}

/** prints set partitioning constraint to file stream */
static
void setpartconsPrint(
   SCIP*            scip,               /**< SCIP data structure */
   SETPARTCONS*     setpartcons,        /**< set partitioning constraint object */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   int v;

   assert(setpartcons != NULL);

   if( file == NULL )
      file = stdout;

   /* print coefficients */
   if( setpartcons->nvars == 0 )
      fprintf(file, "0 ");
   for( v = 0; v < setpartcons->nvars; ++v )
   {
      assert(setpartcons->vars[v] != NULL);
      fprintf(file, "+%s ", SCIPvarGetName(setpartcons->vars[v]));
   }

   /* print right hand side */
   fprintf(file, "= 1\n");
}

/** checks constraint for violation only looking at the fixed variables, applies further fixings if possible */
static
RETCODE processFixings(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< set partitioning constraint to be separated */
   Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   Bool*            reduceddom,         /**< pointer to store TRUE, if a domain reduction was found */
   Bool*            addcut,             /**< pointer to store whether this constraint must be added as a cut */
   Bool*            mustcheck           /**< pointer to store whether this constraint must be checked for feasibility */
   )
{
   CONSDATA* consdata;
   SETPARTCONS* setpartcons;

   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(reduceddom != NULL);
   assert(addcut != NULL);
   assert(mustcheck != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   setpartcons = consdata->setpartcons;
   assert(setpartcons != NULL);
   assert(setpartcons->nvars == 0 || setpartcons->vars != NULL);
   assert(0 <= setpartcons->nfixedzeros && setpartcons->nfixedzeros <= setpartcons->nvars);
   assert(0 <= setpartcons->nfixedones && setpartcons->nfixedones <= setpartcons->nvars);

   *addcut = FALSE;
   *mustcheck = FALSE;

   if( setpartcons->nfixedones >= 2 )
   {
      /* the constraint cannot be feasible, because two variables are already fixed to one */
      CHECK_OKAY( SCIPresetConsAge(scip, cons) );
      *cutoff = TRUE;
   }
   else if( setpartcons->nfixedzeros == setpartcons->nvars )
   {
      if( setpartcons->modifiable )
      {
         /* the constraint cannot be feasible with the active variable set, but due to additional pricing,
          * it may be feasible after the next pricing loop -> just insert it as a cut
          */
         *addcut = TRUE;
      }
      else
      {
         /* the constraint cannot be feasible, because all variables are fixed to zero */
         CHECK_OKAY( SCIPresetConsAge(scip, cons) );
         *cutoff = TRUE;
      }
   }
   else if( setpartcons->nfixedones == 1 && setpartcons->nfixedzeros == setpartcons->nvars - 1 )
   {
      /* all variables are fixed, and the constraint is satisfied -> disable constraint if it is not modifiable */
      if( !setpartcons->modifiable )
      {
         CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
      }
   }
   else if( setpartcons->nfixedones == 1 )
   {
      VAR** vars;
      VAR* var;
      Bool fixedonefound;
      Bool fixed;
      int nvars;
      int v;

      assert(setpartcons->nfixedzeros < setpartcons->nvars - 1);

      /* all other variables except the fixed one must be zero, and if the constraint is not modifiable (i.e. if not
       * more variables may appear), it can be disabled after the other variables are fixed to zero
       */
      
      /* check all variables, if fixings to zero can be applied */
      vars = setpartcons->vars;
      nvars = setpartcons->nvars;
      fixedonefound = FALSE;
      fixed = FALSE;
      for( v = 0; v < nvars; ++v )
      {
         var = vars[v];
         assert(!fixedonefound || SCIPisZero(scip, SCIPvarGetLbLocal(var)));
         assert(SCIPisZero(scip, SCIPvarGetUbLocal(var)) || SCIPisEQ(scip, SCIPvarGetUbLocal(var), 1.0));
         if( SCIPvarGetLbLocal(var) < 0.5 )
         {
            if( SCIPvarGetUbLocal(var) > 0.5 )
            {
               CHECK_OKAY( SCIPchgVarUb(scip, var, 0.0) );
               fixed = TRUE;
            }
         }
         else
            fixedonefound = TRUE;
      }
      /* at least one variable must have been unfixed */
      assert(fixedonefound && fixed);

      if( setpartcons->modifiable )
      {
         CHECK_OKAY( SCIPresetConsAge(scip, cons) );
      }
      else
      {
         CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
      }
      *reduceddom = TRUE;
   }
   else if( setpartcons->nfixedzeros == setpartcons->nvars - 1 && !setpartcons->modifiable )
   {
      VAR** vars;
      VAR* var;
      Bool fixed;
      int nvars;
      int v;

      assert(setpartcons->nfixedones == 0);

      /* only one variable is not fixed to zero in an unmodifiable constraint
       * -> this variable can be fixed to one, and the constraint can be disabled
       */
      
      /* search the single variable that can be fixed */
      vars = setpartcons->vars;
      nvars = setpartcons->nvars;
      fixed = FALSE;
      for( v = 0; v < nvars && !fixed; ++v )
      {
         var = vars[v];
         assert(SCIPisZero(scip, SCIPvarGetLbLocal(var)));
         assert(SCIPisZero(scip, SCIPvarGetUbLocal(var)) || SCIPisEQ(scip, SCIPvarGetUbLocal(var), 1.0));
         if( SCIPvarGetUbLocal(var) > 0.5 )
         {
            CHECK_OKAY( SCIPchgVarLb(scip, var, 1.0) );
            fixed = TRUE;
         }
      }
      assert(fixed);

      CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
      *reduceddom = TRUE;
   }
   else
   {
      /* no information can be obtained from the variable's fixings -> we have to check the constraint manually */
      *mustcheck = TRUE;
   }

   return SCIP_OKAY;
}

/** checks constraint for violation, returns TRUE iff constraint is feasible */
static
Bool check(
   SCIP*            scip,               /**< SCIP data structure */
   SETPARTCONS*     setpartcons,        /**< set partitioning constraint to be checked */
   SOL*             sol                 /**< primal CIP solution */
   )
{
   VAR** vars;
   Real solval;
   Real sum;
   int nvars;
   int v;
   
   /* calculate the constraint's activity */
   vars = setpartcons->vars;
   nvars = setpartcons->nvars;
   sum = 0.0;
   assert(SCIPfeastol(scip) < 0.1); /* to make the comparison against 1.1 working */
   for( v = 0; v < nvars && sum < 1.1; ++v )  /* if sum >= 1.1, it is clearly too large */
   {
      assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY);
      solval = SCIPgetSolVal(scip, sol, vars[v]);
      assert(SCIPisFeasGE(scip, solval, 0.0) && SCIPisFeasLE(scip, solval, 1.0));
      sum += solval;
   }
   return SCIPisFeasEQ(scip, sum, 1.0);
}

/** checks constraint for violation, and adds it as a cut if possible */
static
RETCODE separate(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< set partitioning constraint to be separated */
   Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   Bool*            separated,          /**< pointer to store TRUE, if a cut was found */
   Bool*            reduceddom          /**< pointer to store TRUE, if a domain reduction was found */
   )
{
   CONSDATA* consdata;
   SETPARTCONS* setpartcons;
   Bool addcut;
   Bool mustcheck;

   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(separated != NULL);
   assert(reduceddom != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   setpartcons = consdata->setpartcons;
   assert(setpartcons != NULL);
   assert(setpartcons->nvars == 0 || setpartcons->vars != NULL);
   assert(0 <= setpartcons->nfixedzeros && setpartcons->nfixedzeros <= setpartcons->nvars);
   assert(0 <= setpartcons->nfixedones && setpartcons->nfixedones <= setpartcons->nvars);

   /* check constraint for violation only looking at the fixed variables, apply further fixings if possible */
   CHECK_OKAY( processFixings(scip, cons, cutoff, reduceddom, &addcut, &mustcheck) );

   if( mustcheck )
   {
      ROW* row;

      assert(!addcut);

      /* variable's fixings didn't give us any information -> we have to check the constraint */
      row = consdata->row;
      if( row != NULL )
      {
         /* ignore constraints that are in the LP */
         if( !SCIProwIsInLP(row) )
         {         
            Real feasibility;
            
            feasibility = SCIPgetRowLPFeasibility(scip, row);
            addcut = !SCIPisFeasible(scip, feasibility);
         }
      }
      else
         addcut = !check(scip, setpartcons, NULL);

      if( !addcut )
      {
         /* constraint was feasible -> increase age */
         CHECK_OKAY( SCIPincConsAge(scip, cons) );
      }
   }

   if( addcut )
   {
      if( consdata->row == NULL )
      {
         /* convert set partitioning constraint data into LP row */
         CHECK_OKAY( setpartconsToRow(scip, setpartcons, SCIPconsGetName(cons), &consdata->row) );
      }
      assert(consdata->row != NULL);
      assert(!SCIProwIsInLP(consdata->row));
            
      /* insert LP row as cut */
      CHECK_OKAY( SCIPaddCut(scip, consdata->row, 1.0/(setpartcons->nvars+1)) );
      CHECK_OKAY( SCIPresetConsAge(scip, cons) );
      *separated = TRUE;
   }

   return SCIP_OKAY;
}

/** enforces the pseudo solution on the given constraint */
static
RETCODE enforcePseudo(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< set partitioning constraint to be separated */
   Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   Bool*            infeasible,         /**< pointer to store TRUE, if the constraint was infeasible */
   Bool*            reduceddom,         /**< pointer to store TRUE, if a domain reduction was found */
   Bool*            solvelp             /**< pointer to store TRUE, if the LP has to be solved */
   )
{
   Bool addcut;
   Bool mustcheck;

   assert(!SCIPhasActnodeLP(scip));
   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(infeasible != NULL);
   assert(reduceddom != NULL);
   assert(solvelp != NULL);

   /* check constraint for violation only looking at the fixed variables, apply further fixings if possible */
   CHECK_OKAY( processFixings(scip, cons, cutoff, reduceddom, &addcut, &mustcheck) );

   if( mustcheck )
   {
      CONSDATA* consdata;
      SETPARTCONS* setpartcons;

      assert(!addcut);

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);
      setpartcons = consdata->setpartcons;
      assert(setpartcons != NULL);

      if( check(scip, setpartcons, NULL) )
      {
         /* constraint was feasible -> increase age */
         CHECK_OKAY( SCIPincConsAge(scip, cons) );
      }
      else
      {
         /* constraint was infeasible -> reset age */
         CHECK_OKAY( SCIPresetConsAge(scip, cons) );
         *infeasible = TRUE;
      }
   }

   if( addcut )
   {
      /* a cut must be added to the LP -> we have to solve the LP immediately */
      CHECK_OKAY( SCIPresetConsAge(scip, cons) );
      *solvelp = TRUE;
   }

   return SCIP_OKAY;
}




/*
 * Callback methods of constraint handler
 */

static
DECL_CONSFREE(consFreeSetpart)
{
   CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   CHECK_OKAY( conshdlrdataFree(scip, &conshdlrdata) );

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}

static
DECL_CONSDELETE(consDeleteSetpart)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* free LP row and setpart constraint */
   if( (*consdata)->row != NULL )
   {
      CHECK_OKAY( SCIPreleaseRow(scip, &(*consdata)->row) );
   }
   CHECK_OKAY( setpartconsFree(scip, &(*consdata)->setpartcons) );

   /* free constraint data object */
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

static
DECL_CONSTRANS(consTransSetpart)
{
   CONSDATA* sourcedata;
   CONSDATA* targetdata;
   SETPARTCONS* setpartcons;

   /*debugMessage("Trans method of setpart constraints\n");*/

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPstage(scip) == SCIP_STAGE_INITSOLVE);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->row == NULL);  /* in original problem, there cannot be LP rows */
   assert(sourcedata->setpartcons != NULL);

   /* create constraint data for target constraint */
   CHECK_OKAY( SCIPallocBlockMemory(scip, &targetdata) );

   setpartcons = sourcedata->setpartcons;

   CHECK_OKAY( setpartconsCreateTransformed(scip, &targetdata->setpartcons, setpartcons->nvars, setpartcons->vars,
                  setpartcons->local, setpartcons->modifiable, setpartcons->removeable) );
   targetdata->row = NULL;

   /* create target constraint */
   assert(!SCIPconsIsPropagated(sourcecons));
   CHECK_OKAY( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
                  SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
                  SCIPconsIsPropagated(sourcecons)) );

   return SCIP_OKAY;
}

static
DECL_CONSSEPA(consSepaSetpart)
{
   CONSDATA* consdata;
   Bool cutoff;
   Bool separated;
   Bool reduceddom;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   debugMessage("separating %d/%d set partitioning constraints\n", nusefulconss, nconss);

   *result = SCIP_DIDNOTFIND;

   cutoff = FALSE;
   separated = FALSE;
   reduceddom = FALSE;

   /* step 1: check all useful set partitioning constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff && !reduceddom; ++c )
   {
      CHECK_OKAY( separate(scip, conss[c], &cutoff, &separated, &reduceddom) );
   }

   /* step 2: combine set partitioning constraints to get more cuts */
   todoMessage("further cuts of set partitioning constraints");

   /* step 3: if no cuts were found and we are in the root node, check remaining constraints for feasibility */
   if( SCIPgetActDepth(scip) == 0 )
   {
      for( c = nusefulconss; c < nconss && !cutoff && !separated && !reduceddom; ++c )
      {
         CHECK_OKAY( separate(scip, conss[c], &cutoff, &separated, &reduceddom) );
      }
   }

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( separated )
      *result = SCIP_SEPARATED;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


#if 0 /* ZU LANGSAM UND ZU SCHLECHT */
/** if fractional variables exist, chooses a set partitioning constraint C and a subset S of variables in C,
 *  and branches on (i) x(S) == 1, and (ii) x(S) == 0
 */
static
RETCODE branchLP(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLR*        conshdlr,           /**< set partitioning constraint handler */
   CONS**           conss,              /**< active set partitioning constraints */
   int              nconss,             /**< number of active set partitioning constraints */
   RESULT*          result              /**< pointer to store the result SCIP_BRANCHED, if branching was applied */
   )
{
   CONSHDLRDATA* conshdlrdata;
   INTARRAY* varuses;
   CONS* cons;
   CONSDATA* consdata;
   SETPARTCONS* setpartcons;
   VAR** sortcands;
   VAR* var;
   NODE* node;
   Real solval;
   Real sum;
   int nlpcands;
   int bestcand;
   int bestoverlap;
   int thisoverlap;
   int uses;
   int nselcands;
   int c;
   int i;
   int j;

   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(nconss >= 1);
   assert(result != NULL);

   /* get fractional variables */
   CHECK_OKAY( SCIPgetLPBranchCands(scip, NULL, NULL, NULL, &nlpcands) );
   if( nlpcands == 0 )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   varuses = conshdlrdata->varuses;
   assert(varuses != NULL);

   /* find the set partitioning constraint, where the fractional variables overlap the most with other
    * set partitioning constraints
    */
   bestcand = -1;
   bestoverlap = -1;
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      setpartcons = consdata->setpartcons;
      assert(setpartcons != NULL);

      thisoverlap = 0;
      for( i = 0; i < setpartcons->nvars; ++i )
      {
         solval = SCIPgetVarSol(scip, setpartcons->vars[i]);
         assert(SCIPisGE(scip, solval, 0.0) && SCIPisLE(scip, solval, 1.0));
         if( !SCIPisIntegral(scip, solval) )
            thisoverlap += SCIPgetIntarrayVal(scip, varuses, SCIPvarGetIndex(setpartcons->vars[i]));
      }
      if( thisoverlap > bestoverlap )
      {
         bestcand = c;
         bestoverlap = thisoverlap;
      }
   }
   assert(0 <= bestcand && bestcand < nconss);
   cons = conss[bestcand];
   consdata = SCIPconsGetData(cons);
   setpartcons = consdata->setpartcons;

   /* get temporary memory */
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &sortcands, setpartcons->nvars) );
   
   /* order the constraint's variables by their overlap */
   for( i = 0; i < setpartcons->nvars; ++i )
   {
      var = setpartcons->vars[i];
      uses = SCIPgetIntarrayVal(scip, varuses, SCIPvarGetIndex(var));

      for( j = i; j > 0 && uses > SCIPgetIntarrayVal(scip, varuses, SCIPvarGetIndex(sortcands[j-1])); --j )
      {
         sortcands[j] = sortcands[j-1];
      }
      assert(0 <= j && j <= i);
      sortcands[j] = var;
   }

   /* find the candidate set by adding constraint's variables to the set ordered by their overlap until 0.5 is reached;
    * afterwards, check if the last added variable should be member of the set or not
    */
   sum = 0.0;
   solval = 0.0;
   for( nselcands = 0; nselcands < setpartcons->nvars && sum < 0.5; ++nselcands )
   {
      solval = SCIPgetVarSol(scip, sortcands[nselcands]);
      sum += solval;
   }
   assert(sum >= 0.5);
   if( sum - 0.5 > 0.5 - (sum-solval) )
   {
      sum -= solval;
      nselcands--;
   }
   assert(!SCIPisIntegral(scip, sum));

   /* perform the set partitioning branching on the selected variables */
   assert(0 < nselcands && nselcands < setpartcons->nvars);
   
   /* create left child, fix x_i = 0 for all i \in S */
   CHECK_OKAY( SCIPcreateChild(scip, &node) );
   for( i = 0; i < nselcands; ++i )
   {
      if( !SCIPisZero(scip, SCIPvarGetUbLocal(sortcands[i])) )
      {
         CHECK_OKAY( SCIPchgVarUbNode(scip, node, sortcands[i], 0.0) );
      }
   }
   
   /* create right child: add constraint x(S) == 1 */
   CHECK_OKAY( SCIPcreateChild(scip, &node) );

   /* all other variables of the constraint can be set to 0.0 */
   for( i = nselcands; i < setpartcons->nvars; ++i )
   {
      if( !SCIPisZero(scip, SCIPvarGetUbLocal(sortcands[i])) )
      {
         CHECK_OKAY( SCIPchgVarUbNode(scip, node, sortcands[i], 0.0) );
      }
   }

   if( nselcands == 1 )
   {
      /* only one candidate selected: fix it to 1.0 */
      debugMessage("fixing variable <%s> to 1.0 in right child node\n", SCIPvarGetName(sortcands[0]));
      assert(SCIPisZero(scip, SCIPvarGetLbLocal(sortcands[0])));
      CHECK_OKAY( SCIPchgVarLbNode(scip, node, sortcands[0], 1.0) );

      /* if the constraint is unmodifiable, it can be disabled, because all other variables are already set to 0.0 */
      if( !setpartcons->modifiable )
      {
         CHECK_OKAY( SCIPdisableConsNode(scip, node, cons) );
      }
   }
   else if( setpartcons->modifiable )
   {
      CONS* newcons;
      char name[MAXSTRLEN];
         
      /* for a modifiable constraint, add set partitioning constraint x(S) == 1, because this is stronger than
       * the original constraint with adding x(C\S) == 0 */
      sprintf(name, "SpB%lld", SCIPgetNodenum(scip));
      
      CHECK_OKAY( SCIPcreateConsSetpart(scip, &newcons, name, nselcands, sortcands,
                     TRUE, TRUE, FALSE, TRUE, FALSE, TRUE) );
      CHECK_OKAY( SCIPaddConsNode(scip, node, newcons) );
      CHECK_OKAY( SCIPreleaseCons(scip, &newcons) );
   }

   *result = SCIP_BRANCHED;
   
#ifdef DEBUG
   debugMessage("set partitioning branching: nselcands=%d/%d, weight(S)=%g, A={", 
      nselcands, setpartcons->nvars, sum);
   for( i = 0; i < nselcands; ++i )
      printf(" %s[%g]", SCIPvarGetName(sortcands[i]), SCIPgetSolVal(scip, NULL, sortcands[i]));
   printf(" }\n");
#endif
   
   /* free temporary memory */
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &sortcands) );

   return SCIP_OKAY;
}

#else

/** if fractional variables exist, chooses a set S of them and branches on (i) x(S) == 0, and (ii) x(S) >= 1 */
static
RETCODE branchLP(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLR*        conshdlr,           /**< set partitioning constraint handler */
   RESULT*          result              /**< pointer to store the result SCIP_BRANCHED, if branching was applied */
   )
{
   CONSHDLRDATA* conshdlrdata;
   INTARRAY* varuses;
   VAR** lpcands;
   VAR** sortcands;
   VAR* var;
   Real branchweight;
   Real solval;
   int nlpcands;
   int nsortcands;
   int nselcands;
   int uses;
   int i;
   int j;

   todoMessage("use a better set partitioning branching on LP solution (use SOS branching)");

   assert(conshdlr != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   varuses = conshdlrdata->varuses;
   assert(varuses != NULL);

   /* get fractional variables */
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, NULL, NULL, &nlpcands) );
   if( nlpcands == 0 )
      return SCIP_OKAY;

   /* get temporary memory */
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &sortcands, nlpcands) );
   
   /* sort fractional variables by number of uses in enabled set partitioning constraints */
   nsortcands = 0;
   for( i = 0; i < nlpcands; ++i )
   {
      var = lpcands[i];
      uses = SCIPgetIntarrayVal(scip, varuses, SCIPvarGetIndex(var));
      if( uses > 0 )
      {
         for( j = nsortcands; j > 0 && uses > SCIPgetIntarrayVal(scip, varuses, SCIPvarGetIndex(sortcands[j-1])); --j )
         {
            sortcands[j] = sortcands[j-1];
         }
         assert(0 <= j && j <= nsortcands);
         sortcands[j] = var;
         nsortcands++;
      }
   }
   assert(nsortcands <= nlpcands);

   /* if none of the fractional variables is member of a set partitioning constraint,
    * we are not responsible for doing the branching
    */
   if( nsortcands > 0 )
   {
      /* select the first variables from the sorted candidate list, until MAXBRANCHWEIGHT is reached;
       * then choose one less
       */
      branchweight = 0.0;
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
         NODE* node;

         /* perform the set partitioning branching on the selected variables */
         assert(nselcands <= nlpcands);
         
         /* create left child, fix x_i = 0 for all i \in S */
         CHECK_OKAY( SCIPcreateChild(scip, &node) );
         for( i = 0; i < nselcands; ++i )
         {
            CHECK_OKAY( SCIPchgVarUbNode(scip, node, sortcands[i], 0.0) );
         }

         /* create right child: add constraint x(S) >= 1 */
         CHECK_OKAY( SCIPcreateChild(scip, &node) );
         if( nselcands == 1 )
         {
            /* only one candidate selected: fix it to 1.0 */
            debugMessage("fixing variable <%s> to 1.0 in right child node\n", SCIPvarGetName(sortcands[0]));
            CHECK_OKAY( SCIPchgVarLbNode(scip, node, sortcands[0], 1.0) );
         }
         else
         {
            CONS* newcons;
            char name[MAXSTRLEN];
         
            /* add set covering constraint x(S) >= 1 */
            sprintf(name, "SpB%lld", SCIPgetNodenum(scip));

            CHECK_OKAY( SCIPcreateConsSetcover(scip, &newcons, name, nselcands, sortcands,
                           TRUE, TRUE, FALSE, TRUE, FALSE, TRUE) );
            CHECK_OKAY( SCIPaddConsNode(scip, node, newcons) );
            CHECK_OKAY( SCIPreleaseCons(scip, &newcons) );
         }
      
         *result = SCIP_BRANCHED;
         
#ifdef DEBUG
         debugMessage("set partitioning branching: nselcands=%d/%d, weight(S)=%g, A={", nselcands, nlpcands, branchweight);
         for( i = 0; i < nselcands; ++i )
            printf(" %s[%g]", SCIPvarGetName(sortcands[i]), SCIPgetSolVal(scip, NULL, sortcands[i]));
         printf(" }\n");
#endif
      }
   }

   /* free temporary memory */
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &sortcands) );

   return SCIP_OKAY;
}
#endif

/** if unfixed variables exist, chooses a set S of them and creates |S|+1 child nodes:
 *   - for each variable i from S, create child node with x_0 = ... = x_i-1 = 0, x_i = 1
 *   - create an additional child node x_0 = ... = x_n-1 = 0
 */
static
RETCODE branchPseudo(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLR*        conshdlr,           /**< set partitioning constraint handler */
   RESULT*          result              /**< pointer to store the result SCIP_BRANCHED, if branching was applied */
   )
{
   CONSHDLRDATA* conshdlrdata;
   INTARRAY* varuses;
   VAR** pseudocands;
   VAR** branchcands;
   VAR* var;
   NODE* node;
   int* canduses;
   int npseudocands;
   int maxnbranchcands;
   int nbranchcands;
   int uses;
   int i;
   int j;

   todoMessage("use a better set partitioning branching on pseudo solution (use SOS branching)");

   assert(conshdlr != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   varuses = conshdlrdata->varuses;
   assert(varuses != NULL);

   /* get fractional variables */
   CHECK_OKAY( SCIPgetPseudoBranchCands(scip, &pseudocands, &npseudocands) );
   if( npseudocands == 0 )
      return SCIP_OKAY;

   /* choose the maximal number of branching variables */
   maxnbranchcands = conshdlrdata->npseudobranches-1;
   assert(maxnbranchcands >= 1);

   /* get temporary memory */
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &branchcands, maxnbranchcands) );
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &canduses, maxnbranchcands) );
   
   /* sort fractional variables by number of uses in enabled set partitioning constraints */
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

   if( nbranchcands == 0 )
   {
      /* none of the unfixed variables is member of a set partitioning constraint
       * -> we are not responsible for doing the branching
       */
      return SCIP_OKAY;
   }
   
   /* branch on the first part of the sorted candidates:
    * - for each of these variables i, create a child node x_0 = ... = x_i-1 = 0, x_i = 1
    * - create an additional child node x_0 = ... = x_n-1 = 0
    */
   for( i = 0; i < nbranchcands; ++i )
   {            
      /* create child with x_0 = ... = x_i-1 = 0, x_i = 1 */
      CHECK_OKAY( SCIPcreateChild(scip, &node) );
      for( j = 0; j < i; ++j )
      {
         CHECK_OKAY( SCIPchgVarUbNode(scip, node, branchcands[j], 0.0) );
      }
      CHECK_OKAY( SCIPchgVarLbNode(scip, node, branchcands[i], 1.0) );
   }
   /* create child with x_0 = ... = x_n = 0 */
   CHECK_OKAY( SCIPcreateChild(scip, &node) );
   for( i = 0; i < nbranchcands; ++i )
   {
      CHECK_OKAY( SCIPchgVarUbNode(scip, node, branchcands[i], 0.0) );
   }

   *result = SCIP_BRANCHED;

#ifdef DEBUG
   {
      int nchildren;
      CHECK_OKAY( SCIPgetChildren(scip, NULL, &nchildren) );
      debugMessage("branched on pseudo solution: %d children\n", nchildren);
   }
#endif

   /* free temporary memory */
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &canduses) );
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &branchcands) );

   return SCIP_OKAY;
}



static
DECL_CONSENFOLP(consEnfolpSetpart)
{
   Bool cutoff;
   Bool separated;
   Bool reduceddom;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   debugMessage("LP enforcing %d set partitioning constraints\n", nconss);

   *result = SCIP_FEASIBLE;

   cutoff = FALSE;
   separated = FALSE;
   reduceddom = FALSE;

   /* step 1: check all useful set partitioning constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff && !reduceddom; ++c )
   {
      CHECK_OKAY( separate(scip, conss[c], &cutoff, &separated, &reduceddom) );
   }

   if( !cutoff && !separated && !reduceddom )
   {
      /* step 2: if solution is not integral, choose a variable set to branch on */
      CHECK_OKAY( branchLP(scip, conshdlr, result) );
      if( *result != SCIP_FEASIBLE )
         return SCIP_OKAY;
      
      /* step 3: check all obsolete set partitioning constraints for feasibility */
      for( c = nusefulconss; c < nconss && !cutoff && !separated && !reduceddom; ++c )
      {
         CHECK_OKAY( separate(scip, conss[c], &cutoff, &separated, &reduceddom) );
      }
   }

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( separated )
      *result = SCIP_SEPARATED;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}

static
DECL_CONSENFOPS(consEnfopsSetpart)
{
   Bool cutoff;
   Bool infeasible;
   Bool reduceddom;
   Bool solvelp;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   debugMessage("pseudo enforcing %d set partitioning constraints\n", nconss);

   *result = SCIP_FEASIBLE;

   cutoff = FALSE;
   infeasible = FALSE;
   reduceddom = FALSE;
   solvelp = FALSE;

   /* check all set partitioning constraints for feasibility */
   for( c = 0; c < nconss && !cutoff && !reduceddom && !solvelp; ++c )
   {
      CHECK_OKAY( enforcePseudo(scip, conss[c], &cutoff, &infeasible, &reduceddom, &solvelp) );
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
      
      /* at least one constraint is violated by pseudo solution and we didn't find a better way to resolve this:
       * -> branch on pseudo solution
       */
      CHECK_OKAY( branchPseudo(scip, conshdlr, result) );
   }
   
   return SCIP_OKAY;
}

static
DECL_CONSCHECK(consCheckSetpart)
{
   CONSDATA* consdata;
   CONS* cons;
   SETPARTCONS* setpartcons;
   VAR** vars;
   Real solval;
   int nvars;
   int c;
   int v;
   Bool found;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   /* check all set partitioning constraints for feasibility */
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);
      setpartcons = consdata->setpartcons;
      assert(setpartcons != NULL);
      if( checklprows || consdata->row == NULL || !SCIProwIsInLP(consdata->row) )
      {
         if( !check(scip, setpartcons, sol) )
         {
            /* constraint is violated */
            CHECK_OKAY( SCIPresetConsAge(scip, cons) );
            *result = SCIP_INFEASIBLE;
            return SCIP_OKAY;
         }
         else
         {
            CHECK_OKAY( SCIPincConsAge(scip, cons) );
         }
      }
   }
   
   return SCIP_OKAY;
}




/*
 * Presolving
 */

/** deletes all zero-fixed variables */
static
RETCODE setpartconsApplyFixings(
   SCIP*            scip,               /**< SCIP data structure */
   SETPARTCONS*     setpartcons         /**< set partitioning constraint object */
   )
{
   VAR* var;
   int v;

   assert(setpartcons != NULL);

   if( setpartcons->nfixedzeros >= 1 )
   {
      assert(setpartcons->vars != NULL);

      v = 0;
      while( v < setpartcons->nvars )
      {
         var = setpartcons->vars[v];
         if( SCIPisZero(scip, SCIPvarGetUbGlobal(var)) )
         {
            CHECK_OKAY( setpartconsDelCoefPos(scip, setpartcons, v) );
         }
         else
            ++v;
      }
   }

   return SCIP_OKAY;
}

static
DECL_CONSPRESOL(consPresolSetpart)
{
   CONS* cons;
   CONSDATA* consdata;
   SETPARTCONS* setpartcons;
   Bool infeasible;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* process constraints */
   for( c = 0; c < nconss && *result != SCIP_CUTOFF; ++c )
   {
      cons = conss[c];
      assert(cons != NULL);
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);
      setpartcons = consdata->setpartcons;
      assert(setpartcons != NULL);

      if( !setpartcons->changed )
         continue;

      debugMessage("presolving set partitioning constraint <%s>\n", SCIPconsGetName(cons));

      /* remove all variables that are fixed to zero */
      CHECK_OKAY( setpartconsApplyFixings(scip, setpartcons) );

      /* check for infeasibility due to more than one variable fixed to one */
      if( setpartcons->nfixedones >= 2 )
      {
         debugMessage("set partitioning constraint <%s> is infeasible\n", SCIPconsGetName(cons));
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      /* check, if we can fix all other variables to zero, because a variable fixed to one exists
       *  -> if the constraint is not modifiable, we can delete it (which is done later)
       */
      if( setpartcons->nfixedones == 1 )
      {
         VAR* var;
         int v;

         debugMessage("set partitioning constraint <%s> has a variable fixed to 1.0\n", SCIPconsGetName(cons));
         for( v = 0; v < setpartcons->nvars; ++v )
         {
            var = setpartcons->vars[v];
            if( SCIPisZero(scip, SCIPvarGetLbGlobal(var)) && !SCIPisZero(scip, SCIPvarGetUbGlobal(var)) )
            {
               CHECK_OKAY( SCIPfixVar(scip, var, 0.0, &infeasible) );
               assert(!infeasible);
               (*nfixedvars)++;
               *result = SCIP_SUCCESS;
            }
         }
      }

      if( !setpartcons->modifiable )
      {
         /* constraint is redundant, if one variable is fixed to one; the others are already fixed to zero */
         if( setpartcons->nfixedones == 1 )
         {
            assert(setpartcons->nfixedzeros == setpartcons->nvars-1);
            debugMessage("set partitioning constraint <%s> is redundant\n", SCIPconsGetName(cons));
            CHECK_OKAY( SCIPdelCons(scip, cons) );
            (*ndelconss)++;
            *result = SCIP_SUCCESS;
            continue;
         }

         /* check for infeasibility: are all variables fixed to zero? */
         if( setpartcons->nfixedzeros == setpartcons->nvars )
         {
            debugMessage("set partitioning constraint <%s> is infeasible\n", SCIPconsGetName(cons));
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         /* fix variable and delete constraint, if constraint consists only of a single non-fixed variable */
         if( setpartcons->nfixedzeros == setpartcons->nvars - 1 )
         {
            VAR* var;
            Bool found;
            int v;

            /* search unfixed variable */
            found = FALSE;
            for( v = 0; v < setpartcons->nvars && !found; ++v )
            {
               var = setpartcons->vars[v];
               found = !SCIPisZero(scip, SCIPvarGetUbGlobal(var));
            }
            assert(found);
            debugMessage("set partitioning constraint <%s>: fix <%s> == 1\n", SCIPconsGetName(cons), SCIPvarGetName(var));
            CHECK_OKAY( SCIPfixVar(scip, var, 1.0, &infeasible) );
            assert(!infeasible);
            CHECK_OKAY( SCIPdelCons(scip, cons) );
            (*nfixedvars)++;
            (*ndelconss)++;
            *result = SCIP_SUCCESS;
            continue;
         }

         /* aggregate variable and delete constraint, if constraint consists only of two non-fixed variables */
         if( setpartcons->nfixedzeros == setpartcons->nvars - 2 )
         {
            VAR* var;
            VAR* var1;
            VAR* var2;
            int v;

            /* search unfixed variable */
            var1 = NULL;
            var2 = NULL;
            for( v = 0; v < setpartcons->nvars && var2 == NULL; ++v )
            {
               var = setpartcons->vars[v];
               if( !SCIPisZero(scip, SCIPvarGetUbGlobal(var)) )
               {
                  if( var1 == NULL )
                     var1 = var;
                  else
                     var2 = var;
               }
            }
            assert(var1 != NULL && var2 != NULL);
            if( SCIPvarGetStatus(var1) != SCIP_VARSTATUS_AGGREGATED )
            {
               debugMessage("set partitioning constraint <%s>: aggregate <%s> == 1 - <%s>\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var1), SCIPvarGetName(var2));
               CHECK_OKAY( SCIPaggregateVar(scip, var1, var2, -1.0, 1.0, &infeasible) );
               if( infeasible )
               {
                  debugMessage("set partitioning constraint <%s>: infeasible aggregation <%s> == 1 - <%s>\n",
                     SCIPconsGetName(cons), SCIPvarGetName(var1), SCIPvarGetName(var2));
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }
               CHECK_OKAY( SCIPdelCons(scip, cons) );
               (*naggrvars)++;
               (*ndelconss)++;
               *result = SCIP_SUCCESS;
               continue;
            }
            else if( SCIPvarGetStatus(var2) != SCIP_VARSTATUS_AGGREGATED )
            {
               debugMessage("set partitioning constraint <%s>: aggregate <%s> == 1 - <%s>\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var2), SCIPvarGetName(var1));
               CHECK_OKAY( SCIPaggregateVar(scip, var2, var1, -1.0, 1.0, &infeasible) );
               if( infeasible )
               {
                  debugMessage("set partitioning constraint <%s>: infeasible aggregation <%s> == 1 - <%s>\n",
                     SCIPconsGetName(cons), SCIPvarGetName(var1), SCIPvarGetName(var2));
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }
               CHECK_OKAY( SCIPdelCons(scip, cons) );
               (*naggrvars)++;
               (*ndelconss)++;
               *result = SCIP_SUCCESS;
               continue;
            }
         }
      }

      setpartcons->changed = FALSE;
   }
   
   return SCIP_OKAY;
}




/*
 * variable usage counting
 */

static
DECL_CONSENABLE(consEnableSetpart)
{
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* consdata;
   SETPARTCONS* setpartcons;
   int v;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   setpartcons = consdata->setpartcons;
   assert(setpartcons != NULL);

   debugMessage("enabling information method of set partitioning constraint handler\n");

   /* increase the number of uses for each variable in the constraint */
   for( v = 0; v < setpartcons->nvars; ++v )
   {
      CHECK_OKAY( conshdlrdataIncVaruses(scip, conshdlrdata, setpartcons->vars[v]) );
   }

   return SCIP_OKAY;
}

static
DECL_CONSDISABLE(consDisableSetpart)
{
   CONSHDLRDATA* conshdlrdata;
   INTARRAY* varuses;
   CONSDATA* consdata;
   SETPARTCONS* setpartcons;
   int v;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   varuses = conshdlrdata->varuses;
   assert(varuses != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   setpartcons = consdata->setpartcons;
   assert(setpartcons != NULL);

   debugMessage("disabling information method of set partitioning constraint handler\n");

   /* decrease the number of uses for each variable in the constraint */
   for( v = 0; v < setpartcons->nvars; ++v )
   {
      CHECK_OKAY( conshdlrdataDecVaruses(scip, conshdlrdata, setpartcons->vars[v]) );
   }

   return SCIP_OKAY;
}

static
DECL_LINCONSUPGD(linconsUpgdSetpart)
{
   assert(upgdcons != NULL);

   /* check, if linear constraint can be upgraded to set partitioning constraint
    * -> a set partitioning constraint consists only of binary variables with a coefficient of +1.0,
    *    and has left hand side of +1.0, and right hand side of +1.0: x(S) == 1.0
    */
   if( nposbin == nvars && ncoeffspone == nvars && SCIPisEQ(scip, lhs, 1.0) && SCIPisEQ(scip, rhs, 1.0) )
   {
      debugMessage("upgrading constraint <%s> to set partitioning constraint\n", SCIPconsGetName(cons));

      /* create the set partitioning constraint (an automatically upgraded constraint is always unmodifiable) */
      CHECK_OKAY( SCIPcreateConsSetpart(scip, upgdcons, SCIPconsGetName(cons), nvars, vars,
                     SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     local, FALSE, removeable) );
   }

   return SCIP_OKAY;
}




/*
 * Callback methods of event handler
 */

static
DECL_EVENTEXEC(eventExecSetpart)
{
   SETPARTCONS* setpartcons;
   EVENTTYPE eventtype;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   debugMessage("Exec method of bound change event handler for set partitioning constraints %p\n", eventdata);

   setpartcons = (SETPARTCONS*)eventdata;
   assert(setpartcons != NULL);

   eventtype = SCIPeventGetType(event);
   switch( eventtype )
   {
   case SCIP_EVENTTYPE_LBTIGHTENED:
      setpartcons->nfixedones++;
      break;
   case SCIP_EVENTTYPE_LBRELAXED:
      setpartcons->nfixedones--;
      break;
   case SCIP_EVENTTYPE_UBTIGHTENED:
      setpartcons->nfixedzeros++;
      break;
   case SCIP_EVENTTYPE_UBRELAXED:
      setpartcons->nfixedzeros--;
      break;
   default:
      errorMessage("invalid event type");
      abort();
   }
   assert(0 <= setpartcons->nfixedzeros && setpartcons->nfixedzeros <= setpartcons->nvars);
   assert(0 <= setpartcons->nfixedones && setpartcons->nfixedones <= setpartcons->nvars);

   setpartcons->changed = TRUE;

   debugMessage(" -> constraint has %d zero-fixed and %d one-fixed of %d variables\n", 
      setpartcons->nfixedzeros, setpartcons->nfixedones, setpartcons->nvars);

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for set partitioning constraint and includes it in SCIP */
RETCODE SCIPincludeConsHdlrSetpart(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create event handler for bound change events */
   CHECK_OKAY( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
                  NULL, NULL, NULL,
                  NULL, eventExecSetpart,
                  NULL) );

   /* create constraint handler data */
   CHECK_OKAY( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConsHdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
                  CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ,
                  CONSHDLR_NEEDSCONS,
                  consFreeSetpart, NULL, NULL,
                  consDeleteSetpart, consTransSetpart, 
                  consSepaSetpart, consEnfolpSetpart, consEnfopsSetpart, consCheckSetpart, NULL, consPresolSetpart,
                  consEnableSetpart, consDisableSetpart,
                  conshdlrdata) );

   /* include the linear constraint to set partitioning constraint upgrade in the linear constraint handler */
   CHECK_OKAY( SCIPincludeLinconsUpgrade(scip, linconsUpgdSetpart, LINCONSUPGD_PRIORITY) );

   /* set partitioning constraint handler parameters */
   CHECK_OKAY( SCIPaddIntParam(scip,
                  "conshdlr/setpart/npseudobranches", 
                  "number of children created in pseudo branching",
                  &conshdlrdata->npseudobranches, DEFAULT_NPSEUDOBRANCHES, 2, INT_MAX, NULL, NULL) );
   
   return SCIP_OKAY;
}

/** creates and captures a set partitioning constraint */
RETCODE SCIPcreateConsSetpart(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             local,              /**< is set partitioning constraint only valid locally? */
   Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   assert(scip != NULL);

   /* find the set partitioning constraint handler */
   conshdlr = SCIPfindConsHdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("set partitioning constraint handler not found");
      return SCIP_INVALIDCALL;
   }

   /* create the constraint specific data */
   CHECK_OKAY( SCIPallocBlockMemory(scip, &consdata) );
   if( SCIPstage(scip) == SCIP_STAGE_PROBLEM )
   {
      if( local )
      {
         errorMessage("problem constraint cannot be local");
         return SCIP_INVALIDDATA;
      }

      /* create constraint in original problem */
      CHECK_OKAY( setpartconsCreate(scip, &consdata->setpartcons, nvars, vars, modifiable, removeable) );
   }
   else
   {
      /* create constraint in transformed problem */
      CHECK_OKAY( setpartconsCreateTransformed(scip, &consdata->setpartcons, nvars, vars, 
                     local, modifiable, removeable) );
   }
   consdata->row = NULL;

   /* create constraint (propagation is never used for set partitioning constraints) */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, separate, enforce, check, FALSE) );

   return SCIP_OKAY;
}


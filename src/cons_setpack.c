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

/**@file   cons_setpack.c
 * @brief  constraint handler for set packing constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_setpack.h"
#include "cons_setcover.h"
#include "cons_linear.h"


#define CONSHDLR_NAME          "setpack"
#define CONSHDLR_DESC          "set packing constraint"
#define CONSHDLR_SEPAPRIORITY   +700100
#define CONSHDLR_ENFOPRIORITY   +700100
#define CONSHDLR_CHECKPRIORITY  -700100
#define CONSHDLR_SEPAFREQ             4
#define CONSHDLR_PROPFREQ            -1
#define CONSHDLR_NEEDSCONS         TRUE

#define EVENTHDLR_NAME         "setpack"
#define EVENTHDLR_DESC         "bound change event handler for set packing constraints"

#define LINCONSUPGD_PRIORITY    +700100

#define MINBRANCHWEIGHT               0.4  /**< minimum weight of both sets in set packing branching */
#define MAXBRANCHWEIGHT               0.8  /**< maximum weight of both sets in set packing branching */


/** constraint handler data */
struct ConsHdlrData
{
   INTARRAY*        varuses;            /**< number of times a variable is used in the active set packing constraints */
};

/** set packing constraint data */
struct SetpackCons
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
typedef struct SetpackCons SETPACKCONS; /**< set packing constraint data */

/** constraint data for set packing constraints */
struct ConsData
{
   SETPACKCONS*     setpackcons;        /**< set packing constraint data */
   ROW*             row;                /**< LP row, if constraint is already stored in LP row format */
};




/*
 * Local methods
 */

/** creates constaint handler data for set packing constraint handler */
static
RETCODE conshdlrdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);

   CHECK_OKAY( SCIPallocMemory(scip, conshdlrdata) );
   CHECK_OKAY( SCIPcreateIntarray(scip, &(*conshdlrdata)->varuses) );

   return SCIP_OKAY;
}

/** frees constraint handler data for set packing constraint handler */
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

   debugMessage("varuses of <%s>: %d\n", SCIPvarGetName(var), SCIPgetIntarrayVal(scip, varuses, SCIPvarGetIndex(var)));

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

   debugMessage("varuses of <%s>: %d\n", SCIPvarGetName(var), SCIPgetIntarrayVal(scip, varuses, SCIPvarGetIndex(var)));

   return SCIP_OKAY;
}

/** creates event data for variable at given position, and catches events */
static
RETCODE setpackconsCatchEvent(
   SCIP*            scip,               /**< SCIP data structure */
   SETPACKCONS*     setpackcons,        /**< set packing constraint object */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< array position of variable to catch bound change events for */
   )
{
   VAR* var;

   assert(setpackcons != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < setpackcons->nvars);
   assert(setpackcons->vars != NULL);

   var = setpackcons->vars[pos];
   assert(var != NULL);

   /* catch bound change events on variables */
   CHECK_OKAY( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, (EVENTDATA*)setpackcons) );
   
   /* update the fixed variables counters for this variable */
   if( SCIPisEQ(scip, SCIPvarGetUbLocal(var), 0.0) )
      setpackcons->nfixedzeros++;
   else if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), 1.0) )
      setpackcons->nfixedones++;

   return SCIP_OKAY;
}

/** deletes event data for variable at given position, and drops events */
static
RETCODE setpackconsDropEvent(
   SCIP*            scip,               /**< SCIP data structure */
   SETPACKCONS*     setpackcons,        /**< set packing constraint object */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< array position of variable to catch bound change events for */
   )
{
   VAR* var;

   assert(setpackcons != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < setpackcons->nvars);
   assert(setpackcons->vars != NULL);

   var = setpackcons->vars[pos];
   assert(var != NULL);
   
   CHECK_OKAY( SCIPdropVarEvent(scip, var, eventhdlr, (EVENTDATA*)setpackcons) );

   /* update the fixed variables counters for this variable */
   if( SCIPisEQ(scip, SCIPvarGetUbLocal(var), 0.0) )
      setpackcons->nfixedzeros--;
   else if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), 1.0) )
      setpackcons->nfixedones--;

   return SCIP_OKAY;
}

/** catches bound change events and locks rounding for variable at given position in transformed set packing constraint */
static
RETCODE setpackconsLockCoef(
   SCIP*            scip,               /**< SCIP data structure */
   SETPACKCONS*     setpackcons,        /**< set packing constraint object */
   EVENTHDLR*       eventhdlr,          /**< event handler for bound change events, or NULL */
   int              pos                 /**< position of variable in set packing constraint */
   )
{
   VAR* var;
      
   assert(scip != NULL);
   assert(setpackcons != NULL);
   assert(setpackcons->transformed);
   assert(0 <= pos && pos < setpackcons->nvars);

   var = setpackcons->vars[pos];
   
   debugMessage("locking coefficient <%s> in set packing constraint\n", SCIPvarGetName(var));

   if( eventhdlr == NULL )
   {
      /* get event handler for updating set packing constraint activity bounds */
      eventhdlr = SCIPfindEventHdlr(scip, EVENTHDLR_NAME);
      if( eventhdlr == NULL )
      {
         errorMessage("event handler for set packing constraints not found");
         return SCIP_PLUGINNOTFOUND;
      }
   }

   /* catch bound change events on variable */
   assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL);
   CHECK_OKAY( setpackconsCatchEvent(scip, setpackcons, eventhdlr, pos) );

   /* forbid rounding of variable */
   if( !setpackcons->local )
      SCIPvarForbidRoundUp(var);

   return SCIP_OKAY;
}

/** drops bound change events and unlocks rounding for variable at given position in transformed set packing constraint */
static
RETCODE setpackconsUnlockCoef(
   SCIP*            scip,               /**< SCIP data structure */
   SETPACKCONS*     setpackcons,        /**< set packing constraint object */
   EVENTHDLR*       eventhdlr,          /**< event handler for bound change events, or NULL */
   int              pos                 /**< position of variable in set packing constraint */
   )
{
   VAR* var;

   assert(scip != NULL);
   assert(setpackcons != NULL);
   assert(setpackcons->transformed);
   assert(0 <= pos && pos < setpackcons->nvars);

   var = setpackcons->vars[pos];

   debugMessage("unlocking coefficient <%s> in set packing constraint\n", SCIPvarGetName(var));

   if( eventhdlr == NULL )
   {
      /* get event handler for updating set packing constraint activity bounds */
      eventhdlr = SCIPfindEventHdlr(scip, EVENTHDLR_NAME);
      if( eventhdlr == NULL )
      {
         errorMessage("event handler for set packing constraints not found");
         return SCIP_PLUGINNOTFOUND;
      }
   }
   
   /* drop bound change events on variable */
   assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL);
   CHECK_OKAY( setpackconsDropEvent(scip, setpackcons, eventhdlr, pos) );

   /* allow rounding of variable */
   if( !setpackcons->local )
      SCIPvarAllowRoundUp(var);

   return SCIP_OKAY;
}

/** catches bound change events and locks rounding for all variables in transformed set packing constraint */
static
RETCODE setpackconsLockAllCoefs(
   SCIP*            scip,               /**< SCIP data structure */
   SETPACKCONS*     setpackcons        /**< set packing constraint object */
   )
{
   EVENTHDLR* eventhdlr;
   int i;

   assert(scip != NULL);
   assert(setpackcons != NULL);
   assert(setpackcons->transformed);

   /* get event handler for updating set packing constraint activity bounds */
   eventhdlr = SCIPfindEventHdlr(scip, EVENTHDLR_NAME);
   if( eventhdlr == NULL )
   {
      errorMessage("event handler for set packing constraints not found");
      return SCIP_PLUGINNOTFOUND;
   }

   /* lock every single coefficient */
   for( i = 0; i < setpackcons->nvars; ++i )
   {
      CHECK_OKAY( setpackconsLockCoef(scip, setpackcons, eventhdlr, i) );
   }

   return SCIP_OKAY;
}

/** drops bound change events and unlocks rounding for all variables in transformed set packing constraint */
static
RETCODE setpackconsUnlockAllCoefs(
   SCIP*            scip,               /**< SCIP data structure */
   SETPACKCONS*     setpackcons        /**< set packing constraint object */
   )
{
   EVENTHDLR* eventhdlr;
   int i;

   assert(scip != NULL);
   assert(setpackcons != NULL);
   assert(setpackcons->transformed);

   /* get event handler for updating set packing constraint activity bounds */
   eventhdlr = SCIPfindEventHdlr(scip, EVENTHDLR_NAME);
   if( eventhdlr == NULL )
   {
      errorMessage("event handler for set packing constraints not found");
      return SCIP_PLUGINNOTFOUND;
   }

   /* unlock every single coefficient */
   for( i = 0; i < setpackcons->nvars; ++i )
   {
      CHECK_OKAY( setpackconsUnlockCoef(scip, setpackcons, eventhdlr, i) );
   }

   return SCIP_OKAY;
}

/** deletes coefficient at given position from set packing constraint object */
static
RETCODE setpackconsDelCoefPos(
   SCIP*            scip,               /**< SCIP data structure */
   SETPACKCONS*     setpackcons,        /**< set packing constraint object */
   int              pos                 /**< position of coefficient to delete */
   )
{
   VAR* var;

   assert(setpackcons != NULL);
   assert(0 <= pos && pos < setpackcons->nvars);

   var = setpackcons->vars[pos];
   assert(var != NULL);
   assert(setpackcons->transformed ^ (SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL));

   if( setpackcons->transformed )
   {
      /* drop bound change events and unlock the rounding of variable */
      CHECK_OKAY( setpackconsUnlockCoef(scip, setpackcons, NULL, pos) );
   }

   /* move the last variable to the free slot */
   setpackcons->vars[pos] = setpackcons->vars[setpackcons->nvars-1];
   setpackcons->nvars--;

   setpackcons->changed = TRUE;

   return SCIP_OKAY;
}

/** creates a set packing constraint data object */
static
RETCODE setpackconsCreate(
   SCIP*            scip,               /**< SCIP data structure */
   SETPACKCONS**    setpackcons,        /**< pointer to store the set packing constraint */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< variables of the constraint */
   Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   assert(setpackcons != NULL);
   assert(nvars == 0 || vars != NULL);

   CHECK_OKAY( SCIPallocBlockMemory(scip, setpackcons) );
   if( nvars > 0 )
   {
      VAR* var;
      int v;

      CHECK_OKAY( SCIPduplicateBlockMemoryArray(scip, &(*setpackcons)->vars, vars, nvars) );
      (*setpackcons)->varssize = nvars;
      (*setpackcons)->nvars = nvars;
   }
   else
   {
      (*setpackcons)->vars = NULL;
      (*setpackcons)->varssize = 0;
      (*setpackcons)->nvars = 0;
   }
   (*setpackcons)->nfixedzeros = 0;
   (*setpackcons)->nfixedones = 0;
   (*setpackcons)->local = FALSE;
   (*setpackcons)->modifiable = modifiable;
   (*setpackcons)->removeable = removeable;
   (*setpackcons)->transformed = FALSE;
   (*setpackcons)->changed = TRUE;

   return SCIP_OKAY;
}   

/** creates a transformed set packing constraint data object */
static
RETCODE setpackconsCreateTransformed(
   SCIP*            scip,               /**< SCIP data structure */
   SETPACKCONS**    setpackcons,        /**< pointer to store the set packing constraint */
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

   assert(setpackcons != NULL);
   assert(nvars == 0 || vars != NULL);

   CHECK_OKAY( setpackconsCreate(scip, setpackcons, nvars, vars, modifiable, removeable) );
   (*setpackcons)->local = local;
   (*setpackcons)->transformed = TRUE;

   /* transform the variables */
   for( i = 0; i < (*setpackcons)->nvars; ++i )
   {
      var = (*setpackcons)->vars[i];
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
         (*setpackcons)->vars[i] = var;
         assert(var != NULL);
      }
      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL);
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
   }

   /* catch bound change events and lock the rounding of variables */
   CHECK_OKAY( setpackconsLockAllCoefs(scip, *setpackcons) );

   return SCIP_OKAY;
}

/** frees a set packing constraint data */
static
RETCODE setpackconsFree(
   SCIP*            scip,               /**< SCIP data structure */
   SETPACKCONS**    setpackcons        /**< pointer to store the set packing constraint */
   )
{
   assert(setpackcons != NULL);
   assert(*setpackcons != NULL);

   if( (*setpackcons)->transformed )
   {
      /* drop bound change events and unlock the rounding of variables */
      CHECK_OKAY( setpackconsUnlockAllCoefs(scip, *setpackcons) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*setpackcons)->vars, (*setpackcons)->varssize);
   SCIPfreeBlockMemory(scip, setpackcons);

   return SCIP_OKAY;
}

/** creates an LP row from a set packing constraint data object */
static
RETCODE setpackconsToRow(
   SCIP*            scip,               /**< SCIP data structure */
   SETPACKCONS*     setpackcons,        /**< set packing constraint data */
   const char*      name,               /**< name of the constraint */
   ROW**            row                 /**< pointer to an LP row data object */
   )
{
   int v;

   assert(setpackcons != NULL);
   assert(row != NULL);

   CHECK_OKAY( SCIPcreateRow(scip, row, name, 0, NULL, NULL, -SCIPinfinity(scip), 1.0,
                  setpackcons->local, setpackcons->modifiable, setpackcons->removeable) );
   
   for( v = 0; v < setpackcons->nvars; ++v )
   {
      CHECK_OKAY( SCIPaddVarToRow(scip, *row, setpackcons->vars[v], 1.0) );
   }

   return SCIP_OKAY;
}

/** checks constraint for violation only looking at the fixed variables, applies further fixings if possible */
static
RETCODE processFixings(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< set packing constraint to be separated */
   Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   Bool*            reduceddom,         /**< pointer to store whether a domain reduction was found */
   Bool*            addcut,             /**< pointer to store whether this constraint must be added as a cut */
   Bool*            mustcheck           /**< pointer to store whether this constraint must be checked for feasibility */
   )
{
   CONSDATA* consdata;
   SETPACKCONS* setpackcons;

   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(reduceddom != NULL);
   assert(addcut != NULL);
   assert(mustcheck != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   setpackcons = consdata->setpackcons;
   assert(setpackcons != NULL);
   assert(setpackcons->nvars == 0 || setpackcons->vars != NULL);
   assert(0 <= setpackcons->nfixedzeros && setpackcons->nfixedzeros <= setpackcons->nvars);
   assert(0 <= setpackcons->nfixedones && setpackcons->nfixedones <= setpackcons->nvars);

   *addcut = FALSE;
   *mustcheck = FALSE;

   if( setpackcons->nfixedones >= 2 )
   {
      /* the constraint cannot be feasible, because two variables are already fixed to one */
      CHECK_OKAY( SCIPresetConsAge(scip, cons) );
      *cutoff = TRUE;
   }
   else if( setpackcons->nfixedzeros >= setpackcons->nvars - 1 )
   {
      /* constraint is feasible anyway, because all except one of its variables are fixed to zero
       * -> if constraint is not modifiable (i.e. if not more variables may appear), the constraint can be disabled
       */
      if( !setpackcons->modifiable )
      {
         CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
      }
      else
      {
         CHECK_OKAY( SCIPincConsAge(scip, cons) );
      }
   }
   else if( setpackcons->nfixedones == 1 && setpackcons->nfixedzeros == setpackcons->nvars - 1 )
   {
      /* all variables are fixed, and the constraint is satisfied -> disable constraint if it is not modifiable */
      if( !setpackcons->modifiable )
      {
         CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
      }
   }
   else if( setpackcons->nfixedones == 1 )
   {
      VAR** vars;
      VAR* var;
      Bool fixedonefound;
      Bool fixed;
      int nvars;
      int v;

      assert(setpackcons->nfixedzeros < setpackcons->nvars - 1);

      /* all other variables except the fixed one must be zero, and if the constraint is not modifiable (i.e. if not
       * more variables may appear), it can be disabled after the other variables are fixed to zero
       */
      
      /* check all variables, if fixings to zero can be applied */
      vars = setpackcons->vars;
      nvars = setpackcons->nvars;
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

      if( setpackcons->modifiable )
      {
         CHECK_OKAY( SCIPresetConsAge(scip, cons) );
      }
      else
      {
         CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
      }
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
   SETPACKCONS*     setpackcons,        /**< set packing constraint to be checked */
   SOL*             sol                 /**< primal CIP solution */
   )
{
   VAR** vars;
   Real solval;
   Real sum;
   int nvars;
   int v;
   
   /* calculate the constraint's activity */
   vars = setpackcons->vars;
   nvars = setpackcons->nvars;
   sum = 0.0;
   assert(SCIPfeastol(scip) < 0.1); /* to make the comparison against 1.1 working */
   for( v = 0; v < nvars && sum < 1.1; ++v )  /* if sum >= 1.1, it is clearly too large */
   {
      assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY);
      solval = SCIPgetSolVal(scip, sol, vars[v]);
      assert(SCIPisFeasGE(scip, solval, 0.0) && SCIPisFeasLE(scip, solval, 1.0));
      sum += solval;
   }
   return SCIPisFeasLE(scip, sum, 1.0);
}

/** checks constraint for violation, and adds it as a cut if possible */
static
RETCODE separate(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< set packing constraint to be separated */
   Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   Bool*            separated,          /**< pointer to store whether a cut was found */
   Bool*            reduceddom          /**< pointer to store whether a domain reduction was found */
   )
{
   CONSDATA* consdata;
   SETPACKCONS* setpackcons;
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
   setpackcons = consdata->setpackcons;
   assert(setpackcons != NULL);
   assert(setpackcons->nvars == 0 || setpackcons->vars != NULL);
   assert(0 <= setpackcons->nfixedzeros && setpackcons->nfixedzeros <= setpackcons->nvars);
   assert(0 <= setpackcons->nfixedones && setpackcons->nfixedones <= setpackcons->nvars);

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
         addcut = !check(scip, setpackcons, NULL);

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
         /* convert set packing constraint data into LP row */
         CHECK_OKAY( setpackconsToRow(scip, setpackcons, SCIPconsGetName(cons), &consdata->row) );
      }
      assert(consdata->row != NULL);
      assert(!SCIProwIsInLP(consdata->row));
            
      /* insert LP row as cut */
      CHECK_OKAY( SCIPaddCut(scip, consdata->row, 1.0/(setpackcons->nvars+1)) );
      CHECK_OKAY( SCIPresetConsAge(scip, cons) );
      *separated = TRUE;
   }

   return SCIP_OKAY;
}

/** enforces the pseudo solution on the given constraint */
static
RETCODE enforcePseudo(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< set packing constraint to be separated */
   Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   Bool*            infeasible,         /**< pointer to store whether the constraint was infeasible */
   Bool*            reduceddom,         /**< pointer to store whether a domain reduction was found */
   Bool*            solvelp             /**< pointer to store whether the LP has to be solved */
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
      SETPACKCONS* setpackcons;

      assert(!addcut);

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);
      setpackcons = consdata->setpackcons;
      assert(setpackcons != NULL);

      *infeasible = !check(scip, setpackcons, NULL);

      if( !infeasible )
      {
         /* constraint was feasible -> increase age */
         CHECK_OKAY( SCIPincConsAge(scip, cons) );
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
DECL_CONSFREE(consFreeSetpack)
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
DECL_CONSDELETE(consDeleteSetpack)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* free LP row and setpack constraint */
   if( (*consdata)->row != NULL )
   {
      CHECK_OKAY( SCIPreleaseRow(scip, &(*consdata)->row) );
   }
   CHECK_OKAY( setpackconsFree(scip, &(*consdata)->setpackcons) );

   /* free constraint data object */
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

static
DECL_CONSTRANS(consTransSetpack)
{
   CONSDATA* sourcedata;
   CONSDATA* targetdata;
   SETPACKCONS* setpackcons;

   /*debugMessage("Trans method of setpack constraints\n");*/

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPstage(scip) == SCIP_STAGE_INITSOLVE);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->row == NULL);  /* in original problem, there cannot be LP rows */
   assert(sourcedata->setpackcons != NULL);

   /* create constraint data for target constraint */
   CHECK_OKAY( SCIPallocBlockMemory(scip, &targetdata) );

   setpackcons = sourcedata->setpackcons;

   CHECK_OKAY( setpackconsCreateTransformed(scip, &targetdata->setpackcons, setpackcons->nvars, setpackcons->vars,
                  setpackcons->local, setpackcons->modifiable, setpackcons->removeable) );
   targetdata->row = NULL;

   /* create target constraint */
   assert(!SCIPconsIsPropagated(sourcecons));
   CHECK_OKAY( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
                  SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
                  SCIPconsIsPropagated(sourcecons)) );

   return SCIP_OKAY;
}

static
DECL_CONSSEPA(consSepaSetpack)
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

   debugMessage("separating %d/%d set packing constraints\n", nusefulconss, nconss);

   *result = SCIP_DIDNOTFIND;

   cutoff = FALSE;
   separated = FALSE;
   reduceddom = FALSE;

   /* step 1: check all useful set packing constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff && !reduceddom; ++c )
   {
      CHECK_OKAY( separate(scip, conss[c], &cutoff, &separated, &reduceddom) );
   }

   /* step 2: combine set packing constraints to get more cuts */
   todoMessage("further cuts of set packing constraints");

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

/** if fractional variables exist, chooses a set S of them and branches on (i) x(S) == 0, and (ii) x(S) >= 1 */
static
RETCODE branchLP(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLR*        conshdlr,           /**< set packing constraint handler */
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

   todoMessage("use a better set packing branching on LP solution");

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
   
   /* sort fractional variables by number of uses in enabled set packing constraints */
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

   /* if none of the fractional variables is member of a set packing constraint,
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

         /* perform the set packing branching on the selected variables */
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
            sprintf(name, "SPB%lld", SCIPgetNodenum(scip));

            CHECK_OKAY( SCIPcreateConsSetcover(scip, &newcons, name, nselcands, sortcands,
                           TRUE, TRUE, FALSE, TRUE, FALSE, TRUE) );
            CHECK_OKAY( SCIPaddConsNode(scip, node, newcons) );
            CHECK_OKAY( SCIPreleaseCons(scip, &newcons) );
         }
      
         *result = SCIP_BRANCHED;
         
#ifdef DEBUG
         debugMessage("set packing branching: nselcands=%d/%d, weight(S)=%g, A={", nselcands, nlpcands, branchweight);
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

/** if unfixed variables exist, chooses a set S of them and creates |S|+1 child nodes:
 *   - for each variable i from S, create child node with x_0 = ... = x_i-1 = 0, x_i = 1
 *   - create an additional child node x_0 = ... = x_n-1 = 0
 */
static
RETCODE branchPseudo(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLR*        conshdlr,           /**< set packing constraint handler */
   RESULT*          result              /**< pointer to store the result SCIP_BRANCHED, if branching was applied */
   )
{
   CONSHDLRDATA* conshdlrdata;
   INTARRAY* varuses;
   VAR** pseudocands;
   VAR** sortcands;
   VAR* var;
   NODE* node;
   int npseudocands;
   int nsortcands;
   int nbranchcands;
   int uses;
   int i;
   int j;

   todoMessage("use a better set packing branching on pseudo solution");

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

   /* get temporary memory */
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &sortcands, npseudocands) );
   
   /* sort fractional variables by number of uses in enabled set packing constraints */
   nsortcands = 0;
   for( i = 0; i < npseudocands; ++i )
   {
      var = pseudocands[i];
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
   assert(nsortcands <= npseudocands);

   if( nsortcands == 0 )
   {
      /* none of the unfixed variables is member of a set packing constraint
       * -> we are not responsible for doing the branching
       */
      return SCIP_OKAY;
   }
   
   /* branch on the first part of the sorted candidates:
    * - for each of these variables i, create a child node x_0 = ... = x_i-1 = 0, x_i = 1
    * - create an additional child node x_0 = ... = x_n-1 = 0
    */
   nbranchcands = (nsortcands+9)/10;
   assert(nbranchcands >= 1);
   for( i = 0; i < nbranchcands; ++i )
   {            
      /* create child with x_0 = ... = x_i-1 = 0, x_i = 1 */
      CHECK_OKAY( SCIPcreateChild(scip, &node) );
      for( j = 0; j < i; ++j )
      {
         CHECK_OKAY( SCIPchgVarUbNode(scip, node, sortcands[j], 0.0) );
      }
      CHECK_OKAY( SCIPchgVarLbNode(scip, node, sortcands[i], 1.0) );
   }
   /* create child with x_0 = ... = x_n = 0 */
   CHECK_OKAY( SCIPcreateChild(scip, &node) );
   for( i = 0; i < nbranchcands; ++i )
   {
      CHECK_OKAY( SCIPchgVarUbNode(scip, node, sortcands[i], 0.0) );
   }

   *result = SCIP_BRANCHED;

#ifdef DEBUG
   {
      int nchildren;
      CHECK_OKAY( SCIPgetChildren(scip, NULL, &nchildren) );
      debugMessage("branched on set cover constraint in pseudo solution: %d children\n", nchildren);
   }
#endif

   /* free temporary memory */
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &sortcands) );

   return SCIP_OKAY;
}



static
DECL_CONSENFOLP(consEnfolpSetpack)
{
   Bool cutoff;
   Bool separated;
   Bool reduceddom;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   debugMessage("enforcing %d set packing constraints\n", nconss);

   *result = SCIP_FEASIBLE;

   cutoff = FALSE;
   separated = FALSE;
   reduceddom = FALSE;

   /* step 1: check all useful set packing constraints for feasibility */
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
      
      /* step 3: check all obsolete set packing constraints for feasibility */
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
DECL_CONSENFOPS(consEnfopsSetpack)
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

   *result = SCIP_FEASIBLE;

   cutoff = FALSE;
   infeasible = FALSE;
   reduceddom = FALSE;
   solvelp = FALSE;

   /* check all set packing constraints for feasibility */
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
DECL_CONSCHECK(consCheckSetpack)
{
   CONSDATA* consdata;
   CONS* cons;
   SETPACKCONS* setpackcons;
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

   /* check all set packing constraints for feasibility */
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);
      setpackcons = consdata->setpackcons;
      assert(setpackcons != NULL);
      if( checklprows || consdata->row == NULL || !SCIProwIsInLP(consdata->row) )
      {
         if( !check(scip, setpackcons, sol) )
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
RETCODE setpackconsApplyFixings(
   SCIP*            scip,               /**< SCIP data structure */
   SETPACKCONS*     setpackcons         /**< set packing constraint object */
   )
{
   VAR* var;
   int v;

   assert(setpackcons != NULL);

   if( setpackcons->nfixedzeros >= 1 )
   {
      assert(setpackcons->vars != NULL);

      v = 0;
      while( v < setpackcons->nvars )
      {
         var = setpackcons->vars[v];
         if( SCIPisZero(scip, SCIPvarGetUbGlobal(var)) )
         {
            CHECK_OKAY( setpackconsDelCoefPos(scip, setpackcons, v) );
         }
         else
            ++v;
      }
   }

   return SCIP_OKAY;
}

static
DECL_CONSPRESOL(consPresolSetpack)
{
   CONS* cons;
   CONSDATA* consdata;
   SETPACKCONS* setpackcons;
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
      setpackcons = consdata->setpackcons;
      assert(setpackcons != NULL);

      if( !setpackcons->changed )
         continue;

      debugMessage("presolving set packing constraint <%s>\n", SCIPconsGetName(cons));

      /* remove all variables that are fixed to zero */
      CHECK_OKAY( setpackconsApplyFixings(scip, setpackcons) );

      /* check for infeasibility due to more than one variable fixed to one */
      if( setpackcons->nfixedones >= 2 )
      {
         debugMessage("set packing constraint <%s> is infeasible\n", SCIPconsGetName(cons));
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      /* check, if we can fix all other variables to zero, because a variable fixed to one exists
       *  -> if the constraint is not modifiable, we can delete it (which is done later)
       */
      if( setpackcons->nfixedones == 1 )
      {
         VAR* var;
         int v;

         debugMessage("set packing constraint <%s> has a variable fixed to 1.0\n", SCIPconsGetName(cons));
         for( v = 0; v < setpackcons->nvars; ++v )
         {
            var = setpackcons->vars[v];
            if( SCIPisZero(scip, SCIPvarGetLbGlobal(var)) && !SCIPisZero(scip, SCIPvarGetUbGlobal(var)) )
            {
               CHECK_OKAY( SCIPfixVar(scip, var, 0.0) );
               (*nfixedvars)++;
               *result = SCIP_SUCCESS;
            }
         }
      }

      if( !setpackcons->modifiable )
      {
         /* constaint is redundant, if it consists of not more than one non-fixed variable */
         if( setpackcons->nfixedzeros >= setpackcons->nvars - 1 )
         {
            debugMessage("set packing constraint <%s> is redundant\n", SCIPconsGetName(cons));
            CHECK_OKAY( SCIPdelCons(scip, cons) );
            (*ndelconss)++;
            *result = SCIP_SUCCESS;
            continue;
         }
      }

      setpackcons->changed = FALSE;
   }
   
   return SCIP_OKAY;
}




/*
 * variable usage counting
 */

static
DECL_CONSENABLE(consEnableSetpack)
{
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* consdata;
   SETPACKCONS* setpackcons;
   int v;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   setpackcons = consdata->setpackcons;
   assert(setpackcons != NULL);

   debugMessage("enabling information method of set packing constraint handler\n");

   /* increase the number of uses for each variable in the constraint */
   for( v = 0; v < setpackcons->nvars; ++v )
   {
      CHECK_OKAY( conshdlrdataIncVaruses(scip, conshdlrdata, setpackcons->vars[v]) );
   }

   return SCIP_OKAY;
}

static
DECL_CONSDISABLE(consDisableSetpack)
{
   CONSHDLRDATA* conshdlrdata;
   INTARRAY* varuses;
   CONSDATA* consdata;
   SETPACKCONS* setpackcons;
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
   setpackcons = consdata->setpackcons;
   assert(setpackcons != NULL);

   debugMessage("disabling information method of set packing constraint handler\n");

   /* decrease the number of uses for each variable in the constraint */
   for( v = 0; v < setpackcons->nvars; ++v )
   {
      CHECK_OKAY( conshdlrdataDecVaruses(scip, conshdlrdata, setpackcons->vars[v]) );
   }

   return SCIP_OKAY;
}

static
DECL_LINCONSUPGD(linconsUpgdSetpack)
{
   assert(upgdcons != NULL);

   /* check, if linear constraint can be upgraded to set packing constraint
    * -> a set packing constraint consists only of binary variables with a coefficient of +1.0,
    *    and has infinite left hand side, and right hand side of +1.0: x(S) <= 1.0
    */
   if( nposbin == nvars && ncoeffspone == nvars && SCIPisInfinity(scip, -lhs) && SCIPisEQ(scip, rhs, 1.0) )
   {
      debugMessage("upgrading constraint <%s> to set packing constraint\n", SCIPconsGetName(cons));

      /* create the set packing constraint (an automatically upgraded constraint is always unmodifiable) */
      CHECK_OKAY( SCIPcreateConsSetpack(scip, upgdcons, SCIPconsGetName(cons), nvars, vars,
                     SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     local, FALSE, removeable) );
   }

   return SCIP_OKAY;
}




/*
 * Callback methods of event handler
 */

static
DECL_EVENTEXEC(eventExecSetpack)
{
   SETPACKCONS* setpackcons;
   EVENTTYPE eventtype;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   debugMessage("Exec method of bound change event handler for set packing constraints %p\n", eventdata);

   setpackcons = (SETPACKCONS*)eventdata;
   assert(setpackcons != NULL);

   eventtype = SCIPeventGetType(event);
   switch( eventtype )
   {
   case SCIP_EVENTTYPE_LBTIGHTENED:
      setpackcons->nfixedones++;
      break;
   case SCIP_EVENTTYPE_LBRELAXED:
      setpackcons->nfixedones--;
      break;
   case SCIP_EVENTTYPE_UBTIGHTENED:
      setpackcons->nfixedzeros++;
      break;
   case SCIP_EVENTTYPE_UBRELAXED:
      setpackcons->nfixedzeros--;
      break;
   default:
      errorMessage("invalid event type");
      abort();
   }
   assert(0 <= setpackcons->nfixedzeros && setpackcons->nfixedzeros <= setpackcons->nvars);
   assert(0 <= setpackcons->nfixedones && setpackcons->nfixedones <= setpackcons->nvars);

   setpackcons->changed = TRUE;

   debugMessage(" -> constraint has %d zero-fixed and %d one-fixed of %d variables\n", 
      setpackcons->nfixedzeros, setpackcons->nfixedones, setpackcons->nvars);

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for set packing constraint and includes it in SCIP */
RETCODE SCIPincludeConsHdlrSetpack(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create event handler for bound change events */
   CHECK_OKAY( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
                  NULL, NULL, NULL,
                  NULL, eventExecSetpack,
                  NULL) );

   /* create constraint handler data */
   CHECK_OKAY( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConsHdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
                  CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ,
                  CONSHDLR_NEEDSCONS,
                  consFreeSetpack, NULL, NULL,
                  consDeleteSetpack, consTransSetpack, 
                  consSepaSetpack, consEnfolpSetpack, consEnfopsSetpack, consCheckSetpack, NULL, consPresolSetpack,
                  consEnableSetpack, consDisableSetpack,
                  conshdlrdata) );

   /* include the linear constraint to set packing constraint upgrade in the linear constraint handler */
   CHECK_OKAY( SCIPincludeLinconsUpgrade(scip, linconsUpgdSetpack, LINCONSUPGD_PRIORITY) );

   return SCIP_OKAY;
}

/** creates and captures a set packing constraint */
RETCODE SCIPcreateConsSetpack(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             local,              /**< is set packing constraint only valid locally? */
   Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   assert(scip != NULL);

   /* find the set packing constraint handler */
   conshdlr = SCIPfindConsHdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("set packing constraint handler not found");
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
      CHECK_OKAY( setpackconsCreate(scip, &consdata->setpackcons, nvars, vars, modifiable, removeable) );
   }
   else
   {
      /* create constraint in transformed problem */
      CHECK_OKAY( setpackconsCreateTransformed(scip, &consdata->setpackcons, nvars, vars, 
                     local, modifiable, removeable) );
   }
   consdata->row = NULL;

   /* create constraint (propagation is never used for set packing constraints) */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, separate, enforce, check, FALSE) );

   return SCIP_OKAY;
}


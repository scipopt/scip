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
#pragma ident "@(#) $Id: sepa_intobj.c,v 1.18 2005/05/31 17:20:21 bzfpfend Exp $"

/**@file   sepa_intobj.c
 * @brief  integer objective value separator
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/sepa_intobj.h"


#define SEPA_NAME              "intobj"
#define SEPA_DESC              "integer objective value separator"
#define SEPA_PRIORITY              -100
#define SEPA_FREQ                    -1
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define EVENTHDLR_NAME         "intobj"
#define EVENTHDLR_DESC         "objective change event handler for integer objective value separator"




/*
 * Data structures
 */

/** separator data */
struct SepaData
{
   ROW*             objrow;             /**< objective value equality */
   VAR*             objvar;             /**< objective value variable */
};




/*
 * Local methods
 */

/** creates separator data */
static
RETCODE sepadataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   SEPADATA**       sepadata            /**< pointer to store separator data */
   )
{
   assert(sepadata != NULL);

   CHECK_OKAY( SCIPallocMemory(scip, sepadata) );
   (*sepadata)->objrow = NULL;
   (*sepadata)->objvar = NULL;

   return SCIP_OKAY;
}

/** frees separator data */
static
RETCODE sepadataFree(
   SCIP*            scip,               /**< SCIP data structure */
   SEPADATA**       sepadata            /**< pointer to separator data */
   )
{
   assert(sepadata != NULL);
   assert(*sepadata != NULL);
   assert((*sepadata)->objrow == NULL);
   assert((*sepadata)->objvar == NULL);

   SCIPfreeMemory(scip, sepadata);

   return SCIP_OKAY;
}

/** creates the objective value equality and the objective value variable, if not yet existing */
static
RETCODE createObjRow(
   SCIP*            scip,               /**< SCIP data structure */
   SEPADATA*        sepadata            /**< separator data */
   )
{
   assert(sepadata != NULL);

   if( sepadata->objrow == NULL )
   {
      VAR** vars;
      Real obj;
      int nvars;
      int v;

      /* create and add objective value variable */
      if( sepadata->objvar == NULL )
      {
         CHECK_OKAY( SCIPcreateVar(scip, &sepadata->objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
               SCIP_VARTYPE_IMPLINT, FALSE, TRUE, NULL, NULL, NULL, NULL) );
         CHECK_OKAY( SCIPaddVar(scip, sepadata->objvar) );
         CHECK_OKAY( SCIPaddVarLocks(scip, sepadata->objvar, +1, +1) );
      }

      /* get problem variables */
      vars = SCIPgetOrigVars(scip);
      nvars = SCIPgetNOrigVars(scip);

      /* create objective value equality */
      CHECK_OKAY( SCIPcreateEmptyRow(scip, &sepadata->objrow, "objrow", 0.0, 0.0,
                     FALSE, !SCIPallVarsInProb(scip), TRUE) );

      CHECK_OKAY( SCIPcacheRowExtensions(scip, sepadata->objrow) );
      for( v = 0; v < nvars; ++v )
      {
         obj = SCIPvarGetObj(vars[v]);
         if( !SCIPisZero(scip, obj) )
         {
            CHECK_OKAY( SCIPaddVarToRow(scip, sepadata->objrow, vars[v], obj) );
         }
      }
      CHECK_OKAY( SCIPaddVarToRow(scip, sepadata->objrow, sepadata->objvar, -1.0) );
      CHECK_OKAY( SCIPflushRowExtensions(scip, sepadata->objrow) );

      debugMessage("created objective value row: ");
      debug(SCIPprintRow(scip, sepadata->objrow, NULL));
   }

   return SCIP_OKAY;
}




/*
 * Callback methods of separator
 */

/** destructor of separator to free user data (called when SCIP is exiting) */
static
DECL_SEPAFREE(sepaFreeIntobj)
{  /*lint --e{715}*/
   SEPADATA* sepadata;

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   CHECK_OKAY( sepadataFree(scip, &sepadata) );

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** initialization method of separator (called after problem was transformed) */
#define sepaInitIntobj NULL


/** deinitialization method of separator (called before transformed problem is freed) */
static
DECL_SEPAEXIT(sepaExitIntobj)
{  /*lint --e{715}*/
   SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* release objective variable */
   if( sepadata->objvar != NULL )
   {
      CHECK_OKAY( SCIPreleaseVar(scip, &sepadata->objvar) );
   }

   return SCIP_OKAY;
}


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#define sepaInitsolIntobj NULL


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
DECL_SEPAEXITSOL(sepaExitsolIntobj)
{  /*lint --e{715}*/
   SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* release objective row */
   if( sepadata->objrow != NULL )
   {
      CHECK_OKAY( SCIPreleaseRow(scip, &sepadata->objrow) );
   }

   return SCIP_OKAY;
}


/** execution method of separator */
static
DECL_SEPAEXEC(sepaExecIntobj)
{  /*lint --e{715}*/
   SEPADATA* sepadata;
   Real objval;
   Real intobjval;
   Bool infeasible;
   Bool tightened;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* if the objective value may be fractional, we cannot do anything */
   if( !SCIPisObjIntegral(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* if the current objective value is integral, there is no integral objective value cut */
   objval = SCIPgetLPObjval(scip);
   objval = SCIPretransformObj(scip, objval);
   if( SCIPisFeasIntegral(scip, objval) )
      return SCIP_OKAY;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* the objective value is fractional: create the objective value equality, if not yet existing */
   CHECK_OKAY( createObjRow(scip, sepadata) );

   /* adjust the bounds of the objective value variable */
   if( SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE )
   {
      intobjval = SCIPceil(scip, objval);
      CHECK_OKAY( SCIPtightenVarLb(scip, sepadata->objvar, intobjval, &infeasible, &tightened) );
      debugMessage("new objective variable lower bound: <%s>[%g,%g]\n", 
         SCIPvarGetName(sepadata->objvar), SCIPvarGetLbLocal(sepadata->objvar), SCIPvarGetUbLocal(sepadata->objvar));
   }
   else
   {
      intobjval = SCIPfloor(scip, objval);
      CHECK_OKAY( SCIPtightenVarUb(scip, sepadata->objvar, intobjval, &infeasible, &tightened) );
      debugMessage("new objective variable upper bound: <%s>[%g,%g]\n", 
         SCIPvarGetName(sepadata->objvar), SCIPvarGetLbLocal(sepadata->objvar), SCIPvarGetUbLocal(sepadata->objvar));
   }

   /* add the objective value equality as a cut to the LP */
   if( infeasible )
      *result = SCIP_CUTOFF;
   else if( !SCIProwIsInLP(sepadata->objrow) )
   {
      CHECK_OKAY( SCIPaddCut(scip, sepadata->objrow, FALSE) );
      *result = SCIP_SEPARATED;
   }
   else if( tightened )
      *result = SCIP_REDUCEDDOM;
      
   return SCIP_OKAY;
}



/*
 * event handler for objective changes
 */

/** destructor of event handler to free user data (called when SCIP is exiting) */
#define eventFreeIntobj NULL

/** initialization method of event handler (called after problem was transformed) */
static
DECL_EVENTINIT(eventInitIntobj)
{  /*lint --e{715}*/
   CHECK_OKAY( SCIPcatchEvent(scip, SCIP_EVENTTYPE_VARADDED | SCIP_EVENTTYPE_OBJCHANGED, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
DECL_EVENTEXIT(eventExitIntobj)
{  /*lint --e{715}*/
   CHECK_OKAY( SCIPdropEvent(scip, SCIP_EVENTTYPE_VARADDED | SCIP_EVENTTYPE_OBJCHANGED, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
#define eventInitsolIntobj NULL

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
#define eventExitsolIntobj NULL

/** frees specific event data */
#define eventDeleteIntobj NULL

/** execution methode of objective change event handler */
static
DECL_EVENTEXEC(eventExecIntobj)
{  /*lint --e{715}*/
   EVENTHDLRDATA* eventhdlrdata;
   SEPADATA* sepadata;
   VAR* var;
   Real objdelta;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   sepadata = (SEPADATA*)eventhdlrdata;
   assert(sepadata != NULL);

   /* we don't have anything to do, if the objective value equality doesn't yet exist */
   if( sepadata->objrow == NULL )
      return SCIP_OKAY;

   var = SCIPeventGetVar(event);

   switch( SCIPeventGetType(event) )
   {
   case SCIP_EVENTTYPE_VARADDED:
      debugMessage("variable <%s> with obj=%g was added to the problem\n", SCIPvarGetName(var), SCIPvarGetObj(var));
      objdelta = SCIPvarGetObj(var);
      if( !SCIPisZero(scip, objdelta) )
      {
         CHECK_OKAY( SCIPaddVarToRow(scip, sepadata->objrow, var, SCIPvarGetObj(var)) );
      }
      break;

   case SCIP_EVENTTYPE_OBJCHANGED:
      debugMessage("variable <%s> changed objective value from %g to %g\n", 
         SCIPvarGetName(var), SCIPeventGetOldobj(event), SCIPeventGetNewobj(event));
      objdelta = SCIPeventGetNewobj(event) - SCIPeventGetOldobj(event);
      CHECK_OKAY( SCIPaddVarToRow(scip, sepadata->objrow, var, objdelta) );
      break;

   default:
      errorMessage("invalid event type %x\n", SCIPeventGetType(event));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}




/*
 * separator specific interface methods
 */

/** creates the integer objective value separator and includes it in SCIP */
RETCODE SCIPincludeSepaIntobj(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   SEPADATA* sepadata;
   EVENTHDLRDATA* eventhdlrdata;

   /* create intobj separator data */
   CHECK_OKAY( sepadataCreate(scip, &sepadata) );

   /* include separator */
   CHECK_OKAY( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_DELAY,
         sepaFreeIntobj, sepaInitIntobj, sepaExitIntobj, 
         sepaInitsolIntobj, sepaExitsolIntobj, sepaExecIntobj,
         sepadata) );

   /* include event handler for objective change events */
   eventhdlrdata = (EVENTHDLRDATA*)sepadata;
   CHECK_OKAY( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC, 
         eventFreeIntobj, eventInitIntobj, eventExitIntobj, 
         eventInitsolIntobj, eventExitsolIntobj, eventDeleteIntobj, eventExecIntobj,
         eventhdlrdata) );

   /* add intobj separator parameters */
   /* TODO: (optional) add separator specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}

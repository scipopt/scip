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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_logicor.c,v 1.71 2005/01/21 09:16:50 bzfpfend Exp $"

/**@file   cons_logicor.c
 * @brief  constraint handler for logic or constraints
 *         (equivalent to set covering, but algorithms are suited for depth first search)
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_logicor.h"
#include "cons_linear.h"


#define CONSHDLR_NAME          "logicor"
#define CONSHDLR_DESC          "logic or constraints"
#define CONSHDLR_SEPAPRIORITY   +800000
#define CONSHDLR_ENFOPRIORITY   +800000
#define CONSHDLR_CHECKPRIORITY  -800000
#define CONSHDLR_SEPAFREQ             5
#define CONSHDLR_PROPFREQ             1
#define CONSHDLR_EAGERFREQ          100
#define CONSHDLR_MAXPREROUNDS        -1
#define CONSHDLR_NEEDSCONS         TRUE

#define LINCONSUPGD_PRIORITY    +800000

#define EVENTHDLR_NAME         "logicor"
#define EVENTHDLR_DESC         "bound tighten event handler for logic or constraints"

#define CONFLICTHDLR_NAME      "logicor"
#define CONFLICTHDLR_DESC      "conflict handler creating logic or constraints"
#define CONFLICTHDLR_PRIORITY  LINCONSUPGD_PRIORITY


/**@todo make this a parameter setting */
#define AGEINCREASE(n) (1.0 + 0.2*n)


/** constraint handler data */
struct ConshdlrData
{
   EVENTHDLR*       eventhdlr;          /**< event handler for bound tighten events on watched variables */
};

/** logic or constraint data */
struct ConsData
{
   ROW*             row;                /**< LP row, if constraint is already stored in LP row format */
   VAR**            vars;               /**< variables of the constraint */
   int              varssize;           /**< size of vars array */
   int              nvars;              /**< number of variables in the constraint */
   int              watchedvar1;        /**< position of the first watched variable */
   int              watchedvar2;        /**< position of the second watched variable */
   int              filterpos1;         /**< event filter position of first watched variable */
   int              filterpos2;         /**< event filter position of second watched variable */
};




/*
 * Local methods
 */

#if 0
/** installs rounding locks for the given variable in the given logic or constraint */
static
RETCODE lockRounding(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< logic or constraint */
   VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding down may violate the constraint */
   CHECK_OKAY( SCIPlockVarCons(scip, var, cons, TRUE, FALSE) );

   return SCIP_OKAY;
}
#endif

/** removes rounding locks for the given variable in the given logic or constraint */
static
RETCODE unlockRounding(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< logic or constraint */
   VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding down may violate the constraint */
   CHECK_OKAY( SCIPunlockVarCons(scip, var, cons, TRUE, FALSE) );

   return SCIP_OKAY;
}

/** creates constaint handler data for logic or constraint handler */
static
RETCODE conshdlrdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);

   CHECK_OKAY( SCIPallocMemory(scip, conshdlrdata) );

   /* get event handler for catching bound tighten events on watched variables */
   (*conshdlrdata)->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( (*conshdlrdata)->eventhdlr == NULL )
   {
      errorMessage("event handler for logic or constraints not found\n");
      return SCIP_PLUGINNOTFOUND;
   }
   
   return SCIP_OKAY;
}

/** frees constraint handler data for logic or constraint handler */
static
RETCODE conshdlrdataFree(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   SCIPfreeMemory(scip, conshdlrdata);

   return SCIP_OKAY;
}

/** creates a logic or constraint data object */
static
RETCODE consdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to store the logic or constraint data */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars                /**< variables of the constraint */
   )
{
   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);

   CHECK_OKAY( SCIPallocBlockMemory(scip, consdata) );

   (*consdata)->row = NULL;
   if( nvars > 0 )
   {
      CHECK_OKAY( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
      (*consdata)->varssize = nvars;
      (*consdata)->nvars = nvars;
   }
   else
   {
      (*consdata)->vars = NULL;
      (*consdata)->varssize = 0;
      (*consdata)->nvars = 0;
   }
   (*consdata)->watchedvar1 = -1;
   (*consdata)->watchedvar2 = -1;
   (*consdata)->filterpos1 = -1;
   (*consdata)->filterpos2 = -1;

   /* get transformed variables, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      CHECK_OKAY( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );
   }

   return SCIP_OKAY;
}   

/** frees a logic or constraint data */
static
RETCODE consdataFree(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata            /**< pointer to the logic or constraint */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* release the row */
   if( (*consdata)->row != NULL )
   {
      CHECK_OKAY( SCIPreleaseRow(scip, &(*consdata)->row) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** prints logic or constraint to file stream */
static
void consdataPrint(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< logic or constraint data */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   int v;

   assert(consdata != NULL);

   if( file == NULL )
      file = stdout;

   /* print coefficients */
   fprintf(file, "logicor(");
   for( v = 0; v < consdata->nvars; ++v )
   {
      assert(consdata->vars[v] != NULL);
      if( v > 0 )
         fprintf(file, ", ");
      fprintf(file, "<%s>", SCIPvarGetName(consdata->vars[v]));
   }
   fprintf(file, ")\n");
}

/** stores the given variable numbers as watched variables, and updates the event processing */
static
RETCODE switchWatchedvars(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< logic or constraint */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              watchedvar1,        /**< new first watched variable */
   int              watchedvar2         /**< new second watched variable */
   )
{
   CONSDATA* consdata;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(watchedvar1 == -1 || watchedvar1 != watchedvar2);
   assert(watchedvar1 != -1 || watchedvar2 == -1);
   assert(watchedvar1 == -1 || (0 <= watchedvar1 && watchedvar1 < consdata->nvars));
   assert(watchedvar2 == -1 || (0 <= watchedvar2 && watchedvar2 < consdata->nvars));

   /* if one watched variable is equal to the old other watched variable, just switch positions */
   if( watchedvar1 == consdata->watchedvar2 || watchedvar2 == consdata->watchedvar1 )
   {
      int tmp;
      
      tmp = consdata->watchedvar1;
      consdata->watchedvar1 = consdata->watchedvar2;
      consdata->watchedvar2 = tmp;
      tmp = consdata->filterpos1;
      consdata->filterpos1 = consdata->filterpos2;
      consdata->filterpos2 = tmp;
   }
   assert(watchedvar1 == -1 || watchedvar1 != consdata->watchedvar2);
   assert(watchedvar2 == -1 || watchedvar2 != consdata->watchedvar1);

   /* drop events on old watched variables */
   if( consdata->watchedvar1 != -1 && consdata->watchedvar1 != watchedvar1 )
   {
      assert(consdata->filterpos1 != -1);
      CHECK_OKAY( SCIPdropVarEvent(scip, consdata->vars[consdata->watchedvar1],
            SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr, (EVENTDATA*)cons, consdata->filterpos1) );
   }
   if( consdata->watchedvar2 != -1 && consdata->watchedvar2 != watchedvar2 )
   {
      assert(consdata->filterpos2 != -1);
      CHECK_OKAY( SCIPdropVarEvent(scip, consdata->vars[consdata->watchedvar2],
            SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr, (EVENTDATA*)cons, consdata->filterpos2) );
   }

   /* catch events on new watched variables */
   if( watchedvar1 != -1 && watchedvar1 != consdata->watchedvar1 )
   {
      CHECK_OKAY( SCIPcatchVarEvent(scip, consdata->vars[watchedvar1],
            SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr, (EVENTDATA*)cons, &consdata->filterpos1) );
   }
   if( watchedvar2 != -1 && watchedvar2 != consdata->watchedvar2 )
   {
      CHECK_OKAY( SCIPcatchVarEvent(scip, consdata->vars[watchedvar2],
            SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr, (EVENTDATA*)cons, &consdata->filterpos2) );
   }

   /* set the new watched variables */
   consdata->watchedvar1 = watchedvar1;
   consdata->watchedvar2 = watchedvar2;
   
   return SCIP_OKAY;
}

/** deletes coefficient at given position from logic or constraint data */
static
RETCODE delCoefPos(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< logic or constraint */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< position of coefficient to delete */
   )
{
   CONSDATA* consdata;

   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(consdata->vars[pos]));

   /* remove the rounding locks of variable */
   CHECK_OKAY( unlockRounding(scip, cons, consdata->vars[pos]) );

   if( SCIPconsIsTransformed(cons) )
   {
      /* if the position is watched, drop bound tighten events and stop watching the position */
      if( consdata->watchedvar1 == pos )
      {
         CHECK_OKAY( switchWatchedvars(scip, cons, eventhdlr, consdata->watchedvar2, -1) );
      }
      if( consdata->watchedvar2 == pos )
      {
         CHECK_OKAY( switchWatchedvars(scip, cons, eventhdlr, consdata->watchedvar1, -1) );
      }
   }
   assert(pos != consdata->watchedvar1);
   assert(pos != consdata->watchedvar2);

   /* move the last variable to the free slot */
   consdata->vars[pos] = consdata->vars[consdata->nvars-1];
   consdata->nvars--;

   /* if the last variable (that moved) was watched, update the watched position */
   if( consdata->watchedvar1 == consdata->nvars )
      consdata->watchedvar1 = pos;
   if( consdata->watchedvar2 == consdata->nvars )
      consdata->watchedvar2 = pos;

   CHECK_OKAY( SCIPenableConsPropagation(scip, cons) );

   return SCIP_OKAY;
}

/** deletes all zero-fixed variables, checks for variables fixed to one */
static
RETCODE applyFixings(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< logic or constraint */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   Bool*            redundant           /**< returns whether a variable fixed to one exists in the constraint */
   )
{
   CONSDATA* consdata;
   VAR* var;
   int v;

   assert(eventhdlr != NULL);
   assert(redundant != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->vars != NULL);

   *redundant = FALSE;
   v = 0;
   while( v < consdata->nvars )
   {
      var = consdata->vars[v];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);

      if( SCIPvarGetLbGlobal(var) > 0.5 )
      {
         assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(var), 1.0));
         *redundant = TRUE;
         return SCIP_OKAY;
      }
      else if( SCIPvarGetUbGlobal(var) < 0.5 )
      {
         assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 0.0));
         CHECK_OKAY( delCoefPos(scip, cons, eventhdlr, v) );
      }
      else
         ++v;
   }

   debugMessage("after fixings: ");
   debug(consdataPrint(scip, consdata, NULL));

   return SCIP_OKAY;
}

/** analyzes conflicting assignment on given constraint, and adds conflict clause to problem */
static
RETCODE analyzeConflict(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< logic or constraint that detected the conflict */
   )
{
   CONSDATA* consdata;
   int v;

   /* conflict analysis can only be applied in solving stage */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
   CHECK_OKAY( SCIPinitConflictAnalysis(scip) );
   for( v = 0; v < consdata->nvars; ++v )
   {
      CHECK_OKAY( SCIPaddConflictBinvar(scip, consdata->vars[v]) );
   }

   /* analyze the conflict */
   CHECK_OKAY( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** checks constraint for violation only looking at the watched variables, applies fixings if possible */
static
RETCODE processWatchedVars(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< logic or constraint to be processed */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   Bool*            reduceddom,         /**< pointer to store TRUE, if a domain reduction was found */
   Bool*            addcut,             /**< pointer to store whether this constraint must be added as a cut */
   Bool*            mustcheck           /**< pointer to store whether this constraint must be checked for feasibility */
   )
{
   CONSDATA* consdata;
   VAR** vars;
   Longint nbranchings1;
   Longint nbranchings2;
   Longint nbranchings;
   int nvars;
   int watchedvar1;
   int watchedvar2;
   int v;
   Bool infeasible;

   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(reduceddom != NULL);
   assert(addcut != NULL);
   assert(mustcheck != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->watchedvar1 == -1 || consdata->watchedvar1 != consdata->watchedvar2);

   *addcut = FALSE;
   *mustcheck = FALSE;

   debugMessage("processing watched variables of constraint <%s>\n", SCIPconsGetName(cons));

   vars = consdata->vars;
   nvars = consdata->nvars;
   assert(nvars == 0 || vars != NULL);
   watchedvar1 = consdata->watchedvar1;
   watchedvar2 = consdata->watchedvar2;

   /* check watched variables if they are fixed to one */
   if( (watchedvar1 >= 0 && SCIPvarGetLbLocal(vars[watchedvar1]) > 0.5)
      || (watchedvar2 >= 0 && SCIPvarGetLbLocal(vars[watchedvar2]) > 0.5) )
   {
      /* the variable is fixed to one, making the constraint redundant;
       * remember the variable and disable the constraint
       */
      debugMessage(" -> disabling constraint <%s> (watchedvar fixed to 1.0)\n", SCIPconsGetName(cons));
      CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
      return SCIP_OKAY;
   }

   nbranchings1 = LONGINT_MAX;
   nbranchings2 = LONGINT_MAX;
   watchedvar1 = -1;
   watchedvar2 = -1;
   for( v = 0; v < nvars; ++v )
   {
      /* check, if the variable is fixed */
      if( SCIPvarGetUbLocal(vars[v]) > 0.5 )
      {
         if( SCIPvarGetLbLocal(vars[v]) > 0.5 )
         {
            /* the variable is fixed to one, making the constraint redundant;
             * remember the variable and disable the constraint
             */
            debugMessage(" -> disabling constraint <%s> (variable <%s> fixed to 1.0)\n", 
               SCIPconsGetName(cons), SCIPvarGetName(vars[v]));
            CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
            return SCIP_OKAY;
         }

         /* the variable is unfixed and can be used as watched variable */
         nbranchings = SCIPvarGetNBranchingsCurrentRun(vars[v], SCIP_BRANCHDIR_DOWNWARDS);
         if( nbranchings < nbranchings2 )
         {
            if( nbranchings < nbranchings1 )
            {
               watchedvar2 = watchedvar1;
               nbranchings2 = nbranchings1;
               watchedvar1 = v;
               nbranchings1 = nbranchings;
            }
            else
            {
               watchedvar2 = v;
               nbranchings2 = nbranchings;
            }
         }
      }
   }
   assert(nbranchings1 <= nbranchings2);
   assert(watchedvar1 >= 0 || watchedvar2 == -1);

   if( watchedvar1 == -1 )
   {
      /* there is no unfixed variable left -> the constraint is infeasible
       *  - a modifiable constraint must be added as a cut and further pricing must be performed in the LP solving loop
       *  - an unmodifiable constraint is infeasible and the node can be cut off
       */
      assert(watchedvar2 == -1);

      debugMessage(" -> constraint <%s> is infeasible\n", SCIPconsGetName(cons));

      CHECK_OKAY( SCIPresetConsAge(scip, cons) );
      if( SCIPconsIsModifiable(cons) )
         *addcut = TRUE;
      else
      {
         /* use conflict analysis to get a conflict clause out of the conflicting assignment */
         CHECK_OKAY( analyzeConflict(scip, cons) );

         /* mark the node to be cut off */
         *cutoff = TRUE;
      }
   }
   else if( watchedvar2 == -1 )
   {
      /* there is only one unfixed variable:
       * - a modifiable constraint must be checked manually
       * - an unmodifiable constraint is feasible and can be disabled after the remaining variable is fixed to one
       */
      assert(0 <= watchedvar1 && watchedvar1 < nvars);
      assert(SCIPisEQ(scip, SCIPvarGetLbLocal(vars[watchedvar1]), 0.0));
      assert(SCIPisEQ(scip, SCIPvarGetUbLocal(vars[watchedvar1]), 1.0));
      if( SCIPconsIsModifiable(cons) )
         *mustcheck = TRUE;
      else
      {
         /* fixed remaining variable to one and disable constraint */
         debugMessage(" -> single-literal constraint <%s> (fix <%s> to 1.0) at depth %d\n", 
            SCIPconsGetName(cons), SCIPvarGetName(vars[watchedvar1]), SCIPgetDepth(scip));
         CHECK_OKAY( SCIPinferBinvarCons(scip, vars[watchedvar1], TRUE, cons, 0, &infeasible, NULL) );
         assert(!infeasible);
         CHECK_OKAY( SCIPresetConsAge(scip, cons) );
         CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
         *reduceddom = TRUE;
      }
   }
   else
   {
      debugMessage(" -> new watched variables <%s> and <%s> of constraint <%s> are still unfixed\n",
         SCIPvarGetName(vars[watchedvar1]), SCIPvarGetName(vars[watchedvar2]), SCIPconsGetName(cons));

      /* switch to the new watched variables */
      CHECK_OKAY( switchWatchedvars(scip, cons, eventhdlr, watchedvar1, watchedvar2) );

      /* there are at least two unfixed variables -> the constraint must be checked manually */
      *mustcheck = TRUE;

      /* disable propagation of constraint until a watched variable gets fixed */
      CHECK_OKAY( SCIPdisableConsPropagation(scip, cons) );

      /* increase aging counter */
      CHECK_OKAY( SCIPaddConsAge(scip, cons, AGEINCREASE(consdata->nvars)) );
   }

   return SCIP_OKAY;
}

/** checks constraint for violation, returns TRUE iff constraint is feasible */
static
RETCODE checkCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< logic or constraint to be checked */
   SOL*             sol,                /**< primal CIP solution */
   Bool*            violated            /**< pointer to store whether the given solution violates the constraint */
   )
{
   CONSDATA* consdata;
   VAR** vars;
   Real solval;
   Real sum;
   int nvars;
   int v;

   assert(violated != NULL);

   *violated = FALSE;
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   vars = consdata->vars;
   nvars = consdata->nvars;
   
   /* if we should check the current LP or pseudo solution, look for a fixed-to-one variable in order to disable
    * the constraint
    */
   if( sol == NULL )
   {
      for( v = 0; v < nvars; ++v )
      {
         if( SCIPvarGetLbLocal(vars[v]) > 0.5 )
         {
            CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
            return SCIP_OKAY;
         }
      }
   }

   /* calculate the constraint's activity */
   sum = 0.0;
   solval = 0.0;
   for( v = 0; v < nvars && sum < 1.0; ++v )
   {
      assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY);
      solval = SCIPgetSolVal(scip, sol, vars[v]);
      assert(SCIPisFeasGE(scip, solval, 0.0) && SCIPisFeasLE(scip, solval, 1.0));
      sum += solval;
   }

   *violated = SCIPisFeasLT(scip, sum, 1.0);

   return SCIP_OKAY;
}

/** creates an LP row in a logic or constraint data object */
static
RETCODE createRow(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< logic or constraint */
   )
{
   CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);

   CHECK_OKAY( SCIPcreateEmptyRow(scip, &consdata->row, SCIPconsGetName(cons), 1.0, SCIPinfinity(scip),
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
   
   CHECK_OKAY( SCIPaddVarsToRowSameCoef(scip, consdata->row, consdata->nvars, consdata->vars, 1.0) );

   return SCIP_OKAY;
}

/** adds logic or constraint as cut to the LP */
static
RETCODE addCut(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< logic or constraint */
   )
{
   CONSDATA* consdata;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   if( consdata->row == NULL )
   {
      /* convert logic or constraint data into LP row */
      CHECK_OKAY( createRow(scip, cons) );
   }
   assert(consdata->row != NULL);
   assert(!SCIProwIsInLP(consdata->row));

   debugMessage("adding constraint <%s> as cut to the LP\n", SCIPconsGetName(cons));
            
   /* insert LP row as cut */
   CHECK_OKAY( SCIPaddCut(scip, consdata->row, FALSE) );

   return SCIP_OKAY;
}

/** checks constraint for violation, and adds it as a cut if possible */
static
RETCODE separateCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< logic or constraint to be separated */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   Bool*            separated,          /**< pointer to store TRUE, if a cut was found */
   Bool*            reduceddom          /**< pointer to store TRUE, if a domain reduction was found */
   )
{
   Bool addcut;
   Bool mustcheck;

   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(separated != NULL);
   assert(reduceddom != NULL);

   debugMessage("separating constraint <%s>\n", SCIPconsGetName(cons));

   /* update and check the watched variables, if they were changed since last processing */
   if( SCIPconsIsPropagationEnabled(cons) )
   {
      CHECK_OKAY( processWatchedVars(scip, cons, eventhdlr, cutoff, reduceddom, &addcut, &mustcheck) );
   }
   else
   {
      addcut = FALSE;
      mustcheck = TRUE;
   }

   if( mustcheck )
   {
      CONSDATA* consdata;

      assert(!addcut);
      
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* variable's fixings didn't give us any information -> we have to check the constraint */
      if( consdata->row != NULL )
      {
         /* skip constraints already in the LP */
         if( SCIProwIsInLP(consdata->row) )
            return SCIP_OKAY;
         else
         {
            Real feasibility;
            
            assert(!SCIProwIsInLP(consdata->row));
            feasibility = SCIPgetRowLPFeasibility(scip, consdata->row);
            addcut = SCIPisFeasNegative(scip, feasibility);
         }
      }
      else
      {
         CHECK_OKAY( checkCons(scip, cons, NULL, &addcut) );
      }
   }

   if( addcut )
   {
      /* insert LP row as cut */
      CHECK_OKAY( addCut(scip, cons) );
      CHECK_OKAY( SCIPresetConsAge(scip, cons) );
      *separated = TRUE;
   }

   return SCIP_OKAY;
}

/** enforces the pseudo solution on the given constraint */
static
RETCODE enforcePseudo(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< logic or constraint to be separated */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   Bool*            infeasible,         /**< pointer to store TRUE, if the constraint was infeasible */
   Bool*            reduceddom,         /**< pointer to store TRUE, if a domain reduction was found */
   Bool*            solvelp             /**< pointer to store TRUE, if the LP has to be solved */
   )
{
   Bool addcut;
   Bool mustcheck;

   assert(!SCIPhasCurrentNodeLP(scip));
   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(infeasible != NULL);
   assert(reduceddom != NULL);
   assert(solvelp != NULL);

   /* update and check the watched variables, if they were changed since last processing */
   if( SCIPconsIsPropagationEnabled(cons) )
   {
      CHECK_OKAY( processWatchedVars(scip, cons, eventhdlr, cutoff, reduceddom, &addcut, &mustcheck) );
   }
   else
   {
      addcut = FALSE;
      mustcheck = TRUE;
   }

   if( mustcheck )
   {
      Bool violated;

      assert(!addcut);

      CHECK_OKAY( checkCons(scip, cons, NULL, &violated) );
      if( violated )
      {
         /* constraint was infeasible -> reset age */
         CHECK_OKAY( SCIPresetConsAge(scip, cons) );
         *infeasible = TRUE;
      }
   }
   else if( addcut )
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

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
DECL_CONSFREE(consFreeLogicor)
{  /*lint --e{715}*/
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


/** initialization method of constraint handler (called after problem was transformed) */
#define consInitLogicor NULL


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitLogicor NULL


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#define consInitpreLogicor NULL


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreLogicor NULL


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolLogicor NULL


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
DECL_CONSEXITSOL(consExitsolLogicor)
{
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
DECL_CONSDELETE(consDeleteLogicor)
{  /*lint --e{715}*/
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* free LP row and logic or constraint */
   CHECK_OKAY( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
DECL_CONSTRANS(consTransLogicor)
{  /*lint --e{715}*/
   CONSDATA* sourcedata;
   CONSDATA* targetdata;

   /*debugMessage("Trans method of logic or constraints\n");*/

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->row == NULL);  /* in original problem, there cannot be LP rows */

   /* create constraint data for target constraint */
   CHECK_OKAY( consdataCreate(scip, &targetdata, sourcedata->nvars, sourcedata->vars) );

   /* create target constraint */
   CHECK_OKAY( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), SCIPconsIsRemoveable(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
DECL_CONSINITLP(consInitlpLogicor)
{  /*lint --e{715}*/
   int c;

   for( c = 0; c < nconss; ++c )
   {
      if( SCIPconsIsInitial(conss[c]) )
      {
         CHECK_OKAY( addCut(scip, conss[c]) );
      }
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler */
static
DECL_CONSSEPA(consSepaLogicor)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   Bool cutoff;
   Bool separated;
   Bool reduceddom;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   debugMessage("separating %d/%d logic or constraints\n", nusefulconss, nconss);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   separated = FALSE;
   reduceddom = FALSE;

   /* check all useful logic or constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff && !reduceddom; ++c )
   {
      CHECK_OKAY( separateCons(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &separated, &reduceddom) );
   }

   /* combine logic or constraints to get more cuts */
   /**@todo further cuts of logic or constraints */

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( separated )
      *result = SCIP_SEPARATED;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
DECL_CONSENFOLP(consEnfolpLogicor)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   Bool cutoff;
   Bool separated;
   Bool reduceddom;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   debugMessage("LP enforcing %d logic or constraints\n", nconss);

   *result = SCIP_FEASIBLE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   separated = FALSE;
   reduceddom = FALSE;

   /* check all useful logic or constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff && !reduceddom; ++c )
   {
      CHECK_OKAY( separateCons(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &separated, &reduceddom) );
   }

   /* check all obsolete logic or constraints for feasibility */
   for( c = nusefulconss; c < nconss && !cutoff && !separated && !reduceddom; ++c )
   {
      CHECK_OKAY( separateCons(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &separated, &reduceddom) );
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


/** constraint enforcing method of constraint handler for pseudo solutions */
static
DECL_CONSENFOPS(consEnfopsLogicor)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   Bool cutoff;
   Bool infeasible;
   Bool reduceddom;
   Bool solvelp;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   debugMessage("pseudo enforcing %d logic or constraints\n", nconss);

   *result = SCIP_FEASIBLE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   infeasible = FALSE;
   reduceddom = FALSE;
   solvelp = FALSE;

   /* check all logic or constraints for feasibility */
   for( c = 0; c < nconss && !cutoff && !reduceddom && !solvelp && !infeasible; ++c )
   {
      CHECK_OKAY( enforcePseudo(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &infeasible, &reduceddom, &solvelp) );
   }

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;
   else if( solvelp )
      *result = SCIP_SOLVELP;
   else if( infeasible )
      *result = SCIP_INFEASIBLE;
   
   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
DECL_CONSCHECK(consCheckLogicor)
{  /*lint --e{715}*/
   CONS* cons;
   CONSDATA* consdata;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   /* check all logic or constraints for feasibility */
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);
      if( checklprows || consdata->row == NULL || !SCIProwIsInLP(consdata->row) )
      {
         Bool violated;

         CHECK_OKAY( checkCons(scip, cons, sol, &violated) );
         if( violated )
         {
            /* constraint is violated */
            *result = SCIP_INFEASIBLE;
            return SCIP_OKAY;
         }
      }
   }
   
   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
DECL_CONSPROP(consPropLogicor)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   Bool cutoff;
   Bool reduceddom;
   Bool addcut;
   Bool mustcheck;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   reduceddom = FALSE;

   /* propagate all useful logic or constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      CHECK_OKAY( processWatchedVars(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &reduceddom, &addcut, &mustcheck) );
   }

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
DECL_CONSPRESOL(consPresolLogicor)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   CONS* cons;
   CONSDATA* consdata;
   Bool infeasible;
   Bool redundant;
   Bool fixed;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* process constraints */
   for( c = 0; c < nconss && *result != SCIP_CUTOFF; ++c )
   {
      cons = conss[c];
      assert(cons != NULL);
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      debugMessage("presolving logic or constraint <%s>\n", SCIPconsGetName(cons));

      /* force presolving the constraint in the initial round */
      if( nrounds == 0 )
      {
         CHECK_OKAY( SCIPenableConsPropagation(scip, cons) );
      }

      /* remove all variables that are fixed to zero, check redundancy due to fixed-to-one variable */
      CHECK_OKAY( applyFixings(scip, cons, conshdlrdata->eventhdlr, &redundant) );

      /**@todo find pairs of negated variables in constraint: constraint is redundant */
      /**@todo find sets of equal variables in constraint: multiple entries of variable can be replaced by single entry */

      if( redundant )
      {
         debugMessage("logic or constraint <%s> is redundant\n", SCIPconsGetName(cons));
         CHECK_OKAY( SCIPdelCons(scip, cons) );
         (*ndelconss)++;
         *result = SCIP_SUCCESS;
         continue;
      }
      else if( !SCIPconsIsModifiable(cons) )
      {
         /* if unmodifiable constraint has no variables, it is infeasible,
          * if unmodifiable constraint has only one variable, this one can be fixed and the constraint deleted
          */
         if( consdata->nvars == 0 )
         {
            debugMessage("logic or constraint <%s> is infeasible\n", SCIPconsGetName(cons));
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else if( consdata->nvars == 1 )
         {
            debugMessage("logic or constraint <%s> has only one variable not fixed to 0.0\n",
               SCIPconsGetName(cons));
            
            assert(consdata->vars != NULL);
            assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[0]), 0.0));
            assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[0]), 1.0));
            
            CHECK_OKAY( SCIPfixVar(scip, consdata->vars[0], 1.0, &infeasible, &fixed) );
            if( infeasible )
            {
               debugMessage(" -> infeasible fixing\n");
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            assert(fixed);
            (*nfixedvars)++;

            CHECK_OKAY( SCIPdelCons(scip, cons) );
            (*ndelconss)++;
            *result = SCIP_SUCCESS;
            continue;
         }
      }
   }
   
   /**@todo preprocess pairs of logic or constraints */

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
DECL_CONSRESPROP(consRespropLogicor)
{  /*lint --e{715}*/
   CONSDATA* consdata;
   Bool infervarfound;
   int v;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);
   assert(infervar != NULL);
   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   debugMessage("conflict resolving method of logic or constraint handler\n");

   /* the only deductions are variables infered to 1.0 on logic or constraints where all other variables
    * are assigned to zero
    */
   assert(SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE) > 0.5); /* the inference variable must be assigned to one */

   infervarfound = FALSE;
   for( v = 0; v < consdata->nvars; ++v )
   {
      if( consdata->vars[v] != infervar )
      {
         /* the reason variable must have been assigned to zero */
         assert(SCIPvarGetUbAtIndex(consdata->vars[v], bdchgidx, FALSE) < 0.5);
         CHECK_OKAY( SCIPaddConflictBinvar(scip, consdata->vars[v]) );
      }
      else
      {
         assert(!infervarfound);
         infervarfound = TRUE;
      }
   }
   assert(infervarfound);

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
DECL_CONSLOCK(consLockLogicor)
{  /*lint --e{715}*/
   CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* lock every single coefficient */
   for( i = 0; i < consdata->nvars; ++i )
   {
      CHECK_OKAY( SCIPaddVarLocks(scip, consdata->vars[i], nlockspos, nlocksneg) );
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
DECL_CONSACTIVE(consActiveLogicor)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* consdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   debugMessage("activating information for logic or constraint <%s>\n", SCIPconsGetName(cons));
   debug(consdataPrint(scip, consdata, NULL));

   /* catch events on watched variables */
   if( consdata->watchedvar1 != -1 )
   {
      CHECK_OKAY( SCIPcatchVarEvent(scip, consdata->vars[consdata->watchedvar1],
            SCIP_EVENTTYPE_BOUNDTIGHTENED, conshdlrdata->eventhdlr, (EVENTDATA*)cons, &consdata->filterpos1) );
   }
   if( consdata->watchedvar2 != -1 )
   {
      CHECK_OKAY( SCIPcatchVarEvent(scip, consdata->vars[consdata->watchedvar2],
            SCIP_EVENTTYPE_BOUNDTIGHTENED, conshdlrdata->eventhdlr, (EVENTDATA*)cons, &consdata->filterpos2) );
   }

   return SCIP_OKAY;
}


/** constraint deactivation notification method of constraint handler */
static
DECL_CONSDEACTIVE(consDeactiveLogicor)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* consdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   debugMessage("deactivating information for logic or constraint <%s>\n", SCIPconsGetName(cons));
   debug(consdataPrint(scip, consdata, NULL));

   /* drop events on watched variables */
   if( consdata->watchedvar1 != -1 )
   {
      assert(consdata->filterpos1 != -1);
      CHECK_OKAY( SCIPdropVarEvent(scip, consdata->vars[consdata->watchedvar1],
            SCIP_EVENTTYPE_BOUNDTIGHTENED, conshdlrdata->eventhdlr, (EVENTDATA*)cons, consdata->filterpos1) );
   }
   if( consdata->watchedvar2 != -1 )
   {
      assert(consdata->filterpos2 != -1);
      CHECK_OKAY( SCIPdropVarEvent(scip, consdata->vars[consdata->watchedvar2],
            SCIP_EVENTTYPE_BOUNDTIGHTENED, conshdlrdata->eventhdlr, (EVENTDATA*)cons, consdata->filterpos2) );
   }

   return SCIP_OKAY;
}


/** constraint enabling notification method of constraint handler */
#define consEnableLogicor NULL


/** constraint disabling notification method of constraint handler */
#define consDisableLogicor NULL


/** constraint display method of constraint handler */
static
DECL_CONSPRINT(consPrintLogicor)
{
   consdataPrint(scip, SCIPconsGetData(cons), file);

   return SCIP_OKAY;
}




/*
 * upgrading of linear constraints
 */

/** creates and captures a normalized (with all coefficients +1) logic or constraint */
static
RETCODE createNormalizedLogicor(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Real*            vals,               /**< array with coefficients (+1.0 or -1.0) */
   int              mult,               /**< multiplier on the coefficients(+1 or -1) */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   VAR** transvars;
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);
   assert(mult == +1 || mult == -1);

   /* get temporary memory */
   CHECK_OKAY( SCIPallocBufferArray(scip, &transvars, nvars) );

   /* negate positive or negative variables */
   for( v = 0; v < nvars; ++v )
   {
      if( mult * vals[v] > 0.0 )
         transvars[v] = vars[v];
      else
      {
         CHECK_OKAY( SCIPgetNegatedVar(scip, vars[v], &transvars[v]) );
      }
      assert(transvars[v] != NULL);
   }

   /* create the constraint */
   CHECK_OKAY( SCIPcreateConsLogicor(scip, cons, name, nvars, transvars,
         initial, separate, enforce, check, propagate, local, modifiable, removeable) );

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &transvars);

   return SCIP_OKAY;
}

static
DECL_LINCONSUPGD(linconsUpgdLogicor)
{  /*lint --e{715}*/
   assert(upgdcons != NULL);

   /* check, if linear constraint can be upgraded to logic or constraint
    * - logic or constraints consist only of binary variables with a
    *   coefficient of +1.0 or -1.0 (variables with -1.0 coefficients can be negated):
    *        lhs     <= x1 + ... + xp - y1 - ... - yn <= rhs
    * - negating all variables y = (1-Y) with negative coefficients gives:
    *        lhs + n <= x1 + ... + xp + Y1 + ... + Yn <= rhs + n
    * - negating all variables x = (1-X) with positive coefficients and multiplying with -1 gives:
    *        p - rhs <= X1 + ... + Xp + y1 + ... + yn <= p - lhs
    * - logic or constraints have left hand side of +1.0, and right hand side of +infinity: x(S) >= 1.0
    *    -> without negations:  (lhs == 1 - n  and  rhs == +inf)  or  (lhs == -inf  and  rhs = p - 1)
    */
   if( nposbin + nnegbin == nvars && ncoeffspone + ncoeffsnone == nvars
      && ((SCIPisEQ(scip, lhs, 1.0 - ncoeffsnone) && SCIPisInfinity(scip, rhs))
         || (SCIPisInfinity(scip, -lhs) && SCIPisEQ(scip, rhs, ncoeffspone - 1.0))) )
   {
      int mult;

      debugMessage("upgrading constraint <%s> to logic or constraint\n", SCIPconsGetName(cons));
      
      /* check, if we have to multiply with -1 (negate the positive vars) or with +1 (negate the negative vars) */
      mult = SCIPisInfinity(scip, rhs) ? +1 : -1;
      
      /* create the logic or constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      CHECK_OKAY( createNormalizedLogicor(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, mult,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
   }

   return SCIP_OKAY;
}




/*
 * Callback methods of event handler
 */

static
DECL_EVENTEXEC(eventExecLogicor)
{  /*lint --e{715}*/
   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   debugMessage("Exec method of bound tighten event handler for logic or constraints\n");

   CHECK_OKAY( SCIPenableConsPropagation(scip, (CONS*)eventdata) );

   return SCIP_OKAY;
}




/*
 * Callback methods of conflict handler
 */

static
DECL_CONFLICTEXEC(conflictExecLogicor)
{  /*lint --e{715}*/
   CONS* cons;
   char consname[MAXSTRLEN];
   
   assert(conflicthdlr != NULL);
   assert(strcmp(SCIPconflicthdlrGetName(conflicthdlr), CONFLICTHDLR_NAME) == 0);
   assert(conflictvars != NULL || nconflictvars == 0);
   assert(result != NULL);

   /* don't process already resolved conflicts */
   if( resolved )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* create a constraint out of the conflict set */
   sprintf(consname, "cf%d_%lld", SCIPgetNRuns(scip), SCIPgetNConflictClausesFound(scip));
   CHECK_OKAY( SCIPcreateConsLogicor(scip, &cons, consname, nconflictvars, conflictvars, 
         FALSE, TRUE, FALSE, FALSE, TRUE, local, FALSE, TRUE) );
   CHECK_OKAY( SCIPaddConsNode(scip, node, cons) );
   CHECK_OKAY( SCIPreleaseCons(scip, &cons) );

   *result = SCIP_CONSADDED;

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for logic or constraints and includes it in SCIP */
RETCODE SCIPincludeConshdlrLogicor(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create event handler for bound tighten events on watched variables */
   CHECK_OKAY( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, NULL, NULL,
         NULL, eventExecLogicor,
         NULL) );

   /* create conflict handler for logic or constraints */
   CHECK_OKAY( SCIPincludeConflicthdlr(scip, CONFLICTHDLR_NAME, CONFLICTHDLR_DESC, CONFLICTHDLR_PRIORITY,
         NULL, NULL, NULL, conflictExecLogicor,
         NULL) );

   /* create constraint handler data */
   CHECK_OKAY( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, CONSHDLR_NEEDSCONS,
         consFreeLogicor, consInitLogicor, consExitLogicor, 
         consInitpreLogicor, consExitpreLogicor, consInitsolLogicor, consExitsolLogicor,
         consDeleteLogicor, consTransLogicor, 
         consInitlpLogicor, consSepaLogicor, 
         consEnfolpLogicor, consEnfopsLogicor, consCheckLogicor, 
         consPropLogicor, consPresolLogicor, consRespropLogicor, consLockLogicor,
         consActiveLogicor, consDeactiveLogicor,
         consEnableLogicor, consDisableLogicor,
         consPrintLogicor,
         conshdlrdata) );

   /* include the linear constraint to logicor constraint upgrade in the linear constraint handler */
   CHECK_OKAY( SCIPincludeLinconsUpgrade(scip, linconsUpgdLogicor, LINCONSUPGD_PRIORITY) );

   return SCIP_OKAY;
}


/** creates and captures a logic or constraint */
RETCODE SCIPcreateConsLogicor(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   assert(scip != NULL);

   /* find the logicor constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("logic or constraint handler not found\n");
      return SCIP_INVALIDCALL;
   }

   /* create the constraint specific data */
   CHECK_OKAY( consdataCreate(scip, &consdata, nvars, vars) );

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, removeable) );

   return SCIP_OKAY;
}

/** gets the dual solution of the logic or constraint in the current LP */
Real SCIPgetDualsolLogicor(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   )
{
   CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not a logic or constraint\n");
      abort();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualsol(consdata->row);
   else
      return 0.0;
}

/** gets array of variables in logic or constraint */
VAR** SCIPgetVarsLogicor(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   )
{
   CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not a logic or constraint\n");
      abort();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vars;
}

/** gets number of variables in logic or constraint */
int SCIPgetNVarsLogicor(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   )
{
   CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not a logic or constraint\n");
      abort();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nvars;
}


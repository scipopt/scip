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

/**@file   cons_setcover.c
 * @brief  constraint handler for setcover constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_setcover.h"
#include "cons_linear.h"


#define CONSHDLR_NAME          "setcover"
#define CONSHDLR_DESC          "set covering constraint"
#define CONSHDLR_SEPAPRIORITY   +800000
#define CONSHDLR_ENFOPRIORITY   +800000
#define CONSHDLR_CHECKPRIORITY  -800000
#define CONSHDLR_SEPAFREQ             4
#define CONSHDLR_PROPFREQ            -1
#define CONSHDLR_NEEDSCONS         TRUE

#define EVENTHDLR_NAME         "setcovering"
#define EVENTHDLR_DESC         "bound change event handler for set covering constraints"

#define LINCONSUPGD_PRIORITY   +1000000

#define MINBRANCHWEIGHT               0.4  /**< minimum weight of both sets in set covering branching */
#define MAXBRANCHWEIGHT               0.8  /**< minimum weight of both sets in set covering branching */
#define MAXBRANCHACTIVITY             1.1  /**< maximum activity of constraint to check if branching is possible */


/** constraint handler data */
struct ConsHdlrData
{
   INTARRAY*        varuses;            /**< number of times a variable is used in the active set covering constraints */
};

/** set covering constraint data */
struct SetcoverCons
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
};
typedef struct SetcoverCons SETCOVERCONS; /**< set covering constraint data */

/** constraint data for set covering constraints */
struct ConsData
{
   SETCOVERCONS*    setcovercons;       /**< set covering constraint data */
   ROW*             row;                /**< LP row, if constraint is already stored in LP row format */
};




/*
 * Local methods
 */

/** creates constaint handler data for linear constraint handler */
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

/** frees constraint handler data for linear constraint handler */
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

/** creates a set covering constraint data object */
static
RETCODE setcoverconsCreate(
   SCIP*            scip,               /**< SCIP data structure */
   SETCOVERCONS**   setcovercons,       /**< pointer to store the set covering constraint */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< variables of the constraint */
   Bool             local,              /**< is constraint only locally valid? */
   Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   assert(setcovercons != NULL);
   assert(nvars == 0 || vars != NULL);

   CHECK_OKAY( SCIPallocBlockMemory(scip, setcovercons) );
   if( nvars > 0 )
   {
      CHECK_OKAY( SCIPduplicateBlockMemoryArray(scip, &(*setcovercons)->vars, vars, nvars) );
      (*setcovercons)->varssize = nvars;
      (*setcovercons)->nvars = nvars;
   }
   else
   {
      (*setcovercons)->vars = NULL;
      (*setcovercons)->varssize = 0;
      (*setcovercons)->nvars = 0;
   }
   (*setcovercons)->nfixedzeros = 0;
   (*setcovercons)->nfixedones = 0;
   (*setcovercons)->local = local;
   (*setcovercons)->modifiable = modifiable;
   (*setcovercons)->removeable = removeable;
   (*setcovercons)->transformed = FALSE;

   return SCIP_OKAY;
}   

/** creates a transformed set covering constraint data object */
static
RETCODE setcoverconsCreateTransformed(
   SCIP*            scip,               /**< SCIP data structure */
   SETCOVERCONS**   setcovercons,       /**< pointer to store the set covering constraint */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< variables of the constraint */
   Bool             local,              /**< is constraint only locally valid? */
   Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   EVENTHDLR* eventhdlr;
   int v;

   assert(setcovercons != NULL);
   assert(nvars == 0 || vars != NULL);

   CHECK_OKAY( setcoverconsCreate(scip, setcovercons, nvars, vars, local, modifiable, removeable) );
   (*setcovercons)->transformed = TRUE;

   /* find the bound change event handler */
   eventhdlr = SCIPfindEventHdlr(scip, EVENTHDLR_NAME);
   if( eventhdlr == NULL )
   {
      errorMessage("event handler for processing bound change events on set covering constraints not found");
      return SCIP_ERROR;
   }

   for( v = 0; v < (*setcovercons)->nvars; ++v )
   {
      /* use transformed variables in constraint instead original ones */
      if( SCIPvarGetStatus((*setcovercons)->vars[v]) == SCIP_VARSTATUS_ORIGINAL )
      {
         (*setcovercons)->vars[v] = SCIPvarGetTransformed((*setcovercons)->vars[v]);
         assert((*setcovercons)->vars[v] != NULL);
      }
      assert(SCIPvarGetStatus((*setcovercons)->vars[v]) != SCIP_VARSTATUS_ORIGINAL);
      
      /* catch bound change events on variables */
      CHECK_OKAY( SCIPcatchVarEvent(scip, (*setcovercons)->vars[v], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr,
                     (EVENTDATA*)(*setcovercons)) );
   }

   return SCIP_OKAY;
}

/** frees a set covering constraint data */
static
RETCODE setcoverconsFree(
   SCIP*            scip,               /**< SCIP data structure */
   SETCOVERCONS**   setcovercons        /**< pointer to store the set covering constraint */
   )
{
   assert(setcovercons != NULL);
   assert(*setcovercons != NULL);

   if( (*setcovercons)->transformed )
   {
      EVENTHDLR* eventhdlr;
      int v;
      
      /* find the bound change event handler */
      eventhdlr = SCIPfindEventHdlr(scip, EVENTHDLR_NAME);
      if( eventhdlr == NULL )
      {
         errorMessage("event handler for processing bound change events on set covering constraints not found");
         return SCIP_ERROR;
      }
      
      /* drop bound change events on variables */
      for( v = 0; v < (*setcovercons)->nvars; ++v )
      {
         CHECK_OKAY( SCIPdropVarEvent(scip, (*setcovercons)->vars[v], eventhdlr, (EVENTDATA*)(*setcovercons)) );
      }
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*setcovercons)->vars, (*setcovercons)->varssize);
   SCIPfreeBlockMemory(scip, setcovercons);

   return SCIP_OKAY;
}

/** creates an LP row from a set covering constraint data object */
static
RETCODE setcoverconsToRow(
   SCIP*            scip,               /**< SCIP data structure */
   SETCOVERCONS*    setcovercons,       /**< set covering constraint data */
   const char*      name,               /**< name of the constraint */
   ROW**            row                 /**< pointer to an LP row data object */
   )
{
   int v;

   assert(setcovercons != NULL);
   assert(row != NULL);

   CHECK_OKAY( SCIPcreateRow(scip, row, name, 0, NULL, NULL, 1.0, SCIPinfinity(scip),
                  setcovercons->local, setcovercons->modifiable, setcovercons->removeable) );
   
   for( v = 0; v < setcovercons->nvars; ++v )
   {
      CHECK_OKAY( SCIPaddVarToRow(scip, *row, setcovercons->vars[v], 1.0) );
   }

   return SCIP_OKAY;
}




/*
 * Callback methods of constraint handler
 */

static
DECL_CONSFREE(consFreeSetcover)
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
DECL_CONSDELETE(consDeleteSetcover)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* free LP row and setcover constraint */
   if( (*consdata)->row != NULL )
   {
      CHECK_OKAY( SCIPreleaseRow(scip, &(*consdata)->row) );
   }
   CHECK_OKAY( setcoverconsFree(scip, &(*consdata)->setcovercons) );

   /* free constraint data object */
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

static
DECL_CONSTRANS(consTransSetcover)
{
   CONSDATA* sourcedata;
   CONSDATA* targetdata;
   SETCOVERCONS* setcovercons;

   /*debugMessage("Trans method of setcover constraints\n");*/

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPstage(scip) == SCIP_STAGE_INITSOLVE);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->row == NULL);  /* in original problem, there cannot be LP rows */
   assert(sourcedata->setcovercons != NULL);

   /* create constraint data for target constraint */
   CHECK_OKAY( SCIPallocBlockMemory(scip, &targetdata) );

   setcovercons = sourcedata->setcovercons;

   CHECK_OKAY( setcoverconsCreateTransformed(scip, &targetdata->setcovercons, setcovercons->nvars, setcovercons->vars,
                  setcovercons->local, setcovercons->modifiable, setcovercons->removeable) );
   targetdata->row = NULL;

   /* create target constraint */
   assert(!SCIPconsIsPropagated(sourcecons));
   CHECK_OKAY( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
                  SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
                  SCIPconsIsPropagated(sourcecons)) );

   return SCIP_OKAY;
}

/** checks constraint for violation, and adds it as a cut if possible */
static
RETCODE separate(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< set covering constraint to be separated */
   RESULT*          result              /**< pointer to store the result, if one was found */
   )
{
   CONSDATA* consdata;
   SETCOVERCONS* setcovercons;

   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   setcovercons = consdata->setcovercons;
   assert(setcovercons != NULL);
   assert(0 <= setcovercons->nfixedzeros && setcovercons->nfixedzeros <= setcovercons->nvars);
   assert(0 <= setcovercons->nfixedones && setcovercons->nfixedones <= setcovercons->nvars);

   if( setcovercons->nfixedzeros == setcovercons->nvars )
   {
      if( setcovercons->modifiable )
      {
         /* the constraint cannot be feasible with the active variable set, but due to additional pricing,
          * it may be feasible after the next pricing loop -> just insert it as a cut
          */

         /* convert set covering constraint data into LP row */
         CHECK_OKAY( setcoverconsToRow(scip, setcovercons, SCIPconsGetName(cons), &consdata->row) );
         assert(consdata->row != NULL);
         assert(!SCIProwIsInLP(consdata->row));
         
         /* insert LP row as cut */
         CHECK_OKAY( SCIPaddCut(scip, consdata->row, 1.0/(setcovercons->nvars+1)) );
         *result = SCIP_SEPARATED;
      }
      else
      {
         /* the constraint cannot be feasible, because all variables are fixed to zero */
         *result = SCIP_CUTOFF;
      }
      CHECK_OKAY( SCIPresetConsAge(scip, cons) );
   }
   else if( setcovercons->nfixedones >= 1 )
   {
      /* constraint is feasible anyway, because one of its variables if fixed to one -> constraint can be disabled */
      CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
   }
   else
   {
      ROW* row;

      /* we have to check the constraint */
      row = consdata->row;
      if( row != NULL )
      {
         /* ignore constraints that are in the LP */
         if( !SCIProwIsInLP(row) )
         {         
            Real feasibility;
            
            feasibility = SCIPgetRowLPFeasibility(scip, row);
            /*debugMessage("  row feasibility = %g\n", feasibility);*/
            if( !SCIPisFeasible(scip, feasibility) )
            {
               /* insert LP row as cut */
               CHECK_OKAY( SCIPaddCut(scip, row, -feasibility/SCIProwGetNorm(row)/(SCIProwGetNNonz(row)+1)) );
               CHECK_OKAY( SCIPresetConsAge(scip, cons) );
               *result = SCIP_SEPARATED;
            }
            else
            {
               CHECK_OKAY( SCIPincConsAge(scip, cons) );
            }
         }
      }
      else
      {
         VAR** vars;
         Real solval;
         Real sum;
         int nvars;
         int v;

         /* calculate the constraint's activity */
         vars = setcovercons->vars;
         nvars = setcovercons->nvars;
         sum = 0.0;
         for( v = 0; v < nvars && sum < 1.0; ++v )
         {
            assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY);
            solval = SCIPgetSolVal(scip, NULL, vars[v]);
            assert(SCIPisGE(scip, solval, 0.0) || SCIPisLE(scip, solval, 1.0));
            sum += solval;
         }
         
         if( SCIPisFeasLT(scip, sum, 1.0) )
         {
            /* because sum is smaller than 1.0, the constraint itself can be used as cut */
            
            /* convert set covering constraint data into LP row */
            CHECK_OKAY( setcoverconsToRow(scip, setcovercons, SCIPconsGetName(cons), &consdata->row) );
            assert(consdata->row != NULL);
            assert(!SCIProwIsInLP(consdata->row));
#ifndef NDEBUG
            {
               Real rowactivity;
               rowactivity = SCIPgetRowLPActivity(scip, consdata->row);
               assert(SCIPisSumEQ(scip, rowactivity, sum));
            }
#endif
            
            /* insert LP row as cut */
            CHECK_OKAY( SCIPaddCut(scip, consdata->row, (1.0-sum)/(nvars+1)) );
            CHECK_OKAY( SCIPresetConsAge(scip, cons) );
            *result = SCIP_SEPARATED;
         }
         else
         {
            /* constraint was feasible -> increase age */
            CHECK_OKAY( SCIPincConsAge(scip, cons) );
         }
      }
   }
   
   return SCIP_OKAY;
}

static
DECL_CONSSEPA(consSepaSetcover)
{
   CONSDATA* consdata;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   debugMessage("separating %d/%d set covering constraints\n", nusefulconss, nconss);

   *result = SCIP_DIDNOTFIND;

   /* step 1: check all useful set covering constraints for feasibility */
   for( c = 0; c < nusefulconss && *result != SCIP_CUTOFF; ++c )
   {
      CHECK_OKAY( separate(scip, conss[c], result) );
   }

   /* step 2: combine set covering constraints to get more cuts */
   todoMessage("further cuts of set covering constraints");

   return SCIP_OKAY;
}

/** if fractional variables exist, chooses a set S of them and branches on (i) x(S) >= 1, and (ii) x(S) == 0 */
static
RETCODE branchLP(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLR*        conshdlr,           /**< set covering constraint handler */
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

   todoMessage("use a better set covering branching on LP solution");

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
   
   /* sort fractional variables by number of uses in enabled set covering constraints */
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

   if( nsortcands == 0 )
   {
      /* none of the fractional variables is member of a set covering constraint
       * -> we are not responsible for doing the branching
       */
      return SCIP_OKAY;
   }

#if 0
   /* select the first variables from the sorted candidate list, until MINBRANCHWEIGHT is reached */
   branchweight = 0.0;
   for( nselcands = 0; nselcands < nsortcands && branchweight < MINBRANCHWEIGHT; ++nselcands )
   {
      CHECK_OKAY( SCIPgetVarSol(scip, sortcands[nselcands], &solval) );
      assert(SCIPisFeasGE(scip, solval, 0.0) && SCIPisFeasLE(scip, solval, 1.0));
      branchweight += solval;
   }
#else
   /* select the first variables from the sorted candidate list, until MAXBRANCHWEIGHT is reached; then choose one less */
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
#endif

   /* check, if we accumulated at most MAXBRANCHWEIGHT weight */
   if( MINBRANCHWEIGHT <= branchweight && branchweight <= MAXBRANCHWEIGHT )
   {
      NODE* node;

      /* perform the set covering branching on the selected variables */
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
         char name[255];
         
         /* add constraint x(S) >= 1 */
         sprintf(name, "SCB%lld", SCIPgetNodenum(scip));

         CHECK_OKAY( SCIPcreateConsSetcover(scip, &newcons, name, nselcands, sortcands,
                        TRUE, TRUE, FALSE, TRUE, FALSE, TRUE) );
         CHECK_OKAY( SCIPaddConsNode(scip, node, newcons) );
         CHECK_OKAY( SCIPreleaseCons(scip, &newcons) );
      }
      
      *result = SCIP_BRANCHED;
         
#ifdef DEBUG
      debugMessage("set covering branching: nselcands=%d/%d, weight(S)=%g, A={", nselcands, nlpcands, branchweight);
      for( i = 0; i < nselcands; ++i )
      {
         CHECK_OKAY( SCIPgetSolVal(scip, NULL, sortcands[i], &solval) );
         printf(" %s[%g]", SCIPvarGetName(sortcands[i]), solval);
      }
      printf(" }\n");
#endif
   }

   /* free temporary memory */
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &sortcands) );

   return SCIP_OKAY;
}

static
RETCODE branchPseudo(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLR*        conshdlr,           /**< set covering constraint handler */
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

   todoMessage("use a better set covering branching on pseudo solution");

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
   
   /* sort fractional variables by number of uses in enabled set covering constraints */
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
      /* none of the fractional variables is member of a set covering constraint
       * -> we are not responsible for doing the branching
       */
      return SCIP_OKAY;
   }
   
   /* branch on the first part of the sorted candidates:
    * - for each of these variables i, create a child node x_0 = ... = x_i-1 = 0, x_i = 1
    * - create an additional child node x_0 = ... = x_n = 0
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

#if DEBUG
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
DECL_CONSENFOLP(consEnfolpSetcover)
{
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   debugMessage("enforcing %d set covering constraints\n", nconss);

   *result = SCIP_FEASIBLE;

   /* step 1: check all useful set covering constraints for feasibility */
   for( c = 0; c < nusefulconss && *result != SCIP_CUTOFF; ++c )
   {
      CHECK_OKAY( separate(scip, conss[c], result) );
   }
   if( *result != SCIP_FEASIBLE )
      return SCIP_OKAY;

   /* step 2: if solution is not integral, choose a variable set to branch on */
   CHECK_OKAY( branchLP(scip, conshdlr, result) );
   if( *result != SCIP_FEASIBLE )
      return SCIP_OKAY;

   /* step 3: check all obsolete set covering constraints for feasibility */
   for( c = nusefulconss; c < nconss && *result == SCIP_FEASIBLE; ++c )
   {
      CHECK_OKAY( separate(scip, conss[c], result) );
   }

   return SCIP_OKAY;
}

static
DECL_CONSENFOPS(consEnfopsSetcover)
{
   CONSDATA* consdata;
   SETCOVERCONS* setcovercons;
   int bestbranchcons;
   int bestscore;
   int score;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   /* check all set covering constraints for feasibility */
   bestbranchcons = -1;
   bestscore = 0;
   for( c = 0; c < nconss && (*result == SCIP_FEASIBLE || *result == SCIP_INFEASIBLE); ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      setcovercons = consdata->setcovercons;
      assert(setcovercons != NULL);
      
      if( setcovercons->nfixedzeros == setcovercons->nvars )
      {
         if( setcovercons->modifiable )
         {
            /* the constraint is infeasible with the current variable set, but it may be feasible after pricing in
             * additional variables -> the LP has to be solved to price in new variables
             */
            *result = SCIP_SOLVELP;
         }
         else
         {
            /* the constraint cannot be feasible, because all variables are fixed to zero */
            *result = SCIP_CUTOFF;
         }
         CHECK_OKAY( SCIPresetConsAge(scip, conss[c]) );
      }
      else if( setcovercons->nfixedones >= 1 )
      {
         /* constraint is feasible anyway, because one of its variables if fixed to one -> constraint can be disabled */
         CHECK_OKAY( SCIPdisableConsLocal(scip, conss[c]) );
      }
      else if( setcovercons->nfixedzeros == setcovercons->nvars-1 && !setcovercons->modifiable )
      {
         VAR** vars;
         Real solval;
         int nvars;
         int foundvar;
         int v;
         
         /* all except one variable are fixed to zero: fix the remaining variable to one */

         /* search all variables of the constraint until the non-fixed one is found */
         vars = setcovercons->vars;
         nvars = setcovercons->nvars;
         foundvar = -1;
         for( v = 0; v < nvars && foundvar == -1; ++v )
         {
            assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY);
            if( SCIPisEQ(scip, SCIPvarGetUb(vars[v]), 1.0) )
               foundvar = v;
         }
         assert(0 <= foundvar && foundvar < nvars);

         /* fix variable to one */
         CHECK_OKAY( SCIPchgVarLb(scip, vars[foundvar], 1.0) );

         *result = SCIP_REDUCEDDOM;
      }
      else if( *result == SCIP_FEASIBLE )
      {
         VAR** vars;
         Real solval;
         Bool found;
         int nvars;
         int v;
         
         /* search all variables of the constraint until one is found with solution value 1.0 */
         vars = setcovercons->vars;
         nvars = setcovercons->nvars;
         found = FALSE;
         for( v = 0; v < nvars && !found; ++v )
         {
            assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY);
            solval = SCIPgetSolVal(scip, NULL, vars[v]);
            assert(SCIPisEQ(scip, solval, 0.0) || SCIPisEQ(scip, solval, 1.0));
            found = SCIPisEQ(scip, solval, 1.0);
         }
         if( !found )
         {
            /* constraint is violated by pseudo solution */
            CHECK_OKAY( SCIPresetConsAge(scip, conss[c]) );
            *result = SCIP_INFEASIBLE;
         }
      }
   }
   
   if( *result == SCIP_INFEASIBLE )
   {
      /* at least one constraint is violated by pseudo solution and we didn't find a better way to resolve this:
       * -> branch on pseudo solution
       */
      CHECK_OKAY( branchPseudo(scip, conshdlr, result) );
   }

   return SCIP_OKAY;
}

static
DECL_CONSCHECK(consCheckSetcover)
{
   CONSDATA* consdata;
   SETCOVERCONS* setcovercons;
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

   /* check all set covering constraints for feasibility */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      setcovercons = consdata->setcovercons;
      assert(setcovercons != NULL);
      if( checklprows || consdata->row == NULL || !SCIProwIsInLP(consdata->row) )
      {
         /* search all variables of the constraint until one is found with solution value 1.0 */
         vars = setcovercons->vars;
         nvars = setcovercons->nvars;
         found = FALSE;
         for( v = 0; v < nvars && !found; ++v )
         {
            assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY);
            solval = SCIPgetSolVal(scip, sol, vars[v]);
            assert(SCIPisEQ(scip, solval, 0.0) || SCIPisEQ(scip, solval, 1.0));
            found = SCIPisEQ(scip, solval, 1.0);
         }
         if( !found )
         {
            /* constraint is violated */
            *result = SCIP_INFEASIBLE;
            CHECK_OKAY( SCIPresetConsAge(scip, conss[c]) );
            return SCIP_OKAY;
         }
         else
         {
            CHECK_OKAY( SCIPincConsAge(scip, conss[c]) );
         }
      }
   }
   
   return SCIP_OKAY;
}

static
DECL_CONSENABLE(consEnableSetcover)
{
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* consdata;
   SETCOVERCONS* setcovercons;
   int v;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   setcovercons = consdata->setcovercons;
   assert(setcovercons != NULL);

   debugMessage("enabling information method of set covering constraint handler\n");

   /* increase the number of uses for each variable in the constraint */
   for( v = 0; v < setcovercons->nvars; ++v )
   {
      CHECK_OKAY( conshdlrdataIncVaruses(scip, conshdlrdata, setcovercons->vars[v]) );
   }

   return SCIP_OKAY;
}

static
DECL_CONSDISABLE(consDisableSetcover)
{
   CONSHDLRDATA* conshdlrdata;
   INTARRAY* varuses;
   CONSDATA* consdata;
   SETCOVERCONS* setcovercons;
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
   setcovercons = consdata->setcovercons;
   assert(setcovercons != NULL);

   debugMessage("disabling information method of set covering constraint handler\n");

   /* decrease the number of uses for each variable in the constraint */
   for( v = 0; v < setcovercons->nvars; ++v )
   {
      CHECK_OKAY( conshdlrdataDecVaruses(scip, conshdlrdata, setcovercons->vars[v]) );
   }

   return SCIP_OKAY;
}

static
DECL_LINCONSUPGD(linconsUpgdSetcover)
{
   assert(upgdcons != NULL);
   assert(upgraded != NULL);

   /* check, if linear constraint can be upgraded to set covering constraint
    * -> a set covering constraint consists only of binary variables with a coefficient of +1.0,
    *    and has left hand side of +1.0, and infinite right hand side
    */
   if( nposbin == nvars && ncoeffspone == nvars && SCIPisEQ(scip, lhs, 1.0) && SCIPisInfinity(scip, rhs) )
   {
      debugMessage("upgrading constraint <%s> to set covering constraint\n", SCIPconsGetName(cons));

      /* create the set covering constraint (an automatically upgraded constraint is always unmodifiable) */
      CHECK_OKAY( SCIPcreateConsSetcover(scip, upgdcons, SCIPconsGetName(cons), nvars, vars,
                     SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     local, FALSE, removeable) );
      *upgraded = TRUE;
   }
   else
      *upgraded = FALSE;

   return SCIP_OKAY;
}




/*
 * Callback methods of event handler
 */

static
DECL_EVENTEXEC(eventExecSetcover)
{
   SETCOVERCONS* setcovercons;
   EVENTTYPE eventtype;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   /*debugMessage("Exec method of bound change event handler for set covering constraints\n");*/

   setcovercons = (SETCOVERCONS*)eventdata;
   assert(setcovercons != NULL);

   CHECK_OKAY( SCIPeventGetType(event, &eventtype) );
   switch( eventtype )
   {
   case SCIP_EVENTTYPE_LBTIGHTENED:
      setcovercons->nfixedones++;
      break;
   case SCIP_EVENTTYPE_LBRELAXED:
      setcovercons->nfixedones--;
      break;
   case SCIP_EVENTTYPE_UBTIGHTENED:
      setcovercons->nfixedzeros++;
      break;
   case SCIP_EVENTTYPE_UBRELAXED:
      setcovercons->nfixedzeros--;
      break;
   default:
      errorMessage("invalid event type");
      abort();
   }
   assert(0 <= setcovercons->nfixedzeros && setcovercons->nfixedzeros <= setcovercons->nvars);
   assert(0 <= setcovercons->nfixedones && setcovercons->nfixedones <= setcovercons->nvars);

   debugMessage(" -> constraint has %d zero-fixed and %d one-fixed of %d variables\n", 
      setcovercons->nfixedzeros, setcovercons->nfixedones, setcovercons->nvars);

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for set covering constraint and includes it in SCIP */
RETCODE SCIPincludeConsHdlrSetcover(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create event handler for bound change events */
   CHECK_OKAY( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
                  NULL, NULL, NULL,
                  NULL, eventExecSetcover,
                  NULL) );

   /* create constraint handler data */
   CHECK_OKAY( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConsHdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
                  CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ,
                  CONSHDLR_NEEDSCONS,
                  consFreeSetcover, NULL, NULL,
                  consDeleteSetcover, consTransSetcover, 
                  consSepaSetcover, consEnfolpSetcover, consEnfopsSetcover, consCheckSetcover, NULL,
                  consEnableSetcover, consDisableSetcover,
                  conshdlrdata) );

   /* include the linear constraint to set covering constraint upgrade in the linear constraint handler */
   CHECK_OKAY( SCIPincludeLinconsUpgrade(scip, linconsUpgdSetcover, LINCONSUPGD_PRIORITY) );

   return SCIP_OKAY;
}

/** creates and captures a set covering constraint */
RETCODE SCIPcreateConsSetcover(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             local,              /**< is set covering constraint only valid locally? */
   Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   assert(scip != NULL);

   /* find the set covering constraint handler */
   conshdlr = SCIPfindConsHdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("set covering constraint handler not found");
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
      CHECK_OKAY( setcoverconsCreate(scip, &consdata->setcovercons, nvars, vars, local, modifiable, removeable) );
   }
   else
   {
      /* create constraint in transformed problem */
      CHECK_OKAY( setcoverconsCreateTransformed(scip, &consdata->setcovercons, nvars, vars, 
                     local, modifiable, removeable) );
   }
   consdata->row = NULL;

   /* create constraint (propagation is never used for set covering constraints) */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, separate, enforce, check, FALSE) );

   return SCIP_OKAY;
}


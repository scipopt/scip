/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_linear.c
 * @brief  constraint handler for linear constraints
 * @author Tobias Achterberg
 */

/** Linear constraints can be stored in two different ways: as LP rows with
 *  coefficients for the columns, or as LINCONS objects with coefficients for
 *  the variables. The second way is needed to be able to create the constraints
 *  before the solution process started, and before the necessary columns exist.
 *  At the first moment, a linear constraint is separated, it gets converted into
 *  an LP row.
 *
 *  Linear constraints are separated with a high priority, because they are easy
 *  to separate. The cut pool is implemented by adding linear constraints to the
 *  root node, such that it is separated each time, the linear constraints are
 *  separated. A constraint handler, which generates cuts for the pool should have
 *  a lower separation priority than the linear constraint handler, and it should
 *  have a separation frequency that is a multiple of the frequency of the linear
 *  constraint handler. In this way, it can be avoided to separate the same cut
 *  twice, because if a separation run of the handler is always preceded by a
 *  separation of the linear constraints, the pooled cuts are always satisfied.
 *
 *  Linear constraints are enforced and checked with a very low priority. Checking
 *  of (many) linear constraints is much more involved than checking the solution
 *  values for integrality. Because we are separating the linear constraints quite
 *  often, it is only necessary to enforce them for integral solutions. A constraint
 *  handler which generates pool cuts in its enforcing method should have an
 *  enforcing priority smaller than that of the linear constraint handler to avoid
 *  regenerating cuts which are already in the cut pool.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_linear.h"


typedef struct LinCons LINCONS;         /**< externally stored linear constraint */


#define CONSHDLR_NAME          "linear"
#define CONSHDLR_DESC          "Linear constraints of the form  lhs <= a^T x <= rhs"
#define CONSHDLR_SEPAPRIORITY  +1000000
#define CONSHDLR_ENFOPRIORITY  -1000000
#define CONSHDLR_CHCKPRIORITY  -1000000
#define CONSHDLR_NEEDSCONS         TRUE /**< the constraint handler should only be called, if linear constraints exist */


/** externally stored linear constraint */
struct LinCons
{
   VAR**            vars;               /**< variables of constraint entries */
   Real*            vals;               /**< coefficients of constraint entries */
   Real             lhs;                /**< left hand side of row (for ranged rows) */
   Real             rhs;                /**< right hand side of row */
   int              varssize;           /**< size of the vars- and vals-arrays */
   int              nvars;              /**< number of nonzeros in constraint */
   unsigned int     modifiable:1;       /**< is row modifiable during node processing (subject to column generation)? */
};

/** constraint data for linear constraints */
struct ConsData
{
   union
   {
      LINCONS*      lincons;            /**< external linear constraint data structure */
      ROW*          row;                /**< LP row, if constraint is already stored in LP row format */
   } data;
   unsigned int     islprow:1;          /**< is linear constraint already stored as LP row? */
};


/*
 * memory growing methods for dynamically allocated arrays
 */

static
RETCODE ensureVarsSize(                 /**< ensures, that vars and vals arrays can store at least num entries */
   LINCONS*         lincons,            /**< linear constraint data object */
   SCIP*            scip,               /**< SCIP data structure */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(lincons->nvars <= lincons->varssize);
   
   if( num > lincons->varssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      ALLOC_OKAY( SCIPreallocBlockMemoryArray(scip, lincons->vars, lincons->varssize, newsize) );
      ALLOC_OKAY( SCIPreallocBlockMemoryArray(scip, lincons->vals, lincons->varssize, newsize) );
      lincons->varssize = newsize;
   }
   assert(num <= lincons->varssize);

   return SCIP_OKAY;
}

/*
 * local methods
 */

static
RETCODE linconsCreate(                  /**< creates a linear constraint data object */
   LINCONS**        lincons,            /**< pointer to linear constraint data object */
   SCIP*            scip,               /**< SCIP data structure */
   int              nvars,              /**< number of nonzeros in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Real*            vals,               /**< array with coefficients of constraint entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
   Bool             modifiable          /**< is row modifiable during node processing (subject to column generation)? */
   )
{
   assert(lincons != NULL);
   assert(scip != NULL);

   if( SCIPisG(scip, lhs, rhs) )
   {
      char s[255];
      errorMessage("left hand side of linear constraint greater than right hand side");
      sprintf(s, "  (lhs=%f, rhs=%f)", lhs, rhs);
      errorMessage(s);
      return SCIP_INVALIDDATA;
   }

   ALLOC_OKAY( SCIPallocBlockMemory(scip, *lincons) );

   if( nvars > 0 )
   {
      ALLOC_OKAY( SCIPduplicateBlockMemoryArray(scip, (*lincons)->vars, vars, nvars) );
      ALLOC_OKAY( SCIPduplicateBlockMemoryArray(scip, (*lincons)->vals, vals, nvars) );
   }
   else
   {
      (*lincons)->vars = NULL;
      (*lincons)->vals = NULL;
   }

   (*lincons)->lhs = lhs;
   (*lincons)->rhs = rhs;
   (*lincons)->varssize = nvars;
   (*lincons)->nvars = nvars;
   (*lincons)->modifiable = modifiable;

   return SCIP_OKAY;
}

static
void linconsFree(                       /**< frees a linear constraint data object */
   LINCONS**        lincons,            /**< pointer to linear constraint data object */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(lincons != NULL);
   assert(*lincons != NULL);
   assert((*lincons)->varssize >= 0);

   SCIPfreeBlockMemoryArrayNull(scip, (*lincons)->vars, (*lincons)->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, (*lincons)->vals, (*lincons)->varssize);
   SCIPfreeBlockMemory(scip, *lincons);
}

static
RETCODE linconsAddCoef(                 /**< adds coefficient in linear constraint */
   LINCONS*         lincons,            /**< linear constraint data object */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable of constraint entry */
   Real             val                 /**< coefficients of constraint entry */
   )
{
   assert(lincons != NULL);
   assert(scip != NULL);
   assert(var != NULL);

   CHECK_OKAY( ensureVarsSize(lincons, scip, lincons->nvars+1) );
   lincons->vars[lincons->nvars] = var;
   lincons->vals[lincons->nvars] = val;
   lincons->nvars++;

   return SCIP_OKAY;
}

static
Real linconsGetActivity(                /**< calculates the activity for the linear constraint */
   LINCONS*         lincons             /**< linear constraint data object */
   )
{
   Real activity;
   int v;

   assert(lincons != NULL);
   
   activity = 0.0;
   for( v = 0; v < lincons->nvars; ++v )
   {
      activity += SCIPvarGetPrimsol(lincons->vars[v]) * lincons->vals[v];
   }

   return activity;
}

static
Real linconsGetFeasibility(             /**< calculates the feasibility for the linear constraint */
   LINCONS*         lincons             /**< linear constraint data object */
   )
{
   Real activity;

   assert(lincons != NULL);

   activity = linconsGetActivity(lincons);

   return MIN(lincons->rhs - activity, activity - lincons->lhs);
}

static
RETCODE linconsToRow(                   /**< creates an LP row from a linear constraint data object */
   LINCONS*         lincons,            /**< linear constraint data object */
   const char*      name,               /**< name of the constraint */
   SCIP*            scip,               /**< SCIP data structure */
   ROW**            row                 /**< pointer to an LP row data object */
   )
{
   int v;

   assert(lincons != NULL);
   assert(name != NULL);
   assert(scip != NULL);
   assert(row != NULL);

   CHECK_OKAY( SCIPcreateRow(scip, row, name, 0, NULL, NULL, lincons->lhs, lincons->rhs, lincons->modifiable) );
   
   for( v = 0; v < lincons->nvars; ++v )
   {
      CHECK_OKAY( SCIPaddVarToRow(scip, *row, lincons->vars[v], lincons->vals[v]) );
   }

   return SCIP_OKAY;
}

static
RETCODE applyConstraints(               /**< separates violated inequalities; called from Sepa and Enfo methods */
   CONSHDLR*        conshdlr,           /**< linear constraint handler */
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           conss,              /**< array of constraints to process */
   int              nconss,             /**< number of constraints to process */
   Bool*            found               /**< pointer to store information, if a violated constraint has been found */
   )
{
   CONS* cons;
   CONSDATA* consdata;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(nconss == 0 || conss != NULL);
   assert(nconss >= 0);
   assert(found != NULL);

   *found = FALSE;

   debugMessage("applying %d linear constraints at %p\n", nconss, conss);
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert(cons != NULL);

      debugMessage("applying linear constraint <%s>\n", SCIPconsGetName(cons));
      consdata = SCIPconsGetConsdata(cons);
      assert(consdata != NULL);

      if( consdata->islprow )
      {
         ROW* row = consdata->data.row;
         
         if( !SCIProwIsInLP(row) )
         {         
            Real feasibility;

            CHECK_OKAY( SCIPgetRowFeasibility(scip, row, &feasibility) );
            debugMessage("  row feasibility = %g\n", feasibility);
            if( !SCIPisFeasible(scip, feasibility) )
            {
               /* insert LP row as cut */
               CHECK_OKAY( SCIPaddCut(scip, row, -feasibility/SCIProwGetNorm(row)/(SCIProwGetNNonz(row)+1)) );
               *found = TRUE;
            }
         }
      }
      else
      {
         LINCONS* lincons = consdata->data.lincons;
         Real feasibility;
         
         /* check feasibility of linear constraint; if it's an equality, add the cut in any case */
         feasibility = linconsGetFeasibility(lincons);
         debugMessage("  lincons feasibility = %g\n", feasibility);
         if( !SCIPisFeasible(scip, feasibility) || SCIPisEQ(scip, lincons->lhs, lincons->rhs) )
         {
            ROW* row;
            
            /* convert lincons data into LP row */
            CHECK_OKAY( linconsToRow(lincons, SCIPconsGetName(cons), scip, &row) );
#ifndef NDEBUG
            {
               Real rowfeasibility;
               CHECK_OKAY( SCIPgetRowFeasibility(scip, row, &rowfeasibility) );
               assert(SCIPisEQ(scip, rowfeasibility, feasibility));
            }
#endif

            /* free the lincons data and convert consdata to point to the row */
            linconsFree(&consdata->data.lincons, scip);
            consdata->islprow = TRUE;
            consdata->data.row = row;
            
            /* don't release the row, because we need it as data storage */

            /* insert LP row as cut */
            CHECK_OKAY( SCIPaddCut(scip, row, -feasibility/SCIProwGetNorm(row)/(SCIProwGetNNonz(row)+1)) );
            *found = TRUE;
         }
      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods
 */

static
DECL_CONSINIT(SCIPconsFreeLinear)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   infoMessage(SCIPverbLevel(scip), SCIP_VERBLEVEL_FULL, "free linear constraint handler");

   return SCIP_OKAY;
}

static
DECL_CONSINIT(SCIPconsInitLinear)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   infoMessage(SCIPverbLevel(scip), SCIP_VERBLEVEL_FULL, "initialise linear constraint handler");

   return SCIP_OKAY;
}

static
DECL_CONSEXIT(SCIPconsExitLinear)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   infoMessage(SCIPverbLevel(scip), SCIP_VERBLEVEL_FULL, "exit linear constraint handler");

   return SCIP_OKAY;
}

static
DECL_CONSDELE(SCIPconsDeleLinear)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   if( (*consdata)->islprow )
   {
      CHECK_OKAY( SCIPreleaseRow(scip, &(*consdata)->data.row) );
   }
   else
   {
      linconsFree(&(*consdata)->data.lincons, scip);
   }
   SCIPfreeBlockMemory(scip, *consdata);

   return SCIP_OKAY;
}

static
DECL_CONSTRAN(SCIPconsTranLinear)
{
   LINCONS* lincons;

   debugMessage("Tran method of linear constraints\n");

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(SCIPstage(scip) == SCIP_STAGE_INITSOLVE);
   assert(sourcedata != NULL);
   assert(!sourcedata->islprow);  /* in original problem, there cannot be LP rows */
   assert(sourcedata->data.lincons != NULL);
   assert(targetdata != NULL);

   lincons = sourcedata->data.lincons;

   ALLOC_OKAY( SCIPallocBlockMemory(scip, *targetdata) );
   (*targetdata)->islprow = FALSE;

   CHECK_OKAY( linconsCreate(&(*targetdata)->data.lincons, scip, lincons->nvars, lincons->vars, lincons->vals,
                  lincons->lhs, lincons->rhs, lincons->modifiable) );
   
   return SCIP_OKAY;
}

static
DECL_CONSSEPA(SCIPconsSepaLinear)
{
   Bool found;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   debugMessage("Sepa method of linear constraints\n");

   /* check for violated constraints */
   CHECK_OKAY( applyConstraints(conshdlr, scip, conss, nconss, &found) );

   if( found )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

static
DECL_CONSENFO(SCIPconsEnfoLinear)
{
   Bool found;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   debugMessage("Enfo method of linear constraints\n");

   /* check for violated constraints */
   CHECK_OKAY( applyConstraints(conshdlr, scip, conss, nconss, &found) );

   if( found )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

static
DECL_CONSCHCK(SCIPconsChckLinear)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   todoMessage("Chck method of linear constraints");

   return SCIP_OKAY;
}

static
DECL_CONSPROP(SCIPconsPropLinear)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   todoMessage("Prop method of linear constraints");

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

RETCODE SCIPincludeConsHdlrLinear(      /**< creates the handler for linear constraints and includes it in SCIP */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( SCIPincludeConsHdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHCKPRIORITY, CONSHDLR_NEEDSCONS,
                  SCIPconsFreeLinear, SCIPconsInitLinear, SCIPconsExitLinear, 
                  SCIPconsDeleLinear, SCIPconsTranLinear, 
                  SCIPconsSepaLinear, SCIPconsEnfoLinear, SCIPconsChckLinear, SCIPconsPropLinear,
                  NULL) );

   return SCIP_OKAY;
}

RETCODE SCIPcreateConsLinear(           /**< creates and captures a linear constraint */
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              len,                /**< number of nonzeros in the constraint */
   VAR**            var,                /**< array with variables of constraint entries */
   Real*            val,                /**< array with coefficients of constraint entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
   Bool             model,              /**< is constraint necessary for feasibility? */
   Bool             modifiable          /**< is row modifiable during node processing (subject to column generation)? */
   )
{
   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   assert(scip != NULL);

   /* find the linear constraint handler */
   CHECK_OKAY( SCIPfindConsHdlr(scip, CONSHDLR_NAME, &conshdlr) );
   if( conshdlr == NULL )
   {
      errorMessage("Linear constraint handler not found");
      return SCIP_INVALIDCALL;
   }

   /* create the constraint specific data */
   ALLOC_OKAY( SCIPallocBlockMemory(scip, consdata) );
   consdata->islprow = FALSE;
   CHECK_OKAY( linconsCreate(&consdata->data.lincons, scip, len, var, val, lhs, rhs, modifiable) );

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, model) );

   return SCIP_OKAY;
}

RETCODE SCIPcreateConsLPRow(            /**< creates and captures a linear constraint from an LP row, captures the row */
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   ROW*             row,                /**< LP row */
   Bool             model               /**< is constraint necessary for feasibility? */
   )
{
   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   assert(scip != NULL);

   /* find the linear constraint handler */
   CHECK_OKAY( SCIPfindConsHdlr(scip, CONSHDLR_NAME, &conshdlr) );
   if( conshdlr == NULL )
   {
      errorMessage("Linear constraint handler not found");
      return SCIP_INVALIDCALL;
   }

   /* create the constraint specific data */
   ALLOC_OKAY( SCIPallocBlockMemory(scip, consdata) );
   consdata->islprow = TRUE;
   consdata->data.row = row;

   /* capture the row, because we need it as data storage */
   SCIPcaptureRow(scip, row);

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, SCIProwGetName(row), conshdlr, consdata, model) );

   return SCIP_OKAY;
}

RETCODE SCIPconsLinearAddCoef(          /**< adds coefficient in linear constraint */
   CONS*            cons,               /**< constraint data */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable of constraint entry */
   Real             val                 /**< coefficients of constraint entry */
   )
{
   CONSDATA* consdata;

   assert(cons != NULL);
   assert(scip != NULL);
   assert(var != NULL);

   /*debugMessage("adding coefficient %g * <%s> to linear constraint <%s>\n", val, var->name, SCIPconsGetName(cons));*/

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetConsHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not linear");
      return SCIP_INVALIDDATA;
   }
   
   consdata = SCIPconsGetConsdata(cons);
   if( consdata->islprow )
   {
      assert(consdata->data.row != NULL);
      CHECK_OKAY( SCIPaddVarToRow(scip, consdata->data.row, var, val) );
   }
   else
   {
      assert(consdata->data.lincons != NULL);
      CHECK_OKAY( linconsAddCoef(consdata->data.lincons, scip, var, val) );
   }

   return SCIP_OKAY;
}

RETCODE SCIPconsLinearGetLhs(           /**< gets left hand side of linear constraint */
   CONS*            cons,               /**< constraint data */
   SCIP*            scip,               /**< SCIP data structure */
   Real*            lhs                 /**< pointer to store left hand side */
   )
{
   CONSDATA* consdata;

   assert(cons != NULL);
   assert(scip != NULL);
   assert(lhs != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetConsHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not linear");
      return SCIP_INVALIDDATA;
   }
   
   consdata = SCIPconsGetConsdata(cons);
   if( consdata->islprow )
   {
      assert(consdata->data.row != NULL);
      *lhs = SCIProwGetLhs(consdata->data.row);
   }
   else
   {
      assert(consdata->data.lincons != NULL);
      *lhs = consdata->data.lincons->lhs;
   }

   return SCIP_OKAY;
}

RETCODE SCIPconsLinearGetRhs(           /**< gets right hand side of linear constraint */
   CONS*            cons,               /**< constraint data */
   SCIP*            scip,               /**< SCIP data structure */
   Real*            rhs                 /**< pointer to store right hand side */
   )
{
   CONSDATA* consdata;

   assert(cons != NULL);
   assert(scip != NULL);
   assert(rhs != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetConsHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not linear");
      return SCIP_INVALIDDATA;
   }
   
   consdata = SCIPconsGetConsdata(cons);
   if( consdata->islprow )
   {
      assert(consdata->data.row != NULL);
      *rhs = SCIProwGetRhs(consdata->data.row);
   }
   else
   {
      assert(consdata->data.lincons != NULL);
      *rhs = consdata->data.lincons->rhs;
   }

   return SCIP_OKAY;
}

RETCODE SCIPconsLinearChgLhs(           /**< changes left hand side of linear constraint */
   CONS*            cons,               /**< constraint data */
   SCIP*            scip,               /**< SCIP data structure */
   Real             lhs                 /**< new left hand side */
   )
{
   CONSDATA* consdata;

   assert(cons != NULL);
   assert(scip != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetConsHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not linear");
      return SCIP_INVALIDDATA;
   }
   
   consdata = SCIPconsGetConsdata(cons);
   if( consdata->islprow )
   {
      errorMessage("cannot change sides of linear constraint already stored as LP row");
      return SCIP_INVALIDDATA;
   }
   else
   {
      assert(consdata->data.lincons != NULL);
      consdata->data.lincons->lhs = lhs;
   }

   return SCIP_OKAY;
}

RETCODE SCIPconsLinearChgRhs(           /**< changes right hand side of linear constraint */
   CONS*            cons,               /**< constraint data */
   SCIP*            scip,               /**< SCIP data structure */
   Real             rhs                 /**< new right hand side */
   )
{
   CONSDATA* consdata;

   assert(cons != NULL);
   assert(scip != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetConsHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not linear");
      return SCIP_INVALIDDATA;
   }
   
   consdata = SCIPconsGetConsdata(cons);
   if( consdata->islprow )
   {
      errorMessage("cannot change sides of linear constraint already stored as LP row");
      return SCIP_INVALIDDATA;
   }
   else
   {
      assert(consdata->data.lincons != NULL);
      consdata->data.lincons->rhs = rhs;
   }

   return SCIP_OKAY;
}

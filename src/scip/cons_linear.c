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
#define CONSHDLR_PROPFREQ             3
#define CONSHDLR_NEEDSCONS         TRUE /**< the constraint handler should only be called, if linear constraints exist */

#define EVENTHDLR_NAME         "linear"
#define EVENTHDLR_DESC         "bound change event handler for linear constraints"


/** externally stored linear constraint */
struct LinCons
{
   VAR**            vars;               /**< variables of constraint entries */
   Real*            vals;               /**< coefficients of constraint entries */
   EVENTDATA**      eventdatas;         /**< event datas for bound change events of the variables */
   Real             lhs;                /**< left hand side of row (for ranged rows) */
   Real             rhs;                /**< right hand side of row */
   Real             minactivity;        /**< minimal value w.r.t. the variable's bounds for the constraint's activity */
   Real             maxactivity;        /**< maximal value w.r.t. the variable's bounds for the constraint's activity */
   int              varssize;           /**< size of the vars- and vals-arrays */
   int              nvars;              /**< number of nonzeros in constraint */
   unsigned int     modifiable:1;       /**< is row modifiable during node processing (subject to column generation)? */
   unsigned int     transformed:1;      /**< does the linear constraint data belongs to the transformed problem? */
   unsigned int     validactivitybds:1; /**< are the activity bounds minactivity/maxactivity valid? */
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

/** constraint handler data */
struct ConsHdlrData
{
   EVENTHDLR*       eventhdlr;          /**< event handler to process bound change events with */
};

/** event data for bound change event */
struct EventData
{
   LINCONS*         lincons;            /**< linear constraint to process the bound change for */
   int              varpos;             /**< position of variable in vars array */
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
      if( lincons->transformed )
      {
         ALLOC_OKAY( SCIPreallocBlockMemoryArray(scip, lincons->eventdatas, lincons->varssize, newsize) );
      }
      else
         assert(lincons->eventdatas == NULL);
      lincons->varssize = newsize;
   }
   assert(num <= lincons->varssize);

   return SCIP_OKAY;
}

/*
 * local methods
 */

static
RETCODE linconsCatchEvent(              /**< creates event data for variable at given position, and catches events */
   LINCONS*         lincons,            /**< linear constraint data object */
   SCIP*            scip,               /**< SCIP data structure */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< array position of variable to catch bound change events for */
   )
{
   assert(lincons != NULL);
   assert(scip != NULL);
   assert(0 <= pos && pos < lincons->nvars);
   assert(lincons->vars != NULL);
   assert(lincons->vars[pos] != NULL);
   assert(lincons->eventdatas != NULL);
   assert(lincons->eventdatas[pos] == NULL);

   ALLOC_OKAY( SCIPallocBlockMemory(scip, lincons->eventdatas[pos]) );
   lincons->eventdatas[pos]->lincons = lincons;
   lincons->eventdatas[pos]->varpos = pos;

   CHECK_OKAY( SCIPcatchVarEvent(scip, lincons->vars[pos], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, 
                  lincons->eventdatas[pos]) );

   return SCIP_OKAY;
}

static
RETCODE linconsDropEvent(               /**< deletes event data for variable at given position, and drops events */
   LINCONS*         lincons,            /**< linear constraint data object */
   SCIP*            scip,               /**< SCIP data structure */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< array position of variable to catch bound change events for */
   )
{
   assert(lincons != NULL);
   assert(scip != NULL);
   assert(0 <= pos && pos < lincons->nvars);
   assert(lincons->vars[pos] != NULL);
   assert(lincons->eventdatas[pos] != NULL);
   assert(lincons->eventdatas[pos]->lincons == lincons);
   assert(lincons->eventdatas[pos]->varpos == pos);
   
   CHECK_OKAY( SCIPdropVarEvent(scip, lincons->vars[pos], eventhdlr, lincons->eventdatas[pos]) );

   SCIPfreeBlockMemory(scip, lincons->eventdatas[pos]);

   return SCIP_OKAY;
}

static
RETCODE linconsCreate(                  /**< creates a linear constraint data object of the original problem */
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
   int i;

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
   (*lincons)->eventdatas = NULL;

   (*lincons)->lhs = lhs;
   (*lincons)->rhs = rhs;
   (*lincons)->minactivity = SCIP_INVALID;
   (*lincons)->maxactivity = SCIP_INVALID;
   (*lincons)->varssize = nvars;
   (*lincons)->nvars = nvars;
   (*lincons)->modifiable = modifiable;
   (*lincons)->transformed = FALSE;
   (*lincons)->validactivitybds = FALSE;

   return SCIP_OKAY;
}

static
RETCODE linconsCreateTransformed(       /**< creates a linear constraint data object of the transformed problem */
   LINCONS**        lincons,            /**< pointer to linear constraint data object */
   SCIP*            scip,               /**< SCIP data structure */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              nvars,              /**< number of nonzeros in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Real*            vals,               /**< array with coefficients of constraint entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
   Bool             modifiable          /**< is row modifiable during node processing (subject to column generation)? */
   )
{
   int i;

   assert(lincons != NULL);

   CHECK_OKAY( linconsCreate(lincons, scip, nvars, vars, vals, lhs, rhs, modifiable) );

   /* allocate the additional needed eventdatas array */
   assert((*lincons)->eventdatas == NULL);
   ALLOC_OKAY( SCIPallocBlockMemoryArray(scip, (*lincons)->eventdatas, (*lincons)->varssize) );
   (*lincons)->transformed = TRUE;

   /* transform the variables, catch events */
   for( i = 0; i < (*lincons)->nvars; ++i )
   {
      if( SCIPvarGetStatus((*lincons)->vars[i]) == SCIP_VARSTATUS_ORIGINAL )
      {
         (*lincons)->vars[i] = SCIPvarGetTransformed((*lincons)->vars[i]);
         assert((*lincons)->vars[i] != NULL);
      }
      assert(SCIPvarGetStatus((*lincons)->vars[i]) != SCIP_VARSTATUS_ORIGINAL);
      (*lincons)->eventdatas[i] = NULL;
      CHECK_OKAY( linconsCatchEvent(*lincons, scip, eventhdlr, i) );
   }

   return SCIP_OKAY;
}

static
RETCODE linconsFree(                    /**< frees a linear constraint data object */
   LINCONS**        lincons,            /**< pointer to linear constraint data object */
   SCIP*            scip,               /**< SCIP data structure */
   EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(lincons != NULL);
   assert(*lincons != NULL);
   assert((*lincons)->varssize >= 0);

   /* drop events for included variables */
   if( (*lincons)->transformed )
   {
      for( i = 0; i < (*lincons)->nvars; ++i )
      {
         CHECK_OKAY( linconsDropEvent(*lincons, scip, eventhdlr, i) );
      }
      SCIPfreeBlockMemoryArrayNull(scip, (*lincons)->eventdatas, (*lincons)->varssize);
   }
   assert((*lincons)->eventdatas == NULL);

   SCIPfreeBlockMemoryArrayNull(scip, (*lincons)->vars, (*lincons)->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, (*lincons)->vals, (*lincons)->varssize);
   SCIPfreeBlockMemory(scip, *lincons);

   return SCIP_OKAY;
}

static
RETCODE linconsAddCoef(                 /**< adds coefficient in linear constraint */
   LINCONS*         lincons,            /**< linear constraint data object */
   SCIP*            scip,               /**< SCIP data structure */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   VAR*             var,                /**< variable of constraint entry */
   Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(lincons != NULL);
   assert(scip != NULL);
   assert(var != NULL);

   if( lincons->transformed && SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
   {
      var = SCIPvarGetTransformed(var);
      assert(var != NULL);
   }

   assert(lincons->transformed ^ (SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL));

   CHECK_OKAY( ensureVarsSize(lincons, scip, lincons->nvars+1) );
   lincons->vars[lincons->nvars] = var;
   lincons->vals[lincons->nvars] = val;
   lincons->nvars++;

   if( lincons->transformed )
   {
      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL);
      lincons->eventdatas[lincons->nvars] = NULL;
      CHECK_OKAY( linconsCatchEvent(lincons, scip, eventhdlr, lincons->nvars-1) );
   }

   return SCIP_OKAY;
}

static
RETCODE linconsGetActivity(             /**< calculates the activity of the linear constraint for LP solution */
   LINCONS*         lincons,            /**< linear constraint data object */
   SCIP*            scip,               /**< SCIP data structure */
   Real*            activity            /**< pointer to store the activity */
   )
{
   Real solval;
   int v;

   assert(lincons != NULL);
   assert(lincons->transformed);
   assert(activity != NULL);

   *activity = 0.0;
   for( v = 0; v < lincons->nvars; ++v )
   {
      CHECK_OKAY( SCIPgetVarSol(scip, lincons->vars[v], &solval) );
      *activity += lincons->vals[v] * solval;
   }

   return SCIP_OKAY;
}

static
RETCODE linconsGetFeasibility(          /**< calculates the feasibility of the linear constraint for LP solution */
   LINCONS*         lincons,            /**< linear constraint data object */
   SCIP*            scip,               /**< SCIP data structure */
   Real*            feasibility         /**< pointer to store the feasibility */
   )
{
   Real activity;

   assert(lincons != NULL);
   assert(lincons->transformed);
   assert(feasibility != NULL);

   CHECK_OKAY( linconsGetActivity(lincons, scip, &activity) );
   *feasibility = MIN(lincons->rhs - activity, activity - lincons->lhs);

   return SCIP_OKAY;
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
   assert(lincons->transformed);
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
void linconsUpdateChgLb(                /**< updates minimum and maximum activity for a change in lower bound */
   LINCONS*         lincons,            /**< linear constraint data object */
   Real             oldlb,              /**< old lower bound of variable */
   Real             newlb,              /**< new lower bound of variable */
   Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(lincons != NULL);
   assert(lincons->transformed);

   if( lincons->validactivitybds )
   {
      assert(lincons->minactivity < SCIP_INVALID);
      assert(lincons->maxactivity < SCIP_INVALID);
      if( val > 0.0 )
         lincons->minactivity += val * (newlb - oldlb);
      else
         lincons->maxactivity += val * (newlb - oldlb);
   }
}

static
void linconsUpdateChgUb(                /**< updates minimum and maximum activity for a change in upper bound */
   LINCONS*         lincons,            /**< linear constraint data object */
   Real             oldub,              /**< old upper bound of variable */
   Real             newub,              /**< new upper bound of variable */
   Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(lincons != NULL);
   assert(lincons->transformed);

   if( lincons->validactivitybds )
   {
      assert(lincons->minactivity < SCIP_INVALID);
      assert(lincons->maxactivity < SCIP_INVALID);
      if( val > 0.0 )
         lincons->maxactivity += val * (newub - oldub);
      else
         lincons->minactivity += val * (newub - oldub);
   }
}

static
void linconsUpdateAddVar(               /**< updates minimum and maximum activity for variable addition */
   LINCONS*         lincons,            /**< linear constraint data object */
   VAR*             var,                /**< variable of constraint entry */
   Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(lincons != NULL);
   assert(lincons->transformed);

   if( lincons->validactivitybds )
   {
      assert(lincons->minactivity < SCIP_INVALID);
      assert(lincons->maxactivity < SCIP_INVALID);

      linconsUpdateChgLb(lincons, 0.0, SCIPvarGetLb(var), val);
      linconsUpdateChgUb(lincons, 0.0, SCIPvarGetUb(var), val);
   }
}

static
RETCODE consdataGetLhs(                 /**< gets left hand side of linear constraint */
   CONSDATA*        consdata,           /**< constraint data */
   SCIP*            scip,               /**< SCIP data structure */
   Real*            lhs                 /**< pointer to store left hand side */
   )
{
   assert(consdata != NULL);
   assert(scip != NULL);
   assert(lhs != NULL);

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

static
RETCODE consdataGetRhs(                 /**< gets right hand side of linear constraint */
   CONSDATA*        consdata,           /**< constraint data */
   SCIP*            scip,               /**< SCIP data structure */
   Real*            rhs                 /**< pointer to store right hand side */
   )
{
   assert(consdata != NULL);
   assert(scip != NULL);
   assert(rhs != NULL);

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

static
void linconsCalcActivityBounds(         /**< calculates minimum and maximum activity for constraint */
   LINCONS*         lincons             /**< linear constraint data object */
   )
{
   int i;
   
   assert(lincons != NULL);
   assert(lincons->transformed);
   assert(!lincons->validactivitybds);
   assert(lincons->minactivity >= SCIP_INVALID);
   assert(lincons->maxactivity >= SCIP_INVALID);
   
   lincons->validactivitybds = TRUE;
   lincons->minactivity = 0.0;
   lincons->maxactivity = 0.0;
   
   for( i = 0; i < lincons->nvars; ++ i )
      linconsUpdateAddVar(lincons, lincons->vars[i], lincons->vals[i]);
}

static
RETCODE consdataGetActivityBounds(      /**< gets activity bounds for constraint */
   CONSDATA*        consdata,           /**< constraint data */
   SCIP*            scip,               /**< SCIP data structure */
   Real*            minactivity,        /**< pointer to store the minimal activity */
   Real*            maxactivity         /**< pointer to store the maximal activity */
   )
{
   assert(consdata != NULL);
   assert(scip != NULL);
   assert(minactivity != NULL);
   assert(maxactivity != NULL);

   if( consdata->islprow )
   {
      CHECK_OKAY( SCIPgetRowActivityBounds(scip, consdata->data.row, minactivity, maxactivity) );
   }
   else
   {
      LINCONS* lincons;

      lincons = consdata->data.lincons;
      assert(lincons != NULL);
      if( !lincons->validactivitybds )
         linconsCalcActivityBounds(lincons);
      assert(lincons->minactivity < SCIP_INVALID);
      assert(lincons->maxactivity < SCIP_INVALID);
      *minactivity = lincons->minactivity;
      *maxactivity = lincons->maxactivity;
   }

   return SCIP_OKAY;
}

static
RETCODE separateConstraints(            /**< separates violated inequalities; called from Sepa and Enfo methods */
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

   /*debugMessage("separating %d linear constraints at %p\n", nconss, conss);*/
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert(cons != NULL);

      /*debugMessage("separating linear constraint <%s>\n", SCIPconsGetName(cons));*/
      consdata = SCIPconsGetConsData(cons);
      assert(consdata != NULL);

      if( consdata->islprow )
      {
         ROW* row = consdata->data.row;
         
         if( !SCIProwIsInLP(row) )
         {         
            Real feasibility;

            CHECK_OKAY( SCIPgetRowFeasibility(scip, row, &feasibility) );
            /*debugMessage("  row feasibility = %g\n", feasibility);*/
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
         
         /* check feasibility of linear constraint */
         CHECK_OKAY( linconsGetFeasibility(lincons, scip, &feasibility) );
         /*debugMessage("  lincons feasibility = %g\n", feasibility);*/
         if( !SCIPisFeasible(scip, feasibility) )
         {
            CONSHDLRDATA* conshdlrdata;
            ROW* row;
            
            /* convert lincons data into LP row */
            CHECK_OKAY( linconsToRow(lincons, SCIPconsGetName(cons), scip, &row) );
#ifndef NDEBUG
            {
               Real rowfeasibility;
               CHECK_OKAY( SCIPgetRowFeasibility(scip, row, &rowfeasibility) );
               assert(SCIPisRelEQ(scip, rowfeasibility, feasibility));
            }
#endif
            
            /* get the constraint handler data to access the event handler */
            conshdlrdata = SCIPconshdlrGetData(conshdlr);
            assert(conshdlrdata != NULL);
            assert(conshdlrdata->eventhdlr != NULL);

            /* free the lincons data and convert consdata to point to the row */
            CHECK_OKAY( linconsFree(&consdata->data.lincons, scip, conshdlrdata->eventhdlr) );
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

static
RETCODE checkConstraints(               /**< checks pseudo solution for violated inequalities */
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

   /*debugMessage("checking pseudo solution for %d linear constraints at %p\n", nconss, conss);*/
   for( c = 0; c < nconss && !(*found); ++c )
   {
      cons = conss[c];
      assert(cons != NULL);

      /*debugMessage("checking linear constraint <%s>\n", SCIPconsGetName(cons));*/
      consdata = SCIPconsGetConsData(cons);
      assert(consdata != NULL);

      if( consdata->islprow )
      {
         ROW* row = consdata->data.row;
         Real feasibility;

         CHECK_OKAY( SCIPgetRowPseudoFeasibility(scip, row, &feasibility) );
         /*debugMessage("  row feasibility = %g\n", feasibility);*/
         *found |= !SCIPisFeasible(scip, feasibility);
      }
      else
      {
         LINCONS* lincons = consdata->data.lincons;
         Real feasibility;
         
         /* check feasibility of linear constraint */
         CHECK_OKAY( linconsGetFeasibility(lincons, scip, &feasibility) );
         /*debugMessage("  lincons feasibility = %g\n", feasibility);*/
         *found |= !SCIPisFeasible(scip, feasibility);
      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

static
DECL_CONSFREE(consFreeLinear)
{
   CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   /* free the event handler for bound change events */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   CHECK_OKAY( SCIPeventhdlrFree(&conshdlrdata->eventhdlr, scip) );

   /* free constraint handler data */
   freeMemory(conshdlrdata);
   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}

static
DECL_CONSDELE(consDeleLinear)
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
      CONSHDLRDATA* conshdlrdata;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      CHECK_OKAY( linconsFree(&(*consdata)->data.lincons, scip, conshdlrdata->eventhdlr) );
   }
   SCIPfreeBlockMemory(scip, *consdata);

   return SCIP_OKAY;
}

static
DECL_CONSTRAN(consTranLinear)
{
   CONSHDLRDATA* conshdlrdata;
   LINCONS* lincons;

   /*debugMessage("Tran method of linear constraints\n");*/

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(SCIPstage(scip) == SCIP_STAGE_INITSOLVE);
   assert(sourcedata != NULL);
   assert(!sourcedata->islprow);  /* in original problem, there cannot be LP rows */
   assert(sourcedata->data.lincons != NULL);
   assert(targetdata != NULL);

   ALLOC_OKAY( SCIPallocBlockMemory(scip, *targetdata) );
   (*targetdata)->islprow = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   lincons = sourcedata->data.lincons;
   CHECK_OKAY( linconsCreateTransformed(&(*targetdata)->data.lincons, scip, conshdlrdata->eventhdlr, lincons->nvars, 
                  lincons->vars, lincons->vals, lincons->lhs, lincons->rhs, lincons->modifiable) );

   return SCIP_OKAY;
}

static
DECL_CONSSEPA(consSepaLinear)
{
   Bool found;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   /*debugMessage("Sepa method of linear constraints\n");*/

   /* separate violated constraints */
   CHECK_OKAY( separateConstraints(conshdlr, scip, conss, nconss, &found) );

   if( found )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

static
DECL_CONSENLP(consEnlpLinear)
{
   Bool found;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   /*debugMessage("Enlp method of linear constraints\n");*/

   /* check for violated constraints */

   /* LP is processed at current node -> we can add violated linear constraints to the LP */
   CHECK_OKAY( separateConstraints(conshdlr, scip, conss, nconss, &found) );
   
   if( found )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

static
DECL_CONSENPS(consEnpsLinear)
{
   Bool found;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   /*debugMessage("Enps method of linear constraints\n");*/

   /* check for violated constraints */

   /* LP is not processed at current node -> we just have to check pseudo solution for feasibility */
   CHECK_OKAY( checkConstraints(conshdlr, scip, conss, nconss, &found) );
   
   if( found )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

static
DECL_CONSCHCK(consChckLinear)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   todoMessage("Chck method of linear constraints");

   return SCIP_OKAY;
}

static
DECL_CONSPROP(consPropLinear)
{
   CONS* cons;
   CONSDATA* consdata;
   Real minactivity;
   Real maxactivity;
   Real lhs;
   Real rhs;
   Bool redundant;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   /*debugMessage("Prop method of linear constraints\n");*/

   *result = SCIP_DIDNOTFIND;
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert(cons != NULL);
      consdata = SCIPconsGetConsData(cons);
      assert(consdata != NULL);

      CHECK_OKAY( consdataGetActivityBounds(consdata, scip, &minactivity, &maxactivity) );
      CHECK_OKAY( consdataGetLhs(consdata, scip, &lhs) );
      CHECK_OKAY( consdataGetRhs(consdata, scip, &rhs) );

      if( SCIPisG(scip, minactivity, rhs) || SCIPisL(scip, maxactivity, lhs) )
      {
         debugMessage("linear constraint <%s> is infeasible: activitybounds=[%g,%g], sides=[%g,%g]\n",
            SCIPconsGetName(cons), minactivity, maxactivity, lhs, rhs);
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
      else
      {
         redundant = TRUE;
         if( SCIPisL(scip, minactivity, lhs) )
         {
            redundant = FALSE;
            todoMessage("linear constraint propagation for left hand side");
         }
         if( SCIPisG(scip, maxactivity, rhs) )
         {
            redundant = FALSE;
            todoMessage("linear constraint propagation for left hand side");
         }
         if( redundant )
         {
            debugMessage("linear constraint <%s> is redundant: activitybounds=[%g,%g], sides=[%g,%g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, lhs, rhs);
            todoMessage("disable redundant linear constraint");
         }
      }
   }

   return SCIP_OKAY;
}



/*
 * Callback methods of event handler
 */

static
DECL_EVENTEXEC(eventExecLinear)
{
   LINCONS* lincons;
   Real oldbound;
   Real newbound;
   int varpos;
   EVENTTYPE eventtype;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(event != NULL);

   debugMessage("Exec method of bound change event handler for linear constraints\n");

   lincons = eventdata->lincons;
   varpos = eventdata->varpos;
   assert(lincons != NULL);
   assert(0 <= varpos && varpos < lincons->nvars);

   CHECK_OKAY( SCIPeventGetType(event, &eventtype) );
   CHECK_OKAY( SCIPeventGetOldbound(event, &oldbound) );
   CHECK_OKAY( SCIPeventGetNewbound(event, &newbound) );
#ifndef NDEBUG
   {
      VAR* var;
      CHECK_OKAY( SCIPeventGetVar(event, &var) );
      assert(lincons->vars[varpos] == var);
   }
#endif

   debugMessage(" -> eventtype=0x%x, var=<%s>, oldbound=%g, newbound=%g => activity: [%g,%g]", 
      eventtype, SCIPvarGetName(lincons->vars[varpos]), oldbound, newbound, lincons->minactivity, lincons->maxactivity);

   if( (eventtype & SCIP_EVENTTYPE_LBCHANGED) != 0 )
      linconsUpdateChgLb(lincons, oldbound, newbound, lincons->vals[varpos]);
   else
   {
      assert((eventtype & SCIP_EVENTTYPE_UBCHANGED) != 0);
      linconsUpdateChgUb(lincons, oldbound, newbound, lincons->vals[varpos]);
   }

   debug(printf(" -> [%g,%g]\n", lincons->minactivity, lincons->maxactivity));

   return SCIP_OKAY;
}



/*
 * constraint specific interface methods
 */

RETCODE SCIPincludeConsHdlrLinear(      /**< creates the handler for linear constraints and includes it in SCIP */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create constraint handler data */
   ALLOC_OKAY( allocMemory(conshdlrdata) );

   /* create event handler for bound change events */
   CHECK_OKAY( SCIPeventhdlrCreate(&conshdlrdata->eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, 
                  NULL, NULL, NULL,
                  NULL, eventExecLinear,
                  NULL) );

   /* include constraint handler in SCIP */
   CHECK_OKAY( SCIPincludeConsHdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHCKPRIORITY, CONSHDLR_PROPFREQ,
                  CONSHDLR_NEEDSCONS,
                  consFreeLinear, NULL, NULL,
                  consDeleLinear, consTranLinear, 
                  consSepaLinear, consEnlpLinear, consEnpsLinear, consChckLinear, consPropLinear,
                  conshdlrdata) );

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
   if( SCIPstage(scip) == SCIP_STAGE_PROBLEM )
   {
      /* create constraint in original problem */
      CHECK_OKAY( linconsCreate(&consdata->data.lincons, scip, len, var, val, lhs, rhs, modifiable) );
   }
   else
   {
      CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* create constraint in transformed problem */
      CHECK_OKAY( linconsCreateTransformed(&consdata->data.lincons, scip, conshdlrdata->eventhdlr, len, var, val, lhs, rhs,
                     modifiable) );
   }

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
   Real             val                 /**< coefficient of constraint entry */
   )
{
   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   assert(cons != NULL);
   assert(scip != NULL);
   assert(var != NULL);

   /*debugMessage("adding coefficient %g * <%s> to linear constraint <%s>\n", val, var->name, SCIPconsGetName(cons));*/

   conshdlr = SCIPconsGetConsHdlr(cons);
   assert(conshdlr != NULL);

   if( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not linear");
      return SCIP_INVALIDDATA;
   }
   
   consdata = SCIPconsGetConsData(cons);
   if( consdata->islprow )
   {
      assert(consdata->data.row != NULL);
      CHECK_OKAY( SCIPaddVarToRow(scip, consdata->data.row, var, val) );
   }
   else
   {
      CONSHDLRDATA* conshdlrdata;

      assert(consdata->data.lincons != NULL);

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      CHECK_OKAY( linconsAddCoef(consdata->data.lincons, scip, conshdlrdata->eventhdlr, var, val) );
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
   
   CHECK_OKAY( consdataGetLhs(SCIPconsGetConsData(cons), scip, lhs) );

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
   
   CHECK_OKAY( consdataGetRhs(SCIPconsGetConsData(cons), scip, rhs) );

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
   
   consdata = SCIPconsGetConsData(cons);
   if( consdata->islprow )
   {
      todoMessage("introduce side change methods for rows");
      errorMessage("cannot change sides of linear constraint already stored as LP row (not implemented yet)");
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
   
   consdata = SCIPconsGetConsData(cons);
   if( consdata->islprow )
   {
      todoMessage("introduce side change methods for rows");
      errorMessage("cannot change sides of linear constraint already stored as LP row (not implemented yet)");
      return SCIP_INVALIDDATA;
   }
   else
   {
      assert(consdata->data.lincons != NULL);
      consdata->data.lincons->rhs = rhs;
   }

   return SCIP_OKAY;
}

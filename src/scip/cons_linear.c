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


/** externally stored linear constraint */
struct LinCons
{
   VAR**            var;                /**< variables of constraint entries */
   Real*            val;                /**< coefficients of constraint entries */
   Real             lhs;                /**< left hand side of row (for ranged rows) */
   Real             rhs;                /**< right hand side of row */
   int              size;               /**< size of the var- and val-arrays */
   int              len;                /**< number of nonzeros in constraint */
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
 * local methods
 */

static
RETCODE linconsCreate(                  /**< creates a linear constraint data object */
   LINCONS**        lincons,            /**< pointer to linear constraint data object */
   SCIP*            scip,               /**< SCIP data structure */
   int              len,                /**< number of nonzeros in the constraint */
   VAR**            var,                /**< array with variables of constraint entries */
   Real*            val,                /**< array with coefficients of constraint entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs                 /**< right hand side of row */
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

   if( len > 0 )
   {
      ALLOC_OKAY( SCIPduplicateBlockMemoryArray(scip, (*lincons)->var, var, len) );
      ALLOC_OKAY( SCIPduplicateBlockMemoryArray(scip, (*lincons)->val, val, len) );
   }
   else
   {
      (*lincons)->var = NULL;
      (*lincons)->val = NULL;
   }

   (*lincons)->lhs = lhs;
   (*lincons)->rhs = rhs;
   (*lincons)->size = len;
   (*lincons)->len = len;

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
   assert((*lincons)->size >= 0);

   SCIPfreeBlockMemoryArrayNull(scip, (*lincons)->var, (*lincons)->size);
   SCIPfreeBlockMemoryArrayNull(scip, (*lincons)->val, (*lincons)->size);
   SCIPfreeBlockMemory(scip, *lincons);
}

static
Real linconsGetSlack(                   /**< calculates the slack for the linear constraint */
   LINCONS*         lincons             /**< linear constraint data object */
   )
{
   Real slack;
   int v;

   assert(lincons != NULL);
   
   slack = lincons->rhs;
   for( v = 0; v < lincons->len; ++v )
      slack -= SCIPvarGetPrimsol(lincons->var[v]) * lincons->val[v];

   return slack;
}

static
Real linconsGetFeasibility(             /**< calculates the feasibility for the linear constraint */
   LINCONS*         lincons             /**< linear constraint data object */
   )
{
   Real slack;

   assert(lincons != NULL);

   slack = linconsGetSlack(lincons);

   return MIN(slack, lincons->rhs - lincons->lhs - slack);
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

   CHECK_OKAY( SCIPcreateRow(scip, row, name, 0, NULL, NULL, lincons->lhs, lincons->rhs) );
   
   for( v = 0; v < lincons->len; ++v )
   {
      CHECK_OKAY( SCIPaddVarToRow(scip, *row, lincons->var[v], lincons->val[v]) );
   }

   return SCIP_OKAY;
}

static
RETCODE applyConstraints(               /**< separates violated inequalities; called from Sepa and Enfo methods */
   CONSHDLR*        conshdlr,           /**< linear constraint handler */
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           conss,              /**< array of constraints to process */
   int              nconss              /**< number of constraints to process */
   )
{
   CONS* cons;
   CONSDATA* consdata;
   int c;
   Bool found;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(nconss == 0 || conss != NULL);
   assert(nconss >= 0);

   found = FALSE;

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

            feasibility = SCIProwGetFeasibility(row);
            debugMessage("row feasibility = %g\n", feasibility);
            if( !SCIPisFeasible(scip, feasibility) )
            {
               /* insert LP row as cut
                * we don't need to pool these rows, because they are already stored as linear constraints
                */
               CHECK_OKAY( SCIPaddCut(scip, row, -feasibility/SCIProwGetNorm(row)/(SCIProwGetNNonz(row)+1), FALSE) );
               found = TRUE;
            }
         }
      }
      else
      {
         LINCONS* lincons = consdata->data.lincons;
         Real feasibility;
         
         feasibility = linconsGetFeasibility(lincons);
         debugMessage("lincons feasibility = %g\n", feasibility);
         if( !SCIPisFeasible(scip, feasibility) )
         {
            ROW* row;
            
            /* convert lincons data into LP row */
            CHECK_OKAY( linconsToRow(lincons, SCIPconsGetName(cons), scip, &row) );
            assert(SCIPisEQ(scip, SCIProwGetFeasibility(row), feasibility));

            /* capture the row, because we need it as data storage */
            SCIPcaptureRow(scip, row);

            /* free the lincons data and convert consdata to point to the row */
            linconsFree(&consdata->data.lincons, scip);
            consdata->islprow = TRUE;
            consdata->data.row = row;
            
            /* insert LP row as cut
             * we don't need to pool these rows, because they are already stored as linear constraints
             */
            CHECK_OKAY( SCIPaddCut(scip, row, -feasibility/SCIProwGetNorm(row)/(SCIProwGetNNonz(row)+1), FALSE) );
            found = TRUE;
         }
      }
   }

   if( found )
      return SCIP_SEPARATED;
   else
      return SCIP_FEASIBLE;
}


/*
 * Callback methods
 */

DECL_CONSINIT(SCIPconsInitLinear)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   infoMessage(SCIPverbLevel(scip), SCIP_VERBLEVEL_FULL, "initialise linear constraint handler");

   return SCIP_OKAY;
}

DECL_CONSEXIT(SCIPconsExitLinear)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   infoMessage(SCIPverbLevel(scip), SCIP_VERBLEVEL_FULL, "exit linear constraint handler");

   return SCIP_OKAY;
}

DECL_CONSFREE(SCIPconsFreeLinear)
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

   CHECK_OKAY( linconsCreate(&(*targetdata)->data.lincons, scip, lincons->len, lincons->var, lincons->val,
                  lincons->lhs, lincons->rhs) );
   
   return SCIP_OKAY;
}

DECL_CONSSEPA(SCIPconsSepaLinear)
{
   RETCODE retcode;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   debugMessage("Sepa method of linear constraints\n");

   /* call enforcing method */
   CHECK_OKAY( retcode = applyConstraints(conshdlr, scip, conss, nconss) );

   /* if enforcing method returned SCIP_FEASIBLE, then the separation was not successful */
   if( retcode == SCIP_FEASIBLE )
      return SCIP_FAILURE;
   else
   {
      assert(retcode == SCIP_SEPARATED);
      return SCIP_SEPARATED;
   }
}

DECL_CONSENFO(SCIPconsEnfoLinear)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   debugMessage("Enfo method of linear constraints\n");

   /* call enforcing method */
   return applyConstraints(conshdlr, scip, conss, nconss);
}

DECL_CONSCHCK(SCIPconsChckLinear)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   todoMessage("Chck method of linear constraints");

   return SCIP_OKAY;
}

DECL_CONSPROP(SCIPconsPropLinear)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

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
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHCKPRIORITY,
                  SCIPconsInitLinear, SCIPconsExitLinear, SCIPconsFreeLinear,
                  SCIPconsTranLinear, SCIPconsSepaLinear, SCIPconsEnfoLinear, SCIPconsChckLinear, SCIPconsPropLinear,
                  NULL) );

   return SCIP_OKAY;
}

RETCODE SCIPcreateConsLinear(           /**< creates a linear constraint */
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              len,                /**< number of nonzeros in the constraint */
   VAR**            var,                /**< array with variables of constraint entries */
   Real*            val,                /**< array with coefficients of constraint entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
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
   consdata->islprow = FALSE;
   CHECK_OKAY( linconsCreate(&consdata->data.lincons, scip, len, var, val, lhs, rhs) );

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, model) );

   return SCIP_OKAY;
}

RETCODE SCIPcreateConsLPRow(            /**< creates a linear constraint from an LP row and captures the row */
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

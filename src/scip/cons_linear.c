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

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "cons_linear.h"


typedef struct LinCons LINCONS;         /**< externally stored linear constraint */


#define CONSHDLR_NAME "linear"
#define CONSHDLR_DESC "Linear constraints of the form  lhs <= a^T x <= rhs"

/** externally stored linear constraint */
struct LinCons
{
   char*            name;               /**< name of the constraint */
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
   const char*      name,               /**< name of constraint */
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
      sprintf(s, "constraint:%s (lhs=%f, rhs=%f)", name, lhs, rhs);
      errorMessage(s);
      return SCIP_INVALIDDATA;
   }

   ALLOC_OKAY( SCIPallocBlockMemory(scip, *lincons) );

   ALLOC_OKAY( SCIPduplicateBlockMemoryArray(scip, (*lincons)->name, name, strlen(name)+1) );

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
RETCODE linconsFree(                    /**< frees a linear constraint data object */
   LINCONS**        lincons,            /**< pointer to linear constraint data object */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(lincons != NULL);
   assert(*lincons != NULL);
   assert((*lincons)->size >= 0);

   SCIPfreeBlockMemoryArray(scip, (*lincons)->name, strlen((*lincons)->name)+1);
   SCIPfreeBlockMemoryArrayNull(scip, (*lincons)->var, (*lincons)->size);
   SCIPfreeBlockMemoryArrayNull(scip, (*lincons)->val, (*lincons)->size);
   SCIPfreeBlockMemory(scip, *lincons);

   return SCIP_OKAY;
}


/*
 * Callback methods
 */

DECL_CONSINIT(SCIPconsInitLinear)
{
   assert(self != NULL);
   assert(strcmp(SCIPgetConsHdlrName(self), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   infoMessage(SCIPverbLevel(scip), SCIP_VERBLEVEL_FULL, "initialise linear constraint handler");

   return SCIP_OKAY;
}

DECL_CONSEXIT(SCIPconsExitLinear)
{
   assert(self != NULL);
   assert(strcmp(SCIPgetConsHdlrName(self), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   infoMessage(SCIPverbLevel(scip), SCIP_VERBLEVEL_FULL, "exit linear constraint handler");

   return SCIP_OKAY;
}

DECL_CONSFREE(SCIPconsFreeLinear)
{
   assert(self != NULL);
   assert(strcmp(SCIPgetConsHdlrName(self), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   if( (*consdata)->islprow )
   {
      CHECK_OKAY( SCIPreleaseRow(scip, &(*consdata)->data.row) );
   }
   else
   {
      CHECK_OKAY( linconsFree(&(*consdata)->data.lincons, scip) );
   }
   SCIPfreeBlockMemory(scip, *consdata);

   return SCIP_OKAY;
}

DECL_CONSTRAN(SCIPconsTranLinear)
{
   LINCONS* lincons;

   assert(self != NULL);
   assert(strcmp(SCIPgetConsHdlrName(self), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(SCIPstage(scip) == SCIP_STAGE_SOLVING);
   assert(sourcedata != NULL);
   assert(!sourcedata->islprow);  /* in original problem, there cannot be LP rows */
   assert(sourcedata->data.lincons != NULL);
   assert(targetdata != NULL);

   lincons = sourcedata->data.lincons;

   ALLOC_OKAY( SCIPallocBlockMemory(scip, *targetdata) );
   (*targetdata)->islprow = FALSE;

   CHECK_OKAY( linconsCreate(&(*targetdata)->data.lincons, scip, lincons->name, lincons->len, lincons->var, lincons->val,
                  lincons->lhs, lincons->rhs) );
   
   return SCIP_OKAY;
}

DECL_CONSCHCK(SCIPconsChckLinear)
{
   assert(self != NULL);
   assert(strcmp(SCIPgetConsHdlrName(self), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(cons != NULL);
   assert(psol != NULL);

   todoMessage("Chck method of linear constraints");

   return SCIP_OKAY;
}

DECL_CONSPROP(SCIPconsPropLinear)
{
   assert(self != NULL);
   assert(strcmp(SCIPgetConsHdlrName(self), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(cons != NULL);

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
                  SCIPconsInitLinear, SCIPconsExitLinear, SCIPconsFreeLinear,
                  SCIPconsTranLinear, SCIPconsChckLinear, SCIPconsPropLinear,
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

   /* find the linear constraint handler, or create it if not existing */
   CHECK_OKAY( SCIPfindConsHdlr(scip, CONSHDLR_NAME, &conshdlr) );
   if( conshdlr == NULL )
   {
      errorMessage("Linear constraint handler not found");
      return SCIP_INVALIDCALL;
   }

   /* create the constraint specific data */
   ALLOC_OKAY( SCIPallocBlockMemory(scip, consdata) );
   consdata->islprow = FALSE;
   CHECK_OKAY( linconsCreate(&consdata->data.lincons, scip, name, len, var, val, lhs, rhs) );

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, conshdlr, consdata, model) );

   return SCIP_OKAY;
}

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
#define CONSHDLR_ENFOPRIORITY   -800000
#define CONSHDLR_CHECKPRIORITY  -800000
#define CONSHDLR_SEPAFREQ             4
#define CONSHDLR_PROPFREQ            -1
#define CONSHDLR_NEEDSCONS         TRUE

#define LINCONSUPGD_PRIORITY   +1000000

#define MINBRANCHWEIGHT               0.4  /**< minimum weight of both sets in set covering branching */
#define MAXBRANCHWEIGHT               0.9  /**< minimum weight of both sets in set covering branching */
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
   unsigned int     local:1;            /**< is constraint only valid locally? */
   unsigned int     modifiable:1;       /**< is data modifiable during node processing (subject to column generation)? */
};
typedef struct SetcoverCons SETCOVERCONS; /**< set covering constraint data */

/** constraint data for set covering constraints */
struct ConsData
{
   SETCOVERCONS* setcovercons;          /**< set covering constraint data */
   ROW*          row;                   /**< LP row, if constraint is already stored in LP row format */
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

   CHECK_OKAY( SCIPincIntarrayVal(scip, varuses, SCIPvarGetProbIndex(var), +1) );

   debugMessage("varuses of <%s>: %d\n", SCIPvarGetName(var), SCIPgetIntarrayVal(scip, varuses, SCIPvarGetProbIndex(var)));

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

   CHECK_OKAY( SCIPincIntarrayVal(scip, varuses, SCIPvarGetProbIndex(var), -1) );
   assert(SCIPgetIntarrayVal(scip, varuses, SCIPvarGetProbIndex(var)) >= 0);

   debugMessage("varuses of <%s>: %d\n", SCIPvarGetName(var), SCIPgetIntarrayVal(scip, varuses, SCIPvarGetProbIndex(var)));

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
   Bool             modifiable          /**< is constraint modifiable (subject to column generation)? */
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
   (*setcovercons)->local = local;
   (*setcovercons)->modifiable = modifiable;

   return SCIP_OKAY;
}   

/** creates a transformed set covering constraint data object */
static
RETCODE setcoverconsTransform(
   SCIP*            scip,               /**< SCIP data structure */
   SETCOVERCONS**   setcovercons,       /**< pointer to store the set covering constraint */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< variables of the constraint */
   Bool             local,              /**< is constraint only locally valid? */
   Bool             modifiable          /**< is constraint modifiable (subject to column generation)? */
   )
{
   int i;

   assert(setcovercons != NULL);
   assert(nvars == 0 || vars != NULL);

   CHECK_OKAY( setcoverconsCreate(scip, setcovercons, nvars, vars, local, modifiable) );

   /* use transformed variables in constraint instead original ones */
   for( i = 0; i < (*setcovercons)->nvars; ++i )
   {
      if( SCIPvarGetStatus((*setcovercons)->vars[i]) == SCIP_VARSTATUS_ORIGINAL )
      {
         (*setcovercons)->vars[i] = SCIPvarGetTransformed((*setcovercons)->vars[i]);
         assert((*setcovercons)->vars[i] != NULL);
      }
      assert(SCIPvarGetStatus((*setcovercons)->vars[i]) != SCIP_VARSTATUS_ORIGINAL);
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
                  setcovercons->local, setcovercons->modifiable) );
   
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

   CHECK_OKAY( setcoverconsTransform(scip, &targetdata->setcovercons, setcovercons->nvars, setcovercons->vars,
                  setcovercons->local, setcovercons->modifiable) );
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
   Bool*            separated           /**< pointer to store, iff a cut was found */
   )
{
   CONSDATA* consdata;
   ROW* row;

   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(separated != NULL);

   *separated = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   row = consdata->row;
   if( row != NULL )
   {
      /* ignore constraints that are in the LP */
      if( !SCIProwIsInLP(row) )
      {         
         Real feasibility;
            
         CHECK_OKAY( SCIPgetRowFeasibility(scip, row, &feasibility) );
         /*debugMessage("  row feasibility = %g\n", feasibility);*/
         if( !SCIPisFeasible(scip, feasibility) )
         {
            /* insert LP row as cut */
            CHECK_OKAY( SCIPaddCut(scip, row, -feasibility/SCIProwGetNorm(row)/(SCIProwGetNNonz(row)+1)) );
            CHECK_OKAY( SCIPresetConsAge(scip, cons) );
            *separated = TRUE;
         }
         else
         {
            CHECK_OKAY( SCIPincConsAge(scip, cons) );
         }
      }
   }
   else
   {
      SETCOVERCONS* setcovercons;
      VAR** vars;
      Real solval;
      Real sum;
      int nvars;
      int v;

      setcovercons = consdata->setcovercons;
      assert(setcovercons != NULL);
         
      /* calculate the constraint's activity */
      vars = setcovercons->vars;
      nvars = setcovercons->nvars;
      sum = 0.0;
      for( v = 0; v < nvars && sum < 1.0; ++v )
      {
         assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY);
         CHECK_OKAY( SCIPgetSolVal(scip, NULL, vars[v], &solval) );
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
            CHECK_OKAY( SCIPgetRowActivity(scip, consdata->row, &rowactivity) );
            assert(SCIPisSumEQ(scip, rowactivity, sum));
         }
#endif
            
         /* insert LP row as cut */
         CHECK_OKAY( SCIPaddCut(scip, consdata->row, (1.0-sum)/(nvars+1)) );
         CHECK_OKAY( SCIPresetConsAge(scip, cons) );
         *separated = TRUE;
      }
      else
      {
         assert(1 <= v && v <= nvars);
            
         if( SCIPisEQ(scip, SCIPvarGetLb(vars[v-1]), 1.0) )
         {
            /* constraint is always feasible due to the bounds of vars[v-1] -> disable constraint */
            CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
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
   Bool separated;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   debugMessage("separating %d/%d set covering constraints\n", nusefulconss, nconss);

   *result = SCIP_DIDNOTFIND;

   /* step 1: check all useful set covering constraints for feasibility */
   for( c = 0; c < nusefulconss; ++c )
   {
      CHECK_OKAY( separate(scip, conss[c], &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   }

   /* step 2: combine set covering constraints to get more cuts */
   todoMessage("further cuts of set covering constraints");

   return SCIP_OKAY;
}

static
DECL_CONSENFOLP(consEnfolpSetcover)
{
   Bool separated;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   debugMessage("enforcing %d set covering constraints\n", nconss);

   *result = SCIP_FEASIBLE;

   /* step 1: check all useful set covering constraints for feasibility */
   for( c = 0; c < nusefulconss; ++c )
   {
      CHECK_OKAY( separate(scip, conss[c], &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   }

   /* step 2: if solution is not integral, choose a variable set to branch on */
   todoMessage("set covering branching");

   /* step 3: check all obsolete set covering constraints for feasibility */
   for( c = nusefulconss; c < nconss; ++c )
   {
      CHECK_OKAY( separate(scip, conss[c], &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   }


#if 0
   CONSDATA* consdata;
   SETCOVERCONS* setcovercons;
   ROW* row;
   VAR** vars;
   Real solval;
   Real sum;
   Real branchconssum;
   Bool foundcuts;
   int branchcons;
   int branchconsnvars;
   int nvars;
   int v;
   int foundvar;

   /* check all set covering constraints for feasibility and integrality
    * separating violated constraints is preferred, but if no violation could be found,
    * set covering branching is applied on a constraint with smallest activity and
    * smallest number of variables
    */
   foundcuts = FALSE;
   branchcons = -1;
   branchconssum = SCIPinfinity(scip);
   branchconsnvars = INT_MAX;
   nbranchconss = 0;
   for( c = 0; c < nusefulconss || (c < nconss && !foundcuts && branchcons == -1); ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      row = consdata->row;
      if( row != NULL )
      {
         Real activity;
         
         CHECK_OKAY( SCIPgetRowActivity(scip, row, &activity) );
         /*debugMessage("  row activity = %g\n", activity);*/
         if( SCIPisFeasLT(scip, activity, 1.0) )
         {
            /* insert LP row as cut */
            CHECK_OKAY( SCIPaddCut(scip, row, -feasibility/SCIProwGetNorm(row)/(SCIProwGetNNonz(row)+1)) );
            CHECK_OKAY( SCIPresetConsAge(scip, conss[c]) );
            *result = SCIP_SEPARATED;
            foundcuts = TRUE;
         }
         else
         {
            CHECK_OKAY( SCIPincConsAge(scip, conss[c]) );
         }
      }
      else
      {
         setcovercons = consdata->setcovercons;
         assert(setcovercons != NULL);
         
         /* search all variables of the constraint until one is found with solution value 1.0 */
         vars = setcovercons->vars;
         nvars = setcovercons->nvars;
         foundvar = -1;
         sum = 0.0;
         for( v = 0; v < nvars && foundvar == -1; ++v )
         {
            assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY);
            CHECK_OKAY( SCIPgetSolVal(scip, NULL, vars[v], &solval) );
            assert(SCIPisGE(scip, solval, 0.0) || SCIPisLE(scip, solval, 1.0));
            if( SCIPisEQ(scip, solval, 1.0) )
               foundvar = v;
            sum += solval;
         }
         if( foundvar == -1 )
         {
            todoMessage("vvvvvv remove this");
            assert(SCIPisFeasLT(scip, sum, 1.0)); /* ?????????????????? */

            /* no variable with value 1.0 found, and sum of all values is small enough to branch on the constraint */
            CHECK_OKAY( SCIPresetConsAge(scip, conss[c]) );
            
            if( SCIPisFeasLT(scip, sum, 1.0) )
            {
               /* convert set covering constraint data into LP row */
               CHECK_OKAY( setcoverconsToRow(scip, setcovercons, SCIPconsGetName(conss[c]), &consdata->row) );
               assert(consdata->row != NULL);
               assert(!SCIProwIsInLP(consdata->row));
#ifndef NDEBUG
               {
                  Real rowactivity;
                  CHECK_OKAY( SCIPgetRowActivity(scip, consdata->row, &rowactivity) );
                  assert(SCIPisSumEQ(scip, rowactivity, sum));
               }
#endif

               /* insert LP row as cut */
               CHECK_OKAY( SCIPaddCut(scip, consdata->row, (1.0-sum)/(nvars+1)) );
               *result = SCIP_SEPARATED;
               foundcuts = TRUE;
               
               debugMessage("added set covering cut:");
               debug( SCIPprintRow(scip, consdata->row, NULL) );
            }
            else if( !setcovercons->modifiable && !foundcuts )
            {
               if( SCIPisLT(scip, sum, branchconssum) || (SCIPisLE(scip, sum, branchconssum) && nvars < branchconsnvars) )
               {
                  branchcons = c;
                  branchconssum = sum;
                  branchconsnvars = nvars;
               }
               nbranchconss++;
            }
         }
         else
         {
            assert(0 <= foundvar && foundvar < nvars);

            if( SCIPisEQ(scip, SCIPvarGetLb(vars[foundvar]), 1.0) )
            {
               /* constraint is always feasible due to the bounds of vars[v] -> disable constraint */
               CHECK_OKAY( SCIPdisableConsLocal(scip, conss[c]) );
            }
            else
            {
               /* constraint was feasible -> increase age */
               CHECK_OKAY( SCIPincConsAge(scip, conss[c]) );
            }
         }
      }
   }

   if( !foundcuts && nbranchconss > 0 )
   {
      CONS* cons;
      VAR** avars;
      VAR* bvar;
      Real asum;
      int navars;
   
      abort(); /* ????????????????? */

      /* no constraint was violated -> branch on the set covering constraint best suited for branching
       *  - partition variables into two sets A and B
       *  - create two sons:  (i) x(A) >= 1,   (ii) x(A) == 0
       */
      assert(*result == SCIP_FEASIBLE);
      assert(0 <= branchcons && branchcons < nconss);
      
      *result = SCIP_INFEASIBLE;
      
      cons = conss[branchcons];
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);
      setcovercons = consdata->setcovercons;
      assert(setcovercons != NULL);
      vars = setcovercons->vars;
      nvars = setcovercons->nvars;
      
      /* get temporary memory */
      CHECK_OKAY( SCIPcaptureBufferArray(scip, &avars, nvars) );
      asum = 0.0;
      navars = 0;
      
      /* partition the variables into two sets A and B, try to make x(A), x(B) >= MINBRANCHWEIGHT */
      for( v = 0; v < nvars && asum < MINBRANCHWEIGHT; ++v )
      {
         CHECK_OKAY( SCIPgetSolVal(scip, NULL, vars[v], &solval) );
         if( SCIPisPositive(scip, solval) )
         {
            avars[navars] = vars[v];
            asum += solval;
            navars++;
         }
      }
      
      if( asum < MAXBRANCHWEIGHT )
      {
         NODE* node;
         CONS* newcons;
         char name[255];
         
         assert(navars < nvars);
         assert(v < nvars);
         bvar = vars[v];
         assert(bvar != NULL);
         
         /* create left child, add constraint x(A) >= 1, disable old constraint */
         sprintf(name, "%s_A", SCIPconsGetName(cons));
         CHECK_OKAY( SCIPcreateChild(scip, &node) );
         if( navars == 1 )
         {
            debugMessage("fixing variable <%s> to 1.0 in left child node\n", SCIPvarGetName(avars[0]));
            CHECK_OKAY( SCIPchgVarLbNode(scip, node, avars[0], 1.0) );
         }
         else
         {
            CHECK_OKAY( SCIPcreateConsSetcover(scip, &newcons, name, navars, avars,
                           SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                           TRUE, FALSE) );
            CHECK_OKAY( SCIPaddConsNode(scip, node, newcons) );
            CHECK_OKAY( SCIPreleaseCons(scip, &newcons) );
         }
         CHECK_OKAY( SCIPdisableConsNode(scip, node, cons) );
         
         /* create right child, fix x_i = 0 for all i \in A */
         CHECK_OKAY( SCIPcreateChild(scip, &node) );
         for( v = 0; v < navars; ++v )
         {
            CHECK_OKAY( SCIPchgVarUbNode(scip, node, avars[v], 0.0) );
         }
         if( !setcovercons->modifiable && (nvars - navars == 1) )
         {
            /* only one variable in set B -> fix it to 1.0 and disable constraint */
            debugMessage("fixing variable <%s> to 1.0 in right child node\n", SCIPvarGetName(bvar));
            CHECK_OKAY( SCIPchgVarLbNode(scip, node, bvar, 1.0) );
            CHECK_OKAY( SCIPdisableConsNode(scip, node, cons) );
         }
         
         *result = SCIP_BRANCHED;
         
#ifdef DEBUG
         debugMessage("set covering branching: navars=%d/%d, weight(A)=%g/%g, A={", navars, nvars, asum, sum);
         for( v = 0; v < navars; ++v )
         {
            CHECK_OKAY( SCIPgetSolVal(scip, NULL, avars[v], &solval) );
            printf(" %s[%g]", SCIPvarGetName(avars[v]), solval);
         }
         printf(" }\n");
#endif
      }
      
      /* free temporary memory */
      CHECK_OKAY( SCIPreleaseBufferArray(scip, &avars) );
   }

   assert(!foundcuts || (*result == SCIP_SEPARATED));
   assert((foundcuts || nbranchconss > 0) ^ (*result == SCIP_FEASIBLE));
#endif

   return SCIP_OKAY;
}

static
DECL_CONSENFOPS(consEnfopsSetcover)
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

   abort(); /* ???????????????? */

   *result = SCIP_FEASIBLE;

   /* check all set covering constraints for feasibility */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      setcovercons = consdata->setcovercons;
      assert(setcovercons != NULL);

      /* search all variables of the constraint until one is found with solution value 1.0 */
      vars = setcovercons->vars;
      nvars = setcovercons->nvars;
      found = FALSE;
      for( v = 0; v < nvars && !found; ++v )
      {
         assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY);
         CHECK_OKAY( SCIPgetSolVal(scip, NULL, vars[v], &solval) );
         assert(SCIPisEQ(scip, solval, 0.0) || SCIPisEQ(scip, solval, 1.0));
         found = SCIPisEQ(scip, solval, 1.0);
      }
      if( !found )
      {
         VAR** branchvars;
         int nbranchvars;

         CHECK_OKAY( SCIPresetConsAge(scip, conss[c]) );

         /* constraint is violated: branch on pseudo solution */

         /* get temporary memory */
         CHECK_OKAY( SCIPcaptureBufferArray(scip, &branchvars, nvars) );

         /* find the non-fixed variables in the constraint */
         nbranchvars = 0;
         for( v = 0; v < nvars; ++v )
         {
            if( SCIPisEQ(scip, SCIPvarGetUb(vars[v]), 1.0) )
            {
               branchvars[nbranchvars] = vars[v];
               nbranchvars++;
            }
         }

         if( nbranchvars == 0 )
         {
            /* all variables in the constraint are fixed to 0.0
             * -> if constraint is modifiable, we have to price in additional variables
             * -> if constraint is not modifiable (subject to column generation), the node can be cut off
             */
            if( setcovercons->modifiable )
            {
               errorMessage("modifiable infeasible set covering constraint detected; pricing on pseudo solution not implemented yet");
               abort();
            }
            else
               *result = SCIP_CUTOFF;

            debugMessage("set cover constraint <%s> cuts off node\n", SCIPconsGetName(conss[c]));
         }
         else
         {
            NODE* node;

            /* for each non-fixed variable v, create a child node x_0 = ... = x_v-1 = 0, x_v = 1
             * if the constraint is modifiable, we need the additional child x_0 = ... = x_n = 0 which forces
             * the pricing of new variables to fulfil the constraint
             */
            for( v = 0; v < nbranchvars; ++v )
            {
               int u;

               /* create child with x_0 = ... = x_v-1 = 0, x_v = 1 */
               CHECK_OKAY( SCIPcreateChild(scip, &node) );
               for( u = 0; u < v; ++u )
               {
                  CHECK_OKAY( SCIPchgVarUbNode(scip, node, branchvars[u], 0.0) );
               }
               CHECK_OKAY( SCIPchgVarLbNode(scip, node, branchvars[v], 1.0) );
            }
            if( setcovercons->modifiable )
            {
               /* create child with x_0 = ... = x_v = 0 */
               CHECK_OKAY( SCIPcreateChild(scip, &node) );
               for( v = 0; v < nbranchvars; ++v )
               {
                  CHECK_OKAY( SCIPchgVarUbNode(scip, node, branchvars[v], 0.0) );
               }
            }
            *result = SCIP_BRANCHED;

#if DEBUG
            {
               int nchildren;
               CHECK_OKAY( SCIPgetChildren(scip, NULL, &nchildren) );
               debugMessage("branched on set cover constraint in pseudo solution: %d children\n", nchildren);
            }
#endif
         }
         
         /* free temporary memory */
         CHECK_OKAY( SCIPreleaseBufferArray(scip, &branchvars) );

         return SCIP_OKAY;
      }
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
            CHECK_OKAY( SCIPgetSolVal(scip, sol, vars[v], &solval) );
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
                     local, FALSE) );
      *upgraded = TRUE;
   }
   else
      *upgraded = FALSE;

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
   Bool             modifiable          /**< is row modifiable during node processing (subject to column generation)? */
   )
{
   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   assert(scip != NULL);

   /* find the set covering constraint handler */
   CHECK_OKAY( SCIPfindConsHdlr(scip, CONSHDLR_NAME, &conshdlr) );
   if( conshdlr == NULL )
   {
      errorMessage("set covering constraint handler not found");
      return SCIP_INVALIDCALL;
   }

   /* create the constraint specific data */
   CHECK_OKAY( SCIPallocBlockMemory(scip, &consdata) );
   CHECK_OKAY( setcoverconsCreate(scip, &consdata->setcovercons, nvars, vars, local, modifiable) );
   consdata->row = NULL;

   /* create constraint (propagation is never used for set covering constraints) */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, separate, enforce, check, FALSE) );

   return SCIP_OKAY;
}


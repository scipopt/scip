/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2008 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_countsols.c,v 1.14 2008/04/17 17:49:04 bzfpfets Exp $"

/**@file   cons_countsols.c
 * @brief  constraint handler for counting feasible solutions
 * @author Stefan Heinz
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_and.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/cons_countsols.h"
#include "scip/dialog_default.h"
#include "scip/pub_cons.h"
#include "scip/pub_disp.h"

#ifdef WITH_GMP
#include <gmp.h>
typedef mpz_t                Int;
#else
typedef SCIP_Longint         Int;
#endif 

/* constraint handler properties */
#define CONSHDLR_NAME          "countsols"
#define CONSHDLR_DESC          "constraint to count feasible solutions"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -9999999 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -9999999 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS         0 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

/* default parameter settings */
#define DEFAULT_SPARSETEST         TRUE /**< sparse test on or off */
#define DEFAULT_DISCARDSOLS        TRUE /**< is it allowed to discard solutions */
#define DEFAULT_ACTIVE            FALSE /**< is the constraint handler active */

/* default column settings */
#define DISP_SOLS_NAME             "sols"
#define DISP_SOLS_DESC             "number of detected feasible solutions"
#define DISP_SOLS_HEADER           " sols "
#define DISP_SOLS_WIDTH            6
#define DISP_SOLS_PRIORITY         110000
#define DISP_SOLS_POSITION         100000
#define DISP_SOLS_STRIPLINE        TRUE

#define DISP_CUTS_NAME             "feasST"
#define DISP_CUTS_DESC             "number of detected non trivial feasible subtrees"
#define DISP_CUTS_HEADER           "feasST"
#define DISP_CUTS_WIDTH            6
#define DISP_CUTS_PRIORITY         110000
#define DISP_CUTS_POSITION         110000
#define DISP_CUTS_STRIPLINE        TRUE

/** creates and adds a constraint which cuts off the solution from the feasibility
 *  region 
 * 
 *  input:
 *  - scip            : SCIP main data structure
 *  - sol             : solution to cut off 
 *  - conshdlrdata    : constraint handler data 
*/
#define CUTOFF_CONSTRAINT(x) SCIP_RETCODE x (SCIP* scip, SCIP_SOL* sol, SCIP_CONSHDLRDATA* conshdlrdata)


/** constraint handler data */
struct SCIP_ConshdlrData
{
   int                   feasST;             /**< number of non trivial feasible subtrees */
   int                   nDiscardSols ;      /**< number of discard solutions */
   int                   nNonSparseSols;     /**< number of non sparse solutions */
   Int                   nsols;              /**< number of solutions */
   CUTOFF_CONSTRAINT((*cutoffSolution));     /**< method for cutting of a solution */
   SCIP_Bool             active;             /**< constraint handler active */
   SCIP_Bool             discardsols;        /**< allow to discard solutions */
   SCIP_Bool             sparsetest;         /**< allow to check for sparse solutions */

   int                   nvars;              /**< number of variables in problem */
   SCIP_VAR**            vars;               /**< array containing a copy of all variables before presolving */
};


/* 
 * Local methods for handling the <Int> data structure 
 */

/** allocates memory for the value pointer */
static
void allocInt(
   Int*          value                       /**< pointer to the value to allocate memory */
   )
{
#ifdef WITH_GMP
   mpz_init(*value);
#endif
}


/** sets the value pointer to the new value */
static
void setInt(
   Int*          value,                       /**< pointer to the value to initialize */
   int           newvalue                     /**< new value */
   )
{
#ifdef WITH_GMP
   mpz_set_si(*value, newvalue);
#else
   (*value) = newvalue;
#endif
}


/** free memoy */
static
void freeInt(
   Int*          value                      /**< pointer to the value to free */
   )
{
#ifdef WITH_GMP
   mpz_clear(*value);
#endif
}


/** adds one to the given value */
static
void addOne(
   Int*          value                      /**< pointer to the value to increase */
   )
{
#ifdef WITH_GMP
   mpz_add_ui(*value, *value, 1);
#else
   (*value)++;
#endif
}


/** adds the summand to the given value */
static
void addInt(
   Int*          value,                     /**< pointer to the value to increase */
   Int*          summand                    /**< summand to add on */
   )
{
#ifdef WITH_GMP
   mpz_add(*value, *value, *summand);
#else
   (*value) += (*summand);
#endif
}


/** multiplies the factor to the given vakue */
static
void multInt(
   Int*          value,                     /**< pointer to the value to increase */
   SCIP_Longint  factor                     /**< factor to multiply with */
   )
{
#ifdef WITH_GMP
   mpz_mul_ui (*value, *value, factor);
#else
   (*value) *= factor; 
#endif
}


/* method for creating a string out of an Int which is a mpz_t or SCIP_Longint */
static 
void toString(
   Int      value,                          /**< number */
   char**   buffer                          /**< pointer to buffer for storing the string */
   )
{
#ifdef WITH_GMP
   mpz_get_str(*buffer, 10, value);
#else
   sprintf (*buffer, "%"SCIP_LONGINT_FORMAT"", value);
#endif
}


/*
 * Local methods
 */


/** returns whether a given integer variable is unfixed in the local domain */
static
SCIP_Bool varIsUnfixedLocal(
   SCIP_VAR*             var                 /**< integer variable */
   )
{
   assert( var != NULL );
   assert( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS );
   assert( SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var) >= 0.0 );
   
   return ( SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var) > 0.5 );
}


/** creates the constraint handler data */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                      scip,            /**< SCIP data structure */
   SCIP_CONSHDLRDATA**        conshdlrdata     /**< pointer to store constraint handler data */
   )
{

   SCIP_CALL( SCIPallocMemory(scip, conshdlrdata) );
   
   (*conshdlrdata)->feasST = 0;
   (*conshdlrdata)->nDiscardSols = 0;
   (*conshdlrdata)->nNonSparseSols = 0;
   
   allocInt(&(*conshdlrdata)->nsols);
   
   (*conshdlrdata)->cutoffSolution = NULL;

   (*conshdlrdata)->nvars = 0;
   (*conshdlrdata)->vars = NULL;

   return SCIP_OKAY;
}


/** check if the current parameter setting is correct for a save counting process */
static
SCIP_RETCODE checkParameters(
   SCIP*                      scip             /**< SCIP data structure */
   )
{
   SCIP_HEUR** heuristics;
   int nheuristics;

   int h;
   int value;
   
   SCIP_Bool valid;

   assert( scip != NULL );
   
   valid = TRUE;
   
   if( SCIPgetStage(scip) < SCIP_STAGE_PRESOLVED )
   {
      /* check if dual methods are turnred off */
      SCIPgetIntParam(scip, "presolving/dualfix/maxrounds", &value);
      if( value != 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, 
            "The presolver <dualfix> is not turned off! This might cause a wrong counting process.");
      }
   }

   /* check if all heuristics are turned off */
   heuristics = SCIPgetHeurs(scip);
   nheuristics = SCIPgetNHeurs(scip);

   for( h = 0; h < nheuristics; ++h )
   {
      if( SCIPheurGetFreq(heuristics[h]) != -1 )
      {
         valid = FALSE;
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, 
            "The heuristic <%s> is not turned off! This might cause a wrong counting process.\n", 
            SCIPheurGetName(heuristics[h]));
      }
   }
  
   /* check if restart is turned off */
   SCIPgetIntParam(scip,  "presolving/maxrestarts", &value);
   if( value != 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, 
         "The parameter <presolving/maxrestarts> is not 0 (currently %d)! This might cause a wrong counting process.",
         value);
   }
   
   if( !valid && SCIPgetVerbLevel(scip) != SCIP_VERBLEVEL_FULL )
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, 
         "The current parameter setting might cause a wrong counting process.");
   
   return SCIP_OKAY;
}


/** creates and adds a constraints which cuts off the current solution from the
 *  feasibility region in the case there are only binary variables */
static 
CUTOFF_CONSTRAINT(addBinaryCons)
{
   int v;

   SCIP_VAR** consvars;
   SCIP_VAR** vars;
   int nvars;

   SCIP_Real value;
   SCIP_VAR* var;
   SCIP_VAR* negvar;
   SCIP_CONS* cons;
   
   assert( scip != NULL );
   assert( sol != NULL );
   assert( conshdlrdata != NULL );
    
   SCIP_CALL( SCIPgetPseudoBranchCands(scip, &vars, &nvars, NULL) );
   assert( nvars > 0 );
   
   /* allocate buffer memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nvars) );
   
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
    
      assert( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
      assert( varIsUnfixedLocal(var) );

      value = SCIPgetSolVal(scip, sol, var);
      assert( SCIPisIntegral(scip, value) );

      if (value > 0.5)
      {
         negvar = var;
         SCIP_CALL( SCIPgetNegatedVar(scip, negvar, &vars[nvars]));
      }
      else
         consvars[nvars] = var;
   }
    
   /* create constraint */
   SCIP_CALL( SCIPcreateConsSetcover(scip, &cons, "Setcovering created by countsols", nvars, consvars,
         FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE));
   
   /* add and release constraint */
   SCIP_CALL( SCIPaddConsLocal(scip, cons, NULL) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   
   /* free buffer array */
   SCIPfreeBufferArray(scip, &consvars);
  
   return SCIP_OKAY;
}


/** creates and adds a bound disjunction constraints which cuts off the current solution
 *  from the feasibility region; if only binary variables are involved, then a set
 *  covering constraint is created which is a special case of a bound disjunction
 *  constraint */
static 
CUTOFF_CONSTRAINT(addIntegerCons)
{
   int v;

   SCIP_VAR** consvars;
   SCIP_VAR** vars;
   SCIP_Real* bounds;
   SCIP_BOUNDTYPE* boundtypes;
   int nvars;
   int nbinvars = 0;
   int nconsvars;
   SCIP_VAR* var;
   SCIP_Real value;
   SCIP_Longint lb,ub,valueInt;
   
   SCIP_CONS* cons;

   assert( scip != NULL );
   assert( sol != NULL );
   assert( conshdlrdata != NULL );
  
   SCIP_CALL( SCIPgetPseudoBranchCands(scip, &vars, &nvars, NULL) );
   nconsvars = nvars * 2;
   assert( nvars > 0 );

   /* allocate buffer memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nconsvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bounds, nconsvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtypes, nconsvars) );
   
   nconsvars = 0;

   for( v = nvars - 1; v >= 0; --v )
   {
      var = vars[v];
    
      assert( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS );
      assert( varIsUnfixedLocal(var) );

      if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
      {
         ++nbinvars;
         value = SCIPgetSolVal(scip, sol, var);
         assert( SCIPisIntegral(scip, value) );
         
         if (value < 0.5)
         {
            boundtypes[nconsvars] = SCIP_BOUNDTYPE_LOWER;
            bounds[nconsvars] = 1;
         }
         else 
         {
            boundtypes[nconsvars] = SCIP_BOUNDTYPE_UPPER;
            bounds[nconsvars] = 0;
         }
      }
      else
      {
         assert( SCIPisIntegral(scip, SCIPvarGetLbLocal(var)) );
         assert( SCIPisIntegral(scip, SCIPvarGetUbLocal(var)) );
         assert( SCIPisIntegral(scip, SCIPgetSolVal(scip, sol, var)) );
            
         lb = (SCIP_Longint) SCIPfeasCeil(scip, SCIPvarGetLbLocal(var));
         ub = (SCIP_Longint) SCIPfeasCeil(scip, SCIPvarGetUbLocal(var));
         valueInt = (SCIP_Longint) SCIPfeasCeil(scip, SCIPgetSolVal(scip, sol, var));
         
         if (valueInt == lb)
         {
            boundtypes[nconsvars] = SCIP_BOUNDTYPE_LOWER;
            bounds[nconsvars] = lb + 1;
         }
         else if (valueInt == ub)
         {
            boundtypes[nconsvars] = SCIP_BOUNDTYPE_UPPER;
            bounds[nconsvars] = ub - 1;
         }
         else
         {
            boundtypes[nconsvars] = SCIP_BOUNDTYPE_LOWER;
            bounds[nconsvars] = valueInt + 1;
            consvars[nconsvars] = var;
            ++nconsvars;
            boundtypes[nconsvars] = SCIP_BOUNDTYPE_UPPER;
            bounds[nconsvars] = valueInt - 1;
         }
      }
      
      consvars[nconsvars] = var;
      ++nconsvars;
   }
   
   /* check if only binary variables appear in the constraint; if this is the case we
    * create a set covering constraint instead of a bound disjunction constraint */
   if (nvars == nbinvars )
   {
      for (v = nbinvars - 1; v >= 0; --v)
      {
         /* in the case the bound is zero we have use the negated variable */
         if( bounds[v] == 0)
         {
            SCIP_CALL( SCIPgetNegatedVar(scip, consvars[v], &consvars[v]));
         }
      }
    
      SCIP_CALL( SCIPcreateConsSetcover(scip, &cons, "Setcovering created by countsols", nbinvars, consvars,
            FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE));
   }
   else
   {
      SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &cons, "Bounddisjunction created by countsols", 
            nconsvars, consvars, boundtypes, bounds,
            FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
   }
  
   /* add and release constraint locally */
   SCIP_CALL( SCIPaddConsLocal(scip, cons, NULL) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
  
   /* free buffer memory */
   SCIPfreeBufferArray(scip, &consvars);
   SCIPfreeBufferArray(scip, &bounds);
   SCIPfreeBufferArray(scip, &boundtypes);
  
   return SCIP_OKAY;
}


/** counts the number of solutions represented by sol */
static 
SCIP_RETCODE countSparsesol(
   SCIP*                      scip,             /**< SCIP data structure */
   SCIP_SOL*                  sol,              /**< solution */
   SCIP_Bool                  feasible,         /**< bool if solution is feasible */
   SCIP_CONSHDLRDATA*         conshdlrdata,     /**< constraint handler data */
   SCIP_RESULT*               result            /**< pointer to store the result of the checking process */
   )
{
   assert( scip != NULL );
   assert( sol != NULL );
   assert( conshdlrdata != NULL );
   assert( result != NULL );
   
   if( feasible )
   {
      int v;
      
      Int newsols;
      
      SCIP_VAR** vars;
      int nvars;

      SCIP_VAR* var;
      SCIP_Real lb;
      SCIP_Real ub;
      
      SCIPdebugMessage("counts number of solutions represented through the given one\n");
      
      /**@note aggregations and multi aggregations: we do not have to care about these things
       *       since we count solution from the transformed problem and therefore, SCIP does
       *       it for us */
      
      assert( SCIPgetNPseudoBranchCands(scip) != 0 );
      
      allocInt(&newsols);
      setInt(&newsols, 1);
      
      if( SCIPgetNBinVars(scip) == SCIPgetNVars(scip) )
      {
         SCIP_Longint factor = 1;

         assert( SCIPgetNPseudoBranchCands(scip) < 64 );
         multInt(&newsols, (factor <<= SCIPgetNPseudoBranchCands(scip)) );
      }
      else
      {
         SCIP_CALL( SCIPgetPseudoBranchCands(scip, &vars, &nvars, NULL) );
         for( v = 0; v < nvars; ++v )
         {
            var = vars[v];
            lb = SCIPvarGetLbLocal(var);
            ub = SCIPvarGetUbLocal(var);
            
            SCIPdebugMessage("variable <%s> Local Bounds are [%g,%g]\n", SCIPvarGetName(var), lb, ub);
            
            assert( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS );
            assert( SCIPisIntegral(scip, lb) );
            assert( SCIPisIntegral(scip, ub) );
            assert( SCIPisIntegral(scip, ub - lb) );
            assert( SCIPisLT(scip, lb, ub) );
            
            /* the number of integers laying in the interval [lb,ub] is
             *  (ub - lb + 1); to make everything integral we add another
             *  0.5 and cut the fractional part off */
            multInt(&newsols, (SCIP_Longint)(ub - lb + 1.5) );
         }
         
#ifdef DEBUG
         char buffer[SCIP_MAXSTRLEN];
         
         toString(newsols, &buffer);
         SCIPdebugMessage("add %s solutions\n", buffer );
#endif
      }
      
      *result = SCIP_CUTOFF;
      conshdlrdata->feasST++;
      
      addInt(&conshdlrdata->nsols, &newsols);
      freeInt(&newsols);
   }
   else if(!conshdlrdata->discardsols)
   {
      conshdlrdata->cutoffSolution(scip, sol, conshdlrdata);
      addOne(&conshdlrdata->nsols);
      conshdlrdata->nNonSparseSols++;
   }
   else
      conshdlrdata->nDiscardSols++;
   
   return SCIP_OKAY;
}


/** checks if the new solution is feasible for the logicor constraints */
static 
SCIP_Bool checkLogicor(
   SCIP*                      scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,         /**< constraint handler */
   int                        nconss            /**< number of enabled constraints */
   
   ) 
{
   /**@note the logicor constraints are not fully propagated; therefore, we have to check
    *       them by hand if they are satisfied or not; if a constraint is satisfied we
    *       delete it locally from the branch and bound tree. */
   
   int c, v;
   SCIP_CONS** conss;
   int nvars;
   SCIP_VAR** vars;
   SCIP_Bool satisfied;

   SCIPdebugMessage("check logicor %d contraints\n", nconss);
  
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr),"logicor") == 0 );
   assert( nconss == SCIPconshdlrGetNEnabledConss(conshdlr) );
   
   conss = SCIPconshdlrGetConss(conshdlr);
   assert( conss != NULL );
   
   c = SCIPconshdlrGetNActiveConss(conshdlr) - 1;

   for( ; c >= 0 && nconss > 0; --c )
   {
      SCIPdebugMessage("logicor contraint %d\n", c);
    

      if( !SCIPconsIsEnabled(conss[c]) )
         continue;
    
      nconss--;

      nvars = SCIPgetNVarsLogicor(scip, conss[c]);
      vars = SCIPgetVarsLogicor(scip, conss[c]);
    
      /* calculate the constraint's activity */
      satisfied = FALSE;
      for( v = 0; v < nvars && !satisfied; ++v )
      {
         assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY);

         if( !varIsUnfixedLocal(vars[v] ) ) 
            satisfied = SCIPvarGetLbLocal(vars[v]) > 0.5;
      }

      if( !satisfied )
      {
         SCIPdebugMessage("constraint <%s> cannot be disabled\n", SCIPconsGetName(conss[c]));
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) ) );
         return FALSE;
      }
      
      /* delete constraint from the problem locally since it is satisfied */
      SCIPdelConsLocal(scip, conss[c]);
   }
   return TRUE;
}


/** checks if the new solution is feasible for the knapsack constraints */
static 
SCIP_Bool checkKnapsack(
   SCIP*                      scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,         /**< constraint handler */
   int                        nconss            /**< number of enabled constraints */
   
   ) 
{
   /**@note the knapsack constraints are not fully propagated; therefore, we have to check
    *       them by hand if they are satisfied or not; if a constraint is satisfied we
    *       delete it locally from the branch and bound tree. */
   
   int c, v;

   SCIP_CONS** conss;

   int nvars;
   SCIP_VAR** vars;
   SCIP_Longint* weights;
   SCIP_Longint capacity;
   SCIP_Real capa;

   SCIPdebugMessage("check knapsack %d contraints\n", nconss);
  
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr),"knapsack") == 0 );
   assert( nconss == SCIPconshdlrGetNEnabledConss(conshdlr) );
  
   conss = SCIPconshdlrGetConss(conshdlr);
   assert( conss != NULL );
  

   c = SCIPconshdlrGetNActiveConss(conshdlr) - 1;

   for( ; c >= 0 && nconss > 0; --c )
   {
      SCIPdebugMessage("knapsack contraint %d\n", c);
    
      if( !SCIPconsIsEnabled(conss[c]) )
         continue;
    
      nconss--;
    
      nvars = SCIPgetNVarsKnapsack(scip, conss[c]);
      vars = SCIPgetVarsKnapsack(scip, conss[c]);
      capacity = SCIPgetCapacityKnapsack(scip, conss[c]);
      weights = SCIPgetWeightsKnapsack(scip,conss[c]);
    
      SCIPdebugMessage("knapsack capacity = %"SCIP_LONGINT_FORMAT"\n", capacity);
    
      capa = capacity + 0.1;
    
      for( v = nvars - 1; v >= 0 && capa >= 0 ; --v )
      {
         SCIPdebug( SCIP_CALL( SCIPprintVar( scip, vars[v], NULL) ) );
         SCIPdebugMessage("weight = %"SCIP_LONGINT_FORMAT" :\n", weights[v]);
         assert( SCIPvarIsIntegral(vars[v]) );
         
         /* the weights should be greater or equal to zero */
         assert( weights[v] > -0.5 );
         assert( weights[v] >= 0);
      
         if ( !varIsUnfixedLocal(vars[v]) ) 
         {
            /* variables is fixed locally; therefore, subtract fixed variable value multiplied by
             * the weight; */
            capa -= weights[v] * SCIPvarGetLbLocal(vars[v]); 
         }
         else if (weights[v] > 0.5) 
         {
            /*  variable is unfixed and weight is greater than 0; therefore, subtract upper bound
             *  value multiplied by the weight */
            capa -= weights[v] * SCIPvarGetUbLocal(vars[v]); 
         }
      }
      
      if( SCIPisFeasLT(scip, capa, 0.0) )
      {
         SCIPdebugMessage("constraint %s cannot be disabled\n", SCIPconsGetName(conss[c]));
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) ) );
         return FALSE;
      }
    
      /* delete constraint from the problem locally since it is satisfied */
      SCIPdelConsLocal(scip, conss[c]);
   }
   return TRUE;
}


/** checks if the new solution is feasible for the bounddisjunction constraints */
static 
SCIP_Bool checkBounddisjunction(
   SCIP*                      scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,         /**< constraint handler */
   int                        nconss            /**< number of enabled constraints */
   ) 
{
   /**@note the bounddisjunction constraints are not fully propagated; therefore, we have to check
    *       them by hand if they are satisfied or not; if a constraint is satisfied we
    *       delete it locally from the branch and bound tree */
  
   int c, v;

   SCIP_CONS** conss;

   int nvars;
   SCIP_VAR** vars;
   SCIP_BOUNDTYPE* boundtypes;
   SCIP_Real* bounds;
   SCIP_Bool satisfied;
  
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr),"bounddisjunction") == 0 );
   assert( nconss == SCIPconshdlrGetNEnabledConss(conshdlr) );
  
   conss = SCIPconshdlrGetConss(conshdlr);
   assert( conss != NULL );
  

   c = SCIPconshdlrGetNActiveConss(conshdlr) - 1;
   for( ; c >= 0 && nconss > 0; --c )
   {
      if( !SCIPconsIsEnabled(conss[c]) )
         continue;
    
      nconss--;
      satisfied = FALSE;
    
      nvars = SCIPgetNVarsBounddisjunction(scip, conss[c]);
      vars = SCIPgetVarsBounddisjunction(scip, conss[c]);
    
      boundtypes = SCIPgetBoundtypesBounddisjunction(scip, conss[c]);
      bounds = SCIPgetBoundsBounddisjunction(scip, conss[c]);
    
      for( v = nvars-1; v >= 0 && !satisfied; --v )
      {
         SCIPdebug( SCIPprintVar(scip, vars[v], NULL) );
         assert( SCIPvarGetType(vars[v]) != SCIP_VARTYPE_CONTINUOUS );
      
         /* variable should be in right bounds to delete constraint */
         if (boundtypes[v] == SCIP_BOUNDTYPE_LOWER )
            satisfied = SCIPisFeasGE(scip, SCIPvarGetLbLocal(vars[v]), bounds[v]);
         else
         {
            assert( boundtypes[v] == SCIP_BOUNDTYPE_UPPER );
            satisfied = SCIPisFeasLE(scip, SCIPvarGetUbLocal(vars[v]), bounds[v]);
         }
      }
    
      if (!satisfied)
      {
         SCIPdebugMessage("constraint %s cannot be disabled\n", SCIPconsGetName(conss[c]));
         SCIPdebug(SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) ) );
         return FALSE;
      }
    
      /* delete constraint from the problem locally since it is satisfied */
      SCIPdelConsLocal(scip, conss[c]);
   }
   return TRUE;
}


/** checks if the new solution is feasible for the varbound constraints */
static
SCIP_Bool checkVarbound(
   SCIP*                      scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,         /**< constraint handler */
   int                        nconss            /**< number of enabled constraints */
   
   ) 
{
   /**@note the varbound constraints are not fully propagated; therefore, we have to check
    *       them by hand if they are satisfied or not; if a constraint is satisfied we
    *       delete it locally from the branch and bound tree. */
   
   int c;

   SCIP_CONS** conss;

   SCIP_VAR* var_x;
   SCIP_VAR* var_y;
   SCIP_Real lhs, rhs, coef;

   SCIPdebugMessage("check varbound %d contraints\n", nconss);
  
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr),"varbound") == 0 );
   assert( nconss == SCIPconshdlrGetNEnabledConss(conshdlr) );
  
   conss = SCIPconshdlrGetConss(conshdlr);
   assert( conss != NULL );
  

   c = SCIPconshdlrGetNActiveConss(conshdlr) - 1;

   for( ; c >= 0 && nconss > 0; --c )
   {
      SCIPdebugMessage("varbound contraint %d\n", c);
    
      if( !SCIPconsIsEnabled(conss[c]) )
         continue;
    
      nconss--;
    
      var_x = SCIPgetVarVarbound(scip, conss[c]);
      var_y = SCIPgetVbdvarVarbound(scip, conss[c]);

      assert (SCIPvarGetType(var_y) != SCIP_VARTYPE_CONTINUOUS);
      
      coef = SCIPgetVbdcoefVarbound(scip, conss[c]);
      lhs = SCIPgetLhsVarbound(scip, conss[c]);
      rhs = SCIPgetRhsVarbound(scip, conss[c]);
      
      /* variables y is fixed locally; therefore, subtract fixed variable value multiplied by
       * the coefficient; */
      if (SCIPisGT(scip, SCIPvarGetUbLocal(var_x), rhs - SCIPvarGetUbLocal(var_y) * coef ) 
         || !SCIPisGE(scip, SCIPvarGetLbLocal(var_x), lhs - SCIPvarGetLbLocal(var_y) * coef ))
      {
         SCIPdebugMessage("constraint %s cannot be disabled\n", SCIPconsGetName(conss[c]));
         SCIPdebug(SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) ) );
         SCIPdebugMessage("%s\t lb: %lf\t ub: %lf\n",SCIPvarGetName(var_x), SCIPvarGetLbLocal(var_x), SCIPvarGetUbLocal(var_x));
         SCIPdebugMessage("%s\t lb: %lf\t ub: %lf\n",SCIPvarGetName(var_y), SCIPvarGetLbLocal(var_y), SCIPvarGetUbLocal(var_y));
         
        return FALSE;
      }
      
      /* delete constraint from the problem locally since it is satisfied */
      SCIPdelConsLocal(scip, conss[c]);
   }
   return TRUE;
}


/** check if the current node initializes a non trivial feasible subtree */
static 
SCIP_RETCODE checkFeasSubtree(
   SCIP* scip,                         /**< SCIP main data structure */
   SCIP_SOL* sol,                      /**< solution to check */
   SCIP_Bool* feasible                 /**< pointer to store the result of the check */
   )
{
   int h;

   SCIP_CONSHDLR** conshdlrs;
   int nconshdlrs;

   SCIP_CONSHDLR* conshdlr;
   int nconss;

   SCIPdebugMessage("check if the sparse solution is feasible\n");

   assert( scip != NULL );
   assert( sol != NULL );
   assert( feasible != NULL );
   
   assert( SCIPgetNPseudoBranchCands(scip) != 0 );
   
   *feasible = FALSE;
  
   nconshdlrs = SCIPgetNConshdlrs(scip) - 1;
   conshdlrs = SCIPgetConshdlrs(scip);
   assert (conshdlrs != NULL);

   /* check each constraint handler if there are constraints which are not enabled */
   for (h = nconshdlrs ;  h >= 0 ; --h )
   {
      conshdlr = conshdlrs[h];
      assert( conshdlr != NULL );
      
      nconss = SCIPconshdlrGetNEnabledConss(conshdlr);
      
      /* skip this constraints handler */
      if( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 )
         continue;

      if( nconss > 0 )
      {
         SCIPdebugMessage("constraint handler %s has %d active constraint(s)\n",
            SCIPconshdlrGetName(conshdlr), nconss );

         if (strcmp(SCIPconshdlrGetName(conshdlr), "logicor") == 0)
         {
            if (checkLogicor(scip, conshdlr, nconss) == FALSE )
            {
               SCIPdebugMessage("a <logicor> constraint cannot be disabled\n");
               return SCIP_OKAY;
            }
         }
         else if (strcmp(SCIPconshdlrGetName(conshdlr), "knapsack") == 0)
         {
            if (checkKnapsack(scip, conshdlr, nconss) == FALSE )
            {
               SCIPdebugMessage("a <knapsack> constraint cannot be disabled\n");
               return SCIP_OKAY;
            }
         }
         else if (strcmp(SCIPconshdlrGetName(conshdlr), "bounddisjunction") == 0)
         {
            if (checkBounddisjunction(scip, conshdlr, nconss) == FALSE )
            {
               SCIPdebugMessage("a <bounddisjunction> constraint cannot be disabled\n");
               return SCIP_OKAY;
            }
         }
         else if (strcmp(SCIPconshdlrGetName(conshdlr), "varbound") == 0)
         {
            if (checkVarbound(scip, conshdlr, nconss) == FALSE )
            {
               SCIPdebugMessage("a <varbound> constraint cannot be disabled\n");
               return SCIP_OKAY;
            }
         }
         else
         {
            SCIPdebugMessage("sparse solution is infeasible since the following constraint (and maybe more) is(/are) enabled\n");
            SCIPdebug( SCIP_CALL( SCIPprintCons(scip, SCIPconshdlrGetConss(conshdlr)[0], NULL) ) );
            return SCIP_OKAY;
         }
      }
   }
   
   *feasible = TRUE;
   SCIPdebugMessage("sparse solution is feasible\n");
   
   return SCIP_OKAY;
}


/** check the given solution */
static 
SCIP_RETCODE checkSolution(
   SCIP*                    scip,            /**< SCIP data structure */
   SCIP_SOL*                sol,             /**< solution to add */
   SCIP_CONSHDLRDATA*       conshdlrdata,    /**< constraint handler data */
   SCIP_RESULT*             result           /**< pointer to store the result of the checking process */
   )
{
   SCIP_Bool feasible;

   SCIPdebugMessage("start to add sparse solution\n");

   assert( scip != NULL );
   assert( sol != NULL );
   assert( conshdlrdata != NULL );
   assert( result != NULL );
   
   /* the solution should not be found through a heuristic since in this case the
    * informations of SCIP are not valid for this solution */
   /**@todo it might be not necessary to check this assert since we can check in generale
       all solutions of feasibility independently of the origin; however, the locally fixed
       technique does only work if the solution comes from the branch and bound tree; in
       case the solution comes from a heuristic we should try to sequentially fix the
       variables in the branch and bound tree and check after every fixing if all
       constraints are disabled; at the point where all constraints are disabled the
       unfixed variables are stars; */
   assert( SCIPgetNVars(scip) != 0);
   assert( SCIPsolGetHeur(sol) == NULL);

#ifdef SCIP_DEBUG
   SCIP_VAR* var;
   SCIP_VAR** vars;
   int v;
   int nvars;

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      SCIPdebugMessage("variables <%s> Local Bounds are [%g,%g] Global Bounds are [%g,%g]\n",
         SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
   }
#endif
   
   /* setting result to infeasible since we reject any solution; however,
    * if the solution passes the sparse test the result is set to
    * SCIP_CUTOFF which cuts off the subtree initialized through the current
    * node */
   *result = SCIP_INFEASIBLE;
   
   /* check if solution is completely fixed */
   if( SCIPgetNPseudoBranchCands(scip) == 0 )
   {
      addOne(&conshdlrdata->nsols);
      conshdlrdata->nNonSparseSols++;
   }
   else if( conshdlrdata->sparsetest )
   {
      SCIP_CALL( checkFeasSubtree(scip, sol, &feasible) ) ;
      SCIP_CALL( countSparsesol(scip, sol, feasible, conshdlrdata, result) );
   }
   
   assert( *result == SCIP_INFEASIBLE || *result == SCIP_CUTOFF );
   SCIPdebugMessage("result is %s\n", *result == SCIP_INFEASIBLE ? "SCIP_INFEASIBLE" : "SCIP_CUTOFF" );

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static 
SCIP_DECL_CONSFREE(consFreeCountsols)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(&conshdlrdata != NULL);

   /* free conshdlrdata */
   freeInt(&conshdlrdata->nsols);
   
   SCIPfreeMemory(scip, &conshdlrdata);
   SCIPconshdlrSetData(conshdlr, NULL);
   
   return SCIP_OKAY;
}

/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitCountsols)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL );

   if( conshdlrdata->active )
   {
      int v;
      conshdlrdata->nvars = SCIPgetNVars(scip);

      SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->vars, conshdlrdata->nvars) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &conshdlrdata->vars, SCIPgetVars(scip), conshdlrdata->nvars) );

      /* capture all variables */
      for( v = 0; v < conshdlrdata->nvars; ++v )
      {
         SCIP_CALL( SCIPcaptureVar(scip, conshdlrdata->vars[v]) );
      }
   }

   return SCIP_OKAY;
}

/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitCountsols)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL );
   
   if( conshdlrdata->vars != NULL )
   {
      int v;

      /* release all variables */
      for( v = 0; v < conshdlrdata->nvars; ++v )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &conshdlrdata->vars[v]) );
      }
            
      SCIPfreeMemoryArrayNull(scip, &conshdlrdata->vars);
      conshdlrdata->nvars = 0;
   }      
   
   return SCIP_OKAY;
}

/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#define consInitpreCountsols NULL

/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreCountsols NULL


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The constraint handler may use this call to initialize its branch and bound specific data.
 */
static 
SCIP_DECL_CONSINITSOL(consInitsolCountsols)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   
   assert( SCIPgetStage(scip) == SCIP_STAGE_SOLVING );

   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL );

   /* check if the problem is binary */
   if( SCIPgetNBinVars(scip) == SCIPgetNVars(scip) )
      conshdlrdata->cutoffSolution = &addBinaryCons;
   else
      conshdlrdata->cutoffSolution = &addIntegerCons;
   
   return SCIP_OKAY;
}

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#define consExitsolCountsols NULL

/** frees specific constraint data */
#define consDeleteCountsols NULL

/** transforms constraint data into data belonging to the transformed problem */
#define consTransCountsols NULL

/** LP initialization method of constraint handler */
#define consInitlpCountsols NULL

/** separation method of constraint handler for LP solutions */
#define consSepalpCountsols NULL

/** separation method of constraint handler for arbitrary primal solutions */
#define consSepasolCountsols NULL

/** constraint enforcing method of constraint handler for LP solutions */
static 
SCIP_DECL_CONSENFOLP(consEnfolpCountsols)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIPdebugMessage("method SCIP_DECL_CONSENFOLP(consEnfolpCountsols)\n");
   
   assert( scip != NULL );
   assert( conshdlr != NULL );   
   assert( nconss == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if( conshdlrdata->active )
   {
      if( !solinfeasible )
      {
         SCIP_SOL* sol;
         
         SCIP_CALL( SCIPcreateLPSol(scip, &sol, NULL ) );

         SCIP_CALL( checkSolution(scip, sol, conshdlrdata, result) );
         SCIP_CALL( SCIPfreeSol(scip, &sol) );
      }
      else
         *result = SCIP_INFEASIBLE;
   }
   else
      *result = SCIP_FEASIBLE;
   
   assert( !conshdlrdata->active || *result == SCIP_INFEASIBLE || *result == SCIP_CUTOFF );
   
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsCountsols)
{ /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   
   SCIPdebugMessage("method SCIP_DECL_CONSENFOPS(consEnfopsCountsols)\n");

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( nconss == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

      
   if( conshdlrdata->active )
   {
      if( !solinfeasible )
      {
         SCIP_SOL* sol;
         
         SCIP_CALL( SCIPcreatePseudoSol(scip, &sol, NULL ) );
         
         SCIP_CALL(checkSolution(scip, sol, conshdlrdata, result) );
         SCIP_CALL( SCIPfreeSol(scip, &sol) );
      }
      else
         *result = SCIP_INFEASIBLE;
   }
   else
      *result = SCIP_FEASIBLE;
   
   assert( !conshdlrdata->active || *result == SCIP_INFEASIBLE || *result == SCIP_CUTOFF );
   
   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static 
SCIP_DECL_CONSCHECK(consCheckCountsols)
{  /*lint --e{715}*/  
   /**@todo solutions which come from scip_ckeck should be ignored since it is not clear who
    *       generated these solution; later we should analyze this problem */
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIPdebugMessage("method SCIP_DECL_CONSCHECK(consCheckCountsols)\n");

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if( conshdlrdata->active )
   {
      SCIPwarningMessage("a solution comes in over <SCIP_DECL_CONSCHECK(consCheckCountsols)>; right now these solutions are ignored\n");
      *result = SCIP_INFEASIBLE;
   }
   else
      *result = SCIP_FEASIBLE;
      
   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
#define consPropCountsols NULL

/** presolving method of constraint handler */
#define consPresolCountsols NULL

/** propagation conflict resolving method of constraint handler */
#define consRespropCountsols NULL


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockCountsols)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL );
   
   if( conshdlrdata->vars != NULL )
   {
      int v;
      
      for( v = 0; v < conshdlrdata->nvars; ++v )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, conshdlrdata->vars[v], +1, +1) );
      }
   }
   
   return SCIP_OKAY;
}

/** constraint activation notification method of constraint handler */
#define consActiveCountsols NULL

/** constraint deactivation notification method of constraint handler */
#define consDeactiveCountsols NULL

/** constraint enabling notification method of constraint handler */
#define consEnableCountsols NULL

/** constraint disabling notification method of constraint handler */
#define consDisableCountsols NULL

/** constraint display method of constraint handler */
#define consPrintCountsols NULL


/*
 * Interface methods
 */

/** dialog execution method for the count command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecCount)
{  /*lint --e{715}*/
   SCIP_RETCODE retcode;
   SCIP_HEUR** heuristics;
   int* heuristicfreqs;    
   int nheuristics;        
   SCIP_Bool heuristicsoff;
   SCIP_Bool active;

   int maxrestarts;        
   int maxroundsdualfix;
   
   int h;
   char parametername[256];

   SCIP_Bool valid;
   SCIP_Longint nsols;
   int displayprimalbound;
   int displaygap;
   int displaysols;
   int displayfeasST;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );
   SCIPdialogMessage(scip, NULL, "\n");

   switch( SCIPgetStage(scip) )
   {
   case SCIP_STAGE_INIT:
      SCIPdialogMessage(scip, NULL, "no problem exists\n");
      break;
      
   case SCIP_STAGE_PROBLEM:
      /* activate constraint handler cons_countsols */
      SCIP_CALL( SCIPgetBoolParam(scip, "constraints/"CONSHDLR_NAME"/active", &active) );
      if( !active )
      {
         SCIP_CALL( SCIPsetBoolParam(scip, "constraints/"CONSHDLR_NAME"/active", TRUE) );
      }
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
      /* turn off dual methods */
      if( SCIPfindPresol(scip, "dualfix") != NULL )
      {
         SCIP_CALL( SCIPgetIntParam(scip, "presolving/dualfix/maxrounds", &maxroundsdualfix) );
         
         if( maxroundsdualfix != 0 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, 
               "set parameter <presolving/dualfix/maxrounds> to 0\n");
            SCIP_CALL( SCIPsetIntParam(scip, "presolving/dualfix/maxrounds", 0) );
         }
      }

      /* presolve problem */
      SCIP_CALL( SCIPpresolve(scip) );
      
      /* reset parametername setting for dual methods */
      if( SCIPfindPresol(scip, "dualfix") != NULL )
      {
         if( maxroundsdualfix != 0 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, 
               "reset parameter <presolving/dualfix/maxrounds> to %d\n", maxroundsdualfix);
            SCIP_CALL( SCIPsetIntParam(scip,  "presolving/dualfix/maxrounds", maxroundsdualfix) );
         }
      }
   case SCIP_STAGE_PRESOLVED:
      /* reset activity status of constraint handler cons_countsols */
      if( !active )
      {
         SCIP_CALL( SCIPsetBoolParam(scip, "constraints/"CONSHDLR_NAME"/active", FALSE) );
      }
   case SCIP_STAGE_SOLVING:
      /* check if the problem contains continuous variables */
      if( SCIPgetNContVars(scip) != 0 )
      {   
         SCIPdialogMessage(scip, NULL, "invalid problem for counting; there are continuous variables (after presolving)\n"); 
         break;
      }
      
      /* set parameter setting for a correct count and store the old setting */
      nheuristics = SCIPgetNHeurs(scip);
         
      SCIP_CALL( SCIPallocMemoryArray(scip, &heuristicfreqs, nheuristics) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &heuristics, SCIPgetHeurs(scip), nheuristics) );
      heuristicsoff = TRUE;

      for( h = 0; h < nheuristics; ++h )
      {
         sprintf(parametername, "heuristics/%s/freq", SCIPheurGetName(heuristics[h]));
         heuristicfreqs[h] = SCIPheurGetFreq(heuristics[h]);
         
         if( heuristicfreqs[h] != -1 )
         {
            heuristicsoff = FALSE;
            SCIP_CALL( SCIPsetIntParam(scip, parametername, -1) );
         }
         assert( SCIPheurGetFreq(heuristics[h]) == -1 );
      }
      
      if( !heuristicsoff )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, 
            "set all heuristic calling frequencies to -1\n");
      }

      /* turn off restarts */
      SCIP_CALL( SCIPgetIntParam(scip, "presolving/maxrestarts", &maxrestarts) );

      if( maxrestarts != 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, 
            "set parameter <presolving/maxrestarts> to 0\n");
         SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrestarts", 0) );
      }

      /* turn off primal bound and gap column */
      SCIP_CALL( SCIPgetIntParam(scip, "display/primalbound/active", &displayprimalbound) );
      if( displayprimalbound != 0 )
         SCIP_CALL( SCIPsetIntParam(scip, "display/primalbound/active", 0) );
      SCIP_CALL( SCIPgetIntParam(scip, "display/gap/active", &displaygap) );
      if( displaygap != 0 )
         SCIP_CALL( SCIPsetIntParam(scip, "display/gap/active", 0) );
      
      /* turn on sols and feasST column */
      SCIP_CALL( SCIPgetIntParam(scip, "display/sols/active", &displaysols) );
      if( displayprimalbound != 2 )
         SCIP_CALL( SCIPsetIntParam(scip, "display/sols/active", 2) );
      SCIP_CALL( SCIPgetIntParam(scip, "display/feasST/active", &displayfeasST) );
      if( displayprimalbound != 2 )
         SCIP_CALL( SCIPsetIntParam(scip, "display/feasST/active", 2) );
      
      /* find the countsols constraint handler */
      assert( SCIPfindConshdlr(scip, CONSHDLR_NAME) != NULL );
      
      retcode =  SCIPcount(scip);
      
      valid = FALSE;
      nsols = SCIPgetNCountedSols(scip, &valid);
      
      if( valid )
         SCIPdialogMessage(scip, NULL, "Feasible Solutions : %"SCIP_LONGINT_FORMAT"", nsols);
      else
      {
         char* buffer;
         int buffersize = SCIP_MAXSTRLEN;
         int requiredsize;

         SCIP_CALL( SCIPallocBufferArray(scip, &buffer, buffersize) );
         SCIPgetNCountedSolsstr(scip, &buffer, buffersize, &requiredsize);
         
         if( requiredsize > buffersize )
         {
            SCIP_CALL( SCIPreallocBufferArray(scip, &buffer, requiredsize) );
            SCIPgetNCountedSolsstr(scip, &buffer, buffersize, &requiredsize);
            
         }

         assert( buffersize >= requiredsize );
         SCIPdialogMessage(scip, NULL, "Feasible Solutions : %s", buffer);
         
         SCIPfreeBufferArray(scip, &buffer);
      }

      SCIPdialogMessage(scip, NULL, " (%d non-trivial feasible subtrees)\n", SCIPgetNCountedFeasSubtrees(scip));

      *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

      if( !heuristicsoff )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, 
            "reset heuristic frequencies\n");

         /* reset parameters setting for heuristics */
         for( h = 0; h < nheuristics; ++h )
         {
            sprintf(parametername, "heuristics/%s/freq", SCIPheurGetName(heuristics[h]));
            SCIP_CALL( SCIPsetIntParam(scip, parametername, heuristicfreqs[h]) );
         }
      }
      
      /* reset parametername setting for restarts */
      if( maxrestarts != 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, 
            "reset parameter <presolving/maxrestarts> to %d\n", maxrestarts);
         SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrestarts", maxrestarts) );
      }
      
      /* reset display columns */
      if( displayprimalbound != 0 )
         SCIP_CALL( SCIPsetIntParam(scip, "display/primalbound/active", displayprimalbound) );
      if( displaygap != 0 )
         SCIP_CALL( SCIPsetIntParam(scip, "display/gap/active", displaygap) );
      
      /* reset sols and feasST column */
      if( displaysols != 2 )
         SCIP_CALL( SCIPsetIntParam(scip, "display/sols/active", displaysols) );
      if( displayfeasST != 2 )
         SCIP_CALL( SCIPsetIntParam(scip, "display/feasST/active", displayfeasST) );

      /* free array */
      SCIPfreeMemoryArrayNull(scip, &heuristicfreqs);
      SCIPfreeMemoryArrayNull(scip, &heuristics);
      
      /* evaluate retcode */
      SCIP_CALL( retcode );
      break;
      
   case SCIP_STAGE_SOLVED:
      SCIPdialogMessage(scip, NULL, "problem is already solved\n");
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
   default:
      SCIPerrorMessage("invalid SCIP stage\n");
      return SCIP_INVALIDCALL;
   }
   
   SCIPdialogMessage(scip, NULL, "\n");
   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}


/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(dispOutputSols)
{  /*lint --e{715}*/
   SCIP_CONSHDLR* conshdlr;
   SCIP_Longint sols;
   SCIP_Bool valid;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_SOLS_NAME) == 0);
   assert(scip != NULL);
   
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );
   assert( SCIPconshdlrGetNConss(conshdlr) == 0 );
   
   sols = SCIPgetNCountedSols(scip, &valid);
   
   if( !valid )
   {
      SCIPmessageFPrintInfo(file, "ToMany");
   }
   else
   {
      SCIPdispLongint(file, sols, DISP_SOLS_WIDTH);
   }
   
   return SCIP_OKAY;
}


/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(dispOutputFeasSubtrees)
{ /*lint --e{715}*/
   SCIP_CONSHDLR* conshdlr;
   
   assert(disp != NULL);
   assert(scip != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_CUTS_NAME) == 0);
   
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );
   assert( SCIPconshdlrGetNConss(conshdlr) == 0 );
   
   SCIPdispInt(file, SCIPgetNCountedFeasSubtrees(scip), DISP_CUTS_WIDTH);
   
   return SCIP_OKAY;
}


/** creates the handler for countsols constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrCountsols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   /* create countsol constraint handler data */
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_DIALOG* root;
   SCIP_DIALOG* dialog;
   
   /* create constraint handler specific data here */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata) );
   
   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, 
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         consFreeCountsols, consInitCountsols, consExitCountsols, 
         consInitpreCountsols, consExitpreCountsols, consInitsolCountsols, consExitsolCountsols,
         consDeleteCountsols, consTransCountsols, consInitlpCountsols,
         consSepalpCountsols, consSepasolCountsols, consEnfolpCountsols, consEnfopsCountsols, consCheckCountsols, 
         consPropCountsols, consPresolCountsols, consRespropCountsols, consLockCountsols,
         consActiveCountsols, consDeactiveCountsols, 
         consEnableCountsols, consDisableCountsols,
         consPrintCountsols,
         conshdlrdata) );

   /* add countsols constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, 
         "constraints/"CONSHDLR_NAME"/sparsetest", 
         "should the sparse solution test be turned on?",
         &conshdlrdata->sparsetest, FALSE, DEFAULT_SPARSETEST, NULL, NULL));
   SCIP_CALL( SCIPaddBoolParam(scip, 
         "constraints/"CONSHDLR_NAME"/discardsols", 
         "is it allowed to discard solutions?",
         &conshdlrdata->discardsols, FALSE, DEFAULT_DISCARDSOLS, NULL, NULL));
   SCIP_CALL( SCIPaddBoolParam(scip, 
         "constraints/"CONSHDLR_NAME"/active", 
         "is the constraint handler active?",
         &conshdlrdata->active, FALSE, DEFAULT_ACTIVE, NULL, NULL));
   
   /* add dialog entry for counting */
   root = SCIPgetRootDialog(scip);
   if( root == NULL )
   {
      SCIP_CALL( SCIPcreateRootDialog(scip, &root) );
   }
   assert( root != NULL );
   
   if( !SCIPdialogHasEntry(root, "count") )
   {
      SCIP_CALL( SCIPcreateDialog(scip, &dialog, SCIPdialogExecCount, NULL, NULL,
            "count", "count number of feasible solutions", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* include display column */
   SCIP_CALL( SCIPincludeDisp(scip, DISP_SOLS_NAME, DISP_SOLS_DESC, DISP_SOLS_HEADER, SCIP_DISPSTATUS_OFF, 
         NULL, NULL, NULL, NULL, NULL, dispOutputSols, 
         NULL, DISP_SOLS_WIDTH, DISP_SOLS_PRIORITY, DISP_SOLS_POSITION, DISP_SOLS_STRIPLINE) );
   SCIP_CALL( SCIPincludeDisp(scip, DISP_CUTS_NAME, DISP_CUTS_DESC, DISP_CUTS_HEADER, SCIP_DISPSTATUS_OFF, 
         NULL, NULL, NULL, NULL, NULL, dispOutputFeasSubtrees, 
         NULL, DISP_CUTS_WIDTH, DISP_CUTS_PRIORITY, DISP_CUTS_POSITION, DISP_CUTS_STRIPLINE) );
   
   return SCIP_OKAY;
}


/* execute counting */
extern
SCIP_RETCODE SCIPcount(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Bool active;

   /* activate constraint handler cons_countsols */
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/"CONSHDLR_NAME"/active", &active) );
   if( !active )
   {
      SCIP_CALL( SCIPsetBoolParam(scip, "constraints/"CONSHDLR_NAME"/active", TRUE) );
   }

   /* check if the parameter setting allows a valid counting process */
   SCIP_CALL( checkParameters(scip) );
   
   /* start the solving process */
   SCIP_CALL( SCIPsolve(scip) );
   
   /* reset activity status of constraint handler cons_countsols */
   if( !active )
   {
      SCIP_CALL( SCIPsetBoolParam(scip, "constraints/"CONSHDLR_NAME"/active", FALSE) );
   }
   
   return SCIP_OKAY;
}


/** returns number of feasible solutions found as SCIP_Longint; if the number does not fit into 
 *  a SCIP_Longint the valid flag is set to FALSE */
SCIP_Longint SCIPgetNCountedSols(
   SCIP*                 scip,              /**< SCIP data structure */
   SCIP_Bool*            valid              /**< pointer to store if the return value is valid */             
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   
   /* find the countsols constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

#ifdef WITH_GMP
   *valid = FALSE;
   if( 0 != mpz_fits_slong_p(conshdlrdata->nsols) )
      (*valid) = TRUE;
   
   return mpz_get_si(conshdlrdata->nsols);
#else
   *valid = TRUE;
   return conshdlrdata->nsols;
#endif
}


/** puts the number of counted solutions in the given char* buffer */
void SCIPgetNCountedSolsstr(
   SCIP*                 scip,               /**< SCIP data structure */
   char**                buffer,             /**< buffer to store the number for counted solutions */
   int                   buffersize,         /**< buffer size */
   int*                  requiredsize        /**< pointer to store the required size */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   
   /* find the countsols constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

#ifdef WITH_GMP
   *requiredsize = mpz_sizeinbase( conshdlrdata->nsols, 10 );
   toString(conshdlrdata->nsols, buffer);
#else
   if( conshdlrdata->nsols < pow(10.0, buffersize) )
   {
      toString(conshdlrdata->nsols, buffer);
      *requiredsize = strlen(*buffer);
   }
   else
      *requiredsize = 21;
   
#endif
}


/** returns number of counted non trivial feasible subtrees */
SCIP_Longint SCIPgetNCountedFeasSubtrees(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   
   /* find the countsols constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );
   
   return conshdlrdata->feasST;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop_vbounds.c
 * @ingroup PROPAGATORS
 * @brief  variable upper and lower bound propagator
 * @author Stefan Heinz
 * @author Jens Schulz
 *
 * This propagator uses the variable lower and upper bounds of a variable to reduce variable domains. We (implicitly)
 * create a graph for the variable lower and upper bounds. 
 *
 * 1) Graph construction
 *
 *    For each variable we create a node and for each variable lower (upper) bound we insert an arc (directed) from the
 *    variable which influences the lower (upper) bound of the other variable
 *
 * 2) Create a topological sorted variable array 
 *
 *    This graph is used to create two (almost) topological sorted variable array. One w.r.t. the variable lower bounds
 *    and the other w.r.t. the variable upper bounds. Topological sorted means, a variable which influences the lower
 *    (upper) bound of another variable y is located before y in the corresponding variable array. Note, that in general
 *    a topological sort is not unique.
 *
 * 3) Propagation
 *  
 *    The topological sorted lower and upper bound arrays are used to propagate the variable lower or upper bounds of
 *    the corresponding variables.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/prop_vbounds.h"


#define PROP_NAME              "vbounds"
#define PROP_DESC              "propagates variable upper and lower bounds"
#define PROP_PRIORITY           2000000
#define PROP_FREQ                     1
#define PROP_DELAY                FALSE      /**< should propagation method be delayed, if other propagators found reductions? */

#define EVENTHDLR_NAME         "vbounds"
#define EVENTHDLR_DESC         "bound change event handler for vbounds propagator"



/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   SCIP_VAR**            vars;               /**< array of involved variables */
   SCIP_HASHMAP*         varHashmap;         /**< mapping a variable to its posiotion in the variable array */    
   SCIP_VAR**            lbvars;             /**< topological sorted variables with respect to the variable lower bound */
   SCIP_VAR**            ubvars;             /**< topological sorted variables with respect to the variable upper bound */
   int                   nvars;              /**< number of involved variables */
   int                   nlbvars;            /**< number of variables in variable lower bound array */
   int                   nubvars;            /**< number of variables in variable upper bound array */
   unsigned int          sorted:1;           /**< is the variable array topological sorted */
   unsigned int          propagated:1;       /**< is the lower and upper bound variable array already propagated? */
};


/*
 * Hash map callback methods
 */

/** hash key retrieval function for variables */
static
SCIP_DECL_HASHGETKEY(hashGetKeyVar)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff the indices of both variables are equal */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqVar)
{  /*lint --e{715}*/
   if ( key1 == key2 )
      return TRUE;
   return FALSE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValVar)
{  /*lint --e{715}*/
   assert( SCIPvarGetIndex((SCIP_VAR*) key) >= 0 );
   return (unsigned int) SCIPvarGetIndex((SCIP_VAR*) key);
}


/*
 * Local methods
 */

/** perform depth-first-search from the given variable using the variable lower or upper bounds of the variable */
static
SCIP_RETCODE depthFirstSearch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to start the depth-first-search  */
   SCIP_HASHMAP*         varHashmap,         /**< mapping a variable to its posiotion in the (used) variable array, or NULL */    
   SCIP_VAR**            usedvars,           /**< array of variables which are involved in the propagation, or NULL */
   int*                  nusedvars,          /**< number of variables which are involved in the propagation, or NULL */
   SCIP_HASHTABLE*       connected,          /**< hash table storing if a node was already visited */
   SCIP_VAR**            sortedvars,         /**< array that will contain the topological sorted variables */
   int*                  nsortedvars,        /**< pointer to store the number of already collectes variables in the sorted variables array */
   SCIP_Bool             lowerbound          /**< depth-first-search with respect to the variable lower bounds, otherwise variable upper bound */
   )
{
   SCIP_VAR** vbvars;
   SCIP_VAR* vbvar;
   SCIP_Real scalar;
   SCIP_Real constant;
   int nvbvars;
   int v;

   assert(scip != NULL);
   assert(var != NULL);
   assert(varHashmap == NULL || (varHashmap != NULL && usedvars != NULL && nusedvars != NULL));
   assert(sortedvars != NULL);
   assert(nsortedvars != NULL);
   assert(*nsortedvars >= 0);
   assert(SCIPvarGetProbindex(var) > -1);
   assert(SCIPhashtableExists(connected, var));

   /* mark variable as visited, remove variable from hash table */
   SCIP_CALL( SCIPhashtableRemove(connected, var) );

   if( lowerbound )
   {
      /* get variable lower bounds */
      vbvars = SCIPvarGetVlbVars(var);
      nvbvars = SCIPvarGetNVlbs(var);
   }
   else
   {
      /* get variable upper bounds */
      vbvars = SCIPvarGetVubVars(var);
      nvbvars = SCIPvarGetNVubs(var);
   }
   
   SCIPdebugMessage("variable <%s> has %d variable %s bounds\n", SCIPvarGetName(var), nvbvars, 
      lowerbound ? "lower" : "upper");

   for( v = 0; v < nvbvars; ++v )
   {
      vbvar = vbvars[v];
      assert(vbvar != NULL);
      
      /* transform variable bound variable to an active variable if possible */
      SCIP_CALL( SCIPvarGetProbvarSum(&vbvar, &scalar, &constant) );
      
      /* we could not resolve the variable bound variable to one active variable, therefore, ignore this variable bound */
      if( !SCIPvarIsActive(vbvar) )
         continue;
      
      /* insert variable bound variable into the hash table since they are involve in the later propagation */
      if( varHashmap != NULL && !SCIPhashmapExists(varHashmap, vbvar) )
      {
         SCIPdebugMessage("insert variable <%s> with position %d into the hash map\n", SCIPvarGetName(vbvar), *nusedvars);
         SCIP_CALL( SCIPhashmapInsert(varHashmap, vbvar, (void*) (size_t)(*nusedvars)) );
         usedvars[*nusedvars] =  vbvar;
         (*nusedvars)++;
      }
      
      /* check if the variable bound variable was already visited */
      if( SCIPhashtableExists(connected, vbvar) )
      {
         /* recursively call depth-first-search */
         SCIP_CALL( depthFirstSearch(scip, vbvar, varHashmap, usedvars, nusedvars, connected, sortedvars, nsortedvars, lowerbound) );
      }
   }

   /* store variable in the sorted variable array */
   sortedvars[(*nsortedvars)] = var;
   
   (*nsortedvars)++;
   
   return SCIP_OKAY;
}

/** catches events for variables */
static
SCIP_RETCODE catchEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   int v;

   assert(propdata != NULL);

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   /**@todo try to be more precisely be selecting the variables events 
    *
    * - in case of the lower bound variables we have to watch the LBCHANGED if the variable bound coefficient is
    *   positive and UBCHANGED if the variable bound coefficient is negative 
    * - in case of the upper bound variables we have to watch the UBCHANGED if the variable bound coefficient is
    *   positive and LBCHANGED if the variable bound coefficient is negative 
    */

   /* catch for eached involved variable nound change events */
   for( v = 0; v < propdata->nvars; ++v )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, propdata->vars[v], SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, 
            eventhdlr, (SCIP_EVENTDATA*)propdata, NULL) );
   }
   
   return SCIP_OKAY;
}

/** drops events for variables */
static
SCIP_RETCODE dropEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   int v;

   assert(propdata != NULL);

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   /* drop events on involved variables */
   for( v = 0; v < propdata->nvars; ++v )
   {
      assert(SCIPvarIsTransformed(propdata->vars[v]));
      SCIP_CALL( SCIPdropVarEvent(scip, propdata->vars[v], SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, 
            eventhdlr, (SCIP_EVENTDATA*)propdata, -1) );
   }
   
   return SCIP_OKAY;
}

/** performs propagation of variables lower and upper bounds */
static
SCIP_RETCODE propagateVbounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< vbounds propagator */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_VAR** vars;
   SCIP_VAR** vbvars;
   SCIP_VAR* var;
   SCIP_VAR* vbvar;
   SCIP_Real* coefs;
   SCIP_Real* constants;
   SCIP_Real coef;
   SCIP_Real constant;
   SCIP_Real newbound;
   int inferinfo;
   int nvars;
   int nvbvars;
   int n;
   int v;
   int nchgbds;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;
   
   assert(scip != NULL);
   assert(prop != NULL);
   assert(result != NULL);

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   
   (*result) = SCIP_DIDNOTRUN;

   if( propdata->propagated )
      return SCIP_OKAY;
   
   nchgbds = 0;
   nvars = propdata->nlbvars;
   
   if( nvars > 0 )
   {
      vars = propdata->lbvars;
      assert(vars != NULL);
    
      SCIPdebugMessage("run vbounds (lower) propagator over %d variables\n", nvars);
      
      /* try tighten lower bounds by traversing topological sorted variables from left to right */
      for( v = 0; v < nvars; ++v )
      {
         assert(vars[v] != NULL);
         
         /* get next variable of the topological sorted graph */
         var = vars[v];
         assert(var != NULL );
         
         /* get current lower bound as initialization of new lower bound */
         newbound = SCIPvarGetLbLocal(var);
         inferinfo = 0;
         
         /* get the variable lower bound informations for the current variable */
         vbvars = SCIPvarGetVlbVars(var);
         coefs = SCIPvarGetVlbCoefs(var);
         constants = SCIPvarGetVlbConstants(var);
         nvbvars = SCIPvarGetNVlbs(var);

         /* loop over all variable lower bounds; a variable lower bound has the form: x >= b*y + d*/
         for( n = 0; n < nvbvars; ++n )
         {
            vbvar = vbvars[n];
            coef = coefs[n];
            constant = constants[n];
            
            /* transform variable bound variable to an active variable if possible */
            SCIP_CALL( SCIPvarGetProbvarSum(&vbvar, &coef, &constant) );
         
            if( !SCIPvarIsActive(vbvar) )
               continue;
         
            if( SCIPisPositive(scip, coef) )
            {
               /* if b > 0 => x >= b*lb(y) + d */ 
               if( SCIPisGT(scip, coef*SCIPvarGetLbLocal(vbvar) + constant, newbound) )
               {
                  assert(SCIPvarGetProbindex(vbvar) > -1);
                  newbound = coef*SCIPvarGetLbLocal(vbvar) + constant;
               
                  assert(SCIPhashmapExists(propdata->varHashmap, vbvar));
                  inferinfo = (int)(size_t)SCIPhashmapGetImage(propdata->varHashmap, vbvar) + 1;
               }
            }
            else
            {
               /* if b < 0 => x >= b*ub(y) + d */ 
               if( SCIPisGT(scip, coef*SCIPvarGetUbLocal(vbvar) + constant, newbound) )
               {
                  assert(SCIPvarGetProbindex(vbvar) > -1);
                  newbound = coef*SCIPvarGetUbLocal(vbvar) + constant;

                  assert(SCIPhashmapExists(propdata->varHashmap, vbvar));
                  inferinfo = -((int)(size_t)SCIPhashmapGetImage(propdata->varHashmap, vbvar) + 1);
               }
            }
         }
         
         if(inferinfo != 0)
         {
            SCIPdebugMessage(" -> variable <%s> => lower bound candidate is <%.15g>\n", SCIPvarGetName(SCIPgetVars(scip)[ABS(inferinfo)-1]), newbound);
         }
         
         SCIP_CALL( SCIPinferVarLbProp(scip, var, newbound, prop, inferinfo, &infeasible, &tightened) );
         
         if( infeasible )
         {
            SCIPdebugMessage(" -> lower bound tightening lead to infeasiility\n");
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;      
         } 
         
         if( tightened )
         {
            SCIPdebugMessage(" -> tightened lower bound to <%g> due the %s bound of variable <%s>\n", 
               newbound, inferinfo > 0 ? "lower" : "upper", SCIPvarGetName(SCIPgetVars(scip)[ABS(inferinfo)-1]));
            assert(inferinfo != 0);
            nchgbds++;
         }
         else if( inferinfo != 0 )
         {
            SCIPdebugMessage(" -> new bound was not good enough\n");
         }
      }
   }
   
   nvars = propdata->nubvars;

   if( nvars > 0 )
   {
      vars = propdata->ubvars;
      assert(vars != NULL);

      SCIPdebugMessage("run vbounds (upper) propagator over %d variables\n", nvars);

      /* try to tighten upper bounds by traversing topological sorted variables from right to left */
      for( v = 0; v < nvars; ++v )
      {
         assert(vars[v] != NULL);

         var = vars[v];
         assert(var != NULL);
      
         /* get current upper bound and initialize new upper bound */
         newbound = SCIPvarGetUbLocal(var);
         inferinfo = 0;
      
         /* loop over successor variables to find a better upper bound */
         vbvars = SCIPvarGetVubVars(var);
         coefs = SCIPvarGetVubCoefs(var);
         constants = SCIPvarGetVubConstants(var);
         nvbvars = SCIPvarGetNVubs(var);

         /* loop of the entering arcs of the current node */
         for( n = 0; n < nvbvars; ++n )
         {
            vbvar = vbvars[n];
            coef = coefs[n];
            constant = constants[n];

            /* transform variable bound variable to an active variable if possible */
            SCIP_CALL( SCIPvarGetProbvarSum(&vbvar, &coef, &constant) );

            if( !SCIPvarIsActive(vbvar) )
               continue;

            if( SCIPisPositive(scip, coef) )
            {
               if( SCIPisLT(scip, coef*SCIPvarGetUbLocal(vbvar) + constant, newbound) )
               {
                  /* if b > 0 => x <= b*ub(y) + d */ 
                  assert(SCIPvarGetProbindex(vbvar) > -1);
                  newbound = coef*SCIPvarGetUbLocal(vbvar) + constant;

                  assert(SCIPhashmapExists(propdata->varHashmap, vbvar));
                  inferinfo = (int)(size_t)SCIPhashmapGetImage(propdata->varHashmap, vbvar) + 1;
               }
            }
            else
            {
               if( SCIPisLT(scip, coef*SCIPvarGetLbLocal(vbvar) + constant, newbound) )
               {
                  /* if b < 0 => x <= b*lb(y) + d */ 
                  assert(SCIPvarGetProbindex(vbvar) > -1);
                  newbound = coef*SCIPvarGetLbLocal(vbvar) + constant;

                  assert(SCIPhashmapExists(propdata->varHashmap, vbvar));
                  inferinfo = -((int)(size_t)SCIPhashmapGetImage(propdata->varHashmap, vbvar) + 1);
               }
            }
         }
      
         if(inferinfo != 0)
         {
            SCIPdebugMessage(" -> variable <%s> => upper bound candidate is <%.15g>\n", SCIPvarGetName(SCIPgetVars(scip)[ABS(inferinfo)-1]), newbound);
         }

         /* try new upper bound */
         SCIP_CALL( SCIPinferVarUbProp(scip, var, newbound, prop, inferinfo, &infeasible, &tightened) );
      
         if( infeasible )
         {
            SCIPdebugMessage(" -> upper bound tightening lead to infeasiility\n");
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;      
         } 

         if( tightened )
         {
            SCIPdebugMessage(" -> tightened upper bound to <%g> due the %s bound of variable <%s>\n", 
               newbound, inferinfo < 0 ? "lower" : "upper", SCIPvarGetName(SCIPgetVars(scip)[ABS(inferinfo)-1]));
            assert(inferinfo != 0);
            nchgbds++;
         }
         else if( inferinfo != 0 )
         {
            SCIPdebugMessage(" -> new bound was not good enough\n");
         }
      }
   }   

   /* mark lower and upper bound variable array as propagated */
   propdata->propagated = TRUE;
   
   SCIPdebugMessage("tightened %d variable bounds\n", nchgbds);

   if( nchgbds > 0 )
      (*result) = SCIP_REDUCEDDOM;
   else
      (*result) = SCIP_DIDNOTFIND;
   
   return SCIP_OKAY;
}

/*
 * Callback methods of propagator
 */

/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeVbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   /* free propagator data */
   propdata = SCIPpropGetData(prop);

   SCIPfreeMemory(scip, &propdata);
   SCIPpropSetData(prop, NULL);
   
   return SCIP_OKAY;
}


/** initialization method of propagator (called after problem was transformed) */
#define propInitVbounds NULL

/** deinitialization method of propagator (called before transformed problem is freed) */
#define propExitVbounds NULL

/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
static
SCIP_DECL_PROPINITSOL(propInitsolVbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_VAR** vars;
   int nvars;
   int v;

   SCIPdebugMessage("initialize prop_vbounds propagator for problem <%s>\n", SCIPgetProbName(scip));

   /* free propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   
   /* allocate memory for the arrays of the propdata */
   SCIP_CALL( SCIPallocMemoryArray(scip, &propdata->vars, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &propdata->lbvars, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &propdata->ubvars, nvars) );

   /* create hash table for storing the involved variables */
   assert(propdata->nvars == 0);
   SCIP_CALL( SCIPhashmapCreate(&propdata->varHashmap, SCIPblkmem(scip), SCIPcalcHashtableSize(nvars)) );
   
   /* create the topological sorted variable array with respect to the variable lower bounds */
   assert(propdata->nlbvars == 0);
   SCIP_CALL( SCIPcreateTopoSortedVars(scip, vars, nvars, propdata->varHashmap, propdata->vars, &propdata->nvars, 
         propdata->lbvars, &propdata->nlbvars, TRUE) );

   /* create the topological sorted variable array with respect to the variable upper bounds */
   assert(propdata->nubvars == 0);
   SCIP_CALL( SCIPcreateTopoSortedVars(scip, vars, nvars, propdata->varHashmap, propdata->vars, &propdata->nvars, 
         propdata->ubvars, &propdata->nubvars, FALSE) );

   /* capture all variables */
   for( v = 0; v < propdata->nvars; ++v )
   {
      SCIP_CALL( SCIPcaptureVar(scip, propdata->vars[v]) );
   }
   for( v = 0; v < propdata->nlbvars; ++v )
   {
      SCIP_CALL( SCIPcaptureVar(scip, propdata->lbvars[v]) );
   }
   for( v = 0; v < propdata->nubvars; ++v )
   {
      SCIP_CALL( SCIPcaptureVar(scip, propdata->ubvars[v]) );
   }

   /* catch variable events */
   SCIP_CALL( catchEvents(scip, propdata) );

   propdata->propagated = FALSE;
   propdata->sorted = TRUE;

   return SCIP_OKAY;
}


/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
static
SCIP_DECL_PROPEXITSOL(propExitsolVbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   int v;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* drop all variable events */
   SCIP_CALL( dropEvents(scip, propdata) );

   /* release all variables */
   for( v = 0; v < propdata->nvars; ++v )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &propdata->vars[v]) );
   }
   for( v = 0; v < propdata->nlbvars; ++v )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &propdata->lbvars[v]) );
   }
   for( v = 0; v < propdata->nubvars; ++v )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &propdata->ubvars[v]) );
   }

   /* free hash map */
   SCIPhashmapFree(&propdata->varHashmap);

   /* free varbounds array */
   SCIPfreeMemoryArrayNull(scip, &propdata->lbvars);
   SCIPfreeMemoryArrayNull(scip, &propdata->ubvars);
   SCIPfreeMemoryArrayNull(scip, &propdata->vars);

   propdata->nvars = 0;
   propdata->nlbvars = 0;
   propdata->nubvars = 0;

   return SCIP_OKAY;
}


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecVbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   
   /* perform variable lower and upper bound propagation */
   SCIP_CALL( propagateVbounds(scip, prop, result) );
   
   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropVbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_VAR* var;
   int idx;

   assert(inferinfo != 0);
   assert(ABS(inferinfo) - 1 < SCIPgetNVars(scip));
   
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   
   idx = ABS(inferinfo) - 1;
   assert(idx >= 0);
   assert(idx < propdata->nvars);
   
   var = propdata->vars[idx];

   switch( boundtype )
   {
   case SCIP_BOUNDTYPE_LOWER:
      if( inferinfo > 0 )
      {
         SCIP_CALL( SCIPaddConflictLb(scip, var, bdchgidx ) );
      }
      else
      {
         SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx ) );
      }
      break;
   case SCIP_BOUNDTYPE_UPPER:
      if( inferinfo < 0 )
      {
         SCIP_CALL( SCIPaddConflictLb(scip, var, bdchgidx ) );
      }
      else
      {
         SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx ) );
      }
      break;
   default:
      SCIPerrorMessage("invalid bound type <%d>\n", boundtype);
      SCIPABORT();
   }

   (*result) = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/*
 * Event Handler
 */

/** execution methode of bound change event handler */
static
SCIP_DECL_EVENTEXEC(eventExecVbound)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = (SCIP_PROPDATA*)eventdata;
   assert(propdata != NULL);

   propdata->propagated = FALSE;

   return SCIP_OKAY;
}



/*
 * propagator specific interface methods
 */

/** creates the vbounds propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropVbounds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   
   /* create pseudoobj propagator data */
   SCIP_CALL( SCIPallocMemory(scip, &propdata) );
   propdata->vars = NULL;
   propdata->lbvars = NULL;
   propdata->ubvars = NULL;
   propdata->nvars = 0;
   propdata->nlbvars = 0;
   propdata->nubvars = 0;
   propdata->propagated = FALSE;
   propdata->sorted = FALSE;

   /* include propagator */
   SCIP_CALL( SCIPincludeProp(scip, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY,
         propFreeVbounds, propInitVbounds, propExitVbounds, 
         propInitsolVbounds, propExitsolVbounds, propExecVbounds, propRespropVbounds,
         propdata) );

   /* include event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC, 
         NULL, NULL, NULL, NULL, NULL, NULL, eventExecVbound, NULL) );

   return SCIP_OKAY;
}

/** create a topological sorted variable array of the given variables and stores if (needed) the involved variables into
 *  the corresponding variable array and hash map
 *
 * @note: for all arrays and the hash map (if needed) you need to allocate enough memory before calling this method 
 */
SCIP_RETCODE SCIPcreateTopoSortedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variable which we want sort */
   int                   nvars,              /**< number of variables */
   SCIP_HASHMAP*         varHashmap,         /**< mapping a variable to its posiotion in the (used) variable array, or NULL */    
   SCIP_VAR**            usedvars,           /**< array of variables which are involved in the propagation, or NULL */
   int*                  nusedvars,          /**< number of variables which are involved in the propagation, or NULL */
   SCIP_VAR**            topovars,           /**< array where the topological sorted variables are stored */
   int*                  ntopovars,          /**< pointer to store the number of topological sorted variables */
   SCIP_Bool             lowerbound          /**< topological sorted with respect to the variable lower bounds, otherwise variable upper bound */
   )
{
   SCIP_VAR** sortedvars;
   SCIP_VAR** vbvars;
   SCIP_VAR* var;
   SCIP_HASHTABLE* connected;
   int nvbvars;
   int nsortedvars;
   int i;
   int v;
   
   assert(scip != NULL);   
   assert(vars != NULL || nvars == 0);
   assert(varHashmap == NULL || (varHashmap != NULL && usedvars != NULL && nusedvars != NULL));
   assert(topovars != NULL);
   assert(ntopovars != NULL);
   
   SCIPdebugMessage("create topological sorted variable array with respect to variables %s bounds\n", 
      lowerbound ? "lower" : "upper");

   if( nvars == 0 )
      return SCIP_OKAY;

   assert(vars != NULL);
   
   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &sortedvars, nvars) );
   
   /* create hash table for variables whcih are (still) connected */
   SCIP_CALL( SCIPhashtableCreate(&connected, SCIPblkmem(scip), SCIPcalcHashtableSize(nvars), hashGetKeyVar, hashKeyEqVar, hashKeyValVar, NULL) );
   
   /* detect isolated variables; mark all variables which have at least one entering or leaving arc as connected */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert(var != NULL);
      
      if( lowerbound )
      {
         /* get variable lower bounds */
         vbvars = SCIPvarGetVlbVars(var);
         nvbvars = SCIPvarGetNVlbs(var);
      }
      else
      {
         /* get variable upper bounds */
         vbvars = SCIPvarGetVubVars(var);
         nvbvars = SCIPvarGetNVubs(var);
      }
      
      if( nvbvars > 0 && !SCIPhashtableExists(connected, var) )
      {
         SCIP_CALL( SCIPhashtableInsert(connected, var) );
      }

      for( i = 0; i < nvbvars; ++i )
      {
         /* there is a leaving arc, hence, the variable is connected */  
         assert(vbvars[i] != NULL);
         if( !SCIPhashtableExists(connected, vbvars[i]) )
         {
            SCIP_CALL( SCIPhashtableSafeInsert(connected, vbvars[i]) );
         }
      }
   }

   /* loop over all "connected" variable and find for each connected component a "almost" topological sorted version */
   for( v = 0; v < nvars; ++v )
   {
      if( SCIPhashtableExists(connected, vars[v]) )
      {
         SCIPdebugMessage("start depth-first-search with variable <%s>\n", SCIPvarGetName(vars[v]));
         
         /* use depth first search to get a "almost" topological sorted variables for the connected component which
          * includes vars[v] */
         nsortedvars = 0;
         SCIP_CALL( depthFirstSearch(scip, vars[v], varHashmap, usedvars, nusedvars, connected, 
               sortedvars, &nsortedvars, lowerbound) );
         
         SCIPdebugMessage("detected connected component of size <%d>\n", nsortedvars);
            
         /* */
         if( nsortedvars == 1 )
         {
            assert(!SCIPhashtableExists(connected, vars[v]));
            SCIP_CALL( SCIPhashtableInsert(connected, vars[v]) );
         }
         else
         {
            /* copy and capture variables to make sure, the variables are not deleted */
            for( i = 0; i < nsortedvars; ++i )
            {
               topovars[(*ntopovars)] = sortedvars[i];
               (*ntopovars)++;
            }
         }
      }
   }
   
   assert(*ntopovars <= nvars);
   SCIPdebugMessage("topological sorted array contains %d of %d variables (variable %s bound)\n", 
      *ntopovars, nvars, lowerbound ? "lower" : "upper");
   
   /* free hash table */
   SCIPhashtableFree(&connected);

   /* free buffer memory */
   SCIPfreeBufferArray(scip, &sortedvars);
   
   return SCIP_OKAY;
}

/** returns TRUE if the propagater has the status that all variable lower and upper bounds are propgated */
SCIP_Bool SCIPisPropagatedVbounds(
   SCIP*                 scip                 /**< SCIP data structure */
   )
{
   SCIP_PROP* prop;
   SCIP_PROPDATA* propdata;
   
   prop = SCIPfindProp(scip, PROP_NAME);
   assert(prop != NULL);

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   return propdata->propagated;
}

/** performs propagation of variables lower and upper bounds */
SCIP_RETCODE SCIPexecPropVbounds(
   SCIP*                 scip                 /**< SCIP data structure */
   )
{
   SCIP_PROP* prop;
   SCIP_RESULT result;
   
   prop = SCIPfindProp(scip, PROP_NAME);
   assert(prop != NULL);

   /* perform variable lower and upper bound propagation */
   SCIP_CALL( propagateVbounds(scip, prop, &result) );

   return SCIP_OKAY;
}

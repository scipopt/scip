/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop_vbounds.c
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
#include <string.h>

#include "scip/prop_vbounds.h"

#define PROP_NAME              "vbounds"
#define PROP_DESC              "propagates variable upper and lower bounds"
#define PROP_TIMING             SCIP_PROPTIMING_BEFORELP
#define PROP_PRIORITY           2000000 /**< propagator priority */ 
#define PROP_FREQ                     1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define PROP_PRESOL_PRIORITY          0 /**< priority of the presolving method (>= 0: before, < 0: after constraint handlers); combined with presolvers */
#define PROP_PRESOL_DELAY          TRUE /**< should presolving be delay, if other presolvers found reductions?  */
#define PROP_PRESOL_MAXROUNDS         0 /**< maximal number of presolving rounds the presolver participates in (-1: no
                                         *   limit) */

#define EVENTHDLR_NAME         "vbounds"
#define EVENTHDLR_DESC         "bound change event handler for for vbounds propagator"

#define DEFAULT_USEBDWIDENING      TRUE      /**< should bound widening be used to initialize conflict analysis? */

/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   SCIP_VAR**            vars;               /**< array of involved variables */
   SCIP_HASHMAP*         varHashmap;         /**< mapping a variable to its position in the variable array */    
   SCIP_VAR**            lbvars;             /**< topological sorted variables with respect to the variable lower bound */
   SCIP_VAR**            ubvars;             /**< topological sorted variables with respect to the variable upper bound */
   SCIP_EVENTTYPE*       lbeventtypes;       /**< event types of variables belonging to variable lower bounds */ 
   SCIP_EVENTTYPE*       ubeventtypes;       /**< event types of variables belonging to variable upper bounds */ 
   int                   nvars;              /**< number of involved variables */
   int                   neventvars;         /**< number of variables which are triggered by an event */
   int                   sizevars;           /**< size of the variable array vars */
   int                   nlbvars;            /**< number of variables in variable lower bound array */
   int                   nubvars;            /**< number of variables in variable upper bound array */
   SCIP_Bool             lbpropagated;       /**< is the lower bound variable array already propagated? */
   SCIP_Bool             ubpropagated;       /**< is the upper bound variable array already propagated? */
   SCIP_Bool             usebdwidening;      /**< should bound widening be used to initialize conflict analysis? */
};


/** inference information */
struct InferInfo
{
   union
   {
      struct
      {
         unsigned int    pos:31;             /**< position of the variable which forced that propagation */
         unsigned int    boundtype:1;        /**< bound type which was the reason (0: lower, 1: upper) */
      } asbits;
      int                asint;              /**< inference information as a single int value */
   } val;
};
typedef struct InferInfo INFERINFO;

/** converts an integer into an inference information */
static
INFERINFO intToInferInfo(
   int                   i                   /**< integer to convert */
   )
{
   INFERINFO inferinfo;

   inferinfo.val.asint = i;

   return inferinfo;
}

/** converts an inference information into an int */
static
int inferInfoToInt(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return inferinfo.val.asint;
}

/** returns the propagation rule stored in the inference information */
static
SCIP_BOUNDTYPE inferInfoGetBoundtype(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   assert((SCIP_BOUNDTYPE)inferinfo.val.asbits.boundtype == SCIP_BOUNDTYPE_LOWER 
      || (SCIP_BOUNDTYPE)inferinfo.val.asbits.boundtype == SCIP_BOUNDTYPE_UPPER);
   return (SCIP_BOUNDTYPE)inferinfo.val.asbits.boundtype;
}

/** returns the position stored in the inference information */
static
int inferInfoGetPos(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return inferinfo.val.asbits.pos;
}

/** constructs an inference information out of a propagation rule and a position number */
static
INFERINFO getInferInfo(
   int                   pos,                /**< position of the variable which forced that propagation */
   SCIP_BOUNDTYPE        boundtype           /**< propagation rule that deduced the value */
   )
{
   INFERINFO inferinfo;

   assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);
   assert((int)boundtype >= 0 && (int)boundtype <= 1); /*lint !e685 !e568q*/
   
   inferinfo.val.asbits.pos = pos; /*lint !e732*/
   inferinfo.val.asbits.boundtype = boundtype; /*lint !e641*/
   
   return inferinfo;
}


/*
 * Local methods
 */

/** reset propagation data */
static
void resetPropdata(
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   propdata->vars = NULL;
   propdata->varHashmap = NULL;
   propdata->lbvars = NULL;
   propdata->ubvars = NULL;
   propdata->lbeventtypes = NULL;
   propdata->ubeventtypes = NULL;
   propdata->nvars = 0;
   propdata->neventvars = 0;
   propdata->nlbvars = 0;
   propdata->nubvars = 0;
   propdata->lbpropagated = TRUE;
   propdata->ubpropagated = TRUE;
}

/** gets the requested variables bounds */
static
void getVariableBounds(
   SCIP_VAR*             var,                /**< variable to get the variable bounds from */
   SCIP_VAR***           vbvars,             /**< pointer to store the variable bound array */
   int*                  nvbvars,            /**< pointer to store the number of variable bounds */
   SCIP_Bool             lowerbound          /**< variable lower bounds? (otherwise variable upper bound) */
   )
{
   if( lowerbound )
   {
      /* get variable lower bounds */
      (*vbvars) = SCIPvarGetVlbVars(var);
      (*nvbvars) = SCIPvarGetNVlbs(var);
   }
   else
   {
      /* get variable upper bounds */
      (*vbvars) = SCIPvarGetVubVars(var);
      (*nvbvars) = SCIPvarGetNVubs(var);
   }
}

/** perform depth-first-search from the given variable using the variable lower or upper bounds of the variable */
static
SCIP_RETCODE depthFirstSearch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to start the depth-first-search  */
   SCIP_HASHMAP*         varPosMap,          /**< mapping a variable to its position in the (used) variable array, or NULL */    
   SCIP_VAR**            usedvars,           /**< array of variables which are involved in the propagation, or NULL */
   int*                  nusedvars,          /**< number of variables which are involved in the propagation, or NULL */
   SCIP_HASHTABLE*       connected,          /**< hash table storing if a node was already visited */
   SCIP_VAR**            sortedvars,         /**< array that will contain the topological sorted variables */
   int*                  nsortedvars,        /**< pointer to store the number of already collects variables in the sorted variables array */
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
   assert(varPosMap == NULL || (varPosMap != NULL && usedvars != NULL && nusedvars != NULL));
   assert(sortedvars != NULL);
   assert(nsortedvars != NULL);
   assert(*nsortedvars >= 0);
   assert(SCIPvarGetProbindex(var) > -1);
   assert(SCIPhashtableExists(connected, var));

   /* mark variable as visited, remove variable from hash table */
   SCIP_CALL( SCIPhashtableRemove(connected, var) );

   /* get variable bounds */
   getVariableBounds(var, &vbvars, &nvbvars, lowerbound);
   
   SCIPdebugMessage("variable <%s> has %d variable %s bounds\n", SCIPvarGetName(var), nvbvars, 
      lowerbound ? "lower" : "upper");

   for( v = 0; v < nvbvars; ++v )
   {
      vbvar = vbvars[v];
      assert(vbvar != NULL);
      
      scalar = 1.0;
      constant = 0.0;

      /* transform variable bound variable to an active variable if possible */
      SCIP_CALL( SCIPvarGetProbvarSum(&vbvar, &scalar, &constant) );
      
      /* we could not resolve the variable bound variable to one active variable, therefore, ignore this variable bound */
      if( !SCIPvarIsActive(vbvar) )
         continue;
      
      /* insert variable bound variable into the hash table since they are involved in later propagation */
      if( varPosMap != NULL && !SCIPhashmapExists(varPosMap, vbvar) )
      {
         SCIPdebugMessage("insert variable <%s> with position %d into the hash map\n", SCIPvarGetName(vbvar), *nusedvars);
         SCIP_CALL( SCIPhashmapInsert(varPosMap, vbvar, (void*)(size_t)(*nusedvars)) );
         usedvars[*nusedvars] =  vbvar;
         (*nusedvars)++;
      }
      
      /* check if the variable bound variable was already visited */
      if( SCIPhashtableExists(connected, vbvar) )
      {
         /* recursively call depth-first-search */
         SCIP_CALL( depthFirstSearch(scip, vbvar, varPosMap, usedvars, nusedvars, connected, sortedvars, nsortedvars, lowerbound) );
      }
   }

   /* store variable in the sorted variable array */
   sortedvars[(*nsortedvars)] = var;
   (*nsortedvars)++;
   
   /* insert variable bound variable into the hash table since they are involve in the later propagation */ 
   if( varPosMap != NULL && !SCIPhashmapExists(varPosMap, var) ) 
   { 
      SCIPdebugMessage("insert variable <%s> with position %d into the hash map\n", SCIPvarGetName(var), *nusedvars); 
      SCIP_CALL( SCIPhashmapInsert(varPosMap, var, (void*) (size_t)(*nusedvars)) ); 
      usedvars[*nusedvars] =  var; 
      (*nusedvars)++; 
   } 

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
   SCIP_VAR** vbvars;
   SCIP_VAR* vbvar;
   SCIP_Real* coefs;
   SCIP_Real coef;
   SCIP_Real constant;
   int nvbvars;
   int n;
   SCIP_VAR* var;
   int idx;
   int v;

   assert(propdata != NULL);

   propdata->neventvars = propdata->nvars;

   /* setup arrays of event types lbeventtype and ubeventtype */
   if( propdata->nlbvars > 0 )
   {
      /* we watch the LBCHANGED if the variable bound coefficient is positive and 
       * we watch the UBCHANGED if the variable bound coefficient is negative
       */
      
      assert(propdata->lbeventtypes == NULL);
      SCIP_CALL( SCIPallocMemoryArray(scip, &propdata->lbeventtypes, propdata->neventvars) );
      BMSclearMemoryArray(propdata->lbeventtypes, propdata->neventvars);
      
      for( v = 0; v < propdata->nlbvars; ++v )
      {
         var = propdata->lbvars[v];

         vbvars  = SCIPvarGetVlbVars(var);
         coefs   = SCIPvarGetVlbCoefs(var);
         nvbvars = SCIPvarGetNVlbs(var);

         /* loop over all variable lower bounds; a variable lower bound has the form: x >= b*y + d*/
         for( n = 0; n < nvbvars; ++n )
         {
            vbvar = vbvars[n];
            coef = coefs[n];
            constant = 1.0;
                        
            /* transform variable bound variable to an active variable if possible */
            SCIP_CALL( SCIPvarGetProbvarSum(&vbvar, &coef, &constant) );
         
            if( !SCIPvarIsActive(vbvar) )
               continue;

            assert(SCIPhashmapExists(propdata->varHashmap, vbvar));
            idx = (int)(size_t)SCIPhashmapGetImage(propdata->varHashmap, vbvar);
            assert(idx < propdata->neventvars);
            
            if( coef > 0.0 )
            {
               /* change in lower bound of y may lead to a propagation for x */
               propdata->lbeventtypes[idx] = propdata->lbeventtypes[idx] | SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_VARFIXED;
            }
            else
            {
               /* change in upper bound of y may lead to a propagation for x */
               propdata->lbeventtypes[idx] = propdata->lbeventtypes[idx] | SCIP_EVENTTYPE_UBCHANGED | SCIP_EVENTTYPE_VARFIXED;
            }
         }
      }
   }
   
   if( propdata->nubvars > 0 )
   {
      /* we watch the UBCHANGED if the variable bound coefficient is positive and
       * we watch the LBCHANGED if the variable bound coefficient is negative
       */

      assert(propdata->ubeventtypes == NULL);
      SCIP_CALL( SCIPallocMemoryArray(scip, &propdata->ubeventtypes, propdata->neventvars) );
      BMSclearMemoryArray(propdata->ubeventtypes, propdata->neventvars);

      for( v = 0; v < propdata->nubvars; ++v )
      {
         var = propdata->ubvars[v];

         vbvars  = SCIPvarGetVubVars(var);
         coefs   = SCIPvarGetVubCoefs(var);
         nvbvars = SCIPvarGetNVubs(var);

         /* loop over all variable upper bounds; a variable upper bound has the form: x <= b*y + d*/
         for( n = 0; n < nvbvars; ++n )
         {
            vbvar = vbvars[n];
            coef = coefs[n];
            constant = 1.0;
                        
            /* transform variable bound variable to an active variable if possible */
            SCIP_CALL( SCIPvarGetProbvarSum(&vbvar, &coef, &constant) );
         
            if( !SCIPvarIsActive(vbvar) )
               continue;

            assert(SCIPhashmapExists(propdata->varHashmap, vbvar));
            idx = (int)(size_t)SCIPhashmapGetImage(propdata->varHashmap, vbvar);
            assert(idx < propdata->neventvars);
            
            if( coef > 0.0 )
            {
               /* change in upper bound of y may lead to a propagation for x */
               propdata->ubeventtypes[idx] = propdata->ubeventtypes[idx] | SCIP_EVENTTYPE_UBCHANGED | SCIP_EVENTTYPE_VARFIXED;
            }
            else
            {
               /* change in lower bound of y may lead to a propagation for x */
               propdata->ubeventtypes[idx] = propdata->ubeventtypes[idx] | SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_VARFIXED;
            }
         }
      }
   }
   
   /* catch variable events according to computed event types */
   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   for( v = 0; v < propdata->neventvars; ++v )
   {
      if( propdata->lbeventtypes != NULL && propdata->lbeventtypes[v] != SCIP_EVENTTYPE_DISABLED )
      {
         SCIP_CALL( SCIPcatchVarEvent(scip, propdata->vars[v], propdata->lbeventtypes[v], eventhdlr, (SCIP_EVENTDATA*)(&propdata->lbpropagated), NULL) );
      }

      if( propdata->ubeventtypes != NULL && propdata->ubeventtypes[v] != SCIP_EVENTTYPE_DISABLED )
      {
         SCIP_CALL( SCIPcatchVarEvent(scip, propdata->vars[v], propdata->ubeventtypes[v], eventhdlr, (SCIP_EVENTDATA*)(&propdata->ubpropagated), NULL) );
      }
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
   assert((propdata->nlbvars == 0) == (propdata->lbeventtypes == NULL));
   assert((propdata->nubvars == 0) == (propdata->ubeventtypes == NULL));

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   for( v = 0; v < propdata->neventvars; ++v )
   {
      if( propdata->lbeventtypes != NULL && propdata->lbeventtypes[v] != SCIP_EVENTTYPE_DISABLED )
      {
         SCIP_CALL( SCIPdropVarEvent(scip, propdata->vars[v], propdata->lbeventtypes[v], eventhdlr, (SCIP_EVENTDATA*)(&propdata->lbpropagated), -1) );
      }
      if( propdata->ubeventtypes != NULL && propdata->ubeventtypes[v] != SCIP_EVENTTYPE_DISABLED )
      {
         SCIP_CALL( SCIPdropVarEvent(scip, propdata->vars[v], propdata->ubeventtypes[v], eventhdlr, (SCIP_EVENTDATA*)(&propdata->ubpropagated), -1) );
      }
   }
  
   return SCIP_OKAY;
}

/** resolves a propagation by adding the variable which implied that bound change */
static
SCIP_RETCODE resolvePropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   INFERINFO             inferinfo,          /**< inference information */
   SCIP_BDCHGIDX*        bdchgidx            /**< the index of the bound change, representing the point of time where the change took place */
   )
{
   SCIP_VAR* var;
   SCIP_BOUNDTYPE boundtype;
   int pos;

   assert(propdata != NULL);
   
   boundtype = inferInfoGetBoundtype(inferinfo);
   assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);
   
   pos = inferInfoGetPos(inferinfo);
   assert(pos >= 0);
   assert(pos < propdata->nvars);
   
   var = propdata->vars[pos];
   
   SCIPdebugMessage(" -> add %s bound of variable <%s> as reason\n", 
      boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", SCIPvarGetName(var));
   
   switch( boundtype )
   {
   case SCIP_BOUNDTYPE_LOWER:
   {
      SCIP_CALL( SCIPaddConflictLb(scip, var, bdchgidx ) );
      break;
   }
   case SCIP_BOUNDTYPE_UPPER:
      SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx ) );
      break;
   default:
      SCIPerrorMessage("invalid bound type <%d>\n", boundtype);
      SCIPABORT();
   }

   return SCIP_OKAY;
}

/** relax lower bound of give variable as long as the given inference upper bound leads still to a cutoff and add that
 *  bound change to the conflict set
 */
static
SCIP_RETCODE relaxInfervarLowerbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable for which the upper bound should be relaxed */
   SCIP_Real             inferub,            /**< upper bound which lead to infeasibility */
   SCIP_Real*            newlb               /**< pointer to store the reached relaxed lower bound */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_BDCHGIDX* bdchgidx;
   int nbdchgs;

   /* get number of bound changes */
   nbdchgs = SCIPvarGetNBdchgInfosLb(var);
   bdchgidx = NULL;
   
   assert(nbdchgs > 0 || SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetLbGlobal(var)));
   assert(nbdchgs == 0 || SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPbdchginfoGetNewbound(SCIPvarGetBdchgInfoLb(var, nbdchgs-1))));

   SCIPdebugMessage("variable <%s>[%.15g,%.15g]: nbdchgs %d try to relax lower bound to %.15g\n", 
      SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), nbdchgs, inferub);

   /* try to relax lower bound */
   while( nbdchgs > 0 )
   {
      bdchginfo = SCIPvarGetBdchgInfoLb(var, nbdchgs-1);
      assert(SCIPisEQ(scip, *newlb, SCIPbdchginfoGetNewbound(bdchginfo)));

      SCIPdebugMessage("lower bound change %d oldbd=%.15g, newbd=%.15g, depth=%d, pos=%d, redundant=%u\n", 
         nbdchgs, SCIPbdchginfoGetOldbound(bdchginfo), SCIPbdchginfoGetNewbound(bdchginfo), 
         SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo), SCIPbdchginfoIsRedundant(bdchginfo));

      /* check if the old lower bound is sufficient to prove infeasibility; in case the inference upper bound is
       * greater equal to the next possible relaxed lower bound, then we have to break since in this case the inference
       * upper bound does not lead to a cutoff anymore
       */
      if( SCIPisGE(scip, inferub, SCIPbdchginfoGetOldbound(bdchginfo)) )
         break;
      
      SCIPdebugMessage("***** relaxed lower bound of inference variable <%s> from <%g> to <%g>\n", 
         SCIPvarGetName(var), SCIPbdchginfoGetNewbound(bdchginfo), SCIPbdchginfoGetOldbound(bdchginfo));

      bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
      *newlb = SCIPbdchginfoGetOldbound(bdchginfo);
      nbdchgs--;
   }
      
   /* if the nbdchgs is zero then the local bound matches the global bound, therefore bdchgidx equal to NULL represents
    * the right time point and SCIP finds out that this bound is redundant since it is global
    */
   SCIPdebugMessage("add lower bound of bound change info %d to conflict set\n", nbdchgs);
   SCIP_CALL( SCIPaddConflictLb(scip, var, bdchgidx) );
   
   SCIPdebugMessage("relaxed lower bound to %.15g\n", *newlb);

   return SCIP_OKAY;
}

/** relax lower bound of give variable as long as the given inference bound leads still to a cutoff and add that bound
 *  change to the conflict set
 */
static
SCIP_RETCODE relaxVbdvarLowerbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable for which the upper bound should be relaxed */
   SCIP_Real             coef,               /**< variable bound coefficient */
   SCIP_Real             constant,           /**< variable bound constant */
   SCIP_Real             bound,              /**< bound to exceed */
   SCIP_Real*            inferbound          /**< pointer to store the relaxed inferbound bound */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_BDCHGIDX* bdchgidx;
   int nbdchgs;

   /* get number of bound changes */
   nbdchgs = SCIPvarGetNBdchgInfosLb(var);
   bdchgidx = NULL;
 
   SCIPdebugMessage("variable <%s>[%.15g,%.15g], coef=%.15g, constant=%.15g: nbdchgs %d try to relax lower bound to %.15g\n", 
      SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), coef, constant, nbdchgs, bound);
 
   assert(nbdchgs > 0 || SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetLbGlobal(var)));
   assert(nbdchgs == 0 || SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPbdchginfoGetNewbound(SCIPvarGetBdchgInfoLb(var, nbdchgs-1))));

   /* try to relax lower bound */
   while( nbdchgs > 0 )
   {
      bdchginfo = SCIPvarGetBdchgInfoLb(var, nbdchgs-1);

      SCIPdebugMessage("lower bound change %d oldbd=%.15g, newbd=%.15g, depth=%d, pos=%d, redundant=%u\n", 
         nbdchgs, SCIPbdchginfoGetOldbound(bdchginfo), SCIPbdchginfoGetNewbound(bdchginfo), 
         SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo), SCIPbdchginfoIsRedundant(bdchginfo));
      
      /* check if the old lower bound is sufficient to prove infeasibility; in case the inference bound is greater
       * equal to the next possible relaxed lower bound, then we have to break since in this case the inference bound
       * does not lead to a cutoff anymore
       */
      if( SCIPisGE(scip, bound, coef * SCIPbdchginfoGetOldbound(bdchginfo) + constant) )
         break;
      
      SCIPdebugMessage("***** relaxed lower bound of vbound variable <%s> from <%g> to <%g>\n", 
         SCIPvarGetName(var), SCIPbdchginfoGetNewbound(bdchginfo), SCIPbdchginfoGetOldbound(bdchginfo));

      bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
      *inferbound =  coef * SCIPbdchginfoGetOldbound(bdchginfo) + constant;
      nbdchgs--;
   }
      
   /* if the nbdchgs is zero then the local bound matches the global bound, therefore bdchgidx equal to NULL represents
    * the right time point and SCIP finds out that this bound is redundant since it is global
    */
   SCIPdebugMessage("add lower bound of bound change info %d to conflict set\n", nbdchgs); 
   SCIP_CALL( SCIPaddConflictLb(scip, var, bdchgidx) );
   
   return SCIP_OKAY;
}

/** relax upper bound of give variable as long as the given inference lower bound leads still to a cutoff and add that
 *  bound change to the conflict set
 */
static
SCIP_RETCODE relaxInfervarUpperbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable for which the upper bound should be relaxed */
   SCIP_Real             inferlb,            /**< lower bound which lead to infeasibility */
   SCIP_Real*            newub               /**< pointer to store the reached relaxed upper bound */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_BDCHGIDX* bdchgidx;
   int nbdchgs;

   /* get number of bound changes */
   nbdchgs = SCIPvarGetNBdchgInfosUb(var);
   bdchgidx = NULL;
      
   assert(nbdchgs > 0 || SCIPisEQ(scip, SCIPvarGetUbLocal(var), SCIPvarGetUbGlobal(var)));
   assert(nbdchgs == 0 || SCIPisEQ(scip, SCIPvarGetUbLocal(var), SCIPbdchginfoGetNewbound(SCIPvarGetBdchgInfoUb(var, nbdchgs-1))));

   SCIPdebugMessage("variable <%s>[%.15g,%.15g]: nbdchgs %d try to relax upper bound to %.15g\n", 
      SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), nbdchgs, inferlb);
   
   /* try to relax upper bound */
   while( nbdchgs > 0 )
   {
      bdchginfo = SCIPvarGetBdchgInfoUb(var, nbdchgs-1);
      assert(SCIPisEQ(scip, *newub, SCIPbdchginfoGetNewbound(bdchginfo)));

      SCIPdebugMessage("upper bound change %d oldbd=%.15g, newbd=%.15g, depth=%d, pos=%d, redundant=%u\n", 
         nbdchgs, SCIPbdchginfoGetOldbound(bdchginfo), SCIPbdchginfoGetNewbound(bdchginfo), 
         SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo), SCIPbdchginfoIsRedundant(bdchginfo));
   
      /* check if the old upper bound is sufficient to prove infeasibility; in case the inference lower bound is less
       * equal to the next possible relaxed upper bound, then we have to break since in this case the inference lower bound
       * does not lead to a cutoff anymore
       */
      if( SCIPisLE(scip, inferlb, SCIPbdchginfoGetOldbound(bdchginfo)) )
         break;
         
      SCIPdebugMessage("***** relaxed upper bound of inference variable <%s> from <%g> to <%g>\n", 
         SCIPvarGetName(var), SCIPbdchginfoGetNewbound(bdchginfo), SCIPbdchginfoGetOldbound(bdchginfo));
      
      bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
      *newub = SCIPbdchginfoGetOldbound(bdchginfo);
      nbdchgs--;
   }
   
   /* if the nbdchgs is zero then the local bound matches the global bound, therefore bdchgidx equal to NULL represents
    * the right time point and SCIP finds out that this bound is redundant since it is global
    */
   SCIPdebugMessage("add upper bound of bound change info %d to conflict set\n", nbdchgs);
   SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx) );
   
   SCIPdebugMessage("relaxed upper bound to %.15g\n", *newub);

   return SCIP_OKAY;
}

/** relax upper bound of give variable as long as the given inference bound leads still to a cutoff and add that bound
 *  change to the conflict set
 */
static
SCIP_RETCODE relaxVbdvarUpperbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable for which the upper bound should be relaxed */
   SCIP_Real             coef,               /**< variable bound coefficient */
   SCIP_Real             constant,           /**< variable bound constant */
   SCIP_Real             bound,              /**< bound to exceed */
   SCIP_Real*            inferbound          /**< pointer to store the relaxed inferbound bound */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_BDCHGIDX* bdchgidx;
   int nbdchgs;

   /* get number of bound changes */
   nbdchgs = SCIPvarGetNBdchgInfosUb(var);
   bdchgidx = NULL;
   
   SCIPdebugMessage("variable <%s>[%.15g,%.15g], coef=%.15g, constant=%.15g: nbdchgs %d try to relax upper bound to %.15g\n", 
      SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), coef, constant, nbdchgs, bound);
   
   assert(nbdchgs > 0 || SCIPisEQ(scip, SCIPvarGetUbLocal(var), SCIPvarGetUbGlobal(var)));
   assert(nbdchgs == 0 || SCIPisEQ(scip, SCIPvarGetUbLocal(var), SCIPbdchginfoGetNewbound(SCIPvarGetBdchgInfoUb(var, nbdchgs-1))));

   /* try to relax upper bound */
   while( nbdchgs > 0 )
   {
      bdchginfo = SCIPvarGetBdchgInfoUb(var, nbdchgs-1);

      SCIPdebugMessage("upper bound change %d oldbd=%.15g, newbd=%.15g, depth=%d, pos=%d, redundant=%u\n", 
         nbdchgs, SCIPbdchginfoGetOldbound(bdchginfo), SCIPbdchginfoGetNewbound(bdchginfo), 
         SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo), SCIPbdchginfoIsRedundant(bdchginfo));
      
      /* check if the old upper bound is sufficient to prove infeasibility; in case the inference bound is greater
       * equal to the next possible relaxed upper bound, then we have to break since in this case the inference bound
       * does not lead to a cutoff anymore
       */
      if( SCIPisGE(scip, bound, coef * SCIPbdchginfoGetOldbound(bdchginfo) + constant) )
         break;
         
      SCIPdebugMessage("***** relaxed upper bound of vbound variable <%s> from <%g> to <%g>\n", 
         SCIPvarGetName(var), SCIPbdchginfoGetNewbound(bdchginfo), SCIPbdchginfoGetOldbound(bdchginfo));

      bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
      *inferbound =  coef * SCIPbdchginfoGetOldbound(bdchginfo) + constant;
      nbdchgs--;
   }
      
   /* if the nbdchgs is zero then the local bound matches the global bound, therefore bdchgidx equal to NULL represents
    * the right time point and SCIP finds out that this bound is redundant since it is global
    */
   SCIPdebugMessage("add upper bound of bound change info %d to conflict set\n", nbdchgs);
   SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx) );

   return SCIP_OKAY;
}

/** relaxes bound of give variable as long as the given inference bound leads still to a cutoff and add that bound
 *  change to the conflict set
 */
static
SCIP_RETCODE relaxVbdvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable for which the upper bound should be relaxed */
   SCIP_BOUNDTYPE        boundtype,          /**< boundtype used for the variable bound variable */
   SCIP_Real             coef,               /**< variable bound coefficient */
   SCIP_Real             constant,           /**< variable bound constant */
   SCIP_Real             bound,              /**< bound to exceed */
   SCIP_Real*            inferbound          /**< pointer to store the relaxed inferbound bound */
   )
{
   if( boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      SCIP_CALL( relaxVbdvarLowerbound(scip, var, coef, constant, bound, inferbound) );
   }
   else
   {
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);
      SCIP_CALL( relaxVbdvarUpperbound(scip, var, coef, constant, bound, inferbound) );
   }

   return SCIP_OKAY;
}


/** analyzes an infeasibility which was reached by updating the lower bound of the inference variable above its upper
 *  bound
 */
static
SCIP_RETCODE analyzeConflictLowerbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */   
   SCIP_VAR*             infervar,           /**< variable which lead to a cutoff */
   SCIP_Real             inferlb,            /**< lower bound which lead to infeasibility */
   INFERINFO             inferinfo,          /**< inference information */
   SCIP_Real             coef,               /**< inference variable bound coefficient used */
   SCIP_Real             constant            /**< inference variable bound constant used */
   )
{
   assert(scip != NULL);
   assert(propdata != NULL);
   assert(infervar != NULL);
   assert(SCIPisEQ(scip, SCIPvarGetUbLocal(infervar), SCIPvarGetUbAtIndex(infervar, NULL, FALSE)));
   assert(SCIPisEQ(scip, SCIPvarGetUbAtIndex(infervar, NULL, TRUE), SCIPvarGetUbAtIndex(infervar, NULL, FALSE)));
   assert(SCIPisGT(scip, inferlb, SCIPvarGetUbLocal(infervar)));
   assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);

   /* check if conflict analysis is applicable */
   if( !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   if( propdata->usebdwidening )
   {
      SCIP_VAR* vbdvar;
      SCIP_Real newub;
      SCIP_Real previnferlb;
      int pos;
      
      pos = inferInfoGetPos(inferinfo);
      assert(pos >= 0);
      assert(pos < propdata->nvars);
      
      vbdvar = propdata->vars[pos];
      newub = SCIPvarGetUbLocal(infervar); 
      previnferlb = inferlb;
      
      SCIPdebugMessage("try to create conflict using bound widening order: inference variable, variable bound variable\n");

      /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
      SCIP_CALL( SCIPinitConflictAnalysis(scip) );

      /* try to relax inference variable upper bound bounds */
      SCIP_CALL( relaxInfervarUpperbound(scip, infervar, inferlb, &newub) );
      
      /* try to relax variable bound variable */
      SCIP_CALL( relaxVbdvar(scip, vbdvar, inferInfoGetBoundtype(inferinfo), coef, constant, newub, &previnferlb) );

      /* analyze the conflict */
      SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );

      /* if the upper bound relaxation was successful we try to create another conflict by relaxing the bound of the
       * variable bound variable first 
       */
      if( SCIPisGT(scip, newub, SCIPvarGetUbLocal(infervar)) ) 
      {
         SCIPdebugMessage("try to create conflict using bound widening order: variable bound variable, inference variable\n");

         /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
         SCIP_CALL( SCIPinitConflictAnalysis(scip) );

         newub = SCIPvarGetUbLocal(infervar); 
         
         /* try to relax variable bound variable */
         SCIP_CALL( relaxVbdvar(scip, vbdvar, inferInfoGetBoundtype(inferinfo), coef, constant, newub, &inferlb) );
         
         /* continue conflict analysis only if we improved the inference lower bound; otherwise we end up with previous
          * conflict set 
          */
         if( SCIPisLT(scip, inferlb, previnferlb ) )
         {
            /* try to relax inference variable upper bound bounds */
            SCIP_CALL( relaxInfervarUpperbound(scip, infervar, inferlb, &newub) );
         
            /* analyze the conflict */
            SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );
         }
      }
   }
   else
   {
      /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
      SCIP_CALL( SCIPinitConflictAnalysis(scip) );

      /* add upper bound of the variable for which we tried to change the lower bound */
      SCIP_CALL( SCIPaddConflictUb(scip, infervar, NULL) );
      
      /* add (correct) bound of the variable which let to the new lower bound */
      SCIP_CALL( resolvePropagation(scip, propdata, inferinfo, NULL) );

      /* analyze the conflict */
      SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );
   }
   
   return SCIP_OKAY;
}

/** analyzes an infeasibility which was reached by updating the upper bound of the inference variable below its lower
 *  bound
 */
static
SCIP_RETCODE analyzeConflictUpperbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */   
   SCIP_VAR*             infervar,           /**< variable which lead to a cutoff */
   SCIP_Real             inferub,            /**< upper bound which lead to infeasibility */
   INFERINFO             inferinfo,          /**< inference information */
   SCIP_Real             coef,               /**< inference variable bound coefficient used */
   SCIP_Real             constant            /**< inference variable bound constant used */
   )
{
   assert(scip != NULL);
   assert(propdata != NULL);
   assert(infervar != NULL);
   assert(SCIPisEQ(scip, SCIPvarGetLbLocal(infervar), SCIPvarGetLbAtIndex(infervar, NULL, FALSE)));
   assert(SCIPisEQ(scip, SCIPvarGetLbAtIndex(infervar, NULL, TRUE), SCIPvarGetLbAtIndex(infervar, NULL, FALSE)));
   assert(SCIPisLT(scip, inferub, SCIPvarGetLbLocal(infervar)));
   assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);

   /* check if conflict analysis is applicable */
   if( !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   if( propdata->usebdwidening )
   {
      SCIP_VAR* vbdvar;
      SCIP_Real newlb;
      SCIP_Real previnferub;
      int pos;
      
      pos = inferInfoGetPos(inferinfo);
      assert(pos >= 0);
      assert(pos < propdata->nvars);
      
      vbdvar = propdata->vars[pos];
      newlb = SCIPvarGetLbLocal(infervar); 
      previnferub = inferub;

      SCIPdebugMessage("try to create conflict using bound widening order: inference variable, variable bound variable\n");

      /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
      SCIP_CALL( SCIPinitConflictAnalysis(scip) );
      
      /* try to relax inference variable lower bound bounds */
      SCIP_CALL( relaxInfervarLowerbound(scip, infervar, inferub, &newlb) );
      
      /* try to relax variable bound variable */
      SCIP_CALL( relaxVbdvar(scip, vbdvar, inferInfoGetBoundtype(inferinfo), -coef, -constant, -newlb, &previnferub) );

      /* analyze the conflict */
      SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );

      /* if the lower bound relaxation was successful we try to create another conflict by relaxing the bound of the
       * variable bound variable first 
       */
      if( SCIPisLT(scip, newlb, SCIPvarGetLbLocal(infervar)) ) 
      {
         SCIPdebugMessage("try to create conflict using bound widening order: variable bound variable, inference variable\n");

         /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
         SCIP_CALL( SCIPinitConflictAnalysis(scip) );

         newlb = SCIPvarGetLbLocal(infervar); 
         
         /* try to relax variable bound variable */
         SCIP_CALL( relaxVbdvar(scip, vbdvar, inferInfoGetBoundtype(inferinfo), -coef, -constant, -newlb, &inferub) );

         /* continue conflict analysis only if we improved the inference upper bound w.r.t. the previous conflict
          * analysis run; otherwise we end up with previous conflict set
          */
         if( SCIPisGT(scip, inferub, previnferub ) )
         {
            /* try to relax inference variable upper bound bounds */
            SCIP_CALL( relaxInfervarUpperbound(scip, infervar, inferub, &newlb) );

            /* analyze the conflict */
            SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );
         }
      }

   }
   else
   {
      /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
      SCIP_CALL( SCIPinitConflictAnalysis(scip) );
      
      /* add lower bound of the variable for which we tried to change the upper bound */
      SCIP_CALL( SCIPaddConflictLb(scip, infervar, NULL) );
      
      /* add (correct) bound of the variable which let to the new upper  bound */
      SCIP_CALL( resolvePropagation(scip, propdata, inferinfo, NULL) );

      /* analyze the conflict */
      SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );
   }
   
   return SCIP_OKAY;
}

/** find position of the given variable in the variable array; if it does not exist yet it gets added to the end of the
 *  array 
 */
static
SCIP_RETCODE getVarPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< vbounds propagator data */
   SCIP_VAR*             var,                /**< variable */
   int*                  pos                 /**< pointer to store position in array */
   )
{
   assert(scip != NULL);
   assert(propdata != NULL);
   assert(var != NULL);
   assert(pos != NULL);

   /* get position of vbvar in variable arrays */
   if( SCIPhashmapExists(propdata->varHashmap, var) )
      *pos = (int)(size_t)SCIPhashmapGetImage(propdata->varHashmap, var);
   else
   {
      /* ensure array size */
      if( propdata->sizevars <= propdata->nvars )
      {
         propdata->sizevars = SCIPcalcMemGrowSize(scip, propdata->nvars + 1);
         SCIP_CALL( SCIPreallocMemoryArray(scip, &(propdata->vars), propdata->sizevars) );
      }
      assert(propdata->sizevars > propdata->nvars);

      propdata->vars[propdata->nvars] = var;
      *pos = propdata->nvars;
      propdata->nvars++;

      /* capture variable to ensure the existence */
      SCIP_CALL( SCIPcaptureVar(scip, var) );
      
      /* insert variable bound variable into the hash table since they are involve in propagation */ 
      SCIP_CALL( SCIPhashmapInsert(propdata->varHashmap, var, (void*)(size_t)*pos) );
   }
   
   return SCIP_OKAY;
}

/** performs propagation of variables lower and upper bounds */
static
SCIP_RETCODE propagateVbounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< vbounds propagator */
   SCIP_Bool             force,              /**< should domain changes be forced */
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
   SCIP_Real bestcoef;
   SCIP_Real bestconstant;
   SCIP_Real newbound;
   INFERINFO inferinfo;
   int nvars;
   int nvbvars;
   int pos;
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
   
   if( propdata->lbpropagated && propdata->ubpropagated )
      return SCIP_OKAY;

   nchgbds = 0;
   nvars = propdata->nlbvars;
   
   if( nvars > 0 && !propdata->lbpropagated )
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
         bestcoef = 1.0;
         bestconstant = 0.0;
         inferinfo = getInferInfo(v, SCIP_BOUNDTYPE_UPPER);

         SCIPdebugMessage("try to improve lower bound of variable <%s> (current loc=[%.15g,%.15g])\n",
            SCIPvarGetName(var), newbound, SCIPvarGetUbLocal(var));
         
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
               SCIP_Real candbound;
               SCIP_Real lb;

               lb =  SCIPvarGetLbLocal(vbvar);

               /* ignore variable bound variables with a lower bound of minus infinity */
               if( SCIPisInfinity(scip, -lb) )
                  continue;

               /* compute candidate bound; if b > 0 => x >= b*lb(y) + d */ 
               candbound =  coef * lb + constant;

               /* check if candidate bound is better */
               if( SCIPisGT(scip, candbound, newbound) )
               {
                  assert(SCIPvarGetProbindex(vbvar) > -1);

                  newbound = candbound;
                  bestcoef = coef;
                  bestconstant = constant;
               
                  SCIPdebugMessage(" -> new lower bound candidate <%.15g> due to lower bound of variable <%s> (n=%d)\n",
                     newbound, SCIPvarGetName(vbvar), n);
                  SCIPdebugMessage("         newlb >= %.15g * [%.15g,%.15g] + %.15g\n", 
                     coef, SCIPvarGetLbLocal(vbvar), SCIPvarGetUbLocal(vbvar), constant);

                  /* get position of vbvar in variable arrays */
                  SCIP_CALL( getVarPos(scip, propdata, vbvar, &pos) );

                  /* construct infer info */
                  inferinfo = getInferInfo(pos, SCIP_BOUNDTYPE_LOWER);
               }
            }
            else
            {
               SCIP_Real candbound;
               SCIP_Real ub;

               ub =  SCIPvarGetUbLocal(vbvar);

               /* ignore variable bound variables with an upper bound of infinity */
               if( SCIPisInfinity(scip, ub) )
                  continue;

               /* compute candidate bound; if b < 0 => x >= b*ub(y) + d */ 
               candbound =  coef * ub + constant;

               /* check if candidate bound is better */
               if( SCIPisGT(scip, candbound, newbound) )
               {
                  assert(SCIPvarGetProbindex(vbvar) > -1);

                  newbound = candbound;
                  bestcoef = coef;
                  bestconstant = constant;

                  SCIPdebugMessage(" -> new lower bound candidate <%.15g> due to upper bound of variable <%s> (n=%d)\n",
                     newbound, SCIPvarGetName(vbvar), n);
                  SCIPdebugMessage("         newlb >= %.15g * [%.15g,%.15g] + %.15g\n", 
                     coef, SCIPvarGetLbLocal(vbvar), SCIPvarGetUbLocal(vbvar), constant);
                  
                  /* get position of vbvar in variable arrays */
                  SCIP_CALL( getVarPos(scip, propdata, vbvar, &pos) );

                  /* construct infer info */
                  inferinfo = getInferInfo(pos, SCIP_BOUNDTYPE_UPPER);
               }
            }
         }
         
         /* try to tighten lower bound */
         SCIP_CALL( SCIPinferVarLbProp(scip, var, newbound, prop, inferInfoToInt(inferinfo), force, &infeasible, &tightened) );
         
         if( infeasible )
         {
            /* the infeasible results comes from the fact that the new lower bound lies above the current upper bound */
            assert(SCIPisGT(scip, newbound, SCIPvarGetUbLocal(var)));
               
            SCIPdebugMessage(" -> variable <%s> => variable <%s> lower bound candidate is <%.15g>\n", 
               SCIPvarGetName(propdata->vars[inferInfoGetPos(inferinfo)]), SCIPvarGetName(var), newbound);
            
            SCIPdebugMessage(" -> lower bound tightening lead to infeasibility\n");
            
            /* analyzes an infeasibility via conflict analysis */
            SCIP_CALL( analyzeConflictLowerbound(scip, propdata, var, newbound, inferinfo, bestcoef, bestconstant) );
            *result = SCIP_CUTOFF;

            return SCIP_OKAY;      
         } 
         
         if( tightened )
         {
            SCIPdebugMessage(" -> tightened lower bound to <%g> due the %s bound of variable <%s>\n", 
               newbound, inferInfoGetBoundtype(inferinfo) == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", 
               SCIPvarGetName(propdata->vars[inferInfoGetPos(inferinfo)]));
            nchgbds++;
         }
      }
      
      /* mark lower bound variable array as propagated */
      propdata->lbpropagated = TRUE;
   }
   
   nvars = propdata->nubvars;

   if( nvars > 0 && !propdata->ubpropagated)
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
         bestcoef = 1.0;
         bestconstant = 0.0;
         inferinfo = getInferInfo(v, SCIP_BOUNDTYPE_UPPER);

         SCIPdebugMessage("try to improve upper bound of variable <%s> (current loc=[%.15g,%.15g])\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), newbound);
         
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
               SCIP_Real candbound;
               SCIP_Real ub;

               ub = SCIPvarGetUbLocal(vbvar);

               /* ignore variable bound variables with an upper bound of infinity */
               if( SCIPisInfinity(scip, ub) )
                  continue;
               
               /* compute candidate for new bound; if b > 0 => x <= b*ub(y) + d */ 
               candbound = coef * ub + constant;
               
               /* check if the candidate is better */
               if(  SCIPisLT(scip, candbound, newbound) )
               {
                  assert(SCIPvarGetProbindex(vbvar) > -1);

                  newbound = candbound;
                  bestcoef = coef;
                  bestconstant = constant;

                  SCIPdebugMessage(" -> new upper bound candidate <%.15g> due to upper bound of variable <%s> (n=%d)\n",
                     newbound, SCIPvarGetName(vbvar), n);
                  SCIPdebugMessage("         newub <= %.15g * [%.15g,%.15g] + %.15g\n", 
                     coef, SCIPvarGetLbLocal(vbvar), SCIPvarGetUbLocal(vbvar), constant);
                  
                  /* get position of vbvar in variable arrays */
                  SCIP_CALL( getVarPos(scip, propdata, vbvar, &pos) );

                  /* construct infer info */
                  inferinfo = getInferInfo(pos, SCIP_BOUNDTYPE_UPPER);
               }
            }
            else
            {
               SCIP_Real candbound;
               SCIP_Real lb;

               lb = SCIPvarGetLbLocal(vbvar);

               /* ignore variable bound variables with a lower bound of minus infinity */
               if( SCIPisInfinity(scip, -lb) )
                  continue;
               
               /* compute candidate bound; if b < 0 => x <= b*lb(y) + d */ 
               candbound = coef * lb + constant;

               /* check if candidate bound is better */
               if( SCIPisLT(scip, candbound, newbound) )
               {
                  assert(SCIPvarGetProbindex(vbvar) > -1);

                  newbound = candbound;
                  bestcoef = coef;
                  bestconstant = constant;
                  
                  SCIPdebugMessage(" -> new upper bound candidate <%.15g> due to lower bound of variable <%s> (n=%d)\n",
                     newbound, SCIPvarGetName(vbvar), n);
                  SCIPdebugMessage("         newub <= %.15g * [%.15g,%.15g] + %.15g\n", 
                     coef, SCIPvarGetLbLocal(vbvar), SCIPvarGetUbLocal(vbvar), constant);
                  
                  /* get position of vbvar in variable arrays */
                  SCIP_CALL( getVarPos(scip, propdata, vbvar, &pos) );

                  /* construct infer info */
                  inferinfo = getInferInfo(pos, SCIP_BOUNDTYPE_LOWER);
               }
            }
         }
      
         /* try to tighten upper bound */
         SCIP_CALL( SCIPinferVarUbProp(scip, var, newbound, prop, inferInfoToInt(inferinfo), force, &infeasible, &tightened) );
      
         if( infeasible )
         {
            /* the infeasible results from the fact that the new upper bound lies below the current lower bound */
            assert(SCIPisLT(scip, newbound, SCIPvarGetLbLocal(var)));

            SCIPdebugMessage(" -> variable <%s> => variable <%s> upper bound candidate is <%.15g>\n", 
               SCIPvarGetName(propdata->vars[inferInfoGetPos(inferinfo)]), SCIPvarGetName(var), newbound);

            SCIPdebugMessage(" -> upper bound tightening lead to infeasibility\n");
            
            /* analyzes an infeasibility via conflict analysis */
            SCIP_CALL( analyzeConflictUpperbound(scip, propdata, var, newbound, inferinfo, bestcoef, bestconstant) );
            *result = SCIP_CUTOFF;
            
            return SCIP_OKAY;      
         } 

         if( tightened )
         {
            SCIPdebugMessage(" -> tightened upper bound to <%g> due the %s bound of variable <%s>\n", 
               newbound, inferInfoGetBoundtype(inferinfo) == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", 
               SCIPvarGetName(propdata->vars[inferInfoGetPos(inferinfo)]));
            nchgbds++;
         }
      }
      
      /* mark upper bound variable array as propagated */
      propdata->ubpropagated = TRUE;
   }   

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

/** copy method for propagator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PROPCOPY(propCopyVbounds)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* call inclusion method of propagator */
   SCIP_CALL( SCIPincludePropVbounds(scip) );
 
   return SCIP_OKAY;
}

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


/** presolving initialization method of propagator (called when presolving is about to begin) */
#define propInitpreVbounds NULL


/** presolving deinitialization method of propagator (called after presolving has been finished) */
#define propExitpreVbounds NULL


/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
static
SCIP_DECL_PROPINITSOL(propInitsolVbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_VAR** vars;
   int nvars;
   int v;

   SCIPdebugMessage("initialize vbounds propagator for problem <%s>\n", SCIPgetProbName(scip));

   /* free propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   
   if( nvars == 0 )
      return SCIP_OKAY;

   /* allocate memory for the arrays of the propdata */
   SCIP_CALL( SCIPallocMemoryArray(scip, &propdata->vars, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &propdata->lbvars, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &propdata->ubvars, nvars) );

   /* store size of the variable array */
   propdata->sizevars = nvars;
   
   /* create hash table for storing the involved variables */
   assert(propdata->nvars == 0);
   SCIP_CALL( SCIPhashmapCreate(&propdata->varHashmap, SCIPblkmem(scip), SCIPcalcHashtableSize(5 * nvars)) );
   
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

   /* catch variable events */
   SCIP_CALL( catchEvents(scip, propdata) );

   if( propdata->nlbvars > 0 )
      propdata->lbpropagated = FALSE;
   
   if( propdata->nubvars > 0 )
      propdata->ubpropagated = FALSE;

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

   /* free hash map */
   if( propdata->varHashmap != NULL )
      SCIPhashmapFree(&propdata->varHashmap);
   
   /* free array */
   SCIPfreeMemoryArrayNull(scip, &propdata->lbeventtypes);
   SCIPfreeMemoryArrayNull(scip, &propdata->ubeventtypes);
   SCIPfreeMemoryArrayNull(scip, &propdata->lbvars);
   SCIPfreeMemoryArrayNull(scip, &propdata->ubvars);
   SCIPfreeMemoryArrayNull(scip, &propdata->vars);

   /* reset propagation data */
   resetPropdata(propdata);

   return SCIP_OKAY;
}


/** presolving method of propagator */
#define propPresolVbounds NULL


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecVbounds)
{  /*lint --e{715}*/
   
   /* perform variable lower and upper bound propagation */
   SCIP_CALL( propagateVbounds(scip, prop, FALSE, result) );
   
   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropVbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   
   SCIPdebugMessage("explain %s bound change of variable <%s>\n", 
      boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", SCIPvarGetName(infervar));
   
   SCIP_CALL( resolvePropagation(scip, propdata, intToInferInfo(inferinfo), bdchgidx) );

   (*result) = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/*
 * Event Handler
 */

/** execution method of bound change event handler */
static
SCIP_DECL_EVENTEXEC(eventExecVbound)
{  /*lint --e{715}*/
   SCIP_Bool* propagated;

   propagated = (SCIP_Bool*)eventdata;
   assert(propagated != NULL);

   (*propagated) = FALSE;

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
   
   /*  reset propagation data */
   resetPropdata(propdata);

   /* include propagator */
   SCIP_CALL( SCIPincludeProp(scip, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS, PROP_PRESOL_DELAY,
         propCopyVbounds,
         propFreeVbounds, propInitVbounds, propExitVbounds, propInitpreVbounds, propExitpreVbounds, 
         propInitsolVbounds, propExitsolVbounds, propPresolVbounds, propExecVbounds, propRespropVbounds,
         propdata) );

   /* include event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, eventExecVbound, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/"PROP_NAME"/usebdwidening", "should bound widening be used to initialize conflict analysis?",
         &propdata->usebdwidening, FALSE, DEFAULT_USEBDWIDENING, NULL, NULL) );
   
   
   return SCIP_OKAY;
}

/** create a topological sorted variable array of the given variables and stores if (needed) the involved variables into
 *  the corresponding variable array and hash map
 *
 * @note: for all arrays and the hash map (if requested) you need to allocate enough memory before calling this method 
 */
SCIP_RETCODE SCIPcreateTopoSortedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variable which we want sort */
   int                   nvars,              /**< number of variables */
   SCIP_HASHMAP*         varPosMap,          /**< mapping a variable to its position in the (used) variable array, or NULL */    
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
   int hashsize;
   int i;
   int v;
   
   assert(scip != NULL);   
   assert(vars != NULL || nvars == 0);
   assert(varPosMap == NULL || (varPosMap != NULL && usedvars != NULL && nusedvars != NULL));
   assert(topovars != NULL);
   assert(ntopovars != NULL);
   
   SCIPdebugMessage("create topological sorted variable array with respect to variables %s bounds\n", 
      lowerbound ? "lower" : "upper");

   if( nvars == 0 )
      return SCIP_OKAY;

   assert(vars != NULL);
   
   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &sortedvars, nvars) );
   
   hashsize = SCIPcalcHashtableSize(5 * nvars);

   /* create hash table for variables which are (still) connected */
   SCIP_CALL( SCIPhashtableCreate(&connected, SCIPblkmem(scip), hashsize, SCIPvarGetHashkey, SCIPvarIsHashkeyEq, SCIPvarGetHashkeyVal, NULL) );
   
   /* detect isolated variables; mark all variables which have at least one entering or leaving arc as connected */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert(var != NULL);
      
      if( !SCIPvarIsActive(var) )
         continue;

      /* get variable bounds */
      getVariableBounds(var, &vbvars, &nvbvars, lowerbound);
      
      if( nvbvars > 0 && !SCIPhashtableExists(connected, var) )
      {
         SCIP_CALL( SCIPhashtableInsert(connected, var) );
      }

      for( i = 0; i < nvbvars; ++i )
      {
         if( !SCIPvarIsActive(vbvars[i]) )
            continue;

         /* there is a leaving arc, hence, the variable/node  is connected */  
         assert(vbvars[i] != NULL);
         if( !SCIPhashtableExists(connected, vbvars[i]) )
         {
            SCIP_CALL( SCIPhashtableInsert(connected, vbvars[i]) );
         }
      }
   }

   /* loop over all "connected" variable and find for each connected component a "almost" topological sorted version */
   for( v = 0; v < nvars; ++v )
   {
      if( SCIPhashtableExists(connected, vars[v]) )
      {
         int nsortedvars;

         SCIPdebugMessage("start depth-first-search with variable <%s>\n", SCIPvarGetName(vars[v]));
         
         /* use depth first search to get a "almost" topological sorted variables for the connected component which
          * includes vars[v]
          */
         nsortedvars = 0;
         SCIP_CALL( depthFirstSearch(scip, vars[v], varPosMap, usedvars, nusedvars, connected, sortedvars, &nsortedvars, lowerbound) );
         
         SCIPdebugMessage("detected connected component of size <%d>\n", nsortedvars);
        
         /* copy variables */
         for( i = 0; i < nsortedvars; ++i )
         {
            topovars[(*ntopovars)] = sortedvars[i];
            (*ntopovars)++;
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

/** returns TRUE if the propagator has the status that all variable lower and upper bounds are propagated */
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

   return propdata->lbpropagated && propdata->ubpropagated;
}

/** performs propagation of variables lower and upper bounds */
SCIP_RETCODE SCIPexecPropVbounds(
   SCIP*                 scip,                /**< SCIP data structure */
   SCIP_Bool             force,               /**< should domain changes be forced */
   SCIP_RESULT*          result               /**< pointer to store result */
   )
{
   SCIP_PROP* prop;
   
   prop = SCIPfindProp(scip, PROP_NAME);
   assert(prop != NULL);

   /* perform variable lower and upper bound propagation */
   SCIP_CALL( propagateVbounds(scip, prop, force, result) );

   return SCIP_OKAY;
}

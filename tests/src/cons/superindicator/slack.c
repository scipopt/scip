/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   initlp.c
 * @brief  unit test for checking behaviour of the initlp callback
 *
 * TODO: how to specify the path of the file?
 * What should the test actually test?
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"



#define CONSHDLR_NAME                        "superindicator"

/** type of bound relaxation in SCIPsolveSlack() method */
enum SCIP_SlackType
{
   SCIP_SLACKTYPE_LOWER = 0,                 /**< lower bound */
   SCIP_SLACKTYPE_UPPER = 1,                 /**< upper bound */
   SCIP_SLACKTYPE_BOTH  = 2                  /**< both bounds */
};
typedef enum SCIP_SlackType SCIP_SLACKTYPE;

/** find the position of a variable in an array of variables; returns -1 if not found */
static
int findvarpos(
   SCIP_VAR*             var,                /**< variable to find */
   SCIP_VAR**            vars,               /**< array of variables to be searched */
   int                   nvars               /**< size of the array to be searched */
   )
{
   int c;

   assert(var != NULL);
   assert(vars != NULL);

   for( c = nvars-1; c >= 0; c-- )
   {
      if( vars[c] == var )
         return c;
   }

   return -1;
}

/** copies an array of variables and objective values */
static
SCIP_RETCODE copyVarsAndObjVals(
   SCIP_VAR**            vars,               /**< source array of variables */
   SCIP_VAR**            copiedvars,         /**< target array of variables */
   SCIP_Real*            objvals,            /**< target array of the obj. values */
   int                   nvars               /**< number of variables to copy */
   )
{
   int c;

   assert(vars != NULL);
   assert(copiedvars != NULL);
   assert(objvals != NULL);

   for( c = nvars-1; c >= 0; c-- )
   {
      copiedvars[c] = vars[c];
      assert(copiedvars[c] != NULL);

      objvals[c] = SCIPvarGetObj(copiedvars[c]);
   }

   return SCIP_OKAY;
}

/** solve a slack model minizing the violation of constraints and/or variable bounds */
/**@todo allocating the memory in the correct way */
static
SCIP_RETCODE SCIPsolveSlack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS***          consgroups,         /**< array of constraint groups */
   int                   nconsgroups,        /**< number of constraint groups */
   int*                  consgroupsizes,     /**< array of sizes of constraint groups (may be NULL if all constraint
                                              *   groups contain exactly one constraint) */
   int*                  consgroupprios,     /**< array of priorities of constraint groups */
   SCIP_Real*            consgroupweights,   /**< array of weights of for the slack of each constraint group (may be
                                              *   NULL for unit weights) */
   SCIP_VAR***           vargroups,          /**< array of variable groups on which the bounds has to be relaxed */
   int                   nvargroups,         /**< number of variable groups */
   int*                  vargroupsizes,      /**< array of sizes of variable groups on which the bounds has to be relaxed
                                              *   (may be NULL if all constraintgroups contain exactly one constraint) */
   SCIP_SLACKTYPE**      varslacktypes,      /**< array of variable bound slack types (may be NULL if both bounds are
                                              *   relaxed) */
   int*                  vargroupprios,      /**< array of priorities of variable groups */
   SCIP_Real*            vargroupweights,    /**< array of weights of for the slack of each variable bounds group (may be
                                              *   NULL for unit weights) */
   SCIP_Real*            consgroupsslacks,   /**< buffer array for minimal constraint group relaxations ready to hold
                                              *   nconsgroups + nvargroups elements (has to be allocated for nconsgroups + nvargroups
                                              *   reals or may be NULL. the order correspond first to the array congroups
                                              *   and then the array vargroups) */
   SCIP_SOL**            sol,                /**< solution feasible for the problem relaxed according to
                                              *   consgroupsslacks (may be NULL) */
   SCIP_Real*            consviols,          /**< buffer array for violations of each constraint in the returned
                                              *   solution (has to be allocated for the number of constraints
                                              *   and variables to be relaxed or may be NULL.  will return first
                                              *   the violation for the constraints and then the one for the relaxed bounds) */
   SCIP_Bool             optorig,            /**< shall solution be optimal w.r.t. the original objective function? */
   SCIP_STATUS*          status              /**< pointer to store the solving status of the last problem solved */
   )
{
   SCIP_CONS*** supindconsgroups;
   SCIP_CONS*** supindvarboundgroups;
   SCIP_CONS*** copiedconsgroups;
   SCIP_CONS** varboundconss;
   SCIP_SOL** sols;
   SCIP_VAR*** copiedvargroups;
   SCIP_VAR** origvars;
   SCIP_VAR** stockbinvars;
   SCIP_VAR** allbinvars;
   SCIP_VAR** vars;
   SCIP_Real** solsvals;
   SCIP_Real* origobjvals;
   SCIP_Real* stockweights;
   SCIP_Real previousobjval;
   SCIP_Bool vargroupisempty;
   SCIP_Bool consgroupisempty;
   SCIP_Bool origprobmin;
   SCIP_Bool lastround;
   int* copiedconsgroupprios;
   int* copiedvargroupprios;
   int* consgroupindices;
   int* vargroupindices;
   int nsols;
   int nvars;
   int norigvars;
   int maxprio;
   int currprio;
   int currconsgroup;
   int currvarboundgroup;
   int c;
   int count;
#ifndef NDEBUG
   int minprio;
#endif

   assert(scip != NULL);
   assert(status != NULL);
   assert(nconsgroups >= 0);
   assert(nvargroups >= 0);

   copiedvargroups = NULL;
   varboundconss = NULL;
   supindvarboundgroups = NULL;
   copiedconsgroups = NULL;
   supindconsgroups = NULL;
   copiedvargroupprios = NULL;
   vargroupindices = NULL;
   consgroupindices = NULL;
   copiedconsgroupprios = NULL;

   /* testing whether we have correct set to work on */
   if( SCIPgetStage(scip) != SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("method <SCIPsolveSlack> can only be called in problem stage\n");
      return SCIP_INVALIDCALL;
   }

   if( nconsgroups + nvargroups == 0 )
   {
      SCIPdebugMessage("no relaxed bounds or constraints: solving the original problem\n");
      SCIP_CALL( SCIPsolve(scip) );
      *status = SCIPgetStatus(scip);
      return SCIP_OKAY;
   }

   if( nconsgroups == 0)
   {
      assert(consgroups == NULL);
      assert(consgroupprios == NULL);
      assert(consgroupweights == NULL);
      assert(consgroupsizes == NULL);

      consgroupisempty = TRUE;
   }
   else
   {
      assert(consgroups != NULL);
      assert(consgroupprios != NULL);
      assert(consgroupweights != NULL);
      assert(consgroupsizes != NULL);

      consgroupisempty = FALSE;
   }

   if( nvargroups == 0)
   {
      assert(vargroups == NULL);
      assert(vargroupprios == NULL);
      assert(vargroupweights == NULL);
      assert(vargroupsizes == NULL);

      vargroupisempty = TRUE;
   }
   else
   {
      assert(vargroups != NULL);
      assert(vargroupprios != NULL);
      assert(vargroupweights != NULL);
      assert(vargroupsizes != NULL);

      vargroupisempty = FALSE;
   }

   /* The function should not be called if we don´t want to relax any bound or constraint */
   if(vargroupisempty && consgroupisempty)
   {
      SCIPerrorMessage("method <SCIPsolveSlack> has to be given at least one group of constraint to be relaxed.\n");
      return SCIP_INVALIDCALL;
   }

   /* check whether we have a binary variable on which we want to relax the bounds */
   for(c = 0; c < nvargroups; ++c)
   {
      int d;
      for(d = 0; d < vargroupsizes[c]; d++)
      {
         if( SCIP_VARTYPE_BINARY == SCIPvarGetType(vargroups[c][d]) )
         {
            SCIPerrorMessage("method <SCIPsolveSlack> cannot relaxe binary variables.\n");
            return SCIP_INVALIDCALL;
         }
      }
   }

   /* trick to erase the array of original solution at the end of each solving round *
    *@todo if moved to scip.c need to use the SCIPchgIntParam() instead of the current function */
   SCIP_CALL( SCIPsetIntParam(scip, "limits/maxorigsol", 0) );

   /* checking the objective sense of the original problem and changing it when necessary */
   if( SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE )
      origprobmin = TRUE;
   else
   {
      SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );
      origprobmin = FALSE;
   }

   /* making a copy of the objective value of each variable of the original problem with a copy of their corresponding variables */
   norigvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocMemoryArray(scip, &origobjvals, norigvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &origvars, norigvars) );
   assert(origvars != NULL);
   assert(norigvars != 0);

   /* copying the original objective function to solve the original problem */
   if(optorig)
   {
      SCIP_CALL( copyVarsAndObjVals(SCIPgetVars(scip), origvars, origobjvals, norigvars) );
   }

   /* copying the priority arrays to work on it, creating the indices arrays and sorting them in a decreasing order of their priority */
   if( !consgroupisempty )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &copiedconsgroupprios, consgroupprios, nconsgroups) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &consgroupindices, nconsgroups) );

      for( c = 0; c < nconsgroups; ++c)
      {
         assert(copiedconsgroupprios[c] == consgroupprios[c]);
         consgroupindices[c] = c;
      }

      SCIPsortDownIntInt(copiedconsgroupprios, consgroupindices, nconsgroups);
   }

   if( !vargroupisempty )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &copiedvargroupprios, vargroupprios, nvargroups) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &vargroupindices, nvargroups) );

      for( c = 0; c < nvargroups; ++c)
      {
         assert(copiedvargroupprios[c] == vargroupprios[c]);
         vargroupindices[c]=c;
      }

      SCIPsortDownIntInt(copiedvargroupprios, vargroupindices, nvargroups);
   }

   if( vargroupisempty || ( !consgroupisempty && (copiedconsgroupprios[0] > copiedvargroupprios[0] )))
      maxprio = copiedconsgroupprios[0];
   else
      maxprio = copiedvargroupprios[0];

#ifndef NDEBUG
   if( vargroupisempty || ( !consgroupisempty && copiedconsgroupprios[nconsgroups-1] < copiedvargroupprios[nvargroups-1] ))
      minprio = copiedconsgroupprios[nconsgroups-1];
   else
      minprio = copiedvargroupprios[nvargroups-1];
#endif

   /* remember all the binvars that will be created to return value at the end of the function */
   SCIP_CALL( SCIPallocMemoryArray(scip, &allbinvars, nconsgroups + nvargroups) );

   /* create the superindicator constraints for the member of the consgroups array */
   count = 0;
   if( !consgroupisempty )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &copiedconsgroups, consgroups, nconsgroups) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(supindconsgroups), nconsgroups) );

      for( c = 0; c < nconsgroups; ++c)
      {
         SCIP_VAR* binvar;
         SCIP_VAR* negbinvar;
         char name[SCIP_MAXSTRLEN];
         int d;

         if( consgroupsizes[c] == 0 )
         {
            SCIPdebugMessage("group of constraint #%i is of size 0\n",c+1);
            continue;
         }

         /* allocating the space for the group pointer */
         SCIP_CALL( SCIPallocMemoryArray(scip, &(supindconsgroups[c]), consgroupsizes[c]) );

         /* creating the binary variable with its objective corresponding to the assigned weight and name cooresponding
          * to the one of the first constraint of the group
          * NOTE: adding the constraint at that point already does not creat any problem if we set it the objective value 0 */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_bin", SCIPconsGetName(copiedconsgroups[c][0]) );
         SCIP_CALL( SCIPcreateVar(scip, &binvar, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
               TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPchgVarBranchPriority(scip, binvar, 1) );
         SCIP_CALL( SCIPaddVar(scip, binvar) );

         /* NOTE: we are stocking the binary variables in the same order as given in argument so that the returned values
          * correspond to the array of group of constraint given in argument
          */
         allbinvars[count] = binvar;
         ++count;

         /* get negated variable, since we want to minimize the number of violated constraints */
         SCIP_CALL( SCIPgetNegatedVar(scip, binvar, &negbinvar) );

         for(d = 0; d < consgroupsizes[c]; ++d)
         {
            SCIP_CONS* cons;

            cons = copiedconsgroups[c][d];
            assert(cons != NULL);
            SCIP_CALL( SCIPcaptureCons(scip, cons) );

            /* avoid to construct a superindicator constraint of a superindicator constraint */
            if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0 )
            {
               SCIPerrorMessage("trying to construct a superindicator constraint of a superindicator constraint<%s>\n",SCIPconsGetName(cons));
               return SCIP_INVALIDCALL;
            }

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_supind", SCIPconsGetName(cons) );

            /* constructing the superindicator constraint and capturing (hence have to release it latter) */
            SCIP_CALL( SCIPcreateConsSuperindicator(scip, &supindconsgroups[c][d], name, negbinvar, cons,
                  FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
            SCIP_CALL( SCIPcaptureCons(scip, supindconsgroups[c][d]) );

            /* don´t realease the constraint  as it is not yet added to scip and we will use it later */
            assert(supindconsgroups[c][d] != 0);

            /* remove the constraint which is now replaced by a superindicator constraint */
            SCIP_CALL( SCIPdelCons(scip, cons) );

         }

         /* release the binary variable */
         SCIP_CALL( SCIPreleaseVar(scip, &binvar) );
      }
   }

   /* create the superindicator constraint for the relaxed variable bounds */
   if( !vargroupisempty )
   {
      int countvar;

      /* prepare array for later compilation of return value */
      if(consviols != NULL)
      {
         countvar = 0;

         for(c = 0; c < nvargroups; ++c)
         {
            countvar += vargroupsizes[c];
         }
         SCIP_CALL( SCIPallocMemoryArray(scip, &varboundconss, countvar) );
      }

      /* duplicate the var groups to avoid any change in the original var groups */
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &copiedvargroups, vargroups, nvargroups) );

      /* allocating the memory for the superindvargroups pointer */
      SCIP_CALL( SCIPallocMemoryArray(scip, &(supindvarboundgroups), nvargroups) );

      countvar = 0;
      /* create the superindicator constraints for the variable bounds present in the vargroups array */
      for( c = 0; c < nvargroups; ++c)
      {
         SCIP_VAR* binvar;
         SCIP_VAR* negbinvar;
         SCIP_Real one;
         char name[SCIP_MAXSTRLEN];
         int d;

         assert(vargroupsizes[c] > 0);
         one = 1;

         /* allocating the space for the group pointer */
         SCIP_CALL( SCIPallocMemoryArray(scip, &(supindvarboundgroups[c]), vargroupsizes[c]) );

         /* creating the binary variable with its objective corresponding to the assigned weight and name cooresponding
          * to the one of the first constraint of the group
          * NOTE: adding the bin var constraint at that point already does not creat any problem if we set its objective value to 0 */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_bound_bin", SCIPvarGetName(copiedvargroups[c][0]) );
         SCIP_CALL( SCIPcreateVar(scip, &binvar, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
               TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPchgVarBranchPriority(scip, binvar, 1) );

         /* NOTE: it might be good to add it only later */
         SCIP_CALL( SCIPaddVar(scip, binvar) );
         allbinvars[count] = binvar;
         ++count;

         /* get negated variable, since we want to minimize the number of violated constraints */
         SCIP_CALL( SCIPgetNegatedVar(scip, binvar, &negbinvar) );

         for(d = 0; d < vargroupsizes[c]; ++d)
         {
            SCIP_CONS* cons;
            SCIP_Real lbound;
            SCIP_Real ubound;

            ubound = SCIPvarGetUbOriginal(vargroups[c][d]);
            lbound = SCIPvarGetLbOriginal(vargroups[c][d]);

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_bound_cons", SCIPvarGetName(copiedvargroups[c][0]) );

            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 1, &(vargroups[c][d]), &one, lbound, ubound,
                  FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
            assert(cons != NULL);
            SCIP_CALL( SCIPcaptureCons(scip, cons) );

            if(consviols != NULL)
            {
               varboundconss[countvar] = cons;
               ++countvar;
            }

            switch ( varslacktypes[c][d] )
            {
            case SCIP_SLACKTYPE_LOWER:
               SCIPchgVarLb(scip, vargroups[c][d], -SCIPinfinity(scip));
               break;

            case SCIP_SLACKTYPE_UPPER:
               SCIPchgVarUb(scip, vargroups[c][d], SCIPinfinity(scip));
               break;

            case SCIP_SLACKTYPE_BOTH:
               SCIPchgVarUb(scip, vargroups[c][d], SCIPinfinity(scip));
               SCIPchgVarLb(scip, vargroups[c][d], -SCIPinfinity(scip));
               break;

            default:
               SCIPerrorMessage("invalid slack type\n");
               return SCIP_INVALIDDATA;
            }

            /* constructing the superindicator constraint (releasing the other constraint*/
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_bound_cons_supind", SCIPvarGetName(copiedvargroups[c][0]) );

            SCIP_CALL ( SCIPcreateConsSuperindicator(scip, &supindvarboundgroups[c][d], name, negbinvar, cons,
                  FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

            SCIP_CALL( SCIPcaptureCons(scip, supindvarboundgroups[c][d]) );

            /* remove the constraint which is now replaced by a superindicator constraint */
            SCIP_CALL( SCIPdelCons(scip, cons) );

         }

         /* releasing the variables will add the binvar later */
         SCIP_CALL( SCIPreleaseVar(scip, &binvar) );
      }
   }

   /* setting the initial values for the indices */
   currprio = maxprio;
   currconsgroup = 0;
   currvarboundgroup = 0;
   lastround = FALSE;
   nsols = 0;

   /* initializing the pointers correctly for sols ans sols vals */
   sols = NULL;
   solsvals = NULL;

   SCIP_CALL( SCIPallocMemoryArray(scip, &stockbinvars, nconsgroups + nvargroups) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &stockweights, nconsgroups + nvargroups) );

   /* resetting all objective value to 0 for all variables */
   /* NOTE : doing that at the beginning only so that at each step we have the golbal objective value displayed */
   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);
   for( c = 0; c < nvars; ++c)
   {
      SCIP_VAR* var;

      var =  vars[c];
      assert(var != NULL);

      SCIP_CALL( SCIPchgVarObj(scip, var, 0) );
   }

   previousobjval = 0;
   /* solving the problem step by step */
   do {
      SCIP_CONS* cons;
      SCIP_Real objval;
      char name[SCIP_MAXSTRLEN];

      SCIPdebugMessage("processing groups of priority <%i>\n", currprio);

      assert(currprio >= minprio && (currconsgroup < nconsgroups || currvarboundgroup < nvargroups));
      count = 0;

      /* adding the groups of constraint with the current priority */
      for( ; currconsgroup < nconsgroups && copiedconsgroupprios[currconsgroup] == currprio; ++currconsgroup )
      {
         SCIP_VAR* binvar;

         SCIPdebugMessage("adding group of constraint <%i> to SCIP\n", consgroupindices[currconsgroup]);

         binvar = SCIPgetBinaryVarSuperindicator(supindconsgroups[consgroupindices[currconsgroup]][0]);
         SCIP_CALL( SCIPgetNegatedVar(scip, binvar, &binvar) );

         /* NOTE: the groups of superindicator constraint are not in the same order as the priority array they correspond to */
         SCIP_CALL( SCIPchgVarObj(scip, binvar, consgroupweights[consgroupindices[currconsgroup]]) );

         stockbinvars[count] = binvar;
         stockweights[count] = consgroupweights[consgroupindices[currconsgroup]];

         /* add each constraint */
         for( c = 0; c < consgroupsizes[consgroupindices[currconsgroup]]; ++c )
         {
            SCIPdebugMessage("constraint <%s> added to SCIP\n",
               SCIPconsGetName(supindconsgroups[consgroupindices[currconsgroup]][c]));
            SCIP_CALL( SCIPaddCons(scip, supindconsgroups[consgroupindices[currconsgroup]][c]) );
         }

         ++count;
      }

      /* adding the groups of variable bound with the current priority */
      for( ; currvarboundgroup < nvargroups && copiedvargroupprios[currvarboundgroup] == currprio; ++currvarboundgroup )
      {
         SCIP_VAR* binvar;

         SCIPdebugMessage("adding group of variable bound <%i> to SCIP\n", vargroupindices[currvarboundgroup]);

         binvar = SCIPgetBinaryVarSuperindicator(supindvarboundgroups[vargroupindices[currvarboundgroup]][0]);
         SCIP_CALL( SCIPgetNegatedVar(scip, binvar, &binvar) );

         /* the groups of superindicator constraints are not in the same order as the priority array they correspond to */
         SCIP_CALL( SCIPchgVarObj(scip, binvar, vargroupweights[vargroupindices[currvarboundgroup]]) );

         stockbinvars[count] = binvar;
         stockweights[count] = vargroupweights[vargroupindices[currvarboundgroup]];

         /* add each constraint */
         for( c = 0; c < vargroupsizes[vargroupindices[currvarboundgroup]]; ++c )
         {
            SCIPdebugMessage("variable bound <%s> added to SCIP\n",
               SCIPconsGetName(supindvarboundgroups[vargroupindices[currvarboundgroup]][c]));
            SCIP_CALL( SCIPaddCons(scip, supindvarboundgroups[vargroupindices[currvarboundgroup]][c]) );
         }

         ++count;
      }

      /* give the modified previous solution to SCIP before solving it */
      if( nsols != 0 )
      {
         SCIP_Bool stored;

         for(c = 0; c < nsols; ++c)
         {
            /* the solution has been messed up.*/
            SCIP_CALL( SCIPaddSolFree(scip, &(sols[c]), &stored) );
         }
         SCIPfreeMemoryArray(scip, &sols);
      }

      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "launching the solving round for sets of priority %i.\n", currprio);
      /* setting the next priority value */
      if( currconsgroup < nconsgroups || currvarboundgroup < nvargroups)
      {
         assert(currprio > minprio);
         if( currconsgroup >= nconsgroups )
         {
            assert(currprio > copiedvargroupprios[currvarboundgroup]);
            currprio = copiedvargroupprios[currvarboundgroup];
         }
         else if( currvarboundgroup >= nvargroups )
         {
            assert( currprio > copiedconsgroupprios[currconsgroup]);
            currprio = copiedconsgroupprios[currconsgroup];
         }
         else  if( copiedconsgroupprios[currconsgroup] >= copiedvargroupprios[currvarboundgroup] )
         {
            assert(currprio > copiedconsgroupprios[currconsgroup]);
            currprio = copiedconsgroupprios[currconsgroup];
         }
         else
         {
            assert(currprio > copiedvargroupprios[currvarboundgroup]);
            currprio = copiedvargroupprios[currvarboundgroup];
         }
      }
      else
      {
         assert(currconsgroup == nconsgroups  && currvarboundgroup == nvargroups );
         assert(currprio == minprio);
         lastround = TRUE;
      }

      if(lastround && !optorig)
      {
         SCIP_CALL( SCIPsetIntParam(scip, "limits/maxorigsol", 10) );
      }
      SCIP_CALL( SCIPsolve(scip) );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,"end of the solving round.\n");

      *status = SCIPgetStatus(scip);
      /* we can have an infeasible problem if some constraint are not relaxable */
      if( *status <= SCIP_STATUS_BESTSOLLIMIT )
      {
         goto TERMINATE1;
      }

      nsols = SCIPgetNSols(scip);
      sols = SCIPgetSols(scip);

      /* can end here if there is no more scip instance to be solved */
      if(lastround && !optorig)
         break;

      /* remembering the optimal objective value of the current stage (the first solution in the solution array is always optimal) */
      objval = SCIPgetSolOrigObj(scip, *sols);
      /* we will only use the optimal solutions */
      for(c = 0; c < nsols; ++c)
      {
         assert( SCIPgetSolOrigObj(scip, sols[c]) >= objval);
         if( objval != SCIPgetSolOrigObj(scip, sols[c]))
         {
            /* change nsols to (c+1)-1 as only the previous one was optimal */
            nsols = c;
            break;
         }
      }

      SCIPdebugMessage("number of optimal solutions is  %i.\n", nsols);

      /* correcting the array of sols to be feasible for the next set of constraints to be added */
      SCIP_CALL( SCIPreallocMemoryArray(scip, &solsvals, nsols) );
      nvars = SCIPgetNOrigVars(scip);
      vars = SCIPgetOrigVars(scip);

      for(c = 0; c < nsols; ++c)
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(solsvals[c]), nvars) );
         SCIP_CALL( SCIPgetSolVals(scip, sols[c], nvars, vars, solsvals[c]) );
      }

      /* setting the value for the next set of constraint to make the current solution feasible in the next round.
       * cannot be done in problem stage (need to check the constraint) */
      if( !lastround )
      {
         int d;
         int e;
         /* update the array of solutions for each constraints of each group to be added in the next loop */
         for( c = currconsgroup; c < nconsgroups && copiedconsgroupprios[c] == currprio; ++c )
         {
            SCIPdebugMessage("updating the solutions for the group of constraints %i of priority %i.\n",
               consgroupindices[c], currprio);

            /* for each solution, we check the binvar for the whole group */
            for(d = 0; d < nsols; ++d)
            {
               SCIP_VAR* binvar;
               SCIP_RESULT result;
               int varidx;

               binvar = SCIPgetBinaryVarSuperindicator(supindconsgroups[consgroupindices[c]][0]);
               SCIP_CALL( SCIPgetNegatedVar(scip, binvar, &binvar) );

               assert( SCIPvarIsActive(binvar) );

               varidx = SCIPvarGetProbindex(binvar);
               assert(findvarpos(binvar, vars, nvars) == varidx);

               solsvals[d][varidx] = 0.0;

               for( e = 0; e < consgroupsizes[consgroupindices[c]]; ++e )
               {
                  /* checking the feasibility of each constraint for the given solution */
                  SCIP_CALL( SCIPcheckCons(scip, SCIPgetSlackConsSuperindicator(supindconsgroups[consgroupindices[c]][e]), sols[d], TRUE, TRUE, FALSE, &result) );
                  SCIPdebugMessage("checking the solution for the constraint <%s>.\n",
                     SCIPconsGetName(SCIPgetSlackConsSuperindicator(supindconsgroups[consgroupindices[c]][e])));

                  if( result == SCIP_INFEASIBLE )
                  {
                     solsvals[d][varidx] = 1.0;
                     break;
                  }
               }
            }
         }

         /* update the array of solution for the next variable bounds to be added to the problem */
         for( c = currvarboundgroup; c < nvargroups &&  copiedvargroupprios[c] == currprio; ++c )
         {
            SCIPdebugMessage("updating the solutions for the group of variable bounds <%i>.\n", vargroupindices[c]);

            for(d = 0; d < nsols; ++d)
            {
               SCIP_VAR* binvar;
               SCIP_RESULT result;
               int varidx;

               binvar =  SCIPgetBinaryVarSuperindicator(supindvarboundgroups[vargroupindices[c]][0]);
               SCIP_CALL( SCIPgetNegatedVar(scip, binvar, &binvar) );

               assert( SCIPvarIsActive(binvar) );

               varidx = SCIPvarGetProbindex(binvar);
               assert(findvarpos(binvar, vars, nvars) == varidx);

               solsvals[d][varidx] = 0.0;

               for( e = 0; e < vargroupsizes[vargroupindices[c]]; ++e )
               {
                  /* checking the feasibility of each variable bound for the given solution */
                  SCIP_CALL( SCIPcheckCons(scip, SCIPgetSlackConsSuperindicator(supindvarboundgroups[vargroupindices[c]][e]),
                        sols[d], TRUE, TRUE, FALSE, &result) );
                  SCIPdebugMessage("checking the solution for the constraint <%s>.\n",
                     SCIPconsGetName(SCIPgetSlackConsSuperindicator(supindvarboundgroups[vargroupindices[c]][e])));
                  if( result == SCIP_INFEASIBLE )
                  {
                     solsvals[d][varidx] = 1.0;
                     break;
                  }
               }
            }
         }
      }

      assert(!lastround || optorig);
      /* must be called before creating any new constraint and after having worked on the solutions.
       * destroy my solution arrays (the vals and the valid flag) */
      SCIP_CALL( SCIPfreeTransform(scip) );

      /* creating the sols and setting their values */
      SCIP_CALL( SCIPallocMemoryArray(scip, &sols, nsols) );
      for( c = 0; c < nsols; ++c)
      {
         SCIP_CALL( SCIPcreateOrigSol(scip, &(sols[c]), NULL) );
         SCIP_CALL( SCIPsetSolVals(scip, sols[c], nvars, vars, solsvals[c]) );
         SCIPfreeMemoryArray(scip, &solsvals[c]);
      }
      /* creating the constraint corresponding to the last priority to ensure the optimality of
       * the relaxation of the current priority set*/
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "priority_set_%i", currprio);

      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, count, stockbinvars, stockweights, objval-previousobjval, objval-previousobjval,
            FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      previousobjval = objval;

      SCIPfreeMemoryArray(scip, &solsvals);
   }
   while (!lastround);

   if( optorig )
   {
      SCIP_Bool stored;

      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "solving the realaxed original problem with original objective function\n");

      /* adding the last set of solutions */
      for(c = 0; c < nsols; ++c)
      {
         SCIP_CALL( SCIPaddSolFree(scip, &(sols[c]), &stored) );
      }
      SCIP_CALL( SCIPreallocMemoryArray(scip, &sols, 1) );

      /* change all the objective value to 0 */
      nvars = SCIPgetNVars(scip);
      vars = SCIPgetVars(scip);

      for( c = 0; c < nvars; ++c)
      {
         SCIP_CALL( SCIPchgVarObj(scip, vars[c], 0) );
      }

      /* change the objective sense and the objective value accordingly to the original problem */
      if(!origprobmin)
      {
         SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );
      }

      for(c = 0; c < norigvars ; ++c)
      {
         SCIP_VAR* var;

         var = origvars[c];
         SCIP_CALL( SCIPchgVarObj(scip, var, origobjvals[c]) );
      }

      SCIP_CALL( SCIPsolve(scip) );

      *status = SCIPgetStatus(scip);

      /* we can have an infeasible problem if some constraint are not relaxable */
      if( *status <= SCIP_STATUS_BESTSOLLIMIT )
      {
         goto TERMINATE2;
      }
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "end of the solving of the original problem\n");
   }

   /* setting the returned value (only work on the best solution for now on)*/
   sols[0] = SCIPgetBestSol(scip);
   if(sol != NULL)
      *sol = sols[0];

   if(consgroupsslacks != NULL)
   {
      for( c = 0; c < nconsgroups + nvargroups; ++c)
      {
         /* NOTE: I cannot check whether the memory was correctly allocated here */
         consgroupsslacks[c] = SCIPgetSolVal(scip, sols[0], allbinvars[c]);
      }
   }

   if(consviols != NULL)
   {
      SCIP_RESULT result;
      int d;
      int countvar;

      count = 0;
      countvar = 0;
      for(c = 0; c < nconsgroups; ++c)
      {
         for(d = 0; d < consgroupsizes[c]; ++d)
         {
            SCIP_CALL( SCIPcheckCons(scip, consgroups[c][d], sols[0], TRUE, TRUE, FALSE, &result) );
            assert(result == SCIP_FEASIBLE || result == SCIP_INFEASIBLE);
            if(result == SCIP_FEASIBLE)
               consviols[count] = 0;
            else
               consviols[count] = 1;
            ++count;
         }
      }

      for(c = 0; c < nvargroups; ++c)
      {
         for(d = 0; d < vargroupsizes[c]; ++d)
         {
            SCIP_CALL( SCIPcheckCons(scip, varboundconss[countvar], sols[0], TRUE, TRUE, FALSE, &result) );
            assert(result == SCIP_FEASIBLE || result == SCIP_INFEASIBLE);
            if(result == SCIP_FEASIBLE)
               consviols[count] = 0;
            else
               consviols[count] = 1;
            ++count;
            ++countvar;
         }
      }
   }

   *status = SCIPgetStatus(scip);

   /* freeing the allocated memory for the pointer */
 TERMINATE2:
   if(optorig)
      SCIPfreeMemoryArray(scip, &sols);

 TERMINATE1:
   if(!vargroupisempty && consviols != NULL)
      SCIPfreeMemoryArray(scip, &varboundconss);

   if( !consgroupisempty )
   {
      for( c = 0; c < nconsgroups; ++c)
      {
         int d;
         for(d = 0; d < consgroupsizes[c]; ++d )
         {
            SCIP_CALL( SCIPreleaseCons(scip, &(supindconsgroups[c][d])) );
         }
         SCIPfreeMemoryArray(scip, &(supindconsgroups[c]));
      }

      SCIPfreeMemoryArray(scip, &supindconsgroups);
      SCIPfreeMemoryArray(scip, &copiedconsgroups);
      SCIPfreeMemoryArray(scip, &copiedconsgroupprios);
      SCIPfreeMemoryArray(scip, &consgroupindices);
   }

   if( !vargroupisempty )
   {
      for( c = 0; c < nvargroups; ++c )
      {
         int d;
         for( d = 0; d < vargroupsizes[c]; ++d )
         {
            SCIP_CALL( SCIPreleaseCons(scip, &(supindvarboundgroups[c][d])) );
         }
         SCIPfreeMemoryArray(scip, &(supindvarboundgroups[c]));
      }

      SCIPfreeMemoryArray(scip, &supindvarboundgroups);
      SCIPfreeMemoryArray(scip, &copiedvargroups);
      SCIPfreeMemoryArray(scip, &copiedvargroupprios);
      SCIPfreeMemoryArray(scip, &vargroupindices);
   }

   SCIPfreeMemoryArray(scip, &origobjvals);
   SCIPfreeMemoryArray(scip, &origvars);
   SCIPfreeMemoryArray(scip, &stockbinvars);
   SCIPfreeMemoryArray(scip, &stockweights);
   SCIPfreeMemoryArray(scip, &allbinvars);

   return SCIP_OKAY;
}

/** runs unit test for SCIPsolveSlack() method */
static
SCIP_RETCODE testslack(
   SCIP*                 scip,               /* SCIP data structure */
   int                   wsizeconsgroups,    /* wanted size for the consgroups */
   int                   wsizevargroups,     /* wanted size for the vargroups */
   SCIP_Bool             returnvalue,        /* do we return any value */
   SCIP_Bool             optorig             /* do we solve the original problem */
   )
{
   SCIP_CONS*** consgroups;
   int nconsgroups;
   int* consgroupsizes;
   int* consgroupprios;
   SCIP_Real* consgroupweights;
   SCIP_VAR*** vargroups;
   int nvargroups;
   int* vargroupsizes;
   SCIP_SLACKTYPE** varslacktypes;
   int* vargroupprios;
   SCIP_Real* vargroupweights;
   /* return arguments */
   SCIP_SOL** sols;
   SCIP_STATUS status;
   /* working argument */
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_VAR** usedvars;
   int nusedvars;
   int nvars;
   int nconss;
   int c;
   int d;
   int count;
   SCIP_Real* consgroupsslacks;
   SCIP_SOL* sol;
   SCIP_Real* consviols;

   SCIP_CALL( SCIPallocMemory(scip, &sols) );

   /* initialize the memory needed */

   nconss = SCIPgetNConss(scip);
   nvars = SCIPgetNVars(scip);

   conss = SCIPgetConss(scip);
   vars = SCIPgetVars(scip);

   SCIP_CALL( SCIPallocMemoryArray(scip, &usedvars, nvars) );

   /* sort the array of variable and only keep the non linear one */
   nusedvars = 0;
   for (c = 0; c < nvars; ++c)
   {
      if ( SCIP_VARTYPE_BINARY != SCIPvarGetType(vars[c]) )
      {
         usedvars[nusedvars] = vars[c];
         nusedvars++;
      }
   }

   /* transform the value to create some consgroups and vargroups. */
   /* want ceil(nvars/wsizevargroups) and for cons as well */
   nvargroups = (nusedvars + wsizevargroups -1) / wsizevargroups;
   nconsgroups = (nconss + wsizeconsgroups -1) / wsizeconsgroups;

   SCIP_CALL( SCIPallocMemoryArray(scip, &consgroups, nconsgroups) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consgroupsizes, nconsgroups) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consgroupprios, nconsgroups) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consgroupweights, nconsgroups) );

   if ( nvargroups > 0)
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &vargroups, nvargroups) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &vargroupsizes, nvargroups) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &vargroupprios, nvargroups) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &vargroupweights, nvargroups) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &varslacktypes, nvargroups) );
   }
   else
   {
      vargroups = NULL;
      vargroupsizes = NULL;
      vargroupprios = NULL;
      vargroupweights = NULL;
      varslacktypes = NULL;
   }
   SCIP_CALL( SCIPallocMemoryArray(scip, &consgroupsslacks, nconsgroups + nvargroups) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consviols, nconss + nvars) );

   count = 0;
   for (c = 0; c < nconsgroups ; ++c)
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &consgroups[c], wsizeconsgroups) );

      for (d = 0; d < wsizeconsgroups && count < nconss; ++d)
      {
         consgroups[c][d] = conss[count];
         ++count;
      }
      consgroupsizes[c] = d;
      consgroupprios[c] = c;
      consgroupweights[c] = 1;
   }

   count = 0;
   for (c = 0; c < nvargroups ; ++c)
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &vargroups[c], wsizevargroups) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &varslacktypes[c], wsizevargroups) );

      for (d = 0; d < wsizevargroups && count < nusedvars; ++d)
      {
         vargroups[c][d] = usedvars[count];
         varslacktypes[c][d] = SCIP_SLACKTYPE_BOTH;
         ++count;
      }
      vargroupsizes[c] = d;
      vargroupprios[c] = c;
      vargroupweights[c] = 3;
   }

   if (returnvalue)
   {
      SCIP_CALL( SCIPsolveSlack(scip, consgroups, nconsgroups, consgroupsizes, consgroupprios, consgroupweights, vargroups, nvargroups,
            vargroupsizes, varslacktypes, vargroupprios, vargroupweights, consgroupsslacks, &sol, consviols, optorig, &status) );
   }
   else
   {
      SCIP_CALL( SCIPsolveSlack(scip, consgroups, nconsgroups, consgroupsizes, consgroupprios, consgroupweights, vargroups, nvargroups,
            vargroupsizes, varslacktypes, vargroupprios, vargroupweights, NULL, NULL, NULL, optorig, &status) );
   }

   /* free the memory */
   for (c = 0; c < nconsgroups ; ++c)
   {
      SCIPfreeMemoryArray(scip, &consgroups[c]);
   }

   for (c = 0; c < nvargroups ; ++c)
   {
      SCIPfreeMemoryArray(scip, &vargroups[c]);
      SCIPfreeMemoryArray(scip, &varslacktypes[c]);
   }

   SCIPfreeMemory(scip, &sols);

   SCIPfreeMemoryArray(scip, &consgroups);
   SCIPfreeMemoryArray(scip, &consgroupsizes);
   SCIPfreeMemoryArray(scip, &consgroupprios);
   SCIPfreeMemoryArray(scip, &consgroupweights);

   if (nvargroups > 0)
   {
      SCIPfreeMemoryArray(scip, &vargroups);
      SCIPfreeMemoryArray(scip, &vargroupsizes);
      SCIPfreeMemoryArray(scip, &vargroupprios);
      SCIPfreeMemoryArray(scip, &vargroupweights);
      SCIPfreeMemoryArray(scip, &varslacktypes);
   }

   SCIPfreeMemoryArray(scip, &usedvars);
   SCIPfreeMemoryArray(scip, &consgroupsslacks);
   SCIPfreeMemoryArray(scip, &consviols);

   return SCIP_OKAY;
}

#include "include/scip_test.h"

static SCIP* scip;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!");
}

TestSuite(slack, .init = setup, .fini = teardown);

Test(slack, feasible)
{
   // read a feasible model
   SCIP_CALL( SCIPreadProb(scip, "src/cons/superindicator/misc07.mps.gz", NULL) );

   // call testslack
   //SCIP_CALL( testslack(scip, 13, 17, TRUE, TRUE) );
   //SCIP_CALL( testslack(scip, 1, 1, TRUE, TRUE) );
}

Test(slack, infeasible)
{
   // read a infeasible model
   SCIP_CALL( SCIPreadProb(scip, "src/cons/superindicator/gen_inf.mps.gz", NULL) );
   // call testslack
   SCIP_CALL( testslack(scip, 13, 17, TRUE, TRUE) );
}

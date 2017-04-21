/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   benders.c
 * @brief  methods for Benders' decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/clock.h"
#include "scip/paramset.h"
#include "scip/lp.h"
#include "scip/prob.h"
#include "scip/pricestore.h"
#include "scip/scip.h"
#include "scip/benders.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"

#include "scip/struct_benders.h"
#include "scip/struct_benderscut.h"

#include "scip/benderscut.h"

/* Local methods */

/** adds the auxiliary variables to the Benders' decomposition master problem */
static
SCIP_RETCODE addAuxiliaryVariablesToMaster(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< benders */
   )
{
   SCIP_VAR* auxiliaryvar;
   char varname[SCIP_MAXSTRLEN];    /* the name of the auxiliary variable */
   int i;

   for( i = 0; i < SCIPbendersGetNSubproblems(benders); i++ )
   {
      /* if no optimality cuts have been added for this subproblem, then the auxiliary variable will be created and
       * added */
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "auxiliaryvar_%d", i );
      SCIP_CALL( SCIPcreateVarBasic(scip, &auxiliaryvar, varname, -SCIPinfinity(scip), SCIPinfinity(scip), 1.0,
            SCIP_VARTYPE_CONTINUOUS) );

      SCIP_CALL( SCIPaddVar(scip, auxiliaryvar) );

      benders->auxiliaryvars[i] = auxiliaryvar;

      SCIP_CALL( SCIPreleaseVar(scip, &auxiliaryvar) );
   }

   return SCIP_OKAY;
}



/** compares two benders w. r. to their priority */
SCIP_DECL_SORTPTRCOMP(SCIPbendersComp)
{  /*lint --e{715}*/
   return ((SCIP_BENDERS*)elem2)->priority - ((SCIP_BENDERS*)elem1)->priority;
}

/** comparison method for sorting benders w.r.t. to their name */
SCIP_DECL_SORTPTRCOMP(SCIPbendersCompName)
{
   return strcmp(SCIPbendersGetName((SCIP_BENDERS*)elem1), SCIPbendersGetName((SCIP_BENDERS*)elem2));
}

/** method to call, when the priority of a benders was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdBendersPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetBendersPriority() to mark the benderss unsorted */
   SCIP_CALL( SCIPsetBendersPriority(scip, (SCIP_BENDERS*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** copies the given benders to a new scip */
SCIP_RETCODE SCIPbendersCopyInclude(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_SET*             set,                /**< SCIP_SET of SCIP to copy to */
   SCIP_Bool*            valid               /**< was the copying process valid? */
   )
{
   int i;

   assert(benders != NULL);
   assert(set != NULL);
   assert(valid != NULL);
   assert(set->scip != NULL);

   if( benders->benderscopy != NULL )
   {
      SCIPsetDebugMsg(set, "including benders %s in subscip %p\n", SCIPbendersGetName(benders), (void*)set->scip);
      SCIP_CALL( benders->benderscopy(set->scip, benders, valid) );
   }

   /* calling the copy method for the Benders' cuts */
   SCIPbendersSortBenderscuts(benders);
   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutCopyInclude(benders->benderscuts[i], set) );
   }

   return SCIP_OKAY;
}

/** creates a Benders' decomposition structure
 *  To use the Benders' decomposition for solving a problem, it first has to be activated with a call to SCIPactivateBenders().
 */
SCIP_RETCODE SCIPbendersCreate(
   SCIP_BENDERS**        benders,            /**< pointer to Benders' decomposition data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of Benders' decomposition */
   const char*           desc,               /**< description of Benders' decomposition */
   int                   priority,           /**< priority of the Benders' decomposition */
   int                   nsubproblems,       /**< the number subproblems used in this decomposition */
   SCIP_Bool             cutlp,              /**< should Benders' cuts be generated for LP solutions */
   SCIP_Bool             cutpseudo,          /**< should Benders' cuts be generated for pseudo solutions */
   SCIP_Bool             cutrelax,           /**< should Benders' cuts be generated for relaxation solutions */
   SCIP_DECL_BENDERSCOPY ((*benderscopy)),   /**< copy method of benders or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BENDERSFREE ((*bendersfree)),   /**< destructor of Benders' decomposition */
   SCIP_DECL_BENDERSINIT ((*bendersinit)),   /**< initialize Benders' decomposition */
   SCIP_DECL_BENDERSEXIT ((*bendersexit)),   /**< deinitialize Benders' decomposition */
   SCIP_DECL_BENDERSINITPRE((*bendersinitpre)),/**< presolving initialization method for Benders' decomposition */
   SCIP_DECL_BENDERSEXITPRE((*bendersexitpre)),/**< presolving deinitialization method for Benders' decomposition */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol)),/**< solving process initialization method of Benders' decomposition */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol)),/**< solving process deinitialization method of Benders' decomposition */
   SCIP_DECL_BENDERSGETMASTERVAR((*bendersgetmastervar)),/**< returns the master variable for a given subproblem variable */
   SCIP_DECL_BENDERSEXEC ((*bendersexec)),   /**< the execution method of the Benders' decomposition algorithm */
   SCIP_DECL_BENDERSSOLVESUB((*benderssolvesub)),/**< the solving method for the Benders' decomposition subproblems */
   SCIP_DECL_BENDERSPOSTSOLVE((*benderspostsolve)),/**< called after the subproblems are solved. */
   SCIP_DECL_BENDERSFREESUB((*bendersfreesub)),/**< the freeing method for the Benders' decomposition subproblems */
   SCIP_BENDERSDATA*     bendersdata         /**< Benders' decomposition data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(benders != NULL);
   assert(name != NULL);
   assert(desc != NULL);

   SCIP_ALLOC( BMSallocMemory(benders) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*benders)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*benders)->desc, desc, strlen(desc)+1) );
   (*benders)->priority = priority;
   (*benders)->nsubproblems = nsubproblems;
   (*benders)->cutlp = cutlp;
   (*benders)->cutpseudo = cutpseudo;
   (*benders)->cutrelax = cutrelax;
   (*benders)->benderscopy = benderscopy;
   (*benders)->bendersfree = bendersfree;
   (*benders)->bendersinit = bendersinit;
   (*benders)->bendersexit = bendersexit;
   (*benders)->bendersinitpre = bendersinitpre;
   (*benders)->bendersexitpre = bendersexitpre;
   (*benders)->bendersinitsol = bendersinitsol;
   (*benders)->bendersexitsol = bendersexitsol;
   (*benders)->bendersgetmastervar = bendersgetmastervar;
   (*benders)->bendersexec = bendersexec;
   (*benders)->benderssolvesub = benderssolvesub;
   (*benders)->benderspostsolve = benderspostsolve;
   (*benders)->bendersfreesub = bendersfreesub;
   (*benders)->bendersdata = bendersdata;
   SCIP_CALL( SCIPclockCreate(&(*benders)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*benders)->bendersclock, SCIP_CLOCKTYPE_DEFAULT) );
   (*benders)->ncalls = 0;
   (*benders)->noptcutsfound = 0;
   (*benders)->nfeascutsfound = 0;
   (*benders)->initialized = FALSE;
   (*benders)->addedsubprobs = 0;

   (*benders)->benderscuts = NULL;
   (*benders)->nbenderscuts = 0;
   (*benders)->benderscutssize = 0;
   (*benders)->benderscutssorted = FALSE;
   (*benders)->benderscutsnamessorted = FALSE;

   /* allocating memory for the subproblems array */
   SCIP_ALLOC( BMSallocMemoryArray(&(*benders)->subproblems, (*benders)->nsubproblems) );

   /* allocating memory for the auxiliary variables array */
   SCIP_ALLOC( BMSallocMemoryArray(&(*benders)->auxiliaryvars, (*benders)->nsubproblems) );

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of benders <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*benders)->priority, FALSE, priority, INT_MIN/4, INT_MAX/4,
                  paramChgdBendersPriority, (SCIP_PARAMDATA*)(*benders)) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutlp", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "should Benders' cuts be generated for LP solutions?");
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*benders)->cutlp, FALSE, cutlp, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutpseudo", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "should Benders' cuts be generated for pseudo solutions?");
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*benders)->cutpseudo, FALSE, cutpseudo, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutrelax", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "should Benders' cuts be generated for relaxation solutions?");
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*benders)->cutrelax, FALSE, cutrelax, NULL, NULL) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** calls destructor and frees memory of Benders' decomposition */
SCIP_RETCODE SCIPbendersFree(
   SCIP_BENDERS**        benders,            /**< pointer to Benders' decomposition data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(benders != NULL);
   assert(*benders != NULL);
   assert(!(*benders)->initialized);
   assert(set != NULL);

   /* call destructor of Benders' decomposition */
   if( (*benders)->bendersfree != NULL )
   {
      SCIP_CALL( (*benders)->bendersfree(set->scip, *benders) );
   }

   /* freeing the Benders' cuts */
   for( i = 0; i < (*benders)->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutFree(&((*benders)->benderscuts[i]), set) );
   }
   BMSfreeMemoryArrayNull(&(*benders)->benderscuts);


   BMSfreeMemoryArray(&(*benders)->auxiliaryvars);
   BMSfreeMemoryArray(&(*benders)->subproblems);

   SCIPclockFree(&(*benders)->bendersclock);
   SCIPclockFree(&(*benders)->setuptime);
   BMSfreeMemoryArray(&(*benders)->name);
   BMSfreeMemoryArray(&(*benders)->desc);
   BMSfreeMemory(benders);

   return SCIP_OKAY;
}

/** initializes Benders' decomposition */
SCIP_RETCODE SCIPbendersInit(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   if( benders->initialized )
   {
      SCIPerrorMessage("Benders' decomposition <%s> already initialized\n", benders->name);
      return SCIP_INVALIDCALL;
   }

   if( set->misc_resetstat )
   {
      SCIPclockReset(benders->setuptime);
      SCIPclockReset(benders->bendersclock);

      benders->ncalls = 0;
      benders->noptcutsfound = 0;
      benders->nfeascutsfound = 0;
   }

   if( benders->bendersinit != NULL )
   {
      /* start timing */
      SCIPclockStart(benders->setuptime, set);

      SCIP_CALL( benders->bendersinit(set->scip, benders) );

      /* stop timing */
      SCIPclockStop(benders->setuptime, set);
   }

   /* initialising the Benders' cuts */
   SCIPbendersSortBenderscuts(benders);
   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutInit(benders->benderscuts[i], set) );
   }

   benders->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of Benders' decomposition */
SCIP_RETCODE SCIPbendersExit(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   if( !benders->initialized )
   {
      SCIPerrorMessage("Benders' decomposition <%s> not initialized\n", benders->name);
      return SCIP_INVALIDCALL;
   }

   if( benders->bendersexit != NULL )
   {
      /* start timing */
      SCIPclockStart(benders->setuptime, set);

      SCIP_CALL( benders->bendersexit(set->scip, benders) );

      /* stop timing */
      SCIPclockStop(benders->setuptime, set);
   }

   /* calling the exit method for the Benders' cuts */
   SCIPbendersSortBenderscuts(benders);
   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutExit(benders->benderscuts[i], set) );
   }

   benders->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs the Benders' decomposition that the presolving process is being started */
SCIP_RETCODE SCIPbendersInitpre(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(benders != NULL);
   assert(set != NULL);
   assert(stat != NULL);

   /* adding the auxiliary variables to the master problem */
   SCIP_CALL( addAuxiliaryVariablesToMaster(set->scip, benders) );

   /* call presolving initialization method of Benders' decomposition */
   if( benders->bendersinitpre != NULL )
   {
      /* start timing */
      SCIPclockStart(benders->setuptime, set);

      SCIP_CALL( benders->bendersinitpre(set->scip, benders) );

      /* stop timing */
      SCIPclockStop(benders->setuptime, set);
   }

   return SCIP_OKAY;
}


/** informs the Benders' decomposition that the presolving process has completed */
SCIP_RETCODE SCIPbendersExitpre(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(benders != NULL);
   assert(set != NULL);
   assert(stat != NULL);

   /* call presolving  deinitialization method of Benders' decomposition */
   if( benders->bendersexitpre != NULL )
   {
      /* start timing */
      SCIPclockStart(benders->setuptime, set);

      SCIP_CALL( benders->bendersexitpre(set->scip, benders) );

      /* stop timing */
      SCIPclockStop(benders->setuptime, set);
   }



   return SCIP_OKAY;
}

/** informs Benders' decomposition that the branch and bound process is being started */
SCIP_RETCODE SCIPbendersInitsol(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   /* call solving process initialization method of Benders' decomposition */
   if( benders->bendersinitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(benders->setuptime, set);

      SCIP_CALL( benders->bendersinitsol(set->scip, benders) );

      /* stop timing */
      SCIPclockStop(benders->setuptime, set);
   }

   /* calling the initsol method for the Benders' cuts */
   SCIPbendersSortBenderscuts(benders);
   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutInitsol(benders->benderscuts[i], set) );
   }

   return SCIP_OKAY;
}

/** informs Benders' decomposition that the branch and bound process data is being freed */
SCIP_RETCODE SCIPbendersExitsol(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of Benders' decomposition */
   if( benders->bendersexitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(benders->setuptime, set);

      SCIP_CALL( benders->bendersexitsol(set->scip, benders) );

      /* stop timing */
      SCIPclockStop(benders->setuptime, set);
   }

   /* calling the exitsol method for the Benders' cuts */
   SCIPbendersSortBenderscuts(benders);
   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutExitsol(benders->benderscuts[i], set) );
   }

   return SCIP_OKAY;
}

/** solves the subproblem using the current master problem solution. */
/*  TODO: consider allowing the possibility to pass solution information back from the subproblems instead of the scip
 *  instance. This would allow the use of different solvers for the subproblems, more importantly allowing the use of an
 *  LP solver for LP subproblems. */
SCIP_RETCODE SCIPbendersExec(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_RESULT*          result,             /**< result of the pricing process */
   SCIP_Bool             check               /**< is the execution method called as a check. i.e. no cuts are required */
   )
{
   SCIP_BENDERSCUT** benderscuts;
   int nsubproblems;
   int nbenderscuts;
   int i;
   int j;
   SCIP_Bool infeasible = FALSE;

   assert(benders != NULL);
   assert(result != NULL);
   assert(benders->bendersexec != NULL);

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   SCIP_CALL( benders->bendersexec(set->scip, benders) );

   *result = SCIP_DIDNOTRUN;
   /* solving each of the subproblems for Benders decomposition */
   /* TODO: ensure that the each of the subproblems solve and update the parameters with the correct return values */
   for( i = 0; i < nsubproblems; i++ )
   {
      SCIP_Bool subinfeas = FALSE;

      SCIP_CALL( SCIPbendersSolveSubproblem(benders, set, sol, i, &subinfeas) );

      infeasible = infeasible || subinfeas;

      /* if the subproblems are being solved as part of the conscheck, then we break once an infeasibility is found. The
       * result pointer is set to infeasible and the execution is halted. */
      if( check && infeasible )
      {
         (*result) = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   }

   if( check )
   {
      /* if the subproblems are being solved as part of conscheck, then the results flag must be returned after the solving
       * has completed. No cut is generated during conscheck. */
      if( infeasible )
         (*result) = SCIP_INFEASIBLE;
      else
         (*result) = SCIP_FEASIBLE;
   }
   else
   {
      /* the result flag is set to FEASIBLE as default. If a cut is added, then this is changed to CONS_ADDED */
      (*result) = SCIP_FEASIBLE;

      /* Generating cuts for the subproblems. */
      benderscuts = SCIPbendersGetBenderscuts(benders);
      nbenderscuts = SCIPbendersGetNBenderscuts(benders);

      /* This is done in two loops. The first is by subproblem and the second is by cut type. */
      for( i = 0; i < nsubproblems; i++ )
      {
         for( j = 0; j < nbenderscuts; j++ )
         {
            assert(benderscuts[j] != NULL);

            SCIP_CALL( SCIPbenderscutExec(benderscuts[j], set, benders, i, result) );
         }
      }
   }

   /* calling the post-solve call back for the Benders' decomposition algorithm. This allows the user to work directly
    * with the solved subproblems and the master problem */
   if( benders->benderspostsolve != NULL )
      SCIP_CALL( benders->benderspostsolve(set->scip, benders, infeasible) );

   /* freeing the subproblems after the cuts are generated */
   for( i = 0; i < benders->nsubproblems; i++ )
      SCIP_CALL( SCIPbendersFreeSubproblem(benders, set, i) );

   /* increment the number of calls to the Benders' decomposition subproblem solve */
   benders->ncalls++;

   return SCIP_OKAY;
}

/** solves the subproblems. */
extern
SCIP_RETCODE SCIPbendersSolveSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnum,            /**< the subproblem number */
   SCIP_Bool*            infeasible          /**< returns whether the current subproblem is infeasible */
   )
{
   assert(benders != NULL);
   assert(benders->benderssolvesub != NULL);
   assert(probnum >= 0 && probnum < benders->nsubproblems);

   SCIP_CALL( benders->benderssolvesub(set->scip, benders, sol, probnum, infeasible) );

   return SCIP_OKAY;
}

/** frees the subproblems. */
SCIP_RETCODE SCIPbendersFreeSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   probnum             /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(benders->bendersfreesub != NULL);
   assert(probnum >= 0 && probnum < benders->nsubproblems);

   SCIP_CALL( benders->bendersfreesub(set->scip, benders, probnum) );

   return SCIP_OKAY;
}

/** returns the variable for the given variable of the subproblem. This provides a call back for the mapping between the
 * master and subproblems */
SCIP_VAR* SCIPbendersGetMasterVar(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< the copy of the master variable in the subproblem */
   )
{
   assert(benders != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(benders->bendersgetmastervar != NULL);

   return benders->bendersgetmastervar(set->scip, benders, var);
}

/** gets user data of Benders' decomposition */
SCIP_BENDERSDATA* SCIPbendersGetData(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->bendersdata;
}

/** sets user data of Benders' decomposition; user has to free old data in advance! */
void SCIPbendersSetData(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_BENDERSDATA*     bendersdata         /**< new Benders' decomposition user data */
   )
{
   assert(benders != NULL);

   benders->bendersdata = bendersdata;
}

/** sets copy callback of benders */
void SCIPbendersSetCopy(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSCOPY ((*benderscopy))    /**< copy callback of benders */
   )
{
   assert(benders != NULL);

   benders->benderscopy = benderscopy;
}

/** sets destructor callback of benders */
void SCIPbendersSetFree(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_DECL_BENDERSFREE ((*bendersfree))    /**< destructor of benders */
   )
{
   assert(benders != NULL);

   benders->bendersfree = bendersfree;
}

/** sets initialization callback of benders */
void SCIPbendersSetInit(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_DECL_BENDERSINIT((*bendersinit))     /**< initialize benders */
   )
{
   assert(benders != NULL);

   benders->bendersinit = bendersinit;
}

/** sets deinitialization callback of benders */
void SCIPbendersSetExit(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_DECL_BENDERSEXIT((*bendersexit))     /**< deinitialize benders */
   )
{
   assert(benders != NULL);

   benders->bendersexit = bendersexit;
}

/** sets presolving initialization callback of benders */
void SCIPbendersSetInitpre(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_DECL_BENDERSINITPRE((*bendersinitpre))     /**< initialize benders */
   )
{
   assert(benders != NULL);

   benders->bendersinitpre = bendersinitpre;
}

/** sets presolving deinitialization callback of benders */
void SCIPbendersSetExitpre(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_DECL_BENDERSEXITPRE((*bendersexitpre))     /**< deinitialize benders */
   )
{
   assert(benders != NULL);

   benders->bendersexitpre = bendersexitpre;
}

/** sets solving process initialization callback of benders */
void SCIPbendersSetInitsol(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol))/**< solving process initialization callback of benders */
   )
{
   assert(benders != NULL);

   benders->bendersinitsol = bendersinitsol;
}

/** sets solving process deinitialization callback of benders */
void SCIPbendersSetExitsol(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol))/**< solving process deinitialization callback of benders */
   )
{
   assert(benders != NULL);

   benders->bendersexitsol = bendersexitsol;
}

/** sets post-solve callback of benders */
void SCIPbendersSetPostsolve(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_DECL_BENDERSPOSTSOLVE((*benderspostsolve))/**< solving process deinitialization callback of benders */
   )
{
   assert(benders != NULL);

   benders->benderspostsolve = benderspostsolve;
}

/** gets name of Benders' decomposition */
const char* SCIPbendersGetName(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->name;
}

/** gets description of Benders' decomposition */
const char* SCIPbendersGetDesc(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->desc;
}

/** gets priority of Benders' decomposition */
int SCIPbendersGetPriority(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->priority;
}

/** sets priority of Benders' decomposition */
void SCIPbendersSetPriority(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the Benders' decomposition */
   )
{
   assert(benders != NULL);
   assert(set != NULL);

   benders->priority = priority;
   set->benderssorted = FALSE;
}

/** gets the number of subproblems for the Benders' decomposition */
int SCIPbendersGetNSubproblems(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition data structure */
   )
{
   assert(benders != NULL);

   return benders->nsubproblems;
}

/** returns the SCIP instance for a given subproblem */
SCIP* SCIPbendersSubproblem(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < benders->nsubproblems);

   return benders->subproblems[probnumber];
}

/** gets the number of times, the benders was called and tried to find a variable with negative reduced costs */
int SCIPbendersGetNCalls(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->ncalls;
}

/** gets the number of optimality cuts found by the collection of Benders' decomposition subproblems */
int SCIPbendersGetNOptCutsFound(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->noptcutsfound;
}

/** gets the number of feasibility cuts found by the collection of Benders' decomposition subproblems */
int SCIPbendersGetNFeasCutsFound(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->nfeascutsfound;
}

/** gets time in seconds used in this benders for setting up for next stages */
SCIP_Real SCIPbendersGetSetupTime(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return SCIPclockGetTime(benders->setuptime);
}

/** gets time in seconds used in this benders */
SCIP_Real SCIPbendersGetTime(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return SCIPclockGetTime(benders->bendersclock);
}

/** enables or disables all clocks of \p benders, depending on the value of the flag */
void SCIPbendersEnableOrDisableClocks(
   SCIP_BENDERS*         benders,            /**< the benders for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the benders be enabled? */
   )
{
   assert(benders != NULL);

   SCIPclockEnableOrDisable(benders->setuptime, enable);
   SCIPclockEnableOrDisable(benders->bendersclock, enable);
}

/** is Benders' decomposition initialized? */
SCIP_Bool SCIPbendersIsInitialized(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->initialized;
}

/** are Benders' cuts generated from the LP solutions? */
SCIP_Bool SCIPbendersCutLP(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->cutlp;
}

/** are Benders' cuts generated from the pseudo solutions? */
SCIP_Bool SCIPbendersCutPseudo(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->cutpseudo;
}

/** are Benders' cuts generated from the relaxation solutions? */
SCIP_Bool SCIPbendersCutRelaxation(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->cutrelax;
}

/** Adds a subproblem to the Benders' decomposition data */
void SCIPbendersAddSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP*                 subproblem          /**< subproblem to be added to the data storage */
   )
{
   assert(benders != NULL);
   assert(subproblem != NULL);
   assert(benders->subproblems != NULL);
   assert(benders->addedsubprobs + 1 <= benders->nsubproblems);

   benders->subproblems[benders->addedsubprobs] = subproblem;

   benders->addedsubprobs++;
}

/** returns the auxiliary variable for the given subproblem */
SCIP_VAR* SCIPbendersGetAuxiliaryVar(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->auxiliaryvars[probnumber];
}

/** sets the sorted flags in the Benders' decomposition */
void SCIPbendersSetBenderscutsSorted(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition structure */
   SCIP_Bool             sorted              /**< the value to set the sorted flag to */
   )
{
   assert(benders != NULL);

   benders->benderscutssorted = sorted;
   benders->benderscutsnamessorted = sorted;
}

/** inserts a Benders' cut into the Benders' cuts list */
SCIP_RETCODE SCIPbendersIncludeBenderscut(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BENDERSCUT*      benderscut          /**< Benders' cut */
   )
{
   assert(benders != NULL);
   assert(benderscut != NULL);

   if( benders->nbenderscuts >= benders->benderscutssize )
   {
      benders->benderscutssize = SCIPsetCalcMemGrowSize(set, benders->nbenderscuts+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&benders->benderscuts, benders->benderscutssize) );
   }
   assert(benders->nbenderscuts < benders->benderscutssize);

   benders->benderscuts[benders->nbenderscuts] = benderscut;
   benders->nbenderscuts++;
   benders->benderscutssorted = FALSE;

   return SCIP_OKAY;
}

/** returns the Benders' cut of the given name, or NULL if not existing */
SCIP_BENDERSCUT* SCIPfindBenderscut(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   const char*           name                /**< name of Benderscut' decomposition */
   )
{
   int i;

   assert(benders != NULL);
   assert(name != NULL);

   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      if( strcmp(SCIPbenderscutGetName(benders->benderscuts[i]), name) == 0 )
         return benders->benderscuts[i];
   }

   return NULL;
}

/** returns the array of currently available Benders' cuts; active benders are in the first slots of the array */
SCIP_BENDERSCUT** SCIPbendersGetBenderscuts(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   if( !benders->benderscutssorted )
   {
      SCIPsortPtr((void**)benders->benderscuts, SCIPbenderscutComp, benders->nbenderscuts);
      benders->benderscutssorted = TRUE;
      benders->benderscutsnamessorted = FALSE;
   }

   return benders->benderscuts;
}

/** returns the number of currently available Benders' cuts */
int SCIPbendersGetNBenderscuts(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->nbenderscuts;
}

/** sets the priority of a Benders' decomposition */
SCIP_RETCODE SCIPbendersSetBenderscutPriority(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' cut */
   int                   priority            /**< new priority of the Benders' decomposition */
   )
{
   assert(benders != NULL);
   assert(benderscut != NULL);

   benderscut->priority = priority;
   benders->benderscutssorted = FALSE;

   return SCIP_OKAY;
}

/** sorts benders cuts by priorities */
void SCIPbendersSortBenderscuts(
   SCIP_BENDERS*         benders             /**< benders */
   )
{
   assert(benders != NULL);

   if( !benders->benderscutssorted )
   {
      SCIPsortPtr((void**)benders->benderscuts, SCIPbenderscutComp, benders->nbenderscuts);
      benders->benderscutssorted = TRUE;
      benders->benderscutsnamessorted = FALSE;
   }
}

/** sorts benders cuts by name */
void SCIPbendersSortBenderscutsName(
   SCIP_BENDERS*         benders             /**< benders */
   )
{
   assert(benders != NULL);

   if( !benders->benderscutsnamessorted )
   {
      SCIPsortPtr((void**)benders->benderscuts, SCIPbenderscutCompName, benders->nbenderscuts);
      benders->benderscutssorted = FALSE;
      benders->benderscutsnamessorted = TRUE;
   }
}

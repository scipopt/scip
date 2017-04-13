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
   assert(benders != NULL);
   assert(set != NULL);
   assert(valid != NULL);
   assert(set->scip != NULL);

   if( benders->benderscopy != NULL )
   {
      SCIPsetDebugMsg(set, "including benders %s in subscip %p\n", SCIPbendersGetName(benders), (void*)set->scip);
      SCIP_CALL( benders->benderscopy(set->scip, benders, valid) );
   }
   return SCIP_OKAY;
}

/** creates a variable benders
 *  To use the variable benders for solving a problem, it first has to be activated with a call to SCIPactivateBenders().
 */
SCIP_RETCODE SCIPbendersCreate(
   SCIP_BENDERS**        benders,            /**< pointer to variable benders data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of variable benders */
   const char*           desc,               /**< description of variable benders */
   int                   priority,           /**< priority of the variable benders */
   int                   nsubproblems,       /**< the number subproblems used in this decomposition */
   SCIP_DECL_BENDERSCOPY ((*benderscopy)),   /**< copy method of benders or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BENDERSFREE ((*bendersfree)),   /**< destructor of variable benders */
   SCIP_DECL_BENDERSINIT ((*bendersinit)),   /**< initialize variable benders */
   SCIP_DECL_BENDERSEXIT ((*bendersexit)),   /**< deinitialize variable benders */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol)),/**< solving process initialization method of variable benders */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol)),/**< solving process deinitialization method of variable benders */
   SCIP_DECL_BENDERSGETMASTERVAR((*bendersgetmastervar)),/**< returns the master variable for a given subproblem variable */
   SCIP_DECL_BENDERSEXEC ((*bendersexec)),   /**< the execution method of the Benders' decomposition algorithm */
   SCIP_BENDERSDATA*     bendersdata         /**< variable benders data */
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
   (*benders)->benderscopy = benderscopy;
   (*benders)->bendersfree = bendersfree;
   (*benders)->bendersinit = bendersinit;
   (*benders)->bendersexit = bendersexit;
   (*benders)->bendersinitsol = bendersinitsol;
   (*benders)->bendersexitsol = bendersexitsol;
   (*benders)->bendersgetmastervar = bendersgetmastervar;
   (*benders)->bendersexec = bendersexec;
   (*benders)->bendersdata = bendersdata;
   SCIP_CALL( SCIPclockCreate(&(*benders)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*benders)->bendersclock, SCIP_CLOCKTYPE_DEFAULT) );
   (*benders)->ncalls = 0;
   (*benders)->noptcutsfound = 0;
   (*benders)->nfeascutsfound = 0;
   (*benders)->initialized = FALSE;

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of benders <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*benders)->priority, FALSE, priority, INT_MIN/4, INT_MAX/4,
                  paramChgdBendersPriority, (SCIP_PARAMDATA*)(*benders)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** calls destructor and frees memory of variable benders */
SCIP_RETCODE SCIPbendersFree(
   SCIP_BENDERS**        benders,            /**< pointer to variable benders data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(benders != NULL);
   assert(*benders != NULL);
   assert(!(*benders)->initialized);
   assert(set != NULL);

   /* call destructor of variable benders */
   if( (*benders)->bendersfree != NULL )
   {
      SCIP_CALL( (*benders)->bendersfree(set->scip, *benders) );
   }

   SCIPclockFree(&(*benders)->bendersclock);
   SCIPclockFree(&(*benders)->setuptime);
   BMSfreeMemoryArray(&(*benders)->name);
   BMSfreeMemoryArray(&(*benders)->desc);
   BMSfreeMemory(benders);

   return SCIP_OKAY;
}

/** initializes variable benders */
SCIP_RETCODE SCIPbendersInit(
   SCIP_BENDERS*         benders,            /**< variable benders */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(benders != NULL);
   assert(set != NULL);

   if( benders->initialized )
   {
      SCIPerrorMessage("variable benders <%s> already initialized\n", benders->name);
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
   benders->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of variable benders */
SCIP_RETCODE SCIPbendersExit(
   SCIP_BENDERS*         benders,            /**< variable benders */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(benders != NULL);
   assert(set != NULL);

   if( !benders->initialized )
   {
      SCIPerrorMessage("variable benders <%s> not initialized\n", benders->name);
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
   benders->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs variable benders that the branch and bound process is being started */
SCIP_RETCODE SCIPbendersInitsol(
   SCIP_BENDERS*         benders,            /**< variable benders */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(benders != NULL);
   assert(set != NULL);

   /* call solving process initialization method of variable benders */
   if( benders->bendersinitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(benders->setuptime, set);

      SCIP_CALL( benders->bendersinitsol(set->scip, benders) );

      /* stop timing */
      SCIPclockStop(benders->setuptime, set);
   }

   return SCIP_OKAY;
}

/** informs variable benders that the branch and bound process data is being freed */
SCIP_RETCODE SCIPbendersExitsol(
   SCIP_BENDERS*         benders,            /**< variable benders */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(benders != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of variable benders */
   if( benders->bendersexitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(benders->setuptime, set);

      SCIP_CALL( benders->bendersexitsol(set->scip, benders) );

      /* stop timing */
      SCIPclockStop(benders->setuptime, set);
   }

   return SCIP_OKAY;
}

/** solves the subproblem using the current master problem solution. */
SCIP_RETCODE SCIPbendersExec(
   SCIP_BENDERS*         benders,            /**< variable benders */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RESULT*          result              /**< result of the pricing process */
   )
{
   int i;

   assert(benders != NULL);
   assert(result != NULL);
   assert(benders->bendersexec != NULL);

   *result = SCIP_DIDNOTRUN;
   /* solving each of the subproblems for Benders decomposition */
   /* TODO: ensure that the each of the subproblems solve and update the parameters with the correct return values */
   for( i = 0; i < benders->nsubproblems; i++ )
      SCIP_CALL( benders->bendersexec(set->scip, benders, i) );

   /* TODO: generate cuts. A call back may be useful here. */

   return SCIP_OKAY;
}

/** returns the variable for the given variable of the subproblem. This provides a call back for the mapping between the
 * master and subproblems */
SCIP_VAR* SCIPbendersGetMasterVar(
   SCIP_BENDERS*         benders,            /**< variable benders */
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

/** gets user data of variable benders */
SCIP_BENDERSDATA* SCIPbendersGetData(
   SCIP_BENDERS*         benders             /**< variable benders */
   )
{
   assert(benders != NULL);

   return benders->bendersdata;
}

/** sets user data of variable benders; user has to free old data in advance! */
void SCIPbendersSetData(
   SCIP_BENDERS*         benders,            /**< variable benders */
   SCIP_BENDERSDATA*     bendersdata         /**< new variable benders user data */
   )
{
   assert(benders != NULL);

   benders->bendersdata = bendersdata;
}

/** sets copy callback of benders */
void SCIPbendersSetCopy(
   SCIP_BENDERS*         benders,            /**< variable benders */
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

/** sets master variables mapping function callback of the Benders' decomposition algorithm */
void SCIPbendersSetGetmastervar(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_DECL_BENDERSGETMASTERVAR((*bendersgetmastervar))/**< the master variable to subproblem var mapping function */
   )
{
   assert(benders != NULL);

   benders->bendersgetmastervar = bendersgetmastervar;
}

/** sets the subproblem solving callback of the Benders' decomposition algorithm */
void SCIPbendersSetExec(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_DECL_BENDERSEXEC ((*bendersexec))    /**< the subproblem solving method of the Benders' decomposition algorithm */
   )
{
   assert(benders != NULL);

   benders->bendersexec = bendersexec;
}

/** gets name of variable benders */
const char* SCIPbendersGetName(
   SCIP_BENDERS*         benders             /**< variable benders */
   )
{
   assert(benders != NULL);

   return benders->name;
}

/** gets description of variable benders */
const char* SCIPbendersGetDesc(
   SCIP_BENDERS*         benders             /**< variable benders */
   )
{
   assert(benders != NULL);

   return benders->desc;
}

/** gets priority of variable benders */
int SCIPbendersGetPriority(
   SCIP_BENDERS*         benders             /**< variable benders */
   )
{
   assert(benders != NULL);

   return benders->priority;
}

/** sets priority of variable benders */
void SCIPbendersSetPriority(
   SCIP_BENDERS*         benders,            /**< variable benders */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the variable benders */
   )
{
   assert(benders != NULL);
   assert(set != NULL);

   benders->priority = priority;
   set->benderssorted = FALSE;
}

/** gets the number of times, the benders was called and tried to find a variable with negative reduced costs */
int SCIPbendersGetNCalls(
   SCIP_BENDERS*         benders             /**< variable benders */
   )
{
   assert(benders != NULL);

   return benders->ncalls;
}

/** gets the number of optimality cuts found by the collection of Benders' decomposition subproblems */
int SCIPbendersGetNOptCutsFound(
   SCIP_BENDERS*         benders             /**< variable benders */
   )
{
   assert(benders != NULL);

   return benders->noptcutsfound;
}

/** gets the number of feasibility cuts found by the collection of Benders' decomposition subproblems */
int SCIPbendersGetNFeasCutsFound(
   SCIP_BENDERS*         benders             /**< variable benders */
   )
{
   assert(benders != NULL);

   return benders->nfeascutsfound;
}

/** gets time in seconds used in this benders for setting up for next stages */
SCIP_Real SCIPbendersGetSetupTime(
   SCIP_BENDERS*         benders             /**< variable benders */
   )
{
   assert(benders != NULL);

   return SCIPclockGetTime(benders->setuptime);
}

/** gets time in seconds used in this benders */
SCIP_Real SCIPbendersGetTime(
   SCIP_BENDERS*         benders             /**< variable benders */
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

/** is variable benders initialized? */
SCIP_Bool SCIPbendersIsInitialized(
   SCIP_BENDERS*         benders             /**< variable benders */
   )
{
   assert(benders != NULL);

   return benders->initialized;
}



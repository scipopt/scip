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
//#define SCIP_DEBUG
//#define SCIP_MOREDEBUG
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
#include "scip/cons_linear.h"

#include "scip/struct_benders.h"

/* Defaults for parameters */
#define SCIP_DEFAULT_TRANSFERCUTS          TRUE  /** Should Benders' cuts generated in LNS heuristics be transferred to the main SCIP instance? */
#define SCIP_DEFAULT_CUTSASCONSS           TRUE  /** Should the transferred cuts be added as constraints? */
#define SCIP_DEFAULT_MIPCHECKFREQ             5  /** the number of iterations that the MIP is checked, -1 for always. */
#define SCIP_DEFAULT_LNSCHECK              TRUE  /** should the Benders' decomposition be used in LNS heuristics */
#define SCIP_DEFAULT_LNSMAXDEPTH             -1  /** the maximum depth at which the LNS check is performed */
#define SCIP_DEFAULT_SUBPROBFRAC            1.0  /** the fraction of subproblems that are solved in each iteration */

#define AUXILIARYVAR_NAME     "##bendersauxiliaryvar" /** the name for the Benders' auxiliary variables in the master problem */

/* Local methods */

/** A workaround for GCG. This is a temp vardata that is set for the auxiliary variables */
struct SCIP_VarData
{
   int                   vartype;             /**< the variable type. In GCG this indicates whether the variable is a
                                                   master problem or subproblem variable. */
};


/** adds the auxiliary variables to the Benders' decomposition master problem */
static
SCIP_RETCODE addAuxiliaryVariablesToMaster(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< benders */
   )
{
   SCIP_VAR* auxiliaryvar;
   SCIP_VARDATA* vardata;
   char varname[SCIP_MAXSTRLEN];    /* the name of the auxiliary variable */
   int i;

   /* this is a workaround for GCG. GCG expects that the variable has vardata when added. So a dummy vardata is created */
   SCIP_CALL( SCIPallocBlockMemory(scip, &vardata) );
   vardata->vartype = -1;

   for( i = 0; i < SCIPbendersGetNSubproblems(benders); i++ )
   {
      /* if no optimality cuts have been added for this subproblem, then the auxiliary variable will be created and
       * added */
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s_%d", AUXILIARYVAR_NAME, i );
      SCIP_CALL( SCIPcreateVarBasic(scip, &auxiliaryvar, varname, -SCIPinfinity(scip), SCIPinfinity(scip), 1.0,
            SCIP_VARTYPE_CONTINUOUS) );

      SCIPvarSetData(auxiliaryvar, vardata);

      SCIP_CALL( SCIPaddVar(scip, auxiliaryvar) );

      benders->auxiliaryvars[i] = auxiliaryvar;

      SCIP_CALL( SCIPreleaseVar(scip, &auxiliaryvar) );
   }

   SCIPfreeBlockMemory(scip, &vardata);

   return SCIP_OKAY;
}

/* sets the subproblem objective value array to -infinity */
static
void resetSubproblemObjectiveValue(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition structure */
   )
{
   SCIP* subproblem;
   int nsubproblems;
   int i;

   assert(benders != NULL);

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   for( i = 0; i < nsubproblems; i++ )
   {
      subproblem = SCIPbendersSubproblem(benders, i);
      SCIPbendersSetSubprobObjval(benders, SCIPinfinity(subproblem), i);
   }
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
   SCIPsetBendersPriority(scip, (SCIP_BENDERS*)paramdata, SCIPparamGetInt(param)); /*lint !e740*/

   return SCIP_OKAY;
}

/** copies the given benders to a new scip */
SCIP_RETCODE SCIPbendersCopyInclude(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_SET*             sourceset,          /**< SCIP_SET of SCIP to copy from */
   SCIP_SET*             targetset,          /**< SCIP_SET of SCIP to copy to */
   SCIP_Bool*            valid               /**< was the copying process valid? */
   )
{
   SCIP_BENDERS* targetbenders;  /* the copy of the Benders' decomposition struct in the target set */

   assert(benders != NULL);
   assert(targetset != NULL);
   assert(valid != NULL);
   assert(targetset->scip != NULL);

   (*valid) = FALSE;

   if( benders->benderscopy != NULL && targetset->benders_copybenders )
   {
      SCIPsetDebugMsg(targetset, "including benders %s in subscip %p\n", SCIPbendersGetName(benders), (void*)targetset->scip);
      SCIP_CALL( benders->benderscopy(targetset->scip, benders, valid) );

      /* if the copy was valid, then the Benders cuts are copied. */
      if( (*valid) )
      {
         targetbenders = SCIPsetFindBenders(targetset, SCIPbendersGetName(benders));

         /* storing the pointer to the source scip instance */
         targetbenders->sourcescip = sourceset->scip;

         /* the flag is set to indicate that the Benders' decomposition is a copy */
         targetbenders->iscopy = TRUE;
      }
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
   SCIP_DECL_BENDERSGETVAR((*bendersgetvar)),/**< returns the master variable for a given subproblem variable */
   SCIP_DECL_BENDERSCREATESUB((*benderscreatesub)),/**< creates a Benders' decomposition subproblem */
   SCIP_DECL_BENDERSPRESUBSOLVE((*benderspresubsolve)),/**< called prior to the subproblem solving loop */
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

   /* Checking whether the benderssolvesub and the bendersfreesub are both implemented or both are not implemented */
   if( (benderssolvesub == NULL && bendersfreesub != NULL) || (benderssolvesub != NULL && bendersfreesub == NULL) )
   {
      SCIPerrorMessage("Benders' decomposition <%s> requires that both bendersSolvesub%s and bendersFreesub%s are \
         implemented or neither\n", name, name, name);
      return SCIP_INVALIDCALL;
   }

   SCIP_ALLOC( BMSallocMemory(benders) );
   BMSclearMemory(benders);
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*benders)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*benders)->desc, desc, strlen(desc)+1) );
   (*benders)->priority = priority;
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
   (*benders)->bendersgetvar = bendersgetvar;
   (*benders)->benderscreatesub = benderscreatesub;
   (*benders)->benderspresubsolve = benderspresubsolve;
   (*benders)->benderssolvesub = benderssolvesub;
   (*benders)->benderspostsolve = benderspostsolve;
   (*benders)->bendersfreesub = bendersfreesub;
   (*benders)->bendersdata = bendersdata;
   SCIP_CALL( SCIPclockCreate(&(*benders)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*benders)->bendersclock, SCIP_CLOCKTYPE_DEFAULT) );

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of benders <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*benders)->priority, FALSE, priority, INT_MIN/4, INT_MAX/4,
                  paramChgdBendersPriority, (SCIP_PARAMDATA*)(*benders)) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutlp", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
        "should Benders' cuts be generated for LP solutions?", &(*benders)->cutlp, FALSE, cutlp, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutpseudo", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
        "should Benders' cuts be generated for pseudo solutions?", &(*benders)->cutpseudo, FALSE, cutpseudo, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutrelax", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
        "should Benders' cuts be generated for relaxation solutions?", &(*benders)->cutrelax, FALSE, cutrelax, NULL, NULL) ); /*lint !e740*/

   /* These parameters are left for the user to decide in a settings file. This departs from the usual SCIP convention
    * where the settings available at the creation of the plugin can be set in the function call. */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/transfercuts", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
        "Should Benders' cuts from LNS heuristics be transferred to the main SCIP instance?", &(*benders)->transfercuts,
        FALSE, SCIP_DEFAULT_TRANSFERCUTS, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/mipcheckfreq", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
        "The frequency at which the MIP subproblems are checked, -1 for always", &(*benders)->mipcheckfreq, FALSE,
        SCIP_DEFAULT_MIPCHECKFREQ, -1, INT_MAX, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/lnscheck", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
        "Should Benders' decomposition be used in LNS heurisics?", &(*benders)->lnscheck, FALSE, SCIP_DEFAULT_LNSCHECK,
        NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/lnsmaxdepth", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
        "maximal depth level at which the LNS check is performed (-1: no limit)", &(*benders)->lnsmaxdepth, TRUE,
        SCIP_DEFAULT_LNSMAXDEPTH, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutsasconss", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
        "Should the transferred cuts be added as constraints?", &(*benders)->cutsasconss, FALSE,
        SCIP_DEFAULT_CUTSASCONSS, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/subprobfrac", name);
   SCIP_CALL( SCIPsetAddRealParam(set, messagehdlr, blkmem, paramname,
        "The fraction of subproblems that are solved in each iteration.", &(*benders)->subprobfrac, FALSE,
        SCIP_DEFAULT_SUBPROBFRAC, 0.0, 1.0, NULL, NULL) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** calls destructor and frees memory of Benders' decomposition */
SCIP_RETCODE SCIPbendersFree(
   SCIP_BENDERS**        benders,            /**< pointer to Benders' decomposition data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(benders != NULL);
   assert(*benders != NULL);
   assert(!(*benders)->initialized);
   assert(set != NULL);

   /* call destructor of Benders' decomposition */
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


/** creates the subproblems and registers it with the Benders' decomposition struct */
static
SCIP_RETCODE createSubproblems(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int nsubproblems;
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   /* if the subproblems have already been created, then they will not be created again. This is the case if the
    * transformed problem has been freed and then retransformed. The subproblems should only be created when the problem
    * is first transformed. */
   if( benders->subprobscreated )
      return SCIP_OKAY;

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   /* creating all subproblems */
   for( i = 0; i < nsubproblems; i++ )
   {
      /* calling the create subproblem call back method */
      SCIP_CALL( benders->benderscreatesub(set->scip, benders, i) );
   }

   benders->subprobscreated = TRUE;

   return SCIP_OKAY;

}


/** initializes Benders' decomposition */
SCIP_RETCODE SCIPbendersInit(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
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
      benders->ncutsfound = 0;
      benders->ntransferred = 0;
   }

   /* start timing */
   SCIPclockStart(benders->setuptime, set);

   /* creates the subproblems and sets up the probing mode for LP subproblems. This function calls the benderscreatesub
    * callback. */
   SCIP_CALL( createSubproblems(benders, set) );

   if( benders->bendersinit != NULL )
      SCIP_CALL( benders->bendersinit(set->scip, benders) );

   benders->initialized = TRUE;

   /* stop timing */
   SCIPclockStop(benders->setuptime, set);

   return SCIP_OKAY;
}


/** calls exit method of Benders' decomposition */
SCIP_RETCODE SCIPbendersExit(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(benders != NULL);
   assert(set != NULL);

   if( !benders->initialized )
   {
      SCIPerrorMessage("Benders' decomposition <%s> not initialized\n", benders->name);
      return SCIP_INVALIDCALL;
   }

   /* start timing */
   SCIPclockStart(benders->setuptime, set);

   if( benders->bendersexit != NULL )
      SCIP_CALL( benders->bendersexit(set->scip, benders) );

   benders->initialized = FALSE;

   /* stop timing */
   SCIPclockStop(benders->setuptime, set);

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

   if( !benders->iscopy )
   {
      /* adding the auxiliary variables to the master problem */
      SCIP_CALL( addAuxiliaryVariablesToMaster(set->scip, benders) );
   }

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

   return SCIP_OKAY;
}

/** informs Benders' decomposition that the branch and bound process data is being freed */
SCIP_RETCODE SCIPbendersExitsol(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
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

   return SCIP_OKAY;
}

/** activates benders such that it is called in LP solving loop */
SCIP_RETCODE SCIPbendersActivate(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nsubproblems        /**< the number subproblems used in this decomposition */
   )
{
   int i;

   assert(benders != NULL);
   assert(set != NULL);
   assert(set->stage == SCIP_STAGE_INIT || set->stage == SCIP_STAGE_PROBLEM);

   if( !benders->active )
   {
      benders->active = TRUE;
      set->nactivebenders++;
      set->benderssorted = FALSE;

      benders->nsubproblems = nsubproblems;

      /* allocating memory for the subproblems arrays */
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subproblems, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->auxiliaryvars, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subprobobjval, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->bestsubprobobjval, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subprobislp, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->mastervarscont, benders->nsubproblems) );

      for( i = 0; i < benders->nsubproblems; i++ )
      {
         benders->subproblems[i] = NULL;
         benders->auxiliaryvars[i] = NULL;
         benders->subprobobjval[i] = SCIPsetInfinity(set);
         benders->bestsubprobobjval[i] = SCIPsetInfinity(set);
         benders->subprobislp[i] = FALSE;
         benders->mastervarscont[i] = FALSE;
      }
   }

   return SCIP_OKAY;
}

/** deactivates benders such that it is no longer called in LP solving loop */
void SCIPbendersDeactivate(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(benders != NULL);
   assert(set != NULL);
   assert(set->stage == SCIP_STAGE_INIT || set->stage == SCIP_STAGE_PROBLEM);

   if( benders->active )
   {
#ifndef NDEBUG
      /* checking whether the auxiliary variables and subproblems are all NULL */
      int nsubproblems;
      int i;

      nsubproblems = SCIPbendersGetNSubproblems(benders);

      for( i = 0; i < nsubproblems; i++ )
      {
         assert(benders->auxiliaryvars[i] == NULL);
         assert(benders->subproblems[i] == NULL);
      }
#endif

      benders->active = FALSE;
      set->nactivebenders--;
      set->benderssorted = FALSE;

      /* freeing the memory allocated during the activation of the Benders' decomposition */
      BMSfreeMemoryArray(&benders->mastervarscont);
      BMSfreeMemoryArray(&benders->subprobislp);
      BMSfreeMemoryArray(&benders->bestsubprobobjval);
      BMSfreeMemoryArray(&benders->subprobobjval);
      BMSfreeMemoryArray(&benders->auxiliaryvars);
      BMSfreeMemoryArray(&benders->subproblems);
   }
}

/** returns whether the given Benders decomposition is in use in the current problem */
SCIP_Bool SCIPbendersIsActive(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition structure */
   )
{
   assert(benders != NULL);

   return benders->active;
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
   SCIP_Bool*            infeasible,         /**< is the master problem infeasible with respect to the Benders' cuts? */
   SCIP_BENDERSENFOTYPE  type,               /**< the type of solution being enforced */
   SCIP_Bool             checkint            /**< should the integer solution be checked by the subproblems */
   )
{
   int nsubproblems;
   int subproblemcount;
   int subprobssolved;
   int nbenderscuts;
   int nsolveloops;     /* the number of times the subproblems are solved. An additional loop is required when the subproblem in integer */
   int numnotopt;
   int numtocheck;
   int i;
   int j;
   int l;
   SCIP_Bool optimal;
   SCIP_Bool allchecked;      /* flag to indicate whether all subproblems have been checked */
   int nchecked;              /* the number of subproblems that have been checked */
   SCIP_Bool* subisinfeas;
   SCIP_Bool onlylpcheck;     /* should only the LP be checked in the presence of integer subproblems */

   /* start timing */
   SCIPclockStart(benders->bendersclock, set);

   (*infeasible) = FALSE;

   /* It is assumed that the problem is optimal, until a subproblem is found not to be optimal. However, not all
    * subproblems could be checked in each iteration. As such, it is not possible to state that the problem is optimal
    * if not all subproblems are checked. Situations where this may occur is when a subproblem is a MIP and only the LP
    * is solved. Also, in a distributed computation, then it may be adventageous to only solve some subproblems before
    * resolving the master problem. As such, for a problem to be optimal, then (optimal && allchecked) == TRUE */
   optimal = TRUE;
   nchecked = 0;

   assert(benders != NULL);
   assert(result != NULL);

   /* if the Benders' decomposition is called from a sub-scip, it is assumed that this is an LNS heuristic. As such, the
    * check is not performed and the solution is assumed to be feasible */
   if( benders->iscopy
      && (!benders->lnscheck
         || (benders->lnsmaxdepth > -1 && SCIPgetDepth(benders->sourcescip) > benders->lnsmaxdepth)) )
   {
      (*result) = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

   /* it is not necessary to check all primal solutions by solving the Benders' decomposition subproblems.
    * Only the improving solutions are checked to improve efficiency of the algorithm.
    * If the solution is non-improving, the result FEASIBLE is returned. While this may be incorrect w.r.t to the
    * Benders' subproblems, this solution will never be the optimal solution. A non-improving solution may be used
    * within LNS primal heuristics. If this occurs, the improving solution, if found, will be checked by the solving
    * the Benders' decomposition subproblems.
    * TODO: Add a parameter to control this behaviour.*/
   if( checkint && SCIPsetIsFeasLE(set, SCIPgetPrimalbound(set->scip), SCIPgetSolOrigObj(set->scip, sol)) )
   {
      (*result) = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

#if 0
   /* This is part of testing!!!
    * A flag `onlylpcheck` has been added to force only the LP's to be checked during a LNS heuristic solve. As such,
    * the integer subproblems need to be checked after the completion of the LNS heuristic.
    * If the integer subproblems are checked during the LNS solve, then the Benders' subproblem solves would not be
    * required for the solutions generated from the LNS heuristic. */
   if( sol != NULL && SCIPsolGetHeur(sol) != NULL && strcmp(SCIPheurGetName(SCIPsolGetHeur(sol)), HEUR_RINS) == 0 && benders->lnscheck )
   {
      (*result) = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }
#endif

   /* when Benders' is used in the LNS heuristics, only the LP of the master/subproblems is checked, i.e. no integer
    * cuts are generated. This assumes that in the LNS heuristics the subproblems are all LPs. */
   onlylpcheck = benders->iscopy && benders->lnscheck;

   nsubproblems = SCIPbendersGetNSubproblems(benders);
   if( benders->ncalls == 0 || type == SCIP_BENDERSENFOTYPE_CHECK || onlylpcheck )
      numtocheck = nsubproblems;
   else
      numtocheck = (int) SCIPsetCeil(set, (SCIP_Real) nsubproblems*benders->subprobfrac);
   benders->firstchecked = benders->lastchecked;

   /* allocating memory for the infeasible subproblem array */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &subisinfeas, nsubproblems) );

   for( i = 0; i < nsubproblems; i++ )
      subisinfeas[i] = FALSE;

   /* sets the stored objective function values of the subproblems to infinity */
   resetSubproblemObjectiveValue(benders);

   SCIPdebugMessage("Performing the subproblem solving process. Number of subproblems to check %d\n", numtocheck);

   if( benders->benderspresubsolve != NULL )
      SCIP_CALL( benders->benderspresubsolve(set->scip, benders) );

   *result = SCIP_DIDNOTRUN;

   /* by default the number of solve loops is 1. This is the case is all subproblems are LP or the user has defined a
    * benderssolvesub callback. If there is a subproblem that is not an LP, then 2 solve loops are performed. The first
    * loop is the LP solving loop, the second solves the subproblem to integer optimality. */
   nsolveloops = 1;

   /* the result flag is set to FEASIBLE as default. If a cut is added, then this is changed to CONS_ADDED. If the
    * problem is found to be infeasible due to the auxiliary variable values, then this is changed to INFEASIBLE. */
   (*result) = SCIP_FEASIBLE;

   for( l = 0; l < nsolveloops; l++ )
   {
      SCIPdebugMessage("Benders' decomposition - solve loop %d\n", l);
      numnotopt = 0;
      subproblemcount = 0;

      if( type == SCIP_BENDERSENFOTYPE_CHECK && sol == NULL )
      {
         /* This if statement doesn't make sense to me. Need to determine what happens in the cut generation stage if
          * this is true */
         (*infeasible) = TRUE;
      }
      else
      {
         /* solving each of the subproblems for Benders decomposition */
         /* TODO: ensure that the each of the subproblems solve and update the parameters with the correct return values */
         i = benders->firstchecked;
         //for( i = 0; i < nsubproblems; i++ )
         while( subproblemcount < nsubproblems && numnotopt < numtocheck )
         {
            SCIP_Bool subinfeas = FALSE;
            SCIP_Bool lpsub = SCIPbendersSubprobIsLP(benders, i);
            SCIP_Bool solvesub = TRUE;

            /* for the second solving loop, if the problem is an LP, it is not solved again. If the problem is a MIP,
             * then the subproblem objective function value is set to infinity. However, if the subproblem is proven
             * infeasible from the LP, then the IP loop is not performed. */
            if( l > 0 )
            {
               if( lpsub || subisinfeas[i] )
                  solvesub = FALSE;
               else
                  SCIPbendersSetSubprobObjval(benders, SCIPinfinity(SCIPbendersSubproblem(benders, i)), i);
            }

            if( solvesub )
            {
               SCIP_CALL( SCIPbendersExecSubproblemSolve(benders, set, sol, i, l, FALSE, &subinfeas, type) );

#ifdef SCIP_DEBUG
               if( type == LP )
               {
                  SCIPdebugMessage("LP: Subproblem %d (%f < %f)\n", i, SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i),
                     SCIPbendersGetSubprobObjval(benders, i));
               }
#endif
               (*infeasible) = (*infeasible) || subinfeas;
               subisinfeas[i] = subinfeas;

               /* if the subproblems are being solved as part of the conscheck, then we break once an infeasibility is found.
                * The result pointer is set to (*infeasible) and the execution is halted. */
               if( checkint )
               {
                  /* if the subproblem is feasible, then it is necessary to update the value of the auxiliary variable to the
                   * objective function value of the subproblem. */
                  if( !subinfeas )
                  {
                     SCIP_Bool subproboptimal;

                     SCIP_CALL( SCIPbendersCheckAuxiliaryVar(benders, set, sol, i, &subproboptimal) );

                     if( lpsub || benders->benderssolvesub != NULL || l > 0 || onlylpcheck )
                        optimal = optimal && subproboptimal;

#ifdef SCIP_DEBUG
                     if( lpsub || l > 0 )
                     {
                        if( subproboptimal )
                        {
                           SCIPdebugMessage("Subproblem %d is Optimal (%f >= %f)\n", i,
                              SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i), SCIPbendersGetSubprobObjval(benders, i));
                        }
                        else
                        {
                           SCIPdebugMessage("Subproblem %d is NOT Optimal (%f < %f)\n", i,
                              SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i), SCIPbendersGetSubprobObjval(benders, i));
                        }
                     }
#endif

                     /* only increment the checked count if the subproblem is not an LP, or the solve loop is the MIP
                      * solving loop. Hence, the LP are solved once and the MIPs are solved twice */
                     if( (l == 0 && lpsub) || (l > 0 && !lpsub) || onlylpcheck )
                        nchecked++;


                     if( !subproboptimal )
                     {
                        numnotopt++;
                        assert(numnotopt <= nsubproblems);
                     }
                  }
                  else
                  {
                     numnotopt++;
                     assert(numnotopt <= nsubproblems);
                  }
               }
            }

            subproblemcount++;
            i++;
            if( i >= nsubproblems )
               i = 0;
            benders->lastchecked = i;
         }
      }

      subprobssolved = subproblemcount;


      /* Generating cuts for the subproblems. */
      /* TODO: The cut generating loop will be added for later merge requests */
   }

   allchecked = (nchecked == nsubproblems);

#ifndef NDEBUG
   if( (*result) == SCIP_CONSADDED )
   {
      SCIPdebugMessage("Benders decomposition: Cut added\n");
   }
#endif

   if( checkint && (type == SCIP_BENDERSENFOTYPE_CHECK || (*result) != SCIP_CONSADDED) )
   {
      /* if the subproblems are being solved as part of conscheck, then the results flag must be returned after the solving
       * has completed. */
      if( (*infeasible) || !allchecked )
         (*result) = SCIP_INFEASIBLE;
      else
      {
         (*result) = SCIP_FEASIBLE;

         /* if the subproblems are not infeasible, but they are also not optimal, then the check return the result of
          * feasible and the flag of infeasible. */
         (*infeasible) = !optimal;
      }
   }
#if 0
   /* The else branch existed when the if statement checked for type == CHECK. I think that the else is not needed.
    * Keeping it here to check whether it was needed  */
   else if( type == PSEUDO )
   {
      if( (*infeasible) || !(optimal && allchecked) )
         (*result) = SCIP_INFEASIBLE;
      else
         (*result) = SCIP_FEASIBLE;
   }
#endif

   /* calling the post-solve call back for the Benders' decomposition algorithm. This allows the user to work directly
    * with the solved subproblems and the master problem */
   if( benders->benderspostsolve != NULL )
   {
      printf("checkint %d result %d\n", checkint, (*result));
      SCIP_CALL( benders->benderspostsolve(set->scip, benders, sol, (*infeasible)) );
   }

   /* freeing the subproblems after the cuts are generated */
   i = benders->firstchecked;
   subproblemcount = 0;
   while( subproblemcount < subprobssolved )
      //for( i = 0; i < benders->nsubproblems; i++ )
   {
      SCIP_CALL( SCIPbendersFreeSubproblem(benders, set, i) );

      subproblemcount++;
      i++;
      if( i >= nsubproblems )
         i = 0;
   }

   /* increment the number of calls to the Benders' decomposition subproblem solve */
   benders->ncalls++;

   /* end timing */
   SCIPclockStop(benders->bendersclock, set);

   /* freeing memory */
   SCIPsetFreeBufferArray(set, &subisinfeas);

   return SCIP_OKAY;
}

/** solves the subproblems. */
SCIP_RETCODE SCIPbendersExecSubproblemSolve(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnum,            /**< the subproblem number */
   int                   solveloop,          /**< the solve loop iteration. The first iter is for LP, the second for IP */
   SCIP_Bool             enhancement,        /**< is the solve performed as part of and enhancement? */
   SCIP_Bool*            infeasible,         /**< returns whether the current subproblem is infeasible */
   SCIP_BENDERSENFOTYPE  type                /**< the enforcement type calling this function */
   )
{
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
   assert(benders->bendersfreesub != NULL || (benders->bendersfreesub == NULL && benders->benderssolvesub == NULL));
   assert(probnum >= 0 && probnum < benders->nsubproblems);

   if( benders->bendersfreesub != NULL )
      SCIP_CALL( benders->bendersfreesub(set->scip, benders, probnum) );
   else
   {
   }

   return SCIP_OKAY;
}

/** checks the auxiliary variable value for optimality */
SCIP_RETCODE SCIPbendersCheckAuxiliaryVar(
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber,         /**< the number of the pricing problem */
   SCIP_Bool*            optimal             /**< flag to indicate whether the current subproblem is optimal for the master */
   )
{
   SCIP_Real auxiliaryvarval;
   SCIP_Real soltol;

   assert(benders != NULL);
   assert(set != NULL);
   assert(probnumber >= 0 && probnumber < benders->nsubproblems);

   (*optimal) = FALSE;

   auxiliaryvarval = SCIPbendersGetAuxiliaryVarVal(benders, set, sol, probnumber);

   SCIP_CALL( SCIPsetGetRealParam(set, "benders/solutiontol", &soltol) );

   SCIPsetDebugMsg(set, "Subproblem %d - Auxiliary Variable: %g Subproblem Objective: %g\n", probnumber, auxiliaryvarval,
      SCIPbendersGetSubprobObjval(benders, probnumber));

   if( SCIPsetIsFeasGE(set, auxiliaryvarval + soltol, SCIPbendersGetSubprobObjval(benders, probnumber)) )
      (*optimal) = TRUE;

   return SCIP_OKAY;
}

/** returns the value of the auxiliary variable value in a master problem solution */
SCIP_Real SCIPbendersGetAuxiliaryVarVal(
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber          /**< the number of the pricing problem */
   )
{
   SCIP_VAR* auxiliaryvar;

   assert(benders != NULL);
   assert(set != NULL);

   auxiliaryvar = SCIPbendersGetAuxiliaryVar(benders, probnumber);
   assert(auxiliaryvar != NULL);

   return SCIPgetSolVal(set->scip, sol, auxiliaryvar);
}

/** returns the corresponding master or subproblem variable for the given variable.
 * This provides a call back for the variable mapping between the master and subproblems */
SCIP_RETCODE SCIPbendersGetVar(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< the variable for which the corresponding variable is desired */
   SCIP_VAR**            mappedvar,          /**< the variable that is mapped to var */
   int                   probnumber          /**< the problem number for the desired variable, -1 for the master problem */
   )
{
   assert(benders != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(mappedvar != NULL);
   assert(benders->bendersgetvar != NULL);

   (*mappedvar) = NULL;

   /* if the variable name matches the auxiliary variable, then the master variable is returned as NULL */
   if( strstr(SCIPvarGetName(var), AUXILIARYVAR_NAME) != NULL )
      return SCIP_OKAY;

   SCIP_CALL( benders->bendersgetvar(set->scip, benders, var, mappedvar, probnumber) );

   return SCIP_OKAY;
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

/** sets destructor callback of Benders' decomposition */
void SCIPbendersSetFree(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSFREE ((*bendersfree))    /**< destructor of Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersfree = bendersfree;
}

/** sets initialization callback of Benders' decomposition */
void SCIPbendersSetInit(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSINIT((*bendersinit))     /**< initialize the Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersinit = bendersinit;
}

/** sets deinitialization callback of Benders' decomposition */
void SCIPbendersSetExit(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSEXIT((*bendersexit))     /**< deinitialize the Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersexit = bendersexit;
}

/** sets presolving initialization callback of Benders' decomposition */
void SCIPbendersSetInitpre(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSINITPRE((*bendersinitpre))/**< initialize presolving for Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersinitpre = bendersinitpre;
}

/** sets presolving deinitialization callback of Benders' decomposition */
void SCIPbendersSetExitpre(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSEXITPRE((*bendersexitpre))/**< deinitialize presolving for Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersexitpre = bendersexitpre;
}

/** sets solving process initialization callback of Benders' decomposition */
void SCIPbendersSetInitsol(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol))/**< solving process initialization callback of Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersinitsol = bendersinitsol;
}

/** sets solving process deinitialization callback of Benders' decomposition */
void SCIPbendersSetExitsol(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol))/**< solving process deinitialization callback of Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersexitsol = bendersexitsol;
}

/** sets the pre subproblem solve callback of Benders' decomposition */
void SCIPbendersSetPresubsolve(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSPRESUBSOLVE((*benderspresubsolve))/**< called prior to the subproblem solving loop */
   )
{
   assert(benders != NULL);

   benders->benderspresubsolve = benderspresubsolve;
}

/** sets solve callback of Benders' decomposition */
void SCIPbendersSetSolvesub(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSSOLVESUB((*benderssolvesub))/**< solving method for a Benders' decomposition subproblem */
   )
{
   assert(benders != NULL);

   benders->benderssolvesub = benderssolvesub;
}

/** sets post-solve callback of Benders' decomposition */
void SCIPbendersSetPostsolve(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSPOSTSOLVE((*benderspostsolve))/**< solving process deinitialization callback of Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->benderspostsolve = benderspostsolve;
}

/** sets free subproblem callback of Benders' decomposition */
void SCIPbendersSetFreesub(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSFREESUB((*bendersfreesub))/**< the freeing callback for the subproblem */
   )
{
   assert(benders != NULL);

   benders->bendersfreesub = bendersfreesub;
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
int SCIPbendersGetNCutsFound(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->ncutsfound;
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

/** returns all auxiliary variables */
SCIP_VAR** SCIPbendersGetAuxiliaryVars(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->auxiliaryvars;
}
/** stores the objective function value of the subproblem for use in cut generation */
void SCIPbendersSetSubprobObjval(
   SCIP_BENDERS*         benders,            /**< variable benders */
   SCIP_Real             objval,             /**< the objective function value for the subproblem */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   /* updating the best objval */
   if( objval < benders->bestsubprobobjval[probnumber] )
      benders->bestsubprobobjval[probnumber] = objval;

   benders->subprobobjval[probnumber] = objval;
}

/** returns the objective function value of the subproblem for use in cut generation */
SCIP_Real SCIPbendersGetSubprobObjval(
   SCIP_BENDERS*         benders,            /**< variable benders */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->subprobobjval[probnumber];
}

/* returns whether the subproblem is an LP. This means that the dual solution can be trusted. */
SCIP_Bool SCIPbendersSubprobIsLP(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->subprobislp[probnumber];
}

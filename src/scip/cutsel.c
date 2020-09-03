/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cutsel.c
 * @ingroup OTHER_CFILES
 * @brief  methods for cut selectors
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

//TODO: remove unused headers
#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/clock.h"
#include "scip/stat.h"
#include "scip/visual.h"
#include "scip/paramset.h"
#include "scip/tree.h"
#include "scip/reopt.h"
#include "scip/lp.h"
#include "scip/scip.h"
#include "scip/cutsel.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"

#include "scip/struct_cutsel.h"
#include "scip/struct_scip.h"

/** internal method for creating a cut selector */
static
SCIP_RETCODE doCutselCreate(
   SCIP_CUTSEL**         cutsel,             /**< pointer to store cut selector */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of cut selector */
   const char*           desc,               /**< description of cut selector */
   int                   priority,           /**< priority of the cut selector */
   SCIP_DECL_CUTSELCOPY ((*cutselcopy)),     /**< copy method of cut selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CUTSELFREE ((*cutselfree)),     /**< destructor of cut selector */
   SCIP_DECL_CUTSELINIT ((*cutselinit)),     /**< initialize cut selector */
   SCIP_DECL_CUTSELEXIT ((*cutselexit)),     /**< deinitialize cut selector */
   SCIP_DECL_CUTSELINITSOL((*cutselinitsol)),/**< solving process initialization method of cut selector */
   SCIP_DECL_CUTSELEXITSOL((*cutselexitsol)),/**< solving process deinitialization method of cut selector */
   SCIP_DECL_CUTSELSELECT((*cutselselect)),  /**< cut selection method */
   SCIP_CUTSELDATA*     cutseldata           /**< cut selector data */
   )
{
   //char paramname[SCIP_MAXSTRLEN];
   //char paramdesc[SCIP_MAXSTRLEN];

   assert(cutsel != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(cutselselect != NULL);

   SCIP_ALLOC( BMSallocMemory(cutsel) );
   BMSclearMemory(*cutsel);

   SCIP_ALLOC( BMSduplicateMemoryArray(&(*cutsel)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*cutsel)->desc, desc, strlen(desc)+1) );
   (*cutsel)->priority = priority;
   (*cutsel)->cutselcopy = cutselcopy;
   (*cutsel)->cutselfree = cutselfree;
   (*cutsel)->cutselinit = cutselinit;
   (*cutsel)->cutselexit = cutselexit;
   (*cutsel)->cutselinitsol = cutselinitsol;
   (*cutsel)->cutselexitsol = cutselexitsol;
   (*cutsel)->cutselselect = cutselselect;
   (*cutsel)->cutseldata = cutseldata;
   (*cutsel)->initialized = FALSE;
   /* create clocks */
   SCIP_CALL( SCIPclockCreate(&(*cutsel)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*cutsel)->cutseltime, SCIP_CLOCKTYPE_DEFAULT) );

   // TODO
   ///* add parameters */
   //(void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "cutselection/%s/priority", name);
   //(void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of cut selection rule <%s>", name);
   //SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
   //               &(*cutsel)->priority, FALSE, priority, INT_MIN/4, INT_MAX/2,
   //               paramChgdCutselPriority, (SCIP_PARAMDATA*)(*cutsel)) ); /*lint !e740*/

   return SCIP_OKAY;
}


/** creates a cut selector */
SCIP_RETCODE SCIPcutselCreate(
   SCIP_CUTSEL**         cutsel,             /**< pointer to store cut selector */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of cut selector */
   const char*           desc,               /**< description of cut selector */
   int                   priority,           /**< priority of the cut selector in standard mode */
   SCIP_DECL_CUTSELCOPY ((*cutselcopy)),     /**< copy method of cut selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CUTSELFREE ((*cutselfree)),     /**< destructor of cut selector */
   SCIP_DECL_CUTSELINIT ((*cutselinit)),     /**< initialize cut selector */
   SCIP_DECL_CUTSELEXIT ((*cutselexit)),     /**< deinitialize cut selector */
   SCIP_DECL_CUTSELINITSOL((*cutselinitsol)),/**< solving process initialization method of cut selector */
   SCIP_DECL_CUTSELEXITSOL((*cutselexitsol)),/**< solving process deinitialization method of cut selector */
   SCIP_DECL_CUTSELSELECT((*cutselselect)),  /**< cut selection method */
   SCIP_CUTSELDATA*     cutseldata           /**< cut selector data */
   )
{
   assert(cutsel != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(cutselselect != NULL);

   SCIP_CALL_FINALLY( doCutselCreate(cutsel, set, messagehdlr, blkmem, name, desc, priority,
      cutselcopy, cutselfree, cutselinit, cutselexit, cutselinitsol, cutselexitsol, cutselselect,
      cutseldata), /* TODO: (void) SCIPcutselFree(cutsel, set)*/ NULL );

   return SCIP_OKAY;
}

/** gets name of cut selector */
const char* SCIPcutselGetName(
   SCIP_CUTSEL*         cutsel             /**< cut selector */
   )
{
   assert(cutsel != NULL);

   return cutsel->name;
}

/** calls cut selectors to select cuts */
SCIP_RETCODE SCIPcutselsSelect(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW**            cuts,               /**< array with cuts to select from */
   int                   ncuts,              /**< length of cuts */
   int                   nforcedcuts,        /**< number of forced cuts at start of given array */
   SCIP_Bool             root,               /**< are we at the root node? */
   int                   maxnselectedcuts,   /**< maximum number of cuts to be selected */
   int*                  nselectedcuts       /**< pointer to return number of selected cuts */
   )
{
   int i;
   SCIP_RESULT result = SCIP_DIDNOTFIND;

   /* sort the cut selectors by priority */
   //TODO: SCIPsetSortCutselectors(set);

   printf("We have <%d> cut selectors! \n", set->ncutsels);

   /* try all cut selectors until one succeeds */
   for( i = 0; i < set->ncutsels && result == SCIP_DIDNOTFIND; ++i )
   {
      SCIP_CUTSEL* cutsel;
      cutsel = set->cutsels[i];

      SCIP_CALL( cutsel->cutselselect(set->scip, cutsel, cuts, ncuts, nforcedcuts, root, maxnselectedcuts,
               nselectedcuts, &result) );
      assert(*nselectedcuts <= maxnselectedcuts);
      assert(result == SCIP_SUCCESS || result == SCIP_DIDNOTFIND);
   }

   return SCIP_OKAY;
}

/** copies the given cut selector to a new scip */
SCIP_RETCODE SCIPcutselCopyInclude(
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
    assert(cutsel != NULL);
    assert(set != NULL);
    assert(set->scip != NULL);

    if( cutsel->cutselcopy != NULL )
    {
        SCIPsetDebugMsg(set, "including cut selector %s in subscip %p\n", SCIPcutselGetName(cutsel), (void*)set->scip);
        SCIP_CALL( cutsel->cutselcopy(set->scip, cutsel) );
    }
    return SCIP_OKAY;
}

/** sets copy method of cut selector */
void SCIPcutselSetCopy(
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_DECL_CUTSELCOPY  ((*cutselcopy))  /**< copy method of cut selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
    assert(cutsel != NULL);

    cutsel->cutselcopy = cutselcopy;
}

/** gets user data of cut selector */
SCIP_CUTSELDATA* SCIPcutselGetData(
   SCIP_CUTSEL*          cutsel              /**< cut selector */
   )
{
   assert(cutsel != NULL);

   return cutsel->cutseldata;
}

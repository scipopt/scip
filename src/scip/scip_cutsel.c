//TODO remove every appearance of the word 'node'
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

/**@file   scip_cutsel.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for cut selector plugins
 * @author Mark Turner
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/debug.h"
#include "scip/cutsel.h"
#include "scip/pub_message.h"
#include "scip/scip_cutsel.h"
#include "scip/set.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"

/** creates a cut selector and includes it in SCIP.
 *
 *  @note method has all cut selector callbacks as arguments and is thus changed every time a new
 *        callback is added in future releases; consider using SCIPincludeCutselBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_RETCODE SCIPincludeCutsel(
   SCIP*                 scip,               /**< SCIP data structure */
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
   SCIP_CUTSEL* cutsel;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeCutsel", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether cut selector is already present */
   if( SCIPfindCutsel(scip, name) != NULL )
   {
      SCIPerrorMessage("cut selector <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPcutselCreate(&cutsel, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, priority,
         cutselcopy, cutselfree, cutselinit, cutselexit, cutselinitsol, cutselexitsol,
         cutselselect, cutseldata) );
   SCIP_CALL( SCIPsetIncludeCutsel(scip->set, cutsel) );

   return SCIP_OKAY;
}

/** Creates a cut selector and includes it in SCIP with its most fundamental callbacks. All non-fundamental
 *  (or optional) callbacks as, e.g., init and exit callbacks, will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see SCIPsetCutselCopy(), SCIPsetCutselFree(),
 *  SCIPsetCutselInit(), SCIPsetCutselExit(), SCIPsetCutselInitsol(), and SCIPsetCutselExitsol()
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeCutsel() instead
 */
SCIP_RETCODE SCIPincludeCutselBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTSEL**         cutsel,             /**< reference to a cut selector, or NULL */
   const char*           name,               /**< name of cut selector */
   const char*           desc,               /**< description of cut selector */
   int                   priority,           /**< priority of the cut selector in standard mode */
   SCIP_DECL_CUTSELSELECT((*cutselselect)),  /**< cut selection method */
   SCIP_CUTSELDATA*      cutseldata          /**< cut selector data */
   )
{
    SCIP_CUTSEL* cutselptr;

    SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeCutselBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

    /* check whether cut selector is already present */
    if( SCIPfindCutsel(scip, name) != NULL )
    {
        SCIPerrorMessage("cut selector <%s> already included.\n", name);
        return SCIP_INVALIDDATA;
    }

    SCIP_CALL( SCIPcutselCreate(&cutselptr, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, priority,
                                 NULL, NULL, NULL, NULL, NULL, NULL,
                                 cutselselect, cutseldata) );
    SCIP_CALL( SCIPsetIncludeCutsel(scip->set, cutselptr) );

    if( cutsel != NULL )
        *cutsel = cutselptr;

    return SCIP_OKAY;
}

/** sets copy method of cut selector */
SCIP_RETCODE SCIPsetCutselCopy(
        SCIP*                 scip,               /**< SCIP data structure */
        SCIP_CUTSEL*          cutsel,             /**< cut selector */
        SCIP_DECL_CUTSELCOPY  ((*cutselcopy))  /**< copy method of cut selector or NULL if you don't want to copy your plugin into sub-SCIPs */
)
{
    SCIP_CALL( SCIPcheckStage(scip, "SCIPsetCutselCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

    assert(cutsel != NULL);

    SCIPcutselSetCopy(cutsel, cutselcopy);

    return SCIP_OKAY;
}


/** returns the cut selector of the given name, or NULL if not existing */
SCIP_CUTSEL* SCIPfindCutsel(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of node selector */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindCutsel(scip->set, name);
}


#ifdef TODO // adapt the code below
/** Creates a node selector and includes it in SCIP with its most fundamental callbacks. All non-fundamental
 *  (or optional) callbacks as, e.g., init and exit callbacks, will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see SCIPsetNodeselCopy(), SCIPsetNodeselFree(),
 *  SCIPsetNodeselInit(), SCIPsetNodeselExit(), SCIPsetNodeselInitsol(), and SCIPsetNodeselExitsol()
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeNodesel() instead
 */
SCIP_RETCODE SCIPincludeNodeselBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL**        nodesel,            /**< reference to a node selector, or NULL */
   const char*           name,               /**< name of node selector */
   const char*           desc,               /**< description of node selector */
   int                   stdpriority,        /**< priority of the node selector in standard mode */
   int                   memsavepriority,    /**< priority of the node selector in memory saving mode */
   SCIP_DECL_NODESELSELECT((*nodeselselect)),/**< node selection method */
   SCIP_DECL_NODESELCOMP ((*nodeselcomp)),   /**< node comparison method */
   SCIP_NODESELDATA*     nodeseldata         /**< node selector data */
   )
{
   SCIP_NODESEL* nodeselptr;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeNodeselBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether node selector is already present */
   if( SCIPfindNodesel(scip, name) != NULL )
   {
      SCIPerrorMessage("node selector <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPnodeselCreate(&nodeselptr, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, stdpriority, memsavepriority,
         NULL, NULL, NULL, NULL, NULL, NULL,
         nodeselselect, nodeselcomp, nodeseldata) );
   SCIP_CALL( SCIPsetIncludeNodesel(scip->set, nodeselptr) );

   if( nodesel != NULL )
      *nodesel = nodeselptr;

   return SCIP_OKAY;
}

/** sets copy method of node selector */
SCIP_RETCODE SCIPsetNodeselCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELCOPY ((*nodeselcopy))    /**< copy method of node selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNodeselCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(nodesel != NULL);

   SCIPnodeselSetCopy(nodesel, nodeselcopy);

   return SCIP_OKAY;
}

/** sets destructor method of node selector */
SCIP_RETCODE SCIPsetNodeselFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELFREE ((*nodeselfree))    /**< destructor of node selector */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNodeselFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(nodesel != NULL);

   SCIPnodeselSetFree(nodesel, nodeselfree);

   return SCIP_OKAY;
}

/** sets initialization method of node selector */
SCIP_RETCODE SCIPsetNodeselInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELINIT ((*nodeselinit))    /**< initialize node selector */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNodeselInit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(nodesel != NULL);

   SCIPnodeselSetInit(nodesel, nodeselinit);

   return SCIP_OKAY;
}

/** sets deinitialization method of node selector */
SCIP_RETCODE SCIPsetNodeselExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELEXIT ((*nodeselexit))    /**< deinitialize node selector */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNodeselExit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(nodesel != NULL);

   SCIPnodeselSetExit(nodesel, nodeselexit);

   return SCIP_OKAY;
}

/** sets solving process initialization method of node selector */
SCIP_RETCODE SCIPsetNodeselInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELINITSOL ((*nodeselinitsol))/**< solving process initialization method of node selector */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNodeselInitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(nodesel != NULL);

   SCIPnodeselSetInitsol(nodesel, nodeselinitsol);

   return SCIP_OKAY;
}

/** sets solving process deinitialization method of node selector */
SCIP_RETCODE SCIPsetNodeselExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELEXITSOL ((*nodeselexitsol))/**< solving process deinitialization method of node selector */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNodeselExitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(nodesel != NULL);

   SCIPnodeselSetExitsol(nodesel, nodeselexitsol);

   return SCIP_OKAY;
}

/** returns the array of currently available node selectors */
SCIP_NODESEL** SCIPgetNodesels(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nodesels;
}

/** returns the number of currently available node selectors */
int SCIPgetNNodesels(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nnodesels;
}

/** sets the priority of a node selector */
SCIP_RETCODE SCIPsetNodeselPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   int                   priority            /**< new priority of the node selector */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPnodeselSetPriority(nodesel, scip->set, priority);

   return SCIP_OKAY;
}

/** returns the currently used node selector */
SCIP_NODESEL* SCIPgetNodesel(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetGetNodesel(scip->set, scip->stat);
}
#endif

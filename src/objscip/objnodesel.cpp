/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: objnodesel.cpp,v 1.7 2005/01/18 09:26:49 bzfpfend Exp $"

/**@file   objnodesel.cpp
 * @brief  C++ wrapper for node selectors
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objnodesel.h"




/*
 * Data structures
 */

/** node selector data */
struct NodeselData
{
   scip::ObjNodesel* objnodesel;        /**< node selector object */
   Bool             deleteobject;       /**< should the node selector object be deleted when node selector is freed? */
};




/*
 * Callback methods of node selector
 */

/** destructor of node selector to free user data (called when SCIP is exiting) */
static
DECL_NODESELFREE(nodeselFreeObj)
{  /*lint --e{715}*/
   NODESELDATA* nodeseldata;

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);
   assert(nodeseldata->objnodesel != NULL);

   /* call virtual method of nodesel object */
   CHECK_OKAY( nodeseldata->objnodesel->scip_free(scip, nodesel) );

   /* free nodesel object */
   if( nodeseldata->deleteobject )
      delete nodeseldata->objnodesel;

   /* free nodesel data */
   delete nodeseldata;
   SCIPnodeselSetData(nodesel, NULL);
   
   return SCIP_OKAY;
}


/** initialization method of node selector (called after problem was transformed) */
static
DECL_NODESELINIT(nodeselInitObj)
{  /*lint --e{715}*/
   NODESELDATA* nodeseldata;

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);
   assert(nodeseldata->objnodesel != NULL);

   /* call virtual method of nodesel object */
   CHECK_OKAY( nodeseldata->objnodesel->scip_init(scip, nodesel) );

   return SCIP_OKAY;
}


/** deinitialization method of node selector (called before transformed problem is freed) */
static
DECL_NODESELEXIT(nodeselExitObj)
{  /*lint --e{715}*/
   NODESELDATA* nodeseldata;

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);
   assert(nodeseldata->objnodesel != NULL);

   /* call virtual method of nodesel object */
   CHECK_OKAY( nodeseldata->objnodesel->scip_exit(scip, nodesel) );

   return SCIP_OKAY;
}


/** node selection method of node selector */
static
DECL_NODESELSELECT(nodeselSelectObj)
{  /*lint --e{715}*/
   NODESELDATA* nodeseldata;

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);
   assert(nodeseldata->objnodesel != NULL);

   /* call virtual method of nodesel object */
   CHECK_OKAY( nodeseldata->objnodesel->scip_select(scip, nodesel, selnode) );

   return SCIP_OKAY;
}


/** node comparison method of node selector */
static
DECL_NODESELCOMP(nodeselCompObj)
{  /*lint --e{715}*/
   NODESELDATA* nodeseldata;

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);
   assert(nodeseldata->objnodesel != NULL);

   /* call virtual method of nodesel object */
   return nodeseldata->objnodesel->scip_comp(scip, nodesel, node1, node2);
}




/*
 * node selector specific interface methods
 */

/** creates the node selector for the given node selector object and includes it in SCIP */
RETCODE SCIPincludeObjNodesel(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjNodesel* objnodesel,        /**< node selector object */
   Bool             deleteobject        /**< should the node selector object be deleted when node selector is freed? */
   )
{
   NODESELDATA* nodeseldata;

   /* create node selector data */
   nodeseldata = new NODESELDATA;
   nodeseldata->objnodesel = objnodesel;
   nodeseldata->deleteobject = deleteobject;

   /* include node selector */
   CHECK_OKAY( SCIPincludeNodesel(scip, objnodesel->scip_name_, objnodesel->scip_desc_, 
         objnodesel->scip_stdpriority_, objnodesel->scip_memsavepriority_, objnodesel->scip_lowestboundfirst_,
         nodeselFreeObj, nodeselInitObj, nodeselExitObj, nodeselSelectObj, nodeselCompObj,
         nodeseldata) );

   return SCIP_OKAY;
}

/** returns the nodesel object of the given name, or NULL if not existing */
scip::ObjNodesel* SCIPfindObjNodesel(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of node selector */
   )
{
   NODESEL* nodesel;
   NODESELDATA* nodeseldata;

   nodesel = SCIPfindNodesel(scip, name);
   if( nodesel == NULL )
      return NULL;

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   return nodeseldata->objnodesel;
}
   
/** returns the nodesel object for the given node selector */
scip::ObjNodesel* SCIPgetObjNodesel(
   SCIP*            scip,               /**< SCIP data structure */
   NODESEL*         nodesel             /**< node selector */
   )
{
   NODESELDATA* nodeseldata;

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   return nodeseldata->objnodesel;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: objnodesel.cpp,v 1.4 2004/04/27 15:50:01 bzfpfend Exp $"

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

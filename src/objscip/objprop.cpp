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
#pragma ident "@(#) $Id: objprop.cpp,v 1.2 2004/12/14 12:08:01 bzfpfend Exp $"

/**@file   objprop.cpp
 * @brief  C++ wrapper for propagators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objprop.h"




/*
 * Data structures
 */

/** propagator data */
struct PropData
{
   scip::ObjProp*   objprop;            /**< propagator object */
   Bool             deleteobject;       /**< should the propagator object be deleted when propagator is freed? */
};




/*
 * Callback methods of propagator
 */

/** destructor of propagator to free user data (called when SCIP is exiting) */
static
DECL_PROPFREE(propFreeObj)
{  /*lint --e{715}*/
   PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->objprop != NULL);

   /* call virtual method of prop object */
   CHECK_OKAY( propdata->objprop->scip_free(scip, prop) );

   /* free prop object */
   if( propdata->deleteobject )
      delete propdata->objprop;

   /* free prop data */
   delete propdata;
   SCIPpropSetData(prop, NULL);
   
   return SCIP_OKAY;
}


/** initialization method of propagator (called after problem was transformed) */
static
DECL_PROPINIT(propInitObj)
{  /*lint --e{715}*/
   PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->objprop != NULL);

   /* call virtual method of prop object */
   CHECK_OKAY( propdata->objprop->scip_init(scip, prop) );

   return SCIP_OKAY;
}


/** deinitialization method of propagator (called before transformed problem is freed) */
static
DECL_PROPEXIT(propExitObj)
{  /*lint --e{715}*/
   PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->objprop != NULL);

   /* call virtual method of prop object */
   CHECK_OKAY( propdata->objprop->scip_exit(scip, prop) );

   return SCIP_OKAY;
}


/** execution method of propagator */
static
DECL_PROPEXEC(propExecObj)
{  /*lint --e{715}*/
   PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->objprop != NULL);

   /* call virtual method of prop object */
   CHECK_OKAY( propdata->objprop->scip_exec(scip, prop, result) );

   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator */
static
DECL_PROPRESPROP(propRespropObj)
{  /*lint --e{715}*/
   PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->objprop != NULL);

   /* call virtual method of prop object */
   CHECK_OKAY( propdata->objprop->scip_resprop(scip, prop, infervar, inferinfo, boundtype, bdchgidx, result) );

   return SCIP_OKAY;
}




/*
 * propagator specific interface methods
 */

/** creates the propagator for the given propagator object and includes it in SCIP */
RETCODE SCIPincludeObjProp(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjProp*   objprop,            /**< propagator object */
   Bool             deleteobject        /**< should the propagator object be deleted when propagator is freed? */
   )
{
   PROPDATA* propdata;

   /* create propagator data */
   propdata = new PROPDATA;
   propdata->objprop = objprop;
   propdata->deleteobject = deleteobject;

   /* include propagator */
   CHECK_OKAY( SCIPincludeProp(scip, objprop->scip_name_, objprop->scip_desc_, 
         objprop->scip_priority_, objprop->scip_freq_,
         propFreeObj, propInitObj, propExitObj, propExecObj, propRespropObj,
         propdata) );

   return SCIP_OKAY;
}

/** returns the prop object of the given name, or NULL if not existing */
scip::ObjProp* SCIPfindObjProp(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of propagator */
   )
{
   PROP* prop;
   PROPDATA* propdata;

   prop = SCIPfindProp(scip, name);
   if( prop == NULL )
      return NULL;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   return propdata->objprop;
}
   
/** returns the prop object for the given propagator */
scip::ObjProp* SCIPgetObjProp(
   SCIP*            scip,               /**< SCIP data structure */
   PROP*            prop                /**< propagator */
   )
{
   PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   return propdata->objprop;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   objiisfinder.cpp
 * @brief  C++ wrapper for IIS finders
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objiisfinder.h"

/*
 * Data structures
 */

/** iis finder data */
struct SCIP_IISfinderData
{
   scip::ObjIISfinder*   objiisfinder;       /**< iis finder object */
   SCIP_Bool             deleteobject;       /**< should the iis finder object be deleted when iis finder is freed? */
};

/*
 * Callback methods of iis finder
 */

extern "C"
{

/** copy method for iis finder plugins (called when SCIP copies plugins) */
static
SCIP_DECL_IISFINDERCOPY(iisfinderCopyObj)
{  /*lint --e{715}*/
   SCIP_IISFINDERDATA* iisfinderdata;

   assert(scip != NULL);

   iisfinderdata = SCIPiisfinderGetData(iisfinder);
   assert(iisfinderdata != NULL);
   assert(iisfinderdata->objiisfinder != NULL);
   assert(iisfinderdata->objiisfinder->scip_ != scip);

   if( iisfinderdata->objiisfinder->iscloneable() )
   {
      scip::ObjIISfinder* newobjiisfinder;
      newobjiisfinder = dynamic_cast<scip::ObjIISfinder*> (iisfinderdata->objiisfinder->clone(scip));

      /* call include method of iis finder object */
      SCIP_CALL( SCIPincludeObjIISfinder(scip, newobjiisfinder, TRUE) );
   }

   return SCIP_OKAY;
}

/** destructor of iis finder to free user data (called when SCIP is exiting) */
static
SCIP_DECL_IISFINDERFREE(iisfinderFreeObj)
{  /*lint --e{715}*/
   SCIP_IISFINDERDATA* iisfinderdata;

   iisfinderdata = SCIPiisfinderGetData(iisfinder);
   assert(iisfinderdata != NULL);
   assert(iisfinderdata->objiisfinder != NULL);
   assert(iisfinderdata->objiisfinder->scip_ == scip);

   /* call virtual method of iisfinder object */
   SCIP_CALL( iisfinderdata->objiisfinder->scip_free(scip, iisfinder) );

   /* free iisfinder object */
   if( iisfinderdata->deleteobject )
      delete iisfinderdata->objiisfinder;

   /* free iisfinder data */
   delete iisfinderdata;
   SCIPiisfinderSetData(iisfinder, NULL); /*lint !e64*/

   return SCIP_OKAY;
}

/** iis finder execution method of iisfinder */
static
SCIP_DECL_IISFINDEREXEC(iisfinderExecObj)
{  /*lint --e{715}*/
   SCIP_IISFINDERDATA* iisfinderdata;

   iisfinderdata = SCIPiisfinderGetData(iisfinder);
   assert(iisfinderdata != NULL);
   assert(iisfinderdata->objiisfinder != NULL);

   /* call virtual method of iisfinder object */
   SCIP_CALL( iisfinderdata->objiisfinder->scip_exec(iis, iisfinder, result) );

   return SCIP_OKAY;
}
}

/*
 * iis finder specific interface methods
 */

/** creates the iis finder for the given iis finder object and includes it in SCIP */
SCIP_RETCODE SCIPincludeObjIISfinder(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjIISfinder*   objiisfinder,       /**< iis finder object */
   SCIP_Bool             deleteobject        /**< should the iis finder object be deleted when iis finder is freed? */
   )
{
   SCIP_IISFINDERDATA* iisfinderdata;

   assert(scip != NULL);
   assert(objiisfinder != NULL);

   /* create iis finder data */
   iisfinderdata = new SCIP_IISFINDERDATA;
   iisfinderdata->objiisfinder = objiisfinder;
   iisfinderdata->deleteobject = deleteobject;

   /* include iis finder */
   SCIP_CALL( SCIPincludeIISfinder(scip, objiisfinder->scip_name_, objiisfinder->scip_desc_,
         objiisfinder->scip_priority_,
         iisfinderCopyObj,
         iisfinderFreeObj,
         iisfinderExecObj,
         iisfinderdata) ); /*lint !e429*/

   return SCIP_OKAY; /*lint !e429*/
}

/** returns the iis finder object of the given name, or 0 if not existing */
scip::ObjIISfinder* SCIPfindObjIISfinder(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of iis finder */
   )
{
   SCIP_IISFINDER * iisfinder;
   SCIP_IISFINDERDATA* iisfinderdata;

   iisfinder = SCIPfindIISfinder(scip, name);
   if( iisfinder == NULL )
      return 0;

   iisfinderdata = SCIPiisfinderGetData(iisfinder);
   assert(iisfinderdata != NULL);

   return iisfinderdata->objiisfinder;
}

/** returns the iis finder object for the given iis finder */
scip::ObjIISfinder* SCIPgetObjIISfinder(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_IISFINDER*       iisfinder           /**< iis finder */
   )
{
   SCIP_IISFINDERDATA* iisfinderdata;

   assert(scip != NULL);
   iisfinderdata = SCIPiisfinderGetData(iisfinder);
   assert(iisfinderdata != NULL);

   return iisfinderdata->objiisfinder;
}

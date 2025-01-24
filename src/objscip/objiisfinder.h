/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

/**@file   objiisfinder.h
 * @brief  C++ wrapper for iis finders
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJIISFINDER_H__
#define __SCIP_OBJIISFINDER_H__

#include <cstring>
#include <utility>

#include "scip/scip.h"
#include "objscip/objcloneable.h"

namespace scip
{

/** @brief C++ wrapper for iis finders
 *
 *  This class defines the interface for iis finders implemented in C++.
 *
 *  - \ref IISFINDER "Instructions for implementing an iis finder"
 *  - \ref IISFINDERS "List of available iis finders"
 *  - \ref type_iisfinder.h "Corresponding C interface"
 */
class ObjIISfinder : public ObjCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the iis finder */
   char* scip_name_;

   /** description of the iis finder */
   char* scip_desc_;

   /** priority of the iis finder */
   const int scip_priority_;

   /** default constructor */
   ObjIISfinder(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of iis finder */
      const char*        desc,               /**< description of iis finder */
      int                priority            /**< priority of the iis finder */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_priority_(priority)
   {
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** copy constructor */
   ObjIISfinder(const ObjIISfinder& o) : ObjIISfinder(o.scip_, o.scip_name_, o.scip_desc_, o.scip_priority_) {}

   /** move constructor */
   ObjIISfinder(ObjIISfinder&& o) : scip_(o.scip_), scip_name_(0), scip_desc_(0), scip_priority_(o.scip_priority_)
   {
      std::swap(scip_name_, o.scip_name_);
      std::swap(scip_desc_, o.scip_desc_);
   }

   /** destructor */
   virtual ~ObjIISfinder()
   {
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjIISfinder& operator=(const ObjIISfinder& o) = delete;

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjIISfinder& operator=(ObjIISfinder&& o) = delete;

   /** destructor of iis finder to free user data (called when SCIP is exiting)
    *
    *  @see SCIP_DECL_IISFINDERFREE(x) in @ref type_iisfinder.h
    */
   virtual SCIP_DECL_IISFINDERFREE(scip_free)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** iis finder execution method of iis finder
    *
    *  @see SCIP_DECL_IISFINDEREXEC(x) in @ref type_iisfinder.h
    */
   virtual SCIP_DECL_IISFINDEREXEC(scip_exec) = 0;
};

} /* namespace scip */



/** creates the iis finder for the given iis finder object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is responsible for deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyIISfinder* myiisfinder = new MyIISfinder(...);
 *       SCIP_CALL( SCIPincludeObjIISfinder(scip, &myiisfinder, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete myiisfinder;    // delete iisfinder AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjIISfinder(scip, new MyIISfinder(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyIISfinder is called here
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeObjIISfinder(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjIISfinder*   objiisfinder,       /**< iis finder object */
   SCIP_Bool             deleteobject        /**< should the iis finder object be deleted when iis finder is freed? */
   );

/** returns the iis finder object of the given name, or 0 if not existing */
SCIP_EXPORT
scip::ObjIISfinder* SCIPfindObjIISfinder(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of iis finder */
   );

/** returns the iis finder object for the given iis finder */
SCIP_EXPORT
scip::ObjIISfinder* SCIPgetObjIISfinder(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_IISFINDER*       iisfinder           /**< iis finder */
   );

#endif

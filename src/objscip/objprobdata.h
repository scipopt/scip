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
#pragma ident "@(#) $Id: objprobdata.h,v 1.4 2004/04/15 10:41:25 bzfpfend Exp $"

/**@file   objprobdata.h
 * @brief  C++ wrapper for user problem data
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __OBJPROBDATA_H__
#define __OBJPROBDATA_H__


extern "C" 
{
#include "scip.h"
}


namespace scip
{

/** C++ wrapper object for user problem data */
class ObjProbData
{
public:
   /** default constructor */
   ObjProbData()
   {
   }

   /** destructor */
   virtual ~ObjProbData()
   {
   }

   /** destructor of user problem data to free original user data (called when original problem is freed)
    *
    *  If the "deleteobject" flag in the SCIPcreateObjProb() method was set to TRUE, this method is not needed,
    *  because all the work to delete the user problem data can be done in the destructor of the user problem
    *  data object. If the "deleteobject" flag was set to FALSE, and the user problem data object stays alive
    *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
    *  longer needed.
    */
   virtual RETCODE scip_delorig(
      SCIP*         scip                /**< SCIP data structure */
      )
   {
      return SCIP_OKAY;
   }

   /** creates user data of transformed problem by transforming the original user problem data
    *  (called when problem solving starts)
    *
    *  The user has two possibilities to implement this method:
    *   1. Return the pointer to the original problem data object as pointer to the transformed problem data object.
    *      The user may modify some internal attributes, but he has to make sure, that these modifications are
    *      reversed in the scip_deltrans() method, such that the original problem data is restored. In this case,
    *      he should set *deleteobject to FALSE, because the problem data must not be destructed by SCIP after the
    *      solving process is terminated.
    *   2. Call the copy constructor of the problem data object and return the created copy as transformed problem
    *      data object. In this case, he probably wants to set *deleteobject to TRUE, thus letting SCIP call the
    *      destructor of the object if the transformed problem data is no longer needed.
    */
   virtual RETCODE scip_trans(
      SCIP*         scip,               /**< SCIP data structure */
      ObjProbData** objprobdata,        /**< pointer to store the transformed problem data object */
      Bool*         deleteobject        /**< pointer to store whether SCIP should delete the object after solving */
      )
   {
      assert(objprobdata != NULL);
      assert(deleteobject != NULL);

      /* the default implementation just copies the pointer to the problem data object;
       * SCIP must not destruct the transformed problem data object, because the original problem data must stay alive
       */
      *objprobdata = this;
      *deleteobject = FALSE;

      return SCIP_OKAY;
   }      

   /** destructor of user problem data to free transformed user data (called when transformed problem is freed)
    *
    *  If the "*deleteobject" flag in the scip_trans() method was set to TRUE, this method is not needed,
    *  because all the work to delete the user problem data can be done in the destructor of the user problem
    *  data object. If the "*deleteobject" flag was set to FALSE, and the user problem data object stays alive
    *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
    *  longer needed.
    */
   virtual RETCODE scip_deltrans(
      SCIP*         scip                /**< SCIP data structure */
      )
   {
      return SCIP_OKAY;
   }

};

} /* namespace scip */


   
/** creates empty problem, initializes all solving data structures, and sets the user problem data to point to the
 *  given user data object
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       MyProbData* myprobdata = new MyProbData(...);
 *       CHECK_OKAY( SCIPcreateObjProb(scip, "probname", &myprobdata, FALSE) );
 *       ... // solve the problem
 *       CHECK_OKAY( SCIPfreeProb(scip) );
 *       delete myprobdata;    // delete probdata AFTER SCIPfreeProb() !
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfreeProb() call:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       CHECK_OKAY( SCIPcreateObjProb(scip, "probname", new MyProbData(...), TRUE) );
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );  // problem is freed and destructor of MyProbData is called here
 */
extern
RETCODE SCIPcreateObjProb(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< problem name */
   scip::ObjProbData* objprobdata,      /**< user problem data object */
   Bool             deleteobject        /**< should the user problem data object be deleted when problem is freed? */
   );

/** gets user problem data object
 *  Warning! This method should only be called after a problem was created with SCIPcreateObjProb().
 *  Otherwise, a segmentation fault may arise, or an undefined pointer is returned.
 */
extern
scip::ObjProbData* SCIPgetObjProbData(
   SCIP*            scip                /**< SCIP data structure */
   );

#endif

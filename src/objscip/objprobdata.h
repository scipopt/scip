/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: objprobdata.h,v 1.1 2003/12/08 11:51:03 bzfpfend Exp $"

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

   /** destructor of user problem data to free user data (called when problem is freed)
    *  If the "deleteobject" flag in the SCIPcreateObjProb() method was set to TRUE, this method is not needed,
    *  because all the work to delete the user problem data can be done in the destructor of the user problem
    *  data object. If the "deleteobject" flag was set to FALSE, and the user problem data object stays alive
    *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
    *  longer needed.
    */
   virtual RETCODE scip_delete(
      SCIP*         scip                /**< SCIP data structure */
      )
   {
      return SCIP_OKAY;
   }
   
};

} /* namespace scip */


   
/** creates empty problem, initializes all solving data structures, and sets the user problem data to point to the
 *  given user data object
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

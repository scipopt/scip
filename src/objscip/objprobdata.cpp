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
#pragma ident "@(#) $Id: objprobdata.cpp,v 1.1 2003/12/08 11:51:03 bzfpfend Exp $"

/**@file   objprobdata.cpp
 * @brief  C++ wrapper for user problem data
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objprobdata.h"




/*
 * Data structures
 */

/** user problem data */
struct ProbData
{
   scip::ObjProbData* objprobdata;      /**< user problem data object */
   Bool             deleteobject;       /**< should the user problem data object be deleted when problem is freed? */
};




/*
 * Callback methods of user problem data
 */

/** frees user data of original problem (called when the original problem is freed) */
static
DECL_PROBDELORIG(probDelorigObj)
{  /*lint --e{715}*/
   assert(probdata != NULL);
   assert(*probdata != NULL);
   assert((*probdata)->objprobdata != NULL);

   /* call virtual method of probdata object */
   CHECK_OKAY( (*probdata)->objprobdata->scip_delete(scip) );

   /* free probdata object */
   if( (*probdata)->deleteobject )
      delete (*probdata)->objprobdata;

   /* free probdata data */
   delete *probdata;
   *probdata = NULL;
   
   return SCIP_OKAY;
}


/** creates user data of transformed problem by transforming the original user problem data
 *  (called when problem solving starts)
 */
#define probTransObj NULL


/** frees user data of transformed problem (called when the transformed problem is freed) */
#define probDeltransObj NULL




/*
 * user problem data specific interface methods
 */

/** creates empty problem, initializes all solving data structures, and sets the user problem data to point to the
 *  given user data object
 */
RETCODE SCIPcreateObjProb(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< problem name */
   scip::ObjProbData* objprobdata,      /**< user problem data object */
   Bool             deleteobject        /**< should the user problem data object be deleted when problem is freed? */
   )
{
   PROBDATA* probdata;

   /* create user problem data */
   probdata = new PROBDATA;
   probdata->objprobdata = objprobdata;
   probdata->deleteobject = deleteobject;

   /* create problem */
   CHECK_OKAY( SCIPcreateProb(scip, name, probDelorigObj, probTransObj, probDeltransObj, probdata) );

   return SCIP_OKAY;
}

/** gets user problem data object
 *  Warning! This method should only be called after a problem was created with SCIPcreateObjProb().
 *  Otherwise, a segmentation fault may arise, or an undefined pointer is returned.
 */
scip::ObjProbData* SCIPgetObjProbData(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   PROBDATA* probdata;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->objprobdata;
}


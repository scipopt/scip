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
#pragma ident "@(#) $Id: objvardata.cpp,v 1.3 2005/01/18 09:26:50 bzfpfend Exp $"

/**@file   objvardata.cpp
 * @brief  C++ wrapper for user variable data
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objvardata.h"




/*
 * Data structures
 */

/** user variable data */
struct VarData
{
   scip::ObjVardata* objvardata;        /**< user variable data object */
   Bool             deleteobject;       /**< should the user variable data object be deleted when variable is freed? */
};




/*
 * Callback methods of user variable data
 */

/** frees user data of original variable (called when the original variable is freed) */
static
DECL_VARDELORIG(varDelorigObj)
{  /*lint --e{715}*/
   assert(vardata != NULL);
   assert(*vardata != NULL);
   assert((*vardata)->objvardata != NULL);

   /* call virtual method of vardata object */
   CHECK_OKAY( (*vardata)->objvardata->scip_delorig(scip, var) );

   /* free vardata object */
   if( (*vardata)->deleteobject )
      delete (*vardata)->objvardata;

   /* free vardata data */
   delete *vardata;
   *vardata = NULL;
   
   return SCIP_OKAY;
}


/** creates user data of transformed variable by transforming the original user variable data
 *  (called after variable was transformed)
 */
static
DECL_VARTRANS(varTransObj)
{  /*lint --e{715}*/
   scip::ObjVardata* objvardata;
   Bool deleteobject;

   assert(sourcedata != NULL);
   assert(sourcedata->objvardata != NULL);
   assert(targetdata != NULL);
   assert(*targetdata == NULL);

   /* call virtual method of vardata object */
   CHECK_OKAY( sourcedata->objvardata->scip_trans(scip, targetvar, &objvardata, &deleteobject) );

   /* create transformed user variable data */
   *targetdata = new VARDATA;
   (*targetdata)->objvardata = objvardata;
   (*targetdata)->deleteobject = deleteobject;

   return SCIP_OKAY;
}


/** frees user data of transformed variable (called when the transformed variable is freed) */
static
DECL_VARDELTRANS(varDeltransObj)
{  /*lint --e{715}*/
   assert(vardata != NULL);
   assert(*vardata != NULL);
   assert((*vardata)->objvardata != NULL);

   /* call virtual method of vardata object */
   CHECK_OKAY( (*vardata)->objvardata->scip_deltrans(scip, var) );

   /* free vardata object */
   if( (*vardata)->deleteobject )
      delete (*vardata)->objvardata;

   /* free vardata data */
   delete *vardata;
   *vardata = NULL;
   
   return SCIP_OKAY;
}





/*
 * user variable data specific interface methods
 */

/** create and capture problem variable and associates the given variable data with the variable;
 *  if variable is of integral type, fractional bounds are automatically rounded
 */
RETCODE SCIPcreateObjVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            var,                /**< pointer to variable object */
   const char*      name,               /**< name of variable, or NULL for automatic name creation */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype,            /**< type of variable */
   Bool             initial,            /**< should var's column be present in the initial root LP? */
   Bool             removeable,         /**< is var's column removeable from the LP (due to aging or cleanup)? */
   scip::ObjVardata* objvardata,        /**< user variable data object */
   Bool             deleteobject        /**< should the user variable data object be deleted when variable is freed? */
   )
{
   VARDATA* vardata;

   /* create user variable data */
   vardata = new VARDATA;
   vardata->objvardata = objvardata;
   vardata->deleteobject = deleteobject;

   /* create variable */
   CHECK_OKAY( SCIPcreateVar(scip, var, name, lb, ub, obj, vartype, initial, removeable, 
         varDelorigObj, varTransObj, varDeltransObj, vardata) );

   return SCIP_OKAY;
}

/** gets user variable data object for given problem variable
 *  Warning! This method should only be called after a variable was created with SCIPcreateObjVar().
 *  Otherwise, a segmentation fault may arise, or an undefined pointer is returned.
 */
scip::ObjVardata* SCIPgetObjVardata(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< problem variable */
   )
{
   VARDATA* vardata;

   vardata = SCIPgetVarData(scip, var);
   assert(vardata != NULL);

   return vardata->objvardata;
}


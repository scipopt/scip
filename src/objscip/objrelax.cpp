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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: objrelax.cpp,v 1.5 2005/02/07 14:08:25 bzfpfend Exp $"

/**@file   objrelax.cpp
 * @brief  C++ wrapper for relaxators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objrelax.h"




/*
 * Data structures
 */

/** relaxator data */
struct RelaxData
{
   scip::ObjRelax*  objrelax;           /**< relaxator object */
   Bool             deleteobject;       /**< should the relaxator object be deleted when relaxator is freed? */
};




/*
 * Callback methods of relaxator
 */

/** destructor of relaxator to free user data (called when SCIP is exiting) */
static
DECL_RELAXFREE(relaxFreeObj)
{  /*lint --e{715}*/
   RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(relaxdata->objrelax != NULL);

   /* call virtual method of relax object */
   CHECK_OKAY( relaxdata->objrelax->scip_free(scip, relax) );

   /* free relax object */
   if( relaxdata->deleteobject )
      delete relaxdata->objrelax;

   /* free relax data */
   delete relaxdata;
   SCIPrelaxSetData(relax, NULL);
   
   return SCIP_OKAY;
}


/** initialization method of relaxator (called after problem was transformed) */
static
DECL_RELAXINIT(relaxInitObj)
{  /*lint --e{715}*/
   RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(relaxdata->objrelax != NULL);

   /* call virtual method of relax object */
   CHECK_OKAY( relaxdata->objrelax->scip_init(scip, relax) );

   return SCIP_OKAY;
}


/** deinitialization method of relaxator (called before transformed problem is freed) */
static
DECL_RELAXEXIT(relaxExitObj)
{  /*lint --e{715}*/
   RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(relaxdata->objrelax != NULL);

   /* call virtual method of relax object */
   CHECK_OKAY( relaxdata->objrelax->scip_exit(scip, relax) );

   return SCIP_OKAY;
}


/** solving process initialization method of relaxator (called when branch and bound process is about to begin) */
static
DECL_RELAXINITSOL(relaxInitsolObj)
{  /*lint --e{715}*/
   RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(relaxdata->objrelax != NULL);

   /* call virtual method of relax object */
   CHECK_OKAY( relaxdata->objrelax->scip_initsol(scip, relax) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of relaxator (called before branch and bound process data is freed) */
static
DECL_RELAXEXITSOL(relaxExitsolObj)
{  /*lint --e{715}*/
   RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(relaxdata->objrelax != NULL);

   /* call virtual method of relax object */
   CHECK_OKAY( relaxdata->objrelax->scip_exitsol(scip, relax) );

   return SCIP_OKAY;
}


/** execution method of relaxator */
static
DECL_RELAXEXEC(relaxExecObj)
{  /*lint --e{715}*/
   RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(relaxdata->objrelax != NULL);

   /* call virtual method of relax object */
   CHECK_OKAY( relaxdata->objrelax->scip_exec(scip, relax, result) );

   return SCIP_OKAY;
}




/*
 * relaxator specific interface methods
 */

/** creates the relaxator for the given relaxator object and includes it in SCIP */
RETCODE SCIPincludeObjRelax(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjRelax*  objrelax,           /**< relaxator object */
   Bool             deleteobject        /**< should the relaxator object be deleted when relaxator is freed? */
   )
{
   RELAXDATA* relaxdata;

   /* create relaxator data */
   relaxdata = new RELAXDATA;
   relaxdata->objrelax = objrelax;
   relaxdata->deleteobject = deleteobject;

   /* include relaxator */
   CHECK_OKAY( SCIPincludeRelax(scip, objrelax->scip_name_, objrelax->scip_desc_, 
         objrelax->scip_priority_, objrelax->scip_freq_,
         relaxFreeObj, relaxInitObj, relaxExitObj, 
         relaxInitsolObj, relaxExitsolObj, relaxExecObj,
         relaxdata) );

   return SCIP_OKAY;
}

/** returns the relax object of the given name, or NULL if not existing */
scip::ObjRelax* SCIPfindObjRelax(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of relaxator */
   )
{
   RELAX* relax;
   RELAXDATA* relaxdata;

   relax = SCIPfindRelax(scip, name);
   if( relax == NULL )
      return NULL;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->objrelax;
}
   
/** returns the relax object for the given relaxator */
scip::ObjRelax* SCIPgetObjRelax(
   SCIP*            scip,               /**< SCIP data structure */
   RELAX*           relax               /**< relaxator */
   )
{
   RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->objrelax;
}

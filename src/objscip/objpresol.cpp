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
#pragma ident "@(#) $Id: objpresol.cpp,v 1.8 2005/01/18 09:26:49 bzfpfend Exp $"

/**@file   objpresol.cpp
 * @brief  C++ wrapper for presolvers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objpresol.h"




/*
 * Data structures
 */

/** presolver data */
struct PresolData
{
   scip::ObjPresol* objpresol;          /**< presolver object */
   Bool             deleteobject;       /**< should the presolver object be deleted when presolver is freed? */
};




/*
 * Callback methods of presolver
 */

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
DECL_PRESOLFREE(presolFreeObj)
{  /*lint --e{715}*/
   PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   assert(presoldata->objpresol != NULL);

   /* call virtual method of presol object */
   CHECK_OKAY( presoldata->objpresol->scip_free(scip, presol) );

   /* free presol object */
   if( presoldata->deleteobject )
      delete presoldata->objpresol;

   /* free presol data */
   delete presoldata;
   SCIPpresolSetData(presol, NULL);
   
   return SCIP_OKAY;
}


/** initialization method of presolver (called after problem was transformed) */
static
DECL_PRESOLINIT(presolInitObj)
{  /*lint --e{715}*/
   PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   assert(presoldata->objpresol != NULL);

   /* call virtual method of presol object */
   CHECK_OKAY( presoldata->objpresol->scip_init(scip, presol) );

   return SCIP_OKAY;
}


/** deinitialization method of presolver (called before transformed problem is freed) */
static
DECL_PRESOLEXIT(presolExitObj)
{  /*lint --e{715}*/
   PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   assert(presoldata->objpresol != NULL);

   /* call virtual method of presol object */
   CHECK_OKAY( presoldata->objpresol->scip_exit(scip, presol) );

   return SCIP_OKAY;
}


/** presolving initialization method of presolver (called when presolving is about to begin) */
static
DECL_PRESOLINITPRE(presolInitpreObj)
{  /*lint --e{715}*/
   PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   assert(presoldata->objpresol != NULL);

   /* call virtual method of presol object */
   CHECK_OKAY( presoldata->objpresol->scip_initpre(scip, presol, result) );

   return SCIP_OKAY;
}


/** presolving deinitialization method of presolver (called after presolving has been finished) */
static
DECL_PRESOLEXITPRE(presolExitpreObj)
{  /*lint --e{715}*/
   PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   assert(presoldata->objpresol != NULL);

   /* call virtual method of presol object */
   CHECK_OKAY( presoldata->objpresol->scip_exitpre(scip, presol, result) );

   return SCIP_OKAY;
}


/** execution method of presolver */
static
DECL_PRESOLEXEC(presolExecObj)
{  /*lint --e{715}*/
   PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   assert(presoldata->objpresol != NULL);

   /* call virtual method of presol object */
   CHECK_OKAY( presoldata->objpresol->scip_exec(scip, presol, nrounds,
         nnewfixedvars, nnewaggrvars, nnewchgvartypes, nnewchgbds, nnewholes,
         nnewdelconss, nnewupgdconss, nnewchgcoefs, nnewchgsides,
         nfixedvars, naggrvars, nchgvartypes, nchgbds, naddholes,
         ndelconss, nupgdconss, nchgcoefs, nchgsides, result) );

   return SCIP_OKAY;
}




/*
 * presolver specific interface methods
 */

/** creates the presolver for the given presolver object and includes it in SCIP */
RETCODE SCIPincludeObjPresol(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjPresol* objpresol,          /**< presolver object */
   Bool             deleteobject        /**< should the presolver object be deleted when presolver is freed? */
   )
{
   PRESOLDATA* presoldata;

   /* create presolver data */
   presoldata = new PRESOLDATA;
   presoldata->objpresol = objpresol;
   presoldata->deleteobject = deleteobject;

   /* include presolver */
   CHECK_OKAY( SCIPincludePresol(scip, objpresol->scip_name_, objpresol->scip_desc_, 
         objpresol->scip_priority_, objpresol->scip_maxrounds_,
         presolFreeObj, presolInitObj, presolExitObj, 
         presolInitpreObj, presolExitpreObj, presolExecObj,
         presoldata) );

   return SCIP_OKAY;
}

/** returns the presol object of the given name, or NULL if not existing */
scip::ObjPresol* SCIPfindObjPresol(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of presolver */
   )
{
   PRESOL* presol;
   PRESOLDATA* presoldata;

   presol = SCIPfindPresol(scip, name);
   if( presol == NULL )
      return NULL;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   return presoldata->objpresol;
}
   
/** returns the presol object for the given presolver */
scip::ObjPresol* SCIPgetObjPresol(
   SCIP*            scip,               /**< SCIP data structure */
   PRESOL*          presol              /**< presolver */
   )
{
   PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   return presoldata->objpresol;
}

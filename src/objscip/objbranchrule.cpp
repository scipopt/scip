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
#pragma ident "@(#) $Id: objbranchrule.cpp,v 1.7 2004/09/21 12:14:15 bzfpfend Exp $"

/**@file   objbranchrule.cpp
 * @brief  C++ wrapper for branching rules
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objbranchrule.h"




/*
 * Data structures
 */

/** branching rule data */
struct BranchruleData
{
   scip::ObjBranchrule* objbranchrule;  /**< branching rule object */
   Bool             deleteobject;       /**< should the branching rule object be deleted when branching rule is freed? */
};




/*
 * Callback methods of branching rule
 */

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
DECL_BRANCHFREE(branchFreeObj)
{  /*lint --e{715}*/
   BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->objbranchrule != NULL);

   /* call virtual method of branchrule object */
   CHECK_OKAY( branchruledata->objbranchrule->scip_free(scip, branchrule) );

   /* free branchrule object */
   if( branchruledata->deleteobject )
      delete branchruledata->objbranchrule;

   /* free branchrule data */
   delete branchruledata;
   SCIPbranchruleSetData(branchrule, NULL);
   
   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
static
DECL_BRANCHINIT(branchInitObj)
{  /*lint --e{715}*/
   BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->objbranchrule != NULL);

   /* call virtual method of branchrule object */
   CHECK_OKAY( branchruledata->objbranchrule->scip_init(scip, branchrule) );

   return SCIP_OKAY;
}


/** deinitialization method of branching rule (called before transformed problem is freed) */
static
DECL_BRANCHEXIT(branchExitObj)
{  /*lint --e{715}*/
   BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->objbranchrule != NULL);

   /* call virtual method of branchrule object */
   CHECK_OKAY( branchruledata->objbranchrule->scip_exit(scip, branchrule) );

   return SCIP_OKAY;
}


/** branching execution method for fractional LP solutions */
static
DECL_BRANCHEXECLP(branchExeclpObj)
{  /*lint --e{715}*/
   BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->objbranchrule != NULL);

   /* call virtual method of branchrule object */
   CHECK_OKAY( branchruledata->objbranchrule->scip_execlp(scip, branchrule, allowaddcons, result) );

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static
DECL_BRANCHEXECPS(branchExecpsObj)
{  /*lint --e{715}*/
   BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->objbranchrule != NULL);

   /* call virtual method of branchrule object */
   CHECK_OKAY( branchruledata->objbranchrule->scip_execps(scip, branchrule, allowaddcons, result) );

   return SCIP_OKAY;
}




/*
 * branching rule specific interface methods
 */

/** creates the branching rule for the given branching rule object and includes it in SCIP */
RETCODE SCIPincludeObjBranchrule(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjBranchrule* objbranchrule,  /**< branching rule object */
   Bool             deleteobject        /**< should the branching rule object be deleted when branching rule is freed? */
   )
{
   BRANCHRULEDATA* branchruledata;

   /* create branching rule data */
   branchruledata = new BRANCHRULEDATA;
   branchruledata->objbranchrule = objbranchrule;
   branchruledata->deleteobject = deleteobject;

   /* include branching rule */
   CHECK_OKAY( SCIPincludeBranchrule(scip, objbranchrule->scip_name_, objbranchrule->scip_desc_, 
                  objbranchrule->scip_priority_, objbranchrule->scip_maxdepth_,
                  branchFreeObj, branchInitObj, branchExitObj, branchExeclpObj, branchExecpsObj,
                  branchruledata) );

   return SCIP_OKAY;
}

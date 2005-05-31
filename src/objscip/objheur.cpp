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
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: objheur.cpp,v 1.14 2005/05/31 17:20:08 bzfpfend Exp $"

/**@file   objheur.cpp
 * @brief  C++ wrapper for primal heuristics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objheur.h"




/*
 * Data structures
 */

/** primal heuristic data */
struct HeurData
{
   scip::ObjHeur*   objheur;            /**< primal heuristic object */
   Bool             deleteobject;       /**< should the primal heuristic object be deleted when heuristic is freed? */
};




/*
 * Callback methods of primal heuristic
 */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
DECL_HEURFREE(heurFreeObj)
{  /*lint --e{715}*/
   HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->objheur != NULL);

   /* call virtual method of heur object */
   CHECK_OKAY( heurdata->objheur->scip_free(scip, heur) );

   /* free heur object */
   if( heurdata->deleteobject )
      delete heurdata->objheur;

   /* free heur data */
   delete heurdata;
   SCIPheurSetData(heur, NULL);
   
   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
DECL_HEURINIT(heurInitObj)
{  /*lint --e{715}*/
   HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->objheur != NULL);

   /* call virtual method of heur object */
   CHECK_OKAY( heurdata->objheur->scip_init(scip, heur) );

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
DECL_HEUREXIT(heurExitObj)
{  /*lint --e{715}*/
   HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->objheur != NULL);

   /* call virtual method of heur object */
   CHECK_OKAY( heurdata->objheur->scip_exit(scip, heur) );

   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
DECL_HEURINITSOL(heurInitsolObj)
{  /*lint --e{715}*/
   HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->objheur != NULL);

   /* call virtual method of heur object */
   CHECK_OKAY( heurdata->objheur->scip_initsol(scip, heur) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
DECL_HEUREXITSOL(heurExitsolObj)
{  /*lint --e{715}*/
   HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->objheur != NULL);

   /* call virtual method of heur object */
   CHECK_OKAY( heurdata->objheur->scip_exitsol(scip, heur) );

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
DECL_HEUREXEC(heurExecObj)
{  /*lint --e{715}*/
   HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->objheur != NULL);

   /* call virtual method of heur object */
   CHECK_OKAY( heurdata->objheur->scip_exec(scip, heur, result) );

   return SCIP_OKAY;
}




/*
 * primal heuristic specific interface methods
 */

/** creates the primal heuristic for the given primal heuristic object and includes it in SCIP */
RETCODE SCIPincludeObjHeur(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjHeur*   objheur,            /**< primal heuristic object */
   Bool             deleteobject        /**< should the primal heuristic object be deleted when heuristic is freed? */
   )
{
   HEURDATA* heurdata;

   /* create primal heuristic data */
   heurdata = new HEURDATA;
   heurdata->objheur = objheur;
   heurdata->deleteobject = deleteobject;

   /* include primal heuristic */
   CHECK_OKAY( SCIPincludeHeur(scip, objheur->scip_name_, objheur->scip_desc_, objheur->scip_dispchar_,
         objheur->scip_priority_, objheur->scip_freq_, objheur->scip_freqofs_, objheur->scip_maxdepth_,
         objheur->scip_pseudonodes_, objheur->scip_duringplunging_, objheur->scip_afternode_,
         heurFreeObj, heurInitObj, heurExitObj, 
         heurInitsolObj, heurExitsolObj, heurExecObj,
         heurdata) );

   return SCIP_OKAY;
}

/** returns the heur object of the given name, or NULL if not existing */
scip::ObjHeur* SCIPfindObjHeur(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of primal heuristic */
   )
{
   HEUR* heur;
   HEURDATA* heurdata;

   heur = SCIPfindHeur(scip, name);
   if( heur == NULL )
      return NULL;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   return heurdata->objheur;
}
   
/** returns the heur object for the given primal heuristic */
scip::ObjHeur* SCIPgetObjHeur(
   SCIP*            scip,               /**< SCIP data structure */
   HEUR*            heur                /**< primal heuristic */
   )
{
   HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   return heurdata->objheur;
}

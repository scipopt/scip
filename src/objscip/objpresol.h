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
#pragma ident "@(#) $Id: objpresol.h,v 1.3 2003/12/08 13:24:53 bzfpfend Exp $"

/**@file   objpresol.h
 * @brief  C++ wrapper for presolvers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __OBJPRESOL_H__
#define __OBJPRESOL_H__


extern "C" 
{
#include "scip.h"
}


namespace scip
{

/** C++ wrapper object for presolvers */
class ObjPresol
{
public:
   /** name of the presolver */
   const char* const scip_name_;
   
   /** description of the presolver */
   const char* const scip_desc_;
   
   /** default priority of the presolver */
   const int scip_priority_;

   /** default constructor */
   ObjPresol(
      const char*   name,               /**< name of presolver */
      const char*   desc,               /**< description of presolver */
      int           priority            /**< priority of the presolver */
      )
      : scip_name_(name),
        scip_desc_(desc),
        scip_priority_(priority)
   {
   }

   /** destructor of presolver to free user data (called when SCIP is exiting) */
   virtual RETCODE scip_free(
      SCIP*         scip,               /**< SCIP data structure */
      PRESOL*       presol              /**< the presolver itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** initialization method of presolver (called when problem solving starts) */
   virtual RETCODE scip_init(
      SCIP*         scip,               /**< SCIP data structure */
      PRESOL*       presol              /**< the presolver itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** deinitialization method of presolver (called when problem solving exits) */
   virtual RETCODE scip_exit(
      SCIP*         scip,               /**< SCIP data structure */
      PRESOL*       presol              /**< the presolver itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** execution method of presolver
    *
    *  The presolver should go through the variables and constraints and tighten the domains or
    *  constraints. Each tightening should increase the given total numbers of changes.
    *
    *  possible return values for *result:
    *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
    *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
    *  - SCIP_SUCCESS    : the presolver found a reduction
    *  - SCIP_DIDNOTFIND : the presolver searched, but did not find a presolving change
    *  - SCIP_DIDNOTRUN  : the presolver was skipped
    */
   virtual RETCODE scip_exec(
      SCIP*         scip,               /**< SCIP data structure */
      PRESOL*       presol,             /**< the presolver itself */
      int           nrounds,            /**< no. of presolving rounds already done */
      int           nnewfixedvars,      /**< no. of variables fixed since last call to presolver */
      int           nnewaggrvars,       /**< no. of variables aggregated since last call to presolver */
      int           nnewchgvartypes,    /**< no. of variable type changes since last call to presolver */
      int           nnewchgbds,         /**< no. of variable bounds tightend since last call to presolver */
      int           nnewholes,          /**< no. of domain holes added since last call to presolver */
      int           nnewdelconss,       /**< no. of deleted constraints since last call to presolver */
      int           nnewupgdconss,      /**< no. of upgraded constraints since last call to presolver */
      int           nnewchgcoefs,       /**< no. of changed coefficients since last call to presolver */
      int           nnewchgsides,       /**< no. of changed left or right hand sides since last call to presolver */
      int*          nfixedvars,         /**< pointer to count total number of variables fixed of all presolvers */
      int*          naggrvars,          /**< pointer to count total number of variables aggregated of all presolvers */
      int*          nchgvartypes,       /**< pointer to count total number of variable type changes of all presolvers */
      int*          nchgbds,            /**< pointer to count total number of variable bounds tightend of all presolvers */
      int*          naddholes,          /**< pointer to count total number of domain holes added of all presolvers */
      int*          ndelconss,          /**< pointer to count total number of deleted constraints of all presolvers */
      int*          nupgdconss,         /**< pointer to count total number of upgraded constraints of all presolvers */
      int*          nchgcoefs,          /**< pointer to count total number of changed coefficients of all presolvers */
      int*          nchgsides,          /**< pointer to count total number of changed sides of all presolvers */
      RESULT*       result              /**< pointer to store the result of the presolving call */
      ) = 0;
};

} /* namespace scip */


   
/** creates the presolver for the given presolver object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       MyPresol* mypresol = new MyPresol(...);
 *       CHECK_OKAY( SCIPincludeObjPresol(scip, &mypresol, FALSE) );
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );
 *       delete mypresol;    // delete presol AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       CHECK_OKAY( SCIPincludeObjPresol(scip, new MyPresol(...), TRUE) );
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );  // destructor of MyPresol is called here
 */
extern
RETCODE SCIPincludeObjPresol(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjPresol* objpresol,          /**< presolver object */
   Bool             deleteobject        /**< should the presolver object be deleted when presolver is freed? */
   );

#endif

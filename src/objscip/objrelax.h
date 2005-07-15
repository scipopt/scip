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
#pragma ident "@(#) $Id: objrelax.h,v 1.10 2005/07/15 17:20:03 bzfpfend Exp $"

/**@file   objrelax.h
 * @brief  C++ wrapper for relaxators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJRELAX_H__
#define __SCIP_OBJRELAX_H__


extern "C" 
{
#include "scip/scip.h"
}


namespace scip
{

/** C++ wrapper object for relaxators */
class ObjRelax
{
public:
   /** name of the relaxator */
   const char* const scip_name_;
   
   /** description of the relaxator */
   const char* const scip_desc_;
   
   /** default priority of the relaxator (negative: call after LP, non-negative: call before LP) */
   const int scip_priority_;

   /** frequency for calling relaxator */
   const int scip_freq_;

   /** default constructor */
   ObjRelax(
      const char*   name,               /**< name of relaxator */
      const char*   desc,               /**< description of relaxator */
      int           priority,           /**< priority of the relaxator (negative: after LP, non-negative: before LP) */
      int           freq                /**< frequency for calling relaxator */
      )
      : scip_name_(name),
        scip_desc_(desc),
        scip_priority_(priority),
        scip_freq_(freq)
   {
   }

   /** destructor */
   virtual ~ObjRelax()
   {
   }

   /** destructor of relaxator to free user data (called when SCIP is exiting) */
   virtual RETCODE scip_free(
      SCIP*         scip,               /**< SCIP data structure */
      RELAX*        relax               /**< the relaxator itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** initialization method of relaxator (called after problem was transformed) */
   virtual RETCODE scip_init(
      SCIP*         scip,               /**< SCIP data structure */
      RELAX*        relax               /**< the relaxator itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** deinitialization method of relaxator (called before transformed problem is freed) */
   virtual RETCODE scip_exit(
      SCIP*         scip,               /**< SCIP data structure */
      RELAX*        relax               /**< the relaxator itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** solving process initialization method of relaxator (called when branch and bound process is about to begin)
    *
    *  This method is called when the presolving was finished and the branch and bound process is about to begin.
    *  The relaxator may use this call to initialize its branch and bound specific data.
    *
    */
   virtual RETCODE scip_initsol(
      SCIP*         scip,               /**< SCIP data structure */
      RELAX*        relax               /**< the relaxator itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** solving process deinitialization method of relaxator (called before branch and bound process data is freed)
    *
    *  This method is called before the branch and bound process is freed.
    *  The relaxator should use this call to clean up its branch and bound data.
    */
   virtual RETCODE scip_exitsol(
      SCIP*         scip,               /**< SCIP data structure */
      RELAX*        relax               /**< the relaxator itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** execution method of relaxator
    *
    *  The method is called in the node processing loop. It solves the current subproblem's relaxation.
    *  Like the LP relaxation, the relaxator should only operate on COLUMN variables.
    *
    *  possible return values for *result (if more than one applies, the first in the list should be used):
    *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
    *  - SCIP_CONSADDED  : an additional constraint was generated
    *  - SCIP_REDUCEDDOM : a variable's domain was reduced
    *  - SCIP_SEPARATED  : a cutting plane was generated
    *  - SCIP_SUCCESS    : the relaxator solved the relaxation and should not be called again on the same relaxation
    *  - SCIP_SUSPENDED  : the relaxator interrupted its solving process to wait for additional input (e.g. cutting
    *                      planes); however, it is able to continue the solving in order to improve the dual bound
    *  - SCIP_DIDNOTRUN  : the relaxator was skipped
    */
   virtual RETCODE scip_exec(
      SCIP*         scip,               /**< SCIP data structure */
      RELAX*        relax,              /**< the relaxator itself */
      RESULT*       result              /**< pointer to store the result of the relaxation call */
      ) = 0;
};

} /* namespace scip */


   
/** creates the relaxator for the given relaxator object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       MyRelax* myrelax = new MyRelax(...);
 *       CHECK_OKAY( SCIPincludeObjRelax(scip, &myrelax, FALSE) );
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );
 *       delete myrelax;    // delete relax AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       CHECK_OKAY( SCIPincludeObjRelax(scip, new MyRelax(...), TRUE) );
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );  // destructor of MyRelax is called here
 */
extern
RETCODE SCIPincludeObjRelax(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjRelax*  objrelax,           /**< relaxator object */
   Bool             deleteobject        /**< should the relaxator object be deleted when relaxator is freed? */
   );

/** returns the relax object of the given name, or NULL if not existing */
extern
scip::ObjRelax* SCIPfindObjRelax(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of relaxator */
   );

/** returns the relax object for the given relaxator */
extern
scip::ObjRelax* SCIPgetObjRelax(
   SCIP*            scip,               /**< SCIP data structure */
   RELAX*           relax               /**< relaxator */
   );

#endif

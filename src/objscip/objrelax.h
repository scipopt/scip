/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   objrelax.h
 * @brief  C++ wrapper for relaxation handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJRELAX_H__
#define __SCIP_OBJRELAX_H__

#include <cstring>

#include "scip/scip.h"
#include "objscip/objcloneable.h"

namespace scip
{

/** @brief C++ wrapper for relaxation handlers
 *
 *  This class defines the interface for relaxation handlers implemented in C++. Note that there is a pure virtual
 *  function (this function has to be implemented). This function is: scip_exec().
 *
 *  - \ref RELAX "Instructions for implementing a relaxation handler"
 *  - \ref type_relax.h "Corresponding C interface"
 */
class ObjRelax : public ObjCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the relaxator */
   char* scip_name_;
   
   /** description of the relaxator */
   char* scip_desc_;
   
   /** default priority of the relaxator (negative: call after LP, non-negative: call before LP) */
   const int scip_priority_;

   /** frequency for calling relaxator */
   const int scip_freq_;

   /** default constructor */
   ObjRelax(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of relaxator */
      const char*        desc,               /**< description of relaxator */
      int                priority,           /**< priority of the relaxator (negative: after LP, non-negative: before LP) */
      int                freq                /**< frequency for calling relaxator */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_priority_(priority),
        scip_freq_(freq)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** destructor */
   virtual ~ObjRelax()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** destructor of relaxator to free user data (called when SCIP is exiting) */
   virtual SCIP_RETCODE scip_free(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_RELAX*        relax               /**< the relaxator itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** initialization method of relaxator (called after problem was transformed) */
   virtual SCIP_RETCODE scip_init(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_RELAX*        relax               /**< the relaxator itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** deinitialization method of relaxator (called before transformed problem is freed) */
   virtual SCIP_RETCODE scip_exit(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_RELAX*        relax               /**< the relaxator itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** solving process initialization method of relaxator (called when branch and bound process is about to begin)
    *
    *  This method is called when the presolving was finished and the branch and bound process is about to begin.
    *  The relaxator may use this call to initialize its branch and bound specific data.
    *
    */
   virtual SCIP_RETCODE scip_initsol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_RELAX*        relax               /**< the relaxator itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** solving process deinitialization method of relaxator (called before branch and bound process data is freed)
    *
    *  This method is called before the branch and bound process is freed.
    *  The relaxator should use this call to clean up its branch and bound data.
    */
   virtual SCIP_RETCODE scip_exitsol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_RELAX*        relax               /**< the relaxator itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** execution method of relaxator
    *
    *  The method is called in the node processing loop. It solves the current subproblem's relaxation.
    *  Like the LP relaxation, the relaxator should only operate on COLUMN variables.
    *
    *  possible return values for *result (if more than one applies, the first in the list should be used):
    *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
    *  - SCIP_CONSADDED  : an additional constraint was generated, and the relaxator should not be called again on the
    *                      same relaxation
    *  - SCIP_REDUCEDDOM : a variable's domain was reduced, and the relaxator should not be called again on the same
    *                      relaxation
    *  - SCIP_SEPARATED  : a cutting plane was generated, and the relaxator should not be called again on the same
    *                      relaxation
    *  - SCIP_SUCCESS    : the relaxator solved the relaxation and should not be called again on the same relaxation
    *  - SCIP_SUSPENDED  : the relaxator interrupted its solving process to wait for additional input (e.g. cutting
    *                      planes); however, it is able to continue the solving in order to improve the dual bound
    *  - SCIP_DIDNOTRUN  : the relaxator was skipped
    */
   virtual SCIP_RETCODE scip_exec(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_RELAX*        relax,              /**< the relaxator itself */
      SCIP_Real*         lowerbound,         /**< pointer to store a lowerbound for the current node */
      SCIP_RESULT*       result              /**< pointer to store the result of the relaxation call */
      ) = 0;
};

} /* namespace scip */


   
/** creates the relaxator for the given relaxator object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyRelax* myrelax = new MyRelax(...);
 *       SCIP_CALL( SCIPincludeObjRelax(scip, &myrelax, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete myrelax;    // delete relax AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjRelax(scip, new MyRelax(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyRelax is called here
 */
extern
SCIP_RETCODE SCIPincludeObjRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjRelax*  objrelax,           /**< relaxator object */
   SCIP_Bool             deleteobject        /**< should the relaxator object be deleted when relaxator is freed? */
   );

/** returns the relax object of the given name, or 0 if not existing */
extern
scip::ObjRelax* SCIPfindObjRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of relaxator */
   );

/** returns the relax object for the given relaxator */
extern
scip::ObjRelax* SCIPgetObjRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax               /**< relaxator */
   );

#endif

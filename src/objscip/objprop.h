/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   objprop.h
 * @brief  C++ wrapper for propagators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJPROP_H__
#define __SCIP_OBJPROP_H__

#include <cstring>

#include "scip/scip.h"
#include "objscip/objcloneable.h"

namespace scip
{

/** C++ wrapper object for propagators */
class ObjProp : public ObjCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;
   
   /** name of the propagator */
   char* scip_name_;
   
   /** description of the propagator */
   char* scip_desc_;
   
   /** default priority of the propagator */
   const int scip_priority_;

   /** frequency for calling propagator */
   const int scip_freq_;

   /** should propagator be delayed, if other propagators found reductions? */
   const SCIP_Bool scip_delay_;

   /** default constructor */
   ObjProp(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of propagator */
      const char*        desc,               /**< description of propagator */
      int                priority,           /**< priority of the propagator */
      int                freq,               /**< frequency for calling propagator */
      SCIP_Bool          delay               /**< should propagator be delayed, if other propagators found reductions? */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_priority_(priority),
        scip_freq_(freq),
        scip_delay_(delay)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** destructor */
   virtual ~ObjProp()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** destructor of propagator to free user data (called when SCIP is exiting) */
   virtual SCIP_RETCODE scip_free(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_PROP*         prop                /**< the propagator itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** initialization method of propagator (called after problem was transformed) */
   virtual SCIP_RETCODE scip_init(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_PROP*         prop                /**< the propagator itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** deinitialization method of propagator (called before transformed problem is freed) */
   virtual SCIP_RETCODE scip_exit(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_PROP*         prop                /**< the propagator itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** solving process initialization method of propagator (called when branch and bound process is about to begin)
    *
    *  This method is called when the presolving was finished and the branch and bound process is about to begin.
    *  The propagator may use this call to initialize its branch and bound specific data.
    *
    */
   virtual SCIP_RETCODE scip_initsol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_PROP*         prop                /**< the propagator itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** solving process deinitialization method of propagator (called before branch and bound process data is freed)
    *
    *  This method is called before the branch and bound process is freed.
    *  The propagator should use this call to clean up its branch and bound data.
    */
   virtual SCIP_RETCODE scip_exitsol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_PROP*         prop                /**< the propagator itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** execution method of propagator
    *
    *  Searches for domain propagations. The method is called in the node processing loop.
    *
    *  possible return values for *result:
    *  - SCIP_CUTOFF     : the current node is infeasible for the current domains
    *  - SCIP_REDUCEDDOM : at least one domain reduction was found
    *  - SCIP_DIDNOTFIND : the propagator searched, but did not find a domain reduction
    *  - SCIP_DIDNOTRUN  : the propagator was skipped
    *  - SCIP_DELAYED    : the propagator was skipped, but should be called again
    */
   virtual SCIP_RETCODE scip_exec(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_PROP*         prop,               /**< the propagator itself */
      SCIP_RESULT*       result              /**< pointer to store the result of the propagation call */
      ) = 0;

   /** propagation conflict resolving method of propagator
    *
    *  This method is called during conflict analysis. If the propagator wants to support conflict analysis,
    *  it should call SCIPinferVarLbProp() or SCIPinferVarUbProp() in domain propagation instead of SCIPchgVarLb() or
    *  SCIPchgVarUb() in order to deduce bound changes on variables.
    *  In the SCIPinferVarLbProp() and SCIPinferVarUbProp() calls, the propagator provides a pointer to itself
    *  and an integer value "inferinfo" that can be arbitrarily chosen.
    *  The propagation conflict resolving method must then be implemented, to provide the "reasons" for the bound
    *  changes, i.e. the bounds of variables at the time of the propagation, that forced the propagator to set the
    *  conflict variable's bound to its current value. It can use the "inferinfo" tag to identify its own propagation
    *  rule and thus identify the "reason" bounds. The bounds that form the reason of the assignment must then be provided
    *  by calls to SCIPaddConflictLb() and SCIPaddConflictUb() in the propagation conflict resolving method.
    *
    *  See the description of the propagation conflict resulving method of constraint handlers for further details.
    *
    *  input:
    *  - scip            : SCIP main data structure
    *  - prop            : the propagator itself
    *  - infervar        : the conflict variable whose bound change has to be resolved
    *  - inferinfo       : the user information passed to the corresponding SCIPinferVarLbProp() or SCIPinferVarUbProp() call
    *  - boundtype       : the type of the changed bound (lower or upper bound)
    *  - bdchgidx        : the index of the bound change, representing the point of time where the change took place
    *
    *  output:
    *  - result          : pointer to store the result of the propagation conflict resolving call
    *
    *  possible return values for *result:
    *  - SCIP_SUCCESS    : the conflicting bound change has been successfully resolved by adding all reason bounds
    *  - SCIP_DIDNOTFIND : the conflicting bound change could not be resolved and has to be put into the conflict set
    */
   virtual SCIP_RETCODE scip_resprop(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_PROP*         prop,               /**< the propagator itself */
      SCIP_VAR*          infervar,           /**< the conflict variable whose bound change has to be resolved */
      int                inferinfo,          /**< the user information passed to the corresponding SCIPinferVarLbProp() or
                                              *    SCIPinferVarUbProp() call */
      SCIP_BOUNDTYPE     boundtype,          /**< the type of the changed bound (lower or upper bound) */
      SCIP_BDCHGIDX*     bdchgidx,           /**< the index of the bound change, representing the point of time where the
                                              *   change took place */
      SCIP_RESULT*       result              /**< pointer to store the result of the propagation call */
      ) = 0;
};

} /* namespace scip */


   
/** creates the propagator for the given propagator object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyProp* myprop = new MyProp(...);
 *       SCIP_CALL( SCIPincludeObjProp(scip, &myprop, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete myprop;    // delete prop AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjProp(scip, new MyProp(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyProp is called here
 */
extern
SCIP_RETCODE SCIPincludeObjProp(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjProp*        objprop,            /**< propagator object */
   SCIP_Bool             deleteobject        /**< should the propagator object be deleted when propagator is freed? */
   );

/** returns the prop object of the given name, or 0 if not existing */
extern
scip::ObjProp* SCIPfindObjProp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of propagator */
   );

/** returns the prop object for the given propagator */
extern
scip::ObjProp* SCIPgetObjProp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop                /**< propagator */
   );

#endif

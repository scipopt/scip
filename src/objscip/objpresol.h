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

/**@file   objpresol.h
 * @brief  C++ wrapper for presolvers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJPRESOL_H__
#define __SCIP_OBJPRESOL_H__

#include <cstring>

#include "scip/scip.h"
#include "objscip/objcloneable.h"

namespace scip
{

/** @brief C++ wrapper for presolvers
 *
 *  This class defines the interface for presolvers implemented in C++. Note that there is a pure virtual
 *  function (this function has to be implemented). This function is: scip_exec().
 *
 *  - \ref PRESOL "Instructions for implementing a presolver"
 *  - \ref PRESOLVERS "List of available presolvers"
 *  - \ref type_presol.h "Corresponding C interface"
 */
class ObjPresol : public ObjCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the presolver */
   char* scip_name_;
   
   /** description of the presolver */
   char* scip_desc_;
   
   /** default priority of the presolver */
   const int scip_priority_;

   /** default maximal number of presolving rounds the presolver participates in (-1: no limit) */
   const int scip_maxrounds_;

   /** should presolver be delayed, if other presolvers found reductions? */
   const SCIP_Bool scip_delay_;

   /** default constructor */
   ObjPresol(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of presolver */
      const char*        desc,               /**< description of presolver */
      int                priority,           /**< priority of the presolver */
      int                maxrounds,          /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
      SCIP_Bool          delay               /**< should presolver be delayed, if other presolvers found reductions? */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_priority_(priority),
        scip_maxrounds_(maxrounds),
        scip_delay_(delay)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** destructor */
   virtual ~ObjPresol()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** destructor of presolver to free user data (called when SCIP is exiting) */
   virtual SCIP_RETCODE scip_free(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_PRESOL*       presol              /**< the presolver itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** initialization method of presolver (called after problem was transformed) */
   virtual SCIP_RETCODE scip_init(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_PRESOL*       presol              /**< the presolver itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** deinitialization method of presolver (called before transformed problem is freed) */
   virtual SCIP_RETCODE scip_exit(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_PRESOL*       presol              /**< the presolver itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** presolving initialization method of presolver (called when presolving is about to begin)
    *
    *  possible return values for *result:
    *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
    *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
    *  - SCIP_FEASIBLE   : no infeasibility nor unboundness could be found
    */
   virtual SCIP_RETCODE scip_initpre(
      SCIP*              scip,               /**< SCIP data structure */   
      SCIP_PRESOL*       presol,             /**< presolver */
      SCIP_Bool          isunbounded,        /**< was unboundedness already detected */
      SCIP_Bool          isinfeasible,       /**< was infeasibility already detected */
      SCIP_RESULT*       result              /**< pointer to store the result of the callback method */
      )
   {  /*lint --e{715}*/
      assert(result != NULL);

      *result = SCIP_FEASIBLE;

      return SCIP_OKAY;
   }
   
   /** presolving deinitialization method of presolver (called after presolving has been finished)
    *
    *  possible return values for *result:
    *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
    *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
    *  - SCIP_FEASIBLE   : no infeasibility nor unboundness could be found
    */
   virtual SCIP_RETCODE scip_exitpre(
      SCIP*              scip,               /**< SCIP data structure */   
      SCIP_PRESOL*       presol,             /**< presolver */
      SCIP_Bool          isunbounded,        /**< was unboundedness already detected */
      SCIP_Bool          isinfeasible,       /**< was infeasibility already detected */
      SCIP_RESULT*       result              /**< pointer to store the result of the callback method */
      )
   {  /*lint --e{715}*/
      assert(result != NULL);

      *result = SCIP_FEASIBLE;

      return SCIP_OKAY;
   }

   /** execution method of presolver
    *
    *  The presolver should go through the variables and constraints and tighten the domains or
    *  constraints. Each tightening should increase the given total numbers of changes.
    *
    *  @note the counters state the changes since the last call including the changes of this presolver during its last
    *        last call
    *
    *  possible return values for *result:
    *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
    *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
    *  - SCIP_SUCCESS    : the presolver found a reduction
    *  - SCIP_DIDNOTFIND : the presolver searched, but did not find a presolving change
    *  - SCIP_DIDNOTRUN  : the presolver was skipped
    *  - SCIP_DELAYED    : the presolver should be called again after all (none delayed) wants
    */
   virtual SCIP_RETCODE scip_exec(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_PRESOL*       presol,             /**< the presolver itself */
      int                nrounds,            /**< no. of presolving rounds already done */
      int                nnewfixedvars,      /**< no. of variables fixed since last call to presolver */
      int                nnewaggrvars,       /**< no. of variables aggregated since last call to presolver */
      int                nnewchgvartypes,    /**< no. of variable type changes since last call to presolver */
      int                nnewchgbds,         /**< no. of variable bounds tightend since last call to presolver */
      int                nnewholes,          /**< no. of domain holes added since last call to presolver */
      int                nnewdelconss,       /**< no. of deleted constraints since last call to presolver */
      int                nnewaddconss,       /**< no. of added constraints since last call to presolver */
      int                nnewupgdconss,      /**< no. of upgraded constraints since last call to presolver */
      int                nnewchgcoefs,       /**< no. of changed coefficients since last call to presolver */
      int                nnewchgsides,       /**< no. of changed left or right hand sides since last call to presolver */
      int*               nfixedvars,         /**< pointer to count total number of variables fixed of all presolvers */
      int*               naggrvars,          /**< pointer to count total number of variables aggregated of all presolvers */
      int*               nchgvartypes,       /**< pointer to count total number of variable type changes of all presolvers */
      int*               nchgbds,            /**< pointer to count total number of variable bounds tightend of all presolvers */
      int*               naddholes,          /**< pointer to count total number of domain holes added of all presolvers */
      int*               ndelconss,          /**< pointer to count total number of deleted constraints of all presolvers */
      int*               naddconss,          /**< pointer to count total number of added constraints of all presolvers */
      int*               nupgdconss,         /**< pointer to count total number of upgraded constraints of all presolvers */
      int*               nchgcoefs,          /**< pointer to count total number of changed coefficients of all presolvers */
      int*               nchgsides,          /**< pointer to count total number of changed sides of all presolvers */
      SCIP_RESULT*       result              /**< pointer to store the result of the presolving call */
      ) = 0;
};

} /* namespace scip */


   
/** creates the presolver for the given presolver object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyPresol* mypresol = new MyPresol(...);
 *       SCIP_CALL( SCIPincludeObjPresol(scip, &mypresol, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete mypresol;    // delete presol AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjPresol(scip, new MyPresol(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyPresol is called here
 */
extern
SCIP_RETCODE SCIPincludeObjPresol(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjPresol*      objpresol,          /**< presolver object */
   SCIP_Bool             deleteobject        /**< should the presolver object be deleted when presolver is freed? */
   );

/** returns the presol object of the given name, or 0 if not existing */
extern
scip::ObjPresol* SCIPfindObjPresol(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of presolver */
   );

/** returns the presol object for the given presolver */
extern
scip::ObjPresol* SCIPgetObjPresol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol              /**< presolver */
   );

#endif

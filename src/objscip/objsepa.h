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

/**@file   objsepa.h
 * @brief  C++ wrapper for cut separators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJSEPA_H__
#define __SCIP_OBJSEPA_H__

#include <cstring>

#include "scip/scip.h"
#include "objscip/objcloneable.h"

namespace scip
{

/** C++ wrapper object for cut separators */
class ObjSepa : public ObjCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the cut separator */
   char* scip_name_;
   
   /** description of the cut separator */
   char* scip_desc_;
   
   /** default priority of the cut separator */
   const int scip_priority_;

   /** frequency for calling separator */
   const int scip_freq_;

   /** maximal relative distance from current node's dual bound to primal bound compared to best node's dual bound for applying
    *  separation (0.0: only on current best node, 1.0: on all nodes)
    */
   const SCIP_Real scip_maxbounddist_;

   /** does the separator use a secondary SCIP instance? */
   const SCIP_Bool scip_usessubscip_;

   /** should separator be delayed, if other separators found cuts? */
   const SCIP_Bool scip_delay_;

   /** default constructor */
   ObjSepa(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of cut separator */
      const char*        desc,               /**< description of cut separator */
      int                priority,           /**< priority of the cut separator */
      int                freq,               /**< frequency for calling separator */
      SCIP_Real          maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for applying separation */
      SCIP_Bool          usessubscip,        /**< does the separator use a secondary SCIP instance? */
      SCIP_Bool          delay               /**< should separator be delayed, if other separators found cuts? */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_priority_(priority),
        scip_freq_(freq),
        scip_maxbounddist_(maxbounddist),
        scip_usessubscip_(usessubscip),
        scip_delay_(delay)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** destructor */
   virtual ~ObjSepa()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** destructor of cut separator to free user data (called when SCIP is exiting) */
   virtual SCIP_RETCODE scip_free(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_SEPA*         sepa                /**< the cut separator itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** initialization method of cut separator (called after problem was transformed) */
   virtual SCIP_RETCODE scip_init(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_SEPA*         sepa                /**< the cut separator itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** deinitialization method of cut separator (called before transformed problem is freed) */
   virtual SCIP_RETCODE scip_exit(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_SEPA*         sepa                /**< the cut separator itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** solving process initialization method of separator (called when branch and bound process is about to begin)
    *
    *  This method is called when the presolving was finished and the branch and bound process is about to begin.
    *  The separator may use this call to initialize its branch and bound specific data.
    */
   virtual SCIP_RETCODE scip_initsol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_SEPA*         sepa                /**< the cut separator itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** solving process deinitialization method of separator (called before branch and bound process data is freed)
    *
    *  This method is called before the branch and bound process is freed.
    *  The separator should use this call to clean up its branch and bound data.
    */
   virtual SCIP_RETCODE scip_exitsol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_SEPA*         sepa                /**< the cut separator itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** LP solution separation method of separator
    *
    *  Searches for cutting planes that separate the current LP solution. The method is called in the LP solving loop,
    *  which means that a valid LP solution exists.
    *
    *  possible return values for *result:
    *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
    *  - SCIP_SEPARATED  : a cutting plane was generated
    *  - SCIP_REDUCEDDOM : no cutting plane was generated, but a variable's domain was reduced
    *  - SCIP_CONSADDED  : no cutting plane or domain reduction, but an additional constraint was generated
    *  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
    *  - SCIP_DIDNOTRUN  : the separator was skipped
    *  - SCIP_DELAYED    : the separator was skipped, but should be called again
    */
   virtual SCIP_RETCODE scip_execlp(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_SEPA*         sepa,               /**< the cut separator itself */
      SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
      )
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
   
   /** arbitrary primal solution separation method of separator
    *
    *  Searches for cutting planes that separate the given primal solution. The method is called outside the LP solution
    *  loop (e.g., by a relaxator or a primal heuristic), which means that there is no valid LP solution.
    *
    *  possible return values for *result:
    *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
    *  - SCIP_SEPARATED  : a cutting plane was generated
    *  - SCIP_REDUCEDDOM : no cutting plane was generated, but a variable's domain was reduced
    *  - SCIP_CONSADDED  : no cutting plane or domain reduction, but an additional constraint was generated
    *  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
    *  - SCIP_DIDNOTRUN  : the separator was skipped
    *  - SCIP_DELAYED    : the separator was skipped, but should be called again
    */
   virtual SCIP_RETCODE scip_execsol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_SEPA*         sepa,               /**< the cut separator itself */
      SCIP_SOL*          sol,                /**< primal solution that should be separated */
      SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
      )
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
};

} /* namespace scip */


   
/** creates the cut separator for the given cut separator object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MySepa* mysepa = new MySepa(...);
 *       SCIP_CALL( SCIPincludeObjSepa(scip, &mysepa, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete mysepa;    // delete sepa AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjSepa(scip, new MySepa(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MySepa is called here
 */
extern
SCIP_RETCODE SCIPincludeObjSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjSepa*        objsepa,            /**< cut separator object */
   SCIP_Bool             deleteobject        /**< should the cut separator object be deleted when cut separator is freed? */
   );

/** returns the sepa object of the given name, or 0 if not existing */
extern
scip::ObjSepa* SCIPfindObjSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of cut separator */
   );

/** returns the sepa object for the given cut separator */
extern
scip::ObjSepa* SCIPgetObjSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa                /**< cut separator */
   );

#endif

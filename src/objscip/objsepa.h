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
#pragma ident "@(#) $Id: objsepa.h,v 1.12 2005/02/08 14:22:28 bzfpfend Exp $"

/**@file   objsepa.h
 * @brief  C++ wrapper for cut separators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __OBJSEPA_H__
#define __OBJSEPA_H__


extern "C" 
{
#include "scip.h"
}


namespace scip
{

/** C++ wrapper object for cut separators */
class ObjSepa
{
public:
   /** name of the cut separator */
   const char* const scip_name_;
   
   /** description of the cut separator */
   const char* const scip_desc_;
   
   /** default priority of the cut separator */
   const int scip_priority_;

   /** frequency for calling separator */
   const int scip_freq_;

   /** should separator be delayed, if other separators found cuts? */
   const Bool scip_delay_;

   /** default constructor */
   ObjSepa(
      const char*   name,               /**< name of cut separator */
      const char*   desc,               /**< description of cut separator */
      int           priority,           /**< priority of the cut separator */
      int           freq,               /**< frequency for calling separator */
      Bool          delay               /**< should separator be delayed, if other separators found cuts? */
      )
      : scip_name_(name),
        scip_desc_(desc),
        scip_priority_(priority),
        scip_freq_(freq),
        scip_delay_(delay)
   {
   }

   /** destructor */
   virtual ~ObjSepa()
   {
   }

   /** destructor of cut separator to free user data (called when SCIP is exiting) */
   virtual RETCODE scip_free(
      SCIP*         scip,               /**< SCIP data structure */
      SEPA*         sepa                /**< the cut separator itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** initialization method of cut separator (called after problem was transformed) */
   virtual RETCODE scip_init(
      SCIP*         scip,               /**< SCIP data structure */
      SEPA*         sepa                /**< the cut separator itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** deinitialization method of cut separator (called before transformed problem is freed) */
   virtual RETCODE scip_exit(
      SCIP*         scip,               /**< SCIP data structure */
      SEPA*         sepa                /**< the cut separator itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** solving process initialization method of separator (called when branch and bound process is about to begin)
    *
    *  This method is called when the presolving was finished and the branch and bound process is about to begin.
    *  The separator may use this call to initialize its branch and bound specific data.
    */
   virtual RETCODE scip_initsol(
      SCIP*         scip,               /**< SCIP data structure */
      SEPA*         sepa                /**< the cut separator itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** solving process deinitialization method of separator (called before branch and bound process data is freed)
    *
    *  This method is called before the branch and bound process is freed.
    *  The separator should use this call to clean up its branch and bound data.
    */
   virtual RETCODE scip_exitsol(
      SCIP*         scip,               /**< SCIP data structure */
      SEPA*         sepa                /**< the cut separator itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** execution method of separator
    *
    *  Searches for cutting planes. The method is called in the LP solving loop.
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
   virtual RETCODE scip_exec(
      SCIP*         scip,               /**< SCIP data structure */
      SEPA*         sepa,               /**< the cut separator itself */
      RESULT*       result              /**< pointer to store the result of the separation call */
      ) = 0;
};

} /* namespace scip */


   
/** creates the cut separator for the given cut separator object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       MySepa* mysepa = new MySepa(...);
 *       CHECK_OKAY( SCIPincludeObjSepa(scip, &mysepa, FALSE) );
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );
 *       delete mysepa;    // delete sepa AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       CHECK_OKAY( SCIPincludeObjSepa(scip, new MySepa(...), TRUE) );
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );  // destructor of MySepa is called here
 */
extern
RETCODE SCIPincludeObjSepa(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjSepa*   objsepa,            /**< cut separator object */
   Bool             deleteobject        /**< should the cut separator object be deleted when cut separator is freed? */
   );

/** returns the sepa object of the given name, or NULL if not existing */
extern
scip::ObjSepa* SCIPfindObjSepa(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of cut separator */
   );

/** returns the sepa object for the given cut separator */
extern
scip::ObjSepa* SCIPgetObjSepa(
   SCIP*            scip,               /**< SCIP data structure */
   SEPA*            sepa                /**< cut separator */
   );

#endif

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
#pragma ident "@(#) $Id: objbranchrule.h,v 1.15 2005/01/18 09:26:49 bzfpfend Exp $"

/**@file   objbranchrule.h
 * @brief  C++ wrapper for branching rules
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __OBJBRANCHRULE_H__
#define __OBJBRANCHRULE_H__


#include <cassert>

extern "C" 
{
#include "scip.h"
}


namespace scip
{

/** C++ wrapper object for branching rules */
class ObjBranchrule
{
public:
   /** name of the branching rule */
   const char* const scip_name_;
   
   /** description of the branching rule */
   const char* const scip_desc_;
   
   /** default priority of the branching rule */
   const int scip_priority_;

   /** default maximal depth for applying the branching rule */
   const int scip_maxdepth_;

   /** default maximal relative distance from current node's dual bound to primal bound
    *  compared to best node's dual bound for applying branching rule
    *  (0.0: only on current best node, 1.0: on all nodes)
    */
   const Real scip_maxbounddist_;

   /** default constructor */
   ObjBranchrule(
      const char*   name,               /**< name of branching rule */
      const char*   desc,               /**< description of branching rule */
      int           priority,           /**< priority of the branching rule */
      int           maxdepth,           /**< maximal depth level, up to which this branching rule should be used (or -1) */
      Real          maxbounddist        /**< maximal relative distance from current node's dual bound to primal bound
                                         *   compared to best node's dual bound for applying branching rule
                                         *   (0.0: only on current best node, 1.0: on all nodes) */
      )
      : scip_name_(name),
        scip_desc_(desc),
        scip_priority_(priority),
        scip_maxdepth_(maxdepth),
        scip_maxbounddist_(maxbounddist)
   {
   }

   /** destructor */
   virtual ~ObjBranchrule()
   {
   }

   /** destructor of branching rule to free user data (called when SCIP is exiting) */
   virtual RETCODE scip_free(
      SCIP*         scip,               /**< SCIP data structure */
      BRANCHRULE*   branchrule          /**< the branching rule itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** initialization method of branching rule (called after problem was transformed) */
   virtual RETCODE scip_init(
      SCIP*         scip,               /**< SCIP data structure */
      BRANCHRULE*   branchrule          /**< the branching rule itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** deinitialization method of branching rule (called before transformed problem is freed) */
   virtual RETCODE scip_exit(
      SCIP*         scip,               /**< SCIP data structure */
      BRANCHRULE*   branchrule          /**< the branching rule itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
   virtual RETCODE scip_initsol(
      SCIP*         scip,               /**< SCIP data structure */
      BRANCHRULE*   branchrule          /**< the branching rule itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
   virtual RETCODE scip_exitsol(
      SCIP*         scip,               /**< SCIP data structure */
      BRANCHRULE*   branchrule          /**< the branching rule itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** branching execution method for fractional LP solutions
    *
    *  possible return values for *result (if more than one applies, the first in the list should be used):
    *  - SCIP_CUTOFF     : the current node was detected to be infeasible
    *  - SCIP_CONSADDED  : an additional constraint (e.g. a conflict clause) was generated; this result code must not be
    *                      returned, if allowaddcons is FALSE
    *  - SCIP_REDUCEDDOM : a domain was reduced that rendered the current LP solution infeasible
    *  - SCIP_SEPARATED  : a cutting plane was generated
    *  - SCIP_BRANCHED   : branching was applied
    *  - SCIP_DIDNOTRUN  : the branching rule was skipped
    */
   virtual RETCODE scip_execlp(
      SCIP*         scip,               /**< SCIP data structure */
      BRANCHRULE*   branchrule,         /**< the branching rule itself */
      Bool          allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
      RESULT*       result              /**< pointer to store the result of the branching call */
      )
   {
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
   
   /** branching execution method for not completely fixed pseudo solutions
    *
    *  possible return values for *result (if more than one applies, the first in the list should be used):
    *  - SCIP_CUTOFF     : the current node was detected to be infeasible
    *  - SCIP_CONSADDED  : an additional constraint (e.g. a conflict clause) was generated; this result code must not be
    *                      returned, if allowaddcons is FALSE
    *  - SCIP_REDUCEDDOM : a domain was reduced that rendered the current pseudo solution infeasible
    *  - SCIP_BRANCHED   : branching was applied
    *  - SCIP_DIDNOTRUN  : the branching rule was skipped
    */
   virtual RETCODE scip_execps(
      SCIP*         scip,               /**< SCIP data structure */
      BRANCHRULE*   branchrule,         /**< the branching rule itself */
      Bool          allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
      RESULT*       result              /**< pointer to store the result of the branching call */
      )
   {
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
};

} /* namespace scip */


   
/** creates the branching rule for the given branching rule object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       MyBranchrule* mybranchrule = new MyBranchrule(...);
 *       CHECK_OKAY( SCIPincludeObjBranchrule(scip, &mybranchrule, FALSE) );
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );
 *       delete mybranchrule;    // delete branchrule AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       CHECK_OKAY( SCIPincludeObjBranchrule(scip, new MyBranchrule(...), TRUE) );
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );  // destructor of MyBranchrule is called here
 */
extern
RETCODE SCIPincludeObjBranchrule(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjBranchrule* objbranchrule,  /**< branching rule object */
   Bool             deleteobject        /**< should the branching rule object be deleted when branching rule is freed? */
   );

/** returns the branchrule object of the given name, or NULL if not existing */
extern
scip::ObjBranchrule* SCIPfindObjBranchrule(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of branching rule */
   );

/** returns the branchrule object for the given branching rule */
extern
scip::ObjBranchrule* SCIPgetObjBranchrule(
   SCIP*            scip,               /**< SCIP data structure */
   BRANCHRULE*      branchrule          /**< branching rule */
   );

#endif

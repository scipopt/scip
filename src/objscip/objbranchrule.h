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
#pragma ident "@(#) $Id: objbranchrule.h,v 1.3 2003/12/08 11:51:03 bzfpfend Exp $"

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

   /** default constructor */
   ObjBranchrule(
      const char*   name,               /**< name of branching rule */
      const char*   desc,               /**< description of branching rule */
      int           priority            /**< priority of the branching rule */
      )
      : scip_name_(name),
        scip_desc_(desc),
        scip_priority_(priority)
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
   
   /** initialization method of branching rule (called when problem solving starts) */
   virtual RETCODE scip_init(
      SCIP*         scip,               /**< SCIP data structure */
      BRANCHRULE*   branchrule          /**< the branching rule itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** deinitialization method of branching rule (called when problem solving exits) */
   virtual RETCODE scip_exit(
      SCIP*         scip,               /**< SCIP data structure */
      BRANCHRULE*   branchrule          /**< the branching rule itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** branching execution method for fractional LP solutions */
   virtual RETCODE scip_execlp(
      SCIP*         scip,               /**< SCIP data structure */
      BRANCHRULE*   branchrule,         /**< the branching rule itself */
      RESULT*       result              /**< pointer to store the result of the branching call */
      )
   {
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
   
   /** branching execution method for not completely fixed pseudo solutions */
   virtual RETCODE scip_execps(
      SCIP*         scip,               /**< SCIP data structure */
      BRANCHRULE*   branchrule,         /**< the branching rule itself */
      RESULT*       result              /**< pointer to store the result of the branching call */
      )
   {
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
};

} /* namespace scip */


   
/** creates the branching rule for the given branching rule object and includes it in SCIP */
extern
RETCODE SCIPincludeObjBranchrule(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjBranchrule* objbranchrule,  /**< branching rule object */
   Bool             deleteobject        /**< should the branching rule object be deleted when branching rule is freed? */
   );

#endif

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
#pragma ident "@(#) $Id: objpricer.h,v 1.1 2003/11/27 17:48:44 bzfpfend Exp $"

/**@file   objpricer.h
 * @brief  C++ wrapper for variable pricers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __OBJPRICER_H__
#define __OBJPRICER_H__


extern "C" 
{
#include "scip.h"
}


namespace scip
{

/** C++ wrapper object for variable pricers */
class ObjPricer
{
public:
   /** name of the variable pricer */
   const char* const scip_name_;
   
   /** description of the variable pricer */
   const char* const scip_desc_;
   
   /** default priority of the variable pricer */
   const int scip_priority_;

   /** default constructor */
   ObjPricer(
      const char*   name,               /**< name of variable pricer */
      const char*   desc,               /**< description of variable pricer */
      int           priority            /**< priority of the variable pricer */
      )
      : scip_name_(name),
        scip_desc_(desc),
        scip_priority_(priority)
   {
   }

   /** destructor of variable pricer to free user data (called when SCIP is exiting) */
   virtual RETCODE scip_free(
      SCIP*         scip,               /**< SCIP data structure */
      PRICER*       pricer              /**< the variable pricer itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** initialization method of variable pricer (called when problem solving starts) */
   virtual RETCODE scip_init(
      SCIP*         scip,               /**< SCIP data structure */
      PRICER*       pricer              /**< the variable pricer itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** deinitialization method of variable pricer (called when problem solving exits) */
   virtual RETCODE scip_exit(
      SCIP*         scip,               /**< SCIP data structure */
      PRICER*       pricer              /**< the variable pricer itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** reduced cost pricing method of variable pricer for feasible LPs */
   virtual RETCODE scip_redcost(
      SCIP*         scip,               /**< SCIP data structure */
      PRICER*       pricer              /**< the variable pricer itself */
      ) = 0;
   
   /** farkas pricing method of variable pricer for infeasible LPs */
   virtual RETCODE scip_farkas(
      SCIP*         scip,               /**< SCIP data structure */
      PRICER*       pricer              /**< the variable pricer itself */
      )
   {
      return SCIP_OKAY;
   }
};

} /* namespace scip */


   
/** creates the variable pricer for the given variable pricer object and includes it in SCIP */
RETCODE SCIPincludeObjPricer(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjPricer* objpricer           /**< variable pricer object */
   );

#endif

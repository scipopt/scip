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
#pragma ident "@(#) $Id: objpresol.h,v 1.2 2003/12/08 11:51:03 bzfpfend Exp $"

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
   
   /** execution method of presolver */
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


   
/** creates the presolver for the given presolver object and includes it in SCIP */
extern
RETCODE SCIPincludeObjPresol(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjPresol* objpresol,          /**< presolver object */
   Bool             deleteobject        /**< should the presolver object be deleted when presolver is freed? */
   );

#endif

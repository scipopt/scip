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
#pragma ident "@(#) $Id: objsepa.h,v 1.2 2003/12/08 11:51:04 bzfpfend Exp $"

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

   /** default constructor */
   ObjSepa(
      const char*   name,               /**< name of cut separator */
      const char*   desc,               /**< description of cut separator */
      int           priority,           /**< priority of the cut separator */
      int           freq                /**< frequency for calling separator */
      )
      : scip_name_(name),
        scip_desc_(desc),
        scip_priority_(priority),
        scip_freq_(freq)
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
   
   /** initialization method of cut separator (called when problem solving starts) */
   virtual RETCODE scip_init(
      SCIP*         scip,               /**< SCIP data structure */
      SEPA*         sepa                /**< the cut separator itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** deinitialization method of cut separator (called when problem solving exits) */
   virtual RETCODE scip_exit(
      SCIP*         scip,               /**< SCIP data structure */
      SEPA*         sepa                /**< the cut separator itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** execution method of separator */
   virtual RETCODE scip_exec(
      SCIP*         scip,               /**< SCIP data structure */
      SEPA*         sepa,               /**< the cut separator itself */
      RESULT*       result              /**< pointer to store the result of the separation call */
      ) = 0;
};

} /* namespace scip */


   
/** creates the cut separator for the given cut separator object and includes it in SCIP */
extern
RETCODE SCIPincludeObjSepa(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjSepa*   objsepa,            /**< cut separator object */
   Bool             deleteobject        /**< should the cut separator object be deleted when cut separator is freed? */
   );

#endif

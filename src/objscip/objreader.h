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
#pragma ident "@(#) $Id: objreader.h,v 1.1 2003/11/28 10:05:47 bzfpfend Exp $"

/**@file   objreader.h
 * @brief  C++ wrapper for file readers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __OBJREADER_H__
#define __OBJREADER_H__


extern "C" 
{
#include "scip.h"
}


namespace scip
{

/** C++ wrapper object for file readers */
class ObjReader
{
public:
   /** name of the file reader */
   const char* const scip_name_;
   
   /** description of the file reader */
   const char* const scip_desc_;
   
   /** file extension that reader processes */
   const char* const scip_extension_;

   /** default constructor */
   ObjReader(
      const char*   name,               /**< name of file reader */
      const char*   desc,               /**< description of file reader */
      const char*   extension           /**< file extension that reader processes */
      )
      : scip_name_(name),
        scip_desc_(desc),
        scip_extension_(extension)
   {
   }

   /** destructor of file reader to free user data (called when SCIP is exiting) */
   virtual RETCODE scip_free(
      SCIP*         scip,               /**< SCIP data structure */
      READER*       reader              /**< the file reader itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** problem reading method of reader */
   virtual RETCODE scip_read(
      SCIP*         scip,               /**< SCIP data structure */
      READER*       reader,             /**< the file reader itself */
      const char*   filename,           /**< full path and name of file to read, or NULL if stdin should be used */
      RESULT*       result              /**< pointer to store the result of the file reading call */
      ) = 0;
};

} /* namespace scip */


   
/** creates the file reader for the given file reader object and includes it in SCIP */
RETCODE SCIPincludeObjReader(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjReader* objreader           /**< file reader object */
   );

#endif

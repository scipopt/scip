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
#pragma ident "@(#) $Id: objreader.h,v 1.9 2005/02/14 13:35:38 bzfpfend Exp $"

/**@file   objreader.h
 * @brief  C++ wrapper for file readers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __OBJREADER_H__
#define __OBJREADER_H__


extern "C" 
{
#include "scip/scip.h"
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

   /** destructor */
   virtual ~ObjReader()
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
   
   /** problem reading method of reader
    *
    *  possible return values for *result:
    *  - SCIP_SUCCESS    : the reader read the file correctly and created an appropritate problem
    *  - SCIP_DIDNOTRUN  : the reader is not responsible for given input file
    *
    *  If the reader detected an error in the input file, it should return with RETCODE SCIP_READERR or SCIP_NOFILE.
    */
   virtual RETCODE scip_read(
      SCIP*         scip,               /**< SCIP data structure */
      READER*       reader,             /**< the file reader itself */
      const char*   filename,           /**< full path and name of file to read, or NULL if stdin should be used */
      RESULT*       result              /**< pointer to store the result of the file reading call */
      ) = 0;
};

} /* namespace scip */


   
/** creates the file reader for the given file reader object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       MyReader* myreader = new MyReader(...);
 *       CHECK_OKAY( SCIPincludeObjReader(scip, &myreader, FALSE) );
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );
 *       delete myreader;    // delete reader AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       CHECK_OKAY( SCIPincludeObjReader(scip, new MyReader(...), TRUE) );
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );  // destructor of MyReader is called here
 */
extern
RETCODE SCIPincludeObjReader(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjReader* objreader,          /**< file reader object */
   Bool             deleteobject        /**< should the reader object be deleted when reader is freed? */
   );

/** returns the reader object of the given name, or NULL if not existing */
extern
scip::ObjReader* SCIPfindObjReader(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of file reader */
   );

/** returns the reader object for the given file reader */
extern
scip::ObjReader* SCIPgetObjReader(
   SCIP*            scip,               /**< SCIP data structure */
   READER*          reader              /**< file reader */
   );

#endif

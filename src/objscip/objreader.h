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

/**@file   objreader.h
 * @brief  C++ wrapper for file readers and writers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJREADER_H__
#define __SCIP_OBJREADER_H__

#include <cstring>

#include "scip/scip.h"
#include "objscip/objcloneable.h"

namespace scip
{

/** @brief C++ wrapper for file readers and writers
 *
 *  This class defines the interface for file readers and writers implemented in C++.
 *
 *  - \ref READER "Instructions for implementing a file reader and writer"
 *  - \ref FILEREADERS "List of available file readers and writers"
 */
class ObjReader : public ObjCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the file reader */
   char* scip_name_;
   
   /** description of the file reader */
   char* scip_desc_;
   
   /** file extension that reader processes */
   char* scip_extension_;

   /** default constructor */
   ObjReader(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of file reader */
      const char*        desc,               /**< description of file reader */
      const char*        extension           /**< file extension that reader processes */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_extension_(0)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_extension_, extension, std::strlen(extension)+1) );
   }

   /** destructor */
   virtual ~ObjReader()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
      SCIPfreeMemoryArray(scip_, &scip_extension_);
   }

   /** destructor of file reader to free user data (called when SCIP is exiting) */
   virtual SCIP_RETCODE scip_free(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_READER*       reader              /**< the file reader itself */
      )
   {  /*lint --e{715}*/
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
   virtual SCIP_RETCODE scip_read(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_READER*       reader,             /**< the file reader itself */
      const char*        filename,           /**< full path and name of file to read, or NULL if stdin should be used */
      SCIP_RESULT*       result              /**< pointer to store the result of the file reading call */
      )
   {  /*lint --e{715}*/

      /* set result pointer to indicate that the reading was not performed */
      assert(result != NULL);
      (*result) = SCIP_DIDNOTRUN;

      return SCIP_OKAY;
   }

   /** problem writing method of reader; NOTE: if the parameter "genericnames" is TRUE, then
    *  SCIP already set all variable and constraint names to generic names; therefore, this
    *  method should always use SCIPvarGetName() and SCIPconsGetName(); 
    *
    *  possible return values for *result:
    *  - SCIP_SUCCESS    : the reader read the file correctly and created an appropritate problem
    *  - SCIP_DIDNOTRUN  : the reader is not responsible for given input file
    *
    *  If the reader detected an error in the writing to the file stream, it should return
    *  with RETCODE SCIP_WRITEERROR.
    */
   virtual SCIP_RETCODE scip_write(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_READER*       reader,             /**< the file reader itself */
      FILE*              file,               /**< output file, or NULL if standard output should be used */
      const char*        name,               /**< problem name */
      SCIP_PROBDATA*     probdata,           /**< user problem data set by the reader */
      SCIP_Bool          transformed,        /**< TRUE iff problem is the transformed problem */

      SCIP_OBJSENSE      objsense,           /**< objective sense */
      SCIP_Real          objscale,           /**< scalar applied to objective function; external objective value is
                                              *   extobj = objsense * objscale * (intobj + objoffset) */
      SCIP_Real          objoffset,          /**< objective offset from bound shifting and fixing */
      SCIP_VAR**         vars,               /**< array with active variables ordered binary, integer, implicit, 
                                              *   continuous */
      int                nvars,              /**< number of mutable variables in the problem */
      int                nbinvars,           /**< number of binary variables */
      int                nintvars,           /**< number of general integer variables */
      int                nimplvars,          /**< number of implicit integer variables */
      int                ncontvars,          /**< number of continuous variables */
      SCIP_VAR**         fixedvars,          /**< array with fixed and aggregated variables */
      int                nfixedvars,         /**< number of fixed and aggregated variables in the problem */
      int                startnvars,         /**< number of variables existing when problem solving started */
      SCIP_CONS**        conss,              /**< array with constraints of the problem */
      int                nconss,             /**< number of constraints in the problem */
      int                maxnconss,          /**< maximum number of constraints existing at the same time */
      int                startnconss,        /**< number of constraints existing when problem solving started */
      SCIP_Bool          genericnames,       /**< using generic variable and constraint names? */
      SCIP_RESULT*       result              /**< pointer to store the result of the file reading call */
      )
   {  /*lint --e{715}*/

      /* set result pointer to indicate that the writing was not performed */
      assert(result != NULL);
      (*result) = SCIP_DIDNOTRUN;

      return SCIP_OKAY;
   }
};
   
} /* namespace scip */


   
/** creates the file reader for the given file reader object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyReader* myreader = new MyReader(...);
 *       SCIP_CALL( SCIPincludeObjReader(scip, &myreader, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete myreader;    // delete reader AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjReader(scip, new MyReader(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyReader is called here
 */
extern
SCIP_RETCODE SCIPincludeObjReader(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjReader*      objreader,          /**< file reader object */
   SCIP_Bool             deleteobject        /**< should the reader object be deleted when reader is freed? */
   );

/** returns the reader object of the given name, or 0 if not existing */
extern
scip::ObjReader* SCIPfindObjReader(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of file reader */
   );

/** returns the reader object for the given file reader */
extern
scip::ObjReader* SCIPgetObjReader(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader              /**< file reader */
   );

#endif

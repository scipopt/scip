/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader.h
 * @brief  interface for input file readers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __READER_H__
#define __READER_H__


#include <stdio.h>


typedef struct Reader READER;               /**< reader data structure */
typedef struct ReaderData READERDATA;       /**< reader specific data */


/** destructor of reader to free user data (called when SCIP is exiting)
 *
 *  input:
 *    reader          : the reader itself
 *    scip            : SCIP main data structure
 */
#define DECL_READERFREE(x) RETCODE x (READER* reader, SCIP* scip)

/** problem reading method of reader
 *
 *  input:
 *    reader          : the reader itself
 *    scip            : SCIP main data structure
 *    filename        : full path and name of file to read, or NULL if stdin should be used
 *    result          : pointer to store the result of the propagation call
 *
 *  possible return values for *result:
 *    SCIP_SUCCESS    : the reader read the file correctly
 *    SCIP_DIDNOTRUN  : the reader is not responsible for given input file
 *
 *  If the reader detected an error in the input file, it should return with SCIP_READERR or SCIP_NOFILE.
 */
#define DECL_READERREAD(x) RETCODE x (READER* reader, SCIP* scip, const char* filename, RESULT* result)



#include "scip.h"
#include "retcode.h"
#include "result.h"


extern
RETCODE SCIPreaderCreate(               /**< creates a reader */
   READER**         reader,             /**< pointer to store reader */
   const char*      name,               /**< name of reader */
   const char*      desc,               /**< description of reader */
   const char*      extension,          /**< file extension that reader processes */
   DECL_READERFREE((*readerfree)),      /**< destructor of reader */
   DECL_READERREAD((*readerread)),      /**< read method */
   READERDATA*      readerdata          /**< reader data */
   );

extern
RETCODE SCIPreaderFree(                 /**< frees memory of reader */
   READER**         reader,             /**< pointer to reader data structure */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPreaderRead(                 /**< reads problem data from file with given reader or returns SCIP_DIDNOTRUN */
   READER*          reader,             /**< reader */
   SCIP*            scip,               /**< SCIP data structure */   
   const char*      filename,           /**< name of the input file */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

extern
const char* SCIPreaderGetName(          /**< gets name of reader */
   READER*          reader              /**< reader */
   );

extern
READERDATA* SCIPreaderGetData(          /**< gets user data of reader */
   READER*          reader              /**< reader */
   );

extern
void SCIPreaderSetData(                 /**< sets user data of reader; user has to free old data in advance! */
   READER*          reader,             /**< reader */
   READERDATA*      readerdata          /**< new reader user data */
   );


#endif

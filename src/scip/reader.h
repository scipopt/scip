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
#pragma ident "@(#) $Id: reader.h,v 1.12 2003/11/24 12:12:43 bzfpfend Exp $"

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
 *  - scip            : SCIP main data structure
 *  - reader          : the reader itself
 */
#define DECL_READERFREE(x) RETCODE x (SCIP* scip, READER* reader)

/** problem reading method of reader
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - reader          : the reader itself
 *  - filename        : full path and name of file to read, or NULL if stdin should be used
 *  - result          : pointer to store the result of the propagation call
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : the reader read the file correctly
 *  - SCIP_DIDNOTRUN  : the reader is not responsible for given input file
 *
 *  If the reader detected an error in the input file, it should return with SCIP_READERR or SCIP_NOFILE.
 */
#define DECL_READERREAD(x) RETCODE x (SCIP* scip, READER* reader, const char* filename, RESULT* result)



#include "scip.h"
#include "retcode.h"
#include "result.h"
#include "set.h"


/** creates a reader */
extern
RETCODE SCIPreaderCreate(
   READER**         reader,             /**< pointer to store reader */
   const char*      name,               /**< name of reader */
   const char*      desc,               /**< description of reader */
   const char*      extension,          /**< file extension that reader processes */
   DECL_READERFREE  ((*readerfree)),    /**< destructor of reader */
   DECL_READERREAD  ((*readerread)),    /**< read method */
   READERDATA*      readerdata          /**< reader data */
   );

/** frees memory of reader */
extern
RETCODE SCIPreaderFree(
   READER**         reader,             /**< pointer to reader data structure */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** reads problem data from file with given reader or returns SCIP_DIDNOTRUN */
extern
RETCODE SCIPreaderRead(
   READER*          reader,             /**< reader */
   const SET*       set,                /**< global SCIP settings */
   const char*      filename,           /**< name of the input file */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** gets user data of reader */
extern
READERDATA* SCIPreaderGetData(
   READER*          reader              /**< reader */
   );

/** sets user data of reader; user has to free old data in advance! */
extern
void SCIPreaderSetData(
   READER*          reader,             /**< reader */
   READERDATA*      readerdata          /**< new reader user data */
   );

/** gets name of reader */
extern
const char* SCIPreaderGetName(
   READER*          reader              /**< reader */
   );

/** gets description of reader */
extern
const char* SCIPreaderGetDesc(
   READER*          reader              /**< reader */
   );

/** gets file extension of reader */
extern
const char* SCIPreaderGetExtension(
   READER*          reader              /**< reader */
   );


#endif

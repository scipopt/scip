/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
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
 *  possible return values:
 *    SCIP_OKAY       : normal termination
 *    neg. values     : error codes
 */
#define DECL_READERFREE(x) RETCODE x (READER* reader, SCIP* scip)

/** initialization method of reader (called at problem creation)
 *  possible return values:
 *    SCIP_OKAY   : normal termination
 *    neg. values : error codes
 */
#define DECL_READERINIT(x) RETCODE x (READER* reader, SCIP* scip)

/** deinitialization method of reader (called at problem destruction)
 *  possible return values:
 *    SCIP_OKAY   : normal termination
 *    neg. values : error codes
 */
#define DECL_READEREXIT(x) RETCODE x (READER* reader, SCIP* scip)

/** problem reading method of reader from file, if filename is NULL, read from stdin
 *  possible return values:
 *    SCIP_OKAY   : normal termination
 *    neg. values : error codes
 */
#define DECL_READERREAD(x) RETCODE x (READER* reader, SCIP* scip, const char* filename)



#include "scip.h"
#include "retcode.h"


extern
RETCODE SCIPreaderCreate(               /**< creates a reader */
   READER**         reader,             /**< pointer to store reader */
   const char*      name,               /**< name of reader */
   const char*      desc,               /**< description of reader */
   const char*      extension,          /**< file extension that reader processes */
   DECL_READERFREE((*readerfree)),      /**< destructor of reader */
   DECL_READERINIT((*readerinit)),      /**< initialise reader */
   DECL_READEREXIT((*readerexit)),      /**< deinitialise reader */
   DECL_READERREAD((*readerread)),      /**< read method */
   READERDATA*      readerdata          /**< reader data */
   );

extern
RETCODE SCIPreaderFree(                 /**< frees memory of reader */
   READER**         reader,             /**< pointer to reader data structure */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPreaderInit(                 /**< initializes reader */
   READER*          reader,             /**< reader */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPreaderExit(                 /**< deinitializes reader */
   READER*          reader,             /**< reader */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPreaderRead(                 /**< reads problem data from file with given reader or returns SCIP_DIDNOTRUN */
   READER*          reader,             /**< reader */
   SCIP*            scip,               /**< SCIP data structure */   
   const char*      filename            /**< name of the input file */
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

extern
Bool SCIPreaderIsInitialized(           /**< is reader initialized? */
   READER*          reader              /**< reader */
   );


#endif

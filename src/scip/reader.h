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
#pragma ident "@(#) $Id: reader.h,v 1.21 2005/02/14 13:35:49 bzfpfend Exp $"

/**@file   reader.h
 * @brief  internal methods for input file readers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __READER_H__
#define __READER_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_reader.h"
#include "scip/pub_reader.h"



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
   SET*             set                 /**< global SCIP settings */
   );

/** reads problem data from file with given reader or returns SCIP_DIDNOTRUN */
extern
RETCODE SCIPreaderRead(
   READER*          reader,             /**< reader */
   SET*             set,                /**< global SCIP settings */
   const char*      filename,           /**< name of the input file */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );


#endif

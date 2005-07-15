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
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_reader.h,v 1.7 2005/07/15 17:20:15 bzfpfend Exp $"

/**@file   pub_reader.h
 * @brief  public methods for input file readers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_READER_H__
#define __SCIP_PUB_READER_H__


#include "scip/def.h"
#include "scip/type_reader.h"



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

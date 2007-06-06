/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2007 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_reader.h,v 1.11 2007/06/06 11:25:23 bzfpfend Exp $"

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
SCIP_READERDATA* SCIPreaderGetData(
   SCIP_READER*          reader              /**< reader */
   );

/** sets user data of reader; user has to free old data in advance! */
extern
void SCIPreaderSetData(
   SCIP_READER*          reader,             /**< reader */
   SCIP_READERDATA*      readerdata          /**< new reader user data */
   );

/** gets name of reader */
extern
const char* SCIPreaderGetName(
   SCIP_READER*          reader              /**< reader */
   );

/** gets description of reader */
extern
const char* SCIPreaderGetDesc(
   SCIP_READER*          reader              /**< reader */
   );

/** gets file extension of reader */
extern
const char* SCIPreaderGetExtension(
   SCIP_READER*          reader              /**< reader */
   );


#endif

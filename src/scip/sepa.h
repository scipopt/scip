/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: sepa.h,v 1.28 2004/09/23 15:46:33 bzfpfend Exp $"

/**@file   sepa.h
 * @brief  internal methods for separators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SEPA_H__
#define __SEPA_H__


#include "def.h"
#include "memory.h"
#include "type_retcode.h"
#include "type_result.h"
#include "type_set.h"
#include "type_stat.h"
#include "type_sepastore.h"
#include "type_scip.h"
#include "type_sepa.h"
#include "pub_sepa.h"



/** creates a separator */
extern
RETCODE SCIPsepaCreate(
   SEPA**           sepa,               /**< pointer to separator data structure */
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory for parameter settings */
   const char*      name,               /**< name of separator */
   const char*      desc,               /**< description of separator */
   int              priority,           /**< priority of separator */
   int              freq,               /**< frequency for calling separator */
   DECL_SEPAFREE    ((*sepafree)),      /**< destructor of separator */
   DECL_SEPAINIT    ((*sepainit)),      /**< initialize separator */
   DECL_SEPAEXIT    ((*sepaexit)),      /**< deinitialize separator */
   DECL_SEPAEXEC    ((*sepaexec)),      /**< execution method of separator */
   SEPADATA*        sepadata            /**< separator data */
   );

/** calls destructor and frees memory of separator */
extern
RETCODE SCIPsepaFree(
   SEPA**           sepa,               /**< pointer to separator data structure */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** initializes separator */
extern
RETCODE SCIPsepaInit(
   SEPA*            sepa,               /**< separator */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** calls exit method of separator */
extern
RETCODE SCIPsepaExit(
   SEPA*            sepa,               /**< separator */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** calls execution method of separator */
extern
RETCODE SCIPsepaExec(
   SEPA*            sepa,               /**< separator */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   SEPASTORE*       sepastore,          /**< separation storage */
   int              depth,              /**< depth of current node */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of separator */
extern
void SCIPsepaSetPriority(
   SEPA*            sepa,               /**< separator */
   SET*             set,                /**< global SCIP settings */
   int              priority            /**< new priority of the separator */
   );

#endif

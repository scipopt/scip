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
#pragma ident "@(#) $Id: sepa.h,v 1.39 2005/07/15 17:20:17 bzfpfend Exp $"

/**@file   sepa.h
 * @brief  internal methods for separators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_H__
#define __SCIP_SEPA_H__


#include "scip/def.h"
#include "scip/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_sepastore.h"
#include "scip/type_sepa.h"
#include "scip/pub_sepa.h"



/** creates a separator */
extern
RETCODE SCIPsepaCreate(
   SEPA**           sepa,               /**< pointer to separator data structure */
   SET*             set,                /**< global SCIP settings */
   BLKMEM*          blkmem,             /**< block memory for parameter settings */
   const char*      name,               /**< name of separator */
   const char*      desc,               /**< description of separator */
   int              priority,           /**< priority of separator (>= 0: before, < 0: after constraint handlers) */
   int              freq,               /**< frequency for calling separator */
   Bool             delay,              /**< should separator be delayed, if other separators found cuts? */
   DECL_SEPAFREE    ((*sepafree)),      /**< destructor of separator */
   DECL_SEPAINIT    ((*sepainit)),      /**< initialize separator */
   DECL_SEPAEXIT    ((*sepaexit)),      /**< deinitialize separator */
   DECL_SEPAINITSOL ((*sepainitsol)),   /**< solving process initialization method of separator */
   DECL_SEPAEXITSOL ((*sepaexitsol)),   /**< solving process deinitialization method of separator */
   DECL_SEPAEXEC    ((*sepaexec)),      /**< execution method of separator */
   SEPADATA*        sepadata            /**< separator data */
   );

/** calls destructor and frees memory of separator */
extern
RETCODE SCIPsepaFree(
   SEPA**           sepa,               /**< pointer to separator data structure */
   SET*             set                 /**< global SCIP settings */
   );

/** initializes separator */
extern
RETCODE SCIPsepaInit(
   SEPA*            sepa,               /**< separator */
   SET*             set                 /**< global SCIP settings */
   );

/** calls exit method of separator */
extern
RETCODE SCIPsepaExit(
   SEPA*            sepa,               /**< separator */
   SET*             set                 /**< global SCIP settings */
   );

/** informs separator that the branch and bound process is being started */
extern
RETCODE SCIPsepaInitsol(
   SEPA*            sepa,               /**< separator */
   SET*             set                 /**< global SCIP settings */
   );

/** informs separator that the branch and bound process data is being freed */
extern
RETCODE SCIPsepaExitsol(
   SEPA*            sepa,               /**< separator */
   SET*             set                 /**< global SCIP settings */
   );

/** calls execution method of separator */
extern
RETCODE SCIPsepaExec(
   SEPA*            sepa,               /**< separator */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   SEPASTORE*       sepastore,          /**< separation storage */
   int              depth,              /**< depth of current node */
   Bool             execdelayed,        /**< execute separator even if it is marked to be delayed */
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

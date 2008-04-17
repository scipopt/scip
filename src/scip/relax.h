/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2008 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: relax.h,v 1.16 2008/04/17 17:49:16 bzfpfets Exp $"

/**@file   relax.h
 * @brief  internal methods for relaxators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_RELAX_H__
#define __SCIP_RELAX_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_relax.h"
#include "scip/pub_relax.h"



/** creates a relaxator */
extern
SCIP_RETCODE SCIPrelaxCreate(
   SCIP_RELAX**          relax,              /**< pointer to relaxator data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of relaxator */
   const char*           desc,               /**< description of relaxator */
   int                   priority,           /**< priority of the relaxator (negative: after LP, non-negative: before LP) */
   int                   freq,               /**< frequency for calling relaxator */
   SCIP_DECL_RELAXFREE   ((*relaxfree)),     /**< destructor of relaxator */
   SCIP_DECL_RELAXINIT   ((*relaxinit)),     /**< initialize relaxator */
   SCIP_DECL_RELAXEXIT   ((*relaxexit)),     /**< deinitialize relaxator */
   SCIP_DECL_RELAXINITSOL((*relaxinitsol)),  /**< solving process initialization method of relaxator */
   SCIP_DECL_RELAXEXITSOL((*relaxexitsol)),  /**< solving process deinitialization method of relaxator */
   SCIP_DECL_RELAXEXEC   ((*relaxexec)),     /**< execution method of relaxator */
   SCIP_RELAXDATA*       relaxdata           /**< relaxator data */
   );

/** calls destructor and frees memory of relaxator */
extern
SCIP_RETCODE SCIPrelaxFree(
   SCIP_RELAX**          relax,              /**< pointer to relaxator data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** initializes relaxator */
extern
SCIP_RETCODE SCIPrelaxInit(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls exit method of relaxator */
extern
SCIP_RETCODE SCIPrelaxExit(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs relaxator that the branch and bound process is being started */
extern
SCIP_RETCODE SCIPrelaxInitsol(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs relaxator that the branch and bound process data is being freed */
extern
SCIP_RETCODE SCIPrelaxExitsol(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls execution method of relaxator */
extern
SCIP_RETCODE SCIPrelaxExec(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   int                   depth,              /**< depth of current node */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of relaxator */
extern
void SCIPrelaxSetPriority(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the relaxator */
   );

/** returns whether the relaxation was completely solved at the current node */
extern
SCIP_Bool SCIPrelaxIsSolved(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   );

#endif

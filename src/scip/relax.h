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
#pragma ident "@(#) $Id: relax.h,v 1.8 2005/02/14 13:35:49 bzfpfend Exp $"

/**@file   relax.h
 * @brief  internal methods for relaxators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __RELAX_H__
#define __RELAX_H__


#include "scip/def.h"
#include "scip/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_relax.h"
#include "scip/pub_relax.h"



/** creates a relaxator */
extern
RETCODE SCIPrelaxCreate(
   RELAX**          relax,              /**< pointer to relaxator data structure */
   SET*             set,                /**< global SCIP settings */
   BLKMEM*          blkmem,             /**< block memory for parameter settings */
   const char*      name,               /**< name of relaxator */
   const char*      desc,               /**< description of relaxator */
   int              priority,           /**< priority of the relaxator (negative: after LP, non-negative: before LP) */
   int              freq,               /**< frequency for calling relaxator */
   DECL_RELAXFREE   ((*relaxfree)),     /**< destructor of relaxator */
   DECL_RELAXINIT   ((*relaxinit)),     /**< initialize relaxator */
   DECL_RELAXEXIT   ((*relaxexit)),     /**< deinitialize relaxator */
   DECL_RELAXINITSOL((*relaxinitsol)),  /**< solving process initialization method of relaxator */
   DECL_RELAXEXITSOL((*relaxexitsol)),  /**< solving process deinitialization method of relaxator */
   DECL_RELAXEXEC   ((*relaxexec)),     /**< execution method of relaxator */
   RELAXDATA*       relaxdata           /**< relaxator data */
   );

/** calls destructor and frees memory of relaxator */
extern
RETCODE SCIPrelaxFree(
   RELAX**          relax,              /**< pointer to relaxator data structure */
   SET*             set                 /**< global SCIP settings */
   );

/** initializes relaxator */
extern
RETCODE SCIPrelaxInit(
   RELAX*           relax,              /**< relaxator */
   SET*             set                 /**< global SCIP settings */
   );

/** calls exit method of relaxator */
extern
RETCODE SCIPrelaxExit(
   RELAX*           relax,              /**< relaxator */
   SET*             set                 /**< global SCIP settings */
   );

/** informs relaxator that the branch and bound process is being started */
extern
RETCODE SCIPrelaxInitsol(
   RELAX*           relax,              /**< relaxator */
   SET*             set                 /**< global SCIP settings */
   );

/** informs relaxator that the branch and bound process data is being freed */
extern
RETCODE SCIPrelaxExitsol(
   RELAX*           relax,              /**< relaxator */
   SET*             set                 /**< global SCIP settings */
   );

/** calls execution method of relaxator */
extern
RETCODE SCIPrelaxExec(
   RELAX*           relax,              /**< relaxator */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   int              depth,              /**< depth of current node */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of relaxator */
extern
void SCIPrelaxSetPriority(
   RELAX*           relax,              /**< relaxator */
   SET*             set,                /**< global SCIP settings */
   int              priority            /**< new priority of the relaxator */
   );

/** returns whether the relaxation was completely solved at the current node */
extern
Bool SCIPrelaxIsSolved(
   RELAX*           relax,              /**< relaxator */
   STAT*            stat                /**< dynamic problem statistics */
   );

#endif

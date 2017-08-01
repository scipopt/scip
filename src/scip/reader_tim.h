/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_tim.h
 * @ingroup FILEREADERS
 * @brief  (extended) TIM file reader
 * @author Thorsten Koch
 * @author Tobias Achterberg
 *
 * This reader allows to parse and write TIM files with linear and quadratic constraints and objective,
 * special ordered sets of type 1 and 2, indicators on linear constraints, and semicontinuous variables.
 * For writing, linear (general and specialized), indicator, quadratic, second order cone, and
 * special ordered set constraints are supported.
 *
 * See http://en.wikipedia.org/wiki/TIM_%28format%29 for a description.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_TIM_H__
#define __SCIP_READER_TIM_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the tim file reader into SCIP
 *
 *  @ingroup FileReaderIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeReaderTim(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup FILEREADERS
 *
 * @{
 */

/** reads problem from file */
EXTERN
SCIP_RETCODE SCIPreadTim(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader,             /**< the file reader itself */
   const char*           filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*          result              /**< pointer to store the result of the file reading call */
   );

/** writes problem to file */
EXTERN
SCIP_RETCODE SCIPwriteTim(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL if standard output should be used */
   const char*           name,               /**< problem name */
   SCIP_Bool             transformed,        /**< TRUE iff problem is the transformed problem */
   SCIP_OBJSENSE         objsense,           /**< objective sense */
   SCIP_Real             objscale,           /**< scalar applied to objective function; external objective value is
                                              * extobj = objsense * objscale * (intobj + objoffset) */
   SCIP_Real             objoffset,          /**< objective offset from bound shifting and fixing */
   SCIP_VAR**            vars,               /**< array with active variables ordered binary, integer, implicit, continuous */
   int                   nvars,              /**< number of active variables in the problem */
   int                   nbinvars,           /**< number of binary variables */
   int                   nintvars,           /**< number of general integer variables */
   int                   nimplvars,          /**< number of implicit integer variables */
   int                   ncontvars,          /**< number of continuous variables */
   SCIP_CONS**           conss,              /**< array with constraints of the problem */
   int                   nconss,             /**< number of constraints in the problem */
   SCIP_RESULT*          result              /**< pointer to store the result of the file writing call */
   );

/* @} */

/*
 * Interface methods for the cor and sto files
 */

/* returns the number of stages */
EXTERN
int SCIPtimGetNStage(
   SCIP*                 scip                /**< SCIP data structure */
   );

/* returns the name for a given stage */
EXTERN
const char* SCIPtimGetStageName(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   stagenum            /**< the number of the requested stage */
   );

/* returns the number for a given stage */
EXTERN
int SCIPtimFindStage(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           stage               /**< the name of the requested stage */
   );

/* returns the array of variables for a given stage */
EXTERN
SCIP_VAR** SCIPtimGetStageVars(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   stagenum            /**< the number of the requested stage */
   );

/* returns an array of constraints for a given stage */
EXTERN
SCIP_CONS** SCIPtimGetStageConss(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   stagenum            /**< the number of the requested stage */
   );

/* returns the number of variables for a given stage */
EXTERN
int SCIPtimGetStageNVars(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   stagenum            /**< the number of the requested stage */
   );

/* returns the number of constraints for a given stage */
EXTERN
int SCIPtimGetStageNConss(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   stagenum            /**< the number of the requested stage */
   );

/* returns the number of scenarios for a given stage */
EXTERN
int SCIPtimGetStageNScenarios(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   stagenum            /**< the number of the requested stage */
   );

/* returns the number of scenarios for a given stage */
EXTERN
void SCIPtimSetStageNScenarios(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   stagenum,           /**< the number of the requested stage */
   int                   nscenarios          /**< the number of scenarios to set the stage parameter to */
   );

#ifdef __cplusplus
}
#endif

#endif

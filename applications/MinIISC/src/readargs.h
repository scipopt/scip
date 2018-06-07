/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   readargs.h
 * @brief  read comand line arguments
 * @author Marc Pfetsch
 */

#ifndef __READARGS_H__
#define __READARGS_H__

#include <scip/scip.h>
#include <string.h>

/** get problem name
 *
 *  Returns 0 if maxsize is not sufficient and 1 otherwise.
 */
extern
int getProblemName(
   const char*           filename,           /**< file name */
   char*                 probname,           /**< name of problem (output) */
   int                   maxsize             /**< max length of problem name */
   );

/** read comand line arguments */
extern
SCIP_RETCODE readArguments(
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array with shell parameters */
   const char**          filename,           /**< file name from arguments */
   const char**          settingsname,       /**< name of settings file */
   SCIP_Real*            timelimit,          /**< time limit read from arguments */
   SCIP_Real*            memlimit,           /**< memory limit read from arguments */
   SCIP_Longint*         nodelimit,          /**< node limit read from arguments */
   int*                  dispfreq            /**< display frequency */
   );

#endif

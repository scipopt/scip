/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef __SCIP_PB_MESSAGE_PB_H__
#define __SCIP_PB_MESSAGE_PB_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** message handler data */
struct SCIP_MessagehdlrData
{
   SCIP_Bool             comment;            /**< should output be printed as PB comments */
};

/** creates default PBSCIP message handler */
extern
SCIP_RETCODE SCIPcreateMessagehdlrPbscip(
   SCIP_MESSAGEHDLR**    messagehdlr,        /**< pointer to message handler */
   SCIP_Bool             buffered,           /**< should the output be buffered */
   const char*           filename,           /**< name of log file, or NULL (stdout) */
   SCIP_Bool             quiet               /**< should screen messages be suppressed? */
   );

#ifdef __cplusplus
}
#endif

#endif

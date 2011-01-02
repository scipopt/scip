/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: reader_pip.h,v 1.3 2011/01/02 11:10:43 bzfheinz Exp $"

/**@file   reader_pip.h
 * @ingroup FILEREADERS 
 * @brief  file reader for polynomial mixed-integer programs in PIP format
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_PIP_H__
#define __SCIP_READER_PIP_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the pip file reader into SCIP */
extern
SCIP_RETCODE SCIPincludeReaderPip(
   SCIP*                 scip                /**< SCIP data structure */
   );

/* reads problem from file */
extern
SCIP_RETCODE SCIPreadPip(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_READER*       reader,             /**< the file reader itself */
   const char*        filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*       result              /**< pointer to store the result of the file reading call */
   );

#ifdef __cplusplus
}
#endif

#endif

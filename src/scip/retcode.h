/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   retcode.h
 * @brief  return codes for SCIP methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __RETCODE_H__
#define __RETCODE_H__

#include <stdio.h>

/** return codes for SCIP methods: negative return codes are errors */
enum Retcode
{
   SCIP_OKAY          =   0,            /**< normal termination */
   SCIP_ERROR         =  -1,            /**< unspecified error */
   SCIP_NOMEMORY      =  -2,            /**< insufficient memory error */
   SCIP_READERR       =  -3,            /**< file read error */
   SCIP_NOFILE        =  -4,            /**< file not found error */
   SCIP_LPERROR       =  -5,            /**< error in LP solver */
   SCIP_NOPROBLEM     =  -6,            /**< no problem exists */
   SCIP_INVALIDCALL   =  -7,            /**< method cannot be called at this time in solution process */
   SCIP_INVALIDDATA   =  -8,            /**< error in input data */
   SCIP_INVALIDRESULT =  -9             /**< method returned an invalid result code */
};
typedef enum Retcode RETCODE;           /**< return code for SCIP method */



extern
void SCIPretcodePrint(                  /**< prints error message for return code */
   FILE*            errout,             /**< file stream to write error message */
   RETCODE          retcode             /**< SCIP return code causing the error */
   );

#endif

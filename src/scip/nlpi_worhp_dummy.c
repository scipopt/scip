/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    nlpi_worhp_dummy.c
 * @ingroup DEFPLUGINS_NLPI
 * @brief   dummy WORHP NLP interface
 * @author  Benjamin Mueller
 *
 * This file contains dummy implementations of the interface methods for the Worhp interface.
 * It is used when SCIP is build without Worhp.
 */

#include "scip/nlpi_worhp.h"

/** create solver interface for Worhp solver and includes it into SCIP, if Worhp is available */  /*lint -e{715}*/
SCIP_RETCODE SCIPincludeNlpSolverWorhp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             useip               /**< TRUE for using Interior Point, FALSE for SQP */
   )
{  /*lint --e{715}*/
   assert(scip != NULL);

   return SCIP_OKAY;
}

/** gets string that identifies Worhp (version number) */
const char* SCIPgetSolverNameWorhp(void)
{
   return "WORHP";
}

/** gets string that describes Worhp (version number) */
const char* SCIPgetSolverDescWorhp(void)
{
   return "this is WORHP";
}

/** returns whether Worhp is available, i.e., whether it has been linked in */
SCIP_Bool SCIPisWorhpAvailableWorhp(void)
{
   return FALSE;
}

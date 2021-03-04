/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    nlpi_worhp_dummy.c
 * @ingroup NLPIS
 * @brief   dummy WORHP NLP interface
 * @author  Benjamin Mueller
 */

#include "nlpi/nlpi_worhp.h"

/** create solver interface for Worhp solver and includes it into SCIP, if Worhp is available */
SCIP_RETCODE SCIPincludeNlpSolverWorhp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             useip               /**< TRUE for using Interior Point, FALSE for SQP */
   )
{
   return SCIP_OKAY;
} /*lint !e715*/

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

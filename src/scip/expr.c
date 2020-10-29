/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   expr.c
 * @brief  functions for algebraic expressions
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 */

#include "scip/scip_expr.h"
#include "scip/expr.h"
#include "nlpi/nlpi_ipopt.h" /* for LAPACK */

/** variable mapping data passed on during copying expressions when copying SCIP instances */
typedef struct
{
   SCIP_HASHMAP*         varmap;             /**< SCIP_HASHMAP mapping variables of the source SCIP to corresponding variables of the target SCIP */
   SCIP_HASHMAP*         consmap;            /**< SCIP_HASHMAP mapping constraints of the source SCIP to corresponding constraints of the target SCIP */
   SCIP_Bool             global;             /**< should a global or a local copy be created */
   SCIP_Bool             valid;              /**< indicates whether every variable copy was valid */
} COPY_MAPVAR_DATA;

/** printing to file data */
struct SCIP_ExprPrintData
{
   FILE*                 file;               /**< file to print to */
   SCIP_EXPRITER*        iterator;           /**< iterator to use */
   SCIP_Bool             closefile;          /**< whether file need to be closed when finished printing */
   SCIP_HASHMAP*         leaveexprs;         /**< hashmap storing leave (no children) expressions */
   SCIP_EXPRPRINT_WHAT   whattoprint;        /**< flags that indicate what to print for each expression */
};

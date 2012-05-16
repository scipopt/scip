/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_paramset.h
 * @ingroup PUBLICMETHODS
 * @brief  public methods for handling parameter settings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_PARAMSET_H__
#define __SCIP_PUB_PARAMSET_H__


#include "scip/def.h"
#include "scip/type_paramset.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** returns type of parameter */
extern
SCIP_PARAMTYPE SCIPparamGetType(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns name of parameter */
extern
const char* SCIPparamGetName(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns description of parameter */
extern
const char* SCIPparamGetDesc(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns locally defined parameter specific data */
extern
SCIP_PARAMDATA* SCIPparamGetData(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns locally defined parameter specific data */
extern
SCIP_Bool SCIPparamIsAdvanced(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns value of SCIP_Bool parameter */
extern
SCIP_Bool SCIPparamGetBool(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns default value of SCIP_Bool parameter */
extern
SCIP_Bool SCIPparamGetBoolDefault(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns value of int parameter */
extern
int SCIPparamGetInt(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns minimal value of int parameter */
extern
int SCIPparamGetIntMin(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns maximal value of int parameter */
extern
int SCIPparamGetIntMax(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns default value of int parameter */
extern
int SCIPparamGetIntDefault(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns value of SCIP_Longint parameter */
extern
SCIP_Longint SCIPparamGetLongint(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns minimal value of longint parameter */
extern
SCIP_Longint SCIPparamGetLongintMin(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns maximal value of longint parameter */
extern
SCIP_Longint SCIPparamGetLongintMax(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns default value of SCIP_Longint parameter */
extern
SCIP_Longint SCIPparamGetLongintDefault(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns value of SCIP_Real parameter */
extern
SCIP_Real SCIPparamGetReal(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns minimal value of real parameter */
extern
SCIP_Real SCIPparamGetRealMin(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns maximal value of real parameter */
extern
SCIP_Real SCIPparamGetRealMax(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns default value of SCIP_Real parameter */
extern
SCIP_Real SCIPparamGetRealDefault(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns value of char parameter */
extern
char SCIPparamGetChar(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns default value of char parameter */
extern
char SCIPparamGetCharDefault(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns value of string parameter */
extern
char* SCIPparamGetString(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns default value of String parameter */
extern
char* SCIPparamGetStringDefault(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns whether the parameter is on its default setting */
extern
SCIP_Bool SCIPparamIsDefault(
   SCIP_PARAM*           param               /**< parameter */
   );

#ifdef __cplusplus
}
#endif

#endif

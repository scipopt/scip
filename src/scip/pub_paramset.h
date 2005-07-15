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
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_paramset.h,v 1.8 2005/07/15 17:20:15 bzfpfend Exp $"

/**@file   pub_paramset.h
 * @brief  public methods for handling parameter settings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_PARAMSET_H__
#define __SCIP_PUB_PARAMSET_H__


#include "scip/def.h"
#include "scip/type_paramset.h"
#include "scip/type_scip.h"



/** returns type of parameter */
extern
PARAMTYPE SCIPparamGetType(
   PARAM*           param               /**< parameter */
   );

/** returns name of parameter */
extern
const char* SCIPparamGetName(
   PARAM*           param               /**< parameter */
   );

/** returns description of parameter */
extern
const char* SCIPparamGetDesc(
   PARAM*           param               /**< parameter */
   );

/** returns locally defined parameter specific data */
extern
PARAMDATA* SCIPparamGetData(
   PARAM*           param               /**< parameter */
   );

/** returns value of Bool parameter */
extern
Bool SCIPparamGetBool(
   PARAM*           param               /**< parameter */
   );

/** returns default value of Bool parameter */
extern
Bool SCIPparamGetBoolDefault(
   PARAM*           param               /**< parameter */
   );

/** returns value of int parameter */
extern
int SCIPparamGetInt(
   PARAM*           param               /**< parameter */
   );

/** returns minimal value of int parameter */
extern
int SCIPparamGetIntMin(
   PARAM*           param               /**< parameter */
   );

/** returns maximal value of int parameter */
extern
int SCIPparamGetIntMax(
   PARAM*           param               /**< parameter */
   );

/** returns default value of int parameter */
extern
int SCIPparamGetIntDefault(
   PARAM*           param               /**< parameter */
   );

/** returns value of Longint parameter */
extern
Longint SCIPparamGetLongint(
   PARAM*           param               /**< parameter */
   );

/** returns minimal value of longint parameter */
extern
Longint SCIPparamGetLongintMin(
   PARAM*           param               /**< parameter */
   );

/** returns maximal value of longint parameter */
extern
Longint SCIPparamGetLongintMax(
   PARAM*           param               /**< parameter */
   );

/** returns default value of Longint parameter */
extern
Longint SCIPparamGetLongintDefault(
   PARAM*           param               /**< parameter */
   );

/** returns value of Real parameter */
extern
Real SCIPparamGetReal(
   PARAM*           param               /**< parameter */
   );

/** returns minimal value of real parameter */
extern
Real SCIPparamGetRealMin(
   PARAM*           param               /**< parameter */
   );

/** returns maximal value of real parameter */
extern
Real SCIPparamGetRealMax(
   PARAM*           param               /**< parameter */
   );

/** returns default value of Real parameter */
extern
Real SCIPparamGetRealDefault(
   PARAM*           param               /**< parameter */
   );

/** returns value of char parameter */
extern
char SCIPparamGetChar(
   PARAM*           param               /**< parameter */
   );

/** returns default value of char parameter */
extern
char SCIPparamGetCharDefault(
   PARAM*           param               /**< parameter */
   );

/** returns value of string parameter */
extern
char* SCIPparamGetString(
   PARAM*           param               /**< parameter */
   );

/** returns default value of String parameter */
extern
char* SCIPparamGetStringDefault(
   PARAM*           param               /**< parameter */
   );

/** sets value of Bool parameter */
RETCODE SCIPparamSetBool(
   PARAM*           param,              /**< parameter */
   SCIP*            scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   Bool             value               /**< new value of the parameter */
   );

/** sets value of int parameter */
extern
RETCODE SCIPparamSetInt(
   PARAM*           param,              /**< parameter */
   SCIP*            scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   int              value               /**< new value of the parameter */
   );

/** sets value of Longint parameter */
extern
RETCODE SCIPparamSetLongint(
   PARAM*           param,              /**< parameter */
   SCIP*            scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   Longint          value               /**< new value of the parameter */
   );

/** sets value of Real parameter */
extern
RETCODE SCIPparamSetReal(
   PARAM*           param,              /**< parameter */
   SCIP*            scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   Real             value               /**< new value of the parameter */
   );

/** sets value of char parameter */
extern
RETCODE SCIPparamSetChar(
   PARAM*           param,              /**< parameter */
   SCIP*            scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   char             value               /**< new value of the parameter */
   );

/** sets value of string parameter */
extern
RETCODE SCIPparamSetString(
   PARAM*           param,              /**< parameter */
   SCIP*            scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   const char*      value               /**< new value of the parameter */
   );

/** returns whether the parameter is on its default setting */
extern
Bool SCIPparamIsDefault(
   PARAM*           param               /**< parameter */
   );

/** sets the parameter to its default setting */
extern
RETCODE SCIPparamSetToDefault(
   PARAM*           param,              /**< parameter */
   SCIP*            scip                /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   );


#endif

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
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: intervalarith.h,v 1.2 2005/01/18 09:26:47 bzfpfend Exp $"

/**@file   intervalarith.h
 * @brief  interval arithmetics for provable bounds
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __INTERVALARITH_H__
#define __INTERVALARITH_H__


#include "def.h"



/** interval given by infimum and supremum */
struct Interval
{
   Real             inf;                /**< infimum (lower bound) of interval */
   Real             sup;                /**< supremum (upper bound) of interval */
};
typedef struct Interval INTERVAL;



/*
 * Interval arithmetic operations
 */

/** stores given value as interval */
extern
void SCIPintervalSet(
   INTERVAL*        resultant,          /**< interval to store value into */
   Real             value               /**< value to store */
   );

/** stores given infimum and supremum as interval */
extern
void SCIPintervalSetBounds(
   INTERVAL*        resultant,          /**< interval to store value into */
   Real             inf,                /**< value to store as infimum */
   Real             sup                 /**< value to store as supremum */
   );

/** adds operand1 and operand2 and stores result in resultant */
extern
void SCIPintervalAdd(
   INTERVAL*        resultant,          /**< resultant interval of operation */
   INTERVAL         operand1,           /**< first operand of operation */
   INTERVAL         operand2            /**< second operand of operation */
   );

/** substracts operand2 from operand1 and stores result in resultant */
extern
void SCIPintervalSub(
   INTERVAL*        resultant,          /**< resultant interval of operation */
   INTERVAL         operand1,           /**< first operand of operation */
   INTERVAL         operand2            /**< second operand of operation */
   );

/** multiplies operand1 with operand2 and stores result in resultant */
extern
void SCIPintervalMul(
   INTERVAL*        resultant,          /**< resultant interval of operation */
   INTERVAL         operand1,           /**< first operand of operation */
   INTERVAL         operand2            /**< second operand of operation */
   );

/** returns infimum of interval */
extern
Real SCIPintervalGetInf(
   INTERVAL         interval            /**< interval */
   );

/** returns supremum of interval */
extern
Real SCIPintervalGetSup(
   INTERVAL         interval            /**< interval */
   );


#endif

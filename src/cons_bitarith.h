/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_bitarith.h
 * @brief  constraint handler for bitarith constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CONS_BITARITH_H__
#define __CONS_BITARITH_H__


/** type of arithmetic bit operation */
#define SCIP_NBITARITHTYPES   4         /**< number of different arithmetic bit operations defined */
enum BitarithType
{
   SCIP_BITARITHTYPE_ADD    = 0,        /**< truncated addition:    z  =  (x + y) % len(z)     */
   SCIP_BITARITHTYPE_SHL    = 1,        /**< truncated shift left:  z  =  (x << y) % len(z)    */
   SCIP_BITARITHTYPE_EQ     = 2,        /**< equality operator:     z <-> (x == y)             */
   SCIP_BITARITHTYPE_NOT    = 3,        /**< bitwise not operator:  z  =  ~x                   */
};
typedef enum BitarithType BITARITHTYPE;


#include "scip.h"


/** creates the handler for bitarith constraints and includes it in SCIP */
extern
RETCODE SCIPincludeConsHdlrBitarith(
   SCIP*            scip                /**< SCIP data structure */
   );

/** creates and captures a bitarith constraint */
extern
RETCODE SCIPcreateConsBitarith(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   BITARITHTYPE     bitarithtype,       /**< type of arithmetic bit operation */
   CONS*            operand1,           /**< bitvar constraint: first (left) operand in operation (x) */
   CONS*            operand2,           /**< bitvar constraint: second (right) operand in operation (y) */
   CONS*            resultant,          /**< bitvar constraint: result of operation (z) */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   Bool             removeable          /**< should the constraint be removed from the LP due to aging or cleanup? */
   );

/** gets number of operands for given arithmetic operation */
extern
int SCIPgetArityBitarith(
   BITARITHTYPE     bitarithtype        /**< type of arithmetic bit operation */
   );

#endif

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

/**@file   cons_bitvar.h
 * @brief  constraint handler for bitvar constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CONS_BITVAR_H__
#define __CONS_BITVAR_H__


#include "scip.h"


/** creates the handler for bitvar constraints and includes it in SCIP */
extern
RETCODE SCIPincludeConsHdlrBitvar(
   SCIP*            scip                /**< SCIP data structure */
   );

/** creates and captures a bitvar constraint
 *  Warning! Either the bitvar should be short, or the objective value should be zero, because the objective
 *  value of the most significant bit in the variable would be 2^(nbits-1)*obj
 */
extern
RETCODE SCIPcreateConsBitvar(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nbits,              /**< number of bits in the bitvar */
   Real             obj,                /**< objective value of bitvar variable */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             removeable          /**< should the constraint be removed from the LP due to aging or cleanup? */
   );

/** creates and captures a constant bitvar constraint with constant given as bit vector
 *  Warning! Either the bitvar should be short, or the objective value should be zero, because the objective
 *  value of the most significant bit in the variable would be 2^(nbits-1)*obj
 */
extern
RETCODE SCIPcreateConsBitconst(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nbits,              /**< number of bits in the bitvar */
   Real             obj,                /**< objective value of bitvar constant */
   Bool*            fixedbits           /**< array with the fixed bit values of constant (least significant first) */
   );

/** creates and captures a constant bitvar constraint with constant parsed from a string
 *  Warning! Either the bitvar should be short, or the objective value should be zero, because the objective
 *  value of the most significant bit in the variable would be 2^(nbits-1)*obj
 */
extern
RETCODE SCIPcreateConsBitconstString(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nbits,              /**< number of bits in the bitvar */
   Real             obj,                /**< objective value of bitvar constant */
   const char*      cstring             /**< constant given as a string */
   );

/** gets number of bits in bitvar */
extern
int SCIPgetNBitsBitvar(
   CONS*            cons                /**< bitvar constraint */
   );

/** gets array with bits of bitvar, sorted least significant bit first */
extern
VAR** SCIPgetBitsBitvar(
   CONS*            cons                /**< bitvar constraint */
   );

/** gets variable for single bit in bitvar */
extern
VAR* SCIPgetBitBitvar(
   CONS*            cons,               /**< bitvar constraint */
   int              bit                 /**< bit number to get variable for */
   );

/** gets number of words in bitvar */
extern
int SCIPgetNWordsBitvar(
   CONS*            cons                /**< bitvar constraint */
   );

/** gets array with words of bitvar, sorted least significant word first */
extern
VAR** SCIPgetWordsBitvar(
   CONS*            cons                /**< bitvar constraint */
   );

/** gets variable for single word in bitvar */
extern
VAR* SCIPgetWordBitvar(
   CONS*            cons,               /**< bitvar constraint */
   int              word                /**< word number to get variable for */
   );

/** gets number of bits in a given word of a bitvar */
extern
int SCIPgetNWordBitsBitvar(
   CONS*            cons,               /**< bitvar constraint */
   int              word                /**< word number */
   );

/** returns the number of bits of the given word */
extern
int SCIPgetWordSizeBitvar(
   CONS*            cons,               /**< bitvar constraint */
   int              word                /**< word number */
   );

/** returns the number of different values the given word can store (2^#bits) */
extern
int SCIPgetWordPowerBitvar(
   CONS*            cons,               /**< bitvar constraint */
   int              word                /**< word number */
   );

/** returns whether the bitvar was created as a constant */
extern
Bool SCIPisConstBitvar(
   CONS*            cons                /**< bitvar constraint */
   );

#endif

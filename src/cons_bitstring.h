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

/**@file   cons_bitstring.h
 * @brief  constraint handler for bitstring constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CONS_BITSTRING_H__
#define __CONS_BITSTRING_H__


#include "scip.h"


/** creates the handler for bitstring constraints and includes it in SCIP */
extern
RETCODE SCIPincludeConsHdlrBitstring(
   SCIP*            scip                /**< SCIP data structure */
   );

/** creates and captures a bitstring constraint
 *  Warning! Either the bitstring should be short, or the objective value should be zero, because the objective
 *  value of the most significant bit in the string would be 2^(nbits-1)*obj
 */
extern
RETCODE SCIPcreateConsBitstring(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nbits,              /**< number of bits in the bitstring */
   Real             obj,                /**< objective value of bitstring variable */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             removeable          /**< should the constraint be removed from the LP due to aging or cleanup? */
   );

#endif

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_randnumgen.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for random numbers
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_RANDNUMGEN_H__
#define __SCIP_SCIP_RANDNUMGEN_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons).
 * Additionally, the internal "set.h" is included, such that the defines in set.h are
 * available in optimized mode.
 */
#ifdef NDEBUG
#include "scip/struct_scip.h"
#include "scip/struct_stat.h"
#include "scip/set.h"
#include "scip/tree.h"
#include "scip/misc.h"
#include "scip/var.h"
#include "scip/cons.h"
#include "scip/solve.h"
#include "scip/debug.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup RandomNumbers
 *
 * @{
 */

/** creates and initializes a random number generator
 *
 *  @note The initial seed is changed using SCIPinitializeRandomSeed()
 */
EXTERN
SCIP_RETCODE SCIPcreateRandom(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN**     randnumgen,         /**< random number generator */
   unsigned int          initialseed,        /**< initial random seed */
   SCIP_Bool             useglobalseed       /**< should SCIP's global seed be used to initialise the supplied seed? */
   );

/** frees a random number generator */
EXTERN
void SCIPfreeRandom(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN**     randnumgen          /**< random number generator */
   );

/** initializes a random number generator with a given seed
 *
 *  @note The seed is changed using SCIPinitializeRandomSeed()
 */
EXTERN
void SCIPsetRandomSeed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   unsigned int          seed                /**< new random seed */
   );


/** modifies an initial seed value with the global shift of random seeds */
EXTERN
unsigned int SCIPinitializeRandomSeed(
   SCIP*                 scip,               /**< SCIP data structure */
   unsigned int          initialseedvalue    /**< initial seed value to be modified */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif

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

/**@file   pricer.h
 * @brief  datastructures and methods for variable pricers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PRICER_H__
#define __PRICER_H__


typedef struct Pricer PRICER;           /**< variable pricer data */


/** destructor of variable pricer to free user data (called when SCIP is exiting)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    pricer          : the variable pricer itself
 */
#define DECL_PRICERFREE(x) RETCODE x (SCIP* scip, PRICER* pricer)

/** initialization method of variable pricer (called when problem solving starts)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    pricer          : the variable pricer itself
 */
#define DECL_PRICERINIT(x) RETCODE x (SCIP* scip, PRICER* pricer)

/** deinitialization method of variable pricer (called when problem solving exits)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    pricer          : the variable pricer itself
 */
#define DECL_PRICEREXIT(x) RETCODE x (SCIP* scip, PRICER* pricer)



#endif

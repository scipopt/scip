/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   domain.h
 * @brief  datastructures and methods for managing domains of variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __DOMAIN_H__
#define __DOMAIN_H__

#define SCIP_MAXCOL 0x08000             /**< maximal number of columns; 31 bits available to store column number */

enum BoundType
{
   SCIP_BOUNDTYPE_LOWER = 0,            /**< lower bound */
   SCIP_BOUNDTYPE_UPPER = 1             /**< upper bound */
};
typedef enum BoundType BOUNDTYPE;

typedef struct Domain DOMAIN;           /**< datastructures for storing domains of variables */
typedef struct Holelist HOLELIST;       /**< list of holes in a domain of an integer variable */
typedef struct DomainChg DOMAINCHG;     /**< changes in domains of variables */
typedef struct BoundChg BOUNDCHG;       /**< changes in bounds of variables */
typedef struct HoldChg HOLECHG;         /**< changes in holelist of variables */

#endif

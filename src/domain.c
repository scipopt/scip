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

/**@file   domain.c
 * @brief  datastructures and methods for managing domains of variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "domain.h"

struct Domain
{
   double*          lb;                 /**< lower bounds of variables */
   double*          ub;                 /**< upper bounds of variables */
   HOLELIST*        holelist;           /**< list of holes for each integer variable */
};

struct Holelist
{
   int              first;              /**< first value of hole */
   int              last;               /**< last value of hole */
   HOLELIST*        next;               /**< next hole in list */
};

struct DomainChg
{
   BOUNDCHG*        boundchg;           /**< changes in bounds of variables */
   HOLECHG*         holechg;            /**< changes in hole lists */
};

struct BoundChg
{
   unsigned int     boundtype:1;        /**< type of bound: lower or upper bound */
   unsigned int     col:31;             /**< column number */
   double           newbound;           /**< new value for bound */
   double           oldbound;           /**< old value for bound */
};

struct HoleChg
{
   HOLELIST**       ptr;                /**< changed list pointer */
   HOLELIST*        newlist;            /**< new value of list pointer */
   HOLELIST*        oldlist;            /**< old value of list pointer */
};

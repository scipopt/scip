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

/**@file   lp.h
 * @brief  LP management datastructures and methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __LP_H__
#define __LP_H__

#include "memory.h"
#include "lpi.h"


typedef struct Basis BASIS;             /**< compressed LP basis information */
typedef struct Row ROW;                 /**< row of an LP */
typedef struct RowList ROWLIST;         /**< list of LP rows */


extern
void SCIPunpackBasis(                   /**< unpack row and column basis status from a packed basis object */
   const BASIS*    basis,               /**< pointer to basis data */
   int             ncol,                /**< number of columns in the basis */
   int             nrow,                /**< number of rows in the basis */
   BASESTAT*       cstat,               /**< buffer for storing basis status of columns in unpacked format */
   BASESTAT*       rstat                /**< buffer for storing basis status of rows in unpacked format */
   );

extern
BASIS* SCIPcreateBasis(                 /**< creates a basis object */
   MEMHDR*         memhdr,              /**< block memory header */
   int             ncol,                /**< number of columns in the basis */
   int             nrow,                /**< number of rows in the basis */
   const BASESTAT* cstat,               /**< basis status of columns in unpacked format */
   const BASESTAT* rstat,               /**< basis status of rows in unpacked format */
   const double*   dnorm                /**< optional dual norm vector (or NULL) */
   );

extern
void SCIPfreeBasis(                     /**< free memory of basis data */
   MEMHDR*         memhdr,              /**< block memory header */
   BASIS**         basis,               /**< pointer to basis data */
   int             ncol,                /**< number of columns in the basis */
   int             nrow                 /**< number of rows in the basis */
   );

extern
void SCIPuseBasis(                      /**< announce, that the basis is referenced one more time */
   BASIS*          basis                /**< basis data */
   );

extern
void SCIPdiscardBasis(                  /**< announce, that basis is referenced one time less, free it if no longer used */
   MEMHDR*         memhdr,              /**< block memory header */
   BASIS**         basis,               /**< pointer to basis data */
   int             ncol,                /**< number of columns in the basis */
   int             nrow                 /**< number of rows in the basis */
   );


#endif

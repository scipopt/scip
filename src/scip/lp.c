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

/**@file   lp.c
 * @brief  LP management datastructures and methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "lp.h"
#include "bitencode.h"


typedef DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET DUALPACKETSIZE
typedef SINGLEPACKET ROWPACKET;         /* each row needs one bit of information (basic/nonbasic) */
#define ROWS_PER_PACKET SINGLEPACKETSIZE

struct Basis
{
   int             numuses;             /**< number of times, this warmstart information is referenced */
   COLPACKET*      packcstat;           /**< column basis status in compressed form */
   ROWPACKET*      packrstat;           /**< row basis status in compressed form */
   double*         dnorm;               /**< dual norms of variables */
};

struct Row
{
   int*            ind;                 /**< column indices of row entries */
   double*         val;                 /**< coefficients of row entries */
   double          rhs;                 /**< right hand side of row */
   double          epsilon;             /**< maximal normed violation of row */
   double          eucnorm;             /**< euclidean norm of row vector */
   double          maxval;              /**< maximal absolute value of row vector */
   int             len;                 /**< number of nonzeros in row */
   unsigned int    equality:1;          /**< TRUE iff row is an equality, FALSE iff row is a less or equal inequality */
   unsigned int    inLP:1;              /**< TRUE iff row is included in actual LP */
   unsigned int    insertLP:1;          /**< TRUE iff row should be inserted to LP */
   unsigned int    removeLP:1;          /**< TRUE iff row should be removed from LP */
   unsigned int    sorted:1;            /**< TRUE iff column indices are sorted in increasing order */
};

struct RowList
{
   ROW*            row;                 /**< pointer to this row */
   ROWLIST*        next;                /**< pointer to next rowlist entry */
};


static 
int colpacketNum(                       /**< returns the number of packets needed to store column packet information */
   int ncol                             /**< number of columns to store */
   )
{
   return (ncol+COLS_PER_PACKET-1)/COLS_PER_PACKET;
}

static 
int rowpacketNum(                       /**< returns the number of packets needed to store row packet information */
   int nrow                             /**< number of rows to store */
   )
{
   return (nrow+ROWS_PER_PACKET-1)/ROWS_PER_PACKET;
}

static
void packBasis(                         /**< store row and column basis status in a packed basis object */
   BASIS*          basis,               /**< pointer to basis data */
   int             ncol,                /**< number of columns in the basis */
   int             nrow,                /**< number of rows in the basis */
   const BASESTAT* cstat,               /**< basis status of columns in unpacked format */
   const BASESTAT* rstat                /**< basis status of rows in unpacked format */
   )
{
   assert(basis != NULL);
   assert(basis->packcstat != NULL);
   assert(basis->packrstat != NULL);
   assert(sizeof(BASESTAT) == sizeof(int));

   SCIPencodeDualBit((const int*)cstat, basis->packcstat, ncol);
   SCIPencodeSingleBit((const int*)rstat, basis->packrstat, nrow);
}

void SCIPunpackBasis(                   /**< unpack row and column basis status from a packed basis object */
   const BASIS*    basis,               /**< pointer to basis data */
   int             ncol,                /**< number of columns in the basis */
   int             nrow,                /**< number of rows in the basis */
   BASESTAT*       cstat,               /**< buffer for storing basis status of columns in unpacked format */
   BASESTAT*       rstat                /**< buffer for storing basis status of rows in unpacked format */
   )
{
   assert(basis != NULL);
   assert(basis->packcstat != NULL);
   assert(basis->packrstat != NULL);
   assert(sizeof(BASESTAT) == sizeof(int));

   SCIPdecodeDualBit(basis->packcstat, (int*)cstat, ncol);
   SCIPdecodeSingleBit(basis->packrstat, (int*)rstat, nrow);
}

BASIS* SCIPcreateBasis(                 /**< creates a basis object */
   MEMHDR*         memhdr,              /**< block memory header */
   int             ncol,                /**< number of columns in the basis */
   int             nrow,                /**< number of rows in the basis */
   const BASESTAT* cstat,               /**< basis status of columns in unpacked format */
   const BASESTAT* rstat,               /**< basis status of rows in unpacked format */
   const double*   dnorm                /**< optional dual norm vector (or NULL) */
   )
{
   BASIS* basis;

   CHECK_NULL( allocBlockMemory(memhdr, basis) );
   CHECK_NULL( allocBlockMemoryArray(memhdr, basis->packcstat, colpacketNum(ncol)) );
   CHECK_NULL( allocBlockMemoryArray(memhdr, basis->packrstat, rowpacketNum(nrow)) );
   packBasis(basis, ncol, nrow, cstat, rstat);
   if( dnorm != NULL )
      CHECK_NULL( duplicateBlockMemoryArray(memhdr, basis->dnorm, dnorm, ncol) );

   return basis;
}

void SCIPfreeBasis(                     /**< free memory of basis data */
   MEMHDR*         memhdr,              /**< block memory header */
   BASIS**         basis,               /**< pointer to basis data */
   int             ncol,                /**< number of columns in the basis */
   int             nrow                 /**< number of rows in the basis */
   )
{
   assert(memhdr != NULL);
   assert(basis != NULL);
   assert(*basis != NULL);
   assert((*basis)->numuses == 0);

   freeBlockMemoryArray(memhdr, (*basis)->packcstat, colpacketNum(ncol));
   freeBlockMemoryArray(memhdr, (*basis)->packrstat, rowpacketNum(nrow));
   freeBlockMemoryArrayNull(memhdr, (*basis)->dnorm, ncol);
   freeBlockMemory(memhdr, *basis);
}

void SCIPuseBasis(                      /**< announce, that the basis is referenced one more time */
   BASIS*          basis                /**< basis data */
   )
{
   if( basis != NULL )
      basis->numuses++;
}

void SCIPdiscardBasis(                  /**< announce, that basis is referenced one time less, free it if no longer used */
   MEMHDR*         memhdr,              /**< block memory header */
   BASIS**         basis,               /**< pointer to basis data */
   int             ncol,                /**< number of columns in the basis */
   int             nrow                 /**< number of rows in the basis */
   )
{
   assert(memhdr != NULL);
   assert(basis != NULL);

   if( *basis != NULL )
   {
      assert((*basis)->numuses > 0);
      (*basis)->numuses--;
      if( (*basis)->numuses == 0 )
         SCIPfreeBasis(memhdr, basis, ncol, nrow);
   }
}


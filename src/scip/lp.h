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

#include "lpi.h"


typedef struct Row ROW;                 /**< row of an LP */
typedef struct RowList ROWLIST;         /**< list of LP rows */

struct Lp                               /**< actual LP data */
{
   LPI*             lpi;                /**< LP solver interface */
   ROW*             rows;               /**< array with actual LP rows in correct order */
   int              rowsize;            /**< available slots in row vector */
   int              ncols;              /**< actual number of LP columns */
   int              nrows;              /**< actual number of LP rows (number of used slots in row vector) */
};
typedef struct Lp LP;


extern
ROW* SCIPcreateRow(                     /**< creates an LP row */
   MEM*             mem,                /**< block memory buffers */
   int              len,                /**< number of nonzeros in the row */
   int*             ind,                /**< column indices of row entries */
   double*          val,                /**< coefficients of row entries */
   Bool             equality,           /**< TRUE iff row is an equality, FALSE iff row is a less or equal inequality */
   double           rhs,                /**< right hand side of row */
   double           epsilon             /**< maximal normed violation of row */
   );

extern
void SCIPfreeRow(                       /**< frees an LP row */
   MEM*             mem,                /**< block memory buffers */
   ROW**            row                 /**< pointer to LP row */
   );

extern
void SCIPcaptureRow(                    /**< increases usage counter of LP row */
   ROW*             row                 /**< LP row */
   );

extern
void SCIPreleaseRow(                    /**< decreases usage counter of LP row, and frees memory if necessary */
   MEM*             mem,                /**< block memory buffers */
   ROW**            row                 /**< pointer to LP row */
   );


#endif

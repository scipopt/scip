/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2008 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: reader_col.h,v 1.3 2009/07/31 11:37:13 bzfwinkm Exp $"

/**@file   reader_col.h
 * @brief  COL file reader
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_COL_H__
#define __SCIP_READER_COL_H__

#include "scip/scip.h"
#include "tclique/tclique.h"
#include "probdata_coloring.h"
#include "scip/cons_linear.h"
#include "scip/type_cons.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** bubble sort of two joint arrays of int, sorted s.t. the second array is in non-increasing order */
extern
void COLORreaderBubbleSortIntInt(
   int*                  values,             /**< int array to be permuted in the same way */
   int*                  keys,               /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** bubble sort of two joint arrays of int and reals, sorted s.t. the second array is in non-increasing order */
extern
void COLORreaderBubbleSortIntReal(
   int*                  values,             /**< int array to be permuted in the same way */
   SCIP_Real*            keys,               /**< SCIP_Real array to be sorted */
   int                   len                 /**< length of arrays */
   );


/** creates the initial LP */
extern
SCIP_RETCODE COLORreaderCreateSetsForUncoveredNodes(
   SCIP*                 scip,               /**< SCIP data structure */   
   TCLIQUE_GRAPH*        graph               /**< pointer to graph data structure */
   );


/** includes the col file reader into SCIP */
extern
SCIP_RETCODE SCIPincludeReaderCol(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif

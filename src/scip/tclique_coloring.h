/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                        This file is part of the program                   */
/*                    TCLIQUE --- Algorithm for Maximum Cliques              */
/*                                                                           */
/*    Copyright (C) 1996-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: tclique_coloring.h,v 1.3 2005/04/29 12:56:43 bzfpfend Exp $"

/**@file   tclique_coloring.h
 * @brief  coloring part of algorithm for maximum cliques
 * @author Ralf Borndoerfer
 * @author Zoltan Kormos
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TCLIQUE_COLORING_H__
#define __TCLIQUE_COLORING_H__

#include "scip/memory.h"
#include "scip/tclique_graph.h"


typedef struct _ITV{
   int inf;
   int sup;
} ITV;

typedef struct _LIST_ITV{
   ITV itv;
   struct _LIST_ITV *next;
} LIST_ITV;

typedef struct _NBC{
   int satdeg;
   LIST_ITV *lcitv;
} NBC;




/** colors the positive weighted nodes of a given set of nodes V with the lowest possible number of colors and 
 *  finds a clique in the graph induced by V, an upper bound and an apriori bound for further branching steps
 */
extern
WEIGHT tcliqueColoring( 
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   CHKMEM*          mem,                /**< block memory */
   int*             V,                  /**< nodes for branching */ 
   int              nV,		        /**< number of nodes for branching */
   NBC*             gsd,                /**< neighbour color information of all nodes */
   BOOL*            iscolored,          /**< coloring status of all nodes */
   WEIGHT*          apbound,            /**< pointer to store apriori bound of nodes for branching */ 
   int*             clique,             /**< buffer for storing the clique */
   int*             nclique,            /**< pointer to store number of nodes in the clique */
   WEIGHT*          weightclique        /**< pointer to store the weight of the clique */
   );

#endif

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   termsepadefs.h
 * @brief  includes definitions data structures for terminal separator based methods for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/



#ifndef APPLICATIONS_STP_SRC_TERMSEPADEFS_H_
#define APPLICATIONS_STP_SRC_TERMSEPADEFS_H_

#include "scip/scip.h"
#include "graph.h"


#ifdef __cplusplus
extern "C" {
#endif

/** terminal separator tree bottleneck structure */
typedef struct terminal_separator_tree_bottleneck TBOTTLENECK;


/** separator data needed to build component */
typedef struct terminial_component_initializer
{
   SCIP_Real*            nodes_bdist;        /**< bottleneck computation distance for each node, always reset to -1.0 */
   const int*            sepaterms;          /**< separator terminals NON OWNED */
   int                   sourceterm;         /**< source terminal NOTE: we eliminate the associated sub-graph! */
   int                   nsepatterms;        /**< size of separator */
   int                   ncomponentnodes;    /**< NOTE: possibly overestimate */
   int                   componentnumber;    /**< number of component (0,1,...)*/
   int                   ngraphnodes;        /**< number of nodes of underlying graph, not counting degree 0 nodes */
   int                   maxncompchecks;     /**< maximum number of components to check */
   int                   maxsepasize;        /**< maximum allowed size of a separator */
   SCIP_Bool             rootcompIsProcessed;/**< already processed root component? */
} COMPBUILDER;


#include "stpvector.h"

/** (extended) terminal component */
typedef struct terminal_separator_component
{
   COMPBUILDER*          builder;            /**< initializer; NON-OWNED */
   GRAPH*                subgraph;           /**< graph for (extended) component */
   TBOTTLENECK*          subsolbottleneck;   /**< tree bottleneck on sub-solution */
   int*                  subsolution;        /**< primal solution for (extended) component (CONNECTED/UNKNOWN) */
   int*                  nodemap_orgToSub;   /**< map */
   int*                  nodemap_subToOrg;   /**< map */
   int*                  edgemap_subToOrg;   /**< map */
   int*                  nodes_mark;         /**< marker for nodes of component */
   STP_Vectype(int)      bfsqueue;           /**< queue for BFS */
   SCIP_Real             subprimalobj;
   int                   subnnodes;
   int                   subnedges;
} TERMCOMP;



#ifdef __cplusplus
}
#endif


#endif /* APPLICATIONS_STP_SRC_TERMSEPADEFS_H_ */

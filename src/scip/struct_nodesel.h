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
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_nodesel.h,v 1.5 2005/01/18 09:26:56 bzfpfend Exp $"

/**@file   struct_nodesel.h
 * @brief  datastructures for node selectors and node priority queues
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_NODESEL_H__
#define __STRUCT_NODESEL_H__


#include "def.h"
#include "type_tree.h"
#include "type_nodesel.h"



/** node priority queue data structure;
 *  the fields lowerboundnode, lowerbound, nlowerbounds and validlowerbound are only used for node selection rules,
 *  that don't store the lowest bound node in the first slot of the queue
 */
struct NodePQ
{
   Real             lowerboundsum;      /**< sum of lower bounds of all nodes in the queue */
   Real             lowerbound;         /**< minimal lower bound value of all nodes in the queue */
   NODESEL*         nodesel;            /**< node selector used for sorting the nodes in the queue */
   NODE**           slots;              /**< array of element slots */
   NODE*            lowerboundnode;     /**< node with minimal lower bound, or NULL if not available */
   int              len;                /**< number of used element slots */
   int              size;               /**< total number of available element slots */
   int              nlowerbounds;       /**< number of nodes in the queue with minimal lower bound (0 if invalid) */
   Bool             validlowerbound;    /**< is lower bound value valid? */
};

/** node selector */
struct Nodesel
{
   char*            name;               /**< name of node selector */
   char*            desc;               /**< description of node selector */
   DECL_NODESELFREE ((*nodeselfree));   /**< destructor of node selector */
   DECL_NODESELINIT ((*nodeselinit));   /**< initialize node selector */
   DECL_NODESELEXIT ((*nodeselexit));   /**< deinitialize node selector */
   DECL_NODESELSELECT((*nodeselselect));/**< node selection method */
   DECL_NODESELCOMP ((*nodeselcomp));   /**< node comparison method */
   NODESELDATA*     nodeseldata;        /**< node selector data */
   int              stdpriority;        /**< priority of the node selector in standard mode */
   int              memsavepriority;    /**< priority of the node selector in memory saving mode */
   Bool             lowestboundfirst;   /**< does node comparison sorts w.r.t. lower bound as primal criterion? */
   Bool             initialized;        /**< is node selector initialized? */
};


#endif

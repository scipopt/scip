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
#pragma ident "@(#) $Id: branch.h,v 1.1 2005/03/10 17:11:17 bzfpfend Exp $"

/**@file   branch.h
 * @brief  branch and bound part of algorithm for maximum cliques
 * @author Ralf Borndoerfer
 * @author Zoltan Kormos
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __BRANCH_H__
#define __BRANCH_H__

#include "tclique/graph.h"
#include "tclique/coloring.h"


/* */
#define TCLIQUE_USRCALLBACK(x) BOOL x (int* mwc, int nmwc, WEIGHT weightmwc, void* usrdata)

/** finds maximum weight clique */
void maxClique(
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   TCLIQUE_USRCALLBACK ((*usrcallback)),/**< user function to call on every new solution */
   void*            usrdata,            /**< user data to pass to user callback function */
   int*             mwc, 	        /**< pointer to store nodes of the maximum weight clique */
   int*             nmwc,	        /**< pointer to store number of nodes in the maximum weight clique */
   WEIGHT*          weightmwc,          /**< pointer to store weight of the maximum weight clique */
   WEIGHT           maxweight,          /**< maximum weight of branching nodes in level 0; 0 if not used 
                                         ** (for cliques with at least one fractional node) */
   int              maxtreenodes 	/**< maximum number of nodes of b&b tree */
   );

#endif

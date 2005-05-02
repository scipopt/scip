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
#pragma ident "@(#) $Id: tclique_branch.h,v 1.4 2005/05/02 15:55:30 bzfpfend Exp $"

/**@file   tclique_branch.h
 * @brief  branch and bound part of algorithm for maximum cliques
 * @author Ralf Borndoerfer
 * @author Zoltan Kormos
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TCLIQUE_BRANCH_H__
#define __TCLIQUE_BRANCH_H__

#include "scip/tclique_graph.h"
#include "scip/tclique_coloring.h"


/** user callback method which is called whenever a feasible clique was found
 *  input:
 *   - usrdata      : user data given to tcliqueMaxClique()
 *   - cliquenodes  : array with nodes of the clique
 *   - ncliquenodes : number of nodes in the clique
 *   - cliqueweight : weight of the clique
 *  output:
 *   - minweight    : new minimal weight for feasible cliques
 *   - acceptsol    : setting TRUE makes clique the new best clique, and updates minweight
 *   - stopsolving  : setting TRUE aborts the search for cliques
 */
#define TCLIQUE_USRCALLBACK(x) void x (void* usrdata, int* cliquenodes, int ncliquenodes, WEIGHT cliqueweight, \
      WEIGHT* minweight, BOOL* acceptsol, BOOL* stopsolving)

/** finds maximum weight clique */
extern
void tcliqueMaxClique(
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   TCLIQUE_USRCALLBACK ((*usrcallback)),/**< user function to call on every new solution */
   void*            usrdata,            /**< user data to pass to user callback function */
   int*             maxcliquenodes,     /**< pointer to store nodes of the maximum weight clique */
   int*             nmaxcliquenodes,    /**< pointer to store number of nodes in the maximum weight clique */
   WEIGHT*          maxcliqueweight,    /**< pointer to store weight of the maximum weight clique */
   WEIGHT           maxfirstnodeweight, /**< maximum weight of branching nodes in level 0; 0 if not used 
                                         *   for cliques with at least one fractional node) */
   WEIGHT           minweight,          /**< lower bound for weight of generated cliques */
   int              maxtreenodes 	/**< maximum number of nodes of b&b tree */
   );

#endif

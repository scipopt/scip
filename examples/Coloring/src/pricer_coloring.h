/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pricer_coloring.h
 * @brief  variable pricer for the vertex coloring problem
 * @author Gerald Gamrath
 *
 * This file implements the pricer for the coloring algorithm.
 *
 * It computes maximal stable sets in the current graph whose corresponding variables can improve
 * the current LP solution.  This is done by computing a maximum weighted stable set in the current
 * graph with dual-variables of the node constraints as weights. A variable can improve the
 * solution, if the weight of the corresponding stable set is larger than 1, since it then has
 * negative reduced costs, which are given by (1 - weight of the set).
 *
 * The pricer first tries to compute such a stable set using a a greedy-method. If it fails, the tclique-algorithm is
 * used on the complementary graph. This is a branch-and-bound based algorithm for maximal cliques,
 * included in SCIP.  In this case, not only the best solution is added to the LP, but also all other
 * stable sets found during the branch-and-bound process that could improve the current LP solution
 * are added, limited to a maximal number that can be changed by a parameter.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRICER_COLORING__
#define __SCIP_PRICER_COLORING__

#include "probdata_coloring.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the healthcare variable pricer and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludePricerColoring(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the way, the pricer handles variables with negative reduced costs found during the tclique-algorithm
    if onlybest is true, only the best n variables are added to the lp, while onlybest = false means, that 
    the first n variables with negative reduced costs are added
    Here, n is the value set by setNVarsCreatedPerRound */
extern
void COLORpricerUseOnlyBestStableSets(
   SCIP*                 scip,               /**< SCIP data structure */  
   SCIP_Bool             onlybest            /**< true, if only the best vars should be used */
   );


/* sets, whether the pricing should use the greedy-method */
extern
void COLORpricerUseGreedy(
   SCIP*                 scip,               /**< SCIP data structure */  
   SCIP_Bool             usegreedy           /**< true, if the greedy should be used */
   );


/* sets whether the pricing should use the tclique-method for finding new variables */
extern
void COLORpricerUseTclique(
   SCIP*                 scip,               /**< SCIP data structure */  
   SCIP_Bool             usetclique          /**< true, if the tclique-algorithm should be used */
   );


/* sets the number of variables that the pricer is allowed to create in one round of pricing */
extern
void COLORpricerSetNVarsCreatedPerRound(
   SCIP*                 scip,               /**< SCIP data structure */  
   int                   nvars               /**< maximal number of variables that should be created in one round */
   );

#ifdef __cplusplus
}
#endif

#endif

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_multilinear.h
 * @ingroup SEPARATORS
 * @brief  multilinear separator
 * @author Matthias Walter
 *
 * Separator for cutting planes that are valid for the multilinear polytope \f$ P(G) \f$ of a hypergraph
 * \f$ G = (V,E) \f$, defined as the convex hull of all vectors \f$ z \in \{0,1\}^{V + E} \f$ that
 * satisfy \f$ z_e = \prod_{v \in e} z_v \f$ for all \f$ e \in E \f$. In other words, the variable associated with each
 * (hyper-)edge is equal to the product of the variables associated with the vertices.
 *
 * The hypergraph \f$ G \f$ is computed in the first separation round in which this separator is called. The edges
 * correspond to all AND constraints (see \ref cons_and.h), where \f$ z_e \f$ is the resultant \f$ r \f$ of
 * \f$ r = x_1 \land x_2 \land \dotsb \land x_n \f$ and its incident vertices correspond to the \f$ x_i \f$.
 * Moreover, also (nonlinear) product expressions (see \ref expr_product.h) are scanned. For these, the involved
 * variables are not necessarily binary, i.e., we have \f$ \ell_i \leq x_i \leq u_i \f$, and there might be an
 * additional constant scaling factor \f$ c \f$. Such expressions \f$ r = c \prod\limits_{i=1}^n x_i \f$ are only taken
 * into account if \f$ \ell_i \geq 0 \f$ for all \f$ i \f$ since in this case the cut can be applied to the space
 * \f$ (r',x') \f$ with \f$ r' = \frac{ r }{ c \cdot u_1 \cdot u_2 \cdot \dotsc \cdot u_n } \f$ and
 * \f$ x'_i = \frac{ x_i }{ u_i } \f$ for all \f$ i \f$.
 *
 * ### Flower inequalities ###
 *
 * The implemented cuts are <em>k-flower inequalities</em> for \f$ k=1,2 \f$. Such a cut is determined by a
 * <em>base edge</em> \f$ e \f$ and \f$ k \f$ adjacent edges \f$ f_1, f_2, \dotsc, f_k \f$ that satisfy
 * \f$ f_i \cap e \neq \emptyset \f$ but are disjoint in \f$ e \f$, i.e., \f$ f_i \cap f_j \cap e = \emptyset \f$ for
 * all \f$ i \neq j \f$ holds. The set \f$ R \f$ of <em>remaining nodes</em> is defined as
 * \f$ R := e \setminus \bigcup_{i=1}^k f_i \f$. The inequality reads
 * \f[ z_e + \sum\limits_{i=1}^k (1-z_{f_i}) + \sum\limits_{v \in R} (1-z_v) \geq 1. \f]
 * Validity follows from the fact that the left-hand side is a sum of nonnegative binary terms, and can thus only
 * be violated if (in particular) \f$ z_{f_i} = 1 \f$ for \f$ i=1,2,\dotsc,k \f$ and \f$ z_v = 1 \f$ for all
 * \f$ v \in R \f$ holds. This, however, implies \f$ z_v = 1 \f$ for all \f$ v \in e \f$ and thus \f$ z_e = 1 \f$.
 *
 * Separation can either be done in time \f$ \mathcal{O}( |E|^{k+1} ) \f$ by enumeration or as follows:
 * Replacing an adjacent edge \f$ f_i \f$ by another edge \f$ f_i' \f$ with the same <em>overlap</em> with the base
 * edge, i.e., \f$ f'_i \cap e = f_i \cap e \f$ improves the violation if \f$ 1-z_{f'_i} < 1 - z_{f_i} \f$.
 * Consequently, we compute the set of all <em>overlap sets</em> \f$ \{ e \cap f \mid e,f \in E \} \f$
 * (see \ref hypergraph.h) and compute \f$ \min \{ 1-z_e \mid e \in E : U \subseteq e \} \f$ for all such \f$ U \f$.
 * If the number of overlap sets incident to an edge is constant (say, if \f$ |e| \f$ is constant), then the running
 * time reduces to \f$ \mathcal{O} ( |E| ) \f$.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_MULTILINEAR_H__
#define __SCIP_SEPA_MULTILINEAR_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the multilinear separator and includes it in SCIP
 *
 * @ingroup SeparatorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeSepaMultilinear(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   WeightSpaceVertex.h
 * @brief  Weight space vertex
 * @author Timo Strunk
 *
 * Data structure storing combinatorial and geometric information about a vertex of the weight space polyhedron
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef CLASS_WEIGHTSPACEVERTEX
#define CLASS_WEIGHTSPACEVERTEX

#include <vector>
#include <set>

#include "scip/def.h"

#undef GCC_VERSION
#include "lemon/list_graph.h"

/** data structure for a vertex of the weight space polyhedron */
class WeightSpaceVertex
{
 public:
   /** creates inital vertex */
   WeightSpaceVertex(
      std::vector< const std::vector<SCIP_Real>* >     incident_facets,  /**< incident nondominated points */
      std::vector<SCIP_Real>*                          weight,           /**< weight vector */
      SCIP_Real                                        weighted_objval   /**< weighted objective value */
      );

   /** creates a new point between obsolete and adjacent non obsolete point */
   WeightSpaceVertex(
      const WeightSpaceVertex*          obsolete,           /**< vertex cut off by new solution */
      const WeightSpaceVertex*          adjacent,           /**< adjacent non obsolete vertex */
      const std::vector<SCIP_Real>*     new_facet           /**< new solution cutting off the obsolete vertex */
      );

   /** creates dummy point */
   WeightSpaceVertex();

   /** destructor */
   ~WeightSpaceVertex();

   /** whether this vertex and argument vertex are neighbours in the 1-skeleton*/
   bool isNeighbour(
      const WeightSpaceVertex*          vertex              /**< another weight space vertex */
      ) const;

   /** returns the weighted objective value */
   SCIP_Real getWeightedObjectiveValue() const;

   /** returns the weight */
   const std::vector<SCIP_Real>* getWeight() const;

   /** returns the set of facets defining the vertex */
   const std::vector< const std::vector<SCIP_Real>* >* getFacets() const;

   /** returns the graph node associated with the vertex */
   lemon::ListGraph::Node getNode() const ;

   /** sets the graph node associated with the vertex */
   void setNode(
      lemon::ListGraph::Node  node /**< corresponding node in skeleton graph */
      );

   bool isCorner() const;
   void updateFacet(const std::vector<SCIP_Real>* facet);

/** writes weight space vertex to an output stream */
   void print(
      std::ostream&                   os        /** stream the vector should be written to*/
      ) const;

 private:
   unsigned int                                   nobjs_;                       /**< number of objectives */
   std::vector< const std::vector<SCIP_Real>* >   incident_facets_;             /**< defining facets of form (w,a)*coeffs >= 0 */
   std::vector<SCIP_Real>*                        weight_;                      /**< weight vector */
   SCIP_Real                                      weighted_objective_value_;    /**< weighted objective value */
   lemon::ListGraph::Node                         node_;                        /**< associated graph node */

   /** determines the set of facets based on other vertices' facets */
   void joinFacets(
      const WeightSpaceVertex*          obsolete,           /**< vertex cut off by new solution */
      const WeightSpaceVertex*          adjacent,           /**< adjacent non obsolete vertex */
      const std::vector<SCIP_Real>*     new_facet           /**< new solution cutting off the obsolete vertex */
      );

   /** calculates the weight w and the weighted objective value a
    *  based on w and a for the obsolete and the adjacent vertex */
   void calculate_weight(
      const WeightSpaceVertex*          obsolete,           /**< vertex cut off by new solution */
      const WeightSpaceVertex*          adjacent,           /**< adjacent non obsolete vertex */
      const std::vector<SCIP_Real>*     new_facet           /**< new solution cutting off the obsolete vertex */
   );

};

/** writes weight space vertex to an output stream */
std::ostream& operator<<(
   std::ostream&                   os,       /** stream the vector should be written to*/
   const WeightSpaceVertex         v         /** vertex that should be written */
   );

#endif

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
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
 * @desc  Data structure storing combinatorial and geometric information about a vertex of the weight space polyhedron
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef CLASS_WEIGHTSPACEVERTEX
#define CLASS_WEIGHTSPACEVERTEX

#include <vector>
#include <set>

#include "scip/def.h"
#include "lemon/list_graph.h"

/** data structure for a vertex of the weight space polyhedron */
class WeightSpaceVertex
{
 public:
   /** creates initial corner point */
   WeightSpaceVertex(
      unsigned int                      dimension,          /**< dimension of the weight space */
      const std::vector<SCIP_Real>*     first_sol,          /**< first solution defining the initial facette */
      unsigned int                      nonzeroDimension    /**< axis where the vertex is not zero */
      );

   /** creates a new point between obsolete and adjacent non obsolete point */
   WeightSpaceVertex(
      const WeightSpaceVertex*          obsolete,           /**< vertex cut off by new solution */
      const WeightSpaceVertex*          adjacent,           /**< adjacent non obsolete vertex */
      const std::vector<SCIP_Real>*     new_sol             /**< new solution cutting off the obsolete vertex */
      );

   /** creates dummy point */
   WeightSpaceVertex();

   /** destructor */
   ~WeightSpaceVertex();

   /** whether this represents a corner of the weight space */
   bool isCorner() const;

   /** whether this and point are neighbours in the 1-skeleton*/
   bool isNeighbour(
      const WeightSpaceVertex*          point               /**< another weight space vertex */
      );

   /** changes the adjacent nondominated cost vector of a corner vertex */
   void updateNondomPoint(
      const std::vector<SCIP_Real>*     new_nondom_point    /**< nondom point replacing the old one */
      );

   /** returns the weighted objective value */
   SCIP_Real getWeightedObjectiveValue() const ;

   /** returns the weight */
   const std::vector<SCIP_Real>* getWeight() const ;
   
   /** returns the set of axes where the weight is greater than 0 */
   const std::vector<unsigned int> * getNonZeroDimensions() const ;
   
   /** returns the set of nondominated points that define the vertex */
   const std::vector< const std::vector<SCIP_Real>* >* getIncidentSolutions() const ;
   
   /** returns the graph node associated with the vertex */
   lemon::ListGraph::Node getNode() const ;
   
   /** sets the graph node associated with the vertex */
   void setNode(
      lemon::ListGraph::Node  node /**< corresponding node in skeleton graph */
      );

 private:
   unsigned int                                   nobjs_;                       /**< number of objectives */
   std::vector< const std::vector<SCIP_Real>* >   incident_solutions_;          /**< incident nondominated points */
   std::vector<unsigned int>                      nonzero_dimensions_;          /**< axes where weight is > 0 */
   std::vector<SCIP_Real>*                        weight_;                      /**< weight vector */
   SCIP_Real                                      weighted_objective_value_;    /**< weighted objective value */
   lemon::ListGraph::Node                         node_;                        /**< associated graph node */

   /** sets incident solutions to intersection of given points' incident solutions */
   void joinIncidentSolutions(
      const WeightSpaceVertex*          obsolete,           /**< vertex cut off by new solution */
      const WeightSpaceVertex*          adjacent,           /**< adjacent non obsolete vertex */
      const std::vector<SCIP_Real>*     new_sol             /**< new solution cutting off the obsolete vertex */
      );

   /** sets nonzero dimensions to union of given points nonzero dimensions */
   void joinNonzeroDimensions(
      const WeightSpaceVertex*          obsolete,           /**< vertex cut off by new solution */
      const WeightSpaceVertex*          adjacent            /**< adjacent non obsolete vertex */
      );

/** calculates the weight w and the weighted objective value a 
 *  based on w and a for the obsolete and the adjacent vertex */
   void calculate_weight(
      const WeightSpaceVertex*          obsolete,           /**< vertex cut off by new solution */
      const WeightSpaceVertex*          adjacent,           /**< adjacent non obsolete vertex */      
      const std::vector<SCIP_Real>*     new_sol             /**< new solution cutting off the obsolete vertex */
   );

};

#endif

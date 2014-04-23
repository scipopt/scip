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

/**@file   Skeleton.h
 * @brief  Weight space polyhedron
 * @author Timo Strunk
 *
 * @desc   This class represents the lifted weight space polyhedron.  It supplies weights for the solver to test.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef CLASS_WEIGHTGRAPH
#define CLASS_WEIGHTGRAPH

#include <map>
#include <set>
#include <vector>
#include <queue>

#include "scip/scip.h"
#include "scip/def.h"
#include "lemon/list_graph.h"

class WeightSpaceVertex;

/** the 1-skeleton of the weight space polyhedron */
class Skeleton
{
 public:
   /** constructor */
   Skeleton(
      SCIP*              scip                /**< SCIP solver */
      );

   /** destructor */
   ~Skeleton();

   /** whether there is an untested weight left */
   bool hasNextWeight();

   /** get the next untested weight */
   const std::vector<SCIP_Real> * nextWeight();

   /** determines whether the solution is a new pareto optimum, if so, adds it to the polyhedron */
   bool checkSolution(
      const std::vector<SCIP_Real>*     cost_vector    /**< cost vector that is a candidate to be a nondominated point */
      );

   /** get number of vertices added in last checkSolution call*/
   int getNNewVertices() const;

   /** get number of vertices processed in last checkSolution call*/
   int getNProcessedVertices() const;

 private:
   SCIP*                                          scip_;              /**< SCIP solver */
   std::set<lemon::ListGraph::Node>               untested_nodes_;    /**< nodes for which no weighted run has been
								        * performed yet */
   lemon::ListGraph                               graph_;             /**< the graph structure */
   lemon::ListGraph::NodeMap<WeightSpaceVertex*>  vertex_map_;        /**< map from graph nodes to polygon vertex data */
   lemon::ListGraph::Node                         last_node_;         /**< last tested node */
   std::vector<WeightSpaceVertex*>                vertices_;          /**< list of all generated vertices */
   int                                            n_new_vertices_;    /**< number of vertices added in last 
								       *   checkSolution call*/
   int                                            n_proc_vertices_;   /**< number of vertices processed in last 
								       *   checkSolution call*/

   /* data structures for temporary use in update step*/
   const std::vector<SCIP_Real>*                  new_nondom_point_;  /**< new nondominated cost vector */
   std::vector<WeightSpaceVertex*>*               new_vertices_;      /**< new generated vertices */
   std::queue<lemon::ListGraph::Node>*            unscanned_nodes_;   /**< nodes left to scan for obsolecity */
   std::set<lemon::ListGraph::Node>*              obsolete_nodes_;    /**< nodes identified as obsolete */
   std::vector<lemon::ListGraph::Edge>*           cut_edges_;         /**< edges from obsolete to nonobsolete nodes */

   /** init the polyhedron with the first solution */
   void init(
      const std::vector<SCIP_Real>*     first_nondom_point  /**< non dominated point defining the first facet */ 
      );

   /** wether the new solution makes a given weight space vertex obsolete */
   bool makesObsolete(
      const std::vector<SCIP_Real>*     cost_vector,        /**< cost vector of a solution */
      const WeightSpaceVertex*          vertex              /**< vertex that might be obsolete */ 
      );

   /** updates the polyhedron with the new solution*/
   void addFacet(
      const std::vector<SCIP_Real>*     nondom_point        /**< new nondominated cost vector */
      );

   /** tests all the neighbours of an obsolete node for obsolecity and then removes the node */
   void processObsoleteNode(
      lemon::ListGraph::Node            obs_node            /**< a skeleton node marked as obsolete */
      );

   /** apply changes calculated by add facet */
   void updateGraph();

   /** creates a new node between an obsolete and a non obsolete node */
   void makeIntermediaryPoint( 
      lemon::ListGraph::Edge            cut_edge            /**< an edge between an obsolete and a non-obsolete node */
      );

   /** special method dealing with obsolete nodes that are also corners of the weight space */
   void updateCorner(
      WeightSpaceVertex*                obsolete_vertex     /**< an obsolete vertex that is a corner */
      );

   /** adds a graph node corresponding to new_vertex */
   lemon::ListGraph::Node addNode(
      WeightSpaceVertex*                new_vertex          /**< new vertex not yet represented in the graph */
      );

   /** returns true if data is combinatorically consistent */
   bool graphIsValid() const;
};

#endif

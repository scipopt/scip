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

/**@file   Skeleton.cpp
 * @brief  Weight space polyhedron
 * @author Timo Strunk
 *
 * @desc   This class represents the lifted weight space polyhedron.  It supplies weights for the solver to test.
 * It uses the lemon graph library to store the 1-skeleton of the polyhedron.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <set>

#include "Skeleton.h"
#include "reader_mop.h"
#include "Objectives.h"
#include "WeightSpaceVertex.h"
#include "main.h"

/** constructor */
Skeleton::Skeleton(
   SCIP*                 scip                /**< SCIP solver */
   )
   : scip_(scip),
     graph_(),
     vertex_map_(graph_)
{

}

/** destructor */
Skeleton::~Skeleton()
{
   for( std::vector<WeightSpaceVertex*>::iterator it = vertices_.begin();
        it != vertices_.end();
        ++it )
   {
      delete *it;
   }
}

/** whether there is an untested weight left */
bool Skeleton::hasNextWeight()
{
   return !untested_nodes_.empty();
}

/** get the next untested weight */
const std::vector<SCIP_Real>* Skeleton::nextWeight()
{
   last_node_ = *untested_nodes_.begin();
   untested_nodes_.erase(last_node_);

   assert(vertex_map_[last_node_] != NULL);

   return (vertex_map_[last_node_]->getWeight());
}

/** returns true and updates the polyhedron if the solution is a new pareto optimum */
bool Skeleton::checkSolution(
   const std::vector<SCIP_Real>*        cost_vector    /**< cost vector that is a candidate to be a nondominated point */
   )
{
   bool                  result;

   if( lemon::ListGraph::NodeIt(graph_) == lemon::INVALID )
   {
      /* cost_vector is the first nondom_point */
      init(cost_vector);
      result = true;
      n_new_vertices_ = cost_vector->size();
      n_proc_vertices_ = 0;
   }
   else
   {
      /* cost_vector is a nondominated point if and only if it has a better weighted objective value
       * than all previous nondominated points with respect to the weight it optimizes */
      result = makesObsolete(cost_vector, vertex_map_[last_node_]);
      if(result)
      {
         addFacet(cost_vector);
      }
      else
      {
         n_new_vertices_ = 0;
         n_proc_vertices_ = 1;
      }
   }
   return result;
}

/** init the polyhedron with the first solution */
void Skeleton::init(
   const std::vector<SCIP_Real>*        first_nondom_point  /**< non dominated point defining the first facet */ 
   )
{
   int                        nobjs;
   WeightSpaceVertex*         vertex;
   lemon::ListGraph::Node     node;

   new_vertices_ = new std::vector<WeightSpaceVertex*>();
   nobjs = SCIPgetProbData(scip_)->objectives->getNObjs();

   /* create a node and vertex for every corner of the weight space and a complete graph between them */
   for( int i = 0; i < nobjs; i++ )
   {
      vertex = new WeightSpaceVertex(nobjs, first_nondom_point, i);
      node = addNode(vertex);

      for( std::vector<WeightSpaceVertex*>::iterator it = new_vertices_->begin();
           it != new_vertices_->end();
           ++it )
      {
         graph_.addEdge(node, (**it).getNode());
      }

      new_vertices_->push_back(vertex);
      vertices_.push_back(vertex);
   }

   assert(graphIsValid());

   delete new_vertices_;
}

/** weather the new solution makes a given point obsolete */
bool Skeleton::makesObsolete(
   const std::vector<SCIP_Real>*        cost_vector,        /**< cost vector of a solution */
   const WeightSpaceVertex*             vertex              /**< vertex that might be obsolete */ 
   )
{
   SCIP_Real                       weighted_objective_value;
   const std::vector<SCIP_Real>*   weight;

   assert(cost_vector != NULL);
   assert(vertex != NULL);

   weighted_objective_value = 0.;
   weight = vertex->getWeight();

   /* calculate the weighted objective value of cost_vector */
   for( unsigned int j = 0; j < cost_vector->size(); j++ )
   {
      weighted_objective_value += (*weight)[j] * cost_vector->at(j);
   }

   /* compare weighted objective values */
   if( SCIPgetObjsense(scip_) == SCIP_OBJSENSE_MAXIMIZE )
   {
      return SCIPisGT(
         scip_,
         weighted_objective_value,
         vertex->getWeightedObjectiveValue()
         );
   }
   else
   {
      return SCIPisLT(
         scip_,
         weighted_objective_value,
         vertex->getWeightedObjectiveValue()
         );
   }
}

/** updates the polyhedron with the new solution */
void Skeleton::addFacet(
   const std::vector<SCIP_Real>*        nondom_point        /**< new nondominated cost vector */
   )
{
   lemon::ListGraph::Node     obs_node;

   new_nondom_point_ = nondom_point;
   new_vertices_     = new std::vector<WeightSpaceVertex*>();
   unscanned_nodes_  = new std::queue<lemon::ListGraph::Node>();
   obsolete_nodes_   = new std::set<lemon::ListGraph::Node>();
   cut_edges_        = new std::vector<lemon::ListGraph::Edge>();

   obsolete_nodes_->insert(last_node_);
   unscanned_nodes_->push(last_node_);

   /* scan all obsolete nodes */
   while( !unscanned_nodes_->empty() )
   {
      obs_node = unscanned_nodes_->front();
      unscanned_nodes_->pop();
      processObsoleteNode( obs_node);
   }

   updateGraph();

   delete cut_edges_;
   delete obsolete_nodes_;
   delete unscanned_nodes_;
   delete new_vertices_;
}

/** tests all the neighbours of an obsolete node for obsolecity and then removes the node */
void Skeleton::processObsoleteNode(
   lemon::ListGraph::Node               obs_node            /**< a skeleton node marked as obsolete */
   )
{
   lemon::ListGraph::Node     neighbour;
   WeightSpaceVertex*         adjacent_vertex;

   /* iterate over all neighbours of obs_node */
   for( lemon::ListGraph::IncEdgeIt e(graph_, obs_node); e != lemon::INVALID; ++e )
   {
      neighbour = graph_.oppositeNode(obs_node, e);
      adjacent_vertex = vertex_map_[neighbour];

      if( obsolete_nodes_->find(neighbour) == obsolete_nodes_->end() )
      {
         /* neighbour is not yet marked as obsolete */
         if( makesObsolete(new_nondom_point_, adjacent_vertex) )
         {
            /* neighbour is obsolete, add it to the queue */
            obsolete_nodes_->insert(neighbour);
            unscanned_nodes_->push(neighbour);
            untested_nodes_.erase(neighbour);
         }
         else
         {
            /* neighbour is not obsolete, edge is in cut */
            cut_edges_->push_back(e);
         }
      }
   }
}

void Skeleton::updateGraph()
{
   WeightSpaceVertex*    obsolete_vertex;

   n_proc_vertices_ = obsolete_nodes_->size();
   n_new_vertices_ = cut_edges_->size();

   for( std::vector<lemon::ListGraph::Edge>::iterator it = cut_edges_->begin();
        it != cut_edges_->end();
        ++it )
   {
      makeIntermediaryPoint(*it);
   }

   for( std::set<lemon::ListGraph::Node>::iterator it = obsolete_nodes_->begin();
        it != obsolete_nodes_->end();
        ++it )
   {
      obsolete_vertex = vertex_map_[*it];
      if( obsolete_vertex->isCorner() )
      {
         updateCorner(obsolete_vertex);
      }
      graph_.erase(*it);
   }
   assert(graphIsValid());
}

/** creates a new node between an obsolete and a non obsolete node */
void Skeleton::makeIntermediaryPoint( 
   lemon::ListGraph::Edge               cut_edge            /**< an edge between an obsolete and a non-obsolete node */
   )
{
   lemon::ListGraph::Node     new_node;
   lemon::ListGraph::Node     obsolete_node;
   lemon::ListGraph::Node     adjacent_node;
   WeightSpaceVertex*         obsolete_vertex;
   WeightSpaceVertex*         adjacent_vertex;
   WeightSpaceVertex*         new_vertex;

   /* find out which end of the edge is the obsolete node */
   if( obsolete_nodes_->find(graph_.u(cut_edge)) == obsolete_nodes_->end() )
   {
      /* u is not the obsolete node */
      adjacent_node = graph_.u(cut_edge);
      obsolete_node = graph_.v(cut_edge);
   }
   else
   {
      adjacent_node = graph_.v(cut_edge);
      obsolete_node = graph_.u(cut_edge);
   }
   obsolete_vertex = vertex_map_[obsolete_node];
   adjacent_vertex = vertex_map_[adjacent_node];

   new_vertex = new WeightSpaceVertex(
      obsolete_vertex,
      adjacent_vertex,
      new_nondom_point_);
   vertices_.push_back(new_vertex);
   new_node = addNode(new_vertex);
   graph_.addEdge(new_node, adjacent_node);

   /* add a graph edge for all pairs of combinatorially adjacent new vertices */
   for( std::vector<WeightSpaceVertex*>::iterator pit = new_vertices_->begin();
        pit < new_vertices_->end(); ++pit )
   {
      if( (**pit).isNeighbour(new_vertex) )
      {
         graph_.addEdge((**pit).getNode(), new_node);
      }
   }
   new_vertices_->push_back(new_vertex);
}

/** special method dealing with obsolete nodes that are also corners of the weight space */
void Skeleton::updateCorner(
   WeightSpaceVertex*    obsolete_vertex     /**< an obsolete vertex that is a corner */
   )
{
   lemon::ListGraph::Node new_node;
   WeightSpaceVertex*  adjacent_vertex;

   bool is_tested_node_ = (obsolete_vertex->getNode() == last_node_);
   new_node = graph_.addNode();

   /* marry node and vertex */
   obsolete_vertex->setNode(new_node);
   vertex_map_[new_node] = obsolete_vertex;

   obsolete_vertex->updateNondomPoint(new_nondom_point_);

   /* insert new edges */
   for( std::vector<WeightSpaceVertex*>::iterator pit = new_vertices_->begin();
        pit < new_vertices_->end();
        ++pit )
   {
      adjacent_vertex = *pit;
      if( adjacent_vertex->isNeighbour(obsolete_vertex) )
      {
         graph_.addEdge(adjacent_vertex->getNode(), obsolete_vertex->getNode());
      }
   }
   /* mark updated corner as new vertex */
   new_vertices_->push_back(obsolete_vertex);
   if( !is_tested_node_ )
   {
      untested_nodes_.insert(new_node);
      ++n_new_vertices_;
   }
}

/** adds a graph node corresponding to new_vertex */
lemon::ListGraph::Node Skeleton::addNode(
   WeightSpaceVertex*    new_vertex          /**< new vertex not yet represented in the graph */
   )
{
   lemon::ListGraph::Node result;

   result = graph_.addNode();
   untested_nodes_.insert(result);

   /* marry node and vertex */
   new_vertex->setNode(result);
   vertex_map_[result] = new_vertex;

   return result;
}

/** returns true if all graph edges are actual polyhedron edges */
bool Skeleton::graphIsValid() const
{
   WeightSpaceVertex*    u;
   WeightSpaceVertex*    v;

   for( lemon::ListGraph::EdgeIt e(graph_); e != lemon::INVALID; ++e )
   {
      u = vertex_map_[graph_.u(e)];
      v = vertex_map_[graph_.v(e)];

      if( !(u->isNeighbour(v)) )
      {
         return false;
      }
   }
   return true;
}

/** get number of vertices added in last checkSolution call*/
int Skeleton::getNNewVertices() const
{
   return n_new_vertices_;
}

/** get number of vertices processed in last checkSolution call*/
int Skeleton::getNProcessedVertices() const
{
   return n_proc_vertices_;
}

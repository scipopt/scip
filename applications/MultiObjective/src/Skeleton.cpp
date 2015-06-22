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

/**@file   Skeleton.cpp
 * @brief  Weight space polyhedron
 * @author Timo Strunk
 *
 * This class represents the lifted weight space polyhedron.  It supplies weights for the solver to test.
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
     vertex_map_(graph_),
     last_returned_node_(lemon::INVALID),
     n_new_nodes_(0),
     n_proc_nodes_(0),
     new_facet_(NULL),
     new_vertices_(NULL),
     unscanned_nodes_(NULL),
     obsolete_nodes_(NULL),
     cut_edges_(NULL)
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

   for( std::vector< const std::vector<SCIP_Real>* >::iterator
           it = facets_.begin();
        it != facets_.end();
        ++it )
   {
      delete *it;
   }
}

/** initialize the polyhedron with the first solution
 * by creating a node and vertex for every corner of the weight space and a complete graph between them */
void Skeleton::init(
   const std::vector<SCIP_Real>*                   cost_vector,  /**< cost vector of first solution */
   std::vector< const std::vector<SCIP_Real>* >*   cost_rays     /**< list of known unbounded cost rays */
   )
{
   int                        nobjs = SCIPgetProbData(scip_)->objectives->getNObjs();
   WeightSpaceVertex*         vertex;
   lemon::ListGraph::Node     node;

   createInitialFacets(cost_vector);

   new_vertices_ = new std::vector<WeightSpaceVertex*>();

   for( int i = 0; i < nobjs; i++ )
   {
      vertex = createCorner(i);
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

   delete new_vertices_;

   if( cost_rays != NULL )
   {
      addPrimalRays(cost_rays);
   }

   assert(graphIsValid());
}

/** create all facets defining the inital weight space polyhedron */
void Skeleton::createInitialFacets(
   const std::vector<SCIP_Real>*                   cost_vector   /**< cost vector of first solution */
   )
{
   int nobjs = SCIPgetProbData(scip_)->objectives->getNObjs();

   /* create weight space boundary facets */
   std::vector<SCIP_Real>*    facet;
   for( int i = 0; i < nobjs; ++i )
   {
      facet = new std::vector<SCIP_Real>(nobjs + 1, 0.);
      (*facet)[i] = 1.;
      facets_.push_back(facet);
   }

   /* create first nondom point based facet */
   facets_.push_back(createFacetFromCost(cost_vector));
}

/** create corner vertex of initial weight space polyhedron */
WeightSpaceVertex* Skeleton::createCorner(
   int                                  index               /**< index where weight is 1 */
   )
{
   int nobjs = SCIPgetProbData(scip_)->objectives->getNObjs();

   /* corner is defined by one nondom point and p-1 weight space boundaries*/
   std::vector< const std::vector<SCIP_Real>* > defining_facets(facets_.begin(), facets_.begin() + nobjs);
   const std::vector<SCIP_Real>*                first_nondom_facet = facets_[nobjs];
   defining_facets[index] = first_nondom_facet;

   /* corner has weight (0,..,0,1,0,...,0) */
   std::vector<SCIP_Real>*                      weight             = new std::vector<SCIP_Real>(nobjs, 0.);
   (*weight)[index] = 1.;

   return new WeightSpaceVertex(defining_facets, weight, (*first_nondom_facet)[index]);
}

/** whether there is an untested weight left */
bool Skeleton::hasNextWeight()
{
   return !untested_nodes_.empty();
}

/** get the next untested weight */
const std::vector<SCIP_Real>* Skeleton::nextWeight()
{
   assert( !untested_nodes_.empty() );

   last_returned_node_ = *(untested_nodes_.begin());
   untested_nodes_.erase(last_returned_node_);

   assert(vertex_map_[last_returned_node_] != NULL);

   return (vertex_map_[last_returned_node_]->getWeight());
}

/** returns true and updates the polyhedron if cost vector is a new nondominated point.
 * It is a nondominated point if and only if it has a better weighted objective value
 * than all previous nondominated points with respect to the weight it optimizes */
bool Skeleton::isExtremal(
   const std::vector<SCIP_Real>*        cost_vector    /**< potential new nondominated point */
   )
{
   /* graph must be initialized */
   assert( lemon::ListGraph::NodeIt(graph_) != lemon::INVALID );

   const std::vector<SCIP_Real>* facet = createFacetFromCost(cost_vector);

   bool result = isMakingObsolete(facet, vertex_map_[last_returned_node_], false);

   if(result)
   {
      addFacet(facet);
   }
   else
   {
      n_new_nodes_ = 0;
      n_proc_nodes_ = 1;
      delete facet;
   }

   return result;
}

/** like isExtremal but check all vertices for obsolecity (not just the last returned one)*/
bool Skeleton::isExtremalThorough(
   const std::vector<SCIP_Real>*     cost_vector         /**< cost vector that is a candidate to be a nondominated point */
   )
{
   /* graph must be initialized */
   assert( lemon::ListGraph::NodeIt(graph_) != lemon::INVALID );

   bool result = false;
   const std::vector<SCIP_Real>* facet = createFacetFromCost(cost_vector);

   last_returned_node_ = findObsoleteNode(facet);

   if( last_returned_node_ != lemon::INVALID )
   {
      result = isMakingObsolete(facet, vertex_map_[last_returned_node_], false);
   }

   if(result)
   {
      addFacet(facet);
   }
   else
   {
      n_new_nodes_ = 0;
      n_proc_nodes_ = 1;
      delete facet;
   }

   return result;
}

/** adds a weight space constraint after finding a primal ray with unbounded weighted objective */
void Skeleton::addPrimalRay(
   const std::vector<SCIP_Real>*        cost_ray            /**< cost vector of the unbounded primal ray */
   )
{
   /* graph must be initialized */
   assert( lemon::ListGraph::NodeIt(graph_) != lemon::INVALID );

   const std::vector<SCIP_Real>* facet = createFacetFromRay(cost_ray);
   addFacet(facet);
}

/** like addPrimalRay but check all vertices for obsolecity (not just the last returned one)*/
void Skeleton::addPrimalRayThorough(
   const std::vector<SCIP_Real>*        cost_ray            /**< cost vector of the unbounded primal ray */
   )
{
   /* graph must be initialized */
   assert( lemon::ListGraph::NodeIt(graph_) != lemon::INVALID );

   const std::vector<SCIP_Real>* facet = createFacetFromRay(cost_ray);
   last_returned_node_ = findObsoleteNode(facet);

   if( last_returned_node_ != lemon::INVALID )
   {
      addFacet(facet);
   }
}

/** add multiple rays with unbounded weighted objective */
void Skeleton::addPrimalRays(
   const std::vector< const std::vector<SCIP_Real>* >* cost_rays /**< rays to add */
   )
{
   /* variables for caching update statistics */
   int n_new_nodes_tmp = n_new_nodes_;
   int n_proc_nodes_tmp = n_proc_nodes_;

   for( std::vector< const std::vector<SCIP_Real>* >::const_iterator
           it = cost_rays->begin();
        it != cost_rays->end();
        ++it )
   {
      addPrimalRayThorough(*it);
      n_new_nodes_tmp += n_new_nodes_;
      n_proc_nodes_tmp += n_proc_nodes_;
   }

   n_new_nodes_ = n_new_nodes_tmp;
   n_proc_nodes_ = n_proc_nodes_tmp;
}

/** returns a node corresponding to a vertex made obsolete by facet or INVALID */
lemon::ListGraph::Node Skeleton::findObsoleteNode(
   const std::vector<SCIP_Real>* facet
  )
{
   lemon::ListGraph::Node result = lemon::INVALID;
   for( lemon::ListGraph::NodeIt it(graph_); it != lemon::INVALID; ++it )
   {
      if( isMakingObsolete(facet, vertex_map_[it]) )
      {
         result = it;
         break;
      }
   }

   return result;
}

/** wether the new facet makes a given point obsolete */
bool Skeleton::isMakingObsolete(
   const std::vector<SCIP_Real>*        facet,              /**< facet coefficient vector */
   const WeightSpaceVertex*             vertex,             /**< vertex that might be obsolete */
   bool                                 strict              /**< no tolerance for slight obsolecity */
   )
{
   assert(facet != NULL);
   assert(vertex != NULL);

   int                           nobjs  = SCIPgetProbData(scip_)->objectives->getNObjs();
   const std::vector<SCIP_Real>* weight = vertex->getWeight();

   /* calculate the left hand side value of the vertex with regard to the facet inequality */
   SCIP_Real                     lhs    = (*facet)[nobjs] * vertex->getWeightedObjectiveValue();

   for( int j = 0; j < nobjs; ++j )
   {
      lhs += (*weight)[j] * (*facet)[j];
   }

   /* negative left hand side means the inequality is not satisfied */
   if( strict )
   {
      return lhs < 0;
   }
   else
   {
      return SCIPisLT(scip_, lhs, 0);
   }
}

/** updates the polyhedron with the new solution */
void Skeleton::addFacet(
   const std::vector<SCIP_Real>*        facet        /**< new nondominated cost vector */
   )
{
   assert( lemon::ListGraph::NodeIt(graph_) != lemon::INVALID );

   lemon::ListGraph::Node     obs_node;

   new_facet_ = facet;
   facets_.push_back(facet);

   /* update some stats */
   n_proc_nodes_ = 0;
   n_new_nodes_  = 0;
   untested_nodes_.erase(last_returned_node_);

   /* nodes made obsolete by new facet */
   obsolete_nodes_   = new std::set<lemon::ListGraph::Node>();
   obsolete_nodes_->insert(last_returned_node_);

   /* obsolete nodes whose incident edges have not yet been scanned */
   unscanned_nodes_  = new std::queue<lemon::ListGraph::Node>();
   unscanned_nodes_->push(last_returned_node_);

   /* edges connecting obsolete and non-obsolete nodes */
   cut_edges_        = new std::vector<lemon::ListGraph::Edge>();

   /* scan all obsolete nodes */
   while( !unscanned_nodes_->empty() )
   {
      obs_node = unscanned_nodes_->front();
      unscanned_nodes_->pop();
      scanNode(obs_node);
      ++n_proc_nodes_;
   }

   /* apply changes */
   updateGraph();

   delete cut_edges_;
   delete unscanned_nodes_;
   delete obsolete_nodes_;


   assert(graphIsValid());
}

/** tests all the neighbours of an obsolete node for obsolecity and then removes the node */
void Skeleton::scanNode(
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
         if( isMakingObsolete(new_facet_, adjacent_vertex, true) )
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

/** apply changes calculated by add facet */
void Skeleton::updateGraph()
{
   new_vertices_     = new std::vector<WeightSpaceVertex*>();

   createNewVertices();

   createNewEdges();

   for( std::set<lemon::ListGraph::Node>::iterator it = obsolete_nodes_->begin();
        it != obsolete_nodes_->end();
        ++it )
   {
      graph_.erase(*it);
   }

   delete new_vertices_;
}

/** calculate new vertices from obsolete vertices and add them to the graph */
void Skeleton::createNewVertices()
{
   WeightSpaceVertex*    obsolete_vertex;

   for( std::vector<lemon::ListGraph::Edge>::iterator it = cut_edges_->begin();
        it != cut_edges_->end();
        ++it )
   {
      makeIntermediateVertex(*it);
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
   }
}

/** calculate and add edges between all pairs of combinatorially adjacent new vertices */
void Skeleton::createNewEdges()
{
   for( std::vector<WeightSpaceVertex*>::iterator qit = new_vertices_->begin();
        qit < new_vertices_->end(); ++qit )
   {
      for( std::vector<WeightSpaceVertex*>::iterator pit = new_vertices_->begin();
           pit < qit; ++pit )
      {
         if( (**pit).isNeighbour(*qit) )
         {
            graph_.addEdge((**pit).getNode(), (**qit).getNode());
         }
      }
   }
}

/** creates a new node between an obsolete and a non obsolete node */
void Skeleton::makeIntermediateVertex(
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

   /* create vertex */
   new_vertex = new WeightSpaceVertex(
      obsolete_vertex,
      adjacent_vertex,
      new_facet_);
   vertices_.push_back(new_vertex);
   new_vertices_->push_back(new_vertex);

   /* update graph */
   new_node = addNode(new_vertex);
   graph_.addEdge(new_node, adjacent_node);
}

/** special method dealing with obsolete nodes that are also corners of the weight space */
void Skeleton::updateCorner(
   WeightSpaceVertex*    obsolete_vertex     /**< an obsolete vertex that is a corner */
   )
{
   lemon::ListGraph::Node new_node;

   bool untested = (obsolete_vertex->getNode() != last_returned_node_);
   new_node = addNode(obsolete_vertex, untested);

   obsolete_vertex->updateFacet(new_facet_);

   /* mark updated corner as new vertex */
   new_vertices_->push_back(obsolete_vertex);
}

/** adds a graph node corresponding to new_vertex */
lemon::ListGraph::Node Skeleton::addNode(
   WeightSpaceVertex*    new_vertex,         /**< new vertex not yet represented in the graph */
   bool                  mark_untested       /**< set true if weight should be tested at some point */
   )
{
   lemon::ListGraph::Node result = graph_.addNode();

   /* marry node and vertex */
   new_vertex->setNode(result);
   vertex_map_[result] = new_vertex;

   if( mark_untested )
   {
      untested_nodes_.insert(result);
      ++n_new_nodes_;
   }

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

/** get number of vertices added in last isExtremal call*/
int Skeleton::getNNewVertices() const
{
   return n_new_nodes_;
}

/** get number of vertices processed in last isExtremal call*/
int Skeleton::getNProcessedVertices() const
{
   return n_proc_nodes_;
}

/** returns facet vector corresponding to cost vector */
const std::vector<SCIP_Real>* Skeleton::createFacetFromCost(
   const std::vector<SCIP_Real>*        cost_vector    /**< cost vector of a solution */
   ) const
{
   std::vector<SCIP_Real>* result = new std::vector<SCIP_Real>(*cost_vector);
   result->push_back(-1.);
   return result;
}

/** returns facet vector corresponding to an unbounded cost ray */
const std::vector<SCIP_Real>* Skeleton::createFacetFromRay(
   const std::vector<SCIP_Real>*        cost_ray      /**< cost vector of a primal ray */
   ) const
{
   std::vector<SCIP_Real>* result = new std::vector<SCIP_Real>(*cost_ray);
   result->push_back(0.);
   return result;
}

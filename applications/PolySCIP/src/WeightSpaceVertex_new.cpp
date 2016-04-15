/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*        This file is part of the program PolySCIP                          */
/*                                                                           */
/*    Copyright (C) 2012-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  PolySCIP is distributed under the terms of the ZIB Academic License.     */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with PolySCIP; see the file LICENCE.                               */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   WeightSpaceVertex.cpp
 * @brief  Weight space vertex
 * @author Sebastian Schenker, Timo Strunk
 *
 * Data structure storing combinatorial and geometric information about a vertex of the weight space polyhedron
 */

#include "WeightSpaceVertex.h"
#include "scip/def.h"

#include <algorithm> // std::copy, std::set_intersection
#include <iterator>  // std::ostream_iterator, std::inserter
#include <numeric>   // std::inner_product;
#include <stdio.h>

using std::set;
using std::vector;

/** creates inital vertex */
WeightSpaceVertex::WeightSpaceVertex(
   const vector<unsigned>& incident_facet_inds,  /**< indices of incident facets */
   const vector<SCIP_Real>* weight,                /**< weight vector */
   SCIP_Real weighted_objval                       /**< weighted objective value */
  ) : nObjs_(weight.size()),
      facet_indices_(incident_facet_inds.begin(),incident_facet_inds.end()),
      weight_(weight),
      weighted_obj_val_(weighted_objval)
{
  assert (nObjs_ > 0);
  assert (facet_indices_.size() == nObjs_);
}

/** creates a new point between obsolete and adjacent non obsolete point */
WeightSpaceVertex::WeightSpaceVertex(
   const WeightSpaceVertex& obs,         /**< vertex cut off by new solution */
   const WeightSpaceVertex& adj,         /**< adjacent non obsolete vertex */
   const vector<SCIP_Real>& new_facet,   /**< new solution cutting off the obsolete vertex */
   bool inclinationToAdj                 /**< whether to incline computed weight slightly towards adjacent vertex */
				     ) 
  : nObjs_(obs.getNObjs())
{
   joinFacets(obsolete, adjacent, new_facet);
   calculate_weight(obsolete, adjacent, new_facet);
   assert (weight_ != NULL);
}

/** destructor */
WeightSpaceVertex::~WeightSpaceVertex() {
  if (weight_ != NULL)
    delete weight_;
}

unsigned WeightSpaceVertex::getNObjs() const {
  return nObjs_;
}

/** whether this and point are neighbours in the 1-skeleton */
bool WeightSpaceVertex::isNeighbour(const WeightSpaceVertex& point) const {
  set<unsigned> common_facets;
  set<unsigned>* otherFacets = point.getFacets();

  std::set_intersection(facet_indices_.begin(),
			facet_indices_.end(),
			otherFacets->begin(),
			otherFacets->end(),
			std::inserter( common_facets, common_facets.begin() ) );

  return common_facets.size() == nobjs_-1;
}

/** returns the weighted objective value */
SCIP_Real WeightSpaceVertex::getWeightedObjVal() const {
  return weighted_obj_val_;
}

/** returns the weight vector */
const vector<SCIP_Real>* WeightSpaceVertex::getWeight() const {
  return weight_;
}

/** return i-th element of weight vector */
SCIP_Real WeightSpaceVertex::getWeight(unsigned i) const {
  return (*weight_)[i];
}

/** returns the set of indices of facets that define the vertex */
const set<unsigned>* WeightSpaceVertex::getFacets() const {
   return &facet_indices_;
}

/** returns the graph node associated with the vertex */
lemon::ListGraph::Node WeightSpaceVertex::getNode() const {
   return node_;
}

/** sets the graph node associated with the vertex */
void WeightSpaceVertex::setNode(lemon::ListGraph::Node node /**< corresp. node in skeleton graph */ 
				) {
  node_ = node;
}

/** set facet indices to intersection of obsolet and adjacent vertex plus new facet */
void WeightSpaceVertex::joinFacets(
   const WeightSpaceVertex& obs,           /**< vertex cut off by new solution */
   const WeightSpaceVertex& adj,           /**< adjacent non obsolete vertex */
   unsigned new_facet_index                     /**< new solution cutting off the obsolete vertex */
   )
{
  const set<unsigned>* obsFacets = obs.getFacets(),
    adjFacets = adj.getFacets();
  assert (!obsFacets.empty());
  assert (!adjFacets.empty());
  assert (facet_indices_.empty());
   
  set_intersection(obsFacets->begin(),
		   obsFacets->end(),
		   adjFacets->begin(),
		   adjFacets->end(),
		   std::inserter( facet_indices_, facet_indices_.begin() ) );
  facet_indices.insert(new_facet);
  assert (facet_indices_.size() == nObjs_);
}

/** calculates the weight w and the weighted objective value a
 *  based on w and a for the obsolete and the adjacent vertex
 *  by solving the equation system
 *  I) (w,a)         = h * (w_obs, a_obs) + (1 - h) * (w_adj, a_adj)
 *  II)(w,a) * facet = 0
 *  through insertion of I) into II) and solving for h
 */
void  WeightSpaceVertex::calculate_weight(
      const WeightSpaceVertex& obs,      /**< vertex cut off by new solution */
      const WeightSpaceVertex& adj,      /**< adjacent non obsolete vertex */
      const vector<SCIP_Real>* new_facet,     /**< new solution cutting off the obsolete vertex */
      bool inclinationToAdj                   /**< whether to incline computed weight slightly towards adjacent vertex */
					  ) {
  const vector<SCIP_Real>* obsWeight = obs.getWeight(),
    adjWeight = adj.getWeight();
  assert (obsWeight != NULL);
  assert (adjWeight != NULL);
  assert (new_facet != NULL);
  assert (obsWeight->size() == nObjs_);
  assert (adjWeight->size() == nObjs_);
  assert (new_facet->size() == nObjs_+1);
  
  SCIP_Real lhs_obs = std::inner_product(obsWeight->begin(),
					 obsWeight->end(),
					 new_facet->begin(),
					 obs.getWeightedObjVal()*new_facet->back()),
    lhs_adj = std::inner_product(adjWeight->begin(),
				 adjWeight->end(),
				 new_facet->begin(),
				 adj.getWeightedObjVal()*new_facet->back()),
    factor_adj = lhs_obs / (lhs_obs - lhs_adj);
  
  if (inclinationToAdj)    /* incline weight slightly towards adjacent solution to */
    factor_adj += 0.00001; /* overcome possible unboundedness when checking vertex */
   
  SCIP_Real factor_obs = 1. - factor_adj;

  weighted_objective_value_ = factor_obs*obs.getWeightedObjVal() + factor_adj*adj.getWeightedObjVal();

  weight_ = new std::vector<SCIP_Real>(nObjs_);
  for (unsigned i=0; i<nObjs_; ++i)
    (*weight_)[i] = factor_obs*obs.getWeight(i) + factor_adj*abj.getWeight(i);
}

// bool WeightSpaceVertex::isCorner() const {
// int n_point_facets = 0;
//   for( std::vector< const std::vector<SCIP_Real>* >::const_iterator it = incident_facets_.begin();
//        it != incident_facets_.end();
//        ++it )
//   {
//     if( (**it)[nobjs_] == -1. )
//     {
//       ++n_point_facets;
//     }
//   }

//   return  n_point_facets < 2;
// }

// void WeightSpaceVertex::updateFacet(const std::vector<SCIP_Real>* facet) {
//   assert( isCorner() );
//   assert( weight_->size() == facet->size()-1 );
   
//    for (unsigned i=0; i<nobjs_; ++i) {
//      if ((*incident_facets_[i])[nobjs_] == -1.) {
//        incident_facets_[i] = facet;
//        std::sort(incident_facets_.begin(), incident_facets_.end());
       
//        weighted_objective_value_ = std::inner_product(weight_->begin(), weight_->end(),
// 						      facet->begin(), 0.0);
       
//        return;
//      }
//    }
// }

/** writes weight space vertex to an output stream */
void WeightSpaceVertex::print(std::ostream& os, const vector< vector<SCIP_Real>* >& facets) const {
  set<unsigned>* indices = getFacets();
  vector<SCIP_Real>* weight = getWeight();

  os << "incident facets: " << std::endl;
  std::ostream_iterator<SCIP_Real> os_it(os, " ");
  for (auto it=indices->begin(); it!=indices->end(); ++it) {
    os << "[ ";
    std::copy(facets[*it]->begin(), facets[*it]->end(), os_it);
    os << "]\n";
  }
  os << "weight = [ ";
  std::copy(weight->begin(), weight->end(),os_it);
  os << "]\n";
  os << "weighted objective value = " << getWeightedObjVal() << std::endl;
}


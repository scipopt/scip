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

/**@file   WeightSpaceVertex.cpp
 * @brief  Weight space vertex
 * @author Timo Strunk
 *
 * Data structure storing combinatorial and geometric information about a vertex of the weight space polyhedron
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <iostream>

#include "scip/def.h"

#include "main.h"
#include "WeightSpaceVertex.h"

/** creates inital vertex */
WeightSpaceVertex::WeightSpaceVertex(
   std::vector< const std::vector<SCIP_Real>* >     incident_facets,  /**< incident nondominated points */
   std::vector<SCIP_Real>*                          weight,           /**< weight vector */
   SCIP_Real                                        weighted_objval   /**< weighted objective value */
  ) : nobjs_(weight->size()),
      incident_facets_(incident_facets),
      weight_(weight),
      weighted_objective_value_(weighted_objval)
{
   sort(incident_facets_.begin(), incident_facets_.end());

   assert( nobjs_ > 0 );
}

/** creates a new point between obsolete and adjacent non obsolete point */
WeightSpaceVertex::WeightSpaceVertex(
   const WeightSpaceVertex*          obsolete,           /**< vertex cut off by new solution */
   const WeightSpaceVertex*          adjacent,           /**< adjacent non obsolete vertex */
   const std::vector<SCIP_Real>*     new_facet           /**< new solution cutting off the obsolete vertex */
  ) : nobjs_(obsolete->nobjs_)
{
   joinFacets(obsolete, adjacent, new_facet);
   calculate_weight(obsolete, adjacent, new_facet);

   assert( weight_->size() == nobjs_ );
}

/** default constructor */
WeightSpaceVertex::WeightSpaceVertex()
   : nobjs_(0),
     weight_(NULL)
{
}

/** destructor */
WeightSpaceVertex::~WeightSpaceVertex()
{
   if( weight_ != NULL )
   {
      delete weight_;
   }
}

/** whether this and point are neighbours in the 1-skeleton */
bool WeightSpaceVertex::isNeighbour(const WeightSpaceVertex* point) const
{
   std::vector<const std::vector<SCIP_Real> * > common_facets;
   std::vector<const std::vector<SCIP_Real> * >::iterator cit;

   common_facets.resize(nobjs_);
   cit = set_intersection(
      point->incident_facets_.begin(),
      point->incident_facets_.end(),
      incident_facets_.begin(),
      incident_facets_.end(),
      common_facets.begin());

   return cit - common_facets.begin() >= nobjs_ - 1;
}

/** returns the weighted objective value */
SCIP_Real WeightSpaceVertex::getWeightedObjectiveValue() const
{
   return weighted_objective_value_;
}

/** returns the weight */
const std::vector<SCIP_Real> * WeightSpaceVertex::getWeight() const
{
   return weight_;
}

/** returns the set of nondominated points that define the vertex */
const std::vector< const std::vector<SCIP_Real>* >* WeightSpaceVertex::getFacets() const
{
   return &incident_facets_;
}

/** returns the graph node associated with the vertex */
lemon::ListGraph::Node WeightSpaceVertex::getNode() const
{
   return node_;
}

/** sets the graph node associated with the vertex */
void WeightSpaceVertex::setNode(
   lemon::ListGraph::Node     node /**< corresponding node in skeleton graph */
   )
{
   node_ = node;
}

/** sets incident solutions to intersection of given points incident solutions */
void WeightSpaceVertex::joinFacets(
   const WeightSpaceVertex*             obsolete,           /**< vertex cut off by new solution */
   const WeightSpaceVertex*             adjacent,           /**< adjacent non obsolete vertex */
   const std::vector<SCIP_Real>*        new_facet           /**< new solution cutting off the obsolete vertex */
   )
{
   const std::vector<const std::vector<SCIP_Real>* >&    facetsObs = obsolete->incident_facets_;
   const std::vector<const std::vector<SCIP_Real>* >&    facetsAdj = adjacent->incident_facets_;

   std::vector<const std::vector<SCIP_Real> *>::iterator cit;

   incident_facets_.resize(nobjs_ - 1);
   cit = set_intersection(
      facetsObs.begin(),
      facetsObs.end(),
      facetsAdj.begin(),
      facetsAdj.end(),
      incident_facets_.begin());
   incident_facets_.push_back(new_facet);
   sort(incident_facets_.begin(), incident_facets_.end());
}

/** calculates the weight w and the weighted objective value a
 *  based on w and a for the obsolete and the adjacent vertex
 *  by solving the equation system
 *  I) (w,a)         = h * (w_obs, a_obs) + (1 - h) * (w_adj, a_adj)
 *  II)(w,a) * facet = 0
 *  through insertion of I) into II) and solving for h
 */
/** sets nonzero dimensions to union of given points nonzero dimensions */
void  WeightSpaceVertex::calculate_weight(
      const WeightSpaceVertex*          obsolete,           /**< vertex cut off by new solution */
      const WeightSpaceVertex*          adjacent,           /**< adjacent non obsolete vertex */
      const std::vector<SCIP_Real>*     new_facet           /**< new solution cutting off the obsolete vertex */
  )
{
   const std::vector<SCIP_Real>* weight_obs = obsolete->getWeight();
   const std::vector<SCIP_Real>* weight_adj = adjacent->getWeight();
   SCIP_Real wov_obs = obsolete->getWeightedObjectiveValue();
   SCIP_Real wov_adj = adjacent->getWeightedObjectiveValue();

   assert( weight_obs != NULL );
   assert( weight_adj != NULL );
   assert( new_facet != NULL );
   assert( weight_obs->size() == nobjs_ );
   assert( weight_adj->size() == nobjs_ );
   assert( new_facet->size()  == nobjs_ + 1 );

   SCIP_Real lhs_obs = wov_obs * (*new_facet)[nobjs_];
   SCIP_Real lhs_adj = wov_adj * (*new_facet)[nobjs_];

   for( unsigned int i = 0; i < nobjs_; ++i )
   {
     lhs_obs += (*weight_obs)[i] * (*new_facet)[i];
     lhs_adj += (*weight_adj)[i] * (*new_facet)[i];
   }

   SCIP_Real factor_adj = lhs_obs / (lhs_obs - lhs_adj);
   SCIP_Real factor_obs = 1. - factor_adj;


   weight_ = new std::vector<SCIP_Real>(nobjs_);
   weighted_objective_value_ = factor_obs * wov_obs + factor_adj * wov_adj;

   for( unsigned int i = 0; i < nobjs_; ++i )
   {
     (*weight_)[i] += factor_obs * (*weight_obs)[i] + factor_adj * (*weight_adj)[i];
   }

}

bool WeightSpaceVertex::isCorner() const
{
  int n_point_facets = 0;
  for( std::vector< const std::vector<SCIP_Real>* >::const_iterator it = incident_facets_.begin();
       it < incident_facets_.end();
       ++it )
  {
    if( (**it)[nobjs_] == -1. )
    {
      ++n_point_facets;
    }
  }

  return  n_point_facets < 2;
}

void WeightSpaceVertex::updateFacet(const std::vector<SCIP_Real>* facet)
{
   assert( isCorner() );

   for( unsigned int i = 0; i < nobjs_; ++i )
   {
      if( (*incident_facets_[i])[nobjs_] == -1. )
      {
         incident_facets_[i] = facet;
         weighted_objective_value_ = 0.;
         for( unsigned int i = 0; i < nobjs_; ++i )
         {
            weighted_objective_value_ += (*weight_)[i] * (*facet)[i];
         }
         return;
      }
   }

   assert( weight_->size() == nobjs_ );
}

/** writes weight space vertex to an output stream */
void WeightSpaceVertex::print(
   std::ostream&                   os        /**< stream the vector should be written to*/
   ) const
{
   os << "facets" << std::endl;
   for( std::vector< const std::vector<SCIP_Real>* >::const_iterator it = incident_facets_.begin();
        it != incident_facets_.end();
        ++it )
   {
      os << **it << std::endl;
   }
   os << "coords" << std::endl;
   os << *weight_ << ", " << weighted_objective_value_;
}

/** writes weight space vertex to an output stream */
std::ostream& operator<<(
   std::ostream&                   os,       /**< stream the vector should be written to*/
   const WeightSpaceVertex         v         /**< vertex that should be written */
   )
{
   v.print(os);
   return os;
}

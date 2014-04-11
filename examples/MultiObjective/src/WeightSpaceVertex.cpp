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

/**@file   WeightSpaceVertex.cpp
 * @brief  Weight space vertex
 * @author Timo Strunk
 * 
 * @desc  Data structure storing combinatorial and geometric information about a vertex of the weight space polyhedron
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "WeightSpaceVertex.h"
#include <stdio.h>
#include "scip/def.h"
#include <iostream>
#include "main.h"

/** creates initial corner point */
WeightSpaceVertex::WeightSpaceVertex(
   unsigned int                      dimension,          /**< dimension of the weight space */
   const std::vector<SCIP_Real>*     first_sol,          /**< first solution defining the initial facette */
   unsigned int                      nonzeroDimension    /**< axis where the vertex is not zero */
   )
   : nobjs_(dimension)
{
   weight_ = new std::vector<SCIP_Real>(dimension,0.);
   incident_solutions_.push_back(first_sol);
   nonzero_dimensions_.push_back(nonzeroDimension);
   (*weight_)[nonzeroDimension] = 1.;
   weighted_objective_value_ = (*first_sol)[nonzeroDimension];
}

/** creates a new point between obsolete and adjacent non obsolete point */
WeightSpaceVertex::WeightSpaceVertex(
   const WeightSpaceVertex*        obsolete,           /**< vertex cut off by new solution */
   const WeightSpaceVertex*        adjacent,           /**< adjacent non obsolete vertex */
   const std::vector<SCIP_Real>*   new_sol             /**< new solution cutting off the obsolete vertex */
   )
   : nobjs_(obsolete->nobjs_)
{
   weight_ = new std::vector<SCIP_Real>(nobjs_,0.);

   joinNonzeroDimensions(obsolete, adjacent);
   joinIncidentSolutions(obsolete, adjacent, new_sol);

   assert(incident_solutions_.size() == nonzero_dimensions_.size());

   calculate_weight(obsolete, adjacent, new_sol);
}

/** creates dummy point */
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

/** whether this represents a corner of the weight space */
bool WeightSpaceVertex::isCorner() const 
{
   return nonzero_dimensions_.size()==1;
}

/** whether this and point are neighbours in the 1-skeleton*/
bool WeightSpaceVertex::isNeighbour(const WeightSpaceVertex* point)
{
   std::vector<unsigned int> nonzero_difference;
   std::vector<unsigned int>::iterator dit; 
   unsigned int dims_difference;
   std::vector<const std::vector<SCIP_Real> * > incidence_difference;
   std::vector<const std::vector<SCIP_Real> * >::iterator cit;
   unsigned int cell_difference;

   nonzero_difference.resize(nobjs_);
   dit = set_symmetric_difference(
      point->nonzero_dimensions_.begin(),
      point->nonzero_dimensions_.end(), 
      nonzero_dimensions_.begin(), 
      nonzero_dimensions_.end(), 
      nonzero_difference.begin());
   dims_difference = (unsigned int)(dit - nonzero_difference.begin());

   incidence_difference.resize(2 * nobjs_);
   cit = set_symmetric_difference(
      point->incident_solutions_.begin(),
      point->incident_solutions_.end(), 
      incident_solutions_.begin(), 
      incident_solutions_.end(),
      incidence_difference.begin());
   cell_difference = (unsigned int)(cit - incidence_difference.begin());

   return cell_difference + dims_difference == 2;
}

/** changes the adjacent nondominated cost vector of a corner vertex */
void WeightSpaceVertex::updateNondomPoint(
   const std::vector<SCIP_Real>*   new_nondom_point    /**< nondom point replacing the old one */
   )
{
   assert(nonzero_dimensions_.size() == 1);
   assert(incident_solutions_.size() == 1);

   incident_solutions_.clear();
   incident_solutions_.push_back(new_nondom_point);
   weighted_objective_value_ = (*new_nondom_point)[nonzero_dimensions_[0]];

   assert(incident_solutions_.size() == 1);
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

/** returns the set of axes where the weight is greater than 0 */
const std::vector<unsigned int> * WeightSpaceVertex::getNonZeroDimensions() const 
{
   return &nonzero_dimensions_;
}
  
/** returns the set of nondominated points that define the vertex */
const std::vector< const std::vector<SCIP_Real>* >* WeightSpaceVertex::getIncidentSolutions() const
{
   return &incident_solutions_;
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

/** sets nonzero dimensions to union of given points nonzero dimensions */
void WeightSpaceVertex::joinNonzeroDimensions(
   const WeightSpaceVertex*             obsolete,           /**< vertex cut off by new solution */
   const WeightSpaceVertex*             adjacent            /**< adjacent non obsolete vertex */
   )
{
   const std::vector<unsigned int>* nonzeroObs;
   const std::vector<unsigned int>* nonzeroAdj;
   std::vector<unsigned int>::iterator dit;

   nonzeroObs = obsolete->getNonZeroDimensions();
   nonzeroAdj = adjacent->getNonZeroDimensions();
   nonzero_dimensions_.resize(nobjs_);
   dit = set_union( 
      nonzeroObs->begin(), 
      nonzeroObs->end(), 
      nonzeroAdj->begin(),
      nonzeroAdj->end(), 
      nonzero_dimensions_.begin());
   nonzero_dimensions_.resize(dit - nonzero_dimensions_.begin());
   sort(nonzero_dimensions_.begin(), nonzero_dimensions_.end());
}

/** sets incident solutions to intersection of given points incident solutions */
void WeightSpaceVertex::joinIncidentSolutions(
   const WeightSpaceVertex*             obsolete,           /**< vertex cut off by new solution */
   const WeightSpaceVertex*             adjacent,           /**< adjacent non obsolete vertex */
   const std::vector<SCIP_Real>*        new_sol             /**< new solution cutting off the obsolete vertex */
   )
{
   const std::vector<const std::vector<SCIP_Real>* >* incidentObs;
   const std::vector<const std::vector<SCIP_Real>* >* incidentAdj;
   std::vector<const std::vector<SCIP_Real> *>::iterator cit;

   incidentObs = obsolete->getIncidentSolutions();
   incidentAdj = adjacent->getIncidentSolutions();
   incident_solutions_.resize(std::max(incidentObs->size(), incidentAdj->size()));
   cit = set_intersection(
      incidentObs->begin(), 
      incidentObs->end(), 
      incidentAdj->begin(),
      incidentAdj->end(), 
      incident_solutions_.begin());
   incident_solutions_.resize(cit - incident_solutions_.begin());
   incident_solutions_.push_back(new_sol);
   sort(incident_solutions_.begin(), incident_solutions_.end());
}

/** calculates the weight w and the weighted objective value a 
 *  based on w and a for the obsolete and the adjacent vertex
 *  by solving the equation system
 *  I) (a,w) = h * (a_obs, w_obs) + (1 - h) * (a_adj, w_adj)
 *  II) a = w * y_new
 *  through insertion of I) into II) and solving for h
 */
void WeightSpaceVertex::calculate_weight(
      const WeightSpaceVertex*          obsolete,           /**< vertex cut off by new solution */
      const WeightSpaceVertex*          adjacent,           /**< adjacent non obsolete vertex */      
      const std::vector<SCIP_Real>*     new_sol             /**< new solution cutting off the obsolete vertex */
  )
{
   const std::vector<SCIP_Real>* weight_obs = obsolete->getWeight();
   const std::vector<SCIP_Real>* weight_adj = adjacent->getWeight();
   SCIP_Real wov_obs = obsolete->getWeightedObjectiveValue();
   SCIP_Real wov_adj = adjacent->getWeightedObjectiveValue();

   SCIP_Real excess_obs = wov_obs;
   SCIP_Real excess_adj = wov_adj;

   for( unsigned int i = 0; i < nobjs_; ++i )
   {
     excess_obs -= (*weight_obs)[i] * (*new_sol)[i];
     excess_adj -= (*weight_adj)[i] * (*new_sol)[i];
   }

   SCIP_Real factor_obs = excess_adj / (excess_adj - excess_obs);
   SCIP_Real factor_adj = 1 - factor_obs;

   weighted_objective_value_ = factor_obs * wov_obs + factor_adj * wov_adj;

   for( unsigned int i = 0; i < nobjs_; ++i )
   {
     (*weight_)[i] += factor_obs * (*weight_obs)[i] + factor_adj * (*weight_adj)[i];
   }
}



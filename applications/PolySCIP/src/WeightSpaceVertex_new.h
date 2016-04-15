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

/**@file   WeightSpaceVertex.h
 * @brief  Weight space vertex
 * @author Sebastian Schenker, Timo Strunk
 *
 * Data structure storing combinatorial and geometric information about a vertex of the weight space polyhedron
 */

#ifndef CLASS_WEIGHTSPACEVERTEX
#define CLASS_WEIGHTSPACEVERTEX

#undef GCC_VERSION /* lemon/core.h redefines GCC_VERSION additionally to scip/def.h */
#include "lemon/list_graph.h"
#include "scip/def.h"

#include <set>
#include <vector>

/** data structure for a vertex of the weight space polyhedron */
class WeightSpaceVertex
{
 public:
   /** creates inital vertex */
  WeightSpaceVertex(
		    const std::vector<unsigned>& incident_facet_inds, /**< indices of incident facets */
		    const std::vector<SCIP_Real>* weight,          /**< corresponding weight vector */
		    SCIP_Real weighted_objval                      /**< weighted objective value */
		    );

  /** creates a new point between obsolete and adjacent non obsolete point */
  /* WeightSpaceVertex( */
  /*     const WeightSpaceVertex& obsolete,        /\**< vertex cut off by new solution *\/ */
  /*     const WeightSpaceVertex& adjacent,        /\**< adjacent non obsolete vertex *\/ */
  /*     const std::vector<SCIP_Real>& new_facet,  /\**< new solution cutting off the obsolete vertex *\/ */
  /*     bool inclinationToAdj                     /\**< whether to incline weight slightly towards adjacent vertex *\/ */
  /*     ); */

  /** destructor */
  ~WeightSpaceVertex();
  
  /** whether this vertex and argument vertex are neighbours in the 1-skeleton*/
  bool isNeighbour(
		    const WeightSpaceVertex*          vertex              /**< another weight space vertex */
      ) const;

   unsigned getNObjs() const;

   /** returns the weighted objective value */
   SCIP_Real getWeightedObjVal() const;

   /** returns the weight vector */
   const std::vector<SCIP_Real>* getWeight() const;
   
   /** returns i-th element of weight vector */
   SCIP_Real getWeight(unsigned i) const;

   /** returns the set of indices of facets defining the vertex */
   const std::set<unsigned>* getFacets() const;

   /* /\** returns the graph node associated with the vertex *\/ */
   /* lemon::ListGraph::Node getNode() const ; */

   /* /\** sets the graph node associated with the vertex *\/ */
   /* void setNode( */
   /*    lemon::ListGraph::Node  node /\**< corresponding node in skeleton graph *\/ */
   /*    ); */

   //bool isCorner() const;
   //void updateFacet(const std::vector<SCIP_Real>* facet);

/** writes weight space vertex to an output stream */
   void print(
	      std::ostream& os,        /** stream the vector should be written to*/
	      const std::vector< std::vector<SCIP_Real>* >& all_facets, /** all facets of the weight space polyhedron */
	      ) const;

 private:
   unsigned nObjs_;                      /**< number of objectives */
   std::set<unsigned> facet_indices_;    /**< indices of defining facets of form (w,a)*coeffs >= 0 */
   std::vector<SCIP_Real>* weight_;      /**< weight vector */
   SCIP_Real weighted_obj_val_;          /**< weighted objective value */
   //lemon::ListGraph::Node node_;         /**< associated graph node */

   /* /\** set facet indices to intersection of obsolet and adjacent vertex plus new facet *\/ */
   /* void joinFacets ( */
   /*    const WeightSpaceVertex&          obsolete,         /\**< vertex cut off by new solution *\/ */
   /*    const WeightSpaceVertex&          adjacent,         /\**< adjacent non obsolete vertex *\/ */
   /*    unsigned                          new_facet_index   /\**< new solution cutting off the obsolete vertex *\/ */
   /*    ); */

   /* /\** calculates the weight w and the weighted objective value a */
   /*  *  based on w and a for the obsolete and the adjacent vertex *\/ */
   /* void calculate_weight( */
   /*    const WeightSpaceVertex& obs,              /\**< vertex cut off by new solution *\/ */
   /*    const WeightSpaceVertex& adj,              /\**< adjacent non obsolete vertex *\/ */
   /*    const std::vector<SCIP_Real>* new_facet,   /\**< new solution cutting off the obsolete vertex *\/ */
   /*    bool inclinationToAdj = false              /\**< whether to incline computed weight slightly towards adjacent vertex *\/ */
   /* ); */

};

#endif

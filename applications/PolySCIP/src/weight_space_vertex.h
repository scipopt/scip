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

/** * @brief  Class representing a weight space vertex
 *
 * Data structure storing combinatorial and geometric information
 * about a vertex of the weight space polyhedron. A weight space
 * vertex is represented by a weight 'w' and an weighted objective
 * value 'a' and is a vertex of the (partial) weight space polyhedron
 * P = {(w,a) \in \Lambda \times R : w \cdot y >= a \forall y \in Y'}
 * where Y' is the (current) set of non-dominated points and \Lambda
 * is the set of normalized weights
 */

#ifndef POLYSCIP_SRC_WEIGHT_SPACE_VERTEX_H_INCLUDED
#define POLYSCIP_SRC_WEIGHT_SPACE_VERTEX_H_INCLUDED

#include <memory> // std::shared_ptr
#include <vector>

#include "polyscip.h"
#include "weight_space_facet.h"

namespace polyscip {
  
  /** Vertex of the weight space polyhedron. */
  class WeightSpaceVertex {

  public:
    using FacetContainer = std::vector<std::shared_ptr<const WeightSpaceFacet>>;

    /** Creates a vertex of the (partial) weight space polyhedron.
     *  @param incident_facets Facets defining of the weight space polyhedron defining the vertex
     *  @param weight Corresponding weight
     *  @param weighted_obj_val Corresponding maximal weight objective val in weight space polyhedron
     */
    explicit WeightSpaceVertex(FacetContainer incident_facets, 
			       std::shared_ptr<const Polyscip::WeightType> weight,         
			       Polyscip::ValueType weighted_obj_);
    
    /** Destructor */
    ~WeightSpaceVertex();
  
    /** Checks adjacency of two vertices.
     *  @param other_vertex another weight space vertex
     *  @return true if other_vertex is adjacent, false otherwise
     */
    //    bool isAdjacent(std::shared_ptr<const WeightSpaceVertex> other_vertex) const;

    /** Returns associated weighted objective value.
     *  @return weighted objective value
     */
    Polyscip::ValueType getWeightedObjVal() const;

    /** Checks whether weight of vertex corresponds to unit weight 
     *  @param index index of 1 in unit weight
     *  @return true if weight of vertex is unit weight with 1 at index; false otherwise
     */
    bool hasUnitWeight(Polyscip::WeightType::size_type index) const;
   
    /** Prints weight space vertex information to standard output. 
     *  @param printFacets if true, then defining facets are printed
     */
    void print(bool printFacets = false) const;

  private:
    FacetContainer incident_facets_;                     /**< incident facets */
    std::shared_ptr<const Polyscip::WeightType> weight_; /**< used weight */
    Polyscip::ValueType weighted_obj_val_;               /**< corresponding weighted objective value */

  };

}

#endif // POLYSCIP_SRC_WEIGHT_SPACE_VERTEX_H_INCLUDED

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

/**@file   weight_space_facet.h
 * @brief  Weight space facet class declarations
 * @author Sebastian Schenker
 *
 * Data structure representing a facet of the (partial) weight space
 * polyhedron P={(w,a) : w \cdot y >= a \forall y \in Y} where Y is
 * the (current) set of non-dominated points. A facet (coeffs,rhs) is
 * represented by coefficients 'coeffs' and a right hand side 'rhs'
 * yielding an inequality of the form coeffs \cdot w >= rhs 
 */

#ifndef POLYSCIP_SRC_WEIGHT_SPACE_FACET_H_INCLUDED
#define POLYSCIP_SRC_WEIGHT_SPACE_FACET_H_INCLUDED

#include <vector>

namespace polyscip {

  /** @brief Facet of the (partial) weight space polyhedron. */
  class WeightSpaceFacet {

  public: 
    using ValueType = SCIP_Real;
    using CoeffContainer = std::vector<ValueType>;
    
    /** @brief Creates a facet coeffs \cdot >= rhs of the (partial)
     *  weight space polyhedron.  
     *  @param coeffs Left hand side coefficients of the facet inequality
     *  @param rhs Right hand side value of the facet inequality
     */
    WeightSpaceFacet(const CoeffContainer& coeffs, ValueType rhs);
    
    /** @brief Destructor */
    ~WeightSpaceFacet();

    /** @brief Prints facet information to standard output.
     */
    void print() const;

  private:
    CoeffContainer coeffs_;  /**< left hand side coefficients of the facet inequality */
    ValueType rhs_;          /**< right hand side value of the facet inequality */

  };

}

#endif // POLYSCIP_SRC_WEIGHT_SPACE_FACET_H_INCLUDED

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

/**@file   polyscip.h
 * @brief  PolySCIP solver class
 * @author Sebastian Schenker
 *
 * The PolySCIP solver class.
 */

#ifndef POLYSCIP_SRC_POLYSCIP_H_INCLUDED
#define POLYSCIP_SRC_POLYSCIP_H_INCLUDED

#include "scip/def.h"

#include <vector>

namespace polyscip {
  
  class Polyscip {
  public:
    using ValueType = SCIP_Real;               /**< type for computed values */
    using PointType = std::vector<ValueType>;  /**< type for points in outcome space */
    using RayType = std::vector<ValueType>;     /**< type for rays in outcome space */
    using WeightType = std::vector<ValueType>; /**< type for weights */
    
  private: 
    /** returns true if point is a new non-dominated point; otherwise false */
    bool isNewNondomPoint(const std::vector<SCIP_Real>* point, /**< potential new nondominated point */			  double comp_val = 0.0              /**< value used for checking inequality */ 
			  ) const;

  };


}

#endif // POLYSCIP_SRC_POLYSCIP_H_INCLUDED

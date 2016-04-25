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

/** @brief  PolySCIP solver class
 *
 * The PolySCIP solver class.
 */

#ifndef POLYSCIP_SRC_POLYSCIP_H_INCLUDED
#define POLYSCIP_SRC_POLYSCIP_H_INCLUDED

#include <utility> // std::pair
#include <vector>

#include "scip/def.h"

namespace polyscip {
  
  class Polyscip {
  public:
    using ValueType = SCIP_Real;               /**< type for computed values */
    using PointType = std::vector<ValueType>;  /**< type for points in outcome space; 
						  needs to support: begin(), operator[] */
    using RayType = std::vector<ValueType>;     /**< type for rays in outcome space 
						   needs to support: size() */
    using WeightType = std::vector<ValueType>; /**< type for weights
						  needs to support: at(), size() */
    using PointContainer = std::vector<std::pair<PointType,WeightType>>;
    using RayContainer = std::vector<std::pair<RayType,WeightType>>; /**< Container type for 
									  computed rays
									Needs to support: empty() */
  private: 
    PointContainer supported_nondom_points_;
    RayContainer unbounded_nondom_rays_;

  };


}

#endif // POLYSCIP_SRC_POLYSCIP_H_INCLUDED

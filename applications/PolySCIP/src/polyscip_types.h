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

/** @brief PolySCIP types
 *
 * Types used for PolySCIP.
 */

#ifndef POLYSCIP_SRC_POLYSCIP_TYPES_H_INCLUDED
#define POLYSCIP_SRC_POLYSCIP_TYPES_H_INCLUDED

#include <string>
#include <iostream>
#include <ostream>
#include <utility>
#include <vector>
#include "scip/def.h"

namespace polyscip {

    /**< type for computed values */
    using ValueType = SCIP_Real;
    /**< type for points, rays in outcome space */
    using OutcomeType = std::vector<ValueType>;
    /**< type for solutions in feasible space */
    using SolType = std::vector< std::pair<std::string, ValueType> >;
    /**< type for weights */
    using WeightType = std::vector<ValueType>;
    /**< A result comprises of a solution/ray in feasible space and corresponding
     * non-dominated point in objective space */
    using Result = std::pair<SolType, OutcomeType>;
    using ResultContainer = std::vector<Result>;
}

#endif //POLYSCIP_SRC_POLYSCIP_TYPES_H_INCLUDED
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

/**
 * @brief PolySCIP types
 * @author Sebastian Schenker
 *
 * Types used for PolySCIP solver.
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

    /**< Type for computed values */
    using ValueType = SCIP_Real;
    /**< Type for points, rays in outcome space */
    using OutcomeType = std::vector<ValueType>;
    /**< Type for solutions in feasible space */
    using SolType = std::vector< std::pair<std::string, ValueType> >;
    /**< Type for weights vectors*/
    using WeightType = std::vector<ValueType>;
    /**< A result comprises of a solution/ray in feasible space and corresponding
     * non-dominated outcome in outcome space */
    using Result = std::pair<SolType, OutcomeType>;
    /**< Container for results */
    using ResultContainer = std::vector<Result>;
}

#endif //POLYSCIP_SRC_POLYSCIP_TYPES_H_INCLUDED
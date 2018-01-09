/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*        This file is part of the program PolySCIP                          */
/*                                                                           */
/*    Copyright (C) 2012-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  PolySCIP is distributed under the terms of the ZIB Academic License.     */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with PolySCIP; see the file LICENCE.                               */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**
 * @file polyscip_types.h
 * @brief General types used for PolySCIP
 * @author Sebastian Schenker
 *
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

    using ValueType = SCIP_Real; ///< Type for computed values
    using OutcomeType = std::vector<ValueType>; ///< Type for points, rays in objective space
    using SolType = std::vector< std::pair<std::string, ValueType> >; ///< Type for solutions in feasible space
    using WeightType = std::vector<ValueType>; ///< Type for weight vectors
    using Result = std::pair<SolType, OutcomeType>; ///< A result comprises of a solution/ray in feasible space and a corresponding outcome in objective space
    using ResultContainer = std::vector<Result>; ///< Container type for results

}

#endif //POLYSCIP_SRC_POLYSCIP_TYPES_H_INCLUDED
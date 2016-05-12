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
#include <vector>
#include "scip/def.h"

namespace polyscip {

    /**< type for computed values */
    using ValueType = SCIP_Real;
    /**< type for points, rays in outcome space; needs to support: begin(), size(), operator[] */
    using OutcomeType = std::vector<ValueType>;
    /**< type for weights; needs to support: at(), size() */
    using WeightType = std::vector<ValueType>;

    /** General print function
    * @param container Container to be printed
    * @param description Corresponding description
    * @param os Output stream to print to
    */
    template <typename Container>
    void print(const Container& container,
               std::string description,
               std::ostream& os = std::cout) {
        os << description << "[ ";
        for (const auto& elem : container)
            os << elem << " ";
        os << "]\n";
    };
}

#endif //POLYSCIP_SRC_POLYSCIP_TYPES_H_INCLUDED
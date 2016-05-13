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

/** @brief  Global available functions
 *
 * Some available template functions
 *
 */

#ifndef POLYSCIP_SRC_GLOBAL_FUNCTIONS_H_INCLUDED
#define POLYSCIP_SRC_GLOBAL_FUNCTIONS_H_INCLUDED

#include <iostream>
#include <ostream>

namespace polyscip {

    namespace global {
        /** Allows to get a field of tuple via corresponding field of enum class via std::get
        */
        template<typename E>
        constexpr typename std::underlying_type<E>::type
        toField(E enumerator) noexcept {
            return static_cast<typename std::underlying_type<E>::type>(enumerator);
        };

        /** General print function
        * @param container Container to be printed
        * @param description Corresponding description
        * @param os Output stream to print to
        */
        template<typename Container>
        void print(const Container &container,
                   std::string description,
                   std::ostream &os = std::cout) {
            os << description << "[ ";
            for (const auto &elem : container)
                os << elem << " ";
            os << "]\n";
        };
    }
}
#endif //POLYSCIP_SRC_GLOBAL_FUNCTIONS_H_INCLUDED
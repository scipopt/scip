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
 * @file global_functions.h
 * @brief  Global helper functions
 * @author Sebastian Schenker
 *
 */

#ifndef POLYSCIP_SRC_GLOBAL_FUNCTIONS_H_INCLUDED
#define POLYSCIP_SRC_GLOBAL_FUNCTIONS_H_INCLUDED

#include <functional>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace polyscip {

    namespace global {

        /**
         * Based on code by Stephan T. Lavavej at http://channel9.msdn.com/Series/
         * C9-Lectures-Stephan-T-Lavavej-Core-C-/STLCCSeries6
         */
        namespace impl_own_stl {
            /**
             * Helper function
             * @tparam T
             * @tparam Args
             * @param args
             * @return
             */
            template<typename T, typename ... Args>
            std::unique_ptr<T> make_unique_helper(std::false_type, Args&&... args) {
                return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
            }

            /**
             * Helper function
             * @tparam T
             * @tparam Args
             * @param args
             * @return
             */
            template<typename T, typename ... Args>
            std::unique_ptr<T> make_unique_helper(std::true_type, Args&&... args) {
                static_assert(std::extent<T>::value == 0,
                              "make_unique<T[N]>() is forbidden, please use make_unique<T[]>(),");
                typedef typename std::remove_extent<T>::type U;
                return std::unique_ptr<T>(new U[sizeof...(Args)]{std::forward<Args>(args)...});
            }
        }

        /**
         * @details std::make_unique did not get into the C++11 standard, so we provide it ourselves
         * until installed compiler can be expected to fully support C++14
         * @tparam T
         * @tparam Args
         * @param args
         * @return
         */
        template<typename T, typename ... Args>
        std::unique_ptr<T> make_unique(Args&&... args) {
            return impl_own_stl::make_unique_helper<T>(
                    std::is_array<T>(),std::forward<Args>(args)... );
        }

        /**
         * Print function
         * @param container Container to be printed
         * @param prefix String printed before Container
         * @param suffix String printed after Container
         * @param os Output stream to write to
         * @param negate Indicates whether elements in Container shall be negated
         * @param prec Used precision for output stream
         */
        template<typename Container>
        void print(const Container &container,
                   const std::string& prefix = "",
                   const std::string& suffix = "",
                   std::ostream &os = std::cout,
                   bool negate = false,
                   int prec = 6) {
            os << std::setprecision(prec) << prefix;
            for (const auto &elem : container)
                if (negate) {
                    os << -elem << " ";
                }
                else {
                    os << elem << " ";
                }
            os << suffix;
        }

        /**
         * For conversion between two scalar numeric types where a value might be narrowed.
         * Taken from Stroustroup - "The C++ programming language" 4th edition, page 299
         * @param v Value to convert
         * @return Value in narrowed type
         */
        template<typename Target, typename Source>
        Target narrow_cast(Source v) {
            auto r = static_cast<Target>(v); // convert the value to the target type
            if (static_cast<Source>(r)!=v)
                throw std::runtime_error("narrow_cast<>() failed\n");
            return r;
        }

    }
}
#endif //POLYSCIP_SRC_GLOBAL_FUNCTIONS_H_INCLUDED
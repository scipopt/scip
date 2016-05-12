/* Based on code by Stephan T. Lavavej at http://channel9.msdn.com/Series/
 * C9-Lectures-Stephan-T-Lavavej-Core-C-/STLCCSeries6
 */

/** @brief Implementation of make_unique template
 *
 * Unfortunately, make_unique did not get into the C++11 standard,
 * so we provide it ourselves until installed compiler can be expected to fully
 * support C++14.
 */

#ifndef POLYSCIP_SRC_OWN_MAKE_UNIQUE_H_INCLUDED
#define POLYSCIP_SRC_OWN_MAKE_UNIQUE_H_INCLUDED

#include <memory>
#include <utility>
#include <type_traits>

namespace own_stl {
    namespace impl_own_stl {
    template<typename T, typename ... Args>
            std::unique_ptr<T> make_unique_helper(std::false_type, Args&&... args) {
            return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
        }

        template<typename T, typename ... Args>
        std::unique_ptr<T> make_unique_helper(std::true_type, Args&&... args) {
            static_assert(std::extent<T>::value == 0,
                    "make_unique<T[N]>() is forbidden, please use make_unique<T[]>(),");
            typedef typename std::remove_extent<T>::type U;
            return std::unique_ptr<T>(new U[sizeof...(Args)]{std::forward<Args>(args)...});
        }
    }

    template<typename T, typename ... Args>
    std::unique_ptr<T> make_unique(Args&&... args) {
        return impl_own_stl::make_unique_helper<T>(
            std::is_array<T>(),std::forward<Args>(args)... );
    }
}
#endif //POLYSCIP_SRC_OWN_MAKE_UNIQUE_H_INCLUDED
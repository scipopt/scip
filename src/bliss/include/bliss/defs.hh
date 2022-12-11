#pragma once

/*
  Copyright (c) 2003-2021 Tommi Junttila
  Released under the GNU Lesser General Public License version 3.

  This file is part of bliss.

  bliss is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, version 3 of the License.

  bliss is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with bliss.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cassert>
#include <cstdarg>

//SV https://en.cppreference.com/w/cpp/language/operator_alternative
#include <ciso646>

/** \file
 * \brief Some common definitions.
 */

#define BLISS_VERSION "0.77"
#define BLISS_VERSION_MAJOR 0
#define BLISS_VERSION_MINOR 77

namespace bliss {

/** \brief The version number of bliss. */
static const char * const version = "0.77";



#if defined(BLISS_DEBUG)
#define BLISS_CONSISTENCY_CHECKS
#define BLISS_EXPENSIVE_CONSISTENCY_CHECKS
#endif


#if defined(BLISS_CONSISTENCY_CHECKS)
/* Force a check that the found automorphisms are valid */
#define BLISS_VERIFY_AUTOMORPHISMS
#endif


#if defined(BLISS_CONSISTENCY_CHECKS)
/* Force a check that the generated partitions are equitable */
#define BLISS_VERIFY_EQUITABLEDNESS
#endif

} // namespace bliss


/*! \mainpage Outline
 *
 * This is the C++ API documentation of bliss,
 * produced by running <a href="http://www.doxygen.org">doxygen</a> in
 * the source directory.
 *
 * The algorithms and data structures used in bliss,
 * the graph file format, as well as the compilation process
 * can be found at the
 * <a href="https://users.aalto.fi/tjunttil/bliss">bliss web site</a>.
 *
 * The C++ language API is the main API to bliss.
 * It basically consists of the public methods in the classes
 * * bliss::Graph and
 * * bliss::Digraph.
 *
 * For an example of its use,
 * see the \ref executable "source of the bliss executable".
 *
 * \section capi_sec The C language API
 *
 * The C language API is given in the file bliss_C.h.
 * It is currently only a subset of the C++ API,
 * so consider using the C++ API whenever possible.
 */

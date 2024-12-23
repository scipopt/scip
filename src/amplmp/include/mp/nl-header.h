/**
 Common definitions for NL reader and NL writer.

 Copyright (C) 2014 - 2023 AMPL Optimization Inc.

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Authors: Victor Zverovich, Gleb Belov
 */

#ifndef NLHEADER_H
#define NLHEADER_H

#include <cstddef>
#include <algorithm>
#include <array>
#include <string>
#include <cstring>

#include "mp/nl-header-c.h"

/// Namespace mp.
namespace mp {


/** Information about an optimization problem.
 *  Full documentation on the NL format:
 *  technical report "Writing .nl Files"
 *  (https://ampl.github.io/nlwrite.pdf.)
 */
struct NLProblemInfo : NLProblemInfo_C {
  /// Construct default
  NLProblemInfo() {
    std::memset(this, 0, sizeof(*this));
  }

  /** Returns the number of integer variables (includes binary variables). */
	int num_integer_vars() const {
		return num_linear_binary_vars + num_linear_integer_vars +
				num_nl_integer_vars_in_both + num_nl_integer_vars_in_cons +
				num_nl_integer_vars_in_objs;
	}

	/** Returns the number of continuous variables. */
	int num_continuous_vars() const { return num_vars - num_integer_vars(); }

	/** Returns the total number of common expressions. */
	int num_common_exprs() const {
		return num_common_exprs_in_both + num_common_exprs_in_cons +
				num_common_exprs_in_objs + num_common_exprs_in_single_cons +
				num_common_exprs_in_single_objs;
	}
};


enum {
	/** Maximum number of options reserved for AMPL use in NL and SOL formats. */
	MAX_AMPL_OPTIONS = 9
};

enum {
  /// internal constant
	VBTOL_OPTION_INDEX = 1,
  /// internal constant
	USE_VBTOL_FLAG     = 3
};

/// Namespace arith
namespace arith {

/** Floating-point arithmetic kind. */
enum Kind {

	/** Unknown floating-point arithmetic. */
  UNKNOWN = NL_ARITH_UNKNOWN,

	/**
		\rst
		Standard `IEEE-754 floating point
		<http://en.wikipedia.org/wiki/IEEE_floating_point>`_ - little endian.
		\endrst
	 */
  IEEE_LITTLE_ENDIAN = NL_ARITH_IEEE_LITTLE_ENDIAN,

	/** Standard IEEE-754 floating point - big endian. */
  IEEE_BIG_ENDIAN = NL_ARITH_IEEE_BIG_ENDIAN,

	/**
		\rst
		`IBM floating point
		<http://en.wikipedia.org/wiki/IBM_Floating_Point_Architecture>`_.
		\endrst
	 */
  IBM = NL_ARITH_IBM,

	/** VAX floating point (legacy). */
  VAX = NL_ARITH_VAX,

	/** Cray floating point. */
  CRAY = NL_ARITH_CRAY,

	/** Last floating point. */
  LAST = NL_ARITH_LAST
};

/// Returns floating-point arithmetic kind used on the current system.
Kind GetKind();

/// IsIEEE
inline bool IsIEEE(arith::Kind k) {
	return k == IEEE_LITTLE_ENDIAN || k == IEEE_BIG_ENDIAN;
}
}  // namespace arith


/// NL info
struct NLInfo : NLInfo_C {

  /** Input/output format */
  using Format = NLFormat;

  enum {
    /**
      Text format. The text format is fully portable meaning that an .nl file
      can be written on a machine of one architecture and then read on a
      machine of a different architecture.
     */
    TEXT = NL_FORMAT_TEXT,

    /**
      Binary format. The binary format is not generally portable and should
      normally be used on a single machine.
     */
    BINARY = NL_FORMAT_BINARY
  };

  /// Construct default
  NLInfo() {
    format = BINARY;
    prob_name = "nl_instance";
    num_ampl_options = 3;
    ampl_vbtol = 0.0;
    arith_kind = NL_ARITH_IEEE_LITTLE_ENDIAN;
    flags = WANT_OUTPUT_SUFFIXES;
    std::fill(ampl_options, ampl_options + MAX_AMPL_OPTIONS, 0);
    std::array<long, 3> opt_default {1, 1, 0};
    std::copy(opt_default.begin(), opt_default.end(), ampl_options);
  }
};

/**
	\rst
	An NL `header <http://en.wikipedia.org/wiki/Header_(computing)>`_
	which contains information about problem dimensions, such as the number of
	variables and constraints, and the input format.

	Base class: `mp::ProblemInfo`
	\endrst
 */
struct NLHeader : NLProblemInfo, NLInfo {
};

}  // namespace mp

#endif // NLHEADER_H

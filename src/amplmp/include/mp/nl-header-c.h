/**
 * Base C structs for NL header.
 *
 * Used both in C and C++.
 *
 */
#ifndef NLHEADERC_H
#define NLHEADERC_H

#include <stddef.h>

#ifdef __cplusplus  // Can be used from C++
extern "C" {
#endif


/** Information about an optimization problem.
 *  Full documentation on the NL format:
 *  technical report "Writing .nl Files"
 *  (https://ampl.github.io/nlwrite.pdf.)
 */
typedef struct NLProblemInfo_C {
  /** Total number of variables. */
  int num_vars;

  /**
    Number of algebraic constraints including ranges and equality constraints.
    It doesn't include logical constraints.
   */
  int num_algebraic_cons;

  /** Total number of objectives. */
  int num_objs;

  /** Number of ranges (constraints with -Infinity < LHS < RHS < Infinity). */
  int num_ranges;

  /**
    Number of equality constraints or -1 if unknown (AMPL prior to 19970627).
   */
  int num_eqns;

  /** Number of logical constraints. */
  int num_logical_cons;

  /** Number of random variables. */
  int num_rand_vars;

  /** Number of random defined variables. */
  int num_rand_common_exprs;

  /** Number of random constraints. */
  int num_rand_cons;

  /** Number of random objectives. */
  int num_rand_objs;

  /** Number of random calls. */
  int num_rand_calls;

  /** Number of stages. */
  int num_stages;

  // Nonlinear and complementarity information
  // -----------------------------------------

  /** Total number of nonlinear constraints. */
  int num_nl_cons;

  /** Total number of nonlinear objectives. */
  int num_nl_objs;

  /** Total number of complementarity conditions. */
  int num_compl_conds;

  /** Number of nonlinear complementarity conditions. */
  int num_nl_compl_conds;

  /** Number of complementarities involving double inequalities. */
  int num_compl_dbl_ineqs;

  /** Number of complemented variables with a nonzero lower bound. */
  int num_compl_vars_with_nz_lb;

  // Information about network constraints
  // -------------------------------------

  /** Number of nonlinear network constraints. */
  int num_nl_net_cons;

  /** Number of linear network constraints. */
  int num_linear_net_cons;

  // Information about nonlinear variables
  // -------------------------------------

  /**
    Number of nonlinear variables in constraints including nonlinear
    variables in both constraints and objectives.
   */
  int num_nl_vars_in_cons;

  /**
    Number of nonlinear variables in objectives including nonlinear
    variables in both constraints and objectives.
   */
  int num_nl_vars_in_objs;

  /** Number of nonlinear variables in both constraints and objectives. */
  int num_nl_vars_in_both;

  // Miscellaneous
  // -------------

  /** Number of linear network variables (arcs). */
  int num_linear_net_vars;

  /** Number of functions. */
  int num_funcs;

  // Information about discrete variables
  // ------------------------------------

  /** Number of linear binary variables. */
  int num_linear_binary_vars;

  /** Number of linear non-binary integer variables. */
  int num_linear_integer_vars;

  /**
    Number of integer nonlinear variables in both constraints and objectives.
   */
  int num_nl_integer_vars_in_both;

  /** Number of integer nonlinear variables just in constraints. */
  int num_nl_integer_vars_in_cons;

  /** Number of integer nonlinear variables just in objectives. */
  int num_nl_integer_vars_in_objs;

  // Information about nonzeros
  // --------------------------

  /** Number of nonzeros in constraints' Jacobian. */
  size_t num_con_nonzeros;

  /** Number of nonzeros in all objective gradients. */
  size_t num_obj_nonzeros;

  // Information about names.
  // Does not have to be filled for NLWriter2.
  // -----------------------

  /** Length of longest constraint or objective name
   *  if names are available. */
  int max_con_name_len;

  /** Length of longest variable name if names are available. */
  int max_var_name_len;

  // Information about common expressions
  // ------------------------------------

  /**
    Number of common expressions that appear both in constraints
    and objectives.
   */
  int num_common_exprs_in_both;

  /**
    Number of common expressions that appear in multiple constraints
    and don't appear in objectives.
   */
  int num_common_exprs_in_cons;

  /**
    Number of common expressions that appear in multiple objectives
    and don't appear in constraints.
   */
  int num_common_exprs_in_objs;

  /**
    Number of common expressions that only appear in a single constraint
    and don't appear in objectives.
   */
  int num_common_exprs_in_single_cons;

  /**
    Number of common expressions that only appear in a single objective
    and don't appear in constraints.
   */
  int num_common_exprs_in_single_objs;

} NLProblemInfo_C;


enum {
  /** Maximum number of options reserved for AMPL use in NL and SOL formats. */
  MAX_AMPL_OPTIONS = 9
};


enum {
  VBTOL_OPTION_INDEX = 1,
  USE_VBTOL_FLAG     = 3
};


/** Floating-point arithmetic kind. */
typedef enum NLArithKind {

  /** Unknown floating-point arithmetic. */
  NL_ARITH_UNKNOWN = 0,

  /**
    \rst
    Standard `IEEE-754 floating point
    <http://en.wikipedia.org/wiki/IEEE_floating_point>`_ - little endian.
    \endrst
   */
  NL_ARITH_IEEE_LITTLE_ENDIAN = 1,

  /** Standard IEEE-754 floating point - big endian. */
  NL_ARITH_IEEE_BIG_ENDIAN = 2,

  /**
    \rst
    `IBM floating point
    <http://en.wikipedia.org/wiki/IBM_Floating_Point_Architecture>`_.
    \endrst
   */
  NL_ARITH_IBM = 3,

  /** VAX floating point (legacy). */
  NL_ARITH_VAX = 4,

  /** Cray floating point. */
  NL_ARITH_CRAY = 5,

  /** Last floating point. */
  NL_ARITH_LAST = NL_ARITH_CRAY
} NLArithKind;


/** Input/output format */
typedef int NLFormat;

enum {
  /**
    Text format. The text format is fully portable meaning that an .nl file
    can be written on a machine of one architecture and then read on a
    machine of a different architecture.
   */
  NL_FORMAT_TEXT = 0,

  /**
    Binary format. The binary format is not generally portable and should
    normally be used on a single machine.
   */
  NL_FORMAT_BINARY = 1
};


/** NL Flags. */
enum {
  /** Flag that specifies whether to write output suffixes to a .sol file. */
  WANT_OUTPUT_SUFFIXES = 1
};


/**
 * NL technical info
 */
typedef struct NLInfo_C {

  /** Input/output format. */
  NLFormat format;

  /** The number of options reserved for AMPL use. */
  int num_ampl_options;

  /**
    Values of options reserved for AMPL use.
    Leave the default values if not using AMPL.
   */
  long ampl_options[MAX_AMPL_OPTIONS];

  /**
    Extra info for writing a solution reserved for AMPL use.
    Leave the default value if not using AMPL.
   */
  double ampl_vbtol;

  /**
   * Problem name.
   */
  const char* prob_name;

  /**
    \rst
    Floating-point arithmetic kind used with binary format to check
    if an .nl file is written using a compatible representation of
    floating-point numbers. It is not used with the text format and normally
    set to `NL_ARITH_UNKNOWN` there.
    \endrst
   */
  int arith_kind;

  /**
    \rst
    Flags. Can be either 0 or `mp::NLHeader::WANT_OUTPUT_SUFFIXES`.
    \endrst
   */
  int flags;
} NLInfo_C;


/**
  \rst
  An NL `header <http://en.wikipedia.org/wiki/Header_(computing)>`_
  which contains information about problem dimensions, such as the number of
  variables and constraints, and the input format.
  \endrst
 */
typedef struct NLHeader_C {
  NLProblemInfo_C pi;
  NLInfo_C nli;
} NLHeader_C;


/// Default NLHeader_C
NLHeader_C MakeNLHeader_C_Default(void);


#ifdef __cplusplus
}  // extern "C"
#endif

#endif // NLHEADERC_H

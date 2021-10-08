/*
 Common definitions

 Copyright (C) 2014 - 2016 AMPL Optimization Inc

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

 Author: Victor Zverovich
 */

#ifndef MP_COMMON_H_
#define MP_COMMON_H_

#include "mp/error.h"  // for MP_ASSERT
#include <cstddef>     // for std::size_t

/** The mp namespace. */
namespace mp {

/** Variable information */
namespace var {
/** Variable type. */
enum Type {
  CONTINUOUS, /**< A continuous variable. */
  INTEGER     /**< An integer variable. */
};
}

/** Objective information. */
namespace obj {
/** Objective type. */
enum Type {
  MIN = 0, /**< A minimization objective. */
  MAX = 1  /**< A maximization objective. */
};
}

/** Function information. */
namespace func {
/** Function type. */
enum Type {
  /** A numeric function. */
  NUMERIC  = 0,
  /** A symbolic function accepting numeric and string arguments. */
  SYMBOLIC = 1
};
}

/** Complementarity constraint information. */
class ComplInfo {
 private:
  int flags_;

 public:
  // Flags for the constructor.
  enum {
    /** Constraint upper bound is  infinity (finite variable lower bound). */
    INF_UB = 1,
    /** Constraint lower bound is -infinity (finite variable upper bound). */
    INF_LB = 2
  };

  /**
    \rst
    Constructs a `ComplInfo` object.
    \endrst
   */
  explicit ComplInfo(int flags) : flags_(flags) {
    MP_ASSERT((flags & ~(INF_UB | INF_LB)) == 0,
              "invalid complementarity flags");
  }

  /** Constraint lower bound. */
  double con_lb() const {
    return (flags_ & INF_LB) != 0 ?
          -std::numeric_limits<double>::infinity() : 0;
  }

  /** Constraint upper bound. */
  double con_ub() const {
    return (flags_ & INF_UB) != 0 ? std::numeric_limits<double>::infinity() : 0;
  }
};

/** Suffix information. */
namespace suf {
/** Suffix kind. */
enum Kind {
  VAR     =    0,  /**< Applies to variables. */
  CON     =    1,  /**< Applies to constraints. */
  OBJ     =    2,  /**< Applies to objectives. */
  PROBLEM =    3   /**< Applies to problems. */
};

// SV disabled following line as not compilable with ancient MSVC 10.0
// constexpr int KIND_MASK = Kind::VAR | Kind::CON | Kind::OBJ | Kind ::PROBLEM;
// Suffix flags.
enum {
  FLOAT   =    4,  /**< Suffix values are floating-point numbers. */
  IODECL  =    8,  /**< Declare an INOUT suffix. */
  OUTPUT  = 0x10,  /**< Output suffix: return values from a solver. */
  INPUT   = 0x20,  /**< Input suffix: values were passed to a solver. */
  OUTONLY = 0x40   /**< Output only: reject as an input value. */
};
}

namespace internal {
enum {
  SUFFIX_KIND_MASK = 3,  // Mask for suffix kinds.
  NUM_SUFFIX_KINDS = 4   // The number of suffix kinds.
};
}

namespace sol {
/** Solution status. */
enum Status {
  NOT_CHECKED = -200,
  /** Unknown status. */
  UNKNOWN     =  -1,

  /**
    An optimal solution found for an optimization problem or a feasible
    solution found for a satisfaction problem.
   */
  SOLVED      =   0,

  /** Solution returned but it can be non-optimal or even infeasible. */
  UNCERTAIN   = 100,

  /** Problem is infeasible. */
  INFEASIBLE  = 200,

  /** Problem is unbounded. */
  UNBOUNDED   = 300,

  /** Stopped by a limit, e.g. on iterations or time. */
  LIMIT       = 400,

  /** A solver error. */
  FAILURE     = 500,

  /** Interrupted by the user. */
  INTERRUPTED = 600
};
}

/** Expression information. */
namespace expr {

/**
  \rst
  Expression kind.
  Both AMPL-like and mathematical notation is given for each expression in the
  descriptions below as in :math:`\mathrm{abs}(x) = |x|`, unless they are
  identical such as :math:`\mathrm{sin}(x)` or there is no standard
  mathematical notation.
  \endrst
  */
enum Kind {
  /** An unknown expression. */
  UNKNOWN = 0,

  /** The first expression kind other than the unknown expression kind. */
  FIRST_EXPR,

  /**
    \rst
    The first numeric expression kind. Numeric expression kinds are in
    the range [`~mp::expr::FIRST_NUMERIC`, `~mp::expr::LAST_NUMERIC`].
    \endrst
   */
  FIRST_NUMERIC = FIRST_EXPR,

  /**
    \rst
    A number such as 42 or -1.23e-4.
    \endrst
   */
  NUMBER = FIRST_NUMERIC,

  /**
    \rst
    The first reference expression kind. Reference expression kinds are in
    the range [`~mp::expr::FIRST_REFERENCE`, `~mp::expr::LAST_REFERENCE`].
    \endrst
   */
  FIRST_REFERENCE,

  /** A reference to a variable. */
  VARIABLE = FIRST_REFERENCE,

  /** A reference to a common expression. */
  COMMON_EXPR,

  /** The last reference expression kind. */
  LAST_REFERENCE = COMMON_EXPR,

  /**
    \rst
    The first unary numeric expression kind. Unary numeric expression kinds
    are in the range [`~mp::expr::FIRST_UNARY`, `~mp::expr::LAST_UNARY`].
    \endrst
   */
  FIRST_UNARY,

  /**
    \rst
    A unary minus, :math:`-x`.
    \endrst
   */
  MINUS = FIRST_UNARY,

  /**
    \rst
    The absolute value function, :math:`\mathrm{abs}(x) = |x|`.
    \endrst
   */
  ABS,

  /**
    \rst
    The floor function, :math:`\mathrm{floor}(x) = \lfloor x \rfloor`.
    \endrst
   */
  FLOOR,

  /**
    \rst
    The ceiling function, :math:`\mathrm{ceil}(x) = \lceil x \rceil`.
    \endrst
   */
  CEIL,

  /**
    \rst
    The square root function, :math:`\mathrm{sqrt}(x) = \sqrt{x}`.
    \endrst
   */
  SQRT,

  /**
    \rst
    Squaring: :math:`x \mathop{\verb!^!} 2 = x^2`.
    \endrst
   */
  POW2,

  /**
    \rst
    The natural exponential function, :math:`\mathrm{exp}(x) = e^x`.
    \endrst
   */
  EXP,

  /**
    \rst
    The natural logarithmic function, :math:`\mathrm{log}(x) = \mathrm{ln}(x)`.
    \endrst
   */
  LOG,

  /**
    \rst
    The base 10 logarithmic function,
    :math:`\mathrm{log10}(x) = \mathrm{log}_{10}(x)`.
    \endrst
   */
  LOG10,

  /**
    \rst
    Sine, :math:`\mathrm{sin}(x)`.
    \endrst
   */
  SIN,

  /**
    \rst
    Hyperbolic sine, :math:`\mathrm{sinh}(x)`.
    \endrst
   */
  SINH,

  /**
    \rst
    Cosine, :math:`\mathrm{cos}(x)`.
    \endrst
   */
  COS,

  /**
    \rst
    Hyperbolic cosine, :math:`\mathrm{cosh}(x)`.
    \endrst
   */
  COSH,

  /**
    \rst
    Tangent, :math:`\mathrm{tan}(x)`.
    \endrst
   */
  TAN,

  /**
    \rst
    Hyperbolic tangent, :math:`\mathrm{tan}(x)`.
    \endrst
   */
  TANH,

  /**
    \rst
    Inverse sine, :math:`\mathrm{asin}(x) = \mathrm{sin}^{-1}(x)`.
    \endrst
   */
  ASIN,

  /**
    \rst
    Inverse hyperbolic sine, :math:`\mathrm{asinh}(x) = \mathrm{sinh}^{-1}(x)`.
    \endrst
   */
  ASINH,

  /**
    \rst
    Inverse cosine, :math:`\mathrm{acos}(x) = \mathrm{cos}^{-1}(x)`.
    \endrst
   */
  ACOS,

  /**
    \rst
    Inverse hyperbolic cosine,
    :math:`\mathrm{acosh}(x) = \mathrm{cosh}^{-1}(x)`.
    \endrst
   */
  ACOSH,

  /**
    \rst
    Inverse tangent, :math:`\mathrm{atan}(x) = \mathrm{tan}^{-1}(x)`.
    \endrst
   */
  ATAN,

  /**
    \rst
    Inverse hyperbolic tangent,
    :math:`\mathrm{atanh}(x) = \mathrm{tanh}^{-1}(x)`.
    \endrst
   */
  ATANH,


  /** The last unary numeric expression kind. */
  LAST_UNARY = ATANH,

  /**
    \rst
    The first binary expression kind. Binary expression kinds are in
    the range [`~mp::expr::FIRST_BINARY`, `~mp::expr::LAST_BINARY`].
    \endrst
   */
  FIRST_BINARY,

  /**
    \rst
    Addition, :math:`x + y`.
    \endrst
   */
  ADD = FIRST_BINARY,

  /**
    \rst
    Subtraction, :math:`x - y`.
    \endrst
   */
  SUB,

  /**
    \rst
    The :math:`\mathrm{less}` operation,
    :math:`x \mathop{\rm less} y = \mathrm{max}(x - y, 0)`.
    \endrst
   */
  LESS,

  /**
    \rst
    Multiplication, :math:`x * y = x y`.
    \endrst
   */
  MUL,

  /**
    \rst
    Division, :math:`x / y`.
    \endrst
   */
  DIV,

  /**
    \rst
    Truncated division, :math:`x \mathop{\rm div} y = \mathrm{trunc}(x / y)`.
    \endrst
   */
  TRUNC_DIV,

  /**
    \rst
    The modulo operation, :math:`x \mathop{\rm mod} y`.
    \endrst
   */
  MOD,

  /**
    \rst
    Exponentiation, :math:`x \mathop{\verb!^!} y = x^y`.
    \endrst
   */
  POW,

  /**
    \rst
    Exponentiation with a constant base, :math:`a^x`.
    \endrst
   */
  POW_CONST_BASE,

  /**
    \rst
    Exponentiation with a constant exponent :math:`x^a`.
    \endrst
   */
  POW_CONST_EXP,

  /**
    \rst
    Inverse tangent, :math:`\mathrm{atan2}(y, x) = \mathrm{tan}^{-1}(y/x)`.
    \endrst
   */
  ATAN2,

  /**
    \rst
    The function :math:`\mathrm{precision}(x, n)` which returns :math:`x`
    rounded to :math:`n` significant decimal digits.
    \endrst
   */
  PRECISION,

  /**
    \rst
    The function :math:`\mathrm{round}(x, n)` which returns :math:`x`
    rounded to :math:`n` digits past decimal point.
    \endrst
   */
  ROUND,

  /**
    \rst
    The function :math:`\mathrm{trunc}(x, n)` which returns :math:`x`
    truncated to :math:`n` digits past decimal point.
    \endrst
   */
  TRUNC,

  /** The last binary numeric expression kind. */
  LAST_BINARY = TRUNC,

  /**
    \rst
    An if-then-else expression,
    :math:`\mathrm{if}\;c\;\mathrm{then}\;e_1\;[\mathrm{else}\;e_2]`,
    where :math:`c` is a logical expression representing condition, while
    :math:`e_1` and :math:`e_2` are numeric expressions. The expression
    evaluates to :math:`e_1` if :math:`c` is true and to :math:`e_2` otherwise.
    If the else clause is omitted, :math:`e_2` is assumed to be zero.
    \endrst
   */
  IF,

  /**
    \rst
    A piecewise-linear term,
    :math:`\verb|<<|b_1, ..., b_n; s_1, ..., s_{n + 1}\verb|>> | r`,
    where :math:`b_i` are breakpoints, :math:`s_i` are slopes and :math:`r` is
    a `reference <mp::expr::FIRST_REFERENCE>`.
    \endrst
   */
  PLTERM,

  /**
    \rst
    A function call, :math:`f(e_1, ..., e_n)`, where :math:`f` is a function
    name and :math:`e_i` are numeric or string expressions.
    \endrst
   */
  CALL,

  /**
    \rst
    The first iterated expression kind. Iterated expression kinds are in
    the range [`~mp::expr::FIRST_ITERATED`, `~mp::expr::LAST_ITERATED`].

    The term "iterated" in the context of operators and expressions comes
    from the article `AMPL: A Mathematical Programming Language
    <http://www.ampl.com/REFS/amplmod.pdf>`_ and is used to denote operators
    indexed over sets.
    \endrst
   */
  FIRST_ITERATED,

  /**
    \rst
    A vararg expression, :math:`\mathrm{min}` or :math:`\mathrm{max}`.
    Vararg expression kinds are in the range
    [`~mp::expr::FIRST_VARARG`, `~mp::expr::LAST_VARARG`].
    \endrst
   */
  FIRST_VARARG = FIRST_ITERATED,

  /**
    \rst
    Minimum, :math:`\mathrm{min}(e_1, ..., e_n) = \min_{i=1,...,n} e_i`.
    \endrst
   */
  MIN = FIRST_VARARG,

  /**
    \rst
    Maximum, :math:`\mathrm{max}(e_1, ..., e_n) = \max_{i=1,...,n} e_i`.
    \endrst
   */
  MAX,

  /** The last vararg expression kind. */
  LAST_VARARG = MAX,

  /**
    \rst
    Summation, :math:`\mathrm{sum}(e_1, ..., e_n) = \sum_{i=1}^n e_i`.
    \endrst
   */
  SUM,

  /**
    \rst
    A :math:`\mathrm{numberof}` expression,
    :math:`\mathrm{numberof}\;e_0\;\mathrm{in}\;(e_1, ..., e_n)`, which
    evaluates to the number of times the value of :math:`e_0` appears among the
    values of :math:`e_1, ..., e_n`.
    \endrst
   */
  NUMBEROF,

  /** The last iterated expression kind. */
  LAST_ITERATED = NUMBEROF,

  /**
    \rst
    A symbolic :math:`\mathrm{numberof}` expression.
    :math:`\mathrm{numberof}\;s_0\;\mathrm{in}\;(s_1, ..., s_n)`, which
    evaluates to the number of times the value of :math:`s_0` appears among the
    values of :math:`s_1, ..., s_n`.
    \endrst
   */
  NUMBEROF_SYM,

  /**
    \rst
    A :math:`\mathrm{count}` expression, :math:`\mathrm{count}(l_1, ..., l_n)`,
    where :math:`l_i` are logical expressions. This expression evaluates to
    the number of :math:`l_i` whose values are true.
    \endrst
   */
  COUNT,

  /** The last numeric expression kind. */
  LAST_NUMERIC = COUNT,

  /**
    \rst
    The first logical expression kind. Logical expression kinds are in
    the range [`~mp::expr::FIRST_LOGICAL`, `~mp::expr::LAST_LOGICAL`].
    \endrst
   */
  FIRST_LOGICAL,

  /**
    \rst
    A Boolean (logical) constant, true or false.
    \endrst
   */
  BOOL = FIRST_LOGICAL,

  /**
    \rst
    A logical not, :math:`!l`, where :math:`l` is a logical expression.
    \endrst
   */
  NOT,

  /**
    \rst
    The first binary logical expression kind.
    Binary logical expression kinds are in the range
    [`~mp::expr::FIRST_BINARY_LOGICAL`, `~mp::expr::LAST_BINARY_LOGICAL`].
    \endrst
   */
  FIRST_BINARY_LOGICAL,

  /**
    \rst
    Logical or, :math:`l_1` || :math:`l_2`.
    \endrst
   */
  OR = FIRST_BINARY_LOGICAL,

  /**
    \rst
    Logical and, :math:`l_1` && :math:`l_2`.
    \endrst
   */
  AND,

  /**
    \rst
    If and only if, :math:`l_1` <==> :math:`l_2`.
    \endrst
   */
  IFF,

  /** The last binary logical expression kind. */
  LAST_BINARY_LOGICAL = IFF,

  /**
    \rst
    The first relational expression kind. Relational expression kinds are in
    the range [`~mp::expr::FIRST_RELATIONAL`, `~mp::expr::LAST_RELATIONAL`].
    \endrst
   */
  FIRST_RELATIONAL,

  /**
    \rst
    Less than, :math:`e_1` < :math:`e_2`.
    \endrst
   */
  LT = FIRST_RELATIONAL,

  /**
    \rst
    Less or equal to, :math:`e_1` <= :math:`e_2`.
    \endrst
   */
  LE,

  /**
    \rst
    Equal to, :math:`e_1` = :math:`e_2`.
    \endrst
   */
  EQ,

  /**
    \rst
    Greater or equal to, :math:`e_1` >= :math:`e_2`.
    \endrst
   */
  GE,

  /**
    \rst
    Greater than, :math:`e_1` > :math:`e_2`.
    \endrst
   */
  GT,

  /**
    \rst
    Not equal to, :math:`e_1` != :math:`e_2`.
    \endrst
   */
  NE,

  /** The last relational expression kind. */
  LAST_RELATIONAL = NE,

  /**
    \rst
    The first logical count expression kind.
    Logical count expression kinds are in the range
    [`~mp::expr::FIRST_LOGICAL_COUNT`, `~mp::expr::LAST_LOGICAL_COUNT`].
    \endrst
   */
  FIRST_LOGICAL_COUNT,

  /**
    \rst
    An :math:`\mathrm{atleast}` expression,
    :math:`\mathrm{atleast}\;e\;(l_1, ..., l_n)`, where :math:`e` is a numeric
    expression and :math:`l_i` are logical expressions. It evaluates to true if
    at least :math:`e` expressions :math:`l_i` are true.
    \endrst
   */
  ATLEAST = FIRST_LOGICAL_COUNT,

  /**
    \rst
    An :math:`\mathrm{atmost}` expression,
    :math:`\mathrm{atmost}\;e\;(l_1, ..., l_n)`, where :math:`e` is a numeric
    expression and :math:`l_i` are logical expressions. It evaluates to true if
    at most :math:`e` expressions :math:`l_i` are true.
    \endrst
   */
  ATMOST,

  /**
    \rst
    An :math:`\mathrm{exactly}` expression,
    :math:`\mathrm{exactly}\;e\;(l_1, ..., l_n)`, where :math:`e` is a numeric
    expression and :math:`l_i` are logical expressions. It evaluates to true if
    exactly :math:`e` expressions :math:`l_i` are true.
    \endrst
   */
  EXACTLY,

  /**
    \rst
    The negation of an :math:`\mathrm{atleast}` expression,
    :math:`!\mathrm{atleast}\;e\;(l_1, ..., l_n)`.
    \endrst
   */
  NOT_ATLEAST,

  /**
    \rst
    The negation of an :math:`\mathrm{atmost}` expression,
    :math:`!\mathrm{atmost}\;e\;(l_1, ..., l_n)`.
    \endrst
   */
  NOT_ATMOST,

  /**
    \rst
    The negation of an :math:`\mathrm{exactly}` expression,
    :math:`!\mathrm{exactly}\;e\;(l_1, ..., l_n)`.
    \endrst
   */
  NOT_EXACTLY,

  /** The last logical count expression kind. */
  LAST_LOGICAL_COUNT = NOT_EXACTLY,

  /**
    \rst
    An implication expression,
    :math:`c\;\verb|==>|\;l_1\;[\mathrm{else}\;l_2]`,
    where :math:`c` is a logical expression representing condition, while
    :math:`l_1` and :math:`l_2` are logical expressions. The expression
    evaluates to :math:`l_1` if :math:`c` is true and to :math:`l_2` otherwise.
    If the else clause is omitted, :math:`l_2` is assumed to be true.
    \endrst
   */
  IMPLICATION,

  /**
    \rst
    The first iterated logical expression kind.
    Iterated logical expression kinds are in the range
    [`~mp::expr::FIRST_ITERATED_LOGICAL`, `~mp::expr::LAST_ITERATED_LOGICAL`].
    \endrst
   */
  FIRST_ITERATED_LOGICAL,

  /**
    \rst
    An :math:`\mathrm{exists}` expression,
    :math:`\mathrm{exists}(l_1, ..., l_n)`, where :math:`l_i` are logical
    expressions. It evaluates to true if at least one :math:`l_i` is true.
    \endrst
   */
  EXISTS = FIRST_ITERATED_LOGICAL,

  /**
    \rst
    A :math:`\mathrm{forall}` expression,
    :math:`\mathrm{forall}(l_1, ..., l_n)`, where :math:`l_i` are logical
    expressions. It evaluates to true if all :math:`l_i` are true.
    \endrst
   */
  FORALL,

  /** The last iterated logical expression kind. */
  LAST_ITERATED_LOGICAL = FORALL,

  /**
    \rst
    The first pairwise expression kind. Pairwise expression kinds are in the
    range [`~mp::expr::FIRST_PAIRWISE`, `~mp::expr::LAST_PAIRWISE`].
    \endrst
   */
  FIRST_PAIRWISE,

  /**
    \rst
    An alldifferent expression, :math:`\mathrm{alldiff}(e_1, ..., e_n)`,
    where :math:`e_i` are numeric expressions. It evaluates to true if all
    :math:`e_i` take different values.
    \endrst
   */
  ALLDIFF = FIRST_PAIRWISE,

  /**
    \rst
    The negation of an alldifferent expression,
    :math:`!\mathrm{alldiff}(e_1, ..., e_n)`.
    \endrst
   */
  NOT_ALLDIFF,

  /** The last pairwise expression kind. */
  LAST_PAIRWISE = NOT_ALLDIFF,

  /** The last logical expression kind. */
  LAST_LOGICAL = LAST_PAIRWISE,

  /** A string such as "abc". */
  STRING,

  /**
    \rst
    A symbolic if-then-else expression.
    :math:`\mathrm{if}\;c\;\mathrm{then}\;e_1\;[\mathrm{else}\;e_2]`,
    where :math:`c` is a logical expression representing condition, while
    :math:`e_1` and :math:`e_2` are numeric or string expressions.
    The expression evaluates to :math:`e_1` if :math:`c` is true and to
    :math:`e_2` otherwise. If :math:`e_2` is omitted, it is assumed to be zero.
    \endrst
   */
  IFSYM,

  /** The last expression kind. */
  LAST_EXPR = IFSYM
};

/**
  \rst
  Returns the string representation of the given expression kind.
  Expressions of different kinds can have identical strings.
  For example, `~mp::expr::POW`, `~mp::expr::POW_CONST_BASE` and
  `~mp::expr::POW_CONST_EXP` all have the same representation "^".
  \endrst
 */
const char *str(expr::Kind kind);

/** Returns the NL opcode for the given expression kind. */
int nl_opcode(expr::Kind kind);
}  // namespace expr

#define MP_CONST_DISPATCH(call) static_cast<const Impl*>(this)->call
#define MP_DISPATCH(call) static_cast<Impl*>(this)->call
#define MP_DISPATCH_STATIC(call) Impl::call

namespace internal {

// Suppresses warnings about unused variables.
inline void Unused(...) {}

// Returns true if ExprType is of kind k.
template <typename ExprType>
inline bool Is(expr::Kind k) {
  int kind = k;
  // If FIRST_KIND == LAST_KIND, then a decent optimizing compiler simplifies
  // this to kind == ExprType::FIRST_KIND (checked with GCC 4.8.2).
  // No need to do it ourselves.
  return ExprType::FIRST_KIND <= kind && kind <= ExprType::LAST_KIND;
}

inline bool IsValid(expr::Kind kind) {
  return kind >= expr::UNKNOWN && kind <= expr::LAST_EXPR;
}

// Expression information.
class ExprInfo {
 private:
  static const ExprInfo INFO[];

  friend int expr::nl_opcode(expr::Kind kind);
  friend const char *expr::str(expr::Kind kind);

 public:
  int opcode;
  const char *str;
};

// Maximum NL opcode.
enum { MAX_OPCODE = 82 };

class OpCodeInfo {
 private:
  static const OpCodeInfo INFO[MAX_OPCODE + 1];

 public:
  expr::Kind kind;
  expr::Kind first_kind;  // First member of a kind.

  friend const OpCodeInfo &GetOpCodeInfo(int opcode);
};

inline const OpCodeInfo &GetOpCodeInfo(int opcode) {
  MP_ASSERT(opcode >= 0 && opcode <= MAX_OPCODE, "invalid opcode");
  return OpCodeInfo::INFO[opcode];
}
}  // namespace internal

inline int expr::nl_opcode(expr::Kind kind) {
  MP_ASSERT(internal::IsValid(kind), "invalid expression kind");
  return internal::ExprInfo::INFO[kind].opcode;
}

inline const char *expr::str(expr::Kind kind) {
  MP_ASSERT(internal::IsValid(kind), "invalid expression kind");
  return internal::ExprInfo::INFO[kind].str;
}

/** Information about an optimization problem. */
struct ProblemInfo {
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

  /** Returns the number of integer variables (includes binary variable). */
  int num_integer_vars() const {
    return num_linear_binary_vars + num_linear_integer_vars +
        num_nl_integer_vars_in_both + num_nl_integer_vars_in_cons +
        num_nl_integer_vars_in_objs;
  }

  /** Returns the number of continuous variables. */
  int num_continuous_vars() const { return num_vars - num_integer_vars(); }

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
  std::size_t num_con_nonzeros;

  /** Number of nonzeros in all objective gradients. */
  std::size_t num_obj_nonzeros;

  // Information about names
  // -----------------------

  /** Length of longest constraint name if names are available. */
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
    Number of common expressions that appear in multiple objective
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

  /** Returns the total number of common expressions. */
  int num_common_exprs() const {
    return num_common_exprs_in_both + num_common_exprs_in_cons +
        num_common_exprs_in_objs + num_common_exprs_in_single_cons +
        num_common_exprs_in_single_objs;
  }
};
/* SV disabled enum classes as not compilable with ancient MSVC 10.0
enum class IISStatus {
  non = 0,
  low = 1,
  fix = 2,
  upp = 3,
  mem = 4,
  pmem = 5,
  plow = 6,
  pupp = 7,
  bug = 8
};

enum class BasicStatus {
  none= 0,
  bas = 1,
  sup = 2,
  low = 3,
  upp = 4,
  equ = 5,
  btw = 6
};
*/
}  // namespace mp

#endif  // MP_COMMON_H_

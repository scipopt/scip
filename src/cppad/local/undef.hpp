/* $Id: undef.hpp 2212 2011-11-25 06:22:10Z bradbell $ */
# ifndef CPPAD_UNDEF_INCLUDED
# define CPPAD_UNDEF_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

//      CPPAD_ASSERT_FIRST_CALL_NOT_PARALLEL used by CPPAD_USER_ATOMIC
//      CPPAD_ASSERT_KNOWN needed by cppad_ipopt
# undef CPPAD_ASSERT_NARG_NRES
# undef CPPAD_ASSERT_NOT_PARALLEL
//      CPPAD_ASSERT_UNKNOWN needed by cppad_ipopt
//      CPPAD_BEGIN_NAMESPACE needed by cppad_ipopt
//      CPPAD_BOOL_BINARY in user api
//      CPPAD_BOOL_UNARY in user api
# undef CPPAD_BOOSTVECTOR
# undef CPPAD_COND_EXP
# undef CPPAD_COND_EXP_BASE_REL
# undef CPPAD_COND_EXP_REL
# undef CPPAD_CPPADVECTOR
//      CPPAD_DISCRETE_FUNCTION
//      CPPAD_END_NAMESPACE needed by cppad_ipopt
# undef CPPAD_FOLD_AD_VALUED_BINARY_OPERATION
# undef CPPAD_FOLD_AD_VALUED_BINARY_OPERATOR
# undef CPPAD_FOLD_ASSIGNMENT_OPERATOR
# undef CPPAD_FOLD_BOOL_VALUED_BINARY_OPERATOR
# undef CPPAD_FOR_JAC_SWEEP_TRACE
# undef CPPAD_FORWARD0SWEEP
# undef CPPAD_FORWARD0SWEEP_TRACE
# undef CPPAD_FORWARD_SWEEP_TRACE
# undef CPPAD_GETTIMEOFDAY
//      CPPAD_HASH_TABLE_SIZE used by test_more/optimize.cpp
# undef CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
# undef CPPAD_MAX_NUM_CAPACITIES
# undef CPPAD_MAX_NUM_CAPACITY
//      CPPAD_MAX_NUM_THREADS in user api
# undef CPPAD_MIN_DOUBLE_CAPACITY
# undef CPPAD_NUL
//      CPPAD_NULL is in user api
# undef CPPAD_OP_CODE_TYPE
//      CPPAD_PACKAGE_STRING in user api
# undef CPPAD_REVERSE_SWEEP_TRACE
# undef CPPAD_REV_HES_SWEEP_TRACE
# undef CPPAD_REV_JAC_SWEEP_TRACE
# undef CPPAD_SIZE_T_SAME_UNSIGNED_INT
# undef CPPAD_STANDARD_MATH_UNARY
# undef CPPAD_STANDARD_MATH_UNARY_AD
# undef CPPAD_STDVECTOR
//      CPPAD_TAPE_ADDR_TYPE in user api
# undef CPPAD_TAPE_ID_TYPE 
//      CPPAD_TEST_VECTOR in user api
# undef CPPAD_TRACE_CAPACITY
# undef CPPAD_TRACE_THEAD
# undef CPPAD_TRACE_THREAD
//      CPPAD_TRACK_COUNT in deprecated api
# undef CPPAD_TRACK_DEBUG
//      CPPAD_TRACK_DEL_VEC in deprecated api
//      CPPAD_TRACK_EXTEND in deprecated api
//      CPPAD_TRACK_NEW_VEC in deprecated api
# undef CPPAD_TYPE_ADDR_TYPE
# undef CPPAD_USE_FORWARD0SWEEP
//      CPPAD_USER_ATOMIC
# undef CPPAD_USER_MACRO
# undef CPPAD_USER_MACRO_ONE
# undef CPPAD_USER_MACRO_TWO
# undef CPPAD_VEC_AD_COMPUTED_ASSIGNMENT

# endif

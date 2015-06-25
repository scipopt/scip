/* $Id: configure.hpp.in 3064 2013-12-28 18:01:30Z bradbell $ */
# ifndef CPPAD_CONFIGURE_INCLUDED
# define CPPAD_CONFIGURE_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-13 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Eclipse Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin configure$$
$spell
	CppAD
$$

$section Preprocessor Definitions Used by CppAD$$
$index configure, symbol$$
$index symbol, configure$$

$head Preprocessor Symbols$$
$index preprocessor, symbol$$
$index symbol, processor$$
All of the preprocessor symbols used by CppAD begin with
$code CPPAD_$$
(there are some deprecated symbols that begin with $code CppAD$$).

$end
*/

/*!
\defgroup configure_hpp configure.hpp
\{
 \file configure.hpp
Replacement for config.h so that all preprocessor symbols begin with CPPAD_ 
*/

/*!
\def CPPAD_PACKAGE_STRING
cppad-yyyymmdd as a C string where yyyy is year, mm is month, and dd is day.
*/
# define CPPAD_PACKAGE_STRING "cppad-20140000.3"

/*!
def CPPAD_HAS_NULLPTR
Does this compiler support the a c++11 null-pointer constant nullptr
(true = 1, false = 0).
*/
# define CPPAD_HAS_NULLPTR 0
 
 /*!
def CPPAD_INTERNAL_SPARSE_SET
is the internal representation used for sparse vectors of std::set<size_t>
either sparse_set or sparse_list).
*/
# define CPPAD_INTERNAL_SPARSE_SET sparse_list

/*!
\def CPPAD_IMPLICIT_CTOR_FROM_ANY_TYPE
If this symbol is one, an implicit constor of AD<Base> is defined
where the argument has any type.
Otherwise this constructor is explicit.
*/
# define CPPAD_IMPLICIT_CTOR_FROM_ANY_TYPE 0

/*!
\def CPPAD_BOOSTVECTOR
If this symbol is one, and _MSC_VER is not defined,
we are using boost vector for CPPAD_TESTVECTOR.
It this symbol is zero, 
we are not using boost vector for CPPAD_TESTVECTOR.
*/
# define CPPAD_BOOSTVECTOR 0

/*!
\def CPPAD_CPPADVECTOR
If this symbol is one, 
we are using CppAD vector for CPPAD_TESTVECTOR.
It this symbol is zero, 
we are not using CppAD vector for CPPAD_TESTVECTOR.
*/
# define CPPAD_CPPADVECTOR 1

/*!
\def CPPAD_STDVECTOR
If this symbol is one, 
we are using standard vector for CPPAD_TESTVECTOR.
It this symbol is zero, 
we are not using standard vector for CPPAD_TESTVECTOR.
*/
# define CPPAD_STDVECTOR 0

/*!
\def CPPAD_EIGENVECTOR
If this symbol is one, 
we are using Eigen vector for CPPAD_TESTVECTOR.
If this symbol is zero, 
we are not using Eigen vector for CPPAD_TESTVECTOR.
*/
# define CPPAD_EIGENVECTOR 0

/*!
\def CPPAD_HAS_GETTIMEOFDAY
If this symbol is one, and _MSC_VER is not defined,
this system supports the gettimeofday funcgtion.
Otherwise, this smybol should be zero.
*/
# define CPPAD_HAS_GETTIMEOFDAY 1

/*!
\def CPPAD_SIZE_T_SAME_UNSIGNED_INT 
If this symbol is one, the type size_t is the same as the type unsigned int,
otherwise this symbol is zero.
*/
# define CPPAD_SIZE_T_SAME_UNSIGNED_INT 1

/*!
\def CPPAD_TAPE_ADDR_TYPE
Is the type used to store address on the tape. If not size_t, then 
<code>sizeof(CPPAD_TAPE_ADDR_TYPE) <= sizeof( size_t )</code>
to conserve memory.
This type must support \c std::numeric_limits,
the \c <= operator,
and conversion to \c size_t.
Make sure that the type chosen returns true for is_pod<CPPAD_TAPE_ADDR_TYPE>
in pod_vector.hpp.
This type is later defined as \c addr_t in the CppAD namespace.
*/
# define CPPAD_TAPE_ADDR_TYPE size_t

/*!
\def CPPAD_TAPE_ID_TYPE
Is the type used to store tape identifiers. If not size_t, then 
<code>sizeof(CPPAD_TAPE_ID_TYPE) <= sizeof( size_t )</code>
to conserve memory.
This type must support \c std::numeric_limits,
the \c <= operator,
and conversion to \c size_t.
Make sure that the type chosen returns true for is_pod<CPPAD_TAPE_ID_TYPE>
in pod_vector.hpp.
This type is later defined as \c tape_id_t in the CppAD namespace.
*/
# define CPPAD_TAPE_ID_TYPE size_t
 
/*!
\def CPPAD_MAX_NUM_THREADS
Specifies the maximum number of threads that CppAD can support
(must be greater than or equal four).

The user may define CPPAD_MAX_NUM_THREADS before including any of the CppAD
header files.  If it is not yet defined,
*/
# ifndef CPPAD_MAX_NUM_THREADS
# define CPPAD_MAX_NUM_THREADS 48
# endif

/*! \} */
# endif

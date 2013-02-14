/* $Id: configure.hpp 2285 2012-02-14 13:13:19Z bradbell $ */
# ifndef CPPAD_CONFIGURE_INCLUDED
# define CPPAD_CONFIGURE_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-12 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

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
 \file configure.hpp
Replacement for config.h so that all preprocessor symbols begin with CPPAD_ 
*/

/*!
\def CPPAD_PACKAGE_STRING
cppad-yyyymmdd as a C string where yyyy is year, mm is month, and dd is day.
*/
# define CPPAD_PACKAGE_STRING "cppad-20120101.3"

/*!
\def CPPAD_BOOSTVECTOR
If this symbol is one, and _MSC_VER is not defined,
we are using boost vector for CPPAD_TEST_VECTOR.
It this symbol is zero, 
we are not using boost vector for CPPAD_TEST_VECTOR.
*/
# define CPPAD_BOOSTVECTOR 0

/*!
\def CPPAD_CPPADVECTOR
If this symbol is one, 
we are using CppAD vector for CPPAD_TEST_VECTOR.
It this symbol is zero, 
we are not using CppAD vector for CPPAD_TEST_VECTOR.
*/
# define CPPAD_CPPADVECTOR 1

/*!
\def CPPAD_STDVECTOR
If this symbol is one, 
we are using standard vector for CPPAD_TEST_VECTOR.
It this symbol is zero, 
we are not using standard vector for CPPAD_TEST_VECTOR.
*/
# define CPPAD_STDVECTOR 0

/*!
\def CPPAD_GETTIMEOFDAY
If this symbol is one, and _MSC_VER is not defined,
this system supports the gettimeofday funcgtion.
Otherwise, this smybol should be zero.
*/
# define CPPAD_GETTIMEOFDAY 1

/*!
\def CPPAD_SIZE_T_SAME_UNSIGNED_INT 
If this symgole is one, sizeof(size_t) == sizeof(unsigned int),
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

# endif

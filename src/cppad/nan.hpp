/* $Id: nan.hpp 2085 2011-09-01 14:54:04Z bradbell $ */
# ifndef CPPAD_NAN_INCLUDED
# define CPPAD_NAN_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
/*
$begin nan$$
$spell
	hasnan
	cppad
	hpp
	CppAD
	isnan
	bool
	const
$$

$section Obtain Nan or Determine if a Value is Nan$$

$index isnan$$
$index hasnan$$
$index nan$$

$head Syntax$$
$codei%# include <cppad/nan.hpp>
%$$
$icode%s% = nan(%z%)
%$$
$icode%b% = isnan(%s%)
%$$
$icode%b% = hasnan(%v%)%$$

$head Purpose$$
It obtain and check for the value not a number $code nan$$.
The IEEE standard specifies that a floating point value $icode a$$ 
is $code nan$$ if and only if the following returns true
$codei%
	%a% != %a%
%$$ 
Some systems do not get this correct, so we also use the fact that
zero divided by zero should result in a $code nan$$.
To be specific, if a value is not equal to itself or 
if it is equal to zero divided by zero, it is considered to be a $code nan$$.

$head Include$$
The file $code cppad/nan.hpp$$ is included by $code cppad/cppad.hpp$$
but it can also be included separately with out the rest of 
the $code CppAD$$ routines.

$subhead Macros$$
$index macro, nan$$
$index macro,  isnan$$
$index nan, macro$$
$index isnan, macro$$
Some C++ compilers use preprocessor symbols called $code nan$$ 
and $code isnan$$.
These preprocessor symbols will no longer be defined after 
this file is included. 

$head nan$$
This routine returns a $code nan$$ with the same type as $icode z$$.

$subhead z$$
The argument $icode z$$ has prototype
$codei%
	const %Scalar% &%z% 
%$$
and its value is zero
(see $cref/Scalar/nan/Scalar/$$ for the definition of $icode Scalar$$).

$subhead s$$
The return value $icode s$$ has prototype
$codei%
	%Scalar% %s%
%$$
It is the value $code nan$$ for this floating point type.

$head isnan$$
This routine determines if a scalar value is $code nan$$.

$subhead s$$
The argument $icode s$$ has prototype
$codei%
	const %Scalar% %s%
%$$

$subhead b$$
The return value $icode b$$ has prototype
$codei%
	bool %b%
%$$
It is true if the value $icode s$$ is $code nan$$.

$head hasnan$$
This routine determines if a 
$cref/SimpleVector/$$ has an element that is $code nan$$.

$subhead v$$
The argument $icode v$$ has prototype
$codei%
	const %Vector% &%v%
%$$
(see $cref/Vector/nan/Vector/$$ for the definition of $icode Vector$$).

$subhead b$$
The return value $icode b$$ has prototype
$codei%
	bool %b%
%$$
It is true if the vector $icode v$$ has a $code nan$$.

$head Scalar$$
The type $icode Scalar$$ must support the following operations;
$table
$bold Operation$$ $cnext $bold Description$$  $rnext
$icode%a% / %b%$$ $cnext
	division operator (returns a $icode Scalar$$ object)
$rnext
$icode%a% == %b%$$ $cnext
	equality operator (returns a $code bool$$ object)
$rnext
$icode%a% != %b%$$ $cnext
	not equality operator (returns a $code bool$$ object)
$tend
Note that the division operator will be used with $icode a$$ and $italic b$$
equal to zero. For some types (e.g. $code int$$) this may generate
an exception. No attempt is made to catch any such exception.

$head Vector$$
The type $icode Vector$$ must be a $cref/SimpleVector/$$ class with
elements of type $icode Scalar$$.

$head Parallel Mode$$
$index parallel, user_atomic$$
$index user_atomic, parallel$$
For each type $icode Scalar$$,
the first call to
$codei%
	%b% = isnan(%s%)
%$$
must not be $cref/parallel/ta_in_parallel/$$ execution mode; 
see $code isnan$$ in $cref/parallel_ad/parallel_ad/isnan/$$.

$children%
	example/nan.cpp
%$$
$head Example$$
The file $cref/nan.cpp/$$
contains an example and test of this routine.
It returns true if it succeeds and false otherwise.

$end
*/

# include <cstddef>
# include <cppad/local/cppad_assert.hpp>

// needed before one can use CPPAD_ASSERT_FIRST_CALL_NOT_PARALLEL
# include <cppad/thread_alloc.hpp>

# ifdef nan
# undef nan
# endif
# ifdef isnan
# undef isnan
# endif

namespace CppAD { // BEGIN CppAD namespace

template <class Scalar>
inline Scalar nan(const Scalar &zero)
{	return zero / zero;
}

template <class Scalar>
inline bool isnan(const Scalar &s)
{	CPPAD_ASSERT_FIRST_CALL_NOT_PARALLEL;	
	static Scalar scalar_nan = nan( Scalar(0) );	
	return (s != s) | (s == scalar_nan);
}

template <class Vector>
bool hasnan(const Vector &v)
{
	typedef typename Vector::value_type Scalar;

	bool found_nan;
	size_t i;
	i   = v.size();
	found_nan = false;
	while(i--)
		found_nan |= isnan(v[i]);
	return found_nan;
}

} // End CppAD namespace

# endif

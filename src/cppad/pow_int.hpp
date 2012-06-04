/* $Id: pow_int.hpp 1370 2009-05-31 05:31:50Z bradbell $ */
# ifndef CPPAD_POW_INT_INCLUDED
# define CPPAD_POW_INT_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-08 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
-------------------------------------------------------------------------------
$begin pow_int$$
$spell
	cppad.hpp
	CppAD
	namespace
	const
$$

$index pow, integer$$
$index exponent, integer$$
$index integer, pow$$

$section The Integer Power Function$$

$head Syntax$$
$code # include <cppad/pow_int.h>$$
$pre
$$
$syntax%%z% = pow(%x%, %y%)%$$

$head Purpose$$
Determines the value of the power function 
$latex \[
	{\rm pow} (x, y) = x^y
\] $$
for integer exponents $italic n$$ 
using multiplication and possibly division to compute the value.
The other CppAD $cref/pow/$$ function may use logarithms and exponentiation 
to compute derivatives of the same value
(which will not work if $italic x$$ is less than or equal zero).

$head Include$$
The file $code cppad/pow_int.h$$ is included by $code cppad/cppad.hpp$$
but it can also be included separately with out the rest of 
the $code CppAD$$ routines.
Including this file defines
this version of the $code pow$$ within the $code CppAD$$ namespace.

$head x$$
The argument $italic x$$ has prototype
$syntax%
	const %Type% &%x%
%$$ 

$head y$$
The argument $italic y$$ has prototype
$syntax%
	const int &%y%
%$$ 

$head z$$
The result $italic z$$ has prototype
$syntax%
	%Type% %z%
%$$

$head Type$$
The type $italic Type$$ must support the following operations
where $italic a$$ and $italic b$$ are $italic Type$$ objects
and $italic i$$ is an $code int$$:
$table
$bold Operation$$  $pre  $$
	$cnext $bold Description$$ 
	$cnext $bold Result Type$$ 
$rnext
$syntax%%Type% %a%(%i%)%$$ 
	$cnext construction of a $italic Type$$ object from an $code int$$
	$cnext $italic Type$$
$rnext
$syntax%%a% * %b%$$ 
	$cnext binary multiplication of $italic Type$$ objects
	$cnext $italic Type$$
$rnext
$syntax%%a% / %b%$$ 
	$cnext binary division of $italic Type$$ objects
	$cnext $italic Type$$
$tend

$head Operation Sequence$$
The $italic Type$$ operation sequence used to calculate $italic z$$ is 
$xref/glossary/Operation/Independent/independent/1/$$
of $italic x$$.


$end
-------------------------------------------------------------------------------
*/

namespace CppAD { 

	template <class Type>
	inline Type pow (const Type &x, const int &n)
	{
		Type p(1);
		int n2 = n / 2;

		if( n == 0 )
			return p;
		if( n < 0 )
			return p / pow(x, -n);
		if( n == 1 )
			return x;

		// p = (x^2)^(n/2)
		p = pow( x * x , n2 );

		// n is even case
		if( n % 2 == 0 )
			return p;

		// n is odd case
		return p * x;
	}

}

# endif 

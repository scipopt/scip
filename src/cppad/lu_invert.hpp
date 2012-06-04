/* $Id: lu_invert.hpp 1370 2009-05-31 05:31:50Z bradbell $ */
# ifndef CPPAD_LU_INVERT_INCLUDED
# define CPPAD_LU_INVERT_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-06 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin LuInvert$$
$escape #$$
$spell
	cppad.hpp
	Lu
	Cpp
	jp
	ip
	const
	namespace
	typename
	etmp
$$

$index LuInvert$$
$index linear, invert Lu equation$$
$index equation, Lu invert$$

$section Invert an LU Factored Equation$$

$pre
$$

$head Syntax$$ $code# include <cppad/lu_invert.hpp>$$
$pre
$$
$syntax%LuInvert(%ip%, %jp%, %LU%, %X%)%$$


$head Description$$
Solves the matrix equation $syntax%%A% * %X% = %B%$$ 
using an LU factorization computed by $xref/LuFactor/$$.

$head Include$$
The file $code cppad/lu_invert.hpp$$ is included by $code cppad/cppad.hpp$$
but it can also be included separately with out the rest of 
the $code CppAD$$ routines.

$head Matrix Storage$$
All matrices are stored in row major order.
To be specific, if $latex Y$$ is a vector
that contains a $latex p$$ by $latex q$$ matrix,
the size of $latex Y$$ must be equal to $latex  p * q $$ and for
$latex i = 0 , \ldots , p-1$$,
$latex j = 0 , \ldots , q-1$$,
$latex \[
	Y_{i,j} = Y[ i * q + j ]
\] $$

$head ip$$
The argument $italic ip$$ has prototype
$syntax%
	const %SizeVector% &%ip%
%$$
(see description for $italic SizeVector$$ in
$xref/LuFactor/SizeVector/LuFactor/$$ specifications).
The size of $italic ip$$ is referred to as $italic n$$ in the
specifications below.
The elements of $italic ip$$ determine
the order of the rows in the permuted matrix.

$head jp$$
The argument $italic jp$$ has prototype
$syntax%
	const %SizeVector% &%jp%
%$$
(see description for $italic SizeVector$$ in
$xref/LuFactor/SizeVector/LuFactor/$$ specifications).
The size of $italic jp$$ must be equal to $italic n$$.
The elements of $italic jp$$ determine
the order of the columns in the permuted matrix.

$head LU$$
The argument $italic LU$$ has the prototype
$syntax%
	const %FloatVector% &%LU%
%$$
and the size of $italic LU$$ must equal $latex n * n$$
(see description for $italic FloatVector$$ in
$xref/LuFactor/FloatVector/LuFactor/$$ specifications).

$subhead L$$
We define the lower triangular matrix $italic L$$ in terms of $italic LU$$.
The matrix $italic L$$ is zero above the diagonal
and the rest of the elements are defined by
$syntax%
	%L%(%i%, %j%) = %LU%[ %ip%[%i%] * %n% + %jp%[%j%] ]
%$$
for $latex i = 0 , \ldots , n-1$$ and $latex j = 0 , \ldots , i$$.

$subhead U$$
We define the upper triangular matrix $italic U$$ in terms of $italic LU$$.
The matrix $italic U$$ is zero below the diagonal,
one on the diagonal,
and the rest of the elements are defined by
$syntax%
	%U%(%i%, %j%) = %LU%[ %ip%[%i%] * %n% + %jp%[%j%] ]
%$$
for $latex i = 0 , \ldots , n-2$$ and $latex j = i+1 , \ldots , n-1$$.

$subhead P$$
We define the permuted matrix $italic P$$ in terms of 
the matrix $italic L$$ and the matrix $italic U$$ 
by $syntax%%P% = %L% * %U%$$.

$subhead A$$
The matrix $italic A$$, 
which defines the linear equations that we are solving, is given by
$syntax%
	%P%(%i%, %j%) = %A%[ %ip%[%i%] * %n% + %jp%[%j%] ]
%$$
(Hence 
$italic LU$$ contains a permuted factorization of the matrix $italic A$$.)


$head X$$
The argument $italic X$$ has prototype
$syntax%
	%FloatVector% &%X%
%$$
(see description for $italic FloatVector$$ in
$xref/LuFactor/FloatVector/LuFactor/$$ specifications).
The matrix $italic X$$
must have the same number of rows as the matrix $italic A$$.
The input value of $italic X$$ is the matrix $italic B$$ and the 
output value solves the matrix equation $syntax%%A% * %X% = %B%$$.


$children%
	example/lu_invert.cpp%
	omh/lu_invert_hpp.omh
%$$
$head Example$$
The file $xref/lu_solve.hpp/$$ is a good example usage of 
$code LuFactor$$ with $code LuInvert$$.
The file 
$xref/LuInvert.cpp/$$
contains an example and test of using $code LuInvert$$ by itself.
It returns true if it succeeds and false otherwise.

$head Source$$
The file $cref/lu_invert.hpp/$$ contains the
current source code that implements these specifications.

$end
--------------------------------------------------------------------------
*/
// BEGIN PROGRAM
# include <cppad/local/cppad_assert.hpp>
# include <cppad/check_simple_vector.hpp>
# include <cppad/check_numeric_type.hpp>

namespace CppAD { // BEGIN CppAD namespace

// LuInvert
template <typename SizeVector, typename FloatVector>
void LuInvert(
	const SizeVector  &ip, 
	const SizeVector  &jp, 
	const FloatVector &LU,
	FloatVector       &B )
{	size_t k; // column index in X
	size_t p; // index along diagonal in LU
	size_t i; // row index in LU and X

	typedef typename FloatVector::value_type Float;

	// check numeric type specifications
	CheckNumericType<Float>();

	// check simple vector class specifications
	CheckSimpleVector<Float, FloatVector>();
	CheckSimpleVector<size_t, SizeVector>();

	Float etmp;
	
	size_t n = ip.size();
	CPPAD_ASSERT_KNOWN(
		jp.size() == n,
		"Error in LuInvert: jp must have size equal to n * n"
	);
	CPPAD_ASSERT_KNOWN(
		LU.size() == n * n,
		"Error in LuInvert: Lu must have size equal to n * m"
	);
	size_t m = B.size() / n;
	CPPAD_ASSERT_KNOWN(
		B.size() == n * m,
		"Error in LuSolve: B must have size equal to a multiple of n"
	);

	// temporary storage for reordered solution
	FloatVector x(n);

	// loop over equations
	for(k = 0; k < m; k++)
	{	// invert the equation c = L * b
		for(p = 0; p < n; p++)
		{	// solve for c[p]
			etmp = B[ ip[p] * m + k ] / LU[ ip[p] * n + jp[p] ];
			B[ ip[p] * m + k ] = etmp;
			// subtract off effect on other variables
			for(i = p+1; i < n; i++)
				B[ ip[i] * m + k ] -=
					etmp * LU[ ip[i] * n + jp[p] ];
		}

		// invert the equation x = U * c
		p = n;
		while( p > 0 )
		{	--p;
			etmp       = B[ ip[p] * m + k ];
			x[ jp[p] ] = etmp;
			for(i = 0; i < p; i++ )
				B[ ip[i] * m + k ] -= 
					etmp * LU[ ip[i] * n + jp[p] ];
		}

		// copy reordered solution into B
		for(i = 0; i < n; i++)
			B[i * m + k] = x[i];
	}
	return;
}
} // END CppAD namespace 
// END PROGRAM
# endif

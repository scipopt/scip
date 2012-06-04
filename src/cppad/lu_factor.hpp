/* $Id: lu_factor.hpp 1370 2009-05-31 05:31:50Z bradbell $ */
# ifndef CPPAD_LU_FACTOR_INCLUDED
# define CPPAD_LU_FACTOR_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-08 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin LuFactor$$
$escape #$$
$spell
	cppad.hpp
	Cpp
	Geq
	Lu
	bool
	const
	ip
	jp
	namespace
	std
	typename
$$

$index LuFactor$$
$index linear, Lu factor equation$$
$index equation, Lu factor$$
$index determinant, Lu factor$$
$index solve, Lu factor$$

$section LU Factorization of A Square Matrix$$

$pre
$$

$head Syntax$$ $code# include <cppad/lu_factor.hpp>$$
$pre
$$
$syntax%%sign% = LuFactor(%ip%, %jp%, %LU%)%$$


$head Description$$
Computes an LU factorization of the matrix $italic A$$ 
where $italic A$$ is a square matrix.

$head Include$$
The file $code cppad/lu_factor.hpp$$ is included by $code cppad/cppad.hpp$$
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

$head sign$$
The return value $italic sign$$ has prototype
$syntax%
	int %sign%
%$$
If $italic A$$ is invertible, $italic sign$$ is plus or minus one
and is the sign of the permutation corresponding to the row ordering
$italic ip$$ and column ordering $italic jp$$.
If $italic A$$ is not invertible, $italic sign$$ is zero.

$head ip$$
The argument $italic ip$$ has prototype
$syntax%
	%SizeVector% &%ip%
%$$
(see description of $xref/LuFactor/SizeVector/SizeVector/$$ below).
The size of $italic ip$$ is referred to as $italic n$$ in the
specifications below.
The input value of the elements of $italic ip$$ does not matter.
The output value of the elements of $italic ip$$ determine
the order of the rows in the permuted matrix.

$head jp$$
The argument $italic jp$$ has prototype
$syntax%
	%SizeVector% &%jp%
%$$
(see description of $xref/LuFactor/SizeVector/SizeVector/$$ below).
The size of $italic jp$$ must be equal to $italic n$$.
The input value of the elements of $italic jp$$ does not matter.
The output value of the elements of $italic jp$$ determine
the order of the columns in the permuted matrix.

$head LU$$
The argument $italic LU$$ has the prototype
$syntax%
	%FloatVector% &%LU%
%$$
and the size of $italic LU$$ must equal $latex n * n$$
(see description of $xref/LuFactor/FloatVector/FloatVector/$$ below).

$subhead A$$
We define $italic A$$ as the matrix corresponding to the input 
value of $italic LU$$.

$subhead P$$
We define the permuted matrix $italic P$$ in terms of $italic A$$ by
$syntax%
	%P%(%i%, %j%) = %A%[ %ip%[%i%] * %n% + %jp%[%j%] ]
%$$

$subhead L$$
We define the lower triangular matrix $italic L$$ in terms of the 
output value of $italic LU$$.
The matrix $italic L$$ is zero above the diagonal
and the rest of the elements are defined by
$syntax%
	%L%(%i%, %j%) = %LU%[ %ip%[%i%] * %n% + %jp%[%j%] ]
%$$
for $latex i = 0 , \ldots , n-1$$ and $latex j = 0 , \ldots , i$$.

$subhead U$$
We define the upper triangular matrix $italic U$$ in terms of the
output value of $italic LU$$.
The matrix $italic U$$ is zero below the diagonal,
one on the diagonal,
and the rest of the elements are defined by
$syntax%
	%U%(%i%, %j%) = %LU%[ %ip%[%i%] * %n% + %jp%[%j%] ]
%$$
for $latex i = 0 , \ldots , n-2$$ and $latex j = i+1 , \ldots , n-1$$.

$subhead Factor$$
If the return value $italic sign$$ is non-zero,
$syntax%
	%L% * %U% = %P%
%$$
If the return value of $italic sign$$ is zero,
the contents of $italic L$$ and $italic U$$ are not defined. 

$subhead Determinant$$
$index determinant$$
If the return value $italic sign$$ is zero,
the determinant of $italic A$$ is zero.
If $italic sign$$ is non-zero,
using the output value of $italic LU$$
the determinant of the matrix $italic A$$ is equal to
$syntax%
%sign% * %LU%[%ip%[0], %jp%[0]] * %...% * %LU%[%ip%[%n%-1], %jp%[%n%-1]] 
%$$

$head SizeVector$$
The type $italic SizeVector$$ must be a $xref/SimpleVector/$$ class with
$xref/SimpleVector/Elements of Specified Type/elements of type size_t/$$.
The routine $xref/CheckSimpleVector/$$ will generate an error message
if this is not the case.

$head FloatVector$$
The type $italic FloatVector$$ must be a 
$xref/SimpleVector//simple vector class/$$.
The routine $xref/CheckSimpleVector/$$ will generate an error message
if this is not the case.

$head Float$$
This notation is used to denote the type corresponding
to the elements of a $italic FloatVector$$.
The type $italic Float$$ must satisfy the conditions
for a $xref/NumericType/$$ type.
The routine $xref/CheckNumericType/$$ will generate an error message
if this is not the case.
In addition, the following operations must be defined for any pair
of $italic Float$$ objects $italic x$$ and $italic y$$:

$table
$bold Operation$$ $cnext $bold Description$$  $rnext
$syntax%log(%x%)%$$ $cnext
	returns the logarithm of $italic x$$ as a $italic Float$$ object
$tend

$head AbsGeq$$
Including the file $code lu_factor.hpp$$ defines the template function 
$syntax%
	template <typename %Float%>
	bool AbsGeq<%Float%>(const %Float% &%x%, const %Float% &%y%)
%$$
in the $code CppAD$$ namespace.
This function returns true if the absolute value of 
$italic x$$ is greater than or equal the absolute value of $italic y$$. 
It is used by $code LuFactor$$ to choose the pivot elements.
This template function definition uses the operator 
$code <=$$ to obtain the absolute value for $italic Float$$ objects. 
If this operator is not defined for your use of $italic Float$$,
you will need to specialize this template so that it works for your
use of $code LuFactor$$.
$pre

$$
Complex numbers do not have the operation $code <=$$ defined.
The specializations
$syntax%
bool AbsGeq< std::complex<float> > 
	(const std::complex<float> &%x%, const std::complex<float> &%y%)
bool AbsGeq< std::complex<double> > 
	(const std::complex<double> &%x%, const std::complex<double> &%y%)
%$$ 
are define by including $code lu_factor.hpp$$
These return true if the sum of the square of the real and imaginary parts
of $italic x$$ is greater than or equal the 
sum of the square of the real and imaginary parts of $italic y$$. 

$children%
	example/lu_factor.cpp%
	omh/lu_factor_hpp.omh
%$$
$head Example$$
The file 
$xref/LuFactor.cpp/$$
contains an example and test of using $code LuFactor$$ by itself.
It returns true if it succeeds and false otherwise.
$pre

$$
The file $xref/lu_solve.hpp/$$ provides a useful example usage of 
$code LuFactor$$ with $code LuInvert$$.

$head Source$$
The file $cref/lu_factor.hpp/$$ contains the 
current source code that implements these specifications.

$end
--------------------------------------------------------------------------
*/
// BEGIN PROGRAM

# include <complex>
# include <vector>

# include <cppad/local/cppad_assert.hpp>
# include <cppad/check_simple_vector.hpp>
# include <cppad/check_numeric_type.hpp>

namespace CppAD { // BEGIN CppAD namespace

// AbsGeq
template <typename Float>
inline bool AbsGeq(const Float &x, const Float &y)
{	Float xabs = x;
	if( xabs <= Float(0) )
		xabs = - xabs;
	Float yabs = y;
	if( yabs <= Float(0) )
		yabs = - yabs;
	return xabs >= yabs;
}
inline bool AbsGeq(
	const std::complex<double> &x, 
	const std::complex<double> &y)
{	double xsq = x.real() * x.real() + x.imag() * x.imag();
	double ysq = y.real() * y.real() + y.imag() * y.imag();

	return xsq >= ysq;
}
inline bool AbsGeq(
	const std::complex<float> &x, 
	const std::complex<float> &y)
{	float xsq = x.real() * x.real() + x.imag() * x.imag();
	float ysq = y.real() * y.real() + y.imag() * y.imag();

	return xsq >= ysq;
}

// Lines that are different from code in cppad/local/lu_ratio.hpp end with //
template <class SizeVector, class FloatVector>                          //
int LuFactor(SizeVector &ip, SizeVector &jp, FloatVector &LU)           //
{	
	// type of the elements of LU                                   //
	typedef typename FloatVector::value_type Float;                 //

	// check numeric type specifications
	CheckNumericType<Float>();

	// check simple vector class specifications
	CheckSimpleVector<Float, FloatVector>();
	CheckSimpleVector<size_t, SizeVector>();

	size_t  i, j;          // some temporary indices
	const Float zero( 0 ); // the value zero as a Float object
	size_t  imax;          // row index of maximum element
	size_t  jmax;          // column indx of maximum element
	Float    emax;         // maximum absolute value
	size_t  p;             // count pivots
	int     sign;          // sign of the permutation
	Float   etmp;          // temporary element
	Float   pivot;         // pivot element

	// -------------------------------------------------------
	size_t n = ip.size();
	CPPAD_ASSERT_KNOWN(
		jp.size() == n,
		"Error in LuFactor: jp must have size equal to n"
	);
	CPPAD_ASSERT_KNOWN(
		LU.size() == n * n,
		"Error in LuFactor: LU must have size equal to n * m"
	);
	// -------------------------------------------------------

	// initialize row and column order in matrix not yet pivoted
	for(i = 0; i < n; i++)
	{	ip[i] = i;
		jp[i] = i;
	}
	// initialize the sign of the permutation
	sign = 1;
	// ---------------------------------------------------------

	// Reduce the matrix P to L * U using n pivots
	for(p = 0; p < n; p++)
	{	// determine row and column corresponding to element of 
		// maximum absolute value in remaining part of P
		imax = jmax = n;
		emax = zero;
		for(i = p; i < n; i++)
		{	for(j = p; j < n; j++)
			{	CPPAD_ASSERT_UNKNOWN(
					(ip[i] < n) & (jp[j] < n)
				);
				etmp = LU[ ip[i] * n + jp[j] ];

				// check if maximum absolute value so far
				if( AbsGeq (etmp, emax) )
				{	imax = i;
					jmax = j;
					emax = etmp;
				}
			}
		}
		CPPAD_ASSERT_KNOWN( 
		(imax < n) & (jmax < n) ,
		"LuFactor can't determine an element with "
		"maximum absolute value.\n"
		"Perhaps original matrix contains not a number or infinity.\n" 
		"Perhaps your specialization of AbsGeq is not correct."
		);
		if( imax != p )
		{	// switch rows so max absolute element is in row p
			i        = ip[p];
			ip[p]    = ip[imax];
			ip[imax] = i;
			sign     = -sign;
		}
		if( jmax != p )
		{	// switch columns so max absolute element is in column p
			j        = jp[p];
			jp[p]    = jp[jmax];
			jp[jmax] = j;
			sign     = -sign;
		}
		// pivot using the max absolute element
		pivot   = LU[ ip[p] * n + jp[p] ];

		// check for determinant equal to zero
		if( pivot == zero )
		{	// abort the mission
			return   0;
		}

		// Reduce U by the elementary transformations that maps 
		// LU( ip[p], jp[p] ) to one.  Only need transform elements
		// above the diagonal in U and LU( ip[p] , jp[p] ) is
		// corresponding value below diagonal in L.
		for(j = p+1; j < n; j++)
			LU[ ip[p] * n + jp[j] ] /= pivot;

		// Reduce U by the elementary transformations that maps 
		// LU( ip[i], jp[p] ) to zero. Only need transform elements 
		// above the diagonal in U and LU( ip[i], jp[p] ) is 
		// corresponding value below diagonal in L.
		for(i = p+1; i < n; i++ )
		{	etmp = LU[ ip[i] * n + jp[p] ];
			for(j = p+1; j < n; j++)
			{	LU[ ip[i] * n + jp[j] ] -= 
					etmp * LU[ ip[p] * n + jp[j] ];
			} 
		}
	}
	return sign;
}
} // END CppAD namespace 
// END PROGRAM
# endif

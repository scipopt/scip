/* $Id: lu_solve.hpp 2057 2011-08-11 14:07:11Z bradbell $ */
# ifndef CPPAD_LU_SOLVE_INCLUDED
# define CPPAD_LU_SOLVE_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin LuSolve$$
$escape #$$
$spell
	cppad.hpp
	det
	exp
	Leq
	typename
	bool
	const
	namespace
	std
	Geq
	Lu
	CppAD
	signdet
	logdet
$$

$index LuSolve$$
$index linear, equation$$
$index equation, linear$$
$index determinant, Lu$$
$index solve, linear equation$$

$section Compute Determinant and Solve Linear Equations$$

$pre
$$

$head Syntax$$ $code# include <cppad/lu_solve.hpp>$$
$pre
$$
$syntax%%signdet% = LuSolve(%n%, %m%, %A%, %B%, %X%, %logdet%)%$$


$head Description$$
Use an LU factorization of the matrix $italic A$$ to
compute its determinant 
and solve for $italic X$$ in the linear of equation
$latex \[
	A * X = B
\] $$
where $italic A$$ is an 
$italic n$$ by $italic n$$ matrix,
$italic X$$ is an 
$italic n$$ by $italic m$$ matrix, and
$italic B$$ is an $latex n x m$$ matrix.

$head Include$$
The file $code cppad/lu_solve.hpp$$ is included by $code cppad/cppad.hpp$$
but it can also be included separately with out the rest of 
the $code CppAD$$ routines.

$head Factor and Invert$$
This routine is an easy to user interface to
$xref/LuFactor/$$ and $xref/LuInvert/$$ for computing determinants and
solutions of linear equations.
These separate routines should be used if
one right hand side $italic B$$
depends on the solution corresponding to another
right hand side (with the same value of $italic A$$).
In this case only one call to $code LuFactor$$ is required
but there will be multiple calls to $code LuInvert$$.


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

$head signdet$$
The return value $italic signdet$$ is a $code int$$ value
that specifies the sign factor for the determinant of $italic A$$.
This determinant of $italic A$$ is zero if and only if $italic signdet$$
is zero. 

$head n$$
The argument $italic n$$ has type $code size_t$$ 
and specifies the number of rows in the matrices
$italic A$$,
$italic X$$,
and $italic B$$.
The number of columns in $italic A$$ is also equal to $italic n$$.

$head m$$
The argument $italic m$$ has type $code size_t$$ 
and specifies the number of columns in the matrices
$italic X$$
and $italic B$$.
If $italic m$$ is zero,
only the determinant of $italic A$$ is computed and
the matrices $italic X$$ and $italic B$$ are not used.

$head A$$
The argument $italic A$$ has the prototype
$syntax%
	const %FloatVector% &%A%
%$$
and the size of $italic A$$ must equal $latex n * n$$
(see description of $xref/LuSolve/FloatVector/FloatVector/$$ below).
This is the $latex n$$ by $italic n$$ matrix that 
we are computing the determinant of 
and that defines the linear equation.

$head B$$
The argument $italic B$$ has the prototype
$syntax%
	const %FloatVector% &%B%
%$$
and the size of $italic B$$ must equal $latex n * m$$
(see description of $xref/LuSolve/FloatVector/FloatVector/$$ below).
This is the $latex n$$ by $italic m$$ matrix that 
defines the right hand side of the linear equations.
If $italic m$$ is zero, $italic B$$ is not used.

$head X$$
The argument $italic X$$ has the prototype
$syntax%
	%FloatVector% &%X%
%$$
and the size of $italic X$$ must equal $latex n * m$$
(see description of $xref/LuSolve/FloatVector/FloatVector/$$ below).
The input value of $italic X$$ does not matter.
On output, the elements of $italic X$$ contain the solution
of the equation we wish to solve
(unless $italic signdet$$ is equal to zero).
If $italic m$$ is zero, $italic X$$ is not used.

$head logdet$$
The argument $italic logdet$$ has prototype
$syntax%
	%Float% &%logdet%
%$$
On input, the value of $italic logdet$$ does not matter.
On output, it has been set to the 
log of the determinant of $italic A$$ 
(but not quite).
To be more specific,
the determinant of $italic A$$ is given by the formula
$syntax%
	%det% = %signdet% * exp( %logdet% )
%$$
This enables $code LuSolve$$ to use logs of absolute values
in the case where $italic Float$$ corresponds to a real number.

$head Float$$
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

$head FloatVector$$
The type $italic FloatVector$$ must be a $xref/SimpleVector/$$ class with
$xref/SimpleVector/Elements of Specified Type/elements of type Float/$$.
The routine $xref/CheckSimpleVector/$$ will generate an error message
if this is not the case.

$head LeqZero$$
Including the file $code lu_solve.hpp$$ defines the template function 
$syntax%
	template <typename %Float%>
	bool LeqZero<%Float%>(const %Float% &%x%)
%$$
in the $code CppAD$$ namespace.
This function returns true if $italic x$$ is less than or equal to zero
and false otherwise.
It is used by $code LuSolve$$ to avoid taking the log of
zero (or a negative number if $italic Float$$ corresponds to real numbers).
This template function definition assumes that the operator 
$code <=$$ is defined for $italic Float$$ objects. 
If this operator is not defined for your use of $italic Float$$,
you will need to specialize this template so that it works for your
use of $code LuSolve$$.
$pre

$$
Complex numbers do not have the operation or $code <=$$ defined.
In addition, in the complex case, 
one can take the log of a negative number.
The specializations
$syntax%
	bool LeqZero< std::complex<float> > (const std::complex<float> &%x%)
	bool LeqZero< std::complex<double> >(const std::complex<double> &%x%)
%$$ 
are defined by including $code lu_solve.hpp$$.
These return true if $italic x$$ is zero and false otherwise.

$head AbsGeq$$
Including the file $code lu_solve.hpp$$ defines the template function 
$syntax%
	template <typename %Float%>
	bool AbsGeq<%Float%>(const %Float% &%x%, const %Float% &%y%)
%$$
If the type $italic Float$$ does not support the $code <=$$ operation
and it is not $code std::complex<float>$$ or $code std::complex<double>$$,
see the documentation for $code AbsGeq$$ in $xref/LuFactor/AbsGeq/LuFactor/$$. 

$children%
	example/lu_solve.cpp%
	omh/lu_solve_hpp.omh
%$$
$head Example$$
The file 
$xref/LuSolve.cpp/$$
contains an example and test of using this routine.
It returns true if it succeeds and false otherwise.

$head Source$$
The file $cref/lu_solve.hpp/$$ contains the
current source code that implements these specifications.

$end
--------------------------------------------------------------------------
*/
// BEGIN PROGRAM
# include <complex>
# include <vector>

// link exp for float and double cases
# include <cppad/base_require.hpp>

# include <cppad/local/cppad_assert.hpp>
# include <cppad/check_simple_vector.hpp>
# include <cppad/check_numeric_type.hpp>
# include <cppad/lu_factor.hpp>
# include <cppad/lu_invert.hpp>

namespace CppAD { // BEGIN CppAD namespace

// LeqZero
template <typename Float>
inline bool LeqZero(const Float &x)
{	return x <= Float(0); }
inline bool LeqZero( const std::complex<double> &x )
{	return x == std::complex<double>(0); }
inline bool LeqZero( const std::complex<float> &x )
{	return x == std::complex<float>(0); }

// LuSolve
template <typename Float, typename FloatVector>
int LuSolve(
	size_t             n      ,
	size_t             m      , 
	const FloatVector &A      , 
	const FloatVector &B      , 
	FloatVector       &X      , 
	Float        &logdet      )
{	
	// check numeric type specifications
	CheckNumericType<Float>();

	// check simple vector class specifications
	CheckSimpleVector<Float, FloatVector>();

	size_t        p;       // index of pivot element (diagonal of L)
	int     signdet;       // sign of the determinant
	Float     pivot;       // pivot element

	// the value zero
	const Float zero(0);

	// pivot row and column order in the matrix
	std::vector<size_t> ip(n);
	std::vector<size_t> jp(n);

	// -------------------------------------------------------
	CPPAD_ASSERT_KNOWN(
		A.size() == n * n,
		"Error in LuSolve: A must have size equal to n * n"
	);
	CPPAD_ASSERT_KNOWN(
		B.size() == n * m,
		"Error in LuSolve: B must have size equal to n * m"
	);
	CPPAD_ASSERT_KNOWN(
		X.size() == n * m,
		"Error in LuSolve: X must have size equal to n * m"
	);
	// -------------------------------------------------------

	// copy A so that it does not change
	FloatVector Lu(A);

	// copy B so that it does not change
	X = B;

	// Lu factor the matrix A
	signdet = LuFactor(ip, jp, Lu);

	// compute the log of the determinant
	logdet  = Float(0);
	for(p = 0; p < n; p++)
	{	// pivot using the max absolute element
		pivot   = Lu[ ip[p] * n + jp[p] ];

		// check for determinant equal to zero
		if( pivot == zero )
		{	// abort the mission
			logdet = Float(0);
			return   0;
		}

		// update the determinant
		if( LeqZero ( pivot ) )
		{	logdet += log( - pivot );
			signdet = - signdet;
		}
		else	logdet += log( pivot );

	}

	// solve the linear equations
	LuInvert(ip, jp, Lu, X);

	// return the sign factor for the determinant
	return signdet;
}
} // END CppAD namespace 
// END PROGRAM
# endif

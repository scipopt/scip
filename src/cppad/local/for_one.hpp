/* $Id: for_one.hpp 1598 2009-12-03 02:59:52Z bradbell $ */
# ifndef CPPAD_FOR_ONE_INCLUDED
# define CPPAD_FOR_ONE_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-09 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin ForOne$$
$spell
	dy
	typename
	Taylor
	const
$$


$index partial, first order driver$$
$index first, order partial driver$$
$index driver, first order partial$$

$index easy, partial$$
$index driver, easy partial$$
$index partial, easy$$

$section First Order Partial Derivative: Driver Routine$$

$head Syntax$$
$syntax%%dy% = %f%.ForOne(%x%, %j%)%$$


$head Purpose$$
We use $latex F : B^n \rightarrow B^m$$ to denote the
$xref/glossary/AD Function/AD function/$$ corresponding to $italic f$$.
The syntax above sets $italic dy$$ to the
partial of $latex F$$ with respect to $latex x_j$$; i.e.,
$latex \[
dy 
= \D{F}{ x_j } (x) 
= \left[ 
	\D{ F_0 }{ x_j } (x) , \cdots , \D{ F_{m-1} }{ x_j } (x) 
\right]
\] $$

$head f$$
The object $italic f$$ has prototype
$syntax%
	ADFun<%Base%> %f%
%$$
Note that the $xref/ADFun/$$ object $italic f$$ is not $code const$$
(see $xref/ForOne/ForOne Uses Forward/ForOne Uses Forward/$$ below).

$head x$$
The argument $italic x$$ has prototype
$syntax%
	const %Vector% &%x%
%$$
(see $xref/ForOne/Vector/Vector/$$ below)
and its size 
must be equal to $italic n$$, the dimension of the
$xref/seq_property/Domain/domain/$$ space for $italic f$$.
It specifies
that point at which to evaluate the partial derivative.

$head j$$
The argument $italic j$$ has prototype
$syntax%
	size_t %j%
%$$
an is less than $italic n$$,
$xref/seq_property/Domain/domain/$$ space for $italic f$$.
It specifies the component of $italic F$$ 
for which we are computing the partial derivative.

$head dy$$
The result $italic dy$$ has prototype
$syntax%
	%Vector% %dy%
%$$
(see $xref/ForOne/Vector/Vector/$$ below)
and its size is $latex m$$, the dimension of the
$xref/seq_property/Range/range/$$ space for $italic f$$.
The value of $italic dy$$ is the partial of $latex F$$ with respect to
$latex x_j$$ evaluated at $italic x$$; i.e.,
for $latex i = 0 , \ldots , m - 1$$
$latex \[.
	dy[i] = \D{ F_i }{ x_j } ( x )
\] $$


$head Vector$$
The type $italic Vector$$ must be a $xref/SimpleVector/$$ class with
$xref/SimpleVector/Elements of Specified Type/elements of type/$$
$italic Base$$.
The routine $xref/CheckSimpleVector/$$ will generate an error message
if this is not the case.

$head ForOne Uses Forward$$
After each call to $xref/Forward/$$,
the object $italic f$$ contains the corresponding 
$xref/glossary/Taylor Coefficient/Taylor coefficients/$$.
After $code ForOne$$,
the previous calls to $xref/Forward/$$ are undefined.

$head Example$$
$children%
	example/for_one.cpp
%$$
The routine 
$xref/ForOne.cpp//ForOne/$$ is both an example and test.
It returns $code true$$, if it succeeds and $code false$$ otherwise.

$end
-----------------------------------------------------------------------------
*/

//  BEGIN CppAD namespace
namespace CppAD {

template <typename Base>
template <typename Vector>
Vector ADFun<Base>::ForOne(const Vector &x, size_t j)
{	size_t j1;

	size_t n = Domain();
	size_t m = Range();

	// check Vector is Simple Vector class with Base type elements
	CheckSimpleVector<Base, Vector>();

	CPPAD_ASSERT_KNOWN(
		x.size() == n,
		"ForOne: Length of x not equal domain dimension for f"
	); 
	CPPAD_ASSERT_KNOWN(
		j < n,
		"ForOne: the index j is not less than domain dimension for f"
	);

	// point at which we are evaluating the second partials
	Forward(0, x);

	// direction in which are are taking the derivative
	Vector dx(n);
	for(j1 = 0; j1 < n; j1++)
		dx[j1] = Base(0);
	dx[j] = Base(1);

	// dimension the return value
	Vector dy(m);

	// compute the return value
	dy = Forward(1, dx);

	return dy;
}

} // END CppAD namespace

# endif

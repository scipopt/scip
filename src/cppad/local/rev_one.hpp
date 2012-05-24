/* $Id: rev_one.hpp 1598 2009-12-03 02:59:52Z bradbell $ */
# ifndef CPPAD_REV_ONE_INCLUDED
# define CPPAD_REV_ONE_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-09 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin RevOne$$
$spell
	dw
	Taylor
	const
$$


$index derivative, first order driver$$
$index first, order derivative driver$$
$index driver, first order derivative$$

$index easy, derivative$$
$index driver, easy derivative$$
$index derivative, easy$$

$section First Order Derivative: Driver Routine$$

$head Syntax$$
$syntax%%dw% = %f%.RevOne(%x%, %i%)%$$


$head Purpose$$
We use $latex F : B^n \rightarrow B^m$$ to denote the
$xref/glossary/AD Function/AD function/$$ corresponding to $italic f$$.
The syntax above sets $italic dw$$ to the
derivative of $latex F_i$$ with respect to $latex x$$; i.e.,
$latex \[
dw =
F_i^{(1)} (x) 
= \left[ 
	\D{ F_i }{ x_0 } (x) , \cdots , \D{ F_i }{ x_{n-1} } (x) 
\right]
\] $$

$head f$$
The object $italic f$$ has prototype
$syntax%
	ADFun<%Base%> %f%
%$$
Note that the $xref/ADFun/$$ object $italic f$$ is not $code const$$
(see $xref/RevOne/RevOne Uses Forward/RevOne Uses Forward/$$ below).

$head x$$
The argument $italic x$$ has prototype
$syntax%
	const %Vector% &%x%
%$$
(see $xref/RevOne/Vector/Vector/$$ below)
and its size 
must be equal to $italic n$$, the dimension of the
$xref/seq_property/Domain/domain/$$ space for $italic f$$.
It specifies
that point at which to evaluate the derivative.

$head i$$
The index $italic i$$ has prototype
$syntax%
	size_t %i%
%$$
and is less than $latex m$$, the dimension of the
$xref/seq_property/Range/range/$$ space for $italic f$$.
It specifies the
component of $latex F$$ that we are computing the derivative of.

$head dw$$
The result $italic dw$$ has prototype
$syntax%
	%Vector% %dw%
%$$
(see $xref/RevOne/Vector/Vector/$$ below)
and its size is $italic n$$, the dimension of the
$xref/seq_property/Domain/domain/$$ space for $italic f$$.
The value of $italic dw$$ is the derivative of $latex F_i$$ 
evaluated at $italic x$$; i.e.,
for $latex j = 0 , \ldots , n - 1 $$ 
$latex \[.
	dw[ j ] = \D{ F_i }{ x_j } ( x )
\] $$

$head Vector$$
The type $italic Vector$$ must be a $xref/SimpleVector/$$ class with
$xref/SimpleVector/Elements of Specified Type/elements of type/$$
$italic Base$$.
The routine $xref/CheckSimpleVector/$$ will generate an error message
if this is not the case.

$head RevOne Uses Forward$$
After each call to $xref/Forward/$$,
the object $italic f$$ contains the corresponding 
$xref/glossary/Taylor Coefficient/Taylor coefficients/$$.
After $code RevOne$$,
the previous calls to $xref/Forward/$$ are undefined.

$head Example$$
$children%
	example/rev_one.cpp
%$$
The routine 
$xref/RevOne.cpp//RevOne/$$ is both an example and test.
It returns $code true$$, if it succeeds and $code false$$ otherwise.

$end
-----------------------------------------------------------------------------
*/

//  BEGIN CppAD namespace
namespace CppAD {

template <typename Base>
template <typename Vector>
Vector ADFun<Base>::RevOne(const Vector  &x, size_t i)
{	size_t i1;

	size_t n = Domain();
	size_t m = Range();

	// check Vector is Simple Vector class with Base type elements
	CheckSimpleVector<Base, Vector>();

	CPPAD_ASSERT_KNOWN(
		x.size() == n,
		"RevOne: Length of x not equal domain dimension for f"
	); 
	CPPAD_ASSERT_KNOWN(
		i < m,
		"RevOne: the index i is not less than range dimension for f"
	);

	// point at which we are evaluating the derivative
	Forward(0, x);

	// component which are are taking the derivative of
	Vector w(m);
	for(i1 = 0; i1 < m; i1++)
		w[i1] = 0.;
	w[i] = Base(1);

	// dimension the return value
	Vector dw(n);

	// compute the return value
	dw = Reverse(1, w);

	return dw;
}

} // END CppAD namespace

# endif

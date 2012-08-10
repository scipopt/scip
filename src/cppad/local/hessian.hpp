/* $Id: hessian.hpp 1598 2009-12-03 02:59:52Z bradbell $ */
# ifndef CPPAD_HESSIAN_INCLUDED
# define CPPAD_HESSIAN_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-09 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin Hessian$$
$spell
	hes
	typename
	Taylor
	HesLuDet
	const
$$

$index Hessian, driver$$
$index second, derivative$$
$index driver, Hessian$$

$section Hessian: Easy Driver$$

$head Syntax$$
$syntax%%hes% = %f%.Hessian(%x%, %w%)
%$$
$syntax%%hes% = %f%.Hessian(%x%, %l%)
%$$


$head Purpose$$
We use $latex F : B^n \rightarrow B^m$$ to denote the
$xref/glossary/AD Function/AD function/$$ corresponding to $italic f$$.
The syntax above sets $italic hes$$ to the Hessian 
The syntax above sets $italic h$$ to the Hessian 
$latex \[
	hes = \dpow{2}{x} \sum_{i=1}^m w_i F_i (x) 
\] $$
The routine $cref/sparse_hessian/$$ may be faster in the case
where the Hessian is sparse.

$head f$$
The object $italic f$$ has prototype
$syntax%
	ADFun<%Base%> %f%
%$$
Note that the $xref/ADFun/$$ object $italic f$$ is not $code const$$
(see $xref/Hessian/Hessian Uses Forward/Hessian Uses Forward/$$ below).

$head x$$
The argument $italic x$$ has prototype
$syntax%
	const %Vector% &%x%
%$$
(see $xref/Hessian/Vector/Vector/$$ below)
and its size 
must be equal to $italic n$$, the dimension of the
$xref/seq_property/Domain/domain/$$ space for $italic f$$.
It specifies
that point at which to evaluate the Hessian.

$head l$$
If the argument $italic l$$ is present, it has prototype
$syntax%
	size_t %l%
%$$
and is less than $italic m$$, the dimension of the
$xref/seq_property/Range/range/$$ space for $italic f$$.
It specifies the component of $italic F$$
for which we are evaluating the Hessian.
To be specific, in the case where the argument $italic l$$ is present,
$latex \[
	w_i = \left\{ \begin{array}{ll}
		1 & i = l \\
		0 & {\rm otherwise}
	\end{array} \right.
\] $$

$head w$$
If the argument $italic w$$ is present, it has prototype
$syntax%
	const %Vector% &%w%
%$$
and size $latex m$$.
It specifies the value of $latex w_i$$ in the expression 
for $italic h$$.

$head hes$$
The result $italic hes$$ has prototype
$syntax%
	%Vector% %hes%
%$$
(see $xref/Hessian/Vector/Vector/$$ below)
and its size is $latex n * n$$.
For $latex j = 0 , \ldots , n - 1 $$ 
and $latex \ell = 0 , \ldots , n - 1$$
$latex \[
	hes [ j * n + \ell ] = \DD{ w^{\rm T} F }{ x_j }{ x_\ell } ( x )
\] $$

$head Vector$$
The type $italic Vector$$ must be a $xref/SimpleVector/$$ class with
$xref/SimpleVector/Elements of Specified Type/elements of type/$$
$italic Base$$.
The routine $xref/CheckSimpleVector/$$ will generate an error message
if this is not the case.

$head Hessian Uses Forward$$
After each call to $xref/Forward/$$,
the object $italic f$$ contains the corresponding 
$xref/glossary/Taylor Coefficient/Taylor coefficients/$$.
After $code Hessian$$,
the previous calls to $xref/Forward/$$ are undefined.

$head Example$$
$children%
	example/hessian.cpp%
	example/hes_lagrangian.cpp
%$$
The routines 
$cref/Hessian.cpp/$$ and
$cref/HesLagrangian.cpp/$$
are examples and tests of $code Hessian$$.
They return $code true$$, if they succeed and $code false$$ otherwise.


$end
-----------------------------------------------------------------------------
*/

//  BEGIN CppAD namespace
namespace CppAD {

template <typename Base>
template <typename Vector>
Vector ADFun<Base>::Hessian(const Vector &x, size_t l)
{	size_t i, m = Range();
	CPPAD_ASSERT_KNOWN(
		l < m,
		"Hessian: index i is not less than range dimension for f"
	);

	Vector w(m);
	for(i = 0; i < m; i++)
		w[i] = Base(0);
	w[l] = Base(1);

	return Hessian(x, w);
}


template <typename Base>
template <typename Vector>
Vector ADFun<Base>::Hessian(const Vector &x, const Vector &w)
{	size_t j;
	size_t k;

	size_t n = Domain();

	// check Vector is Simple Vector class with Base type elements
	CheckSimpleVector<Base, Vector>();

	CPPAD_ASSERT_KNOWN(
		x.size() == n,
		"Hessian: length of x not equal domain dimension for f"
	); 
	CPPAD_ASSERT_KNOWN(
		w.size() == Range(),
		"Hessian: length of w not equal range dimension for f"
	);

	// point at which we are evaluating the Hessian
	Forward(0, x);

	// define the return value
	Vector hes(n * n);

	// direction vector for calls to forward
	Vector u(n);
	for(j = 0; j < n; j++)
		u[j] = Base(0);


	// location for return values from Reverse
	Vector ddw(n * 2);

	// loop over forward directions
	for(j = 0; j < n; j++)
	{	// evaluate partials of entire function w.r.t. j-th coordinate
		u[j] = Base(1);
		Forward(1, u);
		u[j] = Base(0);

		// evaluate derivative of partial corresponding to F_i
		ddw = Reverse(2, w);

		// return desired components
		for(k = 0; k < n; k++)
			hes[k * n + j] = ddw[k * 2 + 1];
	}
		
	return hes;
}

} // END CppAD namespace

# endif

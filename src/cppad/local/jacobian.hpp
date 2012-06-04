/* $Id: jacobian.hpp 1598 2009-12-03 02:59:52Z bradbell $ */
# ifndef CPPAD_JACOBIAN_INCLUDED
# define CPPAD_JACOBIAN_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-09 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin Jacobian$$
$spell
	jac
	typename
	Taylor
	Jacobian
	DetLu
	const
$$

$index Jacobian, driver$$
$index first, derivative$$
$index driver, Jacobian$$

$section Jacobian: Driver Routine$$

$head Syntax$$
$syntax%%jac% = %f%.Jacobian(%x%)%$$


$head Purpose$$
We use $latex F : B^n \rightarrow B^m$$ to denote the
$xref/glossary/AD Function/AD function/$$ corresponding to $italic f$$.
The syntax above sets $italic jac$$ to the
Jacobian of $italic F$$ evaluated at $italic x$$; i.e.,
$latex \[
	jac = F^{(1)} (x)
\] $$

$head f$$
The object $italic f$$ has prototype
$syntax%
	ADFun<%Base%> %f%
%$$
Note that the $xref/ADFun/$$ object $italic f$$ is not $code const$$
(see $cref/Forward or Reverse/Jacobian/Forward or Reverse/$$ below).

$head x$$
The argument $italic x$$ has prototype
$syntax%
	const %Vector% &%x%
%$$
(see $xref/Jacobian/Vector/Vector/$$ below)
and its size 
must be equal to $italic n$$, the dimension of the
$xref/seq_property/Domain/domain/$$ space for $italic f$$.
It specifies
that point at which to evaluate the Jacobian.

$head jac$$
The result $italic jac$$ has prototype
$syntax%
	%Vector% %jac%
%$$
(see $xref/Jacobian/Vector/Vector/$$ below)
and its size is $latex m * n$$; i.e., the product of the
$xref/seq_property/Domain/domain/$$
and
$xref/seq_property/Range/range/$$
dimensions for $italic f$$.
For $latex i = 0 , \ldots , m - 1 $$ 
and $latex j = 0 , \ldots , n - 1$$
$latex \[.
	jac[ i * n + j ] = \D{ F_i }{ x_j } ( x )
\] $$


$head Vector$$
The type $italic Vector$$ must be a $xref/SimpleVector/$$ class with
$xref/SimpleVector/Elements of Specified Type/elements of type/$$
$italic Base$$.
The routine $xref/CheckSimpleVector/$$ will generate an error message
if this is not the case.

$head Forward or Reverse$$
This will use order zero Forward mode and either
order one Forward or order one Reverse to compute the Jacobian
(depending on which it estimates will require less work).
After each call to $xref/Forward/$$,
the object $italic f$$ contains the corresponding 
$xref/glossary/Taylor Coefficient/Taylor coefficients/$$.
After each call to $code Jacobian$$,
the previous calls to $xref/Forward/$$ are unspecified.

$head Example$$
$children%
	example/jacobian.cpp
%$$
The routine 
$xref/Jacobian.cpp//Jacobian/$$ is both an example and test.
It returns $code true$$, if it succeeds and $code false$$ otherwise.

$end
-----------------------------------------------------------------------------
*/

//  BEGIN CppAD namespace
namespace CppAD {

template <typename Base, typename Vector>
void JacobianFor(ADFun<Base> &f, const Vector &x, Vector &jac)
{	size_t i;
	size_t j;

	size_t m = f.Domain();
	size_t n = f.Range();

	// check Vector is Simple Vector class with Base type elements
	CheckSimpleVector<Base, Vector>();

	CPPAD_ASSERT_UNKNOWN( x.size()   == f.Domain() );
	CPPAD_ASSERT_UNKNOWN( jac.size() == f.Range() * f.Domain() );

	// argument and result for forward mode calculations
	Vector u(m);
	Vector v(n);

	// initialize all the components
	for(j = 0; j < m; j++)
		u[j] = Base(0);
	
	// loop through the different coordinate directions
	for(j = 0; j < m; j++)
	{	// set u to the j-th coordinate direction
		u[j] = Base(1);

		// compute the partial of f w.r.t. this coordinate direction
		v = f.Forward(1, u);

		// reset u to vector of all zeros
		u[j] = Base(0);

		// return the result
		for(i = 0; i < n; i++)
			jac[ i * m + j ] = v[i];
	}
}
template <typename Base, typename Vector>
void JacobianRev(ADFun<Base> &f, const Vector &x, Vector &jac)
{	size_t i;
	size_t j;

	size_t m = f.Domain();
	size_t n = f.Range();

	CPPAD_ASSERT_UNKNOWN( x.size()   == f.Domain() );
	CPPAD_ASSERT_UNKNOWN( jac.size() == f.Range() * f.Domain() );

	// argument and result for reverse mode calculations
	Vector u(m);
	Vector v(n);

	// initialize all the components
	for(i = 0; i < n; i++)
		v[i] = Base(0);
	
	// loop through the different coordinate directions
	for(i = 0; i < n; i++)
	{	if( f.Parameter(i) )
		{	// return zero for this component of f
			for(j = 0; j < m; j++)
				jac[ i * m + j ] = Base(0);
		}
		else
		{ 
			// set v to the i-th coordinate direction
			v[i] = Base(1);

			// compute the derivative of this component of f
			u = f.Reverse(1, v);

			// reset v to vector of all zeros
			v[i] = Base(0);

			// return the result
			for(j = 0; j < m; j++)
				jac[ i * m + j ] = u[j];
		}
	}
}

template <typename Base>
template <typename Vector>
Vector ADFun<Base>::Jacobian(const Vector &x)
{	size_t i;
	size_t m = Domain();
	size_t n = Range();

	CPPAD_ASSERT_KNOWN(
		x.size() == m,
		"Jacobian: length of x not equal domain dimension for F"
	); 

	// point at which we are evaluating the Jacobian
	Forward(0, x);

	// work factor for forward mode
	size_t workForward = m; 

	// work factor for reverse mode
	size_t workReverse = 0;
	for(i = 0; i < n; i++)
	{	if( ! Parameter(i) )
			++workReverse;
	}

	// choose the method with the least work
	Vector jac( n * m );
	if( workForward <= workReverse )
		JacobianFor(*this, x, jac);
	else	JacobianRev(*this, x, jac);

	return jac;
}

} // END CppAD namespace

# endif

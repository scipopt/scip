/* $Id: sparse_hessian.hpp 2006 2011-07-12 17:35:39Z bradbell $ */
# ifndef CPPAD_SPARSE_HESSIAN_INCLUDED
# define CPPAD_SPARSE_HESSIAN_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin sparse_hessian$$
$spell
	CppAD
	valarray
	std
	Bool
	hes
	const
	Taylor
$$

$section Sparse Hessian: Easy Driver$$
$index SparseHessian$$
$index hessian, sparse$$

$head Syntax$$
$codei%%hes% = %f%.SparseHessian(%x%, %w%)
%$$
$codei%%hes% = %f%.SparseHessian(%x%, %w%, %p%)%$$

$head Purpose$$
We use $latex F : \R^n \rightarrow \R^m$$ do denote the
$cref/AD function/glossary/AD Function/$$
corresponding to $icode f$$. 
The syntax above sets $icode hes$$ to the Hessian 
$latex \[
	hes = \dpow{2}{x} \sum_{i=1}^m w_i F_i (x) 
\] $$
This routine assumes  that the matrix $icode hes$$ is sparse
and uses this assumption sparse to reduce the amount of computation necessary.
One should use speed tests to verify that results are computed faster
than when using the routine $cref/Hessian/$$.

$head f$$
The object $icode f$$ has prototype
$codei%
	ADFun<%Base%> %f%
%$$
Note that the $cref/ADFun/$$ object $icode f$$ is not $code const$$
(see $cref/Uses Forward/sparse_hessian/Uses Forward/$$ below).

$head x$$
The argument $icode x$$ has prototype
$codei%
	const %VectorBase%& %x%
%$$
(see $cref/VectorBase/sparse_hessian/VectorBase/$$ below)
and its size 
must be equal to $icode n$$, the dimension of the
$cref/domain/seq_property/Domain/$$ space for $icode f$$.
It specifies
that point at which to evaluate the Hessian.

$head w$$
The argument $icode w$$ has prototype
$codei%
	const %VectorBase%& %w%
%$$
and size $latex m$$.
It specifies the value of $latex w_i$$ in the expression 
for $icode hes$$.
The more components of $latex w$$ that are identically zero,
the more sparse the resulting Hessian may be (and hence the more efficient
the calculation of $icode hes$$ may be).

$head p$$
The argument $icode p$$ is optional and has prototype
$codei%
	const %VectorSet%& %p%
%$$
(see $cref/VectorSet/sparse_hessian/VectorSet/$$ below)
If it has elements of type $code bool$$,
its size is $latex n * n$$.
If it has elements of type $code std::set<size_t>$$,
its size is $latex n$$ and all its set elements are between
zero and $latex n - 1$$.
It specifies a
$cref/sparsity pattern/glossary/Sparsity Pattern/$$ 
for the Hessian $icode hes$$.
$pre

$$
If this sparsity pattern does not change between calls to 
$codei SparseHessian$$, it should be faster to calculate $icode p$$ once and
pass this argument to $codei SparseHessian$$.
In addition,
if you specify $icode p$$, CppAD will use the same
type of sparsity representation
(vectors of $code bool$$ or vectors of $code std::set<size_t>$$)
for its internal calculations.
Otherwise, the representation
for the internal calculations is unspecified.

$head hes$$
The result $icode hes$$ has prototype
$codei%
	%VectorBase% %hes%
%$$
and its size is $latex n * n$$.
For $latex j = 0 , \ldots , n - 1 $$ 
and $latex \ell = 0 , \ldots , n - 1$$
$latex \[
	hes [ j * n + \ell ] = \DD{ w^{\rm T} F }{ x_j }{ x_\ell } ( x )
\] $$

$head VectorBase$$
The type $icode VectorBase$$ must be a $cref/SimpleVector/$$ class with
$cref/elements of type/SimpleVector/Elements of Specified Type/$$
$icode Base$$.
The routine $cref/CheckSimpleVector/$$ will generate an error message
if this is not the case.

$head VectorSet$$
The type $icode VectorSet$$ must be a $xref/SimpleVector/$$ class with
$xref/SimpleVector/Elements of Specified Type/elements of type/$$
$code bool$$ or $code std::set<size_t>$$;
see $cref/sparsity pattern/glossary/Sparsity Pattern/$$ for a discussion
of the difference.
The routine $xref/CheckSimpleVector/$$ will generate an error message
if this is not the case.

$subhead Restrictions$$
If $icode VectorSet$$ has elements of $code std::set<size_t>$$,
then $icode%p%[%i%]%$$ must return a reference (not a copy) to the 
corresponding set.
According to section 26.3.2.3 of the 1998 C++ standard,
$code std::valarray< std::set<size_t> >$$ does not satisfy
this condition. 

$head Uses Forward$$
After each call to $cref/Forward/$$,
the object $icode f$$ contains the corresponding 
$cref/Taylor coefficients/glossary/Taylor Coefficient/$$.
After $code SparseHessian$$,
the previous calls to $xref/Forward/$$ are undefined.

$head Example$$
$children%
	example/sparse_hessian.cpp
%$$
The routine
$cref/sparse_hessian.cpp/$$ 
is examples and tests of $code sparse_hessian$$.
It return $code true$$, if it succeeds and $code false$$ otherwise.


$end
-----------------------------------------------------------------------------
*/
# include <cppad/local/std_set.hpp>

CPPAD_BEGIN_NAMESPACE
/*!
\file sparse_hessian.hpp
Sparse Hessian driver routine and helper functions.
*/

/*!
Private helper function for SparseHessian(x, w, p).

All of the description in the public member function SparseHessian(x, w, p)
applies.

\param set_type
is a \c bool value. This argument is used to dispatch to the proper source
code depending on the value of \c VectorSet::value_type.

\param x
See \c SparseHessian(x, w, p).

\param w
See \c SparseHessian(x, w, p).

\param p
See \c SparseHessian(x, w, p).

\param hes
is the return value for the corresponding call to \c SparseHessian(x, w, p).
On input, it must have size equal to the domain times range dimension
for this ADFun<Base> object.
On return, it will contain the Hessian.
*/

template <class Base>
template <class VectorBase, class VectorSet>
void ADFun<Base>::SparseHessianCase(
	bool               set_type         ,
	const VectorBase&  x                ,
	const VectorBase&  w                ,
	const VectorSet&   p                ,
	VectorBase&        hes              )
{
	typedef CppAD::vector<size_t> SizeVector;
	typedef CppAD::vectorBool     VectorBool;
	size_t j, k, l;

	size_t n = Domain();

	// check Vector is Simple Vector class with bool elements
	CheckSimpleVector<bool, VectorSet>();

	// check Vector is Simple Vector class with Base type elements
	CheckSimpleVector<Base, VectorBase>();

	CPPAD_ASSERT_KNOWN(
		x.size() == n,
		"SparseHessian: size of x not equal domain dimension for f"
	); 
	CPPAD_ASSERT_KNOWN(
		w.size() == Range(),
		"SparseHessian: size of w not equal range dimension for f"
	);
	CPPAD_ASSERT_KNOWN(
		p.size() == n * n,
		"SparseHessian: using bool values and size of p "
		"not equal square of domain dimension for f"
	);
	
	// initial coloring
	SizeVector color(n);
	for(j = 0; j < n; j++)
		color[j] = j;

	// See GreedyPartialD2Coloring Algorithm Section 3.6.2 of
	// Graph Coloring in Optimization Revisited by
	// Assefaw Gebremedhin, Fredrik Maane, Alex Pothen
	VectorBool forbidden(n);
	for(j = 0; j < n; j++)
	{	// initial all colors as ok for this column
		for(k = 0; k < n; k++)
			forbidden[k] = false;
		// for each row that is connected to column j
		for(l = 0; l < n; l++) if( p[l * n + j] )
		{	// for each column that is connected to row l
			for(k = 0; k < n; k++) if( p[l * n + k] & (j != k) )	
				forbidden[ color[k] ] = true;
		}
		k = 0;
		while( forbidden[k] && k < n )
		{	k++;
			CPPAD_ASSERT_UNKNOWN( k < n );
		}
		color[j] = k;
	}
	size_t n_color = 1;
	for(k = 0; k < n; k++) n_color = std::max(n_color, color[k] + 1);

	// some values
	const Base zero(0);
	const Base one(1);

	// point at which we are evaluating the Hessian
	Forward(0, x);

	// initialize the return value
	for(j = 0; j < n; j++)
		for(k = 0; k < n; k++)
			hes[j * n + k] = zero;

	// direction vector for calls to forward
	VectorBase u(n);

	// location for return values from Reverse
	VectorBase ddw(n * 2);

	// loop over colors
	size_t c;
	for(c = 0; c < n_color; c++)
	{	// determine all the colums with this color
		for(j = 0; j < n; j++)
		{	if( color[j] == c )
				u[j] = one;
			else	u[j] = zero;
		}
		// call forward mode for all these columns at once
		Forward(1, u);

		// evaluate derivative of w^T * F'(x) * u
		ddw = Reverse(2, w);

		// set the corresponding components of the result
		for(j = 0; j < n; j++) if( color[j] == c )
		{	for(l = 0; l < n; l++) 
				if( p[ l * n + j ] )
					hes[l * n + j] = ddw[l * 2 + 1];
		}
	}
}
/*!
Private helper function for SparseHessian(x, w, p).

All of the description in the public member function SparseHessian(x, w, p)
applies.

\param set_type
is a \c std::set<size_t> value.
This argument is used to dispatch to the proper source
code depending on the value of \c VectorSet::value_type.

\param x
See \c SparseHessian(x, w, p).

\param w
See \c SparseHessian(x, w, p).

\param p
See \c SparseHessian(x, w, p).

\param hes
is the return value for the corresponding call to \c SparseHessian(x, w, p).
On input, it must have size equal to the domain times range dimension
for this ADFun<Base> object.
On return, it will contain the Hessian.
*/

template <class Base>
template <class VectorBase, class VectorSet>
void ADFun<Base>::SparseHessianCase(
	const std::set<size_t>&  set_type         ,
	const VectorBase&        x                ,
	const VectorBase&        w                ,
	const VectorSet&         p                ,
	VectorBase&              hes              )
{
	typedef CppAD::vector<size_t> SizeVector;
	typedef CppAD::vectorBool     VectorBool;
	size_t j, k, l;

	size_t n = Domain();

	// check VectorSet is Simple Vector class with sets for elements
	CheckSimpleVector<std::set<size_t>, VectorSet>(
		one_element_std_set<size_t>(), two_element_std_set<size_t>()
	);

	// check Vector is Simple Vector class with Base type elements
	CheckSimpleVector<Base, VectorBase>();

	CPPAD_ASSERT_KNOWN(
		x.size() == n,
		"SparseHessian: size of x not equal domain dimension for f"
	); 
	CPPAD_ASSERT_KNOWN(
		w.size() == Range(),
		"SparseHessian: size of w not equal range dimension for f"
	);
	CPPAD_ASSERT_KNOWN(
		p.size() == n,
		"SparseHessian: using set values and size of p "
		"not equal domain dimension for f"
	);
	
	// initial coloring
	SizeVector color(n);
	for(j = 0; j < n; j++)
		color[j] = j;

	// See GreedyPartialD2Coloring Algorithm Section 3.6.2 of
	// Graph Coloring in Optimization Revisited by
	// Assefaw Gebremedhin, Fredrik Maane, Alex Pothen
	VectorBool forbidden(n);
	std::set<size_t>::const_iterator itr_k, itr_l;
	for(j = 0; j < n; j++)
	{	// initial all colors as ok for this column
		for(k = 0; k < n; k++)
			forbidden[k] = false;

		// for each row that is connected to column j
		itr_l = p[j].begin();
		while( itr_l != p[j].end() )
		{	l = *itr_l++;
			// for each column that is connected to row l
			// (same as each row that is connect to column l)
			itr_k = p[l].begin();
			while( itr_k != p[l].end() )
			{	k = *itr_k++;
				forbidden[ color[k] ] = (j != k);
			}
		}
		k = 0;
		while( forbidden[k] && k < n )
		{	k++;
			CPPAD_ASSERT_UNKNOWN( k < n );
		}
		color[j] = k;
	}
	size_t n_color = 1;
	for(k = 0; k < n; k++) n_color = std::max(n_color, color[k] + 1);

	// some values
	const Base zero(0);
	const Base one(1);

	// point at which we are evaluating the Hessian
	Forward(0, x);

	// initialize the return value
	for(j = 0; j < n; j++)
		for(k = 0; k < n; k++)
			hes[j * n + k] = zero;

	// direction vector for calls to forward
	VectorBase u(n);

	// location for return values from Reverse
	VectorBase ddw(n * 2);

	// loop over colors
	size_t c;
	for(c = 0; c < n_color; c++)
	{	// determine all the colums with this color
		for(j = 0; j < n; j++)
		{	if( color[j] == c )
				u[j] = one;
			else	u[j] = zero;
		}
		// call forward mode for all these columns at once
		Forward(1, u);

		// evaluate derivative of w^T * F'(x) * u
		ddw = Reverse(2, w);

		// set the corresponding components of the result
		for(j = 0; j < n; j++) if( color[j] == c )
		{	itr_l = p[j].begin();
			while( itr_l != p[j].end() )
			{	l = *itr_l++;
				hes[l * n + j] = ddw[l * 2 + 1];
			}
		}
	}
}

/*!
Compute a sparse Hessian.

The C++ source code coresponding to this operation is
\verbatim
	hes = SparseHessian(x, w)
\endverbatim


\tparam Base
is the base type for the recording that is stored in this
ADFun<Base object.

\tparam VectorBase
is a simple vector class with elements of the \a Base.

\tparam VectorSet
is a simple vector class with elements of type
\c bool or \c std::set<size_t>.

\param x
is a vector specifing the point at which to compute the Hessian.

\param w
The Hessian is computed for a weighted sum of the components
of the function corresponding to this ADFun<Base> object.
The argument \a w specifies the weights for each component.
It must have size equal to the range dimension for this ADFun<Base> object.

\param p
is a sparsity pattern for the Hessian.

\return
Will be a vector of size \c n * n containing the Hessian of 
at the point specified by \a x
(where \c n is the domain dimension for this ADFun<Base> object).
*/
template <class Base>
template <class VectorBase, class VectorSet>
VectorBase ADFun<Base>::SparseHessian(
	const VectorBase& x, const VectorBase& w, const VectorSet& p
)
{	size_t n = Domain();
	VectorBase hes(n * n);

	typedef typename VectorSet::value_type Set_type;
	SparseHessianCase(Set_type(), x, w, p, hes);

	return hes;
}

/*!
Compute a sparse Hessian

The C++ source code coresponding to this operation is
\verbatim
	hes = SparseHessian(x, w)
\endverbatim


\tparam Base
is the base type for the recording that is stored in this
ADFun<Base object.

\tparam VectorBase
is a simple vector class with elements of the \a Base.

\param x
is a vector specifing the point at which to compute the Hessian.

\param w
The Hessian is computed for a weighted sum of the components
of the function corresponding to this ADFun<Base> object.
The argument \a w specifies the weights for each component.
It must have size equal to the range dimension for this ADFun<Base> object.

\return
Will be a vector of size \c n * n containing the Hessian of 
at the point specified by \a x
(where \c n is the domain dimension for this ADFun<Base> object).
*/
template <class Base>
template <class VectorBase>
VectorBase ADFun<Base>::SparseHessian(const VectorBase &x, const VectorBase &w)
{	size_t i, j, k;
	typedef CppAD::vectorBool VectorBool;

	size_t m = Range();
	size_t n = Domain();

	// determine the sparsity pattern p for Hessian of w^T F
	VectorBool r(n * n);
	for(j = 0; j < n; j++)
	{	for(k = 0; k < n; k++)
			r[j * n + k] = false;
		r[j * n + j] = true;
	}
	ForSparseJac(n, r);
	//
	VectorBool s(m);
	for(i = 0; i < m; i++)
		s[i] = w[i] != 0;
	VectorBool p = RevSparseHes(n, s);

	// compute sparse Hessian
	VectorBase hes(n * n);
	bool set_type = true; // only used to dispatch compiler to proper case
	SparseHessianCase(set_type, x, w, p, hes);

	return hes;
}

CPPAD_END_NAMESPACE
# endif

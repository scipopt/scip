/* $Id: sparse_jacobian.hpp 2006 2011-07-12 17:35:39Z bradbell $ */
# ifndef CPPAD_SPARSE_JACOBIAN_INCLUDED
# define CPPAD_SPARSE_JACOBIAN_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin sparse_jacobian$$
$spell
	valarray
	std
	CppAD
	Bool
	jac
	Jacobian
	const
	Taylor
$$

$section Sparse Jacobian: Easy Driver$$
$index SparseJacobian$$
$index jacobian, sparse$$

$head Syntax$$
$codei%%jac% = %f%.SparseJacobian(%x%)
%$$
$codei%%jac% = %f%.SparseJacobian(%x%, %p%)%$$

$head Purpose$$
We use $latex F : \R^n \rightarrow \R^m$$ do denote the
$cref/AD function/glossary/AD Function/$$
corresponding to $icode f$$. 
The syntax above sets $icode jac$$ to the Jacobian 
$latex \[
	jac = F^{(1)} (x) 
\] $$
This routine assumes
that the matrix $latex F^{(1)} (x) \in \R^{m \times n}$$ is sparse
and uses this assumption to reduce the amount of computation necessary.
One should use speed tests (e.g. $cref/speed_test/$$)
to verify that results are computed faster
than when using the routine $cref/Jacobian/$$.

$head f$$
The object $icode f$$ has prototype
$codei%
	ADFun<%Base%> %f%
%$$
Note that the $cref/ADFun/$$ object $icode f$$ is not $code const$$
(see $cref/Uses Forward/sparse_jacobian/Uses Forward/$$ below).

$head x$$
The argument $icode x$$ has prototype
$codei%
	const %VectorBase%& %x%
%$$
(see $cref/VectorBase/sparse_jacobian/VectorBase/$$ below)
and its size 
must be equal to $icode n$$, the dimension of the
$cref/domain/seq_property/Domain/$$ space for $icode f$$.
It specifies
that point at which to evaluate the Jacobian.

$head p$$
The argument $icode p$$ is optional and has prototype
$syntax%
	const %VectorSet%& %p%
%$$
(see $cref/VectorSet/sparse_jacobian/VectorSet/$$ below).
If it has elements of type $code bool$$,
its size is $latex m * n$$.
If it has elements of type $code std::set<size_t>$$,
its size is $latex m$$ and all its set elements are between
zero and $latex n - 1$$.
It specifies a 
$cref/sparsity pattern/glossary/Sparsity Pattern/$$ 
for the Jacobian $latex F^{(1)} (x)$$.
$pre

$$
If this sparsity pattern does not change between calls to 
$codei SparseJacobian$$, it should be faster to calculate $icode p$$ once 
(using $cref/ForSparseJac/$$ or $cref/RevSparseJac/$$)
and then pass $icode p$$ to $codei SparseJacobian$$.
In addition,
if you specify $icode p$$, CppAD will use the same
type of sparsity representation 
(vectors of $code bool$$ or vectors of $code std::set<size_t>$$)
for its internal calculations.
Otherwise, the representation
for the internal calculations is unspecified.

$head jac$$
The result $icode jac$$ has prototype
$codei%
	%VectorBase% %jac%
%$$
and its size is $latex m * n$$.
For $latex i = 0 , \ldots , m - 1$$,
and $latex j = 0 , \ldots , n - 1 $$ 
$latex \[
	jac [ i * n + j ] = \D{ F_i }{ x_j } (x)
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
The routine $cref/CheckSimpleVector/$$ will generate an error message
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
After $code SparseJacobian$$,
the previous calls to $xref/Forward/$$ are undefined.

$head Example$$
$children%
	example/sparse_jacobian.cpp
%$$
The routine
$cref/sparse_jacobian.cpp/$$
is examples and tests of $code sparse_jacobian$$.
It return $code true$$, if it succeeds and $code false$$ otherwise.

$end
-----------------------------------------------------------------------------
*/
# include <cppad/local/std_set.hpp>

CPPAD_BEGIN_NAMESPACE
/*!
\file sparse_jacobian.hpp
Sparse Jacobian driver routine and helper functions.
*/

/*!
Private helper function for SparseJacobian(x, p).

All of the description in the public member function SparseJacobian(x, p)
applies.

\param set_type
is a \c bool value. This argument is used to dispatch to the proper souce
code depending on the value of \c VectorSet::value_type.

\param x
See \c SparseJacobian(x, p).

\param p
See \c SparseJacobian(x, p).

\param jac
is the return value for the corresponding call to \c SparseJacobian(x, p).
On input, it must have size equal to the domain times range dimension
for this ADFun<Base> object.
On return, it will contain the Jacobian.
*/
template <class Base>
template <class VectorBase, class VectorSet>
void ADFun<Base>::SparseJacobianCase(
	bool               set_type        ,
	const VectorBase&  x               , 
	const VectorSet&   p               ,
	VectorBase&        jac             )
{
	typedef CppAD::vector<size_t> SizeVector;
	typedef CppAD::vectorBool     VectorBool;
	size_t i, j, k;

	size_t m = Range();
	size_t n = Domain();

	// some values
	const Base zero(0);
	const Base one(1);

	// check VectorSet is Simple Vector class with bool elements
	CheckSimpleVector<bool, VectorSet>();

	// check VectorBase is Simple Vector class with Base type elements
	CheckSimpleVector<Base, VectorBase>();

	CPPAD_ASSERT_KNOWN(
		x.size() == n,
		"SparseJacobian: size of x not equal domain dimension for f"
	); 
	CPPAD_ASSERT_KNOWN(
		p.size() == m * n,
		"SparseJacobian: using bool values and size of p "
		" not equal range dimension times domain dimension for f"
	); 
	CPPAD_ASSERT_UNKNOWN(jac.size() == m * n); 

	// point at which we are evaluating the Jacobian
	Forward(0, x);

	// initialize the return value
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			jac[i * n + j] = zero;

	if( n <= m )
	{	// use forward mode ----------------------------------------
	
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
			for(i = 0; i < m; i++) if( p[i * n + j] )
			{	// for each column that is connected to row i
				for(k = 0; k < n; k++)
					if( p[i * n + k] & (j != k) )	
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
		for(k = 0; k < n; k++) 
			n_color = std::max(n_color, color[k] + 1);

		// direction vector for calls to forward
		VectorBase dx(n);

		// location for return values from Reverse
		VectorBase dy(m);

		// loop over colors
		size_t c;
		for(c = 0; c < n_color; c++)
		{	// determine all the colums with this color
			for(j = 0; j < n; j++)
			{	if( color[j] == c )
					dx[j] = one;
				else	dx[j] = zero;
			}
			// call forward mode for all these columns at once
			dy = Forward(1, dx);

			// set the corresponding components of the result
			for(j = 0; j < n; j++) if( color[j] == c )
			{	for(i = 0; i < m; i++) 
					if( p[ i * n + j ] )
						jac[i * n + j] = dy[i];
			}
		}
	}
	else
	{	// use reverse mode ----------------------------------------
	
		// initial coloring
		SizeVector color(m);
		for(i = 0; i < m; i++)
			color[i] = i;

		// See GreedyPartialD2Coloring Algorithm Section 3.6.2 of
		// Graph Coloring in Optimization Revisited by
		// Assefaw Gebremedhin, Fredrik Maane, Alex Pothen
		VectorBool forbidden(m);
		for(i = 0; i < m; i++)
		{	// initial all colors as ok for this row
			for(k = 0; k < m; k++)
				forbidden[k] = false;
			// for each column that is connected to row i
			for(j = 0; j < n; j++) if( p[i * n + j] )
			{	// for each row that is connected to column j
				for(k = 0; k < m; k++)
					if( p[k * n + j] & (i != k) )	
						forbidden[ color[k] ] = true;
			}
			k = 0;
			while( forbidden[k] && k < m )
			{	k++;
				CPPAD_ASSERT_UNKNOWN( k < n );
			}
			color[i] = k;
		}
		size_t n_color = 1;
		for(k = 0; k < m; k++) 
			n_color = std::max(n_color, color[k] + 1);

		// weight vector for calls to reverse
		VectorBase w(m);

		// location for return values from Reverse
		VectorBase dw(n);

		// loop over colors
		size_t c;
		for(c = 0; c < n_color; c++)
		{	// determine all the rows with this color
			for(i = 0; i < m; i++)
			{	if( color[i] == c )
					w[i] = one;
				else	w[i] = zero;
			}
			// call reverse mode for all these rows at once
			dw = Reverse(1, w);

			// set the corresponding components of the result
			for(i = 0; i < m; i++) if( color[i] == c )
			{	for(j = 0; j < n; j++) 
					if( p[ i * n + j ] )
						jac[i * n + j] = dw[j];
			}
		}
	}
}

/*!
Private helper function for SparseJacobian(x, p).

All of the description in the public member function SparseJacobian(x, p)
applies.

\param set_type
is a \c std::set<size_t> value.
This argument is used to dispatch to the proper souce
code depending on the vlaue of \c VectorSet::value_type.

\param x
See \c SparseJacobian(x, p).

\param p
See \c SparseJacobian(x, p).

\param jac
is the return value for the corresponding call to \c SparseJacobian(x, p).
On input it must have size equalt to the domain times range dimension
for this ADFun<Base> object.
On return, it will contain the Jacobian.
*/
template <class Base>
template <class VectorBase, class VectorSet>
void ADFun<Base>::SparseJacobianCase(
	const std::set<size_t>&  set_type        ,
	const VectorBase&        x               , 
	const VectorSet&         p               ,
	VectorBase&              jac             )
{
	typedef CppAD::vector<size_t> SizeVector;
	typedef CppAD::vectorBool     VectorBool;
	size_t i, j, k;

	size_t m = Range();
	size_t n = Domain();

	// some values
	const Base zero(0);
	const Base one(1);

	// check VectorSet is Simple Vector class with sets for elements
	CheckSimpleVector<std::set<size_t>, VectorSet>(
		one_element_std_set<size_t>(), two_element_std_set<size_t>()
	);

	// check VectorBase is Simple Vector class with Base type elements
	CheckSimpleVector<Base, VectorBase>();

	CPPAD_ASSERT_KNOWN(
		x.size() == n,
		"SparseJacobian: size of x not equal domain dimension for f"
	); 
	CPPAD_ASSERT_KNOWN(
		p.size() == m,
		"SparseJacobian: using sets and size of p "
		"not equal range dimension for f"
	); 
	CPPAD_ASSERT_UNKNOWN(jac.size() == m * n); 

	// point at which we are evaluating the Jacobian
	Forward(0, x);

	// initialize the return value
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			jac[i * n + j] = zero;

	// create a copy of the transpose sparsity pattern
	VectorSet q(n);
	std::set<size_t>::const_iterator itr_i, itr_j;
	for(i = 0; i < m; i++)
	{	itr_j = p[i].begin();
		while( itr_j != p[i].end() )
		{	j = *itr_j++;
			q[j].insert(i);
		}
	}	

	if( n <= m )
	{	// use forward mode ----------------------------------------
	
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

			// for each row connected to column j
			itr_i = q[j].begin();
			while( itr_i != q[j].end() )
			{	i = *itr_i++;
				// for each column connected to row i
				itr_j = p[i].begin();
				while( itr_j != p[i].end() )
				{	// if this is not j, forbid it
					k = *itr_j++;
					forbidden[ color[k] ] = (k != j);
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
		for(k = 0; k < n; k++) 
			n_color = std::max(n_color, color[k] + 1);

		// direction vector for calls to forward
		VectorBase dx(n);

		// location for return values from Reverse
		VectorBase dy(m);

		// loop over colors
		size_t c;
		for(c = 0; c < n_color; c++)
		{	// determine all the colums with this color
			for(j = 0; j < n; j++)
			{	if( color[j] == c )
					dx[j] = one;
				else	dx[j] = zero;
			}
			// call forward mode for all these columns at once
			dy = Forward(1, dx);

			// set the corresponding components of the result
			for(j = 0; j < n; j++) if( color[j] == c )
			{	itr_i = q[j].begin();
				while( itr_i != q[j].end() )
				{	i = *itr_i++;
					jac[i * n + j] = dy[i];
				}
			}
		}
	}
	else
	{	// use reverse mode ----------------------------------------
	
		// initial coloring
		SizeVector color(m);
		for(i = 0; i < m; i++)
			color[i] = i;

		// See GreedyPartialD2Coloring Algorithm Section 3.6.2 of
		// Graph Coloring in Optimization Revisited by
		// Assefaw Gebremedhin, Fredrik Maane, Alex Pothen
		VectorBool forbidden(m);
		for(i = 0; i < m; i++)
		{	// initial all colors as ok for this row
			for(k = 0; k < m; k++)
				forbidden[k] = false;

			// for each column connected to row i
			itr_j = p[i].begin();
			while( itr_j != p[i].end() )
			{	j = *itr_j++;	
				// for each row connected to column j
				itr_i = q[j].begin();
				while( itr_i != q[j].end() )
				{	// if this is not i, forbid it
					k = *itr_i++;
					forbidden[ color[k] ] = (k != i);
				}
			}
			k = 0;
			while( forbidden[k] && k < m )
			{	k++;
				CPPAD_ASSERT_UNKNOWN( k < n );
			}
			color[i] = k;
		}
		size_t n_color = 1;
		for(k = 0; k < m; k++) 
			n_color = std::max(n_color, color[k] + 1);

		// weight vector for calls to reverse
		VectorBase w(m);

		// location for return values from Reverse
		VectorBase dw(n);

		// loop over colors
		size_t c;
		for(c = 0; c < n_color; c++)
		{	// determine all the rows with this color
			for(i = 0; i < m; i++)
			{	if( color[i] == c )
					w[i] = one;
				else	w[i] = zero;
			}
			// call reverse mode for all these rows at once
			dw = Reverse(1, w);

			// set the corresponding components of the result
			for(i = 0; i < m; i++) if( color[i] == c )
			{	itr_j = p[i].begin();
				while( itr_j != p[i].end() )
				{	j = *itr_j++;
					jac[i * n + j] = dw[j];
				}
			}
		}
	}
}

/*!
Compute a sparse Jacobian.

The C++ source code corresponding to this operation is
\verbatim
	jac = SparseJacobian(x, p)
\endverbatim

\tparam Base
is the base type for the recording that is stored in this
ADFun<Base object.

\tparam VectorBase
is a simple vector class with elements of type \a Base.

\tparam VectorSet
is a simple vector class with elements of type 
\c bool or \c std::set<size_t>.

\param x
is a vector specifing the point at which to compute the Jacobian.

\param p
is the sparsity pattern for the Jacobian that we are calculating.

\return
Will be a vector if size \c m * n containing the Jacobian at the
specified point (in row major order).
*/
template <class Base>
template <class VectorBase, class VectorSet>
VectorBase ADFun<Base>::SparseJacobian(
	const VectorBase& x, const VectorSet& p
)
{	size_t m = Range();
	size_t n = Domain();
	VectorBase jac(m * n);

	typedef typename VectorSet::value_type Set_type;
	SparseJacobianCase( Set_type(), x, p, jac);

	return jac;
} 

/*!
Compute a sparse Jacobian.

The C++ source code corresponding to this operation is
\verbatim
	jac = SparseJacobian(x)
\endverbatim

\tparam Base
is the base type for the recording that is stored in this
ADFun<Base object.

\tparam VectorBase
is a simple vector class with elements of the \a Base.

\param x
is a vector specifing the point at which to compute the Jacobian.

\return
Will be a vector of size \c m * n containing the Jacobian at the
specified point (in row major order).
*/
template <class Base>
template <class VectorBase>
VectorBase ADFun<Base>::SparseJacobian( const VectorBase& x )
{	typedef CppAD::vectorBool   VectorBool;

	size_t m = Range();
	size_t n = Domain();

	// sparsity pattern for Jacobian
	VectorBool p(m * n);

	// return result
	VectorBase jac(m * n);

	if( n <= m )
	{	size_t j, k;

		// use forward mode 
		VectorBool r(n * n);
		for(j = 0; j < n; j++)
		{	for(k = 0; k < n; k++)
				r[j * n + k] = false;
			r[j * n + j] = true;
		}
		p = ForSparseJac(n, r);
	}
	else
	{	size_t i, k;

		// use reverse mode 
		VectorBool s(m * m);
		for(i = 0; i < m; i++)
		{	for(k = 0; k < m; k++)
				s[i * m + k] = false;
			s[i * m + i] = true;
		}
		p = RevSparseJac(m, s);
	}
	bool set_type = true; // only used to dispatch compiler to proper case
	SparseJacobianCase(set_type, x, p, jac);
	return jac;
}


CPPAD_END_NAMESPACE
# endif

/* $Id: rev_sparse_jac.hpp 2006 2011-07-12 17:35:39Z bradbell $ */
# ifndef CPPAD_REV_SPARSE_JAC_INCLUDED
# define CPPAD_REV_SPARSE_JAC_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin RevSparseJac$$
$spell
	std
	VecAD
	var
	Jacobian
	Jac
	const
	Bool
	Dep
	proportional
$$

$section Jacobian Sparsity Pattern: Reverse Mode$$ 

$index RevSparseJac$$
$index reverse, sparse Jacobian$$
$index sparse, reverse Jacobian$$
$index pattern, reverse Jacobian$$

$head Syntax$$
$icode%r% = %F%.RevSparseJac(%p%, %s%)%$$

$head Purpose$$
We use $latex F : \R^n \rightarrow R^m$$ to denote the
$xref/glossary/AD Function/AD function/$$ corresponding to $icode f$$.
For a fixed $latex p \times m$$ matrix $latex S$$,
the Jacobian of $latex S * F( x )$$
with respect to $latex x$$ is
$latex \[
	J(x) = S * F^{(1)} ( x )
\] $$
Given a
$xref/glossary/Sparsity Pattern/sparsity pattern/$$ 
for $latex S$$,
$code RevSparseJac$$ returns a sparsity pattern for the $latex J(x)$$.

$head f$$
The object $icode f$$ has prototype
$codei%
	ADFun<%Base%> %f%
%$$

$head x$$
the sparsity pattern is valid for all values of the independent 
variables in $latex x \in \R^n$$
(even if it has $cref/CondExp/$$ or $cref/VecAD/$$ operations).

$head p$$
The argument $icode p$$ has prototype
$codei%
	size_t %p%
%$$
It specifies the number of rows in
$latex S \in \R^{p \times m}$$ and the 
Jacobian $latex J(x) \in \R^{p \times n}$$. 

$head s$$
The argument $icode s$$ has prototype
$codei%
	const %VectorSet%& %s%
%$$
(see $xref/RevSparseJac/VectorSet/VectorSet/$$ below).
If it has elements of type $code bool$$,
its size is $latex p * m$$.
If it has elements of type $code std::set<size_t>$$,
its size is $icode p$$ and all its set elements are between
zero and $latex m - 1$$.
It specifies a 
$xref/glossary/Sparsity Pattern/sparsity pattern/$$ 
for the matrix $icode S$$.

$head r$$
The return value $icode r$$ has prototype
$codei%
	%VectorSet% %r%
%$$
(see $xref/RevSparseJac/VectorSet/VectorSet/$$ below).
If it has elements of type $code bool$$,
its size is $latex p * n$$.
If it has elements of type $code std::set<size_t>$$,
its size is $icode p$$.
It specifies a 
$xref/glossary/Sparsity Pattern/sparsity pattern/$$ 
for the matrix $latex J(x)$$.

$head VectorSet$$
The type $icode VectorSet$$ must be a $xref/SimpleVector/$$ class with
$xref/SimpleVector/Elements of Specified Type/elements of type/$$
$code bool$$ or $code std::set<size_t>$$;
see $cref/sparsity pattern/glossary/Sparsity Pattern/$$ for a discussion
of the difference.

$head Entire Sparsity Pattern$$
Suppose that $latex p = m$$ and
$latex S$$ is the $latex m \times m$$ identity matrix.
In this case, 
the corresponding value for $icode r$$ is a 
sparsity pattern for the Jacobian $latex J(x) = F^{(1)} ( x )$$.

$head Example$$
$children%
	example/rev_sparse_jac.cpp
%$$
The file
$xref/RevSparseJac.cpp/$$
contains an example and test of this operation.
It returns true if it succeeds and false otherwise.

$end
-----------------------------------------------------------------------------
*/

# include <cppad/local/std_set.hpp>

CPPAD_BEGIN_NAMESPACE
/*!
\file rev_sparse_jac.hpp
Reverse mode Jacobian sparsity patterns.
*/

/*!
Calculate Jacobian sparsity patterns using reverse mode.

The C++ source code corresponding to this operation is
\verbatim
	s = f.RevSparseJac(q, r)
\endverbatim

\tparam Base
is the base type for this recording.

\tparam VectorSet
is a simple vector class with elements of type \c bool.

\param p
is the number of rows in the matrix \f$ S \f$.

\param s
is a sparsity pattern for the matrix \f$ S \f$.

\param r
the input value of \a r must be a vector with size \c p*n
where \c n is the number of independent variables
corresponding to the operation sequence stored in \a play. 
The input value of the components of \c r does not matter.
On output, \a r is the sparsity pattern for the matrix
\f[
	J(x) = S * F^{(1)} (x)
\f]
where \f$ F \f$ is the function corresponding to the operation sequence
and \a x is any argument value.

\param total_num_var
is the total number of variable in this recording.

\param dep_taddr
maps dependendent variable index
to the corresponding variable in the tape.

\param ind_taddr
maps independent variable index
to the corresponding variable in the tape.

\param play
is the recording that defines the function we are computing the sparsity 
pattern for.
*/

template <class Base, class VectorSet> 
void RevSparseJacBool(
	size_t                 p                , 
	const VectorSet&       s                ,
	VectorSet&             r                ,
	size_t                 total_num_var    ,
	CppAD::vector<size_t>& dep_taddr        ,
	CppAD::vector<size_t>& ind_taddr        ,
	CppAD::player<Base>&   play             )
{
	// temporary indices
	size_t i, j;

	// check VectorSet is Simple Vector class with bool elements
	CheckSimpleVector<bool, VectorSet>();

	// range and domain dimensions for F
	size_t m = dep_taddr.size();
	size_t n = ind_taddr.size();

	CPPAD_ASSERT_KNOWN(
		p > 0,
		"RevSparseJac: p (first argument) is not greater than zero"
	);

	CPPAD_ASSERT_KNOWN(
		s.size() == p * m,
		"RevSparseJac: s (second argument) length is not equal to\n"
		"p (first argument) times range dimension for ADFun object."
	);

	// vector of sets that will hold the results
	sparse_pack    var_sparsity;
	var_sparsity.resize(total_num_var, p);

	// The sparsity pattern corresponding to the dependent variables
	for(i = 0; i < m; i++)
	{	CPPAD_ASSERT_UNKNOWN( dep_taddr[i] < total_num_var );

		for(j = 0; j < p; j++) if( s[ i * m + j ] )
			var_sparsity.add_element( dep_taddr[i], j );
	}

	// evaluate the sparsity patterns
	RevJacSweep(
		n,
		total_num_var,
		&play,
		var_sparsity
	);

	// return values corresponding to dependent variables
	CPPAD_ASSERT_UNKNOWN( r.size() == p * n );
	for(j = 0; j < n; j++)
	{	CPPAD_ASSERT_UNKNOWN( ind_taddr[j] == (j+1) );

		// ind_taddr[j] is operator taddr for j-th independent variable
		CPPAD_ASSERT_UNKNOWN( play.GetOp( ind_taddr[j] ) == InvOp );

		// extract the result from var_sparsity
		for(i = 0; i < p; i++) 
			r[ i * n + j ] = false;
		CPPAD_ASSERT_UNKNOWN( var_sparsity.end() == p );
		var_sparsity.begin(j+1);
		i = var_sparsity.next_element();
		while( i < p )
		{	r[ i * n + j ] = true;
			i              = var_sparsity.next_element();
		}
	}
}

/*!
\file rev_sparse_jac.hpp
Reverse mode Jacobian sparsity patterns.
*/

/*!
Calculate Jacobian sparsity patterns using reverse mode.

The C++ source code corresponding to this operation is
\verbatim
	s = f.RevSparseJac(q, r)
\endverbatim

\tparam Base
is the base type for this recording.

\tparam VectorSet
is a simple vector class with elements of type \c std::set<size_t>.

\param p
is the number of rows in the matrix \f$ S \f$.

\param s
is a sparsity pattern for the matrix \f$ S \f$.

\param r
the input value of \a r must be a vector with size \a p.
On input, each element of \a r mus be an empty set.
On output, \a r is the sparsity pattern for the matrix
\f[
	J(x) = S * F^{(1)} (x)
\f]
where \f$ F \f$ is the function corresponding to the operation sequence
and \a x is any argument value.

\param total_num_var
is the total number of variable in this recording.

\param dep_taddr
maps dependendent variable index
to the corresponding variable in the tape.

\param ind_taddr
maps independent variable index
to the corresponding variable in the tape.

\param play
is the recording that defines the function we are computing the sparsity 
pattern for.
*/

template <class Base, class VectorSet> 
void RevSparseJacSet(
	size_t                 p                , 
	const VectorSet&       s                ,
	VectorSet&             r                ,
	size_t                 total_num_var    ,
	CppAD::vector<size_t>& dep_taddr        ,
	CppAD::vector<size_t>& ind_taddr        ,
	CppAD::player<Base>&   play             )
{
	// temporary indices
	size_t i, j;
	std::set<size_t>::const_iterator itr;

	// check VectorSet is Simple Vector class with sets for elements
	CheckSimpleVector<std::set<size_t>, VectorSet>(
		one_element_std_set<size_t>(), two_element_std_set<size_t>()
	);

	// range and domain dimensions for F
	size_t m = dep_taddr.size();
	size_t n = ind_taddr.size();

	CPPAD_ASSERT_KNOWN(
		p > 0,
		"RevSparseJac: p (first argument) is not greater than zero"
	);

	CPPAD_ASSERT_KNOWN(
		s.size() == p,
		"RevSparseJac: s (second argument) length is not equal to "
		"p (first argument)."
	);

	// vector of sets that will hold the results
	sparse_set     var_sparsity;
	var_sparsity.resize(total_num_var, p);

	// The sparsity pattern corresponding to the dependent variables
	for(i = 0; i < p; i++)
	{	itr = s[i].begin();
		while(itr != s[i].end())
		{	j = *itr++; 
			CPPAD_ASSERT_KNOWN(
				j < m,
				"RevSparseJac: an element of the set s[i] "
				"has value greater than or equal m."
			);
			CPPAD_ASSERT_UNKNOWN( dep_taddr[j] < total_num_var );
			var_sparsity.add_element( dep_taddr[j], i );
		}
	}

	// evaluate the sparsity patterns
	RevJacSweep(
		n,
		total_num_var,
		&play,
		var_sparsity
	);

	// return values corresponding to dependent variables
	CPPAD_ASSERT_UNKNOWN( r.size() == p );
	for(j = 0; j < n; j++)
	{	CPPAD_ASSERT_UNKNOWN( ind_taddr[j] == (j+1) );

		// ind_taddr[j] is operator taddr for j-th independent variable
		CPPAD_ASSERT_UNKNOWN( play.GetOp( ind_taddr[j] ) == InvOp );

		// extract result from rev_hes_sparsity
		// and add corresponding elements to sets in r
		CPPAD_ASSERT_UNKNOWN( var_sparsity.end() == p );
		var_sparsity.begin(j+1);
		i = var_sparsity.next_element();
		while( i < p )
		{	r[i].insert(j);
			i = var_sparsity.next_element();
		}
	}
}

/*!
User API for Jacobian sparsity patterns using reverse mode.

The C++ source code corresponding to this operation is
\verbatim
	s = f.RevSparseJac(q, r)
\endverbatim

\tparam Base
is the base type for this recording.

\tparam VectorSet
is a simple vector with elements of type \c bool.
or \c std::set<size_t>.

\param p
is the number of rows in the matrix \f$ S \f$.

\param s
is a sparsity pattern for the matrix \f$ S \f$.

\return
If \c VectorSet::value_type is \c bool,
the return value \c r is a vector with size \c p*n
where \c n is the number of independent variables
corresponding to the operation sequence stored in \c f. 
If \c VectorSet::value_type is \c std::set<size_t>,
the return value \c r is a vector of sets with size \c p
and with all its elements between zero and \c n - 1.
The value of \a r is the sparsity pattern for the matrix
\f[
	J(x) = S * F^{(1)} (x)
\f]
where \f$ F \f$ is the function corresponding to the operation sequence
and \a x is any argument value.
*/
template <class Base>
template <class VectorSet>
VectorSet ADFun<Base>::RevSparseJac(
	size_t              p      , 
	const VectorSet&    s      )
{
	VectorSet r;
	typedef typename VectorSet::value_type Set_type;

	RevSparseJacCase(
		Set_type()    ,
		p             ,
		s             ,
		r
	);
	return r;
}

/*!
Private helper function for \c RevSparseJac(p, s).

All of the description in the public member function \c RevSparseJac(p, s)
applies.

\param set_type
is a \c bool value.
This arugment is used to dispatch to the proper source code
depending on the value of \c VectorSet::value_type.

\param p
See \c RevSparseJac(p, s).

\param s
See \c RevSparseJac(p, s).

\param r
is the return value for the corresponding call to RevSparseJac(p, s);
*/

template <class Base>
template <class VectorSet>
void ADFun<Base>::RevSparseJacCase(
	bool                set_type          ,
	size_t              p                 ,
	const VectorSet&    s                 ,
	VectorSet&          r                 )
{	size_t n = Domain();

	// dimension of the result vector
	r.resize( p * n );

	// store results in r
	RevSparseJacBool(
		p              ,
		s              ,
		r              ,
		total_num_var_ ,
		dep_taddr_     ,
		ind_taddr_     ,
		play_
	);
}

/*!
Private helper function for \c RevSparseJac(p, s).

All of the description in the public member function \c RevSparseJac(p, s)
applies.

\param set_type
is a \c std::set<size_t> object.
This arugment is used to dispatch to the proper source code
depending on the value of \c VectorSet::value_type.

\param p
See \c RevSparseJac(p, s).

\param s
See \c RevSparseJac(p, s).

\param r
is the return value for the corresponding call to RevSparseJac(p, s);
*/

template <class Base>
template <class VectorSet>
void ADFun<Base>::RevSparseJacCase(
	const std::set<size_t>&      set_type          ,
	size_t                       p                 ,
	const VectorSet&             s                 ,
	VectorSet&                   r                 )
{	// dimension of the result vector
	r.resize( p );

	// store results in r
	RevSparseJacSet(
		p              ,
		s              ,
		r              ,
		total_num_var_ ,
		dep_taddr_     ,
		ind_taddr_     ,
		play_
	);
}

CPPAD_END_NAMESPACE
# endif

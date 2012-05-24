/* $Id: rev_sparse_hes.hpp 2023 2011-07-25 17:44:12Z bradbell $ */
# ifndef CPPAD_REV_SPARSE_HES_INCLUDED
# define CPPAD_REV_SPARSE_HES_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin RevSparseHes$$
$spell
	std
	VecAD
	Jacobian
	Jac
	Hessian
	Hes
	const
	Bool
	Dep
	proportional
	var
$$

$section Hessian Sparsity Pattern: Reverse Mode$$ 

$index RevSparseHes$$
$index reverse, sparse Hessian$$
$index sparse, reverse Hessian$$
$index pattern, reverse Hessian$$

$head Syntax$$
$icode%h% = %f%.RevSparseHes(%q%, %s%)%$$


$head Purpose$$
We use $latex F : \R^n \rightarrow \R^m$$ to denote the
$xref/glossary/AD Function/AD function/$$ corresponding to $icode f$$.
We define the matrix $latex H(x) \in \R^{q \times n}$$
as the partial with respect to $latex x$$,
of the partial with respect to $latex u$$ (at $latex u = 0$$),
of $latex S * F[ x + R * u ]$$ where
$latex R \in \R^{n \times q}$$ and $latex S \in \R^{1 \times m}$$; i.e.,
$latex \[
	H(x)  =  R^\T (S * F)^{(2)} ( x )
\] $$
Given a
$xref/glossary/Sparsity Pattern/sparsity pattern/$$ 
for the matrix $latex R$$ and the vector $latex S$$,
$code RevSparseHes$$ returns a sparsity pattern for the $latex H(x)$$.

$head f$$
The object $icode f$$ has prototype
$codei%
	const ADFun<%Base%> %f%
%$$

$head x$$
the sparsity pattern is valid for all values of the independent
variables in $latex x \in \R^n$$
(even if it has $cref/CondExp/$$ or $cref/VecAD/$$ operations).

$head q$$
The argument $icode q$$ has prototype
$codei%
	size_t %q%
%$$
It specifies the number of columns in $latex R \in \R^{n \times q}$$
and the number of rows in 
$latex H(x) \in \R^{q \times n}$$.
It must be the same value as in the previous $xref/ForSparseJac/$$ call 
$codei%
	%f%.ForSparseJac(%q%, %r%)
%$$

$head r$$
The argument $icode r$$ in the previous call
$codei%
	%f%.ForSparseJac(%q%, %r%)
%$$
is a $xref/glossary/Sparsity Pattern/sparsity pattern/$$
for the matrix $latex R$$ above.
The type of the elements of
$cref/VectorSet/RevSparseHes/VectorSet/$$ must be the 
same as the type of the elements of $icode r$$.

$head s$$
The argument $icode s$$ has prototype
$codei%
	const %VectorSet%& %s%
%$$
(see $xref/RevSparseHes/VectorSet/VectorSet/$$ below)
If it has elements of type $code bool$$,
its size is $latex m$$.
If it has elements of type $code std::set<size_t>$$,
its size is one and all the elements of $icode%s%[0]%$$
are between zero and $latex m - 1$$.
It specifies a 
$xref/glossary/Sparsity Pattern/sparsity pattern/$$ 
for the vector $icode S$$.

$head h$$
The result $icode h$$ has prototype
$codei%
	%VectorSet%& %h%
%$$
(see $xref/RevSparseHes/VectorSet/VectorSet/$$ below).
If it has elements of type $code bool$$,
its size is $latex q * n$$.
If it has elements of type $code std::set<size_t>$$,
its size is $latex q$$.
It specifies a 
$xref/glossary/Sparsity Pattern/sparsity pattern/$$ 
for the matrix $latex H(x)$$.

$head VectorSet$$
The type $icode VectorSet$$ must be a $cref/SimpleVector/$$ class with
$xref/SimpleVector/Elements of Specified Type/elements of type/$$
$code bool$$ or $code std::set<size_t>$$;
see $cref/sparsity pattern/glossary/Sparsity Pattern/$$ for a discussion
of the difference.
The type of the elements of
$cref/VectorSet/RevSparseHes/VectorSet/$$ must be the 
same as the type of the elements of $icode r$$.

$head Entire Sparsity Pattern$$
Suppose that $latex q = n$$ and
$latex R \in \R^{n \times n}$$ is the $latex n \times n$$ identity matrix.
Further suppose that the $latex S$$ is the $th k$$
$xref/glossary/Elementary Vector/elementary vector/$$; i.e.
$latex \[
S_j = \left\{ \begin{array}{ll}
	1  & {\rm if} \; j = k 
	\\
	0  & {\rm otherwise}
\end{array} \right. 
\] $$
In this case,
the corresponding value $icode h$$ is a 
sparsity pattern for the Hessian matrix
$latex F_k^{(2)} (x) \in \R^{n \times n}$$.

$head Example$$
$children%
	example/rev_sparse_hes.cpp
%$$
The file
$xref/RevSparseHes.cpp/$$
contains an example and test of this operation.
It returns true if it succeeds and false otherwise.

$end
-----------------------------------------------------------------------------
*/
# include <algorithm>
# include <cppad/local/pod_vector.hpp>
# include <cppad/local/std_set.hpp>

CPPAD_BEGIN_NAMESPACE
/*!
\file rev_sparse_hes.hpp
Reverse mode Hessian sparsity patterns.
*/

/*!
Calculate Hessian sparsity patterns using reverse mode.

The C++ source code corresponding to this operation is
\verbatim
	h = f.RevSparseHes(q, s)
\endverbatim

\tparam Base
is the base type for this recording.

\tparam VectorSet
is a simple vector with elements of type \c bool.

\tparam Sparsity
is either \c sparse_pack or \c sparse_set. 

\param q
is the value of \a q in the 
by the previous call of the form 
\verbatim
	f.ForSparseJac(q, r)
\endverbatim
The value \c r in this call is a sparsity pattern for the matrix \f$ R \f$.

\param s
is a vector with size \c m that specifies the sparsity pattern
for the vector \f$ S \f$,
where \c m is the number of dependent variables
corresponding to the operation sequence stored in \a play. 

\param h
the input value of \a h must be a vector with size \c q*n.
The input value of its elements does not matter.
On output, \a h is the sparsity pattern for the matrix
\f[
	H(x) = R^T (S * F)^{(2)} (x)
\f]
where \f$ F \f$ is the function corresponding to the operation sequence
and \a x is any argument value.

\param total_num_var
is the total number of variables in this recording.

\param dep_taddr
maps dependendent variable index
to the corresponding variable in the tape.

\param ind_taddr
maps independent variable index
to the corresponding variable in the tape.

\param play
is the recording that defines the function we are computing the sparsity 
pattern for.

\param for_jac_sparsity
is a vector of sets containing the 
the forward Jacobian sparsity pattern corresponding to 
$latex R$$ for all of the variables on the tape. 
*/

template <class Base, class VectorSet, class Sparsity>
void RevSparseHesBool(
	size_t                    q                 ,
	const VectorSet&          s                 ,
	VectorSet&                h                 ,
	size_t                    total_num_var     ,
	CppAD::vector<size_t>&    dep_taddr         ,
	CppAD::vector<size_t>&    ind_taddr         ,
	CppAD::player<Base>&      play              ,
	Sparsity&                 for_jac_sparsity  )
{
	// temporary indices
	size_t i, j;

	// check Vector is Simple VectorSet class with bool elements
	CheckSimpleVector<bool, VectorSet>();

	// range and domain dimensions for F
	size_t m = dep_taddr.size();
	size_t n = ind_taddr.size();

	CPPAD_ASSERT_KNOWN(
		q == for_jac_sparsity.end(),
		"RevSparseHes: q (first argument) is not equal to its value\n"
		"in the previous call to ForSparseJac with this ADFun object."
	);
	CPPAD_ASSERT_KNOWN(
		s.size() == m,
		"RevSparseHes: s (second argument) length is not equal to\n"
		"range dimension for ADFun object."
	);

	// Array that will hold reverse Jacobian dependency flag.
	// Initialize as true for the dependent variables.
	pod_vector<bool> RevJac;
	RevJac.extend(total_num_var);	
	for(i = 0; i < total_num_var; i++)
		RevJac[i] = false;
	for(i = 0; i < m; i++)
	{	CPPAD_ASSERT_UNKNOWN( dep_taddr[i] < total_num_var );
		RevJac[ dep_taddr[i] ] = s[i];
	}


	// vector of sets that will hold reverse Hessain values
	Sparsity       rev_hes_sparsity;
	rev_hes_sparsity.resize(total_num_var, q);

	// compute the Hessian sparsity patterns
	RevHesSweep(
		n,
		total_num_var,
		&play,
		for_jac_sparsity, 
		RevJac.data(),
		rev_hes_sparsity
	);

	// return values corresponding to independent variables
	CPPAD_ASSERT_UNKNOWN( h.size() == n * q );

	// j is index corresponding to reverse mode partial
	for(j = 0; j < n; j++)
	{	CPPAD_ASSERT_UNKNOWN( ind_taddr[j] < total_num_var );

		// ind_taddr[j] is operator taddr for j-th independent variable
		CPPAD_ASSERT_UNKNOWN( ind_taddr[j] == j + 1 );
		CPPAD_ASSERT_UNKNOWN( play.GetOp( ind_taddr[j] ) == InvOp );

		// extract the result from rev_hes_sparsity
		for(i = 0; i < q; i++) 
			h[ i * n + j ] = false;
		CPPAD_ASSERT_UNKNOWN( rev_hes_sparsity.end() == q );
		rev_hes_sparsity.begin(j + 1);
		i = rev_hes_sparsity.next_element();
		while( i < q )
		{	h[ i * n + j ] = true;
			i = rev_hes_sparsity.next_element();
		}
	}

	return;
}

/*!
Calculate Hessian sparsity patterns using reverse mode.

The C++ source code corresponding to this operation is
\verbatim
	h = f.RevSparseHes(q, s)
\endverbatim

\tparam Base
is the base type for this recording.

\tparam VectorSet
is a simple vector with elements of type \c std::set<size_t>.

\tparam Sparsity
is either \c sparse_pack or \c sparse_set. 

\param q
is the value of \a q in the 
by the previous call of the form 
\verbatim
	f.ForSparseJac(q, r)
\endverbatim
The value \c r in this call is a sparsity pattern for the matrix \f$ R \f$.

\param s
is a vector with size \c m that specifies the sparsity pattern
for the vector \f$ S \f$,
where \c m is the number of dependent variables
corresponding to the operation sequence stored in \a play. 

\param h
the input value of \a h must be a vector with size \a q.
On input, each element of \a h must be an empty set.
The input value of its elements does not matter.
On output, \a h is the sparsity pattern for the matrix
\f[
	H(x) = R^T (S * F)^{(2)} (x)
\f]
where \f$ F \f$ is the function corresponding to the operation sequence
and \a x is any argument value.

\param total_num_var
is the total number of variables in this recording.

\param dep_taddr
maps dependendent variable index
to the corresponding variable in the tape.

\param ind_taddr
maps independent variable index
to the corresponding variable in the tape.

\param play
is the recording that defines the function we are computing the sparsity 
pattern for.

\param for_jac_sparsity
is a vector of sets containing the 
the forward Jacobian sparsity pattern corresponding to 
$latex R$$ for all of the variables on the tape. 
*/

template <class Base, class VectorSet, class Sparsity>
void RevSparseHesSet(
	size_t                    q                 ,
	const VectorSet&          s                 ,
	VectorSet&                h                 ,
	size_t                    total_num_var     ,
	CppAD::vector<size_t>&    dep_taddr         ,
	CppAD::vector<size_t>&    ind_taddr         ,
	CppAD::player<Base>&      play              ,
	Sparsity&                 for_jac_sparsity  )
{
	// temporary indices
	size_t i, j;
	std::set<size_t>::const_iterator itr;

	// check VectorSet is Simple Vector class with sets for elements
	CheckSimpleVector<std::set<size_t>, VectorSet>(
		one_element_std_set<size_t>(), two_element_std_set<size_t>()
	);

	// range and domain dimensions for F
# ifndef NDEBUG
	size_t m = dep_taddr.size();
# endif
	size_t n = ind_taddr.size();

	CPPAD_ASSERT_KNOWN(
		q == for_jac_sparsity.end(),
		"RevSparseHes: q (first argument) is not equal to its value\n"
		"in the previous call to ForSparseJac with this ADFun object."
	);
	CPPAD_ASSERT_KNOWN(
		s.size() == 1,
		"RevSparseHes: s (second argument) length is not equal to one."
	);

	// Array that will hold reverse Jacobian dependency flag.
	// Initialize as true for the dependent variables.
	pod_vector<bool> RevJac;
	RevJac.extend(total_num_var);	
	for(i = 0; i < total_num_var; i++)
		RevJac[i] = false;
	itr = s[0].begin();
	while( itr != s[0].end() )
	{	i = *itr++;
		CPPAD_ASSERT_KNOWN(
			i < m,
			"RevSparseHes: an element of the set s[0] has value "
			"greater than or equal m"
		);
		CPPAD_ASSERT_UNKNOWN( dep_taddr[i] < total_num_var );
		RevJac[ dep_taddr[i] ] = true;
	}


	// vector of sets that will hold reverse Hessain values
	Sparsity       rev_hes_sparsity;
	rev_hes_sparsity.resize(total_num_var, q);

	// compute the Hessian sparsity patterns
	RevHesSweep(
		n,
		total_num_var,
		&play,
		for_jac_sparsity, 
		RevJac.data(),
		rev_hes_sparsity
	);

	// return values corresponding to independent variables
	// j is index corresponding to reverse mode partial
	CPPAD_ASSERT_UNKNOWN( h.size() == q );
	for(j = 0; j < n; j++)
	{	CPPAD_ASSERT_UNKNOWN( ind_taddr[j] < total_num_var );
		CPPAD_ASSERT_UNKNOWN( ind_taddr[j] == j + 1 );
		CPPAD_ASSERT_UNKNOWN( play.GetOp( ind_taddr[j] ) == InvOp );

		// extract the result from rev_hes_sparsity
		// and add corresponding elements to result sets in h
		CPPAD_ASSERT_UNKNOWN( rev_hes_sparsity.end() == q );
		rev_hes_sparsity.begin(j+1);
		i = rev_hes_sparsity.next_element();
		while( i < q )
		{	h[i].insert(j);
			i = rev_hes_sparsity.next_element();
		}
	}

	return;
}



/*!
User API for Hessian sparsity patterns using reverse mode.

The C++ source code corresponding to this operation is
\verbatim
	h = f.RevSparseHes(q, r)
\endverbatim

\tparam Base
is the base type for this recording.

\tparam VectorSet
is a simple vector with elements of type \c bool
or \c std::set<size_t>.

\param q
is the value of \a q in the 
by the previous call of the form 
\verbatim
	f.ForSparseJac(q, r, packed)
\endverbatim
The value \c r in this call is a sparsity pattern for the matrix \f$ R \f$.
The type of the element of \c r for the previous call to \c ForSparseJac
must be the same as the type of the elements of \c s.

\param s
is a vector with size \c m that specifies the sparsity pattern
for the vector \f$ S \f$,
where \c m is the number of dependent variables
corresponding to the operation sequence stored in \a play. 

\return
is a vector with size \c q*n.
containing a sparsity pattern for the matrix
\f[
	H(x) = R^T ( S * F)^{(2)} (x)
\f]
where \f$ F \f$ is the function corresponding to the operation sequence
and \a x is any argument value.
*/

template <class Base>
template <class VectorSet>
VectorSet ADFun<Base>::RevSparseHes(size_t q,  const VectorSet& s)
{	VectorSet h;
	typedef typename VectorSet::value_type Set_type;

	// Should check to make sure q is same as in previous call to
	// forward sparse Jacobian.
	RevSparseHesCase(
		Set_type()    ,
		q             ,
		s             ,
		h
	);

	return h;
}

/*!
Private helper function for RevSparseHes(q, s).

All of the description in the public member function RevSparseHes(q, s)
applies.

\param set_type
is a \c bool value. This argument is used to dispatch to the proper source
code depending on the vlaue of \c VectorSet::value_type.

\param q
See \c RevSparseHes(q, s).

\param s
See \c RevSparseHes(q, s).

\param h
is the return value for the corresponging call to \c RevSparseJac(q, s).
*/
template <class Base>
template <class VectorSet>
void ADFun<Base>::RevSparseHesCase(
	bool              set_type         ,
	size_t            q                ,  
	const VectorSet&  s                ,
	VectorSet&        h                )
{	size_t n = Domain(); 	
	h.resize(q * n );

	CPPAD_ASSERT_KNOWN( 
		for_jac_sparse_pack_.n_set() > 0,
		"RevSparseHes: previous stored call to ForSparseJac did not "
		"use bool for the elements of r."
	);
	CPPAD_ASSERT_UNKNOWN( for_jac_sparse_set_.n_set() == 0 );
	CPPAD_ASSERT_UNKNOWN( for_jac_sparse_pack_.n_set() == total_num_var_ );
	
	// use sparse_pack for the calculation
	CppAD::RevSparseHesBool( 
		q                        ,
		s                        ,
		h                        ,
		total_num_var_           ,
		dep_taddr_               ,
		ind_taddr_               ,
		play_                    ,
		for_jac_sparse_pack_ 
	);
}

/*!
Private helper function for RevSparseHes(q, s).

All of the description in the public member function RevSparseHes(q, s)
applies.

\param set_type
is a \c std::set<size_t> value. 
This argument is used to dispatch to the proper source
code depending on the vlaue of \c VectorSet::value_type.

\param q
See \c RevSparseHes(q, s).

\param s
See \c RevSparseHes(q, s).

\param h
is the return value for the corresponging call to \c RevSparseJac(q, s).
*/
template <class Base>
template <class VectorSet>
void ADFun<Base>::RevSparseHesCase(
	const std::set<size_t>&   set_type         ,
	size_t                    q                ,  
	const VectorSet&          s                ,
	VectorSet&                h                )
{	h.resize(q);

	CPPAD_ASSERT_KNOWN( 
		for_jac_sparse_set_.n_set() > 0,
		"RevSparseHes: previous stored call to ForSparseJac did not "
		"use std::set<size_t> for the elements of r."
	);
	CPPAD_ASSERT_UNKNOWN( for_jac_sparse_pack_.n_set() == 0 );
	CPPAD_ASSERT_UNKNOWN( for_jac_sparse_set_.n_set() == total_num_var_ );
	
	// use sparse_pack for the calculation
	CppAD::RevSparseHesSet( 
		q                        ,
		s                        ,
		h                        ,
		total_num_var_           ,
		dep_taddr_               ,
		ind_taddr_               ,
		play_                    ,
		for_jac_sparse_set_ 
	);
}

CPPAD_END_NAMESPACE
# endif

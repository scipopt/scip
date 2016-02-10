/* $Id: forward.hpp 3098 2014-02-18 03:19:39Z bradbell $ */
# ifndef CPPAD_FORWARD_INCLUDED
# define CPPAD_FORWARD_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-14 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Eclipse Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin Forward$$

$section Forward Mode$$

$childtable%
	omh/forward.omh%
	cppad/local/cap_taylor.hpp%
	cppad/local/num_skip.hpp%
	example/forward.cpp%
	example/forward_mul.cpp
%$$

$end
-----------------------------------------------------------------------------
*/

// documened after Forward but included here so easy to see
# include <cppad/local/cap_taylor.hpp>
# include <cppad/local/num_skip.hpp>

namespace CppAD { // BEGIN_CPPAD_NAMESPACE
/*!
\defgroup forward_hpp forward.hpp
\{
\file forward.hpp
User interface to forward mode computations
*/

/*!
Compute arbitrary order forward mode Taylor coefficieints.

\tparam Base
The type used during the forward mode computations; i.e., the corresponding
recording of operations used the type \c AD<Base>.

\tparam Vector
is a Simple Vector class with eleements of type \c Base.

\param p
is the hightest order for this forward mode computation; i.e., 
after this calculation there will be <code>p+1</code>
Taylor coefficients per variables.

\param x_p
contains Taylor coefficients for the independent variables.
The size of \a x_p must either be \c n or <code>n*(p+1)</code>,
We define <code>q = p + 1 - x_p.size() / n</code>. 
The Taylor coefficients of order k, for
k = q, ... , p are calculated.

\param s
Is the stream where output corresponding to \c PriOp operations will written.
*/

template <typename Base>
template <typename Vector>
Vector ADFun<Base>::Forward(
	size_t p                    , 
	const Vector& x_p           , 
	std::ostream& s             )
{	// temporary indices
	size_t i, j, k;

	// number of independent variables
	size_t n = ind_taddr_.size();

	// number of dependent variables
	size_t m = dep_taddr_.size();

	// check Vector is Simple Vector class with Base type elements
	CheckSimpleVector<Base, Vector>();

	CPPAD_ASSERT_KNOWN(
		size_t(x_p.size()) == n || size_t(x_p.size()) == n*(p+1),
		"Forward: x_p.size() is not equal n or n*(p+1)"
	);
	size_t n_order = size_t(x_p.size()) / n;
	CPPAD_ASSERT_KNOWN(
		p <= taylor_per_var_ || n_order == p + 1,
		"The number of Taylor coefficient currently stored in this ADFun\n"
		"is less than p and x_p.size() != n*(p+1)."
	);  

	// check if the taylor_ matrix needs more columns
	if( taylor_col_dim_ <= p )
		capacity_taylor(p + 1);
	CPPAD_ASSERT_UNKNOWN( taylor_col_dim_ > p );

	// set the p-th order taylor_ coefficients for independent variables
	for(j = 0; j < n; j++)
	{	CPPAD_ASSERT_UNKNOWN( ind_taddr_[j] < total_num_var_ );

		// ind_taddr_[j] is operator taddr for j-th independent variable
		CPPAD_ASSERT_UNKNOWN( play_.GetOp( ind_taddr_[j] ) == InvOp );

		// It is also variable taddr for j-th independent variable
		if( n_order ==  1 )
			taylor_[ind_taddr_[j] * taylor_col_dim_ + p] = x_p[j];
		else for(k = 0; k < n_order; k++)
			taylor_[ind_taddr_[j] * taylor_col_dim_ + k] = 
				x_p[j * n_order + k];
	}

	// evaluate the derivatives
	size_t q = (p + 1) - n_order;
	if( p == 0 )
	{
# if CPPAD_USE_FORWARD0SWEEP
		compare_change_ = forward0sweep(s, true,
			n, total_num_var_, &play_, taylor_col_dim_, taylor_.data(),
			cskip_op_
		);
# else
		compare_change_ = forward_sweep(s, true, q,
			p, n, total_num_var_, &play_, taylor_col_dim_, taylor_.data(),
			cskip_op_
		);
# endif
	}
	else if( q == 0 )
	{	compare_change_ = forward_sweep(s, true, q,
			p, n, total_num_var_, &play_, taylor_col_dim_, taylor_.data(),
			cskip_op_
		);
	}
	else
	{	forward_sweep(s, true, q,
			p, n, total_num_var_, &play_, taylor_col_dim_, taylor_.data(),
			cskip_op_
		);
	}

	// return Taylor coefficients for dependent variables
	Vector y_p;
	if( n_order == 1 )
	{	y_p.resize(m);
		for(i = 0; i < m; i++)
		{	CPPAD_ASSERT_UNKNOWN( dep_taddr_[i] < total_num_var_ );
			y_p[i] = taylor_[dep_taddr_[i] * taylor_col_dim_ + p];
		}
	}
	else
	{	y_p.resize(m * n_order );
		for(i = 0; i < m; i++)	
		{	for(k = 0; k < n_order; k++)
				y_p[ i * n_order + k] = 
					taylor_[ dep_taddr_[i] * taylor_col_dim_ + k ]; 
		}
	}
# ifndef NDEBUG
	if( check_for_nan_ )
	{	bool ok = true;
		if( p == 0 && n_order == 1 )
			ok = ! hasnan(y_p);
		else if( n_order != 1 )
		{	for(i = 0; i < m; i++)
			// on MS Visual Studio 2012, CppAD required in front of isnan ?
			ok &= ! CppAD::isnan( y_p[ i * n_order + 0 ] );
		} 
		CPPAD_ASSERT_KNOWN(ok,
			"y_p = f.Forward(p, x): has a zero order Taylor coefficient "
			"with the value nan."
		);  
		if( p != 0 && n_order == 1 )
			ok = ! hasnan(y_p);
		else if( n_order != 1 )
		{	for(i = 0; i < m; i++)
			{	for(k = 1; k < n_order; k++)
					// Studio 2012, CppAD required in front of isnan ?
					ok &= ! CppAD::isnan( y_p[ i * n_order + k ] );
			}
		}
		CPPAD_ASSERT_KNOWN(ok,
		"y_p = f.Forward(p, x): has a non-zero order Taylor coefficient\n"
		"with the value nan (but zero order coefficients are not nan)."
		);
	}
# endif


	// now we have p + 1  taylor_ coefficients per variable
	taylor_per_var_ = p + 1;

	return y_p;
}

/*! \} */
} // END_CPPAD_NAMESPACE
# endif

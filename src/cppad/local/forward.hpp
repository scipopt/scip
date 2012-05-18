/* $Id: forward.hpp 2114 2011-10-14 18:59:16Z bradbell $ */
# ifndef CPPAD_FORWARD_INCLUDED
# define CPPAD_FORWARD_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin Forward$$

$section Forward Mode$$

$childtable%
	omh/forward.omh%
	cppad/local/cap_taylor.hpp%
	example/forward.cpp
%$$

$end
-----------------------------------------------------------------------------
*/

// documened after Forward but included here so easy to see
# include <cppad/local/cap_taylor.hpp>

CPPAD_BEGIN_NAMESPACE
/*!
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
is the order for this forward mode computation; i.e., the number
of Taylor coefficient is <code>p+1</code>.

\param x_p
Is the \c p th order Taylor coefficient vector for the independent variables.

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
	size_t i, j;

	// number of independent variables
	size_t n = ind_taddr_.size();

	// number of dependent variables
	size_t m = dep_taddr_.size();

	// check Vector is Simple Vector class with Base type elements
	CheckSimpleVector<Base, Vector>();

	CPPAD_ASSERT_KNOWN(
		x_p.size() == n,
		"Second argument to Forward does not have length equal to\n"
		"the dimension of the domain for the corresponding ADFun."
	);
	CPPAD_ASSERT_KNOWN(
		p <= taylor_per_var_,
		"The number of taylor_ coefficient currently stored\n"
		"in this ADFun object is less than p."
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
		taylor_[ind_taddr_[j] * taylor_col_dim_ + p] = x_p[j];
	}

	// evaluate the derivatives
	if( p == 0 )
	{
# if CPPAD_USE_FORWARD0SWEEP
		compare_change_ = forward0sweep(s, true,
			n, total_num_var_, &play_, taylor_col_dim_, taylor_.data()
		);
# else
		compare_change_ = forward_sweep(s, true,
			p, n, total_num_var_, &play_, taylor_col_dim_, taylor_.data()
		);
# endif
	}
	else forward_sweep(s, false,
		p, n, total_num_var_, &play_, taylor_col_dim_, taylor_.data()
	);

	// return the p-th order taylor_ coefficients for dependent variables
	Vector y_p(m);
	for(i = 0; i < m; i++)
	{	CPPAD_ASSERT_UNKNOWN( dep_taddr_[i] < total_num_var_ );
		y_p[i] = taylor_[dep_taddr_[i] * taylor_col_dim_ + p];
	}
# ifndef NDEBUG
	if( hasnan(y_p) )
	{	if( p == 0 )
		{	CPPAD_ASSERT_KNOWN(false,
				"y = f.Forward(0, x): has a nan in y."
			);  
		}
		else
		{	CPPAD_ASSERT_KNOWN(false,
				"y_p = f.Forward(p, x_p): has a nan in y_p for p > 0, "
				"but not for p = 0."
			);
		}
	}
# endif


	// now we have p + 1  taylor_ coefficients per variable
	taylor_per_var_ = p + 1;

	return y_p;
}

CPPAD_END_NAMESPACE
# endif

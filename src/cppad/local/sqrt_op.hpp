/* $Id: sqrt_op.hpp 1641 2010-02-01 16:39:45Z bradbell $ */
# ifndef CPPAD_SQRT_OP_INCLUDED
# define CPPAD_SQRT_OP_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-10 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */


CPPAD_BEGIN_NAMESPACE
/*!
\file sqrt_op.hpp
Forward and reverse mode calculations for z = sqrt(x).
*/


/*!
Compute forward mode Taylor coefficient for result of op = SqrtOp.

The C++ source code corresponding to this operation is
\verbatim
	z = sqrt(x)
\endverbatim

\copydetails forward_unary1_op
*/
template <class Base>
inline void forward_sqrt_op(
	size_t j           ,
	size_t i_z         ,
	size_t i_x         ,
	size_t nc_taylor   , 
	Base*  taylor      )
{	
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(SqrtOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( NumRes(SqrtOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( i_x < i_z );
	CPPAD_ASSERT_UNKNOWN( j < nc_taylor );

	// Taylor coefficients corresponding to argument and result
	Base* x = taylor + i_x * nc_taylor;
	Base* z = taylor + i_z * nc_taylor;

	size_t k;
	if( j == 0 )
		z[j] = sqrt( x[0] );
	else
	{
		z[j] = Base(0);
		for(k = 1; k < j; k++)
			z[j] -= Base(k) * z[k] * z[j-k];
		z[j] /= Base(j);
		z[j] += x[j] / Base(2);
		z[j] /= z[0];
	}
}

/*!
Compute zero order forward mode Taylor coefficient for result of op = SqrtOp.

The C++ source code corresponding to this operation is
\verbatim
	z = sqrt(x)
\endverbatim

\copydetails forward_unary1_op_0
*/
template <class Base>
inline void forward_sqrt_op_0(
	size_t i_z         ,
	size_t i_x         ,
	size_t nc_taylor   , 
	Base*  taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(SqrtOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( NumRes(SqrtOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( i_x < i_z );
	CPPAD_ASSERT_UNKNOWN( 0 < nc_taylor );

	// Taylor coefficients corresponding to argument and result
	Base* x = taylor + i_x * nc_taylor;
	Base* z = taylor + i_z * nc_taylor;

	z[0] = sqrt( x[0] );
}
/*!
Compute reverse mode partial derivatives for result of op = SqrtOp.

The C++ source code corresponding to this operation is
\verbatim
	z = sqrt(x)
\endverbatim

\copydetails reverse_unary1_op
*/

template <class Base>
inline void reverse_sqrt_op(
	size_t      d            ,
	size_t      i_z          ,
	size_t      i_x          ,
	size_t      nc_taylor    , 
	const Base* taylor       ,
	size_t      nc_partial   ,
	Base*       partial      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(SqrtOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( NumRes(SqrtOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( i_x < i_z );
	CPPAD_ASSERT_UNKNOWN( d < nc_taylor );
	CPPAD_ASSERT_UNKNOWN( d < nc_partial );

	// Taylor coefficients and partials corresponding to argument
	Base* px       = partial + i_x * nc_partial;

	// Taylor coefficients and partials corresponding to result
	const Base* z  = taylor  + i_z * nc_taylor;
	Base* pz       = partial + i_z * nc_partial;

	// number of indices to access
	size_t j = d;
	size_t k;
	while(j)
	{	// scale partial w.r.t. z[j]
		pz[j]   /= z[0];

		pz[0]   -= pz[j] * z[j];
		px[j]   += pz[j] / Base(2);
		for(k = 1; k < j; k++)
			pz[k]   -= pz[j] * z[j-k];
		--j;
	}
	px[0] += pz[0] / (Base(2) * z[0]);
}

CPPAD_END_NAMESPACE
# endif

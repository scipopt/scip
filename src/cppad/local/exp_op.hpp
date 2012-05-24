/* $Id: exp_op.hpp 1641 2010-02-01 16:39:45Z bradbell $ */
# ifndef CPPAD_EXP_OP_INCLUDED
# define CPPAD_EXP_OP_INCLUDED

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
\file exp_op.hpp
Forward and reverse mode calculations for z = exp(x).
*/


/*!
Forward mode Taylor coefficient for result of op = ExpOp.

The C++ source code corresponding to this operation is
\verbatim
	z = exp(x)
\endverbatim

\copydetails forward_unary1_op
*/
template <class Base>
inline void forward_exp_op(
	size_t j           ,
	size_t i_z         ,
	size_t i_x         ,
	size_t nc_taylor   , 
	Base*  taylor      )
{	
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(ExpOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( NumRes(ExpOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( i_x < i_z );
	CPPAD_ASSERT_UNKNOWN( j < nc_taylor );

	// Taylor coefficients corresponding to argument and result
	Base* x = taylor + i_x * nc_taylor;
	Base* z = taylor + i_z * nc_taylor;

	size_t k;
	if( j == 0 )
		z[0] = exp( x[0] );
	else
	{
		z[j] = x[1] * z[j-1];
		for(k = 2; k <= j; k++)
			z[j] += Base(k) * x[k] * z[j-k];
		z[j] /= Base(j);
	}
}

/*!
Zero order forward mode Taylor coefficient for result of op = ExpOp.

The C++ source code corresponding to this operation is
\verbatim
	z = exp(x)
\endverbatim

\copydetails forward_unary1_op_0
*/
template <class Base>
inline void forward_exp_op_0(
	size_t i_z         ,
	size_t i_x         ,
	size_t nc_taylor   , 
	Base*  taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(ExpOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( NumRes(ExpOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( i_x < i_z );
	CPPAD_ASSERT_UNKNOWN( 0 < nc_taylor );

	// Taylor coefficients corresponding to argument and result
	Base* x = taylor + i_x * nc_taylor;
	Base* z = taylor + i_z * nc_taylor;

	z[0] = exp( x[0] );
}
/*!
Reverse mode partial derivatives for result of op = ExpOp.

The C++ source code corresponding to this operation is
\verbatim
	z = exp(x)
\endverbatim

\copydetails reverse_unary1_op
*/

template <class Base>
inline void reverse_exp_op(
	size_t      d            ,
	size_t      i_z          ,
	size_t      i_x          ,
	size_t      nc_taylor    , 
	const Base* taylor       ,
	size_t      nc_partial   ,
	Base*       partial      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(ExpOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( NumRes(ExpOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( i_x < i_z );
	CPPAD_ASSERT_UNKNOWN( d < nc_taylor );
	CPPAD_ASSERT_UNKNOWN( d < nc_partial );

	// Taylor coefficients and partials corresponding to argument
	const Base* x  = taylor  + i_x * nc_taylor;
	Base* px       = partial + i_x * nc_partial;

	// Taylor coefficients and partials corresponding to result
	const Base* z  = taylor  + i_z * nc_taylor;
	Base* pz       = partial + i_z * nc_partial;

	// lopp through orders in reverse
	size_t j, k;
	j = d;
	while(j)
	{	// scale partial w.r.t z[j]
		pz[j] /= Base(j);

		for(k = 1; k <= j; k++)
		{	px[k]   += pz[j] * Base(k) * z[j-k]; 	
			pz[j-k] += pz[j] * Base(k) * x[k];
		}
		--j;
	}
	px[0] += pz[0] * z[0];
}

CPPAD_END_NAMESPACE
# endif

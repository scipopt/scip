/* $Id: atan_op.hpp 1641 2010-02-01 16:39:45Z bradbell $ */
# ifndef CPPAD_ATAN_OP_INCLUDED
# define CPPAD_ATAN_OP_INCLUDED

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
\file atan_op.hpp
Forward and reverse mode calculations for z = atan(x).
*/


/*!
Forward mode Taylor coefficient for result of op = AtanOp.

The C++ source code corresponding to this operation is
\verbatim
	z = atan(x)
\endverbatim
The auxillary result is
\verbatim
	y = 1 + x * x
\endverbatim
The value of y, and its derivatives, are computed along with the value
and derivatives of z.

\copydetails forward_unary2_op
*/
template <class Base>
inline void forward_atan_op(
	size_t j           ,
	size_t i_z         ,
	size_t i_x         ,
	size_t nc_taylor   , 
	Base*  taylor      )
{	
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(AtanOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( NumRes(AtanOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( i_x + 1 < i_z );
	CPPAD_ASSERT_UNKNOWN( j < nc_taylor );

	// Taylor coefficients corresponding to argument and result
	Base* x = taylor + i_x * nc_taylor;
	Base* z = taylor + i_z * nc_taylor;
	Base* b = z      -       nc_taylor;  // called y in documentation

	size_t k;
	if( j == 0 )
	{	z[j] = atan( x[0] );
		b[j] = Base(1) + x[0] * x[0];
	}
	else
	{
		b[j] = Base(2) * x[0] * x[j];
		z[j] = Base(0);
		for(k = 1; k < j; k++)
		{	b[j] += x[k] * x[j-k];
			z[j] -= Base(k) * z[k] * b[j-k];
		}
		z[j] /= Base(j);
		z[j] += x[j];
		z[j] /= b[0];
	}
}

/*!
Zero order forward mode Taylor coefficient for result of op = AtanOp.

The C++ source code corresponding to this operation is
\verbatim
	z = atan(x)
\endverbatim
The auxillary result is
\verbatim
	y = 1 + x * x
\endverbatim
The value of y is computed along with the value of z.

\copydetails forward_unary2_op_0
*/
template <class Base>
inline void forward_atan_op_0(
	size_t i_z         ,
	size_t i_x         ,
	size_t nc_taylor   , 
	Base*  taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(AtanOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( NumRes(AtanOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( i_x + 1 < i_z );
	CPPAD_ASSERT_UNKNOWN( 0 < nc_taylor );

	// Taylor coefficients corresponding to argument and result
	Base* x = taylor + i_x * nc_taylor;
	Base* z = taylor + i_z * nc_taylor;
	Base* b = z      -       nc_taylor; // called y in documentation

	z[0] = atan( x[0] );
	b[0] = Base(1) + x[0] * x[0];
}
/*!
Reverse mode partial derivatives for result of op = AtanOp.

The C++ source code corresponding to this operation is
\verbatim
	z = atan(x)
\endverbatim
The auxillary result is
\verbatim
	y = 1 + x * x
\endverbatim
The value of y is computed along with the value of z.

\copydetails reverse_unary2_op
*/

template <class Base>
inline void reverse_atan_op(
	size_t      d            ,
	size_t      i_z          ,
	size_t      i_x          ,
	size_t      nc_taylor    , 
	const Base* taylor       ,
	size_t      nc_partial   ,
	Base*       partial      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(AtanOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( NumRes(AtanOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( i_x + 1 < i_z );
	CPPAD_ASSERT_UNKNOWN( d < nc_taylor );
	CPPAD_ASSERT_UNKNOWN( d < nc_partial );

	// Taylor coefficients and partials corresponding to argument
	const Base* x  = taylor  + i_x * nc_taylor;
	Base* px       = partial + i_x * nc_partial;

	// Taylor coefficients and partials corresponding to first result
	const Base* z  = taylor  + i_z * nc_taylor;
	Base* pz       = partial + i_z * nc_partial;

	// Taylor coefficients and partials corresponding to auxillary result
	const Base* b  = z  - nc_taylor; // called y in documentation
	Base* pb       = pz - nc_partial;

	// number of indices to access
	size_t j = d;
	size_t k;
	while(j)
	{	// scale partials w.r.t z[j] and b[j]
		pz[j] /= b[0];
		pb[j] *= Base(2);

		pb[0] -= pz[j] * z[j]; 
		px[j] += pz[j] + pb[j] * x[0];
		px[0] += pb[j] * x[j];

		// more scaling of partials w.r.t z[j]
		pz[j] /= Base(j);

		for(k = 1; k < j; k++)
		{	pb[j-k] -= pz[j] * Base(k) * z[k];
			pz[k]   -= pz[j] * Base(k) * b[j-k];
			px[k]   += pb[j] * x[j-k];
		}
		--j;
	}
	px[0] += pz[0] / b[0] + pb[0] * Base(2) * x[0];
}

CPPAD_END_NAMESPACE
# endif

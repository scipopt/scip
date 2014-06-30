/* $Id: sub_op.hpp 2910 2013-10-07 13:27:58Z bradbell $ */
# ifndef CPPAD_SUB_OP_INCLUDED
# define CPPAD_SUB_OP_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-13 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Eclipse Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

namespace CppAD { // BEGIN_CPPAD_NAMESPACE
/*!
\defgroup sub_op_hpp sub_op.hpp
\{
\file sub_op.hpp
Forward and reverse mode calculations for z = x - y.
*/

// --------------------------- Subvv -----------------------------------------
/*!
Compute forward mode Taylor coefficients for result of op = SubvvOp.

The C++ source code corresponding to this operation is
\verbatim
	z = x - y
\endverbatim
In the documentation below,
this operations is for the case where both x and y are variables
and the argument \a parameter is not used.

\copydetails forward_binary_op
*/

template <class Base>
inline void forward_subvv_op(
	size_t        q           , 
	size_t        p           , 
	size_t        i_z         ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        nc_taylor   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(SubvvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(SubvvOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < i_z );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[1]) < i_z );
	CPPAD_ASSERT_UNKNOWN( p < nc_taylor );
	CPPAD_ASSERT_UNKNOWN( q <= p );

	// Taylor coefficients corresponding to arguments and result
	Base* x = taylor + arg[0] * nc_taylor;
	Base* y = taylor + arg[1] * nc_taylor;
	Base* z = taylor + i_z    * nc_taylor;

	for(size_t d = q; d <= p; d++)
		z[d] = x[d] - y[d];
}

/*!
Compute zero order forward mode Taylor coefficients for result of op = SubvvOp.

The C++ source code corresponding to this operation is
\verbatim
	z = x - y
\endverbatim
In the documentation below,
this operations is for the case where both x and y are variables
and the argument \a parameter is not used.

\copydetails forward_binary_op_0
*/

template <class Base>
inline void forward_subvv_op_0(
	size_t        i_z         ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        nc_taylor   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(SubvvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(SubvvOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < i_z );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[1]) < i_z );

	// Taylor coefficients corresponding to arguments and result
	Base* x = taylor + arg[0] * nc_taylor;
	Base* y = taylor + arg[1] * nc_taylor;
	Base* z = taylor + i_z    * nc_taylor;

	z[0] = x[0] - y[0];
}

/*!
Compute reverse mode partial derivatives for result of op = SubvvOp.

The C++ source code corresponding to this operation is
\verbatim
	z = x - y
\endverbatim
In the documentation below,
this operations is for the case where both x and y are variables
and the argument \a parameter is not used.

\copydetails reverse_binary_op
*/

template <class Base>
inline void reverse_subvv_op(
	size_t        d           , 
	size_t        i_z         ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        nc_taylor   ,
	const Base*   taylor      ,
	size_t        nc_partial  ,
	Base*         partial     )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(SubvvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(SubvvOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < i_z );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[1]) < i_z );
	CPPAD_ASSERT_UNKNOWN( d < nc_taylor );
	CPPAD_ASSERT_UNKNOWN( d < nc_partial );

	// Partial derivatives corresponding to arguments and result
	Base* px = partial + arg[0] * nc_partial;
	Base* py = partial + arg[1] * nc_partial;
	Base* pz = partial + i_z    * nc_partial;

	// number of indices to access
	size_t i = d + 1;
	while(i)
	{	--i;
		px[i] += pz[i];
		py[i] -= pz[i];
	}
}

// --------------------------- Subpv -----------------------------------------
/*!
Compute forward mode Taylor coefficients for result of op = SubpvOp.

The C++ source code corresponding to this operation is
\verbatim
	z = x - y
\endverbatim
In the documentation below,
this operations is for the case where x is a parameter and y is a variable.

\copydetails forward_binary_op
*/

template <class Base>
inline void forward_subpv_op(
	size_t        q           , 
	size_t        p           , 
	size_t        i_z         ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        nc_taylor   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(SubpvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(SubpvOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[1]) < i_z );
	CPPAD_ASSERT_UNKNOWN( p < nc_taylor );
	CPPAD_ASSERT_UNKNOWN( q <= p );

	// Taylor coefficients corresponding to arguments and result
	Base* y = taylor + arg[1] * nc_taylor;
	Base* z = taylor + i_z    * nc_taylor;

	// Paraemter value
	Base x = parameter[ arg[0] ];
	if( q == 0 )
	{	z[0] = x - y[0];
		q++;
	}
	for(size_t d = q; d <= p; d++)
		z[d] = - y[d];
}
/*!
Compute zero order forward mode Taylor coefficient for result of op = SubpvOp.

The C++ source code corresponding to this operation is
\verbatim
	z = x - y
\endverbatim
In the documentation below,
this operations is for the case where x is a parameter and y is a variable.

\copydetails forward_binary_op_0
*/

template <class Base>
inline void forward_subpv_op_0(
	size_t        i_z         ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        nc_taylor   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(SubpvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(SubpvOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[1]) < i_z );

	// Paraemter value
	Base x = parameter[ arg[0] ];

	// Taylor coefficients corresponding to arguments and result
	Base* y = taylor + arg[1] * nc_taylor;
	Base* z = taylor + i_z    * nc_taylor;

	z[0] = x - y[0];
}

/*!
Compute reverse mode partial derivative for result of op = SubpvOp.

The C++ source code corresponding to this operation is
\verbatim
	z = x - y
\endverbatim
In the documentation below,
this operations is for the case where x is a parameter and y is a variable.

\copydetails reverse_binary_op
*/

template <class Base>
inline void reverse_subpv_op(
	size_t        d           , 
	size_t        i_z         ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        nc_taylor   ,
	const Base*   taylor      ,
	size_t        nc_partial  ,
	Base*         partial     )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(SubvvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(SubvvOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[1]) < i_z );
	CPPAD_ASSERT_UNKNOWN( d < nc_taylor );
	CPPAD_ASSERT_UNKNOWN( d < nc_partial );

	// Partial derivatives corresponding to arguments and result
	Base* py = partial + arg[1] * nc_partial;
	Base* pz = partial + i_z    * nc_partial;

	// number of indices to access
	size_t i = d + 1;
	while(i)
	{	--i;
		py[i] -= pz[i];
	}
}

// --------------------------- Subvp -----------------------------------------
/*!
Compute forward mode Taylor coefficients for result of op = SubvvOp.

The C++ source code corresponding to this operation is
\verbatim
	z = x - y
\endverbatim
In the documentation below,
this operations is for the case where x is a variable and y is a parameter.

\copydetails forward_binary_op
*/

template <class Base>
inline void forward_subvp_op(
	size_t        q           , 
	size_t        p           , 
	size_t        i_z         ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        nc_taylor   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(SubvpOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(SubvpOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < i_z );
	CPPAD_ASSERT_UNKNOWN( p < nc_taylor );
	CPPAD_ASSERT_UNKNOWN( q <= p );

	// Taylor coefficients corresponding to arguments and result
	Base* x = taylor + arg[0] * nc_taylor;
	Base* z = taylor + i_z    * nc_taylor;

	// Parameter value
	Base y = parameter[ arg[1] ];
	if( q == 0 )
	{	z[0] = x[0] - y;
		q++;
	}
	for(size_t d = q; d <= p; d++)
		z[d] = x[d];
}

/*!
Compute zero order forward mode Taylor coefficients for result of op = SubvvOp.

The C++ source code corresponding to this operation is
\verbatim
	z = x - y
\endverbatim
In the documentation below,
this operations is for the case where x is a variable and y is a parameter.

\copydetails forward_binary_op_0
*/

template <class Base>
inline void forward_subvp_op_0(
	size_t        i_z         ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        nc_taylor   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(SubvpOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(SubvpOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < i_z );

	// Parameter value
	Base y = parameter[ arg[1] ];

	// Taylor coefficients corresponding to arguments and result
	Base* x = taylor + arg[0] * nc_taylor;
	Base* z = taylor + i_z    * nc_taylor;

	z[0] = x[0] - y;
}

/*!
Compute reverse mode partial derivative for result of op = SubvpOp.

The C++ source code corresponding to this operation is
\verbatim
	z = x - y
\endverbatim
In the documentation below,
this operations is for the case where x is a variable and y is a parameter.

\copydetails reverse_binary_op
*/

template <class Base>
inline void reverse_subvp_op(
	size_t        d           , 
	size_t        i_z         ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        nc_taylor   ,
	const Base*   taylor      ,
	size_t        nc_partial  ,
	Base*         partial     )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(SubvpOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(SubvpOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < i_z );
	CPPAD_ASSERT_UNKNOWN( d < nc_taylor );
	CPPAD_ASSERT_UNKNOWN( d < nc_partial );

	// Partial derivatives corresponding to arguments and result
	Base* px = partial + arg[0] * nc_partial;
	Base* pz = partial + i_z    * nc_partial;

	// number of indices to access
	size_t i = d + 1;
	while(i)
	{	--i;
		px[i] += pz[i];
	}
}

/*! \} */
} // END_CPPAD_NAMESPACE
# endif

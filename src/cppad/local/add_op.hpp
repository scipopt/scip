/* $Id: add_op.hpp 1986 2011-06-18 20:33:17Z bradbell $ */
# ifndef CPPAD_ADD_OP_INCLUDED
# define CPPAD_ADD_OP_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

CPPAD_BEGIN_NAMESPACE
/*!
\file add_op.hpp
Forward and reverse mode calculations for z = x + y.
*/

// --------------------------- Addvv -----------------------------------------
/*!
Compute forward mode Taylor coefficients for result of op = AddvvOp.

The C++ source code corresponding to this operation is
\verbatim
	z = x + y
\endverbatim
In the documentation below,
this operations is for the case where both x and y are variables
and the argument \a parameter is not used.

\copydetails forward_binary_op
*/

template <class Base>
inline void forward_addvv_op(
	size_t        d           , 
	size_t        i_z         ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        nc_taylor   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(AddvvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(AddvvOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < i_z );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[1]) < i_z );
	CPPAD_ASSERT_UNKNOWN( d < nc_taylor );

	// Taylor coefficients corresponding to arguments and result
	Base* x = taylor + arg[0] * nc_taylor;
	Base* y = taylor + arg[1] * nc_taylor;
	Base* z = taylor + i_z    * nc_taylor;

	z[d] = x[d] + y[d];
}

/*!
Compute zero order forward mode Taylor coefficients for result of op = AddvvOp.

The C++ source code corresponding to this operation is
\verbatim
	z = x + y
\endverbatim
In the documentation below,
this operations is for the case where both x and y are variables
and the argument \a parameter is not used.

\copydetails forward_binary_op_0
*/

template <class Base>
inline void forward_addvv_op_0(
	size_t        i_z         ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        nc_taylor   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(AddvvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(AddvvOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < i_z );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[1]) < i_z );

	// Taylor coefficients corresponding to arguments and result
	Base* x = taylor + arg[0] * nc_taylor;
	Base* y = taylor + arg[1] * nc_taylor;
	Base* z = taylor + i_z    * nc_taylor;

	z[0] = x[0] + y[0];
}

/*!
Compute reverse mode partial derivatives for result of op = AddvvOp.

The C++ source code corresponding to this operation is
\verbatim
	z = x + y
\endverbatim
In the documentation below,
this operations is for the case where both x and y are variables
and the argument \a parameter is not used.

\copydetails reverse_binary_op
*/

template <class Base>
inline void reverse_addvv_op(
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
	CPPAD_ASSERT_UNKNOWN( NumArg(AddvvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(AddvvOp) == 1 );
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
		py[i] += pz[i];
	}
}

// --------------------------- Addpv -----------------------------------------
/*!
Compute forward mode Taylor coefficients for result of op = AddpvOp.

The C++ source code corresponding to this operation is
\verbatim
	z = x + y
\endverbatim
In the documentation below,
this operations is for the case where x is a parameter and y is a variable.

\copydetails forward_binary_op
*/

template <class Base>
inline void forward_addpv_op(
	size_t        d           , 
	size_t        i_z         ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        nc_taylor   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(AddpvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(AddpvOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[1]) < i_z );
	CPPAD_ASSERT_UNKNOWN( d < nc_taylor );

	// Taylor coefficients corresponding to arguments and result
	Base* y = taylor + arg[1] * nc_taylor;
	Base* z = taylor + i_z    * nc_taylor;

# if CPPAD_USE_FORWARD0SWEEP
	CPPAD_ASSERT_UNKNOWN( d > 0 );
	z[d] = y[d];
# else
	// Paraemter value
	Base x = parameter[ arg[0] ];
	if( d == 0 )
		z[d] = x + y[d];
	else	z[d] = y[d];
# endif
}
/*!
Compute zero order forward mode Taylor coefficient for result of op = AddpvOp.

The C++ source code corresponding to this operation is
\verbatim
	z = x + y
\endverbatim
In the documentation below,
this operations is for the case where x is a parameter and y is a variable.

\copydetails forward_binary_op_0
*/

template <class Base>
inline void forward_addpv_op_0(
	size_t        i_z         ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        nc_taylor   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(AddpvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(AddpvOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[1]) < i_z );

	// Paraemter value
	Base x = parameter[ arg[0] ];

	// Taylor coefficients corresponding to arguments and result
	Base* y = taylor + arg[1] * nc_taylor;
	Base* z = taylor + i_z    * nc_taylor;

	z[0] = x + y[0];
}

/*!
Compute reverse mode partial derivative for result of op = AddpvOp.

The C++ source code corresponding to this operation is
\verbatim
	z = x + y
\endverbatim
In the documentation below,
this operations is for the case where x is a parameter and y is a variable.

\copydetails reverse_binary_op
*/

template <class Base>
inline void reverse_addpv_op(
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
	CPPAD_ASSERT_UNKNOWN( NumArg(AddvvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(AddvvOp) == 1 );
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
		py[i] += pz[i];
	}
}


CPPAD_END_NAMESPACE
# endif

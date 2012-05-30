/* $Id: discrete_op.hpp 1986 2011-06-18 20:33:17Z bradbell $ */
# ifndef CPPAD_DISCRETE_OP_INCLUDED
# define CPPAD_DISCRETE_OP_INCLUDED

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
\file discrete_op.hpp
Zero order forward mode for z = f(x) where f is piecewise constant.
*/


/*!
Compute zero order forward mode Taylor coefficient for result of op = DisOp.

The C++ source code corresponding to this operation is
\verbatim
	z = f(x)
\endverbatim
where f is a piecewise constant function (and it's derivative is always
calculated as zero).

\tparam Base
base type for the operator; i.e., this operation was recorded
using AD< \a Base > and computations by this routine are done using type 
\a Base .

\param i_z
variable index corresponding to the result for this operation; 
i.e. the row index in \a taylor corresponding to z. 

\param arg
\a arg[0]
\n
is the index, in the order of the discrete functions defined by the user, 
for this discrete function.
\n
\n
\a arg[1]
variable index corresponding to the argument for this operator;
i.e. the row index in \a taylor corresponding to x.

\param nc_taylor
number of colums in the matrix containing all the Taylor coefficients.

\param taylor
\b Input: \a taylor [ \a arg[1] * \a nc_taylor + 0 ] 
is the zero order Taylor coefficient corresponding to x. 
\n
\b Output: \a taylor [ \a i_z * \a nc_taylor + 0 ] 
is the zero order Taylor coefficient corresponding to z. 

\par Checked Assertions where op is the unary operator with one result:
\li NumArg(op) == 2
\li NumRes(op) == 1
\li \a arg[1] < \a i_z 
\li \a 0 < \a nc_taylor
*/
template <class Base>
inline void forward_dis_op_0(
	size_t        i_z         ,
	const addr_t* arg         ,
	size_t        nc_taylor   , 
	Base*         taylor      )
{	
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(DisOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(DisOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[1]) < i_z );
	CPPAD_ASSERT_UNKNOWN( 0 < nc_taylor );

	// Taylor coefficients corresponding to argument and result
	Base* x = taylor + arg[1] * nc_taylor;
	Base* z = taylor + i_z * nc_taylor;

	z[0]  = discrete<Base>::eval(arg[0], x[0]);
}


CPPAD_END_NAMESPACE
# endif

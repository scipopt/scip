/* $Id: comp_op.hpp 1986 2011-06-18 20:33:17Z bradbell $ */
# ifndef CPPAD_COMP_OP_INCLUDED
# define CPPAD_COMP_OP_INCLUDED

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
\file comp_op.hpp
Zero order forward mode check how man comparisons changed.
*/

/*!
Zero order forward mode execution of op = CompOp.

The C++ source code corresponding to this operation is
\verbatim
	left Cop right
\endverbatim
where Cop is one of the following:
<, <=, == , >=, >, !=.

\tparam Base
base type for the operator; i.e., this operation was recorded
using AD< \a Base > and computations by this routine are done using type 
\a Base.

\param count
If the comparision has the same result as when t operation seqeunce was
recorded, \a count is not affected.
Otherwise it is incremented by one.

\param arg
\n
\a arg[0]
is static cast to size_t from the enum type
\verbatim
	enum CompareOp {
		CompareLt, 
		CompareLe, 
		CompareEq, 
		CompareGe, 
		CompareGt, 
		CompareNe
	}
\endverbatim
for this operation.
\n
\n 
\a arg[1] & 1 
\n
If this expression is true, 
the result of the comparison during taping it true.
Othwise the result if false.
\n
\n
\a arg[1] & 2
\n
if this expression is true, left is a variable, otherwise it is a parameter.
\n
\n
\a arg[1] & 4
\n
if this expression is true, right is a variable, otherwise it is a parameter.
\n

\param num_par
is the lenght of the \a parameter vector.
This value is only used for checking assumptions mentioned below
(and is not used at all when NDEBUG is defined).

\param parameter
Vector of parameters corresponding to the tape.
If left is a parameter, \a parameter[ arg[2] ] is its value.
If right is a parameter, \a parameter[ arg[3] ] is its value.

\param nc_taylor
number of columns in the matrix containing the Taylor coefficients.

\param taylor
Matrix of Taylor coefficients.
If left is a variable, \a taylor[ arg[2] * nc_taylor + 0 ] is its value.
If right is a variable, \a taylor[ arg[3] * nc_taylor + 0 ] is its value.


\par Checked Assertions where op is a binary operator:
\li NumArg(ComOp) == 4
\li NumRes(ComOp) == 0
\li size_t(arg[0]) <= static_cast<size_t>( CompareNe )
\li arg[1] != 0 (either left or right is a variable)
\li if left is a parameter, \a arg[2] < \a num_par
\li if right is a parameter, \a arg[3] < \a num_par

*/
template <class Base>
inline void forward_comp_op_0(
	size_t&        count       ,
	const addr_t*  arg         , 
	size_t         num_par     ,
	const Base*    parameter   ,
	size_t         nc_taylor   ,
	Base*          taylor      )
{	bool result;
	Base left;
	Base right;

	CPPAD_ASSERT_UNKNOWN( NumArg(ComOp) == 4 );
	CPPAD_ASSERT_UNKNOWN( NumRes(ComOp) == 0 );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) <= static_cast<size_t> (CompareNe) );
	CPPAD_ASSERT_UNKNOWN( arg[1] != 0 );

	// result of comparision during recording
	result = (arg[1] & 1) == 1;

	// value of left operand for this forward sweep
	if( arg[1] & 2 )
		left = taylor[ arg[2] * nc_taylor + 0 ];
	else
	{	CPPAD_ASSERT_UNKNOWN( size_t(arg[2]) < num_par );
		left = parameter[ arg[2] ];
	}

	// value of right operand for this forward sweep.
	if( arg[1] & 4 )
		right = taylor[ arg[3] * nc_taylor + 0 ];
	else
	{	CPPAD_ASSERT_UNKNOWN( size_t(arg[3]) < num_par );
		right = parameter[ arg[3] ];
	}
	switch( CompareOp( arg[0] ) )
	{	case CompareLt:
		count += ( result != LessThanZero(left - right) );
		break;

		case CompareLe:
		count += ( result != LessThanOrZero(left - right) );
		break;

		case CompareEq:
		count += ( result != (left == right) );
		break;

		case CompareGe:
		count += ( result != GreaterThanOrZero(left - right) );
		break;

		case CompareGt:
		count += ( result != GreaterThanZero(left - right) );
		break;

		case CompareNe:
		count += ( result != (left != right) );
		break;

		default:
		CPPAD_ASSERT_UNKNOWN(0);
	}
	return;
}
CPPAD_END_NAMESPACE
# endif

/* $Id: sub_eq.hpp 1986 2011-06-18 20:33:17Z bradbell $ */
# ifndef CPPAD_SUB_EQ_INCLUDED
# define CPPAD_SUB_EQ_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

//  BEGIN CppAD namespace
namespace CppAD {

template <class Base>
AD<Base>& AD<Base>::operator -= (const AD<Base> &right)
{	ADTape<Base> *tape = AD<Base>::tape_ptr();
	size_t tape_id = 0;
	if( tape != CPPAD_NULL )
		tape_id = tape->id_;

	// id_ setting for parameters cannot match 0
	bool var_left  = id_       == tape_id;
	bool var_right = right.id_ == tape_id;
	CPPAD_ASSERT_KNOWN(
		Parameter(*this) || var_left ,
		"-=: left operand is a variable for a different thread"
	);
	CPPAD_ASSERT_KNOWN(
		Parameter(right) || var_right ,
		"-=: right operand is a variable for a different thread"
	);

	Base left;
	left    = value_;
	value_ -= right.value_;

	if( var_left )
	{	if( var_right )
		{	// this = variable - variable
			CPPAD_ASSERT_UNKNOWN( NumRes(SubvvOp) == 1 );
			CPPAD_ASSERT_UNKNOWN( NumArg(SubvvOp) == 2 );

			// put operand addresses in tape
			tape->Rec_.PutArg(taddr_, right.taddr_);
			// put operator in the tape
			taddr_ = tape->Rec_.PutOp(SubvvOp);
			// make this a variable
			CPPAD_ASSERT_UNKNOWN( id_ == tape_id );
		}
		else if( IdenticalZero( right.value_ ) )
		{	// this = variable - 0
		}
		else
		{	// this = variable - parameter
			CPPAD_ASSERT_UNKNOWN( NumRes(SubvpOp) == 1 );
			CPPAD_ASSERT_UNKNOWN( NumArg(SubvpOp) == 2 );

			// put operand addresses in tape
			addr_t p = tape->Rec_.PutPar(right.value_);
			tape->Rec_.PutArg(taddr_, p);
			// put operator in the tape
			taddr_ = tape->Rec_.PutOp(SubvpOp);
			// make this a variable
			CPPAD_ASSERT_UNKNOWN( id_ == tape_id );
		}
	}
	else if( var_right  )
	{	// this = parameter - variable
		CPPAD_ASSERT_UNKNOWN( NumRes(SubpvOp) == 1 );
		CPPAD_ASSERT_UNKNOWN( NumArg(SubpvOp) == 2 );

		// put operand addresses in tape
		addr_t p = tape->Rec_.PutPar(left);
		tape->Rec_.PutArg(p, right.taddr_);
		// put operator in the tape
		taddr_ = tape->Rec_.PutOp(SubpvOp);
		// make this a variable
		id_ = tape_id;
	}
	return *this;
}

CPPAD_FOLD_ASSIGNMENT_OPERATOR(-=)

} // END CppAD namespace

# endif 

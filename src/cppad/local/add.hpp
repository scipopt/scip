/* $Id: add.hpp 1986 2011-06-18 20:33:17Z bradbell $ */
# ifndef CPPAD_ADD_INCLUDED
# define CPPAD_ADD_INCLUDED

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
AD<Base> operator + (const AD<Base> &left , const AD<Base> &right)
{	ADTape<Base> *tape = AD<Base>::tape_ptr();
	size_t tape_id = 0;
	if( tape != CPPAD_NULL )
		tape_id = tape->id_;

	// id_ setting for parameters cannot match 0
	bool var_left  = left.id_  == tape_id;
	bool var_right = right.id_ == tape_id;
	CPPAD_ASSERT_KNOWN(
		Parameter(left) || var_left ,
		"+: left operand is a variable for a different thread"
	);
	CPPAD_ASSERT_KNOWN(
		Parameter(right) || var_right ,
		"+: right operand is a variable for a different thread"
	);

	AD<Base> result;
	result.value_  = left.value_ + right.value_;
	CPPAD_ASSERT_UNKNOWN( Parameter(result) );

	if( var_left )
	{	if( var_right )
		{	// result = variable + variable
			CPPAD_ASSERT_UNKNOWN( NumRes(AddvvOp) == 1 );
			CPPAD_ASSERT_UNKNOWN( NumArg(AddvvOp) == 2 );

			// put operand addresses in tape
			tape->Rec_.PutArg(left.taddr_, right.taddr_);
			// put operator in the tape
			result.taddr_ = tape->Rec_.PutOp(AddvvOp);
			// make result a variable
			result.id_ = tape_id;
		}
		else if( IdenticalZero(right.value_) )
		{	// result = variable + 0
			result.make_variable(left.id_, left.taddr_);
		}
		else
		{	// result = variable  + parameter
			//        = parameter + variable
			CPPAD_ASSERT_UNKNOWN( NumRes(AddpvOp) == 1 );
			CPPAD_ASSERT_UNKNOWN( NumArg(AddpvOp) == 2 );

			// put operand addresses in tape
			addr_t p = tape->Rec_.PutPar(right.value_);
			tape->Rec_.PutArg(p, left.taddr_);
			// put operator in the tape
			result.taddr_ = tape->Rec_.PutOp(AddpvOp);
			// make result a variable
			result.id_ = tape_id;
		}
	}
	else if( var_right )
	{	if( IdenticalZero(left.value_) )
		{	// result = 0 + variable
			result.make_variable(right.id_, right.taddr_);
		}
		else
		{	// result = parameter + variable
			CPPAD_ASSERT_UNKNOWN( NumRes(AddpvOp) == 1 );
			CPPAD_ASSERT_UNKNOWN( NumArg(AddpvOp) == 2 );

			// put operand addresses in tape
			addr_t p = tape->Rec_.PutPar(left.value_);
			tape->Rec_.PutArg(p, right.taddr_);
			// put operator in the tape
			result.taddr_ = tape->Rec_.PutOp(AddpvOp);
			// make result a variable
			result.id_ = tape_id;
		}
	}
	return result;
}

// convert other cases into the case above
CPPAD_FOLD_AD_VALUED_BINARY_OPERATOR(+)

} // END CppAD namespace

# endif 

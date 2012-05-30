/* $Id: pow.hpp 2057 2011-08-11 14:07:11Z bradbell $ */
# ifndef CPPAD_POW_INCLUDED
# define CPPAD_POW_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin pow$$
$spell
	Vec
	std
	namespace
	CppAD
	const
$$

$index pow, AD$$
$index exponent, AD function$$

$section The AD Power Function$$

$head Syntax$$
$syntax%%z% = pow(%x%, %y%)%$$


$head Purpose$$
Determines the value of the power function which is defined by
$latex \[
	{\rm pow} (x, y) = x^y
\] $$
This version of the $code pow$$ function may use
logarithms and exponentiation to compute derivatives.
This will not work if $italic x$$ is less than or equal zero.
If the value of $italic y$$ is an integer, 
the $cref/pow_int/$$ function is used to compute this value 
using only multiplication (and division if $italic y$$ is negative). 
(This will work even if $italic x$$ is less than or equal zero.)

$head x$$
The argument $italic x$$ has the following prototype
$syntax%
	const %Type% &%x%
%$$
where $italic Type$$ is
$syntax%VecAD<%Base%>::reference%$$,
$syntax%AD<%Base%>%$$,
$syntax%%Base%$$,
$syntax%double%$$,
or
$syntax%int%$$.

$head y$$
The argument $italic y$$ has the following prototype
$syntax%
	const %Type% &%y%
%$$
where $italic Type$$ is
$syntax%VecAD<%Base%>::reference%$$,
$syntax%AD<%Base%>%$$,
$syntax%%Base%$$,
$syntax%double%$$,
or
$syntax%int%$$.

$head z$$
The result $italic z$$ has prototype
$syntax%
	AD<%Base%> %z%
%$$

$head Standard Types$$
A definition for the $code pow$$ function is included
in the CppAD namespace for the case where both $italic x$$
and $italic y$$ have the same type and that type is
$code float$$ or $code double$$.
 
$head Operation Sequence$$
This is an AD of $italic Base$$
$xref/glossary/Operation/Atomic/atomic operation/1/$$
and hence is part of the current
AD of $italic Base$$
$xref/glossary/Operation/Sequence/operation sequence/1/$$.

$head Example$$
$children%
	example/pow.cpp%
	example/pow_int.cpp
%$$
The files
$xref/Pow.cpp/$$, $cref/pow_int.cpp/$$
contain an examples and tests of this function.   
They returns true if they succeed and false otherwise.

$end
-------------------------------------------------------------------------------
*/

//  BEGIN CppAD namespace
namespace CppAD {
 
// case where x and y are AD<Base> -----------------------------------------
template <class Base> AD<Base> 
pow(const AD<Base> &x, const AD<Base> &y)
{	ADTape<Base> *tape = AD<Base>::tape_ptr();
	size_t tape_id = 0;
	if( tape != CPPAD_NULL )
		tape_id = tape->id_;

	// id_ setting for parameters cannot match 0
	bool var_x = x.id_  == tape_id;
	bool var_y = y.id_ == tape_id;
	CPPAD_ASSERT_KNOWN(
		Parameter(x) || var_x ,
		"pow: first argument is a variable for a different thread"
	);
	CPPAD_ASSERT_KNOWN(
		Parameter(y) || var_y ,
		"pow: second argument is a variable for a different thread"
	);

	AD<Base> result;
	result.value_  = pow(x.value_, y.value_);
	CPPAD_ASSERT_UNKNOWN( Parameter(result) );

	if( var_x )
	{	if( var_y )
		{	// result = variable^variable
			CPPAD_ASSERT_UNKNOWN( NumRes(PowvvOp) == 3 );
			CPPAD_ASSERT_UNKNOWN( NumArg(PowvvOp) == 2 );

			// put operand addresses in tape
			tape->Rec_.PutArg(x.taddr_, y.taddr_);

			// put operator in the tape
			result.taddr_ = tape->Rec_.PutOp(PowvvOp);

			// make result a variable
			result.id_ = tape_id;
		}
		else if( IdenticalZero( y.value_ ) )
		{	// result = variable^0
		}
		else
		{	// result = variable^parameter 
			CPPAD_ASSERT_UNKNOWN( NumRes(PowvpOp) == 3 );
			CPPAD_ASSERT_UNKNOWN( NumArg(PowvpOp) == 2 );

			// put operand addresses in tape
			addr_t p = tape->Rec_.PutPar(y.value_);
			tape->Rec_.PutArg(x.taddr_, p);

			// put operator in the tape
			result.taddr_ = tape->Rec_.PutOp(PowvpOp);

			// make result a variable
			result.id_ = tape_id;
		}
	}
	else if( var_y )
	{	if( IdenticalZero(x.value_) )
		{	// result = 0^variable 
		}
		else
		{	// result = variable^parameter 
			CPPAD_ASSERT_UNKNOWN( NumRes(PowpvOp) == 3 );
			CPPAD_ASSERT_UNKNOWN( NumArg(PowpvOp) == 2 );

			// put operand addresses in tape
			addr_t p = tape->Rec_.PutPar(x.value_);
			tape->Rec_.PutArg(p, y.taddr_);

			// put operator in the tape
			result.taddr_ = tape->Rec_.PutOp(PowpvOp);

			// make result a variable
			result.id_ = tape_id;
		}
	}
	return result;
}
// =========================================================================
// Fold operations in same way as CPPAD_FOLD_AD_VALUED_BINARY_OPERATION(Op)
// -------------------------------------------------------------------------
// Operations with VecAD_reference<Base> and AD<Base> only

template <class Base> AD<Base>
pow(const AD<Base> &x, const VecAD_reference<Base> &y)
{	return pow(x, y.ADBase()); }

template <class Base> AD<Base> 
pow(const VecAD_reference<Base> &x, const VecAD_reference<Base> &y) 
{	return pow(x.ADBase(), y.ADBase()); }

template <class Base> AD<Base>
pow(const VecAD_reference<Base> &x, const AD<Base> &y)
{	return pow(x.ADBase(), y); }
// -------------------------------------------------------------------------
// Operations with Base

template <class Base> AD<Base>
pow(const Base &x, const AD<Base> &y)
{	return pow(AD<Base>(x), y); }

template <class Base> AD<Base>
pow(const Base &x, const VecAD_reference<Base> &y)
{	return pow(AD<Base>(x), y.ADBase()); }

template <class Base> AD<Base>
pow(const AD<Base> &x, const Base &y)
{	return pow(x, AD<Base>(y)); }

template <class Base> AD<Base>
pow(const VecAD_reference<Base> &x, const Base &y)
{	return pow(x.ADBase(), AD<Base>(y)); }
// -------------------------------------------------------------------------
// Operations with double

template <class Base> AD<Base>
pow(const double &x, const AD<Base> &y)
{	return pow(AD<Base>(x), y); }

template <class Base> AD<Base>
pow(const double &x, const VecAD_reference<Base> &y)
{	return pow(AD<Base>(x), y.ADBase()); }

template <class Base> AD<Base>
pow(const AD<Base> &x, const double &y)
{	return pow(x, AD<Base>(y)); }

template <class Base> AD<Base>
pow(const VecAD_reference<Base> &x, const double &y)
{	return pow(x.ADBase(), AD<Base>(y)); }
// -------------------------------------------------------------------------
// Special case to avoid ambuigity when Base is double

inline AD<double>
pow(const double &x, const AD<double> &y)
{	return pow(AD<double>(x), y); }

inline AD<double>
pow(const double &x, const VecAD_reference<double> &y)
{	return pow(AD<double>(x), y.ADBase()); }

inline AD<double>
pow(const AD<double> &x, const double &y)
{	return pow(x, AD<double>(y)); }

inline AD<double>
pow(const VecAD_reference<double> &x, const double &y)
{	return pow(x.ADBase(), AD<double>(y)); }

// =========================================================================
// Fold operations for the cases where x is an int, 
// but let cppad/pow_int.hpp handle the cases where y is an int.
// -------------------------------------------------------------------------
template <class Base> AD<Base> pow
(const int &x, const VecAD_reference<Base> &y)
{	return pow(AD<Base>(x), y.ADBase()); }

template <class Base> AD<Base> pow
(const int &x, const AD<Base> &y)
{	return pow(AD<Base>(x), y); }

} // END CppAD namespace

# endif 

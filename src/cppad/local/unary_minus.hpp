/* $Id: unary_minus.hpp 1369 2009-05-31 01:31:48Z bradbell $ */
# ifndef CPPAD_UNARY_MINUS_INCLUDED
# define CPPAD_UNARY_MINUS_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-06 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin UnaryMinus$$
$spell
	Vec
	const
	inline
$$

$index unary, AD minus operator$$
$index AD, unary minus operator$$
$index minus, AD unary operator$$
$index -, AD unary operator$$

$section AD Unary Minus Operator$$

$head Syntax$$

$syntax%%y% = - %x%$$


$head Purpose$$
Computes the negative of $italic x$$.

$head Base$$
The operation in the syntax above must be supported for the case where
the operand is a $code const$$ $italic Base$$ object.

$head x$$
The operand $italic x$$ has one of the following prototypes
$syntax%
	const AD<%Base%>               &%x%
	const VecAD<%Base%>::reference &%x%
%$$

$head y$$
The result $italic y$$ has type
$syntax%
	AD<%Base%> %y%
%$$
It is equal to the negative of the operand $italic x$$.

$head Operation Sequence$$
This is an AD of $italic Base$$
$xref/glossary/Operation/Atomic/atomic operation/1/$$
and hence is part of the current
AD of $italic Base$$
$xref/glossary/Operation/Sequence/operation sequence/1/$$.

$head Derivative$$
If $latex f$$ is a 
$xref/glossary/Base Function/Base function/$$,
$latex \[
	\D{[ - f(x) ]}{x} = - \D{f(x)}{x}
\] $$

$head Example$$
$children%
	example/unary_minus.cpp
%$$
The file
$xref/UnaryMinus.cpp/$$
contains an example and test of this operation.

$end
-------------------------------------------------------------------------------
*/

//  BEGIN CppAD namespace
namespace CppAD {

// Broken g++ compiler inhibits declaring unary minus a member or friend
template <class Base>
inline AD<Base> AD<Base>::operator - (void) const 
{	// should make a more efficient version by adding unary minus to
	// Operator.h (some day)
	AD<Base> result(0);

	result  -= *this;

	return result;
}


template <class Base>
inline AD<Base> operator - (const VecAD_reference<Base> &right) 
{	return - right.ADBase(); }

}
//  END CppAD namespace


# endif

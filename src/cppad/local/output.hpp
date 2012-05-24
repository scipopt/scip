/* $Id: output.hpp 2091 2011-09-02 19:47:24Z bradbell $ */
# ifndef CPPAD_OUTPUT_INCLUDED
# define CPPAD_OUTPUT_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin Output$$
$spell
	VecAD
	std
	ostream
	const
$$

$index <<, AD output$$
$index AD, stream output$$
$index output, AD$$
$index stream, AD output$$
$index write, AD$$

$section AD Output Stream Operator$$ 

$head Syntax$$
$syntax%%os% << %x%$$


$head Purpose$$
Writes the $italic Base$$ value_, corresponding to $italic x$$,
to the output stream $italic os$$.

$head os$$
The operand $italic os$$ has prototype
$syntax%
	std::ostream &%os%
%$$

$head x$$
The operand $italic x$$ has one of the following prototypes
$syntax%
	const AD<%Base%>               &%x%
	const VecAD<%Base%>::reference &%x%
%$$

$head Result$$
The result of this operation can be used as a reference to $italic os$$.
For example, if the operand $italic y$$ has prototype
$syntax%
	AD<%Base%> %y%
%$$
then the syntax
$syntax%
	%os% << %x% << %y%
%$$
will output the value corresponding to $italic x$$
followed by the value corresponding to $italic y$$. 

$head Operation Sequence$$
The result of this operation is not an
$xref/glossary/AD of Base/AD of Base/$$ object.
Thus it will not be recorded as part of an
AD of $italic Base$$
$xref/glossary/Operation/Sequence/operation sequence/1/$$.

$head Example$$
$children%
	example/output.cpp
%$$
The file
$xref/Output.cpp/$$
contains an example and test of this operation.
It returns true if it succeeds and false otherwise.

$end
------------------------------------------------------------------------------
*/

namespace CppAD { 

	template <class Base>
	CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
	std::ostream& operator << 
	(std::ostream &os, const AD<Base> &x)
	{ 	return (os << x.value_); }

	template <class Base>
	CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
	std::ostream& operator << 
	(std::ostream &os, const VecAD_reference<Base> &x)
	{ 	return (os << x.ADBase()); }

}

# endif

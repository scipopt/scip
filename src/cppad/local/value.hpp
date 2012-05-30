/* $Id: value.hpp 2091 2011-09-02 19:47:24Z bradbell $ */
# ifndef CPPAD_VALUE_INCLUDED
# define CPPAD_VALUE_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin Value$$
$spell
	const
$$

$index Value$$

$index Base, from AD$$
$index AD, convert to Base$$
$index convert, AD to Base$$

$section Convert From an AD Type to its Base Type$$

$head Syntax$$
$syntax%%b% = Value(%x%)%$$


$head Purpose$$
Converts from an AD type to the corresponding
$xref/glossary/Base Type/base type/$$.

$head x$$
The argument $italic x$$ has prototype
$syntax%
	const AD<%Base%> &%x%
%$$

$head b$$
The return value $italic b$$ has prototype
$syntax%
	%Base% %b%
%$$

$head Operation Sequence$$
The result of this operation is not an 
$xref/glossary/AD of Base/AD of Base/$$ object.
Thus it will not be recorded as part of an
AD of $italic Base$$
$xref/glossary/Operation/Sequence/operation sequence/1/$$.

$head Restriction$$
If the argument $italic x$$ is a 
$xref/glossary/Variable/variable/$$ its dependency information 
would not be included in the $code Value$$ result (see above).
For this reason,
the argument $italic x$$ must be a $xref/glossary/Parameter/parameter/$$; i.e.,
it cannot depend on the current
$cref/independent variables/glossary/Tape/Independent Variable/$$.

$head Example$$
$children%
	example/value.cpp
%$$
The file
$xref/Value.cpp/$$
contains an example and test of this operation.

$end
-------------------------------------------------------------------------------
*/

//  BEGIN CppAD namespace
namespace CppAD {

template <class Base>
CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
Base Value(const AD<Base> &x) 
{	Base result;

	CPPAD_ASSERT_KNOWN(
		Parameter(x) ,
		"Value: argument is a variable (not a parameter)" 
	);
		

	result = x.value_;

	return result;
}

}
//  END CppAD namespace


# endif

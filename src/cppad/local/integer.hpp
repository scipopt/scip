/* $Id: integer.hpp 2091 2011-09-02 19:47:24Z bradbell $ */
# ifndef CPPAD_INTEGER_INCLUDED
# define CPPAD_INTEGER_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
------------------------------------------------------------------------------
$begin Integer$$
$spell
	std
	VecAD
	CppAD
	namespace
	const
	bool
$$

$index Integer$$

$index AD, convert to integer$$
$index convert, AD to integer$$

$section Convert From AD to Integer$$

$head Syntax$$
$syntax%%i% = Integer(%x%)%$$


$head Purpose$$
Converts from an AD type to the corresponding integer value.

$head i$$
The result $italic i$$ has prototype
$syntax%
	int %i%
%$$

$head x$$

$subhead Real Types$$
If the argument $italic x$$ has either of the following prototypes:
$syntax%
	const float                %%  &%x%
	const double               %%  &%x%
%$$
the fractional part is dropped to form the integer value.
For example, if $italic x$$ is 1.5, $italic i$$ is 1.
In general, if $latex x \geq 0$$, $italic i$$ is the 
greatest integer less than or equal $italic x$$.
If $latex x \leq 0$$, $italic i$$ is the 
smallest integer greater than or equal $italic x$$.

$subhead Complex Types$$
If the argument $italic x$$ has either of the following prototypes:
$syntax%
	const std::complex<float>  %%  &%x%
	const std::complex<double> %%  &%x%
%$$
The result $italic i$$ is given by
$syntax%
	%i% = Integer(%x%.real())
%$$

$subhead AD Types$$
If the argument $italic x$$ has either of the following prototypes:
$syntax%
	const AD<%Base%>               &%x%
	const VecAD<%Base%>::reference &%x%
%$$
$italic Base$$ must support the $code Integer$$ function and
the conversion has the same meaning as for $italic Base$$.

$head Operation Sequence$$
The result of this operation is not an
$xref/glossary/AD of Base/AD of Base/$$ object.
Thus it will not be recorded as part of an
AD of $italic Base$$
$xref/glossary/Operation/Sequence/operation sequence/1/$$.

$head Example$$
$children%
	example/integer.cpp
%$$
The file
$xref/Integer.cpp/$$
contains an example and test of this operation.

$end 
------------------------------------------------------------------------------
*/


namespace CppAD {

	template <class Base>
	CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
	int Integer(const AD<Base> &x)
	{	return Integer(x.value_); }

	template <class Base>
	CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
	int Integer(const VecAD_reference<Base> &x)
	{	return Integer( x.ADBase() ); }
}
# endif


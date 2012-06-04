/* $Id: near_equal_ext.hpp 2091 2011-09-02 19:47:24Z bradbell $ */
# ifndef CPPAD_NEAR_EQUAL_EXT_INCLUDED
# define CPPAD_NEAR_EQUAL_EXT_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
/*
$begin NearEqualExt$$
$spell
	cout
	endl
	Microsoft
	std
	Cpp
	namespace
	const
	bool
$$

$section Compare AD and Base Objects for Nearly Equal$$

$index NearEqual, AD with Base$$

$head Syntax$$
$syntax%%b% = NearEqual(%x%, %y%, %r%, %a%)%$$


$head Purpose$$
The routine $xref/NearEqual/$$ determines if two objects of
the same type are nearly.
This routine is extended to the case where one object can have type
$italic Type$$ while the other can have type
$syntax%AD<%Type%>%$$ or
$syntax%AD< std::complex<%Type%> >%$$.

$head x$$
The arguments $italic x$$ 
has one of the following possible prototypes:
$syntax%
	const %Type%                     &%x%
	const AD<%Type%>                 &%x%
	const AD< std::complex<%Type%> > &%x%
%$$

$head y$$
The arguments $italic y$$ 
has one of the following possible prototypes:
$syntax%
	const %Type%                     &%y%
	const AD<%Type%>                 &%y%
	const AD< std::complex<%Type%> > &%x%
%$$


$head r$$
The relative error criteria $italic r$$ has prototype
$syntax%
	const %Type% &%r%
%$$
It must be greater than or equal to zero.
The relative error condition is defined as:
$latex \[
	\frac{ | x - y | } { |x| + |y| } \leq r
\] $$

$head a$$
The absolute error criteria $italic a$$ has prototype
$syntax%
	const %Type% &%a%
%$$
It must be greater than or equal to zero.
The absolute error condition is defined as:
$latex \[
	| x - y | \leq a
\] $$

$head b$$
The return value $italic b$$ has prototype
$syntax%
	bool %b%
%$$
If either $italic x$$ or $italic y$$ is infinite or not a number, 
the return value is false.
Otherwise, if either the relative or absolute error 
condition (defined above) is satisfied, the return value is true.
Otherwise, the return value is false.

$head Type$$
The type $italic Type$$ must be a
$xref/NumericType/$$.
The routine $xref/CheckNumericType/$$ will generate
an error message if this is not the case.
If $italic a$$ and $italic b$$ have type $italic Type$$,
the following operation must be defined 
$table
$bold Operation$$     $cnext 
	$bold Description$$ $rnext
$syntax%%a% <= %b%$$  $cnext 
	less that or equal operator (returns a $code bool$$ object)
$tend

$head Operation Sequence$$
The result of this operation is not an
$xref/glossary/AD of Base/AD of Base/$$ object.
Thus it will not be recorded as part of an
AD of $italic Base$$
$xref/glossary/Operation/Sequence/operation sequence/1/$$.

$head Example$$
$children%
	example/near_equal_ext.cpp
%$$
The file $xref/NearEqualExt.cpp/$$ contains an example
and test of this extension of $xref/NearEqual/$$.
It return true if it succeeds and false otherwise.

$end

*/
// BEGIN CppAD namespace
namespace CppAD {
// ------------------------------------------------------------------------

// fold into base type and then use <cppad/near_equal.hpp>
template <class Base>
CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
bool NearEqual(
const AD<Base> &x, const AD<Base> &y, const Base &r, const Base &a)
{	return NearEqual(x.value_, y.value_, r, a);
}

template <class Base>
CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
bool NearEqual(
const Base &x, const AD<Base> &y, const Base &r, const Base &a)
{	return NearEqual(x, y.value_, r, a);
}

template <class Base>
CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
bool NearEqual(
const AD<Base> &x, const Base &y, const Base &r, const Base &a)
{	return NearEqual(x.value_, y, r, a);
}

// fold into AD type and then use cases above
template <class Base>
CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
bool NearEqual(
	const VecAD_reference<Base> &x, const VecAD_reference<Base> &y, 
	const Base &r, const Base &a)
{	return NearEqual(x.ADBase(), y.ADBase(), r, a);
}
template <class Base>
CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
bool NearEqual(const VecAD_reference<Base> &x, const AD<Base> &y, 
	const Base &r, const Base &a)
{	return NearEqual(x.ADBase(), y, r, a);
}
template <class Base>
CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
bool NearEqual(const VecAD_reference<Base> &x, const Base &y, 
	const Base &r, const Base &a)
{	return NearEqual(x.ADBase(), y, r, a);
}
template <class Base>
CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
bool NearEqual(const AD<Base> &x, const VecAD_reference<Base> &y, 
	const Base &r, const Base &a)
{	return NearEqual(x, y.ADBase(), r, a);
}
template <class Base>
CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
bool NearEqual(const Base &x, const VecAD_reference<Base> &y, 
	const Base &r, const Base &a)
{	return NearEqual(x, y.ADBase(), r, a);
}

} // END CppAD namespace

# endif

/* $Id: ad_copy.hpp 1886 2011-02-20 03:17:03Z bradbell $ */
# ifndef CPPAD_AD_COPY_INCLUDED
# define CPPAD_AD_COPY_INCLUDED

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

$begin ad_copy$$
$spell
	Vec
	const
$$

$index AD, copy$$
$index AD, assignment$$

$index constructor, AD$$
$index assignment, AD$$

$index convert, to AD$$
$index Base, convert to AD$$
$index double, convert to AD$$
$index VecAD, convert to AD$$

$section AD Copy Constructor and Assignment Operator$$

$head Syntax$$

$subhead Constructor$$
$syntax%AD<%Base%> %y%(%x%)
%$$
$syntax%AD<%Base%> %y% = %x%$$

$subhead Assignment$$
$syntax%%y% = %x%$$

$head Purpose$$
The constructor creates a new $syntax%AD<%Base%>%$$ object $italic y$$ 
and the assignment operator uses an existing $italic y$$.
In either case, 
$italic y$$ has the same value as $italic x$$,
and the same dependence on the 
$cref/independent variables/glossary/Tape/Independent Variable/$$
($italic y$$ is a 
$xref/glossary/Variable/variable/$$ if and only if $italic x$$ is a variable).

$head x$$
The argument $italic x$$ has prototype
$syntax%
	const %Type% &%x%
%$$ 
where $italic Type$$ is
$syntax%VecAD<%Base%>::reference%$$,
$syntax%AD<%Base%>%$$,
$italic Base$$, or
$code double$$.

$head y$$
The target $italic y$$ has prototype
$syntax%
	AD<%Base%> &%y%
%$$

$head Example$$
$children%
	example/copy_ad.cpp%
	example/copy_base.cpp%
	example/eq.cpp
	
%$$
The following files contain examples and tests of these operations.
Each test returns true if it succeeds and false otherwise.
$table
$rref CopyAD.cpp$$
$rref CopyBase.cpp$$
$rref Eq.cpp$$
$tend

$end
------------------------------------------------------------------------------
*/

//  BEGIN CppAD namespace
namespace CppAD {

/*
// use default copy constructor and assignment operator
// because they may be optimized better by the compiler
template <class Base>
inline AD<Base>::AD(const AD &x) 
{
	value_   = x.value_;
	id_     = x.id_;
	taddr_  = x.taddr_;

	return;
}
template <class Base>
inline AD<Base>& AD<Base>::operator=(const AD<Base> &right)
{	value_   = right.value_;
	id_     = right.id_;
	taddr_  = right.taddr_;

	// check that all variables are parameters while tape is empty
	CPPAD_ASSERT_UNKNOWN(
		Parameter(*this) || (tape_this() != CPPAD_NULL)
	);

	return *this;
}
*/

// Constructor and assignment from Base type
// Initilaize id_ to CPPAD_MAX_NUM_THREADS, so that following conditions hold 
// id_ != 0 , id_ % CPPAD_MAX_NUM_THREADS == 0, id_ != any recording tape id.
// taddr_ is not used, set anyway to avoid compile warning.
template <class Base>
inline AD<Base>::AD(const Base &b) 
: value_(b)
, id_(CPPAD_MAX_NUM_THREADS)
, taddr_(0)
{ }	
template <class Base>
inline AD<Base>& AD<Base>::operator=(const Base &b)
{	value_ = b;
	id_    = CPPAD_MAX_NUM_THREADS;

	// check that this is a parameter
	CPPAD_ASSERT_UNKNOWN( Parameter(*this) );

	return *this;
}

// constructor and assignment from VecAD_reference<Base>
template <class Base>
inline AD<Base>::AD(const VecAD_reference<Base> &x)
{	*this = x.ADBase(); }
template <class Base>
inline AD<Base>& AD<Base>::operator=(const VecAD_reference<Base> &x)
{	return *this = x.ADBase(); }

// constructor and assignment from any other type 
// taddr_ is not used, set anyway to avoid compile warning.
template <class Base>
template <class T>
inline AD<Base>::AD(const T &t) 
: value_(Base(t))
, id_(CPPAD_MAX_NUM_THREADS)
, taddr_(0)
{ }
template <class Base>
template <class T>
inline AD<Base>& AD<Base>::operator=(const T &t)
{	return *this = Base(t); }


} // END CppAD namespace

# endif

/* $Id: default.hpp 1886 2011-02-20 03:17:03Z bradbell $ */
# ifndef CPPAD_DEFAULT_INCLUDED
# define CPPAD_DEFAULT_INCLUDED

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

$begin Default$$
$spell
$$

$section AD Default Constructor$$

$index AD, default construct$$
$index construct, default$$
$index default, constructor$$

$head Syntax$$
$syntax%AD<%Base%> %x%;%$$


$head Purpose$$
Constructs an AD object with an unspecified value.
Directly after this construction,
the object is a $xref/glossary/Parameter/parameter/$$.


$head Example$$
$children%
	example/default.cpp
%$$
The file
$xref/Default.cpp/$$
contains an example and test of this operation.
It returns true if it succeeds and false otherwise.

$end
------------------------------------------------------------------------------
*/
//  BEGIN CppAD namespace
namespace CppAD {

// default constructor
//
// Initilaize id_ to CPPAD_MAX_NUM_THREADS, so that following conditions hold 
// id_ != 0 , id_ % CPPAD_MAX_NUM_THREADS == 0, id_ != any recording tape id.
// taddr_ is not used, set anyway to avoid compile warning.
// value_ = 0 so complex of AD will work, see comments in not_complex_ad.cpp.
template <class Base>
inline AD<Base>::AD(void) 
: value_(Base(0))
, id_(CPPAD_MAX_NUM_THREADS)
, taddr_(0)
{ }	

} // END CppAD namespace

# endif

/* $Id: epsilon.hpp 2057 2011-08-11 14:07:11Z bradbell $ */
# ifndef CPPAD_EPSILON_INCLUDED
# define CPPAD_EPSILON_INCLUDED
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
$begin epsilon$$
$spell
	eps
	CppAD
	namespace
	const
$$
$index AD, epsilon$$
$index epsilon, AD$$

$section Machine Epsilon For AD Types$$

$head Syntax$$
$icode%eps% = epsilon< AD<%Base%> >()%$$

$head Purpose$$
Obtain the value of machine epsilon corresponding
to the type $codei%AD<%Base%>%$$.

$head eps$$
The result $icode eps$$ has prototype
$codei%
	AD<%Base%> eps
%$$

$head Example$$
$children%
	example/epsilon.cpp
%$$
The file
$cref/epsilon.cpp/$$
contains an example and test of this operation.

$end 
------------------------------------------------------------------------------
*/

namespace CppAD {

	template <class ADtype>
	inline ADtype epsilon(void)
	{	typedef typename ADtype::value_type value_type;
		return ADtype ( epsilon<value_type>() ); 
	}

}
# endif

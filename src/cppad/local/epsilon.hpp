/* $Id: epsilon.hpp 2506 2012-10-24 19:36:49Z bradbell $ */
# ifndef CPPAD_EPSILON_INCLUDED
# define CPPAD_EPSILON_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-12 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Eclipse Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
------------------------------------------------------------------------------
$begin epsilon$$
$spell
	std
	eps
	CppAD
	namespace
	const
$$
$index AD, epsilon$$
$index epsilon, AD$$

$section Machine Epsilon For AD Types$$

$head Deprecated$$
$index deprecated, epsilon$$
This routine has bee deprecated.
You should use the $cref/numeric_limits/limits/$$ $code epsilon$$ instead.

$head Syntax$$
$icode%eps% = epsilon<%Float%>()%$$

$head Purpose$$
Obtain the value of machine epsilon corresponding
to the type $icode%Float%$$.

$head Float$$
this type can either be $codei%AD<%Base%>%$$,
or it can be $icode Base$$ for any $codei%AD<%Base%>%$$ type.

$head eps$$
The result $icode eps$$ has prototype
$codei%
	%Float% eps
%$$

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

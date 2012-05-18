/* $Id: var2par.hpp 2091 2011-09-02 19:47:24Z bradbell $ */
# ifndef CPPAD_VAR2PAR_INCLUDED
# define CPPAD_VAR2PAR_INCLUDED

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

$begin Var2Par$$
$spell
	var
	const
$$

$index Var2Par$$
$index variable, convert to parameter$$
$index parameter, convert from variable$$
$index value_, obtain during taping$$

$section Convert an AD Variable to a Parameter$$

$head Syntax$$
$syntax%%y% = Var2Par(%x%)%$$


$head Purpose$$
Returns a 
$cref/parameter/glossary/Parameter/$$ $italic y$$
with the same value as the 
$cref/variable/glossary/Variable/$$ $italic x$$.

$head x$$
The argument $italic x$$ has prototype
$syntax%
	const AD<%Base%> &x
%$$ 
The argument $italic x$$ may be a variable or parameter.


$head y$$
The result $italic y$$ has prototype
$syntax%
	AD<%Base%> &y
%$$ 
The return value $italic y$$ will be a parameter. 


$head Example$$
$children%
	example/var_2par.cpp
%$$
The file
$xref/Var2Par.cpp/$$
contains an example and test of this operation.
It returns true if it succeeds and false otherwise.

$end
------------------------------------------------------------------------------
*/

//  BEGIN CppAD namespace
namespace CppAD {

template <class Base>
CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
AD<Base> Var2Par(const AD<Base> &x)
{	AD<Base> y(x.value_);
	return y;
}


template <class Base>
CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
AD<Base> Var2Par(const VecAD_reference<Base> &x)
{	AD<Base> y(x.ADBase()); 
	y.id_ = 0;
}


} // END CppAD namespace

# endif

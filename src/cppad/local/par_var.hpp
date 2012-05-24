/* $Id: par_var.hpp 2091 2011-09-02 19:47:24Z bradbell $ */
# ifndef CPPAD_PAR_VAR_INCLUDED
# define CPPAD_PAR_VAR_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
/*
---------------------------------------------------------------------------

$begin ParVar$$
$spell
	VecAD
	const
	bool
$$

$index Parameter$$
$index Variable$$
$section Is an AD Object a Parameter or Variable$$

$head Syntax$$
$syntax%%b% = Parameter(%x%)%$$
$pre
$$
$syntax%%b% = Variable(%x%)%$$


$head Purpose$$
Determine if $italic x$$ is a 
$xref/glossary/Parameter/parameter/$$ or 
$xref/glossary/Variable/variable/$$. 

$head x$$
The argument $italic x$$ has prototype
$syntax%
	const AD<%Base%>    &%x%
	const VecAD<%Base%> &%x%
%$$

$head b$$
The return value $italic b$$ has prototype
$syntax%
	bool %b%
%$$
The return value for $code Parameter$$ ($code Variable$$)
is true if and only if $italic x$$ is a parameter (variable).
Note that a $cref/VecAD<Base>/VecAD/$$ object
is a variable if any element of the vector depends on the independent
variables.

$head Operation Sequence$$
The result of this operation is not an
$xref/glossary/AD of Base/AD of Base/$$ object.
Thus it will not be recorded as part of an
AD of $italic Base$$
$xref/glossary/Operation/Sequence/operation sequence/1/$$.

$head Example$$
$children%
	example/par_var.cpp
%$$
The file
$xref/ParVar.cpp/$$
contains an example and test of these functions.
It returns true if it succeeds and false otherwise.

$end
-----------------------------------------------------------------------------
*/

namespace CppAD {
	// Parameter
	template <class Base>
	CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
	bool Parameter(const AD<Base> &x)
	{	size_t thread = x.id_ % CPPAD_MAX_NUM_THREADS;
		return x.id_ != * AD<Base>::id_handle(thread); 
	}

	template <class Base>
	CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
	bool Parameter(const VecAD<Base> &x)
	{	size_t thread = x.id_ % CPPAD_MAX_NUM_THREADS;
		return x.id_ != * AD<Base>::id_handle(thread); 
	}

	// Variable
	template <class Base>
	CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
	bool Variable(const AD<Base> &x)
	{	size_t thread = x.id_ % CPPAD_MAX_NUM_THREADS;
		return x.id_ == * AD<Base>::id_handle(thread); 
	}

	template <class Base>
	CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
	bool Variable(const VecAD<Base> &x)
	{	size_t thread = x.id_ % CPPAD_MAX_NUM_THREADS;
		return x.id_ == * AD<Base>::id_handle(thread); 
	}
} 
// END CppAD namespace


# endif

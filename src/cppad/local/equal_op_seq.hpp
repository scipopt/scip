/* $Id: equal_op_seq.hpp 2091 2011-09-02 19:47:24Z bradbell $ */
# ifndef CPPAD_EQUAL_OP_SEQ_INCLUDED
# define CPPAD_EQUAL_OP_SEQ_INCLUDED

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
$begin EqualOpSeq$$ 
$spell
	Op
	const
	bool
$$

$index EqualOpSeq$$
$index operation, equal sequence$$
$index sequence, equal operation$$
$index equal, operation sequence$$

$section Check if Two Value are Identically Equal$$

$head Syntax$$
$syntax%%b% = EqualOpSeq(%x%, %y%)%$$

$head Purpose$$
Determine if two $italic x$$ and $italic y$$ are identically equal; i.e.,
not only is $icode%x% == %y%$$ true, but
if they are $cref/variables/glossary/Variable/$$,
they correspond have the same
$cref/operation sequence/glossary/Operation/Sequence/$$.

$head Motivation$$
Sometimes it is useful to cache information
and only recalculate when a function's arguments change.
In the case of AD variables, 
it may be important not only when the argument values are equal,
but when they are related to the
$cref/independent variables/glossary/Tape/Independent Variable/$$ 
by the same operation sequence.
After the assignment
$syntax%
	%y% = %x%
%$$
these two AD objects would not only have equal values,
but would also correspond to the same operation sequence.

$head x$$
The argument $italic x$$ has prototype
$syntax%
	const AD<%Base%> &%x%
%$$

$head y$$
The argument $italic y$$ has prototype
$syntax%
	const AD<%Base%> &%y%
%$$

$head b$$
The result $italic b$$ has prototype
$syntax%
	bool %b%
%$$
The result is true if and only if one of the following cases holds:

$list number$$
Both $italic x$$ and $italic y$$ are variables 
and correspond to the same operation sequence.
$lnext
Both $italic x$$ and $italic y$$ are parameters,
$italic Base$$ is an AD type,
and $syntax%EqualOpSeq( Value(%x%) , Value(%y%) )%$$ is true.
$lnext
Both $italic x$$ and $italic y$$ are parameters,
$italic Base$$ is not an AD type,
and $syntax%%x% == %y%%$$ is true.
$lend


$head Example$$
$children%
	example/equal_op_seq.cpp
%$$
The file
$xref/EqualOpSeq.cpp/$$
contains an example and test of $code EqualOpSeq$$.
It returns true if it succeeds and false otherwise.


$end 
------------------------------------------------------------------------------
*/


namespace CppAD {
	template <class Base>
	CPPAD_INLINE_FRIEND_TEMPLATE_FUNCTION
	bool EqualOpSeq(const AD<Base> &x, const AD<Base> &y)
	{	
		if( Parameter(x) )
		{	if( Parameter(y) )
				return EqualOpSeq(x.value_, y.value_);
			else	return false;
		}
		else if( Parameter(y) )
			return false;

		return (x.taddr_ == y.taddr_);
	}
		
}

# endif

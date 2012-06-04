/* $Id: compute_assign.hpp 1369 2009-05-31 01:31:48Z bradbell $ */
# ifndef CPPAD_COMPUTE_ASSIGN_INCLUDED
# define CPPAD_COMPUTE_ASSIGN_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-07 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
-------------------------------------------------------------------------------
$begin compute_assign$$
$spell
	Op
	VecAD
	const
$$

$section AD Computed Assignment Operators$$

$index assignment, operator$$
$index operator, assignment$$

$index +, computed assignment$$
$index add, computed assignment$$
$index plus, computed assignment$$

$index -, computed assignment$$
$index subtract, computed assignment$$
$index minus, computed assignment$$

$index *, computed assignment$$
$index multiply, computed assignment$$
$index times, computed assignment$$

$index /, computed assignment$$
$index divide, computed assignment$$

$head Syntax$$
$syntax%%x% %Op% %y%$$

$head Purpose$$
Performs computed assignment operations 
where either $italic x$$ has type 
$syntax%AD<%Base%>%$$.

$head Op$$
The operator $italic Op$$ is one of the following
$table
$bold Op$$  $cnext $bold Meaning$$ $rnext
$code +=$$   $cnext $italic x$$ is assigned $italic x$$ plus $italic y$$ $rnext
$code -=$$   $cnext $italic x$$ is assigned $italic x$$ minus $italic y$$ $rnext
$code *=$$   $cnext $italic x$$ is assigned $italic x$$ times $italic y$$ $rnext
$code /=$$   $cnext $italic x$$ is assigned $italic x$$ divided by $italic y$$ 
$tend

$head Base$$
The type $italic Base$$ is determined by the operand $italic x$$.

$head x$$
The operand $italic x$$ has the following prototype
$syntax%
	AD<%Base%> &%x%
%$$

$head y$$
The operand $italic y$$ has the following prototype
$syntax%
	const %Type% &%y%
%$$
where $italic Type$$ is
$syntax%VecAD<%Base%>::reference%$$,
$syntax%AD<%Base%>%$$,
$italic Base$$, or
$code double$$.

$head Result$$
$index assignment, multiple$$
$index multiple, assignment$$
The result of this assignment 
can be used as a reference to $italic x$$.
For example, if $italic z$$ has the following type
$syntax%
	AD<%Base%> %z%
%$$
then the syntax
$syntax%
	%z% = %x% += %y%
%$$
will compute $italic x$$ plus $italic y$$
and then assign this value to both $italic x$$ and $italic z$$. 


$head Operation Sequence$$
This is an $cref/atomic/glossary/Operation/Atomic/$$
$cref/AD of Base/glossary/AD of Base/$$ operation
and hence it is part of the current
AD of $italic Base$$ 
$cref/operation sequence/glossary/Operation/Sequence/$$.

$children%
	example/add_eq.cpp%
	example/sub_eq.cpp%
	example/mul_eq.cpp%
	example/div_eq.cpp
%$$

$head Example$$
The following files contain examples and tests of these functions.
Each test returns true if it succeeds and false otherwise.
$table
$rref AddEq.cpp$$
$rref SubEq.cpp$$
$rref MulEq.cpp$$
$rref DivEq.cpp$$
$tend

$head Derivative$$
If $latex f$$ and $latex g$$ are 
$xref/glossary/Base Function/Base functions/$$

$subhead Addition$$
$latex \[
	\D{[ f(x) + g(x) ]}{x} = \D{f(x)}{x} + \D{g(x)}{x}
\] $$

$subhead Subtraction$$
$latex \[
	\D{[ f(x) - g(x) ]}{x} = \D{f(x)}{x} - \D{g(x)}{x}
\] $$

$subhead Multiplication$$
$latex \[
	\D{[ f(x) * g(x) ]}{x} = g(x) * \D{f(x)}{x} + f(x) * \D{g(x)}{x}
\] $$

$subhead Division$$
$latex \[
	\D{[ f(x) / g(x) ]}{x} = 
		[1/g(x)] * \D{f(x)}{x} - [f(x)/g(x)^2] * \D{g(x)}{x}
\] $$

$end 
-----------------------------------------------------------------------------
*/
# include <cppad/local/add_eq.hpp>
# include <cppad/local/sub_eq.hpp>
# include <cppad/local/mul_eq.hpp>
# include <cppad/local/div_eq.hpp>

# endif

/* $Id: ad_binary.hpp 1369 2009-05-31 01:31:48Z bradbell $ */
# ifndef CPPAD_AD_BINARY_INCLUDED
# define CPPAD_AD_BINARY_INCLUDED

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
$begin ad_binary$$
$spell
	Op
	VecAD
	const
$$

$section AD Binary Arithmetic Operators$$

$index binary, operator$$
$index operator, binary$$

$index +, binary operator$$
$index add, binary operator$$
$index plus, binary operator$$

$index -, binary operator$$
$index subtract, binary operator$$
$index minus, binary operator$$

$index *, binary operator$$
$index multiply, binary operator$$
$index times, binary operator$$

$index /, binary operator$$
$index divide, binary operator$$


$head Syntax$$
$syntax%%z% = %x% %Op% %y%$$

$head Purpose$$
Performs arithmetic operations where either $italic x$$ or $italic y$$
has type 
$syntax%AD<%Base%>%$$ or
$cref%VecAD<Base>::reference%VecAD%VecAD<Base>::reference%$$. 

$head Op$$
The operator $italic Op$$ is one of the following
$table
$bold Op$$  $cnext $bold Meaning$$ $rnext
$code +$$   $cnext $italic z$$ is $italic x$$ plus $italic y$$ $rnext
$code -$$   $cnext $italic z$$ is $italic x$$ minus $italic y$$ $rnext
$code *$$   $cnext $italic z$$ is $italic x$$ times $italic y$$ $rnext
$code /$$   $cnext $italic z$$ is $italic x$$ divided by $italic y$$ 
$tend

$head Base$$
The type $italic Base$$ is determined by the operand that
has type $syntax%AD<%Base%>%$$ or $syntax%VecAD<%Base%>::reference%$$.

$head x$$
The operand $italic x$$ has the following prototype
$syntax%
	const %Type% &%x%
%$$
where $italic Type$$ is
$syntax%VecAD<%Base%>::reference%$$,
$syntax%AD<%Base%>%$$,
$italic Base$$, or
$code double$$.

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


$head z$$
The result $italic z$$ has the following prototype
$syntax%
	%Type% %z%
%$$
where $italic Type$$ is
$syntax%AD<%Base%>%$$.

$head Operation Sequence$$
This is an $cref/atomic/glossary/Operation/Atomic/$$
$cref/AD of Base/glossary/AD of Base/$$ operation
and hence it is part of the current
AD of $italic Base$$ 
$cref/operation sequence/glossary/Operation/Sequence/$$.

$children%
	example/add.cpp%
	example/sub.cpp%
	example/mul.cpp%
	example/div.cpp
%$$

$head Example$$
The following files contain examples and tests of these functions.
Each test returns true if it succeeds and false otherwise.
$table
$rref Add.cpp$$
$rref Sub.cpp$$
$rref Mul.cpp$$
$rref Div.cpp$$
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
# include <cppad/local/add.hpp>
# include <cppad/local/sub.hpp>
# include <cppad/local/mul.hpp>
# include <cppad/local/div.hpp>

# endif

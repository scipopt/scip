/* $Id: base_std_math.hpp 2240 2011-12-31 05:33:55Z bradbell $ */
# ifndef CPPAD_BASE_STD_MATH_INCLUDED
# define CPPAD_BASE_STD_MATH_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/* 
$begin base_std_math$$
$spell
	std
	acos
	asin
	atan
	cos
	exp
	sqrt
	const
	CppAD
	namespace
$$

$section Base Type Requirements for Standard Math Functions$$
$index math, base require$$
$index base, math require$$
$index require, base math$$

$head Purpose$$
These definitions are required for the user's code to use the type
$codei%AD<%Base%>%$$:

$head Unary Standard Math$$ 
$index math, base unary$$
$index base, unary math$$
$index unary, base math$$
The type $icode Base$$ must support the following functions
unary standard math functions (in the CppAD namespace):
$table
$bold Syntax$$ $cnext $bold Result$$ 
$rnext
$icode%y% = abs(%x%)%$$  $cnext absolute value     $rnext
$icode%y% = acos(%x%)%$$ $cnext inverse cosine     $rnext
$icode%y% = asin(%x%)%$$ $cnext inverse sine       $rnext
$icode%y% = atan(%x%)%$$ $cnext inverse tangent    $rnext
$icode%y% = cos(%x%)%$$  $cnext cosine             $rnext
$icode%y% = cosh(%x%)%$$ $cnext hyperbolic cosine  $rnext
$icode%y% = exp(%x%)%$$  $cnext exponential        $rnext
$icode%y% = log(%x%)%$$  $cnext natural logarithm  $rnext
$icode%y% = sin(%x%)%$$  $cnext sine               $rnext
$icode%y% = sinh(%x%)%$$ $cnext hyperbolic sine    $rnext
$icode%y% = sqrt(%x%)%$$ $cnext square root        $rnext
$icode%y% = tan(%x%)%$$  $cnext tangent           
$tend
where the arguments and return value have the prototypes
$codei%
	const %Base%& %x%
	%Base%        %y%
%$$


$head CPPAD_STANDARD_MATH_UNARY$$
$index CPPAD_STANDARD_MATH_UNARY$$
The macro invocation, within the CppAD namespace,
$codei%
	CPPAD_STANDARD_MATH_UNARY(%Base%, %Fun%)
%$$
defines the syntax
$codei%
	%y% = CppAD::%Fun%(%x%)
%$$
This macro uses the functions $codei%std::%Fun%$$ which
must be defined and have the same prototype as $codei%CppAD::%Fun%$$.
For example, 
$cref/float/base_float.hpp/Unary Standard Math/$$,
$cref/double/base_double.hpp/Unary Standard Math/$$, and
$cref/complex valid unary math/base_complex.hpp/Valid Unary Math/$$.

$head sign$$
$index sign, base require$$
$index base, sign require$$
$index require, base sign$$
The type $icode Base$$ must support the syntax
$codei%
	%y% = CppAD::sign(%x%)
%$$
which computes
$latex \[
z = \left\{ \begin{array}{ll} 
	+1 & {\rm if} \; x > 0 \\
	 0 & {\rm if} \; x = 0 \\
	-1 & {\rm if} \; x < 0
\end{array} \right.
\] $$
The arguments $icode x$$ has prototype
$codei%
	const %Base%& %x%
%$$
and the return value $icode z$$ has prototype
$codei%
	%Base% %z%
%$$
For example, see
$cref/float/base_float.hpp/sign/$$,
$cref/double/base_double.hpp/sign/$$.
Note that, if ordered comparisons are not defined for the type $icode Base$$,
the $code code sign$$ function should generate an assert if it is used; see
$cref/complex invalid unary math/base_complex.hpp/Invalid Unary Math/$$.

$head pow$$
$index pow, base require$$
$index base, pow require$$
$index require, base pow$$
The type $icode Base$$ must support the syntax
$codei%
	%z% = CppAD::pow(%x%, %y%)
%$$
which computes $latex z = x^y$$.
The arguments $icode x$$ and $italic y$$ have prototypes
$codei%
	const %Base%& %x%
	const %Base%& %y%
%$$
and the return value $icode z$$ has prototype
$codei%
	%Base% %z%
%$$
For example, see
$cref/float/base_float.hpp/pow/$$,
$cref/double/base_double.hpp/pow/$$, and
$cref/complex/base_complex.hpp/pow/$$.

$end
-------------------------------------------------------------------------------
*/

# include <cmath>

CPPAD_BEGIN_NAMESPACE

/*!
\file base_std_math.hpp
Defintions that aid meeting Base type requirements for standard math functions.
*/

/*!
\def CPPAD_STANDARD_MATH_UNARY(Type, Fun)
This macro defines the function
\verbatim
	y = CppAD:Fun(x)
\endverbatim
where the argument \c x and return value \c y have type \c Type
using the corresponding function <code>std::Fun</code>.
*/
# define CPPAD_STANDARD_MATH_UNARY(Type, Fun) \
	inline Type Fun(const Type& x)           \
	{	return std::Fun(x); }

CPPAD_END_NAMESPACE

# endif

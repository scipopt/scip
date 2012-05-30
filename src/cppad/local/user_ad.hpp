/* $Id: user_ad.hpp 2057 2011-08-11 14:07:11Z bradbell $ */
# ifndef CPPAD_USER_AD_INCLUDED
# define CPPAD_USER_AD_INCLUDED

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

$begin AD$$
$spell
	std
	bool
	cos
	Cpp
$$

$section AD Objects$$

$index AD, object$$

$head Purpose$$
The sections listed below describe the operations 
that are available to $xref/glossary/AD of Base/AD of Base/$$ objects.
These objects are used to $cref/tape/glossary/Tape/$$
an AD of $italic Base$$
$xref/glossary/Operation/Sequence/operation sequence/1/$$.
This operation sequence can
be transferred to an $xref/ADFun/$$ object where it
can be used to evaluate the corresponding 
function and derivative values.

$head Base Type Requirements$$
$index Base, require$$
The $italic Base$$ requirements are provided by the CppAD package 
for the following base types:
$code float$$, 
$code double$$,
$code std::complex<float>$$, 
$code std::complex<double>$$.
Otherwise, see $cref/base_require/$$.


$childtable%
	cppad/local/default.hpp%
	cppad/local/ad_copy.hpp%
	cppad/local/convert.hpp%
	cppad/local/ad_valued.hpp%
	cppad/local/bool_valued.hpp%
	cppad/local/vec_ad.hpp%
	cppad/base_require.hpp
%$$

$end
---------------------------------------------------------------------------
*/

# include <cppad/local/default.hpp>
# include <cppad/local/ad_copy.hpp>
# include <cppad/local/convert.hpp>
# include <cppad/local/vec_ad.hpp>
# include <cppad/local/ad_valued.hpp>
# include <cppad/local/bool_valued.hpp>

# endif

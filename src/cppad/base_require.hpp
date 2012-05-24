// $Id: base_require.hpp 2057 2011-08-11 14:07:11Z bradbell $
# ifndef CPPAD_BASE_REQUIRE_INCLUDED
# define CPPAD_BASE_REQUIRE_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin base_require$$
$spell
	eps
	std
	isnan
	Lt
	Le
	Eq
	Ge
	Gt
	cppad.hpp
	namespace
	acos
	asin
	atan
	cos
	sqrt
	optimizations
	bool
	const
	CppAD
	enum
	Lt
	Le
	Eq
	Ge 
	Gt
	inline
	Op
	std
	CondExp
	erf
$$

$index Base, require$$
$index require, Base type$$
$index type, Base require$$

$section AD<Base> Requirements for Base Type$$

$head Syntax$$
$code include <cppad/base_require.hpp>$$

$head Warning$$
This is a preliminary version of these specifications
and it is subject to change in future versions of CppAD.

$head Purpose$$
This section lists the requirements that the type
$icode Base$$ so that the type $codei%AD<%Base%>%$$ can be used.

$subhead Standard Base Types$$
In the case where $icode Base$$ is 
$code float$$, 
$code double$$,
$code std::complex<float>$$, 
$code std::complex<double>$$,
or $codei%AD<%Other%>%$$,
these requirements are provided by including the file
$code cppad/cppad.hpp$$.

$head Include Order$$
If you are linking a non-standard base type to CppAD,
you must first include the file $code cppad/base_require.hpp/$$,
then provide the specifications below, 
and then include the file $code cppad/cppad.hpp$$.

$head Numeric Type$$
The type $icode Base$$ must support all the operations for a 
$cref/NumericType/$$.

$head Integer$$
$index Integer, base require$$
$index base, Integer require$$
$index require, base Integer$$
The type $icode Base$$ must support the syntax
$codei%
	%i% = CppAD::Integer(%x%)
%$$
which converts $icode x$$ to an $code int$$.
The argument $icode x$$ has prototype
$codei%
	const %Base%& %x%
%$$
and the return value $icode i$$ has prototype
$codei%
	int %i%
%$$

$subhead Suggestion$$
In many cases, the $icode Base$$ version of the $code Integer$$ function
can be defined by
$codei%
namespace CppAD {
	inline int Integer(const %Base%& x)
	{	return static_cast<int>(x); }
}
%$$
For example, see
$cref/float/base_float.hpp/Integer/$$,
$cref/double/base_double.hpp/Integer/$$, and
$cref/complex/base_complex.hpp/Integer/$$.

$head epsilon$$
$index epsilon, machine base$$
$index machine, epsilon base$$
$index base, machine epsilon$$
The template specialization
$codei%
	%Base% CppAD::epsilon<%Base%>()
%$$
must return machine epsilon for the type $icode Base$$.
For example, see
$cref/float/base_float.hpp/epsilon/$$,
$cref/float/base_double.hpp/epsilon/$$, and
$cref/float/base_complex.hpp/epsilon/$$,

$head isnan$$
If $icode Base$$ defines the $code isnan$$ function,
you may also have to provide a definition in the CppAD namespace
(to avoid a function ambiguity).
For example, 
$cref/complex isnan/base_complex.hpp/isnan/$$.

$childtable%
	cppad/local/base_cond_exp.hpp%
	omh/base_identical.omh%
	omh/base_ordered.omh%
	cppad/local/base_std_math.hpp%
	cppad/local/base_float.hpp%
	cppad/local/base_double.hpp%
	cppad/local/base_complex.hpp%
	example/base_adolc.hpp
%$$

$end
*/

// definitions that must come before base implementations
# include <cppad/error_handler.hpp>
# include <cppad/local/define.hpp>
# include <cppad/local/cppad_assert.hpp>

// grouping documentation by feature
# include <cppad/local/base_cond_exp.hpp>
# include <cppad/local/base_std_math.hpp>

// must define template AD version of epsilon before the base cases
# include <cppad/local/epsilon.hpp>

// base cases that come with CppAD
# include <cppad/local/base_float.hpp>
# include <cppad/local/base_double.hpp>
# include <cppad/local/base_complex.hpp>

# endif

/* $Id: time_test.hpp 2215 2011-11-29 14:02:53Z bradbell $ */
# ifndef CPPAD_TIME_TEST_INCLUDED
# define CPPAD_TIME_TEST_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin time_test$$
$spell
	gettimeofday
	vec
	cppad.hpp
	Microsoft
	namespace
	std
	const
	cout
	ctime
	ifdef
	const
	endif
	cpp
$$

$index time_test$$
$index test, speed$$
$index test, time$$

$section Determine Amount of Time to Execute a Test$$

$head Syntax$$
$codei%# include <cppad/time_test.hpp>
%$$
$icode%time% = time_test(%test%, %time_min%)%$$

$head Purpose$$
The $code time_test$$ function executes a timing test
and reports the amount of wall clock time for execution.

$head Motivation$$
It is important to separate small calculation units
and test them individually.
This way individual changes can be tested in the context of the
routine that they are in.
On many machines, accurate timing of a very short execution 
sequences is not possible.
In addition,
there may be set up and tear down time for a test that
we do not really want included in the timing.
For this reason $code time_test$$
automatically determines how many times to 
repeat the section of the test that we wish to time.

$head Include$$
The file $code cppad/time_test.hpp$$ defines the 
$code time_test$$ function.
This file is included by $code cppad/cppad.hpp$$
and it can also be included separately with out the rest of 
the $code CppAD$$ routines.

$head test$$
The $code time_test$$ argument $icode test$$ is a function,
or function object, that supports the syntax
$codei%
	%test%(%repeat%)
%$$
and its return value is $code void$$.

$subhead repeat$$
The $icode test$$ argument $icode repeat$$ has prototype
$codei%
	size_t %repeat%
%$$
It specifies the number of times to repeat the test.

$head time_min$$
The argument $icode time_min$$ has prototype
$codei%
	double %time_min%
%$$
It specifies the minimum amount of time in seconds
that the $icode test$$ routine should take.
The $icode repeat$$ argument to $icode test$$ is increased
until this amount of execution time (or more) is reached.

$head time$$
The return value $icode time$$ has prototype
$codei%
	double %time%
%$$ 
and is the number of wall clock seconds that it took
to execute $icode test$$ divided by the value used for $icode repeat$$.

$head Timing$$
The routine $cref elapsed_seconds$$ will be used to determine the
amount of time it took to execute the test.

$children%
	cppad/elapsed_seconds.hpp%
	speed/example/time_test.cpp
%$$
$head Example$$
The routine $cref time_test.cpp$$ is an example and test
of $code time_test$$.

$end
-----------------------------------------------------------------------
*/

# include <cstddef>
# include <cmath>
# include <cppad/elapsed_seconds.hpp>
# include <cppad/local/define.hpp>

CPPAD_BEGIN_NAMESPACE
/*!
\file time_test.hpp
\brief Function that preforms one timing test (for speed of execution).
*/

/*!
Preform one wall clock execution timing test.

\tparam Test
Either the type <code>void (*)(size_t)</code> or a function object
type that supports the same syntax. 

\param test
The function, or function object, that supports the operation
<code>test(repeat)</code> where \c repeat is the number of times
to repeat the tests operaiton that is being timed.

\param time_min
is the minimum amount of time that \c test should take to preform
the repetitions of the operation being timed.
*/
template <class Test>
double time_test(Test test, double time_min )
{
	size_t repeat = 0;
	double s0     = elapsed_seconds();
	double s1     = s0;
	while( s1 - s0 < time_min )
	{	repeat = std::max(size_t(1), 2 * repeat);
		s0     = elapsed_seconds();
		test(repeat);
		s1     = elapsed_seconds();
	}
	double time = (s1 - s0) / double(repeat);
	return time;
}

CPPAD_END_NAMESPACE

// END PROGRAM
# endif

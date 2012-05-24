/* $Id: memory_leak.hpp 2233 2011-12-20 19:34:24Z bradbell $ */
# ifndef CPPAD_MEMORY_LEAK_INCLUDED
# define CPPAD_MEMORY_LEAK_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
/*
$begin memory_leak$$
$spell
	alloc
	hpp
	bool
	inuse
$$

$section Memory Leak Detection$$
$index memory_leak$$
$index leak, memory$$
$index check, memory leak$$

$head Syntax$$
$icode%flag% = %memory_leak()%$$

$head flag$$
The return value $icode flag$$ has prototype
$codei%
	bool %flag%
%$$

$head Purpose$$
This routine checks that the are no memory leaks 
caused by improper use of $cref thread_alloc$$ memory allocator.
The deprecated memory allocator $cref/TrackNewDel/$$ is also checked.
Memory errors in the deprecated allocator $cref omp_alloc$$ are
reported as being in $code thread_alloc$$.

$head inuse$$
It is assumed that, when $code memory_leak$$ is called,
there should not be any memory
$cref/inuse/ta_inuse/$$ or $cref omp_inuse$$ for any thread.
If there is, a message is printed and $code memory_leak$$ returns false.

$head available$$
It is assumed that, when $code memory_leak$$ is called,
there should not be any memory
$cref/available/ta_available/$$ or $cref omp_available$$ for any thread;
i.e., it all has been returned to the system.
If there is memory still available for a thread,
$code memory_leak$$ returns false. 

$head TRACK_COUNT$$
It is assumed that, when $code memory_leak$$ is called,
$cref/TrackCount/TrackNewDel/TrackCount/$$ will return a zero value.
If it returns a non-zero value, 
$code memory_leak$$ returns false.

$head Error Message$$
If this routine returns false, it prints a message
to standard output describing the condition before returning false.

$head Multi-Threading$$
This routine cannot be used in $cref/in_parallel/ta_in_parallel/$$
execution mode.

$end
*/
# include <iostream>
# include <cppad/local/define.hpp>
# include <cppad/omp_alloc.hpp>
# include <cppad/thread_alloc.hpp>
# include <cppad/track_new_del.hpp>

CPPAD_BEGIN_NAMESPACE
/*!
\file memory_leak.hpp
File that implements a memory check at end of a CppAD program
*/

/*!
Function that checks 
allocator \c thread_alloc for misuse that results in memory leaks.
The deprecated routines 
in track_new_del.hpp and omp_alloc.hpp are also checked.

\return
returns \c true, if no error is detected and \c false otherwise.

\par
If an error is detected, diagnostic information is printed to standard
output.
*/
inline bool memory_leak(void)
{
	bool leak = false;
	size_t thread;
	using std::cout;
	using std::endl;
	using CppAD::thread_alloc;
	using CppAD::omp_alloc;
	// --------------------------------------------------------------------
	// check thead_alloc
	CPPAD_ASSERT_KNOWN(
		! thread_alloc::in_parallel(),
		"attempt to use thread_leak in parallel execution mode."
	);
	CPPAD_ASSERT_KNOWN(
		! thread_alloc::in_parallel(),
		"attempt to use thread_leak in parallel execution mode."
	);
	CPPAD_ASSERT_KNOWN(
		thread_alloc::num_threads() == 1,
		"attempt to use thread_leak while num_threads > 1."
	);
	for(thread = 0; thread < CPPAD_MAX_NUM_THREADS; thread++)
	{
		// check that no memory is currently in use for this thread
		size_t num_bytes = thread_alloc::inuse(thread);
		if( num_bytes != 0 )
		{	leak = true;
			cout << "thread_alloc::inuse(thread) = ";
			cout << num_bytes << endl;
		}
		// check that no memory is currently available for this thread
		num_bytes = thread_alloc::available(thread);
		if( num_bytes != 0 )
		{	leak = true;
			cout << "thread_alloc::available(thread) = ";
			cout << num_bytes << endl;
		}
	}
	// ----------------------------------------------------------------------
	// check track_new_del
	if( CPPAD_TRACK_COUNT() != 0 )
	{	leak = true;
		CppAD::TrackElement::Print();
	}
	return leak;
}

CPPAD_END_NAMESPACE
# endif

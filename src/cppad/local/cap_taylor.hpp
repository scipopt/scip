/* $Id: cap_taylor.hpp 2506 2012-10-24 19:36:49Z bradbell $ */
# ifndef CPPAD_CAP_TAYLOR_INCLUDED
# define CPPAD_CAP_TAYLOR_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-12 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Eclipse Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin capacity_taylor$$
$spell
	var
	taylor_
$$

$index Forward, capacity$$
$index capacity_taylor$$
$index capacity, Forward$$
$index control, memory$$
$index memory, control$$

$section Controlling Taylor Coefficients Memory Allocation$$

$head Syntax$$
$icode%f%.capacity_taylor(%c%)%$$

$subhead See Also$$
$cref seq_property$$

$head Purpose$$
The Taylor coefficients calculated by $cref Forward$$ mode calculations
are retained in an $cref ADFun$$ object for subsequent use during 
$cref Reverse$$ mode or higher order Forward mode calculations.
To be specific, a call to $cref/Forward/ForwardAny/$$ with the syntax
$codei%
        %y_p% = %f%.Forward(%p%, %x_p%)
%$$
uses the lower order Taylor coefficients and 
computes the $th p$$ order Taylor coefficients for all
the variables in the operation sequence corresponding to $icode f$$.
You can determine the number of variables in the operation sequence
using the $cref/size_var/seq_property/size_var/$$ function.
The $code capacity_taylor$$ operation allows you to control that 
amount of memory that is retained by an AD function object 
(to hold $code Forward$$ results for subsequent calculations).

$head f$$
The object $icode f$$ has prototype
$codei%
	ADFun<%Base%> %f%
%$$

$head c$$
The argument $icode c$$ has prototype
$codei%
	size_t %c%
%$$
It specifies the number of Taylor coefficients that are allocated for
each variable in the AD operation sequence corresponding to $icode f$$.


$subhead Pre-Allocating Memory$$
If you plan to make calls to $code Forward$$ with the maximum value of 
$icode p$$ equal to $icode q$$,
it should be faster to pre-allocate memory for these calls using
$codei%
	%f%.capacity_taylor(%c%)
%$$
with $icode c$$ equal to $latex q + 1$$.
If you do no do this, $code Forward$$ will automatically allocate memory
and will copy the results to a larger buffer, when necessary.
$pre

$$
Note that each call to $cref Dependent$$ frees the old memory
connected to the function object and sets the corresponding 
taylor capacity to zero.

$subhead Freeing Memory$$
If you no longer need the Taylor coefficients of order $icode q$$
and higher (that are stored in $icode f$$), 
you can reduce the memory allocated to $icode f$$ using
$codei%
	%f%.capacity_taylor(%c%)
%$$
with $icode c$$ equal to $icode q$$.
Note that, if $cref ta_hold_memory$$ is true, this memory is not actually 
returned to the system, but rather held for future use by the same thread. 

$head Original State$$
If $icode f$$ is $cref/constructed/FunConstruct/$$ with the syntax
$codei%
	ADFun<%Base%> %f%(%x%, %y%)
%$$,
there is an implicit call to $code Forward$$ with $icode p$$ equal to zero
and $icode x_p$$ equal to 
the value of the
$cref/independent variables/glossary/Tape/Independent Variable/$$ 
when the AD operation sequence was recorded.

$children%
	example/capacity_taylor.cpp
%$$
$head Example$$
The file 
$cref capacity_taylor.cpp$$
contains an example and test of these operations.
It returns true if it succeeds and false otherwise.

$end
-----------------------------------------------------------------------------
*/

// BEGIN CppAD namespace
namespace CppAD {

template <typename Base>
void ADFun<Base>::capacity_taylor(size_t c)
{	// temporary indices
	size_t i, j, p;

	if( c == taylor_col_dim_ )
		return;

	if( c == 0 )
	{	taylor_.free();
		taylor_per_var_ = 0;
		taylor_col_dim_ = 0;
		return;
	}
	
	// Allocate new matrix with requested number of columns 
	size_t new_len   = c * total_num_var_;
	pod_vector<Base> new_taylor;
	new_taylor.extend(new_len);

	// number of columns to copy
	p = std::min(taylor_per_var_, c);

	// copy the old data into the new matrix
	CPPAD_ASSERT_UNKNOWN( (taylor_per_var_ == 0) | (taylor_.size() != 0) );
	for(i = 0; i < total_num_var_; i++)
	{	for(j = 0; j < p; j++)
		{	new_taylor[i * c + j]  = taylor_[i * taylor_col_dim_ + j];
		}
	}

	// replace taylor_ by new_taylor
	taylor_.swap(new_taylor);
	taylor_col_dim_ = c;
	taylor_per_var_ = p;

	// note that the destructor for new_taylor will free the old taylor memory
	return;
}

} // END CppAD namespace
	

# endif

/* $Id: dependent.hpp 2991 2013-10-22 16:25:15Z bradbell $ */
# ifndef CPPAD_DEPENDENT_INCLUDED
# define CPPAD_DEPENDENT_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-13 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Eclipse Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
/*
$begin Dependent$$
$spell 
	alloc
	num
	taylor_
	ADvector
	const
$$

$spell
$$

$section Stop Recording and Store Operation Sequence$$

$index ADFun, operation sequence$$
$index operation, sequence store$$
$index sequence, operation store$$
$index recording, stop$$
$index tape, stop recording$$
$index Dependent$$

$head Syntax$$
$icode%f%.Dependent(%x%, %y%)%$$

$head Purpose$$
Stop recording and the AD of $icode Base$$
$cref/operation sequence/glossary/Operation/Sequence/$$
that started with the call
$codei%
	Independent(%x%)
%$$
and store the operation sequence in $icode f$$.
The operation sequence defines an 
$cref/AD function/glossary/AD Function/$$
$latex \[
	F : B^n \rightarrow B^m
\] $$
where $latex B$$ is the space corresponding to objects of type $icode Base$$.
The value $latex n$$ is the dimension of the 
$cref/domain/seq_property/Domain/$$ space for the operation sequence.
The value $latex m$$ is the dimension of the 
$cref/range/seq_property/Range/$$ space for the operation sequence
(which is determined by the size of $icode y$$).

$head f$$
The object $icode f$$ has prototype
$codei%
	ADFun<%Base%> %f%
%$$
The AD of $icode Base$$ operation sequence is stored in $icode f$$; i.e.,
it becomes the operation sequence corresponding to $icode f$$.
If a previous operation sequence was stored in $icode f$$,
it is deleted. 

$head x$$
The argument $icode x$$ 
must be the vector argument in a previous call to
$cref Independent$$.
Neither its size, or any of its values, are allowed to change
between calling
$codei%
	Independent(%x%)
%$$
and 
$codei%
	%f%.Dependent(%x%, %y%)
%$$.

$head y$$
The vector $icode y$$ has prototype
$codei%
	const %ADvector% &%y%
%$$
(see $cref/ADvector/FunConstruct/$$ below).
The length of $icode y$$ must be greater than zero
and is the dimension of the range space for $icode f$$.

$head ADvector$$
The type $icode ADvector$$ must be a $cref SimpleVector$$ class with
$cref/elements of type/SimpleVector/Elements of Specified Type/$$
$codei%AD<%Base%>%$$.
The routine $cref CheckSimpleVector$$ will generate an error message
if this is not the case.

$head Taping$$
The tape,
that was created when $codei%Independent(%x%)%$$ was called, 
will stop recording.
The AD operation sequence will be transferred from
the tape to the object $icode f$$ and the tape will then be deleted.

$head Forward$$
No $cref Forward$$ calculation is preformed during this operation.
Thus, directly after this operation,
$codei%
	%f%.size_taylor()
%$$ 
is zero (see $cref size_taylor$$).

$head Parallel Mode$$
$index parallel, Dependent$$
$index Dependent, parallel$$
The call to $code Independent$$,
and the corresponding call to
$codei%
	ADFun<%Base%> %f%( %x%, %y%)
%$$
or 
$codei%
	%f%.Dependent( %x%, %y%)
%$$
or $cref abort_recording$$,
must be preformed by the same thread; i.e.,
$cref/thread_alloc::thread_num/ta_thread_num/$$ must be the same.

$head Example$$
The file
$cref fun_check.cpp$$ 
contains an example and test of this operation.
It returns true if it succeeds and false otherwise.

$end
----------------------------------------------------------------------------
*/


// BEGIN CppAD namespace
namespace CppAD {

/*!
Determine the \c tape corresponding to this exeuction thread and then use
<code>Dependent(tape, y)</code> to store this tapes recording in a function.

\param y [in]
The dependent variable vector for the corresponding function.
*/
template <typename Base>
template <typename ADvector>
void ADFun<Base>::Dependent(const ADvector &y)
{	ADTape<Base>* tape = AD<Base>::tape_ptr();
	CPPAD_ASSERT_KNOWN(
		tape != CPPAD_NULL,
		"Can't store current operation sequence in this ADFun object"
		"\nbecause there is no active tape (for this thread)."
	);

	// code above just determines the tape and checks for errors
	Dependent(tape, y);
}


/*!
Determine the \c tape corresponding to this exeuction thread and then use
<code>Dependent(tape, y)</code> to store this tapes recording in a function.

\param x [in]
The independent variable vector for this tape. This informaiton is
also stored in the tape so a check is done to make sure it is correct
(if NDEBUG is not defined).

\param y [in]
The dependent variable vector for the corresponding function.
*/
template <typename Base>
template <typename ADvector>
void ADFun<Base>::Dependent(const ADvector &x, const ADvector &y)
{
	CPPAD_ASSERT_KNOWN(
		x.size() > 0,
		"Dependent: independent variable vector has size zero."
	);
	CPPAD_ASSERT_KNOWN(
		Variable(x[0]),
		"Dependent: independent variable vector has been changed."
	);
	ADTape<Base> *tape = AD<Base>::tape_ptr(x[0].tape_id_);
	CPPAD_ASSERT_KNOWN(
		tape->size_independent_ == size_t( x.size() ),
		"Dependent: independent variable vector has been changed."
	);
# ifndef NDEBUG
	size_t i, j;
	for(j = 0; j < size_t(x.size()); j++)
	{	CPPAD_ASSERT_KNOWN(
		size_t(x[j].taddr_) == (j+1),
		"ADFun<Base>: independent variable vector has been changed."
		);
		CPPAD_ASSERT_KNOWN(
		x[j].tape_id_ == x[0].tape_id_,
		"ADFun<Base>: independent variable vector has been changed."
		);
	}
	for(i = 0; i < size_t(y.size()); i++)
	{	CPPAD_ASSERT_KNOWN(
		CppAD::Parameter( y[i] ) | (y[i].tape_id_ == x[0].tape_id_) ,
		"ADFun<Base>: dependent vector contains a variable for"
		"\na different tape (thread) than the independent variables."
		);
	}
# endif

	// code above just determines the tape and checks for errors
	Dependent(tape, y);
}

/*!
Replace the floationg point operations sequence for this function object.

\param tape
is a tape that contains the new floating point operation sequence
for this function.
After this operation, all memory allocated for this tape is deleted.

\param y
The dependent variable vector for the function being stored in this object.
*/

template <typename Base>
template <typename ADvector>
void ADFun<Base>::Dependent(ADTape<Base> *tape, const ADvector &y)
{
	size_t   m = y.size();
	size_t   n = tape->size_independent_;
	size_t   i, j;
	size_t   y_taddr;

	// check ADvector is Simple Vector class with AD<Base> elements
	CheckSimpleVector< AD<Base>, ADvector>();

	CPPAD_ASSERT_KNOWN(
		y.size() > 0,
		"ADFun operation sequence dependent variable size is zero size"
	); 

	// set total number of variables in tape, parameter flag, 
	// make a tape copy of dependent variables that are parameters, 
	// and store tape address for each dependent variable
	CPPAD_ASSERT_UNKNOWN( NumRes(ParOp) == 1 );
	dep_parameter_.resize(m);
	dep_taddr_.resize(m);
	for(i = 0; i < m; i++)
	{	dep_parameter_[i] = CppAD::Parameter(y[i]);
		if( dep_parameter_[i] )
			y_taddr = tape->RecordParOp( y[i].value_ );
		else	y_taddr = y[i].taddr_;

		CPPAD_ASSERT_UNKNOWN( y_taddr > 0 );
		dep_taddr_[i] = y_taddr;
	}

	// put an EndOp at the end of the tape
	tape->Rec_.PutOp(EndOp);

	// total number of variables on the tape
	total_num_var_ = tape->Rec_.num_rec_var();

	// conditional skip vector
	cskip_op_.clear();
	cskip_op_.resize( tape->Rec_.num_rec_op() );
	for(i = 0; i < cskip_op_.size(); i++)
		cskip_op_[i] = false;

	// now that each dependent variable has a place in the tape,
	// and there is a EndOp at the end of the tape, we can transfer the 
	// recording to the player and and erase the tape.
	play_.get(tape->Rec_);

	// now we can delete the tape
	AD<Base>::tape_manage(tape_manage_delete);

	// total number of varables in this recording 
	CPPAD_ASSERT_UNKNOWN( total_num_var_ == play_.num_rec_var() );

	// used to determine if there is an operation sequence in *this
	CPPAD_ASSERT_UNKNOWN( total_num_var_ > 0 );

	// free old buffers
	for_jac_sparse_pack_.resize(0, 0);
	for_jac_sparse_set_.resize(0,0);

	// initial row and column dimensions
	taylor_.erase();
	taylor_per_var_   = 0;
	taylor_col_dim_   = 0;

	// set tape address 
	ind_taddr_.resize(n);
	CPPAD_ASSERT_UNKNOWN(
		n < total_num_var_
	);
	for(j = 0; j < n; j++)
	{	CPPAD_ASSERT_UNKNOWN( play_.GetOp(j+1) == InvOp );
		ind_taddr_[j] = j+1;
	}

}

} // END CppAD namespace

# endif

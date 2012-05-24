/* $Id: recorder.hpp 2082 2011-08-31 17:50:58Z bradbell $ */
# ifndef CPPAD_RECORDER_INCLUDED
# define CPPAD_RECORDER_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
# include <cppad/local/hash_code.hpp>
# include <cppad/local/pod_vector.hpp>

CPPAD_BEGIN_NAMESPACE
/*!
\file recorder.hpp
File used to define the recorder class.
*/

/*!
Class used to store an operation sequence while it is being recorded
(the operation sequence is copied to the player class for playback).

\tparam Base
This is an AD< \a Base > operation sequence recording; i.e.,
it records operations of type AD< \a Base >.
*/
template <class Base>
class recorder {
	friend class player<Base>;

private:
	/// offset for this thread in the static hash table
	const size_t thread_offset_;

	/// Number of variables in the recording.
	size_t    num_rec_var_;

	/// The operators in the recording.
	pod_vector<CPPAD_OP_CODE_TYPE> rec_op_;

	/// The VecAD indices in the recording.
	pod_vector<addr_t> rec_vecad_ind_;

	/// The argument indices in the recording
	pod_vector<addr_t> rec_op_arg_;

	/// The parameters in the recording.
	/// Note that Base may not be plain old data, so use false in consructor.
	pod_vector<Base> rec_par_;

	/// Character strings ('\\0' terminated) in the recording.
	pod_vector<char> rec_text_;
// ---------------------- Public Functions -----------------------------------
public:
	/// Default constructor
	recorder(void) : 
	thread_offset_( thread_alloc::thread_num() * CPPAD_HASH_TABLE_SIZE ) ,
	num_rec_var_(0)                                      ,
	rec_op_( std::numeric_limits<addr_t>::max() )        ,
	rec_vecad_ind_( std::numeric_limits<addr_t>::max() ) ,
	rec_op_arg_( std::numeric_limits<addr_t>::max() )    ,
	rec_par_( std::numeric_limits<addr_t>::max() )       ,
	rec_text_( std::numeric_limits<addr_t>::max() )
	{ }

	/// Destructor
	~recorder(void)
	{ }

	/*!
	Erase all information in recording.

	Erases the operation sequence store in this recording 
	(the operation sequence is empty after this operation).
	The buffers used to store the current recording are returned
	to the system (so as to conserve on memory).
	*/
	void Erase(void)
	{	num_rec_var_  = 0;
		rec_op_.erase();
		rec_vecad_ind_.erase();
		rec_op_arg_.erase();
		rec_par_.erase();
		rec_text_.erase();
	}
	/// Start recording the next operator in the operation sequence.
	inline size_t PutOp(OpCode op);
	/// Add a value to the end of the current vector of VecAD indices.
	inline size_t PutVecInd(size_t vec_ind);
	/// Find or add a parameter to the current vector of parameters.
	inline size_t PutPar(const Base &par);
	/// Put one operation argument index in the recording
	inline void PutArg(addr_t arg0); 
	/// Put two operation argument index in the recording
	inline void PutArg(addr_t arg0, addr_t arg1); 
	/// Put three operation argument index in the recording
	inline void PutArg(addr_t arg0, addr_t arg1, addr_t arg2); 
	/// Put four operation argument index in the recording
	inline void PutArg(addr_t arg0, addr_t arg1, addr_t arg2, addr_t arg3); 
	/// Put five operation argument index in the recording
	inline void PutArg(addr_t arg0, addr_t arg1, addr_t arg2, addr_t arg3,
		addr_t arg4);
	/// Put six operation argument index in the recording
	inline void PutArg(addr_t arg0, addr_t arg1, addr_t arg2, addr_t arg3,
		addr_t arg4, addr_t arg5);

	/// Put a character string in the text for this recording.
	inline size_t PutTxt(const char *text);

	/// Number of variables currently stored in the recording.
	size_t num_rec_var(void) const
	{	return num_rec_var_; }

	/// Approximate amount of memory used by the recording 
	size_t Memory(void) const
	{	return rec_op_.capacity()        * sizeof(CPPAD_OP_CODE_TYPE) 
		     + rec_vecad_ind_.capacity() * sizeof(size_t)
		     + rec_op_arg_.capacity()    * sizeof(addr_t)
		     + rec_par_.capacity()       * sizeof(Base)
		     + rec_text_.capacity()      * sizeof(char);
	}
};

/*!
Start recording the next operator in the operation sequence.

This sets the op code for the next operation in this recording.
This call must be followed by putting the corresponding 
\verbatim
	NumArg(op)
\endverbatim
argument indices in the recording.

\param op
Is the op code corresponding to the the operation that is being
recorded. 

\return
The return value is the index of the primary (last) variable 
corresponding to the result of this operation. 
The number of variables corresponding to the operation is given by
\verbatim
	NumRes(op)
\endverbatim
With each call to PutOp 
the return index increases by the number of variables corresponding
to this call to PutOp.
This index starts at zero after the default constructor
and after each call to Erase.
*/
template <class Base>
inline size_t recorder<Base>::PutOp(OpCode op)
{	size_t i    = rec_op_.extend(1);
	rec_op_[i]  = static_cast<CPPAD_OP_CODE_TYPE>(op);
	CPPAD_ASSERT_UNKNOWN( rec_op_.size() == i + 1 );

	// first operator should be a BeginOp and NumRes( BeginOp ) > 0
	num_rec_var_ += NumRes(op);
	CPPAD_ASSERT_UNKNOWN( num_rec_var_ > 0 );

	return num_rec_var_ - 1;
}

/*!
Add a value to the end of the current vector of VecAD indices.

For each VecAD vector, this routine is used to store the length
of the vector followed by the parameter index corresponding to each
value in the vector.
This value for the elements of the VecAD vector corresponds to the
beginning of the operation sequence.

\param vec_ind
is the index to be palced at the end of the vector of VecAD indices.

\return
is the index in the vector of VecAD indices corresponding to this value.
This index starts at zero after the recorder default constructor
and after each call to Erase.
It increments by one for each call to PutVecInd..
*/
template <class Base>
inline size_t recorder<Base>::PutVecInd(size_t vec_ind)
{	size_t i          = rec_vecad_ind_.extend(1);
	rec_vecad_ind_[i] = vec_ind;
	CPPAD_ASSERT_UNKNOWN( rec_vecad_ind_.size() == i + 1 );

	return i;
}

/*!
Find or add a parameter to the current vector of parameters.

\param par
is the parameter to be found or placed in the vector of parameters.

\return
is the index in the parameter vector corresponding to this parameter value.
This value is not necessarily placed at the end of the vector
(because values that are identically equal may be reused).
*/
template <class Base>
size_t recorder<Base>::PutPar(const Base &par)
{	static size_t   hash_table[CPPAD_HASH_TABLE_SIZE * CPPAD_MAX_NUM_THREADS];
	size_t          i;
	size_t          code;

	CPPAD_ASSERT_UNKNOWN( 
		thread_offset_ / CPPAD_HASH_TABLE_SIZE
		== 
		thread_alloc::thread_num() 
	);

	// get hash code for this value
	code = static_cast<size_t>( hash_code(par) );
	CPPAD_ASSERT_UNKNOWN( code < CPPAD_HASH_TABLE_SIZE );

	// If we have a match, return the parameter index
	i = hash_table[code + thread_offset_];
	if( i < rec_par_.size() && IdenticalEqualPar(rec_par_[i], par) )
			return i;
	
	// place a new value in the table
	i           = rec_par_.extend(1);
	rec_par_[i] = par;
	CPPAD_ASSERT_UNKNOWN( rec_par_.size() == i + 1 );

	// make the hash code point to this new value
	hash_table[code + thread_offset_] = i;

	// return the parameter index
	return i;
}
// -------------------------- PutArg --------------------------------------
/*!
Prototype for putting operation argument indices in the recording.

The following syntax
\verbatim
	rec.PutArg(arg0)
	rec.PutArg(arg0, arg1)
	.
	.
	.
	rec.PutArg(arg0, arg1, ..., arg5)
\endverbatim
places the values passed to PutArg at the current end of the
operation argument indices for the recording.
\a arg0 comes before \a arg1, etc.
The proper number of operation argument indices 
corresponding to the operation code op is given by
\verbatim
	NumArg(op)
\endverbatim
The number of the operation argument indices starts at zero
after the default constructor and each call to Erase.
It increases by the number of indices placed by each call to PutArg.
*/
inline void prototype_put_arg(void)
{	// This routine should not be called
	CPPAD_ASSERT_UNKNOWN(false);
}
/*!
Put one operation argument index in the recording

\param arg0
The operation argument index

\copydetails prototype_put_arg 
*/
template <class Base>
inline void recorder<Base>::PutArg(addr_t arg0)
{ 
	size_t i       = rec_op_arg_.extend(1);
	rec_op_arg_[i] =  static_cast<addr_t>( arg0 );
	CPPAD_ASSERT_UNKNOWN( rec_op_arg_.size() == i + 1 );
}
/*!
Put two operation argument index in the recording

\param arg0
First operation argument index.

\param arg1
Second operation argument index.

\copydetails prototype_put_arg 
*/
template <class Base>
inline void recorder<Base>::PutArg(addr_t arg0, addr_t arg1)
{ 
	size_t i         = rec_op_arg_.extend(2);
	rec_op_arg_[i++] =  static_cast<addr_t>( arg0 );
	rec_op_arg_[i]   =  static_cast<addr_t>( arg1 );
	CPPAD_ASSERT_UNKNOWN( rec_op_arg_.size() == i + 1 );
}
/*!
Put three operation argument index in the recording

\param arg0
First operation argument index.

\param arg1
Second operation argument index.

\param arg2
Third operation argument index.

\copydetails prototype_put_arg 
*/
template <class Base>
inline void recorder<Base>::PutArg(addr_t arg0, addr_t arg1, addr_t arg2)
{ 
	size_t i         = rec_op_arg_.extend(3);
	rec_op_arg_[i++] =  static_cast<addr_t>( arg0 );
	rec_op_arg_[i++] =  static_cast<addr_t>( arg1 );
	rec_op_arg_[i]   =  static_cast<addr_t>( arg2 );
	CPPAD_ASSERT_UNKNOWN( rec_op_arg_.size() == i + 1 );
}
/*!
Put four operation argument index in the recording

\param arg0
First operation argument index.

\param arg1
Second operation argument index.

\param arg2
Third operation argument index.

\param arg3
Fourth operation argument index.

\copydetails prototype_put_arg 
*/
template <class Base>
inline void recorder<Base>::PutArg(addr_t arg0, addr_t arg1, addr_t arg2,
	addr_t arg3)
{ 
	size_t i         = rec_op_arg_.extend(4);
	rec_op_arg_[i++] =  static_cast<addr_t>( arg0 );
	rec_op_arg_[i++] =  static_cast<addr_t>( arg1 );
	rec_op_arg_[i++] =  static_cast<addr_t>( arg2 );
	rec_op_arg_[i]   =  static_cast<addr_t>( arg3 );
	CPPAD_ASSERT_UNKNOWN( rec_op_arg_.size() == i + 1 );

}
/*!
Put five operation argument index in the recording

\param arg0
First operation argument index.

\param arg1
Second operation argument index.

\param arg2
Third operation argument index.

\param arg3
Fourth operation argument index.

\param arg4
Fifth operation argument index.

\copydetails prototype_put_arg 
*/
template <class Base>
inline void recorder<Base>::PutArg(addr_t arg0, addr_t arg1, addr_t arg2,
	addr_t arg3, addr_t arg4)
{ 
	size_t i         = rec_op_arg_.extend(5);
	rec_op_arg_[i++] =  static_cast<addr_t>( arg0 );
	rec_op_arg_[i++] =  static_cast<addr_t>( arg1 );
	rec_op_arg_[i++] =  static_cast<addr_t>( arg2 );
	rec_op_arg_[i++] =  static_cast<addr_t>( arg3 );
	rec_op_arg_[i]   =  static_cast<addr_t>( arg4 );
	CPPAD_ASSERT_UNKNOWN( rec_op_arg_.size() == i + 1 );

}
/*!
Put six operation argument index in the recording

\param arg0
First operation argument index.

\param arg1
Second operation argument index.

\param arg2
Third operation argument index.

\param arg3
Fourth operation argument index.

\param arg4
Fifth operation argument index.

\param arg5
Sixth operation argument index.

\copydetails prototype_put_arg 
*/
template <class Base>
inline void recorder<Base>::PutArg(addr_t arg0, addr_t arg1, addr_t arg2, 
	addr_t arg3, addr_t arg4, addr_t arg5)
{ 
	size_t i         = rec_op_arg_.extend(6);
	rec_op_arg_[i++] =  static_cast<addr_t>( arg0 );
	rec_op_arg_[i++] =  static_cast<addr_t>( arg1 );
	rec_op_arg_[i++] =  static_cast<addr_t>( arg2 );
	rec_op_arg_[i++] =  static_cast<addr_t>( arg3 );
	rec_op_arg_[i++] =  static_cast<addr_t>( arg4 );
	rec_op_arg_[i]   =  static_cast<addr_t>( arg5 );
	CPPAD_ASSERT_UNKNOWN( rec_op_arg_.size() == i + 1 );
}
// --------------------------------------------------------------------------
/*!
Put a character string in the text for this recording.

\param text
is a '\\0' terminated character string that is to be put in the
vector of characters corresponding to this recording.
The terminator '\\0' will be included.

\return
is the offset with in the text vector for this recording at which
the character string starts.
*/
template <class Base>
inline size_t recorder<Base>::PutTxt(const char *text)
{
	// determine length of the text including terminating '\0'
	size_t n = 0;
	while( text[n] != '\0' )
		n++;
	CPPAD_ASSERT_UNKNOWN( n <= 1000 ); 
	n++;
	CPPAD_ASSERT_UNKNOWN( text[n-1] == '\0' );

	// copy text including terminating '\0'
	size_t i = rec_text_.extend(n); 
	size_t j;
	for(j = 0; j < n; j++)
		rec_text_[i + j] = text[j];
	CPPAD_ASSERT_UNKNOWN( rec_text_.size() == i + n );

	return i;
}
// -------------------------------------------------------------------------


CPPAD_END_NAMESPACE
# endif

/* $Id: player.hpp 1986 2011-06-18 20:33:17Z bradbell $ */
# ifndef CPPAD_PLAYER_INCLUDED
# define CPPAD_PLAYER_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

CPPAD_BEGIN_NAMESPACE
/*!
\file player.hpp
File used to define the player class.
*/


/*!
Class used to store and play back an operation sequence recording.

\tparam Base
These were AD< \a Base > operations when recorded. Operations during playback
are done using the type  \a Base .
*/
template <class Base>
class player {

// -------------- Variables that define the recording -----------------------
private:
	/// Number of variables in the recording.
	size_t    num_rec_var_;

	/// The operators in the recording.
	pod_vector<CPPAD_OP_CODE_TYPE> rec_op_;

	/// Number of VecAD vectors in the recording
	size_t    num_rec_vecad_vec_;

	/// The VecAD indices in the recording.
	pod_vector<addr_t> rec_vecad_ind_;

	/// The operation argument indices in the recording
	pod_vector<addr_t> rec_op_arg_;

	/// The parameters in the recording.
	/// Note that Base may not be plain old data, so use false in consructor.
	pod_vector<Base> rec_par_;

	/// Character strings ('\\0' terminated) in the recording.
	pod_vector<char> rec_text_;


// --------------- Functions used to create and maniplate a recording -------
public:
	/// Default constructor
	player(void) : 
	num_rec_var_(0)                                      ,
	rec_op_( std::numeric_limits<addr_t>::max() )        ,
	rec_vecad_ind_( std::numeric_limits<addr_t>::max() ) ,
	rec_op_arg_( std::numeric_limits<addr_t>::max() )    ,
	rec_par_( std::numeric_limits<addr_t>::max() )       ,
	rec_text_( std::numeric_limits<addr_t>::max() )
	{ }

	/// Destructor
	~player(void)
	{ }

	// ===============================================================
	// Begin two functions with idential code but different argument types.
	/*!  
 	Moving an operation sequence from a recorder to a player
 
	\param rec
	the object that was used to record the operation sequence.
 	*/
	void operator=(const recorder<Base> &rec)
	{	size_t i;

		// Var
		num_rec_var_        = rec.num_rec_var_;

		// Op
		rec_op_             = rec.rec_op_;

		// VecInd
		rec_vecad_ind_      = rec.rec_vecad_ind_;

		// Arg
		rec_op_arg_         = rec.rec_op_arg_;

		// Par
		rec_par_            = rec.rec_par_;

		// Txt
		rec_text_           = rec.rec_text_;

		// set the number of VecAD vectors
		num_rec_vecad_vec_ = 0;
		for(i = 0; i < rec_vecad_ind_.size(); i += rec_vecad_ind_[i] + 1)
			num_rec_vecad_vec_++;
		// rec_vecad_ind_ contains size of each VecAD followed by
		// the parameter indices used to iniialize it.
		CPPAD_ASSERT_UNKNOWN( i == rec_vecad_ind_.size() );
	}

	/*!  
 	Copying an operation sequence from one player to another
 
	\param play
	the object that contains the operatoion sequence to copy.
 	*/
	void operator=(const player& play)
	{	
		// Var
		num_rec_var_        = play.num_rec_var_;

		// Op
		rec_op_             = play.rec_op_;

		// VecInd
		num_rec_vecad_vec_  = play.num_rec_vecad_vec_;
		rec_vecad_ind_      = play.rec_vecad_ind_;

		// Arg
		rec_op_arg_         = play.rec_op_arg_;

		// Par
		rec_par_            = play.rec_par_;

		// Txt
		rec_text_           = play.rec_text_;
	}
	// End two functions with idential code but different argument types.
	// ===============================================================

	/// Erase all information in an operation sequence player.
	void Erase(void)
	{	
		num_rec_var_       = 0;
		num_rec_vecad_vec_ = 0;

		rec_op_.erase();
		rec_vecad_ind_.erase();
		rec_op_arg_.erase();
		rec_par_.erase();
		rec_text_.erase();
	}

// ------------------ Old method of palying back a recording -----------------
public:
	/*! 
	\brief 
	Fetch an operator from the recording.

	\return 
	the i-th operator in the recording.

	\param i
	the index of the operator in recording
	*/
	OpCode GetOp (size_t i) const
	{	return OpCode(rec_op_[i]); }

	/*! 
	\brief 
	Fetch a VecAD index from the recording.

	\return 
	the i-th VecAD index in the recording.

	\param i
	the index of the VecAD index in recording
	*/
	size_t GetVecInd (size_t i) const
	{	return rec_vecad_ind_[i]; }

	/*! 
	\brief 
	Fetch a parameter from the recording.

	\return 
	the i-th parameter in the recording.

	\param i
	the index of the parameter in recording
	*/
	Base GetPar(size_t i) const
	{	return rec_par_[i]; }

	/*! 
	\brief 
	Fetch entire parameter vector from the recording.

	\return 
	the entire parameter vector.

	*/
	const Base* GetPar(void) const
	{	return rec_par_.data(); }

	/*! 
	\brief 
	Fetch a '\\0' terminated string from the recording.

	\return 
	the beginning of the string.

	\param i
	the index where the string begins. 
	*/
	const char *GetTxt(size_t i) const
	{	CPPAD_ASSERT_UNKNOWN(i < rec_text_.size() );
		return rec_text_.data() + i;
	}
	
	/*!
	\brief
	Replace an argument index in the recording.

	\param i
	is the index, in argument indices, that is to be replaced.

	\param value
	is the new normal index value.
	*/
	void ReplaceInd(size_t i, size_t value)
	{	rec_op_arg_[i] =  static_cast<addr_t>( value ); }

	/// Fetch number of variables in the recording.
	size_t num_rec_var(void) const
	{	return num_rec_var_; }

	/// Fetch number of operators in the recording.
	size_t num_rec_op(void) const
	{	return rec_op_.size(); }

	/// Fetch number of VecAD indices in the recording.
	size_t num_rec_vecad_ind(void) const
	{	return rec_vecad_ind_.size(); }

	/// Fetch number of VecAD vectors in the recording
	size_t num_rec_vecad_vec(void) const
	{	return num_rec_vecad_vec_; }

	/// Fetch number of argument indices in the recording.
	size_t num_rec_op_arg(void) const
	{	return rec_op_arg_.size(); }

	/// Fetch number of parameters in the recording.
	size_t num_rec_par(void) const
	{	return rec_par_.size(); }

	/// Fetch number of characters (representing strings) in the recording.
	size_t num_rec_text(void) const
	{	return rec_text_.size(); }

	/// Fetch a rough measure of amount of memory used to store recording
	/// (just lengths, not capacities). 
	size_t Memory(void) const
	{	return rec_op_.size()        * sizeof(OpCode) 
		     + rec_op_arg_.size()    * sizeof(addr_t)
		     + rec_par_.size()       * sizeof(Base)
		     + rec_text_.size()      * sizeof(char)
		     + rec_vecad_ind_.size() * sizeof(addr_t)
		;
	}

// ------------- Variables used for new methog of playing back a recording ---
private:
	/// Current operator
	OpCode    op_;

	/// Index in recording corresponding to current operator
	size_t    op_index_;

	/// Current offset of the argument indices in rec_op_arg_ 
	size_t    op_arg_;

	/// Index for primary (last) variable corresponding to current operator
	size_t    var_index_;

// ----------- Functions used in new method for palying back a recording ---
public:
	/*!
	Start a play back of the recording during a forward sweep.

	Use repeated calls to next_forward to play back one operator at a time.

	\param op
	The input value of op does not matter. Its output value is the
	first operator in the recording; i.e., BeginOp.

	\param op_arg
	The input value of *op_arg does not matter. Its output value is the
	beginning of the vector of argument indices for the first operation;
	i.e., 0

	\param op_index
	The input value of op_index does not matter. Its output value
	is the index of the next first operator in the recording; i.e., 0.

	\param var_index
	The input value of var_index does not matter. Its output value is the
	index of the primary (last) result corresponding to the the first
	operator (which must be a BeginOp); i.e., 0.
	*/
	void start_forward(
	OpCode& op, const addr_t*& op_arg, size_t& op_index, size_t& var_index)
	{
		op        = op_          = OpCode( rec_op_[0] ); 
		op_arg_   = 0;
		op_arg    = rec_op_arg_.data();
		op_index  = op_index_    = 0;
		var_index = var_index_   = 0;

		CPPAD_ASSERT_UNKNOWN( op_  == BeginOp );
		CPPAD_ASSERT_NARG_NRES(op_, 0, 1);

		return;
	}

	/*!
	Fetch the next operator during a forward sweep.

	Use start_forward to initialize to the first operator; i.e.,
	the BeginOp at the beginning of the recording. 

	\param op
	The input value of op does not matter. Its output value is the
	next operator in the recording.
	For speed, \c next_forward does not check for the special case
	of <tt>op == CSumOp</tt>. In this case, the other return values
	from \c next_forward must be corrected by a call to \c forward_csum.

	\param op_arg
	The input value of *op_arg does not matter. Its output value is the
	beginning of the vector of argument indices for this operation.

	\param op_index
	The input value of op_index does not matter. Its output value
	is the index of the next operator in the recording. Thus the ouput
	value following the previous call to start_forward is one. In addition,
	the output value increases by one with each call to next_forward. 

	\param var_index
	The input value of var_index does not matter. Its output value is the
	index of the primary (last) result corresponding to the operator op.
	*/

	void next_forward(
	OpCode& op, const addr_t*& op_arg, size_t& op_index, size_t& var_index)
	{	using CppAD::NumRes;
		using CppAD::NumArg;

		// index for the next operator 
		op_index    = ++op_index_;

		// first argument for next operator 
		op_arg_    += NumArg(op_);                   // index
		op_arg      = op_arg_ + rec_op_arg_.data();  // pointer

		// next operator
		op          = op_         = OpCode( rec_op_[ op_index_ ] );

		// index for last result for next operator
		var_index   = var_index_ += NumRes(op);

		CPPAD_ASSERT_UNKNOWN( op_arg_ + NumArg(op) <= rec_op_arg_.size() );
		CPPAD_ASSERT_UNKNOWN( var_index_  < num_rec_var_ );
	}
	/*!
	Correct \c next_forward return values when <tt>op == CSumOp</tt>.

	\param op
	The input value of op must be the return value from the previous
	call to \c next_forward and must be \c CSumOp.

	\param op_arg
	The input value of *op_arg must be the return value from the 
	previous call to \c next_forward. Its output value is the
	beginning of the vector of argument indices for this operation.

	\param op_index
	The input value of op_index does must be the return value from the
	previous call to \c next_forward. Its output value
	is the index of this operator in the recording. 

	\param var_index
	The input value of var_index must be the return value from the
	previous call to \c next_forward. Its output value is the
	index of the primary (last) result corresponding to this.
	*/
	void forward_csum(
	OpCode& op, const addr_t*& op_arg, size_t& op_index, size_t& var_index)
	{	using CppAD::NumRes;
		using CppAD::NumArg;
		CPPAD_ASSERT_UNKNOWN( op == CSumOp );
		CPPAD_ASSERT_UNKNOWN( NumArg(CSumOp) == 0 );
		CPPAD_ASSERT_UNKNOWN(
		op_arg[0] + op_arg[1] == op_arg[ 3 + op_arg[0] + op_arg[1] ]
		);
		/*
		The only thing that really needs fixing is op_arg_.
		Actual number of arugments for this operator is
			op_arg[0] + op_arg[1] + 4.
 		We must change op_arg_ so that when you add NumArg(CSumOp)
		you get first argument for next operator in sequence.
		*/
		op_arg_    += op_arg[0] + op_arg[1] + 4;

		CPPAD_ASSERT_UNKNOWN( op_arg_ + NumArg(op) <= rec_op_arg_.size() );
		CPPAD_ASSERT_UNKNOWN( var_index_  < num_rec_var_ );
	}
	/*!
	Get a non-constant version of op_arg returned by previous next_forward

	\return
	The return value is equal to the return value of op_arg 
	corresponding to the previous call to next_forward.
	*/
	addr_t* forward_non_const_arg(void)
	{	return op_arg_ + rec_op_arg_.data(); }

	/*!
	Start a play back of the recording during a reverse sweep.

	Use repeated calls to next_reverse to play back one operator at a time.

	\param op
	The input value of op does not matter. Its output value is the
	last operator in the recording; i.e., EndOp.

	\param op_arg
	The input value of *op_arg does not matter. Its output value is the
	beginning of the vector of argument indices for the last operation;
	(there are no arguments for the last operation so \a op_arg is invalid).

	\param op_index
	The input value of op_index does not matter. Its output value
	is the index of the last operator in the recording.

	\param var_index
	The input value of var_index does not matter. Its output value is the
	index of the primary (last) result corresponding to the the last
	operator (which must be a EndOp).
	(there are no results for the last operation 
	so \a var_index is invalid).
	*/

	void start_reverse(
	OpCode& op, const addr_t*& op_arg, size_t& op_index, size_t& var_index)
	{
		op_arg_     = rec_op_arg_.size();                // index
		op_arg      = op_arg_ + rec_op_arg_.data();      // pointer

		op_index    = op_index_   = rec_op_.size() - 1; 
		var_index   = var_index_  = num_rec_var_ - 1;

		op          = op_         = OpCode( rec_op_[ op_index_ ] );
		CPPAD_ASSERT_UNKNOWN( op_ == EndOp );
		CPPAD_ASSERT_NARG_NRES(op, 0, 0);
		return;
	}

	/*!
	Fetch the next operator during a reverse sweep.

	Use start_reverse to initialize to reverse play back.
	The first call to next_reverse (after start_reverse) will give the 
        last operator in the recording.

	\param op
	The input value of op does not matter. Its output value is the
	next operator in the recording (in reverse order).
	The last operator sets op equal to EndOp.

	\param op_arg
	The input value of *op_arg does not matter. Its output value is the
	beginning of the vector of argument indices for this operation.
	The last operator sets op_arg equal to the beginning of the 
	argument indices for the entire recording.
	For speed, \c next_reverse does not check for the special case
	of <tt>op == CSumOp</tt>. In this case, the other return values
	from \c next_reverse must be corrected by a call to \c reverse_csum.


	\param op_index
	The input value of op_index does not matter. Its output value
	is the index of this operator in the recording. Thus the output
	value following the previous call to start_reverse is equal to 
	the number of variables in the recording minus one.
	In addition, the output value decreases by one with each call to
	next_reverse.
	The last operator sets op_index equal to 0.

	\param var_index
	The input value of var_index does not matter. Its output value is the
	index of the primary (last) result corresponding to the operator op.
	The last operator sets var_index equal to 0 (corresponding to BeginOp
	at beginning of operation sequence).
	*/

	void next_reverse(
	OpCode& op, const addr_t*& op_arg, size_t& op_index, size_t& var_index)
	{	using CppAD::NumRes;
		using CppAD::NumArg;

		// index of the last result for the next operator
		CPPAD_ASSERT_UNKNOWN( var_index_ >= NumRes(op_) );
		var_index   = var_index_ -= NumRes(op_);

		// next operator
		CPPAD_ASSERT_UNKNOWN( op_index_  > 0 );
		op_index    = --op_index_;                                  // index
		op          = op_         = OpCode( rec_op_[ op_index_ ] ); // value

		// first argument for next operator
		CPPAD_ASSERT_UNKNOWN( op_arg_ >= NumArg(op)  );
		op_arg_    -= NumArg(op);                            // index
		op_arg      = op_arg_ + rec_op_arg_.data();          // pointer
	}
	/*!
	Correct \c next_reverse return values when <tt>op == CSumOp</tt>.

	\param op
	The input value of op must be the return value from the previous
	call to \c next_reverse and must be \c CSumOp.

	\param op_arg
	The input value of *op_arg must be the return value from the 
	previous call to \c next_reverse. Its output value is the
	beginning of the vector of argument indices for this operation.

	\param op_index
	The input value of op_index must be the return value from the
	previous call to \c next_reverse. Its output value
	is the index of the this operator in the recording. 

	\param var_index
	The input value of var_index must be the return value from the 
	previous call to \c next_reverse. Its output value is the
	index of the primary (last) result corresponding to this operator.
	*/

	void reverse_csum(
	OpCode& op, const addr_t*& op_arg, size_t& op_index, size_t& var_index)
	{	using CppAD::NumRes;
		using CppAD::NumArg;
		CPPAD_ASSERT_UNKNOWN( op == CSumOp );
		CPPAD_ASSERT_UNKNOWN( NumArg(CSumOp) == 0 );
		/*
		The things needs fixing are op_arg_ and op_arg. Currently, 
		op_arg points first arugment for the previous operator.
		*/
		--op_arg;
		op_arg_    -= (op_arg[0] + 4);
		op_arg      = op_arg_ + rec_op_arg_.data();

		CPPAD_ASSERT_UNKNOWN(
		op_arg[0] + op_arg[1] == op_arg[ 3 + op_arg[0] + op_arg[1] ]
		);
		CPPAD_ASSERT_UNKNOWN( op_index_  < rec_op_.size() );
		CPPAD_ASSERT_UNKNOWN( op_arg_ + NumArg(op) <= rec_op_arg_.size() );
		CPPAD_ASSERT_UNKNOWN( var_index_  < num_rec_var_ );
	}

};

CPPAD_END_NAMESPACE
# endif

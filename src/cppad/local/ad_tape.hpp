/* $Id: ad_tape.hpp 2038 2011-07-31 23:56:44Z bradbell $ */
# ifndef CPPAD_AD_TAPE_INCLUDED
# define CPPAD_AD_TAPE_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
# include <cppad/local/define.hpp>

CPPAD_BEGIN_NAMESPACE

/*!
Class used to hold tape that records AD<Base> operations.

\tparam Base
An <tt>AD<Base></tt> object is used to recording <tt>AD<Base></tt> operations.
*/

template <class Base>
class ADTape {
	// Friends =============================================================

	// classes -------------------------------------------------------------
	friend class AD<Base>;
	friend class ADFun<Base>;
	friend class discrete<Base>;
	friend class user_atomic<Base>;
	friend class VecAD<Base>;
	friend class VecAD_reference<Base>;

	// functions -----------------------------------------------------------
	// PrintFor
	friend void PrintFor <Base> (
		const AD<Base>&    flag   ,
		const char*        before ,
		const AD<Base>&    var    , 
		const char*        after
	);
	// CondExpOp
	friend AD<Base> CondExpOp <Base> (
		enum CompareOp  cop          ,
		const AD<Base> &left         , 
		const AD<Base> &right        , 
		const AD<Base> &trueCase     , 
		const AD<Base> &falseCase 
	);
	// pow
	friend AD<Base> pow <Base>
		(const AD<Base> &x, const AD<Base> &y);
	// Parameter
	friend bool Parameter     <Base> 
		(const AD<Base> &u);
	// Variable
	friend bool Variable      <Base> 
		(const AD<Base> &u);
	// operators -----------------------------------------------------------
	// arithematic binary operators
	friend AD<Base> operator + <Base>
		(const AD<Base> &left, const AD<Base> &right);
	friend AD<Base> operator - <Base>
		(const AD<Base> &left, const AD<Base> &right);
	friend AD<Base> operator * <Base>
		(const AD<Base> &left, const AD<Base> &right);
	friend AD<Base> operator / <Base>
		(const AD<Base> &left, const AD<Base> &right);

	// comparison operators
	friend bool operator < <Base>
		(const AD<Base> &left, const AD<Base> &right);
	friend bool operator <= <Base>
		(const AD<Base> &left, const AD<Base> &right);
	friend bool operator > <Base>
		(const AD<Base> &left, const AD<Base> &right);
	friend bool operator >= <Base>
		(const AD<Base> &left, const AD<Base> &right);
	friend bool operator == <Base>
		(const AD<Base> &left, const AD<Base> &right);
	friend bool operator != <Base>
		(const AD<Base> &left, const AD<Base> &right);
	// ======================================================================

// --------------------------------------------------------------------------
private:
	// ----------------------------------------------------------------------
	// private data
	/*!
	Unique identifier for this tape.  It is always greater than 
	CPPAD_MAX_NUM_THREADS, and different for every tape (even ones that have 
	been deleted). In addition, id_ % CPPAD_MAX_NUM_THREADS is the thread 
	number for this tape. Set by Independent and effectively const
	*/
	size_t                       id_;
	/// Number of independent variables in this tapes reconding.
	/// Set by Independent and effectively const
	size_t         size_independent_;
	/// This is where the information is recorded.
	recorder<Base>              Rec_;
	// ----------------------------------------------------------------------
	// private functions
	//
	// add a parameter to the tape
	size_t RecordParOp(const Base &x);
	
	// see CondExp.h
	void RecordCondExp(
		enum CompareOp  cop           ,
		AD<Base>       &returnValue   ,
		const AD<Base> &left          ,
		const AD<Base> &right         ,
		const AD<Base> &trueCase      ,
		const AD<Base> &falseCase
	);

	// see Compare.h
	void RecordCompare(
		enum CompareOp  cop      ,
		bool            result   ,
		const AD<Base> &left     ,
		const AD<Base> &right
	);

	// place a VecAD object in the tape
	size_t AddVec(
		size_t                   length,
		const pod_vector<Base>&  data
	);

public:
	// constructor
	ADTape(size_t id) : id_(id)
	{ }

	// destructor
	~ADTape(void)
	{	Rec_.Erase(); }

	// public function only used by CppAD::Independent
	template <typename VectorADBase>
	void Independent(VectorADBase &u);

};
// ---------------------------------------------------------------------------
// Private functions
//

/*!
Place a parameter in the tape.

On rare occations it is necessary to place a parameter in the tape; e.g.,
when it is one of the dependent variabes.

\param z
value of the parameter that we are placing in the tape.

\return 
variable index (for this recording) correpsonding to the parameter.

\par Wish List
All these operates are preformed in \c Rec_, so we should
move this routine from <tt>ADTape<Base></tt> to <tt>recorder<Base></tt>.
*/
template <class Base>
size_t ADTape<Base>::RecordParOp(const Base &z)
{	size_t z_taddr;
	size_t ind;
	CPPAD_ASSERT_UNKNOWN( NumRes(ParOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( NumArg(ParOp) == 1 );
	z_taddr = Rec_.PutOp(ParOp);
	ind     = Rec_.PutPar(z);
	Rec_.PutArg(ind);

	return z_taddr;
}

/*!
Put initialization for a VecAD<Base> object in the tape.

This routine should be called once for each VecAD object when just
before it changes from a parameter to a variable.

\param length
size of the <tt>VecAD<Base></tt> object.

\param data
initial values for the <tt>VecAD<Base></tt> object
(values before it becomes a variable).

\par Wish List
All these operates are preformed in \c Rec_, so we should
move this routine from <tt>ADTape<Base></tt> to <tt>recorder<Base></tt>.
*/
template <class Base>
size_t ADTape<Base>::AddVec(size_t length, const pod_vector<Base>& data)
{	CPPAD_ASSERT_UNKNOWN( length > 0 );
	size_t i;
	size_t value_index;

	// store the length in VecInd
	size_t start = Rec_.PutVecInd(length);

	// store indices of the values in VecInd 
	for(i = 0; i < length; i++)
	{
		value_index = Rec_.PutPar( data[i] );
		Rec_.PutVecInd( value_index );
	}
 
	// return the taddr of the length (where the vector starts)
	return start;
}

CPPAD_END_NAMESPACE

# endif

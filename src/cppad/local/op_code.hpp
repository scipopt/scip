/* $Id: op_code.hpp 2240 2011-12-31 05:33:55Z bradbell $ */
# ifndef CPPAD_OP_CODE_INCLUDED
# define CPPAD_OP_CODE_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
# include <string>
# include <sstream>
# include <iomanip>

# include <cppad/local/define.hpp>
# include <cppad/local/cppad_assert.hpp>

// needed before one can use CPPAD_ASSERT_FIRST_CALL_NOT_PARALLEL
# include <cppad/thread_alloc.hpp>

CPPAD_BEGIN_NAMESPACE
/*! 
\file op_code.hpp
Defines the OpCode enum type and functions related to it.

*/


/*!
Type used to distinguish different AD< \a Base > atomic operations.

Each of the operators ends with the characters Op. Ignoring the Op at the end,
the operators appear in alphabetical order. Binary operation where both
operands have type AD< \a Base > use the following convention for thier endings:
\verbatim
    Ending  Left-Operand  Right-Operand
      pvOp     parameter       variable  
      vpOp      variable      parameter  
      vvOp      variable       variable  
\endverbatim
For example, AddpvOp represents the addition operator where the left
operand is a parameter and the right operand is a variable.
*/
enum OpCode {
	AbsOp,    //  abs(variable)
	AcosOp,   // asin(variable)
	AddpvOp,  //      parameter  + variable
	AddvvOp,  //      variable   + variable
	AsinOp,   // asin(variable)
	AtanOp,   // atan(variable)
	BeginOp,  // used to mark the beginning of the tape
	CExpOp,   // CondExp(cop, left, right, trueCase, falseCase)
	ComOp,    // Compare(cop, result, left, right)
	CosOp,    //  cos(variable)
	CoshOp,   // cosh(variable)
	CSumOp,   // Cummulative summation (has variable number of arguments)
	DisOp,    //  discrete::eval(index, variable)
	DivpvOp,  //      parameter  / variable
	DivvpOp,  //      variable   / parameter
	DivvvOp,  //      variable   / variable
	EndOp,    //  used to mark the end of the tape
	ExpOp,    //  exp(variable)
	InvOp,    //                             independent variable
	LdpOp,    //    z[parameter]
	LdvOp,    //    z[variable]
	LogOp,    //  log(variable)
	MulpvOp,  //      parameter  * variable
	MulvvOp,  //      variable   * variable
	ParOp,    //      parameter
	PowpvOp,  //  pow(parameter,   variable)
	PowvpOp,  //  pow(variable,    parameter)
	PowvvOp,  //  pow(variable,    variable)
	PriOp,    //  PrintFor(text, parameter or variable, parameter or variable)
	SignOp,   // sign(variable)
	SinOp,    //  sin(variable)
	SinhOp,   // sinh(variable)
	SqrtOp,   // sqrt(variable)
	StppOp,   //    z[parameter] = parameter
	StpvOp,   //    z[parameter] = variable
	StvpOp,   //    z[variable]  = parameter
	StvvOp,   //    z[variable]  = variable
	SubpvOp,  //      parameter  - variable
	SubvpOp,  //      variable   - parameter
	SubvvOp,  //      variable   - variable
	TanOp,    //  tan(variable)
	TanhOp,   //  tan(variable)
	// user atomic operation codes (note yet implemented)
	UserOp,   //  start of a user atomic operaiton
	UsrapOp,  //  this user atomic argument is a parameter
	UsravOp,  //  this user atomic argument is a variable
	UsrrpOp,  //  this user atomic result is a parameter
	UsrrvOp   //  this user atomic result is a variable
};

/*!
Table containing number of arguments for the corresponding operator.

The i-th element in this table specifes the number of arguments stored for each
occurance of the operator that is the i-th value in the OpCode enum type.
For example, for the first three OpCode enum values we have
\verbatim
OpCode   j   NumArgTable[j]  Meaning
AbsOp    0                1  index of variable we are taking absolute value of
AcosOp   1                1  index of variable we are taking cosine of
AddpvOp  1                2  indices of parameter and variable we are adding
\endverbatim
Note that the meaning of the arguments depends on the operator.
*/
const size_t NumArgTable[] = {
	1, // AbsOp
	1, // AcosOp
	2, // AddpvOp
	2, // AddvvOp
	1, // AsinOp
	1, // AtanOp
	0, // BeginOp
	6, // CExpOp
	4, // ComOp
	1, // CosOp
	1, // CoshOp
	0, // CSumOp   (actually has a variable number of arguments, not zero)
	2, // DisOp
	2, // DivpvOp
	2, // DivvpOp
	2, // DivvvOp
	0, // EndOp
	1, // ExpOp
	0, // InvOp
	3, // LdpOp
	3, // LdvOp
	1, // LogOp
	2, // MulpvOp
	2, // MulvvOp
	1, // ParOp
	2, // PowpvOp
	2, // PowvpOp
	2, // PowvvOp
	5, // PriOp
	1, // SignOp
	1, // SinOp
	1, // SinhOp
	1, // SqrtOp
	3, // StppOp
	3, // StpvOp
	3, // StvpOp
	3, // StvvOp
	2, // SubpvOp
	2, // SubvpOp
	2, // SubvvOp
	1, // TanOp
	1, // TanhOp
	4, // UserOp
	1, // UsrapOp
	1, // UsravOp
	1, // UsrrpOp
	0  // UsrrvOp
};

/*!
Fetch the number of arguments for a specified operator.

\return
Number of arguments corresponding to the specified operator.

\param op 
Operator for which we are fetching the number of arugments.

- Check that \c UsrrvOp is the last defined operator.
- Check that \c CPPAD_OP_CODE_TYPE can support all the operator codes.
- Check that \c op is a valid operator value.
*/
inline size_t NumArg( OpCode op)
{	CPPAD_ASSERT_FIRST_CALL_NOT_PARALLEL;
# ifndef NDEBUG
	// only do these checks once to save time
	static bool first = true;
	if( first )
	{	CPPAD_ASSERT_UNKNOWN( size_t(UsrrvOp) == 
			sizeof(NumArgTable) / sizeof(NumArgTable[0]) - 1
		);
		CPPAD_ASSERT_UNKNOWN( size_t(UsrrvOp) <=
			std::numeric_limits<CPPAD_OP_CODE_TYPE>::max()
		);
		first = false;
	}
	// do this check every time
	CPPAD_ASSERT_UNKNOWN( size_t(op) <= size_t(UsrrvOp) );
# endif

	return NumArgTable[(size_t) op];
}

/*!
Number of variables resulting from the corresponding operation.

The i-th element in this table specifes the number of varibles for each
occurance of the operator that is the i-th value in the OpCode enum type.
For example, for the first three OpCode enum values we have
\verbatim
OpCode   j   NumResTable[j]  Meaning
AbsOp    0                1  variable that is the result of the absolute value
AcosOp   1                2  acos(x) and sqrt(1-x*x) are required for this op
AddpvOp  1                1  variable that is the result of the addition
\endverbatim
*/
// alphabetical order (ignoring the Op at the end)
const size_t NumResTable[] = {
	1, // AbsOp
	2, // AcosOp
	1, // AddpvOp
	1, // AddvvOp
	2, // AsinOp
	2, // AtanOp
	1, // BeginOp  offsets first variable to have index one (not zero)
	1, // CExpOp
	0, // ComOp
	2, // CosOp
	2, // CoshOp
	1, // CSumOp
	1, // DisOp
	1, // DivpvOp
	1, // DivvpOp
	1, // DivvvOp
	0, // EndOp
	1, // ExpOp
	1, // InvOp
	1, // LdpOp
	1, // LdvOp
	1, // LogOp
	1, // MulpvOp
	1, // MulvvOp
	1, // ParOp
	3, // PowpvOp
	3, // PowvpOp
	3, // PowvvOp
	0, // PriOp
	1, // SignOp
	2, // SinOp
	2, // SinhOp
	1, // SqrtOp
	0, // StppOp
	0, // StpvOp
	0, // StvpOp
	0, // StvvOp
	1, // SubpvOp
	1, // SubvpOp
	1, // SubvvOp
	2, // TanOp
	2, // TanhOp
	0, // UserOp
	0, // UsrapOp
	0, // UsravOp
	0, // UsrrpOp
	1, // UsrrvOp
	0  // Last entry not used: avoids warning by g++ 4.3.2 when pycppad builds
};

/*!
Fetch the number of variables resulting from the specified operation.

\return
number of variables resulting from the specified operator.

\param op 
Operator for which we are fetching the number of result variables.
*/
inline size_t NumRes(OpCode op)
{	// check ensuring conversion to size_t is as expected
	CPPAD_ASSERT_UNKNOWN( size_t(UsrrvOp) == 
		sizeof(NumResTable) / sizeof(NumResTable[0]) - 2
	);
	// this test ensures that all indices are within the table
	CPPAD_ASSERT_UNKNOWN( size_t(op) <= size_t(UsrrvOp) );

	return NumResTable[(size_t) op];
}

/*!
Prints a single field corresponding to an operator.

A specified leader is printed in front of the value
and then the value is left justified in the following width character.

\tparam Type
is the type of the value we are printing.

\param os
is the stream that we are printing to.

\param leader
are characters printed before the value.

\param value
is the value being printed.

\param width
is the number of character to print the value in.
If the value does not fit in the width, the value is replace
by width '*' characters.
*/
template <class Type>
void printOpField(
	std::ostream      &os , 
	const char *   leader ,
	const Type     &value , 
	size_t          width )
{
	std::ostringstream buffer;
	std::string        str;

	// first print the leader
	os << leader;

	// print the value into an internal buffer
	buffer << std::setw(width) << value;
	str = buffer.str();

	// length of the string
	size_t len = str.size();
	if( len > width )
	{	size_t i;
		for(i = 0; i < width-1; i++)
			os << str[i];
		os << "*";
		return;
	}

	// count number of spaces at begining
	size_t nspace = 0; 
	while(str[nspace] == ' ' && nspace < len)
		nspace++;

	// left justify the string
	size_t i = nspace;
	while( i < len )
		os << str[i++];

	i = width - len + nspace;
	while(i--)
		os << " "; 
}

/*!
Prints a single operator, its operands, and the corresponding result values.

\tparam Base
Is the base type for these AD< \a Base > operations.

\tparam Value
Determines the type of the values that we are printing.

\param os
is the output stream that the information is printed on.

\param Rec
Is the entire recording for the tape that this operator is in.

\param i_var
is the index for the variable corresponding to the result of this operation
(ignored if NumRes(op) == 0).

\param op
The operator code (OpCode) for this operation.

\param ind
is the vector of argument indices for this operation
(must have NumArg(op) elements).

\param nfz
is the number of forward sweep calculated values of type Value
that correspond to this operation
(ignored if NumRes(op) == 0).

\param fz
points to the first forward calculated value
that correspond to this operation
(ignored if NumRes(op) == 0).

\param nrz
is the number of reverse sweep calculated values of type Value
that correspond to this operation
(ignored if NumRes(op) == 0).

\param rz
points to the first reverse calculated value
that correspond to this operation
(ignored if NumRes(op) == 0).
*/
template <class Base, class Value>
void printOp(
	std::ostream          &os     , 
	const player<Base>   *Rec     ,
	size_t                 i_var  , 
	OpCode                 op     ,
	const addr_t          *ind    ,
	size_t                 nfz    ,
	const  Value          *fz     ,
	size_t                 nrz    ,
	const  Value          *rz     )
{	size_t i;
	
	CPPAD_ASSERT_KNOWN(
		! thread_alloc::in_parallel() ,
		"cannot print trace of AD operations in parallel mode"
	);
	static const char *CompareOpName[] = 
		{ "Lt", "Le", "Eq", "Ge", "Gt", "Ne" };
	static const char *OpName[] = {
		"Abs"   ,
		"Acos"  ,
		"Addpv" ,
		"Addvv" ,
		"Asin"  ,
		"Atan"  ,
		"Begin" ,
		"CExp"  ,
		"Com"   ,
		"Cos"   ,
		"Cosh"  ,
		"CSum"  ,
		"Dis"   ,
		"Divpv" ,
		"Divvp" ,
		"Divvv" ,
		"End"   ,
		"Exp"   ,
		"Inv"   ,
		"Ldp"   ,
		"Ldv"   ,
		"Log"   ,
		"Mulpv" ,
		"Mulvv" ,
		"Par"   ,
		"Powpv" ,
		"Powvp" ,
		"Powvv" ,
		"Pri"   ,
		"Sign"  ,
		"Sin"   ,
		"Sinh"  ,
		"Sqrt"  ,
		"Stpp"  ,
		"Stpv"  ,
		"Stvp"  ,
		"Stvv"  ,
		"Subpv" ,
		"Subvp" ,
		"Subvv" ,
		"Tan"   ,
		"Tanh"  ,
		"User"  ,
		"Usrap" ,
		"Usrav" ,
		"Usrrp" ,
		"Usrrv"
	};
	CPPAD_ASSERT_UNKNOWN( 
		size_t(UsrrvOp) == sizeof(OpName) / sizeof(OpName[0]) - 1
	);

	// print operator
	printOpField(os,  "i=",      i_var, 5);
	if( op == CExpOp )
	{	printOpField(os, "op=", OpName[op], 4); 
		printOpField(os, "", CompareOpName[ ind[0] ], 3);
	}
	else if( op == ComOp )
	{	printOpField(os, "op=", OpName[op], 3); 
		printOpField(os, "", CompareOpName[ ind[0] ], 4);
	}
	else	printOpField(os, "op=", OpName[op], 7); 

	// print other fields
	size_t ncol = 5;
	switch( op )
	{
		case CSumOp:
		/*
		ind[0] = number of addition variables in summation
		ind[1] = number of subtraction variables in summation
		ind[2] = index of parameter that initializes summation
		ind[3], ... , ind[2+ind[0]] = index for positive variables
		ind[3+ind[0]], ..., ind[2+ind[0]+ind[1]] = negative variables 
		ind[3+ind[0]+ind[1]] = ind[0] + ind[1]
		*/
		CPPAD_ASSERT_UNKNOWN( ind[3+ind[0]+ind[1]] == ind[0]+ind[1] );
		printOpField(os, " pr=", Rec->GetPar(ind[2]), ncol);
		for(i = 0; i < ind[0]; i++)
			 printOpField(os, " +v=", ind[3+i], ncol);
		for(i = 0; i < ind[1]; i++)
			 printOpField(os, " -v=", ind[3+ind[0]+i], ncol);
		break;

		case LdpOp:
		CPPAD_ASSERT_UNKNOWN( NumArg(op) == 3 );
		printOpField(os, "off=", ind[0], ncol);
		printOpField(os, "idx=", ind[1], ncol);
		break;

		case LdvOp:
		CPPAD_ASSERT_UNKNOWN( NumArg(op) == 3 );
		printOpField(os, "off=", ind[0], ncol);
		printOpField(os, "  v=", ind[1], ncol);
		break;

		case StppOp:
		CPPAD_ASSERT_UNKNOWN( NumArg(op) == 3 );
		printOpField(os, "off=", ind[0], ncol);
		printOpField(os, "idx=", ind[1], ncol);
		printOpField(os, " pr=", Rec->GetPar(ind[2]), ncol);
		break;

		case StpvOp:
		CPPAD_ASSERT_UNKNOWN( NumArg(op) == 3 );
		printOpField(os, "off=", ind[0], ncol);
		printOpField(os, "idx=", ind[1], ncol);
		printOpField(os, " vr=", ind[2], ncol);
		break;

		case StvpOp:
		CPPAD_ASSERT_UNKNOWN( NumArg(op) == 3 );
		printOpField(os, "off=", ind[0], ncol);
		printOpField(os, " vl=", ind[1], ncol);
		printOpField(os, " pr=", Rec->GetPar(ind[2]), ncol);
		break;

		case StvvOp:
		CPPAD_ASSERT_UNKNOWN( NumArg(op) == 3 );
		printOpField(os, "off=", ind[0], ncol);
		printOpField(os, " vl=", ind[1], ncol);
		printOpField(os, " vr=", ind[2], ncol);
		break;

		case AddvvOp:
		case DivvvOp:
		case MulvvOp:
		case PowvvOp:
		case SubvvOp:
		CPPAD_ASSERT_UNKNOWN( NumArg(op) == 2 );
		printOpField(os, " vl=", ind[0], ncol);
		printOpField(os, " vr=", ind[1], ncol);
		break;

		case AddpvOp:
		case SubpvOp:
		case MulpvOp:
		case PowpvOp:
		case DivpvOp:
		CPPAD_ASSERT_UNKNOWN( NumArg(op) == 2 );
		printOpField(os, " pl=", Rec->GetPar(ind[0]), ncol);
		printOpField(os, " vr=", ind[1], ncol);
		break;

		case DivvpOp:
		case PowvpOp:
		case SubvpOp:
		CPPAD_ASSERT_UNKNOWN( NumArg(op) == 2 );
		printOpField(os, " vl=", ind[0], ncol);
		printOpField(os, " pr=", Rec->GetPar(ind[1]), ncol);
		break;

		case AbsOp:
		case AcosOp:
		case AsinOp:
		case AtanOp:
		case CosOp:
		case CoshOp:
		case ExpOp:
		case LogOp:
		case SignOp:
		case SinOp:
		case SinhOp:
		case SqrtOp:
		case UsravOp:
		case TanOp:
		case TanhOp:
		CPPAD_ASSERT_UNKNOWN( NumArg(op) == 1 );
		printOpField(os, "  v=", ind[0], ncol);
		break;

		case ParOp:
		case UsrapOp:
		case UsrrpOp:
		CPPAD_ASSERT_UNKNOWN( NumArg(op) == 1 );
		printOpField(os, "  p=", Rec->GetPar(ind[0]), ncol);
		break;

		case UserOp:
		CPPAD_ASSERT_UNKNOWN( NumArg(op) == 4 );
		{	const char* name = user_atomic<Base>::name(ind[0]);
			printOpField(os, " f=",   name, ncol);
			printOpField(os, " i=", ind[1], ncol);
			printOpField(os, " n=", ind[2], ncol);
			printOpField(os, " m=", ind[3], ncol);
		}
		break;

		case PriOp:
		CPPAD_ASSERT_NARG_NRES(op, 5, 0);
		if( ind[0] & 1 )
			printOpField(os, " v=", ind[1], ncol);
		else	printOpField(os, " p=", Rec->GetPar(ind[1]), ncol);
		os << "before=\"" << Rec->GetTxt(ind[2]) << "\"";
		if( ind[0] & 2 )
			printOpField(os, " v=", ind[3], ncol);
		else	printOpField(os, " p=", Rec->GetPar(ind[3]), ncol);
		os << "after=\"" << Rec->GetTxt(ind[4]) << "\"";
		break;

		case BeginOp:
		case EndOp:
		case InvOp:
		case UsrrvOp:
		CPPAD_ASSERT_UNKNOWN( NumArg(op) == 0 );
		break;

		case DisOp:
		CPPAD_ASSERT_UNKNOWN( NumArg(op) == 2 );
		{	const char* name = discrete<Base>::name(ind[0]);
			printOpField(os, " f=", name, ncol);
			printOpField(os, " x=", ind[1], ncol);
		}
		break;
	

		case CExpOp:
		CPPAD_ASSERT_UNKNOWN(ind[1] != 0);
		CPPAD_ASSERT_UNKNOWN( NumArg(op) == 6 );
		if( ind[1] & 1 )
			printOpField(os, " vl=", ind[2], ncol);
		else	printOpField(os, " pl=", Rec->GetPar(ind[2]), ncol);
		if( ind[1] & 2 )
			printOpField(os, " vr=", ind[3], ncol);
		else	printOpField(os, " pr=", Rec->GetPar(ind[3]), ncol);
		if( ind[1] & 4 )
			printOpField(os, " vt=", ind[4], ncol);
		else	printOpField(os, " pt=", Rec->GetPar(ind[4]), ncol);
		if( ind[1] & 8 )
			printOpField(os, " vf=", ind[5], ncol);
		else	printOpField(os, " pf=", Rec->GetPar(ind[5]), ncol);
		break;

		case ComOp:
		CPPAD_ASSERT_UNKNOWN(ind[1] != 0);
		CPPAD_ASSERT_UNKNOWN( NumArg(op) == 4 );
		if( ind[1] & 1 )
			printOpField(os, "res=", 1, ncol);
		else	printOpField(os, "res=", 0, ncol);
		if( ind[1] & 2 )
			printOpField(os, " vl=", ind[2], ncol);
		else	printOpField(os, " pl=", Rec->GetPar(ind[2]), ncol);
		if( ind[1] & 4 )
			printOpField(os, " vr=", ind[3], ncol);
		else	printOpField(os, " pr=", Rec->GetPar(ind[3]), ncol);
		break;

		default:
		CPPAD_ASSERT_UNKNOWN(0);
	}
	size_t k;
	if( NumRes(op) > 0 && (op != BeginOp) )
	{ 
		for(k = 0; k < nfz; k++)
			std::cout << "| fz[" << k << "]=" << fz[k];
		for(k = 0; k < nrz; k++)
			std::cout << "| rz[" << k << "]=" << rz[k];
	}
	std::cout << std::endl;
}

CPPAD_END_NAMESPACE
# endif

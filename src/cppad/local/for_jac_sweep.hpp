/* $Id: for_jac_sweep.hpp 2240 2011-12-31 05:33:55Z bradbell $ */
# ifndef CPPAD_FOR_JAC_SWEEP_INCLUDED
# define CPPAD_FOR_JAC_SWEEP_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
# include <set>
# include <cppad/local/pod_vector.hpp>

CPPAD_BEGIN_NAMESPACE
/*!
\file for_jac_sweep.hpp
Compute Forward mode Jacobian sparsity patterns.
*/

/*!
\def CPPAD_FOR_JAC_SWEEP_TRACE
This value is either zero or one. 
Zero is the normal operational value.
If it is one, a trace of every for_jac_sweep computation is printed.
*/
# define CPPAD_FOR_JAC_SWEEP_TRACE 0

/*!
Given the sparsity pattern for the independent variables,
ForJacSweep computes the sparsity pattern for all the other variables.

\tparam Base
base type for the operator; i.e., this operation sequence was recorded
using AD< \a Base > and computations by this routine are done using type 
\a Base.

\tparam Vector_set
is the type used for vectors of sets. It can be either
\c sparse_pack or \c sparse_set.

\param n
is the number of independent variables on the tape.

\param numvar
is the total number of variables on the tape; i.e.,
\a play->num_rec_var().

\param play
The information stored in \a play
is a recording of the operations corresponding to a function
\f[
	F : {\bf R}^n \rightarrow {\bf R}^m
\f]
where \f$ n \f$ is the number of independent variables
and \f$ m \f$ is the number of dependent variables.
The object \a play is effectly constant.
It is not declared const because while playing back the tape
the object \a play holds information about the currentl location
with in the tape and this changes during playback.

\param var_sparsity
\b Input: For j = 1 , ... , \a n, 
the sparsity pattern for the independent variable with index (j-1)
corresponds to the set with index j in \a var_sparsity.
\n
\n
\b Output: For i = \a n + 1 , ... , \a numvar - 1,
the sparsity pattern for the variable with index i on the tape
corresponds to the set with index i in \a var_sparsity.

\par Checked Assertions:
\li numvar == var_sparsity.n_set()
\li numvar == play->num_rec_var()
*/

template <class Base, class Vector_set>
void ForJacSweep(
	size_t                n            ,
	size_t                numvar       ,
	player<Base>*         play         ,
	Vector_set&           var_sparsity )
{
	OpCode           op;
	size_t         i_op;
	size_t        i_var;

	const addr_t*   arg = 0;

	size_t            i, j, k;

	// check numvar argument
	CPPAD_ASSERT_UNKNOWN( play->num_rec_var()  == numvar );
	CPPAD_ASSERT_UNKNOWN( var_sparsity.n_set() == numvar );

	// length of the parameter vector (used by CppAD assert macros)
	const size_t num_par = play->num_rec_par();

	// cum_sparsity accumulates sparsity pattern a cummulative sum
	size_t limit = var_sparsity.end();

	// vecad_sparsity contains a sparsity pattern from each VecAD object
	// to all the other variables.
	// vecad_ind maps a VecAD index (the beginning of the
	// VecAD object) to its from index in vecad_sparsity
	size_t num_vecad_ind   = play->num_rec_vecad_ind();
	size_t num_vecad_vec   = play->num_rec_vecad_vec();
	Vector_set  vecad_sparsity;
	vecad_sparsity.resize(num_vecad_vec, limit);
	pod_vector<size_t> vecad_ind;
	if( num_vecad_vec > 0 )
	{	size_t length;
		vecad_ind.extend(num_vecad_ind);
		j             = 0;
		for(i = 0; i < num_vecad_vec; i++)
		{	// length of this VecAD
			length   = play->GetVecInd(j);
			// set to proper index for this VecAD
			vecad_ind[j] = i; 
			for(k = 1; k <= length; k++)
				vecad_ind[j+k] = num_vecad_vec; // invalid index
			// start of next VecAD
			j       += length + 1;
		}
		CPPAD_ASSERT_UNKNOWN( j == play->num_rec_vecad_ind() );
	}

	// work space used by UserOp.
	typedef std::set<size_t> size_set;
	const size_t user_q = limit; // maximum element plus one
	size_set::iterator set_itr;  // iterator for a standard set
	size_set::iterator set_end;  // end of iterator sequence
	vector< size_set > user_r;   // sparsity pattern for the argument x
	vector< size_set > user_s;   // sparisty pattern for the result y
	size_t user_index = 0;       // indentifier for this user_atomic operation
	size_t user_id    = 0;       // user identifier for this call to operator
	size_t user_i     = 0;       // index in result vector
	size_t user_j     = 0;       // index in argument vector
	size_t user_m     = 0;       // size of result vector
	size_t user_n     = 0;       // size of arugment vector
	// next expected operator in a UserOp sequence
	enum { user_start, user_arg, user_ret, user_end } user_state = user_start;

# if CPPAD_FOR_JAC_SWEEP_TRACE
	std::cout << std::endl;
	CppAD::vector<bool> z_value(limit);
# endif

	// skip the BeginOp at the beginning of the recording
	play->start_forward(op, arg, i_op, i_var);
	CPPAD_ASSERT_UNKNOWN( op == BeginOp );
	while(op != EndOp)
	{
		// this op
		play->next_forward(op, arg, i_op, i_var);
		CPPAD_ASSERT_UNKNOWN( (i_op > n)  | (op == InvOp) );  
		CPPAD_ASSERT_UNKNOWN( (i_op <= n) | (op != InvOp) );  

		// rest of information depends on the case
		switch( op )
		{
			case AbsOp:
			CPPAD_ASSERT_NARG_NRES(op, 1, 1);
			forward_sparse_jacobian_unary_op(
				i_var, arg[0], var_sparsity
			);
			break;
			// -------------------------------------------------

			case AddvvOp:
			CPPAD_ASSERT_NARG_NRES(op, 2, 1);
			forward_sparse_jacobian_binary_op(
				i_var, arg, var_sparsity
			);
			break;
			// -------------------------------------------------

			case AddpvOp:
			CPPAD_ASSERT_NARG_NRES(op, 2, 1);
			forward_sparse_jacobian_unary_op(
				i_var, arg[1], var_sparsity
			);
			break;
			// -------------------------------------------------

			case AcosOp:
			// sqrt(1 - x * x), acos(x)
			CPPAD_ASSERT_NARG_NRES(op, 1, 2);
			forward_sparse_jacobian_unary_op(
				i_var, arg[0], var_sparsity
			);
			break;
			// -------------------------------------------------

			case AsinOp:
			// sqrt(1 - x * x), asin(x)
			CPPAD_ASSERT_NARG_NRES(op, 1, 2);
			forward_sparse_jacobian_unary_op(
				i_var, arg[0], var_sparsity
			);
			break;
			// -------------------------------------------------

			case AtanOp:
			// 1 + x * x, atan(x)
			CPPAD_ASSERT_NARG_NRES(op, 1, 2);
			forward_sparse_jacobian_unary_op(
				i_var, arg[0], var_sparsity
			);
			break;
			// -------------------------------------------------

			case CSumOp:
			// CSumOp has a variable number of arguments and
			// next_forward thinks it one has one argument.
			// we must inform next_forward of this special case.
			play->forward_csum(op, arg, i_op, i_var);
			forward_sparse_jacobian_csum_op(
				i_var, arg, var_sparsity
			);
			break;
			// -------------------------------------------------

			case CExpOp:
			forward_sparse_jacobian_cond_op(
				i_var, arg, num_par, var_sparsity
			);
			break;
			// ---------------------------------------------------

			case ComOp:
			CPPAD_ASSERT_NARG_NRES(op, 4, 0);
			CPPAD_ASSERT_UNKNOWN( arg[1] > 1 );
			break;
			// --------------------------------------------------

			case CosOp:
			// sin(x), cos(x)
			CPPAD_ASSERT_NARG_NRES(op, 1, 2);
			forward_sparse_jacobian_unary_op(
				i_var, arg[0], var_sparsity
			);
			break;
			// ---------------------------------------------------

			case CoshOp:
			// sinh(x), cosh(x)
			CPPAD_ASSERT_NARG_NRES(op, 1, 2);
			forward_sparse_jacobian_unary_op(
				i_var, arg[0], var_sparsity
			);
			break;
			// -------------------------------------------------

			case DisOp:
			CPPAD_ASSERT_NARG_NRES(op, 2, 1);
			var_sparsity.clear(i_var);
			break;
			// -------------------------------------------------

			case DivvvOp:
			CPPAD_ASSERT_NARG_NRES(op, 2, 1);
			forward_sparse_jacobian_binary_op(
				i_var, arg, var_sparsity
			);
			break;
			// -------------------------------------------------

			case DivpvOp:
			CPPAD_ASSERT_NARG_NRES(op, 2, 1);
			forward_sparse_jacobian_unary_op(
				i_var, arg[1], var_sparsity
			);
			break;
			// -------------------------------------------------

			case DivvpOp:
			CPPAD_ASSERT_NARG_NRES(op, 2, 1);
			forward_sparse_jacobian_unary_op(
				i_var, arg[0], var_sparsity
			);
			break;
			// -------------------------------------------------

			case EndOp:
			CPPAD_ASSERT_NARG_NRES(op, 0, 0);
			break;
			// -------------------------------------------------

			case ExpOp:
			CPPAD_ASSERT_NARG_NRES(op, 1, 1);
			forward_sparse_jacobian_unary_op(
				i_var, arg[0], var_sparsity
			);
			break;
			// -------------------------------------------------

			case InvOp:
			CPPAD_ASSERT_NARG_NRES(op, 0, 1);
			// sparsity pattern is already defined
			break;
			// -------------------------------------------------

			case LdpOp:
			forward_sparse_load_op(
				op,
				i_var,
				arg,
				num_vecad_ind,
				vecad_ind.data(),
				var_sparsity,
				vecad_sparsity
			);
			break;
			// -------------------------------------------------

			case LdvOp:
			forward_sparse_load_op(
				op,
				i_var,
				arg,
				num_vecad_ind,
				vecad_ind.data(),
				var_sparsity,
				vecad_sparsity
			);
			break;
			// -------------------------------------------------

			case LogOp:
			CPPAD_ASSERT_NARG_NRES(op, 1, 1);
			forward_sparse_jacobian_unary_op(
				i_var, arg[0], var_sparsity
			);
			break;
			// -------------------------------------------------

			case MulvvOp:
			CPPAD_ASSERT_NARG_NRES(op, 2, 1);
			forward_sparse_jacobian_binary_op(
				i_var, arg, var_sparsity
			);
			break;
			// -------------------------------------------------

			case MulpvOp:
			CPPAD_ASSERT_NARG_NRES(op, 2, 1);
			forward_sparse_jacobian_unary_op(
				i_var, arg[1], var_sparsity
			);
			break;
			// -------------------------------------------------

			case ParOp:
			CPPAD_ASSERT_NARG_NRES(op, 1, 1);
			var_sparsity.clear(i_var);
			break;
			// -------------------------------------------------

			case PowvpOp:
			CPPAD_ASSERT_NARG_NRES(op, 2, 3);
			forward_sparse_jacobian_unary_op(
				i_var, arg[0], var_sparsity
			);
			break;
			// -------------------------------------------------

			case PowpvOp:
			CPPAD_ASSERT_NARG_NRES(op, 2, 3);
			forward_sparse_jacobian_unary_op(
				i_var, arg[1], var_sparsity
			);
			break;
			// -------------------------------------------------

			case PowvvOp:
			CPPAD_ASSERT_NARG_NRES(op, 2, 3);
			forward_sparse_jacobian_binary_op(
				i_var, arg, var_sparsity
			);
			break;
			// -------------------------------------------------

			case PriOp:
			CPPAD_ASSERT_NARG_NRES(op, 5, 0);
			break;
			// -------------------------------------------------

			case SignOp:
			CPPAD_ASSERT_NARG_NRES(op, 1, 1);
			forward_sparse_jacobian_unary_op(
				i_var, arg[0], var_sparsity
			);
			break;
			// -------------------------------------------------

			case SinOp:
			// cos(x), sin(x)
			CPPAD_ASSERT_NARG_NRES(op, 1, 2);
			forward_sparse_jacobian_unary_op(
				i_var, arg[0], var_sparsity
			);
			break;
			// -------------------------------------------------

			case SinhOp:
			// cosh(x), sinh(x)
			CPPAD_ASSERT_NARG_NRES(op, 1, 2);
			forward_sparse_jacobian_unary_op(
				i_var, arg[0], var_sparsity
			);
			break;
			// -------------------------------------------------

			case SqrtOp:
			CPPAD_ASSERT_NARG_NRES(op, 1, 1);
			forward_sparse_jacobian_unary_op(
				i_var, arg[0], var_sparsity
			);
			break;
			// -------------------------------------------------

			case StppOp:
			CPPAD_ASSERT_NARG_NRES(op, 3, 0);
			// storing a parameter does not affect vector sparsity
			break;
			// -------------------------------------------------

			case StpvOp:
			forward_sparse_store_op(
				op,
				arg,
				num_vecad_ind,
				vecad_ind.data(),
				var_sparsity,
				vecad_sparsity
			);
			break;
			// -------------------------------------------------

			case StvpOp:
			CPPAD_ASSERT_NARG_NRES(op, 3, 0);
			// storing a parameter does not affect vector sparsity
			break;
			// -------------------------------------------------

			case StvvOp:
			forward_sparse_store_op(
				op,
				arg,
				num_vecad_ind,
				vecad_ind.data(),
				var_sparsity,
				vecad_sparsity
			);
			break;
			// -------------------------------------------------

			case SubvvOp:
			CPPAD_ASSERT_NARG_NRES(op, 2, 1);
			forward_sparse_jacobian_binary_op(
				i_var, arg, var_sparsity
			);
			break;
			// -------------------------------------------------

			case SubpvOp:
			CPPAD_ASSERT_NARG_NRES(op, 2, 1);
			forward_sparse_jacobian_unary_op(
				i_var, arg[1], var_sparsity
			);
			break;
			// -------------------------------------------------

			case SubvpOp:
			CPPAD_ASSERT_NARG_NRES(op, 2, 1);
			forward_sparse_jacobian_unary_op(
				i_var, arg[0], var_sparsity
			);
			break;
			// -------------------------------------------------

			case TanOp:
			// tan(x)^2, tan(x)
			CPPAD_ASSERT_NARG_NRES(op, 1, 2);
			forward_sparse_jacobian_unary_op(
				i_var, arg[0], var_sparsity
			);
			break;
			// -------------------------------------------------

			case TanhOp:
			// tanh(x)^2, tanh(x)
			CPPAD_ASSERT_NARG_NRES(op, 1, 2);
			forward_sparse_jacobian_unary_op(
				i_var, arg[0], var_sparsity
			);
			break;
			// -------------------------------------------------

			case UserOp:
			// start or end an atomic operation sequence
			CPPAD_ASSERT_UNKNOWN( NumRes( UserOp ) == 0 );
			CPPAD_ASSERT_UNKNOWN( NumArg( UserOp ) == 4 );
			if( user_state == user_start )
			{	user_index = arg[0];
				user_id    = arg[1];
				user_n     = arg[2];
				user_m     = arg[3];
				if(user_r.size() < user_n )
					user_r.resize(user_n);
				if(user_s.size() < user_m )
					user_s.resize(user_m);
				user_j     = 0;
				user_i     = 0;
				user_state = user_arg;
			}
			else
			{	CPPAD_ASSERT_UNKNOWN( user_state == user_end );
				CPPAD_ASSERT_UNKNOWN( user_index == size_t(arg[0]) );
				CPPAD_ASSERT_UNKNOWN( user_id    == size_t(arg[1]) );
				CPPAD_ASSERT_UNKNOWN( user_n     == size_t(arg[2]) );
				CPPAD_ASSERT_UNKNOWN( user_m     == size_t(arg[3]) );
				user_state = user_start;
			}
			break;

			case UsrapOp:
			// parameter argument in an atomic operation sequence
			CPPAD_ASSERT_UNKNOWN( user_state == user_arg );
			CPPAD_ASSERT_UNKNOWN( user_j < user_n );
			CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < num_par );
			// parameters have an empty sparsity pattern
			user_r[user_j].clear();
			++user_j;
			if( user_j == user_n )
			{	// call users function for this operation
				user_atomic<Base>::for_jac_sparse(user_index, user_id,
					user_n, user_m, user_q, user_r, user_s
				);
				user_state = user_ret;
			}
			break;

			case UsravOp:
			// variable argument in an atomic operation sequence
			CPPAD_ASSERT_UNKNOWN( user_state == user_arg );
			CPPAD_ASSERT_UNKNOWN( user_j < user_n );
			CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) <= i_var );
			// set user_r[user_j] to sparsity pattern for variable arg[0]
			user_r[user_j].clear();
			var_sparsity.begin(arg[0]);
			i = var_sparsity.next_element();
			while( i < user_q )
			{	user_r[user_j].insert(i);
				i = var_sparsity.next_element();
			}
			++user_j;
			if( user_j == user_n )
			{	// call users function for this operation
				user_atomic<Base>::for_jac_sparse(user_index, user_id,
					user_n, user_m, user_q, user_r, user_s
				);
				user_state = user_ret;
			}
			break;

			case UsrrpOp:
			// parameter result in an atomic operation sequence
			CPPAD_ASSERT_UNKNOWN( user_state == user_ret );
			CPPAD_ASSERT_UNKNOWN( user_i < user_m );
			user_i++;
			if( user_i == user_m )
				user_state = user_end;
			break;

			case UsrrvOp:
			// variable result in an atomic operation sequence
			CPPAD_ASSERT_UNKNOWN( user_state == user_ret );
			CPPAD_ASSERT_UNKNOWN( user_i < user_m );
			// It might be faster if we add set union to var_sparsity
			// where one of the sets is not in var_sparsity
			set_itr = user_s[user_i].begin();
			set_end = user_s[user_i].end();
			while( set_itr != set_end )
				var_sparsity.add_element(i_var, *set_itr++);
			user_i++;
			if( user_i == user_m )
				user_state = user_end;
			break;
			// -------------------------------------------------

			default:
			CPPAD_ASSERT_UNKNOWN(0);
		}
# if CPPAD_FOR_JAC_SWEEP_TRACE
		// value for this variable
		for(j = 0; j < limit; j++)
			z_value[j] = false;
		var_sparsity.begin(i_var);
		j = var_sparsity.next_element();
		while( j < limit )
		{	z_value[j] = true;
			j = var_sparsity.next_element();
		}
		printOp(
			std::cout,
			play,
			i_var,
			op,
			arg,
			1,
			&z_value,
			0,
			(CppAD::vector<bool> *) CPPAD_NULL
		);
	}
	std::cout << std::endl;
# else
	}
# endif
	CPPAD_ASSERT_UNKNOWN( i_var + 1 == play->num_rec_var() );

	return;
}

CPPAD_END_NAMESPACE

// preprocessor symbols that are local to this file
# undef CPPAD_FOR_JAC_SWEEP_TRACE

# endif

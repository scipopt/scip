/* $Id: reverse_sweep.hpp 2240 2011-12-31 05:33:55Z bradbell $ */
# ifndef CPPAD_REVERSE_SWEEP_INCLUDED
# define CPPAD_REVERSE_SWEEP_INCLUDED

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
\file reverse_sweep.hpp
Compute derivatives of arbitrary order Taylor coefficients.
*/

/*!
\def CPPAD_REVERSE_SWEEP_TRACE
This value is either zero or one. 
Zero is the normal operational value.
If it is one, a trace of every reverse_sweep computation is printed.
*/
# define CPPAD_REVERSE_SWEEP_TRACE 0

/*!
Compute derivative of arbitrary order forward mode Taylor coefficients.

\tparam Base
base type for the operator; i.e., this operation sequence was recorded
using AD< \a Base > and computations by this routine are done using type 
\a Base.

\param d
is the highest order Taylor coefficients that 
we are computing the derivative of.

\param n
is the number of independent variables on the tape.

\param numvar
is the total number of variables on the tape.
This is also equal to the number of rows in the matrix \a Taylor; i.e.,
Rec->num_rec_var().

\param Rec
The information stored in \a Rec
is a recording of the operations corresponding to the function
\f[
	F : {\bf R}^n \rightarrow {\bf R}^m
\f]
where \f$ n \f$ is the number of independent variables and
\f$ m \f$ is the number of dependent variables.
We define the function 
\f$ G : {\bf R}^{n \times d} \rightarrow {\bf R} \f$ by
\f[
G( u ) = \frac{1}{d !} \frac{ \partial^d }{ \partial t^d } 
\left[ 
	\sum_{i=1}^m w_i  F_i ( u^{(0)} + u^{(1)} t + \cdots + u^{(d)} t^d )
\right]_{t=0}
\f]
Note that the scale factor  1 / a d  converts 
the \a d-th partial derivative to the \a d-th order Taylor coefficient.
This routine computes the derivative of \f$ G(u) \f$
with respect to all the Taylor coefficients
\f$ u^{(k)} \f$ for \f$ k = 0 , ... , d \f$.
The vector \f$ w \in {\bf R}^m \f$, and
value of \f$ u \in {\bf R}^{n \times d} \f$
at which the derivative is computed,
are defined below.
\n
\n
The object \a Rec is effectly constant.
There is an exception to this,
while palying back the tape
the object \a Rec holds information about the current location
with in the tape and this changes during palyback. 

\param J
Is the number of columns in the coefficient matrix \a Taylor.
This must be greater than or equal \a d + 1.

\param Taylor
For i = 1 , ... , \a numvar, and for k = 0 , ... , \a d,
\a Taylor [ i * J + k ]
is the k-th order Taylor coefficient corresponding to 
variable with index i on the tape.
The value \f$ u \in {\bf R}^{n \times d} \f$,
at which the derivative is computed,
is defined by
\f$ u_j^{(k)} \f$ = \a Taylor [ j * J + k ]
for j = 1 , ... , \a n, and for k = 0 , ... , \a d.

\param K
Is the number of columns in the partial derivative matrix \a Partial.
It must be greater than or equal \a d + 1.

\param Partial
\b Input:
The last \f$ m \f$ rows of \a Partial are inputs.
The vector \f$ v \f$, used to define \f$ G(u) \f$,
is specified by these rows. 
For i = 0 , ... , m - 1, \a Partial [ ( \a numvar - m + i ) * K + d ] = v_i.
For i = 0 , ... , m - 1 and for k = 0 , ... , d - 1, 
\a Partial [ ( \a numvar - m + i ) * K + k ] = 0.
\n
\n
\b Temporary:
For i = n+1 , ... , \a numvar - 1 and for k = 0 , ... , d, 
the value of \a Partial [ i * K + k ] is used for temporary work space
and its output value is not defined. 
\n
\n
\b Output:
For j = 1 , ... , n and for k = 0 , ... , d, 
\a Partial [ j * K + k ] 
is the partial derivative of \f$ G( u ) \f$ with 
respect to \f$ u_j^{(k)} \f$.

\par Assumptions
The first operator on the tape is a BeginOp,
and the next \a n operators are InvOp operations for the 
corresponding independent variables.
*/
template <class Base>
void ReverseSweep(
	size_t                d,
	size_t                n,
	size_t                numvar,
	player<Base>*         Rec,
	size_t                J,
	const Base*           Taylor,
	size_t                K,
	Base*                 Partial
)
{
	OpCode           op;
	size_t         i_op;
	size_t        i_var;

	const addr_t*   arg = 0;

	// check numvar argument
	CPPAD_ASSERT_UNKNOWN( Rec->num_rec_var() == numvar );
	CPPAD_ASSERT_UNKNOWN( numvar > 0 );

	// length of the parameter vector (used by CppAD assert macros)
	const size_t num_par = Rec->num_rec_par();

	// pointer to the beginning of the parameter vector
	const Base* parameter = 0;
	if( num_par > 0 )
		parameter = Rec->GetPar();

	// work space used by UserOp.
	const size_t user_k  = d;    // order of this forward mode calculation
	const size_t user_k1 = d+1;  // number of orders for this calculation
	vector<size_t> user_ix;      // variable indices for argument vector
	vector<Base> user_tx;        // argument vector Taylor coefficients
	vector<Base> user_ty;        // result vector Taylor coefficients
	vector<Base> user_px;        // partials w.r.t argument vector
	vector<Base> user_py;        // partials w.r.t. result vector
	size_t user_index = 0;       // indentifier for this user_atomic operation
	size_t user_id    = 0;       // user identifier for this call to operator
	size_t user_i     = 0;       // index in result vector
	size_t user_j     = 0;       // index in argument vector
	size_t user_m     = 0;       // size of result vector
	size_t user_n     = 0;       // size of arugment vector
	// next expected operator in a UserOp sequence
	enum { user_start, user_arg, user_ret, user_end } user_state = user_end;

	// temporary indices
	size_t j, ell;

	// Initialize
	Rec->start_reverse(op, arg, i_op, i_var);
	CPPAD_ASSERT_UNKNOWN( op == EndOp );
# if CPPAD_REVERSE_SWEEP_TRACE
	std::cout << std::endl;
# endif
	while(op != BeginOp )
	{	// next op
		Rec->next_reverse(op, arg, i_op, i_var);
# ifndef NDEBUG
		if( i_op <= n )
		{	CPPAD_ASSERT_UNKNOWN((op == InvOp) | (op == BeginOp));
		}
		else	CPPAD_ASSERT_UNKNOWN((op != InvOp) & (op != BeginOp));
# endif

		// rest of informaiton depends on the case
# if CPPAD_REVERSE_SWEEP_TRACE
		size_t       i_tmp  = i_var;
		const Base*  Z_tmp  = Taylor + i_var * J;
		const Base*  pZ_tmp = Partial + i_var * K;

		printOp(
			std::cout, 
			Rec,
			i_tmp,
			op, 
			arg,
			d + 1, 
			Z_tmp, 
			d + 1, 
			pZ_tmp 
		);
# endif

		switch( op )
		{

			case AbsOp:
			reverse_abs_op(
				d, i_var, arg[0], J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------

			case AddvvOp:
			reverse_addvv_op(
				d, i_var, arg, parameter, J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------

			case AddpvOp:
			CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < num_par );
			reverse_addpv_op(
				d, i_var, arg, parameter, J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------

			case AcosOp:
                        // sqrt(1 - x * x), acos(x)
			CPPAD_ASSERT_UNKNOWN( i_var < numvar );
			reverse_acos_op(
				d, i_var, arg[0], J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------

			case AsinOp:
                        // sqrt(1 - x * x), asin(x)
			CPPAD_ASSERT_UNKNOWN( i_var < numvar );
			reverse_asin_op(
				d, i_var, arg[0], J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------

			case AtanOp:
                        // 1 + x * x, atan(x)
			CPPAD_ASSERT_UNKNOWN( i_var < numvar );
			reverse_atan_op(
				d, i_var, arg[0], J, Taylor, K, Partial
			);
			break;
			// -------------------------------------------------

			case BeginOp:
			CPPAD_ASSERT_NARG_NRES(op, 0, 1);
			break;
			// --------------------------------------------------

			case CSumOp:
			// CSumOp has a variable number of arguments and
			// next_reverse thinks it one has one argument.
			// We must inform next_reverse of this special case.
			Rec->reverse_csum(op, arg, i_op, i_var);
			reverse_csum_op(
				d, i_var, arg, K, Partial
			);
			// end of a cummulative summation
			break;
			// -------------------------------------------------

			case CExpOp:
			reverse_cond_op(
				d, 
				i_var, 
				arg, 
				num_par, 
				parameter, 
				J, 
				Taylor,
				K, 
				Partial
			);
			break;
			// --------------------------------------------------

			case ComOp:
			break;
			// --------------------------------------------------

			case CosOp:
			CPPAD_ASSERT_UNKNOWN( i_var < numvar );
			reverse_cos_op(
				d, i_var, arg[0], J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------

			case CoshOp:
			CPPAD_ASSERT_UNKNOWN( i_var < numvar );
			reverse_cosh_op(
				d, i_var, arg[0], J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------

			case DisOp:
			// Derivative of discrete operation is zero so no
			// contribution passes through this operation. 
			break;
			// --------------------------------------------------

			case DivvvOp:
			reverse_divvv_op(
				d, i_var, arg, parameter, J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------

			case DivpvOp:
			CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < num_par );
			reverse_divpv_op(
				d, i_var, arg, parameter, J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------

			case DivvpOp:
			CPPAD_ASSERT_UNKNOWN( size_t(arg[1]) < num_par );
			reverse_divvp_op(
				d, i_var, arg, parameter, J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------

			case ExpOp:
			reverse_exp_op(
				d, i_var, arg[0], J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------
			case LdpOp:
			reverse_load_op(
				op, d, i_var, arg, J, Taylor, K, Partial
			);
			break;
			// -------------------------------------------------

			case LdvOp:
			reverse_load_op(
				op, d, i_var, arg, J, Taylor, K, Partial
			);
			break;
			// -------------------------------------------------

			case InvOp:
			break;
			// --------------------------------------------------

			case LogOp:
			reverse_log_op(
				d, i_var, arg[0], J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------

			case MulvvOp:
			reverse_mulvv_op(
				d, i_var, arg, parameter, J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------

			case MulpvOp:
			CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < num_par );
			reverse_mulpv_op(
				d, i_var, arg, parameter, J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------

			case ParOp:
			break;
			// --------------------------------------------------

			case PowvpOp:
			CPPAD_ASSERT_UNKNOWN( size_t(arg[1]) < num_par );
			reverse_powvp_op(
				d, i_var, arg, parameter, J, Taylor, K, Partial
			);
			break;
			// -------------------------------------------------

			case PowpvOp:
			CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < num_par );
			reverse_powpv_op(
				d, i_var, arg, parameter, J, Taylor, K, Partial
			);
			break;
			// -------------------------------------------------

			case PowvvOp:
			reverse_powvv_op(
				d, i_var, arg, parameter, J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------

			case PriOp:
			// no result so nothing to do
			break;
			// --------------------------------------------------

			case SignOp:
			CPPAD_ASSERT_UNKNOWN( i_var < numvar );
			reverse_sign_op(
				d, i_var, arg[0], J, Taylor, K, Partial
			);
			break;
			// -------------------------------------------------

			case SinOp:
			CPPAD_ASSERT_UNKNOWN( i_var < numvar );
			reverse_sin_op(
				d, i_var, arg[0], J, Taylor, K, Partial
			);
			break;
			// -------------------------------------------------

			case SinhOp:
			CPPAD_ASSERT_UNKNOWN( i_var < numvar );
			reverse_sinh_op(
				d, i_var, arg[0], J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------

			case SqrtOp:
			reverse_sqrt_op(
				d, i_var, arg[0], J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------

			case StppOp:
			break;
			// --------------------------------------------------

			case StpvOp:
			break;
			// -------------------------------------------------

			case StvpOp:
			break;
			// -------------------------------------------------

			case StvvOp:
			break;
			// --------------------------------------------------

			case SubvvOp:
			reverse_subvv_op(
				d, i_var, arg, parameter, J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------

			case SubpvOp:
			CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < num_par );
			reverse_subpv_op(
				d, i_var, arg, parameter, J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------

			case SubvpOp:
			CPPAD_ASSERT_UNKNOWN( size_t(arg[1]) < num_par );
			reverse_subvp_op(
				d, i_var, arg, parameter, J, Taylor, K, Partial
			);
			break;
			// -------------------------------------------------

			case TanOp:
			CPPAD_ASSERT_UNKNOWN( i_var < numvar );
			reverse_tan_op(
				d, i_var, arg[0], J, Taylor, K, Partial
			);
			break;
			// -------------------------------------------------

			case TanhOp:
			CPPAD_ASSERT_UNKNOWN( i_var < numvar );
			reverse_tanh_op(
				d, i_var, arg[0], J, Taylor, K, Partial
			);
			break;
			// --------------------------------------------------

			case UserOp:
			// start or end an atomic operation sequence
			CPPAD_ASSERT_UNKNOWN( NumRes( UserOp ) == 0 );
			CPPAD_ASSERT_UNKNOWN( NumArg( UserOp ) == 4 );
			if( user_state == user_end )
			{	user_index = arg[0];
				user_id    = arg[1];
				user_n     = arg[2];
				user_m     = arg[3];
				if(user_ix.size() < user_n)
					user_ix.resize(user_n);
				if(user_tx.size() < user_n * user_k1)
				{	user_tx.resize(user_n * user_k1);
					user_px.resize(user_n * user_k1);
				}
				if(user_ty.size() < user_m * user_k1)
				{	user_ty.resize(user_m * user_k1);
					user_py.resize(user_m * user_k1);
				}
				user_j     = user_n;
				user_i     = user_m;
				user_state = user_ret;
			}
			else
			{	CPPAD_ASSERT_UNKNOWN( user_state == user_start );
				CPPAD_ASSERT_UNKNOWN( user_index == size_t(arg[0]) );
				CPPAD_ASSERT_UNKNOWN( user_id    == size_t(arg[1]) );
				CPPAD_ASSERT_UNKNOWN( user_n     == size_t(arg[2]) );
				CPPAD_ASSERT_UNKNOWN( user_m     == size_t(arg[3]) );
				user_state = user_end;

				// call users function for this operation
				user_atomic<Base>::reverse(user_index, user_id,
					user_k, user_n, user_m, user_tx, user_ty,
					user_px, user_py
				);
				for(j = 0; j < user_n; j++) if( user_ix[j] > 0 )
				{	for(ell = 0; ell < user_k1; ell++)
						Partial[user_ix[j] * K + ell] +=
							user_px[j * user_k1 + ell];
				}
			}
			break;

			case UsrapOp:
			// parameter argument in an atomic operation sequence
			CPPAD_ASSERT_UNKNOWN( user_state == user_arg );
			CPPAD_ASSERT_UNKNOWN( 0 < user_j && user_j <= user_n );
			CPPAD_ASSERT_UNKNOWN( NumArg(op) == 1 );
			CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < num_par );
			--user_j;
			user_ix[user_j] = 0;
			user_tx[user_j * user_k1 + 0] = parameter[ arg[0]];
			for(ell = 1; ell < user_k1; ell++)
				user_tx[user_j * user_k1 + ell] = Base(0.);
			
			if( user_j == 0 )
				user_state = user_start;
			break;

			case UsravOp:
			// variable argument in an atomic operation sequence
			CPPAD_ASSERT_UNKNOWN( user_state == user_arg );
			CPPAD_ASSERT_UNKNOWN( 0 < user_j && user_j <= user_n );
			CPPAD_ASSERT_UNKNOWN( NumArg(op) == 1 );
			CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) <= i_var );
			CPPAD_ASSERT_UNKNOWN( 0 < arg[0] );
			--user_j;
			user_ix[user_j] = arg[0];
			for(ell = 0; ell < user_k1; ell++)
				user_tx[user_j*user_k1 + ell] = Taylor[ arg[0] * J + ell];
			if( user_j == 0 )
				user_state = user_start;
			break;

			case UsrrpOp:
			// parameter result in an atomic operation sequence
			CPPAD_ASSERT_UNKNOWN( user_state == user_ret );
			CPPAD_ASSERT_UNKNOWN( 0 < user_i && user_i <= user_m );
			CPPAD_ASSERT_UNKNOWN( NumArg(op) == 1 );
			CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < num_par );
			--user_i;
			for(ell = 0; ell < user_k1; ell++)
			{	user_py[user_i * user_k1 + ell] = Base(0.);
				user_ty[user_i * user_k1 + ell] = Base(0.);
			}
			user_ty[user_i * user_k1 + 0] = parameter[ arg[0] ];
			if( user_i == 0 )
				user_state = user_arg;
			break;

			case UsrrvOp:
			// variable result in an atomic operation sequence
			CPPAD_ASSERT_UNKNOWN( user_state == user_ret );
			CPPAD_ASSERT_UNKNOWN( 0 < user_i && user_i <= user_m );
			--user_i;
			for(ell = 0; ell < user_k1; ell++)
			{	user_py[user_i * user_k1 + ell] =
						Partial[i_var * K + ell];
				user_ty[user_i * user_k1 + ell] =
						Taylor[i_var * J + ell];
			}
			if( user_i == 0 )
				user_state = user_arg;
			break;
			// ------------------------------------------------------------

			default:
			CPPAD_ASSERT_UNKNOWN(false);
		}
	}
# if CPPAD_REVERSE_SWEEP_TRACE
	std::cout << std::endl;
# endif
	// values corresponding to BeginOp
	CPPAD_ASSERT_UNKNOWN( i_op == 0 );
	CPPAD_ASSERT_UNKNOWN( i_var == 0 );
}

CPPAD_END_NAMESPACE

// preprocessor symbols that are local to this file
# undef CPPAD_REVERSE_SWEEP_TRACE

# endif

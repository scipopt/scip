/* $Id: store_op.hpp 1986 2011-06-18 20:33:17Z bradbell $ */
# ifndef CPPAD_STORE_OP_INCLUDED
# define CPPAD_STORE_OP_INCLUDED

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
\file store_op.hpp
Changing the current value of a VecAD element.
*/

/*!
Zero order forward mode implementation of op = StppOp.

\copydetails forward_store_op_0
*/
template <class Base>
inline void forward_store_pp_op_0(
	size_t         i_z         ,
	const addr_t*  arg         , 
	size_t         num_par     ,
	size_t         nc_taylor   ,
	Base*          taylor      ,
	size_t         nc_combined ,
	bool*          variable    ,
	size_t*        combined    )
{	size_t i_vec = arg[1];

	// Because the index is a parameter, this indexing error should be 
	// caught and reported to the user when the tape is recording.
	CPPAD_ASSERT_UNKNOWN( i_vec < combined[ arg[0] - 1 ] );

	CPPAD_ASSERT_UNKNOWN( variable != CPPAD_NULL );
	CPPAD_ASSERT_UNKNOWN( combined != CPPAD_NULL );
	CPPAD_ASSERT_UNKNOWN( NumArg(StppOp) == 3 );
	CPPAD_ASSERT_UNKNOWN( NumRes(StppOp) == 0 );
	CPPAD_ASSERT_UNKNOWN( 0 < arg[0] );
	CPPAD_ASSERT_UNKNOWN( arg[0] + i_vec < nc_combined );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[2]) < num_par );

	variable[ arg[0] + i_vec ] = false;
	combined[ arg[0] + i_vec ] = arg[2];
}

/*!
Zero order forward mode implementation of op = StpvOp.

\copydetails forward_store_op_0
*/
template <class Base>
inline void forward_store_pv_op_0(
	size_t         i_z         ,
	const addr_t*  arg         , 
	size_t         num_par     ,
	size_t         nc_taylor   ,
	Base*          taylor      ,
	size_t         nc_combined ,
	bool*          variable    ,
	size_t*        combined    )
{	size_t i_vec = arg[1];

	// Because the index is a parameter, this indexing error should be 
	// caught and reported to the user when the tape is recording.
	CPPAD_ASSERT_UNKNOWN( i_vec < combined[ arg[0] - 1 ] );

	CPPAD_ASSERT_UNKNOWN( variable != CPPAD_NULL );
	CPPAD_ASSERT_UNKNOWN( combined != CPPAD_NULL );
	CPPAD_ASSERT_UNKNOWN( NumArg(StpvOp) == 3 );
	CPPAD_ASSERT_UNKNOWN( NumRes(StpvOp) == 0 );
	CPPAD_ASSERT_UNKNOWN( 0 < arg[0] );
	CPPAD_ASSERT_UNKNOWN( arg[0] + i_vec < nc_combined );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[2]) <= i_z );

	variable[ arg[0] + i_vec ] = true;
	combined[ arg[0] + i_vec ] = arg[2];
}

/*!
Zero order forward mode implementation of op = StvpOp.

\copydetails forward_store_op_0
*/
template <class Base>
inline void forward_store_vp_op_0(
	size_t         i_z         ,
	const addr_t*  arg         , 
	size_t         num_par     ,
	size_t         nc_taylor   ,
	Base*          taylor      ,
	size_t         nc_combined ,
	bool*          variable    ,
	size_t*        combined    )
{	
	CPPAD_ASSERT_UNKNOWN( size_t(arg[1]) <= i_z );
	size_t i_vec = Integer( taylor[ arg[1] * nc_taylor + 0 ] );
	CPPAD_ASSERT_KNOWN( 
		i_vec < combined[ arg[0] - 1 ] ,
		"VecAD: index during zero order forward sweep is out of range"
	);

	CPPAD_ASSERT_UNKNOWN( variable != CPPAD_NULL );
	CPPAD_ASSERT_UNKNOWN( combined != CPPAD_NULL );
	CPPAD_ASSERT_UNKNOWN( NumArg(StvpOp) == 3 );
	CPPAD_ASSERT_UNKNOWN( NumRes(StvpOp) == 0 );
	CPPAD_ASSERT_UNKNOWN( 0 < arg[0] );
	CPPAD_ASSERT_UNKNOWN( arg[0] + i_vec < nc_combined );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[2]) < num_par );

	variable[ arg[0] + i_vec ] = false;
	combined[ arg[0] + i_vec ] = arg[2];
}

/*!
Zero order forward mode implementation of op = StvvOp.

\copydetails forward_store_op_0
*/
template <class Base>
inline void forward_store_vv_op_0(
	size_t         i_z         ,
	const addr_t*  arg         , 
	size_t         num_par     ,
	size_t         nc_taylor   ,
	Base*          taylor      ,
	size_t         nc_combined ,
	bool*          variable    ,
	size_t*        combined    )
{	
	CPPAD_ASSERT_UNKNOWN( size_t(arg[1]) <= i_z );
	size_t i_vec = Integer( taylor[ arg[1] * nc_taylor + 0 ] );
	CPPAD_ASSERT_KNOWN( 
		i_vec < combined[ arg[0] - 1 ] ,
		"VecAD: index during zero order forward sweep is out of range"
	);

	CPPAD_ASSERT_UNKNOWN( variable != CPPAD_NULL );
	CPPAD_ASSERT_UNKNOWN( combined != CPPAD_NULL );
	CPPAD_ASSERT_UNKNOWN( NumArg(StvpOp) == 3 );
	CPPAD_ASSERT_UNKNOWN( NumRes(StvpOp) == 0 );
	CPPAD_ASSERT_UNKNOWN( 0 < arg[0] );
	CPPAD_ASSERT_UNKNOWN( arg[0] + i_vec < nc_combined );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[2]) <= i_z );

	variable[ arg[0] + i_vec ] = true;
	combined[ arg[0] + i_vec ] = arg[2];
}

/*!
Forward mode sparsity operations for StpvOp and StvvOp

\copydetails sparse_store_op
*/
template <class Vector_set>
inline void forward_sparse_store_op(
	OpCode              op             ,
	const addr_t*       arg            , 
	size_t              num_combined   ,
	const size_t*       combined       ,
	Vector_set&         var_sparsity   ,
	Vector_set&         vecad_sparsity )
{
	CPPAD_ASSERT_UNKNOWN( NumArg(op) == 3 );
	CPPAD_ASSERT_UNKNOWN( NumRes(op) == 0 );
	CPPAD_ASSERT_UNKNOWN( 0 < arg[0] );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < num_combined );
	size_t i_v = combined[ arg[0] - 1 ];
	CPPAD_ASSERT_UNKNOWN( i_v < vecad_sparsity.n_set() );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[2]) < var_sparsity.n_set() );

	vecad_sparsity.binary_union(i_v, i_v, arg[2], var_sparsity);

	return;
}

/*!
Reverse mode sparsity operations for StpvOp and StvvOp

This routine is given the sparsity patterns for
G(v[x], y , w , u ... ) and it uses them to compute the 
sparsity patterns for  
\verbatim
	H(y , w , u , ... ) = G[ v[x], y , w , u , ... ]
\endverbatim

\copydetails sparse_store_op
*/
template <class Vector_set>
inline void reverse_sparse_jacobian_store_op(
	OpCode             op              ,
	const addr_t*      arg             , 
	size_t             num_combined    ,
	const size_t*      combined        ,
	Vector_set&        var_sparsity    ,
	Vector_set&        vecad_sparsity  )
{
	CPPAD_ASSERT_UNKNOWN( NumArg(op) == 3 );
	CPPAD_ASSERT_UNKNOWN( NumRes(op) == 0 );
	CPPAD_ASSERT_UNKNOWN( 0 < arg[0] );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < num_combined );
	size_t i_v = combined[ arg[0] - 1 ];
	CPPAD_ASSERT_UNKNOWN( i_v < vecad_sparsity.n_set() );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[2]) < var_sparsity.n_set() );

	var_sparsity.binary_union(arg[2], arg[2], i_v, vecad_sparsity);

	return;
}

/*!
Reverse mode sparsity operations for StpvOp and StvvOp

This routine is given the sparsity patterns for
G(v[x], y , w , u ... )
and it uses them to compute the sparsity patterns for 
\verbatim
	H(y , w , u , ... ) = G[ v[x], y , w , u , ... ]
\endverbatim

\copydetails sparse_store_op

\param var_jacobian
\a var_jacobian[ \a arg[2] ] 
is false (true) if the Jacobian of G with respect to y is always zero 
(may be non-zero).

\param vecad_jacobian
\a vecad_jacobian[i_v] 
is false (true) if the Jacobian with respect to x is always zero 
(may be non-zero).
On input, it corresponds to the function G,
and on output it corresponds to the function H.
*/
template <class Vector_set>
inline void reverse_sparse_hessian_store_op(
	OpCode             op           ,
	const addr_t*      arg          , 
	size_t             num_combined ,
	const size_t*      combined     ,
	Vector_set&        var_sparsity ,
	Vector_set&        vecad_sparsity ,
	bool*              var_jacobian   ,
	bool*              vecad_jacobian )
{
	CPPAD_ASSERT_UNKNOWN( NumArg(op) == 3 );
	CPPAD_ASSERT_UNKNOWN( NumRes(op) == 0 );
	CPPAD_ASSERT_UNKNOWN( 0 < arg[0] );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < num_combined );
	size_t i_v = combined[ arg[0] - 1 ];
	CPPAD_ASSERT_UNKNOWN( i_v < vecad_sparsity.n_set() );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[2]) < var_sparsity.n_set() );

	var_sparsity.binary_union(arg[2], arg[2], i_v, vecad_sparsity);

	var_jacobian[ arg[2] ] |= vecad_jacobian[i_v];

	return;
}

CPPAD_END_NAMESPACE
# endif

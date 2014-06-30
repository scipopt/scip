/* $Id: load_op.hpp 2920 2013-10-11 12:40:52Z bradbell $ */
# ifndef CPPAD_LOAD_OP_INCLUDED
# define CPPAD_LOAD_OP_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-13 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Eclipse Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */


namespace CppAD { // BEGIN_CPPAD_NAMESPACE
/*!
\defgroup load_op_hpp load_op.hpp
\{
\file load_op.hpp
Setting a variable so that it corresponds to current value of a VecAD element.
*/

/*!
Shared documentation for zero order forward mode implementation of 
op = LdpOp or LdvOp (not called).

The C++ source code corresponding to this operation is
\verbatim
	z = y[x]
\endverbatim
where y is a VecAD<Base> vector and x is an AD<Base> index. 
We define the index corresponding to y[x] by
\verbatim
	i_y_x = combined[ arg[0] + i_vec ]
\endverbatim
where i_vec is defined under the heading \a arg[1] below:

\tparam Base
base type for the operator; i.e., this operation was recorded
using AD< \a Base > and computations by this routine are done using type 
\a Base.

\param i_z
is the AD variable index corresponding to the variable z.

\param arg
\n
\a arg[0]
is the offset of this VecAD vector relative to the beginning 
of the \a combined VecAD array.
\n
\n 
\a arg[1] 
\n
If this is the LdpOp operation 
(the index x is a parameter), i_vec is defined by
\verbatim
	i_vec = arg[1]
\endverbatim
If this is the LdvOp operation 
(the index x is a variable), i_vec is defined by
\verbatim
	i_vec = floor( taylor[ arg[1] * nc_taylor + 0 ] )
\endverbatim
where floor(c) is the greatest integer less that or equal c.
\n
\a arg[2]
\b Input: The input value of \a arg[2] does not matter.
\n
\b Output: 
If y[x] is a parameter, \a arg[2] is set to zero 
(which is not a valid variable index).
If y[x] is a variable, 
\a arg[2] is set to the variable index corresponding to y[x]; i.e.  i_y_x.

\param num_par
is the number of parameters in \a parameter.

\param parameter
If y[x] is a parameter, \a parameter [ i_y_x ] is its value.

\param nc_taylor
number of columns in the matrix containing the Taylor coefficients.

\param taylor
\b Input: in LdvOp case, \a taylor[ arg[1] * nc_taylor + 0 ]
is used to compute the index in the definition of i_vec above
\n
\b Input: if y[x] is a variable, \a taylor[ i_y_x * nc_taylor + 0 ]
is the zero order Taylor coefficient for y[x].
\n
\b Output: \a taylor[ i_z * nc_taylor + 0 ]
is the zero order Taylor coefficient for the variable z.

\param nc_combined
is the total number of elements in the combined VecAD array.

\param variable
If \a variable [ \a arg[0] + i_vec ] is true,
y[x] is a variable.  Otherwise it is a parameter.

\param combined
\b Input: \a combined[ \a arg[0] - 1 ] 
is the number of elements in the VecAD vector containing this element.
\n
\b Input: \a combined[ \a arg[0] + i_vec ]
if y[x] is a variable, i_y_x
is its index in the Taylor coefficient array \a taylor.
Otherwise, i_y_x is its index in parameter array \a parameter.

\par Check User Errors
\li In the LdvOp case check that the index is with in range; i.e.
i_vec < combined[ \a arg[0] - 1 ] 

\par Checked Assertions 
\li combined != CPPAD_NULL
\li variable != CPPAD_NULL
\li NumArg(LdpOp) == 3
\li NumRes(LdpOp) == 1
\li 0 <  \a arg[0]
\li \a arg[0] + i_vec < nc_combined
\li i_vec < combined[ \a arg[0] - 1 ] 
\li if y[x] is a parameter, i_y_x < num_par
\li if y[x] is a variable, i_y_x < i_z
\li if x is a variable (LpvOp case), arg[1] < i_z
*/
template <class Base>
inline void forward_load_op_0(
	size_t         i_z         ,
	addr_t*        arg         , 
	size_t         num_par     ,
	const Base*    parameter   ,
	size_t         nc_taylor   ,
	Base*          taylor      ,
	size_t         nc_combined ,
	const bool*    variable    ,
	const size_t*  combined    )
{
	// This routine is only for documentaiton, it should not be used
	CPPAD_ASSERT_UNKNOWN( false );
}
/*!
Shared documentation for sparsity operations corresponding to 
op = LdpOp or LdvOp (not called).

<!-- define sparse_load_op -->
The C++ source code corresponding to this operation is
\verbatim
	z = v[x]
\endverbatim
where v is a VecAD<Base> vector and x is an AD<Base> index. 

\tparam Vector_set
is the type used for vectors of sets. It can be either
\c sparse_pack, \c sparse_set, or \c sparse_list.

\param op
is the code corresponding to this operator; i.e., LdpOp or LdvOp
(only used for error checking).

\param i_z
is the AD variable index corresponding to the variable z; i.e.,
the set with index \a i_z in \a var_sparsity is the sparsity pattern
correpsonding to z.

\param arg
\n
\a arg[0]
is the offset corresponding to this VecAD vector in the VecAD combined array.

\param num_combined
is the total number of elements in the VecAD combinded array.

\param combined
is the VecAD combined array.
\n
\n
\a combined[ \a arg[0] - 1 ]
is the index of the set corresponding to the vector v  in \a vecad_sparsity.
We use the notation i_v for this value; i.e.,
\verbatim
	i_v = combined[ \a arg[0] - 1 ]
\endverbatim

\param var_sparsity
The set with index \a i_z in \a var_sparsity is the sparsity pattern for z.
This is an output for forward mode operations,
and an input for reverse mode operations.

\param vecad_sparsity
The set with index \a i_v is the sparsity pattern for the vector v.
This is an input for forward mode operations.
For reverse mode operations,
the sparsity pattern for z is added to the sparsity pattern for v.

\par Checked Assertions 
\li NumArg(op) == 3
\li NumRes(op) == 1
\li 0         <  \a arg[0]
\li \a arg[0] < \a num_combined
\li i_v       < \a vecad_sparsity.n_set()
<!-- end sparse_load_op -->
*/
template <class Vector_set>
inline void sparse_load_op(
	OpCode              op             ,
	size_t              i_z            ,
	const addr_t*        arg           , 
	size_t              num_combined   ,
	const size_t*       combined       ,
	Vector_set&         var_sparsity   ,
	Vector_set&         vecad_sparsity )
{
	// This routine is only for documentaiton, it should not be used
	CPPAD_ASSERT_UNKNOWN( false );
}


/*!
Zero order forward mode implementation of op = LdpOp.

\copydetails forward_load_op_0
*/
template <class Base>
inline void forward_load_p_op_0(
	size_t         i_z         ,
	addr_t*        arg         , 
	size_t         num_par     ,
	const Base*    parameter   ,
	size_t         nc_taylor   ,
	Base*          taylor      ,
	size_t         nc_combined ,
	const bool*    variable    ,
	const size_t*  combined    )
{	size_t i_vec = arg[1];

	// Because the index is a parameter, this indexing error should be
	// caught and reported to the user at an when the tape is recording.
	CPPAD_ASSERT_UNKNOWN( i_vec < combined[ arg[0] - 1 ] );


	CPPAD_ASSERT_UNKNOWN( variable != CPPAD_NULL );
	CPPAD_ASSERT_UNKNOWN( combined != CPPAD_NULL );
	CPPAD_ASSERT_UNKNOWN( NumArg(LdpOp) == 3 );
	CPPAD_ASSERT_UNKNOWN( NumRes(LdpOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( 0 < arg[0] );
	CPPAD_ASSERT_UNKNOWN( arg[0] + i_vec < nc_combined );

	size_t combined_index = arg[0] + i_vec;
	size_t i_y_x          = combined[ combined_index ];	
	Base* z  = taylor + i_z * nc_taylor;
	if( variable[ combined_index ] )
	{	CPPAD_ASSERT_UNKNOWN( i_y_x < i_z );
		Base* y_x = taylor + i_y_x * nc_taylor;
		arg[2]    = i_y_x;
		z[0]      = y_x[0];
	}
	else
	{	CPPAD_ASSERT_UNKNOWN( i_y_x < num_par );
		Base y_x  = parameter[i_y_x];
		arg[2]    = 0;
		z[0]      = y_x;
	}
}

/*!
Zero order forward mode implementation of op = LdvOp.

\copydetails forward_load_op_0
*/
template <class Base>
inline void forward_load_v_op_0(
	size_t         i_z         ,
	addr_t*        arg         , 
	size_t         num_par     ,
	const Base*    parameter   ,
	size_t         nc_taylor   ,
	Base*          taylor      ,
	size_t         nc_combined ,
	const bool*    variable    ,
	const size_t*  combined    )
{
	CPPAD_ASSERT_UNKNOWN( variable != CPPAD_NULL );
	CPPAD_ASSERT_UNKNOWN( combined != CPPAD_NULL );
	CPPAD_ASSERT_UNKNOWN( NumArg(LdvOp) == 3 );
	CPPAD_ASSERT_UNKNOWN( NumRes(LdvOp) == 1 );
	CPPAD_ASSERT_UNKNOWN( 0 < arg[0] );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[1]) < i_z );

	size_t i_vec = Integer( taylor[ arg[1] * nc_taylor + 0 ] );
	CPPAD_ASSERT_KNOWN( 
		i_vec < combined[ arg[0] - 1 ] ,
		"VecAD: index during zero order forward sweep is out of range"
	);
	CPPAD_ASSERT_UNKNOWN( arg[0] + i_vec < nc_combined );

	size_t combined_index = arg[0] + i_vec;
	size_t i_y_x          = combined[ combined_index ];	


	Base* z  = taylor + i_z * nc_taylor;
	if( variable[ combined_index ] )
	{	CPPAD_ASSERT_UNKNOWN( i_y_x < i_z );
		Base* y_x = taylor + i_y_x * nc_taylor;
		arg[2]    = i_y_x;
		z[0]      = y_x[0];
	}
	else
	{	CPPAD_ASSERT_UNKNOWN( i_y_x < num_par );
		Base y_x  = parameter[i_y_x];
		arg[2]    = 0;
		z[0]      = y_x;
	}
}

/*!
Forward mode, except for zero order, for op = LdpOp or op = LdvOp

The C++ source code corresponding to this operation is
\verbatim
	z = y[x]
\endverbatim
where y is a VecAD<Base> vector and x is an AD<Base> or Base index. 

\tparam Base
base type for the operator; i.e., this operation was recorded
using AD< \a Base > and computations by this routine are done using type 
\a Base.

\param op
is the code corresponding to this operator; i.e., LdpOp or LdvOp
(only used for error checking).

\param q
is the lowest order of the Taylor coefficient that we are computing.

\param p
is the highest order of the Taylor coefficient that we are computing.

\param i_z
is the AD variable index corresponding to the variable z.

\param arg
\a arg[2]
If y[x] is a parameter, <code>arg[2]</code> is zero 
(which is not a valid variable index).
If y[x] is a variable, 
<code>arg[2]</code> is the variable index corresponding to y[x].

\param nc_taylor
number of columns in the matrix containing the Taylor coefficients.

\param taylor
\b Input: if y[x] is a variable, 
<code>taylor[ arg[2] * nc_taylor + k ]</code>
for k = 0 , ... , p,
is the k-order Taylor coefficient corresponding to y[x].
\n
\b Output: <code>taylor[ i_z * nc_taylor + d ]</code>
for k = q , ... , p,
is the k-order Taylor coefficient for the variable z.

\par Checked Assertions 
\li NumArg(op) == 3
\li NumRes(op) == 1
\li p < nc_taylor
\li 0 < q <= p 
\li size_t(arg[2]) < i_z
*/
template <class Base>
inline void forward_load_op(
	OpCode         op          ,
	size_t         q           ,
	size_t         p           ,
	size_t         i_z         ,
	const addr_t*  arg         , 
	size_t         nc_taylor   ,
	Base*          taylor      )
{

	CPPAD_ASSERT_UNKNOWN( NumArg(op) == 3 );
	CPPAD_ASSERT_UNKNOWN( NumRes(op) == 1 );
	CPPAD_ASSERT_UNKNOWN( p < nc_taylor );
	CPPAD_ASSERT_UNKNOWN( 0 < q && q <= p );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[2]) < i_z );

	Base* z      = taylor + i_z * nc_taylor;
	if( arg[2] > 0 )
	{	Base* y_x = taylor + arg[2] * nc_taylor;
		for(size_t d = q; d <= p; d++)
			z[d] = y_x[d];
	}
	else
	{	for(size_t d = q; d <= p; d++)
			z[d] = Base(0);
	}
}

/*!
Reverse mode for op = LdpOp or LdvOp.

The C++ source code corresponding to this operation is
\verbatim
	z = y[x]
\endverbatim
where y is a VecAD<Base> vector and x is an AD<Base> or Base index. 

This routine is given the partial derivatives of a function 
G(z , y[x] , w , u ... )
and it uses them to compute the partial derivatives of 
\verbatim
	H( y[x] , w , u , ... ) = G[ z( y[x] ) , y[x] , w , u , ... ]
\endverbatim

\tparam Base
base type for the operator; i.e., this operation was recorded
using AD< \a Base > and computations by this routine are done using type 
\a Base.

\param op
is the code corresponding to this operator; i.e., LdpOp or LdvOp
(only used for error checking).

\param d
highest order the Taylor coefficient that we are computing the partial
derivative with respect to.

\param i_z
is the AD variable index corresponding to the variable z.

\param arg
\a arg[2]
If y[x] is a parameter, \a arg[2] is zero 
(which is not a valid variable index).
If y[x] is a variable, 
\a arg[2] is the variable index corresponding to y[x].

\param nc_taylor
number of columns in the matrix containing the Taylor coefficients
(not used).

\param taylor
matrix of Taylor coefficients (not used).

\param nc_partial
number of colums in the matrix containing all the partial derivatives
(not used if \a arg[2] is zero).

\param partial
If \a arg[2] is zero, y[x] is a parameter
and no values need to be modified; i.e., \a partial is not used.
Otherwise, y[x] is a variable and:
\n
\n
\a partial [ \a i_z * \a nc_partial + k ] 
for k = 0 , ... , \a d
is the partial derivative of G
with respect to the k-th order Taylor coefficient for z.
\n
\n
If \a arg[2] is not zero,
\a partial [ \a arg[2] * \a nc_partial + k ]
for k = 0 , ... , \a d
is the partial derivative with respect to 
the k-th order Taylor coefficient for x.
On input, it corresponds to the function G,
and on output it corresponds to the the function H. 

\par Checked Assertions 
\li NumArg(op) == 3
\li NumRes(op) == 1
\li d < nc_taylor
\li size_t(arg[2]) < i_z
*/
template <class Base>
inline void reverse_load_op(
	OpCode         op          ,
	size_t         d           ,
	size_t         i_z         ,
	const addr_t*  arg         , 
	size_t         nc_taylor   ,
	const Base*    taylor      ,
	size_t         nc_partial  ,
	Base*          partial     )
{

	CPPAD_ASSERT_UNKNOWN( NumArg(op) == 3 );
	CPPAD_ASSERT_UNKNOWN( NumRes(op) == 1 );
	CPPAD_ASSERT_UNKNOWN( d < nc_taylor );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[2]) < i_z );

	if( arg[2] > 0 )
	{
		Base* pz   = partial + i_z    * nc_partial;
		Base* py_x = partial + arg[2] * nc_partial;
		size_t j = d + 1;
		while(j--)
			py_x[j]   += pz[j];
	}
}


/*!
Forward mode sparsity operations for LdpOp and LdvOp

\copydetails sparse_load_op
*/
template <class Vector_set>
inline void forward_sparse_load_op(
	OpCode             op             ,
	size_t             i_z            ,
	const addr_t*      arg            , 
	size_t             num_combined   ,
	const size_t*      combined       ,
	Vector_set&        var_sparsity   ,
	Vector_set&        vecad_sparsity )
{
	CPPAD_ASSERT_UNKNOWN( NumArg(op) == 3 );
	CPPAD_ASSERT_UNKNOWN( NumRes(op) == 1 );
	CPPAD_ASSERT_UNKNOWN( 0 < arg[0] );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < num_combined );
	size_t i_v = combined[ arg[0] - 1 ];
	CPPAD_ASSERT_UNKNOWN( i_v < vecad_sparsity.n_set() );

	var_sparsity.assignment(i_z, i_v, vecad_sparsity);

	return;
}


/*!
Reverse mode Jacobian sparsity operations for LdpOp and LdvOp

\copydetails sparse_load_op
*/
template <class Vector_set>
inline void reverse_sparse_jacobian_load_op(
	OpCode             op             ,
	size_t             i_z            ,
	const addr_t*      arg            , 
	size_t             num_combined   ,
	const size_t*      combined       ,
	Vector_set&        var_sparsity   ,
	Vector_set&        vecad_sparsity )
{
	CPPAD_ASSERT_UNKNOWN( NumArg(op) == 3 );
	CPPAD_ASSERT_UNKNOWN( NumRes(op) == 1 );
	CPPAD_ASSERT_UNKNOWN( 0 < arg[0] );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < num_combined );
	size_t i_v = combined[ arg[0] - 1 ];
	CPPAD_ASSERT_UNKNOWN( i_v < vecad_sparsity.n_set() );

	vecad_sparsity.binary_union(i_v, i_v, i_z, var_sparsity);

	return;
}


/*!
Reverse mode Hessian sparsity operations for LdpOp and LdvOp

This routine is given the sparsity patterns for 
G(z , v[x] , w , u ... )
and it uses them to compute the sparsity patterns for
\verbatim
	H( v[x] , w , u , ... ) = G[ z( v[x] ) , v[x] , w , u , ... ]
\endverbatim

<!-- replace sparse_load_op -->
The C++ source code corresponding to this operation is
\verbatim
	z = v[x]
\endverbatim
where v is a VecAD<Base> vector and x is an AD<Base> index. 

\tparam Vector_set
is the type used for vectors of sets. It can be either
\c sparse_pack, \c sparse_set, or \c sparse_list.

\param op
is the code corresponding to this operator; i.e., LdpOp or LdvOp
(only used for error checking).

\param i_z
is the AD variable index corresponding to the variable z; i.e.,
the set with index \a i_z in \a var_sparsity is the sparsity pattern
correpsonding to z.

\param arg
\n
\a arg[0]
is the offset corresponding to this VecAD vector in the VecAD combined array.

\param num_combined
is the total number of elements in the VecAD combinded array.

\param combined
is the VecAD combined array.
\n
\n
\a combined[ \a arg[0] - 1 ]
is the index of the set corresponding to the vector v  in \a vecad_sparsity.
We use the notation i_v for this value; i.e.,
\verbatim
	i_v = combined[ \a arg[0] - 1 ]
\endverbatim

\param var_sparsity
The set with index \a i_z in \a var_sparsity is the sparsity pattern for z.
This is an output for forward mode operations,
and an input for reverse mode operations.

\param vecad_sparsity
The set with index \a i_v is the sparsity pattern for the vector v.
This is an input for forward mode operations.
For reverse mode operations,
the sparsity pattern for z is added to the sparsity pattern for v.

\par Checked Assertions 
\li NumArg(op) == 3
\li NumRes(op) == 1
\li 0         <  \a arg[0]
\li \a arg[0] < \a num_combined
\li i_v       < \a vecad_sparsity.n_set()
<!-- end sparse_load_op -->

\param var_jacobian
\a var_jacobian[i_z] 
is false (true) if the Jacobian of G with respect to z is always zero 
(many be non-zero).

\param vecad_jacobian
\a vecad_jacobian[i_v] 
is false (true) if the Jacobian with respect to x is always zero 
(may be non-zero).
On input, it corresponds to the function G,
and on output it corresponds to the function H.

*/
template <class Vector_set>
inline void reverse_sparse_hessian_load_op(
	OpCode             op             ,
	size_t             i_z            ,
	const addr_t*      arg            , 
	size_t             num_combined   ,
	const size_t*      combined       ,
	Vector_set&        var_sparsity   ,
	Vector_set&        vecad_sparsity ,
	bool*              var_jacobian   ,
	bool*              vecad_jacobian )
{
	CPPAD_ASSERT_UNKNOWN( NumArg(op) == 3 );
	CPPAD_ASSERT_UNKNOWN( NumRes(op) == 1 );
	CPPAD_ASSERT_UNKNOWN( 0 < arg[0] );
	CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < num_combined );
	size_t i_v = combined[ arg[0] - 1 ];
	CPPAD_ASSERT_UNKNOWN( i_v < vecad_sparsity.n_set() );

	vecad_sparsity.binary_union(i_v, i_v, i_z, var_sparsity);

	vecad_jacobian[i_v] |= var_jacobian[i_z];

	return;
}


/*! \} */
} // END_CPPAD_NAMESPACE
# endif

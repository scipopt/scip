/* $Id: bender_quad.hpp 1643 2010-02-04 05:45:21Z bradbell $ */
# ifndef CPPAD_BENDER_QUAD_INCLUDED
# define CPPAD_BENDER_QUAD_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-10 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
/*
$begin BenderQuad$$
$spell
	cppad.hpp
	BAvector
	gx
	gxx
	CppAD
	Fy
	dy
	Jacobian
	ADvector
	const
	dg
	ddg
$$


$index Hessian, Bender$$
$index Jacobian, Bender$$
$index BenderQuad$$
$section Computing Jacobian and Hessian of Bender's Reduced Objective$$

$head Syntax$$
$syntax%%
%# include <cppad/cppad.hpp>
BenderQuad(%x%, %y%, %fun%, %g%, %gx%, %gxx%)%$$  

$head See Also$$
$cref/opt_val_hes/$$

$head Problem$$
The type $cref/ADvector/BenderQuad/ADvector/$$ cannot be determined
form the arguments above 
(currently the type $italic ADvector$$ must be 
$syntax%CPPAD_TEST_VECTOR<%Base%>%$$.)
This will be corrected in the future by requiring $italic Fun$$
to define $syntax%%Fun%::vector_type%$$ which will specify the
type $italic ADvector$$.

$head Purpose$$
We are given the optimization problem
$latex \[
\begin{array}{rcl}
	{\rm minimize} & F(x, y) & {\rm w.r.t.} \; (x, y) \in \R^n \times \R^m
\end{array}
\] $$
that is convex with respect to $latex y$$.
In addition, we are given a set of equations $latex H(x, y)$$
such that 
$latex \[
	H[ x , Y(x) ] = 0 \;\; \Rightarrow \;\; F_y [ x , Y(x) ] = 0
\] $$
(In fact, it is often the case that $latex H(x, y) = F_y (x, y)$$.)
Furthermore, it is easy to calculate a Newton step for these equations; i.e.,
$latex \[
	dy = - [ \partial_y H(x, y)]^{-1} H(x, y)
\] $$
The purpose of this routine is to compute the 
value, Jacobian, and Hessian of the reduced objective function
$latex \[
	G(x) = F[ x , Y(x) ]
\] $$
Note that if only the value and Jacobian are needed, they can be
computed more quickly using the relations
$latex \[
	G^{(1)} (x) = \partial_x F [x, Y(x) ]
\] $$ 

$head x$$
The $code BenderQuad$$ argument $italic x$$ has prototype
$syntax%
	const %BAvector% &%x%
%$$
(see $xref/BenderQuad/BAvector/BAvector/$$ below)
and its size must be equal to $italic n$$.
It specifies the point at which we evaluating 
the reduced objective function and its derivatives.


$head y$$
The $code BenderQuad$$ argument $italic y$$ has prototype
$syntax%
	const %BAvector% &%y%
%$$
and its size must be equal to $italic m$$.
It must be equal to $latex Y(x)$$; i.e., 
it must solve the problem in $latex y$$ for this given value of $latex x$$
$latex \[
\begin{array}{rcl}
	{\rm minimize} & F(x, y) & {\rm w.r.t.} \; y \in \R^m
\end{array}
\] $$

$head fun$$
The $code BenderQuad$$ object $italic fun$$ 
must support the member functions listed below.
The $syntax%AD<%Base%>%$$ arguments will be variables for
a tape created by a call to $cref%Independent%$$ from $code BenderQuad$$
(hence they can not be combined with variables corresponding to a 
different tape). 

$subhead fun.f$$
The $code BenderQuad$$ argument $italic fun$$ supports the syntax
$syntax%
	%f% = %fun%.f(%x%, %y%)
%$$
The $syntax%%fun%.f%$$ argument $italic x$$ has prototype
$syntax%
	const %ADvector% &%x%
%$$
(see $xref/BenderQuad/ADvector/ADvector/$$ below)
and its size must be equal to $italic n$$.
The $syntax%%fun%.f%$$ argument $italic y$$ has prototype
$syntax%
	const %ADvector% &%y%
%$$
and its size must be equal to $italic m$$.
The $syntax%%fun%.f%$$ result $italic f$$ has prototype
$syntax%
	%ADvector% %f%
%$$
and its size must be equal to one.
The value of $italic f$$ is
$latex \[
	f = F(x, y)
\] $$.

$subhead fun.h$$
The $code BenderQuad$$ argument $italic fun$$ supports the syntax
$syntax%
	%h% = %fun%.h(%x%, %y%)
%$$
The $syntax%%fun%.h%$$ argument $italic x$$ has prototype
$syntax%
	const %ADvector% &%x%
%$$
and its size must be equal to $italic n$$.
The $syntax%%fun%.h%$$ argument $italic y$$ has prototype
$syntax%
	const %BAvector% &%y%
%$$
and its size must be equal to $italic m$$.
The $syntax%%fun%.h%$$ result $italic h$$ has prototype
$syntax%
	%ADvector% %h%
%$$
and its size must be equal to $italic m$$.
The value of $italic h$$ is
$latex \[
	h = H(x, y)
\] $$.

$subhead fun.dy$$
The $code BenderQuad$$ argument $italic fun$$ supports the syntax
$syntax%
	%dy% = %fun%.dy(%x%, %y%, %h%)

%x%
%$$
The $syntax%%fun%.dy%$$ argument $italic x$$ has prototype
$syntax%
	const %BAvector% &%x%
%$$
and its size must be equal to $italic n$$.
Its value will be exactly equal to the $code BenderQuad$$ argument 
$italic x$$ and values depending on it can be stored as private objects
in $italic f$$ and need not be recalculated.
$syntax%

%y%
%$$
The $syntax%%fun%.dy%$$ argument $italic y$$ has prototype
$syntax%
	const %BAvector% &%y%
%$$
and its size must be equal to $italic m$$.
Its value will be exactly equal to the $code BenderQuad$$ argument 
$italic y$$ and values depending on it can be stored as private objects
in $italic f$$ and need not be recalculated.
$syntax%

%h%
%$$
The $syntax%%fun%.dy%$$ argument $italic h$$ has prototype
$syntax%
	const %ADvector% &%h%
%$$
and its size must be equal to $italic m$$.
$syntax%

%dy%
%$$
The $syntax%%fun%.dy%$$ result $italic dy$$ has prototype
$syntax%
	%ADvector% %dy%
%$$
and its size must be equal to $italic m$$.
The return value $italic dy$$ is given by
$latex \[
	dy = - [ \partial_y H (x , y) ]^{-1} h
\] $$
Note that if $italic h$$ is equal to $latex H(x, y)$$,
$latex dy$$ is the Newton step for finding a zero
of $latex H(x, y)$$ with respect to $latex y$$;
i.e., 
$latex y + dy$$ is an approximate solution for the equation
$latex H (x, y + dy) = 0$$. 

$head g$$
The argument $italic g$$ has prototype
$syntax%
	%BAvector% &%g%
%$$
and has size one.
The input value of its element does not matter.
On output,
it contains the value of $latex G (x)$$; i.e.,
$latex \[
	g[0] = G (x)
\] $$

$head gx$$
The argument $italic gx$$ has prototype
$syntax%
	%BAvector% &%gx%
%$$
and has size $latex n $$.
The input values of its elements do not matter.
On output,
it contains the Jacobian of $latex G (x)$$; i.e.,
for $latex j = 0 , \ldots , n-1$$, 
$latex \[
	gx[ j ] = G^{(1)} (x)_j
\] $$

$head gxx$$
The argument $italic gx$$ has prototype
$syntax%
	%BAvector% &%gxx%
%$$
and has size $latex n \times n$$.
The input values of its elements do not matter.
On output,
it contains the Hessian of $latex G (x)$$; i.e.,
for $latex i = 0 , \ldots , n-1$$, and
$latex j = 0 , \ldots , n-1$$, 
$latex \[
	gxx[ i * n + j ] = G^{(2)} (x)_{i,j}
\] $$

$head BAvector$$
The type $italic BAvector$$ must be a 
$xref/SimpleVector/$$ class. 
We use $italic Base$$ to refer to the type of the elements of 
$italic BAvector$$; i.e.,
$syntax%
	%BAvector%::value_type
%$$

$head ADvector$$
The type $italic ADvector$$ must be a 
$xref/SimpleVector/$$ class with elements of type 
$syntax%AD<%Base%>%$$; i.e.,
$syntax%
	%ADvector%::value_type
%$$
must be the same type as
$syntax%
	AD< %BAvector%::value_type >
%$$.


$head Example$$
$children%
	example/bender_quad.cpp
%$$
The file
$xref/BenderQuad.cpp/$$
contains an example and test of this operation.   
It returns true if it succeeds and false otherwise.


$end
-----------------------------------------------------------------------------
*/

namespace CppAD { // BEGIN CppAD namespace

template <class BAvector, class Fun>
void BenderQuad(
	const BAvector   &x     , 
	const BAvector   &y     , 
	Fun               fun   , 
	BAvector         &g     ,
	BAvector         &gx    ,
	BAvector         &gxx   )
{	// determine the base type
	typedef typename BAvector::value_type Base;

	// check that BAvector is a SimpleVector class
	CheckSimpleVector<Base, BAvector>();

	// declare the ADvector type
	typedef CPPAD_TEST_VECTOR< AD<Base> > ADvector;

	// size of the x and y spaces
	size_t n = x.size();
	size_t m = y.size();

	// check the size of gx and gxx
	CPPAD_ASSERT_KNOWN(
		g.size() == 1,
		"BenderQuad: size of the vector g is not equal to 1"
	);
	CPPAD_ASSERT_KNOWN(
		gx.size() == n,
		"BenderQuad: size of the vector gx is not equal to n"
	);
	CPPAD_ASSERT_KNOWN(
		gxx.size() == n * n,
		"BenderQuad: size of the vector gxx is not equal to n * n"
	);

	// some temporary indices
	size_t i, j;

	// variable versions x
	ADvector vx(n);
	for(j = 0; j < n; j++)
		vx[j] = x[j];
	
	// declare the independent variables
	Independent(vx);

	// evaluate h = H(x, y) 
	ADvector h(m);
	h = fun.h(vx, y);

	// evaluate dy (x) = Newton step as a function of x through h only
	ADvector dy(m);
	dy = fun.dy(x, y, h);

	// variable version of y
	ADvector vy(m);
	for(j = 0; j < m; j++)
		vy[j] = y[j] + dy[j];

	// evaluate G~ (x) = F [ x , y + dy(x) ] 
	ADvector gtilde(1);
	gtilde = fun.f(vx, vy);

	// AD function object that corresponds to G~ (x)
	// We will make heavy use of this tape, so optimize it
	ADFun<Base> Gtilde;
	Gtilde.Dependent(vx, gtilde); 
	Gtilde.optimize();

	// value of G(x)
	g = Gtilde.Forward(0, x);

	// initial forward direction vector as zero
	BAvector dx(n);
	for(j = 0; j < n; j++)
		dx[j] = Base(0);

	// weight, first and second order derivative values
	BAvector dg(1), w(1), ddw(2 * n);
	w[0] = 1.;


	// Jacobian and Hessian of G(x) is equal Jacobian and Hessian of Gtilde
	for(j = 0; j < n; j++)
	{	// compute partials in x[j] direction
		dx[j] = Base(1);
		dg    = Gtilde.Forward(1, dx);
		gx[j] = dg[0];

		// restore the dx vector to zero
		dx[j] = Base(0);

		// compute second partials w.r.t x[j] and x[l]  for l = 1, n
		ddw = Gtilde.Reverse(2, w);
		for(i = 0; i < n; i++)
			gxx[ i * n + j ] = ddw[ i * 2 + 1 ];
	}

	return;
}
	
} // END CppAD namespace

# endif

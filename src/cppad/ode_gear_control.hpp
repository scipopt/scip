/* $Id: ode_gear_control.hpp 2057 2011-08-11 14:07:11Z bradbell $ */
# ifndef CPPAD_ODE_GEAR_CONTROL_INCLUDED
# define CPPAD_ODE_GEAR_CONTROL_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin OdeGearControl$$
$spell
	cppad.hpp
	CppAD
	xf
	xi
	smin
	smax
	eabs
	ef
	maxabs
	nstep
	tf
	sini
	erel
	dep
	const
	tb
	ta
	exp
$$

$index OdeGearControl$$
$index control, Ode Gear$$
$index error, Gear Ode$$
$index differential, Ode Gear control$$
$index equation, Ode Gear control$$

 
$section An Error Controller for Gear's Ode Solvers$$

$head Syntax$$
$syntax%# include <cppad/ode_gear_control.hpp>
%$$
$syntax%%xf% = OdeGearControl(%F%, %M%, %ti%, %tf%, %xi%,
	%smin%, %smax%, %sini%, %eabs%, %erel%, %ef% , %maxabs%, %nstep% )%$$


$head Purpose$$
Let $latex \R$$ denote the real numbers
and let $latex f : \R \times \R^n \rightarrow \R^n$$ be a smooth function.
We define $latex X : [ti , tf] \rightarrow \R^n$$ by 
the following initial value problem:
$latex \[
\begin{array}{rcl}
	X(ti)  & = & xi    \\
	X'(t)  & = & f[t , X(t)] 
\end{array}
\] $$
The routine $cref/OdeGear/$$ is a stiff multi-step method that
can be used to approximate the solution to this equation.
The routine $code OdeGearControl$$ sets up this multi-step method
and controls the error during such an approximation.

$head Include$$
The file $code cppad/ode_gear_control.hpp$$ 
is included by $code cppad/cppad.hpp$$
but it can also be included separately with out the rest of 
the $code CppAD$$ routines.

$head Notation$$
The template parameter types $xref/OdeGearControl/Scalar/Scalar/$$ and
$xref/OdeGearControl/Vector/Vector/$$ are documented below.

$head xf$$
The return value $italic xf$$ has the prototype
$syntax%
	%Vector% %xf%
%$$
and the size of $italic xf$$ is equal to $italic n$$
(see description of $xref/OdeGear/Vector/Vector/$$ below).
It is the approximation for $latex X(tf)$$.

$head Fun$$
The class $italic Fun$$ 
and the object $italic F$$ satisfy the prototype
$syntax%
	%Fun% &%F%
%$$
This must support the following set of calls
$syntax%
	%F%.Ode(%t%, %x%, %f%)
	%F%.Ode_dep(%t%, %x%, %f_x%)
%$$

$subhead t$$
The argument $italic t$$ has prototype
$syntax%
	const %Scalar% &%t%
%$$
(see description of $xref/OdeGear/Scalar/Scalar/$$ below). 

$subhead x$$
The argument $italic x$$ has prototype
$syntax%
	const %Vector% &%x%
%$$
and has size $italic N$$
(see description of $xref/OdeGear/Vector/Vector/$$ below). 

$subhead f$$
The argument $italic f$$ to $syntax%%F%.Ode%$$ has prototype
$syntax%
	%Vector% &%f%
%$$
On input and output, $italic f$$ is a vector of size $italic N$$
and the input values of the elements of $italic f$$ do not matter.
On output,
$italic f$$ is set equal to $latex f(t, x)$$
(see $italic f(t, x)$$ in $xref/OdeGear/Purpose/Purpose/$$). 

$subhead f_x$$
The argument $italic f_x$$ has prototype
$syntax%
	%Vector% &%f_x%
%$$
On input and output, $italic f_x$$ is a vector of size $latex N * N$$
and the input values of the elements of $italic f_x$$ do not matter.
On output, 
$latex \[
	f\_x [i * n + j] = \partial_{x(j)} f_i ( t , x )
\] $$ 

$subhead Warning$$
The arguments $italic f$$, and $italic f_x$$
must have a call by reference in their prototypes; i.e.,
do not forget the $code &$$ in the prototype for 
$italic f$$ and $italic f_x$$.

$head M$$
The argument $italic M$$ has prototype
$syntax%
	size_t %M%
%$$
It specifies the order of the multi-step method; i.e.,
the order of the approximating polynomial
(after the initialization process).
The argument $italic M$$ must greater than or equal one.

$head ti$$
The argument $italic ti$$ has prototype
$syntax%
	const %Scalar% &%ti%
%$$
It specifies the initial time for the integration of 
the differential equation.

$head tf$$
The argument $italic tf$$ has prototype
$syntax%
	const %Scalar% &%tf%
%$$
It specifies the final time for the integration of 
the differential equation.

$head xi$$
The argument $italic xi$$ has prototype
$syntax%
	const %Vector% &%xi%
%$$
and size $italic n$$.
It specifies value of $latex X(ti)$$.

$head smin$$
The argument $italic smin$$ has prototype
$syntax%
	const %Scalar% &%smin%
%$$
The minimum value of $latex T[M] -  T[M-1]$$ in a call to $code OdeGear$$
will be $latex smin$$ except for the last two calls where it may be
as small as $latex smin / 2$$.
The value of $italic smin$$ must be less than or equal $italic smax$$.

$head smax$$
The argument $italic smax$$ has prototype
$syntax%
	const %Scalar% &%smax%
%$$
It specifies the maximum step size to use during the integration; 
i.e., the maximum value for $latex T[M] - T[M-1]$$ 
in a call to $code OdeGear$$.

$head sini$$
The argument $italic sini$$ has prototype
$syntax%
	%Scalar% &%sini%
%$$
The value of $italic sini$$ is the minimum 
step size to use during initialization of the multi-step method; i.e.,
for calls to $code OdeGear$$ where $latex m < M$$.
The value of $italic sini$$ must be less than or equal $italic smax$$
(and can also be less than $italic smin$$).

$head eabs$$
The argument $italic eabs$$ has prototype
$syntax%
	const %Vector% &%eabs%
%$$
and size $italic n$$.
Each of the elements of $italic eabs$$ must be 
greater than or equal zero.
It specifies a bound for the absolute
error in the return value $italic xf$$ as an approximation for $latex X(tf)$$.
(see the 
$xref/OdeGearControl/Error Criteria Discussion/error criteria discussion/$$ 
below). 

$head erel$$
The argument $italic erel$$ has prototype
$syntax%
	const %Scalar% &%erel%
%$$
and is greater than or equal zero.
It specifies a bound for the relative 
error in the return value $italic xf$$ as an approximation for $latex X(tf)$$
(see the 
$xref/OdeGearControl/Error Criteria Discussion/error criteria discussion/$$ 
below). 

$head ef$$
The argument value $italic ef$$ has prototype
$syntax%
	%Vector% &%ef%
%$$
and size $italic n$$.
The input value of its elements does not matter.
On output, 
it contains an estimated bound for the 
absolute error in the approximation $italic xf$$; i.e.,
$latex \[
	ef_i > | X( tf )_i - xf_i |
\] $$

$head maxabs$$
The argument $italic maxabs$$ is optional in the call to $code OdeGearControl$$.
If it is present, it has the prototype
$syntax%
	%Vector% &%maxabs%
%$$
and size $italic n$$.
The input value of its elements does not matter.
On output, 
it contains an estimate for the 
maximum absolute value of $latex X(t)$$; i.e.,
$latex \[
	maxabs[i] \approx \max \left\{ 
		| X( t )_i | \; : \;  t \in [ti, tf] 
	\right\}
\] $$

$head nstep$$
The argument $italic nstep$$ has the prototype
$syntax%
	%size_t% &%nstep%
%$$
Its input value does not matter and its output value
is the number of calls to $xref/OdeGear/$$
used by $code OdeGearControl$$.

$head Error Criteria Discussion$$
The relative error criteria $italic erel$$ and
absolute error criteria $italic eabs$$ are enforced during each step of the
integration of the ordinary differential equations.
In addition, they are inversely scaled by the step size so that
the total error bound is less than the sum of the error bounds.
To be specific, if $latex \tilde{X} (t)$$ is the approximate solution
at time $latex t$$, 
$italic ta$$ is the initial step time,
and $italic tb$$ is the final step time,
$latex \[
\left| \tilde{X} (tb)_j  - X (tb)_j \right| 
\leq 
\frac{tf - ti}{tb - ta}
\left[ eabs[j] + erel \;  | \tilde{X} (tb)_j | \right] 
\] $$
If $latex X(tb)_j$$ is near zero for some $latex tb \in [ti , tf]$$,
and one uses an absolute error criteria $latex eabs[j]$$ of zero,
the error criteria above will force $code OdeGearControl$$
to use step sizes equal to 
$xref/OdeGearControl/smin/smin/$$
for steps ending near $latex tb$$.
In this case, the error relative to $italic maxabs$$ can be judged after
$code OdeGearControl$$ returns.
If $italic ef$$ is to large relative to $italic maxabs$$, 
$code OdeGearControl$$ can be called again 
with a smaller value of $italic smin$$.

$head Scalar$$
The type $italic Scalar$$ must satisfy the conditions
for a $xref/NumericType/$$ type.
The routine $xref/CheckNumericType/$$ will generate an error message
if this is not the case.
In addition, the following operations must be defined for 
$italic Scalar$$ objects $italic a$$ and $italic b$$:

$table
$bold Operation$$ $cnext $bold Description$$  $rnext
$syntax%%a% <= %b%$$ $cnext
	returns true (false) if $italic a$$ is less than or equal 
	(greater than) $italic b$$.
$rnext
$syntax%%a% == %b%$$ $cnext
	returns true (false) if $italic a$$ is equal to $italic b$$.
$rnext
$syntax%log(%a%)%$$ $cnext
	returns a $italic Scalar$$ equal to the logarithm of $italic a$$
$rnext
$syntax%exp(%a%)%$$ $cnext
	returns a $italic Scalar$$ equal to the exponential of $italic a$$
$tend


$head Vector$$
The type $italic Vector$$ must be a $xref/SimpleVector/$$ class with
$xref/SimpleVector/Elements of Specified Type/elements of type Scalar/$$.
The routine $xref/CheckSimpleVector/$$ will generate an error message
if this is not the case.

$head Example$$
$children%
	example/ode_gear_control.cpp
%$$
The file
$xref/OdeGearControl.cpp/$$
contains an example and test a test of using this routine.
It returns true if it succeeds and false otherwise.

$head Theory$$
Let $latex e(s)$$ be the error as a function of the
step size $latex s$$ and suppose that there is a constant
$latex K$$ such that $latex e(s) = K s^m$$.
Let $latex a$$ be our error bound.
Given the value of $latex e(s)$$, a step of size $latex \lambda s$$
would be ok provided that
$latex \[
\begin{array}{rcl}
	a  & \geq & e( \lambda s ) (tf - ti) / ( \lambda s ) \\
	a  & \geq & K \lambda^m s^m (tf - ti) / ( \lambda s ) \\
	a  & \geq & \lambda^{m-1} s^{m-1} (tf - ti) e(s) / s^m \\
	a  & \geq & \lambda^{m-1} (tf - ti) e(s) / s           \\
	\lambda^{m-1} & \leq & \frac{a}{e(s)} \frac{s}{tf - ti}
\end{array}
\] $$
Thus if the right hand side of the last inequality is greater 
than or equal to one, the step of size $latex s$$ is ok. 

$head Source Code$$
The source code for this routine is in the file
$code cppad/ode_gear_control.hpp$$.

$end
--------------------------------------------------------------------------
*/

// link exp and log for float and double
# include <cppad/base_require.hpp>

# include <cppad/ode_gear.hpp>

namespace CppAD { // Begin CppAD namespace

template <class Scalar, class Vector, class Fun>
Vector OdeGearControl(
	Fun             &F     , 
	size_t           M     ,
	const Scalar    &ti    , 
	const Scalar    &tf    , 
	const Vector    &xi    , 
	const Scalar    &smin  , 
	const Scalar    &smax  , 
	Scalar          &sini  , 
	const Vector    &eabs  , 
	const Scalar    &erel  , 
	Vector          &ef    ,
	Vector          &maxabs,
	size_t          &nstep ) 
{
	// check simple vector class specifications
	CheckSimpleVector<Scalar, Vector>();

	// dimension of the state space
	size_t n = xi.size();

	CPPAD_ASSERT_KNOWN(
		M >= 1,
		"Error in OdeGearControl: M is less than one"
	);
	CPPAD_ASSERT_KNOWN(
		smin <= smax,
		"Error in OdeGearControl: smin is greater than smax"
	);
	CPPAD_ASSERT_KNOWN(
		sini <= smax,
		"Error in OdeGearControl: sini is greater than smax"
	);
	CPPAD_ASSERT_KNOWN(
		eabs.size() == n,
		"Error in OdeGearControl: size of eabs is not equal to n"
	);
	CPPAD_ASSERT_KNOWN(
		maxabs.size() == n,
		"Error in OdeGearControl: size of maxabs is not equal to n"
	);

	// some constants
	const Scalar zero(0);
	const Scalar one(1);
	const Scalar one_plus( Scalar(3) / Scalar(2) );
	const Scalar two(2);
	const Scalar ten(10);

	// temporary indices
	size_t i, k;

	// temporary Scalars
	Scalar step, sprevious, lambda, axi, a, root, r;

	// vectors of Scalars
	Vector T  (M + 1);
	Vector X( (M + 1) * n );
	Vector e(n);
	Vector xf(n);

	// initial integer values
	size_t m = 1;
	nstep    = 0;

	// initialize T
	T[0] = ti;

	// initialize X, ef, maxabs
	for(i = 0; i < n; i++) 
	for(i = 0; i < n; i++)
	{	X[i] = xi[i];
		ef[i] = zero;
		X[i]  = xi[i];
		if( zero <= xi[i] )
			maxabs[i] = xi[i];
		else	maxabs[i] = - xi[i];

	}  

	// initial step size
	step = smin;

	while( T[m-1] < tf )
	{	sprevious = step;

	 	// check maximum
		if( smax <= step )
			step = smax;

		// check minimum
		if( m < M )
		{	if( step <= sini )
				step = sini;
		}
		else	if( step <= smin )
				step = smin;

		// check if near the end
		if( tf <= T[m-1] + one_plus * step )
			T[m] = tf;
		else	T[m] = T[m-1] + step;

		// try using this step size
		nstep++;
		OdeGear(F, m, n, T, X, e);
		step = T[m] - T[m-1];

		// compute value of lambda for this step
		lambda = Scalar(10) *  sprevious / step;
		for(i = 0; i < n; i++)
		{	axi = X[m * n + i];
			if( axi <= zero )
				axi = - axi;
			a  = eabs[i] + erel * axi;
			if( e[i] > zero )
			{	if( m == 1 )
					root = (a / e[i]) / ten;
				else
				{	r = ( a / e[i] ) * step / (tf - ti);
					root = exp( log(r) / Scalar(m-1) ); 
				}
				if( root <= lambda )
					lambda = root;
			}
		}

		bool advance;
		if( m == M )
			advance = one <= lambda || step <= one_plus * smin;
		else	advance = one <= lambda || step <= one_plus * sini; 


		if( advance )
		{	// accept the results of this time step
			CPPAD_ASSERT_UNKNOWN( m <= M );
			if( m == M )
			{	// shift for next step
				for(k = 0; k < m; k++)
				{	T[k] = T[k+1];
					for(i = 0; i < n; i++)
						X[k*n + i] = X[(k+1)*n + i];
				}
			}
			// update ef and maxabs
			for(i = 0; i < n; i++)
			{	ef[i] = ef[i] + e[i];
				axi = X[m * n + i];
				if( axi <= zero )
					axi = - axi;
				if( axi > maxabs[i] )
					maxabs[i] = axi;
			}
			if( m != M )
				m++;  // all we need do in this case
		}

		// new step suggested by error criteria 
		step = std::min(lambda , ten) * step / two;
	}
	for(i = 0; i < n; i++)
		xf[i] = X[(m-1) * n + i];

	return xf;
}

} // End CppAD namespace 

# endif

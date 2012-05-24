/* $Id: ode_err_control.hpp 2057 2011-08-11 14:07:11Z bradbell $ */
# ifndef CPPAD_ODE_ERR_CONTROL_INCLUDED
# define CPPAD_ODE_ERR_CONTROL_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin OdeErrControl$$
$spell
	cppad.hpp
	nstep
	maxabs
	exp
	scur
	CppAD
	xf
	tf
	xi
	smin
	smax
	eabs
	erel
	ef
	ta
	tb
	xa
	xb
	const
	eb
$$

$index OdeErrControl$$
$index ODE, control error$$
$index control, ODE error$$
$index error, control ODE$$
$index differential, ODE error control$$
$index equation, ODE error control$$

 
$section An Error Controller for ODE Solvers$$

$head Syntax$$
$code # include <cppad/ode_err_control.hpp>$$
$pre
$$
$syntax%%xf% = OdeErrControl(%method%, %ti%, %tf%, %xi%,
	%smin%, %smax%, %scur%, %eabs%, %erel%, %ef% , %maxabs%, %nstep% )%$$


$head Description$$
Let $latex \R$$ denote the real numbers
and let $latex F : \R \times \R^n \rightarrow \R^n$$ be a smooth function.
We define $latex X : [ti , tf] \rightarrow \R^n$$ by 
the following initial value problem:
$latex \[
\begin{array}{rcl}
	X(ti)  & = & xi    \\
	X'(t)  & = & F[t , X(t)] 
\end{array}
\] $$
The routine $code OdeErrControl$$ can be used to adjust the step size
used an arbitrary integration methods in order to be as fast as possible
and still with in a requested error bound.

$head Include$$
The file $code cppad/ode_err_control.hpp$$ is included by 
$code cppad/cppad.hpp$$
but it can also be included separately with out the rest of 
the $code CppAD$$ routines.

$head Notation$$
The template parameter types $xref/OdeErrControl/Scalar/Scalar/$$ and
$xref/OdeErrControl/Vector/Vector/$$ are documented below.

$head xf$$
The return value $italic xf$$ has the prototype
$syntax%
	%Vector% %xf%
%$$
(see description of $xref/OdeErrControl/Vector/Vector/$$ below). 
and the size of $italic xf$$ is equal to $italic n$$.
If $italic xf$$ contains not a number $cref/nan/$$,
see the discussion of $cref/step/OdeErrControl/Method/Nan/$$.

$head Method$$
The class $italic Method$$
and the object $italic method$$ satisfy the following syntax
$syntax%
	%Method% &%method%
%$$
The object $italic method$$ must support $code step$$ and 
$code order$$ member functions defined below:

$subhead step$$
The syntax
$syntax%
	%method%.step(%ta%, %tb%, %xa%, %xb%, %eb%)
%$$
executes one step of the integration method. 
$syntax%

%ta%
%$$
The argument $italic ta$$ has prototype
$syntax%
	const %Scalar% &%ta%
%$$
It specifies the initial time for this step in the 
ODE integration.
(see description of $xref/OdeErrControl/Scalar/Scalar/$$ below). 
$syntax%

%tb%
%$$
The argument $italic tb$$ has prototype
$syntax%
	const %Scalar% &%tb%
%$$
It specifies the final time for this step in the 
ODE integration.
$syntax%

%xa%
%$$
The argument $italic xa$$ has prototype
$syntax%
	const %Vector% &%xa%
%$$
and size $italic n$$.
It specifies the value of $latex X(ta)$$.
(see description of $xref/OdeErrControl/Vector/Vector/$$ below). 
$syntax%

%xb%
%$$
The argument value $italic xb$$ has prototype
$syntax%
	%Vector% &%xb%
%$$
and size $italic n$$.
The input value of its elements does not matter.
On output, 
it contains the approximation for $latex X(tb)$$ that the method obtains.
$syntax%

%eb%
%$$
The argument value $italic eb$$ has prototype
$syntax%
	%Vector% &%eb%
%$$
and size $italic n$$.
The input value of its elements does not matter.
On output, 
it contains an estimate for the error in the approximation $italic xb$$.
It is assumed (locally) that the error bound in this approximation 
nearly equal to $latex K (tb - ta)^m$$ 
where $italic K$$ is a fixed constant and $italic m$$
is the corresponding argument to $code CodeControl$$.

$subhead Nan$$
If any element of the vector $italic eb$$ or $italic xb$$ are
not a number $code nan$$,
the current step is considered to large.
If this happens with the current step size equal to $italic smin$$, 
$code OdeErrControl$$ returns with $italic xf$$ and $italic ef$$ as vectors
of $code nan$$.

$subhead order$$
If $italic m$$ is $code size_t$$,
the object $italic method$$ must also support the following syntax 
$syntax%
	%m% = %method%.order()
%$$
The return value $italic m$$ is the order of the error estimate;
i.e., there is a constant K such that if $latex ti \leq ta \leq tb \leq tf$$,
$latex \[
	| eb(tb) | \leq K | tb - ta |^m
\] $$
where $italic ta$$, $italic tb$$, and $italic eb$$ are as in 
$syntax%%method%.step(%ta%, %tb%, %xa%, %xb%, %eb%)%$$


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
The step size during a call to $italic method$$ is defined as
the corresponding value of $latex tb - ta$$.
If $latex tf - ti \leq smin$$,
the integration will be done in one step of size $italic tf - ti$$.
Otherwise,
the minimum value of $italic tb - ta$$ will be $latex smin$$
except for the last two calls to $italic method$$ where it may be
as small as $latex smin / 2$$.

$head smax$$
The argument $italic smax$$ has prototype
$syntax%
	const %Scalar% &%smax%
%$$
It specifies the maximum step size to use during the integration; 
i.e., the maximum value for $latex tb - ta$$ in a call to $italic method$$.
The value of $italic smax$$ must be greater than or equal $italic smin$$.

$head scur$$
The argument $italic scur$$ has prototype
$syntax%
	%Scalar% &%scur%
%$$
The value of $italic scur$$ is the suggested next step size,
based on error criteria, to try in the next call to $italic method$$.
On input it corresponds to the first call to $italic method$$,
in this call to $code OdeErrControl$$ (where $latex ta = ti$$).
On output it corresponds to the next call to $italic method$$,
in a subsequent call to $code OdeErrControl$$ (where $italic ta = tf$$).

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
$xref/OdeErrControl/Error Criteria Discussion/error criteria discussion/$$ 
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
$xref/OdeErrControl/Error Criteria Discussion/error criteria discussion/$$ 
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
If on output $italic ef$$ contains not a number $code nan$$,
see the discussion of $cref/step/OdeErrControl/Method/Nan/$$.

$head maxabs$$
The argument $italic maxabs$$ is optional in the call to $code OdeErrControl$$.
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
The argument $italic nstep$$ is optional in the call to $code OdeErrControl$$.
If it is present, it has the prototype
$syntax%
	%size_t% &%nstep%
%$$
Its input value does not matter and its output value
is the number of calls to $syntax%%method%.step%$$ 
used by $code OdeErrControl$$.

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
the error criteria above will force $code OdeErrControl$$
to use step sizes equal to 
$xref/OdeErrControl/smin/smin/$$
for steps ending near $latex tb$$.
In this case, the error relative to $italic maxabs$$ can be judged after
$code OdeErrControl$$ returns.
If $italic ef$$ is to large relative to $italic maxabs$$, 
$code OdeErrControl$$ can be called again 
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
	example/ode_err_control.cpp%
	example/ode_err_maxabs.cpp
%$$
The files
$xref/OdeErrControl.cpp/$$
and
$xref/OdeErrMaxabs.cpp/$$
contain examples and tests of using this routine.
They return true if they succeed and false otherwise.

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
$code cppad/ode_err_control.hpp$$.

$end
--------------------------------------------------------------------------
*/

// link exp and log for float and double
# include <cppad/base_require.hpp>

# include <cppad/local/cppad_assert.hpp>
# include <cppad/check_simple_vector.hpp>
# include <cppad/nan.hpp>

namespace CppAD { // Begin CppAD namespace

template <typename Scalar, typename Vector, typename Method>
Vector OdeErrControl(
	Method          &method, 
	const Scalar    &ti    , 
	const Scalar    &tf    , 
	const Vector    &xi    , 
	const Scalar    &smin  , 
	const Scalar    &smax  , 
	Scalar          &scur  ,
	const Vector    &eabs  , 
	const Scalar    &erel  , 
	Vector          &ef    ,
	Vector          &maxabs,
	size_t          &nstep ) 
{
	// check simple vector class specifications
	CheckSimpleVector<Scalar, Vector>();

	size_t n = xi.size();

	CPPAD_ASSERT_KNOWN(
		smin <= smax,
		"Error in OdeErrControl: smin > smax"
	);
	CPPAD_ASSERT_KNOWN(
		eabs.size() == n,
		"Error in OdeErrControl: size of eabs is not equal to n"
	);
	CPPAD_ASSERT_KNOWN(
		maxabs.size() == n,
		"Error in OdeErrControl: size of maxabs is not equal to n"
	);
	size_t m = method.order();
	CPPAD_ASSERT_KNOWN(
		m > 1,
		"Error in OdeErrControl: m is less than or equal one"
	);

	bool    ok;
	bool    minimum_step;
	size_t  i;
	Vector xa(n), xb(n), eb(n), nan_vec(n);

	// initialization
	Scalar zero(0);
	Scalar one(1);
	Scalar two(2);
	Scalar three(3);
	Scalar m1(m-1);
	Scalar ta = ti;
	for(i = 0; i < n; i++)
	{	nan_vec[i] = nan(zero);
		ef[i]      = zero;
		xa[i]      = xi[i];
		if( zero <= xi[i] )
			maxabs[i] = xi[i];
		else	maxabs[i] = - xi[i];

	}  
	nstep = 0;

	Scalar tb, step, lambda, axbi, a, r, root;
	while( ! (ta == tf) )
	{	// start with value suggested by error criteria
		step = scur;

		// check maximum
		if( smax <= step )
			step = smax;

		// check minimum
		minimum_step = step <= smin;
		if( minimum_step )
			step = smin;

		// check if near the end
		if( tf <= ta + step * three / two )
			tb = tf;
		else	tb = ta + step;

		// try using this step size
		nstep++;
		method.step(ta, tb, xa, xb, eb);
		step = tb - ta;

		// check if this steps error estimate is ok
		ok = ! (hasnan(xb) || hasnan(eb));
		if( (! ok) && minimum_step )
		{	ef = nan_vec;
			return nan_vec;
		}

		// compute value of lambda for this step
		lambda = Scalar(10) * scur / step;
		for(i = 0; i < n; i++)
		{	if( zero <= xb[i] )
				axbi = xb[i];
			else	axbi = - xb[i];
			a    = eabs[i] + erel * axbi;
			if( ! (eb[i] == zero) )
			{	r = ( a / eb[i] ) * step / (tf - ti);
				root = exp( log(r) / m1 ); 
				if( root <= lambda )
					lambda = root;
			}
		}
		if( ok && ( one <= lambda || step <= smin * three / two) )
		{	// this step is within error limits or 
			// close to the minimum size
			ta = tb;
			for(i = 0; i < n; i++)
			{	xa[i] = xb[i];
				ef[i] = ef[i] + eb[i];
				if( zero <= xb[i] )
					axbi = xb[i];
				else	axbi = - xb[i];
				if( axbi > maxabs[i] )
					maxabs[i] = axbi;
			}
		}
		if( ! ok )
		{	// decrease step an see if method will work this time
			scur = step / two;
		}
		else if( ! (ta == tf) )
		{	// step suggested by the error criteria is not used
			// on the last step because it may be very small.
			scur = lambda * step / two;
		}
	}
	return xa;
}

template <typename Scalar, typename Vector, typename Method>
Vector OdeErrControl(
	Method          &method, 
	const Scalar    &ti    , 
	const Scalar    &tf    , 
	const Vector    &xi    , 
	const Scalar    &smin  , 
	const Scalar    &smax  , 
	Scalar          &scur  ,
	const Vector    &eabs  , 
	const Scalar    &erel  , 
	Vector          &ef    ) 
{	Vector maxabs(xi.size());
	size_t nstep;
	return OdeErrControl(
	method, ti, tf, xi, smin, smax, scur, eabs, erel, ef, maxabs, nstep
	);
}

template <typename Scalar, typename Vector, typename Method>
Vector OdeErrControl(
	Method          &method, 
	const Scalar    &ti    , 
	const Scalar    &tf    , 
	const Vector    &xi    , 
	const Scalar    &smin  , 
	const Scalar    &smax  , 
	Scalar          &scur  ,
	const Vector    &eabs  , 
	const Scalar    &erel  , 
	Vector          &ef    ,
	Vector          &maxabs) 
{	size_t nstep;
	return OdeErrControl(
	method, ti, tf, xi, smin, smax, scur, eabs, erel, ef, maxabs, nstep
	);
}

} // End CppAD namespace 

# endif

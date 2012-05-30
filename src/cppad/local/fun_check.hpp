/* $Id: fun_check.hpp 1598 2009-12-03 02:59:52Z bradbell $ */
# ifndef CPPAD_FUN_CHECK_INCLUDED
# define CPPAD_FUN_CHECK_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-09 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin FunCheck$$
$spell
	exp
	bool
	const
	Taylor
$$

$index FunCheck$$
$index check, ADFun$$
$index ADFun, check$$

$section Check an ADFun Sequence of Operations$$

$head Syntax$$
$syntax%%ok% = FunCheck(%f%, %g%, %x%, %r%, %a%)%$$
$pre
$$
$bold See Also$$
$xref/CompareChange/$$


$head Purpose$$
We use $latex F : B^n \rightarrow B^m$$ to denote the
$xref/glossary/AD Function/AD function/$$ corresponding to $italic f$$.
We use $latex G : B^n \rightarrow B^m$$ to denote the
function corresponding to the C++ function object $italic g$$.
This routine check if
$latex \[
	F(x) = G(x)
\]$$
If $latex F(x) \neq G(x)$$, the
$xref/glossary/Operation/Sequence/operation sequence/1/$$
corresponding to $italic f$$ does not represents the algorithm used
by $italic g$$ to calculate values for $latex G$$
(see $xref/FunCheck/Discussion/Discussion/$$ below).

$head f$$
The $code FunCheck$$ argument $italic f$$ has prototype
$syntax%
	ADFun<%Base%> %f%
%$$
Note that the $xref/ADFun/$$ object $italic f$$ is not $code const$$
(see $xref/FunCheck/FunCheck Uses Forward/Forward/$$ below).

$head g$$
The $code FunCheck$$ argument $italic g$$ has prototype
$syntax%
	%Fun% &%g%
%$$
($italic Fun$$ is defined the properties of $italic g$$).
The C++ function object $italic g$$ supports the syntax
$syntax%
	%y% = %g%(%x%)
%$$
which computes $latex y = G(x)$$.

$subhead x$$
The $italic g$$ argument $italic x$$ has prototype
$syntax%
	const %Vector% &%x%
%$$
(see $xref/FunCheck/Vector/Vector/$$ below)
and its size 
must be equal to $italic n$$, the dimension of the
$xref/seq_property/Domain/domain/$$ space for $italic f$$.

$head y$$
The $italic g$$ result $italic y$$ has prototype
$syntax%
	%Vector% %y%
%$$
and its value is $latex G(x)$$.
The size of $italic y$$ 
is equal to $italic m$$, the dimension of the
$xref/seq_property/Range/range/$$ space for $italic f$$.

$head x$$
The $code FunCheck$$ argument $italic x$$ has prototype
$syntax%
	const %Vector% &%x%
%$$
and its size 
must be equal to $italic n$$, the dimension of the
$xref/seq_property/Domain/domain/$$ space for $italic f$$.
This specifies that point at which to compare the values
calculated by $italic f$$ and $italic G$$. 

$head r$$
The $code FunCheck$$ argument $italic r$$ has prototype
$syntax%
	const %Base% &%r%
%$$
It specifies the relative error the element by element
comparison of the value of $latex F(x)$$ and $latex G(x)$$.

$head a$$
The $code FunCheck$$ argument $italic a$$ has prototype
$syntax%
	const %Base% &%a%
%$$
It specifies the absolute error the element by element
comparison of the value of $latex F(x)$$ and $latex G(x)$$.

$head ok$$
The $code FunCheck$$ result $italic ok$$ has prototype
$syntax%
	bool %ok%
%$$
It is true, if for $latex i = 0 , \ldots , m-1$$ 
either the relative error bound is satisfied
$latex \[
| F_i (x) - G_i (x) | 
\leq 
r ( | F_i (x) | + | G_i (x) | ) 
\] $$
or the absolute error bound is satisfied
$latex \[
	| F_i (x) - G_i (x) | \leq a
\] $$
It is false if for some $latex (i, j)$$ neither
of these bounds is satisfied.

$head Vector$$
The type $italic Vector$$ must be a $xref/SimpleVector/$$ class with
$xref/SimpleVector/Elements of Specified Type/elements of type/$$
$italic Base$$.
The routine $xref/CheckSimpleVector/$$ will generate an error message
if this is not the case.

$head FunCheck Uses Forward$$
After each call to $xref/Forward/$$,
the object $italic f$$ contains the corresponding 
$xref/glossary/Taylor Coefficient/Taylor coefficients/$$.
After $code FunCheck$$,
the previous calls to $xref/Forward/$$ are undefined.

$head Discussion$$
Suppose that the algorithm corresponding to $italic g$$ contains
$syntax%
	if( %x% >= 0 )
		%y% = exp(%x%)
	else	%y% = exp(-%x%)
%$$ 
where $italic x$$ and $italic y$$ are $syntax%AD<double>%$$ objects.
It follows that the 
AD of $code double$$ $xref/glossary/Operation/Sequence/operation sequence/1/$$
depends on the value of $italic x$$.
If the sequence of operations stored in $italic f$$ corresponds to 
$italic g$$ with $latex x \geq 0$$, 
the function values computed using $italic f$$ when $latex x < 0$$
will not agree with the function values computed by $latex g$$.
This is because the operation sequence corresponding to $italic g$$ changed
(and hence the object $italic f$$ does not represent the function
$latex G$$ for this value of $italic x$$).
In this case, you probably want to re-tape the calculations
performed by $italic g$$ with the
$cref/independent variables/glossary/Tape/Independent Variable/$$ 
equal to the values in $italic x$$ 
(so AD operation sequence properly represents the algorithm
for this value of independent variables).


$head Example$$
$children%
	example/fun_check.cpp
%$$
The file
$xref/FunCheck.cpp/$$
contains an example and test of this function.   
It returns true if it succeeds and false otherwise.

$end
---------------------------------------------------------------------------
*/

namespace CppAD {
	template <class Base, class Fun, class Vector>
	bool FunCheck(
		ADFun<Base>  &f , 
		Fun          &g , 
		const Vector &x , 
		const Base   &r ,
		const Base   &a )
	{	bool ok = true;
	
		size_t m   = f.Range();
		Vector yf  = f.Forward(0, x); 
		Vector yg  = g(x);
	
		size_t i;
		for(i = 0; i < m; i++)
			ok  &= NearEqual(yf[i], yg[i], r, a);
		return ok;
	}
}

# endif

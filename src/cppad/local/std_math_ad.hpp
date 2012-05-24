/* $Id: std_math_ad.hpp 2057 2011-08-11 14:07:11Z bradbell $ */
# ifndef CPPAD_STD_MATH_AD_INCLUDED
# define CPPAD_STD_MATH_AD_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
/*
-------------------------------------------------------------------------------
$begin std_math_ad$$
$spell
	Vec
	std
	atan
	const
	acos
	asin
	atan
	cos
	exp
	sqrt
	CppAD
	namespace
	tanh
$$

$index standard, AD math unary$$
$index math, AD unary$$
$index unary, AD math$$

$index acos, AD$$
$index asin, AD$$
$index atan, AD$$
$index cos, AD$$
$index cosh, AD$$
$index exp, AD$$
$index log, AD$$
$index log10, AD$$
$index sin, AD$$
$index sinh, AD$$
$index sqrt, AD$$
$index tan, AD$$
$index tanh, AD$$

$section AD Standard Math Unary Functions$$

$head Syntax$$
$syntax%%y% = %fun%(%x%)%$$


$head Purpose$$
Evaluates the one argument standard math function 
$italic fun$$ where its argument is an 
$xref/glossary/AD of Base/AD of/$$ $italic Base$$ object.

$head x$$
The argument $italic x$$ has one of the following prototypes
$syntax%
	const AD<%Base%>               &%x%
	const VecAD<%Base%>::reference &%x%
%$$

$head y$$
The result $italic y$$ has prototype
$syntax%
	AD<%Base%> %y%
%$$


$head Operation Sequence$$
Most of these functions are AD of $italic Base$$
$xref/glossary/Operation/Atomic/atomic operations/1/$$.
In all cases,
The AD of $italic Base$$
operation sequence used to calculate $italic y$$ is 
$xref/glossary/Operation/Independent/independent/1/$$
of $italic x$$.

$head fun$$ 
A definition of $italic fun$$ is included 
for each of the following functions:
$code acos$$,
$code asin$$,
$code atan$$,
$code cos$$,
$code cosh$$,
$code exp$$,
$code log$$,
$code log10$$,
$code sin$$,
$code sinh$$,
$code sqrt$$,
$code tan$$,
$code tanh$$.


$head Examples$$
The following files
contain examples and tests of these functions.   
Each test returns true if it succeeds and false otherwise.
$children%
	example/acos.cpp%
	example/asin.cpp%
	example/atan.cpp%
	example/cos.cpp%
	example/cosh.cpp%
	example/exp.cpp%
	example/log.cpp%
	example/log_10.cpp%
	example/sin.cpp%
	example/sinh.cpp%
	example/sqrt.cpp%
	example/tan.cpp%
	example/tanh.cpp
%$$
$table
$rref Acos.cpp$$
$rref Asin.cpp$$
$rref Atan.cpp$$
$rref Cos.cpp$$
$rref Cosh.cpp$$
$rref Exp.cpp$$
$rref Log.cpp$$
$rref Log10.cpp$$
$rref Sin.cpp$$
$rref Sinh.cpp$$
$rref Sqrt.cpp$$
$rref Tan.cpp$$
$rref Tanh.cpp$$
$tend


$head Derivatives$$
Each of these functions satisfy a standard math function differential equation.
Calculating derivatives using this differential equation 
is discussed for 
both $xref/ForwardTheory/Standard Math Functions/forward/$$
and $xref/ReverseTheory/Standard Math Functions/reverse/$$ mode.
The exact form of the differential equation
for each of these functions is listed below:

$subhead acos$$
$latex \[
\begin{array}{lcr}
	\D{[ {\rm acos} (x) ]}{x} & = & - (1 - x * x)^{-1/2}
\end{array}
\] $$

$subhead asin$$
$latex \[
\begin{array}{lcr}
	\D{[ {\rm asin} (x) ]}{x} & = & (1 - x * x)^{-1/2}
\end{array}
\] $$

$subhead atan$$
$latex \[
\begin{array}{lcr}
        \D{[ {\rm atan} (x) ]}{x} & = & \frac{1}{1 + x^2}
\end{array}
\] $$

$subhead cos$$
$latex \[
\begin{array}{lcr}
        \D{[ \cos (x) ]}{x} & = & - \sin (x)  \\
        \D{[ \sin (x) ]}{x} & = & \cos (x)
\end{array}
\] $$

$subhead cosh$$
$latex \[
\begin{array}{lcr}
        \D{[ \cosh (x) ]}{x} & = & \sinh (x)  \\
        \D{[ \sin (x) ]}{x}  & = & \cosh (x)
\end{array}
\] $$

$subhead exp$$
$latex \[
\begin{array}{lcr}
        \D{[ \exp (x) ]}{x} & = & \exp (x)
\end{array}
\] $$

$subhead log$$
$latex \[
\begin{array}{lcr}
        \D{[ \log (x) ]}{x} & = & \frac{1}{x}
\end{array}
\] $$

$subhead log10$$
This function is special in that it's derivatives are calculated
using the relation
$latex \[
\begin{array}{lcr}
        {\rm log10} (x) & = & \log(x) / \log(10)
\end{array}
\] $$

$subhead sin$$
$latex \[
\begin{array}{lcr}
        \D{[ \sin (x) ]}{x} & = & \cos (x) \\
        \D{[ \cos (x) ]}{x} & = & - \sin (x) 
\end{array}
\] $$

$subhead sinh$$
$latex \[
\begin{array}{lcr}
        \D{[ \sinh (x) ]}{x} & = & \cosh (x)   \\
        \D{[ \cosh (x) ]}{x} & = & \sinh (x)
\end{array}
\] $$

$subhead sqrt$$
$latex \[
\begin{array}{lcr}
        \D{[ {\rm sqrt} (x) ]}{x} & = & \frac{1}{2 {\rm sqrt} (x) }
\end{array}
\] $$

$subhead tan$$
$latex \[
\begin{array}{lcr}
        \D{[ \tan (x) ]}{x} & = & 1 + \tan (x)^2
\end{array}
\] $$

$subhead tanh$$
$latex \[
\begin{array}{lcr}
        \D{[ \tanh (x) ]}{x} & = & 1 - \tanh (x)^2
\end{array}
\] $$

$end
-------------------------------------------------------------------------------
*/

# define CPPAD_STANDARD_MATH_UNARY_AD(Name, Op)                   \
    template <class Base>                                         \
    inline AD<Base> Name(const AD<Base> &x)                       \
    {   return x.Name(); }                                        \
    template <class Base>                                         \
    inline AD<Base> AD<Base>::Name (void) const                   \
    {                                                             \
        AD<Base> result;                                          \
        result.value_ = CppAD::Name(value_);                      \
        CPPAD_ASSERT_UNKNOWN( Parameter(result) );                \
                                                                  \
        if( Variable(*this) )                                     \
        {   CPPAD_ASSERT_UNKNOWN( NumRes(Op) <= 2 );              \
            CPPAD_ASSERT_UNKNOWN( NumArg(Op) == 1 );              \
            ADTape<Base> *tape = tape_this();                     \
            tape->Rec_.PutArg(taddr_);                            \
            result.taddr_ = tape->Rec_.PutOp(Op);                 \
            result.id_    = tape->id_;                            \
        }                                                         \
        return result;                                            \
    }                                                             \
    template <class Base>                                         \
    inline AD<Base> Name(const VecAD_reference<Base> &x)          \
    {   return Name( x.ADBase() ); }

//  BEGIN CppAD namespace
namespace CppAD {

     CPPAD_STANDARD_MATH_UNARY_AD(acos, AcosOp)
     CPPAD_STANDARD_MATH_UNARY_AD(asin, AsinOp)
     CPPAD_STANDARD_MATH_UNARY_AD(atan, AtanOp)
     CPPAD_STANDARD_MATH_UNARY_AD(cos, CosOp)
     CPPAD_STANDARD_MATH_UNARY_AD(cosh, CoshOp)
     CPPAD_STANDARD_MATH_UNARY_AD(exp, ExpOp)
     CPPAD_STANDARD_MATH_UNARY_AD(log, LogOp)
     CPPAD_STANDARD_MATH_UNARY_AD(sin, SinOp)
     CPPAD_STANDARD_MATH_UNARY_AD(sinh, SinhOp)
     CPPAD_STANDARD_MATH_UNARY_AD(sqrt, SqrtOp)
     CPPAD_STANDARD_MATH_UNARY_AD(tan, TanOp)
     CPPAD_STANDARD_MATH_UNARY_AD(tanh, TanhOp)

     // log10
     template <class Base>
     inline AD<Base> log10(const AD<Base> &x)
     {    return CppAD::log(x) / CppAD::log( Base(10) ); }
     template <class Base>
     inline AD<Base> log10(const VecAD_reference<Base> &x)
     {    return CppAD::log(x.ADBase()) / CppAD::log( Base(10) ); }
}

# undef CPPAD_STANDARD_MATH_UNARY_AD

# endif 

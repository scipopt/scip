/* $Id: user_atomic.hpp 2210 2011-11-24 21:56:53Z bradbell $ */
# ifndef CPPAD_USER_ATOMIC_INCLUDED
# define CPPAD_USER_ATOMIC_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin user_atomic$$
$spell
	hes
	std
	Jacobian
	jac
	Tvector
	afun
	vx
	vy
	bool
	namespace
	CppAD
	const
	Taylor
	tx
	ty
	px
	py
$$

$section User Defined Atomic AD Functions$$
$index CPPAD_USER_ATOMIC$$
$index atomic, user function$$
$index user, atomic function$$
$index operation, user atomic$$
$index function, user atomic$$

$head Syntax$$

$subhead Define Function$$
$codei%CPPAD_USER_ATOMIC(%afun%, %Tvector%, %Base%, 
	%forward%, %reverse%, %for_jac_sparse%, %rev_jac_sparse%, %rev_hes_sparse%
)
%$$

$subhead Use Function$$
$icode%afun%(%id%, %ax%, %ay%)
%$$

$subhead Callback Routines$$
$icode%ok% = %forward%(%id%, %k%, %n%, %m%, %vx%, %vy%, %tx%, %ty%)
%$$
$icode%ok% = %reverse%(%id%, %k%, %n%, %m%, %tx%, %ty%, %px%, %py%)
%$$
$icode%ok% = %for_jac_sparse%(%id%, %n%, %m%, %q%, %r%, %s%)
%$$
$icode%ok% = %rev_jac_sparse%(%id%, %n%, %m%, %q%, %r%, %s%)
%$$
$icode%ok% = %rev_hes_sparse%(%id%, %n%, %m%, %q%, %r%, %s%, %t%, %u%, %v%)
%$$

$subhead Free Static Memory$$
$codei%user_atomic<%Base%>::clear()%$$

$head Purpose$$
In some cases, the user knows how to compute the derivative
of a function 
$latex \[
	y = f(x) \; {\rm where} \; f : B^n \rightarrow B^m 
\] $$
more efficiently than by coding it $codei%AD<%Base%>%$$ 
$cref/atomic/glossary/Operation/Atomic/$$ operations
and letting CppAD do the rest.
In this case, $code CPPAD_USER_ATOMIC$$ can be used
add the user code for $latex f(x)$$, and its derivatives,
to the set of $codei%AD<%Base%>%$$ atomic operations. 

$head Partial Implementation$$
The routines 
$cref/forward/user_atomic/forward/$$,
$cref/reverse/user_atomic/reverse/$$,
$cref/for_jac_sparse/user_atomic/for_jac_sparse/$$,
$cref/rev_jac_sparse/user_atomic/rev_jac_sparse/$$, and
$cref/rev_hes_sparse/user_atomic/rev_hes_sparse/$$,
must be defined by the user.
The $icode forward$$ the routine, 
for the case $icode%k% = 0%$$, must be implemented.
Functions with the correct prototype,
that just return $icode false$$, 
can be used for the other cases 
(unless they are required by your calculations). 
For example, you need not implement
$icode forward$$ for the case $icode%k% == 2%$$ until you require
forward mode calculation of second derivatives.

$head CPPAD_USER_ATOMIC$$
$index CPPAD_USER_ATOMIC$$
The macro 
$codei%
CPPAD_USER_ATOMIC(%afun%, %Tvector%, %Base%, 
	%forward%, %reverse%, %for_jac_sparse%, %rev_jac_sparse%, %rev_hes_sparse%
)
%$$ 
defines the $codei%AD<%Base%>%$$ routine $icode afun$$.
This macro can be placed within a namespace 
(not the $code CppAD$$ namespace) 
but must be outside of any routine.

$subhead Tvector$$
The macro argument $icode Tvector$$ must be a
$cref/simple vector template class/SimpleVector/$$.
It determines the type of vectors used as arguments to the routine 
$icode afun$$.

$subhead Base$$
The macro argument $icode Base$$ specifies the 
$cref/base type/base_require/$$
corresponding to $codei%AD<%Base>%$$ operation sequences.
Calling the routine $icode afun$$ will add the operator defined
by this macro to an $codei%AD<%Base>%$$ operation sequence.

$head ok$$
For all routines documented below,
the return value $icode ok$$ has prototype
$codei%
	bool %ok%
%$$
If it is $code true$$, the corresponding evaluation succeeded,
otherwise it failed.

$head id$$
For all routines documented below,
the argument $icode id$$ has prototype
$codei%
	size_t %id%
%$$
Its value in all other calls is the same as in the corresponding
call to $icode afun$$.
It can be used to store and retrieve extra information about
a specific call to $icode afun$$. 

$head k$$
For all routines documented below, the argument $icode k$$ has prototype
$codei%
	size_t %k%
%$$
The value $icode%k%$$ is the order of the Taylor coefficient that
we are evaluating ($cref/forward/user_atomic/forward/$$)
or taking the derivative of ($cref/reverse/user_atomic/reverse/$$).

$head n$$
For all routines documented below, 
the argument $icode n$$ has prototype
$codei%
	size_t %n%
%$$
It is the size of the vector $icode ax$$ in the corresponding call to
$icode%afun%(%id%, %ax%, %ay%)%$$; i.e.,
the dimension of the domain space for $latex y = f(x)$$.

$head m$$
For all routines documented below, the argument $icode m$$ has prototype
$codei%
	size_t %m%
%$$
It is the size of the vector $icode ay$$ in the corresponding call to
$icode%afun%(%id%, %ax%, %ay%)%$$; i.e.,
the dimension of the range space for $latex y = f(x)$$.

$head tx$$
For all routines documented below, 
the argument $icode tx$$ has prototype
$codei%
	const CppAD::vector<%Base%>& %tx%
%$$
and $icode%tx%.size() >= (%k% + 1) * %n%$$.
For $latex j = 0 , \ldots , n-1$$ and $latex \ell = 0 , \ldots , k$$,
we use the Taylor coefficient notation
$latex \[
\begin{array}{rcl}
	x_j^\ell & = & tx [ j * ( k + 1 ) + \ell ]
	\\
	X_j (t) & = & x_j^0 + x_j^1 t^1 + \cdots + x_j^k t^k
\end{array}
\] $$
If $icode%tx%.size() > (%k% + 1) * %n%$$,
the other components of $icode tx$$ are not specified and should not be used.
Note that superscripts represent an index for $latex x_j^\ell$$
and an exponent for $latex t^\ell$$.
Also note that the Taylor coefficients for $latex X(t)$$ correspond
to the derivatives of $latex X(t)$$ at $latex t = 0$$ in the following way:
$latex \[
	x_j^\ell = \frac{1}{ \ell ! } X_j^{(\ell)} (0)
\] $$

$head ty$$
In calls to $cref/forward/user_atomic/forward/$$, 
the argument $icode ty$$ has prototype
$codei%
	CppAD::vector<%Base%>& %ty%
%$$
while in calls to $cref/reverse/user_atomic/reverse/$$ it has prototype
$codei%
	const CppAD::vector<%Base%>& %ty%
%$$
For all calls, $icode%tx%.size() >= (%k% + 1) * %m%$$.
For $latex i = 0 , \ldots , m-1$$ and $latex \ell = 0 , \ldots , k$$,
we use the Taylor coefficient notation
$latex \[
\begin{array}{rcl}
	y_i^\ell & = & ty [ i * ( k + 1 ) + \ell ]
	\\
	Y_i (t)  & = & y_i^0 + y_i^1 t^1 + \cdots + y_i^k t^k + o ( t^k )
\end{array}
\] $$
where $latex o( t^k ) / t^k \rightarrow 0$$ as $latex t \rightarrow 0$$.
If $icode%ty%.size() > (%k% + 1) * %m%$$,
the other components of $icode ty$$ are not specified and should not be used.
Note that superscripts represent an index for $latex y_j^\ell$$
and an exponent for $latex t^\ell$$.
Also note that the Taylor coefficients for $latex Y(t)$$ correspond
to the derivatives of $latex Y(t)$$ at $latex t = 0$$ in the following way:
$latex \[
	y_j^\ell = \frac{1}{ \ell ! } Y_j^{(\ell)} (0)
\] $$

$head afun$$
The macro argument $icode afun$$,
is the name of the AD function corresponding to this atomic
operation (as it is used in the source code).
CppAD uses the other functions,
where the arguments are vectors with elements of type $icode Base$$,
to implement the function 
$codei%
	%afun%(%id%, %ax%, %ay%)
%$$
where the argument are vectors with elements of type $codei%AD<%Base%>%$$.

$subhead ax$$
The $icode afun$$ argument $icode ax$$ has prototype
$codei%
	const %Tvector%< AD<%Base%> >& %ax%
%$$
It is the argument vector $latex x \in B^n$$ 
at which the $codei%AD<%Base%>%$$ version of 
$latex y = f(x)$$ is to be evaluated.
The dimension of the domain space for $latex y = f (x)$$
is specified by $cref/n/user_atomic/n/$$ $codei%= %ax%.size()%$$,
which must be greater than zero.

$subhead ay$$
The $icode afun$$ result $icode ay$$ has prototype
$codei%
	%Tvector%< AD<%Base%> >& %ay%
%$$
The input values of its elements do not matter.
Upon return, it is the $codei%AD<%Base%>%$$ version of the 
result vector $latex y = f(x)$$.
The dimension of the range space for $latex y = f (x)$$
is specified by $cref/m/user_atomic/m/$$ $codei%= %ay%.size()%$$,
which must be greater than zero.

$subhead Parallel Mode$$
$index parallel, user_atomic$$
$index user_atomic, parallel$$
The first call to 
$codei%
	%afun%(%id%, %ax%, %ay%)
%$$
must not be in $cref/parallel/ta_in_parallel/$$ mode.

$head forward$$
The macro argument $icode forward$$ is a
user defined function
$codei%
	%ok% = %forward%(%id%, %k%, %n%, %m%, %vx%, %vy%, %tx%, %ty%)
%$$
that computes results during a $cref/forward/Forward/$$ mode sweep.
For this call, we are given the Taylor coefficients in $icode tx$$ 
form order zero through $icode k$$,
and the Taylor coefficients in  $icode ty$$ with order less than $icode k$$.
The $icode forward$$ routine computes the 
$icode k$$ order Taylor coefficients for $latex y$$ using the definition
$latex Y(t) = f[ X(t) ]$$. 
For example, for $latex i = 0 , \ldots , m-1$$,
$latex \[
\begin{array}{rcl}
y_i^0 & = & Y(0) 
        = f_i ( x^0 )
\\
y_i^1 & = & Y^{(1)} ( 0 ) 
        = f_i^{(1)} ( x^0 ) X^{(1)} ( 0 ) 
        =  f_i^{(1)} ( x^0 ) x^1 
\end{array}
\] $$
Then, for $latex i = 0 , \ldots , m-1$$, it sets
$latex \[
	ty [ i * (k + 1) + k ] = y_i^k
\] $$
The other components of $icode ty$$ must be left unchanged.

$subhead Usage$$
This routine is used,
with $icode%vx%.size() > 0%$$ and $icode%k% == 0%$$,
by calls to $icode afun$$.
It is used,
with $icode%vx%.size() = 0%$$ and
$icode k$$ equal to the order of the derivative begin computed,
by calls to $cref/forward/ForwardAny/$$.

$subhead vx$$
The $icode forward$$ argument $icode vx$$ has prototype
$codei%
	const CppAD::vector<bool>& %vx%
%$$
The case $icode%vx%.size() > 0%$$ occurs 
once for each call to $icode afun$$,
during the call,
and before any of the other callbacks corresponding to that call.
Hence such a call can be used to cache information attached to 
the corresponding $icode id$$
(such as the elements of $icode vx$$).
If $icode%vx%.size() > 0%$$ then
$icode%k% == 0%$$, 
$icode%vx%.size() >= %n%$$, and
for $latex j = 0 , \ldots , n-1$$,
$icode%vx%[%j%]%$$ is true if and only if
$icode%ax%[%j%]%$$ is a $cref/variable/glossary/Variable/$$.
$pre

$$
If $icode%vx%.size() == 0%$$, 
then $icode%vy%.size() == 0%$$ and neither of these vectors
should be used.

$subhead vy$$
The $icode forward$$ argument $icode vy$$ has prototype
$codei%
	CppAD::vector<bool>& %vy%
%$$
If $icode%vy%.size() == 0%$$, it should not be used.
Otherwise, 
$icode%k% == 0%$$ and $icode%vy%.size() >= %m%$$.
The input values of the elements of $icode vy$$ do not matter.
Upon return, for $latex j = 0 , \ldots , m-1$$,
$icode%vy%[%i%]%$$ is true if and only if
$icode%ay%[%j%]%$$ is a variable.
(CppAD uses $icode vy$$ to reduce the necessary computations.)

$head reverse$$
The macro argument $icode reverse$$
is a user defined function
$codei%
	%ok% = %reverse%(%id%, %k%, %n%, %m%, %tx%, %ty%, %px%, %py%)
%$$
that computes results during a $cref/reverse/Reverse/$$ mode sweep. 
The input value of the vectors $icode tx$$ and $icode ty$$
contain Taylor coefficient, up to order $icode k$$,
for $latex X(t)$$ and $latex Y(t)$$ respectively.
We use the $latex \{ x_j^\ell \}$$ and $latex \{ y_i^\ell \}$$
to denote these Taylor coefficients where the implicit range indices are
$latex i = 0 , \ldots , m-1$$,
$latex j = 0 , \ldots , n-1$$,
$latex \ell = 0 , \ldots , k$$.
Using the calculations done by $cref/forward/user_atomic/forward/$$,
the Taylor coefficients $latex \{ y_i^\ell \}$$ are a function of the Taylor
coefficients for $latex \{ x_j^\ell \}$$; i.e., given $latex y = f(x)$$
we define the function
$latex F : B^{n \times (k+1)} \rightarrow B^{m \times (k+1)}$$ by
$latex \[
y_i^\ell =  F_i^\ell ( \{ x_j^\ell \} )
\] $$
We use $latex G : B^{m \times (k+1)} \rightarrow B$$
to denote an arbitrary scalar valued function of the Taylor coefficients for 
$latex Y(t)$$ and write  $latex z = G( \{ y_i^\ell \} )$$.
The $code reverse$$ routine
is given the derivative of $latex z$$ with respect to
$latex \{ y_i^\ell \}$$ and computes its derivative with respect
to $latex \{ x_j^\ell \}$$.

$subhead Usage$$
This routine is used,
with $icode%k% + 1%$$ equal to the order of the derivative being calculated,
by calls to $cref/reverse/reverse_any/$$.

$subhead py$$
The $icode reverse$$ argument $icode py$$ has prototype
$codei%
	const CppAD::vector<%Base%>& %py%
%$$
and $icode%py%.size() >= (%k% + 1) * %m%$$.
For $latex i = 0 , \ldots , m-1$$ and $latex \ell = 0 , \ldots , k$$,
$latex \[
	py[ i * (k + 1 ) + \ell ] = \partial G / \partial y_i^\ell
\] $$
If $icode%py%.size() > (%k% + 1) * %m%$$,
the other components of $icode py$$ are not specified and should not be used.

$subhead px$$
We define the function 
$latex \[
H ( \{ x_j^\ell \} ) = G[ F( \{ x_j^\ell \} ) ]
\] $$
The $icode reverse$$ argument $icode px$$ has prototype
$codei%
	CppAD::vector<%Base%>& %px%
%$$
and $icode%px%.size() >= (%k% + 1) * %n%$$.
The input values of the elements of $icode px$$ do not matter.
Upon return,
for $latex j = 0 , \ldots , n-1$$ and $latex \ell = 0 , \ldots , k$$,
$latex \[
\begin{array}{rcl}
px [ j * (k + 1) + \ell ] & = & \partial H / \partial x_j^\ell
\\
& = & 
( \partial G / \partial \{ y_i^\ell \} ) 
	( \partial \{ y_i^\ell \} / \partial x_j^\ell )
\\
& = & 
\sum_{i=0}^{m-1} \sum_{\ell=0}^k
( \partial G / \partial y_i^\ell ) ( \partial y_i^\ell / \partial x_j^\ell )
\\
& = & 
\sum_{i=0}^{m-1} \sum_{\ell=0}^k
py[ i * (k + 1 ) + \ell ] ( \partial F_i^\ell / \partial x_j^\ell )
\end{array}
\] $$
If $icode%px%.size() > (%k% + 1) * %n%$$,
the other components of $icode px$$ are not specified and should not be used.

$head for_jac_sparse$$
The macro argument $icode for_jac_sparse$$
is a user defined function
$codei%
	%ok% = %for_jac_sparse%(%id%, %n%, %m%, %q%, %r%, %s%)
%$$
that is used to compute results during a forward Jacobian sparsity sweep.
For a fixed $latex n \times q$$ matrix $latex R$$,
the Jacobian of $latex f( x + R * u)$$ with respect to $latex u \in B^q$$ is
$latex \[
	S(x) = f^{(1)} (x) * R
\] $$
Given a $cref/sparsity pattern/glossary/Sparsity Pattern/$$ for $latex R$$,
$icode for_jac_sparse$$ computes a sparsity pattern for $latex S(x)$$.

$subhead Usage$$
This routine is used by calls to $cref/ForSparseJac/$$.

$subhead q$$
The $icode for_jac_sparse$$ argument $icode q$$ has prototype
$codei%
     size_t %q%
%$$
It specifies the number of columns in 
$latex R \in B^{n \times q}$$ and the Jacobian 
$latex S(x) \in B^{m \times q}$$. 

$subhead r$$
The $icode for_jac_sparse$$ argument $icode r$$ has prototype
$codei%
     const CppAD::vector< std::set<size_t> >& %r%
%$$
and $icode%r%.size() >= %n%$$.
For $latex j = 0 , \ldots , n-1$$,
all the elements of $icode%r%[%j%]%$$ are between
zero and $icode%q%-1%$$ inclusive.
This specifies a sparsity pattern for the matrix $latex R$$.

$subhead s$$
The $icode for_jac_sparse$$ return value $icode s$$ has prototype
$codei%
	CppAD::vector< std::set<size_t> >& %s%
%$$
and $icode%s%.size() >= %m%%$$.
The input values of its sets do not matter. Upon return 
for $latex i = 0 , \ldots , m-1$$,
all the elements of $icode%s%[%i%]%$$ are between
zero and $icode%q%-1%$$ inclusive.
This represents a sparsity pattern for the matrix $latex S(x)$$.

$head rev_jac_sparse$$
The macro argument $icode rev_jac_sparse$$
is a user defined function
$codei%
	%ok% = %rev_jac_sparse%(%id%, %n%, %m%, %q%, %r%, %s%)
%$$
that is used to compute results during a reverse Jacobian sparsity sweep.
For a fixed $latex q \times m$$ matrix $latex S$$,
the Jacobian of $latex S * f( x )$$ with respect to $latex x \in B^n$$ is
$latex \[
	R(x) = S * f^{(1)} (x)
\] $$
Given a $cref/sparsity pattern/glossary/Sparsity Pattern/$$ for $latex S$$,
$icode rev_jac_sparse$$ computes a sparsity pattern for $latex R(x)$$.

$subhead Usage$$
This routine is used by calls to $cref/RevSparseJac/$$
and to $cref/optimize/$$.


$subhead q$$
The $icode rev_jac_sparse$$ argument $icode q$$ has prototype
$codei%
     size_t %q%
%$$
It specifies the number of rows in 
$latex S \in B^{q \times m}$$ and the Jacobian 
$latex R(x) \in B^{q \times n}$$. 

$subhead s$$
The $icode rev_jac_sparse$$ argument $icode s$$ has prototype
$codei%
     const CppAD::vector< std::set<size_t> >& %s%
%$$
and $icode%s%.size() >= %m%$$.
For $latex i = 0 , \ldots , m-1$$, 
all the elements of $icode%s%[%i%}%$$
are between zero and $icode%q%-1%$$ inclusive.
This specifies a sparsity pattern for the matrix $latex S^\T$$.

$subhead r$$
The $icode rev_jac_sparse$$ return value $icode r$$ has prototype
$codei%
	CppAD::vector< std::set<size_t> >& %r%
%$$
and $icode%r%.size() >= %n%$$.
The input values of its sets do not matter.
Upon return for $latex j = 0 , \ldots , n-1$$,
all the elements of $icode%r%[%j%]%$$
are between zero and $icode%q%-1%$$ inclusive.
This represents a sparsity pattern for the matrix $latex R(x)^\T$$.

$head rev_hes_sparse$$
The macro argument $icode rev_hes_sparse$$
is a user defined function
$codei%
	%ok% = %rev_hes_sparse%(%id%, %n%, %m%, %q%, %r%, %s%, %t%, %u%, %v%)
%$$
There is an unspecified scalar valued function 
$latex g : B^m \rightarrow B$$.
Given a sparsity pattern for $latex R$$
and information about the function $latex z = g(y)$$,
this routine computes the sparsity pattern for
$latex \[
	V(x) = (g \circ f)^{(2)}( x ) R
\] $$

$subhead Usage$$
This routine is used by calls to $cref/RevSparseHes/$$.

$subhead q$$
The $icode rev_hes_sparse$$ argument $icode q$$ has prototype
$codei%
     size_t %q%
%$$
It specifies the number of columns in the sparsity patterns.

$subhead r$$
The $icode rev_hes_sparse$$ argument $icode r$$ has prototype
$codei%
     const CppAD::vector< std::set<size_t> >& %r%
%$$
and $icode%r%.size() >= %n%$$.
For $latex j = 0 , \ldots , n-1$$,
all the elements of $icode%r%[%j%]%$$ are between
zero and $icode%q%-1%$$ inclusive.
This specifies a sparsity pattern for the matrix $latex R \in B^{n \times q}$$.

$subhead s$$
The $icode rev_hes_sparse$$ argument $icode s$$ has prototype
$codei%
     const CppAD::vector<bool> >& %s%
%$$
and $icode%s%.size() >= %m%$$.
This specifies a sparsity pattern for the matrix 
$latex S(x) = g^{(1)} (y) \in B^{1 \times m}$$.

$subhead t$$
The $icode rev_hes_sparse$$ argument $icode t$$ has prototype
$codei%
     CppAD::vector<bool> >& %t%
%$$
and $icode%t%.size() >= %n%$$.
The input values of its elements do not matter.
Upon return it represents a sparsity pattern for the matrix 
$latex T(x) = (g \circ f)^{(1)} (x) \in B^{1 \times n}$$.

$subhead u$$
The $icode rev_hes_sparse$$ argument $icode u$$ has prototype
$codei%
     const CppAD::vector< std::set<size_t> >& %u%
%$$
and $icode%u%.size() >= %m%$$.
For $latex i = 0 , \ldots , m-1$$,
all the elements of $icode%u%[%i%]%$$
are between zero and $icode%q%-1%$$ inclusive.
This specifies a sparsity pattern
for the matrix $latex U(x) \in B^{m \times q}$$ defined by
$latex \[
\begin{array}{rcl}
U(x) 
& = & 
\partial_u \{ \partial_y g[ y + f^{(1)} (x) R u ] \}_{u=0}
\\
& = & 
\partial_u \{ g^{(1)} [ y + f^{(1)} (x) R u ] \}_{u=0}
\\
& = &
g^{(2)} (y) f^{(1)} (x) R
\end{array}
\] $$

$subhead v$$
The $icode rev_hes_sparse$$ argument $icode v$$ has prototype
$codei%
     CppAD::vector< std::set<size_t> >& %v%
%$$
and $icode%v%.size() >= %n%$$.
The input values of its elements do not matter.
Upon return, for $latex j = 0, \ldots , n-1$$,
all the elements of $icode%v%[%j%]%$$
are between zero and $icode%q%-1%$$ inclusive.
This represents a sparsity pattern for the matrix 
$latex V(x) \in B^{n \times q}$$ defined by
$latex \[
\begin{array}{rcl}
V(x) 
& = & 
\partial_u [ \partial_x (g \circ f) ( x + R u )  ]_{u=0}
\\
& = &
\partial_u [ (g \circ f)^{(1)}( x + R u )  ]_{u=0}
\\
& = &
(g \circ f)^{(2)}( x ) R
\end{array}
\] $$

$head clear$$
User atomic functions hold onto static work space in order to
increase speed by avoiding system memory allocation calls.
The function call $codei%
	user_atomic<%Base%>::clear()
%$$ 
makes to work space $cref/available/ta_available/$$ to
for other uses by the same thread.
This should be called when you are done using the 
user atomic functions for a specific value of $icode Base$$.

$children%
	example/user_tan.cpp%
	example/mat_mul.cpp
%$$
$head Example$$

$subhead Tangent Function$$
The file $cref/user_tan.cpp/$$ contains an example and test
implementation of the tangent function as a user atomic operation.

$subhead Matrix Multiplication$$
The file  $cref/mat_mul.cpp/$$ contains an example and test
implementation of matrix multiplication a a user atomic operation.

$end
------------------------------------------------------------------------------
*/
# include <set>
# include <cppad/local/cppad_assert.hpp>

// needed before one can use CPPAD_ASSERT_FIRST_CALL_NOT_PARALLEL
# include <cppad/thread_alloc.hpp>

CPPAD_BEGIN_NAMESPACE
/*!
\file user_atomic.hpp
user defined atomic operations.
*/

/*!
\def CPPAD_USER_ATOMIC(afun, Tvector, 
	forward, reverse, for_jac_sparse, rev_jac_sparse, rev_hes_sparse 
)
Defines the function <tt>afun(id, ax, ay)</tt>  
where \c id is \c ax and \c ay are vectors with <tt>AD<Base></tt> elements.

\par Tvector
the Simple Vector template class for this function.

\par Base
the base type for the atomic operation.

\par afun 
name of the CppAD defined function that corresponding to this operation.
Note that \c afun, preceeded by a pound sign, 
is a version of \c afun with quotes arround it.

\par forward
name of the user defined function that computes corresponding
results during forward mode.

\par reverse
name of the user defined function that computes corresponding
results during reverse mode.

\par for_jac_sparse
name of the user defined routine that computes corresponding
results during forward mode jacobian sparsity sweeps.

\par rev_jac_sparse
name of the user defined routine that computes corresponding
results during reverse mode jacobian sparsity sweeps.

\par rev_hes_sparse
name of the user defined routine that computes corresponding
results during reverse mode Hessian sparsity sweeps.

\par memory allocation
Note that user_atomic is used as a static object, so its objects
do note get deallocated until the program terminates. 
*/

# define CPPAD_USER_ATOMIC(                                           \
     afun            ,                                                \
     Tvector         ,                                                \
     Base            ,                                                \
	forward         ,                                                \
     reverse         ,                                                \
     for_jac_sparse  ,                                                \
     rev_jac_sparse  ,                                                \
     rev_hes_sparse                                                   \
)                                                                     \
inline void afun (                                                    \
     size_t                               id ,                        \
     const Tvector< CppAD::AD<Base> >&    ax ,                        \
     Tvector< CppAD::AD<Base> >&          ay                          \
)                                                                     \
{	CPPAD_ASSERT_FIRST_CALL_NOT_PARALLEL;                            \
	static CppAD::user_atomic<Base> fun(                             \
          #afun          ,                                            \
          forward        ,                                            \
          reverse        ,                                            \
          for_jac_sparse ,                                            \
          rev_jac_sparse ,                                            \
          rev_hes_sparse                                              \
     );                                                               \
     fun.ad(id, ax, ay);                                              \
}

/*!
Class that actually implements the <tt>afun(id, ax, ay)</tt> calls.

A new user_atomic object is generated each time the user invokes
the CPPAD_USER_ATOMIC macro; see static object in that macro.
*/
template <class Base>
class user_atomic {
	/// type for user routine that computes forward mode results
	typedef bool (*F) (
		size_t                  id ,
		size_t                   k , 
		size_t                   n ,
		size_t                   m ,
		const vector<bool>&     vx ,
		vector<bool>&           vy ,
		const vector<Base>&     tx , 
		vector<Base>&           ty
	);
	/// type for user routine that computes reverse mode results
	typedef bool (*R) (
		size_t                  id ,
		size_t                   k , 
		size_t                   n , 
		size_t                   m , 
		const vector<Base>&     tx , 
		const vector<Base>&     ty ,
		vector<Base>&           px ,
		const vector<Base>&     py
	);
	/// type for user routine that computes forward mode Jacobian sparsity
	typedef bool (*FJS) (
		size_t                           id ,
		size_t                            n ,
		size_t                            m ,
		size_t                            q ,
		const vector< std::set<size_t> >& r ,
		vector< std::set<size_t>  >&      s
	);
	/// type for user routine that computes reverse mode Jacobian sparsity
	typedef bool (*RJS) (
		size_t                           id ,
		size_t                            n ,
		size_t                            m ,
		size_t                            q ,
		vector< std::set<size_t> >&       r ,
		const vector< std::set<size_t> >& s
	);
	/// type for user routine that computes reverse mode Hessian sparsity
	typedef bool (*RHS) (
		size_t                           id ,
		size_t                            n ,
		size_t                            m ,
		size_t                            q ,
		const vector< std::set<size_t> >& r ,
		const vector<bool>&               s ,
		vector<bool>&                     t ,
		const vector< std::set<size_t> >& u ,
		vector< std::set<size_t> >&       v 
	);
private:
	/// users name for the AD version of this atomic operation
	const std::string     name_;
	/// user's implementation of forward mode
	const F                  f_;
	/// user's implementation of reverse mode
	const R                  r_;
	/// user's implementation of forward jacobian sparsity calculations
	const FJS              fjs_;
	/// user's implementation of reverse jacobian sparsity calculations
	const RJS              rjs_;
	/// user's implementation of reverse Hessian sparsity calculations
	const RHS              rhs_;
	/// index of this object in the vector of all objects in this class
	const size_t         index_;

	/// temporary work space used to avoid memory allocation/deallocation
	/// extra information to be passed to the functions
	vector<bool>            vx_;
	vector<bool>            vy_;
	vector<Base>             x_;
	vector<Base>             y_;

	/// List of all objects in this class.
	static std::vector<user_atomic *>& List(void)
	{	CPPAD_ASSERT_FIRST_CALL_NOT_PARALLEL;
		static std::vector<user_atomic *> list;
		return list;
	}
public:
	/*!
 	Constructor called for each invocation of CPPAD_USER_ATOMIC.

	Put this object in the list of all objects for this class and set
	the constant private data name_, f_, r_, and index_.

	\param afun
	is the user's name for the AD version of this atomic operation.

	\param f
	user routine that does forward mode calculations for this operation.

	\param r
	user routine that does reverse mode calculations for this operation.

	\param fjs
	user routine that does forward Jacobian sparsity calculations.

	\param rjs
	user routine that does reverse Jacobian sparsity calculations.

	\param rhs
	user routine that does reverse Hessian sparsity calculations.

	\par
	This constructor can not be used in parallel mode because
	the object it constructs is static. In addition, it changes the
	static object \c List.
	*/
	user_atomic(const char* afun, F f, R r, FJS fjs, RJS rjs, RHS rhs) : 
	name_(afun)
	, f_(f)
	, r_(r)
	, fjs_(fjs)
	, rjs_(rjs)
	, rhs_(rhs)
	, index_( List().size() )
	{	CPPAD_ASSERT_KNOWN(
			! thread_alloc::in_parallel() ,
			"First call to the function *afun is in parallel mode."
		);
		List().push_back(this);
	}

	/*!
 	Implement the user call to <tt>afun(id, ax, ay)</tt>.

	\tparam ADVector
	A simple vector class with elements of type <code>AD<Base></code>.

	\param id
	extra information vector that is just passed through by CppAD,
	and possibly used by user's routines.

	\param ax
	is the argument vector for this call,
	<tt>ax.size()</tt> determines the number of arguments.

	\param ay
	is the result vector for this call,
	<tt>ay.size()</tt> determines the number of results.

	This routine is not \c const because it may modify the works
	space vectors \c x_ and \c y_.
 	*/
	template <class ADVector>
	void ad(size_t id, const ADVector& ax, ADVector& ay)
	{	size_t i, j, k;
		size_t n = ax.size();
		size_t m = ay.size();
# ifndef NDEBUG
		bool ok;
		std::string msg = "user_atomoc: ";
# endif
		//
		if( x_.size() < n )
		{	x_.resize(n);
			vx_.resize(n);
		}
		if( y_.size() < m )
		{	y_.resize(m);
			vy_.resize(m);
		}
		// 
		// Determine if we are going to have to tape this operation
		size_t tape_id     = 0;
		ADTape<Base>* tape = CPPAD_NULL;
		for(j = 0; j < n; j++)
		{	x_[j]   = ax[j].value_;
			vx_[j]  = Variable( ax[j] );
			if ( (tape == CPPAD_NULL) & vx_[j] )
			{	tape    = ax[j].tape_this();
				tape_id = ax[j].id_;
			}
# ifndef NDEBUG
			ok = (tape_id == 0) | Parameter(ax[j]) | (tape_id == ax[j].id_);
			if( ! ok )
			{	msg = msg + name_ + 
				": ax contains variables from different OpenMP threads.";
				CPPAD_ASSERT_KNOWN(false, msg.c_str());
			}
# endif
		}
		// Use zero order forward mode to compute values
		k  = 0;
# if NDEBUG
		f_(id, k, n, m, vx_, vy_, x_, y_);  
# else
		ok = f_(id, k, n, m, vx_, vy_, x_, y_);  
		if( ! ok )
		{	msg = msg + name_ + ": ok returned false from "
				"zero order forward mode calculation.";
			CPPAD_ASSERT_KNOWN(false, msg.c_str());
		}
# endif
		// pass back values
		for(i = 0; i < m; i++)
		{	ay[i].value_ = y_[i];

			// initialize entire vector as a constant (not on tape)
			ay[i].id_    = CPPAD_MAX_NUM_THREADS;
			ay[i].taddr_ = 0;
		}
		// if tape is not null, ay is on the tape
		if( tape != CPPAD_NULL )
		{
			// Note the actual number of results is m
			CPPAD_ASSERT_UNKNOWN( NumRes(UserOp) == 0 );
			CPPAD_ASSERT_UNKNOWN( NumArg(UserOp) == 4 );

			// Begin operators corresponding to one user_atomic operation.
			// Put function index, domain size, range size, and id in tape
			tape->Rec_.PutArg(index_, id, n, m);
			tape->Rec_.PutOp(UserOp);
			// n + m operators follow for this one atomic operation

			// Now put the information for the argument vector in the tape
			CPPAD_ASSERT_UNKNOWN( NumRes(UsravOp) == 0 );
			CPPAD_ASSERT_UNKNOWN( NumRes(UsrapOp) == 0 );
			CPPAD_ASSERT_UNKNOWN( NumArg(UsravOp) == 1 );
			CPPAD_ASSERT_UNKNOWN( NumArg(UsrapOp) == 1 );
			for(j = 0; j < n; j++)
			{	if( vx_[j] )
				{	// information for an argument that is a variable
					tape->Rec_.PutArg(ax[j].taddr_);
					tape->Rec_.PutOp(UsravOp);
				}
				else
				{	// information for an arugment that is parameter
					addr_t p = tape->Rec_.PutPar(ax[j].value_);
					tape->Rec_.PutArg(p);
					tape->Rec_.PutOp(UsrapOp);
				}
			}

			// Now put the information for the results in the tape
			CPPAD_ASSERT_UNKNOWN( NumArg(UsrrpOp) == 1 );
			CPPAD_ASSERT_UNKNOWN( NumRes(UsrrpOp) == 0 );
			CPPAD_ASSERT_UNKNOWN( NumArg(UsrrvOp) == 0 );
			CPPAD_ASSERT_UNKNOWN( NumRes(UsrrvOp) == 1 );
			for(i = 0; i < m; i++)
			{	if( vy_[i] )
				{	ay[i].taddr_ = tape->Rec_.PutOp(UsrrvOp);
					ay[i].id_    = tape_id;
				}
				else
				{	addr_t p = tape->Rec_.PutPar(ay[i].value_);
					tape->Rec_.PutArg(p);
					tape->Rec_.PutOp(UsrrpOp);
				}
			}

			// Put a duplicate UserOp at end of UserOp sequence
			tape->Rec_.PutArg(index_, id, n, m);
			tape->Rec_.PutOp(UserOp);
		} 
		return;
	}

	/// Name corresponding to a user_atomic object
	static const char* name(size_t index)
	{	return List()[index]->name_.c_str(); }
	/*!
 	Link from forward mode sweep to users routine.

	\param index
	index for this function in the list of all user_atomic objects

	\param id
	extra information vector that is just passed through by CppAD,
	and possibly used by user's routines.

	\param k
	order for this forward mode calculation.

	\param n
	domain space size for this calcualtion.

	\param m
	range space size for this calculation.

	\param tx
	Taylor coefficients corresponding to \c x for this calculation.

	\param ty
	Taylor coefficient corresponding to \c y for this calculation

	See the forward mode in user's documentation for user_atomic 
 	*/
	static void forward(
		size_t                index , 
		size_t                   id ,
		size_t                    k ,
		size_t                    n , 
		size_t                    m , 
		const vector<Base>&      tx ,
		vector<Base>&            ty )
	{	
# ifdef _OPENMP
		vector<bool> empty(0);
# else
		static vector<bool> empty(0);
# endif
		
		CPPAD_ASSERT_UNKNOWN( tx.size() >= n * k );
		CPPAD_ASSERT_UNKNOWN( ty.size() >= m * k );

		CPPAD_ASSERT_UNKNOWN(index < List().size() );
		user_atomic* op = List()[index];

		bool ok = op->f_(id, k, n, m, empty, empty, tx, ty);
		if( ! ok )
		{	std::stringstream ss;
			ss << k;
			std::string msg = "user_atomic: ";
			msg = msg + op->name_ + ": ok returned false from " + ss.str()
			    + " order forward mode calculation";
			CPPAD_ASSERT_KNOWN(false, msg.c_str());
		}
	}


	/*!
 	Link from reverse mode sweep to users routine.

	\param index
	index in the list of all user_atomic objects
	corresponding to this function.


	\param id
	extra information vector that is just passed through by CppAD,
	and possibly used by user's routines.

	\param k
	order for this forward mode calculation.

	\param n
	domain space size for this calcualtion.

	\param m
	range space size for this calculation.

	\param tx
	Taylor coefficients corresponding to \c x for this calculation.

	\param ty
	Taylor coefficient corresponding to \c y for this calculation

	\param px
	Partials w.r.t. the \c x Taylor coefficients.

	\param py
	Partials w.r.t. the \c y Taylor coefficients.

	See reverse mode documentation for user_atomic 
 	*/
	static void reverse(
		size_t               index , 
		size_t                  id ,
		size_t                   k ,
		size_t                   n , 
		size_t                   m , 
		const vector<Base>&     tx ,
		const vector<Base>&     ty ,
		vector<Base>&           px ,
		const vector<Base>&     py )
	{
		CPPAD_ASSERT_UNKNOWN(index < List().size() );
		CPPAD_ASSERT_UNKNOWN( tx.size() >= n * k );
		CPPAD_ASSERT_UNKNOWN( px.size() >= n * k );
		CPPAD_ASSERT_UNKNOWN( ty.size() >= m * k );
		CPPAD_ASSERT_UNKNOWN( py.size() >= m * k );
		user_atomic* op = List()[index];

		bool ok = op->r_(id, k, n, m, tx, ty, px, py);
		if( ! ok )
		{	std::stringstream ss;
			ss << k;
			std::string msg = "user_atomic: ";
			msg = op->name_ + ": ok returned false from " + ss.str() 
			    + " order reverse mode calculation";
			CPPAD_ASSERT_KNOWN(false, msg.c_str());
		}
	}

	/*!
 	Link from forward Jacobian sparsity sweep to users routine.

	\param index
	index in the list of all user_atomic objects
	corresponding to this function.

	\param id
	extra information vector that is just passed through by CppAD,
	and possibly used by user's routines.

	\param n
	domain space size for this calcualtion.

	\param m
	range space size for this calculation.

	\param q
	is the column dimension for the Jacobian sparsity partterns.

	\param r
	is the Jacobian sparsity pattern for the argument vector x

	\param s
	is the Jacobian sparsity pattern for the result vector y
	*/
	static void for_jac_sparse(
		size_t                            index ,
		size_t                               id ,
		size_t                                n , 
		size_t                                m , 
		size_t                                q ,
		const vector< std::set<size_t> >&     r ,
		vector< std::set<size_t> >&           s )
	{
		CPPAD_ASSERT_UNKNOWN(index < List().size() );
		CPPAD_ASSERT_UNKNOWN( r.size() >= n );
		CPPAD_ASSERT_UNKNOWN( s.size() >= m );
		user_atomic* op = List()[index];

		bool ok = op->fjs_(id, n, m, q, r, s);
		if( ! ok )
		{	std::string msg = "user_atomic: ";
			msg = msg + op->name_ 
			    + ": ok returned false from for_jac_sparse calculation";
			CPPAD_ASSERT_KNOWN(false, msg.c_str());
		}
	}

	/*!
 	Link from reverse Jacobian sparsity sweep to users routine.

	\param index
	index in the list of all user_atomic objects
	corresponding to this function.

	\param id
	extra information vector that is just passed through by CppAD,
	and possibly used by user's routines.

	\param n
	domain space size for this calcualtion.

	\param m
	range space size for this calculation.

	\param q
	is the row dimension for the Jacobian sparsity partterns.

	\param r
	is the Jacobian sparsity pattern for the argument vector x

	\param s
	is the Jacobian sparsity pattern for the result vector y
	*/
	static void rev_jac_sparse(
		size_t                            index ,
		size_t                               id ,
		size_t                                n , 
		size_t                                m , 
		size_t                                q ,
		vector< std::set<size_t> >&           r ,
		const vector< std::set<size_t> >&     s )
	{
		CPPAD_ASSERT_UNKNOWN(index < List().size() );
		CPPAD_ASSERT_UNKNOWN( r.size() >= n );
		CPPAD_ASSERT_UNKNOWN( s.size() >= m );
		user_atomic* op = List()[index];

		bool ok = op->rjs_(id, n, m, q, r, s);
		if( ! ok )
		{	std::string msg = "user_atomic: ";
			msg = msg + op->name_ 
			    + ": ok returned false from rev_jac_sparse calculation";
			CPPAD_ASSERT_KNOWN(false, msg.c_str());
		}
	}

	/*!
 	Link from reverse Hessian sparsity sweep to users routine.

	\param index
	index in the list of all user_atomic objects
	corresponding to this function.

	\param id
	extra information vector that is just passed through by CppAD,
	and possibly used by user's routines.

	\param n
	domain space size for this calcualtion.

	\param m
	range space size for this calculation.

	\param q
	is the column dimension for the sparsity partterns.

	\param r
	is the forward Jacobian sparsity pattern w.r.t the argument vector x

	\param s
	is the reverse Jacobian sparsity pattern w.r.t the result vector y.

	\param t
	is the reverse Jacobian sparsity pattern w.r.t the argument vector x.

	\param u
	is the Hessian sparsity pattern w.r.t the result vector y.

	\param v
	is the Hessian sparsity pattern w.r.t the argument vector x.
	*/
	static void rev_hes_sparse(
		size_t                            index ,
		size_t                               id ,
		size_t                                n , 
		size_t                                m , 
		size_t                                q ,
		vector< std::set<size_t> >&           r ,
		const vector<bool>&                   s ,
		vector<bool>&                         t ,
		const vector< std::set<size_t> >&     u ,
		vector< std::set<size_t> >&           v )
	{
		CPPAD_ASSERT_UNKNOWN(index < List().size() );
		CPPAD_ASSERT_UNKNOWN( r.size() >= n );
		CPPAD_ASSERT_UNKNOWN( s.size() >= m );
		CPPAD_ASSERT_UNKNOWN( t.size() >= n );
		CPPAD_ASSERT_UNKNOWN( u.size() >= m );
		CPPAD_ASSERT_UNKNOWN( v.size() >= n );
		user_atomic* op = List()[index];

		bool ok = op->rhs_(id, n, m, q, r, s, t, u, v);
		if( ! ok )
		{	std::string msg = "user_atomic: ";
			msg = msg + op->name_ 
			    + ": ok returned false from rev_jac_sparse calculation";
			CPPAD_ASSERT_KNOWN(false, msg.c_str());
		}
	}

	/// Free static CppAD::vector memory used by this class (work space)
	static void clear(void)
	{	size_t i = List().size();
		while(i--)
		{	user_atomic* op = List()[i];
			op->vx_.resize(0);
			op->vy_.resize(0);
			op->x_.resize(0);
			op->y_.resize(0);
		}
		return;
	}
};

CPPAD_END_NAMESPACE
# endif

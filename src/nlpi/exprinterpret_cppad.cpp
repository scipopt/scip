/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    exprinterpret_cppad.cpp
 * @brief   methods to interpret (evaluate) an expression tree "fast" using CppAD
 * @ingroup EXPRINTS
 * @author  Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "nlpi/exprinterpret.h"
#include "scip/def.h"
#include "scip/intervalarith.h"
#include "scip/pub_expr.h"
#include "scip/scip_expr.h"
#include "scip/expr_pow.h"
#include "scip/expr_exp.h"
#include "scip/expr_log.h"
#include "nlpi/expr_varidx.h"

#include <cmath>
#include <cstring>
#include <algorithm>
#include <vector>
using std::vector;

/* Turn off lint warning "747: Significant prototype coercion" and "732: Loss of sign".
 * The first warning is generated for expressions like t[0], where t is a vector, since 0 is an integer constant, but a
 * size_t is expected (usually long unsigned). The second is generated for expressions like t[n], where n is an
 * integer. Both code pieces are likely to be correct. It seems to be impossible to inhibit these messages for
 * vector<*>::operator[] only. */
/*lint --e{747,732}*/

/* defining NO_CPPAD_USER_ATOMIC disables the use of our own implementation of derivatives of power operators
 * via CppAD's user-atomic function feature
 * our customized implementation should give better results (tighter intervals) for the interval data type
 */
/* #define NO_CPPAD_USER_ATOMIC */

/* fallback to non-thread-safe version if C++ is too old to have std::atomic */
#if __cplusplus < 201103L && defined(SCIP_THREADSAFE)
#undef SCIP_THREADSAFE
#endif

/** sign of a value (-1 or +1)
 * 
 * 0.0 has sign +1
 */
#define SIGN(x) ((x) >= 0.0 ? 1.0 : -1.0)

/* CppAD needs to know a fixed upper bound on the number of threads at compile time.
 * It is wise to set it to a power of 2, so that if the tape id overflows, it is likely to start at 0 again, which avoids difficult to debug errors.
 */
#ifndef CPPAD_MAX_NUM_THREADS
#ifdef SCIP_THREADSAFE
#define CPPAD_MAX_NUM_THREADS 64
#else
#define CPPAD_MAX_NUM_THREADS 1
#endif
#endif

/* disable -Wshadow warnings for upcoming includes of CppAD if using some old GCC
 * -Wshadow was too strict with some versions of GCC 4 (https://stackoverflow.com/questions/2958457/gcc-wshadow-is-too-strict)
 */
#ifdef __GNUC__
#if __GNUC__ == 4
#pragma GCC diagnostic ignored "-Wshadow"
#endif
#endif

#include <cppad/cppad.hpp>
#include <cppad/utility/error_handler.hpp>

/* CppAD is not thread-safe by itself, but uses some static datastructures
 * To run it in a multithreading environment, a special CppAD memory allocator that is aware of the multiple threads has to be used.
 * This allocator requires to know the number of threads and a thread number for each thread.
 * To implement this, we follow the team_pthread example of CppAD, which uses pthread's thread-specific data management.
 */
#ifdef SCIP_THREADSAFE

#include <atomic>

/** currently registered number of threads */
static std::atomic_size_t ncurthreads{0};
static thread_local int thread_number{-1};

/** CppAD callback function that indicates whether we are running in parallel mode */
static
bool in_parallel(void)
{
   return ncurthreads > 1;
}

/** CppAD callback function that returns the number of the current thread
 *
 * assigns a new number to the thread if new
 */
static
size_t thread_num(void)
{
   size_t threadnum;

   /* if no thread_number for this thread yet, then assign a new thread number to the current thread
    */
   if( thread_number == -1 )
   {
      thread_number = static_cast<int>(ncurthreads.fetch_add(1, std::memory_order_relaxed));
   }

   threadnum = static_cast<size_t>(thread_number);

   return threadnum;
}

/** sets up CppAD's datastructures for running in multithreading mode
 *
 *  It must be called once before multithreading is started.
 *  For GCC-compatible compilers, this will happen automatically.
 */
extern "C" SCIP_EXPORT char SCIPexprintCppADInitParallel(void);
#ifdef __GNUC__
__attribute__((constructor))
#endif
char SCIPexprintCppADInitParallel(void)
{
   CppAD::thread_alloc::parallel_setup(CPPAD_MAX_NUM_THREADS, in_parallel, thread_num);
   CppAD::parallel_ad<double>();

   return 0;
}

#if !defined(__GNUC__)
/** a dummy variable that is initialized to the result of init_parallel
 *
 *  The purpose is to make sure that init_parallel() is called before any multithreading is started.
 */
static char init_parallel_return = SCIPexprintCppADInitParallel();
#endif

#endif // SCIP_THREADSAFE

using CppAD::AD;

/** expression specific interpreter data */
struct SCIP_ExprIntData
{
public:
   /** constructor */
   SCIP_ExprIntData()
      : val(0.0),
        need_retape(true),
        need_retape_always(false),
        userevalcapability(SCIP_EXPRINTCAPABILITY_ALL),
        hesrowidxs(NULL),
        hescolidxs(NULL),
        hesvalues(NULL),
        hesnnz(0)
   { }

   /** destructor */
   ~SCIP_ExprIntData()
   { }/*lint --e{1540}*/

   /** gives position of index in varidxs vector */
   int getVarPos(
      int                varidx              /**< variable index to look for */
      ) const
   {
      // varidxs is sorted, so can use binary search functions
      assert(std::binary_search(varidxs.begin(), varidxs.end(), varidx));
      return std::lower_bound(varidxs.begin(), varidxs.end(), varidx) - varidxs.begin();
   }

   vector< int >         varidxs;            /**< variable indices used in expression (unique and sorted) */
   vector< AD<double> >  X;                  /**< vector of dependent variables (same size as varidxs) */
   vector< AD<double> >  Y;                  /**< result vector (size 1) */
   CppAD::ADFun<double>  f;                  /**< the function to evaluate as CppAD object */

   vector<double>        x;                  /**< current values of dependent variables (same size as varidxs) */
   double                val;                /**< current function value */
   bool                  need_retape;        /**< will retaping be required for the next point evaluation? */
   bool                  need_retape_always; /**< will retaping be always required? */
   SCIP_EXPRINTCAPABILITY userevalcapability; /**< (intersection of) capabilities of evaluation rountines of user expressions */

   int*                  hesrowidxs;         /**< row indices of Hessian sparsity */
   int*                  hescolidxs;         /**< column indices of Hessian sparsity */
   SCIP_Real*            hesvalues;          /**< values of Hessian */
   int                   hesnnz;             /**< number of nonzeros in Hessian */
};

#ifndef NO_CPPAD_USER_ATOMIC

/** computes sparsity of jacobian for a univariate function during a forward sweep
 *
 *  For a 1 x q matrix R, we have to return the sparsity pattern of the 1 x q matrix S(x) = f'(x) * R.
 *  Since f'(x) is dense, the sparsity of S will be the sparsity of R.
 */
static
bool univariate_for_sparse_jac(
   size_t                     q,             /**< number of columns in R */
   const CppAD::vector<bool>& r,             /**< sparsity of R, columnwise */
   CppAD::vector<bool>&       s              /**< vector to store sparsity of S, columnwise */
   )
{
   assert(r.size() == q);
   assert(s.size() == q);

   s = r;

   return true;
}

/** Computes sparsity of jacobian during a reverse sweep
 *
 *  For a q x 1 matrix R, we have to return the sparsity pattern of the q x 1 matrix S(x) = R * f'(x).
 *  Since f'(x) is dense, the sparsity of S will be the sparsity of R.
 */
static
bool univariate_rev_sparse_jac(
   size_t                     q,             /**< number of rows in R */
   const CppAD::vector<bool>& r,             /**< sparsity of R, rowwise */
   CppAD::vector<bool>&       s              /**< vector to store sparsity of S, rowwise */
   )
{
   assert(r.size() == q);
   assert(s.size() == q);

   s = r;

   return true;
}

/** computes sparsity of hessian during a reverse sweep
 *
 *  Assume V(x) = (g(f(x)))'' R  with f(x) = x^p for a function g:R->R and a matrix R.
 *  we have to specify the sparsity pattern of V(x) and T(x) = (g(f(x)))'.
 */ /*lint -e715*/
static
bool univariate_rev_sparse_hes(
   const CppAD::vector<bool>& vx,            /**< indicates whether argument is a variable, or empty vector */
   const CppAD::vector<bool>& s,             /**< sparsity pattern of S = g'(y) */
   CppAD::vector<bool>&  t,                  /**< vector to store sparsity pattern of T(x) = (g(f(x)))' */
   size_t                q,                  /**< number of columns in R, U, and V */
   const CppAD::vector<bool>& r,             /**< sparsity pattern of R */
   const CppAD::vector<bool>& u,             /**< sparsity pattern of U(x) = g''(f(x)) f'(x) R */
   CppAD::vector<bool>&  v                   /**< vector to store sparsity pattern of V(x) = (g(f(x)))'' R */
   )
{  /*lint --e{439,715}*/  /* @todo take vx into account */
   assert(r.size() == q);
   assert(s.size() == 1);
   assert(t.size() == 1);
   assert(u.size() == q);
   assert(v.size() == q);

   // T(x) = g'(f(x)) * f'(x) = S * f'(x), and f' is not identically 0
   t[0] = s[0];

   // V(x) = g''(f(x)) f'(x) f'(x) R + g'(f(x)) f''(x) R466
   //      = f'(x) U + S f''(x) R, with f'(x) and f''(x) not identically 0
   v = u;
   if( s[0] )
      for( size_t j = 0; j < q; ++j )
         if( r[j] )
            v[j] = true;

   return true;
}


/** Automatic differentiation of x -> x^p, p>=2 integer, as CppAD user-atomic function.
 *
 *  This class implements forward and reverse operations for the function x -> x^p for use within CppAD.
 *  While CppAD would implement integer powers as a recursion of multiplications, we still use pow functions as they allow us to avoid overestimation in interval arithmetics.
 *
 *  @todo treat the exponent as a (variable) argument to the function, with the assumption that we never differentiate w.r.t. it (this should make the approach threadsafe again)
 */
template<class Type>
class atomic_posintpower : public CppAD::atomic_base<Type>
{
public:
   atomic_posintpower()
   : CppAD::atomic_base<Type>("posintpower"),
     exponent(0)
   {
      /* indicate that we want to use bool-based sparsity pattern */
      this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);
   }

private:
   /** exponent value for next call to forward or reverse */
   int exponent;

   /** stores exponent value corresponding to next call to forward or reverse
    *
    * how is this supposed to be threadsafe? (we use only one global instantiation of this class)
    * TODO according to the CppAD 2018 docu, using this function is deprecated; what is the modern way to do this?
    */
   virtual void set_old(size_t id)
   {
      exponent = (int) id;
   }

   /** forward sweep of positive integer power
    *
    * Given the taylor coefficients for x, we have to compute the taylor coefficients for f(x),
    * that is, given tx = (x, x', x'', ...), we compute the coefficients ty = (y, y', y'', ...)
    * in the taylor expansion of f(x) = x^p.
    * Thus, y   = x^p
    *           = tx[0]^p,
    *       y'  = p * x^(p-1) * x'
    *           = p * tx[0]^(p-1) * tx[1],
    *       y'' = 1/2 * p * (p-1) * x^(p-2) * x'^2 + p * x^(p-1) * x''
    *           = 1/2 * p * (p-1) * tx[0]^(p-2) * tx[1]^2 + p * tx[0]^(p-1) * tx[2]
    */
   bool forward(
      size_t                     q,          /**< lowest order Taylor coefficient that we are evaluating */
      size_t                     p,          /**< highest order Taylor coefficient that we are evaluating */
      const CppAD::vector<bool>& vx,         /**< indicates whether argument is a variable, or empty vector */
      CppAD::vector<bool>&       vy,         /**< vector to store which function values depend on variables, or empty vector */
      const CppAD::vector<Type>& tx,         /**< values for taylor coefficients of x */
      CppAD::vector<Type>&       ty          /**< vector to store taylor coefficients of y */
      )
   {
      assert(exponent > 1);
      assert(tx.size() >= p+1);
      assert(ty.size() >= p+1);
      assert(q <= p);

      if( vx.size() > 0 )
      {
         assert(vx.size() == 1);
         assert(vy.size() == 1);
         assert(p == 0);

         vy[0] = vx[0];
      }

      if( q == 0 /* q <= 0 && 0 <= p */ )
      {
         ty[0] = CppAD::pow(tx[0], exponent);
      }

      if( q <= 1 && 1 <= p )
      {
         ty[1] = CppAD::pow(tx[0], exponent-1) * tx[1];
         ty[1] *= double(exponent);
      }

      if( q <= 2 && 2 <= p )
      {
         if( exponent > 2 )
         {
            // ty[2] = 1/2 * exponent * (exponent-1) * pow(tx[0], exponent-2) * tx[1] * tx[1] + exponent * pow(tx[0], exponent-1) * tx[2];
            ty[2]  = CppAD::pow(tx[0], exponent-2) * tx[1] * tx[1];
            ty[2] *= (exponent-1) / 2.0;
            ty[2] += CppAD::pow(tx[0], exponent-1) * tx[2];
            ty[2] *= exponent;
         }
         else
         {
            assert(exponent == 2);
            // ty[2] = 1/2 * exponent * tx[1] * tx[1] + exponent * tx[0] * tx[2];
            ty[2]  = tx[1] * tx[1] + 2.0 * tx[0] * tx[2];
         }
      }

      /* higher order derivatives not implemented */
      if( p > 2 )
         return false;

      return true;
   }

   /** reverse sweep of positive integer power
    *
    * Assume y(x) is a function of the taylor coefficients of f(x) = x^p for x, i.e.,
    *   y(x) = [ x^p, p * x^(p-1) * x', p * (p-1) * x^(p-2) * x'^2 + p * x^(p-1) * x'', ... ].
    * Then in the reverse sweep we have to compute the elements of \f$\partial h / \partial x^[l], l = 0, ..., k,\f$
    * where x^[l] is the l'th taylor coefficient (x, x', x'', ...) and h(x) = g(y(x)) for some function g:R^k -> R.
    * That is, we have to compute
    *\f$
    * px[l] = \partial h / \partial x^[l] = (\partial g / \partial y) * (\partial y / \partial x^[l])
    *       = \sum_{i=0}^k (\partial g / \partial y_i) * (\partial y_i / \partial x^[l])
    *       = \sum_{i=0}^k py[i] * (\partial y_i / \partial x^[l])
    * \f$
    *
    * For k = 0, this means
    *\f$
    * px[0] = py[0] * (\partial y_0 / \partial x^[0])
    *       = py[0] * (\partial x^p / \partial x)
    *       = py[0] * p * tx[0]^(p-1)
    *\f$
    *
    * For k = 1, this means
    * \f$
    * px[0] = py[0] * (\partial y_0 / \partial x^[0]) + py[1] * (\partial y_1 / \partial x^[0])
    *       = py[0] * (\partial x^p / \partial x)     + py[1] * (\partial (p * x^(p-1) * x') / \partial x)
    *       = py[0] * p * tx[0]^(p-1)                 + py[1] * p * (p-1) * tx[0]^(p-2) * tx[1]
    * px[1] = py[0] * (\partial y_0 / \partial x^[1]) + py[1] * (\partial y_1 / \partial x^[1])
    *       = py[0] * (\partial x^p / \partial x')    + py[1] * (\partial (p * x^(p-1) x') / \partial x')
    *       = py[0] * 0                               + py[1] * p * tx[0]^(p-1)
    * \f$
    */ /*lint -e715*/
   bool reverse(
      size_t                     p,          /**< highest order Taylor coefficient that we are evaluating */
      const CppAD::vector<Type>& tx,         /**< values for taylor coefficients of x */
      const CppAD::vector<Type>& ty,         /**< values for taylor coefficients of y */
      CppAD::vector<Type>&       px,         /**< vector to store partial derivatives of h(x) = g(y(x)) w.r.t. x */
      const CppAD::vector<Type>& py          /**< values for partial derivatives of g(x) w.r.t. y */
      )
   { /*lint --e{715}*/
      assert(exponent > 1);
      assert(px.size() >= p+1);
      assert(py.size() >= p+1);
      assert(tx.size() >= p+1);

      switch( p )
      {
         case 0:
            // px[0] = py[0] * exponent * pow(tx[0], exponent-1);
            px[0]  = py[0] * CppAD::pow(tx[0], exponent-1);
            px[0] *= exponent;
            break;

         case 1:
            // px[0] = py[0] * exponent * pow(tx[0], exponent-1) + py[1] * exponent * (exponent-1) * pow(tx[0], exponent-2) * tx[1];
            px[0]  = py[1] * tx[1] * CppAD::pow(tx[0], exponent-2);
            px[0] *= exponent-1;
            px[0] += py[0] * CppAD::pow(tx[0], exponent-1);
            px[0] *= exponent;
            // px[1] = py[1] * exponent * pow(tx[0], exponent-1);
            px[1]  = py[1] * CppAD::pow(tx[0], exponent-1);
            px[1] *= exponent;
            break;

         default:
            return false;
      }

      return true;
   }

   using CppAD::atomic_base<Type>::for_sparse_jac;

   /** computes sparsity of jacobian during a forward sweep
    *
    * For a 1 x q matrix R, we have to return the sparsity pattern of the 1 x q matrix S(x) = f'(x) * R.
    * Since f'(x) is dense, the sparsity of S will be the sparsity of R.
    */
   bool for_sparse_jac(
      size_t                     q,          /**< number of columns in R */
      const CppAD::vector<bool>& r,          /**< sparsity of R, columnwise */
      CppAD::vector<bool>&       s           /**< vector to store sparsity of S, columnwise */
      )
   {
      return univariate_for_sparse_jac(q, r, s);
   }

   using CppAD::atomic_base<Type>::rev_sparse_jac;

   /** computes sparsity of jacobian during a reverse sweep
    *
    *  For a q x 1 matrix R, we have to return the sparsity pattern of the q x 1 matrix S(x) = R * f'(x).
    *  Since f'(x) is dense, the sparsity of S will be the sparsity of R.
    */
   bool rev_sparse_jac(
      size_t                     q,          /**< number of rows in R */
      const CppAD::vector<bool>& r,          /**< sparsity of R, rowwise */
      CppAD::vector<bool>&       s           /**< vector to store sparsity of S, rowwise */
      )
   {
      return univariate_rev_sparse_jac(q, r, s);
   }

   using CppAD::atomic_base<Type>::rev_sparse_hes;

   /** computes sparsity of hessian during a reverse sweep
    *
    *  Assume V(x) = (g(f(x)))'' R  with f(x) = x^p for a function g:R->R and a matrix R.
    *  we have to specify the sparsity pattern of V(x) and T(x) = (g(f(x)))'.
    */
   bool rev_sparse_hes(
      const CppAD::vector<bool>&   vx,       /**< indicates whether argument is a variable, or empty vector */
      const CppAD::vector<bool>&   s,        /**< sparsity pattern of S = g'(y) */
      CppAD::vector<bool>&         t,        /**< vector to store sparsity pattern of T(x) = (g(f(x)))' */
      size_t                       q,        /**< number of columns in R, U, and V */
      const CppAD::vector<bool>&   r,        /**< sparsity pattern of R */
      const CppAD::vector<bool>&   u,        /**< sparsity pattern of U(x) = g''(f(x)) f'(x) R */
      CppAD::vector<bool>&         v         /**< vector to store sparsity pattern of V(x) = (g(f(x)))'' R */
      )
   {
      return univariate_rev_sparse_hes(vx, s, t, q, r, u, v);
   }
};

/** power function with natural exponents */
template<class Type>
static
void posintpower(
   const vector<Type>&   in,                 /**< vector which first argument is base */
   vector<Type>&         out,                /**< vector where to store result in first argument */
   size_t                exponent            /**< exponent */
   )
{
   static atomic_posintpower<typename Type::value_type> pip;
   pip(in, out, exponent);
}

#else

/** power function with natural exponents */
template<class Type>
void posintpower(
   const vector<Type>&   in,                 /**< vector which first argument is base */
   vector<Type>&         out,                /**< vector where to store result in first argument */
   size_t                exponent            /**< exponent */
   )
{
   out[0] = pow(in[0], (int)exponent);
}

#endif


#ifndef NO_CPPAD_USER_ATOMIC

/** Automatic differentiation of x -> sign(x)abs(x)^p, p>=1, as CppAD user-atomic function.
 *
 *  This class implements forward and reverse operations for the function x -> sign(x)abs(x)^p for use within CppAD.
 *  While we otherwise would have to use discontinuous sign and abs functions, our own implementation allows to provide
 *  a continuously differentiable function.
 *
 *  @todo treat the exponent as a (variable) argument to the function, with the assumption that we never differentiate w.r.t. it (this should make the approach threadsafe again)
 */
template<class Type>
class atomic_signpower : public CppAD::atomic_base<Type>
{
public:
   atomic_signpower()
   : CppAD::atomic_base<Type>("signpower"),
     exponent(0.0)
   {
      /* indicate that we want to use bool-based sparsity pattern */
      this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);
   }

private:
   /** exponent for use in next call to forward or reverse */
   SCIP_Real exponent;

   /** stores exponent corresponding to next call to forward or reverse
    *
    * How is this supposed to be threadsafe? (we use only one global instantiation of this class)
    * TODO according to the CppAD 2018 docu, using this function is deprecated; what is the modern way to do this?
    */
   virtual void set_old(size_t id)
   {
      exponent = SCIPgetExponentExprPow((SCIP_EXPR*)(void*)id);
   }

   /** forward sweep of signpower
    *
    * Given the taylor coefficients for x, we have to compute the taylor coefficients for f(x),
    * that is, given tx = (x, x', x'', ...), we compute the coefficients ty = (y, y', y'', ...)
    * in the taylor expansion of f(x) = sign(x)abs(x)^p.
    * Thus, y   = sign(x)abs(x)^p
    *           = sign(tx[0])abs(tx[0])^p,
    *       y'  = p * abs(x)^(p-1) * x'
    *           = p * abs(tx[0])^(p-1) * tx[1],
    *       y'' = 1/2 * p * (p-1) * sign(x) * abs(x)^(p-2) * x'^2 + p * abs(x)^(p-1) * x''
    *           = 1/2 * p * (p-1) * sign(tx[0]) * abs(tx[0])^(p-2) * tx[1]^2 + p * abs(tx[0])^(p-1) * tx[2]
    */
   bool forward(
      size_t                      q,         /**< lowest order Taylor coefficient that we are evaluating */
      size_t                      p,         /**< highest order Taylor coefficient that we are evaluating */
      const CppAD::vector<bool>&  vx,        /**< indicates whether argument is a variable, or empty vector */
      CppAD::vector<bool>&        vy,        /**< vector to store which function values depend on variables, or empty vector */
      const CppAD::vector<Type>&  tx,        /**< values for taylor coefficients of x */
      CppAD::vector<Type>&        ty         /**< vector to store taylor coefficients of y */
      )
   {
      assert(exponent > 0.0);
      assert(tx.size() >= p+1);
      assert(ty.size() >= p+1);
      assert(q <= p);

      if( vx.size() > 0 )
      {
         assert(vx.size() == 1);
         assert(vy.size() == 1);
         assert(p == 0);

         vy[0] = vx[0];
      }

      if( q == 0 /* q <= 0 && 0 <= p */ )
      {
         ty[0] = SIGN(tx[0]) * pow(REALABS(tx[0]), exponent);
      }

      if( q <= 1 && 1 <= p )
      {
            ty[1] = pow(REALABS(tx[0]), exponent - 1.0) * tx[1];
            ty[1] *= exponent;
      }

      if( q <= 2 && 2 <= p )
      {
         if( exponent != 2.0 )
         {
            ty[2]  = SIGN(tx[0]) * pow(REALABS(tx[0]), exponent - 2.0) * tx[1] * tx[1];
            ty[2] *= (exponent - 1.0) / 2.0;
            ty[2] += pow(REALABS(tx[0]), exponent - 1.0) * tx[2];
            ty[2] *= exponent;
         }
         else
         {
            // y'' = 2 (1/2 * sign(x) * x'^2 + |x|*x'') = sign(tx[0]) * tx[1]^2 + 2 * abs(tx[0]) * tx[2]
            ty[2]  = SIGN(tx[0]) * tx[1] * tx[1];
            ty[2] += 2.0 * REALABS(tx[0]) * tx[2];
         }
      }

      /* higher order derivatives not implemented */
      if( p > 2 )
         return false;

      return true;
   }

   /** reverse sweep of signpower
    *
    * Assume y(x) is a function of the taylor coefficients of f(x) = sign(x)|x|^p for x, i.e.,
    *   y(x) = [ f(x), f'(x), f''(x), ... ].
    * Then in the reverse sweep we have to compute the elements of \f$\partial h / \partial x^[l], l = 0, ..., k,\f$
    * where x^[l] is the l'th taylor coefficient (x, x', x'', ...) and h(x) = g(y(x)) for some function g:R^k -> R.
    * That is, we have to compute
    *\f$
    * px[l] = \partial h / \partial x^[l] = (\partial g / \partial y) * (\partial y / \partial x^[l])
    *       = \sum_{i=0}^k (\partial g / \partial y_i) * (\partial y_i / \partial x^[l])
    *       = \sum_{i=0}^k py[i] * (\partial y_i / \partial x^[l])
    *\f$
    *
    * For k = 0, this means
    *\f$
    * px[0] = py[0] * (\partial y_0 / \partial x^[0])
    *       = py[0] * (\partial f(x) / \partial x)
    *       = py[0] * p * abs(tx[0])^(p-1)
    * \f$
    *
    * For k = 1, this means
    *\f$
    * px[0] = py[0] * (\partial y_0  / \partial x^[0]) + py[1] * (\partial y_1   / \partial x^[0])
    *       = py[0] * (\partial f(x) / \partial x)     + py[1] * (\partial f'(x) / \partial x)
    *       = py[0] * p * abs(tx[0])^(p-1)             + py[1] * p * (p-1) * abs(tx[0])^(p-2) * sign(tx[0]) * tx[1]
    * px[1] = py[0] * (\partial y_0  / \partial x^[1]) + py[1] * (\partial y_1 / \partial x^[1])
    *       = py[0] * (\partial f(x) / \partial x')    + py[1] * (\partial f'(x) / \partial x')
    *       = py[0] * 0                                + py[1] * p * abs(tx[0])^(p-1)
    * \f$
    */
   bool reverse(
      size_t                      p,         /**< highest order Taylor coefficient that we are evaluating */
      const CppAD::vector<Type>&  tx,        /**< values for taylor coefficients of x */
      const CppAD::vector<Type>&  ty,        /**< values for taylor coefficients of y */
      CppAD::vector<Type>&        px,        /**< vector to store partial derivatives of h(x) = g(y(x)) w.r.t. x */
      const CppAD::vector<Type>&  py         /**< values for partial derivatives of g(x) w.r.t. y */
      )
   { /*lint --e{715}*/
      assert(exponent > 1);
      assert(px.size() >= p+1);
      assert(py.size() >= p+1);
      assert(tx.size() >= p+1);

      switch( p )
      {
      case 0:
         // px[0] = py[0] * p * pow(abs(tx[0]), p-1);
         px[0]  = py[0] * pow(REALABS(tx[0]), exponent - 1.0);
         px[0] *= p;
         break;

      case 1:
         if( exponent != 2.0 )
         {
            // px[0] = py[0] * p * abs(tx[0])^(p-1) + py[1] * p * (p-1) * abs(tx[0])^(p-2) * sign(tx[0]) * tx[1]
            px[0]  = py[1] * tx[1] * pow(REALABS(tx[0]), exponent - 2.0) * SIGN(tx[0]);
            px[0] *= exponent - 1.0;
            px[0] += py[0] * pow(REALABS(tx[0]), exponent - 1.0);
            px[0] *= exponent;
            // px[1] = py[1] * p * abs(tx[0])^(p-1)
            px[1]  = py[1] * pow(REALABS(tx[0]), exponent - 1.0);
            px[1] *= exponent;
         }
         else
         {
            // px[0] = py[0] * 2.0 * abs(tx[0]) + py[1] * 2.0 * sign(tx[0]) * tx[1]
            px[0]  = py[1] * tx[1] * SIGN(tx[0]);
            px[0] += py[0] * REALABS(tx[0]);
            px[0] *= 2.0;
            // px[1] = py[1] * 2.0 * abs(tx[0])
            px[1]  = py[1] * REALABS(tx[0]);
            px[1] *= 2.0;
         }
         break;

      default:
         return false;
      }

      return true;
   }

   using CppAD::atomic_base<Type>::for_sparse_jac;

   /** computes sparsity of jacobian during a forward sweep
    *
    * For a 1 x q matrix R, we have to return the sparsity pattern of the 1 x q matrix S(x) = f'(x) * R.
    * Since f'(x) is dense, the sparsity of S will be the sparsity of R.
    */
   bool for_sparse_jac(
      size_t                     q,          /**< number of columns in R */
      const CppAD::vector<bool>& r,          /**< sparsity of R, columnwise */
      CppAD::vector<bool>&       s           /**< vector to store sparsity of S, columnwise */
      )
   {
      return univariate_for_sparse_jac(q, r, s);
   }

   using CppAD::atomic_base<Type>::rev_sparse_jac;

   /** computes sparsity of jacobian during a reverse sweep
    *
    *  For a q x 1 matrix R, we have to return the sparsity pattern of the q x 1 matrix S(x) = R * f'(x).
    *  Since f'(x) is dense, the sparsity of S will be the sparsity of R.
    */
   bool rev_sparse_jac(
      size_t                     q,          /**< number of rows in R */
      const CppAD::vector<bool>& r,          /**< sparsity of R, rowwise */
      CppAD::vector<bool>&       s           /**< vector to store sparsity of S, rowwise */
      )
   {
      return univariate_rev_sparse_jac(q, r, s);
   }

   using CppAD::atomic_base<Type>::rev_sparse_hes;

   /** computes sparsity of hessian during a reverse sweep
    *
    * Assume V(x) = (g(f(x)))'' R  with f(x) = sign(x)abs(x)^p for a function g:R->R and a matrix R.
    * we have to specify the sparsity pattern of V(x) and T(x) = (g(f(x)))'.
    */
   bool rev_sparse_hes(
      const CppAD::vector<bool>& vx,         /**< indicates whether argument is a variable, or empty vector */
      const CppAD::vector<bool>& s,          /**< sparsity pattern of S = g'(y) */
      CppAD::vector<bool>&       t,          /**< vector to store sparsity pattern of T(x) = (g(f(x)))' */
      size_t                     q,          /**< number of columns in S and R */
      const CppAD::vector<bool>& r,          /**< sparsity pattern of R */
      const CppAD::vector<bool>& u,          /**< sparsity pattern of U(x) = g''(f(x)) f'(x) R */
      CppAD::vector<bool>&       v           /**< vector to store sparsity pattern of V(x) = (g(f(x)))'' R */
      )
   {
      return univariate_rev_sparse_hes(vx, s, t, q, r, u, v);
   }

};

/** template for evaluation for signpower operator */
template<class Type>
static
void evalSignPower(
   Type&                 resultant,          /**< resultant */
   const Type&           arg,                /**< operand */
   SCIP_EXPR*            expr                /**< expression that holds the exponent */
   )
{
   vector<Type> in(1, arg);
   vector<Type> out(1);

   static atomic_signpower<typename Type::value_type> sp;
   sp(in, out, (size_t)(void*)expr);

   resultant = out[0];
   return;
}

#else

/** template for evaluation for signpower operator
 *
 *  Only implemented for real numbers, thus gives error by default.
 */
template<class Type>
static
void evalSignPower(
   Type&                 resultant,          /**< resultant */
   const Type&           arg,                /**< operand */
   SCIP_EXPR*            expr                /**< expression that holds the exponent */
   )
{  /*lint --e{715}*/
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "evalSignPower()",
      "Error: SignPower not implemented for this value type"
      );
}

/** specialization of signpower evaluation for real numbers */
template<>
void evalSignPower(
   CppAD::AD<double>&    resultant,          /**< resultant */
   const CppAD::AD<double>& arg,             /**< operand */
   SCIP_EXPR*            expr                /**< expression that holds the exponent */
   )
{
   AD<double> adzero(0.);
   SCIP_Real exponent;

   exponent = SCIPgetExponentExprPow(expr);
   assert(exponent >= 1.0);

   if( EPSISINT(exponent, 0.0) )
   {
      // TODO or use powintpower()?
      resultant = CppAD::CondExpGe(arg, adzero, pow(arg, (int)exponent), -pow(-arg, (int)exponent));
   }
   else
   {
      /* pow(0,fractional>1) is not differential in the CppAD world (https://github.com/coin-or/CppAD/discussions/93)
       * this works around this by evaluating pow(eps,exponent) in this case
       */
      resultant = CppAD::CondExpEq(arg, adzero, pow(arg+std::numeric_limits<SCIP_Real>::epsilon(), exponent)-pow(std::numeric_limits<SCIP_Real>::epsilon(), exponent),
         CppAD::CondExpGe(arg, adzero, pow(arg, exponent), -pow(-arg, exponent)));
   }
}

#endif

#if !1  // FIXME
#ifndef NO_CPPAD_USER_ATOMIC

template<class Type>
SCIP_RETCODE exprEvalUser(
   SCIP_EXPR* expr,
   Type* x,
   Type& funcval,
   Type* gradient,
   Type* hessian
   )
{
   return SCIPexprEvalUser(expr, x, &funcval, gradient, hessian); /*lint !e429*/
}

/** Automatic differentiation of user expression as CppAD user-atomic function.
 *
 * This class implements forward and reverse operations for a function given by a user expression for use within CppAD.
 */
template<class Type>
class atomic_userexpr : public CppAD::atomic_base<Type>
{
public:
   atomic_userexpr()
   : CppAD::atomic_base<Type>("userexpr"),
     expr(NULL)
   {
      /* indicate that we want to use bool-based sparsity pattern */
      this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);
   }

private:
   /** user expression */
   SCIP_EXPR* expr;

   /** stores user expression corresponding to next call to forward or reverse
    *
    * how is this supposed to be threadsafe? (we use only one global instantiation of this class)
    * TODO according to the CppAD 2018 docu, using this function is deprecated; what is the modern way to do this?
    */
   virtual void set_old(size_t id)
   {
      expr = (SCIP_EXPR*)(void*)id;
      assert(SCIPexprGetOperator(expr) == SCIP_EXPR_USER);
   }

   /** forward sweep of userexpr
    *
    * We follow http://www.coin-or.org/CppAD/Doc/atomic_forward.xml
    *   Note, that p and q are interchanged!
    *
    * For a scalar variable t, let
    *   Y(t) = f(X(t))
    *   X(t) = x^0 + x^1 t^1 + ... + x^p t^p
    * where for x^i the i an index, while for t^i the i is an exponent.
    * Thus, x^k = 1/k! X^(k) (0),   where X^(k)(.) denotes the k-th derivative.
    *
    * Next, let y^k = 1/k! Y^(k)(0) be the k'th taylor coefficient of Y. Thus,
    *   y^0 = Y^(0)(0)     =     Y(0)   = f(X(0)) = f(x^0)
    *   y^1 = Y^(1)(0)     =     Y'(0)  = f'(X(0)) * X'(0) = f'(x^0) * x^1
    *   y^2 = 1/2 Y^(2)(0) = 1/2 Y''(0) = 1/2 X'(0) * f''(X(0)) X'(0) + 1/2 * f'(X(0)) * X''(0) = 1/2 x^1 * f''(x^0) * x^1 + f'(x^0) * x^2
    *
    * As x^k = (tx[k], tx[(p+1)+k], tx[2*(p+1)+k], ..., tx[n*(p+1)+k], we get
    *   ty[0] = y^0 = f(x^0) = f(tx[{1..n}*(p+1)])
    *   ty[1] = y^1 = f'(x^0) * tx[{1..n}*(p+1)+1] = sum(i=1..n, grad[i] * tx[i*(p+1)+1]),  where grad = f'(x^0)
    *   ty[2] = 1/2 sum(i,j=1..n, x[i*(p+1)+1] * x[j*(p+1)+q] * hessian[i,j]) + sum(i=1..n, grad[i] * x[i*(p+1)+2])
    */
   bool forward(
      size_t                      q,            /**< lowest order Taylor coefficient that we are evaluating */
      size_t                      p,            /**< highest order Taylor coefficient that we are evaluating */
      const CppAD::vector<bool>&  vx,           /**< indicates whether argument is a variable, or empty vector */
      CppAD::vector<bool>&        vy,           /**< vector to store which function values depend on variables, or empty vector */
      const CppAD::vector<Type>&  tx,           /**< values for taylor coefficients of x */
      CppAD::vector<Type>&        ty            /**< vector to store taylor coefficients of y */
   )
   {
      assert(expr != NULL);
      assert(ty.size() == p+1);
      assert(q <= p);

      size_t n = tx.size() / (p+1);
      assert(n == (size_t)SCIPexprGetNChildren(expr)); /*lint !e571*/
      assert(n >= 1);

      if( vx.size() > 0 )
      {
         assert(vx.size() == n);
         assert(vy.size() == 1);
         assert(p == 0);

         /* y_0 is a variable if at least one of the x_i is a variable */
         vy[0] = false;
         for( size_t i = 0; i < n; ++i )
            if( vx[i] )
            {
               vy[0] = true;
               break;
            }
      }

      Type* x = new Type[n];
      for( size_t i = 0; i < n; ++i )
         x[i] = tx[i * (p+1) + 0];  /*lint !e835*/

      if( SCIPcallExprEval(scip, expr, x, &ty[0]) != SCIP_OKAY )
      {
         delete[] x;
         return false;
      }

      Type* gradient = NULL;
      Type* hessian = NULL;

      if( q <= 2 && 1 <= p )
         gradient = new Type[n];
      if( q <= 2 && 2 <= p )
         hessian = new Type[n*n];

      if( exprEvalUser(expr, x, ty[0], gradient, hessian) != SCIP_OKAY )
      {
         delete[] x;
         delete[] gradient;
         delete[] hessian;
         return false;
      }

      if( gradient != NULL )
      {
         ty[1] = 0.0;
         for( size_t i = 0; i < n; ++i )
            ty[1] += gradient[i] * tx[i * (p+1) + 1];
      }

      if( hessian != NULL )
      {
         assert(gradient != NULL);

         ty[2] = 0.0;
         for( size_t i = 0; i < n; ++i )
         {
            for( size_t j = 0; j < n; ++j )
               ty[2] += 0.5 * hessian[i*n+j] * tx[i * (p+1) + 1] * tx[j * (p+1) + 1];

            ty[2] += gradient[i] * tx[i * (p+1) + 2];
         }
      }

      delete[] x;
      delete[] gradient;
      delete[] hessian;

      /* higher order derivatives not implemented */
      if( p > 2 )
         return false;

      return true;
   }

   /** reverse sweep of userexpr
    *
    * We follow http://www.coin-or.org/CppAD/Doc/atomic_reverse.xml
    *   Note, that there q is our p.
    *
    * For a scalar variable t, let
    *   Y(t) = f(X(t))
    *   X(t) = x^0 + x^1 t^1 + ... + x^p t^p
    * where for x^i the i an index, while for t^i the i is an exponent.
    * Thus, x^k = 1/k! X^(k) (0),   where X^(k)(.) denotes the k-th derivative.
    *
    * Next, let y^k = 1/k! Y^(k)(0) be the k'th taylor coefficient of Y. Thus,
    *   Y(t) = y^0 + y^1 t^1 + y^2 t^2 + ...
    * y^0, y^1, ... are the taylor coefficients of f(x).
    *
    * Further, let F(x^0,..,x^p) by given as F^k(x) = y^k. Thus,
    *   F^0(x) = y^0 = Y^(0)(0)   = f(x^0)
    *   F^1(x) = y^1 = Y^(1)(0)   = f'(x^0) * x^1
    *   F^2(x) = y^2 = 1/2 Y''(0) = 1/2 x^1 f''(x^0) x^1 + f'(x^0) x^2
    *
    * Given functions G: R^(p+1) -> R and H: R^(n*(p+1)) -> R, where H(x^0, x^1, .., x^p) = G(F(x^0,..,x^p)),
    * we have to return the value of \f$\partial H / \partial x^l, l = 0..p,\f$ in px. Therefor,
    * \f$
    *  px^l = \partial H / \partial x^l
    *       = sum(k=0..p, (\partial G / \partial y^k) * (\partial y^k / \partial x^l)
    *       = sum(k=0..p, py[k] * (\partial F^k / \partial x^l)
    * \f$
    *
    * For p = 0, this means
    * \f$
    *  px^0 = py[0] * \partial F^0 / \partial x^0
    *       = py[0] * \partial f(x^0) / \partial x^0
    *       = py[0] * f'(x^0)
    * \f$
    *
    * For p = 1, this means (for l = 0):
    * \f[
    *  px^0 = py[0] * \partial F^0    / \partial x^0 + py[1] * \partial F^1 / \partial x^0
    *       = py[0] * \partial f(x^0) / \partial x^0 + py[1] * \partial (f'(x^0) * x^1) / \partial x^0
    *       = py[0] * f'(x^0)                         + py[1] * f''(x^0) * x^1
    * \f]
    * and (for l=1):
    * \[
    *  px^1 = py[0] * \partial F^0    / \partial x^1 + py[1] * \partial F^1 / \partial x^1
    *       = py[0] * \partial f(x^0) / \partial x^1 + py[1] * \partial (f'(x^0) * x^1) / \partial x^0
    *       = py[0] * 0                               + py[1] * f'(x^0)
    * \f]
    *
    * As x^k = (tx[k], tx[(p+1)+k], tx[2*(p+1)+k], ..., tx[n*(p+1)+k] and
    *   px^k = (px[k], px[(p+1)+k], px[2*(p+1)+k], ..., px[n*(p+1)+k], we get
    * for p = 0:
    *   px[i] = (px^0)_i = py[0] * grad[i]
    * for p = 1:
    *   px[i*2+0] = (px^0)_i = py[0] * grad[i] + py[1] * sum(j, hessian[j,i] * tx[j*2+1])
    *   px[i*2+1] = (px^1)_i = py[1] * grad[i]
    */
   bool reverse(
      size_t                      p,            /**< highest order Taylor coefficient that we are evaluating */
      const CppAD::vector<Type>&  tx,           /**< values for taylor coefficients of x */
      const CppAD::vector<Type>&  ty,           /**< values for taylor coefficients of y */
      CppAD::vector<Type>&        px,           /**< vector to store partial derivatives of h(x) = g(y(x)) w.r.t. x */
      const CppAD::vector<Type>&  py            /**< values for partial derivatives of g(x) w.r.t. y */
      )
   {
      assert(expr != NULL);
      assert(px.size() == tx.size());
      assert(py.size() == p+1);

      size_t n = tx.size() / (p+1);
      assert(n == (size_t)SCIPexprGetNChildren(expr)); /*lint !e571*/
      assert(n >= 1);

      Type* x = new Type[n];
      Type funcval;
      Type* gradient = new Type[n];
      Type* hessian = NULL;

      if( p == 1 )
         hessian = new Type[n*n];

      for( size_t i = 0; i < n; ++i )
         x[i] = tx[i * (p+1) + 0]; /*lint !e835*/

      if( exprEvalUser(expr, x, funcval, gradient, hessian) != SCIP_OKAY )
      {
         delete[] x;
         delete[] gradient;
         delete[] hessian;
         return false;
      }

      switch( p )
      {
      case 0:
         // px[j] = (px^0)_j = py[0] * grad[j]
         for( size_t i = 0; i < n; ++i )
            px[i] = py[0] * gradient[i];
         break;

      case 1:
         //  px[i*2+0] = (px^0)_i = py[0] * grad[i] + py[1] * sum(j, hessian[j,i] * tx[j*2+1])
         //  px[i*2+1] = (px^1)_i = py[1] * grad[i]
         assert(hessian != NULL);
         for( size_t i = 0; i < n; ++i )
         {
            px[i*2+0] = py[0] * gradient[i]; /*lint !e835*/
            for( size_t j = 0; j < n; ++j )
               px[i*2+0] += py[1] * hessian[i+n*j] * tx[j*2+1]; /*lint !e835*/

            px[i*2+1] = py[1] * gradient[i];
         }
         break;

      default:
         return false;
      }

      return true;
   } /*lint !e715*/

   using CppAD::atomic_base<Type>::for_sparse_jac;

   /** computes sparsity of jacobian during a forward sweep
    * For a 1 x q matrix R, we have to return the sparsity pattern of the 1 x q matrix S(x) = f'(x) * R.
    * Since we assume f'(x) to be dense, the sparsity of S will be the sparsity of R.
    */
   bool for_sparse_jac(
      size_t                     q,  /**< number of columns in R */
      const CppAD::vector<bool>& r,  /**< sparsity of R, columnwise */
      CppAD::vector<bool>&       s   /**< vector to store sparsity of S, columnwise */
      )
   {
      assert(expr != NULL);
      assert(s.size() == q);

      size_t n = r.size() / q;
      assert(n == (size_t)SCIPexprGetNChildren(expr)); /*lint !e571*/

      // sparsity for S(x) = f'(x) * R
      for( size_t j = 0; j < q; j++ )
      {
         s[j] = false;
         for( size_t i = 0; i < n; i++ )
            s[j] |= (bool)r[i * q + j]; /*lint !e1786*/
      }

      return true;
   }

   using CppAD::atomic_base<Type>::rev_sparse_jac;

   /** computes sparsity of jacobian during a reverse sweep
    * For a q x 1 matrix S, we have to return the sparsity pattern of the q x 1 matrix R(x) = S * f'(x).
    * Since we assume f'(x) to be dense, the sparsity of R will be the sparsity of S.
    */
   bool rev_sparse_jac(
      size_t                     q,  /**< number of rows in R */
      const CppAD::vector<bool>&       rt, /**< sparsity of R, rowwise */
      CppAD::vector<bool>& st  /**< vector to store sparsity of S, rowwise */
      )
   {
      assert(expr != NULL);
      assert(rt.size() == q);

      size_t n = st.size() / q;
      assert(n == (size_t)SCIPexprGetNChildren(expr));

      // sparsity for S(x)^T = f'(x)^T * R^T
      for( size_t j = 0; j < q; j++ )
         for( size_t i = 0; i < n; i++ )
            st[i * q + j] = rt[j];

      return true;
   }

   using CppAD::atomic_base<Type>::rev_sparse_hes;

   /** computes sparsity of hessian during a reverse sweep
    * Assume V(x) = (g(f(x)))'' R  for a function g:R->R and a matrix R.
    * we have to specify the sparsity pattern of V(x) and T(x) = (g(f(x)))'.
    */
   bool rev_sparse_hes(
      const CppAD::vector<bool>&              vx, /**< indicates whether argument is a variable, or empty vector */
      const CppAD::vector<bool>&              s,  /**< sparsity pattern of S = g'(y) */
      CppAD::vector<bool>&                    t,  /**< vector to store sparsity pattern of T(x) = (g(f(x)))' */
      size_t                                  q,  /**< number of columns in S and R */
      const CppAD::vector<bool>& r,  /**< sparsity pattern of R */
      const CppAD::vector<bool>& u,  /**< sparsity pattern of U(x) = g''(f(x)) f'(x) R */
      CppAD::vector<bool>&      v   /**< vector to store sparsity pattern of V(x) = (g(f(x)))'' R */
   )
   {
      assert(expr != NULL);
      size_t n = vx.size();
      assert((size_t)SCIPexprGetNChildren(expr) == n);
      assert(s.size() == 1);
      assert(t.size() == n);
      assert(r.size() == n * q);
      assert(u.size() == q);
      assert(v.size() == n * q);

      size_t i, j, k;

      // sparsity for T(x) = S(x) * f'(x)
      for( i = 0; i < n; ++i )
         t[i] = s[0];

      // V(x) = f'(x)^T * g''(y) * f'(x) * R  +  g'(y) * f''(x) * R
      // U(x) = g''(y) * f'(x) * R
      // S(x) = g'(y)

      // back propagate the sparsity for U
      for( j = 0; j < q; j++ )
         for( i = 0; i < n; i++ )
            v[ i * q + j] = u[j];

      // include forward Jacobian sparsity in Hessian sparsity
      // sparsity for g'(y) * f''(x) * R  (Note f''(x) is assumed to be dense)
      if( s[0] )
         for( j = 0; j < q; j++ )
            for( i = 0; i < n; i++ )
               for( k = 0; k < n; ++k )
                  v[ i * q + j] |= (bool) r[ k * q + j];

      return true;
   }

};

template<class Type>
static
void evalUser(
   Type&                 resultant,          /**< resultant */
   const Type*           args,               /**< operands */
   SCIP_EXPR*            expr                /**< expression that holds the user expression */
   )
{
   assert( args != 0 );
   vector<Type> in(args, args + SCIPexprGetNChildren(expr));
   vector<Type> out(1);

   static atomic_userexpr<typename Type::value_type> u;
   u(in, out, (size_t)(void*)expr);

   resultant = out[0];
   return;
}

#else

template<class Type>
static
void evalUser(
   Type&                 resultant,          /**< resultant */
   const Type*           args,               /**< operands */
   SCIP_EXPR*            expr                /**< expression that holds the user expression */
   )
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "evalUser()",
      "Error: user expressions in CppAD not possible without CppAD user atomic facility"
      );
}

#endif
#endif

/** integer power operation for arbitrary integer exponents */
template<class Type>
static
void evalIntPower(
   Type&                 resultant,          /**< resultant */
   const Type&           arg,                /**< operand */
   const int             exponent            /**< exponent */
   )
{
   if( exponent > 1 )
   {
      vector<Type> in(1, arg);
      vector<Type> out(1);

      posintpower(in, out, exponent);

      resultant = out[0];
      return;
   }

   if( exponent < -1 )
   {
      vector<Type> in(1, arg);
      vector<Type> out(1);

      posintpower(in, out, -exponent);

      resultant = Type(1.0)/out[0];
      return;
   }

   if( exponent == 1 )
   {
      resultant = arg;
      return;
   }

   if( exponent == 0 )
   {
      resultant = Type(1.0);
      return;
   }

   assert(exponent == -1);
   resultant = Type(1.0)/arg;
}

/** CppAD compatible evaluation of an expression for given arguments and parameters */
template<class Type>
static
SCIP_RETCODE eval(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRINTDATA*     exprintdata,        /**< interpreter data for root expression */
   const vector<Type>&   x,                  /**< values of variables */
   Type&                 val                 /**< buffer to store expression value */
   )
{
   Type* buf = NULL;

   assert(expr != NULL);

   // TODO this should iterate instead of using recursion
   //   but the iterdata wouldn't work to hold Type at the moment
   //   they could hold Type*, but then we need to alloc small portions all the time
   //   or we have a big Type-array outside and point to it in iterdata

   if( SCIPisExprVaridx(scip, expr) )
   {
      val = x[exprintdata->getVarPos(SCIPgetIndexExprVaridx(expr))];
      return SCIP_OKAY;
   }
   if( SCIPisExprValue(scip, expr) )
   {
      val = SCIPgetValueExprValue(expr);
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &buf, SCIPexprGetNChildren(expr)) );
   for( int i = 0; i < SCIPexprGetNChildren(expr); ++i )
   {
      SCIP_CALL( eval(scip, SCIPexprGetChildren(expr)[i], exprintdata, x, buf[i]) );
   }

   if( SCIPisExprSum(scip, expr) )
   {
      val = SCIPgetConstantExprSum(expr);
      for( int i = 0; i < SCIPexprGetNChildren(expr); ++i )
         val += SCIPgetCoefsExprSum(expr)[i] * buf[i];
   }
   else if( SCIPisExprProduct(scip, expr) )
   {
      val = SCIPgetCoefExprProduct(expr);
      for( int i = 0; i < SCIPexprGetNChildren(expr); ++i )
         val *= buf[i];
   }
   else if( SCIPisExprPower(scip, expr) )
   {
      SCIP_Real exponent = SCIPgetExponentExprPow(expr);
      if( EPSISINT(exponent, 0.0) )
         evalIntPower(val, buf[0], (int)SCIPgetExponentExprPow(expr));
      else if( exponent == 0.5 )
         val = sqrt(buf[0]);
      else if( exponent < 1.0 )
         val = CppAD::pow(buf[0], SCIPgetExponentExprPow(expr));
      else
      {
         // workaround bug where CppAD claims pow(x,fractional>0) is nondiff at x=0
         // https://github.com/coin-or/CppAD/discussions/93#discussioncomment-327876
         AD<double> adzero(0.);
         val = CppAD::CondExpEq(buf[0], adzero, pow(buf[0]+std::numeric_limits<SCIP_Real>::epsilon(), exponent)-pow(std::numeric_limits<SCIP_Real>::epsilon(), exponent),
            pow(buf[0], exponent));
      }
   }
   else if( SCIPisExprSignpower(scip, expr) )
   {
      evalSignPower(val, buf[0], expr);
   }
   else if( SCIPisExprExp(scip, expr) )
   {
      val = exp(buf[0]);
   }
   else if( SCIPisExprLog(scip, expr) )
   {
      val = log(buf[0]);
   }
   else if( strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), "sin") == 0 )
   {
      val = sin(buf[0]);
   }
   else if( strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), "cos") == 0 )
   {
      val = cos(buf[0]);
   }
   else if( strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), "erf") == 0 )
   {
      val = erf(buf[0]);
   }
   else if( strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), "abs") == 0 )
   {
      val = abs(buf[0]);
   }
   else if( strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), "entropy") == 0 )
   {
      /* -x*log(x) if x > 0, else 0
       * https://coin-or.github.io/CppAD/doc/cond_exp.cpp.htm suggest to use 0 for the x=0 case
       * but then derivatives at 0 are 0, while they are actually infinite (see expr_entropy.c:bwdiff)
       * so we use -sqrt(x) for the x=0 case, as this also has value 0 and first derivative -inf at 0
       */
      val = CppAD::CondExpGt(buf[0], AD<double>(0.), -buf[0] * log(buf[0]), -sqrt(buf[0]));
   }
   else
   {
      // FIXME evalUser(val, buf, expr);
      SCIPerrorMessage("exprint using expr derivative callbacks not yet implemented, (hdlr %s)\n", SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)));
      SCIPfreeBufferArrayNull(scip, &buf);
      return SCIP_ERROR;
   }

   SCIPfreeBufferArray(scip, &buf);

   return SCIP_OKAY;
}

/** replacement for CppAD's default error handler
 *
 *  In debug mode, CppAD gives an error when an evaluation contains a nan.
 *  We do not want to stop execution in such a case, since the calling routine should check for nan's and decide what to do.
 *  Since we cannot ignore this particular error, we ignore all.
 *  @todo find a way to check whether the error corresponds to a nan and communicate this back
 */
static
void cppaderrorcallback(
   bool                  known,              /**< is the error from a known source? */
   int                   line,               /**< line where error occured */
   const char*           file,               /**< file where error occured */
   const char*           cond,               /**< error condition */
   const char*           msg                 /**< error message */
   )
{
   SCIPdebugMessage("ignore CppAD error from %sknown source %s:%d: msg: %s exp: %s\n", known ? "" : "un", file, line, msg, cond);
}

/* install our error handler */
static CppAD::ErrorHandler errorhandler(cppaderrorcallback);

/** gets name and version of expression interpreter */
const char* SCIPexprintGetName(void)
{
   return CPPAD_PACKAGE_STRING;
}

/** gets descriptive text of expression interpreter */
const char* SCIPexprintGetDesc(void)
{
   return "Algorithmic Differentiation of C++ algorithms developed by B. Bell (www.coin-or.org/CppAD)";
}

/** gets capabilities of expression interpreter (using bitflags) */
SCIP_EXPRINTCAPABILITY SCIPexprintGetCapability(
   void
   )
{
   return SCIP_EXPRINTCAPABILITY_FUNCVALUE | SCIP_EXPRINTCAPABILITY_GRADIENT | SCIP_EXPRINTCAPABILITY_HESSIAN;
}

/** creates an expression interpreter object */
SCIP_RETCODE SCIPexprintCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT**        exprint             /**< buffer to store pointer to expression interpreter */
   )
{
   assert(exprint != NULL);

   *exprint = (SCIP_EXPRINT*)1u;  /* some code checks that a non-NULL pointer is returned here, even though it may not point anywhere */

   return SCIP_OKAY;
}

/** frees an expression interpreter object */
SCIP_RETCODE SCIPexprintFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT**        exprint             /**< expression interpreter that should be freed */
   )
{
   assert(exprint != NULL);

   *exprint = NULL;

   return SCIP_OKAY;
}

/** compiles an expression and returns interpreter-specific data for expression
 *
 * @attention the expression is assumed to use varidx expressions but no var expressions
 */
SCIP_RETCODE SCIPexprintCompile(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            rootexpr,           /**< expression */
   SCIP_EXPRINTDATA**    exprintdata         /**< buffer to store pointer to compiled data */
   )
{
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;

   assert(rootexpr != NULL);
   assert(exprintdata != NULL);

   if( *exprintdata == NULL )
   {
      *exprintdata = new SCIP_EXPRINTDATA();
      assert(*exprintdata != NULL);
   }
   else
   {
      (*exprintdata)->need_retape = true;
   }

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, rootexpr, SCIP_EXPRITER_DFS, FALSE) );

   std::set<int> varidxs;
   for( expr = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
   {
      /* cannot handle var-expressions in exprint so far, should be varidx expressions */
      assert(!SCIPisExprVar(scip, expr));

      if( SCIPisExprVaridx(scip, expr) )
         varidxs.insert(SCIPgetIndexExprVaridx(expr));

      //if( strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), "xyz") == 0 )
      //   (*exprintdata)->need_retape_always = true;

      // FIXME data->userevalcapability &= SCIPexprGetUserEvalCapability(expr);
   }

   SCIPfreeExpriter(&it);

   (*exprintdata)->varidxs.reserve(varidxs.size());
   (*exprintdata)->varidxs.insert((*exprintdata)->varidxs.begin(), varidxs.begin(), varidxs.end());

   size_t n = (*exprintdata)->varidxs.size();
   (*exprintdata)->X.resize(n);
   (*exprintdata)->x.resize(n);
   (*exprintdata)->Y.resize(1);

   return SCIP_OKAY;
}

/** frees interpreter data for expression */
SCIP_RETCODE SCIPexprintFreeData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRINTDATA**    exprintdata         /**< pointer to pointer to compiled data to be freed */
   )
{
   assert( exprintdata != NULL);
   assert(*exprintdata != NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &(*exprintdata)->hesrowidxs, (*exprintdata)->hesnnz);
   SCIPfreeBlockMemoryArrayNull(scip, &(*exprintdata)->hescolidxs, (*exprintdata)->hesnnz);
   SCIPfreeBlockMemoryArrayNull(scip, &(*exprintdata)->hesvalues, (*exprintdata)->hesnnz);

   delete *exprintdata;
   *exprintdata = NULL;

   return SCIP_OKAY;
}

/** gives the capability to evaluate an expression by the expression interpreter
 *
 * In cases of user-given expressions, higher order derivatives may not be available for the user-expression,
 * even if the expression interpreter could handle these. This method allows to recognize that, e.g., the
 * Hessian for an expression is not available because it contains a user expression that does not provide
 * Hessians.
 */
SCIP_EXPRINTCAPABILITY SCIPexprintGetExprCapability(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRINTDATA*     exprintdata         /**< interpreter-specific data for expression */
   )
{
   assert(exprintdata != NULL);

   return exprintdata->userevalcapability;
}/*lint !e715*/

/** evaluates an expression tree */
SCIP_RETCODE SCIPexprintEval(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRINTDATA*     exprintdata,        /**< interpreter-specific data for expression */
   SCIP_Real*            varvals,            /**< values of variables */
   SCIP_Real*            val                 /**< buffer to store value of expression */
   )
{
   assert(expr    != NULL);
   assert(exprintdata != NULL);
   assert(varvals != NULL);
   assert(val     != NULL);

   size_t n = exprintdata->varidxs.size();

   if( n == 0 )
   {
      SCIP_CALL( SCIPevalExpr(scip, expr, NULL, 0L) );
      exprintdata->val = *val = SCIPexprGetEvalValue(expr);
      return SCIP_OKAY;
   }

   if( exprintdata->need_retape_always || exprintdata->need_retape )
   {
      for( size_t i = 0; i < n; ++i )
      {
         int idx = exprintdata->varidxs[i];
         exprintdata->X[i] = varvals[idx];
         exprintdata->x[i] = varvals[idx];  /* need this for a following grad or hessian eval with new_x = false */
      }

      CppAD::Independent(exprintdata->X);

      SCIP_CALL( eval(scip, expr, exprintdata, exprintdata->X, exprintdata->Y[0]) );

      exprintdata->f.Dependent(exprintdata->X, exprintdata->Y);

      exprintdata->val = Value(exprintdata->Y[0]);
      SCIPdebugMessage("Eval retaped and computed value %g\n", exprintdata->val);

      // the following is required if the gradient shall be computed by a reverse sweep later
      // exprintdata->val = exprintdata->f.Forward(0, exprintdata->x)[0];

      exprintdata->need_retape = false;
   }
   else
   {
      assert(exprintdata->x.size() >= n);
      for( size_t i = 0; i < n; ++i )
         exprintdata->x[i] = varvals[exprintdata->varidxs[i]];

      exprintdata->val = exprintdata->f.Forward(0, exprintdata->x)[0];
      SCIPdebugMessage("Eval used forward sweep to compute value %g\n", exprintdata->val);
   }

   *val = exprintdata->val;

   return SCIP_OKAY;
}

/** computes value and gradient of an expression tree */
SCIP_RETCODE SCIPexprintGrad(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRINTDATA*     exprintdata,        /**< interpreter-specific data for expression */
   SCIP_Real*            varvals,            /**< values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable values changed since last call to a point evaluation routine? */
   SCIP_Real*            val,                /**< buffer to store expression value */
   SCIP_Real*            gradient            /**< buffer to store expression gradient */
   )
{
   assert(expr     != NULL);
   assert(exprintdata != NULL);
   assert(varvals  != NULL || new_varvals == FALSE);
   assert(val      != NULL);
   assert(gradient != NULL);

   if( new_varvals )
   {
      SCIP_CALL( SCIPexprintEval(scip, exprint, expr, exprintdata, varvals, val) );
   }
   else
      *val = exprintdata->val;

   size_t n = exprintdata->varidxs.size();

   if( n == 0 )
      return SCIP_OKAY;

   vector<double> jac(exprintdata->f.Jacobian(exprintdata->x));

   for( size_t i = 0; i < n; ++i )
      gradient[exprintdata->varidxs[i]] = jac[i];

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "Grad for ");
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\nx    = ");
   for( size_t i = 0; i < n; ++i )
   {
      SCIPinfoMessage(scip, NULL, "\t %g", exprintdata->x[i]);
   }
   SCIPinfoMessage(scip, NULL, "\ngrad = ");
   for( size_t i = 0; i < n; ++i )
   {
      SCIPinfoMessage(scip, NULL, "\t %g", jac[i]);
   }
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   return SCIP_OKAY;
}

/** gives sparsity pattern of lower-triangular part of hessian
 *
 * Since the AD code might need to do a forward sweep, you should pass variable values in here.
 *
 * Result will have (*colidxs)[i] <= (*rowidixs)[i] for i=0..*nnz.
 */
SCIP_RETCODE SCIPexprintHessianSparsity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRINTDATA*     exprintdata,        /**< interpreter-specific data for expression */
   SCIP_Real*            varvals,            /**< values of variables */
   int**                 rowidxs,            /**< buffer to return array with row indices of Hessian elements */
   int**                 colidxs,            /**< buffer to return array with column indices of Hessian elements */
   int*                  nnz                 /**< buffer to return length of arrays */
   )
{
   assert(expr != NULL);
   assert(exprintdata != NULL);
   assert(varvals != NULL);
   assert(rowidxs != NULL);
   assert(colidxs != NULL);
   assert(nnz != NULL);

   if( exprintdata->hesrowidxs == NULL )
   {
      assert(exprintdata->hescolidxs == NULL);
      assert(exprintdata->hesvalues == NULL);
      assert(exprintdata->hesnnz == 0);

      size_t n = exprintdata->varidxs.size();
      if( n == 0 )
      {
         *nnz = 0;
         return SCIP_OKAY;
      }

      size_t nn = n*n;

      if( exprintdata->need_retape_always )
      {
         // pretend dense
         // @todo can we do something better here, e.g., by looking at the expression tree by ourself?

         exprintdata->hesnnz = (n * (n+1))/2;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &exprintdata->hesrowidxs, exprintdata->hesnnz) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &exprintdata->hescolidxs, exprintdata->hesnnz) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &exprintdata->hesvalues, exprintdata->hesnnz) );

         int k = 0;
         for( size_t i = 0; i < n; ++i )
            for( size_t j = 0; j <= i; ++j )
            {
               exprintdata->hesrowidxs[k] = exprintdata->varidxs[i];
               exprintdata->hescolidxs[k] = exprintdata->varidxs[j];
               k++;
            }

#ifdef SCIP_DEBUG
         SCIPinfoMessage(scip, NULL, "HessianSparsity for ");
         SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
         SCIPinfoMessage(scip, NULL, ": all elements, due to discontinuities\n");
#endif

         *rowidxs = exprintdata->hesrowidxs;
         *colidxs = exprintdata->hescolidxs;
         *nnz = exprintdata->hesnnz;

         return SCIP_OKAY;
      }

      if( exprintdata->need_retape )
      {
         SCIP_Real val;
         SCIP_CALL( SCIPexprintEval(scip, exprint, expr, exprintdata, varvals, &val) );
      }  /*lint !e438*/

      SCIPdebugMessage("calling ForSparseJac\n");

      vector<bool> r(nn, false);
      for( size_t i = 0; i < n; ++i )
         r[i*n+i] = true;
      (void) exprintdata->f.ForSparseJac(n, r); // need to compute sparsity for Jacobian first

      SCIPdebugMessage("calling RevSparseHes\n");

      // TODO check whether CppADs sparse hessian can be useful
      vector<bool> s(1, true);
      vector<bool> sparsehes(exprintdata->f.RevSparseHes(n, s));

      for( size_t i = 0; i < nn; ++i )
         if( sparsehes[i] )
         {
            size_t row = i / n;
            size_t col = i % n;
            if( col > row )
               continue;
            ++exprintdata->hesnnz;
         }

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &exprintdata->hesrowidxs, exprintdata->hesnnz) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &exprintdata->hescolidxs, exprintdata->hesnnz) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &exprintdata->hesvalues, exprintdata->hesnnz) );

      int j = 0;
      for( size_t i = 0; i < nn; ++i )
         if( sparsehes[i] )
         {
            size_t row = i / n;
            size_t col = i % n;
            if( col > row )
               continue;
            exprintdata->hesrowidxs[j] = exprintdata->varidxs[row];
            exprintdata->hescolidxs[j] = exprintdata->varidxs[col];
            ++j;
         }

#ifdef SCIP_DEBUG
      SCIPinfoMessage(scip, NULL, "HessianSparsity for ");
      SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, ":");
      for( int i = 0; i < exprintdata->hesnnz; ++i )
      {
         SCIPinfoMessage(scip, NULL, " (%d,%d)", exprintdata->hesrowidxs[i], exprintdata->hescolidxs[i]);
      }
      SCIPinfoMessage(scip, NULL, "\n");
#endif
   }

   *rowidxs = exprintdata->hesrowidxs;
   *colidxs = exprintdata->hescolidxs;
   *nnz = exprintdata->hesnnz;

   return SCIP_OKAY;
}

/** computes value and hessian of an expression
 *
 * Returned arrays rowidxs and colidxs and number of elements nnz are the same as given by SCIPexprintHessianSparsity().
 * Returned array hessianvals will contain the corresponding Hessian elements.
 */
SCIP_RETCODE SCIPexprintHessian(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRINTDATA*     exprintdata,        /**< interpreter-specific data for expression */
   SCIP_Real*            varvals,            /**< values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable values changed since last call to an evaluation routine? */
   SCIP_Real*            val,                /**< buffer to store function value */
   int**                 rowidxs,            /**< buffer to return array with row indices of Hessian elements */
   int**                 colidxs,            /**< buffer to return array with column indices of Hessian elements */
   SCIP_Real**           hessianvals,        /**< buffer to return array with Hessian elements */
   int*                  nnz                 /**< buffer to return length of arrays */
   )
{
   assert(expr != NULL);
   assert(exprintdata != NULL);

   if( exprintdata->hesrowidxs == NULL )
   {
      /* setup sparsity if not done yet */
      int dummy1;
      int* dummy2;

      assert(exprintdata->hescolidxs == NULL);
      assert(exprintdata->hesvalues == NULL);
      assert(exprintdata->hesnnz == 0);

      SCIP_CALL( SCIPexprintHessianSparsity(scip, exprint, expr, exprintdata, varvals, &dummy2, &dummy2, &dummy1) );

      new_varvals = FALSE;
   }

   if( new_varvals )
   {
      SCIP_CALL( SCIPexprintEval(scip, exprint, expr, exprintdata, varvals, val) );
   }
   else
      *val = exprintdata->val;

   size_t n = exprintdata->varidxs.size();

   if( n > 0 )
   {
      // TODO check whether CppADs sparse hessian can be useful
      /* this one uses reverse mode */
      vector<double> hess(exprintdata->f.Hessian(exprintdata->x, 0));

      for( int i = 0; i < exprintdata->hesnnz; ++i )
      {
         exprintdata->hesvalues[i] = hess[exprintdata->getVarPos(exprintdata->hesrowidxs[i]) * n + exprintdata->getVarPos(exprintdata->hescolidxs[i])];
      }
   }

#ifdef SCIP_DEBUG
      SCIPinfoMessage(scip, NULL, "Hessian for ");
      SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\nat x = ");
      for( size_t i = 0; i < n; ++i )
      {
         SCIPinfoMessage(scip, NULL, "\t %g", exprintdata->x[i]);
      }
      SCIPinfoMessage(scip, NULL, "\nis ");
      for( int i = 0; i < exprintdata->hesnnz; ++i )
      {
         SCIPinfoMessage(scip, NULL, " (%d,%d)=%g", exprintdata->hesrowidxs[i], exprintdata->hescolidxs[i], exprintdata->hesvalues[i]);
      }
      SCIPinfoMessage(scip, NULL, "\n");
#endif

   *rowidxs = exprintdata->hesrowidxs;
   *colidxs = exprintdata->hescolidxs;
   *hessianvals = exprintdata->hesvalues;
   *nnz = exprintdata->hesnnz;

   return SCIP_OKAY;
}

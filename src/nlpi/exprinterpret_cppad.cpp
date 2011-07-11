/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    exprinterpret_cppad.cpp
 * @brief   methods to interpret (evaluate) an expression tree "fast" using CppAD
 * @ingroup EXPRINTS
 * @author  Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "nlpi/pub_expr.h"
#include "nlpi/exprinterpret.h"

#include <cmath>
#include <vector>
using std::vector;

/** sign of a value (-1 or +1)
 * 
 * 0.0 has sign +1
 */
#define SIGN(x) ((x) >= 0.0 ? 1.0 : -1.0)

/* in order to use intervals as operands in CppAD,
 * we need to include the intervalarith.hpp very early and require the interval operations to be in the CppAD namespace */
#define SCIPInterval_NAMESPACE CppAD
#include "nlpi/intervalarith.h"

SCIP_Real SCIPInterval_NAMESPACE::SCIPInterval::infinity = SCIP_DEFAULT_INFINITY;
using SCIPInterval_NAMESPACE::SCIPInterval;

#include <cppad/cppad.hpp>
#ifndef CPPAD_PACKAGE_STRING
#include <cppad/config.h>
#define CPPAD_PACKAGE_STRING PACKAGE_STRING
#endif
#include <cppad/declare.hpp>
#include <cppad/error_handler.hpp>

/* Brad recomends using the discrete function feature of CppAD for sign, since it avoids the need for retaping
 * It can be used since it's derivative is almost everywhere 0.0 */

/* sign as function for double */
double sign(const double &x)
{
   return SIGN(x);
}
/* discrete CppAD function sign(double) for use in eval */
CPPAD_DISCRETE_FUNCTION(double, sign)

/* sign as function for SCIPInterval
 * this time outside of the CppAD namespace
 */
SCIPInterval sign(const SCIPInterval& x)
{
   SCIPInterval resultant;

   SCIPintervalSign(&resultant, x);

   return resultant;
}
/* discrete CppAD function sign(SCIPInterval) for use in eval */
CPPAD_DISCRETE_FUNCTION(SCIPInterval, sign)

/** defintion of CondExpOp for SCIPInterval (required by CppAD) */
inline
SCIPInterval CondExpOp(
   enum CppAD::CompareOp cop,
   const SCIPInterval&   left,
   const SCIPInterval&   right,
   const SCIPInterval&   trueCase,
   const SCIPInterval&   falseCase)
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "SCIPInterval CondExpOp(...)",
      "Error: cannot use CondExp with an interval type"
   );

   return SCIPInterval();
}

/** another function that returns whether two intervals are the same (required by CppAD) */
inline
bool EqualOpSeq(
   const SCIPInterval&   x,                  /**< first operand */
   const SCIPInterval&   y                   /**< second operand */
   )
{
   return x == y;
}

/** another function required by CppAD */
inline
bool IdenticalPar(
   const SCIPInterval&   x                   /**< operand */
   )
{
   return true;
}

/** returns whether the interval equals [0,0] */
inline
bool IdenticalZero(
   const SCIPInterval&   x                   /**< operand */
   )
{
   return (x == 0.0);
}

/** returns whether the interval equals [1,1] */
inline
bool IdenticalOne(
   const SCIPInterval&   x                   /**< operand */
   )
{
   return (x == 1.0);
}

/** yet another function that checks whether two intervals are equal */
inline
bool IdenticalEqualPar(
   const SCIPInterval&   x,                  /**< first operand */
   const SCIPInterval&   y                   /**< second operand */
   )
{
   return (x == y);
}

/** greater than zero not defined for intervals */
inline
bool GreaterThanZero(
   const SCIPInterval&   x                   /**< operand */
   )
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "GreaterThanZero(x)",
      "Error: cannot use GreaterThanZero with interval"
   );

   return false;
}

/** greater than or equal zero not defined for intervals */
inline
bool GreaterThanOrZero(
   const SCIPInterval&   x                   /**< operand */
   )
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__ ,
      "GreaterThanOrZero(x)",
      "Error: cannot use GreaterThanOrZero with interval"
   );

   return false;
}

/** less than not defined for intervals */
inline
bool LessThanZero(
   const SCIPInterval&   x                   /**< operand */
)
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "LessThanZero(x)",
      "Error: cannot use LessThanZero with interval"
   );

   return false;
}

/** less than or equal not defined for intervals */
inline
bool LessThanOrZero(
   const SCIPInterval&   x                   /**< operand */
)
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "LessThanOrZero(x)",
      "Error: cannot use LessThanOrZero with interval"
   );

   return false;
}

/** conversion to integers not defined for intervals */
inline
int Integer(
   const SCIPInterval&   x                   /**< operand */
)
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "Integer(x)",
      "Error: cannot use Integer with interval"
   );

   return 0;
}

/** printing of an interval (required by CppAD) */
inline
std::ostream& operator<<(std::ostream& out, const SCIP_INTERVAL& x)
{
   out << '[' << x.inf << ',' << x.sup << ']';
   return out;
}

using CppAD::AD;

/** expression interpreter */
struct SCIP_ExprInt
{
   BMS_BLKMEM*           blkmem;             /**< block memory data structure */
};

/** expression specific interpreter data */
class SCIP_ExprIntData
{
public:
   /* constructor */
   SCIP_ExprIntData()
   : need_retape(true), int_need_retape(true), need_retape_always(false), blkmem(NULL), root(NULL)
   { }

   /* destructor */
   ~SCIP_ExprIntData()
   { }

   vector< AD<double> >  X;                  /**< vector of dependent variables */
   vector< AD<double> >  Y;                  /**< result vector */ 
   CppAD::ADFun<double>  f;                  /**< the function to evaluate as CppAD object */

   vector<double>        x;                  /**< current values of dependent variables */
   double                val;                /**< current function value */
   bool                  need_retape;        /**< will retaping be required for the next point evaluation? */

   vector< AD<SCIPInterval> > int_X;         /**< interval vector of dependent variables */
   vector< AD<SCIPInterval> > int_Y;         /**< interval result vector */
   CppAD::ADFun<SCIPInterval> int_f;         /**< the function to evaluate on intervals as CppAD object */

   vector<SCIPInterval>  int_x;              /**< current interval values of dependent variables */
   SCIPInterval          int_val;            /**< current interval function value */
   bool                  int_need_retape;    /**< will retaping be required for the next interval evaluation? */

   bool                  need_retape_always; /**< will retaping be always required? */

   BMS_BLKMEM*           blkmem;             /**< block memory used to allocate expresstion tree */
   SCIP_EXPR*            root;               /**< copy of expression tree; @todo do we really need to make a copy? */
};

/** template for evaluation for signpower operator
 * only implemented for real numbers, thus gives error by default
 */
template<class Type>
void evalSignPower(
   Type&                 resultant,          /**< resultant */
   Type&                 arg1,               /**< first operand */
   Type&                 arg2                /**< second operand */
   )
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "evalSignPower()",
      "Error: SignPower not implemented for this value type"
   );
}

/** specialization of signpower evaluation for real numbers
 */
template<>
void evalSignPower(
   CppAD::AD<double>&    resultant,          /**< resultant */
   CppAD::AD<double>&    arg1,               /**< first operand */
   CppAD::AD<double>&    arg2                /**< second operand */
   )
{
   if( arg1 == 0.0 )
      resultant = 0.0;
   else if( arg1 > 0.0 )
      resultant =  pow( arg1, arg2);
   else
      resultant = -pow(-arg1, arg2);
}

/** template for evaluation for minimum operator
 * only implemented for real numbers, thus gives error by default
 */
template<class Type>
void evalMin(
   Type&                 resultant,          /**< resultant */
   Type&                 arg1,               /**< first operand */
   Type&                 arg2                /**< second operand */
   )
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "evalMin()",
      "Error: Min not implemented for this value type"
   );
}

/** specialization of minimum evaluation for real numbers
 */
template<>
void evalMin(
   CppAD::AD<double>&    resultant,          /**< resultant */
   CppAD::AD<double>&    arg1,               /**< first operand */
   CppAD::AD<double>&    arg2                /**< second operand */
   )
{
   resultant = MIN(arg1, arg2);
}
/** template for evaluation for maximum operator
 * only implemented for real numbers, thus gives error by default
 */
template<class Type>
void evalMax(
   Type&                 resultant,          /**< resultant */
   Type&                 arg1,               /**< first operand */
   Type&                 arg2                /**< second operand */
   )
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "evalMax()",
      "Error: Max not implemented for this value type"
   );
}

/** specialization of maximum evaluation for real numbers
 */
template<>
void evalMax(
   CppAD::AD<double>&    resultant,          /**< resultant */
   CppAD::AD<double>&    arg1,               /**< first operand */
   CppAD::AD<double>&    arg2                /**< second operand */
   )
{
   resultant = MAX(arg1, arg2);
}


/** template for evaluation for square operator
 * default is to multiply arg with itself
 */
template<class Type>
void evalSquare(
   Type&                 resultant,          /**< resultant */
   Type&                 arg                 /**< operand */
   )
{
   resultant = arg * arg;
}

#if 0 /* @todo find out how to implement a working version of evalSquare for intervals in CppAD */
/** specialization of square evaluation for intervals
 * for intervals, we can get tighter results if we do not just multiply the argument with itself
 */
template<>
void evalSquare(
   CppAD::AD<SCIPInterval>& resultant,          /**< resultant */
   CppAD::AD<SCIPInterval>& arg                 /**< operand */
   )
{
   SCIPInterval result;

   SCIPintervalSquare(SCIPInterval::infinity, &result, Value(arg));

   resultant = result;
}
#endif

/** template for evaluation for square-root operator
 * default is to use the standard sqrt-function
 */
template<class Type>
void evalSqrt(
   Type&                 resultant,          /**< resultant */
   Type&                 arg                 /**< operand */
   )
{
   resultant = sqrt(arg);
}

/** specialization of square-root operator for numbers
 * we perturb the function a little bit so that it's derivatives are defined in 0.0
 */
template<>
void evalSqrt(
   CppAD::AD<double>&    resultant,          /**< resultant */
   CppAD::AD<double>&    arg                 /**< operand */
   )
{
   resultant = sqrt(arg + 1e-20) - 1e-10;
}

/** template for function that sets a value to NaN
 * default is to set it to NAN, if available, and to log(-1.0) otherwise
 */
template<class Type>
void setToNaN(
   Type&                 resultant           /**< resultant */
   )
{
#ifdef NAN
   resultant = NAN;
#else
   resultant = log(-1.0);
#endif
}

/** specialization of setNaN for intervals
 * for intervals, we set the interval to the empty interval
 */
template<>
void setToNaN(
   SCIPInterval&         resultant           /**< resultant */
   )
{
   SCIPintervalSetEmpty(&resultant);
}

/** CppAD compatible evaluation of an expression for given arguments and parameters */
template<class Type>
SCIP_RETCODE eval(
   SCIP_EXPR*            expr,               /**< expression */
   const vector<Type>&   x,                  /**< values of variables */
   SCIP_Real*            param,              /**< values of parameters */
   Type&                 val                 /**< buffer to store expression value */
   )
{
   Type* buf;
   
   assert(expr != NULL);

   /* todo use SCIP_MAXCHILD_ESTIMATE as in expression.c */

   buf = NULL;
   if( SCIPexprGetNChildren(expr) )
   {
      if( BMSallocMemoryArray(&buf, SCIPexprGetNChildren(expr)) == NULL )
         return SCIP_NOMEMORY;

      for( int i = 0; i < SCIPexprGetNChildren(expr); ++i )
         SCIP_CALL( eval(SCIPexprGetChildren(expr)[i], x, param, buf[i]) );
   }

   switch(SCIPexprGetOperator(expr))
   {
      case SCIP_EXPR_VARIDX:
         assert(SCIPexprGetOpIndex(expr) < (int)x.size());
         val = x[SCIPexprGetOpIndex(expr)];
         break;

      case SCIP_EXPR_CONST:
         val = SCIPexprGetOpReal(expr);
         break;

      case SCIP_EXPR_PARAM:
         assert(param != NULL);
         val = param[SCIPexprGetOpIndex(expr)];
         break;

      case SCIP_EXPR_PLUS:
         val = buf[0] + buf[1];
         break;

      case SCIP_EXPR_MINUS:
         val = buf[0] - buf[1];
         break;

      case SCIP_EXPR_MUL:
         val = buf[0] * buf[1];
         break;

      case SCIP_EXPR_DIV:
         val = buf[0] / buf[1];
         break;

      case SCIP_EXPR_SQUARE:
         evalSquare(val, buf[0]);
         break;

      case SCIP_EXPR_SQRT:
         evalSqrt(val, buf[0]);
         break;

      case SCIP_EXPR_POWER:
         val = pow(buf[0], buf[1]);
         break;

      case SCIP_EXPR_EXP:
         val = exp(buf[0]);
         break;

      case SCIP_EXPR_LOG:
         val = log(buf[0]);
         break;

      case SCIP_EXPR_SIN:
         val = sin(buf[0]);
         break;

      case SCIP_EXPR_COS:
         val = cos(buf[0]);
         break;

      case SCIP_EXPR_TAN:
         val = tan(buf[0]);
         break;
#if 0
      case SCIP_EXPR_ERF:
         val = erf(buf[0]);
         break;

      case SCIP_EXPR_ERFI:
         return SCIP_ERROR;
#endif
      case SCIP_EXPR_MIN:
         evalMin(val, buf[0], buf[1]);
         break;

      case SCIP_EXPR_MAX:
         evalMax(val, buf[0], buf[1]);
         break;

      case SCIP_EXPR_ABS:
         val = abs(buf[0]);
         break;

      case SCIP_EXPR_SIGN:
         val = sign(buf[0]);
         break;

      case SCIP_EXPR_SIGNPOWER:
         evalSignPower(val, buf[0], buf[1]);
         break;

      case SCIP_EXPR_INTPOWER:
         val = pow(buf[0], SCIPexprGetIntPowerExponent(expr));
         break;

      case SCIP_EXPR_SUM:
         val = 0.0;
         for (int i = 0; i < SCIPexprGetNChildren(expr); ++i)
            val += buf[i];
         break;

      case SCIP_EXPR_PRODUCT:
         val = 1.0;
         for (int i = 0; i < SCIPexprGetNChildren(expr); ++i)
            val *= buf[i];
         break;

      case SCIP_EXPR_LINEAR:
      {
         SCIP_Real* coefs;

         coefs = SCIPexprGetLinearCoefs(expr);
         assert(coefs != NULL || SCIPexprGetNChildren(expr) == 0);

         val = SCIPexprGetLinearConstant(expr);
         for (int i = 0; i < SCIPexprGetNChildren(expr); ++i)
            val += coefs[i] * buf[i];
         break;
      }
      
      case SCIP_EXPR_QUADRATIC:
      {
         SCIP_QUADELEM* quadelems;
         int nquadelems;
         
         nquadelems = SCIPexprGetNQuadElements(expr);
         quadelems  = SCIPexprGetQuadElements(expr);
         assert(quadelems != NULL || nquadelems == 0);
         
         val = 0.0;
         for (int i = nquadelems; i > 0; --i, ++quadelems)
         {
            if( quadelems->idx1 == quadelems->idx2 )
            {
               Type tmp;
               evalSquare(tmp, buf[quadelems->idx1]);
               val += quadelems->coef * tmp;
            }
            else
            {
               val += quadelems->coef * buf[quadelems->idx1] * buf[quadelems->idx2];
            }
         }
         
         break;
      }

      case SCIP_EXPR_POLYNOMIAL:
      {
         SCIP_EXPRDATA_MONOMIAL** monomials;
         Type childval;
         Type monomialval;
         SCIP_Real exponent;
         int nmonomials;
         int nfactors;
         int* childidxs;
         SCIP_Real* exponents;
         int i;
         int j;

         val = SCIPexprGetPolynomialConstant(expr);

         nmonomials = SCIPexprGetNMonomials(expr);
         monomials  = SCIPexprGetMonomials(expr);

         for( i = 0; i < nmonomials; ++i )
         {
            nfactors  = SCIPexprGetMonomialNFactors(monomials[i]);
            childidxs = SCIPexprGetMonomialChildIndices(monomials[i]);
            exponents = SCIPexprGetMonomialExponents(monomials[i]);
            monomialval  = SCIPexprGetMonomialCoef(monomials[i]);

            for( j = 0; j < nfactors; ++j )
            {
               assert(childidxs[j] >= 0);
               assert(childidxs[j] <  SCIPexprGetNChildren(expr));

               childval = buf[childidxs[j]];
               exponent = exponents[j];

               /* cover some special exponents separately to avoid calling expensive pow function */
               if( exponent == 0.0 )
                  continue;
               if( exponent == 1.0 )
               {
                  monomialval *= childval;
                  continue;
               }
               if( exponent == 2.0 )
               {
                  Type tmp;
                  evalSquare(tmp, childval);
                  monomialval *= tmp;
                  continue;
               }
               if( exponent == 0.5 )
               {
                  Type tmp;
                  evalSqrt(tmp, childval);
                  monomialval *= tmp;
                  continue;
               }
               if( exponent == -1.0 )
               {
                  monomialval /= childval;
                  continue;
               }
               if( exponent == -2.0 )
               {
                  Type tmp;
                  evalSquare(tmp, childval);
                  monomialval /= tmp;
                  continue;
               }
               monomialval *= pow(childval, exponent);
            }

            val += monomialval;
         }

         break;
      }

      default:
         return SCIP_ERROR;
   }

   BMSfreeMemoryArrayNull(&buf);

   return SCIP_OKAY;
}

/** analysis an expression tree whether it requires retaping on every evaluation
 * this may be the case if the evaluation sequence depends on values of operands (e.g., in case of abs, sign, signpower, ...)
 */
bool needAlwaysRetape(SCIP_EXPR* expr)
{
   assert(expr != NULL);
   assert(SCIPexprGetChildren(expr) != NULL || SCIPexprGetNChildren(expr) == 0);

   for( int i = 0; i < SCIPexprGetNChildren(expr); ++i )
   {
      if( needAlwaysRetape(SCIPexprGetChildren(expr)[i]) )
         return true;
   }

   switch( SCIPexprGetOperator(expr) )
   {
      case SCIP_EXPR_MIN:
      case SCIP_EXPR_MAX:
      case SCIP_EXPR_ABS:
      case SCIP_EXPR_SIGNPOWER:
         return true;

      default: ;
   }

   return false;
}

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
   return SCIP_EXPRINTCAPABILITY_FUNCVALUE | SCIP_EXPRINTCAPABILITY_INTFUNCVALUE |
      SCIP_EXPRINTCAPABILITY_GRADIENT | SCIP_EXPRINTCAPABILITY_INTGRADIENT |
      SCIP_EXPRINTCAPABILITY_HESSIAN;
}

/** creates an expression interpreter object */
SCIP_RETCODE SCIPexprintCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPRINT**        exprint             /**< buffer to store pointer to expression interpreter */
)
{
   assert(blkmem  != NULL);
   assert(exprint != NULL);
   
   if( BMSallocMemory(exprint) == NULL )
      return SCIP_NOMEMORY;
   
   (*exprint)->blkmem = blkmem;

   return SCIP_OKAY;
}

/** frees an expression interpreter object */
SCIP_RETCODE SCIPexprintFree(
   SCIP_EXPRINT**        exprint             /**< expression interpreter that should be freed */
)
{
   assert( exprint != NULL);
   assert(*exprint != NULL);
   
   BMSfreeMemory(exprint);

   return SCIP_OKAY;
}

/** compiles an expression tree and stores compiled data in expression tree */
SCIP_RETCODE SCIPexprintCompile(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree                /**< expression tree */
)
{
   assert(tree    != NULL);
   
   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   if (!data)
   {
      data = new SCIP_EXPRINTDATA();
      assert( data != NULL );
      SCIPexprtreeSetInterpreterData(tree, data);
      SCIPdebugMessage("set interpreter data in tree %p to %p\n", (void*)tree, (void*)data);
   }
   else
   {
      data->need_retape     = true;
      data->int_need_retape = true;
   }

   int n = SCIPexprtreeGetNVars(tree);

   data->X.resize(n);
   data->x.resize(n);
   data->Y.resize(1);

   data->int_X.resize(n);
   data->int_x.resize(n);
   data->int_Y.resize(1);

   if( data->root != NULL )
   {
      SCIPexprFreeDeep(exprint->blkmem, &data->root);
   }

   SCIP_EXPR* root = SCIPexprtreeGetRoot(tree);
   
   SCIP_CALL( SCIPexprCopyDeep(exprint->blkmem, &data->root, root) );

   data->need_retape_always = needAlwaysRetape(SCIPexprtreeGetRoot(tree));

   data->blkmem = exprint->blkmem;

   return SCIP_OKAY;
}

/** frees interpreter data */
SCIP_RETCODE SCIPexprintFreeData(
   SCIP_EXPRINTDATA**    interpreterdata     /**< interpreter data that should freed */
)
{
   assert( interpreterdata != NULL);
   assert(*interpreterdata != NULL);

   if( (*interpreterdata)->root != NULL )
      SCIPexprFreeDeep((*interpreterdata)->blkmem, &(*interpreterdata)->root);   

   delete *interpreterdata;
   *interpreterdata = NULL; 

   return SCIP_OKAY;
}

/** notify expression interpreter that a new parameterization is used
 * this probably causes retaping by AD algorithms
 */
SCIP_RETCODE SCIPexprintNewParametrization(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree                /**< expression tree */
)
{
   assert(exprint != NULL);
   assert(tree    != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   if( data != NULL )
   {
      data->need_retape     = true;
      data->int_need_retape = true;
   }
	 
   return SCIP_OKAY;
}

/** evaluates an expression tree */
SCIP_RETCODE SCIPexprintEval(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables */
   SCIP_Real*            val                 /**< buffer to store value */
)
{
   SCIP_EXPRINTDATA* data;

   assert(exprint != NULL);
   assert(tree    != NULL);
   assert(varvals != NULL);
   assert(val     != NULL);

   data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);
   assert(SCIPexprtreeGetNVars(tree) == (int)data->X.size());
   assert(SCIPexprtreeGetRoot(tree)  != NULL);

   int n = SCIPexprtreeGetNVars(tree);

   if( data->need_retape_always || data->need_retape )
   {
      for( int i = 0; i < n; ++i )
      {
         data->X[i] = varvals[i];
         data->x[i] = varvals[i];
      }

      CppAD::Independent(data->X);

      if( data->root != NULL )
         SCIP_CALL( eval(data->root, data->X, SCIPexprtreeGetParamVals(tree), data->Y[0]) );
      else
         data->Y[0] = 0.0;

      data->f.Dependent(data->X, data->Y);

      data->val = Value(data->Y[0]);
      SCIPdebugMessage("Eval retaped and computed value %g\n", data->val);

      data->need_retape = false;
   }
   else
   {
      assert((int)data->x.size() >= n);
      for( int i = 0; i < n; ++i )
         data->x[i] = varvals[i];

      data->val = data->f.Forward(0, data->x)[0];
      SCIPdebugMessage("Eval used foward sweep to compute value %g\n", data->val);
   }

   *val = data->val;

   return SCIP_OKAY;
}

/** evaluates an expression tree on intervals */
extern
SCIP_RETCODE SCIPexprintEvalInt(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        varvals,            /**< interval values of variables */
   SCIP_INTERVAL*        val                 /**< buffer to store interval value of expression */
)
{
   SCIP_EXPRINTDATA* data;

   assert(exprint != NULL);
   assert(tree    != NULL);
   assert(varvals != NULL);
   assert(val     != NULL);

   data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);
   assert(SCIPexprtreeGetNVars(tree) == (int)data->int_X.size());
   assert(SCIPexprtreeGetRoot(tree)  != NULL);

   int n = SCIPexprtreeGetNVars(tree);

   SCIPInterval::infinity = infinity;

   if( data->int_need_retape || data->need_retape_always )
   {
      for( int i = 0; i < n; ++i )
      {
         data->int_X[i] = varvals[i];
         data->int_x[i] = varvals[i];
      }

      CppAD::Independent(data->int_X);

      if( data->root != NULL )
         SCIP_CALL( eval(data->root, data->int_X, SCIPexprtreeGetParamVals(tree), data->int_Y[0]) );
      else
         data->int_Y[0] = 0.0;

      data->int_f.Dependent(data->int_X, data->int_Y);

      data->int_val = Value(data->int_Y[0]);

      data->int_need_retape = false;
   }
   else
   {
      assert((int)data->int_x.size() >= n);
      for( int i = 0; i < n; ++i )
         data->int_x[i] = varvals[i];

      data->int_val = data->int_f.Forward(0, data->int_x)[0];
   }

   *val = data->int_val;

   return SCIP_OKAY;
}

/** computes value and gradient of an expression tree */
SCIP_RETCODE SCIPexprintGrad(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable values changed since last call to a point evaluation routine? */
   SCIP_Real*            val,                /**< buffer to store expression value */
   SCIP_Real*            gradient            /**< buffer to store expression gradient, need to have length at least SCIPexprtreeGetNVars(tree) */
)
{
   assert(exprint  != NULL);
   assert(tree     != NULL);
   assert(varvals  != NULL || new_varvals == FALSE);
   assert(val      != NULL);
   assert(gradient != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);

   if( new_varvals )
   {
      SCIP_CALL( SCIPexprintEval(exprint, tree, varvals, val) );
   }
   else
      *val = data->val;

   int n = SCIPexprtreeGetNVars(tree);

   vector<double> jac(data->f.Jacobian(data->x));

   for( int i = 0; i < n; ++i )
      gradient[i] = jac[i];

#ifdef SCIP_DEBUG
   SCIPdebugMessage("Grad for "); SCIPexprtreePrint(tree, NULL, NULL, NULL); printf("\n");
   SCIPdebugMessage("x    ="); for (int i = 0; i < n; ++i) printf("\t %g", data->x[i]); printf("\n");
   SCIPdebugMessage("grad ="); for (int i = 0; i < n; ++i) printf("\t %g", gradient[i]); printf("\n");
#endif

   return SCIP_OKAY;
}

/** computes interval value and interval gradient of an expression tree */
SCIP_RETCODE SCIPexprintGradInt(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        varvals,            /**< interval values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable interval values changed since last call to an interval evaluation routine? */
   SCIP_INTERVAL*        val,                /**< buffer to store expression interval value */
   SCIP_INTERVAL*        gradient            /**< buffer to store expression interval gradient, need to have length at least SCIPexprtreeGetNVars(tree) */
)
{
   assert(exprint  != NULL);
   assert(tree     != NULL);
   assert(varvals  != NULL || new_varvals == false);
   assert(val      != NULL);
   assert(gradient != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);

   if (new_varvals)
      SCIP_CALL( SCIPexprintEvalInt(exprint, tree, infinity, varvals, val) );
   else
      *val = data->int_val;

   int n = SCIPexprtreeGetNVars(tree);

   vector<SCIPInterval> jac(data->int_f.Jacobian(data->int_x));

   for (int i = 0; i < n; ++i)
      gradient[i] = jac[i];

#ifdef SCIP_DEBUG
   SCIPdebugMessage("GradInt for "); SCIPexprtreePrint(tree, NULL, NULL, NULL); printf("\n");
   SCIPdebugMessage("x    ="); for (int i = 0; i < n; ++i) printf("\t [%g,%g]", SCIPintervalGetInf(data->int_x[i]), SCIPintervalGetSup(data->int_x[i])); printf("\n");
   SCIPdebugMessage("grad ="); for (int i = 0; i < n; ++i) printf("\t [%g,%g]", SCIPintervalGetInf(gradient[i]), SCIPintervalGetSup(gradient[i])); printf("\n");
#endif

   return SCIP_OKAY;
}

/** gives sparsity pattern of hessian
 * NOTE: this function might be replaced later by something nicer 
 * Since the AD code might need to do a forward sweep, you should pass variable values in here.
 */
SCIP_RETCODE SCIPexprintHessianSparsityDense(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables */
   SCIP_Bool*            sparsity            /**< buffer to store sparsity pattern of Hessian, sparsity[i+n*j] indicates whether entry (i,j) is nonzero in the hessian */
)
{
   assert(exprint  != NULL);
   assert(tree     != NULL);
   assert(varvals  != NULL);
   assert(sparsity != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);

   int n = SCIPexprtreeGetNVars(tree);
   int nn = n*n;

   if( data->need_retape_always )
   {
      // @todo can we do something better here, e.g., by looking at the expression tree by ourself?

      for( int i = 0; i < nn; ++i )
         sparsity[i] = TRUE;

#ifdef SCIP_DEBUG
      SCIPdebugMessage("HessianSparsityDense for "); SCIPexprtreePrint(tree, NULL, NULL, NULL); printf("\n");
      SCIPdebugMessage("sparsity = all elements, due to discontinuouities\n");
#endif

      return SCIP_OKAY;
   }

   if( data->need_retape )
   {
      SCIP_Real val;
      SCIP_CALL( SCIPexprintEval(exprint, tree, varvals, &val) );
   }

   vector<bool> r(nn, false);
   for (int i = 0; i < n; ++i)
      r[i*n+i] = true;
   data->f.ForSparseJac(n, r); // need to compute sparsity for Jacobian first

   vector<bool> s(1, true);
   vector<bool> sparsehes(data->f.RevSparseHes(n, s));

   for( int i = 0; i < nn; ++i )
      sparsity[i] = sparsehes[i];

#ifdef SCIP_DEBUG
   SCIPdebugMessage("HessianSparsityDense for "); SCIPexprtreePrint(tree, NULL, NULL, NULL); printf("\n");
   SCIPdebugMessage("sparsity ="); for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) if (sparsity[i*n+j]) printf(" (%d,%d)", i, j); printf("\n");
#endif

   return SCIP_OKAY;
}

/** computes value and dense hessian of an expression tree
 * the full hessian is computed (lower left and upper right triangle)
 */
SCIP_RETCODE SCIPexprintHessianDense(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable values changed since last call to an evaluation routine? */
   SCIP_Real*            val,                /**< buffer to store function value */
   SCIP_Real*            hessian             /**< buffer to store hessian values, need to have size at least n*n */
)
{
   assert(exprint != NULL);
   assert(tree    != NULL);
   assert(varvals != NULL || new_varvals == FALSE);
   assert(val     != NULL);
   assert(hessian != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);

   if( new_varvals )
   {
      SCIP_CALL( SCIPexprintEval(exprint, tree, varvals, val) );
   }
   else
      *val = data->val;

   int n = SCIPexprtreeGetNVars(tree);

   vector<double> hess(data->f.Hessian(data->x, 0));

   int nn = n*n;
   for (int i = 0; i < nn; ++i)
      hessian[i] = hess[i];

#ifdef SCIP_DEBUG
   SCIPdebugMessage("HessianDense for "); SCIPexprtreePrint(tree, NULL, NULL, NULL); printf("\n");
   SCIPdebugMessage("x    ="); for (int i = 0; i < n; ++i) printf("\t %g", data->x[i]); printf("\n");
   SCIPdebugMessage("hess ="); for (int i = 0; i < n*n; ++i) printf("\t %g", hessian[i]); printf("\n");
#endif

   return SCIP_OKAY;
}

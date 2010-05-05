/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: exprinterpret_cppad.cpp,v 1.3 2010/05/05 16:20:13 bzfviger Exp $"

/**@file   exprinterpret_cppad.cpp
 * @brief  methods to interpret (evaluate) an expression tree "fast" using CppAD
 * @author Stefan Vigerske
 * 
 * @todo can allow MIN and MAX if we retape for every new x
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "nlpi/expression.h"
#include "nlpi/exprinterpret.h"

/** sign of a value (-1 or +1)
 * 
 * 0.0 has sign +1
 */
#define SIGN(x) ((x) >= 0.0 ? 1.0 : -1.0)

#include <cppad/cppad.hpp>
#include <cppad/config.h>  // to get PACKAGE_STRING define
#include <cppad/declare.hpp>
#include <cppad/error_handler.hpp>

/* sign function, not using macro */
inline double sign(const double& x) { return SIGN(x); }

/* more sign function definitions for CppAD */
namespace CppAD
{
template <class Base>
inline AD<Base> sign(const AD<Base> &x) { return x / x.Abs(); }  //TODO FIXME: this looks awful

template <class Base>
inline AD<Base> sign(const VecAD_reference<Base> &x) { return sign( x.ADBase() ); }
}

#include <vector>

using std::vector;
using CppAD::AD;

struct SCIP_ExprInt
{
   BMS_BLKMEM*           blkmem;             /**< block memory data structure */
};

class SCIP_ExprIntData
{
public:
   /* constructor */
   SCIP_ExprIntData()
   : need_retape(true), need_retape_always(false), blkmem(NULL), root(NULL)
   { }

   /* destructor */
   ~SCIP_ExprIntData()
   { }

   vector< AD<double> >  X;                  /**< vector of dependent variables */
   vector< AD<double> >  Y;                  /**< result vector */ 
   CppAD::ADFun<double>  f;                  /**< the function to evaluate as CppAD object */

   vector<double>        x;                  /**< current values of dependent variables */
   double                val;                /**< current function value */

   bool                  need_retape;        /**< will retaping be required for the next evaluation? */
   bool                  need_retape_always; /**< will retaping be always required? */

   BMS_BLKMEM*           blkmem;             /**< block memory used to allocate expresstion tree */
   SCIP_EXPR*            root;               /**< copy of expression tree; @todo do we really need to make a copy? */
};

template<class Type>
SCIP_RETCODE eval(SCIP_EXPR* expr, const vector<Type>& x, SCIP_Real* param, Type& val)
{
   Type *buf;
   if( BMSallocMemoryArray(&buf, SCIPexprGetNChildren(expr)) == NULL )
      return SCIP_NOMEMORY;

   for( int i = 0; i < SCIPexprGetNChildren(expr); ++i )
      SCIP_CALL( eval(SCIPexprGetChildren(expr)[i], x, param, buf[i]) );

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
         val = buf[0] * buf[0];
         break;

      case SCIP_EXPR_SQRT:
         val = sqrt(buf[0]);
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

      case SCIP_EXPR_ERF:
         val = erf(buf[0]);
         break;

      case SCIP_EXPR_ERFI:
         return SCIP_ERROR;

      case SCIP_EXPR_MIN:
         return SCIP_ERROR;
         //			val = MIN(buf[0], buf[1]);
         break;

      case SCIP_EXPR_MAX:
         return SCIP_ERROR;
         //			val = MAX(buf[0], buf[1]);
         break;

      case SCIP_EXPR_ABS:
         val = abs(buf[0]);
         break;

      case SCIP_EXPR_SIGN:
         val = sign(buf[0]);
         break;

      case SCIP_EXPR_SIGNPOWER:
         if( buf[0] == 0.0 )
            val = 0.0;
         else if( buf[0] > 0.0 )
            val =  pow( buf[0], buf[1]);
         else
            val = -pow(-buf[0], buf[1]);
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

      default:
         return SCIP_ERROR;
   }

   BMSfreeMemoryArray(&buf);

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
      case SCIP_EXPR_SIGN:
      case SCIP_EXPR_SIGNPOWER:
         return true;

      default: ;
   }

   return false;
}

/** gets name and version of expression interpreter */
const char* SCIPexprintGetName(void)
{
   return PACKAGE_STRING;
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
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree                /** expression tree */
)
{
   assert(tree    != NULL);
   
   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   if (!data)
   {
      data = new SCIP_EXPRINTDATA();
      assert( data != NULL );
      SCIPexprtreeSetInterpreterData(tree, data);
      SCIPdebugMessage("set interpreter data in tree %p to %p\n", tree, data);
   }
   else
   {
      data->need_retape     = true;
   }

   int n = SCIPexprtreeGetNVars(tree);

   data->X.resize(n);
   data->x.resize(n);
   data->Y.resize(1);

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
   SCIP_EXPRINTDATA**    interpreterdata     /** interpreter data that should freed */
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
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree                /** expression tree */
)
{
   assert(exprint != NULL);
   assert(tree    != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   if( data != NULL )
      data->need_retape = true;

   return SCIP_OKAY;
}

/** evaluates an expression tree */
SCIP_RETCODE SCIPexprintEval(
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree,               /** expression tree */
   SCIP_Real*            varvals,            /** values of variables */
   SCIP_Real*            val                 /** buffer to store value */
)
{
   assert(exprint != NULL);
   assert(tree    != NULL);
   assert(varvals != NULL);
   assert(val     != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
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

/** gets number of nonzeros in gradient of expression tree */
SCIP_RETCODE SCIPexprintGetNGradPattern(
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree,               /** expression tree */
   int*                  nnz                 /** buffer to store number of nonzeros */
)
{
   SCIPerrorMessage("SCIPexprintGetNGradPattern not implemented for CppAD");
   return SCIP_ERROR;
}

/** gets sparsity pattern of expression trees gradient */
SCIP_RETCODE SCIPexprintGetGradPattern(
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree,               /** expression tree */
   int*                  gradidx             /** buffer to store gradient indices */
)
{
   SCIPerrorMessage("SCIPexprintGetGradPattern not implemented for CppAD");
   return SCIP_ERROR;
}

/** computes value and gradient of an expression tree */
SCIP_RETCODE SCIPexprintGrad(
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree,               /** expression tree */
   SCIP_Real*            varvals,            /** values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /** have variable values changed since last call to an evaluation routine? */
   SCIP_Real*            val,                /** buffer to store value */
   SCIP_Real*            gradval             /** buffer to store gradient values */
)
{
   SCIPerrorMessage("SCIPexprintEvalGrad not implemented for CppAD");
   return SCIP_ERROR;
}

/** computes value and dense gradient of an expression tree */
SCIP_RETCODE SCIPexprintGradDense(
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree,               /** expression tree */
   SCIP_Real*            varvals,            /** values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /** have variable values changed since last call to an evaluation routine? */
   SCIP_Real*            val,                /** buffer to store value */
   SCIP_Real*            gradient            /** buffer to store gradient */
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
   SCIPdebugMessage("GradDense for "); SCIPexprtreePrint(tree, NULL); printf("\n");
   SCIPdebugMessage("x    ="); for (int i = 0; i < n; ++i) printf("\t %g", data->x[i]); printf("\n");
   SCIPdebugMessage("grad ="); for (int i = 0; i < n; ++i) printf("\t %g", gradient[i]); printf("\n");
#endif

   return SCIP_OKAY;
}

/** gives sparsity pattern of hessian
 * NOTE: this function might be replaced later by something nicer 
 * Since the AD code might need to do a forward sweep, you should pass variable values in here.
 */
SCIP_RETCODE SCIPexprintHessianSparsityDense(
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree,               /** expression tree */
   SCIP_Real*            varvals,            /** values of variables */
   SCIP_Bool*            sparsity            /** buffer to store sparsity pattern of Hessian, sparsity[i+n*j] indicates whether entry (i,j) is nonzero in the hessian */ 
)
{
   assert(exprint  != NULL);
   assert(tree     != NULL);
   assert(varvals  != NULL);
   assert(sparsity != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);

   if( data->need_retape )
   {
      SCIP_Real val;
      SCIP_CALL( SCIPexprintEval(exprint, tree, varvals, &val) );
   }

   int n = SCIPexprtreeGetNVars(tree);
   int nn = n*n;

   vector<bool> r(nn, false);
   for (int i = 0; i < n; ++i)
      r[i*n+i] = true;
   data->f.ForSparseJac(n, r); // need to compute sparsity for Jacobian first

   vector<bool> s(1, true);
   vector<bool> sparsehes(data->f.RevSparseHes(n, s));

   for( int i = 0; i < nn; ++i )
      sparsity[i] = sparsehes[i];

#ifdef SCIP_DEBUG
   SCIPdebugMessage("HessianSparsityDense for "); SCIPexprtreePrint(tree, NULL); printf("\n");
   SCIPdebugMessage("sparsity ="); for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) if (sparsity[i*n+j]) printf(" (%d,%d)", i, j); printf("\n");
#endif

   return SCIP_OKAY;
}

/** computes value and dense hessian of an expression tree
 * the full hessian is computed (lower left and upper right triangle)
 */
SCIP_RETCODE SCIPexprintHessianDense(
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree,               /** expression tree */
   SCIP_Real*            varvals,            /** values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /** have variable values changed since last call to an evaluation routine? */
   SCIP_Real*            val,                /** buffer to store function value */
   SCIP_Real*            hessian             /** buffer to store hessian values, need to have size at least n*n */
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

   //	int k = 0, l = 0;
   //	for (int i = 0; i < n; ++i)
   //	{ // row
   //		for (int j = 0; j <= i; ++j, ++k, ++l)
   //		{ // col
   //			assert(l < n*n);
   //			assert(k < n*(n+1)/2);
   //			hessian[k] = hess[l];
   //		}
   //		l += (n-i)-1;
   //	}
   //	assert(l == n*n);
   //	assert(k == n*(n+1)/2);

#ifdef SCIP_DEBUG
   SCIPdebugMessage("HessianDense for "); SCIPexprtreePrint(tree, NULL); printf("\n");
   SCIPdebugMessage("x    ="); for (int i = 0; i < n; ++i) printf("\t %g", data->x[i]); printf("\n");
   SCIPdebugMessage("hess ="); for (int i = 0; i < n*n; ++i) printf("\t %g", hessian[i]); printf("\n");
#endif

   return SCIP_OKAY;
}

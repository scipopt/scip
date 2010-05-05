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
#pragma ident "@(#) $Id: expression.c,v 1.5 2010/05/05 17:53:12 bzfviger Exp $"

/**@file   expression.c
 * @brief  methods for expressions and expression trees
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdarg.h>
#include <string.h>

#include "nlpi/expression.h"
#include "nlpi/struct_expression.h"
#include "nlpi/exprinterpret.h"

#define SCIP_EXPRESSION_MAXCHILDEST 20       /* estimate on maximal number of children */

/** sign of a value (-1 or +1)
 * 0.0 has sign +1
 */
#define SIGN(x) ((x) >= 0.0 ? 1.0 : -1.0)

/** signature of an expression (pointwise) evaluation function
 * The function should return nan, inf, or -inf in result if the function is undefined for the given arguments.
 *   
 * - opdata    operand data
 * - nargs     number of arguments
 * - argvals   values of arguments 
 * - varvals   values for variables 
 * - paramvals values for parameters 
 * - result    buffer where to store result of evaluation
 */
#define SCIP_DECL_EVAL(x) SCIP_RETCODE x (SCIP_EXPROPDATA  opdata, int nargs, SCIP_Real* argvals, SCIP_Real* varvals, SCIP_Real* paramvals, SCIP_Real* result)

/* element in table of expression operands */
struct SCIPexprOpTableElement
{
   const char*           name;               /**< name of operand (used for printing) */
   int                   nargs;              /**< number of arguments (negative if not fixed) */
   SCIP_DECL_EVAL        ((*eval));          /**< evaluation function */
};

static SCIP_DECL_EVAL( SCIPexprevalPushVar )
{
   assert(result  != NULL);
   assert(varvals != NULL);
   
   *result = varvals[opdata.intval];
   
   return SCIP_OKAY;
}

static SCIP_DECL_EVAL( SCIPexprevalPushValue )
{
   assert(result  != NULL);
   
   *result = opdata.dbl;
   
   return SCIP_OKAY;
}

static SCIP_DECL_EVAL( SCIPexprevalPushParameter )
{
   assert(result    != NULL);
   assert(paramvals != NULL );
   
   *result = paramvals[opdata.intval];
   
   return SCIP_OKAY;
}

static SCIP_DECL_EVAL( SCIPexprevalPlus )
{
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = argvals[0] + argvals[1];
   
   return SCIP_OKAY;
}

static SCIP_DECL_EVAL( SCIPexprevalMinus )
{
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = argvals[0] - argvals[1];
   
   return SCIP_OKAY;
}

static SCIP_DECL_EVAL( SCIPexprevalMult )
{
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = argvals[0] * argvals[1];

   return SCIP_OKAY;
}

static SCIP_DECL_EVAL( SCIPexprevalDiv )
{
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = argvals[0] / argvals[1];
   
   return SCIP_OKAY;
}

static SCIP_DECL_EVAL( SCIPexprevalSqr )
{
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = argvals[0] * argvals[0];
   
   return SCIP_OKAY;
}

static SCIP_DECL_EVAL( SCIPexprevalSqrt )
{
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = sqrt(argvals[0]);
   
   return SCIP_OKAY;
}

static SCIP_DECL_EVAL( SCIPexprevalPower )
{
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = pow(argvals[0], argvals[1]);

   return SCIP_OKAY;
}

static SCIP_DECL_EVAL( SCIPexprevalExp )
{
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = exp(argvals[0]);
   
   return SCIP_OKAY;
}

static SCIP_DECL_EVAL( SCIPexprevalLog )
{
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = log(argvals[0]);
   
   return SCIP_OKAY;
}

static SCIP_DECL_EVAL( SCIPexprevalSin )
{
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = sin(argvals[0]);

   return SCIP_OKAY;
}

static
SCIP_DECL_EVAL( SCIPexprevalCos )
{
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = cos(argvals[0]);
   
   return SCIP_OKAY;
}

static
SCIP_DECL_EVAL( SCIPexprevalTan )
{
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = tan(argvals[0]);
   
   return SCIP_OKAY;
}

static
SCIP_DECL_EVAL( SCIPexprevalErf )
{
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = erf(argvals[0]);
   
   return SCIP_OKAY;
}

static
SCIP_DECL_EVAL( SCIPexprevalErfi )
{
   assert(result  != NULL);
   assert(argvals != NULL);

   /* @TODO implement erfi evaluation */
   SCIPerrorMessage("erfi not implemented");
   
   return SCIP_ERROR;
}

static
SCIP_DECL_EVAL( SCIPexprevalMin )
{
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = MIN(argvals[0], argvals[1]);
   
   return SCIP_OKAY;
}

static
SCIP_DECL_EVAL( SCIPexprevalMax )
{
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = MAX(argvals[0], argvals[1]);

   return SCIP_OKAY;
}

static
SCIP_DECL_EVAL( SCIPexprevalAbs )
{
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = ABS(argvals[0]);
   
   return SCIP_OKAY;
}

static
SCIP_DECL_EVAL( SCIPexprevalSign )
{
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = SIGN(argvals[0]);
   
   return SCIP_OKAY;
}

static
SCIP_DECL_EVAL( SCIPexprevalSignPower )
{
   assert(result  != NULL);
   assert(argvals != NULL);

   if( argvals[0] > 0 )
     *result =  pow( argvals[0], argvals[1]);
   else
     *result = -pow(-argvals[0], argvals[1]);

   return SCIP_OKAY;
}

static
SCIP_DECL_EVAL( SCIPexprevalIntPower )
{
   int       n;

   assert(result  != NULL);
   assert(argvals != NULL);

   n = opdata.intval;

   if( n == 0 )
   {
      *result = 1.0;
      return SCIP_OKAY;
   }

   if( n > 0 )
   {
      if( n == 1 )
      {
         *result = 1.0;
         return SCIP_OKAY;
      }

      *result = SIGN(argvals[0]) * pow(ABS(argvals[0]), n);
      return SCIP_OKAY;
   }

   if( n == -1 )
   {
      *result = 1.0 / argvals[0];
      return SCIP_OKAY;
   }
   *result = SIGN(argvals[0]) / pow(ABS(argvals[0]), -n);

   return SCIP_OKAY;
}

static
SCIP_DECL_EVAL( SCIPexprevalSum )
{
   int i;
   
   assert(result  != NULL);
   assert(argvals != NULL);
   
   *result = 0.0;
   for( i = 0; i < nargs; ++i )
      *result += argvals[i];
   
   return SCIP_OKAY;
}

static 
SCIP_DECL_EVAL( SCIPexprevalProduct )
{
   int i;
   
   assert(result  != NULL);
   assert(argvals != NULL);
   
   *result = 1.0;
   for( i = 0; i < nargs; ++i )
      *result *= argvals[i];
   
   return SCIP_OKAY;
}

#define SCIPexprevalLinear         NULL
#define SCIPexprevalSignomial      NULL
#define SCIPexprevalQuadratic      NULL
#define SCIPexprevalPolynom        NULL

/** table containing for each operand the name, the number of children, and some evaluation functions */
struct SCIPexprOpTableElement SCIPexprOpTable[] =
{
   {NULL,-1,NULL},
   { "variable",          0, SCIPexprevalPushVar       },
   { "constant",          0, SCIPexprevalPushValue     },
   { "parameter",         0, SCIPexprevalPushParameter },
   {NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL}, 
   { "plus",              2, SCIPexprevalPlus          },
   { "minus",             2, SCIPexprevalMinus         },
   { "mul",               2, SCIPexprevalMult          },
   { "div",               2, SCIPexprevalDiv           },
   { "sqr",               1, SCIPexprevalSqr           },
   { "sqrt",              1, SCIPexprevalSqrt          },
   { "power",             2, SCIPexprevalPower         },
   { "exp",               1, SCIPexprevalExp           },
   { "log",               1, SCIPexprevalLog           },
   { "sin",               1, SCIPexprevalSin           },
   { "cos",               1, SCIPexprevalCos           },
   { "tan",               1, SCIPexprevalTan           },
   { "erf",               1, SCIPexprevalErf           },
   { "erfi",              1, SCIPexprevalErfi          },
   { "min",               2, SCIPexprevalMin           },
   { "max",               2, SCIPexprevalMax           },
   { "abs",               1, SCIPexprevalAbs           },
   { "sign",              1, SCIPexprevalSign          },
   { "signpower",         2, SCIPexprevalSignPower     },
   { "intpower",          1, SCIPexprevalIntPower      },
   {NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},
   {NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},
   {NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},
   {NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},
   {NULL,-1,NULL},{NULL,-1,NULL},{NULL,-1,NULL},
   { "sum",              -2, SCIPexprevalSum           },
   { "prod",             -2, SCIPexprevalProduct       },
   { "linear",           -2, SCIPexprevalLinear        },
   { "signomial",        -2, SCIPexprevalSignomial     },
   { "quadratic",        -2, SCIPexprevalQuadratic     },
   { "polynom",          -2, SCIPexprevalPolynom       }
};

/** gives the name of an operand as string */
const char* SCIPexpropGetName(
   SCIP_EXPROP           op                  /**< expression operand */
)
{
   assert(op < SCIP_EXPR_LAST);

   return SCIPexprOpTable[op].name;
}

/** gives the number of children of a simple operand */
int SCIPexpropGetNChildren(
   SCIP_EXPROP           op                  /**< expression operand */
)
{
   assert(op < SCIP_EXPR_LAST);

   return SCIPexprOpTable[op].nargs;
}

/** creates an expression */
SCIP_RETCODE SCIPexprCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           expr,               /**< pointer to buffer for expression address */
   SCIP_EXPROP           op,                 /**< operand of expression */
   ...                                       /**< arguments of operand */
)
{
   va_list         ap;
   SCIP_EXPR**     children;
   SCIP_EXPROPDATA opdata;
   
   assert(blkmem != NULL);
   assert(expr   != NULL);

   switch( op )
   {
      case SCIP_EXPR_VARIDX:
      case SCIP_EXPR_PARAM:
         va_start( ap, op );
         opdata.intval = va_arg( ap, int );
         va_end( ap );
         
         assert( opdata.intval >= 0 );
         
         SCIP_CALL( SCIPexprCreateDirect( blkmem, expr, op, 0, NULL, opdata ) );
         break;
         
      case SCIP_EXPR_CONST:
         va_start(ap, op );
         opdata.dbl = va_arg( ap, SCIP_Real );
         va_end( ap );
         
         SCIP_CALL( SCIPexprCreateDirect( blkmem, expr, op, 0, NULL, opdata ) );
         break;

      /* operands with two children */
      case SCIP_EXPR_PLUS     :
      case SCIP_EXPR_MINUS    :
      case SCIP_EXPR_MUL      :
      case SCIP_EXPR_DIV      :
      case SCIP_EXPR_POWER    :
      case SCIP_EXPR_MIN      :
      case SCIP_EXPR_MAX      :
      case SCIP_EXPR_SIGNPOWER:
         if( BMSallocBlockMemoryArray(blkmem, &children, 2 ) == NULL )
            return SCIP_NOMEMORY;
         
         va_start(ap, op );
         children[0] = va_arg( ap, SCIP_EXPR* );
         children[1] = va_arg( ap, SCIP_EXPR* );
         va_end( ap );
         opdata.data = NULL; /* to avoid compiler warning about use of uninitialised value */
         
         SCIP_CALL( SCIPexprCreateDirect( blkmem, expr, op, 2, children, opdata ) );
         break;

      /* operands with one child */
      case SCIP_EXPR_SQUARE:
      case SCIP_EXPR_SQRT  :
      case SCIP_EXPR_EXP   :
      case SCIP_EXPR_LOG   :
      case SCIP_EXPR_SIN   :
      case SCIP_EXPR_COS   :
      case SCIP_EXPR_TAN   :
      case SCIP_EXPR_ABS   :
      case SCIP_EXPR_SIGN  :
         if( BMSallocBlockMemoryArray(blkmem, &children, 1 ) == NULL )
            return SCIP_NOMEMORY;
         
         va_start(ap, op );
         children[0] = va_arg( ap, SCIP_EXPR* );
         va_end( ap );
         opdata.data = NULL; /* to avoid compiler warning about use of uninitialised value */
         
         SCIP_CALL( SCIPexprCreateDirect( blkmem, expr, op, 1, children, opdata ) );
         break;
         
      case SCIP_EXPR_INTPOWER:
         if( BMSallocBlockMemoryArray(blkmem, &children, 1 ) == NULL )
            return SCIP_NOMEMORY;

         va_start(ap, op );
         children[0] = va_arg( ap, SCIP_EXPR* );
         opdata.intval = va_arg( ap, int);
         va_end( ap );

         SCIP_CALL( SCIPexprCreateDirect( blkmem, expr, op, 1, children, opdata ) );
         break;

      /* complex operands */
      case SCIP_EXPR_SUM      :
      case SCIP_EXPR_PRODUCT  :
      case SCIP_EXPR_LINEAR   :
      case SCIP_EXPR_SIGNOMIAL:
      case SCIP_EXPR_QUADRATIC:
      case SCIP_EXPR_POLYNOM  :
         return SCIP_ERROR;
         /* @TODO implement */
         
      default:
         SCIPerrorMessage("unknown operand: %d\n", op);
         return SCIP_ERROR;
   }
   
   return SCIP_OKAY;
}

/** creates an expression
 * Note, that the expression is allocated but for the children only the pointer is copied.
 */
SCIP_RETCODE SCIPexprCreateDirect(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           expr,               /**< pointer to buffer for expression address */
   SCIP_EXPROP           op,                 /**< operand of expression */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children,           /**< children */
   SCIP_EXPROPDATA       opdata              /**< operand data */
)
{
   assert(blkmem != NULL);
   assert(expr   != NULL);
   assert(children != NULL || nchildren == 0);
   assert(children == NULL || nchildren >  0);
   
   if( BMSallocBlockMemory(blkmem, expr) == NULL )
      return SCIP_NOMEMORY;
   
   (*expr)->op        = op;
   (*expr)->nchildren = nchildren;
   (*expr)->children  = children;
   (*expr)->data      = opdata;
   
   return SCIP_OKAY;
}

/** copies an expression including its children */
SCIP_RETCODE SCIPexprCopyDeep(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           targetexpr,         /**< buffer to store pointer to copied expression */
   SCIP_EXPR*            sourceexpr          /**< expression to copy */
)
{
   assert(blkmem     != NULL);
   assert(targetexpr != NULL);
   assert(sourceexpr != NULL);

   if( BMSduplicateBlockMemory(blkmem, targetexpr, sourceexpr) == NULL )
      return SCIP_NOMEMORY;
   
   if( sourceexpr->nchildren )
   {
      int i;
      
      /* alloc memory for children expressions */
      if( BMSallocBlockMemoryArray(blkmem, &(*targetexpr)->children, sourceexpr->nchildren) == NULL )
         return SCIP_NOMEMORY;

      /* copy children expressions */
      for( i = 0; i < sourceexpr->nchildren; ++i )
      {
         SCIP_CALL( SCIPexprCopyDeep(blkmem, &(*targetexpr)->children[i], sourceexpr->children[i]) );
      }
   }
   else
   {
      assert((*targetexpr)->children == NULL); /* otherwise, sourceexpr->children was not NULL, which is wrong */
   }

   return SCIP_OKAY;
}

/** frees an expression including its children */
void SCIPexprFreeDeep(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           expr                /**< pointer to expression to free */
)
{
   assert(blkmem != NULL);
   assert(expr   != NULL);
   assert(*expr  != NULL);
   
   if( (*expr)->nchildren )
   {
      int i;
      
      assert( (*expr)->children != NULL );
      
      for( i = 0; i < (*expr)->nchildren; ++i )
      {
         SCIPexprFreeDeep(blkmem, &(*expr)->children[i]);
         assert((*expr)->children[i] == NULL);
      }

      BMSfreeBlockMemoryArray(blkmem, &(*expr)->children, (*expr)->nchildren);
   }
   else
   {
      assert( (*expr)->children == NULL );
   }
   
   BMSfreeBlockMemory(blkmem, expr);
}

/** gives operator of expression */
SCIP_EXPROP SCIPexprGetOperator(
   SCIP_EXPR*            expr                /**< expression */
)
{
   assert(expr != NULL);
   
   return expr->op;
}

/** gives number of children of an expression */
int SCIPexprGetNChildren(
   SCIP_EXPR*            expr                /**< expression */
)
{
   assert(expr != NULL);
   
   return expr->nchildren;
}

/** gives pointer to array with children of an expression */
SCIP_EXPR** SCIPexprGetChildren(
   SCIP_EXPR*            expr                /**< expression */
)
{
   assert(expr != NULL);
   
   return expr->children;
}

/** gives index belonging to a SCIP_EXPR_VARIDX or SCIP_EXPR_PARAM operand */
int SCIPexprGetOpIndex(
   SCIP_EXPR*            expr                /**< expression */
)
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_VARIDX || expr->op == SCIP_EXPR_PARAM);
   
   return expr->data.intval;
}

/** gives real belonging to a SCIP_EXPR_CONST operand */ 
SCIP_Real SCIPexprGetOpReal(
   SCIP_EXPR* expr                           /**< expression */
)
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_CONST);
   
   return expr->data.dbl;
}

/** gives void* belonging to a complex operand */
void* SCIPexprGetOpData(
   SCIP_EXPR*            expr                /**< expression */
)
{
   assert(expr != NULL);
   assert(expr->op >= 64); /* only complex operands store their data as void* */
   
   return expr->data.data;
}

/** gives exponent belonging to a SCIP_EXPR_INTPOWER operand */
int SCIPexprGetIntPowerExponent(
   SCIP_EXPR*            expr                /**< expression */
)
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_INTPOWER);

   return expr->data.intval;
}

/** indicates whether the expression contains a SCIP_EXPR_PARAM */
SCIP_Bool SCIPexprHasParam(
   SCIP_EXPR*            expr                /**< expression */
)
{
   int i;

   assert(expr != NULL);

   if( expr->op == SCIP_EXPR_PARAM )
      return TRUE;

   for( i = 0; i < expr->nchildren; ++i )
      if( SCIPexprHasParam(expr->children[i]) )
         return TRUE;

   return FALSE;
}

/** gets maximal degree of expression, or 65535 if not a polynom */
SCIP_RETCODE SCIPexprGetMaxDegree(
   SCIP_EXPR*            expr,               /**< expression */
   int*                  maxdegree           /**< buffer to store maximal degree */
)
{
   int child1;
   int child2;

   assert(expr      != NULL);
   assert(maxdegree != NULL);

   switch (expr->op)
   {
      case SCIP_EXPR_VARIDX:
         *maxdegree = 1;
         break;

      case SCIP_EXPR_CONST:         
      case SCIP_EXPR_PARAM:
         *maxdegree = 0;
         break;

      case SCIP_EXPR_PLUS:
      case SCIP_EXPR_MINUS:
      {
         assert(expr->children[0] != NULL);
         assert(expr->children[1] != NULL);
         
         SCIP_CALL( SCIPexprGetMaxDegree(expr->children[0], &child1) );
         SCIP_CALL( SCIPexprGetMaxDegree(expr->children[1], &child2) );
         
         *maxdegree = MAX(child1, child2);
         break;
      }

      case SCIP_EXPR_MUL:
      {
         assert(expr->children[0] != NULL);
         assert(expr->children[1] != NULL);
         
         SCIP_CALL( SCIPexprGetMaxDegree(expr->children[0], &child1) );
         SCIP_CALL( SCIPexprGetMaxDegree(expr->children[1], &child2) );
         
         *maxdegree = child1 + child2;
         break;
      }

      case SCIP_EXPR_DIV:
      {
         assert(expr->children[0] != NULL);
         assert(expr->children[1] != NULL);
         
         SCIP_CALL( SCIPexprGetMaxDegree(expr->children[0], &child1) );
         SCIP_CALL( SCIPexprGetMaxDegree(expr->children[1], &child2) );
         
         /* if not division by constant, then it is not a polynom */
         *maxdegree = (child2 != 0) ? 65535 : child1;
         break;
      }

      case SCIP_EXPR_SQUARE:
      {
         assert(expr->children[0] != NULL);
         
         SCIP_CALL( SCIPexprGetMaxDegree(expr->children[0], &child1) );
         
         *maxdegree = 2 * child1;
         break;
      }

      case SCIP_EXPR_SQRT:
      {
         assert(expr->children[0] != NULL);
         
         SCIP_CALL( SCIPexprGetMaxDegree(expr->children[0], &child1) );
         
         /* if not squareroot of constant, then no polynomial */
         *maxdegree = (child1 != 0) ? 65535 : 0;
         break;
      }

      case SCIP_EXPR_POWER:
      {
         SCIP_Real val;
         
         assert(expr->children[0] != NULL);
         assert(expr->children[1] != NULL);

         SCIP_CALL( SCIPexprGetMaxDegree(expr->children[0], &child1) );
         SCIP_CALL( SCIPexprGetMaxDegree(expr->children[1], &child2) );

         /* constant ^ constant has degree 0 */
         if( child1 == 0 && child2 == 0 )
         {
            *maxdegree = 0;
            break;
         }
         
         /* non-polynomial ^ non-constant is not a polynom */
         if( child1 >= 65535 || child2 > 0 || SCIPexprHasParam(expr->children[1]) )
         {
            *maxdegree = 65535;
            break;
         }

         /* so it is polynomial ^ constant
          * let's see whether the constant is integral */
         SCIP_CALL( SCIPexprEval(expr->children[1], NULL, NULL, &val) );

         if( val == 0.0 ) /* polynomial ^ 0 == 0 */
            *maxdegree = 0;
         else if( val > 0.0 && floor(val) == val ) /* natural exponent gives polynom again */ 
            *maxdegree = child1 * (int)floor(val);
         else /* negative or nonintegral exponent does not give polynom */ 
            *maxdegree = 65535;

         break;
      }

      case SCIP_EXPR_EXP:
      case SCIP_EXPR_LOG:
      case SCIP_EXPR_SIN:
      case SCIP_EXPR_COS:
      case SCIP_EXPR_TAN:
      case SCIP_EXPR_ERF:
      case SCIP_EXPR_ERFI:
      case SCIP_EXPR_ABS:
      case SCIP_EXPR_SIGN:
      {
         assert(expr->children[0] != NULL);
         
         SCIP_CALL( SCIPexprGetMaxDegree(expr->children[0], &child1) );
         
         /* if argument is not a constant, then no polynom, otherwise it is a constant */
         *maxdegree = (child1 != 0) ? 65535 : 0;
         break;
      }

      case SCIP_EXPR_MIN:
      case SCIP_EXPR_MAX:
      case SCIP_EXPR_SIGNPOWER:
      {
         assert(expr->children[0] != NULL);
         assert(expr->children[1] != NULL);

         SCIP_CALL( SCIPexprGetMaxDegree(expr->children[0], &child1) );
         SCIP_CALL( SCIPexprGetMaxDegree(expr->children[1], &child2) );

         /* if any of the operands is not constant, then it is no polynom */
         *maxdegree = (child1 != 0 || child2 != 0) ? 65535 : 0;
         break;
      }

      case SCIP_EXPR_INTPOWER:
      {
         assert(expr->children[0] != NULL);

         SCIP_CALL( SCIPexprGetMaxDegree(expr->children[0], &child1) );

         /* constant ^ integer or something ^ 0 has degree 0 */
         if( child1 == 0 || expr->data.intval == 0 )
         {
            *maxdegree = 0;
            break;
         }

         /* non-polynomial ^ integer  or  something ^ negative  is not a polynom */
         if( child1 >= 65535 || expr->data.intval < 0 )
         {
            *maxdegree = 65535;
            break;
         }

         /* so it is polynomial ^ natural, which gives a polynom again */
         *maxdegree = child1 * expr->data.intval;

         break;
      }

      default:
         return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** evaluates an expression */
SCIP_RETCODE SCIPexprEval(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real*            varvals,            /**< values for variables, can be NULL if the expression is constant */
   SCIP_Real*            param,              /**< values for parameters, can be NULL if the expression is not parameterized */
   SCIP_Real*            val                 /**< buffer to store value */
)
{
   int i;
   SCIP_Real  staticbuf[SCIP_EXPRESSION_MAXCHILDEST];
   SCIP_Real* buf;

   /* if many children, get large enough memory to store argument values */
   if( expr->nchildren > SCIP_EXPRESSION_MAXCHILDEST )
   {
      if( BMSallocMemoryArray(&buf, expr->nchildren) == NULL )
         return SCIP_NOMEMORY;
   }
   else
   {
      buf = staticbuf;
   }

   /* evaluate children */
   for( i = 0; i < expr->nchildren; ++i )
   {
      SCIP_CALL( SCIPexprEval(expr->children[i], varvals, param,  &buf[i]) );
   }

   /* evaluate this expression */
   assert( SCIPexprOpTable[expr->op].eval != NULL );
   SCIP_CALL( SCIPexprOpTable[expr->op].eval(expr->data, expr->nchildren, buf, varvals, param, val) );

   /* free memory, if allocated before */
   if( staticbuf != buf )
   {
      BMSfreeMemoryArray(&buf);
   }

   return SCIP_OKAY;
}

/** prints an expression */
void SCIPexprPrint(
   SCIP_EXPR*            expr,               /**< expression */
   FILE*                 file,               /**< file for printing, or NULL for stdout */
   const char**          varnames,           /**< names of variables, or NULL for default names */
   const char**          paramnames          /**< names of parameters, or NULL for default names */
)
{
   assert( expr != NULL );

   switch( expr->op )
   {
      case SCIP_EXPR_VARIDX:
         if( varnames != NULL )
         {
            assert(varnames[expr->data.intval] != NULL);
            SCIPmessageFPrintInfo(file, "<%s>", varnames[expr->data.intval]);
         }
         else
         {
            SCIPmessageFPrintInfo(file, "<var%d>", expr->data.intval);
         }
         break;
         
      case SCIP_EXPR_PARAM:
         if( paramnames != NULL )
         {
            assert(paramnames[expr->data.intval] != NULL);
            SCIPmessageFPrintInfo(file, "<%s>", varnames[expr->data.intval]);
         }
         else
         {
            SCIPmessageFPrintInfo(file, "<param%d>", expr->data.intval );
         }
         break;
         
      case SCIP_EXPR_CONST:
         SCIPmessageFPrintInfo(file, "%lf", expr->data.dbl );
         break;

      case SCIP_EXPR_INTPOWER:
         SCIPmessageFPrintInfo(file, "%s(", SCIPexprOpTable[expr->op].name);
         SCIPexprPrint(expr->children[0], file, varnames, paramnames);
         SCIPmessageFPrintInfo(file, ", %d)", expr->data.intval);
         break;

      default:
      {
         int i;
         
         SCIPmessageFPrintInfo(file, "%s(", SCIPexprOpTable[expr->op].name);
         
         for( i = 0; i < expr->nchildren; ++i )
         {
            SCIPexprPrint(expr->children[i], file, varnames, paramnames);
            if( i + 1 < expr->nchildren )
            {
               SCIPmessageFPrintInfo(file, ", ");
            }
         }

         SCIPmessageFPrintInfo(file, ")");
      }
   }
}

/** creates an expression tree */
SCIP_RETCODE SCIPexprtreeCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPRTREE**       tree,               /**< buffer to store address of created expression tree */
   SCIP_EXPR*            root,               /**< pointer to root expression, not copied deep !, can be NULL */
   int                   nvars,              /**< number of variables in variable mapping */
   int                   nparams,            /**< number of parameters in expression */
   SCIP_Real*            params              /**< values for parameters, or NULL (if NULL but nparams > 0, then params is initialized with zeros) */
)
{
   assert(blkmem != NULL);
   assert(tree   != NULL);

   if( BMSallocBlockMemory(blkmem, tree) == NULL )
      return SCIP_NOMEMORY;

   (*tree)->blkmem    = blkmem;
   (*tree)->root      = root;
   (*tree)->nvars     = nvars;
   (*tree)->vars      = NULL;
   (*tree)->nparams   = nparams;
   (*tree)->interpreterdata = NULL;
   
   if( params != NULL )
   {
      assert(nparams > 0);
      if( BMSduplicateBlockMemoryArray(blkmem, &(*tree)->params, params, nparams) == NULL )
         return SCIP_NOMEMORY;
   }
   else if( nparams > 0 )
   {
      if( BMSallocBlockMemoryArray(blkmem, &(*tree)->params, nparams) == NULL )
         return SCIP_NOMEMORY;
      BMSclearMemoryArray((*tree)->params, nparams);
   }
   else
   {
      assert(nparams == 0);
      (*tree)->params = NULL;
   }

   return SCIP_OKAY;
}

/** copies an expression tree */
SCIP_RETCODE SCIPexprtreeCopy(
   BMS_BLKMEM*           blkmem,             /**< block memory that should be used in new expression tree */
   SCIP_EXPRTREE**       targettree,         /**< buffer to store address of copied expression tree */
   SCIP_EXPRTREE*        sourcetree          /**< expression tree to copy */
)
{
   assert(blkmem     != NULL);
   assert(targettree != NULL);
   assert(sourcetree != NULL);

   /* copy expression tree "header" */
   if( BMSduplicateBlockMemory(blkmem, targettree, sourcetree) == NULL )
      return SCIP_NOMEMORY;
   
   /* we may have a new block memory; and we do not want to keep the others interpreter data */
   (*targettree)->blkmem          = blkmem;
   (*targettree)->interpreterdata = NULL;
   
   /* copy variables, if any */
   if( sourcetree->vars != NULL )
   {
      assert(sourcetree->nvars > 0);
      
      if( BMSduplicateBlockMemoryArray(blkmem, &(*targettree)->vars, sourcetree->vars, sourcetree->nvars) == NULL )
         return SCIP_NOMEMORY;
   }
   
   /* copy parameters, if any */
   if( sourcetree->params != NULL )
   {
      assert(sourcetree->nparams > 0);
      
      if( BMSduplicateBlockMemoryArray(blkmem, &(*targettree)->params, sourcetree->params, sourcetree->nparams) == NULL )
         return SCIP_NOMEMORY;
   }

   /* copy expression */
   SCIP_CALL( SCIPexprCopyDeep(blkmem, &(*targettree)->root, sourcetree->root) );
   
   return SCIP_OKAY;
}

/** frees an expression tree */
SCIP_RETCODE SCIPexprtreeFree(
   SCIP_EXPRTREE**       tree                /**< pointer to expression tree that is freed */
)
{
   assert( tree != NULL);
   assert(*tree != NULL);
   
   if( (*tree)->interpreterdata )
   {
      SCIP_CALL( SCIPexprintFreeData(&(*tree)->interpreterdata) );
      assert((*tree)->interpreterdata == NULL);
   }
   
   SCIPexprFreeDeep((*tree)->blkmem, &(*tree)->root);
   assert((*tree)->root == NULL);
   
   BMSfreeBlockMemoryArrayNull((*tree)->blkmem, &(*tree)->vars,   (*tree)->nvars  );
   BMSfreeBlockMemoryArrayNull((*tree)->blkmem, &(*tree)->params, (*tree)->nparams);

   BMSfreeBlockMemory((*tree)->blkmem, tree);
   
   return SCIP_OKAY;
}

/** returns root expression of an expression tree */
SCIP_EXPR* SCIPexprtreeGetRoot(
   SCIP_EXPRTREE*        tree                /**< expression tree */
)
{
   assert(tree != NULL);
   
   return tree->root;
}

/** returns number of variables in expression tree */
int SCIPexprtreeGetNVars(
   SCIP_EXPRTREE*        tree                /**< expression tree */
)
{
   assert(tree != NULL);
   
   return tree->nvars;
}

/** returns values of parameters or NULL if none */
SCIP_Real* SCIPexprtreeGetParamVals(
   SCIP_EXPRTREE*        tree                /**< expression tree */
)
{
   assert(tree != NULL);
   
   return tree->params;
}

/** gets data of expression tree interpreter
 * @return NULL if not set
 */
SCIP_EXPRINTDATA* SCIPexprtreeGetInterpreterData(
   SCIP_EXPRTREE*        tree                /**< expression tree */
)
{
   assert(tree != NULL);
   
   return tree->interpreterdata;
}

/** indicates whether there are parameterized constants (SCIP_EXPR_PARAM) in expression tree */
SCIP_Bool SCIPexprtreeHasParam(
   SCIP_EXPRTREE*        tree                /**< expression tree */
)
{
   assert(tree != NULL);

   return SCIPexprHasParam(tree->root);
}

/** Gives maximal degree of expression in expression tree.
 * If constant expression, gives 0,
 * if linear expression, gives 1,
 * if polynomial expression, gives its maximal degree,
 * otherwise (nonpolynomial nonconstant expressions) gives at least 65535.
 */
SCIP_RETCODE SCIPexprtreeGetMaxDegree(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   int*                  maxdegree           /**< buffer to store maximal degree */
)
{
   assert(tree != NULL);
   
   SCIP_CALL( SCIPexprGetMaxDegree(tree->root, maxdegree) );

   return SCIP_OKAY;
}

/** sets data of expression tree interpreter */
void SCIPexprtreeSetInterpreterData(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_EXPRINTDATA*     interpreterdata     /**< expression interpreter data */
)
{
   assert(tree != NULL);
   assert(interpreterdata != NULL);
   assert(tree->interpreterdata == NULL);

   tree->interpreterdata = interpreterdata;
}

/** evaluates an expression tree */
SCIP_RETCODE SCIPexprtreeEval(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values for variables */
   SCIP_Real*            val                 /**< buffer to store expression tree value */
)
{
   assert(tree    != NULL);
   assert(varvals != NULL || tree->nvars == 0);
   assert(val     != NULL);

   SCIP_CALL( SCIPexprEval(tree->root, varvals, tree->params, val) );
   
   return SCIP_OKAY;
}

/** prints an expression tree */
void SCIPexprtreePrint(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   FILE*                 file,               /**< file for printing, or NULL for stdout */
   const char**          varnames,           /**< names of variables, or NULL for default names */
   const char**          paramnames          /**< names of parameters, or NULL for default names */
)
{
   assert(tree != NULL);

   SCIPexprPrint(tree->root, file, varnames, paramnames);
}

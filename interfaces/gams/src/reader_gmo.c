/* Copyright (C) GAMS Development and others 2011
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Author: Stefan Vigerske
 */

/**@file   reader_gmo.c
 * @ingroup FILEREADERS
 * @brief  GMO file reader
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <math.h>
#include <string.h>

/* dos compiler does not know PI */
#ifndef M_PI
#define M_PI           3.141592653589793238462643
#endif

#include "reader_gmo.h"

#include "gmomcc.h"
#include "gevmcc.h"
#include "gdxcc.h"
#include "optcc.h"
#include "GamsCompatibility.h"
#include "GamsNLinstr.h"

#include "scip/cons_linear.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_indicator.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "scip/heur_subnlp.h"
#include "scip/dialog_default.h"
#include "nlpi/struct_expr.h"

#define READER_NAME             "gmoreader"
#define READER_DESC             "Gams Control file reader (using GMO API)"
#define READER_EXTENSION        "dat"


/*
 * Data structures
 */

/** data for gmo reader */
struct SCIP_ReaderData
{
   gmoHandle_t           gmo;                /**< GAMS model object */
   gevHandle_t           gev;                /**< GAMS environment */
   int                   mipstart;           /**< how to handle initial variable levels */
   char*                 indicatorfile;      /**< name of GAMS options file that contains definitions on indicators */
};

/** problem data */
struct SCIP_ProbData
{
   int                   nvars;              /**< number of variables in vars array */
   SCIP_VAR**            vars;               /**< SCIP variables as corresponding to GMO variables */
   SCIP_VAR*             objvar;             /**< SCIP variable used to model objective function */
   SCIP_VAR*             objconst;           /**< SCIP variable used to model objective constant */
};

/*
 * Callback methods of probdata
 */

/** frees user data of original problem (called when the original problem is freed)
 */
static
SCIP_DECL_PROBDELORIG(probdataDelOrigGmo)
{
   int i;

   assert((*probdata)->vars != NULL || (*probdata)->nvars > 0);

   for( i = 0; i < (*probdata)->nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->vars[i]) );
   }
   SCIPfreeMemoryArray(scip, &(*probdata)->vars);

   if( (*probdata)->objvar != NULL )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->objvar) );
   }

   if( (*probdata)->objconst != NULL )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->objconst) );
   }

   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed)
 */
#if 0
static
SCIP_DECL_PROBTRANS(probdataTransGmo)
{

}
#else
#define probdataTransGmo NULL
#endif

/** frees user data of transformed problem (called when the transformed problem is freed)
 */
#if 0
static
SCIP_DECL_PROBDELTRANS(probdataDelTransGmo)
{
   return SCIP_OKAY;
}
#else
#define probdataDelTransGmo NULL
#endif


/** solving process initialization method of transformed data (called before the branch and bound process begins)
 */
#if 0
static
SCIP_DECL_PROBINITSOL(probdataInitSolGmo)
{

}
#else
#define probdataInitSolGmo NULL
#endif

/** solving process deinitialization method of transformed data (called before the branch and bound data is freed)
 */
#if 0
static
SCIP_DECL_PROBEXITSOL(probdataExitSolGmo)
{
   return SCIP_OKAY;
}
#else
#define probdataExitSolGmo NULL
#endif

/** copies user data of source SCIP for the target SCIP
 */
#if 0
static
SCIP_DECL_PROBCOPY(probdataCopyGmo)
{

}
#else
#define probdataCopyGmo NULL
#endif

/*
 * Local methods of reader
 */


/** ensures that an array of variables has at least a given length */
static
SCIP_RETCODE ensureVarsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to variables array */
   int*                  varssize,           /**< pointer to length of variables array */
   int                   nvars               /**< desired minimal length of array */
   )
{
   assert(scip != NULL);
   assert(vars != NULL);
   assert(varssize != NULL);

   if( nvars < *varssize )
      return SCIP_OKAY;

   *varssize = SCIPcalcMemGrowSize(scip, nvars);
   SCIP_CALL( SCIPreallocBufferArray(scip, vars, *varssize) );
   assert(nvars <= *varssize);

   return SCIP_OKAY;
}

/** adds new children to a linear expression */
static
SCIP_RETCODE exprLinearAdd(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< pointer to store resulting expression */
   int                   nchildren,          /**< number of children to add */
   SCIP_Real*            coefs,              /**< coefficients of additional children */
   SCIP_EXPR**           children,           /**< additional children expressions */
   SCIP_Real             constant            /**< constant to add */
)
{
   SCIP_Real* data;

   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetOperator(expr) == SCIP_EXPR_LINEAR);
   assert(nchildren >= 0);
   assert(coefs != NULL || nchildren == 0);
   assert(children != NULL || nchildren == 0);

   data = (SCIP_Real*)SCIPexprGetOpData(expr);
   assert(data != NULL);

   /* handle simple case of adding a constant */
   if( nchildren == 0 )
   {
      data[SCIPexprGetNChildren(expr)] += constant;

      return SCIP_OKAY;
   }

   /* add new children to expr's children array */
   SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &expr->children, expr->nchildren, expr->nchildren + nchildren) );
   BMScopyMemoryArray(&expr->children[expr->nchildren], children, nchildren);  /*lint !e866*/

   /* add constant and new coefs to expr's data array */
   SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &data, expr->nchildren + 1, expr->nchildren + nchildren + 1) );
   data[expr->nchildren + nchildren] = data[expr->nchildren] + constant;
   BMScopyMemoryArray(&data[expr->nchildren], coefs, nchildren); /*lint !e866*/
   expr->data.data = (void*)data;

   expr->nchildren += nchildren;

   return SCIP_OKAY;
}

/** frees a linear expression, but not its children */
static
void exprLinearFree(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR**           expr                /**< linear expression to free */
   )
{
   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetOperator(*expr) == SCIP_EXPR_LINEAR);

   BMSfreeBlockMemoryArray(blkmem, (SCIP_Real**)&(*expr)->data.data, (*expr)->nchildren + 1); /*lint !e866*/
   BMSfreeBlockMemoryArray(blkmem, &(*expr)->children, (*expr)->nchildren);

   BMSfreeBlockMemory(blkmem, expr);
}

/** creates an expression from the addition of two given expression, with coefficients, and a constant */
static
SCIP_RETCODE exprAdd(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR**           expr,               /**< pointer to store resulting expression */
   SCIP_Real             coef1,              /**< coefficient of first term */
   SCIP_EXPR*            term1,              /**< expression of first term */
   SCIP_Real             coef2,              /**< coefficient of second term */
   SCIP_EXPR*            term2,              /**< expression of second term */
   SCIP_Real             constant            /**< constant term to add */
)
{
   assert(blkmem != NULL);
   assert(expr != NULL);

   if( term1 != NULL && SCIPexprGetOperator(term1) == SCIP_EXPR_CONST )
   {
      constant += coef1 * SCIPexprGetOpReal(term1);
      SCIPexprFreeDeep(blkmem, &term1);
   }

   if( term2 != NULL && SCIPexprGetOperator(term2) == SCIP_EXPR_CONST )
   {
      constant += coef2 * SCIPexprGetOpReal(term2);
      SCIPexprFreeDeep(blkmem, &term2);
   }

   if( term1 == NULL && term2 == NULL )
   {
      SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_CONST, constant) );
      return SCIP_OKAY;
   }

   if( term1 != NULL && SCIPexprGetOperator(term1) == SCIP_EXPR_LINEAR && coef1 != 1.0 )
   {
      /* multiply coefficients and constant of linear expression term1 by coef1 */
      SCIP_Real* coefs;
      int i;

      coefs = SCIPexprGetLinearCoefs(term1);
      assert(coefs != NULL);

      for( i = 0; i < SCIPexprGetNChildren(term1); ++i )
         coefs[i] *= coef1;

      SCIP_CALL( exprLinearAdd(blkmem, term1, 0, NULL, NULL, (coef1-1.0) * SCIPexprGetLinearConstant(term1)) );

      coef1 = 1.0;
   }

   if( term2 != NULL && SCIPexprGetOperator(term2) == SCIP_EXPR_LINEAR && coef2 != 1.0 )
   {
      /* multiply coefficients and constant of linear expression term2 by coef2 */
      SCIP_Real* coefs;
      int i;

      coefs = SCIPexprGetLinearCoefs(term2);
      assert(coefs != NULL);

      for( i = 0; i < SCIPexprGetNChildren(term2); ++i )
         coefs[i] *= coef2;

      SCIP_CALL( exprLinearAdd(blkmem, term2, 0, NULL, NULL, (coef2-1.0) * SCIPexprGetLinearConstant(term2)) );

      coef2 = 1.0;
   }

   if( term1 == NULL || term2 == NULL )
   {
      if( term1 == NULL )
      {
         term1 = term2;
         coef1 = coef2;
         /* term2 = NULL; */
      }
      if( constant != 0.0 || coef1 != 1.0 )
      {
         if( SCIPexprGetOperator(term1) == SCIP_EXPR_LINEAR )
         {
            assert(coef1 == 1.0);

            /* add constant to existing linear expression */
            SCIP_CALL( exprLinearAdd(blkmem, term1, 0, NULL, NULL, constant) );
            *expr = term1;
         }
         else
         {
            /* create new linear expression for coef1 * term1 + constant */
            SCIP_CALL( SCIPexprCreateLinear(blkmem, expr, 1, &term1, &coef1, constant) );
         }
      }
      else
      {
         assert(constant == 0.0);
         assert(coef1 == 1.0);
         *expr = term1;
      }

      return SCIP_OKAY;
   }

   if( SCIPexprGetOperator(term1) == SCIP_EXPR_LINEAR && SCIPexprGetOperator(term2) == SCIP_EXPR_LINEAR )
   {
      assert(coef1 == 1.0);
      assert(coef2 == 1.0);

      SCIP_CALL( exprLinearAdd(blkmem, term1, SCIPexprGetNChildren(term2), SCIPexprGetLinearCoefs(term2), SCIPexprGetChildren(term2), SCIPexprGetLinearConstant(term2) + constant) );
      exprLinearFree(blkmem, &term2);

      *expr = term1;

      return SCIP_OKAY;
   }

   if( SCIPexprGetOperator(term2) == SCIP_EXPR_LINEAR )
   {
      /* if only term2 is linear, then swap */
      SCIP_EXPR* tmp;

      tmp = term2;
      assert(coef2 == 1.0);

      term2 = term1;
      coef2 = coef1;
      term1 = tmp;
      coef1 = 1.0;
   }

   if( SCIPexprGetOperator(term1) == SCIP_EXPR_LINEAR )
   {
      /* add coef2*term2 as extra child to linear expression term1 */
      assert(coef1 == 1.0);

      SCIP_CALL( exprLinearAdd(blkmem, term1, 1, &coef2, &term2, constant) );
      *expr = term1;

      return SCIP_OKAY;
   }

   /* both terms are not linear, then create new linear term for sum */
   {
      SCIP_Real coefs[2];
      SCIP_EXPR* children[2];

      coefs[0] = coef1;
      coefs[1] = coef2;
      children[0] = term1;
      children[1] = term2;

      SCIP_CALL( SCIPexprCreateLinear(blkmem, expr, 2, children, coefs, constant) );
   }

   return SCIP_OKAY;
}

/** creates an expression tree from given GAMS nonlinear instructions */
static
SCIP_RETCODE makeExprtree(
   SCIP*                 scip,               /**< SCIP data structure */
   gmoHandle_t           gmo,                /**< GAMS Model Object */
   int                   codelen,
   int*                  opcodes,
   int*                  fields,
   SCIP_Real*            constants,
   SCIP_EXPRTREE**       exprtree            /**< buffer where to store expression tree */
)
{
   SCIP_PROBDATA* probdata;
   BMS_BLKMEM*   blkmem;
   SCIP_HASHMAP* var2idx;
   SCIP_EXPR**   stack;
   int           stackpos;
   int           stacksize;
   SCIP_VAR**    vars;
   int           nvars;
   int           varssize;
   int           pos;
   SCIP_EXPR*    e;
   SCIP_EXPR*    term1;
   SCIP_EXPR*    term2;
   GamsOpCode    opcode;
   int           address;
   int           varidx;
   int           nargs;

   assert(scip != NULL);
   assert(opcodes != NULL);
   assert(fields != NULL);
   assert(constants != NULL);
   assert(exprtree != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);
   assert(probdata->vars != NULL);

   blkmem = SCIPblkmem(scip);

   stackpos = 0;
   stacksize = 20;
   SCIP_CALL( SCIPallocBufferArray(scip, &stack, stacksize) );

   nvars = 0;
   varssize = 10;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, varssize) );
   SCIP_CALL( SCIPhashmapCreate(&var2idx, blkmem, SCIPgetNVars(scip)) );

   nargs = -1;

   for( pos = 0; pos < codelen; ++pos )
   {
      opcode = (GamsOpCode)opcodes[pos];
      address = fields[pos]-1;

      SCIPdebugMessage("%s: ", GamsOpCodeName[opcode]);

      e = NULL;

      switch( opcode )
      {
         case nlNoOp: /* no operation */
         case nlStore: /* store row */
         case nlHeader:
         {
            SCIPdebugPrintf("ignored\n");
            break;
         }

         case nlPushV: /* push variable */
         {
            address = gmoGetjSolver(gmo, address);
            SCIPdebugPrintf("push variable %d = <%s>\n", address, SCIPvarGetName(probdata->vars[address]));

            if( !SCIPhashmapExists(var2idx, probdata->vars[address]) )
            {
               /* add variable to list of variables */
               SCIP_CALL( ensureVarsSize(scip, &vars, &varssize, nvars+1) );
               assert(nvars < varssize);
               vars[nvars] = probdata->vars[address];
               varidx = nvars;
               ++nvars;
               SCIP_CALL( SCIPhashmapInsert(var2idx, (void*)vars[varidx], (void*)(size_t)varidx) );
            }
            else
            {
               varidx = (int)(size_t)SCIPhashmapGetImage(var2idx, (void*)probdata->vars[address]);
               assert(varidx >= 0);
               assert(varidx < nvars);
               assert(vars[varidx] == probdata->vars[address]);
            }
            SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_VARIDX, varidx) );
            break;
         }

         case nlPushI: /* push constant */
         {
            SCIPdebugPrintf("push constant %g\n", constants[address]);
            SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_CONST, constants[address]) );
            break;
         }

         case nlPushZero: /* push zero */
         {
            SCIPdebugPrintf("push constant zero\n");

            SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_CONST, 0.0) );
            break;
         }

         case nlAdd : /* add */
         {
            SCIPdebugPrintf("add\n");
            assert(stackpos >= 2);
            term1 = stack[stackpos-1];
            --stackpos;
            term2 = stack[stackpos-1];
            --stackpos;

            SCIP_CALL( exprAdd(blkmem, &e, 1.0, term1, 1.0, term2, 0.0) );

            break;
         }

         case nlAddV: /* add variable */
         {
            address = gmoGetjSolver(gmo, address);
            SCIPdebugPrintf("add variable %d = <%s>\n", address, SCIPvarGetName(probdata->vars[address]));

            assert(stackpos >= 1);
            term1 = stack[stackpos-1];
            --stackpos;

            if( !SCIPhashmapExists(var2idx, probdata->vars[address]) )
            {
               /* add variable to list of variables */
               SCIP_CALL( ensureVarsSize(scip, &vars, &varssize, nvars+1) );
               assert(nvars < varssize);
               vars[nvars] = probdata->vars[address];
               varidx = nvars;
               ++nvars;
               SCIP_CALL( SCIPhashmapInsert(var2idx, (void*)vars[varidx], (void*)(size_t)varidx) );
            }
            else
            {
               varidx = (int)(size_t)SCIPhashmapGetImage(var2idx, (void*)probdata->vars[address]);
               assert(varidx >= 0);
               assert(varidx < nvars);
               assert(vars[varidx] == probdata->vars[address]);
            }

            SCIP_CALL( SCIPexprCreate(blkmem, &term2, SCIP_EXPR_VARIDX, varidx) );
            SCIP_CALL( exprAdd(blkmem, &e, 1.0, term1, 1.0, term2, 0.0) );

            break;
         }

         case nlAddI: /* add immediate */
         {
            SCIPdebugPrintf("add constant %g\n", constants[address]);

            assert(stackpos >= 1);
            term1 = stack[stackpos-1];
            --stackpos;

            SCIP_CALL( exprAdd(blkmem, &e, 1.0, term1, 1.0, NULL, constants[address]) );

            break;
         }

         case nlSub: /* substract */
         {
            SCIPdebugPrintf("substract\n");

            assert(stackpos >= 2);
            term1 = stack[stackpos-1];
            --stackpos;
            term2 = stack[stackpos-1];
            --stackpos;

            SCIP_CALL( exprAdd(blkmem, &e, 1.0, term2, -1.0, term1, 0.0) );

            break;
         }

         case nlSubV: /* subtract variable */
         {
            address = gmoGetjSolver(gmo, address);
            SCIPdebugPrintf("subtract variable %d = <%s>\n", address, SCIPvarGetName(probdata->vars[address]));

            assert(stackpos >= 1);
            term1 = stack[stackpos-1];
            --stackpos;

            if( !SCIPhashmapExists(var2idx, probdata->vars[address]) )
            {
               /* add variable to list of variables */
               SCIP_CALL( ensureVarsSize(scip, &vars, &varssize, nvars+1) );
               assert(nvars < varssize);
               vars[nvars] = probdata->vars[address];
               varidx = nvars;
               ++nvars;
               SCIP_CALL( SCIPhashmapInsert(var2idx, (void*)vars[varidx], (void*)(size_t)varidx) );
            }
            else
            {
               varidx = (int)(size_t)SCIPhashmapGetImage(var2idx, (void*)probdata->vars[address]);
               assert(varidx >= 0);
               assert(varidx < nvars);
               assert(vars[varidx] == probdata->vars[address]);
            }

            SCIP_CALL( SCIPexprCreate(blkmem, &term2, SCIP_EXPR_VARIDX, varidx) );
            SCIP_CALL( exprAdd(blkmem, &e, 1.0, term1, -1.0, term2, 0.0) );

            break;
         }

         case nlSubI: /* subtract immediate */
         {
            SCIPdebugPrintf("subtract constant %g\n", constants[address]);

            assert(stackpos >= 1);
            term1 = stack[stackpos-1];
            --stackpos;

            SCIP_CALL( exprAdd(blkmem, &e, 1.0, term1, 1.0, NULL, -constants[address]) );

            break;
         }

         case nlMul: /* multiply */
         {
            SCIPdebugPrintf("subtract\n");

            assert(stackpos >= 2);
            term1 = stack[stackpos-1];
            --stackpos;
            term2 = stack[stackpos-1];
            --stackpos;

            SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_MUL, term2, term1) );
            break;
         }

         case nlMulV: /* multiply variable */
         {
            address = gmoGetjSolver(gmo, address);
            SCIPdebugPrintf("multiply variable %d = <%s>\n", address, SCIPvarGetName(probdata->vars[address]));

            assert(stackpos >= 1);
            term1 = stack[stackpos-1];
            --stackpos;

            if( !SCIPhashmapExists(var2idx, probdata->vars[address]) )
            {
               /* add variable to list of variables */
               SCIP_CALL( ensureVarsSize(scip, &vars, &varssize, nvars+1) );
               assert(nvars < varssize);
               vars[nvars] = probdata->vars[address];
               varidx = nvars;
               ++nvars;
               SCIP_CALL( SCIPhashmapInsert(var2idx, (void*)vars[varidx], (void*)(size_t)varidx) );
            }
            else
            {
               varidx = (int)(size_t)SCIPhashmapGetImage(var2idx, (void*)probdata->vars[address]);
               assert(varidx >= 0);
               assert(varidx < nvars);
               assert(vars[varidx] == probdata->vars[address]);
            }

            SCIP_CALL( SCIPexprCreate(blkmem, &term2, SCIP_EXPR_VARIDX, varidx) );
            SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_MUL, term1, term2) );
            break;
         }

         case nlMulI: /* multiply immediate */
         {
            SCIPdebugPrintf("multiply constant %g\n", constants[address]);

            assert(stackpos >= 1);
            term1 = stack[stackpos-1];
            --stackpos;

            SCIP_CALL( exprAdd(blkmem, &e, constants[address], term1, 1.0, NULL, 0.0) );

            break;
         }

         case nlMulIAdd:
         {
            SCIPdebugPrintf("multiply constant %g and add\n", constants[address]);

            assert(stackpos >= 2);
            term1 = stack[stackpos-1];
            --stackpos;
            term2 = stack[stackpos-1];
            --stackpos;

            SCIP_CALL( exprAdd(blkmem, &e, constants[address], term1, 1.0, term2, 0.0) );

            break;
         }

#if 1
         case nlDiv: /* divide */
         {
            SCIPdebugPrintf("divide\n");

            assert(stackpos >= 2);
            term1 = stack[stackpos-1];
            --stackpos;
            term2 = stack[stackpos-1];
            --stackpos;

            SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_DIV, term2, term1) );
            break;
         }
#elif 1
         case nlDiv: /* divide */
         {
            SCIP_EXPRDATA_MONOMIAL* monomial;
            SCIP_EXPR* children[2];
            double exps[2] = { -1.0, 1.0 };

            SCIPdebugPrintf("divide\n");

            assert(stackpos >= 2);
            children[0] = stack[stackpos-1];
            --stackpos;
            children[1] = stack[stackpos-1];
            --stackpos;

            SCIP_CALL( SCIPexprCreateMonomial(blkmem, &monomial, 1.0, 2, NULL, exps) );
            SCIP_CALL( SCIPexprCreatePolynomial(blkmem, &e, 2, children, 1, &monomial, 0.0, FALSE) );
            break;
         }
#else
         case nlDiv: /* divide */
         {
            SCIPdebugPrintf("divide\n");

            assert(stackpos >= 2);
            term1 = stack[stackpos-1];
            --stackpos;
            term2 = stack[stackpos-1];
            --stackpos;

            SCIP_CALL( SCIPexprCreate(blkmem, &term1, SCIP_EXPR_INTPOWER, term1, -1) );

            SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_MUL, term2, term1) );
            break;
         }
#endif

         case nlDivV: /* divide variable */
         {
            address = gmoGetjSolver(gmo, address);
            SCIPdebugPrintf("divide variable %d = <%s>\n", address, SCIPvarGetName(probdata->vars[address]));

            assert(stackpos >= 1);
            term1 = stack[stackpos-1];
            --stackpos;

            if( !SCIPhashmapExists(var2idx, probdata->vars[address]) )
            {
               /* add variable to list of variables */
               SCIP_CALL( ensureVarsSize(scip, &vars, &varssize, nvars+1) );
               assert(nvars < varssize);
               vars[nvars] = probdata->vars[address];
               varidx = nvars;
               ++nvars;
               SCIP_CALL( SCIPhashmapInsert(var2idx, (void*)vars[varidx], (void*)(size_t)varidx) );
            }
            else
            {
               varidx = (int)(size_t)SCIPhashmapGetImage(var2idx, (void*)probdata->vars[address]);
               assert(varidx >= 0);
               assert(varidx < nvars);
               assert(vars[varidx] == probdata->vars[address]);
            }

            SCIP_CALL( SCIPexprCreate(blkmem, &term2, SCIP_EXPR_VARIDX, varidx) );
            SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_DIV, term1, term2) );
            break;
         }

         case nlDivI: /* divide immediate */
         {
            SCIPdebugPrintf("divide constant %g\n", constants[address]);
            assert(constants[address] != 0.0);

            assert(stackpos >= 1);
            term1 = stack[stackpos-1];
            --stackpos;

            SCIP_CALL( exprAdd(blkmem, &e, 1.0/constants[address], term1, 1.0, NULL, 0.0) );

            break;
         }

         case nlUMin: /* unary minus */
         {
            SCIPdebugPrintf("negate\n");

            assert(stackpos >= 1);
            term1 = stack[stackpos-1];
            --stackpos;

            SCIP_CALL( exprAdd(blkmem, &e, -1.0, term1, 1.0, NULL, 0.0) );

            break;
         }

         case nlUMinV: /* unary minus variable */
         {
            SCIP_Real minusone;

            address = gmoGetjSolver(gmo, address);
            SCIPdebugPrintf("push negated variable %d = <%s>\n", address, SCIPvarGetName(probdata->vars[address]));

            if( !SCIPhashmapExists(var2idx, probdata->vars[address]) )
            {
               /* add variable to list of variables */
               SCIP_CALL( ensureVarsSize(scip, &vars, &varssize, nvars+1) );
               assert(nvars < varssize);
               vars[nvars] = probdata->vars[address];
               varidx = nvars;
               ++nvars;
               SCIP_CALL( SCIPhashmapInsert(var2idx, (void*)vars[varidx], (void*)(size_t)varidx) );
            }
            else
            {
               varidx = (int)(size_t)SCIPhashmapGetImage(var2idx, (void*)probdata->vars[address]);
               assert(varidx >= 0);
               assert(varidx < nvars);
               assert(vars[varidx] == probdata->vars[address]);
            }

            minusone = -1.0;
            SCIP_CALL( SCIPexprCreate(blkmem, &term1, SCIP_EXPR_VARIDX, varidx) );
            SCIP_CALL( SCIPexprCreateLinear(blkmem, &e, 1, &term1, &minusone, 0.0) );

            break;
         }

         case nlFuncArgN:
         {
            SCIPdebugPrintf("set number of arguments = %d\n", address);
            nargs = address;
            break;
         }

         case nlCallArg1:
         case nlCallArg2:
         case nlCallArgN:
         {
            GamsFuncCode func;

            SCIPdebugPrintf("call function ");

            func = (GamsFuncCode)(address+1); /* here the shift by one was not a good idea */

            switch( func )
            {
               case fnmin:
               {
                  SCIPdebugPrintf("min\n");

                  assert(stackpos >= 2);
                  term1 = stack[stackpos-1];
                  --stackpos;
                  term2 = stack[stackpos-1];
                  --stackpos;

                  SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_MIN, term1, term2) );
                  break;
               }

               case fnmax:
               {
                  SCIPdebugPrintf("max\n");

                  assert(stackpos >= 2);
                  term1 = stack[stackpos-1];
                  --stackpos;
                  term2 = stack[stackpos-1];
                  --stackpos;

                  SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_MAX, term1, term2) );
                  break;
               }

               case fnsqr:
               {
                  SCIPdebugPrintf("square\n");

                  assert(stackpos >= 1);
                  term1 = stack[stackpos-1];
                  --stackpos;

                  SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_SQUARE, term1) );
                  break;
               }

               case fnexp:
               {
                  SCIPdebugPrintf("exp\n");

                  assert(stackpos >= 1);
                  term1 = stack[stackpos-1];
                  --stackpos;

                  SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_EXP, term1) );
                  break;
               }

               case fnlog:
               {
                  SCIPdebugPrintf("log\n");

                  assert(stackpos >= 1);
                  term1 = stack[stackpos-1];
                  --stackpos;

                  SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_LOG, term1) );
                  break;
               }

               case fnlog10:
               case fnsllog10:
               case fnsqlog10:
               {
                  SCIPdebugPrintf("log10 = ln * 1/ln(10)\n");

                  assert(stackpos >= 1);
                  term1 = stack[stackpos-1];
                  --stackpos;

                  SCIP_CALL( SCIPexprCreate(blkmem, &term2, SCIP_EXPR_LOG, term1) );
                  SCIP_CALL( SCIPexprCreate(blkmem, &term1, SCIP_EXPR_CONST, 1.0/log(10.0)) );
                  SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_MUL, term2, term1) );
                  break;
               }

               case fnlog2:
               {
                  SCIPdebugPrintf("log2 = ln * 1/ln(2)\n");

                  assert(stackpos >= 1);
                  term1 = stack[stackpos-1];
                  --stackpos;

                  SCIP_CALL( SCIPexprCreate(blkmem, &term2, SCIP_EXPR_LOG, term1) );
                  SCIP_CALL( SCIPexprCreate(blkmem, &term1, SCIP_EXPR_CONST, 1.0/log(2.0)) );
                  SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_MUL, term2, term1) );
                  break;
               }

               case fnsqrt:
               {
                  SCIPdebugPrintf("sqrt\n");

                  assert(stackpos >= 1);
                  term1 = stack[stackpos-1];
                  --stackpos;

                  SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_SQRT, term1) );
                  break;
               }
#if 0
               case fncos:
               {
                  SCIPdebugPrintf("cos\n");

                  assert(stackpos >= 1);
                  term1 = stack[stackpos-1];
                  --stackpos;

                  SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_COS, term1) );
                  break;
               }

               case fnsin:
               {
                  SCIPdebugPrintf("sin\n");

                  assert(stackpos >= 1);
                  term1 = stack[stackpos-1];
                  --stackpos;

                  SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_SIN, term1) );
                  break;
               }

               case fntan:
               {
                  SCIPdebugPrintf("tan\n");

                  assert(stackpos >= 1);
                  term1 = stack[stackpos-1];
                  --stackpos;

                  SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_TAN, term1) );
                  break;
               }
#endif
               case fnpower: /* x ^ y */
               case fnrpower: /* x ^ y */
               case fncvpower: /* constant ^ x */
               case fnvcpower: /* x ^ constant */
               {
                  SCIPdebugPrintf("power\n");

                  assert(stackpos >= 2);
                  term1 = stack[stackpos-1];
                  --stackpos;
                  term2 = stack[stackpos-1];
                  --stackpos;

                  assert(func != fncvpower || SCIPexprGetOperator(term2) == SCIP_EXPR_CONST);
                  assert(func != fnvcpower || SCIPexprGetOperator(term1) == SCIP_EXPR_CONST);

                  if( SCIPexprGetOperator(term1) == SCIP_EXPR_CONST )
                  {
                     /* use intpower if exponent is an integer constant, otherwise use realpower */
                     if( SCIPisIntegral(scip, SCIPexprGetOpReal(term1)) )
                     {
                        SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_INTPOWER, term2, (int)SCIPexprGetOpReal(term1)) );
                        SCIPexprFreeDeep(blkmem, &term1);
                     }
                     else
                     {
                        SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_REALPOWER, term2, SCIPexprGetOpReal(term1)) );
                        SCIPexprFreeDeep(blkmem, &term1);
                     }
                  }
                  else
                  {
                     /* term2^term1 = exp(log(term2)*term1) */
                     SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_LOG, term2) );
                     SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_MUL, e, term1) );
                     SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_EXP, e) );
                  }

                  break;
               }

               case fnsignpower: /* sign(x)*abs(x)^c */
               {
                  SCIPdebugPrintf("signpower\n");

                  assert(stackpos >= 2);
                  term1 = stack[stackpos-1];
                  --stackpos;
                  term2 = stack[stackpos-1];
                  --stackpos;

                  if( SCIPexprGetOperator(term1) == SCIP_EXPR_CONST )
                  {
                     SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_SIGNPOWER, term2, SCIPexprGetOpReal(term1)) );
                     SCIPexprFreeDeep(blkmem, &term1);
                  }
                  else
                  {
                     SCIPerrorMessage("signpower with non-constant exponent not supported.\n");
                     return SCIP_ERROR;
                  }

                  break;
               }

               case fnpi:
               {
                  SCIPdebugPrintf("pi\n");

                  SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_CONST, M_PI) );
                  break;
               }

               case fndiv:
               case fndiv0:
               {
                  SCIPdebugPrintf("divide\n");

                  assert(stackpos >= 2);
                  term1 = stack[stackpos-1];
                  --stackpos;
                  term2 = stack[stackpos-1];
                  --stackpos;

                  SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_DIV, term2, term1) );
                  break;
               }

               case fnslrec:
               case fnsqrec: /* 1/x */
               {
                  SCIPdebugPrintf("reciprocal\n");

                  assert(stackpos >= 1);
                  term1 = stack[stackpos-1];
                  --stackpos;

                  SCIP_CALL( SCIPexprCreate(blkmem, &term2, SCIP_EXPR_CONST, 1.0) );
                  SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_DIV, term2, term1) );
                  break;
               }

               case fnabs:
               {
                  SCIPdebugPrintf("abs\n");

                  assert(stackpos >= 1);
                  term1 = stack[stackpos-1];
                  --stackpos;

                  SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_ABS, term1) );
                  break;
               }

               case fnpoly: /* univariate polynomial */
               {
                  SCIPdebugPrintf("univariate polynomial of degree %d\n", nargs-1);
                  assert(nargs >= 0);
                  switch( nargs )
                  {
                     case 0:
                     {
                        term1 = stack[stackpos-1];
                        --stackpos;

                        SCIPexprFreeDeep(blkmem, &term1);
                        SCIP_CALL( SCIPexprCreate(blkmem, &e, SCIP_EXPR_CONST, 0.0) );
                        break;
                     }

                     case 1: /* "constant" polynomial */
                     {
                        e = stack[stackpos-1];
                        --stackpos;

                        /* delete variable of polynomial */
                        SCIPexprFreeDeep(blkmem, &stack[stackpos-1]);
                        --stackpos;

                        break;
                     }

                     default: /* polynomial is at least linear */
                     {
                        SCIP_EXPRDATA_MONOMIAL** monomials;
                        SCIP_Real exponent;
                        SCIP_Real constant;
                        int nmonomials;
                        int zero;

                        nmonomials = nargs-1;
                        SCIP_CALL( SCIPallocBufferArray(scip, &monomials, nargs-1) );

                        zero = 0;
                        constant = 0.0;
                        for( ; nargs > 0; --nargs )
                        {
                           assert(stackpos > 0);

                           term1 = stack[stackpos-1];
                           assert(SCIPexprGetOperator(term1) == SCIP_EXPR_CONST);

                           if( nargs > 1 )
                           {
                              exponent = (SCIP_Real)(nargs-1);
                              SCIP_CALL( SCIPexprCreateMonomial(blkmem, &monomials[nargs-2], SCIPexprGetOpReal(term1), 1, &zero, &exponent) );
                           }
                           else
                              constant = SCIPexprGetOpReal(term1);

                           SCIPexprFreeDeep(blkmem, &term1);
                           --stackpos;
                        }

                        assert(stackpos > 0);
                        term1 = stack[stackpos-1];
                        --stackpos;

                        SCIP_CALL( SCIPexprCreatePolynomial(blkmem, &e, 1, &term1, nmonomials, monomials, constant, FALSE) );
                        SCIPfreeBufferArray(scip, &monomials);
                     }
                  }
                  nargs = -1;
                  break;
               }

               /* @todo some of these we could also support */
               case fnerrf:
               case fnceil: case fnfloor: case fnround:
               case fnmod: case fntrunc: case fnsign:
               case fnarctan: case fndunfm:
               case fndnorm: case fnerror: case fnfrac: case fnerrorl:
               case fnfact /* factorial */:
               case fnunfmi /* uniform random number */:
               case fnncpf /* fischer: sqrt(x1^2+x2^2+2*x3) */:
               case fnncpcm /* chen-mangasarian: x1-x3*ln(1+exp((x1-x2)/x3))*/:
               case fnentropy /* x*ln(x) */: case fnsigmoid /* 1/(1+exp(-x)) */:
               case fnboolnot: case fnbooland:
               case fnboolor: case fnboolxor: case fnboolimp:
               case fnbooleqv: case fnrelopeq: case fnrelopgt:
               case fnrelopge: case fnreloplt: case fnrelople:
               case fnrelopne: case fnifthen:
               case fnedist /* euclidian distance */:
               case fncentropy /* x*ln((x+d)/(y+d))*/:
               case fngamma: case fnloggamma: case fnbeta:
               case fnlogbeta: case fngammareg: case fnbetareg:
               case fnsinh: case fncosh: case fntanh:
               case fnncpvusin /* veelken-ulbrich */:
               case fnncpvupow /* veelken-ulbrich */:
               case fnbinomial:
               case fnarccos:
               case fnarcsin: case fnarctan2 /* arctan(x2/x1) */:
               default :
               {
                  SCIPdebugPrintf("nr. %d - unsupported. Error.\n", (int)func);
                  SCIPinfoMessage(scip, NULL, "Error: GAMS function %s not supported.\n", GamsFuncCodeName[func]);
                  return SCIP_READERROR;
               }
            } /*lint !e788*/
            break;
         }

         case nlEnd: /* end of instruction list */
         default:
         {
            SCIPinfoMessage(scip, NULL, "Error: GAMS opcode %s not supported.\n", GamsOpCodeName[opcode]);
            return SCIP_READERROR;
         }
      } /*lint !e788*/

      if( e != NULL )
      {
         if( stackpos >= stacksize )
         {
            stacksize = SCIPcalcMemGrowSize(scip, stackpos+1);
            SCIP_CALL( SCIPreallocBufferArray(scip, &stack, stacksize) );
         }
         assert(stackpos < stacksize);
         stack[stackpos] = e;
         ++stackpos;
      }
   }

   /* there should be exactly one element on the stack, which will be the root of our expression tree */
   assert(stackpos == 1);

   SCIP_CALL( SCIPexprtreeCreate(blkmem, exprtree, stack[0], nvars, 0, NULL) );
   SCIP_CALL( SCIPexprtreeSetVars(*exprtree, nvars, vars) );

   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &stack);
   SCIPhashmapFree(&var2idx);

   return SCIP_OKAY;
}

/** creates a SCIP problem from a GMO */
SCIP_RETCODE SCIPcreateProblemReaderGmo(
   SCIP*                 scip,               /**< SCIP data structure */
   gmoRec_t*             gmo,                /**< GAMS Model Object */
   const char*           indicatorfile,      /**< name of file with indicator specification, or NULL */
   int                   mipstart            /**< how to pass initial variable levels from GMO to SCIP */
)
{
   char buffer[GMS_SSSIZE];
   gevHandle_t gev;
   SCIP_Bool objnonlinear;
   SCIP_VAR** vars;
   SCIP_Real minprior;
   SCIP_Real maxprior;
   int i;
   SCIP_Real* coefs = NULL;
   int* indices = NULL;
   int* nlflag;
   SCIP_VAR** consvars = NULL;
   SCIP_VAR** quadvars1;
   SCIP_VAR** quadvars2;
   SCIP_Real* quadcoefs;
   int* qrow;
   int* qcol;
   SCIP_CONS* con;
   int numSos1, numSos2, nzSos;
   SCIP_PROBDATA* probdata;
   int* opcodes;
   int* fields;
   SCIP_Real* constants;
   int nindics;
   int* indicrows;
   int* indiccols;
   int* indiconvals;
   int indicidx;
   size_t namemem;
   SCIP_RETCODE rc = SCIP_OKAY;
   
   assert(scip != NULL);
   assert(gmo != NULL);

   gev = gmoEnvironment(gmo);
   assert(gev != NULL);

   /* we want a real objective function, if it is linear, otherwise keep the GAMS single-variable-objective? */
   gmoObjReformSet(gmo, 1);
   gmoObjStyleSet(gmo, (int)gmoObjType_Fun);
#if 0
   if( gmoObjNLNZ(gmo) > 0 )
      gmoObjStyleSet(gmo, gmoObjType_Var);
   objnonlinear = FALSE;
#else
   objnonlinear = gmoObjStyle(gmo) == (int)gmoObjType_Fun && gmoObjNLNZ(gmo) > 0;
#endif

   /* we want to start indexing at 0 */
   gmoIndexBaseSet(gmo, 0);

   /* we want GMO to use SCIP's value for infinity */
   gmoPinfSet(gmo,  SCIPinfinity(scip));
   gmoMinfSet(gmo, -SCIPinfinity(scip));

   /* create SCIP problem */
   SCIP_CALL( SCIPallocMemory(scip, &probdata) );
   BMSclearMemory(probdata);

   (void) gmoNameInput(gmo, buffer);
   SCIP_CALL( SCIPcreateProb(scip, buffer,
      probdataDelOrigGmo, probdataTransGmo, probdataDelTransGmo,
      probdataInitSolGmo, probdataExitSolGmo, probdataCopyGmo,
      probdata) );

   /* initialize QMaker, if nonlinear */
   if( gmoNLNZ(gmo) > 0 || objnonlinear )
      gmoUseQSet(gmo, 1);

   /* get data on indicator constraints from options object */
   nindics = 0;
   indicrows = NULL;
   indiccols = NULL;
   indiconvals = NULL;
#ifndef WITH_GAMS
   if( indicatorfile != NULL && *indicatorfile != '\0' )
   {
      optHandle_t opt;
      int itype;

      if( !optCreate(&opt, buffer, sizeof(buffer)) )
      {
         SCIPerrorMessage("*** Could not create optionfile handle: %s\n", buffer);
         return SCIP_ERROR;
      }

#if GMOAPIVERSION < 13
      (void) gevGetStrOpt(gev, gevNameSysDir, buffer);
      strcat(buffer, "optscip.def");
      if( optReadDefinition(opt, buffer) )
#else
      if( optReadDefinitionFromPChar(opt, (char*)"indic indicator\ngeneral group 1 1 Dot options and indicators") )
#endif
      {
         for( i = 1; i <= optMessageCount(opt); ++i )
         {
            optGetMessage(opt, i, buffer, &itype);
            if( itype <= (int) optMsgFileLeave || itype == (int) optMsgUserError )
               gevLogStat(gev, buffer);
         }
         optClearMessages(opt);
         return SCIP_ERROR;
      }

      (void) optReadParameterFile(opt, indicatorfile);
      for( i = 1; i <= optMessageCount(opt); ++i )
      {
         optGetMessage(opt, i, buffer, &itype);
         if( itype <= (int) optMsgFileLeave || itype == (int) optMsgUserError )
            gevLogStat(gev, buffer);
      }
      optClearMessages(opt);

      if( optIndicatorCount(opt, &i) > 0 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &indicrows, gmoM(gmo)) );
         SCIP_CALL( SCIPallocBufferArray(scip, &indiccols, gmoM(gmo)) );
         SCIP_CALL( SCIPallocBufferArray(scip, &indiconvals, gmoM(gmo)) );
         if( gmoGetIndicatorMap(gmo, opt, 1, &nindics, indicrows, indiccols, indiconvals) != 0 )
         {
            SCIPerrorMessage("failed to get indicator mapping\n");
            return SCIP_ERROR;
         }
      }

      (void) optFree(&opt);
   }
#endif
   assert(indicrows != NULL || nindics == 0);
   assert(indiccols != NULL || nindics == 0);
   assert(indiconvals != NULL || nindics == 0);

   namemem = (gmoN(gmo) + gmoM(gmo)) * sizeof(char*);

   probdata->nvars = gmoN(gmo);
   SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->vars, probdata->nvars) ); /*lint !e666*/
   vars = probdata->vars;
   
   /* compute range of variable priorities */ 
   minprior = SCIPinfinity(scip);
   maxprior = 0.0;
   if( gmoPriorOpt(gmo) && gmoNDisc(gmo) > 0 )
   {
      for (i = 0; i < gmoN(gmo); ++i)
      {
         if( gmoGetVarTypeOne(gmo, i) == (int) gmovar_X )
            continue; /* GAMS forbids branching priorities for continuous variables */
         if( gmoGetVarPriorOne(gmo,i) < minprior )
            minprior = gmoGetVarPriorOne(gmo,i);
         if( gmoGetVarPriorOne(gmo,i) > maxprior )
            maxprior = gmoGetVarPriorOne(gmo,i);
      }
   }

   /* get objective functions coefficients */
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, gmoN(gmo)+1) ); /* +1 if we have to transform the objective into a constraint */
   if( !objnonlinear )
   {
      if( gmoObjStyle(gmo) == (int) gmoObjType_Fun )
         (void) gmoGetObjVector(gmo, coefs, NULL);
      else
         coefs[gmoObjVar(gmo)] = 1.0;
   }
   else
   {
      BMSclearMemoryArray(coefs, gmoN(gmo));
   }
   
   /* add variables */
   for( i = 0; i < gmoN(gmo); ++i )
   {
      SCIP_VARTYPE vartype;
      SCIP_Real lb, ub;
      lb = gmoGetVarLowerOne(gmo, i);
      ub = gmoGetVarUpperOne(gmo, i);
      switch( gmoGetVarTypeOne(gmo, i) )
      {
         case gmovar_SC:
            lb = 0.0;
            /*lint -fallthrough*/
         case gmovar_X:
         case gmovar_S1:
         case gmovar_S2:
            vartype = SCIP_VARTYPE_CONTINUOUS;
            break;
         case gmovar_B:
            vartype = SCIP_VARTYPE_BINARY;
            break;
         case gmovar_SI:
            lb = 0.0;
            /*lint -fallthrough*/
         case gmovar_I:
            vartype = SCIP_VARTYPE_INTEGER;
            break;
         default:
            SCIPerrorMessage("Unknown variable type.\n");
            return SCIP_INVALIDDATA;
      }
      if( gmoDict(gmo) )
      {
         (void) gmoGetVarNameOne(gmo, i, buffer);
         if( nindics == 0 )
            namemem += strlen(buffer) + 1;
      }
      else
         sprintf(buffer, "x%d", i);
      SCIP_CALL( SCIPcreateVar(scip, &vars[i], buffer, lb, ub, coefs[i], vartype, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, vars[i]) );
      SCIPdebugMessage("added variable ");
      SCIPdebug( SCIPprintVar(scip, vars[i], NULL) );
      
      if( gmoPriorOpt(gmo) && minprior < maxprior && gmoGetVarTypeOne(gmo, i) != (int) gmovar_X )
      {
         /* in GAMS: higher priorities are given by smaller .prior values
            in SCIP: variables with higher branch priority are always preferred to variables with lower priority in selection of branching variable
            thus, we scale the values from GAMS to lie between 0 (lowest prior) and 1000 (highest prior)
         */
         int branchpriority = (int)(1000.0 / (maxprior - minprior) * (maxprior - gmoGetVarPriorOne(gmo, i)));
         SCIP_CALL( SCIPchgVarBranchPriority(scip, vars[i], branchpriority) );
      }
   }
   
   /* setup bound disjunction constraints for semicontinuous/semiinteger variables by saying x <= 0 or x >= gmoGetVarLower */
   if( gmoGetVarTypeCnt(gmo, (int) gmovar_SC) || gmoGetVarTypeCnt(gmo, (int) gmovar_SI) )
   {
      SCIP_BOUNDTYPE bndtypes[2];
      SCIP_Real      bnds[2];
      SCIP_VAR*      bndvars[2];
      SCIP_CONS*     cons;
      char           name[SCIP_MAXSTRLEN];
      
      bndtypes[0] = SCIP_BOUNDTYPE_UPPER;
      bndtypes[1] = SCIP_BOUNDTYPE_LOWER;
      bnds[0] = 0;

      for( i = 0; i < gmoN(gmo); ++i )
      {
         if( gmoGetVarTypeOne(gmo, i) != (int) gmovar_SC && gmoGetVarTypeOne(gmo, i) != (int) gmovar_SI )
            continue;
         
         bndvars[0] = vars[i];
         bndvars[1] = vars[i];
         bnds[1] = gmoGetVarLowerOne(gmo, i);
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "semi%s_%s", gmoGetVarTypeOne(gmo, i) == (int) gmovar_SC ? "con" : "int", SCIPvarGetName(vars[i]));
         SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &cons, name, 2, bndvars, bndtypes, bnds,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIPdebugMessage("added constraint ");
         SCIPdebug( SCIPprintCons(scip, cons, NULL) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }
   
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, gmoN(gmo)+1) ); /* +1 if we have to transform the objective into a constraint */

   /* setup SOS constraints */
   gmoGetSosCounts(gmo, &numSos1, &numSos2, &nzSos);
   if( nzSos > 0 )
   {
      int numSos;
      int* sostype;
      int* sosbeg;
      int* sosind;
      double* soswt;
      int j, k;
      
      numSos = numSos1 + numSos2;
      SCIP_CALL( SCIPallocBufferArray(scip, &sostype, numSos) );
      SCIP_CALL( SCIPallocBufferArray(scip, &sosbeg,  numSos+1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &sosind,  nzSos) );
      SCIP_CALL( SCIPallocBufferArray(scip, &soswt,   nzSos) );

      (void) gmoGetSosConstraints(gmo, sostype, sosbeg, sosind, soswt);
      
      for( i = 0; i < numSos; ++i )
      {
         for( j = sosbeg[i], k = 0; j < sosbeg[i+1]; ++j, ++k )
         {
            consvars[k] = vars[sosind[j]];
            assert(gmoGetVarTypeOne(gmo, sosind[j]) == (sostype[i] == 1 ? (int) gmovar_S1 : (int) gmovar_S2));
         }
         
         sprintf(buffer, "sos%d", i);
         if( sostype[i] == 1 )
         {
            SCIP_CALL( SCIPcreateConsSOS1(scip, &con, buffer, k, consvars, &soswt[sosbeg[i]], TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
         }
         else
         {
            assert(sostype[i] == 2);
            SCIP_CALL( SCIPcreateConsSOS2(scip, &con, buffer, k, consvars, &soswt[sosbeg[i]], TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
         }
         
         SCIP_CALL( SCIPaddCons(scip, con) );
         SCIPdebugMessage("added constraint ");
         SCIPdebug( SCIPprintCons(scip, con, NULL) );
         SCIP_CALL( SCIPreleaseCons(scip, &con) );
      }
      
      SCIPfreeBufferArray(scip, &sostype);
      SCIPfreeBufferArray(scip, &sosbeg);
      SCIPfreeBufferArray(scip, &sosind);
      SCIPfreeBufferArray(scip, &soswt);
   }
   
   /* setup regular constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &indices, gmoN(gmo)) );
   indicidx = 0;
   
   /* alloc some memory, if nonlinear */
   if( gmoNLNZ(gmo) > 0 || objnonlinear )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &nlflag, gmoN(gmo)) );

      SCIP_CALL( SCIPallocBufferArray(scip, &quadvars1, gmoMaxQNZ(gmo)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &quadvars2, gmoMaxQNZ(gmo)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &quadcoefs, gmoMaxQNZ(gmo)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &qrow, gmoMaxQNZ(gmo)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &qcol, gmoMaxQNZ(gmo)) );

      SCIP_CALL( SCIPallocBufferArray(scip, &opcodes, gmoNLCodeSizeMaxRow(gmo)+1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &fields, gmoNLCodeSizeMaxRow(gmo)+1) );
      SCIP_CALL( SCIPduplicateBufferArray(scip, &constants, gmoPPool(gmo), gmoNLConst(gmo)) );

      /* translate special GAMS constants into SCIP variants (gmo does not seem to do this...) */
      for( i = 0; i < gmoNLConst(gmo); ++i )
      {
         if( constants[i] == GMS_SV_PINF )
            constants[i] =  SCIPinfinity(scip);
         else if( constants[i] == GMS_SV_MINF )
            constants[i] = -SCIPinfinity(scip);
         else if( constants[i] == GMS_SV_EPS )
            constants[i] = 0.0;
         else if( constants[i] == GMS_SV_UNDEF || constants[i] == GMS_SV_NA || constants[i] == GMS_SV_NAINT || constants[i] == GMS_SV_ACR )
         {
            SCIPwarningMessage(scip, "Constant %e in nonlinear expressions constants pool cannot be handled by SCIP.\n");
            constants[i] = SCIP_INVALID;
         }
         else if( constants[i] <= -SCIPinfinity(scip) )
            constants[i] = -SCIPinfinity(scip);
         else if( constants[i] >=  SCIPinfinity(scip) )
            constants[i] =  SCIPinfinity(scip);
      }
   }
   else
   {
      nlflag = NULL;

      quadvars1 = NULL;
      quadvars2 = NULL;
      quadcoefs = NULL;
      qrow = NULL;
      qcol = NULL;

      opcodes = NULL;
      fields = NULL;
      constants = NULL;
   }

   for( i = 0; i < gmoM(gmo); ++i )
   {
      double lhs;
      double rhs;
      switch( gmoGetEquTypeOne(gmo, i) )
      {
         case gmoequ_E:
            lhs = rhs = gmoGetRhsOne(gmo, i);
            break;
         case gmoequ_G:
            lhs = gmoGetRhsOne(gmo, i);
            rhs = SCIPinfinity(scip);
            break;
         case gmoequ_L:
            lhs = -SCIPinfinity(scip);
            rhs = gmoGetRhsOne(gmo, i);
            break;
         case gmoequ_N:
            lhs = -SCIPinfinity(scip);
            rhs =  SCIPinfinity(scip);
            break;
         case gmoequ_X:
            SCIPerrorMessage("External functions not supported by SCIP.\n");
            return SCIP_INVALIDDATA;
         case gmoequ_C:
            SCIPerrorMessage("Conic constraints not supported by SCIP interface yet.\n");
            return SCIP_INVALIDDATA;
         case gmoequ_B:
            SCIPerrorMessage("Logic constraints not supported by SCIP interface yet.\n");
            return SCIP_INVALIDDATA;
         default:
            SCIPerrorMessage("unknown equation type.\n");
            return SCIP_INVALIDDATA;
      }

      if( gmoDict(gmo) )
      {
         (void) gmoGetEquNameOne(gmo, i, buffer);
         if( nindics == 0 )
            namemem += strlen(buffer) + 1;
      }
      else
         sprintf(buffer, "e%d", i);

      con = NULL;
      switch( gmoGetEquOrderOne(gmo, i) )
      {
         case gmoorder_L:
         {
            /* linear constraint */
            int j, nz, nlnz;
            (void) gmoGetRowSparse(gmo, i, indices, coefs, NULL, &nz, &nlnz);
            assert(nlnz == 0);

            for( j = 0; j < nz; ++j )
               consvars[j] = vars[indices[j]];

            /* create indicator constraint, if we are at one */
            if( indicidx < nindics && indicrows[indicidx] == i ) /*lint !e613*/
            {
               SCIP_VAR* binvar;

               binvar = vars[indiccols[indicidx]]; /*lint !e613*/
               if( SCIPvarGetType(binvar) != SCIP_VARTYPE_BINARY )
               {
                  SCIPerrorMessage("Indicator variable <%s> is not of binary type.\n", SCIPvarGetName(binvar));
                  return SCIP_ERROR;
               }

               assert(indiconvals[indicidx] == 0 || indiconvals[indicidx] == 1); /*lint !e613*/
               if( indiconvals[indicidx] == 0 ) /*lint !e613*/
               {
                  SCIP_CALL( SCIPgetNegatedVar(scip, binvar, &binvar) );
               }

               if( !SCIPisInfinity(scip, rhs) )
               {
                  SCIP_CALL( SCIPcreateConsIndicator(scip, &con, buffer, binvar, nz, consvars, coefs, rhs,
                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

                  if( !SCIPisInfinity(scip, -lhs) )
                  {
                     SCIP_CALL( SCIPaddCons(scip, con) );
                     SCIPdebugMessage("added constraint ");
                     SCIPdebug( SCIPprintCons(scip, con, NULL) );
                     SCIP_CALL( SCIPreleaseCons(scip, &con) );
                     con = NULL;
                  }
               }
               if( !SCIPisInfinity(scip, -lhs) )
               {
                  for( j = 0; j < nz; ++j )
                     coefs[j] = -coefs[j];
                  SCIP_CALL( SCIPcreateConsIndicator(scip, &con, buffer, binvar, nz, consvars, coefs, -lhs,
                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
               }

               ++indicidx;
            }
            else
            {
               SCIP_CALL( SCIPcreateConsLinear(scip, &con, buffer, nz, consvars, coefs, lhs, rhs,
                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
            }

            break;
         }
         
         case gmoorder_Q:
         {
            /* quadratic constraint */
            int j, nz, nlnz, qnz;
            
            assert(qcol != NULL);
            assert(qrow != NULL);
            assert(quadcoefs != NULL);
            assert(quadvars1 != NULL);
            assert(quadvars2 != NULL);

            (void) gmoGetRowSparse(gmo, i, indices, coefs, NULL, &nz, &nlnz);
            for( j = 0; j < nz; ++j )
               consvars[j] = vars[indices[j]];
            
            qnz = gmoGetRowQNZOne(gmo,i);
            (void) gmoGetRowQ(gmo, i, qcol, qrow, quadcoefs);
            for( j = 0; j < qnz; ++j )
            {
               assert(qcol[j] >= 0);
               assert(qrow[j] >= 0);
               assert(qcol[j] < gmoN(gmo));
               assert(qrow[j] < gmoN(gmo));
               quadvars1[j] = vars[qcol[j]];
               quadvars2[j] = vars[qrow[j]];
               if( qcol[j] == qrow[j] )
                  quadcoefs[j] /= 2.0; /* for some strange reason, the coefficients on the diagonal are multiplied by 2 in GMO. */
            }

            SCIP_CALL( SCIPcreateConsQuadratic(scip, &con, buffer, nz, consvars, coefs, qnz, quadvars1, quadvars2, quadcoefs, lhs, rhs,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
            break;
         }

         case gmoorder_NL:
         {
            /* nonlinear constraint */
            int j, nz, nlnz, linnz;
            int codelen;
            SCIP_EXPRTREE* exprtree;

            assert(nlflag != NULL);

            (void) gmoGetRowSparse(gmo, i, indices, coefs, nlflag, &nz, &nlnz);
            linnz = 0;
            for( j = 0; j < nz; ++j )
            {
               if( !nlflag[j] )
               {
                  consvars[linnz] = vars[indices[j]];
                  coefs[linnz] = coefs[j];
                  ++linnz;
               }
            }

            (void) gmoDirtyGetRowFNLInstr(gmo, i, &codelen, opcodes, fields);
            rc = makeExprtree(scip, gmo, codelen, opcodes, fields, constants, &exprtree);
            if( rc == SCIP_READERROR )
            {
               SCIPinfoMessage(scip, NULL, "Error processing nonlinear instructions of equation %s.\n", buffer);
               goto TERMINATE;
            }
            SCIP_CALL( rc );

            SCIP_CALL( SCIPcreateConsNonlinear(scip, &con, buffer, linnz, consvars, coefs, 1, &exprtree, NULL, lhs, rhs,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

            SCIP_CALL( SCIPexprtreeFree(&exprtree) );
            break;
         }

         default:
            SCIPerrorMessage("Unexpected equation order.\n");
            return SCIP_INVALIDDATA;
      }
      
      assert(con != NULL);
      SCIP_CALL( SCIPaddCons(scip, con) );      
      SCIPdebugMessage("added constraint ");
      SCIPdebug( SCIPprintCons(scip, con, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &con) );

      /* @todo do something about this */
      if( indicidx < nindics && indicrows[indicidx] == i ) /*lint !e613*/
      {
         SCIPerrorMessage("Only linear constraints can be indicatored, currently.\n");
         return SCIP_ERROR;
      }
   }
   
   if( objnonlinear )
   {
      /* make constraint out of nonlinear objective function */
      int j, nz, nlnz, qnz;
      double lhs, rhs;
      
      assert(gmoGetObjOrder(gmo) == (int) gmoorder_L || gmoGetObjOrder(gmo) == (int) gmoorder_Q || gmoGetObjOrder(gmo) == (int) gmoorder_NL);

      SCIP_CALL( SCIPcreateVar(scip, &probdata->objvar, "xobj", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->objvar) );
      SCIPdebugMessage("added objective variable ");
      SCIPdebug( SCIPprintVar(scip, probdata->objvar, NULL) );

      if( gmoGetObjOrder(gmo) != (int) gmoorder_NL )
      {
         assert(qcol != NULL);
         assert(qrow != NULL);
         assert(quadcoefs != NULL);
         assert(quadvars1 != NULL);
         assert(quadvars2 != NULL);

         (void) gmoGetObjSparse(gmo, indices, coefs, NULL, &nz, &nlnz);
         for( j = 0; j < nz; ++j )
            consvars[j] = vars[indices[j]];

         consvars[nz] = probdata->objvar;
         coefs[nz] = -1.0;
         ++nz;

         qnz = gmoObjQNZ(gmo);
         (void) gmoGetObjQ(gmo, qcol, qrow, quadcoefs);
         for( j = 0; j < qnz; ++j )
         {
            assert(qcol[j] >= 0);
            assert(qrow[j] >= 0);
            assert(qcol[j] < gmoN(gmo));
            assert(qrow[j] < gmoN(gmo));
            quadvars1[j] = vars[qcol[j]];
            quadvars2[j] = vars[qrow[j]];
            if( qcol[j] == qrow[j] )
               quadcoefs[j] /= 2.0; /* for some strange reason, the coefficients on the diagonal are multiplied by 2 in GMO */
         }

         if( gmoSense(gmo) == (int) gmoObj_Min )
         {
            lhs = -SCIPinfinity(scip);
            rhs = -gmoObjConst(gmo);
         }
         else
         {
            lhs = -gmoObjConst(gmo);
            rhs = SCIPinfinity(scip);
         }

         SCIP_CALL( SCIPcreateConsQuadratic(scip, &con, "objective", nz, consvars, coefs, qnz, quadvars1, quadvars2, quadcoefs, lhs, rhs,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
      }
      else
      {
         SCIP_Real objfactor;
         int linnz;
         int codelen;
         SCIP_EXPRTREE* exprtree;

         assert(nlflag != NULL);

         (void) gmoGetObjSparse(gmo, indices, coefs, nlflag, &nz, &nlnz);
         linnz = 0;
         for( j = 0; j < nz; ++j )
         {
            if( !nlflag[j] )
            {
               coefs[linnz] = coefs[j];
               consvars[linnz] = vars[indices[j]];
               ++linnz;
            }
         }

         consvars[linnz] = probdata->objvar;
         coefs[linnz] = -1.0;
         ++linnz;

         objfactor = -1.0 / gmoObjJacVal(gmo);

         (void) gmoDirtyGetObjFNLInstr(gmo, &codelen, opcodes, fields);
         rc = makeExprtree(scip, gmo, codelen, opcodes, fields, constants, &exprtree);
         if( rc == SCIP_READERROR )
         {
            SCIPinfoMessage(scip, NULL, "Error processing nonlinear instructions of objective %s.\n", gmoGetObjName(gmo, buffer));
            goto TERMINATE;
         }
         SCIP_CALL( rc );

         if( gmoSense(gmo) == (int) gmoObj_Min )
         {
            lhs = -SCIPinfinity(scip);
            rhs = -gmoObjConst(gmo);
         }
         else
         {
            lhs = -gmoObjConst(gmo);
            rhs = SCIPinfinity(scip);
         }

         SCIP_CALL( SCIPcreateConsNonlinear(scip, &con, "objective", linnz, consvars, coefs, 1, &exprtree, &objfactor, lhs, rhs,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPexprtreeFree(&exprtree) );
      }

      SCIP_CALL( SCIPaddCons(scip, con) );
      SCIPdebugMessage("added objective constraint ");
      SCIPdebug( SCIPprintCons(scip, con, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &con) );
   }
   else if( !SCIPisZero(scip, gmoObjConst(gmo)) )
   {
      /* handle constant term in linear objective by adding a fixed variable */
      SCIP_CALL( SCIPcreateVar(scip, &probdata->objconst, "objconst", 1.0, 1.0, gmoObjConst(gmo), SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->objconst) );
      SCIPdebugMessage("added variable for objective constant: ");
      SCIPdebug( SCIPprintVar(scip, probdata->objconst, NULL) );
   }

   if( gmoSense(gmo) == (int) gmoObj_Max )
      SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );
   
   /* set objective limit, if enabled */
   if( gevGetIntOpt(gev, gevUseCutOff) )
   {
      SCIP_CALL( SCIPsetObjlimit(scip, gevGetDblOpt(gev, gevCutOff)) );
   }

   /* handle initial solution values */
   switch( mipstart )
   {
      case 0 :
      {
         /* don't pass any initial values to SCIP */
         break;
      }

      case 2:
      case 3:
      {
         /* pass all initial values to SCIP and let SCIP check feasibility (2) or repair (3)
          * NOTE: mipstart=3 does not work as expected: heur_completesol does not run if values for all vars are given and for all integer variables integral values are given
          */
         SCIP_SOL* sol;
         SCIP_Real* vals;
         SCIP_Bool stored;

         if( mipstart == 2 )
         {
            /* with this, SCIP will only check feasibility */
            SCIP_CALL( SCIPcreateOrigSol(scip, &sol, NULL) );
         }
         else
         {
            /* with this, SCIP will try to find a feasible solution close by to the initial values */
            SCIP_CALL( SCIPcreatePartialSol(scip, &sol, NULL) );
         }

         SCIP_CALL( SCIPallocBufferArray(scip, &vals, gmoN(gmo)) );
         (void) gmoGetVarL(gmo, vals);

         SCIP_CALL( SCIPsetSolVals(scip, sol, gmoN(gmo), probdata->vars, vals) );

         /* if we have extra variable for objective, then need to set its value too */
         if( probdata->objvar != NULL )
         {
            double objval;
            int numErr;
            (void) gmoEvalFuncObj(gmo, vals, &objval, &numErr);
            if( numErr == 0 )
            {
               SCIP_CALL( SCIPsetSolVal(scip, sol, probdata->objvar, objval) );
            }
         }

         /* if we have extra variable for objective constant, then need to set its value to 1.0 here too */
         if( probdata->objconst != NULL )
         {
            SCIP_CALL( SCIPsetSolVal(scip, sol, probdata->objconst, 1.0) );
         }

         SCIP_CALL( SCIPaddSolFree(scip, &sol, &stored) );
         assert(stored);

         SCIPfreeBufferArray(scip, &vals);

         if( mipstart == 3 )
         {
            SCIPinfoMessage(scip, NULL, "Passed partial solution with values for all variables to SCIP.");
         }

         break;
      }

      case 1:
      case 4:
      {
         /* pass some initial value to SCIP and let SCIP complete solution */
         SCIP_SOL* sol;
         SCIP_Bool stored;
         double tryint = 0.0;
         int nknown;

         if( mipstart == 4 )
            tryint = gevGetDblOpt(gev, gevTryInt);

         SCIP_CALL( SCIPcreatePartialSol(scip, &sol, NULL) );

         nknown = 0;
         for( i = 0; i < gmoN(gmo); ++i )
         {
            if( mipstart == 1 && (gmoGetVarTypeOne(gmo, i) == gmovar_B || gmoGetVarTypeOne(gmo, i) == gmovar_I || gmoGetVarTypeOne(gmo, i) == gmovar_SI) )
            {
               /* 1: set all integer variables */
               SCIP_CALL( SCIPsetSolVal(scip, sol, probdata->vars[i], gmoGetVarLOne(gmo, i)) );
               ++nknown;
            }

            if( mipstart == 4 && (gmoGetVarTypeOne(gmo, i) == gmovar_B || gmoGetVarTypeOne(gmo, i) == gmovar_I || gmoGetVarTypeOne(gmo, i) == gmovar_SI) )
            {
               /* 4: set only integer variables with level close to an integral value, closeness decided by tryint */
               SCIP_Real val;

               val = gmoGetVarLOne(gmo, i);
               if( fabs(round(val)-val) <= tryint )
               {
                  SCIP_CALL( SCIPsetSolVal(scip, sol, probdata->vars[i], val) );
                  ++nknown;
               }
            }
         }

         /* if we have extra variable for objective constant, then can set its value to 1.0 here too */
         if( probdata->objconst != NULL )
         {
            SCIP_CALL( SCIPsetSolVal(scip, sol, probdata->objconst, 1.0) );
         }

         SCIP_CALL( SCIPaddSolFree(scip, &sol, &stored) );
         assert(stored);

         SCIPinfoMessage(scip, NULL, "Passed partial solution with values for %d variables (%.1f%%) to SCIP.", nknown, 100.0*(double)nknown/gmoN(gmo));

         break;
      }

      default:
      {
         SCIPwarningMessage(scip, "Setting mipstart = %d not supported. Ignored.\n", mipstart);
      }
   }

   if( namemem > 1024 * 1024 && nindics == 0 )
   {
      namemem <<= 1;  /* transformed problem has copy of names, so duplicate estimate */
      SCIPinfoMessage(scip, NULL, "Space for names approximately %0.2f MB. Use statement '<modelname>.dictfile=0;' to turn dictionary off.\n", namemem/(1024.0*1024.0));
   }

TERMINATE:
   SCIPfreeBufferArrayNull(scip, &coefs);
   SCIPfreeBufferArrayNull(scip, &indices);
   SCIPfreeBufferArrayNull(scip, &consvars);
   SCIPfreeBufferArrayNull(scip, &nlflag);
   SCIPfreeBufferArrayNull(scip, &quadvars1);
   SCIPfreeBufferArrayNull(scip, &quadvars2);
   SCIPfreeBufferArrayNull(scip, &quadcoefs);
   SCIPfreeBufferArrayNull(scip, &qrow);
   SCIPfreeBufferArrayNull(scip, &qcol);
   SCIPfreeBufferArrayNull(scip, &opcodes);
   SCIPfreeBufferArrayNull(scip, &fields);
   SCIPfreeBufferArrayNull(scip, &constants);
   SCIPfreeBufferArrayNull(scip, &indicrows);
   SCIPfreeBufferArrayNull(scip, &indiccols);
   SCIPfreeBufferArrayNull(scip, &indiconvals);

   /* deinitialize QMaker, if nonlinear */
   if( gmoNLNZ(gmo) > 0 || objnonlinear )
      gmoUseQSet(gmo, 0);

   return rc;
}

/** check solution for feasibility and resolves by NLP solver, if necessary and possible */
static
SCIP_RETCODE checkAndRepairSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            solvals,            /**< solution to check */
   SCIP_Real*            objval,             /**< objective value corresponding to solvals */
   SCIP_Bool             resolvenlp,         /**< whether NLP resolving is allowed */
   SCIP_Bool*            success             /**< to store whether solution is feasible or could be made feasible */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_HEUR* heursubnlp;
   SCIP_SOL* sol;

   assert(scip    != NULL);
   assert(solvals != NULL);
   assert(success != NULL);
   assert(objval  != NULL);

   if( SCIPisTransformed(scip) )
   {
      /* cannot create solutions in SOLVED stage */
      SCIP_CALL( SCIPfreeTransform(scip) );
   }

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPsetSolVals(scip, sol, probdata->nvars, probdata->vars, solvals) );
   if( probdata->objvar != NULL )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, probdata->objvar, *objval) );
   }
   if( probdata->objconst != NULL )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, probdata->objconst, 1.0) );
   }

   SCIP_CALL( SCIPcheckSolOrig(scip, sol, success, FALSE, FALSE) );

   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   if( *success || !resolvenlp )
      return SCIP_OKAY;

   /* assert that we checked already that resolving is possible and makes sense */
   assert(SCIPgetNContVars(scip) > 0);
   assert(SCIPgetNNlpis(scip) > 0);

   heursubnlp = SCIPfindHeur(scip, "subnlp");
   assert(heursubnlp != NULL);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "Attempt solving NLP from original problem with fixed discrete variables.\n");

   /* create transformed problem and recreate sol in transformed problem, so subnlp heuristic can return result in it */
   SCIP_CALL( SCIPtransformProb(scip) );
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPsetSolVals(scip, sol, probdata->nvars, probdata->vars, solvals) );
   if( probdata->objvar != NULL )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, probdata->objvar, *objval) );
   }
   if( probdata->objconst != NULL )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, probdata->objconst, 1.0) );
   }

   SCIP_CALL( SCIPresolveSolHeurSubNlp(scip, heursubnlp, sol, success, 100LL, 10.0) );

   if( *success )
   {
      SCIP_CALL( SCIPcheckSolOrig(scip, sol, success, SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_FULL, FALSE) );

      if( *success )
      {
         SCIP_CALL( SCIPgetSolVals(scip, sol, probdata->nvars, probdata->vars, solvals) );
         *objval = SCIPgetSolOrigObj(scip, sol);

         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "NLP solution is feasible, objective value = %.15e.\n", *objval);
      }
      else
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "NLP solution still not feasible, objective value = %.15e.\n", SCIPgetSolOrigObj(scip, sol));
      }
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "Failed to resolve NLP.\n");
   }

   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   return SCIP_OKAY;
}

/** stores solve information (solution, statistics) in a GMO */
static
SCIP_RETCODE writeGmoSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   gmoHandle_t           gmo                 /**< GAMS Model Object */
)
{
   SCIP_PROBDATA* probdata;
   int nrsol;
   SCIP_Real dualbound;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);
   assert(probdata->vars != NULL);

   nrsol = SCIPgetNSols(scip);

   switch( SCIPgetStatus(scip) )
   {
      default:
      case SCIP_STATUS_UNKNOWN: /* the solving status is not yet known */
         gmoSolveStatSet(gmo, (int) gmoSolveStat_SystemErr);
         gmoModelStatSet(gmo, (int) gmoModelStat_ErrorNoSolution);
         nrsol = 0;
         break;
      case SCIP_STATUS_USERINTERRUPT: /* the user interrupted the solving process (by pressing Ctrl-C) */
         gmoSolveStatSet(gmo, (int) gmoSolveStat_User);
         gmoModelStatSet(gmo, nrsol > 0 ? (gmoNDisc(gmo) ? (int) gmoModelStat_Integer : (int) gmoModelStat_Feasible) : (int) gmoModelStat_NoSolutionReturned);
         break;
      case SCIP_STATUS_NODELIMIT:      /* the solving process was interrupted because the node limit was reached */
      case SCIP_STATUS_STALLNODELIMIT: /* the solving process was interrupted because the node limit was reached */
      case SCIP_STATUS_TOTALNODELIMIT:
         gmoSolveStatSet(gmo, (int) gmoSolveStat_Iteration);
         gmoModelStatSet(gmo, nrsol > 0 ? (gmoNDisc(gmo) ? (int) gmoModelStat_Integer : (int) gmoModelStat_Feasible) : (int) gmoModelStat_NoSolutionReturned);
         break;
      case SCIP_STATUS_TIMELIMIT: /* the solving process was interrupted because the time limit was reached */
      case SCIP_STATUS_MEMLIMIT:  /* the solving process was interrupted because the memory limit was reached */
         gmoSolveStatSet(gmo, (int) gmoSolveStat_Resource);
         gmoModelStatSet(gmo, nrsol > 0 ? (gmoNDisc(gmo) ? (int) gmoModelStat_Integer : (int) gmoModelStat_Feasible) : (int) gmoModelStat_NoSolutionReturned);
         break;
      case SCIP_STATUS_GAPLIMIT: /* the solving process was interrupted because the gap limit was reached */
         gmoSolveStatSet(gmo, (int) gmoSolveStat_Normal);
         gmoModelStatSet(gmo, nrsol > 0 ? (SCIPgetGap(scip) > 0.0 ? (gmoNDisc(gmo) ? (int) gmoModelStat_Integer : (int) gmoModelStat_Feasible) : (int) gmoModelStat_OptimalGlobal): (int) gmoModelStat_NoSolutionReturned);
         break;
      case SCIP_STATUS_SOLLIMIT: /* the solving process was interrupted because the solution limit was reached */
      case SCIP_STATUS_BESTSOLLIMIT: /* the solving process was interrupted because the solution improvement limit was reached */
         gmoSolveStatSet(gmo, (int) gmoSolveStat_Resource);
         gmoModelStatSet(gmo, nrsol > 0 ? (gmoNDisc(gmo) ? (int) gmoModelStat_Integer : (int) gmoModelStat_Feasible) : (int) gmoModelStat_NoSolutionReturned);
         break;
      case SCIP_STATUS_OPTIMAL: /* the problem was solved to optimality, an optimal solution is available */
         gmoSolveStatSet(gmo, (int) gmoSolveStat_Normal);
         gmoModelStatSet(gmo, (int) gmoModelStat_OptimalGlobal);
         break;
      case SCIP_STATUS_INFEASIBLE: /* the problem was proven to be infeasible */
         gmoSolveStatSet(gmo, (int) gmoSolveStat_Normal);
         gmoModelStatSet(gmo, (int) gmoModelStat_InfeasibleNoSolution);
         nrsol = 0;
         break;
      case SCIP_STATUS_UNBOUNDED: /* the problem was proven to be unbounded */
         gmoSolveStatSet(gmo, (int) gmoSolveStat_Normal);
         gmoModelStatSet(gmo, nrsol > 0 ? (int) gmoModelStat_Unbounded : (int) gmoModelStat_UnboundedNoSolution);
         break;
      case SCIP_STATUS_INFORUNBD: /* the problem was proven to be either infeasible or unbounded */
         gmoSolveStatSet(gmo, (int) gmoSolveStat_Normal);
         gmoModelStatSet(gmo, (int) gmoModelStat_NoSolutionReturned);
         nrsol = 0;
         break;
   }

   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING || SCIPgetStage(scip) == SCIP_STAGE_SOLVED )
      dualbound = SCIPgetDualbound(scip);
   else
      dualbound = gmoValNA(gmo);
   gmoSetHeadnTail(gmo, (int) gmoTmipbest, dualbound);
   gmoSetHeadnTail(gmo, (int) gmoTmipnod,  (double) SCIPgetNNodes(scip));
   gmoSetHeadnTail(gmo, (int) gmoHresused, SCIPgetSolvingTime(scip));
   gmoSetHeadnTail(gmo, (int) gmoHiterused, (double) SCIPgetNLPIterations(scip));
   gmoSetHeadnTail(gmo, (int) gmoHdomused, 0.0);

   /* dump all solutions, if more than one found and parameter is set */
   if( nrsol > 1)
   {
      char* indexfilename;

      SCIP_CALL( SCIPgetStringParam(scip, "gams/dumpsolutions", &indexfilename) );
#ifndef WITH_GAMS
      if( indexfilename != NULL && indexfilename[0] )
      {
         char buffer[SCIP_MAXSTRLEN];
         gdxHandle_t gdx;
         int rc;

         if( !gdxCreate(&gdx, buffer, sizeof(buffer)) )
         {
            SCIPerrorMessage("failed to load GDX I/O library: %s\n", buffer);
            return SCIP_OKAY;
         }

         SCIPinfoMessage(scip, NULL, "\nDumping %d alternate solutions:\n", nrsol-1);
         /* create index GDX file */
         if( gdxOpenWrite(gdx, indexfilename, "SCIP DumpSolutions Index File", &rc) == 0 )
         {
            rc = gdxGetLastError(gdx);
            (void) gdxErrorStr(gdx, rc, buffer);
            SCIPerrorMessage("problem writing GDX file %s: %s\n", indexfilename, buffer);
         }
         else
         {
            gdxStrIndexPtrs_t keys;
            gdxStrIndex_t     keysX;
            gdxValues_t       vals;
            SCIP_Real* collev;
            int sloc;
            int i;

            /* create index file */
            GDXSTRINDEXPTRS_INIT(keysX, keys);
            (void) gdxDataWriteStrStart(gdx, "index", "Dumpsolutions index", 1, (int) dt_set, 0);
            for( i = 1; i < nrsol; ++i)
            {
               (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "soln_scip_p%d.gdx", i);
               (void) gdxAddSetText(gdx, buffer, &sloc);
               (void) SCIPsnprintf(keys[0], GMS_SSSIZE, "file%d", i);
               vals[GMS_VAL_LEVEL] = sloc;
               (void) gdxDataWriteStr(gdx, (const char**)keys, vals);
            }
            (void) gdxDataWriteDone(gdx);
            (void) gdxClose(gdx);

            SCIP_CALL( SCIPallocBufferArray(scip, &collev, gmoN(gmo)) );

            /* create point files */
            for( i = 1; i < nrsol; ++i)
            {
               (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "soln_scip_p%d.gdx", i);

               SCIP_CALL( SCIPgetSolVals(scip, SCIPgetSols(scip)[i], probdata->nvars, probdata->vars, collev) );
               (void) gmoSetVarL(gmo, collev);
               if( gmoUnloadSolutionGDX(gmo, buffer, 0, 1, 0) )
               {
                  SCIPerrorMessage("Problems creating point file %s\n", buffer);
               }
               else
               {
                  SCIPinfoMessage(scip, NULL, "Created point file %s\n", buffer);
               }
            }

            SCIPfreeBufferArray(scip, &collev);
         }

         (void) gdxFree(&gdx);
      }

      SCIP_CALL( SCIPgetStringParam(scip, "gams/dumpsolutionsmerged", &indexfilename) );
      if( indexfilename != NULL && indexfilename[0] )
      {
         int solnvarsym;

         if( gmoCheckSolPoolUEL(gmo, "soln_scip_p", &solnvarsym) )
         {
            SCIPerrorMessage("Solution pool scenario label 'soln_scip_p' contained in model dictionary. Cannot dump merged solutions pool.\n");
         }
         else
         {
            void* handle;

            handle = gmoPrepareSolPoolMerge(gmo, indexfilename, nrsol-1, "soln_scip_p");
            if( handle != NULL )
            {
               SCIP_Real* collev;
               int k, i;

               SCIP_CALL( SCIPallocBufferArray(scip, &collev, gmoN(gmo)) );
               for ( k = 0; k < solnvarsym; k++ )
               {
                  gmoPrepareSolPoolNextSym(gmo, handle);
                  for( i = 1; i < nrsol; ++i )
                  {
                     SCIP_CALL( SCIPgetSolVals(scip, SCIPgetSols(scip)[i], probdata->nvars, probdata->vars, collev) );
                     (void) gmoSetVarL(gmo, collev);
                     if( gmoUnloadSolPoolSolution (gmo, handle, i-1) )
                     {
                        SCIPerrorMessage("Problems unloading solution point %d symbol %d\n", i, k);
                     }
                  }
               }
               if( gmoFinalizeSolPoolMerge(gmo, handle) )
               {
                  SCIPerrorMessage("Problems finalizing merged solution pool\n");
               }
            }
            else
            {
               SCIPerrorMessage("Problems preparing merged solution pool\n");
            }
         }
      }
#endif
   }

   /* pass best solution to GMO, if any */
   if( nrsol > 0 )
   {
      SCIP_SOL* sol;
      SCIP_Real* collev;
      SCIP_Real primalbound;

      sol = SCIPgetBestSol(scip);
      assert(sol != NULL);

      SCIP_CALL( SCIPallocBufferArray(scip, &collev, gmoN(gmo)) );
      SCIP_CALL( SCIPgetSolVals(scip, sol, probdata->nvars, probdata->vars, collev) );

#if GMOAPIVERSION < 12
      {
         SCIP_Real* lambda;
         int i;

         SCIP_CALL( SCIPallocBufferArray(scip, &lambda, gmoM(gmo)) );
         for( i = 0; i < gmoM(gmo); ++i )
            lambda[i] = gmoValNA(gmo);

         /* this also sets the gmoHobjval attribute to the level value of GAMS' objective variable */
         gmoSetSolution2(gmo, collev, lambda);

         SCIPfreeBufferArray(scip, &lambda);
      }
#else
      (void) gmoSetSolutionPrimal(gmo, collev);
#endif
      primalbound = SCIPgetPrimalbound(scip);

      SCIPfreeBufferArray(scip, &collev);

      /* if we have an MINLP, check if best solution is really feasible and try to repair and check other solutions otherwise */
      if( gmoObjNLNZ(gmo) != 0 || gmoNLNZ(gmo) != 0 )
      {
         SCIP_Bool success;

         /* check whether best solution is feasible in original problem */
         SCIP_CALL( SCIPcheckSolOrig(scip, sol, &success, FALSE, FALSE) );

         /* look at all solutions, try to repair if not feasible, and keep a feasible one with best objective value */
         if( !success )
         {
            SCIP_Real mainsoltime;
            SCIP_CLOCK* resolveclock;
            int origmaxorigsol;
            int orignlpverblevel;
            int origmaxpresolrounds;
            char origsolvetracefile[SCIP_MAXSTRLEN];
            /* SCIP_Real origfeastol; */
            SCIP_Bool resolvenlp;
            SCIP_Real** solvals;
            SCIP_Real* objvals;
            int nsols;
            int s;

            /* check whether we can run an NLP solver */
            SCIP_CALL( SCIPgetBoolParam(scip, "gams/resolvenlp", &resolvenlp) );
            if( resolvenlp && SCIPgetStage(scip) == SCIP_STAGE_SOLVING && !SCIPisNLPEnabled(scip) )
            {
               SCIPdebugMessage("NLP is disabled: cannot do resolves\n");
               resolvenlp = FALSE;
            }
            if( resolvenlp && SCIPgetNContVars(scip) == 0 )
            {
               SCIPdebugMessage("transformed SCIP problem has no continuous variables: cannot do resolves\n");
               resolvenlp = FALSE;
            }
            if( resolvenlp && SCIPgetNNlpis(scip) == 0 )
            {
               SCIPdebugMessage("no NLP solver: cannot do resolves\n");
               resolvenlp = FALSE;
            }
            if( resolvenlp && SCIPfindHeur(scip, "subnlp") == NULL )
            {
               SCIPdebugMessage("no NLP heuristic available: cannot do resolves\n");
               resolvenlp = FALSE;
            }

            /* try SCIP solutions (limited by limits/maxsol or limits/maxorigsol) */
            nsols = SCIPgetNSols(scip);
            SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nsols) );
            SCIP_CALL( SCIPallocBufferArray(scip, &objvals, nsols) );
            for( s = 0; s < nsols; ++s )
            {
               SCIP_CALL( SCIPallocBufferArray(scip, &solvals[s], gmoN(gmo)) ); /*lint !e866*/
               SCIP_CALL( SCIPgetSolVals(scip, SCIPgetSols(scip)[s], probdata->nvars, probdata->vars, solvals[s]) );
               objvals[s] = SCIPgetSolOrigObj(scip, SCIPgetSols(scip)[s]);
            }

            /* adapt some parameter values for possible resolve's */
            if( resolvenlp )
            {
               char* tmp;

               /* don't store solutions in original problem, so they don't get in a way when transforming for resolve */
               SCIP_CALL( SCIPgetIntParam(scip, "limits/maxorigsol", &origmaxorigsol) );
               SCIP_CALL( SCIPsetIntParam(scip, "limits/maxorigsol", 0) );

               SCIP_CALL( SCIPgetIntParam(scip, "heuristics/subnlp/nlpverblevel", &orignlpverblevel) );
               SCIP_CALL( SCIPsetIntParam(scip, "heuristics/subnlp/nlpverblevel", 1) );

               /* origfeastol = SCIPfeastol(scip); */
               /* SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol", origfeastol / 100.0) ); */

               SCIP_CALL( SCIPgetIntParam(scip, "heuristics/subnlp/maxpresolverounds", &origmaxpresolrounds) );
               SCIP_CALL( SCIPsetIntParam(scip, "heuristics/subnlp/maxpresolverounds", 0) );

               SCIP_CALL( SCIPgetStringParam(scip, "gams/solvetrace/file", &tmp) );
               strcpy(origsolvetracefile, tmp);
               SCIP_CALL( SCIPsetStringParam(scip, "gams/solvetrace/file", "") );
            }

            SCIP_CALL( SCIPcreateClock(scip, &resolveclock) );
            SCIP_CALL( SCIPstartClock(scip, resolveclock) );

            for( s = 0; s < nsols; ++s )
            {
               SCIPinfoMessage(scip, NULL, "Checking feasibility of solution #%0.2d with reported objective value %.15e.\n", s, objvals[s]);

               SCIP_CALL( checkAndRepairSol(scip, solvals[s], &objvals[s], resolvenlp, &success) );

               if( success )
                  break;
            }

            SCIP_CALL( SCIPstopClock(scip, resolveclock) );

            /* add time for checks and NLP resolves to reported solving time */
            mainsoltime = gmoGetHeadnTail(gmo, (int) gmoHresused);
            gmoSetHeadnTail(gmo, (int) gmoHresused, mainsoltime + SCIPgetClockTime(scip, resolveclock));

            SCIP_CALL( SCIPfreeClock(scip, &resolveclock) );

            if( success )
            {
               /* store updated solution in GMO */
#if GMOAPIVERSION < 12
               {
                  SCIP_Real* lambda;
                  int i;

                  SCIP_CALL( SCIPallocBufferArray(scip, &lambda, gmoM(gmo)) );
                  for( i = 0; i < gmoM(gmo); ++i )
                     lambda[i] = gmoValNA(gmo);

                  /* this also sets the gmoHobjval attribute to the level value of GAMS' objective variable */
                  gmoSetSolution2(gmo, solvals[s], lambda);

                  SCIPfreeBufferArray(scip, &lambda);
               }
#else
               (void) gmoSetSolutionPrimal(gmo, solvals[s]);
#endif
               /* update reevaluated objective value */
               objvals[s] = gmoGetHeadnTail(gmo, (int) gmoHobjval);
               primalbound = objvals[s];

               SCIPinfoMessage(scip, NULL, "Solution #%0.2d feasible. Reevaluated objective value = %.15e.\n", s, objvals[s]);

               SCIPinfoMessage(scip, NULL, "\nStatus update:\n");
               SCIPinfoMessage(scip, NULL, "Solving Time (sec) : %.2f\n", gmoGetHeadnTail(gmo, (int) gmoHresused));
               SCIPinfoMessage(scip, NULL, "Primal Bound       : %+.14e\n", objvals[s]);
               if( gmoGetHeadnTail(gmo, (int) gmoTmipbest) != gmoValNA(gmo) ) /*lint !e777*/
               {
                  SCIPinfoMessage(scip, NULL, "Dual Bound         : %+.14e\n", dualbound);
                  SCIPinfoMessage(scip, NULL, "Gap                : ");

                  if( SCIPisEQ(scip, objvals[s], dualbound) )
                     SCIPinfoMessage(scip, NULL, "%.2f %%\n", 0.0);
                  else if( SCIPisZero(scip, dualbound)
                     || SCIPisZero(scip, objvals[s])
                     || (dualbound == gmoValNA(gmo))  /*lint !e777*/
                     || SCIPisInfinity(scip, REALABS(objvals[s]))
                     || SCIPisInfinity(scip, REALABS(dualbound))
                     || objvals[s] * dualbound < 0.0 )
                     SCIPinfoMessage(scip, NULL, "infinite\n");
                  else
                     SCIPinfoMessage(scip, NULL, "%.2f %%\n", REALABS((objvals[s] - dualbound)/MIN(REALABS(dualbound),REALABS(objvals[s])))); /*lint !e666*/
               }
            }
            else
            {
               SCIPinfoMessage(scip, NULL, "None of %d SCIP solutions could be made feasible.\n", nsols);
            }

            /* restore original parameter values */
            if( resolvenlp )
            {
               SCIP_CALL( SCIPsetIntParam(scip, "limits/maxorigsol", origmaxorigsol) ); /*lint !e644*/
               SCIP_CALL( SCIPsetIntParam(scip, "heuristics/subnlp/nlpverblevel", orignlpverblevel) ); /*lint !e644*/
               /* SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol", origfeastol) ); */
               SCIP_CALL( SCIPsetIntParam(scip, "heuristics/subnlp/maxpresolverounds", origmaxpresolrounds) ); /*lint !e644*/
               SCIP_CALL( SCIPsetStringParam(scip, "gams/solvetrace/file", origsolvetracefile) );
            }

            for( s = 0; s < nsols; ++s )
            {
               SCIPfreeBufferArray(scip, &solvals[s]);
            }
            SCIPfreeBufferArray(scip, &solvals);
            SCIPfreeBufferArray(scip, &objvals);
         }

         /* update model status */
         if( !success )
         {
            /* couldn't get a feasible solution, report intermediate infeasible */
            gmoModelStatSet(gmo, (int) gmoModelStat_InfeasibleIntermed);
         }
         else if( !SCIPisEQ(scip, primalbound, dualbound) )
         {
            /* feasible, but gap not closed, so only local optimum */
            gmoModelStatSet(gmo, gmoNDisc(gmo) ? (int) gmoModelStat_Integer : (int) gmoModelStat_Feasible);
         }
         else
         {
            /* feasible and gap closed, so solved globally */
            gmoModelStatSet(gmo, (int) gmoModelStat_OptimalGlobal);
         }

      }
   }

   if( gmoModelType(gmo) == (int) gmoProc_cns )
      switch( gmoModelStat(gmo) )
      {
         case gmoModelStat_OptimalGlobal:
         case gmoModelStat_OptimalLocal:
         case gmoModelStat_Feasible:
         case gmoModelStat_Integer:
            gmoModelStatSet(gmo, (int) gmoModelStat_Solved);
      } /*lint !e744*/

   return SCIP_OKAY;
}

/*
 * Callback methods of reader
 */


/** copy method for reader plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_READERCOPY(readerCopyGmo)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of gmo reader not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
 
   return SCIP_OKAY;
}
#else
#define readerCopyGmo NULL
#endif

/** destructor of reader to free user data (called when SCIP is exiting) */
#if 1
static
SCIP_DECL_READERFREE(readerFreeGmo)
{
   SCIP_READERDATA* readerdata;

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

#if 0 /* TODO should do this if we created gmo */
   if( readerdata->gmo != NULL )
   {
      /* write solution file */
      gmoUnloadSolutionLegacy(readerdata->gmo);

      gmoFree(&readerdata->gmo);
      gevFree(&readerdata->gev);

      gmoLibraryUnload();
      gevLibraryUnload();
   }
#endif

   SCIPfreeMemory(scip, &readerdata);

   return SCIP_OKAY;
}
#else
#define readerFreeGmo NULL
#endif


/** problem reading method of reader */
#if 1
static
SCIP_DECL_READERREAD(readerReadGmo)
{
   SCIP_READERDATA* readerdata;
   gmoHandle_t gmo;
   gevHandle_t gev;
   char buffer[1024];
   
   *result = SCIP_DIDNOTRUN;

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   if( readerdata->gmo == NULL )
   {
      /* initialize GMO and GEV libraries */
      if( !gmoCreate(&readerdata->gmo, buffer, sizeof(buffer)) || !gevCreate(&readerdata->gev, buffer, sizeof(buffer)) )
      {
         SCIPerrorMessage(buffer);
         return SCIP_ERROR;
      }

      gmo = readerdata->gmo;
      gev = readerdata->gev;

      /* load control file */
      if( gevInitEnvironmentLegacy(gev, filename) )
      {
         SCIPerrorMessage("Could not load control file %s\n", filename);
         (void) gmoFree(&gmo);
         (void) gevFree(&gev);
         return SCIP_READERROR;
      }

      if( gmoRegisterEnvironment(gmo, gev, buffer) )
      {
         SCIPerrorMessage("Error registering GAMS Environment: %s\n", buffer);
         (void) gmoFree(&gmo);
         (void) gevFree(&gev);
         return SCIP_ERROR;
      }

      if( gmoLoadDataLegacy(gmo, buffer) )
      {
         SCIPerrorMessage("Could not load model data.\n");
         (void) gmoFree(&gmo);
         (void) gevFree(&gev);
         return SCIP_READERROR;
      }
   }

   SCIP_CALL( SCIPcreateProblemReaderGmo(scip, readerdata->gmo, readerdata->indicatorfile, readerdata->mipstart) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}
#else
#define readerReadGmo NULL
#endif


#if 0
/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteGmo)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of gmo reader not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define readerWriteGmo NULL
#endif

/*
 * Constructs SCIP problem from the one in GMO.
 */

#define DIALOG_READGAMS_NAME                 "readgams"
#define DIALOG_READGAMS_DESC                 "initializes SCIP problem to the one stored in a GAMS modeling object"
#define DIALOG_READGAMS_ISSUBMENU            FALSE

/** copy method for dialog plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_DIALOGCOPY(dialogCopyReadGams)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of ReadGams dialog not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dialogCopyReadGams NULL
#endif

/** destructor of dialog to free user data (called when the dialog is not captured anymore) */
#if 0
static
SCIP_DECL_DIALOGFREE(dialogFreeReadGams)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of ReadGams dialog not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dialogFreeReadGams NULL
#endif

/** description output method of dialog */
#if 0
static
SCIP_DECL_DIALOGDESC(dialogDescReadGams)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of ReadGams dialog not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dialogDescReadGams NULL
#endif


/** execution method of dialog */
static
SCIP_DECL_DIALOGEXEC(dialogExecReadGams)
{  /*lint --e{715}*/
   SCIP_READERDATA* readerdata;

   assert(dialoghdlr != NULL);
   assert(dialog != NULL);
   assert(scip != NULL);

   readerdata = (SCIP_READERDATA*)SCIPdialogGetData(dialog);
   assert(readerdata != NULL);

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   if( readerdata->gmo == NULL )
   {
      SCIPerrorMessage("GMO not initialized, cannot setup GAMS problem\n");
   }
   else
   {
      /* free previous problem and solution data, if existing */
      SCIP_CALL( SCIPfreeProb(scip) );

      SCIP_CALL( SCIPcreateProblemReaderGmo(scip, readerdata->gmo, readerdata->indicatorfile, readerdata->mipstart) );
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
      "\noriginal problem has %d variables (%d bin, %d int, %d cont) and %d constraints\n",
      SCIPgetNOrigVars(scip), SCIPgetNOrigBinVars(scip), SCIPgetNOrigIntVars(scip), SCIPgetNOrigContVars(scip),
      SCIPgetNConss(scip));

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}


/*
 * Writing solution information to GMO dialog
 */

#define DIALOG_WRITEGAMSSOL_NAME             "gamssol"
#define DIALOG_WRITEGAMSSOL_DESC             "writes solution information into GAMS Modeling Object"
#define DIALOG_WRITEGAMSSOL_ISSUBMENU        FALSE

/** copy method for dialog plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_DIALOGCOPY(dialogCopyWriteGamsSol)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of WriteGamsSol dialog not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dialogCopyWriteGamsSol NULL
#endif

/** destructor of dialog to free user data (called when the dialog is not captured anymore) */
#if 0
static
SCIP_DECL_DIALOGFREE(dialogFreeWriteGamsSol)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of WriteGamsSol dialog not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dialogFreeWriteGamsSol NULL
#endif

/** description output method of dialog */
#if 0
static
SCIP_DECL_DIALOGDESC(dialogDescWriteGamsSol)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of WriteGamsSol dialog not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dialogDescWriteGamsSol NULL
#endif


/** execution method of dialog */
static
SCIP_DECL_DIALOGEXEC(dialogExecWriteGamsSol)
{  /*lint --e{715}*/
   SCIP_READERDATA* readerdata;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   readerdata = (SCIP_READERDATA*) SCIPdialogGetData(dialog);
   assert(readerdata != NULL);

   SCIP_CALL( writeGmoSolution(scip, readerdata->gmo) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}


/*
 * Loading GAMS and user option file dialog
 */

#define DIALOG_SETTINGSLOADGAMS_NAME         "loadgams"
#define DIALOG_SETTINGSLOADGAMS_DESC         "loads GAMS settings and SCIP option file specified in GAMS model"
#define DIALOG_SETTINGSLOADGAMS_ISSUBMENU    FALSE

/** copy method for dialog plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_DIALOGCOPY(dialogCopySettingsLoadGams)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of SettingsLoadGams dialog not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dialogCopySettingsLoadGams NULL
#endif

/** destructor of dialog to free user data (called when the dialog is not captured anymore) */
#if 0
static
SCIP_DECL_DIALOGFREE(dialogFreeSettingsLoadGams)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of SettingsLoadGams dialog not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dialogFreeSettingsLoadGams NULL
#endif

/** description output method of dialog */
#if 0
static
SCIP_DECL_DIALOGDESC(dialogDescSettingsLoadGams)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of SettingsLoadGams dialog not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dialogDescSettingsLoadGams NULL
#endif


/** execution method of dialog */
static
SCIP_DECL_DIALOGEXEC(dialogExecSettingsLoadGams)
{  /*lint --e{715}*/
   assert(scip       != NULL);
   assert(dialog     != NULL);
   assert(dialoghdlr != NULL);

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIP_CALL( SCIPreadParamsReaderGmo(scip) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/*
 * reader specific interface methods
 */

/** includes the gmo file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderGmo(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_DIALOG* dialog;
   SCIP_DIALOG* parentdialog;

   /* create gmo reader data */
   SCIP_CALL( SCIPallocMemory(scip, &readerdata) );
   BMSclearMemory(readerdata);
   
   /* include gmo reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerCopyGmo,
         readerFreeGmo, readerReadGmo, readerWriteGmo,
         readerdata) );

   SCIP_CALL( SCIPaddStringParam(scip, "gams/dumpsolutions",
      "name of solutions index gdx file for writing all alternate solutions",
      NULL, FALSE, "", NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, "gams/dumpsolutionsmerged",
      "name of gdx file for writing all alternate solutions into a single file",
      NULL, FALSE, "", NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "gams/resolvenlp",
      "whether to resolve MINLP with fixed discrete variables if best solution violates some constraints",
      NULL, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "gams/mipstart",
      "how to handle initial variable levels",
      &readerdata->mipstart, FALSE, 2, 0, 4, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, "gams/indicatorfile",
      "name of GAMS options file that contains definitions on indicators",
      &readerdata->indicatorfile, FALSE, "", NULL, NULL) );

   /* get parent dialog "write" */
   if( SCIPdialogFindEntry(SCIPgetRootDialog(scip), "write", &parentdialog) != 1 )
   {
      SCIPerrorMessage("sub menu \"write\" not found\n");
      return SCIP_PLUGINNOTFOUND;
   }
   assert(parentdialog != NULL);

   /* create, include, and release dialog */
   if( !SCIPdialogHasEntry(SCIPgetRootDialog(scip), DIALOG_READGAMS_NAME) )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            dialogCopyReadGams, dialogExecReadGams, dialogDescReadGams, dialogFreeReadGams,
            DIALOG_READGAMS_NAME, DIALOG_READGAMS_DESC, DIALOG_READGAMS_ISSUBMENU, (SCIP_DIALOGDATA*)readerdata) );
      SCIP_CALL( SCIPaddDialogEntry(scip, SCIPgetRootDialog(scip), dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }


   /* create, include, and release dialog */
   if( !SCIPdialogHasEntry(parentdialog, DIALOG_WRITEGAMSSOL_NAME) )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            dialogCopyWriteGamsSol, dialogExecWriteGamsSol, dialogDescWriteGamsSol, dialogFreeWriteGamsSol,
            DIALOG_WRITEGAMSSOL_NAME, DIALOG_WRITEGAMSSOL_DESC, DIALOG_WRITEGAMSSOL_ISSUBMENU, (SCIP_DIALOGDATA*)readerdata) );
      SCIP_CALL( SCIPaddDialogEntry(scip, parentdialog, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }


   /* get parent dialog "set" */
   if( SCIPdialogFindEntry(SCIPgetRootDialog(scip), "set", &parentdialog) != 1 )
   {
      SCIPerrorMessage("sub menu \"set\" not found\n");
      return SCIP_PLUGINNOTFOUND;
   }
   assert(parentdialog != NULL);

   /* create, include, and release dialog */
   if( !SCIPdialogHasEntry(parentdialog, DIALOG_SETTINGSLOADGAMS_NAME) )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            dialogCopySettingsLoadGams, dialogExecSettingsLoadGams, dialogDescSettingsLoadGams, dialogFreeSettingsLoadGams,
            DIALOG_SETTINGSLOADGAMS_NAME, DIALOG_SETTINGSLOADGAMS_DESC, DIALOG_SETTINGSLOADGAMS_ISSUBMENU, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, parentdialog, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set gams */
   if( !SCIPdialogHasEntry(parentdialog, "gams") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL, SCIPdialogExecMenu, NULL, NULL,
            "gams", "change parameters for GAMS interface", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, parentdialog, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   return SCIP_OKAY;
}

/** sets the GMO object to use in reader
 * If GMO is set in reader, then reader does not read from file when executed, but sets up problem from GMO
 */
void SCIPsetGMOReaderGmo(
   SCIP*                 scip,               /**< SCIP data structure */
   gmoRec_t*             gmo                 /**< GMO object, or NULL to reset to default behaviour */
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;

   assert(scip != NULL);
   assert(SCIPgetStage(scip) == SCIP_STAGE_INIT);

   reader = SCIPfindReader(scip, READER_NAME);
   assert(reader != NULL);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   readerdata->gmo = gmo;
   readerdata->gev = gmo != NULL ? (gevHandle_t)gmoEnvironment(readerdata->gmo) : NULL;
}

/** passes GAMS options to SCIP and initiates reading of user options file, if given in GMO */
SCIP_RETCODE SCIPreadParamsReaderGmo(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;
   gmoHandle_t gmo;
   gevHandle_t gev;

   assert(scip != NULL);

   reader = SCIPfindReader(scip, READER_NAME);
   assert(reader != NULL);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   assert(readerdata->gmo != NULL);
   assert(readerdata->gev != NULL);

   gmo = readerdata->gmo;
   gev = readerdata->gev;

   if( gevGetIntOpt(gev, gevNodeLim) > 0 )
   {
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", (long long)gevGetIntOpt(gev, gevNodeLim)) );
   }
   SCIP_CALL( SCIPsetRealParam(scip, "limits/time",   gevGetDblOpt(gev, gevResLim)) );
   SCIP_CALL( SCIPsetRealParam(scip, "limits/gap",    gevGetDblOpt(gev, gevOptCR)) );
   if( !SCIPisInfinity(scip, gevGetDblOpt(gev, gevOptCA)) )
   {
      SCIP_CALL( SCIPsetRealParam(scip, "limits/absgap", gevGetDblOpt(gev, gevOptCA)) );
   }
   else
   {
      SCIPwarningMessage(scip, "Value for optca = %g >= value for infinity. Setting solution limit to 1 instead.\n", gevGetDblOpt(gev, gevOptCA));
      SCIP_CALL( SCIPsetIntParam(scip, "limits/solutions", 1) );
   }
   if( gevGetDblOpt(gev, gevWorkSpace) > 0.0 )
   {
      SCIP_CALL( SCIPsetRealParam(scip, "limits/memory", gevGetDblOpt(gev, gevWorkSpace)) );
   }
   SCIP_CALL( SCIPsetIntParam(scip, "lp/threads", gevThreads(gev)) );
#if 0
   SCIP_CALL( SCIPsetIntParam(scip, "timing/clocktype", 2) ); /* wallclock time */
#endif

   /* if log is not kept, then can also set SCIP verblevel to 0 */
   if( gevGetIntOpt(gev, gevLogOption) == 0 )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) );
   }

#ifdef _WIN32
   if( !gevGetIntOpt(gev, gevIDEFlag) )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "display/width", 80) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/lpavgiterations/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/maxdepth/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/time/active", 2) );
   }
#endif

   /* enable column on number of branching on nonlinear variables, if any */
   if( gmoNLNZ(gmo) > 0 || (gmoObjStyle(gmo) == (int) gmoObjType_Fun && gmoObjNLNZ(gmo) > 0) )
   {
      /* enable column on number of branching on continuous variables */
      SCIP_CALL( SCIPsetIntParam(scip, "display/nexternbranchcands/active", 2) );
   }
   /* make sure column on number of branching on fractional variables is shown, if any */
   if( gmoNDisc(gmo) > 0 )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "display/nfrac/active", 2) );
   }

   /* don't print reason why start solution is infeasible, per default */
   SCIP_CALL( SCIPsetBoolParam(scip, "misc/printreason", FALSE) );

   if( gmoOptFile(gmo) > 0 )
   {
      char optfilename[1024];
      SCIP_RETCODE ret;

      (void) gmoNameOptFile(gmo, optfilename);
      SCIPinfoMessage(scip, NULL, "\nreading option file %s\n", optfilename);
      ret = SCIPreadParams(scip, optfilename);
      if( ret != SCIP_OKAY )
      {
         SCIPwarningMessage(scip, "Reading of optionfile %s failed with SCIP return code <%d>!\n", optfilename, ret);
      }
   }

   return SCIP_OKAY;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_nl.c
 * @brief  AMPL .nl file reader
 * @author Stefan Vigerske
 *
 * The code (including some comments) for this reader is based on OSnl2osil.cpp,
 * the nl-reader of the Optimization Services Project (https://projects.coin-or.org/OS).
 *
 * The code for SOS reading is based on the AMPL/Bonmin interface (https://projects.coin-or.org/Bonmin).
 *
 * For an incomplete documentation on how to hook up a solver to AMPL, see http://www.ampl.com/REFS/HOOKING/index.html.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "reader_nl.h"
#include "scip/cons_linear.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"

/* enable the following define to enable recognition of curvature suffix
 * functions of constraints for which this suffix is set are then handled via user-expressions
 */
/* #define CHECKCURVSUFFIX */

/* we need the ABS define from ASL later */
#undef ABS

/* ASL includes */
#include "nlp.h"
/* #include "getstub.h" */
#include "opcode.hd"
#include "asl.h"

#undef filename

#define READER_NAME             "nlreader"
#define READER_DESC             "AMPL .nl file reader"
#define READER_EXTENSION        "nl"


/*
 * Data structures
 */

struct cgrad;

/** problem data */
struct SCIP_ProbData
{
   ASL*                  asl;                /**< ASL data structure */
   SCIP_VAR**            vars;               /**< variables */
   int                   nvars;              /**< number of variables */
   SCIP_Real*            fullx;              /**< scratch memory to store full point during evaluation of AMPL userexpr */
};

/** user expression data */
struct SCIP_UserExprData
{
   SCIP_PROBDATA*        probdata;           /**< problem data structure */
   int                   considx;            /**< index for constraint */
   SCIP_EXPRCURV         curvature;          /**< curvature information */
};

/*
 * Local methods
 */

static
SCIP_DECL_USEREXPREVAL(SCIPuserexprEvalAmpl)
{
   ASL* asl;
   fint nerror = 0;
   cgrad* cg;
   int i;

   assert(data != NULL);
   assert(data->probdata->fullx != NULL);
   assert(funcvalue != NULL);

   asl = data->probdata->asl;

   xunknown();

   i = 0;
   for( cg = Cgrad[data->considx]; cg != NULL; cg = cg->next )
      data->probdata->fullx[cg->varno] = argvals[i++];

   /* function value */
   *funcvalue = conival(data->considx, data->probdata->fullx, &nerror);
   if( nerror != 0 )
   {
      *funcvalue = SCIP_INVALID;
      return SCIP_OKAY;
   }

   /* gradient */
   if( gradient != NULL )
   {
      congrd(data->considx, data->probdata->fullx, gradient, &nerror);

      if( nerror != 0 )
      {
         *funcvalue = SCIP_INVALID;
         return SCIP_OKAY;
      }
   }

   /* TODO implement Hessian */
   assert(hessian == NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_USEREXPRINTEVAL(SCIPuserexprIntEvalAmpl)
{
   ASL* asl;
   cgrad* cg;
   int i;
   int p;
   SCIP_Real extrval;

   assert(data != NULL);
   assert(data->probdata->fullx != NULL);
   assert(funcvalue != NULL);

   asl = data->probdata->asl;

   SCIPintervalSetEntire(infinity, funcvalue);

   if( gradient != NULL )
   {
      for( i = 0; i < nargs; ++i )
         SCIPintervalSetEntire(infinity, &gradient[i]);
   }

   if( hessian != NULL )
   {
      for( i = 0; i < nargs * nargs; ++i )
         SCIPintervalSetEntire(infinity, &hessian[i]);
      for( i = 0; i < nargs; ++i )
      {
         if( data->curvature == SCIP_EXPRCURV_CONVEX )
            SCIPintervalSetBounds(&gradient[i*nargs+i], 0.0, infinity);
         else
            SCIPintervalSetBounds(&gradient[i*nargs+i], -infinity, 0.0);
      }
   }

   if( nargs > 10 )
      return SCIP_OKAY;

   /* we evaluate all corner points and look for
    * a maximal value for a convex function and
    * a minimal value for a concave function
    */
   if( data->curvature == SCIP_EXPRCURV_CONVEX )
      extrval = -infinity;
   else
      extrval =  infinity;

   for( p = 0; p < (1<<nargs); ++p )
   {
      SCIP_Real cornerval;
      fint nerror = 0;
      int q = p;
      for( cg = Cgrad[data->considx]; cg != NULL; cg = cg->next )
      {
         /* if bound at infinity, then we cannot evaluate, so stop */
         if( argvals->inf <= -infinity || argvals->sup >= infinity )
            break;
         data->probdata->fullx[cg->varno] = (q % 2 ? argvals->inf : argvals->sup);
         q = q >> 1;
      }
      /* if found variable at infinity, then stop */
      if( cg != NULL )
         break;

      cornerval = conival(data->considx, data->probdata->fullx, &nerror);

      /* if cannot evaluate, then stop */
      if( nerror != 0 )
         break;

      if( data->curvature == SCIP_EXPRCURV_CONVEX )
      {
         if( cornerval > extrval )
            extrval = cornerval;
      }
      else
      {
         if( cornerval < extrval )
            extrval = cornerval;
      }
   }

   if( p == (1<<nargs) )
   {
      /* all points could be evaluated, so can update funcvalue interval
       * as we have used floating point arithmetic, add a small safety margin
       */
      if( data->curvature == SCIP_EXPRCURV_CONVEX )
         funcvalue->sup = extrval + 1e-9 * fabs(extrval);
      else
         funcvalue->inf = extrval - 1e-9 * fabs(extrval);
   }

   /* TODO find other extreme value by solving a convex box-constrained optimization problem
    * in 1D, could also use bisection */

   return SCIP_OKAY;
}

static
SCIP_DECL_USEREXPRCURV(SCIPuserexprCurvAmpl)
{
   *result = data->curvature;
   return SCIP_OKAY;
}

static
SCIP_DECL_USEREXPRCOPYDATA(SCIPuserexprCopyAmpl)
{
   SCIP_ALLOC( BMSduplicateBlockMemory(blkmem, datatarget, datasource) );

   return SCIP_OKAY;
}

static
SCIP_DECL_USEREXPRFREEDATA(SCIPuserexprFreeAmpl)
{
   BMSfreeBlockMemory(blkmem, &data);
}

#define CHECKSUFFIX(sufname, idx, defaultval) \
   ( ( suf_ ## sufname != NULL && (suf_ ## sufname)->u.i[idx] > 0 ) ? \
     ((suf_ ## sufname)->u.i[idx] == 1) : \
     defaultval )

static
SCIP_RETCODE setupVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   ASL* asl;
   int lower;
   int upper;
   int i;
   SufDesc* suf_initial;
   SufDesc* suf_removable;

   assert(scip != NULL);
   assert(probdata != NULL);
   assert(probdata->asl != NULL);
   assert(probdata->vars == NULL);

   asl = probdata->asl;

   suf_initial = suf_get("initial", ASL_Sufkind_var);
   if( suf_initial != NULL && suf_initial->u.i == NULL )
      suf_initial = NULL;

   suf_removable = suf_get("removable", ASL_Sufkind_var);
   if( suf_removable != NULL && suf_removable->u.i == NULL )
      suf_removable = NULL;

   SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->vars, n_var) );
   probdata->nvars = n_var;

   /* first the nonlinear variables */
   /* welcome to the world of the ASL API */

   lower = 0;
   upper = nlvb - nlvbi;
   for( i = lower; i < upper; ++i ) /* continuous and in an objective and in a constraint */
   {
      SCIP_CALL( SCIPcreateVar(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_CONTINUOUS,
         CHECKSUFFIX(initial, i, TRUE), CHECKSUFFIX(removable, i, FALSE), NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }

   lower = nlvb - nlvbi;
   upper = nlvb;
   for( i = lower; i < upper; ++i ) /* integer and in an objective and in a constraint */
   {
       SCIP_CALL( SCIPcreateVar(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_INTEGER,
          CHECKSUFFIX(initial, i, TRUE), CHECKSUFFIX(removable, i, FALSE), NULL, NULL, NULL, NULL, NULL) );
       SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }

   lower = nlvb;
   /* upper = nlvb + (nlvc - (nlvb + nlvci)); */
   upper = nlvc - nlvci;
   for( i = lower; i < upper; ++i ) /* continuous and just in constraints */
   {
      SCIP_CALL( SCIPcreateVar(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_CONTINUOUS,
         CHECKSUFFIX(initial, i, TRUE), CHECKSUFFIX(removable, i, FALSE), NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }

   lower = nlvc - nlvci;
   upper = nlvc;
   for( i = lower; i < upper; ++i ) /* integer and just in constraints */
   {
      SCIP_CALL( SCIPcreateVar(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_INTEGER,
         CHECKSUFFIX(initial, i, TRUE), CHECKSUFFIX(removable, i, FALSE), NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }

   lower = nlvc;
   /* upper = nlvc + ( nlvo - (nlvc + nlvoi) ); */
   upper = nlvo - nlvoi;
   for( i = lower; i < upper; ++i) /* continuous and just in objectives */
   {
      SCIP_CALL( SCIPcreateVar(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_CONTINUOUS,
         CHECKSUFFIX(initial, i, TRUE), CHECKSUFFIX(removable, i, FALSE), NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }

   lower = nlvo - nlvoi;
   upper = nlvo ;
   for( i = lower; i < upper; ++i ) /* integer and just in objectives */
   {
      SCIP_CALL( SCIPcreateVar(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_INTEGER,
         CHECKSUFFIX(initial, i, TRUE), CHECKSUFFIX(removable, i, FALSE), NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }

   /* now the other variables */

   lower = MAX(nlvc, nlvo);
   upper = MAX(nlvc, nlvo) + nwv;
   for( i = lower; i < upper; ++i ) /* linear arc variables */
   {
      SCIP_CALL( SCIPcreateVar(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_CONTINUOUS,
         CHECKSUFFIX(initial, i, TRUE), CHECKSUFFIX(removable, i, FALSE), NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }


   lower = MAX(nlvc, nlvo) + nwv;
   /* upper = MAX(nlvc, nlvo) + nwv + (n_var - (MAX(nlvc, nlvo) + niv + nbv + nwv) ); */
   upper = n_var -  niv - nbv;
   for( i = lower; i < upper; ++i ) /* other linear variables */
   {
      SCIP_CALL( SCIPcreateVar(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_CONTINUOUS,
         CHECKSUFFIX(initial, i, TRUE), CHECKSUFFIX(removable, i, FALSE), NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }

   lower = n_var -  niv - nbv;
   upper = n_var -  niv ;
   for( i = lower; i < upper; ++i ) /* linear and binary */
   {
      SCIP_CALL( SCIPcreateVar(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_BINARY,
         CHECKSUFFIX(initial, i, TRUE), CHECKSUFFIX(removable, i, FALSE), NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }

   lower = n_var -  niv;
   upper = n_var;
   for( i = lower; i < upper; ++i ) /* linear and integer */
   {
      SCIP_CALL( SCIPcreateVar(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_INTEGER,
         CHECKSUFFIX(initial, i, TRUE), CHECKSUFFIX(removable, i, FALSE), NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }

   /* end of variables -- thank goodness!!! */

   return SCIP_OKAY;
}


/** transforms AMPL expression tree into SCIP expression */
static
SCIP_RETCODE walkExpression(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           scipexpr,           /**< buffer to store pointer to created expression */
   ASL*                  asl,                /**< AMPL solver library structure */
   expr*                 amplexpr,           /**< root node of AMPL expression */
   int*                  exprvaridx,         /**< array with index of problem variables in expression graph */
   int*                  nexprvars,          /**< number of variables in currently processed expression so far */
   int                   nvars,              /**< total number of variables in problem (and length of exprvaridx array) */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occurred */
   )
{
   int opnum;
   int i;

   assert(scip != NULL);
   assert(scipexpr != NULL);
   assert(amplexpr != NULL);
   assert(exprvaridx != NULL || nvars == 0);
   assert(nexprvars != NULL);
   assert(doingfine != NULL);

   opnum = Intcast amplexpr->op;

   switch( opnum )
   {
      case OPNUM:
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), scipexpr, SCIP_EXPR_CONST, ((expr_n*)amplexpr)->v) );
         break;
      }

      case OPVARVAL:
      {
         int varidx;

         varidx = (expr_v*)amplexpr - ((ASL_fg*)asl)->I.var_e_;  /*lint !e826*/

         /* treat the common expression or defined variables */
         if( varidx >= n_var )
         {
            /* process common expression, see also http://www.gerad.ca/~orban/drampl/def-vars.html */
            expr* nlpart;
            linpart* L;
            int n_lin;

            if( varidx - n_var < ncom0 )
            {
               struct cexp* common;

               common = ((const ASL_fg *) asl)->I.cexps_ + varidx - n_var;  /*lint !e826*/

               nlpart = common->e;
               L = common->L;
               n_lin = common->nlin;
            }
            else
            {
               struct cexp1* common;

               common = ((const ASL_fg *) asl)->I.cexps1_ + varidx - n_var - ncom0;  /*lint !e826*/

               nlpart = common->e;
               L = common->L;
               n_lin = common->nlin;
            }

            /* get nonlinear expression corresponding to defined variable */
            SCIP_CALL( walkExpression(scip, scipexpr, asl, nlpart, exprvaridx, nexprvars, nvars, doingfine) );
            if( !*doingfine )
               break;

            if( n_lin > 0 )
            {
               /* add linear part of defined variable */
               SCIP_Real* coefs;
               SCIP_EXPR** children;

               SCIP_CALL( SCIPallocBufferArray(scip, &coefs, n_lin+1) );
               SCIP_CALL( SCIPallocBufferArray(scip, &children, n_lin+1) );
               for( i = 0; i < n_lin; ++i )
               {
                  varidx = ((size_t) (L[i].v.rp) - (size_t) ((ASL_fg*)asl)->I.var_e_) / sizeof (expr_v);  /*lint !e826 !e713*/

                  /* assign index to variable, if we see it the first time */
                  if( exprvaridx[varidx] == -1 )  /*lint !e613*/
                  {
                     exprvaridx[varidx] = *nexprvars;  /*lint !e613*/
                     ++*nexprvars;
                  }

                  SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &children[i], SCIP_EXPR_VARIDX, exprvaridx[varidx]) );  /*lint !e613 */

                  coefs[i] = L[i].fac;
               }
               coefs[i] = 1.0;
               children[i] = *scipexpr;

               SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), scipexpr, n_lin+1, children, coefs, 0.0) );

               SCIPfreeBufferArray(scip, &coefs);
               SCIPfreeBufferArray(scip, &children);
            }
         }
         else
         {
            varidx = amplexpr->a;
            assert(varidx <= nvars);

            /* assign index to variable, if we see it the first time */
            if( exprvaridx[varidx] == -1 )  /*lint !e613*/
            {
               exprvaridx[varidx] = *nexprvars;  /*lint !e613*/
               ++*nexprvars;
            }

            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), scipexpr, SCIP_EXPR_VARIDX, exprvaridx[varidx]) );  /*lint !e613*/
         }

         break;
      }

      case OPPLUS:
      case OPMINUS:
      case OPMULT:
      case OPDIV:
      {
         SCIP_EXPR* child1;
         SCIP_EXPR* child2;
         SCIP_EXPROP operand;

         SCIP_CALL( walkExpression(scip, &child1, asl, amplexpr->L.e, exprvaridx, nexprvars, nvars, doingfine) );
         if( !*doingfine )
            break;

         SCIP_CALL( walkExpression(scip, &child2, asl, amplexpr->R.e, exprvaridx, nexprvars, nvars, doingfine) );
         if( !*doingfine )
         {
            SCIPexprFreeDeep(SCIPblkmem(scip), &child1);
            break;
         }

         switch( opnum )
         {
            case OPPLUS:
               operand = SCIP_EXPR_PLUS;
               break;
            case OPMINUS:
               operand = SCIP_EXPR_MINUS;
               break;
            case OPMULT:
               operand = SCIP_EXPR_MUL;
               break;
            case OPDIV:
               operand = SCIP_EXPR_DIV;
               break;
         }  /*lint !e744*/
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), scipexpr, operand, child1, child2) );  /*lint !e644*/

         break;
      }

      case OPSUMLIST:
      {
         SCIP_EXPR** children;
         int nchildren;

         nchildren = amplexpr->R.ep - amplexpr->L.ep;

         SCIP_CALL( SCIPallocBufferArray(scip, &children, nchildren) );

         for( i = 0; i < nchildren; ++i )
         {
            SCIP_CALL( walkExpression(scip, &children[i], asl, amplexpr->L.ep[i], exprvaridx, nexprvars, nvars, doingfine) );
            if( !*doingfine )
               break;
         }
         if( !*doingfine )
         {
            for( --i; i >= 0; --i )
               SCIPexprFreeDeep(SCIPblkmem(scip), &children[i]);

            SCIPfreeBufferArray(scip, &children);
            break;
         }

         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), scipexpr, SCIP_EXPR_SUM, nchildren, children) );

         SCIPfreeBufferArray(scip, &children);

         break;
      }

      case OPUMINUS: /* negate */
      {
         SCIP_Real minusone;

         SCIP_CALL( walkExpression(scip, scipexpr, asl, amplexpr->L.e, exprvaridx, nexprvars, nvars, doingfine) );
         if( !*doingfine )
            break;

         minusone = -1.0;
         SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), scipexpr, 1, scipexpr, &minusone, 0.0) );

         break;
      }

      case OPPOW: /* general power expr^expr */
      {
         SCIP_EXPR* child1;
         SCIP_EXPR* child2;

         SCIP_CALL( walkExpression(scip, &child1, asl, amplexpr->L.e, exprvaridx, nexprvars, nvars, doingfine) );
         if( !*doingfine )
            break;

         SCIP_CALL( walkExpression(scip, &child2, asl, amplexpr->R.e, exprvaridx, nexprvars, nvars, doingfine) );
         if( !*doingfine )
         {
            SCIPexprFreeDeep(SCIPblkmem(scip), &child1);
            break;
         }

         if( SCIPexprGetOperator(child2) == SCIP_EXPR_CONST )
         {
            /* expr^number is intpower or realpower */
            if( SCIPisIntegral(scip, SCIPexprGetOpReal(child2)) )
            {
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), scipexpr, SCIP_EXPR_INTPOWER, child1, (int)SCIPround(scip, SCIPexprGetOpReal(child2))) );
            }
            else
            {
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), scipexpr, SCIP_EXPR_REALPOWER, child1, SCIPexprGetOpReal(child2)) );
            }
         }
         else if( SCIPexprGetOperator(child1) == SCIP_EXPR_CONST )
         {
            /* number^arg2 is exp(arg2 * ln(number)) */
            if( SCIPexprGetOpReal(child1) < 0.0 )
            {
               SCIPerrorMessage("Negative base in OPPOW expression with nonconstant exponent not allowed in nonlinear expression.\n");
               SCIPexprFreeDeep(SCIPblkmem(scip), &child1);
               SCIPexprFreeDeep(SCIPblkmem(scip), &child2);
               *doingfine = FALSE;
               break;
            }
            else
            {
               SCIP_EXPR* tmp;

               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp, SCIP_EXPR_CONST, log(SCIPexprGetOpReal(child1))) );
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp, SCIP_EXPR_MUL, tmp, child2) );
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), scipexpr, SCIP_EXPR_EXP, tmp) );
            }
         }
         else
         {
            /* arg1^arg2 is exp(arg2 * ln(arg1)) */
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &child1, SCIP_EXPR_LOG, child1) );
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &child2, SCIP_EXPR_MUL, child1, child2) );
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), scipexpr, SCIP_EXPR_EXP, child2) );
         }

         break;
      }

      case OP1POW: /* power expr^number */
      {
         SCIP_CALL( walkExpression(scip, scipexpr, asl, amplexpr->L.e, exprvaridx, nexprvars, nvars, doingfine) );
         if( !*doingfine )
            break;

         /* expr^number is intpower or realpower */
         if( SCIPisIntegral(scip, amplexpr->R.en->v) )
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), scipexpr, SCIP_EXPR_INTPOWER, *scipexpr, (int)SCIPround(scip, amplexpr->R.en->v)) );
         }
         else
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), scipexpr, SCIP_EXPR_REALPOWER, *scipexpr, amplexpr->R.en->v) );
         }

         break;
      }

      case OPCPOW: /* power number^expr */
      {
         SCIP_EXPR* tmp;

         /* number^expr is exp(expr * ln(number)) */
         if( amplexpr->L.en->v < 0.0 )
         {
            SCIPerrorMessage("Negative base in OPCPOW expression with nonconstant exponent not allowed in nonlinear expression.\n");
            *doingfine = FALSE;
            break;
         }

         SCIP_CALL( walkExpression(scip, scipexpr, asl, amplexpr->R.e, exprvaridx, nexprvars, nvars, doingfine) );
         if( !*doingfine )
            break;

         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp, SCIP_EXPR_CONST, log(amplexpr->L.en->v)) );
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp, SCIP_EXPR_MUL, tmp, *scipexpr) );
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), scipexpr, SCIP_EXPR_EXP, tmp) );

         break;
      }

      case OP2POW: /* square */
      case OP_log:
      case OP_sqrt:
      case OP_exp:
      case ABS:
      {
         SCIP_EXPROP operand;

         SCIP_CALL( walkExpression(scip, scipexpr, asl, amplexpr->L.e, exprvaridx, nexprvars, nvars, doingfine) );
         if( !*doingfine )
            break;

         switch( opnum )
         {
            case OP2POW:
               operand = SCIP_EXPR_SQUARE;
               break;
            case OP_log:
               operand = SCIP_EXPR_LOG;
               break;
            case OP_sqrt:
               operand = SCIP_EXPR_SQRT;
               break;
            case OP_exp:
               operand = SCIP_EXPR_EXP;
               break;
            case ABS:
               operand = SCIP_EXPR_ABS;
               break;
            default:
               return SCIP_ERROR;
         }

         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), scipexpr, operand, *scipexpr) );

         break;
      }

      case MINLIST:
      case MAXLIST:
      {
         expr_va* amplexpr_va;
         const de* d;
         SCIP_EXPR* arg;
         SCIP_EXPROP operand;

         operand = opnum == MINLIST ? SCIP_EXPR_MIN : SCIP_EXPR_MAX;
         amplexpr_va = (expr_va*)amplexpr;

         arg = NULL;
         *scipexpr = NULL;
         for( d = amplexpr_va->L.d; d->e; ++d )
         {
            if( *scipexpr == NULL )
            {
               SCIP_CALL( walkExpression(scip, scipexpr, asl, d->e, exprvaridx, nexprvars, nvars, doingfine) );
               if( !*doingfine )
               {
                  SCIPexprFreeDeep(SCIPblkmem(scip), scipexpr);
                  break;
               }
            }
            else
            {
               SCIP_CALL( walkExpression(scip, &arg, asl, d->e, exprvaridx, nexprvars, nvars, doingfine) );
               if( !*doingfine )
               {
                  SCIPexprFreeDeep(SCIPblkmem(scip), &arg);
                  break;
               }
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), scipexpr, operand, *scipexpr, arg) );
            }
         }
         assert(*scipexpr != NULL); /* empty list?? */

         break;
      }

      case OP_cos:
         SCIPerrorMessage("AMPL operand number OP_cos not supported so far.\n");
         *doingfine = FALSE;
         break;

      case OP_sin:
         SCIPerrorMessage("AMPL operand number OP_sin not supported so far.\n");
         *doingfine = FALSE;
         break;

      default:
         SCIPerrorMessage("AMPL operand number %d not supported so far.\n", opnum);
         *doingfine = FALSE;
         break;
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE setupObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_Bool*            success             /**< indicates whether objective was successfully setup */
   )
{
   ASL* asl;
   fint* rowqp;
   fint* colqp;
   real* delsqp;
   int nqpterms;
   struct ograd* og;
   SCIP_Real objconstant;

   assert(scip != NULL);
   assert(probdata != NULL);
   assert(probdata->asl != NULL);
   assert(probdata->vars != NULL || probdata->nvars == 0);
   assert(success != NULL);

   *success = TRUE;

   asl = probdata->asl;
   assert(n_obj >= 0);

   if( n_obj > 1 )
   {
      SCIPwarningMessage(scip, "Multiple objectives not supported by SCIP! Ignoring all other than for first one.\n");
   }
   if( n_obj == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPsetObjsense(scip, objtype[0] == 1 ? SCIP_OBJSENSE_MAXIMIZE : SCIP_OBJSENSE_MINIMIZE) );

   /* if objective is nonlinear, check if quadratic */
   nqpterms = nlo > 0 ? nqpcheck(0, &rowqp, &colqp, &delsqp) : 0;

   objconstant = 0.0;
   /* for nonlinear and quadratic constraints, objconst(0) is already part of that expression, somehow
    * for a linear constraint, we may have a nonlinear expression that consists of the constant only ...
    */
   if( nqpterms >= 0 && ((ASL_fg*)asl)->I.obj_de_->e != NULL )  /*lint !e826*/
   {
      assert(((ASL_fg*)asl)->I.obj_de_->e->op == f_OPNUM || (Intcast ((ASL_fg*)asl)->I.obj_de_->e->op) == OPNUM); /*lint !e826*/
      objconstant += ((expr_n*)((ASL_fg*)asl)->I.obj_de_->e)->v;  /*lint !e826*/
   }

   if( nqpterms == 0 )
   {
      /* linear objective */
      SCIPdebugMessage("objective is linear\n");  /*lint !e534*/

      for( og = Ograd[0]; og; og = og->next )
      {
         assert(og->varno >= 0);
         assert(og->varno < probdata->nvars);

         SCIP_CALL( SCIPchgVarObj(scip, probdata->vars[og->varno], og->coef) );
      }

      /* handle objective constant by adding a fixed variable for it */
      if( !SCIPisZero(scip, objconstant) )
      {
         SCIP_VAR* objconstvar;

         SCIP_CALL( SCIPcreateVarBasic(scip, &objconstvar, "objconstant", objconstant, objconstant, 1.0, SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL( SCIPaddVar(scip, objconstvar) );
         SCIP_CALL( SCIPreleaseVar(scip, &objconstvar) );
      }
   }
   else
   {
      SCIP_CONS* obj;
      SCIP_VAR* objvar;
      SCIP_Real minusone;

      minusone = -1.0;

      SCIP_CALL( SCIPcreateVarBasic(scip, &objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, objvar) );

      if( nqpterms > 0 )
      {
         /* quadratic */
         int i;
         int j;

         SCIPdebugMessage("objective is quadratic\n");  /*lint !e534*/

         assert(rowqp != NULL);  /*lint !e644*/
         assert(colqp != NULL);  /*lint !e644*/
         assert(delsqp != NULL); /*lint !e644*/

         SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &obj, obj_name(0), 1, &objvar, &minusone, 0, NULL, NULL, NULL,
            objtype[0] == 1 ? -objconstant : -SCIPinfinity(scip),
            objtype[0] == 1 ? SCIPinfinity(scip) : -objconstant) );

         /* add quadratic coefficients of objective */
         for( i = 0; i < n_var; ++i )
            for( j = colqp[i]; j < colqp[i+1]; ++j, ++delsqp )
               if( !SCIPisZero(scip, *delsqp) )
               {
                  SCIP_CALL( SCIPaddBilinTermQuadratic(scip, obj, probdata->vars[i], probdata->vars[rowqp[j]], 0.5 * *delsqp) );
               }
         assert(colqp[n_var] == nqpterms);

         /* add linear part of objective */
         for( og = Ograd[0]; og; og = og->next )
         {
            assert(og->varno >= 0);
            assert(og->varno < probdata->nvars);

            if( !SCIPisZero(scip, og->coef) )
            {
               SCIP_CALL( SCIPaddLinearVarQuadratic(scip, obj, probdata->vars[og->varno], og->coef) );
            }
         }
      }
      else
      {
         SCIP_EXPRTREE* exprtree;
         SCIP_EXPR* objexpr;
         SCIP_VAR** exprvars;
         int* exprvaridx;
         int nexprvars;
         int i;

         SCIPdebugMessage("objective is nonlinear\n");  /*lint !e534*/

         SCIP_CALL( SCIPallocBufferArray(scip, &exprvaridx, probdata->nvars) );
         for( i = 0; i < probdata->nvars; ++i )
            exprvaridx[i] = -1;
         nexprvars = 0;

         SCIP_CALL( walkExpression(scip, &objexpr, probdata->asl, ((ASL_fg*)asl)->I.obj_de_->e, exprvaridx, &nexprvars, probdata->nvars, success) );  /*lint !e826*/

         if( !*success )
         {
            SCIPfreeBufferArray(scip, &exprvaridx);
            SCIP_CALL( SCIPreleaseVar(scip, &objvar) );

            return SCIP_OKAY;
         }

         /* assemble array exprvars with SCIP_VAR*'s */
         SCIP_CALL( SCIPallocBufferArray(scip, &exprvars, nexprvars) );
         for( i = 0; i < probdata->nvars; ++i )
         {
            assert(exprvaridx[i] < nexprvars );

            if( exprvaridx[i] >= 0 )
               exprvars[exprvaridx[i]] = probdata->vars[i];  /*lint !e613*/
         }

         /* create expression tree */
         SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, objexpr, nexprvars, 0, NULL) );
         SCIP_CALL( SCIPexprtreeSetVars(exprtree, nexprvars, exprvars) );

         /* nonlinear nonquadratic */
         SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &obj, obj_name(0), 1, &objvar, &minusone, 1, &exprtree, NULL,
            objtype[0] == 1 ? -objconstant : -SCIPinfinity(scip),
            objtype[0] == 1 ? SCIPinfinity(scip) : -objconstant) );

         /* add linear part of objective */
         for( og = Ograd[0]; og; og = og->next )
         {
            assert(og->varno >= 0);
            assert(og->varno < probdata->nvars);

            if( !SCIPisZero(scip, og->coef) )
            {
               SCIP_CALL( SCIPaddLinearVarNonlinear(scip, obj, probdata->vars[og->varno], og->coef) );
            }
         }

         SCIP_CALL( SCIPexprtreeFree(&exprtree) );
         SCIPfreeBufferArray(scip, &exprvars);
         SCIPfreeBufferArray(scip, &exprvaridx);
      }

      /* add objective constraint and forget */
      SCIP_CALL( SCIPaddCons(scip, obj) );
      SCIP_CALL( SCIPreleaseCons(scip, &obj) );

      /* forget objective variable */
      SCIP_CALL( SCIPreleaseVar(scip, &objvar) );
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE setupConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_Bool*            success
   )
{
   SCIP_VAR** exprvars;
   int* exprvaridx;
   int nexprvars;
   SCIP_CONS* cons;
   ASL* asl;
   fint* rowqp;
   fint* colqp;
   real* delsqp;
   int nqpterms;
   struct cgrad* cg;
   int c;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SufDesc* suf_initial = NULL;
   SufDesc* suf_separate = NULL;
   SufDesc* suf_enforce = NULL;
   SufDesc* suf_check = NULL;
   SufDesc* suf_propagate = NULL;
   SufDesc* suf_dynamic = NULL;
   SufDesc* suf_removable = NULL;
   SufDesc* suf_curvature = NULL;
   SCIP_Bool initial;
   SCIP_Bool separate;
   SCIP_Bool enforce;
   SCIP_Bool check;
   SCIP_Bool propagate;
   SCIP_Bool dynamic;
   SCIP_Bool removable;

   assert(scip != NULL);
   assert(probdata != NULL);
   assert(probdata->asl != NULL);
   assert(probdata->vars != NULL || probdata->nvars == 0);
   assert(success != NULL);

   *success = TRUE;

   asl = probdata->asl;
   assert(n_con >= 0);

   /* get AMPL suffix values */
   suf_initial = suf_get("initial", ASL_Sufkind_con);
   if( suf_initial != NULL && suf_initial->u.i == NULL )
      suf_initial = NULL;
   suf_separate = suf_get("separate", ASL_Sufkind_con);
   if( suf_separate != NULL && suf_separate->u.i == NULL )
      suf_separate = NULL;
   suf_enforce = suf_get("enforce", ASL_Sufkind_con);
   if( suf_enforce != NULL && suf_enforce->u.i == NULL )
      suf_enforce = NULL;
   suf_check = suf_get("check", ASL_Sufkind_con);
   if( suf_check != NULL && suf_check->u.i == NULL )
      suf_check = NULL;
   suf_propagate = suf_get("propagate", ASL_Sufkind_con);
   if( suf_propagate != NULL && suf_propagate->u.i == NULL )
      suf_propagate = NULL;
   suf_dynamic = suf_get("dynamic", ASL_Sufkind_con);
   if( suf_dynamic != NULL && suf_dynamic->u.i == NULL )
      suf_dynamic = NULL;
   suf_removable = suf_get("removable", ASL_Sufkind_con);
   if( suf_removable != NULL && suf_removable->u.i == NULL )
      suf_removable = NULL;

#ifdef CHECKCURVSUFFIX
   suf_curvature = suf_get("curvature", ASL_Sufkind_con);
   if( suf_curvature != NULL )
   {
      if( suf_curvature->u.i != NULL )
      {
         SCIPinfoMessage(scip, NULL, "Found curvature suffix for %d constraints in .nl file.\n", (&asl->i.n_var_)[ASL_Sufkind_con & ASL_Sufkind_mask]);
      }
      else
         suf_curvature = NULL;
   }
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &exprvaridx, probdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &exprvars, probdata->nvars) );

   for( c = 0; c < n_con && *success; ++c )
   {
      /* nqpcheck() need to be before LUrhs is processed, as it puts the constant of the quad. term into the sides */
      nqpterms = c < nlc ? nqpcheck(-(c+1), &rowqp, &colqp, &delsqp) : 0;

      lhs = LUrhs[2*c];
      if( SCIPisInfinity(scip, -lhs) )
         lhs = -SCIPinfinity(scip);

      rhs = LUrhs[2*c+1];
      if( SCIPisInfinity(scip, rhs) )
         rhs = SCIPinfinity(scip);

      if( suf_initial != NULL && suf_initial->u.i[c] > 0 )
      {
         assert(suf_initial->u.i[c] == 1 || suf_initial->u.i[c] == 2);
         initial = (suf_initial->u.i[c] == 1);
      }
      else
         initial = TRUE;

      if( suf_separate != NULL && suf_separate->u.i[c] > 0 )
      {
         assert(suf_separate->u.i[c] == 1 || suf_separate->u.i[c] == 2);
         separate = (suf_separate->u.i[c] == 1);
      }
      else
         separate = TRUE;

      if( suf_enforce != NULL && suf_enforce->u.i[c] > 0 )
      {
         assert(suf_enforce->u.i[c] == 1 || suf_enforce->u.i[c] == 2);
         enforce = (suf_enforce->u.i[c] == 1);
      }
      else
         enforce = TRUE;

      if( suf_check != NULL && suf_check->u.i[c] > 0 )
      {
         assert(suf_check->u.i[c] == 1 || suf_check->u.i[c] == 2);
         check = (suf_check->u.i[c] == 1);
      }
      else
         check = TRUE;

      if( suf_propagate != NULL && suf_propagate->u.i[c] > 0 )
      {
         assert(suf_propagate->u.i[c] == 1 || suf_propagate->u.i[c] == 2);
         propagate = (suf_propagate->u.i[c] == 1);
      }
      else
         propagate = TRUE;

      if( suf_dynamic != NULL && suf_dynamic->u.i[c] > 0 )
      {
         assert(suf_dynamic->u.i[c] == 1 || suf_dynamic->u.i[c] == 2);
         dynamic = (suf_dynamic->u.i[c] == 1);
      }
      else
         dynamic = FALSE;

      if( suf_removable != NULL && suf_removable->u.i[c] > 0 )
      {
         assert(suf_removable->u.i[c] == 1 || suf_removable->u.i[c] == 2);
         removable = (suf_removable->u.i[c] == 1);
      }
      else
         removable = FALSE;

      if( nqpterms == 0 )
      {
         /* linear */
         SCIPdebugMessage("constraint %d (%s) is linear\n", c, con_name(c));  /*lint !e534*/

         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, con_name(c), 0, NULL, NULL, lhs, rhs,
            initial, separate, enforce, check, propagate, FALSE, FALSE, dynamic, removable, FALSE) );

         /* add linear coefficients */
         for( cg = Cgrad[c]; cg; cg = cg->next )
            if( !SCIPisZero(scip, cg->coef) )
            {
               SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->vars[cg->varno], cg->coef) );
            }

      }
      else if( suf_curvature != NULL && suf_curvature->u.i[c] > 0 )
      {
         SCIP_USEREXPRDATA* userexprdata;
         SCIP_EXPRTREE* exprtree;
         SCIP_EXPR** childexprs;
         SCIP_EXPR* consexpr;
         int i;

         SCIPdebugMessage("constraint %d (%s) has suffix %d, handle as nonlinear userexpr\n", c, con_name(c), suf_curvature->u.i[c]);  /*lint !e534*/
         assert(suf_curvature->u.i[c] == 1 || suf_curvature->u.i[c] == 2);

         /* TODO would be good to treat linear variables separately */
         nexprvars = 0;
         for( cg = Cgrad[c]; cg != NULL; cg = cg->next )
            exprvars[nexprvars++] = probdata->vars[cg->varno];

         SCIP_CALL( SCIPallocBufferArray(scip, &childexprs, nexprvars) );
         for( i = 0; i < nexprvars; ++i )
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &childexprs[i], SCIP_EXPR_VARIDX, i) );
         }

         SCIP_CALL( SCIPallocBlockMemory(scip, &userexprdata) );
         userexprdata->probdata = probdata;
         userexprdata->considx = c;
         userexprdata->curvature = (suf_curvature->u.i[c] == 1 ? SCIP_EXPRCURV_CONVEX : SCIP_EXPRCURV_CONCAVE);

         SCIP_CALL( SCIPexprCreateUser(SCIPblkmem(scip), &consexpr, nexprvars, childexprs, userexprdata,
            SCIP_EXPRINTCAPABILITY_FUNCVALUE | SCIP_EXPRINTCAPABILITY_INTFUNCVALUE | SCIP_EXPRINTCAPABILITY_GRADIENT | SCIP_EXPRINTCAPABILITY_INTGRADIENT,
            SCIPuserexprEvalAmpl, SCIPuserexprIntEvalAmpl,
            SCIPuserexprCurvAmpl, NULL, NULL,
               SCIPuserexprCopyAmpl, SCIPuserexprFreeAmpl, NULL) );

         SCIPfreeBufferArray(scip, &childexprs);

         /* create expression tree */
         SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, consexpr, nexprvars, 0, NULL) );
         SCIP_CALL( SCIPexprtreeSetVars(exprtree, nexprvars, exprvars) );

         /* create nonlinear constraint */
         SCIP_CALL( SCIPcreateConsNonlinear(scip, &cons, con_name(c), 0, NULL, NULL, 1, &exprtree, NULL, lhs, rhs,
            initial, separate, enforce, check, propagate, FALSE, FALSE, dynamic, removable, FALSE) );

         /* constraint has copy of tree, so free */
         SCIP_CALL( SCIPexprtreeFree(&exprtree) );

         if( probdata->fullx == NULL )
         {
            /* alloc scratch memory for userexpr function evaluations */
            SCIP_CALL( SCIPallocClearMemoryArray(scip, &probdata->fullx, probdata->nvars) );
         }
      }
      else if( nqpterms > 0 )
      {
         /* quadratic */
         int i;
         int j;

         assert(rowqp != NULL);  /*lint !e644*/
         assert(colqp != NULL);  /*lint !e644*/
         assert(delsqp != NULL); /*lint !e644*/

         SCIPdebugMessage("constraint %d (%s) is quadratic\n", c, con_name(c));  /*lint !e534*/

         SCIP_CALL( SCIPcreateConsQuadratic(scip, &cons, con_name(c), 0, NULL, NULL, 0, NULL, NULL, NULL, lhs, rhs,
            initial, separate, enforce, check, propagate, FALSE, FALSE, dynamic, removable) );

         /* add quadratic coefficients of constraint */
         for( i = 0; i < n_var; ++i )
            for( j = colqp[i]; j < colqp[i+1]; ++j, ++delsqp )
               if( !SCIPisZero(scip, *delsqp) )
               {
                  SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, probdata->vars[i], probdata->vars[rowqp[j]], 0.5 * *delsqp) );
               }
         assert(colqp[n_var] == nqpterms);

         /* add linear part of constraint */
         for( cg = Cgrad[c]; cg; cg = cg->next )
            if( !SCIPisZero(scip, cg->coef) )
            {
               SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, probdata->vars[cg->varno], cg->coef) );
            }
      }
      else
      {
         /* nonlinear non-quadratic */
         SCIP_EXPRTREE* exprtree;
         SCIP_EXPR* consexpr;
         int i;

         SCIPdebugMessage("constraint %d (%s) is nonlinear\n", c, con_name(c));  /*lint !e534*/

         for( i = 0; i < probdata->nvars; ++i )
            exprvaridx[i] = -1;
         nexprvars = 0;

         SCIP_CALL( walkExpression(scip, &consexpr, probdata->asl, ((ASL_fg*)asl)->I.con_de_[c].e, exprvaridx, &nexprvars, probdata->nvars, success) );  /*lint !e826*/

         if( !*success )
         {
            SCIPfreeBufferArray(scip, &exprvaridx);
            SCIPfreeBufferArray(scip, &exprvars);

            break;
         }

         /* assemble array exprvars with SCIP_VAR*'s */
         for( i = 0; i < probdata->nvars; ++i )
         {
            assert(exprvaridx[i] < nexprvars );

            if( exprvaridx[i] >= 0 )
               exprvars[exprvaridx[i]] = probdata->vars[i];  /*lint !e613*/
         }

         /* create expression tree */
         SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, consexpr, nexprvars, 0, NULL) );
         SCIP_CALL( SCIPexprtreeSetVars(exprtree, nexprvars, exprvars) );

         /* expression trees without variables raise assertions in CppAD, workaround here for now */
         if( nexprvars > 0 )
         {
            SCIP_CALL( SCIPcreateConsNonlinear(scip, &cons, con_name(c), 0, NULL, NULL, 1, &exprtree, NULL, lhs, rhs,
               initial, separate, enforce, check, propagate, FALSE, FALSE, dynamic, removable, FALSE) );

            /* add linear part of constraint */
            for( cg = Cgrad[c]; cg; cg = cg->next )
               if( !SCIPisZero(scip, cg->coef) )
               {
                  SCIP_CALL( SCIPaddLinearVarNonlinear(scip, cons, probdata->vars[cg->varno], cg->coef) );
               }
         }
         else
         {
            SCIP_Real val;

            SCIP_CALL( SCIPexprtreeEval(exprtree, NULL, &val) );

            if( !SCIPisInfinity(scip, -lhs) )
               lhs -= val;

            if( !SCIPisInfinity(scip,  rhs) )
               rhs -= val;

            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, con_name(c), 0, NULL, NULL, lhs, rhs,
               initial, separate, enforce, check, propagate, FALSE, FALSE, dynamic, removable, FALSE) );

            /* add linear part of constraint */
            for( cg = Cgrad[c]; cg; cg = cg->next )
            {
               SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->vars[cg->varno], cg->coef) );
            }
         }

         SCIP_CALL( SCIPexprtreeFree(&exprtree) );
      }

      assert(cons != NULL);
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   SCIPfreeBufferArray(scip, &exprvaridx);
   SCIPfreeBufferArray(scip, &exprvars);

   return SCIP_OKAY;
}

static
SCIP_RETCODE setupSOS(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_Bool*            success
   )
{
   SCIP_CONS* cons;
   SCIP_VAR** vars;
   int varssize;
   char name[20];
   ASL* asl;
   int copri[2] = {0,0};
   int* starts = NULL;
   int* indices = NULL;
   char* types = NULL;
   SCIP_Real* weights = NULL;
   int* priorities = NULL;
   int num;
   int numNz;
   int i;
   int j;
   int k;

   assert(scip != NULL);
   assert(probdata != NULL);
   assert(probdata->asl != NULL);
   assert(probdata->vars != NULL || probdata->nvars == 0);
   assert(success != NULL);

   *success = TRUE;

   asl = probdata->asl;

   num = suf_sos(0 /*flags*/, &numNz, &types, &priorities, copri, &starts, &indices, &weights);

   if( num == 0 )
      return SCIP_OKAY;
   assert(num > 0);

   varssize = 10;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, varssize) );

   for( i = 0; i < num; ++i )
   {
      if( starts[i+1] - starts[i] > varssize )
      {
         varssize = SCIPcalcMemGrowSize(scip, starts[i+1] - starts[i]);
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, varssize) );
      }

      for( j = starts[i], k = 0; j < starts[i+1]; ++j, ++k )
         vars[k] = probdata->vars[indices[j]];

      sprintf(name, "sos%d", i);
      switch( types[i] )
      {
         case '1' :
         {
            SCIP_CALL( SCIPcreateConsSOS1(scip, &cons, name, k, vars, &weights[starts[i]], TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
            break;
         }

         case '2' :
         {
            SCIP_CALL( SCIPcreateConsSOS2(scip, &cons, name, k, vars, &weights[starts[i]], TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
            break;
         }

         default :
         {
            SCIPerrorMessage("SOS type %c not supported by SCIP.\n", types[i]);
            *success = FALSE;
            break;
         }
      }

      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

static
SCIP_RETCODE setupInitialSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   SCIP_Bool stored;
   SCIP_SOL* sol;
   ASL* asl;
   int i;

   assert(scip != NULL);
   assert(probdata != NULL);

   asl = probdata->asl;

   /* if no initial guess available, then do nothing */
   if( X0 == NULL )
      return SCIP_OKAY;

   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );

   for( i = 0; i < n_var; ++i )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, probdata->vars[i], X0[i]) );
   }

   SCIP_CALL( SCIPaddSolFree(scip, &sol, &stored) );

   return SCIP_OKAY;
}

/*
 * Callback methods of probdata
 */

/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdataDelOrigNl)
{
   int i;

   assert((*probdata)->asl != NULL);
   assert((*probdata)->vars != NULL || (*probdata)->nvars == 0);

   ASL_free(&(*probdata)->asl);

   for( i = 0; i < (*probdata)->nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->vars[i]) );
   }
   SCIPfreeMemoryArrayNull(scip, &(*probdata)->vars);

   SCIPfreeMemoryArrayNull(scip, &(*probdata)->fullx);

   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}

/*
 * Callback methods of reader
 */

#if SCIP_DISABLED_CODE /* TODO: implement, if one finds use for it */
/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyNl)
{  /*lint --e{715}*/
   assert(scip != NULL);

   SCIP_CALL( SCIPincludeReaderNl(scip) );

   return SCIP_OKAY;
}
#endif

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadNl)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_Bool success;
   const char* filebasename;
   FILE* nl;
   ASL* asl;
   SufDecl suftable[15];

   assert(scip != NULL);
   assert(reader != NULL);
   assert(filename != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* initialize ASL */
   asl = ASL_alloc(ASL_read_fg);

   /* make ASL aware of SOS suffixes */
   suftable[0].name = (char*)"ref";
   suftable[0].kind = ASL_Sufkind_var | ASL_Sufkind_real;
   suftable[0].nextra = 0;

   suftable[1].name = (char*)"sos";
   suftable[1].kind = ASL_Sufkind_var;
   suftable[1].nextra = 0;

   suftable[2].name = (char*)"sos";
   suftable[2].kind = ASL_Sufkind_con;
   suftable[2].nextra = 0;

   suftable[3].name = (char*)"sosno";
   suftable[3].kind = ASL_Sufkind_var | ASL_Sufkind_real;
   suftable[3].nextra = 0;

   suftable[4].name = (char*)"priority";
   suftable[4].kind = ASL_Sufkind_var;
   suftable[4].nextra = 0;

   /* make ASL aware of curvature suffix */
   suftable[5].name = (char*)"curvature";
   suftable[5].kind = ASL_Sufkind_con;
   suftable[5].nextra = 0;

   /* make ASL aware of suffixes for constraint flags */
   suftable[6].name = (char*)"initial";  /* should constraint be in initial LP? */
   suftable[6].kind = ASL_Sufkind_con;
   suftable[6].nextra = 0;

   suftable[7].name = (char*)"separate";  /* should constraint be in separated? */
   suftable[7].kind = ASL_Sufkind_con;
   suftable[7].nextra = 0;

   suftable[8].name = (char*)"enforce";  /* should constraint be enforced? */
   suftable[8].kind = ASL_Sufkind_con;
   suftable[8].nextra = 0;

   suftable[9].name = (char*)"check";  /* should constraint be checked? */
   suftable[9].kind = ASL_Sufkind_con;
   suftable[9].nextra = 0;

   suftable[10].name = (char*)"propagate";  /* should constraint be propagated? */
   suftable[10].kind = ASL_Sufkind_con;
   suftable[10].nextra = 0;

   suftable[11].name = (char*)"dynamic";  /* should constraint be subject to aging? */
   suftable[11].kind = ASL_Sufkind_con;
   suftable[11].nextra = 0;

   suftable[12].name = (char*)"removable";  /* should  the relaxation be removed from the LP due to aging or cleanup? */
   suftable[12].kind = ASL_Sufkind_con;
   suftable[12].nextra = 0;

   /* make ASL aware of suffixes for variable flags */
   suftable[13].name = (char*)"initial";  /* should variable be in initial LP? */
   suftable[13].kind = ASL_Sufkind_var;
   suftable[13].nextra = 0;

   suftable[14].name = (char*)"removable";  /* is var's column removable from the LP (due to aging or cleanup)? */
   suftable[14].kind = ASL_Sufkind_var;
   suftable[14].nextra = 0;

   suf_declare(suftable, 15);

   /* let ASL read .nl file
    * jac0dim will do exit(1) if the file is not found
    */
   nl = jac0dim((char*)filename, (fint)strlen(filename));

   if( nl == NULL )
   {
      SCIPerrorMessage("error processing <%s> file\n", filename);
      return SCIP_READERROR;
   }

   want_derivs = 0;
   want_xpi0 = 1;
   (void) qp_read(nl, 0);

   *result = SCIP_DIDNOTFIND;

   SCIPdebugMessage("number of nonzeros    = %d\n", nzc);     /*lint !e534*/
   SCIPdebugMessage("number of variables   = %d\n", n_var);   /*lint !e534*/
   SCIPdebugMessage("number of constraints = %d\n", n_con);   /*lint !e534*/
   SCIPdebugMessage("number of objectives  = %d\n", n_obj);   /*lint !e534*/
   SCIPdebugMessage("number of ranges      = %d\n", nranges); /*lint !e534*/
   SCIPdebugMessage("number of equations   = %d\n", n_eqn);   /*lint !e534*/

   SCIP_CALL( SCIPallocMemory(scip, &probdata) );
   BMSclearMemory(probdata);

   probdata->asl = asl;

   /* initialize empty SCIP problem */
   filebasename = strrchr(filename, '/');
   if( filebasename == NULL )
      filebasename = filename;
   else
      ++filebasename;
   SCIP_CALL( SCIPcreateProbBasic(scip, filebasename) );
   SCIP_CALL( SCIPsetProbData(scip, probdata) );
   SCIP_CALL( SCIPsetProbDelorig(scip, probdataDelOrigNl) );

   SCIP_CALL( setupVariables(scip, probdata) );
   SCIP_CALL( setupInitialSolution(scip, probdata) );

   SCIP_CALL( setupObjective(scip, probdata, &success) );
   if( !success )
      return SCIP_READERROR;

   SCIP_CALL( setupConstraints(scip, probdata, &success) );
   if( !success )
      return SCIP_READERROR;

   SCIP_CALL( setupSOS(scip, probdata, &success) );
   if( !success )
      return SCIP_READERROR;

   *result = SCIP_SUCCESS;

   if( probdata->fullx != NULL )
   {
      /* If fullx is allocated, then we have a user expression, so we need function values and gradients from ASL.
       * When reading the .nl file with qp_read(), we cannot call function/gradient evaluations.
       * A call to qp_opify() is supposed to restore this behavior, but doesn't seem to do this for quadratic constraints (functionvalue=0).
       * Thus, we reread the .nl file again with fg_read().
       * Note: To get Hessians, we would need to use pfgh_read().
       * Alternatively, we could try forbidding user expressions for quadratic constraints (doesn't seem useful).
       */
      ASL_free(&probdata->asl);

      asl = ASL_alloc(ASL_read_fg);
      probdata->asl = asl;
      nl = jac0dim((char*)filename, (fint)strlen(filename));

      want_derivs = 1; /* we want derivatives */
      want_xpi0 = 1; /* we want the initial guess (primal only) */
      asl->i.congrd_mode = 1; /* ask for compact gradient */
      (void) fg_read(nl, ASL_return_read_err);
   }

   return SCIP_OKAY;
}


/*
 * Writing AMPL solution file dialog
 */

#define DIALOG_WRITEAMPLSOL_NAME             "amplsol"
#define DIALOG_WRITEAMPLSOL_DESC             "writes AMPL solution file"
#define DIALOG_WRITEAMPLSOL_ISSUBMENU        FALSE

/** execution method of dialog */
static
SCIP_DECL_DIALOGEXEC(dialogExecWriteAmplSol)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(dialoghdlr != NULL);
   assert(dialog != NULL);

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   if( SCIPgetStage(scip) == SCIP_STAGE_INIT || SCIPgetStage(scip) == SCIP_STAGE_FREE )
   {
      SCIPerrorMessage("No AMPL problem read, cannot write AMPL solution then.\n");
      return SCIP_OKAY;
   }

   /* currently, can pass NULL as reader
    * in the future, the reader should be changed to store its data in the readerdata
    * instead of the problem data
    */
   SCIP_CALL( SCIPwriteAmplSolReaderNl(scip, NULL) );

   return SCIP_OKAY;
}


/*
 * reader specific interface methods
 */

/** includes the AMPL .nl file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderNl(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;
   SCIP_DIALOG* dialog;
   SCIP_DIALOG* parentdialog;

   readerdata = NULL;
   reader = NULL;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );
   assert(reader != NULL);

   /* set non fundamental callbacks via setter functions */
   /* SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyNl) ); */
   /* SCIP_CALL( SCIPsetReaderFree(scip, reader, readerFreeNl) ); */
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadNl) );

   if( SCIPgetRootDialog(scip) != NULL )
   {
      /* get parent dialog "write" */
      if( SCIPdialogFindEntry(SCIPgetRootDialog(scip), "write", &parentdialog) != 1 )
      {
         SCIPerrorMessage("sub menu \"write\" not found\n");
         return SCIP_PLUGINNOTFOUND;
      }
      assert(parentdialog != NULL);

      /* create, include, and release dialog */
      if( !SCIPdialogHasEntry(parentdialog, DIALOG_WRITEAMPLSOL_NAME) )
      {
         SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL, dialogExecWriteAmplSol, NULL, NULL,
            DIALOG_WRITEAMPLSOL_NAME, DIALOG_WRITEAMPLSOL_DESC, DIALOG_WRITEAMPLSOL_ISSUBMENU, NULL) );
         SCIP_CALL( SCIPaddDialogEntry(scip, parentdialog, dialog) );
         SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, "ASL", "AMPL Solver Library developed by D. Gay (www.netlib.com/ampl)") );

   return SCIP_OKAY;
}


/** writes AMPL solution file */
SCIP_RETCODE SCIPwriteAmplSolReaderNl(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          nlreader            /**< AMPL .nl file reader */
   )
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_Real* x;
   SCIP_Real* y;
   const char* msg;
   ASL* asl;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   if( probdata == NULL || probdata->asl == NULL )
   {
      SCIPerrorMessage("No AMPL problem read, cannot write AMPL solution then.\n");
      return SCIP_OKAY;
   }

   asl = probdata->asl;

   /* AMPL solve status codes are at http://www.ampl.com/NEW/statuses.html
    *     number   string       interpretation
    *    0 -  99   solved       optimal solution found
    *  100 - 199   solved?      optimal solution indicated, but error likely
    *  200 - 299   infeasible   constraints cannot be satisfied
    *  300 - 399   unbounded    objective can be improved without limit
    *  400 - 499   limit        stopped by a limit that you set (such as on iterations)
    *  500 - 599   failure      stopped by an error condition in the solver routines
    */

   switch( SCIPgetStatus(scip) )
   {
      case SCIP_STATUS_UNKNOWN:
         msg = "unknown";
         solve_result_num = 500;
         break;
      case SCIP_STATUS_USERINTERRUPT:
         msg = "user interrupt";
         solve_result_num = 450;
         break;
      case SCIP_STATUS_NODELIMIT:
         msg = "node limit reached";
         solve_result_num = 400;
         break;
      case SCIP_STATUS_TOTALNODELIMIT:
         msg = "total node limit reached";
         solve_result_num = 401;
         break;
      case SCIP_STATUS_STALLNODELIMIT:
         msg = "stall node limit reached";
         solve_result_num = 402;
         break;
      case SCIP_STATUS_TIMELIMIT:
         msg = "time limit reached";
         solve_result_num = 403;
         break;
      case SCIP_STATUS_MEMLIMIT:
         msg = "memory limit reached";
         solve_result_num = 404;
         break;
      case SCIP_STATUS_GAPLIMIT:
         msg = "gap limit reached";
         solve_result_num = 405;
         break;
      case SCIP_STATUS_SOLLIMIT:
         msg = "solution limit reached";
         solve_result_num = 406;
         break;
      case SCIP_STATUS_BESTSOLLIMIT:
         msg = "solution improvement limit reached";
         solve_result_num = 407;
         break;
      case SCIP_STATUS_OPTIMAL:
         msg = "optimal solution found";
         solve_result_num = 0;
         break;
      case SCIP_STATUS_INFEASIBLE:
         msg = "infeasible";
         solve_result_num = 200;
         break;
      case SCIP_STATUS_UNBOUNDED:
         msg = "unbounded";
         solve_result_num = 300;
         break;
      case SCIP_STATUS_INFORUNBD:
         msg = "infeasible or unbounded";
         solve_result_num = 299;
         break;
      default:
         solve_result_num = 500;
         SCIPerrorMessage("invalid status code <%d>\n", SCIPgetStatus(scip));
         return SCIP_INVALIDDATA;
   }

   /* get best primal solution */
   x = NULL;
   if( SCIPgetBestSol(scip) != NULL )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &x, probdata->nvars) );
      SCIP_CALL( SCIPgetSolVals(scip, SCIPgetBestSol(scip), probdata->nvars, probdata->vars, x) );
   }

   /* if the problem is an LP and presolving was turned off, then we try to return the vector of dual multipliers */
   y = NULL;
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVED && !SCIPhasPerformedPresolve(scip) && SCIPgetNVars(scip) == SCIPgetNContVars(scip) )
   {
      SCIP_CONSHDLR* linconshdlr;
      SCIP_CONS** conss;
      int nconss;
      int c;

      linconshdlr = SCIPfindConshdlr(scip, "linear");
      assert(linconshdlr != NULL);

      conss = SCIPgetConss(scip);
      nconss = SCIPgetNConss(scip);
      assert(conss != NULL);
      assert(nconss >= 0);

      SCIP_CALL( SCIPallocBufferArray(scip, &y, nconss) );

      for( c = 0; c < SCIPgetNConss(scip); ++c )
      {
         SCIP_CONS* transcons;

         /* dual solution is created by LP solver and therefore only available for linear constraints */
         SCIP_CALL( SCIPgetTransformedCons(scip, conss[c], &transcons) );
         if( transcons == NULL || SCIPconsGetHdlr(transcons) != linconshdlr )
         {
            SCIPfreeBufferArray(scip, &y);
            y = NULL;
            break;
         }

         y[c] = SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE
            ? SCIPgetDualsolLinear(scip, transcons)
            : -SCIPgetDualsolLinear(scip, transcons);
         assert(y[c] != SCIP_INVALID); /*lint !e777*/
      }
   }

   write_sol((char*)msg, x, y, NULL);

   SCIPfreeBufferArrayNull(scip, &x);
   SCIPfreeBufferArrayNull(scip, &y);

   return SCIP_OKAY;
}


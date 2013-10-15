/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cmain.c
 * @brief  AMPL .nl file reader
 * @author Stefan Vigerske
 *
 * The code (including some comments) for this reader is based on OSnl2osil.cpp,
 * the nl-reader of the Optimization Services Project (https://projects.coin-or.org/OS).
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

#if 0
/** reader data */
struct SCIP_ReaderData
{
};
#endif

/** problem data */
struct SCIP_ProbData
{
   ASL*                  asl;                /**< ASL data structure */

   SCIP_VAR**            vars;
   int                   nvars;
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */

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

   assert(scip != NULL);
   assert(probdata != NULL);
   assert(probdata->asl != NULL);
   assert(probdata->vars == NULL);

   asl = probdata->asl;

   SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->vars, n_var) );
   probdata->nvars = n_var;

   /* first the nonlinear variables */
   /* welcome to the world of the ASL API */

   lower = 0;
   upper = nlvb - nlvbi;
   for( i = lower; i < upper; ++i ) /* continuous and in an objective and in a constraint */
   {
      SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }

   lower = nlvb - nlvbi;
   upper = nlvb;
   for( i = lower; i < upper; ++i ) /* integer and in an objective and in a constraint */
   {
       SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_INTEGER) );
       SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }

   lower = nlvb;
   /* upper = nlvb + (nlvc - (nlvb + nlvci)); */
   upper = nlvc - nlvci;
   for( i = lower; i < upper; ++i ) /* continuous and just in constraints */
   {
      SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }

   lower = nlvc - nlvci;
   upper = nlvc;
   for( i = lower; i < upper; ++i ) /* integer and just in constraints */
   {
      SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_INTEGER) );
      SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }

   lower = nlvc;
   /* upper = nlvc + ( nlvo - (nlvc + nlvoi) ); */
   upper = nlvo - nlvoi;
   for( i = lower; i < upper; ++i) /* continuous and just in objectives */
   {
      SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }

   lower = nlvo - nlvoi;
   upper = nlvo ;
   for( i = lower; i < upper; ++i ) /* integer and just in objectives */
   {
      SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_INTEGER) );
      SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }

   /* now the other variables */

   lower = MAX(nlvc, nlvo);
   upper = MAX(nlvc, nlvo) + nwv;
   for( i = lower; i < upper; ++i ) /* linear arc variables */
   {
      SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }


   lower = MAX(nlvc, nlvo) + nwv;
   /* upper = MAX(nlvc, nlvo) + nwv + (n_var - (MAX(nlvc, nlvo) + niv + nbv + nwv) ); */
   upper = n_var -  niv - nbv;
   for( i = lower; i < upper; ++i ) /* other linear variables */
   {
      SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }

   lower = n_var -  niv - nbv;
   upper = n_var -  niv ;
   for( i = lower; i < upper; ++i ) /* linear and binary */
   {
      SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_BINARY) );
      SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }

   lower = n_var -  niv;
   upper = n_var;
   for( i = lower; i < upper; ++i ) /* linear and integer */
   {
      SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->vars[i], var_name(i), LUv[2*i], LUv[2*i+1], 0.0, SCIP_VARTYPE_INTEGER) );
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
            /* process common expression
             * see also http://www.gerad.ca/~orban/drampl/def-vars.html
             */
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

      case MAXLIST:
      {
         int nchildren;

         if( (amplexpr->L.ep == NULL) != (amplexpr->R.ep == NULL) )
         {
            /* what the ...? this seems to be possible, but the non-NULL operand seems to point to garbage
             * so assume that we have no operand
             */
            nchildren = 0;
         }
         else
         {
            nchildren = amplexpr->R.ep - amplexpr->L.ep;
         }

         assert(nchildren >= 0);
         switch( nchildren )
         {
            case 0:
            {
               SCIPwarningMessage(scip, NULL, "MAXLIST with 0 operands. Assuming max(empty list) = 0.\n");
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), scipexpr, SCIP_EXPR_CONST, 0.0) );
               break;
            }

            case 1:
            {
               SCIP_CALL( walkExpression(scip, scipexpr, asl, amplexpr->L.ep != NULL ? amplexpr->L.ep[0] : amplexpr->R.ep[0], exprvaridx, nexprvars, nvars, doingfine) );
               break;
            }

            default:
            {
               SCIP_EXPR* arg1;
               SCIP_EXPR* arg2;

               SCIP_CALL( walkExpression(scip, &arg1, asl, amplexpr->L.ep[0], exprvaridx, nexprvars, nvars, doingfine) );
               if( !*doingfine )
                  break;

               SCIP_CALL( walkExpression(scip, &arg2, asl, amplexpr->L.ep[1], exprvaridx, nexprvars, nvars, doingfine) );
               if( !*doingfine )
               {
                  SCIPexprFreeDeep(SCIPblkmem(scip), &arg1);
                  break;
               }

               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), scipexpr, SCIP_EXPR_MAX, arg1, arg2) );

               for( i = 2; i < nchildren; ++i )
               {
                  SCIP_CALL( walkExpression(scip, &arg1, asl, amplexpr->L.ep[i], exprvaridx, nexprvars, nvars, doingfine) );
                  if( !*doingfine )
                  {
                     SCIPexprFreeDeep(SCIPblkmem(scip), scipexpr);
                     break;
                  }

                  SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), scipexpr, SCIP_EXPR_MAX, *scipexpr, arg1) );
               }
               break;
            }
         }

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
    * for a linear constraint, we may have a nonlinear expression that consists of the constant only...
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

   assert(scip != NULL);
   assert(probdata != NULL);
   assert(probdata->asl != NULL);
   assert(probdata->vars != NULL || probdata->nvars == 0);
   assert(success != NULL);

   *success = TRUE;

   asl = probdata->asl;
   assert(n_con >= 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &exprvaridx, probdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &exprvars, probdata->nvars) );

   for( c = 0; c < n_con && *success; ++c )
   {
      nqpterms = c < nlc ? nqpcheck(-(c+1), &rowqp, &colqp, &delsqp) : 0;
      if( nqpterms == 0 )
      {
         /* linear */
         SCIPdebugMessage("constraint %d (%s) is linear\n", c, con_name(c));  /*lint !e534*/

         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, con_name(c), 0, NULL, NULL, LUrhs[2*c], LUrhs[2*c+1]) );

         /* add linear coefficients */
         for( cg = Cgrad[c]; cg; cg = cg->next )
            if( !SCIPisZero(scip, cg->coef) )
            {
               SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->vars[cg->varno], cg->coef) );
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

         SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, con_name(c), 0, NULL, NULL, 0, NULL, NULL, NULL, LUrhs[2*c], LUrhs[2*c+1]) );

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
            SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, con_name(c), 0, NULL, NULL, 1, &exprtree, NULL, LUrhs[2*c], LUrhs[2*c+1]) );

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
            SCIP_Real lhs;
            SCIP_Real rhs;

            SCIP_CALL( SCIPexprtreeEval(exprtree, NULL, &val) );

            if( !SCIPisInfinity(scip, -LUrhs[2*c]) )
               lhs = LUrhs[2*c] - val;
            else
               lhs = -SCIPinfinity(scip);

            if( !SCIPisInfinity(scip, LUrhs[2*c+1]) )
               rhs = LUrhs[2*c+1] - val;
            else
               rhs = SCIPinfinity(scip);

            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, con_name(c), 0, NULL, NULL, lhs, rhs) );

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

   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}

/*
 * Callback methods of reader
 */

#if 0
/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyNl)
{  /*lint --e{715}*/
   assert(scip != NULL);

   SCIP_CALL( SCIPincludeReaderNl(scip) );

   return SCIP_OKAY;
}
#endif

#if 0
/** destructor of reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreeNl)
{  /*lint --e{715}*/
   SCIP_READERDATA* readerdata;

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   SCIPfreeMemory(scip, &readerdata);

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

   assert(scip != NULL);
   assert(reader != NULL);
   assert(filename != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* initialize ASL */
   asl = ASL_alloc(ASL_read_fg);

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

   SCIP_CALL( setupObjective(scip, probdata, &success) );
   if( !success )
      return SCIP_READERROR;

   SCIP_CALL( setupConstraints(scip, probdata, &success) );
   if( !success )
      return SCIP_READERROR;

   *result = SCIP_SUCCESS;

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
   const char* msg;
   ASL* asl;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   if( probdata == NULL || probdata->asl == NULL )
   {
      SCIPerrorMessage("No AMPL problem read, cannot write AMPL solution then.\n");
      return SCIP_OKAY;
   }

   switch( SCIPgetStatus(scip) )
   {
      case SCIP_STATUS_UNKNOWN:
         msg = "unknown";
         break;
      case SCIP_STATUS_USERINTERRUPT:
         msg = "user interrupt";
         break;
      case SCIP_STATUS_NODELIMIT:
         msg = "node limit reached";
         break;
      case SCIP_STATUS_TOTALNODELIMIT:
         msg = "total node limit reached";
         break;
      case SCIP_STATUS_STALLNODELIMIT:
         msg = "stall node limit reached";
         break;
      case SCIP_STATUS_TIMELIMIT:
         msg = "time limit reached";
         break;
      case SCIP_STATUS_MEMLIMIT:
         msg = "memory limit reached";
         break;
      case SCIP_STATUS_GAPLIMIT:
         msg = "gap limit reached";
         break;
      case SCIP_STATUS_SOLLIMIT:
         msg = "solution limit reached";
         break;
      case SCIP_STATUS_BESTSOLLIMIT:
         msg = "solution improvement limit reached";
         break;
      case SCIP_STATUS_OPTIMAL:
         msg = "optimal solution found";
         break;
      case SCIP_STATUS_INFEASIBLE:
         msg = "infeasible";
         break;
      case SCIP_STATUS_UNBOUNDED:
         msg = "unbounded";
         break;
      case SCIP_STATUS_INFORUNBD:
         msg = "infeasible or unbounded";
         break;
      default:
         SCIPerrorMessage("invalid status code <%d>\n", SCIPgetStatus(scip));
         return SCIP_INVALIDDATA;
   }

   x = NULL;
   if( SCIPgetBestSol(scip) != NULL )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &x, probdata->nvars) );
      SCIP_CALL( SCIPgetSolVals(scip, SCIPgetBestSol(scip), probdata->nvars, probdata->vars, x) );
   }

   asl = probdata->asl;
   write_sol((char*)msg, x, NULL, NULL);

   SCIPfreeBufferArrayNull(scip, &x);

   return SCIP_OKAY;
}


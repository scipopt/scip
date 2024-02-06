/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlhdlr_signomial.c
 * @ingroup DEFPLUGINS_NLHDLR
 * @brief  signomial nonlinear handler
 * @author Liding Xu
 */

#ifdef SCIP_SIGCUT_DEBUG
#ifndef SCIP_SIG_DEBUG
#define SCIP_SIG_DEBUG
#endif
#endif

#ifdef SCIP_SIG_DEBUG
#define SCIP_DEBUG
#endif

#include "scip/nlhdlr_signomial.h"
#include "scip/cons_nonlinear.h"
#include "scip/pub_misc_rowprep.h"
#include "scip/pub_nlhdlr.h"
#include "scip/scip_expr.h"
#include "scip/expr_var.h"
#include "scip/expr_pow.h"

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME                    "signomial"
#define NLHDLR_DESC                    "handler for signomial expressions"
#define NLHDLR_DETECTPRIORITY          30
#define NLHDLR_ENFOPRIORITY            30

/* handler specific parameters */
#define NLHDLR_MAXNUNDERVARS           14    /**< maximum number of variables when underestimating a concave power function (maximum: 14) */
#define NLHDLR_MINCUTSCALE             1e-5  /**< minimum scale factor when scaling a cut (minimum: 1e-6) */


/** nonlinear handler expression data
 *
 * A signomial expression admits the form \f$ cx^a = y \f$, where \f$ y \f$ is an auxiliary variable representing this
 * expression and all \f$ x \f$ are positive.
 *
 * The natural formulation of the expression is defined as \f$ x^a = t = y/c \f$, where \f$ t \f$ is a
 * non-existant slack variable denoting \f$ y/c \f$. The variables in \f$x\f$ with positive exponents form positive
 * variables \f$ u \f$, and the associated exponents form positive exponents \f$ f \f$. The variables in \f$ x \f$ with
 * negative exponents and \f$ t \f$  form negative variables \f$ v \f$, and the associated exponents form negative
 * exponents \f$ g \f$. Let \f$ s =  \max(|f|,|g|) \f$ be a normalization constant, where \f$ |\cdot| \f$ denotes the L1 norm. Apply a scaling step:
 * Dividing the entries of \f$ f \f$  by \f$ s \f$, and dividing the entries of \f$ g \f$  by \f$ s \f$ as well. Then \f$ x^a = t \f$ has
 * a reformulation \f$ u^f = v^g \f$, where \f$ u^f, v^g \$ are two concave power functions.
 */
struct SCIP_NlhdlrExprData
{
   SCIP_Real             coef;               /**< coefficient \f$c\f$ */
   SCIP_EXPR**           factors;            /**< expression factors representing \f$x\f$ */
   int                   nfactors;           /**< number of factors */
   int                   nvars;              /**< number of variables \f$(x,y)\f$ */
   SCIP_Real*            exponents;          /**< exponents \f$e\f$ */
   int                   nposvars;           /**< number of positive variables \f$u\f$ */
   int                   nnegvars;           /**< number of negative variables  \f$v\f$ */
   SCIP_Bool*            signs;              /**< indicators for sign of variables after reformulation, TRUE for positive, FALSE for negative */
   SCIP_Real*            refexponents;       /**< exponents of \f$(x,t)\f$ after reformulation (always positive) */
   SCIP_Bool             isstorecapture;     /**< are all variables already got? */

   /* working parameters will be modified after getting all variables */
   SCIP_VAR**            vars;               /**< variables \f$(x,y)\f$ */
   SCIP_INTERVAL*        intervals;          /**< intervals storing lower and upper bounds of variables \f$(x,y)\f$ */
   SCIP_Real*            box;                /**< the upper/lower bounds of variables, used in SCIPcomputeFacetVertexPolyhedralNonlinear() */
   SCIP_Real*            xstar;              /**< the values of variables, used in SCIPcomputeFacetVertexPolyhedralNonlinear() */
   SCIP_Real*            facetcoefs;         /**< the coefficients of variables, returned by SCIPcomputeFacetVertexPolyhedralNonlinear() */
};

/** nonlinear handler data */
struct SCIP_NlhdlrData
{
   /* parameters */
   int                   maxnundervars;      /**< maximum number of variables in underestimating a concave power function */
   SCIP_Real             mincutscale;        /**< minimum scale factor when scaling a cut */
};

/** data struct to be passed on to vertexpoly-evalfunction (see SCIPcomputeFacetVertexPolyhedralNonlinear) */
typedef struct
{
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata;     /**< expression data */
   SCIP_Bool             sign;               /**< the sign of variables in the reformulated constraint for vertexpoly-evalfunction */
   int                   nsignvars;          /**< the number of variables in the reformulated constraint for vertexpoly-evalfunction */
   SCIP*                 scip;               /**< SCIP data structure */
} VERTEXPOLYFUN_EVALDATA;


/*
 * Local methods
 */

#ifdef SCIP_SIGCUT_DEBUG

/** print the information on a signomial term */
static
SCIP_RETCODE printSignomial(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata      /**< expression data */
   )
{
   assert(expr != NULL);

   /* print overall statistics and the expression */
   SCIPdebugMsg(scip, " #all variables: %d, #positive exponent variables: %d, #negative exponent variables: %d, auxvar: %s \n expr: ",
      nlhdlrexprdata->nvars, nlhdlrexprdata->nposvars, nlhdlrexprdata->nnegvars, SCIPvarGetName(SCIPgetExprAuxVarNonlinear(expr)) );
   SCIPprintExpr(scip, expr, NULL);

   /* if variables are detected, we can print more information */
   if( !nlhdlrexprdata->isstorecapture )
   {
      SCIPinfoMessage(scip, NULL, "\n");
      return SCIP_OKAY;
   }

   /* print the natural formulation of the expression */
   SCIPinfoMessage(scip, NULL, ". natural formulation c x^a = y: %.2f", nlhdlrexprdata->coef);
   for( int i = 0; i < nlhdlrexprdata->nvars - 1; i++ )
   {
      SCIPinfoMessage(scip, NULL, "%s^%.2f", SCIPvarGetName(nlhdlrexprdata->vars[i]), nlhdlrexprdata->exponents[i]);
   }
   SCIPinfoMessage(scip, NULL, " = %s", SCIPvarGetName(nlhdlrexprdata->vars[nlhdlrexprdata->nvars - 1]));

   /* print the reformulation of the expression */
   SCIPinfoMessage(scip, NULL, ". Reformulation u^f = v^g: ");
   if( nlhdlrexprdata->nposvars == 0 )
   {
      SCIPinfoMessage(scip, NULL, "1");
   }
   else
   {
      for( int i = 0; i < nlhdlrexprdata->nvars; i++ )
         if( nlhdlrexprdata->signs[i] )
         {
            SCIPinfoMessage(scip, NULL, "%s^%.2f", SCIPvarGetName(nlhdlrexprdata->vars[i]), nlhdlrexprdata->refexponents[i]);
         }
   }
   SCIPinfoMessage(scip, NULL, " = ");
   if( nlhdlrexprdata->nnegvars == 0 )
   {
      SCIPinfoMessage(scip, NULL, "1");
   }
   else
   {
      for( int i = 0; i < nlhdlrexprdata->nvars; i++ )
      {
         if( !nlhdlrexprdata->signs[i] )
         {
            if( i == nlhdlrexprdata->nvars - 1 )
            {
               SCIPinfoMessage(scip, NULL, "(%s/%.2f)^%.2f", SCIPvarGetName(nlhdlrexprdata->vars[i]), nlhdlrexprdata->coef, nlhdlrexprdata->refexponents[i]);
            }
            else
            {
               SCIPinfoMessage(scip, NULL, "%s^%.2f", SCIPvarGetName(nlhdlrexprdata->vars[i]), nlhdlrexprdata->refexponents[i]);
            }
         }
      }
   }
   SCIPinfoMessage(scip, NULL, "\n");

   return SCIP_OKAY;
}

#endif

/** free the memory of expression data */
static
void freeExprDataMem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata,     /**< expression data */
   SCIP_Bool             ispartial           /**< free the partially allocated memory or the fully allocated memory? */
   )
{

   SCIPfreeBlockMemoryArrayNull(scip, &(*nlhdlrexprdata)->factors, (*nlhdlrexprdata)->nfactors);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->exponents, (*nlhdlrexprdata)->nfactors);
   if( !ispartial )
   {
      int nvars = (*nlhdlrexprdata)->nvars;
      SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->signs, nvars);
      SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->refexponents, nvars);
      SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->vars, nvars);
      SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->intervals, nvars);
      SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->box, 2*nvars);  /*lint !e647*/
      SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->xstar, nvars);
      SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->facetcoefs, nvars);
   }
   SCIPfreeBlockMemory(scip, nlhdlrexprdata);
   *nlhdlrexprdata = NULL;
}

/** reforms a rowprep to a standard form for nonlinear handlers
 *
 * The input rowprep is of the form  \f$ a_u u + a_v v + b \f$.
 * If in the overestimate mode, then we relax \f$ t \le x^a \f$, i.e., \f$ 0 \le u^f - v^g \f$. This implies that \f$ t \f$ is in \f$ v = (v',t) \f$.
 * Therefore, the valid inequality is \f$ 0 \le a_u u + a_v v + b \f$. As overestimate mode requires that \f$ t \f$ is in the left hand side,
 * the coefficients of \f$ t \f$ must be made negative while keeping the sign of the inequality, we can show that \f$ a_t \le 0 \f$, so it suffices
 * to divide the both sides by \f$ -a_t \ge 0\f$, which yields  \f$ t \le (a_u u + a_v' v' + b) / -a_t \f$.
 * A rowprep in standard form only contains an estimator of the expression and no auxvar.
 * If in the underestimate mode, then we relax \f$ x^a \le t \f$, i.e., \f$ u^f - v^g \le 0 \f$. This implies that \f$ t \f$ is in \f$ v = (v',t) \f$.
 * Therefore, the valid inequality is \f$ a_u u + a_v v + b \le 0 \f$. As overestimate mode requires that \f$ t \f$ is in the left hand side,
 * the coefficients of \f$ t \f$ must be made negative while keeping the sign of the inequality, we can show that \f$ a_t \le 0 \f$, so it suffices
 * to divide the both sides by \f$ -a_t \ge 0 \f$, which yields \f$ (a_u u + a_v' v' + b) / -a_t \le t \f$.
 * A rowprep in standard form only contains an estimator of the expression and no auxvar.
 */
static
void reformRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< expression data */
   SCIP_ROWPREP*         rowprep,            /**< cut to be reformulated */
   SCIP_Real             mincutscale,        /**< min scaling factor for the cut in rowprep */
   SCIP_Bool*            success             /**< pointer to store whether the reformulating was successful */
)
{
   int i;
   int nvars;
   SCIP_Real coefauxvar;
   SCIP_Real scale;
   SCIP_Real* coefs;
   SCIP_VAR* auxvar;
   SCIP_VAR** vars;

   assert(rowprep != NULL);
   assert(nlhdlrexprdata != NULL);

   coefs = SCIProwprepGetCoefs(rowprep);
   vars = SCIProwprepGetVars(rowprep);
   nvars = SCIProwprepGetNVars(rowprep);

   /* find auxvar's cut coefficient and set it to zero */
   auxvar = nlhdlrexprdata->vars[nlhdlrexprdata->nfactors];
   coefauxvar = 1.0;
   for( i = 0; i < nvars; i++ )
   {
      if( vars[i] == auxvar )
      {
         coefauxvar = coefs[i];
         coefs[i] = 0.0;
         break;
      }
   }

   if( SCIPisZero(scip, coefauxvar) || coefauxvar >= 0.0 )
   {
      *success = FALSE;
   }
   else
   {
      /* the reformation scales the cut so that coefficients and constant are divided by the absolute value of coefauxvar */
      assert(coefauxvar < 0.0);
      scale = -1.0 / coefauxvar;
      for( i = 0; i < nvars; i++ )
      {
         if( vars[i] == auxvar )
            continue;
         coefs[i] *= scale;
      }
      /* set the side to scale*side by adding -side and adding scale*side */
      SCIProwprepAddSide(rowprep, SCIProwprepGetSide(rowprep) * (-1.0 + scale));

      *success = SCIPisGT(scip, scale, mincutscale);
   }
}

/** store and capture variables associated with the expression and its subexpressions */
static
SCIP_RETCODE storeCaptureVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata      /**< expression data */
   )
{
   int c;

   assert(!nlhdlrexprdata->isstorecapture);

   /* get and capture variables \f$x\f$ */
   for( c = 0; c < nlhdlrexprdata->nfactors; ++c )
   {
      nlhdlrexprdata->vars[c] = SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->factors[c]);
      assert(nlhdlrexprdata->vars[c] != NULL);

      SCIP_CALL( SCIPcaptureVar(scip, nlhdlrexprdata->vars[c]) );
   }

   /* get and capture variable \f$y\f$ */
   nlhdlrexprdata->vars[c] = SCIPgetExprAuxVarNonlinear(expr);
   assert(nlhdlrexprdata->vars[c] != NULL);
   SCIP_CALL( SCIPcaptureVar(scip, nlhdlrexprdata->vars[c]) );

   nlhdlrexprdata->isstorecapture = TRUE;

   return SCIP_OKAY;
}

/** get bounds of variables x,t and check whether they are box constrained signomial variables */
static
void checkSignomialBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< expression data */
   SCIP_Bool*            isboxsignomial      /**< buffer to store whether variables are box constrained signomial variables */
   )
{
   int c;
   SCIP_Real powinf;
   SCIP_Real powsup;
   SCIP_Real productinf = 1;
   SCIP_Real productsup = 1;

   assert(nlhdlrexprdata->isstorecapture);

   *isboxsignomial = FALSE;

   /* get bounds of x */
   for( c = 0; c < nlhdlrexprdata->nfactors; c++ )
   {
      SCIP_Real inf = SCIPvarGetLbLocal(nlhdlrexprdata->vars[c]);
      SCIP_Real sup = SCIPvarGetUbLocal(nlhdlrexprdata->vars[c]);

      /* if the bounds of the variable are not positive and finite, or (bounds equal) then the expression is not a signomial */
      if( !SCIPisPositive(scip, inf) || !SCIPisPositive(scip, sup) || SCIPisInfinity(scip, sup) || SCIPisEQ(scip, inf, sup) )
         return;

      nlhdlrexprdata->intervals[c].inf = inf;
      nlhdlrexprdata->intervals[c].sup = MAX(sup, inf + 0.1);
      powinf = pow(inf, nlhdlrexprdata->exponents[c]);
      powsup = pow(sup, nlhdlrexprdata->exponents[c]);
      productinf *= MIN(powinf, powsup);
      productsup *= MAX(powinf, powsup);
   }

   /* compute bounds of t by bounds of x */
   nlhdlrexprdata->intervals[c].inf = productinf;
   nlhdlrexprdata->intervals[c].sup = productsup;

   *isboxsignomial = TRUE;
}

/** evaluate expression at solution w.r.t. auxiliary variables */
static
SCIP_DECL_VERTEXPOLYFUN(nlhdlrExprEvalPower)
{
   int i;
   int j;
   SCIP_Real val;
   VERTEXPOLYFUN_EVALDATA* evaldata;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata;

   assert(args != NULL);

   evaldata = (VERTEXPOLYFUN_EVALDATA *)funcdata;

   assert(evaldata != NULL);
   assert(nargs == evaldata->nsignvars);

   nlhdlrexprdata = evaldata->nlhdlrexprdata;

#ifdef SCIP_MORE_DEBUG
   SCIPdebugMsg(evaldata->scip, "eval vertexpolyfun\n");
#endif
   val = 1.0;
   for( i = 0, j = 0; i < nlhdlrexprdata->nvars; ++i )
   {
      if( nlhdlrexprdata->signs[i] != evaldata->sign )
         continue;
      /* the reformulated exponent of args[j] is found */
      val *= pow(args[j], nlhdlrexprdata->refexponents[i]);
      j++;
   }

   return val;
}

/** determine whether a power function \f$ w^h \f$ is special and add an overunderestimator or underestimator to a given rowprep
 *
 * \f$ w^h \f$ is special, if all variables are fixed, or it is a constant to estimate, a univariate power to estimate,
 * or a bivariate power to underestimate. The estimator is multiplied by the multiplier and stored in the rowprep.
 */
static
SCIP_RETCODE estimateSpecialPower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nonlinear handler expression data */
   SCIP_Bool             sign,               /**< the sign of variables of the power function */
   SCIP_Real             multiplier,         /**< the multiplier of the estimator */
   SCIP_Bool             overestimate,       /**< whether overestimate or underestimator the power function */
   SCIP_SOL*             sol,                /**< solution \f$ w' \f$  to use in estimation */
   SCIP_ROWPREP*         rowprep,            /**< rowprep where to store estimator */
   SCIP_Bool*            isspecial,          /**< buffer to store whether this function is a special case */
   SCIP_Bool*            success             /**< buffer to store whether successful */
   )
{
   int i;
   SCIP_Bool allfixed = TRUE;
   int nsignvars;

   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(rowprep != NULL);
   assert(isspecial != NULL);
   assert(success != NULL);

   *success = FALSE;
   /* check whether all variables are fixed */
   for( i = 0; i < nlhdlrexprdata->nvars; ++i )
   {
      if( nlhdlrexprdata->signs[i] != sign )
         continue;

      if( !SCIPisRelEQ(scip, nlhdlrexprdata->intervals[i].inf, nlhdlrexprdata->intervals[i].sup) )
      {
         allfixed = FALSE;
         break;
      }
   }

   /* allfixed is a special case */
   if( allfixed )
   {
      /* return a constant */
      SCIP_Real funcval = 1.0;
      SCIP_Real scale;
      SCIP_Real val;

      for( i = 0; i < nlhdlrexprdata->nvars; ++i )
      {
         if( nlhdlrexprdata->signs[i] != sign )
            continue;

         scale = i == (nlhdlrexprdata->nvars - 1) ? nlhdlrexprdata->coef : 1.0;
         val = SCIPgetSolVal(scip, sol, nlhdlrexprdata->vars[i]) / scale;
         funcval *= pow(val, nlhdlrexprdata->refexponents[i]);
      }
      SCIProwprepAddConstant(rowprep, multiplier * funcval);
      *isspecial = TRUE;

      return SCIP_OKAY;
   }

   /* number of variables in the power function */
   nsignvars = sign ? nlhdlrexprdata->nposvars : nlhdlrexprdata->nnegvars;

   /* if the power function has no more than 2 variables, this a special case */
   *isspecial = ( nsignvars <= 1 ) || ( nsignvars == 2  && !overestimate );
   if( !*isspecial )
      return SCIP_OKAY;

   if( nsignvars == 0 )
   {
      /* constant case */
      SCIProwprepAddConstant(rowprep, multiplier);
      *success = TRUE;
   }
   else if( nsignvars == 1 )
   {
      /* univariate case, w^h */
      for( i = 0; i < nlhdlrexprdata->nvars; ++i )
      {
         if( nlhdlrexprdata->signs[i] == sign )
         {
            /* the variable w is found */
            SCIP_Real scale = i == (nlhdlrexprdata->nvars - 1) ? nlhdlrexprdata->coef : 1;
            SCIP_VAR* var = nlhdlrexprdata->vars[i];
            SCIP_Real refexponent = nlhdlrexprdata->refexponents[i];
            if( refexponent == 1.0 )
            {
               /* h = 1, a univariate linear function. Only rescale, no need for linearization */
               SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, multiplier / scale) );
            }
            else
            {
               /* a univariate power function */
               SCIP_Real facetconstant;
               SCIP_Real facetcoef;
               SCIP_Real val = SCIPgetSolVal(scip, sol, var) / scale;
               /* local (using bounds) depends on whether we under- or overestimate */
               SCIP_Bool islocal = !overestimate;
               SCIPestimateRoot(scip, refexponent, overestimate, nlhdlrexprdata->intervals[i].inf, nlhdlrexprdata->intervals[i].sup,
                  val, &facetconstant, &facetcoef, &islocal, success);
               SCIProwprepAddConstant(rowprep,  multiplier * facetconstant);
               SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, multiplier * facetcoef / scale) );
            }
         }
      }
      *success = TRUE;
   }
   else if( nsignvars == 2 && !overestimate )
   {
      /* bivariate case, f(w) = w^h = f_0(w_0) f_1(w_1). The convex under envelope is the maxmimum of
       * two affine functions, each of which is determined by three extreme points the box bound.
       * One affine function is supported by three lower left points; the other affine function is
       * supported by three upper right points. For a point xstar in the box, its corresponding affine function can be determined by
       * which region (upper right or lower left half space) the point is in. Thus, we can determine the region, and use the
       * associated three extreme point to interpolate the affine function.
       */
      SCIP_Bool isupperright;
      SCIP_Real dw0, dw1;
      SCIP_Real f0w0l, f0w0u, f1w1l, f1w1u;
      SCIP_Real fw0lw1u, fw0uw1l;
      SCIP_Real facetconstant;
      SCIP_Real facetcoefs[2] = {0.0, 0.0};
      SCIP_VAR* vars[2] = {NULL, NULL};
      SCIP_Real refexponents[2] = {0.0, 0.0};
      SCIP_Real xstar[2] = {0.0, 0.0};;
      SCIP_Real scale[2] = {0.0, 0.0};;
      SCIP_Real box[4] = {0.0, 0.0, 0.0, 0.0};
      int j = 0;

      /* get refexponents, xstar, scale, and box */
      for( i = 0; i < nlhdlrexprdata->nvars; ++i )
      {
         if( nlhdlrexprdata->signs[i] != sign )
            continue;
         box[2 * j] = nlhdlrexprdata->intervals[i].inf;
         box[2 * j + 1] = nlhdlrexprdata->intervals[i].sup;
         refexponents[j] = nlhdlrexprdata->refexponents[i];
         scale[j] = i == (nlhdlrexprdata->nvars - 1) ? nlhdlrexprdata->coef : 1;
         vars[j] = nlhdlrexprdata->vars[i];
         xstar[j] = SCIPgetSolVal(scip, sol, nlhdlrexprdata->vars[j]) / scale[j];
         j++;
      }

      /* compute the box length*/
      dw0 = box[1] - box[0];
      dw1 = box[3] - box[2];

      /* determine the location (upper right or lower left half sapce) of the xstar.
       * (dw1, dw0) is the direction vector to the upper right half space.
       */
      isupperright = ( (xstar[0] - box[0]) * dw1 + (xstar[1] - box[3]) * dw0 ) > 0.0;

      /* compute function values of f_0, f_1 at vertices */
      f0w0l = pow(box[0], refexponents[0]);
      f0w0u = pow(box[1], refexponents[0]);
      f1w1l = pow(box[2], refexponents[1]);
      f1w1u = pow(box[3], refexponents[1]);
      fw0lw1u = f0w0l * f1w1u;
      fw0uw1l = f0w0u * f1w1l;
      if( isupperright )
      {
         /* underestimator: fw0uw1u + (fw0uw1u - fw0lw1u) / (dw0) * (x0 - w0u) + (fw0uw1u - fw0uw1l) / (dw1) * (x1 - w1u) */
         SCIP_Real fw0uw1u = f0w0u * f1w1u;
         facetcoefs[0] = (fw0uw1u - fw0lw1u) / dw0;
         facetcoefs[1] = (fw0uw1u - fw0uw1l) / dw1;
         facetconstant = fw0uw1u  - facetcoefs[0] * box[1] - facetcoefs[1] * box[3];
      }
      else
      {
         /* underestimator: fw0lw1l + (fw0uw1l - fw0lw1l) / (dw0) * (x0 - w0l) + (fw0lw1u- fw0lw1l) / (dw1) * (x1 - w1l) */
         SCIP_Real fw0lw1l = f0w0l * f1w1l;
         facetcoefs[0] = (fw0uw1l - fw0lw1l) / dw0;
         facetcoefs[1] = (fw0lw1u - fw0lw1l) / dw1;
         facetconstant = fw0lw1l  - facetcoefs[0] * box[0] - facetcoefs[1] * box[2];
      }
      SCIProwprepAddConstant(rowprep,  multiplier * facetconstant);
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, vars[0], multiplier * facetcoefs[0] / scale[0]) );
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, vars[1], multiplier * facetcoefs[1] / scale[1]) );
      *success = TRUE;
   }

   return SCIP_OKAY;
}

/** adds an underestimator for a multivariate concave power function \f$ w^h \f$ to a given rowprep
 *
 * Calls \ref SCIPcomputeFacetVertexPolyhedralNonlinear() for \f$ w^h \f$  and
 * box set to local bounds of auxiliary variables. The estimator is multiplied
 * by the multiplier and stored in the rowprep.
 */
static
SCIP_RETCODE underEstimatePower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nonlinear handler expression data */
   SCIP_Bool             sign,               /**< the sign of variables of the power function */
   SCIP_Real             multiplier,         /**< the multiplier of the estimator */
   SCIP_SOL*             sol,                /**< solution \f$ w' \f$ to use*/
   SCIP_Real             targetvalue,        /**< a target value to achieve; if not reachable, then can give up early */
   SCIP_ROWPREP*         rowprep,            /**< rowprep where to store estimator */
   SCIP_Bool*            success             /**< buffer to store whether successful */
   )
{
   int i;
   int j;
   int nsignvars;
   SCIP_Real facetconstant;
   SCIP_Real scale;
   SCIP_Real* box;
   SCIP_Real* facetcoefs;
   SCIP_Real* xstar;
   VERTEXPOLYFUN_EVALDATA evaldata;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(rowprep != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* number of variables of sign */
   nsignvars = sign ? nlhdlrexprdata->nposvars : nlhdlrexprdata->nnegvars;

   /* data structure to evaluate the power function */
   evaldata.nlhdlrexprdata = nlhdlrexprdata;
   evaldata.sign = sign;
   evaldata.nsignvars = nsignvars;
   evaldata.scip = scip;

   /* compute box constraints and reference value of variables*/
   xstar = nlhdlrexprdata->xstar;
   box = nlhdlrexprdata->box;
   j = 0;
   for( i = 0; i < nlhdlrexprdata->nvars; ++i )
   {
      if( nlhdlrexprdata->signs[i] != sign )
         continue;

      box[2 * j] = nlhdlrexprdata->intervals[i].inf;
      box[2 * j + 1] = nlhdlrexprdata->intervals[i].sup;
      scale = i == (nlhdlrexprdata->nvars - 1) ? nlhdlrexprdata->coef : 1.0;
      xstar[j] = SCIPgetSolVal(scip, sol, nlhdlrexprdata->vars[i]) / scale;
      j++;
   }

   /* find a facet of the underestimator */
   facetcoefs = nlhdlrexprdata->facetcoefs;
   SCIP_CALL( SCIPcomputeFacetVertexPolyhedralNonlinear(scip, conshdlr, FALSE, nlhdlrExprEvalPower, (void*)&evaldata, xstar, box,
      nsignvars, targetvalue, success, facetcoefs, &facetconstant) );

   if( !*success )
   {
      SCIPdebugMsg(scip, "failed to compute facet of convex hull\n");
      return SCIP_OKAY;
   }

   /* set coefficients in the rowprep */
   SCIProwprepAddConstant(rowprep, multiplier * facetconstant);
   j = 0;
   for( i = 0; i < nlhdlrexprdata->nvars; ++i )
   {
      if( nlhdlrexprdata->signs[i] != sign )
         continue;
      scale = i == (nlhdlrexprdata->nvars - 1) ? nlhdlrexprdata->coef : 1.0;
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, nlhdlrexprdata->vars[i], multiplier * facetcoefs[j] / scale) );
      j++;
   }

   return SCIP_OKAY;
}

/** adds an overestimator for a concave power function \f$ w^h \f$ to a given rowprep
 *
 * The estimator is multiplied by the multiplier and stored in the rowprep.
 */
static
SCIP_RETCODE overEstimatePower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nonlinear handler expression data */
   SCIP_Bool             sign,               /**< the sign of variables of the power function */
   SCIP_Real             multiplier,         /**< the multiplier of the estimator */
   SCIP_SOL*             sol,                /**< solution to use */
   SCIP_ROWPREP*         rowprep,            /**< rowprep where to store estimator */
   SCIP_Bool*            success             /**< buffer to store whether successful */
   )
{
   int i;
   SCIP_Real facetcoef;
   SCIP_Real facetconstant;
   SCIP_Real funcval;
   SCIP_Real refexponent;
   SCIP_Real sumrefexponents;
   SCIP_Real scale;
   SCIP_Real val;
   SCIP_VAR* var;
   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(rowprep != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* compute the value and the sum of reformulated exponents of the power function */
   sumrefexponents = 0;
   funcval = 1;
   for( i = 0; i < nlhdlrexprdata->nvars; ++i )
   {
      if( nlhdlrexprdata->signs[i] != sign )
         continue;
      scale = i == (nlhdlrexprdata->nvars - 1) ? nlhdlrexprdata->coef : 1.0;
      val = SCIPgetSolVal(scip, sol, nlhdlrexprdata->vars[i]) / scale;
      val = MAX(val, 0.1);
      refexponent = nlhdlrexprdata->refexponents[i];
      sumrefexponents += refexponent;
      funcval *= pow(val, refexponent);
   }

   /* overestimate by gradient cut: w'^h + h w'^(h - 1)(w - w') */
   facetconstant = (1 - sumrefexponents) * funcval;
   SCIProwprepAddConstant(rowprep,  multiplier * facetconstant);
   for( i = 0; i < nlhdlrexprdata->nvars; ++i )
   {
      if( nlhdlrexprdata->signs[i] != sign )
         continue;
      scale = i == (nlhdlrexprdata->nvars - 1) ? nlhdlrexprdata->coef : 1.0;
      var = nlhdlrexprdata->vars[i];
      val = SCIPgetSolVal(scip, sol, var) / scale;
      val = MAX(val, 0.1);
      facetcoef = nlhdlrexprdata->refexponents[i] * funcval / val;
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, multiplier * facetcoef / scale) );
   }

   *success = TRUE;

   return SCIP_OKAY;
}

/*
 * Callback methods of nonlinear handler
 */

/** nonlinear handler estimation callback */
static
SCIP_DECL_NLHDLRESTIMATE(nlhdlrEstimateSignomial)
{ /*lint --e{715}*/
   SCIP_Bool isboxsignomial;
   SCIP_Bool isspecial;
   SCIP_Bool undersign;
   SCIP_Bool oversign;
   int i;
   int nundervars;
   SCIP_Real undermultiplier;
   SCIP_Real overmultiplier;
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_ROWPREP* rowprep;
   SCIP_Real targetunder;
   SCIP_Real scale;

   assert(conshdlr != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(rowpreps != NULL);
   *success = FALSE;
   *addedbranchscores = FALSE;

   /* store and capture the vars of an expression, if the vars are not stored and captured yet */
   if( !nlhdlrexprdata->isstorecapture )
   {
      SCIP_CALL( storeCaptureVars(scip, expr, nlhdlrexprdata) );
   }

   /* check whether all variables have finite positive bounds, which is necessary for the expression to be a signomial term */
   /* TODO consider allowing 0.0 lower bounds for u variables (and handle gradients at 0.0) */
   isboxsignomial = FALSE;
   checkSignomialBounds(scip, nlhdlrexprdata, &isboxsignomial);

   if( !isboxsignomial )
      return SCIP_OKAY;

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);

   /* determine the sign of variables for over- and underestimators, and the multiplier for estimators in the rowprep */
   if( overestimate )
   {
      /* relax t <= x^a, i.e., 0 <= u^f - v^g, overestimate u^f, underestimate v^g */
      oversign = TRUE;
      undersign = FALSE;
      overmultiplier = 1.0;
      undermultiplier = -1.0;
      nundervars = nlhdlrexprdata->nnegvars;
   }
   else
   {
      /* relax x^a <= t, i.e., u^f - v^g <= 0, underestimate u^f, overestimate v^g */
      oversign = FALSE;
      undersign = TRUE;
      overmultiplier = -1.0;
      undermultiplier = 1.0;
      nundervars = nlhdlrexprdata->nposvars;
   }

   /* compute the value of the overestimator, which is a target value for computing the underestimator efficiently */
   targetunder = 1.0;
   for( i = 0; i < nlhdlrexprdata->nvars; i++ )
   {
      if( nlhdlrexprdata->signs[i] == oversign )
      {
         scale = i == (nlhdlrexprdata->nvars - 1) ? nlhdlrexprdata->coef : 1.0;
         targetunder *= pow(SCIPgetSolVal(scip, sol, nlhdlrexprdata->vars[i]) / scale, nlhdlrexprdata->refexponents[i]);
      }
   }

#ifdef SCIP_SIGCUT_DEBUG
   SCIP_Real targetover = 1.0;

   /* print information on estimators */
   SCIP_CALL( printSignomial(scip, expr, nlhdlrexprdata));
   SCIPinfoMessage(scip, NULL, " Auxvalue: %f, targetvalue: %f, %sestimate.", auxvalue, targetvalue, overestimate ? "over" : "under");

   targetunder = 1.0;
   for( i = 0; i < nlhdlrexprdata->nvars; i++ )
   {
      if( nlhdlrexprdata->signs[i] == undersign )
      {
         scale = i == (nlhdlrexprdata->nvars - 1) ? nlhdlrexprdata->coef : 1.0;
         targetover *= pow(SCIPgetSolVal(scip, sol, nlhdlrexprdata->vars[i]) / scale, nlhdlrexprdata->refexponents[i]);
      }
   }
   SCIPinfoMessage(scip, NULL, " Underestimate: %s, overestimate: %s.",
      undersign ? "positive" : "negative", oversign ? "positive" : "negative");
   SCIPinfoMessage(scip, NULL, " Undervalue (targetover): %f, overvalue (targetunder): %f.", targetover, targetunder);
#endif

   /* create a rowprep and allocate memory for it */
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, TRUE) );
   SCIP_CALL( SCIPensureRowprepSize(scip, rowprep, nlhdlrexprdata->nvars + 1) );

   /* only underestimate a concave function, if the number of variables is below a given threshold */
   if( nundervars <= nlhdlrdata->maxnundervars )
   {
      /* compute underestimator */
      isspecial = FALSE;
      /* check and compute the special case */
      SCIP_CALL( estimateSpecialPower(scip, nlhdlrexprdata, undersign, undermultiplier, FALSE, sol, rowprep, &isspecial, success) );
      if( !isspecial )
      {
         SCIP_CALL( underEstimatePower(scip, conshdlr, nlhdlr, nlhdlrexprdata, undersign, undermultiplier, sol, targetunder, rowprep, success) );
      }
   }

   if( *success )
   {
      /* compute overestimator */
      isspecial = FALSE;
      /* check and compute the special case */
      SCIP_CALL( estimateSpecialPower(scip, nlhdlrexprdata, oversign, overmultiplier, TRUE, sol, rowprep, &isspecial, success) );
      if( !isspecial )
      {
         SCIP_CALL( overEstimatePower(scip, nlhdlrexprdata, oversign, overmultiplier, sol, rowprep, success) );
      }
   }

   if( *success )
   {
      reformRowprep(scip, nlhdlrexprdata, rowprep, nlhdlrdata->mincutscale, success);
      if( *success )
      {
         SCIP_CALL( SCIPsetPtrarrayVal(scip, rowpreps, 0, rowprep) );
      }
   }

   if( !*success )
      SCIPfreeRowprep(scip, &rowprep);

   return SCIP_OKAY;
}

/** callback to detect structure in expression tree */
static
SCIP_DECL_NLHDLRDETECT(nlhdlrDetectSignomial)
{ /*lint --e{715}*/

   assert(expr != NULL);
   assert(enforcing != NULL);
   assert(participating != NULL);

   /* for now, we do not get active if separation is already (or does not need to be) provided */
   if( (*enforcing & SCIP_NLHDLR_METHOD_SEPABOTH) == SCIP_NLHDLR_METHOD_SEPABOTH )
   {
      return SCIP_OKAY;
   }

   /* check for product expressions with more than one child */
   if( SCIPisExprProduct(scip, expr) && SCIPexprGetNChildren(expr) >= 2 )
   {
      int c;
      int nf = SCIPexprGetNChildren(expr);
      int nvars = nf + 1;
      SCIP_Bool ismultilinear = TRUE;

      /* create expression data for the nonlinear handler */
      SCIP_CALL( SCIPallocClearBlockMemory(scip, nlhdlrexprdata) );
      (*nlhdlrexprdata)->nfactors = nf;
      (*nlhdlrexprdata)->nvars = nvars;

      /* allocate memory for expression data */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->factors, nf) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->exponents, nf) );

      /* get monomial information */
      SCIP_CALL( SCIPgetExprMonomialData(scip, expr, &((*nlhdlrexprdata)->coef), (*nlhdlrexprdata)->exponents, (*nlhdlrexprdata)->factors) );

      /* skip multilinear terms, since we wouldn't do better than expr_product */
      for( c = 0; c < nf; c++ )
      {
         if( !SCIPisEQ(scip, (*nlhdlrexprdata)->exponents[c], 1.0) )
         {
            ismultilinear = FALSE;
            break;
         }
      }

      if( ismultilinear )
      {
         /* if multilinear, we free memory of the expression data and do nothing */
         freeExprDataMem(scip, nlhdlrexprdata, TRUE);
         return SCIP_OKAY;
      }
      else
      {
         SCIP_Real normalize;
         SCIP_Real sumlexponents = 0;
         SCIP_Real sumrexponents = 1;
         int nposvars = 0;

         /* allocate more memory for expression data */
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->signs, nvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->refexponents, nvars) );
         SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(*nlhdlrexprdata)->vars, nvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->intervals, nvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->xstar, nvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->facetcoefs, nvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->box, 2 * nvars) );

         (*nlhdlrexprdata)->isstorecapture = FALSE;

         /* detect more information for the reformulation: we first compute the sum
          * of positive and negative exponents and determine the sign indicators
          */
         for( c = 0; c < nf; c++ )
         {
            /* capture sub expressions */
            SCIPcaptureExpr((*nlhdlrexprdata)->factors[c]);
            if( (*nlhdlrexprdata)->exponents[c] > 0.0 )
            {
               /* add a positive exponent */
               sumlexponents += (*nlhdlrexprdata)->exponents[c];
               (*nlhdlrexprdata)->signs[c] = TRUE;
               nposvars++;
            }
            else
            {
               /* subtract a negative exponent */
               sumrexponents -= (*nlhdlrexprdata)->exponents[c];
               (*nlhdlrexprdata)->signs[c] = FALSE;
            }
            /* set null to working variables, meaning that they are not stored yet */
         }
         (*nlhdlrexprdata)->signs[nf] = FALSE;
         (*nlhdlrexprdata)->nposvars = nposvars;
         (*nlhdlrexprdata)->nnegvars = nf - nposvars + 1;

         /* compute the normalization constant */
         normalize = MAX(sumlexponents, sumrexponents);
         /* normalize positive and negative exponents */
         for( c = 0; c < nf; c++ )
         {
            if( (*nlhdlrexprdata)->signs[c] )
               (*nlhdlrexprdata)->refexponents[c] = (*nlhdlrexprdata)->exponents[c] / normalize;
            else
               (*nlhdlrexprdata)->refexponents[c] = -(*nlhdlrexprdata)->exponents[c] / normalize;
         }
         (*nlhdlrexprdata)->refexponents[nf] = 1.0 / normalize;

         /* tell children that we will use their auxvar and use its activity for estimation */
         for( c = 0; c < nf; c++ )
         {
            SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, (*nlhdlrexprdata)->factors[c], TRUE, FALSE, TRUE, TRUE) );
         }
      }
   }

   if( *nlhdlrexprdata != NULL )
   {
      *participating = SCIP_NLHDLR_METHOD_SEPABOTH;
   }

#ifdef SCIP_SIGCUT_DEBUG
   if( *participating )
   {
      SCIPdebugMsg(scip, "scip depth: %d, step: %d, expr pointer: %p, expr data pointer: %p, detected expr: total vars (exps) %d ",
         SCIPgetSubscipDepth(scip), SCIPgetStage(scip), (void *)expr, (void *)nlhdlrexprdata, (*nlhdlrexprdata)->nfactors);
      SCIPprintExpr(scip, expr, NULL);
      SCIPinfoMessage(scip, NULL, " participating: %d\n", *participating);
   }
#endif

   return SCIP_OKAY;
}

/** auxiliary evaluation callback of nonlinear handler */
static SCIP_DECL_NLHDLREVALAUX(nlhdlrEvalauxSignomial)
{ /*lint --e{715}*/
   int c;
   SCIP_Real val;
   SCIP_VAR* var;

   *auxvalue = nlhdlrexprdata->coef;
   for( c = 0; c < nlhdlrexprdata->nfactors; ++c )
   {
      assert(nlhdlrexprdata->factors[c] != NULL);

      var = SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->factors[c]);
      assert(var != NULL);

      val = SCIPgetSolVal(scip, sol, var);
      if( SCIPisPositive(scip, val) )
      {
         *auxvalue *= pow(val, nlhdlrexprdata->exponents[c]);
      }
      else
      {
         *auxvalue = SCIP_INVALID;
         break;
      }
   }

   return SCIP_OKAY;
}

/** nonlinear handler copy callback */
static
SCIP_DECL_NLHDLRCOPYHDLR(nlhdlrCopyhdlrSignomial)
{ /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPnlhdlrGetName(sourcenlhdlr), NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeNlhdlrSignomial(targetscip) );

   return SCIP_OKAY;
}

/** callback to free data of handler */
static
SCIP_DECL_NLHDLRFREEHDLRDATA(nlhdlrFreehdlrDataSignomial)
{ /*lint --e{715}*/
   assert(nlhdlrdata != NULL);

   SCIPfreeBlockMemory(scip, nlhdlrdata);

   return SCIP_OKAY;
}

/** callback to free expression specific data */
static
SCIP_DECL_NLHDLRFREEEXPRDATA(nlhdlrFreeExprDataSignomial)
{ /*lint --e{715}*/
   int c;

   /* release expressions */
   for( c = 0; c < (*nlhdlrexprdata)->nfactors; c++ )
   {
      SCIP_CALL( SCIPreleaseExpr(scip, &(*nlhdlrexprdata)->factors[c]) );
   }

   /* release variables */
   if( (*nlhdlrexprdata)->isstorecapture )
   {
      for( c = 0; c < (*nlhdlrexprdata)->nvars; c++ )
      {
         if( (*nlhdlrexprdata)->vars[c] != NULL )
         {
            SCIP_CALL( SCIPreleaseVar(scip, &(*nlhdlrexprdata)->vars[c]) );
         }
      }
   }

   /* free memory */
   freeExprDataMem(scip, nlhdlrexprdata, FALSE);

   return SCIP_OKAY;
}


/*
 * nonlinear handler specific interface methods
 */

/** includes signomial nonlinear handler in nonlinear constraint handler */
SCIP_RETCODE
SCIPincludeNlhdlrSignomial(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_NLHDLR* nlhdlr;

   assert(scip != NULL);

   /* create nonlinear handler specific data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &nlhdlrdata));
   BMSclearMemory(nlhdlrdata);

   SCIP_CALL( SCIPincludeNlhdlrNonlinear(scip, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_DETECTPRIORITY,
      NLHDLR_ENFOPRIORITY, nlhdlrDetectSignomial, nlhdlrEvalauxSignomial, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIPnlhdlrSetCopyHdlr(nlhdlr, nlhdlrCopyhdlrSignomial);
   SCIPnlhdlrSetFreeHdlrData(nlhdlr, nlhdlrFreehdlrDataSignomial);
   SCIPnlhdlrSetFreeExprData(nlhdlr, nlhdlrFreeExprDataSignomial);
   SCIPnlhdlrSetSepa(nlhdlr, NULL, NULL, nlhdlrEstimateSignomial, NULL);

   /* parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "nlhdlr/" NLHDLR_NAME "/maxnundervars",
         "maximum number of variables when underestimating a concave power function",
         &nlhdlrdata->maxnundervars, TRUE, NLHDLR_MAXNUNDERVARS, 2, 14, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "nlhdlr/" NLHDLR_NAME "/mincutscale",
         "minimum scale factor when scaling a cut",
         &nlhdlrdata->mincutscale, TRUE, NLHDLR_MINCUTSCALE, 1e-6, 1e6, NULL, NULL) );

   return SCIP_OKAY;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
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


#define SCIP_SIGCUT_DEBUG_
#define SCIP_SIG_DEBUG

#ifdef SCIP_SIGCUT_DEBUG
#ifndef SCIP_SIG_DEBUG_
#define SCIP_SIG_DEBUG
#endif
#endif

#ifdef SCIP_SIG_DEBUG
#define SCIP_DEBUG
#endif


#include <string.h>

#include "nlhdlr_signomial.h"
#include "scip/cons_nonlinear.h"
#include "scip/struct_misc.h"
#include "scip/pub_misc_rowprep.h"
#include "scip/pub_nlhdlr.h"
#include "scip/scip_expr.h"
#include "scip/expr_var.h"


/* fundamental nonlinear handler properties */
#define NLHDLR_NAME                    "signomial"
#define NLHDLR_DESC                    "signomial handler for expressions"
#define NLHDLR_DETECTPRIORITY          30     
#define NLHDLR_ENFOPRIORITY            30

/* handler specific parameters */
#define NLHDLR_MAXNUNDERVARS           14    /**< maximum number of variables when underestimating a concave power function (maximum: 14) */
#define NLHDLR_MINCUTSCALE             1e-5  /**< minimum scale factor when scaling a cut (minimum: 1e-6) */


/**
 * nonlinear handler expression data
 * 
 * A signomial expression admits the form \f$ cx^e = y \f$, where \f$ y \f$ is an auxiliary variable representing this
 * expression. The natural fomrulation of the expression is defined as \f$ x^e = t = y/c \f$, where \f$ t \f$ is a
 * non-existant slack variable denoting \f$ y/c \f$. The variables in $x$ with positive exponents form positive
 * variables \f$ u \f$, and the associated exponents form positive exponents \f$ f\ f$. The variables in \f$ x \f$ with
 * negative exponents and \f$ t \f$  form negative variables \f$ v \f$, and the associated exponents form negative
 * exponents \f$ g \f$. Let \$ s =  \max(|f|,|g|) \$ be a normalization constant, where \$ || \$ denotes the L1 norm. Apply a scaling step:
 * Dividing the entries of \f$ f \f$  by \f$ s \f$, and dividing the entries of \f$ g \f$  by \f$ s \f$ as well. Then \f$ x^e = t \f$ has a reformulation 
 * \f$ u^f = v^g \f$, where \f$ u^f, v^g \$ are two concave power functions.
 */
struct SCIP_NlhdlrExprData
{
   SCIP_EXPR*            expr;               /**< expression */
   SCIP_Real             coef;               /**< coeffcient \f$c\f$ */
   SCIP_EXPR**           factors;            /**< expression factors representing \f$x\f$ */
   int                   nfactors;           /**< number of factors */
   int                   nvars;              /**< number of variables \f$(x,y)\f$ */
   SCIP_Real*            exponents;          /**< exponents \f$e\f$ */
   int                   nposvars;           /**< number of positive variables \f$u\f$ */
   int                   nnegvars;           /**< number of negative variables  \f$v\f$ */
   SCIP_Bool*            signs;              /**< indicators for sign of variables after reformulation, TRUE for positive, FALSE for negative */
   SCIP_Real*            refexponents;       /**< exponents of \f$(x,t)\f$ after reformulation */
   SCIP_Bool             isgetvars;          /**< are all variables already got? */

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
   int                   nexprs;             /**< the number of registered signomial expressions */

   /* parameters */
   int                   maxnundervars;      /**< maximum numbver of variables in underestimating a concave power function */
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
   SCIP *scip,
   SCIP_NLHDLREXPRDATA *nlhdlrexprdata
   )
{
   assert(nlhdlrexprdata->expr != NULL);

   /* print overall statistics and the expression */
   SCIPdebugMsg(scip, " #all variables: %d, #positive exponent variables: %d, #negative exponent variables: %d, auxvar: %s \n expr: ",
      nlhdlrexprdata->nvars, nlhdlrexprdata->nposvars, nlhdlrexprdata->nnegvars, SCIPvarGetName(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->expr)) );
   SCIPprintExpr(scip, nlhdlrexprdata->expr, NULL);

   /* if variables are detected, we can print more information */
   if( !nlhdlrexprdata->isgetvars )
   {
      SCIPinfoMessage(scip, NULL, "\n");
      return SCIP_OKAY;
   }

   /* print the natural formulation of the expression */   
   SCIPinfoMessage(scip, NULL, ". natural formulation c x^a = y: %.2f", nlhdlrexprdata->coef);
   for( int i = 0; i < nlhdlrexprdata->nvars - 1; i++ )
      SCIPinfoMessage(scip, NULL, "%s^%.2f", SCIPvarGetName(nlhdlrexprdata->vars[i]), nlhdlrexprdata->exponents[i]);
   SCIPinfoMessage(scip, NULL, " = %s", SCIPvarGetName(nlhdlrexprdata->vars[nlhdlrexprdata->nvars - 1]));

   /* print the  reformulation of the expression */  
   SCIPinfoMessage(scip, NULL, ". Reformulation u^f = v^g: ");
   if( nlhdlrexprdata->nposvars == 0 )
   {
      SCIPinfoMessage(scip, NULL, "1");
   }
   else
   {
      for( int i = 0; i < nlhdlrexprdata->nvars; i++ )
      {
         if( nlhdlrexprdata->signs[i] )
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
               SCIPinfoMessage(scip, NULL, "(%s/%.2f)^%.2f", SCIPvarGetName(nlhdlrexprdata->vars[i]), nlhdlrexprdata->coef, nlhdlrexprdata->refexponents[i]);
            else
               SCIPinfoMessage(scip, NULL, "%s^%.2f", SCIPvarGetName(nlhdlrexprdata->vars[i]), nlhdlrexprdata->refexponents[i]);
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
      SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->signs, (*nlhdlrexprdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->refexponents, (*nlhdlrexprdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->vars, (*nlhdlrexprdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->intervals, (*nlhdlrexprdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->box, 2 * (*nlhdlrexprdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->xstar, (*nlhdlrexprdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->facetcoefs, (*nlhdlrexprdata)->nvars);
   }
   SCIPfreeBlockMemory(scip, nlhdlrexprdata);
   *nlhdlrexprdata = NULL;
}

/** reform rowprep to a standard form for nonlinear handlers.
 * A rowprep in standard form only contains an estimator of the expression and no auxvar.
 */
static
SCIP_Real reformRowprep(
   SCIP*                 scip,
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,
   SCIP_ROWPREP*         rowprep,
   SCIP_Bool*            success
)
{
   assert(rowprep != NULL);
   assert(nlhdlrexprdata != NULL);

   SCIP_Real* coefs =  SCIProwprepGetCoefs(rowprep);
   SCIP_VAR** vars = SCIProwprepGetVars(rowprep);
   int nvars = SCIProwprepGetNVars(rowprep);

   /* find auxvar's cut coefficient and set it to zero */
   SCIP_VAR* auxvar = nlhdlrexprdata->vars[nlhdlrexprdata->nfactors];
   SCIP_Real coefauxvar = 0;
   for( int i = 0; i < nvars; i++ )
   {  
      if( vars[i] != auxvar )
         continue;
      coefauxvar = coefs[i];
      coefs[i] = 0;
   }

   if( SCIPisZero(scip, coefauxvar) )
   {
      *success = FALSE;
      return 0.0;
   }

   /* the reformation scales the cut, such that coefficients and constant are dividing by the absolute value of coefauxvar */
   coefauxvar = 1 / fabs(coefauxvar);
   
   for( int i = 0; i < nvars; i++ )
   {  
      if( vars[i] == auxvar )
         continue;
      coefs[i] *= coefauxvar;
   }

   rowprep->side *= coefauxvar;

   return coefauxvar;
}

/** get variables associated with the expression and its subexpressions */
static
SCIP_RETCODE getVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata      /**< expression data */
   )
{
   int c;

   assert(!nlhdlrexprdata->isgetvars);

   /* get and capture variables \f$x\f$ */
   for( c = 0; c < nlhdlrexprdata->nfactors; ++c )
   {
      nlhdlrexprdata->vars[c] = NULL;
      nlhdlrexprdata->vars[c] = SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->factors[c]);
      assert(nlhdlrexprdata->vars[c] != NULL);
      SCIP_CALL( SCIPcaptureVar(scip, nlhdlrexprdata->vars[c]) );
   }

   /* get and capture variable \f$y\f$ */
   nlhdlrexprdata->vars[c] = NULL;
   nlhdlrexprdata->vars[c] = SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->expr);
   assert(nlhdlrexprdata->vars[c] != NULL);
   SCIP_CALL( SCIPcaptureVar(scip, nlhdlrexprdata->vars[c]) );

   nlhdlrexprdata->isgetvars = TRUE;
   return SCIP_OKAY;
}

/** get bounds of variables $x,t$ and check whether they are box constrained signomial variables */
static
void getCheckBds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< expression data */
   SCIP_Bool*            isboxsignomial      /**< buffer to store whether variables are box constrained signomial variables */
   )
{
   SCIP_Real productinf = 1;
   SCIP_Real productsup = 1;
   int c;

   assert(nlhdlrexprdata->isgetvars);

   *isboxsignomial = FALSE;

   /* get bounds of \f$x\f$ */
   for( c = 0; c < nlhdlrexprdata->nfactors; c++ )
   {
      SCIP_Real inf = SCIPvarGetLbLocal(nlhdlrexprdata->vars[c]);
      SCIP_Real sup = SCIPvarGetUbLocal(nlhdlrexprdata->vars[c]);

      /* if the bounds of the variable are not positive and finite, then the expression is not a signomial */
      if( !SCIPisPositive(scip, inf) || !SCIPisPositive(scip, sup) || SCIPisInfinity(scip, sup) )
         return;
      nlhdlrexprdata->intervals[c].inf = inf;
      nlhdlrexprdata->intervals[c].sup = sup;
      SCIP_Real powinf = pow(inf, nlhdlrexprdata->exponents[c]);
      SCIP_Real powsup = pow(sup, nlhdlrexprdata->exponents[c]);
      productinf *= fmin(powinf, powsup);
      productsup *= fmax(powinf, powsup);
   }

   /* compute bounds of \f$t\f$  by bounds of \f$x\f$ */
   nlhdlrexprdata->intervals[c].inf = productinf;
   nlhdlrexprdata->intervals[c].sup = productsup;

   #ifdef SCIP_SIG_DEBUG_
      if( !SCIPisEQ(scip, productinf, SCIPvarGetLbLocal(nlhdlrexprdata->vars[c])) ){
         SCIPdebugMsg(scip, "%f %f \n", productinf, SCIPvarGetLbLocal(nlhdlrexprdata->vars[c]));
      }
      if( !SCIPisEQ(scip, productsup, SCIPvarGetUbLocal(nlhdlrexprdata->vars[c])) ){
         SCIPdebugMsg(scip, "%f %f \n", productsup, SCIPvarGetUbLocal(nlhdlrexprdata->vars[c]));
      }
   #endif

   *isboxsignomial = TRUE;
}

/** evaluate expression at solution w.r.t. auxiliary variables */
static 
SCIP_DECL_VERTEXPOLYFUN(nlhdlrExprEvalPower)
{
   VERTEXPOLYFUN_EVALDATA *evaldata = (VERTEXPOLYFUN_EVALDATA *)funcdata;

   assert(args != NULL);
   assert(nargs == evaldata->nsignvars);
   assert(evaldata != NULL);

   SCIP_NLHDLREXPRDATA *nlhdlrexprdata = evaldata->nlhdlrexprdata;

#ifdef SCIP_MORE_DEBUG
   SCIPdebugMsg(evaldata->scip, "eval vertexpolyfun at\n");
#endif
   SCIP_Real val = 1;
   int j = 0;
   for ( int i = 0 ; i < nlhdlrexprdata->nvars; ++i )
   {
      if( nlhdlrexprdata->signs[i] != evaldata->sign ) 
         continue;
      /* the reformulated exponent of argrs[j] is found */
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
   SCIP*                 scip,                               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA   *nlhdlrexprdata,                    /**< nonlinear handler expression data */
   SCIP_Bool             sign,                               /**< the sign of variables of the power function */     
   SCIP_Real             multiplier,                         /**< the mulitplier of the estimator */ 
   SCIP_Bool             overestimate,                       /**< whether overestimate or underestimator the power function */ 
   SCIP_SOL*             sol,                                /**< solution \f$ w' \f$  to use in estimation */
   SCIP_ROWPREP*         rowprep,                            /**< rowprep where to store estimator */
   SCIP_Bool*            isspecial,                          /**< buffer to store whether this function is a special case */ 
   SCIP_Bool*            success                             /**< buffer to store whether successful */
   )
{
   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(rowprep != NULL);
   assert(isspecial != NULL);
   assert(success != NULL);

   *success = FALSE;
   SCIP_Bool allfixed = TRUE;

   /* check whether all variables are fixed */
   for( int i = 0; i < nlhdlrexprdata->nvars; ++i )
   {
      if( nlhdlrexprdata->signs[i] != sign )
         continue;

      if( !SCIPisRelEQ(scip, nlhdlrexprdata->intervals[i].inf, nlhdlrexprdata->intervals[i].sup) )
         allfixed = FALSE;
   }


   /* allfixed is a special case */
   if( allfixed )
   {
      /* SCIPcomputeFacetVertexPolyhedralNonlinear prints a warning and does not succeed if all is fixed */
      SCIPdebugMsg(scip, "all variables fixed, skip estimate\n");
      *isspecial = TRUE;
      return SCIP_OKAY;
   }

   /* number of variables in the power function */
   int nsignvars = sign ? nlhdlrexprdata->nposvars : nlhdlrexprdata->nnegvars;

   /* if the power function has no more than 2 variables, this a special case */
   *isspecial = ( nsignvars <= 1 ) || ( nsignvars == 2  && !overestimate );
   if( !*isspecial )
      return SCIP_OKAY;
   
   if( nsignvars == 0 )
   {
      /* constant case */
      SCIProwprepAddConstant(rowprep, multiplier * 1);
      *success = TRUE;
   }
   else if( nsignvars == 1 )
   {
      /* univariate case, \f$ w^h \f$ */
      for( int i = 0; i < nlhdlrexprdata->nvars; ++i )
      {
         if( nlhdlrexprdata->signs[i] == sign )
         {
            /* the variable \f$ w \f$ is found */
            SCIP_Real scale = i == (nlhdlrexprdata->nvars - 1) ? nlhdlrexprdata->coef : 1;
            SCIP_VAR* var = nlhdlrexprdata->vars[i];
            SCIP_Real refexponent = nlhdlrexprdata->refexponents[i];
            if( SCIPisEQ(scip, refexponent, 1.0) )
            {
               /* \f$ h = 1 \f$, a univariate linear function. Only rescale, no need for linearization */
               SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, multiplier / scale));
            }
            else
            {
               /* a univariate power function */
               SCIP_Real facetconstant;
               SCIP_Real facetcoef;
               SCIP_Real val = SCIPgetSolVal(scip, sol, var) / scale;
               if( overestimate )
               {
                  /* overestimate by gradient cut: \f$ w'^h + h  w'^(h-1) (w - w') =  (1- h) w'^h + h  w'^(h-1) w  \f$ */
                  SCIP_Real funcval = pow(val, refexponent);
                  facetconstant = (1 - refexponent) * funcval;
                  facetcoef = refexponent * funcval / val;
               }
               else
               {
                  /* underestimate by secant approximation: \f$ ((sup^h - inf^h) / (sup - inf)) (w - inf) + inf^h \f$ */
                  SCIP_Real inf = nlhdlrexprdata->intervals[i].inf;
                  SCIP_Real sup = nlhdlrexprdata->intervals[i].sup;
                  SCIP_Real powinf = pow(inf, refexponent);
                  facetcoef = (pow(sup, refexponent) - powinf) / (sup - inf);
                  facetconstant = powinf - facetcoef * inf;
               }
               SCIProwprepAddConstant(rowprep,  multiplier * facetconstant);
               SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, multiplier * facetcoef / scale));
            }
         }
      }
      *success = TRUE;
   }
   else if( nsignvars == 2 && !overestimate ){
      /* bivariate case, \f$ f(w) = w^h = f_0(w_0) f_1(w_1)  \f$ */
      SCIP_VAR* vars[2] = {NULL, NULL}; 
      SCIP_Real refexponents[2] = {0.0, 0.0};
      SCIP_Real xstar[2] = {0.0, 0.0};;
      SCIP_Real scale[2] = {0.0, 0.0};; 
      SCIP_Real box[4] = {0.0, 0.0, 0.0, 0.0};
      int j = 0;
      /* get refexponents, xtar, scale, and box */
      for( int i = 0; i < nlhdlrexprdata->nvars; ++i )
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
      SCIP_Real dw0 = box[1] - box[0];
      SCIP_Real dw1 = box[3] - box[2];
      /* determine the location (upper right or lower left half sapce) of the xtar. 
       * \f$ (dw1, dw0) \f$ is the direction vector to the upper right half space. */
      SCIP_Bool isupperright = ( (xstar[0] - box[0]) * dw1 + (xstar[1] - box[3]) * dw0 ) > 0;
      /* compute function values of \f$ f_0, f_1 \f$ at vetices */
      SCIP_Real f0w0l = pow(box[0], refexponents[0]);
      SCIP_Real f0w0u = pow(box[1], refexponents[0]);
      SCIP_Real f1w1l = pow(box[2], refexponents[1]);
      SCIP_Real f1w1u = pow(box[3], refexponents[1]);
      SCIP_Real fw0lw1u = f0w0l * f1w1u;
      SCIP_Real fw0uw1l = f0w0u * f1w1l;
      SCIP_Real facetcoefs[2];
      SCIP_Real facetconstant;
      if( isupperright )
      {
         /* underestimator: \f$ fw0uw1u + (fw0uw1u - fw0lw1u) / (dw0) * (x0 - w0u) + (fw0uw1u - fw0uw1l) / (dw1) * (x1 - w1u) \f$ */
         SCIP_Real fw0uw1u = f0w0u * f1w1u;
         facetcoefs[0] = (fw0uw1u - fw0lw1u) / dw0;
         facetcoefs[1] = (fw0uw1u - fw0uw1l) / dw1;
         facetconstant = fw0uw1u  - facetcoefs[0] * box[1] - facetcoefs[1] * box[3];
      }
      else
      {
         /* underestimator: \f$ fw0lw1l + (fw0uw1l - fw0lw1l) / (dw0) * (x0 - w0l) + (fw0lw1u- fw0lw1l) / (dw1) * (x1 - w1l) \f$ */
         SCIP_Real fw0lw1l = f0w0l * f1w1l;
         facetcoefs[0] = (fw0uw1l - fw0lw1l) / dw0;
         facetcoefs[1] = (fw0lw1u - fw0lw1l) / dw1;
         facetconstant = fw0lw1l  - facetcoefs[0] * box[0] - facetcoefs[1] * box[2];
      }
      SCIProwprepAddConstant(rowprep,  multiplier * facetconstant);
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, vars[0], multiplier * facetcoefs[0] / scale[0]));
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, vars[1], multiplier * facetcoefs[1] / scale[1]));
      *success = TRUE;
   }

   return SCIP_OKAY;
}

/** adds an under estimator for a concave power concave function \f$ w^h \f$ to a given rowprep
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
   SCIP_Real             multiplier,         /**< the mulitplier of the estimator */ 
   SCIP_SOL*             sol,                /**< solution \f$ w' \f$ to use*/
   SCIP_Real             targetvalue,        /**< a target value to achieve; if not reachable, then can give up early */
   SCIP_ROWPREP*         rowprep,            /**< rowprep where to store estimator */
   SCIP_Bool*            success             /**< buffer to store whether successful */
   )
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(rowprep != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* number of variables of sign */
   int nsignvars = sign ? nlhdlrexprdata->nposvars : nlhdlrexprdata->nnegvars;

   /* data struture to evaluate the power function */
   VERTEXPOLYFUN_EVALDATA evaldata;
   evaldata.nlhdlrexprdata = nlhdlrexprdata;
   evaldata.sign = sign;
   evaldata.nsignvars = nsignvars;
   evaldata.scip = scip;

   SCIP_Real* xstar = nlhdlrexprdata->xstar;
   SCIP_Real* box = nlhdlrexprdata->box;
   /* compute box constraints and reference value of variables*/
   int j = 0;
   for( int i = 0; i < nlhdlrexprdata->nvars; ++i )
   {
      if( nlhdlrexprdata->signs[i] != sign )
         continue;

      box[2 * j] = nlhdlrexprdata->intervals[i].inf;
      box[2 * j + 1] = nlhdlrexprdata->intervals[i].sup;
      SCIP_Real scale = i == (nlhdlrexprdata->nvars - 1) ? nlhdlrexprdata->coef : 1;
      xstar[j] = SCIPgetSolVal(scip, sol, nlhdlrexprdata->vars[i]) / scale;
      j++;
   }

   SCIP_Real* facetcoefs = nlhdlrexprdata->facetcoefs;
   SCIP_Real facetconstant;
   /* find a facet of the underestimator */
   SCIP_CALL( SCIPcomputeFacetVertexPolyhedralNonlinear(scip, conshdlr, FALSE, nlhdlrExprEvalPower, (void*)&evaldata, xstar, box,
      nsignvars, targetvalue, success, facetcoefs, &facetconstant));

   if( !*success )
   {
      SCIPdebugMsg(scip, "failed to compute facet of convex hull\n");
      return SCIP_OKAY;
   }

   /* set coefficients in the rowprep */
   SCIProwprepAddConstant(rowprep, multiplier * facetconstant);
   j = 0;
   for( int i = 0; i < nlhdlrexprdata->nvars; ++i )
   {
      if( nlhdlrexprdata->signs[i] != sign )
         continue;
      SCIP_Real scale = i == (nlhdlrexprdata->nvars - 1) ? nlhdlrexprdata->coef : 1;
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, nlhdlrexprdata->vars[i], multiplier * facetcoefs[j] / scale));
      j++;
   }

   return SCIP_OKAY;
}

/** adds an over estimator for a concave power concave function \f$ w^h \f$ to a given rowprep
 *
 * The estimator is multiplied by the multiplier and stored in the rowprep.
 */
static SCIP_RETCODE overEstimatePower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nonlinear handler expression data */
   SCIP_Bool             sign,               /**< the sign of variables of the power function */ 
   SCIP_Real             multiplier,         /**< the mulitplier of the estimator */ 
   SCIP_SOL*             sol,                /**< solution to use*/
   SCIP_ROWPREP*         rowprep,            /**< rowprep where to store estimator */
   SCIP_Bool*            success             /**< buffer to store whether successful */
   )
{
   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(rowprep != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* compute the value and the sum of reformulated exponents of the power function */
   SCIP_Real sumrefexponents = 0;
   SCIP_Real funcval = 1;
   for( int i = 0; i < nlhdlrexprdata->nvars; ++i )
   {
      if( nlhdlrexprdata->signs[i] != sign )
         continue;
      SCIP_Real scale = i == (nlhdlrexprdata->nvars - 1) ? nlhdlrexprdata->coef : 1;
      SCIP_Real val = SCIPgetSolVal(scip, sol, nlhdlrexprdata->vars[i]) / scale;
      SCIP_Real refexponent = nlhdlrexprdata->refexponents[i];
      sumrefexponents += refexponent;
      funcval *= pow(val, refexponent);
   }

   /* overestimate by gradient cut:  \f$ w'^h + h w'^(h - 1)(w - w')  \f$ */
   SCIP_Real facetconstant = (1 - sumrefexponents) * funcval;
   SCIProwprepAddConstant(rowprep,  multiplier * facetconstant);
   for( int i = 0; i < nlhdlrexprdata->nvars; ++i )
   {
      if( nlhdlrexprdata->signs[i] != sign )
         continue;
      SCIP_Real scale = i == (nlhdlrexprdata->nvars - 1) ? nlhdlrexprdata->coef : 1;
      SCIP_VAR* var = nlhdlrexprdata->vars[i];
      SCIP_Real val = SCIPgetSolVal(scip, sol, var) / scale;
      SCIP_Real facetcoef = nlhdlrexprdata->refexponents[i] * funcval / val;
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, multiplier * facetcoef / scale));
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

   assert(conshdlr != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(rowpreps != NULL);
   assert(expr == nlhdlrexprdata->expr);

   *success = FALSE;
   *addedbranchscores = FALSE;
	
   /* get variables*/
   if( !nlhdlrexprdata->isgetvars )
      SCIP_CALL(getVars(scip, nlhdlrexprdata));

   /* check whether the expression is a signomial term */
   SCIP_Bool isboxsignomial = FALSE;
   getCheckBds(scip, nlhdlrexprdata, &isboxsignomial);

   if( !isboxsignomial )
      return SCIP_OKAY;

   SCIP_NLHDLRDATA *nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   SCIP_Bool undersign;
   SCIP_Bool oversign;
   SCIP_Real undermultiplier;
   SCIP_Real overmultiplier;
   int nundervars;

   /* determine the sign of variables for over and under estimators, and the multiplier for estimators in the rowprep */
   if( overestimate )
   {
      /* relax \f$ t \le x^a \f$, i.e., \f$ 0 \le u^f - v^g \f$, overestimate \f$ u^f \$f, underestimate \f$ v^g \f$ */
      oversign = TRUE;
      undersign = FALSE;
      overmultiplier = 1;
      undermultiplier = -1;
      nundervars = nlhdlrexprdata->nnegvars;
   }
   else
   {
      /* relax \f$ x^a \le t \f$, i.e.,  \f$ u^f - v^g \le 0 \f$, underestimate \f$ u^f \f$,  overestimate \f$ v^g \f$ */
      oversign = FALSE;
      undersign = TRUE;
      overmultiplier = -1;
      undermultiplier = 1;
      nundervars = nlhdlrexprdata->nposvars;
   }

   /* compute the value of the overestimator, which is a target value for computing the underestimator efficiently */
   SCIP_Real targetunder = 1;
   for( int i = 0; i < nlhdlrexprdata->nvars; i++ )
   {
      if( nlhdlrexprdata->signs[i] == oversign )
      {
         SCIP_Real scale = i == (nlhdlrexprdata->nvars - 1) ? nlhdlrexprdata->coef : 1;
         targetunder *= pow(SCIPgetSolVal(scip, sol, nlhdlrexprdata->vars[i]) / scale, nlhdlrexprdata->refexponents[i]);
      }
   }
#ifdef SCIP_SIGCUT_DEBUG
   /* print information on estimators */
   SCIP_CALL( printSignomial(scip, nlhdlrexprdata));
   SCIPinfoMessage(scip, NULL, " Auxvalue: %f, targetvalue: %f, %sestimate.", auxvalue, targetvalue, overestimate ? "over" : "under");
   SCIP_Real targetover = 1;
   for ( int i = 0; i < nlhdlrexprdata->nvars; i++ )
   {
      if( nlhdlrexprdata->signs[i] == undersign )
      {
         SCIP_Real scale = i == (nlhdlrexprdata->nvars - 1) ? nlhdlrexprdata->coef : 1;
         targetover *= pow(SCIPgetSolVal(scip, sol, nlhdlrexprdata->vars[i]) / scale, nlhdlrexprdata->refexponents[i]);
      }
   }
   SCIPinfoMessage(scip, NULL, " Underestimate: %s, overestimate: %s.", 
      undersign ? "positive" : "negative", oversign ? "positive" : "negative");
   SCIPinfoMessage(scip, NULL, " Undervalue (targetover): %f, overvalue (targetunder): %f.", targetover, targetunder);
#endif

   /* create a rowprep and allocate memory for it */
   SCIP_ROWPREP *rowprep;
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, TRUE));
   SCIP_CALL( SCIPensureRowprepSize(scip, rowprep, nlhdlrexprdata->nvars + 1));

   SCIP_Bool isspecial;
   /* only underesimate a concave function, if the number of variables satisfy the criteria */
   if( nundervars <= nlhdlrdata->maxnundervars )
   {
      /* compute underestimator */
      isspecial = FALSE;
      /* check and compute the special case */
      SCIP_CALL( estimateSpecialPower(scip, nlhdlrexprdata, undersign, undermultiplier, FALSE, sol, rowprep, &isspecial, success));
      if( !isspecial )
         SCIP_CALL( underEstimatePower(scip, conshdlr, nlhdlr, nlhdlrexprdata, undersign, undermultiplier, sol, targetunder, rowprep, success));
   }


   if( *success )
   {
      /* compute overestimator */
      isspecial = FALSE;
      #ifdef SCIP_SIGCUT_DEBUG
      SCIP_Real prevconstant =  SCIProwprepGetSide(rowprep);
      #endif
      /* check and compute the special case */
      SCIP_CALL( estimateSpecialPower(scip, nlhdlrexprdata, oversign, overmultiplier, TRUE, sol, rowprep, &isspecial, success));
      if( !isspecial )
         SCIP_CALL( overEstimatePower(scip, nlhdlrexprdata, oversign, overmultiplier, sol, rowprep, success));

   }


   if( *success )
   {
      SCIP_Real scale = reformRowprep(scip, nlhdlrexprdata, rowprep, success);
      *success = *success && SCIPisGT(scip, scale, nlhdlrdata->mincutscale);
      if( *success )
         SCIP_CALL( SCIPsetPtrarrayVal(scip, rowpreps, 0, rowprep));
   }
   if(!*success)
      SCIPfreeRowprep(scip, &rowprep);

   return SCIP_OKAY;
}


/** callback to detect structure in expression tree */
static 
SCIP_DECL_NLHDLRDETECT(nlhdlrDetectSignomial)
{ /*lint --e{715}*/
   SCIP_NLHDLRDATA *nlhdlrdata;
   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);

   assert(nlhdlrdata != NULL);
   assert(expr != NULL);
   assert(enforcing != NULL);
   assert(participating != NULL);

   /* for now, we do not care about separation if it is not required */
   if( (*enforcing & SCIP_NLHDLR_METHOD_SEPABOTH) == SCIP_NLHDLR_METHOD_SEPABOTH )
   {
      return SCIP_OKAY;
   }

   /* check for product expressions with more than one children */
   if( SCIPisExprProduct(scip, expr) && SCIPexprGetNChildren(expr) >= 2 )
   {
      int nf = SCIPexprGetNChildren(expr);
      int nvars = nf + 1;

      /* create expression data for the nonlinear handler */
      SCIP_CALL( SCIPallocClearBlockMemory(scip, nlhdlrexprdata));
      (*nlhdlrexprdata)->nfactors = nf;
      (*nlhdlrexprdata)->nvars = nvars;

      SCIP_EXPR *** ptrfactors = &(*nlhdlrexprdata)->factors;
      SCIP_Real ** ptrexponents = &(*nlhdlrexprdata)->exponents;
      /* allocat memory for expression data */
      SCIPallocBlockMemoryArray(scip, ptrfactors, nf);
      SCIPallocBlockMemoryArray(scip, ptrexponents, nf);

      /* get monomial information */
      SCIP_CALL( SCIPgetExprMonomialData(scip, expr, &((*nlhdlrexprdata)->coef), *ptrexponents, *ptrfactors) );

      /* skip multilinear terms */
      SCIP_Bool ismultilinear = FALSE;
      if( nf >= 2)
      {
         ismultilinear = TRUE;
         for ( int c = 0; c < nf; c++ )
         {
            if( !SCIPisEQ(scip, (*nlhdlrexprdata)->exponents[c], 1.0) )
            {
               ismultilinear = FALSE;
               break;
            }
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
         /* allocate more memory for expression data */
         SCIP_Bool ** ptrsigns =  &(*nlhdlrexprdata)->signs;
         SCIP_Real ** ptrrefexponents  = &(*nlhdlrexprdata)->refexponents;
         SCIP_VAR *** ptrvars = &(*nlhdlrexprdata)->vars;
         SCIP_INTERVAL ** ptrintervals = &(*nlhdlrexprdata)->intervals;
         SCIP_Real ** ptrxstar =  &(*nlhdlrexprdata)->xstar;
         SCIP_Real ** ptrfacetcoefs = &(*nlhdlrexprdata)->facetcoefs;
         SCIP_Real ** ptrbox = &(*nlhdlrexprdata)->box;
         SCIPallocBlockMemoryArray(scip, ptrsigns, nvars);
         SCIPallocBlockMemoryArray(scip, ptrrefexponents, nvars);
         SCIPallocBlockMemoryArray(scip, ptrvars, nvars);
         SCIPallocBlockMemoryArray(scip, ptrintervals, nvars);
         SCIPallocBlockMemoryArray(scip, ptrxstar, nvars);
         SCIPallocBlockMemoryArray(scip, ptrfacetcoefs, nvars);
         SCIPallocBlockMemoryArray(scip, ptrbox, 2 * nvars);

         (*nlhdlrexprdata)->isgetvars = FALSE;
         (*nlhdlrexprdata)->expr = expr;
         SCIPcaptureExpr(expr);

         /* detect more information in reformulation,  and we first 
          * compute the sum of positive and negative exponents and determine the sign indicators 
          */
         SCIP_Real sumlexponents = 0;
         SCIP_Real sumrexponents = 1;
         int nposvars = 0;
         for (int c = 0; c < nf; c++)
         {
            /* capture sub expressions */
            SCIPcaptureExpr((*nlhdlrexprdata)->factors[c]);
            if( SCIPisPositive(scip, (*nlhdlrexprdata)->exponents[c]) )
            {
               /* add a positive exponent */
               sumlexponents += (*nlhdlrexprdata)->exponents[c];
               (*nlhdlrexprdata)->signs[c] = TRUE;
               nposvars++;
            }
            else
            {
               /* minus a negative exponent */
               sumrexponents -= (*nlhdlrexprdata)->exponents[c];
               (*nlhdlrexprdata)->signs[c] = FALSE;
            }
         }
         (*nlhdlrexprdata)->signs[nf] = FALSE;
         (*nlhdlrexprdata)->nposvars = nposvars;
         (*nlhdlrexprdata)->nnegvars = nf - nposvars + 1;

         /* compute the normalization constant */
         SCIP_Real normalize = fmax(sumlexponents, sumrexponents);
         /* normalize positive and negative exponents */
         for (int c = 0; c < nf; c++)
         {
            if( (*nlhdlrexprdata)->signs[c] )
               (*nlhdlrexprdata)->refexponents[c] = (*nlhdlrexprdata)->exponents[c] / normalize;
            else
               (*nlhdlrexprdata)->refexponents[c] = -(*nlhdlrexprdata)->exponents[c] / normalize;
         }
         (*nlhdlrexprdata)->refexponents[nf] = 1 / normalize;

         /* tell children that we will use their auxvar and use its activity for both estimate and domain propagation */
         for( int c = 0; c < nf; c++ )
            SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, (*nlhdlrexprdata)->factors[c], TRUE, FALSE, TRUE, TRUE));

         nlhdlrdata->nexprs++;
      }
   }

   if( *nlhdlrexprdata != NULL )
   {
      /* we want to join separation, if not disabled by parameter */
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

   *auxvalue = nlhdlrexprdata->coef;
   for ( int c = 0; c < nlhdlrexprdata->nfactors; ++c )
   {
      assert(nlhdlrexprdata->factors[c] != NULL);
      SCIP_VAR* var = SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->factors[c]);
      assert(var != NULL);
      SCIP_Real val = SCIPgetSolVal(scip, sol, var);
      *auxvalue *= pow(val, nlhdlrexprdata->exponents[c]);
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

   SCIP_CALL (SCIPincludeNlhdlrSignomial(targetscip));

   return SCIP_OKAY;
}

/** callback to free data of handler */
static 
SCIP_DECL_NLHDLRFREEHDLRDATA(nlhdlrFreehdlrDataSignomial)
{ /*lint --e{715}*/
   assert(nlhdlrdata != NULL);
   assert((*nlhdlrdata)->nexprs == 0);

   SCIPfreeBlockMemory(scip, nlhdlrdata);

   return SCIP_OKAY;
}

/** callback to free expression specific data */
static
SCIP_DECL_NLHDLRFREEEXPRDATA(nlhdlrFreeExprDataSignomial)
{ /*lint --e{715}*/
   SCIP_NLHDLRDATA *nlhdlrdata;
   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);

   assert(expr != NULL);

   /* release expressions */
   SCIP_CALL(SCIPreleaseExpr(scip, &(*nlhdlrexprdata)->expr));
   for( int c = 0; c < (*nlhdlrexprdata)->nfactors; c++ )
      SCIP_CALL(SCIPreleaseExpr(scip, &(*nlhdlrexprdata)->factors[c]));

   /* release variables */
   if( (*nlhdlrexprdata)->isgetvars )
   {
      for( int c = 0; c < (*nlhdlrexprdata)->nvars; c++ ) 
         SCIP_CALL(SCIPreleaseVar(scip, &(*nlhdlrexprdata)->vars[c]));
   }

   /* free memory */
   freeExprDataMem(scip, nlhdlrexprdata, FALSE);

   nlhdlrdata->nexprs--;

   return SCIP_OKAY;
}


/*
 * nonlinear handler specific interface methods
 */

/** includes signomial nonlinear handler to consexpr */
SCIP_RETCODE 
SCIPincludeNlhdlrSignomial(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NLHDLRDATA * nlhdlrdata;
   SCIP_NLHDLR *nlhdlr;

   assert(scip != NULL);

   /* create nonlinear handler specific data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &nlhdlrdata));
   BMSclearMemory(nlhdlrdata);

   SCIP_CALL( SCIPincludeNlhdlrNonlinear(scip, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_DETECTPRIORITY,
      NLHDLR_ENFOPRIORITY, nlhdlrDetectSignomial, nlhdlrEvalauxSignomial, nlhdlrdata));
   assert(nlhdlr != NULL);

   SCIPnlhdlrSetCopyHdlr(nlhdlr, nlhdlrCopyhdlrSignomial);
   SCIPnlhdlrSetFreeHdlrData(nlhdlr, nlhdlrFreehdlrDataSignomial);
   SCIPnlhdlrSetFreeExprData(nlhdlr, nlhdlrFreeExprDataSignomial);
   SCIPnlhdlrSetSepa(nlhdlr, NULL, NULL, nlhdlrEstimateSignomial, NULL);

   /* parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "nlhdlr/" NLHDLR_NAME "/maxnundervars",
         "maximum number of variables when underestimating a concave power function", 
         &nlhdlrdata->maxnundervars, TRUE, NLHDLR_MAXNUNDERVARS, 2, 14, NULL, NULL));
   SCIP_CALL( SCIPaddRealParam(scip, "nlhdlr/" NLHDLR_NAME "/mincutscale",
         "minimum scale factor when scaling a cut", 
         &nlhdlrdata->mincutscale, TRUE, NLHDLR_MINCUTSCALE, 1e-6, 1e6, NULL, NULL));

   return SCIP_OKAY;
}

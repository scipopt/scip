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

/**@file   reader_zpl.c
 * @brief  ZIMPL model file reader
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/reader_zpl.h"

#ifdef WITH_ZIMPL

#include <assert.h>
#include <unistd.h>
#include <string.h>

#include "scip/cons_exactlp.h"
#include "scip/cons_linear.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "scip/cons_indicator.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_nonlinear.h"
#include "scip/intervalarith.h"
#include "scip/pub_misc.h"

#ifdef WITH_EXACTSOLVE
/* @Note: gmp.h has to be included before zimpl/mme.h */
#include "gmp.h"
#endif

/* @Note: Due to dependencies we need the following order. */
/* include the ZIMPL headers necessary to define the LP and MINLP construction interface */
#include "zimpl/bool.h"
#include "zimpl/ratlptypes.h"
#include "zimpl/mme.h"

#include "zimpl/numb.h"
#include "zimpl/bound.h"
#include "zimpl/mono.h"
#include "zimpl/term.h"

#include "zimpl/xlpglue.h"
#include "zimpl/zimpllib.h"

#define READER_NAME             "zplreader"
#define READER_DESC             "file reader for ZIMPL model files"
#define READER_EXTENSION        "zpl"

#ifdef WITH_EXACTSOLVE
#define LARGEBOUND              1e+06
#endif

/*
 * LP construction interface of ZIMPL
 */

/* we only support ZIMPL with a version higher than 3.2.0 */
#if (ZIMPL_VERSION >= 320)

/* ZIMPL does not support user data in callbacks - we have to use static variables */
struct
SCIP_ReaderData
{
   SCIP*                 scip;               /**< scip data structure */
   SCIP_SOL*             sol;                /**< primal solution candidate */
   SCIP_Bool             valid;              /**< is the primal solution candidate valid */
   SCIP_Bool             branchpriowarning;  /**< store if the waring regarding fractional value for the branching
                                              *   priority was already posted */
   SCIP_Bool             readerror;          /**< was a reading be discovered */
#ifdef WITH_EXACTSOLVE
   SCIP_CONSHDLRDATA*    conshdlrdata;       /**< exactlp constraint handler data */
   SCIP_OBJSENSE         objsense;           /**< objective sense */
   SCIP_VAR**            vars;               /**< variables in the order they are added to the problem (var->index) */
   int                   nvars;              /**< number of variables */
   int                   varssize;           /**< size of variable specific arrays */
   int                   ninfbounds;         /**< number of variables with infinite bound in safe dual bounding method */
   int                   ninfintbounds;      /**< number of integer variables with infinite bound in safe db method */
   int                   nlargebounds;       /**< number of variables with large bound in safe dual bounding method */
   mpq_t*                obj;                /**< objective function values of variables */
   mpq_t*                lb;                 /**< lower bounds of variables */
   mpq_t*                ub;                 /**< upper bounds of variables */
   int                   nconss;             /**< number of constraints */
   int                   consssize;          /**< size of constraints specific arrays */
   int                   nsplitconss;        /**< number of constraints we would have to be split for a FP-relaxation */
   mpq_t*                lhs;                /**< left hand sides of constraints */
   mpq_t*                rhs;                /**< right hand sides of constraints */
   int                   nnonz;              /**< number of nonzero elements in the constraint matrix */
   int                   nonzsize;           /**< size of non-zero entries specific arrays */
   int                   nintegral;          /**< number of integral nonzero elements in the constraint matrix */
   int*                  beg;                /**< start index of each constraint in ind and val array */
   int*                  len;                /**< number of nonzeros in val array corresponding to constraint */
   int*                  ind;                /**< variable indices (var->probindex) of constraint matrix entries */
   mpq_t*                val;                /**< values of nonzero constraint matrix entries (and some zeros) */
   mpq_t                 minabsval;          /**< minimum absolute nonzero constraint matrix, lhs, or rhs entry */
   mpq_t                 maxabsval;          /**< maximum absolute nonzero constraint matrix, lhs, or rhs entry */
   SCIP_Bool             objneedscaling;     /**< do objective values need scaling because some aren't FP-representable? */
#endif
};

/** Allocate storage for the mathematical program instance generated by ZIMPL. xlp_alloc() is the first xlpglue routine
 *  that will be called by ZIMPL. The user_data pointer may hold an arbitray value.
 */
Lps* xlp_alloc(
   const char*           name,               /**< name of the problem */
   Bool                  need_startval,      /**< does ZIMPL provides a primal solution candidate */
   void*                 user_data           /**< user data which was previously passed to ZIMPL */
   )
{  /*lint --e{715}*/
   SCIP* scip;
   SCIP_READERDATA* readerdata;
   SCIP_Bool usestartsol;

   readerdata = (SCIP_READERDATA*)user_data;
   assert(readerdata != NULL);

   scip = readerdata->scip;
   assert(scip != NULL);

   /* create problem */
   SCIP_CALL_ABORT( SCIPcreateProb(scip, name, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* check if are interested in the primal solution candidate */
   SCIP_CALL_ABORT( SCIPgetBoolParam(scip, "reading/zplreader/usestartsol", &usestartsol) );

   if( usestartsol )
   {
      /* create primal solution */
      SCIP_CALL_ABORT( SCIPcreateSol(scip, &readerdata->sol, NULL) );
      readerdata->valid = TRUE;
   }

   /* return the reader data pointer to receive it all other ZIMPL call backs */
   return (Lps*) readerdata;
}

/** free storage for mathematical program. xlp_free() is the last xlpglue routine that will be called by Zimpl */
void xlp_free(
   Lps*                  data                /**< pointer to reader data */
   )
{  /*lint --e{715}*/
   /* nothing to be done here */
}

/** does there already exists a constraint with the given name? */ 
Bool xlp_conname_exists(
   const Lps*            data,               /**< pointer to reader data */
   const char*           name                /**< constraint name to check */
   )
{
   SCIP_READERDATA* readerdata;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

   /* check if constraint with the given name already exists */
   return (SCIPfindCons(readerdata->scip, name) != NULL);
}


#ifdef WITH_EXACTSOLVE
/** store a linear constraint in the exactlp constraint matrix and checks whether it has to be split into lhs and rhs
 *  constraints for building an FP-relaxation
 */
static
Bool storeLinearConstraint(
   Lps*                  data,               /**< pointer to reader data */
   ConType               type,               /**< constraint type (LHS, RHS, EQUAL, RANGE, etc) */
   const Term*           term                /**< term to use */
   )
{
   SCIP_READERDATA* readerdata;
   Bool split;
   mpq_t absval;
   int i;

   readerdata = (SCIP_READERDATA*)data;

   assert(readerdata != NULL);
   assert(readerdata->nconss >= 0);
   assert(readerdata->nnonz >= 0);
   assert(readerdata->beg[readerdata->nconss] == readerdata->nnonz);
   assert(readerdata->len[readerdata->nconss] == 0);

#ifdef CREATEEXACTLPCONS_OUT
   SCIPdebugMessage("zimpl reader: store linear cons with <%d> vars in matrix, nconss<%d>, next matrixbeg<%d>, next matrixlen<%d>\n",
      term_get_elements(term), readerdata->nconss, readerdata->beg[readerdata->nconss],
      readerdata->len[readerdata->nconss]);
#endif

   /* store coefficients of constraint */
   split = FALSE;
   mpq_init(absval);
   for( i = 0; i < term_get_elements(term); i++ )
   {
      SCIP_VAR* scipvar;
      int varidx;

      /* coefficient of variable */
      assert(!numb_equal(mono_get_coeff(term_get_element(term, i)), numb_zero()));
      assert(mono_is_linear(term_get_element(term, i)));
      numb_get_mpq(mono_get_coeff(term_get_element(term, i)), readerdata->val[readerdata->nnonz]);

      /* when we use an FP-relaxation for calculating dual bounds, we have to split the row into lhs and rhs part if
       *  - lhs and rhs of constraints are not -inf and inf, respectively, and
       *  - current coefficient is not FP-representable
       */
      if( type == CON_RANGE || type == CON_EQUAL )
         split = split || !mpqIsReal(readerdata->scip, readerdata->val[readerdata->nnonz]);

      /* index of variable */
      scipvar = (SCIP_VAR*)mono_get_var(term_get_element(term, i), 0);
      varidx = SCIPvarGetIndex(scipvar);
      assert(SCIPvarGetIndex(readerdata->vars[varidx]) == varidx);
      assert(readerdata->vars[varidx] == scipvar);
      readerdata->ind[readerdata->nnonz] = varidx;

#ifdef CREATEEXACTLPCONS_OUT
      {
         char s[SCIP_MAXSTRLEN];

         gmp_snprintf(s, "   i=%d: current nnonz<%d>, var<%s>, approx coef<%g>, stored exact coef<%Qd> (split:%d)\n",
            readerdata->nnonz, i, SCIPvarGetName(scipvar), numb_todbl(mono_get_coeff(term_get_element(term, i))),
            readerdata->val[readerdata->nnonz], split);
         SCIPdebugMessage(s);
      }
#endif

      /* update number of integral nonzero coefficients */
      if( mpqIsIntegral(readerdata->val[readerdata->nnonz]) )
         readerdata->nintegral++;

      /* update minimum and maximal absolute nonzero constraint matrix, lhs, or rhs entry */
      mpq_abs(absval, readerdata->val[readerdata->nnonz]);
      if( mpq_cmp(absval, readerdata->maxabsval) > 0 )
         mpq_set(readerdata->maxabsval, absval);
      if( mpq_cmp(absval, readerdata->minabsval) < 0 )
         mpq_set(readerdata->minabsval, absval);

      /* update number of nonzero coefficients */
      readerdata->nnonz++;
   }
   mpq_clear(absval);

   /* update number of constraints */
   readerdata->nconss++;

   /* initialize start index and length of next constraint */
   readerdata->beg[readerdata->nconss] = readerdata->nnonz;
   readerdata->len[readerdata->nconss] = 0;
   assert(readerdata->beg[readerdata->nconss-1] <= readerdata->beg[readerdata->nconss]);

   /* set length of constraint */
   readerdata->len[readerdata->nconss-1] = readerdata->beg[readerdata->nconss] - readerdata->beg[readerdata->nconss-1];

   return split;
}
#endif

/** method creates a constraint and is called directly from ZIMPL
 *
 *  @note this method is used by ZIMPL from version 3.00; 
 */
Bool xlp_addcon_term(
   Lps*                  data,               /**< pointer to reader data */
   const char*           name,               /**< constraint name */
   ConType               type,               /**< constraint type (LHS, RHS, EQUAL, RANGE, etc) */
   const Numb*           lhs,                /**< left hand side */
   const Numb*           rhs,                /**< right hand side */
   unsigned int          flags,              /**< special constraint flags, see ratlptypes.h */
   const Term*           term                /**< term to use */
   )
{
   SCIP* scip;
   SCIP_READERDATA* readerdata;
   SCIP_CONS* cons;
   SCIP_Real sciplhs;
   SCIP_Real sciprhs;
   SCIP_Bool initial;
   SCIP_Bool separate;
   SCIP_Bool enforce;
   SCIP_Bool check;
   SCIP_Bool propagate;
   SCIP_Bool local;
   SCIP_Bool modifiable;
   SCIP_Bool dynamic;
   SCIP_Bool removable;
   SCIP_Bool usercut;
   SCIP_Bool lazycut;
#ifdef WITH_EXACTSOLVE
   SCIP_Bool lhsgiven;
   SCIP_Bool rhsgiven;
   mpq_t absval;
#endif

   int  i;
   int  maxdegree;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

   scip = readerdata->scip;
   assert(scip != NULL);

#ifdef WITH_EXACTSOLVE
   /* reallocate and initialize constraint information; ranged constraints might be splitted into two constraints */
   assert(readerdata->nconss <= readerdata->consssize);
   if( readerdata->nconss + 3 > readerdata->consssize )
   {
      readerdata->consssize = MAX(2 * readerdata->consssize, readerdata->consssize + 3);

      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip, &readerdata->beg, readerdata->consssize) );
      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip, &readerdata->len, readerdata->consssize) );
      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip, &readerdata->lhs, readerdata->consssize) );
      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip, &readerdata->rhs, readerdata->consssize) );
      for( i = readerdata->nconss; i < readerdata->consssize; ++i )
      {
         mpq_init(readerdata->lhs[i]);
         mpq_init(readerdata->rhs[i]);
      }
   }

   /* reallocate and initialize matrix information; ranged constraints might be splitted into two constraints */
   assert(readerdata->nnonz <= readerdata->nonzsize);
   if( readerdata->nnonz + (2 * term_get_elements(term)) > readerdata->nonzsize )
   {
      readerdata->nonzsize = MAX(2 * readerdata->nonzsize, readerdata->nonzsize + (2 * term_get_elements(term)));
      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip, &readerdata->val, readerdata->nonzsize) );
      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip, &readerdata->ind, readerdata->nonzsize) );
      for( i = readerdata->nnonz; i < readerdata->nonzsize; ++i )
      {
         mpq_init(readerdata->val[i]);
      }
   }
#endif

#ifdef WITH_EXACTSOLVE
   lhsgiven = FALSE;
   rhsgiven = FALSE;

   /* get exact lhs and rhs */
   switch( type )
   {
   case CON_FREE:
      mpq_set(readerdata->lhs[readerdata->nconss], *negInfinity(readerdata->conshdlrdata));
      mpq_set(readerdata->rhs[readerdata->nconss], *posInfinity(readerdata->conshdlrdata));
      break;
   case CON_LHS:
      numb_get_mpq(lhs, readerdata->lhs[readerdata->nconss]);
      mpq_set(readerdata->rhs[readerdata->nconss], *posInfinity(readerdata->conshdlrdata));
      lhsgiven = TRUE;
      break;
   case CON_RHS:
      mpq_set(readerdata->lhs[readerdata->nconss], *negInfinity(readerdata->conshdlrdata));
      numb_get_mpq(rhs, readerdata->rhs[readerdata->nconss]);
      rhsgiven = TRUE;
      break;
   case CON_RANGE:
      numb_get_mpq(lhs, readerdata->lhs[readerdata->nconss]);
      numb_get_mpq(rhs, readerdata->rhs[readerdata->nconss]);
      lhsgiven = TRUE;
      rhsgiven = TRUE;
      break;
   case CON_EQUAL:
      numb_get_mpq(lhs, readerdata->lhs[readerdata->nconss]);
      numb_get_mpq(rhs, readerdata->rhs[readerdata->nconss]);
      assert(mpq_equal(readerdata->lhs[readerdata->nconss], readerdata->rhs[readerdata->nconss]) != 0);
      lhsgiven = TRUE;
      rhsgiven = TRUE;
      break;
   default:
      SCIPwarningMessage(scip, "invalid constraint type <%d> in ZIMPL callback xlp_addcon()\n", type);
      numb_get_mpq(lhs, readerdata->lhs[readerdata->nconss]);
      numb_get_mpq(rhs, readerdata->rhs[readerdata->nconss]);
      readerdata->readerror = TRUE;
      break;
   }

   /* update minimum and maximum absolute nonzero entry of constraint matrix, lhs and rhs vector */
   mpq_init(absval);
   if( lhsgiven && mpq_sgn(readerdata->lhs[readerdata->nconss]) != 0 )
   {
      mpq_abs(absval, readerdata->lhs[readerdata->nconss]);
      if( mpq_cmp(absval, readerdata->maxabsval) > 0 )
         mpq_set(readerdata->maxabsval, absval);
      if( mpq_cmp(absval, readerdata->minabsval) < 0 )
         mpq_set(readerdata->minabsval, absval);
   }
   if( rhsgiven && mpq_sgn(readerdata->rhs[readerdata->nconss]) != 0 )
   {
      mpq_abs(absval, readerdata->rhs[readerdata->nconss]);
      if( mpq_cmp(absval, readerdata->maxabsval) > 0 )
         mpq_set(readerdata->maxabsval, absval);
      if( mpq_cmp(absval, readerdata->minabsval) < 0 )
         mpq_set(readerdata->minabsval, absval);
   }
   mpq_clear(absval);
#else
   /* get double-precision lhs and rhs */
   switch( type )
   {
   case CON_FREE:
      sciplhs = -SCIPinfinity(scip);
      sciprhs = SCIPinfinity(scip);
      break;
   case CON_LHS:
      sciplhs = (SCIP_Real)numb_todbl(lhs);
      sciprhs = SCIPinfinity(scip);
      break;
   case CON_RHS:
      sciplhs = -SCIPinfinity(scip);
      sciprhs = (SCIP_Real)numb_todbl(rhs);
      break;
   case CON_RANGE:
      sciplhs = (SCIP_Real)numb_todbl(lhs);
      sciprhs = (SCIP_Real)numb_todbl(rhs);
      break;
   case CON_EQUAL:
      sciplhs = (SCIP_Real)numb_todbl(lhs);
      sciprhs = (SCIP_Real)numb_todbl(rhs);
      assert(sciplhs == sciprhs);  /*lint !e777*/
      break;
   default:
      SCIPwarningMessage(scip, "invalid constraint type <%d> in ZIMPL callback xlp_addcon()\n", type);
      sciplhs = (SCIP_Real)numb_todbl(lhs);
      sciprhs = (SCIP_Real)numb_todbl(rhs);
      readerdata->readerror = TRUE;
      break;
   }
#endif

   cons = NULL;

   /* default values */
   initial = TRUE;
   separate = TRUE;
   propagate = TRUE;
   enforce = TRUE;
   check = TRUE;
   removable = FALSE;
   local = FALSE;
   modifiable = FALSE;
   dynamic = FALSE;

   usercut = (flags & LP_FLAG_CON_SEPAR) != 0;
   lazycut = (flags & LP_FLAG_CON_CHECK) != 0;

   /* evaluate constraint flags */
   if( usercut && lazycut )
   {
      initial = FALSE;
      separate = TRUE;
      check = TRUE;
   }
   else if( usercut )
   {
      initial = FALSE;
      separate = TRUE;
      check = FALSE;
   }
   else if( lazycut )
   {
      initial = FALSE;
      separate = FALSE;
      check = TRUE;
   }

   maxdegree = term_get_degree(term);

   if (maxdegree <= 1)
   {
      /* if the constraint gives an indicator constraint */
      if ( flags & LP_FLAG_CON_INDIC )
      {
         Bool lhsIndCons = FALSE;  /* generate lhs form for indicator constraints */
         Bool rhsIndCons = FALSE;  /* generate rhs form for indicator constraints */

#ifdef WITH_EXACTSOLVE
         SCIPerrorMessage("xpl_addcon_term: exact version for indicator constraints not supported\n");
         return TRUE;
#endif

         /* currently indicator constraints can only handle "<=" constraints */
         switch( type )
         {
         case CON_LHS:
            lhsIndCons = TRUE;
            break;
         case CON_RHS:
            rhsIndCons = TRUE;
            break;
         case CON_RANGE:
         case CON_EQUAL:
            lhsIndCons = TRUE;
            rhsIndCons = TRUE;
            break;
         case CON_FREE:
            /*lint -fallthrough*/
         default:
            SCIPerrorMessage("invalid constraint type <%d> in ZIMPL callback xlp_addcon()\n", type);
            readerdata->readerror = TRUE;
            break;
         }

         /* insert lhs form of indicator */
         if ( lhsIndCons )
         {
            SCIP_CALL_ABORT( SCIPcreateConsIndicator(scip, &cons, name, NULL, 0, NULL, NULL, -sciplhs,
                  initial, separate, enforce, check, propagate, local, dynamic, removable, FALSE) );
            SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );

            for( i = 0; i < term_get_elements(term); i++ )
            {
               SCIP_VAR* scipvar;
               SCIP_Real scipval;
               const Mono* mono = term_get_element(term, i);
               MFun mfun;

               scipvar = (SCIP_VAR*)mono_get_var(mono, 0);

               /* check whether variable is the binary variable */
               mfun = mono_get_function(mono);
               if (mfun == MFUN_TRUE || mfun == MFUN_FALSE)
               {
                  scipvar = (SCIP_VAR*)mono_get_var(mono, 0);
                  SCIP_CALL_ABORT( SCIPsetBinaryVarIndicator(scip, cons, scipvar) );
               }
               else
               {
                  assert(!numb_equal(mono_get_coeff(mono), numb_zero()));
                  assert(mono_is_linear(mono));

                  scipval = -numb_todbl(mono_get_coeff(mono));
                  SCIP_CALL_ABORT( SCIPaddVarIndicator(scip, cons, scipvar, scipval) );
               }
            }
         }

         /* insert rhs form of indicator */
         if ( rhsIndCons )
         {
            SCIP_CALL_ABORT( SCIPcreateConsIndicator(scip, &cons, name, NULL, 0, NULL, NULL, sciprhs,
                  initial, separate, enforce, check, propagate, local, dynamic, removable, FALSE) );
            SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );

            for( i = 0; i < term_get_elements(term); i++ )
            {
               SCIP_VAR* scipvar;
               SCIP_Real scipval;
               const Mono* mono = term_get_element(term, i);
               MFun mfun;

               scipvar = (SCIP_VAR*)mono_get_var(mono, 0);

               /* check whether variable is the binary variable */
               mfun = mono_get_function(mono);
               if (mfun == MFUN_TRUE || mfun == MFUN_FALSE)
               {
                  scipvar = (SCIP_VAR*)mono_get_var(mono, 0);
                  SCIP_CALL_ABORT( SCIPsetBinaryVarIndicator(scip, cons, scipvar) );
               }
               else
               {
                  assert(!numb_equal(mono_get_coeff(mono), numb_zero()));
                  assert(mono_is_linear(mono));

                  scipval = numb_todbl(mono_get_coeff(mono));
                  SCIP_CALL_ABORT( SCIPaddVarIndicator(scip, cons, scipvar, scipval) );
               }
            }
         }
      }
      else
      {
#ifdef WITH_EXACTSOLVE
         Bool split;

         /* store constraint in the exactlp constraint matrix */
         split = storeLinearConstraint(data, type, term);

         /* update number of linear constraints that need to be split in case of an FP relaxation */
         if( split )
            readerdata->nsplitconss++;

         /* for an FP-relaxation, split constraint if necessary */
         if( SCIPuseFPRelaxation(scip) && split )
         {
            assert(type == CON_RANGE || type == CON_EQUAL);

            /* change constraint just stored to an rhs constraint */
            mpq_set(readerdata->lhs[readerdata->nconss-1], *negInfinity(readerdata->conshdlrdata));

            /* store the constraint again as lhs constraint */
            numb_get_mpq(lhs, readerdata->lhs[readerdata->nconss]);
            mpq_set(readerdata->rhs[readerdata->nconss], *posInfinity(readerdata->conshdlrdata));
            split = storeLinearConstraint(data, type, term);
            assert(split);
         }
#else
         SCIP_CALL_ABORT( SCIPcreateConsLinear(scip, &cons, name, 0, NULL, NULL, sciplhs, sciprhs,
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
         SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );

         for( i = 0; i < term_get_elements(term); i++ )
         {
            SCIP_VAR* scipvar;
            SCIP_Real scipval;

            assert(!numb_equal(mono_get_coeff(term_get_element(term, i)), numb_zero()));
            assert(mono_is_linear(term_get_element(term, i)));

            scipvar = (SCIP_VAR*)mono_get_var(term_get_element(term, i), 0);
            scipval = numb_todbl(mono_get_coeff(term_get_element(term, i)));

            SCIP_CALL_ABORT( SCIPaddCoefLinear(scip, cons, scipvar, scipval) );
         }
#endif
      }
   }
   else if (maxdegree == 2)
   {
      int        nlinvars;
      int        nquadterms;
      SCIP_VAR** linvars;
      SCIP_VAR** quadvar1;
      SCIP_VAR** quadvar2;
      SCIP_Real* lincoefs;
      SCIP_Real* quadcoefs;
      Mono*      monom;

#ifdef WITH_EXACTSOLVE
      SCIPerrorMessage("xpl_addcon_term: exact version for degree == 2 not supported\n");
      return TRUE;
#endif

      nlinvars   = 0;
      nquadterms = 0;

      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &linvars,   term_get_elements(term)) );
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &quadvar1,  term_get_elements(term)) );
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &quadvar2,  term_get_elements(term)) );
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &lincoefs,  term_get_elements(term)) );
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &quadcoefs, term_get_elements(term)) );

      for( i = 0; i < term_get_elements(term); ++i )
      {
         monom = term_get_element(term, i);
         assert(!numb_equal(mono_get_coeff(monom), numb_zero()));
         assert(mono_get_degree(monom) <= 2);
         assert(mono_get_degree(monom) > 0);
         if (mono_get_degree(monom) == 1)
         {
            linvars [nlinvars] = (SCIP_VAR*)mono_get_var(monom, 0);
            lincoefs[nlinvars] = numb_todbl(mono_get_coeff(monom));
            ++nlinvars;
         }
         else
         {
            assert(mono_get_degree(monom) == 2);
            quadvar1 [nquadterms] = (SCIP_VAR*)mono_get_var(monom, 0);
            quadvar2 [nquadterms] = (SCIP_VAR*)mono_get_var(monom, 1);
            quadcoefs[nquadterms] = numb_todbl(mono_get_coeff(monom));
            ++nquadterms;
         }
      }

      SCIP_CALL_ABORT( SCIPcreateConsQuadratic(scip, &cons, name, nlinvars, linvars, lincoefs, nquadterms, quadvar1, quadvar2, quadcoefs, sciplhs, sciprhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );
      SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );

      SCIPfreeBufferArray(scip, &linvars);
      SCIPfreeBufferArray(scip, &quadvar1);
      SCIPfreeBufferArray(scip, &quadvar2);
      SCIPfreeBufferArray(scip, &lincoefs);
      SCIPfreeBufferArray(scip, &quadcoefs);
   }
   else
   {
      SCIP_VAR** polyvars;
      int        npolyvars;
      int        polyvarssize;
      SCIP_HASHMAP* polyvarmap;
      SCIP_VAR** vars;
      int        nvars;
      int        varssize;
      SCIP_HASHMAP* varmap;
      SCIP_EXPRDATA_MONOMIAL** simplemonomials;
      int        nsimplemonomials;
      int        simplemonomialssize;
      SCIP_EXPR** extramonomials;
      SCIP_Real* extracoefs;
      int        nextramonomials;
      int        extramonomialssize;
      Mono*      monomial;
      SCIP_Bool  fail;
      int varpos;
      int j;

#ifdef WITH_EXACTSOLVE
      SCIPerrorMessage("xpl_addcon_term: exact version for degree > 2 not supported\n");
      return TRUE;
#endif

      fail = FALSE;

      vars = NULL;
      nvars = 0;
      varssize = 0;
      varmap = NULL;

      polyvars = NULL;
      npolyvars = 0;
      polyvarssize = 0;

      simplemonomials = NULL;
      nsimplemonomials = 0;
      simplemonomialssize = 0;

      extramonomials = NULL;
      extracoefs = NULL;
      nextramonomials = 0;
      extramonomialssize = 0;

      SCIP_CALL_ABORT( SCIPhashmapCreate(&varmap, SCIPblkmem(scip), SCIPcalcMemGrowSize(scip, 10)) );
      SCIP_CALL_ABORT( SCIPhashmapCreate(&polyvarmap, SCIPblkmem(scip), SCIPcalcMemGrowSize(scip, 10)) );

      for( i = 0; i < term_get_elements(term); ++i )
      {
         monomial = term_get_element(term, i);
         assert(monomial != NULL);
         assert(!numb_equal(mono_get_coeff(monomial), numb_zero()));
         assert(mono_get_degree(monomial) > 0);

         if( mono_get_function(monomial) == MFUN_NONE )
         {
            /* nonlinear monomial without extra function around it */
            SCIP_Real one;

            one = 1.0;

            /* create SCIP monomial */
            if( simplemonomialssize == 0 )
            {
               simplemonomialssize = SCIPcalcMemGrowSize(scip, 1);
               SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &simplemonomials, simplemonomialssize) );
            }
            else if( simplemonomialssize < nsimplemonomials + 1 )
            {
               simplemonomialssize = SCIPcalcMemGrowSize(scip, nsimplemonomials+1);
               SCIP_CALL_ABORT( SCIPreallocBufferArray(scip, &simplemonomials, simplemonomialssize) );
            }
            assert(simplemonomials != NULL);
            SCIP_CALL_ABORT( SCIPexprCreateMonomial(SCIPblkmem(scip), &simplemonomials[nsimplemonomials], numb_todbl(mono_get_coeff(monomial)), 0, NULL, NULL) );

            for( j = 0; j < mono_get_degree(monomial); ++j )
            {
               /* get variable index in polyvars; add to polyvars if not existing yet */
               if( !SCIPhashmapExists(polyvarmap, (void*)mono_get_var(monomial, j)) )  /*lint !e826*/
               {
                  if( polyvarssize == 0 )
                  {
                     polyvarssize = SCIPcalcMemGrowSize(scip, 1);
                     SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &polyvars, polyvarssize) );
                  }
                  else if( polyvarssize < npolyvars + 1 )
                  {
                     polyvarssize = SCIPcalcMemGrowSize(scip, npolyvars+1);
                     SCIP_CALL_ABORT( SCIPreallocBufferArray(scip, &polyvars, polyvarssize) );
                  }
                  assert(polyvars != NULL);

                  polyvars[npolyvars] = (SCIP_VAR*)mono_get_var(monomial, j);  /*lint !e826*/
                  ++npolyvars;
                  varpos = npolyvars-1;
                  SCIP_CALL_ABORT( SCIPhashmapInsert(polyvarmap, (void*)mono_get_var(monomial, j), (void*)(size_t)varpos) );  /*lint !e826*/
               }
               else
               {
                  varpos = (int)(size_t)SCIPhashmapGetImage(polyvarmap, (void*)mono_get_var(monomial, j));  /*lint !e826*/
               }
               assert(polyvars != NULL);
               assert(polyvars[varpos] == (SCIP_VAR*)mono_get_var(monomial, j));

               SCIP_CALL_ABORT( SCIPexprAddMonomialFactors(SCIPblkmem(scip), simplemonomials[nsimplemonomials], 1, &varpos, &one) );
            }
            SCIPexprMergeMonomialFactors(simplemonomials[nsimplemonomials], 0.0);

            ++nsimplemonomials;
         }
         else
         {
            /* nonlinear monomial with extra function around it, put into new expression */
            SCIP_EXPR** children;
            SCIP_EXPR* expr;
            SCIP_EXPROP op;

            switch( mono_get_function(monomial) )
            {
            case MFUN_SQRT:
               op = SCIP_EXPR_SQRT;
               break;
            case MFUN_LOG:
               op = SCIP_EXPR_LOG;
               break;
            case MFUN_EXP:
               op = SCIP_EXPR_EXP;
               break;
            default:
               SCIPerrorMessage("ZIMPL function %d not supported\n", mono_get_function(monomial));
               fail = TRUE;
               break;
            }  /*lint !e788*/
            if( fail )
               break;

            if( extramonomialssize == 0 )
            {
               extramonomialssize = SCIPcalcMemGrowSize(scip, 1);
               SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &extramonomials, extramonomialssize) );
               SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &extracoefs,  extramonomialssize) );
            }
            else if( extramonomialssize < nextramonomials + 1 )
            {
               extramonomialssize = SCIPcalcMemGrowSize(scip, nextramonomials+1);
               SCIP_CALL_ABORT( SCIPreallocBufferArray(scip, &extramonomials, extramonomialssize) );
               SCIP_CALL_ABORT( SCIPreallocBufferArray(scip, &extracoefs,  extramonomialssize) );
            }
            assert(extracoefs != NULL);
            assert(extramonomials != NULL);
            extracoefs[nextramonomials] = numb_todbl(mono_get_coeff(monomial));

            /* create children expressions */
            SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &children, mono_get_degree(monomial)) );
            for( j = 0; j < mono_get_degree(monomial); ++j )
            {
               /* get variable index in vars; add to vars if not existing yet */
               if( !SCIPhashmapExists(varmap, (void*)mono_get_var(monomial, j)) )  /*lint !e826*/
               {
                  if( varssize == 0 )
                  {
                     varssize = SCIPcalcMemGrowSize(scip, 1);
                     SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &vars, varssize) );
                  }
                  else if( varssize < nvars + 1 )
                  {
                     varssize = SCIPcalcMemGrowSize(scip, nvars+1);
                     SCIP_CALL_ABORT( SCIPreallocBufferArray(scip, &vars, varssize) );
                  }
                  assert(vars != NULL);

                  vars[nvars] = (SCIP_VAR*)mono_get_var(monomial, j);  /*lint !e826*/
                  ++nvars;
                  varpos = nvars-1;
                  SCIP_CALL_ABORT( SCIPhashmapInsert(varmap, (void*)mono_get_var(monomial, j), (void*)(size_t)varpos) );  /*lint !e826*/
               }
               else
               {
                  varpos = (int)(size_t)SCIPhashmapGetImage(varmap, (void*)mono_get_var(monomial, j));  /*lint !e826*/
               }
               assert(vars != NULL);
               assert(vars[varpos] == (SCIP_VAR*)mono_get_var(monomial, j));

               SCIP_CALL_ABORT( SCIPexprCreate(SCIPblkmem(scip), &children[j], SCIP_EXPR_VARIDX, varpos) );
            }

            /* create expression for product of variables */
            SCIP_CALL_ABORT( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_PRODUCT, mono_get_degree(monomial), children) );
            /* create expression for function of product of variables */
            SCIP_CALL_ABORT( SCIPexprCreate(SCIPblkmem(scip), &extramonomials[nextramonomials], op, expr) );  /*lint !e644*/

            ++nextramonomials;
         }
      }

      if( !fail )
      {
         SCIP_EXPRTREE* exprtree;
         SCIP_EXPR* polynomial;
         SCIP_EXPR** children;
         int nchildren;

         assert(polyvars != NULL || npolyvars == 0);

         nchildren = npolyvars + nextramonomials;
         SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &children, nchildren) );
         /* add polynomial variables to vars
          * create children expressions for polynomial variables
          */
         for( i = 0; i < npolyvars; ++i )
         {
            /* get variable index in vars; add to vars if not existing yet */
            if( !SCIPhashmapExists(varmap, (void*)polyvars[i]) )  /*lint !e613*/
            {
               if( varssize == 0 )
               {
                  varssize = SCIPcalcMemGrowSize(scip, 1);
                  SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &vars, varssize) );
               }
               else if( varssize < nvars + 1 )
               {
                  varssize = SCIPcalcMemGrowSize(scip, nvars+1);
                  SCIP_CALL_ABORT( SCIPreallocBufferArray(scip, &vars, varssize) );
               }
               assert(vars != NULL);

               vars[nvars] = polyvars[i];  /*lint !e613*/
               ++nvars;
               varpos = nvars-1;
               SCIP_CALL_ABORT( SCIPhashmapInsert(varmap, (void*)polyvars[i], (void*)(size_t)varpos) );  /*lint !e613*/
            }
            else
            {
               varpos = (int)(size_t)SCIPhashmapGetImage(varmap, (void*)polyvars[i]);  /*lint !e613*/
            }
            assert(vars[varpos] == polyvars[i]);  /*lint !e613*/

            SCIP_CALL_ABORT( SCIPexprCreate(SCIPblkmem(scip), &children[i], SCIP_EXPR_VARIDX, varpos) );  /*lint !e866*/
         }

         /* add simple monomials as additional children */
         BMScopyMemoryArray(&children[npolyvars], extramonomials, nextramonomials);  /*lint !e866*/

         assert(extracoefs     != NULL || nextramonomials == 0);
         assert(extramonomials != NULL || nextramonomials == 0);

         /* create polynomial expression including simple monomials */
         SCIP_CALL_ABORT( SCIPexprCreatePolynomial(SCIPblkmem(scip), &polynomial, nchildren, children, nsimplemonomials, simplemonomials, 0.0, FALSE) );
         /* add extra monomials */
         for( i = 0; i < nextramonomials; ++i )
         {
            SCIP_EXPRDATA_MONOMIAL* monomialdata;
            int childidx;
            SCIP_Real exponent;

            childidx = npolyvars + i;
            exponent = 1.0;
            SCIP_CALL_ABORT( SCIPexprCreateMonomial(SCIPblkmem(scip), &monomialdata, extracoefs[i], 1, &childidx, &exponent) );  /*lint !e613*/
            SCIP_CALL_ABORT( SCIPexprAddMonomials(SCIPblkmem(scip), polynomial, 1, &monomialdata, FALSE) );
         }

         SCIPfreeBufferArray(scip, &children);

         /* create expression tree */
         SCIP_CALL_ABORT( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, polynomial, nvars, 0, NULL) );
         SCIP_CALL_ABORT( SCIPexprtreeSetVars(exprtree, nvars, vars) );

         /* create constraint */
         SCIP_CALL_ABORT( SCIPcreateConsNonlinear(scip, &cons, name, 0, NULL, NULL, 1, &exprtree, NULL, sciplhs, sciprhs,
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
         SCIP_CALL_ABORT( SCIPexprtreeFree(&exprtree) );
         SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );
      }

      /* free memory */
      SCIPhashmapFree(&varmap);
      SCIPfreeBufferArrayNull(scip, &vars);
      SCIPhashmapFree(&polyvarmap);
      SCIPfreeBufferArrayNull(scip, &polyvars);
      SCIPfreeBufferArrayNull(scip, &simplemonomials);
      SCIPfreeBufferArrayNull(scip, &extramonomials);
      SCIPfreeBufferArrayNull(scip, &extracoefs);

      if( fail )
         return TRUE;
   }

   if( cons != NULL )
   {
      SCIP_CALL_ABORT( SCIPreleaseCons(scip, &cons) );
   }

   return FALSE;
}

/** method adds a variable; is called directly by ZIMPL */
Var* xlp_addvar(
   Lps*                  data,               /**< pointer to reader data */
   const char*           name,               /**< variable name */
   VarClass              usevarclass,        /**< variable type */
   const Bound*          lower,              /**< lower bound */
   const Bound*          upper,              /**< upper bound */
   const Numb*           priority,           /**< branching priority */
   const Numb*           startval            /**< start value for the variable within in the start solution */
   )
{  /*lint --e{715}*/
   SCIP* scip;
   SCIP_READERDATA* readerdata;
   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_VARTYPE vartype;
   SCIP_Bool initial;
   SCIP_Bool removable;
   SCIP_Bool dynamiccols;
   int branchpriority;
   Var* zplvar;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

   scip = readerdata->scip;
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPgetBoolParam(scip, "reading/zplreader/dynamiccols", &dynamiccols) );

#ifdef WITH_EXACTSOLVE
   /* reallocate and initialize variable specific information */
   assert(readerdata->nvars <= readerdata->varssize);
   if( readerdata->nvars == readerdata->varssize )
   {
      int i;

      readerdata->varssize *= 2;
      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip_, &readerdata->vars, readerdata->varssize) );
      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip_, &readerdata->lb, readerdata->varssize) );
      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip_, &readerdata->ub, readerdata->varssize) );
      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip_, &readerdata->obj, readerdata->varssize) );
      for( i = readerdata->nvars; i < readerdata->varssize; ++i )
      {
         mpq_init(readerdata->lb[i]);
         mpq_init(readerdata->ub[i]);
         mpq_init(readerdata->obj[i]);
         mpq_set_si(readerdata->obj[i], 0, 1);
      }
   }
   assert(readerdata->nvars < readerdata->varssize);

   /* get exact lower bounds for exactlp constraint handler and safe FP-values for FP-problem */
   switch( bound_get_type(lower) )
   {
   case BOUND_VALUE:
      numb_get_mpq(bound_get_value(lower), readerdata->lb[readerdata->nvars]);
      lb = mpqGetRealRelax(scip, readerdata->lb[readerdata->nvars], GMP_RNDD);
      break;
   case BOUND_INFTY:
      mpq_set(readerdata->lb[readerdata->nvars], *posInfinity(readerdata->conshdlrdata));
      lb = SCIPinfinity(scip);
      break;
   case BOUND_MINUS_INFTY:
      mpq_set(readerdata->lb[readerdata->nvars], *negInfinity(readerdata->conshdlrdata));
      lb = -SCIPinfinity(scip);
      break;
   case BOUND_ERROR:
   default:
      SCIPerrorMessage("invalid lower bound type <%d> in ZIMPL reader\n", bound_get_type(lower));
      mpq_set_si(readerdata->lb[readerdata->nvars], 0, 1);
      lb = 0.0;
      break;
   }

   /* get exact upper bounds for exactlp constraint handler and safe FP-values for FP-problem */
   switch( bound_get_type(upper) )
   {
   case BOUND_VALUE:
      numb_get_mpq(bound_get_value(upper), readerdata->ub[readerdata->nvars]);
      ub = mpqGetRealRelax(scip, readerdata->ub[readerdata->nvars], GMP_RNDU);
      break;
   case BOUND_INFTY:
      mpq_set(readerdata->ub[readerdata->nvars], *posInfinity(readerdata->conshdlrdata));
      ub = SCIPinfinity(scip);
      break;
   case BOUND_MINUS_INFTY:
      mpq_set(readerdata->ub[readerdata->nvars], *negInfinity(readerdata->conshdlrdata));
      ub = -SCIPinfinity(scip);
      break;
   case BOUND_ERROR:
   default:
      SCIPerrorMessage("invalid upper bound type <%d> in ZIMPL reader\n", bound_get_type(upper));
      mpq_set_si(readerdata->ub[readerdata->nvars], 0, 1);
      ub = 0.0;
      break;
   }
#else
   switch( bound_get_type(lower) )
   {
   case BOUND_VALUE:
      lb = (SCIP_Real)numb_todbl(bound_get_value(lower));
      break;
   case BOUND_INFTY:
      lb = SCIPinfinity(scip);
      break;
   case BOUND_MINUS_INFTY:
      lb = -SCIPinfinity(scip);
      break;
   case BOUND_ERROR:
   default:
      SCIPerrorMessage("invalid lower bound type <%d> in ZIMPL reader\n", bound_get_type(lower));
      lb = 0.0;
      break;
   }

   switch( bound_get_type(upper) )
   {
   case BOUND_VALUE:
      ub = (SCIP_Real)numb_todbl(bound_get_value(upper));
      break;
   case BOUND_INFTY:
      ub = SCIPinfinity(scip);
      break;
   case BOUND_MINUS_INFTY:
      ub = -SCIPinfinity(scip);
      break;
   case BOUND_ERROR:
   default:
      SCIPerrorMessage("invalid upper bound type <%d> in ZIMPL reader\n", bound_get_type(upper));
      ub = 0.0;
      break;
   }
#endif

   /* get variable type */
   switch( usevarclass )
   {
   case VAR_CON:
      vartype = SCIP_VARTYPE_CONTINUOUS;
      break;
   case VAR_INT:
      vartype = SCIP_VARTYPE_INTEGER;
      break;
   case VAR_IMP:
      vartype = SCIP_VARTYPE_IMPLINT;
      break;
   default:
      SCIPwarningMessage(scip, "invalid variable class <%d> in ZIMPL callback xlp_addvar()\n", usevarclass);
      vartype = SCIP_VARTYPE_CONTINUOUS;
      readerdata->readerror = TRUE;
      break;
   }

#ifdef WITH_EXACTSOLVE
   /* update number of variables with infinite or large bounds and integer variables with infinite bounds */
   if( bound_get_type(lower) == BOUND_MINUS_INFTY || bound_get_type(upper) == BOUND_INFTY  )
   {
      readerdata->ninfbounds++;
      if( vartype != SCIP_VARTYPE_CONTINUOUS )
         readerdata->ninfintbounds++;
   }
   else if( lb <= -LARGEBOUND || ub >= LARGEBOUND )
      readerdata->nlargebounds++;
#endif

   initial = !dynamiccols;
   removable = dynamiccols;

   /* create variable */
   SCIP_CALL_ABORT( SCIPcreateVar(scip, &var, name, lb, ub, 0.0, vartype, initial, removable, NULL, NULL, NULL, NULL, NULL) );

   /* add variable to the problem; we are releasing the variable later */
   SCIP_CALL_ABORT( SCIPaddVar(scip, var) );

   if( !numb_equal(priority, numb_unknown()) )
   {
      if( numb_is_int(priority) )
	 branchpriority = numb_toint(priority);
      else
      {
	 if( !readerdata->branchpriowarning )
	 {
	    SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL,
	       "ZIMPL reader: fractional branching priorities in input - rounding down to integer values\n");
	    readerdata->branchpriowarning = TRUE;
	 }
	 branchpriority = (int)numb_todbl(priority);
      }

      /* change the branching priority of the variable */
      SCIP_CALL_ABORT( SCIPchgVarBranchPriority(scip, var, branchpriority) );
   }

   /* check if we are willing to except a primal solution candidate */
   if( readerdata->valid )
   {
      /* if the number is unknown we have no valid primal solution candidate */
      if( numb_equal(startval, numb_unknown()) )
      {
         SCIPdebugMessage("primal solution candidate contains an unknown value for variable <%s>(%g)\n", 
            SCIPvarGetName(var), (SCIP_Real)numb_todbl(startval));
         readerdata->valid = FALSE;
      }
      else
      {
         assert(readerdata->sol != NULL);
         SCIPdebugMessage("change solution solution <%p>: <%s> = <%g>\n", 
            readerdata->sol, SCIPvarGetName(var), (SCIP_Real)numb_todbl(startval));

         /* set value within the primal solution candidate */
         SCIP_CALL_ABORT( SCIPsetSolVal(scip, readerdata->sol, var, (SCIP_Real)numb_todbl(startval)) );
      }
   }

   /* copy the variable pointer before we release the variable */
   zplvar = (Var*)var;

#ifdef WITH_EXACTSOLVE
   /* store variable for exactlp constraint handler */
   readerdata->vars[readerdata->nvars] = var;
   readerdata->nvars++;

#ifdef CREATEEXACTLPCONS_OUT
   SCIPdebugMessage("zimpl reader: added new variable");
   SCIPprintVar(scip_, var, NULL);
#endif
#endif

   /* release variable */
   SCIP_CALL_ABORT( SCIPreleaseVar(scip, &var) );

   return zplvar;
}

/** add a SOS constraint. Add a given a Zimpl term as an SOS constraint to the mathematical program */
Bool xlp_addsos_term(
   Lps*                  data,               /**< pointer to reader data */
   const char*           name,               /**< constraint name */
   SosType               type,               /**< SOS type */
   const Numb*           priority,           /**< priority */
   const Term*           term                /**< terms indicating sos */
   )
{
   /*lint --e{715}*/
   SCIP* scip;
   SCIP_READERDATA* readerdata;
   SCIP_CONS* cons;
   SCIP_Bool initial;
   SCIP_Bool separate;
   SCIP_Bool enforce;
   SCIP_Bool check;
   SCIP_Bool propagate;
   SCIP_Bool local;
   SCIP_Bool dynamic;
   SCIP_Bool removable;
   int i;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

#ifdef WITH_EXACTSOLVE
   SCIPerrorMessage("xlp_addsos_termr: exact version not supported.\n");
   readerdata->readerror = TRUE;
   return FALSE;
#endif

   scip = readerdata->scip;
   assert(scip != NULL);

   switch( type )
   {
   case SOS_TYPE1:
      initial = TRUE;
      separate = TRUE;
      enforce = TRUE;
      check = enforce;
      propagate = TRUE;
      local = FALSE;
      dynamic = FALSE;
      removable = dynamic;

      SCIP_CALL_ABORT( SCIPcreateConsSOS1(scip, &cons, name, 0, NULL, NULL,
            initial, separate, enforce, check, propagate, local, dynamic, removable, FALSE) );
      SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );

      for( i = 0; i < term_get_elements(term); i++ )
      {
         SCIP_VAR* var;
         SCIP_Real weight;

         assert( mono_is_linear(term_get_element(term, i)) );

         var = (SCIP_VAR*) mono_get_var(term_get_element(term, i), 0);
         weight = numb_todbl(mono_get_coeff(term_get_element(term, i)));

         SCIP_CALL_ABORT( SCIPaddVarSOS1(scip, cons, var, weight) );
      }
      SCIP_CALL_ABORT( SCIPreleaseCons(scip, &cons) );
      break;
   case SOS_TYPE2:
      initial = TRUE;
      separate = TRUE;
      enforce = TRUE;
      check = enforce;
      propagate = TRUE;
      local = FALSE;
      dynamic = FALSE;
      removable = dynamic;

      SCIP_CALL_ABORT( SCIPcreateConsSOS2(scip, &cons, name, 0, NULL, NULL, 
            initial, separate, enforce, check, propagate, local, dynamic, removable, FALSE) );
      SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );
      for( i = 0; i < term_get_elements(term); i++ )
      {
         SCIP_VAR* var;
         SCIP_Real weight;

         assert( mono_is_linear(term_get_element(term, i)) );

         var = (SCIP_VAR*) mono_get_var(term_get_element(term, i), 0);
         weight = numb_todbl(mono_get_coeff(term_get_element(term, i)));

         SCIP_CALL_ABORT( SCIPaddVarSOS2(scip, cons, var, weight) );
      }
      SCIP_CALL_ABORT( SCIPreleaseCons(scip, &cons) );
      break;
   case SOS_ERR:
      /*lint -fallthrough*/
   default:
      SCIPerrorMessage("invalid SOS type <%d> in ZIMPL callback xlp_addsos_term()\n", type);
      readerdata->readerror = TRUE;
      break;
   }

   return FALSE;
}

/** returns the variable name */
const char* xlp_getvarname(
   const Lps*            data,               /**< pointer to reader data */
   const Var*            var                 /**< variable */
   )
{
#ifndef NDEBUG
   SCIP* scip;
   SCIP_READERDATA* readerdata;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

   scip = readerdata->scip;
   assert(scip != NULL);
#endif

   return SCIPvarGetName((SCIP_VAR*)var);
}

/** return variable type */
VarClass xlp_getclass(
   const Lps*            data,               /**< pointer to reader data */
   const Var*            var                 /**< variable */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_VAR* scipvar;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

   scipvar = (SCIP_VAR*)var;
   switch( SCIPvarGetType(scipvar) )
   {
   case SCIP_VARTYPE_BINARY:
   case SCIP_VARTYPE_INTEGER:
      return VAR_INT;
   case SCIP_VARTYPE_IMPLINT:
      return VAR_IMP;
   case SCIP_VARTYPE_CONTINUOUS:
      return VAR_CON;
   default:
      SCIPerrorMessage("invalid SCIP variable type <%d> in ZIMPL callback xlp_getclass()\n", SCIPvarGetType(scipvar));
      readerdata->readerror = TRUE;
      break;
   }

   return VAR_CON;
}

/** returns lower bound */
Bound* xlp_getlower(
   const Lps*            data,               /**< pointer to reader data */
   const Var*            var                 /**< variable */
   )
{
   SCIP* scip;
   SCIP_READERDATA* readerdata;
   SCIP_VAR* scipvar;
   SCIP_Real lb;
   char s[SCIP_MAXSTRLEN];
   BoundType boundtype;
   Numb* numb;
   Bound* bound;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

#ifdef WITH_EXACTSOLVE
   SCIPerrorMessage("xlp_getlower: exact version not supported.\n");
   readerdata->readerror = TRUE;
   return NULL;
#endif

   scip = readerdata->scip;
   assert(scip != NULL);

   scipvar = (SCIP_VAR*)var;
   assert(scipvar != NULL);

   /* collect lower bound */
   lb = SCIPvarGetLbGlobal(scipvar);
   numb = NULL;

   /* check if lower bound is infinity */
   if( SCIPisInfinity(scip, -lb) )
      boundtype = BOUND_MINUS_INFTY;
   else if( SCIPisInfinity(scip, lb) )
      boundtype = BOUND_INFTY;
   else
   {
      boundtype = BOUND_VALUE;

      /* create double form string */
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "%.20f", lb);
      numb = numb_new_ascii(s);
   }

   /* create bound */
   bound = bound_new(boundtype, numb);

   if( numb != NULL )
      numb_free(numb);

   return bound;
}

/** returns upper bound */
Bound* xlp_getupper(
   const Lps*            data,               /**< pointer to reader data */
   const Var*            var                 /**< variable */
   )
{
   SCIP* scip;
   SCIP_READERDATA* readerdata;
   SCIP_VAR* scipvar;
   SCIP_Real ub;
   char s[SCIP_MAXSTRLEN];
   BoundType boundtype;
   Numb* numb;
   Bound* bound;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

#ifdef WITH_EXACTSOLVE
   SCIPerrorMessage("xlp_getupper: exact version not supported.\n");
   readerdata->readerror = TRUE;
   return NULL;
#endif

   scip = readerdata->scip;
   assert(scip != NULL);

   scipvar = (SCIP_VAR*)var;
   assert(scipvar != NULL);

   /* collect upper bound */
   ub = SCIPvarGetUbGlobal(scipvar);
   numb = NULL;

   /* check if upper bound is infinity */
   if( SCIPisInfinity(scip, -ub) )
      boundtype = BOUND_MINUS_INFTY;
   else if( SCIPisInfinity(scip, ub) )
      boundtype = BOUND_INFTY;
   else
   {
      boundtype = BOUND_VALUE;
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "%.20f", ub);
      numb = numb_new_ascii(s);
   }

   /* create ZIMPL bound */
   bound = bound_new(boundtype, numb);

   if (numb != NULL)
      numb_free(numb);

   return bound;
}

/* set the name of the objective function */
void xlp_objname(
   Lps*                  data,               /**< pointer to reader data */
   const char*           name                /**< name of the objective function */
   )
{  /*lint --e{715}*/
   /* nothing to be done */
}

/* set the name of the objective function */
void xlp_setdir(
   Lps*                  data,               /**< pointer to reader data */
   Bool                  minimize            /**<True if the problem should be minimized, False if it should be maximized  */
   )
{
   SCIP* scip;
   SCIP_READERDATA* readerdata;
   SCIP_OBJSENSE objsense;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

   scip = readerdata->scip;
   assert(scip != NULL);

   objsense = (minimize ? SCIP_OBJSENSE_MINIMIZE : SCIP_OBJSENSE_MAXIMIZE);
   SCIP_CALL_ABORT( SCIPsetObjsense(scip, objsense) );

#ifdef WITH_EXACTSOLVE
   readerdata->objsense = objsense;
#endif
}

/** changes objective coefficient of a variable */
void xlp_addtocost(
   Lps*                  data,               /**< pointer to reader data */
   Var*                  var,                /**< variable */
   const Numb*           cost                /**< objective coefficient */
   )
{
   SCIP* scip;
   SCIP_READERDATA* readerdata;
   SCIP_VAR* scipvar;
   SCIP_Real scipval;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

   scip = readerdata->scip;
   assert(scip != NULL);

   scipvar = (SCIP_VAR*)var;

#ifdef WITH_EXACTSOLVE
   {
      mpq_t scipvalmpq;
      int varidx;

      varidx = SCIPvarGetIndex(scipvar);
      assert(SCIPvarGetIndex(readerdata->vars[varidx]) == varidx);
      assert(readerdata->vars[varidx] == scipvar);

      mpq_init(scipvalmpq);

      /* get exact objective coefficient used in exactlp constraint handler and approximate FP-value used in FP-problem */
      numb_get_mpq(cost, scipvalmpq);
      scipval = mpqGetRealApprox(scip, scipvalmpq);

#ifdef CREATEEXACTLPCONS_OUT
      SCIPdebugMessage("zimpl reader: change obj<%g> of var: add<%g> as approx", SCIPvarGetObj(scipvar), scipval);
      gmp_printf(" (<%Qd> as exact)", scipvalmpq);
#endif

      /* check if we need to scale the objective function to FP-representable values. if we want to work on an
       * FP-relaxation we have to ensure that the objective coefficient stored in scip for the variable is identical to
       * the one in the exactlp constraint handler, otherwise the dual bounding methods are not reliable
       */
      if( SCIPuseFPRelaxation(scip) && !readerdata->objneedscaling )
      {
         mpq_t tmp;

         mpq_init(tmp);
         mpq_set_d(tmp, scipval);

         if( !mpq_equal(scipvalmpq, tmp) )
            readerdata->objneedscaling = TRUE;

         mpq_clear(tmp);
      }

      /* change objective coefficient of variable in exactlp constraint handler and in FP problem */
      mpq_add(readerdata->obj[varidx], readerdata->obj[varidx], scipvalmpq);
      SCIP_CALL_ABORT( SCIPchgVarObj(scip, scipvar, SCIPvarGetObj(scipvar) + scipval) );

#ifdef CREATEEXACTLPCONS_OUT
      SCIPprintVar(scip_, scipvar, NULL);
#endif

      mpq_clear(scipvalmpq);
   }
#else
   /* get objective coefficient */
   scipval = numb_todbl(cost);

   /* change objective coefficient of variable */
   SCIP_CALL_ABORT( SCIPchgVarObj(scip, scipvar, SCIPvarGetObj(scipvar) + scipval) );
#endif
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyZpl)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderZpl(scip) );
 
   return SCIP_OKAY;
}


/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeZpl NULL


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadZpl)
{  /*lint --e{715}*/
   SCIP_READERDATA* readerdata;
   SCIP_RETCODE retcode;
   char oldpath[SCIP_MAXSTRLEN];
   char buffer[SCIP_MAXSTRLEN];
   char compextension[SCIP_MAXSTRLEN];
   char namewithoutpath[SCIP_MAXSTRLEN];
   char* path;
   char* name;
   char* extension;
   char* compression;
   char* paramstr;

   SCIP_Bool changedir;
   int i;

#ifdef WITH_EXACTSOLVE
   SCIP_CONS* cons;
   SCIP_VAR* vartmp;
   SCIP_Bool initial;
   SCIP_Bool separate;
   SCIP_Bool enforce;
   SCIP_Bool check;
   SCIP_Bool propagate;
   SCIP_Bool local;
   SCIP_Bool modifiable;
   SCIP_Bool dynamic;
   SCIP_Bool removable;
   SCIP_Bool stickingatnode;
   char consname[SCIP_MAXSTRLEN];
   mpq_t lbtmp;
   mpq_t ubtmp;
   mpq_t objtmp;
   int idx;
   int nswitch;
   int c;
#endif

   SCIP_CALL( SCIPgetBoolParam(scip, "reading/zplreader/changedir", &changedir) );

   SCIP_CALL( SCIPallocBuffer(scip, &readerdata) );

   readerdata->scip = scip;
   readerdata->sol = NULL;
   readerdata->valid = FALSE;
   readerdata->branchpriowarning = FALSE;
   readerdata->readerror = FALSE;

#ifdef WITH_EXACTSOLVE
   /* initialize exactlp relevant values of reader data */
   readerdata->conshdlrdata = SCIPconshdlrGetData(SCIPfindConshdlr(scip, "exactlp"));
   assert(readerdata->conshdlrdata != NULL);
   readerdata->objsense = SCIP_OBJSENSE_MINIMIZE;
   readerdata->nvars = 0;
   readerdata->ninfbounds = 0;
   readerdata->ninfintbounds = 0;
   readerdata->nlargebounds = 0;
   readerdata->nconss = 0;
   readerdata->consssize = 0;
   readerdata->nsplitconss = 0;
   readerdata->nnonz = 0;
   readerdata->nonzsize = 0;
   readerdata->nintegral = 0;
   mpq_init(readerdata->minabsval);
   mpq_init(readerdata->maxabsval);
   mpq_set(readerdata->minabsval, *posInfinity(readerdata->conshdlrdata));
   mpq_set(readerdata->maxabsval, *negInfinity(readerdata->conshdlrdata));
   readerdata->objneedscaling = FALSE;

   /* allocate and initialize variable specific information */
   readerdata->varssize = 1024;
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &readerdata->vars, readerdata->varssize) );
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &readerdata->lb, readerdata->varssize) );
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &readerdata->ub, readerdata->varssize) );
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &readerdata->obj, readerdata->varssize) );
   for( i = 0; i < readerdata->varssize; ++i )
   {
      mpq_init(readerdata->lb[i]);
      mpq_init(readerdata->ub[i]);
      mpq_init(readerdata->obj[i]);
      mpq_set_si(readerdata->obj[i], 0, 1);
   }

   /* allocate and initialize constraint specific information */
   readerdata->consssize = 1024;
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &readerdata->beg, readerdata->consssize) );
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &readerdata->len, readerdata->consssize) );
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &readerdata->lhs, readerdata->consssize) );
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &readerdata->rhs, readerdata->consssize) );
   for( i = 0; i < readerdata->consssize; ++i )
   {
      mpq_init(readerdata->lhs[i]);
      mpq_init(readerdata->rhs[i]);
   }

   /* set initialize start index and length of first constraint */
   readerdata->beg[readerdata->nconss] = 0;
   readerdata->len[readerdata->nconss] = 0;

   /* allocate and initialize matrix specific information */
   readerdata->nonzsize = 1024;
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &readerdata->val, readerdata->nonzsize) );
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &readerdata->ind, readerdata->nonzsize) );
   for( i = 0; i < readerdata->nonzsize; ++i )
   {
      mpq_init(readerdata->val[i]);
   }
#endif

   path = NULL;
   oldpath[0] = '\0';
   if( changedir )
   {
      /* change to the directory of the ZIMPL file, s.t. paths of data files read by the ZIMPL model are relative to
       * the location of the ZIMPL file
       */
      (void)strncpy(buffer, filename, SCIP_MAXSTRLEN-1);
      buffer[SCIP_MAXSTRLEN-1] = '\0';
      SCIPsplitFilename(buffer, &path, &name, &extension, &compression);
      if( compression != NULL )
         (void) SCIPsnprintf(compextension, SCIP_MAXSTRLEN, ".%s", compression);
      else
         *compextension = '\0';
      (void) SCIPsnprintf(namewithoutpath, SCIP_MAXSTRLEN, "%s.%s%s", name, extension, compextension);
      if( getcwd(oldpath, SCIP_MAXSTRLEN) == NULL )
      {
         SCIPerrorMessage("error getting the current path\n");
         return SCIP_READERROR;
      }
      if( path != NULL )
      {
         if( chdir(path) != 0 )
         {
            SCIPerrorMessage("error changing to directory <%s>\n", path);
            return SCIP_NOFILE;
         }
      }
      filename = namewithoutpath;
   }

   /* get current path for output */
   if( SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_NORMAL )
   {
      char currentpath[SCIP_MAXSTRLEN];
      if( getcwd(currentpath, SCIP_MAXSTRLEN) == NULL )
      {
         SCIPerrorMessage("error getting the current path\n");
         return SCIP_READERROR;
      }
      /* an extra blank line should be printed separately since the buffer message handler only handle up to one line
       *  correctly */
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n");
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "base directory for ZIMPL parsing: <%s>\n", currentpath);
      /* an extra blank line should be printed separately since the buffer message handler only handle up to one line
       *  correctly */
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n");
   }

   /* get the parameter string */
   SCIP_CALL( SCIPgetStringParam(scip, "reading/zplreader/parameters", &paramstr) );
   if( strcmp(paramstr, "-") == 0 )
   {
      /* call ZIMPL parser without arguments */
#ifdef WITH_REDUCEDSOLVE
      if( !zpl_read(filename, FALSE, (void*)readerdata) )
         readerdata->readerror = TRUE;
#else
      if( !zpl_read(filename, TRUE, (void*)readerdata) )
         readerdata->readerror = TRUE;
#endif
   }
   else
   {
      char dummy[2] = "x";
      char** argv;
      int argc;
      int p;
      int len;

      len = (int) strlen(paramstr);
      SCIP_CALL( SCIPallocBufferArray(scip, &argv, len+1) );
      argv[0] = dummy; /* argument 0 is irrelevant */
      argc = 1;
      p = 0;
      while( p < len )
      {
         int arglen;

         /* process next argument */
         SCIP_CALL( SCIPallocBufferArray(scip, &argv[argc], len+1) );  /*lint !e866*/
         arglen = 0;

         /* skip spaces */
         while( p < len && paramstr[p] == ' ' )
            p++;

         /* process characters */
         while( p < len && paramstr[p] != ' ' )
         {
            switch( paramstr[p] )
            {
            case '"':
               p++;
               /* read characters as they are until the next " */
               while( p < len && paramstr[p] != '"' )
               {
                  argv[argc][arglen] = paramstr[p];
                  arglen++;
                  p++;
               }
               p++; /* skip final " */
               break;
            case '\\':
               /* read next character as it is */
               p++;
               argv[argc][arglen] = paramstr[p];
               arglen++;
               p++;
               break;
            default:
               argv[argc][arglen] = paramstr[p];
               arglen++;
               p++;
               break;
            }
         }
         argv[argc][arglen] = '\0';

         /* check for empty argument */
         if( arglen == 0 )
         {
            SCIPfreeBufferArray(scip, &argv[argc]);
         }
         else
            argc++;
      }

      /* append file name as last argument */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &argv[argc], filename, (int) strlen(filename)+1) );  /*lint !e866*/
      argc++;

      /* display parsed arguments */
      if( SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_FULL )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "ZIMPL arguments:\n");
         for( i = 1; i < argc; ++i )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "%d: <%s>\n", i, argv[i]);
         }
      }

      /* call ZIMPL parser with arguments */
#ifdef WITH_REDUCEDSOLVE
      if( !zpl_read_with_args(argv, argc, FALSE, (void*)readerdata) )
         readerdata->readerror = TRUE;
#else
      if( !zpl_read_with_args(argv, argc, TRUE, (void*)readerdata) )
         readerdata->readerror = TRUE;
#endif
      /* free argument memory */
      for( i = argc - 1; i >= 1; --i )
      {
         SCIPfreeBufferArray(scip, &argv[i]);
      }
      SCIPfreeBufferArray(scip, &argv);
   }

   if( changedir )
   {
      /* change directory back to old path */
      if( path != NULL )
      {
         if( chdir(oldpath) != 0 )
         {
            SCIPwarningMessage(scip, "error changing back to directory <%s>\n", oldpath);
         }
      }
   }

#ifdef WITH_EXACTSOLVE
   /* create exactlp constraint from variable, constraint and matrix information stored in reader data */
#ifdef CREATEEXACTLPCONS_OUT
   SCIPdebugMessage("zimpl reader: after adding all vars sort them by probindex instead of index\n");
#endif

   /* sort variable specific information by variable's probindex instead of index */
   for( c = 0; c < readerdata->nconss; ++c )
   {
      for( i = readerdata->beg[c]; i < readerdata->beg[c] + readerdata->len[c]; ++i )
      {
         assert(readerdata->ind[i] == SCIPvarGetIndex(readerdata->vars[readerdata->ind[i]]));
         readerdata->ind[i] = SCIPvarGetProbindex(readerdata->vars[readerdata->ind[i]]);
      }
   }

   mpq_init(lbtmp);
   mpq_init(ubtmp);
   mpq_init(objtmp);
   idx = 0;
   nswitch = 0;
   for( i = 0; i < readerdata->nvars-1+nswitch; ++i )
   {
      int probindex;

      assert(idx <= i);
      assert(idx < readerdata->nvars);
      assert(nswitch <= readerdata->nvars);

      probindex = SCIPvarGetProbindex(readerdata->vars[idx]);
      assert(probindex >= idx);

#ifdef CREATEEXACTLPCONS_OUT
      SCIPdebugMessage("i=%d:  on idxpos<%d> var<%s> with probidx<%d> | on probidxpos<%d> var<%s> with probidx<%d>\n", i,
         idx, SCIPvarGetName(readerdata->vars[idx]), SCIPvarGetProbindex(readerdata->vars[idx]),
         probindex, SCIPvarGetName(readerdata->vars[probindex]), SCIPvarGetProbindex(readerdata->vars[probindex]) );
#endif

      /* position of current variable is incorrect */
      if( idx != probindex )
      {
         /* swap current variable with the one that is located at the current variable's correct position */
         mpq_set(lbtmp, readerdata->lb[probindex]);
         mpq_set(ubtmp, readerdata->ub[probindex]);
         mpq_set(objtmp, readerdata->obj[probindex]);
         vartmp = readerdata->vars[probindex];

         mpq_set(readerdata->lb[probindex], readerdata->lb[idx]);
         mpq_set(readerdata->ub[probindex], readerdata->ub[idx]);
         mpq_set(readerdata->obj[probindex], readerdata->obj[idx]);
         readerdata->vars[probindex] = readerdata->vars[idx];

         mpq_set(readerdata->lb[idx], lbtmp);
         mpq_set(readerdata->ub[idx], ubtmp);
         mpq_set(readerdata->obj[idx], objtmp);
         readerdata->vars[idx] = vartmp;

         nswitch++;
#ifdef CREATEEXACTLPCONS_OUT
         SCIPdebugMessage("   swap: idxpos<%d> var<%s> with probidx<%d> | on probidxpos<%d> var<%s> with probidx<%d>\n",
            idx, SCIPvarGetName(readerdata->vars[idx]), SCIPvarGetProbindex(readerdata->vars[idx]),
            probindex, SCIPvarGetName(readerdata->vars[probindex]), SCIPvarGetProbindex(readerdata->vars[probindex]) );
#endif
      }
      /* position of current variable is correct */
      else
      {
#ifdef CREATEEXACTLPCONS_OUT
         SCIPdebugMessage("NO swap: idxpos<%d> var<%s> with probidx<%d> | on probidxpos<%d> var<%s> with probidx<%d>\n",
            idx, SCIPvarGetName(readerdata->vars[idx]), SCIPvarGetProbindex(readerdata->vars[idx]),
            probindex, SCIPvarGetName(readerdata->vars[probindex]), SCIPvarGetProbindex(readerdata->vars[probindex]) );
#endif

         /* move to variable on next position */
         idx++;
      }
   }

   mpq_clear(objtmp);
   mpq_clear(ubtmp);
   mpq_clear(lbtmp);

#ifndef NDEBUG
   for( i = 0; i < readerdata->nvars; ++i )
   {
      assert(SCIPvarGetProbindex(readerdata->vars[i]) == i);

      assert(SCIPisInfinity(scip, SCIPvarGetLbLocal(readerdata->vars[i]))
         || SCIPisInfinity(scip, -SCIPvarGetLbLocal(readerdata->vars[i]))
         || mpqGetRealRelax(scip, readerdata->lb[i], GMP_RNDD) == SCIPvarGetLbLocal(readerdata->vars[i]));

      assert(SCIPisInfinity(scip, SCIPvarGetUbLocal(readerdata->vars[i]))
         || SCIPisInfinity(scip, -SCIPvarGetUbLocal(readerdata->vars[i]))
         || mpqGetRealRelax(scip, readerdata->ub[i], GMP_RNDU) == SCIPvarGetUbLocal(readerdata->vars[i]));

      assert(SCIPisInfinity(scip, SCIPvarGetObj(readerdata->vars[i]))
         || SCIPisInfinity(scip, -SCIPvarGetObj(readerdata->vars[i]))
         || mpqGetRealApprox(scip, readerdata->obj[i]) == SCIPvarGetObj(readerdata->vars[i]));
   }
#endif

   /* create exactlp constraint */
   initial = TRUE;
   separate = TRUE;
   enforce = TRUE;
   check = enforce;
   propagate = TRUE;
   local = FALSE;
   modifiable = FALSE;
   dynamic = FALSE;
   removable = dynamic;
   stickingatnode = FALSE;
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "exactlp");

   SCIP_CALL_ABORT( SCIPcreateConsExactlp(scip, &cons, consname, readerdata->objsense, readerdata->nvars,
         readerdata->ninfbounds, readerdata->ninfintbounds, readerdata->nlargebounds, readerdata->obj,readerdata->lb,
         readerdata->ub, readerdata->nconss, readerdata->nsplitconss, readerdata->lhs, readerdata->rhs, readerdata->nnonz,
         readerdata->nintegral, readerdata->beg, readerdata->len, readerdata->ind, readerdata->val, readerdata->minabsval,
         readerdata->maxabsval, readerdata->objneedscaling, initial, separate, enforce, check, propagate, local,
         modifiable, dynamic, removable, stickingatnode) );
   SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );
   if( cons != NULL )
   {
      SCIP_CALL_ABORT( SCIPreleaseCons(scip, &cons) );
   }

#ifdef CREATEEXACTLPCONS_OUT
   SCIPdebugMessage("zimpl reader: released exactlp constraint\n");
#endif
#endif /* end of WITH_EXACTSOLVE */

   if( readerdata->valid )
   {
      SCIP_Bool stored;

      assert(readerdata->sol != NULL);

      stored = FALSE;

#ifdef WITH_EXACTSOLVE
      /* current exact version of SCIP only supports primal solutions found within the exactlp constraint handler */
      SCIP_CALL( SCIPfreeSol(scip, &readerdata->sol) );
#else
      /* add primal solution to solution candidate storage, frees the solution afterwards */
      SCIP_CALL( SCIPaddSolFree(scip, &readerdata->sol, &stored) );
#endif
      if( stored )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "ZIMPL starting solution candidate accepted\n");
      }
   }

   *result = SCIP_SUCCESS;

#ifdef WITH_EXACTSOLVE
   /* free matrix specific information */
   assert(readerdata->nnonz <= readerdata->nonzsize);
   for( i = 0; i < readerdata->nonzsize; ++i )
   {
      mpq_clear(readerdata->val[i]);
   }
   mpq_clear(readerdata->minabsval);
   mpq_clear(readerdata->maxabsval);
   SCIPfreeMemoryArray(scip, &readerdata->ind);
   SCIPfreeMemoryArray(scip, &readerdata->val);

   /* free constraint specific information */
   assert(readerdata->nconss <= readerdata->consssize);
   for( i = 0; i < readerdata->consssize; ++i )
   {
      mpq_clear(readerdata->rhs[i]);
      mpq_clear(readerdata->lhs[i]);
   }
   SCIPfreeMemoryArray(scip, &readerdata->rhs);
   SCIPfreeMemoryArray(scip, &readerdata->lhs);
   SCIPfreeMemoryArray(scip, &readerdata->len);
   SCIPfreeMemoryArray(scip, &readerdata->beg);

   /* free variable specific information */
   assert(readerdata->nvars <= readerdata->varssize);
   for( i = 0; i < readerdata->varssize; ++i )
   {
      mpq_clear(readerdata->obj[i]);
      mpq_clear(readerdata->ub[i]);
      mpq_clear(readerdata->lb[i]);
   }
   SCIPfreeMemoryArray(scip, &readerdata->obj);
   SCIPfreeMemoryArray(scip, &readerdata->ub);
   SCIPfreeMemoryArray(scip, &readerdata->lb);
   SCIPfreeMemoryArray(scip, &readerdata->vars);
#endif

   /* evaluate if a reading error occurred */
   if( readerdata->readerror )
      retcode = SCIP_READERROR;
   else
      retcode = SCIP_OKAY;

   /* free primal solution candidate */
   if( readerdata->sol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &readerdata->sol) );
   }

   /* free reader data */
   SCIPfreeBuffer(scip, &readerdata);

   return retcode;
}


/** problem writing method of reader */
#define readerWriteZpl NULL

#endif
#endif



/*
 * reader specific interface methods
 */

/** includes the zpl file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderZpl(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
#ifdef WITH_ZIMPL
#if (ZIMPL_VERSION >= 320)
   SCIP_READERDATA* readerdata;
   char extcodename[SCIP_MAXSTRLEN];

   /* create zpl reader data */
   readerdata = NULL;

   /* include zpl reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerCopyZpl, readerFreeZpl, readerReadZpl, readerWriteZpl, readerdata) );

   /* add zpl reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/zplreader/dynamiccols", "should columns be added and removed dynamically to the LP?",
         NULL, FALSE, FALSE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/zplreader/changedir", "should the current directory be changed to that of the ZIMPL file before parsing?",
         NULL, FALSE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/zplreader/usestartsol", "should ZIMPL starting solutions be forwarded to SCIP?",
         NULL, FALSE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddStringParam(scip,
         "reading/zplreader/parameters", "additional parameter string passed to the ZIMPL parser (or - for no additional parameters)",
         NULL, FALSE, "-", NULL, NULL) );

   (void) SCIPsnprintf(extcodename, SCIP_MAXSTRLEN, "ZIMPL %d.%d.%d", 
      ZIMPL_VERSION/100, (ZIMPL_VERSION%100)/10, ZIMPL_VERSION%10);
   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, extcodename, "Zuse Institute Mathematical Programming Language developed by T. Koch (zimpl.zib.de)"));
#else
   SCIPwarningMessage(scip, "SCIP does only support ZIMPL 3.2.0 and higher. Please update your ZIMPL version %d.%d.%d\n",
      ZIMPL_VERSION/100, (ZIMPL_VERSION%100)/10, ZIMPL_VERSION%10);
#endif
#endif

   return SCIP_OKAY;
}

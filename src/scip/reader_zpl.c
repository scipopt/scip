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

/**@file   reader_zpl.c
 * @ingroup FILEREADERS 
 * @brief  ZIMPL model file reader
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define CREATEEXACTLPCONS_OUT /** uncomment to get more debug messages (SCIP_DEBUG) about creation of exactlp constraints; 
//                               *  also in cons_exactlp.c */

#include "scip/reader_zpl.h"

#ifdef WITH_ZIMPL

#include <unistd.h>
#include <string.h>
#include <assert.h>

#ifdef EXACTSOLVE
#include <gmp.h> /* gmp.h has to be included before zimpl/mme.h */
#include "scip/cons_exactlp.h"
#include "scip/intervalarith.h"
#endif

#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "scip/cons_indicator.h"
#include "scip/cons_quadratic.h"
#include "scip/pub_misc.h"

/* include the ZIMPL headers necessary to define the LP and MINLP construction interface */
#include "zimpl/bool.h"
#include "zimpl/ratlptypes.h"
#include "zimpl/mme.h"
#include "zimpl/xlpglue.h"
#include "zimpl/zimpllib.h"

/** ZIMPL_VERSION is defined by ZIMPL version 3.00 and higher. ZIMPL 3.00 made same changes in the interface to SCIP. */
#if (ZIMPL_VERSION >= 300)
#include "zimpl/mono.h"
#endif

#define READER_NAME             "zplreader"
#define READER_DESC             "file reader for ZIMPL model files"
#define READER_EXTENSION        "zpl"

/*
 * LP construction interface of ZIMPL
 */

/* ZIMPL does not support user data in callbacks - we have to use a static variables */
static SCIP* scip_ = NULL;
static SCIP_Real* startvals_ = NULL;
static SCIP_VAR** startvars_ = NULL;
static int startvalssize_ = 0;
static int nstartvals_ = 0;
static SCIP_Bool issuedbranchpriowarning_ = FALSE;
static SCIP_Bool readerror_ = FALSE;

#ifdef EXACTSOLVE
static SCIP_CONSHDLRDATA* conshdlrdata_ = NULL;
static SCIP_VAR** vars_ = NULL;  
static SCIP_OBJSENSE objsense_;
static mpq_t* lb_ = NULL;
static mpq_t* ub_ = NULL;
static mpq_t* obj_ = NULL;
static int nvars_ = 0;
static int ninfbounds_ = 0;
static int ninfintbounds_ = 0;
static int nlargebounds_ = 0;
static int varssize_ = 0;

static int* beg_ = NULL;
static int* len_ = NULL;
static mpq_t* lhs_ = NULL;
static mpq_t* rhs_ = NULL;
static int nconss_ = 0;
static int nsplitconss_ = 0;
static int consssize_ = 0;

static mpq_t* val_ = NULL;
static int* ind_ = NULL;
static int nnonz_ = 0;
static int nintegral_ = 0;
static int nonzsize_ = 0;
static mpq_t minabsval_;
static mpq_t maxabsval_;
static SCIP_Bool objneedscaling_ = FALSE;
#define LARGEBOUND              1e+06
#endif


void xlp_alloc(const char* name, Bool need_startval)
{  /*lint --e{715}*/
   /* create problem */
   SCIP_CALL_ABORT( SCIPcreateProb(scip_, name, NULL, NULL, NULL, NULL, NULL, NULL) );
}

void xlp_free(void)
{   
   /* nothing to be done here */
}

void xlp_stat(void)
{
   /* nothing to be done here */
}

void xlp_scale(void)
{
   /* nothing to be done here */
}

void xlp_write(FILE* fp, LpFormat format, const char* title)
{  /*lint --e{715}*/
   /* nothing to be done here */
}

void xlp_transtable(FILE* fp, LpFormat format)
{  /*lint --e{715}*/
   /* nothing to be done here */
}

void xlp_orderfile(FILE* fp, LpFormat format)
{  /*lint --e{715}*/
   /* nothing to be done here */
}

void xlp_mstfile(FILE* fp, LpFormat format)
{  /*lint --e{715}*/
   /* nothing to be done here */
}

Bool xlp_conname_exists(const char* conname)
{  /*lint --e{715}*/
   return (SCIPfindCons(scip_, conname) != NULL);
}

/** ZIMPL_VERSION is defined by ZIMPL version 3.00 and higher. ZIMPL 3.00 made same changes in the interface to SCIP. */
#if (ZIMPL_VERSION >= 300)

/** method creates a linear constraint and is called directly from ZIMPL 
 *
 *  @note this method is used by ZIMPL from version 3.00; 
 */
#ifndef EXACTSOLVE
Bool xlp_addcon_term(
   const char*           name,               /**< constraint name */
   ConType               type,               /**< constraint type */
   const Numb*           lhs,                /**< left hand side */
   const Numb*           rhs,                /**< right hand side */
   unsigned int          flags,              /**< special constraint flags */
   const Term*           term      
   )
{
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
   int  i;
   int  maxdegree;
   
   switch( type )
   {
   case CON_FREE:
      sciplhs = -SCIPinfinity(scip_);
      sciprhs = SCIPinfinity(scip_);
      break;
   case CON_LHS:
      sciplhs = (SCIP_Real)numb_todbl(lhs);
      sciprhs = SCIPinfinity(scip_);
      break;
   case CON_RHS:
      sciplhs = -SCIPinfinity(scip_);
      sciprhs = (SCIP_Real)numb_todbl(rhs);
      break;
   case CON_RANGE:
      sciplhs = (SCIP_Real)numb_todbl(lhs);
      sciprhs = (SCIP_Real)numb_todbl(rhs);
      break;
   case CON_EQUAL:
      sciplhs = (SCIP_Real)numb_todbl(lhs);
      sciprhs = (SCIP_Real)numb_todbl(rhs);
      assert(sciplhs == sciprhs);
      break;
   default:
      SCIPwarningMessage("invalid constraint type <%d> in ZIMPL callback xlp_addcon()\n", type);
      sciplhs = (SCIP_Real)numb_todbl(lhs);
      sciprhs = (SCIP_Real)numb_todbl(rhs);
      readerror_ = TRUE;
      break;
   }

   initial = !((flags & LP_FLAG_CON_SEPAR) != 0);
   separate = TRUE;
   enforce = TRUE;
   check = enforce;
   propagate = TRUE;
   local = FALSE;
   modifiable = FALSE;
   dynamic = ((flags & LP_FLAG_CON_SEPAR) != 0);
   removable = dynamic;
   maxdegree = term_get_degree(term);

   if (maxdegree <= 1)
   {
      /* if the constraint gives an indicator constraint */
      if ( flags & LP_FLAG_CON_INDIC )
      {
         Bool lhsIndCons = FALSE;  /* generate lhs form for indicator constraints */
         Bool rhsIndCons = FALSE;  /* generate rhs form for indicator constraints */

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
         default:
            SCIPwarningMessage("invalid constraint type <%d> in ZIMPL callback xlp_addcon()\n", type);
            readerror_ = TRUE;
            break;
         }

         /* insert lhs form of indicator */
         if ( lhsIndCons )
         {
            SCIP_CALL_ABORT( SCIPcreateConsIndicator(scip_, &cons, name, NULL, 0, NULL, NULL, -sciplhs,
                  initial, separate, enforce, check, propagate, local, dynamic, removable, FALSE) );
            SCIP_CALL_ABORT( SCIPaddCons(scip_, cons) );

            for (i = 0; i < term_get_elements(term); i++)
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
                  SCIP_CALL( SCIPsetBinaryVarIndicator(scip_, cons, scipvar) );
               }
               else
               {
                  assert(!numb_equal(mono_get_coeff(mono), numb_zero()));
                  assert(mono_is_linear(mono));

                  scipval = -numb_todbl(mono_get_coeff(mono));
                  SCIP_CALL_ABORT( SCIPaddVarIndicator(scip_, cons, scipvar, scipval) );
               }
            }
         }

         /* insert rhs form of indicator */
         if ( rhsIndCons )
         {
            SCIP_CALL_ABORT( SCIPcreateConsIndicator(scip_, &cons, name, NULL, 0, NULL, NULL, sciprhs,
                  initial, separate, enforce, check, propagate, local, dynamic, removable, FALSE) );
            SCIP_CALL_ABORT( SCIPaddCons(scip_, cons) );

            for (i = 0; i < term_get_elements(term); i++)
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
                  SCIP_CALL( SCIPsetBinaryVarIndicator(scip_, cons, scipvar) );
               }
               else
               {
                  assert(!numb_equal(mono_get_coeff(mono), numb_zero()));
                  assert(mono_is_linear(mono));

                  scipval = numb_todbl(mono_get_coeff(mono));
                  SCIP_CALL_ABORT( SCIPaddVarIndicator(scip_, cons, scipvar, scipval) );
               }
            }
         }
      }
      else
      {
         SCIP_CALL_ABORT( SCIPcreateConsLinear(scip_, &cons, name, 0, NULL, NULL, sciplhs, sciprhs,
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
         SCIP_CALL_ABORT( SCIPaddCons(scip_, cons) );

         for (i = 0; i < term_get_elements(term); i++)
         {
            SCIP_VAR* scipvar;
            SCIP_Real scipval;

            assert(!numb_equal(mono_get_coeff(term_get_element(term, i)), numb_zero()));
            assert(mono_is_linear(term_get_element(term, i)));

            scipvar = (SCIP_VAR*)mono_get_var(term_get_element(term, i), 0);
            scipval = numb_todbl(mono_get_coeff(term_get_element(term, i)));

            SCIP_CALL_ABORT( SCIPaddCoefLinear(scip_, cons, scipvar, scipval) );
         }
      }
   }
   else if (maxdegree == 2)
   {
      int        n_linvar   = 0;
      int        n_quadterm = 0;
      SCIP_VAR** linvar;
      SCIP_VAR** quadvar1;
      SCIP_VAR** quadvar2;
      SCIP_Real* lincoeff;
      SCIP_Real* quadcoeff;
      Mono*      monom;
  	  
      SCIP_CALL( SCIPallocBufferArray(scip_, &linvar,    term_get_elements(term)) );
      SCIP_CALL( SCIPallocBufferArray(scip_, &quadvar1,  term_get_elements(term)) );
      SCIP_CALL( SCIPallocBufferArray(scip_, &quadvar2,  term_get_elements(term)) );
      SCIP_CALL( SCIPallocBufferArray(scip_, &lincoeff,  term_get_elements(term)) );
      SCIP_CALL( SCIPallocBufferArray(scip_, &quadcoeff, term_get_elements(term)) );
  	  
      for (i = 0; i < term_get_elements(term); ++i)
      {
         monom = term_get_element(term, i);
         assert(!numb_equal(mono_get_coeff(monom), numb_zero()));
         assert(mono_get_degree(monom) <= 2);
         assert(mono_get_degree(monom) > 0);
         if (mono_get_degree(monom) == 1)
         {
            linvar  [n_linvar] = (SCIP_VAR*)mono_get_var(monom, 0);
            lincoeff[n_linvar] = numb_todbl(mono_get_coeff(monom));
            ++n_linvar;
         }
         else
         {
            assert(mono_get_degree(monom) == 2);
            quadvar1 [n_quadterm] = (SCIP_VAR*)mono_get_var(monom, 0);
            quadvar2 [n_quadterm] = (SCIP_VAR*)mono_get_var(monom, 1);
            quadcoeff[n_quadterm] = numb_todbl(mono_get_coeff(monom));
            ++n_quadterm;
         }
      }

      SCIP_CALL_ABORT( SCIPcreateConsQuadratic(scip_, &cons, name, n_linvar, linvar, lincoeff, n_quadterm, quadvar1, quadvar2, quadcoeff, sciplhs, sciprhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );
      SCIP_CALL_ABORT( SCIPaddCons(scip_, cons) );
      
      SCIPfreeBufferArray(scip_, &linvar);
      SCIPfreeBufferArray(scip_, &quadvar1);
      SCIPfreeBufferArray(scip_, &quadvar2);
      SCIPfreeBufferArray(scip_, &lincoeff);
      SCIPfreeBufferArray(scip_, &quadcoeff);
   }
   else
   {
      SCIPerrorMessage("xpl_addcon_term for degree > 2 not implemented\n");
      return TRUE;
   }

   SCIP_CALL_ABORT( SCIPreleaseCons(scip_, &cons) );
   
   return FALSE;
}
#else
/** method stores a linear constraint in the constraints matrix and checks whether it has to be split into lhs and rhs 
 *  constraints for building an FP relaxation
 */
static
Bool storeConstraint(
   ConType               type,               /**< constraint type */
   const Term*           term
   )
{
   Bool split;
   mpq_t absval;
   int i;

   assert(nconss_ >= 0);
   assert(nnonz_ >= 0);
   assert(beg_[nconss_] == nnonz_);
   assert(len_[nconss_] == 0);
  
#ifdef CREATEEXACTLPCONS_OUT
   SCIPdebugMessage("zimpl reader: store linear cons with <%d> vars in matrix, nconss<%d>, next matrixbeg<%d>, next matrixlen<%d>\n", 
      term_get_elements(term), nconss_, beg_[nconss_], len_[nconss_]);
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
      numb_get_mpq(mono_get_coeff(term_get_element(term, i)), val_[nnonz_]);

      /* when we use an FP relaxation for calculating dual bounds, we have to split the row into lhs and rhs part if
       *  - lhs and rhs of constraints are not -inf and inf, respectively, and
       *  - current coefficient is not FP representable
       */
      if( type == CON_RANGE || type == CON_EQUAL )
         split = split || !mpqIsReal(scip_, val_[nnonz_]);

      /* index of variable */
      scipvar = (SCIP_VAR*)mono_get_var(term_get_element(term, i), 0);
      varidx = SCIPvarGetIndex(scipvar);
      assert(SCIPvarGetIndex(vars_[varidx]) == varidx);
      assert(vars_[varidx] == scipvar);
      ind_[nnonz_] = varidx;

#ifdef CREATEEXACTLPCONS_OUT
      gmp_printf("   i=%d: current nnonz<%d>, var<%s>, approx coef<%g>, stored exact coef<%Qd> (split:%d)\n", 
         nnonz_, i, SCIPvarGetName(scipvar), numb_todbl(mono_get_coeff(term_get_element(term, i))), val_[nnonz_], split);
#endif

      /* update number of inetgral nonzero coefficients */
      if( mpqIsIntegral(val_[nnonz_]) )
         nintegral_++;

      /* update minimum and maximal absolute nonzero constraint matrix, lhs, or rhs entry */
      mpq_abs(absval, val_[nnonz_]);
      if( mpq_cmp(absval, maxabsval_) > 0 )
         mpq_set(maxabsval_, absval);
      if( mpq_cmp(absval, minabsval_) < 0 )
         mpq_set(minabsval_, absval);
 
      /* update number of nonzero coefficients */
      nnonz_++;

   }
   mpq_clear(absval);


   /* update number of constraints */
   nconss_++;
      
   /* initialize start index and length of next constraint */
   beg_[nconss_] = nnonz_;  
   len_[nconss_] = 0;
   assert(beg_[nconss_-1] <= beg_[nconss_]);

   /* set length of constraint */
   len_[nconss_-1] = beg_[nconss_] - beg_[nconss_-1];

   return split;
}

/** method creates a linear constraint and is called directly from ZIMPL 
 *
 *  @note this method is used by ZIMPL from version 2.10; 
 */
Bool xlp_addcon_term(
   const char*           name,               /**< constraint name */
   ConType               type,               /**< constraint type */
   const Numb*           lhs,                /**< left hand side */
   const Numb*           rhs,                /**< right hand side */
   unsigned int          flags,              /**< special constraint flags */
   const Term*           term      
   )
{
   SCIP_Bool lhsgiven;
   SCIP_Bool rhsgiven;
   mpq_t absval;
   int  maxdegree;
   int  i;

   /* reallocate and initialize constraint information; ranged constraints might be splitted into two constraints */ 
   assert(nconss_ <= consssize_);
   if( nconss_ + 3 > consssize_ )
   {
      consssize_ = MAX(2 * consssize_, consssize_ + 3);

      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip_, &beg_, consssize_) );
      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip_, &len_, consssize_) );
      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip_, &lhs_, consssize_) );
      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip_, &rhs_, consssize_) );
      for( i = nconss_; i < consssize_; ++i )
      {
         mpq_init(lhs_[i]);
         mpq_init(rhs_[i]);
      }
   }

   /* reallocate and initialize matrix information; ranged constraints might be splitted into two constraints */ 
   assert(nnonz_ <= nonzsize_);
   if( nnonz_ + (2 * term_get_elements(term)) > nonzsize_ )
   {
      nonzsize_ = MAX(2 * nonzsize_, nonzsize_ + (2 * term_get_elements(term)));
      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip_, &val_, nonzsize_) );
      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip_, &ind_, nonzsize_) );
      for( i = nnonz_; i < nonzsize_; ++i )
      {
         mpq_init(val_[i]);
      }
   }   

   lhsgiven = FALSE;
   rhsgiven = FALSE;
   mpq_init(absval);

   /* get exact lhs and rhs for exactlp constraint handler */
   switch( type )
   {
   case CON_FREE:
      mpq_set(lhs_[nconss_], *negInfinity(conshdlrdata_));
      mpq_set(rhs_[nconss_], *posInfinity(conshdlrdata_));
      break;
   case CON_LHS:
      numb_get_mpq(lhs, lhs_[nconss_]);
      mpq_set(rhs_[nconss_], *posInfinity(conshdlrdata_));
      lhsgiven = TRUE;
      break;
   case CON_RHS:
      mpq_set(lhs_[nconss_], *negInfinity(conshdlrdata_));
      numb_get_mpq(rhs, rhs_[nconss_]);
      rhsgiven = TRUE;
      break;
   case CON_RANGE:
      numb_get_mpq(lhs, lhs_[nconss_]);
      numb_get_mpq(rhs, rhs_[nconss_]);
      lhsgiven = TRUE;
      rhsgiven = TRUE;
      break;
   case CON_EQUAL:
      numb_get_mpq(lhs, lhs_[nconss_]);
      numb_get_mpq(rhs, rhs_[nconss_]);
      assert(mpq_equal(lhs_[nconss_], rhs_[nconss_]) != 0);
      lhsgiven = TRUE;
      rhsgiven = TRUE;
      break;
   default:
      SCIPwarningMessage("invalid constraint type <%d> in ZIMPL callback xlp_addcon()\n", type);
      numb_get_mpq(lhs, lhs_[nconss_]);
      numb_get_mpq(rhs, rhs_[nconss_]);
      readerror_ = TRUE;
      break;
   }

   /* update minimum and maximal absolute nonzero constraint matrix, lhs, or rhs entry */
   if( lhsgiven && mpq_sgn(lhs_[nconss_]) != 0 )
   {
      mpq_abs(absval, lhs_[nconss_]);
      if( mpq_cmp(absval, maxabsval_) > 0 )
         mpq_set(maxabsval_, absval);
      if( mpq_cmp(absval, minabsval_) < 0 )
         mpq_set(minabsval_, absval);
   }
   if( rhsgiven && mpq_sgn(rhs_[nconss_]) != 0 )
   {
      mpq_abs(absval, rhs_[nconss_]);
      if( mpq_cmp(absval, maxabsval_) > 0 )
         mpq_set(maxabsval_, absval);
      if( mpq_cmp(absval, minabsval_) < 0 )
         mpq_set(minabsval_, absval);
   }

   mpq_clear(absval);

   if( (flags & LP_FLAG_CON_SEPAR) != 0 )
   {
      SCIPwarningMessage("xlp_addcon_term(): given nondefault CON_SEPAR flag ignored for exact solving process\n");
   }

   maxdegree = term_get_degree(term);
   if (maxdegree <= 1)
   {
      /* if the constraint gives an indicator constraint */
      if ( flags & LP_FLAG_CON_INDIC )
      {
         SCIPerrorMessage("xpl_addcon_term: exact version for indicator constraints not supported\n");
         return TRUE;
      }
      else
      {
         Bool split;

         split = storeConstraint(type, term);

         /* update number of linear constraints that need to be split in case of an FP relaxation */
         if( split )
            nsplitconss_++;

         if( SCIPuseFPRelaxation(scip_) && split )
         {
            assert(type == CON_RANGE || type == CON_EQUAL);

            /* change constraint just stored to an rhs constraints */
            mpq_set(lhs_[nconss_-1], *negInfinity(conshdlrdata_));

            /* store the constraint again as lhs constraint */
            numb_get_mpq(lhs, lhs_[nconss_]);
            mpq_set(rhs_[nconss_], *posInfinity(conshdlrdata_));
            split = storeConstraint(type, term);
            assert(split);
         }
      }
   }
   else if (maxdegree == 2)
   {
      SCIPerrorMessage("xpl_addcon_term: exact version for degree == 2 not supported\n");
      return TRUE;
   }
   else
   {
      SCIPerrorMessage("xpl_addcon_term for degree > 2 not implemented\n");
      return TRUE;
   }

   return FALSE;
}

#endif


#else /* else of #if (ZIMPL_VERSION >= 300) */


/** method creates a linear constraint and is called directly from ZIMPL
 *
 *  @note this method is used by ZIMPL up to version 2.09; 
 */
#ifndef EXACTSOLVE
Con* xlp_addcon(
   const char*           name,               /**< constraint name */
   ConType               type,               /**< constraint type */
   const Numb*           lhs,                /**< left hand side */
   const Numb*           rhs,                /**< right hand side */
   unsigned int          flags               /**< special constraint flags */
   )
{
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
   Con* zplcon;

   switch( type )
   {
   case CON_FREE:
      sciplhs = -SCIPinfinity(scip_);
      sciprhs = SCIPinfinity(scip_);
      break;
   case CON_LHS:
      sciplhs = (SCIP_Real)numb_todbl(lhs);
      sciprhs = SCIPinfinity(scip_);
      break;
   case CON_RHS:
      sciplhs = -SCIPinfinity(scip_);
      sciprhs = (SCIP_Real)numb_todbl(rhs);
      break;
   case CON_RANGE:
      sciplhs = (SCIP_Real)numb_todbl(lhs);
      sciprhs = (SCIP_Real)numb_todbl(rhs);
      break;
   case CON_EQUAL:
      sciplhs = (SCIP_Real)numb_todbl(lhs);
      sciprhs = (SCIP_Real)numb_todbl(rhs);
      assert(sciplhs == sciprhs);/*lint !e777 */
      break;
   default:
      SCIPwarningMessage("invalid constraint type <%d> in ZIMPL callback xlp_addcon()\n", type);
      sciplhs = (SCIP_Real)numb_todbl(lhs);
      sciprhs = (SCIP_Real)numb_todbl(rhs);
      readerror_ = TRUE;
      break;
   }

   initial = !((flags & LP_FLAG_CON_SEPAR) != 0);
   separate = TRUE;
   enforce = TRUE;
   check = enforce;
   propagate = TRUE;
   local = FALSE;
   modifiable = FALSE;
   dynamic = ((flags & LP_FLAG_CON_SEPAR) != 0);
   removable = dynamic;

   SCIP_CALL_ABORT( SCIPcreateConsLinear(scip_, &cons, name, 0, NULL, NULL, sciplhs, sciprhs,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
   zplcon = (Con*)cons; /* this is ugly, because our CONS-pointer will be released; but in this case we know that the
                         * CONS will not be destroyed by SCIPreleaseCons() */
   SCIP_CALL_ABORT( SCIPaddCons(scip_, cons) );
   SCIP_CALL_ABORT( SCIPreleaseCons(scip_, &cons) );

   return zplcon;
}
#else
Con* xlp_addcon(
   const char*           name,               /**< constraint name */
   ConType               type,               /**< constraint type */
   const Numb*           lhs,                /**< left hand side */
   const Numb*           rhs,                /**< right hand side */
   unsigned int          flags               /**< special constraint flags */
   )
{
   SCIPerrorMessage("xlp_addcon: exact version not supported.\n");
   readerror_ = TRUE;
   return NULL;
}
#endif

/** adds coefficient/variable to linear constraint 
 *
 *  @note this method is used by ZIMPL up to version 2.09; 
 */
#ifndef EXACTSOLVE
void xlp_addtonzo(
   Var*            var,                /**< variable to add */
   Con*            con,                /**< constraint */
   const Numb*     numb                /**< variable coefficient */
   )
{
   SCIP_CONS* scipcons;
   SCIP_VAR* scipvar;
   SCIP_Real scipval;

   scipcons = (SCIP_CONS*)con;
   scipvar = (SCIP_VAR*)var;
   scipval = numb_todbl(numb);

   SCIP_CALL_ABORT( SCIPaddCoefLinear(scip_, scipcons, scipvar, scipval) );
}
#else
void xlp_addtonzo(
   Var*            var,                /**< variable to add */
   Con*            con,                /**< constraint */
   const Numb*     numb                /**< variable coefficient */
   )
{
   SCIPerrorMessage("xlp_addtonzo: exact version not supported.\n");
   readerror_ = TRUE;
}
#endif

#endif /* end of #if (ZIMPL_VERSION >= 300) */


/** method adds an variable; is called directly by ZIMPL */
#ifndef EXACTSOLVE
Var* xlp_addvar(
   const char*           name,               /**< variable name */
   VarClass              usevarclass,        /**< variable type */
   const Bound*          lower,              /**< lower bound */
   const Bound*          upper,              /**< upper bound */
   const Numb*           priority,           /**< branching priority */
   const Numb*           startval            /**< */
   )
{  /*lint --e{715}*/
   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_VARTYPE vartype;
   SCIP_Bool initial;
   SCIP_Bool removable;
   SCIP_Bool dynamiccols;
   SCIP_Bool usestartsol;

   Var* zplvar;
   int branchpriority;

   SCIP_CALL_ABORT( SCIPgetBoolParam(scip_, "reading/zplreader/dynamiccols", &dynamiccols) );
   SCIP_CALL_ABORT( SCIPgetBoolParam(scip_, "reading/zplreader/usestartsol", &usestartsol) );

   switch( bound_get_type(lower) )
   {
   case BOUND_VALUE:
      lb = (SCIP_Real)numb_todbl(bound_get_value(lower));
      break;
   case BOUND_INFTY:
      lb = SCIPinfinity(scip_);
      break;
   case BOUND_MINUS_INFTY:
      lb = -SCIPinfinity(scip_);
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
      ub = SCIPinfinity(scip_);
      break;
   case BOUND_MINUS_INFTY:
      ub = -SCIPinfinity(scip_);
      break;
   case BOUND_ERROR:
   default:
      SCIPerrorMessage("invalid upper bound type <%d> in ZIMPL reader\n", bound_get_type(upper));
      ub = 0.0;
      break;
   }

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
      SCIPwarningMessage("invalid variable class <%d> in ZIMPL callback xlp_addvar()\n", usevarclass);
      vartype = SCIP_VARTYPE_CONTINUOUS;
      readerror_ = TRUE;
      break;
   }
   initial = !dynamiccols;
   removable = dynamiccols;

   if( numb_is_int(priority) )
      branchpriority = numb_toint(priority);
   else
   {
      if( !issuedbranchpriowarning_ )
      {
         SCIPverbMessage(scip_, SCIP_VERBLEVEL_MINIMAL, NULL,
            "ZIMPL reader: fractional branching priorities in input - rounding down to integer values\n");
         issuedbranchpriowarning_ = TRUE;
      }
      branchpriority = (int)numb_todbl(priority);
   }

   SCIP_CALL_ABORT( SCIPcreateVar(scip_, &var, name, lb, ub, 0.0, vartype, initial, removable, NULL, NULL, NULL, NULL) );

   zplvar = (Var*)var; /* this is ugly, because our VAR-pointer will be released; but in this case we know that the VAR will not be
                          destroyed by SCIPreleaseVar() */
   SCIP_CALL_ABORT( SCIPaddVar(scip_, var) );
   SCIP_CALL_ABORT( SCIPchgVarBranchPriority(scip_, var, branchpriority) );

   if( usestartsol )
   {
      if( nstartvals_ >= startvalssize_ )
      {
         startvalssize_ *= 2;
         SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip_, &startvals_, startvalssize_) );
         SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip_, &startvars_, startvalssize_) );
      }
      assert( nstartvals_ < startvalssize_ );
      startvals_[nstartvals_] = (SCIP_Real)numb_todbl(startval);
      startvars_[nstartvals_] = var;
      nstartvals_++;
   }

   return zplvar;
}
#else
Var* xlp_addvar(
   const char*           name,               /**< variable name */
   VarClass              usevarclass,        /**< variable type */
   const Bound*          lower,              /**< lower bound */
   const Bound*          upper,              /**< upper bound */
   const Numb*           priority,           /**< branching priority */
   const Numb*           startval            /**< */
   )
{  /*lint --e{715}*/
   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_VARTYPE vartype;
   SCIP_Bool initial;
   SCIP_Bool removable;
   SCIP_Bool dynamiccols;
   SCIP_Bool usestartsol;

   Var* zplvar;
   int branchpriority;
   int i;

   SCIP_CALL_ABORT( SCIPgetBoolParam(scip_, "reading/zplreader/dynamiccols", &dynamiccols) );
   SCIP_CALL_ABORT( SCIPgetBoolParam(scip_, "reading/zplreader/usestartsol", &usestartsol) );

   /* reallocate and initialize variable specific information */ 
   assert(nvars_ <= varssize_);
   if( nvars_ == varssize_ )
   {
      varssize_ *= 2;
      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip_, &vars_, varssize_) );
      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip_, &lb_, varssize_) );
      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip_, &ub_, varssize_) );
      SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip_, &obj_, varssize_) );
      for( i = nvars_; i < varssize_; ++i )
      {
         mpq_init(lb_[i]);
         mpq_init(ub_[i]);
         mpq_init(obj_[i]);
         mpq_set_si(obj_[i], 0, 1); 
      }
   }
   assert(nvars_ < varssize_);

   /* get exact bounds for exactlp constraint handler and safe FP values for FP problem */
   SCIPdebugMessage("zimpl reader: store lower bound of variable\n");
   switch( bound_get_type(lower) )
   {
   case BOUND_VALUE:
      numb_get_mpq(bound_get_value(lower), lb_[nvars_]);
      lb = mpqGetRealRelax(scip_, lb_[nvars_], GMP_RNDD);
      break;
   case BOUND_INFTY:
      mpq_set(lb_[nvars_], *posInfinity(conshdlrdata_));
      lb = SCIPinfinity(scip_);
      break;
   case BOUND_MINUS_INFTY:
      mpq_set(lb_[nvars_], *negInfinity(conshdlrdata_));
      lb = -SCIPinfinity(scip_);
      break;
   case BOUND_ERROR:
   default:
      SCIPerrorMessage("invalid lower bound type <%d> in ZIMPL reader\n", bound_get_type(lower));
      mpq_set_si(lb_[nvars_], 0, 1);
      lb = 0.0;
      break;
   }

   SCIPdebugMessage("zimpl reader: store upper bound of variable\n");
   switch( bound_get_type(upper) )
   {
   case BOUND_VALUE:
      numb_get_mpq(bound_get_value(upper), ub_[nvars_]);
      ub = mpqGetRealRelax(scip_, ub_[nvars_], GMP_RNDU);
      break;
   case BOUND_INFTY:
      mpq_set(ub_[nvars_], *posInfinity(conshdlrdata_));
      ub = SCIPinfinity(scip_);
      break;
   case BOUND_MINUS_INFTY:
      mpq_set(ub_[nvars_], *negInfinity(conshdlrdata_));
      ub = -SCIPinfinity(scip_);
      break;
   case BOUND_ERROR:
   default:
      SCIPerrorMessage("invalid upper bound type <%d> in ZIMPL reader\n", bound_get_type(upper));
      mpq_set_si(ub_[nvars_], 0, 1);
      ub = 0.0;
      break;
   }

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
      SCIPwarningMessage("invalid variable class <%d> in ZIMPL callback xlp_addvar()\n", usevarclass);
      vartype = SCIP_VARTYPE_CONTINUOUS;
      readerror_ = TRUE;
      break;
   }

   /* update of variables with infinite or large bounds and integer variables with infinite bounds */
   if( bound_get_type(lower) == BOUND_MINUS_INFTY || bound_get_type(upper) == BOUND_INFTY  ) 
   {
      ninfbounds_++;
      if( vartype != SCIP_VARTYPE_CONTINUOUS ) 
         ninfintbounds_++;
   }
   else if( lb <= -LARGEBOUND || ub >= LARGEBOUND ) 
      nlargebounds_++;

   /* get branching priority of variable and additional information */
   initial = !dynamiccols;
   removable = dynamiccols;
   if( numb_is_int(priority) )
      branchpriority = numb_toint(priority);
   else
   {
      if( !issuedbranchpriowarning_ )
      {
         SCIPverbMessage(scip_, SCIP_VERBLEVEL_MINIMAL, NULL,
            "ZIMPL reader: fractional branching priorities in input - rounding down to integer values\n");
         issuedbranchpriowarning_ = TRUE;
      }
      branchpriority = (int)numb_todbl(priority);
   }

   /* add variable to FP problem */
   SCIP_CALL_ABORT( SCIPcreateVar(scip_, &var, name, lb, ub, 0.0, vartype, initial, removable, NULL, NULL, NULL, NULL) );
   assert(SCIPvarGetIndex(var) == nvars_);
   zplvar = (Var*)var; /* this is ugly, because our VAR-pointer will be released; but in this case we know that the VAR will not be
                          destroyed by SCIPreleaseVar() */
   SCIP_CALL_ABORT( SCIPaddVar(scip_, var) );
   SCIP_CALL_ABORT( SCIPchgVarBranchPriority(scip_, var, branchpriority) );

   if( usestartsol )
   {
      if( nstartvals_ >= startvalssize_ )
      {
         startvalssize_ *= 2;
         SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip_, &startvals_, startvalssize_) );
         SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip_, &startvars_, startvalssize_) );
      }
      assert( nstartvals_ < startvalssize_ );
      startvals_[nstartvals_] = (SCIP_Real)numb_todbl(startval);
      startvars_[nstartvals_] = var;
      nstartvals_++;
   }

#ifdef CREATEEXACTLPCONS_OUT
   SCIPdebugMessage("zimpl reader: added new variable");
   SCIPprintVar(scip_, var, NULL);
#endif

   /* store variable for exactlp constraint handler */
   vars_[nvars_] = var;
   nvars_++;

   return zplvar;
}
#endif


/** ZIMPL_VERSION is defined by ZIMPL version 3.00 and higher. ZIMPL 3.00 made same changes in the interface to SCIP. */
#if (ZIMPL_VERSION >= 300)

#ifndef EXACTSOLVE
Bool xlp_addsos_term(
   const char*           name,               /**< constraint name */
   SosType               type,               /**< SOS type */
   const Numb*           priority,           /**< priority */
   const Term*           term                /**< terms indicating sos */
   )
{
   /*lint --e{715}*/
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

      SCIP_CALL_ABORT( SCIPcreateConsSOS1(scip_, &cons, name, 0, NULL, NULL,
            initial, separate, enforce, check, propagate, local, dynamic, removable, FALSE) );
      SCIP_CALL_ABORT( SCIPaddCons(scip_, cons) );

      for (i = 0; i < term_get_elements(term); i++)
      {
         SCIP_VAR* var;
         SCIP_Real weight;

         assert( mono_is_linear(term_get_element(term, i)) );

         var = (SCIP_VAR*) mono_get_var(term_get_element(term, i), 0);
         weight = numb_todbl(mono_get_coeff(term_get_element(term, i)));

	 SCIP_CALL_ABORT( SCIPaddVarSOS1(scip_, cons, var, weight) );
      }
      SCIP_CALL_ABORT( SCIPreleaseCons(scip_, &cons) );
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

      SCIP_CALL_ABORT( SCIPcreateConsSOS2(scip_, &cons, name, 0, NULL, NULL, 
            initial, separate, enforce, check, propagate, local, dynamic, removable, FALSE) );
      SCIP_CALL_ABORT( SCIPaddCons(scip_, cons) );
      for (i = 0; i < term_get_elements(term); i++)
      {
         SCIP_VAR* var;
         SCIP_Real weight;

         assert( mono_is_linear(term_get_element(term, i)) );

         var = (SCIP_VAR*) mono_get_var(term_get_element(term, i), 0);
         weight = numb_todbl(mono_get_coeff(term_get_element(term, i)));

	 SCIP_CALL_ABORT( SCIPaddVarSOS2(scip_, cons, var, weight) );
      }
      SCIP_CALL_ABORT( SCIPreleaseCons(scip_, &cons) );
      break;
   case SOS_ERR:
   default:
      SCIPwarningMessage("invalid SOS type <%d> in ZIMPL callback xlp_addsos_term()\n", type);
      readerror_ = TRUE;
      break;
   }

   return FALSE;
}
#else
Bool xlp_addsos_term(
   const char*           name,               /**< constraint name */
   SosType               type,               /**< SOS type */
   const Numb*           priority,           /**< priority */
   const Term*           term                /**< terms indicating sos */
   )
{
   SCIPerrorMessage("xlp_addsos_term: exact version not supported.\n");
   readerror_ = TRUE;
   return TRUE;
}
#endif

#else /* else of #if (ZIMPL_VERSION >= 300) */

/** method adds SOS constraints */
#ifndef EXACTSOLVE
Sos* xlp_addsos(
   const char*           name,               /**< constraint name */
   SosType               type,               /**< SOS type */
   const Numb*           priority            /**< priority */
   )
{  /*lint --e{715}*/
   SCIP_CONS* cons;
   SCIP_Bool initial;
   SCIP_Bool separate;
   SCIP_Bool enforce;
   SCIP_Bool check;
   SCIP_Bool propagate;
   SCIP_Bool local;
   SCIP_Bool dynamic;
   SCIP_Bool removable;
   Sos* zplsos = NULL;

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

      SCIP_CALL_ABORT( SCIPcreateConsSOS1(scip_, &cons, name, 0, NULL, NULL,
            initial, separate, enforce, check, propagate, local, dynamic, removable, FALSE) );
      zplsos = (Sos*)cons; /* this is ugly, because our CONS-pointer will be released; but in this case we know that the CONS will not be
			      destroyed by SCIPreleaseCons() */
      SCIP_CALL_ABORT( SCIPaddCons(scip_, cons) );
      SCIP_CALL_ABORT( SCIPreleaseCons(scip_, &cons) );
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

      SCIP_CALL_ABORT( SCIPcreateConsSOS2(scip_, &cons, name, 0, NULL, NULL, 
            initial, separate, enforce, check, propagate, local, dynamic, removable, FALSE) );
      zplsos = (Sos*)cons; /* this is ugly, because our CONS-pointer will be released; but in this case we know that the CONS will not be
			      destroyed by SCIPreleaseCons() */
      SCIP_CALL_ABORT( SCIPaddCons(scip_, cons) );
      SCIP_CALL_ABORT( SCIPreleaseCons(scip_, &cons) );
      break;
   case SOS_ERR:
   default:
      SCIPwarningMessage("invalid SOS type <%d> in ZIMPL callback xlp_addsos()\n", type);
      readerror_ = TRUE;
      break;
   }

   return zplsos;
}
#else
Sos* xlp_addsos(
   const char*           name,               /**< constraint name */
   SosType               type,               /**< SOS type */
   const Numb*           priority            /**< priority */
   )
{  /*lint --e{715}*/
   SCIPerrorMessage("xlp_addsos: exact version not supported.\n");
   readerror_ = TRUE;
   return NULL;
}
#endif

/** adds a variable to a SOS constraint */
#ifndef EXACTSOLVE
void xlp_addtosos(
   Sos*                  sos,                /**< SOS constraint */ 
   Var*                  var,                /**< variable to add */
   const Numb*           weight              /**< weight of the variable */
   )
{
   /* this function should maybe get the type of the sos constraint, this
      would simplify the code below. */
   
   SCIP_CONS* scipcons;
   SCIP_VAR* scipvar;
   
   scipcons = (SCIP_CONS*)sos;
   scipvar = (SCIP_VAR*)var;
   
   if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(scipcons)), "SOS1") == 0 )
      SCIP_CALL_ABORT( SCIPaddVarSOS1(scip_, scipcons, scipvar, numb_todbl(weight)) );
   else
      SCIP_CALL_ABORT( SCIPaddVarSOS2(scip_, scipcons, scipvar, numb_todbl(weight)) );
}
#else
void xlp_addtosos(
   Sos*                  sos,                /**< SOS constraint */ 
   Var*                  var,                /**< variable to add */
   const Numb*           weight              /**< weight of the variable */
   )
{
   SCIPerrorMessage("xlp_addtosos: exact version not supported.\n");
   readerror_ = TRUE;
}
#endif

#endif /* end of #if (ZIMPL_VERSION >= 300) */


/** ZIMPL_VERSION is defined by ZIMPL version 3.00 and higher. ZIMPL 3.00 made same changes in the interface to SCIP. */
#if (ZIMPL_VERSION >= 300)

/** retuns the variable name */
const char* xlp_getvarname(
   const Var*            var                 /**< variable */
   )
{
   return SCIPvarGetName((SCIP_VAR*)var);
}

#endif /* end of #if (ZIMPL_VERSION >= 300) */


/** return variable type */
VarClass xlp_getclass(
   const Var*            var                 /**< variable */
   )
{
   SCIP_VAR* scipvar;

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
      SCIPwarningMessage("invalid SCIP variable type <%d> in ZIMPL callback xlp_getclass()\n", SCIPvarGetType(scipvar));
      readerror_ = TRUE;
      break;
   }

   return VAR_CON;
}

/** returns lower bound */
#ifndef EXACTSOLVE
Bound* xlp_getlower(
   const Var*            var                 /**< variable */
   )
{
   SCIP_VAR* scipvar;
   SCIP_Real lb;
   char s[SCIP_MAXSTRLEN];
   BoundType boundtype;
   Numb* numb;
   Bound* bound;

   scipvar = (SCIP_VAR*)var;
   lb = SCIPvarGetLbGlobal(scipvar);
   (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "%.20f", lb);
   numb = numb_new_ascii(s); /* ????? isn't there a method numb_new_dbl()? */
   if( SCIPisInfinity(scip_, -lb) )
      boundtype = BOUND_MINUS_INFTY;
   else if( SCIPisInfinity(scip_, lb) )
      boundtype = BOUND_INFTY;
   else
      boundtype = BOUND_VALUE;
   bound = bound_new(boundtype, numb);
   numb_free(numb);

   return bound;
}
#else 
Bound* xlp_getlower(
   const Var*            var                 /**< variable */
   )
{
   SCIPerrorMessage("xlp_getlower: exact version not supported yet.\n");
   readerror_ = TRUE;
   return NULL;
}
#endif

/** returns upper bound */
#ifndef EXACTSOLVE
Bound* xlp_getupper(
   const Var*            var                 /**< variable */
   )
{
   SCIP_VAR* scipvar;
   SCIP_Real ub;
   char s[SCIP_MAXSTRLEN];
   BoundType boundtype;
   Numb* numb = NULL;
   Bound* bound;

   scipvar = (SCIP_VAR*)var;
   ub = SCIPvarGetUbGlobal(scipvar);
   if( SCIPisInfinity(scip_, -ub) )
      boundtype = BOUND_MINUS_INFTY;
   else if( SCIPisInfinity(scip_, ub) )
      boundtype = BOUND_INFTY;
   else
   {
      boundtype = BOUND_VALUE;
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "%.20f", ub);
      numb = numb_new_ascii(s); /* ????? isn't there a method numb_new_dbl()? */
   }
   bound = bound_new(boundtype, numb);

   if (numb != NULL)
      numb_free(numb);

   return bound;
}
#else 
Bound* xlp_getupper(
   const Var*            var                 /**< variable */
   )
{
   SCIPerrorMessage("xlp_getupper: exact version not supported yet.\n");
   readerror_ = TRUE;
   return NULL;
}
#endif


void xlp_objname(const char* name)
{  /*lint --e{715}*/
   /* nothing to be done */
}

#ifndef EXACTSOLVE
void xlp_setdir(Bool minimize)
{
   SCIP_OBJSENSE objsense;

   objsense = (minimize ? SCIP_OBJSENSE_MINIMIZE : SCIP_OBJSENSE_MAXIMIZE);
   SCIP_CALL_ABORT( SCIPsetObjsense(scip_, objsense) );
}
#else
void xlp_setdir(Bool minimize)
{
   objsense_ = (minimize ? SCIP_OBJSENSE_MINIMIZE : SCIP_OBJSENSE_MAXIMIZE);
   SCIP_CALL_ABORT( SCIPsetObjsense(scip_, objsense_) );
}

#endif

/** changes objective coefficient of a variable */
#ifndef EXACTSOLVE
void xlp_addtocost(
   Var*            var,                /**< variable */
   const Numb*     cost                /**< objective coefficient */
   )
{
   SCIP_VAR* scipvar;
   SCIP_Real scipval;

   scipvar = (SCIP_VAR*)var;
   scipval = numb_todbl(cost);

   SCIP_CALL_ABORT( SCIPchgVarObj(scip_, scipvar, SCIPvarGetObj(scipvar) + scipval) );
}
#else 
void xlp_addtocost(
   Var*            var,                /**< variable */
   const Numb*     cost                /**< objective coefficient */
   )
{
   SCIP_VAR* scipvar;
   SCIP_Real scipval;
   mpq_t scipvalmpq;
   int varidx;
   
   scipvar = (SCIP_VAR*)var;
   varidx = SCIPvarGetIndex(scipvar);
   assert(SCIPvarGetIndex(vars_[varidx]) == varidx);
   assert(vars_[varidx] == scipvar);

   mpq_init(scipvalmpq);

   /* get exact objective coefficient used in exactlp constraint handler and approximate FP value used in FP problem */
   numb_get_mpq(cost, scipvalmpq);
   scipval = mpqGetRealApprox(scip_, scipvalmpq);

#ifdef CREATEEXACTLPCONS_OUT
   SCIPdebugMessage("zimpl reader: change obj<%g> of var: add<%g> as approx", SCIPvarGetObj(scipvar), scipval);
   gmp_printf(" (<%Qd> as exact)", scipvalmpq);
#endif

   /* if we want to work on an FP relaxation we have to ensure that obj coefficient stored in scip for the variable 
    * is identical to the one in the exactlp constraint handler, otherwise the dual bounding methods are not reliable
    */
   if( SCIPuseFPRelaxation(scip_) && !objneedscaling_ )
   {
      mpq_t tmp;

      mpq_init(tmp);
      mpq_set_d(tmp, scipval);

      if( !mpq_equal(scipvalmpq, tmp) )
         objneedscaling_ = TRUE;
      
      mpq_clear(tmp);
   }
   
   /* change objective coefficient of variable in exactlp constraint handler and in FP problem */
   mpq_add(obj_[varidx], obj_[varidx], scipvalmpq);
   SCIP_CALL_ABORT( SCIPchgVarObj(scip_, scipvar, SCIPvarGetObj(scipvar) + scipval) );

#ifdef CREATEEXACTLPCONS_OUT
   SCIPprintVar(scip_, scipvar, NULL);
#endif

   mpq_clear(scipvalmpq);

}
#endif

Bool xlp_presolve(void)
{
   /* nothing to be done */
   return TRUE;
}

Bool xlp_hassos(void)
{
   return TRUE;
}

#if (ZIMPL_VERSION >= 300)
#else
Bool xlp_concheck(const Con* con)
{  /*lint --e{715}*/
   return TRUE;
}
#endif


/*
 * Callback methods of reader
 */

/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeZpl NULL


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadZpl)
{  /*lint --e{715}*/
   char oldpath[SCIP_MAXSTRLEN];
   char buffer[SCIP_MAXSTRLEN];
   char compextension[SCIP_MAXSTRLEN];
   char namewithoutpath[SCIP_MAXSTRLEN];
   char* path;
   char* name;
   char* extension;
   char* compression;
   char* paramstr;
   SCIP_SOL* startsol;

   int i;
   SCIP_Bool changedir;
   SCIP_Bool usestartsol;
   SCIP_Bool success;
   
#ifdef EXACTSOLVE
   SCIP_CONS* cons;
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

   SCIP_VAR* vartmp;
   mpq_t lbtmp;
   mpq_t ubtmp;
   mpq_t objtmp;
   int idx;
   int nswitch;

   int c;
   assert(nconss_ == 0);
#endif

   SCIP_CALL( SCIPgetBoolParam(scip, "reading/zplreader/changedir", &changedir) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/zplreader/usestartsol", &usestartsol) );

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
         return SCIP_PARSEERROR;
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
         return SCIP_PARSEERROR;
      }
      /* an extra blank line should be printed separately since the buffer message handler only handle up to one line
       *  correctly */
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n");
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "base directory for ZIMPL parsing: <%s>\n", currentpath);
      /* an extra blank line should be printed separately since the buffer message handler only handle up to one line
       *  correctly */
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n");
   }

   /* set static variables (ZIMPL callbacks do not support user data) */
   scip_ = scip;
   issuedbranchpriowarning_ = FALSE;
   readerror_ = FALSE;

   if( usestartsol )
   {
      startvalssize_ = 1024;
      SCIP_CALL_ABORT( SCIPallocMemoryArray(scip_,&startvals_,startvalssize_) );
      SCIP_CALL_ABORT( SCIPallocMemoryArray(scip_,&startvars_,startvalssize_) );
   }

#ifdef EXACTSOLVE
   conshdlrdata_ = SCIPconshdlrGetData(SCIPfindConshdlr(scip, "exactlp"));
   assert(conshdlrdata_ != NULL);

   /* allocate and initialize variable specific information */ 
   varssize_ = 1024;
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip_, &vars_, varssize_) );
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip_, &lb_, varssize_) );
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip_, &ub_, varssize_) );
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip_, &obj_, varssize_) );
   for( i = 0; i < varssize_; ++i )
   {
      mpq_init(lb_[i]);
      mpq_init(ub_[i]);
      mpq_init(obj_[i]);
      mpq_set_si(obj_[i], 0, 1); 
   }

   /* allocate and initialize constraint specific information */ 
   consssize_ = 1024;
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip_, &beg_, consssize_) );
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip_, &len_, consssize_) );
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip_, &lhs_, consssize_) );
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip_, &rhs_, consssize_) );
   for( i = 0; i < consssize_; ++i )
   {
      mpq_init(lhs_[i]);
      mpq_init(rhs_[i]);
   }

   /* set initialize start index and length of first constraint */
   beg_[nconss_] = 0;  
   len_[nconss_] = 0;

   /* allocate and initialize matrix specific information */ 
   nonzsize_ = 1024; 
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip_, &val_, nonzsize_) );
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip_, &ind_, nonzsize_) );
   for( i = 0; i < nonzsize_; ++i )
   {
      mpq_init(val_[i]);
   }

   /* initialize minimum and maximal absolute nonzero constraint matrix, lhs, or rhs entry */
   mpq_init(minabsval_);
   mpq_init(maxabsval_);
   mpq_set(minabsval_, *posInfinity(conshdlrdata_));
   mpq_set(maxabsval_, *negInfinity(conshdlrdata_));
#endif

   /* get the parameter string */
   SCIP_CALL( SCIPgetStringParam(scip, "reading/zplreader/parameters", &paramstr) );
   if( strcmp(paramstr, "-") == 0 )
   {
      /* call ZIMPL parser without arguments */
#ifdef REDUCEDSOLVE
      if( !zpl_read(filename, FALSE) )
         readerror_ = TRUE;
#else
      if( !zpl_read(filename, TRUE) )
         readerror_ = TRUE;
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
         SCIP_CALL( SCIPallocBufferArray(scip, &argv[argc], len+1) );
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
      SCIP_CALL( SCIPduplicateBufferArray(scip, &argv[argc], filename, (int) strlen(filename)+1) );
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
#ifdef REDUCEDSOLVE
      if( !zpl_read_with_args(argv, argc, FALSE) )
         readerror_ = TRUE;
#else
      if( !zpl_read_with_args(argv, argc, TRUE) )
         readerror_ = TRUE;
#endif
      /* free argument memory */
      for( i = argc-1; i >= 1; --i )
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
            SCIPwarningMessage("error changing back to directory <%s>\n", oldpath);
         }
      }
   }

#ifdef EXACTSOLVE
   /* create exactlp constraint from variable, constraint and matrix information stored */
   
#ifdef CREATEEXACTLPCONS_OUT
   SCIPdebugMessage("zimpl reader: after adding all vars sort them by probindex instead of index\n");
#endif

   /* sort variable specific information by variable's probindex instead of index */
   for( c = 0; c < nconss_; ++c )
   {
      for( i = beg_[c]; i < beg_[c] + len_[c]; ++i )
      {
         assert(ind_[i] == SCIPvarGetIndex(vars_[ind_[i]]));
         ind_[i] = SCIPvarGetProbindex(vars_[ind_[i]]);
      }
   }

   mpq_init(lbtmp);
   mpq_init(ubtmp);
   mpq_init(objtmp);
   idx = 0;
   nswitch = 0;

   for( i = 0; i < nvars_-1+nswitch; ++i )
   {
      int probindex;

      assert(idx <= i);
      assert(idx < nvars_);
      assert(nswitch <= nvars_);

      probindex = SCIPvarGetProbindex(vars_[idx]);
      assert(probindex >= idx);

#ifdef CREATEEXACTLPCONS_OUT
      SCIPdebugMessage("i=%d:  on idxpos<%d> var<%s> with probidx<%d> | on probidxpos<%d> var<%s> with probidx<%d>\n", i, 
         idx, SCIPvarGetName(vars_[idx]), SCIPvarGetProbindex(vars_[idx]), 
         probindex, SCIPvarGetName(vars_[probindex]), SCIPvarGetProbindex(vars_[probindex]) );
#endif

      /* position of current variable is incorrect */      
      if( idx != probindex )
      {
         /* swap current variable with the one that is located at the current variable's correct position */  
         mpq_set(lbtmp, lb_[probindex]);
         mpq_set(ubtmp, ub_[probindex]);
         mpq_set(objtmp, obj_[probindex]);
         vartmp = vars_[probindex];

         mpq_set(lb_[probindex], lb_[idx]);
         mpq_set(ub_[probindex], ub_[idx]);
         mpq_set(obj_[probindex], obj_[idx]);
         vars_[probindex] = vars_[idx];

         mpq_set(lb_[idx], lbtmp);
         mpq_set(ub_[idx], ubtmp);
         mpq_set(obj_[idx], objtmp);
         vars_[idx] = vartmp;

         nswitch++;
#ifdef CREATEEXACTLPCONS_OUT
         SCIPdebugMessage("   swap: idxpos<%d> var<%s> with probidx<%d> | on probidxpos<%d> var<%s> with probidx<%d>\n", 
            idx, SCIPvarGetName(vars_[idx]), SCIPvarGetProbindex(vars_[idx]), 
            probindex, SCIPvarGetName(vars_[probindex]), SCIPvarGetProbindex(vars_[probindex]) );
#endif
      }
      /* position of current variable is correct */      
      else
      {
#ifdef CREATEEXACTLPCONS_OUT
         SCIPdebugMessage("NO swap: idxpos<%d> var<%s> with probidx<%d> | on probidxpos<%d> var<%s> with probidx<%d>\n", 
            idx, SCIPvarGetName(vars_[idx]), SCIPvarGetProbindex(vars_[idx]), 
            probindex, SCIPvarGetName(vars_[probindex]), SCIPvarGetProbindex(vars_[probindex]) );
#endif

         /* move to variable on next position */
         idx++;
      }
   }

   mpq_clear(objtmp);
   mpq_clear(ubtmp);
   mpq_clear(lbtmp);

#ifndef NDEBUG
   for( i = 0; i < nvars_; ++i )
   {
      assert(SCIPvarGetProbindex(vars_[i]) == i);
 
      assert(SCIPisInfinity(scip_, SCIPvarGetLbLocal(vars_[i])) || SCIPisInfinity(scip_, -SCIPvarGetLbLocal(vars_[i]))
         || mpqGetRealRelax(scip_, lb_[i], GMP_RNDD) == SCIPvarGetLbLocal(vars_[i]));

      assert(SCIPisInfinity(scip_, SCIPvarGetUbLocal(vars_[i])) || SCIPisInfinity(scip_, -SCIPvarGetUbLocal(vars_[i]))
         || mpqGetRealRelax(scip_, ub_[i], GMP_RNDU) == SCIPvarGetUbLocal(vars_[i]));

      assert(SCIPisInfinity(scip_, SCIPvarGetObj(vars_[i])) || SCIPisInfinity(scip_, -SCIPvarGetObj(vars_[i]))
         || mpqGetRealApprox(scip_, obj_[i]) == SCIPvarGetObj(vars_[i]));
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
   

   SCIP_CALL_ABORT( SCIPcreateConsExactlp(scip, &cons, consname, objsense_, nvars_, ninfbounds_, ninfintbounds_, 
         nlargebounds_, obj_, lb_, ub_, nconss_, nsplitconss_, lhs_, rhs_, nnonz_, nintegral_, beg_, len_, ind_, val_, 
         minabsval_, maxabsval_, objneedscaling_, initial, separate, enforce, check, propagate, local, modifiable, dynamic,
         removable, stickingatnode) );
   SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );
   SCIP_CALL_ABORT( SCIPreleaseCons(scip, &cons) );
#ifdef CREATEEXACTLPCONS_OUT
   SCIPdebugMessage("zimpl reader: released exactlp constraint\n");
#endif
#endif   

   if( usestartsol )
   {
      /* transform the problem such that adding primal solutions is possible */
      SCIP_CALL( SCIPtransformProb(scip) );
#ifdef EXACTSOLVE
#endif
      SCIP_CALL( SCIPcreateSol(scip, &startsol, NULL) );
      for( i = 0; i < nstartvals_; i++ )
      {
         SCIP_CALL( SCIPsetSolVal(scip, startsol, startvars_[i], startvals_[i]) );
#ifndef EXACTSOLVE
         SCIP_CALL( SCIPreleaseVar(scip, &startvars_[i]) ); /* variables are still needed */
#endif
      }
   
      success = FALSE;
      /* current exact version of SCIP only supports primal solutions found within the exactlp constraint handler */
      if( !SCIPisExactSolve(scip) )
      {
         SCIP_CALL( SCIPtrySolFree(scip, &startsol, TRUE, TRUE, TRUE, &success) );
      }
      else
      {
         SCIP_CALL( SCIPfreeSol(scip, &startsol) );
      }
      if( success && SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_FULL )
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "ZIMPL starting solution accepted\n");

      SCIPfreeMemoryArray(scip_,&startvals_);
      SCIPfreeMemoryArray(scip_,&startvars_);
      nstartvals_ = 0;
      startvalssize_ = 0;
   }
   
#ifdef EXACTSOLVE
   /* free matrix specific information */ 
   assert(nnonz_ <= nonzsize_);
   for( i = 0; i < nonzsize_; ++i )
   {
      mpq_clear(val_[i]);
   }
   mpq_clear(minabsval_);
   mpq_clear(maxabsval_);
   SCIPfreeMemoryArray(scip_, &ind_);
   SCIPfreeMemoryArray(scip_, &val_);
   nnonz_ = 0;
   nintegral_ = 0;
   nonzsize_ = 0;

   /* free constraint specific information */ 
   assert(nconss_ <= consssize_);
   for( i = 0; i < consssize_; ++i )
   {
      mpq_clear(rhs_[i]);
      mpq_clear(lhs_[i]);
   }
   SCIPfreeMemoryArray(scip_, &rhs_);
   SCIPfreeMemoryArray(scip_, &lhs_);
   SCIPfreeMemoryArray(scip_, &len_);
   SCIPfreeMemoryArray(scip_, &beg_);
   nconss_ = 0;
   nsplitconss_ = 0;
   consssize_ = 0;

   /* free variable specific information */ 
   assert(nvars_ <= varssize_);
   for( i = 0; i < varssize_; ++i )
   {
      mpq_clear(obj_[i]);
      mpq_clear(ub_[i]);
      mpq_clear(lb_[i]);
      if( i < nvars_ )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &vars_[i]) );
      } 
   }
   SCIPfreeMemoryArray(scip_, &obj_);
   SCIPfreeMemoryArray(scip_, &ub_);
   SCIPfreeMemoryArray(scip_, &lb_);
   SCIPfreeMemoryArray(scip_, &vars_);
   nvars_ = 0;
   ninfbounds_ = 0;
   nlargebounds_ = 0;
   varssize_ = 0;
#endif

   *result = SCIP_SUCCESS;

   if( readerror_ )
      return SCIP_PARSEERROR;
   else
      return SCIP_OKAY;
}


/** problem writing method of reader */
#define readerWriteZpl NULL


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
   SCIP_READERDATA* readerdata;

   /* create zpl reader data */
   readerdata = NULL;

   /* include zpl reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerFreeZpl, readerReadZpl, readerWriteZpl, readerdata) );

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
#endif

   return SCIP_OKAY;
}

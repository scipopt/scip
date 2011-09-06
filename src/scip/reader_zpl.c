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

/**@file   reader_zpl.c
 * @ingroup FILEREADERS 
 * @brief  ZIMPL model file reader
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/reader_zpl.h"

#ifdef WITH_ZIMPL

#include <unistd.h>
#include <string.h>
#include <assert.h>

#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "scip/cons_indicator.h"
#include "scip/cons_quadratic.h"
#include "scip/pub_misc.h"

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

/*
 * LP construction interface of ZIMPL
 */

/* ZIMPL does not support user data in callbacks - we have to use static variables */
static SCIP* scip_ = NULL;
static SCIP_Real* startvals_ = NULL;
static SCIP_VAR** startvars_ = NULL;
static int startvalssize_ = 0;
static int nstartvals_ = 0;
static SCIP_Bool issuedbranchpriowarning_ = FALSE;
static SCIP_Bool readerror_ = FALSE;

void xlp_alloc(const char* name, Bool need_startval)
{  /*lint --e{715}*/
   /* create problem */
   SCIP_CALL_ABORT( SCIPcreateProb(scip_, name, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );
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


/** method creates a linear constraint and is called directly from ZIMPL 
 *
 *  @note this method is used by ZIMPL from version 3.00; 
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
      assert(sciplhs == sciprhs);  /*lint !e777*/
      break;
   default:
      SCIPwarningMessage("invalid constraint type <%d> in ZIMPL callback xlp_addcon()\n", type);
      sciplhs = (SCIP_Real)numb_todbl(lhs);
      sciprhs = (SCIP_Real)numb_todbl(rhs);
      readerror_ = TRUE;
      break;
   }

   cons = NULL;
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
               /*lint -fallthrough*/
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
                  SCIP_CALL_ABORT( SCIPsetBinaryVarIndicator(scip_, cons, scipvar) );
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
                  SCIP_CALL_ABORT( SCIPsetBinaryVarIndicator(scip_, cons, scipvar) );
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
  	  
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip_, &linvar,    term_get_elements(term)) );
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip_, &quadvar1,  term_get_elements(term)) );
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip_, &quadvar2,  term_get_elements(term)) );
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip_, &lincoeff,  term_get_elements(term)) );
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip_, &quadcoeff, term_get_elements(term)) );
  	  
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

   if( cons != NULL )
   {
      SCIP_CALL_ABORT( SCIPreleaseCons(scip_, &cons) );
   }
   
   return FALSE;
}

/** method adds an variable; is called directly by ZIMPL */
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

   SCIP_CALL_ABORT( SCIPcreateVar(scip_, &var, name, lb, ub, 0.0, vartype, initial, removable, NULL, NULL, NULL, NULL, NULL) );
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
      /*lint -fallthrough*/
   default:
      SCIPwarningMessage("invalid SOS type <%d> in ZIMPL callback xlp_addsos_term()\n", type);
      readerror_ = TRUE;
      break;
   }

   return FALSE;
}

/** returns the variable name */
const char* xlp_getvarname(
   const Var*            var                 /**< variable */
   )
{
   return SCIPvarGetName((SCIP_VAR*)var);
}

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

/** returns upper bound */
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

void xlp_objname(const char* name)
{  /*lint --e{715}*/
   /* nothing to be done */
}

void xlp_setdir(Bool minimize)
{
   SCIP_OBJSENSE objsense;

   objsense = (minimize ? SCIP_OBJSENSE_MINIMIZE : SCIP_OBJSENSE_MAXIMIZE);
   SCIP_CALL_ABORT( SCIPsetObjsense(scip_, objsense) );
}

/** changes objective coefficient of a variable */
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

Bool xlp_presolve(void)
{
   /* nothing to be done */
   return TRUE;
}

Bool xlp_hassos(void)
{
   return TRUE;
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

   /* set static variables (ZIMPL callbacks do not support user data) */
   scip_ = scip;
   issuedbranchpriowarning_ = FALSE;
   readerror_ = FALSE;

   if( usestartsol )
   {
      startvalssize_ = 1024;
      SCIP_CALL_ABORT( SCIPallocMemoryArray(scip_, &startvals_, startvalssize_) );
      SCIP_CALL_ABORT( SCIPallocMemoryArray(scip_, &startvars_, startvalssize_) );
   }

   /* get the parameter string */
   SCIP_CALL( SCIPgetStringParam(scip, "reading/zplreader/parameters", &paramstr) );
   if( strcmp(paramstr, "-") == 0 )
   {
      /* call ZIMPL parser without arguments */
      if( !zpl_read(filename, TRUE) )
         readerror_ = TRUE;
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
      if( !zpl_read_with_args(argv, argc, TRUE) )
         readerror_ = TRUE;

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
            SCIPwarningMessage("error changing back to directory <%s>\n", oldpath);
         }
      }
   }

   if( usestartsol )
   {
      /* if read failed, transformProb might fail also, due to lack of a problem */
      if( nstartvals_ > 0 && !readerror_ )
      {
         /* transform the problem such that adding primal solutions is possible */
         SCIP_CALL( SCIPtransformProb(scip) );
         SCIP_CALL( SCIPcreateSol(scip, &startsol, NULL) );
         for( i = 0; i < nstartvals_; i++ )
         {
            SCIP_CALL( SCIPsetSolVal(scip, startsol, startvars_[i], startvals_[i]) );
            SCIP_CALL( SCIPreleaseVar(scip, &startvars_[i]) );
         }
   
         success = FALSE;
         SCIP_CALL( SCIPtrySolFree(scip, &startsol, FALSE, TRUE, TRUE, TRUE, &success) );
         if( success && SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_FULL )
            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "ZIMPL starting solution accepted\n");
      }
      
      SCIPfreeMemoryArray(scip_, &startvals_);
      SCIPfreeMemoryArray(scip_, &startvars_);
      nstartvals_ = 0;
      startvalssize_ = 0;
   }

   *result = SCIP_SUCCESS;

   if( readerror_ )
      return SCIP_READERROR;
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
   char extcodename[100];

   /* create zpl reader data */
   readerdata = NULL;

   /* include zpl reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerCopyZpl,
         readerFreeZpl, readerReadZpl, readerWriteZpl,
         readerdata) );

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

   (void) SCIPsnprintf(extcodename, sizeof(extcodename), "ZIMPL %d.%d.%d", ZIMPL_VERSION/100, (ZIMPL_VERSION%100)/10, ZIMPL_VERSION%10);  /*lint !e845*/
   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, extcodename, "Zuse Institute Mathematical Programming Language developed by T. Koch (zimpl.zib.de)"));
#endif

   return SCIP_OKAY;
}

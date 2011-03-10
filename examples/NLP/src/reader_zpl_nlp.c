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
 * @brief  ZIMPL model file reader
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/reader_zpl.h"

#ifdef WITH_ZIMPL

#include <unistd.h>
#include <string.h>
#include <assert.h>

#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "cons_polynomial.h"

/* include the ZIMPL headers necessary to define the LP construction interface */
#include "zimpl/bool.h"
#include "zimpl/ratlptypes.h"
#include "zimpl/mme.h"
#include "zimpl/mono.h"
#include "zimpl/xlpglue.h"

extern
Bool zpl_read(const char* filename);
extern
Bool zpl_read_with_args(int argc, char** argv);



#define READER_NAME             "zplreader"
#define READER_DESC             "file reader for ZIMPL model files"
#define READER_EXTENSION        "zpl"




/*
 * LP construction interface of ZIMPL
 */

/* ZIMPL does not support user data in callbacks - we have to use a static variables */
static SCIP* scip_ = NULL;
static SCIP_Bool issuedbranchpriowarning_ = FALSE;
static SCIP_Bool readerror_ = FALSE;

void xlp_alloc(const char* name, Bool need_startval)
{
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
{
   /* nothing to be done here */
}

void xlp_transtable(FILE* fp, LpFormat format)
{
   /* nothing to be done here */
}

void xlp_orderfile(FILE* fp, LpFormat format)
{
   /* nothing to be done here */
}

void xlp_mstfile(FILE* fp, LpFormat format)
{
   /* nothing to be done here */
}

void xlp_sosfile(FILE* fp, LpFormat format)
{
   /* nothing to be done here */
}

Bool xlp_conname_exists(const char* conname)
{
   return (SCIPfindCons(scip_, conname) != NULL);
}

Con* xlp_addcon(const char* name, ConType type, const Numb* lhs, const Numb* rhs, unsigned int flags)
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

   SCIP_CALL_ABORT( SCIPcreateConsLinear(scip_, &cons, name, 0, NULL, NULL, sciplhs, sciprhs,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
   zplcon = (Con*)cons; /* this is ugly, because our CONS-pointer will be released; but in this case we know that the CONS will not be
                           destroyed by SCIPreleaseCons() */
   SCIP_CALL_ABORT( SCIPaddCons(scip_, cons) );
   SCIP_CALL_ABORT( SCIPreleaseCons(scip_, &cons) );

   return zplcon;
}

Con* xlp_addcon_term(const char* name, ConType type, const Numb* lhs, const Numb* rhs, unsigned int flags, const Term* term)
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
   int  i;
   
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

   if (term_is_linear(term))
   {
      SCIP_CALL_ABORT( SCIPcreateConsLinear(scip_, &cons, name, 0, NULL, NULL, sciplhs, sciprhs,
                          initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
      zplcon = (Con*)cons; /* this is ugly, because our CONS-pointer will be released; but in this case we know that the CONS will not be
                           destroyed by SCIPreleaseCons() */
      SCIP_CALL_ABORT( SCIPaddCons(scip_, cons) );

      for(i = 0; i < term_get_elements(term); i++)
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
   else
   {
      Polynomial* polynomial;
      Monomial**  monomials;

      SCIP_CALL_ABORT( SCIPallocBufferArray(scip_, &monomials, term_get_elements(term)) );

      term_print(stdout, term, FALSE);
      
      for(i = 0; i < term_get_elements(term); i++)
      {
         const Mono*     mono = term_get_element(term, i);
         const MonoElem* e;
         int             k;
         SCIP_VAR**      vars;
         SCIP_Real*      power;
         SCIP_Real       coef = numb_todbl(mono_get_coeff(mono));

         SCIP_CALL_ABORT( SCIPallocBufferArray(scip_, &vars, mono_get_count(mono)) );
         SCIP_CALL_ABORT( SCIPallocBufferArray(scip_, &power, mono_get_count(mono)) );
#if 0
         for(i = 0; i < mono_get_count(mono); i++)
         {
            vars [i] = mono_get_var(mono, i);
            power[i] = numb_todbl(mono_get_power(mono, i));
         }
#else
         for(e = &mono->first, k = 0; e != NULL; e = e->next, k++)
         {
            vars [k] = (SCIP_VAR*)entry_get_var(e->entry);
            power[k] = numb_todbl(e->power);
         }
#endif
         monomialCreate(scip_, &monomials[i], mono->count, vars, power, coef);

         SCIPfreeBufferArray(scip_, &vars);
         SCIPfreeBufferArray(scip_, &power);
      }
      SCIP_CALL_ABORT( polynomialCreate(scip_, &polynomial, term_get_elements(term), monomials) );

      SCIP_CALL_ABORT( SCIPcreateConsPolynomial(scip_, &cons, name, polynomial, sciplhs, sciprhs,
                          initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
      zplcon = (Con*)cons; /* this is ugly, because our CONS-pointer will be released; but in this case we know that the CONS will not be
                           destroyed by SCIPreleaseCons() */

      SCIP_CALL_ABORT( SCIPaddCons(scip_, cons) );
      releasePolynomial(scip_, &polynomial);

      for(i = 0; i < term_get_elements(term); i++)
      {
         releaseMonomial(scip_, &monomials[i]);
      }
      SCIPfreeBufferArray(scip_, &monomials);
   }

   SCIP_CALL_ABORT( SCIPreleaseCons(scip_, &cons) );
   
   return zplcon;
}

Var* xlp_addvar(const char* name, VarClass usevarclass, const Bound* lower, const Bound* upper, const Numb* priority, const Numb* startval)
{
   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_VARTYPE vartype;
   SCIP_Bool initial;
   SCIP_Bool removable;
   SCIP_Bool dynamiccols;
   Var* zplvar;
   int branchpriority;

   SCIP_CALL_ABORT( SCIPgetBoolParam(scip_, "reading/zplreader/dynamiccols", &dynamiccols) );

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
   SCIP_CALL_ABORT( SCIPreleaseVar(scip_, &var) );

   return zplvar;
}

Sos* xlp_addsos(const char* name, SosType type, const Numb* priority)
{
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
   Sos* zplsos;

   switch( type )
   {
   case SOS_TYPE1:
      break;
   case SOS_TYPE2:
      SCIPwarningMessage("SOS type 2 is not supported by SCIP\n");
      readerror_ = TRUE;
      break;
   case SOS_ERR:
   default:
      SCIPwarningMessage("invalid SOS type <%d> in ZIMPL callback xlp_addsos()\n", type);
      readerror_ = TRUE;
      break;
   }

   initial = TRUE;
   separate = TRUE;
   enforce = TRUE;
   check = enforce;
   propagate = TRUE;
   local = FALSE;
   modifiable = FALSE;
   dynamic = FALSE;
   removable = dynamic;

   SCIP_CALL_ABORT( SCIPcreateConsSetpack(scip_, &cons, name, 0, NULL,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
   zplsos = (Sos*)cons; /* this is ugly, because our CONS-pointer will be released; but in this case we know that the CONS will not be
                           destroyed by SCIPreleaseCons() */
   SCIP_CALL_ABORT( SCIPaddCons(scip_, cons) );
   SCIP_CALL_ABORT( SCIPreleaseCons(scip_, &cons) );

   return zplsos;
}

void xlp_addtosos(Sos* sos, Var* var, const Numb* weight)
{
   SCIP_CONS* scipcons;
   SCIP_VAR* scipvar;

   scipcons = (SCIP_CONS*)sos;
   scipvar = (SCIP_VAR*)var;

   SCIP_CALL_ABORT( SCIPaddCoefSetppc(scip_, scipcons, scipvar) );
}

VarClass xlp_getclass(const Var* var)
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

Bound* xlp_getlower(const Var* var)
{
   SCIP_VAR* scipvar;
   SCIP_Real lb;
   char s[SCIP_MAXSTRLEN];
   BoundType boundtype;
   Numb* numb;
   Bound* bound;

   scipvar = (SCIP_VAR*)var;
   lb = SCIPvarGetLbGlobal(scipvar);
   SCIPsnprintf(s, SCIP_MAXSTRLEN, "%.20f", lb);
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

Bound* xlp_getupper(const Var* var)
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
      SCIPsnprintf(s, SCIP_MAXSTRLEN, "%.20f", ub);
      numb = numb_new_ascii(s); /* ????? isn't there a method numb_new_dbl()? */
   }
   bound = bound_new(boundtype, numb);

   if (numb != NULL)
      numb_free(numb);

   return bound;
}

void xlp_objname(const char* name)
{
   /* nothing to be done */
}

void xlp_setdir(Bool minimize)
{
   SCIP_OBJSENSE objsense;

   objsense = (minimize ? SCIP_OBJSENSE_MINIMIZE : SCIP_OBJSENSE_MAXIMIZE);
   SCIP_CALL_ABORT( SCIPsetObjsense(scip_, objsense) );
}

void xlp_addtonzo(Var* var, Con* con, const Numb* numb)
{
   SCIP_CONS* scipcons;
   SCIP_VAR* scipvar;
   SCIP_Real scipval;

   scipcons = (SCIP_CONS*)con;
   scipvar = (SCIP_VAR*)var;
   scipval = numb_todbl(numb);

   SCIP_CALL_ABORT( SCIPaddCoefLinear(scip_, scipcons, scipvar, scipval) );
}

void xlp_addtocost(Var* var, const Numb* cost)
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

Bool xlp_concheck(const Con* con)
{
   return TRUE;
}



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
   SCIP_Bool changedir;

   SCIP_CALL( SCIPgetBoolParam(scip, "reading/zplreader/changedir", &changedir) );

   path = NULL;
   oldpath[0] = '\0';
   if( changedir )
   {
      /* change to the directory of the ZIMPL file, s.t. paths of data files read by the ZIMPL model are relative to
       * the location of the ZIMPL file
       */
      strncpy(buffer, filename, SCIP_MAXSTRLEN-1);
      buffer[SCIP_MAXSTRLEN-1] = '\0';
      SCIPsplitFilename(buffer, &path, &name, &extension, &compression);
      if( compression != NULL )
         SCIPsnprintf(compextension, SCIP_MAXSTRLEN, ".%s", compression);
      else
         *compextension = '\0';
      SCIPsnprintf(namewithoutpath, SCIP_MAXSTRLEN, "%s.%s%s", name, extension, compextension);
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
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\nbase directory for ZIMPL parsing: <%s>\n\n", currentpath);
   }

   /* set static variables (ZIMPL callbacks do not support user data) */
   scip_ = scip;
   issuedbranchpriowarning_ = FALSE;
   readerror_ = FALSE;

   /* get the parameter string */
   SCIP_CALL( SCIPgetStringParam(scip, "reading/zplreader/parameters", &paramstr) );
   if( strcmp(paramstr, "-") == 0 )
   {
      /* call ZIMPL parser without arguments */
      if( !zpl_read(filename) )
         readerror_ = TRUE;
   }
   else
   {
      char dummy[2] = "x";
      char** argv;
      int argc;
      int p;
      int len;
      int i;

      len = strlen(paramstr);
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
      SCIP_CALL( SCIPduplicateBufferArray(scip, &argv[argc], filename, strlen(filename)+1) );
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
      if( !zpl_read_with_args(argc, argv) )
         readerror_ = TRUE;

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

   *result = SCIP_SUCCESS;

   if( readerror_ )
      return SCIP_READERROR;
   else
      return SCIP_OKAY;
}

#endif




/*
 * reader specific interface methods
 */

/** includes the zpl file reader in SCIP */
SCIP_DECL_INCLUDEPLUGIN(SCIPincludeReaderZpl)
{
#ifdef WITH_ZIMPL
   SCIP_READERDATA* readerdata;

   /* create zpl reader data */
   readerdata = NULL;
   
   /* include zpl reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         SCIPincludeReaderZpl,
         readerFreeZpl, readerReadZpl, readerdata) );

   /* add zpl reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/zplreader/dynamiccols", "should columns be added and removed dynamically to the LP?",
         NULL, FALSE, FALSE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, 
         "reading/zplreader/changedir", "should the current directory be changed to that of the ZIMPL file before parsing?",
         NULL, FALSE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddStringParam(scip,
         "reading/zplreader/parameters", "additional parameter string passed to the ZIMPL parser (or - for no additional parameters)",
         NULL, FALSE, "-", NULL, NULL) );
#endif

   return SCIP_OKAY;
}

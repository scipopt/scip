/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: reader_zpl.c,v 1.10 2005/09/27 11:20:53 bzfpfend Exp $"

/**@file   reader_zpl.c
 * @brief  ZIMPL model file reader
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <unistd.h>
#include <string.h>
#include <assert.h>

#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "scip/reader_zpl.h"


#ifdef WITH_ZIMPL

/* include the ZIMPL headers necessary to define the LP construction interface */
#include "zimpl/bool.h"
#include "zimpl/ratlptypes.h"
#include "zimpl/mme.h"
#include "zimpl/xlpglue.h"

extern
Bool zpl_read(const char* filename);



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

void xlp_alloc(const char* name)
{
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

void xlp_write(FILE* fp, LpFormat format)
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
   SCIP_Bool removeable;
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
   removeable = dynamic;

   SCIP_CALL_ABORT( SCIPcreateConsLinear(scip_, &cons, name, 0, NULL, NULL, sciplhs, sciprhs,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removeable) );
   zplcon = (Con*)cons; /* this is ugly, because our CONS-pointer will be released; but in this case we know that the CONS will not be
                           destroyed by SCIPreleaseCons() */
   SCIP_CALL_ABORT( SCIPaddCons(scip_, cons) );
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
   SCIP_Bool removeable;
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
   case VAR_BIN:
      vartype = SCIP_VARTYPE_BINARY;
      break;
   default:
      SCIPwarningMessage("invalid variable class <%d> in ZIMPL callback xlp_addvar()\n", usevarclass);
      vartype = SCIP_VARTYPE_CONTINUOUS;
      readerror_ = TRUE;
      break;
   }
   initial = !dynamiccols;
   removeable = dynamiccols;

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

   SCIP_CALL_ABORT( SCIPcreateVar(scip_, &var, name, lb, ub, 0.0, vartype, initial, removeable, NULL, NULL, NULL, NULL) );
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
   SCIP_Bool removeable;
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
   removeable = dynamic;

   SCIP_CALL_ABORT( SCIPcreateConsSetpack(scip_, &cons, name, 0, NULL,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removeable) );
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
      return VAR_BIN;
   case SCIP_VARTYPE_INTEGER:
      return VAR_INT;
   case SCIP_VARTYPE_IMPLINT:
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
   sprintf(s, "%.20f", lb);
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
   Numb* numb;
   Bound* bound;

   scipvar = (SCIP_VAR*)var;
   ub = SCIPvarGetUbGlobal(scipvar);
   sprintf(s, "%.20f", ub);
   numb = numb_new_ascii(s); /* ????? isn't there a method numb_new_dbl()? */
   if( SCIPisInfinity(scip_, -ub) )
      boundtype = BOUND_MINUS_INFTY;
   else if( SCIPisInfinity(scip_, ub) )
      boundtype = BOUND_INFTY;
   else
      boundtype = BOUND_VALUE;
   bound = bound_new(boundtype, numb);
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

   SCIP_CALL_ABORT( SCIPchgVarObj(scip_, scipvar, scipval) );
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

/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeZpl NULL


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadZpl)
{  /*lint --e{715}*/
   char oldpath[SCIP_MAXSTRLEN];
   char buffer[SCIP_MAXSTRLEN];
   char namewithoutpath[SCIP_MAXSTRLEN];
   char compextension[SCIP_MAXSTRLEN];
   char* path;
   char* name;
   char* extension;
   char* compression;

   /* change to the directory of the ZIMPL file, s.t. paths of data files read by the ZIMPL model are relative to
    * the location of the ZIMPL file
    */
   strncpy(buffer, filename, SCIP_MAXSTRLEN-1);
   buffer[SCIP_MAXSTRLEN-1] = '\0';
   SCIPsplitFilename(buffer, &path, &name, &extension, &compression);
   if( compression != NULL )
      sprintf(compextension, ".%s", compression);
   else
      *compextension = '\0';
   sprintf(namewithoutpath, "%s.%s%s", name, extension, compextension);
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

   /* set static variables (ZIMPL callbacks do not support user data) */
   scip_ = scip;
   issuedbranchpriowarning_ = FALSE;
   readerror_ = FALSE;

   /* call ZIMPL parser */
   if( !zpl_read(namewithoutpath) )
      readerror_ = TRUE;

   /* change directory back to old path */
   if( path != NULL )
   {
      if( chdir(oldpath) != 0 )
      {
         SCIPwarningMessage("error changing back to directory <%s>\n", oldpath);
      }
   }

   *result = SCIP_SUCCESS;

   if( readerror_ )
      return SCIP_PARSEERROR;
   else
      return SCIP_OKAY;
}

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
         readerFreeZpl, readerReadZpl, readerdata) );

   /* add zpl reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/zplreader/dynamiccols", "should columns be added and removed dynamically to the LP?",
         NULL, FALSE, NULL, NULL) );
#endif

   return SCIP_OKAY;
}

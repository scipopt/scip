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
#pragma ident "@(#) $Id: reader_zpl.c,v 1.1 2005/08/26 13:01:24 bzfpfend Exp $"

/**@file   reader_zpl.c
 * @brief  ZIMPL model file reader
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "scip/reader_zpl.h"


#define READER_NAME             "zplreader"
#define READER_DESC             "zpl file reader"
#define READER_EXTENSION        "zpl"




/*
 * Data structures
 */

/* TODO: (optional) fill in the necessary reader data */

/** data for zpl reader */
struct SCIP_ReaderData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of ZIMPL
 */

typedef double Numb;
typedef double Bound;
typedef unsigned int Bool;
typedef int LpFormat;
typedef int ConType;
typedef int VarClass;
typedef void Sos;
typedef void Var;
typedef void Con;

#define ZIMPL_DECL_XLP_ALLOC(x) void x(const char* name);
#define ZIMPL_DECL_XLP_FREE(x) void x(void);
#define ZIMPL_DECL_XLP_STAT(x) void x(void);
#define ZIMPL_DECL_XLP_SCALE(x) void x(void);
#define ZIMPL_DECL_XLP_WRITE(x) void x(FILE* fp, LpFormat format);
#define ZIMPL_DECL_XLP_TRANSTABLE(x) void x(FILE* fp, LpFormat format);
#define ZIMPL_DECL_XLP_ORDERFILE(x) void x(FILE* fp, LpFormat format);
#define ZIMPL_DECL_XLP_MSTFILE(x) void x(FILE* fp, LpFormat format);
#define ZIMPL_DECL_XLP_SOSFILE(x) void x(FILE* fp, LpFormat format);
#define ZIMPL_DECL_XLP_CONNAME_EXISTS(x) Bool x(const char* conname);
#define ZIMPL_DECL_XLP_ADDCON(x) Con* x(const char* name, ConType type, const Numb* lhs, const Numb* rhs, unsigned int flags);
#define ZIMPL_DECL_XLP_ADDVAR(x) Var* x(const char* name, VarClass usevarclass, const Bound* lower, const Bound* upper, const Numb* priority, const Numb* startval);
#define ZIMPL_DECL_XLP_ADDSOS(x) Sos* x(const char* name, SosType type, const Numb* priority);
#define ZIMPL_DECL_XLP_ADDTOSOS(x) void x(Sos* sos, Var* var, const Numb* weight);
#define ZIMPL_DECL_XLP_GETCLASS(x) VarClass x(const Var* var);
#define ZIMPL_DECL_XLP_GETLOWER(x) Bound* x(const Var* var);
#define ZIMPL_DECL_XLP_GETUPPER(x) Bound* x(const Var* var);
#define ZIMPL_DECL_XLP_OBJNAME(x) void x(const char* name);
#define ZIMPL_DECL_XLP_SETDIR(x) void x(Bool minimize);
#define ZIMPL_DECL_XLP_ADDTONZO(x) void x(Var* var, Con* con, const Numb* numb);
#define ZIMPL_DECL_XLP_ADDTOCOST(x) void x(Var* var, const Numb* cost);
#define ZIMPL_DECL_XLP_PRESOLVE(x) void x(void);
#define ZIMPL_DECL_XLP_HASSOS(x) Bool x(void);

/* ZIMPL does not support user data in callbacks - we have to use a static variables */
static SCIP* scip_ = NULL;
static SCIP_Bool issuedbranchpriowarning_ = FALSE;
static SCIP_Bool readerror_ = FALSE;

static
ZIMPL_DECL_XLP_ALLOC(zimplXlpAlloc)
{
   /* create problem */
   SCIP_CALL_ABORT( SCIPcreateProb(scip, name, NULL, NULL, NULL, NULL, NULL, NULL) );
}

static
ZIMPL_DECL_XLP_FREE(zimplXlpFree)
{
   /* nothing to be done */
}

static
ZIMPL_DECL_XLP_STAT(zimplXlpStat)
{
   /* nothing to be done */
}

static
ZIMPL_DECL_XLP_SCALE(zimplXlpScale)
{
   /* nothing to be done */
}

static
ZIMPL_DECL_XLP_WRITE(zimplXlpWrite)
{
   /* nothing to be done */
}

static
ZIMPL_DECL_XLP_TRANSTABLE(zimplXlpTranstable)
{
   /* nothing to be done */
}

static
ZIMPL_DECL_XLP_ORDERFILE(zimplXlpOrderfile)
{
   /* nothing to be done */
}

static
ZIMPL_DECL_XLP_MSTFILE(zimplXlpMstfile)
{
   /* nothing to be done */
}

static
ZIMPL_DECL_XLP_SOSFILE(zimplXlpSosfile)
{
   /* nothing to be done */
}

static
ZIMPL_DECL_XLP_CONNAME_EXISTS(zimplXlpConname_Exists)
{
   return (SCIPfindCons(scip_, conname) != NULL);
}

static
ZIMPL_DECL_XLP_ADDCON(zimplXlpAddcon)
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

static
ZIMPL_DECL_XLP_ADDVAR(zimplXlpAddvar)
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

   SCIP_CALL( SCIPgetBoolParam(scip_, "reading/zplreader/dynamiccols", &dynamiccols) );

   lb = (SCIP_Real)numb_todbl(bound_get_value(lower));
   ub = (SCIP_Real)numb_todbl(bound_get_value(upper));
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
         SCIP_CALL_ABORT( SCIPverbMessage(scip_, SCIP_VERBLEVEL_MINIMAL, NULL,
               "ZIMPL reader: fractional branching priorities in input - rounding down to integer values\n") );
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

static
ZIMPL_DECL_XLP_ADDSOS(zimplXlpAddsos)
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

   return zplcon;
}

static
ZIMPL_DECL_XLP_ADDTOSOS(zimplXlpAddtosos)
{
   SCIP_CONS* scipcons;
   SCIP_VAR* scipvar;

   scipcons = (SCIP_CONS*)sos;
   scipvar = (SCIP_VAR*)var;

   SCIP_CALL_ABORT( SCIPaddCoefSetppc(scip_, scipcons, scipvar) );
}

static
ZIMPL_DECL_XLP_GETCLASS(zimplXlpGetclass)
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

static
ZIMPL_DECL_XLP_GETLOWER(zimplXlpGetlower)
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

static
ZIMPL_DECL_XLP_GETUPPER(zimplXlpGetupper)
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

static
ZIMPL_DECL_XLP_OBJNAME(zimplXlpObjname)
{
   /* nothing to be done */
}

static
ZIMPL_DECL_XLP_SETDIR(zimplXlpSetdir)
{
   SCIP_OBJSENSE objsense;

   objsense = (minimize ? SCIP_OBJSENSE_MINIMIZE : SCIP_OBJSENSE_MAXIMIZE);
   SCIP_CALL_ABORT( SCIPsetObjsense(scip_, objsense) );
}

static
ZIMPL_DECL_XLP_ADDTONZO(zimplXlpAddtonzo)
{
   SCIP_CONS* scipcons;
   SCIP_VAR* scipvar;
   SCIP_REAL scipval;

   scipcons = (SCIP_CONS*)con;
   scipvar = (SCIP_VAR*)var;
   scipval = numb_todbl(numb);

   SCIP_CALL_ABORT( SCIPaddCoefLinear(scip_, scipcons, scipvar, scipval) );
}

static
ZIMPL_DECL_XLP_ADDTOCOST(zimplXlpAddtocost)
{
   SCIP_VAR* scipvar;
   SCIP_REAL scipval;

   scipvar = (SCIP_VAR*)var;
   scipval = numb_todbl(cost);

   SCIP_CALL_ABORT( SCIPchgVarObj(scip_, scipvar, scipval) );
}

static
ZIMPL_DECL_XLP_PRESOLVE(zimplXlpPresolve)
{
   /* nothing to be done */
}

static
ZIMPL_DECL_XLP_HASSOS(zimplXlpHassos)
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
   ZIMPL* zimpl;
   SCIP_Bool error;

   /* set static variables (ZIMPL callbacks do not support user data) */
   scip_ = scip;
   issuedbranchpriowarning_ = FALSE;
   readerror_ = FALSE;

   /* create ZIMPL data */
   SCIP_ALLOC( zimpl = zimpl_create() );
   
   /* install ZIMPL callbacks */
   zimpl_install_xlp_callbacks(zimpl, zimplXlpAlloc, zimplXlpFree, zimplXlpStat, zimplXlpScale,
      zimplXlpWrite, zimplXlpTranstable, zimplXlpOrderfile, zimplXlpMstfile, zimplXlpSosfile,
      zimplXlpConname_Exists, zimplXlpAddcon, zimplXlpAddvar, zimplXlpAddsos, zimplXlpAddtosos,
      zimplXlpGetclass, zimplXlpGetlower, zimplXlpGetupper, zimplXlpObjname, zimplXlpSetdir,
      zimplXlpAddtonzo, zimplXlpAddtocost, zimplXlpPresolve, zimplXlpHassos);

   /* call ZIMPL parser */
   error = (SCIP_Bool)zimpl_exec_prog(zimpl, filename);
   if( error || readerror_ )
   {
      SCIPwarningMessage("error reading ZIMPL file <%s>\n", filename);
      return SCIP_READERROR;
   }
   
   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}




/*
 * reader specific interface methods
 */

/** includes the zpl file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderZpl(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
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

   return SCIP_OKAY;
}

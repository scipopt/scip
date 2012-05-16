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

/**@file   solex.c
 * @brief  methods for storing exact primal CIP solutions
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/scip.h"
#include "scip/set.h"
#include "scip/misc.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/solex.h"
#include "scip/cons_exactlp.h"

#ifndef NDEBUG
#include "scip/struct_solex.h"
#endif

#ifdef WITH_EXACTSOLVE
#include "gmp.h"

/** clears solution arrays of exact primal CIP solution */
static
SCIP_RETCODE solexClearArrays(
   SCIP_SOLEX*           sol                 /**< exact primal CIP solution */
   )
{
   assert(sol != NULL);

   SCIP_CALL( SCIPmpqarrayClear(sol->vals) );
   SCIP_CALL( SCIPboolarrayClear(sol->valid) );

   return SCIP_OKAY;
}

/** sets value of variable in the solution's array */
static
SCIP_RETCODE solexSetArrayVal(
   SCIP_SOLEX*           sol,                /**< exact primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   const mpq_t           val                 /**< value to set variable to */
   )
{
   int idx;

   assert(sol != NULL);

   idx = SCIPvarGetIndex(var);

   /* mark the variable valid */
   SCIP_CALL( SCIPboolarraySetVal(sol->valid, set, idx, TRUE) );

   /* set the value in the solution array */
   SCIP_CALL( SCIPmpqarraySetVal(sol->vals, set, idx, val) );

   return SCIP_OKAY;
}

/** returns the value of the variable in the given solution */
static
void solexGetArrayVal(
   SCIP_SOLEX*           sol,                /**< primal CIP solution */
   SCIP_VAR*             var,                /**< problem variable */
   mpq_t                 val                 /**< pointer to store value of variable */
   )
{
   int idx;

   assert(sol != NULL);

   idx = SCIPvarGetIndex(var);

   /* check, if the variable's value is valid */
   if( SCIPboolarrayGetVal(sol->valid, idx) )
   {
      SCIPmpqarrayGetVal(sol->vals, idx, val);
   }
   else
   {
#ifndef NDEBUG
      mpq_t tmp;

      mpq_init(tmp);
      SCIPmpqarrayGetVal(sol->vals, idx, tmp);
      assert(mpq_sgn(tmp) == 0);
#endif

      /* return the variable's value corresponding to the origin */
      switch( sol->solorigin )
      {
      case SCIP_SOLORIGIN_ORIGINAL:
      case SCIP_SOLORIGIN_ZERO:
         mpq_set_d(val, 0.0);
         break;

      case SCIP_SOLORIGIN_LPSOL:
      case SCIP_SOLORIGIN_PSEUDOSOL:
      case SCIP_SOLORIGIN_UNKNOWN:
         SCIPerrorMessage("solution origin not supported in exact version <%d>\n", sol->solorigin);
         SCIPABORT();
         break;

      default:
         SCIPerrorMessage("unknown solution origin <%d>\n", sol->solorigin);
         SCIPABORT();
         break;
      }
   }
}

/** prints given value and takes care of whether string requires more than SCIP_MAXSTRLEN characters */
static
void solexPrintValue(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   mpq_t                 value               /**< given value to print */
   )
{
   char s[SCIP_MAXSTRLEN];
   int n;

   n = gmp_snprintf(s, SCIP_MAXSTRLEN, "%20Qd", value);
   if( n >= SCIP_MAXSTRLEN )
   {
      char* bigs;

      if( SCIPallocMemorySize(scip, &bigs, n+1) != SCIP_OKAY )
      {
         SCIPmessagePrintInfo(SCIPgetMessagehdlr(scip), "string too long\n");
      }
      else
      {
#ifndef NDEBUG
         int m;
         m = gmp_snprintf(bigs, n+1, "%20Qd", value);
         assert(m == n);
#else
         gmp_snprintf(bigs, n+1, "%20Qd", value);
#endif
         SCIPmessagePrintInfo(SCIPgetMessagehdlr(scip), bigs);
         SCIPfreeMemory(scip, &bigs);
      }
   }
   else
   {
      SCIPmessagePrintInfo(SCIPgetMessagehdlr(scip), s);
   }
}

/** creates exact primal CIP solution, initialized to zero */
SCIP_RETCODE SCIPsolexCreate(
   SCIP_SOLEX**          sol,                /**< pointer to exact primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, sol) );
   SCIP_CALL( SCIPmpqarrayCreate(&(*sol)->vals, blkmem) );
   SCIP_CALL( SCIPboolarrayCreate(&(*sol)->valid, blkmem) );
   (*sol)->heur = heur;
   (*sol)->solorigin = SCIP_SOLORIGIN_ZERO;

   mpq_init((*sol)->obj);
   mpq_set_d((*sol)->obj, 0.0);

   (*sol)->primalindex = -1;

   return SCIP_OKAY;
}

/** frees exact primal CIP solution */
SCIP_RETCODE SCIPsolexFree(
   SCIP_SOLEX**          sol,                /**< pointer to exact primal CIP solution */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(sol != NULL);
   assert(*sol != NULL);

   SCIP_CALL( solexClearArrays(*sol) );
   SCIP_CALL( SCIPmpqarrayFree(&(*sol)->vals) );
   SCIP_CALL( SCIPboolarrayFree(&(*sol)->valid) );
   BMSfreeBlockMemory(blkmem, sol);

   return SCIP_OKAY;
}

/** sets value of variable in exact primal CIP solution */
SCIP_RETCODE SCIPsolexSetVal(
   SCIP_SOLEX*           sol,                /**< exact primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable to add to solution */
   const mpq_t           obj,                /**< objective value of variable */
   const mpq_t           lb,                 /**< global lower bound of variable */
   const mpq_t           val                 /**< solution value of variable */
   )
{
   mpq_t oldval;
#ifdef SCIP_DEBUG
   char s[SCIP_MAXSTRLEN];
#endif

   assert(sol != NULL);
   assert(var != NULL);

#ifdef SCIP_DEBUG
   gmp_snprintf(s, SCIP_MAXSTRLEN, "setting value of <%s> in solution %p to %Qd\n", SCIPvarGetName(var), (void*)sol, val);
   SCIPdebugMessage(s);
#endif

   mpq_init(oldval);

   /* we want to store only values for non-fixed variables (LOOSE or COLUMN); others have to be transformed */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( sol->solorigin == SCIP_SOLORIGIN_ORIGINAL )
      {
         solexGetArrayVal(sol, var, oldval);

         if( mpq_equal(val, oldval) == 0 )
         {
            mpq_t prod;

            SCIP_CALL( solexSetArrayVal(sol, set, var, val) );

            /* update objective */
            mpq_init(prod);
            mpq_mul(prod, obj, oldval);
            mpq_sub(sol->obj, sol->obj, prod);
            mpq_mul(prod, obj, val);
            mpq_add(sol->obj, sol->obj, prod);

            mpq_clear(prod);
         }
      }
      else
      {
         SCIP_CALL( SCIPsolexSetVal(sol, set, SCIPvarGetTransVar(var), obj, lb, val) );
      }
      break;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      assert(sol->solorigin != SCIP_SOLORIGIN_ORIGINAL);

      solexGetArrayVal(sol, var, oldval);

      if( mpq_equal(val, oldval) == 0 )
      {
         mpq_t prod;

         SCIP_CALL( solexSetArrayVal(sol, set, var, val) );

         /* update objective */
         mpq_init(prod);
         mpq_mul(prod, obj, oldval);
         mpq_sub(sol->obj, sol->obj, prod);
         mpq_mul(prod, obj, val);
         mpq_add(sol->obj, sol->obj, prod);

         mpq_clear(prod);
      }
      break;

   case SCIP_VARSTATUS_FIXED:
      assert(sol->solorigin != SCIP_SOLORIGIN_ORIGINAL);

      if( mpq_equal(val, lb) == 0 )
      {
#ifdef SCIP_DEBUG
         gmp_snprintf(s, SCIP_MAXSTRLEN, "cannot set solution value for variable <%s> fixed to %Qd to different value %Qd\n",
            SCIPvarGetName(var), lb, val);
         SCIPerrorMessage(s);
#endif
         SCIPABORT();
      }
      break;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      /** @todo exiptodo: presolving extension
       *  - implement this if exact version of SCIP supports aggregated variables
       */
      SCIPerrorMessage("cannot set solution value for aggregated variable\n");
      SCIPABORT();
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot set solution value for multiple aggregated variable\n");
      SCIPABORT();
      break;

   case SCIP_VARSTATUS_NEGATED:
      /** @todo exiptodo: presolving extension
       *  - implement this if exact version of SCIP supports negated variables
       */
      SCIPerrorMessage("cannot set solution value for negated variable\n");
      SCIPABORT();
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
   }

   mpq_clear(oldval);

   return SCIP_OKAY;
}

/** returns value of variable in exact primal CIP solution */
void SCIPsolexGetVal(
   SCIP_SOLEX*           sol,                /**< exact primal CIP solution */
   SCIP_VAR*             var,                /**< variable to get value for */
   mpq_t                 val                 /**< pointer to store value of variable */
   )
{
   assert(sol != NULL);
   assert(var != NULL);

   /* only values for non fixed variables (LOOSE or COLUMN) are stored; others have to be transformed */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( sol->solorigin == SCIP_SOLORIGIN_ORIGINAL )
         solexGetArrayVal(sol, var, val);
      else
         SCIPsolexGetVal(sol, SCIPvarGetTransVar(var), val);
      break;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      assert(sol->solorigin != SCIP_SOLORIGIN_ORIGINAL);
      solexGetArrayVal(sol, var, val);
      break;

   case SCIP_VARSTATUS_FIXED:
      /** @todo exiptodo: presolving extension
       *  - implement this if exact version of SCIP supports fixed variables
       */
      SCIPerrorMessage("cannot get solution value of fixed variable\n");
      SCIPABORT();
      break;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      /** @todo exiptodo: presolving extension
       *  - implement this if exact version of SCIP supports aggregated variables
       */
      SCIPerrorMessage("cannot get solution value of aggregated variable\n");
      SCIPABORT();
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      /** @todo exiptodo: presolving extension
       *  - implement this if exact version of SCIP supports multiaggregated variables
       */
      SCIPerrorMessage("cannot get solution value of multiaggregated variable\n");
      SCIPABORT();
      break;

   case SCIP_VARSTATUS_NEGATED:
      /** @todo exiptodo: presolving extension
       *  - implement this if exact version of SCIP supports negated variables
       */
      SCIPerrorMessage("cannot get solution value of negated variable\n");
      SCIPABORT();
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
   }
}

/** gets objective value of exact primal CIP solution in transformed problem */
void SCIPsolexGetObj(
   SCIP_SOLEX*           sol,                /**< primal CIP solution */
   mpq_t                 obj                 /**< pointer to store objective value of solution */
   )
{
   assert(sol != NULL);

   /* for original solutions, sol->obj contains the external objective value */
   if( sol->solorigin == SCIP_SOLORIGIN_ORIGINAL )
   {
      /** @todo exiptodo: heuristics extension
       *  - implement this if exact version of SCIP supports getting objective value of original solutions */
      SCIPerrorMessage("cannot get objectiv value of original solution\n");
      SCIPABORT();
   }
   else
   {
      mpq_set(obj, sol->obj);
   }
}

/** @todo exiptodo: heuristics extension
 *  - if heuristics are allowed to create and add exact solutions, then we need to implement SCIPsolexCheck(). currently,
 *    we do not need it because only cons_exactlp.c adds solutions and they are valid for sure.
 *  - should it be implemented here or in cons_exactlp.c?
 */

/** returns whether the given exact solutions in transformed space are equal */
SCIP_Bool SCIPsolexsAreEqual(
   SCIP_SOLEX*           sol1,               /**< first exact primal CIP solution */
   SCIP_SOLEX*           sol2,               /**< second exact primal CIP solution */
   SCIP_PROB*            prob                /**< transformed problem data */
   )
{
   SCIP_Bool areequal;
   mpq_t val1;
   mpq_t val2;
   int v;

   assert(sol1 != NULL);
   assert(sol2 != NULL);
   assert(sol1->solorigin != SCIP_SOLORIGIN_ORIGINAL);
   assert(sol2->solorigin != SCIP_SOLORIGIN_ORIGINAL);
   assert(prob != NULL);

   areequal = TRUE;

   if( mpq_equal(sol1->obj, sol2->obj) != 0 )
      return FALSE;

   mpq_init(val1);
   mpq_init(val2);

   for( v = 0; v < prob->nvars && areequal; ++v )
   {
      SCIPsolexGetVal(sol1, prob->vars[v], val1);
      SCIPsolexGetVal(sol2, prob->vars[v], val2);
      if( mpq_equal(sol1->obj, sol2->obj) != 0 )
         areequal = FALSE;
   }

   mpq_clear(val2);
   mpq_clear(val1);

   return areequal;
}


/** outputs non-zero elements of exact solution to file stream */
SCIP_RETCODE SCIPsolexPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOLEX*           sol,                /**< primal CIP solution */
   SCIP_PROB*            prob,               /**< problem data (original or transformed) */
   SCIP_PROB*            transprob,          /**< transformed problem data or NULL (to display priced variables) */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
#if 0
   SCIP_CONS** conss;
#endif
   mpq_t solval;
   mpq_t obj;
   int v;

   assert(sol != NULL);
   assert(prob != NULL);
   assert(prob->transformed || transprob != NULL);

   assert(SCIPisExactSolve(scip));
#if 0
   conss = SCIPgetConss(scip);
   assert(conss != NULL);
#endif
   assert(SCIPgetNConss(scip) == 1);

   mpq_init(solval);
   mpq_init(obj);

   /* display variables of problem data */
   for( v = 0; v < prob->nfixedvars; ++v )
   {
      assert(prob->fixedvars[v] != NULL);

      SCIPsolexGetVal(sol, prob->fixedvars[v], solval);
      if( printzeros || mpq_sgn(solval) != 0 )
      {
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%-32s ", SCIPvarGetName(prob->fixedvars[v]));
         solexPrintValue(scip, file, solval);
#if 0
         /** @todo exiptodo: presolving extension
          *  - in addition, print objective coefficients of the variables (as in the inexact version)
          *  - the problem is that currently SCIPvarGetObjExactlp() always returns the obj coef for the transformed
          *    problem, i.e. after objective scaling
          */
         SCIPvarGetObjExactlp(conss[0], prob->fixedvars[v], obj);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, " \t(obj: ");
         solexPrintValue(scip, file, obj);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, ")");
#endif
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "\n");
      }
   }
   for( v = 0; v < prob->nvars; ++v )
   {
      assert(prob->vars[v] != NULL);

      SCIPsolexGetVal(sol, prob->vars[v], solval);

      if( printzeros || mpq_sgn(solval) != 0 )
      {
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%-32s ", SCIPvarGetName(prob->vars[v]));
         solexPrintValue(scip, file, solval);
#if 0
         SCIPvarGetObjExactlp(conss[0], prob->vars[v], obj);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, " \t(obj: ");
         solexPrintValue(scip, file, obj);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, ")");
#endif
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "\n");
      }
   }

   /* display additional priced variables (if given problem data is original problem) */
   if( !prob->transformed && sol->solorigin != SCIP_SOLORIGIN_ORIGINAL )
   {
      assert(transprob != NULL);
      for( v = 0; v < transprob->nfixedvars; ++v )
      {
         assert(transprob->fixedvars[v] != NULL);

         if( SCIPvarIsTransformedOrigvar(transprob->fixedvars[v]) )
            continue;

         SCIPsolexGetVal(sol, transprob->fixedvars[v], solval);
         if( printzeros || mpq_sgn(solval) != 0 )
         {
            SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%-32s ", SCIPvarGetName(transprob->fixedvars[v]));
            solexPrintValue(scip, file, solval);
#if 0
            SCIPvarGetObjExactlp(conss[0], transprob->fixedvars[v], obj);
            SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, " \t(obj: ");
            solexPrintValue(scip, file, obj);
            SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, ")");
#endif
            SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "\n");
         }
      }
      for( v = 0; v < transprob->nvars; ++v )
      {
         assert(transprob->vars[v] != NULL);

         if( SCIPvarIsTransformedOrigvar(transprob->vars[v]) )
            continue;

         SCIPsolexGetVal(sol, transprob->vars[v], solval);
         if( printzeros || mpq_sgn(solval) != 0 )
         {
            SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%-32s ", SCIPvarGetName(transprob->vars[v]));
            solexPrintValue(scip, file, solval);
#if 0
            SCIPvarGetObjExactlp(conss[0], transprob->vars[v], obj);
            SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, " \t(obj: ");
            solexPrintValue(scip, file, obj);
            SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, ")");
#endif
            SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "\n");
         }
      }
   }

   mpq_clear(obj);
   mpq_clear(solval);

   return SCIP_OKAY;
}

/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPsolexGetOrigin

/** gets origin of exact solution */
SCIP_SOLORIGIN SCIPsolexGetOrigin(
   SCIP_SOLEX*           sol                 /**< exact primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->solorigin;
}

#endif

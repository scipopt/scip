/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_nlp.c,v 1.23 2009/09/04 14:37:17 bzfberth Exp $"

/**@file    heur_nlp.c
 * @ingroup PRIMALHEURISTICS
 * @brief   NLP local search primal heuristic
 * @author  Stefan Vigerske
 * 
 * @todo catch changes on global bound changes (e.g., due to new incumbents) and propagate to nlpi
 * @todo set cutoff or similar
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_nlp.h"
#include "scip/nlpi.h"
#ifdef WITH_IPOPT
#include "scip/nlpi_ipopt.h"
#endif

#include "scip/cons_linear.h"
#include "scip/cons_varbound.h"
#include "scip/cons_quadratic.h"

#define HEUR_NAME             "nlp"
#define HEUR_DESC             "primal heuristic that performs a local search in a QCP after fixing integer variables"
#define HEUR_DISPCHAR         'Q'
#define HEUR_PRIORITY         -2000000
#define HEUR_FREQ             10
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_NLPI*            nlpi;               /**< NLP solver interface */
   SCIP_NLPSTATISTICS*   nlpstatistics;      /**< statistics from NLP solver */
   
   int                   nvars;              /**< number of variables in NLP */
   SCIP_VAR**            var_nlp2scip;       /**< mapping variables in NLP to SCIP variables */
   SCIP_HASHMAP*         var_scip2nlp;       /**< mapping variables in SCIP to NLP variables */
                         
   int                   ndiscrvars;         /**< number of discrete variables */
   int*                  discrvars;          /**< indices of discrete variables */
                         
   int                   nvarbndconss;       /**< number of variable bound constraints */
   SCIP_CONS**           varbndconss;        /**< variable bound constraints */
                         
   SCIP_SOL*             startcand;          /**< candidate for start point for heuristic */
   SCIP_Real             startcandviol;      /**< violation of start point candidate w.r.t. constraint that reported this candidate */
                         
   int                   nlpverblevel;       /**< verbosity level of NLP solver */
   int                   nlpiterlimit;       /**< iteration limit of NLP solver; 0 for off */
   SCIP_Real             nlptimelimit;       /**< time limit of NLP solver; 0 for off */
   SCIP_Real             resolvetolfactor;   /**< factor for feasiblity tolerance when resolving NLP due to disagreement of feasibility */
   SCIP_Bool             resolvefromscratch; /**< whether a resolve of an NLP due to disagreement of feasibility should be from the original starting point or the infeasible solution */
                         
   SCIP_Longint          iterused;           /**< number of iterations used so far */
   int                   iteroffset;         /**< number of iterations added to the contingent of the total number of iterations */
   SCIP_Real             iterquot;           /**< contingent of NLP iterations in relation to the number of nodes in SCIP */
   int                   itermin;            /**< minimal number of iterations required to start local search */
                         
   SCIP_Bool             varboundexplicit;  /**< whether variable bound constraints should be handled explicitly before solving an NLP instead of adding them as linear constraints to to the NLP */
};


/*
 * Local methods
 */

#ifdef WITH_IPOPT
/** adds linear constraints to NLP */
static
SCIP_RETCODE addLinearConstraints(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_CONSHDLR*        linconshdlr         /**< constraint handler for linear constraints */
   )
{
   SCIP_CONS**      conss;
   SCIP_CONS**      useconss;
   SCIP_Real*       lhs;
   SCIP_Real*       rhs;
   SCIP_Real*       coeff;
   int*             rowoffset;
   int*             colindex;

   int              nnz;
   int              nconss;
   int              nuseconss;
   int              i;
   int              j;
   int              k;

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(linconshdlr != NULL);
   
   nconss = SCIPconshdlrGetNConss(linconshdlr);
   conss  = SCIPconshdlrGetConss(linconshdlr);
   
   if( !nconss )
      return SCIP_OKAY;
   
   SCIP_CALL( SCIPallocBufferArray(scip, &useconss, nconss) );
   
   nnz = 0;
   nuseconss = 0;
   for( i = 0; i < nconss; ++i )
   {
      for( k = 0; k < SCIPgetNVarsLinear(scip, conss[i]); ++k )
      {
         if( SCIPvarGetType(SCIPgetVarsLinear(scip, conss[i])[k]) > SCIP_VARTYPE_INTEGER )
         {
            nnz += SCIPgetNVarsLinear(scip, conss[i]);
            useconss[nuseconss] = conss[i];
            ++nuseconss;
            break;
         }
      }
   }
   
   if( !nuseconss )
   {
      SCIPfreeBufferArray(scip, &useconss);
      return SCIP_OKAY;
   }
   
   SCIP_CALL( SCIPallocBufferArray(scip, &lhs, nuseconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhs, nuseconss) );
   /* three arrays to store linear constraints matrix as sparse matrix */ 
   SCIP_CALL( SCIPallocBufferArray(scip, &rowoffset, nuseconss+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colindex, nnz) ); /* array to hold column indices of all linear constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &coeff, nnz) ); /* array to hold coefficients of all linear constraints */
   
   j = 0;
   for( i = 0; i < nuseconss; ++i )
   {
      SCIP_VAR** vars;
      SCIP_Real* vals;
      int nvars;

      vars = SCIPgetVarsLinear(scip, useconss[i]);
      vals = SCIPgetValsLinear(scip, useconss[i]);
      nvars = SCIPgetNVarsLinear(scip, useconss[i]);
      assert(vars != NULL);
      assert(vals != NULL);

      lhs[i] = SCIPgetLhsLinear(scip, useconss[i]);
      rhs[i] = SCIPgetRhsLinear(scip, useconss[i]);
      rowoffset[i] = j;

      /* copy coefficients (vals) into long coeff array starting at position j */
      BMScopyMemoryArray(&coeff[j], vals, nvars);

      for( k = 0; k < nvars; ++k, ++j )
      {
         assert(SCIPhashmapExists(heurdata->var_scip2nlp, vars[k]));
         colindex[j] = (int) (size_t) SCIPhashmapGetImage(heurdata->var_scip2nlp, vars[k]);
      }
   }
   rowoffset[nuseconss] = nnz;
   assert(j == nnz);
   
   SCIP_CALL( SCIPnlpiAddConstraints(scip, heurdata->nlpi, nuseconss,
         lhs, rhs, rowoffset, colindex, coeff,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );
   
   SCIPfreeBufferArray(scip, &coeff);
   SCIPfreeBufferArray(scip, &colindex);
   SCIPfreeBufferArray(scip, &rowoffset);
   SCIPfreeBufferArray(scip, &rhs);
   SCIPfreeBufferArray(scip, &lhs);
   SCIPfreeBufferArray(scip, &useconss);
   
   return SCIP_OKAY;
}

/** collect variable bound constraints */
static
SCIP_RETCODE collectVarBoundConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_CONSHDLR*        varbndconshdlr      /**< constraint handler for variable bound constraints */
   )
{
   SCIP_CONS**           conss;              /* all varbound constraints */
                                            
   int                   nconss;             /* total number of varbound constraints */
   int                   nconss4nlp;         /* number of varbound constraints we have to add to NLP because vbdvar is only implicit integer */
   int                   i;                 

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(varbndconshdlr != NULL);
   
   assert(heurdata->nvarbndconss == 0);
   assert(heurdata->varbndconss == NULL);

   nconss = SCIPconshdlrGetNConss(varbndconshdlr);
   if( !nconss )
      return SCIP_OKAY;
   
   conss = SCIPconshdlrGetConss(varbndconshdlr);
   
   if( heurdata->varboundexplicit )
   {
      /* count for how many constraint the Vbdvar is only implicit integer, so we need to add constraint to NLP */
      nconss4nlp = 0;
      for( i = 0; i < nconss; ++i )
      {
         if( SCIPvarGetType(SCIPgetVbdvarVarbound(scip, conss[i])) > SCIP_VARTYPE_INTEGER )
            nconss4nlp++;
         else
         {
            SCIP_CALL( SCIPcaptureCons(scip, conss[i]) );
         }
      }
   }
   else
      nconss4nlp = nconss;

   
   if( nconss4nlp )
   {
      SCIP_Real*       lhs;
      SCIP_Real*       rhs;
      int*             rowoffset;
      int*             colindex;
      SCIP_Real*       coeff;
      int              j;
      int              k;

      heurdata->nvarbndconss = nconss - nconss4nlp;
      if( heurdata->nvarbndconss )
         SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->varbndconss, heurdata->nvarbndconss) );

      SCIP_CALL( SCIPallocBufferArray(scip, &lhs, nconss4nlp) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rhs, nconss4nlp) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rowoffset, nconss4nlp + 1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &colindex, 2*nconss4nlp) );
      SCIP_CALL( SCIPallocBufferArray(scip, &coeff, 2*nconss4nlp) );
      
      j = 0;
      k = 0;
      
      for( i = 0; i < nconss; ++i )
      {
         /* treat constraint explicitly via boundchanges */
         if( heurdata->varboundexplicit && SCIPvarGetType(SCIPgetVbdvarVarbound(scip, conss[i])) <= SCIP_VARTYPE_INTEGER )
         { 
            assert(j < heurdata->nvarbndconss);
            heurdata->varbndconss[j] = conss[i];
            ++j;
            continue;
         }        
         /* else: add constraint to NLP */
         
         /* varbound constraints: lhs <= x + c * y <= rhs */
         lhs[k] = SCIPgetLhsVarbound(scip, conss[i]);
         rhs[k] = SCIPgetRhsVarbound(scip, conss[i]);
         rowoffset[k] = 2*k;
         
         coeff[2*k] = 1.0;
         assert(SCIPhashmapExists(heurdata->var_scip2nlp, SCIPgetVarVarbound(scip, conss[i])));
         colindex[2*k] = (int) (size_t) SCIPhashmapGetImage(heurdata->var_scip2nlp, SCIPgetVarVarbound(scip, conss[i]));
         
         coeff[2*k+1] = SCIPgetVbdcoefVarbound(scip, conss[i]);
         assert(SCIPhashmapExists(heurdata->var_scip2nlp, SCIPgetVbdvarVarbound(scip, conss[i])));
         colindex[2*k+1] = (int) (size_t) SCIPhashmapGetImage(heurdata->var_scip2nlp, SCIPgetVbdvarVarbound(scip, conss[i]));

         ++k;
      }

      rowoffset[k] = 2 * nconss4nlp;
      assert(j == heurdata->nvarbndconss);
      assert(k == nconss4nlp);
      
      SCIP_CALL( SCIPnlpiAddConstraints(scip, heurdata->nlpi, nconss4nlp,
            lhs, rhs, rowoffset, colindex, coeff,
            NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );
      
      SCIPfreeBufferArray(scip, &coeff);
      SCIPfreeBufferArray(scip, &colindex);
      SCIPfreeBufferArray(scip, &rowoffset);
      SCIPfreeBufferArray(scip, &rhs);
      SCIPfreeBufferArray(scip, &lhs);
   }
   else
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &heurdata->varbndconss, conss, nconss) );
      heurdata->nvarbndconss = nconss;
   }
   
   return SCIP_OKAY;
}

/** adds quadratic constraints to NLP */
static
SCIP_RETCODE addQuadraticConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_CONSHDLR*        quadconshdlr        /**< constraint handler for quadratic constraints */
   )
{   
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(quadconshdlr != NULL);
   
   if( !SCIPconshdlrGetNConss(quadconshdlr) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPconsInitNlpiQuadratic(scip, quadconshdlr, heurdata->nlpi, SCIPconshdlrGetNConss(quadconshdlr), 
         SCIPconshdlrGetConss(quadconshdlr), heurdata->var_scip2nlp) );
   
   return SCIP_OKAY;
}

/** sets up NLP from constraints in SCIP */
static
SCIP_RETCODE setupNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data structure */
   )
{
   SCIP_Real*            varlb;
   SCIP_Real*            varub;
   SCIP_Real*            objcoeff;
   int*                  objvar;
   int                   i;
   int                   cnt;                /* counter on discrete variables */
   int                   nconshdlrs;
   SCIP_CONSHDLR**       conshdlrs;

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(heurdata->nlpi == NULL);
   assert(heurdata->var_nlp2scip == NULL);
   assert(heurdata->var_scip2nlp == NULL);

#ifdef WITH_IPOPT
   SCIP_CALL( SCIPcreateNlpSolverIpopt(scip, &heurdata->nlpi) );
#else
   SCIPerrorMessage("No NLP solver available. Cannot setup NLP.\n");
   return SCIP_ERROR;
#endif
   SCIP_CALL( SCIPnlpiInit(scip, heurdata->nlpi, "nlp") );
   
   SCIP_CALL( SCIPnlpiSetIntPar(scip, heurdata->nlpi, SCIP_NLPPAR_VERBLEVEL, heurdata->nlpverblevel) );
   if( heurdata->nlptimelimit )
   {
      SCIP_CALL( SCIPnlpiSetRealPar(scip, heurdata->nlpi, SCIP_NLPPAR_TILIM, heurdata->nlptimelimit) );
   }

   /* add variables to NLP solver; capture variables */
   heurdata->nvars = SCIPgetNVars(scip);
   heurdata->ndiscrvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &heurdata->var_nlp2scip, SCIPgetVars(scip), heurdata->nvars) );
   SCIP_CALL( SCIPhashmapCreate(&heurdata->var_scip2nlp, SCIPblkmem(scip), heurdata->nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->discrvars, heurdata->ndiscrvars) );
   
   SCIP_CALL( SCIPallocBufferArray(scip, &varlb, heurdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varub, heurdata->nvars) );
   
   cnt = 0; 
   for( i = 0; i < heurdata->nvars; ++i )
   {
      varlb[i] = SCIPvarGetLbGlobal(heurdata->var_nlp2scip[i]);
      varub[i] = SCIPvarGetUbGlobal(heurdata->var_nlp2scip[i]);
      
      SCIP_CALL( SCIPcaptureVar(scip, heurdata->var_nlp2scip[i]) );

      SCIP_CALL( SCIPhashmapInsert(heurdata->var_scip2nlp, heurdata->var_nlp2scip[i], (void*)(size_t)i) );

      if( SCIPvarGetType(heurdata->var_nlp2scip[i]) <= SCIP_VARTYPE_INTEGER) /* binary or integer */
      {
         heurdata->discrvars[cnt] = i;
         ++cnt;
      }
   }
   
   SCIP_CALL( SCIPnlpiAddVars(scip, heurdata->nlpi, heurdata->nvars, varlb, varub, NULL, NULL) );
   
   SCIPfreeBufferArray(scip, &varub);
   SCIPfreeBufferArray(scip, &varlb);
      
   /* add (linear) objective to NLP solver */
   SCIP_CALL( SCIPallocBufferArray(scip, &objcoeff, heurdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &objvar, heurdata->nvars) );

   cnt = 0;
   for( i = 0; i < heurdata->nvars; ++i )
   {
      if( SCIPvarGetObj(heurdata->var_nlp2scip[i]) )
      {
         objcoeff[cnt] = SCIPvarGetObj(heurdata->var_nlp2scip[i]) * SCIPgetObjsense(scip);
         objvar[cnt]   = i;
         ++cnt;
      }
   }
   SCIP_CALL( SCIPnlpiSetObjective(scip, heurdata->nlpi, cnt, objvar, objcoeff, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0.0) );
   
   SCIPfreeBufferArray(scip, &objvar);
   SCIPfreeBufferArray(scip, &objcoeff);
   
   /* add all constraints which can be handled */   
   nconshdlrs = SCIPgetNConshdlrs(scip);
   conshdlrs  = SCIPgetConshdlrs(scip);
   for( i = 0; i < nconshdlrs; ++i )
   {
      assert(SCIPconshdlrGetNConss(conshdlrs[i]) >= 0);
      if( SCIPconshdlrGetNConss(conshdlrs[i]) == 0 )
         continue;
      
      if( strcmp(SCIPconshdlrGetName(conshdlrs[i]), "linear") == 0 )
      {
         SCIP_CALL( addLinearConstraints(scip, heurdata, conshdlrs[i]) );
      }
      else if( strcmp(SCIPconshdlrGetName(conshdlrs[i]), "quadratic") == 0 )
      {
         SCIP_CALL( addQuadraticConstraints(scip, heurdata, conshdlrs[i]) );
      }
      else if( strcmp(SCIPconshdlrGetName(conshdlrs[i]), "varbound") == 0 )
      {
         SCIP_CALL( collectVarBoundConstraints(scip, heurdata, conshdlrs[i]) );
      }
      else if( strcmp(SCIPconshdlrGetName(conshdlrs[i]), "logicor") == 0 || strcmp(SCIPconshdlrGetName(conshdlrs[i]), "setppc") == 0 
         || strcmp(SCIPconshdlrGetName(conshdlrs[i]), "knapsack") == 0 )
      {
         /* skip because combinatorial part is fixed in NLP */ 
         SCIPdebugMessage("skip adding logicor, setppc, and knapsack constraints to NLP\n"); 
      }
      else
      { 
         /* @TODO any other constraints to consider here? */
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "skip addition of %d constraints of type %s to NLP\n", 
            SCIPconshdlrGetNConss(conshdlrs[i]), SCIPconshdlrGetName(conshdlrs[i])); 
      }      
   }
   
   SCIP_CALL( SCIPnlpStatisticsCreate(scip, &heurdata->nlpstatistics) );

   return SCIP_OKAY;
}
#endif

/** for a fixation of discrete variables, applies the variable bound constraints to the NLP */
static
SCIP_RETCODE applyVarBoundConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data structure */
   )
{
   SCIP_VAR*             var;
   SCIP_Real*            varlb;
   SCIP_Real*            varub;
   SCIP_HASHMAP*         varmap;
   SCIP_CONS*            cons;
   int*                  varidx;
   void*                 idx;

   SCIP_Real             shift;
   int                   varcnt;             /* number of variables for which we have bounds to change */
   int                   i;
  
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(heurdata->nvarbndconss == 0 || heurdata->varbndconss != NULL);
   
   if( heurdata->nvarbndconss == 0 )
      return SCIP_OKAY;
   
   SCIP_CALL( SCIPallocBufferArray(scip, &varlb, heurdata->nvarbndconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varub, heurdata->nvarbndconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varidx, heurdata->nvarbndconss) );
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scip), heurdata->nvarbndconss) );
   varcnt = 0;
   
   for( i = 0; i < heurdata->nvarbndconss; ++i )
   {
      cons = heurdata->varbndconss[i];
      var = SCIPgetVarVarbound(scip, cons);
      assert(cons != NULL);
      assert(var != NULL);

      /* integer variables have been fixed in heurExecNlp already, so we skip them here */
      if( SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
         continue;
      /* There Should be an assert here!!! ?????????????????????????? */

      idx = SCIPhashmapGetImage(varmap, var);
      if( idx == NULL )
      { 
         /* variable appeared first time */
         assert(SCIPhashmapExists(heurdata->var_scip2nlp, var));
         varidx[varcnt] = (int)(size_t)SCIPhashmapGetImage(heurdata->var_scip2nlp, var);
         varlb[varcnt] = SCIPgetLhsVarbound(scip, cons);
         varub[varcnt] = SCIPgetRhsVarbound(scip, cons);
         
         shift = SCIPgetVbdcoefVarbound(scip, cons) * SCIPgetSolVal(scip, heurdata->startcand, SCIPgetVbdvarVarbound(scip, cons));
         if( !SCIPisInfinity(scip, -varlb[varcnt]) )
            varlb[varcnt] -= shift;
         if( !SCIPisInfinity(scip,  varub[varcnt]) )
            varub[varcnt] -= shift;
         
         varlb[varcnt] = MAX(varlb[varcnt],SCIPvarGetLbGlobal(var));
         varub[varcnt] = MIN(varub[varcnt],SCIPvarGetUbGlobal(var));
         
         SCIP_CALL( SCIPhashmapInsert(varmap, var, (void*)(size_t)(varcnt+1)) );
         
         SCIPdebugMessage("%s: var %s at %d now bounded in [%g, %g] due to %s = %g\n",
            SCIPconsGetName(cons), SCIPvarGetName(var), varidx[varcnt],
            varlb[varcnt], varub[varcnt], SCIPvarGetName(SCIPgetVbdvarVarbound(scip, cons)),
            SCIPgetSolVal(scip, heurdata->startcand, SCIPgetVbdvarVarbound(scip, cons)) );
         
         ++varcnt;
      }
      else
      {
         /* variable appeared before and was stored at position (int)idx - 1 */
         SCIP_Real lhs;
         SCIP_Real rhs;
         int idx_;
         
         assert(idx_ < varcnt);
      
         lhs = SCIPgetLhsVarbound(scip, cons);
         rhs = SCIPgetRhsVarbound(scip, cons);
         idx_ = ((int)(size_t)idx) - 1;
         
         shift = SCIPgetVbdcoefVarbound(scip, cons) * SCIPgetSolVal(scip, heurdata->startcand, SCIPgetVbdvarVarbound(scip, cons));
         if( !SCIPisInfinity(scip, -lhs) )
            lhs -= shift;
         if( !SCIPisInfinity(scip,  rhs) )
            rhs -= shift;
      
         varlb[idx_] = MAX(varlb[idx_],lhs);
         varub[idx_] = MIN(varub[idx_],rhs);
   
         SCIPdebugMessage("%s: var %s at %d now bounded in [%g, %g] due to %s = %g  [updated]\n",
            SCIPconsGetName(cons), SCIPvarGetName(var), varidx[idx_],
            varlb[idx_], varub[idx_], SCIPvarGetName(SCIPgetVbdvarVarbound(scip, cons)),
            SCIPgetSolVal(scip, heurdata->startcand, SCIPgetVbdvarVarbound(scip, cons)) );
      }      
   }
   
   SCIP_CALL( SCIPnlpiChgVarBounds(scip, heurdata->nlpi, varcnt, varidx, varlb, varub) );

   SCIPhashmapFree(&varmap);   
   SCIPfreeBufferArray(scip, &varidx);
   SCIPfreeBufferArray(scip, &varub);
   SCIPfreeBufferArray(scip, &varlb);

   return SCIP_OKAY;
}

/** free NLP data structure */
static
SCIP_RETCODE destroyNLP(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_HEURDATA*        heurdata            /**< heuristic data structure */
   )
{
   int i;

   assert(scip != NULL);
   assert(heurdata != NULL);
  
   assert(heurdata->nlpi != NULL);
   assert(heurdata->nlpstatistics != NULL);
   
   for( i = 0; i < heurdata->nvars; ++i )
      SCIP_CALL( SCIPreleaseVar(scip, &heurdata->var_nlp2scip[i]) );
   
   SCIPfreeMemoryArray(scip, &heurdata->var_nlp2scip);
   SCIPhashmapFree(&heurdata->var_scip2nlp);
   
   SCIPfreeMemoryArray(scip, &heurdata->discrvars);

   SCIP_CALL( SCIPnlpiFree(scip, &heurdata->nlpi) );
   assert(heurdata->nlpi == NULL);
   
   SCIPnlpStatisticsFree(scip, &heurdata->nlpstatistics);
   assert(heurdata->nlpstatistics == NULL);
   
   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeNlp)
{
   SCIP_HEURDATA* heurdata;
   assert(scip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->nlpi == NULL);
   assert(heurdata->var_nlp2scip == NULL);
   assert(heurdata->var_scip2nlp == NULL);
   assert(heurdata->discrvars == NULL);
   assert(heurdata->startcand == NULL);
   
   SCIPfreeMemoryNull(scip, &heurdata);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
#define heurInitNlp NULL


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#define heurExitNlp NULL


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolNlp)
{
   SCIP_Bool havenlp; /* first assume all continuous variables are linear */
   
   assert(scip != NULL);
   assert(heur != NULL);
   
   if( SCIPheurGetFreq(heur) < 0 )
      return SCIP_OKAY;

   havenlp = FALSE;      

   /* do not build NLP if there are no nonlinear continuous or impl. integer variables */
   if( SCIPgetNContVars(scip) > 0 || SCIPgetNImplVars(scip) > 0 )
   {  
      /* check if we have nonlinear continuous or impl. integer variables */
      SCIP_CONSHDLR* conshdlr;
      
      conshdlr = SCIPfindConshdlr(scip, "quadratic");
      if( conshdlr && SCIPconshdlrGetNConss(conshdlr) )
      {
         int i;
         int j;
         for( i = 0; !havenlp && i < SCIPconshdlrGetNConss(conshdlr); ++i )
         {
            SCIP_VAR** quadvars;
            SCIP_CONS* cons;
            int nquadvars;

            cons = SCIPconshdlrGetConss(conshdlr)[i];
            assert(cons != NULL);

            nquadvars = SCIPgetNQuadVarsQuadratic(scip, cons);
            quadvars = SCIPgetQuadVarsQuadratic(scip, cons);
            assert(nquadvars == 0 || quadvars != NULL);

            for( j = 0; !havenlp && j < nquadvars; ++j )
               if( SCIPvarGetType(quadvars[j]) == SCIP_VARTYPE_IMPLINT || SCIPvarGetType(quadvars[j]) == SCIP_VARTYPE_CONTINUOUS )
                  havenlp = TRUE;
            /* Are implicit integers really valid here? ???????????????? */
         }
      }
   }

   if( havenlp )
   {
      SCIP_HEURDATA* heurdata;

      /* get heuristic's data */
      heurdata = SCIPheurGetData(heur);
      assert( heurdata != NULL );

#ifdef WITH_IPOPT
      setupNLP(scip, heurdata);
#else
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "No NLP solver available. NLP local search heuristic will not run.\n");
#endif

   }
   else
   {
      SCIPdebugMessage("No nonlinear continuous variables. NLP local search heuristic will not run.\n");
   }
   
   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolNlp)
{
   SCIP_HEURDATA* heurdata;
   assert(scip != NULL);
   assert(heur != NULL);

   /* get heuristic's data */  
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   
   if( heurdata->nlpi )
   {
      SCIP_CALL( destroyNLP(scip, heurdata) );
   }
   
   if( heurdata->nvarbndconss )
   {
      int i;

      assert(heurdata->varbndconss != NULL);
      for( i = 0; i < heurdata->nvarbndconss; ++i )
      {
         assert(heurdata->varbndconss[i] != NULL);
         SCIP_CALL( SCIPreleaseCons(scip, &heurdata->varbndconss[i]) );
      }

      SCIPfreeMemoryArray(scip, &heurdata->varbndconss);
      heurdata->nvarbndconss = 0;
   }
   
   if( heurdata->startcand )
   {
      SCIP_CALL( SCIPfreeSol(scip, &heurdata->startcand) );
   }

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecNlp)
{  
   SCIP_HEURDATA* heurdata;
   SCIP_Real*     startpoint;
   SCIP_Longint   itercontingent;
   SCIP_Real      timelimit;

   assert(scip != NULL);
   assert(heur != NULL);
  
   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   
   *result = SCIP_DIDNOTRUN;
   
   /* probably do not have continuous or implicit integer variables in nonlinear functions */
   if( heurdata->nlpi == NULL )
      return SCIP_OKAY;
   /* Strange comment. Isn't that the case for MIP, too? ????????????????????? */
   
   if( heurdata->startcand == NULL )
   {
      /* only call heuristic if optimal LP solution is available */
      if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      {
         SCIPdebugMessage("skip NLP heuristic because no start candidate given and no LP solution available\n");
         return SCIP_OKAY;
      }

      /* only call heuristic, if there are no fractional variables */
      if( SCIPgetNLPBranchCands(scip) > 0 )
      {
         SCIPdebugMessage("skip NLP heuristic because no start candidate given and current LP solution is fractional\n");
         return SCIP_OKAY;
      }
   }

   itercontingent = (SCIP_Longint)(heurdata->iterquot * SCIPgetNNodes(scip));

   /* weight by previous success of heuristic */
   itercontingent = (SCIP_Longint)(itercontingent * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   itercontingent += heurdata->iteroffset; 
   itercontingent -= heurdata->iterused;

   if( itercontingent < heurdata->itermin )
   {
      SCIPdebugMessage("skip NLP heuristic; contingent=%"SCIP_LONGINT_FORMAT"; minimal number of iterations=%d; success ratio=%g\n", 
         itercontingent, heurdata->itermin, (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
      return SCIP_OKAY;
   }

   if( heurdata->nlpiterlimit > 0 )
      itercontingent = MIN(itercontingent, heurdata->nlpiterlimit);

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
   {
      timelimit -= SCIPgetSolvingTime(scip);
      if( timelimit < 10.0 )
      {
         SCIPdebugMessage("skip NLP heuristic; only %g seconds time left\n", timelimit);
         return SCIP_OKAY;
      }
      if( heurdata->nlptimelimit > 0 )
         timelimit = MIN(heurdata->nlptimelimit, timelimit);

      SCIP_CALL( SCIPnlpiSetRealPar(scip, heurdata->nlpi, SCIP_NLPPAR_TILIM, timelimit) );
   }

   *result = SCIP_DIDNOTFIND;
   
   SCIP_CALL( SCIPallocBufferArray(scip, &startpoint, heurdata->nvars) );

   SCIP_CALL( SCIPgetSolVals(scip, heurdata->startcand, heurdata->nvars, heurdata->var_nlp2scip, startpoint) );
   SCIP_CALL( SCIPnlpiSetInitialGuess(scip, heurdata->nlpi, startpoint) );

   /* fix discrete variables */
   if( heurdata->ndiscrvars > 0 )
   { 
      SCIP_Real* discrfix;
      int i;
      
      SCIP_CALL( SCIPallocBufferArray(scip, &discrfix, heurdata->ndiscrvars) );

      for( i = 0; i < heurdata->ndiscrvars; ++i )
      {
         discrfix[i] = startpoint[heurdata->discrvars[i]];

         /* round fractional variables to the nearest integer */
         if( !SCIPisIntegral(scip, discrfix[i]) )
            discrfix[i] = SCIPfrac(scip, discrfix[i]) > 0.5 ? SCIPceil(scip, discrfix[i]) : SCIPfloor(scip, discrfix[i]);

         /* adjust value to the global bounds of the corresponding SCIP variable */
         discrfix[i] = MAX(discrfix[i], SCIPvarGetLbGlobal(heurdata->var_nlp2scip[heurdata->discrvars[i]]));
         discrfix[i] = MIN(discrfix[i], SCIPvarGetUbGlobal(heurdata->var_nlp2scip[heurdata->discrvars[i]]));
      }
      /* apply fixings */
      SCIP_CALL( SCIPnlpiChgVarBounds(scip, heurdata->nlpi, heurdata->ndiscrvars, heurdata->discrvars, discrfix, discrfix) );

      SCIPfreeBufferArray(scip, &discrfix);
   }
   SCIP_CALL( applyVarBoundConstraints(scip, heurdata) );
   
   SCIPfreeBufferArray(scip, &startpoint);   

   if( heurdata->startcand != NULL )
      SCIPfreeSol(scip, &heurdata->startcand);
      
   SCIP_CALL( SCIPnlpiSetIntPar(scip, heurdata->nlpi, SCIP_NLPPAR_ITLIM, itercontingent) );
   SCIPdebugMessage("start NLP solve with iteration limit %"SCIP_LONGINT_FORMAT" and timelimit %g\n", itercontingent, timelimit);
   
   SCIP_CALL( SCIPnlpiSolve(scip, heurdata->nlpi) );

   SCIPdebugMessage("NLP solver returned with termination status %d and solution status %d\n", 
      SCIPnlpiGetTermstat(scip, heurdata->nlpi), SCIPnlpiGetSolstat(scip, heurdata->nlpi));
   
   if( SCIPnlpiGetTermstat(scip, heurdata->nlpi) >= SCIP_NLPITERMSTAT_MEMERR )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, 
         "NLP solver returned with bad termination status %d. Will not run NLP heuristic again for this run.\n",  
         SCIPnlpiGetTermstat(scip, heurdata->nlpi));
      SCIP_CALL( destroyNLP(scip, heurdata) );
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPnlpiGetStatistics(scip, heurdata->nlpi, heurdata->nlpstatistics) );

   heurdata->iterused += SCIPnlpStatisticsGetNIterations(heurdata->nlpstatistics);
   SCIPdebugMessage("NLP solver used %d iterations and %g seconds\n", 
      SCIPnlpStatisticsGetNIterations(heurdata->nlpstatistics), SCIPnlpStatisticsGetTotalTime(heurdata->nlpstatistics));

   if( SCIPnlpiGetSolstat(scip, heurdata->nlpi) <= SCIP_NLPSOLSTAT_FEASIBLE )
   {
      SCIP_Real* primals;
      SCIP_SOL*  sol;
      SCIP_Bool  stored; 
      SCIP_Bool  feasible;

      SCIP_CALL( SCIPnlpiGetSolution(scip, heurdata->nlpi, &primals) );
      assert(primals != NULL);

      SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
      SCIP_CALL( SCIPsetSolVals(scip, sol, heurdata->nvars, heurdata->var_nlp2scip, primals) );

      SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, TRUE, TRUE, &feasible) );

      /* add solution */
      if( feasible )
      { 
         SCIP_CALL( SCIPaddSol(scip, sol, &stored) );
         if( stored )
         {
            SCIPdebugMessage("SCIP stored solution\n");
            *result = SCIP_FOUNDSOL;
         }
      }
      else if( heurdata->resolvetolfactor < 1.0 )
      { 
         /* resolve with tighter tolerances */
         SCIPdebugMessage("solution reported by NLP solver not feasible for SCIP, resolve with feasibility tolerance %g\n", heurdata->resolvetolfactor*SCIPfeastol(scip));
         SCIP_CALL( SCIPnlpiSetRealPar(scip, heurdata->nlpi, SCIP_NLPPAR_FEASTOL, heurdata->resolvetolfactor*SCIPfeastol(scip)) );

         if( !heurdata->resolvefromscratch )
         {
            SCIP_CALL( SCIPnlpiSetInitialGuess(scip, heurdata->nlpi, primals) );
         }

         SCIP_CALL( SCIPnlpiSolve(scip, heurdata->nlpi) );
         SCIP_CALL( SCIPnlpiGetStatistics(scip, heurdata->nlpi, heurdata->nlpstatistics) );
         heurdata->iterused += SCIPnlpStatisticsGetNIterations(heurdata->nlpstatistics);

         if( SCIPnlpiGetSolstat(scip, heurdata->nlpi) <= SCIP_NLPSOLSTAT_FEASIBLE )
         { 
            /* still feasible, hope that SCIP accepts it now */
            SCIP_CALL( SCIPnlpiGetSolution(scip, heurdata->nlpi, &primals) );
            assert(primals != NULL);

            SCIP_CALL( SCIPsetSolVals(scip, sol, heurdata->nvars, heurdata->var_nlp2scip, primals) );

            /* ???????? What is this call good for? The value of feasible is not used! */
#ifdef SCIP_DEBUG
            SCIP_CALL( SCIPcheckSolOrig(scip, sol, &feasible, TRUE, TRUE) );
#endif
            SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, TRUE, &stored) );
            if( stored )
            {
               SCIPdebugMessage("SCIP stored solution\n");
               *result = SCIP_FOUNDSOL;
            }
            else
            {
               SCIPdebugMessage("NLP solver said feasible after tigthening tolerances, but SCIP still did not store solution\n");
            }
         }
         
         /* reset to original tolerance */
         SCIP_CALL( SCIPnlpiSetRealPar(scip, heurdata->nlpi, SCIP_NLPPAR_FEASTOL, SCIPfeastol(scip)) );
      }

      SCIP_CALL( SCIPfreeSol(scip, &sol) );
   }
   
   /* TODO reset time and iterlimit in nlp solver? */
   
   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the NLP local search primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurNlp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create Nlp primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   BMSclearMemory(heurdata);
   
   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, heurFreeNlp, heurInitNlp, heurExitNlp, heurInitsolNlp, heurExitsolNlp, heurExecNlp,
         heurdata) );

   /* add Nlp primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/"HEUR_NAME"/nlpverblevel", "verbosity level of NLP solver", 
         &heurdata->nlpverblevel, FALSE, 0, 0, 2, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/"HEUR_NAME"/nlpiterlimit", "iteration limit of NLP solver; 0 to use solver default", 
         &heurdata->nlpiterlimit, FALSE, 0, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/nlptimelimit", "time limit of NLP solver; 0 to use solver default", 
         &heurdata->nlptimelimit, FALSE, 0, 0, SCIPinfinity(scip), NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/resolvetolfactor", "if SCIP does not accept a solution which the NLP solver thinks is feasible, the feasibility tolerance is reduced by this factor and the NLP resolved (set to 1. to turn off resolve", 
         &heurdata->resolvetolfactor, FALSE, 0.01, 0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/resolvefromscratch", "whether a resolve of an NLP due to disagreement of feasibility should be from the original starting point or the infeasible solution", 
         &heurdata->resolvefromscratch, FALSE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/"HEUR_NAME"/iteroffset", "number of iterations added to the contingent of the total number of iterations", 
         &heurdata->iteroffset, FALSE, 500, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/iterquotient", "contingent of NLP iterations in relation to the number of nodes in SCIP", 
         &heurdata->iterquot,  FALSE, 0.1, 0, SCIPinfinity(scip), NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/"HEUR_NAME"/itermin", "contingent of NLP iterations in relation to the number of nodes in SCIP", 
         &heurdata->itermin,   FALSE, 300, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/varboundexplicit", "whether variable bound constraints should be handled explicitly before solving NLP instead of adding them to the NLP", 
         &heurdata->varboundexplicit, FALSE, TRUE, NULL, NULL) );

   return SCIP_OKAY;
}

/** updates the starting point for the NLP heuristic
 * 
 * Is called by a constraint handler that handles nonlinear constraints when a check on feasibility of a solution fails.
 */
SCIP_RETCODE SCIPheurNlpUpdateStartpoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< NLP heuristic */
   SCIP_SOL*             solcand,            /**< solution candidate */
   SCIP_Real             violation           /**< constraint violation of solution candidate */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_Bool      takenew;
   
   assert(scip != NULL);
   assert(heur != NULL);
   assert(solcand != NULL);
   assert(SCIPisPositive(scip, violation));
   
   if( SCIPgetStage(scip) >= SCIP_STAGE_SOLVED )
      return SCIP_OKAY;
   
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   
   /* probably have no continuous variables */
   if( !heurdata->nlpi )
      return SCIP_OKAY;
   
   takenew = (heurdata->startcand == NULL);
   if( heurdata->startcand == NULL || SCIPisGT(scip, heurdata->startcandviol, violation) ||
      SCIPisRelGT(scip, SCIPgetSolTransObj(scip, heurdata->startcand)*SCIPgetObjsense(scip), SCIPgetSolTransObj(scip, solcand)*SCIPgetObjsense(scip)) )
   {
      if( heurdata->startcand != NULL )
      {
         SCIP_CALL( SCIPfreeSol(scip, &heurdata->startcand) );
      }
      SCIP_CALL( SCIPcreateSolCopy(scip, &heurdata->startcand, solcand) );
      SCIP_CALL( SCIPunlinkSol(scip, heurdata->startcand) );
      heurdata->startcandviol = violation;
   }

   return SCIP_OKAY;
}

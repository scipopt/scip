/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_nlp.c,v 1.11 2009/08/06 13:09:02 bzfviger Exp $"

/**@file    heur_nlp.c
 * @ingroup PRIMALHEURISTICS
 * @brief   NLP local search primal heuristic
 * @author  Stefan Vigerske
 */

/* @TODO catch changes on global bound changes (e.g., due to new incumbents) and propagate to nlpi
 * @TODO set cutoff or similar
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
#ifdef WITH_SOC3
#include "cons_soc3.h"
#include "cons_soc.h"
#endif

#define HEUR_NAME             "nlp"
#define HEUR_DESC             "primal heuristic that performs a local search in an NLP after fixing integer variables"
#define HEUR_DISPCHAR         'X'  /* @TODO: need a nice letter */
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
   SCIP_NLPI*     nlpi;                  /**< NLP solver interface */
   SCIP_NLPSTATISTICS* nlpstatistics;    /**< statistics from NLP solver */
   
   int            nvars;                 /**< number of variables in NLP */
   SCIP_VAR**     var_nlp2scip;          /**< mapping variables in NLP to SCIP variables */
   SCIP_HASHMAP*  var_scip2nlp;          /**< mapping variables in SCIP to NLP variables */
   
   int            ndiscrvars;            /**< number of discrete variables */
   int*           discrvars;             /**< indices of discrete variables */
   
   int            nvarbndconss;          /**< number of variable bound constraints */
   SCIP_CONS**    varbndconss;           /**< variable bound constraints */
   
   SCIP_SOL*      startcandidate;        /**< candidate for start point for heuristic */
   SCIP_Real      startcand_violation;   /**< violation of start point candidate w.r.t. constraint that reported this candidate */
   
   int            nlpverblevel;          /**< verbosity level of NLP solver */
   int            nlpiterlimit;          /**< iteration limit of NLP solver; 0 for off */
   SCIP_Real      nlptimelimit;          /**< time limit of NLP solver; 0 for off */
   SCIP_Real      resolvetolfactor;      /**< factor for feasiblity tolerance when resolving NLP due to disagreement of feasibility */
   
   SCIP_Longint   iterused;              /**< number of iterations used so far */
   int            iteroffset;            /**< number of iterations added to the contingent of the total number of iterations */
   SCIP_Real      iterquot;              /**< contingent of NLP iterations in relation to the number of nodes in SCIP */
   int            itermin;               /**< minimal number of iterations required to start local search */
   
   SCIP_Bool      varbound_explicit;     /**< whether variable bound constraints should be handled explicitly before solving an NLP instead of adding them as linear constraints to to the NLP */
};


/*
 * Local methods
 */
static SCIP_RETCODE addLinearConstraints(SCIP* scip, SCIP_HEUR* heur, SCIP_CONSHDLR* linconshdlr);
static SCIP_RETCODE collectVarBoundConstraints(SCIP* scip, SCIP_HEUR* heur, SCIP_CONSHDLR* varbndconshdlr);
static SCIP_RETCODE addQuadraticConstraints(SCIP* scip, SCIP_HEUR* heur, SCIP_CONSHDLR* quadconshdlr);
#ifdef WITH_SOC3
static SCIP_RETCODE addSOCConstraints(SCIP* scip, SCIP_HEUR* heur, SCIP_CONSHDLR* socconshdlr);
static SCIP_RETCODE addSOC3Constraints(SCIP* scip, SCIP_HEUR* heur, SCIP_CONSHDLR* soc3conshdlr);
#endif

/** sets up NLP from constraints in SCIP */
static
SCIP_RETCODE setupNLP(
   SCIP*       scip,  /**< SCIP data structure */
   SCIP_HEUR*  heur   /**< Heuristic */
   )
{
   SCIP_HEURDATA*   heurdata;
   SCIP_Real*       varlb;
   SCIP_Real*       varub;
   SCIP_Real*       objcoeff;
   int*             objvar;
   int              i, cnt;
   int              nconshdlrs;
   SCIP_CONSHDLR**  conshdlrs;

   assert(scip != NULL);
   assert(heur != NULL);
   
   heurdata = SCIPheurGetData(heur);
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
   if (heurdata->nlptimelimit)
     SCIP_CALL( SCIPnlpiSetRealPar(scip, heurdata->nlpi, SCIP_NLPPAR_TILIM, heurdata->nlptimelimit) );
   
   /* add variables to NLP solver; capture variables */
   heurdata->nvars      = SCIPgetNVars(scip);
   heurdata->ndiscrvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &heurdata->var_nlp2scip, SCIPgetVars(scip), heurdata->nvars) );
   SCIP_CALL( SCIPhashmapCreate(&heurdata->var_scip2nlp, SCIPblkmem(scip), heurdata->nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->discrvars, heurdata->ndiscrvars) );
   
   SCIP_CALL( SCIPallocBufferArray(scip, &varlb, heurdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varub, heurdata->nvars) );
   
   cnt = 0; /* counter on discrete variables */
   for (i = 0; i < heurdata->nvars; ++i)
   {
      varlb[i] = SCIPvarGetLbGlobal(heurdata->var_nlp2scip[i]);
      varub[i] = SCIPvarGetUbGlobal(heurdata->var_nlp2scip[i]);
      
      SCIP_CALL( SCIPcaptureVar(scip, heurdata->var_nlp2scip[i]) );

      SCIP_CALL( SCIPhashmapInsert(heurdata->var_scip2nlp, heurdata->var_nlp2scip[i], (void*)(size_t)i) );

      if (SCIPvarGetType(heurdata->var_nlp2scip[i]) <= SCIP_VARTYPE_INTEGER) /* binary or integer */
         heurdata->discrvars[cnt++] = i;
   }
   
   SCIP_CALL( SCIPnlpiAddVars(scip, heurdata->nlpi, heurdata->nvars, varlb, varub, NULL, NULL) );
   
   SCIPfreeBufferArray(scip, &varlb);
   SCIPfreeBufferArray(scip, &varub);
   
   /* add (linear) objective to NLP solver */
   SCIP_CALL( SCIPallocBufferArray(scip, &objcoeff, heurdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &objvar,   heurdata->nvars) );
   cnt = 0;
   for (i = 0; i < heurdata->nvars; ++i)
   {
      if (SCIPvarGetObj(heurdata->var_nlp2scip[i]))
      {
         objcoeff[cnt] = SCIPvarGetObj(heurdata->var_nlp2scip[i]) * SCIPgetObjsense(scip);
         objvar[cnt]   = i;
         ++cnt;
      }
   }
   SCIP_CALL( SCIPnlpiSetObjective(scip, heurdata->nlpi, cnt, objvar, objcoeff, 0, NULL, NULL, NULL, NULL, NULL, 0.) );
   
   SCIPfreeBufferArray(scip, &objcoeff);
   SCIPfreeBufferArray(scip, &objvar);
   
   /* add as many constraints as I can understand */
   
   nconshdlrs = SCIPgetNConshdlrs(scip);
   conshdlrs  = SCIPgetConshdlrs(scip);
   for (i = 0; i < nconshdlrs; ++i)
   {
      if (!SCIPconshdlrGetNConss(conshdlrs[i]))
         continue;
      
      if (strcmp(SCIPconshdlrGetName(conshdlrs[i]), "linear") == 0)
         SCIP_CALL( addLinearConstraints(scip, heur, conshdlrs[i]) );
      else if (strcmp(SCIPconshdlrGetName(conshdlrs[i]), "quadratic") == 0)
         SCIP_CALL( addQuadraticConstraints(scip, heur, conshdlrs[i]) );
      else if (strcmp(SCIPconshdlrGetName(conshdlrs[i]), "varbound") == 0)
         SCIP_CALL( collectVarBoundConstraints(scip, heur, conshdlrs[i]) );
      else if (strcmp(SCIPconshdlrGetName(conshdlrs[i]), "logicor") == 0)
         { SCIPdebugMessage("skip adding logic or constraints to NLP\n"); } /* skip because combinatorial part is fixed in NLP */ 
      else if (strcmp(SCIPconshdlrGetName(conshdlrs[i]), "setppc") == 0)
         { SCIPdebugMessage("skip adding setppc constraints to NLP\n"); }   /* skip because combinatorial part is fixed in NLP */ 
      else if (strcmp(SCIPconshdlrGetName(conshdlrs[i]), "knapsack") == 0)
         { SCIPdebugMessage("skip adding knapsack constraints to NLP\n"); } /* skip because combinatorial part is fixed in NLP */ 
#ifdef WITH_SOC3
      else if (strcmp(SCIPconshdlrGetName(conshdlrs[i]), "soc" ) == 0)
         SCIP_CALL( addSOCConstraints (scip, heur, conshdlrs[i]) );
      else if (strcmp(SCIPconshdlrGetName(conshdlrs[i]), "SOC3") == 0)
         SCIP_CALL( addSOC3Constraints(scip, heur, conshdlrs[i]) );
#endif
      else if (SCIPconshdlrGetNConss(conshdlrs[i]))
         { SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "skip addition of %d constraints of type %s to NLP\n", SCIPconshdlrGetNConss(conshdlrs[i]), SCIPconshdlrGetName(conshdlrs[i])); }
      /* @TODO any other constraints to consider here? */
   }
   
   SCIP_CALL( SCIPnlpStatisticsCreate(scip, &heurdata->nlpstatistics) );

   return SCIP_OKAY;
}

/** adds linear constraints to NLP */
static
SCIP_RETCODE addLinearConstraints(
   SCIP*          scip,          /**< SCIP data structure */ 
   SCIP_HEUR*     heur,          /**< NLP heuristic */
   SCIP_CONSHDLR* linconshdlr    /**< constraint handler for linear constraints */
   )
{
   SCIP_HEURDATA*   heurdata;
   int              i, j, k;
   int              nconss;
   SCIP_CONS**      conss;
   int              nnz;
   SCIP_Real*       lhs;
   SCIP_Real*       rhs;
   int*             rowoffset;
   int*             colindex;
   SCIP_Real*       coeff;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(linconshdlr != NULL);
   
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   
   nconss = SCIPconshdlrGetNConss(linconshdlr);
   conss  = SCIPconshdlrGetConss(linconshdlr);
   
   if (!nconss)
      return SCIP_OKAY;
   
   nnz = 0;
   for (i = 0; i < nconss; ++i)
      nnz += SCIPgetNVarsLinear(scip, conss[i]);
   
   SCIP_CALL( SCIPallocBufferArray(scip, &lhs, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhs, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowoffset, nconss+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colindex, nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coeff, nnz) );
   
   j = 0;
   for (i = 0; i < nconss; ++i)
   {
      lhs[i] = SCIPgetLhsLinear(scip, conss[i]);
      rhs[i] = SCIPgetRhsLinear(scip, conss[i]);
      rowoffset[i] = j;
      memcpy(&coeff[j], SCIPgetValsLinear(scip, conss[i]), SCIPgetNVarsLinear(scip, conss[i]) * sizeof(SCIP_Real));
      for (k = 0; k < SCIPgetNVarsLinear(scip, conss[i]); ++k, ++j)
      {
         assert(SCIPhashmapExists(heurdata->var_scip2nlp, SCIPgetVarsLinear(scip, conss[i])[k]));
         colindex[j] = (int) (size_t) SCIPhashmapGetImage(heurdata->var_scip2nlp, SCIPgetVarsLinear(scip, conss[i])[k]);
      }
   }
   rowoffset[nconss] = nnz;
   assert(j == nnz);
   
   SCIP_CALL( SCIPnlpiAddConstraints(scip, heurdata->nlpi, nconss,
      lhs, rhs,
      rowoffset, colindex, coeff,
      NULL, NULL, NULL, NULL, NULL,
      NULL, NULL
#ifdef WITH_DSL
      , NULL
#endif
      ) );
   
   SCIPfreeBufferArray(scip, &lhs);
   SCIPfreeBufferArray(scip, &rhs);
   SCIPfreeBufferArray(scip, &rowoffset);
   SCIPfreeBufferArray(scip, &colindex);
   SCIPfreeBufferArray(scip, &coeff);
   
   return SCIP_OKAY;
}

/** collect variable bound constraints */
static
SCIP_RETCODE collectVarBoundConstraints(
   SCIP*          scip,           /**< SCIP data structure */
   SCIP_HEUR*     heur,           /**< NLP heuristic */
   SCIP_CONSHDLR* varbndconshdlr  /**< constraint handler for variable bound constraints */
   )
{
   SCIP_HEURDATA*   heurdata;
   SCIP_CONS**      conss;       /* all varbound constraints */
   int              nconss;      /* total number of varbound constraints */
   int              nconss4nlp;  /* number of varbound constraints we have to add to NLP because vbdvar is only implicit integer */
   int              i;
   
   assert(scip != NULL);
   assert(heur != NULL);
   assert(varbndconshdlr != NULL);
   
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->nvarbndconss == 0);
   assert(heurdata->varbndconss == NULL);

   nconss = SCIPconshdlrGetNConss(varbndconshdlr);
   if (!nconss)
      return SCIP_OKAY;
   
   conss = SCIPconshdlrGetConss(varbndconshdlr);
   
   if (heurdata->varbound_explicit)
   {
      /* count for how many constraint the Vbdvar is only implicit integer, so we need to add constraint to NLP */
      nconss4nlp = 0;
      for (i = 0; i < nconss; ++i)
      {
         if (SCIPvarGetType(SCIPgetVbdvarVarbound(scip, conss[i])) > SCIP_VARTYPE_INTEGER)
            nconss4nlp++;
         else
            SCIP_CALL( SCIPcaptureCons(scip, conss[i]) );
      }
   }
   else
   {
      nconss4nlp = nconss;
   }
   
   if (nconss4nlp)
   {
      SCIP_Real*       lhs;
      SCIP_Real*       rhs;
      int*             rowoffset;
      int*             colindex;
      SCIP_Real*       coeff;
      int              j, k;

      heurdata->nvarbndconss = nconss - nconss4nlp;
      if (heurdata->nvarbndconss)
         SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->varbndconss, heurdata->nvarbndconss) );

      SCIP_CALL( SCIPallocBufferArray(scip, &lhs, nconss4nlp) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rhs, nconss4nlp) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rowoffset, nconss4nlp + 1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &colindex, 2 * nconss4nlp) );
      SCIP_CALL( SCIPallocBufferArray(scip, &coeff, 2 * nconss4nlp) );
      
      i = 0; j = 0; k = 0;
      for (i = 0; i < nconss; ++i)
      {
         if (heurdata->varbound_explicit && SCIPvarGetType(SCIPgetVbdvarVarbound(scip, conss[i])) <= SCIP_VARTYPE_INTEGER)
         { /* treat constraint explicitly via boundchanges */
            assert(j < heurdata->nvarbndconss);
            heurdata->varbndconss[j] = conss[i];
            ++j;
            continue;
         }
         
         /* else: add constraint to NLP */
         
         /* varbound constraints: lhs <= x + c * y <= rhs */
         lhs[k]       = SCIPgetLhsVarbound(scip, conss[i]);
         rhs[k]       = SCIPgetRhsVarbound(scip, conss[i]);
         rowoffset[k] = 2*k;
         
         coeff[2*k]    = 1.0;
         assert(SCIPhashmapExists(heurdata->var_scip2nlp, SCIPgetVarVarbound(scip, conss[i])));
         colindex[2*k] = (int) (size_t) SCIPhashmapGetImage(heurdata->var_scip2nlp, SCIPgetVarVarbound(scip, conss[i]));
         
         coeff[2*k+1]    = SCIPgetVbdcoefVarbound(scip, conss[i]);
         assert(SCIPhashmapExists(heurdata->var_scip2nlp, SCIPgetVbdvarVarbound(scip, conss[i])));
         colindex[2*k+1] = (int) (size_t) SCIPhashmapGetImage(heurdata->var_scip2nlp, SCIPgetVbdvarVarbound(scip, conss[i]));

         ++k;
      }
      rowoffset[k] = 2 * nconss4nlp;
      assert(j == heurdata->nvarbndconss);
      assert(k == nconss4nlp);
      
      SCIP_CALL( SCIPnlpiAddConstraints(scip, heurdata->nlpi, nconss4nlp,
         lhs, rhs,
         rowoffset, colindex, coeff,
         NULL, NULL, NULL, NULL, NULL,
         NULL, NULL
#ifdef WITH_DSL
         , NULL
#endif
         ) );
      
      SCIPfreeBufferArray(scip, &lhs);
      SCIPfreeBufferArray(scip, &rhs);
      SCIPfreeBufferArray(scip, &rowoffset);
      SCIPfreeBufferArray(scip, &colindex);
      SCIPfreeBufferArray(scip, &coeff);
   }
   else
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &heurdata->varbndconss, conss, nconss) );
      heurdata->nvarbndconss = nconss;
   }
   
   return SCIP_OKAY;
}

/** for a fixation of discrete variables, applies the variable bound constraints to the NLP */
static
SCIP_RETCODE applyVarBoundConstraints(
   SCIP*       scip,  /**< SCIP data structure */
   SCIP_HEUR*  heur   /**< NLP heuristic */
   )
{
   SCIP_HEURDATA*   heurdata;
   int              i;
   SCIP_VAR*        var;
   SCIP_Real*       varlb;
   SCIP_Real*       varub;
   int*             varidx;
   int              varcnt;
   SCIP_HASHMAP*    varmap;
   SCIP_Real        shift;
   SCIP_CONS*       cons;
   void*            idx;
  
   assert(scip != NULL);
   assert(heur != NULL);
  
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->nvarbndconss == 0 || heurdata->varbndconss != NULL);
   
   if (heurdata->nvarbndconss == 0)
      return SCIP_OKAY;
   
   SCIP_CALL( SCIPallocBufferArray(scip, &varlb,  heurdata->nvarbndconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varub,  heurdata->nvarbndconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varidx, heurdata->nvarbndconss) );
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scip), heurdata->nvarbndconss) );
   varcnt = 0; /* how many variable where we have bounds to change so far */
   
   for (i = 0; i < heurdata->nvarbndconss; ++i)
   {
      cons = heurdata->varbndconss[i];
      var = SCIPgetVarVarbound(scip, cons);
      /* variable bounds should been fixed already, do not release again */
      if (SCIPvarGetType(var) <= SCIP_VARTYPE_INTEGER)
         continue;
      
      idx = SCIPhashmapGetImage(varmap, var);
      if (idx == NULL)
      { /* variable appeared first time */
         assert(SCIPhashmapExists(heurdata->var_scip2nlp, var));
         varidx[varcnt] = (int) (size_t) SCIPhashmapGetImage(heurdata->var_scip2nlp, var);
         varlb[varcnt] = SCIPgetLhsVarbound(scip, cons);
         varub[varcnt] = SCIPgetRhsVarbound(scip, cons);
         
         shift = SCIPgetVbdcoefVarbound(scip, cons) * SCIPgetSolVal(scip, heurdata->startcandidate, SCIPgetVbdvarVarbound(scip, cons));
         if (!SCIPisInfinity(scip, -varlb[varcnt]))
            varlb[varcnt] -= shift;
         if (!SCIPisInfinity(scip,  varub[varcnt]))
            varub[varcnt] -= shift;
         
         if (varlb[varcnt] < SCIPvarGetLbGlobal(var))
            varlb[varcnt] = SCIPvarGetLbGlobal(var);
         if (varub[varcnt] > SCIPvarGetUbGlobal(var))
            varub[varcnt] = SCIPvarGetUbGlobal(var);
         
         SCIP_CALL( SCIPhashmapInsert(varmap, var, (void*)(size_t)(varcnt+1)) );
         
         SCIPdebugMessage("%s: var %s at %d now bounded in [%g, %g] due to %s = %g\n",
            SCIPconsGetName(cons), SCIPvarGetName(var), varidx[varcnt],
            varlb[varcnt], varub[varcnt], SCIPvarGetName(SCIPgetVbdvarVarbound(scip, cons)),
            SCIPgetSolVal(scip, heurdata->startcandidate, SCIPgetVbdvarVarbound(scip, cons)) );
         
         ++varcnt;
      }
      else
      { /* variable appeared before and was stored at position (int)idx - 1 */
         SCIP_Real lhs = SCIPgetLhsVarbound(scip, cons);
         SCIP_Real rhs = SCIPgetRhsVarbound(scip, cons);
         int idx_ = ((int)(size_t)idx) - 1;
         
         assert(idx_ < varcnt);
         
         shift = SCIPgetVbdcoefVarbound(scip, cons) * SCIPgetSolVal(scip, heurdata->startcandidate, SCIPgetVbdvarVarbound(scip, cons));
         if (!SCIPisInfinity(scip, -lhs))
            lhs -= shift;
         if (!SCIPisInfinity(scip,  rhs))
            rhs -= shift;
         
         if (lhs > varlb[idx_])
            varlb[idx_] = lhs;
         if (rhs < varub[idx_])
            varub[idx_] = rhs;
   
         SCIPdebugMessage("%s: var %s at %d now bounded in [%g, %g] due to %s = %g  [updated]\n",
            SCIPconsGetName(cons), SCIPvarGetName(var), varidx[idx_],
            varlb[idx_], varub[idx_], SCIPvarGetName(SCIPgetVbdvarVarbound(scip, cons)),
            SCIPgetSolVal(scip, heurdata->startcandidate, SCIPgetVbdvarVarbound(scip, cons)) );
      }      
   }
   
   SCIP_CALL( SCIPnlpiChgVarBounds(scip, heurdata->nlpi, varcnt, varidx, varlb, varub) );
   
   SCIPfreeBufferArray(scip, &varlb);
   SCIPfreeBufferArray(scip, &varub);
   SCIPfreeBufferArray(scip, &varidx);
   SCIPhashmapFree(&varmap);

   return SCIP_OKAY;
}

/** adds quadratic constraints to NLP */
static
SCIP_RETCODE addQuadraticConstraints(
   SCIP*          scip,         /**< SCIP data structure */
   SCIP_HEUR*     heur,         /**< NLP heuristic */
   SCIP_CONSHDLR* quadconshdlr  /**< constraint handler for quadratic constraints */
   )
{
   SCIP_HEURDATA*   heurdata;
   
   assert(scip != NULL);
   assert(heur != NULL);
   assert(quadconshdlr != NULL);
   
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   
   if (!SCIPconshdlrGetNConss(quadconshdlr))
      return SCIP_OKAY;

   SCIP_CALL( SCIPconsInitnlpiQuadratic(scip, quadconshdlr, heurdata->nlpi, SCIPconshdlrGetNConss(quadconshdlr), SCIPconshdlrGetConss(quadconshdlr), heurdata->var_scip2nlp) );
   
   return SCIP_OKAY;
}

#ifdef WITH_SOC3
static
SCIP_RETCODE addSOCConstraints(
   SCIP*          scip,          /**< SCIP data structure */
   SCIP_HEUR*     heur,          /**< NLP heuristic */
   SCIP_CONSHDLR* socconshdlr    /**< constraint handler for SOC constraints */
   )
{
   SCIP_HEURDATA*   heurdata;
   
   assert(scip != NULL);
   assert(heur != NULL);
   assert(socconshdlr != NULL);
   
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   
   if (!SCIPconshdlrGetNConss(socconshdlr))
      return SCIP_OKAY;

   SCIP_CALL( SCIPconsInitnlpiSOC(scip, socconshdlr, heurdata->nlpi, SCIPconshdlrGetNConss(socconshdlr), SCIPconshdlrGetConss(socconshdlr), heurdata->var_scip2nlp) );
   
   return SCIP_OKAY; 
}

static
SCIP_RETCODE addSOC3Constraints(
   SCIP*          scip,          /**< SCIP data structure */
   SCIP_HEUR*     heur,          /**< NLP heuristic */
   SCIP_CONSHDLR* soc3conshdlr   /**< constraint handler for SOC3 constraints */
   )
{
   SCIP_HEURDATA*   heurdata;
   
   assert(scip != NULL);
   assert(heur != NULL);
   assert(soc3conshdlr != NULL);
   
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   
   if (!SCIPconshdlrGetNConss(soc3conshdlr))
      return SCIP_OKAY;

   SCIP_CALL( SCIPconsInitnlpiSOC3(scip, soc3conshdlr, heurdata->nlpi, SCIPconshdlrGetNConss(soc3conshdlr), SCIPconshdlrGetConss(soc3conshdlr), heurdata->var_scip2nlp) );
   
   return SCIP_OKAY; 
}
#endif

/** frees NLP */
static
SCIP_RETCODE destroyNLP(
   SCIP*        scip,  /**< SCIP data structure */ 
   SCIP_HEUR*   heur   /**< NLP heuristic */
   )
{
   int i;
   SCIP_HEURDATA* heurdata;
   assert(scip != NULL);
   assert(heur != NULL);
  
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->nlpi != NULL);
   assert(heurdata->nlpstatistics != NULL);
   
   for (i = 0; i < heurdata->nvars; ++i)
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
   assert(heurdata->startcandidate == NULL);
   
   SCIPfreeMemoryNull(scip, &heurdata);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
#if 0
static
SCIP_DECL_HEURINIT(heurInitNlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Nlp primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitNlp NULL
#endif


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_HEUREXIT(heurExitNlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Nlp primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitNlp NULL
#endif


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolNlp)
{
   SCIP_Bool havenlp = FALSE; /* first assume all continuous variables are linear */
   
   assert(scip != NULL);
   assert(heur != NULL);
   
   if (SCIPheurGetFreq(heur) < 0)
      return SCIP_OKAY;
      
   /* do not build NLP if there are no nonlinear continuous or impl. integer variables
    */
   if (SCIPgetNContVars(scip) || SCIPgetNImplVars(scip))
   {  /* check if we have nonlinear continuous variables */
      SCIP_CONSHDLR* conshdlr;
      
      conshdlr = SCIPfindConshdlr(scip, "quadratic");
      if (conshdlr && SCIPconshdlrGetNConss(conshdlr))
      {
         SCIP_CONS* cons;
         int        i, j;
         for (i = 0; !havenlp && i < SCIPconshdlrGetNConss(conshdlr); ++i)
         {
            cons = SCIPconshdlrGetConss(conshdlr)[i];
            for (j = 0; !havenlp && j < SCIPgetNQuadVarsQuadratic(cons); ++j)
               if (SCIPvarGetType(SCIPgetQuadVarsQuadratic(cons)[j]) > SCIP_VARTYPE_INTEGER)
                  havenlp = TRUE;
         }
      }
      
#ifdef WITH_SOC3
      if (!havenlp)
      {
         conshdlr = SCIPfindConshdlr(scip, "soc");
         if (conshdlr && SCIPconshdlrGetNConss(conshdlr))
            havenlp = TRUE;  /* @TODO check if some SOC constraint has a continuous variable */
      }
#endif
   }

   if (havenlp)
   {
#ifndef WITH_IPOPT
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "No NLP solver available. NLP local search heuristic will not run.\n");
      return SCIP_OKAY;
#endif
      setupNLP(scip, heur);
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
  
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   
   if (heurdata->nlpi)
      SCIP_CALL( destroyNLP(scip, heur) );
   
   if (heurdata->nvarbndconss)
   {
      int i;
      assert(heurdata->varbndconss != NULL);
      for (i = 0; i < heurdata->nvarbndconss; ++i)
      {
         assert(heurdata->varbndconss[i] != NULL);
         SCIP_CALL( SCIPreleaseCons(scip, &heurdata->varbndconss[i]) );
      }
      SCIPfreeMemoryArray(scip, &heurdata->varbndconss);
      heurdata->nvarbndconss = 0;
   }
   
   if (heurdata->startcandidate)
      SCIP_CALL( SCIPfreeSol(scip, &heurdata->startcandidate) );

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
  
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   
   *result = SCIP_DIDNOTRUN;
   
   /* probably do not have continuous or implicit integer variables */
   if (!heurdata->nlpi)
      return SCIP_OKAY;
   
   if (!heurdata->startcandidate)
   {
      int nlpcands;
      
      /* only call heuristic if optimal LP solution is available */
      if (SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL)
      {
         SCIPdebugMessage("skip NLP heuristic because no start candidate given and no LP solution available\n");
         return SCIP_OKAY;
      }

      /* only call heuristic, if there are no fractional variables */
      SCIP_CALL( SCIPgetLPBranchCands(scip, NULL, NULL, NULL, &nlpcands, NULL) );
      if (nlpcands)
      {
         SCIPdebugMessage("skip NLP heuristic because no start candidate given and current LP solution is fractional\n");
         return SCIP_OKAY;
      }
   }

   itercontingent  = (SCIP_Longint)(heurdata->iterquot * SCIPgetNNodes(scip));
   /* weight by previous success of heuristic */
   itercontingent  = (SCIP_Longint)(itercontingent * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   itercontingent += heurdata->iteroffset; 
   itercontingent -= heurdata->iterused;
   if (itercontingent < heurdata->itermin)
   {
      SCIPdebugMessage("skip NLP heuristic; contingent=%"SCIP_LONGINT_FORMAT"; minimal number of iterations=%d; success ratio=%g\n", itercontingent, heurdata->itermin, (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
      return SCIP_OKAY;
   }
   if (heurdata->nlpiterlimit)
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
      if (heurdata->nlptimelimit)
         timelimit = MIN(heurdata->nlptimelimit, timelimit);
      SCIP_CALL( SCIPnlpiSetRealPar(scip, heurdata->nlpi, SCIP_NLPPAR_TILIM, timelimit) );
   }
   
   SCIP_CALL( SCIPallocBufferArray(scip, &startpoint, heurdata->nvars) );

   SCIP_CALL( SCIPgetSolVals(scip, heurdata->startcandidate, heurdata->nvars, heurdata->var_nlp2scip, startpoint) );
   SCIP_CALL( SCIPnlpiSetInitialGuess(scip, heurdata->nlpi, startpoint) );

   if (heurdata->ndiscrvars)
   { /* fix discrete variables */
      int i;
      SCIP_Real* discrfix;
      
      SCIP_CALL( SCIPallocBufferArray(scip, &discrfix, heurdata->ndiscrvars) );
      for (i = 0; i < heurdata->ndiscrvars; ++i)
      {
         discrfix[i] = startpoint[heurdata->discrvars[i]];
         if (!SCIPisIntegral(scip, discrfix[i]))
            discrfix[i] = SCIPfrac(scip, discrfix[i]) > .5 ? SCIPceil(scip, discrfix[i]) : SCIPfloor(scip, discrfix[i]);
         if (discrfix[i] < SCIPvarGetLbGlobal(heurdata->var_nlp2scip[heurdata->discrvars[i]]))
            discrfix[i] = SCIPvarGetLbGlobal(heurdata->var_nlp2scip[heurdata->discrvars[i]]);
         else if (discrfix[i] > SCIPvarGetUbGlobal(heurdata->var_nlp2scip[heurdata->discrvars[i]]))
            discrfix[i] = SCIPvarGetUbGlobal(heurdata->var_nlp2scip[heurdata->discrvars[i]]);
      }
      SCIP_CALL( SCIPnlpiChgVarBounds(scip, heurdata->nlpi, heurdata->ndiscrvars, heurdata->discrvars, discrfix, discrfix) );
      SCIPfreeBufferArray(scip, &discrfix);
   }
   SCIP_CALL( applyVarBoundConstraints(scip, heur) );
   
   SCIPfreeBufferArray(scip, &startpoint);   
   if (heurdata->startcandidate)
      SCIPfreeSol(scip, &heurdata->startcandidate);
      
   SCIP_CALL( SCIPnlpiSetIntPar(scip, heurdata->nlpi, SCIP_NLPPAR_ITLIM, itercontingent) );
   SCIPdebugMessage("start NLP solve with iteration limit %"SCIP_LONGINT_FORMAT" and timelimit %g\n", itercontingent, timelimit);
   
   SCIP_CALL( SCIPnlpiSolve(scip, heurdata->nlpi) );

   SCIPdebugMessage("NLP solver returned with termination status %d and solution status %d\n", SCIPnlpiGetTermstat(scip, heurdata->nlpi), SCIPnlpiGetSolstat(scip, heurdata->nlpi));
   
   if (SCIPnlpiGetTermstat(scip, heurdata->nlpi) >= SCIP_NLPITERMSTAT_MEMERR)
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "NLP solver returned with bad termination status %d. Will not run NLP heuristic again for this run.\n",  SCIPnlpiGetTermstat(scip, heurdata->nlpi));
      SCIP_CALL( destroyNLP(scip, heur) );
      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPnlpiGetStatistics(scip, heurdata->nlpi, heurdata->nlpstatistics) );

   heurdata->iterused += SCIPnlpStatisticsGetNIterations(heurdata->nlpstatistics);
   SCIPdebugMessage("NLP solver used %d iterations and %g seconds\n", SCIPnlpStatisticsGetNIterations(heurdata->nlpstatistics), SCIPnlpStatisticsGetTotalTime(heurdata->nlpstatistics));

   if (SCIPnlpiGetSolstat(scip, heurdata->nlpi) <= SCIP_NLPSOLSTAT_FEASIBLE)
   {
      SCIP_Real* primals = NULL;
      SCIP_SOL*  sol;
      SCIP_Bool  stored; 
      SCIP_Bool  feasible;

      SCIP_CALL( SCIPnlpiGetSolution(scip, heurdata->nlpi, &primals) );
      assert(primals);

      SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
      SCIP_CALL( SCIPsetSolVals(scip, sol, heurdata->nvars, heurdata->var_nlp2scip, primals) );

/*#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPcheckSolOrig(scip, sol, &feasible, TRUE, TRUE) );
#endif*/

      SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, TRUE, TRUE, &feasible) );
      if (feasible)
      { /* add solution */ 
         SCIP_CALL( SCIPaddSolFree(scip, &sol, &stored) );
         if (stored)
         {
            SCIPdebugMessage("SCIP stored solution\n");
            *result = SCIP_FOUNDSOL;
         }
      }
      else if (heurdata->resolvetolfactor < 1.0)
      { /* resolve with tighter tolerances */
         SCIPdebugMessage("solution reported by NLP solver not feasible for SCIP, resolve with feasibility tolerance %g\n", heurdata->resolvetolfactor*SCIPfeastol(scip));
         SCIP_CALL( SCIPnlpiSetRealPar(scip, heurdata->nlpi, SCIP_NLPPAR_FEASTOL, heurdata->resolvetolfactor*SCIPfeastol(scip)) );
         SCIP_CALL( SCIPnlpiSetInitialGuess(scip, heurdata->nlpi, primals) );
         SCIP_CALL( SCIPnlpiSolve(scip, heurdata->nlpi) );
         SCIP_CALL( SCIPnlpiGetStatistics(scip, heurdata->nlpi, heurdata->nlpstatistics) );
         heurdata->iterused += SCIPnlpStatisticsGetNIterations(heurdata->nlpstatistics);
         if (SCIPnlpiGetSolstat(scip, heurdata->nlpi) <= SCIP_NLPSOLSTAT_FEASIBLE)
         { /* still feasible, hope that SCIP accepts it now */
            SCIP_CALL( SCIPnlpiGetSolution(scip, heurdata->nlpi, &primals) );
            assert(primals);
            SCIP_CALL( SCIPsetSolVals(scip, sol, heurdata->nvars, heurdata->var_nlp2scip, primals) );
#ifdef SCIP_DEBUG
            SCIP_CALL( SCIPcheckSolOrig(scip, sol, &feasible, TRUE, TRUE) );
#endif
            SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, TRUE, &stored) );
            if (stored)
            {
               SCIPdebugMessage("SCIP stored solution\n");
               *result = SCIP_FOUNDSOL;
            }
            else
            {
               SCIPdebugMessage("NLP solver said feasible after tigthening tolerances, but SCIP still did not store solution\n");
            }
         }
         else
            SCIP_CALL( SCIPfreeSol(scip, &sol) );
         
         /* reset to original tolerance */
         SCIP_CALL( SCIPnlpiSetRealPar(scip, heurdata->nlpi, SCIP_NLPPAR_FEASTOL, SCIPfeastol(scip)) );
      }
      else
         SCIP_CALL( SCIPfreeSol(scip, &sol) );
   }
   else
      *result = SCIP_DIDNOTFIND;
   
   /* TODO reset time and iterlimit in nlp solver? */
   
   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the Nlp primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurNlp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create Nlp primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   memset(heurdata, 0, sizeof(SCIP_HEURDATA));
   
   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING,
         heurFreeNlp, heurInitNlp, heurExitNlp, 
         heurInitsolNlp, heurExitsolNlp, heurExecNlp,
         heurdata) );

   /* add Nlp primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/"HEUR_NAME"/nlpverblevel",      "verbosity level of NLP solver",                                                  &heurdata->nlpverblevel, FALSE, 0,   0, 2,                  NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/"HEUR_NAME"/nlpiterlimit",      "iteration limit of NLP solver; 0 to use solver default",                         &heurdata->nlpiterlimit, FALSE, 0,   0, INT_MAX,            NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/nlptimelimit",      "time limit of NLP solver; 0 to use solver default",                              &heurdata->nlptimelimit, FALSE, 0,   0, SCIPinfinity(scip), NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/resolvetolfactor",  "if SCIP does not accept a solution which the NLP solver thinks is feasible, the feasibility tolerance is reduced by this factor and the NLP resolved (set to 1. to turn off resolve", &heurdata->resolvetolfactor, TRUE, 0.005, 0, 1., NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/"HEUR_NAME"/iteroffset",        "number of iterations added to the contingent of the total number of iterations", &heurdata->iteroffset,   FALSE, 500, 0, INT_MAX,            NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/iterquotient",      "contingent of NLP iterations in relation to the number of nodes in SCIP",        &heurdata->iterquot,     FALSE, 0.1, 0, SCIPinfinity(scip), NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/"HEUR_NAME"/itermin",           "contingent of NLP iterations in relation to the number of nodes in SCIP",        &heurdata->itermin,      FALSE, 300, 0, INT_MAX,            NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/varbound_explicit", "whether variable bound constraints should be handled explicitly before solving NLP instead of adding them to the NLP", &heurdata->varbound_explicit, FALSE, TRUE, NULL, NULL) );

   return SCIP_OKAY;
}

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
   
   if (SCIPgetStage(scip) >= SCIP_STAGE_SOLVED)
      return SCIP_OKAY;
   
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   
   /* probably have no continuous variables */
   if (!heurdata->nlpi)
      return SCIP_OKAY;
   
   takenew = !heurdata->startcandidate;
   if (!heurdata->startcandidate ||
      SCIPisGT(scip, heurdata->startcand_violation, violation) ||
      SCIPisRelGT(scip, SCIPgetSolTransObj(scip, heurdata->startcandidate)*SCIPgetObjsense(scip), SCIPgetSolTransObj(scip, solcand)*SCIPgetObjsense(scip)) )
   {
      if (heurdata->startcandidate)
         SCIP_CALL( SCIPfreeSol(scip, &heurdata->startcandidate) );
      SCIP_CALL( SCIPcreateSolCopy(scip, &heurdata->startcandidate, solcand) );
      SCIP_CALL( SCIPunlinkSol(scip, heurdata->startcandidate) );
      heurdata->startcand_violation = violation;
   }

   return SCIP_OKAY;
}

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

/**@file    heur_nlp.c
 * @ingroup PRIMALHEURISTICS
 * @brief   NLP local search primal heuristic
 * @author  Stefan Vigerske
 * 
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
#include "scip/cons_soc.h"
#ifdef WITH_SOC3
#include "cons_soc3.h"
#endif
#ifdef WITH_NL
#include "cons_nl_C.h"
#include "expression.h"
#endif
#ifdef WITH_SIGNPOWER
#include "cons_signpower.h"
#endif
#ifdef WITH_UNIVARDEFINITE
#include "cons_univardefinite.h"
#endif

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
   SCIP_Bool             triedsetupnlp;      /**< whether we have tried to setup an NLP */ 
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for global bound change events */
   
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
                         
   SCIP_Bool             varboundexplicit;   /**< whether variable bound constraints should be handled explicitly before solving an NLP instead of adding them as linear constraints to to the NLP */
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
   int*             nlininds;
   int**            lininds;
   SCIP_Real**      linvals;

   int              nconss;
   int              nuseconss;
   int              i;
   int              k;

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(linconshdlr != NULL);
   
   nconss = SCIPconshdlrGetNConss(linconshdlr);
   conss  = SCIPconshdlrGetConss(linconshdlr);
   
   if( !nconss )
      return SCIP_OKAY;
   
   SCIP_CALL( SCIPallocBufferArray(scip, &useconss, nconss) );
   
   /* count the number of constraints to add to NLP
    * constraints with only integer variables are not added to the NLP, since all its variables will be fixed later */
   nuseconss = 0;
   for( i = 0; i < nconss; ++i )
   {
      if( SCIPconsIsLocal(conss[i]) )
         continue;
      for( k = 0; k < SCIPgetNVarsLinear(scip, conss[i]); ++k )
      {
         if( SCIPvarGetType(SCIPgetVarsLinear(scip, conss[i])[k]) > SCIP_VARTYPE_INTEGER )
         {
            useconss[nuseconss] = conss[i];
            ++nuseconss;
            break;
         }
      }
   }
   
   if( !nuseconss )
   {  /* no linear constraints to add to NLP, so we are done */
      SCIPfreeBufferArray(scip, &useconss);
      return SCIP_OKAY;
   }
   
   SCIP_CALL( SCIPallocBufferArray(scip, &lhs, nuseconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhs, nuseconss) );
   /* arrays to store linear constraints */ 
   SCIP_CALL( SCIPallocBufferArray(scip, &nlininds, nuseconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lininds,  nuseconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals,  nuseconss) );
   
   for( i = 0; i < nuseconss; ++i )
   {
      SCIP_VAR** vars;
      SCIP_Real* vals;
      int nvars;

      vars  = SCIPgetVarsLinear(scip, useconss[i]);
      vals  = SCIPgetValsLinear(scip, useconss[i]);
      nvars = SCIPgetNVarsLinear(scip, useconss[i]);
      assert(vars != NULL);
      assert(vals != NULL);

      lhs[i] = SCIPgetLhsLinear(scip, useconss[i]);
      rhs[i] = SCIPgetRhsLinear(scip, useconss[i]);
      
      SCIP_CALL( SCIPallocBufferArray(scip, &lininds[i], nvars) );
      nlininds[i] = nvars;
      linvals[i]  = vals;
     
      /* fill column indices of variables (vars) in conss[i] with variable indices in NLP */
      for( k = 0; k < nvars; ++k )
      {
         assert(SCIPhashmapExists(heurdata->var_scip2nlp, vars[k]));
         lininds[i][k] = (int) (size_t) SCIPhashmapGetImage(heurdata->var_scip2nlp, vars[k]);
      }
   }
   
   SCIP_CALL( SCIPnlpiAddConstraints(scip, heurdata->nlpi, nuseconss,
         lhs, rhs, nlininds, lininds, linvals,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );
   
   for( i = 0; i < nuseconss; ++i )
   {
      SCIPfreeBufferArray(scip, &lininds[i]);
   }
   
   SCIPfreeBufferArray(scip, &lininds);
   SCIPfreeBufferArray(scip, &linvals);
   SCIPfreeBufferArray(scip, &nlininds);
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
   int                   nconsslocal;        /* number of local varbound constraints */
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
   
   /* count for how many global constraint the Vbdvar is only implicit integer, so we need to add constraint to NLP */
   nconss4nlp = 0;
   nconsslocal = 0;
   for( i = 0; i < nconss; ++i )
   {
      if( SCIPconsIsLocal(conss[i]) )
      {
         ++nconsslocal;
      }
      else if( SCIPvarGetType(SCIPgetVbdvarVarbound(scip, conss[i])) > SCIP_VARTYPE_INTEGER )
      {  /* will add constraint to NLP soon */
         nconss4nlp++;
      }
      else
      {  /* will handle constraint explicitely before solving NLP */
         SCIP_CALL( SCIPcaptureCons(scip, conss[i]) );
      }
   }

   if( nconsslocal == 0 && nconss4nlp != 0 )
   {
      SCIP_Real*       lhs;
      SCIP_Real*       rhs;
      int*             nlininds;
      int**            lininds;
      SCIP_Real**      linvals;
      int              j;
      int              k;

      /* number of and space for varbound constraints that are handled explicitely */
      heurdata->nvarbndconss = nconss - nconss4nlp - nconsslocal;
      if( heurdata->nvarbndconss )
         SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->varbndconss, heurdata->nvarbndconss) );

      /* memory to store those varbound constraints that are added to the NLP as linear constraint */
      SCIP_CALL( SCIPallocBufferArray(scip, &lhs,      nconss4nlp) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rhs,      nconss4nlp) );
      SCIP_CALL( SCIPallocBufferArray(scip, &nlininds, nconss4nlp) );
      SCIP_CALL( SCIPallocBufferArray(scip, &lininds,  nconss4nlp) );
      SCIP_CALL( SCIPallocBufferArray(scip, &linvals,  nconss4nlp) );
      
      j = 0; /* the number of explicitely handled varbound constraints passed so far */
      k = 0; /* the number of varbound constraints for the NLP passed so far */
      
      for( i = 0; i < nconss; ++i )
      {
         /* skip local constraints */
         if( SCIPconsIsLocal(conss[i]) )
            continue;
         
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
         
         nlininds[k] = 2;
         SCIP_CALL( SCIPallocBufferArray(scip, &lininds[k], 2) );
         SCIP_CALL( SCIPallocBufferArray(scip, &linvals[k], 2) );
         
         /* add x term to coefficient matrix */
         linvals[k][0] = 1.0;
         assert(SCIPhashmapExists(heurdata->var_scip2nlp, SCIPgetVarVarbound(scip, conss[i])));
         lininds[k][0] = (int) (size_t) SCIPhashmapGetImage(heurdata->var_scip2nlp, SCIPgetVarVarbound(scip, conss[i]));
         
         /* add c * y term to coefficient matrix */
         linvals[k][1] = SCIPgetVbdcoefVarbound(scip, conss[i]);
         assert(SCIPhashmapExists(heurdata->var_scip2nlp, SCIPgetVbdvarVarbound(scip, conss[i])));
         lininds[k][1] = (int) (size_t) SCIPhashmapGetImage(heurdata->var_scip2nlp, SCIPgetVbdvarVarbound(scip, conss[i]));

         ++k;
      }

      assert(j == heurdata->nvarbndconss);
      assert(k == nconss4nlp);
      
      SCIP_CALL( SCIPnlpiAddConstraints(scip, heurdata->nlpi, nconss4nlp,
            lhs, rhs, nlininds, lininds, linvals,
            NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );
      
      for( k = 0; k < nconss4nlp; ++k )
      {
         SCIPfreeBufferArray(scip, &lininds[k]);
         SCIPfreeBufferArray(scip, &linvals[k]);
      }
      
      SCIPfreeBufferArray(scip, &nlininds);
      SCIPfreeBufferArray(scip, &lininds);
      SCIPfreeBufferArray(scip, &linvals);
      SCIPfreeBufferArray(scip, &rhs);
      SCIPfreeBufferArray(scip, &lhs);
   }
   else
   {  /* all varbound constraints are global and handled explicitely */
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

   /* let the constraint handler for quadratic constraints fill the NLP with its constraints */
   SCIP_CALL( SCIPconsInitNlpiQuadratic(scip, quadconshdlr, heurdata->nlpi, SCIPconshdlrGetNConss(quadconshdlr), 
         SCIPconshdlrGetConss(quadconshdlr), heurdata->var_scip2nlp) );
   
   return SCIP_OKAY;
}

/** adds second order cone constraints to NLP */
static
SCIP_RETCODE addSOCConstraints(
   SCIP*          scip,          /**< SCIP data structure */
   SCIP_HEURDATA* heurdata,      /**< heuristic data structure */
   SCIP_CONSHDLR* socconshdlr    /**< constraint handler for SOC constraints */
   )
{
   assert(scip        != NULL);
   assert(heurdata    != NULL);
   assert(socconshdlr != NULL);
   
   if (!SCIPconshdlrGetNConss(socconshdlr))
      return SCIP_OKAY;

   SCIP_CALL( SCIPconsInitNlpiSOC(scip, socconshdlr, heurdata->nlpi, SCIPconshdlrGetNConss(socconshdlr), SCIPconshdlrGetConss(socconshdlr), heurdata->var_scip2nlp) );
   
   return SCIP_OKAY; 
}

#ifdef WITH_SOC3
/** adds three dimensional second order cone constraints to NLP */
static
SCIP_RETCODE addSOC3Constraints(
   SCIP*          scip,          /**< SCIP data structure */
   SCIP_HEURDATA* heurdata,      /**< heuristic data structure */
   SCIP_CONSHDLR* soc3conshdlr   /**< constraint handler for SOC3 constraints */
   )
{
   assert(scip         != NULL);
   assert(heurdata     != NULL);
   assert(soc3conshdlr != NULL);
   
   if (!SCIPconshdlrGetNConss(soc3conshdlr))
      return SCIP_OKAY;

   SCIP_CALL( SCIPconsInitnlpiSOC3(scip, soc3conshdlr, heurdata->nlpi, SCIPconshdlrGetNConss(soc3conshdlr), SCIPconshdlrGetConss(soc3conshdlr), heurdata->var_scip2nlp) );
   
   return SCIP_OKAY; 
}
#endif

#ifdef WITH_NL
/** adds nonlinear constraints to NLP */
static
SCIP_RETCODE addNonlinearConstraints(
   SCIP*          scip,          /**< SCIP data structure */
   SCIP_HEURDATA* heurdata,      /**< heuristic data structure */
   SCIP_CONSHDLR* nlconshdlr     /**< constraint handler for nonlinear constraints */
   )
{
   assert(scip       != NULL);
   assert(heurdata   != NULL);
   assert(nlconshdlr != NULL);

   if (!SCIPconshdlrGetNConss(nlconshdlr))
      return SCIP_OKAY;

   SCIP_CALL( SCIPconsInitnlpiNonlinear(scip, nlconshdlr, heurdata->nlpi, SCIPconshdlrGetNConss(nlconshdlr), SCIPconshdlrGetConss(nlconshdlr), heurdata->var_scip2nlp) );

   return SCIP_OKAY;
}
#endif

#ifdef WITH_SIGNPOWER
/** adds nonlinear constraints to NLP */
static
SCIP_RETCODE addSignpowerConstraints(
   SCIP*          scip,          /**< SCIP data structure */
   SCIP_HEURDATA* heurdata,      /**< heuristic data structure */
   SCIP_CONSHDLR* sgnpowconshdlr /**< constraint handler for signpower constraints */
   )
{
   assert(scip           != NULL);
   assert(heurdata       != NULL);
   assert(sgnpowconshdlr != NULL);

   if (!SCIPconshdlrGetNConss(sgnpowconshdlr))
      return SCIP_OKAY;

   SCIP_CALL( SCIPconsInitNlpiSignpower(scip, sgnpowconshdlr, heurdata->nlpi, SCIPconshdlrGetNConss(sgnpowconshdlr), SCIPconshdlrGetConss(sgnpowconshdlr), heurdata->var_scip2nlp) );

   return SCIP_OKAY;
}
#endif

#ifdef WITH_UNIVARDEFINITE
/** adds univariate definite constraints to NLP */
static
SCIP_RETCODE addUnivardefiniteConstraints(
   SCIP*          scip,          /**< SCIP data structure */
   SCIP_HEURDATA* heurdata,      /**< heuristic data structure */
   SCIP_CONSHDLR* univardefiniteconshdlr    /**< constraint handler for univariate definite constraints */
   )
{
   assert(scip        != NULL);
   assert(heurdata    != NULL);
   assert(univardefiniteconshdlr != NULL);

   if (!SCIPconshdlrGetNConss(univardefiniteconshdlr))
      return SCIP_OKAY;

   SCIP_CALL( SCIPconsInitNlpiUnivardefinite(scip, univardefiniteconshdlr, heurdata->nlpi, SCIPconshdlrGetNConss(univardefiniteconshdlr), SCIPconshdlrGetConss(univardefiniteconshdlr), heurdata->var_scip2nlp) );

   return SCIP_OKAY;
}
#endif

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

   /* create an instance of our NLP solver; so far we have only IPOPT */
#ifdef WITH_IPOPT
   SCIP_CALL( SCIPcreateNlpSolverIpopt(scip, &heurdata->nlpi) );
#else
   SCIPerrorMessage("No NLP solver available. Cannot setup NLP.\n");
   return SCIP_ERROR;
#endif
   SCIP_CALL( SCIPnlpiInit(scip, heurdata->nlpi, "nlp") );
   
   /* set some parameters of NLP solver */ 
   SCIP_CALL( SCIPnlpiSetIntPar(scip, heurdata->nlpi, SCIP_NLPPAR_VERBLEVEL, heurdata->nlpverblevel) );
   if( heurdata->nlptimelimit )
   {
      SCIP_CALL( SCIPnlpiSetRealPar(scip, heurdata->nlpi, SCIP_NLPPAR_TILIM, heurdata->nlptimelimit) );
   }

   /* collect and capture variables, collect discrete variables, collect bounds
    * assign variable indices to variables */
   heurdata->nvars = SCIPgetNVars(scip);
   heurdata->ndiscrvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &heurdata->var_nlp2scip, SCIPgetVars(scip), heurdata->nvars) );
   SCIP_CALL( SCIPhashmapCreate(&heurdata->var_scip2nlp, SCIPblkmem(scip), heurdata->nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->discrvars, heurdata->ndiscrvars) );
   
   SCIP_CALL( SCIPallocBufferArray(scip, &varlb, heurdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varub, heurdata->nvars) );
   
   cnt = 0; /* counts number of discrete variables passed so far */
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
   
   /* add variables to NLP solver */
   SCIP_CALL( SCIPnlpiAddVars(scip, heurdata->nlpi, heurdata->nvars, varlb, varub, NULL, NULL) );
   
   SCIPfreeBufferArray(scip, &varub);
   SCIPfreeBufferArray(scip, &varlb);

   /* collect objective coefficients for minimization objective */
   SCIP_CALL( SCIPallocBufferArray(scip, &objcoeff, heurdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &objvar, heurdata->nvars) );

   cnt = 0; /* number of nonzeros in objective passed so far */
   for( i = 0; i < heurdata->nvars; ++i )
   {
      if( SCIPvarGetObj(heurdata->var_nlp2scip[i]) )
      {
         /* NLPI understands only minimization problems, so we turn maximization problems into minimization problems */ 
         objcoeff[cnt] = SCIPvarGetObj(heurdata->var_nlp2scip[i]) * (int)SCIPgetObjsense(scip);
         objvar[cnt]   = i;
         ++cnt;
      }
   }
   /* add (linear) objective to NLP solver */
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
         /* skip combinatorial constraints, since all their variables will be fixed when the NLP is solved */ 
         SCIPdebugMessage("skip adding logicor, setppc, and knapsack constraints to NLP\n"); 
      }
      else if (strcmp(SCIPconshdlrGetName(conshdlrs[i]), "soc" ) == 0)
      {
         SCIP_CALL( addSOCConstraints (scip, heurdata, conshdlrs[i]) );
      }
#ifdef WITH_SOC3
      else if (strcmp(SCIPconshdlrGetName(conshdlrs[i]), "SOC3") == 0)
      {
         SCIP_CALL( addSOC3Constraints(scip, heurdata, conshdlrs[i]) );
      }
#endif
#ifdef WITH_NL
      else if (strcmp(SCIPconshdlrGetName(conshdlrs[i]), "nonlinear" ) == 0)
      {
         SCIP_CALL( addNonlinearConstraints (scip, heurdata, conshdlrs[i]) );
      }
#endif
#ifdef WITH_SIGNPOWER
      else if (strcmp(SCIPconshdlrGetName(conshdlrs[i]), "signpower" ) == 0)
      {
         SCIP_CALL( addSignpowerConstraints (scip, heurdata, conshdlrs[i]) );
      }
#endif
#ifdef WITH_UNIVARDEFINITE
      else if (strcmp(SCIPconshdlrGetName(conshdlrs[i]), "univardefinite" ) == 0)
      {
         SCIP_CALL( addUnivardefiniteConstraints (scip, heurdata, conshdlrs[i]) );
      }
#endif
      else
      {
         /* @TODO any other constraints to consider here? */
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "skip addition of %d constraints of type %s to NLP\n", 
            SCIPconshdlrGetNConss(conshdlrs[i]), SCIPconshdlrGetName(conshdlrs[i])); 
      }
   }
   
   /* initialize data structure for NLP solve statistics */
   SCIP_CALL( SCIPnlpStatisticsCreate(scip, &heurdata->nlpstatistics) );
   
   /** catch variable global bounds change events */
   for( i = 0; i < heurdata->nvars; ++i )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, heurdata->var_nlp2scip[i], SCIP_EVENTTYPE_GBDCHANGED, heurdata->eventhdlr, (SCIP_EVENTDATA*)heurdata, NULL) );
   }

   return SCIP_OKAY;
}
#endif

/** process variable global bound change event */
static
SCIP_DECL_EVENTEXEC(processVarEvent)
{
   SCIP_HEURDATA* heurdata;
   SCIP_VAR*      var;
   int            nlpidx;
   SCIP_Real      lb, ub;
   
   assert(scip      != NULL);
   assert(event     != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);

   heurdata = (SCIP_HEURDATA*)eventdata;
   assert(heurdata  != NULL);

   var = SCIPeventGetVar(event);
   assert(var != NULL);
   
   assert(SCIPhashmapExists(heurdata->var_scip2nlp, var));
   nlpidx = (int) (size_t) SCIPhashmapGetImage(heurdata->var_scip2nlp, var);

   lb = SCIPvarGetLbGlobal(var);
   ub = SCIPvarGetUbGlobal(var);
   SCIP_CALL( SCIPnlpiChgVarBounds(scip, heurdata->nlpi, 1, &nlpidx, &lb, &ub) );

   return SCIP_OKAY;
}

/** Checks if an NLP should be setup and if positive, then sets up the NLP.
 * 
 * Looks at the current problem and if it has continuous variables in nonlinear constraints.
 * If so and if there is an NLP solver available, then sets up the NLP.
 */
static
SCIP_RETCODE checkCIPandSetupNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur                /**< NLP heuristic */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_Bool      havenlp;

   assert(scip != NULL);
   assert(heur != NULL);
   
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   
   if( heurdata->nlpi != NULL )
   {
      assert(heurdata->triedsetupnlp == TRUE);
      SCIPdebugMessage("already have NLP; skip additional setup.\n");
      return SCIP_OKAY;
   }
   
   heurdata->triedsetupnlp = TRUE;
   
   assert(SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED);
   assert(SCIPgetStage(scip) < SCIP_STAGE_SOLVED);
   
   havenlp = FALSE; /* first assume all continuous variables are linear */

   /* do not build NLP if there are no nonlinear continuous or impl. integer variables */
   if( SCIPgetNContVars(scip) > 0 || SCIPgetNImplVars(scip) > 0 )
   {
      /* check if we have nonlinear continuous or impl. integer variables in the quadratic constraints */
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
            
            if( SCIPconsIsLocal(cons) )
               continue;

            nquadvars = SCIPgetNQuadVarsQuadratic(scip, cons);
            quadvars = SCIPgetQuadVarsQuadratic(scip, cons);
            assert(nquadvars == 0 || quadvars != NULL);

            for( j = 0; !havenlp && j < nquadvars; ++j )
               if( SCIPvarGetType(quadvars[j]) == SCIP_VARTYPE_IMPLINT || SCIPvarGetType(quadvars[j]) == SCIP_VARTYPE_CONTINUOUS )
                  havenlp = TRUE;
         }
      }

      if (!havenlp)
      {
         conshdlr = SCIPfindConshdlr(scip, "soc");
         if (conshdlr && SCIPconshdlrGetNConss(conshdlr))
            havenlp = TRUE;  /* @TODO check if some SOC constraint has a continuous variable */
      }

#ifdef WITH_NL
      if (!havenlp)
      {
         conshdlr = SCIPfindConshdlr(scip, "nonlinear");
         if (conshdlr && SCIPconshdlrGetNConss(conshdlr))
         {
            SCIP_CONS*     cons;
            SCIP_EXPRTREE* exprtree;
            int            i, j;
            for (i = 0; !havenlp && i < SCIPconshdlrGetNConss(conshdlr); ++i)
            {
               cons = SCIPconshdlrGetConss(conshdlr)[i];
               
               if( SCIPconsIsLocal(cons) )
                  continue;

               exprtree = SCIPgetExprtreeNonlinear(cons);
               if (!exprtree)
                  continue;
               for (j = 0; !havenlp && j < SCIPexprtreeGetNVars(exprtree); ++j)
                  if (SCIPvarGetType(SCIPexprtreeGetVars(exprtree)[j]) > SCIP_VARTYPE_IMPLINT)
                     havenlp = TRUE;
            }
         }
      }
#endif

#ifdef WITH_SIGNPOWER
      if (!havenlp)
      {
         conshdlr = SCIPfindConshdlr(scip, "signpower");
         if (conshdlr && SCIPconshdlrGetNConss(conshdlr))
         {
            SCIP_CONS*     cons;
            SCIP_VAR*      x;
            int            i;
            for( i = 0; !havenlp && i < SCIPconshdlrGetNConss(conshdlr); ++i )
            {
               cons     = SCIPconshdlrGetConss(conshdlr)[i];
               
               if( SCIPconsIsLocal(cons) )
                  continue;

               x        = SCIPgetNonlinearVarSignpower(scip, cons);
               if( x == NULL )
                  continue;
               if( SCIPvarGetType(x) > SCIP_VARTYPE_IMPLINT )
                  havenlp = TRUE;
            }
         }
      }
#endif

#ifdef WITH_UNIVARDEFINITE
      if (!havenlp)
      {
         conshdlr = SCIPfindConshdlr(scip, "univardefinite");
         if (conshdlr && SCIPconshdlrGetNConss(conshdlr))
         {
            SCIP_CONS*     cons;
            SCIP_VAR*      x;
            int            i;
            for( i = 0; !havenlp && i < SCIPconshdlrGetNConss(conshdlr); ++i )
            {
               cons     = SCIPconshdlrGetConss(conshdlr)[i];
               
               if( SCIPconsIsLocal(cons) )
                  continue;

               x        = SCIPgetNonlinearVarUnivardefinite(scip, cons);
               if( x == NULL )
                  continue;
               if( SCIPvarGetType(x) > SCIP_VARTYPE_IMPLINT )
                  havenlp = TRUE;
            }
         }
      }
#endif
   }
   
   if( !havenlp )
   { /* no suitable nonlinearities -> forget about NLP and return */
      SCIPdebugMessage("No nonlinear continuous variables. NLP local search heuristic will not run.\n");
      return SCIP_OKAY;
   }

   /* setup NLP (if NLP solver available) */
#ifdef WITH_IPOPT
   SCIP_CALL( setupNLP(scip, heurdata) );
#else
   SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "No NLP solver available. NLP local search heuristic will not run.\n");
   havenlp = FALSE;
#endif
   
   return SCIP_OKAY;   
}

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

      /* integer variables have been fixed in heurExecNlp already
       * thus, we should not touch (and maybe unfix) these variables here
       * further, we believe that the given startpoint satisfies the varbound constraints, so that the current fixation is feasible 
       * @todo however, as Timo suggests, an assert would be nice here; unfortunately, this is not trivial to implement since there are no SCIPnlpiGetVarLb/Ub functions */
      if( SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
         continue;

      /* check if we passed already a varbound constraint on variable var */ 
      idx = SCIPhashmapGetImage(varmap, var);
      if( idx == NULL )
      { 
         /* variable appeared first time */
         assert(SCIPhashmapExists(heurdata->var_scip2nlp, var));
         varidx[varcnt] = (int)(size_t)SCIPhashmapGetImage(heurdata->var_scip2nlp, var);

         /* this constraint can only be handled explicitely, if the Vbdvar is discrete, i.e., fixed for an NLP solve */
         assert(SCIPvarGetType(SCIPgetVbdvarVarbound(scip, cons)) <= SCIP_VARTYPE_INTEGER);
         
         /* compute bounds on var determined by varbound constraint */
         varlb[varcnt] = SCIPgetLhsVarbound(scip, cons);
         varub[varcnt] = SCIPgetRhsVarbound(scip, cons);
       
         shift = SCIPgetVbdcoefVarbound(scip, cons) * SCIPgetSolVal(scip, heurdata->startcand, SCIPgetVbdvarVarbound(scip, cons));
         if( !SCIPisInfinity(scip, -varlb[varcnt]) )
            varlb[varcnt] -= shift;
         if( !SCIPisInfinity(scip,  varub[varcnt]) )
            varub[varcnt] -= shift;
         
         if (SCIPvarGetLbGlobal(var) > varlb[varcnt])
            varlb[varcnt] = SCIPvarGetLbGlobal(var);
         if (SCIPvarGetUbGlobal(var) < varub[varcnt])
            varub[varcnt] = SCIPvarGetUbGlobal(var);
         if( varlb[varcnt] > varub[varcnt] )
            varlb[varcnt] = varub[varcnt];
         
         /* remember that a bound change for variable var is now stored at position varbnd */
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

         idx_ = ((int)(size_t)idx) - 1;
         assert(idx_ < varcnt);

         /* compute bounds on var determined by this varbound constraint */
         lhs = SCIPgetLhsVarbound(scip, cons);
         rhs = SCIPgetRhsVarbound(scip, cons);

         shift = SCIPgetVbdcoefVarbound(scip, cons) * SCIPgetSolVal(scip, heurdata->startcand, SCIPgetVbdvarVarbound(scip, cons));
         if( !SCIPisInfinity(scip, -lhs) )
            lhs -= shift;
         if( !SCIPisInfinity(scip,  rhs) )
            rhs -= shift;

         /* possibly tighten previously stored bound change on variable var with newly computed bounds */ 
         varlb[idx_] = MAX(varlb[idx_],lhs);
         varub[idx_] = MIN(varub[idx_],rhs);
         if( varlb[idx_] > varub[idx_] )
            varlb[idx_] = varub[idx_];
   
         SCIPdebugMessage("%s: var %s at %d now bounded in [%g, %g] due to %s = %g  [updated]\n",
            SCIPconsGetName(cons), SCIPvarGetName(var), varidx[idx_],
            varlb[idx_], varub[idx_], SCIPvarGetName(SCIPgetVbdvarVarbound(scip, cons)),
            SCIPgetSolVal(scip, heurdata->startcand, SCIPgetVbdvarVarbound(scip, cons)) );
      }      
   }
   
   /* apply bound changes on variables in varbound constraint to NLP */
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
   {
      SCIP_CALL( SCIPdropVarEvent(scip, heurdata->var_nlp2scip[i], SCIP_EVENTTYPE_GBDCHANGED, heurdata->eventhdlr, (SCIP_EVENTDATA*)heurdata, -1) );
      SCIP_CALL( SCIPreleaseVar(scip, &heurdata->var_nlp2scip[i]) );
   }
   
   SCIPfreeMemoryArray(scip, &heurdata->var_nlp2scip);
   SCIPhashmapFree(&heurdata->var_scip2nlp);
   
   SCIPfreeMemoryArray(scip, &heurdata->discrvars);

   SCIP_CALL( SCIPnlpiFree(scip, &heurdata->nlpi) );
   assert(heurdata->nlpi == NULL);
   
   SCIPnlpStatisticsFree(scip, &heurdata->nlpstatistics);
   assert(heurdata->nlpstatistics == NULL);
   
   return SCIP_OKAY;
}


/** main procedure of the NLP heuristic */
SCIP_RETCODE SCIPapplyNlpHeur(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                       */
   SCIP_RESULT*          result,             /**< result data structure                                          */
   SCIP_SOL*             refpoint,           /**< point to take fixation of discrete variables from, and startpoint for NLP solver; if NULL, then LP solution is used */
   SCIP_Longint          itercontingent,     /**< iteration limit for NLP solver                                 */
   SCIP_Real             timelimit,          /**< time limit for NLP solver                                      */
   SCIP_Longint*         iterused            /**< buffer to store number of iterations used by NLP solver, or NULL if not of interest */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_Real*     startpoint;

   assert(scip != NULL);
   assert(heur != NULL);
  
   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* try to setup NLP if not tried before */
   if( heurdata->nlpi == NULL && !heurdata->triedsetupnlp )
   {
      SCIP_CALL( checkCIPandSetupNLP(scip, heur) );
   }

   /* not initialized */
   if( heurdata->nlpi == NULL )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   assert(heurdata->nvars > 0);
   assert(heurdata->var_nlp2scip != NULL);
   
   if( iterused != NULL )
      *iterused = 0;

   /* get starting values (=refpoint, if not NULL; otherwise LP solution (or pseudo solution)) */
   SCIP_CALL( SCIPallocBufferArray(scip, &startpoint, heurdata->nvars) );
   SCIP_CALL( SCIPgetSolVals(scip, refpoint, heurdata->nvars, heurdata->var_nlp2scip, startpoint) );

   /* fix discrete variables to values in startpoint */
   if( heurdata->ndiscrvars > 0 )
   {
      SCIP_Real* discrfix;
      int i;

      SCIP_CALL( SCIPallocBufferArray(scip, &discrfix, heurdata->ndiscrvars) );

      for( i = 0; i < heurdata->ndiscrvars; ++i )
      {
         discrfix[i] = startpoint[heurdata->discrvars[i]];

         /* only apply to integer feasible points */
         if( !SCIPisFeasIntegral(scip, discrfix[i]) )
         {
            if( refpoint || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
            {
               SCIPdebugMessage("skip NLP heuristic because start candidate not integer feasible: var %d is %g\n", i, discrfix[i]);
               SCIPfreeBufferArray(scip, &startpoint);
               *result = SCIP_DIDNOTRUN;
               return SCIP_OKAY;
            }
            /* otherwise we desperately wanna run the NLP heur, so we continue and round what we have */
         }
         /* if we do not really have a startpoint, then we should take care that we do not fix variables to very large values
          * set to 0 here and project on bounds below */ 
         if( ABS(discrfix[i]) > 1E+10 && !refpoint && SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
            discrfix[i] = 0;

         /* round fractional variables to the nearest integer */
         discrfix[i] = SCIPfrac(scip, discrfix[i]) > 0.5 ? SCIPceil(scip, discrfix[i]) : SCIPfloor(scip, discrfix[i]);

         /* adjust value to the global bounds of the corresponding SCIP variable */
         discrfix[i] = MAX(discrfix[i], SCIPvarGetLbGlobal(heurdata->var_nlp2scip[heurdata->discrvars[i]]));
         discrfix[i] = MIN(discrfix[i], SCIPvarGetUbGlobal(heurdata->var_nlp2scip[heurdata->discrvars[i]]));

         assert(discrfix[i] == SCIPfloor(scip, discrfix[i]));
         assert(discrfix[i] == SCIPceil(scip, discrfix[i]));

         startpoint[heurdata->discrvars[i]] = discrfix[i];
      }
      /* apply fixings */
      SCIP_CALL( SCIPnlpiChgVarBounds(scip, heurdata->nlpi, heurdata->ndiscrvars, heurdata->discrvars, discrfix, discrfix) );

      SCIPfreeBufferArray(scip, &discrfix);
   }
   /* apply those variable bound constraints that we can apply explicitely */
   SCIP_CALL( applyVarBoundConstraints(scip, heurdata) );
   
   /* set time and iteration limit for NLP solver */
   SCIP_CALL( SCIPnlpiSetIntPar(scip, heurdata->nlpi, SCIP_NLPPAR_ITLIM, (int)itercontingent) );
   SCIP_CALL( SCIPnlpiSetRealPar(scip, heurdata->nlpi, SCIP_NLPPAR_TILIM, timelimit) );

   /* pass initial guess to NLP solver, if we have one; otherwise clear previous guess and let NLP solver choose */
   if( refpoint || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIP_CALL( SCIPnlpiSetInitialGuess(scip, heurdata->nlpi, startpoint) );
   }
   else
   {
      SCIP_CALL( SCIPnlpiSetInitialGuess(scip, heurdata->nlpi, NULL) );
   }
   
   SCIPfreeBufferArray(scip, &startpoint);

   /* let the NLP solver do its magic */
   SCIPdebugMessage("start NLP solve with iteration limit %"SCIP_LONGINT_FORMAT" and timelimit %g\n", itercontingent, timelimit);
   SCIP_CALL( SCIPnlpiSolve(scip, heurdata->nlpi) );

   SCIPdebugMessage("NLP solver returned with termination status %d and solution status %d\n", 
      SCIPnlpiGetTermstat(scip, heurdata->nlpi), SCIPnlpiGetSolstat(scip, heurdata->nlpi));
   
   if( SCIPnlpiGetTermstat(scip, heurdata->nlpi) >= SCIP_NLPTERMSTAT_MEMERR )
   {  /* oops, something did not go well at all */
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, 
         "NLP solver returned with bad termination status %d. Will not run NLP heuristic again for this run.\n",  
         SCIPnlpiGetTermstat(scip, heurdata->nlpi));
      SCIP_CALL( destroyNLP(scip, heurdata) );
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPnlpiGetStatistics(scip, heurdata->nlpi, heurdata->nlpstatistics) );

   if( iterused != NULL )
      *iterused += SCIPnlpStatisticsGetNIterations(heurdata->nlpstatistics);
   SCIPdebugMessage("NLP solver used %d iterations and %g seconds\n", 
      SCIPnlpStatisticsGetNIterations(heurdata->nlpstatistics), SCIPnlpStatisticsGetTotalTime(heurdata->nlpstatistics));

   if( SCIPnlpiGetSolstat(scip, heurdata->nlpi) <= SCIP_NLPSOLSTAT_FEASIBLE )
   {  /* NLP solver claims, it found a feasible (maybe even optimal) solution */
      SCIP_Real* primals;
      SCIP_SOL*  sol;
      SCIP_Bool  stored;
      SCIP_Bool  feasible;

      /* get solution from NLP solver and create a SCIP_SOL out of it */
      SCIP_CALL( SCIPnlpiGetSolution(scip, heurdata->nlpi, &primals) );
      assert(primals != NULL);

      SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
      SCIP_CALL( SCIPsetSolVals(scip, sol, heurdata->nvars, heurdata->var_nlp2scip, primals) );

      SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, TRUE, TRUE, &feasible) );

      if( feasible )
      {  /* SCIP agrees that solution is feasible (yippi!), so we add it */ 
         SCIP_CALL( SCIPaddSol(scip, sol, &stored) );
         if( stored )
         {
            SCIPdebugMessage("SCIP stored solution with value %g\n", SCIPgetSolOrigObj(scip, sol));
            *result = SCIP_FOUNDSOL;
         }
      }
      else if( heurdata->resolvetolfactor < 1.0 )
      {
         /* resolve with tighter tolerances */
         SCIPdebugMessage("solution reported by NLP solver not feasible for SCIP, resolve with feasibility tolerance %g\n", heurdata->resolvetolfactor*SCIPfeastol(scip));
         SCIP_CALL( SCIPnlpiSetRealPar(scip, heurdata->nlpi, SCIP_NLPPAR_FEASTOL, heurdata->resolvetolfactor*SCIPfeastol(scip)) );

         /* set almost-feasible solution as starting point for NLP solver, if user allows us */
         if( !heurdata->resolvefromscratch )
         {
            SCIP_CALL( SCIPnlpiSetInitialGuess(scip, heurdata->nlpi, primals) );
         }

         /* solve again */
         SCIP_CALL( SCIPnlpiSolve(scip, heurdata->nlpi) );
         SCIP_CALL( SCIPnlpiGetStatistics(scip, heurdata->nlpi, heurdata->nlpstatistics) );
         if( iterused != NULL )
            *iterused += SCIPnlpStatisticsGetNIterations(heurdata->nlpstatistics);

         if( SCIPnlpiGetSolstat(scip, heurdata->nlpi) <= SCIP_NLPSOLSTAT_FEASIBLE )
         { 
            /* NLP solver claims feasible again, hope that SCIP accepts it now */
            SCIP_CALL( SCIPnlpiGetSolution(scip, heurdata->nlpi, &primals) );
            assert(primals != NULL);

            SCIP_CALL( SCIPsetSolVals(scip, sol, heurdata->nvars, heurdata->var_nlp2scip, primals) );

#ifdef SCIP_DEBUG
            /* print the infeasibilities to stdout */
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
   
   /* TODO: reset time and iterlimit in nlp solver? */

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
   SCIP_HEURDATA* heurdata;
   
   assert(scip != NULL);
   assert(heur != NULL);
   
   /* do not setup NLP if heuristic is never called by SCIP */
   if( SCIPheurGetFreq(heur) < 0 )
      return SCIP_OKAY;

   /* try to setup NLP; fails if there are no nonlinear continuous variables or there is no NLP solver */
   SCIP_CALL( checkCIPandSetupNLP(scip, heur) );

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   
   if( heurdata->nlpi == NULL )
      return SCIP_OKAY;
   
   /* if the heuristic is called at the root node, we want to be called directly after the initial root LP solve */
   if( SCIPheurGetFreqofs(heur) == 0 )
      SCIPheurSetTimingmask(heur, SCIP_HEURTIMING_DURINGLPLOOP | HEUR_TIMING);

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

   SCIPheurSetTimingmask(heur, HEUR_TIMING);
   
   heurdata->triedsetupnlp = FALSE;

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecNlp)
{  /*lint --e{666,715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_Longint   itercontingent;
   SCIP_Real      timelimit;
   SCIP_Longint   iterused;

   assert(scip != NULL);
   assert(heur != NULL);
  
   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   
   /* obviously, we did not do anything yet */
   *result = SCIP_DIDNOTRUN;
   
   /* InitsolNlp decided that we do not need an NLP solver
    * probably because we do not have nonlinear continuous or implicit integer variables */
   if( heurdata->nlpi == NULL )
      return SCIP_OKAY;
   
   if( heurdata->startcand == NULL )
   {  /* if no start candidate is given, we like to consider the LP solution of the current node */
      /* at least if we are not called the first time, we call the heuristic only if an optimal LP solution is available 
       * if we are called the first time and the LP is unbounded, then we are quite desperate and still give the NLP a try */
      if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      {
         if( SCIPgetNNodes(scip) > 1 || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_UNBOUNDEDRAY )
         {
            *result = SCIP_DELAYED;
            SCIPdebugMessage("NLP heuristic delayed because no start candidate given and no LP solution available\n");
            return SCIP_OKAY;
         }
         else
         {
            SCIPdebugMessage("LP is unbounded in root node, so we are quite desperate; run NLP heuristic and pray\n");
         }
      }
      else if( SCIPgetNLPBranchCands(scip) > 0 )
      { /* only call heuristic, if there are no fractional variables */
         *result = SCIP_DELAYED;
         SCIPdebugMessage("NLP heuristic delayed because no start candidate given and current LP solution is fractional\n");
         return SCIP_OKAY;
      }
   }
   else
   {
      SCIPdebugMessage("have startcand from heur %s\n", SCIPsolGetHeur(heurdata->startcand) ? SCIPheurGetName(SCIPsolGetHeur(heurdata->startcand)) : "NULL");
   }

   /* compute the contingent on number of iterations that the NLP solver is allowed to use
      we make it depending on the current number of processed nodes */  
   itercontingent = (SCIP_Longint)(heurdata->iterquot * SCIPgetNNodes(scip));

   /* weight by previous success of heuristic */
   itercontingent = (SCIP_Longint)(itercontingent * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   /* add the fixed offset */
   itercontingent += heurdata->iteroffset;
   /* substract the number of iterations used so far */
   itercontingent -= heurdata->iterused;

   if( itercontingent < heurdata->itermin )
   {  /* not enough iterations left to start NLP solver */
      SCIPdebugMessage("skip NLP heuristic; contingent=%"SCIP_LONGINT_FORMAT"; minimal number of iterations=%d; success ratio=%g\n", 
         itercontingent, heurdata->itermin, (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
      return SCIP_OKAY;
   }

   /* enforce user given iteration limit, if given */
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
      /* enforce user given time limit, if given */
      if( heurdata->nlptimelimit > 0 )
         timelimit = MIN(heurdata->nlptimelimit, timelimit);

      SCIP_CALL( SCIPnlpiSetRealPar(scip, heurdata->nlpi, SCIP_NLPPAR_TILIM, timelimit) );
   }
   else if( heurdata->nlptimelimit > 0 )
   {  /* enforce user given time limit */
      SCIP_CALL( SCIPnlpiSetRealPar(scip, heurdata->nlpi, SCIP_NLPPAR_TILIM, heurdata->nlptimelimit) );
   }

   /* so far we have not found any solution, but now we are willing to search for one */
   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPapplyNlpHeur(scip, heur, result, heurdata->startcand, itercontingent, timelimit, &iterused) );
   heurdata->iterused += iterused;
   
   /* forget startcand */
   if( heurdata->startcand != NULL )
      SCIP_CALL( SCIPfreeSol(scip, &heurdata->startcand) );
   
   /* reset timing, if it was changed temporary (at the root node) */
   if( heurtiming != HEUR_TIMING )
      SCIPheurSetTimingmask(heur, HEUR_TIMING);

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

   /* include variable event handler */
   SCIP_CALL( SCIPincludeEventhdlr(scip, HEUR_NAME, "propagates a global bound change to the NLP",
      NULL, NULL, NULL, NULL, NULL, NULL, processVarEvent, NULL) );
   heurdata->eventhdlr = SCIPfindEventhdlr(scip, HEUR_NAME);

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, heurFreeNlp, heurInitNlp, heurExitNlp, heurInitsolNlp, heurExitsolNlp, heurExecNlp,
         heurdata) );

   /* add Nlp primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/"HEUR_NAME"/nlpverblevel",
         "verbosity level of NLP solver",
         &heurdata->nlpverblevel, FALSE, 0, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/"HEUR_NAME"/nlpiterlimit",
         "iteration limit of NLP solver; 0 to use solver default",
         &heurdata->nlpiterlimit, FALSE, 0, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/nlptimelimit",
         "time limit of NLP solver; 0 to use solver default",
         &heurdata->nlptimelimit, FALSE, 0.0, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/resolvetolfactor",
         "if SCIP does not accept a NLP feasible solution, resolve NLP with feas. tolerance reduced by this factor (set to 1.0 to turn off resolve)",
         &heurdata->resolvetolfactor, TRUE, 0.01, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/resolvefromscratch",
         "should the NLP resolve be started from the original starting point or the infeasible solution?",
         &heurdata->resolvefromscratch, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/"HEUR_NAME"/iteroffset",
         "number of iterations added to the contingent of the total number of iterations",
         &heurdata->iteroffset, FALSE, 500, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/iterquotient",
         "contingent of NLP iterations in relation to the number of nodes in SCIP",
         &heurdata->iterquot,  FALSE, 0.1, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/"HEUR_NAME"/itermin",
         "contingent of NLP iterations in relation to the number of nodes in SCIP",
         &heurdata->itermin,   FALSE, 300, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/varboundexplicit",
         "should variable bound constraints be handled explicitly before solving the NLP instead of adding them to the NLP?",
         &heurdata->varboundexplicit, TRUE, TRUE, NULL, NULL) );

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
   
   assert(scip != NULL);
   assert(heur != NULL);
   assert(solcand != NULL);
   assert(SCIPisPositive(scip, violation));
   
   /* game is over already; no more interest in starting points */
   if( SCIPgetStage(scip) >= SCIP_STAGE_SOLVED )
      return SCIP_OKAY;
   
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   
   /* we do not have a NLP, so we also do not need a starting point */
   if( heurdata->nlpi == NULL )
      return SCIP_OKAY;
   
   /* if we have no point yet, or the new point has a lower constraint violation, or it has a better objective function value,
    * then take the new point */
   if( heurdata->startcand == NULL || SCIPisGT(scip, heurdata->startcandviol, violation) ||
      SCIPisRelGT(scip, SCIPgetSolTransObj(scip, heurdata->startcand)*(int)SCIPgetObjsense(scip), SCIPgetSolTransObj(scip, solcand)*(int)SCIPgetObjsense(scip)) )
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

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
#pragma ident "@(#) $Id: heur_nlp.c,v 1.76 2010/08/30 16:50:07 bzfwinkm Exp $"

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
#include "nlpi/nlpi.h"

#include "scip/cons_linear.h"
#include "scip/cons_varbound.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_indicator.h"

#define HEUR_NAME             "nlp"
#define HEUR_DESC             "primal heuristic that performs a local search in an NLP after fixing integer variables"
#define HEUR_DISPCHAR         'Q'
#define HEUR_PRIORITY         -2000000
#define HEUR_FREQ             1
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
   SCIP_NLPIPROBLEM*     nlpiprob;           /**< problem in NLP solver interface */
   SCIP_NLPSTATISTICS*   nlpstatistics;      /**< statistics from NLP solver */
   SCIP_Bool             triedsetupnlp;      /**< whether we have tried to setup an NLP */ 
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for global bound change events */
   SCIP_DECL_HEURNLPNLPIINIT((**nlpiinits)); /**< user-given NLPI initialization methods */
   SCIP_DECL_HEURNLPHAVECONS((**haveconss)); /**< user-given NLPI information methods */
   int                   nlpiinitssize;      /**< size of nlpiinits array */
   int                   nnlpiinits;         /**< number of NLPI initialization methods */
   
   int                   nvars;              /**< number of variables in NLP */
   SCIP_VAR**            var_nlp2scip;       /**< mapping variables in NLP to SCIP variables */
   SCIP_HASHMAP*         var_scip2nlp;       /**< mapping variables in SCIP to NLP variables */
                         
   int                   ndiscrvars;         /**< number of discrete variables */
   int*                  discrvars;          /**< indices of discrete variables */
                         
   int                   nvarbndconss;       /**< number of variable bound constraints */
   SCIP_CONS**           varbndconss;        /**< variable bound constraints */
   
   SCIP_CONSHDLR*        conshdlrindicator;  /**< constraint handler for indicator constraints */
                         
   SCIP_SOL*             startcand;          /**< candidate for start point for heuristic */
   SCIP_Real             startcandviol;      /**< violation of start point candidate w.r.t. constraint that reported this candidate */
                         
   int                   nlpverblevel;       /**< verbosity level of NLP solver */
   int                   nlpiterlimit;       /**< iteration limit of NLP solver; 0 for off */
   SCIP_Real             nlptimelimit;       /**< time limit of NLP solver; 0 for off */
   SCIP_Real             resolvetolfactor;   /**< factor for feasiblity tolerance when resolving NLP due to disagreement of feasibility */
   SCIP_Bool             resolvefromscratch; /**< whether a resolve of an NLP due to disagreement of feasibility should be from the original starting point or the infeasible solution */
   char*                 nlpsolver;          /**< name of NLP solver to use */
   char*                 nlpoptfile;         /**< name of NLP solver specific option file */
   SCIP_Bool             names;              /**< whether to pass variable and constraint names to the NLP */
                         
   SCIP_Longint          iterused;           /**< number of iterations used so far */
   int                   iteroffset;         /**< number of iterations added to the contingent of the total number of iterations */
   SCIP_Real             iterquot;           /**< contingent of NLP iterations in relation to the number of nodes in SCIP */
   int                   itermin;            /**< minimal number of iterations required to start local search */
   SCIP_Bool             runalways;          /**< whether to run NLP heuristic always (independent of iteroffset,iterquot,itermin) */
                         
   SCIP_Bool             varboundexplicit;   /**< whether variable bound constraints should be handled explicitly before solving an NLP instead of adding them as linear constraints to to the NLP */
};


/*
 * Local methods
 */

/** adds linear constraints to NLP */
static
SCIP_RETCODE addLinearConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< problem in NLP solver */
   SCIP_HASHMAP*         varsmap,            /**< mapping SCIP variables to variable indices in NLP */
   SCIP_HASHMAP*         conssmap,           /**< hashmap where to add mapping from SCIP constraints to row indices in NLP, or NULL */
   int*                  nlpconsscounter,    /**< counter where to add number of constraints added to NLP, or NULL; need to be not-NULL if conssmap != NULL */
   SCIP_CONSHDLR*        linconshdlr,        /**< constraint handler for linear constraints */
   SCIP_Bool             onlysubnlp,         /**< whether to add only constraints that are relevant for the NLP obtained by fixing all discrete variables in the CIP */
   SCIP_Bool             names               /**< whether to gives variable and constraint names to NLPI */
   )
{
   SCIP_CONS**      conss;
   SCIP_CONS**      useconss;
   SCIP_Real*       lhs;
   SCIP_Real*       rhs;
   const char**     consnames;
   int*             nlininds;
   int**            lininds;
   SCIP_Real**      linvals;

   int              nconss;
   int              nuseconss;
   int              i;
   int              k;

   assert(scip != NULL);
   assert(nlpi != NULL);
   assert(nlpiprob != NULL);
   assert(varsmap != NULL);
   assert(linconshdlr != NULL);
   assert(conssmap == NULL || nlpconsscounter != NULL);
   
   nconss = SCIPconshdlrGetNConss(linconshdlr);
   conss  = SCIPconshdlrGetConss(linconshdlr);
   
   if( !nconss )
      return SCIP_OKAY;
   
   SCIP_CALL( SCIPallocBufferArray(scip, &useconss, nconss) );
   
   /* count the number of constraints to add to NLP
    * if onlysubnlp, then constraints with only integer variables are not added to the NLP, since all its variables will be fixed later */
   nuseconss = 0;
   for( i = 0; i < nconss; ++i )
   {
      if( SCIPconsIsLocal(conss[i]) )
         continue;
      if( onlysubnlp )
      {
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
      else
      {
         useconss[nuseconss] = conss[i];
         ++nuseconss;
      }
   }
   
   if( !nuseconss )
   {  /* no linear constraints to add to NLP, so we are done */
      SCIPfreeBufferArray(scip, &useconss);
      return SCIP_OKAY;
   }
   
   SCIP_CALL( SCIPallocBufferArray(scip, &lhs, nuseconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhs, nuseconss) );
   if( names )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &consnames, nuseconss) );
   }
   else
      consnames = NULL;
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
      if( names )
         consnames[i] = SCIPconsGetName(useconss[i]);
      
      SCIP_CALL( SCIPallocBufferArray(scip, &lininds[i], nvars) );
      nlininds[i] = nvars;
      linvals[i]  = vals;
     
      /* fill column indices of variables (vars) in conss[i] with variable indices in NLP */
      for( k = 0; k < nvars; ++k )
      {
         assert(SCIPhashmapExists(varsmap, vars[k]));
         lininds[i][k] = (int) (size_t) SCIPhashmapGetImage(varsmap, vars[k]);
      }
      
      if( conssmap != NULL )
      {
         SCIP_CALL( SCIPhashmapInsert(conssmap, useconss[i], (void*)(size_t)*nlpconsscounter) );
      }
      if( nlpconsscounter != NULL )
         ++*nlpconsscounter;
   }
   
   SCIP_CALL( SCIPnlpiAddConstraints(nlpi, nlpiprob, nuseconss,
         lhs, rhs, nlininds, lininds, linvals,
         NULL, NULL, NULL, NULL, consnames) );
   
   for( i = 0; i < nuseconss; ++i )
   {
      SCIPfreeBufferArray(scip, &lininds[i]);
   }
   
   SCIPfreeBufferArray(scip, &lininds);
   SCIPfreeBufferArray(scip, &linvals);
   SCIPfreeBufferArray(scip, &nlininds);
   SCIPfreeBufferArray(scip, &rhs);
   SCIPfreeBufferArray(scip, &lhs);
   SCIPfreeBufferArrayNull(scip, &consnames);
   SCIPfreeBufferArray(scip, &useconss);
   
   return SCIP_OKAY;
}

/** adds linear constraints to NLP */
static
SCIP_RETCODE addLogicOrConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< problem in NLP solver */
   SCIP_HASHMAP*         varsmap,            /**< mapping SCIP variables to variable indices in NLP */
   SCIP_HASHMAP*         conssmap,           /**< hashmap where to add mapping from SCIP constraints to row indices in NLP, or NULL */
   int*                  nlpconsscounter,    /**< counter where to add number of constraints added to NLP, or NULL; need to be not-NULL if conssmap != NULL */
   SCIP_CONSHDLR*        logicorconshdlr,    /**< constraint handler for logic or constraints */
   SCIP_Bool             names               /**< whether to gives variable and constraint names to NLPI */
   )
{
   SCIP_CONS**           conss;              /* all logicor constraints */
   int                   nconss;             /* total number of logicor constraints */
   int                   nconss4nlp;
   int                   i, j, k;
   SCIP_Real*            lhs;
   const char**          consnames;
   int*                  nlininds;
   int**                 lininds;
   SCIP_Real**           linvals;

   assert(scip != NULL);
   assert(nlpi != NULL);
   assert(nlpiprob != NULL);
   assert(varsmap != NULL);
   assert(logicorconshdlr != NULL);
   assert(conssmap == NULL || nlpconsscounter != NULL);

   nconss = SCIPconshdlrGetNConss(logicorconshdlr);
   if( !nconss )
      return SCIP_OKAY;

   conss = SCIPconshdlrGetConss(logicorconshdlr);

   nconss4nlp = 0;
   for( i = 0; i < nconss; ++i )
   {
      if( SCIPconsIsLocal(conss[i]) )
         continue;
      nconss4nlp++;
   }

   if( nconss4nlp == 0 )
   { /* all logic or constraints are local, nothing to add to NLP */
      return SCIP_OKAY;
   }

   /* memory to store those logic or constraints that are added to the NLP as linear constraint */
   SCIP_CALL( SCIPallocBufferArray(scip, &lhs,      nconss4nlp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlininds, nconss4nlp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lininds,  nconss4nlp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals,  nconss4nlp) );
   if( names )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &consnames, nconss4nlp) );
   }
   else
      consnames = NULL;

   k = 0; /* the number of logic or constraints for the NLP passed so far */
   for( i = 0; i < nconss; ++i )
   {
      /* skip local constraints */
      if( SCIPconsIsLocal(conss[i]) )
         continue;

      /* logic or constraints: 1 == sum_j x_j */
      lhs[k] = 1.0;

      nlininds[k] = SCIPgetNVarsLogicor(scip, conss[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &lininds[k], nlininds[k]) );
      SCIP_CALL( SCIPallocBufferArray(scip, &linvals[k], nlininds[k]) );
      if( names )
         consnames[k] = SCIPconsGetName(conss[i]);

      for( j = 0; j < nlininds[k]; ++j )
      {
         linvals[k][j] = 1.0;

         if( !SCIPhashmapExists(varsmap, SCIPgetVarsLogicor(scip, conss[i])[j]) )
            break;
         lininds[k][j] = (int) (size_t) SCIPhashmapGetImage(varsmap, SCIPgetVarsLogicor(scip, conss[i])[j]);
      }

      if( j < nlininds[k] )
      {
         SCIPfreeBufferArray(scip, &lininds[k]);
         SCIPfreeBufferArray(scip, &linvals[k]);
         --nconss4nlp;
         continue;
      }
      
      if( conssmap != NULL )
      {
         SCIP_CALL( SCIPhashmapInsert(conssmap, conss[i], (void*)(size_t)*nlpconsscounter) );
      }
      if( nlpconsscounter != NULL )
         ++*nlpconsscounter;

      ++k;
   }
   assert(k == nconss4nlp);

   SCIP_CALL( SCIPnlpiAddConstraints(nlpi, nlpiprob, nconss4nlp,
      lhs, lhs, nlininds, lininds, linvals,
      NULL, NULL, NULL, NULL, consnames) );

   for( k = 0; k < nconss4nlp; ++k )
   {
      SCIPfreeBufferArray(scip, &lininds[k]);
      SCIPfreeBufferArray(scip, &linvals[k]);
   }

   SCIPfreeBufferArray(scip, &nlininds);
   SCIPfreeBufferArray(scip, &lininds);
   SCIPfreeBufferArray(scip, &linvals);
   SCIPfreeBufferArray(scip, &lhs);
   SCIPfreeBufferArrayNull(scip, &consnames);

   return SCIP_OKAY;
}

/** add setppc constraints */
static
SCIP_RETCODE addSetppcConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< problem in NLP solver */
   SCIP_HASHMAP*         varsmap,            /**< mapping SCIP variables to variable indices in NLP */
   SCIP_HASHMAP*         conssmap,           /**< hashmap where to add mapping from SCIP constraints to row indices in NLP, or NULL */
   int*                  nlpconsscounter,    /**< counter where to add number of constraints added to NLP, or NULL; need to be not-NULL if conssmap != NULL */
   SCIP_CONSHDLR*        setppcconshdlr,     /**< constraint handler for setppc constraints */
   SCIP_Bool             names               /**< whether to gives variable and constraint names to NLPI */
   )
{
   SCIP_CONS**           conss;              /* all setppc constraints */
   int                   nconss;             /* total number of setpcc constraints */
   int                   nconss4nlp;
   int                   i, j, k;
   SCIP_Real*            lhs;
   SCIP_Real*            rhs;
   const char**          consnames;
   int*                  nlininds;
   int**                 lininds;
   SCIP_Real**           linvals;

   assert(scip != NULL);
   assert(nlpi != NULL);
   assert(nlpiprob != NULL);
   assert(varsmap != NULL);
   assert(setppcconshdlr != NULL);
   assert(conssmap == NULL || nlpconsscounter != NULL);

   nconss = SCIPconshdlrGetNConss(setppcconshdlr);
   if( !nconss )
      return SCIP_OKAY;

   conss = SCIPconshdlrGetConss(setppcconshdlr);

   nconss4nlp = 0;
   for( i = 0; i < nconss; ++i )
   {
      if( SCIPconsIsLocal(conss[i]) )
         continue;
      nconss4nlp++;
   }

   if( nconss4nlp == 0 )
   { /* all setppc constraints are local, nothing to add to NLP */
      return SCIP_OKAY;
   }

   /* memory to store those setppc constraints that are added to the NLP as linear constraint */
   SCIP_CALL( SCIPallocBufferArray(scip, &lhs,      nconss4nlp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhs,      nconss4nlp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlininds, nconss4nlp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lininds,  nconss4nlp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals,  nconss4nlp) );
   if( names )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &consnames, nconss4nlp) );
   }
   else
      consnames = NULL;

   k = 0; /* the number of setppc constraints for the NLP passed so far */
   for( i = 0; i < nconss; ++i )
   {
      /* skip local constraints */
      if( SCIPconsIsLocal(conss[i]) )
         continue;

      /* setppc constraints: 1 ~ sum_j x_j */
      switch( SCIPgetTypeSetppc(scip, conss[i]) )
      {
         case SCIP_SETPPCTYPE_PARTITIONING:
            lhs[k] = 1.0;
            rhs[k] = 1.0;
            break;

         case SCIP_SETPPCTYPE_PACKING:
            lhs[k] = -SCIPinfinity(scip);
            rhs[k] = 1.0;
            break;

         case SCIP_SETPPCTYPE_COVERING:
            lhs[k] = 1.0;
            rhs[k] = SCIPinfinity(scip);
            break;
      }

      nlininds[k] = SCIPgetNVarsSetppc(scip, conss[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &lininds[k], nlininds[k]) );
      SCIP_CALL( SCIPallocBufferArray(scip, &linvals[k], nlininds[k]) );

      for( j = 0; j < nlininds[k]; ++j )
      {
         linvals[k][j] = 1.0;

         if( !SCIPhashmapExists(varsmap, SCIPgetVarsSetppc(scip, conss[i])[j]) )
            break;
         lininds[k][j] = (int) (size_t) SCIPhashmapGetImage(varsmap, SCIPgetVarsSetppc(scip, conss[i])[j]);
      }

      if( j < nlininds[k] )
      {
         SCIPfreeBufferArray(scip, &lininds[k]);
         SCIPfreeBufferArray(scip, &linvals[k]);
         --nconss4nlp;
         continue;
      }
      
      if( names )
         consnames[k] = SCIPconsGetName(conss[i]);

      if( conssmap != NULL )
      {
         SCIP_CALL( SCIPhashmapInsert(conssmap, conss[i], (void*)(size_t)*nlpconsscounter) );
      }
      if( nlpconsscounter != NULL )
         ++*nlpconsscounter;

      ++k;
   }
   assert(k == nconss4nlp);

   SCIP_CALL( SCIPnlpiAddConstraints(nlpi, nlpiprob, nconss4nlp,
      lhs, rhs, nlininds, lininds, linvals,
      NULL, NULL, NULL, NULL, consnames) );

   for( k = 0; k < nconss4nlp; ++k )
   {
      SCIPfreeBufferArray(scip, &lininds[k]);
      SCIPfreeBufferArray(scip, &linvals[k]);
   }

   SCIPfreeBufferArray(scip, &nlininds);
   SCIPfreeBufferArray(scip, &lininds);
   SCIPfreeBufferArray(scip, &linvals);
   SCIPfreeBufferArray(scip, &lhs);
   SCIPfreeBufferArray(scip, &rhs);
   SCIPfreeBufferArrayNull(scip, &consnames);

   return SCIP_OKAY;
}

/** add knapsack constraints */
static
SCIP_RETCODE addKnapsackConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< problem in NLP solver */
   SCIP_HASHMAP*         varsmap,            /**< mapping SCIP variables to variable indices in NLP */
   SCIP_HASHMAP*         conssmap,           /**< hashmap where to add mapping from SCIP constraints to row indices in NLP, or NULL */
   int*                  nlpconsscounter,    /**< counter where to add number of constraints added to NLP, or NULL; need to be not-NULL if conssmap != NULL */
   SCIP_CONSHDLR*        knapsackconshdlr,   /**< constraint handler for knapsack constraints */
   SCIP_Bool             names               /**< whether to gives variable and constraint names to NLPI */
   )
{
   SCIP_CONS**           conss;              /* all knapsack constraints */
   int                   nconss;             /* total number of knapsack constraints */
   int                   nconss4nlp;
   int                   i, j, k;
   SCIP_Real*            lhs;
   SCIP_Real*            rhs;
   const char**          consnames;
   int*                  nlininds;
   int**                 lininds;
   SCIP_Real**           linvals;

   assert(scip != NULL);
   assert(nlpi != NULL);
   assert(nlpiprob != NULL);
   assert(varsmap != NULL);
   assert(knapsackconshdlr != NULL);
   assert(conssmap == NULL || nlpconsscounter != NULL);

   nconss = SCIPconshdlrGetNConss(knapsackconshdlr);
   if( !nconss )
      return SCIP_OKAY;

   conss = SCIPconshdlrGetConss(knapsackconshdlr);

   nconss4nlp = 0;
   for( i = 0; i < nconss; ++i )
   {
      if( SCIPconsIsLocal(conss[i]) )
         continue;
      nconss4nlp++;
   }

   if( nconss4nlp == 0 )
   { /* all knapsack constraints are local, nothing to add to NLP */
      return SCIP_OKAY;
   }

   /* memory to store those knapsack constraints that are added to the NLP as linear constraint */
   SCIP_CALL( SCIPallocBufferArray(scip, &lhs,      nconss4nlp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhs,      nconss4nlp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlininds, nconss4nlp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lininds,  nconss4nlp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals,  nconss4nlp) );
   if( names )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &consnames, nconss4nlp) );
   }
   else
      consnames = NULL;

   k = 0; /* the number of knapsack constraints for the NLP passed so far */
   for( i = 0; i < nconss; ++i )
   {
      /* skip local constraints */
      if( SCIPconsIsLocal(conss[i]) )
         continue;

      lhs[k] = -SCIPinfinity(scip);
      rhs[k] = (SCIP_Real)SCIPgetCapacityKnapsack(scip, conss[i]);

      nlininds[k] = SCIPgetNVarsKnapsack(scip, conss[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &lininds[k], nlininds[k]) );
      SCIP_CALL( SCIPallocBufferArray(scip, &linvals[k], nlininds[k]) );
      
      if( names )
         consnames[k] = SCIPconsGetName(conss[i]);

      for( j = 0; j < nlininds[k]; ++j )
      {
         linvals[k][j] = (SCIP_Real)SCIPgetWeightsKnapsack(scip, conss[i])[j];

         if( !SCIPhashmapExists(varsmap, SCIPgetVarsKnapsack(scip, conss[i])[j]) )
            break;
         lininds[k][j] = (int) (size_t) SCIPhashmapGetImage(varsmap, SCIPgetVarsKnapsack(scip, conss[i])[j]);
      }

      if( j < nlininds[k] )
      {
         SCIPfreeBufferArray(scip, &lininds[k]);
         SCIPfreeBufferArray(scip, &linvals[k]);
         --nconss4nlp;
         continue;
      }

      if( conssmap != NULL )
      {
         SCIP_CALL( SCIPhashmapInsert(conssmap, conss[i], (void*)(size_t)*nlpconsscounter) );
      }
      if( nlpconsscounter != NULL )
         ++*nlpconsscounter;

      ++k;
   }
   assert(k == nconss4nlp);

   SCIP_CALL( SCIPnlpiAddConstraints(nlpi, nlpiprob, nconss4nlp,
      lhs, rhs, nlininds, lininds, linvals,
      NULL, NULL, NULL, NULL, consnames) );

   for( k = 0; k < nconss4nlp; ++k )
   {
      SCIPfreeBufferArray(scip, &lininds[k]);
      SCIPfreeBufferArray(scip, &linvals[k]);
   }

   SCIPfreeBufferArray(scip, &nlininds);
   SCIPfreeBufferArray(scip, &lininds);
   SCIPfreeBufferArray(scip, &linvals);
   SCIPfreeBufferArray(scip, &lhs);
   SCIPfreeBufferArray(scip, &rhs);
   SCIPfreeBufferArrayNull(scip, &consnames);

   return SCIP_OKAY;
}

/** collect variable bound constraints */
static
SCIP_RETCODE collectVarBoundConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< problem in NLP solver */
   SCIP_HASHMAP*         varsmap,            /**< mapping SCIP variables to variable indices in NLP */
   SCIP_HASHMAP*         conssmap,           /**< hashmap where to add mapping from SCIP constraints to row indices in NLP, or NULL */
   int*                  nlpconsscounter,    /**< counter where to add number of constraints added to NLP, or NULL; need to be not-NULL if conssmap != NULL */
   int*                  nexplvbndconss,     /**< buffer to store number of variable bound constraints that need to be handled explicitely by the user, or NULL */
   SCIP_CONS***          explvbndconss,      /**< buffer to store array of variable bound constraints that need to be handled explicitely by the user, or NULL */
   SCIP_CONSHDLR*        varbndconshdlr,     /**< constraint handler for variable bound constraints */
   SCIP_Bool             names               /**< whether to gives variable and constraint names to NLPI */
   )
{
   SCIP_CONS**           conss;              /* all varbound constraints */

   int                   nconss;             /* total number of varbound constraints */
   int                   nconss4nlp;         /* number of varbound constraints we have to add to NLP because vbdvar is only implicit integer */
   int                   nconsslocal;        /* number of local varbound constraints */
   int                   i;                 

   assert(scip != NULL);
   assert(nlpi != NULL);
   assert(nlpiprob != NULL);
   assert(varsmap != NULL);
   assert(conssmap == NULL || nlpconsscounter != NULL);
   assert((explvbndconss == NULL) == (nexplvbndconss == NULL));
   assert(varbndconshdlr != NULL);
   
   nconss = SCIPconshdlrGetNConss(varbndconshdlr);
   if( !nconss )
      return SCIP_OKAY;
   
   conss = SCIPconshdlrGetConss(varbndconshdlr);
   
   /* count for how many global constraint the Vbdvar is only implicit integer, so we need to add constraint to NLP */
   nconss4nlp = 0;
   nconsslocal = 0;
   if( nexplvbndconss != NULL )
      *nexplvbndconss = 0;
   for( i = 0; i < nconss; ++i )
   {
      if( SCIPconsIsLocal(conss[i]) )
      {
         ++nconsslocal;
      }
      else if( explvbndconss == NULL || SCIPvarGetType(SCIPgetVbdvarVarbound(scip, conss[i])) > SCIP_VARTYPE_INTEGER )
      {  /* will add constraint to NLP soon */
         ++nconss4nlp;
      }
      else
      {  /* will handle constraint explicitely before solving NLP */
         ++*nexplvbndconss;
      }
   }

   if( nconsslocal == 0 && nconss4nlp != 0 )
   {
      SCIP_Real*       lhs;
      SCIP_Real*       rhs;
      const char**     consnames;
      int*             nlininds;
      int**            lininds;
      SCIP_Real**      linvals;
      int              j;
      int              k;

      /* number of and space for varbound constraints that are handled explicitely */
      assert(nconss == (nexplvbndconss ? *nexplvbndconss : 0) + nconss4nlp + nconsslocal);
      if( nexplvbndconss != NULL && *nexplvbndconss != 0 )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, explvbndconss, *nexplvbndconss) );
      }

      /* memory to store those varbound constraints that are added to the NLP as linear constraint */
      SCIP_CALL( SCIPallocBufferArray(scip, &lhs,      nconss4nlp) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rhs,      nconss4nlp) );
      SCIP_CALL( SCIPallocBufferArray(scip, &nlininds, nconss4nlp) );
      SCIP_CALL( SCIPallocBufferArray(scip, &lininds,  nconss4nlp) );
      SCIP_CALL( SCIPallocBufferArray(scip, &linvals,  nconss4nlp) );
      if( names )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &consnames, nconss4nlp) );
      }
      else
         consnames = NULL;
      
      j = 0; /* the number of explicitely handled varbound constraints passed so far */
      k = 0; /* the number of varbound constraints for the NLP passed so far */
      
      for( i = 0; i < nconss; ++i )
      {
         /* skip local constraints */
         if( SCIPconsIsLocal(conss[i]) )
            continue;
         
         /* treat constraint explicitly via boundchanges */
         if( explvbndconss != NULL && SCIPvarGetType(SCIPgetVbdvarVarbound(scip, conss[i])) <= SCIP_VARTYPE_INTEGER )
         {
            assert(j < *nexplvbndconss);
            (*explvbndconss)[j] = conss[i];
            ++j;
            continue;
         }
         /* else: add constraint to NLP */
         
         /* varbound constraints: lhs <= x + c * y <= rhs */
         lhs[k] = SCIPgetLhsVarbound(scip, conss[i]);
         rhs[k] = SCIPgetRhsVarbound(scip, conss[i]);
         if( names )
            consnames[k] = SCIPconsGetName(conss[i]);
         
         nlininds[k] = 2;
         SCIP_CALL( SCIPallocBufferArray(scip, &lininds[k], 2) );
         SCIP_CALL( SCIPallocBufferArray(scip, &linvals[k], 2) );
         
         /* add x term to coefficient matrix */
         linvals[k][0] = 1.0;
         assert(SCIPhashmapExists(varsmap, SCIPgetVarVarbound(scip, conss[i])));
         lininds[k][0] = (int) (size_t) SCIPhashmapGetImage(varsmap, SCIPgetVarVarbound(scip, conss[i]));
         
         /* add c * y term to coefficient matrix */
         linvals[k][1] = SCIPgetVbdcoefVarbound(scip, conss[i]);
         assert(SCIPhashmapExists(varsmap, SCIPgetVbdvarVarbound(scip, conss[i])));
         lininds[k][1] = (int) (size_t) SCIPhashmapGetImage(varsmap, SCIPgetVbdvarVarbound(scip, conss[i]));

         ++k;
         
         if( conssmap != NULL )
         {
            SCIP_CALL( SCIPhashmapInsert(conssmap, conss[i], (void*)(size_t)*nlpconsscounter) );
         }
         if( nlpconsscounter != NULL )
            ++*nlpconsscounter;
      }

      assert((nexplvbndconss == NULL && j == 0) || j == *nexplvbndconss);
      assert(k == nconss4nlp);
      
      SCIP_CALL( SCIPnlpiAddConstraints(nlpi, nlpiprob, nconss4nlp,
            lhs, rhs, nlininds, lininds, linvals,
            NULL, NULL, NULL, NULL, consnames) );
      
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
      SCIPfreeBufferArrayNull(scip, &consnames);
   }
   else
   {  /* all varbound constraints are global and handled explicitely */
      assert(explvbndconss  != NULL);
      assert(nexplvbndconss != NULL);
      SCIP_CALL( SCIPduplicateMemoryArray(scip, explvbndconss, conss, nconss) );
      *nexplvbndconss = nconss;
   }
   
   return SCIP_OKAY;
}

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
   SCIP_CALL( SCIPnlpiChgVarBounds(heurdata->nlpi, heurdata->nlpiprob, 1, &nlpidx, &lb, &ub) );

   return SCIP_OKAY;
}

/* sets up NLPI problem for heuristic */
static
SCIP_RETCODE setupNlp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur                /**< NLP heuristic */
   )
{
   int                   i;
   int                   cnt;                /* counter on discrete variables */
   SCIP_VAR*             var;
   SCIP_HEURDATA*        heurdata;

   assert(scip != NULL);
   assert(heur != NULL);
   
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->nlpi == NULL);
   assert(heurdata->var_nlp2scip == NULL);
   assert(heurdata->var_scip2nlp == NULL);

   if( SCIPgetNNlpis(scip) == 0 )
   {
      SCIPerrorMessage("No NLP solver available. Cannot setup NLP.\n");
      return SCIP_ERROR;
   }

   /* select an NLP solver; default selection is to take the solver with highest priority */
   assert(heurdata->nlpsolver != NULL);
   if( heurdata->nlpsolver[0] == '\0' )
   {
      heurdata->nlpi = SCIPgetNlpis(scip)[SCIPgetNNlpis(scip)-1];
   }
   else
   {
      heurdata->nlpi = SCIPfindNlpi(scip, heurdata->nlpsolver);
      if( heurdata->nlpi == NULL )
      {
         SCIPerrorMessage("NLP solver <%s> not found.\n", heurdata->nlpsolver);
         return SCIP_PARAMETERWRONGVAL;
      }
   }
   
   if( !heurdata->varboundexplicit )
   {
      heurdata->nvarbndconss = 0;
      heurdata->varbndconss = NULL;
   }
   
   SCIP_CALL( SCIPconstructNlpiProblemNlpHeur(scip, heur, heurdata->nlpi, &heurdata->nlpiprob,
      &heurdata->nvars, NULL, &heurdata->var_scip2nlp, NULL, 
      heurdata->varboundexplicit ? &heurdata->nvarbndconss : NULL,
      heurdata->varboundexplicit ? &heurdata->varbndconss : NULL,
      TRUE, heurdata->names) );
   assert(heurdata->nlpiprob != NULL);
      
   /* set some parameters of NLP solver */ 
   if( heurdata->nlptimelimit )
   {
      SCIP_CALL( SCIPnlpiSetRealPar(heurdata->nlpi, heurdata->nlpiprob, SCIP_NLPPAR_TILIM, heurdata->nlptimelimit) );
   }
   if( heurdata->nlpoptfile != NULL && *heurdata->nlpoptfile != '\0' )
   {
      SCIP_CALL( SCIPnlpiSetStringPar(heurdata->nlpi, heurdata->nlpiprob, SCIP_NLPPAR_OPTFILE, heurdata->nlpoptfile) );
   }

   /* collect and capture variables, collect discrete variables, collect bounds
    * assign variable indices to variables */
   heurdata->ndiscrvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->var_nlp2scip, heurdata->nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->discrvars, heurdata->ndiscrvars) );
   
   cnt = 0; /* counts number of discrete variables passed so far */
   for( i = 0; i < heurdata->nvars; ++i )
   {
      var = SCIPgetVars(scip)[i];
      
      heurdata->var_nlp2scip[i] = var;
      assert( (int)(size_t)SCIPhashmapGetImage(heurdata->var_scip2nlp, var) == i);
      
      SCIP_CALL( SCIPcaptureVar(scip, var) );

      if( SCIPvarGetType(var) <= SCIP_VARTYPE_INTEGER) /* binary or integer */
      {
         heurdata->discrvars[cnt] = i;
         ++cnt;
      }
      
      /** catch variable global bounds change events */
      SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_GBDCHANGED, heurdata->eventhdlr, (SCIP_EVENTDATA*)heurdata, NULL) );
   }
   
   /* catch explicit variable bound constraints */
   for( i = 0; i < heurdata->nvarbndconss; ++i )
   {
      SCIP_CALL( SCIPcaptureCons(scip, heurdata->varbndconss[i]) );
   }
   
   /* initialize data structure for NLP solve statistics */
   SCIP_CALL( SCIPnlpStatisticsCreate(&heurdata->nlpstatistics) );

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
   int            i;

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
      for( i = 0; !havenlp && i < heurdata->nnlpiinits; ++i )
      {
         SCIP_CALL( (*heurdata->haveconss[i])(scip, TRUE, &havenlp) );
      }
   }
   
   if( !havenlp )
   { /* no suitable nonlinearities -> forget about NLP and return */
      SCIPdebugMessage("No nonlinear continuous variables. NLP local search heuristic will not run.\n");
      return SCIP_OKAY;
   }

   /* setup NLP (if NLP solver available) */
   if( SCIPgetNNlpis(scip) > 0 )
   {
      SCIP_CALL( setupNlp(scip, heur) );
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "No NLP solver available. NLP local search heuristic will not run.\n");
      havenlp = FALSE;
   }
   
   return SCIP_OKAY;   
}

/** for a fixation of discrete variables, applies the variable bound constraints to the NLP */
static
SCIP_RETCODE applyVarBoundConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_SOL*             refpoint            /**< point to take fixation of discrete variables from */
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
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scip), SCIPcalcHashtableSize(5 * heurdata->nvarbndconss)) );
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
       
         shift = SCIPgetVbdcoefVarbound(scip, cons) * SCIPgetSolVal(scip, refpoint, SCIPgetVbdvarVarbound(scip, cons));
         if( !SCIPisInfinity(scip, -varlb[varcnt]) )
            varlb[varcnt] -= shift;
         if( !SCIPisInfinity(scip,  varub[varcnt]) )
            varub[varcnt] -= shift;
         
         if( SCIPvarGetLbGlobal(var) > varlb[varcnt] )
            varlb[varcnt] = SCIPvarGetLbGlobal(var);
         if( SCIPvarGetUbGlobal(var) < varub[varcnt] )
            varub[varcnt] = SCIPvarGetUbGlobal(var);
         if( varlb[varcnt] > varub[varcnt] )
            varlb[varcnt] = varub[varcnt];
         
         /* remember that a bound change for variable var is now stored at position varbnd */
         SCIP_CALL( SCIPhashmapInsert(varmap, var, (void*)(size_t)(varcnt+1)) );
#if 0
         SCIPdebugMessage("%s: var %s at %d now bounded in [%g, %g] due to %s = %g\n",
            SCIPconsGetName(cons), SCIPvarGetName(var), varidx[varcnt],
            varlb[varcnt], varub[varcnt], SCIPvarGetName(SCIPgetVbdvarVarbound(scip, cons)),
            SCIPgetSolVal(scip, refpoint, SCIPgetVbdvarVarbound(scip, cons)) );
#endif
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

         shift = SCIPgetVbdcoefVarbound(scip, cons) * SCIPgetSolVal(scip, refpoint, SCIPgetVbdvarVarbound(scip, cons));
         if( !SCIPisInfinity(scip, -lhs) )
            lhs -= shift;
         if( !SCIPisInfinity(scip,  rhs) )
            rhs -= shift;

         /* possibly tighten previously stored bound change on variable var with newly computed bounds */ 
         varlb[idx_] = MAX(varlb[idx_],lhs);
         varub[idx_] = MIN(varub[idx_],rhs);
         if( varlb[idx_] > varub[idx_] )
            varlb[idx_] = varub[idx_];
#if 0
         SCIPdebugMessage("%s: var %s at %d now bounded in [%g, %g] due to %s = %g  [updated]\n",
            SCIPconsGetName(cons), SCIPvarGetName(var), varidx[idx_],
            varlb[idx_], varub[idx_], SCIPvarGetName(SCIPgetVbdvarVarbound(scip, cons)),
            SCIPgetSolVal(scip, refpoint, SCIPgetVbdvarVarbound(scip, cons)) );
#endif
      }
   }
   
   /* apply bound changes on variables in varbound constraint to NLP */
   SCIP_CALL( SCIPnlpiChgVarBounds(heurdata->nlpi, heurdata->nlpiprob, varcnt, varidx, varlb, varub) );

   SCIPhashmapFree(&varmap);   
   SCIPfreeBufferArray(scip, &varidx);
   SCIPfreeBufferArray(scip, &varub);
   SCIPfreeBufferArray(scip, &varlb);

   return SCIP_OKAY;
}

/** for a fixation of discrete variables, applies the indicator constraints to the NLP
 * if binaries of indicator constraints are fixed to 1, fixes the corresponding slack variable to 0 */
static
SCIP_RETCODE applyIndicatorConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_SOL*             refpoint            /**< point to take fixation of discrete variables from */
   )
{
   SCIP_CONS* cons;
   int        nconss;
   int        c;
   int        cnt;
   int*       varidx;
   SCIP_Real* varlb;
   SCIP_Real* varub;
   SCIP_VAR*  var;
   SCIP_Real  scalar;
   SCIP_Real  offset;
   
   assert(scip     != NULL);
   assert(heurdata != NULL);
   
   if( heurdata->conshdlrindicator == NULL )
      return SCIP_OKAY;
   
   nconss = SCIPconshdlrGetNConss(heurdata->conshdlrindicator);
   if( nconss == 0 )
      return SCIP_OKAY;
   
   SCIP_CALL( SCIPallocBufferArray(scip, &varidx, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varlb,  nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varub,  nconss) );

   cnt = 0;
   for( c = 0; c < nconss; ++c )
   {
      cons = SCIPconshdlrGetConss(heurdata->conshdlrindicator)[c];
      assert(cons != NULL);
      
      var = SCIPgetSlackVarIndicator(cons);
      scalar = 1.0;
      offset = 0.0;
      SCIP_CALL( SCIPvarGetProbvarSum(&var, &scalar, &offset) );

      varidx[cnt] = (int)(size_t)SCIPhashmapGetImage(heurdata->var_scip2nlp, var);
      /* if slack variable is multiaggregated or fixed, then we skip this propagation */
      if( varidx[cnt] < 0 || scalar == 0.0 )
         continue;

      /* if the binary variable indicates that the slack should be 0.0, set it this way,
       * otherwise slack is allowed to vary within global bounds (resets possible fixings from previous call to heuristic)
       */
      if( SCIPgetSolVal(scip, refpoint, SCIPgetBinaryVarIndicator(cons)) > 0.5 )
      {
         varlb[cnt] = -offset / scalar;
         varub[cnt] = varlb[cnt];

         SCIPdebugMessage("%s: var <%s> at %d now fixed to 0.0 due to <%s> = %g\n",
            SCIPconsGetName(cons), SCIPvarGetName(SCIPgetSlackVarIndicator(cons)), varidx[cnt],
            SCIPvarGetName(SCIPgetBinaryVarIndicator(cons)),
            SCIPgetSolVal(scip, refpoint, SCIPgetBinaryVarIndicator(cons)) );
      }
      else
      {
         varlb[cnt] = (SCIPvarGetLbGlobal(SCIPgetSlackVarIndicator(cons)) - offset) / scalar;
         varub[cnt] = (SCIPvarGetUbGlobal(SCIPgetSlackVarIndicator(cons)) - offset) / scalar;
      }
      ++cnt;
   }

   /* apply bound changes on variables in varbound constraint to NLP */
   SCIP_CALL( SCIPnlpiChgVarBounds(heurdata->nlpi, heurdata->nlpiprob, cnt, varidx, varlb, varub) );

   SCIPfreeBufferArray(scip, &varidx);
   SCIPfreeBufferArray(scip, &varlb);
   SCIPfreeBufferArray(scip, &varub);

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

   SCIP_CALL( SCIPnlpiFreeProblem(heurdata->nlpi, &heurdata->nlpiprob) );
   assert(heurdata->nlpiprob == NULL);

   heurdata->nlpi = NULL;
   
   SCIPnlpStatisticsFree(&heurdata->nlpstatistics);
   assert(heurdata->nlpstatistics == NULL);

   if( heurdata->nvarbndconss )
   {
      assert(heurdata->varbndconss != NULL);
      for( i = 0; i < heurdata->nvarbndconss; ++i )
      {
         assert(heurdata->varbndconss[i] != NULL);
         SCIP_CALL( SCIPreleaseCons(scip, &heurdata->varbndconss[i]) );
      }

      SCIPfreeMemoryArray(scip, &heurdata->varbndconss);
      heurdata->nvarbndconss = 0;
   }

   return SCIP_OKAY;
}

/** sets up NLP from constraints in SCIP */
SCIP_RETCODE SCIPconstructNlpiProblemNlpHeur(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_NLPI*            nlpi,               /**< NLP solver interface for which the NLP should be setup */
   SCIP_NLPIPROBLEM**    nlpiproblem,        /**< buffer to store pointer to new NLPI problem */
   int*                  nnlpvars,           /**< buffer to store number of variables in NLPI problem (should be SCIPgetNVars(scip)), or NULL */
   int*                  nnlpconss,          /**< buffer to store number of constraints in NLPI problem, or NULL */
   SCIP_HASHMAP**        varsmap,            /**< buffer to store pointer to mapping from SCIP variables to variable indicies in NLPI_PROBLEM, or NULL if not needed */
   SCIP_HASHMAP**        conssmap,           /**< buffer to store pointer to mapping from SCIP constraints to row indices in NLPI_PROBLEM, or NULL if not needed */
   int*                  nexplvbdconss,      /**< buffer to store number of variable bound constraints that should be handled explicitely by the user, or NULL if all should be added to NLPI (makes only sense if onlysubnlp == TRUE) */
   SCIP_CONS***          explvbdconss,       /**< buffer to store array of variable bound constraint that should be handled explicitely by the user, or NULL if all should be added to NLPI (makes only sense if onlysubnlp == TRUE) */
   SCIP_Bool             onlysubnlp,         /**< whether to add only constraints that are relevant for the NLP obtained by fixing all discrete variables in the CIP */
   SCIP_Bool             names               /**< whether to gives variable and constraint names to NLPI */
   )
{
   SCIP_HEURDATA*        heurdata;
   int                   nvars;
   SCIP_VAR**            vars;
   SCIP_Real*            varlb;
   SCIP_Real*            varub;
   const char**          varnames;
   SCIP_HASHMAP*         var_scip2nlp;
   SCIP_Real*            objcoeff;
   int*                  objvar;
   int                   i;
   int                   cnt;                /* counter on discrete variables */
   SCIP_CONSHDLR*        conshdlr;
   int                   mynnlpconss;
   
   assert(scip != NULL);
   assert(nlpi != NULL);
   assert(nlpiproblem != NULL);
   assert(nexplvbdconss == NULL || onlysubnlp); /* if nvbdconss != NULL, then onlysubnlp should be TRUE */
   assert(explvbdconss == NULL  || onlysubnlp); /* if  vbdconss != NULL, then onlysubnlp should be TRUE */

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* get SCIP variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* create an NLPI problem */
   SCIP_CALL( SCIPnlpiCreateProblem(nlpi, nlpiproblem, "nlp") );
   assert(*nlpiproblem != NULL);

   /* set value for infinity */
   SCIP_CALL( SCIPnlpiSetRealPar(nlpi, *nlpiproblem, SCIP_NLPPAR_INFINITY, SCIPinfinity(scip)) );

   /* init mapping from SCIP variables to NLPI indices */
   SCIP_CALL( SCIPhashmapCreate(&var_scip2nlp, SCIPblkmem(scip), SCIPcalcHashtableSize(5 * nvars)) );
   
   /* alloc space for variable lower bound, upper bounds, and names */
   SCIP_CALL( SCIPallocBufferArray(scip, &varlb, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varub, nvars) );
   if( names )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &varnames, nvars) );
   }
   else
   {
      varnames = NULL;
   }
   
   /* collect variable lower bound, upper bounds, and names */
   for( i = 0; i < nvars; ++i )
   {
      varlb[i] = SCIPvarGetLbGlobal(vars[i]);
      varub[i] = SCIPvarGetUbGlobal(vars[i]);
      if( varnames != NULL )
         varnames[i] = SCIPvarGetName(vars[i]);
      
      SCIP_CALL( SCIPhashmapInsert(var_scip2nlp, vars[i], (void*)(size_t)i) );
   }
   
   /* add variables to NLP solver */
   SCIP_CALL( SCIPnlpiAddVars(nlpi, *nlpiproblem, nvars, varlb, varub, varnames) );
   
   SCIPfreeBufferArray(scip, &varub);
   SCIPfreeBufferArray(scip, &varlb);
   SCIPfreeBufferArrayNull(scip, &varnames);

   /* collect objective coefficients for minimization objective */
   SCIP_CALL( SCIPallocBufferArray(scip, &objcoeff, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &objvar, nvars) );

   cnt = 0; /* number of nonzeros in objective passed so far */
   for( i = 0; i < nvars; ++i )
   {
      if( SCIPvarGetObj(vars[i]) )
      {
         /* the transformed problem in SCIP is always a minimization problem,
          * thus we do not have to take care of the objective sense here */
         assert( !SCIPvarIsOriginal(vars[i]) );
         objcoeff[cnt] = SCIPvarGetObj(vars[i]);
         objvar[cnt]   = i;
         ++cnt;
      }
   }
   /* add (linear) objective to NLP solver */
   SCIP_CALL( SCIPnlpiSetObjective(nlpi, *nlpiproblem, cnt, objvar, objcoeff, 0, NULL, NULL, NULL, SCIPgetTransObjoffset(scip)) );
   
   SCIPfreeBufferArray(scip, &objvar);
   SCIPfreeBufferArray(scip, &objcoeff);
   
   /* initialize hashmap for constraints if user wants to have one */
   if( conssmap != NULL )
   {
      SCIP_CALL( SCIPhashmapCreate(conssmap, SCIPblkmem(scip), SCIPcalcHashtableSize(5 * SCIPgetNConss(scip))) );
   }
   mynnlpconss = 0;
   
   /* add linear constraints to NLP solver */
   conshdlr = SCIPfindConshdlr(scip, "linear");
   if( conshdlr != NULL )
   {
      SCIP_CALL( addLinearConstraints(scip, nlpi, *nlpiproblem, var_scip2nlp, conssmap ? *conssmap : NULL, &mynnlpconss, conshdlr, onlysubnlp, names) );
   }

   /* add varbound constraints to NLP solver */
   conshdlr = SCIPfindConshdlr(scip, "varbound");
   if( conshdlr != NULL )
   {
      SCIP_CALL( collectVarBoundConstraints(scip, nlpi, *nlpiproblem, var_scip2nlp, conssmap ? *conssmap : NULL, &mynnlpconss, nexplvbdconss, explvbdconss, conshdlr, names) );
   }
   
   if( !onlysubnlp )
   {
      /* add also logicor constraints to NLP */
      conshdlr = SCIPfindConshdlr(scip, "logicor");
      if( conshdlr != NULL )
      {
         SCIP_CALL( addLogicOrConstraints(scip, nlpi, *nlpiproblem, var_scip2nlp, conssmap ? *conssmap : NULL, &mynnlpconss, conshdlr, names) );
      }

      /* add also setppc constraints to NLP */
      conshdlr = SCIPfindConshdlr(scip, "setppc");
      if( conshdlr != NULL )
      {
         SCIP_CALL( addSetppcConstraints(scip, nlpi, *nlpiproblem, var_scip2nlp, conssmap ? *conssmap : NULL, &mynnlpconss, conshdlr, names) );
      }

      /* add also knapsack constraints to NLP solver */
      conshdlr = SCIPfindConshdlr(scip, "knapsack");
      if( conshdlr != NULL )
      {
         SCIP_CALL( addKnapsackConstraints(scip, nlpi, *nlpiproblem, var_scip2nlp, conssmap ? *conssmap : NULL, &mynnlpconss, conshdlr, names) );
      }
   }
   
   /** call user given NLPI initialization functions */
   for( i = 0; i < heurdata->nnlpiinits; ++i )
   {
      SCIP_CALL( (*heurdata->nlpiinits[i])(scip, nlpi, *nlpiproblem, var_scip2nlp, conssmap ? *conssmap : NULL, &mynnlpconss, onlysubnlp, names) );
   }
   
   if( nnlpvars != NULL )
      *nnlpvars = nvars;
   
   if( varsmap != NULL )
      *varsmap = var_scip2nlp;
   else
      SCIPhashmapFree(&var_scip2nlp);
   
   if( nnlpconss != NULL )
      *nnlpconss = mynnlpconss;
   
   return SCIP_OKAY;
}

/** main procedure of the NLP heuristic */
SCIP_RETCODE SCIPapplyNlpHeur(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                       */
   SCIP_RESULT*          result,             /**< result data structure                                          */
   SCIP_SOL*             refpoint,           /**< point to take fixation of discrete variables from, and startpoint for NLP solver; if NULL, then LP solution is used */
   SCIP_Longint          itercontingent,     /**< iteration limit for NLP solver, or -1 for default of NLP heuristic */
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
               SCIPfreeBufferArray(scip, &discrfix);
               SCIPfreeBufferArray(scip, &startpoint);
               if( SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
               {
                  SCIP_CALL( destroyNLP(scip, heurdata) );
                  heurdata->triedsetupnlp = FALSE;
               }
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
      SCIP_CALL( SCIPnlpiChgVarBounds(heurdata->nlpi, heurdata->nlpiprob, heurdata->ndiscrvars, heurdata->discrvars, discrfix, discrfix) );

      SCIPfreeBufferArray(scip, &discrfix);
   }
   /* apply those variable bound constraints that we can apply explicitely */
   SCIP_CALL( applyVarBoundConstraints(scip, heurdata, refpoint) );

   /* apply indicator constraints */
   SCIP_CALL( applyIndicatorConstraints(scip, heurdata, refpoint) );

   /* set time and iteration limit for NLP solver */
   if( itercontingent == -1 && heurdata->nlpiterlimit > 0 )
      itercontingent = heurdata->nlpiterlimit;

   if( itercontingent > 0 )
   {
      SCIP_CALL( SCIPnlpiSetIntPar(heurdata->nlpi, heurdata->nlpiprob, SCIP_NLPPAR_ITLIM, (int)itercontingent) );
   }

   SCIP_CALL( SCIPnlpiSetRealPar(heurdata->nlpi, heurdata->nlpiprob, SCIP_NLPPAR_TILIM, timelimit) );

   SCIP_CALL( SCIPnlpiSetIntPar(heurdata->nlpi, heurdata->nlpiprob, SCIP_NLPPAR_VERBLEVEL, heurdata->nlpverblevel) );

   /* pass initial guess to NLP solver, if we have one; otherwise clear previous guess and let NLP solver choose */
   if( refpoint || (SCIPgetStage(scip) == SCIP_STAGE_SOLVING && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL) )
   {
      SCIP_CALL( SCIPnlpiSetInitialGuess(heurdata->nlpi, heurdata->nlpiprob, startpoint) );
   }
   else
   {
      SCIP_CALL( SCIPnlpiSetInitialGuess(heurdata->nlpi, heurdata->nlpiprob, NULL) );
   }
   
   SCIPfreeBufferArray(scip, &startpoint);

   *result = SCIP_DIDNOTFIND;

   /* let the NLP solver do its magic */
   SCIPdebugMessage("start NLP solve with iteration limit %"SCIP_LONGINT_FORMAT" and timelimit %g\n", itercontingent, timelimit);
   SCIP_CALL( SCIPnlpiSolve(heurdata->nlpi, heurdata->nlpiprob) );

   SCIPdebugMessage("NLP solver returned with termination status %d and solution status %d\n",
      SCIPnlpiGetTermstat(heurdata->nlpi, heurdata->nlpiprob), SCIPnlpiGetSolstat(heurdata->nlpi, heurdata->nlpiprob));
   
   if( SCIPnlpiGetTermstat(heurdata->nlpi, heurdata->nlpiprob) >= SCIP_NLPTERMSTAT_MEMERR )
   {  /* oops, something did not go well at all */
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, 
         "NLP solver returned with bad termination status %d. Will not run NLP heuristic again for this run.\n",  
         SCIPnlpiGetTermstat(heurdata->nlpi, heurdata->nlpiprob));
      SCIP_CALL( destroyNLP(scip, heurdata) );
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPnlpiGetStatistics(heurdata->nlpi, heurdata->nlpiprob, heurdata->nlpstatistics) );

   if( iterused != NULL )
      *iterused += SCIPnlpStatisticsGetNIterations(heurdata->nlpstatistics);
   SCIPdebugMessage("NLP solver used %d iterations and %g seconds\n", 
      SCIPnlpStatisticsGetNIterations(heurdata->nlpstatistics), SCIPnlpStatisticsGetTotalTime(heurdata->nlpstatistics));

   if( SCIPnlpiGetSolstat(heurdata->nlpi, heurdata->nlpiprob) <= SCIP_NLPSOLSTAT_FEASIBLE )
   {  /* NLP solver claims, it found a feasible (maybe even optimal) solution */
      SCIP_Real* primals;
      SCIP_SOL*  sol;
      SCIP_Bool  stored;
      SCIP_Bool  feasible;

      /* get solution from NLP solver and create a SCIP_SOL out of it */
      SCIP_CALL( SCIPnlpiGetSolution(heurdata->nlpi, heurdata->nlpiprob, &primals) );
      assert(primals != NULL);

      SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
      SCIP_CALL( SCIPsetSolVals(scip, sol, heurdata->nvars, heurdata->var_nlp2scip, primals) );

      SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, TRUE, TRUE, &feasible) );

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
         SCIP_CALL( SCIPnlpiSetRealPar(heurdata->nlpi, heurdata->nlpiprob, SCIP_NLPPAR_FEASTOL, heurdata->resolvetolfactor*SCIPfeastol(scip)) );

         /* set almost-feasible solution as starting point for NLP solver, if user allows us */
         if( !heurdata->resolvefromscratch )
         {
            SCIP_CALL( SCIPnlpiSetInitialGuess(heurdata->nlpi, heurdata->nlpiprob, primals) );
         }

         /* solve again */
         SCIP_CALL( SCIPnlpiSolve(heurdata->nlpi, heurdata->nlpiprob) );
         SCIPdebugMessage("NLP solver returned with termination status %d and solution status %d\n",
            SCIPnlpiGetTermstat(heurdata->nlpi, heurdata->nlpiprob), SCIPnlpiGetSolstat(heurdata->nlpi, heurdata->nlpiprob));
         SCIP_CALL( SCIPnlpiGetStatistics(heurdata->nlpi, heurdata->nlpiprob, heurdata->nlpstatistics) );
         if( iterused != NULL )
            *iterused += SCIPnlpStatisticsGetNIterations(heurdata->nlpstatistics);

         if( SCIPnlpiGetSolstat(heurdata->nlpi, heurdata->nlpiprob) <= SCIP_NLPSOLSTAT_FEASIBLE )
         { 
            /* NLP solver claims feasible again, hope that SCIP accepts it now */
            SCIP_CALL( SCIPnlpiGetSolution(heurdata->nlpi, heurdata->nlpiprob, &primals) );
            assert(primals != NULL);

            SCIP_CALL( SCIPsetSolVals(scip, sol, heurdata->nvars, heurdata->var_nlp2scip, primals) );

#ifdef SCIP_DEBUG
            /* print the infeasibilities to stdout */
            SCIP_CALL( SCIPcheckSolOrig(scip, sol, &feasible, TRUE, TRUE) );
#endif
            SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, TRUE, &stored) );
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
         SCIP_CALL( SCIPnlpiSetRealPar(heurdata->nlpi, heurdata->nlpiprob, SCIP_NLPPAR_FEASTOL, SCIPfeastol(scip)) );
      }

      SCIP_CALL( SCIPfreeSol(scip, &sol) );
   }

   /* if the heuristic was applied before solving has started, then destroy NLP, since EXITSOL may not be called */
   if( SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
   {
      SCIP_CALL( destroyNLP(scip, heurdata) );
      heurdata->triedsetupnlp = FALSE;
   }

   /* TODO: reset time and iterlimit in nlp solver? */

   return SCIP_OKAY;   
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyNlp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurNlp(scip) );
 
   return SCIP_OKAY;
}

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
   
   SCIPfreeMemoryArrayNull(scip, &heurdata->nlpiinits);
   SCIPfreeMemoryArrayNull(scip, &heurdata->haveconss);

   SCIPfreeMemoryNull(scip, &heurdata);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitNlp)
{
   SCIP_HEURDATA* heurdata;
   
   assert(scip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   
   heurdata->conshdlrindicator = SCIPfindConshdlr(scip, "indicator");

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitNlp)
{
   SCIP_HEURDATA* heurdata;
   
   assert(scip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->startcand == NULL);
   
   heurdata->conshdlrindicator = NULL;

   return SCIP_OKAY;
}

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

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->nlpi == NULL);

   /* try to setup NLP */
   SCIP_CALL( checkCIPandSetupNLP(scip, heur) );
   
   /* it's ok if NLP was not setup (e.g., because there are no nonlinear continuous variables or there is no NLP solver) */
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

   if( heurdata->startcand != NULL )
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
            SCIPdebugMessage("NLP heuristic delayed because no start candidate given and no LP solution available; LP status = %d\n", SCIPgetLPSolstat(scip));
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

   if( !heurdata->runalways )
   {
      /* check if enough nodes have been processed so that we want to run the heuristic again */

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
   }
   else
   {
      itercontingent = -1;
   }

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

      SCIP_CALL( SCIPnlpiSetRealPar(heurdata->nlpi, heurdata->nlpiprob, SCIP_NLPPAR_TILIM, timelimit) );
   }
   else if( heurdata->nlptimelimit > 0 )
   {  /* enforce user given time limit */
      SCIP_CALL( SCIPnlpiSetRealPar(heurdata->nlpi, heurdata->nlpiprob, SCIP_NLPPAR_TILIM, heurdata->nlptimelimit) );
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
         NULL,
         NULL, NULL, NULL, NULL, NULL, NULL, processVarEvent, NULL) );
   heurdata->eventhdlr = SCIPfindEventhdlr(scip, HEUR_NAME);
   
   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, 
         heurCopyNlp,
         heurFreeNlp, heurInitNlp, heurExitNlp, heurInitsolNlp, heurExitsolNlp, heurExecNlp,
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

   SCIP_CALL( SCIPaddStringParam(scip, "heuristics/"HEUR_NAME"/nlpoptfile",
         "name of an NLP solver specific options file",
         &heurdata->nlpoptfile, TRUE, "", NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, "heuristics/"HEUR_NAME"/nlpsolver",
         "name of an NLP solver to use (empty value means to use solver with highest priority)",
         &heurdata->nlpsolver, FALSE, "", NULL, NULL) );

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

   SCIP_CALL( SCIPaddBoolParam (scip, "heuristics/"HEUR_NAME"/runalways",
         "whether to run NLP heuristic always if starting point available (does not use iteroffset,iterquot,itermin)",
         &heurdata->runalways, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/varboundexplicit",
         "should variable bound constraints be handled explicitly before solving the NLP instead of adding them to the NLP?",
         &heurdata->varboundexplicit, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/names",
         "should variable and constraint names passed to the NLP solver?",
         &heurdata->names, FALSE, FALSE, NULL, NULL) );
   
   return SCIP_OKAY;
}

/** includes an NLPI initialization method into the NLP heuristic
 * can be used by constraint handlers to register a function that inserts their constraints into an NLPI */
SCIP_RETCODE SCIPincludeHeurNlpNlpiInit(
   SCIP*                   scip,               /**< SCIP data structure */
   SCIP_DECL_HEURNLPHAVECONS((*havecons)),     /**< method to call for checking if potential constraints for the NLP are present */
   SCIP_DECL_HEURNLPNLPIINIT((*nlpiinit)),     /**< method to call for initializing NLP */
   const char*             conshdlrname        /**< name of the constraint handler */
   )
{
   SCIP_HEUR*     heur;
   SCIP_HEURDATA* heurdata;
   int            i;
   
   assert(scip != NULL);
   assert(havecons != NULL);
   assert(nlpiinit != NULL);
   assert(conshdlrname != NULL );
   
   /* find the NLP heuristic handler */
   heur = SCIPfindHeur(scip, HEUR_NAME);
   if( heur == NULL )
   {
      SCIPerrorMessage("NLP heuristic not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   for( i = heurdata->nnlpiinits - 1; i >= 0; --i )
   {
      if( heurdata->nlpiinits[i] == nlpiinit )
      {
#ifdef SCIP_DEBUG
         SCIPwarningMessage("Try to add already known initialization method %p for constraint handler <%s>.\n", nlpiinit, conshdlrname);
#endif
         return SCIP_OKAY;
      }
   }

   /* insert NLPI initialization method into heuristic data */
   assert(heurdata->nnlpiinits <= heurdata->nlpiinitssize);
   if( heurdata->nnlpiinits+1 > heurdata->nlpiinitssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, heurdata->nnlpiinits+1);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &heurdata->nlpiinits, newsize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &heurdata->haveconss, newsize) );
      heurdata->nlpiinitssize = newsize;
   }
   assert(heurdata->nnlpiinits+1 <= heurdata->nlpiinitssize);
   
   heurdata->nlpiinits[heurdata->nnlpiinits] = nlpiinit;
   heurdata->haveconss[heurdata->nnlpiinits] = havecons;
   heurdata->nnlpiinits++;

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

   /* too early or the game is over already: no more interest in starting points */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;
   
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   
   /* we do not have a NLP, so we also do not need a starting point */
   if( heurdata->nlpi == NULL )
      return SCIP_OKAY;
   
   SCIPdebugMessage("consider solution candidate with violation %g and objective %g\n", violation, SCIPgetSolTransObj(scip, solcand));

   /* if we have no point yet, or the new point has a lower constraint violation, or it has a better objective function value,
    * then take the new point */
   if( heurdata->startcand == NULL || SCIPisGT(scip, heurdata->startcandviol, violation) ||
      SCIPisRelGT(scip, SCIPgetSolTransObj(scip, heurdata->startcand), SCIPgetSolTransObj(scip, solcand)) )
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

/** gets NLPI interface used by NLP heuristic, or NULL if none */
SCIP_NLPI* SCIPgetNlpiHeurNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur                /**< heuristic data structure                                       */
   )
{
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   return heurdata->nlpi;
}

/** gets the NLP that is used by the NLP heuristic as NLPI problem, or NULL if not constructed */
SCIP_NLPIPROBLEM* SCIPgetNlpiProblemHeurNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur                /**< heuristic data structure                                       */
   )
{
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   return heurdata->nlpiprob;
}

/** gets number of variables in NLP
 * usually, this equals SCIPgetNVars()
 */
int SCIPgetGetNVarsHeurNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur                /**< heuristic data structure                                       */
   )
{
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   return heurdata->nvars;
}

/** gets mapping of NLP variables to SCIP variables */
SCIP_VAR** SCIPgetNlpVarsHeurNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur                /**< heuristic data structure                                       */
   )
{
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   return heurdata->var_nlp2scip;
}

/** gets mapping of SCIP variables to NLP variables */
SCIP_HASHMAP* SCIPgetVarMappingHeurNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur                /**< heuristic data structure                                       */
   )
{
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   return heurdata->var_scip2nlp;
}

/** gets nlpi initialization callbacks */
void SCIPgetNlpiSetupDataHeurNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                       */
   SCIP_DECL_HEURNLPNLPIINIT((***nlpiinits)),/**< buffer to store pointer to array of NLPIINIT callbacks, or NULL */
   SCIP_DECL_HEURNLPHAVECONS((***haveconss)),/**< buffer to store pointer to array of HAVECONS callbacks, or NULL */
   int*                  nnlpiinits          /**< buffer to store number of NLPIINIT and HAVECONS callbacks, or NULL */
   )
{
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   if( nlpiinits != NULL )
      *nlpiinits = heurdata->nlpiinits;
   if( haveconss != NULL )
      *haveconss = heurdata->haveconss;
   if( nnlpiinits != NULL )
      *nnlpiinits = heurdata->nnlpiinits;
}

/** gets current start candidate of heuristic, or NULL if none;
 * if forget is TRUE and return is non-NULL, then its the obligation of the user to free the returned solution */
SCIP_SOL* SCIPgetStartcandHeurNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                       */
   SCIP_Bool             forget              /**< should the heuristic forget this starting candidate?           */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_SOL*      startcand;

   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   
   startcand = heurdata->startcand;
   
   if( forget && startcand != NULL )
   {
      heurdata->startcand = NULL;
   }
   
   return startcand;
}

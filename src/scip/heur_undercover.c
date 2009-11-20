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
#pragma ident "@(#) $Id: heur_undercover.c,v 1.9 2009/11/20 17:02:28 bzfgleix Exp $"

/**@file   heur_undercover.c
 * @ingroup PRIMALHEURISTICS
 * @brief  undercover primal heuristic for MIQCPs

 * @author Ambros Gleixner
 *
 * @todo: implement ppc strategy constraint 'v'iolation and 'b'ranching status
 * @todo: parameter tuning for subproblem and ppc problem
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>


#include "scip/heur_undercover.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_setppc.h"
#include "scip/cons_soc.h"
#ifdef WITH_UNIVARDEFINITE
#include "cons_univardefinite.h"
#endif
#ifdef WITH_CONSBRANCHNL
#include "cons_branchnonlinear.h"
#endif
#include "scip/heur_nlp.h"

#define HEUR_NAME             "undercover"
#define HEUR_DESC             "solves a linearization of an MIQCP determined by a set covering approach"
#define HEUR_DISPCHAR         'u'
#define HEUR_PRIORITY         -1110000
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE

#define DEFAULT_MINNODES      (SCIP_Longint)500/* minimum number of nodes to regard in the subproblem                 */
#define DEFAULT_MAXNODES      (SCIP_Longint)5000/* maximum number of nodes to regard in the subproblem                */
#define DEFAULT_MINIMPROVE    0.01           /* factor by which heuristic should at least improve the incumbent       */
#define DEFAULT_NODESOFS      (SCIP_Longint)500/* number of nodes added to the contingent of the total nodes          */
#define DEFAULT_NODESQUOT     0.1            /* subproblem nodes in relation to nodes of the original problem         */
#define DEFAULT_PPCOBJQUOT    1.0            /* additional penalty factor for fixing continuous variables             */
#define DEFAULT_DOMRED        1.0            /* reduce domain of selected variables by this factor around LP value    */
#define DEFAULT_LOCKSROUNDING TRUE           /* shall LP values for integer vars be rounded according to locks?       */
#define DEFAULT_ONLYCONVEXIFY FALSE          /* should we only fix/dom.red. variables creating nonconvexity?          */
#define DEFAULT_ONLYATINTEGER FALSE          /* should heuristic be called only at integer feasible nodes?            */
#define DEFAULT_POSTNLP       TRUE           /* should the NLP heuristic be called to polish a feasible solution?     */
#define DEFAULT_PPCSTRAT      'u'            /* strategy for finding a ppc solution                                   */
#define PPCSTRATS             "bcdlmtuv"     /* strategies for finding a ppc solution                                 */


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem                 */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem                 */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes          */
   SCIP_Longint          nusednodes;         /**< nodes already used by heuristic in earlier calls                    */
   SCIP_Real             minimprove;         /**< factor by which heuristic should at least improve the incumbent     */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem       */
   SCIP_Real             ppcobjquot;         /**< additional penalty factor for fixing continuous variables           */
   SCIP_Real             domred;             /**< reduce domain of selected variables by this factor around LP value  */
   SCIP_Bool             locksrounding;      /**< shall LP values for integer vars be rounded according to locks?     */
   SCIP_Bool             onlyconvexify;      /**< should we only fix/dom.red. variables creating nonconvexity?        */
   SCIP_Bool             onlyatinteger;      /**< should heuristic be called only at integer feasible nodes?          */
   SCIP_Bool             postnlp;            /**< should the NLP heuristic be called to polish a feasible solution?   */
   SCIP_Bool             run;                /**< should heuristic run, i.e. are nonlinear constraints present?       */
   char                  ppcstrat;           /**< strategy for finding a ppc solution                                 */
};




/*
 * Local methods
 */

/** creates a set covering/packing problem to determine a number of variables to be fixed */
static
SCIP_RETCODE createPpcProblem(
   SCIP*                 scip,               /**< original SCIP data structure                                        */
   SCIP*                 ppcscip,            /**< SCIP data structure for the set covering problem                    */
   SCIP_VAR**            ppcvars,            /**< the variables of the set covering problem                           */
   char                  strat,              /**< strategy for finding a ppc solution                                 */
   SCIP_Real             objquot,            /**< additional penalty factor for fixing continuous variables           */
   SCIP_Bool             local,              /**< shall locally fixed variables be treated as constants               */
   SCIP_Bool             onlyconvexify,      /**< should we only fix/dom.red. variables creating nonconvexity?        */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully       */
   )
{
   SCIP_VAR** vars;
   SCIP_CONSHDLR* conshdlr;

   SCIP_Real* conscounter;
   SCIP_Real* termcounter;
   SCIP_Real conobjfac;
   SCIP_Real intobjfac;


   int nvars;
   int ncontvars;
   int i;

   char name[SCIP_MAXSTRLEN];

   assert( scip != NULL );
   assert( ppcscip != NULL );
   assert( ppcvars != NULL );
   assert( strchr(PPCSTRATS, strat) != NULL );
   assert( SCIPisFeasPositive(scip, objquot) && SCIPisFeasPositive(scip, 1.0/objquot) );
   assert( success != NULL );

   *success = FALSE;

   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, &ncontvars) );

   /* get name of the original problem and add the string "_undercoverppc" */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_undercoverppc", SCIPgetProbName(scip));

   /* create the ppc problem */
   SCIP_CALL( SCIPcreateProb(ppcscip, name, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* create a binary variable in the ppc problem for each variable in the original problem (fix it or not?) */
   SCIPallocBufferArray(scip, &conscounter, nvars);
   SCIPallocBufferArray(scip, &termcounter, nvars);
   for( i = 0; i < nvars; ++i )
   {
      /* get name of the original variable and add the string "_ppc" */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_ppc", SCIPvarGetName(vars[i]));

      SCIP_CALL( SCIPcreateVar(ppcscip, &ppcvars[i], name, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY,
            TRUE, FALSE, NULL, NULL, NULL, NULL) );

      SCIP_CALL( SCIPaddVar(ppcscip, ppcvars[i]) );

      conscounter[i] = 0.0;
      termcounter[i] = 0.0;
   }

   /* go through all quadratic and bilinear terms in the original problem */
   conshdlr = SCIPfindConshdlr(scip, "quadratic");
   if( conshdlr != NULL )
   {
      for( i = 0; i < SCIPconshdlrGetNConss(conshdlr); ++i )
      {
         SCIP_CONS* quadcons;

         SCIP_Real nterms;
         SCIP_Bool infeas;
         SCIP_Bool fixed;
         SCIP_Bool nottofix;
         int t;

         quadcons = SCIPconshdlrGetConss(conshdlr)[i];
         assert( quadcons != NULL );

         /* TODO: remove this dirty hack for counting constraint containment */
         nterms = SCIPgetNQuadVarsQuadratic(scip, quadcons) + SCIPgetNBilinTermsQuadratic(scip, quadcons);
         nterms *= 2.0;

         /* fix variable in quadratic terms to linearize/convexify it */
         for( t = 0; t < SCIPgetNQuadVarsQuadratic(scip, quadcons); ++t )
         {
            SCIP_VAR* quadvar;
            int probindex;

            quadvar = SCIPgetQuadVarsQuadratic(scip, quadcons)[t];
            assert( quadvar != NULL );

            /* if constraints with inactive variables are present, we will have difficulty creating the subscip later */
            probindex = SCIPvarGetProbindex(quadvar);
            if( probindex == -1 )
            {
               SCIPdebugMessage("undercover heuristic detected constraint <%s> with inactive variables\n", SCIPconsGetName(quadcons));
               return SCIP_OKAY;
            }

            /* if the variable is fixed in the original problem, the term is already linear: nothing to do */
            /* if the variable has zero coefficient in the original problem, the term is already linear: nothing to do */
            nottofix = local && SCIPisFeasEQ(scip, SCIPvarGetLbLocal(quadvar), SCIPvarGetUbLocal(quadvar));
            nottofix = nottofix || SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(quadvar), SCIPvarGetUbGlobal(quadvar));
            nottofix = nottofix || SCIPisZero(scip, SCIPgetSqrCoefsQuadVarsQuadratic(scip, quadcons)[t]);

            /* if we want to convexify only and coefficient and bounds are accordingly: nothing to do */
            nottofix = nottofix || (onlyconvexify && SCIPisInfinity(scip, -SCIPgetLhsQuadratic(scip, quadcons)) && SCIPgetSqrCoefsQuadVarsQuadratic(scip, quadcons)[t] >= 0);
            nottofix = nottofix || (onlyconvexify && SCIPisInfinity(scip, SCIPgetRhsQuadratic(scip, quadcons)) && SCIPgetSqrCoefsQuadVarsQuadratic(scip, quadcons)[t] <= 0);

            if( nottofix )
               continue;

            SCIP_CALL( SCIPfixVar(ppcscip, ppcvars[probindex], 1.0, &infeas, &fixed) );
            assert( !infeas );
            assert( fixed );

            conscounter[probindex] += 1.0/nterms;
            ++(termcounter[probindex]);

            SCIPdebugMessage("undercover heuristic: fixing var %s in set covering problem to 1.\n", SCIPvarGetName(ppcvars[probindex]));
         }

         /* create set covering constraints: fix at least one variable in each bilinear term */
         for( t = 0; t < SCIPgetNBilinTermsQuadratic(scip, quadcons); ++t )
         {
            SCIP_VAR* bilinvar1;
            SCIP_VAR* bilinvar2;
            SCIP_VAR* ppcconsvars[2];
            SCIP_CONS* ppccons;
            int probindex1;
            int probindex2;

            bilinvar1 = SCIPgetBilinVars1Quadratic(scip, quadcons)[t];
            bilinvar2 = SCIPgetBilinVars2Quadratic(scip, quadcons)[t];
            assert( bilinvar1 != NULL );
            assert( bilinvar2 != NULL );

            /* if constraints with inactive variables are present, we will have difficulty creating the subscip later */
            probindex1 = SCIPvarGetProbindex(bilinvar1);
            probindex2 = SCIPvarGetProbindex(bilinvar2);
            if( probindex1 == -1 || probindex2 == -1 )
            {
               SCIPdebugMessage("undercover heuristic detected constraint <%s> with inactive variables\n", SCIPconsGetName(quadcons));
               return SCIP_OKAY;
            }

            /* if one of the variables is fixed in the original problem, the term is already linear: nothing to do */
            /* if the term has zero coefficient in the original problem, it is already linear: nothing to do */
            nottofix = local && SCIPisFeasEQ(scip, SCIPvarGetLbLocal(bilinvar1), SCIPvarGetUbLocal(bilinvar1));
            nottofix = nottofix || (local && SCIPisFeasEQ(scip, SCIPvarGetLbLocal(bilinvar2), SCIPvarGetUbLocal(bilinvar2)));
            nottofix = nottofix || SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(bilinvar1), SCIPvarGetUbGlobal(bilinvar1));
            nottofix = nottofix || SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(bilinvar2), SCIPvarGetUbGlobal(bilinvar2));
            nottofix = nottofix || SCIPisZero(scip, SCIPgetBilinCoefsQuadratic(scip, quadcons)[t]);

            if( nottofix )
               continue;

            /* get name of the original constraint and add the string "_bilin.." */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_bilin%d", SCIPconsGetName(quadcons), t);

            ppcconsvars[0] = ppcvars[probindex1];
            ppcconsvars[1] = ppcvars[probindex2];

            SCIP_CALL( SCIPcreateConsSetcover(ppcscip, &ppccons, name, 2, ppcconsvars,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

            if( ppccons == NULL )
            {
               SCIPdebugMessage("failed to create set covering constraint %s: terminating undercover heuristic\n", name);
               return SCIP_OKAY;
            }

            SCIP_CALL( SCIPaddCons(ppcscip, ppccons) );
            SCIP_CALL( SCIPreleaseCons(ppcscip, &ppccons) );

            conscounter[probindex1] += 1.0/nterms;
            conscounter[probindex2] += 1.0/nterms;
            ++(termcounter[probindex1]);
            ++(termcounter[probindex2]);
         }

         for( t = 0; t < nvars; ++t )
            conscounter[t] = SCIPceil(scip, conscounter[t]);
      }
   }

   /* go through all SOC constraints in the original problem */
   conshdlr = SCIPfindConshdlr(scip, "soc");
   if( conshdlr != NULL && !onlyconvexify )
   {
      for( i = 0; i < SCIPconshdlrGetNConss(conshdlr); ++i )
      {
         SCIP_CONS* soccons;
         SCIP_CONS* ppccons;
         SCIP_VAR** soclhsvars;
         SCIP_VAR* socrhsvar;
         SCIP_VAR** ppcconsvars;

         SCIP_Bool nottofix;
         int ntofix;
         int probindex;
         int t;

         soccons = SCIPconshdlrGetConss(conshdlr)[i];
         assert( soccons != NULL );

         SCIP_CALL( SCIPallocBufferArray(ppcscip, &ppcconsvars, SCIPgetNLhsVarsSOC(scip, soccons) + 1) );
         ntofix = 0;

         /* right hand side variable */
         socrhsvar = SCIPgetRhsVarSOC(scip, soccons);
         assert( socrhsvar != NULL );

         /* if constraints with inactive variables are present, we will have difficulty creating the subscip later */
         probindex = SCIPvarGetProbindex(socrhsvar);
         if( probindex == -1 )
         {
            SCIPdebugMessage("undercover heuristic detected constraint <%s> with inactive variables\n", SCIPconsGetName(soccons));
            return SCIP_OKAY;
         }

         /* if the variable is fixed in the original problem, the term is already linear: nothing to do */
         /* if the term has zero coefficient in the original problem, it is already linear: nothing to do */
         nottofix = local && SCIPisFeasEQ(scip, SCIPvarGetLbLocal(socrhsvar), SCIPvarGetUbLocal(socrhsvar));
         nottofix = nottofix || SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(socrhsvar), SCIPvarGetUbGlobal(socrhsvar));
         nottofix = nottofix || SCIPisZero(scip, SCIPgetRhsCoefSOC(scip, soccons));

         if( !nottofix )
         {
            SCIP_CALL( SCIPgetNegatedVar(ppcscip, ppcvars[probindex], &ppcconsvars[ntofix]) );
            ++ntofix;

            conscounter[probindex] += 1.0/nvars;
            ++termcounter[probindex];
         }

         /* left hand side variables */
         soclhsvars = SCIPgetLhsVarsSOC(scip, soccons);
         assert( soclhsvars != NULL );

         for( t = 0; t < SCIPgetNLhsVarsSOC(scip, soccons); ++t )
         {
            assert( soclhsvars[t] != NULL );

            /* if constraints with inactive variables are present, we will have difficulty creating the subscip later */
            probindex = SCIPvarGetProbindex(soclhsvars[t]);
            if( probindex == -1 )
            {
               SCIPdebugMessage("undercover heuristic detected constraint <%s> with inactive variables\n", SCIPconsGetName(soccons));
               return SCIP_OKAY;
            }

            /* if the variable is fixed in the original problem, the term is already linear: nothing to do */
            /* if the term has zero coefficient in the original problem, it is already linear: nothing to do */
            nottofix = local && SCIPisFeasEQ(scip, SCIPvarGetLbLocal(soclhsvars[t]), SCIPvarGetUbLocal(soclhsvars[t]));
            nottofix = nottofix || SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(soclhsvars[t]), SCIPvarGetUbGlobal(soclhsvars[t]));
            nottofix = nottofix || (SCIPgetLhsCoefsSOC(scip, soccons) != NULL && SCIPisZero(scip, SCIPgetLhsCoefsSOC(scip, soccons)[t]));

            if( nottofix )
               continue;

            SCIP_CALL( SCIPgetNegatedVar(ppcscip, ppcvars[probindex], &ppcconsvars[ntofix]) );
            ++ntofix;

            conscounter[probindex] += 1.0/nvars;
            ++termcounter[probindex];
         }

         if( ntofix > 0 )
         {
            /* get name of the original constraint and add the string "_soc" */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_soc", SCIPconsGetName(soccons));

            SCIP_CALL( SCIPcreateConsSetpack(ppcscip, &ppccons, name, ntofix, ppcconsvars,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

            if( ppccons == NULL )
            {
               SCIPdebugMessage("failed to create set packing constraint %s: terminating undercover heuristic\n", name);
               SCIPfreeBufferArray(ppcscip, &ppcconsvars);
               return SCIP_OKAY;
            }

            SCIP_CALL( SCIPaddCons(ppcscip, ppccons) );
            SCIP_CALL( SCIPreleaseCons(ppcscip, &ppccons) );

            for( t = 0; t < nvars; ++t )
               conscounter[t] = SCIPceil(scip, conscounter[t]);
         }

         SCIPfreeBufferArray(ppcscip, &ppcconsvars);
      }
   }

#ifdef WITH_UNIVARDEFINITE
   /* go through all "univariate definite" constraints in the original problem */
   conshdlr = SCIPfindConshdlr(scip, "univardefinite");
   if( conshdlr != NULL )
   {
      for( i = 0; i < SCIPconshdlrGetNConss(conshdlr); ++i )
      {
         SCIP_CONS* uvdcons;
         SCIP_VAR* uvdvar;

         SCIP_Bool infeas;
         SCIP_Bool fixed;
         SCIP_Bool nottofix;
         int probindex;

         uvdcons = SCIPconshdlrGetConss(conshdlr)[i];
         assert( uvdcons != NULL );

         uvdvar = SCIPgetNonlinearVarUnivardefinite(scip, uvdcons);
         assert( uvdvar != NULL );

         /* if constraints with inactive variables are present, we will have difficulty creating the subscip later */
         probindex = SCIPvarGetProbindex(uvdvar);
         if( probindex == -1 || SCIPvarGetProbindex(SCIPgetLinearVarUnivardefinite(scip, uvdcons)) )
         {
            SCIPdebugMessage("undercover heuristic detected constraint <%s> with inactive variables\n", SCIPconsGetName(uvdcons));
            return SCIP_OKAY;
         }

         /* if the variable is fixed in the original problem, the term is already linear: nothing to do */
         /* if the variable has zero coefficient in the original problem, the term is already linear: nothing to do */
         nottofix = local && SCIPisFeasEQ(scip, SCIPvarGetLbLocal(uvdvar), SCIPvarGetUbLocal(uvdvar));
         nottofix = nottofix || SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(uvdvar), SCIPvarGetUbGlobal(uvdvar));
         nottofix = nottofix || SCIPisZero(scip, SCIPgetCoefNonlinearUnivardefinite(scip, uvdcons));

         /* if we want to convexify only and function and bounds are accordingly: nothing to do */
         nottofix = nottofix || (onlyconvexify && SCIPisNonlinearFunctionConvexUnivardefinite(scip, uvdcons) && SCIPisInfinity(scip, -SCIPgetLhsUnivardefinite(scip, uvdcons)));
         nottofix = nottofix || (onlyconvexify && !SCIPisNonlinearFunctionConvexUnivardefinite(scip, uvdcons) && SCIPisInfinity(scip, SCIPgetRhsUnivardefinite(scip, uvdcons)));

         if( nottofix )
            continue;

         SCIP_CALL( SCIPfixVar(ppcscip, ppcvars[probindex], 1.0, &infeas, &fixed) );
         assert( !infeas );
         assert( fixed );

         SCIPdebugMessage("undercover heuristic: fixing var %s in set covering problem to 1.\n", SCIPvarGetName(ppcvars[probindex]));
      }
   }
#endif

   /* set obj value of ppc variables according to strat */
   conobjfac = objquot > 1.0 ? objquot : 1.0;
   intobjfac = objquot > 1.0 ? 1.0 : 1.0/objquot;

   switch( strat )
   {
   case 'c':
      /* number of influenced nonlinear constraints */
      for( i = 0; i < nvars - ncontvars; ++i )
      {
         SCIP_CALL( SCIPchgVarObj(ppcscip, ppcvars[i], intobjfac*conscounter[i]) );
      }
      for( ; i < nvars; ++i )
      {
         SCIP_CALL( SCIPchgVarObj(ppcscip, ppcvars[i], conobjfac*conscounter[i]) );
      }
      break;
   case 'd':
      /* domain size */
      for( i = 0; i < nvars - ncontvars; ++i )
      {
         SCIP_CALL( SCIPchgVarObj(ppcscip, ppcvars[i],
               intobjfac*(local ? SCIPvarGetUbLocal(vars[i]) - SCIPvarGetLbLocal(vars[i]) : SCIPvarGetUbGlobal(vars[i]) - SCIPvarGetLbGlobal(vars[i]))) );
      }
      for( ; i < nvars; ++i )
      {
         SCIP_CALL( SCIPchgVarObj(ppcscip, ppcvars[i],
               conobjfac*(local ? SCIPvarGetUbLocal(vars[i]) - SCIPvarGetLbLocal(vars[i]) : SCIPvarGetUbGlobal(vars[i]) - SCIPvarGetLbGlobal(vars[i]))) );
      }
      break;
   case 'l':
      /* number of locks */
      for( i = 0; i < nvars - ncontvars; ++i )
      {
         SCIP_CALL( SCIPchgVarObj(ppcscip, ppcvars[i], intobjfac*(SCIPvarGetNLocksDown(vars[i]) + SCIPvarGetNLocksUp(vars[i]) + 1.0)) );
      }
      for( ; i < nvars; ++i )
      {
         SCIP_CALL( SCIPchgVarObj(ppcscip, ppcvars[i], conobjfac*(SCIPvarGetNLocksDown(vars[i]) + SCIPvarGetNLocksUp(vars[i]) + 1.0)) );
      }
      break;
   case 'm':
      /* min(up locks, down locks) */
      for( i = 0; i < nvars - ncontvars; ++i )
      {
         SCIP_CALL( SCIPchgVarObj(ppcscip, ppcvars[i], intobjfac*(MIN(SCIPvarGetNLocksDown(vars[i]), SCIPvarGetNLocksUp(vars[i])) + 1.0)) );
      }
      for( ; i < nvars; ++i )
      {
         SCIP_CALL( SCIPchgVarObj(ppcscip, ppcvars[i], conobjfac*(MIN(SCIPvarGetNLocksDown(vars[i]), SCIPvarGetNLocksUp(vars[i])) + 1.0)) );
      }
      break;
   case 't':
      /* number of influenced nonlinear terms */
      for( i = 0; i < nvars - ncontvars; ++i )
      {
         SCIP_CALL( SCIPchgVarObj(ppcscip, ppcvars[i], intobjfac*termcounter[i]) );
      }
      for( ; i < nvars; ++i )
      {
         SCIP_CALL( SCIPchgVarObj(ppcscip, ppcvars[i], conobjfac*termcounter[i]) );
      }
      break;
   case 'u':
      /* unit penalties for variables of same type */
      for( i = 0; i < nvars - ncontvars; ++i )
      {
         SCIP_CALL( SCIPchgVarObj(ppcscip, ppcvars[i], intobjfac) );
      }
      for( ; i < nvars; ++i )
      {
         SCIP_CALL( SCIPchgVarObj(ppcscip, ppcvars[i], conobjfac) );
      }
      break;
   default:
      SCIPerrorMessage("undercover heuristic: invalid ppc strategy <%c>\n", strat);
      SCIPfreeBufferArray(scip, &termcounter);
      SCIPfreeBufferArray(scip, &conscounter);
      return SCIP_ERROR;
   }

   SCIPfreeBufferArray(scip, &termcounter);
   SCIPfreeBufferArray(scip, &conscounter);
   *success = TRUE;
   return SCIP_OKAY;
}


/** creates a subproblem for subscip by fixing the variables with ppcsolvals[.] == 1 */
static
SCIP_RETCODE createSubProblem(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem                         */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                                */
   SCIP_Real*            ppcsolvals,         /**< ppcsolvals[i] == 1 if var. i should be fixed/dom.red. in subproblem */
   SCIP_Real             domred,             /**< reduce domain of selected variables by this factor around LP value */
   SCIP_Bool             locksrounding,      /**< shall LP values for integer vars be rounded according to locks? */
   SCIP_Bool             local,              /**< shall local LP rows be copied and local bounds be used? */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully  */
   )
{
   SCIP_VAR** vars;                          /* original scip variables                    */
   SCIP_ROW** rows;                          /* original scip rows                         */
   SCIP_HASHMAP* varmap;

   char name[SCIP_MAXSTRLEN];

   int nrows;
   int nvars;

   int fixingcounter;
   int roundedfixingcounter;
   int i;

   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* get name of the original problem and add the string "_undercoversub" */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_undercoversub", SCIPgetProbName(scip));

   /* create the subproblem */
   SCIP_CALL( SCIPcreateProb(subscip, name, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(subscip), nvars) );

   /* create the variables of the subproblem */
   fixingcounter = 0;
   roundedfixingcounter = 0;
   for( i = 0; i < nvars; ++i )
   {
      assert(SCIPvarGetProbindex(vars[i]) == i);

      /* iff ppcsolvals[i] == 1, variable is fixed/dom.red. in the subproblem, otherwise it is just copied */
      if( SCIPisEQ(scip, ppcsolvals[i], 1.0) )
      {
         SCIP_Real lpsolval;

         /* get the current LP solution value */
         lpsolval = SCIPvarGetLPSol(vars[i]);

         /* fix integer variables to rounded lpsolval */
         if( SCIPvarGetType(vars[i]) != SCIP_VARTYPE_CONTINUOUS )
         {
            if( SCIPisFeasFracIntegral(scip, lpsolval) )
               ++roundedfixingcounter;

            if( locksrounding && SCIPvarGetNLocksDown(vars[i]) != SCIPvarGetNLocksUp(vars[i]) )
            {
               lpsolval = SCIPvarGetNLocksDown(vars[i]) < SCIPvarGetNLocksUp(vars[i]) ? SCIPfeasFloor(scip, lpsolval) : SCIPfeasCeil(scip, lpsolval);
            }
            else
            {
               lpsolval = SCIPfeasFrac(scip, lpsolval) <= 0.5 ? SCIPfeasFloor(scip, lpsolval) : SCIPfeasCeil(scip, lpsolval);
            }
            assert( lpsolval >= SCIPvarGetLbGlobal(vars[i]) );
            assert( lpsolval <= SCIPvarGetUbGlobal(vars[i]) );
         }

         if( SCIPisEQ(scip, domred, 1.0) )
         {
            /* create fixed variable */
            SCIP_CALL( SCIPcreateVar(subscip, &subvars[i], SCIPvarGetName(vars[i]),
                  lpsolval, lpsolval, SCIPvarGetObj(vars[i]), SCIPvarGetType(vars[i]),
                  SCIPvarIsInitial(vars[i]), SCIPvarIsRemovable(vars[i]), NULL, NULL, NULL, NULL) );
         }
         else
         {
            SCIP_Real domsize;

            /* get the current domain size */
            domsize = local ? SCIPvarGetUbLocal(vars[i]) - SCIPvarGetLbLocal(vars[i]) : SCIPvarGetUbGlobal(vars[i]) - SCIPvarGetLbGlobal(vars[i]);
            domsize = MAX(0.0, domsize);
            if( SCIPisInfinity(scip, domsize) )
               domsize = 1.0;
            /* compute half of the new domain size */
            domsize *= 1.0 - domred;
            domsize *= 0.5;

            /* create variable with reduced domain */
            SCIP_CALL( SCIPcreateVar(subscip, &subvars[i], SCIPvarGetName(vars[i]),
                  MAX(lpsolval - domsize, local ? SCIPvarGetLbLocal(vars[i]) : SCIPvarGetLbGlobal(vars[i])),
                  MIN(lpsolval + domsize, local ? SCIPvarGetUbLocal(vars[i]) : SCIPvarGetUbGlobal(vars[i])),
                  SCIPvarGetObj(vars[i]), SCIPvarGetType(vars[i]), SCIPvarIsInitial(vars[i]), SCIPvarIsRemovable(vars[i]), NULL, NULL, NULL, NULL) );
         }
         ++fixingcounter;
      }
      else 
      {
         SCIP_CALL( SCIPcreateVar(subscip, &subvars[i], SCIPvarGetName(vars[i]),
               local ? SCIPvarGetLbLocal(vars[i]) : SCIPvarGetLbGlobal(vars[i]), local ? SCIPvarGetUbLocal(vars[i]) : SCIPvarGetUbGlobal(vars[i]),
               SCIPvarGetObj(vars[i]), SCIPvarGetType(vars[i]), SCIPvarIsInitial(vars[i]), SCIPvarIsRemovable(vars[i]), NULL, NULL, NULL, NULL) );
      }

      SCIP_CALL( SCIPaddVar(subscip, subvars[i]) );
      SCIP_CALL( SCIPhashmapInsert(varmap, vars[i], subvars[i]) );
   }

   SCIPdebugMessage("undercover heuristic fixed %d variables (%d integer variables to rounded LP value)\n", fixingcounter, roundedfixingcounter);

   /* abort, if all variables were fixed to their current LP value */
   if( fixingcounter == nvars && roundedfixingcounter == 0)
   {
      *success = FALSE;
      SCIPhashmapFree(&varmap);
      return SCIP_OKAY;
   }

   /* get the rows and their number */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) ); 

   SCIPdebugMessage("undercover heuristic copying %s lp rows\n", local ? "local" : "global");

   /* copy all rows to linear constraints */
   for( i = 0; i < nrows; i++ )
   {
      SCIP_CONS* cons;
      SCIP_VAR** consvars;
      SCIP_COL** cols;
      SCIP_Real constant;
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real* vals;
      int nnonz;
      int j;
          
      /* ignore rows that are only locally valid */
      if( !local && SCIProwIsLocal(rows[i]) )
         continue;

      /* get the row's data */
      constant = SCIProwGetConstant(rows[i]);
      lhs = SCIProwGetLhs(rows[i]) - constant;
      rhs = SCIProwGetRhs(rows[i]) - constant;
      vals = SCIProwGetVals(rows[i]);
      nnonz = SCIProwGetNNonz(rows[i]);
      cols = SCIProwGetCols(rows[i]);
      
      assert( lhs <= rhs );
      
      /* allocate memory array to be filled with the corresponding subproblem variables */
      SCIP_CALL( SCIPallocBufferArray(subscip, &consvars, nnonz) );
      for( j = 0; j < nnonz; j++ ) 
         consvars[j] = subvars[SCIPvarGetProbindex(SCIPcolGetVar(cols[j]))];

      /* create a new linear constraint and add it to the subproblem */
      SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, SCIProwGetName(rows[i]), nnonz, consvars, vals, lhs, rhs,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
      
      /* free temporary memory */
      SCIPfreeBufferArray(subscip, &consvars);
   }
   
   /* copy quadratic and soc constraints */
   for( i = 0; i < 2; ++i )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONS* cons;
      SCIP_CONS* conscopy;
      SCIP_Bool succeed;
      int c;

      conshdlr = (i == 0) ? SCIPfindConshdlr(scip, "quadratic") : SCIPfindConshdlr(scip, "soc");
      if( conshdlr == NULL )
         continue;

      SCIPdebugMessage("undercover heuristic attempting to copy %d %s constraints\n", SCIPconshdlrGetNConss(conshdlr), SCIPconshdlrGetName(conshdlr));

      for( c = 0; c < SCIPconshdlrGetNConss(conshdlr); ++c )
      {
         cons = SCIPconshdlrGetConss(conshdlr)[c];
         assert( cons != NULL );

         SCIP_CALL( SCIPcopyCons(subscip, &conscopy, NULL, conshdlr, scip, cons, varmap,
             SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
             SCIPconsIsPropagated(cons), TRUE, SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
             FALSE, &succeed) );

         if( succeed )
         {
            SCIP_CALL( SCIPaddCons(subscip, conscopy) );
            SCIP_CALL( SCIPreleaseCons(subscip, &conscopy) );
         }
         else
         {
            SCIPdebugMessage("failed to copy constraint %s\n", SCIPconsGetName(cons));
         }
      }
   }

#ifdef WITH_UNIVARDEFINITE
   /* copy "univariate definite" constraints */
   if( SCIPfindConshdlr(scip, "univardefinite") != NULL )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONS* cons;
      SCIP_CONS* conscopy;
      SCIP_Bool succeed;
      int c;

      conshdlr = SCIPfindConshdlr(scip, "univardefinite");
      assert( conshdlr != NULL );

      SCIPdebugMessage("undercover heuristic attempting to copy %d %s constraints\n", SCIPconshdlrGetNConss(conshdlr), SCIPconshdlrGetName(conshdlr));

      for( c = 0; c < SCIPconshdlrGetNConss(conshdlr); ++c )
      {
         cons = SCIPconshdlrGetConss(conshdlr)[c];
         assert( cons != NULL );

         SCIP_CALL( SCIPcopyCons(subscip, &conscopy, NULL, conshdlr, scip, cons, varmap,
             SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
             SCIPconsIsPropagated(cons), TRUE, SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
             FALSE, &succeed) );

         if( succeed )
         {
            SCIP_CALL( SCIPaddCons(subscip, conscopy) );
            SCIP_CALL( SCIPreleaseCons(subscip, &conscopy) );
         }
         else
         {
            SCIPdebugMessage("failed to copy constraint %s\n", SCIPconsGetName(cons));
         }
      }
   }
#endif

   SCIPhashmapFree(&varmap);

   *success = TRUE;
   return SCIP_OKAY;
}


/** solve ppc problem */
static
SCIP_RETCODE solvePpcProblem(
   SCIP*                 ppcscip,            /**< SCIP data structure for the ppc problem */
   SCIP_VAR**            ppcvars,            /**< variables of the ppc problem */
   SCIP_Real*            ppcsolvals,         /**< best feasible solution values */
   SCIP_Real             timelimit,          /**< time limit */
   SCIP_Real             memorylimit,        /**< memory limit */
   SCIP_Bool*            success             /**< feasible solution found? */
   )
{
#ifdef NDEBUG
   SCIP_RETCODE retstat;
#endif

#ifdef SCIP_DEBUG
   SCIP_Real totalpenalty;
   int nfixed;
   int i;
#endif

   assert( ppcscip != NULL );
   assert( ppcvars != NULL );
   assert( ppcsolvals != NULL );
   assert( timelimit > 0.0 );
   assert( memorylimit > 0.0 );
   assert( success != NULL );

   *success = FALSE;

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(ppcscip, "misc/catchctrlc", FALSE) );
 
   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(ppcscip, "display/verblevel", 0) );
 
   /* set limits for the ppc problem */
   SCIP_CALL( SCIPsetRealParam(ppcscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(ppcscip, "limits/memory", memorylimit) );

   /* forbid recursive call of heuristics solving subMIPs */
#if 0
   SCIP_CALL( SCIPsetIntParam(ppcscip, "heuristics/undercover/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(ppcscip, "heuristics/rens/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(ppcscip, "heuristics/crossover/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(ppcscip, "heuristics/oneopt/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(ppcscip, "heuristics/rins/freq", -1) ); 
   SCIP_CALL( SCIPsetIntParam(ppcscip, "heuristics/localbranching/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(ppcscip, "heuristics/mutation/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(ppcscip, "heuristics/dins/freq", -1) );
#endif

#ifdef SCIP_DEBUG
   /* for debugging, enable MIP output */
   SCIP_CALL( SCIPsetIntParam(ppcscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(ppcscip, "display/freq", 100000000) );
#endif

   /* presolve ppc problem */
#ifdef NDEBUG
   /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
    * Hence in optimized mode, the return code is catched and a warning is printed, only in debug mode, SCIP will stop.
    */
   retstat = SCIPpresolve(ppcscip);
   if( retstat != SCIP_OKAY )
   { 
      SCIPwarningMessage("error while presolving ppc problem in undercover heuristic: subSCIP terminated with code <%d>\n", retstat);
      return SCIP_OKAY;
   }
#else
   SCIP_CALL( SCIPpresolve(ppcscip) );
#endif

   SCIPdebugMessage("undercover presolved ppc problem: %d vars, %d cons\n", SCIPgetNVars(ppcscip), SCIPgetNConss(ppcscip));

   /* solve ppc problem */
#ifdef NDEBUG
   retstat = SCIPsolve(ppcscip);
   if( retstat != SCIP_OKAY )
   { 
      SCIPwarningMessage("error while solving ppc problem in undercover heuristic: subSCIP terminated with code <%d>\n", retstat);
   }
#else
   SCIP_CALL( SCIPsolve(ppcscip) );
#endif

   /* check, whether a solution was found and save the best in ppcsolvals */
   if( SCIPgetNSols(ppcscip) == 0 )
      return SCIP_OKAY;

   assert( SCIPgetBestSol(ppcscip) != NULL );

   SCIP_CALL( SCIPgetSolVals(ppcscip, SCIPgetBestSol(ppcscip), SCIPgetNOrigVars(ppcscip), ppcvars, ppcsolvals) );

#ifdef SCIP_DEBUG
   nfixed = 0;
   totalpenalty = 0.0;
   for( i = 0; i < SCIPgetNOrigVars(ppcscip); ++i )
   {
      if( ppcsolvals[i] > 0.5 )
         ++nfixed;
      totalpenalty += SCIPvarGetObj(ppcvars[i]);
   }

   SCIPdebugMessage("undercover found a ppc solution: %d/%d variables fixed, normalized penalty=%f\n",
         nfixed, SCIPgetNOrigVars(ppcscip), SCIPgetSolOrigObj(ppcscip, SCIPgetBestSol(ppcscip))/totalpenalty);
#else
   SCIPdebugMessage("undercover found a ppc solution\n");
#endif

   *success = TRUE;
   return SCIP_OKAY;
}


/** solve subproblem */
static
SCIP_RETCODE solveSubProblem(
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem */
   SCIP_Bool             checklinear,        /**< should we check if subscip is linear after presolve? */
   SCIP_Real             timelimit,          /**< time limit */
   SCIP_Real             memorylimit,        /**< memory limit */
   SCIP_Longint          nstallnodes,        /**< number of stalling nodes for the subproblem */
   SCIP_Bool*            success             /**< feasible solution found? */
   )
{
#ifdef NDEBUG
   SCIP_RETCODE retstat;
#endif

   assert( subscip != NULL );
   assert( timelimit > 0.0 );
   assert( memorylimit > 0.0 );
   assert( success != NULL );

   *success = FALSE;

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );
 
   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
 
   /* set limits for the sub problem */
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", nstallnodes) );

   /* forbid recursive call of heuristics solving subMIPs */
#if 0
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/undercover/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/rens/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/crossover/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/oneopt/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/rins/freq", -1) ); 
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/localbranching/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/mutation/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/dins/freq", -1) );
#endif

#ifdef SCIP_DEBUG
   /* for debugging, enable MIP output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100) );
#endif

   /* SCIP_CALL( SCIPreadParams(subscip, "undercoversub.set") ); */

   /* presolve subproblem */
#ifdef NDEBUG
   /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
    * Hence in optimized mode, the return code is catched and a warning is printed, only in debug mode, SCIP will stop.
    */
   retstat = SCIPpresolve(subscip);
   if( retstat != SCIP_OKAY )
   { 
      SCIPwarningMessage("error while presolving subMIQCP in undercover heuristic: subSCIP terminated with code <%d>\n", retstat);
      return SCIP_OKAY;
   }
#else
   SCIP_CALL( SCIPpresolve(subscip) );
#endif

   SCIPdebugMessage("undercover presolved subMIQCP: %d vars, %d cons\n", SCIPgetNVars(subscip), SCIPgetNConss(subscip));

   /* check for nonlinear constraints */
   if( SCIPgetStatus(subscip) != SCIP_STATUS_INFEASIBLE && checklinear )
   {
      SCIP_CONSHDLR* conshdlr;

      conshdlr = SCIPfindConshdlr(subscip, "quadratic");
      *success = conshdlr == NULL || SCIPconshdlrGetNConss(conshdlr) == 0;

      conshdlr = SCIPfindConshdlr(subscip, "soc");
      *success = *success && (conshdlr == NULL || SCIPconshdlrGetNConss(conshdlr) == 0);

#ifdef WITH_UNIVARDEFINITE
      conshdlr = SCIPfindConshdlr(subscip, "univardefinite");
      *success = *success && (conshdlr == NULL || SCIPconshdlrGetNConss(conshdlr) == 0);
#endif

      assert( *success );
      if( !(*success) )
      {
         SCIPwarningMessage("undercover presolved subMIQCP contains nonlinear constraints; terminating heuristic\n");
         return SCIP_OKAY;
      }
      else
         *success = FALSE;
   }

   /* solve subproblem */
#ifdef NDEBUG
   retstat = SCIPsolve(subscip);
   if( retstat != SCIP_OKAY )
   { 
      SCIPwarningMessage("error while solving subMIQCP in undercover heuristic: subSCIP terminated with code <%d>\n", retstat);
   }
#else
   SCIP_CALL( SCIPsolve(subscip) );
#endif

   /* check, whether a solution was found */
   if( SCIPgetNSols(subscip) == 0 )
      return SCIP_OKAY;

   *success = TRUE;
   return SCIP_OKAY;
}


/** copy the solution of the subproblem to newsol */
static
SCIP_RETCODE copySol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   SCIP_SOL**            newsol,             /**< solution to the original problem                    */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */
   int        nvars;
        
   assert( scip != NULL );
   assert( subscip != NULL );
   assert( subvars != NULL );
   assert( subsol != NULL );
   assert( newsol != NULL );
   assert( *newsol != NULL );

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   assert( nvars == SCIPgetNOrigVars(subscip) );  
 
   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );
   SCIP_CALL( SCIPsetSolVals(scip, *newsol, nvars, vars, subsolvals) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}


/** main procedure of the undercover heuristic */
SCIP_RETCODE SCIPapplyUndercover(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                       */
   SCIP_RESULT*          result,             /**< result data structure                                          */
   SCIP_Real             timelimit,          /**< time limit                                                     */        
   SCIP_Real             memorylimit,        /**< memory limit                                                   */
   char                  ppcstrat,           /**< strategy for finding a ppc solution                            */
   SCIP_Real             ppcobjquot,         /**< additional penalty factor for fixing continuous variables      */
   SCIP_Real             domred,             /**< reduce domain of selected variables by this factor around LP value */
   SCIP_Bool             locksrounding,      /**< shall LP values for integer vars be rounded according to locks? */
   SCIP_Bool             onlyconvexify,      /**< should we only fix/dom.red. variables creating nonconvexity?   */
   SCIP_Real             minimprove,         /**< factor by which heuristic should at least improve the incumbent*/
   SCIP_Longint          nstallnodes,        /**< number of stalling nodes for the subproblem                    */
   SCIP_Bool             postnlp             /**< shall NLP heuristic be called when a feas. solution was found? */
   )
{
   SCIP_VAR** vars;                          /* original problem's variables */

   SCIP* ppcscip;
   SCIP_VAR** ppcvars;                       /* ppc problem's variables */
   SCIP_Real* ppcsolvals;                    /* solution to ppc problem */
  
   SCIP* subscip;                            /* SCIP data strucutre for solving subMIQCP */
   SCIP_VAR** subvars;                       /* subMIQCP's variables */

   SCIP_Real subtimelimit;                   /* time limit for subMIQCP */
   SCIP_Bool success;

   int nvars;
   int i;

   assert( scip != NULL );
   assert( heur != NULL );
   assert( result != NULL );
   assert( *result == SCIP_DIDNOTFIND );

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* initializing ppc problem */
   SCIP_CALL( SCIPcreate(&ppcscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(ppcscip) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ppcvars, nvars) ); 

   /* create ppc problem */
   SCIPdebugMessage("undercover heuristic creating ppc problem\n");
   success = FALSE;
   SCIP_CALL( createPpcProblem(scip, ppcscip, ppcvars, ppcstrat, ppcobjquot, TRUE, FALSE, &success) );

   /* solve ppc problem */
   if( success )
   {
      SCIPdebugMessage("undercover heuristic solving ppc problem\n");
      SCIP_CALL( SCIPallocBufferArray(scip, &ppcsolvals, nvars) ); 
      SCIP_CALL( solvePpcProblem(ppcscip, ppcvars, ppcsolvals, timelimit, memorylimit, &success) );
   }

   subtimelimit = timelimit - SCIPgetTotalTime(ppcscip); /* ???????????????? */

   /* free memory from ppc problem */
   SCIP_CALL( SCIPfreeTransform(ppcscip) );
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(ppcscip, &(ppcvars[i])) );
   }
   SCIPfreeBufferArray(scip, &ppcvars);
   SCIP_CALL( SCIPfree(&ppcscip) );

   if( !success )
   {
      SCIPdebugMessage("undercover heuristic terminating: problems creating/solving ppc problem\n");
      return SCIP_OKAY;
   }

   if( subtimelimit < 10.0 )
   {
      SCIPdebugMessage("undercover heuristic terminating: subtimelimit=%f\n", subtimelimit);
      return SCIP_OKAY;
   }

   /* initializing subMIQCP */
   SCIP_CALL( SCIPcreate(&subscip) );

#ifdef WITH_CONSBRANCHNL
   SCIP_CALL( SCIPincludeConshdlrBranchNonlinear(subscip) );
#endif
   
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) ); 

#ifdef WITH_UNIVARDEFINITE
   SCIP_CALL( SCIPincludeConshdlrUnivardefinite(subscip) );
#endif
   
   /* create subMIQCP */
   SCIPdebugMessage("undercover heuristic creating subMIQCP\n");
   success = FALSE;
   SCIP_CALL( createSubProblem(scip, subscip, subvars, ppcsolvals, domred, locksrounding, TRUE, &success) );

   /* if there is already a solution, add an objective cutoff */
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_Real cutoff;

      assert( !SCIPisInfinity(scip, SCIPgetUpperbound(scip)) );   

      if( SCIPisInfinity(scip, -SCIPgetLowerbound(scip)) )
         cutoff = (SCIPgetUpperbound(scip) >= 0 ? 1.0 - minimprove : 1.0 + minimprove)*SCIPgetUpperbound(scip);
      else
         cutoff = (1.0 - minimprove)*SCIPgetUpperbound(scip) + minimprove*SCIPgetLowerbound(scip);

      cutoff = MIN(SCIPgetUpperbound(scip) - SCIPsumepsilon(scip), cutoff);
      SCIP_CALL( SCIPsetObjlimit(subscip, cutoff) );
   }

   if( success )
   {
      SCIP_SOL** subsols;
      SCIP_SOL* newsol;
      int nsubsols;

      /* solve subMIQCP */
      SCIPdebugMessage("undercover heuristic solving subMIQCP\n");
      SCIP_CALL( solveSubProblem(subscip, !onlyconvexify && SCIPisEQ(scip, domred, 1.0), subtimelimit, memorylimit, nstallnodes, &success) );

      /* check, whether a solution was found;
       * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
       */
      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      assert( success || nsubsols == 0 );
      assert( !success || nsubsols > 0 );

      /* create new solution for the original problem */
      SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );

      success = FALSE;
      for( i = 0; i < nsubsols && !success; ++i )
      {
         /* try to add new solution to scip */
         SCIP_CALL( copySol(scip, subscip, subvars, subsols[i], &newsol, &success) );
         SCIP_CALL( SCIPtrySol(scip, newsol, TRUE, TRUE, TRUE, &success) );
      }

#ifdef SCIP_DEBUG
      if( success )
      {
         assert( i >= 1 );
         if( i == 1 )
         {
            SCIPdebugMessage("undercover heuristic found %d solutions in subMIQCP; best solution feasible in original problem\n", nsubsols);
         }
         else
         {
            SCIPdebugMessage("undercover heuristic found %d solutions in subMIQCP; %d%s best solution feasible in original problem\n",
                  nsubsols, i, i == 2 ? "nd" : (i == 3 ? "rd" : "th"));
         }
      }
      else
      {
         SCIPdebugMessage("undercover heuristic found no solutions in subMIQCP\n");
      }
#endif

      if( success )
      {
         *result = SCIP_FOUNDSOL;

         /* call NLP heuristic for post optimization */
         if( postnlp )
         {
            SCIP_HEUR* nlpheur;
            SCIP_Real nlptimelimit;
            SCIP_RESULT nlpresult;

            nlpheur = SCIPfindHeur(scip, "nlp");

            if( nlpheur != NULL )
            {
               nlptimelimit = subtimelimit - SCIPgetTotalTime(subscip); /* ??????????????????? */

               if( nlptimelimit >= 10.0 )
               {
                  SCIPdebugMessage("undercover heuristic calling NLP heuristic for post optimization\n");

                  SCIP_CALL( SCIPapplyNlpHeur(scip, nlpheur, &nlpresult, newsol, INT_MAX, nlptimelimit, NULL) );

                  SCIPdebugMessage("NLP heuristic called by undercover %ssuccessfully\n", nlpresult == SCIP_FOUNDSOL ? "" : "un");
               }
            }
         }
      }

      /* free solution */
      SCIP_CALL( SCIPfreeSol(scip, &newsol) );
   }
   else
   {
      SCIPdebugMessage("undercover heuristic failed creating subMIQCP; terminating\n");
   }

   /* free remaining memory */
   SCIP_CALL( SCIPfreeTransform(subscip) );
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(subscip, &(subvars[i])) );
   }
   SCIPfreeBufferArray(scip, &subvars);
   SCIP_CALL( SCIPfree(&subscip) );
   SCIPfreeBufferArray(scip, &ppcsolvals);

   return SCIP_OKAY;
}



/*
 * Callback methods of primal heuristic
 */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeUndercover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitUndercover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* initialize data */
   heurdata->nusednodes = 0;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#define heurExitUndercover NULL


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolUndercover)
{  /*lint --e{715}*/
   SCIP_CONSHDLR* conshdlr;
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* look for nonlinear constraints */
   heurdata->run = TRUE;

   conshdlr = SCIPfindConshdlr(scip, "quadratic");
   if( conshdlr != NULL && SCIPconshdlrGetNConss(conshdlr) > 0 )
      return SCIP_OKAY;

   conshdlr = SCIPfindConshdlr(scip, "soc");
   if( conshdlr != NULL && SCIPconshdlrGetNConss(conshdlr) > 0 )
      return SCIP_OKAY;

#ifdef WITH_UNIVARDEFINITE
   conshdlr = SCIPfindConshdlr(scip, "univardefinite");
   if( conshdlr != NULL && SCIPconshdlrGetNConss(conshdlr) > 0 )
      return SCIP_OKAY;
#endif

   SCIPdebugMessage("undercover heuristic will not run for <%s> (no known nonlinear constraints present)\n", SCIPgetProbName(scip));
   heurdata->run = FALSE;
   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolUndercover NULL


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecUndercover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;                  /* heuristic's data                    */
   SCIP_Real timelimit;                      /* timelimit for the subproblem        */
   SCIP_Real memorylimit;
   SCIP_Longint nstallnodes;                 /* number of stalling nodes for the subproblem */

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );
   /* assert( SCIPhasCurrentNodeLP(scip) ); this leads to an assert e.g. in bell3a ?????????????????????? */

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* do not call heuristic, if no nonlinear constraints are present */
   if( !(heurdata->run) )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
 
   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIPdebugMessage("undercover heuristic delayed (no optimal LP solution)\n");
      *result = SCIP_DELAYED;
      return SCIP_OKAY;
   }

   /* only call heuristic, if LP solution is integral */
   if( heurdata->onlyatinteger )
   {
      SCIP_VAR** vars;
      int nvars;
      int ncontvars;
      int i;

      SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, &ncontvars) );

      for( i = 0; i < nvars - ncontvars; ++i )
      {
         if( !SCIPisFeasIntegral(scip, SCIPvarGetLPSol(vars[i])) )
         {
            SCIPdebugMessage("undercover heuristic not running (LP solution fractional)\n");
            *result = SCIP_DIDNOTRUN;
            return SCIP_OKAY;
         }
      }
   }

   *result = SCIP_DIDNOTRUN;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));
   
   /* reward heuristic if it succeeded often */
   nstallnodes = (SCIP_Longint)(nstallnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur) + 1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nstallnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-MIP as 100 nodes */
   nstallnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nstallnodes -= heurdata->nusednodes;
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);
   nstallnodes = MAX(nstallnodes, 1);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nstallnodes < heurdata->minnodes )
   {
      SCIPdebugMessage("skipping undercover heuristic: nstallnodes=%"SCIP_LONGINT_FORMAT", minnodes=%"SCIP_LONGINT_FORMAT"\n", nstallnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(scip, memorylimit) )   
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
   if( timelimit < 10.0 || memorylimit <= 0.0 )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIPdebugMessage("calling undercover heuristic for <%s>\n", SCIPgetProbName(scip));

   SCIP_CALL( SCIPapplyUndercover(scip, heur, result, timelimit, memorylimit,
         heurdata->ppcstrat, heurdata->ppcobjquot, heurdata->domred, heurdata->locksrounding, heurdata->onlyconvexify,
         heurdata->minimprove, nstallnodes, heurdata->postnlp) );

   return SCIP_OKAY;
}





/*
 * primal heuristic specific interface methods
 */

/** creates the undercover primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurUndercover(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create undercover primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING,
         heurFreeUndercover, heurInitUndercover, heurExitUndercover, 
         heurInitsolUndercover, heurExitsolUndercover, heurExecUndercover,
         heurdata) );

   /* add undercover primal heuristic parameters */
 
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
 
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/maxnodes",
         "maximum number of nodes to regard in the subproblem",

         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
 
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minimprove",
         "factor by which the heuristic should at least improve the incumbent  ",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/locksrounding",
         "shall LP values for integer vars be rounded according to locks?",
         &heurdata->locksrounding, TRUE, DEFAULT_LOCKSROUNDING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/onlyconvexify",
         "should we only fix/dom.red. variables creating nonconvexity?",
         &heurdata->onlyconvexify, TRUE, DEFAULT_ONLYCONVEXIFY, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/onlyatinteger",
         "should heuristic be called only at integer feasible nodes?",
         &heurdata->onlyatinteger, TRUE, DEFAULT_ONLYATINTEGER, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/postnlp",
         "should the NLP heuristic be called to polish a feasible solution?",
         &heurdata->postnlp, TRUE, DEFAULT_POSTNLP, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/"HEUR_NAME"/ppcstrategy",
         "strategy for the ppc problem ('b'ranching status, influenced nonlinear 'c'onstraints/'t'erms, 'd'omain size, 'l'ocks, 'm'in of up/down locks, 'u'nit penalties, constraint 'v'iolation)",
         &heurdata->ppcstrat, TRUE, DEFAULT_PPCSTRAT, PPCSTRATS, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/ppcobjquot",
         "additional penalty factor for fixing continuous variables",
         &heurdata->ppcobjquot, TRUE, DEFAULT_PPCOBJQUOT, SCIPfeastol(scip), 1.0/SCIPfeastol(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/domainreduction",
         "reduce domain of selected variables by this factor around (possibly rounded) LP value (1.0 to fix)",
         &heurdata->domred, TRUE, DEFAULT_DOMRED, 0.0, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}

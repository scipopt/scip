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

/**@file   heur_undercover.c
 * @ingroup PRIMALHEURISTICS
 * @brief  undercover primal heuristic for MIQCPs
 *
 * @author Timo Berthold
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
#include "scip/cons_and.h"
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
#ifdef WITH_SIGNPOWER
#include "cons_signpower.h"
#endif
#include "scip/heur_nlp.h"

#define HEUR_NAME             "undercover"
#define HEUR_DESC             "solves a linearization of an MIQCP determined by a set covering approach"
#define HEUR_DISPCHAR         'u'
#define HEUR_PRIORITY         -1110000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE

#define DEFAULT_MINNODES      (SCIP_Longint)500/**< minimum number of nodes to regard in the subproblem                 */
#define DEFAULT_MAXNODES      (SCIP_Longint)500/**< maximum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINIMPROVE    0.01           /**< factor by which heuristic should at least improve the incumbent       */
#define DEFAULT_NODESOFS      (SCIP_Longint)500/**< number of nodes added to the contingent of the total nodes          */
#define DEFAULT_NODESQUOT     0.1            /**< subproblem nodes in relation to nodes of the original problem         */
#define DEFAULT_PPCOBJQUOT    1.0            /**< additional penalty factor for fixing continuous variables             */
#define DEFAULT_DOMRED        1.0            /**< reduce domain of selected variables by this factor around LP value    */
#define DEFAULT_LOCKSROUNDING TRUE           /**< shall LP values for integer vars be rounded according to locks?       */
#define DEFAULT_ONLYCONVEXIFY FALSE          /**< should we only fix/dom.red. variables creating nonconvexity?          */
#define DEFAULT_ONLYATINTEGER FALSE          /**< should heuristic be called only at integer feasible nodes?            */
#define DEFAULT_GLOBALBOUNDS  FALSE          /**< should global bounds on variables be used instead of local bounds at focus node? */
#define DEFAULT_POSTNLP       TRUE           /**< should the NLP heuristic be called to polish a feasible solution?     */
#define DEFAULT_BEFORECUTS    TRUE           /**< should undercover called at root node before cut separation?          */
#define DEFAULT_FIXANDPROP    TRUE           /**< should undercover fix consecutively and propagate fixings?            */
#define DEFAULT_BACKTRACK     TRUE           /**< use one level of backtracking if infeasibility is encountered?        */
#define DEFAULT_PPCSTRAT      'u'            /**< strategy for finding a ppc solution                                   */
#define PPCSTRATS             "bcdlmtuv"     /**< strategies for finding a ppc solution                                 */


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
   SCIP_Bool             globalbounds;       /**< should global bounds on variables be used instead of local bounds at focus node? */
   SCIP_Bool             postnlp;            /**< should the NLP heuristic be called to polish a feasible solution?   */
   SCIP_Bool             beforecuts;         /**< should undercover be called at root node before cut separation?     */
   SCIP_Bool             fixandprop;         /**< should undercover fix consecutively and propagate fixings?          */
   SCIP_Bool             backtrack;          /**< use one level of backtracking if infeasibility is encountered?      */
   SCIP_Bool             run;                /**< should heuristic run, i.e. are nonlinear constraints present?       */
   char                  ppcstrat;           /**< strategy for finding a ppc solution                                 */
};


/*
 * Local methods
 */


static
/* determines, whether a term is already linear, because the variable is fixed or the coefficient is zero */
SCIP_Bool termIsLinear(
   SCIP*                 scip,               /**< SCIP data structure                 */
   SCIP_VAR*             var,                /**< variable to check                    */
   SCIP_Real             coeff,              /**< coefficient to check                 */
   SCIP_Bool             local               /**< should only local bounds be checked? */
   )
{
   SCIP_Bool islinear;

   /* if the variable has zero coefficient in the original problem, the term is linear */
   islinear = SCIPisZero(scip, coeff);

   /* if the variable is fixed in the original problem, the term is linear */
   if( local )
      islinear = islinear || SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
   else
      islinear = islinear || SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));

   return islinear;
}

static
/** determines, whether a term is convex */
SCIP_Bool termIsConvex(
   SCIP*                 scip,               /**< original SCIP data structure      */
   SCIP_Real             lhs,                /**< left hand side of the constraint  */
   SCIP_Real             rhs,                /**< right hand side of the constraint */
   SCIP_Bool             sign                /**< signature of the term             */
   )
{
   if( sign )
      return SCIPisInfinity(scip, -lhs);
   else
      return SCIPisInfinity(scip, rhs);
}

static
/** increases count by one and marks that this field has already been increased */
void  incConsCounter(
   int*                  counter,            /**< array to count in how many constraints a variable appears */
   SCIP_Bool*            marker,             /**< was this variable already counted for this constraint?    */
   int                   idx                 /**< problem index of the variable                             */
   )
{
   ++counter[idx];
   marker[idx] = TRUE;
   return;
}

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

   SCIP_Bool* consmarker;
   int* conscounter;
   int* termcounter;
   SCIP_Real conobjfac;
   SCIP_Real intobjfac;

   int nvars;
   int ncontvars;
   int i;
   int t;
   int probindex;

   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(ppcscip != NULL);
   assert(ppcvars != NULL);
   assert(strchr(PPCSTRATS, strat) != NULL);
   assert(SCIPisFeasPositive(scip, objquot) && SCIPisFeasPositive(scip, 1.0/objquot));
   assert(success != NULL);

   *success = FALSE;

   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, &ncontvars) );

   /* get name of the original problem and add the string "_undercoverppc" */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_undercoverppc", SCIPgetProbName(scip));

   /* create the ppc problem */
   SCIP_CALL( SCIPcreateProb(ppcscip, name, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* allocate and initialize to zero arrays for weighted objectives */
   SCIPallocBufferArray(scip, &consmarker, nvars);
   SCIPallocBufferArray(scip, &conscounter, nvars);
   SCIPallocBufferArray(scip, &termcounter, nvars);
   BMSclearMemoryArray(conscounter, nvars);
   BMSclearMemoryArray(termcounter, nvars);

   /* create a binary variable in the ppc problem for each variable in the original problem (fix it or not?) */
   for( i = 0; i < nvars; ++i )
   {
      /* get name of the original variable and add the string "_ppc" */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_ppc", SCIPvarGetName(vars[i]));

      SCIP_CALL( SCIPcreateVar(ppcscip, &ppcvars[i], name, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY,
            TRUE, FALSE, NULL, NULL, NULL, NULL) );

      SCIP_CALL( SCIPaddVar(ppcscip, ppcvars[i]) );
   }

   /* go through all and constraints in the original problem */
   conshdlr = SCIPfindConshdlr(scip, "and");
   if( conshdlr != NULL )
   {
      for( i = 0; i < SCIPconshdlrGetNConss(conshdlr); ++i )
      {
         SCIP_CONS* andcons;
         SCIP_CONS* ppccons;
         SCIP_VAR** andvars;
         SCIP_VAR** ppcconsvars;

         int ntofix;

         andcons = SCIPconshdlrGetConss(conshdlr)[i];
         assert(andcons != NULL);

         SCIP_CALL( SCIPallocBufferArray(ppcscip, &ppcconsvars, SCIPgetNVarsAnd(scip, andcons)) );
         ntofix = 0;
         BMSclearMemoryArray(consmarker, nvars);

         /* if constraints with inactive variables are present, we will have difficulty creating the subscip later */
         if( SCIPvarGetProbindex(SCIPgetResultantAnd(scip, andcons)) == -1 )
         {
            SCIPdebugMessage("undercover heuristic detected constraint <%s> with inactive variables\n", SCIPconsGetName(andcons));
            SCIPfreeBufferArray(ppcscip, &ppcconsvars);
            goto TERMINATE;
         }

         /* and variables */
         andvars = SCIPgetVarsAnd(scip, andcons);
         assert(andvars != NULL);

         for( t = 0; t < SCIPgetNVarsAnd(scip, andcons); ++t )
         {
            assert(andvars[t] != NULL);

            /* if constraints with inactive variables are present, we will have difficulty creating the subscip later */
            probindex = SCIPvarGetProbindex(andvars[t]);
            if( probindex == -1 )
            {
               SCIPdebugMessage("undercover heuristic detected constraint <%s> with inactive variables\n", SCIPconsGetName(andcons));
               SCIPfreeBufferArray(ppcscip, &ppcconsvars);
               goto TERMINATE;
            }

            /* if variable is already fixed: nothing to do; if already fixed to 0 entire constraint is already linear */
            if( termIsLinear(scip, andvars[t], 1.0, local) )
            {
               if( (local && SCIPisFeasZero(scip, SCIPvarGetLbLocal(andvars[t]))) || (!local && SCIPisFeasZero(scip, SCIPvarGetLbGlobal(andvars[t]))) )
               {
                  ntofix = 0;
                  break;
               }
               else
                  continue;
            }

            SCIP_CALL( SCIPgetNegatedVar(ppcscip, ppcvars[probindex], &ppcconsvars[ntofix]) );
            ++ntofix;

            ++termcounter[probindex];
            if( !consmarker[probindex] )
               incConsCounter(conscounter, consmarker, probindex);
         }

         if( ntofix > 0 )
         {
            /* get name of the original constraint and add the string "_and" */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_and", SCIPconsGetName(andcons));

            SCIP_CALL( SCIPcreateConsSetpack(ppcscip, &ppccons, name, ntofix, ppcconsvars,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

            if( ppccons == NULL )
            {
               SCIPdebugMessage("failed to create set packing constraint %s: terminating undercover heuristic\n", name);
               SCIPfreeBufferArray(ppcscip, &ppcconsvars);
               goto TERMINATE;
            }

            SCIP_CALL( SCIPaddCons(ppcscip, ppccons) );
            SCIP_CALL( SCIPreleaseCons(ppcscip, &ppccons) );
         }

         SCIPfreeBufferArray(ppcscip, &ppcconsvars);
      }
   }

   /* go through all quadratic and bilinear terms in the original problem */
   conshdlr = SCIPfindConshdlr(scip, "quadratic");
   if( conshdlr != NULL )
   {
      for( i = 0; i < SCIPconshdlrGetNConss(conshdlr); ++i )
      {
         SCIP_CONS* quadcons;
         SCIP_Real* sqrcoefs;
         SCIP_Real* bilincoefs;

         SCIP_Bool infeas;
         SCIP_Bool fixed;
         
         /* TODO: if onlyconvexify, then use SCIPisConvexQuadratic and SCIPisConcaveQuadratic to decide here if the whole constraint does not need to be considered */

         /* get coefficient arrays of constraint */
         quadcons = SCIPconshdlrGetConss(conshdlr)[i];
         assert(quadcons != NULL);
         sqrcoefs = SCIPgetSqrCoefsQuadVarsQuadratic(scip, quadcons);
         bilincoefs = SCIPgetBilinCoefsQuadratic(scip, quadcons);
         BMSclearMemoryArray(consmarker, nvars);

         /* fix variable in quadratic terms to linearize/convexify it */
         for( t = 0; t < SCIPgetNQuadVarsQuadratic(scip, quadcons); ++t )
         {
            SCIP_VAR* quadvar;
         
            quadvar = SCIPgetQuadVarsQuadratic(scip, quadcons)[t];
            assert(quadvar != NULL);

            /* if constraints with inactive variables are present, we will have difficulty creating the subscip later */
            probindex = SCIPvarGetProbindex(quadvar);
            if( probindex == -1 )
            {
               SCIPdebugMessage("undercover heuristic detected constraint <%s> with inactive variables\n", SCIPconsGetName(quadcons));
               goto TERMINATE;
            }

            if( termIsLinear(scip, quadvar, sqrcoefs[t], local) )
               continue;     

            if( onlyconvexify && termIsConvex(scip, SCIPgetLhsQuadratic(scip, quadcons), SCIPgetRhsQuadratic(scip, quadcons), sqrcoefs[t] >= 0) )
               continue;

            SCIP_CALL( SCIPfixVar(ppcscip, ppcvars[probindex], 1.0, &infeas, &fixed) );
            assert(!infeas);
            assert(fixed);

            ++termcounter[probindex];
            if( !consmarker[probindex] )
               incConsCounter(conscounter, consmarker, probindex);

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
            assert(bilinvar1 != NULL);
            assert(bilinvar2 != NULL);

            /* if constraints with inactive variables are present, we will have difficulty creating the subscip later */
            probindex1 = SCIPvarGetProbindex(bilinvar1);
            probindex2 = SCIPvarGetProbindex(bilinvar2);
            if( probindex1 == -1 || probindex2 == -1 )
            {
               SCIPdebugMessage("undercover heuristic detected constraint <%s> with inactive variables\n", SCIPconsGetName(quadcons));
               goto TERMINATE;
            }

            /* if the term is linear, because one of the variables is fixed or the coefficient is zero, continue */
            if( termIsLinear(scip, bilinvar1, bilincoefs[t], local) || termIsLinear(scip, bilinvar2, bilincoefs[t], local) )
               continue;

            /* get name of the original constraint and add the string "_bilin.." */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_bilin%d", SCIPconsGetName(quadcons), t);

            ppcconsvars[0] = ppcvars[probindex1];
            ppcconsvars[1] = ppcvars[probindex2];

            SCIP_CALL( SCIPcreateConsSetcover(ppcscip, &ppccons, name, 2, ppcconsvars,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

            if(ppccons == NULL)
            {
               SCIPdebugMessage("failed to create set covering constraint %s: terminating undercover heuristic\n", name);
               goto TERMINATE;
            }

            SCIP_CALL( SCIPaddCons(ppcscip, ppccons) );
            SCIP_CALL( SCIPreleaseCons(ppcscip, &ppccons) );

            /* increase appearance counters for both variables */
            ++termcounter[probindex1];
            if( !consmarker[probindex1] )
               incConsCounter(conscounter, consmarker, probindex1);

            ++termcounter[probindex2];
            if( !consmarker[probindex2] )
               incConsCounter(conscounter, consmarker, probindex2);
         }        
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

         int ntofix;

         soccons = SCIPconshdlrGetConss(conshdlr)[i];
         assert(soccons != NULL);

         SCIP_CALL( SCIPallocBufferArray(ppcscip, &ppcconsvars, SCIPgetNLhsVarsSOC(scip, soccons) + 1) );
         ntofix = 0;
         BMSclearMemoryArray(consmarker, nvars);

         /* right hand side variable */
         socrhsvar = SCIPgetRhsVarSOC(scip, soccons);
         assert(socrhsvar != NULL);

         /* if constraints with inactive variables are present, we will have difficulty creating the subscip later */
         probindex = SCIPvarGetProbindex(socrhsvar);
         if( probindex == -1 )
         {
            SCIPdebugMessage("undercover heuristic detected constraint <%s> with inactive variables\n", SCIPconsGetName(soccons));
            SCIPfreeBufferArray(ppcscip, &ppcconsvars);
            goto TERMINATE;
         }

         if( !termIsLinear(scip, socrhsvar, SCIPgetRhsCoefSOC(scip, soccons), local) )
         {
            SCIP_CALL( SCIPgetNegatedVar(ppcscip, ppcvars[probindex], &ppcconsvars[ntofix]) );
            ++ntofix;

            incConsCounter(conscounter, consmarker, probindex);
            ++termcounter[probindex];
         }

         /* left hand side variables */
         soclhsvars = SCIPgetLhsVarsSOC(scip, soccons);
         assert(soclhsvars != NULL);

         for( t = 0; t < SCIPgetNLhsVarsSOC(scip, soccons); ++t )
         {
            SCIP_Real coef;
            assert(soclhsvars[t] != NULL);

            /* if constraints with inactive variables are present, we will have difficulty creating the subscip later */
            probindex = SCIPvarGetProbindex(soclhsvars[t]);
            if( probindex == -1 )
            {
               SCIPdebugMessage("undercover heuristic detected constraint <%s> with inactive variables\n", SCIPconsGetName(soccons));
               SCIPfreeBufferArray(ppcscip, &ppcconsvars);
               goto TERMINATE;
            }

            coef = SCIPgetLhsCoefsSOC(scip, soccons) == NULL ? 1.0 : SCIPgetLhsCoefsSOC(scip, soccons)[t];
            if( termIsLinear(scip, soclhsvars[t], coef, local) )
               continue;

            SCIP_CALL( SCIPgetNegatedVar(ppcscip, ppcvars[probindex], &ppcconsvars[ntofix]) );
            ++ntofix;

            ++termcounter[probindex];
            if( !consmarker[probindex] )
               incConsCounter(conscounter, consmarker, probindex);
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
               goto TERMINATE;
            }

            SCIP_CALL( SCIPaddCons(ppcscip, ppccons) );
            SCIP_CALL( SCIPreleaseCons(ppcscip, &ppccons) );
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
         SCIP_CONS* ppccons;
         SCIP_VAR* linearuvdvar;
         SCIP_VAR* nonlinearuvdvar;
         SCIP_VAR* ppcconsvars[2];

         uvdcons = SCIPconshdlrGetConss(conshdlr)[i];
         assert(uvdcons != NULL);

         linearuvdvar = SCIPgetLinearVarUnivardefinite(scip, uvdcons);
         assert(linearuvdvar != NULL);

         nonlinearuvdvar = SCIPgetNonlinearVarUnivardefinite(scip, uvdcons);
         assert(nonlinearuvdvar != NULL);

         /* if constraints with inactive variables are present, we will have difficulty creating the subscip later */
         if( SCIPvarGetProbindex(linearuvdvar) == -1 || SCIPvarGetProbindex(nonlinearuvdvar) == -1 )
         {
            SCIPdebugMessage("undercover heuristic detected constraint <%s> with inactive variables\n", SCIPconsGetName(uvdcons));
            goto TERMINATE;
         }

         /* if one of the variables is fixed or has zero coefficient, continue */
         if( termIsLinear(scip, linearuvdvar, SCIPgetCoefLinearUnivardefinite(scip, uvdcons), local)
            || termIsLinear(scip, nonlinearuvdvar, SCIPgetCoefNonlinearUnivardefinite(scip, uvdcons), local) )
            continue; 

         /* if we want to convexify only and function and bounds are accordingly: nothing to do */
         if( onlyconvexify && termIsConvex(scip, SCIPgetLhsUnivardefinite(scip, uvdcons), SCIPgetRhsUnivardefinite(scip, uvdcons), SCIPisNonlinearFunctionConvexUnivardefinite(scip, uvdcons)) )
            continue;

         if( local ? SCIPgetMonotonicityLocalUnivardefinite(scip, uvdcons) : SCIPgetMonotonicityGlobalUnivardefinite(scip, uvdcons) )
         {
            /* get name of the original constraint and add the string "_uvd" */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_uvd", SCIPconsGetName(uvdcons));

            ppcconsvars[0] = ppcvars[SCIPvarGetProbindex(linearuvdvar)];
            ppcconsvars[1] = ppcvars[SCIPvarGetProbindex(nonlinearuvdvar)];

            SCIP_CALL( SCIPcreateConsSetcover(ppcscip, &ppccons, name, 2, ppcconsvars,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

            if( ppccons == NULL )
            {
               SCIPdebugMessage("failed to create set covering constraint %s: terminating undercover heuristic\n", name);
               goto TERMINATE;
            }

            SCIP_CALL( SCIPaddCons(ppcscip, ppccons) );
            SCIP_CALL( SCIPreleaseCons(ppcscip, &ppccons) );
         }
         else
         {
            SCIP_Bool infeas;
            SCIP_Bool fixed;

            SCIP_CALL( SCIPfixVar(ppcscip, ppcvars[SCIPvarGetProbindex(nonlinearuvdvar)], 1.0, &infeas, &fixed) );
            assert(!infeas);
            assert(fixed);

            SCIPdebugMessage("undercover heuristic: fixing var %s in set covering problem to 1.\n", SCIPvarGetName(ppcvars[SCIPvarGetProbindex(nonlinearuvdvar)]));
         }

         ++(conscounter[SCIPvarGetProbindex(linearuvdvar)]);
         ++(conscounter[SCIPvarGetProbindex(nonlinearuvdvar)]);
         ++(termcounter[SCIPvarGetProbindex(linearuvdvar)]);
         ++(termcounter[SCIPvarGetProbindex(nonlinearuvdvar)]);
      }
   }
#endif

#ifdef WITH_SIGNPOWER
   /* go through all "signpower" constraints in the original problem */
   conshdlr = SCIPfindConshdlr(scip, "signpower");
   if( conshdlr != NULL )
   {
      for( i = 0; i < SCIPconshdlrGetNConss(conshdlr); ++i )
      {
         SCIP_CONS* spcons;
         SCIP_CONS* ppccons;
         SCIP_VAR* linearspvar;
         SCIP_VAR* nonlinearspvar;
         SCIP_VAR* ppcconsvars[2];

         spcons = SCIPconshdlrGetConss(conshdlr)[i];
         assert(spcons != NULL);

         linearspvar = SCIPgetLinearVarSignpower(scip, spcons);
         assert(linearspvar != NULL);

         nonlinearspvar = SCIPgetNonlinearVarSignpower(scip, spcons);
         assert(nonlinearspvar != NULL);

         /* if constraints with inactive variables are present, we will have difficulty creating the subscip later */
         if( SCIPvarGetProbindex(linearspvar) == -1 || SCIPvarGetProbindex(nonlinearspvar) == -1 )
         {
            SCIPdebugMessage("undercover heuristic detected constraint <%s> with inactive variables\n", SCIPconsGetName(spcons));
            goto TERMINATE;
         }

         /* if one of the variables is fixed or has zero coefficient, continue */
         if( termIsLinear(scip, linearspvar, SCIPgetCoefLinearSignpower(scip, spcons), local)
            || termIsLinear(scip, nonlinearspvar, 1.0, local) )
            continue;

         /* if we want to convexify only and constraint and bounds are accordingly: nothing to do */
         if( onlyconvexify )
         {
            /* if constraint is x|x|^{n-1} - cz <= rhs with x >= 0, then it's a convex constraint already */
            if( SCIPisInfinity(scip, -SCIPgetLhsSignpower(scip, spcons)) && !SCIPisNegative(scip, local ? SCIPvarGetLbLocal(nonlinearspvar) : SCIPvarGetLbGlobal(nonlinearspvar)) )
               continue;
            /* if constraint is lhs <= x|x|^{n-1} - cz with x <= 0, then it's a convex constraint already */
            if( SCIPisInfinity(scip,  SCIPgetRhsSignpower(scip, spcons)) && !SCIPisPositive(scip, local ? SCIPvarGetUbLocal(nonlinearspvar) : SCIPvarGetUbGlobal(nonlinearspvar)) )
               continue;
         }

         /* need to fix either linear or nonlinear variable to get rid of nonlinearity */

         /* get name of the original constraint and add the string "_sp" */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_sp", SCIPconsGetName(spcons));

         ppcconsvars[0] = ppcvars[SCIPvarGetProbindex(linearspvar)];
         ppcconsvars[1] = ppcvars[SCIPvarGetProbindex(nonlinearspvar)];

         SCIP_CALL( SCIPcreateConsSetcover(ppcscip, &ppccons, name, 2, ppcconsvars,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

         if( ppccons == NULL )
         {
            SCIPdebugMessage("failed to create set covering constraint %s: terminating undercover heuristic\n", name);
            goto TERMINATE;
         }

         SCIP_CALL( SCIPaddCons(ppcscip, ppccons) );
         SCIP_CALL( SCIPreleaseCons(ppcscip, &ppccons) );

         ++(conscounter[SCIPvarGetProbindex(linearspvar)]);
         ++(conscounter[SCIPvarGetProbindex(nonlinearspvar)]);
         ++(termcounter[SCIPvarGetProbindex(linearspvar)]);
         ++(termcounter[SCIPvarGetProbindex(nonlinearspvar)]);
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
      SCIPfreeBufferArray(scip, &consmarker);
      return SCIP_ERROR;
   }

   *success = TRUE;

 TERMINATE: 
   SCIPfreeBufferArray(scip, &termcounter); 
   SCIPfreeBufferArray(scip, &conscounter);
   SCIPfreeBufferArray(scip, &consmarker);
 
   return SCIP_OKAY;
}

static
/** calculate the reduced bounds of a variable in the subMIQCP */
void calculateBounds(
   SCIP*                 scip, 
   SCIP_VAR*             var, 
   SCIP_Real*            lb, 
   SCIP_Real*            ub,
   SCIP_Real             domred,
   int*                  roundedfixingcounter,
   SCIP_Bool             locksrounding
   )
{  
   SCIP_Real lpsolval;
   SCIP_Real domsize;

   /* get the current LP solution value */
   lpsolval = SCIPvarGetLPSol(var);

   /* fix integer variables to rounded lpsolval */
   if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
   {
      if( SCIPisFeasIntegral(scip, lpsolval) )
         ++roundedfixingcounter;
      else if( locksrounding )
      {
         if( SCIPvarGetNLocksDown(var) < SCIPvarGetNLocksUp(var) ) 
            lpsolval = SCIPfeasFloor(scip, lpsolval);
         else if( SCIPvarGetNLocksDown(var) > SCIPvarGetNLocksUp(var) )
            lpsolval = SCIPfeasCeil(scip, lpsolval);
         else 
            lpsolval = SCIPfeasFrac(scip, lpsolval) <= 0.5 ? SCIPfeasFloor(scip, lpsolval) : SCIPfeasCeil(scip, lpsolval);
      }
      else
      {
         if( SCIPfeasFrac(scip, lpsolval) < 0.5)
            lpsolval = SCIPfeasFloor(scip, lpsolval);
         else if( SCIPfeasFrac(scip, lpsolval) > 0.5 )
            lpsolval = SCIPfeasCeil(scip, lpsolval);
         else
            lpsolval = SCIPvarGetNLocksDown(var) < SCIPvarGetNLocksUp(var) ? SCIPfeasFloor(scip, lpsolval) : SCIPfeasCeil(scip, lpsolval);
      }

      assert(SCIPisFeasIntegral(scip, lpsolval));
      assert(SCIPvarGetLbGlobal(var) <= lpsolval && lpsolval <= SCIPvarGetUbGlobal(var));
   }

   /* get the current domain size */
   domsize = ub-lb;
   domsize = MAX(0.0, domsize);
   if( SCIPisInfinity(scip, domsize) )
      domsize = 1.0;

   /* compute half of the new domain size */
   domsize *= 1.0 - domred;
   domsize *= 0.5;

   /* move lpsolval into the domain. It may be slightly outside due to numerical issues,
      or the domain may have been infered due to propagation */
   lpsolval = MIN(lpsolval,*ub);
   lpsolval = MAX(lpsolval,*lb);

   *lb = MAX(*lb, lpsolval - domsize);
   *ub = MIN(*ub, lpsolval + domsize);

   assert(*lb <= *ub);
}


/** calculates up to four alternative values for backtracking, if fixing the variable failed. 
 * The alternatives are the two bounds of the variable, and the averages of the bounds and the fixing value.
 * For infinite bounds, fixval +/- abs(fixval) will be used instead. 
 */
static
void calculateAlternatives(
   SCIP*                 scip, 
   SCIP_VAR*             var, 
   SCIP_Real             fixval, 
   SCIP_Real*            alternatives, 
   int*                  nalternatives
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   /* for binary variables, there is only one possible alternative value for backtracking */
   if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
   {
      assert(SCIPisFeasEQ(scip, fixval, 0.0) || SCIPisFeasEQ(scip, fixval, 1.0));
      alternatives[0] = 1.0 - fixval;
      *nalternatives = 1;
      return;
   }
   
   *nalternatives = 0;
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   
   assert(!SCIPisEQ(scip, lb, ub));
 
   /* use the lower bound as alternative value, if it is infinite, use x'-|x'| */  
   if( SCIPisInfinity(scip, -lb) )
   { 
      /* lower bound is infinite, use x'-|x'|. If x' is zero, use -1.0 instead */
      if( SCIPisFeasZero(scip, fixval) )
         alternatives[*nalternatives] = -1.0;
      else
         alternatives[*nalternatives] = fixval - ABS(fixval);

      assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPisFeasIntegral(scip,alternatives[*nalternatives]));
      (*nalternatives)++;     
   }
   else if( !SCIPisEQ(scip, lb, fixval) )
   {
      /* lower bound is finite. Use it as alternative backtracking value */
      alternatives[*nalternatives] = lb;
      (*nalternatives)++;
   }
   
   /* use the upper bound as alternative value, if it is infinite, use x'+|x'| instead */  
   if( SCIPisInfinity(scip, ub) )
   {
      /* upper bound is infinite, use x'+|x'|. If x' is zero, use 1.0 instead */
      if( SCIPisZero(scip, fixval) )
         alternatives[*nalternatives] = 1.0;
      else
         alternatives[*nalternatives] = fixval + ABS(fixval);

      assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPisFeasIntegral(scip,alternatives[*nalternatives]));
      (*nalternatives)++;     
   }
   else if( !SCIPisEQ(scip, ub, fixval) )
   {
      /* upper bound is finite. Use it as alternative backtracking value */
      alternatives[*nalternatives] = ub;
      (*nalternatives)++;
   }
   assert(*nalternatives == 1 || *nalternatives == 2);

   /* use the average of x' and lower bound as alternative value, if this is not equal to any of the other values */  
   if( !SCIPisInfinity(scip, -lb) && !SCIPisEQ(scip, lb, fixval) 
      && (SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || !SCIPisEQ(scip, lb, fixval-1)) )
   {
      alternatives[*nalternatives] = (lb+fixval)/2.0;
      
      /* round up for discrete variables */
      if(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS)
      {
         alternatives[*nalternatives] = SCIPceil(scip, alternatives[*nalternatives]);
         assert(!SCIPisEQ(scip, alternatives[*nalternatives], fixval));
      }         
      (*nalternatives)++; 
   }

   /* use the average of x' and upper bound as alternative value, if this is not equal to any of the other values */  
   if( !SCIPisInfinity(scip, ub) && !SCIPisEQ(scip, ub, fixval) 
      && (SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || !SCIPisEQ(scip, ub, fixval+1)) )
   {
      alternatives[*nalternatives] = (ub+fixval)/2.0;
      
      if(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS)
      {
         alternatives[*nalternatives] = SCIPfloor(scip, alternatives[*nalternatives]);
         assert(!SCIPisEQ(scip, alternatives[*nalternatives], fixval));
      }         
      (*nalternatives)++; 
   }  
   
#ifndef NDEBUG
   {
      int i;
      int j;
      for( i=0; i < *nalternatives; ++i)
      {
         assert(!SCIPisEQ(scip, alternatives[i], fixval));
         for( j=i+1; j < *nalternatives; ++j)
            assert(!SCIPisEQ(scip, alternatives[i], alternatives[j]));               
      }
   }
#endif     

}

/** creates a subproblem for subSCIP by fixing the variables with ppcsolvals[.] == 1 */
static
SCIP_RETCODE createSubProblem(
   SCIP*                 scip,               /**< original SCIP data structure                                        */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem                              */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                                     */
   SCIP_Real*            ppcsolvals,         /**< ppcsolvals[i] == 1 if var. i should be fixed/dom.red. in subproblem */
   SCIP_Real             domred,             /**< reduce domain of selected variables by this factor around LP value  */
   SCIP_Bool             fixandprop,         /**< should undercover fix consecutively and propagate fixings?          */
   SCIP_Bool             backtrack,          /**< use one level of backtracking if infeasibility is encountered?      */
   SCIP_Bool             locksrounding,      /**< shall LP values for integer vars be rounded according to locks?     */
   SCIP_Bool             local,              /**< shall local LP rows be copied and local bounds be used?             */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully       */
   )
{
   SCIP_VAR** vars;                          /* original SCIP variables */
   SCIP_HASHMAP* varmap;

   char name[SCIP_MAXSTRLEN];
   int nvars;
   int fixingcounter;
   int roundedfixingcounter;
   int i;

   *success = FALSE;

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
      SCIP_Real lb;
      SCIP_Real ub;
      
      assert(SCIPvarGetProbindex(vars[i]) == i);

      lb = local ? SCIPvarGetLbLocal(vars[i]) :  SCIPvarGetLbGlobal(vars[i]);
      ub = local ? SCIPvarGetUbLocal(vars[i]) :  SCIPvarGetLbGlobal(vars[i]);
     
      /* iff ppcsolvals[i] == 1, variable is fixed/dom.red. in the subproblem */
      if( SCIPisEQ(scip, ppcsolvals[i], 1.0) && !fixandprop )
      {
         calculateBounds(scip, vars [i], &lb, &ub, domred, &roundedfixingcounter, locksrounding);  

         SCIPdebugMessage("%s %s variable <%s> to [%g, %g] in subMIQCP\n", SCIPisEQ(scip, lb, ub) ? "fix" : "restrict", 
            SCIPvarGetType(vars[i]) == SCIP_VARTYPE_CONTINUOUS ? "continuous" : "discrete",
            SCIPvarGetName(vars[i]), SCIPvarGetLbGlobal(subvars[i]), SCIPvarGetUbGlobal(subvars[i]));
         
         ++fixingcounter;         
      }
      
      /* variable is created. domain reductions were only performed when we want to fix all variables at a time, without intermediate propagation */
      SCIP_CALL( SCIPcreateVar(subscip, &subvars[i], SCIPvarGetName(vars[i]), lb, ub, SCIPvarGetObj(vars[i]), 
            SCIPvarGetType(vars[i]), SCIPvarIsInitial(vars[i]), SCIPvarIsRemovable(vars[i]), NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(subscip, subvars[i]) );
      SCIP_CALL( SCIPhashmapInsert(varmap, vars[i], subvars[i]) );
   }

   SCIPdebugMessage("undercover heuristic fixed %d variables (%d integer variables to rounded LP value)\n", fixingcounter, roundedfixingcounter);

   /* fix-and-propagate loop */
   if( fixandprop )
   {
      int nbacktracks;
      int nalternatives;
      SCIP_Real* alternatives;

      /* start probing */
      SCIP_CALL( SCIPstartProbing(scip) );
      
      /* if we want to use global bounds and are not at the root node, reset bounds */
      for( i = 0; i < nvars; ++i )
      {
         if( SCIPgetDepth(scip) > 0 && !local )
         {
            SCIP_CALL( SCIPchgVarLbProbing(scip, vars[i], SCIPvarGetLbGlobal(vars[i])) );
            SCIP_CALL( SCIPchgVarUbProbing(scip, vars[i], SCIPvarGetUbGlobal(vars[i])) );
         }
      } 

      /* we try at most four alternative values in backtracking */
      SCIP_CALL( SCIPallocBufferArray(scip, &alternatives, 4) );
    
      fixingcounter = 0;
      roundedfixingcounter = 0;

      /* for each variable, which has to be fixed: fix and propagate */   
      for( i = nvars-1; i >= 0; --i )
      {
         assert(SCIPvarGetProbindex(vars[i]) == i);

         /* iff ppcsolvals[i] == 1, variable is fixed/dom.red. in the subproblem */
         if( SCIPisEQ(scip, ppcsolvals[i], 1.0)  )
         {
            SCIP_Real lb;
            SCIP_Real ub;
            SCIP_Real oldlb;
            SCIP_Real oldub;

            SCIP_Bool infeasible;
            SCIP_Longint ndomredsfound;
  
            lb = SCIPvarGetLbLocal(vars[i]);
            ub = SCIPvarGetUbLocal(vars[i]);

            /* if the variable was already fixed, e.g. by propagation, continue */
            if( SCIPisEQ(scip, lb, ub) )
               continue;
            
            oldlb = lb;
            oldub = ub;
            nalternatives = 0;

            /* create next probing node */
            SCIP_CALL( SCIPnewProbingNode(scip) );

            /* calculate the new bounds of the variable */
            calculateBounds(scip, vars[i], &lb, &ub, domred, &roundedfixingcounter, locksrounding);  
            if( backtrack )
            {
               /* Backtracking without fixing is not supported yet */
               if( domred != 1.0 )
               {
                  SCIPerrorMessage("Backtracking for domain reduction not implemented yet.\n");
                  return SCIP_INVALIDDATA;
               }
               else
               {
                  /* calculate between one and four alternative fixing values for the case of backtracking */
                  assert(SCIPisEQ(scip,lb,ub));
                  calculateAlternatives(scip, vars[i], lb, alternatives, &nalternatives);
                  assert(nalternatives >= 1);
               }
            }

            nbacktracks = 0;
            /* backtracking loop */
            do
            {
               /* if we are in backtracking, try next alternative value */
               if( nbacktracks > 0 )
               {          
                  lb = alternatives[nbacktracks-1];
                  ub = alternatives[nbacktracks-1];                  
               }

               /* change variable bounds in probing, subMIQCP will follow later */
               if( SCIPisLbBetter(scip, lb, oldlb, oldub) )
               {
                  SCIP_CALL( SCIPchgVarLbProbing(scip, vars[i], lb) );
               }
               if( SCIPisUbBetter(scip, ub, oldlb, oldub) )
               {
                  SCIP_CALL( SCIPchgVarUbProbing(scip, vars[i], ub) );
               }
               
               SCIPdebugMessage("tentatively %s variable <%s> to [%g, %g] for probing\n", SCIPisEQ(scip, lb, ub) ? "fix" : "restrict",                   
                  SCIPvarGetName(vars[i]), lb, ub);

               /* propagate the bound change */       
               SCIP_CALL( SCIPpropagateProbing(scip, 0, &infeasible, &ndomredsfound) );
               SCIPdebugMessage("  --> propagation reduced %lld further domains\n", ndomredsfound);
               
               /* if propagation led to a cutoff, backtrack or abort */
               if( infeasible )
               {  
                  if( backtrack && nbacktracks < nalternatives )
                  {
                     SCIPdebugMessage("  --> cutoff detected - backtracking\n");
                     SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip)-1) );
                     SCIP_CALL( SCIPnewProbingNode(scip) );
                     nbacktracks++;
                  }
                  else 
                  {
                     SCIPfreeBufferArray(scip, &alternatives);
                     SCIP_CALL( SCIPendProbing(scip) );
                     SCIPdebugMessage("  --> cutoff detected - abort\n");
                     goto TERMINATE;
                  }
               }
               else 
               {
                  /* if the fixing did not lead to a cutoff, transfer it to the subMIQCP */
                  SCIPdebugMessage("  --> finally %s %s variable <%s> to [%g, %g] in subMIQCP\n", SCIPisEQ(scip, lb, ub) ? "fix" : "restrict", 
                  SCIPvarGetType(vars[i]) == SCIP_VARTYPE_CONTINUOUS ? "continuous" : "discrete",
                  SCIPvarGetName(vars[i]), lb, ub);
                  SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i], lb) );
                  SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i], ub) );
                  break;
               }
            }
            while( backtrack && nbacktracks <= nalternatives );

            ++fixingcounter;         
         }
      }
      SCIPfreeBufferArray(scip, &alternatives);
      SCIP_CALL( SCIPendProbing(scip) );

      SCIPdebugMessage("undercover heuristic fixed %d variables (%d integer variables to rounded LP value) during probing\n", fixingcounter, roundedfixingcounter);

   }

   /* abort, if nothing was fixed or all variables were fixed to their current LP value */
   /* CAN this happen??????? The first condition means, that we have a MIP, 
    * the second, that the LP solution was integral
    * if( fixingcounter == 0 || (fixingcounter == nvars && roundedfixingcounter == 0) )
    *   goto TERMINATE;
    */

   /* copy all constraints */
   for( i = 0; i < SCIPgetNConshdlrs(scip); ++i )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONS* cons;
      SCIP_CONS* conscopy;
      SCIP_Bool succeed;
      int c;

      conshdlr = SCIPgetConshdlrs(scip)[i];

      SCIPdebugMessage("undercover heuristic attempting to copy %d %s constraints\n", SCIPconshdlrGetNConss(conshdlr), SCIPconshdlrGetName(conshdlr));

      /* loop over all constraint handlers and copy all their constraints */
      for( c = 0; c < SCIPconshdlrGetNConss(conshdlr); ++c )
      {
         cons = SCIPconshdlrGetConss(conshdlr)[c];
         assert(cons != NULL);

         SCIP_CALL( SCIPcopyCons(subscip, &conscopy, NULL, conshdlr, scip, cons, varmap,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
               SCIPconsIsPropagated(cons), TRUE, SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
               FALSE, &succeed) );

         /* some constraint handlers may not have a copy constructor. Since we only need a relaxation of the problem, we can proceed */
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

   *success = TRUE;

 TERMINATE:
   SCIPhashmapFree(&varmap);

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

   assert(ppcscip != NULL);
   assert(ppcvars != NULL);
   assert(ppcsolvals != NULL);
   assert(timelimit > 0.0);
   assert(memorylimit > 0.0);
   assert(success != NULL);

   *success = FALSE;

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(ppcscip, "misc/catchctrlc", FALSE) );
 
   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(ppcscip, "display/verblevel", 0) );
 
   /* set limits for the ppc problem */
   SCIP_CALL( SCIPsetRealParam(ppcscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(ppcscip, "limits/memory", memorylimit) );

   /* forbid recursive call of heuristics solving subMIPs */
   SCIP_CALL( SCIPsetIntParam(ppcscip, "heuristics/undercover/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(ppcscip, "heuristics/rens/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(ppcscip, "heuristics/crossover/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(ppcscip, "heuristics/oneopt/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(ppcscip, "heuristics/rins/freq", -1) ); 
   SCIP_CALL( SCIPsetIntParam(ppcscip, "heuristics/localbranching/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(ppcscip, "heuristics/mutation/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(ppcscip, "heuristics/dins/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(ppcscip, "separating/rapidlearning/freq", -1) );

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

   assert(SCIPgetBestSol(ppcscip) != NULL);

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
   SCIP_Longint          nodelimit,          /**< node limit */
   SCIP_Longint          nstallnodes         /**< number of stalling nodes for the subproblem */
   )
{
#ifdef NDEBUG
   SCIP_RETCODE retstat;
#endif

   assert(subscip != NULL);
   assert(timelimit > 0.0);
   assert(memorylimit > 0.0);

   SCIP_CALL( SCIPsetHeuristicsAggressive(subscip) );
#ifndef WITH_UNIVARDEFINITE
   SCIP_CALL( SCIPreadParams(subscip, "settings/emphasis/feasibility.set") );
   SCIP_CALL( SCIPreadParams(subscip, "settings/presolving/fast.set") );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", 100) );
#else
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", nstallnodes) );
#endif

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );
 
   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
 
   /* set limits for the sub problem */
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nodelimit) );

   /* forbid recursive call of undercover heuristic */
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/undercover/freq", -1) );

#ifdef SCIP_DEBUG
   /* for debugging, enable MIP output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100) );
#endif

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

#ifndef NDEBUG
   /* check for nonlinear constraints */
   if( SCIPgetStatus(subscip) != SCIP_STATUS_INFEASIBLE && checklinear )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_Bool  islinear;

      islinear = TRUE;

      conshdlr = SCIPfindConshdlr(subscip, "and");
      islinear = islinear && (conshdlr == NULL || SCIPconshdlrGetNConss(conshdlr) == 0);

      conshdlr = SCIPfindConshdlr(subscip, "quadratic");
      islinear = islinear && (conshdlr == NULL || SCIPconshdlrGetNConss(conshdlr) == 0);

      conshdlr = SCIPfindConshdlr(subscip, "soc");
      islinear = islinear && (conshdlr == NULL || SCIPconshdlrGetNConss(conshdlr) == 0);

      conshdlr = SCIPfindConshdlr(subscip, "univardefinite");
      islinear = islinear && (conshdlr == NULL || SCIPconshdlrGetNConss(conshdlr) == 0);

      conshdlr = SCIPfindConshdlr(subscip, "signpower");
      islinear = islinear && (conshdlr == NULL || SCIPconshdlrGetNConss(conshdlr) == 0);

      assert(islinear);
   }
#endif

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
        
   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);
   assert(newsol != NULL);
   assert(*newsol != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   assert(nvars == SCIPgetNOrigVars(subscip));
 
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
   SCIP_Longint          nodelimit,          /**< node limit                                                     */
   char                  ppcstrat,           /**< strategy for finding a ppc solution                            */
   SCIP_Real             ppcobjquot,         /**< additional penalty factor for fixing continuous variables      */
   SCIP_Real             domred,             /**< reduce domain of selected variables by this factor around LP value */
   SCIP_Bool             fixandprop,         /**< should undercover fix consecutively and propagate fixings?          */
   SCIP_Bool             backtrack,          /**< use one level of backtracking if infeasibility is encountered?      */
   SCIP_Bool             locksrounding,      /**< shall LP values for integer vars be rounded according to locks? */
   SCIP_Bool             onlyconvexify,      /**< should we only fix/dom.red. variables creating nonconvexity?   */
   SCIP_Bool             globalbounds,       /**< should global bounds on variables be used instead of local bounds at focus node? */
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

   SCIP_Bool success;

   int nvars;
   int i;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(result != NULL);
   assert(*result == SCIP_DIDNOTFIND);

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* initializing ppc problem */
   SCIP_CALL( SCIPcreate(&ppcscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(ppcscip) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ppcvars, nvars) ); 
   SCIP_CALL( SCIPallocBufferArray(scip, &ppcsolvals, nvars) ); 
   /* create ppc problem */
   SCIPdebugMessage("undercover heuristic creating ppc problem\n");
   success = FALSE;
   SCIP_CALL( createPpcProblem(scip, ppcscip, ppcvars, ppcstrat, ppcobjquot, !globalbounds, onlyconvexify, &success) );
   
   if( !success )
   {
      SCIPdebugMessage("undercover heuristic terminating: problems creating ppc problem\n");
      goto TERMINATEPPC;
   }

   /* solve ppc problem */
   SCIPdebugMessage("undercover heuristic solving ppc problem\n");
   SCIP_CALL( solvePpcProblem(ppcscip, ppcvars, ppcsolvals, timelimit, memorylimit, &success) );

   if( !success )
   {  
      SCIPdebugMessage("undercover heuristic terminating: problems solving ppc problem\n");
      goto TERMINATEPPC;
   }

   timelimit -= SCIPgetTotalTime(ppcscip);
   if( timelimit < 10.0 )
   { 
      SCIPdebugMessage("undercover heuristic terminating: subtimelimit=%f\n", timelimit);
      goto TERMINATEPPC;
   }

   /* initializing subMIQCP */
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) ); 

#ifdef WITH_CONSBRANCHNL
   SCIP_CALL( SCIPincludeConshdlrBranchNonlinear(subscip) );
#endif
   
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

#ifdef WITH_UNIVARDEFINITE
   SCIP_CALL( SCIPincludeConshdlrUnivardefinite(subscip) );
#endif
#ifdef WITH_SIGNPOWER
   SCIP_CALL( SCIPincludeConshdlrSignpower(subscip) );
#endif

   /* create subMIQCP */
   SCIPdebugMessage("undercover heuristic creating subMIQCP\n");
   success = FALSE;
   SCIP_CALL( createSubProblem(scip, subscip, subvars, ppcsolvals, domred, fixandprop, backtrack, locksrounding, !globalbounds, &success) );

   if( !success )
   { 
      SCIPdebugMessage("undercover heuristic terminating: problems creating subproblem\n");
      goto TERMINATE;
   }

   /* if there is already a solution, add an objective cutoff */
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_Real cutoff;

      assert(!SCIPisInfinity(scip, SCIPgetUpperbound(scip)));

      if( SCIPisInfinity(scip, -SCIPgetLowerbound(scip)) )
         cutoff = (SCIPgetUpperbound(scip) >= 0 ? 1.0 - minimprove : 1.0 + minimprove)*SCIPgetUpperbound(scip);
      else
         cutoff = (1.0 - minimprove)*SCIPgetUpperbound(scip) + minimprove*SCIPgetLowerbound(scip);

      cutoff = MIN(SCIPgetUpperbound(scip) - SCIPsumepsilon(scip), cutoff);
      SCIP_CALL( SCIPsetObjlimit(subscip, cutoff) );
   }
      
   /* solve subMIQCP */
   SCIPdebugMessage("undercover heuristic solving subMIQCP\n");
   SCIP_CALL( solveSubProblem(subscip, !onlyconvexify && SCIPisEQ(scip, domred, 1.0), timelimit, memorylimit, nodelimit, nstallnodes) );
 
   success = FALSE;
   
   /* check, whether a solution was found */
   if( SCIPgetNSols(subscip) > 0 )
   {
      SCIP_SOL** subsols;
      SCIP_SOL* newsol;
      int nsubsols;

      /* check, whether a solution was found;
       * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
       */
      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      assert(subsols != NULL);

      /* create new solution for the original problem */
      SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );

      for( i = 0; i < nsubsols && !success; ++i )
      {
         /* try to add new solution to scip */
         SCIP_CALL( copySol(scip, subscip, subvars, subsols[i], &newsol, &success) );
         SCIP_CALL( SCIPtrySol(scip, newsol, TRUE, TRUE, TRUE, &success) );
      }

      if( success )
      {   
         assert(i >= 1);
         *result = SCIP_FOUNDSOL;

         SCIPdebugMessage("undercover heuristic found %d solutions in subMIQCP; solution %d feasible in original problem\n", nsubsols, i);
       
         /* call NLP heuristic for post optimization */
         if( postnlp )
         {
            SCIP_HEUR* nlpheur;
            SCIP_RESULT nlpresult;

            nlpheur = SCIPfindHeur(scip, "nlp");

            if( nlpheur != NULL )
            {
               timelimit -= SCIPgetTotalTime(subscip); 
               if( timelimit >= 10.0 )
               {
                  SCIPdebugMessage("undercover heuristic calling NLP heuristic for post optimization\n");

                  SCIP_CALL( SCIPapplyNlpHeur(scip, nlpheur, &nlpresult, newsol, INT_MAX, timelimit, NULL) );

                  SCIPdebugMessage("NLP heuristic called by undercover %ssuccessfully\n", nlpresult == SCIP_FOUNDSOL ? "" : "un");
               }
            }
         }
      }
      /* free solution */
      SCIP_CALL( SCIPfreeSol(scip, &newsol) );
   }
   
   SCIPdebugMessage("Terminate undercover heuristic%s.\n", success ? "" : "---did not find a feasible solution");
   
   /* free memory from subproblem */
 TERMINATE:
   for( i = nvars-1; i > 0; --i )
   {
      SCIP_CALL( SCIPreleaseVar(subscip, &(subvars[i])) );
   }
   SCIPfreeBufferArray(scip, &subvars);
   SCIP_CALL( SCIPfree(&subscip) );

   /* free memory from ppc problem */      
 TERMINATEPPC:
   SCIPfreeBufferArray(scip, &ppcsolvals);
   for( i = nvars-1; i > 0; --i )
   {
      SCIP_CALL( SCIPreleaseVar(ppcscip, &(ppcvars[i])) );
   }
   SCIPfreeBufferArray(scip, &ppcvars);
   SCIP_CALL( SCIPfree(&ppcscip) );
     
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

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

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

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

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
   SCIP_CONSHDLR* conshdlrand;        /* constraint handler for and constraints or NULL                      */
   SCIP_CONSHDLR* conshdlrquad;       /* constraint handler for quadratic constraints or NULL                */
   SCIP_CONSHDLR* conshdlrsoc;        /* constraint handler for second order cone constraints or NULL        */
   SCIP_CONSHDLR* conshdlrunivar;     /* constraint handler for univariate definite constraints or NULL      */
   SCIP_CONSHDLR* conshdlrsignpower;  /* constraint handler for signpower constraints or NULL      */
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* if the heuristic is called at the root node, we may want to be called directly after the initial root LP solve */
   if( heurdata->beforecuts && SCIPheurGetFreqofs(heur) == 0 )
      SCIPheurSetTimingmask(heur, SCIP_HEURTIMING_DURINGLPLOOP);

   /* look for nonlinear constraints */
   conshdlrand = SCIPfindConshdlr(scip, "and");
   conshdlrquad = SCIPfindConshdlr(scip, "quadratic");
   conshdlrsoc  = SCIPfindConshdlr(scip, "soc");
   conshdlrunivar = SCIPfindConshdlr(scip, "univardefinite");
   conshdlrsignpower = SCIPfindConshdlr(scip, "signpower");

   /* only run heuristic, when there is at least on nonlinear constraint */
   heurdata->run = FALSE;
   heurdata->run =  heurdata->run || (conshdlrand != NULL && SCIPconshdlrGetNConss(conshdlrand) > 0); 
   heurdata->run =  heurdata->run || (conshdlrquad != NULL && SCIPconshdlrGetNConss(conshdlrquad) > 0); 
   heurdata->run =  heurdata->run || (conshdlrsoc != NULL && SCIPconshdlrGetNConss(conshdlrsoc) > 0);
   heurdata->run =  heurdata->run || (conshdlrunivar != NULL && SCIPconshdlrGetNConss(conshdlrunivar) > 0);
   heurdata->run =  heurdata->run || (conshdlrsignpower != NULL && SCIPconshdlrGetNConss(conshdlrsignpower) > 0);

   if( !heurdata->run )
   {
      SCIPdebugMessage("undercover heuristic will not run for <%s> (no known nonlinear constraints present)\n", SCIPgetProbName(scip));
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolUndercover)
{
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   SCIPheurSetTimingmask(heur, HEUR_TIMING);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecUndercover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;                  /* heuristic's data                    */
   SCIP_Real timelimit;                      /* timelimit for the subproblem        */
   SCIP_Real memorylimit;
   SCIP_Longint nstallnodes;                 /* number of stalling nodes for the subproblem */

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);
   /* TODO: assert(SCIPhasCurrentNodeLP(scip)); this leads to an assert e.g. in bell3a ?????????????????????? */

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

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

   /* only call undercover once at the root */
   if( SCIPgetDepth(scip) == 0 && SCIPheurGetNCalls(heur) > 0 )
      return SCIP_OKAY;
   
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

   /* leave some time for NLP heuristic */
   if( heurdata->postnlp && timelimit >= 100.0 )
      timelimit -= 11.0;

   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(scip, memorylimit) )   
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
   if( timelimit < 10.0 || memorylimit <= 0.0 )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* reset timing, if it was changed temporary (at the root node) */
   if( heurtiming != HEUR_TIMING )
      SCIPheurSetTimingmask(heur, HEUR_TIMING);

   SCIPdebugMessage("calling undercover heuristic for <%s>\n", SCIPgetProbName(scip));

   SCIP_CALL( SCIPapplyUndercover(scip, heur, result, timelimit, memorylimit, heurdata->maxnodes, heurdata->ppcstrat, 
         heurdata->ppcobjquot, heurdata->domred, heurdata->fixandprop, heurdata->backtrack, heurdata->locksrounding, heurdata->onlyconvexify, 
         heurdata->globalbounds, heurdata->minimprove, nstallnodes, heurdata->postnlp) );

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
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, -1.0, 1.0, NULL, NULL) );

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

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/globalbounds",
         "should global bounds on variables be used instead of local bounds at focus node?",
         &heurdata->globalbounds, TRUE, DEFAULT_GLOBALBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/beforecuts",
         "should the heuristic be called at root node before cut separation?",
         &heurdata->beforecuts, TRUE, DEFAULT_BEFORECUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/fixandprop",
         "should undercover fix consecutively and propagate fixings?",
         &heurdata->fixandprop, TRUE, DEFAULT_FIXANDPROP, NULL, NULL) );

SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/"HEUR_NAME"/backtrack", 
         "use one level of backtracking if infeasibility is encountered?",
         &heurdata->backtrack, FALSE, DEFAULT_BACKTRACK, NULL, NULL) );
 
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

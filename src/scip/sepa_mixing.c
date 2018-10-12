/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_mixing.c
 * @brief  mixing inequality separator
 * @author Weikun Chen
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/pub_implics.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_sepa.h"
#include "scip/pub_var.h"
#include "scip/scip_branch.h"
#include "scip/scip_cut.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sepa.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_var.h"
#include "scip/sepa_mixing.h"
#include <string.h>


#define SEPA_NAME              "mixing"
#define SEPA_DESC              "mixing inequality separator"
#define SEPA_PRIORITY               -50
#define SEPA_FREQ                    10
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */
#define DEFAULT_USELOACLBOUNDS    TRUE  /**< should local bounds be used? */

/** separator-specific data for the implied bounds separator */
struct SCIP_SepaData
{
   SCIP_Bool             uselocalbounds;  /**< should local bounds be used? */
};

/*
 * local methods
 */

/** adds the given cut */
static
SCIP_RETCODE addCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_Real*            cutcoefs,           /**< coefficients of active variables in cut */
   int*                  cutinds,            /**< problem indices of variables in cut */
   int                   cutnnz,             /**< number of non-zeros in cut */
   SCIP_Real             cutrhs,             /**< right hand side of cut */
   SCIP_Bool             cutislocal,         /**< is the cut only locally valid? */
   SCIP_Bool*            cutoff,             /**< whether a cutoff has been detected */
   int*                  ncuts               /**< pointer to update number of cuts added */
   )
{
   SCIP_VAR** vars;
   SCIP_ROW* cut;
   char cutname[SCIP_MAXSTRLEN];
   int v;

   assert(ncuts != NULL);
   assert(cutoff != NULL);
   *cutoff = FALSE;
   /* get active problem variables */
   vars = SCIPgetVars(scip);

   /* construct cut name */
   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "mix%d_x%d", SCIPgetNLPs(scip), *ncuts);

   /* create empty cut */
   SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), cutrhs,
                                       cutislocal, FALSE, TRUE) );

   /* cache the row extension and only flush them if the cut gets added */
   SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

   /* collect all non-zero coefficients */
   for( v = 0; v < cutnnz; ++v )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, cut, vars[cutinds[v]], cutcoefs[v]) );
   }
   /* flush all changes before adding the cut */
   SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
   /* set cut rank */
   SCIProwChgRank(cut, 1);
#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, " -> found cut");
      SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
#endif
   if( cutislocal )
   {
      SCIP_CALL( SCIPaddRow(scip, cut, FALSE, cutoff) );
   }
   else
   {
      SCIP_CALL( SCIPaddPoolCut(scip, cut) );
   }
   (*ncuts)++;
   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   return SCIP_OKAY;
}

/** searches and adds implied bound cuts that are violated by the given solution value array */
static
SCIP_RETCODE separateCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real*            solvals,            /**< array with solution values of all problem variables */
   int                   ncontvars,          /**< number of continuous variables */
   SCIP_Bool*            cutoff,             /**< whether a cutoff has been detected */
   int*                  ncuts               /**< pointer to store the number of generated cuts */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_VAR* var;
   SCIP_VAR** vars;
   SCIP_Real* vlbmixcoefs;
   SCIP_Real* vlbmixsols;
   SCIP_Real* vubmixcoefs;
   SCIP_Real* vubmixsols;
   SCIP_Real* cutcoefs;
   SCIP_Real cutrhs;
   int* vlbmixinds;
   int* vubmixinds;
   int* cutinds;
   int* vlbmixsigns;
   int* vubmixsigns;
   int vlbmixsize;
   int vubmixsize;
   int cutnnz;

   int firstnonbinvars;
   int nvars;
   int i;
   int k; 
   sepadata = SCIPsepaGetData(sepa);
   *cutoff = FALSE;
   *ncuts = 0;
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   firstnonbinvars = SCIPgetNBinVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &vlbmixcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vlbmixsols, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vlbmixinds, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vlbmixsigns, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vubmixcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vubmixsols, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vubmixinds, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vubmixsigns, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutinds, nvars) );

   for( i=firstnonbinvars; i<nvars; i++ )
   {
      if( SCIPvarGetProbindex(vars[i]) < 0 )
         continue;
      cutnnz = 0;
      cutrhs = 0;
      vlbmixsize = 0;
      var=vars[i];
//      if( strchr( SCIPvarGetName(var), 'I') != NULL )
//         continue;
      assert( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY );
      SCIP_VAR** vlbvars;
      SCIP_Real* vlbcoefs;
      SCIP_Real* vlbconsts;
      SCIP_VAR** vubvars;
      SCIP_Real* vubcoefs;
      SCIP_Real* vubconsts;
      SCIP_Real maxabscoef;
      SCIP_Real curmaxabscoef;
      SCIP_Real activity;
      SCIP_Real lastcoef;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Bool islocallb;
      SCIP_Bool islocalub;
      SCIP_Bool cutislocal; /* Is it a local cut or global cut */
      int maxabsind;
      int maxabssign;
      int nvlb;
      int nvub;
      int j;
      vlbvars = SCIPvarGetVlbVars(var);
      vlbcoefs = SCIPvarGetVlbCoefs(var);
      vlbconsts = SCIPvarGetVlbConstants(var);
      nvlb = SCIPvarGetNVlbs(var);
      islocallb = FALSE;
      if( nvlb == 0 )
         goto VUB;
      maxabscoef = 0.0;
      maxabsind = -1;
      lb = SCIPvarGetLbGlobal(var);
      if( sepadata->uselocalbounds && lb < SCIPvarGetLbLocal(var) )
      {
         /* This is a lcoal cut */
         islocallb = TRUE;
         lb = SCIPvarGetLbLocal(var);
      }

      /* obtain a new lower bound if possible */
      for( j=0; j < nvlb; j++ )
      {
         SCIP_Real tmplb;
         if( SCIPvarIsBinary(vlbvars[j]) && (SCIPvarGetProbindex(vlbvars[j]) >= 0) )
         {
            tmplb = (vlbcoefs[j] > 0) ? vlbconsts[j] : (vlbconsts[j] + vlbcoefs[j]);
            if( tmplb >= lb )
            {
               islocallb = FALSE;
               lb = tmplb;
            }
         }
      }

      if( SCIPisFeasLT(scip, SCIPvarGetUbLocal(var), lb) )
      {
         /*Never happen: Already done by propagation*/
      }

      for( j=0; j < nvlb; j++  )
      {
         SCIP_Real maxactivity;
         SCIP_Real coef;
         SCIP_Real constant;
         if( SCIPvarIsBinary(vlbvars[j]) && (SCIPvarGetProbindex(vlbvars[j]) >= 0) )
         {
            maxactivity = (vlbcoefs[j] > 0) ? (vlbconsts[j] + vlbcoefs[j]) : vlbconsts[j];
            if( SCIPisFeasLE(scip, maxactivity, lb) )
            {
               /* this implied bounds constraint is redundant */
               continue;
            }
            if( vlbcoefs[j] > 0 )
            {
               coef = maxactivity - lb;
               constant = 0;
               vlbmixsigns[vlbmixsize] = 0;
            }
            else
            {
               coef = lb - maxactivity;
               constant = maxactivity - lb;
               vlbmixsigns[vlbmixsize] = 1;
            }

            vlbmixcoefs[vlbmixsize] = REALABS(coef);
            vlbmixinds[vlbmixsize] = SCIPvarGetProbindex(vlbvars[j]);
            vlbmixsols[vlbmixsize] = (!vlbmixsigns[vlbmixsize]) ? SCIPvarGetLPSol(vlbvars[j]) : (1-SCIPvarGetLPSol(vlbvars[j]));
            if( maxabscoef < vlbmixcoefs[vlbmixsize] )
            {
               maxabscoef = vlbmixcoefs[vlbmixsize];
               maxabsind = vlbmixinds[vlbmixsize];
               maxabssign = vlbmixsigns[vlbmixsize];
            }
            vlbmixsize += 1;
         }
      }
      if( vlbmixsize == 0 )
         goto VUB;
      SCIPsortDownRealRealIntInt(vlbmixsols, vlbmixcoefs, vlbmixinds,  vlbmixsigns, vlbmixsize);
      /* the nonbinary variables */
      activity = -(SCIPvarGetLPSol(var) - lb);
      cutcoefs[cutnnz] = -1;
      cutinds[cutnnz] = SCIPvarGetProbindex(var);
      cutrhs = cutrhs - lb;
      cutnnz++;
      lastcoef = 0;
      for( j=0; j < vlbmixsize; j++  )
      {
         SCIP_Real solval;
         solval = vlbmixsols[j];
         if( activity + solval*(maxabscoef-lastcoef) < 0.0 || SCIPisFeasZero(scip, solval) )
            break;
         else
         {
            if( SCIPisLT(scip, vlbmixcoefs[j], lastcoef) )
               continue;
            else
            {
               activity += (vlbmixcoefs[j]-lastcoef)*solval;

               cutrhs -= vlbmixsigns[j]*(vlbmixcoefs[j] -lastcoef);
               cutcoefs[cutnnz] = !vlbmixsigns[j] ? (vlbmixcoefs[j] - lastcoef) : (lastcoef - vlbmixcoefs[j]);
               cutinds[cutnnz] = vlbmixinds[j];
               cutnnz++;
               lastcoef = vlbmixcoefs[j];
            }
         }
      }
      /* adding the variable with maximal coefficient */
      if( lastcoef + 1e-6 < maxabscoef )
      {
         cutrhs -= maxabssign*(maxabscoef - lastcoef);
         cutcoefs[cutnnz] = !maxabssign ? maxabscoef - lastcoef : lastcoef - maxabscoef;
         cutinds[cutnnz] = maxabsind;
         cutnnz++;
      }
      if( SCIPisEfficacious(scip, activity) )
      {
         SCIP_CALL( addCut(scip, sepa, cutcoefs, cutinds, cutnnz, cutrhs, islocallb, cutoff, ncuts) );
      }

VUB:
      cutnnz = 0;
      cutrhs = 0;
      vubmixsize = 0;

      vubvars = SCIPvarGetVubVars(var);
      vubcoefs = SCIPvarGetVubCoefs(var);
      vubconsts = SCIPvarGetVubConstants(var);
      nvub = SCIPvarGetNVubs(var);
      islocalub = FALSE;
      if( nvub == 0 )
         goto CONFLICT;
      maxabscoef = 0.0;
      maxabsind = -1;
      ub = SCIPvarGetUbGlobal(var);
      if( sepadata->uselocalbounds && ub > SCIPvarGetUbLocal(var) )
      {
         /* This is a lcoal cut */
         islocalub = TRUE;
         ub = SCIPvarGetUbLocal(var);
      }

      /* obtain a new upper bound if possible */
      for( j=0; j < nvub; j++ )
      {
         SCIP_Real tmpub;
         if( SCIPvarIsBinary(vubvars[j]) && (SCIPvarGetProbindex(vubvars[j]) >= 0) )
         {
            tmpub = (vubcoefs[j] < 0) ? vubconsts[j] : (vubconsts[j] + vubcoefs[j]);
            if( tmpub <= ub )
            {
               islocalub = FALSE;
               ub = tmpub;
            }
         }
      }


      for( j=0; j < nvub; j++  )
      {
         SCIP_Real minactivity;
         SCIP_Real coef;
         SCIP_Real constant;
         if( SCIPvarIsBinary(vubvars[j]) && (SCIPvarGetProbindex(vubvars[j]) >= 0) )
         {
            minactivity = (vubcoefs[j] < 0) ? (vubconsts[j] + vubcoefs[j]) : vubconsts[j];
            if( SCIPisFeasLE(scip, ub, minactivity) )
            {
               /* this implied bounds constraint is redundant */
               continue;
            }
            if( vubcoefs[j] > 0 )
            {
               coef = ub - minactivity;
               constant = minactivity;
               vubmixsigns[vubmixsize] = 1;
            }
            else
            {
               coef = minactivity - ub;
               constant = ub;
               vubmixsigns[vubmixsize] = 0;
            }

            vubmixcoefs[vubmixsize] = REALABS(coef);
            vubmixinds[vubmixsize] = SCIPvarGetProbindex(vubvars[j]);
            vubmixsols[vubmixsize] = (!vubmixsigns[vubmixsize])? SCIPvarGetLPSol(vubvars[j]):1-SCIPvarGetLPSol(vubvars[j]);
            if( maxabscoef < vubmixcoefs[vubmixsize] )
            {
               maxabscoef = vubmixcoefs[vubmixsize];
               maxabsind = vubmixinds[vubmixsize];
               maxabssign = vubmixsigns[vubmixsize];
            }
            vubmixsize += 1;
         }
      }
      if( vubmixsize == 0 )
         goto CONFLICT;
      SCIPsortDownRealRealIntInt(vubmixsols, vubmixcoefs, vubmixinds,  vubmixsigns, vubmixsize);
      /* the nonbinary variables */
      activity = SCIPvarGetLPSol(var) - ub;
      cutcoefs[cutnnz] = 1;
      cutinds[cutnnz] = SCIPvarGetProbindex(var);
      cutrhs = ub;
      cutnnz++;
      lastcoef = 0;
      for( j=0; j < vubmixsize; j++  )
      {
         SCIP_Real solval;
         solval = vubmixsols[j];
         if( activity + solval*(maxabscoef-lastcoef) < 0.0 || SCIPisFeasZero(scip, solval) )
            break;
         else
         {
            if( SCIPisLT(scip, vubmixcoefs[j], lastcoef) )
               continue;
            else
            {
               activity += (vubmixcoefs[j]-lastcoef)*solval;

               cutrhs -= vubmixsigns[j]*(vubmixcoefs[j] -lastcoef);
               cutcoefs[cutnnz] = !vubmixsigns[j] ? vubmixcoefs[j] - lastcoef : lastcoef - vubmixcoefs[j];
               cutinds[cutnnz] = vubmixinds[j];
               cutnnz++;
               lastcoef = vubmixcoefs[j];
            }
         }
      }
      /* adding the variable with maximal coefficient */
      if( lastcoef + 1e-6 < maxabscoef )
      {
         cutrhs -= maxabssign*(maxabscoef - lastcoef);
         cutcoefs[cutnnz] = !maxabssign ? maxabscoef - lastcoef : lastcoef - maxabscoef;
         cutinds[cutnnz] = maxabsind;
         cutnnz++;
      }
      if( SCIPisEfficacious(scip, activity) )
      {
         SCIP_CALL( addCut(scip, sepa, cutcoefs, cutinds, cutnnz, cutrhs, islocalub, cutoff, ncuts) );
      }
CONFLICT:

      if( (vlbmixsize == 0) || (vubmixsize == 0) )
         continue;
      cutislocal = islocallb || islocalub;
      for( j=0; j<vlbmixsize; j++ )
      {
         SCIP_Real solval;
         solval = vlbmixsols[j];
         if( !SCIPisEfficacious(scip, solval+vubmixsols[0]-1) )
            break;
         for( k=0; k<vubmixsize; k++ )
         {
            if( SCIPisEfficacious(scip, solval+vubmixsols[k]-1) )
            {
               SCIP_Real tmp;
               tmp = lb + vlbmixcoefs[j] + vubmixcoefs[k] - ub;
               if( SCIPisEfficacious(scip, tmp) )
               {
                  cutnnz = 2;
                  cutrhs = 1.0;
                  cutcoefs[0] = vlbmixsigns[j] ? -1.0 : 1.0;
                  cutcoefs[1] = vubmixsigns[k] ? -1.0 : 1.0;
                  cutinds[0] = vlbmixinds[j];
                  cutinds[1] = vubmixinds[k];
                  cutrhs = vlbmixsigns[j] ? (cutrhs - 1.0) : cutrhs;
                  cutrhs = vubmixsigns[k] ? (cutrhs - 1.0) : cutrhs;
                  SCIP_CALL( addCut(scip, sepa, cutcoefs, cutinds, cutnnz, cutrhs, cutislocal, cutoff, ncuts) );
               }
            }
            else
               break;
         }
      }
   }
   SCIPfreeBufferArray(scip, &vlbmixcoefs);
   SCIPfreeBufferArray(scip, &vlbmixsols);
   SCIPfreeBufferArray(scip, &vlbmixinds);
   SCIPfreeBufferArray(scip, &vlbmixsigns);
   SCIPfreeBufferArray(scip, &vubmixcoefs);
   SCIPfreeBufferArray(scip, &vubmixsols);
   SCIPfreeBufferArray(scip, &vubmixinds);
   SCIPfreeBufferArray(scip, &vubmixsigns);
   SCIPfreeBufferArray(scip, &cutcoefs);
   SCIPfreeBufferArray(scip, &cutinds);

   return SCIP_OKAY;
}


/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyMixing)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaMixing(scip) );

   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeMixing)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* get separation data and free it */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   SCIPfreeBlockMemory(scip, &sepadata);

   /* reset data pointer to NULL */
   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpMixing)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_Real* solvals;
   SCIP_Bool cutoff;
   int ncontvars;
   int nvars;
   int nfracs;
   int ncuts;

   assert(sepa != NULL);
   assert(scip != NULL);

   *result = SCIP_DIDNOTRUN;

   /* gets active problem variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, &ncontvars) );
   if( ncontvars == 0 )
      return SCIP_OKAY;


   /* get solution values for all variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars) );
   SCIP_CALL( SCIPgetVarSols(scip, nvars, vars, solvals) );

   /* call the cut separation */
   SCIP_CALL( separateCuts(scip, sepa, NULL, solvals, ncontvars, &cutoff, &ncuts) );

   /* adjust result code */
   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if ( ncuts > 0 )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_DIDNOTFIND;

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &solvals);

   return SCIP_OKAY;
}


/** arbitrary primal solution separation method of separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolMixing)
{  /*lint --e{715}*/
#if 0
   SCIP_VAR** vars;
   SCIP_VAR** fracvars;
   SCIP_Real* solvals;
   SCIP_Real* fracvals;
   SCIP_Bool cutoff;
   int nvars;
   int nbinvars;
   int nfracs;
   int ncuts;
   int i;

   assert(sepa != NULL);
   assert(scip != NULL);

   *result = SCIP_DIDNOTRUN;

   /* gets active problem variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, NULL, NULL, NULL) );
   if( nbinvars == 0 )
      return SCIP_OKAY;

   /* get solution values for all variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars) );
   SCIP_CALL( SCIPgetSolVals(scip, sol, nvars, vars, solvals) );

   /* get binary problem variables that are fractional in given solution */
   SCIP_CALL( SCIPallocBufferArray(scip, &fracvars, nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fracvals, nbinvars) );
   nfracs = 0;
   for( i = 0; i < nbinvars; ++i )
   {
      if( !SCIPisFeasIntegral(scip, solvals[i]) )
      {
         fracvars[nfracs] = vars[i];
         fracvals[nfracs] = solvals[i];
         nfracs++;
      }
   }

   /* call the cut separation */
   ncuts = 0;
   cutoff = FALSE;

   if( nfracs > 0 )
   {
      SCIP_CALL( separateCuts(scip, sepa, sol, solvals, ncontvars, &cutoff, &ncuts) );
   }

   /* adjust result code */
   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if ( ncuts > 0 )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_DIDNOTFIND;

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &fracvals);
   SCIPfreeBufferArray(scip, &fracvars);
   SCIPfreeBufferArray(scip, &solvals);
#endif
   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the mixing separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaMixing(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create mixing separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   assert(sepadata != NULL);

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpMixing, sepaExecsolMixing,
         sepadata) );
   assert(sepa != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyMixing) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeMixing) );

   /* add separator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/mixing/uselocalbounds",
         "should local bounds be used?",
         &sepadata->uselocalbounds, TRUE, DEFAULT_USELOACLBOUNDS, NULL, NULL) );

   return SCIP_OKAY;
}

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
 * @brief  mixing/star inequality separator
 * @author Weikun Chen
 *
 * For an overview see:
 * Atamturk, A., Nemhauser, G.L. and Savelsbergh, M.W., 2000@n
 * The mixed vertex packing problem.@n
 * Mathematical Programming, 89(1), 35-53.
 *
 * Some remarks:
 * - Besides the mixing inequality, we also add the conflict inequality.
 * - Currently, we found the performance is really bad on the neos-565672 instances.
 *   The reason is that, with the mixing inequality, SCIP spends a lot of time on cutting planes generations.
 * - Currently, we do not consider sparsity of the cuts as we only care about finding a most violated cut.
 * - Besides the most violated cut we consider, we also add an additional variable to make the cut be the strongest one,
 *   even the additional variable does not contribute any violation.
 *
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
#include "scip/scip_tree.h"
#include "scip/sepa_mixing.h"
#include <string.h>


#define SEPA_NAME              "mixing"
#define SEPA_DESC              "mixing inequality separator"
#define DEFAULT_MAXROUNDS            -1 /**< maximal number of mixing separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        -1 /**< maximal number of mixing separation rounds in the root node (-1: unlimited) */
#define SEPA_PRIORITY               -50
#define SEPA_FREQ                    10
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */
#define DEFAULT_USELOACLBOUNDS    FALSE /**< should local bounds be used? */
#define DEFAULT_ISCUTSONINTS      FALSE /**< should general integer variables be used to generate cuts? */

/** separator-specific data for the implied bounds separator */
struct SCIP_SepaData
{
   SCIP_Bool             uselocalbounds;     /**< should local bounds be used? */
   SCIP_Bool             iscutsonints;       /**< should general integer variables be used to generate cuts? */
   int                   maxrounds;          /**< maximal number of mixing separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of mixing separation rounds in the root node (-1: unlimited) */
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
   assert(cutcoefs != NULL);
   assert(cutinds != NULL);
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
      /* local cuts we add to the sepastore */
      SCIP_CALL( SCIPaddRow(scip, cut, FALSE, cutoff) );
   }
   else
   {
      SCIP_CALL( SCIPaddPoolCut(scip, cut) );
   }
   (*ncuts)++;
   /* release the row */
   SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   return SCIP_OKAY;
}

/** searches and adds mixing cuts that are violated by the given solution value array */
static
SCIP_RETCODE separateCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
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
   int firstvars;
   int nvars;
   int i;
   int k;
   assert(cutoff != NULL);
   assert(ncuts != NULL);

   *cutoff = FALSE;
   *ncuts = 0;
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   /* get the index of the first considered variable */
   if( sepadata->iscutsonints )
   {
      /* generate cuts based on all nonbinary variabls */
      firstvars = SCIPgetNBinVars(scip);
   }
   else
   {
      /* only generate cuts based on continuous variables */
      firstvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) + SCIPgetNImplVars(scip);
   }
   /* allocate temporary memory */
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

   for( i=firstvars; i<nvars; i++ )
   {
      SCIP_VAR** vlbvars;
      SCIP_Real* vlbcoefs;
      SCIP_Real* vlbconsts;
      SCIP_VAR** vubvars;
      SCIP_Real* vubcoefs;
      SCIP_Real* vubconsts;
      SCIP_Real maxabscoef;
      SCIP_Real activity;
      SCIP_Real lastcoef;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Bool islocallb;  /* Is it a local lower bound or global lower bound */
      SCIP_Bool islocalub;  /* Is it a local upper bound or global upper bound */
      SCIP_Bool cutislocal; /* Is it a local cut or global cut */
      int maxabsind;
      int maxabssign;
      int nvlb;
      int nvub;
      int j;

      if( SCIPvarGetProbindex(vars[i]) < 0 )
         continue;
      cutnnz = 0;
      cutrhs = 0;
      vlbmixsize = 0;
      var=vars[i];
      assert( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY );
      /* get variable lower bounds information */
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
      if( SCIPisFeasEQ(scip, SCIPvarGetUbLocal(var), SCIPvarGetLPSol(var)) )
         goto VUB;

      /* obtain a new lower (and global) bound if possible */
      for( j=0; j < nvlb; j++ )
      {
         SCIP_Real tmplb;
         if( SCIPvarIsBinary(vlbvars[j]) && (SCIPvarGetProbindex(vlbvars[j]) >= 0) )
         {
            tmplb = (vlbcoefs[j] > 0) ? vlbconsts[j] : (vlbconsts[j] + vlbcoefs[j]);
            if( tmplb >= lb )
            {
               /* variable bounds information is global, as a result, the cut should be global at this point */
               islocallb = FALSE;
               lb = tmplb;
            }
         }
      }

      if( SCIPisFeasLT(scip, SCIPvarGetUbLocal(var), lb) )
      {
         /*Never happen: already done by propagation*/
      }

      /* extract the useful variable bounds information (binary and nonredundant) */
      for( j=0; j < nvlb; j++  )
      {
         /* consider only active and binary variable */
         if( SCIPvarIsBinary(vlbvars[j]) && (SCIPvarGetProbindex(vlbvars[j]) >= 0) )
         {
            SCIP_Real maxactivity;
            SCIP_Real coef;
            maxactivity = (vlbcoefs[j] > 0) ? (vlbconsts[j] + vlbcoefs[j]) : vlbconsts[j];
            coef = 0.0;
            if( SCIPisFeasLE(scip, maxactivity, lb) )
            {
               /* this variable bound constraint is redundant */
               continue;
            }
            if( vlbcoefs[j] > 0 )
            {
               coef = maxactivity - lb;
               vlbmixsigns[vlbmixsize] = 0.0;
            }
            else
            {
               coef = lb - maxactivity;
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
      /* stop if no variable lower bounds information exists */
      if( vlbmixsize == 0 )
         goto VUB;
      /* stop if the current solution value is larger than the maxial coefficient */
      if( SCIPisFeasGT(scip, (SCIPvarGetLPSol(var) - lb), maxabscoef) )
         goto VUB;
      /* sort the lp solutions decreasingly */
      SCIPsortDownRealRealIntInt(vlbmixsols, vlbmixcoefs, vlbmixinds,  vlbmixsigns, vlbmixsize);
      /* add the continuous variable */
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
         /* stop if we could not find a violated cut */
         if( activity + solval*(maxabscoef-lastcoef) < 0.0 || SCIPisFeasZero(scip, solval) )
            break;
         else
         {
            /* skip if we have already add a variable with bigger coefficients */
            if( SCIPisLE(scip, vlbmixcoefs[j], lastcoef) )
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
      /* adding the variable with maximal coefficient to make sure the cut is stronger enough */
      if( SCIPisGT(scip, maxabscoef, lastcoef) )
      {
         cutrhs -= maxabssign*(maxabscoef - lastcoef);
         cutcoefs[cutnnz] = !maxabssign ? maxabscoef - lastcoef : lastcoef - maxabscoef;
         cutinds[cutnnz] = maxabsind;
         cutnnz++;
      }
      if( SCIPisEfficacious(scip, activity) )
      {
         /* add the cut if the violtion is good enough */
         SCIP_CALL( addCut(scip, sepa, cutcoefs, cutinds, cutnnz, cutrhs, islocallb, cutoff, ncuts) );
      }

VUB:
      cutnnz = 0;
      cutrhs = 0;
      vubmixsize = 0;

      /* get variable upper bounds information */
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
      if( SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetLPSol(var)) )
         goto CONFLICT;
      if( sepadata->uselocalbounds && ub > SCIPvarGetUbLocal(var) )
      {
         /* This is a lcoal cut */
         islocalub = TRUE;
         ub = SCIPvarGetUbLocal(var);
      }

      /* obtain a new upper bound (and global) if possible */
      for( j=0; j < nvub; j++ )
      {
         SCIP_Real tmpub;
         if( SCIPvarIsBinary(vubvars[j]) && (SCIPvarGetProbindex(vubvars[j]) >= 0) )
         {
            tmpub = (vubcoefs[j] < 0) ? vubconsts[j] : (vubconsts[j] + vubcoefs[j]);
            if( tmpub <= ub )
            {
               /* variable bounds information is global, as a result, the cut should be global at this point */
               islocalub = FALSE;
               ub = tmpub;
            }
         }
      }


      /* extract the useful variable bounds information (binary and nonredundant) */
      for( j=0; j < nvub; j++  )
      {
         /* consider only active and binary variable */
         if( SCIPvarIsBinary(vubvars[j]) && (SCIPvarGetProbindex(vubvars[j]) >= 0) )
         {
            SCIP_Real minactivity;
            SCIP_Real coef;
            minactivity = (vubcoefs[j] < 0) ? (vubconsts[j] + vubcoefs[j]) : vubconsts[j];
            coef = 0.0;
            if( SCIPisFeasLE(scip, ub, minactivity) )
            {
               /* this variable bound constraint is redundant */
               continue;
            }
            if( vubcoefs[j] > 0 )
            {
               coef = ub - minactivity;
               vubmixsigns[vubmixsize] = 1;
            }
            else
            {
               coef = minactivity - ub;
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
      /* stop if no variable upper bounds information exists */
      if( vubmixsize == 0 )
         goto CONFLICT;
      /* stop if the current solution value is larger than the maxial coefficient */
      if( SCIPisFeasGT(scip, (ub - SCIPvarGetLPSol(var)), maxabscoef) )
         goto CONFLICT;
      /* sort the lp solutions decreasingly */
      SCIPsortDownRealRealIntInt(vubmixsols, vubmixcoefs, vubmixinds,  vubmixsigns, vubmixsize);
      /* add the continuous variables */
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
            if( SCIPisLE(scip, vubmixcoefs[j], lastcoef) )
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
      if( SCIPisGT(scip, maxabscoef, lastcoef) )
      {
         cutrhs -= maxabssign*(maxabscoef - lastcoef);
         cutcoefs[cutnnz] = !maxabssign ? maxabscoef - lastcoef : lastcoef - maxabscoef;
         cutinds[cutnnz] = maxabsind;
         cutnnz++;
      }
      if( SCIPisEfficacious(scip, activity) )
      {
         /* add the cut if the violtion is good enough */
         SCIP_CALL( addCut(scip, sepa, cutcoefs, cutinds, cutnnz, cutrhs, islocalub, cutoff, ncuts) );
      }

      /* combine the variable lower bounds information and upper bounds information together to generate cuts */
CONFLICT:
      /* stop if no useful variable lower (or upper) bounds information exists */
      if( (vlbmixsize == 0) || (vubmixsize == 0) )
         continue;
      cutislocal = islocallb || islocalub;
      for( j=0; j<vlbmixsize; j++ )
      {
         SCIP_Real solval;
         solval = vlbmixsols[j];
         /* stop if no violated cut exists */
         if( !SCIPisEfficacious(scip, solval+vubmixsols[0]-1) )
            break;
         for( k=0; k<vubmixsize; k++ )
         {
            /* add the cut if the violation is good enough */
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
   /* free temporary memory */
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
   SCIP_SEPADATA* sepadata;
   SCIP_Bool cutoff;
   int nbinvars;
   int nvars;
   int ncuts;
   int ncalls;
   int depth;

   assert(sepa != NULL);
   assert(scip != NULL);

   *result = SCIP_DIDNOTRUN;
   depth = SCIPgetDepth(scip);
   ncalls = SCIPsepaGetNCallsAtNode(sepa);
   sepadata = SCIPsepaGetData(sepa);
   /* only call the mixing cut separator a given number of times at each node */
   if( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;

   /* gets numver of active problem variables and number of binary variables */
   SCIP_CALL( SCIPgetVarsData(scip, NULL, &nvars, &nbinvars, NULL, NULL, NULL) );
   /* if all the active problem variables are binary, stop  */
   if( nvars == nbinvars )
      return SCIP_OKAY;



   /* call the cut separation */
   SCIP_CALL( separateCuts(scip, sepa, NULL, &cutoff, &ncuts) );

   /* adjust result code */
   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if ( ncuts > 0 )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_DIDNOTFIND;


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
         sepaExeclpMixing, NULL,
         sepadata) );
   assert(sepa != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyMixing) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeMixing) );

   /* add separator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/mixing/uselocalbounds",
         "should local bounds be used?", &sepadata->uselocalbounds, TRUE,
         DEFAULT_USELOACLBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/mixing/iscutsonints",
         "hould general integer variables be used to generate cuts?", &sepadata->iscutsonints, TRUE,
         DEFAULT_ISCUTSONINTS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/mixing/maxrounds",
         "maximal number of mixing separation rounds per node (-1: unlimited)",
         &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/mixing/maxroundsroot",
         "maximal number of mixing separation rounds in the root node (-1: unlimited)",
         &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   return SCIP_OKAY;
}

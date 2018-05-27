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

/**@file   scip_nonlinear.c
 * @brief  public methods for nonlinear functions
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Robert Lion Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <ctype.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/
#endif


#include "lpi/lpi.h"
#include "nlpi/exprinterpret.h"
#include "nlpi/nlpi.h"
#include "scip/benders.h"
#include "scip/benderscut.h"
#include "scip/branch.h"
#include "scip/branch_nodereopt.h"
#include "scip/clock.h"
#include "scip/compr.h"
#include "scip/concsolver.h"
#include "scip/concurrent.h"
#include "scip/conflict.h"
#include "scip/conflictstore.h"
#include "scip/cons.h"
#include "scip/cons_linear.h"
#include "scip/cutpool.h"
#include "scip/cuts.h"
#include "scip/debug.h"
#include "scip/def.h"
#include "scip/dialog.h"
#include "scip/dialog_default.h"
#include "scip/disp.h"
#include "scip/event.h"
#include "scip/heur.h"
#include "scip/heur_ofins.h"
#include "scip/heur_reoptsols.h"
#include "scip/heur_trivialnegation.h"
#include "scip/heuristics.h"
#include "scip/history.h"
#include "scip/implics.h"
#include "scip/interrupt.h"
#include "scip/lp.h"
#include "scip/mem.h"
#include "scip/message_default.h"
#include "scip/misc.h"
#include "scip/nlp.h"
#include "scip/nodesel.h"
#include "scip/paramset.h"
#include "scip/presol.h"
#include "scip/presolve.h"
#include "scip/pricer.h"
#include "scip/pricestore.h"
#include "scip/primal.h"
#include "scip/prob.h"
#include "scip/prop.h"
#include "scip/reader.h"
#include "scip/relax.h"
#include "scip/reopt.h"
#include "scip/retcode.h"
#include "scip/scipbuildflags.h"
#include "scip/scipcoreplugins.h"
#include "scip/scipgithash.h"
#include "scip/sepa.h"
#include "scip/sepastore.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/solve.h"
#include "scip/stat.h"
#include "scip/syncstore.h"
#include "scip/table.h"
#include "scip/tree.h"
#include "scip/var.h"
#include "scip/visual.h"
#include "xml/xml.h"

#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nonlinear.h"
#include "scip/scip_numerics.h"
#include "scip/scip_prob.h"

#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_nlp.h"
#include "scip/pub_var.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif

/** computes coefficients of linearization of a square term in a reference point */
void SCIPaddSquareLinearization(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             sqrcoef,            /**< coefficient of square term */
   SCIP_Real             refpoint,           /**< point where to linearize */
   SCIP_Bool             isint,              /**< whether corresponding variable is a discrete variable, and thus linearization could be moved */
   SCIP_Real*            lincoef,            /**< buffer to add coefficient of linearization */
   SCIP_Real*            linconstant,        /**< buffer to add constant of linearization */
   SCIP_Bool*            success             /**< buffer to set to FALSE if linearization has failed due to large numbers */
   )
{
   assert(scip != NULL);
   assert(lincoef != NULL);
   assert(linconstant != NULL);
   assert(success != NULL);

   if( sqrcoef == 0.0 )
      return;

   if( SCIPisInfinity(scip, REALABS(refpoint)) )
   {
      *success = FALSE;
      return;
   }

   if( !isint || SCIPisIntegral(scip, refpoint) )
   {
      SCIP_Real tmp;

      /* sqrcoef * x^2  ->  tangent in refpoint = sqrcoef * 2 * refpoint * (x - refpoint) */

      tmp = sqrcoef * refpoint;

      if( SCIPisInfinity(scip, 2.0 * REALABS(tmp)) )
      {
         *success = FALSE;
         return;
      }

      *lincoef += 2.0 * tmp;
      tmp *= refpoint;
      *linconstant -= tmp;
   }
   else
   {
      /* sqrcoef * x^2 -> secant between f=floor(refpoint) and f+1 = sqrcoef * (f^2 + ((f+1)^2 - f^2) * (x-f))
       * = sqrcoef * (-f*(f+1) + (2*f+1)*x)
       */
      SCIP_Real f;
      SCIP_Real coef;
      SCIP_Real constant;

      f = SCIPfloor(scip, refpoint);

      coef     =  sqrcoef * (2.0 * f + 1.0);
      constant = -sqrcoef * f * (f + 1.0);

      if( SCIPisInfinity(scip, REALABS(coef)) || SCIPisInfinity(scip, REALABS(constant)) )
      {
         *success = FALSE;
         return;
      }

      *lincoef     += coef;
      *linconstant += constant;
   }
}

/** computes coefficients of secant of a square term */
void SCIPaddSquareSecant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             sqrcoef,            /**< coefficient of square term */
   SCIP_Real             lb,                 /**< lower bound on variable */
   SCIP_Real             ub,                 /**< upper bound on variable */
   SCIP_Real             refpoint,           /**< point for which to compute value of linearization */
   SCIP_Real*            lincoef,            /**< buffer to add coefficient of secant */
   SCIP_Real*            linconstant,        /**< buffer to add constant of secant */
   SCIP_Bool*            success             /**< buffer to set to FALSE if secant has failed due to large numbers or unboundedness */
   )
{
   SCIP_Real coef;
   SCIP_Real constant;

   assert(scip != NULL);
   assert(!SCIPisInfinity(scip,  lb));
   assert(!SCIPisInfinity(scip, -ub));
   assert(SCIPisLE(scip, lb, ub));
   assert(SCIPisLE(scip, lb, refpoint));
   assert(SCIPisGE(scip, ub, refpoint));
   assert(lincoef != NULL);
   assert(linconstant != NULL);
   assert(success != NULL);

   if( sqrcoef == 0.0 )
      return;

   if( SCIPisInfinity(scip, -lb) || SCIPisInfinity(scip, ub) )
   {
      /* unboundedness */
      *success = FALSE;
      return;
   }

   /* sqrcoef * x^2 -> sqrcoef * (lb * lb + (ub*ub - lb*lb)/(ub-lb) * (x-lb)) = sqrcoef * (lb*lb + (ub+lb)*(x-lb))
    *  = sqrcoef * ((lb+ub)*x - lb*ub)
    */
   coef     =  sqrcoef * (lb + ub);
   constant = -sqrcoef * lb * ub;
   if( SCIPisInfinity(scip, REALABS(coef)) || SCIPisInfinity(scip, REALABS(constant)) )
   {
      *success = FALSE;
      return;
   }

   *lincoef     += coef;
   *linconstant += constant;
}

/** computes coefficients of linearization of a bilinear term in a reference point */
void SCIPaddBilinLinearization(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bilincoef,          /**< coefficient of bilinear term */
   SCIP_Real             refpointx,          /**< point where to linearize first  variable */
   SCIP_Real             refpointy,          /**< point where to linearize second variable */
   SCIP_Real*            lincoefx,           /**< buffer to add coefficient of first  variable in linearization */
   SCIP_Real*            lincoefy,           /**< buffer to add coefficient of second variable in linearization */
   SCIP_Real*            linconstant,        /**< buffer to add constant of linearization */
   SCIP_Bool*            success             /**< buffer to set to FALSE if linearization has failed due to large numbers */
   )
{
   SCIP_Real constant;

   assert(scip != NULL);
   assert(lincoefx != NULL);
   assert(lincoefy != NULL);
   assert(linconstant != NULL);
   assert(success != NULL);

   if( bilincoef == 0.0 )
      return;

   if( SCIPisInfinity(scip, REALABS(refpointx)) || SCIPisInfinity(scip, REALABS(refpointy)) )
   {
      *success = FALSE;
      return;
   }

   /* bilincoef * x * y ->  bilincoef * (refpointx * refpointy + refpointy * (x - refpointx) + refpointx * (y - refpointy))
    *                    = -bilincoef * refpointx * refpointy + bilincoef * refpointy * x + bilincoef * refpointx * y
    */

   constant = -bilincoef * refpointx * refpointy;

   if( SCIPisInfinity(scip, REALABS(bilincoef * refpointx)) || SCIPisInfinity(scip, REALABS(bilincoef * refpointy))
      || SCIPisInfinity(scip, REALABS(constant)) )
   {
      *success = FALSE;
      return;
   }

   *lincoefx    += bilincoef * refpointy;
   *lincoefy    += bilincoef * refpointx;
   *linconstant += constant;
}

/** computes coefficients of McCormick under- or overestimation of a bilinear term */
void SCIPaddBilinMcCormick(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bilincoef,          /**< coefficient of bilinear term */
   SCIP_Real             lbx,                /**< lower bound on first variable */
   SCIP_Real             ubx,                /**< upper bound on first variable */
   SCIP_Real             refpointx,          /**< reference point for first variable */
   SCIP_Real             lby,                /**< lower bound on second variable */
   SCIP_Real             uby,                /**< upper bound on second variable */
   SCIP_Real             refpointy,          /**< reference point for second variable */
   SCIP_Bool             overestimate,       /**< whether to compute an overestimator instead of an underestimator */
   SCIP_Real*            lincoefx,           /**< buffer to add coefficient of first  variable in linearization */
   SCIP_Real*            lincoefy,           /**< buffer to add coefficient of second variable in linearization */
   SCIP_Real*            linconstant,        /**< buffer to add constant of linearization */
   SCIP_Bool*            success             /**< buffer to set to FALSE if linearization has failed due to large numbers */
   )
{
   SCIP_Real constant;
   SCIP_Real coefx;
   SCIP_Real coefy;

   assert(scip != NULL);
   assert(!SCIPisInfinity(scip,  lbx));
   assert(!SCIPisInfinity(scip, -ubx));
   assert(!SCIPisInfinity(scip,  lby));
   assert(!SCIPisInfinity(scip, -uby));
   assert(SCIPisLE(scip, lbx, ubx));
   assert(SCIPisLE(scip, lby, uby));
   assert(SCIPisLE(scip, lbx, refpointx));
   assert(SCIPisGE(scip, ubx, refpointx));
   assert(SCIPisLE(scip, lby, refpointy));
   assert(SCIPisGE(scip, uby, refpointy));
   assert(lincoefx != NULL);
   assert(lincoefy != NULL);
   assert(linconstant != NULL);
   assert(success != NULL);

   if( bilincoef == 0.0 )
      return;

   if( overestimate )
      bilincoef = -bilincoef;

   if( SCIPisRelEQ(scip, lbx, ubx) && SCIPisRelEQ(scip, lby, uby) )
   {
      /* both x and y are mostly fixed */
      SCIP_Real cand1;
      SCIP_Real cand2;
      SCIP_Real cand3;
      SCIP_Real cand4;

      coefx = 0.0;
      coefy = 0.0;

      /* estimate x * y by constant */
      cand1 = lbx * lby;
      cand2 = lbx * uby;
      cand3 = ubx * lby;
      cand4 = ubx * uby;

      /* take most conservative value for underestimator */
      if( bilincoef < 0.0 )
         constant = bilincoef * MAX( MAX(cand1, cand2), MAX(cand3, cand4) );
      else
         constant = bilincoef * MIN( MIN(cand1, cand2), MIN(cand3, cand4) );
   }
   else if( bilincoef > 0.0 )
   {
      /* either x or y is not fixed and coef > 0.0 */
      if( !SCIPisInfinity(scip, -lbx) && !SCIPisInfinity(scip, -lby) &&
         (SCIPisInfinity(scip,  ubx) || SCIPisInfinity(scip,  uby)
            || (uby - refpointy) * (ubx - refpointx) >= (refpointy - lby) * (refpointx - lbx)) )
      {
         if( SCIPisRelEQ(scip, lbx, ubx) )
         {
            /* x*y = lbx * y + (x-lbx) * y >= lbx * y + (x-lbx) * lby >= lbx * y + min{(ubx-lbx) * lby, 0 * lby} */
            coefx    =  0.0;
            coefy    =  bilincoef * lbx;
            constant =  bilincoef * (lby < 0.0 ? (ubx-lbx) * lby : 0.0);
         }
         else if( SCIPisRelEQ(scip, lby, uby) )
         {
            /* x*y = lby * x + (y-lby) * x >= lby * x + (y-lby) * lbx >= lby * x + min{(uby-lby) * lbx, 0 * lbx} */
            coefx    =  bilincoef * lby;
            coefy    =  0.0;
            constant =  bilincoef * (lbx < 0.0 ? (uby-lby) * lbx : 0.0);
         }
         else
         {
            coefx    =  bilincoef * lby;
            coefy    =  bilincoef * lbx;
            constant = -bilincoef * lbx * lby;
         }
      }
      else if( !SCIPisInfinity(scip, ubx) && !SCIPisInfinity(scip, uby) )
      {
         if( SCIPisRelEQ(scip, lbx, ubx) )
         {
            /* x*y = ubx * y + (x-ubx) * y >= ubx * y + (x-ubx) * uby >= ubx * y + min{(lbx-ubx) * uby, 0 * uby} */
            coefx    =  0.0;
            coefy    =  bilincoef * ubx;
            constant =  bilincoef * (uby > 0.0 ? (lbx - ubx) * uby : 0.0);
         }
         else if( SCIPisRelEQ(scip, lby, uby) )
         {
            /* x*y = uby * x + (y-uby) * x >= uby * x + (y-uby) * ubx >= uby * x + min{(lby-uby) * ubx, 0 * ubx} */
            coefx    =  bilincoef * uby;
            coefy    =  0.0;
            constant =  bilincoef * (ubx > 0.0 ? (lby - uby) * ubx : 0.0);
         }
         else
         {
            coefx    =  bilincoef * uby;
            coefy    =  bilincoef * ubx;
            constant = -bilincoef * ubx * uby;
         }
      }
      else
      {
         *success = FALSE;
         return;
      }
   }
   else
   {
      /* either x or y is not fixed and coef < 0.0 */
      if( !SCIPisInfinity(scip,  ubx) && !SCIPisInfinity(scip, -lby) &&
         (SCIPisInfinity(scip, -lbx) || SCIPisInfinity(scip,  uby)
            || (ubx - lbx) * (refpointy - lby) <= (uby - lby) * (refpointx - lbx)) )
      {
         if( SCIPisRelEQ(scip, lbx, ubx) )
         {
            /* x*y = ubx * y + (x-ubx) * y <= ubx * y + (x-ubx) * lby <= ubx * y + max{(lbx-ubx) * lby, 0 * lby} */
            coefx    =  0.0;
            coefy    =  bilincoef * ubx;
            constant =  bilincoef * (lby < 0.0 ? (lbx - ubx) * lby : 0.0);
         }
         else if( SCIPisRelEQ(scip, lby, uby) )
         {
            /* x*y = lby * x + (y-lby) * x <= lby * x + (y-lby) * ubx <= lby * x + max{(uby-lby) * ubx, 0 * ubx} */
            coefx    =  bilincoef * lby;
            coefy    =  0.0;
            constant =  bilincoef * (ubx > 0.0 ? (uby - lby) * ubx : 0.0);
         }
         else
         {
            coefx    =  bilincoef * lby;
            coefy    =  bilincoef * ubx;
            constant = -bilincoef * ubx * lby;
         }
      }
      else if( !SCIPisInfinity(scip, -lbx) && !SCIPisInfinity(scip, uby) )
      {
         if( SCIPisRelEQ(scip, lbx, ubx) )
         {
            /* x*y = lbx * y + (x-lbx) * y <= lbx * y + (x-lbx) * uby <= lbx * y + max{(ubx-lbx) * uby, 0 * uby} */
            coefx    =  0.0;
            coefy    =  bilincoef * lbx;
            constant =  bilincoef * (uby > 0.0 ? (ubx - lbx) * uby : 0.0);
         }
         else if( SCIPisRelEQ(scip, lby, uby) )
         {
            /* x*y = uby * x + (y-uby) * x <= uby * x + (y-uby) * lbx <= uby * x + max{(lby-uby) * lbx, 0 * lbx} */
            coefx    =  bilincoef * uby;
            coefy    =  0.0;
            constant =  bilincoef * (lbx < 0.0 ? (lby - uby) * lbx : 0.0);
         }
         else
         {
            coefx    =  bilincoef * uby;
            coefy    =  bilincoef * lbx;
            constant = -bilincoef * lbx * uby;
         }
      }
      else
      {
         *success = FALSE;
         return;
      }
   }

   if( SCIPisInfinity(scip, REALABS(coefx)) || SCIPisInfinity(scip, REALABS(coefy))
      || SCIPisInfinity(scip, REALABS(constant)) )
   {
      *success = FALSE;
      return;
   }

   if( overestimate )
   {
      coefx    = -coefx;
      coefy    = -coefy;
      constant = -constant;
   }

   SCIPdebugMsg(scip, "%.20g * x[%.20g,%.20g] * y[%.20g,%.20g] %c= %.20g * x %+.20g * y %+.20g\n", bilincoef, lbx, ubx,
      lby, uby, overestimate ? '<' : '>', coefx, coefy, constant);

   *lincoefx    += coefx;
   *lincoefy    += coefy;
   *linconstant += constant;
}


/** computes coefficients of linearization of a bilinear term in a reference point when given a linear inequality
 *  involving only the variables of the bilinear term
 *
 *  @note the formulas are extracted from "Convex envelopes of bivariate functions through the solution of KKT systems"
 *        by Marco Locatelli
 */
void SCIPcomputeBilinEnvelope1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bilincoef,          /**< coefficient of bilinear term */
   SCIP_Real             lbx,                /**< lower bound on first variable */
   SCIP_Real             ubx,                /**< upper bound on first variable */
   SCIP_Real             refpointx,          /**< reference point for first variable */
   SCIP_Real             lby,                /**< lower bound on second variable */
   SCIP_Real             uby,                /**< upper bound on second variable */
   SCIP_Real             refpointy,          /**< reference point for second variable */
   SCIP_Bool             overestimate,       /**< whether to compute an overestimator instead of an underestimator */
   SCIP_Real             xcoef,              /**< x coefficient of linear inequality; must be in {-1,0,1} */
   SCIP_Real             ycoef,              /**< y coefficient of linear inequality */
   SCIP_Real             constant,           /**< constant of linear inequality */
   SCIP_Real* RESTRICT   lincoefx,           /**< buffer to store coefficient of first  variable in linearization */
   SCIP_Real* RESTRICT   lincoefy,           /**< buffer to store coefficient of second variable in linearization */
   SCIP_Real* RESTRICT   linconstant,        /**< buffer to store constant of linearization */
   SCIP_Bool* RESTRICT   success             /**< buffer to store whether linearization was successful */
   )
{
   SCIP_Real xs[2] = {lbx, ubx};
   SCIP_Real ys[2] = {lby, uby};
   SCIP_Real minx;
   SCIP_Real maxx;
   SCIP_Real miny;
   SCIP_Real maxy;
   SCIP_Real tmp;
   SCIP_Real mj;
   SCIP_Real qj;
   SCIP_Real xj;
   SCIP_Real yj;
   SCIP_Real vx;
   SCIP_Real vy;
   int n;
   int i;

   assert(scip != NULL);
   assert(!SCIPisInfinity(scip,  lbx));
   assert(!SCIPisInfinity(scip, -ubx));
   assert(!SCIPisInfinity(scip,  lby));
   assert(!SCIPisInfinity(scip, -uby));
   assert(SCIPisLE(scip, lbx, ubx));
   assert(SCIPisLE(scip, lby, uby));
   assert(SCIPisLE(scip, lbx, refpointx));
   assert(SCIPisGE(scip, ubx, refpointx));
   assert(SCIPisLE(scip, lby, refpointy));
   assert(SCIPisGE(scip, uby, refpointy));
   assert(lincoefx != NULL);
   assert(lincoefy != NULL);
   assert(linconstant != NULL);
   assert(success != NULL);
   assert(xcoef == 0.0 || xcoef == -1.0 || xcoef == 1.0); /*lint !e777*/
   assert(ycoef != SCIP_INVALID && ycoef != 0.0); /*lint !e777*/
   assert(constant != SCIP_INVALID); /*lint !e777*/

   *success = FALSE;
   *lincoefx = SCIP_INVALID;
   *lincoefy = SCIP_INVALID;
   *linconstant = SCIP_INVALID;

   /* reference point does not satisfy linear inequality */
   if( SCIPisFeasGT(scip, xcoef * refpointx - ycoef * refpointy - constant, 0.0) )
      return;

   /* compute minimal and maximal bounds on x and y for accepting the reference point */
   minx = lbx + 0.01 * (ubx-lbx);
   maxx = ubx - 0.01 * (ubx-lbx);
   miny = lby + 0.01 * (uby-lby);
   maxy = uby - 0.01 * (uby-lby);

   /* check whether the reference point is in [minx,maxx]x[miny,maxy] */
   if( SCIPisLE(scip, refpointx, minx) || SCIPisGE(scip, refpointx, maxx)
      || SCIPisLE(scip, refpointy, miny) || SCIPisGE(scip, refpointy, maxy) )
      return;

   /* always consider xy without the bilinear coefficient */
   if( bilincoef < 0.0 )
      overestimate = !overestimate;

   /* we use same notation as in "Convex envelopes of bivariate functions through the solution of KKT systems", 2016 */
   mj = xcoef / ycoef;
   qj = -constant / ycoef;

   /* mj > 0 => underestimate; mj < 0 => overestimate */
   if( SCIPisNegative(scip, mj) != overestimate )
      return;

   /* get the corner point that satisfies the linear inequality xcoef*x <= ycoef*y + constant */
   if( !overestimate )
   {
      ys[0] = uby;
      ys[1] = lby;
   }

   vx = SCIP_INVALID;
   vy = SCIP_INVALID;
   n = 0;
   for( i = 0; i < 2; ++i )
   {
      SCIP_Real activity = xcoef * xs[i] - ycoef * ys[i] - constant;
      if( SCIPisLE(scip, activity, 0.0) )
      {
         /* corner point is satisfies inequality */
         vx = xs[i];
         vy = ys[i];
      }
      else if( SCIPisFeasGT(scip, activity, 0.0) )
         /* corner point is clearly cut off */
         ++n;
   }

   /* skip if no corner point satisfies the inequality or if no corner point is cut off (that is, all corner points satisfy the inequality almost [1e-9..1e-6]) */
   if( n != 1 || vx == SCIP_INVALID || vy == SCIP_INVALID ) /*lint !e777*/
      return;

   tmp = mj*(refpointx - vx) + vy - refpointy;
   if( SCIPisZero(scip, tmp) )
      return;

   /* (xj,yj) is the projection onto the line xcoef*x = ycoef*y + constant */
   xj = (refpointx*(vy - qj) - vx*(refpointy - qj)) / tmp;
   yj = mj * xj + qj;

   assert(SCIPisFeasEQ(scip, xcoef*xj - ycoef*yj - constant, 0.0));

   /* check whether the projection is in [minx,maxx] x [miny,maxy]; this avoids numerical difficulties when the
    * projection is close to the variable bounds
    */
   if( SCIPisLE(scip, xj, minx) || SCIPisGE(scip, xj, maxx) || SCIPisLE(scip, yj, miny) || SCIPisGE(scip, yj, maxy) )
      return;

   assert(vy - mj*vx - qj != 0.0);

   *lincoefy = (mj*SQR(xj) - 2.0*mj*vx*xj - qj*vx + vx*vy) / (vy - mj*vx - qj);
   *lincoefx = 2.0*mj*xj + qj - mj*(*lincoefy);
   *linconstant = -mj*SQR(xj) - (*lincoefy)*qj;

   /* consider the bilinear coefficient */
   *lincoefx *= bilincoef;
   *lincoefy *= bilincoef;
   *linconstant *= bilincoef;

   /* cut needs to be tight at (vx,vy) and (xj,yj); otherwise we consider the cut to be numerically bad */
   *success = SCIPisFeasEQ(scip, (*lincoefx)*vx + (*lincoefy)*vy + (*linconstant), bilincoef*vx*vy)
      && SCIPisFeasEQ(scip, (*lincoefx)*xj + (*lincoefy)*yj + (*linconstant), bilincoef*xj*yj);

#ifndef NDEBUG
   {
      SCIP_Real activity = (*lincoefx)*refpointx + (*lincoefy)*refpointy + (*linconstant);

      /* cut needs to under- or overestimate the bilinear term at the reference point */
      if( bilincoef < 0.0 )
         overestimate = !overestimate;

      if( overestimate )
         assert(SCIPisFeasGE(scip, activity, bilincoef*refpointx*refpointy));
      else
         assert(SCIPisFeasLE(scip, activity, bilincoef*refpointx*refpointy));
   }
#endif
}

/** helper function to compute the convex envelope of a bilinear term when two linear inequalities are given; we
 *  use the same notation and formulas as in Locatelli 2016
 */
static
void computeBilinEnvelope2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             x,                  /**< reference point for x */
   SCIP_Real             y,                  /**< reference point for y */
   SCIP_Real             mi,                 /**< coefficient of x in the first linear inequality */
   SCIP_Real             qi,                 /**< constant in the first linear inequality */
   SCIP_Real             mj,                 /**< coefficient of x in the second linear inequality */
   SCIP_Real             qj,                 /**< constant in the second linear inequality */
   SCIP_Real* RESTRICT   xi,                 /**< buffer to store x coordinate of the first point */
   SCIP_Real* RESTRICT   yi,                 /**< buffer to store y coordinate of the first point */
   SCIP_Real* RESTRICT   xj,                 /**< buffer to store x coordinate of the second point */
   SCIP_Real* RESTRICT   yj,                 /**< buffer to store y coordinate of the second point */
   SCIP_Real* RESTRICT   xcoef,              /**< buffer to store the x coefficient of the envelope */
   SCIP_Real* RESTRICT   ycoef,              /**< buffer to store the y coefficient of the envelope */
   SCIP_Real* RESTRICT   constant            /**< buffer to store the constant of the envelope */
   )
{
   assert(xi != NULL);
   assert(yi != NULL);
   assert(xj != NULL);
   assert(yj != NULL);
   assert(xcoef != NULL);
   assert(ycoef != NULL);
   assert(constant != NULL);

   if( SCIPisEQ(scip, mi, mj) )
   {
      *xi = (x + mi * y - qi) / (2.0*mi);
      *yi = mi*(*xi) + qi;
      *xj = (*xi) + (qi - qj)/ (2.0*mi);
      *yj = mj * (*xj) + qj;
      *ycoef = (*xi) + (qi - qj) / (4.0*mi); /* note that this is wrong in Locatelli 2016 */
      *xcoef = 2.0*mi*(*xi) - mi * (*ycoef) + qi;
      *constant = -mj*SQR(*xj) - (*ycoef) * qj;
   }
   else if( mi > 0.0 )
   {
      assert(mj > 0.0);

      *xi = (y + SQRT(mi*mj)*x - qi) / (REALABS(mi) + SQRT(mi*mj));
      *yi = mi*(*xi) + qi;
      *xj = (y + SQRT(mi*mj)*x - qj) / (REALABS(mj) + SQRT(mi*mj));
      *yj = mj*(*xj) + qj;
      *ycoef = (2.0*mj*(*xj) + qj - 2.0*mi*(*xi) - qi) / (mj - mi);
      *xcoef = 2.0*mj*(*xj) + qj - mj*(*ycoef);
      *constant = -mj*SQR(*xj) - (*ycoef) * qj;
   }
   else
   {
      assert(mi < 0.0 && mj < 0.0);

      /* apply variable transformation x = -x in case for overestimation */
      computeBilinEnvelope2(scip, -x, y, -mi, qi, -mj, qj, xi, yi, xj, yj, xcoef, ycoef, constant);

      /* revert transformation; multiply cut by -1 and change -x by x */
      *xi = -(*xi);
      *xj = -(*xj);
      *ycoef = -(*ycoef);
      *constant = -(*constant);
   }
}

/** computes coefficients of linearization of a bilinear term in a reference point when given two linear inequality
 *  involving only the variables of the bilinear term
 *
 *  @note the formulas are extracted from "Convex envelopes of bivariate functions through the solution of KKT systems"
 *        by Marco Locatelli
 *
 */
void SCIPcomputeBilinEnvelope2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bilincoef,          /**< coefficient of bilinear term */
   SCIP_Real             lbx,                /**< lower bound on first variable */
   SCIP_Real             ubx,                /**< upper bound on first variable */
   SCIP_Real             refpointx,          /**< reference point for first variable */
   SCIP_Real             lby,                /**< lower bound on second variable */
   SCIP_Real             uby,                /**< upper bound on second variable */
   SCIP_Real             refpointy,          /**< reference point for second variable */
   SCIP_Bool             overestimate,       /**< whether to compute an overestimator instead of an underestimator */
   SCIP_Real             xcoef1,             /**< x coefficient of linear inequality; must be in {-1,0,1} */
   SCIP_Real             ycoef1,             /**< y coefficient of linear inequality */
   SCIP_Real             constant1,          /**< constant of linear inequality */
   SCIP_Real             xcoef2,             /**< x coefficient of linear inequality; must be in {-1,0,1} */
   SCIP_Real             ycoef2,             /**< y coefficient of linear inequality */
   SCIP_Real             constant2,          /**< constant of linear inequality */
   SCIP_Real* RESTRICT   lincoefx,           /**< buffer to store coefficient of first  variable in linearization */
   SCIP_Real* RESTRICT   lincoefy,           /**< buffer to store coefficient of second variable in linearization */
   SCIP_Real* RESTRICT   linconstant,        /**< buffer to store constant of linearization */
   SCIP_Bool* RESTRICT   success             /**< buffer to store whether linearization was successful */
   )
{
   SCIP_Real mi, mj, qi, qj, xi, xj, yi, yj;
   SCIP_Real xcoef, ycoef, constant;
   SCIP_Real minx, maxx, miny, maxy;

   assert(scip != NULL);
   assert(!SCIPisInfinity(scip,  lbx));
   assert(!SCIPisInfinity(scip, -ubx));
   assert(!SCIPisInfinity(scip,  lby));
   assert(!SCIPisInfinity(scip, -uby));
   assert(SCIPisLE(scip, lbx, ubx));
   assert(SCIPisLE(scip, lby, uby));
   assert(SCIPisLE(scip, lbx, refpointx));
   assert(SCIPisGE(scip, ubx, refpointx));
   assert(SCIPisLE(scip, lby, refpointy));
   assert(SCIPisGE(scip, uby, refpointy));
   assert(lincoefx != NULL);
   assert(lincoefy != NULL);
   assert(linconstant != NULL);
   assert(success != NULL);
   assert(xcoef1 != 0.0 && xcoef1 != SCIP_INVALID); /*lint !e777*/
   assert(ycoef1 != SCIP_INVALID && ycoef1 != 0.0); /*lint !e777*/
   assert(constant1 != SCIP_INVALID); /*lint !e777*/
   assert(xcoef2 != 0.0 && xcoef2 != SCIP_INVALID); /*lint !e777*/
   assert(ycoef2 != SCIP_INVALID && ycoef2 != 0.0); /*lint !e777*/
   assert(constant2 != SCIP_INVALID); /*lint !e777*/

   *success = FALSE;
   *lincoefx = SCIP_INVALID;
   *lincoefy = SCIP_INVALID;
   *linconstant = SCIP_INVALID;

   /* reference point does not satisfy linear inequalities */
   if( SCIPisFeasGT(scip, xcoef1 * refpointx - ycoef1 * refpointy - constant1, 0.0)
      || SCIPisFeasGT(scip, xcoef2 * refpointx - ycoef2 * refpointy - constant2, 0.0) )
      return;

   /* compute minimal and maximal bounds on x and y for accepting the reference point */
   minx = lbx + 0.01 * (ubx-lbx);
   maxx = ubx - 0.01 * (ubx-lbx);
   miny = lby + 0.01 * (uby-lby);
   maxy = uby - 0.01 * (uby-lby);

   /* check the reference point is in the interior of the domain */
   if( SCIPisLE(scip, refpointx, minx) || SCIPisGE(scip, refpointx, maxx)
      || SCIPisLE(scip, refpointy, miny) || SCIPisFeasGE(scip, refpointy, maxy) )
      return;

   /* the sign of the x-coefficients of the two inequalities must be different; otherwise the convex or concave
    * envelope can be computed via SCIPcomputeBilinEnvelope1 for each inequality separately
    */
   if( (xcoef1 > 0) == (xcoef2 > 0) )
      return;

   /* always consider xy without the bilinear coefficient */
   if( bilincoef < 0.0 )
      overestimate = !overestimate;

   /* we use same notation as in "Convex envelopes of bivariate functions through the solution of KKT systems", 2016 */
   mi = xcoef1 / ycoef1;
   qi = -constant1 / ycoef1;
   mj = xcoef2 / ycoef2;
   qj = -constant2 / ycoef2;

   /* mi, mj > 0 => underestimate; mi, mj < 0 => overestimate */
   if( SCIPisNegative(scip, mi) != overestimate || SCIPisNegative(scip, mj) != overestimate )
      return;

   /* compute cut according to Locatelli 2016 */
   computeBilinEnvelope2(scip, refpointx, refpointy, mi, qi, mj, qj, &xi, &yi, &xj, &yj, &xcoef, &ycoef, &constant);
   assert(SCIPisEQ(scip, mi*xi + qi, yi));
   assert(SCIPisEQ(scip, mj*xj + qj, yj));

   /* it might happen that (xi,yi) = (xj,yj) if the two lines intersect */
   if( SCIPisEQ(scip, xi, xj) && SCIPisEQ(scip, yi, yj) )
      return;

   /* check whether projected points are in the interior */
   if( SCIPisLE(scip, xi, minx) || SCIPisGE(scip, xi, maxx) || SCIPisLE(scip, yi, miny) || SCIPisGE(scip, yi, maxy) )
      return;
   if( SCIPisLE(scip, xj, minx) || SCIPisGE(scip, xj, maxx) || SCIPisLE(scip, yj, miny) || SCIPisGE(scip, yj, maxy) )
      return;

   *lincoefx = bilincoef * xcoef;
   *lincoefy = bilincoef * ycoef;
   *linconstant = bilincoef * constant;

   /* cut needs to be tight at (vx,vy) and (xj,yj) */
   *success = SCIPisFeasEQ(scip, (*lincoefx)*xi + (*lincoefy)*yi + (*linconstant), bilincoef*xi*yi)
      && SCIPisFeasEQ(scip, (*lincoefx)*xj + (*lincoefy)*yj + (*linconstant), bilincoef*xj*yj);

#ifndef NDEBUG
   {
      SCIP_Real activity = (*lincoefx)*refpointx + (*lincoefy)*refpointy + (*linconstant);

      /* cut needs to under- or overestimate the bilinear term at the reference point */
      if( bilincoef < 0.0 )
         overestimate = !overestimate;

      if( overestimate )
         assert(SCIPisFeasGE(scip, activity, bilincoef*refpointx*refpointy));
      else
         assert(SCIPisFeasLE(scip, activity, bilincoef*refpointx*refpointy));
   }
#endif
}

/** creates an NLP relaxation and stores it in a given NLPI problem; the function computes for each variable which the
 *  number of non-linearly occurrences and stores it in the nlscore array
 *
 *  @note the first row corresponds always to the cutoff row (even if cutoffbound is SCIPinfinity(scip))
 **/
SCIP_RETCODE SCIPcreateNlpiProb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLROW**          nlrows,             /**< nonlinear rows */
   int                   nnlrows,            /**< total number of nonlinear rows */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< empty nlpi problem */
   SCIP_HASHMAP*         var2idx,            /**< empty hash map to store mapping between variables and indices in nlpi
                                              *   problem */
   SCIP_Real*            nlscore,            /**< array to store the score of each nonlinear variable (NULL if not
                                              *   needed) */
   SCIP_Real             cutoffbound,        /**< cutoff bound */
   SCIP_Bool             setobj,             /**< should the objective function be set? */
   SCIP_Bool             onlyconvex          /**< filter only for convex constraints */
   )
{
   SCIP_EXPRTREE** exprtrees;
   int** exprvaridxs;
   SCIP_QUADELEM** quadelems;
   int* nquadelems;
   SCIP_Real** linvals;
   int** lininds;
   int* nlininds;
   SCIP_Real* lhss;
   SCIP_Real* rhss;
   const char** names;
   SCIP_VAR** vars;
   int nvars;
   SCIP_Real* lbs;
   SCIP_Real* ubs;
   SCIP_Real* objvals = NULL;
   int* objinds = NULL;
   const char** varnames;
   int nobjinds;
   int nconss;
   int i;

   assert(nlpiprob != NULL);
   assert(var2idx != NULL);
   assert(nlrows != NULL);
   assert(nnlrows > 0);
   assert(nlpi != NULL);

   SCIPdebugMsg(scip, "call SCIPcreateConvexNlpNlobbt() with cutoffbound %g\n", cutoffbound);

   if( nlscore != NULL )
   {
      BMSclearMemoryArray(nlscore, SCIPgetNVars(scip));
   }
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   nconss = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &exprtrees, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &exprvaridxs, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nquadelems, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lininds, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlininds, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &names, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lhss, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhss, nnlrows + 1) );

   if( setobj )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &objvals, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &objinds, nvars) );
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &lbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varnames, nvars) );

   /* create a unique mapping between variables and {0,..,nvars-1} */
   nobjinds = 0;
   for( i = 0; i < nvars; ++i )
   {
      assert(vars[i] != NULL);
      SCIP_CALL( SCIPhashmapInsert(var2idx, (void*)vars[i], (void*)(size_t)i) );

      lbs[i] = SCIPvarGetLbLocal(vars[i]);
      ubs[i] = SCIPvarGetUbLocal(vars[i]);
      varnames[i] = SCIPvarGetName(vars[i]);

      /* collect non-zero objective coefficients */
      if( setobj && !SCIPisZero(scip, SCIPvarGetObj(vars[i])) )
      {
         assert(objvals != NULL);
         assert(objinds != NULL);

         objvals[nobjinds] = SCIPvarGetObj(vars[i]);
         objinds[nobjinds] = i;
         ++nobjinds;
      }
   }

   /* add variables */
   SCIP_CALL( SCIPnlpiAddVars(nlpi, nlpiprob, nvars, lbs, ubs, varnames) );
   SCIPfreeBufferArray(scip, &varnames);
   SCIPfreeBufferArray(scip, &ubs);
   SCIPfreeBufferArray(scip, &lbs);

   /* set the objective function */
   if( setobj )
   {
      if( nobjinds > 0 )
      {
         SCIP_CALL( SCIPnlpiSetObjective(nlpi, nlpiprob, nobjinds, objinds, objvals, 0, NULL, NULL, NULL, 0.0) );
      }

      SCIPfreeBufferArray(scip, &objinds);
      SCIPfreeBufferArray(scip, &objvals);
   }

   /* add row for cutoff bound even if cutoffbound == SCIPinfinity() */
   lhss[nconss] = -SCIPinfinity(scip);
   rhss[nconss] = cutoffbound;
   names[nconss] = "objcutoff";
   lininds[nconss] = NULL;
   linvals[nconss] = NULL;
   nlininds[nconss] = 0;
   nquadelems[nconss] = 0;
   quadelems[nconss] = NULL;
   exprtrees[nconss] = NULL;
   exprvaridxs[nconss] = NULL;

   SCIP_CALL( SCIPallocBufferArray(scip, &lininds[nconss], nvars) ); /*lint !e866*/
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals[nconss], nvars) ); /*lint !e866*/

   for( i = 0; i < nvars; ++i )
   {
      if( !SCIPisZero(scip, SCIPvarGetObj(vars[i])) )
      {
         linvals[nconss][nlininds[nconss]] = SCIPvarGetObj(vars[i]);
         lininds[nconss][nlininds[nconss]] = i;
         ++nlininds[nconss];
      }
   }
   ++nconss;

   /* add convex nonlinear rows to NLPI problem */
   for( i = 0; i < nnlrows; ++i )
   {
      SCIP_Bool userhs;
      SCIP_Bool uselhs;
      int k;

      assert(nlrows[i] != NULL);

      uselhs = FALSE;
      userhs = FALSE;

      /* check curvature together with constraint sides of a nonlinear row */
      if( SCIPnlrowGetNQuadElems(nlrows[i]) == 0 && SCIPnlrowGetExprtree(nlrows[i]) == NULL )
      {
         uselhs = TRUE;
         userhs = TRUE;
      }
      else
      {
         if( (!onlyconvex || SCIPnlrowGetCurvature(nlrows[i]) == SCIP_EXPRCURV_CONVEX)
            && !SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrows[i])) )
            userhs = TRUE;
         if( (!onlyconvex || SCIPnlrowGetCurvature(nlrows[i]) == SCIP_EXPRCURV_CONCAVE)
            && !SCIPisInfinity(scip, SCIPnlrowGetLhs(nlrows[i])) )
            uselhs = TRUE;
      }

      if( !uselhs && !userhs )
         continue;

      lhss[nconss] = uselhs ? SCIPnlrowGetLhs(nlrows[i]) - SCIPnlrowGetConstant(nlrows[i]) : -SCIPinfinity(scip);
      rhss[nconss] = userhs ? SCIPnlrowGetRhs(nlrows[i]) - SCIPnlrowGetConstant(nlrows[i]) :  SCIPinfinity(scip);
      names[nconss] = SCIPnlrowGetName(nlrows[i]);
      nlininds[nconss] = 0;
      lininds[nconss] = NULL;
      linvals[nconss] = NULL;
      nquadelems[nconss] = 0;
      quadelems[nconss] = NULL;
      exprtrees[nconss] = NULL;
      exprvaridxs[nconss] = NULL;

      /* copy linear part */
      if( SCIPnlrowGetNLinearVars(nlrows[i]) > 0 )
      {
         SCIP_VAR* var;

         nlininds[nconss] = SCIPnlrowGetNLinearVars(nlrows[i]);

         SCIP_CALL( SCIPallocBufferArray(scip, &lininds[nconss], nlininds[nconss]) ); /*lint !e866*/
         SCIP_CALL( SCIPallocBufferArray(scip, &linvals[nconss], nlininds[nconss]) ); /*lint !e866*/

         for( k = 0; k < nlininds[nconss]; ++k )
         {
            var = SCIPnlrowGetLinearVars(nlrows[i])[k];
            assert(var != NULL);
            assert(SCIPhashmapExists(var2idx, (void*)var));

            lininds[nconss][k] = (int)(size_t)SCIPhashmapGetImage(var2idx, (void*)var);
            assert(var == vars[lininds[nconss][k]]);
            linvals[nconss][k] = SCIPnlrowGetLinearCoefs(nlrows[i])[k];
         }
      }

      /* copy quadratic part */
      if( SCIPnlrowGetNQuadElems(nlrows[i]) > 0 )
      {
         SCIP_QUADELEM quadelem;
         SCIP_VAR* var1;
         SCIP_VAR* var2;

         nquadelems[nconss] = SCIPnlrowGetNQuadElems(nlrows[i]);
         SCIP_CALL( SCIPallocBufferArray(scip, &quadelems[nconss], nquadelems[nconss]) ); /*lint !e866*/

         for( k = 0; k < nquadelems[nconss]; ++k )
         {
            quadelem = SCIPnlrowGetQuadElems(nlrows[i])[k];

            var1 = SCIPnlrowGetQuadVars(nlrows[i])[quadelem.idx1];
            assert(var1 != NULL);
            assert(SCIPhashmapExists(var2idx, (void*)var1));

            var2 = SCIPnlrowGetQuadVars(nlrows[i])[quadelem.idx2];
            assert(var2 != NULL);
            assert(SCIPhashmapExists(var2idx, (void*)var2));

            quadelems[nconss][k].coef = quadelem.coef;
            quadelems[nconss][k].idx1 = (int)(size_t)SCIPhashmapGetImage(var2idx, (void*)var1);
            quadelems[nconss][k].idx2 = (int)(size_t)SCIPhashmapGetImage(var2idx, (void*)var2);

            /* expr.c assumes that the indices are ordered */
            if( quadelems[nconss][k].idx1 > quadelems[nconss][k].idx2 )
            {
               SCIPswapInts(&quadelems[nconss][k].idx1, &quadelems[nconss][k].idx2);
            }
            assert(quadelems[nconss][k].idx1 <= quadelems[nconss][k].idx2);

            /* update nlscore */
            if( nlscore != NULL )
            {
               ++nlscore[quadelems[nconss][k].idx1];
               if( quadelems[nconss][k].idx1 != quadelems[nconss][k].idx2 )
                  ++nlscore[quadelems[nconss][k].idx2];
            }
         }
      }

      /* copy expression tree */
      if( SCIPnlrowGetExprtree(nlrows[i]) != NULL )
      {
         SCIP_VAR* var;

         /* note that we don't need to copy the expression tree here since only the mapping between variables in the
          * tree and the corresponding indices change; this mapping is stored in the exprvaridxs array
          */
         exprtrees[nconss] = SCIPnlrowGetExprtree(nlrows[i]);

         SCIP_CALL( SCIPallocBufferArray(scip, &exprvaridxs[nconss], SCIPexprtreeGetNVars(exprtrees[nconss])) ); /*lint !e866*/

         for( k = 0; k < SCIPexprtreeGetNVars(exprtrees[nconss]); ++k )
         {
            var = SCIPexprtreeGetVars(exprtrees[nconss])[k];
            assert(var != NULL);
            assert(SCIPhashmapExists(var2idx, (void*)var));

            exprvaridxs[nconss][k] = (int)(size_t)SCIPhashmapGetImage(var2idx, (void*)var);

            /* update nlscore */
            if( nlscore != NULL )
               ++nlscore[exprvaridxs[nconss][k]];
         }
      }

      ++nconss;
   }
   assert(nconss > 0);

   /* pass all constraint information to nlpi */
   SCIP_CALL( SCIPnlpiAddConstraints(nlpi, nlpiprob, nconss, lhss, rhss, nlininds, lininds, linvals, nquadelems,
         quadelems, exprvaridxs, exprtrees, names) );

   /* free memory */
   for( i = nconss - 1; i > 0; --i )
   {
      if( exprtrees[i] != NULL )
      {
         assert(exprvaridxs[i] != NULL);
         SCIPfreeBufferArray(scip, &exprvaridxs[i]);
      }

      if( nquadelems[i] > 0 )
      {
         assert(quadelems[i] != NULL);
         SCIPfreeBufferArray(scip, &quadelems[i]);
      }

      if( nlininds[i] > 0 )
      {
         assert(linvals[i] != NULL);
         assert(lininds[i] != NULL);
         SCIPfreeBufferArray(scip, &linvals[i]);
         SCIPfreeBufferArray(scip, &lininds[i]);
      }
   }
   /* free row for cutoff bound even if objective is 0 */
   SCIPfreeBufferArray(scip, &linvals[i]);
   SCIPfreeBufferArray(scip, &lininds[i]);

   SCIPfreeBufferArray(scip, &rhss);
   SCIPfreeBufferArray(scip, &lhss);
   SCIPfreeBufferArray(scip, &names);
   SCIPfreeBufferArray(scip, &nlininds);
   SCIPfreeBufferArray(scip, &lininds);
   SCIPfreeBufferArray(scip, &linvals);
   SCIPfreeBufferArray(scip, &nquadelems);
   SCIPfreeBufferArray(scip, &quadelems);
   SCIPfreeBufferArray(scip, &exprvaridxs);
   SCIPfreeBufferArray(scip, &exprtrees);

   return SCIP_OKAY;
}

/** updates bounds of each variable and the cutoff row in the nlpiproblem */
SCIP_RETCODE SCIPupdateNlpiProb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< nlpi problem representing the convex NLP relaxation */
   SCIP_HASHMAP*         var2nlpiidx,        /**< mapping between variables and nlpi indices */
   SCIP_VAR**            nlpivars,           /**< array containing all variables of the nlpi */
   int                   nlpinvars,          /**< total number of nlpi variables */
   SCIP_Real             cutoffbound         /**< new cutoff bound */
   )
{
   SCIP_Real* lbs;
   SCIP_Real* ubs;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int* inds;
   int i;

   SCIPdebugMsg(scip, "call SCIPupdateConvexNlpNlobbt()\n");

   /* update variable bounds */
   SCIP_CALL( SCIPallocBufferArray(scip, &lbs, nlpinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubs, nlpinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, nlpinvars) );

   for( i = 0; i < nlpinvars; ++i )
   {
      assert(nlpivars[i] != NULL);
      assert(SCIPhashmapExists(var2nlpiidx, (void*)nlpivars[i]));

      lbs[i] = SCIPvarGetLbLocal(nlpivars[i]);
      ubs[i] = SCIPvarGetUbLocal(nlpivars[i]);
      inds[i] = (int)(uintptr_t)SCIPhashmapGetImage(var2nlpiidx, (void*)nlpivars[i]);
      assert(inds[i] >= 0 && inds[i] < nlpinvars);
   }

   SCIP_CALL( SCIPnlpiChgVarBounds(nlpi, nlpiprob, nlpinvars, inds, lbs, ubs) );

   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &ubs);
   SCIPfreeBufferArray(scip, &lbs);

   /* update cutoff row */
   lhs = -SCIPinfinity(scip);
   rhs = cutoffbound;
   i = 0;

   SCIP_CALL( SCIPnlpiChgConsSides(nlpi, nlpiprob, 1, &i, &lhs, &rhs) );

   return SCIP_OKAY;
}

/** adds linear rows to the NLP relaxation */
SCIP_RETCODE SCIPaddNlpiProbRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< nlpi problem */
   SCIP_HASHMAP*         var2idx,            /**< empty hash map to store mapping between variables and indices in nlpi
                                              *   problem */
   SCIP_ROW**            rows,               /**< rows to add */
   int                   nrows               /**< total number of rows to add */
   )
{
   const char** names;
   SCIP_Real* lhss;
   SCIP_Real* rhss;
   SCIP_Real** linvals;
   int** lininds;
   int* nlininds;
   int i;

   assert(nlpi != NULL);
   assert(nlpiprob != NULL);
   assert(var2idx != NULL);
   assert(nrows == 0 || rows != NULL);

   SCIPdebugMsg(scip, "call SCIPaddConvexNlpRowsNlobbt() with %d rows\n", nrows);

   if( nrows <= 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &names, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lhss, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhss, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lininds, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlininds, nrows) );

   for( i = 0; i < nrows; ++i )
   {
      int k;

      assert(rows[i] != NULL);
      assert(SCIProwGetNNonz(rows[i]) <= SCIPgetNVars(scip));

      names[i] = SCIProwGetName(rows[i]);
      lhss[i] = SCIProwGetLhs(rows[i]) - SCIProwGetConstant(rows[i]);
      rhss[i] = SCIProwGetRhs(rows[i]) - SCIProwGetConstant(rows[i]);
      nlininds[i] = SCIProwGetNNonz(rows[i]);
      linvals[i] = SCIProwGetVals(rows[i]);
      lininds[i] = NULL;

      SCIP_CALL( SCIPallocBufferArray(scip, &lininds[i], SCIProwGetNNonz(rows[i])) ); /*lint !e866*/

      for( k = 0; k < SCIProwGetNNonz(rows[i]); ++k )
      {
         SCIP_VAR* var;

         var = SCIPcolGetVar(SCIProwGetCols(rows[i])[k]);
         assert(var != NULL);
         assert(SCIPhashmapExists(var2idx, (void*)var));

         lininds[i][k] = (int)(size_t)SCIPhashmapGetImage(var2idx, (void*)var);
         assert(lininds[i][k] >= 0 && lininds[i][k] < SCIPgetNVars(scip));
      }
   }

   /* pass all linear rows to the nlpi */
   SCIP_CALL( SCIPnlpiAddConstraints(nlpi, nlpiprob, nrows, lhss, rhss, nlininds, lininds, linvals, NULL,
         NULL, NULL, NULL, names) );

   /* free memory */
   for( i = nrows - 1; i >= 0; --i )
   {
      SCIPfreeBufferArray(scip, &lininds[i]);
   }
   SCIPfreeBufferArray(scip, &nlininds);
   SCIPfreeBufferArray(scip, &lininds);
   SCIPfreeBufferArray(scip, &linvals);
   SCIPfreeBufferArray(scip, &rhss);
   SCIPfreeBufferArray(scip, &lhss);
   SCIPfreeBufferArray(scip, &names);

   return SCIP_OKAY;
}

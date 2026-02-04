/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_certificate.c
 * @brief  public methods for certified solving
 * @author Leon Eifler
 * @author Ambros Gleixner
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <ctype.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>
#ifndef _WIN32
#include <strings.h> /*lint --e{766}*/
#endif


#include "lpi/lpi.h"
#include "scip/exprinterpret.h"
#include "scip/nlpi.h"
#include "scip/benders.h"
#include "scip/benderscut.h"
#include "scip/branch.h"
#include "scip/branch_nodereopt.h"
#include "scip/certificate.h"
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
#include "scip/lpexact_bounding.h"
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
#include "scip/sepastoreexact.h"
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

#include "scip/scip_certificate.h"
#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nlp.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_var.h"

#include "scip/pub_cons.h"
#include "scip/pub_fileio.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_sol.h"
#include "scip/pub_var.h"
#include "scip/pub_lpexact.h"
#include "scip/struct_certificate.h"
#include "scip/struct_lpexact.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif

/** returns the sense of an inequality */
static
char getInequalitySense(
   SCIP_Bool             isgreaterthan       /**< is the inequality a greater than inequality? */
   )
{
   return isgreaterthan ? 'G' : 'L';
}

/** returns whether certificate output is activated
 *
 *  @todo add a flag set->certificate_enabled to store the return value of this method for easier and faster access
 */
SCIP_Bool SCIPisCertified(
   SCIP*                 scip                /**< certificate information */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(scip->stat != NULL);

   if( !(scip->set->exact_enable) )
      return FALSE;
   else if( scip->set->stage == SCIP_STAGE_SOLVING )
      return SCIPcertificateIsEnabled(scip->stat->certificate);
   else
      return !(scip->set->certificate_filename[0] == '-' && scip->set->certificate_filename[1] == '\0');
}

/** should the certificate track bound changes?
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Bool SCIPshouldCertificateTrackBounds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   return SCIPisCertified(scip) && scip->set->stage >= SCIP_STAGE_INITSOLVE && !SCIPinProbing(scip);
}

/** prints constraint to certificate
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcertifyCons(
   SCIP*                 scip,               /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the original solution be printed or in transformed space */
   const char*           consname,           /**< name of the constraint */
   const char            sense,              /**< sense of the constraint, i.e., G, L, or E */
   SCIP_RATIONAL*        side,               /**< left/right-hand side */
   int                   len,                /**< number of nonzeros */
   int*                  ind,                /**< index array */
   SCIP_RATIONAL**       val                 /**< coefficient array */
   )
{
   assert(scip != NULL);
   assert(scip->stat != NULL);
   assert(scip->stat->certificate != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcertifyCons", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcertificatePrintCons(scip->stat->certificate, isorigfile, consname, sense, side, len, ind, val) );

   return SCIP_OKAY;
}

/** prints activity bound to proof section
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcertifyActivityVarBound(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           linename,           /**< name of the unsplitting line */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound (upper/lower) */
   SCIP_Real             newbound,           /**< pointer to lower bound on the objective, NULL indicating infeasibility */
   SCIP_Bool             ismaxactivity,      /**< TRUE for maxactivity, FALSE for minactivity */
   SCIP_CONS*            constraint,         /**< the constraint */
   SCIP_VAR*             variable,           /**< the variable */
   SCIP_ROWEXACT*        row,                /**< the  corresponding row, or NULL if constraint has no row representation */
   SCIP_RATIONAL**       vals,               /**< value array */
   SCIP_RATIONAL*        lhs,                /**< lhs of the constraint */
   SCIP_RATIONAL*        rhs,                /**< rhs of the constraint */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars               /**< number of values */
   )
{
   /* It would be more efficient if we could do this all in fp artihmetic. However, this is not trivial because the
    * translations between aggregate variables need to be done exactly.
    */
   SCIP_RATIONAL* newboundex;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcertifyActivityVarBound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &newboundex) );
   SCIPrationalSetReal(newboundex, newbound);
   (void) SCIPcertifyActivityVarBoundExact(scip, linename, boundtype,
      newboundex, ismaxactivity, constraint, variable, row, vals, lhs, rhs, vars, nvars);
   SCIPrationalFreeBuffer(SCIPbuffer(scip), &newboundex);

   return SCIP_OKAY;
}

/** prints activity bound to proof section
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcertifyActivityVarBoundExact(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           linename,           /**< name of the unsplitting line */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound (upper/lower) */
   SCIP_RATIONAL*        newbound,           /**< pointer to lower bound on the objective, NULL indicating infeasibility */
   SCIP_Bool             ismaxactivity,      /**< TRUE for maxactivity, FALSE for minactivity */
   SCIP_CONS*            constraint,         /**< the constraint */
   SCIP_VAR*             variable,           /**< the variable */
   SCIP_ROWEXACT*        row,                /**< the  corresponding row, or NULL if constraint has no row representation */
   SCIP_RATIONAL**       vals,               /**< value array */
   SCIP_RATIONAL*        lhs,                /**< lhs of the constraint */
   SCIP_RATIONAL*        rhs,                /**< rhs of the constraint */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars               /**< number of values */
   )
{
   SCIP_CERTIFICATE* certificate;
   SCIP_Longint res;
   SCIP_RETCODE ret = SCIP_OKAY;
   SCIP_RATIONAL* val;
   SCIP_Bool upperboundcontribution;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcertifyActivityVarBoundExact", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert(scip != NULL);
   assert(scip->stat != NULL);
   assert(scip->stat != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(constraint)), "exactlinear") == 0);

   certificate = scip->stat->certificate;

   switch( variable->varstatus )
   {
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_ORIGINAL:
         SCIPABORT();
         return SCIP_ERROR;
      case SCIP_VARSTATUS_NEGATED:

         SCIPrationalMultReal(newbound, newbound, -1.0);
         assert( SCIPvarGetNegationConstant(variable) == 1 );
         SCIPrationalAddReal(newbound, newbound, 1.0);
         ret = SCIPcertifyActivityVarBoundExact(scip, linename,
               boundtype == SCIP_BOUNDTYPE_UPPER ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER,
               newbound, ismaxactivity, constraint, variable->negatedvar, row, vals, lhs, rhs, vars, nvars);
         SCIPrationalAddReal(newbound, newbound, -1.0);
         SCIPrationalMultReal(newbound, newbound, -1.0);
         return ret;
         break;
      case SCIP_VARSTATUS_AGGREGATED:
         SCIPrationalAddProdReal(newbound, variable->exactdata->aggregate.constant, -1.0);
         SCIPrationalDiv(newbound, newbound, variable->exactdata->aggregate.scalar);
         ret = SCIPcertifyActivityVarBoundExact(scip, linename,
               (boundtype == SCIP_BOUNDTYPE_UPPER) == SCIPrationalIsPositive(variable->exactdata->aggregate.scalar) ? SCIP_BOUNDTYPE_UPPER: SCIP_BOUNDTYPE_LOWER,
               newbound,  ismaxactivity, constraint, variable->data.aggregate.var, row, vals, lhs, rhs, vars, nvars);
         SCIPrationalMult(newbound, newbound, variable->exactdata->aggregate.scalar);
         return ret;
         break;
      case SCIP_VARSTATUS_COLUMN:
         break;
      default:
         SCIPABORT();
         return SCIP_ERROR;
   }

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return ret;

   certificate->indexcounter++;

   if( linename == NULL )
   {
      SCIPcertificatePrintProofMessage(certificate, "ACT_L%d ", certificate->indexcounter - 1);
   }
   else
   {
      SCIPcertificatePrintProofMessage(certificate, "%s ", linename);
   }

   /* find the correct value in the constraint */
   val = NULL;
   for( int i = 0; i < nvars; i++ )
   {
      if( vars[i] == variable )
      {
         val = vals[i];
         break;
      }
   }

   assert(val != NULL);

   /* Do we need an upper bound on the contribution val[i]*x_i? (otherwise a lowerbound) */
   upperboundcontribution = (boundtype == SCIP_BOUNDTYPE_UPPER) == SCIPrationalIsPositive(val);
   SCIPcertificatePrintProofMessage(certificate, "%c ", getInequalitySense(upperboundcontribution));

   /* new bound = -newbound * val for now, we print a second line where we scale with 1/val */
   SCIPrationalMult(newbound, newbound, val);
   SCIPrationalNegate(newbound, newbound);
   SCIP_CALL_ABORT( SCIPcertificatePrintProofRational(certificate, newbound) );

   /* print coeffictent of variable -> val */
   SCIPrationalNegate(newbound, newbound);
   SCIPrationalDiv(newbound, newbound, val);
   SCIPcertificatePrintProofMessage(certificate, " 1 %d ", SCIPvarGetCertificateIndex(variable));

   /* negate val, print it and reset it again */
   SCIPrationalNegate(val, val);
   SCIP_CALL_ABORT( SCIPcertificatePrintProofRational(certificate, val) );
   SCIPrationalNegate(val, val);

   if(row != NULL)
      res = SCIPcertificateGetRowIndex(certificate, row, !ismaxactivity);
   else
      res = SCIPcertificateGetConsIndex(certificate, constraint, lhs, rhs, !ismaxactivity);

   SCIPcertificatePrintProofMessage(certificate, " { lin %d %d -1", nvars, res);

   /* print all other variables with their correct bounds */
   for( int i = 0; i < nvars; i++ )
   {
      SCIP_VAR* ivar;
      bool is_upper_bound;
      SCIP_Longint certificateindex;

      assert(!SCIPrationalIsAbsInfinity(vals[i]));

      ivar = vars[i];
      if( ivar == variable )
         continue;

      is_upper_bound = upperboundcontribution != SCIPrationalIsPositive(vals[i]);

      assert(upperboundcontribution != ismaxactivity);

      certificateindex = is_upper_bound ? SCIPvarGetUbCertificateIndexLocal(ivar) :  SCIPvarGetLbCertificateIndexLocal(ivar);
      SCIPcertificatePrintProofMessage(certificate, " %d ", certificateindex);
      SCIP_CALL( SCIPcertificatePrintProofRational(certificate, vals[i]) );
   }
   SCIPcertificatePrintProofMessage(certificate, " } -1\n");

   /* now scale with 1/val */
   certificate->indexcounter++;
   SCIPcertificatePrintProofMessage(certificate, "ACT_L%d %c ", certificate->indexcounter - 1, getInequalitySense(boundtype == SCIP_BOUNDTYPE_LOWER));
   SCIP_CALL( SCIPcertificatePrintProofRational(certificate, newbound) );
   SCIPcertificatePrintProofMessage(certificate, " 1 %d 1 { lin 1 %d ", SCIPvarGetCertificateIndex(variable), certificate->indexcounter - 2);
   SCIPrationalInvert(val, val);
   SCIPrationalNegate(val, val);
   SCIP_CALL( SCIPcertificatePrintProofRational(certificate, val) );

   /* Return val to its original state: */
   SCIPrationalNegate(val, val);
   SCIPrationalInvert(val, val);
   SCIPcertificatePrintProofMessage(certificate, " } -1\n", SCIPvarGetCertificateIndex(variable), certificate->indexcounter - 2);

   /* if variable is integer, round the new bound */
   if( !SCIPrationalIsAbsInfinity(newbound) && SCIPvarGetType(variable) != SCIP_VARTYPE_CONTINUOUS && !SCIPrationalIsIntegral(newbound) )
   {
      certificate->indexcounter++;

      SCIPcertificatePrintProofMessage(certificate, "ACT_R%d %c ", certificate->indexcounter - 1, getInequalitySense(boundtype == SCIP_BOUNDTYPE_LOWER));
      SCIPrationalRoundInteger(newbound, newbound, boundtype == SCIP_BOUNDTYPE_UPPER ? SCIP_R_ROUND_DOWNWARDS : SCIP_R_ROUND_UPWARDS);

      SCIP_CALL( SCIPcertificatePrintProofRational(certificate, newbound) );

      SCIPcertificatePrintProofMessage(certificate, " 1 %d 1", SCIPvarGetCertificateIndex(variable));
      SCIPcertificatePrintProofMessage(certificate, " { rnd 1 %d 1 } -1\n", certificate->indexcounter - 2);
   }

#ifndef NDEBUG
   certificate->lastinfo->isbound = TRUE;
   certificate->lastinfo->boundtype = boundtype;
   certificate->lastinfo->varindex = SCIPvarGetCertificateIndex(variable);
   certificate->lastinfo->isglobal = FALSE;
   certificate->lastinfo->certificateindex = certificate->indexcounter - 1;
   SCIPrationalSetRational(certificate->lastinfo->boundval, newbound);
#endif
   (void) SCIPcertificateSetLastBoundIndex(certificate, certificate->indexcounter - 1);

   return SCIP_OKAY;
}

/** prints activity conflict to certificate file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcertifyActivityConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_ROWEXACT*        row,                /**< corresponding row, or NULL if constraint does not have representation as row */
   SCIP_RATIONAL*        lhs,                /**< lhs of the constraint */
   SCIP_RATIONAL*        rhs,                /**< rhs of the constraint */
   int                   nvals,              /**< number of values */
   SCIP_RATIONAL**       vals,               /**< value array */
   SCIP_VAR**            vars,               /**< variable array */
   SCIP_RATIONAL*        diff,               /**< difference between min/max activity as lhs/rhs */
   SCIP_Bool             userhs              /**< is rhs or lhs used */
   )
{
   SCIP_CERTIFICATE* certificate;
   SCIP_Longint conscertificateindex;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcertifyActivityConflict", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPisCertified(scip) )
      return SCIP_OKAY;

   certificate = SCIPgetCertificate(scip);

   if( row != NULL )
      conscertificateindex = SCIPcertificateGetRowIndex(SCIPgetCertificate(scip), row, userhs);
   else
      conscertificateindex = SCIPcertificateGetConsIndex(certificate, cons, lhs, rhs, userhs);

   assert(conscertificateindex != LONG_MAX);

   SCIPcertificatePrintProofMessage(certificate, "ActivityConflict%d ", certificate->indexcounter);
   SCIPcertificatePrintProofMessage(certificate, userhs ? "G " : "L ");

   SCIP_CALL( SCIPcertificatePrintProofRational(certificate, diff) );
   SCIPcertificatePrintProofMessage(certificate, " 0 { lin %d %d -1", nvals + 1, conscertificateindex);
   for( int i = 0; i < nvals; i++ )
   {
      SCIP_VAR* var;
      bool is_upper_bound;
      SCIP_Longint certificateindex;
      var = row == NULL ? vars[i] : row->cols[i]->var;
      is_upper_bound = userhs != SCIPrationalIsPositive(vals[i]);
      certificateindex = is_upper_bound ? SCIPvarGetUbCertificateIndexLocal(var) :  SCIPvarGetLbCertificateIndexLocal(var);
      SCIPcertificatePrintProofMessage(certificate, " %d ", certificateindex);
      SCIP_CALL( SCIPcertificatePrintProofRational(certificate, vals[i]) );
   }
   SCIPcertificatePrintProofMessage(certificate, " } -1\n");

   SCIP_CALL( SCIPcertificateUpdateParentData(certificate, SCIPgetCurrentNode(scip), certificate->indexcounter, NULL) );
   certificate->indexcounter++;
   return SCIP_OKAY;
}

/** adds aggregation information to certificate for one row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPaddCertificateAggrInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< agrrrow that results from the aggregation */
   SCIP_ROW**            aggrrows,           /**< array of rows used fo the aggregation */
   SCIP_Real*            weights,            /**< array of weights */
   int                   naggrrows,          /**< length of the arrays */
   SCIP_ROW**            negslackrows,       /**< array of rows that are added implicitly with negative slack */
   SCIP_Real*            negslackweights,    /**< array of negative slack weights */
   int                   nnegslackrows       /**< length of the negative slack array */
   )
{
   assert(scip != NULL);
   assert(scip->stat != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddCertificateAggrInfo", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcertificateNewAggrInfo(scip, aggrrow, aggrrows, weights, naggrrows, negslackrows, negslackweights, nnegslackrows) );

   return SCIP_OKAY;
}

/** stores the active aggregation information in the certificate data structures for a row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPstoreCertificateActiveAggrInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< row that aggregation-info is stored for */
   )
{
   SCIP_CERTIFICATE* certificate;
   SCIP_AGGREGATIONINFO* aggrinfo;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPstoreCertificateActiveAggrInfo", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPisExact(scip) || !SCIPisCertified(scip) )
      return SCIP_OKAY;

   certificate = SCIPgetCertificate(scip);

   assert(certificate != NULL);
   assert(certificate->workingaggrinfo);

   aggrinfo = certificate->aggrinfo[certificate->naggrinfos - 1];
   certificate->workingaggrinfo = FALSE;

   assert(aggrinfo != NULL);

   SCIP_CALL( SCIPhashmapSetImage(certificate->aggrinfohash, (void*) row, (void*) aggrinfo) );

   return SCIP_OKAY;
}

/** frees the active aggregation information
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPfreeCertificateActiveAggrInfo(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CERTIFICATE* certificate;
   SCIP_AGGREGATIONINFO* aggrinfo;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPfreeCertificateActiveAggrInfo", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPisExact(scip) || !SCIPisCertified(scip) )
      return SCIP_OKAY;

   certificate = SCIPgetCertificate(scip);

   assert(certificate != NULL);

   /* if the mirinfo is used it gets tranformed into sparse format, don't free it in that case */
   if( !certificate->workingaggrinfo )
      return SCIP_OKAY;

   aggrinfo = certificate->aggrinfo[certificate->naggrinfos - 1];

   assert(aggrinfo != NULL);

   SCIP_CALL( SCIPcertificateFreeAggrInfo(scip->set, certificate, scip->lp, aggrinfo, NULL) );

   certificate->workingaggrinfo = FALSE;

   return SCIP_OKAY;
}

/** adds aggregation information to certificate for one row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPaddCertificateMirInfo(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->stat != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddCertificateMirInfo", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcertificateNewMirInfo(scip) );

   return SCIP_OKAY;
}

/** stores the active mir information in the certificate data structures for a row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPstoreCertificateActiveMirInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< row that mirinfo is stored for */
   )
{
   SCIP_CERTIFICATE* certificate;
   SCIP_MIRINFO* mirinfo;
   int i;
   int csplit;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPstoreCertificateActiveMirInfo", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPisExact(scip) || !SCIPisCertified(scip) )
      return SCIP_OKAY;

   certificate = SCIPgetCertificate(scip);

   assert(certificate != NULL);
   assert(certificate->workingmirinfo);

   certificate->workingmirinfo = FALSE;

   mirinfo = certificate->mirinfo[certificate->nmirinfos - 1];

   assert(mirinfo != NULL);

   assert(mirinfo->nsplitvars == SCIPgetNVars(scip));

   csplit = 0;

   /* make the mirinfo sparse again */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mirinfo->varinds), mirinfo->nsplitvars) );

   for( i = 0; i < mirinfo->nsplitvars; i++ )
   {
      if( mirinfo->splitcoefficients[i] != 0.0 || mirinfo->upperused[i] || mirinfo->localbdused[i] )
      {
         mirinfo->splitcoefficients[csplit] = mirinfo->splitcoefficients[i];
         mirinfo->upperused[csplit] = mirinfo->upperused[i];
         mirinfo->localbdused[csplit] = mirinfo->localbdused[i];
         mirinfo->varinds[csplit] = i;
         csplit++;
      }
   }

   mirinfo->nsplitvars = csplit;
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(mirinfo->splitcoefficients), SCIPgetNVars(scip), mirinfo->nsplitvars) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(mirinfo->varinds), SCIPgetNVars(scip), mirinfo->nsplitvars) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(mirinfo->upperused), SCIPgetNVars(scip), mirinfo->nsplitvars) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(mirinfo->localbdused), SCIPgetNVars(scip), mirinfo->nsplitvars) );

   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &mirinfo->slackrows, SCIPgetNVars(scip), mirinfo->nslacks) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &mirinfo->slacksign, SCIPgetNVars(scip), mirinfo->nslacks) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &mirinfo->slackcoefficients, SCIPgetNVars(scip), mirinfo->nslacks) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &mirinfo->slackweight, SCIPgetNVars(scip), mirinfo->nslacks) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &mirinfo->slackscale, SCIPgetNVars(scip), mirinfo->nslacks) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &mirinfo->slackusedcoef, SCIPgetNVars(scip), mirinfo->nslacks) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &mirinfo->slackroundeddown, SCIPgetNVars(scip), mirinfo->nslacks) );

   SCIP_CALL( SCIPhashmapSetImage(certificate->mirinfohash, (void*) row, (void*) mirinfo) );

   return SCIP_OKAY;
}

/** print MIR cut to certificate file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcertifyMirCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< row that needs to be certified */
   )
{
   SCIP_CERTIFICATE* certificate;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcertifyMirCut", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPisExact(scip) || !SCIPisCertified(scip) )
      return SCIP_OKAY;

   certificate = SCIPgetCertificate(scip);
   SCIP_CALL( SCIPcertificatePrintMirCut(scip->set, scip->lp, certificate, scip->transprob, row, 'L') );

   return SCIP_OKAY;
}

/** frees the active mir information
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPfreeCertificateActiveMirInfo(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CERTIFICATE* certificate;
   SCIP_MIRINFO* mirinfo;
   int i;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPfreeCertificateActiveMirInfo", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPisExact(scip) || !SCIPisCertified(scip) )
      return SCIP_OKAY;

   certificate = SCIPgetCertificate(scip);

   assert(certificate != NULL);

   if( certificate->nmirinfos == 0 )
      return SCIP_OKAY;

   mirinfo = certificate->mirinfo[certificate->nmirinfos - 1];

   assert(mirinfo != NULL);

   /* if the mirinfo is used it gets tranformed into sparse format, don't free it in that case */
   if( !certificate->workingmirinfo )
      return SCIP_OKAY;

   for(i = 0; i < mirinfo->nslacks; ++i)
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(mirinfo->slackrows[i])) );
   }

   assert(mirinfo->varinds == NULL);

   SCIPrationalFreeBlock(SCIPblkmem(scip), &(mirinfo->rhs));
   SCIPrationalFreeBlock(SCIPblkmem(scip), &(mirinfo->frac));

   SCIPfreeBlockMemoryArray(scip, &(mirinfo->splitcoefficients), SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &(mirinfo->upperused), SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &(mirinfo->localbdused), SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &(mirinfo->slackrows), SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &(mirinfo->slackcoefficients), SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &(mirinfo->slacksign), SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &(mirinfo->slackweight), SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &(mirinfo->slackscale), SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &(mirinfo->slackusedcoef), SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &(mirinfo->slackroundeddown), SCIPgetNVars(scip));
   SCIPfreeBlockMemory(scip, &mirinfo);
   certificate->nmirinfos--;
   certificate->workingmirinfo = FALSE;

   return SCIP_OKAY;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

/**@file   certificate.h
 * @brief  methods for certificate output
 * @author Ambros Gleixner
 * @author Daniel Steffy
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CERTIFICATE_H__
#define __SCIP_CERTIFICATE_H__


#include "scip/def.h"
#include "scip/type_set.h"
#include "scip/type_cuts.h"
#include "scip/type_stat.h"
#include "scip/type_tree.h"
#include "scip/type_certificate.h"
#include "scip/type_sol.h"
#include "scip/type_lpexact.h"
#include "scip/type_cons.h"
#include "scip/type_var.h"
#include "scip/pub_fileio.h"
#include "scip/type_prob.h"
#ifdef SCIP_WITH_GMP
#include "gmp.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** creates certificate data structure */
SCIP_RETCODE SCIPcertificateCreate(
   SCIP_CERTIFICATE**    certificate,        /**< pointer to store the certificate information */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** frees certificate data structure */
void SCIPcertificateFree(
   SCIP_CERTIFICATE**    certificate         /**< pointer to store the certificate information */
   );

/** initializes certificate information and creates files for certificate output */
SCIP_RETCODE SCIPcertificateInit(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** initializes certificate information and creates files for certificate output */
SCIP_RETCODE SCIPcertificateInitTransFile(
   SCIP*                 scip                /**< scip data structure */
   );

/** returns whether the certificate output is activated? */
SCIP_Bool SCIPsetCertificateEnabled(
   SCIP_SET*             set                 /**< SCIP settings */
   );

/** returns current certificate index (return -1 if certificate not active) */
SCIP_Real SCIPcertificateGetFilesize(
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   );

/** returns current certificate index*/
SCIP_Longint SCIPcertificateGetCurrentIndex(
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   );

#ifndef NDEBUG
/** checks if information is consistent with printed certificate line */
SCIP_Bool SCIPcertificateEnsureLastBoundInfoConsistent(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_VAR*             var,                /**< variable that gets changed */
   SCIP_BOUNDTYPE        boundtype,          /**< lb or ub changed? */
   SCIP_Real             newbound,           /**< new bound */
   SCIP_Bool             needsglobal         /**< if the bound needs to be global */
   );
#endif

/** sets the objective function used when printing dual bounds */
SCIP_RETCODE SCIPcertificateSetAndPrintObjective(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the original solution be printed or in transformed space */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Rational**       coefs,              /**< objective function coefficients */
   int                   nvars               /**< number of variables */
   );

/** prints a string to the problem section of the certificate file */
void SCIPcertificatePrintProblemMessage(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the original solution be printed or in transformed space */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a string to the proof section of the certificate file */
void SCIPcertificatePrintProofMessage(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a rational number to the problem section of the certificate file */
void SCIPcertificatePrintProblemRational(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the original solution be printed or in transformed space */
   SCIP_Rational*        val,                /**< Rational to print to the problem*/
   int                   base                /**< The base representation*/
   );

/** prints a rational number to the proof section of the certificate file */
void SCIPcertificatePrintProofRational(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Rational*        val,                /**< Rational to print to the problem*/
   int                   base                /**< The base representation*/
   );

/** prints a comment to the problem section of the certificate file */
void SCIPcertificatePrintProblemComment(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the original solution be printed or in transformed space */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a comment to the proof section of the certificate file */
void SCIPcertificatePrintProofComment(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints variable section header */
void SCIPcertificatePrintVarHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the original solution be printed or in transformed space */
   int                   nvars               /**< number of variables */
   );

/** prints version header */
void SCIPcertificatePrintVersionHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile          /**< should the original solution be printed or in transformed space */
   );

/** prints integer section header */
void SCIPcertificatePrintIntHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the original solution be printed or in transformed space */
   int                   nvars               /**< number of variables */
   );

/** prints constraint section header */
void SCIPcertificatePrintConsHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the original solution be printed or in transformed space */
   int                   nconss,             /**< number of all constraints */
   int                   nboundconss         /**< number of bound constraints */
   );

/** prints derivation section header */
void SCIPcertificatePrintDerHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile          /**< shoud the line be printed to the origfile or the transfile */
   );

/** prints constraint */
void SCIPcertificatePrintCons(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the original solution be printed or in transformed space */
   const char*           consname,           /**< name of the constraint */
   const char            sense,              /**< sense of the constraint, i.e., G, L, or E */
   SCIP_Rational*        side,               /**< left/right-hand side */
   int                   len,                /**< number of nonzeros */
   int*                  ind,                /**< index array */
   SCIP_Rational**       val                 /**< coefficient array */
   );

/** prints verification of row as a MIR cut (viewed as a split cut) */
SCIP_RETCODE SCIPcertificatePrintMirCut(
   SCIP_SET*             set,                /**< SCIP settings */
   SCIP_LP*              lp,                 /**< SCIP lp data structure */
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_PROB*            prob,               /**< SCIP problem data */
   SCIP_ROW*             row,                /**< the row to be printed */
   const char            sense               /**< sense of the constraint, i.e., G, L, or E */
   );

/** create a new node data structure for the current node */
SCIP_RETCODE SCIPcertificateTransAggrrow(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_PROB*            prob,               /**< SCIP problem data */
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate */
   SCIP_AGGRROW*         aggrrow,            /**< agrrrow that results from the aggregation */
   SCIP_ROW*             row,                /**< the cut that we are attempting to prove */
   SCIP_ROW**            aggrrows,           /**< array of rows used fo the aggregation */
   SCIP_Real*            weights,            /**< array of weights */
   int                   naggrrows           /**< length of the arrays */
   );

/** Print cutoff bound for objective value **/
SCIP_RETCODE SCIPcertificatePrintCutoffBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate */
   SCIP_Rational*        bound,              /**< the bound */
   SCIP_Longint*         certificateline     /**< save the line index */
   );

/** create a new node data structure for the current node */
SCIP_RETCODE SCIPcertificatePrintAggrrow(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_LP*              lp,                 /**< SCIP lp data structure */
   SCIP_PROB*            prob,               /**< SCIP problem data */
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate */
   SCIP_AGGRROW*         aggrrow,            /**< agrrrow that results from the aggregation */
   SCIP_ROW**            aggrrows,           /**< array of rows used fo the aggregation */
   SCIP_Real*            weights,            /**< array of weights */
   int                   naggrrows,          /**< length of the arrays */
   SCIP_Bool             local,              /**< true if local bound information can be used */
   unsigned long*        certificateline     /**< pointer to store the certificate line index or NULL */
   );

/** prints a variable bound to the problem section of the certificate file and returns line index */
SCIP_RETCODE SCIPcertificatePrintBoundCons(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the original solution be printed or in transformed space */
   const char*           boundname,          /**< name of the bound constraint */
   SCIP_VAR*             var,                /**< variable to print the bound cons for */
   SCIP_Rational*        boundval,           /**< value of the bound */
   SCIP_Bool             isupper             /**< is it the upper bound? */
   );

/** installs updated node data in parent node */
SCIP_RETCODE SCIPcertificateUpdateParentData(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_NODE*            node,               /**< node data structure */
   SCIP_Longint          fileindex,          /**< index of new bound */
   SCIP_Rational*        newbound            /**< pointer to value of new bound, NULL if infeasible */
   );

/** prints dual bound to proof section */
SCIP_Longint SCIPcertificatePrintDualbound(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   const char*           linename,           /**< name of the unsplitting line */
   SCIP_Rational*        lowerbound,         /**< pointer to lower bound on the objective, NULL indicating infeasibility */
   int                   len,                /**< number of dual multipiers */
   SCIP_Longint*         ind,                /**< index array */
   SCIP_Rational**       val                 /**< array of dual multipliers */
   );

/** Print a dual bound from an exact lp solution */
SCIP_RETCODE SCIPcertificatePrintDualboundExactLP(
   SCIP_CERTIFICATE*     certificate,        /**< scip certificate struct */
   SCIP_LPEXACT*         lpexact,            /**< the exact lp */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_NODE*            node,               /**< the current node */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Bool             usefarkas           /**< should an infeasibility proof be printed? */
   );

/** Print a dual bound from an exact lp solution */
SCIP_RETCODE SCIPcertificatePrintDualboundPseudo(
   SCIP_CERTIFICATE*     certificate,        /**< scip certificate struct */
   SCIP_LPEXACT*         lpexact,            /**< the exact lp */
   SCIP_NODE*            node,               /**< current node */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Bool             lowerchanged,       /**< to the modified indices address a change in lb or ub? */
   int                   modifiedvarindex,   /**< index of modified variable, or -1 */
   SCIP_Longint          boundchangeindex,   /**< index of unprocessed bound change in the certificate, or -1 */
   SCIP_Real             psval               /**< the pseudo obj value */
   );

/** prints the bound that a node inherits from its parent to the certificate */
SCIP_RETCODE SCIPcertificatePrintInheritedBound(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_NODE*            node                /**< node data */
   );

/** update the parent certificate node data when branching, print branching into certificate if not already present */
SCIP_RETCODE SCIPcertificatePrintBranching(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp,                 /**< LP informations */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NODE*            node,               /**< node data */
   SCIP_VAR*             branchvar,          /**< the variable that gets branched on */
   SCIP_BOUNDTYPE        boundtype,          /**< the bounding type */
   SCIP_Real             newbound            /**< the new bound */
   );

/** create a new node data structure for the current node */
SCIP_RETCODE SCIPcertificateNewNodeData(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   );

/** create a new split info structure for the current cut */
SCIP_RETCODE SCIPcertificateNewMirInfo(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** free all aggregation information */
SCIP_RETCODE SCIPcertificateClearAggrinfo(
   SCIP*                 scip                /**< global SCIP data structure */
   );

/** free all mir information */
SCIP_RETCODE SCIPcertificateClearMirinfo(
   SCIP*                 scip                /**< global SCIP data structure */
   );

/** free aggregation information */
SCIP_RETCODE SCIPcertificateFreeAggrInfo(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate structure */
   SCIP_LP*              lp,                 /**< SCIP lp data structure */
   SCIP_AGGREGATIONINFO* aggrinfo,           /**< SCIP aggregation info */
   SCIP_ROW*             row                 /**< new row, that info should be stored for */
   );

/** free mir information for row */
SCIP_RETCODE SCIPcertificateFreeMirInfo(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate structure */
   SCIP_LP*              lp,                 /**< SCIP lp data structure */
   SCIP_MIRINFO*         mirinfo,           /**< SCIP mir info */
   SCIP_ROW*             row                 /**< row that should be freed, or NULL if not needed */
   );

/** create a new aggregation info for a row */
SCIP_RETCODE SCIPcertificateNewAggrInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< agrrrow that results from the aggregation */
   SCIP_ROW**            aggrrows,           /**< array of rows used fo the aggregation */
   SCIP_Real*            weights,            /**< array of weights */
   int                   naggrrows,          /**< length of the arrays */
   SCIP_ROW**            negslackrows,       /**< array of rows that are added implicitly with negative slack */
   SCIP_Real*            negslackweights,    /**< array of negative slack weights */
   int                   nnegslackrows       /**< length of the negative slack array */
   );

/** prints unsplitting information to proof section */
void SCIPcertificatePrintUnsplitting(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_NODE*            node                /**< node data */
   );

/** prints RTP section with lowerbound and upperbound range */
void SCIPcertificatePrintRtpRange(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_Bool             isorigfile,         /**< should the original solution be printed or in transformed space */
   SCIP_Rational*        lowerbound,         /**< pointer to lower bound on the objective */
   SCIP_Rational*        upperbound          /**< pointer to upper bound on the objective */
   );

/** prints the last part of the certificate header (RTP range/sol, ...) */
SCIP_RETCODE SCIPcertificatePrintResult(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             isorigfile,         /**< should the original solution be printed or in transformed space */
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   );

/** prints the last part of the certificate header (RTP range/sol, ...) */
SCIP_RETCODE SCIPcertificateSaveFinalbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   );

/** prints RTP section for infeasibility */
void SCIPcertificatePrintRtpInfeas(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_Bool             isorigfile          /**< should the original solution be printed or in transformed space */
   );

/** prints SOL header and exact solution to certificate file */
void SCIPcertificatePrintSolExact(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** set the node to have its own bound proof */
SCIP_RETCODE SCIPcertificateSetInheritanceData(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_NODE*            node,               /**< node data structure */
   SCIP_Longint          fileindex,          /**< index of new bound */
   SCIP_Rational*        newbound            /**< the inherited bound */
   );

/** sets the last bound index for the certificate */
SCIP_RETCODE SCIPcertificateSetLastBoundIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_Longint          index               /**< index of new bound */
   );

/** returns the last bound index for the certificate */
SCIP_Longint SCIPcertificateGetLastBoundIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CERTIFICATE*     certificate         /**< certificate data structure */
   );

/** returns the index for a row in the certificate */
SCIP_Longint SCIPcertificateGetRowIndex(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_ROWEXACT*        row,                /**< row to consider */
   SCIP_Bool             rhs                 /**< whether we want the index for the rhs or the lhs */
   );

/** prints a proof that boundchange is leads to infeasibility */
SCIP_RETCODE SCIPcertificatePrintCutoffConflictingBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Rational*        lb,                 /**< lower bound */
   SCIP_Rational*        ub,                 /**< upper bound */
   SCIP_Longint          lbindex,            /**< index of the lower bound */
   SCIP_Longint          ubindex             /**< index of the upper bound */
   );

/** prints a proof for a new global bound */
SCIP_RETCODE SCIPcertificatePrintGlobalBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate */
   SCIP_VAR*             var,                /**< variable */
   SCIP_BOUNDTYPE        boundtype,          /**< Whether we have an upper bound or a lower bound */
   SCIP_Rational*        value,              /**< value of the bound */
   SCIP_Longint          certificateindex    /**< index in the certificate */
   );

/** should the certificate track bound changes? */
SCIP_Bool SCIPcertificateShouldTrackBounds(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif

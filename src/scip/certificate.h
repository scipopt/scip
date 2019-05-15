/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   certificate.h
 * @brief  methods for certificate output
 * @author Ambros Gleixner
 * @author Daniel Steffy
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CERTIFICATE_H__
#define __SCIP_CERTIFICATE_H__


#include "scip/def.h"
#include "scip/type_set.h"
#include "scip/type_tree.h"
#include "scip/type_certificate.h"
#include "scip/type_solex.h"
#include "scip/pub_fileio.h"
#include "gmp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates certificate data structure */
extern
SCIP_RETCODE SCIPcertificateCreate(
   SCIP_CERTIFICATE**    certificate,        /**< pointer to store the certificate information */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** frees certificate data structure */
extern
void SCIPcertificateFree(
   SCIP_CERTIFICATE**    certificate         /**< pointer to store the certificate information */
   );

/** initializes certificate information and creates files for certificate output */
extern
SCIP_RETCODE SCIPcertificateInit(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** closes the certificate output files */
extern
void SCIPcertificateExit(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** returns whether the certificate output is activated? */
extern
SCIP_Bool SCIPcertificateIsActive(
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   );

/** returns current certificate file size in MB */
extern
SCIP_Real SCIPcertificateGetFilesize(
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   );

/** sets the objective function used when printing dual bounds */
extern
SCIP_RETCODE SCIPcertificateSetAndPrintObjective(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const mpq_t*          coefs,              /**< objective function coefficients */
   int                   nvars               /**< number of variables */
   );

/** prints a string to the problem section of the certificate file */
extern
void SCIPcertificatePrintProblemMessage(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a string to the proof section of the certificate file */
extern
void SCIPcertificatePrintProofMessage(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a rational number to the problem section of the certificate file */
extern
void SCIPcertificatePrintProblemRational(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const mpq_t           val,                /**< Rational to print to the problem*/
   int                   base                /**< The base representation*/
   );

/** prints a rational number to the proof section of the certificate file */
extern
void SCIPcertificatePrintProofRational(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const mpq_t           val,                /**< Rational to print to the problem*/
   int                   base                /**< The base representation*/
   );

/** prints an integer to the problem section of the certificate file */
extern
void SCIPcertificatePrintProblemInteger(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const mpz_t           val,                /**< Integer to print to the problem*/
   int                   base                /**< The base representation*/
   );

/** prints an integer to the proof section of the certificate file */
extern
void SCIPcertificatePrintProofInteger(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const mpz_t           val,                /**< Integer to print to the problem*/
   int                   base                /**< The base representation*/
   );

/** prints a comment to the problem section of the certificate file */
extern
void SCIPcertificatePrintProblemComment(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a comment to the proof section of the certificate file */
extern
void SCIPcertificatePrintProofComment(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints variable section header */
extern
void SCIPcertificatePrintVarHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   int                   nvars               /**< number of variables */
   );

/** prints version header */
extern
void SCIPcertificatePrintVersionHeader(
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   );

/** prints integer section header */
extern
void SCIPcertificatePrintIntHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   int                   nvars               /**< number of variables */
   );

/** prints constraint section header */
extern
void SCIPcertificatePrintConsHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   int                   nconss,             /**< number of all constraints */
   int                   nboundconss         /**< number of bound constraints */
   );

/** prints derivation section header */
extern
void SCIPcertificatePrintDerHeader(
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   );

/** prints constraint */
extern
void SCIPcertificatePrintCons(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const char*           consname,           /**< name of the constraint */
   const char            sense,              /**< sense of the constraint, i.e., G, L, or E */
   const mpq_t           side,               /**< left/right-hand side */
   int                   len,                /**< number of nonzeros */
   int*                  ind,                /**< index array */
   mpq_t*                val                 /**< coefficient array */
   );

/** prints a variable bound to the problem section of the certificate file and returns line index */
extern
SCIP_RETCODE SCIPcertificatePrintBoundCons(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const char*           boundname,          /**< name of the bound constraint */
   int                   varindex,           /**< index of the variable */
   const mpq_t           boundval,           /**< value of the bound */
   SCIP_Bool             isupper             /**< is it the upper bound? */
   );

/** checks whether variable bound assumption is present; prints it if not; returns index */
extern
SCIP_Longint SCIPcertificatePrintBoundAssumption(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const char*           assumptionname,     /**< name of the bound constraint */
   int                   varindex,           /**< index of the variable */
   const mpq_t           boundval,           /**< value of the bound */
   SCIP_Bool             isupper             /**< is it the upper bound? */
   );

/** prints dual bound to proof section */
extern
SCIP_Longint SCIPcertificatePrintDualbound(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           linename,           /**< name of the unsplitting line */
   const mpq_t*          lowerbound,         /**< pointer to lower bound on the objective, NULL indicating infeasibility */
   int                   len,                /**< number of dual multipiers */
   int*                  ind,                /**< index array */
   const mpq_t*          val                 /**< array of dual multipliers */
   );

/** prints unsplitting information to proof section */
extern
int SCIPcertificatePrintUnsplitting(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   const char*           linename,           /**< name of the unsplitting line */
   const mpq_t*          lowerbound,         /**< pointer to lower bound on the objective, NULL indicating infeasibility */
   int                   derindex_left,      /**< index of the first derivation */
   int                   assumptionindex_left,/**< index of the first unsplitting assumption */
   int                   derindex_right,     /**< index of the second derivation */
   int                   assumptionindex_right/**< index of the second unsplitting assumption */
   );

/** prints RTP section with lowerbound and upperbound range */
extern
void SCIPcertificatePrintRtpRange(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   const mpq_t*          lowerbound,         /**< pointer to lower bound on the objective, NULL if negative infinity */
   const mpq_t*          upperbound          /**< pointer to upper bound on the objective, NULL if positive infinity */
   );


/** prints RTP section for infeasibility */
extern
void SCIPcertificatePrintRtpInfeas(
   SCIP_CERTIFICATE*     certificate         /**< certificate data structure */
   );

/** prints SOL header and exact solution to certificate file */
extern
void SCIPcertificatePrintSolex(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOLEX*           sol                 /**< primal CIP solution */
   );

#ifdef __cplusplus
}
#endif

#endif

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

/**@file   struct_certificate.h
 * @brief  data structures for certificate output
 * @author Ambros Gleixner
 * @author Daniel Steffy
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_CERTIFICATE_H__
#define __SCIP_STRUCT_CERTIFICATE_H__

#include <stdio.h>

#include "scip/def.h"
#include "scip/type_certificate.h"
#include "scip/type_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data structure for hashing bounds of variables in a certificate file */
struct SCIP_CertificateBound
{
   SCIP_Longint          fileindex;          /**< index of this bound in the certificate file */
   int                   varindex;           /**< index of the variable */
   SCIP_Rational*        boundval;           /**< value of the bound */
   SCIP_Bool             isupper;            /**< is it the upper bound? */
};

/** certificate data structure */
struct SCIP_Certificate
{
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< message handler to use */
   SCIP_HASHTABLE*       varboundtable;      /**< hash table for mapping variable bounds to line index in file */
   struct SCIP_CertificateBound* workbound;  /**< temporary memory for hashing bound information */
   BMS_BLKMEM*           blkmem;             /**< SCIP block memory */
   SCIP_Longint          indexcounter;       /**< counter for line indices in file */
   SCIP_Longint          conscounter;        /**< counter for line indices in constraint section */
   SCIP_FILE*            file;               /**< file to store problem definition */
   SCIP_FILE*            derivationfile;     /**< file to store derivations temporarily */
   char*                 derivationfilename; /**< name of the derivation file */
   char*                 objstring;          /**< string for buffering the objective function */
   SCIP_Real             filesize;           /**< size of derivation file in MB */
};

#ifdef __cplusplus
}
#endif

#endif

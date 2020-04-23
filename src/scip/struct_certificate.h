/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
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

struct SCIP_Certnodedata
{
   SCIP_Longint          assumptionindex_self;/**< Line Index where assumption is printed */
   SCIP_Longint          assumptionindex_left;/**< Line Index of left branch assumption */
   SCIP_Longint          derindex_left;      /**< Line Index of derivation assuming assumption left */
   SCIP_Rational*        derbound_left;      /**< Bound of left derivation */
   SCIP_Longint          assumptionindex_right;/**< Line Index of right branch assumption */
   SCIP_Longint          derindex_right;     /**< Line Index of derivation assuming assumption right */
   SCIP_Rational*        derbound_right;     /**< Bound of right derivation */
   SCIP_Longint          derindex_inherit;   /**< Line index of bound inherited from parent */
   SCIP_Rational*        derbound_inherit;   /**< inherited bound */
   unsigned int          leftfilled:1;       /**< Is the left node filled ? */
   unsigned int          leftinfeas:1;       /**< Is the left node infeasible ? */
   unsigned int          rightfilled:1;      /**< Is the node right filled ? */
   unsigned int          rightinfeas:1;      /**< Is the node right infeasible ? */
   unsigned int          inheritedbound:1;   /**< did the node inherit its bound from its parent node? */
};

/** certificate data structure */
struct SCIP_Certificate
{
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< message handler to use */
   SCIP_HASHTABLE*       varboundtable;      /**< hash table for mapping variable bounds to line index in file */
   SCIP_CERTIFICATEBOUND** boundvals;        /**< array to store rationals in varboundtable to avoid memory leak */
   SCIP_Longint          boundvalsize;       /**< size of boundvals array */
   SCIP_Longint          nboundvals;         /**< number of elements in boundvals array */
   SCIP_HASHMAP*         nodedatahash;       /**< Hashmap storing pointer to data of each node */
   SCIP_CERTIFICATEBOUND* workbound;         /**< temporary memory for hashing bound information */
   BMS_BLKMEM*           blkmem;             /**< SCIP block memory */
   SCIP_Longint          indexcounter;       /**< counter for line indices in file */
   SCIP_Longint          conscounter;        /**< counter for line indices in constraint section */
   SCIP_FILE*            file;               /**< file to store problem definition */
   SCIP_FILE*            derivationfile;     /**< file to store derivations temporarily */
   char*                 derivationfilename; /**< name of the derivation file */
   char*                 objstring;          /**< string for buffering the objective function */
   SCIP_Real             filesize;           /**< size of derivation file in MB */
   SCIP_HASHMAP*         rowdatahash;        /**< Hashmap storing mapping between rows and file index */
   SCIP_Rational*        rootbound;          /**< the bound for the root node */
   SCIP_Longint          derindex_root;      /**< index of root bound in certificate */
   SCIP_Bool             rootinfeas;         /**< is the root node infeasible */
   SCIP_Bool             objintegral;        /**< is the objective always integral? copy this so we don't need the prob everywhere */
   SCIP_Rational**       vals;               /**< we maintain an array for solvals so we don't have to reallocate at every bounding call */
   int                   valssize;           /**< the size of the vals array */
};

#ifdef __cplusplus
}
#endif

#endif

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

/**@file   struct_certificate.h
 * @ingroup INTERNALAPI
 * @brief  data structures for certificate output
 * @author Ambros Gleixner
 * @author Daniel Steffy
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_CERTIFICATE_H__
#define __SCIP_STRUCT_CERTIFICATE_H__

#include <stdio.h>

#include "scip/def.h"
#include "scip/pub_fileio.h"
#include "scip/type_certificate.h"
#include "scip/type_var.h"
#include "scip/type_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data structure for hashing bounds of variables in a certificate file */
struct SCIP_CertificateBound
{
   int                   varindex;           /**< index of this bound in the certificate file */
   SCIP_RATIONAL*        boundval;           /**< value of the bound */
   SCIP_BOUNDTYPE        boundtype;          /**< is it the upper bound? */
   SCIP_Bool             isbound;            /**< is the last printed index a bound? if it is not, the other information is not useful */
   SCIP_Bool             isglobal;           /**< is the bound global? */
   SCIP_Longint          certificateindex;   /**< index of the bound in the certificate file */
};

/** data structure for storing necessary information to print verified aggregation of rows */
struct SCIP_AggregationInfo
{
   SCIP_AGGRROW*         aggrrow;            /**< aggregation row to be saved */
   SCIP_ROW**            aggrrows;           /**< array of rows used for the aggregation */
   SCIP_ROW**            negslackrows;       /**< array of rows that are implicitly added (using negative slack) */
   SCIP_Real*            weights;            /**< array of weights */
   SCIP_Real*            negslackweights;    /**< array of weights for the negslackrows */
   SCIP_Real*            substfactor;        /**< factor used in the substition of slack variables (weight)/(1-f0) */
   int                   naggrrows;          /**< length of the aggrrows array */
   int                   nnegslackrows;      /**< length of the negslackrows array */
   SCIP_Longint          fileindex;          /**< index of the aggregated row in the certificate file */
   SCIP_Longint          arpos;              /**< position in the aggrinfo array, so we can access it from the hashmap */
};

/** data structure for certifying MIR cut (splitcoefs, rhs, fractionality f, 1/1-f, scaling factor, which bounds to use) */
struct SCIP_MirInfo
{
   SCIP_Real*            splitcoefficients;  /**< coefficients in the split, saved in the complemented variable space */
   SCIP_Real*            slackcoefficients;  /**< coefficients for integer slacks that enter the split */
   SCIP_Real*            slackweight;     /**< continuous part of integer slack that needs to be accounted for */
   SCIP_Bool*            slackroundeddown;     /**< original part of integer slack that needs to be accounted for */
   SCIP_Real*            slackscale;         /**< original part of integer slack that needs to be accounted for */
   SCIP_Real*            slackusedcoef;      /**< coef that was actually used in the slack subsititution */
   SCIP_ROW**            slackrows;          /**< rows whos integer slack is in the split */
   int*                  varinds;            /**< indices of variables in split */
   int*                  slacksign;          /**< was rhs or lhs used for integer slacks? +1 -> rhs, -1 -> lhs*/
   SCIP_Bool*            upperused;          /**< TRUE if ub was used to complement variable, FALSE if lb was used */
   SCIP_Bool*            localbdused;        /**< TRUE if local bound was used to complement variable, FALSE if global was used */
   int                   nsplitvars;         /**< number of variables in the split */
   int                   nlocalvars;         /**< number of local bounds used in transformation */
   int                   nslacks;            /**< number of integer slacks in the split */
   int                   nrounddownslacks;
   SCIP_RATIONAL*        rhs;                /**< rhs of the split disjunction */
   SCIP_RATIONAL*        frac;               /**< fractionality of the rhs in the mir cut */
   SCIP_Longint          arpos;              /**< position in the mirinfo array, so we can access it from the hashmap */
   SCIP_INTERVAL         onedivoneminusf0;   /**< rounded value of 1/(1-f0) that was used in MIR procedure */
   SCIP_Real             scale;              /**< scaling factor that was used in cut-postprocessing */
   SCIP_Real             unroundedrhs;       /**< we need to save the rhs if we round down integral cuts for certification */
};

struct SCIP_Certnodedata
{
   SCIP_Longint          assumptionindex_self;/**< line index where node's last assumption is printed */
   SCIP_Longint          derindex_self;      /**< line index of node's own bound, initially inherited from parent */
   SCIP_RATIONAL*        derbound_self;      /**< node's own bound, initially inherited from parent */
   SCIP_Longint          assumptionindex_left;/**< line index of left branch assumption */
   SCIP_Longint          derindex_left;      /**< line index of derivation assuming assumption left */
   SCIP_RATIONAL*        derbound_left;      /**< bound of left derivation */
   SCIP_Longint          assumptionindex_right;/**< line index of right branch assumption */
   SCIP_Longint          derindex_right;     /**< line index of derivation assuming assumption right */
   SCIP_RATIONAL*        derbound_right;     /**< bound of right derivation */
   unsigned int          leftfilled:1;       /**< is the data for the left child node set? */
   unsigned int          leftinfeas:1;       /**< is the left node infeasible ? */
   unsigned int          rightfilled:1;      /**< is the data for the right child node set? */
   unsigned int          rightinfeas:1;      /**< is the node right infeasible ? */
   unsigned int          inheritedbound:1;   /**< did the node inherit its bound from its parent node? */
};

/** certificate data structure */
struct SCIP_Certificate
{
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< message handler to use */
   SCIP_HASHMAP*         nodedatahash;       /**< Hashmap storing pointer to data of each node */
   SCIP_HASHMAP*         aggrinfohash;       /**< Hashmap storing aggregation information of rows */
   SCIP_HASHMAP*         mirinfohash;        /**< Hashmap storing split disjunctions */
   SCIP_AGGREGATIONINFO** aggrinfo;          /**< array to store the aggregation info to avoid memory leaks */
   SCIP_MIRINFO**        mirinfo;            /**< array to store the split info to avoid memory leaks */
   SCIP_Longint          aggrinfosize;       /**< size of aggrinfo array */
   SCIP_Longint          naggrinfos;         /**< number of elements in aggrinfo array */
   SCIP_Longint          mirinfosize;        /**< size of mirinfo array */
   SCIP_Longint          nmirinfos;          /**< number of elements in mirinfo array */
   SCIP_CERTIFICATEBOUND* lastinfo;          /**< information on last printed certificate index */
   BMS_BLKMEM*           blkmem;             /**< SCIP block memory */
   SCIP_Longint          indexcounter;       /**< counter for line indices in file */
   SCIP_Longint          indexcounter_ori;   /**< counter for line indices in origial problem vipr file */
   SCIP_Longint          conscounter;        /**< counter for line indices in constraint section */
   SCIP_Longint          lastboundindex;     /**< place to store the last bound index to avoid having to add it to the signature of SCIPvarChgUbLocal, varProcessChgUbLocal */
   SCIP_FILE*            origfile;           /**< file to store original problem definition */
   SCIP_FILE*            transfile;          /**< file to store transformed problem (after presolving) */
   SCIP_FILE*            derivationfile;     /**< file to store derivations temporarily */
   SCIP_Bool             transfile_initialized; /**< boolean to store if the transfile has been initialized */
   char*                 derivationfilename; /**< name of the derivation file */
   char*                 origfilename;       /**< name of the original problem file */
   SCIP_Real             filesize;           /**< size of derivation file in MB */
   SCIP_Real             maxfilesize;        /**< maximum size of derivation file in MB (stop printing if exceeded) */
   SCIP_HASHMAP*         rowdatahash;        /**< Hashmap storing mapping between rows and file index */
   SCIP_RATIONAL*        rootbound;          /**< the bound for the root node */
   SCIP_RATIONAL*        finalbound;         /**< the final dual bound value */
   SCIP_Longint          derindex_root;      /**< index of root bound in certificate */
   SCIP_Bool             rootinfeas;         /**< is the root node infeasible */
   SCIP_Bool             objintegral;        /**< is the objective always integral? copy this so we don't need the prob everywhere */
   SCIP_Bool             workingmirinfo;     /**< true if mirinfo is under construction and not sparsely stored, false otherwise */
   SCIP_Bool             workingaggrinfo;    /**< true if aggrinfo is under construction (last entry not in hashmap), false otherwise */
   SCIP_RATIONAL**       vals;               /**< we maintain an array for solvals so we don't have to reallocate at every bounding call */
   int                   valssize;           /**< the size of the vals array */
};

#ifdef __cplusplus
}
#endif

#endif

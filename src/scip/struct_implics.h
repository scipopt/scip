/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_implics.h
 * @brief  datastructures for implications, variable bounds, and cliques
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_IMPLICS_H__
#define __SCIP_STRUCT_IMPLICS_H__


#include "scip/def.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_implics.h"

#ifdef __cplusplus
extern "C" {
#endif

/** variable bounds of a variable x in the form x <= b*z + d  or  x >= b*z + d */
struct SCIP_VBounds
{
   SCIP_VAR**            vars;               /**< variables z    in variable bounds x <= b*z + d  or  x >= b*z + d */
   SCIP_Real*            coefs;              /**< coefficients b in variable bounds x <= b*z + d  or  x >= b*z + d */
   SCIP_Real*            constants;          /**< constants d    in variable bounds x <= b*z + d  or  x >= b*z + d */
   int                   len;                /**< number of existing variable bounds (used slots in arrays) */
   int                   size;               /**< size of vars, coefs, and constants arrays */
};

/** implications for binary variable x in the form 
 *    x  <= 0  ==>  y <= b or y >= b (stored in arrays[0]) 
 *    x  >= 1  ==>  y <= b or y >= b (stored in arrays[1]) 
 *  implications with    binary variable y are stored at the beginning of arrays (sorted by pointer of y)
 *  implications with nonbinary variable y are stored at the end       of arrays (sorted by pointer of y)
 */
struct SCIP_Implics
{
   SCIP_VAR**            vars[2];            /**< variables y  in implications y <= b or y >= b */
   SCIP_BOUNDTYPE*       types[2];           /**< types        of implications y <= b (SCIP_BOUNDTYPE_UPPER)
                                              *                             or y >= b (SCIP_BOUNDTYPE_LOWER) */
   SCIP_Real*            bounds[2];          /**< bounds b     in implications y <= b or y >= b */
   int*                  ids[2];             /**< unique ids of implications; < 0 iff implication is a shortcut,
                                              *   i.e., it was added as part of the transitive closure of another implication */
   int                   size[2];            /**< size of implvars, implbounds and implvals arrays  for x <= 0 and x >= 1*/
   int                   nimpls[2];          /**< number of all implications                        for x <= 0 and x >= 1 */
   int                   nbinimpls[2];       /**< number of     implications with binary variable y for x <= 0 and x >= 1 */
};

/** single clique, stating that at most one of the binary variables can be fixed to the corresponding value */
struct SCIP_Clique
{
   SCIP_VAR**            vars;               /**< variables in the clique */
   SCIP_Bool*            values;             /**< values of the variables in the clique */
   int                   nvars;              /**< number of variables in the clique */
   int                   size;               /**< size of vars and values arrays */
   unsigned int          id:31;              /**< unique identifier of clique */
   unsigned int          eventsissued:1;     /**< were the IMPLADDED events on the variables already issued? */
};

/** list of cliques for a single variable */
struct SCIP_CliqueList
{
   SCIP_CLIQUE**         cliques[2];         /**< cliques the variable fixed to FALSE/TRUE is member of */
   int                   ncliques[2];        /**< number of cliques the variable fixed to FALSE/TRUE is member of */
   int                   size[2];            /**< size of cliques arrays */
};

/** collection of cliques */
struct SCIP_CliqueTable
{
   SCIP_CLIQUE**         cliques;            /**< cliques stored in the table */
   int                   ncliques;           /**< number of cliques stored in the table */
   int                   size;               /**< size of cliques array */
   int                   ncreatedcliques;    /**< number of ever created cliques */
   int                   ncleanupfixedvars;  /**< number of fixed variables when the last cleanup was performed */
   int                   ncleanupaggrvars;   /**< number of aggregated variables when the last cleanup was performed */
   int                   ncleanupcliques;    /**< number of cliques stored when the last cleanup was performed */
};

#ifdef __cplusplus
}
#endif

#endif

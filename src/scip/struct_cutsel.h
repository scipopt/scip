/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_cutsel.h
 * @ingroup INTERNALAPI
 * @brief  data structures for cut selectors
 * @author Felipe Serrano
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_CUTSEL_H__
#define __SCIP_STRUCT_CUTSEL_H__


#include "scip/def.h"
#include "scip/type_cutsel.h"

#ifdef __cplusplus
extern "C" {
#endif

/** cut selector */
struct SCIP_Cutsel
{
   char*                 name;               /**< name of cut selector */
   char*                 desc;               /**< description of cut selector */
   SCIP_DECL_CUTSELCOPY ((*cutselcopy));     /**< copy method of cut selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CUTSELFREE ((*cutselfree));     /**< destructor of cut selector */
   SCIP_DECL_CUTSELINIT ((*cutselinit));     /**< initialize cut selector */
   SCIP_DECL_CUTSELEXIT ((*cutselexit));     /**< deinitialize cut selector */
   SCIP_DECL_CUTSELINITSOL((*cutselinitsol));/**< solving process initialization method of cut selector */
   SCIP_DECL_CUTSELEXITSOL((*cutselexitsol));/**< solving process deinitialization method of cut selector */
   SCIP_DECL_CUTSELSELECT((*cutselselect));  /**< cut selection method */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this cut selector for the next stages */
   SCIP_CLOCK*           cutseltime;         /**< cut selector execution time */
   SCIP_CUTSELDATA*      cutseldata;         /**< cut selector data */
   int                   priority;           /**< priority of the cut selector */
   SCIP_Bool             initialized;        /**< is cut selector initialized? */
};

#ifdef __cplusplus
}
#endif

#endif

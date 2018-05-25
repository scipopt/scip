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

/**@file   scip_datastructures.c
 * @brief  public methods for data structures
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

#include "scip/scip_datastructures.h"
#include "scip/scip_mem.h"

#include "scip/pub_message.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif

/** creates a dynamic array of real values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPcreateRealarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY**      realarray           /**< pointer to store the real array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPrealarrayCreate(realarray, SCIPblkmem(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of real values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPfreeRealarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY**      realarray           /**< pointer to the real array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPrealarrayFree(realarray) );

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPextendRealarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPrealarrayExtend(realarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, minidx, maxidx) );

   return SCIP_OKAY;
}

/** clears a dynamic real array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPclearRealarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPrealarrayClear(realarray) );

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array
 *
 *  @return  value of entry in dynamic array
 */
SCIP_Real SCIPgetRealarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   idx                 /**< array index to get value for */
   )
{
   assert(scip != NULL);

   return SCIPrealarrayGetVal(realarray, idx);
}

/** sets value of entry in dynamic array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPsetRealarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   idx,                /**< array index to set value for */
   SCIP_Real             val                 /**< value to set array index to */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPrealarraySetVal(realarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, idx, val) );

   return SCIP_OKAY;
}

/** increases value of entry in dynamic array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPincRealarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   idx,                /**< array index to increase value for */
   SCIP_Real             incval              /**< value to increase array index */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPrealarrayIncVal(realarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, idx, incval) );

   return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements
 *
 *  @return the minimal index of all stored non-zero elements
 */
int SCIPgetRealarrayMinIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   )
{
   assert(scip != NULL);

   return SCIPrealarrayGetMinIdx(realarray);
}

/** returns the maximal index of all stored non-zero elements
 *
 *  @return the maximal index of all stored non-zero elements
 */
int SCIPgetRealarrayMaxIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   )
{
   assert(scip != NULL);

   return SCIPrealarrayGetMaxIdx(realarray);
}

/** creates a dynamic array of int values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPcreateIntarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY**       intarray            /**< pointer to store the int array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPintarrayCreate(intarray, SCIPblkmem(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of int values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPfreeIntarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY**       intarray            /**< pointer to the int array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPintarrayFree(intarray) );

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPextendIntarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPintarrayExtend(intarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, minidx, maxidx) );

   return SCIP_OKAY;
}

/** clears a dynamic int array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPclearIntarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPintarrayClear(intarray) );

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array
 *
 *  @return value of entry in dynamic array
 */
int SCIPgetIntarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   idx                 /**< array index to get value for */
   )
{
   assert(scip != NULL);

   return SCIPintarrayGetVal(intarray, idx);
}

/** sets value of entry in dynamic array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPsetIntarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   idx,                /**< array index to set value for */
   int                   val                 /**< value to set array index to */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPintarraySetVal(intarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, idx, val) );

   return SCIP_OKAY;
}

/** increases value of entry in dynamic array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPincIntarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   idx,                /**< array index to increase value for */
   int                   incval              /**< value to increase array index */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPintarrayIncVal(intarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, idx, incval) );

   return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements
 *
 *  @return the minimal index of all stored non-zero elements
 */
int SCIPgetIntarrayMinIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   )
{
   assert(scip != NULL);

   return SCIPintarrayGetMinIdx(intarray);
}

/** returns the maximal index of all stored non-zero elements
 *
 *  @return the maximal index of all stored non-zero elements
 */
int SCIPgetIntarrayMaxIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   )
{
   assert(scip != NULL);

   return SCIPintarrayGetMaxIdx(intarray);
}

/** creates a dynamic array of bool values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPcreateBoolarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY**      boolarray           /**< pointer to store the bool array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPboolarrayCreate(boolarray, SCIPblkmem(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of bool values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPfreeBoolarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY**      boolarray           /**< pointer to the bool array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPboolarrayFree(boolarray) );

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPextendBoolarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPboolarrayExtend(boolarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, minidx, maxidx) );

   return SCIP_OKAY;
}

/** clears a dynamic bool array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPclearBoolarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPboolarrayClear(boolarray) );

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array
 *
 *  @return value of entry in dynamic array at position idx
 */
SCIP_Bool SCIPgetBoolarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   idx                 /**< array index to get value for */
   )
{
   assert(scip != NULL);

   return SCIPboolarrayGetVal(boolarray, idx);
}

/** sets value of entry in dynamic array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPsetBoolarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   idx,                /**< array index to set value for */
   SCIP_Bool             val                 /**< value to set array index to */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPboolarraySetVal(boolarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, idx, val) );

   return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements
 *
 *  @return the minimal index of all stored non-zero elements
 */
int SCIPgetBoolarrayMinIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   assert(scip != NULL);

   return SCIPboolarrayGetMinIdx(boolarray);
}

/** returns the maximal index of all stored non-zero elements
 *
 *  @return the maximal index of all stored non-zero elements
 */
int SCIPgetBoolarrayMaxIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   assert(scip != NULL);

   return SCIPboolarrayGetMaxIdx(boolarray);
}

/** creates a dynamic array of pointers
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPcreatePtrarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY**       ptrarray            /**< pointer to store the int array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPptrarrayCreate(ptrarray, SCIPblkmem(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of pointers
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPfreePtrarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY**       ptrarray            /**< pointer to the int array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPptrarrayFree(ptrarray) );

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPextendPtrarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic int array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPptrarrayExtend(ptrarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, minidx, maxidx) );

   return SCIP_OKAY;
}

/** clears a dynamic pointer array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPclearPtrarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray            /**< dynamic int array */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPptrarrayClear(ptrarray) );

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
void* SCIPgetPtrarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic int array */
   int                   idx                 /**< array index to get value for */
   )
{
   assert(scip != NULL);

   return SCIPptrarrayGetVal(ptrarray, idx);
}

/** sets value of entry in dynamic array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPsetPtrarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic int array */
   int                   idx,                /**< array index to set value for */
   void*                 val                 /**< value to set array index to */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPptrarraySetVal(ptrarray, scip->set->mem_arraygrowinit, scip->set->mem_arraygrowfac, idx, val) );

   return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements
 *
 *  @return the minimal index of all stored non-zero elements
 */
int SCIPgetPtrarrayMinIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   )
{
   assert(scip != NULL);

   return SCIPptrarrayGetMinIdx(ptrarray);
}

/** returns the maximal index of all stored non-zero elements
 *
 *  @return the maximal index of all stored non-zero elements
 */
int SCIPgetPtrarrayMaxIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   )
{
   assert(scip != NULL);

   return SCIPptrarrayGetMaxIdx(ptrarray);
}

/** creates directed graph structure */
SCIP_RETCODE SCIPcreateDigraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIGRAPH**        digraph,            /**< pointer to store the created directed graph */
   int                   nnodes              /**< number of nodes */
   )
{
   assert(scip != NULL);
   assert(digraph != NULL);

   SCIP_CALL( SCIPdigraphCreate(digraph, scip->mem->probmem, nnodes) );

   return SCIP_OKAY;
}

/** copies directed graph structure
 *
 *  The copying procedure uses the memory of the passed SCIP instance. The user must ensure that the digraph lives
 *  as most as long as the SCIP instance.
 *
 *  @note The data in nodedata is copied verbatim. This possibly has to be adapted by the user.
 */
SCIP_RETCODE SCIPcopyDigraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIGRAPH**        targetdigraph,      /**< pointer to store the copied directed graph */
   SCIP_DIGRAPH*         sourcedigraph       /**< source directed graph */
   )
{
   assert(scip != NULL);
   assert(sourcedigraph != NULL);
   assert(targetdigraph != NULL);

   SCIP_CALL( SCIPdigraphCopy(targetdigraph, sourcedigraph, scip->mem->probmem) );

   return SCIP_OKAY;
}

/** creates a disjoint set (union find) structure \p uf for \p ncomponents many components (of size one) */
SCIP_RETCODE SCIPcreateDisjointset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DISJOINTSET**    djset,              /**< disjoint set (union find) data structure */
   int                   ncomponents         /**< number of components */
   )
{
   assert(scip != NULL);
   assert(djset != NULL);

   SCIP_CALL( SCIPdisjointsetCreate(djset, scip->mem->probmem, ncomponents) );

   return SCIP_OKAY;
}

/** frees the disjoint set (union find) data structure */
void SCIPfreeDisjointset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DISJOINTSET**    djset               /**< disjoint set (union find) data structure */
   )
{
   assert(scip != NULL);

   SCIPdisjointsetFree(djset, scip->mem->probmem);
}

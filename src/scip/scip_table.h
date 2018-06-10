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

/**@file   scip_table.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for statistics table plugins
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_TABLE_H__
#define __SCIP_SCIP_TABLE_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_set.h"
#include "scip/type_table.h"

/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons).
 * Additionally, the internal "set.h" is included, such that the defines in set.h are
 * available in optimized mode.
 */
#ifdef NDEBUG
#include "scip/struct_scip.h"
#include "scip/struct_stat.h"
#include "scip/set.h"
#include "scip/tree.h"
#include "scip/misc.h"
#include "scip/var.h"
#include "scip/cons.h"
#include "scip/solve.h"
#include "scip/debug.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicTableMethods
 *
 * @{
 */

/** creates a statistics table and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeTable(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of statistics table */
   const char*           desc,               /**< description of statistics table */
   SCIP_Bool             active,             /**< should the table be activated by default? */
   SCIP_DECL_TABLECOPY   ((*tablecopy)),     /**< copy method of statistics table or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_TABLEFREE   ((*tablefree)),     /**< destructor of statistics table */
   SCIP_DECL_TABLEINIT   ((*tableinit)),     /**< initialize statistics table */
   SCIP_DECL_TABLEEXIT   ((*tableexit)),     /**< deinitialize statistics table */
   SCIP_DECL_TABLEINITSOL ((*tableinitsol)), /**< solving process initialization method of statistics table */
   SCIP_DECL_TABLEEXITSOL ((*tableexitsol)), /**< solving process deinitialization method of statistics table */
   SCIP_DECL_TABLEOUTPUT ((*tableoutput)),   /**< output method */
   SCIP_TABLEDATA*       tabledata,          /**< statistics table data */
   int                   position,           /**< position of statistics table */
   SCIP_STAGE            earlieststage       /**< output of the statistics table is only printed from this stage onwards */
   );

/** returns the statistics table of the given name, or NULL if not existing */
EXTERN
SCIP_TABLE* SCIPfindTable(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of statistics table */
   );

/** returns the array of currently available statistics tables */
EXTERN
SCIP_TABLE** SCIPgetTables(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available statistics tables */
EXTERN
int SCIPgetNTables(
   SCIP*                 scip                /**< SCIP data structure */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif

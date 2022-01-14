/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   stpprioqueue.h
 * @brief  priority queue with integer keys
 * @author Daniel Rehfeldt
 *
 * Implements a (minimum) priority queue with integer keys.
 * NOTE: for efficiency reasons we don't want to give a compare callback,
 * as for example done in the SCIP default priority queue.
 * Also, this implementation is faster than the SCIP default priority queue.
 * todo: if needed for other key type, either use template pattern
 * or intrusive design with macros...
 */


#ifndef APPLICATIONS_STP_SRC_STPPRIOQUEUE_H_
#define APPLICATIONS_STP_SRC_STPPRIOQUEUE_H_

typedef struct stp_priority_queue STP_PQ;

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

extern SCIP_RETCODE stpprioqueue_create(SCIP*, int, STP_PQ**);
extern void         stpprioqueue_free(SCIP*, STP_PQ**);
extern void         stpprioqueue_deleteMin(void**, int*, STP_PQ*);
extern void*        stpprioqueue_deleteMinReturnData(STP_PQ*);
extern const void*  stpprioqueue_peakMinData(const STP_PQ*);
extern int          stpprioqueue_peakMinKey(const STP_PQ*);
extern SCIP_RETCODE stpprioqueue_insert(SCIP*, void*, int, STP_PQ*);
extern void         stpprioqueue_clean(STP_PQ*);
extern SCIP_Bool    stpprioqueue_isClean(const STP_PQ*);

#ifdef __cplusplus
}
#endif

#endif /* APPLICATIONS_STP_SRC_STPPRIOQUEUE_H_ */

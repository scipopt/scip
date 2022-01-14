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

/**@file   mincut.h
 * @brief  Minimum cut routines for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file encompasses minimum cut routines for Steiner tree problems.
 *
 */

#ifndef APPLICATIONS_STP_SRC_MINCUT_H_
#define APPLICATIONS_STP_SRC_MINCUT_H_

#include "graph.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct terminal_separator_storage TERMSEPAS;

extern SCIP_RETCODE mincut_termsepasInit(SCIP*, const GRAPH*, int, int, TERMSEPAS**);
extern void         mincut_termsepasFree(SCIP*, TERMSEPAS**);
extern int          mincut_termsepasGetNall(const TERMSEPAS*);
extern int          mincut_termsepasGetN(const TERMSEPAS*, int);
extern const int*   mincut_termsepasGetFirst(int, TERMSEPAS*, int*, int*);
extern const int*   mincut_termsepasGetNext(int, TERMSEPAS*, int*, int*);
extern int          mincut_termsepasGetSource(const TERMSEPAS*);
extern SCIP_Bool mincut_findTerminalSeparatorsIsPromising(const GRAPH*);
extern SCIP_RETCODE mincut_findTerminalSeparators(SCIP*, SCIP_RANDNUMGEN*, GRAPH*, TERMSEPAS*);
extern SCIP_RETCODE mincut_separateLp(SCIP*, SCIP_CONSHDLR*, SCIP_RANDNUMGEN*, const int*, GRAPH*, int, int*);


#ifdef __cplusplus
}
#endif

#endif /* APPLICATIONS_STP_SRC_MINCUT_H_ */

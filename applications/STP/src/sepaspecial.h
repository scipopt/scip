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

/**@file   sepaspecial.h
 * @brief  Separator for Steiner tree problem contraints beyond flow-balance-directed-cut constraints
 * @author Daniel Rehfeldt
 *
 * This file includes some special separator routines beyond the flow-balance directed cut formulation constraints.
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef APPLICATIONS_STP_SRC_SEPASPECIAL_H_
#define APPLICATIONS_STP_SRC_SEPASPECIAL_H_

#include "scip/scip.h"
#include "graph.h"


#ifdef __cplusplus
extern "C" {
#endif



/** cuts for pseudo ancestor cliques */
typedef struct pseudoancestor_cliques PACLIQUES;

/** cuts for implications between potential terminals */
typedef struct prize_collecting_implications PCIMPLICATION;

/** cuts for implications between non-terminals and terminals */
typedef struct vertex_terminal_implications VTIMPLICATION;


extern SCIP_RETCODE sepaspecial_pacliquesInit(SCIP*, const GRAPH*, PACLIQUES**);
extern void         sepaspecial_pacliquesFree(SCIP*, PACLIQUES**);
extern SCIP_RETCODE sepaspecial_pacliquesSeparate(SCIP*, SCIP_CONSHDLR*, PACLIQUES*, int, int*);
extern int          sepaspecial_pcimplicationsGetNstarts(const PCIMPLICATION*);
extern const int*   sepaspecial_pcimplicationsGetStarts(const PCIMPLICATION*);
extern const int*   sepaspecial_pcimplicationsGetVerts(const PCIMPLICATION*);
extern SCIP_RETCODE sepaspecial_pcimplicationsInit(SCIP*, const GRAPH*, PCIMPLICATION**);
extern void         sepaspecial_pcimplicationsFree(SCIP*, PCIMPLICATION**);
extern SCIP_RETCODE sepaspecial_pcimplicationsSeparate(SCIP*, SCIP_CONSHDLR*, PCIMPLICATION*, int, int*);
extern SCIP_RETCODE sepaspecial_vtimplicationsInit(SCIP*, const GRAPH*, VTIMPLICATION**);
extern void         sepaspecial_vtimplicationsFree(SCIP*, VTIMPLICATION**);
extern SCIP_RETCODE sepaspecial_vtimplicationsSeparate(SCIP*, SCIP_CONSHDLR*, VTIMPLICATION*, int, int*);



#ifdef __cplusplus
}
#endif


#endif /* APPLICATIONS_STP_SRC_SEPASPECIAL_H_ */

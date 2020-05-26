/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   graphdefs.h
 * @brief  includes graph definitions used for Steiner tree problems
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/



#ifndef APPLICATIONS_STP_SRC_GRAPHDEFS_H_
#define APPLICATIONS_STP_SRC_GRAPHDEFS_H_


#define STP_SPG                      0
#define STP_SAP                      1
#define STP_PCSPG                    2
#define STP_RPCSPG                   3
#define STP_NWSPG                    4
#define STP_DCSTP                    5
#define STP_NWPTSPG                  6
#define STP_RSMT                     7
#define STP_OARSMT                   8
#define STP_MWCSP                    9
#define STP_DHCSTP                   10
#define STP_GSTP                     11
#define STP_RMWCSP                   12
#define STP_BRMWCSP                  13

#define EAT_FREE     -1
#define EAT_LAST     -2
#define EAT_HIDE     -3

#define STP_TERM           0        /**< terminal */
#define STP_TERM_NONE     -1        /**< non-terminal */
#define STP_TERM_PSEUDO   -2        /**< pseudo-terminal (for PC/MW variants) */
#define STP_TERM_NONLEAF  -3        /**< non-leaf (pseudo-) terminal (for PC/MW variants) */

#define STP_CENTER_OK    0           /**< do nothing */
#define STP_CENTER_DEG   1           /**< find maximum degree */
#define STP_CENTER_SUM   2           /**< find the minimum distance sum */
#define STP_CENTER_MIN   3           /**< find the minimum largest distance */
#define STP_CENTER_ALL   4           /**< find the minimum distance sum to all knots */

#define TERM2EDGE_NOTERM      -1    /**< for PC/MW: vertex is no terminal */
#define TERM2EDGE_FIXEDTERM   -2    /**< for PC/MW: vertex is fixed terminal; artificial root is also considered a fixed terminal */
#define TERM2EDGE_NONLEAFTERM -3    /**< for PC/MW: vertex is non-leaf terminal */

#define SDSTAR_BASE_UNSET  -1
#define SDSTAR_BASE_KILLED -2

#define STP_DELPSEUDO_MAXGRAD   7
#define STP_DELPSEUDO_MAXNEDGES 21


/* ((((edge) % 2) == 0) ? ((edge) + 1) : ((edge) - 1)) without branch */
#define flipedge(edge) ( ((edge) + 1) - 2 * ((edge) % 2) )
#define flipedge_Uint(edge) ( (((unsigned int) edge) + 1) - 2 * (((unsigned int) edge) % 2) )

#define CONNECT      0
#define UNKNOWN    (-1)
#define FARAWAY            1e15
#define BLOCKED            1e10              /**< used for temporarily blocking an edge */
#define BLOCKED_MINOR      (BLOCKED - 1.0)   /**< used for permanently blocking an edge;
                                              * different from BLOCKED because of weird prize sum in reduce_base.c */

#define EDGE_BLOCKED       0
#define EDGE_MODIFIABLE    1

#define MST_MODE   0
#define FSP_MODE   1
#define BSP_MODE   2

#define Is_term(a)         ((a) >= 0)
#define Is_pseudoTerm(a)   ((a) == STP_TERM_PSEUDO)
#define Is_nonleafTerm(a)  ((a) == STP_TERM_NONLEAF)
#define Is_anyTerm(a)      ((a) >= 0 || (a) == STP_TERM_PSEUDO || (a) == STP_TERM_NONLEAF )
#define Edge_anti(a) ((((a) % 2) > 0) ? (a) - 1 : (a) + 1)

/* stp file format */
#define STP_FILE_MAGIC       0x33d32945
#define STP_FILE_VERSION_MAJOR   1
#define STP_FILE_VERSION_MINOR   0

typedef enum { FF_BEA, FF_STP, FF_PRB, FF_GRD } FILETYPE;


#endif /* APPLICATIONS_STP_SRC_GRAPHDEFS_H_ */

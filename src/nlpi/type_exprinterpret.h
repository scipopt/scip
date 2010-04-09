/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_exprinterpret.h,v 1.1 2010/04/09 20:55:02 bzfviger Exp $"

/**@file   type_exprinterpret.h
 * @brief  type definitions for expression interpreter
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_EXPRINTERPRET_H__
#define __SCIP_TYPE_EXPRINTERPRET_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_ExprInt     SCIP_EXPRINT;      /**< an expression interpreter */
typedef struct SCIP_ExprIntData SCIP_EXPRINTDATA;  /**< data of an expression interpreter */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_TYPE_EXPRINTERPRET_H__ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_message.h,v 1.1 2005/07/15 17:20:15 bzfpfend Exp $"

/**@file   pub_message.h
 * @brief  public methods for message output
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_MESSAGE_H__
#define __SCIP_PUB_MESSAGE_H__


/** sets the user data of the message handler */
extern
RETCODE SCIPmessagehdlrSetData(
   MESSAGEHDLR*     messagehdlr,        /**< message handler; must not be NULL */
   MESSAGEHDLRDATA* messagehdlrdata     /**< new message handler data to attach to the handler */
   );

/** returns the user data of the message handler */
extern
MESSAGEHDLRDATA* SCIPmessagehdlrGetData(
   MESSAGEHDLR*     messagehdlr         /**< message handler */
   );


#endif


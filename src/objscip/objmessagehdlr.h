/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   objmessagehdlr.h
 * @brief  C++ wrapper for message handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJMESSAGEHDLR_H__
#define __SCIP_OBJMESSAGEHDLR_H__

#include "scip/scip.h"
#include "objscip/objcloneable.h"

namespace scip
{

/** C++ wrapper object for message handlers */
class ObjMessagehdlr
{
public:
   /** should the output be buffered up to the next newline? */
   const SCIP_Bool scip_bufferedoutput_;

   /** default constructor */
   ObjMessagehdlr(
      SCIP_Bool          bufferedoutput      /**< should the output be buffered up to the next newline? */
      )
      : scip_bufferedoutput_(bufferedoutput)
   {
   }

   /** destructor */
   virtual ~ObjMessagehdlr()
   {
   }

   /** error message print method of message handler
    *
    *  This method is invoked, if SCIP wants to display an error message to the screen or a file
    */
   virtual void scip_error(
      SCIP_MESSAGEHDLR*  messagehdlr,        /**< the message handler itself */
      FILE*              file,               /**< file stream to print into */
      const char*        msg                 /**< string to output into the file */
      ) = 0;

   /** warning message print method of message handler
    *
    *  This method is invoked, if SCIP wants to display a warning message to the screen or a file
    */
   virtual void scip_warning(
      SCIP_MESSAGEHDLR*  messagehdlr,        /**< the message handler itself */
      FILE*              file,               /**< file stream to print into */
      const char*        msg                 /**< string to output into the file */
      ) = 0;

   /** dialog message print method of message handler
    *
    *  This method is invoked, if SCIP wants to display a dialog message to the screen or a file
    */
   virtual void scip_dialog(
      SCIP_MESSAGEHDLR*  messagehdlr,        /**< the message handler itself */
      FILE*              file,               /**< file stream to print into */
      const char*        msg                 /**< string to output into the file */
      ) = 0;

   /** info message print method of message handler
    *
    *  This method is invoked, if SCIP wants to display an information message to the screen or a file
    */
   virtual void scip_info(
      SCIP_MESSAGEHDLR*  messagehdlr,        /**< the message handler itself */
      FILE*              file,               /**< file stream to print into */
      const char*        msg                 /**< string to output into the file */
      ) = 0;
};

} /* namespace scip */


   
/** creates the message handler for the given message handler object */
extern
SCIP_RETCODE SCIPcreateObjMessagehdlr(
   SCIP_MESSAGEHDLR**    messagehdlr,        /**< pointer to store the message handler */
   scip::ObjMessagehdlr* objmessagehdlr,     /**< message handler object */
   SCIP_Bool             deleteobject        /**< should the message handler object be deleted when message handler is freed? */
   );

/** destroys the message handler that was created by SCIPcreateObjMessagehdlr();
 *  if deleteobject was set to TRUE in SCIPcreateObjMessagehdlr(), the message handler object is deleted
 */
extern
SCIP_RETCODE SCIPfreeObjMessagehdlr(
   SCIP_MESSAGEHDLR**    messagehdlr         /**< pointer to the message handler */
   );

/** returns the message handler object for the given message handler */
extern
scip::ObjMessagehdlr* SCIPgetObjMessagehdlr(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

#endif

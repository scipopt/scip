/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: objreader.cpp,v 1.4 2004/09/23 15:46:30 bzfpfend Exp $"

/**@file   objreader.cpp
 * @brief  C++ wrapper for file readers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objreader.h"




/*
 * Data structures
 */

/** file reader data */
struct ReaderData
{
   scip::ObjReader* objreader;          /**< file reader object */
   Bool             deleteobject;       /**< should the reader object be deleted when reader is freed? */
};




/*
 * Callback methods of file reader
 */

/** destructor of file reader to free user data (called when SCIP is exiting) */
static
DECL_READERFREE(readerFreeObj)
{  /*lint --e{715}*/
   READERDATA* readerdata;

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   assert(readerdata->objreader != NULL);

   /* call virtual method of reader object */
   CHECK_OKAY( readerdata->objreader->scip_free(scip, reader) );

   /* free reader object */
   if( readerdata->deleteobject )
      delete readerdata->objreader;

   /* free reader data */
   delete readerdata;
   SCIPreaderSetData(reader, NULL);
   
   return SCIP_OKAY;
}


/** problem reading method of reader */
static
DECL_READERREAD(readerReadObj)
{  /*lint --e{715}*/
   READERDATA* readerdata;

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   assert(readerdata->objreader != NULL);

   /* call virtual method of reader object */
   CHECK_OKAY( readerdata->objreader->scip_read(scip, reader, filename, result) );

   return SCIP_OKAY;
}




/*
 * file reader specific interface methods
 */

/** creates the file reader for the given file reader object and includes it in SCIP */
RETCODE SCIPincludeObjReader(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjReader* objreader,          /**< file reader object */
   Bool             deleteobject        /**< should the reader object be deleted when reader is freed? */
   )
{
   READERDATA* readerdata;

   /* create file reader data */
   readerdata = new READERDATA;
   readerdata->objreader = objreader;
   readerdata->deleteobject = deleteobject;

   /* include file reader */
   CHECK_OKAY( SCIPincludeReader(scip, objreader->scip_name_, objreader->scip_desc_, objreader->scip_extension_,
         readerFreeObj, readerReadObj,
         readerdata) );

   return SCIP_OKAY;
}

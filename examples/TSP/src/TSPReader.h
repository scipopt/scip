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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: TSPReader.h,v 1.2 2005/03/16 10:46:11 bzfberth Exp $"

/**@file   TSPReader.h
 * @brief  C++ file reader for TSP data files
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TSPREADER_H__
#define __TSPREADER_H__

#include <iostream>
#include <fstream>
#include <string>
#include "gminucut.h"
#include "objscip/objscip.h"

namespace tsp
{

/** SCIP file reader for TSP data files */
class TSPReader : public scip::ObjReader
{
public:

   Bool round_lengths_;
  
   /** default constructor */
   TSPReader(SCIP* scip)
      : scip::ObjReader("tspreader", "file reader for TSP files", "tsp")
   {
      /* add TSP reader parameters */
    SCIPaddBoolParam(scip,
                  "reading/tspreader/round_lengths", "should lenghts of edges be rounded to nearest integer?", 
                  &round_lengths_, TRUE, NULL, NULL);
   }

   /** destructor */
   virtual ~TSPReader()
   {
   }

   /** destructor of file reader to free user data (called when SCIP is exiting) */
   virtual RETCODE scip_free(
      SCIP*         scip,               /**< SCIP data structure */
      READER*       reader              /**< the file reader itself */
      );
   
   /** problem reading method of reader
    *
    *  possible return values for *result:
    *  - SCIP_SUCCESS    : the reader read the file correctly and created an appropritate problem
    *  - SCIP_DIDNOTRUN  : the reader is not responsible for given input file
    *
    *  If the reader detected an error in the input file, it should return with RETCODE SCIP_READERR or SCIP_NOFILE.
    */
   virtual RETCODE scip_read(
      SCIP*         scip,               /**< SCIP data structure */
      READER*       reader,             /**< the file reader itself */
      const char*   filename,           /**< full path and name of file to read, or NULL if stdin should be used */
      RESULT*       result              /**< pointer to store the result of the file reading call */
      );

private:
   
   void TSPReader::getNodesFromFile(
      std::ifstream&      filedata, 
      double*      x_coords, 
      double*      y_coords, 
      GRAPH*        graph
      );

   bool TSPReader::checkValid(
      GRAPH*         graph, 
      std::string        name,
      std::string        type, 
      std::string        edgeweighttype, 
      int           nnodes
      );
};

} /* namespace tsp */

#endif

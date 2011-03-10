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

/**@file   ReaderTSP.h
 * @brief  C++ file reader for TSP data files
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TSPREADER_H__
#define __TSPREADER_H__

#include <iostream>
#include <fstream>
#include <string>
#include "GomoryHuTree.h"
#include "objscip/objscip.h"

namespace tsp
{

/** SCIP file reader for TSP data files */
class ReaderTSP : public scip::ObjReader
{
public:

   SCIP_Bool round_lengths_;
  
   /** default constructor */
   ReaderTSP(SCIP* scip)
      : scip::ObjReader(scip, "tspreader", "file reader for TSP files", "tsp")
   {
      /* add TSP reader parameters */
      SCIPaddBoolParam(scip,
         "reading/tspreader/round_lengths", "should lenghts of edges be rounded to nearest integer?", 
         &round_lengths_, FALSE, TRUE, NULL, NULL);
   }

   /** destructor */
   virtual ~ReaderTSP()
   {
   }

   /** destructor of file reader to free user data (called when SCIP is exiting) */
   virtual SCIP_RETCODE scip_free(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_READER*       reader              /**< the file reader itself */
      );
   
   /** problem reading method of reader
    *
    *  possible return values for *result:
    *  - SCIP_SUCCESS    : the reader read the file correctly and created an appropritate problem
    *  - SCIP_DIDNOTRUN  : the reader is not responsible for given input file
    *
    *  If the reader detected an error in the input file, it should return with RETCODE SCIP_READERR or SCIP_NOFILE.
    */
   virtual SCIP_RETCODE scip_read(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_READER*       reader,             /**< the file reader itself */
      const char*        filename,           /**< full path and name of file to read, or NULL if stdin should be used */
      SCIP_RESULT*       result              /**< pointer to store the result of the file reading call */
      );

   /** problem writing method of reader; NOTE: if the parameter "genericnames" is TRUE, then
    *  SCIP already set all variable and constraint names to generic names; therefore, this
    *  method should always use SCIPvarGetName() and SCIPconsGetName(); 
    *
    *  possible return values for *result:
    *  - SCIP_SUCCESS    : the reader read the file correctly and created an appropritate problem
    *  - SCIP_DIDNOTRUN  : the reader is not responsible for given input file
    *
    *  If the reader detected an error in the writing to the file stream, it should return
    *  with RETCODE SCIP_WRITEERROR.
    */
   virtual SCIP_RETCODE scip_write(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_READER*       reader,             /**< the file reader itself */
      FILE*              file,               /**< output file, or NULL if standard output should be used */
      const char*        name,               /**< problem name */
      SCIP_PROBDATA*     probdata,           /**< user problem data set by the reader */
      SCIP_Bool          transformed,        /**< TRUE iff problem is the transformed problem */

      SCIP_OBJSENSE      objsense,           /**< objective sense */
      SCIP_Real          objscale,           /**< scalar applied to objective function; external objective value is
                                              *   extobj = objsense * objscale * (intobj + objoffset) */
      SCIP_Real          objoffset,          /**< objective offset from bound shifting and fixing */
      SCIP_VAR**         vars,               /**< array with active variables ordered binary, integer, implicit, 
                                              *   continuous */
      int                nvars,              /**< number of mutable variables in the problem */
      int                nbinvars,           /**< number of binary variables */
      int                nintvars,           /**< number of general integer variables */
      int                nimplvars,          /**< number of implicit integer variables */
      int                ncontvars,          /**< number of continuous variables */
      SCIP_VAR**         fixedvars,          /**< array with fixed and aggregated variables */
      int                nfixedvars,         /**< number of fixed and aggregated variables in the problem */
      int                startnvars,         /**< number of variables existing when problem solving started */
      SCIP_CONS**        conss,              /**< array with constraints of the problem */
      int                nconss,             /**< number of constraints in the problem */
      int                maxnconss,          /**< maximum number of constraints existing at the same time */
      int                startnconss,        /**< number of constraints existing when problem solving started */
      SCIP_Bool          genericnames,       /**< using generic variable and constraint names? */
      SCIP_RESULT*       result              /**< pointer to store the result of the file reading call */
      );
   
private:
   
   /** parses the node list */ 
   void getNodesFromFile(
      std::ifstream&     filedata,           /**< filestream containing the data to extract */
      double*            x_coords,           /**< double array to be filled with the x-coordinates of the nodes */
      double*            y_coords,           /**< same for y-coordinates */
      GRAPH*             graph               /**< the graph which is to be generated by the nodes */
      );

   /** method asserting, that the file has had the correct format and everything was set correctly */
   bool checkValid(
      GRAPH*             graph,              /**< the constructed graph, schould not be NULL */ 
      std::string        name,               /**< the name of the file */
      std::string        type,               /**< the type of the problem, should be "TSP" */
      std::string        edgeweighttype,     /**< type of the edgeweights, should be "EUC_2D", "MAX_2D", or "MAN_2D" */
      int                nnodes              /**< dimension of the problem, should at least be one */
      );
   
   /** adds a variable to both halfedges and captures it for usage in the graph */
   SCIP_RETCODE addVarToEdges(
      SCIP*                 scip,               /**< SCIP data structure */
      GRAPHEDGE*            edge,               /**< an edge of the graph */
      SCIP_VAR*             var                 /**< variable corresponding to that edge */
      );

};

} /* namespace tsp */

#endif

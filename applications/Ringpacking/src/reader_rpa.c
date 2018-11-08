/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_rpa.c
 * @brief  Ringpacking problem reader
 * @author Benjamin Mueller
 *
 * This file implements the reader/parser used to read the ringpacking input data. For more details see \ref RINGPACKING_READER.
 *
 * @page RINGPACKING_READER Parsing the input format and creating the problem
 *
 * In the <code>data</code> directory you find a few data files which contain each one ringpacking problem. They have
 * the following structure. In the first line the name of the instance is stated. In the second line you find three
 * integer numbers. The first one gives you the number of different ring types \f$T\f$, the second and third the width
 * and height of the rectangles, respectively. The remaining lines each contain one integer and two floats which
 * together specify one ring type. The integer gives the demand and the floats correspond to the inner and outer radius
 * of the respective type.
 *
 * For parsing that data, we implemented a reader plugin for \SCIP. A reader has several callback methods and at least
 * one interface methods (the one including the reader into \SCIP). For our purpose we only implemented the \ref
 * READERREAD "READERREAD" callback and the interface method which adds the reader plugin to \SCIP.
 *
 * @section RINGPACKING_READERINCLUDE The SCIPincludeReaderRpa() interface method
 *
 * The interface method <code>SCIPincludeReaderRpa()</code> is called to add the reader plugin to \SCIP (see
 * cmain.c). This means \SCIP gets informed that this reader is available for reading input files. Therefore, the
 * function <code>SCIPincludeReader()</code> is called within this method which passes all necessary information of the
 * reader to SCIP. This information includes the name of the reader, a description, and the file extension for which the
 * file reader is in charge. In our case we selected the file extension "rpa". This means that all files which have
 * this file extension are passed to our reader for parsing. Besides these information the call
 * <code>SCIPincludeReader()</code> also passes for each callback of the reader a function pointers
 * (some of them might be NULL pointers). These function pointers are used by \SCIP to run the reader. For more
 * information about all available reader callbacks we refer to the \ref READER "How to add file readers" tutorial. In
 * the remaining section we restrict ourself to the callback <code>READERREAD</code> which is the only one we
 * implemented for the ringpacking example. All other callbacks are not required for this example.
 *
 * @section RINGPACKING_READERREAD The READERREAD callback method
 *
 * The READERREAD callback is in charge of parsing a file and creating the problem. To see the list of arguments this
 * functions gets to see the file type_reader.h in the source of \SCIP. The following arguments are of interest in our
 * case. First of all the \SCIP pointer, the file name, and the SCIP_RESULT pointer. The \SCIP pointer gives us the
 * current environment. The file name states the file which we should open and parse. Last but not least, the SCIP_RESULT
 * pointer is required to tell \SCIP if the parsing process was successfully or not. Note that in type_reader.h you also
 * find a list of allowable result values for the SCIP_RESULT pointer and the <code>SCIP_RETCODE</code> which is the
 * return value of this function.
 *
 * @subsection RINGPACKING_PARSING Parsing the problem
 *
 * The file can be opened and parsed with your favorite methods. In this case we are using the functionality provided by
 * \SCIP since this has some nice side effects. We are using the function SCIPfopen() which can besides standard
 * files also handle files which are packed. To find all files related to the parsing of a file, we refer to the file pub_misc.h
 * in the source of SCIP. Parsing the data out of the file is not that hard. Please look at the code and comments
 * therein for more details.
 *
 * @subsection RINGPACKING_CREATING Creating the problem
 *
 * After parsing the file the final task for the reader is to create the problem. In our case, we pass the collected data
 * to the \ref probdata_rpa.h "main problem data plugin". For this, we use the interface methods
 * SCIPprobdataCreate() which is provided by the
 * problem data plugin (see probdata_rpa.c). After that, the reader sets the result value for the SCIP_RESULT
 * pointer to <code>SCIP_SUCCESS</code> and returns with a proper <code>SCIP_RETCODE</code>.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "reader_rpa.h"
#include "probdata_rpa.h"

/**@name Reader properties
 *
 * @{
 */

#define READER_NAME             "rpareader"
#define READER_DESC             "file reader for ringpacking data format"
#define READER_EXTENSION        "rpa"

/**@} */

/* default values of parameters */
#define DEFAULT_VERIFICATION_NLPTILIMSOFT    1.0       /**< soft time limit for each verification NLP */
#define DEFAULT_VERIFICATION_NLPNODELIMSOFT  1000L     /**< soft node limit for each verification NLP */
#define DEFAULT_VERIFICATION_HEURTILIMSOFT   1.0       /**< soft time limit for heuristic verification */
#define DEFAULT_VERIFICATION_HEURITERLIMSOFT 100       /**< soft iteration limit for each heuristic verification */
#define DEFAULT_VERIFICATION_TOTALTILIMSOFT  1200.0    /**< total time limit for all verification problems during the enumeration */


/**@name Callback methods
 *
 * @{
 */

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadRpa)
{  /*lint --e{715}*/
   SCIP_FILE* file;
   SCIP_Real* rints;
   SCIP_Real* rexts;
   int* demands;
   SCIP_Bool error;
   char name[SCIP_MAXSTRLEN];
   char buffer[SCIP_MAXSTRLEN];
   SCIP_Real width;
   SCIP_Real height;
   SCIP_Real r_int;
   SCIP_Real r_ext;
   int demand;
   int ntypes;
   int nread;
   int lineno;
   int i;

   *result = SCIP_DIDNOTRUN;
   width = -1.0;
   height = -1.0;

   /* open file */
   file = SCIPfopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   lineno = 0;
   sprintf(name, "++ uninitialized ++");
   ntypes = 0;

   /* read problem dimension */
   if( !SCIPfeof(file) )
   {
      /* get next line */
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         return SCIP_READERROR;
      lineno++;

      /* parse instance name line */
      nread = sscanf(buffer, "%s", name);
      if( nread != 1 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         return SCIP_READERROR;
      }

      /* get next line */
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         return SCIP_READERROR;
      lineno++;

      /* parse dimension line */
      nread = sscanf(buffer, "%d %" SCIP_REAL_FORMAT " %" SCIP_REAL_FORMAT "\n", &ntypes, &width, &height);
      if( nread < 3 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         return SCIP_READERROR;
      }
   }

   SCIPdebugMessage("instance name = %s\n", name);
   SCIPdebugMessage("width = %e height = %e\n", MAX(width,height), MIN(width,height));

   /* allocate buffer memory for storing the demands, rints, rexts */
   SCIP_CALL( SCIPallocBufferArray(scip, &demands, ntypes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rints, ntypes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rexts, ntypes) );

   /* ring types */
   r_int = 0.0;
   r_ext = 0.0;
   demand = 0;
   i = 0;
   error = FALSE;

   while( !SCIPfeof(file) && !error )
   {
      /* get next line */
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         break;
      lineno++;

      /* parse the line */
      nread = sscanf(buffer, "%d %" SCIP_REAL_FORMAT " %" SCIP_REAL_FORMAT "\n", &demand, &r_int, &r_ext);
      if( nread == 0 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         error = TRUE;
         break;
      }

      if( r_int > r_ext )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: internal radius is greater than the external one\n");
         error = TRUE;
         break;
      }

      if( demand <= 0 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: demand has to be positive\n");
         error = TRUE;
         break;
      }

      demands[i] = demand;
      rints[i] = r_int;
      rexts[i] = r_ext;
      ++i;

      if( i == ntypes )
         break;
   }

   if( i < ntypes )
   {
      SCIPwarningMessage(scip, "found %d different types of rings, needed %d\n", i, ntypes);
      error = TRUE;
   }

   if( !SCIPisPositive(scip, width) || !SCIPisPositive(scip, height) )
   {
      SCIPwarningMessage(scip, "non-positive width and height = (%f, %f)!\n", width, height);
      error = TRUE;
   }

   if( !error )
   {
      /* sort rings by their external radii */
      SCIPsortDownRealRealInt(rexts, rints, demands, ntypes);

      /* create and set problem data */
      SCIP_CALL( SCIPprobdataCreate(scip, filename, demands, rints, rexts, ntypes, MAX(width,height), MIN(width,height)) );
      SCIP_CALL( SCIPprobdataSetupProblem(scip) );
   }

   (void)SCIPfclose(file);
   SCIPfreeBufferArray(scip, &rints);
   SCIPfreeBufferArray(scip, &rexts);
   SCIPfreeBufferArray(scip, &demands);

   if( error )
      return SCIP_READERROR;

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** includes the rpa file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderRpa(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create ringpacking reader data */
   readerdata = NULL;

   /* include ringpacking reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );
   assert(reader != NULL);

   /* add soft verification parameters */
   SCIP_CALL( SCIPaddRealParam(scip, "ringpacking/verification/nlptilimsoft", "soft time limit for verification NLP",
                               NULL, FALSE, DEFAULT_VERIFICATION_NLPTILIMSOFT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "ringpacking/verification/nlpnodelimsoft",
                                  "soft node limit for verification NLP", NULL, FALSE,
                                  DEFAULT_VERIFICATION_NLPNODELIMSOFT, 0L, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "ringpacking/verification/heurtilimsoft",
                               "soft time limit for heuristic verification", NULL, FALSE,
                               DEFAULT_VERIFICATION_HEURTILIMSOFT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "ringpacking/verification/heuriterlimsoft",
                              "soft iteration limit for heuristic verification", NULL, FALSE,
                              DEFAULT_VERIFICATION_HEURITERLIMSOFT, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "ringpacking/verification/totaltilimsoft",
                               "total time limit for all verification problems during the enumeration", NULL, FALSE,
                               DEFAULT_VERIFICATION_TOTALTILIMSOFT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadRpa) );

   return SCIP_OKAY;
}

/**@} */

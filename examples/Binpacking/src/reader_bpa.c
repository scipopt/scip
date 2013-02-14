/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_bpa.c
 * @brief  Binpacking problem reader file reader
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * This file implements the reader/parser used to read the binpacking input data. For more details see \ref READER.
 *
 * @page READER Parsing the input format and creating the problem
 *
 * In the <code>data</code> directory you find a few data files which contain each one binpacking problem. These data
 * files have the following structure. In the first line the name of the instance is stated. In the second line you find
 * three integer numbers. The first one gives you the capacity \f$\kappa\f$, the second the number of items, and the
 * last integer states the value of a known feasible solution. This means an upper bound on the number of needed
 * bins. The remaining lines give the size for each item.
 *
 * For parsing that data, we implemented a reader plugin for \SCIP. A reader has several callback methods and at least
 * one interface methods (the one including the reader into \SCIP). For our purpose we only implemented the \ref
 * READERREAD "READERREAD" callback and the interface method which adds the reader plugin to \SCIP.
 *
 * @section READERINCLUDE The SCIPincludeReaderBpa() interface method
 *
 * The interface method <code>SCIPincludeReaderBpa()</code> is called to add the reader plugin to \SCIP (see
 * cmain.c). This means \SCIP gets informed that this reader is available for reading input files. Therefore, the
 * function <code>SCIPincludeReader()</code> is called within this method which passes all necessary information of the
 * reader to SCIP. This information includes the name of the reader, a description, and the file extension for which the
 * file reader is in charge. In our case we selected the file extension "bpa". This means that all files which have
 * this file extension are passed to our reader for parsing. Besides these information the call
 * <code>SCIPincludeReader()</code> also passes for each callback of the reader a function pointers
 * (some of them might be NULL pointers). These function
 * pointers are used by \SCIP to run the reader. For more information about all available reader callbacks we refer to
 * the <a href="http://scip.zib.de/doc/html/READER.html">How to add file readers</a> tutorial. In the remaining section
 * we restrict ourself to the callback <code>READERREAD</code> which is the only one we implemented for the binpacking
 * example. All other callbacks are not required for this example.
 *
 * @section READERREAD The READERREAD callback method
 *
 * The READERREAD callback is in charge of parsing a file and creating the problem. To see the list of arguments this
 * functions gets see the file type_reader.h in the source of \SCIP. The following arguments are of interest in our
 * case. First of all the \SCIP pointer, the file name, and the SCIP_RESULT pointer. The \SCIP pointer gives us the
 * current environment. The file name states the file which we should open and parse. Last but not least, the SCIP_RESULT
 * pointer is required to tell \SCIP if the parsing process was successfully or
 * not. Note that in type_reader.h you also find a list of allowable result values for the SCIP_RESULT pointer and the
 * <code>SCIP_RETCODE</code> which is the return value of this function.
 *
 * @subsection PARSING Parsing the problem
 *
 * The file can be opened and parsed with your favorite methods. In this case we are using the functionality provided by
 * \SCIP since this has some nice side effects. We are using the function SCIPfopen() which can besides standard
 * files also handle files which are packed. To find all files related to the parsing of a file, we refer to the file pub_misc.h
 * in the source of SCIP. Parsing the data out of the file is not that hard. Please look at the code and comments
 * therein for more details.
 *
 * @subsection CREATING Creating the problem
 *
 * After parsing the file the final task for the reader is to create the problem. In our case, we pass the collected data
 * to the \ref probdata_binpacking.h "main problem data plugin". For this, we use the interface methods
 * SCIPprobdataCreate() which is provided by the
 * problem data plugin (see probdata_binpacking.c). After that, the reader sets the result value for the SCIP_RESULT
 * pointer to <code>SCIP_SUCCESS</code> and returns with a proper <code>SCIP_RETCODE</code>.
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/cons_setppc.h"

#include "probdata_binpacking.h"
#include "reader_bpa.h"

/**@name Reader properties
 *
 * @{
 */

#define READER_NAME             "bpareader"
#define READER_DESC             "file reader for binpacking data format"
#define READER_EXTENSION        "bpa"

/**@} */


/**@name Callback methods
 *
 * @{
 */

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadBpa)
{  /*lint --e{715}*/
   SCIP_FILE* file;
   SCIP_Longint* weights;
   int* ids;
   SCIP_Bool error;

   char name[SCIP_MAXSTRLEN];
   char format[16];
   char buffer[SCIP_MAXSTRLEN];
   int capacity;
   int nitems;
   int bestsolvalue;
   int nread;
   int weight;
   int nweights;
   int lineno;

   *result = SCIP_DIDNOTRUN;

   /* open file */
   file = SCIPfopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   lineno = 0;

   /* read problem name */
   if( !SCIPfeof(file) )
   {
      /* get next line */
      if( SCIPfgets(buffer, sizeof(buffer), file) == NULL )
         return SCIP_READERROR;
      lineno++;

      /* parse dimension line */
      sprintf(format, "%%%ds\n", SCIP_MAXSTRLEN);
      nread = sscanf(buffer, format, name);
      if( nread == 0 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         return SCIP_READERROR;
      }

      SCIPdebugMessage("problem name <%s>\n", name);
   }

   /* read problem dimension */
   if( !SCIPfeof(file) )
   {
      /* get next line */
      if( SCIPfgets(buffer, sizeof(buffer), file) == NULL )
         return SCIP_READERROR;
      lineno++;

      /* parse dimension line */
      nread = sscanf(buffer, "%d %d %d\n", &capacity, &nitems, &bestsolvalue);
      if( nread < 2 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         return SCIP_READERROR;
      }

      SCIPdebugMessage("capacity = <%d>, number of items = <%d>, best known solution = <%d>\n", capacity, nitems, bestsolvalue);
   }


   /* allocate buffer memory for storing the weights and ids temporary */
   SCIP_CALL( SCIPallocBufferArray(scip, &weights, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ids, nitems) );

   /* pasre weights */
   nweights = 0;
   error = FALSE;

   while( !SCIPfeof(file) && !error )
   {
      /* get next line */
      if( SCIPfgets(buffer, sizeof(buffer), file) == NULL )
         break;
      lineno++;

      /* parse the line */
      nread = sscanf(buffer, "%d\n", &weight);
      if( nread == 0 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         error = TRUE;
         break;
      }

      SCIPdebugMessage("found weight %d <%d>\n", nweights, weight);
      weights[nweights] = weight;
      ids[nweights] = nweights;
      nweights++;

      if( nweights == nitems )
         break;
   }

   if( nweights < nitems )
   {
      SCIPwarningMessage(scip, "set nitems from <%d> to <%d> since the file <%s> only contains <%d> weights\n", nitems, weights, filename, weights);
      nitems = nweights;
   }

   if( !error )
   {
      /* create a new problem in SCIP */
      SCIP_CALL( SCIPprobdataCreate(scip, name, ids, weights, nitems, (SCIP_Longint)capacity) );
   }

   (void)SCIPfclose(file);
   SCIPfreeBufferArray(scip, &ids);
   SCIPfreeBufferArray(scip, &weights);

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

/** includes the bpa file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderBpa(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create binpacking reader data */
   readerdata = NULL;

   /* include binpacking reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );
   assert(reader != NULL);

   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadBpa) );

   return SCIP_OKAY;
}

/**@} */

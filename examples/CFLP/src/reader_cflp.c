/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_cflp.c
 * @brief  Cflp problem reader file reader
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * This file implements the reader/parser used to read the cflp input data. For more details see \ref CFLP_READER.
 *
 * @page CFLP_READER Parsing the input format and creating the problem
 *
 * In the <code>data</code> directory you find a few data files which contain each one cflp problem. These data
 * files have the following structure. In the first line the name of the instance is stated. In the second line you find
 * three integer numbers. The first one gives you the capacity \f$\kappa\f$, the second the number of items, and the
 * last integer states the value of a known feasible solution. This means an upper bound on the number of needed
 * bins. The remaining lines give the size for each item.
 *
 * For parsing that data, we implemented a reader plugin for \SCIP. A reader has several callback methods and at least
 * one interface methods (the one including the reader into \SCIP). For our purpose we only implemented the \ref
 * READERREAD "READERREAD" callback and the interface method which adds the reader plugin to \SCIP.
 *
 * @section CFLP_READERINCLUDE The SCIPincludeReaderCflp() interface method
 *
 * The interface method <code>SCIPincludeReaderCflp()</code> is called to add the reader plugin to \SCIP (see
 * cmain.c). This means \SCIP gets informed that this reader is available for reading input files. Therefore, the
 * function <code>SCIPincludeReader()</code> is called within this method which passes all necessary information of the
 * reader to SCIP. This information includes the name of the reader, a description, and the file extension for which the
 * file reader is in charge. In our case we selected the file extension "cflp". This means that all files which have
 * this file extension are passed to our reader for parsing. Besides these information the call
 * <code>SCIPincludeReader()</code> also passes for each callback of the reader a function pointers
 * (some of them might be NULL pointers). These function
 * pointers are used by \SCIP to run the reader. For more information about all available reader callbacks we refer to
 * the \ref READER "How to add file readers" tutorial. In the remaining section
 * we restrict ourself to the callback <code>READERREAD</code> which is the only one we implemented for the cflp
 * example. All other callbacks are not required for this example.
 *
 * @section CFLP_READERREAD The READERREAD callback method
 *
 * The READERREAD callback is in charge of parsing a file and creating the problem. To see the list of arguments this
 * functions gets see the file type_reader.h in the source of \SCIP. The following arguments are of interest in our
 * case. First of all the \SCIP pointer, the file name, and the SCIP_RESULT pointer. The \SCIP pointer gives us the
 * current environment. The file name states the file which we should open and parse. Last but not least, the SCIP_RESULT
 * pointer is required to tell \SCIP if the parsing process was successfully or
 * not. Note that in type_reader.h you also find a list of allowable result values for the SCIP_RESULT pointer and the
 * <code>SCIP_RETCODE</code> which is the return value of this function.
 *
 * @subsection CFLP_PARSING Parsing the problem
 *
 * The file can be opened and parsed with your favorite methods. In this case we are using the functionality provided by
 * \SCIP since this has some nice side effects. We are using the function SCIPfopen() which can besides standard
 * files also handle files which are packed. To find all files related to the parsing of a file, we refer to the file pub_misc.h
 * in the source of SCIP. Parsing the data out of the file is not that hard. Please look at the code and comments
 * therein for more details.
 *
 * @subsection CFLP_CREATING Creating the problem
 *
 * After parsing the file the final task for the reader is to create the problem. In our case, we pass the collected data
 * to the \ref probdata_cflp.h "main problem data plugin". For this, we use the interface methods
 * SCIPprobdataCreate() which is provided by the
 * problem data plugin (see probdata_cflp.c). After that, the reader sets the result value for the SCIP_RESULT
 * pointer to <code>SCIP_SUCCESS</code> and returns with a proper <code>SCIP_RETCODE</code>.
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/cons_setppc.h"

#include "probdata_cflp.h"
#include "reader_cflp.h"
#include "benders_cflp.h"

/**@name Reader properties
 *
 * @{
 */

#define READER_NAME             "cflpreader"
#define READER_DESC             "file reader for cflp data format"
#define READER_EXTENSION        "txt"


#define DEFAULT_USEBENDERS       TRUE

struct SCIP_ReaderData
{
   SCIP_Bool             usebenders;         /**< should Benders' decomposition be used to solve the problem */
};

/**@} */

/** creates the reader data */
static
SCIP_RETCODE readerdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA**     readerdata          /**< pointer to store the reader data */
   )
{
   assert(scip != NULL);
   assert(readerdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, readerdata) );
   (*readerdata)->usebenders = TRUE;

   return SCIP_OKAY;
}

/** frees the reader data */
static
SCIP_RETCODE readerdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA**     readerdata          /**< pointer to store the reader data */
   )
{
   assert(scip != NULL);
   assert(readerdata != NULL);
   assert(*readerdata != NULL);

   SCIPfreeBlockMemory(scip, readerdata);

   return SCIP_OKAY;
}


/**@name Callback methods
 *
 * @{
 */

/** destructor of reader to free user data (called when SCIP is exiting)*/
static
SCIP_DECL_READERFREE(readerFreeCflp)
{
   SCIP_READERDATA* readerdata;

   assert(scip != NULL);
   assert(reader != NULL);

   readerdata = SCIPreaderGetData(reader);

   SCIP_CALL( readerdataFree(scip, &readerdata) );

   return SCIP_OKAY;
}


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadCflp)
{  /*lint --e{715}*/
   SCIP_READERDATA* readerdata;

   SCIP_FILE* file;
   SCIP_Bool error;

   char name[SCIP_MAXSTRLEN];
   char format[16];
   char buffer[SCIP_MAXSTRLEN];
   SCIP_Real fixedcost;
   SCIP_Real capacity;
   SCIP_Real* demands;
   SCIP_Real** costs;
   int nfacilities;
   int ncustomers;

   SCIP_Real demand;
   SCIP_Real cost;
   int facility;
   int customer;


   int nread;
   int lineno;
   int i;
   int j;

   readerdata = SCIPreaderGetData(reader);

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
   sprintf(name, "++ uninitialized ++");

   /* read problem name */
   if( !SCIPfeof(file) )
   {
      while( lineno < 2 )
      {
         /* get next line */
         if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
            return SCIP_READERROR;
         lineno++;
      }

      /* parse dimension line */
      sprintf(format, "%%%ds\n", SCIP_MAXSTRLEN);
      nread = sscanf(buffer, format, name);
      if( nread == 0 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         return SCIP_READERROR;
      }

      SCIPdebugMsg(scip, "problem name <%s>\n", name);
   }

   capacity = 0;
   nfacilities = 0;
   fixedcost = 0;

   /* read problem dimension */
   if( !SCIPfeof(file) )
   {
      while( lineno < 4 )
      {
         /* get next line */
         if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
            return SCIP_READERROR;
         lineno++;
      }

      /* parse dimension line */
      nread = sscanf(buffer, "%d %lf %lf\n", &nfacilities, &fixedcost, &capacity);
      if( nread < 2 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         return SCIP_READERROR;
      }

      SCIPdebugMsg(scip, "number of facilities = <%d> facility capacity = <%f>, fixed cost = <%f>\n", nfacilities, capacity, fixedcost);
   }

   /* The number of customers is the same as the number of facilities */
   ncustomers = nfacilities;
   capacity = capacity*ncustomers;


   /* allocate buffer memory for storing the demand and the transport cost temporarily */
   SCIP_CALL( SCIPallocBufferArray(scip, &demands, ncustomers) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costs, ncustomers) );

   for( i = 0; i < ncustomers; i++ )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &costs[i], ncustomers) );

      for( j = 0; j < ncustomers; j++ )
         costs[i][j] = 0.0;

      demands[i] = 0.0;
   }

   /* parse demands and costs */
   error = FALSE;

   while( !SCIPfeof(file) && !error )
   {
      /* get next line */
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         break;
      lineno++;

      if( lineno >= 6 )
      {
         /* parse the line */
         nread = sscanf(buffer, "%d %d %lf %lf\n", &facility, &customer, &cost, &demand);
         if( nread == 0 )
         {
            SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
            error = TRUE;
            break;
         }

         SCIPdebugMsg(scip, "(%d, %d) found cost <%e> and demand <%d>\n", facility, customer, cost, demand);
         demands[customer - 1] += demand;
         costs[facility - 1][customer - 1] = cost;
      }
   }

   if( !error )
   {
      if( readerdata->usebenders )
      {
         /* creating the Benders' decomposition structure */
         SCIP_CALL( SCIPincludeBendersCflp(scip, ncustomers) );
      }

      /* create a new problem in SCIP */
      SCIP_CALL( SCIPprobdataCreate(scip, name, costs, demands, capacity, fixedcost, ncustomers, nfacilities,
              ncustomers, readerdata->usebenders) );
   }

   (void)SCIPfclose(file);

   for( i = ncustomers - 1; i >= 0; i-- )
      SCIPfreeBufferArray(scip, &costs[i]);

   SCIPfreeBufferArray(scip, &costs);
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

/** includes the cflp file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderCflp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create cflp reader data */
   readerdata = NULL;

   SCIP_CALL( readerdataCreate(scip, &readerdata) );

   /* include cflp reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );
   assert(reader != NULL);

   SCIP_CALL( SCIPsetReaderFree(scip, reader, readerFreeCflp) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadCflp) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "readers/" READER_NAME "/usebenders", "Should Benders' decomposition be used to solve the problem?",
         &readerdata->usebenders, FALSE, DEFAULT_USEBENDERS, NULL, NULL) );

   return SCIP_OKAY;
}

/**@} */

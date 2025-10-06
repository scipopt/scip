/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_scflp.c
 * @brief  SCFLP reader file reader
 * @author Stephen J. Maher
 *
 * This file implements the reader/parser used to read the CAP input data and builds a SCFLP instance. For more details
 * see \ref SCFLP_READER.
 *
 * @page SCFLP_READER Parsing the input format
 *
 * In the <code>data</code> directory you find a few data files which contain each one CAP problem. These data
 * files have the following structure:
 *    - The first line specifies the number of facilities (n) and the number of customers (m).
 *    - The next n lines specifies for each facility the setup cost and propduction capacity.
 *    - The following lines state the parameters for each of the customers.
 *      - First, the demand of customer is given.
 *      - The next set of lines states the transportation cost from each facility to the customer. Each line contains 7
 *      facilities, so there are ceil(n/7) lines for the transportation costs for each customer.
 *
 * For parsing that data, we implemented a reader plugin for \SCIP. A reader has several callback methods and at least
 * one interface methods (the one including the reader into \SCIP). For our purpose we only implemented the \ref
 * READERREAD callback and the interface method which adds the reader plugin to \SCIP.
 *
 * @section SCFLP_READERINCLUDE The SCIPincludeReaderScflp() interface method
 *
 * The interface method <code>SCIPincludeReaderScflp()</code> is called to add the reader plugin to \SCIP (see
 * cmain.c). This means \SCIP gets informed that this reader is available for reading input files. Therefore, the
 * function <code>SCIPincludeReader()</code> is called within this method which passes all necessary information of the
 * reader to SCIP. This information includes the name of the reader, a description, and the file extension for which the
 * file reader is in charge. In our case we selected the file extension "cap". This means that all files which have
 * this file extension are passed to our reader for parsing. Besides these information the call
 * <code>SCIPincludeReader()</code> also passes for each callback of the reader a function pointers
 * (some of them might be NULL pointers). These function
 * pointers are used by \SCIP to run the reader. For more information about all available reader callbacks we refer to
 * the \ref READER "How to add file readers" tutorial. In the remaining section
 * we restrict ourself to the callback <code>READERREAD</code> which is the only one we implemented for the SCFLP
 * example. All other callbacks are not required for this example.
 *
 * @section SCFLP_READERREAD The READERREAD callback method
 *
 * The READERREAD callback is in charge of parsing a file and creating the problem. To see the list of arguments this
 * functions gets see the file type_reader.h in the source of \SCIP. The following arguments are of interest in our
 * case. First of all the \SCIP pointer, the file name, and the SCIP_RESULT pointer. The \SCIP pointer gives us the
 * current environment. The file name states the file which we should open and parse. Last but not least, the SCIP_RESULT
 * pointer is required to tell \SCIP if the parsing process was successfully or
 * not. Note that in type_reader.h you also find a list of allowable result values for the SCIP_RESULT pointer and the
 * <code>SCIP_RETCODE</code> which is the return value of this function.
 *
 * In the READERREAD callback, the scenarios for the stochastic program are generated. A scenario represents a different
 * realisation of demand for each customer. The scenarios are generated by sampling demands from a normal distribution
 * with a mean given by the deterministic demand and the standard deviation sampled from a uniform distribution with the
 * range [0.1*mean, 0.3*mean]. The number of scenarios can be set using a runtime parameter.
 *
 * @subsection SCFLP_PARSING Parsing the problem
 *
 * The file can be opened and parsed with your favorite methods. In this case we are using the functionality provided by
 * \SCIP since this has some nice side effects. We are using the function SCIPfopen() which can besides standard
 * files also handle files which are packed. To find all files related to the parsing of a file, we refer to the file pub_misc.h
 * in the source of SCIP. Parsing the data out of the file is not that hard. Please look at the code and comments
 * therein for more details.
 *
 * @subsection SCFLP_CREATING Creating the problem
 *
 * After parsing the file the final task for the reader is to create the problem. In our case, we pass the collected data
 * to the \ref probdata_scflp.h "main problem data plugin". For this, we use the interface methods
 * SCIPprobdataCreate() which is provided by the
 * problem data plugin (see probdata_scflp.c). After that, the reader sets the result value for the SCIP_RESULT
 * pointer to <code>SCIP_SUCCESS</code> and returns with a proper <code>SCIP_RETCODE</code>.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "probdata_scflp.h"
#include "reader_scflp.h"

/**@name Reader properties
 *
 * @{
 */

#define READER_NAME             "scflpreader"
#define READER_DESC             "file reader for CAP data format and creates a SCFLP instance"
#define READER_EXTENSION        "cap"


#define DEFAULT_USEBENDERS       TRUE
#define DEFAULT_NUMSCENARIOS      250
#define DEFAULT_RANDOMSEED          1
#define DEFAULT_QUADCOSTS       FALSE   /**< should the problem be formulated with quadratic costs */


#define MAXNCOSTS                   7

struct SCIP_ReaderData
{
   SCIP_Bool             usebenders;         /**< should Benders' decomposition be used to solve the problem */
   int                   nscenarios;         /**< the number of scenarios */
   int                   randomseed;         /**< the random seed used to generate the scenarios */
   SCIP_Bool             quadcosts;          /**< should the problem be formulated with quadratic costs */
};

/**@} */

/** generates a normally distributed random number */
static
SCIP_Real generateGaussianNoise(
   SCIP_RANDNUMGEN*      randomgen,          /**< the random number generator */
   SCIP_Real             mean,               /**< the mean of the normal distribution */
   SCIP_Real             stdDev,             /**< the standard deviation of the normal distribution */
   SCIP_Real*            spare,              /**< the spare value that has been collected */
   SCIP_Bool*            hasspare            /**< whether a spare value exists */
   )
{
   SCIP_Real u;
   SCIP_Real v;
   SCIP_Real s;

   /* if a spare exists, then it is returned */
   if( (*hasspare) )
   {
      (*hasspare) = FALSE;
      return mean + stdDev * (*spare);
   }

   (*hasspare) = TRUE;
   do {
      u = SCIPrandomGetReal(randomgen, 0.0, 1.0);
      v = SCIPrandomGetReal(randomgen, 0.0, 1.0);
      s = u * u + v * v;
   }
   while( (s >= 1.0) || (s == 0.0) );

   s = sqrt(-2.0 * log(s) / s);
   (*spare) = v * s;

   return mean + stdDev * u * s;
}


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
SCIP_DECL_READERFREE(readerFreeScflp)
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
SCIP_DECL_READERREAD(readerReadScflp)
{  /*lint --e{715}*/
   SCIP_READERDATA* readerdata;

   SCIP_RANDNUMGEN* randomgen;

   SCIP_FILE* file;
   SCIP_Bool error;

   char name[SCIP_MAXSTRLEN];
   char buffer[SCIP_MAXSTRLEN];
   SCIP_Real* fixedcost;
   SCIP_Real* capacity;
   SCIP_Real* detdemands;
   SCIP_Real** demands;
   SCIP_Real** costs;
   int nfacilities;
   int ncustomers;
   int nscenarios;

   SCIP_Real facilitycost;
   SCIP_Real facilitycap;
   SCIP_Real demand;
   int facility;
   int customer;
   int ncostlines;   /* this is the number of lines that have facility costs */

   int nread;
   int lineno;
   int i;
   int j;

   readerdata = SCIPreaderGetData(reader);

   /* creating the random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &randomgen, (unsigned int) readerdata->randomseed, FALSE) );

   nscenarios = readerdata->nscenarios;

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
   (void)sprintf(name, "SCFLP");

   nfacilities = 0;
   ncustomers = 0;

   /* read problem name */
   if( !SCIPfeof(file) )
   {
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         return SCIP_READERROR;

      /* parse dimension line */
      nread = sscanf(buffer, " %d %d\n", &nfacilities, &ncustomers);
      if( nread == 0 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         return SCIP_READERROR;
      }

      SCIPdebugMsg(scip, "number of facilities = <%d> number of customers = <%d>\n", nfacilities, ncustomers);
   }

   /* computing the number of customer cost lines */
   ncostlines = (int) SCIPceil(scip, (SCIP_Real)nfacilities/MAXNCOSTS);

   /* allocate buffer memory for storing the demand and the transport cost temporarily */
   SCIP_CALL( SCIPallocBufferArray(scip, &capacity, nfacilities) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedcost, nfacilities) );

   for( i = 0; i < nfacilities; i++ )
   {
      capacity[i] = 0.0;
      fixedcost[i] = 0.0;
   }

   error = FALSE;

   /* read facility data */
   facility = 0;
   while( !SCIPfeof(file) && !error && facility < nfacilities )
   {
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         break;

      /* parse dimension line */
      nread = sscanf(buffer, " %lf %lf\n", &facilitycap, &facilitycost);
      if( nread < 2 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         error = TRUE;
         break;
      }

      SCIPdebugMsg(scip, "facility %d capacity = <%f>, fixed cost = <%f>\n", facility, facilitycap, facilitycost);
      capacity[facility] = facilitycap;
      fixedcost[facility] = facilitycost;
      facility++;
      assert(facility <= nfacilities);
   }

   /* allocate buffer memory for storing the demand and the transport cost temporarily */
   SCIP_CALL( SCIPallocBufferArray(scip, &detdemands, ncustomers) );
   SCIP_CALL( SCIPallocBufferArray(scip, &demands, ncustomers) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costs, ncustomers) );

   /* TODO: convert the 2D arrays into contiguous arrays. This has benefits from a memory management point of view. */
   for( i = 0; i < ncustomers; i++ )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &costs[i], nfacilities) );
      SCIP_CALL( SCIPallocBufferArray(scip, &demands[i], nscenarios) );

      for( j = 0; j < nfacilities; j++ )
         costs[i][j] = 0.0;

      for( j = 0; j < nscenarios; j++ )
         demands[i][j] = 0.0;

      detdemands[i] = 0.0;
   }

   /* parse demands and costs */

   lineno = 0;
   customer = 0;
   facility = 0;
   while( !SCIPfeof(file) && !error )
   {
      /* get next line */
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         break;

      if( lineno == 0 )
      {
         /* parse the line */
         nread = sscanf(buffer, " %lf\n", &demand);
         if( nread == 0 )
         {
            SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
            error = TRUE;
            break;
         }

         SCIPdebugMsg(scip, "customer %d demand <%f>\n", customer, demand);
         detdemands[customer] = demand;
      }
      else if( lineno <= ncostlines )
      {
         SCIP_Real tmpcosts[MAXNCOSTS];
         /* parse the line */
         nread = sscanf(buffer, " %lf %lf %lf %lf %lf %lf %lf\n", &tmpcosts[0], &tmpcosts[1], &tmpcosts[2],
            &tmpcosts[3], &tmpcosts[4], &tmpcosts[5], &tmpcosts[6]);

         if( nread == 0 )
         {
            SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
            error = TRUE;
            break;
         }

         for( i = 0; i < nread; i++ )
         {
            SCIPdebugMsg(scip, "(%d, %d) found cost <%e>\n", customer, facility, tmpcosts[i]);
            costs[customer][facility] = tmpcosts[i] / demand; /*lint !e530*/
            facility++;
         }
      }

      if( lineno == ncostlines )
      {
         lineno = 0;
         facility = 0;
         customer++;
      }
      else
         lineno++;
   }

   /* if there has been no error in reading the data, then we generate the scenarios */
   if( !error )
   {
      SCIP_Real mean;
      SCIP_Real stdDev;
      SCIP_Real spare;
      SCIP_Bool hasspare;

      /* creating the scenario data */
      for( i = 0; i < ncustomers; i++ )
      {
         mean = detdemands[i];
         stdDev = SCIPrandomGetReal(randomgen, 0.1*mean, 0.2*mean);
         hasspare = FALSE;
         for( j = 0; j < nscenarios; j++ )
         {
            demands[i][j] = SCIPround(scip, generateGaussianNoise(randomgen, mean, stdDev, &spare, &hasspare));

            SCIPdebugMsg(scip, "Scenario %d: customer %d demand <%f>\n", j, i, demands[i][j]);
         }
      }

      /* create a new problem in SCIP */
      SCIP_CALL( SCIPprobdataCreate(scip, name, costs, demands, capacity, fixedcost, ncustomers, nfacilities,
              nscenarios, readerdata->usebenders, readerdata->quadcosts) );
   }

   (void)SCIPfclose(file);

   for( i = ncustomers - 1; i >= 0; i-- )
   {
      SCIPfreeBufferArray(scip, &demands[i]);
      SCIPfreeBufferArray(scip, &costs[i]);
   }

   SCIPfreeBufferArray(scip, &costs);
   SCIPfreeBufferArray(scip, &demands);
   SCIPfreeBufferArray(scip, &detdemands);

   SCIPfreeBufferArray(scip, &fixedcost);
   SCIPfreeBufferArray(scip, &capacity);

   /* freeing the random number generator */
   SCIPfreeRandom(scip, &randomgen);

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

/** includes the scflp file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderScflp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create scflp reader data */
   readerdata = NULL;

   SCIP_CALL( readerdataCreate(scip, &readerdata) );

   /* include scflp reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );
   assert(reader != NULL);

   SCIP_CALL( SCIPsetReaderFree(scip, reader, readerFreeScflp) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadScflp) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/" READER_NAME "/usebenders", "Should Benders' decomposition be used to solve the problem?",
         &readerdata->usebenders, FALSE, DEFAULT_USEBENDERS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "reading/" READER_NAME "/numscenarios",
         "the number of scenarios that will be randomly generated",
         &readerdata->nscenarios, FALSE, DEFAULT_NUMSCENARIOS, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "reading/" READER_NAME "/randomseed",
         "the random seed used to generate the scenarios",
         &readerdata->randomseed, FALSE, DEFAULT_RANDOMSEED, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/" READER_NAME "/quadcosts", "should the problem be formulated with quadratic costs",
         &readerdata->quadcosts, FALSE, DEFAULT_QUADCOSTS, NULL, NULL) );

   return SCIP_OKAY;
}

/**@} */

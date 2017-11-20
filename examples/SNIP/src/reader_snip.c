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

/**@file   reader_snip.c
 * @brief  Snip problem reader file reader
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * This file implements the reader/parser used to read the snip input data. For more details see \ref SNIP_READER.
 *
 * @page SNIP_READER Parsing the input format and creating the problem
 *
 * In the <code>data</code> directory you find a few data files which contain each one snip problem. These data
 * files have the following structure. In the first line the name of the instance is stated. In the second line you find
 * three integer numbers. The first one gives you the snipacity \f$\kappa\f$, the second the number of items, and the
 * last integer states the value of a known feasible solution. This means an upper bound on the number of needed
 * bins. The remaining lines give the size for each item.
 *
 * For parsing that data, we implemented a reader plugin for \SCIP. A reader has several callback methods and at least
 * one interface methods (the one including the reader into \SCIP). For our purpose we only implemented the \ref
 * READERREAD "READERREAD" callback and the interface method which adds the reader plugin to \SCIP.
 *
 * @section SNIP_READERINCLUDE The SCIPincludeReaderSnip() interface method
 *
 * The interface method <code>SCIPincludeReaderSnip()</code> is called to add the reader plugin to \SCIP (see
 * cmain.c). This means \SCIP gets informed that this reader is available for reading input files. Therefore, the
 * function <code>SCIPincludeReader()</code> is called within this method which passes all necessary information of the
 * reader to SCIP. This information includes the name of the reader, a description, and the file extension for which the
 * file reader is in charge. In our case we selected the file extension "snip". This means that all files which have
 * this file extension are passed to our reader for parsing. Besides these information the call
 * <code>SCIPincludeReader()</code> also passes for each callback of the reader a function pointers
 * (some of them might be NULL pointers). These function
 * pointers are used by \SCIP to run the reader. For more information about all available reader callbacks we refer to
 * the \ref READER "How to add file readers" tutorial. In the remaining section
 * we restrict ourself to the callback <code>READERREAD</code> which is the only one we implemented for the snip
 * example. All other callbacks are not required for this example.
 *
 * @section SNIP_READERREAD The READERREAD callback method
 *
 * The READERREAD callback is in charge of parsing a file and creating the problem. To see the list of arguments this
 * functions gets see the file type_reader.h in the source of \SCIP. The following arguments are of interest in our
 * case. First of all the \SCIP pointer, the file name, and the SCIP_RESULT pointer. The \SCIP pointer gives us the
 * current environment. The file name states the file which we should open and parse. Last but not least, the SCIP_RESULT
 * pointer is required to tell \SCIP if the parsing process was successfully or
 * not. Note that in type_reader.h you also find a list of allowable result values for the SCIP_RESULT pointer and the
 * <code>SCIP_RETCODE</code> which is the return value of this function.
 *
 * @subsection SNIP_PARSING Parsing the problem
 *
 * The file can be opened and parsed with your favorite methods. In this case we are using the functionality provided by
 * \SCIP since this has some nice side effects. We are using the function SCIPfopen() which can besides standard
 * files also handle files which are packed. To find all files related to the parsing of a file, we refer to the file pub_misc.h
 * in the source of SCIP. Parsing the data out of the file is not that hard. Please look at the code and comments
 * therein for more details.
 *
 * @subsection SNIP_CREATING Creating the problem
 *
 * After parsing the file the final task for the reader is to create the problem. In our case, we pass the collected data
 * to the \ref probdata_snip.h "main problem data plugin". For this, we use the interface methods
 * SCIPprobdataCreate() which is provided by the
 * problem data plugin (see probdata_snip.c). After that, the reader sets the result value for the SCIP_RESULT
 * pointer to <code>SCIP_SUCCESS</code> and returns with a proper <code>SCIP_RETCODE</code>.
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "probdata_snip.h"
#include "reader_snip.h"

/**@name Reader properties
 *
 * @{
 */

#define READER_NAME             "snipreader"
#define READER_DESC             "file reader for snip data format"
#define READER_EXTENSION        "snip"


#define DEFAULT_USEBENDERS      TRUE

struct SCIP_ReaderData
{
   SCIP_Bool             usebenders;         /**< should Benders' decomposition be used to solve the problem */
};

/**@} */

/** reads the gain for traversing an arc that does not have a sensor */
static
SCIP_RETCODE readArcGain(
   SCIP*                 scip,               /**< the SCIP data structure */
   const char*           filename,           /**< the file name for the arcgain data */
   SCIP_Real*            probwosensor,       /**< array to store the probability of detection without a sensor */
   int*                  arcids,             /**< array to store the arc ids */
   SCIP_Bool*            isnode,             /**< flag to indicate whether an index is a node */
   int                   narcswosensor       /**< the number of arcs without sensors */
   )
{
   SCIP_FILE* file;
   SCIP_Bool error;
   int head;
   int tail;
   SCIP_Real cost;
   int nread;
   int i;

   char buffer[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(probwosensor != NULL);

   /* open file */
   file = SCIPfopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   error = FALSE;

   /* read facility data */
   i = 0;
   while( !SCIPfeof(file) && !error && i < narcswosensor )
   {
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         break;

      /* parse dimension line */
      nread = sscanf(buffer, "%d %d %lf\n", &tail, &head, &cost);
      if( nread < 3 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", i, filename, buffer);
         error = TRUE;
         break;
      }

      SCIPdebugMsg(scip, "ARCGAIN - arc %d: (%d, %d) cost = <%f>\n", i, tail, head, cost);
      probwosensor[i] = cost;
      arcids[i] = tail*1000 + head;
      i++;
      assert(i <= narcswosensor);

      isnode[tail] = TRUE;
      isnode[head] = TRUE;
   }

   return SCIP_OKAY;
}

/** reads the gain from traversing an arc that may contain a sensor */
static
SCIP_RETCODE readInterdictionArc(
   SCIP*                 scip,               /**< the SCIP data structure */
   const char*           filename,           /**< the file name for the arcgain data */
   SCIP_Real*            intdictwosensor,    /**< array to store the probability of detection without a sensor */
   SCIP_Real*            intdictwsensor,     /**< array to store the probability of detection with a sensor */
   int*                  intdictarcids,      /**< array to store the arc ids */
   SCIP_Bool*            isnode,             /**< flag to indicate whether an index is a node */
   int                   narcswsensor        /**< the number of arcs without sensors */
   )
{
   SCIP_FILE* file;
   SCIP_Bool error;
   int head;
   int tail;
   SCIP_Real wosensor;
   SCIP_Real wsensor;
   int nread;
   int i;

   char buffer[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(intdictwosensor != NULL);
   assert(intdictwsensor != NULL);

   /* open file */
   file = SCIPfopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   error = FALSE;

   /* read facility data */
   i = 0;
   while( !SCIPfeof(file) && !error && i < narcswsensor )
   {
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         break;

      /* parse dimension line */
      nread = sscanf(buffer, "%d %d %lf %lf\n", &tail, &head, &wosensor, &wsensor);
      if( nread < 4 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", i, filename, buffer);
         error = TRUE;
         break;
      }

      SCIPdebugMsg(scip, "INTERDICTIONARCS - arc %d: (%d, %d) without sensor = <%f> with sensor = <%f>\n", i, tail, head, wosensor,
         wsensor);
      intdictwosensor[i] = wosensor;
      intdictwsensor[i] = wsensor;
      intdictarcids[i] = tail*1000 + head;
      i++;
      assert(i <= narcswsensor);

      isnode[tail] = TRUE;
      isnode[head] = TRUE;
   }

   return SCIP_OKAY;
}

/** reads the scenarios */
static
SCIP_RETCODE readScenarios(
   SCIP*                 scip,               /**< the SCIP data structure */
   const char*           filename,           /**< the file name for the arcgain data */
   SCIP_Real*            scenariocost,       /**< array to store the costs for the scenarios */
   int*                  scenarioarcids,     /**< array to store the arc ids */
   SCIP_Bool*            isnode,             /**< flag to indicate whether an index is a node */
   int                   nscenarios          /**< the number of scenarios */
   )
{
   SCIP_FILE* file;
   SCIP_Bool error;
   int source;
   int sink;
   SCIP_Real cost;
   int nread;
   int i;

   char buffer[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(scenariocost != NULL);
   assert(scenarioarcids != NULL);

   /* open file */
   file = SCIPfopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   error = FALSE;

   /* read facility data */
   i = 0;
   while( !SCIPfeof(file) && !error && i < nscenarios )
   {
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         break;

      /* parse dimension line */
      nread = sscanf(buffer, "%d %d %lf\n", &source, &sink, &cost);
      if( nread < 3 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", i, filename, buffer);
         error = TRUE;
         break;
      }

      SCIPdebugMsg(scip, "SCENARIOS - scenario %d: source = <%d>, sink = <%d>, cost = <%f>\n", i, source, sink, cost);
      scenariocost[i] = cost;
      scenarioarcids[i] = source*1000 + sink;
      i++;
      assert(i <= nscenarios);

      isnode[source] = TRUE;
      isnode[sink] = TRUE;
   }

   return SCIP_OKAY;
}

/** reads the shortest path information */
static
SCIP_RETCODE readShortestPaths(
   SCIP*                 scip,               /**< the SCIP data structure */
   const char*           filename,           /**< the file name for the arcgain data */
   SCIP_Real**           shortestpaths,      /**< array to store the shortest paths for each scenario */
   int                   nscenarios,         /**< the number of scenarios */
   int                   nnodes              /**< the number of nodes */
   )
{
   SCIP_FILE* file;
   SCIP_Bool error;
   SCIP_Real cost;
   int scenario;
   int node;

   char buffer[10000];
   char* splitbuf;

   assert(scip != NULL);
   assert(shortestpaths != NULL);

   /* open file */
   file = SCIPfopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   error = FALSE;

   /* read facility data */
   scenario = 0;
   while( !SCIPfeof(file) && !error && scenario < nscenarios )
   {
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         break;

      /* parse dimension line */
      node = 0;
      splitbuf = strtok(buffer, " ");
      while( splitbuf != NULL && node < nnodes )
      {
         cost = atof(splitbuf);
         SCIPdebugMsg(scip, "SHORTESTPATHS - scenario %d: source = <%d>, cost = <%lf>\n", scenario, node, cost);
         shortestpaths[scenario][node] = cost;
         node++;
         splitbuf = strtok(NULL, " ");
         assert(node <= nnodes);
      }

      scenario++;
      assert(scenario <= nscenarios);
   }

   return SCIP_OKAY;
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

   SCIP_CALL( SCIPallocBlockMemory(scip, &(*readerdata)) );
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

   SCIPfreeBlockMemory(scip, &(*readerdata));

   return SCIP_OKAY;
}


/**@name Callback methods
 *
 * @{
 */

/** destructor of reader to free user data (called when SCIP is exiting)*/
static
SCIP_DECL_READERFREE(readerFreeSnip)
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
SCIP_DECL_READERREAD(readerReadSnip)
{  /*lint --e{715}*/
   SCIP_READERDATA* readerdata;

   SCIP_FILE* file;
   SCIP_Bool error;

   char newfilename[SCIP_MAXSTRLEN];
   char* fromlastslash;
   char* parent;

   char name[SCIP_MAXSTRLEN];
   char buffer[SCIP_MAXSTRLEN];
   char format[16];

   SCIP_Real* probwosensor;
   SCIP_Real* intdictwosensor;
   SCIP_Real* intdictwsensor;
   SCIP_Real* scenariocost;
   SCIP_Real** shortestpaths;
   int* arcids;
   int* intdictarcids;
   int* scenarioarcids;
   int* nodemapping;
   SCIP_Bool* isnode;
   int narcs;
   int nnodes;
   int nsensors;
   int nscenarios;
   int budget;
   int snipnumber;
   int maxnodenum;


   int nread;
   int lineno;
   int i;

   readerdata = SCIPreaderGetData(reader);

   *result = SCIP_DIDNOTRUN;
   error = FALSE;

   fromlastslash = strrchr(filename, '/');

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &parent, filename, strlen(filename) - (strlen(fromlastslash) - 1)) );
   /* open file */
   file = SCIPfopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   lineno = 0;
   sprintf(name, "SNIP");

   /* read problem name */
   if( !SCIPfeof(file) )
   {
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         return SCIP_READERROR;

      /* parse dimension line */
      nread = sscanf(buffer, "%d %d %d %d %d %d %d\n", &snipnumber, &nnodes, &narcs, &nsensors, &nscenarios, &budget,
         &maxnodenum);
      if( nread < 5 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         return SCIP_READERROR;
      }

      SCIPdebugMsg(scip, "nodes = <%d>, arcs = <%d>, sensors = <%d>, scenarios = <%d>, budget = <%d>\n",
         nnodes, narcs, nsensors, nscenarios, budget);
   }

   /* because the maxnodenum is a label and not an index, it must be increased by 1 */
   maxnodenum = maxnodenum + 1;

   /* allocate buffer memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &probwosensor, narcs - nsensors) );
   SCIP_CALL( SCIPallocBufferArray(scip, &intdictwosensor, narcs - nsensors) );
   SCIP_CALL( SCIPallocBufferArray(scip, &intdictwsensor, nsensors) );
   SCIP_CALL( SCIPallocBufferArray(scip, &scenariocost, nscenarios) );
   SCIP_CALL( SCIPallocBufferArray(scip, &shortestpaths, nscenarios) );
   for( i = 0; i < nscenarios; i++ )
      SCIP_CALL( SCIPallocBufferArray(scip, &shortestpaths[i], nnodes) );

   SCIP_CALL( SCIPallocBufferArray(scip, &arcids, narcs - nsensors) );
   SCIP_CALL( SCIPallocBufferArray(scip, &intdictarcids, nsensors) );
   SCIP_CALL( SCIPallocBufferArray(scip, &scenarioarcids, nscenarios) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodemapping, maxnodenum) );
   SCIP_CALL( SCIPallocBufferArray(scip, &isnode, maxnodenum) );

   /* initialising the isnode array */
   for( i = 0; i < maxnodenum; i++ )
      isnode[i] = FALSE;

   /* read arc gain data */
   if( !SCIPfeof(file) )
   {
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         return SCIP_READERROR;

      /* parse dimension line */
      sprintf(format, "%%%ds\n", SCIP_MAXSTRLEN);
      nread = sscanf(buffer, format, name);
      if( nread == 0 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         return SCIP_READERROR;
      }

      (void) SCIPsnprintf(newfilename, SCIP_MAXSTRLEN, "%s%s", parent, name);
      SCIPdebugMsg(scip, "reading data from <%s>\n", newfilename);
      SCIP_CALL( readArcGain(scip, newfilename, probwosensor, arcids, isnode, narcs - nsensors) );
   }

   /* read interdiction arc data */
   if( !SCIPfeof(file) )
   {
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         return SCIP_READERROR;

      /* parse dimension line */
      sprintf(format, "%%%ds\n", SCIP_MAXSTRLEN);
      nread = sscanf(buffer, format, name);
      if( nread == 0 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         return SCIP_READERROR;
      }

      (void) SCIPsnprintf(newfilename, SCIP_MAXSTRLEN, "%s%s", parent, name);
      SCIPdebugMsg(scip, "reading data from <%s>\n", newfilename);
      SCIP_CALL( readInterdictionArc(scip, newfilename, intdictwosensor, intdictwsensor, intdictarcids, isnode, nsensors) );
   }

   /* read scenario data */
   if( !SCIPfeof(file) )
   {
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         return SCIP_READERROR;

      /* parse dimension line */
      sprintf(format, "%%%ds\n", SCIP_MAXSTRLEN);
      nread = sscanf(buffer, format, name);
      if( nread == 0 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         return SCIP_READERROR;
      }

      (void) SCIPsnprintf(newfilename, SCIP_MAXSTRLEN, "%s%s", parent, name);
      SCIPdebugMsg(scip, "reading data from <%s>\n", newfilename);
      SCIP_CALL( readScenarios(scip, newfilename, scenariocost, scenarioarcids, isnode, nscenarios) );
   }

   /* read shortest paths data */
   if( !SCIPfeof(file) )
   {
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         return SCIP_READERROR;

      /* parse dimension line */
      sprintf(format, "%%%ds\n", SCIP_MAXSTRLEN);
      nread = sscanf(buffer, format, name);
      if( nread == 0 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         return SCIP_READERROR;
      }

      (void) SCIPsnprintf(newfilename, SCIP_MAXSTRLEN, "%s%s", parent, name);
      SCIPdebugMsg(scip, "reading data from <%s>\n", newfilename);
      SCIP_CALL( readShortestPaths(scip, newfilename, shortestpaths, nscenarios, nnodes) );
   }



   if( !error )
   {
      SCIP_Real multiplier;
      int nodecount;

      if( snipnumber == 3 )
         multiplier = 0.1;
      else if( snipnumber == 4 )
         multiplier = 0;
      else
      {
         SCIPwarningMessage(scip, "invalid input in file <%s>: SNIP number = <%d>. SNIP number must be 3 or 4\n",
            filename, snipnumber);
         return SCIP_READERROR;
      }

      nodecount = 0;
      for( i = 0; i < maxnodenum; i++ )
      {
         if( isnode[i] )
         {
            nodemapping[i] = nodecount;
            nodecount++;
         }
      }
      assert(nodecount == nnodes);

      SCIPinfoMessage(scip, NULL, "SNIP settings: number = <%d>, budget = <%d>\n", snipnumber, budget);
      SCIPinfoMessage(scip, NULL, "==========================================\n\n");

      /* create a new problem in SCIP */
      SCIP_CALL( SCIPprobdataCreate(scip, name, scenariocost, probwosensor, intdictwosensor, shortestpaths,
            scenarioarcids, arcids, intdictarcids, nodemapping, budget, multiplier, narcs, nnodes, nsensors, nscenarios,
            readerdata->usebenders) );
   }

   (void)SCIPfclose(file);

   /* free buffer memory */
   SCIPfreeBufferArray(scip, &nodemapping);
   SCIPfreeBufferArray(scip, &isnode);
   SCIPfreeBufferArray(scip, &scenarioarcids);
   SCIPfreeBufferArray(scip, &intdictarcids);
   SCIPfreeBufferArray(scip, &arcids);
   for( i = nscenarios - 1; i >= 0; i-- )
      SCIPfreeBufferArray(scip, &shortestpaths[i]);
   SCIPfreeBufferArray(scip, &shortestpaths);
   SCIPfreeBufferArray(scip, &scenariocost);
   SCIPfreeBufferArray(scip, &intdictwsensor);
   SCIPfreeBufferArray(scip, &intdictwosensor);
   SCIPfreeBufferArray(scip, &probwosensor);

   SCIPfreeBlockMemoryArray(scip, &parent, strlen(parent));

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

/** includes the snip file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderSnip(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create snip reader data */
   readerdata = NULL;

   SCIP_CALL( readerdataCreate(scip, &readerdata) );

   /* include snip reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );
   assert(reader != NULL);

   SCIP_CALL( SCIPsetReaderFree(scip, reader, readerFreeSnip) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadSnip) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/" READER_NAME "/usebenders", "Should Benders' decomposition be used to solve the problem?",
         &readerdata->usebenders, FALSE, DEFAULT_USEBENDERS, NULL, NULL) );

   return SCIP_OKAY;
}

/**@} */

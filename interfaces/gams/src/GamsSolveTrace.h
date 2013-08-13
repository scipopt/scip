/* Copyright (C) GAMS Development and others 2011
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Author: Stefan Vigerske
 */

#ifndef GAMSBBTRACE_H_
#define GAMSBBTRACE_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GAMS_solvetrace GAMS_SOLVETRACE; /**< GAMS solve trace data structure */

/** creates a GAMS solve trace data structure and initializes trace file for writing
 * @return 0, if successful; nonzero, if failure
 */
extern
int GAMSsolvetraceCreate(
   GAMS_SOLVETRACE**     solvetrace,         /**< buffer to store pointer of GAMS solve trace data structure */
   const char*           filename,           /**< name of trace file to write */
   const char*           solverid,           /**< solver identifier string */
   int                   optfilenr,          /**< number of options file, or 0 if none */
   const char*           probname,           /**< problem name */
   double                infinity,           /**< solver value for infinity */
   int                   nodefreq,           /**< interval in number of nodes when to write N-lines to trace files, 0 to disable N-lines */
   double                timefreq            /**< interval in seconds when to write T-lines to trace files, 0 to disable T-lines */
);

/** closes trace file and frees GAMS solve trace data structure */
extern
void GAMSsolvetraceFree(
   GAMS_SOLVETRACE**     solvetrace          /**< pointer to GAMS solve trace data structure to be freed */
);

/** adds line to GAMS solve trace file */
extern
void GAMSsolvetraceAddLine(
   GAMS_SOLVETRACE*      solvetrace,         /**< GAMS solve trace data structure */
   long long             nnodes,             /**< number of enumerated nodes so far */
   double                dualbnd,            /**< current dual bound */
   double                primalbnd           /**< current primal bound */
);

/** adds end line to GAMS solve trace file */
extern
void GAMSsolvetraceAddEndLine(
   GAMS_SOLVETRACE*      solvetrace,         /**< GAMS solve trace data structure */
   long long             nnodes,             /**< number of enumerated nodes so far */
   double                dualbnd,            /**< current dual bound */
   double                primalbnd           /**< current primal bound */
);

/** set a new value for infinity */
extern
void GAMSsolvetraceSetInfinity(
   GAMS_SOLVETRACE*      solvetrace,         /**< GAMS solve trace data structure */
   double                infinity            /**< new value for infinity */
);

/** sets starttime to the current time */
extern
void GAMSsolvetraceResetStarttime(
   GAMS_SOLVETRACE*      solvetrace          /**< GAMS solve trace data structure */
);

#ifdef __cplusplus
}
#endif

#endif /* GAMSBBTRACE_H_ */

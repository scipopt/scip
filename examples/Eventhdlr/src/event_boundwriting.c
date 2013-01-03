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

/**@file   event_boundwriting.c
 * @ingroup EVENTHDLR
 * @brief  event handler for writing primal and dual bound for all open nodes
 * @author Michael Winkler
 *
 *
 *  This event handler writes to a specified file at a given frequency the dual bounds of all open nodes and the current
 *  primal bound.
 *
 *  setting "set misc boundwriting freq 1000" will lead to write every 1000 nodes the primal bound and the dual bound
 *  front of all open nodes
 *
 *  setting "set misc boundwriting filename bounds.txt" will write the bounds to the files
 *  ascending from bounds1.txt over bounds2.txt to boundsN.txt were N is the last number for writing bounds (no
 *  filename means to write to standard out)
 *
 *  setting "set misc writesubmipdualbound TRUE" will lead to resolve each open node in a subSCIP until the root in the
 *  subSCIP is solved and as a result will print this resulting dual bound
 *
 *  An output could look as follows (here writesubmipdualbound is set to TRUE):
 *
 *  PB 201
 *
 *  5913 34 192.1 193.5
 *
 *  2884 26 162.1 162.1
 *
 *  The first line above shows the Primalbound. All following lines will show first the node number, second the depth of
 *  this open node, third the dual bound of the open node, and last the dual bound of the root node in a subSCIP of the
 *  resolved node.
 *
 *  @note you can write all bounds to a single file by adding the line
 *  \code
 *  #define ONEFILE
 *  \endcode
 *  and recompiling this example.
 *
 *  @note If you want to get a better human readable format for printing, define
 *  \code
 *  #define LONGSTATS
 *  \endcode
 *  and recompile this example.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <string.h>

#include "event_boundwriting.h"


#define EVENTHDLR_NAME         "boundwriting"
#define EVENTHDLR_DESC         "event handler for writing current primalbound and global dualbound and for all open nodes the dualbound"

#define DEFAULT_FREQ                   0LL   /**< frequency for writing primal and dual bounds */
#define DEFAULT_FILENAME                ""   /**< filename to write to */
#define DEFAULT_WRITESUBMIPDUALBOUND FALSE   /**< write dualbound after solving supmip root for all open node */

/* should the printing be with more information */
//#define LONGSTATS

/* do we want to create for each printing new files or only one */
//#define ONEFILE

/*
 * Data structures
 */

/** LP reading data */
struct SCIP_EventhdlrData
{
   FILE*                  file;
   SCIP_Longint           freq;
   char*                  filename;
   SCIP_Bool              isopen;
   char                   oldfilename[SCIP_MAXSTRLEN];
   int                    filenumber;
   SCIP_Bool              writesubmipdualbound;
};


/*
 * Local methods
 */

/** initializes the reader data */
static
void initEventhdlrdata(
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   assert(eventhdlrdata != NULL);

   eventhdlrdata->filename = NULL;
   eventhdlrdata->file = NULL;
   eventhdlrdata->isopen = FALSE;
   eventhdlrdata->oldfilename[0] = '\0';
   eventhdlrdata->filenumber = 1;
}

/* apply all changes to the submip */
static
SCIP_RETCODE applyDomainChanges(
   SCIP*                 subscip,            /**< scip to apply domain changes */
   SCIP_VAR**            vars,               /**< variables in original scip instance */
   SCIP_Real*            bounds,             /**< bounds which should be applied */
   SCIP_BOUNDTYPE*       boundtypes,         /**< bound types for bounds which should be applied */
   int                   nvars,              /**< number of variables */
   SCIP_HASHMAP*         varmap              /**< hashmap for identifiing the corresponding variables in subscip */
   )
{
   SCIP_VAR* subscipvar;
   int v;

   assert(subscip != NULL);
   assert(vars != NULL);
   assert(bounds != NULL);
   assert(boundtypes != NULL);
   assert(nvars > 0);
   assert(varmap != NULL);

   /* need to do a downwards loop, because ing decisions are collect from bottom to root and if we ed twice on a variable, then the weaker bound is behind the stronger in this array */
   for( v = nvars - 1; v >= 0; --v )
   {
      subscipvar = (SCIP_VAR*) SCIPhashmapGetImage(varmap, vars[v]);
      assert(subscipvar != NULL);

      if( boundtypes[v] == SCIP_BOUNDTYPE_LOWER )
      {
         SCIP_CALL( SCIPchgVarLbGlobal(subscip, subscipvar, bounds[v]) );
      }
      else
      {
         assert(boundtypes[v] == SCIP_BOUNDTYPE_UPPER);
         SCIP_CALL( SCIPchgVarUbGlobal(subscip, subscipvar, bounds[v]) );
      }
   }

   return SCIP_OKAY;
}

/** call writing method */
static
SCIP_RETCODE writeBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< file to write to or NULL */
   SCIP_Bool             writesubmipdualbound/**< write dualbounds of submip roots for all open nodes */
   )
{
   SCIP_NODE** opennodes;
   int nopennodes;
   int n;
   int v;

   assert(scip != NULL);

   nopennodes = -1;

#ifdef LONGSTATS
   SCIPinfoMessage(scip, file, "Status after %"SCIP_LONGINT_FORMAT" processed nodes (%d open)\n", SCIPgetNNodes(scip), SCIPgetNNodesLeft(scip));

   SCIPinfoMessage(scip, file, "Primalbound: %g\n", SCIPgetPrimalbound(scip));
   SCIPinfoMessage(scip, file, "Dualbound: %g\n", SCIPgetDualbound(scip));
#else
   SCIPinfoMessage(scip, file, "PB %g\n", SCIPgetPrimalbound(scip));
#endif

   /* get all open nodes and therefor print all dualbounds */
   for( v = 2; v >= 0; --v )
   {
      SCIP_NODE* node;

      switch( v )
      {
      case 2:
         SCIP_CALL( SCIPgetChildren(scip, &opennodes, &nopennodes) );
         break;
      case 1:
         SCIP_CALL( SCIPgetSiblings(scip, &opennodes, &nopennodes) );
         break;
      case 0:
         SCIP_CALL( SCIPgetLeaves(scip, &opennodes, &nopennodes) );
         break;
      default:
	 assert(0);
	 break;
      }
      assert(nopennodes >= 0);

      /* print all node information */
      for( n = nopennodes - 1; n >= 0 && !SCIPisStopped(scip); --n )
      {
         node = opennodes[n];

         if( writesubmipdualbound )
         {
            SCIP* subscip;
            SCIP_Bool valid;
            SCIP_HASHMAP* varmap;                     /* mapping of SCIP variables to sub-SCIP variables */
            SCIP_VAR** vars;                          /* original problem's variables                    */
            int nvars;
            SCIP_Real submipdb;
	    SCIP_Bool cutoff;

            SCIP_CALL( SCIPcreate(&subscip) );

            SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

            /* create the variable mapping hash map */
            SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(subscip), SCIPcalcHashtableSize(5 * nvars)) );

            submipdb = SCIP_INVALID;
            valid = FALSE;
	    cutoff = FALSE;
            SCIP_CALL( SCIPcopy(scip, subscip, varmap, NULL, "__boundwriting", TRUE, FALSE, TRUE, &valid) );

            if( valid )
            {
               SCIP_VAR** branchvars;
               SCIP_Real* branchbounds;
               SCIP_BOUNDTYPE* boundtypes;
               int nbranchvars;
               int size;

               size = SCIPnodeGetDepth(node);

               /* allocate memory for all branching decisions */
               SCIP_CALL( SCIPallocBufferArray(scip, &branchvars, size) );
               SCIP_CALL( SCIPallocBufferArray(scip, &branchbounds, size) );
               SCIP_CALL( SCIPallocBufferArray(scip, &boundtypes, size) );

               /* we assume that we only have one branching decision at each node */
               SCIPnodeGetAncestorBranchings( node, branchvars, branchbounds, boundtypes, &nbranchvars, size );

               /* check if did not have enough memory */
               if( nbranchvars > size )
               {
                  size = nbranchvars;
                  SCIP_CALL( SCIPallocBufferArray(scip, &branchvars, size) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &branchbounds, size) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &boundtypes, size) );

                  /* now getting all information */
                  SCIPnodeGetAncestorBranchings( node, branchvars, branchbounds, boundtypes, &nbranchvars, size );
               }

               /* apply all changes to the submip */
               SCIP_CALL( applyDomainChanges(subscip, branchvars, branchbounds, boundtypes, nbranchvars, varmap) );

               /* free memory for all branching decisions */
               SCIPfreeBufferArray(scip, &boundtypes);
               SCIPfreeBufferArray(scip, &branchbounds);
               SCIPfreeBufferArray(scip, &branchvars);

	       /* do not abort subproblem on CTRL-C */
	       SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );
	       /* disable output to console */
	       SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
	       /* solve only root node */
	       SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", 1LL) );

	       /* set cutoffbound as objective limit for subscip */
	       SCIP_CALL( SCIPsetObjlimit(subscip, SCIPgetCutoffbound(scip)) );

	       SCIP_CALL( SCIPsolve(subscip) );

	       cutoff = (SCIPgetStatus(subscip) == SCIP_STATUS_INFEASIBLE);
	       submipdb = SCIPgetDualbound(subscip) * SCIPgetTransObjscale(scip) + SCIPgetTransObjoffset(scip);
	    }

#ifdef LONGSTATS
            SCIPinfoMessage(scip, file, "Node %"SCIP_LONGINT_FORMAT" (depth %d): dualbound: %g, nodesubmiprootdualbound: %g %s\n", SCIPnodeGetNumber(node), SCIPnodeGetDepth(node), SCIPgetNodeDualbound(scip, node), submipdb, cutoff ? "(cutoff)" : "");
#else
	    SCIPinfoMessage(scip, file, "%"SCIP_LONGINT_FORMAT" %d %g %g %s\n", SCIPnodeGetNumber(node), SCIPnodeGetDepth(node), SCIPgetNodeDualbound(scip, node), submipdb, cutoff ? "(cutoff)" : "");
#endif

            /* free hash map */
            SCIPhashmapFree(&varmap);

            SCIP_CALL( SCIPfree(&subscip) );
         }
         else
         {
#ifdef LONGSTATS
            SCIPinfoMessage(scip, file, "Node %"SCIP_LONGINT_FORMAT" (depth %d): dualbound: %g\n", SCIPnodeGetNumber(node), SCIPnodeGetDepth(node), SCIPgetNodeDualbound(scip, node));
#else
            SCIPinfoMessage(scip, file, "%"SCIP_LONGINT_FORMAT" %d %g\n", SCIPnodeGetNumber(node), SCIPnodeGetDepth(node), SCIPgetNodeDualbound(scip, node));
#endif
         }
      }
   }

#ifdef LONGSTATS
   SCIPinfoMessage(scip, file, "\n");
#endif

   return SCIP_OKAY;
}

/*
 * Callback methods of event handler
 */

/** copy method for event handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_EVENTCOPY(eventCopyBoundwriting)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   /* call inclusion method of event handler */
   SCIP_CALL( SCIPincludeEventHdlrBoundwriting(scip) );

   return SCIP_OKAY;
}

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeBoundwriting)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(!eventhdlrdata->isopen);
   assert(eventhdlrdata->file == NULL);

   SCIPfreeMemory(scip, &eventhdlrdata);

   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitBoundwriting)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   /* notify SCIP that your event handler wants to react on the event type best solution found */
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitBoundwriting)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   /* notify SCIP that your event handler wants to drop the event type best solution found */
   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, NULL, -1) );

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   if( eventhdlrdata->isopen )
   {
      (void) fclose(eventhdlrdata->file);
      eventhdlrdata->isopen = FALSE;
   }
   eventhdlrdata->file = NULL;
   eventhdlrdata->oldfilename[0] = '\0';
   eventhdlrdata->filenumber = 1;

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecBoundwriting)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(((SCIPeventGetType(event) & SCIP_EVENTTYPE_NODESOLVED) == SCIP_EVENTTYPE_NODEFEASIBLE) || ((SCIPeventGetType(event) & SCIP_EVENTTYPE_NODESOLVED) == SCIP_EVENTTYPE_NODEINFEASIBLE) || ((SCIPeventGetType(event) & SCIP_EVENTTYPE_NODESOLVED) == SCIP_EVENTTYPE_NODEBRANCHED));

   SCIPdebugMessage("exec method of event handler for writing primal- and dualbounds\n");

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

#ifdef ONEFILE
   /* check if we need to open the file */
   if( strlen(eventhdlrdata->filename) > 0 && !eventhdlrdata->isopen )
   {
      assert(eventhdlrdata->file == NULL);
      assert(eventhdlrdata->oldfilename[0] == '\0');

      eventhdlrdata->file = fopen(eventhdlrdata->filename, "w");
      (void)strncpy(eventhdlrdata->oldfilename, eventhdlrdata->filename, SCIP_MAXSTRLEN);

      if( eventhdlrdata->file == NULL )
      {
         SCIPerrorMessage("cannot create file <%s> for writing\n", eventhdlrdata->filename);
         SCIPprintSysError(eventhdlrdata->filename);
         return SCIP_FILECREATEERROR;
      }
      eventhdlrdata->isopen = TRUE;

#ifdef LONGSTATS
      SCIPinfoMessage(scip, eventhdlrdata->file, "Problem: %s (%d Original Constraints, %d Original Variables)\n", SCIPgetProbName(scip), SCIPgetNOrigConss(scip), SCIPgetNOrigVars(scip) );
      SCIPinfoMessage(scip, eventhdlrdata->file, "\t (%d Presolved Constraints, %d Presolved Variables, (%d binary, %d integer, %d implicit integer, %d continuous))\n", SCIPgetNConss(scip), SCIPgetNVars(scip), SCIPgetNBinVars(scip), SCIPgetNIntVars(scip), SCIPgetNImplVars(scip), SCIPgetNContVars(scip));
      SCIPinfoMessage(scip, eventhdlrdata->file, "\n");
#endif
   }
#endif

   /* call writing only at right moment */
   if( eventhdlrdata->freq == 0 || (SCIPgetNNodes(scip) % eventhdlrdata->freq) != 0 )
      return SCIP_OKAY;

#ifndef ONEFILE
   if( strlen(eventhdlrdata->filename) > 0 )
   {
      char name[SCIP_MAXSTRLEN];
      char number[SCIP_MAXSTRLEN];
      char* pch;
      int n;

      assert(eventhdlrdata->file == NULL);
      assert(!eventhdlrdata->isopen);

      if( eventhdlrdata->oldfilename[0] == '\0' )
         (void)strncpy(eventhdlrdata->oldfilename, eventhdlrdata->filename, SCIP_MAXSTRLEN);

      /* find last '.' to append filenumber */
      pch=strrchr(eventhdlrdata->filename,'.');

      assert(eventhdlrdata->filenumber > 0);
      n=sprintf(number, "%"SCIP_LONGINT_FORMAT"", eventhdlrdata->filenumber * eventhdlrdata->freq);
      assert(n > 0);
      assert(n < SCIP_MAXSTRLEN);

      /* if no point is found, extend directly */
      if( pch == NULL )
      {
         (void)strncpy(name, eventhdlrdata->filename, (unsigned int)(SCIP_MAXSTRLEN - n));
         strncat(name, number, (unsigned int)n);
      }
      else
      {
         int len;

         if( pch-(eventhdlrdata->filename) > (SCIP_MAXSTRLEN - n) )
            len = SCIP_MAXSTRLEN - n;
         else
            len = pch-(eventhdlrdata->filename);

         (void)strncpy(name, eventhdlrdata->filename, (unsigned int)len);
         name[len] = '\0';
         strncat(name, number, (unsigned int)n);
	 assert(len+n < SCIP_MAXSTRLEN);
         name[len+n] = '\0';

         if( len + n + strlen(&(eventhdlrdata->filename[len])) < SCIP_MAXSTRLEN )
         {
            strncat(name, &(eventhdlrdata->filename[len]), strlen(&(eventhdlrdata->filename[len])));
            name[strlen(eventhdlrdata->filename)+n] = '\0';
         }
      }

      eventhdlrdata->file = fopen(name, "w");

      if( eventhdlrdata->file == NULL )
      {
         SCIPerrorMessage("cannot create file <%s> for writing\n", eventhdlrdata->filename);
         SCIPprintSysError(eventhdlrdata->filename);
         return SCIP_FILECREATEERROR;
      }
      eventhdlrdata->isopen = TRUE;

#ifdef LONGSTATS
      SCIPinfoMessage(scip, eventhdlrdata->file, "Problem: %s (%d Original Constraints, %d Original Variables)\n", SCIPgetProbName(scip), SCIPgetNOrigConss(scip), SCIPgetNOrigVars(scip) );
      SCIPinfoMessage(scip, eventhdlrdata->file, "\t (%d Active Constraints, %d Active Variables, (%d binary, %d integer, %d implicit integer, %d continuous))\n", SCIPgetNConss(scip), SCIPgetNVars(scip), SCIPgetNBinVars(scip), SCIPgetNIntVars(scip), SCIPgetNImplVars(scip), SCIPgetNContVars(scip));
      SCIPinfoMessage(scip, eventhdlrdata->file, "\n");
#endif
   }
#endif

#ifndef NDEBUG
   /* check the filename did not change during the solving */
   if( strlen(eventhdlrdata->filename) > 0 && eventhdlrdata->isopen )
   {
      char tmp[SCIP_MAXSTRLEN];

      (void)strncpy(tmp, eventhdlrdata->filename, SCIP_MAXSTRLEN);

      /* the name should stay the same */
      assert(strcmp(tmp, eventhdlrdata->oldfilename) == 0);
   }
#endif

   /* call writing method */
   SCIP_CALL( writeBounds(scip, eventhdlrdata->file, eventhdlrdata->writesubmipdualbound) );

#ifndef ONEFILE
   if( strlen(eventhdlrdata->filename) > 0 )
   {
      assert(eventhdlrdata->isopen);

      (void) fclose(eventhdlrdata->file);
      eventhdlrdata->isopen = FALSE;
      eventhdlrdata->file = NULL;
      ++(eventhdlrdata->filenumber);
   }
#endif

   return SCIP_OKAY;
}

/** includes event handler for writing primal- and dualbound for all open nodes */
SCIP_RETCODE SCIPincludeEventHdlrBoundwriting(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   /* create bounds reader data */
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );
   initEventhdlrdata(eventhdlrdata);

   eventhdlr = NULL;
   /* create event handler for events on watched variables */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecBoundwriting, eventhdlrdata) );
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPsetEventhdlrCopy(scip, eventhdlr, eventCopyBoundwriting) );
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeBoundwriting) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitBoundwriting) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitBoundwriting) );

   /* add boundwriting parameters */
   SCIP_CALL( SCIPaddLongintParam(scip,
         "eventhdlr/"EVENTHDLR_NAME"/freq",
         "in which frequency should all bounds be written(0: never)",
         &eventhdlrdata->freq, FALSE, DEFAULT_FREQ, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddStringParam(scip,
         "eventhdlr/"EVENTHDLR_NAME"/filename",
         "filename to write all bounds to",
         &eventhdlrdata->filename, FALSE, DEFAULT_FILENAME, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "eventhdlr/"EVENTHDLR_NAME"/writesubmipdualbound",
         "should the dualbound of the submip-root which was created out of an open node be printed",
         &eventhdlrdata->writesubmipdualbound, FALSE, DEFAULT_WRITESUBMIPDUALBOUND, NULL, NULL) );

   return SCIP_OKAY;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip.c
 * @brief  SCIP callable library
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip.h"
#include "set.h"
#include "mem.h"
#include "prob.h"
#include "stat.h"
#include "solve.h"
#include "lpi.h"


/** SCIP main data structure */
struct Scip
{
   STAGE            stage;              /**< SCIP operation stage */
   SET*             set;                /**< global SCIP settings */
   MEM*             mem;                /**< block memory buffers */
   PROB*            origprob;           /**< original problem data */
   PROB*            transprob;          /**< transformed problem after presolve */
   STAT*            stat;               /**< dynamic problem statistics */
   TREE*            tree;               /**< branch and bound tree */
   LP*              lp;                 /**< LP data */
   PRICE*           price;              /**< storage for priced variables */
   SEPA*            sepa;               /**< storage for separated cuts */
   BRANCHCAND*      branchcand;         /**< storage for branching candidates */
   CUTPOOL*         cutpool;            /**< global cut pool */
   PRIMAL*          primal;             /**< primal data and solution storage */
   EVENTFILTER*     eventfilter;        /**< event filter for global (not variable dependent) events */
   EVENTQUEUE*      eventqueue;         /**< event queue to cache events and process them later (bound change events) */
};



/** checks, if SCIP is in one of the feasible stages */
static
RETCODE checkStage(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      method,             /**< method that was called */
   Bool             init,               /**< may method be called in the INIT stage? */
   Bool             problem,            /**< may method be called in the PROBLEM stage? */
   Bool             initsolve,          /**< may method be called in the INITSOLVE stage? */
   Bool             presolving,         /**< may method be called in the PRESOLVING stage? */
   Bool             solving,            /**< may method be called in the SOLVING stage? */
   Bool             solved,             /**< may method be called in the SOLVED stage? */
   Bool             freesolve           /**< may method be called in the FREESOLVE stage? */
   )
{
   char s[255];

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(scip->mem != NULL);
   assert(method != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      assert(scip->origprob == NULL);
      assert(scip->stat == NULL);
      assert(scip->transprob == NULL);
      assert(scip->lp == NULL);
      assert(scip->price == NULL);
      assert(scip->sepa == NULL);
      assert(scip->branchcand == NULL);
      assert(scip->cutpool == NULL);
      assert(scip->primal == NULL);
      assert(scip->eventfilter == NULL);
      assert(scip->eventqueue == NULL);
      assert(scip->tree == NULL);

      if( !init )
      {
         sprintf(s, "Cannot call method <%s> in initialization stage", method);
         errorMessage(s);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_PROBLEM:
      assert(scip->origprob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob == NULL);
      assert(scip->lp == NULL);
      assert(scip->price == NULL);
      assert(scip->sepa == NULL);
      assert(scip->branchcand == NULL);
      assert(scip->cutpool == NULL);
      assert(scip->primal == NULL);
      assert(scip->eventfilter == NULL);
      assert(scip->eventqueue == NULL);
      assert(scip->tree == NULL);

      if( !problem )
      {
         sprintf(s, "Cannot call method <%s> in problem creation stage", method);
         errorMessage(s);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
      assert(scip->origprob != NULL);
      assert(scip->stat != NULL);
      assert(scip->tree == NULL);

      if( !initsolve )
      {
         sprintf(s, "Cannot call method <%s> in solve initialization stage", method);
         errorMessage(s);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
      assert(scip->origprob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);
      assert(scip->sepa != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->cutpool != NULL);
      assert(scip->primal != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->tree == NULL);

      if( !presolving )
      {
         sprintf(s, "Cannot call method <%s> in presolving stage", method);
         errorMessage(s);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
      assert(scip->origprob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);
      assert(scip->sepa != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->cutpool != NULL);
      assert(scip->primal != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->tree != NULL);

      if( !solving )
      {
         sprintf(s, "Cannot call method <%s> in solving stage", method);
         errorMessage(s);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVED:
      assert(scip->origprob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);
      assert(scip->sepa != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->cutpool != NULL);
      assert(scip->primal != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->tree != NULL);

      if( !solved )
      {
         sprintf(s, "Cannot call method <%s> in problem solved stage", method);
         errorMessage(s);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_FREESOLVE:
      assert(scip->origprob != NULL);
      assert(scip->stat != NULL);

      if( !freesolve )
      {
         sprintf(s, "Cannot call method <%s> in solve deinitialization stage", method);
         errorMessage(s);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   default:
      errorMessage("Unknown SCIP stage");
      return SCIP_ERROR;
   }
}




/*
 * miscellaneous methods
 */

/** returns scip version number */
Real SCIPversion(
   void
   )
{
   return (Real)(SCIP_VERSION)/100.0;
}

/** prints a version information line to a file stream */
void SCIPprintVersion(
   FILE*            file                /**< file stream to write version information, or NULL for stdout */
   )
{
   if( file == NULL )
      file = stdout;

   fprintf(file, "SCIP version %g [precision: %d byte]", SCIPversion(), (int)sizeof(Real));
#ifdef SCIP_BLOCKMEMORY
   fprintf(file, " [memory: block]");
#else
   fprintf(file, " [memory: standard]");
#endif
#ifndef NDEBUG
   fprintf(file, " [mode: debug]");
#else
   fprintf(file, " [mode: optimized]");
#endif
   fprintf(file, " [LP solver: %s]\n", SCIPlpiGetSolverName());
}

/** prints error message and aborts program execution */
void SCIPerror(
   FILE*            errout,             /**< file stream to write error message */
   RETCODE          retcode             /**< SCIP return code causing the error */
   )
{
   fprintf(errout, "SCIP Error (%d): ", retcode);
   SCIPretcodePrint(errout, retcode);
   fprintf(errout, "\n");
   abort();
}




/*
 * general SCIP methods
 */

/** creates and initializes SCIP data structures */
RETCODE SCIPcreate(
   SCIP**           scip                /**< pointer to SCIP data structure */
   )
{
   assert(scip != NULL);

   ALLOC_OKAY( allocMemory(scip) );

   (*scip)->stage = SCIP_STAGE_INIT;

   CHECK_OKAY( SCIPsetCreate(&(*scip)->set, *scip) );
   CHECK_OKAY( SCIPmemCreate(&(*scip)->mem) );
   (*scip)->origprob = NULL;
   (*scip)->stat = NULL;
   (*scip)->transprob = NULL;
   (*scip)->tree = NULL;
   (*scip)->lp = NULL;
   (*scip)->price = NULL;
   (*scip)->sepa = NULL;
   (*scip)->branchcand = NULL;
   (*scip)->cutpool = NULL;
   (*scip)->primal = NULL;
   (*scip)->eventfilter = NULL;
   (*scip)->eventqueue = NULL;
   
   return SCIP_OKAY;
}

/** frees SCIP data structures */
RETCODE SCIPfree(
   SCIP**           scip                /**< pointer to SCIP data structure */
   )
{
   assert(scip != NULL);

   CHECK_OKAY( checkStage(*scip, "SCIPfree", TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   CHECK_OKAY( SCIPfreeProb(*scip) );
   assert((*scip)->stage == SCIP_STAGE_INIT);

   CHECK_OKAY( SCIPmemFree(&(*scip)->mem) );
   CHECK_OKAY( SCIPsetFree(&(*scip)->set) );

   freeMemory(scip);

   return SCIP_OKAY;
}

/** gets verbosity level for message output */
VERBLEVEL SCIPverbLevel(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->verblevel;
}

/** sets verbosity level for message output */
RETCODE SCIPsetVerbLevel(
   SCIP*            scip,               /**< SCIP data structure */
   VERBLEVEL        verblevel           /**< verbosity level for message output */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetVerbLevel", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetSetVerbLevel(scip->set, verblevel);
}

/** returns current stage of SCIP */
STAGE SCIPstage(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return scip->stage;
}




/*
 * SCIP user functionality methods: managing plugins
 */

/** creates a reader and includes it in SCIP */
RETCODE SCIPincludeReader(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of reader */
   const char*      desc,               /**< description of reader */
   const char*      extension,          /**< file extension that reader processes */
   DECL_READERFREE((*readerfree)),      /**< destructor of reader */
   DECL_READERREAD((*readerread)),      /**< read method */
   READERDATA*      readerdata          /**< reader data */
   )
{
   READER* reader;

   CHECK_OKAY( checkStage(scip, "SCIPincludeReader", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPreaderCreate(&reader, name, desc, extension, readerfree, readerread, readerdata) );
   CHECK_OKAY( SCIPsetIncludeReader(scip->set, reader) );
   
   return SCIP_OKAY;
}

/** creates a constraint handler and includes it in SCIP */
RETCODE SCIPincludeConsHdlr(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of constraint handler */
   const char*      desc,               /**< description of constraint handler */
   int              sepapriority,       /**< priority of the constraint handler for separation */
   int              enfopriority,       /**< priority of the constraint handler for constraint enforcing */
   int              chckpriority,       /**< priority of the constraint handler for checking infeasibility */
   int              propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
   Bool             needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
   DECL_CONSFREE((*consfree)),          /**< destructor of constraint handler */
   DECL_CONSINIT((*consinit)),          /**< initialise constraint handler */
   DECL_CONSEXIT((*consexit)),          /**< deinitialise constraint handler */
   DECL_CONSDELE((*consdele)),          /**< free specific constraint data */
   DECL_CONSTRAN((*constran)),          /**< transform constraint data into data belonging to the transformed problem */
   DECL_CONSSEPA((*conssepa)),          /**< separate cutting planes */
   DECL_CONSENLP((*consenlp)),          /**< enforcing constraints for LP solutions */
   DECL_CONSENPS((*consenps)),          /**< enforcing constraints for pseudo solutions */
   DECL_CONSCHCK((*conschck)),          /**< check feasibility of primal solution */
   DECL_CONSPROP((*consprop)),          /**< propagate variable domains */
   CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   CONSHDLR* conshdlr;

   CHECK_OKAY( checkStage(scip, "SCIPincludeConsHdlr", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconshdlrCreate(&conshdlr, name, desc, sepapriority, enfopriority, chckpriority, propfreq, needscons,
                  consfree, consinit, consexit, consdele, constran, conssepa, consenlp, consenps, conschck, consprop,
                  conshdlrdata) );
   CHECK_OKAY( SCIPsetIncludeConsHdlr(scip->set, conshdlr) );
   
   return SCIP_OKAY;
}

/** finds the constraint handler of the given name */
RETCODE SCIPfindConsHdlr(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of constraint handler */
   CONSHDLR**       conshdlr            /**< pointer for storing the constraint handler (returns NULL, if not found) */
   )
{
   assert(name != NULL);
   assert(conshdlr != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPfindConsHdlr", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetFindConsHdlr(scip->set, name, conshdlr) );

   return SCIP_OKAY;
}

/** creates a primal heuristic and includes it in SCIP */
RETCODE SCIPincludeHeur(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of primal heuristic */
   const char*      desc,               /**< description of primal heuristic */
   char             dispchar,           /**< display character of primal heuristic */
   int              priority,           /**< priority of the primal heuristic */
   int              freq,               /**< frequency for calling primal heuristic */
   DECL_HEURFREE((*heurfree)),          /**< destructor of primal heuristic */
   DECL_HEURINIT((*heurinit)),          /**< initialise primal heuristic */
   DECL_HEUREXIT((*heurexit)),          /**< deinitialise primal heuristic */
   DECL_HEUREXEC((*heurexec)),          /**< execution method of primal heuristic */
   HEURDATA*        heurdata            /**< primal heuristic data */
   )
{
   HEUR* heur;

   CHECK_OKAY( checkStage(scip, "SCIPincludeHeur", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPheurCreate(&heur, name, desc, dispchar, priority, freq,
                  heurfree, heurinit, heurexit, heurexec, heurdata) );
   CHECK_OKAY( SCIPsetIncludeHeur(scip->set, heur) );
   
   return SCIP_OKAY;
}

/** finds the primal heuristic of the given name */
RETCODE SCIPfindHeur(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of primal heuristic */
   HEUR**           heur                /**< pointer for storing the primal heuristic (returns NULL, if not found) */
   )
{
   assert(name != NULL);
   assert(heur != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPfindHeur", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetFindHeur(scip->set, name, heur) );

   return SCIP_OKAY;
}

/** creates an event handler and includes it in SCIP */
RETCODE SCIPincludeEventhdlr(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of event handler */
   const char*      desc,               /**< description of event handler */
   DECL_EVENTFREE((*eventfree)),        /**< destructor of event handler */
   DECL_EVENTINIT((*eventinit)),        /**< initialise event handler */
   DECL_EVENTEXIT((*eventexit)),        /**< deinitialise event handler */
   DECL_EVENTDELE((*eventdele)),        /**< free specific event data */
   DECL_EVENTEXEC((*eventexec)),        /**< execute event handler */
   EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   EVENTHDLR* eventhdlr;

   CHECK_OKAY( checkStage(scip, "SCIPincludeEventHdlr", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPeventhdlrCreate(&eventhdlr, name, desc, 
                  eventfree, eventinit, eventexit, eventdele, eventexec,
                  eventhdlrdata) );
   CHECK_OKAY( SCIPsetIncludeEventHdlr(scip->set, eventhdlr) );
   
   return SCIP_OKAY;
}

/** finds the event handler of the given name */
RETCODE SCIPfindEventHdlr(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of event handler */
   EVENTHDLR**      eventhdlr           /**< pointer for storing the event handler (returns NULL, if not found) */
   )
{
   assert(name != NULL);
   assert(eventhdlr != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPfindEventHdlr", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetFindEventHdlr(scip->set, name, eventhdlr) );

   return SCIP_OKAY;
}

/** creates a node selector and includes it in SCIP */
RETCODE SCIPincludeNodesel(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of node selector */
   const char*      desc,               /**< description of node selector */
   DECL_NODESELFREE((*nodeselfree)),    /**< destructor of node selector */
   DECL_NODESELINIT((*nodeselinit)),    /**< initialise node selector */
   DECL_NODESELEXIT((*nodeselexit)),    /**< deinitialise node selector */
   DECL_NODESELSLCT((*nodeselslct)),    /**< node selection method */
   DECL_NODESELCOMP((*nodeselcomp)),    /**< node comparison method */
   NODESELDATA*     nodeseldata,        /**< node selector data */
   Bool             lowestboundfirst    /**< does node comparison sorts w.r.t. lower bound as primal criterion? */
   )
{
   NODESEL* nodesel;

   CHECK_OKAY( checkStage(scip, "SCIPincludeNodesel", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPnodeselCreate(&nodesel, name, desc,
                  nodeselfree, nodeselinit, nodeselexit, nodeselslct, nodeselcomp, nodeseldata, lowestboundfirst) );
   CHECK_OKAY( SCIPsetIncludeNodesel(scip->set, nodesel) );
   
   return SCIP_OKAY;
}

/** creates a branching rule and includes it in SCIP */
RETCODE SCIPincludeBranchrule(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of branching rule */
   const char*      desc,               /**< description of branching rule */
   int              priority,           /**< priority of the branching rule */
   DECL_BRANCHFREE((*branchfree)),      /**< destructor of branching rule */
   DECL_BRANCHINIT((*branchinit)),      /**< initialise branching rule */
   DECL_BRANCHEXIT((*branchexit)),      /**< deinitialise branching rule */
   DECL_BRANCHEXLP((*branchexlp)),      /**< branching execution method for fractional LP solutions */
   DECL_BRANCHEXPS((*branchexps)),      /**< branching execution method for not completely fixed pseudo solutions */
   BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   )
{
   BRANCHRULE* branchrule;

   CHECK_OKAY( checkStage(scip, "SCIPincludeBranchrule", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPbranchruleCreate(&branchrule, name, desc, priority,
                  branchfree, branchinit, branchexit, branchexlp, branchexps, branchruledata) );
   CHECK_OKAY( SCIPsetIncludeBranchrule(scip->set, branchrule) );
   
   return SCIP_OKAY;
}

/** creates a display column and includes it in SCIP */
RETCODE SCIPincludeDisp(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of display column */
   const char*      desc,               /**< description of display column */
   const char*      header,             /**< head line of display column */
   DECL_DISPFREE((*dispfree)),          /**< destructor of display column */
   DECL_DISPINIT((*dispinit)),          /**< initialise display column */
   DECL_DISPEXIT((*dispexit)),          /**< deinitialise display column */
   DECL_DISPOUTP((*dispoutp)),          /**< output method */
   DISPDATA*        dispdata,           /**< display column data */
   int              width,              /**< width of display column (no. of chars used) */
   int              priority,           /**< priority of display column */
   int              position,           /**< relative position of display column */
   Bool             stripline           /**< should the column be separated with a line from its right neighbour? */
   )
{
   DISP* disp;

   CHECK_OKAY( checkStage(scip, "SCIPincludeDisp", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPdispCreate(&disp, name, desc, header, dispfree, dispinit, dispexit, dispoutp, dispdata,
                  width, priority, position, stripline) );
   CHECK_OKAY( SCIPsetIncludeDisp(scip->set, disp) );
   
   return SCIP_OKAY;
}




/*
 * global problem methods
 */

/** creates empty problem and initializes all solving data structures */
RETCODE SCIPcreateProb(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< problem name */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateProb", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   scip->stage = SCIP_STAGE_PROBLEM;
   
   CHECK_OKAY( SCIPstatCreate(&scip->stat) );
   CHECK_OKAY( SCIPprobCreate(&scip->origprob, name) );
   
   return SCIP_OKAY;
}

/** reads problem from file and initializes all solving data structures */
RETCODE SCIPreadProb(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      filename            /**< problem file name */
   )
{
   RESULT result;
   int i;
   char s[255];

   assert(filename != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPreadProb", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* try all readers until one could read the file */
   result = SCIP_DIDNOTRUN;
   for( i = 0; i < scip->set->nreaders && result == SCIP_DIDNOTRUN; ++i )
   {
      CHECK_OKAY( SCIPreaderRead(scip->set->readers[i], scip, filename, &result) );
   }

   switch( result )
   {
   case SCIP_DIDNOTRUN:
      failureMessage("No reader for input file <%s> available\n", filename);
      return SCIP_READERR;

   case SCIP_SUCCESS:
      assert(scip->origprob != NULL);
      
      sprintf(s, "original problem has %d variables (%d bin, %d int, %d impl, %d cont) and %d constraints",
         scip->origprob->nvars, scip->origprob->nbin, scip->origprob->nint, scip->origprob->nimpl, scip->origprob->ncont,
         scip->origprob->nconss);
      infoMessage(SCIPverbLevel(scip), SCIP_VERBLEVEL_HIGH, s);
#if 0
      printf(" var names :  ");
      SCIPhashtablePrintStatistics(scip->origprob->varnames);
      printf(" cons names:  ");
      SCIPhashtablePrintStatistics(scip->origprob->consnames);
#endif

      return SCIP_OKAY;

   default:
      assert(i < scip->set->nreaders);
      sprintf(s, "invalid result code <%d> from reader <%s> reading file <%s>", 
         result, SCIPreaderGetName(scip->set->readers[i]), filename);
      errorMessage(s);
      return SCIP_READERR;
   }
}

/** frees problem and solution process data */
RETCODE SCIPfreeProb(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeProb", TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   CHECK_OKAY( SCIPfreeSolve(scip) );
   assert(scip->stage == SCIP_STAGE_INIT || scip->stage == SCIP_STAGE_PROBLEM);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      return SCIP_OKAY;

   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPprobFree(&scip->origprob, scip->mem->probmem, scip->set, scip->lp) );
      CHECK_OKAY( SCIPstatFree(&scip->stat) );

      scip->stage = SCIP_STAGE_INIT;

      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

/** sets objective sense of problem */
RETCODE SCIPsetObjsense(
   SCIP*            scip,               /**< SCIP data structure */
   OBJSENSE         objsense            /**< new objective sense */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetObjsense", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   if( objsense != SCIP_OBJSENSE_MAXIMIZE && objsense != SCIP_OBJSENSE_MINIMIZE )
   {
      errorMessage("Invalid objective sense");
      return SCIP_INVALIDDATA;
   }

   SCIPprobSetObjsense(scip->origprob, objsense);
   
   return SCIP_OKAY;
}
   
/** adds variable to the problem */
RETCODE SCIPaddVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to add */
   )
{
   assert(var != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddVar", FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      if( var->varstatus != SCIP_VARSTATUS_ORIGINAL )
      {
         errorMessage("Cannot add transformed variables to original problem");
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPprobAddVar(scip->origprob, scip->mem->probmem, scip->set, var) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
      if( var->varstatus == SCIP_VARSTATUS_ORIGINAL )
      {
         errorMessage("Cannot add original variables to transformed problem");
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPprobAddVar(scip->transprob, scip->mem->solvemem, scip->set, var) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

/** gets variables of the problem */
RETCODE SCIPgetVars(
   SCIP*            scip,               /**< SCIP data structure */
   VAR***           vars,               /**< pointer to store variables array or NULL if not needed */
   int*             nvars,              /**< pointer to store number of variables or NULL if not needed */
   int*             nbin,               /**< pointer to store number of binary variables or NULL if not needed */
   int*             nint,               /**< pointer to store number of integer variables or NULL if not needed */
   int*             nimpl,              /**< pointer to store number of implicit integral vars or NULL if not needed */
   int*             ncont               /**< pointer to store number of continous variables or NULL if not needed */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      if( vars != NULL )
         *vars = scip->origprob->vars;
      if( nvars != NULL )
         *nvars = scip->origprob->nvars;
      if( nbin != NULL )
         *nbin = scip->origprob->nbin;
      if( nint != NULL )
         *nint = scip->origprob->nint;
      if( nimpl != NULL )
         *nimpl = scip->origprob->nimpl;
      if( ncont != NULL )
         *ncont = scip->origprob->ncont;
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      if( vars != NULL )
         *vars = scip->transprob->vars;
      if( nvars != NULL )
         *nvars = scip->transprob->nvars;
      if( nbin != NULL )
         *nbin = scip->transprob->nbin;
      if( nint != NULL )
         *nint = scip->transprob->nint;
      if( nimpl != NULL )
         *nimpl = scip->transprob->nimpl;
      if( ncont != NULL )
         *ncont = scip->transprob->ncont;
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

/** finds variable of given name in the problem */
RETCODE SCIPfindVar(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of variable to find */
   VAR**            var                 /**< pointer to store the variable, returns NULL if not found */
   )
{
   assert(name != NULL);
   assert(var != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPfindVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      *var = SCIPprobFindVar(scip->origprob, name);
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
      *var = SCIPprobFindVar(scip->transprob, name);
      if( *var == NULL )
         *var = SCIPprobFindVar(scip->origprob, name);
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

/** adds global constraint to the problem */
RETCODE SCIPaddCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(cons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddCons", FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPprobAddCons(scip->origprob, scip->mem->probmem, scip->set, cons) );
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPprobAddCons(scip->transprob, scip->mem->solvemem, scip->set, cons) );
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPnodeAddCons(scip->tree->root, scip->mem->solvemem, scip->set, cons) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

/** finds constraint of given name in the problem */
RETCODE SCIPfindCons(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of constraint to find */
   CONS**           cons                /**< pointer to store the constraint, returns NULL if not found */
   )
{
   assert(name != NULL);
   assert(cons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPfindCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      *cons = SCIPprobFindCons(scip->origprob, name);
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
      *cons = SCIPprobFindCons(scip->transprob, name);
      if( *cons == NULL )
         *cons = SCIPprobFindCons(scip->origprob, name);
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}




/*
 * local subproblem methods
 */

/** adds local constraint to the given node (and all of its subnodes) */
RETCODE SCIPaddConsNode(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to add constraint to, or NULL for active node */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(cons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddConsNode", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( node == NULL )
      node = scip->tree->actnode;

   CHECK_OKAY( SCIPnodeAddCons(node, scip->mem->solvemem, scip->set, cons) );
   
   return SCIP_OKAY;
}

/** disables constraint's separation, enforcing, and propagation capabilities at the given node (and all subnodes) */
RETCODE SCIPdisableConsNode(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to add constraint to, or NULL for active node */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(cons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPdisableConsNode", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( node == NULL )
      node = scip->tree->actnode;

   CHECK_OKAY( SCIPnodeDisableCons(node, scip->mem->solvemem, scip->set, cons) );
   
   return SCIP_OKAY;
}




/*
 * solve methods
 */

/** solves problem */
RETCODE SCIPsolve(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   char s[255];
   int i;

   CHECK_OKAY( checkStage(scip, "SCIPsolve", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      /* switch stage to INITSOLVE */
      scip->stage = SCIP_STAGE_INITSOLVE;

      /* mark statistics before solving */
      SCIPstatMark(scip->stat);

      /* init solve data structures */
      CHECK_OKAY( SCIPlpCreate(&scip->lp, scip->mem->solvemem, scip->set, SCIPprobGetName(scip->origprob)) );
      CHECK_OKAY( SCIPpriceCreate(&scip->price) );
      CHECK_OKAY( SCIPsepaCreate(&scip->sepa) );
      CHECK_OKAY( SCIPcutpoolCreate(&scip->cutpool, scip->set->agelimit) );
      CHECK_OKAY( SCIPeventfilterCreate(&scip->eventfilter, scip->mem->solvemem) );
      CHECK_OKAY( SCIPeventqueueCreate(&scip->eventqueue) );

      /* init callback methods */
      CHECK_OKAY( SCIPsetInitCallbacks(scip->set) );

      /* copy problem in solve memory */
      CHECK_OKAY( SCIPprobTransform(scip->origprob, scip->mem->solvemem, scip->set, scip->stat, &scip->transprob) );

      /* activate constraints in the problem */
      CHECK_OKAY( SCIPprobActivate(scip->transprob, scip->set) );

      /* create primal solution storage */
      CHECK_OKAY( SCIPprimalCreate(&scip->primal, scip->mem->solvemem, scip->set, scip->transprob, scip->lp) );

      /* create branching candidate storage */
      CHECK_OKAY( SCIPbranchcandCreate(&scip->branchcand, scip->set, scip->transprob) );

      /* switch stage to PRESOLVING */
      scip->stage = SCIP_STAGE_PRESOLVING;

      /* presolve problem */
      todoMessage("problem presolving");

      /* print presolved problem statistics */
      sprintf(s, "presolved problem has %d variables (%d bin, %d int, %d impl, %d cont) and %d constraints",
         scip->transprob->nvars, scip->transprob->nbin, scip->transprob->nint, scip->transprob->nimpl,
         scip->transprob->ncont, scip->transprob->nconss);
      infoMessage(SCIPverbLevel(scip), SCIP_VERBLEVEL_HIGH, s);

      for( i = 0; i < scip->set->nconshdlrs; ++i )
      {
         int nactiveconss;

         nactiveconss = SCIPconshdlrGetNActiveConss(scip->set->conshdlrs[i]);
         if( nactiveconss > 0 )
         {
            sprintf(s, " %5d constraints of type <%s>", nactiveconss, SCIPconshdlrGetName(scip->set->conshdlrs[i]));
            infoMessage(SCIPverbLevel(scip), SCIP_VERBLEVEL_HIGH, s);
         }
      }

      /* create branch-and-bound tree */
      CHECK_OKAY( SCIPtreeCreate(&scip->tree, scip->mem->solvemem, scip->set, scip->stat, scip->lp, scip->transprob) );
      
      /* switch stage to SOLVING */
      scip->stage = SCIP_STAGE_SOLVING;

      /* fallthrough */

   case SCIP_STAGE_SOLVING:
      /* continue solution process */
      infoMessage(SCIPverbLevel(scip), SCIP_VERBLEVEL_NORMAL, "");
      CHECK_OKAY( SCIPsolveCIP(scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->tree, 
                     scip->lp, scip->price, scip->sepa, scip->branchcand, scip->cutpool, scip->primal,
                     scip->eventfilter, scip->eventqueue) );

      /* detect, whether problem is solved */
      todoMessage("detect, whether problem is solved");
      scip->stage = SCIP_STAGE_SOLVED;

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVED:
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

/** frees all solution process data, only original problem is kept */
RETCODE SCIPfreeSolve(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeSolve", TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_PROBLEM:
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      /* switch stage to FREESOLVE */
      scip->stage = SCIP_STAGE_FREESOLVE;

      /* deactivate the active node */
      CHECK_OKAY( SCIPnodeActivate(NULL, scip->mem->solvemem, scip->set, scip->stat, scip->tree, scip->lp, 
                     scip->branchcand, scip->eventqueue) );

      /* deactivate constraints in the problem */
      CHECK_OKAY( SCIPprobDeactivate(scip->transprob) );

      /* clear the LP */
      CHECK_OKAY( SCIPlpClear(scip->lp, scip->mem->solvemem, scip->set) );

      /* free solution process data */
      CHECK_OKAY( SCIPtreeFree(&scip->tree, scip->mem->solvemem, scip->set, scip->lp) );
      CHECK_OKAY( SCIPbranchcandFree(&scip->branchcand) );
      CHECK_OKAY( SCIPprimalFree(&scip->primal, scip->mem->solvemem) );
      CHECK_OKAY( SCIPprobFree(&scip->transprob, scip->mem->solvemem, scip->set, scip->lp) );

      /* exit callback methods */
      CHECK_OKAY( SCIPsetExitCallbacks(scip->set) );

      /* free solve data structures */
      CHECK_OKAY( SCIPeventfilterFree(&scip->eventfilter, scip->mem->solvemem, scip->set) );
      CHECK_OKAY( SCIPeventqueueFree(&scip->eventqueue) );
      CHECK_OKAY( SCIPcutpoolFree(&scip->cutpool, scip->mem->solvemem, scip->set, scip->lp) );
      CHECK_OKAY( SCIPsepaFree(&scip->sepa) );
      CHECK_OKAY( SCIPpriceFree(&scip->price) );
      CHECK_OKAY( SCIPlpFree(&scip->lp, scip->mem->solvemem, scip->set) );

      clearBlockMemoryNull(scip->mem->solvemem);

      /* reset statistics to the point before solving started */
      SCIPstatReset(scip->stat);

      /* switch stage to PROBLEM */
      scip->stage = SCIP_STAGE_PROBLEM;

      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}




/*
 * variable methods
 */

/** create and capture problem variable */
RETCODE SCIPcreateVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            var,                /**< pointer to variable object */
   const char*      name,               /**< name of column */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype             /**< type of variable */
   )
{
   assert(var != NULL);
   assert(name != NULL);
   assert(lb <= ub);

   CHECK_OKAY( checkStage(scip, "SCIPcreateVar", FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPvarCreate(var, scip->mem->probmem, scip->set, scip->stat, name, lb, ub, obj, vartype) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPvarCreateTransformed(var, scip->mem->solvemem, scip->set, scip->stat, name, lb, ub, obj, vartype) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

/** increases usage counter of variable */
RETCODE SCIPcaptureVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to capture */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcaptureVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   SCIPvarCapture(var);

   return SCIP_OKAY;
}

/** decreases usage counter of variable, and frees memory if necessary */
RETCODE SCIPreleaseVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            var                 /**< pointer to variable */
   )
{
   assert(var != NULL);
   assert(*var != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPreleaseVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPvarRelease(var, scip->mem->probmem, scip->set, scip->lp) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
      if( (*var)->varstatus == SCIP_VARSTATUS_ORIGINAL )
      {
         errorMessage("cannot release original variables while solving the problem");
         return SCIP_INVALIDCALL;
      }
      CHECK_OKAY( SCIPvarRelease(var, scip->mem->solvemem, scip->set, scip->lp) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

/** gets solution value for variable in active node */
RETCODE SCIPgetVarSol(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to get solution value for */
   Real*            solval              /**< pointer to store the solution value */
   )
{
   assert(var != NULL);
   assert(solval != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetVarSol", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   *solval = SCIPvarGetSol(var, scip->tree);

   return SCIP_OKAY;
}

/** changes lower bound of variable in the given node; if possible, adjust bound to integral value */
RETCODE SCIPchgVarLbNode(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to change bound at, or NULL for active node */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarLbNode", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( node == NULL )
      node = scip->tree->actnode;
   
   SCIPvarAdjustLb(var, scip->set, &newbound);

   CHECK_OKAY( SCIPnodeAddBoundchg(node, scip->mem->solvemem, scip->set, scip->stat, scip->lp, scip->tree, 
                  scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER) );
   
   return SCIP_OKAY;
}

/** changes upper bound of variable in the given node; if possible, adjust bound to integral value */
RETCODE SCIPchgVarUbNode(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to change bound at, or NULL for active node */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarUbNode", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( node == NULL )
      node = scip->tree->actnode;
   
   SCIPvarAdjustUb(var, scip->set, &newbound);

   CHECK_OKAY( SCIPnodeAddBoundchg(node, scip->mem->solvemem, scip->set, scip->stat, scip->lp, scip->tree, 
                  scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER) );
   
   return SCIP_OKAY;
}

/** depending on SCIP's stage, changes lower bound of variable in the problem, in preprocessing, or in active node;
 *  if possible, adjust bound to integral value
 */
RETCODE SCIPchgVarLb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarLb", FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIPvarAdjustLb(var, scip->set, &newbound);

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPvarChgLb(var, scip->mem->probmem, scip->set, scip->stat, scip->lp, scip->tree,
                     scip->branchcand, scip->eventqueue, newbound) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
      if( var->varstatus == SCIP_VARSTATUS_ORIGINAL )
      {
         errorMessage("cannot change bounds of original variables while presolving the problem");
         return SCIP_INVALIDCALL;
      }
      CHECK_OKAY( SCIPvarChgLb(var, scip->mem->solvemem, scip->set, scip->stat, scip->lp, scip->tree,
                     scip->branchcand, scip->eventqueue, newbound) );
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
      if( var->varstatus == SCIP_VARSTATUS_ORIGINAL )
      {
         errorMessage("cannot change bounds of original variables while solving the problem");
         return SCIP_INVALIDCALL;
      }
      CHECK_OKAY( SCIPnodeAddBoundchg(scip->tree->actnode, scip->mem->solvemem, scip->set, scip->stat, scip->lp,
                     scip->tree, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

/** depending on SCIP's stage, changes upper bound of variable in the problem, in preprocessing, or in active node;
 *  if possible, adjust bound to integral value
 */
RETCODE SCIPchgVarUb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarUb", FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIPvarAdjustUb(var, scip->set, &newbound);

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPvarChgUb(var, scip->mem->probmem, scip->set, scip->stat, scip->lp, scip->tree,
                     scip->branchcand, scip->eventqueue, newbound) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
      if( var->varstatus == SCIP_VARSTATUS_ORIGINAL )
      {
         errorMessage("cannot change bounds of original variables while presolving the problem");
         return SCIP_INVALIDCALL;
      }
      CHECK_OKAY( SCIPvarChgUb(var, scip->mem->solvemem, scip->set, scip->stat, scip->lp, scip->tree,
                     scip->branchcand, scip->eventqueue, newbound) );
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
      if( var->varstatus == SCIP_VARSTATUS_ORIGINAL )
      {
         errorMessage("cannot change bounds of original variables while solving the problem");
         return SCIP_INVALIDCALL;
      }
      CHECK_OKAY( SCIPnodeAddBoundchg(scip->tree->actnode, scip->mem->solvemem, scip->set, scip->stat, scip->lp,
                     scip->tree, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

/** changes type of variable in the problem */
RETCODE SCIPchgVarType(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   VARTYPE          vartype             /**< new type of variable */
   )
{
   assert(var != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPchgVarType", FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(var->varstatus == SCIP_VARSTATUS_ORIGINAL);
      if( var->probindex >= 0 )
      {
         CHECK_OKAY( SCIPprobChgVarType(scip->origprob, var, vartype) );
      }
      else
      {
         CHECK_OKAY( SCIPvarChgType(var, vartype) );
      }
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
      if( var->varstatus == SCIP_VARSTATUS_ORIGINAL )
      {
         errorMessage("cannot change type of original variables while solving the problem");
         return SCIP_INVALIDCALL;
      }
      if( var->probindex >= 0 )
      {
         CHECK_OKAY( SCIPprobChgVarType(scip->transprob, var, vartype) );
      }
      else
      {
         CHECK_OKAY( SCIPvarChgType(var, vartype) );
      }
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }   
}




/*
 * constraint methods
 */

/** creates and captures a constraint of the given constraint handler
 *  Warning! If a constraint is marked to be checked for feasibility but not to be enforced, a LP or pseudo solution
 *  may be declared feasible even if it violates this particular constraint.
 *  This constellation should only be used, if no LP or pseudo solution can violate the constraint -- e.g. if a
 *  local constraint is redundant due to the variable's local bounds.
 */
RETCODE SCIPcreateCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to constraint */
   const char*      name,               /**< name of constraint */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   CONSDATA*        consdata,           /**< data for this specific constraint */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate           /**< should the constraint be propagated during node processing? */
   )
{
   assert(cons != NULL);
   assert(name != NULL);
   assert(conshdlr != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPcreateCons", FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPconsCreate(cons, scip->mem->probmem, name, conshdlr, consdata, separate, enforce, check, propagate,
                     TRUE) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPconsCreate(cons, scip->mem->solvemem, name, conshdlr, consdata, separate, enforce, check, propagate,
                     FALSE) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

/** increases usage counter of constraint */
RETCODE SCIPcaptureCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to capture */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcaptureCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   SCIPconsCapture(cons);

   return SCIP_OKAY;
}

/** decreases usage counter of constraint, and frees memory if necessary */
RETCODE SCIPreleaseCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons                /**< pointer to constraint */
   )
{
   assert(cons != NULL);
   assert(*cons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPreleaseCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIPconsRelease(cons, scip->mem->probmem, scip->set);
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
      if( SCIPconsIsOriginal(*cons) )
      {
         errorMessage("cannot release original constraint while solving the problem");
         return SCIP_INVALIDCALL;
      }
      SCIPconsRelease(cons, scip->mem->solvemem, scip->set);
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}




/*
 * LP methods
 */

/** checks, whether the LP was solved in the active node */
RETCODE SCIPhasActnodeLP(
   SCIP*            scip,               /**< SCIP data structure */
   Bool*            actnodehaslp        /**< pointer to store whether the active node has LP information */
   )
{
   assert(actnodehaslp != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPhasActnodeLP", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   *actnodehaslp = scip->tree->actnodehaslp;

   return SCIP_OKAY;
}

/** gets solution status of actual LP */
RETCODE SCIPgetLPSolstat(
   SCIP*            scip,               /**< SCIP data structure */
   LPSOLSTAT*       lpsolstat           /**< pointer to store the LP solution status */
   )
{
   assert(lpsolstat != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetLPSolstat", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   *lpsolstat = SCIPlpGetSolstat(scip->lp);

   return SCIP_OKAY;
}

/** gets actual LP columns */
RETCODE SCIPgetLPCols(
   SCIP*            scip,               /**< SCIP data structure */
   COL***           cols,               /**< pointer to store the array of LP columns, or NULL */
   int*             ncols               /**< pointer to store the number of LP columns, or NULL */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetLPCols", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );   

   if( cols != NULL )
      *cols = scip->lp->cols;
   if( ncols != NULL )
      *ncols = scip->lp->ncols;

   return SCIP_OKAY;
}

/** gets actual LP rows */
RETCODE SCIPgetLPRows(
   SCIP*            scip,               /**< SCIP data structure */
   ROW***           rows,               /**< pointer to store the array of LP rows, or NULL */
   int*             nrows               /**< pointer to store the number of LP rows, or NULL */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetLPRows", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );   
   
   if( rows != NULL )
      *rows = scip->lp->rows;
   if( nrows != NULL )
      *nrows = scip->lp->nrows;

   return SCIP_OKAY;
}




/*
 * LP row methods
 */

/** creates and captures an LP row */
RETCODE SCIPcreateRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW**            row,                /**< pointer to row */
   const char*      name,               /**< name of row */
   int              len,                /**< number of nonzeros in the row */
   COL**            col,                /**< array with columns of row entries */
   Real*            val,                /**< array with coefficients of row entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
   Bool             local,              /**< is row only valid locally? */
   Bool             modifiable          /**< is row modifiable during node processing (subject to column generation)? */
   )
{
   assert(row != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPcreateRow", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIProwCreate(row, scip->mem->solvemem, scip->set, scip->stat, scip->lp,
                  name, len, col, val, lhs, rhs, local, modifiable) );

   return SCIP_OKAY;
}

/** increases usage counter of LP row */
RETCODE SCIPcaptureRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< row to capture */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcaptureRow", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   SCIProwCapture(row);

   return SCIP_OKAY;
}

/** decreases usage counter of LP row, and frees memory if necessary */
RETCODE SCIPreleaseRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW**            row                 /**< pointer to LP row */
   )
{
   assert(row != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPreleaseRow", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIProwRelease(row, scip->mem->solvemem, scip->set, scip->lp) );
   
   return SCIP_OKAY;
}

/** forbids roundings of variables in row that may violate row */
RETCODE SCIPforbidRowRounding(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPforbidRowRounding", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   SCIProwForbidRounding(row, scip->set);
   
   return SCIP_OKAY;
}

/** allows roundings of variables in row that may violate row */
RETCODE SCIPallowRowRounding(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPallowRowRounding", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   SCIProwAllowRounding(row, scip->set);
   
   return SCIP_OKAY;
}

/** changes left hand side of LP row */
RETCODE SCIPchgRowLhs(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real             lhs                 /**< new left hand side */
   )
{
   assert(row != NULL);

   errorMessage("sides of row must not be changed (not implemented yet)");
   abort();
#if 0
   CHECK_OKAY( checkStage(scip, "SCIPchgRowLhs", FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIProwChgLhs(row, scip->set, scip->lp, lhs) );
   return SCIP_OKAY;
#endif
}

/** changes right hand side of LP row */
RETCODE SCIPchgRowRhs(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real             rhs                 /**< new right hand side */
   )
{
   assert(row != NULL);

   errorMessage("sides of row must not be changed (not implemented yet)");
   abort();
#if 0
   CHECK_OKAY( checkStage(scip, "SCIPchgRowRhs", FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIProwChgRhs(row, scip->set, scip->lp, rhs) );
   return SCIP_OKAY;
#endif
}

/** resolves variable to columns and adds them with the coefficient to the row */
RETCODE SCIPaddVarToRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   VAR*             var,                /**< problem variable */
   Real             val                 /**< value of coefficient */
   )
{
   assert(row != NULL);
   assert(var != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddVarToRow", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPvarAddToRow(var, scip->mem->solvemem, scip->set, scip->lp, scip->stat, row, val) );
   
   return SCIP_OKAY;
}

/** returns the minimal and maximal activity of a row w.r.t. the column's bounds */
RETCODE SCIPgetRowActivityBounds(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real*            minactivity,        /**< pointer to store the minimal activity, or NULL */
   Real*            maxactivity         /**< pointer to store the maximal activity, or NULL */
   )
{
   assert(row != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetRowActivityBounds", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
   
   CHECK_OKAY( SCIProwGetActivityBounds(row, scip->mem->solvemem, scip->set, scip->lp, minactivity, maxactivity) );

   return SCIP_OKAY;
}

/** gets activity bounds for row after setting variable to zero */
RETCODE SCIPgetRowActivityResiduals(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   VAR*             var,                /**< variable to calculate activity residual for */
   Real             val,                /**< coefficient value of variable in linear constraint */
   Real*            minresactivity,     /**< pointer to store the minimal residual activity */
   Real*            maxresactivity      /**< pointer to store the maximal residual activity */
   )
{
   assert(row != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetRowActivityResiduals", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
   
   CHECK_OKAY( SCIProwGetActivityResiduals(row, scip->mem->solvemem, scip->set, scip->lp, var, val, 
                  minresactivity, maxresactivity) );

   return SCIP_OKAY;
}

/** invalidates activity bounds, such that they are recalculated in next get */
RETCODE SCIPinvalidRowActivityBounds(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPinvalidRowActivityBounds", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
   
   CHECK_OKAY( SCIProwInvalidActivityBounds(row) );

   return SCIP_OKAY;
}

/** returns the activity of a row in the last LP solution */
RETCODE SCIPgetRowActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real*            activity            /**< pointer to store the row's activity */
   )
{
   assert(row != NULL);
   assert(activity != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetRowActivity", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
   
   *activity = SCIProwGetActivity(row, scip->stat);
   
   return SCIP_OKAY;
}

/** returns the feasibility of a row in the last LP solution */
RETCODE SCIPgetRowFeasibility(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real*            feasibility         /**< pointer to store the row's feasibility */
   )
{
   assert(row != NULL);
   assert(feasibility != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetRowFeasibility", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
   
   *feasibility = SCIProwGetFeasibility(row, scip->stat);
   
   return SCIP_OKAY;
}

/** returns the activity of a row for the actual pseudo solution */
RETCODE SCIPgetRowPseudoActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real*            pseudoactivity      /**< pointer to store the row's pseudo activity */
   )
{
   assert(row != NULL);
   assert(pseudoactivity != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetRowPseudoActivity", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIProwGetPseudoActivity(row, scip->mem->solvemem, scip->set, scip->lp, pseudoactivity) );
   
   return SCIP_OKAY;
}

/** returns the feasibility of a row for the actual pseudo solution */
RETCODE SCIPgetRowPseudoFeasibility(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real*            pseudofeasibility   /**< pointer to store the row's pseudo feasibility */
   )
{
   assert(row != NULL);
   assert(pseudofeasibility != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetRowPseudoFeasibility", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIProwGetPseudoFeasibility(row, scip->mem->solvemem, scip->set, scip->lp, pseudofeasibility) );
   
   return SCIP_OKAY;
}

/** returns the activity of a row for the given primal solution */
RETCODE SCIPgetRowSolActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   SOL*             sol,                /**< primal CIP solution */
   Real*            solactivity         /**< pointer to store the row's activity for the solution */
   )
{
   assert(row != NULL);
   assert(solactivity != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetRowSolActivity", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIProwGetSolActivity(row, scip->mem->solvemem, scip->set, scip->stat, sol, solactivity) );
   
   return SCIP_OKAY;
}

/** returns the feasibility of a row for the given primal solution */
RETCODE SCIPgetRowSolFeasibility(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   SOL*             sol,                /**< primal CIP solution */
   Real*            solfeasibility      /**< pointer to store the row's feasibility for the solution */
   )
{
   assert(row != NULL);
   assert(solfeasibility != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetRowSolFeasibility", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIProwGetSolFeasibility(row, scip->mem->solvemem, scip->set, scip->stat, sol, solfeasibility) );
   
   return SCIP_OKAY;
}

/** output row to file stream */
RETCODE SCIPprintRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   assert(row != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPprintRow", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIProwPrint(row, file);

   return SCIP_OKAY;
}




/*
 * cutting plane methods
 */

/** adds cut to separation storage */
RETCODE SCIPaddCut(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             cut,                /**< separated cut */
   Real             score               /**< separation score of cut (the larger, the better the cut) */
   )
{
   assert(cut != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddCut", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   assert(scip->tree->actnode != NULL);

   if( !scip->tree->actnodehaslp )
   {
      errorMessage("cannot add cuts, because node LP is not processed");
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( SCIPsepaAddCut(scip->sepa, scip->mem->solvemem, scip->set, scip->lp,
                  cut, score, (scip->tree->actnode->depth == 0)) );
   
   return SCIP_OKAY;
}

/** if not already existing, adds row to global cut pool */
RETCODE SCIPpoolCut(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< cutting plane to add */
   )
{
   assert(row != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPpoolCut", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPcutpoolAddRow(scip->cutpool, scip->mem->solvemem, scip->set, row) );

   return SCIP_OKAY;
}

/** gets actual number of rows in the global cut pool */
RETCODE SCIPgetPoolsize(
   SCIP*            scip,               /**< SCIP data structure */
   int*             poolsize            /**< pointer to store the number of rows in the global cut pool */
   )
{
   assert(poolsize != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetPoolsize", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *poolsize = SCIPcutpoolGetNCuts(scip->cutpool);
   
   return SCIP_OKAY;
}




/*
 * branching methods
 */

/** gets branching candidates for LP solution branching (fractional variables) */
RETCODE SCIPgetLPBranchCands(
   SCIP*            scip,               /**< SCIP data structure */
   VAR***           lpcands,            /**< pointer to store the array of LP branching candidates, or NULL */
   Real**           lpcandssol,         /**< pointer to store the array of LP candidate solution values, or NULL */
   Real**           lpcandsfrac,        /**< pointer to store the array of LP candidate fractionalities, or NULL */
   int*             nlpcands            /**< pointer to store the number of LP branching candidates, or NULL */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetLPBranchCands", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPbranchcandGetLPCands(scip->branchcand, scip->set, scip->stat, scip->lp,
                  lpcands, lpcandssol, lpcandsfrac, nlpcands) );
   
   return SCIP_OKAY;
}

/** gets branching candidates for pseudo solution branching (nonfixed variables) */
RETCODE SCIPgetPseudoBranchCands(
   SCIP*            scip,               /**< SCIP data structure */
   VAR***           pseudocands,        /**< pointer to store the array of pseudo branching candidates, or NULL */
   int*             npseudocands        /**< pointer to store the number of pseudo branching candidates, or NULL */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetLPBranchCands", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPbranchcandGetPseudoCands(scip->branchcand, scip->set, scip->transprob, pseudocands, npseudocands) );

   return SCIP_OKAY;
}

/** creates a child node of the active node */
RETCODE SCIPcreateChild(
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           node                /**< pointer to node data structure */
   )
{
   assert(node != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPcreateChild", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPnodeCreate(node, scip->mem->solvemem, scip->set, scip->tree) );
   
   return SCIP_OKAY;
}

/** branches on a variable; if solution value x' is fractional, two child nodes are created
 *  (x <= floor(x'), x >= ceil(x')), if solution value is integral, three child nodes are created
 *  (x <= x'-1, x == x', x >= x'+1)
 */
RETCODE SCIPbranchVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to branch on */
   )
{
   Real solval;

   assert(var != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPbranchVar", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( SCIPsetIsFixed(scip->set, var->dom.lb, var->dom.ub) )
   {
      char s[255];
      sprintf(s, "cannot branch on variable <%s> with fixed domain [%g,%g]", var->name, var->dom.lb, var->dom.ub);
      errorMessage(s);
      return SCIP_INVALIDDATA;
   }

   CHECK_OKAY( SCIPtreeBranchVar(scip->tree, scip->mem->solvemem, scip->set, scip->stat, scip->lp, scip->branchcand,
                  scip->eventqueue, var) );

   return SCIP_OKAY;
}

/** calls branching rules to branch on an LP solution */
RETCODE SCIPbranchLP(
   SCIP*            scip,               /**< SCIP data structure */
   RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   )
{
   VAR** lpcands;
   int nlpcands;
   int i;

   assert(result != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPbranchLP", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   /* get branching candidates */
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, NULL, NULL, &nlpcands) );
   if( nlpcands == 0 )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTRUN;

   /* try all branching rules until one succeeded to branch */
   for( i = 0; i < scip->set->nbranchrules && *result == SCIP_DIDNOTRUN; ++i )
   {
      CHECK_OKAY( SCIPbranchruleExecLPSol(scip->set->branchrules[i], scip, result) );
   }

   if( *result == SCIP_DIDNOTRUN )
   {
      /* no branching method succeeded in choosing a branching: just branch on the first fractional variable */
      CHECK_OKAY( SCIPbranchVar(scip, lpcands[0]) );
      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}




/*
 * primal solutions
 */

/** creates a primal solution, initialized to zero */
RETCODE SCIPcreateSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateSol", FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPsolCreate(sol, scip->mem->solvemem, scip->stat, heur) );

   return SCIP_OKAY;
}

/** creates a primal solution, initialized to the actual LP solution */
RETCODE SCIPcreateLPSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateLPSol", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( !scip->tree->actnodehaslp )
   {
      errorMessage("LP solution does not exist");
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( SCIPsolCreateLPSol(sol, scip->mem->solvemem, scip->stat, scip->lp, heur) );

   return SCIP_OKAY;
}

/** creates a primal solution, initialized to the actual pseudo solution */
RETCODE SCIPcreatePseudoSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreatePseudoSol", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPsolCreatePseudoSol(sol, scip->mem->solvemem, scip->stat, scip->tree, heur) );

   return SCIP_OKAY;
}

/** creates a primal solution, initialized to the actual solution */
RETCODE SCIPcreateActSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateActSol", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPsolCreateActSol(sol, scip->mem->solvemem, scip->stat, scip->tree, scip->lp, heur) );

   return SCIP_OKAY;
}

/** frees primal CIP solution */
RETCODE SCIPfreeSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol                 /**< pointer to the solution */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeSol", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsolFree(sol, scip->mem->solvemem) );

   return SCIP_OKAY;
}

/** links a primal solution to the actual LP solution */
RETCODE SCIPlinkLPSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPlinkLPSol", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( !scip->tree->actnodehaslp )
   {
      errorMessage("LP solution does not exist");
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( SCIPsolLinkLPSol(sol, scip->mem->solvemem, scip->stat, scip->lp) );

   return SCIP_OKAY;
}

/** links a primal solution to the actual pseudo solution */
RETCODE SCIPlinkPseudoSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPlinkPseudoSol", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPsolLinkPseudoSol(sol, scip->mem->solvemem, scip->stat, scip->tree) );

   return SCIP_OKAY;
}

/** links a primal solution to the actual LP or pseudo solution */
RETCODE SCIPlinkActSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPlinkActSol", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPsolLinkActSol(sol, scip->mem->solvemem, scip->stat, scip->tree, scip->lp) );

   return SCIP_OKAY;
}

/** clears a primal solution */
RETCODE SCIPclearSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPclearSol", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsolClear(sol, scip->mem->solvemem, scip->stat) );

   return SCIP_OKAY;
}

/** stores solution values of variables in solution's own array */
RETCODE SCIPunlinkSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPunlinkSol", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsolUnlink(sol, scip->mem->solvemem, scip->set, scip->transprob) );

   return SCIP_OKAY;
}

/** sets value of variable in primal CIP solution */
RETCODE SCIPsetSolVal(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   VAR*             var,                /**< variable to add to solution */
   Real             val                 /**< solution value of variable */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetSolVal", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsolSetVal(sol, scip->mem->solvemem, scip->set, scip->stat, var, val) );

   return SCIP_OKAY;
}

/** increases value of variable in primal CIP solution */
RETCODE SCIPincSolVal(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   VAR*             var,                /**< variable to increase solution value for */
   Real             incval              /**< increment for solution value of variable */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPincSolVal", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsolIncVal(sol, scip->mem->solvemem, scip->set, scip->stat, var, incval) );

   return SCIP_OKAY;
}

/** returns value of variable in primal CIP solution, or in actual LP/pseudo solution */
RETCODE SCIPgetSolVal(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution, or NULL for actual LP/pseudo solution */
   VAR*             var,                /**< variable to get value for */
   Real*            solval              /**< pointer to store the solution value */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetSolVal", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( sol != NULL )
   {
      CHECK_OKAY( SCIPsolGetVal(sol, scip->mem->solvemem, scip->set, scip->stat, var, solval) );
   }
   else
   {
      CHECK_OKAY( checkStage(scip, "SCIPgetSolVal(sol==NULL)", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
      *solval = SCIPvarGetSol(var, scip->tree);
   }

   return SCIP_OKAY;
}

/** returns objective value of primal CIP solution, or actual LP/pseudo objective value */
RETCODE SCIPgetSolObj(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution, or NULL for actual LP/pseudo objective value */
   Real*            objval              /**< pointer to store the objective value */
   )
{
   assert(objval != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetSolObj", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( sol != NULL )
      *objval = SCIPprobExternObjval(scip->origprob, SCIPprobExternObjval(scip->transprob, SCIPsolGetObj(sol)));
   else
   {
      CHECK_OKAY( checkStage(scip, "SCIPgetSolObj(sol==NULL)", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
      if( scip->tree->actnodehaslp )
      {
         assert(scip->lp->solved);
         *objval = SCIPprobExternObjval(scip->origprob, SCIPprobExternObjval(scip->transprob, scip->lp->objval));
      }
      else
         *objval = SCIPprobExternObjval(scip->origprob, SCIPprobExternObjval(scip->transprob, scip->tree->actpseudoobjval));
   }
   
   return SCIP_OKAY;
}

/** returns transformed objective value of primal CIP solution, or transformed actual LP/pseudo objective value */
RETCODE SCIPgetSolTransObj(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution, or NULL for actual LP/pseudo objective value */
   Real*            objval              /**< pointer to store the objective value */
   )
{
   assert(objval != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetSolTransObj", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( sol != NULL )
      *objval = SCIPsolGetObj(sol);
   else
   {
      CHECK_OKAY( checkStage(scip, "SCIPgetSolTransObj(sol==NULL)", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
      if( scip->tree->actnodehaslp )
      {
         assert(scip->lp->solved);
         *objval = scip->lp->objval;
      }
      else
         *objval = scip->tree->actpseudoobjval;
   }
   
   return SCIP_OKAY;
}

/** gets node number, where this solution was found */
RETCODE SCIPgetSolNodenum(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   Longint*         nodenum             /**< pointer to store the node number */
   )
{
   assert(nodenum != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetSolNodenum", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   *nodenum = SCIPsolGetNodenum(sol);

   return SCIP_OKAY;
}

/** gets heuristic, that found this solution (or NULL if it's from the tree) */
RETCODE SCIPgetSolHeur(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   HEUR**           heur                /**< pointer to store the heuristic */
   )
{
   assert(heur != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetSolHeur", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   *heur = SCIPsolGetHeur(sol);

   return SCIP_OKAY;
}

/** outputs non-zero original variables of solution to file stream */
RETCODE SCIPprintSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPprintSol", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( file == NULL )
      file = stdout;

   fprintf(file, "objective value:                 %f\n",
      SCIPprobExternObjval(scip->origprob, SCIPprobExternObjval(scip->transprob, SCIPsolGetObj(sol))));
   
   CHECK_OKAY( SCIPsolPrint(sol, scip->mem->solvemem, scip->set, scip->stat, scip->origprob, file) );

   return SCIP_OKAY;
}

/** outputs non-zero transformed variables of solution to file stream */
RETCODE SCIPprintTransSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPprintSolTrans", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( file == NULL )
      file = stdout;

   fprintf(file, "transformed objective value:     %f\n", SCIPsolGetObj(sol));
   
   CHECK_OKAY( SCIPsolPrint(sol, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, file) );
   
   return SCIP_OKAY;
}

/** gets best feasible primal solution found so far */
RETCODE SCIPgetBestSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol                 /**< pointer to store the solution, returns NULL if no solution available */
   )
{
   assert(sol != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetBestSol", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   if( scip->primal->nsols > 0 ) 
   {
      assert(scip->primal->sols != NULL);
      assert(scip->primal->sols[0] != NULL);
      *sol = scip->primal->sols[0];
   }
   else
      *sol = NULL;
   
   return SCIP_OKAY;
}

/** outputs best feasible primal solution found so far to file stream */
RETCODE SCIPprintBestSol(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   SOL* sol;

   CHECK_OKAY( checkStage(scip, "SCIPprintBestSol", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   CHECK_OKAY( SCIPgetBestSol(scip, &sol) );

   if( file == NULL )
      file = stdout;

   if( sol == NULL )
      fprintf(file, "no solution available\n");
   else
      SCIPprintSol(scip, sol, file);
   
   return SCIP_OKAY;
}

/** adds feasible primal solution to solution storage by moving it */
RETCODE SCIPaddSolMove(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol                 /**< pointer to primal CIP solution; is cleared in function call */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddSolMove", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPprimalAddSolMove(scip->primal, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->tree, 
                  scip->lp, scip->eventfilter, sol) );

   return SCIP_OKAY;
}

/** adds feasible primal solution to solution storage by copying it */
RETCODE SCIPaddSolCopy(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal CIP solution */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddSolCopy", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPprimalAddSolCopy(scip->primal, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->tree, 
                  scip->lp, scip->eventfilter, sol) );

   return SCIP_OKAY;
}

/** checks solution for feasibility; if possible, adds it to storage by moving */
RETCODE SCIPtrySolMove(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   Bool             chckintegrality,    /**< has integrality to be checked? */
   Bool             chcklprows,         /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPtrySolMove", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPprimalTrySolMove(scip->primal, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->tree, 
                  scip->lp, scip->eventfilter, sol, chckintegrality, chcklprows, stored) );

   return SCIP_OKAY;
}

/** checks solution for feasibility; if possible, adds it to storage by copying */
RETCODE SCIPtrySolCopy(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal CIP solution */
   Bool             chckintegrality,    /**< has integrality to be checked? */
   Bool             chcklprows,         /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPtrySolCopy", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPprimalTrySolCopy(scip->primal, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->tree, 
                  scip->lp, scip->eventfilter, sol, chckintegrality, chcklprows, stored) );

   return SCIP_OKAY;
}




/*
 * event methods
 */

/** catches a global (not variable dependent) event */
RETCODE SCIPcatchEvent(
   SCIP*            scip,               /**< SCIP data structure */
   EVENTTYPE        eventtype,          /**< event type mask to select events to catch */
   EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler when processing this event */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcatchEvent", FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPeventfilterAdd(scip->eventfilter, scip->mem->solvemem, scip->set, eventtype, eventhdlr, eventdata) );

   return SCIP_OKAY;
}

/** drops a global event (stops to track event) */
RETCODE SCIPdropEvent(
   SCIP*            scip,               /**< SCIP data structure */
   EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler when processing this event */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPdropEvent", FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE) );

   CHECK_OKAY( SCIPeventfilterDel(scip->eventfilter, scip->mem->solvemem, scip->set, eventhdlr, eventdata) );
   
   return SCIP_OKAY;
}

/** catches a domain change event on the given variable */
RETCODE SCIPcatchVarEvent(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to catch event for */
   EVENTTYPE        eventtype,          /**< event type mask to select events to catch */
   EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler when processing this event */
   )
{
   assert(var != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPcatchVarEvent", FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   if( var->varstatus == SCIP_VARSTATUS_ORIGINAL )
   {
      errorMessage("cannot catch events for original variables");
      return SCIP_INVALIDDATA;
   }

   if( (eventtype & SCIP_EVENTTYPE_DOMCHANGED) == 0 )
   {
      errorMessage("event is not a domain change event");
      return SCIP_INVALIDDATA;
   }

   CHECK_OKAY( SCIPeventfilterAdd(var->eventfilter, scip->mem->solvemem, scip->set, eventtype, eventhdlr, eventdata) );

   return SCIP_OKAY;
}

/** drops a domain change event (stops to track event) on the given variable */
RETCODE SCIPdropVarEvent(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to drop event for */
   EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler when processing this event */
   )
{
   assert(var != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPdropVarEvent", FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE) );

   if( var->varstatus == SCIP_VARSTATUS_ORIGINAL )
   {
      errorMessage("cannot drop events for original variables");
      return SCIP_INVALIDDATA;
   }

   CHECK_OKAY( SCIPeventfilterDel(var->eventfilter, scip->mem->solvemem, scip->set, eventhdlr, eventdata) );
   
   return SCIP_OKAY;
}




/*
 * tree methods
 */

/** gets children of active node */
RETCODE SCIPgetChildren(
   SCIP*            scip,               /**< SCIP data structure */
   NODE***          children,           /**< pointer to store children array */
   int*             nchildren           /**< pointer to store number of children */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetChildren", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   *children = scip->tree->children;
   *nchildren = scip->tree->nchildren;
   
   return SCIP_OKAY;
}

/** gets siblings of active node */
RETCODE SCIPgetSiblings(
   SCIP*            scip,               /**< SCIP data structure */
   NODE***          siblings,           /**< pointer to store siblings array */
   int*             nsiblings           /**< pointer to store number of siblings */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetSiblings", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   *siblings = scip->tree->siblings;
   *nsiblings = scip->tree->nsiblings;
   
   return SCIP_OKAY;
}

/** gets the best child of the active node */
RETCODE SCIPgetBestChild(
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           bestchild           /**< pointer to store best child */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetBestChild", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   *bestchild = SCIPtreeGetBestChild(scip->tree, scip->set);
   
   return SCIP_OKAY;
}

/** gets the best sibling of the active node */
RETCODE SCIPgetBestSibling(
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           bestsibling         /**< pointer to store best sibling */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetBestSibling", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   *bestsibling = SCIPtreeGetBestSibling(scip->tree, scip->set);
   
   return SCIP_OKAY;
}

/** gets the best leaf from the node queue */
RETCODE SCIPgetBestLeaf(
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           bestleaf            /**< pointer to store best leaf */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetBestLeaf", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   *bestleaf = SCIPtreeGetBestLeaf(scip->tree);
   
   return SCIP_OKAY;
}

/** gets the best node from the tree (child, sibling, or leaf) */
RETCODE SCIPgetBestNode(
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           bestnode            /**< pointer to store best node */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetBestNode", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   *bestnode = SCIPtreeGetBestNode(scip->tree, scip->set);
   
   return SCIP_OKAY;
}




/*
 * statistic methods
 */

/** gets number of processed nodes, including the active node */
RETCODE SCIPgetNodenum(
   SCIP*            scip,               /**< SCIP data structure */
   Longint*         nodenum             /**< pointer to store the number of processed nodes */
   )
{
   assert(nodenum != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetNodenum", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   *nodenum = scip->stat->nnodes;

   return SCIP_OKAY;
}

/** gets number of nodes left in the tree (children + siblings + leaves) */
RETCODE SCIPgetNNodesLeft(
   SCIP*            scip,               /**< SCIP data structure */
   int*             nnodes              /**< pointer to store the number of processed nodes */
   )
{
   assert(nnodes != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetNNodesLeft", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *nnodes = SCIPtreeGetNNodes(scip->tree);

   return SCIP_OKAY;
}

/** gets total number of simplex iterations used so far */
RETCODE SCIPgetNLPIterations(
   SCIP*            scip,               /**< SCIP data structure */
   int*             lpiterations        /**< pointer to store the iterations */
   )
{
   assert(lpiterations != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetNLPIterations", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *lpiterations = scip->stat->nlpiterations;
   
   return SCIP_OKAY;
}

/** gets depth of active node */
RETCODE SCIPgetActDepth(
   SCIP*            scip,               /**< SCIP data structure */
   int*             actdepth            /**< pointer to store the depth */
   )
{
   assert(actdepth != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetActDepth", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   if( scip->tree->actnode != NULL )
      *actdepth = scip->tree->actnode->depth;
   else
      *actdepth = -1;
   
   return SCIP_OKAY;
}

/** gets maximal depth of all processed nodes */
RETCODE SCIPgetMaxDepth(
   SCIP*            scip,               /**< SCIP data structure */
   int*             maxdepth            /**< pointer to store the depth */
   )
{
   assert(maxdepth != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetMaxDepth", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *maxdepth = scip->stat->maxdepth;
   
   return SCIP_OKAY;
}

/** gets actual plunging depth (succ. times, a child was selected as next node) */
RETCODE SCIPgetPlungeDepth(
   SCIP*            scip,               /**< SCIP data structure */
   int*             plungedepth         /**< pointer to store the depth */
   )
{
   assert(plungedepth != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetPlungeDepth", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   *plungedepth = scip->stat->plungedepth;
   
   return SCIP_OKAY;
}

/** gets total number of active constraints at the current node */
RETCODE SCIPgetNActiveConss(
   SCIP*            scip,               /**< SCIP data structure */
   int*             nactiveconss        /**< pointer to store the number of active constraints */
   )
{
   int h;

   assert(nactiveconss != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetNActiveConss", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   *nactiveconss = 0;
   for( h = 0; h < scip->set->nconshdlrs; ++h )
      *nactiveconss += SCIPconshdlrGetNActiveConss(scip->set->conshdlrs[h]);

   return SCIP_OKAY;
}

/** gets total number of enabled constraints at the current node */
RETCODE SCIPgetNEnabledConss(
   SCIP*            scip,               /**< SCIP data structure */
   int*             nenabledconss       /**< pointer to store the number of enabled constraints */
   )
{
   int h;

   assert(nenabledconss != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetNEnabledConss", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   *nenabledconss = 0;
   for( h = 0; h < scip->set->nconshdlrs; ++h )
      *nenabledconss += SCIPconshdlrGetNEnabledConss(scip->set->conshdlrs[h]);

   return SCIP_OKAY;
}

/** gets dual bound of active node */
RETCODE SCIPgetActDualBound(
   SCIP*            scip,               /**< SCIP data structure */
   Real*            actdualbound        /**< pointer to store the dual bound */
   )
{
   assert(actdualbound != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetActDualBound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *actdualbound = SCIPprobExternObjval(scip->origprob, 
      SCIPprobExternObjval(scip->transprob, SCIPtreeGetActLowerbound(scip->tree)));
   
   return SCIP_OKAY;
}

/** gets dual bound of active node in transformed problem */
RETCODE SCIPgetActTransDualBound(
   SCIP*            scip,               /**< SCIP data structure */
   Real*            actdualbound        /**< pointer to store the dual bound */
   )
{
   assert(actdualbound != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetActTransDualBound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *actdualbound = SCIPtreeGetActLowerbound(scip->tree);
   
   return SCIP_OKAY;
}

/** gets average dual bound of all unprocessed nodes */
RETCODE SCIPgetAvgDualBound(
   SCIP*            scip,               /**< SCIP data structure */
   Real*            avgdualbound        /**< pointer to store the average dual bound */
   )
{
   assert(avgdualbound != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetAvgDualBound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *avgdualbound = SCIPprobExternObjval(scip->origprob, 
      SCIPprobExternObjval(scip->transprob, SCIPtreeGetAvgLowerbound(scip->tree, scip->primal->upperbound)));
   
   return SCIP_OKAY;
}

/** gets average dual bound of all unprocessed nodes in transformed problem */
RETCODE SCIPgetAvgTransDualBound(
   SCIP*            scip,               /**< SCIP data structure */
   Real*            avgdualbound        /**< pointer to store the average dual bound */
   )
{
   assert(avgdualbound != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetAvgTransDualBound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *avgdualbound = SCIPtreeGetAvgLowerbound(scip->tree, scip->primal->upperbound);
   
   return SCIP_OKAY;
}

/** gets global dual bound */
RETCODE SCIPgetDualBound(
   SCIP*            scip,               /**< SCIP data structure */
   Real*            dualbound           /**< pointer to store the dual bound */
   )
{
   assert(dualbound != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetDualBound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *dualbound = SCIPprobExternObjval(scip->origprob, 
      SCIPprobExternObjval(scip->transprob, SCIPtreeGetLowerbound(scip->tree, scip->set)));
   
   return SCIP_OKAY;
}

/** gets global dual bound in transformed problem */
RETCODE SCIPgetTransDualBound(
   SCIP*            scip,               /**< SCIP data structure */
   Real*            dualbound           /**< pointer to store the dual bound */
   )
{
   assert(dualbound != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetTransDualBound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *dualbound = SCIPtreeGetLowerbound(scip->tree, scip->set);
   
   return SCIP_OKAY;
}

/** gets global primal bound */
RETCODE SCIPgetPrimalBound(
   SCIP*            scip,               /**< SCIP data structure */
   Real*            primalbound         /**< pointer to store the primal bound */
   )
{
   assert(primalbound != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetPrimalBound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *primalbound = SCIPprobExternObjval(scip->origprob, SCIPprobExternObjval(scip->transprob, scip->primal->upperbound));
   
   return SCIP_OKAY;
}

/** gets global primal bound in transformed problem */
RETCODE SCIPgetTransPrimalBound(
   SCIP*            scip,               /**< SCIP data structure */
   Real*            primalbound         /**< pointer to store the primal bound */
   )
{
   assert(primalbound != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetTransPrimalBound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *primalbound = scip->primal->upperbound;
   
   return SCIP_OKAY;
}




/*
 * numeric values and comparisons
 */

/** returns value treated as infinity */
Real SCIPinfinity(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->infinity;
}

/** returns value treated as zero */
Real SCIPepsilon(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->epsilon;
}

/** returns value treated as zero for sums of floating point values */
Real SCIPsumepsilon(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->sumepsilon;
}

/** checks, if values are in range of epsilon */
Bool SCIPisEQ(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsEQ(scip->set, val1, val2);
}

/** checks, if val1 is (more than epsilon) lower than val2 */
Bool SCIPisLT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsLT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than epsilon) greater than val2 */
Bool SCIPisLE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsLE(scip->set, val1, val2);
}

/** checks, if val1 is (more than epsilon) greater than val2 */
Bool SCIPisGT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsGT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than epsilon) lower than val2 */
Bool SCIPisGE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsGE(scip->set, val1, val2);
}

/** checks, if value is in range epsilon of 0.0 */
Bool SCIPisZero(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsZero(scip->set, val);
}

/** checks, if value is greater than epsilon */
Bool SCIPisPositive(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsPositive(scip->set, val);
}

/** checks, if value is lower than -epsilon */
Bool SCIPisNegative(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsNegative(scip->set, val);
}

/** checks, if values are in range of sumepsilon */
Bool SCIPisSumEQ(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumEQ(scip->set, val1, val2);
}

/** checks, if val1 is (more than sumepsilon) lower than val2 */
Bool SCIPisSumLT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumLT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than sumepsilon) greater than val2 */
Bool SCIPisSumLE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumLE(scip->set, val1, val2);
}

/** checks, if val1 is (more than sumepsilon) greater than val2 */
Bool SCIPisSumGT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumGT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than sumepsilon) lower than val2 */
Bool SCIPisSumGE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumGE(scip->set, val1, val2);
}

/** checks, if value is in range sumepsilon of 0.0 */
Bool SCIPisSumZero(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumZero(scip->set, val);
}

/** checks, if value is greater than sumepsilon */
Bool SCIPisSumPositive(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumPositive(scip->set, val);
}

/** checks, if value is lower than -sumepsilon */
Bool SCIPisSumNegative(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumNegative(scip->set, val);
}

/** checks, if values are in range of feasibility tolerance */
Bool SCIPisFeasEQ(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasEQ(scip->set, val1, val2);
}

/** checks, if val1 is (more than feasibility tolerance) lower than val2 */
Bool SCIPisFeasLT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasLT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than feasibility tolerance) greater than val2 */
Bool SCIPisFeasLE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasLE(scip->set, val1, val2);
}

/** checks, if val1 is (more than feasibility tolerance) greater than val2 */
Bool SCIPisFeasGT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasGT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than feasibility tolerance) lower than val2 */
Bool SCIPisFeasGE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasGE(scip->set, val1, val2);
}

/** checks, if value is in range feasibility tolerance of 0.0 */
Bool SCIPisFeasZero(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasZero(scip->set, val);
}

/** checks, if value is greater than feasibility tolerance */
Bool SCIPisFeasPositive(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasPositive(scip->set, val);
}

/** checks, if value is lower than -feasibility tolerance */
Bool SCIPisFeasNegative(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasNegative(scip->set, val);
}

/** checks, if relative difference of values is in range of epsilon */
Bool SCIPisRelEQ(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsRelEQ(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is lower than epsilon */
Bool SCIPisRelLT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsRelLT(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is not greater than epsilon */
Bool SCIPisRelLE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsRelLE(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is greater than epsilon */
Bool SCIPisRelGT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsRelGT(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is not lower than -epsilon */
Bool SCIPisRelGE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsRelGE(scip->set, val1, val2);
}

/** checks, if rel. difference of values is in range of sumepsilon */
Bool SCIPisSumRelEQ(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumRelEQ(scip->set, val1, val2);
}

/** checks, if rel. difference of val1 and val2 is lower than sumepsilon */
Bool SCIPisSumRelLT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumRelLT(scip->set, val1, val2);
}

/** checks, if rel. difference of val1 and val2 is not greater than sumepsilon */
Bool SCIPisSumRelLE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumRelLE(scip->set, val1, val2);
}

/** checks, if rel. difference of val1 and val2 is greater than sumepsilon */
Bool SCIPisSumRelGT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumRelGT(scip->set, val1, val2);
}

/** checks, if rel. difference of val1 and val2 is not lower than -sumepsilon */
Bool SCIPisSumRelGE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumRelGE(scip->set, val1, val2);
}

/** checks, if value is (positive) infinite */
Bool SCIPisInfinity(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against infinity */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsInfinity(scip->set, val);
}

/** rounds value down to the next integer */
Real SCIPfloor(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetFloor(scip->set, val);
}

/** rounds value up to the next integer */
Real SCIPceil(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetCeil(scip->set, val);
}

/** checks, if value is integral within the LP feasibility bounds */
Bool SCIPisIntegral(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsIntegral(scip->set, val);
}

/** checks, if value is non-negative within the LP feasibility bounds */
Bool SCIPisFeasible(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasible(scip->set, val);
}




/*
 * memory management
 */

/** returns block memory to use at the current time */
MEMHDR* SCIPmemhdr(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(scip->mem != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_PROBLEM:
      return scip->mem->probmem;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
      return scip->mem->solvemem;

   default:
      errorMessage("Unknown SCIP stage");
      return NULL;
   }
}

/** calculate memory size for dynamically allocated arrays */
int SCIPcalcMemGrowSize(
   SCIP*            scip,               /**< SCIP data structure */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);

   return SCIPsetCalcMemGrowSize(scip->set, num);
}

/** gets a memory buffer with at least the given size */
RETCODE SCIPcaptureBuffer(
   SCIP*            scip,               /**< SCIP data structure */
   void**           ptr,                /**< pointer to store the buffer */
   int              size                /**< required size in bytes of buffer */
   )
{
   assert(ptr != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPcaptureBuffer", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetCaptureBufferSize(scip->set, ptr, size) );

   return SCIP_OKAY;
}

/** releases a memory buffer */
RETCODE SCIPreleaseBuffer(
   SCIP*            scip,               /**< SCIP data structure */
   void**           ptr,                /**< pointer to store the buffer */
   int              dummysize           /**< used to get a safer define for SCIPreleaseBufferSize/Array */
   )
{
   assert(ptr != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPreleaseBuffer", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetReleaseBufferSize(scip->set, ptr);
   
   return SCIP_OKAY;
}




/*
 * dynamic arrays
 */

/** creates a dynamic array of real values */
RETCODE SCIPcreateRealarray(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY**      realarray           /**< pointer to store the real array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPrealarrayCreate(realarray, SCIPmemhdr(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of real values */
RETCODE SCIPfreeRealarray(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY**      realarray           /**< pointer to the real array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPrealarrayFree(realarray, SCIPmemhdr(scip)) );
   
   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
RETCODE SCIPextendRealarray(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray,          /**< dynamic real array */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPextendRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPrealarrayExtend(realarray, SCIPmemhdr(scip), scip->set, minidx, maxidx) );
   
   return SCIP_OKAY;
}

/** clears a dynamic real array */
RETCODE SCIPclearRealarray(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray           /**< dynamic real array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPclearRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   SCIPrealarrayClear(realarray);
   
   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
Real SCIPgetRealarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray,          /**< dynamic real array */
   int              idx                 /**< array index to get value for */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   return SCIPrealarrayGetVal(realarray, idx);
   
   return SCIP_OKAY;
}

/** sets value of entry in dynamic array */
RETCODE SCIPsetRealarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray,          /**< dynamic real array */
   int              idx,                /**< array index to set value for */
   Real             val                 /**< value to set array index to */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPrealarraySetVal(realarray, SCIPmemhdr(scip), scip->set, idx, val) );
   
   return SCIP_OKAY;
}

/** increases value of entry in dynamic array */
RETCODE SCIPincRealarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray,          /**< dynamic real array */
   int              idx,                /**< array index to increase value for */
   Real             incval              /**< value to increase array index */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPincRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPrealarrayIncVal(realarray, SCIPmemhdr(scip), scip->set, idx, incval) );
   
   return SCIP_OKAY;
}

/** creates a dynamic array of int values */
RETCODE SCIPcreateIntarray(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY**       intarray            /**< pointer to store the int array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPintarrayCreate(intarray, SCIPmemhdr(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of int values */
RETCODE SCIPfreeIntarray(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY**       intarray            /**< pointer to the int array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPintarrayFree(intarray, SCIPmemhdr(scip)) );
   
   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
RETCODE SCIPextendIntarray(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray,           /**< dynamic int array */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPextendIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPintarrayExtend(intarray, SCIPmemhdr(scip), scip->set, minidx, maxidx) );
   
   return SCIP_OKAY;
}

/** clears a dynamic int array */
RETCODE SCIPclearIntarray(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray            /**< dynamic int array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPclearIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   SCIPintarrayClear(intarray);
   
   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
int SCIPgetIntarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray,           /**< dynamic int array */
   int              idx                 /**< array index to get value for */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   return SCIPintarrayGetVal(intarray, idx);
   
   return SCIP_OKAY;
}

/** sets value of entry in dynamic array */
RETCODE SCIPsetIntarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray,           /**< dynamic int array */
   int              idx,                /**< array index to set value for */
   int              val                 /**< value to set array index to */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPintarraySetVal(intarray, SCIPmemhdr(scip), scip->set, idx, val) );
   
   return SCIP_OKAY;
}

/** increases value of entry in dynamic array */
RETCODE SCIPincIntarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray,           /**< dynamic int array */
   int              idx,                /**< array index to increase value for */
   int              incval              /**< value to increase array index */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPincIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPintarrayIncVal(intarray, SCIPmemhdr(scip), scip->set, idx, incval) );
   
   return SCIP_OKAY;
}

/** creates a dynamic array of bool values */
RETCODE SCIPcreateBoolarray(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY**      boolarray           /**< pointer to store the bool array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPboolarrayCreate(boolarray, SCIPmemhdr(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of bool values */
RETCODE SCIPfreeBoolarray(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY**      boolarray           /**< pointer to the bool array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPboolarrayFree(boolarray, SCIPmemhdr(scip)) );
   
   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
RETCODE SCIPextendBoolarray(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPextendBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPboolarrayExtend(boolarray, SCIPmemhdr(scip), scip->set, minidx, maxidx) );
   
   return SCIP_OKAY;
}

/** clears a dynamic bool array */
RETCODE SCIPclearBoolarray(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPclearBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   SCIPboolarrayClear(boolarray);
   
   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
Bool SCIPgetBoolarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int              idx                 /**< array index to get value for */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   return SCIPboolarrayGetVal(boolarray, idx);
   
   return SCIP_OKAY;
}

/** sets value of entry in dynamic array */
RETCODE SCIPsetBoolarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int              idx,                /**< array index to set value for */
   Bool             val                 /**< value to set array index to */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPboolarraySetVal(boolarray, SCIPmemhdr(scip), scip->set, idx, val) );
   
   return SCIP_OKAY;
}




/*
 * debugging methods
 */

#ifndef NDEBUG

/** prints output about used memory */
void SCIPdebugMemory(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   memoryDiagnostic();

   printf("\nProblem Block Memory (%p):\n", scip->mem->probmem);
   blockMemoryDiagnostic(scip->mem->probmem);

   printf("\nSolution Block Memory (%p):\n", scip->mem->solvemem);
   blockMemoryDiagnostic(scip->mem->solvemem);
}

#endif

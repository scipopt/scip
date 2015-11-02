/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file   heur_repair.c
 * @brief  repair primal heuristic
 * @author Gregor Hendel, Thomas Nagel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define SCIP_STATISTIC
#include <assert.h>
#include <string.h>

#include <stdio.h>
#include <stdlib.h>

#include "scip/heur_repair.h"
#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"
#include "cons_varbound.h"

#define HEUR_NAME             "repair"
#define HEUR_DESC             "repair heuristic"
#define HEUR_DISPCHAR         '!'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_FILENAME      "-" /**< file name of a solution to be used as infeasible starting point */
#define DEFAULT_ROUNDIT       TRUE /**< if it is True : fractional variables which are not fractional in the given
                                    *    solution are rounded, if it is FALSE : solvingprocess of this heuristic is stoped*/

/*
 * Data structures
 */


/** primal heuristic data */
struct SCIP_HeurData
{
   char*                 filename;           /**< file name of a solution to be used as infeasible starting point */
   SCIP_Bool             roundit;            /**< if it is True : fractional variables which are not fractional in the given
                                              *    solution are rounded, if it is FALSE : solvingprocess of this heuristic is stoped*/
   int                   subnodes;
   int                   subiters;
   SCIP_Real             subpresoltime;
   int                   runs;

   int                   ninvalidvars;
   int                   norvars;
   SCIP_Real             relinvalidvars;
   int                   ninvalidcons;
   int                   norcons;
   SCIP_Real             relinvalidcons;
};


/*
 * Local methods
 */

/** reads a given SCIP solution file, problem has to be transformed in advance */
static
SCIP_RETCODE readSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,
   const char*           fname               /**< name of the input file */
   )
{
   SCIP_FILE* file;
   SCIP_Bool error;
   SCIP_Bool unknownvariablemessage;
   SCIP_Bool usevartable;
   int lineno;

   assert(scip != NULL);
   assert(fname != NULL);

   SCIP_CALL( SCIPgetBoolParam(scip, "misc/usevartable", &usevartable) );

   if( !usevartable )
   {
      SCIPerrorMessage("Cannot read solution file if vartable is disabled. Make sure parameter 'misc/usevartable' is set to TRUE.\n");
      return SCIP_READERROR;
   }

   /* open input file */
   file = SCIPfopen(fname, "r");
   if( file == NULL )
   {
      if( strcmp(fname, DEFAULT_FILENAME) != 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
               "Warning: Solution to repair not found!\n");
      }
      return SCIP_OKAY;
   }

   /* read the file */
   error = FALSE;
   unknownvariablemessage = FALSE;
   lineno = 0;
   while( !SCIPfeof(file) && !error )
   {
      char buffer[SCIP_MAXSTRLEN];
      char varname[SCIP_MAXSTRLEN];
      char valuestring[SCIP_MAXSTRLEN];
      char objstring[SCIP_MAXSTRLEN];
      SCIP_VAR* var;
      SCIP_Real value;
      int nread;

      /* get next line */
      if( SCIPfgets(buffer, (int) sizeof(buffer), file) == NULL )
         break;
      lineno++;

      /* there are some lines which may proceed the solution information */
      if( strncasecmp(buffer, "solution status:", 16) == 0 || strncasecmp(buffer, "objective value:", 16) == 0 ||
         strncasecmp(buffer, "Log started", 11) == 0 || strncasecmp(buffer, "Variable Name", 13) == 0 ||
         strncasecmp(buffer, "All other variables", 19) == 0 || strncasecmp(buffer, "\n", 1) == 0 ||
         strncasecmp(buffer, "NAME", 4) == 0 || strncasecmp(buffer, "ENDATA", 6) == 0 )    /* allow parsing of SOL-format on the MIPLIB 2003 pages */
         continue;

      /* parse the line */
      nread = sscanf(buffer, "%s %s %s\n", varname, valuestring, objstring);
      if( nread < 2 )
      {
         SCIPerrorMessage("Invalid input line %d in solution file <%s>: <%s>.\n", lineno, fname, buffer);
         error = TRUE;
         break;
      }

      /* find the variable */
      var = SCIPfindVar(scip, varname);
      if( var == NULL )
      {
         if( !unknownvariablemessage )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "unknown variable <%s> in line %d of solution file <%s>\n",
               varname, lineno, fname);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  (further unknown variables are ignored)\n");
            unknownvariablemessage = TRUE;
         }
         continue;
      }

      /* cast the value */
      if( strncasecmp(valuestring, "inv", 3) == 0 )
         continue;
      else if( strncasecmp(valuestring, "+inf", 4) == 0 || strncasecmp(valuestring, "inf", 3) == 0 )
         value = SCIPinfinity(scip);
      else if( strncasecmp(valuestring, "-inf", 4) == 0 )
         value = -SCIPinfinity(scip);
      else
      {
         nread = sscanf(valuestring, "%lf", &value);
         if( nread != 1 )
         {
            SCIPerrorMessage("Invalid solution value <%s> for variable <%s> in line %d of solution file <%s>.\n",
               valuestring, varname, lineno, fname);
            error = TRUE;
            break;
         }
      }

      /* set the solution value of the variable, if not multiaggregated */
      if( SCIPisTransformed(scip) && SCIPvarGetStatus(SCIPvarGetProbvar(var)) == SCIP_VARSTATUS_MULTAGGR )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "ignored solution value for multiaggregated variable <%s>\n", SCIPvarGetName(var));
      }
      else
      {
         SCIP_RETCODE retcode;
         retcode = SCIPsetSolVal(scip, sol, var, value);

         if( retcode == SCIP_INVALIDDATA )
         {
            if( SCIPvarGetStatus(SCIPvarGetProbvar(var)) == SCIP_VARSTATUS_FIXED )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "ignored conflicting solution value for fixed variable <%s>\n",
                  SCIPvarGetName(var));
            }
            else
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "ignored solution value for multiaggregated variable <%s>\n",
                  SCIPvarGetName(var));
            }
         }
         else
         {
            SCIP_CALL( retcode );
         }
      }
   }

   /* close input file */
   SCIPfclose(file);

      return SCIP_OKAY;
}

/** checks if all fractional variables in the given solution are fractional. */
static
SCIP_RETCODE checkCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,
   SCIP_Bool             roundit,
   SCIP_Bool*            success
   )
{
   SCIP_VAR** vars;
   int nvars;
   int nfracvars;
   int nbinvars;
   int nintvars;
   int i;
   *success = TRUE;

   assert(NULL != sol);

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   /* proofs if the candidates are fractional */
   nfracvars = nbinvars + nintvars;
   for (i = 0; i < nfracvars; ++i)
   {
      SCIP_Real value;
      value = SCIPgetSolVal(scip, sol, vars[i]);
      if( !SCIPisFeasIntegral(scip, value) )
      {
         if( roundit )
         {
            SCIP_Real roundedvalue;
            roundedvalue = SCIPround(scip, value);
            SCIP_CALL(SCIPsetSolVal(scip, sol, vars[i], roundedvalue));
         }
         else
         {
            *success = FALSE;
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL,"All variables finished.\n");
            return SCIP_OKAY;
         }
      }
   }
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL,"All variables rounded.\n");
   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_HEUR*            heur,               /**< RENS heuristic structure                            */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;                         /* the original problem's number of variables      */
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}
/* put your local methods here, and declare them static */


/*
 * Callback methods of primal heuristic
 */


/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_HEURCOPY(heurCopyRepair)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of repair primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurCopyRepair NULL
#endif

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeRepair)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);
   SCIPfreeMemory(scip, &heurdata);

   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitRepair)
{  /*lint --e{715}*/

   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);

   heurdata->subiters = -1;
   heurdata->subnodes = -1;
   heurdata->subpresoltime = 0;

   heurdata->ninvalidvars = 0;
   heurdata->norvars = 0;
   heurdata->relinvalidvars = 0;
   heurdata->ninvalidcons = 0;
   heurdata->norcons = 0;
   heurdata->relinvalidcons = 0;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitRepair)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_Real time;
   SCIP_Real relvars;
   SCIP_Real relcons;
   char message[1024];
   int invalids;
   int ninvars;
   int ninvcons;
   int nvars;
   int ncons;
   int iterations;
   int nodes;
   int runs;

   heurdata = SCIPheurGetData(heur);
   invalids = heurdata->ninvalidvars+heurdata->ninvalidcons;
   ninvars = heurdata->ninvalidvars;
   ninvcons = heurdata->ninvalidcons;
   nvars = heurdata->norvars;
   ncons = heurdata->norcons;
   iterations = heurdata->subiters;
   nodes = heurdata->subnodes;
   time = heurdata->subpresoltime;
   runs = heurdata->runs;

   /* TODO remove this from heuristic data */
   heurdata->relinvalidvars = MAX((SCIP_Real)heurdata->norvars, 1.0);
   heurdata->relinvalidvars = heurdata->ninvalidvars/heurdata->relinvalidvars;
   heurdata->relinvalidcons = MAX((SCIP_Real)heurdata->norcons, 1.0);
   heurdata->relinvalidcons = heurdata->ninvalidcons/heurdata->relinvalidcons;

   relvars = heurdata->relinvalidvars;
   relcons = heurdata->relinvalidcons;


   /* Prints all statistic data for an user*/
   SCIPstatistic(
      sprintf(message, "Repair: \n\t total invalids: %d\n\n\t invalid variables: %d\n\t total variables: %d\n\t relative invalid variables: %.2f%%\n", invalids, ninvars, nvars, 100 * relvars);
      sprintf(message, "%s  \n\n\t invalid constraints: %d\n\t total constraints: %d\n\t relative invalid constraints: %.2f%%\n", message, ninvcons, ncons, 100* relcons);
      sprintf(message, "%s  \n\n\t iterations: %d\n\t nodes: %d\n\t number of runs: %d\n\t presolve time: %.2f s\n", message,iterations,nodes,runs,time);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, message);
   )
   return SCIP_OKAY;
}



/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_HEURINITSOL(heurInitsolRepair)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of repair primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitsolRepair NULL
#endif


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_HEUREXITSOL(heurExitsolRepair)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of repair primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitsolRepair NULL
#endif


/** execution method of primal heuristic.
 * Repair needs an incorrect solution, in which all variables are in their bound. */
static
SCIP_DECL_HEUREXEC(heurExecRepair)
{ /*lint --e{715}*/
   SCIP* subscip;
   SCIP_VAR** vars;
   SCIP_VAR** subvars = NULL;
   SCIP_ROW** rows;
   SCIP_SOL* sol;
   SCIP_SOL* subsol;
   /*SCIP_SOL* newsol;*/
   SCIP_HEURDATA* heurdata;
   /*char solfilename[1024];*/
   /*FILE* solfile;*/
   SCIP_RETCODE retcode = SCIP_OKAY;
   SCIP_Real timelimit;
   SCIP_Real memorylimit;
   char probname[SCIP_MAXSTRLEN];
   int i;
   int nbinvars;
   int nintvars;
   int nvars;
   int nrows;
   SCIP_Bool cutoff;
   SCIP_Bool success;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( !SCIPhasCurrentNodeLP(scip) )
      return SCIP_OKAY;

   if( !SCIPisLPConstructed(scip) )
   {
      SCIP_CALL(SCIPconstructLP(scip, &cutoff));

      if( cutoff )
         return SCIP_OKAY;
   }

   /* create zero solution */
   SCIP_CALL(SCIPcreateOrigSol(scip, &sol, heur));

   heurdata = SCIPheurGetData(heur);

   /* use read method to enter solution from a file */
   retcode = readSol(scip, sol, heurdata->filename);

   if( SCIP_NOFILE == retcode )
   {
      if( strcmp(heurdata->filename, DEFAULT_FILENAME) != 0 )
         SCIPwarningMessage(scip, "cannot open file <%s> for reading\n",
               heurdata->filename);

      goto TERMINATE;
   }
   else if( retcode != SCIP_OKAY )
   {
      goto TERMINATE;
   }
   SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "Repair: Solution file read.\n");

   /* checks the integrality of all discrete variable */
   SCIP_CALL( checkCands(scip, sol, heurdata->roundit, &success) );

   if( !success )
      goto TERMINATE;

   *result = SCIP_DIDNOTFIND;

   /* initializes the subscip */
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );
   SCIP_CALL( SCIPcopyParamSettings(scip, subscip) );

   /* get name of the original problem and add the string "_repairsub" */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_repairsub",
         SCIPgetProbName(scip));

   SCIP_CALL( SCIPcreateProb(subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPcreateSol(subscip, &subsol, heur) );

   /*Adds all original variables and adapts their bounds*/
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

   for( i = 0; i < nvars; ++i )
   {
      SCIP_CONS* cons;
      SCIP_Real orlb;
      SCIP_Real orub;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real value;
      SCIP_Real slack;
      SCIP_VARTYPE vartype;
      char varname[1024];
      char slackvarname[1024];
      char consvarname[1024];

      heurdata->norvars++;
      slack = 0.0;
      orlb = SCIPvarGetLbGlobal(vars[i]);
      orub = SCIPvarGetUbGlobal(vars[i]);
      value = SCIPgetSolVal(scip, sol, vars[i]);
      vartype = SCIPvarGetType(vars[i]);
      sprintf(varname, "sub_%s", SCIPvarGetName(vars[i]));

      /* if the value of x is out of bound, sets the scip to an correcting value*/
      if( SCIPisFeasLT(scip, value, orlb) )
      {
         lb = value;
         slack = orlb - value;
      }
      else
      {
         lb = orlb;
      }
      if( SCIPisFeasGT(scip, value, orub) )
      {
         ub = value;
         slack = orub - value;
      }
      else
      {
         ub = orub;
      }
      if( !SCIPisFeasZero(scip, slack) && SCIP_VARTYPE_BINARY == vartype )
         vartype = SCIP_VARTYPE_INTEGER;

      /* Adds the sub representing variable to the subscip. */
      SCIP_CALL(SCIPcreateVarBasic(subscip, &subvars[i], varname, lb, ub, 0, vartype));
      SCIP_CALL(SCIPaddVar(subscip, subvars[i]));
      SCIP_CALL( SCIPsetSolVal(subscip, subsol, subvars[i], value) );

      /* if necessary adds an constraint to represent the original bounds of x.*/
      if( !SCIPisFeasEQ(scip, slack, 0.0) )
      {
         SCIP_VAR* newvar;
         heurdata->ninvalidvars++;
         sprintf(slackvarname, "artificialslack_%s", SCIPvarGetName(vars[i]));
         sprintf(consvarname, "boundcons_%s", SCIPvarGetName(vars[i]));
         SCIP_CALL( SCIPcreateVarBasic(subscip, &newvar, slackvarname, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY));
         SCIP_CALL(SCIPaddVar(subscip, newvar));
         SCIP_CALL( SCIPsetSolVal(subscip, subsol, newvar, 1.0) );
         SCIP_CALL( SCIPcreateConsBasicVarbound(subscip, &cons, consvarname, subvars[i], newvar, slack, orlb, orub) );
         SCIP_CALL( SCIPaddCons(subscip, cons) );
         SCIP_CALL( SCIPreleaseVar(subscip, &newvar) );
         SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
      }
   }

   /* check solution for feasibility regarding the LP rows (SCIPgetRowSolActivity()) */
   rows = SCIPgetLPRows(scip);
   nrows = SCIPgetNLPRows(scip);

   for (i = 0; i < nrows; ++i)
   {
      SCIP_COL** cols;
      SCIP_VAR** consvars;
      SCIP_CONS* cons;
      SCIP_VAR* newvar;
      SCIP_Real* vals;
      SCIP_Real slack;
      char varname[1024];
      int nnonz;
      SCIP_Real constant;
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real rowsolact;
      int j;

      heurdata->norcons++;

      /* gets the values to check the constraint */
      constant = SCIProwGetConstant(rows[i]);
      lhs = SCIPisInfinity(scip, -SCIProwGetLhs(rows[i])) ?
            SCIProwGetLhs(rows[i]) : SCIProwGetLhs(rows[i]) - constant;
      rhs = SCIPisInfinity(scip, SCIProwGetRhs(rows[i])) ?
            SCIProwGetRhs(rows[i]) : SCIProwGetRhs(rows[i]) - constant;
      rowsolact = SCIPgetRowSolActivity(scip, rows[i], sol) - constant;
      vals = SCIProwGetVals(rows[i]);

      assert(SCIPisFeasLE(scip, lhs, rhs));

      nnonz = SCIProwGetNNonz(rows[i]);
      cols = SCIProwGetCols(rows[i]);
      SCIP_CALL(SCIPallocBufferArray(scip, &consvars, nnonz));

      /* translate all variables from the original scip to the subscip with subscip variables. */
      for (j = 0; j < nnonz; ++j)
      {
         int pos;

         pos = SCIPvarGetProbindex(SCIPcolGetVar(cols[j]));
         consvars[j] = subvars[pos];
      }

      SCIP_CALL(
            SCIPcreateConsBasicLinear(subscip, &cons, SCIProwGetName(rows[i]),
                  nnonz, consvars, vals, lhs, rhs));

      /* Sets the slack if its necessary */
      if( SCIPisFeasLT(scip, rowsolact, lhs) )
      {
         slack = lhs - rowsolact;
      }
      else if( SCIPisFeasGT(scip, rowsolact, rhs) )
      {
         slack = rhs - rowsolact;
      }
      else
      {
         slack = 0;
      }

      /*if necessary adds an new artificial slack variable*/
      if( !SCIPisFeasEQ(subscip, slack, 0.0) )
      {
         heurdata->ninvalidcons++;
         sprintf(varname, "artificialslack_%s", SCIProwGetName(rows[i]));
         SCIP_CALL( SCIPcreateVarBasic(subscip, &newvar, varname, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
         SCIP_CALL( SCIPaddVar(subscip, newvar) );
         SCIP_CALL( SCIPsetSolVal(subscip, subsol, newvar, 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(subscip, cons, newvar, slack) );
         SCIP_CALL( SCIPreleaseVar(subscip, &newvar) );
      }
      /*Adds the Constraint and release it.*/
      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
      SCIPfreeBufferArray(scip, &consvars);
   }

   /*Adds the given solution to the subscip. */
   SCIP_CALL( SCIPaddSolFree(subscip, &subsol, &success) );

   if( !success )
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Repair was not good enough.\n");

   /* check whether there is enough time and memory left */
   SCIP_CALL(SCIPgetRealParam(scip, "limits/time", &timelimit));
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL(SCIPgetRealParam(scip, "limits/memory", &memorylimit));

   /* substract the memory already used by the main SCIP and the estimated memory usage of external software */
   if( !SCIPisInfinity(scip, memorylimit) )
   {
      memorylimit -= SCIPgetMemUsed(scip) / 1048576.0;
      memorylimit -= SCIPgetMemExternEstim(scip) / 1048576.0;
   }

   /* abort if no time is left or not enough memory to create a copy of SCIP, including external memory usage */
   if( timelimit <= 0.0
         || memorylimit <= 2.0 * SCIPgetMemExternEstim(scip) / 1048576.0 )
      goto TERMINATE;

   /* set limits for the subproblem */
   /*heurdata->nodelimit = nstallnodes;
   SCIP_CALL(SCIPsetLongintParam(subscip, "limits/nodes", nstallnodes));*/
   SCIP_CALL(SCIPsetRealParam(subscip, "limits/time", timelimit));
   SCIP_CALL(SCIPsetRealParam(subscip, "limits/memory", memorylimit));

   /* forbid recursive call of heuristics and separators solving sub-SCIPs */
   SCIP_CALL(SCIPsetSubscipsOff(subscip, TRUE));

   /* disable output to console */
   SCIP_CALL(SCIPsetIntParam(subscip, "display/verblevel", 0));

#ifdef SCIP_DEBUG
   /* for debugging Repair, enable MIP output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000000) );
#endif

   /* presolve the subproblem */
   retcode = SCIPpresolve(subscip);

   /* errors in solving the subproblem should not kill the overall solving process;
    * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL(retcode);
#endif
      SCIPwarningMessage(scip,
            "Error while presolving subproblem in REPAIR heuristic; sub-SCIP terminated with code <%d>\n",
            retcode);

      /* free */
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL(SCIPfree(&subscip));
      return SCIP_OKAY;
   }

   SCIPwriteOrigProblem(subscip,"stein27_inf_friday.cip",NULL,FALSE);
   /* solve the subproblem */
   retcode = SCIPsolve(subscip);

   if( retcode != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL(retcode);
#endif
      SCIPwarningMessage(scip,
            "Error while solving subproblem in REPAIR heuristic; sub-SCIP terminated with code <%d>\n",
            retcode);

      /* free */
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL(SCIPfree(&subscip));
      return SCIP_OKAY;
   }

   /* print solving statistics of subproblem if we are in SCIP's debug mode */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );
   /*sprintf(solfilename, "%s.rsol", heurdata->filename);
   solfile = fopen(solfilename, "w");
   SCIP_CALL(SCIPprintBestSol(subscip, solfile, 0));

   fclose(solfile);*/

   assert(SCIPgetNSols(subscip) > 0);
   if( SCIPisFeasZero(scip, SCIPgetPrimalbound(subscip)) )
   {
      SCIP_CALL(createNewSol(scip, subscip, subvars, heur, SCIPgetBestSol(subscip), &success));
      /*assert(success || SCIPgetNSols(scip) > 0);*/

      if( success )
         *result = SCIP_FOUNDSOL;
   }

   heurdata->subiters = SCIPgetNLPIterations(subscip);
   heurdata->subnodes = SCIPgetNTotalNodes(subscip);
   heurdata->subpresoltime = SCIPgetPresolvingTime(subscip);
   heurdata->runs = SCIPgetNRuns(subscip);

   /* terminates the solving process  */
TERMINATE:
   SCIPfreeSol(scip, &sol);
   /*SCIPfreeSol(subscip, &subsol);*/
   SCIPfreeBufferArray(scip, &subvars);
   SCIPfree(&subscip);
   return retcode;
}



/*
 * primal heuristic specific interface methods
 */

/** creates the repair primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRepair(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create repair primal heuristic data */
   heurdata = NULL;

   SCIP_CALL( SCIPallocMemory(scip ,&heurdata) );

   heur = NULL;

   /* include primal heuristic */
   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecRepair, heurdata) );

   assert(heur != NULL);
   assert(heurdata != NULL);


   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyRepair) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeRepair) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitRepair) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitRepair) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolRepair) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolRepair) );

   /* add repair primal heuristic parameters */

   heurdata->filename = NULL;
   /* add string parameter for filename containing a solution */
   SCIP_CALL( SCIPaddStringParam(scip, "heuristics/"HEUR_NAME"/filename",
         "file name of a solution to be used as infeasible starting point, [-] if not available",
         &heurdata->filename, FALSE, DEFAULT_FILENAME, NULL, NULL) );

   /* add bool parameter for decision how to deal with unfractional cands */
   SCIP_CALL(
         SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/roundit",
               "True : fractional variables which are not fractional in the given solution are rounded, "
               "FALSE : solvingprocess of this heuristic is stoped. ",
               &heurdata->roundit, FALSE, DEFAULT_ROUNDIT, NULL, NULL));

   return SCIP_OKAY;
}

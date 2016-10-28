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
 * @author Gregor Hendel
 * @author Thomas Nagel
 *
 */

/* author bzfhende: TODO write a short description of what the heuristic does (modify heur_repair.h accordingly) */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define SCIP_STATISTIC
#include <assert.h>
#include <string.h>

#include <stdio.h>
#include <stdlib.h>

#include <libgen.h>

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
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_FILENAME      "-" /**< file name of a solution to be used as infeasible starting point */
#define DEFAULT_ROUNDIT       TRUE /**< if it is True : fractional variables which are not fractional in the given
                                    *    solution are rounded, if it is FALSE : solving process of this heuristic is stopped */
#define DEFAULT_USEOBJFACTOR  FALSE /**< should a scaled objective function for original variables be used in repair subproblem? */
#define DEFAULT_USEVARFIX     TRUE  /**< should variable fixings be used in repair subproblem? */
#define DEFAULT_USESLACKVARS  FALSE /**< should slack variables be used in repair subproblem? */
#define DEFAULT_ALPHA         2.0

/*
 * Data structures
 */


/** primal heuristic data */
struct SCIP_HeurData
{
   char*                 filename;           /**< file name of a solution to be used as infeasible starting point */
   SCIP_Bool             roundit;            /**< if it is True : fractional variables which are not fractional in the given
                                              *    solution are rounded, if it is FALSE : solvingprocess of this heuristic is stoped */
   SCIP_Bool             useobjfactor;       /**< should a scaled objective function for original variables be used in repair subproblem? */
   SCIP_Bool             usevarfix;          /**< should variable fixings be used in repair subproblem? */
   SCIP_Bool             useslackvars;       /**< should slack variables be used in repair subproblem? */


   int                   subnodes;           /**< number of nodes which were necessary to solve the subscip */
   int                   subiters;           /**< contains total number of iterations used in primal and dual simplex and barrier algorithm to solve the subscip */
   SCIP_Real             subpresoltime;      /**< time for presolving the subscip */
   int                   runs;               /**< number of branch and bound runs performed to solve the subscip */

   int                   nviolatedvars;      /**< number of violated vars in the given solution */
   int                   norvars;            /**< number of all vars in the given problem */
   SCIP_Real             relviolatedvars;    /**< relative number of violated vars */
   int                   nviolatedcons;      /**< number of violated cons in the given solution */
   int                   norcons;            /**< number of all cons in the given problem */
   SCIP_Real             relviolatedcons;    /**< relative number of violated cons */

   int                   nvarfixed;          /**< number of all variables fixed in the sub problem */
   SCIP_Real             relvarfixed;        /**< relative number of fixed vars */

   SCIP_Real             orsolval;           /**< value of the solution find by repair, in the original Problem*/
   SCIP_Real             improovedoldsol;    /**< value of the given sol sfter beeing improoved by SCIP */
   SCIP_SOL*             infsol;             /**< infeasible solution to start with */
   SCIP_Real             alpha;              /**< */
};


/** TODO GENERAL TODOS:
 *  - implement a minfixingrate-parameter
 *  - concentrate on the binary and integer variables for the fixing.
 *  - fold the methods extendedProgram and variableFixings into one that does both (and calls subroutines where necessary)
 */

/*
 * Local methods
 */

/** reads a given SCIP solution file, problem has to be transformed in advance */
/* author bzfhende: TODO this is okay, but in the meantime when you merge master back in, there is a new function SCIPreadSol() */
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

   if( strcmp(fname, DEFAULT_FILENAME) == 0 )
   {
      /*SCIPfreeSol(scip, &sol);*/
      SCIPlinkLPSol(scip, sol);
      return SCIP_OKAY;
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
      return SCIP_NOFILE;
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
         SCIPerrorMessage("violated input line %d in solution file <%s>: <%s>.\n", lineno, fname, buffer);
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
            SCIPerrorMessage("violated solution value <%s> for variable <%s> in line %d of solution file <%s>.\n",
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

/**
 * computes a factor to norm the original objectiv uper bound to 1.
 */
static
SCIP_RETCODE getObjectiveFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,               /**< SCIP data structure */
   SCIP_Real*            factor,             /**< SCIP_Real to save the factor for the old objective function*/
   SCIP_Bool*            success             /**< SCIP_Bool: Is the factor real?*/
   )
{
   SCIP_VAR** vars;
   SCIP_Real lprelaxobj;
   SCIP_Real upperbound;
   SCIP_Real objoffset;
   int nvars;
   int i;

   *success = TRUE;
   *factor = 0.0;
   upperbound = 0.0;

   /* IF there is a solution already found DO:*/
   lprelaxobj = SCIPgetLowerbound(scip);

   if( SCIPisInfinity(scip, -lprelaxobj) )
   {
      *factor = 0.0;
      return SCIP_OKAY;
   }



   if( !SCIPisInfinity(scip, SCIPgetUpperbound(scip)) )
   {
      upperbound = SCIPgetUpperbound(scip);
   }
   else
   {
      SCIP_CALL(SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL));

      for (i = 0; i < nvars; ++i)
      {
         upperbound = SCIPvarGetObj(vars[i]);
         if( SCIPisInfinity(scip, upperbound) || SCIPisInfinity(scip, -upperbound) )
         {
            /* TODO fancy diving function */
            *factor = 1 / SCIP_REAL_MAX;
            return SCIP_OKAY;
         }
         else if( SCIPisZero(scip, upperbound) )
         {
            continue;
         }
         /* author bzfhende: TODO I think you should treat infinite bounds explicitly. By the way, it seems that you should
          * exchange the use of lower and upper bounds for obtaining an "antipseudo"-solution-like upperbound */

         else if( SCIPisGT(scip, 0.0, upperbound) )
         {
            *factor += upperbound * SCIPvarGetLbGlobal(vars[i]);
         }
         else
         {
            *factor += upperbound * SCIPvarGetUbGlobal(vars[i]);
         }
      }
   }

   /* Ending-sequence */
   *factor = upperbound - lprelaxobj;
   if( !SCIPisZero(scip, *factor) )
   {
      *factor = 1 / *factor;
   }

   /* sets an offset, which guarantee positive objective values */
   objoffset = -lprelaxobj * (*factor);
   SCIP_CALL( SCIPaddOrigObjoffset(subscip, -objoffset) );

   return SCIP_OKAY;
}

/** */
static
SCIP_Real getPotentialContributed(
   SCIP*             scip,
   SCIP_SOL*         sol,
   SCIP_VAR*         var,
   SCIP_Real         coeffiant,
   int               sgn
   )
{
   if( 0 > sgn * coeffiant )
   {
      if(SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)))
      {
         return SCIPinfinity(scip);
      }
      return coeffiant * (SCIPgetSolVal(scip, sol, var) - SCIPvarGetLbGlobal(var));
   }
   else
   {
      if(SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)))
      {
         return -SCIPinfinity(scip);
      }
      return coeffiant * (SCIPgetSolVal(scip, sol, var) - SCIPvarGetUbGlobal(var));
   }
}





/**
 * tries to fix a variable, if a it is possible, the potentials of rows is adapted.
 */
static
SCIP_Bool tryFixVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< SCIP data structure */
   SCIP_Real*            potential,          /**< Array with all potential values */
   SCIP_Real*            slack,              /**< Array with all potential values */
   SCIP_VAR*             var,                /**< Variable to be fixed? */
   int*                  inftycounter,       /**< counters how many variables have an infinity potential in a row */
   SCIP_HEURDATA*        heurdata            /**< repairs heuristic data */
   )
{
   SCIP_ROW** rows;
   SCIP_COL* col;
   SCIP_Real* vals;
   SCIP_Real alpha = heurdata->alpha;
   int nrows;
   int i;
   int sgn;
   int rowindex;

   if( SCIPisFeasLT(scip, SCIPgetSolVal(scip, sol, var), SCIPvarGetLbGlobal(var)) )
   {
      return FALSE;
   }
   if( SCIPisFeasGT(scip, SCIPgetSolVal(scip, sol, var), SCIPvarGetUbGlobal(var)) )
   {
      return FALSE;
   }

   col = SCIPvarGetCol(var);
   rows = SCIPcolGetRows(col);
   nrows = SCIPcolGetNLPNonz(col);
   vals = SCIPcolGetVals(col);

   /*iterate over nonzero rows*/
   for( i = 0; i < nrows; ++i)
   {
      SCIP_Real contribution;
      rowindex = SCIProwGetLPPos(rows[i]);
      assert(rowindex >= 0);
      if( 0 > rowindex )
      {
         continue;
      }
      if( 0 == slack[rowindex] )
      {
         continue;
      }
      sgn = 1;
      if( 0 > slack[rowindex])
      {
         sgn = -1;
      }

      contribution = getPotentialContributed(scip, sol, var, vals[i], sgn);
      if( !SCIPisInfinity(scip, REALABS(contribution)) )
      {
         potential[rowindex] -= contribution;
      }
      else
      {
         inftycounter[rowindex]--;
      }
      /* author bzfhende: TODO I think this needs to be a -sign because fixing this variable reduces the potential */
      /*if( 0 > sgn * vals[i] )
      {
         assert(!SCIPisInfinity(scip, SCIPvarGetLbGlobal(var)));
         potential[rowindex] += vals[i] * (SCIPgetSolVal(scip, sol, var) - SCIPvarGetLbGlobal(var));
      }
      else
      {
         assert(!SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)));
         potential[rowindex] += vals[i] * (SCIPgetSolVal(scip, sol, var) - SCIPvarGetUbGlobal(var));
      }*/
      /* todo comment negative slack means violated right hand side, but in this case, the potential is always positive because of the single contributions you compute */
      assert(0 <= inftycounter[rowindex]);
      if( 0 == inftycounter[rowindex] && REALABS(potential[rowindex]) < alpha * REALABS(slack[rowindex]) )
      {
         /* author bzfhende: TODO this is hard to read here */

         for( ; i >= 0; --i )
         {
            sgn = 1;
            /*rowindex = SCIProwGetIndex(rows[i]);*/
            if( 0 == slack[rowindex] )
            {
               continue;
            }
            rowindex = SCIProwGetLPPos(rows[i]);
            if( 0 > slack[rowindex])
            {
               sgn = -1;
            }
            contribution = getPotentialContributed(scip, sol, var, vals[i], sgn);
            if( !SCIPisInfinity(scip, REALABS(contribution)) )
            {
               potential[rowindex] += contribution;
            }
            else
            {
               inftycounter[rowindex]++;
            }
         }
         return FALSE;
      }
   }

   return TRUE;
}
/** checks if all integral variables in the given solution are integral. */
static
SCIP_RETCODE checkCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /* author bzfhende: TODO missing documentation */
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

   /* test if the candidates are fractional */
   nfracvars = nbinvars + nintvars;
   for( i = 0; i < nfracvars; ++i)
   {
      SCIP_Real value;
      value = SCIPgetSolVal(scip, sol, vars[i]);
      if( !SCIPisFeasIntegral(scip, value) )
      {
         /* author bzfhende: TODO use more sophisticated rounding that involves SCIPvarGetNLocksUp/Down() */
         if( roundit )
         {
            SCIP_Real roundedvalue;

            if( SCIPvarGetNLocksUp(vars[i]) > SCIPvarGetNLocksDown(vars[i]) )
            {
               roundedvalue = SCIPceil(scip, value - 1.0);
            }
            else
            {
               roundedvalue = SCIPfloor(scip, value + 1.0);
            }

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
   SCIP_HEURDATA*         heurdata,           /**< repairs heurdata                                    */
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
   SCIP_Real  valuetmp;                      /* var to save the original value of the solution temporally*/

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
   valuetmp = SCIPgetSolOrigObj(scip, newsol);
   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, TRUE, TRUE, TRUE, success) );

   if( success )
      heurdata->orsolval = valuetmp;

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/** reads a given SCIP solution file, problem has to be transformed in advance */
static
SCIP_RETCODE extendedProgramm(
   SCIP*                 scip,             /* author bzfhende: TODO missing documentation */
   SCIP_HEUR*            heur,
   SCIP_RESULT*          result
   )
{
   SCIP* subscip = NULL;
   SCIP_VAR** vars;
   SCIP_VAR** subvars = NULL;
   SCIP_ROW** rows;
   SCIP_SOL* sol;
   SCIP_SOL* subsol;
   SCIP_HEURDATA* heurdata;
   SCIP_RETCODE retcode = SCIP_OKAY;
   SCIP_Real timelimit;
   SCIP_Real memorylimit;
   SCIP_Real factor;
   char probname[SCIP_MAXSTRLEN];
   int i;
   int nbinvars;
   int nintvars;
   int nvars;
   int nrows;
   SCIP_Bool success;

   heurdata = SCIPheurGetData(heur);
   sol = heurdata->infsol;

   /* initializes the subscip */
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );
   SCIP_CALL( SCIPcopyParamSettings(scip, subscip) );

   /* use inference branching */
   if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
    }

   /* get name of the original problem and add the string "_repairsub" */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_repairsub",
         SCIPgetProbName(scip));

   SCIP_CALL( SCIPcreateProb(subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPcreateSol(subscip, &subsol, heur) );

   /* gets all original variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

   /* author bzfhende: TODO use new debug message function SCIPdebugMsg(scip, ...) */
   SCIPdebugMessage("\n\n Calling objective factor calculation \n\n");
   if( heurdata->useobjfactor )
   {
      SCIP_CALL( getObjectiveFactor(scip, subscip, &factor, &success) );
   }
   else
   {
      factor = 0.0;
   }

   SCIPdebugMessage("\n\n  Calculated objective factor %15.9f\n\n",factor);

   if( !success )
   {
      factor = 0.0;
      SCIPwarningMessage(scip,"Repair could not set a factor!\n");
   }

   /* Adds all original variables and adapts their bounds */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CONS* cons;
      SCIP_Real lborig;
      SCIP_Real uborig;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real objval;
      SCIP_Real value;
      SCIP_Real slack;
      SCIP_VARTYPE vartype;
      char varname[1024];

      heurdata->norvars++;
      slack = 0.0;
      lborig = SCIPvarGetLbGlobal(vars[i]);
      uborig = SCIPvarGetUbGlobal(vars[i]);
      value = SCIPgetSolVal(scip, sol, vars[i]);
      vartype = SCIPvarGetType(vars[i]);
      objval = SCIPvarGetObj(vars[i])*factor;

      if( SCIPisZero(scip, objval) )
      {
         objval = 0.0;
      }

      sprintf(varname, "sub_%s", SCIPvarGetName(vars[i]));

      /* if the value of x is lower than the lower bound, sets the slack to an correcting value */
      if( SCIPisFeasLT(scip, value, lborig) )
      {
         lb = value;
         slack = lborig - value;
      }
      else
      {
         lb = lborig;
      }

      /* if the value of x is bigger than the upper bound, sets the slack to an correcting value*/
      if( SCIPisFeasGT(scip, value, uborig) )
      {
         ub = value;
         slack = uborig - value;
      }
      else
      {
         ub = uborig;
      }

      /* if an binary variable is out of bound, generalize it to an integer variable */
      if( !SCIPisFeasZero(scip, slack) && SCIP_VARTYPE_BINARY == vartype )
         vartype = SCIP_VARTYPE_INTEGER;

      /* Adds the sub representing variable to the subscip. */
      SCIP_CALL( SCIPcreateVarBasic(subscip, &subvars[i], varname, lb, ub, objval, vartype) );
      SCIP_CALL( SCIPaddVar(subscip, subvars[i]) );
      SCIP_CALL( SCIPsetSolVal(subscip, subsol, subvars[i], value) );

      /* if necessary adds an constraint to represent the original bounds of x.*/
      if( !SCIPisFeasEQ(scip, slack, 0.0) )
      {
         SCIP_VAR* newvar;
         char consvarname[1024];
         char slackvarname[1024];

         sprintf(slackvarname, "artificialslack_%s", SCIPvarGetName(vars[i]));
         sprintf(consvarname, "boundcons_%s", SCIPvarGetName(vars[i]));

         /* initialize and add an artificial slack variable */
         /* author bzfhende: TODO didn't we want to change the variable type to continuous? */
         SCIP_CALL( SCIPcreateVarBasic(subscip, &newvar, slackvarname, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY));
         SCIP_CALL(SCIPaddVar(subscip, newvar));

         /* set the value of the slack variable to 1 to make subsolution feasible */
         SCIP_CALL( SCIPsetSolVal(subscip, subsol, newvar, 1.0) );

         /* adds a linear constraint to represent the old bounds */
         SCIP_CALL( SCIPcreateConsBasicVarbound(subscip, &cons, consvarname, subvars[i], newvar, slack, lborig, uborig) );
         SCIP_CALL( SCIPaddCons(subscip, cons) );
         SCIP_CALL( SCIPreleaseVar(subscip, &newvar) );
         SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

         /* increases the counter for violated vars */
         heurdata->nviolatedvars++;
      }
   }

   /* check solution for feasibility regarding the LP rows (SCIPgetRowSolActivity()) */
   rows = SCIPgetLPRows(scip);
   nrows = SCIPgetNLPRows(scip);

   /* Adds all original constraints and adapts their bounds */
   for( i = 0; i < nrows; ++i )
   {
      SCIP_COL** cols;
      SCIP_VAR** consvars;
      SCIP_CONS* cons;
      SCIP_VAR* newvar;
      SCIP_Real* vals;
      SCIP_Real slack;
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
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nnonz) );

      /* translate all variables from the original scip to the subscip with subscip variables. */
      for( j = 0; j < nnonz; ++j )
      {
         int pos;
         /* author bzfhende: TODO here, you should avoid columns with position -1 */
         pos = SCIPvarGetProbindex(SCIPcolGetVar(cols[j]));
         consvars[j] = subvars[pos];
      }

      /* create a new linear constraint, representing the old one */
      SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &cons, SCIProwGetName(rows[i]),
               nnonz, consvars, vals, lhs, rhs) );

      /* sets the slack if its necessary: negative slack for right hand side violation, positive slack for left hand side violation */
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
         char varname[1024];
         sprintf(varname, "artificialslack_%s", SCIProwGetName(rows[i]));
         /* author bzfhende: TODO vartype should be continuous, right? */
         SCIP_CALL( SCIPcreateVarBasic(subscip, &newvar, varname, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
         SCIP_CALL( SCIPaddVar(subscip, newvar) );
         SCIP_CALL( SCIPsetSolVal(subscip, subsol, newvar, 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(subscip, cons, newvar, slack) );
         SCIP_CALL( SCIPreleaseVar(subscip, &newvar) );
         heurdata->nviolatedcons++;
      }
      /*Adds the Constraint and release it.*/
      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
      SCIPfreeBufferArray(scip, &consvars);
   }

   /*if( 0 >= heurdata->nviolatedcons && 0 >= heurdata->nviolatedvars )
   {
      SCIP_CALL( SCIPaddSolFree(scip, &subsol, &success) );
      if( !success )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Given repair-solution was feasible, but not good enough.\n");
      }
      else
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Given repair-solution was feasible.\n");
      }
      *result = SCIP_FOUNDSOL;
      return SCIP_OKAY;
   }*/

   /*Adds the given solution to the subscip. */
   heurdata->improovedoldsol = SCIPgetSolOrigObj(subscip, subsol);
   SCIP_CALL( SCIPaddSolFree(subscip, &subsol, &success) );

   /* author bzfhende: TODO make this a meaningful debug(!) message, currently, it states 'repair was not good enough' where it should
    * say something like 'subsolution could not be added' */
   if( !success )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Repair was not good enough.\n");
   }

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
   {
      SCIP_Longint nnodes;
      SCIP_CALL(SCIPgetLongintParam(scip, "limits/nodes", &nnodes));
      nnodes = (SCIP_Longint)(0.1 * nnodes);
      SCIP_CALL(SCIPsetRealParam(subscip, "limits/time", timelimit));
      SCIP_CALL(SCIPsetRealParam(subscip, "limits/memory", memorylimit));
      SCIP_CALL(SCIPsetObjlimit(subscip,1.0));
      SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nnodes) );
      SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", MAX(10, nnodes/10)) );
      SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", 3) );
   }
   /*SCIP_CALL(SCIPsetLongintParam(subscip, "limits/nodes", 1));*/

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

   /* solve the subproblem */
   retcode = SCIPsolve(subscip);
   heurdata->improovedoldsol = SCIPgetSolOrigObj(subscip, SCIPgetBestSol(subscip));

   if( retcode != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL(retcode);
#endif
      SCIPwarningMessage(scip,
            "Error while solving subproblem in REPAIR heuristic; sub-SCIP terminated with code <%d>\n",
            retcode);

      /* free */
      for( i = 0; i < nvars; ++i )
      {
         SCIP_CALL( SCIPreleaseVar(subscip, &subvars[i]) );
      }
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL(SCIPfree(&subscip));
      return SCIP_OKAY;
   }

   /* print solving statistics of subproblem if we are in SCIP's debug mode */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );

   assert(SCIPgetNSols(subscip) > 0);
   /* author bzfhende: TODO in the case of continuous slack variables, you would need to check the feasibility of the
    * returned solution by checking all slack variables
    */
   if( SCIPisLT(scip, SCIPgetPrimalbound(subscip), 1.0) )
   {
      SCIP_CALL(createNewSol(heurdata, scip, subscip, subvars, heur, SCIPgetBestSol(subscip), &success));

      if( success )
         *result = SCIP_FOUNDSOL;
   }

   /*if( SCIPisFeasZero(scip, SCIPgetPrimalbound(subscip)) )
   {
      SCIP_CALL(createNewSol(heurdata, scip, subscip, subvars, heur, SCIPgetBestSol(subscip), &success));

      if( success )
         *result = SCIP_FOUNDSOL;
   }
   */

   heurdata->subiters = SCIPgetNLPIterations(subscip);
   heurdata->subnodes = SCIPgetNTotalNodes(subscip);
   heurdata->subpresoltime = SCIPgetPresolvingTime(subscip);
   heurdata->runs = SCIPgetNRuns(subscip);

   /* terminates the solving process  */
TERMINATE:
   SCIPfreeSol(scip, &sol);
   /*SCIPfreeSol(subscip, &subsol);*/
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(subscip, &subvars[i]) );
   }
   SCIPfreeBufferArrayNull(scip, &subvars);
   if( NULL != subscip )
   {
      SCIP_CALL( SCIPfree(&subscip) );
   }

   return SCIP_OKAY;
}

/**
 *
 */
/* author bzfhende: TODO adapt method head according to the SCIP style guidelines */

static
SCIP_RETCODE varFixings(SCIP* scip, SCIP_HEUR* heur, SCIP_RESULT* result)
{
   SCIP* subscip = NULL;
   SCIP_VAR** vars = NULL;
   SCIP_VAR** subvars = NULL;
   SCIP_ROW** rows = NULL;
   SCIP_CONS** subcons = NULL;
   int* nviolatedrows = NULL;
   int* permutation = NULL;
   int* inftycounter = NULL;
   SCIP_SOL* sol = NULL;
   SCIP_SOL* subsol = NULL;
   SCIP_HEURDATA* heurdata = NULL;
   SCIP_Real* potential = NULL;
   SCIP_Real* slack = NULL;
   SCIP_RETCODE retcode = SCIP_OKAY;
   SCIP_Real timelimit;
   SCIP_Real memorylimit;
   SCIP_Real factor;
   char probname[SCIP_MAXSTRLEN];
   int i;
   int nbinvars;
   int nintvars;
   int nvars;
   int nrows;
   SCIP_Bool success;

   heurdata = SCIPheurGetData(heur);
   sol = heurdata->infsol;

   /* initializes the subscip */
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );
   SCIP_CALL( SCIPcopyParamSettings(scip, subscip) );

   /* use inference branching */
   if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
    }

   /* get name of the original problem and add the string "_repairsub" */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_repairsub",
         SCIPgetProbName(scip));

   SCIP_CALL( SCIPcreateProb(subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* get name of the original problem and add the string "_repairsub" */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_repairfixsub_fix",
         SCIPgetProbName(scip));

   SCIP_CALL( SCIPcreateSol(subscip, &subsol, heur) );

   /* Gets all original variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nviolatedrows, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &permutation, nvars) );

   SCIPdebugMessage("\n\n Calling objective factor calculation \n\n");
      if( heurdata->useobjfactor )
      {
         SCIP_CALL( getObjectiveFactor(scip, subscip, &factor, &success) );
      }
      else
      {
         factor = 0.0;
      }

   /* Adds all original variables */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real objval;
      SCIP_Real value;
      SCIP_VARTYPE vartype;
      char varname[1024];

      heurdata->norvars++;
      lb = SCIPvarGetLbGlobal(vars[i]);
      ub = SCIPvarGetUbGlobal(vars[i]);
      value = SCIPgetSolVal(scip, sol, vars[i]);
      vartype = SCIPvarGetType(vars[i]);

      nviolatedrows[i]=0;

      sprintf(varname, "sub_%s", SCIPvarGetName(vars[i]));

      objval = SCIPvarGetObj(vars[i]);
      /* Adds the sub representing variable to the subscip. */
      SCIP_CALL( SCIPcreateVarBasic(subscip, &subvars[i], varname, lb, ub, objval, vartype) );
      SCIP_CALL( SCIPaddVar(subscip, subvars[i]) );
      SCIP_CALL( SCIPsetSolVal(subscip, subsol, subvars[i], value) );

      /* author bzfhende: TODO whenever this happens, you do not need to add the solution because it violates a bound */
      if( SCIPisFeasLT(scip, value, lb) || SCIPisFeasGT(scip,value,ub) )
      {
         heurdata->nviolatedvars++;
      }
   }

   /* check solution for feasibility regarding the LP rows (SCIPgetRowSolActivity()) */
   rows = SCIPgetLPRows(scip);
   nrows = SCIPgetNLPRows(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &potential, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &slack, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subcons, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inftycounter, nrows) );
   /* Adds all original constraints and computes potentials and slacks */
   for (i = 0; i < nrows; ++i)
   {
      SCIP_COL** cols;
      SCIP_VAR** consvars;
      SCIP_Real* vals;
      SCIP_Real constant;
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real rowsolact;
      int nnonz;
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
      potential[i] = 0.0;
      inftycounter[i] = 0;

      assert(SCIPisFeasLE(scip, lhs, rhs));

      nnonz = SCIProwGetNNonz(rows[i]);
      cols = SCIProwGetCols(rows[i]);
      SCIP_CALL(SCIPallocBufferArray(subscip, &consvars, nnonz));

      /* sets the slack if its necessary */
      if( SCIPisFeasLT(scip, rowsolact, lhs) )
      {
         slack[i] = lhs - rowsolact;
         heurdata->nviolatedcons++;
      }
      else if( SCIPisFeasGT(scip, rowsolact, rhs) )
      {
         slack[i] = rhs - rowsolact;
         heurdata->nviolatedcons++;
      }
      else
      {
         slack[i] = 0.0;
      }

      /* translate all variables from the original scip to the subscip with subscip variables. */
      for( j = 0; j < nnonz; ++j )
      {
         SCIP_Real contribution;
         int pos;
         int sgn = 1;

         /* negative slack represents a right hand side violation */
         if( 0 > slack[i] )
         {
            assert(!SCIPisInfinity(scip, rhs));
            sgn = -1;
         }

         pos = SCIPvarGetProbindex(SCIPcolGetVar(cols[j]));
         consvars[j] = subvars[pos];

         /* compute potentials */
         /*if( 0 > sgn * vals[j] )
         {
            assert(!SCIPisInfinity(scip, SCIPvarGetLbGlobal(vars[pos])));
            potential[i] += vals[j] * (SCIPvarGetLbGlobal(vars[pos]) - SCIPgetSolVal(scip, sol, vars[pos]));
         }
         else
         {
            assert(!SCIPisInfinity(scip, SCIPvarGetUbGlobal(vars[pos])));
            potential[i] += vals[j] * (SCIPvarGetUbGlobal(vars[pos]) - SCIPgetSolVal(scip, sol, vars[pos]));
         }*/
         contribution = getPotentialContributed(scip, sol, vars[pos], vals[j], sgn);
         if( !SCIPisInfinity(scip, REALABS(contribution)) )
         {
            potential[i] += contribution;
         }
         else
         {
            inftycounter[i]++;
         }

         if( !SCIPisZero(scip, slack[i]) )
         {
            nviolatedrows[pos]++;
         }
      }


      /* create a new linear constraint, representing the old one */
      SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &subcons[i], SCIProwGetName(rows[i]),
               nnonz, consvars, vals, lhs, rhs) );

      /*Adds the Constraint and release it.*/
      SCIP_CALL( SCIPaddCons(subscip, subcons[i]) );
      SCIP_CALL( SCIPreleaseCons(subscip, &subcons[i]) );
      SCIPfreeBufferArray(subscip, &consvars);
   }

   /* author bzfhende: TODO don't you need an #ifdef here? */
   {
      FILE* probfile;
      char* bfilename;
      char probfilename[SCIP_MAXSTRLEN];

      bfilename = basename(heurdata->filename);

      sprintf(probfilename, "%s.lp", bfilename);

      SCIPdebugMessage("All temp vars initialized");

      probfile = fopen(probfilename, "w");
      /* test if file exists */
      if( NULL != probfile )
      {
         SCIP_CALL(SCIPprintOrigProblem(subscip, probfile, "lp", FALSE));
         fclose(probfile);
      }
      else
      {
         SCIPwarningMessage(scip, "Could not open file <%s> for storing infeasible repair subproblem\n", probfilename);
      }
   }

   /* get the greedy order */
   for( i = 0; i < nvars; ++i )
   {
      permutation[i] = i;
   }
   SCIPsortIntInt(nviolatedrows, permutation, nvars);

   /* loops over variables and greedily fix variables, but preserve the cover property that enough slack is given to violated rows */
   heurdata->nvarfixed = 0;
   for( i = 0; i < nvars; ++i )
   {
      /* author bzfhende: TODO get rid of the printf's */
      printf("Try to fix %s",SCIPvarGetName(vars[permutation[i]]));
      if( tryFixVar(scip, sol, potential, slack, vars[permutation[i]], inftycounter, heurdata) )
      {
         SCIP_Bool infeasible = FALSE;
         SCIP_Bool fixed = TRUE;
         SCIP_CALL( SCIPfixVar(subscip, subvars[permutation[i]], SCIPgetSolVal(scip, sol, vars[permutation[i]]), &infeasible, &fixed) );
         assert(!infeasible && fixed);
         heurdata->nvarfixed++;
         printf("  fixed!!!");
      }
      else
      {
         printf("  not.");
      }
      printf("\n");
    }
   printf("fixings finished\n\n");

   /*todo delete all redundant code, especially all addings of original variables and constraints !!!*/
   if( heurdata->useslackvars )
   {
      /*SCIP_VAR** tmpvars = NULL;*/
      SCIP_ROW** tmprows = NULL;
      int ntmprows;
      /* Gets all original variables */
      /*SCIP_CALL( SCIPgetVarsData(scip, &tmpvars, &ntmpvars, NULL, NULL, NULL, NULL) );*/
      tmprows = SCIPgetLPRows(scip);
      ntmprows = SCIPgetNLPRows(scip);

      /* Adds all original variables and decides which variables should be fixed */
      for( i = 0; i < nvars; ++i )
      {
         SCIP_CONS* cons;
         SCIP_Real lborig;
         SCIP_Real uborig;
         SCIP_Real lb;
         SCIP_Real ub;
         SCIP_Real objval;
         SCIP_Real value;
         SCIP_Real sla;
         SCIP_VARTYPE vartype;
         char varname[1024];
         char slackvarname[1024];
         char consvarname[1024];

         sla = 0.0;
         lborig = SCIPvarGetLbGlobal(subvars[i]);
         uborig = SCIPvarGetUbGlobal(subvars[i]);
         value = SCIPgetSolVal(scip, sol, subvars[i]);
         vartype = SCIPvarGetType(subvars[i]);

         nviolatedrows[i] = 0;

         sprintf(varname, "sub_%s", SCIPvarGetName(subvars[i]));

         /* if the value of x is lower than the lower bound, sets the slack to an correcting value */
         if( SCIPisFeasLT(scip, value, lborig) )
         {
            lb = value;
            sla = lborig - value;
            SCIPchgVarLbGlobal(subscip,subvars[i],lb);
         }
         else
         {
            lb = lborig;
         }

         /* if the value of x is bigger than the upper bound, sets the slack to an correcting value*/
         if( SCIPisFeasGT(scip, value, uborig) )
         {
            ub = value;
            sla = uborig - value;
            SCIPchgVarUbGlobal(subscip,subvars[i],ub);
         }
         else
         {
            ub = uborig;
         }

         objval = SCIPvarGetObj(vars[i])*factor;

         if(SCIPisZero(scip, objval))
         {
            objval = 0.0;
         }
         /* if an binary variable is out of bound, generalize it to an integer variable */
         if( !SCIPisFeasZero(scip, sla) && SCIP_VARTYPE_BINARY == vartype )
         {
            vartype = SCIP_VARTYPE_INTEGER;
            SCIPchgVarType(subscip,subvars[i],vartype, &success);
         }

         /* if necessary adds an constraint to represent the original bounds of x.*/
         if( !SCIPisFeasEQ(scip, sla, 0.0) )
         {
            SCIP_VAR* newvar;
            sprintf(slackvarname, "artificialslack_%s", SCIPvarGetName(vars[i]));
            sprintf(consvarname, "boundcons_%s", SCIPvarGetName(vars[i]));

            /* initialize and add an artificial slack variable */
            if(heurdata->useobjfactor)
            {
               SCIP_CALL( SCIPcreateVarBasic(subscip, &newvar, slackvarname, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY));
            }
            else
            {
               SCIP_CALL( SCIPcreateVarBasic(subscip, &newvar, slackvarname, 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS));
            }
            SCIP_CALL(SCIPaddVar(subscip, newvar));

            /* set the value of the slack variable to 1 to punish the use of it. */
            SCIP_CALL( SCIPsetSolVal(subscip, subsol, newvar, 1.0) );

            /* adds a linear constraint to represent the old bounds */
            SCIP_CALL( SCIPcreateConsBasicVarbound(subscip, &cons, consvarname, subvars[i], newvar, sla, lborig, uborig) );
            SCIP_CALL( SCIPaddCons(subscip, cons) );
            SCIP_CALL( SCIPreleaseVar(subscip, &newvar) );
            SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

            /* increases the counter for violated vars */
            heurdata->nviolatedvars++;
            }
         }
      /* todo iterate over slack*/
      /* Adds all original constraints and decides which variables should be fixed */
      for (i = 0; i < ntmprows; ++i)
      {
         SCIP_VAR* newvar;
         char varname[1024];

         /*if necessary adds an new artificial slack variable*/
         if( !SCIPisFeasEQ(subscip, slack[i], 0.0) )
         {
            sprintf(varname, "artificialslack_%s", SCIProwGetName(tmprows[i]));
            SCIP_CALL( SCIPcreateVarBasic(subscip, &newvar, varname, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
            SCIP_CALL( SCIPaddVar(subscip, newvar) );
            SCIP_CALL( SCIPsetSolVal(subscip, subsol, newvar, 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(subscip, subcons[i], newvar, slack[i]) );
            SCIP_CALL( SCIPreleaseVar(subscip, &newvar) );
         }
      }
   }

   /*Adds the given solution to the subscip. */
   heurdata->improovedoldsol = SCIPgetSolOrigObj(subscip, subsol);
   SCIP_CALL( SCIPaddSolFree(subscip, &subsol, &success) );

   if( !success )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Repair was not good enough.\n");
   }

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
   SCIP_CALL(SCIPsetObjlimit(subscip,1.0));
   /*SCIP_CALL(SCIPsetLongintParam(subscip, "limits/nodes", 1));*/

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
      for( i = 0; i < nvars; ++i )
      {
         SCIP_CALL( SCIPreleaseVar(subscip, &subvars[i]) );
      }
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL(SCIPfree(&subscip));
      return SCIP_OKAY;
   }

   success = FALSE;

   /* author bzfhende: TODO comment this */
   if( SCIPgetBestSol(subscip) != NULL )
   {
      heurdata->improovedoldsol = SCIPgetSolOrigObj(subscip, SCIPgetBestSol(subscip));
      /* print solving statistics of subproblem if we are in SCIP's debug mode */
      SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );

      assert(SCIPgetNSols(subscip) > 0);
      SCIP_CALL(createNewSol(heurdata, scip, subscip, subvars, heur, SCIPgetBestSol(subscip), &success));
      /*assert(success);*/

      *result = SCIP_FOUNDSOL;
   }
   else
   {
      printf("No solution found!\n");
   }

   heurdata->subiters = SCIPgetNLPIterations(subscip);
   heurdata->subnodes = SCIPgetNTotalNodes(subscip);
   heurdata->subpresoltime = SCIPgetPresolvingTime(subscip);
   heurdata->runs = SCIPgetNRuns(subscip);

   /* terminates the solving process  */
TERMINATE:
   SCIPfreeSol(scip, &sol);
   /*SCIPfreeSol(subscip, &subsol);*/
   SCIPfreeBufferArrayNull(scip, &nviolatedrows);
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(subscip, &subvars[i]) );
   }

   /* author bzfhende: TODO buffer arrays must be freed in reverse order of allocation */
   SCIPfreeBufferArrayNull(scip, &subvars);
   SCIPfreeBufferArrayNull(scip, &permutation);
   SCIPfreeBufferArrayNull(scip, &inftycounter);
   SCIPfreeBufferArrayNull(scip, &potential);
   SCIPfreeBufferArrayNull(scip, &slack);
   SCIPfreeBufferArrayNull(scip, &subcons);
   if( NULL != subscip )
   {
      SCIPfree(&subscip);
   }
   printf("repair finished\n");
   return SCIP_OKAY;
}


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
   heurdata->runs = 0;

   heurdata->nviolatedvars = 0;
   heurdata->norvars = 0;
   heurdata->relviolatedvars = 0;
   heurdata->nviolatedcons = 0;
   heurdata->norcons = 0;
   heurdata->relviolatedcons = 0;

   heurdata->nvarfixed = 0;
   heurdata->relvarfixed = -1;

   heurdata->orsolval = SCIP_INVALID;

   heurdata->improovedoldsol = SCIP_REAL_MAX;

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
   SCIP_Real relfixed;
   char solval[1024];
   char message[2048];
   int violateds;
   int ninvars;
   int ninvcons;
   int nvars;
   int ncons;
   int iterations;
   int nodes;
   int runs;

   heurdata = SCIPheurGetData(heur);
   violateds = heurdata->nviolatedvars+heurdata->nviolatedcons;
   ninvars = heurdata->nviolatedvars;
   ninvcons = heurdata->nviolatedcons;
   nvars = heurdata->norvars;
   ncons = heurdata->norcons;
   iterations = heurdata->subiters;
   nodes = heurdata->subnodes;
   time = heurdata->subpresoltime;
   runs = heurdata->runs;

   if( SCIP_INVALID == heurdata->orsolval )
   {
      sprintf(solval,"--");
   }
   else
   {
      sprintf(solval,"%15.9g",heurdata->orsolval);
   }

   heurdata->relviolatedvars = MAX((SCIP_Real)heurdata->norvars, 1.0);
   heurdata->relviolatedvars = heurdata->nviolatedvars/heurdata->relviolatedvars;
   heurdata->relviolatedcons = MAX((SCIP_Real)heurdata->norcons, 1.0);
   heurdata->relviolatedcons = heurdata->nviolatedcons/heurdata->relviolatedcons;

   heurdata->relvarfixed = MAX((SCIP_Real)heurdata->norvars, 1.0);
   heurdata->relvarfixed = heurdata->nvarfixed/heurdata->relvarfixed;
   relvars = heurdata->relviolatedvars;
   relcons = heurdata->relviolatedcons;
   relfixed = heurdata->relvarfixed;


   sprintf(message, "<repair> \n\t total violateds: %d\n\n\t violated variables: %d\n\t total variables: %d\n\t relative violated variables: %.2f%%\n", violateds, ninvars, nvars, 100 * relvars);
   sprintf(message, "%s  \n\n\t violated constraints: %d\n\t total constraints: %d\n\t relative violated constraints: %.2f%%\n", message, ninvcons, ncons, 100* relcons);
   sprintf(message, "%s  \n\n\t fixed variables: %d\n\t relative fixed varibales: %.2f%%\n", message, heurdata->nvarfixed, 100* relfixed);
   sprintf(message, "%s  \n\n\t iterations: %d\n\t nodes: %d\n\t number of runs: %d\n\t presolve time: %.2f s\n", message,iterations,nodes,runs,time);
   sprintf(message, "%s  \n\n\t Value of repairs best solution: %s\n improoved orsolval: %6f \n</repair>\n\n", message,solval, heurdata->improovedoldsol);
   /* prints all statistic data for a user*/
   SCIPstatistic(
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, message);
   )
   return SCIP_OKAY;
}

#ifdef REPAIRWRITEPROB
static
SCIP_RETCODE writeDebugInformation(
   SCIP*                 scip,
   SCIP*                 subscip,
   SCIP_SOL*             subsol,
   SCIP_HEURDATA*        heurdata
   )
{
   SCIPdebugMessage("Print files:\n");
   {
      FILE* solfile;
      FILE* probfile;
      char* bfilename;
      char solfilename[SCIP_MAXSTRLEN];
      char probfilename[SCIP_MAXSTRLEN];

      bfilename = basename(heurdata->filename);

      sprintf(solfilename, "%s.sol", bfilename);
      sprintf(probfilename, "%s.cip", bfilename);

      SCIPdebugMessage("All temp vars initialized");

      solfile = fopen(solfilename, "w");



      /* test if file exists */
      if( NULL != solfile )
      {
         SCIP_CALL(SCIPprintSol(scip, subsol, solfile, TRUE));
         fclose(solfile);
       }
       else
       {
          SCIPwarningMessage(scip, "Could not open file <%s> for storing infeasible repair solution\n", solfilename);
       }

      probfile = fopen(probfilename, "w");
      /* test if file exists */
      if( NULL != probfile )
      {
         SCIP_CALL(SCIPprintOrigProblem(subscip, probfile, "cip", FALSE));
         fclose(probfile);
      }
      else
      {
         SCIPwarningMessage(scip, "Could not open file <%s> for storing infeasible repair subproblem\n", probfilename);
      }
   }
}
#else
#define writeDebugInformation(scip, subscip, subsol, heurdata) SCIP_OKAY
#endif



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
   SCIP_HEURDATA* heurdata = NULL;
   SCIP_RETCODE retcode;
   SCIP_Bool success;

   retcode = SCIP_OKAY;
   heurdata = SCIPheurGetData(heur);
   printf("%s\n",heurdata->filename);

   /* if repair allready run, stop*/
   if(0 < SCIPheurGetNCalls(heur)){
      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }

   /* checks the result pointer*/
   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( !SCIPhasCurrentNodeLP(scip) )
      return SCIP_OKAY;

   if( !SCIPisLPConstructed(scip) )
   {
      SCIP_CALL(SCIPconstructLP(scip, &success));

      if( success )
         return SCIP_OKAY;
   }

   /* create zero solution */
   SCIP_CALL( SCIPcreateOrigSol(scip, &(heurdata->infsol), heur) );

   heurdata = SCIPheurGetData(heur);

   /* use read method to enter solution from a file */
   retcode = readSol(scip, heurdata->infsol, heurdata->filename);

   if( SCIP_NOFILE == retcode )
   {
      if( strcmp(heurdata->filename, DEFAULT_FILENAME) != 0 )
         SCIPwarningMessage(scip, "cannot open file <%s> for reading\n",
               heurdata->filename);

      SCIPfreeSol(scip, &(heurdata->infsol));
      return SCIP_OKAY;
   }
   else if( retcode != SCIP_OKAY )
   {
      SCIPfreeSol(scip, &(heurdata->infsol));
      return SCIP_OKAY;
   }
   SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "Repair: Solution file read.\n");

   /* checks the integrality of all discrete variable */
   SCIP_CALL( checkCands(scip, heurdata->infsol, heurdata->roundit, &success) );
   if( !success )
   {
      SCIPdebugMessage("Hello Termination\n");
      SCIPfreeSol(scip, &(heurdata->infsol));
      return SCIP_OKAY;
   }
   *result = SCIP_DIDNOTFIND;


   if(heurdata->usevarfix)
   {
      retcode = varFixings(scip, heur, result);
   }
   else
   {
      retcode = extendedProgramm(scip, heur, result);
   }
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
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/roundit",
         "True : fractional variables which are not fractional in the given solution are rounded, "
         "FALSE : solving process of this heuristic is stopped. ",
         &heurdata->roundit, FALSE, DEFAULT_ROUNDIT, NULL, NULL));

   /* add bool parameter for decision how the objective function should be */
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/useobjfactor",
         "should a scaled objective function for original variables be used in repair subproblem?",
         &heurdata->useobjfactor, FALSE, DEFAULT_USEOBJFACTOR, NULL, NULL));

   /* add bool parameter for decision if variable fixings should be used */
      SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/usevarfix",
            "should variable fixings be used in repair subproblem?",
            &heurdata->usevarfix, FALSE, DEFAULT_USEVARFIX, NULL, NULL));
   /* add bool parameter for decision how the objective function should be */
      SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/useslackvars",
            "should slack variables be used in repair subproblem?",
            &heurdata->useslackvars, FALSE, DEFAULT_USESLACKVARS, NULL, NULL));

      SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/alpha", "factor for the potential of var fixings",
               &heurdata->alpha, TRUE, DEFAULT_ALPHA, 0.0, 100.00, NULL, NULL) );


   return SCIP_OKAY;
}

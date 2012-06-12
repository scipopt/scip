/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   matscip.c
 * @brief  mex function for the SCIP-MATLAB interface
 * @author Timo Berthold
 * @author Ambros Gleixner
 */

#include "mex.h"

#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <scip/cons_linear.h>

/* standard "main" method for mex interface */
void mexFunction(
   int                   nlhs,               /* number of expected outputs */
   mxArray*              plhs[],             /* array of pointers to output arguments */
   int                   nrhs,               /* number of inputs */
   const mxArray*        prhs[]              /* array of pointers to input arguments */
   )
{
   SCIP* scip;
   SCIP_VAR** vars;
   SCIP_Real* objs;
   SCIP_Real* lhss;
   SCIP_Real* rhss;
   SCIP_Real* lbs;
   SCIP_Real* ubs;
   SCIP_Real* matrix;
   SCIP_Real* bestsol;
   SCIP_Real* objval;
   char* vartypes;
   char objsense[SCIP_MAXSTRLEN];

   int nvars;
   int nconss;
   int stringsize;
   int i;

   if( SCIPmajorVersion() < 2 )
   {
      mexErrMsgTxt("SCIP versions less than 2.0 are not supported\n");
      return;
   }

   /* initialize SCIP */
   SCIP_CALL_ABORT( SCIPcreate(&scip) );

   /* output SCIP information */
   SCIPprintVersion(scip, NULL);

   /* include default SCIP plugins */
   SCIP_CALL_ABORT( SCIPincludeDefaultPlugins(scip) );

   if( nlhs != 2 || nrhs != 8 )
      mexErrMsgTxt("invalid number of parameters. Call as [bestsol, objval] = matscip(matrix, lhs, rhs, obj, lb, ub, vartype, objsense)\n");

   if( mxIsSparse(prhs[0]) )
      mexErrMsgTxt("sparse matrices are not supported yet"); /* ???????? of course this has to change */

   /* get linear constraint coefficient matrix */
   matrix = mxGetPr(prhs[0]);
   if( matrix == NULL )
      mexErrMsgTxt("matrix must not be NULL");
   if( mxGetNumberOfDimensions(prhs[0]) != 2 )
      mexErrMsgTxt("matrix must have exactly two dimensions");

   /* get dimensions of matrix */
   nconss = mxGetM(prhs[0]);
   nvars = mxGetN(prhs[0]);
   assert(nconss > 0);
   assert(nvars > 0);

   /* get left hand sides of linear constraints */
   lhss = mxGetPr(prhs[1]);
   if( mxGetM(prhs[1]) != nconss )
      mexErrMsgTxt("dimension of left hand side vector does not match matrix dimension");
   assert(lhss != NULL);

   /* get right hand sides of linear constraints */
   rhss = mxGetPr(prhs[2]);
   if( mxGetM(prhs[2]) != nconss )
      mexErrMsgTxt("dimension of right hand side vector does not match matrix dimension");
   assert(rhss != NULL);

   /* get objective coefficients */
   objs = mxGetPr(prhs[3]);
   if( mxGetM(prhs[3]) != nvars )
      mexErrMsgTxt("dimension of objective coefficient vector does not match matrix dimension");

   /* get lower bounds of variables */
   lbs = mxGetPr(prhs[4]);
   if( mxGetM(prhs[4]) != nvars )
      mexErrMsgTxt("dimension of lower bound vector does not match matrix dimension");

   /* get upper bounds of variables */
   ubs = mxGetPr(prhs[5]);
   if( mxGetM(prhs[5]) != nvars )
      mexErrMsgTxt("dimension of upper bound vector does not match matrix dimension");

   /* allocate memory for variable type characters */
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &vartypes, nvars+1) );

   /* get variable types */
   if( mxGetString(prhs[6], vartypes, nvars+1)  != 0 )
      mexErrMsgTxt("Error when parsing variable types, maybe a wrong vector dimension?");

   /* get objective sense */
   stringsize = mxGetNumberOfElements(prhs[7]);
   if( stringsize != 3 )
      mexErrMsgTxt("objective sense must be a three character word: \"max\" or \"min\"");
   if( mxGetString(prhs[7], objsense, stringsize+1) != 0)
      mexErrMsgTxt("Error when parsing objective sense string");
   if( strcmp(objsense,"max") != 0 && strcmp(objsense,"min") != 0 )
      mexErrMsgTxt("objective sense must be either \"max\" or \"min\"");

   /* get output parameters */
   plhs[0] = mxCreateDoubleMatrix(nvars, 1, mxREAL);
   bestsol = mxGetPr(plhs[0]);
   plhs[1] = mxCreateDoubleScalar(mxREAL);
   objval  = mxGetPr(plhs[1]);

   /* create SCIP problem */
   SCIP_CALL_ABORT( SCIPcreateProb(scip, "mex_prob", NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* allocate memory for variable array */
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &vars, nvars) );

   /* create variables */
   for( i = 0; i < nvars; ++i)
   {
      SCIP_VARTYPE vartype;
      char varname[SCIP_MAXSTRLEN];

      /* convert vartype character to SCIP vartype */
      if( vartypes[i] == 'i' )
         vartype = SCIP_VARTYPE_INTEGER;
      else if( vartypes[i] == 'b' )
         vartype = SCIP_VARTYPE_BINARY;
      else if( vartypes[i] == 'c' )
         vartype = SCIP_VARTYPE_CONTINUOUS;
      else
         mexErrMsgTxt("unkown variable type");

      /* variables get canonic names x_i */
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "x_%d", i);

      /* create variable object and add it to SCIP */
      SCIP_CALL_ABORT( SCIPcreateVar(scip, &vars[i], varname, lbs[i], ubs[i], objs[i],
            vartype, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      assert(vars[i] != NULL);
      SCIP_CALL_ABORT( SCIPaddVar(scip, vars[i]) );
   }

   /* create linear constraints */
   for( i = 0; i < nconss; ++i )
   {
      SCIP_CONS* cons;
      char consname[SCIP_MAXSTRLEN];
      int j;

      /* constraints get canonic names cons_i */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "cons_%d", i);

      /* create empty linear constraint */
      SCIP_CALL_ABORT( SCIPcreateConsLinear(scip, &cons, consname, 0, NULL, NULL, lhss[i], rhss[i],
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      /* add non-zero coefficients to linear constraint */
      for( j = 0; j < nvars; ++j )
      {
         if( !SCIPisFeasZero(scip, matrix[i+j*nconss]) )
         {
            SCIP_CALL_ABORT( SCIPaddCoefLinear(scip, cons, vars[j], matrix[i+j*nconss]) );
         }
      }

      /* add constraint to SCIP and release it */
      SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );
      SCIP_CALL_ABORT( SCIPreleaseCons(scip, &cons) );
   }

   /* set objective sense in SCIP */
   if( strcmp(objsense,"max") == 0)
   {
      SCIP_CALL_ABORT( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );
   }
   else if( strcmp(objsense,"min") == 0)
   {
      SCIP_CALL_ABORT( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );
   }
   else
      /* this should have been caught earlier when parsing objsense */
      mexErrMsgTxt("unkown objective sense");

   /* solve SCIP problem */
   SCIP_CALL_ABORT( SCIPsolve(scip) );

   /* if SCIP found a solution, pass it back into MATLAB output parameters */
   if( SCIPgetNSols > 0 )
   {
      SCIP_SOL* scipbestsol;

      /* get incumbent solution vector */
      scipbestsol = SCIPgetBestSol(scip);
      assert(scipbestsol != NULL);

      /* get objective value of incumbent solution */
      *objval = SCIPgetSolOrigObj(scip, scipbestsol);
      assert(!SCIPisInfinity(scip, REALABS(*objval)));

      /* copy solution values into output vector */
      for( i = 0; i < nvars; ++i )
         bestsol[i] = SCIPgetSolVal(scip,scipbestsol,vars[i]);
   }

   /* release variables */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL_ABORT( SCIPreleaseVar(scip, &vars[i]) );
   }

   /* free memory for variable arrays */
   SCIPfreeMemoryArray(scip, &vartypes);
   SCIPfreeMemoryArray(scip, &vars);

   /* deinitialize SCIP */
   SCIP_CALL_ABORT( SCIPfree(&scip) );

   /* check for memory leaks */
   BMScheckEmptyMemory();

   return;
}

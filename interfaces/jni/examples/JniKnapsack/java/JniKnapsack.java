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

/**@file   Knapsack.java
 * @brief  main class for JNI knapsack example
 * @author Timo Berthold
 * @author Stefan Heinz
 */

import de.zib.jscip.nativ.*;
import de.zib.jscip.nativ.jni.*;


public class JniKnapsack {

   /** main function */
   public static void main(
      String[]         args                /**< command line arguments */
      )
   {
      /** load SCIP libraries */
      JniScipLibraryLoader.loadLibrary();

      try
      {

	 /**@note C pointer are longs in the JNI interface */
	 long scip;
	 long consLinear;
	 long consKnapsack;
	 long consLogicor;
	 long consSetpart;
	 long consSetpack;
	 long consSetcover;

	 /* create the SCIP environment */
	 JniScip env = new JniScip();

	 /* create the SCIP variable environment */
	 JniScipVar envVar = new JniScipVar();

	 /* create the SCIP  constraint environment */
	 JniScipCons envCons = new JniScipCons();

	 /* create the SCIP knapsack constraint environment */
	 JniScipConsKnapsack envConsKnapsack = new JniScipConsKnapsack();

	 /* create the SCIP linear constraint environment */
	 JniScipConsLinear envConsLinear = new JniScipConsLinear();

	 /* create SCIP instance */
	 scip = env.create();

	 /* print version to stdout */
	 env.printVersion(scip, 0);

	 /* set message handler of SCIP quiet (as a result SCIP does not produce any output to stdout and stderr) */
	 env.setMessagehdlrQuiet(scip, true);

	 /* write all SCIP output to log file */
	 env.setMessagehdlrLogfile(scip, "scip.log");

	 /* include default plugins od f SCIP */
	 env.includeDefaultPlugins(scip);

	 /* create empty problem with name "knapsack" */
	 env.createProbBasic(scip, "knapsack");

	 /* create a knapsack constraint with zero variable and a capacity og 10 (variables are added later) */
	 consKnapsack = envConsKnapsack.createConsBasicKnapsack(scip, "knapsack", 0, null, null, 10);

	 /* create linear constraint with zero variable, left hand side of zero, right hand side of 100 (variables are
	  * added later)
	  */
	 consLinear = envConsLinear.createConsBasicLinear(scip, "linear", 0, null, null, 0.0, 100.0);


         /* create 10 variables*/
         int nvars = 10;
         long vars[] = new long[nvars];

         /* instantiating and setting the variables to be binary */
         for( int i = 0; i < nvars; i++ )
         {
            /* create variable */
            vars[i] = env.createVarBasic(scip, "x" + i, 0.0, 1.0, -1.0, 0);
            System.out.println("created variables <" + envVar.varGetName(vars[i]) + ">");

            /* add variable to problem */
            env.addVar(scip, vars[i]);

            /* add variable to constraints */
            envConsKnapsack.addCoefKnapsack(scip, consKnapsack, vars[i], i + 10);
            envConsLinear.addCoefLinear(scip, consLinear, vars[i], i + 10);

	    System.out.println("variables pointer <" + vars[i] + ">");
	 }

         /* add constraints to problem */
	 env.addCons(scip, consKnapsack);
	 env.addCons(scip, consLinear);


	 int nvals =  envConsLinear.getNVarsLinear(scip, consLinear);
	 double[] vals = envConsLinear.getValsLinear(scip, consLinear);

	 for( int i = 0; i < nvals; ++i )
	 {
	    System.out.println("*** " + vals[i]);
	 }

         /* release constraints */
	 env.releaseCons(scip, consKnapsack);
	 env.releaseCons(scip, consLinear);

	 /* solve problem */
	 env.solve(scip);

	 /* get best solution */
	 long sol = env.getBestSol(scip);

	 int status = env.getStatus(scip);

	 System.out.println("*** " + status + "\n*** Value found: " + env.getPrimalbound(scip));

	 /* in case of an optimal solution was found display solution */
	 if( status == JniScipStatus.SCIP_STATUS_OPTIMAL )
	 {
            /* collect solution values */
	    double[] values = env.getSolVals(scip, sol, nvars, vars);

            /* print solution value */
	    for (int i = 0; i < nvars; ++i )
	    {
	       System.out.println(envVar.varGetName(vars[i]) + " : " + values[i]);
	    }
	 }
      }
      catch (NativeScipException e)
      {
	 System.out.println(e.getMessage());
      }

      try
      {
	 /**@note C pointer are longs in the JNI interface */
	 long scip;

	 /* create the SCIP environment */
	 JniScip env = new JniScip();

	 /* create the SCIP variable environment */
	 JniScipVar envVar = new JniScipVar();

	 scip = env.create();

	 /* include default plugins of SCIP */
	 env.includeDefaultPlugins(scip);

	 env.readProb(scip, "data/test.lp", "");

	 env.readSol(scip, "data/solution.sol");

	 long[] sols = env.getSols(scip);

	 long[] vars = env.getVars(scip);
	 int nvars = env.getNVars(scip);

	 for( int i = 0; i < nvars; ++i )
	 {
	    System.out.println("variable <" + envVar.varGetName(vars[i]) + "> -> <" + env.getSolVal(scip, sols[0], vars[i]) + ">");
	 }

	 env.presolve(scip);

	 env.printBestSol(scip, 0, false);

	 env.free(scip);
      }
      catch (NativeScipException e)
      {
	 System.out.println(e.getMessage());
      }

   }
}

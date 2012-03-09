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

/**@file   cppmain.cpp
 * @brief  main file for p-median pricer example
 * @author Joerg Rambau
 * @author Andreas Tuchscherer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

// standard library includes
#include <stdio.h>
#include <iostream>
#include <vector>

// scip includes
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

// user defined includes
#include "pricer_distance.h"


// namespace usage 
using namespace std;
using namespace scip;


static
SCIP_RETCODE readParams(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< parameter file name, or NULL */
   )
{
   if( filename != NULL )
   {
      if( SCIPfileExists(filename) )
      {
         printf("reading parameter file <%s>\n", filename);
         SCIP_CALL( SCIPreadParams(scip, filename) );
      }
      else
         printf("parameter file <%s> not found - using default parameters\n", filename);
   }
   else if( SCIPfileExists("pmedian.set") )
   {
      printf("reading parameter file <pmedian.set>\n");
      SCIP_CALL( SCIPreadParams(scip, "pmedian.set") );
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE fromCommandLine(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< input file name */
   )
{
   /********************
    * Problem Creation *
    ********************/

   printf("\nread problem <%s>\n", filename);
   printf("============\n\n");
   SCIP_CALL( SCIPreadProb(scip, filename, NULL) );


   /*******************
    * Problem Solving *
    *******************/

   /* solve problem */
   printf("\nsolve problem\n");
   printf("=============\n\n");
   SCIP_CALL( SCIPsolve(scip) );

   printf("\nprimal solution:\n");
   printf("================\n\n");
   SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );


   /**************
    * Statistics *
    **************/

   printf("\nStatistics\n");
   printf("==========\n\n");

   SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   return SCIP_OKAY;
}

// static
// RETCODE interactive(
//    SCIP*                 scip                /**< SCIP data structure */
//    )
// {
//    /* start user interactive mode */
//    SCIP_CALL( SCIPstartInteraction(scip) );

//    return SCIP_OKAY;
// }

static
SCIP_RETCODE runSCIP(
   int                   argc,
   char**                argv
   )
{
   SCIP* scip = NULL;


   /***********************
    * Version information *
    ***********************/

   SCIPprintVersion(NULL);
   printf("\n");


   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include distance pricer */
   ObjPricerDistance* distance_pricer_ptr = new ObjPricerDistance(scip);
   SCIP_CALL(SCIPincludeObjPricer(scip, distance_pricer_ptr, true));
   
   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );


   /**************
    * Parameters *
    **************/
   
   if( argc >= 3 )
   {
      SCIP_CALL( readParams(scip, argv[2]) );
   }
   else
   {
      SCIP_CALL( readParams(scip, NULL) );
   }


   /**************
    * Start SCIP *
    **************/

   if( argc >= 2 )
   {
      SCIP_CALL( fromCommandLine(scip, argv[1]) );
   }
   else
   {
      printf("\n");

      // SCIP_CALL( interactive(scip) );

      /*** read and solve problem ***/

      const int num_points (3);
      const int num_centers(2);

      SCIP_CALL(SCIPcreateProb(scip, "p-median", 0, 0, 0, 0, 0, 0, 0));
      
      SCIP_CONS* cons;
      
      vector<SCIP_CONS*> setpart_consptr_vector(num_points);
      for (int i = 0; i < num_points; ++i)
      {
         SCIP_CALL( SCIPcreateConsSetpart(scip, &cons, "set partitioning", 
               0, 0, 
               true,  /* initial */         
               false, /* separate */        
               true,  /* enforce */         
               true,  /* check */           
               true,  /* propagate */       
               false, /* local */           
               true,  /* modifiable */      
               false, /* dynamic */         
               false, /* removable */       
               false  /* stickingatnode */ ) );
         SCIP_CALL(SCIPaddCons(scip, cons));
         setpart_consptr_vector[i] = cons;
      }

      vector<SCIP_CONS*> setpack_consptr_vector(num_points);
      for (int i = 0; i < num_points; ++i)
      {
         SCIP_CALL( SCIPcreateConsSetpack(scip, &cons, "set packing", 
               0, 0, 
               true,  /* initial */         
               false, /* separate */        
               true,  /* enforce */         
               true,  /* check */           
               true,  /* propagate */       
               false, /* local */           
               true,  /* modifiable */      
               false, /* dynamic */         
               false, /* removable */       
               false  /* stickingatnode */ ) );
         SCIP_CALL(SCIPaddCons(scip, cons));
         setpack_consptr_vector[i] = cons;
      }

      SCIP_CONS* card_cons;
      SCIP_CALL( SCIPcreateConsLinear(scip, &card_cons, "cardinality", 
            0, 0, 0,
            -SCIPinfinity(scip), num_centers,
            true,  /* initial */         
            false, /* separate */        
            true,  /* enforce */         
            true,  /* check */           
            true,  /* propagate */       
            false, /* local */           
            true,  /* modifiable */      
            false, /* dynamic */         
            false, /* removable */       
            false  /* stickingatnode */ ) );
      SCIP_CALL(SCIPaddCons(scip, card_cons));

      distance_pricer_ptr->init(num_points,
				num_centers,
				setpart_consptr_vector,
				setpack_consptr_vector,
				card_cons);

      //SCIP_CALL(distance_pricer_ptr->scip_farkas(scip, 0));

      SCIP_CALL(SCIPactivatePricer(scip, SCIPfindPricer(scip, "Distance_Pricer")));
   }

   SCIP_CALL(SCIPprintOrigProblem(scip, NULL, NULL, FALSE));
   

   /***
       Solve:
   ***/
   SCIP_CALL(SCIPsolve(scip));

   /***
       Statistics:
    ***/
   SCIP_CALL(SCIPprintStatistics(scip, NULL));

   SCIP_CALL(SCIPprintTransProblem(scip, NULL, NULL, FALSE));

   SCIP_CALL(SCIPprintBestSol(scip, NULL, FALSE));

   
   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

int
main(
   int                   argc,
   char**                argv
   )
{
   SCIP_RETCODE retcode;

   retcode = runSCIP(argc, argv);
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode, stderr);
      return -1;
   }

   return 0;
}

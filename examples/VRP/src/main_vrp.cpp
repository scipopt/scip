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

/**@file
 * @brief  main file for VRP pricer example
 * @author Andreas Bley
 * @author Marc Pfetsch
 *
 * We want to solve the vehicle routing problem on a graph \f$G = (V,E)\f$ with
 * \f$V = J \cup \{d\}\f$, where d is the depot and the distances are given by the
 * length function \f$l_e: E \rightarrow R_{\geq 0}\f$.
 *
 * Consider the MIP formulation
 *
 * \f[
 *  \begin{array}[t]{rll}
 *    \min &  \displaystyle \sum_{e \in E} l_e y_e \\
 *         & & \\
 *   s.t.  & -y_e + \sum_{t \in T_k} a^t_e x_t  \leq 0, &  \forall e \in E\\
 *         &  \displaystyle \sum_{t \in T_k} a^t_j x_t = 1, &  \forall j \in J \\
 *         &  y(\delta(j)) = 2, &  \forall j \in J \\
 *         &  y_e \in \{0,1,2\},  & \forall e \in E \\
 *         &  x_t  \in [0,1], & \forall t \in T_k
 *  \end{array}
 * \f]
 *
 * where \f$T_k\f$ is the set of tours visiting at most k customers
 * with repetitions of customers allowed and \f$a^t_e\f$ (\f$a^t_j\f$) counts how often
 * edge e (node j) is traversed in \f$t \in T_k\f$.
 *
 * Examples and the file format are given at https://neo.lcc.uma.es/vrp/vrp-instances/capacitated-vrp-instances/.
 */

/* standard library includes */
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

/* scip includes */
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

/* user defined includes */
#include "pricer_vrp.h"


/* namespace usage */
using namespace std;
using namespace scip;


/** read VRP problem */
static
int read_problem(
   const char*           filename,           /**< filename */
   int&                  num_nodes,          /**< number of nodes in instance */
   int&                  capacity,           /**< capacity in instance */
   vector<int>&          demand,             /**< array of demands of instance */
   vector<vector<int> >& dist                /**< distances between nodes */
   )
{
   const string DIMENSION           = "DIMENSION";
   const string DEMAND_SECTION      = "DEMAND_SECTION";
   const string DEPOT_SECTION       = "DEPOT_SECTION";
   const string EDGE_WEIGHT_TYPE    = "EDGE_WEIGHT_TYPE";
   const string EUC_2D              = "EUC_2D";
   const string EXPLICIT            = "EXPLICIT";
   const string LOWER_DIAG_ROW      = "LOWER_DIAG_ROW";
   const string EDGE_WEIGHT_FORMAT  = "EDGE_WEIGHT_FORMAT";
   const string EDGE_WEIGHT_SECTION = "EDGE_WEIGHT_SECTION";
   const string NODE_COORD_SECTION  = "NODE_COORD_SECTION";
   const string CAPACITY            = "CAPACITY";

   ifstream file(filename);

   num_nodes = -1;
   capacity = -1;

   if ( ! file )
   {
      cerr << "Cannot open file " << filename << endl;
      return 1;
   }

   string edge_weight_type   = "";
   string edge_weight_format = "";
   vector<int> x;
   vector<int> y;

   while ( file )
   {
      //--------------------
      // Read keyword.
      //--------------------
      string key;
      string dummy;
      file >> key;

      if ( key == DIMENSION )
      {
         file >> dummy;
         file >> num_nodes;

         assert( num_nodes >= 0 );
         demand.resize(num_nodes, 0); /*lint !e732 !e747*/
         dist.resize(num_nodes); /*lint !e732 !e747*/
         for (int i = 0; i < num_nodes; ++i)
            dist[i].resize(i, 0);  /*lint !e732 !e747*/
      }

      if ( key == CAPACITY )
      {
         file >> dummy;
         file >> capacity;
      }
      else if ( key == EDGE_WEIGHT_TYPE )
      {
         file >> dummy;
         file >> edge_weight_type;
         if ( edge_weight_type != EUC_2D && edge_weight_type != EXPLICIT )
         {
            cerr << "Wrong " << EDGE_WEIGHT_TYPE << " " << edge_weight_type << endl;
            return 1;
         }
         if ( edge_weight_type == EUC_2D )
         {
            assert( num_nodes >= 0 );
            x.resize(num_nodes, 0); /*lint !e732 !e747*/
            y.resize(num_nodes, 0); /*lint !e732 !e747*/
         }
      }
      else if ( key == EDGE_WEIGHT_FORMAT )
      {
         file >> dummy;
         file >> edge_weight_format;
      }
      else if ( key == EDGE_WEIGHT_FORMAT + ":" )
      {
         file >> edge_weight_format;
      }
      else if ( key == EDGE_WEIGHT_SECTION )
      {
         if ( edge_weight_type   != EXPLICIT || edge_weight_format != LOWER_DIAG_ROW )
         {
            cerr << "Error. Unsupported edge length type." << endl;
            return 1;
         }
         assert( num_nodes >= 0 );
         for (int i = 0; i < num_nodes; ++i)
         {
            for (int j = 0; j < i; ++j)
            {
               int l;
               file >> l;
               dist[i][j] = l;  /*lint !e732 !e747*/
            }
         }
      }
      else if ( key == NODE_COORD_SECTION )
      {
         if ( edge_weight_type != EUC_2D )
         {
            cerr << "Error. Data file contains " << EDGE_WEIGHT_TYPE << " " << edge_weight_type << " and " << NODE_COORD_SECTION << endl;
            return 1;
         }
         assert( num_nodes >= 0 );
         for (int i = 0; i < num_nodes; ++i)
         {
            int j, xi, yi;
            file >> j;
            file >> xi;
            file >> yi;
            if ( j != i+1 )
            {
               cerr << "Error reading " << NODE_COORD_SECTION << endl;
               return 1;
            }
            x[i] = xi; /*lint !e732 !e747*/
            y[i] = yi; /*lint !e732 !e747*/
         }
         for (int i = 0; i < num_nodes; ++i)
         {
            for (int j = 0; j < i; ++j)
            {
               int dx = x[i] - x[j]; /*lint !e732 !e747 !e864*/
               int dy = y[i] - y[j]; /*lint !e732 !e747 !e864*/
               dist[i][j] = int( sqrt((double)dx*dx + dy*dy) + 0.5 ); /*lint !e732 !e747 !e790*/
            }
         }
      }
      else if ( key == DEMAND_SECTION )
      {
         assert( num_nodes >= 0 );
         for (int i = 0; i < num_nodes; ++i)
         {
            int j, d;
            file >> j;
            file >> d;
            if ( j != i+1 )
            {
               cerr << "Error reading " << DEMAND_SECTION << endl;
               return 1;
            }
            demand[i] = d; /*lint !e732 !e747*/
         }
      }
      else if ( key == DEPOT_SECTION )
      {
         for (int i = 0; i != -1 ;)
         {
            file >> i;
            if ( i != -1 && i != 1 )
            {
               cerr << "Error: This file specifies other depots than 1." << endl;
               return 1;
            }
         }
      }
      else 
      {
         (void) getline(file, dummy);
      }
   }

   return 0;
}


//------------------------------------------------------------
static
SCIP_RETCODE execmain(int argc, char** argv)
{
   SCIP* scip = NULL;

   cout << "Solving the vehicle routing problem using SCIP." << endl;
   cout << "Implemented by Andreas Bley." << endl << endl;

   if ( argc != 2 && argc != 3 )
   {
      cerr << "Usage: vrp [-h] datafile" << endl;
      cerr << "Options:" << endl;
      cerr << " -h  Uses hop limit instead of capacity limit for tours."<< endl;
      return SCIP_INVALIDDATA;
   }


   /**********************
    * Setup problem data *
    **********************/

   const char* VRP_PRICER_NAME = "VRP_Pricer";

   vector<vector<int> > dist;
   vector<int> demand;
   int capacity;
   int num_nodes;

   if ( read_problem(argv[argc-1], num_nodes, capacity, demand, dist) )
   {
      cerr << "Error reading data file " << argv[argc-1] << endl;
      return SCIP_READERROR;
   }
   assert( num_nodes >= 0 );
   assert( capacity >= 0 );

   cout << "Number of nodes: " << num_nodes << endl;

   if ( argc == 3 )
   {
      if ( string("-h") != argv[1] )
      {
         cerr << "Unknow option " << argv[2] << endl;
         return SCIP_PARAMETERUNKNOWN;
      }

      int total_demand = 0;
      for (int i = 1; i< num_nodes; ++i)
         total_demand += demand[i]; /*lint !e732 !e747*/

      if( total_demand == 0.0 )
      {
         cerr << "Total demand is zero!" << endl;
         return SCIP_INVALIDDATA;
      }

      capacity = (num_nodes - 1) * capacity / total_demand;
      demand.assign(num_nodes, 1);
      demand[0] = 0; /*lint !e747*/
      cout << "Max customers per tour: " << capacity << endl << endl;
   }
   else
      cout << "Max demand per tour: " << capacity << endl << endl;

   /**************
    * Setup SCIP *
    **************/

   /* initialize SCIP environment */
   SCIP_CALL( SCIPcreate(&scip) );

   /***********************
    * Version information *
    ***********************/

   SCIPprintVersion(scip, NULL);
   SCIPinfoMessage(scip, NULL, "\n");

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* set verbosity parameter */
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 5) );
   /* SCIP_CALL( SCIPsetBoolParam(scip, "display/lpinfo", TRUE) ); */

   /* create empty problem */
   SCIP_CALL( SCIPcreateProb(scip, "VRP", 0, 0, 0, 0, 0, 0, 0) );

   /* add arc-routing variables */
   char var_name[255];
   vector< vector<SCIP_VAR*> > arc_var( num_nodes ); /*lint !e732 !e747*/
   for (size_t i = 0; i < (size_t)num_nodes; ++i)
   {
      arc_var[i].resize(i, (SCIP_VAR*) NULL); /*lint !e732 !e747*/
      for (size_t j = 0; j < i; ++j)
      {
         SCIP_VAR* var;
         (void) SCIPsnprintf(var_name, 255, "E%d_%d", i, j );

         SCIP_CALL( SCIPcreateVar(scip,
                     &var,                   // returns new index
                     var_name,               // name
                     0.0,                    // lower bound
                     2.0,                    // upper bound
                     dist[i][j],             // objective
                     SCIP_VARTYPE_INTEGER,   // variable type
                     TRUE,                   // initial
                     FALSE,                  // forget the rest ...
                     NULL, NULL, NULL, NULL, NULL) );  /*lint !e732 !e747*/
         SCIP_CALL( SCIPaddVar(scip, var) );
         arc_var[i][j] = var; /*lint !e732 !e747*/
      }
   }

   /* add arc-routing - tour constraints */
   char con_name[255];
   vector< vector<SCIP_CONS*> > arc_con( num_nodes );  /*lint !e732 !e747*/
   for (size_t i = 0; i < (size_t)num_nodes; ++i)
   {
      arc_con[i].resize(i, (SCIP_CONS*)NULL); /*lint !e732 !e747*/
      for (size_t j = 0; j < i; ++j)
      {
         SCIP_CONS* con;
         (void) SCIPsnprintf(con_name, 255, "A%d_%d", i, j);
         SCIP_VAR* idx = arc_var[i][j]; /*lint !e732 !e747*/
         SCIP_Real coeff = -1;
         SCIP_CALL( SCIPcreateConsLinear(scip, &con, con_name, 1, &idx, &coeff,
                     -SCIPinfinity(scip),    /* lhs */
                     0.0,                    /* rhs */
                     TRUE,                   /* initial */
                     FALSE,                  /* separate */
                     TRUE,                   /* enforce */
                     TRUE,                   /* check */
                     TRUE,                   /* propagate */
                     FALSE,                  /* local */
                     TRUE,                   /* modifiable */
                     FALSE,                  /* dynamic */
                     FALSE,                  /* removable */
                     FALSE) );               /* stickingatnode */
         SCIP_CALL( SCIPaddCons(scip, con) );
         arc_con[i][j] = con;  /*lint !e732 !e747*/
      }
   }

   /* add arc-routing - degree constraints */
   for (size_t i = 1; i < (size_t)num_nodes; ++i)
   {
      SCIP_CONS* con;
      (void) SCIPsnprintf(con_name, 255, "D%d", i);
      SCIP_CALL( SCIPcreateConsLinear(scip, &con, con_name, 0, 0, 0,
                  2.0,                    /* lhs */
                  2.0,                    /* rhs */
                  TRUE,                   /* initial */
                  FALSE,                  /* separate */
                  TRUE,                   /* enforce */
                  TRUE,                   /* check */
                  TRUE,                   /* propagate */
                  FALSE,                  /* local */
                  FALSE,                  /* modifiable */
                  FALSE,                  /* dynamic */
                  FALSE,                  /* removable */
                  FALSE) );               /* stickingatnode */
      SCIP_CALL( SCIPaddCons(scip, con) );
      for (size_t j = 0; j < (size_t)num_nodes; ++j)
      {
         if ( j != i )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, con, i > j ? arc_var[i][j] : arc_var[j][i], 1.0) ); /*lint !e732 !e747*/
         }
      }
      SCIP_CALL( SCIPreleaseCons(scip, &con) );
   }

   /* add set packing constraints (Node 0 is the depot) */
   vector<SCIP_CONS*> part_con(num_nodes, (SCIP_CONS*)NULL);  /*lint !e732 !e747*/
   for (size_t i = 1; i < (size_t)num_nodes; ++i)
   {
      SCIP_CONS* con = NULL;
      (void) SCIPsnprintf(con_name, 255, "C%d", i);
      SCIP_CALL( SCIPcreateConsLinear( scip, &con, con_name, 0, NULL, NULL,
                                       1.0,                /* lhs */
                                       SCIPinfinity(scip), /* rhs */
                                       TRUE,  /* initial */
                                       FALSE, /* separate */
                                       TRUE,  /* enforce */
                                       TRUE,  /* check */
                                       TRUE,  /* propagate */
                                       FALSE, /* local */
                                       TRUE,  /* modifiable */
                                       FALSE, /* dynamic */
                                       FALSE, /* removable */
                                       FALSE  /* stickingatnode */ ) );
      SCIP_CALL( SCIPaddCons(scip, con) );
      part_con[i] = con;  /*lint !e732 !e747*/
   }

   /* include VRP pricer */
   ObjPricerVRP* vrp_pricer_ptr = new ObjPricerVRP(scip, VRP_PRICER_NAME, num_nodes, capacity, demand, dist,
      arc_var, arc_con, part_con);

   SCIP_CALL( SCIPincludeObjPricer(scip, vrp_pricer_ptr, TRUE) );

   /* activate pricer */
   SCIP_CALL( SCIPactivatePricer(scip, SCIPfindPricer(scip, VRP_PRICER_NAME)) );

   //    SCIP_CALL( SCIPwriteOrigProblem(scip, "vrp_init.lp", "lp", FALSE) );


   /*************
    *  Solve    *
    *************/

   SCIP_CALL( SCIPsolve(scip) );


   /**************
    * Statistics *
    *************/
   SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );



   /********************
    * Deinitialization *
    ********************/

   /* release variables */
   for (size_t i = 0; i < (size_t)num_nodes; ++i)
   {
      if ( i > 0 )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &part_con[i]) );
      }
      for (size_t j = 0; j < i; ++j)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &arc_var[i][j]) );
         SCIP_CALL( SCIPreleaseCons(scip, &arc_con[i][j]) );
      }
   }


   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

int main(int argc, char** argv)
{
   return execmain(argc, argv) != SCIP_OKAY ? 1 : 0;
}

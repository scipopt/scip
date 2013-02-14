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
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file
 * @brief  main file for VRP pricer example
 * @author Andreas Bley
 * @author Marc Pfetsch
 *
 * We want to solve the vehicle routing problem on a graph G = (V,E) with
 * V = J cup {d}, where d is the depot and the distances are given by the
 * length function l_e: E -> R_{<= 0}.
 *
 * Consider the MIP formulation
 *
 *   min  sum_{e in E} l_e y_e
 *   s.t. -y_e + sum_{t in T_k} a^t_e x_t <= 0,   for all e in E
 *               sum_{t in T_k} a^t_j x_t == 1,   for all j in J
 *         y(delta(j))                     == 2,   for all j in J
 *         y_e                       in {0,1,2},   for all e in E
 *                              x_t  in [0,1],     for all t in T_k
 *
 * where T_k is the set of tours visiting at most k customers
 * with repetitions of customers allowed and a^t_e (a^t_j) counts how often
 * edge e (node j) is traversed in t in T_k.
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
int read_problem(
   const char*           filename,           /**< filename */
   int&                  num_nodes,          /**< number of nodes in instance */
   int&                  capacity,           /**< capacity in instance */
   vector<int>&          demand,             /**< array of demands of instance */
   vector<vector<int> >& distance            /**< distances between nodes */
   )
{
   static const string DIMENSION           = "DIMENSION";
   static const string DEMAND_SECTION      = "DEMAND_SECTION";
   static const string DEPOT_SECTION       = "DEPOT_SECTION";
   static const string EDGE_WEIGHT_TYPE    = "EDGE_WEIGHT_TYPE";
   static const string EUC_2D              = "EUC_2D";
   static const string EXPLICIT            = "EXPLICIT";
   static const string LOWER_DIAG_ROW      = "LOWER_DIAG_ROW";
   static const string EDGE_WEIGHT_FORMAT  = "EDGE_WEIGHT_FORMAT";
   static const string EDGE_WEIGHT_SECTION = "EDGE_WEIGHT_SECTION";
   static const string NODE_COORD_SECTION  = "NODE_COORD_SECTION";
   static const string CAPACITY            = "CAPACITY";

   ifstream file(filename);

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

         demand.resize(num_nodes, 0);
         distance.resize(num_nodes);
         for (int i = 0; i < num_nodes; ++i)
            distance[i].resize(i, 0);
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
            x.resize(num_nodes, 0);
            y.resize(num_nodes, 0);
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
         for (int i = 0; i < num_nodes; ++i)
         {
            for (int j = 0; j < i; ++j)
            {
               int l;
               file >> l;
               distance[i][j] = l;
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
            x[i] = xi;
            y[i] = yi;
         }
         for (int i = 0; i < num_nodes; ++i)
         {
            for (int j = 0; j < i; ++j)
            {
               int dx = x[i] - x[j];
               int dy = y[i] - y[j];
               distance[i][j] = int( sqrt((double)dx*dx + dy*dy) + 0.5 );
            }
         }
      }
      else if ( key == DEMAND_SECTION )
      {
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
            demand[i] = d;
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
         getline(file, dummy);
      }
   }

   return 0;
}



//------------------------------------------------------------
int main(int argc, char** argv)
{
   SCIP* scip = NULL;

   cout << "Solving the vehicle routing problem using SCIP." << endl;
   cout << "Implemented by Andreas Bley." << endl << endl;

   if ( argc != 2 && argc != 3 )
   {
      cerr << "Usage: vrp [-h] datafile" << endl;
      cerr << "Options:" << endl;
      cerr << " -h  Uses hop limit instead of capacity limit for tours."<< endl;
      return 1;
   }


   /**********************
    * Setup problem data *
    **********************/

   static const char* VRP_PRICER_NAME = "VRP_Pricer";

   vector<vector<int> >  distance;
   vector<int> demand;
   int capacity;
   int num_nodes;

   if ( read_problem(argv[argc-1], num_nodes, capacity, demand, distance) )
   {
      cerr << "Error reading data file " << argv[argc-1] << endl;
      return 1;
   }

   cout << "Number of nodes: " << num_nodes << endl;

   if ( argc == 3 )
   {
      if ( string("-h") != argv[1] )
      {
         cerr << "Unknow option " << argv[2] << endl;
         return 1;
      }

      int total_demand = 0;
      for (int i = 1; i< num_nodes; ++i)
         total_demand += demand[i];
      capacity = (num_nodes - 1) * capacity / total_demand;
      demand.assign(num_nodes, 1);
      demand[0] = 0;
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
   vector< vector<SCIP_VAR*> > arc_var( num_nodes );
   for (int i = 0; i < num_nodes; ++i)
   {
      arc_var[i].resize(i, (SCIP_VAR*) NULL);
      for (int j = 0; j < i; ++j)
      {
         SCIP_VAR* var;
         SCIPsnprintf(var_name, 255, "E%d_%d", i, j );

         SCIP_CALL( SCIPcreateVar(scip,
                     &var,                   // returns new index
                     var_name,               // name
                     0.0,                    // lower bound
                     2.0,                    // upper bound
                     distance[i][j],         // objective
                     SCIP_VARTYPE_INTEGER,   // variable type
                     true,                   // initial
                     false,                  // forget the rest ...
                     0, 0, 0, 0, 0) );
         SCIP_CALL( SCIPaddVar(scip, var) );
         arc_var[i][j] = var;
      }
   }

   /* add arc-routing - tour constraints */
   char con_name[255];
   vector< vector<SCIP_CONS*> > arc_con( num_nodes );
   for (int i = 0; i < num_nodes; ++i)
   {
      arc_con[i].resize(i, (SCIP_CONS*)NULL);
      for (int j = 0; j < i; ++j)
      {
         SCIP_CONS* con;
         SCIPsnprintf(con_name, 255, "A%d_%d", i, j);
         SCIP_VAR* index = arc_var[i][j];
         SCIP_Real coeff = -1;
         SCIP_CALL( SCIPcreateConsLinear(scip, &con, con_name, 1, &index, &coeff,
                     -SCIPinfinity(scip),    /* lhs */
                     0.0,                    /* rhs */
                     true,                   /* initial */
                     false,                  /* separate */
                     true,                   /* enforce */
                     true,                   /* check */
                     true,                   /* propagate */
                     false,                  /* local */
                     true,                   /* modifiable */
                     false,                  /* dynamic */
                     false,                  /* removable */
                     false) );               /* stickingatnode */
         SCIP_CALL( SCIPaddCons(scip, con) );
         arc_con[i][j] = con;
      }
   }

   /* add arc-routing - degree constraints */
   for (int i = 1; i < num_nodes; ++i)
   {
      SCIP_CONS* con;
      SCIPsnprintf(con_name, 255, "D%d", i);
      SCIP_CALL( SCIPcreateConsLinear(scip, &con, con_name, 0, 0, 0,
                  2.0,                    /* lhs */
                  2.0,                    /* rhs */
                  true,                   /* initial */
                  false,                  /* separate */
                  true,                   /* enforce */
                  true,                   /* check */
                  true,                   /* propagate */
                  false,                  /* local */
                  false,                  /* modifiable */
                  false,                  /* dynamic */
                  false,                  /* removable */
                  false) );               /* stickingatnode */
      SCIP_CALL( SCIPaddCons(scip, con) );
      for (int j = 0; j < num_nodes; ++j)
      {
         if ( j != i )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, con, i > j ? arc_var[i][j] : arc_var[j][i], 1.0) );
         }
      }
   }

   /* add set packing constraints (Node 0 is the depot) */
   vector<SCIP_CONS*> part_con(num_nodes, (SCIP_CONS*)NULL);
   for (int i = 1; i < num_nodes; ++i)
   {
      SCIP_CONS* con = NULL;
      SCIPsnprintf(con_name, 255, "C%d", i);
      SCIP_CALL( SCIPcreateConsLinear( scip, &con, con_name, 0, NULL, NULL,
                                       1.0,                /* lhs */
                                       SCIPinfinity(scip), /* rhs */
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
      SCIP_CALL( SCIPaddCons(scip, con) );
      part_con[i] = con;
   }

   /* include VRP pricer */
   ObjPricerVRP* vrp_pricer_ptr = new ObjPricerVRP(scip, VRP_PRICER_NAME, num_nodes, capacity, demand, distance, 
      arc_var, arc_con, part_con);

   SCIP_CALL( SCIPincludeObjPricer(scip, vrp_pricer_ptr, true) );

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
   //SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );



   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return 0;
}

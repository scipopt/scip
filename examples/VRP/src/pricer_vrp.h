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

/**@file pricer_vrp.h
 * @brief VRP pricer plugin
 * @author Andreas Bley
 * @author Marc Pfetsch
 */

#ifndef __SCIP_PRICER_VRP_H__
#define __SCIP_PRICER_VRP_H__

#include "objscip/objscip.h"
#include "scip/pub_var.h"

#include <vector>
#include <list>

using namespace std;
using namespace scip;


/** pricer class */
class ObjPricerVRP : public ObjPricer /*lint --e{3713}*/
{
public:

   /** Constructs the pricer object with the data needed */
   ObjPricerVRP(
      SCIP*                               scip,        /**< SCIP pointer */
      const char*                         p_name,      /**< name of pricer */
      const int                           p_num_nodes, /**< number of nodes */
      const int                           p_capacity,  /**< vehicle capacity */
      const vector< int >&                p_demand,    /**< demand array */
      const vector< vector<int> >&        p_distance,  /**< matrix of distances */
      const vector< vector<SCIP_VAR*> >&  p_arc_var,   /**< matrix of arc variables */
      const vector< vector<SCIP_CONS*> >& p_arc_con,   /**< matrix of arc constraints */
      const vector<SCIP_CONS* >&          p_part_con   /**< array of partitioning constraints */
      );

   /** Destructs the pricer object. */
   virtual ~ObjPricerVRP();

   /** initialization method of variable pricer (called after problem was transformed) */
   virtual SCIP_DECL_PRICERINIT(scip_init);

   /** reduced cost pricing method of variable pricer for feasible LPs */
   virtual SCIP_DECL_PRICERREDCOST(scip_redcost);

   /** farkas pricing method of variable pricer for infeasible LPs */
   virtual SCIP_DECL_PRICERFARKAS(scip_farkas);

   /** perform pricing */
   SCIP_RETCODE pricing(
      SCIP*              scip,               /**< SCIP data structure */
      bool               isfarkas            /**< whether we perform Farkas pricing */
      ) const;

   /** add tour variable to problem */
   SCIP_RETCODE add_tour_variable(
      SCIP*              scip,               /**< SCIP data structure */
      const list<int>&   tour                /**< list of nodes in tour */
      ) const;

   /** return negative reduced cost tour (uses restricted shortest path dynamic programming algorithm) */
   double find_shortest_tour(
      const vector< vector<double> >& length,     /**< matrix of lengths */
      list<int>&         tour                /**< list of nodes in tour */
      ) const;


protected:

   /** return number of nodes */
   inline int num_nodes() const
   {
      return _num_nodes;
   }

   /** return vehicle capacity */
   inline int capacity() const
   {
      return _capacity;
   }

   /** return demand of node i*/
   inline int demand(
      const int          i                   /**< node */
      ) const
   {
      return _demand[i]; /*lint !e747 !e732*/
   }

   /** return distance between nodes i and j */
   inline double distance(
      const int          i,                  /**< first node */
      const int          j                   /**< second node */
      ) const
   {
      return ( i > j ? _distance[i][j] : _distance[j][i] ); /*lint !e747 !e732*/
   }

   /** return variable corresponding to arc between i and j */
   inline SCIP_VAR* arc_var(
      const int          i,                  /**< first node */
      const int          j                   /**< second node */
      ) const
   {
      return ( i > j ? _arc_var[i][j] : _arc_var[j][i] ); /*lint !e747 !e732*/
   }

   /** return constraint corresponding to arc between i and j */
   inline SCIP_CONS* arc_con(
      const int          i,                  /**< first node */
      const int          j                   /**< second node */
      ) const
   {
      return ( i > j ? _arc_con[i][j] : _arc_con[j][i] ); /*lint !e747 !e732*/
   }

   /** return partitioning constraint for node i */
   inline SCIP_CONS* part_con(
      const int          i                   /**< node */
      ) const
   {
      return _part_con[i]; /*lint !e747 !e732*/
   }

   /** whether edge between node i and j exists */
   inline bool have_edge(
      const int          i,                  /**< first node */
      const int          j                   /**< second node */
      ) const
   {
      /* return whether variable is not fixed to 0 */
      return ( SCIPvarGetUbLocal( arc_var(i, j) ) > 0.5 );
   }


private:

   const int                    _num_nodes;  /**< number of nodes */
   const int                    _capacity;   /**< vehicle capacity */
   const vector< int >          _demand;     /**< demand array */
   const vector< vector<int> >  _distance;   /**< distance matrix */

   vector< vector<SCIP_VAR*> >  _arc_var;    /**< matrix of arc variables */
   vector< vector<SCIP_CONS*> > _arc_con;    /**< matrix of arc constraints */
   vector<SCIP_CONS* >          _part_con;   /**< array of partitioning constraints */
};

#endif

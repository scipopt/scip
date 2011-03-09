/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   
 * @brief  VRP pricer plugin
 * @author Andreas Bley
 */


#ifndef __SCIP_PRICER_VRP_H__
#define __SCIP_PRICER_VRP_H__

#include "objscip/objscip.h"

#include "scip/cons_setppc.h"
#include "scip/cons_linear.h"
#include "scip/var.h"

#include <vector>
#include <list>

using namespace std;
using namespace scip;

/* ------------------------------------------------------------ */
class ObjPricerVRP : public ObjPricer 
{
public:
   

   ObjPricerVRP( 
     SCIP*                               scip,
     const char*                         p_name,
     const int                           p_num_nodes,
     const int                           p_capacity,
     const vector< int >&                p_demand,
     const vector< vector<int> >&        p_distance,
     const vector< vector<SCIP_VAR*> >&  p_arc_var,
     const vector< vector<SCIP_CONS*> >& p_arc_con,   
     const vector<SCIP_CONS* >&          p_part_con 
     );
   
   virtual ~ObjPricerVRP();


   virtual SCIP_RETCODE scip_init   ( SCIP* scip, SCIP_PRICER* pricer );

   virtual SCIP_RETCODE scip_redcost( SCIP* scip, SCIP_PRICER* pricer, SCIP_Real* lowerbound, SCIP_RESULT* result );
  
   virtual SCIP_RETCODE scip_farkas ( SCIP* scip, SCIP_PRICER* pricer );


   SCIP_RETCODE
   add_tour_variable
   ( SCIP*              scip, 
     const list<int>&   tour );
      

   /* return negative reduced cost tour 
      (uses restricted shortest path dynamic programming algorithm) */
   double find_shortest_tour
   ( const vector< vector<double> >& length,
     list<int>&                      tour );
   
protected:
   
   inline int num_nodes() const
   {
      return _num_nodes;
   }

   inline int capacity() const
   {
      return _capacity;
   }

   inline int demand( const int i ) const
   {
      return _demand[i];
   }

   inline double distance( const int i, const int j ) const
   {
      return ( i>j ? _distance[i][j] : _distance[j][i] );
   }

   inline SCIP_VAR* arc_var( const int i, const int j ) const
   {
      return ( i>j ? _arc_var[i][j] : _arc_var[j][i] );
   }

   inline SCIP_CONS* arc_con( const int i, const int j ) const
   {
      return ( i>j ? _arc_con[i][j] : _arc_con[j][i] );
   }

   inline SCIP_CONS* part_con( const int i ) const
   {
      return ( _part_con[i] );
   }
   
   inline bool have_edge( const int i, const int j ) const
   {
      return ( SCIPvarGetUbLP( arc_var( i, j ) ) >= 1 );
   }



private:

   const int                    _num_nodes;
   const int                    _capacity;
   const vector< int >          _demand;
   const vector< vector<int> >  _distance;

   vector< vector<SCIP_VAR*> >  _arc_var;   
   vector< vector<SCIP_CONS*> > _arc_con;   
   vector<SCIP_CONS* >          _part_con;   
};


#endif


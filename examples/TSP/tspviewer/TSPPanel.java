/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
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

/**
 * @author Timo Berthold
 *
 * Panel which can visualize a TSP tour
 */

import javax.swing.JPanel;
import java.awt.Color;
import java.awt.Graphics;

public class TSPPanel extends JPanel 
{
   private double[][] coords;
   private int diam = 8;

   // scales the coordinates of the tour such that it can be drawn into the panel
   private int[][] scaleTour(int xs, int ys)
   {
      int xsize = xs - diam;
      int ysize = ys - diam;

      int nnodes = coords[0].length;
      int[][] tour = new int[2][nnodes];
      double[] borders = new double[4]; // the extreme values of the coords
		
      // initialize 
      for( int i = 0; i < 3; i++ )
      {
         if( i % 2 == 0 )
            borders[i] = Double.MAX_VALUE;
         else
            borders[i] = Double.MIN_VALUE;
      }

      // find the extreme values of the coords
      for( int i = 0; i < nnodes; i++ )
      {
         if( coords[0][i] < borders[0] )
            borders[0] = coords[0][i];
         if( coords[0][i] > borders[1] )
            borders[1] = coords[0][i];
         if( coords[1][i] < borders[2] )
            borders[2] = coords[1][i];
         if( coords[1][i] > borders[3] )
            borders[3] = coords[1][i];
      }
		
      // calculate the scaling factor
      double x_scale = Math.max(borders[1] - borders[0], 1.0);
      double y_scale = Math.max(borders[3] - borders[2], 1.0);
      double scale = Math.min(xsize/x_scale, ysize/y_scale);
		
      // calculate the shift		
      double x_off = 0.0;
      double y_off = 0.0;
      if( xsize/x_scale < ysize/y_scale )
      {
         x_off = diam/2;
         y_off = diam/2 + (ysize - scale*y_scale)/2;
      }
      else
      {
         x_off = diam/2 + (xsize - scale*x_scale)/2;
         y_off = diam/2;
      }
		
      // scale the coordinates
      for( int i = 0; i < nnodes; i++ )
      {
         tour[0][i] =  (int) Math.round(x_off + scale * (coords[0][i] - borders[0]));
         tour[1][i] =  (int) Math.round(ys - (y_off + scale * (coords[1][i] - borders[2])));
      }
		
      return tour;
   }
	
   // overwritten paint method, scales tour and paints it into the panel
   protected void paintComponent(Graphics g) 
   {		
      // clear panel
      g.setColor(Color.white);
      g.fillRect(0, 0, getSize().width, getSize().height);
      if( coords == null )
         return;

      int[][] tour = scaleTour(getSize().width, getSize().height);
      int nnodes = tour[0].length;
		
      // draw red lines connecting successive nodes 
      g.setColor(Color.red);
      for( int i = 0; i < nnodes; i++ )
      {
         g.drawLine(tour[0][i], tour[1][i], tour[0][(i+1) % nnodes], tour[1][(i+1) % nnodes]);
      }

      // draw a yellow circle with orange boundary for each node
      for( int i = 0; i < nnodes; i++ )
      {
         g.setColor(Color.yellow);
         g.fillOval(tour[0][i] - diam/2, tour[1][i] - diam/2, diam, diam);
         g.setColor(Color.orange);
         g.drawOval(tour[0][i] - diam/2, tour[1][i] - diam/2, diam, diam);
      }
   }
	
	
   /**
    * @return Returns the tour.
    */
   public double[][] getTour() {
      return coords;
   }
   /**
    * @param tour The tour to set.
    */
   public void setTour(double[][] coords) {
      this.coords = coords;
   }
}

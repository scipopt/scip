import javax.swing.JPanel;
import java.awt.Color;
import java.awt.Graphics;
/*
 * Created on 01.03.2005
 *
 */

/**
 * @author Timo Berthold
 *
 * Panel which can visualize a TSP tour
 */
public class TSPPanel extends JPanel 
{
	
	private double[][] coords;
	
	// scales the coordinates of the tour such that it can be drawn into the panel
	private int[][] scaleTour(int xs, int ys)
	{
		int xsize = xs - 10;
		int ysize = ys - 10;

		int nnodes = coords[0].length;
		int[][] tour = new int[2][nnodes];
		int i;
		double[] borders = new double[4]; //the extreme values of the coords
		
		//initialize 
		for(i = 0; i < 3; i++)
			if(i % 2 == 0)
				borders[i] = Double.MAX_VALUE;
			else
				borders[i] = Double.MIN_VALUE;
	
		//find the extreme values of the coords
		for( i = 0; i < nnodes; i++ )
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
		double scale = Math.min(xsize/x_scale, ysize/y_scale );
		
		tour = new int[2][nnodes];
		
		// calculate the shift		
		double x_off = 0.0;
		double y_off = 0.0;
		
		if(xsize/x_scale < ysize/y_scale)
		{
			x_off = - borders[0];
			y_off = (ysize - scale*y_scale)/2 - borders[2];
		}
		else
		{
			y_off = - borders[2];
			x_off = (xsize - scale*x_scale)/2 - borders[0];
		}
		
		// scale the coordinates
		for( i = 0; i < nnodes; i++ )
		{
			tour[0][i] =  (int) Math.round( 5 + x_off + scale * ( coords[0][i] - borders[0] ) );
			tour[1][i] =  (int) Math.round( ysize - ( 5 + y_off + scale * ( coords[1][i] - borders[2] ) ) );
		}
		
		return tour;
	}
	
	// overwritten paint method, scales tour and paints it into the panel
	protected void paintComponent(Graphics g) 
	{		
		// clear panel
		g.setColor(Color.white);
		g.fillRect(0,0,getSize().width,getSize().height);
		if(coords == null)
			return;
		int[][] tour = scaleTour(getSize().width,getSize().height);
		int nnodes = tour[0].length;
		int[] x = new int[nnodes];
		int[] y = new int[nnodes];
		int d = 8;
		int i;
		
		for(i = 0; i < nnodes; i++)
		{		
			x[i] = tour[0][i];
			y[i] = tour[1][i];
		}
		
		// draw red lines connecting successive nodes 
		g.setColor(Color.red);
		for(i = 0; i < nnodes; i++)
			g.drawLine(x[i] + 4, y[i] + 4, x[ (i+1) % nnodes] + 4, y[ (i+1) % nnodes ] + 4);
		
		// draw a yellow circle with orange boundary for each node
		for(i = 0; i < nnodes; i++)
		{
			g.setColor(Color.yellow);
			g.fillOval(x[i], y[i], d, d);
			g.setColor(Color.orange);
			g.drawOval(x[i], y[i], d, d);
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

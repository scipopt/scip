import java.awt.Graphics;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.StringTokenizer;
import javax.swing.JFrame;
import java.awt.GridLayout;

/*
 * Created on 02.03.2005
 *
 */

/**
 * @author Timo Berthold
 *
 */
public class TSPViewer extends Thread {
	
	private static File file;
	private static File lock; //lock file forbidding other programms to change the file during it is read
	
	public TSPViewer(String s)
	{
		file = new File(s);	
		lock = new File(s+".lock");
	}
	
	// read in a file
	private static double[][] readFile()
	{
		
		int nnodes;
		int i;
		
		double[][] coords = new double[2][];
		
		try
		{
			String line;
			String token;
			StringTokenizer strtok;
			int count = 0;
			
			// if the file is in use, wait
			while( lock.exists() && count < 10 )
			{
				sleep(500);
				count++;
			}
			if( lock.exists() )
			{
				System.out.println("Cannot read file <temp.tour>, because lockfile <temp.lock> still exists");
				return null;
			}
			
			//create lock file forbidding other programms to use the file while we read
			lock.createNewFile();
			BufferedReader in = new BufferedReader( new FileReader(file) );
			if( (line = in.readLine()) != null)
				nnodes = Integer.parseInt(line);
			else	
				nnodes = 0;
			
			coords[0] = new double[nnodes];
			coords[1] = new double[nnodes];
			
			i = 0;	
			// read out the coordinates
			while ( i < nnodes && (line = in.readLine()) !=  null)		
			{         
				strtok = new StringTokenizer(line);	
				strtok.nextToken();
				coords[0][i] = Double.parseDouble(strtok.nextToken());
				coords[1][i] = Double.parseDouble(strtok.nextToken());		
				
				i++;
			}
			in.close();
			lock.delete();
			
		} 
		catch (Exception e)
		{
			System.err.println("TSP Viewer: TSPFrame: File input error");
		}
		return coords;	
	}
	
	/** 
	 * monitors whether the input file has been changed and if so, repaints.
	 * */
	public void run()
	{
		long updated = file.lastModified();
		
		//Create and set up the frame
		JFrame.setDefaultLookAndFeelDecorated(true);
        JFrame frame = new JFrame("TSP Viewer");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(800,600);
        Graphics g = frame.getGraphics();
        
        //Add the components
        TSPPanel panel = new TSPPanel();             
        frame.getContentPane().setLayout(new GridLayout(1,1));
        frame.getContentPane().add(panel);
        panel.paint(g);
        //Display the window
        frame.setVisible(true);
        updated = -1;
        //keep active as ong as the frame exists
        while(frame != null)
        {
        	long upd = updated;
        	if(file.exists())
        		upd = file.lastModified();
        	// check, whether file has been updated. If so, repaint, otherwise sleep for half a second
        	if( upd != updated )
        	{
        		updated = upd;
        		double[][] tour = readFile();
        		if(tour==null)
        			updated = -1;
        		else
        			panel.setTour(tour);  
        		frame.repaint();
        	} 
        	else
				try 
				{
					sleep(500);
				} 
        		catch (InterruptedException e) 
				{
					e.printStackTrace();
				}
        }
	}
	
	/**
	 * creates one Frame and makes it run
	 * @param args The name of the file to display
	 */
	public static void main(String[] args) 
	{
		TSPViewer window = new TSPViewer(args[0]);
		window.run();
	}
}

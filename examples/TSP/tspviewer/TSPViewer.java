/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**
 * @author Timo Berthold
 * Frame displaying solutions of a TSP and updating from file
 */

import java.awt.Graphics;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.StringTokenizer;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JLabel;
import javax.swing.JButton;
import java.awt.BorderLayout;

public class TSPViewer extends Thread {
	
   private static File file;
   private static File lock; //lock file forbidding other programms to change the file during it is read
   private int actionsteps = 0;
   
   public TSPViewer(String s)
   {
      file = new File(s);	
      lock = new File(s+".lock");
      System.out.println("Watching file <" + file.getPath() + ">");
   }
	
   // read in a file
   private static TSPSolution readFile(int number)
   {
		
      int nnodes = 0;
      int i;
      String name = "";
      double obj = -1.0;
	  
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
            System.out.println("Cannot read file <" + file.getPath() + ">, because lockfile <"
               + lock.getPath() + "> still exists");
            return null;
         }
			
         // create lock file forbidding other programms to use the file while we read
         lock.createNewFile();
         BufferedReader in = new BufferedReader(new FileReader(file));
         
         if( (line = in.readLine()) != null )
         {
            if( line.equals("RESET") )
            {
               lock.delete();
               in.close();
               return null;
            }
            else	
               nnodes = Integer.parseInt(line); 
         }
         if( ((line = in.readLine()) != null ))
            name = new String(line);
         if( ((line = in.readLine()) != null ))
            obj = Double.parseDouble(line);
         
         
         coords[0] = new double[nnodes];
         coords[1] = new double[nnodes];
			
         // read out the coordinates
         i = 0;	
         while( i < nnodes && (line = in.readLine()) != null )
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
      catch( Exception e )
      {
         System.err.println("TSP Viewer: TSPFrame: File input error");
      }

      return new TSPSolution(coords, obj, number, name);	
   }
  
   /** 
    * monitors whether the input file has been changed and if so, repaints.
    */
   public void run()
   {
      long updated = file.lastModified();
	  	
	  
      // create and set up the frame
      JFrame.setDefaultLookAndFeelDecorated(true);
      JFrame frame = new JFrame("TSP Viewer");
      frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
      frame.setSize(800, 600);
      Graphics g = frame.getGraphics();
        
      // add the components
      TSPPanel tsppanel = new TSPPanel();    
      JPanel panel = new JPanel();
      JButton prev = new JButton("<<");
      prev.addActionListener( new ActionListener() {
            public void actionPerformed(ActionEvent e) {
               actionsteps--;
            }
         } );
      prev.setEnabled(false);
      JButton next = new JButton(">>");
      next.addActionListener( new ActionListener() {
            public void actionPerformed(ActionEvent e) {
               actionsteps++;
            }
         } );
      next.setEnabled(false);
      JLabel label = new JLabel("No solution found yet");
      
      panel.add(prev);
      panel.add(next);
      panel.add(label);
      
      frame.getContentPane().setLayout(new BorderLayout());
      frame.getContentPane().add(tsppanel, BorderLayout.CENTER);
      frame.getContentPane().add(panel, BorderLayout.SOUTH);
      tsppanel.paint(g);

      // display the window
      frame.setVisible(true);

      
      LinkedList solutions = new LinkedList();
      int currsol = 0;
      updated = -1;
      int nsols = 0;
      // keep active as long as the frame exists	
      while( frame != null )
      {
         // check, whether file has been updated. If so, repaint
         long upd = updated;
         if( file.exists() )
            upd = file.lastModified();
         if( upd != updated )
         {
            updated = upd;
            TSPSolution sol = readFile(nsols+1);
            
            if( sol == null )
            {
               nsols = 0;
               currsol = 0;
               next.setEnabled(false);
               prev.setEnabled(false);          
               label.setText("no solution found yet");
               tsppanel.setTour(null);
               solutions = new LinkedList();
            }
            else
            {
               solutions.addLast( sol );
               tsppanel.setTour(sol.getTourcoords());
               label.setText(sol.getNumber()+". solution, found by "+sol.getHeur()+" , objective: "+sol.getObjval());
               currsol = nsols;
               nsols++;
               if(nsols > 1)
                  prev.setEnabled(true);
               next.setEnabled(false);
               actionsteps = 0;
            }             
            frame.repaint();
         } 
         // check if action was performed on the buttons. If so, repaint
         else if( actionsteps != 0 )
         {
            currsol = Math.max(0, Math.min(currsol+actionsteps, nsols-1));
            ListIterator it = solutions.listIterator(currsol);
            TSPSolution sol = (TSPSolution) it.next();
            tsppanel.setTour(sol.getTourcoords());
            /*double[][] t = sol.getTourcoords();
              System.out.print("Painting tour number " + sol.getNumber()+": ");
              for(int j = 0; j < t.length; j++)
              System.out.print("("+t[j][0]+"|"+t[j][1]+") - ");
              System.out.print("\n");*/
            label.setText(sol.getNumber()+". solution, found by "+sol.getHeur()+" , objective: "+sol.getObjval());
            if( currsol != 0 )
               prev.setEnabled(true);
            else
               prev.setEnabled(false);
            if( currsol != nsols -1 )
               next.setEnabled(true);
            else
               next.setEnabled(false);
            actionsteps = 0;
            frame.repaint();
         }
         // if nothing has happened, sleep a bit
         else
         {
            try 
            {
               sleep(500);
            } 
            catch( InterruptedException e )
            {
               e.printStackTrace();
            }
         }
      }
   }
	
   /**
    * creates one Frame and makes it run
    * @param args The name of the file to display
    */
   public static void main(String[] args) 
   {
      String filename = "../temp.tour";

      if( args.length >= 1 )
         filename = args[0];

      TSPViewer window = new TSPViewer(filename);
      window.run();
   }
}

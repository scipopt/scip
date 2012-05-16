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
 *
 * Object just storing information about a TSP solution
 */

public class TSPSolution {
	
	private double[][] tourcoords;
    private	double objval;
	private String heur;
	private int number;
	

	/**
	 * @param tourcoords
	 * @param objval
	 * @param heur
	 */
	public TSPSolution(double[][] tourcoords, double objval, int number, String heur) {
		super();
		this.tourcoords = tourcoords;
		this.objval = objval;
		this.heur = heur;
		this.number = number;
	}
	
	
	/**
	 * @return Returns the heur.
	 */
	public String getHeur() {
		return heur;
	}
	/**
	 * @param heur The heur to set.
	 */
	public void setHeur(String heur) {
		this.heur = heur;
	}
	/**
	 * @return Returns the objval.
	 */
	public double getObjval() {
		return objval;
	}
	/**
	 * @param objval The objval to set.
	 */
	public void setObjval(double objval) {
		this.objval = objval;
	}
	/**
	 * @return Returns the tourcoords.
	 */
	public double[][] getTourcoords() {
		return tourcoords;
	}
	/**
	 * @param tourcoords The tourcoords to set.
	 */
	public void setTourcoords(double[][] tourcoords) {
		this.tourcoords = tourcoords;
	}
	
	/**
	 * @return Returns the number.
	 */
	public int getNumber() {
		return number;
	}
	/**
	 * @param number The number to set.
	 */
	public void setNumber(int number) {
		this.number = number;
	}
}

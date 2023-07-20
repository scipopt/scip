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

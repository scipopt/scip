/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the examples to                     */
/*                         An introduction to SCIP                           */
/*                                                                           */
/*    Copyright (C) 2007 Cornelius Schwarz                                   */
/*                                                                           */
/*                  2007 University of Bayreuth                              */
/*                                                                           */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**
 * @file queens.cpp
 * @author Cornelius Schwarz
 * @brief n-queens examlple implementation
 */

#include "queens.hpp"
#include <sstream>
#include "scip_exception.hpp"

using namespace std;
using namespace scipexamples;

/* constructor */
scipexamples::QueensSolver::QueensSolver(int n)
   : _scip(0), _n(n), _vars(n, vector<SCIP_VAR *>(n)), _cons()
{
   // initialize scip
   SCIP_CALL_EXC( SCIPcreate(& _scip) );

   // load default plugins linke separators, heuristics, etc.
   SCIP_CALL_EXC( SCIPincludeDefaultPlugins(_scip) );

   // disable scip output to stdout
   SCIP_CALL_EXC( SCIPsetMessagehdlr(_scip, NULL) );

   // create an empty problem
   SCIP_CALL_EXC( SCIPcreateProb(_scip, "queens", NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   // set the objective sense to maximize, default is minimize
   SCIP_CALL_EXC( SCIPsetObjsense(_scip, SCIP_OBJSENSE_MAXIMIZE) );

   // create a binary variable for every field (i,j) on the chess board
   ostringstream namebuf;
   for (int i = 0; i < _n; ++i)
   {
      for(int j = 0; j < _n; ++j)
      {
	 SCIP_VAR * var;
	 namebuf.str("");
	 namebuf << "x#" << i << "#" << j;

	 // create the SCIP_VAR object
	 SCIP_CALL_EXC( SCIPcreateVar(_scip, & var, namebuf.str().c_str(), 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );

	 // add the SCIP_VAR object to the scip problem
	 SCIP_CALL_EXC( SCIPaddVar(_scip, var) );

	 // storing the SCIP_VAR pointer for later access
	 _vars[i][j] = var;
      }
   }

   // create constraints
   // one queen per row
   for (int i = 0; i < _n; ++i)
   {
      SCIP_CONS * cons;
      namebuf.str("");
      namebuf<<"row_"<<i;

      // create SCIP_CONS object
      // this is an equality since there must be a queen in every row
      SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, & cons, namebuf.str().c_str(), 0, NULL, NULL, 1.0, 1.0,
					  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      // add the vars belonging to field in this row to the constraint
      for (int j = 0; j < _n; ++j)
	 SCIP_CALL_EXC( SCIPaddCoefLinear(_scip, cons, _vars[i][j], 1.0) );

      // add the constraint to scip
      SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );

      // store the constraint for later on
      _cons.push_back(cons);
   }

   // this is the same with the col constraints ( one queen per column)
   for (int j = 0; j < _n; ++j)
   {
      SCIP_CONS * cons;
      namebuf.str("");
      namebuf << "col_" << j;
      SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, & cons, namebuf.str().c_str(), 0, NULL, NULL, 1.0, 1.0,
					  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      for (int i = 0; i < _n; ++i)
	 SCIP_CALL_EXC( SCIPaddCoefLinear(_scip, cons, _vars[i][j], 1.0) );

      SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );
      _cons.push_back(cons);
   }

   // now we create the constraints for the diagonals
   // there is only one queen allowed in each diagonal, but there do not need to be one. Therefore we add a <= 1 constraint
   // in this problem case we can set the lhs to zero instead of -SCIPinfinity
   // diag col down
   for (int j = 0; j < _n; ++j)
   {
      SCIP_CONS * cons;
      namebuf.str("");
      namebuf << "diag_col_down_" << j;
      SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, & cons, namebuf.str().c_str(), 0, NULL, NULL, 0.0, 1.0,
					 TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      for (int i = 0; i < _n - j; ++i)
	 SCIP_CALL_EXC( SCIPaddCoefLinear(_scip, cons, _vars[i][j+i], 1.0) );
      SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );
      _cons.push_back(cons);
   }

   // diag row down
   for (int i = 0; i < _n; ++i)
   {
      SCIP_CONS * cons;
      namebuf.str("");
      namebuf<<"diag_row_down_"<<i;
      SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, & cons, namebuf.str().c_str(), 0, NULL, NULL, 0.0, 1.0,
					  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      for (int j = 0; j < _n - i; ++j)
	 SCIP_CALL_EXC( SCIPaddCoefLinear(_scip, cons, _vars[i+j][j], 1.0) );
      SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );
      _cons.push_back(cons);
   }

   // diag col up
   for (int j = 0; j < _n; ++j)
   {
      SCIP_CONS * cons;
      namebuf.str("");
      namebuf<<"diag_col_up_"<<j;
      SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, & cons, namebuf.str().c_str(), 0, NULL, NULL, 0.0, 1.0,
					  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));

      for (int i = 0; i < _n - j; ++i)
	 SCIP_CALL_EXC( SCIPaddCoefLinear(_scip, cons, _vars[i][_n - j - i - 1], 1.0) );
      SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );
      _cons.push_back(cons);
   }

   // diag row up
   for (int i = 0; i < _n; ++i)
   {
      SCIP_CONS * cons;
      namebuf.str("");
      namebuf<<"diag_row_up"<<i;
      SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, & cons, namebuf.str().c_str(), 0, NULL, NULL, 0.0, 1.0,
					  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      for (int j = 0; j < _n - i; ++j)
	 SCIP_CALL_EXC( SCIPaddCoefLinear(_scip, cons, _vars[i+j][_n - j - 1], 1.0) );
      SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );
      _cons.push_back(cons);
   }
}


/* display the solution */
void scipexamples::QueensSolver::disp(std::ostream& out)
{
   // get the best found solution from scip
   SCIP_SOL * sol = SCIPgetBestSol(_scip);
   out << "solution for " << _n << "-queens:" << endl << endl;

   // when SCIP did not succeed then sol is NULL
   if (sol == NULL)
   {
      out << "no solution found" << endl;
      return;
   }

   for (int i = 0; i < _n; ++i)
   {
      for (int j = 0; j < _n; ++j)
	 out << " ---";
      out << endl;

      for (int j = 0; j < _n; ++j)
      {
	 out << "| ";
	 // we display a D in every field if x[i][j] = 1 which means that a queen will be placed there
	 if ( SCIPgetSolVal(_scip, sol, _vars[i][j]) > 0.5 )
	    out << "D ";
	 else
	    out << "  ";
      }
      out << "|" << endl;
   }
   for (int j = 0; j < _n; ++j)
      out << " ---";
   out << endl;
}


/* destructor */
scipexamples::QueensSolver::~QueensSolver(void)
{
   // since the SCIPcreateVar captures all variables, we have to release them now
   for (int i = 0; i < _n; ++i)
   {
      for(int j = 0; j < _n; ++j)
	 SCIP_CALL_EXC( SCIPreleaseVar(_scip, & _vars[i][j]) );
   }
   _vars.clear();

   // the same for the constraints
   for (vector<SCIP_CONS *>::size_type i = 0; i < _cons.size(); ++i)
      SCIP_CALL_EXC( SCIPreleaseCons(_scip, &_cons[i]) );
   _cons.clear();

   // after releasing all vars and cons we can free the scip problem
   // remember this has allways to be the last call to scip
   SCIP_CALL_EXC( SCIPfree( & _scip) );
}

/* solve the n-queens problem */
void QueensSolver::solve(void)
{
   // this tells scip to start the solution process
   SCIP_CALL_EXC( SCIPsolve(_scip) );
}

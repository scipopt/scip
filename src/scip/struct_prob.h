/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_prob.h,v 1.9 2004/04/27 15:50:05 bzfpfend Exp $"

/**@file   struct_prob.h
 * @brief  datastructures for storing and manipulating the main problem
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_PROB_H__
#define __STRUCT_PROB_H__


#include "def.h"
#include "type_misc.h"
#include "type_var.h"
#include "type_prob.h"
#include "type_cons.h"


/** main problem to solve */
struct Prob
{
   Real             objoffset;          /**< objective offset from bound shifting and fixing (fixed vars result) */
   Real             objlim;             /**< objective limit as external value */
   char*            name;               /**< problem name */
   DECL_PROBDELORIG ((*probdelorig));   /**< frees user data of original problem */
   DECL_PROBTRANS   ((*probtrans));     /**< creates user data of transformed problem by transforming original user data */
   DECL_PROBDELTRANS((*probdeltrans));  /**< frees user data of transformed problem */
   DECL_PROBINITSOL ((*probinitsol));   /**< solving process initialization method of transformed data */
   DECL_PROBEXITSOL ((*probexitsol));   /**< solving process deinitialization method of transformed data */
   PROBDATA*        probdata;           /**< user problem data set by the reader */
   HASHTABLE*       varnames;           /**< hash table storing variable's names */
   VAR**            fixedvars;          /**< array with fixed and aggregated variables */
   VAR**            vars;               /**< array with active variables ordered binary, integer, implicit, continuous */
   HASHTABLE*       consnames;          /**< hash table storing constraints' names */
   CONS**           conss;              /**< array with constraints of the problem */
   int              fixedvarssize;      /**< available slots in fixedvars array */
   int              nfixedvars;         /**< number of fixed and aggregated variables in the problem */
   int              varssize;           /**< available slots in vars array */
   int              nvars;              /**< number of mutable variables in the problem (used slots in vars array) */
   int              nbinvars;           /**< number of binary variables */
   int              nintvars;           /**< number of general integer variables */
   int              nimplvars;          /**< number of implicit integer variables */
   int              ncontvars;          /**< number of continuous variables */
   int              ncolvars;           /**< number of variables with attached column information */
   int              consssize;          /**< available slots in conss array */
   int              nconss;             /**< number of constraints in the problem (number of used slots in conss array) */
   int              maxnconss;          /**< maximum number of constraints existing at the same time */
   int              startnconss;        /**< number of constraints existing when problem solving started */
   OBJSENSE         objsense;           /**< objective sense */
   Bool             objisintegral;      /**< is objective value always integral for feasible solutions? */
   Bool             transformed;        /**< TRUE iff problem is the transformed problem */
};


#endif

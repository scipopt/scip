/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_sepastore.h,v 1.2 2003/12/15 17:45:34 bzfpfend Exp $"

/**@file   struct_sepastore.h
 * @brief  datastructures for storing separated cuts
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_SEPASTORE_H__
#define __STRUCT_SEPASTORE_H__


#include "def.h"
#include "type_lp.h"
#include "type_var.h"
#include "type_sepastore.h"


/** storage for separated cuts
 *  cuts with only one variable are treated specifically as bound changes, and they are not counted for exceeding the
 *  separation limit
 */
struct SepaStore
{
   ROW**            cuts;               /**< array with separated cuts sorted by score */
   Real*            score;              /**< score for each separated cut (e.g. violation/(eucnorm * #nonzeros)) */
   VAR**            bdchgvars;          /**< variables to change bound for */
   Real*            bdchgvals;          /**< new bounds of variables */
   BOUNDTYPE*       bdchgtypes;         /**< type of bounds to change */
   int              cutssize;           /**< size of cuts and score arrays */
   int              ncuts;              /**< number of separated cuts (max. is set->maxsepacuts) */
   int              bdchgssize;         /**< size of bdchgvars, bdchgvals, and bdchgtypes arrays */
   int              nbdchgs;            /**< number of bound changes */
   int              ncutsfound;         /**< total number of cuts (and bound changes) found so far */
   int              ncutsapplied;       /**< total number of cuts (and bound changes) applied to the LPs */
   Bool             initiallp;          /**< is the separation storage currently being filled with the initial LP rows? */
};


#endif

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cmain.c
 * @brief  main file for C compilation
 * @author Tobias Achterberg
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "scip.h"
#include "cons_linear.h"


static const int   nrows = 2;
static const int   ncols = 2;

static const char* row_name[] = { "lin1"  , "lin2"   };
static const int   row_len [] = {        2,        2 };
static const int   row_idx [] = {   0,   1,   0,   1 };
static const Real  row_val [] = { 2.0, 1.0, 1.0, 2.0 };
static const Real  row_rhs [] = {      2.0,      2.0 };

int
main(void)
{
   SCIP* scip = NULL;
   RETCODE retcode;
   COL** cols;
   COL** rowcols;
   Real* rowvals;
   CONS* cons;
   int r;
   int c;
   int i;
   int pos;
   char cname[255];

   printf("SCIP version %g\n", SCIP_VERSION);

   allocMemoryArray(cols, ncols);
   allocMemoryArray(rowcols, ncols);
   allocMemoryArray(rowvals, ncols);

   CHECK_SCIP( SCIPcreate(&scip) );
   CHECK_SCIP( SCIPcreateProb(scip, "test.lp") );

#ifndef NDEBUG
   SCIPdebugMemory(scip);
#endif

   for( c = 0; c < ncols; ++c )
   {
      sprintf(cname, "x%d", c);
      CHECK_SCIP( SCIPcreateCol(scip, &cols[c], cname, 0.0, 10.0, -1.0, SCIP_COLTYPE_INTEGER) );
   }

   pos = 0;
   for( r = 0; r < nrows; ++r )
   {
      for( i = 0; i < row_len[r]; ++i )
      {
         rowcols[i] = cols[row_idx[pos]];
         rowvals[i] = row_val[pos];
         pos++;
      }
      CHECK_SCIP( SCIPconsCreate_Linear(scip, &cons, row_name[r], row_len[r], rowcols, rowvals, row_rhs[r], 0.0, 1e-6,
                     SCIP_ROWTYPE_LESSEQUAL, TRUE, TRUE) );
   }

   CHECK_SCIP( SCIPsolve(scip) );

#ifndef NDEBUG
   SCIPdebugMemory(scip);
#endif

   CHECK_SCIP( SCIPfree(&scip) );

   freeMemoryArray(rowvals);
   freeMemoryArray(rowcols);
   freeMemoryArray(cols);

#ifndef NDEBUG
   memoryCheckEmpty();
#endif

   return 0;
}


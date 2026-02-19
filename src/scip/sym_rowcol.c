/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   sym_rowcol.c
 * @ingroup DEFPLUGINS_SYM
 * @brief  symmetry handler for row and column symmetries
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_linear.h"
#include "scip/cons_orbitope.h"
#include "scip/event_shadowtree.h"
#include "scip/pub_cons.h"
#include "scip/pub_implics.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_sym.h"
#include "scip/pub_var.h"
#include "scip/scip_cons.h"
#include "scip/scip_event.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sym.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"
#include "scip/sym_rowcol.h"
#include "scip/symmetry.h"
#include "scip/symmetry_lexred.h"
#include "scip/symmetry_orbitopal.h"
#include "scip/type_implics.h"
#include "scip/type_set.h"

/* symmetry handler properties */
#define SYM_NAME            "sym_rowcol"
#define SYM_DESC            "symmetry handler for row and column symmetries"
#define SYM_PRIORITY           -100000       /**< priority of try-add function*/
#define SYM_PROPPRIORITY       -100000       /**< priority of propagation method */
#define SYM_PROPTIMING SCIP_PROPTIMING_BEFORELP /**< timing of propagator method */
#define SYM_PROPFREQ                 1       /**< frequence of propagator method */
#define SYM_DELAYPROP            FALSE       /**< Should propagation method be delayed, if other propagators found reductions? */

/* default value of parameters */
#define SYM_DEFAULT_DETECTDOUBLELEX TRUE     /**< Should we check whether the components can be handled by double lex matrices? */
#define SYM_DEFAULT_HANDLESIGNEDORBITOPES TRUE /**< Shall signed orbitopes be handled? */
#define SYM_DEFAULT_USEDYNAMICPROP TRUE      /**< Shall dynamic propagation be used for orbitopes? */


/*
 * Data structures
 */

/** symmetry component data */
struct SCIP_SymCompData
{
   SCIP_CONS**           conss;              /**< static symmetry handling constraints */
   int                   nconss;             /**< number of static symmetry handling constraints */
   int                   maxnconss;          /**< maximum number of constraints conss can hold */
   SCIP_LEXREDDATA*      lexreddata;         /**< container for lexicographic reduction propagation; */
   SCIP_Bool             lexredactive;       /**< whether lexicographic reduction is active on this component */
};

/** symmetry handler data */
struct SCIP_SymhdlrData
{
   SCIP_EVENTHDLR*       shadowtreeeventhdlr; /**< shadow tree event handler */
   SCIP_ORBITOPALREDDATA* orbitopalreddata;  /**< container for orbitopal reduction data */
   SCIP_Bool             detectdoublelex;    /**< Should we check whether the components can be handled by double lex matrices? */
   SCIP_Bool             handlesignedorbitopes;/**< Shall signed orbitopes be handled? */
   SCIP_Bool             usedynamicprop;     /**< Shall dynamic propagation be used for orbitopes? */
};

/*
 * Local methods
 */

/** orbitopal reduction */
static
SCIP_RETCODE addOrbitopesDynamic(
   SCIP*                 scip,               /**< SCIP instance */
   SYM_SYMTYPE           symtype,            /**< symmetry type */
   int**                 perms,              /**< permutations */
   SCIP_VAR**            permvars,           /**< variables the symmetries act on */
   int                   npermvars,          /**< number of variables */
   SCIP_Real*            permvardomaincenter,/**< domain center of variables */
   int                   componentid,        /**< ID of component for which orbitope is added */
   char*                 partialname,        /**< partial name for orbitope constraint */
   int**                 varidxmatrix,       /**< matrix containing variable indices in orbitope matrix */
   int                   nrows,              /**< number of rows of orbitope */
   int                   ncols,              /**< number of columns of orbitope */
   SCIP_CONS***          orbitopeconss,      /**< pointer to store array of conss derived from orbitopes */
   int*                  norbitopeconss,     /**< pointer to store number of conss in orbitopeconss */
   int*                  maxnorbitopeconss,  /**< pointer to store maximum number of conss orbitopeconss can hold */
   SCIP_ORBITOPALREDDATA* orbitopalreddata,  /**< data needed for orbitopal reduction */
   SCIP_LEXREDDATA*      lexreddata,         /**< data needed for lexicographic reduction */
   SCIP_Bool*            lexredactive,       /**< pointer to store whether lexicographic reduction is active */
   SCIP_Bool*            success             /**< pointer to store whether orbitope could be added successfully */
   )
{
   char name[SCIP_MAXSTRLEN];
   int i;
   int j;
   SCIP_Bool ispporbitope;
   SCIP_VAR*** varmatrix;
   SCIP_Bool* pprows;
   int npprows;
   SCIP_ORBITOPETYPE type;

   assert(scip != NULL);
   assert(perms != NULL);
   assert(permvars != NULL);
   assert(npermvars > 0);
   assert(permvardomaincenter != NULL || symtype != SYM_SYMTYPE_SIGNPERM);
   assert(varidxmatrix != NULL);
   assert(nrows > 0);
   assert(ncols > 0);
   assert(orbitopeconss != NULL);
   assert(norbitopeconss != NULL);
   assert(maxnorbitopeconss != NULL);
   assert(orbitopalreddata != NULL);
   assert(lexreddata != NULL);
   assert(lexredactive != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* add linear constraints x_1 >= x_2 >= ... >= x_ncols for single-row orbitopes */
   if( nrows == 1 )
   {
      /* restrict to the packing and partitioning rows */
      SCIP_CONS* cons;
      SCIP_VAR* consvars[2];
      SCIP_Real conscoefs[2] = { -1.0, 1.0 };

      /* for all adjacent column pairs, add linear constraint */
      SCIP_CALL( SCIPensureBlockMemoryArray(scip, orbitopeconss, maxnorbitopeconss, *norbitopeconss + ncols - 1) );
      for( i = 0; i < ncols - 1; ++i )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_1row_comp_%d_col%d", partialname, componentid, i);

         consvars[0] = permvars[varidxmatrix[0][i]];
         consvars[1] = permvars[varidxmatrix[0][i + 1]];

         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, consvars, conscoefs, -SCIPinfinity(scip), 0.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

         SCIP_CALL( SCIPaddCons(scip, cons) );
         (*orbitopeconss)[(*norbitopeconss)++] = cons;
      }

      *success = TRUE;

      return SCIP_OKAY;
   }

   /* for only 2 columns, the the component can be completely handled by lexicographic reduction */
   if( ncols == 2 )
   {
      int* orbisackperm;

      /* If the component is an orbitope with 2 columns, then there is 1 generator of order 2. */
      orbisackperm = perms[0];

      SCIP_CALL( SCIPlexicographicReductionAddPermutation(scip, lexreddata,
            permvars, npermvars, orbisackperm, symtype, permvardomaincenter, TRUE, success) );

      if( *success )
      {
         *lexredactive = TRUE;
         return SCIP_OKAY;
      }
   }

   /* create orbitope variable matrix */
   SCIP_CALL( SCIPallocBufferArray(scip, &varmatrix, nrows) );
   for( i = 0; i < nrows; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &varmatrix[i], ncols) );
      for( j = 0; j < ncols; ++j )
         varmatrix[i][j] = permvars[varidxmatrix[i][j]];
   }

   pprows = NULL;
   SCIP_CALL( SCIPisPackingPartitioningOrbitope(scip, varmatrix, nrows, ncols, &pprows, &npprows, &type) );

   /* does it have at least 3 packing-partitioning rows? */
   ispporbitope = npprows >= 3;  /* (use same magic number as cons_orbitope.c) */

   if( ispporbitope ) /* @todo if it's a pporbitope, we do it statically right now. */
   {
      /* restrict to the packing and partitioning rows */
      SCIP_CONS* cons;
      SCIP_VAR*** ppvarsarrayonlypprows;
      int r;

      assert(pprows != NULL);

      SCIP_CALL( SCIPallocBufferArray(scip, &ppvarsarrayonlypprows, npprows) );

      r = 0;
      for( i = 0; i < nrows; ++i )
      {
         if( pprows[i] )
         {
            assert(r < npprows);
            ppvarsarrayonlypprows[r++] = varmatrix[i];
         }
      }
      assert(r == npprows);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_pp", partialname);
      SCIP_CALL( SCIPcreateConsOrbitope(scip, &cons, name, ppvarsarrayonlypprows, SCIP_ORBITOPETYPE_PACKING,
            npprows, ncols, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(scip, cons) );

      /* check whether we need to resize */
      SCIP_CALL( SCIPensureBlockMemoryArray(scip, orbitopeconss, maxnorbitopeconss, *norbitopeconss + 1) );
      (*orbitopeconss)[(*norbitopeconss)++] = cons;
      *success = TRUE;

      SCIPfreeBufferArray(scip, &ppvarsarrayonlypprows);
   }
   else
   {
      /* use orbitopal reduction for component */
      SCIP_COLUMNORDERING columnordering;
      SCIP_VAR** orbitopevarmatrix;
      int nelem;
      int pos = 0;

      /* variable array */
      nelem = nrows * ncols;
      SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevarmatrix, nelem) );
      for( i = 0; i < nrows; ++i )
      {
         for( j = 0; j < ncols; ++j )
            orbitopevarmatrix[pos++] = varmatrix[i][j];
      }

      /* get column ordering */
      columnordering = SCIPorbitopalReductionGetDefaultColumnOrdering(orbitopalreddata);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_full", partialname);
      SCIP_CALL( SCIPorbitopalReductionAddOrbitope(scip, orbitopalreddata,
            SCIP_ROWORDERING_BRANCHING, columnordering,
            orbitopevarmatrix, nrows, ncols, success) );
      *success = TRUE;

      SCIPfreeBufferArray(scip, &orbitopevarmatrix);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &pprows, nrows);

   for( i = nrows - 1; i >= 0; --i )
   {
      SCIPfreeBufferArray(scip, &varmatrix[i]);
   }
   SCIPfreeBufferArray(scip, &varmatrix);

   return SCIP_OKAY;
}

/** returns whether every variable in each row of an orbitope matrix has same center */
static
SCIP_Bool isEquallyCenteredOrbitope(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_Real*            vardomaincenter,    /**< array containing domain centers of each variable */
   int**                 varidxmatrix,       /**< matrix containing variable indices of orbitope */
   int                   startrow,           /**< row of varidxmatrix in which orbitope begins */
   int                   endrow,             /**< row of varidxmatrix after which orbitope ends */
   int                   startcol,           /**< column of varidxmatrix in which orbitope begins */
   int                   endcol,             /**< column of varidxmatrix after which orbitope ends */
   SCIP_Bool             equalrowcenters     /**< whether rows are centered equally (otherwise, columns) */
   )
{
   int i;
   int j;

   assert(scip != NULL);
   assert(varidxmatrix != NULL);
   assert(0 <= startrow && startrow < endrow);
   assert(0 <= startcol && startcol < endcol);

   /* if we are dealing with unsigned permutations, then all rows have the same center */
   if( vardomaincenter == NULL )
      return TRUE;

   if( equalrowcenters )
   {
      for( i = startrow; i < endrow; ++i )
      {
         for( j = startcol; j < endcol - 1; ++j )
         {
            /* different centers detected */
            if( !SCIPisEQ(scip, vardomaincenter[varidxmatrix[i][j]],
                  vardomaincenter[varidxmatrix[i][j + 1]]) )
               return FALSE;
         }
      }
   }
   else
   {
      for( j = startcol; j < endcol; ++j )
      {
         for( i = startrow; i < endrow - 1; ++i )
         {
            /* different centers detected */
            if( !SCIPisEQ(scip, vardomaincenter[varidxmatrix[i][j]],
                  vardomaincenter[varidxmatrix[i + 1][j]]) )
               return FALSE;
         }
      }
   }

   return TRUE;
}

/** checks whether a proper signed permutation flips a (partial) orbitope column */
static
SCIP_RETCODE hasOrbitopeColumnFlip(
   int**                 varidxmatrix,       /**< matrix containing variable indices of orbitope */
   int                   startrow,           /**< row of varidxmatrix in which orbitope starts */
   int                   endrow,             /**< row of varidxmatrix after which orbitope ends */
   int                   startcol,           /**< column of varidxmatrix in which orbitope starts */
   int                   endcol,             /**< column of varidxmatrix after which orbitope ends */
   int*                  signedperm,         /**< signed permutation to be checked */
   int                   npermvars,          /**< number of variables symmetries act on */
   SCIP_Bool             transposed,         /**< whether the orbitope is transposed in varidxmatrix */
   int*                  flipablerows,       /**< allocated array to store rows admitting a flip */
   int*                  nflipablerows       /**< pointer to store number of flipable rows */
   )
{
   int flipcolumn = -1;
   int i;
   int j;

   assert(varidxmatrix != NULL);
   assert(0 <= startrow && startrow < endrow);
   assert(0 <= startcol && startcol < endcol);
   assert(signedperm != NULL);
   assert(npermvars >= (endrow - startrow) * (endcol - startcol));
   assert(flipablerows != NULL);
   assert(nflipablerows != NULL);

   *nflipablerows = 0;

   /* iterate over matrix and check for flip */
   for( i = startrow; i < endrow; ++i )
   {
      for( j = startcol; j < endcol; ++j )
      {
         if( signedperm[varidxmatrix[i][j]] == npermvars + varidxmatrix[i][j] )
         {
            if( flipcolumn == -1 )
            {
               flipcolumn = transposed ? i : j;
               flipablerows[(*nflipablerows)++] = transposed ? j : i;
            }
            else if( (transposed && flipcolumn != i) || (!transposed && flipcolumn != j) )
            {
               /* variables from two columns are flipped */
               *nflipablerows = 0;
               return SCIP_OKAY;
            }
            else
               flipablerows[(*nflipablerows)++] = transposed ? j : i;
         }
         else if( signedperm[varidxmatrix[i][j]] != varidxmatrix[i][j] )
         {
            /* variable is not flipped */
            *nflipablerows = 0;
            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}

/** handles orbitope action by static or dynamic symmetry handling methods */
static
SCIP_RETCODE handleOrbitope(
   SCIP*                 scip,               /**< SCIP instance */
   SYM_SYMTYPE           symtype,            /**< symmetry type */
   int**                 perms,              /**< permutations */
   SCIP_VAR**            permvars,           /**< variables the permutations act on */
   int                   npermvars,          /**< number of variables */
   SCIP_Real*            permvardomaincenter,/**< domain centers of permvars */
   int                   componentid,        /**< ID of component to which orbitope is added */
   int**                 varidxmatrix,       /**< matrix containing variable indices of orbitope */
   int                   nrows,              /**< number of rows of matrix */
   int                   ncols,              /**< number of columns of matrix */
   char*                 partialname,        /**< partial name to be extended by constraints */
   SCIP_Bool             issigned,           /**< whether the first row of the orbitope can be sign-flipped */
   SCIP_Bool             handlestatically,   /**< whether the orbitope shall be handled statically */
   SCIP_LEXREDDATA*      lexreddata,         /**< data needed for lexicographic reduction */
   SCIP_Bool*            lexredactive,       /**< pointer to store whether lexicographic reduction is active */
   SCIP_ORBITOPALREDDATA* orbitopalreddata,  /**< data for orbitopal reduction */
   SCIP_CONS***          orbitopeconss,      /**< pointer to store array of conss derived from orbitopes */
   int*                  norbitopeconss,     /**< pointer to store number of conss in orbitopeconss */
   int*                  maxnorbitopeconss,  /**< pointer to store maximum number of conss orbitopeconss can hold */
   SCIP_Bool             usedynamicprop,     /**< Shall orbitopes be handled dynamically? */
   SCIP_Bool*            success,            /**< pointer to store whether orbitope could be added successfully */
   int*                  nchgbds             /**< pointer to store number of bound changes (or NULL) */
   )
{
   assert(scip != NULL);
   assert(perms != NULL);
   assert(permvars != NULL);
   assert(npermvars > 0);
   assert(permvardomaincenter != NULL || symtype != SYM_SYMTYPE_SIGNPERM);
   assert(varidxmatrix != NULL);
   assert(nrows > 0);
   assert(ncols > 0);
   assert(lexreddata != NULL);
   assert(lexredactive != NULL);
   assert(orbitopalreddata != NULL);
   assert(orbitopeconss != NULL);
   assert(norbitopeconss != NULL);
   assert(maxnorbitopeconss != NULL);
   assert(success != NULL);
   assert(nchgbds != NULL);

   *success = FALSE;

   if( handlestatically )
   {
      char name[SCIP_MAXSTRLEN];
      SCIP_Real consvals[2] = {-1.0, 1.0};
      SCIP_VAR* consvars[2];
      SCIP_VAR** orbitopevarmatrix;
      SCIP_CONS* cons;
      int nconss;
      int nelem;
      int pos;
      int i;
      int j;

      /* handle orbitope */
      nelem = nrows * ncols;
      SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevarmatrix, nelem) );
      for( i = 0, pos = 0; i < nrows; ++i )
      {
         for( j = 0; j < ncols; ++j )
            orbitopevarmatrix[pos++] = permvars[varidxmatrix[i][j]];
      }

      SCIP_CALL( SCIPorbitopalReductionAddOrbitope(scip, orbitopalreddata, SCIP_ROWORDERING_NONE,
            SCIP_COLUMNORDERING_NONE, orbitopevarmatrix, nrows, ncols, success) );

      /* compute number of constraints to handle signed part of the orbitope */
      nconss = ncols - 1;
      if( issigned )
         nconss += ncols;

      /* create linear constraints */
      SCIP_CALL( SCIPensureBlockMemoryArray(scip, orbitopeconss, maxnorbitopeconss, *norbitopeconss + nconss) );

      /* handle symmetries by enforcing sorted variables in first row */
      for( j = 0; j < ncols - 1; ++j )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_sort_%d", partialname, j);
         consvars[0] = orbitopevarmatrix[j];
         consvars[1] = orbitopevarmatrix[j + 1];

         /* enforce constraints to be in LP since this seems to have a positive impact for orbitopes with cont. vars */
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, consvars, consvals, -SCIPinfinity(scip), 0.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         (*orbitopeconss)[(*norbitopeconss)++] = cons;
         SCIP_CALL( SCIPaddCons(scip, cons) );
      }

      if( issigned )
      {
         SCIP_VAR* var;
         SCIP_Real bound;

         assert(permvardomaincenter != NULL);

         /* the first row is contained in the upper half of the variable domain */
         for( j = 0; j < ncols; ++j )
         {
            var = orbitopevarmatrix[j];
            bound = permvardomaincenter[varidxmatrix[0][j]];

            if( SCIPisLT(scip, SCIPvarGetLbLocal(var), bound) )
            {
               /* improve lower bound */
               SCIP_CALL( SCIPchgVarLb(scip, var, bound) );
               ++(*nchgbds);
            }
         }
      }
      SCIPfreeBufferArray(scip, &orbitopevarmatrix);

      *success = TRUE;
   }

   /* if symmetries have not been handled yet */
   if( !(*success) )
   {
      /* dynamic propagation */
      if( usedynamicprop )
      {
         SCIP_CALL( addOrbitopesDynamic(scip, symtype, perms, permvars, npermvars, permvardomaincenter, componentid,
               partialname, varidxmatrix, nrows, ncols, orbitopeconss, norbitopeconss, maxnorbitopeconss,
               orbitopalreddata, lexreddata, lexredactive, success) );
      }
      /* static variant only for binary variables */
      else
      {
         SCIP_VAR*** orbitopematrix;
         SCIP_CONS* cons;
         int i;
         int j;
         int nbinrows = 0;

         SCIP_CALL( SCIPallocBufferArray(scip, &orbitopematrix, nrows) );
         for( i = 0; i < nrows; ++i )
         {
            /* skip rows without binary variables */
            if( !SCIPvarIsBinary(permvars[varidxmatrix[i][0]]) )
               continue;

            SCIP_CALL( SCIPallocBufferArray(scip, &orbitopematrix[nbinrows], ncols) );
            for( j = 0; j < ncols; ++j )
            {
               assert(SCIPvarIsBinary(permvars[varidxmatrix[i][j]]));
               orbitopematrix[nbinrows][j] = permvars[varidxmatrix[i][j]];
            }
            ++nbinrows;
         }

         if( nbinrows > 0 )
         {
            SCIP_CALL( SCIPcreateConsOrbitope(scip, &cons, partialname, orbitopematrix, SCIP_ORBITOPETYPE_FULL,
                  nbinrows, ncols, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
                  FALSE) );
            SCIP_CALL( SCIPaddCons(scip, cons) );

            /* do not release constraint here - will be done later */
            SCIP_CALL( SCIPensureBlockMemoryArray(scip, orbitopeconss, maxnorbitopeconss, *norbitopeconss + 1) );
            (*orbitopeconss)[(*norbitopeconss)++] = cons;

            *success = TRUE;

            for( i = nbinrows - 1; i >= 0; --i )
            {
               SCIPfreeBufferArray(scip, &orbitopematrix[i]);
            }
         }
         SCIPfreeBufferArray(scip, &orbitopematrix);
      }
   }

   return SCIP_OKAY;
}

/** handles double lex orbitope action by static symmetry handling methods */
static
SCIP_RETCODE handleDoubleLexOrbitope(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_VAR**            permvars,           /**< variables the permutations act on */
   SCIP_Real*            permvardomaincenter,/**< domain center of variables */
   int**                 varidxmatrix,       /**< matrix containing variable indices of orbitope */
   int                   nrows,              /**< number of rows of matrix */
   int                   ncols,              /**< number of columns of matrix */
   char*                 partialname,        /**< partial name to be extended by constraints */
   int                   nsignedrows,        /**< the first number of rows that can be sign-flipped */
   SCIP_ORBITOPALREDDATA* orbitopalreddata,  /**< data for orbitopal reduction */
   SCIP_CONS***          orbitopeconss,      /**< pointer to store array of conss derived from orbitopes */
   int*                  norbitopeconss,     /**< pointer to store number of conss in orbitopeconss */
   int*                  maxnorbitopeconss,  /**< pointer to store maximum number of conss orbitopeconss can hold */
   SCIP_Bool*            success,            /**< pointer to store whether orbitope could be added successfully */
   int*                  nchgbds             /**< pointer to store number of bound changes (or NULL) */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_Real consvals[2] = {-1.0, 1.0};
   SCIP_VAR* consvars[2];
   SCIP_VAR** orbitopevarmatrix;
   SCIP_CONS* cons;
   int nsignedconss;
   int nsortconss;
   int nactiverows;
   int nactrowsprev;
   int nelem;
   int pos;
   int i;
   int j;

   assert(scip != NULL);
   assert(permvars != NULL);
   assert(permvardomaincenter != NULL || nsignedconss == 0);
   assert(varidxmatrix != NULL);
   assert(nrows > 0);
   assert(ncols > 0);
   assert(0 <= nsignedrows && nsignedrows <= nrows);
   assert(orbitopalreddata != NULL);
   assert(orbitopeconss != NULL);
   assert(norbitopeconss != NULL);
   assert(maxnorbitopeconss != NULL);
   assert(success != NULL);
   assert(nchgbds != NULL);

   /* prepare data for orbitope matrices */
   nelem = nrows * ncols;
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevarmatrix, nelem) );

   /* compute number of constraints to handle signed part of the orbitope */
   nsignedconss = 0;
   nsortconss = ncols - 1;
   if( nsignedrows > 0 )
   {
      nactiverows = nsignedrows;
      nactrowsprev = nrows;

      for( j = 0; j < ncols; ++j )
      {
         nsortconss += MAX(nactrowsprev - nactiverows - 1, 0);
         nsignedconss += nactiverows;
         nactrowsprev = nactiverows;

         /* ceil(nactiverows / 2) */
         nactiverows = (int) ((nactiverows + 1) / 2);
      }
      assert(nactiverows >= 1);

      nsortconss += nactiverows - 1;
   }
   else
      nsortconss += nrows - 1;

   /* create linear constraints */
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, orbitopeconss, maxnorbitopeconss,
         *norbitopeconss + nsortconss + nsignedconss) );

   /* sort variables in first row */
   for( j = 0; j < ncols - 1; ++j )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_sort_%d", partialname, j);
      consvars[0] = permvars[varidxmatrix[0][j]];
      consvars[1] = permvars[varidxmatrix[0][j + 1]];

      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, consvars, consvals, -SCIPinfinity(scip), 0.0,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      (*orbitopeconss)[(*norbitopeconss)++] = cons;
      SCIP_CALL( SCIPaddCons(scip, cons) );
   }

   /* handle symmetries by enforcing sorted columns via orbitopal reduction */
   for( i = 0, pos = 0; i < nrows; ++i )
   {
      for( j = 0; j < ncols; ++j )
         orbitopevarmatrix[pos++] = permvars[varidxmatrix[i][j]];
   }

   SCIP_CALL( SCIPorbitopalReductionAddOrbitope(scip, orbitopalreddata,
         SCIP_ROWORDERING_NONE, SCIP_COLUMNORDERING_NONE,
         orbitopevarmatrix, nrows, ncols, success) );

   if( nsignedconss > 0 )
   {
      int k;

      nactiverows = nsignedrows;
      nactrowsprev = nrows;
      for( j = 0; j < ncols; ++j )
      {
         /* ceil(nactiverows / 2) */
         nactiverows = (int) ((nactiverows + 1) / 2);

         /* the second half of active rows can be sorted by linear inequalities */
         for( i = nactiverows; i < nactrowsprev - 1; ++i )
         {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_sort_col_%d_row_%d", partialname, j, i);

            consvars[0] = permvars[varidxmatrix[i][0]];
            consvars[1] = permvars[varidxmatrix[i + 1][0]];

            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, consvars, consvals,
                  -SCIPinfinity(scip), 0.0,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

            (*orbitopeconss)[(*norbitopeconss)++] = cons;
            SCIP_CALL( SCIPaddCons(scip, cons) );
         }

         /* we can also sort the second half of rows by orbitopal reduction */
         if( nactrowsprev - nactiverows > 1 )
         {
            for( k = 0, pos = 0; k < ncols; ++k )
            {
               for( i = nactiverows; i < nactrowsprev; ++i )
                  orbitopevarmatrix[pos++] = permvars[varidxmatrix[i][k]];
            }

            SCIP_CALL( SCIPorbitopalReductionAddOrbitope(scip, orbitopalreddata,
                  SCIP_ROWORDERING_NONE, SCIP_COLUMNORDERING_NONE,
                  orbitopevarmatrix, ncols, nactrowsprev - nactiverows, success) );
         }
         nactrowsprev = nactiverows;

         /* the first half of the active rows are in the upper part of the variable domain */
         for( i = 0; i < nactiverows; ++i )
         {
            SCIP_Real bound;

            consvars[0] = permvars[varidxmatrix[i][j]];
            bound = permvardomaincenter[varidxmatrix[i][j]];

            /* improve lower bound either by changing the bound or a linear constraint */
            if( SCIPisLT(scip, SCIPvarGetLbLocal(consvars[0]), bound) )
            {
               SCIP_CALL( SCIPchgVarLb(scip, consvars[0], bound) );
               ++(*nchgbds);
            }
         }
      }

      /* within the remaining active rows, the rows can be sorted */
      if( nactiverows > 1 )
      {
         for( i = 0; i < nactiverows - 1; ++i )
         {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_sortfirstactive_%d", partialname, i);

            consvars[0] = permvars[varidxmatrix[i][0]];
            consvars[1] = permvars[varidxmatrix[i + 1][0]];

            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, consvars, consvals,
                  -SCIPinfinity(scip), 0.0,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

            (*orbitopeconss)[(*norbitopeconss)++] = cons;
            SCIP_CALL( SCIPaddCons(scip, cons) );
         }

         /* also apply orbitopal reduction */
         for( j = 0, pos = 0; j < ncols; ++j )
         {
            for( i = 0; i < nactiverows; ++i )
               orbitopevarmatrix[pos++] = permvars[varidxmatrix[i][j]];
         }

         SCIP_CALL( SCIPorbitopalReductionAddOrbitope(scip, orbitopalreddata,
               SCIP_ROWORDERING_NONE, SCIP_COLUMNORDERING_NONE,
               orbitopevarmatrix, ncols, nactiverows, success) );
      }
      assert(*norbitopeconss <= *maxnorbitopeconss );
   }
   else
   {
      /* sort first column */
      for( i = 0; i < nrows - 1; ++i )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_sort_col_0_row_%d", partialname, i);

         consvars[0] = permvars[varidxmatrix[i][0]];
         consvars[1] = permvars[varidxmatrix[i + 1][0]];

         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, consvars, consvals,
               -SCIPinfinity(scip), 0.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         (*orbitopeconss)[(*norbitopeconss)++] = cons;
         SCIP_CALL( SCIPaddCons(scip, cons) );
      }

      /* apply orbitopal fixing to row permutations */
      for( j = 0, pos = 0; j < ncols; ++j )
      {
         for( i = 0; i < nrows; ++i )
            orbitopevarmatrix[pos++] = permvars[varidxmatrix[i][j]];
      }
      SCIP_CALL( SCIPorbitopalReductionAddOrbitope(scip, orbitopalreddata,
            SCIP_ROWORDERING_NONE, SCIP_COLUMNORDERING_NONE,
            orbitopevarmatrix, ncols, nrows, success) );
   }

   SCIPfreeBufferArray(scip, &orbitopevarmatrix);

   *success = TRUE;

   return SCIP_OKAY;
}

/** handles double lex matrix */
static
SCIP_RETCODE handleDoublelLexMatrix(
   SCIP*                 scip,               /**< SCIP instance */
   SYM_SYMTYPE           symtype,            /**< symmetry type */
   int**                 perms,              /**< permutations */
   SCIP_VAR**            permvars,           /**< variables the permutations act on */
   int                   npermvars,          /**< number of variables */
   SCIP_Real*            permvardomaincenter,/**< domain center of variables */
   int                   id,                 /**< ID of component that is handled */
   int**                 varidxmatrix,       /**< matrix containing variable indices of double lex matrix */
   int                   nrows,              /**< number of rows of matrix */
   int                   ncols,              /**< number of columns of matrix */
   int*                  rowsbegin,          /**< array indicating where a new row block begins */
   int*                  colsbegin,          /**< array indicating where a new column block begins */
   int                   nrowblocks,         /**< number of row blocks */
   int                   ncolblocks,         /**< number of column blocks */
   int**                 signedperms,        /**< array of proper signed permutations */
   int                   nsignedperms,       /**< number of proper signed permutations */
   SCIP_Bool             handlesignedorbitopes,/**< Shall signed orbitopes be handled? */
   SCIP_Bool             usedynamicprop,     /**< Shall dynamic orbitope propagation be used */
   SCIP_LEXREDDATA*      lexreddata,         /**< data needed for lexicographic reduction */
   SCIP_Bool*            lexredactive,       /**< pointer to store whether lexicographic reduction is active */
   SCIP_ORBITOPALREDDATA* orbitopalreddata,  /**< data for orbitopal reduction */
   SCIP_CONS***          orbitopeconss,      /**< pointer to store array of conss derived from orbitopes */
   int*                  norbitopeconss,     /**< pointer to store number of conss in orbitopeconss */
   int*                  maxnorbitopeconss,  /**< pointer to store maximum number of conss orbitopeconss can hold */
   SCIP_Bool*            success,            /**< pointer to store whether orbitope could be added successfully */
   int*                  nchgbds             /**< pointer to store number of bound changes (or NULL) */
   )
{
   char partialname[SCIP_MAXSTRLEN];
   int** orbitopematrix;
   int* flipableidx;
   int nflipableidx;
   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   SCIP_Bool tmpsuccess;
   int maxdim;
   int i;
   int p;
   int j;

   assert(scip != NULL);
   assert(perms != NULL);
   assert(permvars != NULL);
   assert(permvardomaincenter != NULL || symtype != SYM_SYMTYPE_SIGNPERM);
   assert(varidxmatrix != NULL);
   assert(nrows > 0);
   assert(ncols > 0);
   assert(rowsbegin != NULL);
   assert(colsbegin != NULL);
   assert(nrowblocks > 0);
   assert(ncolblocks > 0);
   assert(signedperms != NULL);
   assert(nsignedperms >= 0);
   assert(lexreddata != NULL);
   assert(lexredactive != NULL);
   assert(orbitopalreddata != NULL);
   assert(orbitopeconss != NULL);
   assert(norbitopeconss != NULL);
   assert(maxnorbitopeconss != NULL);
   assert(success != NULL);
   assert(nchgbds != NULL);

   *success = FALSE;

   maxdim = MAX(nrows, ncols);
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitopematrix, maxdim) );
   for( i = 0; i < maxdim; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &orbitopematrix[i], maxdim) );
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &flipableidx, maxdim) );

   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, 2) );

   /* detect how double-lex matrices can be handled
    *
    * We use the following strategy:
    *
    * (1) In case of 1 column- and row-block, we check whether one of them can be handled
    *     by signed orbitopes. If this is the case, reorder the rows or columns based
    *     on the signed permutations.
    *
    * (2) If there are multiple column- and row-blocks, just handle classical permutation symmetries.
    */

   if( nrowblocks == 1 && ncolblocks == 1 )
   {
      /* Case 1 */
      SCIP_Bool hascolflip = FALSE;
      SCIP_Bool hasrowflip = FALSE;
      SCIP_Bool canusecolorbitope = FALSE;
      SCIP_Bool canuseroworbitope = FALSE;
      int q;

      /* check whether orbitopes can be used to handle column and row swaps (requires equally centered rows/columns) */
      if( isEquallyCenteredOrbitope(scip, permvardomaincenter, varidxmatrix, 0, nrows, 0, ncols, TRUE) )
         canusecolorbitope = TRUE;
      if( isEquallyCenteredOrbitope(scip, permvardomaincenter, varidxmatrix, 0, nrows, 0, ncols, FALSE) )
         canuseroworbitope = TRUE;

      nflipableidx = 0;
      if( handlesignedorbitopes && canusecolorbitope )
      {
         /* check whether the signed permutations flip entries within a single column of the orbitope matrix */
         for( q = 0; q < nsignedperms && nflipableidx == 0; ++q )
         {
            SCIP_CALL( hasOrbitopeColumnFlip(varidxmatrix, 0, nrows, 0, ncols,
                  signedperms[q], npermvars, FALSE, flipableidx, &nflipableidx) );
         }

         /* a signed orbitope can be used to handle the column symmetries */
         if( nflipableidx > 0 )
            hascolflip = TRUE;
      }

      if( handlesignedorbitopes && !hascolflip && canuseroworbitope )
      {
         assert(nflipableidx == 0);

         /* check whether the signed permutations flip entries within a single column of transposed orbitope matrix */
         for( q = 0; q < nsignedperms && nflipableidx == 0; ++q )
         {
            SCIP_CALL( hasOrbitopeColumnFlip(varidxmatrix, 0, nrows, 0, ncols,
                  signedperms[q], npermvars, TRUE, flipableidx, &nflipableidx) );
         }

         /* a signed orbitope can be used to handle the row symmetries */
         if( nflipableidx > 0 )
            hasrowflip = TRUE;
      }

      /* handle row and column symmetries by potentially reordered orbitopes */
      if( hascolflip )
      {
         int isigned = 0;
         int iunsigned;

         /* handle column symmetries by reordering the rows of the orbitope */
         iunsigned = nflipableidx;
         for( i = 0; i < nrows; ++i )
         {
            if( isigned < nflipableidx && flipableidx[isigned] == i )
            {
               for( j = 0; j < ncols; ++j )
                  orbitopematrix[isigned][j] = varidxmatrix[i][j];
               ++isigned;
            }
            else
            {
               for( j = 0; j < ncols; ++j )
                  orbitopematrix[iunsigned][j] = varidxmatrix[i][j];
               ++iunsigned;
            }
         }
         assert(isigned == nflipableidx);
         assert(iunsigned == nrows);

         (void) SCIPsnprintf(partialname, SCIP_MAXSTRLEN, "orbitope_component_%d_doublelex_col_0", id);

         SCIP_CALL( handleDoubleLexOrbitope(scip, permvars, permvardomaincenter, orbitopematrix, nrows, ncols,
               partialname, nflipableidx, orbitopalreddata, orbitopeconss, norbitopeconss, maxnorbitopeconss,
               &tmpsuccess, nchgbds) );
         *success = *success || tmpsuccess;
      }
      else if( hasrowflip )
      {
         int jsigned = 0;
         int junsigned;

         /* handle row symmetries by reordering the columns of the orbitope */
         junsigned = nflipableidx;
         for( j = 0; j < ncols; ++j )
         {
            if( jsigned < nflipableidx && flipableidx[jsigned] == j )
            {
               for( i = 0; i < nrows; ++i )
                  orbitopematrix[jsigned][i] = varidxmatrix[i][j];
               ++jsigned;
            }
            else
            {
               for( i = 0; i < nrows; ++i )
                  orbitopematrix[junsigned][i] = varidxmatrix[i][j];
               ++junsigned;
            }
         }
         assert(jsigned == nflipableidx);
         assert(junsigned == ncols);

         (void) SCIPsnprintf(partialname, SCIP_MAXSTRLEN, "orbitope_component_%d_doublelex_row_0", id);

         SCIP_CALL( handleDoubleLexOrbitope(scip, permvars, permvardomaincenter, orbitopematrix, ncols, nrows,
               partialname, nflipableidx, orbitopalreddata, orbitopeconss, norbitopeconss, maxnorbitopeconss,
               &tmpsuccess, nchgbds) );
         *success = *success || tmpsuccess;
      }
   }

   /* if no symmetries have been handled yet, handle column and row symmetries without signed permutations */
   if( !(*success) )
   {
      /* ensure that we can store orbitope constraints in probdata */
      SCIP_CALL( SCIPensureBlockMemoryArray(scip, orbitopeconss, maxnorbitopeconss,
            *norbitopeconss + nrowblocks + ncolblocks) );

      /* handle column symmetries via original column and row ordering */
      for( p = 0; p < ncolblocks; ++p )
      {
         int jj;
         j = 0;

         /* we can only handle the orbitope if all variables in a row have the same domain center */
         if( !isEquallyCenteredOrbitope(scip, permvardomaincenter, varidxmatrix, 0, nrows,
               colsbegin[p], colsbegin[p + 1], TRUE) )
            continue;

         /* create the orbitope matrix */
         for( i = 0; i < nrows; ++i )
         {
            for( j = 0, jj = colsbegin[p]; jj < colsbegin[p + 1]; ++j, ++jj )
               orbitopematrix[i][j] = varidxmatrix[i][jj];
         }

         (void) SCIPsnprintf(partialname, SCIP_MAXSTRLEN, "orbitope_component_%d_doublelex_col_%d", id, p);

         SCIP_CALL( handleOrbitope(scip, symtype, perms, permvars, npermvars, permvardomaincenter, id, orbitopematrix,
               nrows, j, partialname, FALSE, TRUE, lexreddata, lexredactive, orbitopalreddata,
               orbitopeconss, norbitopeconss, maxnorbitopeconss, usedynamicprop, &tmpsuccess, nchgbds) );
         *success = *success || tmpsuccess;
      }

      /* handle row symmetries via original column and row ordering */
      for( p = 0; p < nrowblocks; ++p )
      {
         int ii;

         /* we can only handle the orbitope if all variables in a row have the same domain center */
         if( !isEquallyCenteredOrbitope(scip, permvardomaincenter, varidxmatrix,
               rowsbegin[p], rowsbegin[p + 1], 0, ncols, FALSE) )
            continue;

         /* create the orbitope matrix */
         for( i = 0, ii = rowsbegin[p]; ii < rowsbegin[p + 1]; ++i, ++ii )
         {
            for( j = 0; j < ncols; ++j )
               orbitopematrix[j][i] = varidxmatrix[ii][j];
         }

         (void) SCIPsnprintf(partialname, SCIP_MAXSTRLEN, "orbitope_component_%d_doublelex_row_%d", id, p);

         SCIP_CALL( handleOrbitope(scip, symtype, perms, permvars, npermvars, permvardomaincenter, id, orbitopematrix,
               ncols, i, partialname, FALSE, TRUE, lexreddata, lexredactive, orbitopalreddata, orbitopeconss,
               norbitopeconss, maxnorbitopeconss, usedynamicprop, &tmpsuccess, nchgbds) );
         *success = *success || tmpsuccess;
      }
   }

   SCIPfreeBufferArray(scip, &consvals);
   SCIPfreeBufferArray(scip, &consvars);
   SCIPfreeBufferArray(scip, &flipableidx);

   for( i = maxdim - 1; i >= 0; --i )
   {
      SCIPfreeBufferArray(scip, &orbitopematrix[i]);
   }
   SCIPfreeBufferArray(scip, &orbitopematrix);

   return SCIP_OKAY;
}

/** tries to handle symmetries of single lex matrices (orbitopes) or double lex matrices */
static
SCIP_RETCODE tryHandleSingleOrDoubleLexMatrices(
   SCIP*                 scip,               /**< SCIP instance */
   SYM_SYMTYPE           symtype,            /**< symmetry type */
   int**                 perms,              /**< permutations */
   int                   nperms,             /**< number of permutations */
   SCIP_VAR**            permvars,           /**< variables the permutations act on */
   int                   npermvars,          /**< number of variables */
   SCIP_Real*            permvardomaincenter,/**< array of centers of variable domains */
   SCIP_Bool             detectsinglelex,    /**< whether single lex matrices shall be detected */
   SCIP_Bool             handlesignedorbitopes,/**< Shall signed orbitopes be handled? */
   SCIP_Bool             usedynamicprop,     /**< Shall dynamic propagation be used for orbitopes? */
   int                   cidx,               /**< index of component */
   int*                  nchgbds,            /**< pointer to store number of bound changes (or NULL) */
   SCIP_LEXREDDATA*      lexreddata,         /**< data needed for lexicographic reduction */
   SCIP_Bool*            lexredactive,       /**< pointer to store whether lexicographic reduction is active */
   SCIP_ORBITOPALREDDATA* orbitopalreddata,  /**< data for orbitopal reduction */
   SCIP_CONS***          orbitopeconss,      /**< pointer to store array of conss derived from orbitopes */
   int*                  norbitopeconss,     /**< pointer to store number of conss in orbitopeconss */
   int*                  maxnorbitopeconss,  /**< pointer to store maximum number of conss orbitopeconss can hold */
   SCIP_Bool*            success             /**< pointer to store whether symmetries are handled */
   )
{
   int** lexmatrix = NULL;
   int* lexrowsbegin = NULL;
   int* lexcolsbegin = NULL;
   int nrows;
   int ncols;
   int nrowmatrices;
   int ncolmatrices;
   int** selectedperms = NULL;
   int i;
   int nonpermidx = -1;
   SCIP_Real percentageunsigned;
   SCIP_Bool isorbitope;
   SCIP_Bool locsuccess = FALSE;
   SCIP_Bool* isproperperm;
   int nselectedperms = 0;

   assert(scip != NULL);
   assert(perms != NULL);
   assert(nperms > 0);
   assert(permvars != NULL);
   assert(npermvars > 0);
   assert(nchgbds != NULL);
   assert(lexreddata != NULL);
   assert(lexredactive != NULL);
   assert(orbitopalreddata != NULL);
   assert(orbitopeconss != NULL);
   assert(norbitopeconss != NULL);
   assert(maxnorbitopeconss != NULL);
   assert(success != NULL);

   *nchgbds = 0;
   *success = FALSE;

   /* get proper permutations of component and possibly store index of one non-permutation */
   SCIP_CALL( SCIPallocBufferArray(scip, &isproperperm, nperms) );
   if( symtype == SYM_SYMTYPE_PERM )
   {
      selectedperms = perms;
      nselectedperms = nperms;
      for( i = 0; i < nperms; ++i )
         isproperperm[i] = TRUE;
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &selectedperms, nperms) );
      for( i = 0; i < nperms; ++i )
      {
         isproperperm[i] = isProperPerm(symtype, perms[i], npermvars);
         if( isproperperm[i] )
            selectedperms[nselectedperms++] = perms[i];
         else
            nonpermidx = i;
      }
   }
   percentageunsigned = (SCIP_Real) nselectedperms / (SCIP_Real) nperms;

   if( nselectedperms == 0 )
   {
      SCIPfreeBufferArray(scip, &isproperperm);
      if( symtype != SYM_SYMTYPE_PERM )
      {
         SCIPfreeBufferArrayNull(scip, &selectedperms);
      }
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPdetectSingleOrDoubleLexMatrices(scip, detectsinglelex, selectedperms, nselectedperms, npermvars,
         &locsuccess, &isorbitope, &lexmatrix, &nrows, &ncols, &lexrowsbegin, &lexcolsbegin,
         &nrowmatrices, &ncolmatrices) );

   /* possibly handle double lex matrix or orbitope */
   if( locsuccess )
   {
      assert(lexmatrix != NULL);
      assert(nrows > 0);
      assert(ncols > 0);

      /* possibly store non-permutation symmetries (in many cases, there is only one) */
      if( nselectedperms != nperms )
      {
         assert(nselectedperms < nperms);

         if( nselectedperms == nperms - 1 )
         {
            selectedperms[0] = selectedperms[nonpermidx];
            nselectedperms = 1;
         }
         else
         {
            nselectedperms = 0;
            for( i = 0; i < nperms; ++i )
            {
               if( !isproperperm[i] )
                  selectedperms[nselectedperms++] = perms[i];
            }
         }
      }

      if( isorbitope )
      {
         int** orbitopematrix;
         char partialname[SCIP_MAXSTRLEN];

         locsuccess = FALSE;

         /* signed permutations can only handle the orbitope if all variables per row have the same domain center */
         if( symtype != (int) SYM_SYMTYPE_PERM )
         {
            if( !isEquallyCenteredOrbitope(scip, permvardomaincenter, lexmatrix, 0, nrows, 0, ncols, TRUE) )
               goto FREEMEMORY;
         }

         (void) SCIPsnprintf(partialname, SCIP_MAXSTRLEN, "orbitope_component_%d", cidx);

         if( handlesignedorbitopes )
         {
            int* flipablerows;
            int nflipablerows = 0;
            int p;

            SCIP_CALL( SCIPallocBufferArray(scip, &flipablerows, nrows) );

            /* check whether the signed permutations flip entries within a single column of the orbitope matrix
             *
             * It is sufficient to find one such signed permutations, because only one row will incorporate
             * information about the sign change.
             */
            for( p = 0; p < nselectedperms && nflipablerows == 0; ++p )
            {
               SCIP_CALL( hasOrbitopeColumnFlip(lexmatrix, 0, nrows, 0, ncols, selectedperms[p], npermvars, FALSE,
                     flipablerows, &nflipablerows) );
            }

            /* possibly flip rows to be able to handle signed orbitopes */
            if( nflipablerows > 0 )
            {
               int j;
               int isigned = 0;
               int iunsigned;

               iunsigned = nflipablerows;
               SCIP_CALL( SCIPallocBufferArray(scip, &orbitopematrix, nrows) );
               for( i = 0; i < nrows; ++i )
               {
                  SCIP_CALL( SCIPallocBufferArray(scip, &orbitopematrix[i], ncols) );

                  if( isigned < nflipablerows && flipablerows[isigned] == i )
                  {
                     for( j = 0; j < ncols; ++j )
                        orbitopematrix[isigned][j] = lexmatrix[i][j];
                     ++isigned;
                  }
                  else
                  {
                     for( j = 0; j < ncols; ++j )
                        orbitopematrix[iunsigned][j] = lexmatrix[i][j];
                     ++iunsigned;
                  }
               }
               assert(isigned == nflipablerows);
               assert(iunsigned == nrows);

               SCIP_CALL( handleOrbitope(scip, symtype, perms, permvars, npermvars, permvardomaincenter, cidx,
                     orbitopematrix, nrows, ncols, partialname, TRUE, TRUE, lexreddata, lexredactive, orbitopalreddata,
                     orbitopeconss, norbitopeconss, maxnorbitopeconss, usedynamicprop, &locsuccess, nchgbds) );

               for( i = nrows - 1; i >= 0; --i )
               {
                  SCIPfreeBufferArray(scip, &orbitopematrix[i]);
               }
               SCIPfreeBufferArray(scip, &orbitopematrix);
            }
            SCIPfreeBufferArray(scip, &flipablerows);
         }

         /* if we have not handled the orbitope yet, handle it as unsigned orbitope and the orbitope is large */
         if( (!locsuccess) && percentageunsigned > 0.8 )
         {
            SCIP_CALL( handleOrbitope(scip, symtype, perms, permvars, npermvars, permvardomaincenter, cidx, lexmatrix,
                  nrows, ncols, partialname, FALSE, FALSE, lexreddata, lexredactive, orbitopalreddata, orbitopeconss,
                  norbitopeconss, maxnorbitopeconss, usedynamicprop, &locsuccess, nchgbds) );
         }
      }
      else
      {
         SCIP_CALL( handleDoublelLexMatrix(scip, symtype, perms, permvars, npermvars, permvardomaincenter, cidx,
               lexmatrix, nrows, ncols, lexrowsbegin, lexcolsbegin, nrowmatrices, ncolmatrices, selectedperms, nselectedperms,
               handlesignedorbitopes, usedynamicprop, lexreddata, lexredactive, orbitopalreddata, orbitopeconss,
               norbitopeconss, maxnorbitopeconss, &locsuccess, nchgbds) );
      }

   FREEMEMORY:
      /* free memory not needed anymore */
      for( i = nrows - 1; i >= 0; --i )
      {
         SCIPfreeBlockMemoryArray(scip, &lexmatrix[i], ncols);
      }
      SCIPfreeBlockMemoryArray(scip, &lexmatrix, nrows);

      if( ncolmatrices > 0 )
      {
         SCIPfreeBlockMemoryArray(scip, &lexcolsbegin, ncolmatrices + 1);
      }

      if( nrowmatrices > 0 )
      {
         SCIPfreeBlockMemoryArray(scip, &lexrowsbegin, nrowmatrices + 1);
      }
   }
   SCIPfreeBufferArray(scip, &isproperperm);
   if( symtype != SYM_SYMTYPE_PERM )
   {
      SCIPfreeBufferArrayNull(scip, &selectedperms);
   }

   if( locsuccess )
      *success = TRUE;

   return SCIP_OKAY;
}

/*
 * Callback methods of symmetry handler
 */

/** addition method for symmetry method handler plugins (tries to add symmetry handling method for given symmetries) */
static
SCIP_DECL_SYMHDLRTRYADD(symhdlrTryaddRowCol)
{  /*lint --e{715}*/
   SCIP_SYMHDLRDATA* symhdlrdata;

   assert(success != NULL);
   assert(symhdlr != NULL);
   assert(scip != NULL);
   assert(perms != NULL);
   assert(nperms >= 0);
   assert(permvars != NULL || npermvars == 0);
   assert(naddedconss != NULL);
   assert(nchgbds != NULL);

   *success = FALSE;
   *nchgbds = 0;

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata != NULL);
   assert(symhdlrdata->shadowtreeeventhdlr != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, symcompdata) );
   (*symcompdata)->conss = NULL;
   (*symcompdata)->nconss = 0;
   (*symcompdata)->maxnconss = 0;
   (*symcompdata)->lexreddata = NULL;
   (*symcompdata)->lexredactive = FALSE;

   SCIP_CALL( SCIPincludeLexicographicReduction(scip, &(*symcompdata)->lexreddata, symhdlrdata->shadowtreeeventhdlr) );
   assert((*symcompdata)->lexreddata != NULL);

   SCIP_CALL( tryHandleSingleOrDoubleLexMatrices(scip, symtype, perms, nperms, permvars, npermvars, permvardomcenter,
         !symhdlrdata->detectdoublelex, symhdlrdata->handlesignedorbitopes, symhdlrdata->usedynamicprop,
         id, nchgbds, (*symcompdata)->lexreddata, &(*symcompdata)->lexredactive, symhdlrdata->orbitopalreddata,
         &(*symcompdata)->conss, &(*symcompdata)->nconss, &(*symcompdata)->maxnconss, success) );

   if( !(*success) )
   {
      assert(*nchgbds == 0);

      SCIP_CALL( SCIPlexicographicReductionReset(scip, (*symcompdata)->lexreddata) );
      SCIP_CALL( SCIPlexicographicReductionFree(scip, &(*symcompdata)->lexreddata) );

      SCIPfreeBlockMemory(scip, symcompdata);
      *symcompdata = NULL;

      return SCIP_OKAY;
   }

   *naddedconss = (*symcompdata)->nconss;

   return SCIP_OKAY;
}

/** solving process deinitialization method of symmetry handler (called before branch and bound process data is freed) */
static
SCIP_DECL_SYMHDLREXITSOL(symhdlrExitsolRowCol)
{ /*lint --e{715}*/
   SCIP_SYMHDLRDATA* symhdlrdata;
   SCIP_SYMCOMPDATA* symdata;
   int s;

   assert(symcomps != NULL || nsymcomps == 0);
   assert(symhdlr != NULL);

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata != NULL);

   for( s = 0; s < nsymcomps; ++s )
   {
      assert(symcomps[s] != NULL);

      symdata = SCIPsymcompGetData(symcomps[s]);
      assert(symdata != NULL);
      assert(symdata->lexreddata != NULL || !symdata->lexredactive);

      SCIP_CALL( SCIPlexicographicReductionReset(scip, symdata->lexreddata) );
   }
   assert(symhdlrdata->orbitopalreddata != NULL);

   SCIP_CALL( SCIPorbitopalReductionReset(scip, symhdlrdata->orbitopalreddata) );

   return SCIP_OKAY;
}

/** deinitialization method of symmetry handler (called before transformed problem is freed) */
static
SCIP_DECL_SYMHDLREXIT(symhdlrExitRowCol)
{  /*lint --e{715}*/
   SCIP_SYMHDLRDATA* symhdlrdata;
   SCIP_SYMCOMPDATA* symdata;
   int s;
   int c;

   assert(symcomps != NULL || nsymcomps == 0);
   assert(symhdlr != NULL);

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata != NULL);

   for( s = 0; s < nsymcomps; ++s )
   {
      assert(symcomps[s] != NULL);

      symdata = SCIPsymcompGetData(symcomps[s]);
      assert(symdata != NULL);
      assert(symdata->lexreddata != NULL || !symdata->lexredactive);

      SCIP_CALL( SCIPlexicographicReductionReset(scip, symdata->lexreddata) );
      SCIP_CALL( SCIPlexicographicReductionFree(scip, &symdata->lexreddata) );

      for( c = 0; c < symdata->nconss; ++c )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &symdata->conss[c]) );
      }
      SCIPfreeBlockMemoryArrayNull(scip, &symdata->conss, symdata->maxnconss);

      SCIPfreeBlockMemory(scip, &symdata);
   }
   assert(symhdlrdata->orbitopalreddata != NULL);

   SCIP_CALL( SCIPorbitopalReductionReset(scip, symhdlrdata->orbitopalreddata) );

   return SCIP_OKAY;
}

/** destructor of symmetry handler to free symmetry handler data (called when SCIP is exiting) */
static
SCIP_DECL_SYMHDLRFREE(symhdlrFreeRowCol)
{  /*lint --e{715}*/
   SCIP_SYMHDLRDATA* symhdlrdata;

   assert(scip != NULL);
   assert(symhdlr != NULL);

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata != NULL);
   assert(symhdlrdata->orbitopalreddata != NULL);

   SCIP_CALL( SCIPorbitopalReductionFree(scip, &symhdlrdata->orbitopalreddata) );

   SCIPfreeBlockMemory(scip, &symhdlrdata);

   return SCIP_OKAY;
}

/** domain propagation method of symmetry handler */
static
SCIP_DECL_SYMHDLRPROP(symhdlrPropRowCol)
{  /*lint --e{715}*/
   SCIP_SYMHDLRDATA* symhdlrdata;
   SCIP_SYMCOMPDATA* symcompdata;
   SCIP_Bool infeasible = FALSE;
   SCIP_Bool didrun = FALSE;
   SCIP_Bool didrunlocal;
   int nredlocal;
   int nreds = 0;
   int s;

   assert(symhdlr != NULL);
   assert(symcomps != NULL || nsymcomps == 0);
   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   /* do not run if we are in the root or not yet solving */
   if( SCIPgetDepth(scip) <= 0 || SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata != NULL);

   /* run orbital reduction */
   SCIP_CALL( SCIPorbitopalReductionPropagate(scip, symhdlrdata->orbitopalreddata,
         &infeasible, &nreds, &didrun) );

   if( infeasible )
      return SCIP_OKAY;

   /* run lexicographic reduction */
   for( s = 0; s < nsymcomps; ++s )
   {
      assert(symcomps[s] != NULL);

      symcompdata = SCIPsymcompGetData(symcomps[s]);
      assert(symcompdata != NULL);

      if( !symcompdata->lexredactive )
         continue;

     SCIP_CALL( SCIPlexicographicReductionPropagate(scip, symcompdata->lexreddata,
            &infeasible, &nredlocal, &didrunlocal) );
      nreds += nredlocal;
      didrun = didrun || didrunlocal;
      if( infeasible )
         return SCIP_OKAY;
   }

   if( infeasible )
      *result = SCIP_CUTOFF;
   else if( nreds > 0 )
      *result = SCIP_REDUCEDDOM;
   else if( didrun )
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** presolving method of symmetry handler */
static
SCIP_DECL_SYMHDLRPRESOL(symhdlrPresolRowCol)
{  /*lint --e{715}*/
   SCIP_SYMCOMPDATA* symdata;
   int s;
   int c;

   assert(result != NULL);
   assert(symcomps != NULL || nsymcomps == 0);

   *result = nsymcomps > 0 ? SCIP_DIDNOTFIND : SCIP_DIDNOTRUN;

   for( s = 0; s < nsymcomps; ++s )
   {
      symdata = SCIPsymcompGetData(symcomps[s]);

      for( c = 0; c < symdata->nconss; ++c )
      {
         SCIP_CALL( SCIPpresolCons(scip, symdata->conss[c], nrounds, presoltiming, nnewfixedvars, nnewaggrvars,
               nnewchgvartypes, nnewchgbds, nnewholes, nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs,
               nnewchgsides, nfixedvars, naggrvars, nchgvartypes, nchgbds, naddholes, ndelconss, naddconss,
               nupgdconss, nchgcoefs, nchgsides, result) );

         /* exit if cutoff or unboundedness has been detected */
         if( *result == SCIP_CUTOFF || *result == SCIP_UNBOUNDED )
         {
            SCIPdebugMsg(scip, "Presolving constraint <%s> detected cutoff or unboundedness.\n",
               SCIPconsGetName(symdata->conss[c]));
            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}

/** symmetry component display method of symmetry handler */
static
SCIP_DECL_SYMHDLRPRINT(symhdlrPrintRowCol)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}

/** include symmetry handler for row and column symmetries */
SCIP_RETCODE SCIPincludeSymhdlrRowCol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SYMHDLRDATA* symhdlrdata = NULL;
   SCIP_EVENTHDLR* eventhdlr;

   assert(scip != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &symhdlrdata) );

   SCIP_CALL( SCIPincludeSymhdlrBasic(scip, SYM_NAME, SYM_DESC, SYM_PRIORITY, SYM_PROPPRIORITY, 0, -1,
         SYM_PROPFREQ, -1, SYM_DELAYPROP, FALSE, 1.0, 1, SYM_PROPTIMING, SCIP_PRESOLTIMING_FAST,
         symhdlrTryaddRowCol, NULL, symhdlrFreeRowCol, NULL, symhdlrExitRowCol,
         NULL, symhdlrExitsolRowCol, NULL, NULL, NULL, symhdlrPropRowCol,
         NULL, symhdlrPresolRowCol, symhdlrPrintRowCol, symhdlrdata) );

   /* include shadow tree event handler if it is not included yet */
   eventhdlr = SCIPfindEventhdlr(scip, "event_shadowtree");
   if( eventhdlr == NULL )
   {
      SCIP_CALL( SCIPincludeEventHdlrShadowTree(scip, &symhdlrdata->shadowtreeeventhdlr) );
      assert(symhdlrdata->shadowtreeeventhdlr != NULL);
   }
   else
      symhdlrdata->shadowtreeeventhdlr = eventhdlr;

   SCIP_CALL( SCIPincludeOrbitopalReduction(scip, &symhdlrdata->orbitopalreddata) );
   assert(symhdlrdata->orbitopalreddata != NULL);

   /* add parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "symmetries/" SYM_NAME "/detectdoublelex",
         "Should we check whether the components can be handled by double lex matrices?",
         &symhdlrdata->detectdoublelex, TRUE, SYM_DEFAULT_DETECTDOUBLELEX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "symmetries/" SYM_NAME "/handlesignedorbitopes",
         "Shall signed orbitopes be handled?",
         &symhdlrdata->handlesignedorbitopes, TRUE, SYM_DEFAULT_HANDLESIGNEDORBITOPES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "symmetries/" SYM_NAME "/usedynamicprop",
         "Shall dynamic orbitope propagation be used?",
         &symhdlrdata->usedynamicprop, TRUE, SYM_DEFAULT_USEDYNAMICPROP, NULL, NULL) );

   return SCIP_OKAY;
}

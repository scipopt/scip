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

/**@file   sol.c
 * @brief  datastructures and methods for storing primal IP solutions
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <limits.h>

#include "sol.h"



/*
 * dynamic memory arrays
 */

static
RETCODE solEnsureValsMem(               /**< resizes vals array to be able to store the given index */
   SOL*             sol,                /**< primal IP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var                 /**< variable to get storage for */
   )
{
   int index;

   assert(sol != NULL);
   assert(sol->vars != NULL || sol->firstindex == -1);
   assert(sol->vals != NULL || sol->firstindex == -1);
   assert(set != NULL);
   assert(var != NULL);

   index = var->index;

   if( index >= sol->firstindex && index < sol->firstindex + sol->nvals )
   {
      int pos;

      pos = index - sol->firstindex;
      assert(0 <= pos && pos < sol->nvals);
      if( sol->vars[pos] == NULL )
      {
         assert(sol->vals[pos] == 0.0);
         sol->vars[pos] = var;
         SCIPvarCapture(var);
      }
   }
   else if( sol->firstindex == -1 )
   {
      if( sol->vars == NULL )
      {
         int newsize;
         
         assert(sol->vals == NULL);
         assert(sol->valssize == 0);
         assert(SCIPsetIsZero(set, sol->obj));
         newsize = SCIPsetCalcMemGrowSize(set, 1);
         ALLOC_OKAY( allocBlockMemoryArray(memhdr, sol->vars, newsize) );
         ALLOC_OKAY( allocBlockMemoryArray(memhdr, sol->vals, newsize) );
         sol->valssize = newsize;
      }
      assert(sol->valssize > 0);
      sol->nvals = 1;
      sol->firstindex = index;      
      sol->vars[0] = var;
      sol->vals[0] = 0.0;
      SCIPvarCapture(var);
   }
   else if( index < sol->firstindex )
   {
      VAR** newvars;
      Real* newvals;
      int indexinc;
      int newnvals;
      int newsize;
      int i;

      indexinc = sol->firstindex - index;
      newnvals = sol->nvals + indexinc;
      newsize = SCIPsetCalcMemGrowSize(set, newnvals);
      ALLOC_OKAY( allocBlockMemoryArray(memhdr, newvars, newsize) );
      ALLOC_OKAY( allocBlockMemoryArray(memhdr, newvals, newsize) );
      copyMemoryArray(&newvars[indexinc], sol->vars, sol->nvals);
      copyMemoryArray(&newvals[indexinc], sol->vals, sol->nvals);
      newvars[0] = var;
      newvals[0] = 0.0;
      for( i = 1; i < indexinc; ++i )
      {
         newvars[i] = NULL;
         newvals[i] = 0.0;
      }
      freeBlockMemoryArray(memhdr, sol->vars, sol->valssize);
      freeBlockMemoryArray(memhdr, sol->vals, sol->valssize);
      sol->vars = newvars;
      sol->vals = newvals;
      sol->valssize = newsize;
      sol->nvals = newnvals;
      sol->firstindex = index;
      SCIPvarCapture(var);
   }
   else
   {
      int newnvals;
      int i;

      assert(index >= sol->firstindex + sol->nvals);
      newnvals = index - sol->firstindex + 1;
      if( newnvals > sol->valssize )
      {
         int newsize;

         newsize = SCIPsetCalcMemGrowSize(set, newnvals);
         ALLOC_OKAY( reallocBlockMemoryArray(memhdr, sol->vars, sol->valssize, newsize) );
         ALLOC_OKAY( reallocBlockMemoryArray(memhdr, sol->vals, sol->valssize, newsize) );
         sol->valssize = newsize;
      }
      for( i = sol->nvals; i < newnvals-1; ++i )
      {
         sol->vars[i] = NULL;
         sol->vals[i] = 0.0;
      }
      sol->vars[newnvals-1] = var;
      sol->vals[newnvals-1] = 0.0;
      sol->nvals = newnvals;
      SCIPvarCapture(var);
   }
   assert(sol->firstindex <= index && index < sol->firstindex + sol->nvals);
   assert(sol->nvals <= sol->valssize);

   return SCIP_OKAY;
}



RETCODE SCIPsolCreate(                  /**< creates primal IP solution */
   SOL**            sol,                /**< pointer to primal IP solution */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(sol != NULL);
   assert(memhdr != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, *sol) );   
   (*sol)->vars = NULL;
   (*sol)->vals = NULL;
   (*sol)->obj = 0.0;
   (*sol)->nvals = 0;
   (*sol)->valssize = 0;
   (*sol)->firstindex = -1;

   return SCIP_OKAY;
}

void SCIPsolFree(                       /**< frees primal IP solution */
   SOL**            sol,                /**< pointer to primal IP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data (or NULL, if it's an original variable) */
   )
{
   assert(sol != NULL);
   assert(*sol != NULL);

   SCIPsolClear(*sol, memhdr, set, lp);

   freeBlockMemoryArrayNull(memhdr, (*sol)->vars, (*sol)->valssize);
   freeBlockMemoryArrayNull(memhdr, (*sol)->vals, (*sol)->valssize);
   freeBlockMemory(memhdr, *sol);
}

void SCIPsolClear(                      /**< clears primal IP solution */
   SOL*             sol,                /**< primal IP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data (or NULL, if it's an original variable) */
   )
{
   int v;

   assert(sol != NULL);

   /* release variables in solution */
   for( v = 0; v < sol->nvals; ++v )
      SCIPvarRelease(&sol->vars[v], memhdr, set, lp);

   sol->obj = 0.0;
   sol->nvals = 0;
   sol->firstindex = -1;
}

RETCODE SCIPsolSetVal(                  /**< sets value of variable in primal IP solution */
   SOL*             sol,                /**< primal IP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< variable to add to solution */
   Real             val                 /**< solution value of variable */
   )
{
   int pos;

   assert(sol != NULL);
   assert(var != NULL);

   CHECK_OKAY( solEnsureValsMem(sol, memhdr, set, var) );

   pos = var->index - sol->firstindex;
   assert(pos >= 0 && pos < sol->nvals);
   assert(sol->vars[pos] == var);

   sol->obj -= var->obj * sol->vals[pos];
   sol->vals[pos] = val;
   sol->obj += var->obj * val;

   return SCIP_OKAY;
}

RETCODE SCIPsolIncVal(                  /**< increases value of variable in primal IP solution */
   SOL*             sol,                /**< primal IP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< variable to add to solution */
   Real             incval              /**< increment for solution value of variable */
   )
{
   int pos;

   assert(sol != NULL);
   assert(var != NULL);

   CHECK_OKAY( solEnsureValsMem(sol, memhdr, set, var) );

   pos = var->index - sol->firstindex;
   assert(pos >= 0 && pos < sol->nvals);
   assert(sol->vars[pos] == var);

   sol->vals[pos] += incval;
   sol->obj += var->obj * incval;

   return SCIP_OKAY;
}

RETCODE SCIPsolCopyLPSol(               /**< copys LP solution to primal IP solution */
   SOL*             sol,                /**< primal IP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(sol != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);

   SCIPsolClear(sol, memhdr, set, lp);

   if( lp->ncols > 0 )
   {
      VAR* firstvar;
      VAR* lastvar;
      VAR* var;
      int firstindex;
      int lastindex;
      int index;
      int pos;
      int c;

      /* find the variables in LP with smallest and largest index */
      firstindex = INT_MAX;
      lastindex = INT_MIN;
      for( c = 0; c < lp->ncols; ++c )
      {
         assert(lp->cols[c] != NULL);
         var = lp->cols[c]->var;
         assert(var != NULL);
         assert(var->varstatus == SCIP_VARSTATUS_COLUMN);
         assert(var->data.col == lp->cols[c]);
         index = var->index;

         if( index < firstindex )
         {
            firstindex = index;
            firstvar = var;
         }
         if( index > lastindex )
         {
            lastindex = index;
            lastvar = var;
         }
      }
      assert(firstvar != NULL && firstvar->index == firstindex);
      assert(lastvar != NULL && lastvar->index == lastindex);

      /* get memory to store all variables of LP */
      CHECK_OKAY( solEnsureValsMem(sol, memhdr, set, firstvar) );
      CHECK_OKAY( solEnsureValsMem(sol, memhdr, set, lastvar) );
      assert(sol->firstindex == firstindex);
      assert(sol->firstindex + sol->nvals > lastindex);
      assert(sol->obj == 0.0);

      /* store the variables of LP */
      for( c = 0; c < lp->ncols; ++c )
      {
         assert(lp->cols[c] != NULL);
         var = lp->cols[c]->var;
         assert(var != NULL);
         assert(var->varstatus == SCIP_VARSTATUS_COLUMN);
         assert(var->data.col == lp->cols[c]);
         pos = var->index - firstindex;

         sol->vars[pos] = var;
         sol->vals[pos] = var->data.col->primsol;
         sol->obj += var->obj * sol->vals[pos];
      }
   }

   return SCIP_OKAY;
}

Real SCIPsolGetVal(                     /**< returns value of variable in primal IP solution */
   SOL*             sol,                /**< primal IP solution */
   VAR*             var                 /**< variable to get value for */
   )
{
   int pos;

   assert(sol != NULL);
   assert(var != NULL);

   pos = var->index - sol->firstindex;

   if( pos < 0 || pos >= sol->nvals)
      return 0.0;
   else
   {
      assert(sol->vals[pos] == 0.0 || sol->vars[pos] == var);
      return sol->vals[pos];
   }
}

void SCIPsolPrint(                      /**< outputs non-zero elements of solution to file stream */
   SOL*             sol,                /**< primal IP solution */
   const SET*       set,                /**< global SCIP settings */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   int i;

   assert(sol != NULL);
   
   if( file == NULL )
      file = stdout;

   fprintf(file, "obj=%g", sol->obj);
   for( i = 0; i < sol->nvals; ++i )
   {
      if( !SCIPsetIsZero(set, sol->vals[i]) )
      {
         assert(sol->vars[i] != NULL);
         fprintf(file, ", %s=%g", sol->vars[i]->name, sol->vals[i]);
      }
   }
   fprintf(file, "\n");
}


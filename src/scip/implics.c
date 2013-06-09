/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   implics.c
 * @brief  methods for implications, variable bounds, and clique tables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>
#include "string.h"

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/event.h"
#include "scip/var.h"
#include "scip/implics.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/debug.h"

#ifndef NDEBUG
#include "scip/struct_implics.h"
#endif



/*
 * methods for variable bounds
 */

/** creates a variable bounds data structure */
static
SCIP_RETCODE vboundsCreate(
   SCIP_VBOUNDS**        vbounds,            /**< pointer to store variable bounds data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(vbounds != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, vbounds) );
   (*vbounds)->vars = NULL;
   (*vbounds)->coefs = NULL;
   (*vbounds)->constants = NULL;
   (*vbounds)->len = 0;
   (*vbounds)->size = 0;

   return SCIP_OKAY;
}

/** frees a variable bounds data structure */
void SCIPvboundsFree(
   SCIP_VBOUNDS**        vbounds,            /**< pointer to store variable bounds data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(vbounds != NULL);

   if( *vbounds != NULL )
   {
      BMSfreeBlockMemoryArrayNull(blkmem, &(*vbounds)->vars, (*vbounds)->size);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*vbounds)->coefs, (*vbounds)->size);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*vbounds)->constants, (*vbounds)->size);
      BMSfreeBlockMemory(blkmem, vbounds);
   }
}

/** ensures, that variable bounds arrays can store at least num entries */
static
SCIP_RETCODE vboundsEnsureSize(
   SCIP_VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(vbounds != NULL);
   
   /* create variable bounds data structure, if not yet existing */
   if( *vbounds == NULL )
   {
      SCIP_CALL( vboundsCreate(vbounds, blkmem) );
   }
   assert(*vbounds != NULL);
   assert((*vbounds)->len <= (*vbounds)->size);

   if( num > (*vbounds)->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(*vbounds)->vars, (*vbounds)->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(*vbounds)->coefs, (*vbounds)->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(*vbounds)->constants, (*vbounds)->size, newsize) );
      (*vbounds)->size = newsize;
   }
   assert(num <= (*vbounds)->size);

   return SCIP_OKAY;
}

/** binary searches the insertion position of the given variable in the vbounds data structure */
static
SCIP_RETCODE vboundsSearchPos(
   SCIP_VBOUNDS*         vbounds,            /**< variable bounds data structure, or NULL */
   SCIP_VAR*             var,                /**< variable to search in vbounds data structure */
   SCIP_Bool             negativecoef,       /**< is coefficient b negative? */
   int*                  insertpos,          /**< pointer to store position where to insert new entry */
   SCIP_Bool*            found               /**< pointer to store whether the same variable was found at the returned pos */
   )
{
   SCIP_Bool exists;
   int pos;

   assert(insertpos != NULL);
   assert(found != NULL);

   /* check for empty vbounds data */
   if( vbounds == NULL )
   {
      *insertpos = 0;
      *found = FALSE;
      return SCIP_OKAY;
   }
   assert(vbounds->len >= 0);

   /* binary search for the given variable */
   exists = SCIPsortedvecFindPtr((void**)vbounds->vars, SCIPvarComp, var, vbounds->len, &pos);

   if( exists )
   {
      /* we found the variable: check if the sign of the coefficient matches */
      assert(var == vbounds->vars[pos]);
      if( (vbounds->coefs[pos] < 0.0) == negativecoef )
      {
         /* the variable exists with the same sign at the current position */
         *insertpos = pos;
         *found = TRUE;
      }
      else if( negativecoef )
      {
         if( pos+1 < vbounds->len && vbounds->vars[pos+1] == var )
         {
            /* the variable exists with the desired sign at the next position */
            assert(vbounds->coefs[pos+1] < 0.0);
            *insertpos = pos+1;
            *found = TRUE;
         }
         else
         {
            /* the negative coefficient should be inserted to the right of the positive one */
            *insertpos = pos+1;
            *found = FALSE;
         }
      }
      else
      {
         if( pos-1 >= 0 && vbounds->vars[pos-1] == var )
         {
            /* the variable exists with the desired sign at the previous position */
            assert(vbounds->coefs[pos-1] > 0.0);
            *insertpos = pos-1;
            *found = TRUE;
         }
         else
         {
            /* the positive coefficient should be inserted to the left of the negative one */
            *insertpos = pos;
            *found = FALSE;
         }
      }
   }
   else
   {
      *insertpos = pos;
      *found = FALSE;
   }

   return SCIP_OKAY;
}

/** adds a variable bound to the variable bounds data structure */
SCIP_RETCODE SCIPvboundsAdd(
   SCIP_VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BOUNDTYPE        vboundtype,         /**< type of variable bound (LOWER or UPPER) */
   SCIP_VAR*             var,                /**< variable z    in x <= b*z + d  or  x >= b*z + d */
   SCIP_Real             coef,               /**< coefficient b in x <= b*z + d  or  x >= b*z + d */
   SCIP_Real             constant,           /**< constant d    in x <= b*z + d  or  x >= b*z + d */
   SCIP_Bool*            added               /**< pointer to store whether the variable bound was added */
   )
{
   int insertpos;
   SCIP_Bool found;

   assert(vbounds != NULL);
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);
   assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS);
   assert(added != NULL);
   assert(!SCIPsetIsZero(set, coef));

   *added = FALSE;

   /* identify insertion position of variable */
   SCIP_CALL( vboundsSearchPos(*vbounds, var, (coef < 0.0), &insertpos, &found) );
   if( found )
   {
      /* the same variable already exists in the vbounds data structure: use the better vbound */
      assert(*vbounds != NULL);
      assert(0 <= insertpos && insertpos < (*vbounds)->len);
      assert((*vbounds)->vars[insertpos] == var);
      assert(((*vbounds)->coefs[insertpos] < 0.0) == (coef < 0.0));

      if( vboundtype == SCIP_BOUNDTYPE_UPPER )
      {
         if( constant + MIN(coef, 0.0) < (*vbounds)->constants[insertpos] + MIN((*vbounds)->coefs[insertpos], 0.0) )
         {
            (*vbounds)->coefs[insertpos] = coef;
            (*vbounds)->constants[insertpos] = constant;
            *added = TRUE;
         }
      }
      else
      {
         if( constant + MAX(coef, 0.0) > (*vbounds)->constants[insertpos] + MAX((*vbounds)->coefs[insertpos], 0.0) )
         {
            (*vbounds)->coefs[insertpos] = coef;
            (*vbounds)->constants[insertpos] = constant;
            *added = TRUE;
         }
      }
   }
   else
   {
      int i;

      /* the given variable does not yet exist in the vbounds */
      SCIP_CALL( vboundsEnsureSize(vbounds, blkmem, set, *vbounds != NULL ? (*vbounds)->len+1 : 1) );
      assert(*vbounds != NULL);
      assert(0 <= insertpos && insertpos <= (*vbounds)->len);
      assert(0 <= insertpos && insertpos < (*vbounds)->size);

      /* insert variable at the correct position */
      for( i = (*vbounds)->len; i > insertpos; --i )
      {
         assert(!SCIPsetIsZero(set, (*vbounds)->coefs[i-1]));
         (*vbounds)->vars[i] = (*vbounds)->vars[i-1];
         (*vbounds)->coefs[i] = (*vbounds)->coefs[i-1];
         (*vbounds)->constants[i] = (*vbounds)->constants[i-1];
      }
      assert(!SCIPsetIsZero(set, coef));
      (*vbounds)->vars[insertpos] = var;
      (*vbounds)->coefs[insertpos] = coef;
      (*vbounds)->constants[insertpos] = constant;
      (*vbounds)->len++;
      *added = TRUE;
   }

   return SCIP_OKAY;
}

/** removes from variable x a variable bound x >=/<= b*z + d with binary or integer z */
SCIP_RETCODE SCIPvboundsDel(
   SCIP_VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             vbdvar,             /**< variable z    in x >=/<= b*z + d */
   SCIP_Bool             negativecoef        /**< is coefficient b negative? */
   )
{
   SCIP_Bool found;
   int pos;
   int i;

   assert(vbounds != NULL);
   assert(*vbounds != NULL);

   /* searches for variable z in variable bounds of x */
   SCIP_CALL( vboundsSearchPos(*vbounds, vbdvar, negativecoef, &pos, &found) );
   if( !found )
      return SCIP_OKAY;

   assert(0 <= pos && pos < (*vbounds)->len);
   assert((*vbounds)->vars[pos] == vbdvar);
   assert(((*vbounds)->coefs[pos] < 0.0) == negativecoef);

   /* removes z from variable bounds of x */
   for( i = pos; i < (*vbounds)->len - 1; i++ )
   {
      (*vbounds)->vars[i] = (*vbounds)->vars[i+1];
      (*vbounds)->coefs[i] = (*vbounds)->coefs[i+1];
      (*vbounds)->constants[i] = (*vbounds)->constants[i+1];
   }
   (*vbounds)->len--;

#ifndef NDEBUG
   SCIP_CALL( vboundsSearchPos(*vbounds, vbdvar, negativecoef, &pos, &found) );
   assert(!found);
#endif

   /* free vbounds data structure, if it is empty */
   if( (*vbounds)->len == 0 )
      SCIPvboundsFree(vbounds, blkmem);

   return SCIP_OKAY;
}

/** reduces the number of variable bounds stored in the given variable bounds data structure */
void SCIPvboundsShrink(
   SCIP_VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   newnvbds            /**< new number of variable bounds */
   )
{
   assert(vbounds != NULL);
   assert(*vbounds != NULL);
   assert(newnvbds <= (*vbounds)->len);

   if( newnvbds == 0 )
      SCIPvboundsFree(vbounds, blkmem);
   else
      (*vbounds)->len = newnvbds;
}




/*
 * methods for implications
 */

#ifndef NDEBUG
/** comparator function for implication variables in the implication data structure */
static
SCIP_DECL_SORTPTRCOMP(compVars)
{  /*lint --e{715}*/
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VARTYPE var1type;
   SCIP_VARTYPE var2type;
   int var1idx;
   int var2idx;

   var1 = (SCIP_VAR*)elem1;
   var2 = (SCIP_VAR*)elem2;
   assert(var1 != NULL);
   assert(var2 != NULL);
   var1type = SCIPvarGetType(var1);
   var2type = SCIPvarGetType(var2);
   var1idx = SCIPvarGetIndex(var1);
   var2idx = SCIPvarGetIndex(var2);

   if( var1type == var2type )
   {
      if( var1idx < var2idx )
         return -1;
      else if( var1idx > var2idx )
         return +1;
      else
         return 0;
   }
   else
   {
      if( var1type == SCIP_VARTYPE_BINARY )
         return -1;
      if( var2type == SCIP_VARTYPE_BINARY )
         return +1;
      else if( var1idx < var2idx )
         return -1;
      else if( var1idx > var2idx )
         return +1;
      else
      {
         SCIPABORT();
         /*lint --e{527}*/
         return 0;
      }
   }
}

/** performs integrity check on implications data structure */
static
void checkImplics(
   SCIP_IMPLICS*         implics,            /**< implications data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_Bool varfixing;

   if( implics == NULL )
      return;

   varfixing = FALSE;
   do
   {
      SCIP_VAR** vars;
      SCIP_BOUNDTYPE* types;
      SCIP_Real* bounds;
      int nimpls;
      int nbinimpls;
      int i;

      vars = implics->vars[varfixing];
      types = implics->types[varfixing];
      bounds = implics->bounds[varfixing];
      nimpls = implics->nimpls[varfixing];
      nbinimpls = implics->nbinimpls[varfixing];

      assert(0 <= nbinimpls && nbinimpls <= nimpls && nimpls <= implics->size[varfixing]);
      assert(nimpls == 0 || vars != NULL);
      assert(nimpls == 0 || types != NULL);
      assert(nimpls == 0 || bounds != NULL);

      if( nbinimpls > 0 )
      {
         /* in case of implication we cannot use SCIPvarIsBinary() to check for binaries since the implication are
          * sorted with respect to variable type; this means first the binary variables (SCIPvarGetType(var) ==
          * SCIP_VARTYPE_BINARY) and second all others;
          */
         assert(SCIPvarGetType(vars[0]) == SCIP_VARTYPE_BINARY);
         assert((types[0] == SCIP_BOUNDTYPE_LOWER) == (bounds[0] > 0.5));
         assert(SCIPsetIsFeasZero(set, bounds[0]) || SCIPsetIsFeasEQ(set, bounds[0], 1.0));
      }

      for( i = 1; i < nbinimpls; ++i )
      {
         int cmp;

         /* in case of implication we cannot use SCIPvarIsBinary() to check for binaries since the implication are
          * sorted with respect to variable type; this means first the binary variables (SCIPvarGetType(var) ==
          * SCIP_VARTYPE_BINARY) and second all others;
          */
         assert(SCIPvarGetType(vars[i]) == SCIP_VARTYPE_BINARY);
         assert((types[i] == SCIP_BOUNDTYPE_LOWER) == (bounds[i] > 0.5));
         assert(SCIPsetIsFeasZero(set, bounds[i]) || SCIPsetIsFeasEQ(set, bounds[i], 1.0));

         cmp = compVars(vars[i-1], vars[i]);
         assert(cmp <= 0);
         assert((cmp == 0) == (vars[i-1] == vars[i]));
         assert(cmp < 0 || (types[i-1] == SCIP_BOUNDTYPE_LOWER && types[i] == SCIP_BOUNDTYPE_UPPER));
      }

      for( i = nbinimpls; i < nimpls; ++i )
      {
         int cmp;

         /* in case of implication we cannot use SCIPvarIsBinary() to check for binaries since the implication are
          * sorted with respect to variable type; this means first the binary variables (SCIPvarGetType(var) ==
          * SCIP_VARTYPE_BINARY) and second all others;
          */
         assert(SCIPvarGetType(vars[i]) != SCIP_VARTYPE_BINARY);

         if( i == 0 )
            continue;

         cmp = compVars(vars[i-1], vars[i]);
         assert(cmp <= 0);
         assert((cmp == 0) == (vars[i-1] == vars[i]));
         assert(cmp < 0 || (types[i-1] == SCIP_BOUNDTYPE_LOWER && types[i] == SCIP_BOUNDTYPE_UPPER));
      }

      varfixing = !varfixing;
   }
   while( varfixing == TRUE );
}
#else
#define checkImplics(implics,set) /**/
#endif

/** creates an implications data structure */
static
SCIP_RETCODE implicsCreate(
   SCIP_IMPLICS**        implics,            /**< pointer to store implications data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(implics != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, implics) );

   (*implics)->vars[0] = NULL;
   (*implics)->types[0] = NULL;
   (*implics)->bounds[0] = NULL;
   (*implics)->ids[0] = NULL;
   (*implics)->size[0] = 0;
   (*implics)->nimpls[0] = 0;
   (*implics)->nbinimpls[0] = 0;

   (*implics)->vars[1] = NULL;
   (*implics)->types[1] = NULL;
   (*implics)->bounds[1] = NULL;
   (*implics)->ids[1] = NULL;
   (*implics)->size[1] = 0;
   (*implics)->nimpls[1] = 0;
   (*implics)->nbinimpls[1] = 0;

   return SCIP_OKAY;
}

/** frees an implications data structure */
void SCIPimplicsFree(
   SCIP_IMPLICS**        implics,            /**< pointer of implications data structure to free */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(implics != NULL);

   if( *implics != NULL )
   {
      BMSfreeBlockMemoryArrayNull(blkmem, &(*implics)->vars[0], (*implics)->size[0]);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*implics)->types[0], (*implics)->size[0]);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*implics)->bounds[0], (*implics)->size[0]);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*implics)->ids[0], (*implics)->size[0]);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*implics)->vars[1], (*implics)->size[1]);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*implics)->types[1], (*implics)->size[1]);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*implics)->bounds[1], (*implics)->size[1]);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*implics)->ids[1], (*implics)->size[1]);
      BMSfreeBlockMemory(blkmem, implics);
   }
}

/** ensures, that arrays for x == 0 or x == 1 in implications data structure can store at least num entries */
static
SCIP_RETCODE implicsEnsureSize(
   SCIP_IMPLICS**        implics,            /**< pointer to implications data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             varfixing,          /**< FALSE if size of arrays for x == 0 has to be ensured, TRUE for x == 1 */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(implics != NULL);
   
   /* create implications data structure, if not yet existing */
   if( *implics == NULL )
   {
      SCIP_CALL( implicsCreate(implics, blkmem) );
   }
   assert(*implics != NULL);
   assert((*implics)->nimpls[varfixing] <= (*implics)->size[varfixing]);

   if( num > (*implics)->size[varfixing] )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(*implics)->vars[varfixing], (*implics)->size[varfixing],
            newsize) ); /*lint !e866*/
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(*implics)->types[varfixing], (*implics)->size[varfixing], 
            newsize) ); /*lint !e866*/
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(*implics)->bounds[varfixing], (*implics)->size[varfixing],
            newsize) ); /*lint !e866*/
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(*implics)->ids[varfixing], (*implics)->size[varfixing],
            newsize) ); /*lint !e866*/
      (*implics)->size[varfixing] = newsize;
   }
   assert(num <= (*implics)->size[varfixing]);

   return SCIP_OKAY;
}

/** gets the positions of the implications y >= l and y <= u in the implications data structure;
 *  if no lower or upper bound implication for y was found, -1 is returned
 */
static
void implicsSearchVar(
   SCIP_IMPLICS*         implics,            /**< implications data structure */
   SCIP_Bool             varfixing,          /**< FALSE if y is searched in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y to search for */
   int*                  poslower,           /**< pointer to store position of y_lower (-1 if not found) */
   int*                  posupper,           /**< pointer to store position of y_upper (-1 if not found) */
   int*                  posadd              /**< pointer to store position of first y entry, or where a new y entry
                                              *   should be placed */
   )
{
   SCIP_Bool found;
   int left;
   int right;
   int pos;

   assert(implics != NULL);
   assert(poslower != NULL);
   assert(posupper != NULL);
   assert(posadd != NULL);

   /* set left and right pointer */
   if( SCIPvarGetType(implvar) == SCIP_VARTYPE_BINARY )
   {
      /* check whether there are implications with binary variable y */
      if( implics->nbinimpls[varfixing] == 0 )
      {
         /* there are no implications with binary variable y */
         *posadd = 0;
         *poslower = -1;
         *posupper = -1;
         return;
      }      
      left = 0;
      right = implics->nbinimpls[varfixing] - 1;
   }
   else
   {
      /* check whether there are implications with non-binary variable y */
      if( implics->nimpls[varfixing] == implics->nbinimpls[varfixing] )
      {
         /* there are no implications with non-binary variable y */
         *posadd = implics->nbinimpls[varfixing];
         *poslower = -1;
         *posupper = -1;
         return;
      }
      left = implics->nbinimpls[varfixing];
      right = implics->nimpls[varfixing] - 1;
   }
   assert(left <= right);

   /* search for the position in the sorted array (via binary search) */
   found = SCIPsortedvecFindPtr((void**)(&(implics->vars[varfixing][left])), SCIPvarComp, (void*)implvar, right-left+1, &pos);

   /* adjust position */
   pos += left;

   if( !found )
   {
      /* y was not found */
      assert(pos >= right || compVars((void*)implics->vars[varfixing][pos], (void*)implvar) > 0);
      assert(pos == left || compVars((void*)implics->vars[varfixing][pos-1], (void*)implvar) < 0);
      *poslower = -1;
      *posupper = -1;
      *posadd = pos;
   }
   else
   {
      /* y was found */
      assert(implvar == implics->vars[varfixing][pos]);

      /* set poslower and posupper */
      if( implics->types[varfixing][pos] == SCIP_BOUNDTYPE_LOWER )
      {
         /* y was found as y_lower (on position middle) */
         *poslower = pos;
         if( pos + 1 < implics->nimpls[varfixing] && implics->vars[varfixing][pos+1] == implvar )
         {  
            assert(implics->types[varfixing][pos+1] == SCIP_BOUNDTYPE_UPPER);
            *posupper = pos + 1;
         }
         else
            *posupper = -1;
         *posadd = pos;
      }
      else
      {
         /* y was found as y_upper (on position pos) */
         *posupper = pos;
         if( pos - 1 >= 0 && implics->vars[varfixing][pos-1] == implvar )
         {  
            assert(implics->types[varfixing][pos-1] == SCIP_BOUNDTYPE_LOWER);
            *poslower = pos - 1;
            *posadd = pos - 1;
         }
         else
         {
            *poslower = -1;
            *posadd = pos;
         }
      }
   }
}

/** returns whether variable y is already contained in implications for x == 0 or x == 1 with the given impltype
 *  y can be contained in structure with y >= b (y_lower) and y <= b (y_upper) 
 */
static
SCIP_Bool implicsSearchImplic(
   SCIP_IMPLICS*         implics,            /**< implications data structure */
   SCIP_Bool             varfixing,          /**< FALSE if y is searched in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y to search for */
   SCIP_BOUNDTYPE        impltype,           /**< type of implication y <=/>= b to search for */
   int*                  poslower,           /**< pointer to store position of y_lower (inf if not found) */
   int*                  posupper,           /**< pointer to store position of y_upper (inf if not found) */
   int*                  posadd              /**< pointer to store correct position (with respect to impltype) to add y */
   )
{
   assert(implics != NULL);
   assert(poslower != NULL);
   assert(posupper != NULL);
   assert(posadd != NULL);

   implicsSearchVar(implics, varfixing, implvar, poslower, posupper, posadd);
   assert(*poslower == -1 || *posupper == -1 || *posupper == (*poslower)+1);
   assert(*poslower == -1 || *posadd == *poslower);
   assert(*poslower >= 0 || *posupper == -1 || *posadd == *posupper);
   assert(0 <= *posadd && *posadd <= implics->nimpls[varfixing]);

   if( impltype == SCIP_BOUNDTYPE_LOWER )
      return (*poslower >= 0);
   else
   {
      if( *poslower >= 0 )
      {
         assert(*posadd == *poslower);
         (*posadd)++;
      }
      return (*posupper >= 0);
   }
}

/** adds an implication x == 0/1 -> y <= b or y >= b to the implications data structure;
 *  the implication must be non-redundant
 */
SCIP_RETCODE SCIPimplicsAdd(
   SCIP_IMPLICS**        implics,            /**< pointer to implications data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool             varfixing,          /**< FALSE if implication for x == 0 has to be added, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   SCIP_Real             implbound,          /**< bound b    in implication y <= b or y >= b */
   SCIP_Bool             isshortcut,         /**< is the implication a shortcut, i.e., added as part of the transitive closure of another implication? */
   SCIP_Bool*            conflict,           /**< pointer to store whether implication causes a conflict for variable x */
   SCIP_Bool*            added               /**< pointer to store whether the implication was added */
   )
{
   int poslower;
   int posupper;
   int posadd;
   SCIP_Bool found;
#ifndef NDEBUG
   int k;
#endif

   assert(implics != NULL);
   assert(*implics == NULL || (*implics)->nbinimpls[varfixing] <= (*implics)->nimpls[varfixing]);
   assert(stat != NULL);
   assert(SCIPvarIsActive(implvar));
   assert(SCIPvarGetStatus(implvar) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(implvar) == SCIP_VARSTATUS_LOOSE); 
   assert((impltype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsFeasGT(set, implbound, SCIPvarGetLbGlobal(implvar)))
      || (impltype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsFeasLT(set, implbound, SCIPvarGetUbGlobal(implvar))));
   assert(conflict != NULL);
   assert(added != NULL);

   SCIPdebugMessage("adding implication to implics %p [%u]: <%s> %s %g\n",
      (void*)*implics, varfixing, SCIPvarGetName(implvar), impltype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", implbound);

   checkImplics(*implics, set);

   *conflict = FALSE;
   *added = FALSE;

   /* check if variable is already contained in implications data structure */
   if( *implics != NULL )
   {
      found = implicsSearchImplic(*implics, varfixing, implvar, impltype, &poslower, &posupper, &posadd);
      assert(-1 <= poslower && poslower < (*implics)->nimpls[varfixing]);
      assert(-1 <= posupper && posupper < (*implics)->nimpls[varfixing]);
      assert(0 <= posadd && posadd <= (*implics)->nimpls[varfixing]);
      assert(poslower == -1 || (*implics)->types[varfixing][poslower] == SCIP_BOUNDTYPE_LOWER);
      assert(posupper == -1 || (*implics)->types[varfixing][posupper] == SCIP_BOUNDTYPE_UPPER);
   }
   else
   {
      found = FALSE;
      poslower = -1;
      posupper = -1;
      posadd = 0;
   }

   if( impltype == SCIP_BOUNDTYPE_LOWER )
   {
      assert(found == (poslower >= 0));

      /* check if y >= b is redundant */
      if( poslower >= 0 && SCIPsetIsFeasLE(set, implbound, (*implics)->bounds[varfixing][poslower]) )
      {
         SCIPdebugMessage(" -> implication is redundant to <%s> >= %g\n", 
            SCIPvarGetName(implvar), (*implics)->bounds[varfixing][poslower]);
         return SCIP_OKAY;
      }

      /* check if y >= b causes conflict for x (i.e. y <= a (with a < b) is also valid) */
      if( posupper >= 0 && SCIPsetIsFeasGT(set, implbound, (*implics)->bounds[varfixing][posupper]) )
      {      
         SCIPdebugMessage(" -> implication is conflicting to <%s> <= %g\n", 
            SCIPvarGetName(implvar), (*implics)->bounds[varfixing][posupper]);
         *conflict = TRUE;
         return SCIP_OKAY;
      }

      *added = TRUE;

      /* check if entry of the same type already exists */
      if( found )
      {
         assert(poslower >= 0);
         assert(posadd == poslower);

         /* add y >= b by changing old entry on poslower */
         assert((*implics)->vars[varfixing][poslower] == implvar);
         assert(SCIPsetIsFeasGT(set, implbound, (*implics)->bounds[varfixing][poslower]));
         (*implics)->bounds[varfixing][poslower] = implbound;
      }
      else
      {
         assert(poslower == -1);

         /* add y >= b by creating a new entry on posadd */
         SCIP_CALL( implicsEnsureSize(implics, blkmem, set, varfixing,
               *implics != NULL ? (*implics)->nimpls[varfixing]+1 : 1) );
         assert(*implics != NULL);

	 if( (*implics)->nimpls[varfixing] - posadd > 0 )
	 {
	    int amount = ((*implics)->nimpls[varfixing] - posadd);

#ifndef NDEBUG
	    for( k = (*implics)->nimpls[varfixing]; k > posadd; k-- )
	    {
	       assert(compVars((void*)(*implics)->vars[varfixing][k-1], (void*)implvar) >= 0);
	    }
#endif
	    BMSmoveMemoryArray(&((*implics)->types[varfixing][posadd+1]), &((*implics)->types[varfixing][posadd]), amount); /*lint !e866*/
	    BMSmoveMemoryArray(&((*implics)->ids[varfixing][posadd+1]), &((*implics)->ids[varfixing][posadd]), amount); /*lint !e866*/
	    BMSmoveMemoryArray(&((*implics)->vars[varfixing][posadd+1]), &((*implics)->vars[varfixing][posadd]), amount); /*lint !e866*/
	    BMSmoveMemoryArray(&((*implics)->bounds[varfixing][posadd+1]), &((*implics)->bounds[varfixing][posadd]), amount); /*lint !e866*/
	 }

         (*implics)->vars[varfixing][posadd] = implvar;
         (*implics)->types[varfixing][posadd] = impltype;
         (*implics)->bounds[varfixing][posadd] = implbound;
         (*implics)->ids[varfixing][posadd] = (isshortcut ? -stat->nimplications : stat->nimplications);
         if( SCIPvarGetType(implvar) == SCIP_VARTYPE_BINARY )
            (*implics)->nbinimpls[varfixing]++;
         (*implics)->nimpls[varfixing]++;
#ifndef NDEBUG
         for( k = posadd-1; k >= 0; k-- )
            assert(compVars((void*)(*implics)->vars[varfixing][k], (void*)implvar) <= 0);
#endif
         stat->nimplications++;
      }
   }
   else
   {
      assert(found == (posupper >= 0));

      /* check if y <= b is redundant */
      if( posupper >= 0 && SCIPsetIsFeasGE(set, implbound, (*implics)->bounds[varfixing][posupper]) )
      {
         SCIPdebugMessage(" -> implication is redundant to <%s> <= %g\n", 
            SCIPvarGetName(implvar), (*implics)->bounds[varfixing][posupper]);
         return SCIP_OKAY;
      }

      /* check if y <= b causes conflict for x (i.e. y >= a (with a > b) is also valid) */
      if( poslower >= 0 && SCIPsetIsFeasLT(set, implbound, (*implics)->bounds[varfixing][poslower]) )
      {      
         SCIPdebugMessage(" -> implication is conflicting to <%s> >= %g\n", 
            SCIPvarGetName(implvar), (*implics)->bounds[varfixing][poslower]);
         *conflict = TRUE;
         return SCIP_OKAY;
      }

      *added = TRUE;

      /* check if entry of the same type already exists */
      if( found )
      {
         assert(posupper >= 0);
         assert(posadd == posupper);

         /* add y <= b by changing old entry on posupper */
         assert((*implics)->vars[varfixing][posupper] == implvar);
         assert(SCIPsetIsFeasLT(set, implbound,(*implics)->bounds[varfixing][posupper]));
         (*implics)->bounds[varfixing][posupper] = implbound;
      }
      else
      {
         /* add y <= b by creating a new entry on posadd */
         assert(posupper == -1);

         SCIP_CALL( implicsEnsureSize(implics, blkmem, set, varfixing,
               *implics != NULL ? (*implics)->nimpls[varfixing]+1 : 1) );
         assert(*implics != NULL);

	 if( (*implics)->nimpls[varfixing] - posadd > 0 )
	 {
	    int amount = ((*implics)->nimpls[varfixing] - posadd);

#ifndef NDEBUG
	    for( k = (*implics)->nimpls[varfixing]; k > posadd; k-- )
	    {
	       assert(compVars((void*)(*implics)->vars[varfixing][k-1], (void*)implvar) >= 0);
	    }
#endif
	    BMSmoveMemoryArray(&((*implics)->types[varfixing][posadd+1]), &((*implics)->types[varfixing][posadd]), amount); /*lint !e866*/
	    BMSmoveMemoryArray(&((*implics)->ids[varfixing][posadd+1]), &((*implics)->ids[varfixing][posadd]), amount); /*lint !e866*/
	    BMSmoveMemoryArray(&((*implics)->vars[varfixing][posadd+1]), &((*implics)->vars[varfixing][posadd]), amount); /*lint !e866*/
	    BMSmoveMemoryArray(&((*implics)->bounds[varfixing][posadd+1]), &((*implics)->bounds[varfixing][posadd]), amount); /*lint !e866*/
	 }

         (*implics)->vars[varfixing][posadd] = implvar;
         (*implics)->types[varfixing][posadd] = impltype;
         (*implics)->bounds[varfixing][posadd] = implbound;
         (*implics)->ids[varfixing][posadd] = (isshortcut ? -stat->nimplications : stat->nimplications);
         if( SCIPvarGetType(implvar) == SCIP_VARTYPE_BINARY )
            (*implics)->nbinimpls[varfixing]++;
         (*implics)->nimpls[varfixing]++;
#ifndef NDEBUG
         for( k = posadd-1; k >= 0; k-- )
            assert(compVars((void*)(*implics)->vars[varfixing][k], (void*)implvar) <= 0);
#endif
         stat->nimplications++;
      }
   }
    
   checkImplics(*implics, set);

   return SCIP_OKAY;
}

/** removes the implication  x <= 0 or x >= 1  ==>  y <= b  or  y >= b  from the implications data structure */
SCIP_RETCODE SCIPimplicsDel(
   SCIP_IMPLICS**        implics,            /**< pointer to implications data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             varfixing,          /**< FALSE if y should be removed from implications for x <= 0, TRUE for x >= 1 */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_BOUNDTYPE        impltype            /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   )
{
   int poslower;
   int posupper; 
   int posadd;
   SCIP_Bool found;

   assert(implics != NULL);
   assert(*implics != NULL);
   assert(implvar != NULL);

   SCIPdebugMessage("deleting implication from implics %p [%u]: <%s> %s x\n",
      (void*)*implics, varfixing, SCIPvarGetName(implvar), impltype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=");

   checkImplics(*implics, set);

   /* searches for y in implications of x */
   found = implicsSearchImplic(*implics, varfixing, implvar, impltype, &poslower, &posupper, &posadd);
   if( !found )
   {
      SCIPdebugMessage(" -> implication was not found\n");
      return SCIP_OKAY;
   }

   assert((impltype == SCIP_BOUNDTYPE_LOWER && poslower >= 0 && posadd == poslower) 
      || (impltype == SCIP_BOUNDTYPE_UPPER && posupper >= 0 && posadd == posupper));
   assert(0 <= posadd && posadd < (*implics)->nimpls[varfixing]);
   assert((SCIPvarGetType(implvar) == SCIP_VARTYPE_BINARY) == (posadd < (*implics)->nbinimpls[varfixing]));
   assert((*implics)->vars[varfixing][posadd] == implvar);
   assert((*implics)->types[varfixing][posadd] == impltype);

   /* removes y from implications of x */
   if( (*implics)->nimpls[varfixing] - posadd > 1 )
   {
      int amount = ((*implics)->nimpls[varfixing] - posadd - 1);

      BMSmoveMemoryArray(&((*implics)->types[varfixing][posadd]), &((*implics)->types[varfixing][posadd+1]), amount); /*lint !e866*/
      BMSmoveMemoryArray(&((*implics)->vars[varfixing][posadd]), &((*implics)->vars[varfixing][posadd+1]), amount); /*lint !e866*/
      BMSmoveMemoryArray(&((*implics)->bounds[varfixing][posadd]), &((*implics)->bounds[varfixing][posadd+1]), amount); /*lint !e866*/
   }
   (*implics)->nimpls[varfixing]--;

   if( SCIPvarGetType(implvar) == SCIP_VARTYPE_BINARY )
   {
      assert(posadd < (*implics)->nbinimpls[varfixing]);
      (*implics)->nbinimpls[varfixing]--;
   }

   /* free implics data structure, if it is empty */
   if( (*implics)->nimpls[0] == 0 && (*implics)->nimpls[1] == 0 )
      SCIPimplicsFree(implics, blkmem);

   checkImplics(*implics, set);

   return SCIP_OKAY;
}

/** returns which implications on given variable y are contained in implications for x == 0 or x == 1 */
void SCIPimplicsGetVarImplics(
   SCIP_IMPLICS*         implics,            /**< implications data structure */
   SCIP_Bool             varfixing,          /**< FALSE if y should be searched in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y to search for */
   SCIP_Bool*            haslowerimplic,     /**< pointer to store whether there exists an implication y >= l */
   SCIP_Bool*            hasupperimplic      /**< pointer to store whether there exists an implication y <= u */
   )
{
   int poslower;
   int posupper;
   int posadd;

   assert(haslowerimplic != NULL);
   assert(hasupperimplic != NULL);

   implicsSearchVar(implics, varfixing, implvar, &poslower, &posupper, &posadd);

   *haslowerimplic = (poslower >= 0);
   *hasupperimplic = (posupper >= 0);
}

/** returns whether an implication y <= b or y >= b is contained in implications for x == 0 or x == 1 */
SCIP_Bool SCIPimplicsContainsImpl(
   SCIP_IMPLICS*         implics,            /**< implications data structure */
   SCIP_Bool             varfixing,          /**< FALSE if y should be searched in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y to search for */
   SCIP_BOUNDTYPE        impltype            /**< type of implication y <=/>= b to search for */
   )
{
   int poslower;
   int posupper;
   int posadd;

   return implicsSearchImplic(implics, varfixing, implvar, impltype, &poslower, &posupper, &posadd);
}




/*
 * methods for cliques
 */

/** creates a clique data structure */
static
SCIP_RETCODE cliqueCreate(
   SCIP_CLIQUE**         clique,             /**< pointer to store clique data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   size,               /**< initial size of clique */
   int                   id                  /**< unique identifier of the clique */
   )
{
   assert(clique != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, clique) );
   if( size > 0 )
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*clique)->vars, size) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*clique)->values, size) );
   }
   else
   {
      (*clique)->vars = NULL;
      (*clique)->values = NULL;
   }
   (*clique)->nvars = 0;
   (*clique)->size = size;
   (*clique)->id = (unsigned int)id;
   (*clique)->eventsissued = FALSE;

   return SCIP_OKAY;
}

/** frees a clique data structure */
static
void cliqueFree(
   SCIP_CLIQUE**         clique,             /**< pointer to store clique data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(clique != NULL);

   if( *clique != NULL )
   {
      BMSfreeBlockMemoryArrayNull(blkmem, &(*clique)->vars, (*clique)->size);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*clique)->values, (*clique)->size);
      BMSfreeBlockMemory(blkmem, clique);
   }
}

/** ensures, that clique arrays can store at least num entries */
static
SCIP_RETCODE cliqueEnsureSize(
   SCIP_CLIQUE*          clique,             /**< clique data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(clique != NULL);
   
   if( num > clique->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &clique->vars, clique->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &clique->values, clique->size, newsize) );
      clique->size = newsize;
   }
   assert(num <= clique->size);

   return SCIP_OKAY;
}

/** returns the position of the given variable/value pair in the clique; returns -1 if variable/value pair is not member
 *  of the clique
 */
int SCIPcliqueSearchVar(
   SCIP_CLIQUE*          clique,             /**< clique data structure */
   SCIP_VAR*             var,                /**< variable to search for */
   SCIP_Bool             value               /**< value of the variable in the clique */
   )
{
   int varidx;
   int left;
   int right;

   assert(clique != NULL);

   varidx = SCIPvarGetIndex(var);
   left = -1;
   right = clique->nvars;
   while( left < right-1 )
   {
      int middle;
      int idx;

      middle = (left+right)/2;
      idx = SCIPvarGetIndex(clique->vars[middle]);
      assert(idx >= 0);
      if( varidx < idx )
         right = middle;
      else if( varidx > idx )
         left = middle;
      else
      {
         assert(var == clique->vars[middle]);

         /* now watch out for the correct value */
         if( clique->values[middle] < value )
         {
            int i;
            for( i = middle+1; i < clique->nvars && clique->vars[i] == var; ++i )
            {
               if( clique->values[i] == value )
                  return i;
            }
            return -1;
         }
         if( clique->values[middle] > value )
         {
            int i;
            for( i = middle-1; i >= 0 && clique->vars[i] == var; --i )
            {
               if( clique->values[i] == value )
                  return i;
            }
            return -1;
         }
         return middle;
      }
   }

   return -1;
}

/** returns whether the given variable/value pair is member of the given clique */
SCIP_Bool SCIPcliqueHasVar(
   SCIP_CLIQUE*          clique,             /**< clique data structure */
   SCIP_VAR*             var,                /**< variable to remove from the clique */
   SCIP_Bool             value               /**< value of the variable in the clique */
   )
{
   return (SCIPcliqueSearchVar(clique, var, value) >= 0);
}

/** adds a single variable to the given clique */
SCIP_RETCODE SCIPcliqueAddVar(
   SCIP_CLIQUE*          clique,             /**< clique data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable to add to the clique */
   SCIP_Bool             value,              /**< value of the variable in the clique */
   SCIP_Bool*            doubleentry,        /**< pointer to store whether the variable and value occurs twice in the clique */
   SCIP_Bool*            oppositeentry       /**< pointer to store whether the variable with opposite value is in the clique */
   )
{
   int pos;
   int i;

   assert(clique != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarIsBinary(var));
   assert(doubleentry != NULL);
   assert(oppositeentry != NULL);

   SCIPdebugMessage("adding variable <%s> == %u to clique %u\n", SCIPvarGetName(var), value, clique->id);

   *doubleentry = FALSE;
   *oppositeentry = FALSE;

   /* allocate memory */
   SCIP_CALL( cliqueEnsureSize(clique, blkmem, set, clique->nvars+1) );

   /* search for insertion position by binary variable, note that first the entries are order after variable index and
    * second after the bool value of the corresponding variable
    */
   (void) SCIPsortedvecFindPtr((void**) clique->vars, SCIPvarComp, var, clique->nvars, &pos);

   assert(pos >= 0 && pos <= clique->nvars);
   /* remember insertion position for later, pos might change */
   i = pos;

   if( pos < clique->nvars )
   {
      const int amount = clique->nvars - pos;

      /* moving elements to correct position */
      BMSmoveMemoryArray(&(clique->vars[pos+1]), &(clique->vars[pos]), amount); /*lint !e866*/
      BMSmoveMemoryArray(&(clique->values[pos+1]), &(clique->values[pos]), amount); /*lint !e866*/
      clique->nvars++;

      /* insertion for a variable with cliquevalue FALSE */
      if( !value )
      {
	 /* find last entry with the same variable and value behind the insertion position */
	 for( ; pos < clique->nvars - 1 && clique->vars[pos + 1] == var && clique->values[pos + 1] == value; ++pos ); /*lint !e722*/

	 /* check if the same variable with other value also exists */
	 if( pos < clique->nvars - 1 && clique->vars[pos + 1] == var )
	 {
	    assert(clique->values[pos + 1] != value);
	    *oppositeentry = TRUE;
	 }

	 /* check if we found the same variable with the same value more than once */
	 if( pos != i )
	    *doubleentry = TRUE;
	 else
	 {
	    /* find last entry with the same variable and different value before the insertion position */
	    for( ; pos > 0 && clique->vars[pos - 1] == var && clique->values[pos - 1] != value; --pos ); /*lint !e722*/

	    /* check if we found the same variable with the same value more than once */
	    if( pos > 0 && clique->vars[pos - 1] == var )
	    {
	       assert(clique->values[pos - 1] == value);

	       *doubleentry = TRUE;
	    }
	    /* if we found the same variable with different value, we need to order them correctly */
	    if( pos != i )
	    {
	       assert(clique->vars[pos] == var);
	       assert(clique->values[pos] != value);

	       clique->values[pos] = value;
	       value = !value;
	    }
	 }
      }
      /* insertion for a variable with cliquevalue TRUE */
      else
      {
	 /* find last entry with the same variable and different value behind the insertion position */
	 for( ; pos < clique->nvars - 1 && clique->vars[pos + 1] == var && clique->values[pos + 1] != value; ++pos ); /*lint !e722*/

	 /* check if the same variable with other value also exists */
	 if( pos < clique->nvars - 1 && clique->vars[pos + 1] == var )
	 {
	    assert(clique->values[pos + 1] == value);
	    *doubleentry = TRUE;
	 }

	 /* check if we found the same variable with different value */
	 if( pos != i )
	 {
	    *oppositeentry = TRUE;

	    /* if we found the same variable with different value, we need to order them correctly */
	    assert(clique->vars[pos] == var);
	    assert(clique->values[pos] != value);

	    clique->values[pos] = value;
	    value = !value;
	 }
	 else
	 {
	    /* find last entry with the same variable and value before the insertion position */
	    for( ; pos > 0 && clique->vars[pos - 1] == var && clique->values[pos - 1] == value; --pos ); /*lint !e722*/

	    if( pos != i )
	       *doubleentry = TRUE;

	    /* check if we found the same variable with different value up front */
	    if( pos > 0 && clique->vars[pos - 1] == var && clique->values[pos - 1] != value )
	       *oppositeentry = TRUE;
	 }
      }
   }
   else
      clique->nvars++;

   clique->vars[i] = var;
   clique->values[i] = value;
   clique->eventsissued = FALSE;

   return SCIP_OKAY;
}

/** removes a single variable from the given clique */
SCIP_RETCODE SCIPcliqueDelVar(
   SCIP_CLIQUE*          clique,             /**< clique data structure */
   SCIP_VAR*             var,                /**< variable to remove from the clique */
   SCIP_Bool             value               /**< value of the variable in the clique */
   )
{
   int pos;

   assert(clique != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarIsBinary(var));

   SCIPdebugMessage("deleting variable <%s> == %u from clique %u\n", SCIPvarGetName(var), value, clique->id);

   /* find variable in clique */
   pos = SCIPcliqueSearchVar(clique, var, value);
   assert(0 <= pos && pos < clique->nvars);
   assert(clique->vars[pos] == var);
   assert(clique->values[pos] == value);

   /* remove entry from clique */
   for( ; pos < clique->nvars-1; ++pos )
   {
      clique->vars[pos] = clique->vars[pos+1];
      clique->values[pos] = clique->values[pos+1];
   }
   clique->nvars--;

   return SCIP_OKAY;
}

/** gets the position of the given clique in the cliques array; returns -1 if clique is not member of cliques array */
static
int cliquesSearchClique(
   SCIP_CLIQUE**         cliques,            /**< array of cliques */
   int                   ncliques,           /**< number of cliques in the cliques array */
   SCIP_CLIQUE*          clique              /**< clique to search for */
   )
{
   unsigned int cliqueid;
   int left;
   int right;

   assert(cliques != NULL || ncliques == 0);
   assert(clique != NULL);

   cliqueid = clique->id;
   left = -1;
   right = ncliques;
   while( left < right-1 )
   {
      unsigned int id;
      int middle;

      assert(cliques != NULL);
      middle = (left+right)/2;
      id = cliques[middle]->id;
      if( cliqueid < id )
         right = middle;
      else if( cliqueid > id )
         left = middle;
      else
      {
         assert(clique == cliques[middle]);
         return middle;
      }
   }

   return -1;
}

#ifndef NDEBUG
/** checks whether clique appears in all clique lists of the involved variables */
static
void cliqueCheck(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   )
{
   int i;

   assert(clique != NULL);

   for( i = 0; i < clique->nvars; ++i )
   {
      SCIP_CLIQUE** cliques;
      int ncliques;
      int pos;

      assert(i == 0 || SCIPvarGetIndex(clique->vars[i-1]) <= SCIPvarGetIndex(clique->vars[i]));
      assert(i == 0 || clique->vars[i-1] != clique->vars[i] || clique->values[i-1] <= clique->values[i]);
      ncliques = SCIPvarGetNCliques(clique->vars[i], clique->values[i]);
      cliques = SCIPvarGetCliques(clique->vars[i], clique->values[i]);
      pos = cliquesSearchClique(cliques, ncliques, clique);
      assert(0 <= pos && pos < ncliques);
      assert(cliques[pos] == clique);
   }
}
#else
#define cliqueCheck(clique) /**/
#endif

/** creates a clique list data structure */
static
SCIP_RETCODE cliquelistCreate(
   SCIP_CLIQUELIST**     cliquelist,         /**< pointer to store clique list data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(cliquelist != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, cliquelist) );
   (*cliquelist)->cliques[0] = NULL;
   (*cliquelist)->cliques[1] = NULL;
   (*cliquelist)->ncliques[0] = 0;
   (*cliquelist)->ncliques[1] = 0;
   (*cliquelist)->size[0] = 0;
   (*cliquelist)->size[1] = 0;

   return SCIP_OKAY;
}

/** frees a clique list data structure */
void SCIPcliquelistFree(
   SCIP_CLIQUELIST**     cliquelist,         /**< pointer to the clique list data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(cliquelist != NULL);

   if( *cliquelist != NULL )
   {
      BMSfreeBlockMemoryArrayNull(blkmem, &(*cliquelist)->cliques[0], (*cliquelist)->size[0]);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*cliquelist)->cliques[1], (*cliquelist)->size[1]);
      BMSfreeBlockMemory(blkmem, cliquelist);
   }
}

/** ensures, that clique list arrays can store at least num entries */
static
SCIP_RETCODE cliquelistEnsureSize(
   SCIP_CLIQUELIST*      cliquelist,         /**< clique list data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             value,              /**< value of the variable for which the clique list should be extended */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(cliquelist != NULL);

   if( num > cliquelist->size[value] )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &cliquelist->cliques[value], cliquelist->size[value], newsize) ); /*lint !e866*/
      cliquelist->size[value] = newsize;
   }
   assert(num <= cliquelist->size[value]);

   return SCIP_OKAY;
}

/** adds a clique to the clique list */
SCIP_RETCODE SCIPcliquelistAdd(
   SCIP_CLIQUELIST**     cliquelist,         /**< pointer to the clique list data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             value,              /**< value of the variable for which the clique list should be extended */
   SCIP_CLIQUE*          clique              /**< clique that should be added to the clique list */
   )
{
   unsigned int id;
   int i;

   assert(cliquelist != NULL);

   /* allocate memory */
   if( *cliquelist == NULL )
   {
      SCIP_CALL( cliquelistCreate(cliquelist, blkmem) );
   }
   assert(*cliquelist != NULL);
   SCIP_CALL( cliquelistEnsureSize(*cliquelist, blkmem, set, value, (*cliquelist)->ncliques[value]+1) );
   assert((*cliquelist)->cliques[value] != NULL);

   SCIPdebugMessage("adding clique %u to cliquelist %p value %u (length: %d)\n", 
      clique->id, (void*)*cliquelist, value, (*cliquelist)->ncliques[value]);
   
   /* insert clique into list, sorted by clique id */
   id = clique->id;
   for( i = (*cliquelist)->ncliques[value]; i > 0 && (*cliquelist)->cliques[value][i-1]->id > id; --i )
      (*cliquelist)->cliques[value][i] = (*cliquelist)->cliques[value][i-1];
   (*cliquelist)->cliques[value][i] = clique;
   (*cliquelist)->ncliques[value]++;

   return SCIP_OKAY;
}

/** removes a clique from the clique list */
SCIP_RETCODE SCIPcliquelistDel(
   SCIP_CLIQUELIST**     cliquelist,         /**< pointer to the clique list data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool             value,              /**< value of the variable for which the clique list should be reduced */
   SCIP_CLIQUE*          clique              /**< clique that should be deleted from the clique list */
   )
{
   int pos;

   assert(cliquelist != NULL);
   assert(*cliquelist != NULL);

   SCIPdebugMessage("deleting clique %u from cliquelist %p value %u (length: %d)\n", 
      clique->id, (void*)*cliquelist, value, (*cliquelist)->ncliques[value]);
   
   pos = cliquesSearchClique((*cliquelist)->cliques[value], (*cliquelist)->ncliques[value], clique);
   assert(0 <= pos && pos < (*cliquelist)->ncliques[value]);
   assert((*cliquelist)->cliques[value][pos] == clique);

   /* remove clique from list */
   for( ; pos < (*cliquelist)->ncliques[value] - 1; ++pos )
      (*cliquelist)->cliques[value][pos] = (*cliquelist)->cliques[value][pos+1];
   (*cliquelist)->ncliques[value]--;

   /* free cliquelist if it is empty */
   if( (*cliquelist)->ncliques[0] == 0 && (*cliquelist)->ncliques[1] == 0 )
      SCIPcliquelistFree(cliquelist, blkmem);

   return SCIP_OKAY;
}

/** returns whether the given clique lists have a non-empty intersection, i.e. whether there is a clique that appears
 *  in both lists
 */
SCIP_Bool SCIPcliquelistsHaveCommonClique(
   SCIP_CLIQUELIST*      cliquelist1,        /**< first clique list data structure */
   SCIP_Bool             value1,             /**< value of first variable */
   SCIP_CLIQUELIST*      cliquelist2,        /**< second clique list data structure */
   SCIP_Bool             value2              /**< value of second variable */
   )
{
   SCIP_CLIQUE** cliques1;
   SCIP_CLIQUE** cliques2;
   int ncliques1;
   int ncliques2;
   int i1;
   int i2;

   if( cliquelist1 == NULL || cliquelist2 == NULL )
      return FALSE;

   ncliques1 = cliquelist1->ncliques[value1];
   cliques1 = cliquelist1->cliques[value1];
   ncliques2 = cliquelist2->ncliques[value2];
   cliques2 = cliquelist2->cliques[value2];

   i1 = 0;
   i2 = 0;

   if( i1 < ncliques1 && i2 < ncliques2 )
   {
      int cliqueid;

      /* make the bigger clique the first one */
      if( ncliques2 > ncliques1 )
      {
         SCIP_CLIQUE** tmpc;
         int tmpi;

         tmpc = cliques1;
         tmpi = ncliques1;
         cliques1 = cliques2;
         ncliques1 = ncliques2;
         cliques2 = tmpc;
         ncliques2 = tmpi;
      }

      /* check whether both clique lists have a same clique */
      while( TRUE )  /*lint !e716*/
      {
	 cliqueid = SCIPcliqueGetId(cliques2[i2]);

	 /* if last item in clique1 has a smaller index than the actual clique in clique2, than cause of increasing order
	  * there will be no same item and we can stop */
	 if( SCIPcliqueGetId(cliques1[ncliques1 - 1]) < cliqueid )
	    break;

	 while( SCIPcliqueGetId(cliques1[i1]) < cliqueid )
	 {
	    ++i1;
	    assert(i1 < ncliques1);
	 }
	 cliqueid = SCIPcliqueGetId(cliques1[i1]);

	 /* if last item in clique2 has a smaller index than the actual clique in clique1, than cause of increasing order
	  * there will be no same item and we can stop */
	 if( SCIPcliqueGetId(cliques2[ncliques2 - 1]) < cliqueid )
	    break;

	 while( SCIPcliqueGetId(cliques2[i2]) < cliqueid )
	 {
	    ++i2;
	    assert(i2 < ncliques2);
	 }
	 if( SCIPcliqueGetId(cliques2[i2]) == cliqueid )
	    return TRUE;
      }
   }
   return FALSE;
}

/** removes all listed entries from the cliques */
void SCIPcliquelistRemoveFromCliques(
   SCIP_CLIQUELIST*      cliquelist,         /**< clique list data structure */
   SCIP_VAR*             var                 /**< active problem variable the clique list belongs to */
   )
{
   assert(SCIPvarIsBinary(var));

   if( cliquelist != NULL )
   {
      int value;

      SCIPdebugMessage("removing variable <%s> from cliques (%d with value 0, %d with value 1)\n",
         SCIPvarGetName(var), cliquelist->ncliques[0], cliquelist->ncliques[1]);

      for( value = 0; value < 2; ++value )
      {
         int i;

         assert(SCIPvarGetCliques(var, (SCIP_Bool)value) == cliquelist->cliques[value]);
         assert(SCIPvarGetNCliques(var, (SCIP_Bool)value) == cliquelist->ncliques[value]);
         for( i = 0; i < cliquelist->ncliques[value]; ++i )
         {
            SCIP_CLIQUE* clique;
            int pos;

            clique = cliquelist->cliques[value][i];
            assert(clique != NULL);

            SCIPdebugMessage(" -> removing variable <%s> == %d from clique %u (size %d)\n",
               SCIPvarGetName(var), value, clique->id, clique->nvars);

            /* binary search the position of the variable in the clique */
            pos = SCIPcliqueSearchVar(clique, var, (SCIP_Bool)value);
            assert(0 <= pos && pos < clique->nvars);
            assert(clique->vars[pos] == var);
            assert(clique->values[pos] == (SCIP_Bool)value);

            /* remove the entry from the clique */
            if( clique->nvars - pos - 1 > 0 )
            {
               BMSmoveMemoryArray(&(clique->vars[pos]), &(clique->vars[pos+1]), clique->nvars - pos - 1); /*lint !e866*/
               BMSmoveMemoryArray(&(clique->values[pos]), &(clique->values[pos+1]), clique->nvars - pos - 1); /*lint !e866*/
            }
            clique->nvars--;

            cliqueCheck(clique);
         }
      }
   }
}

/** creates a clique table data structure */
SCIP_RETCODE SCIPcliquetableCreate(
   SCIP_CLIQUETABLE**    cliquetable         /**< pointer to store clique table data structure */
   )
{
   assert(cliquetable != NULL);

   SCIP_ALLOC( BMSallocMemory(cliquetable) );
   (*cliquetable)->cliques = NULL;
   (*cliquetable)->ncliques = 0;
   (*cliquetable)->size = 0;
   (*cliquetable)->ncreatedcliques = 0;
   (*cliquetable)->ncleanupfixedvars = 0;
   (*cliquetable)->ncleanupaggrvars = 0;
   (*cliquetable)->ncleanupcliques = 0;

   return SCIP_OKAY;
}

/** frees a clique table data structure */
SCIP_RETCODE SCIPcliquetableFree(
   SCIP_CLIQUETABLE**    cliquetable,        /**< pointer to store clique table data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int i;

   assert(cliquetable != NULL);
   assert(*cliquetable != NULL);

   /* free all cliques */
   for( i = 0; i < (*cliquetable)->ncliques; ++i )
   {
      cliqueFree(&(*cliquetable)->cliques[i], blkmem);
   }

   /* free clique table data */
   BMSfreeMemoryArrayNull(&(*cliquetable)->cliques);
   BMSfreeMemory(cliquetable);

   return SCIP_OKAY;
}

/** ensures, that clique table arrays can store at least num entries */
static
SCIP_RETCODE cliquetableEnsureSize(
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(cliquetable != NULL);

   if( num > cliquetable->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&cliquetable->cliques, newsize) );
      cliquetable->size = newsize;
   }
   assert(num <= cliquetable->size);

   return SCIP_OKAY;
}

/** adds a clique to the clique table, using the given values for the given variables;
 *  performs implications if the clique contains the same variable twice
 */
SCIP_RETCODE SCIPcliquetableAdd(
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR**            vars,               /**< binary variables in the clique: at most one can be set to the given value */
   SCIP_Bool*            values,             /**< values of the variables in the clique; NULL to use TRUE for all vars */
   int                   nvars,              /**< number of variables in the clique */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to count the number of performed bound changes, or NULL */
   )
{
   SCIP_CLIQUE* clique;
   int i;

   assert(cliquetable != NULL);
   assert(vars != NULL);

   SCIPdebugMessage("adding clique %d with %d vars to clique table\n", cliquetable->ncliques, nvars);

   /* create the clique data structure */
   SCIP_CALL( cliqueCreate(&clique, blkmem, nvars, cliquetable->ncreatedcliques) );
   cliquetable->ncreatedcliques++;

   /* add clique to clique table */
   SCIP_CALL( cliquetableEnsureSize(cliquetable, set, cliquetable->ncliques+1) );
   cliquetable->cliques[cliquetable->ncliques] = clique;
   cliquetable->ncliques++;

   /* check clique on debugging solution */
   SCIP_CALL( SCIPdebugCheckClique(set, vars, values, nvars) ); /*lint !e506 !e774*/

   /* add the corresponding active problem variables to the clique */
   for( i = 0; i < nvars; ++i )
   {
      /* put the clique into the sorted clique table of the variable */
      SCIP_CALL( SCIPvarAddClique(vars[i], blkmem, set, stat, lp, branchcand, eventqueue,
            values != NULL ? values[i] : TRUE, clique, infeasible, nbdchgs) );
   }

   cliqueCheck(clique);

   return SCIP_OKAY;
}

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashgetkeyClique)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff both keys are equal */
static
SCIP_DECL_HASHKEYEQ(hashkeyeqClique)
{  /*lint --e{715}*/
   SCIP_CLIQUE* clique1;
   SCIP_CLIQUE* clique2;
   int i;

   clique1 = (SCIP_CLIQUE*)key1;
   clique2 = (SCIP_CLIQUE*)key2;
   assert(clique1 != NULL);
   assert(clique2 != NULL);

   if( clique1->nvars != clique2->nvars )
      return FALSE;

   /* the variables are sorted: we can simply check the equality of each pair of variable/values */
   for( i = 0; i < clique1->nvars; ++i )
   {
      if( clique1->vars[i] != clique2->vars[i] || clique1->values[i] != clique2->values[i] )
         return FALSE;
   }

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashkeyvalClique)
{  /*lint --e{715}*/
   SCIP_CLIQUE* clique;
   unsigned int hashval;
   int i;

   clique = (SCIP_CLIQUE*)key;
   hashval = 0;
   for( i = 0; i < clique->nvars; ++i )
   {
      hashval *= 31;
      hashval += (unsigned int)(((size_t)clique->vars[i]) >> 1) + (unsigned int)clique->values[i];
   }

   return hashval;
}

/** removes all empty and single variable cliques from the clique table, and converts all two variable cliques
 *  into implications; removes double entries from the clique table
 */
SCIP_RETCODE SCIPcliquetableCleanup(
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected */
   )
{
   SCIP_HASHTABLE* hashtable;
   int hashtablesize;
   int i;

   assert(cliquetable != NULL);
   assert(stat != NULL);
   assert(infeasible != NULL);

   *infeasible = FALSE;

   /* check if we have anything to do */
   if( stat->npresolfixedvars == cliquetable->ncleanupfixedvars
      && stat->npresolaggrvars == cliquetable->ncleanupaggrvars
      && cliquetable->ncliques == cliquetable->ncleanupcliques )
      return SCIP_OKAY;

   /* delay events */
   SCIP_CALL( SCIPeventqueueDelay(eventqueue) );

   /* create hash table to test for multiple cliques */
   hashtablesize = SCIPcalcHashtableSize(10*cliquetable->ncliques);
   hashtablesize = MAX(hashtablesize, (set->misc_usesmalltables ? SCIP_HASHSIZE_CLIQUES_SMALL : SCIP_HASHSIZE_CLIQUES));
   SCIP_CALL( SCIPhashtableCreate(&hashtable, blkmem, hashtablesize,
         hashgetkeyClique, hashkeyeqClique, hashkeyvalClique, NULL) );

   i = 0;
   while( i < cliquetable->ncliques && !(*infeasible) )
   {
      SCIP_CLIQUE* clique;

      clique = cliquetable->cliques[i];
      if( clique->nvars == 2 )
      {
         /* add the 2-clique as implication (don't use transitive closure; otherwise new cliques can be generated) */
         if( SCIPvarGetType(clique->vars[0]) == SCIP_VARTYPE_BINARY )
         {
            SCIP_CALL( SCIPvarAddImplic(clique->vars[0], blkmem, set, stat, lp, cliquetable, branchcand, eventqueue,
                  clique->values[0], clique->vars[1], clique->values[1] ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER,
                  (SCIP_Real)(!clique->values[1]), FALSE, infeasible, NULL) );
         }
	 else if( SCIPvarGetType(clique->vars[1]) == SCIP_VARTYPE_BINARY )
         {
            SCIP_CALL( SCIPvarAddImplic(clique->vars[1], blkmem, set, stat, lp, cliquetable, branchcand, eventqueue,
                  clique->values[1], clique->vars[0], clique->values[0] ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER,
                  (SCIP_Real)(!clique->values[0]), FALSE, infeasible, NULL) );
         }
         else
         {
            /* in case the variable are not of binary type we have to add the implication as variable bound */

            assert(SCIPvarGetType(clique->vars[0]) != SCIP_VARTYPE_BINARY && SCIPvarIsBinary(clique->vars[0]));

            /* add variable upper or rather variable lower bound on vars[0] */
            if( clique->values[0] )
            {
               SCIP_CALL( SCIPvarAddVub(clique->vars[0], blkmem, set, stat, lp, cliquetable, branchcand, eventqueue,
                     clique->vars[1], clique->values[1] ? -1.0 : 1.0, clique->values[1] ? 1.0 : 0.0, FALSE, infeasible, NULL) );
            }
            else
            {
               SCIP_CALL( SCIPvarAddVlb(clique->vars[0], blkmem, set, stat, lp, cliquetable, branchcand, eventqueue, 
                     clique->vars[1], clique->values[1] ? 1.0 : -1.0, clique->values[1] ? 0.0 : 1.0, FALSE, infeasible, NULL) );
            }
         }
      }
      
      /* check if the clique is already contained in the clique table, or if it is redundant (too small) */
      if( clique->nvars <= 2 || SCIPhashtableExists(hashtable, (void*)clique) )
      {
         int j;

         /* delete the clique from the variables' clique lists */
         for( j = 0; j < clique->nvars; ++j )
         {
            SCIP_CALL( SCIPvarDelCliqueFromList(clique->vars[j], blkmem, clique->values[j], clique) );
         }

         /* free clique and remove it from clique table */
         cliqueFree(&cliquetable->cliques[i], blkmem);
         cliquetable->cliques[i] = cliquetable->cliques[cliquetable->ncliques-1];
         cliquetable->ncliques--;
      }
      else
      {
         SCIP_CALL( SCIPhashtableInsert(hashtable, (void*)clique) );
         i++;
         if( !clique->eventsissued )
         {
            int j;

            /* issue IMPLADDED event on each variable in the clique */
            for( j = 0; j < clique->nvars; ++j )
            {
               SCIP_EVENT* event;

               SCIP_CALL( SCIPeventCreateImplAdded(&event, blkmem, clique->vars[j]) );
               SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, NULL, &event) );
            }
            clique->eventsissued = TRUE;
         }
      }
   }

   /* free hash table */
   SCIPhashtableFree(&hashtable);

   /* remember the number of fixed variables and cliques in order to avoid unnecessary cleanups */
   cliquetable->ncleanupfixedvars = stat->npresolfixedvars;
   cliquetable->ncleanupaggrvars = stat->npresolaggrvars;
   cliquetable->ncleanupcliques = cliquetable->ncliques;

   /* process events */
   SCIP_CALL( SCIPeventqueueProcess(eventqueue, blkmem, set, NULL, lp, branchcand, NULL) );
   
   return SCIP_OKAY;
}


/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPvboundsGetNVbds
#undef SCIPvboundsGetVars
#undef SCIPvboundsGetCoefs
#undef SCIPvboundsGetConstants
#undef SCIPimplicsGetNImpls
#undef SCIPimplicsGetNBinImpls
#undef SCIPimplicsGetVars
#undef SCIPimplicsGetTypes
#undef SCIPimplicsGetBounds
#undef SCIPimplicsGetIds
#undef SCIPcliqueGetNVars
#undef SCIPcliqueGetVars
#undef SCIPcliqueGetValues
#undef SCIPcliqueGetId
#undef SCIPcliquelistGetNCliques
#undef SCIPcliquelistGetCliques
#undef SCIPcliquelistCheck
#undef SCIPcliquetableGetNCliques
#undef SCIPcliquetableGetCliques

/** gets number of variable bounds contained in given variable bounds data structure */
int SCIPvboundsGetNVbds(
   SCIP_VBOUNDS*         vbounds             /**< variable bounds data structure */
   )
{
   return vbounds != NULL ? vbounds->len : 0;
}

/** gets array of variables contained in given variable bounds data structure */
SCIP_VAR** SCIPvboundsGetVars(
   SCIP_VBOUNDS*         vbounds             /**< variable bounds data structure */
   )
{
   return vbounds != NULL ? vbounds->vars : NULL;
}

/** gets array of coefficients contained in given variable bounds data structure */
SCIP_Real* SCIPvboundsGetCoefs(
   SCIP_VBOUNDS*         vbounds             /**< variable bounds data structure */
   )
{
   return vbounds != NULL ? vbounds->coefs : NULL;
}

/** gets array of constants contained in given variable bounds data structure */
SCIP_Real* SCIPvboundsGetConstants(
   SCIP_VBOUNDS*         vbounds             /**< variable bounds data structure */
   )
{
   return vbounds != NULL ? vbounds->constants : NULL;
}

/** gets number of implications for a given binary variable fixing */
int SCIPimplicsGetNImpls(
   SCIP_IMPLICS*         implics,            /**< implication data */
   SCIP_Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   )
{
   return implics != NULL ? implics->nimpls[varfixing] : 0;
}

/** gets number of implications on binary variables for a given binary variable fixing */
int SCIPimplicsGetNBinImpls(
   SCIP_IMPLICS*         implics,            /**< implication data */
   SCIP_Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   )
{
   return implics != NULL ? implics->nbinimpls[varfixing] : 0;
}

/** gets array with implied variables for a given binary variable fixing */
SCIP_VAR** SCIPimplicsGetVars(
   SCIP_IMPLICS*         implics,            /**< implication data */
   SCIP_Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   )
{
   return implics != NULL ? implics->vars[varfixing] : NULL;
}

/** gets array with implication types for a given binary variable fixing */
SCIP_BOUNDTYPE* SCIPimplicsGetTypes(
   SCIP_IMPLICS*         implics,            /**< implication data */
   SCIP_Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   )
{
   return implics != NULL ? implics->types[varfixing] : NULL;
}

/** gets array with implication bounds for a given binary variable fixing */
SCIP_Real* SCIPimplicsGetBounds(
   SCIP_IMPLICS*         implics,            /**< implication data */
   SCIP_Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   )
{
   return implics != NULL ? implics->bounds[varfixing] : NULL;
}

/** Gets array with unique implication identifiers for a given binary variable fixing.
 *  If an implication is a shortcut, i.e., it was added as part of the transitive closure of another implication,
 *  its id is negative, otherwise it is nonnegative.
 */
int* SCIPimplicsGetIds(
   SCIP_IMPLICS*         implics,            /**< implication data */
   SCIP_Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   )
{
   return implics != NULL ? implics->ids[varfixing] : NULL;
}

/** gets number of variables in the cliques */
int SCIPcliqueGetNVars(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   )
{
   assert(clique != NULL);

   return clique->nvars;
}

/** gets array of active problem variables in the cliques */
SCIP_VAR** SCIPcliqueGetVars(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   )
{
   assert(clique != NULL);

   return clique->vars;
}

/** gets array of values of active problem variables in the cliques, i.e. whether the variable is fixed to FALSE or
 *  to TRUE in the clique
 */
SCIP_Bool* SCIPcliqueGetValues(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   )
{
   assert(clique != NULL);

   return clique->values;
}

/** gets unique identifier of the clique */
int SCIPcliqueGetId(
   SCIP_CLIQUE*          clique              /**< clique data structure */
   )
{
   assert(clique != NULL);

   return (int) clique->id;
}
   
/** returns the number of cliques stored in the clique list */
int SCIPcliquelistGetNCliques(
   SCIP_CLIQUELIST*      cliquelist,         /**< clique list data structure */
   SCIP_Bool             value               /**< value of the variable for which the cliques should be returned */
   )
{
   return cliquelist != NULL ? cliquelist->ncliques[value] : 0;
}

/** returns the cliques stored in the clique list, or NULL if the clique list is empty */
SCIP_CLIQUE** SCIPcliquelistGetCliques(
   SCIP_CLIQUELIST*      cliquelist,         /**< clique list data structure */
   SCIP_Bool             value               /**< value of the variable for which the cliques should be returned */
   )
{
   return cliquelist != NULL ? cliquelist->cliques[value] : NULL;
}

/** checks whether variable is contained in all cliques of the cliquelist */
void SCIPcliquelistCheck(
   SCIP_CLIQUELIST*      cliquelist,         /**< clique list data structure */
   SCIP_VAR*             var                 /**< variable, the clique list belongs to */
   )
{
#ifndef NDEBUG
   int value;

   assert(SCIPvarGetNCliques(var, FALSE) == SCIPcliquelistGetNCliques(cliquelist, FALSE));
   assert(SCIPvarGetCliques(var, FALSE) == SCIPcliquelistGetCliques(cliquelist, FALSE));
   assert(SCIPvarGetNCliques(var, TRUE) == SCIPcliquelistGetNCliques(cliquelist, TRUE));
   assert(SCIPvarGetCliques(var, TRUE) == SCIPcliquelistGetCliques(cliquelist, TRUE));

   for( value = 0; value < 2; ++value )
   {
      SCIP_CLIQUE** cliques;
      int ncliques;
      int i;
      
      ncliques = SCIPcliquelistGetNCliques(cliquelist, (SCIP_Bool)value);
      cliques = SCIPcliquelistGetCliques(cliquelist, (SCIP_Bool)value);
      for( i = 0; i < ncliques; ++i )
      {

         SCIP_CLIQUE* clique;
         int pos;

         clique = cliques[i];
         assert(clique != NULL);

         pos = SCIPcliqueSearchVar(clique, var, (SCIP_Bool)value);
         assert(0 <= pos && pos < clique->nvars);
         assert(clique->vars[pos] == var);
         assert(clique->values[pos] == (SCIP_Bool)value);
      }
   }
#endif
}

/** gets the number of cliques stored in the clique table */
int SCIPcliquetableGetNCliques(
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   assert(cliquetable != NULL);

   return cliquetable->ncliques;
}

/** gets the array of cliques stored in the clique table */
SCIP_CLIQUE** SCIPcliquetableGetCliques(
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   assert(cliquetable != NULL);

   return cliquetable->cliques;
}

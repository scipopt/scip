/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: misc.h,v 1.17 2005/01/18 09:26:48 bzfpfend Exp $"

/**@file   misc.h
 * @brief  internal miscellaneous methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __MISC_H__
#define __MISC_H__


#include "def.h"
#include "memory.h"
#include "type_retcode.h"
#include "type_set.h"
#include "type_misc.h"
#include "pub_misc.h"



/*
 * Dynamic Arrays
 */

/** creates a dynamic array of real values */
extern
RETCODE SCIPrealarrayCreate(
   REALARRAY**      realarray,          /**< pointer to store the real array */
   MEMHDR*          memhdr              /**< block memory */
   );

/** creates a copy of a dynamic array of real values */
extern
RETCODE SCIPrealarrayCopy(
   REALARRAY**      realarray,          /**< pointer to store the copied real array */
   MEMHDR*          memhdr,             /**< block memory */
   REALARRAY*       sourcerealarray     /**< dynamic real array to copy */
   );

/** frees a dynamic array of real values */
extern
RETCODE SCIPrealarrayFree(
   REALARRAY**      realarray           /**< pointer to the real array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
RETCODE SCIPrealarrayExtend(
   REALARRAY*       realarray,          /**< dynamic real array */
   SET*             set,                /**< global SCIP settings */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic real array */
extern
RETCODE SCIPrealarrayClear(
   REALARRAY*       realarray           /**< dynamic real array */
   );

/** gets value of entry in dynamic array */
extern
Real SCIPrealarrayGetVal(
   REALARRAY*       realarray,          /**< dynamic real array */
   int              idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
extern
RETCODE SCIPrealarraySetVal(
   REALARRAY*       realarray,          /**< dynamic real array */
   SET*             set,                /**< global SCIP settings */
   int              idx,                /**< array index to set value for */
   Real             val                 /**< value to set array index to */
   );

/** increases value of entry in dynamic array */
extern
RETCODE SCIPrealarrayIncVal(
   REALARRAY*       realarray,          /**< dynamic real array */
   SET*             set,                /**< global SCIP settings */
   int              idx,                /**< array index to increase value for */
   Real             incval              /**< value to increase array index */
   );

/** returns the minimal index of all stored non-zero elements */
extern
int SCIPrealarrayGetMinIdx(
   REALARRAY*       realarray           /**< dynamic real array */
   );

/** returns the maximal index of all stored non-zero elements */
extern
int SCIPrealarrayGetMaxIdx(
   REALARRAY*       realarray           /**< dynamic real array */
   );

/** creates a dynamic array of int values */
extern
RETCODE SCIPintarrayCreate(
   INTARRAY**       intarray,           /**< pointer to store the int array */
   MEMHDR*          memhdr              /**< block memory */
   );

/** creates a copy of a dynamic array of int values */
extern
RETCODE SCIPintarrayCopy(
   INTARRAY**       intarray,           /**< pointer to store the copied int array */
   MEMHDR*          memhdr,             /**< block memory */
   INTARRAY*        sourceintarray      /**< dynamic int array to copy */
   );

/** frees a dynamic array of int values */
extern
RETCODE SCIPintarrayFree(
   INTARRAY**       intarray            /**< pointer to the int array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
RETCODE SCIPintarrayExtend(
   INTARRAY*        intarray,           /**< dynamic int array */
   SET*             set,                /**< global SCIP settings */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic int array */
extern
RETCODE SCIPintarrayClear(
   INTARRAY*        intarray            /**< dynamic int array */
   );

/** gets value of entry in dynamic array */
extern
int SCIPintarrayGetVal(
   INTARRAY*        intarray,           /**< dynamic int array */
   int              idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
extern
RETCODE SCIPintarraySetVal(
   INTARRAY*        intarray,           /**< dynamic int array */
   SET*             set,                /**< global SCIP settings */
   int              idx,                /**< array index to set value for */
   int              val                 /**< value to set array index to */
   );

/** increases value of entry in dynamic array */
extern
RETCODE SCIPintarrayIncVal(
   INTARRAY*        intarray,           /**< dynamic int array */
   SET*             set,                /**< global SCIP settings */
   int              idx,                /**< array index to increase value for */
   int              incval              /**< value to increase array index */
   );

/** returns the minimal index of all stored non-zero elements */
extern
int SCIPintarrayGetMinIdx(
   INTARRAY*        intarray            /**< dynamic int array */
   );

/** returns the maximal index of all stored non-zero elements */
extern
int SCIPintarrayGetMaxIdx(
   INTARRAY*        intarray            /**< dynamic int array */
   );

/** creates a dynamic array of bool values */
extern
RETCODE SCIPboolarrayCreate(
   BOOLARRAY**      boolarray,          /**< pointer to store the bool array */
   MEMHDR*          memhdr              /**< block memory */
   );

/** creates a copy of a dynamic array of bool values */
extern
RETCODE SCIPboolarrayCopy(
   BOOLARRAY**      boolarray,          /**< pointer to store the copied bool array */
   MEMHDR*          memhdr,             /**< block memory */
   BOOLARRAY*       sourceboolarray     /**< dynamic bool array to copy */
   );

/** frees a dynamic array of bool values */
extern
RETCODE SCIPboolarrayFree(
   BOOLARRAY**      boolarray           /**< pointer to the bool array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
RETCODE SCIPboolarrayExtend(
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   SET*             set,                /**< global SCIP settings */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic bool array */
extern
RETCODE SCIPboolarrayClear(
   BOOLARRAY*       boolarray           /**< dynamic bool array */
   );

/** gets value of entry in dynamic array */
extern
Bool SCIPboolarrayGetVal(
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int              idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
extern
RETCODE SCIPboolarraySetVal(
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   SET*             set,                /**< global SCIP settings */
   int              idx,                /**< array index to set value for */
   Bool             val                 /**< value to set array index to */
   );

/** returns the minimal index of all stored non-zero elements */
extern
int SCIPboolarrayGetMinIdx(
   BOOLARRAY*       boolarray           /**< dynamic bool array */
   );

/** returns the maximal index of all stored non-zero elements */
extern
int SCIPboolarrayGetMaxIdx(
   BOOLARRAY*       boolarray           /**< dynamic bool array */
   );

/** creates a dynamic array of pointer values */
extern
RETCODE SCIPptrarrayCreate(
   PTRARRAY**       ptrarray,           /**< pointer to store the ptr array */
   MEMHDR*          memhdr              /**< block memory */
   );

/** creates a copy of a dynamic array of pointer values */
extern
RETCODE SCIPptrarrayCopy(
   PTRARRAY**       ptrarray,           /**< pointer to store the copied ptr array */
   MEMHDR*          memhdr,             /**< block memory */
   PTRARRAY*        sourceptrarray      /**< dynamic ptr array to copy */
   );

/** frees a dynamic array of pointer values */
extern
RETCODE SCIPptrarrayFree(
   PTRARRAY**       ptrarray            /**< pointer to the ptr array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
RETCODE SCIPptrarrayExtend(
   PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   SET*             set,                /**< global SCIP settings */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic pointer array */
extern
RETCODE SCIPptrarrayClear(
   PTRARRAY*        ptrarray            /**< dynamic ptr array */
   );

/** gets value of entry in dynamic array */
extern
void* SCIPptrarrayGetVal(
   PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   int              idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
extern
RETCODE SCIPptrarraySetVal(
   PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   SET*             set,                /**< global SCIP settings */
   int              idx,                /**< array index to set value for */
   void*            val                 /**< value to set array index to */
   );

/** returns the minimal index of all stored non-zero elements */
extern
int SCIPptrarrayGetMinIdx(
   PTRARRAY*        ptrarray            /**< dynamic ptr array */
   );

/** returns the maximal index of all stored non-zero elements */
extern
int SCIPptrarrayGetMaxIdx(
   PTRARRAY*        ptrarray            /**< dynamic ptr array */
   );


#endif

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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_paramset.h,v 1.6 2005/01/21 09:17:08 bzfpfend Exp $"

/**@file   struct_paramset.h
 * @brief  datastructures for handling parameter settings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_PARAMSET_H__
#define __STRUCT_PARAMSET_H__


#include "def.h"
#include "type_misc.h"



/** data for Bool parameters */
struct BoolParam
{
   Bool*            valueptr;           /**< pointer to store the current parameter value, or NULL */
   Bool             curvalue;           /**< stores the current parameter value if it is not stored in *valueptr */
   Bool             defaultvalue;       /**< default value of the parameter */
};
typedef struct BoolParam BOOLPARAM;

/** data for int parameters */
struct IntParam
{
   int*             valueptr;           /**< pointer to store the current parameter value, or NULL */
   int              curvalue;           /**< stores the current parameter value if it is not stored in *valueptr */
   int              defaultvalue;       /**< default value of the parameter */
   int              minvalue;           /**< minimum value for parameter */
   int              maxvalue;           /**< maximum value for parameter */
};
typedef struct IntParam INTPARAM;

/** data for Longint parameters */
struct LongintParam
{
   Longint          curvalue;           /**< stores the current parameter value if it is not stored in *valueptr */
   Longint          defaultvalue;       /**< default value of the parameter */
   Longint          minvalue;           /**< minimum value for parameter */
   Longint          maxvalue;           /**< maximum value for parameter */
   Longint*         valueptr;           /**< pointer to store the current parameter value, or NULL */
};
typedef struct LongintParam LONGINTPARAM;

/** data for Real parameters */
struct RealParam
{
   Real             curvalue;           /**< stores the current parameter value if it is not stored in *valueptr */
   Real             defaultvalue;       /**< default value of the parameter */
   Real             minvalue;           /**< minimum value for parameter */
   Real             maxvalue;           /**< maximum value for parameter */
   Real*            valueptr;           /**< pointer to store the current parameter value, or NULL */
};
typedef struct RealParam REALPARAM;

/** data for char parameters */
struct CharParam
{
   char*            valueptr;           /**< pointer to store the current parameter value, or NULL */
   char*            allowedvalues;      /**< array with possible parameter values, or NULL if not restricted */
   char             curvalue;           /**< stores the current parameter value if it is not stored in *valueptr */
   char             defaultvalue;       /**< default value of the parameter */
};
typedef struct CharParam CHARPARAM;

/** data for char* parameters */
struct StringParam
{
   char**           valueptr;           /**< pointer to store the current parameter value, or NULL */
   char*            curvalue;           /**< stores the current parameter value if it is not stored in *valueptr */
   char*            defaultvalue;       /**< default value of the parameter */
};
typedef struct StringParam STRINGPARAM;

/** single parameter */
struct Param
{
   union
   {
      BOOLPARAM     boolparam;          /**< data for Bool parameters */
      INTPARAM      intparam;           /**< data for int parameters */
      LONGINTPARAM  longintparam;       /**< data for Longint parameters */
      REALPARAM     realparam;          /**< data for Real parameters */
      CHARPARAM     charparam;          /**< data for char parameters */
      STRINGPARAM   stringparam;        /**< data for char* parameters */
   } data;
   char*            name;               /**< name of the parameter */
   char*            desc;               /**< description of the parameter */
   DECL_PARAMCHGD   ((*paramchgd));     /**< change information method of parameter */
   PARAMDATA*       paramdata;          /**< locally defined parameter specific data */
   PARAMTYPE        paramtype;          /**< type of this parameter */
};

/** set of parameters */
struct ParamSet
{
   HASHTABLE*       hashtable;          /**< hash table to store the parameters */
   PARAM**          params;             /**< array with parameters */
   int              nparams;            /**< number of parameters */
   int              paramssize;         /**< size of params array */
};


#endif

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   paramset.h
 * @brief  methods and datastructures for handling parameter settings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PARAMSET_H__
#define __PARAMSET_H__


/** possible parameter types */
enum ParamType
{
   SCIP_PARAMTYPE_BOOL    = 0,           /**< bool values: TRUE or FALSE */
   SCIP_PARAMTYPE_INT     = 1,           /**< integer values */
   SCIP_PARAMTYPE_LONGINT = 2,           /**< long integer values */
   SCIP_PARAMTYPE_REAL    = 3,           /**< real values */
   SCIP_PARAMTYPE_CHAR    = 4,           /**< characters */
   SCIP_PARAMTYPE_STRING  = 5            /**< strings: arrays of characters */
};
typedef enum ParamType PARAMTYPE;

typedef struct Param PARAM;             /**< single parameter */

typedef struct ParamSet PARAMSET;       /**< set of parameters */



#include <math.h>

#include "def.h"
#include "retcode.h"



/** creates parameter set */
extern
RETCODE SCIPparamsetCreate(
   PARAMSET**       paramset            /**< pointer to store the parameter set */
   );

/** frees parameter set */
extern
void SCIPparamsetFree(
   PARAMSET**       paramset,           /**< pointer to the parameter set */
   MEMHDR*          memhdr              /**< block memory */
   );

/** creates a bool parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPparamsetAddBool(
   PARAMSET*        paramset,           /**< parameter set */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   Bool*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   Bool             defaultvalue        /**< default value of the parameter */
   );

/** creates a int parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPparamsetAddInt(
   PARAMSET*        paramset,           /**< parameter set */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   int*             valueptr,           /**< pointer to store the current parameter value, or NULL */
   int              defaultvalue        /**< default value of the parameter */
   );

/** creates a Longint parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPparamsetAddLongint(
   PARAMSET*        paramset,           /**< parameter set */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   Longint*         valueptr,           /**< pointer to store the current parameter value, or NULL */
   Longint          defaultvalue        /**< default value of the parameter */
   );

/** creates a Real parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPparamsetAddReal(
   PARAMSET*        paramset,           /**< parameter set */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   Real*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   Real             defaultvalue        /**< default value of the parameter */
   );

/** creates a char parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPparamsetAddChar(
   PARAMSET*        paramset,           /**< parameter set */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   char*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   char             defaultvalue        /**< default value of the parameter */
   );

/** creates a string parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPparamsetAddString(
   PARAMSET*        paramset,           /**< parameter set */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   char**           valueptr,           /**< pointer to store the current parameter value, or NULL */
   const char*      defaultvalue        /**< default value of the parameter */
   );

/** gets the value of an existing Bool parameter */
extern
RETCODE SCIPparamsetGetBool(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   Bool*            value               /**< pointer to store the parameter */
   );

/** gets the value of an existing int parameter */
extern
RETCODE SCIPparamsetGetInt(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   int*             value               /**< pointer to store the parameter */
   );

/** gets the value of an existing Longint parameter */
extern
RETCODE SCIPparamsetGetLongint(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   Longint*         value               /**< pointer to store the parameter */
   );

/** gets the value of an existing Real parameter */
extern
RETCODE SCIPparamsetGetReal(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   Real*            value               /**< pointer to store the parameter */
   );

/** gets the value of an existing char parameter */
extern
RETCODE SCIPparamsetGetChar(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   char*            value               /**< pointer to store the parameter */
   );

/** gets the value of an existing string parameter */
extern
RETCODE SCIPparamsetGetString(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   char**           value               /**< pointer to store the parameter */
   );

/** changes the value of an existing Bool parameter */
extern
RETCODE SCIPparamsetSetBool(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   Bool             value               /**< new value of the parameter */
   );

/** changes the value of an existing int parameter */
extern
RETCODE SCIPparamsetSetInt(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   int              value               /**< new value of the parameter */
   );

/** changes the value of an existing Longint parameter */
extern
RETCODE SCIPparamsetSetLongint(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   Longint          value               /**< new value of the parameter */
   );

/** changes the value of an existing Real parameter */
extern
RETCODE SCIPparamsetSetReal(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   Real             value               /**< new value of the parameter */
   );

/** changes the value of an existing char parameter */
extern
RETCODE SCIPparamsetSetChar(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   char             value               /**< new value of the parameter */
   );

/** changes the value of an existing string parameter */
extern
RETCODE SCIPparamsetSetString(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   const char*      value               /**< new value of the parameter */
   );

#endif

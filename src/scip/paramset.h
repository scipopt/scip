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
typedef struct ParamData PARAMDATA;     /**< locally defined parameter specific data */
typedef struct ParamSet PARAMSET;       /**< set of parameters */


/** information method for changes in the parameter
 *
 *  Method is called if the parameter was changed through a SCIPparamsetSetXxx() call
 *  (which is called by SCIPsetXxxParam()).
 *  It will not be called, if the parameter was changed directly by changing the value
 *  in the memory location.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    param           : the changed parameter (already set to its new value)
 */
#define DECL_PARAMCHGD(x) RETCODE x (SCIP* scip, PARAM* param)



#include <math.h>

#include "def.h"
#include "retcode.h"



/*
 * Parameter methods
 */

/** returns type of parameter */
extern
PARAMTYPE SCIPparamGetType(
   PARAM*           param               /**< parameter */
   );

/** returns name of parameter */
extern
const char* SCIPparamGetName(
   PARAM*           param               /**< parameter */
   );

/** returns description of parameter */
extern
const char* SCIPparamGetDesc(
   PARAM*           param               /**< parameter */
   );

/** returns locally defined parameter specific data */
extern
PARAMDATA* SCIPparamGetData(
   PARAM*           param               /**< parameter */
   );

/** returns value of Bool parameter */
extern
Bool SCIPparamGetBool(
   PARAM*           param               /**< parameter */
   );

/** returns value of int parameter */
extern
int SCIPparamGetInt(
   PARAM*           param               /**< parameter */
   );

/** returns minimal value of int parameter */
extern
int SCIPparamGetIntMin(
   PARAM*           param               /**< parameter */
   );

/** returns maximal value of int parameter */
extern
int SCIPparamGetIntMax(
   PARAM*           param               /**< parameter */
   );

/** returns value of Longint parameter */
extern
Longint SCIPparamGetLongint(
   PARAM*           param               /**< parameter */
   );

/** returns minimal value of longint parameter */
extern
Longint SCIPparamGetLongintMin(
   PARAM*           param               /**< parameter */
   );

/** returns maximal value of longint parameter */
extern
Longint SCIPparamGetLongintMax(
   PARAM*           param               /**< parameter */
   );

/** returns value of Real parameter */
extern
Real SCIPparamGetReal(
   PARAM*           param               /**< parameter */
   );

/** returns minimal value of real parameter */
extern
Real SCIPparamGetRealMin(
   PARAM*           param               /**< parameter */
   );

/** returns maximal value of real parameter */
extern
Real SCIPparamGetRealMax(
   PARAM*           param               /**< parameter */
   );

/** returns value of char parameter */
extern
char SCIPparamGetChar(
   PARAM*           param               /**< parameter */
   );

/** returns value of string parameter */
extern
char* SCIPparamGetString(
   PARAM*           param               /**< parameter */
   );

/** sets value of Bool parameter */
RETCODE SCIPparamSetBool(
   PARAM*           param,              /**< parameter */
   SCIP*            scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   Bool             value               /**< new value of the parameter */
   );

/** sets value of int parameter */
extern
RETCODE SCIPparamSetInt(
   PARAM*           param,              /**< parameter */
   SCIP*            scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   int              value               /**< new value of the parameter */
   );

/** sets value of Longint parameter */
extern
RETCODE SCIPparamSetLongint(
   PARAM*           param,              /**< parameter */
   SCIP*            scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   Longint          value               /**< new value of the parameter */
   );

/** sets value of Real parameter */
extern
RETCODE SCIPparamSetReal(
   PARAM*           param,              /**< parameter */
   SCIP*            scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   Real             value               /**< new value of the parameter */
   );

/** sets value of char parameter */
extern
RETCODE SCIPparamSetChar(
   PARAM*           param,              /**< parameter */
   SCIP*            scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   char             value               /**< new value of the parameter */
   );

/** sets value of string parameter */
extern
RETCODE SCIPparamSetString(
   PARAM*           param,              /**< parameter */
   SCIP*            scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   const char*      value               /**< new value of the parameter */
   );




/*
 * Parameter set methods
 */

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
   const char*      desc,               /**< description of the parameter */
   Bool*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   Bool             defaultvalue,       /**< default value of the parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a int parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPparamsetAddInt(
   PARAMSET*        paramset,           /**< parameter set */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   int*             valueptr,           /**< pointer to store the current parameter value, or NULL */
   int              defaultvalue,       /**< default value of the parameter */
   int              minvalue,           /**< minimum value for parameter */
   int              maxvalue,           /**< maximum value for parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a Longint parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPparamsetAddLongint(
   PARAMSET*        paramset,           /**< parameter set */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   Longint*         valueptr,           /**< pointer to store the current parameter value, or NULL */
   Longint          defaultvalue,       /**< default value of the parameter */
   Longint          minvalue,           /**< minimum value for parameter */
   Longint          maxvalue,           /**< maximum value for parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a Real parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPparamsetAddReal(
   PARAMSET*        paramset,           /**< parameter set */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   Real*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   Real             defaultvalue,       /**< default value of the parameter */
   Real             minvalue,           /**< minimum value for parameter */
   Real             maxvalue,           /**< maximum value for parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a char parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPparamsetAddChar(
   PARAMSET*        paramset,           /**< parameter set */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   char*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   char             defaultvalue,       /**< default value of the parameter */
   const char*      allowedvalues,      /**< array with possible parameter values, or NULL if not restricted */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a string parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPparamsetAddString(
   PARAMSET*        paramset,           /**< parameter set */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   char**           valueptr,           /**< pointer to store the current parameter value, or NULL */
   const char*      defaultvalue,       /**< default value of the parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
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
   SCIP*            scip,               /**< SCIP data structure */   
   const char*      name,               /**< name of the parameter */
   Bool             value               /**< new value of the parameter */
   );

/** changes the value of an existing int parameter */
extern
RETCODE SCIPparamsetSetInt(
   PARAMSET*        paramset,           /**< parameter set */
   SCIP*            scip,               /**< SCIP data structure */   
   const char*      name,               /**< name of the parameter */
   int              value               /**< new value of the parameter */
   );

/** changes the value of an existing Longint parameter */
extern
RETCODE SCIPparamsetSetLongint(
   PARAMSET*        paramset,           /**< parameter set */
   SCIP*            scip,               /**< SCIP data structure */   
   const char*      name,               /**< name of the parameter */
   Longint          value               /**< new value of the parameter */
   );

/** changes the value of an existing Real parameter */
extern
RETCODE SCIPparamsetSetReal(
   PARAMSET*        paramset,           /**< parameter set */
   SCIP*            scip,               /**< SCIP data structure */   
   const char*      name,               /**< name of the parameter */
   Real             value               /**< new value of the parameter */
   );

/** changes the value of an existing char parameter */
extern
RETCODE SCIPparamsetSetChar(
   PARAMSET*        paramset,           /**< parameter set */
   SCIP*            scip,               /**< SCIP data structure */   
   const char*      name,               /**< name of the parameter */
   char             value               /**< new value of the parameter */
   );

/** changes the value of an existing string parameter */
extern
RETCODE SCIPparamsetSetString(
   PARAMSET*        paramset,           /**< parameter set */
   SCIP*            scip,               /**< SCIP data structure */   
   const char*      name,               /**< name of the parameter */
   const char*      value               /**< new value of the parameter */
   );

/** reads parameters from a file */
RETCODE SCIPparamsetRead(
   PARAMSET*        paramset,           /**< parameter set */
   SCIP*            scip,               /**< SCIP data structure */   
   const char*      filename            /**< file name */
   );

/** writes all parameters in the parameter set to a file */
RETCODE SCIPparamsetWrite(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      filename,           /**< file name, or NULL for stdout */
   Bool             comments            /**< should parameter descriptions be written as comments? */
   );

/** returns the array of parameters */
extern
PARAM** SCIPparamsetGetParams(
   PARAMSET*        paramset            /**< parameter set */
   );

/** returns the number of parameters in the parameter set */
extern
int SCIPparamsetGetNParams(
   PARAMSET*        paramset            /**< parameter set */
   );

#endif

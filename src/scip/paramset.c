/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   paramset.c
 * @brief  methods for handling parameter settings
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/set.h"
#include "scip/paramset.h"

#include "scip/struct_paramset.h"



/*
 * Parameter methods
 */

/** hash key retrieval function for parameters */
static
SCIP_DECL_HASHGETKEY(hashGetKeyParam)
{  /*lint --e{715}*/
   SCIP_PARAM* param;

   param = (SCIP_PARAM*)elem;
   assert(param != NULL);

   return param->name;
}

/** checks parameter value according to the given feasible domain; issues a warning message if value was invalid */
static
SCIP_RETCODE paramCheckBool(
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_Bool             value               /**< value to check */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_BOOL);

   if( value != TRUE && value != FALSE )
   {
      SCIPwarningMessage("Invalid value <%d> for bool parameter <%s>. Must be <0> (FALSE) or <1> (TRUE).\n",
         value, param->name);
      return SCIP_PARAMETERWRONGVAL;
   }

   return SCIP_OKAY;
}

/** checks parameter value according to the given feasible domain; issues a warning message if value was invalid */
static
SCIP_RETCODE paramCheckInt(
   SCIP_PARAM*           param,              /**< parameter */
   int                   value               /**< value to check */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_INT);

   if( value < param->data.intparam.minvalue || value > param->data.intparam.maxvalue )
   {
      SCIPwarningMessage("Invalid value <%d> for int parameter <%s>. Must be in range [%d,%d].\n",
         value, param->name, param->data.intparam.minvalue, param->data.intparam.maxvalue);
      return SCIP_PARAMETERWRONGVAL;
   }
   
   return SCIP_OKAY;
}

/** checks parameter value according to the given feasible domain; issues a warning message if value was invalid */
static
SCIP_RETCODE paramCheckLongint(
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_Longint          value               /**< value to check */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_LONGINT);

   if( value < param->data.longintparam.minvalue || value > param->data.longintparam.maxvalue )
   {
      SCIPwarningMessage("Invalid value <%"SCIP_LONGINT_FORMAT"> for longint parameter <%s>. Must be in range [%"SCIP_LONGINT_FORMAT",%"SCIP_LONGINT_FORMAT"].\n",
         value, param->name, param->data.longintparam.minvalue, param->data.longintparam.maxvalue);
      return SCIP_PARAMETERWRONGVAL;
   }
   
   return SCIP_OKAY;
}

/** checks parameter value according to the given feasible domain; issues a warning message if value was invalid */
static
SCIP_RETCODE paramCheckReal(
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_Real             value               /**< value to check */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_REAL);

   if( value < param->data.realparam.minvalue || value > param->data.realparam.maxvalue )
   {
      SCIPwarningMessage("Invalid real parameter value <%.15g> for parameter <%s>. Must be in range [%.15g,%.15g].\n",
         value, param->name, param->data.realparam.minvalue, param->data.realparam.maxvalue);
      return SCIP_PARAMETERWRONGVAL;
   }

   return SCIP_OKAY;
}

/** checks parameter value according to the given feasible domain; issues a warning message if value was invalid */
static
SCIP_RETCODE paramCheckChar(
   SCIP_PARAM*           param,              /**< parameter */
   char                  value               /**< value to check */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_CHAR);

   if( value == '\b' || value == '\f' || value == '\n' || value == '\r' || value == '\v' )
   {
      SCIPwarningMessage("Invalid char parameter value <%x>.\n", (int)value);
      return SCIP_PARAMETERWRONGVAL;
   }

   if( param->data.charparam.allowedvalues != NULL )
   {
      char* c;

      c = param->data.charparam.allowedvalues;
      while( *c != '\0' && *c != value )
         c++;

      if( *c != value )
      {
         SCIPwarningMessage("Invalid char parameter value <%c>. Must be in set {%s}.\n",
            value, param->data.charparam.allowedvalues);
         return SCIP_PARAMETERWRONGVAL;
      }
   }

   return SCIP_OKAY;
}

/** checks parameter value according to the given feasible domain; issues a warning message if value was invalid */
static
SCIP_RETCODE paramCheckString(
   SCIP_PARAM*           param,              /**< parameter */
   const char*           value               /**< value to check */
   )
{
   unsigned int i;

   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_STRING);

   if( value == NULL )
   {
      SCIPwarningMessage("Cannot assign a NULL string to a string parameter.\n");
      return SCIP_PARAMETERWRONGVAL;
   }

   for( i = 0; i < strlen(value); ++i )
   {
      if( value[i] == '\b' || value[i] == '\f' || value[i] == '\n' || value[i] == '\r' || value[i] == '\v' )
      {
         SCIPwarningMessage("Invalid character <%x> in string parameter at position %d.\n", (int)value[i], i);
         return SCIP_PARAMETERWRONGVAL;
      }
   }

   return SCIP_OKAY;
}

/** writes the parameter to a file */
static
SCIP_RETCODE paramWrite(
   SCIP_PARAM*           param,              /**< parameter */
   FILE*                 file,               /**< file stream to write parameter to, or NULL for stdout  */
   SCIP_Bool             comments,           /**< should parameter descriptions be written as comments? */
   SCIP_Bool             onlychanged         /**< should only the parameters been written, that are changed from default? */
   )
{
   assert(param != NULL);

   /* write parameters at default values only, if the onlychanged flag is not set */
   if( onlychanged && SCIPparamIsDefault(param) )
      return SCIP_OKAY;

   /* write parameter description, bounds, and defaults as comments */
   if( comments )
   {
      SCIPmessageFPrintInfo(file, "# %s\n", param->desc);
      switch( param->paramtype )
      {
      case SCIP_PARAMTYPE_BOOL:
         SCIPmessageFPrintInfo(file, "# [type: bool, range: {TRUE,FALSE}, default: %s]\n",
            param->data.boolparam.defaultvalue ? "TRUE" : "FALSE");
         break;
      case SCIP_PARAMTYPE_INT:
         SCIPmessageFPrintInfo(file, "# [type: int, range: [%d,%d], default: %d]\n", 
            param->data.intparam.minvalue, param->data.intparam.maxvalue, param->data.intparam.defaultvalue);
         break;
      case SCIP_PARAMTYPE_LONGINT:
         SCIPmessageFPrintInfo(file, "# [type: longint, range: [%"SCIP_LONGINT_FORMAT",%"SCIP_LONGINT_FORMAT"], default: %"SCIP_LONGINT_FORMAT"]\n", 
            param->data.longintparam.minvalue, param->data.longintparam.maxvalue, param->data.longintparam.defaultvalue);
         break;
      case SCIP_PARAMTYPE_REAL:
         SCIPmessageFPrintInfo(file, "# [type: real, range: [%.15g,%.15g], default: %.15g]\n",
            param->data.realparam.minvalue, param->data.realparam.maxvalue, param->data.realparam.defaultvalue);
         break;
      case SCIP_PARAMTYPE_CHAR:
         SCIPmessageFPrintInfo(file, "# [type: char, range: {%s}, default: %c]\n",
            param->data.charparam.allowedvalues != NULL ? param->data.charparam.allowedvalues : "all chars",
            param->data.charparam.defaultvalue);
         break;
      case SCIP_PARAMTYPE_STRING:
         SCIPmessageFPrintInfo(file, "# [type: string, default: \"%s\"]\n", param->data.stringparam.defaultvalue);
         break;
      default:
         SCIPerrorMessage("unknown parameter type\n");
         return SCIP_INVALIDDATA;
      }
   }

   /* write parameter value */
   SCIPmessageFPrintInfo(file, "%s = ", param->name);
   switch( param->paramtype )
   {
   case SCIP_PARAMTYPE_BOOL:
      SCIPmessageFPrintInfo(file, "%s\n", SCIPparamGetBool(param) ? "TRUE" : "FALSE");
      break;
   case SCIP_PARAMTYPE_INT:
      SCIPmessageFPrintInfo(file, "%d\n", SCIPparamGetInt(param));
      break;
   case SCIP_PARAMTYPE_LONGINT:
      SCIPmessageFPrintInfo(file, "%"SCIP_LONGINT_FORMAT"\n", SCIPparamGetLongint(param));
      break;
   case SCIP_PARAMTYPE_REAL:
      SCIPmessageFPrintInfo(file, "%.15g\n", SCIPparamGetReal(param));
      break;
   case SCIP_PARAMTYPE_CHAR:
      SCIPmessageFPrintInfo(file, "%c\n", SCIPparamGetChar(param));
      break;
   case SCIP_PARAMTYPE_STRING:
      SCIPmessageFPrintInfo(file, "\"%s\"\n", SCIPparamGetString(param));
      break;
   default:
      SCIPerrorMessage("unknown parameter type\n");
      return SCIP_INVALIDDATA;
   }

   if( comments )
      SCIPmessageFPrintInfo(file, "\n");

   return SCIP_OKAY;
}

/** if a bool parameter exits with the given parameter name it is set to the new value */
static
SCIP_RETCODE paramSetBool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   const char*           paramname,          /**< parameter name */
   SCIP_Bool             value,              /**< new value of the parameter */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_PARAM* param;

   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)paramname);
   if( param != NULL )
   {
      assert(SCIPparamGetType(param) == SCIP_PARAMTYPE_BOOL);
      SCIP_CALL( SCIPparamSetBool(param, scip, value, quiet) );
   }
#ifndef NDEBUG
   else
   {
      SCIPwarningMessage("unknown hard coded bool parameter <%s>\n", paramname);
   }
#endif
   
   return SCIP_OKAY;
}

/** if an integer parameter exits with the given parameter name it is set to the new value */
static
SCIP_RETCODE paramSetInt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   const char*           paramname,          /**< parameter name */
   int                   value,              /**< new value of the parameter */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_PARAM* param;

   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)paramname);
   if( param != NULL )
   {
      assert(SCIPparamGetType(param) == SCIP_PARAMTYPE_INT);
      SCIP_CALL( SCIPparamSetInt(param, scip, value, quiet) );
   }
#ifndef NDEBUG
   else
   {
      SCIPwarningMessage("unknown hard coded int parameter <%s>\n", paramname);
   }
#endif

   return SCIP_OKAY;
}

/** if a long integer parameter exits with the given parameter name it is set to the new value */
static
SCIP_RETCODE paramSetLongint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   const char*           paramname,          /**< parameter name */
   SCIP_Longint          value,              /**< new value of the parameter */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_PARAM* param;

   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)paramname);
   if( param != NULL )
   {
      assert(SCIPparamGetType(param) == SCIP_PARAMTYPE_LONGINT);
      SCIP_CALL( SCIPparamSetLongint(param, scip, value, quiet) );
   }
#ifndef NDEBUG
   else
   {
      SCIPwarningMessage("unknown hard coded Longint parameter <%s>\n", paramname);
   }
#endif

   return SCIP_OKAY;
}

/** if a real parameter exits with the given parameter name it is set to the new value */
static
SCIP_RETCODE paramSetReal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   const char*           paramname,          /**< parameter name */
   SCIP_Real             value,              /**< new value of the parameter */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_PARAM* param;

   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)paramname);
   if( param != NULL )
   {
      assert(SCIPparamGetType(param) == SCIP_PARAMTYPE_REAL);
      SCIP_CALL( SCIPparamSetReal(param, scip, value, quiet) );
   }
#ifndef NDEBUG
   else
   {
      SCIPwarningMessage("unknown hard coded real parameter <%s>\n", paramname);
   }
#endif
   
   return SCIP_OKAY;
}

/** copies value of source Bool parameter to target Bool parameter*/
static
SCIP_RETCODE paramCopyBool(
   SCIP_PARAM*           sourceparam,        /**< source Bool parameter */
   SCIP_PARAM*           targetparam,        /**< target Bool parameter */
   SCIP*                 targetscip          /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   )
{
   SCIP_Bool value;

   assert(sourceparam != NULL);
   assert(targetparam != NULL);

   /* get value of source parameter and copy it to target parameter */
   value = SCIPparamGetBool(sourceparam);
   SCIP_CALL( SCIPparamSetBool(targetparam, targetscip, value, TRUE) );

   return SCIP_OKAY;
}

/** copies value of source int parameter to target int parameter*/
static
SCIP_RETCODE paramCopyInt(
   SCIP_PARAM*           sourceparam,        /**< source int parameter */
   SCIP_PARAM*           targetparam,        /**< target int parameter */
   SCIP*                 targetscip          /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   )
{
   int value;

   assert(sourceparam != NULL);
   assert(targetparam != NULL);

   /* get value of source parameter and copy it to target parameter */
   value = SCIPparamGetInt(sourceparam);
   SCIP_CALL( SCIPparamSetInt(targetparam,targetscip, value, TRUE) );

   return SCIP_OKAY;
}

/** copies value of source longint parameter to target longint parameter*/
static
SCIP_RETCODE paramCopyLongint(
   SCIP_PARAM*           sourceparam,        /**< source longint parameter */
   SCIP_PARAM*           targetparam,        /**< target longint parameter */
   SCIP*                 targetscip          /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   )
{
   SCIP_Longint value;

   assert(sourceparam != NULL);
   assert(targetparam != NULL);

   /* get value of source parameter and copy it to target parameter */
   value = SCIPparamGetLongint(sourceparam);
   SCIP_CALL( SCIPparamSetLongint(targetparam,targetscip, value, TRUE) );

   return SCIP_OKAY;
}

/** copies value of source real parameter to target real parameter*/
static
SCIP_RETCODE paramCopyReal(
   SCIP_PARAM*           sourceparam,        /**< source real parameter */
   SCIP_PARAM*           targetparam,        /**< target real parameter */
   SCIP*                 targetscip          /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   )
{
   SCIP_Real value;

   assert(sourceparam != NULL);
   assert(targetparam != NULL);

   /* get value of source parameter and copy it to target parameter */
   value = SCIPparamGetReal(sourceparam);
   SCIP_CALL( SCIPparamSetReal(targetparam,targetscip, value, TRUE) );

   return SCIP_OKAY;
}

/** copies value of source char parameter to target char parameter*/
static
SCIP_RETCODE paramCopyChar(
   SCIP_PARAM*           sourceparam,        /**< source char parameter */
   SCIP_PARAM*           targetparam,        /**< target char parameter */
   SCIP*                 targetscip          /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   )
{
   char value;

   assert(sourceparam != NULL);
   assert(targetparam != NULL);

   /* get value of source parameter and copy it to target parameter */
   value = SCIPparamGetChar(sourceparam);
   SCIP_CALL( SCIPparamSetChar(targetparam,targetscip, value, TRUE) );

   return SCIP_OKAY;
}

/** copies value of source string parameter to target string parameter*/
static
SCIP_RETCODE paramCopyString(
   SCIP_PARAM*           sourceparam,        /**< source string parameter */
   SCIP_PARAM*           targetparam,        /**< target string parameter */
   SCIP*                 targetscip          /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   )
{
   char* value;

   assert(sourceparam != NULL);
   assert(targetparam != NULL);

   /* get value of source parameter and copy it to target parameter */
   value = SCIPparamGetString(sourceparam);
   SCIP_CALL( SCIPparamSetString(targetparam,targetscip, value, TRUE) );

   return SCIP_OKAY;
}

/** returns type of parameter */
SCIP_PARAMTYPE SCIPparamGetType(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);
   
   return param->paramtype;
}

/** returns name of parameter */
const char* SCIPparamGetName(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);

   return param->name;
}

/** returns description of parameter */
const char* SCIPparamGetDesc(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);

   return param->desc;
}

/** returns locally defined parameter specific data */
SCIP_PARAMDATA* SCIPparamGetData(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);

   return param->paramdata;
}

/** returns locally defined parameter specific data */
SCIP_Bool SCIPparamIsAdvanced(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);

   return param->isadvanced;
}

/** returns value of SCIP_Bool parameter */
SCIP_Bool SCIPparamGetBool(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_BOOL);

   if( param->data.boolparam.valueptr != NULL )
      return *param->data.boolparam.valueptr;
   else
      return param->data.boolparam.curvalue;
}

/** returns default value of SCIP_Bool parameter */
SCIP_Bool SCIPparamGetBoolDefault(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_BOOL);

   return param->data.boolparam.defaultvalue;
}

/** returns value of int parameter */
int SCIPparamGetInt(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_INT);

   if( param->data.intparam.valueptr != NULL )
      return *param->data.intparam.valueptr;
   else
      return param->data.intparam.curvalue;
}

/** returns minimal value of int parameter */
int SCIPparamGetIntMin(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_INT);

   return param->data.intparam.minvalue;
}

/** returns maximal value of int parameter */
int SCIPparamGetIntMax(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_INT);

   return param->data.intparam.maxvalue;
}

/** returns default value of int parameter */
int SCIPparamGetIntDefault(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_INT);

   return param->data.intparam.defaultvalue;
}

/** returns value of SCIP_Longint parameter */
SCIP_Longint SCIPparamGetLongint(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_LONGINT);

   if( param->data.longintparam.valueptr != NULL )
      return *param->data.longintparam.valueptr;
   else
      return param->data.longintparam.curvalue;
}

/** returns minimal value of longint parameter */
SCIP_Longint SCIPparamGetLongintMin(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_LONGINT);

   return param->data.longintparam.minvalue;
}

/** returns maximal value of longint parameter */
SCIP_Longint SCIPparamGetLongintMax(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_LONGINT);

   return param->data.longintparam.maxvalue;
}

/** returns default value of SCIP_Longint parameter */
SCIP_Longint SCIPparamGetLongintDefault(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_LONGINT);

   return param->data.longintparam.defaultvalue;
}

/** returns value of SCIP_Real parameter */
SCIP_Real SCIPparamGetReal(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_REAL);

   if( param->data.realparam.valueptr != NULL )
      return *param->data.realparam.valueptr;
   else
      return param->data.realparam.curvalue;
}

/** returns minimal value of real parameter */
SCIP_Real SCIPparamGetRealMin(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_REAL);

   return param->data.realparam.minvalue;
}

/** returns maximal value of real parameter */
SCIP_Real SCIPparamGetRealMax(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_REAL);

   return param->data.realparam.maxvalue;
}

/** returns default value of SCIP_Real parameter */
SCIP_Real SCIPparamGetRealDefault(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_REAL);

   return param->data.realparam.defaultvalue;
}

/** returns value of char parameter */
char SCIPparamGetChar(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_CHAR);

   if( param->data.charparam.valueptr != NULL )
      return *param->data.charparam.valueptr;
   else
      return param->data.charparam.curvalue;
}

/** returns default value of char parameter */
char SCIPparamGetCharDefault(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_CHAR);

   return param->data.charparam.defaultvalue;
}

/** returns value of string parameter */
char* SCIPparamGetString(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_STRING);

   if( param->data.stringparam.valueptr != NULL )
      return *param->data.stringparam.valueptr;
   else
      return param->data.stringparam.curvalue;
}

/** returns default value of String parameter */
char* SCIPparamGetStringDefault(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_STRING);

   return param->data.stringparam.defaultvalue;
}

/** sets value of SCIP_Bool parameter */
SCIP_RETCODE  SCIPparamSetBool(
   SCIP_PARAM*           param,              /**< parameter */
   SCIP*                 scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   SCIP_Bool             value,              /**< new value of the parameter */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   assert(param != NULL);

   /* check, if value is possible for the parameter */
   SCIP_CALL_QUIET( paramCheckBool(param, value) );

   /* set the parameter's current value */
   if( param->data.boolparam.valueptr != NULL )
      *param->data.boolparam.valueptr = value;
   else
      param->data.boolparam.curvalue = value;

   /* call the parameter's change information method */
   if( param->paramchgd != NULL && scip != NULL )
   {
      SCIP_CALL( param->paramchgd(scip, param) );
   }

   if( !quiet )
   {
      SCIP_CALL( paramWrite(param, NULL, FALSE, TRUE) );
   }

   return SCIP_OKAY;
}

/** sets value of int parameter */
SCIP_RETCODE SCIPparamSetInt(
   SCIP_PARAM*           param,              /**< parameter */
   SCIP*                 scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   int                   value,              /**< new value of the parameter */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   assert(param != NULL);

   /* check, if value is possible for the parameter */
   SCIP_CALL_QUIET( paramCheckInt(param, value) );

   /* set the parameter's current value */
   if( param->data.intparam.valueptr != NULL )
      *param->data.intparam.valueptr = value;
   else
      param->data.intparam.curvalue = value;

   /* call the parameter's change information method */
   if( param->paramchgd != NULL && scip != NULL )
   {
      SCIP_CALL( param->paramchgd(scip, param) );
   }

   if( !quiet )
   {
      SCIP_CALL( paramWrite(param, NULL, FALSE, TRUE) );
   }

   return SCIP_OKAY;
}

/** sets value of SCIP_Longint parameter */
SCIP_RETCODE SCIPparamSetLongint(
   SCIP_PARAM*           param,              /**< parameter */
   SCIP*                 scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   SCIP_Longint          value,              /**< new value of the parameter */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   assert(param != NULL);

   /* check, if value is possible for the parameter */
   SCIP_CALL_QUIET( paramCheckLongint(param, value) );

   /* set the parameter's current value */
   if( param->data.longintparam.valueptr != NULL )
      *param->data.longintparam.valueptr = value;
   else
      param->data.longintparam.curvalue = value;

   /* call the parameter's change information method */
   if( param->paramchgd != NULL && scip != NULL )
   {
      SCIP_CALL( param->paramchgd(scip, param) );
   }

   if( !quiet )
   {
      SCIP_CALL( paramWrite(param, NULL, FALSE, TRUE) );
   }

   return SCIP_OKAY;
}

/** sets value of SCIP_Real parameter */
SCIP_RETCODE SCIPparamSetReal(
   SCIP_PARAM*           param,              /**< parameter */
   SCIP*                 scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   SCIP_Real             value,              /**< new value of the parameter */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   assert(param != NULL);

   /* check, if value is possible for the parameter */
   value = MAX(value, SCIP_REAL_MIN);
   value = MIN(value, SCIP_REAL_MAX);
   SCIP_CALL_QUIET( paramCheckReal(param, value) );

   /* set the parameter's current value */
   if( param->data.realparam.valueptr != NULL )
      *param->data.realparam.valueptr = value;
   else
      param->data.realparam.curvalue = value;

   /* call the parameter's change information method */
   if( param->paramchgd != NULL && scip != NULL )
   {
      SCIP_CALL( param->paramchgd(scip, param) );
   }

   if( !quiet )
   {
      SCIP_CALL( paramWrite(param, NULL, FALSE, TRUE) );
   }

   return SCIP_OKAY;
}

/** sets value of char parameter */
SCIP_RETCODE SCIPparamSetChar(
   SCIP_PARAM*           param,              /**< parameter */
   SCIP*                 scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   char                  value,              /**< new value of the parameter */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   assert(param != NULL);

   /* check, if value is possible for the parameter */
   SCIP_CALL_QUIET( paramCheckChar(param, value) );

   /* set the parameter's current value */
   if( param->data.charparam.valueptr != NULL )
      *param->data.charparam.valueptr = value;
   else
      param->data.charparam.curvalue = value;

   /* call the parameter's change information method */
   if( param->paramchgd != NULL && scip != NULL )
   {
      SCIP_CALL( param->paramchgd(scip, param) );
   }

   if( !quiet )
   {
      SCIP_CALL( paramWrite(param, NULL, FALSE, TRUE) );
   }

   return SCIP_OKAY;
}

/** sets value of string parameter */
SCIP_RETCODE SCIPparamSetString(
   SCIP_PARAM*           param,              /**< parameter */
   SCIP*                 scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   const char*           value,              /**< new value of the parameter */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   assert(param != NULL);

   /* check, if value is possible for the parameter */
   SCIP_CALL_QUIET( paramCheckString(param, value) );

   /* set the parameter's current value */
   if( param->data.stringparam.valueptr != NULL )
   {
      BMSfreeMemoryArrayNull(param->data.stringparam.valueptr);
      BMSduplicateMemoryArray(param->data.stringparam.valueptr, value, strlen(value)+1);
   }
   else
   {
      BMSfreeMemoryArrayNull(&param->data.stringparam.curvalue);
      BMSduplicateMemoryArray(&param->data.stringparam.curvalue, value, strlen(value)+1);
   }

   /* call the parameter's change information method */
   if( param->paramchgd != NULL && scip != NULL )
   {
      SCIP_CALL( param->paramchgd(scip, param) );
   }

   if( !quiet )
   {
      SCIP_CALL( paramWrite(param, NULL, FALSE, TRUE) );
   }
   
   return SCIP_OKAY;
}

/** returns whether the parameter is on its default setting */
SCIP_Bool SCIPparamIsDefault(
   SCIP_PARAM*           param               /**< parameter */
   )
{
   assert(param != NULL);

   switch( param->paramtype )
   {
   case SCIP_PARAMTYPE_BOOL:
      return (SCIPparamGetBool(param) == SCIPparamGetBoolDefault(param));

   case SCIP_PARAMTYPE_INT:
      return (SCIPparamGetInt(param) == SCIPparamGetIntDefault(param));

   case SCIP_PARAMTYPE_LONGINT:
      return (SCIPparamGetLongint(param) == SCIPparamGetLongintDefault(param));

   case SCIP_PARAMTYPE_REAL:
      return EPSZ(SCIPparamGetReal(param) - SCIPparamGetRealDefault(param), 1e-16);

   case SCIP_PARAMTYPE_CHAR:
      return (SCIPparamGetChar(param) == SCIPparamGetCharDefault(param));

   case SCIP_PARAMTYPE_STRING:
      return (strcmp(SCIPparamGetString(param), SCIPparamGetStringDefault(param)) == 0);
      
   default:
      SCIPerrorMessage("unknown parameter type\n");
      SCIPABORT();
      return FALSE; /*lint !e527*/
   }
}

/** sets the parameter to its default setting */
SCIP_RETCODE SCIPparamSetToDefault(
   SCIP_PARAM*           param,              /**< parameter */
   SCIP*                 scip                /**< SCIP data structure, or NULL if paramchgd method should not be called */     
   )
{
   assert(param != NULL);

   switch( param->paramtype )
   {
   case SCIP_PARAMTYPE_BOOL:
      SCIP_CALL( SCIPparamSetBool(param, scip, SCIPparamGetBoolDefault(param), TRUE) );
      break;

   case SCIP_PARAMTYPE_INT:
      SCIP_CALL( SCIPparamSetInt(param, scip, SCIPparamGetIntDefault(param), TRUE) );
      break;

   case SCIP_PARAMTYPE_LONGINT:
      SCIP_CALL( SCIPparamSetLongint(param, scip, SCIPparamGetLongintDefault(param), TRUE) );
      break;

   case SCIP_PARAMTYPE_REAL:
      SCIP_CALL( SCIPparamSetReal(param, scip, SCIPparamGetRealDefault(param), TRUE) );
      break;

   case SCIP_PARAMTYPE_CHAR:
      SCIP_CALL( SCIPparamSetChar(param, scip, SCIPparamGetCharDefault(param), TRUE) );
      break;

   case SCIP_PARAMTYPE_STRING:
      SCIP_CALL( SCIPparamSetString(param, scip, SCIPparamGetStringDefault(param), TRUE) );
      break;
      
   default:
      SCIPerrorMessage("unknown parameter type\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** creates a parameter with name and description, does not set the type specific parameter values themselves */
static
SCIP_RETCODE paramCreate(
   SCIP_PARAM**          param,              /**< pointer to the parameter */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(param != NULL);
   assert(name != NULL);
   assert(desc != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, param) );
   
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*param)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*param)->desc, desc, strlen(desc)+1) );

   (*param)->paramchgd = paramchgd;
   (*param)->paramdata = paramdata;

   return SCIP_OKAY;
}

/** creates a SCIP_Bool parameter, and sets its value to default */
static
SCIP_RETCODE paramCreateBool(
   SCIP_PARAM**          param,              /**< pointer to the parameter */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Bool*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Bool             defaultvalue,       /**< default value of the parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(param != NULL);
   assert(name != NULL);

   SCIP_CALL( paramCreate(param, blkmem, name, desc, paramchgd, paramdata) );

   (*param)->paramtype = SCIP_PARAMTYPE_BOOL;
   (*param)->data.boolparam.valueptr = valueptr;
   (*param)->isadvanced = isadvanced;
   (*param)->data.boolparam.defaultvalue = defaultvalue;

   SCIP_CALL( SCIPparamSetBool(*param, NULL, defaultvalue, TRUE) );

   return SCIP_OKAY;
}

/** creates a int parameter, and sets its value to default */
static
SCIP_RETCODE paramCreateInt(
   SCIP_PARAM**          param,              /**< pointer to the parameter */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   int*                  valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   int                   defaultvalue,       /**< default value of the parameter */
   int                   minvalue,           /**< minimum value for parameter */
   int                   maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(param != NULL);
   assert(name != NULL);

   SCIP_CALL( paramCreate(param, blkmem, name, desc, paramchgd, paramdata) );

   (*param)->paramtype = SCIP_PARAMTYPE_INT;
   (*param)->data.intparam.valueptr = valueptr;
   (*param)->isadvanced = isadvanced;
   (*param)->data.intparam.defaultvalue = defaultvalue;
   (*param)->data.intparam.minvalue = minvalue;
   (*param)->data.intparam.maxvalue = maxvalue;

   SCIP_CALL( SCIPparamSetInt(*param, NULL, defaultvalue, TRUE) );

   return SCIP_OKAY;
}

/** creates a SCIP_Longint parameter, and sets its value to default */
static
SCIP_RETCODE paramCreateLongint(
   SCIP_PARAM**          param,              /**< pointer to the parameter */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Longint*         valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Longint          defaultvalue,       /**< default value of the parameter */
   SCIP_Longint          minvalue,           /**< minimum value for parameter */
   SCIP_Longint          maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(param != NULL);
   assert(name != NULL);

   SCIP_CALL( paramCreate(param, blkmem, name, desc, paramchgd, paramdata) );

   (*param)->paramtype = SCIP_PARAMTYPE_LONGINT;
   (*param)->data.longintparam.valueptr = valueptr;
   (*param)->isadvanced = isadvanced;
   (*param)->data.longintparam.defaultvalue = defaultvalue;
   (*param)->data.longintparam.minvalue = minvalue;
   (*param)->data.longintparam.maxvalue = maxvalue;

   SCIP_CALL( SCIPparamSetLongint(*param, NULL, defaultvalue, TRUE) );

   return SCIP_OKAY;
}

/** creates a SCIP_Real parameter, and sets its value to default */
static
SCIP_RETCODE paramCreateReal(
   SCIP_PARAM**          param,              /**< pointer to the parameter */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Real*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Real             defaultvalue,       /**< default value of the parameter */
   SCIP_Real             minvalue,           /**< minimum value for parameter */
   SCIP_Real             maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(param != NULL);
   assert(name != NULL);

   SCIP_CALL( paramCreate(param, blkmem, name, desc, paramchgd, paramdata) );

   (*param)->paramtype = SCIP_PARAMTYPE_REAL;
   (*param)->data.realparam.valueptr = valueptr;
   (*param)->isadvanced = isadvanced;
   (*param)->data.realparam.defaultvalue = defaultvalue;
   (*param)->data.realparam.minvalue = minvalue;
   (*param)->data.realparam.maxvalue = maxvalue;

   SCIP_CALL( SCIPparamSetReal(*param, NULL, defaultvalue, TRUE) );

   return SCIP_OKAY;
}

/** creates a char parameter, and sets its value to default */
static
SCIP_RETCODE paramCreateChar(
   SCIP_PARAM**          param,              /**< pointer to the parameter */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   char*                 valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   char                  defaultvalue,       /**< default value of the parameter */
   const char*           allowedvalues,      /**< array with possible parameter values, or NULL if not restricted */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(param != NULL);
   assert(name != NULL);

   SCIP_CALL( paramCreate(param, blkmem, name, desc, paramchgd, paramdata) );

   (*param)->paramtype = SCIP_PARAMTYPE_CHAR;
   (*param)->data.charparam.valueptr = valueptr;
   (*param)->isadvanced = isadvanced;
   (*param)->data.charparam.defaultvalue = defaultvalue;
   if( allowedvalues != NULL )
   {
      SCIP_ALLOC( BMSduplicateMemoryArray(&(*param)->data.charparam.allowedvalues, allowedvalues, strlen(allowedvalues)+1) );
   }
   else
      (*param)->data.charparam.allowedvalues = NULL;

   SCIP_CALL( SCIPparamSetChar(*param, NULL, defaultvalue, TRUE) );

   return SCIP_OKAY;
}

/** creates a string parameter, and sets its value to default */
static
SCIP_RETCODE paramCreateString(
   SCIP_PARAM**          param,              /**< pointer to the parameter */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   char**                valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   const char*           defaultvalue,       /**< default value of the parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(param != NULL);
   assert(name != NULL);
   assert(valueptr == NULL || *valueptr == NULL);
   assert(defaultvalue != NULL);

   SCIP_CALL( paramCreate(param, blkmem, name, desc, paramchgd, paramdata) );

   (*param)->paramtype = SCIP_PARAMTYPE_STRING;
   (*param)->data.stringparam.valueptr = valueptr;
   (*param)->isadvanced = isadvanced;
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*param)->data.stringparam.defaultvalue, defaultvalue, strlen(defaultvalue)+1) );
   (*param)->data.stringparam.curvalue = NULL;

   SCIP_CALL( SCIPparamSetString(*param, NULL, defaultvalue, TRUE) );

   return SCIP_OKAY;
}

/** frees a single parameter */
static
void paramFree(
   SCIP_PARAM**          param,              /**< pointer to the parameter */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(param != NULL);
   assert(*param != NULL);

   switch( (*param)->paramtype )
   {
   case SCIP_PARAMTYPE_BOOL:
   case SCIP_PARAMTYPE_INT:
   case SCIP_PARAMTYPE_LONGINT:
   case SCIP_PARAMTYPE_REAL:
      break;
   case SCIP_PARAMTYPE_CHAR:
      BMSfreeMemoryArrayNull(&(*param)->data.charparam.allowedvalues);
      break;
   case SCIP_PARAMTYPE_STRING:
      BMSfreeMemoryArray(&(*param)->data.stringparam.defaultvalue);
      if( (*param)->data.stringparam.valueptr == NULL )
      {
         BMSfreeMemoryArray(&(*param)->data.stringparam.curvalue);
      }
      else
      {
         BMSfreeMemoryArray((*param)->data.stringparam.valueptr);
      }
      break;
   default:
      SCIPerrorMessage("invalid parameter type\n");
      SCIPABORT();
   }

   BMSfreeMemoryArray(&(*param)->name);
   BMSfreeMemoryArray(&(*param)->desc);
   BMSfreeBlockMemory(blkmem, param);
}

/** sets SCIP_Bool parameter according to the value of the given string */
static
SCIP_RETCODE paramParseBool(
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_SET*             set,                /**< global SCIP settings */
   char*                 valuestr            /**< value in string format (may be modified during parse) */
   )
{
   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_BOOL);
   assert(set != NULL);
   assert(valuestr != NULL);

   if( strcasecmp(valuestr, "TRUE") == 0 )
   {
      SCIP_CALL( SCIPparamSetBool(param, set->scip, TRUE, TRUE) );
   }
   else if( strcasecmp(valuestr, "FALSE") == 0 )
   {
      SCIP_CALL( SCIPparamSetBool(param, set->scip, FALSE, TRUE) );
   }
   else
   {
      SCIPerrorMessage("invalid parameter value <%s> for SCIP_Bool parameter <%s>\n", valuestr, param->name);
      return SCIP_READERROR;
   }
   
   return SCIP_OKAY;
}

/** sets int parameter according to the value of the given string */
static
SCIP_RETCODE paramParseInt(
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_SET*             set,                /**< global SCIP settings */
   char*                 valuestr            /**< value in string format (may be modified during parse) */
   )
{
   int value;

   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_INT);
   assert(set != NULL);
   assert(valuestr != NULL);

   if( sscanf(valuestr, "%d", &value) == 1 )
   {
      SCIP_CALL( SCIPparamSetInt(param, set->scip, value, TRUE) );
   }
   else
   {
      SCIPerrorMessage("invalid parameter value <%s> for int parameter <%s>\n", valuestr, param->name);
      return SCIP_READERROR;
   }
   
   return SCIP_OKAY;
}

/** sets SCIP_Longint parameter according to the value of the given string */
static
SCIP_RETCODE paramParseLongint(
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_SET*             set,                /**< global SCIP settings */
   char*                 valuestr            /**< value in string format (may be modified during parse) */
   )
{
   SCIP_Longint value;

   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_LONGINT);
   assert(set != NULL);
   assert(valuestr != NULL);

   if( sscanf(valuestr, "%"SCIP_LONGINT_FORMAT, &value) == 1 )
   {
      SCIP_CALL( SCIPparamSetLongint(param, set->scip, value, TRUE) );
   }
   else
   {
      SCIPerrorMessage("invalid parameter value <%s> for SCIP_Longint parameter <%s>\n", valuestr, param->name);
      return SCIP_READERROR;
   }
   
   return SCIP_OKAY;
}

/** sets SCIP_Real parameter according to the value of the given string */
static
SCIP_RETCODE paramParseReal(
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_SET*             set,                /**< global SCIP settings */
   char*                 valuestr            /**< value in string format (may be modified during parse) */
   )
{
   SCIP_Real value;

   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_REAL);
   assert(set != NULL);
   assert(valuestr != NULL);

   if( sscanf(valuestr, "%"SCIP_REAL_FORMAT, &value) == 1 )
   {
      SCIP_CALL( SCIPparamSetReal(param, set->scip, value, TRUE) );
   }
   else
   {
      SCIPerrorMessage("invalid parameter value <%s> for SCIP_Real parameter <%s>\n", valuestr, param->name);
      return SCIP_READERROR;
   }
   
   return SCIP_OKAY;
}

/** sets Char parameter according to the value of the given string */
static
SCIP_RETCODE paramParseChar(
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_SET*             set,                /**< global SCIP settings */
   char*                 valuestr            /**< value in string format (may be modified during parse) */
   )
{
   char value;

   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_CHAR);
   assert(set != NULL);
   assert(valuestr != NULL);

   if( sscanf(valuestr, "%c", &value) == 1 )
   {
      SCIP_CALL( SCIPparamSetChar(param, set->scip, value, TRUE) );
   }
   else
   {
      SCIPerrorMessage("invalid parameter value <%s> for char parameter <%s>\n", valuestr, param->name);
      return SCIP_READERROR;
   }
   
   return SCIP_OKAY;
}

/** sets String parameter according to the value of the given string */
static
SCIP_RETCODE paramParseString(
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_SET*             set,                /**< global SCIP settings */
   char*                 valuestr            /**< value in string format (may be modified during parse) */
   )
{
   unsigned int len;

   assert(param != NULL);
   assert(param->paramtype == SCIP_PARAMTYPE_STRING);
   assert(set != NULL);
   assert(valuestr != NULL);

   /* check for quotes */
   len = strlen(valuestr);
   if( len <= 1 || valuestr[0] != '"' || valuestr[len-1] != '"' )
   {
      SCIPerrorMessage("invalid parameter value <%s> for string parameter <%s> (string has to be in double quotes)\n",
         valuestr, param->name);
      return SCIP_READERROR;
   }

   /* remove the quotes */
   valuestr[len-1] = '\0';
   valuestr++;
   SCIP_CALL( SCIPparamSetString(param, set->scip, valuestr, TRUE) );
   
   return SCIP_OKAY;
}


/*
 * Parameter set methods
 */

/** creates parameter set */
SCIP_RETCODE SCIPparamsetCreate(
   SCIP_PARAMSET**       paramset,           /**< pointer to store the parameter set */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(paramset != NULL);

   SCIP_ALLOC( BMSallocMemory(paramset) );

   SCIP_CALL( SCIPhashtableCreate(&(*paramset)->hashtable, blkmem, SCIP_HASHSIZE_PARAMS,
         hashGetKeyParam, SCIPhashKeyEqString, SCIPhashKeyValString, NULL) );

   (*paramset)->params = NULL;
   (*paramset)->nparams = 0;
   (*paramset)->paramssize = 0;

   return SCIP_OKAY;
}

/** frees parameter set */
void SCIPparamsetFree(
   SCIP_PARAMSET**       paramset,           /**< pointer to the parameter set */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int i;

   assert(paramset != NULL);
   assert(*paramset != NULL);
   assert((*paramset)->paramssize == 0 || (*paramset)->params != NULL);
   assert((*paramset)->paramssize >= (*paramset)->nparams);

   for( i = 0; i < (*paramset)->nparams; ++i )
   {
      paramFree(&(*paramset)->params[i], blkmem);
   }

   SCIPhashtableFree(&(*paramset)->hashtable);

   BMSfreeMemoryArrayNull(&(*paramset)->params);
   BMSfreeMemory(paramset);
}

/** adds parameter to the parameter set */
static
SCIP_RETCODE paramsetAdd(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP_PARAM*           param               /**< parameter to add */
   )
{
   assert(paramset != NULL);
   assert(param != NULL);

   /* insert the parameter name to the hash table */
   SCIP_CALL( SCIPhashtableSafeInsert(paramset->hashtable, (void*)param) );

   /* ensure, that there is enough space in the params array */
   if( paramset->nparams >= paramset->paramssize )
   {
      paramset->paramssize *= 2;
      paramset->paramssize = MAX(paramset->paramssize, paramset->nparams+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&paramset->params, paramset->paramssize) );
   }
   assert(paramset->nparams < paramset->paramssize);
   
   /* insert parameter in the params array */
   paramset->params[paramset->nparams] = param;
   paramset->nparams++;
   
   return SCIP_OKAY;
}

/** creates a SCIP_Bool parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPparamsetAddBool(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Bool*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Bool             defaultvalue,       /**< default value of the parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   SCIP_PARAM* param;

   assert(paramset != NULL);

   /* create the parameter */
   SCIP_CALL( paramCreateBool(&param, blkmem, name, desc, valueptr, isadvanced, defaultvalue, paramchgd, paramdata) );

   /* add parameter to the parameter set */
   SCIP_CALL( paramsetAdd(paramset, param) );
   
   return SCIP_OKAY;
}

/** creates a int parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPparamsetAddInt(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   int*                  valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   int                   defaultvalue,       /**< default value of the parameter */
   int                   minvalue,           /**< minimum value for parameter */
   int                   maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   SCIP_PARAM* param;

   assert(paramset != NULL);

   /* create the parameter */
   SCIP_CALL( paramCreateInt(&param, blkmem, name, desc, valueptr, isadvanced, defaultvalue, minvalue, maxvalue, 
         paramchgd, paramdata) );

   /* add parameter to the parameter set */
   SCIP_CALL( paramsetAdd(paramset, param) );
   
   return SCIP_OKAY;
}

/** creates a SCIP_Longint parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPparamsetAddLongint(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Longint*         valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Longint          defaultvalue,       /**< default value of the parameter */
   SCIP_Longint          minvalue,           /**< minimum value for parameter */
   SCIP_Longint          maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   SCIP_PARAM* param;

   assert(paramset != NULL);

   /* create the parameter */
   SCIP_CALL( paramCreateLongint(&param, blkmem, name, desc, valueptr, isadvanced, defaultvalue, minvalue, maxvalue, 
         paramchgd, paramdata) );

   /* add parameter to the parameter set */
   SCIP_CALL( paramsetAdd(paramset, param) );
   
   return SCIP_OKAY;
}

/** creates a SCIP_Real parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPparamsetAddReal(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Real*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Real             defaultvalue,       /**< default value of the parameter */
   SCIP_Real             minvalue,           /**< minimum value for parameter */
   SCIP_Real             maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   SCIP_PARAM* param;

   assert(paramset != NULL);

   /* create the parameter */
   SCIP_CALL( paramCreateReal(&param, blkmem, name, desc, valueptr, isadvanced, defaultvalue, minvalue, maxvalue, 
         paramchgd, paramdata) );

   /* add parameter to the parameter set */
   SCIP_CALL( paramsetAdd(paramset, param) );
   
   return SCIP_OKAY;
}

/** creates a char parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPparamsetAddChar(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   char*                 valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   char                  defaultvalue,       /**< default value of the parameter */
   const char*           allowedvalues,      /**< array with possible parameter values, or NULL if not restricted */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   SCIP_PARAM* param;

   assert(paramset != NULL);

   /* create the parameter */
   SCIP_CALL( paramCreateChar(&param, blkmem, name, desc, valueptr, isadvanced, defaultvalue, allowedvalues, 
         paramchgd, paramdata) );

   /* add parameter to the parameter set */
   SCIP_CALL( paramsetAdd(paramset, param) );
   
   return SCIP_OKAY;
}

/** creates a string parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPparamsetAddString(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   char**                valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   const char*           defaultvalue,       /**< default value of the parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   SCIP_PARAM* param;

   assert(paramset != NULL);

   /* create the parameter */
   SCIP_CALL( paramCreateString(&param, blkmem, name, desc, valueptr, isadvanced, defaultvalue, paramchgd, paramdata) );

   /* add parameter to the parameter set */
   SCIP_CALL( paramsetAdd(paramset, param) );
   
   return SCIP_OKAY;
}

/** returns the name of the given parameter type */
static
const char* paramtypeGetName(
   SCIP_PARAMTYPE        paramtype           /**< type of parameter */
   )
{
   static const char* paramtypename[] = {
      "Bool",    /* SCIP_PARAMTYPE_BOOL    = 0 */
      "int",     /* SCIP_PARAMTYPE_INT     = 1 */
      "Longint", /* SCIP_PARAMTYPE_LONGINT = 2 */
      "Real",    /* SCIP_PARAMTYPE_REAL    = 3 */
      "char",    /* SCIP_PARAMTYPE_CHAR    = 4 */
      "string"   /* SCIP_PARAMTYPE_STRING  = 5 */
   };

   return paramtypename[(int)paramtype];
}

/** gets the value of an existing SCIP_Bool parameter */
SCIP_RETCODE SCIPparamsetGetBool(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   const char*           name,               /**< name of the parameter */
   SCIP_Bool*            value               /**< pointer to store the parameter */
   )
{
   SCIP_PARAM* param;

   assert(paramset != NULL);
   assert(value != NULL);

   /* retrieve parameter from hash table */
   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
   {
      SCIPerrorMessage("parameter <%s> unknown\n", name);
      return SCIP_PARAMETERUNKNOWN;
   }
   if( param->paramtype != SCIP_PARAMTYPE_BOOL )
   {
      SCIPerrorMessage("wrong parameter type - parameter <%s> has type <%s> instead of <%s>\n", 
         name, paramtypeGetName(param->paramtype), paramtypeGetName(SCIP_PARAMTYPE_BOOL));
      return SCIP_PARAMETERWRONGTYPE;
   }

   /* get the parameter's current value */
   *value = SCIPparamGetBool(param);

   return SCIP_OKAY;
}

/** gets the value of an existing int parameter */
SCIP_RETCODE SCIPparamsetGetInt(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   const char*           name,               /**< name of the parameter */
   int*                  value               /**< pointer to store the parameter */
   )
{
   SCIP_PARAM* param;

   assert(paramset != NULL);
   assert(value != NULL);

   /* retrieve parameter from hash table */
   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
   {
      SCIPerrorMessage("parameter <%s> unknown\n", name);
      return SCIP_PARAMETERUNKNOWN;
   }
   if( param->paramtype != SCIP_PARAMTYPE_INT )
   {
      SCIPerrorMessage("wrong parameter type - parameter <%s> has type <%s> instead of <%s>\n", 
         name, paramtypeGetName(param->paramtype), paramtypeGetName(SCIP_PARAMTYPE_INT));
      return SCIP_PARAMETERWRONGTYPE;
   }

   /* get the parameter's current value */
   *value = SCIPparamGetInt(param);

   return SCIP_OKAY;
}

/** gets the value of an existing SCIP_Longint parameter */
SCIP_RETCODE SCIPparamsetGetLongint(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   const char*           name,               /**< name of the parameter */
   SCIP_Longint*         value               /**< pointer to store the parameter */
   )
{
   SCIP_PARAM* param;

   assert(paramset != NULL);
   assert(value != NULL);

   /* retrieve parameter from hash table */
   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
   {
      SCIPerrorMessage("parameter <%s> unknown\n", name);
      return SCIP_PARAMETERUNKNOWN;
   }
   if( param->paramtype != SCIP_PARAMTYPE_LONGINT )
   {
      SCIPerrorMessage("wrong parameter type - parameter <%s> has type <%s> instead of <%s>\n", 
         name, paramtypeGetName(param->paramtype), paramtypeGetName(SCIP_PARAMTYPE_LONGINT));
      return SCIP_PARAMETERWRONGTYPE;
   }

   /* get the parameter's current value */
   *value = SCIPparamGetLongint(param);

   return SCIP_OKAY;
}

/** gets the value of an existing SCIP_Real parameter */
SCIP_RETCODE SCIPparamsetGetReal(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   const char*           name,               /**< name of the parameter */
   SCIP_Real*            value               /**< pointer to store the parameter */
   )
{
   SCIP_PARAM* param;

   assert(paramset != NULL);
   assert(value != NULL);

   /* retrieve parameter from hash table */
   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
   {
      SCIPerrorMessage("parameter <%s> unknown\n", name);
      return SCIP_PARAMETERUNKNOWN;
   }
   if( param->paramtype != SCIP_PARAMTYPE_REAL )
   {
      SCIPerrorMessage("wrong parameter type - parameter <%s> has type <%s> instead of <%s>\n", 
         name, paramtypeGetName(param->paramtype), paramtypeGetName(SCIP_PARAMTYPE_REAL));
      return SCIP_PARAMETERWRONGTYPE;
   }

   /* get the parameter's current value */
   *value = SCIPparamGetReal(param);

   return SCIP_OKAY;
}

/** gets the value of an existing char parameter */
SCIP_RETCODE SCIPparamsetGetChar(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   const char*           name,               /**< name of the parameter */
   char*                 value               /**< pointer to store the parameter */
   )
{
   SCIP_PARAM* param;

   assert(paramset != NULL);
   assert(value != NULL);

   /* retrieve parameter from hash table */
   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
   {
      SCIPerrorMessage("parameter <%s> unknown\n", name);
      return SCIP_PARAMETERUNKNOWN;
   }
   if( param->paramtype != SCIP_PARAMTYPE_CHAR )
   {
      SCIPerrorMessage("wrong parameter type - parameter <%s> has type <%s> instead of <%s>\n", 
         name, paramtypeGetName(param->paramtype), paramtypeGetName(SCIP_PARAMTYPE_CHAR));
      return SCIP_PARAMETERWRONGTYPE;
   }

   /* get the parameter's current value */
   *value = SCIPparamGetChar(param);

   return SCIP_OKAY;
}

/** gets the value of an existing string parameter */
SCIP_RETCODE SCIPparamsetGetString(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   const char*           name,               /**< name of the parameter */
   char**                value               /**< pointer to store the parameter */
   )
{
   SCIP_PARAM* param;

   assert(paramset != NULL);
   assert(value != NULL);

   /* retrieve parameter from hash table */
   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
   {
      SCIPerrorMessage("parameter <%s> unknown\n", name);
      return SCIP_PARAMETERUNKNOWN;
   }
   if( param->paramtype != SCIP_PARAMTYPE_STRING )
   {
      SCIPerrorMessage("wrong parameter type - parameter <%s> has type <%s> instead of <%s>\n", 
         name, paramtypeGetName(param->paramtype), paramtypeGetName(SCIP_PARAMTYPE_STRING));
      return SCIP_PARAMETERWRONGTYPE;
   }

   /* get the parameter's current value */
   *value = SCIPparamGetString(param);

   return SCIP_OKAY;
}

/** changes the value of an existing parameter */
SCIP_RETCODE SCIPparamsetSet(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   void*                 value               /**< new value of the parameter */
   )
{
   SCIP_PARAM* param;
   
   assert(paramset != NULL);
   assert(set != NULL);

   /* retrieve parameter from hash table */
   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
   {
      SCIPerrorMessage("parameter <%s> unknown\n", name);
      return SCIP_PARAMETERUNKNOWN;
   }
   
   switch( param->paramtype )
   {
   case SCIP_PARAMTYPE_BOOL:
      /* set the parameter's current value */
      SCIP_CALL( SCIPparamSetBool(param, set->scip, (SCIP_Bool) (size_t) value, TRUE) );
      break;
      
   case  SCIP_PARAMTYPE_INT:
      /* set the parameter's current value */
      SCIP_CALL( SCIPparamSetInt(param, set->scip, (int) (size_t) value, TRUE) );
      break;
      
   case SCIP_PARAMTYPE_LONGINT:
      /* set the parameter's current value */
      SCIP_CALL( SCIPparamSetLongint(param, set->scip, (SCIP_Longint) (size_t) value, TRUE) );
      break;
      
   case SCIP_PARAMTYPE_REAL:
      /* set the parameter's current value */
      SCIP_CALL( SCIPparamSetReal(param, set->scip, (SCIP_Real) (size_t) value, TRUE) );
      break;
      
   case SCIP_PARAMTYPE_CHAR:
      /* set the parameter's current value */
      SCIP_CALL( SCIPparamSetChar(param, set->scip, (char) (size_t) value, TRUE) );
      break;
      
   case SCIP_PARAMTYPE_STRING:
      /* set the parameter's current value */
      SCIP_CALL( SCIPparamSetString(param, set->scip, (char*) value, TRUE) );
      break;
      
   default:
      SCIPerrorMessage("unknown parameter type\n");
      return SCIP_INVALIDDATA;
   }
   
   return SCIP_OKAY;
}

/** changes the value of an existing SCIP_Bool parameter */
SCIP_RETCODE SCIPparamsetSetBool(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   SCIP_Bool             value               /**< new value of the parameter */
   )
{
   SCIP_PARAM* param;

   assert(paramset != NULL);
   assert(set != NULL);

   /* retrieve parameter from hash table */
   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
   {
      SCIPerrorMessage("parameter <%s> unknown\n", name);
      return SCIP_PARAMETERUNKNOWN;
   }
   if( param->paramtype != SCIP_PARAMTYPE_BOOL )
   {
      SCIPerrorMessage("wrong parameter type - parameter <%s> has type <%s> instead of <%s>\n", 
         name, paramtypeGetName(param->paramtype), paramtypeGetName(SCIP_PARAMTYPE_BOOL));
      return SCIP_PARAMETERWRONGTYPE;
   }

   /* set the parameter's current value */
   SCIP_CALL( SCIPparamSetBool(param, set->scip, value, TRUE) );
   
   return SCIP_OKAY;
}

/** changes the value of an existing int parameter */
SCIP_RETCODE SCIPparamsetSetInt(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   int                   value               /**< new value of the parameter */
   )
{
   SCIP_PARAM* param;

   assert(paramset != NULL);
   assert(set != NULL);

   /* retrieve parameter from hash table */
   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
   {
      SCIPerrorMessage("parameter <%s> unknown\n", name);
      return SCIP_PARAMETERUNKNOWN;
   }
   if( param->paramtype != SCIP_PARAMTYPE_INT )
   {
      SCIPerrorMessage("wrong parameter type - parameter <%s> has type <%s> instead of <%s>\n", 
         name, paramtypeGetName(param->paramtype), paramtypeGetName(SCIP_PARAMTYPE_INT));
      return SCIP_PARAMETERWRONGTYPE;
   }

   /* set the parameter's current value */
   SCIP_CALL( SCIPparamSetInt(param, set->scip, value, TRUE) );

   return SCIP_OKAY;
}

/** changes the value of an existing SCIP_Longint parameter */
SCIP_RETCODE SCIPparamsetSetLongint(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   SCIP_Longint          value               /**< new value of the parameter */
   )
{
   SCIP_PARAM* param;

   assert(paramset != NULL);
   assert(set != NULL);

   /* retrieve parameter from hash table */
   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
   {
      SCIPerrorMessage("parameter <%s> unknown\n", name);
      return SCIP_PARAMETERUNKNOWN;
   }
   if( param->paramtype != SCIP_PARAMTYPE_LONGINT )
   {
      SCIPerrorMessage("wrong parameter type - parameter <%s> has type <%s> instead of <%s>\n", 
         name, paramtypeGetName(param->paramtype), paramtypeGetName(SCIP_PARAMTYPE_LONGINT));
      return SCIP_PARAMETERWRONGTYPE;
   }

   /* set the parameter's current value */
   SCIP_CALL( SCIPparamSetLongint(param, set->scip, value, TRUE) );

   return SCIP_OKAY;
}

/** changes the value of an existing SCIP_Real parameter */
SCIP_RETCODE SCIPparamsetSetReal(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   SCIP_Real             value               /**< new value of the parameter */
   )
{
   SCIP_PARAM* param;

   assert(paramset != NULL);
   assert(set != NULL);

   /* retrieve parameter from hash table */
   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
   {
      SCIPerrorMessage("parameter <%s> unknown\n", name);
      return SCIP_PARAMETERUNKNOWN;
   }
   if( param->paramtype != SCIP_PARAMTYPE_REAL )
   {
      SCIPerrorMessage("wrong parameter type - parameter <%s> has type <%s> instead of <%s>\n", 
         name, paramtypeGetName(param->paramtype), paramtypeGetName(SCIP_PARAMTYPE_REAL));
      return SCIP_PARAMETERWRONGTYPE;
   }

   /* set the parameter's current value */
   SCIP_CALL( SCIPparamSetReal(param, set->scip, value, TRUE) );
   
   return SCIP_OKAY;
}

/** changes the value of an existing char parameter */
SCIP_RETCODE SCIPparamsetSetChar(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   char                  value               /**< new value of the parameter */
   )
{
   SCIP_PARAM* param;

   assert(paramset != NULL);
   assert(set != NULL);

   /* retrieve parameter from hash table */
   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
   {
      SCIPerrorMessage("parameter <%s> unknown\n", name);
      return SCIP_PARAMETERUNKNOWN;
   }
   if( param->paramtype != SCIP_PARAMTYPE_CHAR )
   {
      SCIPerrorMessage("wrong parameter type - parameter <%s> has type <%s> instead of <%s>\n", 
         name, paramtypeGetName(param->paramtype), paramtypeGetName(SCIP_PARAMTYPE_CHAR));
      return SCIP_PARAMETERWRONGTYPE;
   }

   /* set the parameter's current value */
   SCIP_CALL( SCIPparamSetChar(param, set->scip, value, TRUE) );

   return SCIP_OKAY;
}

/** changes the value of an existing string parameter */
SCIP_RETCODE SCIPparamsetSetString(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   const char*           value               /**< new value of the parameter */
   )
{
   SCIP_PARAM* param;

   assert(paramset != NULL);
   assert(set != NULL);

   /* retrieve parameter from hash table */
   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
   {
      SCIPerrorMessage("parameter <%s> unknown\n", name);
      return SCIP_PARAMETERUNKNOWN;
   }
   if( param->paramtype != SCIP_PARAMTYPE_STRING )
   {
      SCIPerrorMessage("wrong parameter type - parameter <%s> has type <%s> instead of <%s>\n", 
         name, paramtypeGetName(param->paramtype), paramtypeGetName(SCIP_PARAMTYPE_STRING));
      return SCIP_PARAMETERWRONGTYPE;
   }

   /* set the parameter's current value */
   SCIP_CALL( SCIPparamSetString(param, set->scip, value, TRUE) );

   return SCIP_OKAY;
}

/** parses a parameter file line "paramname = paramvalue" and sets parameter accordingly */
static
SCIP_RETCODE paramsetParse(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP_SET*             set,                /**< global SCIP settings */
   char*                 line                /**< line to parse (is modified during parse, but not freed) */
   )
{
   SCIP_PARAM* param;
   char* paramname;
   char* paramvaluestr;
   char* lastquote;
   SCIP_Bool quoted;

   assert(paramset != NULL);
   assert(line != NULL);

   /* find the start of the parameter name */
   while( *line == ' ' || *line == '\t' || *line == '\r' )
      line++;
   if( *line == '\0' || *line == '\n' || *line == '#' )
      return SCIP_OKAY;
   paramname = line;

   /* find the end of the parameter name */
   while( *line != ' ' && *line != '\t' && *line != '\r' && *line != '\n' && *line != '#' && *line != '\0' && *line != '=' )
      line++;
   if( *line == '=' )
   {
      *line = '\0';
      line++;
   }
   else
   {
      *line = '\0';
      line++;

      /* search for the '=' char in the line */
      while( *line == ' ' || *line == '\t' || *line == '\r' )
         line++;
      if( *line != '=' )
      {
         SCIPerrorMessage("character '=' was expected after the parameter name\n");
         return SCIP_READERROR;
      }
      line++;
   }

   /* find the start of the parameter value string */
   while( *line == ' ' || *line == '\t' || *line == '\r' )
      line++;
   if( *line == '\0' || *line == '\n' || *line == '#' )
   {
      SCIPerrorMessage("parameter value is missing\n");
      return SCIP_READERROR;
   }
   paramvaluestr = line;

   /* find the end of the parameter value string */
   quoted = (*paramvaluestr == '"');
   lastquote = NULL;
   while( (quoted || (*line != ' ' && *line != '\t' && *line != '\r' && *line != '\n' && *line != '#')) && *line != '\0' )
   {
      if( *line == '"' )
         lastquote = line;
      line++;
   }
   if( lastquote != NULL )
      line = lastquote+1;
   if( *line == '#' )
      *line = '\0';
   else if( *line != '\0' )
   {
      /* check, if the rest of the line is clean */
      *line = '\0';
      line++;
      while( *line == ' ' || *line == '\t' || *line == '\r' )
         line++;
      if( *line != '\0' && *line != '\n' && *line != '#' )
      {
         SCIPerrorMessage("additional characters after parameter value\n");
         return SCIP_READERROR;
      }
   }

   /* retrieve parameter from hash table */
   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)paramname);
   if( param == NULL )
   {
      SCIPwarningMessage("unknown parameter <%s>\n", paramname);
      return SCIP_OKAY;
   }

   /* set parameter's value */
   switch( param->paramtype )
   {
   case SCIP_PARAMTYPE_BOOL:
      SCIP_CALL( paramParseBool(param, set, paramvaluestr) );
      break;
   case SCIP_PARAMTYPE_INT:
      SCIP_CALL( paramParseInt(param, set, paramvaluestr) );
      break;
   case SCIP_PARAMTYPE_LONGINT:
      SCIP_CALL( paramParseLongint(param, set, paramvaluestr) );
      break;
   case SCIP_PARAMTYPE_REAL:
      SCIP_CALL( paramParseReal(param, set, paramvaluestr) );
      break;
   case SCIP_PARAMTYPE_CHAR:
      SCIP_CALL( paramParseChar(param, set, paramvaluestr) );
      break;
   case SCIP_PARAMTYPE_STRING:
      SCIP_CALL( paramParseString(param, set, paramvaluestr) );
      break;
   default:
      SCIPerrorMessage("unknown parameter type\n");
      return SCIP_INVALIDDATA;
   }
   
   return SCIP_OKAY;
}

/** reads parameters from a file */
SCIP_RETCODE SCIPparamsetRead(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           filename            /**< file name */
   )
{
   SCIP_RETCODE retcode;
   FILE* file;
   char line[1024];
   int lineno;

   assert(paramset != NULL);
   assert(filename != NULL);

   /* open the file for reading */
   file = fopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   /* read the parameters from the file */
   lineno = 0;
   retcode = SCIP_OKAY;
   while( fgets(line, sizeof(line), file) != NULL && retcode == SCIP_OKAY )
   {
      lineno++;
      retcode = paramsetParse(paramset, set, line);
   }

   /* close input file */
   fclose(file);

   if( retcode == SCIP_READERROR )
   {
      SCIPerrorMessage("input error in file <%s> line %d\n", filename, lineno);
   }
   else
   {
      SCIP_CALL( retcode );
   }
   
   return SCIP_OKAY;
}

/** writes all parameters in the parameter set to a file */
SCIP_RETCODE SCIPparamsetWrite(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   const char*           filename,           /**< file name, or NULL for stdout */
   SCIP_Bool             comments,           /**< should parameter descriptions be written as comments? */
   SCIP_Bool             onlychanged         /**< should only the parameters been written, that are changed from default? */
   )
{
   FILE* file;
   int i;

   assert(paramset != NULL);

   /* open the file for writing */
   if( filename != NULL )
   {
      file = fopen(filename, "w");
      if( file == NULL )
      {
         SCIPerrorMessage("cannot open file <%s> for writing\n", filename);
         SCIPprintSysError(filename);
         return SCIP_FILECREATEERROR;
      }
   }
   else
      file = NULL;

   if( comments )
   {
      /* display the SCIP version as comment in the first line */
#if( SCIP_SUBVERSION == 0 )
         SCIPmessageFPrintInfo(file, "# SCIP version %d.%d.%d\n", 
            SCIP_VERSION/100, (SCIP_VERSION/10) % 10, SCIP_VERSION % 10);
#else
         SCIPmessageFPrintInfo(file, "# SCIP version %d.%d.%d.%d\n", 
            SCIP_VERSION/100, (SCIP_VERSION/10) % 10, SCIP_VERSION % 10, SCIP_SUBVERSION);
#endif
      
      SCIPmessageFPrintInfo(file, "\n");
   }
   
   /* write the parameters to the file */
   for( i = 0; i < paramset->nparams; ++i )
   {
      SCIP_CALL( paramWrite(paramset->params[i], file, comments, onlychanged) );
   }

   /* close output file */
   if( filename != NULL )
   {
      assert(file != NULL);
      fclose(file);
   }

   return SCIP_OKAY;
}

/** installs default values for all parameters */
SCIP_RETCODE SCIPparamsetSetToDefaults(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip                /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   )
{
   int i;

   /* set all parameters to their default values */
   for( i = 0; i < paramset->nparams; ++i )
   {
      SCIP_CALL( SCIPparamSetToDefault(paramset->params[i], scip) );
   }

   return SCIP_OKAY;
}

/** installs default value for a single parameter */
SCIP_RETCODE SCIPparamsetSetToDefault(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   const char*           paramname           /**< name of the parameter */
   )
{
   SCIP_PARAM* param;

   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)paramname);
   
   if( param == NULL )
   {
      SCIPerrorMessage("parameter <%s> unknown\n", paramname);
      return SCIP_PARAMETERUNKNOWN;
   }

   SCIP_CALL( SCIPparamSetToDefault(param, scip) );
   
   return SCIP_OKAY;
}

/** resets all parameters which start with given prefix to their default value */
static
SCIP_RETCODE paramsetSetDefault(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           prefix              /**< prefix of the parameter name */
   )
{
   SCIP_PARAM* param;
   int nparams;
   size_t prefixlen;
   int p;

   nparams = paramset->nparams;
   prefixlen = strlen(prefix);

   for( p = 0; p < nparams; ++p )
   {
      param = paramset->params[p];
      /* check if the prefix matches */
      if( strncmp(param->name, prefix, prefixlen) == 0 )
      {
         SCIP_CALL( SCIPparamSetToDefault(param, scip) );
      }
   }

   return SCIP_OKAY;
}

/** resets parameters changed by other SCIPparamsetSetToHeuristicsXyz functions to their default values */
static
SCIP_RETCODE paramsetSetHeuristicsDefault(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{  /*lint --e{715}*/
   /* reset all parameter which start with "heuristics" in their name to default */
   SCIP_CALL( paramsetSetDefault(paramset, scip, "heuristics") );

   return SCIP_OKAY;
}

/** sets heuristics to aggressive */
static
SCIP_RETCODE paramsetSetHeuristicsAggressive(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_HEUR** heurs;
   SCIP_PARAM* param;
   char paramname[SCIP_MAXSTRLEN];
   int nheurs;
   int i;

   heurs = SCIPgetHeurs(scip);
   nheurs = SCIPgetNHeurs(scip);

   for( i = 0; i < nheurs; ++i )
   {
      const char* heurname;
      heurname = SCIPheurGetName(heurs[i]);

      /* get frequency parameter of heuristic */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/freq", heurname);
      param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)paramname);
      
      if( param != NULL )
      {
         int deffreq;
         int newfreq;

         assert(SCIPparamGetType(param) == SCIP_PARAMTYPE_INT);
         deffreq = SCIPparamGetIntDefault(param);

         /* change frequency to half of the default value, if it is > 0, otherwise set to 20 */
         if( deffreq == -1 || deffreq == 0 )
         {
            newfreq = 20;
         }
         else
         {
            newfreq = (int) SCIPceil(scip, deffreq/2.0);
            newfreq = MAX(newfreq, 1);
         }

         SCIP_CALL( paramSetInt(scip, paramset, paramname, newfreq, quiet) );
      }

      /* LP iteration limits only get increased for heuristics which are activated by default */
      if( SCIPparamGetIntDefault(param) > -1 )
      {
         /* construct (possible) parameter name for LP iteration offset */
         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/maxlpiterofs", heurname);
         param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)paramname);

         if( param != NULL && SCIPparamGetType(param) == SCIP_PARAMTYPE_INT )
         {
            /* set LP iteration offset to 1.5 time the current value */
            SCIP_CALL( paramSetInt(scip, paramset, paramname, (int) (1.5 * SCIPparamGetIntDefault(param)), quiet) );
         }
      
         /* construct (possible) parameter name for LP iteration quotient parameter */
         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/maxlpiterquot", heurname);
         param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)paramname);

         if( param != NULL && SCIPparamGetType(param) == SCIP_PARAMTYPE_REAL )
         {
            /* set LP iteration quotient to 1.5 time the current value */
            SCIP_CALL( paramSetReal(scip, paramset, paramname, 1.5*SCIPparamGetRealDefault(param), quiet) );
         }
      }
   }  

   /* set specific parameters for RENS heuristic */
   SCIP_CALL( paramSetLongint(scip, paramset, "heuristics/rens/nodesofs", (SCIP_Longint)2000, quiet) );
   SCIP_CALL( paramSetReal(scip, paramset, "heuristics/rens/minfixingrate", 0.3, quiet) );
   
   /* set specific parameters for Crossover heuristic */
   SCIP_CALL( paramSetLongint(scip, paramset, "heuristics/crossover/nwaitingnodes", (SCIP_Longint)20, quiet) );
   SCIP_CALL( paramSetBool(scip, paramset, "heuristics/crossover/dontwaitatroot", TRUE, quiet) );
   SCIP_CALL( paramSetReal(scip, paramset, "heuristics/crossover/nodesquot", 0.15, quiet) );
   SCIP_CALL( paramSetReal(scip, paramset, "heuristics/crossover/minfixingrate", 0.5, quiet) );

   return SCIP_OKAY;
}

/** sets heuristics to fast */
static
SCIP_RETCODE paramsetSetHeuristicsFast(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   int i;
#  define NEXPENSIVEHEURFREQS 13
   static const char* const expensiveheurfreqs[NEXPENSIVEHEURFREQS] = {
      "heuristics/coefdiving/freq",
      "heuristics/crossover/freq",
      "heuristics/feaspump/freq",
      "heuristics/fracdiving/freq",
      "heuristics/guideddiving/freq",
      "heuristics/linesearchdiving/freq",
      "heuristics/subnlp/freq",
      "heuristics/objpscostdiving/freq",
      "heuristics/pscostdiving/freq",
      "heuristics/rens/freq",
      "heuristics/rootsoldiving/freq",
      "heuristics/undercover/freq",
      "heuristics/veclendiving/freq"
   };

   SCIP_CALL( paramsetSetHeuristicsDefault(paramset, scip, quiet) );

   /* explicitly turn off expensive heuristics, if included */
   for( i = 0; i < NEXPENSIVEHEURFREQS; ++i )
      if( SCIPhashtableRetrieve(paramset->hashtable, (void*)expensiveheurfreqs[i]) != NULL )
      {
         SCIP_CALL( paramSetInt(scip, paramset, expensiveheurfreqs[i], -1, quiet) );
      }
   
   return SCIP_OKAY;
}

/** turns all heuristics off */
static
SCIP_RETCODE paramsetSetHeuristicsOff(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_HEUR** heurs;
   char paramname[SCIP_MAXSTRLEN];
   int nheurs;
   int i;

   heurs = SCIPgetHeurs(scip);
   nheurs = SCIPgetNHeurs(scip);

   /* no need to reset heuristic parameters, since we switch off of them off anyway */

   for( i = 0; i < nheurs; ++i )
   {
      const char* heurname;
      heurname = SCIPheurGetName(heurs[i]);

      /* get frequency parameter of heuristic */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/freq", heurname);

      SCIP_CALL( paramSetInt(scip, paramset, paramname, -1, quiet) );
   }
   
   return SCIP_OKAY;
}

/** resets all parameters that start "presolving" in their name to their default value; additionally set the parameters
 *  which  might be previous changed by the methods SCIPparamsetSetToPresolving{Off,Fast,Aggressive} to their default value
 */
static
SCIP_RETCODE paramsetSetPresolvingDefault(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{  /*lint --e{715}*/
   SCIP_CONSHDLR** conshdlrs;
   SCIP_PARAM* param;
   char paramname[SCIP_MAXSTRLEN];
   int nconshdlrs;
   int i;

   /* reset all parameter which start with "presolving" in their name to default */
   SCIP_CALL( paramsetSetDefault(paramset, scip, "presolving") );
   
   conshdlrs = SCIPgetConshdlrs(scip);
   nconshdlrs = SCIPgetNConshdlrs(scip);

   /* reset presolving settings for each individual constraint handler */
   for( i = 0; i < nconshdlrs; ++i )
   {
      const char* conshdlrname;
      conshdlrname = SCIPconshdlrGetName(conshdlrs[i]);

      /* reset maxprerounds parameter of presolvers */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/maxprerounds", conshdlrname);
      SCIP_CALL( SCIPparamsetSetToDefault(paramset, scip, paramname) );

      /* reset presolpairwise parameter of constraint handler */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/presolpairwise", conshdlrname);
      param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)paramname);
      if( param != NULL && SCIPparamGetType(param) == SCIP_PARAMTYPE_BOOL )
      {
         SCIP_CALL( SCIPparamsetSetToDefault(paramset, scip, paramname) );
      }
   }

   /* explicitly reset parameters of knapsack constraint handler, if the constraint handler is included */
   SCIP_CALL( SCIPparamsetSetToDefault(paramset, scip, "constraints/knapsack/disaggregation") );
   SCIP_CALL( SCIPparamsetSetToDefault(paramset, scip, "constraints/knapsack/simplifyinequalities") );

   /* explicitly reset parameters of linear constraint handler, if the constraint handler is included */
   SCIP_CALL( SCIPparamsetSetToDefault(paramset, scip, "constraints/linear/simplifyinequalities") );

   return SCIP_OKAY;
}

/** sets presolving to aggressive */
static
SCIP_RETCODE paramsetSetPresolvingAggressive(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_PARAM* param;
   char paramname[SCIP_MAXSTRLEN];
   
   /* reset previous changes on presolving parameters */
   SCIP_CALL( paramsetSetPresolvingDefault(paramset, scip, quiet) );

   /* explicitly change restart parameters */
   SCIP_CALL( paramSetReal(scip, paramset, "presolving/restartfac", 0.03, quiet) );
   SCIP_CALL( paramSetReal(scip, paramset, "presolving/restartminred", 0.06, quiet) );

   /* explicitly change parameters of linear constraint handler, if included */
   SCIP_CALL( paramSetBool(scip, paramset, "constraints/linear/simplifyinequalities", TRUE, quiet) );
   
   /* explicitly change parameters of knapsack constraint handler, if included */
   SCIP_CALL( paramSetBool(scip, paramset, "constraints/knapsack/simplifyinequalities", TRUE, quiet) );
   
   /* explicitly change parameters of presolver boundshift */
   SCIP_CALL( paramSetInt(scip, paramset, "presolving/boundshift/maxrounds", -1, quiet) );
 
   /* explicitly change parameters of probing */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "propagating/probing/maxuseless");
   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)paramname);
   if( param != NULL )
   {
      int defvalue;
    
      assert(SCIPparamGetType(param) == SCIP_PARAMTYPE_INT);
      defvalue = SCIPparamGetIntDefault(param);
    
      SCIP_CALL( paramSetInt(scip, paramset, paramname, (int) (1.5 * defvalue), quiet) );
   }
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "propagating/probing/maxtotaluseless");
   param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)paramname);
   if( param != NULL )
   {
      int defvalue;
    
      assert(SCIPparamGetType(param) == SCIP_PARAMTYPE_INT);
      defvalue = SCIPparamGetIntDefault(param);
    
      SCIP_CALL( paramSetInt(scip, paramset, paramname, (int) (1.5 * defvalue), quiet) );
   }
 
   return SCIP_OKAY;
}

/** sets presolving to fast */
static
SCIP_RETCODE paramsetSetPresolvingFast(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_CONSHDLR** conshdlrs;
   SCIP_PARAM* param;
   char paramname[SCIP_MAXSTRLEN];
   int nconshdlrs;
   int i;

   /* reset previous changes on presolving parameters */
   SCIP_CALL( paramsetSetPresolvingDefault(paramset, scip, quiet) );

   conshdlrs = SCIPgetConshdlrs(scip);
   nconshdlrs = SCIPgetNConshdlrs(scip);

   /* turn off pairwise comparison for each constraint handler that has this feature */
   for( i = 0; i < nconshdlrs; ++i )
   {
      const char* conshdlrname;
      conshdlrname = SCIPconshdlrGetName(conshdlrs[i]);

      /* get presolpairwise parameter of constraint handler */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/presolpairwise", conshdlrname);
      param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)paramname);

      if( param != NULL && SCIPparamGetType(param) == SCIP_PARAMTYPE_BOOL )
      {
         SCIP_CALL( paramSetBool(scip, paramset, paramname, FALSE, quiet) );
      }
   }

   /* explicitly turn off restarts */
   SCIP_CALL( paramSetInt(scip, paramset, "presolving/maxrestarts", 0, quiet) );

   /* turn off probing */
   SCIP_CALL( paramSetInt(scip, paramset, "propagating/probing/maxprerounds", 0, quiet) );

   /* explicitly change parameters of knapsack constraint handler, if the constraint handler is included */
   SCIP_CALL( paramSetBool(scip, paramset, "constraints/knapsack/disaggregation", FALSE, quiet) );
   
   return SCIP_OKAY;
}

/** turns all presolving off */
static
SCIP_RETCODE paramsetSetPresolvingOff(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_PRESOL** presols;
   SCIP_CONSHDLR** conshdlrs;
   char paramname[SCIP_MAXSTRLEN];
   int npresols;
   int nconshdlrs;
   int i;

   /* no need to reset to default presolve settings first, since we disable all presolve's anyway, correct? */

   presols = SCIPgetPresols(scip);
   npresols = SCIPgetNPresols(scip);
   
   /* turn each individual presolver off */
   for( i = 0; i < npresols; ++i )
   {
      const char* presolname;
      presolname = SCIPpresolGetName(presols[i]);
      
      /* get maxrounds parameter of presolvers */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "presolving/%s/maxrounds", presolname);
      
      SCIP_CALL( paramSetInt(scip, paramset, paramname, 0, quiet) );
   }

   conshdlrs = SCIPgetConshdlrs(scip);
   nconshdlrs = SCIPgetNConshdlrs(scip);

   /* turn off presolving for each individual constraint handler */
   for( i = 0; i < nconshdlrs; ++i )
   {
      const char* conshdlrname;
      conshdlrname = SCIPconshdlrGetName(conshdlrs[i]);
      
      /* get maxprerounds parameter of presolvers */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/maxprerounds", conshdlrname);
      
      SCIP_CALL( paramSetInt(scip, paramset, paramname, 0, quiet) );
   }
   
   /* explicitly turn off restarts */
   SCIP_CALL( paramSetInt(scip, paramset, "presolving/maxrestarts", 0, quiet) );

   return SCIP_OKAY;
}

/** reset parameters that may have been changed by other SCIPparamsetSetToSeparatingXyz to their default values */
static
SCIP_RETCODE paramsetSetSeparatingDefault(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{  /*lint --e{715}*/
   SCIP_CONSHDLR** conshdlrs;
   char paramname[SCIP_MAXSTRLEN];
   int nconshdlrs;
   int i;

   /* reset all parameter which start with "separating" in their name to default */
   SCIP_CALL( paramsetSetDefault(paramset, scip, "separating") );

   conshdlrs = SCIPgetConshdlrs(scip);
   nconshdlrs = SCIPgetNConshdlrs(scip);

   /* reset each individual constraint handler separation settings */
   for( i = 0; i < nconshdlrs; ++i )
   {
      const char* conshdlrname;
      conshdlrname = SCIPconshdlrGetName(conshdlrs[i]);

      /* reset separation frequency parameter of constraint handler, if available */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/sepafreq", conshdlrname);
      SCIP_CALL( SCIPparamsetSetToDefault(paramset, scip, paramname) );
      
      /* reset maximal separated cuts in root node of constraint handler, if available */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/maxsepacutsroot", conshdlrname);
      if( SCIPhashtableRetrieve(paramset->hashtable, (void*)paramname) != NULL )
      {
         SCIP_CALL( SCIPparamsetSetToDefault(paramset, scip, paramname) );
      }
   }

   SCIP_CALL( SCIPparamsetSetToDefault(paramset, scip, "constraints/linear/separateall") );
   
   return SCIP_OKAY;
}

/** sets separating to aggressive */
static
SCIP_RETCODE paramsetSetSeparatingAggressive(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_CONSHDLR** conshdlrs;
   SCIP_SEPA** sepas;
   SCIP_PARAM* param;
   char paramname[SCIP_MAXSTRLEN];
   int nconshdlrs;
   int nsepas;
   int i;

   sepas = SCIPgetSepas(scip);
   nsepas = SCIPgetNSepas(scip);

   /* set all separating parameters to default values */
   SCIP_CALL( paramsetSetSeparatingDefault(paramset, scip, quiet) );
   
   /* set separating parameters of all separators */
   for( i = 0; i < nsepas; ++i )
   {
      const char* sepaname;
      sepaname = SCIPsepaGetName(sepas[i]);

      /* intobj separator should stay disabled */
      if( strcmp(sepaname, "intobj") == 0 )
         continue;
         
      /* get frequency parameter of separator */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "separating/%s/freq", sepaname);
      param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)paramname);
      
      if( param != NULL )
      {
         int deffreq;
         int newfreq;

         assert(SCIPparamGetType(param) == SCIP_PARAMTYPE_INT);
         deffreq = SCIPparamGetIntDefault(param);
         
         /* for enabled separators, change frequency to at least every 20th depths and 
          * enable disabled separators
          */
         if( deffreq == -1 )
            newfreq = 0;
         else if( deffreq == 0 )
            newfreq = 20;
         else
            newfreq = MIN(deffreq, 20);
         
         SCIP_CALL( paramSetInt(scip, paramset, paramname, newfreq, quiet) );
      }

      /* get maximum number of rounds in root node */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "separating/%s/maxroundsroot", sepaname);
      param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)paramname);
      
      if( param != NULL )
      {
         int defrounds;

         assert(SCIPparamGetType(param) == SCIP_PARAMTYPE_INT);
         defrounds = SCIPparamGetIntDefault(param);
         
         /* increase the maximum number of rounds in the root node by factor of 1.5 */
         SCIP_CALL( paramSetInt(scip, paramset, paramname, (int) (1.5 * defrounds), quiet) );
      }

      /* get maximum number of cuts per separation in root node */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "separating/%s/maxsepacutsroot", sepaname);
      param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)paramname);
      
      if( param != NULL )
      {
         int defnumber;

         assert(SCIPparamGetType(param) == SCIP_PARAMTYPE_INT);
         defnumber = SCIPparamGetIntDefault(param);
         
         /* increase the maximum number of cut per separation rounds in the root node by factor of 2 */
         SCIP_CALL( paramSetInt(scip, paramset, paramname, 2*defnumber, quiet) );
      }
   }

   conshdlrs = SCIPgetConshdlrs(scip);
   nconshdlrs = SCIPgetNConshdlrs(scip);

   /* set separating parameters of all constraint handlers */
   for( i = 0; i < nconshdlrs; ++i )
   {
      const char* conshdlrname;
      conshdlrname = SCIPconshdlrGetName(conshdlrs[i]);

      /* get separating frequency parameter of constraint handler */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/sepafreq", conshdlrname);
      param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)paramname);
      
      if( param != NULL )
      {
         int deffreq;
         int newfreq;

         assert(SCIPparamGetType(param) == SCIP_PARAMTYPE_INT);
         deffreq = SCIPparamGetIntDefault(param);
         
         /* for constraint handlers with enabled separation, change frequency to at least every 10th depths and 
          * enable disabled separation routines 
          */
         if( deffreq == -1 )
            newfreq = 0;
         else if( deffreq == 0 )
            newfreq = 10;
         else
            newfreq = MIN(deffreq, 10);
         
         SCIP_CALL( paramSetInt(scip, paramset, paramname, newfreq, quiet) );
      }

      /* get maximal separated cuts in root node  of constraint handler */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/maxsepacutsroot", conshdlrname);
      param = (SCIP_PARAM*)SCIPhashtableRetrieve(paramset->hashtable, (void*)paramname);
      
      if( param != NULL )
      {
         int defnumber;
         
         assert(SCIPparamGetType(param) == SCIP_PARAMTYPE_INT);
         defnumber = SCIPparamGetIntDefault(param);
         
         /* change maximal cuts in root node to at least 500 */
         SCIP_CALL( paramSetInt(scip, paramset, paramname, MAX(defnumber, 500), quiet) );
      }
   }

   /* explicitly change general separating parameters */
   SCIP_CALL( paramSetReal(scip, paramset, "separating/minorthoroot", 0.1, quiet) );
   SCIP_CALL( paramSetInt(scip, paramset, "separating/maxroundsrootsubrun", 5, quiet) );
   SCIP_CALL( paramSetInt(scip, paramset, "separating/maxaddrounds", 5, quiet) );
   SCIP_CALL( paramSetInt(scip, paramset, "separating/maxcutsroot", 5000, quiet) );
   SCIP_CALL( paramSetInt(scip, paramset, "separating/poolfreq", 10, quiet) );

   /* explicitly change a separating parameter of the linear constraint handler */
   SCIP_CALL( paramSetBool(scip, paramset, "constraints/linear/separateall", TRUE, quiet) );
   
   /* explicitly change a separating parameter of cmir separator */
   SCIP_CALL( paramSetInt(scip, paramset, "separating/cmir/maxfailsroot", 200, quiet) );

   /* explicitly change separating parameters of mcf separator */
   SCIP_CALL( paramSetInt(scip, paramset, "separating/mcf/maxtestdelta", -1, quiet) );
   SCIP_CALL( paramSetBool(scip, paramset, "separating/mcf/trynegscaling", TRUE, quiet) );

   return SCIP_OKAY;
}

/** sets separating to fast */
static
SCIP_RETCODE paramsetSetSeparatingFast(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   /* reset previous changes on separating parameters */
   SCIP_CALL( paramsetSetSeparatingDefault(paramset, scip, quiet) );

   /* explicitly decrease maxbounddist */
   SCIP_CALL( paramSetReal(scip, paramset, "separating/maxbounddist", 0.0, quiet) );

   /* explicitly turn off expensive separators */
   SCIP_CALL( paramSetInt(scip, paramset, "constraints/and/sepafreq", 0, quiet) );
   SCIP_CALL( paramSetInt(scip, paramset, "separating/cmir/maxroundsroot", 5, quiet) );
   SCIP_CALL( paramSetInt(scip, paramset, "separating/cmir/maxtriesroot", 100, quiet) );
   SCIP_CALL( paramSetInt(scip, paramset, "separating/cmir/maxaggrsroot", 3, quiet) );
   SCIP_CALL( paramSetInt(scip, paramset, "separating/cmir/maxsepacutsroot", 200, quiet) );
   SCIP_CALL( paramSetInt(scip, paramset, "separating/flowcover/freq", -1, quiet) );
   SCIP_CALL( paramSetInt(scip, paramset, "separating/gomory/maxroundsroot", 20, quiet) );
   SCIP_CALL( paramSetInt(scip, paramset, "separating/gomory/maxsepacutsroot", 200, quiet) );
   SCIP_CALL( paramSetInt(scip, paramset, "separating/mcf/freq", -1, quiet) );
   SCIP_CALL( paramSetReal(scip, paramset, "separating/redcost/maxbounddist", 0.0, quiet) );
   SCIP_CALL( paramSetInt(scip, paramset, "separating/strongcg/maxroundsroot", 10, quiet) );
   SCIP_CALL( paramSetInt(scip, paramset, "separating/strongcg/maxsepacutsroot", 200, quiet) );
   
   return SCIP_OKAY;
}

/** turns all cuts off */
static
SCIP_RETCODE paramsetSetSeparatingOff(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_SEPA** sepas;
   SCIP_CONSHDLR** conshdlrs;
   char paramname[SCIP_MAXSTRLEN];
   int nsepas;
   int nconshdlrs;
   int i;

   sepas = SCIPgetSepas(scip);
   nsepas = SCIPgetNSepas(scip);
   
   /* no need to reset separation parameters to default here, since we disable all separators anyway */

   /* turn each individual separator off */
   for( i = 0; i < nsepas; ++i )
   {
      const char* sepaname;
      sepaname = SCIPsepaGetName(sepas[i]);
      
      /* get frequency parameter of separator */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "separating/%s/freq", sepaname);
      SCIP_CALL( paramSetInt(scip, paramset, paramname, -1, quiet) );
   }

   conshdlrs = SCIPgetConshdlrs(scip);
   nconshdlrs = SCIPgetNConshdlrs(scip);

   /* turn off separation for each individual constraint handler */
   for( i = 0; i < nconshdlrs; ++i )
   {
      const char* conshdlrname;
      conshdlrname = SCIPconshdlrGetName(conshdlrs[i]);
      
      /* get separation frequency parameter of constraint handler */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/sepafreq", conshdlrname);
      SCIP_CALL( paramSetInt(scip, paramset, paramname, -1, quiet) );
   }
   
   return SCIP_OKAY;
}

/** sets parameters to 
 *  - SCIP_PARAMSETTING_DEFAULT to use default values (see also SCIPparamsetSetToDefault())
 *  - SCIP_PARAMSETTING_COUNTER to get feasible and "fast" counting process
 *  - SCIP_PARAMSETTING_CPSOLVER to get CP like search (e.g. no LP relaxation)
 *  - SCIP_PARAMSETTING_EASYCIP to solve easy problems fast
 *  - SCIP_PARAMSETTING_FEASIBILITY to detect feasibility fast 
 *  - SCIP_PARAMSETTING_HARDLP to be capable to handle hard LPs
 *  - SCIP_PARAMSETTING_OPTIMALITY to prove optimality fast
 */
SCIP_RETCODE SCIPparamsetSetEmphasis(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMEMPHASIS    paramemphasis,      /**< parameter emphasis */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   /* reset all parameter to default */
   SCIP_CALL( SCIPparamsetSetToDefaults(paramset, scip) );
      
   switch( paramemphasis )
   {
   case SCIP_PARAMSETTING_DEFAULT:
      /* the default values are already set */
      break;

   case SCIP_PARAMEMPHASIS_COUNTER:
      /* avoid logicor upgrade since the logicor constraint handler does not perform full propagation */ 
      SCIP_CALL( paramSetBool(scip, paramset, "constraints/linear/upgrade/logicor", FALSE, quiet) );
   
      /* set priority for inference branching to highest possible value */
      SCIP_CALL( paramSetInt(scip, paramset, "branching/inference/priority", INT_MAX/4, quiet) );
 
      /* set priority for depth first search to highest possible value */
      SCIP_CALL( paramSetInt(scip, paramset, "nodeselection/dfs/stdpriority", INT_MAX/4, quiet) );
      
      /* avoid that the ZIMPL reader transforms the problem before the problem is generated */
      SCIP_CALL( paramSetBool(scip, paramset, "reading/zplreader/usestartsol", FALSE, quiet) );

      /* turn off all heuristics */
      SCIP_CALL( paramsetSetHeuristicsOff(paramset, scip, quiet) );

      /* turn off all separation */
      SCIP_CALL( paramsetSetSeparatingOff(paramset, scip, quiet) );
      
      /* turn off restart */
      SCIP_CALL( paramSetInt(scip, paramset, "presolving/maxrestarts", 0, quiet) );

      /* unlimited propagation round in any branch and bound node */
      SCIP_CALL( paramSetInt(scip, paramset, "propagating/maxrounds", -1, quiet) );
      SCIP_CALL( paramSetInt(scip, paramset, "propagating/maxroundsroot", -1, quiet) );

      /* adjust conflict analysis for depth first search */
      SCIP_CALL( paramSetInt(scip, paramset, "conflict/fuiplevels", 1, quiet) );        
      SCIP_CALL( paramSetBool(scip, paramset, "conflict/dynamic", FALSE, quiet) );
   
      /* prefer binary variables for branching */
      SCIP_CALL( paramSetBool(scip, paramset, "branching/preferbinary", TRUE, quiet) );

      /* turn on aggressive constraint aging */ 
      SCIP_CALL( paramSetInt(scip, paramset, "constraints/agelimit", 1, quiet) );       
      break;

   case SCIP_PARAMEMPHASIS_CPSOLVER:
      /* use only first unique implication point */
      SCIP_CALL( paramSetInt(scip, paramset, "conflict/fuiplevels", 1, quiet) );

      /* conflict constraints should be subject to aging to reduce the number of stored conflicts */
      SCIP_CALL( paramSetBool(scip, paramset, "conflict/dynamic", FALSE, quiet) );

      /* perform a restart after 100 conflict analysis calls were successfully */
      SCIP_CALL( paramSetInt(scip, paramset, "conflict/restartnum", 100, quiet) );

      /* do not check pseudo solution (for performance reasons) */
      SCIP_CALL( paramSetBool(scip, paramset, "constraints/disableenfops", TRUE, quiet) );

      /* turn of LP relaxation */
      SCIP_CALL( paramSetInt(scip, paramset, "lp/solvefreq", -1, quiet) );

      /* set priority for depth first search to highest possible value */
      SCIP_CALL( paramSetInt(scip, paramset, "nodeselection/dfs/stdpriority", INT_MAX/4, quiet) );

      /* turn off pseudo cost propagation  */
      SCIP_CALL( paramSetInt(scip, paramset, "propagating/pseudoobj/freq", -1, quiet) );

      break;

   case SCIP_PARAMEMPHASIS_EASYCIP:
      /* set heuristics to fast, to avoid spending to much time for involved heuristics */
      SCIP_CALL( paramsetSetHeuristicsFast(paramset, scip, quiet) );

      /* set presolving to fast, to avoid spending to much time for involved presolving */
      SCIP_CALL( paramsetSetPresolvingFast(paramset, scip, quiet) );

      /* set separating to fast, to avoid spending to much time for involved separators */
      SCIP_CALL( paramsetSetSeparatingFast(paramset, scip, quiet) );
      
      break;

   case SCIP_PARAMEMPHASIS_FEASIBILITY:
      /* set heuristics aggressive */
      SCIP_CALL( paramsetSetHeuristicsAggressive(paramset, scip, quiet) );
      
      /* reduce the amount of separation rounds and disable most expensive separators */
      SCIP_CALL( paramSetInt(scip, paramset, "separating/maxrounds", 1, quiet) );
      SCIP_CALL( paramSetInt(scip, paramset, "separating/maxroundsroot", 5, quiet) );
      SCIP_CALL( paramSetInt(scip, paramset, "separating/cmir/freq", -1, quiet) );
      SCIP_CALL( paramSetInt(scip, paramset, "separating/flowcover/freq", -1, quiet) );
      SCIP_CALL( paramSetInt(scip, paramset, "separating/mcf/freq", -1, quiet) );
   
      /* set priority for node selection "restartdfs" to be higher as the current used one */
      SCIP_CALL( paramSetInt(scip, paramset, "nodeselection/restartdfs/stdpriority", INT_MAX/4, quiet) );
      break;

   case SCIP_PARAMEMPHASIS_HARDLP:
      /* set heuristics to fast, to avoid heuristics which solve also an LP */
      SCIP_CALL( paramsetSetHeuristicsFast(paramset, scip, quiet) );

      /* set presolving to fast, to avoid spending to much time for involved presolving */
      SCIP_CALL( paramsetSetPresolvingFast(paramset, scip, quiet) );
      
      /* reduce the amount of strong branching */
      SCIP_CALL( paramSetReal(scip, paramset, "branching/relpscost/maxreliable", 1.0, quiet) );
      SCIP_CALL( paramSetInt(scip, paramset, "branching/relpscost/inititer", 10, quiet) );
      
      /* reduce the amount of separation rounds */
      SCIP_CALL( paramSetInt(scip, paramset, "separating/maxrounds", 1, quiet) );
      SCIP_CALL( paramSetInt(scip, paramset, "separating/maxroundsroot", 5, quiet) );

      break;

   case SCIP_PARAMEMPHASIS_OPTIMALITY:
      /* set cuts aggressive */
      SCIP_CALL( paramsetSetSeparatingAggressive(paramset, scip, quiet) );

      /* increase the amount of strong branching */
      SCIP_CALL( paramSetInt(scip, paramset, "branching/fullstrong/maxdepth", 10, quiet) );
      SCIP_CALL( paramSetReal(scip, paramset, "branching/fullstrong/maxbounddist", 0.0, quiet) );
      SCIP_CALL( paramSetReal(scip, paramset, "branching/relpscost/sbiterquot", 1.0, quiet) );
      SCIP_CALL( paramSetInt(scip, paramset, "branching/relpscost/sbiterofs", 1000000, quiet) );

      break;

   default:
      SCIPerrorMessage("the parameter setting <%d> is not allowed for emphasis call\n", paramemphasis);
      return SCIP_INVALIDCALL;
   }  
   return SCIP_OKAY;
}

/** sets parameters to deactivate separators and heuristics that use auxiliary SCIP instances; should be called for
 *  auxiliary SCIP instances to avoid recursion
 */
SCIP_RETCODE SCIPparamsetSetToSubscipsOff(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_HEUR** heurs;
   SCIP_SEPA** sepas;

   char paramname[SCIP_MAXSTRLEN];

   int nheurs;
   int nsepas;
   int i;

   heurs = SCIPgetHeurs(scip);
   nheurs = SCIPgetNHeurs(scip);

   /* disable all heuristics that use auxiliary SCIP instances */
   for( i = 0; i < nheurs; ++i )
   {
      if( SCIPheurUsesSubscip(heurs[i]) )
      {
         const char* heurname;
         heurname = SCIPheurGetName(heurs[i]);

         /* get frequency parameter of heuristic */
         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/freq", heurname);

         SCIP_CALL( paramSetInt(scip, paramset, paramname, -1, quiet) );
      }
   }

   sepas = SCIPgetSepas(scip);
   nsepas = SCIPgetNSepas(scip);

   /* disable all separators that use auxiliary SCIP instances */
   for( i = 0; i < nsepas; ++i )
   {
      if( SCIPsepaUsesSubscip(sepas[i]) )
      {
         const char* sepaname;
         sepaname = SCIPsepaGetName(sepas[i]);
      
         /* get frequency parameter of separator */
         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "separating/%s/freq", sepaname);
         SCIP_CALL( paramSetInt(scip, paramset, paramname, -1, quiet) );
      }
   }
   
   return SCIP_OKAY;
}

/** sets heuristic parameters values to 
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all heuristic parameters 
 *  - SCIP_PARAMSETTING_FAST such that the time spend for heuristic is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the heuristic are called more aggregative
 *  - SCIP_PARAMSETTING_OFF which turn off all heuristics
 */
SCIP_RETCODE SCIPparamsetSetHeuristics(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   switch( paramsetting )
   {
   case SCIP_PARAMSETTING_DEFAULT:
      SCIP_CALL( paramsetSetHeuristicsDefault(paramset, scip, quiet) );
      break;
   case SCIP_PARAMSETTING_OFF:
      SCIP_CALL( paramsetSetHeuristicsOff(paramset, scip, quiet) );
      break;
   case SCIP_PARAMSETTING_FAST:
      SCIP_CALL( paramsetSetHeuristicsFast(paramset, scip, quiet) );
      break;
   case SCIP_PARAMSETTING_AGGRESSIVE:
      SCIP_CALL( paramsetSetHeuristicsAggressive(paramset, scip, quiet) );
      break;
   default:
      SCIPerrorMessage("the parameter setting <%d> is not allowed for heuristics\n", paramsetting);
      return SCIP_INVALIDCALL;
   }
   
   return SCIP_OKAY;
}

/** sets presolving parameters to 
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all presolving parameters 
 *  - SCIP_PARAMSETTING_FAST such that the time spend for presolving is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the presolving is more aggregative
 *  - SCIP_PARAMSETTING_OFF which turn off all presolving
 */
SCIP_RETCODE SCIPparamsetSetPresolving(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   switch( paramsetting )
   {
   case SCIP_PARAMSETTING_DEFAULT:
      SCIP_CALL( paramsetSetPresolvingDefault(paramset, scip, quiet) );
      break;
   case SCIP_PARAMSETTING_OFF:
      SCIP_CALL( paramsetSetPresolvingOff(paramset, scip, quiet) );
      break;
   case SCIP_PARAMSETTING_FAST:
      SCIP_CALL( paramsetSetPresolvingFast(paramset, scip, quiet) );
      break;
   case SCIP_PARAMSETTING_AGGRESSIVE:
      SCIP_CALL( paramsetSetPresolvingAggressive(paramset, scip, quiet) );
      break;
   default:
      SCIPerrorMessage("the parameter setting <%d> is not allowed for presolving\n", paramsetting);
      return SCIP_INVALIDCALL;
   }
   
   return SCIP_OKAY;
}

/** sets separating parameters to 
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all separating parameters 
 *  - SCIP_PARAMSETTING_FAST such that the time spend for separating is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the separating is done more aggregative
 *  - SCIP_PARAMSETTING_OFF which turn off all separating
 */
SCIP_RETCODE SCIPparamsetSetSeparating(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   switch( paramsetting )
   {
   case SCIP_PARAMSETTING_DEFAULT:
      SCIP_CALL( paramsetSetSeparatingDefault(paramset, scip, quiet) );
      break;
   case SCIP_PARAMSETTING_OFF:
      SCIP_CALL( paramsetSetSeparatingOff(paramset, scip, quiet) );
      break;
   case SCIP_PARAMSETTING_FAST:
      SCIP_CALL( paramsetSetSeparatingFast(paramset, scip, quiet) );
      break;
   case SCIP_PARAMSETTING_AGGRESSIVE:
      SCIP_CALL( paramsetSetSeparatingAggressive(paramset, scip, quiet) );
      break;
   default:
      SCIPerrorMessage("the parameter setting <%d> is not allowed for separating\n", paramsetting);
      return SCIP_INVALIDCALL;
   }
   
   return SCIP_OKAY;
}

/** returns the array of parameters */
SCIP_PARAM** SCIPparamsetGetParams(
   SCIP_PARAMSET*        paramset            /**< parameter set */
   )
{
   assert(paramset != NULL);

   return paramset->params;
}

/** returns the number of parameters in the parameter set */
int SCIPparamsetGetNParams(
   SCIP_PARAMSET*        paramset            /**< parameter set */
   )
{
   assert(paramset != NULL);

   return paramset->nparams;
}

/** copies all parameter values of the source parameter set to the corresponding parameters in the target set */
SCIP_RETCODE SCIPparamsetCopyParams(
   SCIP_PARAMSET*          sourceparamset,     /**< source parameter set */
   SCIP_PARAMSET*          targetparamset,     /**< target parameter set */
   SCIP*                   targetscip          /**< target SCIP data structure */
   )
{
   int i;

   assert(sourceparamset != NULL);
   assert(targetparamset != NULL);
   assert(sourceparamset != targetparamset);
   assert(targetscip != NULL);

   assert(sourceparamset->nparams == 0 || sourceparamset->params != NULL);
   assert(targetparamset->nparams == 0 || targetparamset->params != NULL);

   for( i = 0; i < sourceparamset->nparams; ++i )
   { 
      SCIP_PARAM* targetparam;
      const char* paramname;
      assert(sourceparamset->params[i] != NULL);     

      /* find parameter of same name in target scip */
      paramname = SCIPparamGetName(sourceparamset->params[i]);
      targetparam = (SCIP_PARAM*)SCIPhashtableRetrieve(targetparamset->hashtable, (void*)paramname);

      /* if a plugin was not copied, the parameter does not exist in the target SCIP */
      if( targetparam == NULL )
         continue;
      
      assert(SCIPparamGetType(sourceparamset->params[i]) == SCIPparamGetType(targetparam));
      
      /* set value of target parameter to value of source parameter */
      switch( SCIPparamGetType(sourceparamset->params[i]) )
      {
      case SCIP_PARAMTYPE_BOOL:
         SCIP_CALL( paramCopyBool(sourceparamset->params[i], targetparam, targetscip) );
         break;

      case SCIP_PARAMTYPE_INT:
         SCIP_CALL( paramCopyInt(sourceparamset->params[i], targetparam, targetscip) );
         break;

      case SCIP_PARAMTYPE_LONGINT:
         SCIP_CALL( paramCopyLongint(sourceparamset->params[i], targetparam, targetscip) );
         break;
         
      case SCIP_PARAMTYPE_REAL:
         SCIP_CALL( paramCopyReal(sourceparamset->params[i], targetparam, targetscip) );
         break;
         
      case SCIP_PARAMTYPE_CHAR:
         SCIP_CALL( paramCopyChar(sourceparamset->params[i], targetparam, targetscip) );
         break;
         
      case SCIP_PARAMTYPE_STRING:
         /* the vbc parameters are explicitly not copied to avoid that the vbc file of the original SCIP is over
          * written; to avoid that hard coded comparison, each parameter could get a Bool flag which tells if the value
          * of that parameter can be copied */
         if( strncmp(sourceparamset->params[i]->name, "vbc/", 4) != 0 )
         {
            SCIP_CALL( paramCopyString(sourceparamset->params[i], targetparam, targetscip) ); 
         }
         break;
         
      default:
         SCIPerrorMessage("unknown parameter type\n");
         return SCIP_INVALIDDATA;
      }
   }   
  
   return SCIP_OKAY;
}

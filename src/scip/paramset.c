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

/**@file   paramset.c
 * @brief  methods and datastructures for handling parameter settings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <math.h>

#include "memory.h"
#include "misc.h"
#include "message.h"
#include "paramset.h"



/** data for Bool parameters */
struct BoolParam
{
   Bool*            valueptr;           /**< pointer to store the current parameter value, or NULL */
   Bool             actvalue;           /**< stores the actual parameter value if it is not stored in *valueptr */
   Bool             defaultvalue;       /**< default value of the parameter */
};
typedef struct BoolParam BOOLPARAM;

/** data for int parameters */
struct IntParam
{
   int*             valueptr;           /**< pointer to store the current parameter value, or NULL */
   int              actvalue;           /**< stores the actual parameter value if it is not stored in *valueptr */
   int              defaultvalue;       /**< default value of the parameter */
};
typedef struct IntParam INTPARAM;

/** data for Longint parameters */
struct LongintParam
{
   Longint*         valueptr;           /**< pointer to store the current parameter value, or NULL */
   Longint          actvalue;           /**< stores the actual parameter value if it is not stored in *valueptr */
   Longint          defaultvalue;       /**< default value of the parameter */
};
typedef struct LongintParam LONGINTPARAM;

/** data for Real parameters */
struct RealParam
{
   Real*            valueptr;           /**< pointer to store the current parameter value, or NULL */
   Real             actvalue;           /**< stores the actual parameter value if it is not stored in *valueptr */
   Real             defaultvalue;       /**< default value of the parameter */
};
typedef struct RealParam REALPARAM;

/** data for char parameters */
struct CharParam
{
   char*            valueptr;           /**< pointer to store the current parameter value, or NULL */
   char             actvalue;           /**< stores the actual parameter value if it is not stored in *valueptr */
   char             defaultvalue;       /**< default value of the parameter */
};
typedef struct CharParam CHARPARAM;

/** data for char* parameters */
struct StringParam
{
   char**           valueptr;           /**< pointer to store the current parameter value, or NULL */
   char*            actvalue;           /**< stores the actual parameter value if it is not stored in *valueptr */
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
   unsigned int     paramtype:3;        /**< type of this parameter */
};

/** set of parameters */
struct ParamSet
{
   HASHTABLE*       hashtable;          /**< hash table to store the parameters */
   PARAM**          params;             /**< array with parameters */
   int              nparams;            /**< number of parameters */
   int              paramssize;         /**< size of params array */
};




/*
 * local methods
 */

/** hash key retrieval function for parameters */
static
DECL_HASHGETKEY(hashGetKeyParam)
{
   PARAM* param;

   param = (PARAM*)elem;
   assert(param != NULL);

   return param->name;
}

/** creates a Bool parameter, and sets its value to default */
static
RETCODE paramCreateBool(
   PARAM**          param,              /**< pointer to the parameter */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   Bool*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   Bool             defaultvalue        /**< default value of the parameter */
   )
{
   assert(param != NULL);
   assert(name != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, param) );
   
   ALLOC_OKAY( duplicateMemoryArray(&(*param)->name, name, strlen(name)+1) );
   (*param)->paramtype = SCIP_PARAMTYPE_BOOL;
   (*param)->data.boolparam.valueptr = valueptr;
   (*param)->data.boolparam.defaultvalue = defaultvalue;
   if( valueptr == NULL )
      (*param)->data.boolparam.actvalue = defaultvalue;
   else
      *valueptr = defaultvalue;

   return SCIP_OKAY;
}

/** creates a int parameter, and sets its value to default */
static
RETCODE paramCreateInt(
   PARAM**          param,              /**< pointer to the parameter */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   int*             valueptr,           /**< pointer to store the current parameter value, or NULL */
   int              defaultvalue        /**< default value of the parameter */
   )
{
   assert(param != NULL);
   assert(name != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, param) );
   
   ALLOC_OKAY( duplicateMemoryArray(&(*param)->name, name, strlen(name)+1) );
   (*param)->paramtype = SCIP_PARAMTYPE_INT;
   (*param)->data.intparam.valueptr = valueptr;
   (*param)->data.intparam.defaultvalue = defaultvalue;
   if( valueptr == NULL )
      (*param)->data.intparam.actvalue = defaultvalue;
   else
      *valueptr = defaultvalue;

   return SCIP_OKAY;
}

/** creates a Longint parameter, and sets its value to default */
static
RETCODE paramCreateLongint(
   PARAM**          param,              /**< pointer to the parameter */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   Longint*         valueptr,           /**< pointer to store the current parameter value, or NULL */
   Longint          defaultvalue        /**< default value of the parameter */
   )
{
   assert(param != NULL);
   assert(name != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, param) );
   
   ALLOC_OKAY( duplicateMemoryArray(&(*param)->name, name, strlen(name)+1) );
   (*param)->paramtype = SCIP_PARAMTYPE_LONGINT;
   (*param)->data.longintparam.valueptr = valueptr;
   (*param)->data.longintparam.defaultvalue = defaultvalue;
   if( valueptr == NULL )
      (*param)->data.longintparam.actvalue = defaultvalue;
   else
      *valueptr = defaultvalue;

   return SCIP_OKAY;
}

/** creates a Real parameter, and sets its value to default */
static
RETCODE paramCreateReal(
   PARAM**          param,              /**< pointer to the parameter */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   Real*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   Real             defaultvalue        /**< default value of the parameter */
   )
{
   assert(param != NULL);
   assert(name != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, param) );
   
   ALLOC_OKAY( duplicateMemoryArray(&(*param)->name, name, strlen(name)+1) );
   (*param)->paramtype = SCIP_PARAMTYPE_REAL;
   (*param)->data.realparam.valueptr = valueptr;
   (*param)->data.realparam.defaultvalue = defaultvalue;
   if( valueptr == NULL )
      (*param)->data.realparam.actvalue = defaultvalue;
   else
      *valueptr = defaultvalue;

   return SCIP_OKAY;
}

/** creates a char parameter, and sets its value to default */
static
RETCODE paramCreateChar(
   PARAM**          param,              /**< pointer to the parameter */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   char*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   char             defaultvalue        /**< default value of the parameter */
   )
{
   assert(param != NULL);
   assert(name != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, param) );
   
   ALLOC_OKAY( duplicateMemoryArray(&(*param)->name, name, strlen(name)+1) );
   (*param)->paramtype = SCIP_PARAMTYPE_CHAR;
   (*param)->data.charparam.valueptr = valueptr;
   (*param)->data.charparam.defaultvalue = defaultvalue;
   if( valueptr == NULL )
      (*param)->data.charparam.actvalue = defaultvalue;
   else
      *valueptr = defaultvalue;

   return SCIP_OKAY;
}

/** creates a string parameter, and sets its value to default */
static
RETCODE paramCreateString(
   PARAM**          param,              /**< pointer to the parameter */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   char**           valueptr,           /**< pointer to store the current parameter value, or NULL */
   const char*      defaultvalue        /**< default value of the parameter */
   )
{
   assert(param != NULL);
   assert(name != NULL);
   assert(valueptr == NULL || *valueptr == NULL);
   assert(defaultvalue != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, param) );
   
   ALLOC_OKAY( duplicateMemoryArray(&(*param)->name, name, strlen(name)+1) );
   (*param)->paramtype = SCIP_PARAMTYPE_STRING;
   (*param)->data.stringparam.valueptr = valueptr;
   ALLOC_OKAY( duplicateMemoryArray(&(*param)->data.stringparam.defaultvalue, defaultvalue, strlen(defaultvalue)+1) );
   if( valueptr == NULL )
   {
      ALLOC_OKAY( duplicateMemoryArray(&(*param)->data.stringparam.actvalue, defaultvalue, strlen(defaultvalue)+1) );
   }
   else
   {
      ALLOC_OKAY( duplicateMemoryArray(valueptr, defaultvalue, strlen(defaultvalue)+1) );
   }

   return SCIP_OKAY;
}

/** frees a single parameter */
static
void paramFree(
   PARAM**          param,              /**< pointer to the parameter */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(param != NULL);
   assert(*param != NULL);

   switch( (*param)->paramtype )
   {
   case SCIP_PARAMTYPE_STRING:
      freeMemoryArray(&(*param)->data.stringparam.defaultvalue);
      if( (*param)->data.stringparam.valueptr == NULL )
      {
         freeMemoryArray(&(*param)->data.stringparam.actvalue);
      }
      else
      {
         freeMemoryArray(&(*param)->data.stringparam.valueptr);
      }
      break;
   default:
      break;
   }

   freeMemoryArray(&(*param)->name);
   freeBlockMemory(memhdr, param);
}


/*
 * interface methods
 */

/** creates parameter set */
RETCODE SCIPparamsetCreate(
   PARAMSET**       paramset            /**< pointer to store the parameter set */
   )
{
   assert(paramset != NULL);

   ALLOC_OKAY( allocMemory(paramset) );

   CHECK_OKAY( SCIPhashtableCreate(&(*paramset)->hashtable, SCIP_HASHSIZE_PARAMS,
                  hashGetKeyParam, SCIPhashKeyEqString, SCIPhashKeyValString) );

   (*paramset)->params = NULL;
   (*paramset)->nparams = 0;
   (*paramset)->paramssize = 0;

   return SCIP_OKAY;
}

/** frees parameter set */
void SCIPparamsetFree(
   PARAMSET**       paramset,           /**< pointer to the parameter set */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   int i;

   assert(paramset != NULL);
   assert(*paramset != NULL);
   assert((*paramset)->paramssize == 0 || (*paramset)->params != NULL);
   assert((*paramset)->paramssize >= (*paramset)->nparams);

   for( i = 0; i < (*paramset)->nparams; ++i )
   {
      paramFree(&(*paramset)->params[i], memhdr);
   }

   SCIPhashtableFree(&(*paramset)->hashtable, memhdr);

   freeMemoryArrayNull(&(*paramset)->params);
   freeMemory(paramset);
}

/** adds parameter to the parameter set */
static
RETCODE paramsetAdd(
   PARAMSET*        paramset,           /**< parameter set */
   MEMHDR*          memhdr,             /**< block memory */
   PARAM*           param               /**< parameter to add */
   )
{
   assert(paramset != NULL);
   assert(param != NULL);

   /* insert the parameter name to the hash table */
   CHECK_OKAY( SCIPhashtableInsert(paramset->hashtable, memhdr, (void*)param) );

   /* ensure, that there is enough space in the params array */
   if( paramset->nparams >= paramset->paramssize )
   {
      paramset->paramssize *= 2;
      paramset->paramssize = MAX(paramset->paramssize, paramset->nparams+1);
      ALLOC_OKAY( reallocMemoryArray(&paramset->params, paramset->paramssize) );
   }
   assert(paramset->nparams < paramset->paramssize);
   
   /* insert parameter in the params array */
   paramset->params[paramset->nparams] = param;
   paramset->nparams++;
   
   return SCIP_OKAY;
}

/** creates a Bool parameter, sets it to its default value, and adds it to the parameter set */
RETCODE SCIPparamsetAddBool(
   PARAMSET*        paramset,           /**< parameter set */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   Bool*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   Bool             defaultvalue        /**< default value of the parameter */
   )
{
   PARAM* param;

   assert(paramset != NULL);

   /* create the parameter */
   CHECK_OKAY( paramCreateBool(&param, memhdr, name, valueptr, defaultvalue) );

   /* add parameter to the parameter set */
   CHECK_OKAY( paramsetAdd(paramset, memhdr, param) );
   
   return SCIP_OKAY;
}

/** creates a int parameter, sets it to its default value, and adds it to the parameter set */
RETCODE SCIPparamsetAddInt(
   PARAMSET*        paramset,           /**< parameter set */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   int*             valueptr,           /**< pointer to store the current parameter value, or NULL */
   int              defaultvalue        /**< default value of the parameter */
   )
{
   PARAM* param;

   assert(paramset != NULL);

   /* create the parameter */
   CHECK_OKAY( paramCreateInt(&param, memhdr, name, valueptr, defaultvalue) );

   /* add parameter to the parameter set */
   CHECK_OKAY( paramsetAdd(paramset, memhdr, param) );
   
   return SCIP_OKAY;
}

/** creates a Longint parameter, sets it to its default value, and adds it to the parameter set */
RETCODE SCIPparamsetAddLongint(
   PARAMSET*        paramset,           /**< parameter set */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   Longint*         valueptr,           /**< pointer to store the current parameter value, or NULL */
   Longint          defaultvalue        /**< default value of the parameter */
   )
{
   PARAM* param;

   assert(paramset != NULL);

   /* create the parameter */
   CHECK_OKAY( paramCreateLongint(&param, memhdr, name, valueptr, defaultvalue) );

   /* add parameter to the parameter set */
   CHECK_OKAY( paramsetAdd(paramset, memhdr, param) );
   
   return SCIP_OKAY;
}

/** creates a Real parameter, sets it to its default value, and adds it to the parameter set */
RETCODE SCIPparamsetAddReal(
   PARAMSET*        paramset,           /**< parameter set */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   Real*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   Real             defaultvalue        /**< default value of the parameter */
   )
{
   PARAM* param;

   assert(paramset != NULL);

   /* create the parameter */
   CHECK_OKAY( paramCreateReal(&param, memhdr, name, valueptr, defaultvalue) );

   /* add parameter to the parameter set */
   CHECK_OKAY( paramsetAdd(paramset, memhdr, param) );
   
   return SCIP_OKAY;
}

/** creates a char parameter, sets it to its default value, and adds it to the parameter set */
RETCODE SCIPparamsetAddChar(
   PARAMSET*        paramset,           /**< parameter set */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   char*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   char             defaultvalue        /**< default value of the parameter */
   )
{
   PARAM* param;

   assert(paramset != NULL);

   /* create the parameter */
   CHECK_OKAY( paramCreateChar(&param, memhdr, name, valueptr, defaultvalue) );

   /* add parameter to the parameter set */
   CHECK_OKAY( paramsetAdd(paramset, memhdr, param) );
   
   return SCIP_OKAY;
}

/** creates a string parameter, sets it to its default value, and adds it to the parameter set */
RETCODE SCIPparamsetAddString(
   PARAMSET*        paramset,           /**< parameter set */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   char**           valueptr,           /**< pointer to store the current parameter value, or NULL */
   const char*      defaultvalue        /**< default value of the parameter */
   )
{
   PARAM* param;

   assert(paramset != NULL);

   /* create the parameter */
   CHECK_OKAY( paramCreateString(&param, memhdr, name, valueptr, defaultvalue) );

   /* add parameter to the parameter set */
   CHECK_OKAY( paramsetAdd(paramset, memhdr, param) );
   
   return SCIP_OKAY;
}

/** gets the value of an existing Bool parameter */
RETCODE SCIPparamsetGetBool(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   Bool*            value               /**< pointer to store the parameter */
   )
{
   PARAM* param;

   assert(paramset != NULL);
   assert(value != NULL);

   /* retrieve parameter from hash table */
   param = SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
      return SCIP_UNKNOWNPARAMETER;
   if( param->paramtype != SCIP_PARAMTYPE_BOOL )
      return SCIP_WRONGPARAMETERTYPE;

   /* set the actual parameter's value */
   if( param->data.boolparam.valueptr != NULL )
      *value = *param->data.boolparam.valueptr;
   else
      *value = param->data.boolparam.actvalue;

   return SCIP_OKAY;
}

/** gets the value of an existing int parameter */
RETCODE SCIPparamsetGetInt(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   int*             value               /**< pointer to store the parameter */
   )
{
   PARAM* param;

   assert(paramset != NULL);
   assert(value != NULL);

   /* retrieve parameter from hash table */
   param = SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
      return SCIP_UNKNOWNPARAMETER;
   if( param->paramtype != SCIP_PARAMTYPE_INT )
      return SCIP_WRONGPARAMETERTYPE;

   /* set the actual parameter's value */
   if( param->data.intparam.valueptr != NULL )
      *value = *param->data.intparam.valueptr;
   else
      *value = param->data.intparam.actvalue;

   return SCIP_OKAY;
}

/** gets the value of an existing Longint parameter */
RETCODE SCIPparamsetGetLongint(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   Longint*         value               /**< pointer to store the parameter */
   )
{
   PARAM* param;

   assert(paramset != NULL);
   assert(value != NULL);

   /* retrieve parameter from hash table */
   param = SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
      return SCIP_UNKNOWNPARAMETER;
   if( param->paramtype != SCIP_PARAMTYPE_LONGINT )
      return SCIP_WRONGPARAMETERTYPE;

   /* set the actual parameter's value */
   if( param->data.longintparam.valueptr != NULL )
      *value = *param->data.longintparam.valueptr;
   else
      *value = param->data.longintparam.actvalue;

   return SCIP_OKAY;
}

/** gets the value of an existing Real parameter */
RETCODE SCIPparamsetGetReal(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   Real*            value               /**< pointer to store the parameter */
   )
{
   PARAM* param;

   assert(paramset != NULL);
   assert(value != NULL);

   /* retrieve parameter from hash table */
   param = SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
      return SCIP_UNKNOWNPARAMETER;
   if( param->paramtype != SCIP_PARAMTYPE_REAL )
      return SCIP_WRONGPARAMETERTYPE;

   /* set the actual parameter's value */
   if( param->data.realparam.valueptr != NULL )
      *value = *param->data.realparam.valueptr;
   else
      *value = param->data.realparam.actvalue;

   return SCIP_OKAY;
}

/** gets the value of an existing char parameter */
RETCODE SCIPparamsetGetChar(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   char*            value               /**< pointer to store the parameter */
   )
{
   PARAM* param;

   assert(paramset != NULL);
   assert(value != NULL);

   /* retrieve parameter from hash table */
   param = SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
      return SCIP_UNKNOWNPARAMETER;
   if( param->paramtype != SCIP_PARAMTYPE_CHAR )
      return SCIP_WRONGPARAMETERTYPE;

   /* set the actual parameter's value */
   if( param->data.charparam.valueptr != NULL )
      *value = *param->data.charparam.valueptr;
   else
      *value = param->data.charparam.actvalue;

   return SCIP_OKAY;
}

/** gets the value of an existing string parameter */
RETCODE SCIPparamsetGetString(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   char**           value               /**< pointer to store the parameter */
   )
{
   PARAM* param;

   assert(paramset != NULL);
   assert(value != NULL);

   /* retrieve parameter from hash table */
   param = SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
      return SCIP_UNKNOWNPARAMETER;
   if( param->paramtype != SCIP_PARAMTYPE_STRING )
      return SCIP_WRONGPARAMETERTYPE;

   /* set the actual parameter's value */
   if( param->data.stringparam.valueptr != NULL )
      *value = *param->data.stringparam.valueptr;
   else
      *value = param->data.stringparam.actvalue;

   return SCIP_OKAY;
}

/** changes the value of an existing Bool parameter */
RETCODE SCIPparamsetSetBool(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   Bool             value               /**< new value of the parameter */
   )
{
   PARAM* param;

   assert(paramset != NULL);

   /* retrieve parameter from hash table */
   param = SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
      return SCIP_UNKNOWNPARAMETER;
   if( param->paramtype != SCIP_PARAMTYPE_BOOL )
      return SCIP_WRONGPARAMETERTYPE;

   /* set the actual parameter's value */
   if( param->data.boolparam.valueptr != NULL )
      *param->data.boolparam.valueptr = value;
   else
      param->data.boolparam.actvalue = value;

   return SCIP_OKAY;
}

/** changes the value of an existing int parameter */
RETCODE SCIPparamsetSetInt(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   int              value               /**< new value of the parameter */
   )
{
   PARAM* param;

   assert(paramset != NULL);

   /* retrieve parameter from hash table */
   param = SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
      return SCIP_UNKNOWNPARAMETER;
   if( param->paramtype != SCIP_PARAMTYPE_INT )
      return SCIP_WRONGPARAMETERTYPE;

   /* set the actual parameter's value */
   if( param->data.intparam.valueptr != NULL )
      *param->data.intparam.valueptr = value;
   else
      param->data.intparam.actvalue = value;

   return SCIP_OKAY;
}

/** changes the value of an existing Longint parameter */
RETCODE SCIPparamsetSetLongint(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   Longint          value               /**< new value of the parameter */
   )
{
   PARAM* param;

   assert(paramset != NULL);

   /* retrieve parameter from hash table */
   param = SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
      return SCIP_UNKNOWNPARAMETER;
   if( param->paramtype != SCIP_PARAMTYPE_LONGINT )
      return SCIP_WRONGPARAMETERTYPE;

   /* set the actual parameter's value */
   if( param->data.longintparam.valueptr != NULL )
      *param->data.longintparam.valueptr = value;
   else
      param->data.longintparam.actvalue = value;

   return SCIP_OKAY;
}

/** changes the value of an existing Real parameter */
RETCODE SCIPparamsetSetReal(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   Real             value               /**< new value of the parameter */
   )
{
   PARAM* param;

   assert(paramset != NULL);

   /* retrieve parameter from hash table */
   param = SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
      return SCIP_UNKNOWNPARAMETER;
   if( param->paramtype != SCIP_PARAMTYPE_REAL )
      return SCIP_WRONGPARAMETERTYPE;

   /* set the actual parameter's value */
   if( param->data.realparam.valueptr != NULL )
      *param->data.realparam.valueptr = value;
   else
      param->data.realparam.actvalue = value;

   return SCIP_OKAY;
}

/** changes the value of an existing char parameter */
RETCODE SCIPparamsetSetChar(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   char             value               /**< new value of the parameter */
   )
{
   PARAM* param;

   assert(paramset != NULL);

   /* retrieve parameter from hash table */
   param = SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
      return SCIP_UNKNOWNPARAMETER;
   if( param->paramtype != SCIP_PARAMTYPE_CHAR )
      return SCIP_WRONGPARAMETERTYPE;

   /* set the actual parameter's value */
   if( param->data.charparam.valueptr != NULL )
      *param->data.charparam.valueptr = value;
   else
      param->data.charparam.actvalue = value;

   return SCIP_OKAY;
}

/** changes the value of an existing string parameter */
RETCODE SCIPparamsetSetString(
   PARAMSET*        paramset,           /**< parameter set */
   const char*      name,               /**< name of the parameter */
   const char*      value               /**< new value of the parameter */
   )
{
   PARAM* param;

   assert(paramset != NULL);

   /* retrieve parameter from hash table */
   param = SCIPhashtableRetrieve(paramset->hashtable, (void*)name);
   if( param == NULL )
      return SCIP_UNKNOWNPARAMETER;
   if( param->paramtype != SCIP_PARAMTYPE_STRING )
      return SCIP_WRONGPARAMETERTYPE;

   /* set the actual parameter's value */
   if( param->data.stringparam.valueptr != NULL )
   {
      freeMemoryArray(param->data.stringparam.valueptr);
      duplicateMemoryArray(param->data.stringparam.valueptr, value, strlen(value)+1);
   }
   else
   {
      freeMemoryArray(&param->data.stringparam.actvalue);
      duplicateMemoryArray(&param->data.stringparam.actvalue, value, strlen(value)+1);
   }

   return SCIP_OKAY;
}


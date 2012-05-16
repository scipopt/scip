/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_osil.c
 * @brief  OS instance language (OSiL) format file reader
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <math.h>

#include "scip/reader_osil.h"
#include "scip/scip.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/cons_linear.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_nonlinear.h"
#include "xml/xml.h"


#define READER_NAME             "osilreader"
#define READER_DESC             "file reader for OS instance language (OSiL) format"
#define READER_EXTENSION        "osil"

#ifndef M_PI
#define M_PI           3.141592653589793238462643
#endif

#ifndef M_E
#define M_E            2.7182818284590452354
#endif


/*
 * Data structures
 */

typedef enum {
   LINEAR,
   QUADRATIC,
   NONLINEAR

} CONSTYPE;


/*
 * Local methods
 */

/** create variables with bounds and type according to xml data */
static
SCIP_RETCODE readVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       datanode,           /**< XML root node for instance data */
   SCIP_VAR***           vars,               /**< buffer to store pointer to variable array */
   int*                  nvars,              /**< buffer to store number of variables */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occured */
   )
{
   const XML_NODE* variables;
   const XML_NODE* varnode;
   const char* attrval;
   int varssize;

   assert(scip != NULL);
   assert(datanode != NULL);
   assert(vars != NULL);
   assert(nvars != NULL);
   assert(doingfine != NULL);

   *vars = NULL;
   *nvars = 0;

   variables = xmlFindNodeMaxdepth(datanode, "variables", 0, 1);

   if( variables == NULL )
   {
      /* no variables: strange but ok so far */
      return SCIP_OKAY;
   }

   attrval = xmlGetAttrval(variables, "numberOfVariables");
   if( attrval == NULL )
   {
      SCIPerrorMessage("numberOfVariables attribute not found.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   varssize = (int)strtol(attrval, (char**)&attrval, 10);
   if( *attrval != '\0' || varssize < 0 )
   {
      SCIPerrorMessage("Invalid value for numberOfVariables attribute.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }
   assert(varssize >= 0);

   SCIP_CALL( SCIPallocBufferArray(scip, vars, varssize) );

   for( varnode = xmlFirstChild(variables); varnode != NULL; varnode = xmlNextSibl(varnode) )
   {
      const char* varname;
      SCIP_VARTYPE vartype;
      SCIP_Real varlb;
      SCIP_Real varub;
      SCIP_Real semibound;

      if( varssize == *nvars )
      {
         SCIPerrorMessage("More variables than expected.\n");
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* find variable name */
      varname = xmlGetAttrval(varnode, "name");

      /* find variable lower bound (default is 0.0 !) */
      attrval = xmlGetAttrval(varnode, "lb");
      if( attrval == NULL )
         varlb = 0.0;
      else if( strcmp(attrval, "-INF") == 0 )
         varlb = -SCIPinfinity(scip);
      else if( strcmp(attrval, "INF") == 0 )
         varlb = SCIPinfinity(scip);
      else
      {
         varlb = strtod(attrval, (char**)&attrval);
         if( *attrval != '\0' )
         {
            SCIPerrorMessage("error parsing variable lower bound %s\n", attrval);
            *doingfine = FALSE;
            return SCIP_OKAY;
         }
      }

      /* find variable upper bound (default is infinity) */
      attrval = xmlGetAttrval(varnode, "ub");
      if( attrval == NULL )
         varub = SCIPinfinity(scip);
      else if( strcmp(attrval, "-INF") == 0 )
         varub = -SCIPinfinity(scip);
      else if( strcmp(attrval, "INF") == 0 )
         varub = SCIPinfinity(scip);
      else
      {
         varub = strtod(attrval, (char**)&attrval);
         if( *attrval != '\0' )
         {
            SCIPerrorMessage("error parsing variable upper bound %s\n", attrval);
            *doingfine = FALSE;
            return SCIP_OKAY;
         }
      }

      semibound = SCIP_INVALID;

      /* find variable type (default is continuous) */
      attrval = xmlGetAttrval(varnode, "type");
      if( attrval == NULL )
         vartype = SCIP_VARTYPE_CONTINUOUS;
      else switch( *attrval )
      {
      case 'C':
         vartype = SCIP_VARTYPE_CONTINUOUS;
         break;
      case 'B':
         vartype = SCIP_VARTYPE_BINARY;
         if( varub > 1.0 )
            varub = 1.0;
         break;
      case 'I':
         vartype = SCIP_VARTYPE_INTEGER;
         break;
      case 'D':
         vartype = SCIP_VARTYPE_CONTINUOUS;
         if( varlb > 0.0 )
            semibound = varlb;
         varlb = 0.0;
         break;
      case 'J':
         vartype = SCIP_VARTYPE_INTEGER;
         if( varlb > 0.0 )
            semibound = varlb;
         varlb = 0.0;
         break;
      default:
         SCIPerrorMessage("Unsupported variable type: %s\n", attrval);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      if( vartype != SCIP_VARTYPE_CONTINUOUS )
      {
         varlb = SCIPceil(scip, varlb);
         varub = SCIPfloor(scip, varub);
      }

      /* create SCIP variable */
      SCIP_CALL( SCIPcreateVar(scip, &(*vars)[*nvars], varname, varlb, varub, 0.0, vartype, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
      assert((*vars)[*nvars] != NULL);

      SCIP_CALL( SCIPaddVar(scip, (*vars)[*nvars]) );

      /* if variable is actually semicontinuous or semiintegral, create bounddisjunction constraint (var <= 0.0 || var >= semibound) */
      if( semibound != SCIP_INVALID )
      {
         SCIP_CONS* cons;
         SCIP_VAR* consvars[2];
         SCIP_BOUNDTYPE boundtypes[2];
         SCIP_Real bounds[2];
         char name[SCIP_MAXSTRLEN];

         consvars[0] = (*vars)[*nvars];
         consvars[1] = (*vars)[*nvars];

         boundtypes[0] = SCIP_BOUNDTYPE_UPPER;
         boundtypes[1] = SCIP_BOUNDTYPE_LOWER;

         bounds[0] = 0.0;
         bounds[1] = semibound;

         SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_semibound", SCIPvarGetName((*vars)[*nvars]));

         SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &cons, name, 2, consvars, boundtypes, bounds, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }

      ++*nvars;
   }
   if( *nvars < varssize )
   {
      SCIPerrorMessage("Less variables than expected.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** setup linear coefficients and constant of objective and objective sense */
static
SCIP_RETCODE readObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       datanode,           /**< XML root node for instance data */
   SCIP_VAR**            vars,               /**< variables in order of OSiL indices */
   int                   nvars,              /**< number of variables */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occured */
   )
{
   const XML_NODE* objective;
   const XML_NODE* coefnode;
   const char* attrval;

   assert(scip != NULL);
   assert(datanode != NULL);
   assert(vars != NULL || nvars == 0);
   assert(doingfine != NULL);

   /* check for first objective */
   objective = xmlFindNodeMaxdepth(datanode, "obj", 0, 2);

   /* if no objective, then nothing to do here */
   if( objective == NULL )
      return SCIP_OKAY;

   /* objective sense */
   attrval = xmlGetAttrval(objective, "maxOrMin");
   if( attrval == NULL )
   {
      SCIPerrorMessage("Objective sense missing\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }
   else if( strcmp(attrval, "min") == 0 )
   {
      SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE);
   }
   else if( strcmp(attrval, "max") == 0 )
   {
      SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE);
   }
   else
   {
      SCIPerrorMessage("Cannot parse objective sense: %s\n", attrval);
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   /* objective coefficients */
   for( coefnode = xmlFirstChild(objective); coefnode != NULL; coefnode = xmlNextSibl(coefnode) )
   {
      SCIP_Real val;
      int idx;

      attrval = xmlGetAttrval(coefnode, "idx");
      if( attrval == NULL )
      {
         SCIPerrorMessage("No idx in objective coefficient.\n");
         *doingfine = FALSE;
         return SCIP_OKAY;
      }
      idx = (int)strtol(attrval, (char**)&attrval, 10);
      if( *attrval != '\0' )
      {
         SCIPerrorMessage("Error parsing objective coefficient index: %s\n", xmlGetAttrval(coefnode, "idx"));
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      if( idx < 0 || idx >= nvars )
      {
         SCIPerrorMessage("Invalid objective coefficient index: %d\n", idx);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      if( xmlFirstChild(coefnode) == NULL || xmlGetData(xmlFirstChild(coefnode)) == NULL )
      {
         SCIPerrorMessage("No coefficient stored for objective coefficient %d\n", idx);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      attrval = xmlGetData(xmlFirstChild(coefnode));
      val = strtod(attrval, (char**)&attrval);
      if( *attrval != '\0' )
      {
         SCIPerrorMessage("Error parsing objective coefficient value: %s\n", xmlGetData(xmlFirstChild(coefnode)));
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      SCIP_CALL( SCIPchgVarObj(scip, vars[idx], val) );
   }

   /* objective constant: model as fixed variable, if nonzero */
   attrval = xmlGetAttrval(objective, "constant");
   if( attrval != NULL )
   {
      SCIP_Real objconst;

      objconst = strtod(attrval, (char**)&attrval);
      if( *attrval != '\0' )
      {
         SCIPerrorMessage("Error parsing objective constant: %s\n", xmlGetAttrval(objective, "constant"));
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      if( objconst != 0.0 )
      {
         SCIP_VAR* objconstvar;

         SCIP_CALL( SCIPcreateVar(scip, &objconstvar, "objconstvar", objconst, objconst, 1.0, SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, objconstvar) );
         SCIP_CALL( SCIPreleaseVar(scip, &objconstvar) );
      }
   }

   if( xmlNextSibl(objective) != NULL )
   {
      SCIPerrorMessage("Multiple objectives not supported by SCIP.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** setup constraint sides as linear constraints
 *
 * constraints are not added to the problem yet
 */
static
SCIP_RETCODE readConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       datanode,           /**< XML root node for instance data */
   SCIP_CONS***          conss,              /**< buffer to store array of (linear) constraints */
   CONSTYPE**            constypes,          /**< buffer to store type of constraints (will be all LINEAR) */
   int*                  nconss,             /**< buffer to store number of constraints */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occured */
   )
{
   const XML_NODE* constraints;
   const XML_NODE* consnode;
   const char* attrval;
   int consssize;
   char name[20];

   assert(scip != NULL);
   assert(datanode != NULL);
   assert(conss != NULL);
   assert(constypes != NULL);
   assert(nconss != NULL);
   assert(doingfine != NULL);

   *conss = NULL;
   *constypes = NULL;
   *nconss = 0;

   constraints = xmlFindNodeMaxdepth(datanode, "constraints", 0, 1);

   /* if no constraints, then nothing to do here */
   if( constraints == NULL )
      return SCIP_OKAY;

   attrval = xmlGetAttrval(constraints, "numberOfConstraints");
   if( attrval == NULL )
   {
      SCIPerrorMessage("numberOfConstraints attribute not found.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   consssize = (int)strtol(attrval, (char**)&attrval, 10);
   if( *attrval != '\0' || consssize < 0 )
   {
      SCIPerrorMessage("Invalid value for numberOfConstraints attribute.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }
   assert(consssize >= 0);

   SCIP_CALL( SCIPallocBufferArray(scip, conss, consssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, constypes, consssize) );

   for( consnode = xmlFirstChild(constraints); consnode != NULL; consnode = xmlNextSibl(consnode) )
   {
      const char* consname;
      SCIP_Real conslhs;
      SCIP_Real consrhs;

      if( consssize == *nconss )
      {
         SCIPerrorMessage("More constraints than expected.\n");
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* find constraint name */
      consname = xmlGetAttrval(consnode, "name");
      if( consname == NULL )
      {
         SCIPsnprintf(name, 20, "cons%d", *nconss);
         consname = name;
      }

      /* find constraint lower bound (=lhs) (default is -infinity) */
      attrval = xmlGetAttrval(consnode, "lb");
      if( attrval == NULL )
         conslhs = -SCIPinfinity(scip);
      else if( strcmp(attrval, "-INF") == 0 )
         conslhs = -SCIPinfinity(scip);
      else if( strcmp(attrval, "INF") == 0 )
         conslhs = SCIPinfinity(scip);
      else
      {
         conslhs = strtod(attrval, (char**)&attrval);
         if( *attrval != '\0' )
         {
            SCIPerrorMessage("error parsing constraint lower bound %s\n", attrval);
            *doingfine = FALSE;
            return SCIP_OKAY;
         }
      }

      /* find constraint upper bound (=rhs) (default is +infinity) */
      attrval = xmlGetAttrval(consnode, "ub");
      if( attrval == NULL )
         consrhs = SCIPinfinity(scip);
      else if( strcmp(attrval, "-INF") == 0 )
         consrhs = -SCIPinfinity(scip);
      else if( strcmp(attrval, "INF") == 0 )
         consrhs = SCIPinfinity(scip);
      else
      {
         consrhs = strtod(attrval, (char**)&attrval);
         if( *attrval != '\0' )
         {
            SCIPerrorMessage("error parsing constraint upper bound %s\n", attrval);
            *doingfine = FALSE;
            return SCIP_OKAY;
         }
      }

      /* find constraint constant (default is 0.0) and substract from lhs/rhs */
      attrval = xmlGetAttrval(consnode, "constant");
      if( attrval != NULL )
      {
         SCIP_Real consconstant;

         consconstant = strtod(attrval, (char**)&attrval);
         if( *attrval != '\0' )
         {
            SCIPerrorMessage("error parsing constraint constant %s\n", attrval);
            *doingfine = FALSE;
            return SCIP_OKAY;
         }
         if( conslhs > -SCIPinfinity(scip) )
            conslhs -= consconstant;
         if( consrhs <  SCIPinfinity(scip) )
            consrhs -= consconstant;
      }

      /* create SCIP linear constraint */
      SCIP_CALL( SCIPcreateConsLinear(scip, &(*conss)[*nconss], consname, 0, NULL, NULL, conslhs, consrhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      assert((*conss)[*nconss] != NULL);

      (*constypes)[*nconss] = LINEAR;

      ++*nconss;
   }
   if( *nconss < consssize )
   {
      SCIPerrorMessage("Less constraints than expected.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** parse linear coefficients of constraints */
static
SCIP_RETCODE readLinearCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       datanode,           /**< XML root node for instance data */
   SCIP_VAR**            vars,               /**< variables in order of OSiL indices */
   int                   nvars,              /**< number of variables */
   SCIP_CONS**           conss,              /**< constraints in order of OSiL indices */
   CONSTYPE*             constypes,          /**< type of constraints (assumed to be LINEAR) */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occured */
   )
{
   const XML_NODE* lincoef;
   const XML_NODE* startnode;
   const XML_NODE* idxnode;
   const XML_NODE* valnode;
   const XML_NODE* elnode;
   const char* attrval;
   SCIP_Bool rowmajor;
   int* start;
   int* idx;
   SCIP_Real* val;
   int nnz;
   int count;

   assert(scip != NULL);
   assert(datanode != NULL);
   assert(vars != NULL || nvars == 0);
   assert(conss != NULL || nconss == 0);
   assert(constypes != NULL || nconss == 0);
   assert(doingfine != NULL);

   lincoef = xmlFindNodeMaxdepth(datanode, "linearConstraintCoefficients", 0, 1);

   if( lincoef == NULL )
      return SCIP_OKAY;

   attrval = xmlGetAttrval(lincoef, "numberOfValues");
   if( attrval == NULL )
   {
      SCIPerrorMessage("numberOfValues attribute for linearConstraintCoefficients not found.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   nnz = (int)strtol(attrval, (char**)&attrval, 10);
   if( *attrval != '\0' || nnz < 0 )
   {
      SCIPerrorMessage("Invalid value for numberOfValues attribute.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }
   assert(nnz >= 0);

   startnode = xmlFindNodeMaxdepth(lincoef, "start", 0, 1);
   if( startnode == NULL )
   {
      SCIPerrorMessage("start node not found under linearConstraintCoefficients.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   idxnode = xmlFindNodeMaxdepth(lincoef, "rowIdx", 0, 1);
   if( idxnode != NULL )
   {
      if( xmlFindNodeMaxdepth(lincoef, "colIdx", 0, 1) != NULL )
      {
         SCIPerrorMessage("Both rowIdx and colIdx found under linearConstraintCoefficients.\n");
         *doingfine = FALSE;
         return SCIP_OKAY;
      }
      rowmajor = FALSE;
   }
   else
   {
      idxnode = xmlFindNodeMaxdepth(lincoef, "colIdx", 0, 1);
      if( idxnode == NULL )
      {
         SCIPerrorMessage("rowIdx and colIdx not found under linearConstraintCoefficients.\n");
         *doingfine = FALSE;
         return SCIP_OKAY;
      }
      rowmajor = TRUE;
   }

   valnode = xmlFindNodeMaxdepth(lincoef, "value", 0, 1);
   if( valnode == NULL )
   {
      SCIPerrorMessage("value node not found under linearConstraintCoefficients.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   start = NULL;
   idx = NULL;
   val = NULL;

   /* read row or column start indices */
   SCIP_CALL( SCIPallocBufferArray(scip, &start, (rowmajor ? nconss : nvars) + 1) );

   count = 0;
   for( elnode = xmlFirstChild(startnode); elnode != NULL; elnode = xmlNextSibl(elnode), ++count )
   {
      if( strcmp(xmlGetName(elnode), "el") != 0 )
      {
         SCIPerrorMessage("expected <el>-node, but got %s\n", xmlGetName(elnode));
         *doingfine = FALSE;
         goto CLEANUP;
      }
      if( count >= (rowmajor ? nconss : nvars) + 1 )
      {
         SCIPerrorMessage("too many elements under start node, expected %d, got at least %d.\n", (rowmajor ? nconss : nvars) + 1, count + 1);
         *doingfine = FALSE;
         goto CLEANUP;
      }
      if( xmlFirstChild(elnode) == NULL || xmlGetData(xmlFirstChild(elnode)) == NULL )
      {
         SCIPerrorMessage("<el>-node without data.\n");
         *doingfine = FALSE;
         goto CLEANUP;
      }

      start[count] = (int)strtol(xmlGetData(xmlFirstChild(elnode)), (char**)&attrval, 10);

      if( *attrval != '\0' || start[count] < 0 || (start[count] > nnz) )
      {
         SCIPerrorMessage("Invalid value '%s' in <el> node.\n", xmlGetData(elnode));
         *doingfine = FALSE;
         goto CLEANUP;
      }
   }
   if( count != (rowmajor ? nconss : nvars) + 1 )
   {
      SCIPerrorMessage("Got only %d start entries, expected %d many.\n", count, (rowmajor ? nconss : nvars) + 1);
      *doingfine = FALSE;
      goto CLEANUP;
   }

   /* read row or column indices */
   SCIP_CALL( SCIPallocBufferArray(scip, &idx, nnz) );

   count = 0;
   for( elnode = xmlFirstChild(idxnode); elnode != NULL; elnode = xmlNextSibl(elnode), ++count )
   {
      if( strcmp(xmlGetName(elnode), "el") != 0 )
      {
         SCIPerrorMessage("expected <el>-node, but got %s\n", xmlGetName(elnode));
         *doingfine = FALSE;
         goto CLEANUP;
      }
      if( count >= nnz )
      {
         SCIPerrorMessage("too many elements under rowIdx or colIdx node, expected %d, got at least %d.\n", nnz, count + 1);
         *doingfine = FALSE;
         goto CLEANUP;
      }
      if( xmlFirstChild(elnode) == NULL || xmlGetData(xmlFirstChild(elnode)) == NULL )
      {
         SCIPerrorMessage("<el>-node without data.\n");
         *doingfine = FALSE;
         goto CLEANUP;
      }

      idx[count] = (int)strtol(xmlGetData(xmlFirstChild(elnode)), (char**)&attrval, 10);

      if( *attrval != '\0' || idx[count] < 0 || (idx[count] >= (rowmajor ? nvars : nconss)) )
      {
         SCIPerrorMessage("Invalid value '%s' in <el> node.\n", xmlGetData(elnode));
         *doingfine = FALSE;
         goto CLEANUP;
      }
   }
   if( count != nnz )
   {
      SCIPerrorMessage("Got only %d rowIdx or colIdx entries, expected %d many.\n", count, nnz);
      *doingfine = FALSE;
      goto CLEANUP;
   }

   /* read coefficient values */
   SCIP_CALL( SCIPallocBufferArray(scip, &val, nnz) );

   count = 0;
   for( elnode = xmlFirstChild(valnode); elnode != NULL; elnode = xmlNextSibl(elnode), ++count )
   {
      if( strcmp(xmlGetName(elnode), "el") != 0 )
      {
         SCIPerrorMessage("expected <el>-node, but got %s\n", xmlGetName(elnode));
         *doingfine = FALSE;
         goto CLEANUP;
      }
      if( count >= nnz )
      {
         SCIPerrorMessage("too many elements under value node, expected %d, got at least %d.\n", nnz, count + 1);
         *doingfine = FALSE;
         goto CLEANUP;
      }
      if( xmlFirstChild(elnode) == NULL || xmlGetData(xmlFirstChild(elnode)) == NULL )
      {
         SCIPerrorMessage("<el>-node without data.\n");
         *doingfine = FALSE;
         goto CLEANUP;
      }

      val[count] = strtod(xmlGetData(xmlFirstChild(elnode)), (char**)&attrval);

      if( *attrval != '\0' || (val[count] != val[count]) )
      {
         SCIPerrorMessage("Invalid value '%s' in <el> node.\n", xmlGetData(elnode));
         *doingfine = FALSE;
         goto CLEANUP;
      }
   }
   if( count != nnz )
   {
      SCIPerrorMessage("Got only %d value entries, expected %d many.\n", count, nnz);
      *doingfine = FALSE;
      goto CLEANUP;
   }

   /* add coefficients to linear cons */
   if( rowmajor )
   {
      int row;
      int pos;
      for( row = 0; row < nconss; ++row )
      {
         /* these asserts were checked above */
         assert(start[row] >= 0);
         assert(start[row+1] >= 0);
         assert(start[row] < nnz);
         assert(start[row+1] <= nnz);
         for( pos = start[row]; pos < start[row+1]; ++pos )
         {
            /* these asserts were checked above */
            assert(pos >= 0);
            assert(pos < nnz);
            assert(idx[pos] >= 0);
            assert(idx[pos] < nvars);

            assert(constypes[row] == LINEAR);

            SCIP_CALL( SCIPaddCoefLinear(scip, conss[row], vars[idx[pos]], val[pos]) );
         }
      }
   }
   else
   {
      int col;
      int pos;
      for( col = 0; col < nvars; ++col )
      {
         /* these asserts were checked above */
         assert(start[col] >= 0);
         assert(start[col+1] >= 0);
         assert(start[col] <= nnz);
         assert(start[col+1] <= nnz);
         for( pos = start[col]; pos < start[col+1]; ++pos )
         {
            /* these asserts were checked above */
            assert(pos >= 0);
            assert(pos < nnz);
            assert(idx[pos] >= 0);
            assert(idx[pos] < nconss);

            assert(constypes[idx[pos]] == LINEAR);

            SCIP_CALL( SCIPaddCoefLinear(scip, conss[idx[pos]], vars[col], val[pos]) );
         }
      }
   }

 CLEANUP:
   SCIPfreeBufferArrayNull(scip, &start);
   SCIPfreeBufferArrayNull(scip, &idx);
   SCIPfreeBufferArrayNull(scip, &val);

   return SCIP_OKAY;
}

/** read quadratic coefficients of constraints and objective */
static
SCIP_RETCODE readQuadraticCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       datanode,           /**< XML root node for instance data */
   SCIP_VAR**            vars,               /**< variables in order of OSiL indices */
   int                   nvars,              /**< number of variables */
   SCIP_CONS**           conss,              /**< constraints in order of OSiL indices */
   CONSTYPE*             constypes,          /**< type of constraints (assumed to be LINEAR) */
   int                   nconss,             /**< number of constraints */
   SCIP_CONS**           objcons,            /**< buffer to store constraint for nonlinear part of objective function, or to add to if already existing */
   CONSTYPE*             objconstype,        /**< buffer to store type of objective constraint, if created (should be QUADRATIC) */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occured */
)
{
   const XML_NODE* quadcoef;
   const XML_NODE* qterm;
   const char* attrval;
   SCIP_CONS* cons;
   int nqterms;
   int count;
   int considx;
   int varidx1;
   int varidx2;
   SCIP_Real coef;

   assert(scip != NULL);
   assert(datanode != NULL);
   assert(vars != NULL || nvars == 0);
   assert(conss != NULL || nconss == 0);
   assert(constypes != NULL || nconss == 0);
   assert(objcons != NULL);
   assert(doingfine != NULL);

   quadcoef = xmlFindNodeMaxdepth(datanode, "quadraticCoefficients", 0, 1);

   if( quadcoef == NULL )
      return SCIP_OKAY;

   attrval = xmlGetAttrval(quadcoef, "numberOfQuadraticTerms");
   if( attrval == NULL )
   {
      SCIPerrorMessage("numberOfQuadraticTerms attribute for quadraticCoefficients not found.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   nqterms = (int)strtol(attrval, (char**)&attrval, 10);
   if( *attrval != '\0' || nqterms < 0 )
   {
      SCIPerrorMessage("Invalid value for numberOfQuadraticTerms attribute.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }
   assert(nqterms >= 0);

   count = 0;
   for( qterm = xmlFirstChild(quadcoef); qterm != NULL; qterm = xmlNextSibl(qterm), ++count )
   {
      if( strcmp(xmlGetName(qterm), "qTerm") != 0 )
      {
         SCIPerrorMessage("expected <qTerm>-node, but got %s\n", xmlGetName(qterm));
         *doingfine = FALSE;
         return SCIP_OKAY;
      }
      if( count >= nqterms )
      {
         SCIPerrorMessage("too many quadratic terms, expected %d, got at least %d.\n", nqterms, count + 1);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* get constraint index, or -1 for objective */
      attrval = xmlGetAttrval(qterm, "idx");
      if( attrval == NULL )
      {
         SCIPerrorMessage("Missing 'idx' attribute in qTerm.\n");
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      considx = (int)strtol(attrval, (char**)&attrval, 10);
      if( *attrval != '\0' || considx < -1 || considx >= nconss )
      {
         SCIPerrorMessage("Invalid value '%s' in idx attribute of qTerm node.\n", xmlGetAttrval(qterm, "idx"));
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* get index of first variable */
      attrval = xmlGetAttrval(qterm, "idxOne");
      if( attrval == NULL )
      {
         SCIPerrorMessage("Missing 'idxOne' attribute in qTerm.\n");
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      varidx1 = (int)strtol(attrval, (char**)&attrval, 10);
      if( *attrval != '\0' || varidx1 < 0 || varidx1 >= nvars )
      {
         SCIPerrorMessage("Invalid value '%s' in idxOne attribute of qTerm node.\n", xmlGetAttrval(qterm, "idxOne"));
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* get index of second variable */
      attrval = xmlGetAttrval(qterm, "idxTwo");
      if( attrval == NULL )
      {
         SCIPerrorMessage("Missing 'idxTwo' attribute in qTerm.\n");
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      varidx2 = (int)strtol(attrval, (char**)&attrval, 10);
      if( *attrval != '\0' || varidx2 < 0 || varidx2 >= nvars )
      {
         SCIPerrorMessage("Invalid value '%s' in idxTwo attribute of qTerm node.\n", xmlGetAttrval(qterm, "idxTwo"));
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      /* get (optional) coefficient of quadratic term */
      attrval = xmlGetAttrval(qterm, "coef");
      if( attrval != NULL )
      {
         coef = strtod(attrval, (char**)&attrval);
         if( *attrval != '\0' || (coef != coef) )
         {
            SCIPerrorMessage("Invalid value '%s' in coef attribute of qTerm node.\n", xmlGetAttrval(qterm, "coef"));
            *doingfine = FALSE;
            return SCIP_OKAY;
         }
      }
      else
      {
         /* default is 1.0 according to specification */
         coef = 1.0;
      }

      /* skip zero coefficients */
      if( coef == 0.0 )
         continue;

      if( considx == -1 )
      {
         if( *objcons == NULL )
         {
            /* create constraint to hold quadratic part of objective */
            SCIP_VAR* objvar;
            SCIP_Real minusone;

            SCIP_CALL( SCIPcreateVar(scip, &objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, objvar) );

            minusone = -1.0;
            SCIP_CALL( SCIPcreateConsQuadratic(scip, objcons, "objcons", 1, &objvar, &minusone, 0, NULL, NULL, NULL,
               SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE ? -SCIPinfinity(scip) : 0.0,
               SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE ?  SCIPinfinity(scip) : 0.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
            *objconstype = QUADRATIC;

            SCIP_CALL( SCIPreleaseVar(scip, &objvar) );
         }
         cons = *objcons;
         assert(*objconstype == QUADRATIC);
      }
      else if( constypes[considx] == LINEAR )
      {
         /* replace linear constraint by quadratic constraint */
         cons = conss[considx];

         SCIP_CALL( SCIPcreateConsQuadratic(scip, &cons, SCIPconsGetName(cons),
            SCIPgetNVarsLinear(scip, cons), SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons),
            0, NULL, NULL, NULL,
            SCIPgetLhsLinear(scip, cons), SCIPgetRhsLinear(scip, cons),
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );

         SCIP_CALL( SCIPreleaseCons(scip, &conss[considx]) );

         conss[considx] = cons;
         constypes[considx] = QUADRATIC;
      }
      else
      {
         cons = conss[considx];
         assert(constypes[considx] == QUADRATIC);
      }

      SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, vars[varidx1], vars[varidx2], coef) );
   }

   if( count != nqterms )
   {
      SCIPerrorMessage("Got only %d quadratic terms, expected %d many.\n", count, nqterms);
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** transforms OSnL expression tree into SCIP expression */
static
SCIP_RETCODE readExpression(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< buffer to store pointer to created expression */
   const XML_NODE*       node,               /**< root node of expression to be read */
   int*                  exprvaridx,         /**< array with index of problem variables in expression graph */
   int*                  nexprvars,          /**< number of variables in currently processed expression so far */
   int                   nvars,              /**< total number of variables in problem (and length of exprvaridx array) */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occurred */
   )
{
   const char* exprname;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(node != NULL);
   assert(exprvaridx != NULL || nvars == 0);
   assert(nexprvars != NULL);
   assert(doingfine != NULL);

   exprname = xmlGetName(node);
   assert(exprname != NULL);

   *expr = NULL;

   /* zero argument operands */
   if( strcmp(exprname, "variable") == 0 )
   {
      const char* attrval;
      SCIP_Real coef;
      int idx;

      attrval = xmlGetAttrval(node, "idx");
      if( attrval == NULL )
      {
         SCIPerrorMessage("idx attribute required for variable node\n");
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      idx = (int)strtol(attrval, (char**)&attrval, 10);
      if( *attrval != '\0' || idx < 0 || idx >= nvars )
      {
         SCIPerrorMessage("Invalid value '%s' in idx attribute of variable node.\n", xmlGetAttrval(node, "idx"));
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      attrval = xmlGetAttrval(node, "coef");
      if( attrval != NULL )
      {
         coef = strtod(attrval, (char**)&attrval);
         if( *attrval != '\0' || (coef != coef) )
         {
            SCIPerrorMessage("Invalid value '%s' in coef attribute of number node.\n", xmlGetAttrval(node, "coef"));
            *doingfine = FALSE;
            return SCIP_OKAY;
         }
      }
      else
      {
         coef = 1.0;
      }

      /* assign new index to variable, if we see it the first time */
      if( exprvaridx[idx] == -1 )
      {
         exprvaridx[idx] = *nexprvars;
         ++*nexprvars;
      }

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_VARIDX, exprvaridx[idx]) );
      if( coef != 1.0 )
      {
         SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), expr, 1, expr, &coef, 0.0) );
      }

      return SCIP_OKAY;
   }

   if( strcmp(exprname, "number") == 0 )
   {
      const char* attrval;
      SCIP_Real val;

      attrval = xmlGetAttrval(node, "type");
      if( attrval != NULL && (strcmp(attrval, "real") != 0) )
      {
         SCIPerrorMessage("only number's of type real supported\n");
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      attrval = xmlGetAttrval(node, "value");
      if( attrval != NULL )
      {
         val = strtod(attrval, (char**)&attrval);
         if( *attrval != '\0' || (val != val) )
         {
            SCIPerrorMessage("Invalid value '%s' in value attribute of number node.\n", xmlGetAttrval(node, "value"));
            *doingfine = FALSE;
            return SCIP_OKAY;
         }
      }
      else
      {
         /* according to OSnL.xsd, the value attribute is optional
          * I guess the default is the empty string, which should correspond to 0.0
          */
         val = 0.0;
      }

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_CONST, val) );

      return SCIP_OKAY;
   }

   if( strcmp(exprname, "PI") == 0 )
   {
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_CONST, M_PI) );

      return SCIP_OKAY;
   }

   if( strcmp(exprname, "E") == 0 )
   {
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_CONST, M_E) );

      return SCIP_OKAY;
   }


   /* single argument operands */
   if( strcmp(exprname, "negate") == 0 ||
      strcmp(exprname, "abs") == 0 ||
      strcmp(exprname, "squareRoot") == 0 ||
      strcmp(exprname, "square") == 0 ||
      strcmp(exprname, "exp") == 0 ||
      strcmp(exprname, "ln") == 0 ||
      strcmp(exprname, "log10") == 0
      )
   {
      SCIP_EXPR* arg;

      if( xmlFirstChild(node) == NULL || xmlNextSibl(xmlFirstChild(node)) != NULL )
      {
         SCIPerrorMessage("expected exactly 1 children in <%s> node\n", exprname);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      SCIP_CALL( readExpression(scip, &arg, xmlFirstChild(node), exprvaridx, nexprvars, nvars, doingfine) );
      if( !*doingfine )
         return SCIP_OKAY;

      if( strcmp(exprname, "negate") == 0 )
      {
         SCIP_Real minusone;

         minusone = -1.0;
         SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), expr, 1, &arg, &minusone, 0.0) );
      }
      else if( strcmp(exprname, "abs") == 0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_ABS, arg) );
      }
      else if( strcmp(exprname, "squareRoot") == 0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_SQRT, arg) );
      }
      else if( strcmp(exprname, "square") == 0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_SQUARE, arg) );
      }
      else if( strcmp(exprname, "exp") == 0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_EXP, arg) );
      }
      else if( strcmp(exprname, "ln") == 0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_LOG, arg) );
      }
      else /* if( strcmp(exprname, "log10") == 0 ) */
      {
         /* log10(expr) = ln(expr)*1/ln(10) */
         SCIP_EXPR* tmp;

         assert(strcmp(exprname, "log10") == 0);

         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp, SCIP_EXPR_CONST, 1.0/log(10.0)) );
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &arg, SCIP_EXPR_LOG, arg) );
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_MUL, arg, tmp) );
      }

      return SCIP_OKAY;
   }

   /* two argument operands */
   if( strcmp(exprname, "plus") == 0 ||
      strcmp(exprname, "minus") == 0 ||
      strcmp(exprname, "times") == 0 ||
      strcmp(exprname, "divide") == 0 ||
      strcmp(exprname, "power") == 0 ||
      strcmp(exprname, "log") == 0 ||
      strcmp(exprname, "min") == 0 ||
      strcmp(exprname, "max") == 0
     )
   {
      SCIP_EXPR* arg1;
      SCIP_EXPR* arg2;

      if( xmlFirstChild(node) == NULL ||
         xmlNextSibl(xmlFirstChild(node)) == NULL ||
         xmlNextSibl(xmlNextSibl(xmlFirstChild(node))) != NULL )
      {
         SCIPerrorMessage("expected exactly 2 children in <%s> node\n", exprname);
         *doingfine = FALSE;
         return SCIP_OKAY;
      }

      SCIP_CALL( readExpression(scip, &arg1, xmlFirstChild(node), exprvaridx, nexprvars, nvars, doingfine) );
      if( !*doingfine )
         return SCIP_OKAY;

      SCIP_CALL( readExpression(scip, &arg2, xmlNextSibl(xmlFirstChild(node)), exprvaridx, nexprvars, nvars, doingfine) );
      if( !*doingfine )
      {
         SCIPexprFreeDeep(SCIPblkmem(scip), &arg1);
         return SCIP_OKAY;
      }

      if( strcmp(exprname, "plus") == 0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_PLUS, arg1, arg2) );
      }
      else if( strcmp(exprname, "minus") == 0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_MINUS, arg1, arg2) );
      }
      else if( strcmp(exprname, "times") == 0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_MUL, arg1, arg2) );
      }
      else if( strcmp(exprname, "divide") == 0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_DIV, arg1, arg2) );
      }
      else if( strcmp(exprname, "power") == 0 )
      {
         if( SCIPexprGetOperator(arg2) == SCIP_EXPR_CONST )
         {
            /* expr^number is intpower or realpower */
            if( SCIPisIntegral(scip, SCIPexprGetOpReal(arg2)) )
            {
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_INTPOWER, arg1, (int)SCIPround(scip, SCIPexprGetOpReal(arg2))) );
            }
            else
            {
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_REALPOWER, arg1, SCIPexprGetOpReal(arg2)) );
            }
         }
         else if( SCIPexprGetOperator(arg1) == SCIP_EXPR_CONST )
         {
            /* number^arg2 is exp(arg2 * ln(number)) */
            if( SCIPexprGetOpReal(arg1) < 0.0 )
            {
               SCIPerrorMessage("negative base in power with nonconstant exponent not allowed\n");
               SCIPexprFreeDeep(SCIPblkmem(scip), &arg1);
               SCIPexprFreeDeep(SCIPblkmem(scip), &arg2);
               *doingfine = FALSE;
               return SCIP_OKAY;
            }
            else
            {
               SCIP_EXPR* tmp;

               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp, SCIP_EXPR_CONST, log(SCIPexprGetOpReal(arg1))) );
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp, SCIP_EXPR_MUL, tmp, arg2) );
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_EXP, tmp) );
            }
         }
         else
         {
            /* arg1^arg2 is exp(arg2 * ln(arg1)) */
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &arg1, SCIP_EXPR_LOG, arg1) );
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &arg2, SCIP_EXPR_MUL, arg1, arg2) );
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_EXP, arg2) );
         }
      }
      else if( strcmp(exprname, "log") == 0 )
      {
         /* logarithm of arg2 w.r.t. base arg1 = ln(arg2) / ln(arg1) */
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &arg1, SCIP_EXPR_LOG, arg1) );
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &arg2, SCIP_EXPR_LOG, arg2) );
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_DIV, arg2, arg1) );
      }
      else if( strcmp(exprname, "min") == 0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_MIN, arg1, arg2) );
      }
      else /* if( strcmp(exprname, "max") == 0 ) */
      {
         assert(strcmp(exprname, "max") == 0);

         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_MAX, arg1, arg2) );
      }

      return SCIP_OKAY;
   }

   /* arbitrary argument operands */
   if( strcmp(exprname, "sum") == 0 || strcmp(exprname, "product") == 0 )
   {
      const XML_NODE* argnode;
      SCIP_EXPR** args;
      int nargs;
      int argssize;

      if( xmlFirstChild(node) == NULL )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, SCIP_EXPR_CONST, (strcmp(exprname, "sum") == 0) ? 0.0 : 1.0) );

         return SCIP_OKAY;
      }

      argssize = 5;
      SCIP_CALL( SCIPallocBufferArray(scip, &args, argssize) );

      nargs = 0;
      for( argnode = xmlFirstChild(node); argnode != NULL; argnode = xmlNextSibl(argnode), ++nargs )
      {
         if( nargs >= argssize )
         {
            argssize = SCIPcalcMemGrowSize(scip, nargs + 1);
            SCIP_CALL( SCIPreallocBufferArray(scip, &args, argssize) );
         }
         assert(nargs < argssize);

         SCIP_CALL( readExpression(scip, &args[nargs], argnode, exprvaridx, nexprvars, nvars, doingfine) );
         if( !*doingfine )
         {
            assert(args[nargs] == NULL);
            break;
         }
      }

      if( *doingfine )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), expr, (strcmp(exprname, "sum") == 0) ? SCIP_EXPR_SUM : SCIP_EXPR_PRODUCT, nargs, args) );
      }
      else
      {
         for( ; nargs > 0; --nargs )
            SCIPexprFreeDeep(SCIPblkmem(scip), &args[nargs-1]);
      }

      SCIPfreeBufferArray(scip, &args);

      return SCIP_OKAY;
   }

   if( strcmp(exprname, "quadratic") == 0 )
   {
      const char* attrval;
      const XML_NODE* qterm;
      SCIP_QUADELEM* quadelems;
      int nquadelems;
      int quadelemssize;
      int* quadvarsidxs;
      int nquadvars;
      int i;

      quadelemssize = 5;
      SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, quadelemssize) );
      nquadelems = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &quadvarsidxs, nvars) );
      for( i = 0; i < nvars; ++i )
         quadvarsidxs[i] = -1;
      nquadvars = 0;

      for( qterm = xmlFirstChild(node); qterm != NULL; qterm = xmlNextSibl(qterm), ++nquadelems )
      {
         if( strcmp(xmlGetName(qterm), "qpTerm") != 0 )
         {
            SCIPerrorMessage("unexpected '%s' node in <quadratic> expression\n", xmlGetName(qterm));
            *doingfine = FALSE;
            break;
         }

         if( nquadelems >= quadelemssize )
         {
            quadelemssize = SCIPcalcMemGrowSize(scip, nquadelems + 1);
            SCIP_CALL( SCIPreallocBufferArray(scip, &quadelems, quadelemssize) );
         }
         assert(quadelemssize > nquadelems);

         attrval = xmlGetAttrval(qterm, "idxOne");
         if( attrval == NULL )
         {
            SCIPerrorMessage("missing idxOne attribute in qpTerm\n");
            *doingfine = FALSE;
            break;
         }

         quadelems[nquadelems].idx1 = (int)strtol(attrval, (char**)&attrval, 10);
         if( *attrval != '\0' || quadelems[nquadelems].idx1 < 0 || quadelems[nquadelems].idx1 >= nvars )
         {
            SCIPerrorMessage("invalid value '%s' for idxOne attribute of qpTerm\n", xmlGetAttrval(qterm, "idxOne"));
            *doingfine = FALSE;
            break;
         }

         attrval = xmlGetAttrval(qterm, "idxTwo");
         if( attrval == NULL )
         {
            SCIPerrorMessage("missing idxOne attribute in qpTerm\n");
            *doingfine = FALSE;
            break;
         }

         quadelems[nquadelems].idx2 = (int)strtol(attrval, (char**)&attrval, 10);
         if( *attrval != '\0' || quadelems[nquadelems].idx2 < 0 || quadelems[nquadelems].idx2 >= nvars )
         {
            SCIPerrorMessage("invalid value '%s' for idxTwo attribute of qpTerm\n", xmlGetAttrval(qterm, "idxTwo"));
            *doingfine = FALSE;
            break;
         }

         attrval = xmlGetAttrval(qterm, "coef");
         if( attrval != NULL )
         {
            quadelems[nquadelems].coef = strtod(attrval, (char**)&attrval);
            if( *attrval != '\0' || (quadelems[nquadelems].coef != quadelems[nquadelems].coef) )
            {
               SCIPerrorMessage("invalid value '%s' for coef attribute of qpTerm\n", xmlGetAttrval(qterm, "coef"));
               *doingfine = FALSE;
               break;
            }
         }
         else
         {
            quadelems[nquadelems].coef = 1.0;
         }

         /* get index for first variable in quadratic element */
         if( quadvarsidxs[quadelems[nquadelems].idx1] < 0 )
         {
            quadvarsidxs[quadelems[nquadelems].idx1] = nquadvars;
            quadelems[nquadelems].idx1 = nquadvars;

            ++nquadvars;
         }
         else
         {
            quadelems[nquadelems].idx1 = quadvarsidxs[quadelems[nquadelems].idx1];
         }

         /* get index for second variable in quadratic element */
         if( quadvarsidxs[quadelems[nquadelems].idx2] < 0 )
         {
            quadvarsidxs[quadelems[nquadelems].idx2] = nquadvars;
            quadelems[nquadelems].idx2 = nquadvars;

            ++nquadvars;
         }
         else
         {
            quadelems[nquadelems].idx2 = quadvarsidxs[quadelems[nquadelems].idx2];
         }

         /* swap indices if in wrong order */
         if( quadelems[nquadelems].idx1 > quadelems[nquadelems].idx2 )
         {
            int tmp;

            tmp = quadelems[nquadelems].idx1;
            quadelems[nquadelems].idx1 = quadelems[nquadelems].idx2;
            quadelems[nquadelems].idx2 = tmp;
         }
      }

      if( *doingfine )
      {
         SCIP_EXPR** children;

         /* setup array with children expressions corresponding to variables */
         SCIP_CALL( SCIPallocBufferArray(scip, &children, nquadvars) );
         for( i = 0; i < nvars; ++i )
         {
            if( quadvarsidxs[i] == -1 )
               continue;

            /* assign new index to variable, if we see it the first time in this exprtree */
            if( exprvaridx[i] == -1 )
            {
               exprvaridx[i] = *nexprvars;
               ++*nexprvars;
            }

            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &children[quadvarsidxs[i]], SCIP_EXPR_VARIDX, exprvaridx[i]) );
         }

         SCIP_CALL( SCIPexprCreateQuadratic(SCIPblkmem(scip), expr, nquadvars, children, 0.0, NULL, nquadelems, quadelems) );

         SCIPfreeBufferArray(scip, &children);
      }

      SCIPfreeBufferArray(scip, &quadelems);
      SCIPfreeBufferArray(scip, &quadvarsidxs);
   }


   SCIPerrorMessage("Expression operand <%s> not supported by SCIP so far.\n", exprname);
   *doingfine = FALSE;

   return SCIP_OKAY;
}


/** read nonlinear expressions of constraints and objective */
static
SCIP_RETCODE readNonlinearExprs(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       datanode,           /**< XML root node for instance data */
   SCIP_VAR**            vars,               /**< variables in order of OSiL indices */
   int                   nvars,              /**< number of variables */
   SCIP_CONS**           conss,              /**< constraints in order of OSiL indices */
   CONSTYPE*             constypes,          /**< type of constraints (assumed to be LINEAR) */
   int                   nconss,             /**< number of constraints */
   SCIP_CONS**           objcons,            /**< buffer to store constraint for nonlinear part of objective function, or to add to if already existing */
   CONSTYPE*             objconstype,        /**< buffer to store type of objective constraint, if created (should be QUADRATIC) */
   SCIP_Bool*            doingfine           /**< buffer to indicate whether no errors occurred */
)
{
   const XML_NODE* nlexprs;
   const XML_NODE* nlexpr;
   const char* attrval;
   SCIP_EXPRTREE* exprtree;
   SCIP_EXPR* expr;
   SCIP_VAR** exprvars;
   int* exprvaridx;
   int nexprvars;
   int nnlexprs;
   int count;
   int considx;
   int i;

   assert(scip != NULL);
   assert(datanode != NULL);
   assert(vars != NULL || nvars == 0);
   assert(conss != NULL || nconss == 0);
   assert(constypes != NULL || nconss == 0);
   assert(objcons != NULL);
   assert(doingfine != NULL);

   nlexprs = xmlFindNodeMaxdepth(datanode, "nonlinearExpressions", 0, 1);

   if( nlexprs == NULL )
      return SCIP_OKAY;

   attrval = xmlGetAttrval(nlexprs, "numberOfNonlinearExpressions");
   if( attrval == NULL )
   {
      SCIPerrorMessage("numberOfNonlinearExpressions attribute for nonlinearExpressions not found.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }

   nnlexprs = (int)strtol(attrval, (char**)&attrval, 10);
   if( *attrval != '\0' || nnlexprs < 0 )
   {
      SCIPerrorMessage("Invalid value for numberOfNonlinearExpressions attribute.\n");
      *doingfine = FALSE;
      return SCIP_OKAY;
   }
   assert(nnlexprs >= 0);

   /* buffer array to store index of variable in expression graph, or -1 if not present */
   SCIP_CALL( SCIPallocBufferArray(scip, &exprvaridx, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &exprvars, nvars) );

   count = 0;
   for( nlexpr = xmlFirstChild(nlexprs); nlexpr != NULL; nlexpr = xmlNextSibl(nlexpr), ++count )
   {
      if( strcmp(xmlGetName(nlexpr), "nl") != 0 )
      {
         SCIPerrorMessage("expected <nl>-node, but got %s\n", xmlGetName(nlexpr));
         *doingfine = FALSE;
         break;
      }
      if( count >= nnlexprs )
      {
         SCIPerrorMessage("too many nonlinear expressions, expected %d, got at least %d.\n", nnlexprs, count + 1);
         *doingfine = FALSE;
         break;
      }

      /* treat empty expression as 0.0 and continue */
      if( xmlFirstChild(nlexprs) == NULL )
         continue;

      /* get constraint index, or -1 for objective */
      attrval = xmlGetAttrval(nlexpr, "idx");
      if( attrval == NULL )
      {
         SCIPerrorMessage("Missing 'idx' attribute in <nl>.\n");
         *doingfine = FALSE;
         break;
      }

      considx = (int)strtol(attrval, (char**)&attrval, 10);
      if( *attrval != '\0' || considx < -1 || considx >= nconss )
      {
         SCIPerrorMessage("Invalid value '%s' in idx attribute of <nl> node.\n", xmlGetAttrval(nlexpr, "idx"));
         *doingfine = FALSE;
         break;
      }

      expr = NULL;
      nexprvars = 0;
      for( i = 0; i < nvars; ++i )
         exprvaridx[i] = -1;

      /* turn OSiL expression into SCIP expression and assign indices to variables */
      SCIP_CALL( readExpression(scip, &expr, xmlFirstChild(nlexpr), exprvaridx, &nexprvars, nvars, doingfine) );
      if( !*doingfine )
      {
         assert(expr == NULL);
         break;
      }

      /* assemble array with SCIP_VAR*'s */
      for( i = 0; i < nvars; ++i )
      {
         assert(exprvaridx[i] < nexprvars );

         if( exprvaridx[i] >= 0 )
            exprvars[exprvaridx[i]] = vars[i];
      }

      /* create expression tree */
      SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, expr, nexprvars, 0, NULL) );
      SCIP_CALL( SCIPexprtreeSetVars(exprtree, nexprvars, exprvars) );

      /* add expression tree to objective or constraint */
      if( considx == -1 && *objcons == NULL )
      {
         /* create constraint to hold nonlinear part of objective */
         SCIP_VAR* objvar;
         SCIP_Real minusone;
         SCIP_Real one;

         SCIP_CALL( SCIPcreateVar(scip, &objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, objvar) );

         minusone = -1.0;
         one = 1.0;
         SCIP_CALL( SCIPcreateConsNonlinear(scip, objcons, "objcons", 1, &objvar, &minusone, 1, &exprtree, &one,
            SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE ? -SCIPinfinity(scip) : 0.0,
            SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE ?  SCIPinfinity(scip) : 0.0,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         *objconstype = NONLINEAR;

         SCIP_CALL( SCIPreleaseVar(scip, &objvar) );
      }
      else
      {
         SCIP_CONS** cons;
         SCIP_CONS* oldcons;
         CONSTYPE* constype;

         if( considx == -1 )
         {
            cons = objcons;
            constype = objconstype;
         }
         else
         {
            cons = &conss[considx];
            constype = &constypes[considx];
         }
         oldcons = *cons;

         /* replace cons by nonlinear constraint or add to already existing nonlinear constraint */
         switch( *constype )
         {
         case LINEAR:
         {
            SCIP_Real one;

            one = 1.0;
            SCIP_CALL( SCIPcreateConsNonlinear(scip, cons, SCIPconsGetName(*cons),
               SCIPgetNVarsLinear(scip, *cons), SCIPgetVarsLinear(scip, *cons), SCIPgetValsLinear(scip, *cons),
               1, &exprtree, &one,
               SCIPgetLhsLinear(scip, *cons), SCIPgetRhsLinear(scip, *cons),
               SCIPconsIsInitial(*cons), SCIPconsIsSeparated(*cons), SCIPconsIsEnforced(*cons),
               SCIPconsIsChecked(*cons), SCIPconsIsPropagated(*cons), SCIPconsIsLocal(*cons),
               SCIPconsIsModifiable(*cons), SCIPconsIsDynamic(*cons), SCIPconsIsRemovable(*cons), SCIPconsIsStickingAtNode(*cons)) );

            SCIP_CALL( SCIPreleaseCons(scip, &oldcons) );

            break;
         }

         case QUADRATIC:
         {
            SCIP_EXPRTREE* exprtrees[2];
            SCIP_Real exprcoefs[2];

            SCIP_EXPR* quadexpr;
            SCIP_QUADELEM* quadelems;
            SCIP_Real* lincoefs;
            SCIP_EXPR** children;
            SCIP_QUADVARTERM* quadvarterms;
            SCIP_BILINTERM* bilinterms;
            int nquadelems;
            int nquadvars;
            int nbilin;
            int j;

            exprtrees[0] = exprtree;
            exprcoefs[0] = 1.0;

            /* turn quadratic part into expression tree */
            SCIP_CALL( SCIPsortQuadVarTermsQuadratic(scip, *cons) );

            quadvarterms = SCIPgetQuadVarTermsQuadratic(scip, *cons);
            nquadvars = SCIPgetNQuadVarTermsQuadratic(scip, *cons);
            bilinterms = SCIPgetBilinTermsQuadratic(scip, *cons);
            nbilin = SCIPgetNBilinTermsQuadratic(scip, *cons);

            SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, nquadvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &children, nquadvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, nbilin + nquadvars) );
            nquadelems = 0;

            for( i = 0; i < nquadvars; ++i )
            {
               lincoefs[i] = quadvarterms[i].lincoef;
               exprvars[i] = quadvarterms[i].var;

               if( quadvarterms[i].sqrcoef != 0.0 )
               {
                  quadelems[nquadelems].idx1 = i;
                  quadelems[nquadelems].idx2 = i;
                  quadelems[nquadelems].coef = quadvarterms[i].sqrcoef;
                  ++nquadelems;
               }

               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &children[i], SCIP_EXPR_VARIDX, i) );

               for( j = 0; j < quadvarterms[i].nadjbilin; ++j )
               {
                  if( bilinterms[quadvarterms[i].adjbilin[j]].var1 == quadvarterms[i].var )
                  {
                     int otheridx;

                     assert(bilinterms[quadvarterms[i].adjbilin[j]].var2 != quadvarterms[i].var);

                     SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, *cons, bilinterms[quadvarterms[i].adjbilin[j]].var2, &otheridx) );
                     assert(otheridx >= 0);
                     assert(otheridx < nquadvars);

                     quadelems[nquadelems].idx1 = MIN(i, otheridx);
                     quadelems[nquadelems].idx2 = MAX(i, otheridx);
                     quadelems[nquadelems].coef = bilinterms[quadvarterms[i].adjbilin[j]].coef;
                     ++nquadelems;
                  }
               }
            }

            SCIP_CALL( SCIPexprCreateQuadratic(SCIPblkmem(scip), &quadexpr, nquadvars, children, 0.0, lincoefs, nquadelems, quadelems) );
            SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtrees[1], quadexpr, nquadvars, 0, NULL) );
            SCIP_CALL( SCIPexprtreeSetVars(exprtrees[1], nquadvars, exprvars) );
            exprcoefs[1] = 1.0;

            SCIPfreeBufferArray(scip, &lincoefs);
            SCIPfreeBufferArray(scip, &children);
            SCIPfreeBufferArray(scip, &quadelems);

            SCIP_CALL( SCIPcreateConsNonlinear(scip, cons, SCIPconsGetName(*cons),
               SCIPgetNLinearVarsNonlinear(scip, *cons), SCIPgetLinearVarsNonlinear(scip, *cons), SCIPgetLinearCoefsNonlinear(scip, *cons),
               2, exprtrees, exprcoefs,
               SCIPgetLhsNonlinear(scip, *cons), SCIPgetRhsNonlinear(scip, *cons),
               SCIPconsIsInitial(*cons), SCIPconsIsSeparated(*cons), SCIPconsIsEnforced(*cons),
               SCIPconsIsChecked(*cons), SCIPconsIsPropagated(*cons), SCIPconsIsLocal(*cons),
               SCIPconsIsModifiable(*cons), SCIPconsIsDynamic(*cons), SCIPconsIsRemovable(*cons), SCIPconsIsStickingAtNode(*cons)) );

            SCIP_CALL( SCIPreleaseCons(scip, &oldcons) );

            break;
         }

         case NONLINEAR:
         {
            SCIP_Real one;

            one = 1.0;
            SCIP_CALL( SCIPaddExprtreesNonlinear(scip, *cons, 1, &exprtree, &one) );
            break;
         }
         }

         *constype = NONLINEAR;
      }

      SCIP_CALL( SCIPexprtreeFree(&exprtree) );
   }

   SCIPfreeBufferArray(scip, &exprvars);
   SCIPfreeBufferArray(scip, &exprvaridx);

   return SCIP_OKAY;
}

/*
 * Callback methods of reader
 */


/** copy method for reader plugins (called when SCIP copies plugins) */
#if 1
static
SCIP_DECL_READERCOPY(readerCopyOsil)
{  /*lint --e{715}*/
   assert(scip != NULL);

   SCIP_CALL( SCIPincludeReaderOsil(scip) );

   return SCIP_OKAY;
}
#else
#define readerCopyOsil NULL
#endif

/** destructor of reader to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_READERFREE(readerFreeOsil)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of osil reader not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define readerFreeOsil NULL
#endif


/** problem reading method of reader */
#if 1
static
SCIP_DECL_READERREAD(readerReadOsil)
{  /*lint --e{715}*/
   const char* name;
   XML_NODE* start;
   const XML_NODE* header;
   const XML_NODE* data;
   SCIP_RETCODE retcode;
   SCIP_Bool doingfine;
   SCIP_VAR** vars;
   int nvars;
   SCIP_CONS** conss;
   CONSTYPE* constypes;
   int nconss;
   SCIP_CONS* objcons;
   CONSTYPE objconstype;
   int i;

   assert(scip != NULL);
   assert(reader != NULL);
   assert(result != NULL);
   assert(filename != NULL);

   *result = SCIP_DIDNOTRUN;
   retcode = SCIP_READERROR;
   doingfine = TRUE;
   vars = NULL;
   nvars = 0;
   conss = NULL;
   constypes = NULL;
   nconss = 0;
   objcons = NULL;

   /* read OSiL xml file */
   start = xmlProcess(filename);

   if( start == NULL )
   {
      SCIPerrorMessage("Some error occured during parsing the OSiL XML file.\n");
      goto CLEANUP;
   }

   SCIPdebug( xmlShowNode(start) );

   /* parse header to get problem name */
   name = filename;
   header = xmlFindNodeMaxdepth(start, "instanceHeader", 0, 2);
   if( header != NULL )
   {
      const XML_NODE* namenode;

      namenode = xmlFindNodeMaxdepth(header, "name", 0, 2);

      if( namenode != NULL && xmlFirstChild(namenode) != NULL )
         name = xmlGetData(xmlFirstChild(namenode));
   }

   /* create SCIP problem */
   SCIP_CALL( SCIPcreateProb(scip, name, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* process instance data */
   data = xmlFindNodeMaxdepth(start, "instanceData", 0, 2);
   if( data == NULL )
   {
      SCIPerrorMessage("Node <instanceData> not found.\n");
      goto CLEANUP;
   }

   /* read variables */
   SCIP_CALL( readVariables(scip, data, &vars, &nvars, &doingfine) );
   if( !doingfine )
      goto CLEANUP;
   assert(vars != NULL || nvars == 0);

   /* read objective sense, coefficient, and constant */
   SCIP_CALL( readObjective(scip, data, vars, nvars, &doingfine) );
   if( !doingfine )
      goto CLEANUP;

   /* read constraint data (names, constants, lhs/rhs) */
   SCIP_CALL( readConstraints(scip, data, &conss, &constypes, &nconss, &doingfine) );
   if( !doingfine )
      goto CLEANUP;
   assert(conss != NULL || nconss == 0);

   /* read linear coefficients matrix */
   SCIP_CALL( readLinearCoefs(scip, data, vars, nvars, conss, constypes, nconss, &doingfine) );
   if( !doingfine )
      goto CLEANUP;

   /* read quadratic coefficients (turns linear constraints into quadratic ones, may create objcons) */
   SCIP_CALL( readQuadraticCoefs(scip, data, vars, nvars, conss, constypes, nconss, &objcons, &objconstype, &doingfine) );
   if( !doingfine )
      goto CLEANUP;

   /* read nonlinear expressions (turns constraints into nonlinear ones, may create objcons) */
   SCIP_CALL( readNonlinearExprs(scip, data, vars, nvars, conss, constypes, nconss, &objcons, &objconstype, &doingfine) );
   if( !doingfine )
      goto CLEANUP;

   /* add constraints to problem */
   for( i = 0; i < nconss; ++i )
   {
      assert(conss[i] != NULL);
      SCIP_CALL( SCIPaddCons(scip, conss[i]) );
   }
   if( objcons != NULL )
   {
      SCIP_CALL( SCIPaddCons(scip, objcons) );
   }

   *result = SCIP_SUCCESS;
   retcode = SCIP_OKAY;

CLEANUP:
   /* free xml data */
   if( start != NULL )
      xmlFreeNode(start);

   /* free variables */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &vars[i]) );
   }
   SCIPfreeBufferArrayNull(scip, &vars);

   /* free constraints */
   for( i = 0; i < nconss; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &conss[i]) );
   }
   SCIPfreeBufferArrayNull(scip, &conss);
   SCIPfreeBufferArrayNull(scip, &constypes);

   if( objcons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &objcons) );
   }

   return retcode;
}
#else
#define readerReadOsil NULL
#endif


#if 0
/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteOsil)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of osil reader not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define readerWriteOsil NULL
#endif


/*
 * reader specific interface methods
 */

/** includes the osil file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderOsil(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   /* include osil reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerCopyOsil, readerFreeOsil, readerReadOsil, readerWriteOsil, NULL) );

   return SCIP_OKAY;
}

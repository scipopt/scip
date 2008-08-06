/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2008 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: reader_mps.c,v 1.94 2008/08/06 09:20:09 bzfwolte Exp $"

//#define SCIP_DEBUG

/**@file   reader_mps.c
 * @brief  (extended) MPS file reader
 * @author Thorsten Koch
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 * @author Stefan Heinz
 *
 * This reader/writer handles MPS files in extended MPS format, as it
 * is used by CPLEX. In the extended format the limits on variable
 * name lengths and coefficients are considerably relaxed. The columns
 * in the format are then separated by whitespaces.
 *
 * @todo Test for uniqueness of variable and constraint names (after cutting down).
 * @todo Check whether constructing the names for aggregated constraint yields name clashes (aggrXXX).
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>

#include "scip/reader_mps.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"

#define READER_NAME             "mpsreader"
#define READER_DESC             "file reader for MIPs in IBM's Mathematical Programming System format"
#define READER_EXTENSION        "mps"



/*
 * mps reader internal methods
 */

#define MPS_MAX_LINELEN  1024
#define MPS_MAX_NAMELEN   256
#define MPS_MAX_VALUELEN   26
#define MPS_MAX_FIELDLEN   20

#define PATCH_CHAR    '_'
#define BLANK         ' '

enum MpsSection
   {
      MPS_NAME,
      MPS_OBJSEN,
      MPS_OBJNAME,
      MPS_ROWS,
      MPS_USERCUTS,
      MPS_LAZYCONS,
      MPS_COLUMNS,
      MPS_RHS,
      MPS_RANGES,
      MPS_BOUNDS,
      MPS_SOS,
      MPS_ENDATA
   };
typedef enum MpsSection MPSSECTION;

struct MpsInput
{
   MPSSECTION           section;
   SCIP_FILE*           fp;
   int                  lineno;
   SCIP_OBJSENSE        objsense;
   SCIP_Bool            haserror;
   char                 buf[MPS_MAX_LINELEN];
   const char*          f0;
   const char*          f1;
   const char*          f2;
   const char*          f3;
   const char*          f4;
   const char*          f5;
   char                 probname[MPS_MAX_NAMELEN];
   char                 objname [MPS_MAX_NAMELEN];
   SCIP_Bool            isinteger;
   SCIP_Bool            isnewformat;
   SCIP_Bool            semicontwarning;
};
typedef struct MpsInput MPSINPUT;

/* sparse matrix representation */
struct SparseMatrix
{
   SCIP_Real*            values;           /**< matrix element */
   SCIP_VAR**            columns;          /**< corresponding variables */
   const char**          rows;             /**< corresponding constraint names */ 
   int                   nentries;         /**< number of elements in the arrays */
   int                   sentries;         /**< number of slots in the arrays */
};
typedef struct SparseMatrix SPARSEMATRIX;


static
SCIP_RETCODE mpsinputCreate(
   SCIP*                 scip,
   MPSINPUT**            mpsi,
   SCIP_FILE*            fp
   )
{
   assert(mpsi != NULL);
   assert(fp != NULL);

   SCIP_CALL( SCIPallocMemory(scip, mpsi) );

   (*mpsi)->section     = MPS_NAME;
   (*mpsi)->fp          = fp;
   (*mpsi)->lineno      = 0;
   (*mpsi)->objsense    = SCIP_OBJSENSE_MINIMIZE;
   (*mpsi)->haserror    = FALSE;
   (*mpsi)->isinteger   = FALSE;
   (*mpsi)->isnewformat = FALSE;
   (*mpsi)->semicontwarning = FALSE;
   (*mpsi)->buf     [0] = '\0';
   (*mpsi)->probname[0] = '\0';
   (*mpsi)->objname [0] = '\0';
   (*mpsi)->f0          = NULL;
   (*mpsi)->f1          = NULL;
   (*mpsi)->f2          = NULL;
   (*mpsi)->f3          = NULL;
   (*mpsi)->f4          = NULL;
   (*mpsi)->f5          = NULL;

   return SCIP_OKAY;
}

static
void mpsinputFree(
   SCIP*                 scip,
   MPSINPUT**            mpsi
   )
{
   SCIPfreeMemory(scip, mpsi);
}

static
MPSSECTION mpsinputSection(
   const MPSINPUT*       mpsi
   )
{
   assert(mpsi != NULL);

   return mpsi->section;
}

#if 0
static
int mpsinputLineno(
   const MPSINPUT*       mpsi
   )
{
   assert(mpsi != NULL);

   return mpsi->lineno;
}
#endif

static
const char* mpsinputField0(
   const MPSINPUT*       mpsi
   )
{
   assert(mpsi != NULL);

   return mpsi->f0;
}

static
const char* mpsinputField1(
   const MPSINPUT*       mpsi
   )
{
   assert(mpsi != NULL);

   return mpsi->f1;
}

static
const char* mpsinputField2(
   const MPSINPUT*       mpsi
   )
{
   assert(mpsi != NULL);

   return mpsi->f2;
}

static
const char* mpsinputField3(
   const MPSINPUT*       mpsi
   )
{
   assert(mpsi != NULL);

   return mpsi->f3;
}

static
const char* mpsinputField4(
   const MPSINPUT*       mpsi
   )
{
   assert(mpsi != NULL);

   return mpsi->f4;
}

static
const char* mpsinputField5(
   const MPSINPUT*       mpsi
   )
{
   assert(mpsi != NULL);

   return mpsi->f5;
}

#if 0
static
const char* mpsinputProbname(
   const MPSINPUT*       mpsi
   )
{
   assert(mpsi != NULL);

   return mpsi->probname;
}
#endif

static
const char* mpsinputObjname(
   const MPSINPUT*       mpsi
   )
{
   assert(mpsi != NULL);

   return mpsi->objname;
}

static
SCIP_OBJSENSE mpsinputObjsense(
   const MPSINPUT*       mpsi
   )
{
   assert(mpsi != NULL);

   return mpsi->objsense;
}

static
SCIP_Bool mpsinputHasError(
   const MPSINPUT*       mpsi
   )
{
   assert(mpsi != NULL);

   return mpsi->haserror;
}

static
SCIP_Bool mpsinputIsInteger(
   const MPSINPUT*       mpsi
   )
{
   assert(mpsi != NULL);

   return mpsi->isinteger;
}

static
void mpsinputSetSection(
   MPSINPUT*             mpsi,
   MPSSECTION            section
   )
{
   assert(mpsi != NULL);

   mpsi->section = section;
}

static
void mpsinputSetProbname(
   MPSINPUT*             mpsi,
   const char*           probname
   )
{
   assert(mpsi     != NULL);
   assert(probname != NULL);
   assert(strlen(probname) < sizeof(mpsi->probname));

   strncpy(mpsi->probname, probname, MPS_MAX_NAMELEN - 1);
}

static
void mpsinputSetObjname(
   MPSINPUT*             mpsi,
   const char*           objname
   )
{
   assert(mpsi != NULL);
   assert(objname != NULL);
   assert(strlen(objname) < sizeof(mpsi->objname));

   strncpy(mpsi->objname, objname, MPS_MAX_NAMELEN - 1);
}

static
void mpsinputSetObjsense(
   MPSINPUT*             mpsi,
   SCIP_OBJSENSE         sense
   )
{
   assert(mpsi != NULL);

   mpsi->objsense = sense;
}

static
void mpsinputSyntaxerror(
   MPSINPUT*             mpsi
   )
{
   assert(mpsi != NULL);

   SCIPwarningMessage("Syntax error in line %d\n", mpsi->lineno);
   mpsi->section  = MPS_ENDATA;
   mpsi->haserror = TRUE;
}

static
void mpsinputEntryIgnored(
   SCIP*                 scip,               /**< SCIP data structure */
   MPSINPUT*             mpsi,
   const char*           what,
   const char*           what_name,
   const char*           entity,
   const char*           entity_name,
   SCIP_VERBLEVEL        verblevel
   )
{
   assert(mpsi        != NULL);
   assert(what        != NULL);
   assert(what_name   != NULL);
   assert(entity      != NULL);
   assert(entity_name != NULL);

   SCIPverbMessage(scip, verblevel, NULL,
      "Warning line %d: %s \"%s\" for %s \"%s\" ignored\n", mpsi->lineno, what, what_name, entity, entity_name);
}

/* fill the line from \p pos up to column 80 with blanks. */
static
void clearFrom(
   char*                 buf,
   unsigned int          pos
   )
{
   unsigned int i;

   for(i = pos; i < 80; i++)
      buf[i] = BLANK;
   buf[80] = '\0';
}

/* change all blanks inside a field to #PATCH_CHAR. */
static
void patchField(
   char*                 buf,
   int                   beg,
   int                   end
   )
{
   int i;

   while((beg <= end) && (buf[end] == BLANK))
      end--;

   while((beg <= end) && (buf[beg] == BLANK))
      beg++;

   for(i = beg; i <= end; i++)
      if (buf[i] == BLANK)
         buf[i] = PATCH_CHAR;
}

/* read a mps format data line and parse the fields. */
static
SCIP_Bool mpsinputReadLine(
   MPSINPUT*             mpsi
   )
{
   unsigned int len;
   unsigned int i;
   int space;
   char* s;
   SCIP_Bool is_marker;
   SCIP_Bool is_empty;
   char* nexttok;

   do
   {
      mpsi->f0 = mpsi->f1 = mpsi->f2 = mpsi->f3 = mpsi->f4 = mpsi->f5 = 0;
      is_marker = FALSE;
      is_empty = FALSE;

      /* Read until we have a not comment line. */
      do
      {
         if (NULL == SCIPfgets(mpsi->buf, sizeof(mpsi->buf), mpsi->fp))
            return FALSE;
         mpsi->lineno++;
      }
      while(*mpsi->buf == '*');

      /* Normalize line */
      len = strlen(mpsi->buf);

      for(i = 0; i < len; i++)
         if ((mpsi->buf[i] == '\t') || (mpsi->buf[i] == '\n') || (mpsi->buf[i] == '\r'))
            mpsi->buf[i] = BLANK;

      if (len < 80)
         clearFrom(mpsi->buf, len);

      assert(strlen(mpsi->buf) >= 80);

      /* Look for new section */
      if (*mpsi->buf != BLANK)
      {
         mpsi->f0 = SCIPstrtok(&mpsi->buf[0], " ", &nexttok);

         assert(mpsi->f0 != 0);

         mpsi->f1 = SCIPstrtok(NULL, " ", &nexttok);

         return TRUE;
      }

      /* If we decide to use the new format we never revert this decision */
      if (!mpsi->isnewformat)
      {
         /* Test for fixed format comments */
         if ((mpsi->buf[14] == '$') && (mpsi->buf[13] == ' '))
            clearFrom(mpsi->buf, 14);
         else if ((mpsi->buf[39] == '$') && (mpsi->buf[38] == ' '))
            clearFrom(mpsi->buf, 39);

         /* Test for fixed format */
         space = mpsi->buf[12] | mpsi->buf[13]
            | mpsi->buf[22] | mpsi->buf[23]
            | mpsi->buf[36] | mpsi->buf[37] | mpsi->buf[38]
            | mpsi->buf[47] | mpsi->buf[48]
            | mpsi->buf[61] | mpsi->buf[62] | mpsi->buf[63];

         if (space == BLANK)
         {
            /* Now we have space at the right positions.
             * But are there also the non space where they
             * should be ?
             */
            SCIP_Bool number = isdigit(mpsi->buf[24]) || isdigit(mpsi->buf[25])
               || isdigit(mpsi->buf[26]) || isdigit(mpsi->buf[27])
               || isdigit(mpsi->buf[28]) || isdigit(mpsi->buf[29])
               || isdigit(mpsi->buf[30]) || isdigit(mpsi->buf[31])
               || isdigit(mpsi->buf[32]) || isdigit(mpsi->buf[33])
               || isdigit(mpsi->buf[34]) || isdigit(mpsi->buf[35]);

            /* len < 13 is handle ROW lines with embedded spaces
             * in the names correctly
             */
            if (number || len < 13)
            {
               /* We assume fixed format, so we patch possible embedded spaces. */
               patchField(mpsi->buf,  4, 12);
               patchField(mpsi->buf, 14, 22);
               patchField(mpsi->buf, 39, 47);
            }
            else
            {
               if (  mpsi->section == MPS_COLUMNS || mpsi->section == MPS_RHS
                  || mpsi->section == MPS_RANGES  || mpsi->section == MPS_BOUNDS)
                  mpsi->isnewformat = TRUE;
            }
         }
         else
         {
            mpsi->isnewformat = TRUE;
         }
      }
      s = &mpsi->buf[1];

      /* At this point it is not clear if we have a indicator field.
       * If there is none (e.g. empty) f1 will be the first name field.
       * If there is one, f2 will be the first name field.
       *
       * Initially comment marks '$' are only allowed in the beginning
       * of the 2nd and 3rd name field. We test all fields but the first.
       * This makes no difference, since if the $ is at the start of a value
       * field, the line will be errornous anyway.
       */
      do
      {
         if (NULL == (mpsi->f1 = SCIPstrtok(s, " ", &nexttok)))
            break;

         if ((NULL == (mpsi->f2 = SCIPstrtok(NULL, " ", &nexttok))) || (*mpsi->f2 == '$'))
         {
            mpsi->f2 = 0;
            break;
         }
         if (!strcmp(mpsi->f2, "'MARKER'"))
            is_marker = TRUE;

         if ((NULL == (mpsi->f3 = SCIPstrtok(NULL, " ", &nexttok))) || (*mpsi->f3 == '$'))
         {
            mpsi->f3 = 0;
            break;
         }
         if (is_marker)
         {
            if (!strcmp(mpsi->f3, "'INTORG'"))
               mpsi->isinteger = TRUE;
            else if (!strcmp(mpsi->f3, "'INTEND'"))
               mpsi->isinteger = FALSE;
            else
               break; /* unknown marker */
         }
         if (!strcmp(mpsi->f3, "'MARKER'"))
            is_marker = TRUE;

         if ((NULL == (mpsi->f4 = SCIPstrtok(NULL, " ", &nexttok))) || (*mpsi->f4 == '$'))
         {
            mpsi->f4 = 0;
            break;
         }
         if (is_marker)
         {
            if (!strcmp(mpsi->f4, "'INTORG'"))
               mpsi->isinteger = TRUE;
            else if (!strcmp(mpsi->f4, "'INTEND'"))
               mpsi->isinteger = FALSE;
            else
               break; /* unknown marker */
         }
         if ((NULL == (mpsi->f5 = SCIPstrtok(NULL, " ", &nexttok))) || (*mpsi->f5 == '$'))
            mpsi->f5 = 0;
      }
      while(FALSE);

      /* check for empty lines */
      is_empty = (mpsi->f0 == NULL && mpsi->f1 == NULL);
   }
   while(is_marker || is_empty);

   return TRUE;
}

/* Insert \p name as field 1 and shift all other fields up. */
static
void mpsinputInsertName(
   MPSINPUT*             mpsi,
   const char*           name,
   SCIP_Bool             second
   )
{
   assert(mpsi != NULL);
   assert(name != NULL);

   mpsi->f5 = mpsi->f4;
   mpsi->f4 = mpsi->f3;
   mpsi->f3 = mpsi->f2;

   if (second)
      mpsi->f2 = name;
   else
   {
      mpsi->f2 = mpsi->f1;
      mpsi->f1 = name;
   }
}

/* Process NAME section. */
static
SCIP_RETCODE readName(
   MPSINPUT*             mpsi
   )
{
   assert(mpsi != NULL);

   /* This has to be the Line with the NAME section. */
   if (!mpsinputReadLine(mpsi) || mpsinputField0(mpsi) == NULL || strcmp(mpsinputField0(mpsi), "NAME"))
   {
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }

   /* Sometimes the name is omitted. */
   mpsinputSetProbname(mpsi, (mpsinputField1(mpsi) == 0) ? "_MPS_" : mpsinputField1(mpsi));

   /*printf("Problem name   : %s\n", mpsinputProbname(mpsi));*/

   /* This hat to be a new section */
   if (!mpsinputReadLine(mpsi) || (mpsinputField0(mpsi) == NULL))
   {
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }

   if (!strncmp(mpsinputField0(mpsi), "ROWS", 4))
      mpsinputSetSection(mpsi, MPS_ROWS);
   else if (!strncmp(mpsinputField0(mpsi), "USERCUTS", 8))
      mpsinputSetSection(mpsi, MPS_USERCUTS);
   else if (!strncmp(mpsinputField0(mpsi), "LAZYCONS", 8))
      mpsinputSetSection(mpsi, MPS_LAZYCONS);
   else if (!strncmp(mpsinputField0(mpsi), "OBJSEN", 6))
      mpsinputSetSection(mpsi, MPS_OBJSEN);
   else if (!strncmp(mpsinputField0(mpsi), "OBJNAME", 7))
      mpsinputSetSection(mpsi, MPS_OBJNAME);
   else
   {
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/* Process OBJSEN section. This Section is an ILOG extension. */
static
SCIP_RETCODE readObjsen(
   MPSINPUT*             mpsi
   )
{
   assert(mpsi != NULL);

   /* This has to be the Line with MIN or MAX. */
   if (!mpsinputReadLine(mpsi) || (mpsinputField1(mpsi) == NULL))
   {
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }

   if (!strncmp(mpsinputField1(mpsi), "MIN", 3))
      mpsinputSetObjsense(mpsi, SCIP_OBJSENSE_MINIMIZE);
   else if (!strncmp(mpsinputField1(mpsi), "MAX", 3))
      mpsinputSetObjsense(mpsi, SCIP_OBJSENSE_MAXIMIZE);
   else
   {
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }

   /* Look for ROWS, USERCUTS, LAZYCONS, or OBJNAME Section */
   if (!mpsinputReadLine(mpsi) || mpsinputField0(mpsi) == NULL)
   {
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }

   if (!strcmp(mpsinputField0(mpsi), "ROWS"))
      mpsinputSetSection(mpsi, MPS_ROWS);
   else if (!strcmp(mpsinputField0(mpsi), "USERCUTS"))
      mpsinputSetSection(mpsi, MPS_USERCUTS);
   else if (!strcmp(mpsinputField0(mpsi), "LAZYCONS"))
      mpsinputSetSection(mpsi, MPS_LAZYCONS);
   else if (!strcmp(mpsinputField0(mpsi), "OBJNAME"))
      mpsinputSetSection(mpsi, MPS_OBJNAME);
   else
   {
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/* Process OBJNAME section. This Section is an ILOG extension. */
static
SCIP_RETCODE readObjname(
   MPSINPUT*             mpsi
   )
{
   assert(mpsi != NULL);

   /* This has to be the Line with the name. */
   if (!mpsinputReadLine(mpsi) || mpsinputField1(mpsi) == NULL)
   {
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }

   mpsinputSetObjname(mpsi, mpsinputField1(mpsi));

   /* Look for ROWS, USERCUTS, or LAZYCONS Section */
   if (!mpsinputReadLine(mpsi) || mpsinputField0(mpsi) == NULL)
   {
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }
   if (!strcmp(mpsinputField0(mpsi), "ROWS"))
      mpsinputSetSection(mpsi, MPS_ROWS);
   else if (!strcmp(mpsinputField0(mpsi), "USERCUTS"))
      mpsinputSetSection(mpsi, MPS_USERCUTS);
   else if (!strcmp(mpsinputField0(mpsi), "LAZYCONS"))
      mpsinputSetSection(mpsi, MPS_LAZYCONS);
   else
      mpsinputSyntaxerror(mpsi);

   return SCIP_OKAY;
}

/* Process ROWS, USERCUTS, or LAZYCONS section. */
static
SCIP_RETCODE readRows(
   MPSINPUT*             mpsi,
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   while(mpsinputReadLine(mpsi))
   {
      if (mpsinputField0(mpsi) != NULL)
      {
         if (!strcmp(mpsinputField0(mpsi), "ROWS"))
            mpsinputSetSection(mpsi, MPS_ROWS);
         else if (!strcmp(mpsinputField0(mpsi), "USERCUTS"))
            mpsinputSetSection(mpsi, MPS_USERCUTS);
         else if (!strcmp(mpsinputField0(mpsi), "LAZYCONS"))
            mpsinputSetSection(mpsi, MPS_LAZYCONS);
         else if (!strcmp(mpsinputField0(mpsi), "COLUMNS"))
            mpsinputSetSection(mpsi, MPS_COLUMNS);
         else
            mpsinputSyntaxerror(mpsi);

         return SCIP_OKAY;
      }

      if (*mpsinputField1(mpsi) == 'N')
      {
         if (*mpsinputObjname(mpsi) == '\0')
            mpsinputSetObjname(mpsi, mpsinputField2(mpsi));
         else
            mpsinputEntryIgnored(scip, mpsi, "row", mpsinputField2(mpsi), "objective function", "N", SCIP_VERBLEVEL_NORMAL);
      }
      else
      {
         SCIP_CONS* cons;
         SCIP_Bool dynamicrows;
         SCIP_Bool dynamicconss;
         SCIP_Bool initial;
         SCIP_Bool separate;
         SCIP_Bool enforce;
         SCIP_Bool check;
         SCIP_Bool propagate;
         SCIP_Bool local;
         SCIP_Bool modifiable;
         SCIP_Bool dynamic;
         SCIP_Bool removable;

         cons = SCIPfindCons(scip, mpsinputField2(mpsi));
         if( cons != NULL )
            break;

         SCIP_CALL( SCIPgetBoolParam(scip, "reading/mpsreader/dynamicconss", &dynamicconss) );
         SCIP_CALL( SCIPgetBoolParam(scip, "reading/mpsreader/dynamicrows", &dynamicrows) );
         initial = !dynamicrows && (mpsinputSection(mpsi) == MPS_ROWS);
         separate = TRUE;
         enforce = (mpsinputSection(mpsi) != MPS_USERCUTS);
         check = (mpsinputSection(mpsi) != MPS_USERCUTS);
         propagate = TRUE;
         local = FALSE;
         modifiable = FALSE;
         dynamic = dynamicconss;
         removable = dynamicrows || (mpsinputSection(mpsi) == MPS_USERCUTS);

         switch(*mpsinputField1(mpsi))
         {
         case 'G' :
            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, mpsinputField2(mpsi), 0, NULL, NULL, 0.0, SCIPinfinity(scip),
                  initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
            break;
         case 'E' :
            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, mpsinputField2(mpsi), 0, NULL, NULL, 0.0, 0.0,
                  initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
            break;
         case 'L' :
            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, mpsinputField2(mpsi), 0, NULL, NULL, -SCIPinfinity(scip), 0.0,
                  initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
            break;
         default :
            mpsinputSyntaxerror(mpsi);
            return SCIP_OKAY;
         }
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }
   mpsinputSyntaxerror(mpsi);

   return SCIP_OKAY;
}

/* Process COLUMNS section. */
static
SCIP_RETCODE readCols(
   MPSINPUT*             mpsi,
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   char          colname[MPS_MAX_NAMELEN] = { '\0' };
   SCIP_CONS*    cons;
   SCIP_VAR*     var;
   SCIP_Real     val;

   var = NULL;
   while(mpsinputReadLine(mpsi))
   {
      if (mpsinputField0(mpsi) != 0)
      {
         if (strcmp(mpsinputField0(mpsi), "RHS"))
            break;

         /* add the last variable to the problem */
         if( var != NULL )
         {
            SCIP_CALL( SCIPaddVar(scip, var) );
            SCIP_CALL( SCIPreleaseVar(scip, &var) );
         }
         assert(var == NULL);

         mpsinputSetSection(mpsi, MPS_RHS);
         return SCIP_OKAY;
      }
      if (mpsinputField1(mpsi) == NULL || mpsinputField2(mpsi) == NULL || mpsinputField3(mpsi) == NULL)
         break;

      /* new column? */
      if (strcmp(colname, mpsinputField1(mpsi)))
      {
         SCIP_Bool dynamiccols;

         /* add the last variable to the problem */
         if( var != NULL )
         {
            SCIP_CALL( SCIPaddVar(scip, var) );
            SCIP_CALL( SCIPreleaseVar(scip, &var) );
         }
         assert(var == NULL);

         strncpy(colname, mpsinputField1(mpsi), MPS_MAX_NAMELEN - 1);

         SCIP_CALL( SCIPgetBoolParam(scip, "reading/mpsreader/dynamiccols", &dynamiccols) );

         if( mpsinputIsInteger(mpsi) )
         {
            /* for integer variables, default bounds are 0 <= x <= 1, and default cost is 0 */
            SCIP_CALL( SCIPcreateVar(scip, &var, colname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY, !dynamiccols, dynamiccols,
                  NULL, NULL, NULL, NULL) );
         }
         else
         {
            /* for continuous variables, default bounds are 0 <= x, and default cost is 0 */
            SCIP_CALL( SCIPcreateVar(scip, &var, colname, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS,
                  !dynamiccols, dynamiccols, NULL, NULL, NULL, NULL) );
         }
      }
      assert(var != NULL);

      val = atof(mpsinputField3(mpsi));

      if (!strcmp(mpsinputField2(mpsi), mpsinputObjname(mpsi)))
      {
         SCIP_CALL( SCIPchgVarObj(scip, var, val) );
      }
      else
      {
         cons = SCIPfindCons(scip, mpsinputField2(mpsi));
         if( cons == NULL )
            mpsinputEntryIgnored(scip, mpsi, "Column", mpsinputField1(mpsi), "row", mpsinputField2(mpsi), SCIP_VERBLEVEL_FULL);
         else if( !SCIPisZero(scip, val) )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, var, val) );
         }
      }
      if (mpsinputField5(mpsi) != NULL)
      {
         assert(mpsinputField4(mpsi) != NULL);

         val = atof(mpsinputField5(mpsi));

         if (!strcmp(mpsinputField4(mpsi), mpsinputObjname(mpsi)))
         {
            SCIP_CALL( SCIPchgVarObj(scip, var, val) );
         }
         else
         {
            cons = SCIPfindCons(scip, mpsinputField4(mpsi));
            if( cons == NULL )
               mpsinputEntryIgnored(scip, mpsi, "Column", mpsinputField1(mpsi), "row", mpsinputField4(mpsi), SCIP_VERBLEVEL_FULL);
            else if( !SCIPisZero(scip, val) )
            {
               SCIP_CALL( SCIPaddCoefLinear(scip, cons, var, val) );
            }
         }
      }
   }
   mpsinputSyntaxerror(mpsi);

   return SCIP_OKAY;
}

/* Process RHS section. */
static
SCIP_RETCODE readRhs(
   MPSINPUT*             mpsi,
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   char        rhsname[MPS_MAX_NAMELEN] = { '\0' };
   SCIP_CONS*  cons;
   SCIP_Real   lhs;
   SCIP_Real   rhs;
   SCIP_Real   val;

   while(mpsinputReadLine(mpsi))
   {
      if (mpsinputField0(mpsi) != NULL)
      {
         if (!strcmp(mpsinputField0(mpsi), "RANGES"))
            mpsinputSetSection(mpsi, MPS_RANGES);
         else if (!strcmp(mpsinputField0(mpsi), "BOUNDS"))
            mpsinputSetSection(mpsi, MPS_BOUNDS);
         else if (!strcmp(mpsinputField0(mpsi), "SOS"))
            mpsinputSetSection(mpsi, MPS_SOS);
         else if (!strcmp(mpsinputField0(mpsi), "ENDATA"))
            mpsinputSetSection(mpsi, MPS_ENDATA);
         else
            break;
         return SCIP_OKAY;
      }
      if ((mpsinputField2(mpsi) != NULL && mpsinputField3(mpsi) == NULL)
         || (mpsinputField4(mpsi) != NULL && mpsinputField5(mpsi) == NULL))
         mpsinputInsertName(mpsi, "_RHS_", FALSE);

      if (mpsinputField1(mpsi) == NULL || mpsinputField2(mpsi) == NULL || mpsinputField3(mpsi) == NULL)
         break;

      if (*rhsname == '\0')
         strncpy(rhsname, mpsinputField1(mpsi), MPS_MAX_NAMELEN - 1);

      if (!strcmp(rhsname, mpsinputField1(mpsi)))
      {
         cons = SCIPfindCons(scip, mpsinputField2(mpsi));
         if( cons == NULL )
            mpsinputEntryIgnored(scip, mpsi, "RHS", mpsinputField1(mpsi), "row", mpsinputField2(mpsi), SCIP_VERBLEVEL_NORMAL);
         else
         {
            val = atof(mpsinputField3(mpsi));

            /* find out the row sense */
            lhs = SCIPgetLhsLinear(scip, cons);
            rhs = SCIPgetRhsLinear(scip, cons);
            if( SCIPisInfinity(scip, -lhs) )
            {
               /* lhs = -infinity -> lower or equal */
               assert(SCIPisZero(scip, rhs));
               SCIP_CALL( SCIPchgRhsLinear(scip, cons, val) );
            }
            else if( SCIPisInfinity(scip, rhs) )
            {
               /* rhs = +infinity -> greater or equal */
               assert(SCIPisZero(scip, lhs));
               SCIP_CALL( SCIPchgLhsLinear(scip, cons, val) );
            }
            else
            {
               /* lhs > -infinity, rhs < infinity -> equality */
               assert(SCIPisZero(scip, lhs));
               assert(SCIPisZero(scip, rhs));
               SCIP_CALL( SCIPchgLhsLinear(scip, cons, val) );
               SCIP_CALL( SCIPchgRhsLinear(scip, cons, val) );
            }
	    SCIPdebugMessage("RHS <%s> lhs: %g  rhs: %g  val: <%22.12g>\n", mpsinputField2(mpsi), lhs, rhs, val);
         }
         if (mpsinputField5(mpsi) != NULL)
         {
            cons = SCIPfindCons(scip, mpsinputField4(mpsi));
            if( cons == NULL )
               mpsinputEntryIgnored(scip, mpsi, "RHS", mpsinputField1(mpsi), "row", mpsinputField4(mpsi), SCIP_VERBLEVEL_NORMAL);
            else
            {
               val = atof(mpsinputField5(mpsi));

               /* find out the row sense */
               lhs = SCIPgetLhsLinear(scip, cons);
               rhs = SCIPgetRhsLinear(scip, cons);
               if( SCIPisInfinity(scip, -lhs) )
               {
                  /* lhs = -infinity -> lower or equal */
                  assert(SCIPisZero(scip, rhs));
                  SCIP_CALL( SCIPchgRhsLinear(scip, cons, val) );
               }
               else if( SCIPisInfinity(scip, rhs) )
               {
                  /* rhs = +infinity -> greater or equal */
                  assert(SCIPisZero(scip, lhs));
                  SCIP_CALL( SCIPchgLhsLinear(scip, cons, val) );
               }
               else
               {
                  /* lhs > -infinity, rhs < infinity -> equality */
                  assert(SCIPisZero(scip, lhs));
                  assert(SCIPisZero(scip, rhs));
                  SCIP_CALL( SCIPchgLhsLinear(scip, cons, val) );
                  SCIP_CALL( SCIPchgRhsLinear(scip, cons, val) );
               }
	       SCIPdebugMessage("RHS <%s> lhs: %g  rhs: %g  val: <%22.12g>\n", mpsinputField4(mpsi), lhs, rhs, val);
            }
         }
      }
   }
   mpsinputSyntaxerror(mpsi);

   return SCIP_OKAY;
}

/* Process RANGES section */
static
SCIP_RETCODE readRanges(
   MPSINPUT*             mpsi,
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   char        rngname[MPS_MAX_NAMELEN] = { '\0' };
   SCIP_CONS*  cons;
   SCIP_Real   lhs;
   SCIP_Real   rhs;
   SCIP_Real   val;

   while(mpsinputReadLine(mpsi))
   {
      if (mpsinputField0(mpsi) != NULL)
      {
         /*printf("Range name     : %s\n", rngname);*/

         if (!strcmp(mpsinputField0(mpsi), "BOUNDS"))
            mpsinputSetSection(mpsi, MPS_BOUNDS);
         else if (!strcmp(mpsinputField0(mpsi), "SOS"))
            mpsinputSetSection(mpsi, MPS_SOS);
         else if (!strcmp(mpsinputField0(mpsi), "ENDATA"))
            mpsinputSetSection(mpsi, MPS_ENDATA);
         else
            break;
         return SCIP_OKAY;
      }
      if ((mpsinputField2(mpsi) != NULL && mpsinputField3(mpsi) == NULL)
         || (mpsinputField4(mpsi) != NULL && mpsinputField5(mpsi) == NULL))
         mpsinputInsertName(mpsi, "_RNG_", FALSE);

      if (mpsinputField1(mpsi) == NULL || mpsinputField2(mpsi) == NULL || mpsinputField3(mpsi) == NULL)
         break;

      if (*rngname == '\0')
         strncpy(rngname, mpsinputField1(mpsi), MPS_MAX_NAMELEN - 1);

      /* The rules are:
       * Row Sign   LHS             RHS
       * ----------------------------------------
       *  G   +/-   rhs             rhs + |range|
       *  L   +/-   rhs - |range|   rhs
       *  E   +     rhs             rhs + range
       *  E   -     rhs + range     rhs
       * ----------------------------------------
       */
      if (!strcmp(rngname, mpsinputField1(mpsi)))
      {
         cons = SCIPfindCons(scip, mpsinputField2(mpsi));
         if( cons == NULL )
            mpsinputEntryIgnored(scip, mpsi, "Range", mpsinputField1(mpsi), "row", mpsinputField2(mpsi), SCIP_VERBLEVEL_NORMAL);
         else
         {
            val = atof(mpsinputField3(mpsi));

            /* find out the row sense */
            lhs = SCIPgetLhsLinear(scip, cons);
            rhs = SCIPgetRhsLinear(scip, cons);
            if( SCIPisInfinity(scip, -lhs) )
            {
               /* lhs = -infinity -> lower or equal */
               SCIP_CALL( SCIPchgLhsLinear(scip, cons, rhs - REALABS(val)) );
            }
            else if( SCIPisInfinity(scip, rhs) )
            {
               /* rhs = +infinity -> greater or equal */
               SCIP_CALL( SCIPchgRhsLinear(scip, cons, lhs + REALABS(val)) );
            }
            else
            {
               /* lhs > -infinity, rhs < infinity -> equality */
               assert(SCIPisEQ(scip, lhs, rhs));
               if( val >= 0.0 )
               {
                  SCIP_CALL( SCIPchgRhsLinear(scip, cons, rhs + val) );
               }
               else
               {
                  SCIP_CALL( SCIPchgLhsLinear(scip, cons, lhs + val) );
               }
            }
         }
         if (mpsinputField5(mpsi) != NULL)
         {
            cons = SCIPfindCons(scip, mpsinputField4(mpsi));
            if( cons == NULL )
               mpsinputEntryIgnored(scip, mpsi, "Range", mpsinputField1(mpsi), "row", mpsinputField4(mpsi), SCIP_VERBLEVEL_NORMAL);
            else
            {
               val = atof(mpsinputField5(mpsi));

               /* find out the row sense */
               lhs = SCIPgetLhsLinear(scip, cons);
               rhs = SCIPgetRhsLinear(scip, cons);
               if( SCIPisInfinity(scip, -lhs) )
               {
                  /* lhs = -infinity -> lower or equal */
                  SCIP_CALL( SCIPchgLhsLinear(scip, cons, rhs - REALABS(val)) );
               }
               else if( SCIPisInfinity(scip, rhs) )
               {
                  /* rhs = +infinity -> greater or equal */
                  SCIP_CALL( SCIPchgRhsLinear(scip, cons, lhs + REALABS(val)) );
               }
               else
               {
                  /* lhs > -infinity, rhs < infinity -> equality */
                  assert(SCIPisEQ(scip, lhs, rhs));
                  if( val >= 0.0 )
                  {
                     SCIP_CALL( SCIPchgRhsLinear(scip, cons, rhs + val) );
                  }
                  else
                  {
                     SCIP_CALL( SCIPchgLhsLinear(scip, cons, lhs + val) );
                  }
               }
            }
         }
      }
   }
   mpsinputSyntaxerror(mpsi);

   return SCIP_OKAY;
}

/* Process BOUNDS section. */
static
SCIP_RETCODE readBounds(
   MPSINPUT*             mpsi,
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   char        bndname[MPS_MAX_NAMELEN] = { '\0' };
   SCIP_VAR*   var;
   SCIP_Real   val;

   while(mpsinputReadLine(mpsi))
   {
      if (mpsinputField0(mpsi) != 0)
      {
         /*printf("Bound name     : %s\n", bndname);*/

         if (!strcmp(mpsinputField0(mpsi), "SOS"))
            mpsinputSetSection(mpsi, MPS_SOS);
         else if (!strcmp(mpsinputField0(mpsi), "ENDATA"))
            mpsinputSetSection(mpsi, MPS_ENDATA);
	 else
	    break;
         return SCIP_OKAY;
      }
      /* Is the value field used ? */
      if (  !strcmp(mpsinputField1(mpsi), "LO")  /* lower bound given in field 4 */
         || !strcmp(mpsinputField1(mpsi), "UP")  /* upper bound given in field 4 */
         || !strcmp(mpsinputField1(mpsi), "FX")  /* fixed value given in field 4 */
         || !strcmp(mpsinputField1(mpsi), "LI")  /* CPLEX extension: lower bound of integer variable given in field 4 */
         || !strcmp(mpsinputField1(mpsi), "UI")  /* CPLEX extension: upper bound of integer variable given in field 4 */
         || !strcmp(mpsinputField1(mpsi), "SC")) /* CPLEX extension: semi continuous variable */
      {
         if (mpsinputField3(mpsi) != NULL && mpsinputField4(mpsi) == NULL)
            mpsinputInsertName(mpsi, "_BND_", TRUE);
         if( !mpsi->semicontwarning && !strcmp(mpsinputField1(mpsi), "SC") )
         {
            mpsinputEntryIgnored(scip, mpsi, "not supported semi continuous declaration", mpsinputField1(mpsi),
               "variable", mpsinputField3(mpsi), SCIP_VERBLEVEL_NORMAL);
            mpsi->semicontwarning = TRUE;
         }
      }
      else if( !strcmp(mpsinputField1(mpsi), "FR") /* free variable */
         || !strcmp(mpsinputField1(mpsi), "MI")    /* lower bound is minus infinity */
         || !strcmp(mpsinputField1(mpsi), "PL")    /* upper bound is plus infinity */
         || !strcmp(mpsinputField1(mpsi), "BV"))   /* CPLEX extension: binary variable */
      {
         if (mpsinputField2(mpsi) != NULL && mpsinputField3(mpsi) == NULL)
            mpsinputInsertName(mpsi, "_BND_", TRUE);
      }
      else
      {
         mpsinputSyntaxerror(mpsi);
         return SCIP_OKAY;
      }

      if (mpsinputField1(mpsi) == NULL || mpsinputField2(mpsi) == NULL || mpsinputField3(mpsi) == NULL)
         break;

      if (*bndname == '\0')
         strncpy(bndname, mpsinputField2(mpsi), MPS_MAX_NAMELEN - 1);

      /* Only read the first Bound in section */
      if (!strcmp(bndname, mpsinputField2(mpsi)))
      {
         var = SCIPfindVar(scip, mpsinputField3(mpsi));
         if( var == NULL )
            mpsinputEntryIgnored(scip, mpsi, "column", mpsinputField3(mpsi), "bound", bndname, SCIP_VERBLEVEL_NORMAL);
         else
         {
            if( mpsinputField4(mpsi) == NULL )
               val = 0.0;
            else
               val = atof(mpsinputField4(mpsi));

            /* if a bound of a binary variable is given, the variable is converted into an integer variable
             * with default bounds 0 <= x <= infinity
             */
            if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
            {
               if( (mpsinputField1(mpsi)[1] == 'I') /* ILOG extension (Integer Bound) */
                  || (!(mpsinputField1(mpsi)[0] == 'L' && SCIPisEQ(scip, val, 0.0))
                     && !(mpsinputField1(mpsi)[0] == 'U' && SCIPisEQ(scip, val, 1.0))) )
               {
                  assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 0.0));
                  assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(var), 1.0));
                  SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_INTEGER) );
                  SCIP_CALL( SCIPchgVarUb(scip, var, SCIPinfinity(scip)) );
               }
            }

            switch(mpsinputField1(mpsi)[0])
            {
            case 'L':
               if( mpsinputField1(mpsi)[1] == 'I' ) /* ILOG extension (Integer Bound) */
               {
                  SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_INTEGER) );
               }
               SCIP_CALL( SCIPchgVarLb(scip, var, val) );
               break;
            case 'U':
               if( mpsinputField1(mpsi)[1] == 'I' ) /* ILOG extension (Integer Bound) */
               {
                  SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_INTEGER) );
               }
               SCIP_CALL( SCIPchgVarUb(scip, var, val) );
               break;
            case 'F':
               if (mpsinputField1(mpsi)[1] == 'X')
               {
                  SCIP_CALL( SCIPchgVarLb(scip, var, val) );
                  SCIP_CALL( SCIPchgVarUb(scip, var, val) );
               }
               else
               {
                  SCIP_CALL( SCIPchgVarLb(scip, var, -SCIPinfinity(scip)) );
                  SCIP_CALL( SCIPchgVarUb(scip, var, +SCIPinfinity(scip)) );
               }
               break;
            case 'M':
               SCIP_CALL( SCIPchgVarLb(scip, var, -SCIPinfinity(scip)) );
               break;
            case 'P':
               SCIP_CALL( SCIPchgVarUb(scip, var, +SCIPinfinity(scip)) );
               break;
            case 'B' : /* Ilog extension (Binary) */
               SCIP_CALL( SCIPchgVarLb(scip, var, 0.0) );
               SCIP_CALL( SCIPchgVarUb(scip, var, 1.0) );
               SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY) );
               break;
            default:
               mpsinputSyntaxerror(mpsi);
               return SCIP_OKAY;
            }
         }
      }
   }
   mpsinputSyntaxerror(mpsi);

   return SCIP_OKAY;
}


/* Process SOS section.
 *
 * We read the SOS section, which is a nonstandard section introduced by CPLEX.
 *
 * @note Currently do not support the standard way of specifying SOS constraints via markers.
 */
static
SCIP_RETCODE readSOS(
   MPSINPUT*             mpsi,
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Bool initial, separate, enforce, check, propagate;
   SCIP_Bool local, modifiable, dynamic, removable;
   char name[MPS_MAX_NAMELEN] = { '\0' };
   SCIP_CONS*  cons = NULL;
   int consType = -1;
   int cnt = 0;

   /* standard settings for SOS constraints: */
   initial = TRUE;
   separate = FALSE;
   enforce = TRUE;
   check = TRUE;
   propagate = TRUE;
   local = FALSE;
   modifiable = FALSE;
   dynamic = FALSE;
   removable = FALSE;

   /* loop through section */
   while ( mpsinputReadLine(mpsi) )
   {
      int type = -1;

      /* check if next section is found */
      if ( mpsinputField0(mpsi) != NULL )
      {
	 if ( ! strcmp(mpsinputField0(mpsi), "ENDATA") )
            mpsinputSetSection(mpsi, MPS_ENDATA);
	 break;
      }
      if ( mpsinputField1(mpsi) == NULL && mpsinputField2(mpsi) == NULL )
      {
	 SCIPerrorMessage("empty data in a non-comment line.\n");
	 mpsinputSyntaxerror(mpsi);
	 return SCIP_OKAY;
      }

      /* check for new SOS set */
      if ( strcmp(mpsinputField1(mpsi), "S1") == 0 )
	 type = 1;
      if ( strcmp(mpsinputField1(mpsi), "S2") == 0 )
	 type = 2;

      /* add last constraint and create a new one */
      if ( type > 0 )
      {
	 assert( type == 1 || type == 2 );
	 if ( cons != NULL )
	 {
	    /* add last constraint */
	    SCIP_CALL( SCIPaddCons(scip, cons) );
	    SCIPdebugMessage("(line %d) added constraint <%s>: ", mpsi->lineno, SCIPconsGetName(cons));
	    SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
	    SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	 }

	 /* check name */
	 if ( mpsinputField2(mpsi) != NULL )
	    strncpy(name, mpsinputField2(mpsi), MPS_MAX_NAMELEN - 1);
	 else
	 {
	    /* create new name */
	    snprintf(name, MPS_MAX_NAMELEN, "SOS%d", ++cnt);
	 }

	 /* create new SOS constraint */
	 if ( type == 1 )
	 {
	    /* we do not know the name of the constraint */
	    SCIP_CALL( SCIPcreateConsSOS1(scip, &cons, name, 0, NULL, NULL, initial, separate, enforce, check, propagate,
                  local, modifiable, dynamic, removable) );
	 }
	 else
	 {
	    assert( type == 2 );
	    SCIP_CALL( SCIPcreateConsSOS2(scip, &cons, name, 0, NULL, NULL, initial, separate, enforce, check, propagate,
                  local, modifiable, dynamic, removable) );
	 }
	 consType = type;
	 SCIPdebugMessage("created constraint <%s> of type %d.\n", name, type);
	 /* note: we ignore the priorities! */
      }
      else
      {
	 /* otherwise we are in the section given variables */
	 SCIP_VAR* var = NULL;
	 SCIP_Real weight = 0.0;
	 char* endptr;

	 if ( consType != 1 && consType != 2 )
	 {
	    SCIPerrorMessage("missing SOS type specification.\n");
	    mpsinputSyntaxerror(mpsi);
	    return SCIP_OKAY;
	 }

	 /* get variable */
	 var = SCIPfindVar(scip, mpsinputField1(mpsi));
	 if ( var == NULL )
	 {
	    /* ignore unkown variables - we would not know the type anyway */
            mpsinputEntryIgnored(scip, mpsi, "column", mpsinputField1(mpsi), "SOS", name, SCIP_VERBLEVEL_NORMAL);
	 }
	 else
	 {
	    /* get weight */
	    weight = strtod(mpsinputField2(mpsi), &endptr);
	    if ( endptr == mpsinputField2(mpsi) || *endptr != '\0' )
	    {
	       SCIPerrorMessage("weight for variable <%s> not specified.\n", mpsinputField1(mpsi));
	       mpsinputSyntaxerror(mpsi);
	       return SCIP_OKAY;
	    }

	    /* add variable and weight */
	    assert( consType == 1 || consType == 2 );
	    switch (consType)
	    {
	    case 1: SCIP_CALL( SCIPaddVarSOS1(scip, cons, var, weight) ); break;
	    case 2: SCIP_CALL( SCIPaddVarSOS2(scip, cons, var, weight) ); break;
	    default: abort(); // should not happen
	    }
	    SCIPdebugMessage("added variable <%s> with weight %g.\n", SCIPvarGetName(var), weight);
	 }
	 /* check other fields */
	 if ( ( mpsinputField3(mpsi) != NULL && *mpsinputField3(mpsi) != '\0' ) ||
            ( mpsinputField4(mpsi) != NULL && *mpsinputField4(mpsi) != '\0' ) ||
            ( mpsinputField5(mpsi) != NULL && *mpsinputField5(mpsi) != '\0' ) )
	 {
	    SCIPwarningMessage("ignoring data in fields 3-5 <%s> <%s> <%s>.\n",
               mpsinputField3(mpsi), mpsinputField4(mpsi), mpsinputField5(mpsi));
	 }
      }
   }

   if ( cons != NULL )
   {
      /* add last constraint */
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIPdebugMessage("(line %d) added constraint <%s>: ", mpsi->lineno, SCIPconsGetName(cons));
      SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   return SCIP_OKAY;
}



/* Read LP in "MPS File Format".
 *
 *  A specification of the MPS format can be found at
 *
 *  http://plato.asu.edu/ftp/mps_format.txt,
 *  ftp://ftp.caam.rice.edu/pub/people/bixby/miplib/mps_format,
 *
 *  and in the
 *
 *  ILOG CPLEX Reference Manual
 *
 *  This routine should read all valid MPS format files.
 *  What it will not do, is to find all cases where a file is ill formed.
 *  If this happens it may complain and read nothing or read "something".
 */
static
SCIP_RETCODE readMps(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of the input file */
   )
{
   SCIP_FILE* fp;
   MPSINPUT* mpsi;
   SCIP_Bool error;

   assert(scip != NULL);
   assert(filename != NULL);

   fp = SCIPfopen(filename, "r");
   if (fp == NULL)
   {
      char buf[1024];
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      (void) strerror_r(errno, buf, 1024);
      SCIPerrorMessage("%s: %s\n", filename, buf);
      return SCIP_NOFILE;
   }

   SCIP_CALL( mpsinputCreate(scip, &mpsi, fp) );

   SCIP_CALL( readName(mpsi) );

   SCIP_CALL( SCIPcreateProb(scip, mpsi->probname, NULL, NULL, NULL, NULL, NULL, NULL) );

   if (mpsinputSection(mpsi) == MPS_OBJSEN)
   {
      SCIP_CALL( readObjsen(mpsi) );
   }
   if (mpsinputSection(mpsi) == MPS_OBJNAME)
   {
      SCIP_CALL( readObjname(mpsi) );
   }
   while (mpsinputSection(mpsi) == MPS_ROWS
      || mpsinputSection(mpsi) == MPS_USERCUTS
      || mpsinputSection(mpsi) == MPS_LAZYCONS)
   {
      SCIP_CALL( readRows(mpsi, scip) );
   }
   if (mpsinputSection(mpsi) == MPS_COLUMNS)
   {
      SCIP_CALL( readCols(mpsi, scip) );
   }
   if (mpsinputSection(mpsi) == MPS_RHS)
   {
      SCIP_CALL( readRhs(mpsi, scip) );
   }
   if (mpsinputSection(mpsi) == MPS_RANGES)
   {
      SCIP_CALL( readRanges(mpsi, scip) );
   }
   if (mpsinputSection(mpsi) == MPS_BOUNDS)
   {
      SCIP_CALL( readBounds(mpsi, scip) );
   }
   if (mpsinputSection(mpsi) == MPS_SOS)
   {
      SCIP_CALL( readSOS(mpsi, scip) );
   }
   if (mpsinputSection(mpsi) != MPS_ENDATA)
      mpsinputSyntaxerror(mpsi);

   SCIPfclose(fp);

   error = mpsinputHasError(mpsi);

   if( !error )
   {
      SCIP_CALL( SCIPsetObjsense(scip, mpsinputObjsense(mpsi)) );

      /*printf("Objective sense: %s\n", (mpsinputObjsense(mpsi) == SCIP_OBJSENSE_MINIMIZE) ? "Minimize" : "Maximize");*/
   }
   mpsinputFree(scip, &mpsi);

   if( error )
      return SCIP_PARSEERROR;
   else
      return SCIP_OKAY;
}

/*
 * local methods for writing problem
 */

/** hash key retrieval function for variables */
static
SCIP_DECL_HASHGETKEY(hashGetKeyVar)
{
   return elem;
}

/** returns TRUE iff the indices of both variables are equal */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqVar)
{
   if ( key1 == key2 )
      return TRUE;
   return FALSE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValVar)
{
   assert( SCIPvarGetIndex((SCIP_VAR*) key) >= 0 );
   return (unsigned int) SCIPvarGetIndex((SCIP_VAR*) key);
}


/* computes the field width such that the output file is nicely arranged */
static
unsigned int computeFieldWidth(
   unsigned int             width              /**< required width */
   )
{
   width = MAX(8, width);
   return MIN(MPS_MAX_FIELDLEN, width);
}


/* output two strings in columns 1 and 2 with computed widths */
static
void printRecord(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           col1,               /**< column 1 */
   const char*           col2,               /**< column 2 */
   unsigned int          maxnamelen          /**< maximum name length */
   )
{
   unsigned int fieldwidth;
   char format[32];

   assert( scip != NULL );
   assert( col1 != NULL );
   assert( col2 != NULL );
   assert( strlen(col1) < MPS_MAX_NAMELEN );
   assert( strlen(col2) < MPS_MAX_VALUELEN );
   assert( maxnamelen > 0 );

   fieldwidth = computeFieldWidth(maxnamelen);
   snprintf(format, 32," %%-%ds %%%ds ", fieldwidth, MPS_MAX_VALUELEN - 1);

   SCIPinfoMessage(scip, file, format, col1, col2);
}

/* output two strings in columns 1 (width 2) and 2 (width 8) */
static
void printStart(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           col1,               /**< column 1 */
   const char*           col2,               /**< column 2 */
   int                   maxnamelen          /**< maximum name length (-1 if irrelevant) */
   )
{
   unsigned int fieldwidth;
   char format[32];

   assert( scip != NULL );
   assert( col1 != NULL );
   assert( col2 != NULL );
   assert( strlen(col1) <= 2 );
   assert( strlen(col2) < MPS_MAX_NAMELEN );
   assert( maxnamelen == -1 || maxnamelen > 0 );

   if( maxnamelen < 0 )
   {
      /* format does not matter */
      snprintf(format, 32, " %%-2.2s %%-s ");
   }
   else
   {
      fieldwidth = computeFieldWidth((unsigned int) maxnamelen);
      snprintf(format, 32, " %%-2.2s %%-%ds ", fieldwidth);
   }

   SCIPinfoMessage(scip, file, format, col1, col2);
}

/* prints the given data as column entry */
static
void printEntry(
   SCIP*              scip,               /**< SCIP data structure */
   FILE*              file,               /**< output file (or NULL for standard output) */
   const char*        varname,            /**< variable name */
   const char*        consname,           /**< constraint name */
   SCIP_Real          value,              /**< value to display */
   int*               recordcnt,          /**< pointer to store the number of records per line */
   unsigned int       maxnamelen          /**< maximum name length */
   )
{
   char valuestr[MPS_MAX_VALUELEN] = { '\0' };

   assert( scip != NULL );
   assert( recordcnt != NULL );
   assert( *recordcnt >= 0 && *recordcnt < 2 );

   if( !SCIPisZero(scip, value) )
   {
      snprintf(valuestr, MPS_MAX_VALUELEN, "%25.15g", value);

      if ( *recordcnt == 0 )
      {
         /* start new line with an empty first column and the variable name in the second column */
         printStart(scip, file, "", varname, (int) maxnamelen);
         *recordcnt = 0;
      }

      printRecord(scip, file, consname, valuestr, maxnamelen);
      (*recordcnt)++;

      if ( *recordcnt == 2 )
      {
         /* each line can have at most two records */
         SCIPinfoMessage(scip, file, "\n");
         *recordcnt = 0;
      }
   }
}

/* prints the constraint type to file stream */
static
void printRowType(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   const char*           name                /**< constraint name */
   )
{
   char rowtype[2];
   assert( scip != NULL );

   assert( !SCIPisInfinity(scip, -lhs) || !SCIPisInfinity(scip, rhs) );
   assert( SCIPisGT(scip, rhs, lhs) || SCIPisEQ(scip, lhs, rhs) );
   assert( name != NULL );


   if( SCIPisEQ(scip, lhs, rhs) )
      snprintf(rowtype, 2, "%s", "E");
   else
   {
      /* in case the right hand side and the left hand side are not infinity we print a
       * less or equal constraint and put the right hand side in the RHS section and the
       * left hand side (hidden) in the RANGE section */
      if( !SCIPisInfinity(scip, rhs) )
         snprintf(rowtype, 2, "%s", "L");
      else
      {
         assert( !SCIPisInfinity(scip, -lhs) );
         snprintf(rowtype, 2, "%s", "G");
      }
   }

   printStart(scip, file, rowtype, name, -1);
   SCIPinfoMessage(scip, file, "\n");
}


/** initializes the sparse matrix */
static
SCIP_RETCODE initializeMatrix(
   SCIP*              scip,               /**< SCIP data structure */
   SPARSEMATRIX**     matrix,             /**< pointer to sparse matrix containing the entries */
   int                slots               /**< number of slots */
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, matrix) );
   (*matrix)->nentries = 0;
   (*matrix)->sentries = slots;
   SCIP_CALL( SCIPallocBufferArray(scip, &(*matrix)->values, (*matrix)->sentries) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*matrix)->columns, (*matrix)->sentries) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*matrix)->rows, (*matrix)->sentries) );

   return SCIP_OKAY;
}

/** this method takes care that the required capacity is available in the sparse matrix */
static
SCIP_RETCODE checkSparseMatrixCapacity(
   SCIP*                 scip,               /**< SCIP data structure */
   SPARSEMATRIX*         matrix,             /**< sparse matrix for storing the coefficient */
   int                   capacity            /**< needed capacity */
   )
{
   if( matrix->nentries + capacity >= matrix->sentries )
   {
      matrix->sentries = matrix->sentries * 2 + capacity;
      SCIP_CALL( SCIPreallocBufferArray(scip, &matrix->values, matrix->sentries) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &matrix->columns, matrix->sentries) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &matrix->rows, matrix->sentries) );
   }
   return SCIP_OKAY;
}

/** frees the sparse matrix */
static
void freeMatrix(
   SCIP*              scip,               /**< SCIP data structure */
   SPARSEMATRIX*      matrix              /**< sparse matrix to free */
   )
{
   SCIPfreeBufferArray(scip, &matrix->columns);
   SCIPfreeBufferArray(scip, &matrix->rows);
   SCIPfreeBufferArray(scip, &matrix->values);
   
   SCIPfreeBuffer(scip, &matrix);
}



/** transforms given variables, scalars, and constant to the corresponding active variables, scalars, and constant */
static
SCIP_RETCODE getActiveVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< vars array to get active variables for */
   SCIP_Real*            scalars,            /**< scalars a_1, ..., a_n in linear sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c  */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int requiredsize;
   int v;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( scalars != NULL );
   assert( nvars != NULL );
   assert( constant != NULL );

   if( transformed )
   {
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, *nvars, constant, &requiredsize) );

      if( requiredsize > *nvars )
      {
         *nvars = requiredsize;
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, *nvars ) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &scalars, *nvars ) );

         SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, *nvars, constant, &requiredsize) );
         assert( requiredsize <= *nvars );
      }
   }
   else
   {
      for( v = 0; v < *nvars; ++v )
         SCIP_CALL( SCIPvarGetOrigvarSum(&vars[v], &scalars[v], constant) );
   }
   return SCIP_OKAY;
}


/** computes the coefficient for the given variables and linear constraint information */
static
SCIP_RETCODE getLinearCoeffs(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           consname,           /**< name of the constraint */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of coefficients values (or NULL if all coefficient values are 1) */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             transformed,        /**< transformed constraint? */
   SPARSEMATRIX*         matrix,             /**< sparse matrix for storing the coefficient */
   SCIP_Real*            rhs                 /**< pointer to right hand side */
   )
{
   int v;
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   int nactivevars;
   SCIP_Real activeconstant = 0.0;

   assert( scip != NULL );
   assert( nvars == 0 || vars != NULL );
   assert( nvars > 0 || vars == NULL );
   assert( !SCIPisInfinity(scip, *rhs) );
   assert( matrix != NULL );
   
   /* if the variables array contains no variables, then return without
    * doing any thing; The MPS format and LP format do not forbid this
    * situation */
   if( nvars == 0 ) 
      return SCIP_OKAY;

   /* duplicate variable and value array */
   nactivevars = nvars;
   SCIP_CALL( SCIPduplicateBufferArray(scip, &activevars, vars, nactivevars ) );

   if( vals != NULL )
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activevals, vals, nactivevars ) );
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nactivevars) );

      for( v = 0; v < nactivevars; ++v )
         activevals[v] = 1.0;
   }

   /* retransform given variables to active variables */
   SCIP_CALL( getActiveVariables(scip, activevars, activevals, &nactivevars, &activeconstant, transformed) );

   /* copy the (matrix) row into the sparse matrix */
   SCIP_CALL( checkSparseMatrixCapacity(scip, matrix, nactivevars) );
   assert( matrix->nentries + nactivevars < matrix->sentries );
   
   for( v = 0; v < nactivevars; ++v )
   {
      matrix->values[matrix->nentries] = activevals[v];
      matrix->columns[matrix->nentries] = activevars[v];
      matrix->rows[matrix->nentries] = consname;
      matrix->nentries++;
   }
   
   /* adjust right hand side */
   (*rhs) -= activeconstant;

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);

   return SCIP_OKAY;
}


/** check whether given variables are aggregated and put them into an array without duplication */
static
SCIP_RETCODE collectAggregatedVars(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_VAR**         vars,               /**< variable array */
   int                nvars,              /**< number of mutable variables in the problem */
   SCIP_VAR***        aggvars,            /**< pointer to array storing the aggregated variables on output */
   int*               naggvars,           /**< pointer to number of aggregated variables on output */
   int*               saggvars,           /**< pointer to number of slots in aggvars array */
   SCIP_HASHTABLE*    varAggregated       /**< hashtable for checking duplicates */
   )
{
   int v;
   SCIP_VAR* var;
   SCIP_VARSTATUS status;

   assert( scip != NULL );
   assert( aggvars != NULL );
   assert( naggvars != NULL );
   assert( saggvars != NULL );
   assert( varAggregated != NULL );

   /* check variables */
   for (v = 0; v < nvars; ++v)
   {
      var = vars[v];
      status = SCIPvarGetStatus(var);

      /* collect aggregated variables in a list */
      if( status >= SCIP_VARSTATUS_AGGREGATED )
      {
         assert( status == SCIP_VARSTATUS_AGGREGATED || 
            status == SCIP_VARSTATUS_MULTAGGR ||
            status == SCIP_VARSTATUS_NEGATED );

	 if ( ! SCIPhashtableExists(varAggregated, (void*) var) )
	 {
            if( (*saggvars) <= (*naggvars) )
            {
               (*saggvars) *= 2;
               SCIP_CALL( SCIPreallocBufferArray(scip, aggvars, (*saggvars)) );
            }
            assert( (*saggvars) > (*naggvars) );

	    (*aggvars)[*naggvars] = var;
	    (*naggvars)++;
	    SCIP_CALL( SCIPhashtableInsert(varAggregated, (void*) var) );
	 }
      }
   }
   return SCIP_OKAY;
}


/** method check if the variable names are not longer than MPS_MAX_NAMELEN - 1*/
static
SCIP_RETCODE checkVarnames(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_VAR**         vars,               /**< array of variables */
   int                nvars,              /**< number of variables */
   unsigned int*      maxnamelen,         /**< pointer to store rhe maximum name lenght */
   const char***      varnames,           /**< pointer to array of varible names */
   SCIP_HASHMAP**     varnameHashmap      /**< pointer to hash map storing variable, variable name mapping */      
   )
{
   int v;
   int faulty;
   char* varname;
   SCIP_VAR* var;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( maxnamelen != NULL );

   faulty = 0;

   /* allocate memory */
   SCIP_CALL( SCIPhashmapCreate(varnameHashmap, SCIPblkmem(scip), nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, varnames, nvars) );

   /* check if the variable names are not to long */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert( var != NULL );
      
      if( strlen(SCIPvarGetName(var)) >= MPS_MAX_NAMELEN )
      {
         faulty++;
         (*maxnamelen) = MPS_MAX_NAMELEN - 1;
      }
      else
      {
	 size_t l = strlen(SCIPvarGetName(var));
         (*maxnamelen) = MAX(*maxnamelen, l);
      }
 
      SCIP_CALL( SCIPallocBufferArray(scip, &varname, (int) *maxnamelen + 1) );
      snprintf(varname, (*maxnamelen) + 1, "%s", SCIPvarGetName(var) );
      
      /* insert variable with variable name into hash map */
      assert( !SCIPhashmapExists(*varnameHashmap, var) );
      SCIP_CALL( SCIPhashmapInsert(*varnameHashmap, var, (void*) (size_t) varname) );

      (*varnames)[v] = varname;
   }
   
   if( faulty > 0 )
   {
      SCIPwarningMessage("there are %d variable names which have to be cut down to %d characters; LP might be corrupted\n", 
         faulty, MPS_MAX_NAMELEN - 1);
   }
   return SCIP_OKAY;
}

/** method check if the constraint names are not longer than MPS_MAX_NAMELEN - 1 */
static
SCIP_RETCODE checkConsnames(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONS**        origconss,          /**< array of all constraints */
   int                norigconss,         /**< number of all constraints */
   SCIP_CONS***       conss,              /**< array to store only relevant constraints */
   int*               nconss,             /**< number of relevant constraints */
   SCIP_Bool          transformed,        /**< TRUE iff problem is the transformed problem */
   unsigned int*      maxnamelen,         /**< pointer to store rhe maximum name lenght */
   const char***      consnames           /**< pointer to array of constraint names */
   )
{
   int c, i;
   int faulty;
   SCIP_CONS* cons;
   char* consname;

   assert( scip != NULL );
   assert( conss != NULL );
   assert( maxnamelen != NULL );

   faulty = 0;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, conss, norigconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, consnames, norigconss) );

   for( i = 0, c = 0; i < norigconss; ++i )
   {
      size_t l;

      cons = origconss[i];
      assert( cons != NULL );

      /* in case the transformed probelms is written only constraint are
       * posted which are enabled in the current node */
      if( transformed && !SCIPconsIsEnabled(cons) )
         continue;
      
      (*conss)[c] = cons;
      
      if( strlen(SCIPconsGetName(cons)) >= MPS_MAX_NAMELEN )
      {
         faulty++;
         (*maxnamelen) = MPS_MAX_NAMELEN - 1;
      }

      l = strlen(SCIPconsGetName(cons));
      (*maxnamelen) = MAX(*maxnamelen, l);
      
      SCIP_CALL( SCIPallocBufferArray(scip, &consname, (int) *maxnamelen + 1) );
      snprintf(consname, (*maxnamelen) + 1, "%s", SCIPconsGetName(cons) );
      
      (*consnames)[c] = consname;
      c++;
   }
   
   if( faulty > 0 )
   {
      SCIPwarningMessage("there are %d constraint names which have to be cut down to %d characters; LP might be corrupted\n", 
         faulty, MPS_MAX_NAMELEN - 1);
   }
   
   (*nconss) = c;
   
   return SCIP_OKAY;
}


/* outputs the COULMNS section of the MPS format */
static
void printColumnSection(
   SCIP*              scip,               /**< SCIP data structure */
   FILE*              file,               /**<  output file, or NULL if standard output should be used */
   SPARSEMATRIX*      matrix,             /**< sparse matrix containing the entries */
   SCIP_HASHMAP*      varnameHashmap,     /**< map from SCIP_VAR* to variable name */
   unsigned int       maxnamelen          /**< maximum name lenght */
   )
{
   SCIP_Bool intSection;
   SCIP_VAR* var;
   const char* varname;
   SCIP_Real value;
   int v;
   int recordcnt;

   /* sort sparse matrix w.r.t. the variable indices */
   SCIPsortPtrPtrReal((void**) matrix->columns, (void**) matrix->rows, matrix->values, 
         SCIPvarComp, matrix->nentries);

   /* print COLUMNS section */
   SCIPinfoMessage(scip, file, "COLUMNS\n");

   intSection = FALSE;
   
   for( v = 0; v < matrix->nentries; ++v )
   {
      var = matrix->columns[v];
      assert( var != NULL );
         
      if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS && intSection)
      {
         /* end integer section in MPS format */
         printStart(scip, file, "", "INTEND", (int) maxnamelen);
         printRecord(scip, file, "'MARKER'", "", maxnamelen);
         printRecord(scip, file, "'INTEND'", "", maxnamelen);
         SCIPinfoMessage(scip, file, "\n", maxnamelen);
         intSection = FALSE;
      }
      else if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS && !intSection )
      {
         /* start integer section in MPS format */
         printStart(scip, file, "", "INTSTART", (int) maxnamelen);
         printRecord(scip, file, "'MARKER'", "", maxnamelen);
         printRecord(scip, file, "'INTORG'", "", maxnamelen);
         SCIPinfoMessage(scip, file, "\n");
         intSection = TRUE;
      }

      SCIPdebugMessage("create entries for variable <%s>\n", SCIPvarGetName(var));
      
      /* record count; there are at most two records per line */
      recordcnt = 0;
       
      /* get variable name */
      assert ( SCIPhashmapExists(varnameHashmap, var) );
      varname = (const char*) (size_t) SCIPhashmapGetImage(varnameHashmap, var);

      /* output all entries of the same variable */
      do {
         value = matrix->values[v];
         
         /* print record to file */
         printEntry( scip, file, varname, matrix->rows[v], value, &recordcnt, maxnamelen );
         v++;
      }
      while (v < matrix->nentries && var == matrix->columns[v]);
         
      if( recordcnt == 1 )
         SCIPinfoMessage(scip, file, "\n");
         
      v--;
   }
}


/** outputs the right hand side section */
static
void printRhsSection(
   SCIP*              scip,               /**< SCIP data structure */
   FILE*              file,               /**<  output file, or NULL if standard output should be used */
   SCIP_CONS**        conss,              /**< constraint array */
   int                nconss,             /**< number of constraints */
   const char**       consnames,          /**< constraint names */
   SCIP_Real*         rhss,               /**< right hand side array */
   SCIP_Bool          transformed,        /**< TRUE iff problem is the transformed problem */
   unsigned int       maxnamelen          /**< maximum name lenght */
   )
{
   int c;
   int recordcnt;
   SCIP_CONS* cons;

   assert( rhss != NULL );

   SCIPinfoMessage(scip, file, "RHS\n");
   SCIPdebugMessage("start printing RHS section\n");
   
   recordcnt = 0;

   /* take care of the linear constraints */
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert( cons != NULL );

      /* skip all contraints which have a right hand side of infinity */
      if( SCIPisInfinity(scip, rhss[c]) )
         continue;
         
         /* in case the transformed problems is written only constraint are
       * posted which are enabled in the current node; the conss array
       * should only contain relevant constraints since these are collected
       * in the beginning in the methode checkConsnames()  */
      assert( !transformed || SCIPconsIsEnabled(cons) );
      
      assert( conss[c] != NULL );
      assert( consnames[c] != NULL );
      assert( !SCIPisInfinity(scip, rhss[c]) );
      
      printEntry( scip, file, "RHS", consnames[c], rhss[c], &recordcnt, maxnamelen );
   }

   if( recordcnt == 1 )
      SCIPinfoMessage(scip, file, "\n");
}


/** outputs the range section */
static
void printRangeSection(
   SCIP*              scip,               /**< SCIP data structure */
   FILE*              file,               /**<  output file, or NULL if standard output should be used */
   SCIP_CONS**        conss,              /**< constraint array */
   int                nconss,             /**< number of constraints */
   const char**       consnames,          /**< constraint names */
   SCIP_Bool          transformed,        /**< TRUE iff problem is the transformed problem */
   unsigned int       maxnamelen          /**< maximum name lenght */
   )
{   
   int c;
   int recordcnt;
   
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;

   SCIP_CONS* cons;
   SCIP_Real lhs;
   SCIP_Real rhs;
   

   SCIPinfoMessage(scip, file, "RANGES\n");
   SCIPdebugMessage("start printing RANGES section\n");

   recordcnt = 0;
   
   for( c = 0; c < nconss; ++c  )
   {
      cons = conss[c];
      assert( cons != NULL);

      /* in case the transformed problems is written only constraint are
       * posted which are enabled in the current node; the conss array
       * should only contain relevant constraints since these are collected
       * in the beginning in the methode checkConsnames()  */
      assert( !transformed || SCIPconsIsEnabled(cons) );

      assert( consnames[c] != NULL );

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );
      
      conshdlrname = SCIPconshdlrGetName(conshdlr);

      if( strcmp(conshdlrname, "linear") == 0 )
      {
         lhs = SCIPgetLhsLinear(scip, cons);
         rhs = SCIPgetRhsLinear(scip, cons);
      }
      else if( strcmp(conshdlrname, "varbound") == 0 )
      {
         lhs = SCIPgetLhsVarbound(scip, cons);
         rhs = SCIPgetRhsVarbound(scip, cons);
      }
      else
         continue;
      
      if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) && !SCIPisEQ(scip, rhs, lhs) )
      {
         assert( SCIPisGT(scip, rhs, lhs) );
         printEntry( scip, file, "RANGE", consnames[c], rhs - lhs, &recordcnt, maxnamelen );
      }
   }
   if(recordcnt == 1 )
      SCIPinfoMessage(scip, file, "\n");
}


/** output bound section */
static
void printBoundSection(
   SCIP*              scip,               /**< SCIP data structure */
   FILE*              file,               /**<  output file, or NULL if standard output should be used */
   SCIP_VAR**         vars,               /**< active variables */
   int                nvars,              /**< number of active variables */
   SCIP_VAR**         aggvars,            /**< needed aggregated variables */
   int                naggvars,           /**< number of aggregated variables */
   SCIP_Bool          transformed,        /**< TRUE iff problem is the transformed problem */
   const char**       varnames,           /**< array with variable names */
   unsigned int       maxnamelen          /**< maximum name lenght */
   )
{
   int v;
   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Real ub;
   const char* varname;
   char valuestr[MPS_MAX_VALUELEN] = { '\0' };

   assert( scip != NULL );
   assert( vars != NULL );
   
   SCIPinfoMessage(scip, file, "BOUNDS\n");
   SCIPdebugMessage("start printing BOUNDS section\n");

   /* output the active variables */
   for (v = 0; v < nvars; ++v)
   {
      var = vars[v];
      assert( var != NULL );

      /* get variable name */
      varname = varnames[v];
      assert( strncmp(varname, SCIPvarGetName(var), maxnamelen) == 0 );   
      
      if( transformed )
      {
         /* in case the transformed is written only local bounds are posted
          * which are valid in the current node */
         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);
      }
      else
      {
         lb = SCIPvarGetLbOriginal(var);
         ub = SCIPvarGetUbOriginal(var);
      }

      /* take care of binary variables */
      if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
      {
         printStart(scip, file, "BV", "Bound", (int) maxnamelen);
         printRecord(scip, file, varname, "", maxnamelen);
         SCIPinfoMessage(scip, file, "\n");
         continue;
      }

      /* take care of free variables */
      if ( SCIPisInfinity(scip, -lb) && SCIPisInfinity(scip, ub) )
      {
         /* variable is free */
         printStart(scip, file, "FR", "Bound", (int) maxnamelen);
         printRecord(scip, file, varname, "", maxnamelen);
         SCIPinfoMessage(scip, file, "\n");
         continue;
      }

      /* take care of fixed variables */
      if ( SCIPisEQ(scip, lb, ub) )
      {
         /* variable is fixed */
         snprintf(valuestr, MPS_MAX_VALUELEN, "%25.15g", lb);
         printStart(scip, file, "FX", "Bound", (int) maxnamelen);
         printRecord(scip, file, varname, valuestr, maxnamelen);
         SCIPinfoMessage(scip, file, "\n");
         continue;
      }

      /* print lower bound */
      if ( SCIPisInfinity(scip, -lb) )
      {
         assert( !SCIPisInfinity(scip, ub) );
         printStart(scip, file, "MI", "Bound", (int) maxnamelen);
         printRecord(scip, file, varname, "", maxnamelen);
         SCIPinfoMessage(scip, file, "\n");
      }
      else
      {
         if ( SCIPisZero(scip, lb) )
         {
            /* variables are nonnegative by default - so we skip these variables */
            if ( SCIPisInfinity(scip, ub) )
               continue;
            lb = 0.0;
         }
         else
         {
            snprintf(valuestr, MPS_MAX_VALUELEN, "%25.15g", lb);
            printStart(scip, file, "LO", "Bound", (int) maxnamelen);
            printRecord(scip, file, varname, valuestr, maxnamelen);
            SCIPinfoMessage(scip, file, "\n");
         }
      }

      /* print upper bound as far this one is not infinity */
      if( !SCIPisInfinity(scip, ub) )
      {
         snprintf(valuestr, MPS_MAX_VALUELEN, "%25.15g", ub);
         printStart(scip, file, "UP", "Bound", (int) maxnamelen);
         printRecord(scip, file, varname, valuestr, maxnamelen);
         SCIPinfoMessage(scip, file, "\n");
      }
   }
   
   /* output aggregated variables as 'free' */
   for (v = 0; v < naggvars; ++v)
   {
      var = aggvars[v];
      assert( var != NULL );

      /* get variable name */
      varname = varnames[nvars + v];
      assert( strncmp(varname, SCIPvarGetName(var), maxnamelen) == 0 );   
      
      /* variable is free */
      printStart(scip, file, "FR", "Bound", (int) maxnamelen);
      printRecord(scip, file, varname, "", maxnamelen);
      SCIPinfoMessage(scip, file, "\n");
   }
}


/*
 * Callback methods of reader
 */

/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeMps NULL


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadMps)
{  /*lint --e{715}*/
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIP_CALL( readMps(scip, filename) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteMps)
{  /*lint --e{715}*/
   int c, v, i;
   char* namestr;

   SCIP_CONS** relconss;
   int nrelconss;

   SCIP_CONS* cons = NULL;
   const char* consname;
   const char** consnames;

   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;

   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real* rhss;
   SCIP_Real value;

   SCIP_VAR* var = NULL;
   const char* varname;
   const char** varnames;
   
   char valuestr[MPS_MAX_VALUELEN] = { '\0' };

   SCIP_CONS** consSOS1;
   SCIP_CONS** consSOS2;
   int nConsSOS1 = 0;
   int nConsSOS2 = 0;

   SCIP_HASHMAP* varnameHashmap;           /* hash map from SCIP_VAR* to variable name */
   SPARSEMATRIX* matrix;
   
   SCIP_VAR** aggvars;
   int naggvars = 0;
   int saggvars;
   SCIP_HASHTABLE* varAggregatedHash;


   SCIP_VAR** consvars;
   int nconsvars;
   SCIP_Real* vals;
   SCIP_Longint* weights;

   SCIP_Bool needRANGES = FALSE;
   unsigned int maxnamelen = 0;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   
   /* check if the variable names are not too long and build the "variable" -> "variable name" hash map */
   SCIP_CALL( checkVarnames(scip, vars, nvars, &maxnamelen, &varnames, &varnameHashmap) );
 
   /* check if the constraint names are too long and build the constraint names */
   SCIP_CALL( checkConsnames(scip, conss, nconss, &relconss, &nrelconss, transformed, &maxnamelen, &consnames) );

   conss = relconss;
   nconss = nrelconss;

   /* collect SOS constraints in array for later output */
   SCIP_CALL( SCIPallocBufferArray(scip, &consSOS1, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consSOS2, nconss) );

   /* create hashtable for storing aggregated variables */
   saggvars = nvars;
   SCIP_CALL( SCIPallocBufferArray(scip, &aggvars, saggvars) );
   SCIP_CALL( SCIPhashtableCreate(&varAggregatedHash, SCIPblkmem(scip), 1000, hashGetKeyVar, hashKeyEqVar, hashKeyValVar, 
         NULL) );
   
   /* initialize sparse matrix */
   SCIP_CALL( initializeMatrix(scip, &matrix, (nvars * 2)) );
   assert( matrix->sentries >= nvars );

   /* initialize rhs vactor */
   SCIP_CALL( SCIPallocBufferArray(scip, &rhss, nconss) );
   
   /* print statistics as comment to file stream */
   SCIPinfoMessage(scip, file, "* SCIP STATISTICS\n");
   SCIPinfoMessage(scip, file, "*   Problem name     : %s\n", name);
   SCIPinfoMessage(scip, file, "*   Variables        : %d (%d binary, %d integer, %d implicit integer, %d continuous)\n",
      nvars, nbinvars, nintvars, nimplvars, ncontvars);
   SCIPinfoMessage(scip, file, "*   Constraints      : %d\n", nconss);
   SCIPinfoMessage(scip, file, "*   Obj. scale       : %.15g\n", objscale);
   SCIPinfoMessage(scip, file, "*   Obj. offset      : %.15g\n", objoffset);

   /* print NAME of the problem */
   SCIPinfoMessage(scip, file, "%-14s%s\n", "NAME", name);

   /* start ROWS section */
   SCIPinfoMessage(scip, file, "ROWS\n");
   
   /* print row type for the objective function */
   printStart(scip, file, "N", "Obj", -1);
   SCIPinfoMessage(scip, file, "\n");

   /* first fill the matrix with the objective coefficients */
   for( v = 0; v < nvars; ++v )
   {
      /* take care of the objective entry */
      var = vars[v];
      value = SCIPvarGetObj(var);

      if( !SCIPisZero(scip, value) )
      { 
         /* convert maximization problem into minimization since MPS format the objective is to minimize */
         if( objsense == SCIP_OBJSENSE_MAXIMIZE )
            value *= -1.0;
       
         assert( matrix->nentries < matrix->sentries );
         
         matrix->values[matrix->nentries] = value;
         matrix->columns[matrix->nentries] = var;
         matrix->rows[matrix->nentries] = "Obj";
         matrix->nentries++;
      }
   }

   /* loop over all constraints */
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert( cons != NULL);

      /* in case the transformed problems is written only constraint are
       * posted which are enabled in the current node; the conss array
       * should only contain relevant constraints since these are collected
       * in the beginning in the methode checkConsnames()  */
      assert( !transformed || SCIPconsIsEnabled(cons) );

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );
      
      conshdlrname = SCIPconshdlrGetName(conshdlr);

      /* construct constraint name */
      consname = consnames[c];
      assert( 0 == strncmp(consname, SCIPconsGetName(cons), maxnamelen) );
     
      if( strcmp(conshdlrname, "linear") == 0 )
      {
         lhs = SCIPgetLhsLinear(scip, cons);
         rhs = SCIPgetRhsLinear(scip, cons);

         /* there is nothing to do if the left hand side is minus infinity and the right side is infinity */
         if( !SCIPisInfinity(scip, -lhs) || !SCIPisInfinity(scip, rhs) )
         {
            if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) && !SCIPisEQ(scip, lhs, rhs) )
               needRANGES = TRUE;

            /* print row entry */
            printRowType(scip, file, lhs, rhs, consname);

            if( SCIPisInfinity(scip, rhs) )
               rhss[c] = lhs;
            else
               rhss[c] = rhs;
            
            assert( !SCIPisInfinity(scip, rhss[c]) );
            
            /* compute column entries */
            SCIP_CALL( getLinearCoeffs(scip, consname, 
                  SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons),
                  SCIPgetNVarsLinear(scip, cons), transformed, matrix, &rhss[c]) );
         }
      }
      else if( strcmp(conshdlrname, "setppc") == 0 )
      {
         
         /* print row entry */
         switch ( SCIPgetTypeSetppc(scip, cons) )
         {
         case SCIP_SETPPCTYPE_PARTITIONING :
            printRowType(scip, file, 1.0, 1.0, consname);
            break;
         case SCIP_SETPPCTYPE_PACKING :
            printRowType(scip, file, -SCIPinfinity(scip), 1.0, consname);
            break;
         case SCIP_SETPPCTYPE_COVERING :
            printRowType(scip, file, 1.0, SCIPinfinity(scip), consname);
            break;
         }

         rhss[c] = 1.0;
         
         /* compute column entries */
         SCIP_CALL( getLinearCoeffs(scip, consname, SCIPgetVarsSetppc(scip, cons), NULL,
               SCIPgetNVarsSetppc(scip, cons), transformed, matrix, &rhss[c]) );
      }
      else if( strcmp(conshdlrname, "logicor") == 0 )
      {
         /* print row entry */
         printRowType(scip, file, 1.0, SCIPinfinity(scip), consname);
         
         rhss[c] = 1.0;

         /* compute column entries */
         SCIP_CALL( getLinearCoeffs(scip, consname, SCIPgetVarsLogicor(scip, cons), NULL,
               SCIPgetNVarsLogicor(scip, cons), transformed, matrix, &rhss[c]) );
      }
      else if( strcmp(conshdlrname, "knapsack") == 0 )
      {
         /* print row entry */
         printRowType(scip, file, -SCIPinfinity(scip), (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons), consname);

         nconsvars = SCIPgetNVarsKnapsack(scip, cons);
         weights = SCIPgetWeightsKnapsack(scip, cons);

         /* copy Longint array to SCIP_Real array */
         SCIP_CALL( SCIPallocBufferArray(scip, &vals, nconsvars ) );
         for( i = 0; i < nconsvars; ++i )
            vals[i] = weights[i];

         rhss[c] = SCIPgetCapacityKnapsack(scip, cons);

         /* compute column entries */
         SCIP_CALL( getLinearCoeffs(scip, consname, SCIPgetVarsKnapsack(scip, cons), vals,
               nconsvars, transformed, matrix, &rhss[c]) );

         SCIPfreeBufferArray(scip, &vals);
      }
      else if( strcmp(conshdlrname, "varbound") == 0 )
      {
         lhs = SCIPgetLhsVarbound(scip, cons);
         rhs = SCIPgetRhsVarbound(scip, cons);

         /* there is nothing to do if the left hand side is minus infinity and the right side is infinity */
         if( !SCIPisInfinity(scip, -lhs) || !SCIPisInfinity(scip, rhs) )
         {
            if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) && !SCIPisEQ(scip, lhs, rhs) )
               needRANGES = TRUE;

            /* print row entry */
            printRowType(scip, file, lhs, rhs, consname);

            /* allocat memory */
            SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2) );
            SCIP_CALL( SCIPallocBufferArray(scip, &vals, 2) );
      
            consvars[0] = SCIPgetVarVarbound(scip, cons);
            consvars[1] = SCIPgetVbdvarVarbound(scip, cons);
      
            vals[0] = 1.0;
            vals[1] = SCIPgetVbdcoefVarbound(scip, cons);
      
            if( SCIPisInfinity(scip, rhs) )
               rhss[c] = lhs;
            else
               rhss[c] = rhs;
            
            assert( !SCIPisInfinity(scip, rhss[c]) ); 

            /* compute column entries */
            SCIP_CALL( getLinearCoeffs(scip, consname, consvars, vals, 2, transformed, matrix, &rhss[c]) );
            
            SCIPfreeBufferArray(scip, &consvars);
            SCIPfreeBufferArray(scip, &vals);
         }
      }
      else if ( strcmp(conshdlrname, "SOS1") == 0 )
      {
	 /* store constraint */
	 consSOS1[nConsSOS1++] = cons;

         /* check for aggregated variables in SOS1 constraints for later output
          * of aggregations as linear constraints */
         consvars = SCIPgetVarsSOS1(scip, cons);
         nconsvars = SCIPgetNVarsSOS1(scip, cons);

         /* SOS constraint do not have a right hand side; mark this with SCIPinfinity(scip) */
         rhss[c] = SCIPinfinity(scip); 

         SCIP_CALL( collectAggregatedVars(scip, consvars, nconsvars, &aggvars, &naggvars, &saggvars, varAggregatedHash) );
      }
      else if ( strcmp(conshdlrname, "SOS2") == 0 )
      {
	 /* store constraint */
	 consSOS2[nConsSOS2++] = cons;

         /* check for aggregated variables in SOS2 constraints for later output
          * aggregations as linear constraints */
         consvars = SCIPgetVarsSOS2(scip, cons);
         nconsvars = SCIPgetNVarsSOS2(scip, cons);
         
         /* SOS constraint do not have a right hand side; mark this with SCIPinfinity(scip) */
         rhss[c] = SCIPinfinity(scip); 
         
         SCIP_CALL( collectAggregatedVars(scip, consvars, nconsvars, &aggvars, &naggvars, &saggvars, varAggregatedHash) );
      }
      else
      {
         /* unknown constraint type; mark this with SCIPinfinity(scip) */
         rhss[c] = SCIPinfinity(scip); 
         
         SCIPwarningMessage("constraint handler <%s> cannot print requested format\n", conshdlrname );
      }
   }

   /* free hash table */
   SCIPhashtableFree(&varAggregatedHash);

   if( naggvars > 0 )
   {
      /* construct variables name of the needed aggregated variables and
       * the constraint names for the aggregation constraints */

      /* realloc memory */
      SCIP_CALL( SCIPreallocBufferArray(scip, &consnames, nconss + naggvars) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &rhss, nconss + naggvars) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &varnames, nvars + naggvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 1) );

      for (c = 0; c < naggvars; ++c )
      {
	 size_t l;

         /* create variable name */
         var = aggvars[c];
         
         if( strlen(SCIPvarGetName(var)) >= MPS_MAX_NAMELEN )
            maxnamelen = MPS_MAX_NAMELEN - 1;
         else
	 {
	    l = strlen(SCIPvarGetName(var));
            maxnamelen = MAX(maxnamelen, l);
	 }
         
         SCIP_CALL( SCIPallocBufferArray(scip, &namestr, MPS_MAX_NAMELEN) );
         snprintf(namestr,  MPS_MAX_NAMELEN, "%s", SCIPvarGetName(var) );
            
         /* insert variable with variable name into hash map */
         varnames[nvars + c] = namestr;
         assert( !SCIPhashmapExists(varnameHashmap, var) );
         SCIP_CALL( SCIPhashmapInsert(varnameHashmap, var, (void*) (size_t) namestr) );

         /* output row type (it is an equation) */
         SCIP_CALL( SCIPallocBufferArray(scip, &namestr, MPS_MAX_VALUELEN) );
         snprintf(namestr, MPS_MAX_VALUELEN, "aggr%d", c );
         printRowType(scip, file, 1.0, 1.0, namestr);

	 l = strlen(namestr);
         maxnamelen = MAX(maxnamelen, l);
         
         consnames[nconss + c] = namestr;

         consvars[0] = aggvars[c];
         rhss[nconss + c] = 0.0;
         
         /* compute column entries */
         SCIP_CALL( getLinearCoeffs(scip, namestr, consvars, NULL, 1, transformed, matrix, &rhss[nconss + c]) );
         
         /* add the aggregated variables to the sparse matrix */
         SCIP_CALL( checkSparseMatrixCapacity(scip, matrix, 1) );
         matrix->values[matrix->nentries] = -1.0;
         matrix->columns[matrix->nentries] = aggvars[c];
         matrix->rows[matrix->nentries] = namestr;
         matrix->nentries++;
      }
      SCIPfreeBufferArray(scip, &consvars);
   }
   
   /* output COLUMNS section */
   printColumnSection(scip, file, matrix, varnameHashmap, maxnamelen);
   
   /* output RHS section */
   printRhsSection(scip, file, conss, nconss, consnames, rhss, transformed, maxnamelen);
   
   if( needRANGES )
   {
      /* output RANGES section */
      printRangeSection(scip, file, conss, nconss, consnames, transformed, maxnamelen);
   }
   
   /* output BOUNDS section */
   printBoundSection(scip, file, vars, nvars, aggvars, naggvars, transformed, varnames, maxnamelen);


   /* print SOS section */
   if ( nConsSOS1 > 0 || nConsSOS2 > 0 )
   {
      SCIP_Real* sosweights;

      SCIPinfoMessage(scip, file, "SOS\n");
      SCIPdebugMessage("start printing SOS section\n");

      SCIP_CALL( SCIPallocBufferArray(scip, &namestr, MPS_MAX_NAMELEN) );

      /* first output SOS1 constraints */
      for (c = 0; c < nConsSOS1; ++c)
      {
	 cons = consSOS1[c];
	 consvars = SCIPgetVarsSOS1(scip, cons);
	 nconsvars = SCIPgetNVarsSOS1(scip, cons);
	 sosweights = SCIPgetWeightsSOS1(scip, cons);
         snprintf(namestr, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );
         
         printStart(scip, file, "S1", namestr, -1);
	 SCIPinfoMessage(scip, file, "\n");

	 for (v = 0; v < nconsvars; ++v)
	 {
            /* get variable name */
            assert ( SCIPhashmapExists(varnameHashmap, consvars[v]) );
            varname = (const char*) (size_t) SCIPhashmapGetImage(varnameHashmap, consvars[v]);
            
            printStart(scip, file, "", varname, (int) maxnamelen);
            
	    if ( sosweights != NULL )
               snprintf(valuestr, MPS_MAX_VALUELEN, "%25.15g", sosweights[v]);
	    else
               snprintf(valuestr, MPS_MAX_VALUELEN, "%25d ", v);

            SCIPinfoMessage(scip, file, "%25s\n", valuestr);
	 }
      }

      /* next output SOS2 constraints */
      for( c = 0; c < nConsSOS2; ++c )
      {
	 cons = consSOS2[c];
	 consvars = SCIPgetVarsSOS2(scip, cons);
	 nconsvars = SCIPgetNVarsSOS2(scip, cons);
	 sosweights = SCIPgetWeightsSOS2(scip, cons);
         snprintf(namestr, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );

         printStart(scip, file, "S2", namestr, -1);
	 SCIPinfoMessage(scip, file, "\n");

	 for (v = 0; v < nconsvars; ++v)
	 {
            /* get variable name */
            assert ( SCIPhashmapExists(varnameHashmap, consvars[v]) );
            varname = (const char*) (size_t) SCIPhashmapGetImage(varnameHashmap, consvars[v]);

            printStart(scip, file, "", varname, (int) maxnamelen);

	    if ( sosweights != NULL )
               snprintf(valuestr, MPS_MAX_VALUELEN, "%25.15g", sosweights[v]);
	    else
               snprintf(valuestr, MPS_MAX_VALUELEN, "%25d ", v);

	    SCIPinfoMessage(scip, file, "%25s\n", valuestr);
	 }
      }
      SCIPfreeBufferArray(scip, &namestr);
   }

   /* free variable and constraint name array */
   for( v = 0; v < nvars + naggvars; ++v )
      SCIPfreeBufferArray(scip, &varnames[v]);
   SCIPfreeBufferArray(scip, &varnames);

   for( c = 0; c < nconss + naggvars; ++c )
   {
      /* in case the transformed problems is written only constraint are
       * posted which are enabled in the current node; the conss array
       * should only contain relevant constraints since these are collected
       * in the beginning in the methode checkConsnames()  */
      assert( !transformed || SCIPconsIsEnabled(cons) );
      
      SCIPfreeBufferArray(scip, &consnames[c]);
   }
   SCIPfreeBufferArray(scip, &consnames);
   
   /* free matrix data structure */
   freeMatrix(scip, matrix);

   /* free variable hashmap */
   SCIPhashmapFree(&varnameHashmap);

   SCIPfreeBufferArray(scip, &aggvars);
   SCIPfreeBufferArray(scip, &rhss);
   SCIPfreeBufferArray(scip, &relconss);

   /* free buffer arrays for SOS1 and SOS2 */
   SCIPfreeBufferArray(scip, &consSOS1);
   SCIPfreeBufferArray(scip, &consSOS2);

   /* print end of data line */
   SCIPinfoMessage(scip, file, "ENDATA");
   
   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/*
 * mps file reader specific interface methods
 */

/** includes the mps file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderMps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   /* create mps reader data */
   readerdata = NULL;

   /* include mps reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerFreeMps, readerReadMps, readerWriteMps, readerdata) );

   /* add mps reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/mpsreader/dynamicconss", "should model constraints be subject to aging?",
         NULL, FALSE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/mpsreader/dynamiccols", "should columns be added and removed dynamically to the LP?",
         NULL, FALSE, FALSE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/mpsreader/dynamicrows", "should rows be added and removed dynamically to the LP?",
         NULL, FALSE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2007 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: reader_mps.c,v 1.85 2008/03/19 19:04:03 bzfpfets Exp $"

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
#define MPS_MAX_NAMELEN   255
#define MPS_MAX_VALUELEN   25
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
   char                 probname[MPS_MAX_LINELEN];
   char                 objname [MPS_MAX_LINELEN];
   SCIP_Bool            isinteger;
   SCIP_Bool            isnewformat;
   SCIP_Bool            semicontwarning;
};
typedef struct MpsInput MPSINPUT;



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

   strcpy(mpsi->probname, probname);
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

   strcpy(mpsi->objname, objname);
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
   int                   verblevel
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
   char          colname[MPS_MAX_LINELEN] = { '\0' };
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

         strcpy(colname, mpsinputField1(mpsi));

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
   char        rhsname[MPS_MAX_LINELEN] = { '\0' };
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
         strcpy(rhsname, mpsinputField1(mpsi));

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
   char        rngname[MPS_MAX_LINELEN] = { '\0' };
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
         strcpy(rngname, mpsinputField1(mpsi));

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
   char        bndname[MPS_MAX_LINELEN] = { '\0' };
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
         strcpy(bndname, mpsinputField2(mpsi));

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
   char name[SCIP_MAXSTRLEN] = { '\0' };
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
	    strncpy(name, mpsinputField2(mpsi), SCIP_MAXSTRLEN);
	 else
	 {
	    /* create new name */
	    snprintf(name, SCIP_MAXSTRLEN, "SOS%d", ++cnt);
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

   if (NULL == (fp = SCIPfopen(filename, "r")))
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPerrorMessage("%s: %s\n", filename, strerror(errno));
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
   return SCIPvarGetIndex((SCIP_VAR*) key);
}


/* computes the field width such that the output file is nicely arranged */
static
int computeFieldWidth(
   int                   width              /**< required width */
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
   int                   maxnamelen          /**< maximum name lenght */
   )
{
   int fieldwidth;
   char format[16];

   assert( scip != NULL );
   assert( col1 != NULL );
   assert( col2 != NULL );
   assert( strlen(col1) <= MPS_MAX_NAMELEN );
   assert( strlen(col2) <= MPS_MAX_VALUELEN );
   assert( maxnamelen > 0 );

   fieldwidth = computeFieldWidth(maxnamelen);
   sprintf(format, " %%-%ds %%%ds ", fieldwidth, MPS_MAX_VALUELEN);

   SCIPinfoMessage(scip, file, format, col1, col2);
}

/* output two strings in columns 1 (width 2) and 2 (width 8) */
static
void printStart(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           col1,               /**< column 1 */
   const char*           col2,               /**< column 2 */
   int                   maxnamelen          /**< maximum name lenght */
   )
{
   int fieldwidth;
   char format[16];

   assert( scip != NULL );
   assert( col1 != NULL );
   assert( col2 != NULL );
   assert( strlen(col1) <= 2 );
   assert( strlen(col2) <= MPS_MAX_NAMELEN );
   assert( maxnamelen == -1 || maxnamelen > 0 );

   if( maxnamelen == -1 )
   {
      /* format does not matter */
      sprintf(format, " %%-2.2s %%-s ");
   }
   else
   {
      fieldwidth = computeFieldWidth(maxnamelen);
      sprintf(format, " %%-2.2s %%-%ds ", fieldwidth);
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
   int                maxnamelen          /**< maximum name lenght */
   )
{
   char valuestr[MPS_MAX_VALUELEN + 1];

   assert( scip != NULL );
   assert( recordcnt != NULL );
   assert( *recordcnt >= 0 && *recordcnt < 2 );

   if( !SCIPisZero(scip, value) )
   {
      snprintf(valuestr, MPS_MAX_VALUELEN + 1, "%25.15g", value);

      if ( *recordcnt == 0 )
      {
         /* finish last line which already contains two records */
         printStart(scip, file, "", varname, maxnamelen);
         *recordcnt = 0;
      }

      printRecord(scip, file, consname, valuestr, maxnamelen);
      (*recordcnt)++;

      if ( *recordcnt == 2 )
      {
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
      sprintf(rowtype, "%s", "E");
   else
   {
      /* in case the right hand side and the left hand side are not infinity we print a
       * less or equal constraint and put the right hand side in the RHS section and the
       * left hand side (hidden) in the RANGE section */
      if( !SCIPisInfinity(scip, rhs) )
         sprintf(rowtype, "%s", "L");
      else
      {
         assert( !SCIPisInfinity(scip, -lhs) );
         sprintf(rowtype, "%s", "G");
      }
   }

   printStart(scip, file, rowtype, name, -1);
   SCIPinfoMessage(scip, file, "\n");
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

/** computes the coefficient for the given variable and linear constraint information */
static
SCIP_RETCODE getLinearCoeff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variables of interest */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of coefficients values (or NULL if all coefficient values are 1) */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             transformed,        /**< transformed constraint? */
   SCIP_Real*            coeff               /**< pointer to store the coefficient */
   )
{
   int v;
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   int nactivevars;
   SCIP_Real activeconstant = 0.0;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( nvars > 0 );
   assert( coeff != NULL );
   assert( *coeff == 0.0 );

   /* duplicate variable and value array */
   nactivevars = nvars;
   SCIPduplicateBufferArray(scip, &activevars, vars, nactivevars );

   if( vals != NULL )
      SCIPduplicateBufferArray(scip, &activevals, vals, nactivevars );
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nactivevars) );

      for( v = 0; v < nactivevars; ++v )
         activevals[v] = 1.0;
   }

   /* retransform given variables to active variables */
   SCIP_CALL( getActiveVariables(scip, activevars, activevals, &nactivevars, &activeconstant, transformed) );

   /* search for the variables in the active variables array */
   for( v = 0; v < nactivevars; ++v )
   {
      if( activevars[v] == var )
      {
         *coeff = activevals[v];
         break;
      }
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);

   return SCIP_OKAY;
}


/** computes the right side for the given linear constraint information */
static
SCIP_RETCODE getLinearRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of coefficients values (or NULL if all coefficient values are 1) */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             transformed,        /**< transformed constraint? */
   SCIP_Real*            rhs                 /**< pointer to store the right hand side */
   )
{
   int v;
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   int nactivevars;
   SCIP_Real activeconstant = 0.0;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( nvars > 0 );
   assert( rhs != NULL );

   /* duplicate variable and value array */
   nactivevars = nvars;
   SCIPduplicateBufferArray(scip, &activevars, vars, nactivevars );

   if( vals != NULL )
      SCIPduplicateBufferArray(scip, &activevals, vals, nactivevars );
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nactivevars) );

      for( v = 0; v < nactivevars; ++v )
         activevals[v] = 1.0;
   }

   /* retransform given variables to active variables */
   SCIP_CALL( getActiveVariables(scip, activevars, activevals, &nactivevars, &activeconstant, transformed) );

   /* compute rhs */
   *rhs -= activeconstant;

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);

   return SCIP_OKAY;
}


/** check whether given variables are aggregated and put them into an array without duplication */
static
SCIP_RETCODE collectAggregatedVars(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_Bool          transformed,        /**< TRUE iff problem is the transformed problem */
   int                nvars,              /**< number of mutable variables in the problem */
   SCIP_VAR**         vars,               /**< variable array */
   int*               nAggregatedVars,    /**< number of aggregated variables on output */
   SCIP_VAR**         aggregatedVars,     /**< array storing the aggregated variables on output */
   SCIP_HASHTABLE*    varAggregated,      /**< hashtable for checking duplicates */
   SCIP_Real*         aggregatedConstant  /**< constant of aggregation */
   )
{
   int j;
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   int nactivevars;
   SCIP_Real activeconstant = 0.0;

   assert( scip != NULL );
   assert( nAggregatedVars != NULL );
   assert( varAggregated != NULL );
   assert( aggregatedVars != NULL );
   assert( aggregatedConstant != NULL );

   /* reserve space for variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &activevars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nvars) );

   /* check variables */
   for (j = 0; j < nvars; ++j)
   {
      SCIP_VARSTATUS varStatus;
      SCIP_VAR* var;

      var = vars[j];
      varStatus = SCIPvarGetStatus(var);

      /* collect aggregated variables in a list */
      if ( varStatus == SCIP_VARSTATUS_AGGREGATED || varStatus == SCIP_VARSTATUS_MULTAGGR || varStatus == SCIP_VARSTATUS_NEGATED )
      {
	 if ( ! SCIPhashtableExists(varAggregated, (void*) var) )
	 {
	    /* set up list to obtain substitution variables */
	    nactivevars = 1;

	    activevars[0] = var;
	    activevals[0] = 1.0;
	    activeconstant = 0.0;

	    /* retransform given variables to active variables */
	    SCIP_CALL( getActiveVariables(scip, activevars, activevals, &nactivevars, &activeconstant, transformed) );

	    aggregatedVars[*nAggregatedVars] = var;
	    aggregatedConstant[(*nAggregatedVars)++] = activeconstant;
	    SCIP_CALL( SCIPhashtableInsert(varAggregated, (void*) var) );
	 }
      }
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);

   return SCIP_OKAY;
}

/** write variables for aggregated constraints in SOS constraints */
static
SCIP_RETCODE printAggregatedConsVar(
   SCIP*              scip,               /**< SCIP data structure */
   FILE*              file,               /**< output file (or NULL for standard output) */
   SCIP_VAR*          var,                /**< variable */
   const char*        varname,            /**< variable name */
   SCIP_Bool          transformed,        /**< TRUE iff problem is the transformed problem */
   int                nvars,              /**< number of mutable variables in the problem */
   int                nConsSOS1,          /**< number of SOS1 constraints in consSOS1 array */
   SCIP_CONS**        consSOS1,           /**< array of SOS1 constraints */
   int                nConsSOS2,          /**< number of SOS2 constraints in consSOS2 array */
   SCIP_CONS**        consSOS2,           /**< array of SOS1 constraints */
   int                nAggregatedVars,    /**< number of aggregated variables */
   SCIP_VAR**         aggregatedVars,     /**< array storing the aggregated variables */
   int*               recordcnt,          /**< variable counter */
   int                maxnamelen          /**< pointer to store rhe maximum name lenght */
   )
{
   int c, v;

   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   int nactivevars;
   SCIP_Real activeconstant = 0.0;
   char consname[MPS_MAX_NAMELEN + 1];

   assert( scip != NULL );
   assert( consSOS1 != NULL );
   assert( consSOS2 != NULL );

   /* write aggregation constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &activevars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nvars) );

   for (c = 0; c < nAggregatedVars; ++c)
   {
      /* set up list to obtain substitution variables */
      nactivevars = 1;

      activevars[0] = aggregatedVars[c];
      activevals[0] = 1.0;
      activeconstant = 0.0;

      /* retransform given variables to active variables */
      SCIP_CALL( getActiveVariables(scip, activevars, activevals, &nactivevars, &activeconstant, transformed) );

      activevals[nactivevars] = -1.0;
      activevars[nactivevars] = aggregatedVars[c];
      ++nactivevars;

      snprintf(consname, MPS_MAX_NAMELEN, "aggr%d", c);

      /* check for variable */
      for (v = 0; v < nactivevars; ++v)
      {
	 if ( activevars[v] == var )
	    printEntry(scip, file, varname, consname, activevals[v], recordcnt, maxnamelen);
      }
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);

   return SCIP_OKAY;
}

/** method check if the variable names are not longer than LP_MAX_NAMELEN */
static
void checkVarnames(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_VAR**         vars,               /**< array of variables */
   int                nvars,              /**< number of variables */
   unsigned int*      maxnamelen          /**< pointer to store rhe maximum name lenght */
   )
{
   int v;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( maxnamelen != NULL );

   /* check if the variable names are not to long */
   for( v = 0; v < nvars; ++v )
   {
      if( strlen(SCIPvarGetName(vars[v])) > MPS_MAX_NAMELEN )
      {
         SCIPwarningMessage("there is a variable name which has to be cut down to %d characters; LP might be corrupted\n", MPS_MAX_NAMELEN);
         (*maxnamelen) = MPS_MAX_NAMELEN;
         return;
      }

      (*maxnamelen) = MAX(*maxnamelen, strlen(SCIPvarGetName(vars[v])));
   }
}

/** method check if the constraint names are not longer than LP_MAX_NAMELEN */
static
void checkConsnames(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONS**        conss,              /**< array of constraints */
   int                nconss,             /**< number of constraints */
   SCIP_Bool          transformed,        /**< TRUE iff problem is the transformed problem */
   unsigned int*      maxnamelen          /**< pointer to store rhe maximum name lenght */
   )
{
   int c;

   assert( scip != NULL );
   assert( conss != NULL );
   assert( maxnamelen != NULL );

   for( c = 0; c < nconss; ++c )
   {
      if( strlen(SCIPconsGetName(conss[c])) > MPS_MAX_NAMELEN )
      {
         SCIPwarningMessage("there is a constraint name which has to be cut down to %d characters; LP might be corrupted\n",
            MPS_MAX_NAMELEN);

         (*maxnamelen) = MPS_MAX_NAMELEN;
         return;
      }

      (*maxnamelen) = MAX(*maxnamelen, strlen(SCIPconsGetName(conss[c])));
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
   int recordcnt;

   SCIP_CONS* cons = NULL;
   char consname[MPS_MAX_LINELEN + 1];
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;

   int nint;
   SCIP_Real lhs;
   SCIP_Real rhs;

   SCIP_VAR* var = NULL;
   char varname[MPS_MAX_LINELEN + 1];
   SCIP_Real value;
   char valuestr[MPS_MAX_VALUELEN + 1];
   SCIP_Real lb;
   SCIP_Real ub;

   SCIP_CONS** linears;
   SCIP_CONS** setppcs;
   SCIP_CONS** logicors;
   SCIP_CONS** knapsacks;
   SCIP_CONS** varbounds;
   SCIP_CONS** consSOS1;
   SCIP_CONS** consSOS2;
   int nlinears;
   int nsetppcs;
   int nlogicors;
   int nknapsacks;
   int nvarbounds;
   int nConsSOS1 = 0;
   int nConsSOS2 = 0;

   SCIP_VAR** aggregatedVars;
   int nAggregatedVars = 0;
   SCIP_HASHTABLE* varAggregatedHash;
   SCIP_Real* aggregatedConstant;

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

   /* check if the variable names are not too long */
   checkVarnames(scip, vars, nvars, &maxnamelen);
   /* check if the constraint names are too long */
   checkConsnames(scip, conss, nconss, transformed, &maxnamelen);

   /* before we start printing the problem in MPS format we sort the constraints by their type */
   SCIP_CALL( SCIPallocBufferArray(scip, &linears, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &setppcs, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &logicors, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &knapsacks, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varbounds, nconss) );

   /* collect SOS constraints in array for later output */
   SCIP_CALL( SCIPallocBufferArray(scip, &consSOS1, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consSOS2, nconss) );

   nlinears = 0;
   nsetppcs = 0;
   nlogicors = 0;
   nknapsacks = 0;
   nvarbounds = 0;

   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert( cons != NULL);

      /* in case the transformed is written only constraint are posted which are enabled in the current node */
      if( transformed && !SCIPconsIsEnabled(cons) )
         continue;

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);

      if( strcmp(conshdlrname, "linear") == 0 )
      {
         lhs = SCIPgetLhsLinear(scip, cons);
         rhs = SCIPgetRhsLinear(scip, cons);

         /* there nothing to do if the left hand side is minus infinity and the right side is infinity */
         if( !SCIPisInfinity(scip, -lhs) || !SCIPisInfinity(scip, rhs) )
         {
            if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) && !SCIPisEQ(scip, lhs, rhs) )
               needRANGES = TRUE;

            linears[nlinears] = cons;
            nlinears++;
         }
      }
      else if( strcmp(conshdlrname, "setppc") == 0 )
      {
         setppcs[nsetppcs] = cons;
         nsetppcs++;
      }
      else if( strcmp(conshdlrname, "logicor") == 0 )
      {
         logicors[nlogicors] = cons;
         nlogicors++;
      }
      else if( strcmp(conshdlrname, "knapsack") == 0 )
      {
         knapsacks[nknapsacks] = cons;
         nknapsacks++;
      }
      else if( strcmp(conshdlrname, "varbound") == 0 )
      {
         lhs = SCIPgetLhsVarbound(scip, cons);
         rhs = SCIPgetRhsVarbound(scip, cons);

         /* there nothing to do if the left hand side is minus infinity and the right side is infinity */
         if( !SCIPisInfinity(scip, -lhs) || !SCIPisInfinity(scip, rhs) )
         {
            if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) && !SCIPisEQ(scip, lhs, rhs) )
               needRANGES = TRUE;

            varbounds[nvarbounds] = cons;
            nvarbounds++;
         }
      }
      else if ( strcmp(conshdlrname, "SOS1") == 0 )
      {
	 /* store constraint */
	 consSOS1[nConsSOS1++] = cons;
      }
      else if ( strcmp(conshdlrname, "SOS2") == 0 )
      {
	 /* store constraint */
	 consSOS2[nConsSOS2++] = cons;
      }
      else
      {
         SCIPwarningMessage("constraint handler <%s> can not print requested format\n", conshdlrname );
      }
   }

   /* create hashtable for storing aggregated variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &aggregatedVars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &aggregatedConstant, nvars) );
   SCIP_CALL( SCIPhashtableCreate(&varAggregatedHash, SCIPblkmem(scip), 1000, hashGetKeyVar, hashKeyEqVar, hashKeyValVar) );

   /* check for aggregated variables in SOS1 constraints and output aggregations as linear constraints */
   for (c = 0; c < nConsSOS1; ++c)
   {
      cons = consSOS1[c];
      consvars = SCIPgetVarsSOS1(scip, cons);
      nconsvars = SCIPgetNVarsSOS1(scip, cons);

      SCIP_CALL( collectAggregatedVars(scip, transformed, nconsvars, consvars, &nAggregatedVars, aggregatedVars,
				       varAggregatedHash, aggregatedConstant) );
   }

   /* check for aggregated variables in SOS2 constraints and output aggregations as linear constraints */
   for (c = 0; c < nConsSOS2; ++c)
   {
      cons = consSOS2[c];
      consvars = SCIPgetVarsSOS2(scip, cons);
      nconsvars = SCIPgetNVarsSOS2(scip, cons);

      SCIP_CALL( collectAggregatedVars(scip, transformed, nconsvars, consvars, &nAggregatedVars, aggregatedVars,
				       varAggregatedHash, aggregatedConstant) );
   }

   /* now we start writing the problem */

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

   /* print ROWS section */
   SCIPinfoMessage(scip, file, "ROWS\n");

   printStart(scip, file, "N", "Obj", -1);
   SCIPinfoMessage(scip, file, "\n");

   /* take care of the linear constraints */
   for( c = 0; c < nlinears; ++c )
   {
      cons = linears[c];
      assert( cons != NULL);

      lhs = SCIPgetLhsLinear(scip, cons);
      rhs = SCIPgetRhsLinear(scip, cons);

      /* there nothing to do if the left hand side is minus infinity and the right side infinity */
      assert( !SCIPisInfinity(scip, -lhs) || !SCIPisInfinity(scip, rhs) );
      snprintf(consname, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );
      printRowType(scip, file, lhs, rhs, consname);
   }

   /* take care of the setppc constraints */
   for( c = 0; c < nsetppcs; ++c )
   {
      cons = setppcs[c];
      snprintf(consname, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );

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
   }

   /* take care of the logicor constraints */
   for( c = 0; c < nlogicors; ++c )
   {
      cons = logicors[c];
      snprintf(consname, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );
      printRowType(scip, file, 1.0, SCIPinfinity(scip), consname);
   }

   /* take care of the knapsack constraints */
   for( c = 0; c < nknapsacks; ++c )
   {
      cons = knapsacks[c];
      snprintf(consname, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );
      printRowType(scip, file, -SCIPinfinity(scip), SCIPgetCapacityKnapsack(scip, cons), consname);
   }

   /* take care of the varbound constraints */
   for( c = 0; c < nvarbounds; ++c )
   {
      cons = varbounds[c];
      lhs = SCIPgetLhsVarbound(scip, cons);
      rhs = SCIPgetRhsVarbound(scip, cons);

      /* there nothing to do if the left hand side is minus infinity and the right side infinity */
      assert( !SCIPisInfinity(scip, -lhs) || !SCIPisInfinity(scip, rhs) );

      snprintf(consname, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );
      printRowType(scip, file, lhs, rhs, consname);
   }

   /* take care of the aggregated constraints for SOS constraints */
   for (c = 0; c < nAggregatedVars; ++c )
   {
      /* output type */
      snprintf(consname, MPS_MAX_NAMELEN, "aggr%d", c );
      printRowType(scip, file, aggregatedConstant[c], aggregatedConstant[c], consname);
   }

   /* print COLUMNS section */
   SCIPinfoMessage(scip, file, "COLUMNS\n");

   nint = nbinvars + nintvars + nimplvars;

   if ( nint > 0 )
   {
      printStart(scip, file, "", "INTSTART", maxnamelen);
      printRecord(scip, file, "'MARKER'", "", maxnamelen);
      printRecord(scip, file, "'INTORG'", "", maxnamelen);
      SCIPinfoMessage(scip, file, "\n");
   }

   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert( var != NULL );
      assert( (SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS && v < nint) ||
	      (v >= nint && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS) );

      recordcnt = 0;

      SCIPdebugMessage("create entries for variable <%s>\n", SCIPvarGetName(var));
      /* create variable name */
      snprintf(varname, MPS_MAX_NAMELEN, "%s", SCIPvarGetName(var) );

      /* first take care of the objective function */
      value = SCIPvarGetObj(var);

      /* convert maximization problem into minimization since MPS format the objective is to minimize */
      if ( objsense == SCIP_OBJSENSE_MAXIMIZE )
         value *= -1.0;

      printEntry( scip, file, varname, "Obj", value, &recordcnt, maxnamelen );

      /* take care of the linear constraints */
      for( c = 0; c < nlinears; ++c )
      {
         cons = linears[c];
         value = 0.0;
         snprintf(consname, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );

         SCIP_CALL( getLinearCoeff(scip, var, SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons),
				   SCIPgetNVarsLinear(scip, cons), transformed, &value) );

         printEntry( scip, file, varname, consname, value, &recordcnt, maxnamelen );
      }

      /* take care of the setppc constraints */
      for( c = 0; c < nsetppcs; ++c )
      {
         cons = setppcs[c];
         value = 0.0;
         snprintf(consname, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );

         SCIP_CALL( getLinearCoeff(scip, var, SCIPgetVarsSetppc(scip, cons), NULL,
				   SCIPgetNVarsSetppc(scip, cons), transformed, &value) );

         printEntry( scip, file, varname, consname, value, &recordcnt, maxnamelen );
      }

      /* take care of the logicor constraints */
      for( c = 0; c < nlogicors; ++c )
      {
         cons = logicors[c];
         value = 0.0;
         snprintf(consname, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );

         SCIP_CALL( getLinearCoeff(scip, var, SCIPgetVarsLogicor(scip, cons), NULL,
				   SCIPgetNVarsLogicor(scip, cons), transformed, &value) );

         printEntry( scip, file, varname, consname, value, &recordcnt, maxnamelen );
      }

      /* take care of the knapsack constraints */
      for( c = 0; c < nknapsacks; ++c )
      {
         cons = knapsacks[c];
         value = 0.0;
         snprintf(consname, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );

         nconsvars = SCIPgetNVarsKnapsack(scip, cons);
         weights = SCIPgetWeightsKnapsack(scip, cons);

         /* copy Longint array to SCIP_Real array */
         SCIP_CALL( SCIPallocBufferArray(scip, &vals, nconsvars ) );
         for( i = 0; i < nconsvars; ++i )
            vals[i] = weights[i];

         SCIP_CALL( getLinearCoeff(scip, var, SCIPgetVarsKnapsack(scip, cons), vals,
               nconsvars, transformed, &value) );

         printEntry( scip, file, varname, consname, value, &recordcnt, maxnamelen );

         SCIPfreeBufferArray(scip, &vals);
      }

      /* take care of the varbound constraints */
      for( c = 0; c < nvarbounds; ++c )
      {
         cons = varbounds[c];
         value = 0.0;
         snprintf(consname, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );

         /* allocat memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2) );
         SCIP_CALL( SCIPallocBufferArray(scip, &vals, 2) );

         consvars[0] = SCIPgetVarVarbound(scip, cons);
         consvars[1] = SCIPgetVbdvarVarbound(scip, cons);

         vals[0] = 1.0;
         vals[1] = SCIPgetVbdcoefVarbound(scip, cons);

         SCIP_CALL( getLinearCoeff(scip, var, consvars, vals, 2, transformed, &value) );

         printEntry( scip, file, varname, consname, value, &recordcnt, maxnamelen );

         /* free buffer array */
         SCIPfreeBufferArray(scip, &consvars);
         SCIPfreeBufferArray(scip, &vals);
      }

      /* take care of the aggregated constraints needed for SOS constraints */
      SCIP_CALL( printAggregatedConsVar(scip, file, var, varname, nvars, transformed, nConsSOS1, consSOS1, nConsSOS2, consSOS2,
            nAggregatedVars, aggregatedVars, &recordcnt, maxnamelen) );

      if( recordcnt == 1 )
         SCIPinfoMessage(scip, file, "\n");

      /* print end of integer variables */
      if( v + 1 == nint )
      {
         printStart(scip, file, "", "INTEND", maxnamelen);
         printRecord(scip, file, "'MARKER'", "", maxnamelen);
         printRecord(scip, file, "'INTEND'", "", maxnamelen);
         SCIPinfoMessage(scip, file, "\n", maxnamelen);
      }
   }

   for (v = 0; v < nfixedvars; ++v)
   {
      var = fixedvars[v];
      assert( var != NULL );

      SCIPdebugMessage("create entries for variable <%s>\n", SCIPvarGetName(var));

      recordcnt = 0;
      /* create variable name */
      snprintf(varname, MPS_MAX_NAMELEN, "%s", SCIPvarGetName(var) );

      /* take care of the aggregated constraints needed for SOS constraints */
      SCIP_CALL( printAggregatedConsVar(scip, file, var, varname, nvars, transformed, nConsSOS1, consSOS1, nConsSOS2, consSOS2,
            nAggregatedVars, aggregatedVars, &recordcnt, maxnamelen) );

      if( recordcnt == 1 )
         SCIPinfoMessage(scip, file, "\n");
   }

   /* print RHS section */
   SCIPinfoMessage(scip, file, "RHS\n");
   SCIPdebugMessage("start printing RHS section\n");

   recordcnt = 0;

   /* take care of the linear constraints */
   for( c = 0; c < nlinears; ++c )
   {
      cons = linears[c];
      snprintf(consname, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );

      value = SCIPgetRhsLinear(scip, cons);

      if( SCIPisInfinity(scip, value) )
         value = SCIPgetLhsLinear(scip, cons);

      assert( !SCIPisInfinity(scip, value) && !SCIPisInfinity(scip, -value) );

      SCIP_CALL( getLinearRhs(scip, SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons),
			      SCIPgetNVarsLinear(scip, cons), transformed, &value) );

      printEntry( scip, file, "RHS", consname, value, &recordcnt, maxnamelen );
   }

   /* take care of the setppc constraints */
   for( c = 0; c < nsetppcs; ++c )
   {
      cons = setppcs[c];
      snprintf(consname, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );
      value = 1.0;

      SCIP_CALL( getLinearRhs(scip, SCIPgetVarsSetppc(scip, cons), NULL,
            SCIPgetNVarsSetppc(scip, cons), transformed, &value) );

      printEntry( scip, file, "RHS", consname, value, &recordcnt, maxnamelen );
   }

   /* take care of the logicor constraints */
   for( c = 0; c < nlogicors; ++c )
   {
      cons = logicors[c];
      snprintf(consname, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );
      value = 1.0;

      SCIP_CALL( getLinearRhs(scip, SCIPgetVarsLogicor(scip, cons), NULL,
            SCIPgetNVarsLogicor(scip, cons), transformed, &value) );

      printEntry( scip, file, "RHS", consname, value, &recordcnt, maxnamelen );
   }

   /* take care of the knapsack constraints */
   for( c = 0; c < nknapsacks; ++c )
   {
      cons = knapsacks[c];
      snprintf(consname, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );

      value = SCIPgetCapacityKnapsack(scip, cons);;
      nconsvars = SCIPgetNVarsKnapsack(scip, cons);
      weights = SCIPgetWeightsKnapsack(scip, cons);

      /* copy Longint array to SCIP_Real array */
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, nconsvars ) );
      for( i = 0; i < nconsvars; ++i )
         vals[i] = weights[i];

      SCIP_CALL( getLinearRhs(scip, SCIPgetVarsKnapsack(scip, cons), vals,
            nconsvars, transformed, &value) );

      printEntry( scip, file, "RHS", consname, value, &recordcnt, maxnamelen );

      SCIPfreeBufferArray(scip, &vals);
   }

   /* take care of the varbound constraints */
   for( c = 0; c < nvarbounds; ++c )
   {
      cons = varbounds[c];
      value = 0.0;
      snprintf(consname, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );

      /* allocat memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, 2) );

      consvars[0] = SCIPgetVarVarbound(scip, cons);
      consvars[1] = SCIPgetVbdvarVarbound(scip, cons);

      vals[0] = 1.0;
      vals[1] = SCIPgetVbdcoefVarbound(scip, cons);

      value = SCIPgetRhsVarbound(scip, cons);

      if ( SCIPisInfinity(scip, value) )
         value = SCIPgetLhsVarbound(scip, cons);

      assert( !SCIPisInfinity(scip, value) && !SCIPisInfinity(scip, -value) );

      SCIP_CALL( getLinearRhs(scip, consvars, vals, 2, transformed, &value) );

      printEntry( scip, file, "RHS", consname, value, &recordcnt, maxnamelen );

      /* free buffer array */
      SCIPfreeBufferArray(scip, &consvars);
      SCIPfreeBufferArray(scip, &vals);
   }

   /* take care of the aggregated constraints for SOS constraints */
   for (c = 0; c < nAggregatedVars; ++c )
   {
      /* output rhs */
      snprintf(consname, MPS_MAX_NAMELEN, "aggr%d", c);
      printEntry(scip, file, "RHS", consname, -aggregatedConstant[c], &recordcnt, maxnamelen);
   }

   if( recordcnt == 1 )
      SCIPinfoMessage(scip, file, "\n");

   if( needRANGES )
   {
      /* print RANGES section */
      SCIPinfoMessage(scip, file, "RANGES\n");
      SCIPdebugMessage("start printing RANGES section\n");
      recordcnt = 0;

      for( c = 0; c < nlinears; ++c  )
      {
         cons = linears[c];
         snprintf(consname, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );

         lhs = SCIPgetLhsLinear(scip, cons);
         rhs = SCIPgetRhsLinear(scip, cons);

         if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) && !SCIPisEQ(scip, rhs, lhs) )
         {
            assert( SCIPisGT(scip, rhs, lhs) );
            printEntry( scip, file, "RANGE", consname, rhs - lhs, &recordcnt, maxnamelen );
         }
      }

      if(recordcnt == 1 )
         SCIPinfoMessage(scip, file, "\n");
   }

   /* print BOUNDS section */
   SCIPinfoMessage(scip, file, "BOUNDS\n");
   SCIPdebugMessage("start printing BOUNDS section\n");

   for (v = 0; v < nvars; ++v)
   {
      var = vars[v];
      assert( var != NULL );
      snprintf(varname, MPS_MAX_NAMELEN, "%s", SCIPvarGetName(var) );

      if( transformed )
      {
         /* in case the transformed is written only local bounds are posted which are valid in the current node */
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
         printStart(scip, file, "BV", "Bound", maxnamelen);
         printRecord(scip, file, varname, "", maxnamelen);
         SCIPinfoMessage(scip, file, "\n");
         continue;
      }

      /* take care of free variables */
      if ( SCIPisInfinity(scip, -lb) && SCIPisInfinity(scip, ub) )
      {
         /* variable is free */
         printStart(scip, file, "FR", "Bound", maxnamelen);
         printRecord(scip, file, varname, "", maxnamelen);
         SCIPinfoMessage(scip, file, "\n");
         continue;
      }

      /* take care of fixed variables */
      if ( SCIPisEQ(scip, lb, ub) )
      {
         /* variable is fixed */
         snprintf(valuestr, MPS_MAX_VALUELEN+1, "%25.15g", lb);
         printStart(scip, file, "FX", "Bound", maxnamelen);
         printRecord(scip, file, varname, valuestr, maxnamelen);
         SCIPinfoMessage(scip, file, "\n");
         continue;
      }

      /* print lower bound */
      if ( SCIPisInfinity(scip, -lb) )
      {
         assert( !SCIPisInfinity(scip, ub) );
         printStart(scip, file, "MI", "Bound", maxnamelen);
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
            snprintf(valuestr, MPS_MAX_VALUELEN+1, "%25.15g", lb);
            printStart(scip, file, "LO", "Bound", maxnamelen);
            printRecord(scip, file, varname, valuestr, maxnamelen);
            SCIPinfoMessage(scip, file, "\n");
         }
      }

      /* print upper bound as far this one is not infinity */
      if( !SCIPisInfinity(scip, ub) )
      {
         snprintf(valuestr, MPS_MAX_VALUELEN+1, "%25.15g", ub);
         printStart(scip, file, "UP", "Bound", maxnamelen);
         printRecord(scip, file, varname, valuestr, maxnamelen);
         SCIPinfoMessage(scip, file, "\n");
      }
   }
   /* output aggregated variables as 'free' */
   for (v = 0; v < nAggregatedVars; ++v)
   {
      var = aggregatedVars[v];
      assert( var != NULL );
      snprintf(varname, MPS_MAX_NAMELEN, "%s", SCIPvarGetName(var) );

      /* variable is free */
      printStart(scip, file, "FR", "Bound", maxnamelen);
      printRecord(scip, file, varname, "", maxnamelen);
      SCIPinfoMessage(scip, file, "\n");
   }

   /* print SOS section */
   if ( nConsSOS1 > 0 || nConsSOS2 > 0 )
   {
      SCIP_Real* sosweights;

      SCIPinfoMessage(scip, file, "SOS\n");
      SCIPdebugMessage("start printing SOS section\n");

      /* first output SOS1 constraints */
      for (c = 0; c < nConsSOS1; ++c)
      {
	 cons = consSOS1[c];
	 consvars = SCIPgetVarsSOS1(scip, cons);
	 nconsvars = SCIPgetNVarsSOS1(scip, cons);
	 sosweights = SCIPgetWeightsSOS1(scip, cons);
         snprintf(consname, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );

         printStart(scip, file, "S1", consname, -1);
	 SCIPinfoMessage(scip, file, "\n");

	 for (v = 0; v < nconsvars; ++v)
	 {
            snprintf(varname, MPS_MAX_NAMELEN, "%s", SCIPvarGetName(consvars[v]) );
            printStart(scip, file, "", varname, maxnamelen);

	    if ( sosweights != NULL )
               sprintf(valuestr, "%25.15g", sosweights[v]);
	    else
               sprintf(valuestr, "%25d ", v);

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
         snprintf(consname, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );

         printStart(scip, file, "S2", consname, -1);
	 SCIPinfoMessage(scip, file, "\n");

	 for (v = 0; v < nconsvars; ++v)
	 {
            snprintf(varname, MPS_MAX_NAMELEN, "%s", SCIPvarGetName(consvars[v]) );
            printStart(scip, file, "", varname, maxnamelen);

	    if ( sosweights != NULL )
               sprintf(valuestr, "%25.15g", sosweights[v]);
	    else
               sprintf(valuestr, "%25d ", v);

	    SCIPinfoMessage(scip, file, "%25s\n", valuestr);
	 }
      }
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &linears);
   SCIPfreeBufferArray(scip, &setppcs);
   SCIPfreeBufferArray(scip, &logicors);
   SCIPfreeBufferArray(scip, &knapsacks);
   SCIPfreeBufferArray(scip, &varbounds);
   SCIPfreeBufferArray(scip, &aggregatedVars);
   SCIPfreeBufferArray(scip, &aggregatedConstant);
   SCIPhashtableFree(&varAggregatedHash);
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

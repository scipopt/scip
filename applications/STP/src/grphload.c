/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   grphload.c
 * @brief  Methods for loading Steiner problems in .stp format
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 *
 * This file includes methods for reading a Steiner problem in .stp format.
 *
 * A list of all interface methods can be found in grph.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*lint -esym(750,GRPHLOAD_C) -esym(766,errno.h)                              */

#define GRPHLOAD_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>         /* errno */
#include <assert.h>
#include <stdarg.h>        /* message: va_list etc */

#if defined(_MSC_VER)
#include  <io.h>
#endif

#if defined(_WIN32) || defined(_WIN64) || defined(_MSC_VER)
#ifndef R_OK
#define R_OK 1
#endif
#else
#include <unistd.h>        /* R_OK  */
#endif

#include "portab.h"
#include "grph.h"

#define MSG_FATAL    0
#define MSG_ERROR    1
#define MSG_WARN     2
#define MSG_INFO     3
#define MSG_DEBUG    4

#define VERBOSE      MSG_INFO

/* Try to get the maximum length of a path.
 *
 * WARNING: if a path found or build during the scanning process is
 *          longer than defined below, the program will probably
 *          crash, because scanf() will overwrite some memory.
 */
#if defined(PATH_MAX)                     /* Found this on SCO UNIX */
#define MAX_PATH_LEN PATH_MAX
#elif defined(_MAX_PATH)                  /* Watcom C on MSDOS */
#define MAX_PATH_LEN _MAX_PATH
#elif defined(_POSIX_PATH_MAX)            /* Maybe POSIX */
#define MAX_PATH_LEN _POSIX_PATH_MAX
#else
#define MAX_PATH_LEN 1024
#endif

/* Extension separators
 */
#define EXTSEP             '.'

#define MAX_LINE_LEN       1024
#define MAX_KEYWORD_LEN    64
#define MAX_STRING_LEN     256
#define MAX_ARGUMENTS      8

struct key
{
   const char*           keyword;
   int                   sw_code;
   const char*           format;
};

#define KEY_SECTION              0001
#define KEY_EOF                  9999
#define KEY_END                  9998

#define KEY_COMMENT_NAME         1001
#define KEY_COMMENT_DATE         1002
#define KEY_COMMENT_CREATOR      1003
#define KEY_COMMENT_PROBLEM      1004
#define KEY_COMMENT_REMARK       1005

#define KEY_GRAPH_NODES          2001
#define KEY_GRAPH_EDGES          2002
#define KEY_GRAPH_E              2003
#define KEY_GRAPH_A              2004
#define KEY_GRAPH_AA             2005
#define KEY_GRAPH_OBSTACLES      2006
#define KEY_GRAPH_HOPLIMIT       2007

#define KEY_TERMINALS_END        3001
#define KEY_TERMINALS_TERMINALS  3002
#define KEY_TERMINALS_T          3003
#define KEY_TERMINALS_TP         3004
#define KEY_TERMINALS_ROOT       3005
#define KEY_TERMINALS_ROOTP      3006
#define KEY_TERMINALS_TG         3007
#define KEY_TERMINALS_GROUPS     3008
#define KEY_TERMINALS_TR         3009

#define KEY_COORDINATES_DD       4001
#define KEY_COORDINATES_DDD      4002
#define KEY_COORDINATES_DDDD     4003
#define KEY_COORDINATES_DDDDD    4004
#define KEY_COORDINATES_DDDDDD   4005
#define KEY_COORDINATES_DDDDDDD  4006
#define KEY_COORDINATES_DDDDDDDD 4007

#define KEY_COORDINATES_END      4011
#define KEY_COORDINATES_GRID     4012

#define KEY_SOLUTION_VALUE       4021
#define KEY_SOLUTION_DATE        4022
#define KEY_SOLUTION_TIME        4023
#define KEY_SOLUTION_STEINER     4024
#define KEY_SOLUTION_S           4025

#define KEY_PRESOLVE_DATE        5001
#define KEY_PRESOLVE_FIXED       5002
#define KEY_PRESOLVE_LOWER       5003
#define KEY_PRESOLVE_UPPER       5004
#define KEY_PRESOLVE_TIME        5005
#define KEY_PRESOLVE_EA          5006
#define KEY_PRESOLVE_EC          5007
#define KEY_PRESOLVE_ED          5008
#define KEY_PRESOLVE_ES          5009


#define KEY_NODEWEIGHTS_NW       6000
#define KEY_NODEWEIGHTS_END      6001

#define KEY_MAXDEGS_MD           8000

#define KEY_OBSTACLES_RR         9000
#define KEY_OBSTACLES_END        9001

#define KEY_HOPCONS_LIM          10000
#define KEY_HOPCONS_FACTOR       10001

#define KEY_TREE_S               11000

static const struct key keyword_table[] =
   {
      /*
       * *** The keywords MUST be sorted alphabetically ! ***
       */
      {  ".eof",                     KEY_EOF,                    NULL        },
      {  ".section",                 KEY_SECTION,                NULL        },

      {  "comment.creator",          KEY_COMMENT_CREATOR,        "s"         },
      {  "comment.date",             KEY_COMMENT_DATE,           "s"         },
      {  "comment.end",              KEY_END,                    NULL        },
      {  "comment.name",             KEY_COMMENT_NAME,           "s"         },
      {  "comment.problem",          KEY_COMMENT_PROBLEM,        "s"         },
      {  "comment.remark",           KEY_COMMENT_REMARK,         "s"         },

      {  "coordinates.dd",           KEY_COORDINATES_DD,         "nnn"       },
      {  "coordinates.ddd",          KEY_COORDINATES_DDD,        "nnnn"      },
      {  "coordinates.dddd",         KEY_COORDINATES_DDDD,       "nnnnn"     },
      {  "coordinates.ddddd",        KEY_COORDINATES_DDDDD,      "nnnnnn"    },
      {  "coordinates.dddddd",       KEY_COORDINATES_DDDDDD,     "nnnnnnn"   },
      {  "coordinates.ddddddd",      KEY_COORDINATES_DDDDDDD,    "nnnnnnnn"  },
      {  "coordinates.dddddddd",     KEY_COORDINATES_DDDDDDDD,   "nnnnnnnnn" },
      {  "coordinates.end",          KEY_COORDINATES_END,        NULL        },
      {  "coordinates.grid",         KEY_COORDINATES_GRID,       NULL        },

      {  "graph.a",                  KEY_GRAPH_A,                "nnn"       },
      {  "graph.aa",                 KEY_GRAPH_AA,               "nnnn"      },
      {  "graph.e",                  KEY_GRAPH_E,                "nnn"       },
      {  "graph.edges",              KEY_GRAPH_EDGES,            "n"         },
      {  "graph.end",                KEY_END,                    NULL        },
      {  "graph.hoplimit",           KEY_GRAPH_HOPLIMIT,         "n"         },
      {  "graph.nodes",              KEY_GRAPH_NODES,            "n"         },
      {  "graph.obstacles",          KEY_GRAPH_OBSTACLES,        "n"         },

      {  "hopconstraint.limit",      KEY_HOPCONS_LIM,            "n"         },
      {  "hopconstraint.factor",     KEY_HOPCONS_FACTOR,         "nn"        },

      {  "maximumdegrees.end",       KEY_END,                    NULL        },
      {  "maximumdegrees.md",        KEY_MAXDEGS_MD,             "n"         },

      {  "nodeweights.end",          KEY_NODEWEIGHTS_END,        NULL        },
      {  "nodeweights.nw",           KEY_NODEWEIGHTS_NW,         "n"         },

      {  "obstacles.end",            KEY_OBSTACLES_END,          NULL        },
      {  "obstacles.rr",             KEY_OBSTACLES_RR,           "nnnn"      },

      {  "presolve.date",            KEY_PRESOLVE_DATE,          "s"         },
      {  "presolve.ea",              KEY_PRESOLVE_ED,            "nnnn"      },
      {  "presolve.ec",              KEY_PRESOLVE_EC,            "nnn"       },
      {  "presolve.ed",              KEY_PRESOLVE_ED,            "nnn"       },
      {  "presolve.end",             KEY_END,                    NULL        },
      {  "presolve.es",              KEY_PRESOLVE_ES,            "nn"        },
      {  "presolve.fixed",           KEY_PRESOLVE_FIXED,         "n"         },
      {  "presolve.lower",           KEY_PRESOLVE_LOWER,         "n"         },
      {  "presolve.orgnodes",        KEY_EOF,                    "n"         },
      {  "presolve.time",            KEY_PRESOLVE_TIME,          "n"         },
      {  "presolve.upper",           KEY_PRESOLVE_UPPER,         "n"         },

      {  "solution.date",            KEY_SOLUTION_DATE,          "s"         },
      {  "solution.end",             KEY_END,                    NULL        },
      {  "solution.s",               KEY_SOLUTION_S,             "n"         },
      {  "solution.steiner",         KEY_SOLUTION_STEINER,       "n"         },
      {  "solution.time",            KEY_SOLUTION_TIME,          "n"         },
      {  "solution.value",           KEY_SOLUTION_VALUE,         "n"         },

      {  "terminals.end",            KEY_TERMINALS_END,          NULL        },
      {  "terminals.groups",         KEY_TERMINALS_GROUPS,       "n"         },
      {  "terminals.root",           KEY_TERMINALS_ROOT,         "n"         },
      {  "terminals.rootp",          KEY_TERMINALS_ROOTP,        "n"         },
      {  "terminals.t",              KEY_TERMINALS_T,            "n"         },
      {  "terminals.terminals",      KEY_TERMINALS_TERMINALS,    "n"         },
      {  "terminals.tg",             KEY_TERMINALS_TG,           "nn"        },
      {  "terminals.tp",             KEY_TERMINALS_TP,           "nn"        },
      {  "terminals.tr",             KEY_TERMINALS_TR,           "nn"        },

      {  "tree.s",                   KEY_TREE_S,                 NULL        },
   };

struct section
{
   const char*           name;
   const char*           extension;
   const int             flag;
   int                   mark;
};

#define FLAG_OPTIONAL     1
#define FLAG_REQUIRED     2

#define SECTION_MISSING    0
#define SECTION_EXISTEND   1

/* Extension NULL = no separate file possible !
 */
static struct section section_table[] =
   {
      { "",            "stp", FLAG_REQUIRED, SECTION_EXISTEND },

      /*
       * *** The section names MUST be sorted alphabetically ! ***
       */
      { "comment",     NULL,  FLAG_OPTIONAL, SECTION_MISSING },
      { "coordinates", "crd", FLAG_OPTIONAL, SECTION_MISSING },
      { "graph",       "grp", FLAG_REQUIRED, SECTION_MISSING },
      { "maximumdegrees", "mdg", FLAG_OPTIONAL, SECTION_MISSING },
      { "nodeweights", "nwg", FLAG_OPTIONAL, SECTION_MISSING },
      { "obstacles",   "obs", FLAG_OPTIONAL, SECTION_MISSING },
      { "presolve",    "prs", FLAG_OPTIONAL, SECTION_MISSING },
      { "solution",    "slt", FLAG_OPTIONAL, SECTION_MISSING },
      { "terminals",   "trm", FLAG_OPTIONAL, SECTION_MISSING },
      { "tree",        "tre", FLAG_OPTIONAL, SECTION_MISSING },
   };

typedef struct current_file
{
   char                  filename[MAX_PATH_LEN];
   int                   line;
   FILE*                 fp;
   struct                section* section;
} CURF;

typedef union parameter
{
   double                n;  /* Could be long long */
   char                  s[MAX_STRING_LEN];
} PARA;

/*---------------------------------------------------------------------------*/
/*--- Name     : String to Lower                                          ---*/
/*--- Function : Converts a string to lower case.                         ---*/
/*--- Arguments: Pointer to string.                                       ---*/
/*--- Returns  : The argument, but now pointing to a lower case string.   ---*/
/*---------------------------------------------------------------------------*/
static
char* strlower(
   char* s
   )
{
   char* t;

   for( t = s; *s != '\0'; s++ )
      *s = (char)tolower(*s);

   return t;
}

/*---------------------------------------------------------------------------*/
/*--- Name     : Print Message                                            ---*/
/*--- Function : Prints a message on stderr.                              ---*/
/*--- Arguments: Type of message, Info about filename and line number,    ---*/
/*---            printf format string and parameters to be printed.       ---*/
/*--- Returns  : Nothing                                                  ---*/
/*---------------------------------------------------------------------------*/
static void message(
   unsigned int          type,
   const CURF*           curf,
   const char*           msg,
   ...)
{
   va_list params;

   const char* intro[] = { "Fatal Error", "Error      ", "Warning    ",
                           "Info       ", "Debug      "
   };
   const char* header  = "*** %s File \"%s\" line %05d: ";

   assert(type           <  sizeof(intro) / sizeof(char*));
   assert(curf           != NULL);
   assert(curf->filename != NULL);
   assert(curf->line     >= 0);
   assert(msg            != NULL);

   va_start(params, msg);

   if (type <= VERBOSE)
   {
      (void)fprintf(stderr, header, intro[type], curf->filename, curf->line);
      (void)vfprintf(stderr, msg, params);
      (void)fputc('\n', stderr);
   }
   va_end(params);
}

/*---------------------------------------------------------------------------*/
/*--- Name     : Key Compare                                              ---*/
/*--- Function : Compares the key with an element of the list.            ---*/
/*--- Parameter: Pointer to key, pointer to element                       ---*/
/*--- Returns  : <0 : key<elem, =0 : key=elem, >0 : key>elem              ---*/
/*---------------------------------------------------------------------------*/
static int key_cmp(
   const void*           key,
   const void*           elem)
{
   assert(key                                != NULL);
   assert(elem                               != NULL);
   assert(((const struct key*)elem)->keyword != NULL);

   return(strcmp((const char*)key, ((const struct key*)elem)->keyword));
}

/*---------------------------------------------------------------------------*/
/*--- Name     : Section Compare                                          ---*/
/*--- Function : Compares the key with an section name.                   ---*/
/*--- Parameter: Pointer to key, pointer to section                       ---*/
/*--- Returns  : <0 : key<sec, =0 : key=sec, >0 : key>sec                 ---*/
/*---------------------------------------------------------------------------*/
static int sec_cmp(
   const void*          key,
   const void*          section)
{
   assert(key                                    != NULL);
   assert(section                                != NULL);
   assert(((const struct section*)section)->name != NULL);

   return(strcmp((const char*)key, ((const struct section*)section)->name));
}

/*---------------------------------------------------------------------------*/
/*--- Name     : Get Arguments                                            ---*/
/*--- Function : Extract the arguments following a keyword.               ---*/
/*--- Parameter: Current file info, argument format string, input line,   ---*/
/*---            pointer to array of MAX_ARGUMENTS PARA's.                ---*/
/*--- Returns  : 0 for success and < 0 for failure.                       ---*/
/*---------------------------------------------------------------------------*/
static int get_arguments(
   const CURF*           curf,
   const char*           format,
   const char*           s,
   PARA*                 para)
{
   const char* err_missmatch_v = "Wrong Syntax";
   const char* msg_hello_ss    = "get_arguments(\"%s\", \"%s\")";

   int missmatch = FALSE;
   int i;
   int decimal_spaces;
   SCIP_Bool is_negative;
   assert(format != NULL);
   assert(s      != NULL);
   assert(para   != NULL);

   message(MSG_DEBUG, curf, msg_hello_ss, format, s);

   /* We try until we run out of input or have enough arguments.
    */
   while((*s != '\0') && (*format != '\0') && !missmatch)
   {
      missmatch = TRUE;

      switch(*format)
      {
      case 'n' :  /* Numeric */
         /* Go to next digit.
          */
         while((*s != '\0') && !isdigit(*s) && (*s != '.') && (*s != '-') )
         {
            s++;
         }
         /* Something left ?
          */
         if (*s != '\0')
         {
            assert(isdigit(*s) || (*s == '.') || (*s == '-'));

            /* Get it.
             */
            para->n = 0;
            decimal_spaces = -1;
            is_negative = FALSE;
            while(isdigit(*s) || (*s == '.') || (*s == '-'))
            {
               if( *s == '.' )
                  decimal_spaces = 0;
               else if( *s == '-' )
                  is_negative = TRUE;
               else if( decimal_spaces != -1 )
                  para->n = para->n + pow(10.0, (double) -(++decimal_spaces)) * (*s - '0');
               else
                  para->n = para->n * 10 + (*s - '0');
               s++;
            }
            if( is_negative )
               para->n = (-1) * para->n;
            missmatch = FALSE;
         }
         break;
      case 's' :  /* String */
         /* Go to the beginning of the string.
          */
         while((*s != '\0') && (*s != '\"'))
            s++;

         /* Someting left ?
          */
         if (*s != '\0')
         {
            assert(*s == '\"');

            /* Get the String.
             */
            i = 0;
            s++;

            while((*s != '\0') && (*s != '\"') && (i < MAX_STRING_LEN - 1))
               para->s[i++] = *s++;

            para->s[i] = '\0';
            missmatch  = FALSE;
         }
         break;
      case 'b' :  /* Bool */
      case 'd' :  /* Date */
      default  :
         /* CONSTCOND */
         assert(FALSE);
         break;
      }
      if (missmatch)
         message(MSG_ERROR, curf, err_missmatch_v);
      else
      {
         para++;
         format++;

      }
   }
   return((*format == '\0') ? SUCCESS : FAILURE);
}

/*---------------------------------------------------------------------------*/
/*--- Name     : Open File                                                ---*/
/*--- Function : Opens a file, processes the file header and secures it's ---*/
/*---            the right file.                                          ---*/
/*--- Arguments: Info about filename and wanted section (file type).      ---*/
/*--- Returns  : 0 for success with curf->fp and curf->line filled or     ---*/
/*---            or < 0 for failure.                                      ---*/
/*---------------------------------------------------------------------------*/
static int open_file(
   CURF* curf,
   unsigned char main_file
)
{
   const char* err_cantopen_s   = "%s.";
   const char* err_noheader_v   = "Wrong file header.";
   const char* err_nomagic_d    = "Wrong Magic-Number %d.";
   const char* err_wrongtype_ss = "Wrong file type. Found %s, expected %s.";
   const char* msg_version_dd   = "Format Version %d.%d.";

   char         type[4];
   unsigned int magic;
   int          version_major;
   int          version_minor;
   int result = FAILURE;

   assert(curf           != NULL);
   assert(curf->filename != NULL);
   assert(curf->section  != NULL);

   /* Prepare the result.
    */
   curf->line = 1;
   curf->fp   = NULL;

   /* reading in the main file? */
   if( main_file )
   {
      char fillname_gr[MAX_STRING_LEN];

      (void)sprintf(fillname_gr, "%s.%s",
            curf->filename, "gr");

      /* try to open .gr */
      if ((curf->fp = fopen(fillname_gr, "r")) != NULL)
      {
         (void)sprintf(curf->filename, "%s", fillname_gr);
         return(SUCCESS);
      }
      else
      {
         /* try to open .stp */
         char fillname_stp[MAX_STRING_LEN];
         (void) sprintf(fillname_stp, "%s.%s", curf->filename, "stp");

         if ((curf->fp = fopen(fillname_stp, "r")) == NULL)
         {
            message(MSG_FATAL, curf, err_cantopen_s, strerror(errno));
            return result;
         }
         (void)sprintf(curf->filename, "%s", fillname_stp);
      }
   }


   /* Try to open the file...
    */
   if (!main_file && (curf->fp = fopen(curf->filename, "r")) == NULL)
      message(MSG_FATAL, curf, err_cantopen_s, strerror(errno));

   /* Read Header...
    */
   else if (fscanf(curf->fp, "%8x %3s File, STP Format Version %2d.%2d \n",
         &magic, type, &version_major, &version_minor) != 4)
      message(MSG_FATAL, curf, err_noheader_v);

   /* Test Magic...
    */
   else if (magic != STP_MAGIC)
      message(MSG_FATAL, curf, err_nomagic_d, magic);

   /* Did we get the right type of file ?
    */
   else if (strcmp(strlower(type), curf->section->extension))
      message(MSG_FATAL, curf, err_wrongtype_ss, type, curf->section->extension);
   else
   {
      /* Succeeded. Just give a warning if the file has a different
       * version number than the reader and hope it will be ok.
       */
      if ((version_major != VERSION_MAJOR) || (version_minor != VERSION_MINOR))
         message(MSG_WARN, curf, msg_version_dd, version_minor, version_major);

      result = SUCCESS;
   }
   return(result);
}

/*---------------------------------------------------------------------------*/
/*--- Name     : Start a new Section                                      ---*/
/*--- Function : Starts a new section, maybe in another file.             ---*/
/*--- Parameter: Path for files, Basename, pointer to actual file info,   ---*/
/*---            Pointer to where to save actual file info, inputline.    ---*/
/*--- Returns  : 0 for success and < 0 for failure.                       ---*/
/*---------------------------------------------------------------------------*/
static int start_section(
   const char* pathname,
   const char* basename,
   CURF*       curf,
   CURF*       save,
   const char* s)
{
   const char* err_missing_v   = "Section name missing";
   const char* err_badsect_s   = "Unknown section name [%s]";
   const char* err_required_s  = "Can't access required file [%s]";

   CURF            temp;
   char            inclname[MAX_PATH_LEN];
   char            sectname[MAX_KEYWORD_LEN];
   char            dummy   [MAX_KEYWORD_LEN];
   int             tokens;
   int             ret = FAILURE;

   assert(pathname != NULL);
   assert(basename != NULL);
   assert(curf     != NULL);
   assert(save     != NULL);
   assert(s        != NULL);

   sectname[0] = '\0';
   inclname[0] = '\0';

   /* Extract names
    */
   if( (tokens = sscanf(s, "%63s %63s %s", sectname, dummy, inclname)) < 1 )
      message(MSG_FATAL, curf, err_missing_v);
   else
   {
      /* Known section ?
       */
      if( strcmp(strlower(sectname),"comments") == 0 )
         sectname[7] = '\0';

      temp.section = (struct section*)bsearch(strlower(sectname),
         &section_table[1],
         (sizeof(section_table) / sizeof(struct section)) - 1,
         sizeof(struct section), sec_cmp);

      if( temp.section == NULL )
         message(MSG_FATAL, curf, err_badsect_s, sectname);
      else
      {
         /* Is this section in a separate file ?
          */
         if( tokens == 1 || (tokens == 2 && strcmp(strlower(sectname),"tree") == 0) )
         {
            curf->section        = temp.section;
            curf->section->mark |= SECTION_EXISTEND;
            ret                  = SUCCESS;
         }
         else
         {
            (void)sprintf(temp.filename, "%s%s%c%s",
               pathname,
               (*inclname == '\0') ? basename : inclname,
               EXTSEP,
               temp.section->extension);
#if defined(_MSC_VER)
            if (_access(temp.filename, R_OK))
#else
            if (access(temp.filename, R_OK))
#endif
            {
               /* We can't access the include file.
                * If the section is optional, we just ignore
                * the whole thing, otherwise we have a problem.
                */
               if (temp.section->flag & FLAG_REQUIRED)
                  message(MSG_FATAL, curf, err_required_s, temp.filename);
               else
               {
                  temp.section = &section_table[0];
                  ret          = SUCCESS;
               }
            }
            else
            {
               if (!open_file(&temp, FALSE) )
               {
                  *save                = *curf;
                  *curf                = temp;
                  curf->section->mark |= SECTION_EXISTEND;

                  ret                  = SUCCESS;
               }
            }
         }
      }
   }
   return(ret);
}

static
SCIP_RETCODE init_coordinates(
   SCIP* scip,
   GRAPH* g,
   PARA* para,
   double*** coordinates,
   int* grid_dim,
   int* termcount,
   int  dim,
   int  nodes
   )
{
   int i;

   if( *coordinates == NULL )
   {
      assert(g == NULL);
      assert(termcount != NULL);
      assert(grid_dim != NULL);
      assert(*termcount == 0);
      assert(nodes > 0);

      *grid_dim = dim;

      /* allocate memory for the coordinate arrays */
      SCIP_CALL( SCIPallocMemoryArray(scip, coordinates, dim) );

      for( i = 0; i < dim; i++ )
	 SCIP_CALL( SCIPallocMemoryArray(scip, &((*coordinates)[i]), nodes) ); /*lint !e866*/
   }

   for( i = 0; i < dim; i++ )
      (*coordinates)[i][*termcount] = (double)para[i + 1].n;

   (*termcount)++;
   return SCIP_OKAY;
}

static
int get_scale_order(
   SCIP_Real number
   )
{
   int ints;
   int order;
   int trail_zeroes;
   int i;
   int length;
   char s;
   char str_number[SCIP_MAXSTRLEN];
   (void)SCIPsnprintf(str_number, SCIP_MAXSTRLEN, "%f", number);
   length = (int) strlen(str_number);
   if( SCIP_MAXSTRLEN < length )
      length = (int) SCIP_MAXSTRLEN;

   for( i = 0; i < length; i++ )
   {
      if( str_number[length - i - 1] != '0' )
         break;
   }
   trail_zeroes = i;

   for( i = 0; i < length; i++ )
   {
      s = str_number[i];
      if( s == '.' )
         break;
   }
   ints = i;
   order = length - ints - trail_zeroes - 1;

   return order;
}


/* scales coordinates in such a way, that they become integer */
static
SCIP_RETCODE scale_coords(
   double** coordinates,
   int*** scaled_coords,
   int* scale_order,
   int nnodes,
   int grid_dim
   )
{
   int i;
   int j;
   int tmp;
   int scale_factor;
   int max_order = 0;

   assert(coordinates != NULL);
   assert(nnodes > 0);
   assert(grid_dim > 1);

   SCIP_CALL( SCIPallocMemoryArray(scip, scaled_coords, grid_dim) );

   for( i = 0; i < grid_dim; i++ )
      for( j = 0; j < nnodes; j++ )
      {
         tmp = get_scale_order(coordinates[i][j]);

         if( max_order < tmp )
         {
            max_order = tmp;
         }
      }

   *scale_order = max_order;
   scale_factor = (int) pow(10.0, (double) max_order);

   for( i = 0; i < grid_dim; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*scaled_coords)[i]), nnodes) ); /*lint !e866*/
      for( j = 0; j < nnodes; j++ )
      {
         (*scaled_coords)[i][j] = (int) (coordinates[i][j] * scale_factor);
      }
   }
   return SCIP_OKAY;
}

/*---------------------------------------------------------------------------*/
/*--- Name     : Steiner Tree Problem Load                                ---*/
/*--- Function : Reads a file in STP format and parses it.                ---*/
/*--- Parameter: Pointer to filename, Pointer to presolve struct          ---*/
/*---------------------------------------------------------------------------*/
SCIP_RETCODE graph_load(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< pointer to store the graph */
   const char*           file,               /**< file to load */
   PRESOL*               presol              /**< presolving struct */
   )
{
   const char*  err_unknown_s   = "Unknown keyword [%s]";
   const char*  err_include_v   = "Include in included file";
   const char*  err_missing_s   = "Required section %s missing";
   const char*  msg_newsect_s   = "Processing Section %s";
   const char*  msg_keyword_sd  = "Found Keyword \"%s\", code = %d";
   const char*  err_badedge_ddd = "Bad edge %d-%d (%d nodes)";
   const char*  err_badroot_dd  = "Bad root %d (%d nodes)";
   const char*  err_baddeg_dd   = "More degree constraints (%d) than nodes (%d)";
   const char*  msg_finish_dddd = "Nodes: %d  Edges: %d  Terminals: %d  Source=%d\n";

   const char*  endofline = "#;\n\r";
   const char*  separator = " \t:=";

   static CURF  curf_null = { "", 0, NULL, NULL };

   GRAPH*       g       = NULL;
   CURF         curf;
   CURF         save;
   PARA         para    [MAX_ARGUMENTS];
   double       nodeweight;
   char         buffer  [MAX_LINE_LEN];
   char         pathname[MAX_PATH_LEN];
   char         basename[MAX_PATH_LEN];
   char         keyword [MAX_KEYWORD_LEN];
   int          stop_input = FALSE;
   int          ret        = FAILURE;
   char*        s;
   char*        t;
   struct key*  p;
   double**     coordinates = NULL;
   int          i;
   int          head;
   int          tail;
   int          tgroups = 0;
   int          grid_dim = -1;
   int          terms = 0;
   int          nodes = 0;
   int          edges = 0;
   int          nwcount = 0;
   int          degcount = 0;
   int          hoplimit = UNKNOWN;
   int          stp_type = -1;
   int          termcount = 0;
   int          nobstacles = -1;
   int          scale_order = 1;
   int          obstacle_counter = 0;
   int**        scaled_coordinates = NULL;
   int**        obstacle_coords = NULL;
   int          transformed = 0;

   assert(file != NULL);

   /* No section loaded so far.
    */
   for( i = 1; i < (int)(sizeof(section_table) / sizeof(section_table[0])); i++ )
      section_table[i].mark = SECTION_MISSING;

   /* Get the names...
    */
   (void)strcpy(pathname, file);

   /* Did we get a path ?
    */
   if( (s = strrchr(pathname, DIRSEP[0])) == NULL )
   {
      /* No, no path.
       */
      (void)strcpy(basename, pathname);
      pathname[0] = '\0';
   }
   else
   {
      /* Yes, there is a path in front.
       */
      s++;
      (void)strcpy(basename, s);
      *s = '\0';
   }

   /* Strip of the extension
    */
   if ((s = strrchr(basename, EXTSEP)) != NULL)
      *s = '\0';

   /* Build filename
    */
   curf.section     = &section_table[0];
   save             = curf_null;

   (void)sprintf(curf.filename, "%s%s",
      pathname,
      basename);

   /* Open the file...
    */
   if (!open_file(&curf, TRUE))
   {
      /* We read while a file is open...
       */
      while((curf.fp != NULL) && !stop_input)
      {
         /* Read a line.
          */
         if ((s = fgets(buffer, (int) sizeof(buffer), curf.fp)) == NULL)
         {
            /* No more lines available, so we can close the file.
             */
            (void)fclose(curf.fp);

            /* If we are in an included file, we come back to our
             * main file. Otherwise fp_save is NULL and than fp will
             * get NULL so we're finished.
             */
            curf = save;
            save = curf_null;

            continue;
         }
         /* Count line number
          */
         curf.line++;

         /* Find the start of the interesting part of an inputline.
          */
         while(isspace(*s))
            s++;

         /* Find the end of the interesting portion of an input line.
          * Either the start of a comment or the final NL or CR.
          * Since all lines but the last must have at least a NL
          * t is nearly never NULL.
          */
         if ((t = strpbrk(s, endofline)) != NULL)
            *t = '\0';

         /* Is there an interesting part left ?
          */
         if (*s == '\0')
            continue;

         /* Build a keyword of form "sectionname.keyword"
          */
         (void)strcpy(keyword, curf.section->name);

         i          = (int)strlen(keyword);
         keyword[i] = '.';

         for(i++;
             (i < MAX_KEYWORD_LEN - 1) && (isalpha(*s) || (*s == '_'));
             i++, s++)
            keyword[i] = (char)tolower(*s);

         keyword[i] = '\0';

         /* Skip junk following the keyword.
          */
         while((*s != '\0') && (strchr(separator, *s) != NULL))
            s++;

         if( strcmp(keyword,"comments") == 0 )
            keyword[i - 1] = '\0';

         /* Did we know the keyword ?
          */
         p = (struct key*)bsearch(keyword, keyword_table,
            sizeof(keyword_table) / sizeof(struct key),
            sizeof(struct key), key_cmp);

         if (p == NULL)
            message(MSG_ERROR, &curf, err_unknown_s, keyword);
         else
         {
            char newformat[5] = "";
            assert(p != NULL);

            message(MSG_DEBUG, &curf, msg_keyword_sd, p->keyword, p->sw_code);

            /* Yes, so lets get the rest of the line if possible
             */
            if( (stp_type == STP_MWCSP || stp_type == STP_RMWCSP) && p->format != NULL && (p->sw_code == KEY_TERMINALS_T || p->sw_code == KEY_GRAPH_E) )
               strcpy(newformat, "nn");
            else if( stp_type == STP_SAP && p->sw_code == KEY_GRAPH_A )
               strcpy(newformat, "nnnn");

            if( p->format == NULL || !get_arguments(&curf, (const char*)( newformat[0] != '\0' ? newformat : p->format ), s, para) )
            {
               /* Now, what should we do ?
                */
               switch(p->sw_code)
               {
               case KEY_SECTION : /* a new Section starts. */
                  if (save.fp == NULL)
                     stop_input = start_section(pathname, basename,
                        &curf, &save, s);
                  else
                  {
                     message(MSG_FATAL, &curf, err_include_v);
                     stop_input = TRUE;
                  }
                  if (!stop_input)
                     message(MSG_INFO, &curf, msg_newsect_s, curf.section->name);
                  break;
               case KEY_END : /* END found. */
                  curf.section = &section_table[0];
                  break;
               case KEY_TREE_S : /* fall through */
               case KEY_EOF : /* EOF found */
                  ret        = SUCCESS;
                  stop_input = TRUE;

                  /* Test if all required section were found.
                   */
                  for(i = 0; (unsigned)i < sizeof(section_table) / sizeof(struct section); i++)
                  {
                     if ((section_table[i].flag & FLAG_REQUIRED)
                        && !(section_table[i].mark & SECTION_EXISTEND))
                     {
                        message(MSG_FATAL, &curf, err_missing_s, section_table[i].name);
                        ret = FAILURE;
                     }
                  }
                  break;
               case KEY_COMMENT_NAME :
               case KEY_COMMENT_DATE :
               case KEY_COMMENT_CREATOR :
               case KEY_COMMENT_PROBLEM :
                  (void)printf("Problem: [%s]\n", para[0].s);
                  if( strcmp(para[0].s, "SPG") == 0 )
                     stp_type = STP_SPG;
                  else if( strcmp(para[0].s, "PCSPG") == 0
                     || strcmp(para[0].s, "Prize-Collecting Steiner Problem in Graphs") == 0 )
                     stp_type = STP_PCSPG;
                  else if( strcmp(para[0].s, "RPCST") == 0
                     || strcmp(para[0].s, "Rooted Prize-Collecting Steiner Problem in Graphs") == 0 )
                     stp_type = STP_RPCSPG;
                  else if( strcmp(para[0].s, "NWSPG") == 0 )
                     stp_type = STP_NWSPG;
                  else if( strcmp(para[0].s, "DCST") == 0 )
                     stp_type = STP_DCSTP;
                  else if( strcmp(para[0].s, "RSMT") == 0 )
                     stp_type = STP_RSMT;
                  else if( strcmp(para[0].s, "OARSMT") == 0 )
                     stp_type = STP_OARSMT;
                  else if( strcmp(para[0].s, "Maximum Node Weight Connected Subgraph") == 0
                     || strcmp(para[0].s, "MWCS") == 0 )
                     stp_type = STP_MWCSP;
                  else if( strcmp(para[0].s, "Rooted Maximum Node Weight Connected Subgraph") == 0
                     || strcmp(para[0].s, "RMWCS") == 0 )
                     stp_type = STP_RMWCSP;
                  else if( strcmp(para[0].s, "HCDST") == 0 )
                     stp_type = STP_DHCSTP;
                  else if( strcmp(para[0].s, "SAP") == 0 )
                     stp_type = STP_SAP;
                  else if( strcmp(para[0].s, "GSTP") == 0 )
                     stp_type = STP_GSTP;
                  break;
               case KEY_COMMENT_REMARK :
                  (void)printf("Comment: [%s]\n", para[0].s);
                  if( strcmp(para[0].s, "Transformed") == 0 )
                     transformed = 1;
                  break;
               case KEY_GRAPH_NODES :
                  assert(para != NULL);
                  nodes = (int)para[0].n;
                  break;
               case KEY_GRAPH_OBSTACLES :
                  nobstacles = (int)para[0].n;
                  if( nobstacles > 0 )
                     stp_type = STP_OARSMT;
                  break;
               case KEY_GRAPH_HOPLIMIT :
                  hoplimit = (int)para[0].n;
                  stp_type = STP_DHCSTP;
                  break;
               case KEY_GRAPH_EDGES :
                  edges = (int)para[0].n;
                  break;
               case KEY_GRAPH_A :
               case KEY_GRAPH_AA :
               case KEY_GRAPH_E :
                  if( (int)para[0].n > nodes || (int)para[1].n > nodes )
                  {
                     message(MSG_FATAL, &curf, err_badedge_ddd,
                        (int)para[0].n, (int)para[1].n, nodes);
                     ret = FAILURE;
                     break;
                  }

                  if( g == NULL )
                  {
                     if( stp_type == STP_GSTP )
                        SCIP_CALL( graph_init(scip, graph, nodes * 2, edges * 2 + nodes * nodes, 1) );
                     else
                        SCIP_CALL( graph_init(scip, graph, nodes, edges * 2, 1) );

                     g = *graph;
                     assert(g != NULL);
                     assert(g->source == UNKNOWN);
                     for( i = 0; i < nodes; i++ )
                        graph_knot_add(g, -1);

                     if( stp_type == STP_DHCSTP )
                     {
                        assert(hoplimit != UNKNOWN);
                        g->hoplimit = hoplimit;
                     }
                  }

                  if( stp_type == STP_DHCSTP )
                  {
                     tail = (int)para[0].n - 1;
                     head = (int)para[1].n - 1;
                     /* check whether the anti-parallel arc has already been added */
                     for( i = g->inpbeg[head]; i != EAT_LAST; i = g->ieat[i] )
                        if( g->tail[i] == tail )
                           break;
                     if( i == EAT_LAST )
                        graph_edge_add(scip, g, tail, head, (double)para[2].n, FARAWAY);
                     else
                        g->cost[i] = (double)para[2].n;

                  }
                  else if( stp_type == STP_SAP )
                  {
                     graph_edge_add(scip, g, (int)para[0].n - 1, (int)para[1].n - 1, (double)para[2].n, (double)para[3].n);
                  }
                  else if( stp_type == STP_MWCSP || stp_type == STP_RMWCSP )
                  {
                     graph_edge_add(scip, g, (int)para[0].n - 1, (int)para[1].n - 1, 0.0, 0.0);
                  }
                  else
                  {
                     graph_edge_add(scip, g, (int)para[0].n - 1, (int)para[1].n - 1,
                        (double)para[2].n,
                        (p->sw_code == KEY_GRAPH_E)
                        ? (double)para[2].n
                        : (double)para[3].n);
                  }
                  break;
               case KEY_MAXDEGS_MD :
                  assert(g != NULL);
                  assert((int)para[0].n >= 0);

                  if( degcount < nodes )
                  {
                     if( g->maxdeg == NULL )
                     {
                        SCIP_CALL( SCIPallocMemoryArray(scip, &(g->maxdeg), nodes ) );
                        stp_type = STP_DCSTP;
                     }
                     g->maxdeg[degcount++] = (int)para[0].n;
                  }
                  else
                  {
                     message(MSG_FATAL, &curf, err_baddeg_dd,
                        degcount, nodes);
                     ret = FAILURE;
                  }
                  break;
               case KEY_NODEWEIGHTS_NW :
                  nodeweight = (double) para[0].n;
                  assert(g != NULL);
                  assert(presol != NULL);

                  if( stp_type != STP_NWSPG )
                     stp_type = STP_NWSPG;

                  if( Is_term(g->term[nwcount]) )
                     presol->fixed += nodeweight;
                  else
                     /* add node-weight to edge-weights of all incoming edges */
                     for( i = g->inpbeg[nwcount]; i != EAT_LAST; i = g->ieat[i] )
                        g->cost[i] += nodeweight;
                  nwcount++;
                  break;
               case KEY_NODEWEIGHTS_END :
                  curf.section = &section_table[0];
                  break;
               case KEY_OBSTACLES_RR :
                  assert(nobstacles > 0);
                  if( obstacle_coords == NULL )
                  {
                     assert(obstacle_counter == 0);
                     SCIP_CALL( SCIPallocBufferArray(scip, &obstacle_coords, 4) );
                     for( i = 0; i < 4; i++ )
                        SCIP_CALL( SCIPallocBufferArray(scip, &(obstacle_coords[i]), nobstacles) ); /*lint !e866*/
                  }
                  for( i = 0; i < 4; i++ )
                     obstacle_coords[i][obstacle_counter] = (int)para[i].n;
                  obstacle_counter++;
                  break;
               case KEY_OBSTACLES_END :
                  curf.section = &section_table[0];

                  if( obstacle_counter != nobstacles )
                  {
                     message(MSG_FATAL, &curf, "obstacle number does not match coordinates \n");
                     ret = FAILURE;
                     break;
                  }
                  if( scaled_coordinates == NULL )
                  {
                     message(MSG_FATAL, &curf, "coordinates not given \n");
                     ret = FAILURE;
                     break;
                  }
                  assert(g == NULL);

                  // todo fix problem with edges over obstacles
                  message(MSG_FATAL, &curf, "Obstacle avoiding RSMT problems are currently not supported in this format \n");

                  SCIP_CALL( graph_obstgrid_create(scip, graph, scaled_coordinates, obstacle_coords, nodes, grid_dim, nobstacles, scale_order) );
                  g = *graph;
                  if( obstacle_coords != NULL )
                     for( i = 3; i >= 0; i-- )
                        SCIPfreeBufferArrayNull(scip, &(obstacle_coords[i]));
                  SCIPfreeBufferArrayNull(scip, &(obstacle_coords));
                  break;
               case KEY_TERMINALS_END :
                  if( transformed == 0 )
                  {
                     if( stp_type == STP_RMWCSP )
                     {
                        assert(nodes == termcount);
                        if( g != NULL )
                        {
                           SCIP_CALL( graph_pc_2rmw(scip, g) );
                        }
                        else
                        {
                           message(MSG_FATAL, &curf, "graph not initialized \n");
                           ret = FAILURE;
                           break;
                        }
                     }
                     else if( stp_type == STP_MWCSP )
                     {
                        assert(nodes == termcount);
                        if( g != NULL )
                        {
                           SCIP_CALL( graph_pc_2mw(scip, g, g->prize) );
                        }
                        else
                        {
                           message(MSG_FATAL, &curf, "graph not initialized \n");
                           ret = FAILURE;
                           break;
                        }
                     }
                     else if( stp_type == STP_PCSPG )
                     {
                        SCIP_CALL( graph_pc_2pc(scip, g) );
                     }
                     else if( stp_type == STP_RPCSPG )
                     {
                        SCIP_CALL( graph_pc_2rpc(scip, g) );
                     }
                  }
                  curf.section = &section_table[0];
                  break;
               case KEY_TERMINALS_TERMINALS :
                  terms = (int)para[0].n;
                  assert(terms > 0);

                  if( stp_type == STP_MWCSP || stp_type == STP_RMWCSP )
                  {
                     assert(terms == nodes);
                     assert(g != NULL);
                     if( g->prize == NULL )
                        SCIP_CALL( graph_pc_init(scip, g, terms, -1) );
                  }
                  break;
               case KEY_TERMINALS_GROUPS :
                  assert(stp_type == STP_GSTP);
                  tgroups = (int)para[0].n;
                  presol->fixed -= tgroups * 1e+8;
                  for( i = 0; i < tgroups; i++ )
                  {
                     graph_knot_add(g, 0);
                  }
                  break;
               case KEY_TERMINALS_ROOT :
                  assert(g != NULL);

                  if ((int)para[0].n <= nodes)
                  {
                     g->source = (int)para[0].n - 1;
                     graph_knot_chg(g, (int)para[0].n - 1, 0);
                  }
                  else
                  {
                     message(MSG_FATAL, &curf, err_badroot_dd,
                        (int)para[0].n, nodes);
                     ret = FAILURE;
                  }
                  break;
               case KEY_TERMINALS_ROOTP :
                  assert(g != NULL);
                  assert(terms > 0);
                  g->source = (int)para[0].n - 1;
                  graph_knot_chg(g, (int)para[0].n - 1, 0);
                  stp_type = STP_RPCSPG;
                  if( g->prize == NULL )
                     SCIP_CALL( graph_pc_init(scip, g, nodes, -1) );

                  g->prize[(int)para[0].n - 1] = FARAWAY;
                  break;
               case KEY_TERMINALS_T :
                  if( stp_type == STP_MWCSP || stp_type == STP_RMWCSP )
                  {
                     assert(g != NULL);
                     assert(g->prize != NULL);
                     g->prize[(int)para[0].n - 1] = (double)para[1].n;
                     if( SCIPisGT(scip, (double)para[1].n, 0.0) )
                        presol->fixed -= (double)para[1].n;
                     termcount++;
                  }
                  else
                     graph_knot_chg(g, (int)para[0].n - 1, 0);
                  break;
               case KEY_TERMINALS_TG :
                  assert(g != NULL);
                  assert(tgroups > 0);
                  graph_edge_add(scip, g, (int)para[0].n - 1, g->knots - tgroups + (int)para[1].n - 1, 1e+8,  1e+8);
                  assert(Is_term(g->term[g->knots - tgroups + (int)para[1].n - 1]));
                  break;
               case KEY_TERMINALS_TP :
                  assert(g != NULL);
                  graph_knot_chg(g, (int)para[0].n - 1, 0);
                  if( g->prize == NULL )
                  {
                     assert(stp_type != STP_RPCSPG);
                     stp_type = STP_PCSPG;
                     SCIP_CALL( graph_pc_init(scip, g, nodes, -1) );
                  }
                  g->prize[(int)para[0].n - 1] = (double)para[1].n;
                  termcount++;
                  break;
               case KEY_TERMINALS_TR :
                  assert(stp_type == STP_RMWCSP);
                  assert(g != NULL);
                  assert(g->prize != NULL);
                  g->prize[(int)para[0].n - 1] = FARAWAY;
                  presol->fixed -= (double)para[1].n;
                  graph_knot_chg(g, (int)para[0].n - 1, 0);
                  termcount++;
                  break;
               case KEY_COORDINATES_DD :
                  /* in this case coordinates are not needed */
                  if( terms > 0 )
                  {
                     ret        = SUCCESS;
                     stop_input = TRUE;
                     break;
                  }
                  SCIP_CALL( init_coordinates(scip, g, para, &coordinates, &grid_dim, &termcount, 2, nodes) );
                  break;
               case KEY_COORDINATES_DDD :
                  SCIP_CALL( init_coordinates(scip, g, para, &coordinates, &grid_dim, &termcount, 3, nodes) );
                  break;
               case KEY_COORDINATES_DDDD :
                  SCIP_CALL( init_coordinates(scip, g, para, &coordinates, &grid_dim, &termcount, 4, nodes) );
                  break;
               case KEY_COORDINATES_DDDDD :
                  SCIP_CALL( init_coordinates(scip, g, para, &coordinates, &grid_dim, &termcount, 5, nodes) );
                  break;
               case KEY_COORDINATES_DDDDDD :
                  SCIP_CALL( init_coordinates(scip, g, para, &coordinates, &grid_dim, &termcount, 6, nodes) );
                  break;
               case KEY_COORDINATES_DDDDDDD :
                  SCIP_CALL( init_coordinates(scip, g, para, &coordinates, &grid_dim, &termcount, 7, nodes) );
                  break;
               case KEY_COORDINATES_DDDDDDDD :
                  SCIP_CALL( init_coordinates(scip, g, para, &coordinates, &grid_dim, &termcount, 8, nodes) );
                  break;
               case KEY_COORDINATES_END :
                  assert(g == NULL);
                  assert(grid_dim > 1);

                  curf.section = &section_table[0];
                  if( termcount != nodes )
                  {
                     message(MSG_FATAL, &curf, "node number does not match coordinates \n");
                     ret = FAILURE;
                     break;
                  }

                  /* scale all coordinates such that they are integers */
                  SCIP_CALL( scale_coords(coordinates, &scaled_coordinates, &scale_order, nodes, grid_dim) );

                  if( coordinates != NULL )
                     for( i = 0; i < grid_dim; i++ )
                        SCIPfreeMemoryArrayNull(scip, &(coordinates[i]));

                  SCIPfreeMemoryArrayNull(scip, &coordinates);

                  if( stp_type != STP_OARSMT )
                  {
                     SCIP_CALL( graph_grid_create(scip, graph, scaled_coordinates, nodes, grid_dim, scale_order) );
                     g = *graph;
                  }

                  break;
               case KEY_COORDINATES_GRID :
                  break;
               case KEY_PRESOLVE_FIXED :
                  if (presol != NULL)
                     presol->fixed += (double)para[0].n;
                  break;
               case KEY_PRESOLVE_DATE :
                  (void)printf("Found presolve information %s\n",
                     para[0].s);
                  break;
               case KEY_PRESOLVE_LOWER :
                  if (presol != NULL)
                     presol->lower = (double)para[0].n;
                  break;
               case KEY_PRESOLVE_UPPER :
                  if (presol != NULL)
                     presol->upper = (double)para[0].n;
                  break;
               case KEY_PRESOLVE_TIME :
                  if (presol != NULL)
                     presol->time = (int)para[0].n;
                  break;
               case KEY_PRESOLVE_EA :
                  break;
               case KEY_PRESOLVE_EC :
                  break;
               case KEY_PRESOLVE_ED :
                  break;
               case KEY_PRESOLVE_ES :
                  break;
               default :
                  /* CONSTCOND */
                  assert(FALSE);
                  break;
               }
            }
         }
      }
   }
   /* Was there an error in an included file ?
    * If so, close the main file.
    */
   if( save.fp != NULL )
      (void)fclose(save.fp);

   /* Close the actual file anyway. Since we stop at encountering
    * a line with "EOF" on it, this will be the normal case.
    */
   if( curf.fp != NULL )
      (void)fclose(curf.fp);

   if( ret == SUCCESS )
   {
      assert(g != NULL);

      if( g->source == UNKNOWN )
      {
         for( i = 0; i < g->knots; i++ )
            if ((g->term[i] == 0)
               && ((g->source < 0) || (g->grad[i] > g->grad[g->source])))
               g->source = i;
      }

      if( g->stp_type == UNKNOWN )
      {
         if( stp_type != UNKNOWN )
            g->stp_type = stp_type;
         else
            g->stp_type = STP_SPG;
      }

      (void)printf(msg_finish_dddd,
         g->knots, g->edges, g->terms, g->source);

      assert(graph_valid(g));
      return SCIP_OKAY;
   }
   else
   {
      if( obstacle_coords != NULL )
      {
         for( i = 0; i < 4; i++ )
            SCIPfreeBufferArrayNull(scip, &(obstacle_coords[i]));
         SCIPfreeBufferArrayNull(scip, &(obstacle_coords));
      }
      return SCIP_READERROR;
   }

}

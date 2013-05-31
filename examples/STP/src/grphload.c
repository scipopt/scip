/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   Type....: Function                                                      */
/*   File....: grphload.c                                                    */
/*   Name....: Graph File Loader                                             */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*lint -esym(750,GRPHLOAD_C) -esym(766,errno.h)                              */

#define GRPHLOAD_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>         /* errno */
#include <assert.h>
#include <stdarg.h>        /* message: va_list etc */

#include <unistd.h>        /* R_OK  */

#include "portab.h"
#include "crc32.h"
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
 *          longer than defined blow the programm will probably
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
   const char* keyword;
   int         sw_code;
   const char* format;
};

#define KEY_SECTION              0001
#define KEY_EOF                  9999
#define KEY_END                  9998

#define KEY_COMMENT_NAME         1001
#define KEY_COMMENT_DATE         1002
#define KEY_COMMENT_CREATOR      1003
#define KEY_COMMENT_REMARK       1004

#define KEY_GRAPH_NODES          2001
#define KEY_GRAPH_EDGES          2002
#define KEY_GRAPH_E              2003
#define KEY_GRAPH_A              2004

#define KEY_TERMINALS_TERMINALS  3001
#define KEY_TERMINALS_T          3002
#define KEY_TERMINALS_ROOT       3003

#define KEY_COORDINATES_DD       4001
#define KEY_COORDINATES_DDD      4002
#define KEY_COORDINATES_GRID     4003

#define KEY_SOLUTION_VALUE       4001
#define KEY_SOLUTION_DATE        4002
#define KEY_SOLUTION_TIME        4003
#define KEY_SOLUTION_STEINER     4004
#define KEY_SOLUTION_S           4005

#define KEY_PRESOLVE_DATE        5001
#define KEY_PRESOLVE_FIXED       5002
#define KEY_PRESOLVE_LOWER       5003
#define KEY_PRESOLVE_UPPER       5004
#define KEY_PRESOLVE_TIME        5005

static const struct key keyword_table[] =
{
   /*
    * *** The keywords MUST be sorted alphabetically ! ***
    */
   {  ".eof",                KEY_EOF,                 NULL   },
   {  ".section",            KEY_SECTION,             NULL   },

   {  "comment.creator",     KEY_COMMENT_CREATOR,     "s"    },
   {  "comment.date",        KEY_COMMENT_DATE,        "s"    },
   {  "comment.end",         KEY_END,                 NULL   },
   {  "comment.name",        KEY_COMMENT_NAME,        "s"    },
   {  "comment.remark",      KEY_COMMENT_REMARK,      "s"    },

   {  "coordinates.dd",      KEY_COORDINATES_DD,      "nnn"  },
   {  "coordinates.ddd",     KEY_COORDINATES_DDD,     "nnnn" },
   {  "coordinates.end",     KEY_END,                 NULL   },
   {  "coordinates.grid",    KEY_COORDINATES_GRID,    NULL   },

   {  "graph.a",             KEY_GRAPH_A,             "nnnn" },
   {  "graph.e",             KEY_GRAPH_E,             "nnn"  },
   {  "graph.edges",         KEY_GRAPH_EDGES,         "n"    },
   {  "graph.end",           KEY_END,                 NULL   },
   {  "graph.nodes",         KEY_GRAPH_NODES,         "n"    },

   {  "presolve.date",       KEY_PRESOLVE_DATE,       "s"    },
   {  "presolve.end",        KEY_END,                 NULL   },
   {  "presolve.fixed",      KEY_PRESOLVE_FIXED,      "n"    },
   {  "presolve.lower",      KEY_PRESOLVE_LOWER,      "n"    },
   {  "presolve.time",       KEY_PRESOLVE_TIME,       "n"    },
   {  "presolve.upper",      KEY_PRESOLVE_UPPER,      "n"    },

   {  "solution.date",       KEY_SOLUTION_DATE,       "s"    },
   {  "solution.end",        KEY_END,                 NULL   },
   {  "solution.s",          KEY_SOLUTION_S,          "n"    },
   {  "solution.steiner",    KEY_SOLUTION_STEINER,    "n"    },
   {  "solution.time",       KEY_SOLUTION_TIME,       "n"    },
   {  "solution.value",      KEY_SOLUTION_VALUE,      "n"    },

   {  "terminals.end",       KEY_END,                 NULL   },
   {  "terminals.root",      KEY_TERMINALS_ROOT,      "n"    },
   {  "terminals.t",         KEY_TERMINALS_T,         "n"    },
   {  "terminals.terminals", KEY_TERMINALS_TERMINALS, "n"    },
};

struct section
{
   const char* name;
   const char* extension;
   const int   flag;
   int         mark;
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
   { "comment",     NULL,  FLAG_REQUIRED, SECTION_MISSING },
   { "coordinates", "crd", FLAG_OPTIONAL, SECTION_MISSING },
   { "graph",       "grp", FLAG_REQUIRED, SECTION_MISSING },
   { "presolve",    "prs", FLAG_OPTIONAL, SECTION_MISSING },
   { "solution",    "slt", FLAG_OPTIONAL, SECTION_MISSING },
   { "terminals",   "trm", FLAG_REQUIRED, SECTION_MISSING },
};

typedef struct current_file
{
   char   filename[MAX_PATH_LEN];
   int    line;
   FILE*  fp;
   struct section* section;
} CURF;

typedef union parameter
{
   double n;  /* Could be long long */
   char   s[MAX_STRING_LEN];
} PARA;

/*---------------------------------------------------------------------------*/
/*--- Name     : String to Lower                                          ---*/
/*--- Function : Converts a string to lower case.                         ---*/
/*--- Arguments: Pointer to string.                                       ---*/
/*--- Returns  : The argument, but now pointing to a lower case string.   ---*/
/*---------------------------------------------------------------------------*/
static char* strlower(
   char* s)
{
   char* t;

   for(t = s; *s != '\0'; s++)
      *s = (char)tolower(*s);

   return(t);
}

/*---------------------------------------------------------------------------*/
/*--- Name     : Print Message                                            ---*/
/*--- Function : Prints a message on stderr.                              ---*/
/*--- Arguments: Type of message, Info about filename and line number,    ---*/
/*---            printf format string and parameters to be printed.       ---*/
/*--- Returns  : Nothing                                                  ---*/
/*---------------------------------------------------------------------------*/
static void message(
   unsigned int type,
   const CURF*  curf,
   const char*  msg,
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
   const void* key,
   const void* elem)
{
   assert(key                                != NULL);
   assert(elem                               != NULL);
   assert(((const struct key*)elem)->keyword != NULL);
#if 0
   (void)fprintf(stderr, "key [%s] elem [%s]\n",
      (const char*)key, ((const struct key*)elem)->keyword);
#endif
   return(strcmp((const char*)key, ((const struct key*)elem)->keyword));
}

/*---------------------------------------------------------------------------*/
/*--- Name     : Section Compare                                          ---*/
/*--- Function : Compares the key with an section name.                   ---*/
/*--- Parameter: Pointer to key, pointer to section                       ---*/
/*--- Returns  : <0 : key<sec, =0 : key=sec, >0 : key>sec                 ---*/
/*---------------------------------------------------------------------------*/
static int sec_cmp(
   const void* key,
   const void* section)
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
   const CURF* curf,
   const char* format,
   const char* s,
   PARA*       para)
{
   const char* err_missmatch_v = "Wrong Syntax";
   const char* msg_hello_ss    = "get_arguments(\"%s\", \"%s\")";

   int missmatch = FALSE;
   int i;

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
            while((*s != '\0') && !isdigit(*s))
               s++;

            /* Someting left ?
             */
            if (*s != '\0')
            {
               assert(isdigit(*s));

               /* Get it.
                */
               para->n = 0;

               while(isdigit(*s))
               {
                  para->n = para->n * 10 + (*s - '0');
                  s++;
               }
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
   CURF* curf)
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

   /* Try to open the file...
    */
   if ((curf->fp = fopen(curf->filename, "r")) == NULL)
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
   const char* err_duplicate_s = "Duplicate Section [%s]";

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
   if ((tokens = sscanf(s, "%63s %63s %s", sectname, dummy, inclname)) < 1)
      message(MSG_FATAL, curf, err_missing_v);
   else
   {
      /* Known section ?
       */
     temp.section = (struct section*)bsearch(strlower(sectname),
         &section_table[1],
         (sizeof(section_table) / sizeof(struct section)) - 1,
         sizeof(struct section), sec_cmp);

      if (temp.section == NULL)
         message(MSG_FATAL, curf, err_badsect_s, sectname);
      else
      {
         if (temp.section->mark & SECTION_EXISTEND)
            message(MSG_FATAL, curf, err_duplicate_s, sectname);
         else
         {
            /* Is this section in a separate file ?
             */
            if (tokens == 1)
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

               if (access(temp.filename, R_OK))
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
                  if (!open_file(&temp))
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
   }
   return(ret);
}

/*---------------------------------------------------------------------------*/
/*--- Name     : Steiner Tree Problem Load                                ---*/
/*--- Function : Reads a file in STP format and parses it.                ---*/
/*--- Parameter: Pointer to filename, Pointer to presolve struct          ---*/
/*--- Returns  : 0 for success, < 0 for failure                           ---*/
/*---------------------------------------------------------------------------*/
GRAPH* graph_load(
   const char*  file,
   PRESOL*      presol)
{
   const char*  err_unknown_s   = "Unknown keyword [%s]";
   const char*  err_include_v   = "Include in included file";
   const char*  err_missing_s   = "Required section %s missing";
   const char*  msg_newsect_s   = "Processing Section %s";
   const char*  msg_keyword_sd  = "Found Keyword \"%s\", code = %d";
   const char*  err_badedge_ddd = "Bad edge %d-%d (%d nodes)";
   const char*  err_badroot_dd  = "Bad root %d (%d nodes)";
   const char*  msg_finish_dddd = "Knots: %d  Edges: %d  Terminals: %d  Source=%d\n";

   const char*  endofline = "#;\n\r";
   const char*  separator = " \t:=";

   static CURF  curf_null = { "", 0, NULL, NULL };

   GRAPH*       g       = NULL;
   CURF         curf;
   CURF         save;
   PARA         para    [MAX_ARGUMENTS];
   char         buffer  [MAX_LINE_LEN];
   char         pathname[MAX_PATH_LEN];
   char         basename[MAX_PATH_LEN];
   char         keyword [MAX_KEYWORD_LEN];
   int          stop_input = FALSE;
   int          ret        = FAILURE;
   unsigned int crc        = 0;
   char*        s;
   char*        t;
   struct key*  p;
   int          i;
   int          nodes = 0;
   int          edges = 0;
   int          has_coordinates = FALSE;
   int          is_gridgraph    = FALSE;

   assert(file != NULL);

   /* No section loaded so far.
    */
   for(i = 1; i < (int)(sizeof(section_table) / sizeof(section_table[0])); i++)
      section_table[i].mark = SECTION_MISSING;

   /* Get the names...
    */
   (void)strcpy(pathname, file);

   /* Did we get a path ?
    */
   if ((s = strrchr(pathname, DIRSEP[0])) == NULL)
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

   (void)sprintf(curf.filename, "%s%s.%s",
      pathname,
      basename,
      curf.section->extension);

   /* Open the file...
    */
   if (!open_file(&curf))
   {
      /* We read while a file is open...
       */
      while((curf.fp != NULL) && !stop_input)
      {
         /* Read a line.
          */
         if ((s = fgets(buffer, sizeof(buffer), curf.fp)) == NULL)
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

         /* Computation of Checksum
          */
         crc = crc32((unsigned char*)s, crc);

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

         /* Did we know the keyword ?
          */
         p = (struct key*)bsearch(keyword, keyword_table,
            sizeof(keyword_table) / sizeof(struct key),
            sizeof(struct key), key_cmp);

         if (p == NULL)
            message(MSG_ERROR, &curf, err_unknown_s, keyword);
         else
         {
            assert(p != NULL);

            message(MSG_DEBUG, &curf, msg_keyword_sd, p->keyword, p->sw_code);

            /* Yes, so lets get the rest of the line if possible
             */
            if ((p->format == NULL)
               || !get_arguments(&curf, p->format, s, para))
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
                     message(MSG_INFO, &curf, msg_newsect_s,
                        curf.section->name);

                  /* Reset CRC
                   */
                  crc = 0;
                  break;
               case KEY_END : /* END found. */
                  message(MSG_INFO, &curf, "CRC [%X]", crc);

                  curf.section = &section_table[0];
                  break;
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
               case KEY_COMMENT_REMARK :
                  (void)printf("Comment: [%s]\n", para[0].s);
                  break;
               case KEY_GRAPH_NODES :
                  nodes = (int)para[0].n;
                  break;
               case KEY_GRAPH_EDGES :
                  edges = (int)para[0].n;
                  break;
               case KEY_GRAPH_A :
               case KEY_GRAPH_E :
                  if (g == NULL)
                  {
                     g = graph_init(nodes, edges * 2, 1, 0);

                     assert(g != NULL);

                     for(i = 0; i < nodes; i++)
                        graph_knot_add(g, -1, 0, 0);

                     g->source[0] = -1;
                  }
                  if (((int)para[0].n <= nodes) && ((int)para[1].n <= nodes))
                  {
                     graph_edge_add(g, (int)para[0].n - 1, (int)para[1].n - 1,
                        (double)para[2].n,
                        (p->sw_code == KEY_GRAPH_E)
                        ? (double)para[2].n
                        : (double)para[3].n);
                  }
                  else
                  {
                     message(MSG_FATAL, &curf, err_badedge_ddd,
                        (int)para[0].n, (int)para[1].n, nodes);
                     ret = FAILURE;
                  }
                  break;
               case KEY_TERMINALS_TERMINALS :
                  break;
               case KEY_TERMINALS_ROOT :
                  assert(g != NULL);

                  if ((int)para[0].n <= nodes)
                     g->source[0] = (int)para[0].n - 1;
                  else
                  {
                     message(MSG_FATAL, &curf, err_badroot_dd,
                        (int)para[0].n, nodes);
                     ret = FAILURE;
                  }
                  break;
               case KEY_TERMINALS_T :
                  graph_knot_chg(g, (int)para[0].n - 1, 0,
                     NO_CHANGE, NO_CHANGE);
                  break;
               case KEY_COORDINATES_DD :
                  has_coordinates = TRUE;
                  graph_knot_chg(g, (int)para[0].n - 1, NO_CHANGE,
                     (int)para[1].n, (int)para[2].n);
                  break;
               case KEY_COORDINATES_DDD :
                  break;
               case KEY_COORDINATES_GRID :
                  is_gridgraph = TRUE;
                  break;
               case KEY_PRESOLVE_FIXED :
                  if (presol != NULL)
                     presol->fixed = (double)para[0].n;
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
               default :
                  /* CONSTCOND */
                  assert(FALSE);
                  break;
               }
            }
         }
      }
   }
   /* Was there an error in an incuded file ?
    * Than close the main file.
    */
   if (save.fp != NULL)
      (void)fclose(save.fp);

   /* Close the actual file anyway. Since we stop at encountering
    * a line with "EOF" on it, this will be the normal case.
    */
   if (curf.fp != NULL)
      (void)fclose(curf.fp);

   if (ret == SUCCESS)
   {
      assert(g != NULL);

      if (g->source[0] == -1)
      {
         for(i = 0; i < g->knots; i++)
            if ((g->term[i] == 0)
               && ((g->source[0] < 0) || (g->grad[i] > g->grad[g->source[0]])))
               g->source[0] = i;
      }
      graph_flags(g, (has_coordinates ? GRAPH_HAS_COORDINATES : 0)
                      | (is_gridgraph ? GRAPH_IS_GRIDGRAPH    : 0));

      (void)printf(msg_finish_dddd,
         g->knots, g->edges, g->terms, g->source[0]);

      assert(graph_valid(g));
   }
   return((ret == SUCCESS) ? g : NULL);
}
